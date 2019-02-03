#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <unistd.h>
#include <gmp.h>
#include "mpzutil.h"
#include "bitmap.h"
#include "ntutil.h"
#include "cstd.h"

/*
    Copyright 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// This module is a grab-bag of stuff, much of which has nothing to do with
// large integer arithmetic (mpz).  This module should be split up and refined...


#define _abs(a)					((a)<0?-(a):(a))

#define expr_valid_op(op)		((op) == 'E' || (op) == 'e' || (op) == '^' || (op) == '*' || (op) == '+' || (op) == '-')

static mpz_t mpz_util_primorial;
static int *mpz_util_primes;
static mpz_t _mpz_temp;
#if MPZ_MAX_SMALL_PRIME_FACTOR_INDEX < 256
static unsigned char mpz_small_factors[MPZ_MAX_SMALL_INTEGER+1];	// for composite i the i-th entry is the prime index of the smallest proper prime divisor, 0 if none (<= 172 for i <= 2^20)
#else
static unsigned short mpz_small_factors[MPZ_MAX_SMALL_INTEGER+1];	// for composite i the i-th entry is the prime index of a proper prime divisor, 0 if none (which must be <= 309 for i <= 2^22)
#endif

short *mpz_tiny_sqrts[MPZ_MAX_TINY_PRIME+1];
static short *mpz_tiny_sqrt_table;

#define MPZ_PP_TABSIZE						28

static struct {
	unsigned int i0, i1, i2;						// m0 contains product of primes with index in (i0,i1], m covers (i0,i2] and m = m0*m1
	unsigned int unsused;
	unsigned long b;							// any composite below b contains a prime factor with index <= i2
	mpz_t m0, m;
} _pptab[MPZ_PP_TABSIZE];

unsigned long _mpz_randomf (mpz_t o, mpz_t N, mpz_t factors[], mpz_ftree_t factor_trees[], unsigned long w);
int mpz_pollard_rho (mpz_t d, mpz_t n);

static gmp_randstate_t mpz_util_rands;
static int mpz_util_inited;

void mpz_util_clear ()
{
	register int i;
	
	if ( ! mpz_util_inited ) return;
	mpz_clear (_mpz_temp);
	mpz_clear (mpz_util_primorial);
	for ( i = 0 ; i < MPZ_PP_TABSIZE ; i++ ) { mpz_clear (_pptab[i].m0); mpz_clear (_pptab[i].m); }
	mem_free (mpz_tiny_sqrt_table);
	gmp_randclear (mpz_util_rands);
	mpz_util_inited = 0;
}

void mpz_util_init ()
{
	unsigned  i, j, k, n, p;
	unsigned long x;
	short *sp;

	if ( mpz_util_inited ) return;
	mpz_init (_mpz_temp);

	mpz_util_primes = prime_small_primes ();

	// i-th entry of small factors array holds the index of the smallest prime divisor of i for composite i, 0 o.w.	
	for ( i = MPZ_MAX_SMALL_PRIME_FACTOR_INDEX ; i  ; i-- ) {
		p = mpz_util_primes[i];
		n = MPZ_MAX_SMALL_INTEGER / p;
		for ( j = 2 ; j <= n ; j++ ) mpz_small_factors[j*p] = i;
	}	

	// Compute various primorials and prime products
	mpz_init2 (mpz_util_primorial, 94027);					// the big Kahuna - the product of primes < 2^16
	mpz_set_ui (mpz_util_primorial, 1);
	for ( i = 1 ; i <= MPZ_TINY_PRIMES ; i++ ) mpz_mul_ui (mpz_util_primorial, mpz_util_primorial, mpz_util_primes[i]);
	
	// early gcd gaps are smaller
	for ( j = 0 ; j < MPZ_PP_TABSIZE ; j++ ) {
		if ( j < 2 ) {
			n = 64;
		} else if ( j < 4 ) {
			n = 128;
		} else {
			n = 256;
		}
		mpz_init2 (_pptab[j].m0,n/2*16);   mpz_init2 (_pptab[j].m,n*16);
		mpz_set_ui (_pptab[j].m0, 1);  mpz_set_ui (_pptab[j].m, 1);
		_pptab[j].i0 = i-1;
		if ( i+256 < MPZ_GCD_PRIMES ) {
			_pptab[j].i1 = i+n/2-1;			// note top of range is closed
			_pptab[j].i2 = i+n-1;
			_pptab[j].b = (unsigned long)mpz_util_primes[_pptab[j].i2+1]*(unsigned long)mpz_util_primes[_pptab[j].i2+1];
			for ( k = 0 ; k < n/8 ; k++ ) {
				x = (unsigned long)mpz_util_primes[i]*(unsigned long)mpz_util_primes[i+1]*(unsigned long)mpz_util_primes[i+2]*(unsigned long)mpz_util_primes[i+3];
				mpz_mul_ui (_pptab[j].m0, _pptab[j].m0, x);
				i += 4;
			}
			mpz_set (_pptab[j].m, _pptab[j].m0);
			for ( k = 0 ; k < n/8 ; k++ ) {
				x = (unsigned long)mpz_util_primes[i]*(unsigned long)mpz_util_primes[i+1]*(unsigned long)mpz_util_primes[i+2]*(unsigned long)mpz_util_primes[i+3];
				mpz_mul_ui (_pptab[j].m, _pptab[j].m, x);
				i += 4;
			}
		} else {
			_pptab[j].i2 = MPZ_GCD_PRIMES;
			_pptab[j].i1 = _pptab[j].i0 + (_pptab[j].i2-_pptab[j].i0) / 2;
			_pptab[j].b = 0xFFFFFFFFUL;
			while ( i <= _pptab[j].i1 ) mpz_mul_ui (_pptab[j].m0, _pptab[j].m0, mpz_util_primes[i++]);
			mpz_set (_pptab[j].m, _pptab[j].m0);		
			while ( i <= _pptab[j].i2 ) mpz_mul_ui (_pptab[j].m, _pptab[j].m, mpz_util_primes[i++]);
		}
		mpz_mul (mpz_util_primorial, mpz_util_primorial, _pptab[j].m);
	}
	if ( i <= MPZ_GCD_PRIMES ) { err_printf ("Initialization alignment problem computing prime products, only %d small primes used\n", i);  abort(); }
	
	// create table of sqrts mod p for tiny primes
	mpz_tiny_sqrt_table = (short*) mem_alloc(MPZ_TINY_PRIME_SUM*sizeof(*mpz_tiny_sqrt_table));
	for ( i = 0 ; i < MPZ_TINY_PRIME_SUM ; i++ ) mpz_tiny_sqrt_table[i] = -1;
	sp = mpz_tiny_sqrt_table;
	for ( i = 1 ; i <= MPZ_TINY_PRIMES ; i++ ) {
		p = mpz_util_primes[i];
		mpz_tiny_sqrts[p] = sp;  sp += p;
		for ( j = 0 ; j < p ; j++ ) mpz_tiny_sqrts[p][(j*j)%p] = j;
	}
	
	gmp_randinit_default (mpz_util_rands);
	gmp_randseed_ui (mpz_util_rands, cstd_seed());
	mpz_util_inited = 1;
}

unsigned long ui_phi (unsigned long n)
{
	unsigned long p[MPZ_MAX_UI_PP_FACTORS];
	unsigned long h[MPZ_MAX_UI_PP_FACTORS];
	unsigned long m;
	int i, j, w;
	
	if ( !n ) { err_printf ("Attempt to compute phi(0)!\n"); abort(); }
	w = ui_factor(p,h,n);
	for ( m = 1, i = 0 ; i < w ; i++ ) {
		m *= (p[i]-1);
		for ( j = 1 ; j < h[i] ; j++ ) m *= p[i];
	}
	return m;
}

int ui_mu (unsigned long n)									// returns \mu(n) (Moebius function)
{
	unsigned long p[MPZ_MAX_UI_PP_FACTORS];
	unsigned long h[MPZ_MAX_UI_PP_FACTORS];
	int i, w;

	if ( !n ) { err_printf ("Attempt to compute mu(0)!\n"); abort(); }
	if ( n==1 ) return 1;
	w = ui_factor(p,h,n);
	for ( i = 0 ; i < w ; i++ ) if ( h[i] > 1 ) return 0;
	return ( (w&1) ? -1 : 1 );
}

static char prime_to_105[105] = { 0,1,1,0,1,0,0,0,1,0,0,1,0,1,0,0,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,1,0,0,0,0,
	                                                 1,1,0,0,0,0,1,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,0,1,1,0,0,1,0,1,1,0,0,1,0,1,0,0,1,0,0,0,1,0,1,1, };
extern unsigned long _ff_p;
  	
int ui_is_prime (unsigned long p)
{
	unsigned long q[MPZ_MAX_UI_PP_FACTORS];
	unsigned long h[MPZ_MAX_UI_PP_FACTORS];
	int sts;
	
	if ( p <= MPZ_MAX_SMALL_PRIME ) return ui_is_small_prime(p);
	sts = fast_prime_check_log (p);
	if ( sts >= 0 ) return sts;
	if ( ! prime_to_105[p%105] ) return 0;
//printf ("ui_is_prime(%ld)  calling factor with p=%ld\n", p, _ff_p);
	if ( ui_factor(q,h,p) > 1 ) return 0;
	if ( h[0] > 1 ) return 0;
	return 1;
}


// returns the least p for which the p-adic valuation of a and b differ
unsigned long ui_pdiff (unsigned long a, unsigned long b)
{
	unsigned long d, p;
	unsigned long q[MPZ_MAX_UI_PP_FACTORS];
	unsigned long h[MPZ_MAX_UI_PP_FACTORS];
	int w;
	int i;
	
	mpz_util_init();
	if ( a == b ) return 0;
	d = ui_gcd(a,b);
	a /= d;
	b /= d;
	// check small primes before factoring
	for ( i = 1 ; i < 25 ; i++ ) {
		if ( a%mpz_util_primes[i] ) {
			if ( ! (b%mpz_util_primes[i]) ) break;
		} else {
			if ( b%mpz_util_primes[i] ) break;
		}
	}
	if ( i < 25 ) return mpz_util_primes[i];
	w = ui_factor (q, h, a);
	if ( w ) p = q[0]; else p = 0;
	w = ui_factor (q, h, b);
	if ( w && p > q[0] ) p = q[0];
	return p;
}


// note that o and a can overlap
void mpz_parallel_invert (mpz_t o[], mpz_t a[], unsigned n, mpz_t p)
{
	static int init;
	static mpz_t c[MPZ_MAX_INVERTS];
	static mpz_t u, v;
	register unsigned i;
	
	if ( ! init ) {
		for ( i = 0 ; i < MPZ_MAX_INVERTS ; i++ ) mpz_init (c[i]);
		mpz_init (u);  mpz_init (v);
		init = 1;
	}
	if ( ! n ) return;
	if ( n > MPZ_MAX_INVERTS ) { err_printf ("Exceeded MPZ_MAX_INVERTS, %d > %d\n", n, MPZ_MAX_INVERTS);  abort(); }
	mpz_set (c[0], a[0]);
	for ( i = 1 ; i < n ; i++ ) mpz_mulm (c[i], c[i-1], a[i], p);
	if ( ! mpz_invert (u, c[n-1], p) ) { err_printf ("Invert failed in mpz_parallel_invert!\n"); }
	for ( i = n-1 ; i ; i-- ) {
		mpz_mulm (v, c[i-1], u, p);
		mpz_mulm (u, u, a[i], p);
		mpz_set (o[i], v);
	}
	mpz_set (o[0], u);
}

// given n integers a[i], computes o = \prod a[i], destroying a[] in the process
void mpz_mul_set (mpz_t o, mpz_t a[], unsigned long n)
{
	register unsigned long i;
	
	if ( ! n ) { mpz_set_ui(o, 1); return; }
	while ( n >= 2 ) {
		for ( i = 0 ; i < (n+1)/2 ; i++ ) {
			if ( 2*i+1 < n ) {
				mpz_mul (a[i], a[2*i], a[2*i+1]);
			} else {
				mpz_set (a[i], a[2*i]);
			}
		}
		n = i;
	}
	if ( n == 2 ) {
		mpz_mul (o, a[0], a[1]);
	} else {
		mpz_set (o, a[0]);
	}
}

// this simple recursive algorithm uses less space but is seriously slow
void ui_crt_coeff_r (mpz_t M, unsigned long a[], unsigned long m[], int n)
{
	mpz_t A0, A1;
	register unsigned long t;
	register int i,j;

	if ( n == 1 ) { a[0] = 1; mpz_set_ui(M,m[0]); return; }
	if ( n == 2 ) { a[0] = m[1]%m[0]; a[1]=m[0]%m[1]; mpz_set_ui(M,m[0]); mpz_mul_ui(M,M,m[1]); return; }
	j = n>>1;
	mpz_init(A0);  mpz_init(A1);
	ui_crt_coeff_r (A0, a, m, j);
	ui_crt_coeff_r (A1, a+j, m+j, n-j);
	for ( i = 0 ; i < j ; i++ ) { t = mpz_fdiv_ui(A1,m[i]); a[i] *= t; }
	for ( ; i < n ; i++ ) { t = mpz_fdiv_ui(A0,m[i]); a[i] *= t; }
	mpz_mul(M,A0,A1);
	mpz_clear(A0); mpz_clear(A1);
}

// given n relatively prime moduli m[i], computes a[i] = \prod_{j!=i}m[j] mod m[i].  w[] is workspace for (n+lg(n)+3) mpz's
// we could make this work with less space, but it will slow things down and isn't worth the complication
void ui_crt_coeff (unsigned long a[], unsigned long m[], int n, mpz_t *w)
{
	mpz_t *e;
	mpz_t *lm[32];
	int lcnt[32];
	register int j, k;
	
//printf ("crt_coeff inputs "); for ( j = 0 ; j < n ; j++ ) printf ("%ld ", m[j]); puts("");
	if ( n==1 ) { a[0] = 1; return; }
	if ( n==2 ) { a[0] = m[1]%m[0];  a[1] = m[0]%m[1];  return; }
	e = w+3;
	// build product tree level 1
	lcnt[0] = n;  lm[1] = e;  k = 1;
	for ( j = 0 ; j+1 < lcnt[k-1] ; j+=2 ) { mpz_set_ui(w[0], m[j]);  mpz_mul_ui(lm[k][j>>1], w[0], m[j+1]); }
	if ( j < lcnt[k-1] ) { mpz_set_ui (lm[k][j>>1],m[j]); j+= 2; }
	lcnt[k] = j>>1;
	e += lcnt[k];
//gmp_printf ("Level %d nodes: ", k); for ( j = 0 ; j < lcnt[k] ; j++ ) gmp_printf ("%Zd ", lm[k][j]); puts("");
	// build product tree levels 2 and above
	for ( k = 2 ; lcnt[k-1] > 2 ; k++ ) {
		lm[k] = e;
		for ( j = 0 ; j+1 < lcnt[k-1] ; j+= 2 ) mpz_mul(lm[k][j>>1],lm[k-1][j],lm[k-1][j+1]);
		if ( j < lcnt[k-1] ) { mpz_set(lm[k][j>>1],lm[k-1][j]); j+=2; }
		lcnt[k] = j>>1;
		e += lcnt[k];
//gmp_printf ("Level %d nodes: ", k); for ( j = 0 ; j < lcnt[k] ; j++ ) gmp_printf ("%Zd ", lm[k][j]); puts("");
	}
	if ( lcnt[--k] != 2 ) { printf ("Error, %d nodes at top of product tree, expected exactly 2\n", lcnt[k]); abort(); }
	mpz_mod(w[0],lm[k][1],lm[k][0]);  mpz_mod(lm[k][1],lm[k][0],lm[k][1]); mpz_set(lm[k][0],w[0]);
//	mpz_mod(la[k][0],lm[k][1],lm[k][0]);  mpz_mod(la[k][1],lm[k][0],lm[k][1]);		// set nodes at top level to contain reduced complements, we must have k>=1
	e += 2;
//gmp_printf ("Top level %d nodes: %Zd %Zd\n", k, la[k][0], la[k][1]);
	// reduce complements down the tree
	for ( k-- ;  k > 0 ; k-- ) {
		for ( j = 0 ; j+1 < lcnt[k] ; j+=2 ) {
			mpz_mod (w[1],lm[k+1][j>>1],lm[k][j]);							// reduce parent complement
			mpz_mul (w[0],w[1],lm[k][j+1]);  mpz_mod (w[2],w[0],lm[k][j]);		// multiply right sibling into left and reduce (enlarges complement to include sibling)
			mpz_mod (w[1],lm[k+1][j>>1],lm[k][j+1]);						// reduce parent complement
			mpz_mul (w[0],w[1],lm[k][j]);  mpz_mod (lm[k][j+1],w[0],lm[k][j+1]);	// multiply left sibling into right and reduce (ditto)
			mpz_set(lm[k][j],w[2]);
		}
		if ( j < lcnt[k] ) mpz_set(lm[k][j],lm[k+1][j>>1]);						// if no sibling, just copy parent value, its already reduced
//gmp_printf ("Level %d nodes: ", k); for ( j = 0 ; j < lcnt[k] ; j++ ) gmp_printf ("%Zd ", la[k][j]); puts("");
	}
	// reduce complements to level 0, which is the output a[]
	for ( j = 0 ; j+1 < lcnt[k] ; j+= 2 ) {
		mpz_mod_ui (w[1],lm[k+1][j>>1],m[j]);											// reduce parent complement
		mpz_mul_ui (w[0],w[1],m[j+1]);  mpz_mod_ui (w[0],w[0],m[j]); a[j] = mpz_get_ui(w[0]);		// multiply right sibling into left and reduce (enlarges complement to include sibling)
		mpz_mod_ui (w[1],lm[k+1][j>>1],m[j+1]);										// reduce parent complement
		mpz_mul_ui (w[0],w[1],m[j]);  mpz_mod_ui (w[0],w[0],m[j+1]); a[j+1] = mpz_get_ui(w[0]);	// multiply left sibling into right and reduce (ditto)		
	}
	if ( j < lcnt[k] ) a[j] = mpz_get_ui(lm[k+1][j>>1]);
/*
{
long bits,lbits,mbits;
bits = mbits = 0;
for ( k = 1 ; k < levels ; k++ ) {
lbits = 0;
for ( j = 0 ; j < lcnt[k] ; j++ ) {
bits += mpz_sizeinbase(lm[k][j],2);
lbits += mpz_sizeinbase(lm[k][j],2);
}
if ( lbits > mbits ) mbits = lbits;
}
printf ("CRT comp used %ld bytes mpz space, could have used just %ld bytes\n", bits/8, 2*mbits/8);
}
*/
//printf ("crt_coeff outputs "); for ( j = 0 ; j < n ; j++ ) printf ("%ld ", a[j]); puts("");
}


#define MPZ_PM_CACHE_SIZE		100

static struct {
	mpz_t p;
	mpz_t maxp;
	mpz_t endp;
	mpz_t o;
	unsigned long maxbits;
	int n;
} _mpz_pm_cache[MPZ_PM_CACHE_SIZE];
unsigned long _mpz_pm_cache_count;


#define MPZ_TIER_SIZE		256

// multiplies o by the product of all primes in (p,maxp] up to maxbits and updates p to last prime used or > maxp if all used.  return number of primes multiplied.
int mpz_prime_mult (mpz_t o, mpz_t p, mpz_t maxp, unsigned long maxbits)
{
	static int init, warn;
	static mpz_t t1[MPZ_TIER_SIZE];
	static mpz_t t2[MPZ_TIER_SIZE];
	static mpz_t t3[MPZ_TIER_SIZE];
	unsigned long bits;
	register int i, i1, i2, i3, n;
	
	if ( ! init ) {
		for ( i = 0 ; i < MPZ_TIER_SIZE ; i++ ) { mpz_init (t1[i]); mpz_init (t2[i]); mpz_init (t3[i]); }
		for ( i = 0 ; i < MPZ_PM_CACHE_SIZE ; i++ ) {
			mpz_init (_mpz_pm_cache[i].p);
			mpz_init (_mpz_pm_cache[i].endp);
			mpz_init (_mpz_pm_cache[i].maxp);
			mpz_init (_mpz_pm_cache[i].o);
		}
		init = 1;
	}
	for ( i = 0 ; i < _mpz_pm_cache_count ; i++ ) {
		if ( mpz_cmp (p,_mpz_pm_cache[i].p) == 0 && mpz_cmp (maxp, _mpz_pm_cache[i].maxp) == 0 &&
		     maxbits == _mpz_pm_cache[i].maxbits ) {
			mpz_set (p, _mpz_pm_cache[i].endp);
			mpz_set (o, _mpz_pm_cache[i].o);
			return _mpz_pm_cache[i].n;
		}
	}
	if ( _mpz_pm_cache_count < MPZ_PM_CACHE_SIZE ) { 
		i = _mpz_pm_cache_count;
		mpz_set (_mpz_pm_cache[i].p, p);
		mpz_set (_mpz_pm_cache[i].maxp, maxp);
		_mpz_pm_cache[i].maxbits = maxbits;
	} else {
		if ( ! warn ) { err_printf ("Prime product cache full in mpz_prime_mult...\n");  warn = 1; }
	}
	i1 = i2 = i3 = 0;
	bits = mpz_sizeinbase (o, 2);
	for ( n = 0 ; bits < maxbits ; n++ ) {
		mpz_nextprime (p, p);
		if ( mpz_cmp (p, maxp) > 0 ) break;
		if ( i3 == MPZ_TIER_SIZE ) {
			if ( i2 == MPZ_TIER_SIZE ) {
				if ( i1 >= MPZ_TIER_SIZE-2 ) break;
				mpz_mul_set (t1[i1++], t2, i2);
				i2 = 0;
			}
			mpz_mul_set (t2[i2++], t3, i3);
			i3 = 0;
		}
		mpz_set (t3[i3++], p);
		bits += mpz_sizeinbase (p, 2);
	}
	if ( i2 == MPZ_TIER_SIZE ) {
		mpz_mul_set (t1[i1++], t2, i2);
		i2 = 0;
	}
	mpz_mul_set (t2[i2++], t3, i3);
	mpz_mul_set (t1[i1++], t2, i2);
	mpz_mul_set (t3[0], t1, i1);
	mpz_mul (o, o, t3[0]);
	if ( _mpz_pm_cache_count < MPZ_PM_CACHE_SIZE ) {
		i = _mpz_pm_cache_count++;
		mpz_set (_mpz_pm_cache[i].endp, p);
		mpz_set (_mpz_pm_cache[i].o, o);
		_mpz_pm_cache[i].n = n;
	}
	return n;
}

// computes the product of all prime powers <= L
// this function is designed to be called once - it dynamically allocates and frees all memory
void mpz_power_primorial (mpz_t o, unsigned long L)	
{
	mpz_t t1[MPZ_TIER_SIZE];
	mpz_t t2[MPZ_TIER_SIZE];
	mpz_t t3[MPZ_TIER_SIZE];
	mpz_t p, q, t;
	unsigned long roots[33];
	register int i, i1, i2, i3, j1, j2, n;
	
	n = ui_len(L);
	if ( n > 32 ) { err_printf ("primorial %u too large to compute - be reasonable.\n", n);  abort(); }
	for ( i = n ; i >= 2 ; i-- ) roots[i] = (unsigned long) floor(pow((double)L,1.0/(double)i));
	roots[1] = L;
	roots[0] = -1;
	
	// just init tier 3 initially, init tier 1 and 2 variables as needed
	for ( i = 0 ; i < MPZ_TIER_SIZE ; i++ ) mpz_init2 (t3[i], MPZ_TIER_SIZE*n);
	j1 = j2 = 0;

	mpz_init (p);  mpz_init (q);  mpz_init (t);  mpz_set_ui (p, 1);
	i1 = i2 = i3 = 0;
	i = n;
	for ( ;; ) {
		mpz_nextprime (p, p);
		while ( i && mpz_cmp_ui (p,roots[i]) > 0 ) i--;
		if ( ! i ) break;
		mpz_pow_ui (q, p, i);
		if ( i3 == MPZ_TIER_SIZE ) {
			if ( i2 == MPZ_TIER_SIZE ) {
				if ( i1 >= MPZ_TIER_SIZE-2 ) { err_printf ("Insufficient memory in mpz_power_primorial on input %lu - increase MPZ_TIER_SIZE or add a tier\n", L);  abort(); }
				if ( i1 == j1 ) mpz_init (t1[j1++]);
				mpz_mul_set (t1[i1++], t2, i2);
				i2 = 0;
			}
			if ( i2 == j2 ) mpz_init (t2[j2++]);
			mpz_mul_set (t2[i2++], t3, i3);
			i3 = 0;
		}
		mpz_set (t3[i3++], q);
	}
	if ( i2 == MPZ_TIER_SIZE ) {
		if ( i1 == j1 ) mpz_init (t1[j1++]);
		mpz_mul_set (t1[i1++], t2, i2);
		i2 = 0;
	}
	// free memory as we go to keep peak usage down
	if ( i2 == j2 ) mpz_init (t2[j2++]);
	mpz_mul_set (t2[i2++], t3, i3);
	for ( i = 0 ; i < MPZ_TIER_SIZE ; i++ ) mpz_clear (t3[i]);
	if ( i1 == j1 ) mpz_init (t1[j1++]);
	mpz_mul_set (t1[i1++], t2, i2);
	for ( i = 0 ; i < j2 ; i++ ) mpz_clear (t2[i]);
	mpz_mul_set (o, t1, i1);
	for ( i = 0 ; i < j1 ; i++ ) mpz_clear (t1[i]);
	mpz_clear (p);  mpz_clear (q);  mpz_clear (t);
	return;
}


int mpz_remove_small_primes (mpz_t o, mpz_t n, unsigned long exps[], unsigned long maxprimes)
{
	static int init;
	static mpz_t d, x, t;
	int i, w;
	
	if ( ! init ) { mpz_util_init();  mpz_init (d);  mpz_init (x);  mpz_init (t);  init = 1; }
	if ( maxprimes > MPZ_GCD_PRIMES ) { err_printf ("maxprimes value %lu exceeded MAX_GCD_PRIMES = %d in mpz_remove_small_primes\n", maxprimes, MPZ_GCD_PRIMES); abort(); }
	mpz_gcd (d, n, mpz_util_primorial);
	mpz_set (o, n);
	exps[0] = 0;
	w = 0;
	for ( i = 1 ; i <= maxprimes ; i++ ) {
		if ( mpz_divisible_ui_p (d, mpz_util_primes[i]) ) {
			mpz_set_ui (x, mpz_util_primes[i]);
			exps[i] = mpz_remove (o, o, x);
			w++;
		} else {
			exps[i] = 0;
		}
	}
	return w;
}

// Let p be the largest prime factor of n.  If n/p <= L, return n/p, otherwise return 0
unsigned long mpz_nearprime (mpz_t n, unsigned long L)
{
	static int init;
	static mpz_t d, x, t;
	int i;
	
	if ( ! init ) { mpz_util_init();  mpz_init (d);  mpz_init (x);  mpz_init (t);  init = 1; }
	if ( L > MPZ_MAX_GCD_INTEGER ) { err_printf ("cofactor bound too large in mpz_nearprime\n"); abort(); }
	mpz_gcd (d, n, mpz_util_primorial);
	if ( mpz_cmp_ui (d, L) > 0 ) return 0;
	mpz_set (t, n);
	for ( i = 1 ; i <= MPZ_GCD_PRIMES ; i++ ) {
		if ( mpz_divisible_ui_p (d, mpz_util_primes[i]) ) {
			mpz_set_ui (x, mpz_util_primes[i]);
			mpz_remove (t, t, x);
			mpz_remove (d, d, x);
			if ( mpz_cmp_ui(d,1) == 0 ) break;
		}
	}
	mpz_divexact (x, n, t);
	if ( mpz_cmp_ui (x, L) > 0 ) return 0;
	if ( mpz_cmp_ui(t,1)==0 ) return mpz_get_ui(x);
	if ( ! mpz_probab_prime_p (t, 5) ) return 0;
	return mpz_get_ui (x);
}


int mpz_remove_small_squares (mpz_t o, mpz_t n)
{
	static int init;
	static mpz_t d, x, m;
	int i, j;
	
	if ( ! init ) { mpz_util_init();  mpz_init (d);  mpz_init (x);   mpz_init (m);  init = 1; }
	mpz_set (o, n);
	if ( mpz_cmp_ui (n, 1) <= 0 ) return 0;
	mpz_gcd (d, n, mpz_util_primorial);
	mpz_set (m, n);
	for ( i = 1 ; i <= MPZ_GCD_PRIMES ; i++ ) {
		if ( mpz_divisible_ui_p (d, mpz_util_primes[i]) ) {
			mpz_set_ui (x, mpz_util_primes[i]);
			j = mpz_remove (m, m, x);
			if ( j%2 ) mpz_mul (m, m, x);
		}
	}
	mpz_set (o, m);
	return 1;
}


// replaces divisors p^n of n by p for primes p <= MPZ_MAX_GCD_PRIME
int mpz_flatten_small (mpz_t o, mpz_t n)
{
	static int init;
	static mpz_t d, x, m;
	int i;
	
	if ( ! init ) { mpz_util_init();  mpz_init (d);  mpz_init (x);   mpz_init (m);  init = 1; }
	mpz_set (o, n);
	if ( mpz_cmp_ui (n, 1) <= 0 ) return 0;
	mpz_gcd (d, n, mpz_util_primorial);
	mpz_set (m, n);
	for ( i = 1 ; i <= MPZ_GCD_PRIMES ; i++ ) {
		if ( mpz_divisible_ui_p (d, mpz_util_primes[i]) ) {
			mpz_set_ui (x, mpz_util_primes[i]);
			mpz_remove (m, m, x);
			mpz_mul (m, m, x);
		}
	}
	mpz_set (o, m);
	return 1;
}


void mpz_print_factors (mpz_t N)
{
	static int init;
	static mpz_t P;
	unsigned long p[MPZ_MAX_SMALL_FACTORS];
	unsigned long h[MPZ_MAX_SMALL_FACTORS];
	int i, w;
	
	if ( ! init ) { mpz_init (P);  init = 1; }
	w = mpz_factor_small (p, h, P, N, MPZ_MAX_SMALL_FACTORS, 128);
	if ( mpz_cmp_ui (P,1) > 0 ) gmp_printf ("%Zd ", P);
	for ( i = 0 ; i < w ; i++ ) {
		if ( h[i] > 1 ) {
			printf ("%lu^%lu ", p[i], h[i]);
		} else if ( h[i] == 1 ) {
			printf ("%lu ", p[i]);			
		}
	}
	puts ("");
}


#define _remove(n,d,c)			while ( !((n)%(d)) ) { (n)/=(d); (c)++; }
#define _tdiv(n,d,h,p,w)			{ _remove(n,d,h[w]); if ( (h)[(w)] ) { (p)[(w)++] = (d); (h)[(w)] = 0; } }

// returns prime-power factors ordered by prime
int ui_factor (unsigned long p[MPZ_MAX_UI_PP_FACTORS], unsigned long h[MPZ_MAX_UI_PP_FACTORS], unsigned long n)
{
	static int init;
	static mpz_t N, D, D1;
	register unsigned long d, d1;
	register int i, k, w;

	if ( ! init ) { if ( ! mpz_util_inited ) mpz_util_init(); mpz_init (N);  mpz_init (D);  mpz_init (D1); mpz_util_init(); init = 1;}
	if ( n == 0 ) { err_printf ("attempt to factor 0\n");  abort(); }
	w = 0;
	h[w] = 0;
	while ( ! (n&0x1) ) { n >>= 1; h[w]++; }
	if ( h[w] ) { p[w++] = 2;  h[w] = 0; }
	if ( n == 1 ) return w;

	// factor small integers using precomputed factor table
	if ( n <= MPZ_MAX_SMALL_INTEGER ) {
		register int sp, tp;

		while ( (tp = mpz_small_factors[n]) ) {
			sp = mpz_util_primes[tp];
			if ( w && p[w-1] == sp ) { h[w-1]++; } else { p[w] = sp; h[w++] = 1; }
			n /= sp;
		}
		if ( w && p[w-1] == n ) { h[w-1]++; } else { p[w] = n; h[w++] = 1; }
		return w;
	}

	// Clear out all tiny primes via trial division.  Yes, it is worth hardwiring all this - it's much faster.
	_tdiv(n,3,h,p,w); _tdiv(n,5,h,p,w); _tdiv(n,7,h,p,w); _tdiv(n,11,h,p,w); _tdiv(n,13,h,p,w); _tdiv(n,17,h,p,w); _tdiv(n,19,h,p,w); _tdiv(n,23,h,p,w); _tdiv(n,29,h,p,w);
	_tdiv(n,31,h,p,w); _tdiv(n,37,h,p,w); _tdiv(n,41,h,p,w); _tdiv(n,43,h,p,w); _tdiv(n,47,h,p,w); _tdiv(n,53,h,p,w); _tdiv(n,59,h,p,w); _tdiv(n,61,h,p,w); _tdiv(n,67,h,p,w); _tdiv(n,71,h,p,w);
	_tdiv(n,73,h,p,w); _tdiv(n,79,h,p,w); _tdiv(n,83,h,p,w); _tdiv(n,89,h,p,w); _tdiv(n,97,h,p,w); _tdiv(n,101,h,p,w); _tdiv(n,103,h,p,w); _tdiv(n,107,h,p,w); _tdiv(n,109,h,p,w); _tdiv(n,113,h,p,w);
	_tdiv(n,127,h,p,w); _tdiv(n,131,h,p,w); _tdiv(n,137,h,p,w); _tdiv(n,139,h,p,w); _tdiv(n,149,h,p,w); _tdiv(n,151,h,p,w); _tdiv(n,157,h,p,w); _tdiv(n,163,h,p,w); _tdiv(n,167,h,p,w); _tdiv(n,173,h,p,w);
	_tdiv(n,179,h,p,w); _tdiv(n,181,h,p,w); _tdiv(n,191,h,p,w); _tdiv(n,193,h,p,w); _tdiv(n,197,h,p,w); _tdiv(n,199,h,p,w); _tdiv(n,211,h,p,w); _tdiv(n,223,h,p,w); _tdiv(n,227,h,p,w); _tdiv(n,229,h,p,w);
	_tdiv(n,233,h,p,w); _tdiv(n,239,h,p,w); _tdiv(n,241,h,p,w); _tdiv(n,251,h,p,w);
	
	// the code above assumes MPZ_MAX_TINY_PRIME is 251!!!
	
	if ( n == 1 ) return w;
	if ( n <= MPZ_MAX_SMALL_INTEGER ) {
		register int sp, tp;

		while ( (tp = mpz_small_factors[n]) ) {
			sp = mpz_util_primes[tp];
			if ( w && p[w-1] == sp ) { h[w-1]++; } else { p[w] = sp; h[w++] = 1; }
			n /= sp;
		}
		if ( w && p[w-1] == n ) { h[w-1]++; } else { p[w] = n; h[w++] = 1; }
		return w;
	}

	// at this point n > MPZ_MAX_SMALL_INTEGER and contains no primes <= MPZ_MAX_TINY_PRIME
	// time to start whacking it with some gcd's
	mpz_set_ui (N, n);
	for ( i = 0 ; i < MPZ_PP_TABSIZE ; i++ ) {
		mpz_gcd (D, N,_pptab[i].m);
		if ( mpz_cmp_ui(D,1) != 0 ) {
			d = mpz_get_ui(D);
			if ( d <= mpz_util_primes[_pptab[i].i2]  ) {		// sliced off just one prime - this is what we want
				p[w] = d;
				h[w] = 0;
				do { n /= d;  h[w]++; } while ( ! (n%d) );	// need to remove all occurrences of d
				w++;
				mpz_set_ui (N, n);
			} else {
//printf ("multiple primes in gcd span\n");
				mpz_gcd (D1, D, _pptab[i].m0);
				if ( mpz_cmp_ui(D1,1) ) {
					d1 = mpz_get_ui (D1);
					if ( d1 <= mpz_util_primes[_pptab[i].i1] ) {
						p[w] = d1;
						h[w] = 0;
						do { n /= d1;  h[w]++; } while ( ! (n%d1) );	// need to remove all occurrences of d1
						w++;
					} else {
//printf ("trial dividing low i=%d, d1=%lu, (%u,%u)\n", i, d1, mpz_util_primes[_pptab[i].i0],mpz_util_primes[_pptab[i].i1]);
						h[w] = 0;
						for ( k = _pptab[i].i0+1 ; k <= _pptab[i].i1 ; k++ ) _tdiv (n, mpz_util_primes[k], h, p, w);
					}
					d1 = d/d1;
				} else {
					d1 = d;
				}
				if ( d1 > 1 ) {
					if ( d1 < mpz_util_primes[_pptab[i].i2]) {
						p[w] = d1;
						h[w] = 0;
						do { n /= d1;  h[w]++; } while ( ! (n%d1) );	// need to remove all occurrences of d1
						w++;
					} else {
//printf ("trial dividing high i=%d, d1=%lu, (%u,%u)\n", i, d1, mpz_util_primes[_pptab[i].i1],mpz_util_primes[_pptab[i].i2]);
						h[w]  = 0;
						for ( k = _pptab[i].i1+1 ; k <= _pptab[i].i2 ; k++ ) _tdiv (n, mpz_util_primes[k], h, p, w);
					}
				}
			}
			if ( n == 1 ) return w;
		}
		if ( n < _pptab[i].b ) {
//printf ("added remainder n=%d\n", n);
			p[w] = n; h[w++] = 1;
			return w;
		}
		mpz_set_ui (N, n);
	}

	// If the original n's second largest prime was below MPZ_SMALL_PRIME, what remains is prime - this is likely for n in the 32-48 bit range, so check here.
	if  ( mpz_probab_prime_p (N, 10) ) {
		p[w] = n; h[w++] = 1;
		return w;
	}
	
	// resort to bigger guns - this can't happen unless N is a composite with more than 32 bits and no prime factors below MPZ_MAX_GCD_PRIME
	// we could go direct to a single-precision version of pollard rho here, but we don't expect this to happen much anyway
	w += mpz_factor_small(p+w,h+w,D,N,MPZ_MAX_UI_PP_FACTORS-w, 64);
	return w;
}

// returns primes in order
int mpz_factor_small (unsigned long p[], unsigned long h[], mpz_t bigp, mpz_t n, int max_factors, int max_hard_bits)
{
	static int init;
	static mpz_t d, x, m;
	unsigned long pt, ht;
	int i, j, w;
	
	if ( ! init ) { mpz_util_init();  mpz_init (d);  mpz_init (x);   mpz_init (m);  init = 1; }
	if ( ! mpz_sgn (n) ) { err_printf ("Attempt to factor 0");  abort(); }
	mpz_set_ui (bigp, 1);
	if ( mpz_sgn(n) < 0 ) mpz_neg (m, n); else mpz_set (m,n);
	if ( mpz_cmp_ui(m,1) == 0 ) return 0;
	mpz_gcd (d, m, mpz_util_primorial);
	w = 0;
	for ( i = 1 ; i <= MPZ_GCD_PRIMES ; i++ ) {
		if ( mpz_cmp_ui (d, 1) == 0 ) break;
		if ( mpz_divisible_ui_p (d, mpz_util_primes[i]) ) {
			p[w] = mpz_util_primes[i];
			mpz_set_ui (x, p[w]);
			h[w] = mpz_remove (m, m, x);
			w++;
			if ( w == max_factors ) { err_printf ("Exceeded max_factors %d in mpz_factor_small\n", max_factors);  abort(); }
			mpz_remove (d, d, x);
		}
	}
	if ( mpz_cmp_ui (m, 1) == 0 ) return w;
	if ( mpz_sizeinbase (m, 2) > max_hard_bits ) return -1;
	while ( ! mpz_probab_prime_p (m, 5) ) {
		mpz_pollard_rho (d, m);
		if ( ! mpz_fits_ulong_p (d) ) { err_printf ("factor too big to fit in unsigned long!");   return -1; }
		pt = mpz_get_ui (d);
		ht = mpz_remove (m, m, d);
		for ( i = 0 ; i < w ; i++ ) if ( p[i] > pt ) break;
		for ( j = w ; j-1 >= i ; j-- ) { p[j] = p[j-1];  h[j] = h[j-1]; }
		p[j] = pt;
		h[j] = ht;			
		w++;
		if ( mpz_cmp_ui (m,1) == 0 ) return w;
		if ( w == max_factors ) { err_printf ("Exceeded max_factors %d in mpz_factor_small\n", max_factors);  abort(); }
	}
	if ( ! mpz_fits_ulong_p (m) ) {
		mpz_set (bigp, m);
	} else {
		mpz_set_ui (bigp, 1);
		pt = mpz_get_ui (m);
		ht = 1;
		for ( i = 0 ; i < w ; i++ ) if ( p[i] > pt ) break;
		for ( j = w ; j-1 >= i ; j-- ) { p[j] = p[j-1];  h[j] = h[j-1]; }
		p[j] = pt;
		h[j] = ht;			
		w++;
	}
	return w;
}


int mpz_pollard_rho (mpz_t d, mpz_t n)
{
	static int init;
	static mpz_t x, y, z, t, P;
	register unsigned c, i, j, k, m;

	if ( ! init ) { mpz_util_init();  mpz_init (x);  mpz_init (y);  mpz_init(t);  mpz_init (P); mpz_init (z); init = 1; }
	m = 0;
	for ( c = 1 ;; c++ ) {
		mpz_set_ui (x, 2);
		mpz_set_ui (y, 2);
		i = 1;
		k = 2;
		for (;;) {
			mpz_set_ui (P, 1);
			mpz_set (t, x);
			for ( j = 0 ; i < k && j < 20 ; j++ ) {
				m++;
				i++;
				mpz_mul (x, x, x);
				mpz_add_ui (x, x, c);
				mpz_mod (x, x, n);
				mpz_sub (z, y, x);
				mpz_mulm (P, P, z, n);
			}
			mpz_gcd (d, P, n);
			if ( mpz_cmp_ui (d,1) != 0 ) {
				mpz_set (x,t);
				do {
					mpz_mul (x, x, x);
					mpz_add_ui (x, x, c);
					mpz_mod (x, x, n);
					mpz_sub (z, y, x);
					mpz_gcd (d, z, n);
				} while ( mpz_cmp_ui (d,1) == 0 );
				if ( mpz_cmp (d,n) != 0 ) return m;
				break;
			}
			if ( i == k ) {
				mpz_set (y,x);
				k *= 2;
			}
		}
	}
	return 1;
}


// computes the y-coarse part of x
int mpz_coarse_part (mpz_t o, mpz_t x, mpz_t y)
{
	static int init;
	static mpz_t d, z;
	int i;
	
	if ( ! init ) { mpz_util_init();  mpz_init (d);  mpz_init (z);   init = 1; }
	if ( mpz_cmp (x,y) <= 0 ) { mpz_set_ui (o,1);  return 1; }
	mpz_gcd (d, x, mpz_util_primorial);
	mpz_set (o, x);
	for ( i = 1 ; i <= MPZ_GCD_PRIMES ; i++ ) {
		if ( mpz_cmp_ui (y, mpz_util_primes[i]) < 0 ) break;
		if ( mpz_divisible_ui_p (d, mpz_util_primes[i]) ) {
			mpz_set_ui (z, mpz_util_primes[i]);
			mpz_remove (o, o, z);
		}
	}
	if ( i <= MPZ_GCD_PRIMES ) return 1;
	if ( mpz_cmp (o,y) <= 0 ) { mpz_set_ui (o, 1);  return 1; }
	mpz_set_ui (d, MPZ_MAX_GCD_PRIME);
	while ( ! mpz_probab_prime_p (o, 5) ) {
		while ( ! mpz_divisible_p (o, d) ) {
			mpz_nextprime (d, d);
			if ( mpz_cmp (d, y) > 0 ) return 1;
		}
		mpz_remove (o, o, d);
	}
	return 1;
}

#define MAX_MERSENNE			20
static int mersenne_tab[MAX_MERSENNE+1] = { 1, 2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423 };
 
int mpz_eval_expr (mpz_t o, char *expr)
{
	static int init;
	static mpz_t p, x, y, z;
    char *s, *t, op, nextop;
    int i, digits, n;

 	if ( ! init ) { mpz_init (p);  mpz_init (x);  mpz_init (y);  mpz_init (z);  init = 1; }
    	mpz_util_init();	
    s = expr;
    if ( *s == 'D' || *s == 'd' ) {
        digits = atoi (s+1);
		mpz_ui_pow_ui (x, 10, digits);
		do {
			mpz_urandomm (o, mpz_util_rands, x);
		} while ( ! mpz_tstbit (o, 0) );
		if ( mpz_congruent_ui_p (o, 1, 4) ) mpz_mul_2exp (o, o, 2);
		gmp_info_printf ("Random Discriminant: %Zd\n", o);
		return 1;
    }
    if ( *s == 'M' || *s == 'm' ) {
        digits = atoi (s+1);
	for ( i = 0 ; i <= MAX_MERSENNE ; i++ ) if ( digits == mersenne_tab[i] ) break;
	if ( i > MAX_MERSENNE ) err_printf ("Warning, specified Mersenne prime unknown\n");
		mpz_ui_pow_ui (x, 2, digits);
		mpz_sub_ui(o,x,1);
		info_printf ("Mersenne prime: 2^%d-1\n", digits);
		return 1;
    }
    if ( *s == 'R' || *s == 'r' ) {
        digits = atoi (s+1);
		mpz_ui_pow_ui (x, 10, digits);
		mpz_urandomm (o, mpz_util_rands, x);
		gmp_info_printf ("Random Number: %Zd\n", o);
		return 1;
    }
    if ( *s == 'P' || *s == 'p' ) {
        digits = atoi (s+1);
		mpz_ui_pow_ui (x, 10, digits);
		mpz_urandomm (o, mpz_util_rands, x);
		do {
			mpz_nextprime (o, o);
		} while ( ! mpz_probab_prime_p (o, 20) );			// We really want to be sure here so we don't screw up test cases
		gmp_info_printf ("Random Prime: %Zd\n", o);
		return 1;
    }
    if ( *s == 'B' || *s == 'b' ) {
        digits = atoi (s+1);
		mpz_ui_pow_ui (x, 10, digits);
		mpz_urandomm (o, mpz_util_rands, x);
		do {
			mpz_nextprime (o, o);
		} while ( ! mpz_congruent_ui_p (o, 3, 4) || ! mpz_probab_prime_p (o, 20) );			// We really want to be sure here so we don't screw up test cases
		gmp_info_printf ("Random Prime 3mod4: %Zd\n", o);
		return 1;
    }
    if ( *s == 'C' || *s == 'c' ) {
        digits = atoi (s+1);
		mpz_ui_pow_ui (x, 10, digits/2);
		mpz_urandomm (y, mpz_util_rands, x);
		mpz_nextprime (y, y);
		mpz_urandomm (z, mpz_util_rands, x);
		mpz_nextprime (z, z);
		mpz_mul (o, y, z);
		gmp_info_printf ("Random Composite of 2 Primes: %Zd\n", o);
		return 1;
    }
   mpz_set_ui (o, 1);
    mpz_set_ui (x, 1);
    op = '*';
    if ( *s=='-' ) s++;
    while ( expr_valid_op(op) ) {
	    t = s;
	    while ( isdigit (*s) ) s++;
	    if ( s == t ) return 0;
	    nextop = *s;
	    *s++ = '\0';
	    mpz_set (y, x);
	    if ( mpz_set_str (x, t, 0) != 0 ) return 0;
		if ( nextop == '!' || nextop == '#' || nextop == '$' ) {
			if ( nextop == '!' ) {
				n = atoi (t);
				mpz_fac_ui (x, n);
			} else if ( nextop == '#' ) {
				mpz_set_ui (p,1);
				mpz_set_ui (z,1);
				for(;;) {
					mpz_nextprime (p, p);
					if ( mpz_cmp (p, x) > 0 ) break;
					mpz_mul (z, z, p);
				}
				mpz_set (x, z); 
			} else if ( nextop == '$' ) {
				mpz_set_ui (p,1);
				mpz_set_ui (z,1);
				for ( i = 0 ; mpz_cmp_ui (x, i) > 0 ; i++ ) {
					mpz_nextprime (p, p);
					mpz_mul (z, z, p);
				}
				mpz_set (x, z); 
			}
			nextop = *s;
			*s++ = '\0';
		}
	    if ( op == 'E' || op == 'e' || op == '^' ) {
	    	n = atoi (t);
		  	if ( n > 0 ) {
		  		mpz_pow_ui (x, y, n-1);
		  		mpz_mul (o, o, x);
		  	}
	    } else if ( op == '*' ) {
	    	mpz_mul (o, o, x);
	    } else if ( op == '+' ) {
	    	mpz_add (o, o, x);
	    	return 1;
	    } else if ( op == '-' ) {
	    	mpz_sub (o, o, x);
	    	return 1;
	    }
	    op = nextop;
	}
	if ( expr[0]=='-' ) mpz_neg(o,o);
	return 1;
}


static unsigned long mpz_mulm_counter, mpz_powm_counter, mpz_powm_tiny_counter;
static clock_t mpz_counter_reset_time;

void mpz_mulm (mpz_t o, mpz_t a, mpz_t b, mpz_t m)
    { mpz_mul (o, a, b);  mpz_mod (o, o, m); mpz_mulm_counter++; }

void mpz_addm (mpz_t o, mpz_t a, mpz_t b, mpz_t m)
    { mpz_add (o, a, b);  mpz_mod (o, o, m); }

void mpz_subm (mpz_t o, mpz_t a, mpz_t b, mpz_t m)
    { mpz_sub (o, a, b);  mpz_mod (o, o, m); }

// assumes input already mod m
void mpz_negm (mpz_t a, mpz_t m)
	{ if ( mpz_sgn(a) ) mpz_sub (a, m, a); }

void mpz_subm_ui (mpz_t o, mpz_t a, unsigned long b, mpz_t m)
    { mpz_sub_ui (o, a, b);  mpz_mod (o, o, m); }
    
void mpz_set_i (mpz_t o, long n)
{
	if ( n < 0 ) { mpz_set_ui (o, (unsigned long)(-n));  mpz_neg (o, o); } else { mpz_set_ui (o, (unsigned long)n); }
}

void mpq_set_i (mpq_t o, long n)
{
	if ( n < 0 ) { mpq_set_ui (o, (unsigned long)(-n), 1UL);  mpq_neg (o, o); } else { mpq_set_ui (o, (unsigned long)n, 1UL); }
}

// Optimize small powers - 0, 1, and 2 - assume b < m
void mpz_powm_tiny (mpz_t o, mpz_t b, unsigned e, mpz_t m)
{
    switch (e) {
    case 0: mpz_set_ui (o, 1);  return;
    case 1: if ( o != b ) mpz_set (o, b);  return;
    case 2: mpz_mulm (o, b, b, m);  return;
    default: mpz_powm_ui (o, b, e, m); mpz_powm_tiny_counter++;
    }
    return;
}

void mpz_powm_big (mpz_t o, mpz_t b, mpz_t e, mpz_t m)
    { mpz_powm (o, b, e, m);  mpz_powm_counter++; }

void mpz_reset_counters ()
{
    mpz_mulm_counter = 0;
    mpz_powm_counter = 0;
    mpz_powm_tiny_counter = 0;
    mpz_counter_reset_time = clock();
}

void mpz_report_counters ()
{
    info_printf ("Counters: %ld mult, %ld tiny exp, %ld exp, %ld msec\n", mpz_mulm_counter, mpz_powm_tiny_counter, mpz_powm_counter,
    	    	 delta_msecs(mpz_counter_reset_time, clock()));
}


// standard 4-ary exponentiation (fixed 2-bit window)
// p must fit in 32-bits.  This is not checked
unsigned long ui_powm_ui (unsigned long a, unsigned long e, unsigned long p)
{
	register unsigned long c, m;
	unsigned long b[4];
	register int i, j;
	
	if ( ! e ) return 1;
	i = _asm_highbit(e);
	if ( i&1 ) i--;
	m = 3UL<<i;
	b[1] = ui_mod(a,p);
	b[2] = ui_mod(b[1]*b[1],p);
	b[3] = ui_mod(b[1]*b[2],p);
	c = b[(m&e)>>i];
	for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
		c = ui_mod(c*c,p);  c = ui_mod(c*c,p);
		j = (m&e)>>i;
		if ( j ) c = ui_mod(c*b[j],p);
	}
	return c;
}


/*
long i_sqrt_modprime_slow (long n, long p)
{
	static int init;
	static mpz_t a, P;
	
	if ( ! init ) { mpz_init(a); mpz_init(P); init = 1; }
	
	mpz_set_ui(P,p);
	if ( n < 0 ) { 
		mpz_set_ui(a,-n);  mpz_mod(a,a,P); if ( mpz_sgn(a) ) mpz_sub(a,P,a);
	} else {
		mpz_set_ui(a,n);  mpz_mod(a,a,P);
	}
	if ( ! mpz_sqrt_modprime (a,a,P) ) return -1;
	return (long)mpz_get_ui(a);
}
*/

int _mpzutil_nr35_tab[35] =    { 0, 0, 5, 5, 0, 7, 7, 5, 5, 0, 7, 0, 5, 5, 0, 0, 0, 5, 5, 7, 7, 0, 5, 5, 7, 0, 7, 5, 5, 0, 0, 7, 5, 5, 7};

// Algorithm 1.5.1 in [CANT] - standard Tonelli-Shanks n^2 algorithm.  p must be an odd prime and this is not verified.
// Note that ff_sqrt is faster, but this code is useful if you only need one sqrt and p is smallish.  currently only supports p<2^32
 long i_sqrt_modprime (long a, long p)
{
	register long q, x, y, z, b;
	register int i, r, m;

	if ( p >= (1UL<<(ULONG_BITS/2)) ) { printf ("i_sqrt_modprime only supports p < 2^32\n"); abort(); }
	a = i_mod (a, p);
	if ( ! a ) return 0;
	if ( a==1 ) return 1;
	// Get a non-residue, use hardwired tests to catch all but 1/32 of the cases
	if ( (p&3)==3 ) {
		x = p-1;						// -1 is  a non-residue whenever p = 3 mod 4
	} else if ( (p&7)==5 ) {
		x = 2;						// 2 is a non-residue whenever p = 5 mod 8
	} else if ( (p%3)==2 ) {
		x = p-3;						// -3 is a non-residue whenever p = 2 mod 3
	} else {
		x = _mpzutil_nr35_tab[p%35];		// use 5 or 7 when possible
		if ( ! x ) for ( x = 11 ; ui_legendre(x,p) >= 0 ; x += 2 );
	}
	q = (p-1)>>1;
	for ( r = 1 ; !(q&1) ; r++ ) q >>= 1;
	y = ui_powm_ui (x,q,p);
	q = (q-1)>>1;
	z = ui_powm_ui (a,q,p);
	b = ui_mod(a*z,p);
	x= b;
	b = ui_mod(b*z,p);
	while ( b != 1 ) {
		q = b;
		for ( m = 1 ; m <= r; m++ ) {
			q = ui_mod(q*q,p);
			if ( q==1 ) break;
		}
		if ( m > r ) { printf ("Unexpected result: m=%d > r=%d in i_sqrt_modprime?!  a = %ld, p = %ld\n", m, r, a, p);  abort(); }
		if ( m == r ) return -1;
		q = y;
		for ( i = 0 ; i < r-m-1 ; i++ ) q = ui_mod(q*q,p);
		y = ui_mod(q*q,p);
		r = m;
		x = ui_mod(x*q,p);
		b = ui_mod(b*y,p);
	}
	if ( ui_mod(x*x,p) != a ) { printf ("i_sqrt_modprime failed %ld^2 != %ld mod %ld\n", x, a, p); abort(); }
	return x;
}

// Algorithm 1.5.1 in [CANT] - standard Tonelli-Shanks algorithm
// there seems to be a bug in here for very large p
int mpz_sqrt_modprime (mpz_t o, mpz_t a, mpz_t p)
{
	static mpz_t q, x, y, b;
	static int init;
	int i, r, m;

	if ( ! init ) { mpz_init (q);  mpz_init (x);  mpz_init (y);  mpz_init (b);  init = 1; }
	if ( ! mpz_sgn (a) ) { mpz_set_ui (o, 0);  return 1; }
	if ( mpz_cmp_ui (a, 1) == 0 ) { mpz_set_ui (o, 1);  return 1; }
	if ( mpz_divisible_p(a,p) ) { mpz_set_ui(o,0); return 1; }
	mpz_sub_ui (q, p, 1);
	mpz_set_ui (x, 2);
	r = mpz_remove (q, q, x);
	mpz_set_ui(x,3);
	do {
		mpz_add_ui(x,x,2);
	} while ( mpz_jacobi (x, p) != -1 );
	mpz_powm (y, x, q, p);
	mpz_sub_ui (q, q, 1);
	mpz_divexact_ui (q, q, 2);
	mpz_powm (x, a, q, p);
	mpz_mulm (b, a, x, p);
	mpz_mulm (b, b, x, p);
	mpz_mulm (x, a, x, p);
	for (;;) {
		if ( mpz_cmp_ui (b, 1) == 0 ) { mpz_set (o, x);  break; }
		mpz_set (q, b);
		for ( m = 1 ; m <= r; m++ ) {
			mpz_mulm (q, q, q, p);
			if ( mpz_cmp_ui (q, 1) == 0 ) break;
		}
		if ( m > r ) { gmp_printf ("Unexpected result: m=%d > r=%d in mpz_sqrt_modprime?!  a = %Zd, p = %Zd\n", m, r, a, p);  return 0; }
		if ( m == r ) return 0;
		mpz_set (q, y);
		for ( i = 0 ; i < r-m-1 ; i++ ) mpz_mulm (q, q, q, p);
		mpz_mulm (y, q, q, p);
		r = m;
		mpz_mulm (x, x, q, p);
		mpz_mulm (b, b, y, p);
	}
///	dbg_printf ("%Zd is a square root of %Zd mod %Zd\n", o, a, p);
	return 1;
}

// these are horribly inefficient
unsigned long mpz_get_bits_ui (mpz_t a, unsigned long i, unsigned long n)		// gets bits i thru i+n-1 and returns as ui
{
	register unsigned long j, k, x;

/*	y = 0;
	for ( j = i ; j < i+n ; j++ ) {
		if ( mpz_tstbit (a, j) ) y |= (1<<(j-i));
	}
	return y;
*/
	if ( mp_bits_per_limb != ULONG_BITS ) { err_printf ("mpz_get_bits_ui needs mp_bits_per_limb == ULONG_BITS!\n");  abort(); }
	j = i  >> ULONG_BITS_LOG2;
	if ( j >= a->_mp_size ) return 0;
	x = a->_mp_d[j];
	k = i-(j<<ULONG_BITS_LOG2);		// k = i mod ULONG_BITS < ULONG_BITS
	x >>= k;
	k = ULONG_BITS-k;				// k = # bits in x > 0
	if ( k < n ) {
		if ( ++j >= a->_mp_size ) return x;
		x |= (a->_mp_d[j]&((1UL<<(n-k))-1))<<k;	// get n-k more bits
	} else if ( k > n ) {
		x &= (1UL<<n)-1;
	}

//	if ( x != y ) { err_printf ("bug in mpz_get_bits_ui(%Zx(hex),i,n) = %d != %d\n", a, i, n, x, y);  abort(); }

	return x;
}

// Note bits is assumed to contain only the n bits to be set, the higher order bits must be zero
unsigned long mpz_set_bits_ui (mpz_t a, unsigned long i, unsigned long n, unsigned long bits)		// sets bits i thru i+n-1 and returns as ui
{
	register unsigned long j, k, k2, x;
	
/*
	for ( j = i ; j < i+n ; j++ ) {
		if ( bits&(1<<(j-i)) ) {
			mpz_setbit (a, j);
		} else {
			mpz_clrbit (a, j);
		}
	}
	return bits;
*/	
	if ( ! n ) return 0;
	if ( mp_bits_per_limb != ULONG_BITS ) { err_printf ("mpz_get_bits_ui needs mp_bits_per_limb == ULONG_BITS!\n");  abort(); }
	j = i  >> ULONG_BITS_LOG2;
	if ( j >= a->_mp_size ) { mpz_setbit (a, i+n-1); }	// force size adjustment and realloc if required
	k = i-(j<<ULONG_BITS_LOG2);		// k = i mod ULONG_BITS < ULONG_BITS
	x = bits<<k;
	x |= a->_mp_d[j] & ((1UL<<k)-1);		// get k low bits that should remain the same
	k2 = ULONG_BITS-k;				// k2 = # bits in set so far > 0
	if ( k2 < n ) {
		a->_mp_d[j] = x;
		if ( ++j >= a->_mp_size ) { mpz_setbit (a, i+n-1); }	// force size adjustment and realloc if required
		x = bits >>k2;
		x |= (a->_mp_d[j]&~((1UL<<(n-k2))-1));	// get high bits that should remain the same
		a->_mp_d[j] = x;
	} else if ( k2 > n ) {
		x |= a->_mp_d[j] & ~((1UL<<(k+n))-1); // get high bits past k+n bits that should remain the same
		a->_mp_d[j] = x;
	} else {
		a->_mp_d[j] = x;
	}
	while ( a->_mp_size && ! a->_mp_d[a->_mp_size-1] ) a->_mp_size--;
	return bits;
}

void mpz_randomm (mpz_t o, mpz_t m)
{
	mpz_util_init();
	mpz_urandomm (o, mpz_util_rands, m);
}

void mpz_randomb (mpz_t o, int b)
{
	mpz_util_init();
	mpz_urandomb (o, mpz_util_rands, b);
}

unsigned long ui_randomm (unsigned long m)
{
	mpz_util_init(); 
	return gmp_urandomm_ui (mpz_util_rands, m);
}

unsigned long ui_randomb (unsigned long b)
{
	mpz_util_init(); 
	return gmp_urandomb_ui (mpz_util_rands, b);
}

unsigned long ui_pp_div (unsigned long n, unsigned long p)
{
	unsigned long q;
	
	if ( ! n ) return 1;
	q = 1;
	while ( (n%p) == 0 ) {
		n /= p;
		q *= p;
	}
	return q;
}


// return 1 if the integer represented by the prime factorization (p,h,w) <ULONG_MAX contains a divisor in [min,max]
// slow implementation
int ui_divisor_in_interval (unsigned long p[], unsigned long h[], int w, unsigned long min, unsigned long max)
{
	register unsigned long n;
	unsigned long e[MPZ_MAX_FACTORS+1];
	register int i, j;
	
	if ( min > max ) return 0;
	if ( min <= 1 ) return 1;
	for ( i = 0 ; i < w ; i++ ) e[i] = 0;
	for (;;) {
		n = 1;
		for ( j = 0 ; j < w ; j++ ) {
			if ( e[j] ) {
				n *= p[j];
				for ( i = 1 ; i < e[j] ; i++ ) n *= p[j];
			}
		}
//		printf ("n=%d\n", n);
		if ( n >= min && n <= max ) return 1;
		e[0]++;
		for ( i = 0 ; e[i] > h[i] ; i++ ) { if ( i == w-1 ) return 0;  e[i] = 0;  e[i+1]++; }
	}
}

// a is compatible with b if it is composed entirely of primes dividing b
int ui_compatible  (unsigned long a, unsigned long b)
{
	unsigned long d;
	
	while ( a != 1 ) {
		d = ui_gcd(a,b);
		if ( d == 1 ) return 0;
		a /= d;
	}
	return 1;	
}

// a is compatible with b if it is composed entirely of primes dividing b
int mpz_compatible  (mpz_t a, mpz_t b)
{
	static int init;
	static mpz_t t, d;
	
	if ( ! init ) { mpz_init(d);  mpz_init(t);  init = 1; }
	mpz_set (t, a);
	while ( mpz_cmp_ui (t,1) != 0 ) {
		mpz_gcd (d, t, b);
		if ( mpz_cmp_ui(d,1) == 0 ) return 0;
		mpz_divexact(t,t,d);
	}
	return 1;	
}

unsigned long ui_crt (unsigned long a, unsigned long M, unsigned long b, unsigned long N)
{
	unsigned long x;
	
	if ( a > b ) {
		x = ui_inverse (N, M);
		return (x*N*(a-b)+b)%(M*N);
	} else {
		x = ui_inverse (M, N);
		return (x*M*(b-a)+a)%(M*N);
	}
}

unsigned long ui_gcd_ext (unsigned long a, unsigned long b, long *x, long *y)
{
	register unsigned long q, r, r0, r1;
	register long s, t, s0, s1, t0, t1;
	
	if ( a < b ) return (ui_gcd_ext (b, a, y, x));
	if ( b == 0 ) {
		if ( x ) *x = 1;
		if ( y ) *y = 0;
		return a;
	}
	if ( x ) { s0 = 1;  s1 = 0; }
	if ( y ) { t1 = 1;  t0 = 0; }
	r0 = a;
	r1 = b;
	while ( r1 > 0 ) {
		q = r0/r1;
		r = r0 - q*r1;
		r0 = r1;
		r1 = r;
		if ( y ) {
			t = t0 - q*t1;
			t0 = t1;
			t1 = t;
		}
		if ( x ) {
			s = s0 - q*s1;
			s0 = s1;
			s1 = s;
		}
	}
	if ( x ) *x = s0;
	if ( y ) *y = t0;
	return r0;
}

// Assumes gcd(a,b,N) = 1 but does not verify this
unsigned long bach_gcd (long a, long b, unsigned long N)
{
	unsigned long c, g, h;
	unsigned long M;
	
	g = ui_gcd(abs(a),N);
	if ( g == 1 ) return 0;
	if ( ui_gcd(abs(a+b),N) == 1 ) return 1;
	M = N;
	for ( h = g ; h > 1 ; ) {
		do { M /= h; } while ( (M%h) == 0 );
		h = ui_gcd (M, g);			// could use h instead of g?
	}
	// M>1 should now be a factor of N that is relatively prime to both g and N
	if ( M <= 1 || (N%M)!=0 ) { out_printf ("Error M=%ld in bach_gcd\n", M);  abort(); }
	c = (M*ui_inverse(M,N/M)) % N;
	if ( c == 0 ) { out_printf ("Error c == 0 in bach_gcd\n");  abort(); }
	return c;
}




int mpz_eval_term_ui (mpz_t o, unsigned long numvars, unsigned long vars[], unsigned long exps[])
{
	static int init;
	static mpz_t x;
	int i;
	
	if ( ! init ) { mpz_init (x);  init = 1; }
	
	mpz_set_ui (o, 1);
	for ( i = 0 ; i < numvars ; i++ ) {
		if ( exps[i] > 0 ) {
			mpz_ui_pow_ui (x, vars[i], exps[i]);
			mpz_mul (o, o, x);
		}
	}
	return 1;
}


char *ui_term_to_string (char *buf, unsigned long numvars, unsigned long vars[], unsigned long exps[])
{
	char *s;
	int i;
	
	if ( ! numvars ) { strcpy (buf, "1");  return buf; }
	s = buf;
	*s = '\0';
	for ( i = 0 ; i < numvars ; i++ ) {
		if ( exps[i] > 0 ) {
			if ( s > buf ) *s++ = '*';
			sprintf (s, "%lu^%lu", vars[i], exps[i]);
			while ( *s ) s++;
		}
	}
	return buf;
}

unsigned long ui_wt (unsigned long x)
{
	unsigned long i, n;
	
	n = 0;
	for ( i = 1 ; i <= x ; i <<= 1 ) {
		if ( i&x ) n++;
	}
	return n;
}

unsigned long ui_lg_ceil (unsigned long x)
{
	unsigned long i;

	i = _asm_highbit(x);
	if ( x==(1UL<<i) ) return i; else return i+1;
}

unsigned long ui_binexp_cost (unsigned long x)
	{ return ui_len(x) + ui_wt(x) - 2; }

unsigned long ui_get_bits (unsigned long x, unsigned long pos, unsigned long bits)
{
	unsigned long i, m;
	
	m = 0;
	for ( i = pos ; i < pos+bits ; i++ ) m |= (1<<i);
	m &= x;
	return (m >> pos);
}

char *ui_bit_string (char *buf, unsigned long x)
{
	int i;
	
	for ( i = ui_len(x)-1 ; i >= 0 ; i-- ) *buf++ = ( (x & (1<<i)) ? 1 : 0 );
	*buf = '\0';
	return buf;
}

void ui_invert_permutation (unsigned long p2[], unsigned long p1[], unsigned long n)
{
	unsigned long i;
	
	for ( i = 0 ; i < n ; i++ ) p2[p1[i]] = i;
}


double mpz_log2 (mpz_t a)
{
	unsigned long b, n;
	
	b = mpz_sizeinbase(a,2);
	if ( ! b ) return 0.0;			// define lg(0) = 0 for convenience
	if ( b < ULONG_BITS ) {
		n = mpz_get_ui (a);
		return log2((double)n);
	} else {
		n = mpz_get_bits_ui (a, b-(ULONG_BITS-1), ULONG_BITS-2);
		return (double)(b-1) + log2(1.0+(double)n/(double)(1UL<<(ULONG_BITS-2)));
	}
}


int _qsort_mpz_cmp (const void *a, const void *b)
{
	return mpz_cmp (*((mpz_t *)a), *((mpz_t *)b));
}

int ui_qsort_cmp (const void *a, const void *b)
{
	if ( *((unsigned long *)a) < *((unsigned long*)b) ) return -1;
	if ( *((unsigned long *)a) > *((unsigned long*)b) ) return 1;
	return 0;
}

int dbl_qsort_cmp (const void *a, const void *b)
{
	if ( *((double *)a) < *((double*)b) ) return -1;
	if ( *((double *)a) > *((double*)b) ) return 1;
	return 0;
}


void mpz_sort (mpz_t a[], unsigned long n)
{
	qsort (a, n, sizeof(a[0]), _qsort_mpz_cmp);
}


unsigned long mpz_randomf (mpz_t o, mpz_t N, mpz_t factors[], unsigned long w)
{
	return mpz_randomft (o, N, factors, NULL, w);
}


/*
	The algorithm below is an implementation of Bach's algorithm for computing random factored integers.
	See "Analytic Methods in the Analysis and Design of Number-Theoretic Algorithms", Eric Bach, MIT Press 1984
	Returns the factorization of a random integer in the range (N/2,N]
*/

unsigned long mpz_randomft (mpz_t o, mpz_t N, mpz_t factors[], mpz_ftree_t factor_trees[], unsigned long w)
{
	mpz_t d, Q, M;
	mpf_t x;
	mpz_ftree_t t;
	double b, u;
	unsigned long n;
	int i, j, k, m, q;
	
	mpz_init (d);
	mpz_init (M);
	mpz_util_init ();

	if ( w == 0 ) { err_printf ("Factor array too small - increase size!\n");  abort(); }

	// Base case
	if ( mpz_cmp_ui (N, MPZ_MAX_SMALL_INTEGER) <= 0 ) {
		// pick random r in (N/2,N]
		mpz_fdiv_q_2exp (M, N, 1);
		mpz_sub (d, N, M);
		mpz_urandomm (o, mpz_util_rands, d);
		mpz_add (o, o, M);
		mpz_add_ui (o, o, 1);
		// factor r using small factors array
		n = mpz_get_ui (o);
		k = 0;
		while ( n > 1 ) {
			if ( k >= w ) { err_printf ("Factor array too small - increase size!\n");  abort(); }
			q = mpz_small_factors[n];
			if ( ! q ) q = n;
			if ( n%q ) { err_printf ("Small factor array assert failed - %lu not divisible by %d\n", n, q);  abort(); }
			n /= q;
			for ( i = 0 ; i < k ; i++ ) {
				if ( mpz_divisible_ui_p (factors[i], q) ) {
					mpz_mul_ui (factors[i], factors[i], q);
					break;
				}
			}
			if ( i == k ) mpz_set_ui (factors[k++], q);
		}
		mpz_clear (d);
		mpz_clear (M);
		return k;
	}
	
	mpz_init (Q);
	mpf_init (x);
	t  = 0;
	for (;;) {
		if ( factor_trees ) {
			t = mpz_randomfpp (Q, N);
		} else {
			mpz_randompp (Q, N);
		}
		mpz_fdiv_q (M, N, Q);
		k = mpz_randomf (o, M, factors, w-1);
		for ( i = 0 ; i < k ; i++ ) {
			if ( mpz_divisible_p (factors[i], Q) ) {
				mpz_mul (factors[i], factors[i], Q);
				break;
			}
		}
		if ( i == k ) {
			if ( factor_trees ) factor_trees[k] = t;
			mpz_set (factors[k++], Q);
		} else {
			if ( factor_trees ) mpz_clear_ftree (t);
		}
		mpz_set_ui (o, 1);
		for ( i = 0 ; i < k ; i++ ) mpz_mul (o, o, factors[i]);
		b = (mpz_log2 (N) - 1.0) / mpz_log2 (o);
		mpf_urandomb (x, mpz_util_rands, 64);
		u = mpf_get_d (x);
		if ( u <= b ) break;
dbg_printf ("random reject %f > %f\n", u, b);
		if ( factor_trees ) for ( i = 0 ; i < k ; i++ ) mpz_clear_ftree(t);
	}
	// Need to clean-up factorizations of p-1 to make sure we deal with cases where we generated a prime power with
	// the factorization of p^k-1.  Fortunately (p-1)|p^k-1, so all the prime factors of p-1 are present in the factorization
	// of p^k-1, we just need to pick them out.
	if ( factor_trees ) {
dbg_printf ("cleaning up\n");
		for ( i = 0 ; i < k ; i++ ) {
			if ( ! mpz_perfect_power_p (factors[i]) ) continue;
			if ( ! mpz_pp_base (Q, factors[i]) ) { err_printf ("Unable to extract prime power base from a known prime power!");  abort(); }
			mpz_sub_ui (Q, Q, 1);
			m = 0;
			t = factor_trees[i];
			for ( j = 0 ; j < t->w ; j++ ) {
				mpz_gcd (d, Q, t->factors[j]);
				if ( mpz_cmp_ui (d,1) > 0 ) mpz_set (t->factors[m++], d);
			}
			for ( j = m ; j < t->w ; j++ ) mpz_clear (t->factors[j]);
			t->w = m;
			// sanity check the results
			mpz_set_ui (d, 1);
			for ( j = 0 ; j < t->w ; j++ ) mpz_mul (d, d, t->factors[j]);
			if ( mpz_cmp (Q, d) != 0 ) { gmp_err_printf ("Invalid factorization of p-1 after prime power cleanup!  %Zd != %Zd\n", d, Q);  abort(); }
		}
	}
	mpz_clear (d);
	mpz_clear (M);
	mpz_clear (Q);
	mpf_clear (x);
	return k;
}

// Returns a random prime power in the range [2,N]
void mpz_randompp (mpz_t Q, mpz_t N)
{
	static int init;
	static mpz_t M, M2, J, p, d, r;
	static mpf_t x, y;
	double b, u;
	unsigned long j, n;

	if ( ! init ) {
		mpz_init (M);  mpz_init (M2); mpz_init (J);  mpz_init (p);  mpz_init (d);  mpz_init (r);  mpf_init (x);  mpf_init (y); 
		mpz_util_init ();
		init = 1;
	}
	for (;;) {
		// select random j in [1,log2(N)]
		n = (unsigned long)ceil (mpz_log2 (N));
		j = gmp_urandomm_ui (mpz_util_rands, n) + 1;
		mpz_ui_pow_ui (J, 2, j);
		mpz_urandomm (r, mpz_util_rands, J);
		mpz_add (Q, J, r);
		if ( mpz_cmp (Q, N) > 0 ) continue;
		if ( mpz_perfect_power_p (Q) ) {
			if ( ! mpz_pp_base (p, Q) ) continue;
		} else {
			mpz_set (p, Q);
		}
		mpz_gcd (d, p, mpz_util_primorial);
		if ( mpz_cmp_ui (d,1) != 0 && mpz_cmp (d,p) != 0 ) continue;
		if ( ! mpz_probab_prime_p (p, 5) ) continue;
		mpf_urandomb (x, mpz_util_rands, 64);
		u = mpf_get_d (x);
		mpz_fdiv_q (M, N, Q);
		mpz_mul_ui (d, Q, 2);
		mpz_fdiv_q (M2, N, d);
		mpz_sub (M, M, M2);
//err_printf ("#(N/2q,N/q] = %Zd for N = %Zd, q = %Zd\n", M, N, Q);
		mpz_mul (M, M, J);
		mpf_set_z (x, M);
		mpf_set_z (y, N);
		mpf_div (x, x, y);
		b = mpf_get_d (x);
//err_printf ("b = %f, Q = %Zd, p = %Zd, log2(p) = %f, log2(N) = %f\n", b, Q, p, mpz_log2(p), mpz_log2(N));
		b *= mpz_log2(p) / mpz_log2(N);
		if ( u < b ) break;
//err_printf ("random reject %f >= %f\n", u, b);
	}
}

// optimize later
unsigned long ui_pp_base (unsigned long pp)
{
	static int init;
	static mpz_t x, b;
	
	if ( ! init ) { mpz_init (x);  mpz_init (b); }
	mpz_set_ui (x, pp);
	if ( mpz_probab_prime_p (x, 5) ) return mpz_get_ui (x);
	if ( ! mpz_pp_base (b, x) ) return 0;
	return mpz_get_ui (b);
}


/*
	Extracts base from a prime power and returns the exponent
*/
unsigned long mpz_pp_base (mpz_t b, mpz_t q)
{
	static int init;
	static mpz_t p, x, y, z;
	int c, n;
	
	if ( ! init ) { mpz_util_init ();  mpz_init (p);  mpz_init (x);  mpz_init (y);  mpz_init (z);  init = 1; }
	mpz_gcd (p, q, mpz_util_primorial);
	if ( mpz_cmp_ui (p, 1) != 0 ) {
		if ( mpz_cmp_ui(p,MPZ_MAX_SMALL_PRIME) > 0 || ! ui_is_small_prime (mpz_get_ui(p)) ) return 0;		// q divisible by more than 1 small prime
		mpz_fdiv_q (x, q, p);
		n = 1;
		while ( mpz_cmp_ui (x,1) > 0 ) {
			mpz_fdiv_qr (x, z, x, p);
			if ( mpz_cmp_ui (z, 0) != 0 ) return 0;		// q is not a prime power
			n++;
		}
		mpz_set (b, p);
		return n;
	} else {
		// This is not particular efficient, but it rarely gets used.
		if ( mpz_probab_prime_p (q, 5) ) { mpz_set (b, q);  return 1; }
		n = mpz_sizeinbase (q,2)/ui_lg_floor(mpz_util_primes[MPZ_GCD_PRIMES]);
		while ( n >= 0 ) {
			mpz_ui_pow_ui (z, 2, (mpz_sizeinbase (q,2) / n) + 1);		// upper bound
			mpz_ui_pow_ui (y, 2, (mpz_sizeinbase (q,2)-1) / n);			// lower bound
			do {
				mpz_add (x, y, z);
				mpz_tdiv_q_ui (x, x, 2);									// middle
				mpz_pow_ui (p, x, n);
				c = mpz_cmp (p, q);
				if ( ! c ) { if ( ! mpz_probab_prime_p (x, 5) ) return 0;  mpz_set (b, x);  return n; }
				if ( c < 0 ) {
					mpz_add_ui (y, x, 1);
				} else {
					mpz_sub_ui (z, x, 1);
				}
			} while ( mpz_cmp (y,z) <= 0 );
			n--;
		}
	}
	return 0;
}

mpz_ftree_t mpz_randomfp (mpz_t o, mpz_t N)
{
	static int init;
	static mpz_t N2;
	mpz_ftree_t t;
	int i;
	
	if ( ! init ) { mpz_init(N2);  init = 1; }
	mpz_sub_ui (N2, N, 1);
	mpz_fdiv_q_ui (N2, N2, 2);
	t = (struct factor_tree_item_struct *)mem_alloc (sizeof (*t));
	for ( i = 0 ; i < MPZ_MAX_TREE_FACTORS ; i++ ) mpz_init (t->factors[i]);
	for (;;) {
		t->w = mpz_randomf (o, N2, t->factors, MPZ_MAX_TREE_FACTORS);
		mpz_mul_ui (o, o, 2);
		mpz_add_ui (o, o, 1);
		if ( mpz_probab_prime_p (o, 5) ) break;
		out_printf (".");
	}
	mpz_sort (t->factors, t->w);
	for ( i = 0 ; i < t->w ; i++ ) {
		if ( ! mpz_tstbit (t->factors[i], 0) ) break;
	}
	if ( i < t->w  ) {
		mpz_mul_ui (t->factors[i], t->factors[i], 2);
	} else {
		if ( t->w >= MPZ_MAX_TREE_FACTORS ) { err_printf ("Exceeded MPZ_MAX_TREE_FACTORS!");  abort(); }
		for ( i = t->w ; i > 0 ; i-- ) mpz_set (t->factors[i], t->factors[i-1]);
		mpz_set_ui (t->factors[0], 2);
		t->w++;
	}
	return t;	
}


mpz_ftree_t mpz_randomfpp (mpz_t o, mpz_t N)
{
	static int init;
	static mpz_t d, p, N2;
	mpz_ftree_t t;
	int i;
	
	if ( ! init ) { mpz_init (d);  mpz_init (p);  mpz_init(N2);  init = 1; }
	mpz_sub_ui (N2, N, 1);
	mpz_fdiv_q_ui (N2, N2, 2);
	t = (struct factor_tree_item_struct*)mem_alloc (sizeof (*t));
	for ( i = 0 ; i < MPZ_MAX_TREE_FACTORS ; i++ ) mpz_init (t->factors[i]);
	for (;;) {
		t->w = mpz_randomf (o, N2, t->factors, MPZ_MAX_TREE_FACTORS);
		out_printf (".");
		mpz_mul_ui (o, o, 2);
		mpz_add_ui (o, o, 1);
		if ( mpz_perfect_power_p (o) ) {
			if ( ! mpz_pp_base (p, o) ) continue;
		} else {
			mpz_set (p, o);
		}
		mpz_gcd (d, p, mpz_util_primorial);
		if ( mpz_cmp_ui (d,1) != 0 && mpz_cmp (d,p) != 0 ) continue;
		if ( ! mpz_probab_prime_p (p, 5) ) continue;
	}
	mpz_sort (t->factors, t->w);
	for ( i = 0 ; i < t->w ; i++ ) {
		if ( ! mpz_tstbit (t->factors[i], 0) ) break;
	}
	if ( i < t->w  ) {
		mpz_mul_ui (t->factors[i], t->factors[i], 2);
	} else {
		if ( t->w >= MPZ_MAX_TREE_FACTORS ) { err_printf ("Exceeded MPZ_MAX_TREE_FACTORS!");  abort(); }
		for ( i = t->w ; i > 0 ; i-- ) mpz_set (t->factors[i], t->factors[i-1]);
		mpz_set_ui (t->factors[0], 2);
		t->w++;
	}
	return t;	
}


void mpz_clear_ftree (mpz_ftree_t t)
{
	int i;
	
	for ( i = 0 ; i < t->w ; i++ ) {
		if ( t->subtrees[i] ) mpz_clear_ftree (t->subtrees[i]);
	}
	for ( i = 0 ; i < MPZ_MAX_TREE_FACTORS ; i++ ) mpz_clear (t->factors[i]);
}


unsigned long mpz_store_bytes (mpz_t a)
{
	unsigned long size;
	
	size = _ui_ceil_ratio(mpz_sizeinbase(a,2),8);
	if ( size+2 >= (1<<15) ) { gmp_err_printf ("Integer %Zd too large to store!\n", a);  abort(); }
	return size+2;
}


void mpz_store (char *s, mpz_t a)
{
	size_t size, count;
	short cnt;
	
	size = _ui_ceil_ratio(mpz_sizeinbase(a,2),8);
	if ( size+1 >= (1<<15) ) { gmp_err_printf ("Integer %Zd too large to store!\n", a);  abort(); }
	mpz_export (s+sizeof(short), &count, 1, 1, 0, 0, a);
	if ( count >= (1<<15) ) { gmp_err_printf ("mpz_export wrote too many bytes for integer %Zd\n", a);  abort(); }
	cnt = (short)count;
	if ( mpz_sgn(a) < 0 ) cnt = -cnt;
	*((short*)s) = cnt;
	return;	
}


void mpz_retrieve (mpz_t o, char *s)
{
	short cnt;
	int sign;
	
	cnt = *((short*)s);
	sign = 1;
	if ( cnt < 0 ) { sign = -1;  cnt = -cnt; }
	mpz_import (o, (unsigned long)cnt, 1, 1, 0, 0, s+2);
	if ( sign < 0 ) mpz_neg (o, o);
	return;
}


unsigned long ui_eval_expr (char *expr)
{
	char *s;
	unsigned long n;
	int i, j, k;

	for ( s = expr ; *s && *s != 'e' ; s++ );
	if ( *s ) {
		k = atoi(expr);
		j = atoi(s+1);
		for ( i = 0, n = 1 ; i < j ; i++ ) n *= k;
	} else {
		n = atoll (expr);
	}
	for ( s = expr ; *s && *s != '+' && *s != '-' ; s++ );
	if ( *s == '+' ) return n+atol(s+1);
	if ( *s == '-' ) return n-atol(s+1);
	return n;
}

NAF_entry_t ui_NAF_table[2][4] = { {{0,0}, {1,0}, {0,0}, {-1,1}}, {{1,0},{0,1},{-1,1},{0,1}} };
unsigned char NAF_pbits_tab[1024] = {0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 64, 64, 64, 65, 66, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 72, 72, 72, 73, 74, 80, 80, 81, 80, 80, 80, 81, 82, 84, 84, 85, 128, 128, 128, 129, 130, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 136, 136, 136, 137, 138, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 128, 128, 128, 129, 130, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 136, 136, 136, 137, 138, 144, 144, 145, 144, 144, 144, 145, 146, 148, 148, 149, 160, 160, 160, 161, 162, 160, 160, 161, 160, 160, 160, 161, 162, 164, 164, 165, 168, 168, 168, 169, 170, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 64, 64, 64, 65, 66, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 72, 72, 72, 73, 74, 80, 80, 81, 80, 80, 80, 81, 82, 84, 84, 85, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 64, 64, 64, 65, 66, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 72, 72, 72, 73, 74, 80, 80, 81, 80, 80, 80, 81, 82, 84, 84, 85, 128, 128, 128, 129, 130, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 136, 136, 136, 137, 138, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 128, 128, 128, 129, 130, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 136, 136, 136, 137, 138, 144, 144, 145, 144, 144, 144, 145, 146, 148, 148, 149, 160, 160, 160, 161, 162, 160, 160, 161, 160, 160, 160, 161, 162, 164, 164, 165, 168, 168, 168, 169, 170, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 64, 64, 64, 65, 66, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 72, 72, 72, 73, 74, 80, 80, 81, 80, 80, 80, 81, 82, 84, 84, 85, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, };
unsigned char NAF_nbits_tab[1024] = {0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 85, 84, 84, 82, 81, 80, 80, 80, 81, 80, 80, 74, 73, 72, 72, 72, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 66, 65, 64, 64, 64, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 170, 169, 168, 168, 168, 165, 164, 164, 162, 161, 160, 160, 160, 161, 160, 160, 162, 161, 160, 160, 160, 149, 148, 148, 146, 145, 144, 144, 144, 145, 144, 144, 138, 137, 136, 136, 136, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 130, 129, 128, 128, 128, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 138, 137, 136, 136, 136, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 130, 129, 128, 128, 128, 85, 84, 84, 82, 81, 80, 80, 80, 81, 80, 80, 74, 73, 72, 72, 72, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 66, 65, 64, 64, 64, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 85, 84, 84, 82, 81, 80, 80, 80, 81, 80, 80, 74, 73, 72, 72, 72, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 66, 65, 64, 64, 64, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 170, 169, 168, 168, 168, 165, 164, 164, 162, 161, 160, 160, 160, 161, 160, 160, 162, 161, 160, 160, 160, 149, 148, 148, 146, 145, 144, 144, 144, 145, 144, 144, 138, 137, 136, 136, 136, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 130, 129, 128, 128, 128, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 138, 137, 136, 136, 136, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 130, 129, 128, 128, 128, 85, 84, 84, 82, 81, 80, 80, 80, 81, 80, 80, 74, 73, 72, 72, 72, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 66, 65, 64, 64, 64, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, };
unsigned char NAF_c_tab[1024] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, };
	
// n must be less than 2^63, takes about 14 nsecs for n~2^30 and about 32 nsecs for n~2^60 (AMD 2.5GHz Athlon64)
void ui_NAF (unsigned long *pbits, unsigned long *nbits, unsigned long n)
{
	register int c, m;
	
	m=n&0x1FF;
	*pbits = NAF_pbits_tab[m];
	*nbits = NAF_nbits_tab[m];
	if ( n < 128 ) return;
	c = NAF_c_tab[m];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 8;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 8;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 16;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 16;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 24;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 24;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 32;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 32;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 40;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 40;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 48;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 48;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 56;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 56;
	return;
}

	
unsigned long ui_next_prime (unsigned long p)
{
	static int init;
	static mpz_t P;
	
	if ( p < MPZ_MAX_SMALL_PRIME ) {
		mpz_util_init();
		return ui_small_prime(ui_small_prime_index(p)+1);
	}
	if ( ! init ) { mpz_init (P); init = 1; } 
	mpz_set_ui(P,p);
	mpz_nextprime(P,P);
	return mpz_get_ui(P);
}


// Given 0<d<p, finds a solution to x^2+dy^2=p, where p is an odd prime
// variant of Algorithm 1.5.2 in Cohen
int mpz_cornacchia (mpz_t x, mpz_t y, mpz_t d, mpz_t p)
{
	static mpz_t k, t, a, b, c, r, L;
	static int init;
	
	if ( ! init ) { mpz_init(k); mpz_init(t); mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(r); mpz_init(L); init = 1; }
	if ( mpz_sgn(d) <= 0 || mpz_cmp(d,p)>= 0 ) return 0;
	mpz_neg (t, d);
	if ( ! mpz_sqrt_modprime(k, t, p) ) return 0;
	mpz_mul_2exp(t,k,1);
	if ( mpz_cmp(t,p) < 0 ) mpz_sub(b,p,k);			// make sure p/2 < k < p
	else mpz_set(b,k);
	mpz_set(a,p);
	mpz_sqrt(L,p);
	while ( mpz_cmp(b,L) > 0 ) { mpz_mod(r,a,b);  mpz_set(a,b);  mpz_set(b,r); }
	mpz_mul(t,b,b);
	mpz_sub(r,p,t);
	if ( ! mpz_divisible_p(r,d) ) return 0;
	mpz_divexact(c,r,d);
	if ( ! mpz_perfect_square_p(c) ) return 0;
	mpz_sqrt(y,c);
	mpz_set(x,b);
	return 1;
}

// Given 0<d<4p and d congruent to 0 or 3 mod 4, finds a solution to x^2+dy^2=4p, where p is an odd prime
// Algorithm 1.5.3 in Cohen
int mpz_cornacchia4 (mpz_t x, mpz_t y, mpz_t d, mpz_t p)
{
	static mpz_t x0, D, a, b, c, r, L, p4;
	static int init;
	int dm4;

	if ( ! init ) { mpz_init(x0); mpz_init(D); mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(r); mpz_init(L); mpz_init(p4); init = 1; }
	mpz_mul_2exp (p4,p,2); 
	if ( mpz_sgn(d) <= 0 || mpz_cmp(d,p4)>= 0 ) return 0;
	dm4 = mpz_tstbit(d,0) + 2*mpz_tstbit(d,1);
	if ( dm4 && dm4 != 3 ) return 0;
	mpz_sub (D, p, d); while ( mpz_sgn(D) < 0 ) mpz_add(D,D,p);	// D = -d mod p
	if ( ! mpz_sqrt_modprime(x0, D, p) ) return 0;				// x0 = sqrt(D) mod p
	if ( mpz_tstbit(x0,0) ^ mpz_tstbit(d,0) ) mpz_sub(x0,p,x0);		// make sure parity of x0 matches parity of d (=parity of -d)
	mpz_mul_2exp (a,p,1);  mpz_set(b,x0);
	mpz_sqrt(L,p4);
	while ( mpz_cmp(b,L) > 0 ) { mpz_mod(r,a,b);  mpz_set(a,b);  mpz_set(b,r); }
	mpz_mul(c,b,b);
	mpz_sub(r,p4,c);
	if ( ! mpz_divisible_p(r,d) ) return 0;
	mpz_divexact(c,r,d);
	if ( ! mpz_perfect_square_p(c) ) return 0;
	mpz_sqrt(y,c);
	mpz_set(x,b);
	return 1;
}

/*
	Need to be careful of overflow here, it is quite possible that m3 (and/or k3) won't fit in a long, and we will set it to -1 in this case.
	Typically, the caller is only interested in elements of the sequence within a bounded range that does fit in a long, and one can
	easily have k3 small and m3 huge, yielding one element of interest.

	For safety/simplicity we resort to GMP to avoid any overflow issues
*/
int i_aseq_intersection (long *m3, long *k3, long m1, long k1, long m2, long k2)
{
	static mpz_t M, X, Y;
	static int init;
	long d;
	long u, v;
	
	if ( ! init ) { mpz_init(M); mpz_init(X); mpz_init(Y); init = 1; }
	if ( m1 < 0 ) m1 = -m1;
	if ( m2 < 0 ) m2 = -m2;
	d = ui_gcd_ext (m1, m2, &u, &v);
	if ( (k1-k2)%d ) return 0;
	m1 /= d;
	mpz_set_ui(M,m1); mpz_mul_ui(M,M,m2);
	m2 /= d;
	if ( mpz_cmp_ui(M,LONG_MAX) > 0 ) *m3 = -1; else *m3 = mpz_get_ui(M);
	// *k3 = (u*m1*k2 + v*m2*k1) % *m3;
	if ( u < 0 ) { mpz_set_ui(X,-u); mpz_neg(X,X); } else mpz_set_ui(X,u);
	if ( v < 0 ) { mpz_set_ui(Y,-v); mpz_neg(Y,Y); } else mpz_set_ui(Y,v);
	mpz_mul_ui(X,X,m1); mpz_mul_ui(X,X,k2); mpz_mul_ui(Y,Y,m2); mpz_mul_ui(Y,Y,k1);
	mpz_add(X,X,Y); mpz_mod(X,X,M);
	if ( mpz_cmp_ui(X,LONG_MAX) > 0 ) *k3 = -1; else *k3 = mpz_get_ui(X);
	return 1;
}

unsigned long ui_divisor_function (unsigned long n)
{
	unsigned long p[MPZ_MAX_UI_PP_FACTORS];
	unsigned long h[MPZ_MAX_UI_PP_FACTORS];
	unsigned long d;
	int i, w;
	
	if ( !n ) return 0;
	if ( n==1 ) return 1;
	w = ui_factor(p,h,n);
	for ( i = 0, d=1 ; i < w ; i++ ) d *= h[i]+1;
	return d;
}

int mpz_exp_window_parameters (int *window_size, int *max_digit, mpz_t e, int max_window_size)
{
	register int c, d, i, j, k, t;
	int mink, minc, maxd, minmaxd; 
	
	t = mpz_sizeinbase (e, 2)-1;
	minmaxd = 1;
	mink = 1;
	minc = t + mpz_popcount(e) - 1;							// binary method
	for ( k = 2 ; k <= t+1 && k <= max_window_size ; k++ ) {
		j = k;
		c = 0;
		d = mpz_get_bits_ui (e, t-j+1, j);
		while ( ! (d&0x1) && d != 2 ) { d >>= 1;	j--; }
		maxd = d;
		for ( i = t-j ; i >= 0 ; ) {
			for ( ; i >= 0 && ! mpz_tstbit (e,i) ; i-- ) c++;
			if ( i < 0 ) break;
			j = _ui_min (i+1,k);
			d = mpz_get_bits_ui (e, i-j+1, j);
			while ( ! (d&0x1) ) { d >>= 1;	j--; }
			if ( d > maxd ) maxd = d;
			i -= j;
			for ( ; j > 0 ; j-- ) c++;
			c++;
		}
		c += (maxd+1)/2;
		if ( c < minc ) {
			for ( mink = 0 ; (1<<mink) < maxd ; mink++ );
			minc = c;
			minmaxd = maxd;
		}
	}
	*window_size = mink;
	*max_digit = minmaxd;
	return minc;
}


#define MAX_MF		16

// computes positive solutions (x,y) to x^2-dy^2=m with d>m^2, m nonzero, and x < 2^h.  does NOT compute solution to u^2-dv^2=1
int mpz_pell_solver (mpz_t x[MPZ_MAX_PELL_SOLUTIONS], mpz_t y[MPZ_MAX_PELL_SOLUTIONS], long d, long m, int h)
{
	static int init;
	static mpz_t p[3], q[3], t0, t1;
	long P[3],Q[3],a[3], a0;
	long mf[MAX_MF], md[MAX_MF], am;
	register double z;
	register int nf, i, j, k, n, nm1, nm2, sm;

	if ( d <= m*m ) { printf ("Error in mpz_pell_solver, must have d > m^2, d=%ld, m=%ld\n", d, m); abort(); }
	if ( ! init ) { for ( k = 0 ; k < 3 ; k++ ) { mpz_init(p[k]);  mpz_init(q[k]); } mpz_init(t0); mpz_init(t1); init =1; }
	// make a list of square divisors of m and corresponding factors
	sm = ( m < 0 ? -1 : 1 );  am = sm*m;
	mf[0] = am; md[0] = 1;
	n = sqrt(am);
	for ( nf = 1, j = 2 ; j <= n ; j++ ) {
		if ( !(m%(j*j)) ) { if ( nf >= MAX_MF ) { printf ("Exceeded MAX_MF in mpz_pell_solver for m=%ld\n", m); abort(); } md[nf] = j;  mf[nf] = am/(j*j); nf++; }
	}
	P[0] = 0; Q[0] = 1; a[0] = a0 = (long)(z=sqrt(d));  mpz_set_ui(p[0],a0); mpz_set_ui(q[0],1);
	// check for early solution - fixes previous bug
	j = n = 0;
	mpz_mul(t0,p[n],p[n]);  mpz_mul(t1,q[n],q[n]);  mpz_mul_ui(t1,t1,d);  mpz_sub(t0,t0,t1);
	if ( mpz_sgn(t0)==sm ) {
		for ( i = 0 ; i < nf ; i++ ) {
			if ( mpz_cmpabs_ui(t0,mf[i])==0 ) {
				mpz_mul_ui(x[j],p[n],md[i]);  mpz_mul_ui(y[j],q[n],md[i]);
				if ( ++j == MPZ_MAX_PELL_SOLUTIONS ) { printf ("Hit MAX_SOLUTIONS=%d in mpz_pell_solver for d=%ld, m=%ld\n", j, d, m); return j; }
			}
		}
	}
	P[1] = a0; Q[1] = d-a0*a0; if ( ! Q[1] ) return 0;   							// allow for the case that d is actually square, but always return 0 in this situation
	a[1] = (z+P[1]) / Q[1];  mpz_set_ui(p[1],a0*a[1]+1);  mpz_set_ui(q[1],a[1]);  n = 1;
	for ( k = 2 ; ; k++ ) {
		if ( mpz_sgn(t0)==sm ) {
			for ( i = 0 ; i < nf ; i++ ) {
				if ( mpz_cmpabs_ui(t0,mf[i])==0 ) {
					mpz_mul_ui(x[j],p[n],md[i]);  mpz_mul_ui(y[j],q[n],md[i]);
					if ( ++j == MPZ_MAX_PELL_SOLUTIONS ) { printf ("Hit MAX_SOLUTIONS=%d in mpz_pell_solver for d=%ld, m=%ld\n", j, d, m); return j; }
				}
			}
		}
		n = k%3;  nm1=(k-1)%3;  nm2 = (k-2)%3;
		P[n] = a[nm1]*Q[nm1]-P[nm1];
		Q[n] = (d-P[n]*P[n])/Q[nm1];
		a[n] = (z+P[n]) / Q[n];
		mpz_mul_ui(t0,p[nm1],a[n]);  mpz_add(p[n],t0,p[nm2]);
		if ( h && mpz_sizeinbase(p[n],2) > h ) break;
		mpz_mul_ui(t0,q[nm1],a[n]);  mpz_add(q[n],t0,q[nm2]);
		mpz_mul(t0,p[n],p[n]);  mpz_mul(t1,q[n],q[n]);  mpz_mul_ui(t1,t1,d);  mpz_sub(t0,t0,t1);
	}
	return j;
}
