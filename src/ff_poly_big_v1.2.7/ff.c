#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include "gmp.h"
#include "ff.h"
#include "cstd.h"

/*
    Copyright 2007-2013 Andrew V. Sutherland
    See LICENSE file for license details.
*/

ff_t _ff_t1;
ff_t _ff_p;
ff_t _ff_2p;
ff_t _ff_2g;		// generator of Sylow 2-subgroup and necessarily a quadratic non-residue
ff_t _ff_2gi;		// inverse of _ff_2g
ff_t _ff_3g;		// generator of Sylow 3-subgroup
ff_t _ff_negone;
ff_t _ff_negthree;
ff_t _ff_half;
ff_t _ff_third;
ff_t _ff_fourth;
ff_t _ff_fifth;		// fractions >= 1/5 are set only when needed
ff_t _ff_seventh;
ff_t _ff_eleventh;
ff_t _ff_thirteenth;
ff_t _ff_frac[FF_FRACTIONS+1];						// _ff_frac[i] = 1/i for i > 0.  Only set if ff_invert_fractions is called.  These are all set at once, independent of the fixed fractions above.
static unsigned long mod256itab[256];
static ff_t __ff_mont_itab[FF_ITAB_SIZE];				// private copy for single modulus environment

static ff_t _ff_2exp;								// 2^(p+3)/4-1 set only for p=5mod8, used for fast sqrts
unsigned long _ff_p2_m;						// odd part of p-1, p=2^p2_e*p2_m+1
int _ff_p2_e;									// power of 2 dividing p-1
unsigned long _ff_p3_m;						// p=3^p3_e*p3_m+1 when p=1mod3, p=3^p3_3*p3_m-1 when p=2mod3
int _ff_p3_m1mod3;								// true if p3_m is 1 mod 3
int _ff_p3_e;									// power of 3 dividing  p-1 when p=1mod3, power of 3 dividing p+1 when p=2mod3
ff_t _ff_2Sylow_tab[64][2];						// 2Sylow_tab[i][0]=g^(2^i) for 0<=i<=e-1, where g generates the 2-Sylow subgroup.
											// 2Sylow_tab[i][1] = 2Sylow_tab[i][0]*2Sylow_tab[i+1][0]
static ff_t _ff_3Sylow_tab[42][2];					// 3Sylow_tab[i][0]=g^(3^i) for 0<=i<=e-1, where g generates the 3-Sylow subgroup (when it is non-trivial)
											// 3Sylow_tab[i][1] = 3Sylow_tab[i][0]^2
ff_t _ff_binomial_tab[((FF_BINOMIAL_MAX+1)*(FF_BINOMIAL_MAX+2))/2];
int _ff_binomial_top;

ff_t _ff_cbrt_unity;
int _ff_cbrt_setup;

int _ff_p1mod3;

int _ff_sqrt_chain_len, _ff_cbrt_chain_len, _ff_invcbrt_chain_len;
int _ff_sqrt_chain[FF_MAX_CHAIN_LEN], _ff_cbrt_chain[FF_MAX_CHAIN_LEN], _ff_invcbrt_chain[FF_MAX_CHAIN_LEN];

ff_t *_ff_mont_itab = __ff_mont_itab;				// my be repointed to mange multiple moduli
ff_t _ff_mont_R;
ff_t _ff_mont_R2;
unsigned long _ff_mont_pni;


void ff_montgomery_setup (int ring);

void ff_ext_setup(void);

// nr35_tab[p%35] is either 0 or a non-residue mod p for p=1mod4 (if p is 2 or 3 mod 5, then 5 is a non-residue mod p; if p is 3, 5, or 6 mod 7, then 7 is a non-residue mod p).
// nric35_tab[p%35] is, the coefficient k s.t. kp+1=0 mod ff_nr35_tab[p%35]
int _ff_nr35_tab[35] =    { 0, 0, 5, 5, 0, 7, 7, 5, 5, 0, 7, 0, 5, 5, 0, 0, 0, 5, 5, 7, 7, 0, 5, 5, 7, 0, 7, 5, 5, 0, 0, 7, 5, 5, 7};
int _ff_nric35_tab[35] = { 0, 0, 2, 3, 0, 4, 1, 2, 3, 0, 2, 0, 2, 3, 0, 0, 0, 2, 3, 4, 1, 0, 2, 3, 2, 0, 4, 2, 3, 0, 0, 2, 2, 3, 1};

gmp_randstate_t ff_randstate;
int ff_randstate_init;

	
// Note that it is possible to use non-prime p, and we don't want to waste time checking primality anyway.
// We assume that either the caller has checked the primality of p, or is aware that it is working over a ring
void ff_setup_ui (unsigned long p)
{
	static int init = 0;
	unsigned long n;
	
#if  FF_WORDS != 1
	err_printf ("ff_setup_ui only supports single-word Montgomery or native representation\n");
	abort();
#endif
#if FF_MONTGOMERY
	if ( ! init ) { for ( n = 1 ; n < 256 ; n+= 2 ) mod256itab[n] = ff_ui_inverse(n,256); init = 1; }		// This is only done once
#endif

	assert (ULONG_BITS >= 64);
	if ( p && _ff_p == p ) return;
	if ( p > (1UL<<FF_BITS) ) { printf ("ff_setup: p=%lu is too large, FF_BITS = %d\n", p, FF_BITS); abort(); }
	
	// Don't deal with char 2 (or any other even primes, for that matter...)
	if ( !(p&1) ) { printf ("ff_setup: p=%lu must be an odd prime\n", p);  abort(); }

	_ff_p = p;
	_ff_2p = p+p;
#if FF_MONTGOMERY
	ff_montgomery_setup (0);
#endif
	_ff_p2_m = (_ff_p-1)/2;
	_ff_set_ui (_ff_half, _ff_p2_m+1);  _ff_mult(_ff_fourth,_ff_half,_ff_half);
	_ff_neg(_ff_negone, _ff_mont_R);
	_ff_add(_ff_negthree,_ff_negone,_ff_negone);  _ff_addto(_ff_negthree,_ff_negone);
	for ( _ff_p2_e = 1 ; !(_ff_p2_m&0x1) ; _ff_p2_m>>=1, _ff_p2_e++ );
	_ff_cbrt_setup = 0;

	_ff_p1mod3 = ((p%3)==1);
	n = (_ff_p1mod3? (2*p+1UL)/3UL : (p==3?0:(p+1UL)/3UL));
	_ff_set_ui(_ff_third,n);				// note this is set to zero for p=3
	if ( _ff_p1mod3 ) {
		_ff_p3_m = (_ff_p-1)/3;
		for ( _ff_p3_e = 1 ; !(_ff_p3_m%3) ; _ff_p3_m /= 3, _ff_p3_e++ );
		_ff_p3_m1mod3 = ((_ff_p3_m%3)==1);
	} else {
		// this code is meaningless if p=3, but we won't bother checking this because these values won't be used in this case
		_ff_p3_m = (_ff_p+1)/3;
		for ( _ff_p3_e = 1 ; !(_ff_p3_m%3) ; _ff_p3_m /= 3, _ff_p3_e++ );
		_ff_p3_m1mod3 = ((_ff_p3_m%3)==1);
	}
	
	_ff_2exp = _ff_2g = _ff_3g = 0;								// will get set when needed
	_ff_fifth=_ff_seventh=_ff_eleventh=_ff_thirteenth=_ff_frac[1]=0;	// ditto
	_ff_binomial_top = 0;
#if FF_NO_EXT_SETUP != 1
	ff_ext_setup();
#endif
	return;
}

// macros to handle setting ff_t types in raw mode -- needed to handle multiples of p and montgomery values
#if FF_WORDS == 1
#define _ff_raw_set_mpz(z,X)		((z) = mpz_get_ui(X))
#define _ff_raw_get_mpz(Z,x)		(mpz_set_ui(Z,x),Z)
#endif


#if FF_WORDS == 1 && FF_MONTGOMERY

void ff_montgomery_setup (int ring)
{
	register unsigned long p, t;
	register int i;
	
	// precompute -p^{-1} mod B (B is either 2^32 or 2^64), R = B mod p and R^2 mod p
	// Use Jebelean's trick for computing p^{-1} mod B (see HECHECC p. 190 remark 10.4 (ii))
	p = _ff_p;
	t = mod256itab[p&0xFF];
	t = (2*t + (-(p*t*t)))&0xFFFF;
	t = (2*t + (-(p*t*t)))&0xFFFFFFFF;
#if FF_HALF_WORDS == 1
	_ff_mont_pni = (-t)&0xFFFFFFFF;
	t = (1UL<<FF_MONTGOMERY_RBITS)%p;
	_ff_mont_R = (ff_t)t;
	t = (t*t)%p;
	_ff_mont_R2 = (ff_t)t;
#else
	t = (2*t + (-(p*t*t)))&0xFFFFFFFFFFFFFFFF;
	_ff_mont_pni = -t;
	_ff_mont_R = (1UL<<(FF_MONTGOMERY_RBITS-1))%p;
	_ff_mont_R += _ff_mont_R;
	if ( _ff_mont_R > p ) _ff_mont_R -= p;
		
	// to avoid depending on multi-precision arithmetic, we square R mod p by doubling it RBITS times (R = 2^RBITS)
	// probably no slower than a modular reduction anyway
	_ff_mont_R2 = _ff_mont_R;
	for ( i = 0 ; i < FF_MONTGOMERY_RBITS ; i++ ) {
		_ff_mont_R2 += _ff_mont_R2;
		if ( _ff_mont_R2 > p ) _ff_mont_R2 -= p;
	}
#endif
	if ( ring ) return;
	t = 1;
	for ( i = 0 ; i <= 3*FF_MONTGOMERY_RBITS ; i++ ) {
		_ff_mont_itab[i] = t;
		t += t;													// note that _ff_p < 2^(ULONG_BITS-1) so overflow can't occur
		if ( t >= _ff_p ) t -= _ff_p;
	}
//	printf ("Montgomery setup p = %lu, -p^{-1} mod b = %lu, R mod p = %lu, R^2 mod p = %lu\n", _ff_p, _ff_mont_pni, _ff_mont_R, _ff_mont_R2);
}

// This is a streamlined version of Algorithm 11.12 in HECHECC p. 208
// The input is in Montgomery representation: x = [a] = aR mod p
// This functions computes v = [a^{-1}] = a^{-1}R mod p = x^{-1}R^2 mod p
unsigned long ff_montgomery1_invert (unsigned long x)
{
	register unsigned long r, s, t, v;
	register int k;

	if ( _ff_zero(x) ) { printf ("ff_invert: attempt to invert 0,\n"); abort();  return 0; }
	r = x;  s = 1;  t = _ff_p;  v = 0;  k = 0;
	while ( r > 0 ) {
		if ( !(t&0x1UL) ) {
			t >>= 1;  s <<= 1;
		} else if ( !(r&0x1UL) ) {
			r >>= 1;  v <<= 1;
		} else if ( t > r ) {
			t = (t-r)>>1;  v += s;  s <<= 1;
		} else {
			r = (r-t)>>1; s += v;  v <<= 1;
		}
		k++;
//		if ( k > 2*_ff_pbits ) { printf ("ff_invert: k = %d > 2*pbits!  raw input %lu\n", k, x);  abort(); }
	}
	if ( v >= _ff_p ) v -= _ff_p;				
	v = _ff_p - v;														// v cannot be zero if x is invertible (which we assume to be the case)
	// This accomplishes steps 10-12 of Alg 11.12 in one multiplication
	v = ff_montgomery1_mult (v, _ff_mont_itab[3*FF_MONTGOMERY_RBITS-k]);		// Montgomery lookup table holds powers 2^k to save time
#if ! FF_FAST
	if ( ff_montgomery1_mult (v,x) != _ff_mont_R ) {printf ("ff_montgomery1_invert failed, %lu*%lu != 1, raw input: %lu  raw output %lu\n", _ff_get_ui(x), _ff_get_ui(v), x, v); abort(); }
#endif
	return v;
}

#endif

// Algorithm 11.15 of [HECHECC], P. 209 - note that z and x can overlap
void ff_parallel_invert (ff_t z[], ff_t x[], unsigned n)
{
	ff_t c[FF_MAX_PARALLEL_INVERTS];
	register ff_t u, v;
	register unsigned i;

	if ( n > FF_MAX_PARALLEL_INVERTS ) {
		for ( i = 0 ; i+FF_MAX_PARALLEL_INVERTS < n ; i += FF_MAX_PARALLEL_INVERTS ) ff_parallel_invert(z+i,x+i,FF_MAX_PARALLEL_INVERTS);
		z += i;  x += i;  n-= i;
	}
	if ( ! n ) return;
	
	_ff_set (c[0], x[0]);
	for ( i = 1 ; i < n ; i++ ) _ff_mult (c[i], c[i-1], x[i]);
	 _ff_invert (u, c[n-1]);
	for ( i = n-1 ; i > 0 ; i-- ) {
		_ff_mult (v, c[i-1], u);
		_ff_mult (u, u, x[i]);
		_ff_set (z[i], v);
	}
	_ff_set (z[0], u);
}

int ff_parallel_invert_check (ff_t z[], ff_t x[], unsigned n)
{
	ff_t c[FF_MAX_PARALLEL_INVERTS];
	register ff_t u, v;
	register unsigned i;

	if ( n > FF_MAX_PARALLEL_INVERTS ) {
		for ( i = 0 ; i+FF_MAX_PARALLEL_INVERTS < n ; i += FF_MAX_PARALLEL_INVERTS ) ff_parallel_invert(z+i,x+i,FF_MAX_PARALLEL_INVERTS);
		z += i;  x += i;  n-= i;
	}
	if ( ! n ) return 1;
	
	_ff_set (c[0], x[0]);
	for ( i = 1 ; i < n ; i++ ) _ff_mult (c[i], c[i-1], x[i]);
	if ( _ff_zero(c[n-1]) ) return 0;
	 _ff_invert (u, c[n-1]);
	for ( i = n-1 ; i > 0 ; i-- ) {
		_ff_mult (v, c[i-1], u);
		_ff_mult (u, u, x[i]);
		_ff_set (z[i], v);
	}
	_ff_set (z[0], u);
	return 1;
}

static inline void ff_setup_2exp(void) { ff_t t; if ( _ff_2exp ) return;  _ff_set_one(t);  _ff_x2(t);  ff_exp_ui (&_ff_2exp,&t,(_ff_p+3)/4-1); }

#if FF_NO_EXT_SETUP!=1

void ff_setup_2g (void)
{
	ff_t t;
	register int i,n;
	
	if ( _ff_2g ) return;
	_ff_sqrt_chain_len = ff_precompute_exp_chain (_ff_sqrt_chain, (_ff_p2_m-1)/2);

	// When p=3 mod 4 life is very simple, since -1 then generates the 2-Sylow subgroup
	if ( _ff_p2_e == 1 ) { _ff_set(_ff_2g,_ff_negone);  _ff_set(_ff_2gi,_ff_negone); _ff_set(_ff_2Sylow_tab[0][0],_ff_negone);  return; }
	
	// hardwire a few quick tests for non-residues to speed things up.  This should catch all but 1/32 of the primes, on average.
	_ff_set_zero(t);
	if ( ! _ff_p1mod3 ) { _ff_set_i (t,-3); goto SylowSetup; }					// if p is 2 mod 3 then -3 is a non-residue
	if ( (_ff_p&0x7UL) == 5 ) { _ff_set_ui (t, 2); goto SylowSetup; }			// if p is 5 mod 8 then 2 is a non-residue
	n = _ff_nr35_tab[_ff_p%35];
	if ( n ) { _ff_set_ui (t, n); goto SylowSetup; }
	// We only reach this prob 1/32, and the expected cost of the computation below is less than 32 modular reductions
	// Note that legendre is faster here than computing square roots because n is small
	for ( n = 11 ; ; n++ ) if ( ui_legendre (n, _ff_p) < 0 ) break;
	_ff_set_ui(t,n);
	
SylowSetup:

	ff_exp_ui (&_ff_2g, &t, _ff_p2_m);
	_ff_set(_ff_2Sylow_tab[0][0], _ff_2g);
	_ff_set(_ff_2gi,_ff_2g);
	// compute the inverse of _ff_2g as we square up, this is faster than calling _ff_inverse
	for ( i = 1 ; i < _ff_p2_e-1 ; i++ ) { _ff_square(_ff_2Sylow_tab[i][0],_ff_2Sylow_tab[i-1][0]);  ff_mult(_ff_2gi,_ff_2gi,_ff_2Sylow_tab[i][0]); }
	_ff_set(_ff_2Sylow_tab[i][0],_ff_negone);
	ff_negate(_ff_2gi);
	for ( i = 0 ; i < _ff_p2_e-1 ; i++ ) _ff_mult(_ff_2Sylow_tab[i][1], _ff_2Sylow_tab[i][0], _ff_2Sylow_tab[i+1][0]);
	_ff_set(_ff_2Sylow_tab[i][1],_ff_negone);
}

/*
   Extended sqrt algorithm, based on Tonelli-Shanks, computes the square root in F_p^2 if necessary.
   The return value is 0 if the root is in F_p^2, in which case o=a^{-1/2}/z where F_p^2=F_p[z]/(z^2-_ff_2g).

   This algorithm is deterministic (given precomputed 2-Sylow generator) - caller should flip coin to randomize choice of sqrt.
*/
int ff_invsqrt (ff_t o[1], ff_t a[1], int ext)
{
	ff_t b, x, y;
	int sts;

	ff_setup_2g();
	if ( _ff_zero(a[0]) ) { _ff_set_zero (o[0]);  return 1; }			// need to handle zero separately (make sure setup_2g gets called though)

//	ff_exp_ui (&x, a, (_ff_p2_m-1)/2);			// x = a^{(m-1)/2}
	ff_exp_chain (&x, a, _ff_sqrt_chain, _ff_sqrt_chain_len);			// x = a^{(m-1)/2}
	ff_mult (b, a[0], x);
	ff_mult (b, b, x);						// b = a^m is in the 2-Sylow subgroup
	sts = ff_2Sylow_invsqrt(&y,&b, ext);
	if ( ! ext && ! sts ) return 0;
	ff_mult (x, x, y);						// x = a^{(m-1)/2} * b^{-1/2} = a^{(m-1)/2} * a^{-m/2} = a^{-1/2}
#if ! FF_FAST
	_ff_square (y, x);
	if ( ! sts ) ff_mult(y,y,_ff_2g);
	ff_mult(y,y,a[0]);
	if ( ! _ff_one(y) ) { 
		printf ("p=%lu, a=%lu, b=%lu, 2g=%lu, 2gi=%lu\n", _ff_p, _ff_get_ui(a[0]), _ff_get_ui(b), _ff_get_ui(_ff_2g), _ff_get_ui(_ff_2gi));
		if ( sts ) {
			printf ("ff_invsqrt failed, %lu^2 != 1/%lu\n", _ff_get_ui(x),  _ff_get_ui(a[0])); abort();
		} else {
			printf ("ff_sqrt failed, %lu^2*%lu != 1/%lu\n", _ff_get_ui(x), _ff_get_ui(_ff_2g), _ff_get_ui(a[0])); abort();
		}
	}
#endif
	_ff_set (o[0], x); 
	return sts;
}

/*
   Computes a^{-1/2} for a in the 2-Sylow subgroup, computing in F_p^2 if necessary.
   Returns 0 if the root is in F_p^2, in which case o=a^{-1/2}/z as above, otherwise returns 1

   Uses precomputed table of 2e powers of the 2-Sylow generator to reduce running time by a factor of 4 over standard Tonelli-Shanks (still O(e^2)).
   Assumes _ff_p2_e >= 4 and a has order at least 4 (easy cases handle inline by ff.h)
*/
int _ff_2Sylow_invsqrt (ff_t o[1], ff_t a[1], int ext)
{
	register ff_t b, q, t, w1, w2;
	register int j, k, s;
	int sts;

	if ( _ff_p2_e <= 2 ) { printf ("assert _ff_p2_e > 2 failed for p=%ld, _ff_p2_e=%d\n", _ff_p, _ff_p2_e); abort(); }
//	assert ( _ff_p2_e > 2 );

	// set w1 and w2 to the two elements of order 4 in the 2-Sylow subgroup (note we know e > 2 and a has order at least 4 since it is not 1 or -1)
	_ff_set (w1, _ff_2Sylow_tab[_ff_p2_e-2][0]);
	_ff_set (w2, _ff_2Sylow_tab[_ff_p2_e-2][1]);
	_ff_set_one (t);
	_ff_set (b,a[0]);
	sts = 1; j = 0;
	for(;;) {
		_ff_set (q, b);
		for ( s = 2 ; s <= _ff_p2_e ; s++ ) {						// s <= _ff_p2_e is an unnecessary safety check if a is actually in the 2-Sylow
			j=1;
			if ( _ff_equal(q,w1) ) break;
			j=0;
			if ( _ff_equal(q,w2) ) break;
			ff_square (q,q);
		}
		k = _ff_p2_e-s;
#if ! FF_FAST
		if ( k < 0 ) { printf ("Unexpected result: k<0 in ff_2Sylow_invsqrt?!  a = %lu, p = %lu\n", _ff_get_ui(a[0]), _ff_p); abort(); }
#endif	
		if ( !k ) {											// this can only happen the first time through
			if ( ! ext ) return 0;
			sts=0;
			if ( j ) {
				ff_mult(b,b,_ff_2gi);							// clear the low bit rather than adding
				ff_mult(t,t,_ff_2gi);							// setting sts=0 implicitly multiplied t by s^{1/2}, so multiplying by s^{-1} gives s^{-1/2}
			} else {
				ff_mult(b,b,_ff_2g);							// sts=0 multiplied t by s^{1/2}, so we don't need to do anything to t
			}
		} else {
			ff_mult (b, b, _ff_2Sylow_tab[k][j]);					// the product of all the elements S[k] we multiply into b here is b^{-1}, since we terminate with b=1
			ff_mult (t, t, _ff_2Sylow_tab[k-1][j]);					// multiply t by S[k-1]=sqrt(S[k]), the product of all these will be b^{-1/2}
		}
		if ( _ff_equal(b,_ff_negone) ) { ff_mult(t, t, w1); break; }
		if ( _ff_one(b) ) break;
	}
#if ! FF_FAST
	_ff_square (b, t);
	if ( ! sts ) ff_mult(b,b,_ff_2g);
	ff_mult(b,b,a[0]);
	if ( ! _ff_one(b) ) {
		printf ("a=%lu, p=%lu\n", _ff_get_ui(a[0]), _ff_p);
		if ( sts ) {
			printf ("ff_2Sylow_invsqrt failed, %lu^2 *  %lu != 1\n", _ff_get_ui(t), _ff_get_ui(a[0]));
		} else {		
			printf ("ff_2Sylow_invsqrt failed, %lu^2 * %lu *  %lu != 1\n", _ff_get_ui(t), _ff_get_ui(_ff_2g), _ff_get_ui(a[0]));
		}
		abort();
	}
#endif
	_ff_set(o[0],t);
	return sts;
}

// Recursive discrete log algorithm - not currently used as it is  slower than ff_2Sylow_invsqrt for p < 2^64
// Given a and b in the 2 Sylow subgroup, returns a nonnegative integer k < 2^e s.t. _a^k = b, or -1 if no such k exists
int ff_2Sylow_dlog (ff_t a, ff_t b, int e)
{
	register ff_t  x0, x1,y;
	register int d, i, k0, k1;
	ff_t t;
//printf ("2Sylow_dlog(%ld,%ld,%d)\n", _ff_get_ui(a), _ff_get_ui(b), e);
	if ( _ff_one(b) ) return 0;
	if ( _ff_equal(b,_ff_negone) ) return (1<<(e-1));			// non-generic optimization specific to fields
	if ( e < 2 ) return -1;
	if ( e == 2 ) {
		if ( _ff_equal (b,a) ) return 1;
		_ff_mult(y,a,b);
		if ( _ff_one(y) ) return 3;
		return -1;
	}
	d = e/2;
	// compute x0=a^(2^(e-d)), y=b^(2^(e-d)), x1 = a^(2^d)
	_ff_set(x0,a); _ff_set(y,b);
	_ff_set_zero(x1);									// avoid compiler warning
	for ( i = 0 ; i < e-d ; i++ ) {
		if ( i == d ) _ff_set(x1,x0);
		ff_square(x0,x0);
		ff_square(y,y);
	}
	if ( i == d ) _ff_set(x1,x0);
	k0 = ff_2Sylow_dlog (x0,y,d);
//printf ("k0=%d\n", k0);
	if ( k0 < 0 ) return -1;
	k1 = (1<<e)-k0;
	ff_exp_ui(&t, &a, k1);
	_ff_mult(y,b,t);
	k1 = ff_2Sylow_dlog(x1,y,e-d);
//printf("k1=%d\n", k1);
	if ( k1 < 0 ) return -1;
	return (1<<d)*k1+k0;
}

void _ff_setup_3g (void)
{
	register int i;
	register unsigned long n;
	ff_t r,s,t;

	_ff_cbrt_setup = 1;
	if ( ! _ff_p1mod3 ) {
		_ff_cbrt_chain_len = ff_precompute_exp_chain (_ff_cbrt_chain, (2*_ff_p-1)/3);
		_ff_invcbrt_chain_len = ff_precompute_exp_chain (_ff_invcbrt_chain, (_ff_p-2)/3);
		return;
	}
	n = (_ff_p-1)/(3*_ff_p3_m);					// p-1 = 3*m*n

	// we could use cubic reciprocity here
	_ff_set_one(t);
	_ff_x2(t);
	for(;;) { 
		ff_exp_ui (&r,&t, _ff_p3_m);				// exponentiation into 3-Sylow group first
		if ( _ff_one(r) ) { _ff_inc(t); continue; }
		if ( _ff_p3_e == 1 ) break;
		ff_exp_ui (&s, &r, n);					// s = t^((p-1)/3)
		if ( ! _ff_one(s) ) break;					// if s is not 1 then r is a generator of the 3-Sylow (and a cubic nonresidue)
		_ff_inc(t);
	}
	_ff_set (_ff_3g, r);
	_ff_set (_ff_3Sylow_tab[0][0], r);
	_ff_square (_ff_3Sylow_tab[0][1],_ff_3Sylow_tab[0][0]);
	for ( i = 1 ; i < _ff_p3_e ; i++ ) {
		_ff_mult(_ff_3Sylow_tab[i][0],_ff_3Sylow_tab[i-1][0],_ff_3Sylow_tab[i-1][1]);		// tab[i][0] = tab[i-1][0]^3
		_ff_square(_ff_3Sylow_tab[i][1], _ff_3Sylow_tab[i][0]);
	}
	_ff_set(_ff_cbrt_unity,_ff_3Sylow_tab[_ff_p3_e-1][0]);
	if ( _ff_p3_m1mod3 ) {
		_ff_cbrt_chain_len = ff_precompute_exp_chain (_ff_cbrt_chain, (_ff_p3_m-1)/3);
	} else {
		_ff_cbrt_chain_len = ff_precompute_exp_chain (_ff_cbrt_chain, (_ff_p3_m-2)/3);
	}
//printf("cbrt_unity = %d\n", _ff_get_ui(_ff_cbrt_unity));
}

/*
   Tonelli-Shanks cuberoot algorithm

    This algorithm is deterministic (given precomputed 3-Sylow generator) - caller should flip coin to randomize choice of cbrt.
    Computes cube root and optionally the inverse of the cube root (oi may be null)
*/
int ff_cbrt_invcbrt (ff_t o[1], ff_t *oi, ff_t a[1])
{
	ff_t b, x, y;

	ff_setup_3g();							// always do this first to make sure _ff_cbrt_unity is set for caller and so we have a chain for exponentiation
	if ( _ff_zero(a[0]) ) { _ff_set_zero (o[0]);  if ( oi ) { printf ("Attempt to invert zero in ff_cbrt_invcbrt!\n");  abort(); } return 1; }
	
	if ( ! _ff_p1mod3 ) {
		// every element of F_p is a cubic residue
		if ( _ff_p == 3 ) {_ff_set(o[0],a[0]);  if ( oi ) _ff_set(*oi,a[0]);  return 1; }
		if ( ! oi ) { ff_exp_chain (o, a, _ff_cbrt_chain, _ff_cbrt_chain_len); return 1; } // a^((2p-1)/3)
		ff_exp_chain (oi, a, _ff_invcbrt_chain, _ff_invcbrt_chain_len);				 // a^((p-2)/3)
		_ff_square(x,*oi);
		_ff_mult(o[0],a[0],x);
		return 1;
	}

	if ( _ff_p3_m1mod3 ) {
		ff_exp_chain (&x, a, _ff_cbrt_chain, _ff_cbrt_chain_len);					// a^((m-1)/3)
		_ff_square (y, x); _ff_mult (b, x, y);		// b = x^3 = a^(m-1)
		ff_mult (b, b, a[0]);					// b = a^m is in the 3-Sylow subgroup
	} else {
		ff_exp_chain (&x, a, _ff_cbrt_chain, _ff_cbrt_chain_len);					// a^((m-2)/3)
		_ff_square (y, x); _ff_mult (b, x, y);		// b = x^3 = a^(m-2)
		_ff_square (y,a[0]); ff_mult (b, b, y);	// b = a^m is in the 3-Sylow subgroup
	}

	// Check if b is the 3-Sylow generator or its square--this happens with probability 2/3^e.  For 2/3 of the primes this is 2/3, which makes it worth avoiding the call to ff_2Sylow_invcbrt
	if ( _ff_equal(b,_ff_3Sylow_tab[0][0]) || _ff_equal(b,_ff_3Sylow_tab[0][1]) ) return 0;
	if ( ! ff_3Sylow_invcbrt(&y,&b) ) return 0;			// y = a^{-m/3}
	ff_mult(x,x,y);									// x = a^{-1/3} or a^{-2/3}
	if ( _ff_p3_m1mod3 ) {
		if ( oi ) _ff_set(*oi,x);
		_ff_square(y,x); _ff_mult(x,y,a[0]);
	} else {
		if ( oi ) { _ff_square(y,x); _ff_mult(*oi,y,a[0]); }
		ff_mult(x,x,a[0]);
	}
//#if ! FF_FAST
	_ff_square (y, x);  ff_mult(y,y,x);
	if ( ! _ff_equal(y,a[0]) ) { printf ("ff_cbrt failed, %lu^3 = %lu != %lu mod %lu\n", _ff_get_ui(x), _ff_get_ui(y), _ff_get_ui(a[0]), _ff_p); abort(); }
//#endif
	_ff_set(o[0],x);
	return 1;
}

// computes a^{-1/3} for a in the 3-Sylow subgroup, returns 0 if not a quadratic residue
// uses precomputed table of 2e powers of the 3-Sylow generator to reduce running time by a factor of 2 over standard Tonelli-Shanks (still O(e^2)).
int ff_3Sylow_invcbrt (ff_t o[1], ff_t a[1])
{
	register ff_t b, q, q1, t, w1, w2;
	register int j, k, s;
	
	ff_setup_3g();										// always do this first to make sure _ff_cbrt_unity is set for caller	
	// handle easy cases first
	if ( _ff_one(a[0]) ) { _ff_set_one(o[0]);  return 1; }		// use 1 as the cube root of 1
	if ( _ff_p3_e == 1 ) return 0;

	// set w1 and w2 to the two elements of order 3 in the 3-Sylow (i.e. the two non-trivial cube roots of unity)
	_ff_set (w1, _ff_3Sylow_tab[_ff_p3_e-1][0]);
	_ff_set (w2, _ff_3Sylow_tab[_ff_p3_e-1][1]);
	_ff_set_one (t);
	_ff_set (b,a[0]);
	j = 0;									// avoid compiler warning
	do {
		_ff_set (q, b);
		for ( s = 1 ; s < _ff_p3_e+1 ; s++ ) {		// s<e+1 is just a safety check in case a isn't in the 3-Sylow, this could be removed
			j=1;
			if ( _ff_equal(q,w1) ) break;
			j=0;
			if ( _ff_equal(q,w2) ) break;
			_ff_set(q1,q);
			ff_square (q,q);  ff_mult(q,q,q1);
		}
		k = _ff_p3_e-s;
#if ! FF_FAST
		if ( k < 0 ) { printf ("Unexpected result: k<0 in ff_3Sylow_invsqrt?!  a = %lu, p = %lu\n", _ff_get_ui(a[0]), _ff_p);  abort(); }
#endif
		if ( k <= 0 ) return 0;
		ff_mult (b, b, _ff_3Sylow_tab[k][j]);		// the product of all the elements S[k] we multiply into b here is b^{-1}, since we terminate with b=1
		ff_mult (t, t, _ff_3Sylow_tab[k-1][j]);		// multiply t by S[k-1]=cbrt(S[k]), the product of all these will be b^{-1/3}
	} while ( !_ff_one(b) );
#if ! FF_FAST
	_ff_square (b, t); ff_mult(b,b,t);
	ff_mult(b,b,a[0]);
	if ( ! _ff_one(b) ) { printf ("ff_3Sylow_invcbrt failed, %lu^3 *  %lu != 1\n", _ff_get_ui(t), _ff_get_ui(a[0])); abort(); }
#endif
	_ff_set(o[0],t);
	return 1;
}

#endif

// standard 4-ary exponentiation (fixed 2-bit window), o and a may be the same ptr
void ff_exp_ui (ff_t o[1], ff_t a[1], unsigned long e)
{
	register int i, j;
	register ff_t c;
	register unsigned long m;
	ff_t b[4];
	
	if ( ! e ) { _ff_set_one (o[0]);  return; }
//	if ( _ff_one(e) ) { _ff_set (o[0], a[0]);  return; }
	// avoid tests to optimize for e <3 or a==0,1,-1
	i = _asm_highbit(e);
	if ( i&1 ) i--;
	m = 3UL<<i;
	_ff_set (b[1], a[0]);
	_ff_square (b[2],b[1]);
	_ff_mult(b[3],b[2],b[1]);
	_ff_set (c, b[(m&e)>>i]);
	for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
		ff_square(c,c);  ff_square(c,c);
		j = (m&e)>>i;
		if ( j ) ff_mult(c,c,b[j]);
	}
//printf ("%ld^%ld=%ld mod %ld\n", _ff_get_ui(a[0]),e,_ff_get_ui(c),_ff_p);
	_ff_set (o[0], c);
}

int ff_precompute_exp_chain (int chain[FF_MAX_CHAIN_LEN], unsigned long e)
{
	register unsigned long m;
	register int i, j, k;
	
	if ( ! e ) return 0;
	m = e;
	i = ui_len(e)-1;
	for ( j = 0 ; ; j+=2 ) {
		if ( i < 2 ) {
			switch(e) {
			case 1: chain[j]=0; chain[j+1] = 0; break;
			case 2: chain[j]=0; chain[j+1] = 1; break;
			case 3: chain[j]=1; if ( j ) chain[j-1] += 1; chain[j+1] = 0; break;
			default: printf("bug in exp_chain e=%ld\n",m); abort();
			}
			break;
		}
		switch( e>>(i-2) ) {
		case 4: chain[j] = 0; chain[j+1]=2; break;
		case 5: chain[j] = 2; if ( j ) chain[j-1] += 2; chain[j+1]=0; break;
		case 6: chain[j] = 1; if ( j ) chain[j-1] += 1; chain[j+1]=1; break;
		case 7: chain[j] = 3; if ( j ) chain[j-1] += 2; chain[j+1]=0; break;
		default: printf("bug in exp_chain e=%ld, i=%d, e>>(i-2) = %ld\n", m, i, e>>(i-2)); abort();
		}
		i-=2;
		if ( ! i ) break;
		e &= (1UL<<i)-1;  i--;
		chain[j+1]+=1;
		if ( ! e ) { chain[j+1] += i; break; }
		while ( ! (e&(1UL<<i)) ) { chain[j+1]++; i--; }
	}
	k = j+2;
	// verify chain
	e = 0;
	for ( i = 0 ; i < k ; i+= 2 ) {
		e += 2*chain[i]+1;
		e <<= chain[i+1];
	}
	if ( e != m ) { printf ("chain verification failed, %ld != %ld\n", e, m); abort(); }
//printf("Verified chain for e=%ld (p=%ld)\n", e, _ff_p);
	return k;
}


void ff_exp_chain (ff_t o[1], ff_t a[1], int chain[], int len)
{
	register int i, j;
	register ff_t c;
	ff_t b[4];
	
	if ( ! len ) { _ff_set_one(o[0]); return; }
	_ff_set (b[0], a[0]);
	_ff_square (b[3],b[0]);
	_ff_mult(b[1],b[0],b[3]);
	_ff_mult(b[2],b[1],b[3]);
	ff_mult(b[3],b[3],b[2]);
	_ff_set(c, b[chain[0]]);
	for ( j = 0 ; j < chain[1] ; j++ ) ff_square(c,c);
	for ( i = 2 ; i < len ; i+=2 ) {
		ff_mult(c,c,b[chain[i]]);
		for ( j = 0 ; j < chain[i+1] ; j++ ) ff_square(c,c);
	}
	_ff_set (o[0], c);
}

// This function is copied from mpzutil.c to make ff.c/ff.h/asm.h self-contained
unsigned long ff_ui_inverse (unsigned long a, unsigned long m)
{
	register unsigned long q, r, r0, r1;
	register long t, t0, t1;

	if ( a >= m ) a = a%m;
	if ( a == 0 ) return 0;

	t0 = 0;  t1 = 1;
	r0 = m;  r1 = a;
	while ( r1 > 0 ) {
		q = r0 / r1;
		r = r0 - q*r1;
		r0 = r1;  r1 = r;
		t = t0 - q*t1;
		t0 = t1;  t1 = t;
	}
	if ( r0 > 1 ) return 0;
	if ( t0 < 0 ) return m - ((unsigned long)(-t0));
	return (unsigned long)t0;
}

/*
  Fast square-root algorithm using a single exponentiation, requires p = 3 mod 4 or p = 5 mod 8.
  Relies on Thm 8.1 on p. 107 of Bressoud "Factorization and Primality Testing".
  Returns 1 if successful, 0 if a is not a quadratic residue.

  Not currently used, ff_sqrt is effectively just as fast.
	
  Currently only supported in single word implementations.  Overlap is OK.
*/

int ff_fast_sqrt (ff_t o[1], ff_t a[1])
{
	register ff_t  i;
	ff_t temp;
	
#if FF_WORDS > 1
	return 0;
#endif
	if ( _ff_zero(a[0]) ) { _ff_set_zero(o[0]);  return 1; }
	if ( (_ff_p&3)==3 ) {
		ff_exp_ui (&temp, a, (_ff_p+1)>>2);
	} else if ((_ff_p&7)==5 ) {
		ff_exp_ui (&temp, a, (_ff_p+3)>>3);
		_ff_square (i,temp);
		if ( _ff_equal (i, a[0]) ) {
			_ff_set(o[0], temp);
			return 1;
		}
		ff_setup_2exp();
		ff_mult(temp, temp, _ff_2exp);
	} else {
		return 0;
	}
	_ff_square(i,temp);
	if ( ! _ff_equal(i,a[0]) ) return 0;
	_ff_set(o[0],temp);
	return 1;
}

/*
	For p=3 mod 4, we may compute y=x^{1/4} as x^{(p+1)^2/16} since (p+1)^2/16 = 2^2/16 = 1/4 mod p-1.
	This yields y^4 = x^{(p+1)^2/4} = x^{(p-1)/2+1}^2 = x^{((p-1)/2)^2 + (p-1) + 1} = +/- x, 
	If y^4 = -x then x has no fourth root, since only one of x or -x has a square root for p=3 mod 4.
*/
int _ff_fast_fourth_root (ff_t *a, ff_t *x)	// _ff_p must be 3 mod 4, not verified!
{
	static mpz_t E;				// use mpz for exponent computation which may overflow 64 bits
	static unsigned long lastp;
	static int init;
	static int chain_len;
	static int chain[FF_MAX_CHAIN_LEN];
	register ff_t s,t;
	
	if ( ! init ) { mpz_init(E); init = 1; }
	if ( lastp != _ff_p ) {
		mpz_set_ui(E,(_ff_p+1)>>2);  mpz_mul(E,E,E);  mpz_mod_ui(E,E,_ff_p-1);
		chain_len = ff_precompute_exp_chain(chain,mpz_get_ui(E));
		lastp = _ff_p;
	}
	_ff_set(s,*x);						// save value of x in case a and x overlap
	ff_exp_chain(a,x,chain,chain_len);
	ff_square(t,*a); ff_square(t,t);
	return _ff_equal(s,t);
}


/*
	For p=11 mod 12, we may compute y=x^{1/6} as x^{(p+1)/4 * (2p-1)/3} since (p+1)/4 * (2p-1)/3 = 1/6 mod p-1
	This yields y^6 = x^{(p+1)/2 * (2p-1)} = x^{((p-1)/2+1) * 1} = +/- x
	If y^6 = -x then x has no sixth root, since only one of x or -x has a square root for p=3 mod 4.
*/
int _ff_fast_sixth_root (ff_t *a, ff_t *x)	// _ff_p must be 11 mod 12, not verified!
{
	static mpz_t E;
	static unsigned long lastp;
	static int init;
	static int chain_len;
	static int chain[FF_MAX_CHAIN_LEN];
	ff_t s,t;
	
	if ( ! init ) { mpz_init(E); init = 1; }
	if ( lastp != _ff_p ) {
		mpz_set_ui(E,(_ff_p+1)>>2);  mpz_mul_ui(E,E,(2*_ff_p-1)/3);  mpz_mod_ui (E,E,_ff_p-1);
		chain_len = ff_precompute_exp_chain(chain,mpz_get_ui(E));
		lastp = _ff_p;
	}
	_ff_set(s,*x);						// save value of x in case a and x overlap
	ff_exp_chain(a,x,chain,chain_len);
	ff_square(t,*a); ff_mult(t,t,*a); ff_square(t,t);
	return _ff_equal(s,t);
}

/*
	For p=3 mod 4, we may compute y=x^{1/8} as x^{(p+1)^3/64} since (p+1)^3/64 = 2^3/64 = 1/8 mod p-1.
	This yields y^8 = x^{(p+1)^3/8} = x^{(p-1)/2+1}^3 = x^{((p-1)/2)^3 + 3*((p-1)/2)^2 + 3*(p-1)/2 + 1} = +/- x, 
	If y^8 = -x then x has no eighth root, since only one of x or -x has a square root for p=3 mod 4.
*/
int _ff_fast_eighth_root (ff_t *a, ff_t *x)
{
	static mpz_t E;
	static unsigned long lastp;
	static int init;
	static int chain_len;
	static int chain[FF_MAX_CHAIN_LEN];
	register ff_t s,t;
	
	if ( ! init ) { mpz_init(E); init = 1; }
	if ( lastp != _ff_p ) {
		if ( (_ff_p&3) != 3 ) { printf ("ff_eighth_root only supported for p=3 mod 4\n"); abort(); }
		mpz_set_ui(E,(_ff_p+1)>>2);  mpz_mul(E,E,E); mpz_mod_ui(E,E,_ff_p-1);  mpz_mul_ui(E,E,(_ff_p+1)>>2);  mpz_mod_ui(E,E,_ff_p-1);
		chain_len = ff_precompute_exp_chain(chain,mpz_get_ui(E));
		lastp = _ff_p;
	}
	_ff_set(s,*x);						// save value of x in case a and x overlap
	ff_exp_chain(a,x,chain,chain_len);
	ff_square(t,*a); ff_square(t,t); ff_square(t,t);
	return _ff_equal(s,t);
}


/*
	For p=3 mod 4, we may compute y=x^{1/16} as x^{(p+1)^4/256} since (p+1)^4/256 = 2^4/256 = 1/16 mod p-1.
	This yields y^16 = x^{(p+1)^4/16} = x^{(p-1)/2+1}^4 = x^{((p-1)/2)^4 + 4*((p-1)/2)^3 + 6*((p-1)/2)^2 +4*(p-1)/2 + 1} = +/- x, 
	If y^16 = -x then x has no sixteenth root, since only one of x or -x has a square root for p=3 mod 4.
*/
int _ff_fast_sixteenth_root (ff_t *a, ff_t *x)
{
	static mpz_t E;
	static unsigned long lastp;
	static int init;
	static int chain_len;
	static int chain[FF_MAX_CHAIN_LEN];
	register ff_t s,t;
	
	if ( ! init ) { mpz_init(E); init = 1; }
	if ( lastp != _ff_p ) {
		if ( (_ff_p&3) != 3 ) { printf ("ff_eighth_root only supported for p=3 mod 4\n"); abort(); }
		mpz_set_ui(E,(_ff_p+1)>>2);  mpz_mul(E,E,E); mpz_mod_ui(E,E,_ff_p-1); mpz_mul(E,E,E);  mpz_mod_ui(E,E,_ff_p-1);
		chain_len = ff_precompute_exp_chain(chain,mpz_get_ui(E));
		lastp = _ff_p;
	}
	_ff_set(s,*x);						// save value of x in case a and x overlap
	ff_exp_chain(a,x,chain,chain_len);
	ff_square(t,*a); ff_square(t,t); ff_square(t,t);
	return _ff_equal(s,t);
}


void ff_setup_fifth(void)
{
	register unsigned long x;
	
	if ( _ff_p > ULONG_MAX>>2 ) {		// avoid overflow
		x = ff_ui_inverse (5,_ff_p);
	} else {
		switch(_ff_p%5) {
		case 1: x = (4*_ff_p+1)/5; break;
		case 2: x = (2*_ff_p+1)/5; break;
		case 3: x = (3*_ff_p+1)/5; break;
		case 4: x = (_ff_p+1)/5; break;
		default: return;
		}
	}
	_ff_set_ui(_ff_fifth,x);
}


void ff_setup_seventh(void)
{
	register unsigned long x;
	
	if ( _ff_p > ULONG_MAX>>3 ) {		// avoid overflow
		x = ff_ui_inverse (7,_ff_p);
	} else {
		switch(_ff_p%7) {
		case 0: x = 0; break;
		case 1: x = (6*_ff_p+1)/7; break;
		case 2: x = (3*_ff_p+1)/7; break;
		case 3: x = (2*_ff_p+1)/7; break;
		case 4: x = (5*_ff_p+1)/7; break;
		case 5: x = (4*_ff_p+1)/7; break;
		case 6: x = (_ff_p+1)/7; break;
		default: return;
		}
	}
	_ff_set_ui(_ff_seventh,x);
}


void ff_setup_eleventh(void)
{
	register unsigned long x;
	
	if ( _ff_p > ULONG_MAX>>4 ) {		// avoid overflow
		x = ff_ui_inverse (11,_ff_p);
	} else {
		switch(_ff_p%11) {
		case 0: x = 0; break;
		case 1: x = (10*_ff_p+1)/11; break;
		case 2: x = (5*_ff_p+1)/11; break;
		case 3: x = (7*_ff_p+1)/11; break;
		case 4: x = (8*_ff_p+1)/11; break;
		case 5: x = (2*_ff_p+1)/11; break;
		case 6: x = (9*_ff_p+1)/11; break;
		case 7: x = (3*_ff_p+1)/11; break;
		case 8: x = (4*_ff_p+1)/11; break;
		case 9: x = (6*_ff_p+1)/11; break;
		case 10: x = (_ff_p+1)/11; break;
		default: return;
		}
	}
	_ff_set_ui(_ff_eleventh,x);
}

void ff_setup_thirteenth(void)
{
	register unsigned long x;
	
	if ( _ff_p > ULONG_MAX>>4 ) {		// avoid overflow
		x = ff_ui_inverse (13,_ff_p);
	} else {
		switch(_ff_p%13) {
		case 0: x = 0; break;
		case 1: x = (12*_ff_p+1)/13; break;
		case 2: x = (6*_ff_p+1)/13; break;
		case 3: x = (4*_ff_p+1)/13; break;
		case 4: x = (3*_ff_p+1)/13; break;
		case 5: x = (5*_ff_p+1)/13; break;
		case 6: x = (2*_ff_p+1)/13; break;
		case 7: x = (11*_ff_p+1)/13; break;
		case 8: x = (8*_ff_p+1)/13; break;
		case 9: x = (10*_ff_p+1)/13; break;
		case 10: x = (9*_ff_p+1)/13; break;
		case 11: x = (7*_ff_p+1)/13; break;
		case 12: x = (_ff_p+1)/13; break;
		default: return;
		}
	}
	_ff_set_ui(_ff_thirteenth,x);
}

void ff_invert_small_int (ff_t z[], int x)
{
	register ff_t t0;
	unsigned long y;
	
	switch (x) {
	case 2: _ff_set(z[0],_ff_half); break;
	case 3: _ff_set(z[0],_ff_third); break;
	case 4: _ff_square(z[0],_ff_half); break;
	case 5: if ( ! _ff_fifth ) ff_setup_fifth(); _ff_set(z[0],_ff_fifth); break;
	case 6: _ff_mult(z[0],_ff_third,_ff_half); break;
	case 7: if ( ! _ff_seventh ) ff_setup_seventh(); _ff_set(z[0],_ff_seventh); break;
	case 8: _ff_square(t0,_ff_half); _ff_mult(z[0],t0,_ff_half); break;
	case 9: _ff_square(z[0],_ff_third); break;
	case 10: if ( ! _ff_fifth ) ff_setup_fifth(); _ff_mult(z[0],_ff_half,_ff_fifth); break;
	case 11: if ( ! _ff_eleventh ) ff_setup_eleventh(); _ff_set(z[0],_ff_eleventh); break;
	case 12: _ff_square(t0,_ff_half); _ff_mult(z[0],t0,_ff_third); break;
	case 13: if ( ! _ff_thirteenth ) ff_setup_thirteenth(); _ff_set(z[0],_ff_thirteenth); break;
	case 14: if ( ! _ff_seventh ) ff_setup_seventh(); ff_mult(z[0],_ff_half,_ff_seventh); break;
	case 15: if ( ! _ff_fifth ) ff_setup_fifth(); ff_mult(z[0],_ff_third,_ff_fifth); break;
	case 16: _ff_square(t0,_ff_half); _ff_square(z[0],t0); break;
	case 18: _ff_square(t0,_ff_third); _ff_mult(z[0],t0,_ff_half); break;
	case 20: if ( ! _ff_fifth ) ff_setup_fifth(); _ff_square(t0,_ff_half); _ff_mult(z[0],t0,_ff_fifth); break;
	case 22: if ( ! _ff_eleventh ) ff_setup_eleventh(); _ff_mult(z[0],_ff_half,_ff_eleventh); break;
	case 24: _ff_square(t0,_ff_half); ff_mult(t0,t0,_ff_half); _ff_mult(z[0],t0,_ff_third); break;
	case 26: if ( ! _ff_thirteenth ) ff_setup_thirteenth(); _ff_mult(z[0],_ff_half,_ff_thirteenth); break;
	case 27: _ff_square(t0,_ff_third); _ff_mult(z[0],t0,_ff_third); break;
	case 28: if ( ! _ff_seventh ) ff_setup_seventh(); _ff_square(t0,_ff_half); _ff_mult(z[0],t0,_ff_seventh); break;
	case 30: if ( ! _ff_fifth ) ff_setup_fifth(); _ff_mult(t0,_ff_half,_ff_third); _ff_mult(z[0],_ff_fifth,t0); break;
	case 32: _ff_square(t0,_ff_half); ff_square(t0,t0); _ff_mult(z[0],t0,_ff_half); break;
	case 33: if ( ! _ff_eleventh ) ff_setup_eleventh(); _ff_mult(z[0],_ff_third,_ff_eleventh); break;
	case 36: _ff_mult(t0,_ff_half,_ff_third); _ff_square(z[0],t0); break;
	default:
		y = ff_ui_inverse (x, _ff_p);
		_ff_set_ui(z[0],y);
	}
	if ( _ff_zero(z[0]) ) { printf ("Attempt to invert 0 (x=%d) over F_%ld in ff_invert_small\n", x, _ff_p); abort(); }
}

void _ff_setup_fractions (void)
{
	ff_t nums[FF_FRACTIONS+1];
	register ff_t n;
	register int i,j;
	
	if ( ! _ff_zero(_ff_frac[1]) ) return;
	_ff_set_one(_ff_frac[1]);  _ff_set(_ff_frac[2],_ff_half); _ff_set(_ff_frac[3],_ff_third); _ff_mult(_ff_frac[4],_ff_half,_ff_half);
	if ( _ff_p <= 5 ) return;
	_ff_set_ui(n,5);
	j = ( _ff_p <= FF_FRACTIONS ? _ff_p-1 : FF_FRACTIONS );
	for ( i = 5 ; i <= j ; i++ ) { _ff_set(nums[i],n); _ff_inc(n); }
	ff_parallel_invert(_ff_frac+5,nums+5,j-4);
}

// this is slower
void ff_dot_product4 (ff_t z[1], ff_t x[], ff_t y[], int n)
{
	register int i;
	register ff_t t0,t1;

//printf ("dot(%d) x[] = ", n); for ( i = 0 ; i < n ; i++ ) printf ("%ld ", _ff_get_ui(x[i])); puts("");
//printf ("dot(%d) y[] = ", n); for ( i = 0 ; i < n ; i++ ) printf ("%ld ", _ff_get_ui(y[i])); puts("");
#if ! FF_FAST
	if ( _ff_p & 0x7000000000000000UL ) { printf ("p=%ld too large in ff_dot_product\n", _ff_p); abort(); }
#endif

	_ff_set_zero(t0);
	for ( i = 0 ; i <= n-4 ; i+=4 ) {
		_ff_sum_4_mults(t1,x[i],x[i+1],x[i+2],x[i+3],y[i+3],y[i+2],y[i+1],y[i]);
		_ff_addto(t0,t1);
	}
	switch ( n-i ) {
	case 1: _ff_mult(t1,x[i],y[i]); _ff_addto(t0,t1); break;
	case 2: _ff_sum_2_mults(t1,x[i],x[i+1],y[i+1],y[i]); _ff_addto(t0,t1); break;
	case 3: _ff_sum_3_mults(t1,x[i],x[i+1],x[i+2],y[i+2],y[i+1],y[i]); _ff_addto(t0,t1); break;
	}
	_ff_set(z[0],t0);
//printf ("dot(x,y,%d) = %ld\n", _ff_get_ui(t0));
}

// for speed, we compare using the internal Montgomery representation.
// applications should make no assumptions about this comparison, 1 > 1+1 is possible.
int _ff_cmp (const void *a, const void *b)
{
	register ff_t *i, *j;
	i = (ff_t *)a; j = (ff_t *)b;
	if ( *i < *j ) return -1;
	if ( *i > *j ) return 1;
	return 0;
}

// sort list so we can search it.  note that for this purpose, we don't care about *how* its sorted
void ff_organize (ff_t a[], int n)
	{ qsort (a, n, sizeof(ff_t), _ff_cmp); }

int ff_find (ff_t x, ff_t a[], int n)
	{ register ff_t *p;  p=bsearch (&x, a, n, sizeof(ff_t), _ff_cmp);  return (p?p-a:-1); }

// returns a pointer to an array continaing binom(n,0),...,binom(n,n) as finite field elements
ff_t *ff_binomials (int n)
{
	if ( n > _ff_binomial_top ) {
		register int i, j, k;
		
		if ( n > FF_BINOMIAL_MAX ) { printf ("n=%d exceedx FF_BINOMIAL_MAX=%d\n", n, FF_BINOMIAL_MAX); abort(); }
		k = ((_ff_binomial_top+1)*(_ff_binomial_top+2))/2;												// k points to new row
		for ( i = _ff_binomial_top+1 ; i <= n ; i++ ) {
			_ff_set_one(_ff_binomial_tab[k]);
			for ( j = k-i+1 ; j < k ; j++ ) _ff_add(_ff_binomial_tab[j+i],_ff_binomial_tab[j-1],_ff_binomial_tab[j]); 	// j points to row above k
			_ff_set_one(_ff_binomial_tab[k+i]);
			k += i+1;
		}
		_ff_binomial_top = n;
	}
	return _ff_binomial_tab+(n*n+n)/2;
}

#if FF_NO_EXT_SETUP != 1

/*
	Solves the norm equation for the current prime modulus, given a negative discriminant D with |D|<4p
	using a modified cornacchia algorithm (Algorithm 1.5.3 in Cohen).  The returned values of t and v are always positive.
*/
int ff_norm_equation (long *t, long *v, long D)
{
	ff_t Dp, xp;
	register long x0, a, b, c, r, L, d, p4, nD;

	assert (FF_BITS <= 62);
	nD = (unsigned long)(-D);
	p4 = _ff_p<<2;
	if ( D>= 0 || nD >= p4 ) return 0;
	if ( (nD&3) && (nD&3)!=3 ) return 0;
	for ( d = D+(long)_ff_p ; d < 0 ; d += (long)_ff_p );
	_ff_set_ui(Dp,d);
	if ( ! ff_sqrt(&xp,&Dp) ) return 0;
	x0 = _ff_get_ui(xp);
	if ( (x0&1) ^ (nD&1) ) x0 = _ff_p-x0;				// make sure parity of x0 matches parity of D
	a = _ff_p<<1;  b = x0;
	L = (long) sqrt(p4);
	while ( b > L ) { r = ui_mod (a,b);  a = b;  b= r; }		// partial Euclidean algorithm
	x0 = p4-b*b;
	c = x0/nD;
	if ( c*nD != x0 ) return 0;
	r = (long)(sqrt(c)+0.5);
	if ( r*r != c ) return 0;
	*v = (long)r;
	*t = (long)b;
	return 1;
}

#endif

/*
	Computes o[i] = prod_{j != i} a[j], using 2n multiplications. 
	Input and output arrays cannot overlap, output array must have space for
	either n+ceil(log_2(n)) entries or the least power of 2 greater than n, whichever is less.
*/
void ff_complementary_products (ff_t o[], ff_t a[], int n)
{
	ff_t *levels[64];
	int lens[64];
	register ff_t *e, *pchild, *pparent, *echild, t0;
	register int i;
	
	if ( n <= 2 ) { if ( n==2 ) { _ff_set(o[0],a[1]); _ff_set(o[1],a[0]); } else if ( n==1 ) _ff_set_one(o[0]); return; }
	
	// build product tree using a[] as the leaf level (implicitly append 1 to each level to make lengths even)
	levels[0] = a;  lens[0] = n;  e = o;
	for ( i = 0 ; lens[i] > 2 ; i++ ) {
		levels[i+1] = e;  pchild = levels[i];  echild = pchild+lens[i];
		for ( pparent = e ; pchild < echild-1 ; pparent++, pchild+=2 ) _ff_mult(pparent[0],pchild[0],pchild[1]);
		if ( pchild < echild ) { _ff_set(pparent[0],pchild[0]); pparent++; }
		lens[i+1] = pparent-e;   e = pparent;
	}
	if ( lens[i] != 2 ) { err_printf ("Bug in ff_products with n=%d, top level %d doesn't have two nodes!\n", n, i); abort(); }
	pchild = levels[i];
	_ff_set(t0,pchild[0]); _ff_set(pchild[0],pchild[1]); _ff_set(pchild[1],t0);									// swap children at top level,  replacing each by its complement
	for ( i-- ; i ; i-- ) {
		pparent = levels[i+1]; pchild = levels[i];  echild = pchild + lens[i];
		for ( ; pchild < echild-1 ; pparent++, pchild+=2 )
			{ _ff_mult(t0,pparent[0],pchild[1]); _ff_mult(pchild[1],pparent[0],pchild[0]); _ff_set(pchild[0],t0);  }	// each childs complement is the product of its parents complement and its sibling
		if ( pchild < echild ) { _ff_set(pchild[0],pparent[0]); }
	}
	// compute bottom level from right to left, to avoid overwriting parent level in the left half
	pparent = levels[1]+lens[1]-1;
	if ( n&1 ) { _ff_set(o[n-1],pparent[0]);  pparent--;  pchild = a+n-3;  echild = o+n-3; } else { pchild = a+n-2;  echild = o+n-2; }
	for ( ; pchild >= a ; pparent--, pchild -= 2, echild -=2 ) { _ff_mult(echild[1],pparent[0],pchild[0]); _ff_mult(echild[0],pparent[0],pchild[1]); }
}


// implementation of Cornacchia's algorithm to find positive (x,y) s.t. x^2+dy^2 = p  (d>0) using Cohen, Alg. 1.5.2
// returns 0 if no solution exists.
int ff_cornacchia (long *x, long *y, long d)
{
	ff_t k, t;
	long a, b, c, r, L, p;
	
	if ( d <= 0 ) { err_printf ("ff_cornacchia: d=%ld must be positive\n", d); abort(); }
	p = _ff_p;
	_ff_set_i (t, -d);
	if ( ! ff_sqrt (&k, &t) ) return 0;
	b = _ff_get_ui(k);
	if ( b < (p>>1) ) b = p-b;
	a = p;
	L = sqrt(p);
	while ( b > L ) { r = a%b;  a=b;  b=r; }
	r = p-b*b;
	c = r/d;
	if ( c*d != r ) return 0;
	if ( ! (*y=i_sqrt(c)) && c ) return 0;
	*x=b;
	return 1;
}


int _qsort_ff_ui_cmp (const void *a, const void *b)
{
	unsigned long x, y;
	
	x = _ff_get_ui (*((ff_t *)a));  y = _ff_get_ui (*((ff_t *)b));
	if ( x < y ) return -1;
	if ( x > y ) return 1;
	return 0;
}
