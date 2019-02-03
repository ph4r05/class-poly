#ifndef _MPZUTIL_INCLUDE_
#define _MPZUTIL_INCLUDE_

/*
    Copyright 2007-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// This module is a grab-bag of stuff, a lot of which has nothing to do with
// GMP's multi-precision integer arithmetic (mpz).  This module should be split up and refined

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <gmp.h>
#include "prime.h"
#include "ntutil.h"
#include "cstd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MPZ_E						M_E
#define MPZ_GAMMA				0.577215664901532860606512090082		// euler-mascheroni constant
#define MPZ_EGAMMA				1.78107241799019798523650410311		//e^gamma
#define MPZ_SQRT2					M_SQRT2
#define MPZ_LN2					M_LN2
#define MPZ_PI					M_PI
#define MPZ_LOG2E					1.44269504088896340735992468100189
#define MPZ_HALFLOGPI				0.57236494292470008707171367567653

#ifndef ULONG_BITS
#if ULONG_MAX == 18446744073709551615UL
#define ULONG_BITS 64
#else
#error mpzutil requires an unsigned long to hold (at least) 64 bits
#endif
#endif

#ifndef ULONG_BITS_LOG2
#if ULONG_BITS == 64
#define ULONG_BITS_LOG2			6
#else
#error mpzutil requires an unsigned long to hold (at least) 64 bits
#endif
#endif

#define MPZ_MAX_TINY_PRIME			251								// this is hardwired into ui_factor, don't change without also changing ui_factor!  needs to cover to sqrt of MAX_GCD_INTEGER
#define MPZ_MAX_TINY_INTEGER		255
#define MPZ_TINY_PRIMES			54
#define MPZ_TINY_PRIME_SUM			6081							// sum of primes <= MPZ_MAX_TINY_PRIME

// gcd primes are used in gcd-based trial-division factoring
#define MPZ_GCD_PRIMES				6542
#define MPZ_MAX_GCD_PRIME				65521
#define MPZ_MAX_GCD_INTEGER			0xFFFF

// small primes are precomputed and stored in memory, used to support fast enumeration of larger primes (up to the square of  MAX_SMALL_INTEGER), and for fast factoring of small integers
#define MPZ_SMALL_PRIMES				PRIME_SMALL_PRIMES
#define MPZ_MAX_SMALL_PRIME			PRIME_MAX_SMALL_PRIME
#define MPZ_MAX_SMALL_INTEGER			PRIME_MAX_SMALL_INTEGER

#define MPZ_MAX_SMALL_PRIME_FACTOR_INDEX	309//172						// index of the largest prime <= sqrt(MPZ_MAX_SMALL_INTEGER)

#define MPZ_MAX_ENUM_PRIME			PRIME_MAX_ENUM

#define MPZ_MAX_TREE_FACTORS			100
#define MPZ_MAX_INVERTS				200
#define MPZ_MAX_SMALL_FACTORS			256

#define MPZ_MAX_UI_PP_FACTORS			20
#define MAX_UI_PP_FACTORS				MPZ_MAX_UI_PP_FACTORS			// more sensible name
#define MPZ_MAX_FACTORS				64							// we could be a bit more generous here

#define MPZ_MAX_PELL_SOLUTIONS			50

int mpz_pell_solver (mpz_t x[MPZ_MAX_PELL_SOLUTIONS], mpz_t y[MPZ_MAX_PELL_SOLUTIONS], long d, long m, int h);


static inline int mod (int n, int m) { n %= m;  return ( n<0 ? (m + n) : n ); }
static inline int modl (long n, long m) { n %= m;  return ( n<0 ? (m + n) : n ); }

struct factor_tree_item_struct {
	int w;
	mpz_t factors[MPZ_MAX_TREE_FACTORS];							// list of prime power factors
	struct factor_tree_item_struct *subtrees[MPZ_MAX_TREE_FACTORS];	// subtrees factor p_i - 1 where p_i is the base of the ith pp factor
};

typedef struct factor_tree_item_struct *mpz_ftree_t;

#define _ui_ceil_ratio(a,b)		(((a)%(b)) ? (a)/(b)+1 : (a)/(b))
#define _ui_max(a,b)			((a)>(b)?(a):(b))
#define _ui_min(a,b)			((a)>(b)?(b):(a))

#define i_sgn_c(x)		((x)<0?'-':'+')
#define i_abs(x)		((x)<0?-(x):(x))

#define mpz_get_i(z)			(((long)mpz_get_ui(z))*mpz_sgn(z))
#define mylround(x)			((long)floor((x)+0.5))				// assumes x is positive - workaround problem with builtin lround function

static inline double logfac(long n) { return n*log(n) - (double)n + log((double)n*(1.0+4.0*n*(1.0+2.0*n)))/6.0 + MPZ_HALFLOGPI; }		// Ramanujan approximation, accurate to O(1/n^3)

void mpz_util_init();
void mpz_util_clear();

int mpz_eval_expr (mpz_t o, char *expr);
void mpz_mulm (mpz_t o, mpz_t a, mpz_t b, mpz_t m);
void mpz_powm_tiny (mpz_t o, mpz_t b, unsigned e, mpz_t m);
void mpz_powm_big (mpz_t o, mpz_t b, mpz_t e, mpz_t m);
void mpz_addm (mpz_t o, mpz_t a, mpz_t b, mpz_t m);
void mpz_subm (mpz_t o, mpz_t a, mpz_t b, mpz_t m);
void mpz_subm_ui (mpz_t o, mpz_t a, unsigned long b, mpz_t m);
void mpz_negm (mpz_t a, mpz_t m);
void mpz_set_i (mpz_t o, long n);
void mpq_set_i (mpq_t o, long n);
void mpz_parallel_invert (mpz_t o[], mpz_t a[], unsigned n, mpz_t p);
int mpz_eval_term_ui (mpz_t o, unsigned long numvars, unsigned long vars[], unsigned long exps[]);
char *ui_term_to_string (char *buf, unsigned long numvars, unsigned long vars[], unsigned long exps[]);

int mpz_sqrt_modprime (mpz_t o, mpz_t a, mpz_t p);
int mpz_cornacchia (mpz_t x, mpz_t y, mpz_t d, mpz_t p);
int mpz_cornacchia4 (mpz_t x, mpz_t y, mpz_t d, mpz_t p);

int mpz_factor_small (unsigned long p[], unsigned long h[], mpz_t bigp, mpz_t n, int max_factors, int max_hard_bits);
void mpz_print_factors (mpz_t N);
int mpz_remove_small_squares (mpz_t o, mpz_t n);
int mpz_flatten_small (mpz_t o, mpz_t n);  // Attempts to compute the largest square-free divisor of n but doesn't remove squares of primes larger than MPZ_MAX_GCD_PRIME
int mpz_coarse_part (mpz_t o, mpz_t x, mpz_t y);
void mpz_mul_set (mpz_t o, mpz_t *a, unsigned long n);		// note overwrites array specified by a
int mpz_prime_mult (mpz_t o, mpz_t p, mpz_t maxp, unsigned long maxbits);
void mpz_power_primorial (mpz_t o, unsigned long L);
int mpz_compatible (mpz_t a, mpz_t b);

void mpz_randomm (mpz_t o, mpz_t m);
void mpz_randomb (mpz_t o, int b);

// code moved to prime.h, functions here are for backward compatibility
static inline prime_enum_ctx_t *fast_prime_enum_start_w (long start, long end, long window) { return prime_enum_start_w (start, end, window); }
static inline prime_enum_ctx_t *fast_prime_enum_start (long start, long end, int powers) { return prime_enum_start (start, end, powers); }
static inline long fast_prime_enum(prime_enum_ctx_t *ctx) { return prime_enum(ctx); }
static inline long fast_prime_enum_powers (prime_enum_ctx_t *ctx) { return prime_enum_powers (ctx); }
static inline int fast_prime_enum_exp (prime_enum_ctx_t *ctx) { return prime_enum_exp (ctx); }
static inline int fast_prime_enum_base (prime_enum_ctx_t *ctx) { return prime_enum_base (ctx); }
static inline long fast_prime_power_enum (prime_enum_ctx_t *ctx) { return prime_power_enum(ctx); }
static inline void fast_prime_enum_end (prime_enum_ctx_t *ctx) { prime_enum_end(ctx); }
static inline void fast_prime_start_logging (long start, int window) { prime_start_logging (start, window); }
static inline int fast_prime_check_log (long n) { return prime_check_log (n); }
static inline void fast_prime_stop_logging (void) { prime_stop_logging (); }
static inline int ui_is_small_prime (unsigned long p) { return is_small_prime(p); }
static inline unsigned long ui_small_prime (unsigned long n) { return nth_small_prime(n); }
static inline unsigned long ui_small_prime_index (unsigned long p) { return small_prime_pi(p); }
static inline unsigned long ui_primorial (int w) { return primorial (w); }
static inline unsigned long ui_primorial_phi (int w) { return primorial_phi (w); }
static inline unsigned long ui_pi_bound (unsigned long n) { return pi_bound (n); }

int ui_is_prime (unsigned long p);										// checks log, then uses factoring.  for large p (> 2^40) may resort to a probabilistic test

unsigned long ui_randomm(unsigned long m);
unsigned long ui_randomb (unsigned long b);
unsigned long ui_next_prime (unsigned long p);							// this just calls GMP, use fast_prime_enum if you care at all about speed

unsigned long mpz_randomf (mpz_t o, mpz_t N, mpz_t factors[], unsigned long w);
unsigned long mpz_randomft (mpz_t o, mpz_t N,  mpz_t factors[], mpz_ftree_t factor_trees[], unsigned long w);
mpz_ftree_t mpz_randomfp (mpz_t o, mpz_t N);
mpz_ftree_t mpz_randomfpp (mpz_t o, mpz_t N);
void mpz_clear_ftree (mpz_ftree_t t);

void mpz_randompp (mpz_t Q, mpz_t N);
unsigned long mpz_pp_base (mpz_t b, mpz_t q);
unsigned long ui_pp_base (unsigned long pp);
unsigned long ui_pp_div (unsigned long n, unsigned long p);

unsigned long ui_pdiff (unsigned long a, unsigned long b);		// returns least p s.t. the power of p dividing a and b differ
int ui_factor (unsigned long p[], unsigned long h[], unsigned long n);
int ui_divisor_in_interval (unsigned long p[], unsigned long h[], int w, unsigned long min, unsigned long max);
static inline unsigned long ui_flatten (unsigned long n)			// returns the product of the prime divisors of n
	{ unsigned long p[MAX_UI_PP_FACTORS]; unsigned long h[MAX_UI_PP_FACTORS];  int w = ui_factor (p, h, n);  unsigned long d = 1; for ( int i = 0 ; i < w ; i++ ) d *= p[i];  return d; }

int mpz_remove_small_primes (mpz_t o, mpz_t n, unsigned long exps[], unsigned long maxprimes);		// returns # distinct primes removed
unsigned long mpz_nearprime (mpz_t n, unsigned long L);

unsigned long ui_phi (unsigned long n);						// returns \phi(n) = |Z/nZ|^*
int ui_mu (unsigned long n);									// returns \mu(n) (Moebius function)
unsigned long ui_divisor_function (unsigned long n);				// returns d(n), the number of divisors of n (0 if n=0).
unsigned long ui_binomial (int n, int k);						// only handles n <= 28

// Computes m3 and k3 such that {m3*x+k3} = {m1*x+k1} \cup {m2*x+k2} or returns 0 if the intersection is empty
int i_aseq_intersection (long *m3, long *k3, long m1, long k1, long m2, long k2);

long i_sqrt_modprime (long n, long p);

unsigned long mpz_get_bits_ui (mpz_t a, unsigned long i, unsigned long n);		// gets bits i thru i+n-1 and returns as ui
unsigned long mpz_set_bits_ui (mpz_t a, unsigned long i, unsigned long n, unsigned long bits);		// sets bits i thru i+n-1 and returns as ui
double mpz_log2 (mpz_t a);		// approximation guarunteed to be between floor(lg(a)) and lg(a)

void mpz_reset_counters ();
void mpz_report_counters ();

void mpz_sort (mpz_t a[], unsigned long n);

unsigned long ui_gcd_ext (unsigned long a, unsigned long b, long *x, long *y);
unsigned long ui_crt (unsigned long a, unsigned long M, unsigned long b, unsigned long N);	// returns x cong a%M and cong b%N with 0<=x<M*N, assumes gcd(M,N)=1
unsigned long bach_gcd (long a, long b, unsigned long N);
int ui_compatible  (unsigned long a, unsigned long b);

void ui_crt_coeff (unsigned long a[], unsigned long m[], int n, mpz_t *w);	// computes a[i] = \prod_{j!=i}m[j] mod m[i] for i from 0 to n-1,  w is workspace for 2(n+ceil(lg(n))) mpz_t's
void ui_crt_coeff2 (unsigned long a[], unsigned long m[], int n, mpz_t *w);
void ui_crt_coeff_r (mpz_t M, unsigned long a[], unsigned long m[], int n);

int ui_qsort_cmp (const void *a, const void *b);
int dbl_qsort_cmp (const void *a, const void *b);
unsigned long ui_wt (unsigned long x);								// number of bits in binary rep of x
unsigned long ui_lg_ceil (unsigned long x);								// ceil(lg(x))
#define ui_lg_floor(x)	_asm_highbit(x)								// floor(lg(x))
unsigned long ui_binexp_cost (unsigned long x);
unsigned long ui_get_bits (unsigned long x, unsigned long pos, unsigned long bits);
char *ui_bit_string (char *buf, unsigned long x);

unsigned long mpz_store_bytes (mpz_t a);
void mpz_store (char *s, mpz_t a);
void mpz_retrieve (mpz_t o, char *s);

// compute a*x+b
static inline void mpz_linear_ui (mpz_t o, mpz_t x, unsigned long a, unsigned long b)
{
	mpz_mul_ui(o,x,a);  mpz_add_ui(o,o,b);
}

// compute a*x^2+b*x+c
static inline void mpz_quadratic_ui (mpz_t o, mpz_t x, unsigned long a, unsigned long b, unsigned long c)
{
	mpz_mul_ui(o,x,a);  mpz_add_ui(o,o,b);  mpz_mul (o,o,x);  mpz_add_ui(o,o,c);
}

// compute a3*x^3+a2*x^2+a1*a+a0
static inline void mpz_cubic_ui (mpz_t o, mpz_t x, unsigned long a3, unsigned long a2, unsigned long a1, unsigned long a0)
{
	mpz_mul_ui(o,x,a3);  mpz_add_ui(o,o,a2);  mpz_mul (o,o,x);  mpz_add_ui(o,o,a1);  mpz_mul(o,o,x);  mpz_add_ui(o,o,a0);
}

// compute a4*x^4+a3*x^3+a2*x^2+a1*a+a0
static inline void mpz_quartic_ui (mpz_t o, mpz_t x, unsigned long a4, unsigned long a3, unsigned long a2, unsigned long a1, unsigned long a0)
{
	mpz_mul_ui(o,x,a4);  mpz_add_ui(o,o,a3);  mpz_mul (o,o,x);  mpz_add_ui(o,o,a2);  mpz_mul (o,o,x);  mpz_add_ui(o,o,a1);  mpz_mul(o,o,x);  mpz_add_ui(o,o,a0);
}

// counts multiples of k in [Min,Max] - assumes that the least multiple of k greater than Max fits in an unsigned long (this will be true if k and Max are both less than U
static inline unsigned long ui_multiples_in_range (unsigned long k, unsigned long Min, unsigned long Max)
{
	register unsigned long a, b;
	
	a = _ui_ceil_ratio (Min, k)*k;			// least multiple of k >= Min
	b = (Max/k)*k;						// greatest multiple of k <= Max
	if ( b < a ) return 0;				// we could avoid this using signed arithmetic (and probably should, but we could have overflow problems if Max won't fit in a long)
	return (b-a)/k + 1;
}

static inline unsigned long ui_trim (unsigned long x) { while ( x & 0x1 ) x <<= 1; return x; }
unsigned long ui_eval_expr (char *expr);

extern short *mpz_tiny_sqrts[MPZ_MAX_TINY_PRIME+1];

static inline long tiny_sqrt_modprime (long n, long p) { register long x;  if ( p < 2 || p > MPZ_MAX_TINY_PRIME || ! mpz_tiny_sqrts[p] ) return -2;  x = n%p;  if ( x < 0 ) x += p;  return mpz_tiny_sqrts[p][x]; }

typedef struct NAF_entry_struct {
	char digit;
	char c;
} NAF_entry_t;

extern NAF_entry_t ui_NAF_table[2][4];

void ui_NAF (unsigned long *pbits, unsigned long *nbits, unsigned long n);

static inline int ui_NAF_byte (unsigned char *pbits, unsigned char *nbits, unsigned long n, int c)
{
	register int i,j,k;
	register unsigned char pb, nb;
	
	pb=nb=0;
	for ( i = 0 ; i < 8 ; i++ ) { j=n&3; k=ui_NAF_table[c][j].digit; if ( k>0 ) pb|=(1<<i); if ( k < 0 ) nb|=(1<<i); c = ui_NAF_table[c][j].c; n>>=1; }
	*pbits = pb; *nbits=nb;
	return c;
}

// The ppf datatype for prime-power factorizations  should be regarded as
// a transparent (glass box) datatype - functions are free to access the internal
// members and manipulate them directly, provided they conform to
// the convention that the primes are listed in increasing order of prime
// (not prime power)

#define PPF_MAX_FACTORS		MPZ_MAX_UI_PP_FACTORS

typedef struct ppf_struct {
	unsigned long p[PPF_MAX_FACTORS];
	unsigned long h[PPF_MAX_FACTORS];					// making these unsigned longs is a silly legacy - should be changed
	int w;
} ppf_t[1];

// This function does not check for overflow
static unsigned long inline ppf_eval (ppf_t n)
{
	register unsigned long x;
	register int i, j;
	
	x = 1;
	for ( i = 0 ; i < n->w ; i++ ) for ( j = 0 ; j < n->h[i] ; j++ ) x *= n->p[i];
	return x;
}

static inline void ppf_set_one (ppf_t n) { n->w = 0; }

static inline int ppf_factor (ppf_t f, unsigned long n)
	{ return (f->w = ui_factor(f->p, f->h, n)); }

// divides a by b, placing result in a.  Assumes a is divisible by b
static inline void ppf_divexact (ppf_t a, ppf_t b)
{
	register int i, j;
	
	for ( i = j = 0 ; j < a->w && i < b->w ; j++ )  if ( a->p[j]==b->p[i] ) a->h[j] -= b->h[i++];
//	if ( i < b->w ) { printf ("ppf_divexact not exact!"); exit (0); }
	for ( i = j = 0 ; j < a->w ; j++ ) if ( a->h[j] ) { a->p[i] = a->p[j]; a->h[i] = a->h[j]; i++; }		// remove primes whose exponent is now zero
	a->w = i;
}

static inline void ppf_mult (ppf_t m, ppf_t a, ppf_t b)
{
	register int i, j, k;
	
	for ( i = j = k = 0 ; i < a->w || j < b->w ; k++ ) {
		if ( i < a->w && (j>=b->w || a->p[i] <= b->p[j]) ) {
			m->p[k] = a->p[i];
			m->h[k] = a->h[i];
			if ( j < b->w && b->p[j]==a->p[i] ) {
				m->h[k] += b->h[j];
				j++;
			}
			i++;
		} else {
			m->p[k] = b->p[j];
			m->h[k] = b->h[j];
			j++;
		}
	}
	b->w = k;	
}

static inline void ppf_lcm (ppf_t m, ppf_t a, ppf_t b)
{
	register int i, j, k;
	
	for ( i = j = k = 0 ; i < a->w || j < b->w ; k++ ) {
		if ( i < a->w && (j>=b->w || a->p[i] <= b->p[j]) ) {
			m->p[k] = a->p[i];
			m->h[k] = a->h[i];
			if ( j < b->w && b->p[j]==a->p[i] ) {
				if ( b->h[j] > a->h[i] ) m->h[k] = b->h[j];
				j++;
			}
			i++;
		} else {
			m->p[k] = b->p[j];
			m->h[k] = b->h[j];
			j++;
		}
	}
	b->w = k;
}

static inline void ppf_copy(ppf_t a, ppf_t b) { register int i; for ( i = 0 ; i < b->w ; i++ ) { a->p[i] = b->p[i]; a->h[i] = b->h[i]; } a->w = i; }
	
static inline void ppf_split (ppf_t m1, ppf_t m2, ppf_t m)
{
	register int i, j;
	
	m1->w = (m->w+1)/2;
	for ( i = 0 ; i < m1->w ; i++ ) { m1->p[i] = m->p[i]; m1->h[i] = m->h[i]; }
	m2->w = m->w - m1->w;
	for ( j = 0 ; i < m->w ; i++,j++ ) { m2->p[j] = m->p[i]; m2->h[j] = m->h[i]; }
}

// Given an integer n, computes an array A of length m = n or 4n such that A[p mod m] = kron(n,p) for all odd primes p
static inline int qr_mod (char A[], int maxm, int n)
{
	register int i, j, k, m, p, maxk;
	
	m = ( ! (n&3) || !((n-1)&3) ) ? n : n<<2;
	if ( m < 0 ) m = -m;
	if ( m > maxm ) return 0;
	memset (A, 0, m);
	maxk = (int) ui_phi(m);
	i = 2;
	for ( k = p = 0 ; k < maxk || p < m ; i++ ) {
		p = (int) ui_small_prime(i);
		j = p%m;
		if ( ! (m%p) ) { A[j] = 0; continue; }
		if ( ! A[j] ) { A[j] = legendre (n,p); k++; }
	}
	return m;
}

// given arrays A and B with modular constraints modulo mA and mB computs matrix C and returns modulus mC
// such that A[p mod mA] == 1 and B[p mod mB] == 1 iff C[p mod mC] == 1
static inline int qr_mod_merge (char C[], int maxm, char A[], int mA, char B[], int mB)
{
	register int i, j, k, m;
	
	m = (mA*mB)/ui_gcd(mA,mB);
	if ( m > maxm ) return 0;
	for ( i = j = k = 0 ; i < m ; i++,j++,k++ ) {
		if ( j == mA ) j = 0;  if ( k == mB ) k = 0;
		if ( A[j] == 1 && B[k] == 1 ) C[i] = 1; else C[i] = 0;
	}
	return m;
}

#ifdef __cplusplus
}
#endif

#endif
