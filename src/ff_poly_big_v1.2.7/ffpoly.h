#ifndef _FFPOLY_INCLUDE
#define _FFPOLY_INCLUDE

/*
    Copyright 2007-2013 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <assert.h>
#include <gmp.h>
#include "ff.h"
#include "ffext.h"
#include "cstd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FF_POLY_VERSION_STRING		"1.2.7"

typedef struct ff_poly_modulus_struct {
	ff_t *g, *rgi, *w, *w1;					// rgi is the inverse of rev(g) mod x^e
	int d, e;								// d+e is the max degree of a poly that may be reduced
} ff_poly_modulus_t[1];

#define FF_POLY_MAX_DEGREE				1023	// you can exceed this by using ffpolybig.c directly
#define FF_POLY_SMALL_DEGREE			255		// largest small polynomial must be larger than any small crossover values
#define FF_POLY_MULT_SMALL_DEGREE		20		// crossover degree for big mults
#define FF_POLY_SQUARE_SMALL_DEGREE	64		// crossover degree for big squares
#define FF_POLY_MOD_SMALL_DEGREE		255		// crossover degree for mult-based polymods (setup+reduce)

#define FF_POLY_SPLIT					1		// indicates that poly is a product of linear factors over Fp
#define FF_POLY_ONE_ROOT				2		// asks for just one root (even if there are more)
#define FF_POLY_EXACTLY_ONE_ROOT		4		// used in ff_poly_roots_d3 to indicate a root should be returned only when there is exactly one (but the correct count is always returned)

// handy macros - these require integer variables d_f be declared for each poly f
#define _poly_set(a,b)				ff_poly_copy(a,&d_ ## a,b,d_ ## b)
#define _poly_print(a)				ff_poly_print(a,d_ ## a)
#define _poly_neg(c,a)				ff_poly_neg(c,&d_ ## c,a,d_ ## a)
#define _poly_monic(c,a)				ff_poly_monic(c,&d_ ## c,a,d_ ## a)
#define _poly_addto(c,a)				ff_poly_addto(c,&d_ ## c,a,d_ ## a)
#define _poly_add(c,a,b)				ff_poly_add(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_sub(c,a,b)				ff_poly_sub(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_mult(c,a,b)				ff_poly_mult(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_div(q,r,a,b)				ff_poly_div(q,&d_ ## q,r,&d_ ## r,a,d_ ## a,b,d_ ## b)
#define _poly_mod(c,a,b)				ff_poly_mod(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_gcd(c,a,b)				ff_poly_gcd(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_gcdext(c,u,v,a,b)		ff_poly_gcdext(c,&d_ ## c,u,&d_ ## u, v,&d_ ## v,a,d_ ## a,b,d_ ## b)
#define _poly_expmod(c,a,e,b)			ff_poly_exp_mod(c,&d_ ## c,a,&d_ ## a,e,b,&d_ ## b)

static inline ff_t *ff_poly_alloc (int d) { return mem_alloc ((d+1)*sizeof(ff_t)); }
static inline void ff_poly_free (ff_t *f, int d) { mem_free (f); }
static inline void ff_poly_randomize (ff_t f[], int d) { register int i; _ff_random_nz(f[d]); for ( i = 0 ; i < d; i++ ) _ff_random(f[i]); }


static inline void ff_poly_set_mpz (ff_t f[], mpz_t F[], int d)
	{ register int i;  for ( i = 0 ; i <= d ; i++ ) _ff_set_mpz(f[i],F[i]); }

void ff_poly_med_weierstrass (ff_t f[4], ff_t W[5]);
void ff_poly_short_weierstrass (ff_t f[4], ff_t W[5]);

int ff_poly_parse (ff_t f[], int maxd, char *expr);
void ff_poly_print (ff_t f[], int d_f);
int ff_poly_sprint (char *s, ff_t f[], int d_f);

void ff_poly_twist (ff_t g[], ff_t f[], int d);
static inline void ff_poly_twist_mpz (ff_t f[], mpz_t F[], int d)
{
	ff_poly_set_mpz (f, F, d);
	ff_poly_twist (f, f, d);
}

static inline int ff_poly_is_zero (ff_t f[], int d_f) { return d_f == -1 ? 1 : 0; }
static inline int ff_poly_is_one (ff_t f[], int d_f ) { return (!d_f && _ff_one(f[0])) ? 1 : 0; }

int ff_poly_irreducible (ff_t f[], int d_f, int *pnroots);		// pnroots is the number of distinct roots (optional)
int ff_poly_factorization_pattern_and_root (int counts[], ff_t f[], int d_f, ff_t *r);  	// counts[i] = # of irred factors of deg i.  returns total number of irred factors, if r is non-null sets r to a root (if any exist)
static inline int ff_poly_factorization_pattern (int counts[], ff_t f[], int d_f)
	{ return ff_poly_factorization_pattern_and_root (counts, f, d_f, 0); }
static inline int ff_poly_count_factors (ff_t f[], int d_f)		// total number of  irreducible factors (including multiplicity) of f, which need not be monic
	{ return ff_poly_factorization_pattern_and_root (0, f, d_f, 0); }
int ff_poly_factors (ff_t r[], int n[], ff_t f[], int d_f);			// determines the irreducible monic factors of f (which need not be monic).  Returned poly's are implicitly monic (leading 1 omitted), and concatenated in r[], with degrees in n[]
int ff_poly_roots (ff_t *r, ff_t f[], int d_f);					// f must be monic. returns all Fp-roots (not necessarily distinct)
int _ff_poly_distinct_roots (ff_t *r, ff_t f[], int d_f, int flags);	// f must be monic.  FF_POLY_SPLIT flag => f splits into linear factors.  FF_POLY_ONE_ROOT flag => only one root will be computed (if any).  Returns count of distinct roots
static inline int ff_poly_find_root (ff_t *r, ff_t f[], int d_f)	// finds one Fp-root and returns count of distinct roots, or returns 0 if f has no Fp-roots
	{ return _ff_poly_distinct_roots (r,f,d_f,FF_POLY_ONE_ROOT); }
int ff_poly_distinct_roots (ff_t *r, ff_t f[], int d_f);			// f must be monic
int ff_poly_count_roots (ff_t f[], int d_f);					// returns the number of Fp-roots of f, with multiplicity
int ff_poly_count_distinct_roots (ff_t f[], int d_f);			// returns the number of distinct Fp-roots of f
	
int ff_poly_discriminant_nonzero (ff_t f[], int d);			// returns 1 if discriminant is nonzero, 0 otherwise (uses a gcd with derivative, does not actually compute the discriminant)
int ff_poly_gcd (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b);		// output g is not made monic.  returns deg
int ff_poly_gcd_reduce (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b);	// output g is not made monic.  returns deg
int ff_poly_gcdext (ff_t g[], int *pd_g, ff_t u[], int *pd_u, ff_t v[], int *pd_v, ff_t a[], int d_a, ff_t b[], int d_b);	// computes monic extended gcd, returns deg
int ff_poly_inv_mod (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g);
static inline int ff_poly_inv_modulus (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_poly_modulus_t mod)
	{ return ff_poly_inv_mod (h, pd_h, f, d_f, mod->g, mod->d); }
int ff_poly_trivial_gcd (ff_t f[], int d_f, ff_t g[], int d_g);		// very fast code to test whether f and g have a common factor
void _ff_poly_div (ff_t q[], int *pd_q, ff_t r[], int *pd_r, ff_t a[], int d_a, ff_t b[], int d_b, ff_t *work);	// work must hold space for d_a+d_b+2 entries
static inline void ff_poly_div (ff_t q[], int *pd_q, ff_t r[], int *pd_r, ff_t a[], int d_a, ff_t b[], int d_b)	// if division is known to be exact, r and pd_r may be null
	{ ff_t *w = mem_alloc ((d_a+d_b+2)*sizeof(ff_t));  _ff_poly_div (q, pd_q, r, pd_r, a, d_a, b, d_b, w);  mem_free (w); }
void  ff_poly_inverse_mod_xn (ff_t *g, int *d_g, ff_t *f, int d_f, int n, ff_t *work);					// work must hold space for 2n-1 entries
void  ff_poly_log_mod_xn (ff_t *g, int *d_g, ff_t *f, int d_f, int n);
void  ff_poly_exp_mod_xn (ff_t *g, int *d_g, ff_t *f, int d_f, int n);
void ff_poly_antiderivative (ff_t *g, int *d_g, ff_t *f, int d_f);
void ff_poly_fold_mod_xn (ff_t *f, int *pd_f, ff_t *a, int d_a, ff_t *b, int d_b, ff_t *c, int d_c, int n);
void ff_poly_S_mod_xn (ff_t *S, int *pd_S, ff_t A, ff_t B, ff_t *C, int d_C, int n);

void ff_poly_mult_small (ff_t o[], ff_t f[], ff_t g[], int n);				// polys must be the same size (padded to degree n, meaning n+1 coefficients, o must have space for poly of degree 2n
void ff_poly_mult_big (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b);		// requires linking with zn_poly, only called when FF_POLY_BIG is defined
void ff_poly_mult (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b);	// will use ff_poly_mult_big for polys of degree > POLY_SMALL_DEGREE (provided FF_POLY_BIG is defined when ffpoly.c is compiled)

// g must be large enough to hold a poly of degree 2d_f
void ff_poly_square_small (ff_t g[], ff_t f[], int d_f);
void ff_poly_square_big (ff_t g[], ff_t f[], int d_f);					// requires linking with zn_poly
void ff_poly_square (ff_t g[], ff_t f[], int d_f);						// will use ff_poly_square_big for polys of degree > POLY_SMALL_DEGREE (provided FF_POLY_BIG is defined when ffpoly.c is compiled)

// uses ff_poly_from_roots_small for d <= 64 (very fast), otherwise calls ff_poly_from_roots_big if FF_POLY_BIG is defined or ff_poly_from_roots_naive (very slow) if not
void ff_poly_from_roots (ff_t f[], ff_t r[], int d);

// poly mod functions assume g is monic!
void ff_poly_mod_setup_max (ff_poly_modulus_t mod, ff_t *g, int d_g, int max);
static inline void ff_poly_mod_setup (ff_poly_modulus_t mod, ff_t *g, int d_g)
	{ ff_poly_mod_setup_max (mod, g, d_g, 0); }
void ff_poly_mod_clear (ff_poly_modulus_t mod) ;
void ff_poly_mod_reduce (ff_t *h, int *d_h, ff_t *f, int d_f, ff_poly_modulus_t mod);

// important: ff_polymod_small functions assume g is monic, depressed, and has negated coefficients, i.e. g = x^d - g[d-2]*x^(d-2) - g[d-1]x^(d-1) - ... - g[0]
void ff_poly_mod_small (ff_t *h, int *d_h, ff_t *f, int d_f, ff_t *g, int d_g);
void ff_poly_mod_small_inplace (ff_t *f, int d_f, ff_t *g, int d_g);

// these functions do not require the modulus g to be monic or depressed, and the coefficients are not negated
void ff_poly_mod_big (ff_t *h, int *d_h, ff_t *f, int d_f, ff_t *g, int d_g);
void ff_poly_mod (ff_t h[], int *pd_g, ff_t *f, int d_f, ff_t *g, int d_g);

void ff_poly_xn_mod (ff_t g[], int *pd_g, unsigned long e, ff_t f[], int d_f);
void ff_poly_xn_mod_mpz (ff_t g[], int *pd_g, mpz_t e, ff_t f[], int d_f);
void ff_poly_xn_modulus  (ff_t g[], int *pd_g, unsigned long e, ff_poly_modulus_t mod);
void ff_poly_xn_modulus_mpz (ff_t g[], int *pd_g, mpz_t e, ff_poly_modulus_t mod);
void ff_poly_xpan_mod (ff_t g[], int *pd_g, ff_t a, unsigned long e, ff_t f[], int d_f);
void ff_poly_xpan_modulus  (ff_t g[], int *pd_g, ff_t a, unsigned long e, ff_poly_modulus_t mod);
void ff_poly_pow_mod (ff_t g[], int *pd_g, ff_t a[], int d_a, unsigned long e, ff_t f[], int d_f);
void ff_poly_pow_mod_mpz (ff_t g[], int *pd_g, ff_t a[], int d_a, mpz_t e, ff_t f[], int d_f);
void ff_poly_pow_modulus (ff_t g[], int *pd_g, ff_t a[], int d_a, unsigned long e, ff_poly_modulus_t mod);
void ff_poly_pow_modulus_mpz (ff_t g[], int *pd_g, ff_t a[], int d_a, mpz_t e, ff_poly_modulus_t mod);

int ff_poly_sqrt_mod (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g);
int ff_poly_sqrt_modulus (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_poly_modulus_t mod);

void ff_poly_depress_monic_inplace (ff_t s[1], ff_t f[], int d_f);
static inline void ff_poly_depress_monic (ff_t s[1], ff_t g[], ff_t f[], int d_f)
	{ if ( f != g ) { register int i; for ( i = 0 ; i <= d_f ; i++ ) _ff_set(g[i], f[i]); } ff_poly_depress_monic_inplace(s,g,d_f); }

void ff_poly_translate (ff_t g[], int *pd_g, ff_t f[], int d_f, ff_t c);	// g(x)=f(x+c), works in place, currently d_f must be <= FF_BINOMIAL_MAX

int ff_poly_roots_d2 (ff_t r[2], ff_t f[3], int d_f);				// f may be quadratic or linear, needn't be monic
int _ff_poly_roots_d3 (ff_t r[3], ff_t f[4], ff_t *pd, int flags);		// f must be of the form x^3+ax+b, r may be null if only a count is desired, pd is an optional pointer to sqrt(-D/3) (used for 3-torsion)
int _ff_poly_roots_d3_mod3 (ff_t r[3], ff_t f[4]);				// used only when _ff_p==3.  f must be monic
int _ff_poly_roots_d4 (ff_t r[4], ff_t f[5], ff_t *pd, int flags);		// f must be of the form x^3+ax^2+bx+c, r is required, pd is optional
int _ff_poly_roots_d4_mod3 (ff_t r[4], ff_t f[4]);				// used only when _ff_p==3.  f must be monic

// all the ff_poly_g1 functions have been moved to ecurve in smalljac
void ff_poly_xn_mod_d3 (ff_t g[3], unsigned long n, ff_t f[4]);	// computes x^n mod f=x^3+ax+b
void ff_poly_xn_mod_d4 (ff_t g[3], unsigned long n, ff_t f[5]);
void ff_poly_xpan_mod_d3 (ff_t g[3], ff_t a, unsigned long n, ff_t f[4]);
void ff_poly_xpan_mod_d4 (ff_t g[4], ff_t a, unsigned long n, ff_t f[5]);

void ff_poly_xpan_mod_d2 (ff_t g[2], ff_t a, unsigned long n, ff_t f[1]);
void ff_poly_xpan_mod_d2_o (ff_t g[2], ff_t a, unsigned long n, ff_t f[1]);

static inline int ff_poly_degree (ff_t f[], int d_f) { while ( d_f>=0 && _ff_zero(f[d_f]) ) d_f--;  return d_f; }

static inline void ff_poly_copy (ff_t b[], int *pd_b, ff_t a[], int d_a)
{
	register int i;
	
	if ( b!=a ) for ( i = 0 ; i <= d_a ; i++ ) _ff_set(b[i], a[i]);
	if ( pd_b ) *pd_b = d_a;
}

// note that f may have leading zeroes after reversing inplace -- use ff_poly_degree to get the correct degree
static inline void ff_poly_reverse_inplace (ff_t f[], int d)
{
	register ff_t t;
	register int i,j;
	
	if ( d <= 0 ) return;
	j = (d>>1) + (d&1);		// j = ceil(d_f/2)
	for ( i = 0 ; i < j ; i++ ) { _ff_set(t,f[i]);  _ff_set(f[i],f[d-i]);  _ff_set (f[d-i],t); }
}

static inline void ff_poly_reverse (ff_t g[], int *pd_g, ff_t f[], int d_f)
{
	register int i;
	
	if ( g == f ) ff_poly_reverse_inplace (f, d_f); else for ( i = 0 ; i <= d_f ; i++ ) _ff_set (g[i],f[d_f-i]);
	if ( pd_g ) *pd_g = ff_poly_degree (g, d_f);
}

static inline int ff_poly_equal (ff_t a[], int d_a, ff_t b[], int d_b)
{
	register int i;
	
	if ( d_a != d_b ) return 0;
	for ( i = 0 ; i <= d_a ; i++ ) if ( ! _ff_equal(a[i],b[i]) ) return 0;
	return 1;
}

static inline void ff_poly_set_i (ff_t f[], long F[], int d) { register int i;  for ( i = 0 ; i <= d ; i++ ) _ff_set_i(f[i],F[i]); }


// overlap ok
static inline void ff_poly_double (ff_t b[], int *d_b, ff_t a[], int d_a)
{
	register int i;
	
	for ( i = 0 ; i <= d_a ; i++ ) _ff_add(b[i],a[i],a[i]);
	if ( d_b ) *d_b = d_a;
}

// overlap ok
static inline void ff_poly_triple (ff_t b[], int *d_b, ff_t a[], int d_a)
{
	register int i;
	register ff_t t0;
	
	for ( i = 0 ; i <= d_a ; i++ ) { _ff_add(t0,a[i],a[i]); _ff_add(b[i],a[i],t0); }
	if ( d_b ) *d_b = d_a;
}


// overlap ok
static inline void ff_poly_scalar_mult (ff_t c[], int *d_c, ff_t a, ff_t b[], int d_b)
{
	register int i;
	
	for ( i = 0 ; i <= d_b ; i++ ) ff_mult(c[i],a,b[i]);
	if ( d_c ) *d_c = d_b;
}

// overlap not ok, use ff_poly_addto
static inline void ff_poly_add (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b)
{
	register int i,d_c;

	assert (c!=a && c != b);
	if ( d_a < 0 ) { ff_poly_copy (c, pd_c, b, d_b);  return; }
	if ( d_b < 0 ) { ff_poly_copy (c, pd_c, a, d_a);  return; }
	
	d_c = d_a;
	if ( d_b > d_a ) d_c = d_b;
	for ( i = 0 ; i <= d_c ; i++ ) {
		_ff_set_zero (c[i]);
		if ( i <= d_a ) _ff_addto (c[i], a[i]);
		if ( i <= d_b ) _ff_addto (c[i], b[i]);
	}
	*pd_c = ff_poly_degree(c,d_c);
}

// overlap ok
static inline void ff_poly_addto (ff_t c[], int *pd_c, ff_t b[], int d_b)
{
	register int i;

	if ( d_b < 0 ) return;
	for ( i = 0 ; i <= d_b ; i++ ) {
		if ( i > *pd_c ) {
			_ff_set (c[i], b[i]);
			*pd_c = i;
		} else {
			_ff_addto (c[i], b[i]);
		}
	}
	*pd_c = ff_poly_degree(c,*pd_c);
}

// overlap ok
static inline void ff_poly_neg (ff_t c[], int *pd_c, ff_t b[], int d_b)
{
	int i;

	if ( pd_c ) *pd_c = d_b;
	if ( c == b ) { for ( i = 0 ; i <= d_b ; i++ ) ff_negate(c[i]); }
	else { for ( i = 0 ; i <= d_b ; i++ ) _ff_neg(c[i], b[i]); }
}

// overlap not ok - use poly_subfrom instead
static inline void ff_poly_sub (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b)
{
	int d_c;
	int i;

	assert (c!=a && c !=b);
	if ( d_a < 0 ) { ff_poly_neg (c, pd_c, b, d_b);  return; }
	if ( d_b < 0 ) { ff_poly_copy (c, pd_c, a, d_a);  return; }
	
	d_c = d_a;
	if ( d_b > d_a ) d_c = d_b;
	for ( i = 0 ; i <= d_c ; i++ ) {
		_ff_set_zero (c[i]);
		if ( i <= d_a ) _ff_addto (c[i], a[i]);
		if ( i <= d_b ) _ff_subfrom (c[i], b[i]);
	}
	*pd_c = ff_poly_degree(c,d_c);
}

// overlap ok
static inline void ff_poly_subfrom (ff_t c[], int *pd_c, ff_t b[], int d_b)
{
	register int i;

	if ( d_b < 0 ) return;
	for ( i = 0 ; i <= d_b ; i++ ) {
		if ( i > *pd_c ) {
			_ff_neg (c[i], b[i]);
			*pd_c = i;
		} else {
			_ff_subfrom (c[i], b[i]);
		}
	}
	*pd_c = ff_poly_degree(c,*pd_c);
}

// overlap ok
static inline void ff_poly_shift_up (ff_t b[], int *pd_b, ff_t a[], int d_a, int n)
{
	register int i;
	
	for ( i = d_a+n ; i >= n ; i-- ) _ff_set(b[i],a[i-n]);
	while ( i >= 0 ) { _ff_set_zero(b[i]); i--; }
	*pd_b = d_a+n;
}

// overlap ok
static inline void ff_poly_shift_down (ff_t b[], int *pd_b, ff_t a[], int d_a, int n)
{
	register int i;
	
	*pd_b = d_a-n;
	for ( i = 0 ; i <= *pd_b ; i++ ) _ff_set(b[i],a[i+n]);
}

static inline void ff_poly_monic (ff_t c[], int *pd_c, ff_t b[], int d_b)	// c=b ok
{
	register ff_t z;
	int i;
	
	if ( d_b < 0 || _ff_one (b[d_b]) ) { ff_poly_copy (c, pd_c, b, d_b);  return; }
	_ff_invert (z, b[d_b]);
	if ( pd_c ) *pd_c = d_b;
	for ( i = 0 ; i < d_b ; i++ ) ff_mult(c[i],z,b[i]);
	_ff_set_one(c[i]);
}


static inline int ff_poly_x3axb_disc (ff_t d[1], ff_t f[2])					// f(x)=x^3+ax+b, returns D =-4a^3-27b^2
{
	register ff_t t0, t1, t2;
	
	_ff_square(t0,f[0]);  _ff_set_ui(t1,27);  ff_mult(t2,t0,t1);				// t2 = 27f0^2
	_ff_square(t0,f[1]);  _ff_mult(t1,t0,f[1]);  _ff_x2(t1); _ff_x2(t1);		// t1 = 4f1^3
	_ff_addto(t1,t2); _ff_neg(d[0],t1);
	return ! _ff_zero(d[0]);
}

// computes f(x) = x^3+ax+b with a=f1 satisfying y^2=f(x)
static inline void ff_poly_x3axb_from_f1_x_y (ff_t f[3], ff_t f1, ff_t x, ff_t y)
{
	register ff_t t0, t1;
	
	_ff_square(t0,x); _ff_mult(t1,t0,x); _ff_mult(t0,f1,x); _ff_addto(t1,t0);	// x^3 + f1*x
	_ff_square(t0,y); _ff_subfrom(t0,t1);								// b = y^2 - x^3 - a*x
	_ff_set_one(f[3]); _ff_set_zero(f[2]); _ff_set(f[1],f1); _ff_set(f[0],t0);
}

static inline int ff_poly_roots_d1 (ff_t r[1], ff_t f[2])				// f needn't be monic
{
	register ff_t t1;
	if ( _ff_one(f[1]) ) { _ff_set_one(t1); } else { _ff_invert(t1,f[1]); }
	ff_negate(t1); _ff_mult(r[0],t1,f[0]);  return 1;
}

// assumes f monic, char > 3, sets t to translation value: depressed cubic is f(x-t), stored in place
static inline void ff_depress_cubic (ff_t t[1], ff_t f[4])
{
	register ff_t t1, t2, t3, t4;
	
	_ff_mult(t1,_ff_third,f[2]);		// f2/3
	_ff_set(t[0],t1);
	_ff_square(t2,t1);			// f2^2/9
	_ff_add(t3,t2,t2);			// 2f2^2/9
	_ff_sub(t4,t3,f[1]);
	_ff_addto(t2,t3);
	_ff_subfrom(f[1],t2);
	_ff_mult(t2,t1,t4);
	_ff_addto(f[0],t2);
	_ff_set_zero(f[2]);
	// 3M+4A
}

static inline void ff_depress_quartic (ff_t t[1], ff_t f[5])
{
	register ff_t t0, t1, t2, t3, t4;

	_ff_mult(t1,_ff_fourth,f[3]);
	_ff_set(t[0],t1);
	_ff_neg(t0,t1);							// t0 = c
	_ff_square(t2,t0);						// t2 = c^2
	_ff_triple(t3,t2);						// t3 = 3c^2
	_ff_sub (t4,f[2],t3);						// t4 = f2-3c^2
	_ff_sum_2_mults(t1,t0,t2,t4,f[1]);			// t1 = cf1+c^2(f2-3c^2)
	_ff_addto(f[0],t1);						// f0 += cf1+c^2(f2-3c^2)
	_ff_subfrom(t4,t2);						// t4 = f2-4c^2
	_ff_mult(t1,t0,t4);						// t1 = c(f2-4c^2)
	_ff_add(t0,t1,t1);						// t0 = 2c(f2-4c^2)
	_ff_addto(f[1],t0);						// f1 += 2c(f2-4c^2)
	_ff_add(t0,t3,t3);						// t0 = 6c^2
	_ff_subfrom(f[2],t0);						// f2 -= 6c^2
	_ff_set_zero(f[3]);						// f3 = 0
	// 5M+10A
}

static inline int ff_poly_roots_d3 (ff_t r[3], ff_t f[4])			// f must be monic, returns all roots with multiplicity (not nescarilly distinct)
{
	ff_t g[4],t;
	int i, k;

	assert ( _ff_one(f[3]) );
	if ( _ff_zero(f[2]) ) return _ff_poly_roots_d3 (r, f, 0, 0);
	if ( _ff_p == 3 ) return _ff_poly_roots_d3_mod3 (r, f);
	for ( i = 0 ; i <= 3 ; i++ ) _ff_set(g[i],f[i]);
	ff_depress_cubic(&t,g);
	k = _ff_poly_roots_d3 (r, g, 0, 0);
	for ( i = 0 ; i < k ; i++ ) _ff_subfrom(r[i],t); 
	return k;
}

static inline int ff_poly_roots_d4 (ff_t r[4], ff_t f[5])			// f must be monic, returns all roots with multiplicity (not nescarilly distinct)
{
	ff_t g[5], t;
	int i, k;

	assert ( _ff_one(f[4]) );
	if ( _ff_zero(f[3]) ) return _ff_poly_roots_d4 (r, f, 0, 0);
	for ( i = 0 ; i <= 4 ; i++ ) _ff_set(g[i],f[i]);
	ff_depress_quartic(&t,g);
	k = _ff_poly_roots_d4 (r, g, 0, 0);
	for ( i = 0 ; i < k ; i++ ) _ff_subfrom(r[i],t);
	return k;
}


// computes g=f/(x-r), assuming f(r) = 0.  f and g may point to the same place
static inline void ff_poly_remove_root (ff_t g[], ff_t f[], int d, ff_t r[1])
{
	register ff_t r0,t1,t2;
	register int i;

	_ff_set(r0,r[0]);
	_ff_set(t1,f[d]);
	for ( i = d-1 ; i > 0 ; i-- ) {
		_ff_set(t2,f[i]);
		_ff_set(g[i],t1);
		ff_mult(t1,t1,r0);
		_ff_addto (t1,t2);
	}
	_ff_set(g[0],t1);
}

// computes g=f*(x-r).  f is assumed to be monic and of degree > 0.  f and g may point to the same place
static inline void ff_poly_add_root (ff_t g[], ff_t f[], int d, ff_t r[1])
{
	register ff_t r0,t0;
	register int i;

	_ff_neg(r0,r[0]);
	_ff_set(g[d+1],f[d]);
	for ( i = d ; i > 0 ; i-- ) {
		_ff_mult(t0,r0,f[i]);
		_ff_add(g[i],f[i-1],t0);
	}
	_ff_mult(g[0],f[0],r0);
}

// computes g=f/(x-r), f is monic degree 3, f and g may point to the same place
static inline void ff_poly_remove_root_d3 (ff_t g[], ff_t f[], ff_t r[1])
{
	register ff_t t1,t2;
	
	_ff_add(t1,f[2],*r);
	ff_mult(t2,t1,*r);
	_ff_add (g[0],t2,f[1]);
	_ff_set(g[1],t1);
	_ff_set_one(g[2]);
}

// computes g=f/(x-r), f is monic degree 4, f and g may point to the same place
static inline void ff_poly_remove_root_d4 (ff_t g[], ff_t f[], ff_t r[1])
{
	register ff_t r0,t1,t2;

	_ff_set(r0,r[0]);
	_ff_add(t1,f[3],r0);
	_ff_set(t2,f[2]);
	_ff_set(g[2],t1);
	ff_mult(t1,t1,r0);
	_ff_addto (t1,t2);
	_ff_set(t2,f[1]);
	_ff_set(g[1],t1);
	ff_mult(t1,t1,r0);
	_ff_add (g[0],t1,t2);
	_ff_set_one(g[3]);
}

// f must be monic
static inline void ff_poly_remove_root_d2 (ff_t g[], ff_t f[], ff_t r[1])
{
	_ff_add(g[0],f[1],*r);
	_ff_set_one(g[1]);
}

// compute o=f(x)
static inline void ff_poly_eval (ff_t o[1], ff_t f[], int d, ff_t x[1])
{
	register ff_t t,y;
	register int i;
	
	if ( d < 0 ) { _ff_set_zero(o[0]); return; }
	_ff_set (y,f[d]);
	_ff_set (t,*x);
	for ( i = d-1 ; i >= 0 ; i-- ) ff_multadd(y,y,t,f[i]);
	_ff_set(o[0],y);
}

// naive algorithm, uses (d+1)(d+2)/2 mults and d(d+1) adds, overlap not allowed
static inline void ff_poly_square_naive (ff_t o[], ff_t f[], int d)
{
	register ff_t t;
	register int i,j;
	
	for ( i = 0 ; i <= d ; i++ ) { _ff_square(o[2*i],f[i]); if ( i ) _ff_set_zero(o[2*i-1]); }
	for ( i = 0 ; i <= d ; i++ ) for ( j = 0 ; j < i ; j++ ) { _ff_mult(t,f[i],f[j]); _ff_x2(t);  _ff_addto(o[i+j],t); }	
}

void ff_poly_from_roots_small (ff_t f[], ff_t r[], int d);		// d cannot exceed 64

// naive implementation: d(d-1)/2 mults. Both ff_poly_from_roots_small and ff_poly_from_roots_big are much faster
static inline void ff_poly_from_roots_naive (ff_t f[], ff_t r[], int d)
{
	register ff_t t;
	register int i,j;

	_ff_set_one(f[d]);
	_ff_neg(f[d-1],r[d-1]);
	for ( i = d-2 ; i >= 0 ; i-- ) {
		_ff_mult(f[i],r[i],f[i+1]);
		ff_negate(f[i]);
		for ( j = i+1 ; j < d-1 ; j++ ) {
			_ff_mult(t,r[i],f[j+1]);
			_ff_subfrom(f[j],t);
		}
		_ff_subfrom (f[j],r[i]);
	}
}

// g=f is ok
static inline void ff_poly_derivative (ff_t g[], int *pd_g, ff_t f[], int d_f)
{
	register int i;
	register ff_t t;
	
	if ( d_f > 0 ) {
		_ff_set(g[0],f[1]);  _ff_set_one(t);
		for ( i = 1 ; i < d_f ; i++ ) {  _ff_inc(t); _ff_mult(g[i],t,f[i+1]); }
		if ( pd_g ) *pd_g = ff_poly_degree(g,d_f-1);
	} else {
		if ( pd_g ) *pd_g = -1;
	}
}

static inline void ff_poly_mod_xn (ff_t f[], int *pd_f, int n)
{
	if ( *pd_f < n ) return;
	*pd_f = ff_poly_degree(f,n-1);
}

static inline void ff_poly_random (ff_t f[], int d)
	{ register int i;  for ( i = 0 ; i <= d ; i++ ) _ff_random(f[i]); return; }
	
static inline void ff_poly_random_monic (ff_t f[], int d)	// d must be >= 0
	{ register int i;  for ( i = 0 ; i < d ; i++ ) _ff_random(f[i]); _ff_set_one(f[d]); return; }

static inline void ff_poly_random_depressed_monic (ff_t f[], int d)	// d must be >= 1
	{ register int i;  for ( i = 0 ; i < d-1 ; i++ ) _ff_random(f[i]); _ff_set_one(f[d]); _ff_set_zero(f[d-1]); return; }

#ifdef __cplusplus
}
#endif

#endif
