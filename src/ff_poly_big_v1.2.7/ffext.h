#ifndef _FFEXT_INCLUDE_
#define _FFEXT_INCLUDE_

/*
    Copyright 2008-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include "ff.h"

#ifdef __cplusplus
extern "C" {
#endif

// subset of finite field arithmetic supported in degree 2 and 3 extension fields
// elements are represented as polys over F_p[x]/(g(x)) where g(x)=x^2-_ff_2g for degree 2
// and in degree 3, g(x)=x^3-_ff_3g if p=1mod3 else g(x)=x^3-x-_ff_3g

extern int _ff2_cbrt_setup;
extern ff_t _ff2_cbrt_unity[2];
extern ff_t _ff3_zp[3];
extern ff_t _ff3_z2p[3];
extern ff_t _ff3_f[4];

void ff_ext_setup(void);
void _ff3_setup(void);
static inline void ff3_setup(void) { if ( ! _ff3_f[3] ) _ff3_setup(); }
static inline void ff2_setup(void) { if ( ! _ff_2g ) ff_setup_2g(); }


// for convenience
static inline void ff2_random(ff_t o[2]) { _ff_random(o[0]);  _ff_random_nz(o[1]); }
static inline void ff3_random(ff_t o[3]) { do { _ff_random(o[0]);  _ff_random(o[1]);  _ff_random(o[2]); } while ( _ff_zero(o[1]) && _ff_zero(o[2]) ); }
static inline void ff2_set_zero(ff_t o[2]) { _ff_set_zero(o[0]); _ff_set_zero(o[1]); }
static inline void ff3_set_zero(ff_t o[3]) { _ff_set_zero(o[0]); _ff_set_zero(o[1]); _ff_set_zero(o[2]); }
static inline void ff2_set_one(ff_t o[2]) { _ff_set_one(o[0]); _ff_set_zero(o[1]); }
static inline void ff3_set_one(ff_t o[3]) { _ff_set_one(o[0]); _ff_set_zero(o[1]); _ff_set_zero(o[2]); }
static inline void ff2_set_ui (ff_t o[2], unsigned long x) { _ff_set_ui(o[0],x); _ff_set_zero(o[1]); }
static inline void ff3_set_ui (ff_t o[3], unsigned long x) { _ff_set_ui(o[0],x); _ff_set_zero(o[1]); _ff_set_zero (o[2]); }
static inline int ff2_zero(ff_t o[2]) { return ( _ff_zero(o[0]) && _ff_zero(o[1]) ); }
static inline int ff3_zero(ff_t o[3]) { return ( _ff_zero(o[0]) && _ff_zero(o[1]) &&  _ff_zero(o[2]) ); }
static inline int ff2_one(ff_t o[2]) { return ( _ff_one(o[0]) && _ff_zero(o[1]) ); }
static inline int ff3_one(ff_t o[3]) { return ( _ff_one(o[0]) && _ff_zero(o[1]) && _ff_zero(o[2]) ); }
static inline int ff2_neg_one(ff_t o[2]) { return ( _ff_neg_one(o[0]) && _ff_zero(o[1]) ); }
static inline int ff2_pm_one(ff_t o[2]) { return ( _ff_pm_one(o[0]) && _ff_zero(o[1]) ); }
static inline int ff3_neg_one(ff_t o[3]) { return ( _ff_neg_one(o[0]) && _ff_zero(o[1]) &&  _ff_zero(o[2]) ); }
static inline int ff3_pm_one(ff_t o[3]) { return ( _ff_pm_one(o[0]) && _ff_zero(o[1]) && _ff_zero(o[2]) ); }
static inline void ff2_set(ff_t o[2], ff_t a[2]) { _ff_set(o[0],a[0]);  _ff_set(o[1],a[1]);  }
static inline void ff3_set(ff_t o[3], ff_t a[3]) { _ff_set(o[0],a[0]);  _ff_set(o[1],a[1]);  _ff_set(o[2],a[2]);  }
static inline void ff2_set_ff(ff_t o[3], ff_t a) { _ff_set(o[0],a);  _ff_set_zero(o[1]);  }
static inline void ff3_set_ff(ff_t o[3], ff_t a) { _ff_set(o[0],a);  _ff_set_zero(o[1]);  _ff_set_zero(o[2]);  }	
static inline int ff2_equal(ff_t o[3], ff_t a[3]) { return (_ff_equal(a[0],o[0])&&_ff_equal(a[1],o[1])); }
static inline int ff3_equal(ff_t o[3], ff_t a[3]) { return (_ff_equal(a[0],o[0])&&_ff_equal(a[1],o[1])&&_ff_equal(a[2],o[2])); }

static inline void ff2_neg(ff_t o[2], ff_t a[2]) { _ff_neg(o[0],a[0]); _ff_neg(o[1],a[1]); }
static inline void ff3_neg(ff_t o[3], ff_t a[3]) { _ff_neg(o[0],a[0]); _ff_neg(o[1],a[1]); _ff_neg(o[2],a[2]); }
static inline void ff2_negate(ff_t a[2]) { ff_negate(a[0]); ff_negate(a[1]); }
static inline void ff3_negate(ff_t a[3]) { ff_negate(a[0]); ff_negate(a[1]); ff_negate(a[2]); }

// overlap ok (we avoid _ff_add(x,y,x) and _ff_sub(x,y,x)
static inline void ff2_add(ff_t o[2], ff_t a[2], ff_t b[2]) { if ( o == b ) { _ff_addto (o[0], a[0]); _ff_addto(o[1],a[1]); } else {_ff_add(o[0],a[0],b[0]);   _ff_add(o[1],a[1],b[1]); } }
static inline void ff3_add(ff_t o[3], ff_t a[3], ff_t b[3]) { if ( o == b ) {_ff_addto (o[0], a[0]); _ff_addto(o[1],a[1]); _ff_addto(o[2],a[2]); } else { _ff_add(o[0],a[0],b[0]);   _ff_add(o[1],a[1],b[1]);   _ff_add(o[2],a[2],b[2]); } }
static inline void ff2_sub(ff_t o[2], ff_t a[2], ff_t b[2]) { if ( o == b ) { _ff_subfrom(o[0],a[0]); _ff_subfrom(o[1],a[1]); ff2_negate(o); } else { _ff_sub(o[0],a[0],b[0]);   _ff_sub(o[1],a[1],b[1]);  } }
static inline void ff3_sub(ff_t o[3], ff_t a[3], ff_t b[3]) { if ( o == b ) { _ff_subfrom(o[0],a[0]); _ff_subfrom(o[1],a[1]); _ff_subfrom(o[2],a[2]);  ff2_negate(o); } else { _ff_sub(o[0],a[0],b[0]);   _ff_sub(o[1],a[1],b[1]);   _ff_sub(o[2],a[2],b[2]);  } }

#define ff2_norm(o,a) 	ff2_norm_s(o,a,_ff_2g)
static inline void ff2_norm_s (ff_t o[1], ff_t a[2], ff_t s)
{
	register ff_t t1,t2;
	
	_ff_square(t1,a[1]);
	_ff_neg(t2,s);							// it is worth the negation in order to use sum_2_mults
	_ff_sum_2_mults (o[0],a[0],t1,t2,a[0]);		// N(a) = a[0]^2 - s*a[1]^2
}

static inline void ff2_trace (ff_t o[1], ff_t a[2])
	{ ff_add(o[0],a[0],a[0]); }
	
static inline void ff2_bar (ff_t o[2], ff_t a[2])
	{ _ff_set (o[0], a[0]); _ff_neg (o[1], a[1]); }

static inline void ff2_scalar_mult (ff_t o[2], ff_t c, ff_t a[2]) { ff_mult(o[0],c,a[0]); ff_mult(o[1],c,a[1]); }

#define ff2_mult(o,a,b)	ff2_mult_s(o,a,b,_ff_2g)
static inline void ff2_mult_s (ff_t o[2], ff_t a[2], ff_t b[2], ff_t s)	// compute (a[1]z+a[0])(b[1]z+b[0]) mod (z^2-s)
{
	register ff_t t1;

	_ff_mult(t1,a[1],b[1]);
	_ff_sum_2_mults(o[1],a[0],a[1],b[0],b[1]);			// o[1] = a[0]b[1]+a[1]b[0]
	_ff_sum_2_mults(o[0],a[0],s,t1,b[0]);				// o[0] = a[0]^2+a[1]b[1]s
}

// multiplies (b[1]z+b[0])*(z+a) mod (z^2-s)
static inline void ff2_mult_zpa_s (ff_t o[2], ff_t b[2], ff_t a, ff_t s)
{
	register ff_t t1;
	
	_ff_mult(t1,a,b[1]); _ff_addto(t1,b[0]);
	_ff_sum_2_mults(o[0],b[0],b[1],s,a);
	_ff_set(o[1],t1);
}

#define ff2_square(o,a)	ff2_square_s(o,a,_ff_2g)
static inline void ff2_square_s (ff_t o[2], ff_t a[2], ff_t s)
{
	register ff_t t1;

	_ff_square(t1,a[1]);
	_ff_dbl_mult(o[1],a[0],a[1]);
	_ff_sum_2_mults(o[0],a[0],s,t1,a[0]);
}

#define ff2_invert(o,a)	ff2_invert_s(o,a,_ff_2g)
static inline void ff2_invert_s (ff_t o[2], ff_t a[2], ff_t s)
{
	ff_t t;

	ff2_norm_s (&t, a, s);
	_ff_invert (t, t);
	ff2_scalar_mult (o, t, a);
	ff_negate (o[1]);
}

#define ff2_exp_ui(o,a,e)		ff2_exp_ui_s(o,a,e,_ff_2g)
void ff2_exp_ui_s (ff_t o[2], ff_t a[2], unsigned long e, ff_t s);	// computes in F_p[z]/(z^2-s)
void ff2_exp_ui_s_x (ff_t o[2], ff_t a[2], unsigned long e, ff_t s);

// Cube roots in F_p^2 are only relevant when p=2mod3, since otherwise the 3-Sylow of F_p^2 is the 3-Sylow of F_p.
void _ff2_setup_cbrt(void);
static inline void ff2_setup_cbrt(void) { if ( ! _ff2_cbrt_setup ) _ff2_setup_cbrt(); }
int ff2_3Sylow_invcbrt (ff_t o[2], ff_t a[2]);

// Note ff_setup_3g() must be called before using any of the ff3 functions below
// ff3_poly_eval, ff3_sqrt, and ff3_exp_ui do this automatically

static inline void ff3_scalar_mult (ff_t o[3], ff_t c, ff_t a[3]) { ff_mult(o[0],c,a[0]); ff_mult(o[1],c,a[1]); ff_mult(o[2],c,a[2]); }
void ff3_mult (ff_t o[3], ff_t a[3], ff_t b[3]);			// too long to bother inlining
void ff3_square (ff_t o[3], ff_t a[3]);

// exponentiating by p (the Frobenius map) is very fast, potentially only 2M and at most 6M+4A, faster than ff3_mult by a lot
// computes a^p = a[0]+a[1]z^p+a[2]z^{2p} using 2M or 6M+4A, overlap ok
static inline void ff3_exp_p (ff_t o[3], ff_t a[3])
{
//	ff_t t1[3], t2[3];
	register ff_t w0,w1,w2;
	
	if ( _ff_p1mod3 ) {
		// if p=1mod3 we know z^p is a multiple of z and z^2p is a multiple of z^2
		_ff_set(o[0],a[0]);
		ff_mult(o[1],a[1],_ff3_zp[1]);
		ff_mult(o[2],a[2],_ff3_z2p[2]);
	} else {
		_ff_sum_2_mults(w0,a[1],a[2],_ff3_z2p[0],_ff3_zp[0]);
		_ff_sum_2_mults(w1,a[1],a[2],_ff3_z2p[1],_ff3_zp[1]);
		_ff_sum_2_mults(w2,a[1],a[2],_ff3_z2p[2],_ff3_zp[2]);
		_ff_add(o[0],w0,a[0]);
		_ff_set(o[1],w1);
		_ff_set(o[2],w2);
//		ff3_scalar_mult(t1,a[1],_ff3_zp);
//		ff3_scalar_mult(t2,a[2],_ff3_z2p);
//		ff3_add(t1,t1,t2);
//		_ff_addto(t1[0],a[0]);
//		ff3_set(o,t1);
	}
}


extern int _ff3_trace_z2;

// tr(a[0]+a[1]z+a[2]z^2 = 3a[0]+a[1]tr(z)+a[2]tr(z^2) = 3a[0]+a[2]tr(z^2) since tr(z)=0 for z^3-z-s=0 and z^3-s=0
static inline void ff3_trace (ff_t o[1], ff_t a[3])
{
	register ff_t t1,t2;
	
	_ff_add(t1,a[0],a[0]);
	_ff_add(o[0],t1,a[0]);
	if ( ! _ff3_trace_z2 ) return;
	_ff_add(t2,a[2],a[2]);			// we rely on the fact that tr(z^2) is either 0 or 2
	_ff_addto(o[0],t2);	
}

// norm(a)=a*a^p*a^(p^2), uses 16M+17A or 27M+23A, about 2 or 3 ff3_mults.
static inline void ff3_norm (ff_t o[1], ff_t a[3])
{
	ff_t x[3],y[3];
	register ff_t t1;

	ff3_exp_p (x,a);				// x=a^p
	ff3_exp_p (y,x);				// y=a^(p^2)
	ff3_mult(x,x,y);				// x=a^(p^2+p)
	// we know the norm is in F_p, so we only compute the constant coefficient of a*x
	_ff_sum_2_mults(t1,a[1],a[2],x[1],x[2]);
	_ff_sum_2_mults(o[0],a[0],t1,_ff_3g,x[0]);  // o=a0x0+(a2x1+a1x2)s		note that this works for both p=1mod3 or p=2mod3
}

void ff3_minpoly (ff_t f[4], ff_t a[3]);

// these functions evaluate a poly defined over F_p at a point in F_p^2 or F_p^3
void ff2_poly_eval_ff (ff_t y[2], ff_t f[], int d, ff_t x[2]);
void ff3_poly_eval_ff (ff_t y[3], ff_t f[], int d, ff_t x[3]);

// g(x)=f(x+c), works in place, currently d_f must be <= FF_BINOMIAL_MAX
void ff2_poly_translate (ff_t g[], int *pd_g, ff_t f[], int d_f, ff_t c[2]);


// these functions evaluate a poly defined over F_p^2 or F_p^3, computing y = f(x)
// the array f holds all the coefficients, ordered by degree; a total of (2*(d+1) entries in Fp^2, 3*(d+1) entries in Fp^3)
void ff2_poly_eval (ff_t y[2], ff_t f[], int d, ff_t x[2]);
void ff3_poly_eval (ff_t y[3], ff_t f[], int d, ff_t x[3]);

void ff2_nonresidue (ff_t nr[2]);
int ff2_sqrt (ff_t o[2], ff_t a[2]);
int ff3_sqrt (ff_t o[3], ff_t a[3]);

// finds roots of a monic quadratic, returns number of roots found (0 or 2)
// f and r may overlap, note that each element of Fp^2 occupies 2 positions in these arrays (so first root is r[1]*z+r[0] and second is r[3]*z+r[2])
static inline int ff2_poly_roots_d2 (ff_t r[4], ff_t f[4])
{
	ff_t w[2];
	
	if ( ff2_zero (f+2) ) { ff2_neg (w,f);  if ( ! ff2_sqrt (r, w) ) return 0;  ff2_neg (r+2,r); return 2; }
	// w = f[1]^2 - 4*f[0]
	ff2_add (w, f, f); ff2_add (w, w, w); 
	ff2_square (r, f+2);  ff2_sub (w, r, w);
	if ( ! ff2_sqrt (w, w) ) return 0;
	ff2_sub (r, w, f+2);
	ff2_scalar_mult (r, _ff_half, r);
	ff2_sub (r+2, r, w);
	return 2;
}

void ff2_poly_print (ff_t f[], int d_f);

int ff3_trsqrt_zpa_mod_rs (ff_t o[1], ff_t a, ff_t r, ff_t s);			// Computes tr(sqrt(z)) in F_p^3=F_p[z]/(z^3-rz-s).  This is a support function for factoring quartics.

void ff3_exp_ui (ff_t o[3], ff_t a[3], unsigned long e);
void ff3_exp_ui_rs (ff_t o[3], ff_t a[3], unsigned long e, ff_t r, ff_t s);
// inversion not currently supported/needed in degree 3

// computes z^n mod f, where z is our generator for F_p^3
void ff3_zn_mod (ff_t o[3], unsigned long n, ff_t f[4]);			// only looks at f[0] and f[1], assumes f[2]=0 and f[3]=1

// this is a stub that currently just aborts with an unimplemented error
int ff2_poly_discriminant_nonzero (ff_t f[], int d);

// computes d = -4f1^3 - 27f0^2 (returns 0 if d==0, 1 ow)
int ff2_poly_x3axb_disc (ff_t d[2], ff_t f[4]);

// compute d=disc(x^3+f2*x^2+f1*x+f0) = f2^2*(f1^2-f0*f2) - f1^3 in F_3^2 (returns 0 if d==0, 1 ow)
int ff2_poly_disc_char3 (ff_t d[2], ff_t f[6]);

// transforms a Weierstrass equation [a1,a3,a2,a4,a6] into an isomorphic curve y^2=f(x)=x^3+f2*x^2+f1*x+f0
void ff2_poly_med_weierstrass (ff_t f[8], ff_t W[10]);

// transforms a Weierstrass equation [a1,a3,a2,a4,a6] into an isomorphic curve y^2=f(x)=x^3+f1*x+f0 (requires characteristic != 3)
void ff2_poly_short_weierstrass (ff_t f[8], ff_t W[10]);

void ff2_poly_depress_monic_cubic (ff_t g[4], ff_t f[8]);

#ifdef __cplusplus
}
#endif

#endif
