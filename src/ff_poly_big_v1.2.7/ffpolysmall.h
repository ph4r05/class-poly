#ifndef _FFPOLY_SMALL_
#define _FFPOLY_SMALL_

/*
    Copyright 2009-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// Heavily optimized inlines for arithmetic operations on low-degree polynomials over Fp.
// This code takes a long time to compile but the speed up at run-time makes it well worth the wait...

#include "ff.h"
#include "ffpolyfromroots.h"
#include "ffpoly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
	The modular poly functions in ffpolysmall assume the modulus is of the form g(x) = x^n - g_{n-1}x^{n-2} - ... - g1x - g0.
	The coefficients of x^n and x^{n+1} are implicity 1 and 0, and the other coefficients are implicitly negated.
	This does not apply to the gcd functions!
*/
void ff_poly_xpan_mod_d2 (ff_t h[2], ff_t a, unsigned long n, ff_t g[1]);
void ff_poly_xn_mod_d3 (ff_t h[3], unsigned long n, ff_t g[2]);
void ff_poly_xpan_mod_d3 (ff_t h[3], ff_t a, unsigned long n, ff_t g[2]);
void ff_poly_xn_mod_d4 (ff_t h[4], unsigned long n, ff_t g[3]);
void ff_poly_xpan_mod_d4 (ff_t h[4], ff_t a, unsigned long n, ff_t g[3]);
void ff_poly_xn_mod_d5 (ff_t h[5], unsigned long n, ff_t g[4]);
void ff_poly_xpan_mod_d5 (ff_t h[5], ff_t a, unsigned long n, ff_t g[4]);
void ff_poly_xn_mod_d6 (ff_t h[6], unsigned long n, ff_t g[5]);
void ff_poly_xpan_mod_d6 (ff_t h[6], ff_t a, unsigned long n, ff_t g[5]);
void ff_poly_xn_mod_d7 (ff_t h[7], unsigned long n, ff_t g[6]);
void ff_poly_xpan_mod_d7 (ff_t h[7], ff_t a, unsigned long n, ff_t g[6]);
void ff_poly_xn_mod_d8 (ff_t h[8], unsigned long n, ff_t g[7]);
void ff_poly_xpan_mod_d8 (ff_t h[8], ff_t a, unsigned long n, ff_t g[7]);
void ff_poly_xn_mod_d10 (ff_t h[10], unsigned long n, ff_t g[9]);
void ff_poly_xpan_mod_d10 (ff_t h[10], ff_t a, unsigned long n, ff_t g[9]);
void ff_poly_xn_mod_d11 (ff_t h[11], unsigned long n, ff_t g[10]);
void ff_poly_xpan_mod_d11 (ff_t h[11], ff_t a, unsigned long n, ff_t g[10]);
void ff_poly_xn_mod_d12 (ff_t h[12], unsigned long n, ff_t g[11]);
void ff_poly_xpan_mod_d12 (ff_t h[12], ff_t a, unsigned long n, ff_t g[11]);
void ff_poly_xn_mod_d13 (ff_t h[13], unsigned long n, ff_t g[12]);
void ff_poly_xpan_mod_d13 (ff_t h[13], ff_t a, unsigned long n, ff_t g[12]);
void ff_poly_xn_mod_d15 (ff_t h[15], unsigned long n, ff_t g[14]);
void ff_poly_xpan_mod_d15 (ff_t h[15], ff_t a, unsigned long n, ff_t g[14]);
void ff_poly_xn_mod_d17 (ff_t h[17], unsigned long n, ff_t g[16]);
void ff_poly_xpan_mod_d17 (ff_t h[17], ff_t a, unsigned long n, ff_t g[16]);
void ff_poly_xn_mod_d19 (ff_t h[19], unsigned long n, ff_t g[18]);
void ff_poly_xpan_mod_d19 (ff_t h[19], ff_t a, unsigned long n, ff_t g[18]);
void ff_poly_xn_mod_d23 (ff_t h[23], unsigned long n, ff_t g[22]);
void ff_poly_xpan_mod_d23 (ff_t h[23], ff_t a, unsigned long n, ff_t g[22]);
void ff_poly_xn_mod_d29 (ff_t h[29], unsigned long n, ff_t g[28]);
void ff_poly_xpan_mod_d29 (ff_t h[29], ff_t a, unsigned long n, ff_t g[28]);
void ff_poly_xn_mod_d31 (ff_t h[31], unsigned long n, ff_t g[30]);
void ff_poly_xpan_mod_d31 (ff_t h[31], ff_t a, unsigned long n, ff_t g[30]);
void ff_poly_xn_mod_d37 (ff_t h[37], unsigned long n, ff_t g[36]);
void ff_poly_xpan_mod_d37 (ff_t h[37], ff_t a, unsigned long n, ff_t g[36]);
void ff_poly_xn_mod_d41 (ff_t h[41], unsigned long n, ff_t g[40]);
void ff_poly_xpan_mod_d41 (ff_t h[41], ff_t a, unsigned long n, ff_t g[40]);
void ff_poly_xn_mod_d43 (ff_t h[43], unsigned long n, ff_t g[42]);
void ff_poly_xpan_mod_d43 (ff_t h[43], ff_t a, unsigned long n, ff_t g[42]);
void ff_poly_xn_mod_small (ff_t h[], unsigned long n, ff_t g[], int d);
void ff_poly_xpan_mod_small (ff_t h[], ff_t a, unsigned long n, ff_t g[], int d);

void ff_poly_exp_mod_2 (ff_t g[2], ff_t h[2], unsigned long n, ff_t f[2]);		// Computes g = h^n mod f, where f is of the form x^2-f1x-f0 and g is degree 1 (possibly with zero leading coefficients),  h and g may overlap
void ff_poly_exp_mod_4 (ff_t h[4], ff_t f[4], unsigned long n, ff_t g[3]);		// computes h = f^n mod g with deg(f)=3, g=x^4 - g2x^2 - g1x - g0 (note signs)
void ff_poly_compose_mod_4 (ff_t h[4], ff_t f1[4], ff_t f2[4], ff_t g[3]);			// computes h(x) = f1(f2(x)) mod g(x)
void ff_poly_compose_mod_6 (ff_t g[6], ff_t h1[6], ff_t h2[6], ff_t f[5]);			// computes g(x) = h1(h2(x)) mod f(x) with overlap permittted, assumes f = x^6 - f4x^4 - f3x^3 - ... - f0


void ff_poly_mod_small (ff_t *h, int *d_h, ff_t *f, int d_f, ff_t *g, int d_g);		// g must be monic, f need not be
void ff_poly_mod_small_inplace (ff_t *f, int d_f, ff_t *g, int d_g);				// g must be monic, f need not be

void ff_poly_mult_small (ff_t o[], ff_t f[], ff_t g[], int n);
void ff_poly_square_small (ff_t o[], ff_t f[], int n);
void ff_poly_square_mod_small  (ff_t o[], ff_t f[], ff_t g[], int n) ;				// deg(f)=n-1, g = x^n - g_{n-2}x^{n-2} - ... - g0 (note signs)

void ff_poly_gcd_small (ff_t h[], int *d_h, ff_t f[], int d_f, ff_t g[], int d_g);		// computes the gcd of two monic polys
void ff_poly_bgcd (ff_t h[], int *d_h, ff_t f[], int d_f, ff_t g[], int d_g);			// returns gcd of two monic polys known to have non-trivial gcd of bounded degree (specified by d_h, assumed 1 if d_h=null)
static inline void ff_poly_clf (ff_t h[], ff_t f[], int d_f, ff_t g[], int d_g)			// returns a common linear factor of two monic polys known to have gcd of degree 1 (the results are undefined if this is not true!)
    { ff_poly_bgcd(h,0,f,d_f,g,d_g); }

int ff_poly_dbl_root_d3 (ff_t r[1], ff_t f[3]);

// computes cubic resolvent of f=x^4+f2x^2+f1x+f0 and depresses it
static inline void ff_depressed_cubic_resolvent (ff_t t[1], ff_t g[4], ff_t f[5])
{
	register ff_t t1,t2;

	_ff_set_one(g[3]);
	_ff_add(t1,f[2],f[2]);
	_ff_neg(g[2],t1);															// g2 = -2f2
	_ff_add(t1,f[0],f[0]); _ff_x2(t1);
	_ff_square(t2,f[2]);
	_ff_sub(g[1],t2,t1);															// g1 = f2^2-4f0
	_ff_square(g[0],f[1]);														// g0 = f1^2
	ff_depress_cubic (t,g);
	// total 5M+9A
}

static inline void poly_print(ff_t f[], int d)
{
	register int i;
	printf("["); for ( i = d ; i >= 0 ; i-- ) printf ("+%lu*X^%d", _ff_get_ui(f[i]), i); printf ("]");
}

// computes h(x) and nonzero scalar c such that f*h=c mod g (assumes f is invertible mod g = x^4 -g2x^2 - g1x - g0, note signs!)
static inline void ff_poly_invert_mod_4 (ff_t h[4], ff_t c[1], ff_t f[4], ff_t g[3])
{
	register ff_t t0,t1,t2,s0,s1,f0,f1,f2,f3,f32,w0,w1,w2,w3,w4,w5,w6;
	
	_ff_set(f0,f[0]); _ff_set(f1,f[1]); _ff_set(f2,f[2]); _ff_set(f3,f[3]);
	_ff_neg(w0,f[3]); _ff_square(t0,w0); _ff_neg(f32,t0);
	// compute t2x^2+t1x+t0 = f3^2*g - (f3x-f2)*f
	_ff_sum_3_mults(t2,g[2],f1,f2,f2,w0,f32);
	_ff_sum_3_mults(t1,g[1],f0,f1,f2,w0,f32);
	_ff_sum_2_mults(t0,g[0],f0,f2,f32);

	// compute s1x+s0 = t2^2*f - (t2f3x+t2f2-t1f3)*t
	_ff_square(w2,t2);  _ff_square(w1,t1);  _ff_neg(w4,t0);		// w1 = t1^2, w2 = t2^2
	_ff_sum_2_mults(w3,f2,f3,t0,t1);  _ff_neg(w5,t2);
	_ff_sum_3_mults(s1,f1,w3,f3,w1,w5,w2);					// s1 = t2^2f1-t2(f3t0+f2t1)+t1^2f3
	_ff_sum_2_mults(w3,t1,t2,f2,w0);						// w3 = f2t2-f3t1
	_ff_sum_2_mults(s0,f0,w4,w3,w2);						// s0 = t2^2f0 - t0(f2t2-f3t1)

	// compute c = s1^2*t - (s1t2x+s1t1-s0t2)*s
	_ff_sum_2_mults(w4,s0,s1,t1,w5);						// w4 = s1t1-s0t2
	_ff_square(w1,s1); _ff_neg(w6,s0);						// w1 = s1^2
	_ff_sum_2_mults(c[0],w4,t0,w1,w6);					// c = s1^2t0 - s0(s1t1-s0t2);

	// compute h = -s1^2*(f3x-f2) + (s1t2x+s1t1-s0t2)(t2^2+(t2f3x+t2f2-t1f3)(f3x-f2) satisfying h*f = c mod g
	_ff_mult(w5,w2,f32);
	_ff_mult(h[3],s1,w5);								// h3 = -s1t2^2f3^2
	_ff_neg(t0,s0); _ff_mult (h[2],t0,w5);					// h2 = s0t2^2f3^2
	_ff_mult(w5,f2,w3); _ff_sub(w6,w5,w2); _ff_mult(w3,s1,t2);	// w6 = f2(f2t2-f3t1)-t2^2, w3 = s1t2
	_ff_neg(t0,t1); _ff_mult(w2,t0,f32);						// w2 = t1f3^2
	_ff_sum_3_mults(h[1],w1,w3,w2,w4,w6,w0);				// h1 = -s1^2f3 + s1t2((t2f2-t1f3)f2-t2^2) + (s1t2-s0t2)t1f3^2
	_ff_sum_2_mults(h[0],w4,w1,f2,w6);					// h0 = (s1t2-s0t2)((t2f2-t1f3)f2-t2^2) + s1^2f2
	// 36M + 20A (21 redc)
}

/*
	The gcd_linear functions all require that f and g have a common linear factor (their behavior is not defined otherwise).
	They don't check for leading zero coefficients that can arise, and consequently may set h[1]=0. 
	The caller should check this and revert to the more general (but slower) ff_poly_gcd_linear if this happens.
*/

static inline void ff_poly_mgcd_linear_2_2 (ff_t h[2], ff_t f[2], ff_t g[2])				// f monic deg 2, g monic deg 2, f and g are not modified, h may overlap either
{
	_ff_sub(h[1],f[1],g[1]); _ff_sub(h[0],f[0],g[0]);
	// 2A
}

static inline void ff_poly_gcd_linear_2_2 (ff_t h[2], ff_t f[2], ff_t g[2])				// f deg 2, g deg 2, f and g are not modified, h may overlap either
{
	register ff_t t2;
	
	_ff_neg(t2,f[2]);
	_ff_sum_2_mults(h[1],f[1],t2,g[1],g[2]);
	_ff_sum_2_mults(h[0],f[0],t2,g[0],g[2]);
	// 4M+3A
}

static inline void ff_poly_mgcd_linear_3_2 (ff_t h[2], ff_t f[3], ff_t g[2])				// f monic deg 3, g monic deg 2, f and g are not modified, h may overlap either
{
	register ff_t t0,t1,t2;

	_ff_sub(t2,f[2],g[1]);
	_ff_sub(t1,f[1],g[0]); _ff_mult(t0,t2,g[1]); _ff_sub(h[1],t1,t0);
	_ff_mult(t0,t2,g[0]); _ff_sub(h[0],f[0],t0);
	// 2M+4A
}

static inline void ff_poly_mgcd_linear_4_2 (ff_t h[2], ff_t f[4], ff_t g[2])				// f monic deg 4, g monic deg 2, f and g are not modified, h may overlap either
{
	register ff_t s0,s1,t0,t1,t;

	_ff_set(s0,g[0]); _ff_set(s1,g[1]);
	_ff_sub(t1,f[3],s1);  _ff_sub(t0,f[2],s0); _ff_mult(t,t1,s1); _ff_subfrom(t0,t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_sub(h[1],f[1],t); _ff_mult(t,t0,s0); _ff_sub(h[0],f[0],t);
	// 4M+6A
}

static inline void ff_poly_mgcd_linear_6_2 (ff_t h[2], ff_t f[6], ff_t g[2])				// f monic deg 6, g monic deg 2, f and g are not modified, h may overlap either
{
	register ff_t s0,s1,t0,t1,t;

	_ff_set(s0,g[0]); _ff_set(s1,g[1]);
	_ff_sub(t1,f[5],s1);  _ff_sub(t0,f[4],s0); _ff_mult(t,t1,s1); _ff_subfrom(t0,t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[3],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[2],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1);
	_ff_sub(h[1],f[1],t); _ff_mult(t,t0,s0); _ff_sub(h[0],f[0],t);
	// 8M+10A
}

static inline void ff_poly_mgcd_linear_8_2 (ff_t h[2], ff_t f[8], ff_t g[2])				// f monic deg 8, g monic deg 2, f and g are not modified, h may overlap either
{
	register ff_t s0,s1,t0,t1,t;

	_ff_set(s0,g[0]); _ff_set(s1,g[1]);
	_ff_sub(t1,f[7],s1);  _ff_sub(t0,f[6],s0); _ff_mult(t,t1,s1); _ff_subfrom(t0,t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[5],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[4],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[3],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[2],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1);
	_ff_sub(h[1],f[1],t); _ff_mult(t,t0,s0); _ff_sub(h[0],f[0],t);
	// 12M+14A
}

static inline void ff_poly_mgcd_linear_12_2 (ff_t h[2], ff_t f[12], ff_t g[2])				// f monic deg 12, g monic deg 2, f and g are not modified, h may overlap either
{
	register ff_t s0,s1,t0,t1,t;

	_ff_set(s0,g[0]); _ff_set(s1,g[1]);
	_ff_sub(t1,f[11],s1);  _ff_sub(t0,f[10],s0); _ff_mult(t,t1,s1); _ff_subfrom(t0,t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[9],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[8],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[7],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[6],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[5],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[4],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[3],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[2],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1);
	_ff_sub(h[1],f[1],t); _ff_mult(t,t0,s0); _ff_sub(h[0],f[0],t);
	// 20M+18A
}

static inline void ff_poly_mgcd_linear_14_2 (ff_t h[2], ff_t f[12], ff_t g[2])				// f monic deg 14, g monic deg 2, f and g are not modified, h may overlap either
{
	register ff_t s0,s1,t0,t1,t;

	_ff_set(s0,g[0]); _ff_set(s1,g[1]);
	_ff_sub(t1,f[13],s1);  _ff_sub(t0,f[12],s0); _ff_mult(t,t1,s1); _ff_subfrom(t0,t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[11],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[10],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[9],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[8],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[7],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[6],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[5],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[4],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[3],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1); _ff_set(t1,t0); _ff_sub(t0,f[2],t);
	_ff_sum_2_mults(t,t0,t1,s0,s1);
	_ff_sub(h[1],f[1],t); _ff_mult(t,t0,s0); _ff_sub(h[0],f[0],t);
	// 24M+20A
}

static inline void ff_poly_mgcd_linear_n_2 (ff_t h[2], ff_t *f, int d_f, ff_t g[2])		// f monic deg d_f >= 2, g monic deg 2, f and g are not modified, h may overlap either
{
	register int i;
	register ff_t s0, s1, t0, t1, t;
	
	if ( d_f == 2 ) { ff_poly_mgcd_linear_2_2(h,f,g); return; }
	_ff_set(s0,g[0]); _ff_set(s1,g[1]);
	_ff_sub(t1,f[d_f-1],s1);
	_ff_sub(t0,f[d_f-2],s0); _ff_mult(t,t1,s1); _ff_subfrom(t0,t);
	for ( i = d_f-3 ; i > 0 ; i-- ) {
		_ff_sum_2_mults(t,t0,t1,s0,s1);
		_ff_set(t1,t0); _ff_sub(t0,f[i],t);
	}
	_ff_set(h[1],t0);
	_ff_mult(t,t1,s0); _ff_sub(h[0],f[0],t);
	// (2d-4)M  + (2d-2)A (using d-1 redc's)
}

static inline void ff_poly_gcd_linear_4_2 (ff_t h[2], ff_t f[4], ff_t g[3])				// f deg 4 monic, g deg 2 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, s02, t0, t2, t3;
	
	_ff_neg(s0,g[0]); _ff_neg(s1,g[1]); _ff_set(s2,g[2]); ff_mult(s02,s0,s2);
	_ff_mult(t3,f[3],s2); _ff_addto(t3,s1);
	_ff_square(t0,s2); 
	_ff_sum_2_mults(t2,f[2],t3,s1,t0); _ff_addto(t2,s02);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(h[1],f[1],t2,t3,s02,s1,t0);
	_ff_sum_2_mults(h[0],f[0],t2,s0,t0);
	// 11M  + 7A (using 7 redc's)
}

static inline void ff_poly_gcd_linear_6_2 (ff_t h[2], ff_t f[6], ff_t g[3])				// f deg 6 monic, g deg 2 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, s02, t0, t2, t3;
	
	_ff_neg(s0,g[0]); _ff_neg(s1,g[1]); _ff_set(s2,g[2]); ff_mult(s02,s0,s2);
	_ff_mult(t3,f[5],s2); _ff_addto(t3,s1);
	_ff_square(t0,s2); 
	_ff_sum_2_mults(t2,f[4],t3,s1,t0); _ff_addto(t2,s02);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[3],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[2],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(h[1],f[1],t2,t3,s02,s1,t0);
	_ff_sum_2_mults(h[0],f[0],t2,s0,t0);
	// 19M  + 11A (using 11 redc's)
}

static inline void ff_poly_gcd_linear_8_2 (ff_t h[2], ff_t f[8], ff_t g[3])				// f deg 8 monic, g deg 2 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, s02, t0, t2, t3;
	
	_ff_neg(s0,g[0]); _ff_neg(s1,g[1]); _ff_set(s2,g[2]); ff_mult(s02,s0,s2);
	_ff_mult(t3,f[7],s2); _ff_addto(t3,s1);
	_ff_square(t0,s2); 
	_ff_sum_2_mults(t2,f[6],t3,s1,t0); _ff_addto(t2,s02);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[5],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[4],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[3],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[2],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(h[1],f[1],t2,t3,s02,s1,t0);
	_ff_sum_2_mults(h[0],f[0],t2,s0,t0);
	// 27M  + 15A (using 15 redc's)
}

static inline void ff_poly_gcd_linear_12_2 (ff_t h[2], ff_t f[12], ff_t g[3])				// f deg 12 monic, g deg 2 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, s02, t0, t2, t3;
	
	_ff_neg(s0,g[0]); _ff_neg(s1,g[1]); _ff_set(s2,g[2]); ff_mult(s02,s0,s2);
	_ff_mult(t3,f[11],s2); _ff_addto(t3,s1);
	_ff_square(t0,s2); 
	_ff_sum_2_mults(t2,f[10],t3,s1,t0); _ff_addto(t2,s02);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[9],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[8],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[7],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[6],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[5],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[4],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[3],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[2],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(h[1],f[1],t2,t3,s02,s1,t0);
	_ff_sum_2_mults(h[0],f[0],t2,s0,t0);
	// 59M  + 23A (using 23 redc's)
}

static inline void ff_poly_gcd_linear_14_2 (ff_t h[2], ff_t f[14], ff_t g[3])				// f deg 12 monic, g deg 2 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, s02, t0, t2, t3;
	
	_ff_neg(s0,g[0]); _ff_neg(s1,g[1]); _ff_set(s2,g[2]); ff_mult(s02,s0,s2);
	_ff_mult(t3,f[13],s2); _ff_addto(t3,s1);
	_ff_square(t0,s2); 
	_ff_sum_2_mults(t2,f[12],t3,s1,t0); _ff_addto(t2,s02);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[11],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[10],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[9],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[8],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[7],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[6],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[5],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[4],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t3,f[3],t2,t3,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(t2,f[2],t3,t2,s02,s1,t0);
	ff_mult(t0,t0,s2);
	_ff_sum_3_mults(h[1],f[1],t2,t3,s02,s1,t0);
	_ff_sum_2_mults(h[0],f[0],t2,s0,t0);
	// 59M  + 23A (using 23 redc's)
}

static inline void ff_poly_gcd_linear_n_2 (ff_t h[2], ff_t *f, int d_f, ff_t g[3])			// f deg d_f > 2, g deg 2,  neither need be monic, f and g are not modified, h may overlap either
{
	register int i;
	register ff_t s0, s1, s2, t0, t1, t;
	
	_ff_neg(s0,g[0]); _ff_neg(s1,g[1]); _ff_square(s2,g[2]); ff_mult(s0,s0,g[2]);
	_ff_sum_2_mults(t1,s1,g[2],f[d_f-1],f[d_f]);
	_ff_sum_2_mults(t0,s1,s2,f[d_f-2],t1); _ff_addto(t0,s0);
	for ( i = d_f-3 ; i > 0 ; i-- ) {
		ff_mult(s2,s2,g[2]);
		_ff_sum_3_mults(t,s0,s1,s2,f[i],t0,t1);
		_ff_set(t1,t0); _ff_set(t0,t);
	}
	_ff_set(h[1],t0);
	_ff_neg(s0,g[0]);
	_ff_sum_2_mults(h[0],s0,s2,f[0],t1);
	// (4d-4)M  + (2d)A (using 2d redc's)
}

static inline void ff_poly_mgcd_linear_3_3 (ff_t h[2], ff_t f[5], ff_t g[3])				// f monic deg 3, g monic deg 3, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, t0, t1, t2, t;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]);
	_ff_sub(t2,s2,f[2]); _ff_sub(t1,f[1],s1); _ff_sub(t0,f[0],s0);					// negate t2 when setting it here to save doing it below
	
	// t2,t1,t0 is now f mod g.  reduce g =  1,s2,s1,s0 mod -t2,t1,t0
	_ff_mult(t,t2,s2);  _ff_add(s2,t,t1);
	_ff_mult(t,t2,s1);  _ff_add(s1,t,t0);  _ff_sum_2_mults(h[1],s1,s2,t1,t2);
	_ff_mult(t,t2,s0); _ff_sum_2_mults(h[0],t,s2,t0,t2);
	// total 7M+7A (5 redc)
}


static inline void ff_poly_mgcd_linear_4_3 (ff_t h[2], ff_t f[4], ff_t g[3])				// f monic deg 4, g monic deg 3, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, t0, t1, t2, t3, t;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]);
	_ff_sub(t3,f[3],s2);
	_ff_sub(t2,s1,f[2]); _ff_mult(t,t3,s2); _ff_addto(t2,t);						// negate t2 when setting it here to save doing it below
	_ff_sub(t1,f[1],s0); _ff_mult(t,t3,s1); _ff_subfrom(t1,t);
	_ff_mult(t,t3,s0); _ff_sub(t0,f[0],t);
	
	// t2,t1,t0 is now f mod g.  reduce g =  1,s2,s1,s0 mod -t2,t1,t0
	_ff_mult(t,t2,s2);  _ff_add(s2,t,t1);
	_ff_mult(t,t2,s1);  _ff_add(s1,t,t0);  _ff_sum_2_mults(h[1],s1,s2,t1,t2);
	_ff_mult(t,t2,s0); _ff_sum_2_mults(h[0],t,s2,t0,t2);
	// total 10M+10A  (8 redc)
}


static inline void ff_poly_mgcd_linear_6_3 (ff_t h[2], ff_t f[6], ff_t g[3])				// f monic deg 6, g monic deg 3, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, t0, t1, t2, t, u;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]);
	_ff_sub(t2,f[5],s2);
	_ff_sub(t1,f[4],s1); _ff_mult(t,t2,s2); _ff_subfrom(t1,t);
	_ff_sub(t0,f[3],s0); _ff_sum_2_mults(t,t1,t2,s1,s2); _ff_subfrom(t0,t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[2],t);
	_ff_sum_2_mults(t,t1,t2,s0,s1);
	_ff_neg(t2,t0); _ff_mult(u,t1,s0);										// negate t2 when setting it here to save doing it below
	_ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);
	
	// t2,t1,t0 is now f mod g.  reduce g =  1,s2,s1,s0 mod -t2,t1,t0
	_ff_mult(t,t2,s2);  _ff_add(s2,t,t1);
	_ff_mult(t,t2,s1);  _ff_add(s1,t,t0);  _ff_sum_2_mults(h[1],s1,s2,t1,t2);
	_ff_mult(t,t2,s0); _ff_sum_2_mults(h[0],t,s2,t0,t2);
	// total (3d-2)M  (using d+4 redc's)
}


static inline void ff_poly_mgcd_linear_8_3 (ff_t h[2], ff_t f[8], ff_t g[3])				// f monic deg 8, g monic deg 3, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, t0, t1, t2, t, u;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]);
	_ff_sub(t2,f[7],s2);
	_ff_sub(t1,f[6],s1); _ff_mult(t,t2,s2); _ff_subfrom(t1,t);
	_ff_sub(t0,f[5],s0); _ff_sum_2_mults(t,t1,t2,s1,s2); _ff_subfrom(t0,t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[4],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[3],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[2],t);
	_ff_sum_2_mults(t,t1,t2,s0,s1);
	_ff_neg(t2,t0); _ff_mult(u,t1,s0);										// negate t2 when setting it here to save doing it below
	_ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);
	
	// t2,t1,t0 is now f mod g.  reduce g =  1,s2,s1,s0 mod -t2,t1,t0
	_ff_mult(t,t2,s2);  _ff_add(s2,t,t1);
	_ff_mult(t,t2,s1);  _ff_add(s1,t,t0);  _ff_sum_2_mults(h[1],s1,s2,t1,t2);
	_ff_mult(t,t2,s0); _ff_sum_2_mults(h[0],t,s2,t0,t2);
	// total (3d-2)M  (using d+4 redc's)
}


static inline void ff_poly_mgcd_linear_12_3 (ff_t h[2], ff_t f[12], ff_t g[3])				// f monic deg 8, g monic deg 3, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, t0, t1, t2, t, u;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]);
	_ff_sub(t2,f[11],s2);
	_ff_sub(t1,f[10],s1); _ff_mult(t,t2,s2); _ff_subfrom(t1,t);
	_ff_sub(t0,f[9],s0); _ff_sum_2_mults(t,t1,t2,s1,s2); _ff_subfrom(t0,t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[8],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[7],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[6],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[5],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[4],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[3],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[2],t);
	_ff_sum_2_mults(t,t1,t2,s0,s1);
	_ff_neg(t2,t0); _ff_mult(u,t1,s0);										// negate t2 when setting it here to save doing it below
	_ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);
	
	// t2,t1,t0 is now f mod g.  reduce g =  1,s2,s1,s0 mod -t2,t1,t0
	_ff_mult(t,t2,s2);  _ff_add(s2,t,t1);
	_ff_mult(t,t2,s1);  _ff_add(s1,t,t0);  _ff_sum_2_mults(h[1],s1,s2,t1,t2);
	_ff_mult(t,t2,s0); _ff_sum_2_mults(h[0],t,s2,t0,t2);
	// total (3d-2)M  (using d+4 redc's)
}


static inline void ff_poly_mgcd_linear_14_3 (ff_t h[2], ff_t f[14], ff_t g[3])				// f monic deg 8, g monic deg 3, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, t0, t1, t2, t, u;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]);
	_ff_sub(t2,f[13],s2);
	_ff_sub(t1,f[11],s1); _ff_mult(t,t2,s2); _ff_subfrom(t1,t);
	_ff_sub(t0,f[10],s0); _ff_sum_2_mults(t,t1,t2,s1,s2); _ff_subfrom(t0,t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[10],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[9],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[8],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[7],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[6],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[5],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[4],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[3],t);
	_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
	_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[2],t);
	_ff_sum_2_mults(t,t1,t2,s0,s1);
	_ff_neg(t2,t0); _ff_mult(u,t1,s0);										// negate t2 when setting it here to save doing it below
	_ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);
	
	// t2,t1,t0 is now f mod g.  reduce g =  1,s2,s1,s0 mod -t2,t1,t0
	_ff_mult(t,t2,s2);  _ff_add(s2,t,t1);
	_ff_mult(t,t2,s1);  _ff_add(s1,t,t0);  _ff_sum_2_mults(h[1],s1,s2,t1,t2);
	_ff_mult(t,t2,s0); _ff_sum_2_mults(h[0],t,s2,t0,t2);
	// total (3d-2)M  (using d+4 redc's)
}


static inline void ff_poly_mgcd_linear_n_3 (ff_t h[2], ff_t *f, int d_f, ff_t g[3])		// f monic deg d_f >= 3, g monic deg 3, f and g are not modified, h may overlap either
{
	register int i;
	register ff_t s0, s1, s2, t0, t1, t2, t, u;
	
	switch(d_f) {
	case 3: ff_poly_mgcd_linear_3_3(h,f,g); return;
	case 4: ff_poly_mgcd_linear_4_3(h,f,g); return;
	}
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]);
	_ff_sub(t2,f[d_f-1],s2);
	_ff_sub(t1,f[d_f-2],s1); _ff_mult(t,t2,s2); _ff_subfrom(t1,t);
	_ff_sub(t0,f[d_f-3],s0); _ff_sum_2_mults(t,t1,t2,s1,s2); _ff_subfrom(t0,t);
	for ( i = d_f-4 ; i > 1 ; i-- ) {
		_ff_sum_3_mults(t,t0,t1,t2,s0,s1,s2);
		_ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[i],t);
	}
	_ff_sum_2_mults(t,t1,t2,s0,s1);
	_ff_neg(t2,t0); _ff_mult(u,t1,s0);										// negate t2 when setting it here to save doing it below
	_ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);
	
	// t2,t1,t0 is now f mod g.  reduce g =  1,s2,s1,s0 mod -t2,t1,t0
	_ff_mult(t,t2,s2);  _ff_add(s2,t,t1);
	_ff_mult(t,t2,s1);  _ff_add(s1,t,t0);  _ff_sum_2_mults(h[1],s1,s2,t1,t2);
	_ff_mult(t,t2,s0); _ff_sum_2_mults(h[0],t,s2,t0,t2);
	// total (3d-2)M  (using d+4 redc's)
}

static inline void ff_poly_gcd_linear_3_2_reg (ff_t h[2], ff_t f0, ff_t f1, ff_t f2, ff_t f3, ff_t g0, ff_t g1, ff_t g2)					// f and g need not be monic
{
	register ff_t w1, w2, w3;

	_ff_neg(w3,g2);  _ff_sum_2_mults(w1,f2,f3,g1,w3); _ff_mult(w2,f3,w3); _ff_square(w3,g2);
	_ff_sum_3_mults(h[1],f1,g0,g1,w1,w2,w3);
	_ff_sum_2_mults(h[0],f0,g0,w1,w3);
	// 9M+5A (5 redc)
}

static inline void ff_poly_gcd_linear_4_3_reg (ff_t h[2], ff_t f0, ff_t f1, ff_t f2, ff_t f3, ff_t f4, ff_t g0, ff_t g1, ff_t g2, ff_t g3)		// f and g need not be monic
{
	register ff_t w1, w2, w3;
	
	_ff_neg(w3,g3);  _ff_sum_2_mults(w1,f3,f4,g2,w3); _ff_mult(w2,f4,w3); _ff_square(w3,g3);
	_ff_sum_3_mults(f2,f2,g1,g2,w1,w2,w3);
	_ff_sum_3_mults(f1,f1,g0,g1,w1,w2,w3);
	_ff_sum_2_mults(f0,f0,g0,w1,w3);
	ff_poly_gcd_linear_3_2_reg(h,g0,g1,g2,g3,f0,f1,f2);
	// Total of 21M+12A (11 redc)
}

static inline void ff_poly_gcd_linear_5_4_reg (ff_t h[2], ff_t f0, ff_t f1, ff_t f2, ff_t f3, ff_t f4, ff_t f5, ff_t g0, ff_t g1, ff_t g2, ff_t g3, ff_t g4)			// f and g need not be monic
{
	register ff_t w1, w2, w3;
	
	_ff_neg(w3,g4);  _ff_sum_2_mults(w1,f4,f5,g3,w3); _ff_mult(w2,f5,w3); _ff_square(w3,g4);
	_ff_sum_3_mults(f3,f3,g2,g3,w1,w2,w3);
	_ff_sum_3_mults(f2,f2,g1,g2,w1,w2,w3);
	_ff_sum_3_mults(f1,f1,g0,g1,w1,w2,w3);
	_ff_sum_2_mults(f0,f0,g0,w1,w3);
	ff_poly_gcd_linear_4_3_reg(h,g0,g1,g2,g3,g4,f0,f1,f2,f3);
	// Total of 36M+21A (18 redc)
}

static inline void ff_poly_gcd_linear_6_5_reg (ff_t h[2], ff_t f0, ff_t f1, ff_t f2, ff_t f3, ff_t f4, ff_t f5, ff_t f6, ff_t g0, ff_t g1, ff_t g2, ff_t g3, ff_t g4, ff_t g5)			// f and g need not be monic
{
	register ff_t w1, w2, w3;
	
	_ff_neg(w3,g5);  _ff_sum_2_mults(w1,f5,f6,g4,w3); _ff_mult(w2,f6,w3); _ff_square(w3,g5);
	_ff_sum_3_mults(f4,f4,g3,g4,w1,w2,w3);
	_ff_sum_3_mults(f3,f3,g2,g3,w1,w2,w3);
	_ff_sum_3_mults(f2,f2,g1,g2,w1,w2,w3);
	_ff_sum_3_mults(f1,f1,g0,g1,w1,w2,w3);
	_ff_sum_2_mults(f0,f0,g0,w1,w3);
	ff_poly_gcd_linear_5_4_reg(h,g0,g1,g2,g3,g4,g5,f0,f1,f2,f3,f4);
	// Total of 54M+32A (26 redc)
}

static inline void ff_poly_gcd_linear_7_6_reg (ff_t h[2], ff_t f0, ff_t f1, ff_t f2, ff_t f3, ff_t f4, ff_t f5, ff_t f6, ff_t f7, ff_t g0, ff_t g1, ff_t g2, ff_t g3, ff_t g4, ff_t g5, ff_t g6)			// f and g need not be monic
{
	register ff_t w1, w2, w3;
	
	_ff_neg(w3,g6);  _ff_sum_2_mults(w1,f6,f7,g5,w3); _ff_mult(w2,f7,w3); _ff_square(w3,g6);
	_ff_sum_3_mults(f5,f5,g4,g5,w1,w2,w3);
	_ff_sum_3_mults(f4,f4,g3,g4,w1,w2,w3);
	_ff_sum_3_mults(f3,f3,g2,g3,w1,w2,w3);
	_ff_sum_3_mults(f2,f2,g1,g2,w1,w2,w3);
	_ff_sum_3_mults(f1,f1,g0,g1,w1,w2,w3);
	_ff_sum_2_mults(f0,f0,g0,w1,w3);
	ff_poly_gcd_linear_6_5_reg(h,g0,g1,g2,g3,g4,g5,g6,f0,f1,f2,f3,f4,f5);
	// Total of 75M+45A (35 redc)
}

static inline void ff_poly_gcd_linear_8_7_reg (ff_t h[2], ff_t f0, ff_t f1, ff_t f2, ff_t f3, ff_t f4, ff_t f5, ff_t f6, ff_t f7, ff_t f8, ff_t g0, ff_t g1, ff_t g2, ff_t g3, ff_t g4, ff_t g5, ff_t g6, ff_t g7)			// f and g need not be monic
{
	register ff_t w1, w2, w3;
	
	_ff_neg(w3,g6);  _ff_sum_2_mults(w1,f7,f8,g6,w3); _ff_mult(w2,f8,w3); _ff_square(w3,g7);
	_ff_sum_3_mults(f6,f6,g5,g6,w1,w2,w3);
	_ff_sum_3_mults(f5,f5,g4,g5,w1,w2,w3);
	_ff_sum_3_mults(f4,f4,g3,g4,w1,w2,w3);
	_ff_sum_3_mults(f3,f3,g2,g3,w1,w2,w3);
	_ff_sum_3_mults(f2,f2,g1,g2,w1,w2,w3);
	_ff_sum_3_mults(f1,f1,g0,g1,w1,w2,w3);
	_ff_sum_2_mults(f0,f0,g0,w1,w3);
	ff_poly_gcd_linear_6_5_reg(h,g0,g1,g2,g3,g4,g5,g6,f0,f1,f2,f3,f4,f5);
	// Total of 99M+60A (45 redc)
}

// this is currently not used and hasn't been tested
static inline void ff_poly_gcd_linear_4_3 (ff_t h[2], ff_t f[4], ff_t g[4])				// f deg 4 monic, g deg 3 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t t0, t1, t2, s0, ng0, ng1, ng2, g3, g32;

	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng2,g[2]); _ff_set(g3,g[3]); _ff_square(g32,g3);
	_ff_mult(t0,f[3],g3); _ff_addto(t0,ng2);  _ff_mult(s0,t0,ng2);
	_ff_sum_2_mults(t2,f[2],ng1,g3,g32);  _ff_addto(t2,s0);
	_ff_sum_3_mults(t1,f[1],t1,s0,ng0,ng1,g32);
	_ff_sum_2_mults(t0,f[0],t2,ng0,g32);
	ff_poly_gcd_linear_3_2_reg(h,g[0],g[1],g[2],g[3],t0,t1,t2);
	// Total of 19M+14A (11 redc) (could use fewer additions)
}

static inline void ff_poly_gcd_linear_6_3 (ff_t h[2], ff_t f[6], ff_t g[4])				// f deg 6 monic, g deg 3 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, g3, ng0, ng1, t0, t1, t2, t3;
	
	_ff_set(g3,g[3]); _ff_neg(s2,g[2]); _ff_neg(ng1,g[1]); _ff_mult(s1,ng1,g3); _ff_square(t0,g3); _ff_neg(ng0,g[0]); _ff_mult(s0,t0,ng0);
	_ff_mult(t3,f[5],g3); _ff_addto(t3,s2);
	_ff_sum_2_mults(t2,f[4],t3,s2,t0); _ff_addto(t2,s1);
	ff_mult(t0,t0,g3);
	_ff_sum_3_mults(t1,f[3],t2,t3,s1,s2,t0);  _ff_addto(t1,s0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t3,f[2],t1,t2,t3,s0,s1,s2,t0);
	_ff_mult(s0,g3,ng0);
	_ff_sum_3_mults(t2,f[1],t1,t2,s0,ng1,t0);
	_ff_sum_2_mults(t1,f[0],t1,ng0,t0);
	// 21M + 15A (12 redc)
	ff_poly_gcd_linear_3_2_reg(h,g[0],g[1],g[2],g[3],t1,t2,t3);
	// 30M  + 20A (17 redc)
}

static inline void ff_poly_gcd_linear_8_3 (ff_t h[2], ff_t f[8], ff_t g[4])				// f deg 8 monic, g deg 3 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, g3, ng0, ng1, t0, t1, t2, t3;
	
	_ff_set(g3,g[3]); _ff_neg(s2,g[2]); _ff_neg(ng1,g[1]); _ff_mult(s1,ng1,g3); _ff_square(t0,g3); _ff_neg(ng0,g[0]); _ff_mult(s0,t0,ng0);
	_ff_mult(t3,f[7],g3); _ff_addto(t3,s2);
	_ff_sum_2_mults(t2,f[6],t3,s2,t0); _ff_addto(t2,s1);
	ff_mult(t0,t0,g3);
	_ff_sum_3_mults(t1,f[5],t2,t3,s1,s2,t0);  _ff_addto(t1,s0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t3,f[4],t1,t2,t3,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t2,f[3],t3,t1,t2,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t1,f[2],t2,t3,t1,s0,s1,s2,t0);
	_ff_mult(s0,g3,ng0);
	_ff_sum_3_mults(t3,f[1],t2,t3,s0,ng1,t0);
	_ff_sum_2_mults(t2,f[0],t2,ng0,t0);
	ff_poly_gcd_linear_3_2_reg(h,g[0],g[1],g[2],g[3],t2,t3,t1);
	// 39M  + 26A (20 redc)
}

static inline void ff_poly_gcd_linear_12_3 (ff_t h[2], ff_t f[12], ff_t g[4])			// f deg 12 monic, g deg 3 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, g3, ng0, ng1, t0, t1, t2, t3;
	
	_ff_set(g3,g[3]); _ff_neg(s2,g[2]); _ff_neg(ng1,g[1]); _ff_mult(s1,ng1,g3); _ff_square(t0,g3); _ff_neg(ng0,g[0]); _ff_mult(s0,t0,ng0);
	_ff_mult(t3,f[11],g3); _ff_addto(t3,s2);
	_ff_sum_2_mults(t2,f[10],t3,s2,t0); _ff_addto(t2,s1);
	ff_mult(t0,t0,g3);
	_ff_sum_3_mults(t1,f[9],t2,t3,s1,s2,t0);  _ff_addto(t1,s0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t3,f[8],t1,t2,t3,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t2,f[7],t3,t1,t2,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t1,f[6],t2,t3,t1,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t3,f[5],t1,t2,t3,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t2,f[4],t3,t1,t2,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t1,f[3],t2,t3,t1,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t3,f[2],t1,t2,t3,s0,s1,s2,t0);
	_ff_mult(s0,g3,ng0);
	_ff_sum_3_mults(t2,f[1],t1,t2,s0,ng1,t0);
	_ff_sum_2_mults(t1,f[0],t1,ng0,t0);
	ff_poly_gcd_linear_3_2_reg(h,g[0],g[1],g[2],g[3],t1,t2,t3);
}

static inline void ff_poly_gcd_linear_14_3 (ff_t h[2], ff_t f[14], ff_t g[4])			// f deg 14 monic, g deg 3 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, g3, ng0, ng1, t0, t1, t2, t3;
	
	_ff_set(g3,g[3]); _ff_neg(s2,g[2]); _ff_neg(ng1,g[1]); _ff_mult(s1,ng1,g3); _ff_square(t0,g3); _ff_neg(ng0,g[0]); _ff_mult(s0,t0,ng0);
	_ff_mult(t3,f[13],g3); _ff_addto(t3,s2);
	_ff_sum_2_mults(t2,f[12],t3,s2,t0); _ff_addto(t2,s1);
	ff_mult(t0,t0,g3);
	_ff_sum_3_mults(t1,f[11],t2,t3,s1,s2,t0);  _ff_addto(t1,s0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t3,f[10],t1,t2,t3,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t2,f[9],t3,t1,t2,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t1,f[8],t2,t3,t1,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t3,f[7],t1,t2,t3,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t2,f[6],t3,t1,t2,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t1,f[5],t2,t3,t1,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t3,f[4],t1,t2,t3,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t2,f[3],t3,t1,t2,s0,s1,s2,t0);
	_ff_mult(t0,t0,g3);
	_ff_sum_4_mults(t1,f[2],t2,t3,t1,s0,s1,s2,t0);
	_ff_mult(s0,g3,ng0);
	_ff_sum_3_mults(t3,f[1],t2,t3,s0,ng1,t0);
	_ff_sum_2_mults(t2,f[0],t2,ng0,t0);
	ff_poly_gcd_linear_3_2_reg(h,g[0],g[1],g[2],g[3],t2,t3,t1);
}

static inline void ff_poly_gcd_linear_n_3 (ff_t h[2], ff_t *f, int d_f, ff_t g[4])			// f deg > 4 monic, g deg 3 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t s0, s1, s2, g3, ng0, ng1, t0, t1, t2, t3, t;
	register int i;
	
	_ff_set(g3,g[3]); _ff_neg(s2,g[2]); _ff_neg(ng1,g[1]); _ff_mult(s1,ng1,g3); _ff_square(t0,g3); _ff_neg(ng0,g[0]); _ff_mult(s0,t0,ng0);
	_ff_mult(t3,f[d_f-1],g3); _ff_addto(t3,s2);
	_ff_sum_2_mults(t2,f[d_f-2],t3,s2,t0); _ff_addto(t2,s1);
	ff_mult(t0,t0,g3);
	_ff_sum_3_mults(t1,f[d_f-3],t2,t3,s1,s2,t0);  _ff_addto(t1,s0);
	for ( i = d_f-4 ; i > 1 ; i-- ) {
		_ff_mult(t0,t0,g3);
		_ff_sum_4_mults(t,f[i],t1,t2,t3,s0,s1,s2,t0);
		_ff_set(t3,t2); _ff_set(t2,t1); _ff_set(t1,t);
	}
	_ff_mult(s0,g3,ng0);
	_ff_sum_3_mults(t3,f[1],t2,t3,s0,ng1,t0);
	_ff_sum_2_mults(t2,f[0],t2,ng0,t0);
	ff_poly_gcd_linear_3_2_reg(h,g[0],g[1],g[2],g[3],t2,t3,t1);
}

static inline void ff_poly_gcd_linear_6_4 (ff_t h[2], ff_t f[6], ff_t g[5])				// f deg 6 monic, g deg 4 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t ng0, ng1, ng2, ng3, g4, g42, g43, s1, t0, t1, t2, t3;

	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng2,g[2]); _ff_neg(ng3,g[3]);_ff_set(g4,g[4]); _ff_square(g42,g4); _ff_mult(g43,g42,g4); 
	_ff_mult(t1,f[5],g4); _ff_addto(t1,ng3);
	_ff_sum_3_mults(t0,f[4],t1,g4,ng2,ng3,g42);
	_ff_mult(s1,g4,t1);
	_ff_sum_4_mults(t3,f[3],t0,s1,g42,ng1,ng2,ng3,g43);
	_ff_sum_4_mults(t2,f[2],t0,s1,g42,ng0,ng1,ng2,g43);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,g43);
	_ff_sum_2_mults(t0,f[0],t0,ng0,g43);
	// 20M + 16A (9 redc)
	ff_poly_gcd_linear_4_3_reg(h,g[0],g[1],g[2],g[3],g[4],t0,t1,t2,t3);
	// 41M  + 28A (20 redc)
}

static inline void ff_poly_gcd_linear_8_5 (ff_t h[2], ff_t f[8], ff_t g[6])				// f deg 8 monic, g deg 5 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t ng0, ng1, ng2, ng3, ng4, g5, g52, g53, g54, t0, t1, t2, t3, t4, s1, s2, s3;
	
	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng2,g[2]); _ff_neg(ng3,g[3]); _ff_neg(ng4,g[4]); _ff_set(g5,g[5]);
	_ff_square(g52,g5); _ff_mult(g53,g52,g5); _ff_square(g54,g52); _ff_mult(s3,ng3,g5); _ff_mult(s2,ng2,g52); _ff_mult(s1,ng1,g53);
	_ff_mult(t2,f[7],g5); _ff_addto(t2,ng4);
	_ff_sum_2_mults(t1,f[6],t2,ng4,g52); _ff_addto(t1,s3);
	_ff_sum_3_mults(t0,f[5],t1,t2,s3,ng4,g53);  _ff_addto(t0,s2);
	_ff_sum_4_mults(t4,f[4],t0,t1,t2,s2,s3,ng4,g54); _ff_addto(t4,s1);
	_ff_mult(s2,g52,t2); _ff_mult(s1,g5,t1);
	_ff_sum_5_mults(t3,f[3],t0,s1,s2,g53,ng0,ng1,ng2,ng3,g54);
	_ff_sum_4_mults(t2,f[2],t0,s1,s2,ng0,ng1,ng2,g54);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,g54);
	_ff_sum_2_mults(t0,f[0],t0,ng0,g54);
	// 32M + 25A (16 redc)
	ff_poly_gcd_linear_5_4_reg(h,g[0],g[1],g[2],g[3],g[4],g[5],t0,t1,t2,t3,t4);
	// 68M  + 50A (34 redc)
}

static inline void ff_poly_gcd_linear_12_5 (ff_t h[2], ff_t f[12], ff_t g[6])				// f deg 12 monic, g deg 5 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t ng0, ng1, ng2, ng3, ng4, g5, g52, g53, g54, w, t0, t1, t2, t3, t4, s0, s1, s2, s3;

	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng2,g[2]); _ff_neg(ng3,g[3]); _ff_neg(ng4,g[4]); _ff_set(g5,g[5]);
	_ff_square(g52,g5); _ff_mult(g53,g52,g5); _ff_square(g54,g52); _ff_mult(s3,ng3,g5); _ff_mult(s2,ng2,g52); _ff_mult(s1,ng1,g53); _ff_mult(s0,ng0,g54);
	_ff_mult(t1,f[11],g5); _ff_addto(t1,ng4);
	_ff_sum_2_mults(t0,f[10],t1,ng4,g52); _ff_addto(t0,s3);
	_ff_sum_3_mults(t4,f[9],t0,t1,s3,ng4,g53);  _ff_addto(t4,s2);
	_ff_sum_4_mults(t3,f[8],t4,t0,t1,s2,s3,ng4,g54); _ff_addto(t3,s1);
	_ff_mult(w,g54,g5);
	_ff_sum_5_mults(t2,f[7],t3,t4,t0,t1,s1,s2,s3,ng4,w);  _ff_addto(t2,s0);
	ff_mult(w,w,g5);
	_ff_sum_6_mults(t1,f[6],t2,t3,t4,t0,t1,s0,s1,s2,s3,ng4,w);
	ff_mult(w,w,g5);
	_ff_sum_6_mults(t0,f[5],t1,t2,t3,t4,t0,s0,s1,s2,s3,ng4,w);
	ff_mult(w,w,g5);
	_ff_sum_6_mults(t4,f[4],t0,t1,t2,t3,t4,s0,s1,s2,s3,ng4,w);
	_ff_mult(s3,g53,t3); _ff_mult(s2,g52,t2); _ff_mult(s1,g5,t1);
	_ff_sum_5_mults(t3,f[3],t0,s1,s2,s3,ng0,ng1,ng2,ng3,w);
	_ff_sum_4_mults(t2,f[2],t0,s1,s2,ng0,ng1,ng2,w);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,w);
	_ff_sum_2_mults(t0,f[0],t0,ng0,w);
	// 61M + 26A (26 redc)
	ff_poly_gcd_linear_5_4_reg(h,g[0],g[1],g[2],g[3],g[4],g[5],t0,t1,t2,t3,t4);
	// 97M + 47A (44 redc)
}

static inline void ff_poly_gcd_linear_14_5 (ff_t h[2], ff_t f[14], ff_t g[6])				// f deg 12 monic, g deg 5 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t ng0, ng1, ng2, ng3, ng4, g5, g52, g53, g54, w, t0, t1, t2, t3, t4, s0, s1, s2, s3;

	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng2,g[2]); _ff_neg(ng3,g[3]); _ff_neg(ng4,g[4]); _ff_set(g5,g[5]);
	_ff_square(g52,g5); _ff_mult(g53,g52,g5); _ff_square(g54,g52); _ff_mult(s3,ng3,g5); _ff_mult(s2,ng2,g52); _ff_mult(s1,ng1,g53); _ff_mult(s0,ng0,g54);
	_ff_mult(t3,f[13],g5); _ff_addto(t3,ng4);
	_ff_sum_2_mults(t2,f[12],t3,ng4,g52); _ff_addto(t2,s3);
	_ff_sum_3_mults(t1,f[11],t2,t3,s3,ng4,g53);  _ff_addto(t1,s2);
	_ff_sum_4_mults(t0,f[10],t1,t2,t3,s2,s3,ng4,g54); _ff_addto(t0,s1);
	_ff_mult(w,g54,g5);
	_ff_sum_5_mults(t4,f[9],t0,t1,t2,t3,s1,s2,s3,ng4,w);  _ff_addto(t4,s0);
	ff_mult(w,w,g5);
	_ff_sum_6_mults(t3,f[8],t4,t0,t1,t2,t3,s0,s1,s2,s3,ng4,w);
	ff_mult(w,w,g5);
	_ff_sum_6_mults(t2,f[7],t3,t4,t0,t1,t2,s0,s1,s2,s3,ng4,w);
	ff_mult(w,w,g5);
	_ff_sum_6_mults(t1,f[6],t2,t3,t4,t0,t1,s0,s1,s2,s3,ng4,w);
	ff_mult(w,w,g5);
	_ff_sum_6_mults(t0,f[5],t1,t2,t3,t4,t0,s0,s1,s2,s3,ng4,w);
	ff_mult(w,w,g5);
	_ff_sum_6_mults(t4,f[4],t0,t1,t2,t3,t4,s0,s1,s2,s3,ng4,w);
	_ff_mult(s3,g53,t3); _ff_mult(s2,g52,t2); _ff_mult(s1,g5,t1);
	_ff_sum_5_mults(t3,f[3],t0,s1,s2,s3,ng0,ng1,ng2,ng3,w);
	_ff_sum_4_mults(t2,f[2],t0,s1,s2,ng0,ng1,ng2,w);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,w);
	_ff_sum_2_mults(t0,f[0],t0,ng0,w);
	// 75M + 36A (30 redc)
	ff_poly_gcd_linear_5_4_reg(h,g[0],g[1],g[2],g[3],g[4],g[5],t0,t1,t2,t3,t4);
	// 111M  + 51A (48 redc)
}

static inline void ff_poly_gcd_linear_n_5 (ff_t h[2], ff_t *f, int d_f, ff_t g[6])				// f deg > 10 monic, g deg 5 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t ng0, ng1, ng2, ng3, ng4, g5, g52, g53, g54, w, t, t0, t1, t2, t3, t4, s0, s1, s2, s3;
	register int i;

	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng2,g[2]); _ff_neg(ng3,g[3]); _ff_neg(ng4,g[4]); _ff_set(g5,g[5]);
	_ff_square(g52,g5); _ff_mult(g53,g5,g52); _ff_square(g54,g52);
	_ff_mult(s0,ng0,g54); _ff_mult(s1,ng1,g53); _ff_mult(s2,ng2,g52); _ff_mult(s3,ng3,g5);
	_ff_mult(t3,f[d_f-1],g5); _ff_addto(t3,ng4);
	_ff_sum_2_mults(t2,f[d_f-2],t3,ng4,g52); _ff_addto(t2,s3);
	_ff_sum_3_mults(t1,f[d_f-3],t2,t3,s3,ng4,g53); _ff_addto(t1,s2);
	_ff_sum_4_mults(t0,f[d_f-4],t1,t2,t3,s2,s3,ng4,g54); _ff_addto(t0, s1);
	ff_mult(w,g54,g5);
	_ff_sum_5_mults(t4,f[d_f-5],t0,t1,t2,t3,s1,s2,s3,ng4,w); _ff_addto(t4, s0);
	for ( i = d_f-6 ; i > 3 ; i-- ) {
		ff_mult(w,w,g5);
		_ff_sum_6_mults(t,f[i],t4,t0,t1,t2,t3,s0,s1,s2,s3,ng4,w);
		_ff_set(t3,t2); _ff_set(t2,t1); _ff_set(t1,t0); _ff_set(t0,t4); _ff_set(t4,t);
	}
	_ff_mult(s3,g53,t3); _ff_mult(s2,g52,t2); _ff_mult(s1,g5,t1);
	_ff_sum_5_mults(t3,f[3],t0,s1,s2,s3,ng0,ng1,ng2,ng3,w);
	_ff_sum_4_mults(t2,f[2],t0,s1,s2,ng0,ng1,ng2,w);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,w);
	_ff_sum_2_mults(t0,f[0],t0,ng0,w);
	// (7n-24)M (2n+2)redc
	ff_poly_gcd_linear_5_4_reg(h,g[0],g[1],g[2],g[3],g[4],g[5],t0,t1,t2,t3,t4);
	// (7n+12)M (2n+20)redc
}

static inline void ff_poly_gcd_linear_8_6 (ff_t h[2], ff_t f[8], ff_t g[7])				// f deg 8 monic, g deg 6 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t ng0, ng1, ng2, ng3, ng4, ng5, g6, g62, g63, t0, t1, t2, t3, t4, t5, s1;

	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng2,g[2]); _ff_neg(ng3,g[3]); _ff_neg(ng4,g[4]); _ff_neg(ng5,g[5]); _ff_set(g6,g[6]); 
	_ff_square(g62,g6); _ff_mult(g63,g62,g6);
	_ff_mult(t1,f[7],g6); _ff_addto(t1,ng5);
	_ff_mult(s1,g6,t1);
	_ff_sum_3_mults(t0,f[6],t1,g6,ng4,ng5,g62);
	_ff_sum_4_mults(t5,f[5],t0,s1,g62,ng3,ng4,ng5,g63);
	_ff_sum_4_mults(t4,f[4],t0,s1,g62,ng2,ng3,ng4,g63);
	_ff_sum_4_mults(t3,f[3],t0,s1,g62,ng1,ng2,ng3,g63);
	_ff_sum_4_mults(t2,f[2],t0,s1,g62,ng0,ng1,ng2,g63);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,g63);
	_ff_sum_2_mults(t0,f[0],t0,ng0,g63);
	// 28M + 24A (11 redc)
	ff_poly_gcd_linear_6_5_reg(h,g[0],g[1],g[2],g[3],g[4],g[5],g[6],t0,t1,t2,t3,t4,t5);
	// 82M + 56A (35 redc)
}

static inline void ff_poly_gcd_linear_12_6 (ff_t h[2], ff_t f[12], ff_t g[7])			// f deg 12 monic, g deg 6 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t ng0, ng1, ng2, ng3, ng4, ng5, g6, g62, g63, g64, g65, g66, g67, t0, t1, t2, t3, t4, t5, s1, s2, s3, s4, s5;
	
	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng2,g[2]); _ff_neg(ng3,g[3]); _ff_neg(ng4,g[4]); _ff_neg(ng5,g[5]); _ff_set(g6,g[6]); 
	_ff_square(g62,g6); _ff_mult(g63,g62,g6); _ff_square(g64,g62); _ff_mult(g65,g64,g6); _ff_square(g66,g63); _ff_mult(g67,g66,g6);
	_ff_mult(s1,g6,ng4); _ff_mult(s2,g62,ng3); _ff_mult(s3,g63,ng2); _ff_mult(s4,g64,ng1); _ff_mult(s5,g65,ng0);
	_ff_mult(t5,f[11],g6); _ff_addto(t5,ng5);
	_ff_sum_2_mults(t4,f[10],t5,ng5,g62); _ff_addto(t4,s1);
	_ff_sum_3_mults(t3,f[9],t4,t5,s1,ng5,g63); _ff_addto(t3,s2);
	_ff_sum_4_mults(t2,f[8],t3,t4,t5,s2,s1,ng5,g64); _ff_addto(t2,s3);
	_ff_sum_5_mults(t1,f[7],t2,t3,t4,t5,s3,s2,s1,ng5,g65); _ff_addto(t1,s4);
	_ff_sum_6_mults(t0,f[6],t1,t2,t3,t4,t5,s4,s3,s2,s1,ng5,g66); _ff_addto(t0,s5);
	_ff_mult(s1,g6,t1); _ff_mult(s2,g62,t2); _ff_mult(s3,g63,t3); _ff_mult(s4,g64,t4); 
	_ff_sum_7_mults(t5,f[5],t0,s1,s2,s3,s4,s5,t5,ng1,ng2,ng3,ng4,ng5,g67);
	_ff_sum_6_mults(t4,f[4],t0,s1,s2,s3,s4,ng0,ng1,ng2,ng3,ng4,g67);
	_ff_sum_5_mults(t3,f[3],t0,s1,s2,s3,ng0,ng1,ng2,ng3,g67);
	_ff_sum_4_mults(t2,f[2],t0,s1,s2,ng0,ng1,ng2,g67);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,g67);
	_ff_sum_2_mults(t0,f[0],t0,ng0,g67);
	// 63M + 49A (27 redc)
	ff_poly_gcd_linear_6_5_reg(h,g[0],g[1],g[2],g[3],g[4],g[5],g[6],t0,t1,t2,t3,t4,t5);
	// 117M + 81A (53 redc)
}

static inline void ff_poly_gcd_linear_14_6 (ff_t h[2], ff_t f[12], ff_t g[7])			// f deg 14 monic, g deg 6 not nesc monic, f and g are not modified, h may overlap either
{
	register ff_t ng0, ng1, ng2, ng3, ng4, ng5, g6, g62, g63, g64, g65, w, t0, t1, t2, t3, t4, t5, s1, s2, s3, s4, s5;
	
	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng2,g[2]); _ff_neg(ng3,g[3]); _ff_neg(ng4,g[4]); _ff_neg(ng5,g[5]); _ff_set(g6,g[6]); 
	_ff_square(g62,g6); _ff_mult(g63,g62,g6); _ff_square(g64,g62); _ff_mult(g65,g64,g6);
	_ff_mult(s1,g6,ng4); _ff_mult(s2,g62,ng3); _ff_mult(s3,g63,ng2); _ff_mult(s4,g64,ng1); _ff_mult(s5,g65,ng0);
	_ff_mult(t1,f[13],g6); _ff_addto(t1,ng5);
	_ff_sum_2_mults(t0,f[12],t1,ng5,g62); _ff_addto(t0,s1);
	_ff_sum_3_mults(t5,f[11],t0,t1,s1,ng5,g63); _ff_addto(t5,s2);
	_ff_sum_4_mults(t4,f[10],t5,t0,t1,s2,s1,ng5,g64); _ff_addto(t4,s3);
	_ff_sum_5_mults(t3,f[9],t4,t5,t0,t1,s3,s2,s1,ng5,g65); _ff_addto(t3,s4);
	ff_mult(w,g65,g6);
	_ff_sum_6_mults(t2,f[8],t3,t4,t5,t0,t1,s4,s3,s2,s1,ng5,w); _ff_addto(t2,s5);
	ff_mult(w,w,g6);
	_ff_sum_7_mults(t1,f[7],t2,t3,t4,t5,t0,t1,s5,s4,s3,s2,s1,ng5,w);
	ff_mult(w,w,g6);
	_ff_sum_7_mults(t0,f[6],t1,t2,t3,t4,t5,t0,s5,s4,s3,s2,s1,ng5,w);
	ff_mult(w,w,g6);
	_ff_mult(s1,g6,t1); _ff_mult(s2,g62,t2); _ff_mult(s3,g63,t3); _ff_mult(s4,g64,t4); 
	_ff_sum_7_mults(t5,f[5],t0,s1,s2,s3,s4,s5,t5,ng1,ng2,ng3,ng4,ng5,w);
	_ff_sum_6_mults(t4,f[4],t0,s1,s2,s3,s4,ng0,ng1,ng2,ng3,ng4,w);
	_ff_sum_5_mults(t3,f[3],t0,s1,s2,s3,ng0,ng1,ng2,ng3,w);
	_ff_sum_4_mults(t2,f[2],t0,s1,s2,ng0,ng1,ng2,w);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,w);
	_ff_sum_2_mults(t0,f[0],t0,ng0,w);
	// 79M + 61A (31 redc)
	ff_poly_gcd_linear_6_5_reg(h,g[0],g[1],g[2],g[3],g[4],g[5],g[6],t0,t1,t2,t3,t4,t5);
	// 133M + 93A (57 redc)
}

static inline void ff_poly_mgcd_linear_5_5 (ff_t h[2], ff_t f[7], ff_t g[5])				// f monic deg 5, g monic deg 5
{
	register ff_t s0, s1, s2, s3, t0, t1, t2, t3, t4, w1, w2,w3;
	
	_ff_sub(t4,f[4],g[4]);  _ff_sub(t3,f[3],g[3]);  _ff_sub(t2,f[2],g[2]);  _ff_sub(t1,f[1],g[1]);  _ff_sub(t1,f[0],g[0]);  

	// reduce 1 s4 s3 s2 s1 s0 mod t4 t3 t2 t1 t0 here (saves 2 mults versus calling 5_4 code)
	_ff_neg(w1,t4);  _ff_mult(w2,w1,g[4]); _ff_addto(w2,t3);  _ff_square(w3,w1);
	_ff_sum_3_mults(s3,g[3],t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,g[2],t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,g[1],t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,g[0],t0,w2,w3);

	ff_poly_gcd_linear_4_3_reg(h,t0,t1,t2,t3,t4,s0,s1,s2,s3);
	// Total 44M+38A (23 redc)
}

static inline void ff_poly_mgcd_linear_m_4 (ff_t h[2], ff_t *f, int d_f, ff_t g[4])				// f monic deg 4 to 6, g monic deg 4
{
	register ff_t s0, s1, s2, s3, t0, t1, t2, t3, w1, w2, w3, t;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]);
	switch (d_f) {
	case 4:
		_ff_sub(t3,f[3],s3);  _ff_sub(t2,f[2],s2);  _ff_sub(t1,f[1],s1);  _ff_sub(t0,f[0],s0);  break;
	case 5:
		_ff_sub(w1,f[4],s3);
		_ff_sub(t3,f[3],s2); _ff_mult(t,w1,s3); _ff_subfrom(t3,t);
		_ff_sub(t2,f[2],s1); _ff_mult(t,w1,s2); _ff_subfrom(t2,t);
		_ff_sub(t1,f[1],s0); _ff_mult(t,w1,s1); _ff_subfrom(t1,t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 6:
		_ff_sub(w2,f[5],s3);
		_ff_sub(w1,f[4],s2); _ff_mult(t,w2,s3); _ff_subfrom(w1,t);
		_ff_sub(t3,f[3],s1); _ff_sum_2_mults(t,w1,w2,s2,s3); _ff_subfrom(t3,t);
		_ff_sub(t2,f[2],s0); _ff_sum_2_mults(t,w1,w2,s1,s2); _ff_subfrom(t2,t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	default:
		printf ("Bad degree %d in ff_poly_mgcd_linear_m_4\n", d_f); abort();
	}
	
	// reduce 1 s3 s2 s1 s0 mod t3 t2 t1 t0 here (saves 2 mults versus calling 4_3 code)
	_ff_neg(w1,t3);  _ff_mult(w2,w1,s3); _ff_addto(w2,t2);  _ff_square(w3,w1);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_3_2_reg(h,t0,t1,t2,t3,s0,s1,s2);
}


static inline void ff_poly_mgcd_linear_n_4 (ff_t h[2], ff_t *f, int d_f, ff_t g[4])				// f monic deg >= 4, g monic deg 4
{
	register int i;
	register ff_t s0, s1, s2, s3,  t0, t1, t2, t3, w1, w2, w3, t, u;
	
	if ( d_f < 7 ) { ff_poly_mgcd_linear_m_4(h,f,d_f,g); return; }
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]);
	_ff_sub(t3,f[d_f-1],s3);
	_ff_sub(t2,f[d_f-2],s2); _ff_mult(t,t3,s3); _ff_subfrom(t2,t);
	_ff_sub(t1,f[d_f-3],s1); _ff_sum_2_mults(t,t2,t3,s2,s3); _ff_subfrom(t1,t);
	_ff_sub(t0,f[d_f-4],s0); _ff_sum_3_mults(t,t1,t2,t3,s1,s2,s3); _ff_subfrom(t0,t);
	for ( i = d_f-5; i > 2 ; i-- ) {
		_ff_sum_4_mults(t,t0,t1,t2,t3,s0,s1,s2,s3);
		 _ff_set(t3,t2); _ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[i],t);
	}
	_ff_sum_3_mults(u,t1,t2,t3,s0,s1,s2);  _ff_set(t3,t0);
	_ff_sum_2_mults(t,t1,t2,s0,s1);  _ff_sub(t2,f[2],u);
	_ff_mult(u,t1,s0);  _ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);

	// reduce 1 s3 s2 s1 s0 mod t3 t2 t1 t0 here (saves 2 mults versus calling 4_3 code)
	_ff_neg(w1,t3);  _ff_mult(w2,w1,s3); _ff_addto(w2,t2);  _ff_square(w3,w1);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_3_2_reg(h,t0,t1,t2,t3,s0,s1,s2);
}


static inline void ff_poly_mgcd_linear_m_5 (ff_t h[2], ff_t *f, int d_f, ff_t g[5])				// f monic deg 5 to 8, g monic deg 5
{
	register ff_t s0, s1, s2, s3, s4, t0, t1, t2, t3, t4, w1, w2, w3, t;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]);
	switch (d_f) {
	case 5:
		_ff_sub(t4,f[4],s4);  _ff_sub(t3,f[3],s3);  _ff_sub(t2,f[2],s2);  _ff_sub(t1,f[1],s1);  _ff_sub(t0,f[0],s0);  break;
	case 6:
		_ff_sub(w1,f[5],s4);
		_ff_sub(t4,f[4],s3); _ff_mult(t,w1,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s2); _ff_mult(t,w1,s3); _ff_subfrom(t3,t);
		_ff_sub(t2,f[2],s1); _ff_mult(t,w1,s2); _ff_subfrom(t2,t);
		_ff_sub(t1,f[1],s0); _ff_mult(t,w1,s1); _ff_subfrom(t1,t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 7:
		_ff_sub(w2,f[6],s4);
		_ff_sub(w1,f[5],s3); _ff_mult(t,w2,s4); _ff_subfrom(w1,t);
		_ff_sub(t4,f[4],s2); _ff_sum_2_mults(t,w1,w2,s3,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s1); _ff_sum_2_mults(t,w1,w2,s2,s3); _ff_subfrom(t3,t);
		_ff_sub(t2,f[2],s0); _ff_sum_2_mults(t,w1,w2,s1,s2); _ff_subfrom(t2,t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 8:
		_ff_sub(w3,f[7],s4);
		_ff_sub(w2,f[6],s3); _ff_mult(t,w3,s4); _ff_subfrom(w2,t);
		_ff_sub(w1,f[5],s2); _ff_sum_2_mults(t,w2,w3,s3,s4); _ff_subfrom(w1,t);
		_ff_sub(t4,f[4],s1); _ff_sum_3_mults(t,w1,w2,w3,s2,s3,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s0); _ff_sum_3_mults(t,w1,w2,w3,s1,s2,s3); _ff_subfrom(t3,t);
		_ff_sum_3_mults(t,w1,w2,w3,s0,s1,s2); _ff_sub(t2,f[2],t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	default:
		printf ("Bad degree %d in ff_poly_mgcd_linear_m_5\n", d_f); abort();
	}
	
	// reduce 1 s4 s3 s2 s1 s0 mod t4 t3 t2 t1 t0 here (saves 2 mults versus calling 5_4 code)
	_ff_neg(w1,t4);  _ff_mult(w2,w1,s4); _ff_addto(w2,t3);  _ff_square(w3,w1);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_4_3_reg(h,t0,t1,t2,t3,t4,s0,s1,s2,s3);
}


static inline void ff_poly_mgcd_linear_7_5 (ff_t h[2], ff_t f[7], ff_t g[5])					// f monic deg 7, g monic deg 5
{
	register ff_t s0, s1, s2, s3, s4, t0, t1, t2, t3, t4, w1, w2,w3;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]);
	_ff_sub(w2,f[6],s4);
	_ff_sub(w1,f[5],s3); _ff_mult(w3,w2,s4); _ff_subfrom(w1,w3);
	_ff_sub(t4,f[4],s2); _ff_sum_2_mults(w3,w1,w2,s3,s4); _ff_subfrom(t4,w3);
	_ff_sub(t3,f[3],s1); _ff_sum_2_mults(w3,w1,w2,s2,s3); _ff_subfrom(t3,w3);
	_ff_sub(t2,f[2],s0); _ff_sum_2_mults(w3,w1,w2,s1,s2); _ff_subfrom(t2,w3);
	_ff_sum_2_mults(w3,w1,w2,s0,s1); _ff_sub(t1,f[1],w3);
	_ff_mult(w3,w1,s0); _ff_sub(t0,f[0],w3);
	
	// reduce 1 s4 s3 s2 s1 s0 mod t4 t3 t2 t1 t0 here (saves 2 mults versus calling 5_4 code)
	_ff_neg(w1,t4);  _ff_mult(w2,w1,s4); _ff_addto(w2,t3);  _ff_square(w3,w1);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_4_3_reg(h,t0,t1,t2,t3,t4,s0,s1,s2,s3);
	// Total 44M+38A (23 redc)
}

static inline void ff_poly_mgcd_linear_n_5 (ff_t h[2], ff_t *f, int d_f, ff_t g[5])				// f monic deg >= 5, g monic deg 5
{
	register int i;
	register ff_t s0, s1, s2, s3, s4, t0, t1, t2, t3, t4, w1, w2, w3, t, u;
	
	if ( d_f < 9 ) { ff_poly_mgcd_linear_m_5(h,f,d_f,g); return; }
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]);
	_ff_sub(t4,f[d_f-1],s4);
	_ff_sub(t3,f[d_f-2],s3); _ff_mult(t,t4,s4); _ff_subfrom(t3,t);
	_ff_sub(t2,f[d_f-3],s2); _ff_sum_2_mults(t,t3,t4,s3,s4); _ff_subfrom(t2,t);
	_ff_sub(t1,f[d_f-4],s1); _ff_sum_3_mults(t,t2,t3,t4,s2,s3,s4); _ff_subfrom(t1,t);
	_ff_sub(t0,f[d_f-5],s0); _ff_sum_4_mults(t,t1,t2,t3,t4,s1,s2,s3,s4); _ff_subfrom(t0,t);
	for ( i = d_f-6; i > 3 ; i-- ) {
		_ff_sum_5_mults(t,t0,t1,t2,t3,t4,s0,s1,s2,s3,s4);
		_ff_set(t4,t3); _ff_set(t3,t2); _ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[i],t);
	}
	_ff_sum_4_mults(t,t1,t2,t3,t4,s0,s1,s2,s3);  	_ff_set(t4,t0);
	_ff_sum_3_mults(u,t1,t2,t3,s0,s1,s2);  _ff_sub(t3,f[3],t);
	_ff_sum_2_mults(t,t1,t2,s0,s1);  _ff_sub(t2,f[2],u);
	_ff_mult(u,t1,s0);  _ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);

	// reduce 1 s4 s3 s2 s1 s0 mod t4 t3 t2 t1 t0 here (saves 2 mults versus calling 5_4 code)
	_ff_neg(w1,t4); _ff_mult(w2,w1,s4); _ff_addto(w2,t3);  _ff_square(w3,w1);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_4_3_reg(h,t0,t1,t2,t3,t4,s0,s1,s2,s3);
}

static inline void ff_poly_mgcd_linear_m_6 (ff_t h[2], ff_t *f, int d_f, ff_t g[6])				// f monic deg 6 to 10, g monic deg 6
{
	register ff_t s0, s1, s2, s3, s4, s5, t0, t1, t2, t3, t4, t5, w1, w2, w3, w4, t;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]); _ff_set(s5, g[5]);
	switch (d_f) {
	case 6:
		_ff_sub(t5,f[5],s5);  _ff_sub(t4,f[4],s4);  _ff_sub(t3,f[3],s3);  _ff_sub(t2,f[2],s2);  _ff_sub(t1,f[1],s1);  _ff_sub(t0,f[0],s0);  break;
	case 7:
		_ff_sub(w1,f[6],s5);
		_ff_sub(t5,f[5],s4); _ff_mult(t,w1,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s3); _ff_mult(t,w1,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s2); _ff_mult(t,w1,s3); _ff_subfrom(t3,t);
		_ff_sub(t2,f[2],s1); _ff_mult(t,w1,s2); _ff_subfrom(t2,t);
		_ff_sub(t1,f[1],s0); _ff_mult(t,w1,s1); _ff_subfrom(t1,t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 8:
		_ff_sub(w2,f[7],s5);
		_ff_sub(w1,f[6],s4); _ff_mult(t,w2,s5); _ff_subfrom(w1,t);
		_ff_sub(t5,f[5],s3); _ff_sum_2_mults(t,w1,w2,s4,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s2); _ff_sum_2_mults(t,w1,w2,s3,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s1); _ff_sum_2_mults(t,w1,w2,s2,s3); _ff_subfrom(t3,t);
		_ff_sub(t2,f[2],s0); _ff_sum_2_mults(t,w1,w2,s1,s2); _ff_subfrom(t2,t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 9:
		_ff_sub(w3,f[8],s5);
		_ff_sub(w2,f[7],s4); _ff_mult(t,w3,s5); _ff_subfrom(w2,t);
		_ff_sub(w1,f[6],s3); _ff_sum_2_mults(t,w2,w3,s4,s5); _ff_subfrom(w1,t);
		_ff_sub(t5,f[5],s2); _ff_sum_3_mults(t,w1,w2,w3,s3,s4,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s1); _ff_sum_3_mults(t,w1,w2,w3,s2,s3,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s0); _ff_sum_3_mults(t,w1,w2,w3,s1,s2,s3); _ff_subfrom(t3,t);
		_ff_sum_3_mults(t,w1,w2,w3,s0,s1,s2); _ff_sub(t2,f[2],t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 10:
		_ff_sub(w4,f[9],s5);
		_ff_sub(w3,f[8],s4); _ff_mult(t,w4,s5); _ff_subfrom(w3,t);
		_ff_sub(w2,f[7],s3); _ff_sum_2_mults(t,w3,w4,s4,s5); _ff_subfrom(w2,t);
		_ff_sub(w1,f[6],s2); _ff_sum_3_mults(t,w2,w3,w4,s3,s4,s5); _ff_subfrom(w1,t);
		_ff_sub(t5,f[5],s1); _ff_sum_4_mults(t,w1,w2,w3,w4,s2,s3,s4,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s0); _ff_sum_4_mults(t,w1,w2,w3,w4,s1,s2,s3,s4); _ff_subfrom(t4,t);
		_ff_sum_4_mults(t,w1,w2,w3,w4,s0,s1,s2,s3); _ff_sub(t3,f[3],t);
		_ff_sum_3_mults(t,w1,w2,w3,s0,s1,s2); _ff_sub(t2,f[2],t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	default:
		printf ("Invalid degree %d in ff_poly_mgcd_linear_m_6\n", d_f);  abort();
	}

	// reduce 1 s5 s4 s3 s2 s1 s0 mod t5 t4 t3 t2 t1 t0 here (saves 2 mults versus calling 7_6 code)
	_ff_neg(w1,t5);  _ff_mult(w2,w1,s5); _ff_addto(w2,t4);  _ff_square(w3,w1);
	_ff_sum_3_mults(s4,s4,t4,t3,w1,w2,w3);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_5_4_reg(h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
}

static inline void ff_poly_mgcd_linear_n_6 (ff_t h[2], ff_t *f, int d_f, ff_t g[6])				// f monic deg >= 6, g monic deg 6
{
	register int i;
	register ff_t s0, s1, s2, s3, s4, s5, t0, t1, t2, t3, t4, t5, w1, w2, w3, t, u;
	
	if ( d_f < 11 ) { ff_poly_mgcd_linear_m_6(h,f,d_f,g); return; }

	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]); _ff_set(s5,g[5]);
	_ff_sub(t5,f[d_f-1],s5);
	_ff_sub(t4,f[d_f-2],s4); _ff_mult(t,t5,s5); _ff_subfrom(t4,t);
	_ff_sub(t3,f[d_f-3],s3); _ff_sum_2_mults(t,t4,t5,s4,s5); _ff_subfrom(t3,t);
	_ff_sub(t2,f[d_f-4],s2); _ff_sum_3_mults(t,t3,t4,t5,s3,s4,s5); _ff_subfrom(t2,t);
	_ff_sub(t1,f[d_f-5],s1); _ff_sum_4_mults(t,t2,t3,t4,t5,s2,s3,s4,s5); _ff_subfrom(t1,t);
	_ff_sub(t0,f[d_f-6],s0); _ff_sum_5_mults(t,t1,t2,t3,t4,t5,s1,s2,s3,s4,s5); _ff_subfrom(t0,t);
	for ( i = d_f-7 ; i > 4 ; i-- ) {
		_ff_sum_6_mults(t,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4,s5);
		_ff_set(t5,t4); _ff_set(t4,t3); _ff_set(t3,t2); _ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[i],t);
	}
	_ff_sum_5_mults(t,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4); _ff_set(t5,t0);
	_ff_sum_4_mults(u,t1,t2,t3,t4,s0,s1,s2,s3);  _ff_sub(t4,f[4],t);
	_ff_sum_3_mults(t,t1,t2,t3,s0,s1,s2);  	_ff_sub(t3,f[3],u);
	_ff_sum_2_mults(u,t1,t2,s0,s1);  _ff_sub(t2,f[2],t);
	_ff_mult(t,t1,s0);  _ff_sub(t1,f[1],u); _ff_sub(t0,f[0],t);

	// reduce 1 s5 s4 s3 s2 s1 s0 mod t5 t4 t3 t2 t1 t0 here (saves 2 mults versus calling 7_6 code)
	_ff_neg(w1,t5);  _ff_mult(w2,w1,s5); _ff_addto(w2,t4);  _ff_square(w3,w1);
	_ff_sum_3_mults(s4,s4,t4,t3,w1,w2,w3);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_5_4_reg(h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
}

static inline void ff_poly_mgcd_linear_12_7 (ff_t h[2], ff_t f[12], ff_t g[7])				// f monic deg 12, g monic deg 7
{
	register ff_t s0, s1, s2, s3, s4, s5, s6, t0, t1, t2, t3, t4, t5, t6, w1, w2, w3, w4, w5, t;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]); _ff_set(s5, g[5]); _ff_set(s6,g[6]);
	_ff_sub(w5,f[11],s6);
	_ff_sub(w4,f[10],s5); _ff_mult(t,w5,s6); _ff_subfrom(w4,t);
	_ff_sub(w3,f[9],s4); _ff_sum_2_mults(t,w4,w5,s5,s6); _ff_subfrom(w3,t);
	_ff_sub(w2,f[8],s3); _ff_sum_3_mults(t,w3,w4,w5,s4,s5,s6); _ff_subfrom(w2,t);
	_ff_sub(w1,f[7],s2); _ff_sum_4_mults(t,w2,w3,w4,w5,s3,s4,s5,s6); _ff_subfrom(w1,t);
	_ff_sub(t6,f[6],s1); _ff_sum_5_mults(t,w1,w2,w3,w4,w5,s2,s3,s4,s5,s6); _ff_subfrom(t6,t);
	_ff_sub(t5,f[5],s0); _ff_sum_5_mults(t,w1,w2,w3,w4,w5,s1,s2,s3,s4,s5); _ff_subfrom(t5,t);
	 _ff_sum_5_mults(t,w1,w2,w3,w4,w5,s0,s1,s2,s3,s4); _ff_sub(t4,f[4],t);
	_ff_sum_4_mults(t,w1,w2,w3,w4,s0,s1,s2,s3); _ff_sub(t3,f[3],t);
	_ff_sum_3_mults(t,w1,w2,w3,s0,s1,s2); _ff_sub(t2,f[2],t);
	_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
	_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);

	// reduce 1 s6 s5 s4 s3 s2 s1 s0 mod t6 t5 t4 t3 t2 t1 t0 here (saves 2 mults versus calling 7_6 code)
	_ff_neg(w1,t6);  _ff_mult(w2,w1,s6); _ff_addto(w2,t5);  _ff_square(w3,w1);
	_ff_sum_3_mults(s5,s5,t5,t4,w1,w2,w3);
	_ff_sum_3_mults(s4,s4,t4,t3,w1,w2,w3);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_6_5_reg(h,t0,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5);
}


static inline void ff_poly_mgcd_linear_m_7 (ff_t h[2], ff_t *f, int d_f, ff_t g[7])				// f monic deg 7 to 12, g monic deg 7
{
	register ff_t s0, s1, s2, s3, s4, s5, s6, t0, t1, t2, t3, t4, t5, t6, w1, w2, w3, w4, w5, t;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]); _ff_set(s5, g[5]); _ff_set(s6,g[6]);
	switch (d_f) {
	case 7:
		_ff_sub(t6,f[6],s6);  _ff_sub(t5,f[5],s5);  _ff_sub(t4,f[4],s4);  _ff_sub(t3,f[3],s3);  _ff_sub(t2,f[2],s2);  _ff_sub(t1,f[1],s1);  _ff_sub(t0,f[0],s0);  break;
	case 8:
		_ff_sub(w1,f[7],s6);
		_ff_sub(t6,f[6],s5); _ff_mult(t,w1,s6); _ff_subfrom(t6,t);
		_ff_sub(t5,f[5],s4); _ff_mult(t,w1,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s3); _ff_mult(t,w1,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s2); _ff_mult(t,w1,s3); _ff_subfrom(t3,t);
		_ff_sub(t2,f[2],s1); _ff_mult(t,w1,s2); _ff_subfrom(t2,t);
		_ff_sub(t1,f[1],s0); _ff_mult(t,w1,s1); _ff_subfrom(t1,t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 9:
		_ff_sub(w2,f[8],s6);
		_ff_sub(w1,f[7],s5); _ff_mult(t,w2,s6); _ff_subfrom(w1,t);
		_ff_sub(t6,f[6],s4); _ff_sum_2_mults(t,w1,w2,s5,s6); _ff_subfrom(t6,t);
		_ff_sub(t5,f[5],s3); _ff_sum_2_mults(t,w1,w2,s4,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s2); _ff_sum_2_mults(t,w1,w2,s3,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s1); _ff_sum_2_mults(t,w1,w2,s2,s3); _ff_subfrom(t3,t);
		_ff_sub(t2,f[2],s0); _ff_sum_2_mults(t,w1,w2,s1,s2); _ff_subfrom(t2,t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 10:
		_ff_sub(w3,f[9],s6);
		_ff_sub(w2,f[8],s5); _ff_mult(t,w3,s6); _ff_subfrom(w2,t);
		_ff_sub(w1,f[7],s4); _ff_sum_2_mults(t,w2,w3,s5,s6); _ff_subfrom(w1,t);
		_ff_sub(t6,f[6],s3); _ff_sum_3_mults(t,w1,w2,w3,s4,s5,s6); _ff_subfrom(t6,t);
		_ff_sub(t5,f[5],s2); _ff_sum_3_mults(t,w1,w2,w3,s3,s4,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s1); _ff_sum_3_mults(t,w1,w2,w3,s2,s3,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s0); _ff_sum_3_mults(t,w1,w2,w3,s1,s2,s3); _ff_subfrom(t3,t);
		_ff_sum_3_mults(t,w1,w2,w3,s0,s1,s2); _ff_sub(t2,f[2],t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 11:
		_ff_sub(w4,f[10],s6);
		_ff_sub(w3,f[9],s5); _ff_mult(t,w4,s6); _ff_subfrom(w3,t);
		_ff_sub(w2,f[8],s4); _ff_sum_2_mults(t,w3,w4,s5,s6); _ff_subfrom(w2,t);
		_ff_sub(w1,f[7],s3); _ff_sum_3_mults(t,w2,w3,w4,s4,s5,s6); _ff_subfrom(w1,t);
		_ff_sub(t6,f[6],s2); _ff_sum_4_mults(t,w1,w2,w3,w4,s3,s4,s5,s6); _ff_subfrom(t6,t);
		_ff_sub(t5,f[5],s1); _ff_sum_4_mults(t,w1,w2,w3,w4,s2,s3,s4,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s0); _ff_sum_4_mults(t,w1,w2,w3,w4,s1,s2,s3,s4); _ff_subfrom(t4,t);
		_ff_sum_4_mults(t,w1,w2,w3,w4,s0,s1,s2,s3); _ff_sub(t3,f[3],t);
		_ff_sum_3_mults(t,w1,w2,w3,s0,s1,s2); _ff_sub(t2,f[2],t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 12:
		_ff_sub(w5,f[11],s6);
		_ff_sub(w4,f[10],s5); _ff_mult(t,w5,s6); _ff_subfrom(w4,t);
		_ff_sub(w3,f[9],s4); _ff_sum_2_mults(t,w4,w5,s5,s6); _ff_subfrom(w3,t);
		_ff_sub(w2,f[8],s3); _ff_sum_3_mults(t,w3,w4,w5,s4,s5,s6); _ff_subfrom(w2,t);
		_ff_sub(w1,f[7],s2); _ff_sum_4_mults(t,w2,w3,w4,w5,s3,s4,s5,s6); _ff_subfrom(w1,t);
		_ff_sub(t6,f[6],s1); _ff_sum_5_mults(t,w1,w2,w3,w4,w5,s2,s3,s4,s5,s6); _ff_subfrom(t6,t);
		_ff_sub(t5,f[5],s0); _ff_sum_5_mults(t,w1,w2,w3,w4,w5,s1,s2,s3,s4,s5); _ff_subfrom(t5,t);
		 _ff_sum_5_mults(t,w1,w2,w3,w4,w5,s0,s1,s2,s3,s4); _ff_sub(t4,f[4],t);
		_ff_sum_4_mults(t,w1,w2,w3,w4,s0,s1,s2,s3); _ff_sub(t3,f[3],t);
		_ff_sum_3_mults(t,w1,w2,w3,s0,s1,s2); _ff_sub(t2,f[2],t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	default:
		printf ("Invalid degree %d in ff_poly_mgcd_linear_m_7\n", d_f); abort();
	}

	// reduce 1 s6 s5 s4 s3 s2 s1 s0 mod t6 t5 t4 t3 t2 t1 t0 here (saves 2 mults versus calling 7_6 code)
	_ff_neg(w1,t6);  _ff_mult(w2,w1,s6); _ff_addto(w2,t5);  _ff_square(w3,w1);
	_ff_sum_3_mults(s5,s5,t5,t4,w1,w2,w3);
	_ff_sum_3_mults(s4,s4,t4,t3,w1,w2,w3);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_6_5_reg(h,t0,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5);
}

static inline void ff_poly_mgcd_linear_14_7 (ff_t h[2], ff_t f[14],  ff_t g[7])				// f monic deg 14, g monic deg 7
{
	register ff_t s0, s1, s2, s3, s4, s5, s6, t0, t1, t2, t3, t4, t5, t6, w1, w2, w3, t, u;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]); _ff_set(s5,g[5]); _ff_set(s6,g[6]);
	_ff_sub(t6,f[13],s6);
	_ff_sub(t5,f[12],s5); _ff_mult(t,t6,s6); _ff_subfrom(t5,t);
	_ff_sub(t4,f[11],s4); _ff_sum_2_mults(t,t5,t6,s5,s6); _ff_subfrom(t4,t);
	_ff_sub(t3,f[10],s3); _ff_sum_3_mults(t,t4,t5,t6,s4,s5,s6); _ff_subfrom(t3,t);
	_ff_sub(t2,f[9],s2); _ff_sum_4_mults(t,t3,t4,t5,t6,s3,s4,s5,s6); _ff_subfrom(t2,t);
	_ff_sub(t1,f[8],s1); _ff_sum_5_mults(t,t2,t3,t4,t5,t6,s2,s3,s4,s5,s6); _ff_subfrom(t1,t);
	_ff_sub(t0,f[7],s0); _ff_sum_6_mults(t,t1,t2,t3,t4,t5,t6,s1,s2,s3,s4,s5,s6); _ff_subfrom(t0,t);
	_ff_sum_7_mults(t,t0,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5,s6);
	_ff_set(t6,t5); _ff_set(t5,t4); _ff_set(t4,t3); _ff_set(t3,t2); _ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[6],t);		// this shift could be avoided
	_ff_sum_6_mults(t,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5); _ff_set(t6,t0);
	_ff_sum_5_mults(u,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);  _ff_sub(t5,f[5],t);
	_ff_sum_4_mults(t,t1,t2,t3,t4,s0,s1,s2,s3);  	_ff_sub(t4,f[4],u);
	_ff_sum_3_mults(u,t1,t2,t3,s0,s1,s2);  _ff_sub(t3,f[3],t);
	_ff_sum_2_mults(t,t1,t2,s0,s1);  _ff_sub(t2,f[2],u);
	_ff_mult(u,t1,s0);  _ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);

	// reduce 1 s6 s5 s4 s3 s2 s1 s0 mod t6 t5 t4 t3 t2 t1 t0 here (saves 2 mults versus calling 7_6 code)
	_ff_neg(w1,t6);  _ff_mult(w2,w1,s6); _ff_addto(w2,t5);  _ff_square(w3,w1);
	_ff_sum_3_mults(s5,s5,t5,t4,w1,w2,w3);
	_ff_sum_3_mults(s4,s4,t4,t3,w1,w2,w3);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_6_5_reg(h,t0,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5);
}

static inline void ff_poly_mgcd_linear_18_7 (ff_t h[2], ff_t f[18], ff_t g[7])				// f monic deg 18, g monic deg 7
{
	register ff_t s0, s1, s2, s3, s4, s5, s6, t0, t1, t2, t3, t4, t5, t6, w1, w2, w3, t, u;

	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]); _ff_set(s5,g[5]); _ff_set(s6,g[6]);
	_ff_sub(t6,f[17],s6);
	_ff_sub(t5,f[16],s5); _ff_mult(t,t6,s6); _ff_subfrom(t5,t);
	_ff_sub(t4,f[15],s4); _ff_sum_2_mults(t,t5,t6,s5,s6); _ff_subfrom(t4,t);
	_ff_sub(t3,f[14],s3); _ff_sum_3_mults(t,t4,t5,t6,s4,s5,s6); _ff_subfrom(t3,t);
	_ff_sub(t2,f[13],s2); _ff_sum_4_mults(t,t3,t4,t5,t6,s3,s4,s5,s6); _ff_subfrom(t2,t);
	_ff_sub(t1,f[12],s1); _ff_sum_5_mults(t,t2,t3,t4,t5,t6,s2,s3,s4,s5,s6); _ff_subfrom(t1,t);
	_ff_sub(t0,f[11],s0); _ff_sum_6_mults(t,t1,t2,t3,t4,t5,t6,s1,s2,s3,s4,s5,s6); _ff_subfrom(t0,t);
	_ff_sum_7_mults(t,t0,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5,s6); _ff_sub(t6,f[10],t);
	_ff_sum_7_mults(t,t6,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4,s5,s6); _ff_sub(t5,f[9],t);
	_ff_sum_7_mults(t,t5,t6,t0,t1,t2,t3,t4,s0,s1,s2,s3,s4,s5,s6); _ff_sub(t4,f[8],t);
	_ff_sum_7_mults(t,t4,t5,t6,t0,t1,t2,t3,s0,s1,s2,s3,s4,s5,s6); _ff_sub(t3,f[7],t);
	_ff_sum_7_mults(t,t3,t4,t5,t6,t0,t1,t2,s0,s1,s2,s3,s4,s5,s6);
	_ff_set(t2,t4); _ff_set(t4,t6); _ff_set(t6,t1); _ff_set(t1,t3); _ff_set(t3,t5); _ff_set(t5,t0); _ff_sub(t0,f[6],t);	// this shift could be avoided
	_ff_sum_6_mults(t,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5); _ff_set(t6,t0);
	_ff_sum_5_mults(u,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);  _ff_sub(t5,f[5],t);
	_ff_sum_4_mults(t,t1,t2,t3,t4,s0,s1,s2,s3);  	_ff_sub(t4,f[4],u);
	_ff_sum_3_mults(u,t1,t2,t3,s0,s1,s2);  _ff_sub(t3,f[3],t);
	_ff_sum_2_mults(t,t1,t2,s0,s1);  _ff_sub(t2,f[2],u);
	_ff_mult(u,t1,s0);  _ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);

	// reduce 1 s6 s5 s4 s3 s2 s1 s0 mod t6 t5 t4 t3 t2 t1 t0 here (saves 2 mults versus calling 7_6 code)
	_ff_neg(w1,t6);  _ff_mult(w2,w1,s6); _ff_addto(w2,t5);  _ff_square(w3,w1);
	_ff_sum_3_mults(s5,s5,t5,t4,w1,w2,w3);
	_ff_sum_3_mults(s4,s4,t4,t3,w1,w2,w3);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_6_5_reg(h,t0,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5);
}

static inline void ff_poly_mgcd_linear_20_7 (ff_t h[2], ff_t f[20], ff_t g[7])				// f monic deg 20, g monic deg 7
{
	register ff_t s0, s1, s2, s3, s4, s5, s6, t0, t1, t2, t3, t4, t5, t6, w1, w2, w3, t, u;

	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]); _ff_set(s5,g[5]); _ff_set(s6,g[6]);
	_ff_sub(t6,f[19],s6);
	_ff_sub(t5,f[18],s5); _ff_mult(t,t6,s6); _ff_subfrom(t5,t);
	_ff_sub(t4,f[17],s4); _ff_sum_2_mults(t,t5,t6,s5,s6); _ff_subfrom(t4,t);
	_ff_sub(t3,f[16],s3); _ff_sum_3_mults(t,t4,t5,t6,s4,s5,s6); _ff_subfrom(t3,t);
	_ff_sub(t2,f[15],s2); _ff_sum_4_mults(t,t3,t4,t5,t6,s3,s4,s5,s6); _ff_subfrom(t2,t);
	_ff_sub(t1,f[14],s1); _ff_sum_5_mults(t,t2,t3,t4,t5,t6,s2,s3,s4,s5,s6); _ff_subfrom(t1,t);
	_ff_sub(t0,f[13],s0); _ff_sum_6_mults(t,t1,t2,t3,t4,t5,t6,s1,s2,s3,s4,s5,s6); _ff_subfrom(t0,t);
	_ff_sum_7_mults(t,t0,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5,s6); _ff_sub(t6,f[12],t);
	_ff_sum_7_mults(t,t6,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4,s5,s6); _ff_sub(t5,f[11],t);
	_ff_sum_7_mults(t,t5,t6,t0,t1,t2,t3,t4,s0,s1,s2,s3,s4,s5,s6); _ff_sub(t4,f[10],t);
	_ff_sum_7_mults(t,t4,t5,t6,t0,t1,t2,t3,s0,s1,s2,s3,s4,s5,s6); _ff_sub(t3,f[9],t);
	_ff_sum_7_mults(t,t3,t4,t5,t6,t0,t1,t2,s0,s1,s2,s3,s4,s5,s6); _ff_sub(t2,f[8],t);
	_ff_sum_7_mults(t,t2,t3,t4,t5,t6,t0,t1,s0,s1,s2,s3,s4,s5,s6); _ff_sub(t1,f[7],t);
	_ff_sum_7_mults(t,t1,t2,t3,t4,t5,t6,t0,s0,s1,s2,s3,s4,s5,s6); _ff_sub(t0,f[6],t);
	_ff_sum_6_mults(t,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5); _ff_set(t6,t0);
	_ff_sum_5_mults(u,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);  _ff_sub(t5,f[5],t);
	_ff_sum_4_mults(t,t1,t2,t3,t4,s0,s1,s2,s3);  	_ff_sub(t4,f[4],u);
	_ff_sum_3_mults(u,t1,t2,t3,s0,s1,s2);  _ff_sub(t3,f[3],t);
	_ff_sum_2_mults(t,t1,t2,s0,s1);  _ff_sub(t2,f[2],u);
	_ff_mult(u,t1,s0);  _ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);

	// reduce 1 s6 s5 s4 s3 s2 s1 s0 mod t6 t5 t4 t3 t2 t1 t0 here (saves 2 mults versus calling 7_6 code)
	_ff_neg(w1,t6);  _ff_mult(w2,w1,s6); _ff_addto(w2,t5);  _ff_square(w3,w1);
	_ff_sum_3_mults(s5,s5,t5,t4,w1,w2,w3);
	_ff_sum_3_mults(s4,s4,t4,t3,w1,w2,w3);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_6_5_reg(h,t0,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5);
}

static inline void ff_poly_mgcd_linear_n_7 (ff_t h[2], ff_t *f, int d_f, ff_t g[7])				// f monic deg >= 7, g monic deg 7
{
	register int i;
	register ff_t s0, s1, s2, s3, s4, s5, s6, t0, t1, t2, t3, t4, t5, t6, w1, w2, w3, t, u;
	
	if ( d_f < 13 ) { ff_poly_mgcd_linear_m_7(h,f,d_f,g); return; }

	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]); _ff_set(s5,g[5]); _ff_set(s6,g[6]);
	_ff_sub(t6,f[d_f-1],s6);
	_ff_sub(t5,f[d_f-2],s5); _ff_mult(t,t6,s6); _ff_subfrom(t5,t);
	_ff_sub(t4,f[d_f-3],s4); _ff_sum_2_mults(t,t5,t6,s5,s6); _ff_subfrom(t4,t);
	_ff_sub(t3,f[d_f-4],s3); _ff_sum_3_mults(t,t4,t5,t6,s4,s5,s6); _ff_subfrom(t3,t);
	_ff_sub(t2,f[d_f-5],s2); _ff_sum_4_mults(t,t3,t4,t5,t6,s3,s4,s5,s6); _ff_subfrom(t2,t);
	_ff_sub(t1,f[d_f-6],s1); _ff_sum_5_mults(t,t2,t3,t4,t5,t6,s2,s3,s4,s5,s6); _ff_subfrom(t1,t);
	_ff_sub(t0,f[d_f-7],s0); _ff_sum_6_mults(t,t1,t2,t3,t4,t5,t6,s1,s2,s3,s4,s5,s6); _ff_subfrom(t0,t);
	for ( i = d_f-8 ; i > 5 ; i-- ) {
		_ff_sum_7_mults(t,t0,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5,s6);
		_ff_set(t6,t5); _ff_set(t5,t4); _ff_set(t4,t3); _ff_set(t3,t2); _ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[i],t);
	}
	_ff_sum_6_mults(t,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5); _ff_set(t6,t0);
	_ff_sum_5_mults(u,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);  _ff_sub(t5,f[5],t);
	_ff_sum_4_mults(t,t1,t2,t3,t4,s0,s1,s2,s3);  	_ff_sub(t4,f[4],u);
	_ff_sum_3_mults(u,t1,t2,t3,s0,s1,s2);  _ff_sub(t3,f[3],t);
	_ff_sum_2_mults(t,t1,t2,s0,s1);  _ff_sub(t2,f[2],u);
	_ff_mult(u,t1,s0);  _ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);

	// reduce 1 s6 s5 s4 s3 s2 s1 s0 mod t6 t5 t4 t3 t2 t1 t0 here (saves 2 mults versus calling 7_6 code)
	_ff_neg(w1,t6);  _ff_mult(w2,w1,s6); _ff_addto(w2,t5);  _ff_square(w3,w1);
	_ff_sum_3_mults(s5,s5,t5,t4,w1,w2,w3);
	_ff_sum_3_mults(s4,s4,t4,t3,w1,w2,w3);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_6_5_reg(h,t0,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5);
}

static inline void ff_poly_mgcd_linear_m_8 (ff_t h[2], ff_t *f, int d_f, ff_t g[8])				// f monic deg 8 to 14, g monic deg 8
{
	register ff_t s0, s1, s2, s3, s4, s5, s6, s7, t0, t1, t2, t3, t4, t5, t6, t7, w1, w2, w3, w4, w5, w6, t;
	
	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]); _ff_set(s5, g[5]); _ff_set(s6,g[6]); _ff_set(s7,g[7]);
	switch (d_f) {
	case 8:
		_ff_sub(t7,f[7],s7); _ff_sub(t6,f[6],s6);  _ff_sub(t5,f[5],s5);  _ff_sub(t4,f[4],s4);  _ff_sub(t3,f[3],s3);  _ff_sub(t2,f[2],s2);  _ff_sub(t1,f[1],s1);  _ff_sub(t0,f[0],s0);  break;
	case 9:
		_ff_sub(w1,f[8],s7);
		_ff_sub(t7,f[7],s6); _ff_mult(t,w1,s7); _ff_subfrom(t7,t);
		_ff_sub(t6,f[6],s5); _ff_mult(t,w1,s6); _ff_subfrom(t6,t);
		_ff_sub(t5,f[5],s4); _ff_mult(t,w1,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s3); _ff_mult(t,w1,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s2); _ff_mult(t,w1,s3); _ff_subfrom(t3,t);
		_ff_sub(t2,f[2],s1); _ff_mult(t,w1,s2); _ff_subfrom(t2,t);
		_ff_sub(t1,f[1],s0); _ff_mult(t,w1,s1); _ff_subfrom(t1,t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 10:
		_ff_sub(w2,f[9],s7);
		_ff_sub(w1,f[8],s6); _ff_mult(t,w2,s7); _ff_subfrom(w1,t);
		_ff_sub(t7,f[7],s5); _ff_sum_2_mults(t,w1,w2,s6,s7); _ff_subfrom(t7,t);
		_ff_sub(t6,f[6],s4); _ff_sum_2_mults(t,w1,w2,s5,s6); _ff_subfrom(t6,t);
		_ff_sub(t5,f[5],s3); _ff_sum_2_mults(t,w1,w2,s4,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s2); _ff_sum_2_mults(t,w1,w2,s3,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s1); _ff_sum_2_mults(t,w1,w2,s2,s3); _ff_subfrom(t3,t);
		_ff_sub(t2,f[2],s0); _ff_sum_2_mults(t,w1,w2,s1,s2); _ff_subfrom(t2,t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 11:
		_ff_sub(w3,f[10],s7);
		_ff_sub(w2,f[9],s6); _ff_mult(t,w3,s7); _ff_subfrom(w2,t);
		_ff_sub(w1,f[8],s5); _ff_sum_2_mults(t,w2,w3,s6,s7); _ff_subfrom(w1,t);
		_ff_sub(t7,f[7],s4); _ff_sum_3_mults(t,w1,w2,w3,s5,s6,s7); _ff_subfrom(t7,t);
		_ff_sub(t6,f[6],s3); _ff_sum_3_mults(t,w1,w2,w3,s4,s5,s6); _ff_subfrom(t6,t);
		_ff_sub(t5,f[5],s2); _ff_sum_3_mults(t,w1,w2,w3,s3,s4,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s1); _ff_sum_3_mults(t,w1,w2,w3,s2,s3,s4); _ff_subfrom(t4,t);
		_ff_sub(t3,f[3],s0); _ff_sum_3_mults(t,w1,w2,w3,s1,s2,s3); _ff_subfrom(t3,t);
		_ff_sum_3_mults(t,w1,w2,w3,s0,s1,s2); _ff_sub(t2,f[2],t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 12:
		_ff_sub(w4,f[11],s7);
		_ff_sub(w3,f[10],s6); _ff_mult(t,w4,s7); _ff_subfrom(w3,t);
		_ff_sub(w2,f[9],s5); _ff_sum_2_mults(t,w3,w4,s6,s7); _ff_subfrom(w2,t);
		_ff_sub(w1,f[8],s4); _ff_sum_3_mults(t,w2,w3,w4,s5,s6,s7); _ff_subfrom(w1,t);
		_ff_sub(t7,f[7],s3); _ff_sum_4_mults(t,w1,w2,w3,w4,s4,s5,s6,s7); _ff_subfrom(t7,t);
		_ff_sub(t6,f[6],s2); _ff_sum_4_mults(t,w1,w2,w3,w4,s3,s4,s5,s6); _ff_subfrom(t6,t);
		_ff_sub(t5,f[5],s1); _ff_sum_4_mults(t,w1,w2,w3,w4,s2,s3,s4,s5); _ff_subfrom(t5,t);
		_ff_sub(t4,f[4],s0); _ff_sum_4_mults(t,w1,w2,w3,w4,s1,s2,s3,s4); _ff_subfrom(t4,t);
		_ff_sum_4_mults(t,w1,w2,w3,w4,s0,s1,s2,s3); _ff_sub(t3,f[3],t);
		_ff_sum_3_mults(t,w1,w2,w3,s0,s1,s2); _ff_sub(t2,f[2],t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 13:
		_ff_sub(w5,f[12],s7);
		_ff_sub(w4,f[11],s6); _ff_mult(t,w5,s7); _ff_subfrom(w4,t);
		_ff_sub(w3,f[10],s5); _ff_sum_2_mults(t,w4,w5,s6,s7); _ff_subfrom(w3,t);
		_ff_sub(w2,f[9],s4); _ff_sum_3_mults(t,w3,w4,w5,s5,s6,s7); _ff_subfrom(w2,t);
		_ff_sub(w1,f[8],s3); _ff_sum_4_mults(t,w2,w3,w4,w5,s4,s5,s6,s7); _ff_subfrom(w1,t);
		_ff_sub(t7,f[7],s2); _ff_sum_5_mults(t,w1,w2,w3,w4,w5,s3,s4,s5,s6,s7); _ff_subfrom(t7,t);
		_ff_sub(t6,f[6],s1); _ff_sum_5_mults(t,w1,w2,w3,w4,w5,s2,s3,s4,s5,s6); _ff_subfrom(t6,t);
		_ff_sub(t5,f[5],s0); _ff_sum_5_mults(t,w1,w2,w3,w4,w5,s1,s2,s3,s4,s5); _ff_subfrom(t5,t);
		 _ff_sum_5_mults(t,w1,w2,w3,w4,w5,s0,s1,s2,s3,s4); _ff_sub(t4,f[4],t);
		_ff_sum_4_mults(t,w1,w2,w3,w4,s0,s1,s2,s3); _ff_sub(t3,f[3],t);
		_ff_sum_3_mults(t,w1,w2,w3,s0,s1,s2); _ff_sub(t2,f[2],t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	case 14:
		_ff_sub(w6,f[13],s7);
		_ff_sub(w5,f[12],s6); _ff_mult(t,w6,s7); _ff_subfrom(w5,t);
		_ff_sub(w4,f[11],s5); _ff_sum_2_mults(t,w5,w6,s6,s7); _ff_subfrom(w4,t);
		_ff_sub(w3,f[10],s4); _ff_sum_3_mults(t,w4,w5,w6,s5,s6,s7); _ff_subfrom(w3,t);
		_ff_sub(w2,f[9],s3); _ff_sum_4_mults(t,w3,w4,w5,w6,s4,s5,s6,s7); _ff_subfrom(w2,t);
		_ff_sub(w1,f[8],s2); _ff_sum_5_mults(t,w2,w3,w4,w5,w6,s3,s4,s5,s6,s7); _ff_subfrom(w1,t);
		_ff_sub(t7,f[7],s1); _ff_sum_6_mults(t,w1,w2,w3,w4,w5,w6,s2,s3,s4,s5,s6,s7); _ff_subfrom(t7,t);
		_ff_sub(t6,f[6],s0); _ff_sum_6_mults(t,w1,w2,w3,w4,w5,w6,s1,s2,s3,s4,s5,s6); _ff_subfrom(t6,t);
		 _ff_sum_6_mults(t,w1,w2,w3,w4,w5,w6,s0,s1,s2,s3,s4,s5); _ff_sub(t5,f[5],t);
		 _ff_sum_5_mults(t,w1,w2,w3,w4,w5,s0,s1,s2,s3,s4); _ff_sub(t4,f[4],t);
		_ff_sum_4_mults(t,w1,w2,w3,w4,s0,s1,s2,s3); _ff_sub(t3,f[3],t);
		_ff_sum_3_mults(t,w1,w2,w3,s0,s1,s2); _ff_sub(t2,f[2],t);
		_ff_sum_2_mults(t,w1,w2,s0,s1); _ff_sub(t1,f[1],t);
		_ff_mult(t,w1,s0); _ff_sub(t0,f[0],t);
		break;
	default:
		printf ("Bad degree %d in ff_poly_mgcd_linear_m_8\n", d_f); abort();
	}

	// reduce 1 s7 s6 s5 s4 s3 s2 s1 s0 mod t7 t6 t5 t4 t3 t2 t1 t0 here (saves 2 mults versus calling 7_6 code)
	_ff_neg(w1,t7);  _ff_mult(w2,w1,s7); _ff_addto(w2,t6);  _ff_square(w3,w1);
	_ff_sum_3_mults(s6,s6,t6,t5,w1,w2,w3);
	_ff_sum_3_mults(s5,s5,t5,t4,w1,w2,w3);
	_ff_sum_3_mults(s4,s4,t4,t3,w1,w2,w3);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_7_6_reg(h,t0,t1,t2,t3,t4,t5,t6,t7,s0,s1,s2,s3,s4,s5,s6);
}

static inline void ff_poly_mgcd_linear_n_8 (ff_t h[2], ff_t *f, int d_f, ff_t g[8])				// f monic deg >= 7, g monic deg 7
{
	register int i;
	register ff_t s0, s1, s2, s3, s4, s5, s6, s7, t0, t1, t2, t3, t4, t5, t6, t7, w1, w2, w3, t, u;
	
	if ( d_f < 15 ) { ff_poly_mgcd_linear_m_8(h,f,d_f,g); return; }

	_ff_set(s0,g[0]); _ff_set(s1,g[1]); _ff_set(s2,g[2]); _ff_set(s3, g[3]); _ff_set(s4,g[4]); _ff_set(s5,g[5]); _ff_set(s6,g[6]); _ff_set(s7,g[7]);
	_ff_sub(t7,f[d_f-1],s7);
	_ff_sub(t6,f[d_f-2],s6); _ff_mult(t,t7,s7); _ff_subfrom(t6,t);
	_ff_sub(t5,f[d_f-3],s5); _ff_sum_2_mults(t,t6,t7,s6,s7); _ff_subfrom(t5,t);
	_ff_sub(t4,f[d_f-4],s4); _ff_sum_3_mults(t,t5,t6,t7,s5,s6,s7); _ff_subfrom(t4,t);
	_ff_sub(t3,f[d_f-5],s3); _ff_sum_4_mults(t,t4,t5,t6,t7,s4,s5,s6,s7); _ff_subfrom(t3,t);
	_ff_sub(t2,f[d_f-6],s2); _ff_sum_5_mults(t,t3,t4,t5,t6,t7,s3,s4,s5,s6,s7); _ff_subfrom(t2,t);
	_ff_sub(t1,f[d_f-7],s1); _ff_sum_6_mults(t,t2,t3,t4,t5,t6,t7,s2,s3,s4,s5,s6,s7); _ff_subfrom(t1,t);
	_ff_sub(t0,f[d_f-8],s0); _ff_sum_7_mults(t,t1,t2,t3,t4,t5,t6,t7,s1,s2,s3,s4,s5,s6,s7); _ff_subfrom(t0,t);
	for ( i = d_f-9 ; i > 6 ; i-- ) {
		_ff_sum_8_mults(t,t0,t1,t2,t3,t4,t5,t6,t7,s0,s1,s2,s3,s4,s5,s6,s7);
		_ff_set(t7,t6); _ff_set(t6,t5); _ff_set(t5,t4); _ff_set(t4,t3); _ff_set(t3,t2); _ff_set(t2,t1); _ff_set(t1,t0); _ff_sub(t0,f[i],t);
	}
	_ff_sum_7_mults(u,t1,t2,t3,t4,t5,t6,t7,s0,s1,s2,s3,s4,s5,s6); _ff_set(t7,t0);
	_ff_sum_6_mults(t,t1,t2,t3,t4,t5,t6,s0,s1,s2,s3,s4,s5); _ff_sub(t6,f[6],u);
	_ff_sum_5_mults(u,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);  _ff_sub(t5,f[5],t);
	_ff_sum_4_mults(t,t1,t2,t3,t4,s0,s1,s2,s3);  	_ff_sub(t4,f[4],u);
	_ff_sum_3_mults(u,t1,t2,t3,s0,s1,s2);  _ff_sub(t3,f[3],t);
	_ff_sum_2_mults(t,t1,t2,s0,s1);  _ff_sub(t2,f[2],u);
	_ff_mult(u,t1,s0);  _ff_sub(t1,f[1],t); _ff_sub(t0,f[0],u);

	// reduce 1 s7 s6 s5 s4 s3 s2 s1 s0 mod t7 t6 t5 t4 t3 t2 t1 t0 here (saves 2 mults versus calling 7_6 code)
	_ff_neg(w1,t7);  _ff_mult(w2,w1,s7); _ff_addto(w2,t6);  _ff_square(w3,w1);
	_ff_sum_3_mults(s6,s6,t6,t5,w1,w2,w3);
	_ff_sum_3_mults(s5,s5,t5,t4,w1,w2,w3);
	_ff_sum_3_mults(s4,s4,t4,t3,w1,w2,w3);
	_ff_sum_3_mults(s3,s3,t3,t2,w1,w2,w3);
	_ff_sum_3_mults(s2,s2,t2,t1,w1,w2,w3);
	_ff_sum_3_mults(s1,s1,t1,t0,w1,w2,w3);
	_ff_sum_2_mults(s0,s0,t0,w2,w3);

	ff_poly_gcd_linear_7_6_reg(h,t0,t1,t2,t3,t4,t5,t6,t7,s0,s1,s2,s3,s4,s5,s6);
}

// Replaces f by f mod g, with g monic (not verified)
static inline void ff_poly_mod_small_old (ff_t *f, int d_f, ff_t *g, int d_g)
{
	register int j, k, m;
	register ff_t t0;
	
//printf ("Computing "); poly_print(f,d_f); printf (" mod "); poly_print(g,d_g); puts("");
	
	for ( m = d_g-1 ; m >= 0 && _ff_zero(g[m]) ; m-- );
	for ( j = d_f ; j >= d_g ; j-- ) if ( ! _ff_zero(f[j]) ) for ( k = 0 ; k <= m ; k++ ) { _ff_mult(t0,f[j],g[k]); _ff_subfrom(f[j-d_g+k],t0); }

//printf ("Result "); poly_print(f,d_g-1); puts("");
}

// The poly_mod functions assume g is monic, but f need not be

// assumes d_g >= 2, d_f = d_g+1
static inline void ff_poly_mod_delta1 (ff_t *f, ff_t *g, int d_g)
{
	register int i;
	register ff_t t0;
	
	_ff_mult(t0,f[d_g+1],g[d_g-1]); _ff_subfrom(f[d_g],t0);
	for ( i = d_g-1 ; i ; i-- ) { _ff_sum_2_mults_arr(t0,f+d_g,g+i-1); _ff_subfrom(f[i],t0); }
	_ff_mult(t0,f[d_g],g[0]); _ff_subfrom(f[0],t0);
}

// assumes d_g >= 3, d_f = d_g+2
static inline void ff_poly_mod_delta2 (ff_t *f, ff_t *g, int d_g)
{
	register int i;
	register ff_t t0;
	
	_ff_mult(t0,f[d_g+2],g[d_g-1]); _ff_subfrom(f[d_g+1],t0);
	_ff_sum_2_mults_arr(t0,f+d_g+1,g+d_g-2); _ff_subfrom(f[d_g],t0);
	for ( i = d_g-1 ; i>1 ; i-- ) { _ff_sum_3_mults_arr(t0,f+d_g,g+i-2); _ff_subfrom(f[i],t0); }
	_ff_sum_2_mults_arr(t0,f+d_g,g); _ff_subfrom(f[1],t0);
	_ff_mult(t0,f[d_g],g[0]); _ff_subfrom(f[0],t0);
}

// assumes d_g >= 4, d_f = d_g+3
static inline void ff_poly_mod_delta3 (ff_t *f, ff_t *g, int d_g)
{
	register int i;
	register ff_t t0;
	
	_ff_mult(t0,f[d_g+3],g[d_g-1]); _ff_subfrom(f[d_g+2],t0);
	_ff_sum_2_mults_arr(t0,f+d_g+2,g+d_g-2); _ff_subfrom(f[d_g+1],t0);
	_ff_sum_3_mults_arr(t0,f+d_g+1,g+d_g-3); _ff_subfrom(f[d_g],t0);
	for ( i = d_g-1 ; i>2 ; i-- ) { _ff_sum_4_mults_arr(t0,f+d_g,g+i-3); _ff_subfrom(f[i],t0); }
	_ff_sum_3_mults_arr(t0,f+d_g,g); _ff_subfrom(f[2],t0);
	_ff_sum_2_mults_arr(t0,f+d_g,g); _ff_subfrom(f[1],t0);
	_ff_mult(t0,f[d_g],g[0]); _ff_subfrom(f[0],t0);
}

// assumes d_g >= 5, d_f = d_g+4
static inline void ff_poly_mod_delta4 (ff_t *f, ff_t *g, int d_g)
{
	register int i;
	register ff_t t0;
	
	_ff_mult(t0,f[d_g+4],g[d_g-1]); _ff_subfrom(f[d_g+3],t0);
	_ff_sum_2_mults_arr(t0,f+d_g+3,g+d_g-2); _ff_subfrom(f[d_g+2],t0);
	_ff_sum_3_mults_arr(t0,f+d_g+2,g+d_g-3); _ff_subfrom(f[d_g+1],t0);
	_ff_sum_4_mults_arr(t0,f+d_g+1,g+d_g-4); _ff_subfrom(f[d_g],t0);
	for ( i = d_g-1 ; i>3 ; i-- ) { _ff_sum_5_mults_arr(t0,f+d_g,g+i-4); _ff_subfrom(f[i],t0); }
	_ff_sum_4_mults_arr(t0,f+d_g,g); _ff_subfrom(f[3],t0);
	_ff_sum_3_mults_arr(t0,f+d_g,g); _ff_subfrom(f[2],t0);
	_ff_sum_2_mults_arr(t0,f+d_g,g); _ff_subfrom(f[1],t0);
	_ff_mult(t0,f[d_g],g[0]); _ff_subfrom(f[0],t0);
}

// assumes d_g >= 6, d_f = d_g+5
static inline void ff_poly_mod_delta5 (ff_t *f, ff_t *g, int d_g)
{
	register int i;
	register ff_t t0;
	
	_ff_mult(t0,f[d_g+5],g[d_g-1]); _ff_subfrom(f[d_g+4],t0);
	_ff_sum_2_mults_arr(t0,f+d_g+4,g+d_g-2); _ff_subfrom(f[d_g+3],t0);
	_ff_sum_3_mults_arr(t0,f+d_g+3,g+d_g-3); _ff_subfrom(f[d_g+2],t0);
	_ff_sum_4_mults_arr(t0,f+d_g+2,g+d_g-4); _ff_subfrom(f[d_g+1],t0);
	_ff_sum_5_mults_arr(t0,f+d_g+1,g+d_g-5); _ff_subfrom(f[d_g],t0);
	for ( i = d_g-1 ; i>4 ; i-- ) { _ff_sum_6_mults_arr(t0,f+d_g,g+i-5); _ff_subfrom(f[i],t0); }
	_ff_sum_5_mults_arr(t0,f+d_g,g); _ff_subfrom(f[4],t0);
	_ff_sum_4_mults_arr(t0,f+d_g,g); _ff_subfrom(f[3],t0);
	_ff_sum_3_mults_arr(t0,f+d_g,g); _ff_subfrom(f[2],t0);
	_ff_sum_2_mults_arr(t0,f+d_g,g); _ff_subfrom(f[1],t0);
	_ff_mult(t0,f[d_g],g[0]); _ff_subfrom(f[0],t0);
}

static inline void ff_poly_mod_d2 (ff_t *f, int d_f, ff_t *g)
{
	register int i;
	register ff_t t0;
	
	switch (d_f) {
		case 0: case 1: return;
		case 2: _ff_mult(t0,f[2],g[1]); _ff_subfrom(f[1],t0); _ff_mult(t0,f[2],g[0]); _ff_subfrom(f[0],t0); return;
	}
	_ff_mult(t0,f[d_f],g[1]); _ff_subfrom(f[d_f-1],t0);
	for ( i = d_f-2 ; i ; i-- ) { _ff_sum_2_mults_arr(t0,f+i+1,g); _ff_subfrom(f[i],t0); }
	_ff_mult(t0,f[2],g[0]); _ff_subfrom(f[0],t0);
}

static inline void ff_poly_mod_d3 (ff_t *f, int d_f, ff_t *g)
{
	register int i;
	register ff_t t0;
	
	switch (d_f) {
		case 0: case 1: case 2: return;
		case 3: _ff_mult(t0,f[3],g[2]); _ff_subfrom(f[2],t0); _ff_mult(t0,f[3],g[1]); _ff_subfrom(f[1],t0); _ff_mult(t0,f[3],g[0]); _ff_subfrom(f[0],t0); return;
	}
	_ff_mult(t0,f[d_f],g[2]); _ff_subfrom(f[d_f-1],t0);
	_ff_sum_2_mults_arr(t0,f+d_f-1,g+1); _ff_subfrom(f[d_f-2],t0);
	for ( i = d_f-3 ; i>1 ; i-- ) { _ff_sum_3_mults_arr(t0,f+i+1,g); _ff_subfrom(f[i],t0); }
	_ff_sum_2_mults_arr(t0,f+3,g); _ff_subfrom(f[1],t0);
	_ff_mult(t0,f[3],g[0]); _ff_subfrom(f[0],t0);
}

static inline void ff_poly_mod_d4 (ff_t *f, int d_f, ff_t *g)
{
	register int i;
	register ff_t t0;
	
	switch (d_f) {
		case 0: case 1: case 2: case 3: return;
		case 4: _ff_mult(t0,f[4],g[3]); _ff_subfrom(f[3],t0); _ff_mult(t0,f[4],g[2]); _ff_subfrom(f[2],t0); _ff_mult(t0,f[4],g[1]); _ff_subfrom(f[1],t0); _ff_mult(t0,f[4],g[0]); _ff_subfrom(f[0],t0); return;
		case 5: ff_poly_mod_delta1 (f, g, 4); return;
	}
	_ff_mult(t0,f[d_f],g[3]); _ff_subfrom(f[d_f-1],t0);
	_ff_sum_2_mults_arr(t0,f+d_f-1,g+2); _ff_subfrom(f[d_f-2],t0);
	_ff_sum_3_mults_arr(t0,f+d_f-2,g+1); _ff_subfrom(f[d_f-3],t0);
	for ( i = d_f-4 ; i>2 ; i-- ) { _ff_sum_4_mults_arr(t0,f+i+1,g); _ff_subfrom(f[i],t0); }
	_ff_sum_3_mults_arr(t0,f+4,g); _ff_subfrom(f[2],t0);
	_ff_sum_2_mults_arr(t0,f+4,g); _ff_subfrom(f[1],t0);
	_ff_mult(t0,f[4],g[0]); _ff_subfrom(f[0],t0);
}

static inline void ff_poly_mod_d5 (ff_t *f, int d_f, ff_t *g)
{
	register int i;
	register ff_t t0;
	
	switch (d_f) {
		case 0: case 1: case 2: case 3: case 4: return;
		case 5: _ff_mult(t0,f[5],g[4]); _ff_subfrom(f[4],t0); _ff_mult(t0,f[5],g[3]); _ff_subfrom(f[3],t0); _ff_mult(t0,f[5],g[2]); _ff_subfrom(f[2],t0); _ff_mult(t0,f[5],g[1]); _ff_subfrom(f[1],t0); _ff_mult(t0,f[5],g[0]); _ff_subfrom(f[0],t0); return;
		case 6: ff_poly_mod_delta1 (f, g, 5); return;
		case 7: ff_poly_mod_delta2 (f, g, 5); return;
		case 8: ff_poly_mod_delta3 (f, g, 5); return;
	}
	_ff_mult(t0,f[d_f],g[4]); _ff_subfrom(f[d_f-1],t0);
	_ff_sum_2_mults_arr(t0,f+d_f-1,g+3); _ff_subfrom(f[d_f-2],t0);
	_ff_sum_3_mults_arr(t0,f+d_f-2,g+2); _ff_subfrom(f[d_f-3],t0);
	_ff_sum_4_mults_arr(t0,f+d_f-3,g+1); _ff_subfrom(f[d_f-4],t0);
	for ( i = d_f-5 ; i>3 ; i-- ) { _ff_sum_5_mults_arr(t0,f+i+1,g); _ff_subfrom(f[i],t0); }
	_ff_sum_4_mults_arr(t0,f+5,g); _ff_subfrom(f[3],t0);
	_ff_sum_3_mults_arr(t0,f+5,g); _ff_subfrom(f[2],t0);
	_ff_sum_2_mults_arr(t0,f+5,g); _ff_subfrom(f[1],t0);
	_ff_mult(t0,f[5],g[0]); _ff_subfrom(f[0],t0);
}


static inline void ff_poly_mod_d6 (ff_t *f, int d_f, ff_t *g)
{
	register int i;
	register ff_t t0;
	
	switch (d_f) {
		case 0: case 1: case 2: case 3: case 4: case 5: return;
		case 6: _ff_mult(t0,f[6],g[5]); _ff_subfrom(f[5],t0); _ff_mult(t0,f[6],g[4]); _ff_subfrom(f[4],t0); _ff_mult(t0,f[6],g[3]); _ff_subfrom(f[3],t0); _ff_mult(t0,f[6],g[2]); _ff_subfrom(f[2],t0); _ff_mult(t0,f[6],g[1]); _ff_subfrom(f[1],t0); _ff_mult(t0,f[6],g[0]); _ff_subfrom(f[0],t0); return;
		case 7: ff_poly_mod_delta1 (f, g, 6); return;
		case 8: ff_poly_mod_delta2 (f, g, 6); return;
		case 9: ff_poly_mod_delta3 (f, g, 6); return;
		case 10: ff_poly_mod_delta4 (f, g, 6); return;
	}
	_ff_mult(t0,f[d_f],g[5]); _ff_subfrom(f[d_f-1],t0);
	_ff_sum_2_mults_arr(t0,f+d_f-1,g+4); _ff_subfrom(f[d_f-2],t0);
	_ff_sum_3_mults_arr(t0,f+d_f-2,g+3); _ff_subfrom(f[d_f-3],t0);
	_ff_sum_4_mults_arr(t0,f+d_f-3,g+2); _ff_subfrom(f[d_f-4],t0);
	_ff_sum_5_mults_arr(t0,f+d_f-4,g+1); _ff_subfrom(f[d_f-5],t0);
	for ( i = d_f-6 ; i>4 ; i-- ) { _ff_sum_6_mults_arr(t0,f+i+1,g); _ff_subfrom(f[i],t0); }
	_ff_sum_5_mults_arr(t0,f+6,g); _ff_subfrom(f[4],t0);
	_ff_sum_4_mults_arr(t0,f+6,g); _ff_subfrom(f[3],t0);
	_ff_sum_3_mults_arr(t0,f+6,g); _ff_subfrom(f[2],t0);
	_ff_sum_2_mults_arr(t0,f+6,g); _ff_subfrom(f[1],t0);
	_ff_mult(t0,f[6],g[0]); _ff_subfrom(f[0],t0);
}


// computes xf^2 mod g where g=x^n - g[n-2]x^{n-2} - ... - g[0] (note signs!)  and f has degree n-1 (n must be > 0)
static inline void ff_poly_square_mult_x_mod_small (ff_t o[], ff_t f[], ff_t g[], int n) 
{
	ff_t t[FF_POLY_SMALL_DEGREE+1];
	register int i;
	register ff_t t0, t1;

	n--;	// set n to degree of f and output
	_ff_square(t[n+1],f[n]);
	_ff_dbl_mult(t[n],f[n-1],f[n]);
	for ( i=n-2 ; i >= 0 ; i-- ) {
		_ff_sq_coeff (t0, f+i, n-i);
		_ff_sum_mults (t1, t+i+3, g+i+1, n-i-1);
		_ff_add(t[i+1],t0,t1);
	}
	_ff_sq_coeff (t0,f,n-1);
	_ff_sum_mults(t1,t+2,g,n);
	_ff_add(o[n],t0,t1);
	for ( i = n-2 ; i >= 0 ; i-- ) {
		_ff_sq_coeff (t0,f,i);
		_ff_sum_mults(t1,t+1,g,i+2);
		_ff_add(o[i+1],t0,t1);
	}
	_ff_mult(o[0],t[1],g[0]);
}

// For a monic poly f of degree 3, returns a linear multiple g1*x+g0 of gcd(f,f').
// This will be 0 if f has a triple root, and if f has a double root, then it is equal to -g0/g1.
// Uses 3M+9A
static inline void ff_poly_dbl_root_linfac_d3 (ff_t g[2], ff_t f[3])
{
	register ff_t t1, t2;
	
	_ff_triple(t1,f[1]);  _ff_square(t2,f[2]); _ff_subfrom(t2,t1); _ff_add(g[1],t2,t2);
	_ff_triple(t1,f[0]); _ff_triple(t2,t1); _ff_mult(t1,f[1],f[2]); _ff_sub(g[0],t1,t2);
}


/*
	The naming conventions for the inlines below refer to the number of coefficients, not the degree.
	For all the functions below, the output array may equal either input array.
	
	Note that we use far more register variables than we need to (indeed more than may exist), it seems to generate faster code.
	This code was tuned on an Athlon AMD 64 processor, YMMV.
	Loop unwinding is well worthwhile - typically 20% improvement or better.
*/


// 1A
static inline void ff_poly_add_1_2 (ff_t o[], ff_t f[], ff_t g[])
	{ _ff_add(o[0],f[0],g[0]); _ff_set(o[1],g[1]); }
	
// 2M
static inline void ff_poly_mult_1_2 (ff_t o[], ff_t f[], ff_t g[])
	{  _ff_mult(o[1],f[0],g[1]); _ff_mult(o[0],f[0],g[0]); }

// 2A
static inline void ff_poly_add_2_2 (ff_t o[], ff_t f[], ff_t g[])
	{ _ff_add(o[0],f[0],g[0]); _ff_add(o[1],f[1],g[1]); }

// 3M+1A
static inline void ff_poly_square_2o(ff_t o[], ff_t f[])
	{ register ff_t t; _ff_square(o[2],f[1]);  _ff_mult(t,f[0],f[1]); _ff_add(o[1],t,t); _ff_square(o[0],f[0]); }

// It is a bit of mystery why this is faster, but it is.  The order of the operations below has a noticeable impact on the speed (i.e don't mess with it)
static inline void ff_poly_square_2 (ff_t o[], ff_t f[])
	{ register ff_t t1; _ff_dbl_mult(t1,f[0],f[1]);  _ff_square(o[2],f[1]); _ff_square(o[0],f[0]); _ff_set(o[1],t1); }

// reduces f mod g=x^2-g1x-g0 in place (note signs, this is different from ff_poly_mod_d2!)
static inline void ff_poly_mod_2 (ff_t *f, int d_f, ff_t *g)
{
	register int i;
	register ff_t t0;
	
	switch (d_f) {
		case 0: case 1: return;
		case 2: _ff_mult(t0,f[2],g[1]); _ff_addto(f[1],t0); _ff_mult(t0,f[2],g[0]); _ff_addto(f[0],t0); return;
	}
	_ff_mult(t0,f[d_f],g[1]); _ff_addto(f[d_f-1],t0);
	for ( i = d_f-2 ; i ; i-- ) { _ff_sum_2_mults_arr(t0,f+i+1,g); _ff_addto(f[i],t0); }
	_ff_mult(t0,f[2],g[0]); _ff_addto(f[0],t0);
}		

// squares a linear poly modulo a monic quadratic x^2-g1x-g0 (note signs, g1 need not be zero), o may equal f
static inline void ff_poly_square_mod_2 (ff_t o[2], ff_t f[2], ff_t g[2])
{
	register ff_t t1, t2;
	
	_ff_square(t2,f[1]);
	_ff_sum_2_mults_d1(t1,f[0],t2,g[1],f[1]);
	_ff_sum_2_mults (o[0],f[0],t2,g[0],f[0]);
	_ff_set(o[1],t1);
}

// multiplies two linear polys modulo a monic quadratic x^2-g1x-g0 (note signs, g1 need not be zero), o may equal f
static inline void ff_poly_mult_mod_2 (ff_t o[2], ff_t f1[2], ff_t f2[2], ff_t g[2])
{
	register ff_t t1, t2;
	
	_ff_mult(t2,f1[1],f2[1]);
	_ff_sum_3_mults(t1,f1[0],f1[1],t2,g[1],f2[0],f2[1]);
	_ff_sum_2_mults (o[0],f1[0],t2,g[0],f2[0]);
	_ff_set(o[1],t1);
}

// 3M+4A 
static inline void ff_poly_mult_2_2o (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t0, t1, t2;

	_ff_add(t0,f[0],f[1]);  _ff_add(t1,g[0],g[1]);  _ff_mult(t2,t0,t1);
	_ff_mult(o[2],f[1],g[1]);  _ff_mult(o[0],f[0],g[0]);  _ff_subfrom(t2,o[0]); _ff_sub(o[1],t2,o[2]);
}

static inline void ff_poly_mult_2_2 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1;

	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_mult(o[0],f[0],g[0]);  _ff_mult(o[2],f[1],g[1]); _ff_set(o[1],t1);
}


// 2A
static inline void ff_poly_add_2_3 (ff_t o[], ff_t f[], ff_t g[])
	{ _ff_add(o[0],f[0],g[0]); _ff_add(o[1],f[1],g[1]); _ff_set(o[2],g[2]); }
	
// 5M+5A
static inline void ff_poly_mult_2_3o (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t;
	
	_ff_mult(t,f[0],g[2]);
	_ff_mult(o[3],f[1],g[2]);
	ff_poly_mult_2_2o(o,f,g);
	_ff_addto(o[2],t);
}

static inline void ff_poly_mult_2_3 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1,t2;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_2_mults(t2,f[0],f[1],g[1],g[2]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[3],f[1],g[2]);
	_ff_set(o[1],t1); _ff_set(o[2],t2);
}

// 3A
static inline void ff_poly_add_3_3 (ff_t o[], ff_t f[], ff_t g[])
	{ _ff_add(o[0],f[0],g[0]); _ff_add(o[1],f[1],g[1]); _ff_add(o[2],f[2],g[2]); }

// 6M+13A
static inline void ff_poly_mult_3_3o (ff_t o[], ff_t f[], ff_t g[])
{
	ff_t t0, t1, t2, t3, t4;
	
	_ff_add(t0,f[0],f[1]); _ff_add(t1,g[0],g[1]); _ff_mult(t2,t0,t1);
	_ff_add(t0,f[0],f[2]); _ff_add(t1,g[0],g[2]); _ff_mult(t3,t0,t1);
	_ff_add(t0,f[1],f[2]); _ff_add(t1,g[1],g[2]); _ff_mult(t4,t0,t1);
	_ff_mult(o[0],f[0],g[0]);  _ff_mult(o[2],f[1],g[1]); _ff_mult(o[4],f[2],g[2]);  
	_ff_subfrom(t2,o[0]); _ff_sub(o[1],t2,o[2]);
	_ff_subfrom(t4,o[4]); _ff_sub(o[3],t4,o[2]);
	_ff_subfrom(t3,o[4]); _ff_subfrom(t3,o[0]); _ff_addto(o[2],t3);
}

static inline void ff_poly_mult_3_3 (ff_t o[], ff_t f[], ff_t g[])
{
	ff_t t1, t2, t3;
	
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_2_mults(t3,f[1],f[2],g[1],g[2]);
	_ff_mult(o[0],f[0],g[0]);  _ff_mult(o[4],f[2],g[2]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3);
}

// 6M+4A
static inline void ff_poly_square_3o (ff_t o[], ff_t f[])
{
	register ff_t t;

	_ff_square(o[4],f[2]);
	_ff_mult(t,f[1],f[2]);  _ff_add(o[3],t,t);
	_ff_mult(t,f[0],f[2]); _ff_x2(t);
	ff_poly_square_2 (o,f);
	_ff_addto(o[2],t);
}

static inline void ff_poly_square_3 (ff_t o[], ff_t f[])
{
	register ff_t t1,t2,t3;

	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_dbl_mult(t3,f[1],f[2]);
	_ff_square(o[0],f[0]); _ff_square(o[4],f[2]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3);
}


// squares a quadratic modulo a depressed cubic of the form x^3-g1x-g0 (note the negative signs!), o may equal f
static inline void ff_poly_square_mod_3 (ff_t o[3], ff_t f[3], ff_t g[2])
{
	register ff_t t1,t2,t3,t4;

	_ff_square(t4,f[2]);									// x^4 coeff
	_ff_dbl_mult(t3,f[1],f[2]);								// x^3 coeff
	
	_ff_sum_3_mults_d1(t2,f[0],f[1],t4,g[1],f[1],f[2]);			// o[2] = 2f0f2+f1^2+t4g1
	_ff_sum_3_mults_d1(t1,f[0],t3,t4,g[0],g[1],f[1]);			// o[1] = 2f0f1+t3g1+t4g0
	_ff_sum_2_mults(o[0],f[0],t3,g[0],f[0]);					// o[0] = f0^2+t3g0
	_ff_set(o[1],t1); _ff_set(o[2],t2); 
}

// multiplies two quadratics modulo a depressed cubic of the form x^3-g1x-g0 (note the negative signs!), o may equal f1 or f2
static inline void ff_poly_mult_mod_3 (ff_t o[3], ff_t f1[3], ff_t f2[3], ff_t g[2])
{
	register ff_t t1,t2,t3,t4;

	_ff_mult(t4,f1[2],f2[2]);											// x^4 coeff
	_ff_sum_2_mults(t3,f1[1],f1[2],f2[1],f2[2]);							// x^3 coeff
	
	_ff_sum_4_mults(t2,f1[0],f1[1],f1[2],t4,g[1],f2[0],f2[1],f2[2]);			// o[2] = f10f22+f11f21+f12f20+t4g1
	_ff_sum_4_mults(t1,f1[0],f1[1],t3,t4,g[0],g[1],f2[0],f2[1]);				// o[1] = f10f21+f11f20+t3g1+t4g0
	_ff_sum_2_mults(o[0],f1[0],t3,g[0],f2[0]);							// o[0] = f10f20+t3g0
	_ff_set(o[1],t1); _ff_set(o[2],t2);
}


// Computes x*f^2 mod g for f degree 2, g=x^3-g1x-g0 (note signs)
static inline void ff_poly_square_mult_x_mod_3 (ff_t o[3], ff_t f[3], ff_t g[2])
{
	register ff_t t1,t2,t3,t4,t5;

	_ff_square(t5,f[2]);											// t5 = f2^2
	_ff_dbl_mult(t4,f[1],f[2]);										// t4 = f1f2

	_ff_sum_3_mults_d1(t3,f[0],f[1],t5,g[1],f[1],f[2]);					// t3 = 2f0f2+f1^2+t5g1
	_ff_sum_3_mults_d1 (t2,f[0],t4,t5,g[0],g[1],f[1]);					// o[2] = 2f0f1+t4g1+t5g0
	_ff_sum_3_mults (t1,f[0],t3,t4,g[0],g[1],f[0]);						// o[1] = f0^2+t3g1+t4g0 
	_ff_mult(o[0],t3,g[0]);										// o[0] = t3g0
	_ff_set(o[1],t1);  _ff_set(o[2],t2);
}


// 8M+18A
static inline void ff_poly_mult_3_4o (ff_t o[], ff_t f[], ff_t g[])
{
	ff_t u[2],v[2],w[3];
	
	ff_poly_add_1_2(u,f+2,f);
	ff_poly_add_2_2(v,g,g+2);
	ff_poly_mult_2_2(w,u,v);
	ff_poly_mult_1_2(o+4,f+2,g+2);
	ff_poly_mult_2_2(o,f,g);
	_ff_subfrom(w[2],o[2]);
	_ff_subfrom(w[1],o[5]); _ff_sub(o[3],w[1],o[1]);
	_ff_subfrom(w[0],o[4]); _ff_subfrom(w[0],o[0]);
	_ff_addto(o[4],w[2]); _ff_addto (o[2],w[0]);
}

static inline void ff_poly_mult_3_4 (ff_t o[], ff_t f[], ff_t g[])
{
	ff_t t1,t2,t3,t4;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_3_mults(t3,f[0],f[1],f[2],g[1],g[2],g[3]);
	_ff_sum_2_mults(t4,f[1],f[2],g[2],g[3]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[5],f[2],g[3]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4);
}

// 9M+11A
static inline void ff_poly_square_4o (ff_t o[], ff_t f[])
{
	register ff_t t0, t1, t2;
	
	_ff_add(t0,f[0],f[1]);  _ff_add(t1,f[2],f[3]);  _ff_mult(t2,t0,t1);
	_ff_mult(t0,f[0],f[2]); _ff_mult(t1,f[1],f[3]); _ff_subfrom(t2,t0); _ff_subfrom(t2,t1);
	ff_poly_square_2(o+4,f+2);
	ff_poly_square_2(o,f);
	_ff_x2(t0); _ff_addto(o[2],t0);
	_ff_add(o[3],t2,t2);	
	_ff_x2(t1); _ff_addto(o[4],t1);
}

static inline void ff_poly_square_4 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3;

	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_3_mults_s3(o[4],f[1],f[2],f[3]);
	_ff_dbl_mult(o[5],f[2],f[3]);
	_ff_square(o[0],f[0]); _ff_square(o[6],f[3]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3);
}


// squares a cubic mod a depressed quartic of the form x^4-g2x^2-g1x-g0 (note signs!)
static inline void ff_poly_square_mod_4 (ff_t o[4], ff_t f[4], ff_t g[3])
{
	register ff_t t1,t2,t3,t4,t5,t6;

	_ff_square(t6,f[3]);											// t6 = f3^2
	_ff_dbl_mult(t5,f[2],f[3]);										// t5 = f2f3
	_ff_sum_3_mults_d1 (t4,f[1],f[2],t6,g[2],f[2],f[3]);					// t4 = 2f1f3+f2^2+t6g2
	
	_ff_sum_4_mults_d2(t3,f[0],f[1],t5,t6,g[1],g[2],f[2],f[3]);			// o[3] = 2f0f3+2f1f2+t5g2+t6g1
	_ff_sum_5_mults_d1(t2,f[0],f[1],t4,t5,t6,g[0],g[1],g[2],f[1],f[2]);		// o[2] = 2f0f2+f1^2+t4g2+t5g1+t6g0
	_ff_sum_3_mults_d1(t1,f[0],t4,t5,g[0],g[1],f[1]);					// o[1] = 2f0f1+t4g1+t5g0
	_ff_sum_2_mults(o[0],f[0],t4,g[0],f[0]);							// o[0] = f0^2+t4g0
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3);
}


// multiplies two cubics mod a depressed quartic of the form x^4-g2x^2-g1x-g0 (note signs!)
static inline void ff_poly_mult_mod_4 (ff_t o[4], ff_t f1[4], ff_t f2[4], ff_t g[3])
{
	register ff_t t1,t2,t3,t4,t5,t6;

	_ff_mult(t6,f1[3],f2[3]);										// x^6 coeff
	_ff_sum_2_mults(t5,f1[2],f1[3],f2[2],f2[3]);						// x^5 coeff
	_ff_sum_4_mults(t4,f1[1],f1[2],f1[3],t6,g[2],f2[1],f2[2],f2[3]);		// x^4 coeff
	_ff_sum_6_mults(t3,f1[0],f1[1],f1[2],f1[3],t5,t6,g[1],g[2],f2[0],f2[1],f2[2],f2[3]);
	_ff_sum_6_mults(t2,f1[0],f1[1],f1[2],t4,t5,t6,g[0],g[1],g[2],f2[0],f2[1],f2[2]);
	_ff_sum_4_mults(t1,f1[0],f1[1],t4,t5,g[0],g[1],f2[0],f2[1]);
	_ff_sum_2_mults(o[0],f1[0],t4,g[0],f2[0]);
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3);
}


// Computes x*f^2 mod g for f degree 3, g=x^4-g2x^2-g1x-g0 (note signs)
static inline void ff_poly_square_mult_x_mod_4 (ff_t o[4], ff_t f[4], ff_t g[3])
{
	register ff_t t1,t2,t3,t4,t5,t6, t7;

	_ff_square(t7,f[3]);											// t7 = f3^2
	_ff_dbl_mult(t6,f[2],f[3]);										// t6 = f2f3
	_ff_sum_3_mults_d1 (t5,f[1],f[2],t7,g[2],f[2],f[3]);					// t5 = 2f1f3+f2^2+t7g2
	_ff_sum_4_mults_d2(t4,f[0],f[1],t6,t7,g[1],g[2],f[2],f[3]);			// t4 = 2f0f3+2f1f2+t6g2+t7g1
	_ff_sum_5_mults_d1(t3,f[0],f[1],t5,t6,t7,g[0],g[1],g[2],f[1],f[2]);		// o[3] = 2f0f2+f1^2+t5fg2+t6t1+t7g0
	_ff_sum_4_mults_d1 (t2,f[0],t4,t5,t6,g[0],g[1],g[2],f[1]);			// o[2] = 2f0f1+t6g0+t5g1+t4g2
	_ff_sum_3_mults (t1,f[0],t4,t5,g[0],g[1],f[0]);						// o[1] = f0^2+t4g1+t5g0
	_ff_mult(o[0],t4,g[0]);										// o[0] = t4g0
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3);
}


// 9M+24A
static inline void ff_poly_mult_4_4o (ff_t o[], ff_t f[], ff_t g[])
{
	ff_t u[2],v[2],w[3];
	
	ff_poly_add_2_2 (u,f,f+2); ff_poly_add_2_2(v,g,g+2);
	ff_poly_mult_2_2(w,u,v);
	ff_poly_mult_2_2(o+4,f+2,g+2);
	ff_poly_mult_2_2(o,f,g);
	_ff_subfrom(w[2],o[6]); _ff_subfrom(w[2],o[2]);
	_ff_subfrom(w[1],o[5]); _ff_sub(o[3],w[1],o[1]);
	_ff_subfrom(w[0],o[4]); _ff_subfrom(w[0],o[0]);
	_ff_addto(o[4],w[2]); _ff_addto (o[2],w[0]);
}

static inline void ff_poly_mult_4_4 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_3_mults(o[4],f[1],f[2],f[3],g[1],g[2],g[3]);
	_ff_sum_2_mults(o[5],f[2],f[3],g[2],g[3]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[6],f[3],g[3]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3);
}

// 13M+28A
static inline void ff_poly_mult_4_5o (ff_t o[], ff_t f[], ff_t g[])
{
	ff_t u[2],v[3],w[5];

	ff_poly_add_2_2(u,f,f+2);  ff_poly_add_2_2(v,g,g+2); _ff_set(v[2],g[4]);
	ff_poly_mult_2_3(w,u,v);
	ff_poly_mult_2_3(o+4,f+2,g+2);
	ff_poly_mult_2_2(o,f,g);
	_ff_subfrom(w[3],o[7]);
	_ff_subfrom(w[2],o[6]); _ff_subfrom(w[2],o[2]);
	_ff_subfrom(w[1],o[5]); _ff_sub(o[3],w[1],o[1]);
	_ff_subfrom(w[0],o[4]); _ff_subfrom(w[0],o[0]);
	_ff_addto(o[5],w[3]);_ff_addto(o[4],w[2]);
	_ff_addto (o[2],w[0]);
}

// 13M+28A
static inline void ff_poly_mult_4_5 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_4_mults(t4,f[0],f[1],f[2],f[3],g[1],g[2],g[3],g[4]);
	_ff_sum_3_mults(o[5],f[1],f[2],f[3],g[2],g[3],g[4]);
	_ff_sum_2_mults(o[6],f[2],f[3],g[3],g[4]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[7],f[3],g[4]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4);
}

// 14M+17A
static inline void ff_poly_square_5o (ff_t o[], ff_t f[])
{
	ff_t w[4];
	register ff_t t0;
	
	ff_poly_mult_2_3(w,f+3,f);
	ff_poly_square_2(o+6,f+3);
	ff_poly_square_3(o,f);
	_ff_add(t0,w[0],w[0]); _ff_addto(o[3],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[4],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[6],t0);
	_ff_add(o[5],w[2],w[2]);
}

static inline void ff_poly_square_5 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_4_mults_s4(o[5],f[1],f[2],f[3],f[4]);
	_ff_sum_3_mults_s3(o[6],f[2],f[3],f[4]);
	_ff_dbl_mult(o[7],f[3],f[4]);
	_ff_square(o[0],f[0]); _ff_square(o[8],f[4]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4);
}

// computes f^2 mod g  for f degree 4 and g(x) = x^5-g3x^3-g2x^2-g1x-g0 (note signs!)
static inline void ff_poly_square_mod_5 (ff_t o[5], ff_t f[5], ff_t g[4])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8;

	_ff_square(t8,f[4]);													// t8 = f4^2
	_ff_dbl_mult(t7,f[3],f[4]);												// t7 = 2f3f4
	_ff_sum_3_mults_d1 (t6,f[2],f[3],t8,g[3],f[3],f[4]);							// t6 = 2f2f4+f3^2+t8g3
	_ff_sum_4_mults_d2(t5,f[1],f[2],t7,t8,g[2],g[3],f[3],f[4]);					// t5 = 2f1f4+2f2f3+t7g3+t8g2
	
	_ff_sum_6_mults_d2(t4,f[0],f[1],f[2],t6,t7,t8,g[1],g[2],g[3],f[2],f[3],f[4]);		// o[4] = 2f0f4+2f1f3+f2^2+t6g3+t7g2+t8g1
	_ff_sum_6_mults_d2(t3,f[0],f[1],t5,t6,t7,t8,g[0],g[1],g[2],g[3],f[2],f[3]);		// o[3] = 2f0f3+2f1f2+t5g3+t6g2+t7t1+t8g0
	_ff_sum_5_mults_d1(t2,f[0],f[1],t5,t6,t7,g[0],g[1],g[2],f[1],f[2]);				// o[2] = 2f0f2+f1^2+5g2+t6g1+t7g0
	_ff_sum_3_mults_d1(t1,f[0],t5,t6,g[0],g[1],f[1]);							// o[1] = 2f0f1+t5g1+t6g0
	_ff_sum_2_mults(o[0],f[0],t5,g[0],f[0]);									// o[0] = f0^2+t5g0
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4);
}

// Computes x*f^2 mod g for f degree 4, g=x^5-g3x^3-g2x^2-g1x-g0 (note signs)
static inline void ff_poly_square_mult_x_mod_5 (ff_t o[5], ff_t f[5], ff_t g[4])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9;

	_ff_square(t9,f[4]);													// t9 = f4^2
	_ff_dbl_mult(t8,f[3],f[4]);												// t8 = 2f3f4
	_ff_sum_3_mults_d1 (t7,f[2],f[3],t9,g[3],f[3],f[4]);							// t7 = 2f2f4+f3^2+t9g3
	_ff_sum_4_mults_d2(t6,f[1],f[2],t8,t9,g[2],g[3],f[3],f[4]);					// t6 = 2f1f4+2f2f3+t8g3+t9g2
	_ff_sum_6_mults_d2(t5,f[0],f[1],f[2],t7,t8,t9,g[1],g[2],g[3],f[2],f[3],f[4]);		// t5 = 2f0f4+2f1f3+f2^2+t7g3+t8g2+t9g1
	
	_ff_sum_6_mults_d2(t4,f[0],f[1],t6,t7,t8,t9,g[0],g[1],g[2],g[3],f[2],f[3]);		// o[4] = 2f0f3+2f1f2+t6g3+t7g2+t8t1+t9g0
	_ff_sum_6_mults_d1(t3,f[0],f[1],t5,t6,t7,t8,g[0],g[1],g[2],g[3],f[1],f[2]);		// o[3] = 2f0f2+f1^2+t5g3+t6g2+t7g1+t8g0
	_ff_sum_4_mults_d1(t2,f[0],t5,t6,t7,g[0],g[1],g[2],f[1]);					// o[2] = 2f0f1+t5g2+t6g1+t7g0
	_ff_sum_3_mults(t1,f[0],t5,t6,g[0],g[1],f[0]);								// o[1] = f0^2+t5g1+t6g0
	_ff_mult(o[0],t5,g[0]);
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4);
}


/*
// 14M+44A (it is possible to do this with 13M--see Montgomery 567 Karatsuba paper--but it requires many more additions)
static inline void ff_poly_mult_5_5o (ff_t o[], ff_t f[], ff_t g[])
{
	ff_t u[3],v[3],w[5];
	
	_ff_set(u[0],f[2]); ff_poly_add_2_2 (u+1,f,f+3); _ff_set(v[0],g[2]); ff_poly_add_2_2(v+1,g,g+3);		// u = f1+xf0 where f=f1(x)x^2+f0(x), ditto for v and g
	ff_poly_mult_3_3(o+4,f+2,g+2);														
	ff_poly_mult_2_2(o,f,g);																
	ff_poly_mult_3_3_c(w,u,v,o[4]);
	_ff_subfrom(w[4],o[8]); _ff_subfrom(w[4],o[2]);
	_ff_subfrom(w[3],o[7]); _ff_subfrom(w[3],o[1]);
	_ff_subfrom(w[2],o[6]); _ff_sub(o[3],w[2],o[0]);
	_ff_subfrom(w[1],o[5]);
	_ff_addto(o[5],w[4]); _ff_addto(o[4],w[3]);
	_ff_addto (o[2],w[1]);
}
*/

static inline void ff_poly_mult_5_5 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_4_mults(o[5],f[1],f[2],f[3],f[4],g[1],g[2],g[3],g[4]);
	_ff_sum_3_mults(o[6],f[2],f[3],f[4],g[2],g[3],g[4]);
	_ff_sum_2_mults(o[7],f[3],f[4],g[3],g[4]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[8],f[4],g[4]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4);
}

// 17M+49A 
static inline void ff_poly_mult_5_6o (ff_t o[], ff_t f[], ff_t g[])
{
	ff_t u[3],v[3],w[5];

	ff_poly_add_2_2 (u,f+3,f); _ff_set(u[2],f[2]); ff_poly_add_3_3(v,g,g+3);
	ff_poly_mult_3_3(w,u,v);
	ff_poly_mult_2_3(o+6,f+3,g+3);
	ff_poly_mult_3_3(o,f,g);
	_ff_subfrom(w[4],o[4]);
	_ff_subfrom(w[3],o[9]); _ff_subfrom(w[3],o[3]);
	_ff_subfrom(w[1],o[7]); _ff_subfrom(w[1],o[1]);
	_ff_subfrom(w[0],o[6]); _ff_subfrom(w[0],o[0]);
	_ff_subfrom(w[2],o[8]); _ff_sub(o[5],w[2],o[2]);
	_ff_addto(o[7],w[4]); _ff_addto(o[6],w[3]);
	_ff_addto(o[4],w[1]); _ff_addto (o[3],w[0]);
}

static inline void ff_poly_mult_5_6 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4, t5;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_5_mults(t5,f[0],f[1],f[2],f[3],f[4],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_4_mults(o[6],f[1],f[2],f[3],f[4],g[2],g[3],g[4],g[5]);
	_ff_sum_3_mults(o[7],f[2],f[3],f[4],g[3],g[4],g[5]);
	_ff_sum_2_mults(o[8],f[3],f[4],g[4],g[5]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[9],f[4],g[5]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5);
}

// 18M+30A
static inline void ff_poly_square_6o (ff_t o[], ff_t f[])
{
	ff_t w[5];
	register ff_t t0;
	
	ff_poly_mult_3_3(w,f+3,f);
	ff_poly_square_3(o+6,f+3);
	ff_poly_square_3(o,f);
	_ff_add(t0,w[0],w[0]); _ff_addto(o[3],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[4],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[6],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[7],t0);
	_ff_add(o[5],w[2],w[2]);
}

static inline void ff_poly_square_6 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_5_mults_s5(o[6],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_4_mults_s4(o[7],f[2],f[3],f[4],f[5]);
	_ff_sum_3_mults_s3(o[8],f[3],f[4],f[5]);
	_ff_dbl_mult(o[9],f[4],f[5]);
	_ff_square(o[0],f[0]); _ff_square(o[10],f[5]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5);
}

// computes f^2 mod g  for f degree 5 and g(x) = x^6-g4x^4-g3x^3-g2x^2-g1x-g0 (note signs!)
static inline void ff_poly_square_mod_6 (ff_t o[6], ff_t f[6], ff_t g[5])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;

	_ff_square(t10,f[5]);																// t10 = f5^2
	_ff_dbl_mult(t9,f[4],f[5]);																// t9 = 2f4f5
	_ff_sum_3_mults_d1 (t8,f[3],f[4],t10,g[4],f[4],f[5]);										// t8 = 2f3f5+f4^2+t10g4
	_ff_sum_4_mults_d2(t7,f[2],f[3],t9,t10,g[3],g[4],f[4],f[5]);									// t7 = 2f2f5+2f3f4+t9g4+t10g3
	_ff_sum_6_mults_d2(t6,f[1],f[2],f[3],t8,t9,t10,g[2],g[3],g[4],f[3],f[4],f[5]);						// t6 = 2f1f5+2f2f4+f3^2+t8g4+t9g3+t10g2
	
	_ff_sum_7_mults_d3(t5,f[0],f[1],f[2],t7,t8,t9,t10,g[1],g[2],g[3],g[4],f[3],f[4],f[5]);				// o[5] = 2f0f5+2f1f4+2f2f3+t7g4+t8g3+t9g3+t10g1
	_ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t6,t7,t8,t9,t10,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);			// o[4] = 2f0f4+2f1f3+f2^2+t6g4+t7g3+t8g2+t9g1+t10g0
	_ff_sum_6_mults_d2(t3,f[0],f[1],t6,t7,t8,t9,g[0],g[1],g[2],g[3],f[2],f[3]);						// o[3] = 2f0f3+2f1f2+t6g3+t7t2+t8g1+t9g0
	_ff_sum_5_mults_d1(t2,f[0],f[1],t6,t7,t8,g[0],g[1],g[2],f[1],f[2]);								// o[2] = 2f0f2+f1^2+t6g2+t7g1+t8g0
	_ff_sum_3_mults_d1(t1,f[0],t6,t7,g[0],g[1],f[1]);											// o[1] = 2f0f1+t6g1+t7g0
	_ff_sum_2_mults(o[0],f[0],t6,g[0],f[0]);													// o[0] = f0^2+t6g0
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5);
}

// Computes x*f^2 mod g for f degree 5, g=x^6--g4x^4-g3x^3g2x^2-g1x-g0 (note signs)
static inline void ff_poly_square_mult_x_mod_6 (ff_t o[6], ff_t f[6], ff_t g[5])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11;

	_ff_square(t11,f[5]);																// t11 = f5^2
	_ff_dbl_mult(t10,f[4],f[5]);															// t10 = 2f4f5
	_ff_sum_3_mults_d1 (t9,f[3],f[4],t11,g[4],f[4],f[5]);										// t9 = 2f3f5+f4^2+t11g4
	_ff_sum_4_mults_d2(t8,f[2],f[3],t10,t11,g[3],g[4],f[4],f[5]);									// t8 = 2f2f5+2f3f4+t10g4+t11g3
	_ff_sum_6_mults_d2(t7,f[1],f[2],f[3],t9,t10,t11,g[2],g[3],g[4],f[3],f[4],f[5]);						// t7 = 2f1f5+2f2f4+f3^2+t9g4+t10g3+t11g2
	_ff_sum_7_mults_d3(t6,f[0],f[1],f[2],t8,t9,t10,t11,g[1],g[2],g[3],g[4],f[3],f[4],f[5]);				// t6 = 2f0f5+2f1f4+2f2f3+t8g4+t9g3+t10g2+t11g1
	
	_ff_sum_8_mults_d2(t5,f[0],f[1],f[2],t7,t8,t9,t10,t11,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);			// o[5] = 2f0f4+2f1f3+f2^2+t7g4+t8g3+t9g2+t10g1+t11g0
	_ff_sum_7_mults_d2(t4,f[0],f[1],t6,t7,t8,t9,t10,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);				// o[4] = 2f0f3+2f1f2+t6g4+t7g3+t8g2+t9t1+t10g0
	_ff_sum_6_mults_d1(t3,f[0],f[1],t6,t7,t8,t9,g[0],g[1],g[2],g[3],f[1],f[2]);						// o[3] = 2f0f2+f1^2+t6g3+t7g2+t8g1+t9g0
	_ff_sum_4_mults_d1(t2,f[0],t6,t7,t8,g[0],g[1],g[2],f[1]);									// o[2] = 2f0f1+t6g2+t7g1+t8g0
	_ff_sum_3_mults(t1,f[0],t6,t7,g[0],g[1],f[0]);												// o[1] = f0^2+t6g1+t7g0
	_ff_mult(o[0],t6,g[0]);
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5);
}

// multiplies two quintics mod a depressed sextic of the form x^6-g4x^4-g3x^3-g2x^2-g1x-g0 (note signs!)
static inline void ff_poly_mult_mod_6 (ff_t o[5], ff_t f1[6], ff_t f2[6], ff_t g[5])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;

	_ff_mult(t10,f1[5],f2[5]);
	_ff_sum_2_mults(t9,f1[4],f1[5],f2[4],f2[5]);
	_ff_sum_4_mults(t8,f1[3],f1[4],f1[5],t10,g[4],f2[3],f2[4],f2[5]);
	_ff_sum_6_mults(t7,f1[2],f1[3],f1[4],f1[5],t9,t10,g[3],g[4],f2[2],f2[3],f2[4],f2[5]);
	_ff_sum_8_mults(t6,f1[1],f1[2],f1[3],f1[4],f1[5],t8,t9,t10,g[2],g[3],g[4],f2[1],f2[2],f2[3],f2[4],f2[5]);
	_ff_sum_10_mults(t5,f1[0],f1[1],f1[2],f1[3],f1[4],f1[5],t7,t8,t9,t10,g[1],g[2],g[3],g[4],f2[0],f2[1],f2[2],f2[3],f2[4],f2[5]);
	_ff_sum_10_mults(t4,f1[0],f1[1],f1[2],f1[3],f1[4],t6,t7,t8,t9,t10,g[0],g[1],g[2],g[3],g[4],f2[0],f2[1],f2[2],f2[3],f2[4]);
	_ff_sum_8_mults(t3,f1[0],f1[1],f1[2],f1[3],t6,t7,t8,t9,g[0],g[1],g[2],g[3],f2[0],f2[1],f2[2],f2[3]);
	_ff_sum_6_mults(t2,f1[0],f1[1],f1[2],t6,t7,t8,g[0],g[1],g[2],f2[0],f2[1],f2[2]);
	_ff_sum_4_mults(t1,f1[0],f1[1],t6,t7,g[0],g[1],f2[0],f2[1]);
	_ff_sum_2_mults(o[0],f1[0],t6,g[0],f2[0]);
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5);
}

// 18M+59A
static inline void ff_poly_mult_6_6o (ff_t o[], ff_t f[], ff_t g[])
{
	ff_t u[3],v[3],w[5];
	
	ff_poly_add_3_3 (u,f,f+3); ff_poly_add_3_3(v,g,g+3);
	ff_poly_mult_3_3(w,u,v);
	ff_poly_mult_3_3(o+6,f+3,g+3);
	ff_poly_mult_3_3(o,f,g);
	_ff_subfrom(w[4],o[10]); _ff_subfrom(w[4],o[4]);
	_ff_subfrom(w[3],o[9]); _ff_subfrom(w[3],o[3]);
	_ff_subfrom(w[1],o[7]); _ff_subfrom(w[1],o[1]);
	_ff_subfrom(w[0],o[6]); _ff_subfrom(w[0],o[0]);
	_ff_subfrom(w[2],o[8]); _ff_sub(o[5],w[2],o[2]);
	_ff_addto(o[4],w[1]); _ff_addto (o[3],w[0]);
	_ff_addto(o[7],w[4]); _ff_addto (o[6],w[3]);
}

static inline void ff_poly_mult_6_6 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4,t5;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_6_mults(t5,f[0],f[1],f[2],f[3],f[4],f[5],g[0],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_5_mults(o[6],f[1],f[2],f[3],f[4],f[5],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_4_mults(o[7],f[2],f[3],f[4],f[5],g[2],g[3],g[4],g[5]);
	_ff_sum_3_mults(o[8],f[3],f[4],f[5],g[3],g[4],g[5]);
	_ff_sum_2_mults(o[9],f[4],f[5],g[4],g[5]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[10],f[5],g[5]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5);
}

// 22M+71A
static inline void ff_poly_mult_6_7o (ff_t o[], ff_t f[], ff_t g[])
{
	ff_t u[3],v[4],w[6];

	ff_poly_add_3_3(u,f,f+3);  ff_poly_add_3_3(v,g,g+3); _ff_set(v[3],g[6]);
	ff_poly_mult_3_4(w,u,v);
	ff_poly_mult_3_4(o+6,f+3,g+3);
	ff_poly_mult_3_3(o,f,g);
	_ff_subfrom(w[5],o[11]);
	_ff_subfrom(w[4],o[10]); _ff_subfrom(w[4],o[4]);
	_ff_subfrom(w[3],o[9]); _ff_subfrom(w[3],o[3]);
	_ff_subfrom(w[2],o[8]); _ff_sub(o[5],w[2],o[2]);
	_ff_subfrom(w[1],o[7]); _ff_subfrom(w[1],o[1]);
	_ff_subfrom(w[0],o[6]); _ff_subfrom(w[0],o[0]);
	_ff_addto(o[8],w[5]); _ff_addto(o[7],w[4]); 
	_ff_addto(o[6],w[3]); _ff_addto(o[4],w[1]); 
	_ff_addto(o[3],w[0]);
}

static inline void ff_poly_mult_6_7 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4,t5, t6;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_6_mults(t5,f[0],f[1],f[2],f[3],f[4],f[5],g[0],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_6_mults(t6,f[0],f[1],f[2],f[3],f[4],f[5],g[1],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_5_mults(o[7],f[1],f[2],f[3],f[4],f[5],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_4_mults(o[8],f[2],f[3],f[4],f[5],g[3],g[4],g[5],g[6]);
	_ff_sum_3_mults(o[9],f[3],f[4],f[5],g[4],g[5],g[6]);
	_ff_sum_2_mults(o[10],f[4],f[5],g[5],g[6]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[11],f[5],g[6]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6);
}

// 23M+44A
static inline void ff_poly_square_7o (ff_t o[], ff_t f[])
{
	ff_t w[6];
	register ff_t t0;
	
	ff_poly_mult_3_4(w,f+4,f);
	ff_poly_square_3(o+8,f+4);
	ff_poly_square_4(o,f);
	// it is still worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[4],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[5],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[6],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[8],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[9],t0);
	_ff_add(o[7],w[3],w[3]);
}

static inline void ff_poly_square_7 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_6_mults_s6(o[7],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_5_mults_s5(o[8],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_4_mults_s4(o[9],f[3],f[4],f[5],f[6]);
	_ff_sum_3_mults_s3(o[10],f[4],f[5],f[6]);
	_ff_dbl_mult(o[11],f[5],f[6]);
	_ff_square(o[0],f[0]); _ff_square(o[12],f[6]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6);
}

// computes f^2 mod g  for f degree 6 and g(x) = x^7-g5x^5-g4x^4-g3x^3-g2x^2-g1x-g0 (note signs!)
static inline void ff_poly_square_mod_7 (ff_t o[7], ff_t f[7], ff_t g[6])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12;

	_ff_square(t12,f[6]);																// t12 = f6^2
	_ff_dbl_mult(t11,f[5],f[6]);															// t11 = 2f5f6
	_ff_sum_3_mults_d1 (t10,f[4],f[5],t12,g[5],f[5],f[6]);										// t10 = 2f4f6+f5^2+t12g5
	_ff_sum_4_mults_d2(t9,f[3],f[4],t11,t12,g[4],g[5],f[5],f[6]);									// t9 = 2f3f6+2f4f5+t11g5+t12g4
	_ff_sum_6_mults_d2(t8,f[2],f[3],f[4],t10,t11,t12,g[3],g[4],g[5],f[4],f[5],f[6]);					// t8 = 2f2f6+2f3f5+f4^2+t10g5+t11g4+t12g3
	_ff_sum_7_mults_d3(t7,f[1],f[2],f[3],t9,t10,t11,t12,g[2],g[3],g[4],g[5],f[4],f[5],f[6]);				// t7 = 2f1f6+2f2f5+2f3f4+t9g5+t10g4+t11g3+t12g2
	
	_ff_sum_9_mults_d3(t6,f[0],f[1],f[2],f[3],t8,t9,t10,t11,t12,g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5],f[6]);	// o[6] = 2f0f6+2f1f5+2f2f4+f3^2+t8g5+t9g4+t10g3+t11g2+t12g1
	_ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t7,t8,t9,t10,t11,t12,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);	// o[5] = 2f0f5+2f1f4+2f2f3+t7g5+t8g4+t9g3+t10g2+t11g1+t12g0
	_ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t7,t8,t9,t10,t11,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);			// o[4] = 2f0f4+2f1f3+f2^2+t7g4+t8g3+t9g2+t10g1+t11g0
	_ff_sum_6_mults_d2(t3,f[0],f[1],t7,t8,t9,t10,g[0],g[1],g[2],g[3],f[2],f[3]);						// o[3] = 2f0f3+2f1f2+t7t3+t8g2+t9g1+t10g0
	_ff_sum_5_mults_d1(t2,f[0],f[1],t7,t8,t9,g[0],g[1],g[2],f[1],f[2]);								// o[2] = 2f0f2+f1^2+t7g2+t8g1+t9g0
	_ff_sum_3_mults_d1(t1,f[0],t7,t8,g[0],g[1],f[1]);											// o[1] = 2f0f1+t7g1+t8g0
	_ff_sum_2_mults(o[0],f[0],t7,g[0],f[0]);													// o[0] = f0^2+t7g0
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6);
}

// Computes x*f^2 mod g for f degree 6, g=x^7-g5x^5-g4x^4-g3x^3g2x^2-g1x-g0 (note signs)
static inline void ff_poly_square_mult_x_mod_7 (ff_t o[7], ff_t f[7], ff_t g[6])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13;

	_ff_square(t13,f[6]);																	// t13 = f6^2
	_ff_dbl_mult(t12,f[5],f[6]);																// t12 = 2f5f6
	_ff_sum_3_mults_d1 (t11,f[4],f[5],t13,g[5],f[5],f[6]);											// t11 = 2f4f6+f5^2+t13g5
	_ff_sum_4_mults_d2(t10,f[3],f[4],t12,t13,g[4],g[5],f[5],f[6]);										// t10 = 2f3f6+2f4f5+t12g5+t13g4
	_ff_sum_6_mults_d2(t9,f[2],f[3],f[4],t11,t12,t13,g[3],g[4],g[5],f[4],f[5],f[6]);						// t9 = 2f2f6+2f3f5+f4^2+t11g5+t12g4+t13g3
	_ff_sum_7_mults_d3(t8,f[1],f[2],f[3],t10,t11,t12,t13,g[2],g[3],g[4],g[5],f[4],f[5],f[6]);					// t8 = 2f1f6+2f2f5+2f3f4+t10g5+t11g4+t12g3+t13g2
	_ff_sum_9_mults_d3(t7,f[0],f[1],f[2],f[3],t9,t10,t11,t12,t13,g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5],f[6]);	// t7 = 2f0f6+2f1f5+2f2f4+f5^2+t9g5+t10g4+t11g3+t12g2+t13g1
	
	_ff_sum_9_mults_d3(t6,f[0],f[1],f[2],t8,t9,t10,t11,t12,t13,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);	// o[6] = 2f0f5+2f1f4+2f2f3+t8g5+t9g4+t10g3+t11g2+t12g1+t13g0
	_ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t7,t8,t9,t10,t11,t12,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);		// o[5] = 2f0f4+2f1f3+f2^2+t7g5+t8g4+t9g3+t10g2+t11g1+t12g0
	_ff_sum_7_mults_d2(t4,f[0],f[1],t7,t8,t9,t10,t11,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);					// o[4] = 2f0f3+2f1f2+t7g4+t8g3+t9g2+t10g1+t11g0
	_ff_sum_6_mults_d1(t3,f[0],f[1],t7,t8,t9,t10,g[0],g[1],g[2],g[3],f[1],f[2]);							// o[3] = 2f0f2+f1^2+t7g3+t8g2+t9g1+t10g0
	_ff_sum_4_mults_d1(t2,f[0],t7,t8,t9,g[0],g[1],g[2],f[1]);										// o[2] = 2f0f1+t7g2+t8g1+t9g0
	_ff_sum_3_mults(t1,f[0],t7,t8,g[0],g[1],f[0]);													// o[1] = f0^2+t7g1+t8g0
	_ff_mult(o[0],t7,g[0]);																	// o[0] = t7g0
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6);
}

// 24M+85A (it is possible to do this with 23M or even 22M--see Montgomery 567 Karatsuba paper).
static inline void ff_poly_mult_7_7o (ff_t o[], ff_t f[], ff_t g[])
{
	ff_t u[4],v[4],w[7];
	
	ff_poly_add_3_3 (u,f,f+3); _ff_set(u[3],f[6]); ff_poly_add_3_3(v,g,g+3); _ff_set(v[3],g[6]);
	ff_poly_mult_4_4(w,u,v);
	ff_poly_mult_4_4(o+6,f+3,g+3);
	ff_poly_mult_3_3(o,f,g);
	_ff_subfrom(w[6],o[12]);
	_ff_subfrom(w[5],o[11]);
	_ff_subfrom(w[4],o[10]); _ff_subfrom(w[4],o[4]);
	_ff_subfrom(w[3],o[9]); _ff_subfrom(w[3],o[3]);
	_ff_subfrom(w[2],o[8]); _ff_sub(o[5],w[2],o[2]);
	_ff_subfrom(w[1],o[7]); _ff_subfrom(w[1],o[1]);
	_ff_subfrom(w[0],o[6]); _ff_subfrom(w[0],o[0]);
	_ff_addto(o[4],w[1]); _ff_addto (o[3],w[0]);
	_ff_addto(o[7],w[4]); _ff_addto (o[6],w[3]);
	_ff_addto(o[9],w[6]); _ff_addto (o[8],w[5]);
}

static inline void ff_poly_mult_7_7 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4,t5, t6;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_6_mults(t5,f[0],f[1],f[2],f[3],f[4],f[5],g[0],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_7_mults(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6],g[0],g[1],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_6_mults(o[7],f[1],f[2],f[3],f[4],f[5],f[6],g[1],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_5_mults(o[8],f[2],f[3],f[4],f[5],f[6],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_4_mults(o[9],f[3],f[4],f[5],f[6],g[3],g[4],g[5],g[6]);
	_ff_sum_3_mults(o[10],f[4],f[5],f[6],g[4],g[5],g[6]);
	_ff_sum_2_mults(o[11],f[5],f[6],g[5],g[6]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[12],f[6],g[6]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6);
}

static inline void ff_poly_mult_7_8 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4,t5, t6, t7;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_6_mults(t5,f[0],f[1],f[2],f[3],f[4],f[5],g[0],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_7_mults(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6],g[0],g[1],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_7_mults(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],g[1],g[2],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_6_mults(o[8],f[1],f[2],f[3],f[4],f[5],f[6],g[2],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_5_mults(o[9],f[2],f[3],f[4],f[5],f[6],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_4_mults(o[10],f[3],f[4],f[5],f[6],g[4],g[5],g[6],g[7]);
	_ff_sum_3_mults(o[11],f[4],f[5],f[6],g[5],g[6],g[7]);
	_ff_sum_2_mults(o[12],f[5],f[6],g[6],g[7]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[13],f[6],g[7]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7);
}

// 27M+59A
static inline void ff_poly_square_8o (ff_t o[], ff_t f[])
{
	ff_t w[7];
	register ff_t t0;
	
	ff_poly_mult_4_4(w,f+4,f);
	ff_poly_square_4(o+8,f+4);
	ff_poly_square_4(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[4],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[5],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[6],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[8],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[9],t0);
	_ff_add(t0,w[6],w[6]); _ff_addto(o[10],t0);
	_ff_add(o[7],w[3],w[3]);
}


static inline void ff_poly_square_8 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_7_mults_s7(o[8],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_6_mults_s6(o[9],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_5_mults_s5(o[10],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_4_mults_s4(o[11],f[4],f[5],f[6],f[7]);
	_ff_sum_3_mults_s3(o[12],f[5],f[6],f[7]);
	_ff_dbl_mult(o[13],f[6],f[7]);
	_ff_square(o[0],f[0]); _ff_square(o[14],f[7]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7);
}

// computes f^2 mod g  for f degree 7 and g(x) = x^8-g6x^6-g5x^5-g4x^4-g3x^3-g2x^2-g1x-g0 (note signs!)
static inline void ff_poly_square_mod_8 (ff_t o[8], ff_t f[8], ff_t g[7])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14;

	_ff_square(t14,f[7]);																	// t14 = f7^2
	_ff_dbl_mult(t13,f[6],f[7]);																// t13 = 2f6f7
	_ff_sum_3_mults_d1 (t12,f[5],f[6],t14,g[6],f[6],f[7]);											// t12 = 2f5f7+f6^2+t14g6
	_ff_sum_4_mults_d2(t11,f[4],f[5],t13,t14,g[5],g[6],f[6],f[7]);										// t11 = 2f4f7+2f5f6+t13g6+t14g5
	_ff_sum_6_mults_d2(t10,f[3],f[4],f[5],t12,t13,t14,g[4],g[5],g[6],f[5],f[6],f[7]);						// t10 = 2f3f7+2f4f6+f5^2+t12g6+t13g5+t14g4
	_ff_sum_7_mults_d3(t9,f[2],f[3],f[4],t11,t12,t13,t14,g[3],g[4],g[5],g[6],f[5],f[6],f[7]);					// t9 = 2f2f7+2f3f6+2f4f5+t11g6+t12g5+t13g4+t14g3
	_ff_sum_9_mults_d3(t8,f[1],f[2],f[3],f[4],t10,t11,t12,t13,t14,g[2],g[3],g[4],g[5],g[6],f[4],f[5],f[6],f[7]);	// t8 = 2f1f7+2f2f6+2f3f5+f4^2+t10g6+t11g5+t12g4+t13g3+g14g2
	
	_ff_sum_10_mults_d4(t7,f[0],f[1],f[2],f[3],t9,t10,t11,t12,t13,t14,g[1],g[2],g[3],g[4],g[5],g[6],f[4],f[5],f[6],f[7]);		// o[7] = 2f0f7+2f1f6+2f2f5+2f3f4+t9g6+t10g5+t12g4+t12g3+t13g2+t14g1
	_ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t8,t9,t10,t11,t12,t13,t14,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);	// o[6] = 2f0f6+2f1f5+2f2f4+f3^2+t8g6+t9g5+t10g4+t11g3+t12g2+t13g1+t14g0
	_ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t8,t9,t10,t11,t12,t13,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);	// o[5] = 2f0f5+2f1f4+2f2f3+t8g5+t9g4+t10g3+t11g2+t12g1+t13g0
	_ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t8,t9,t10,t11,t12,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);			// o[4] = 2f0f4+2f1f3+f2^2+t8g4+t9g3+t10g2+t11g1+t12g0
	_ff_sum_6_mults_d2(t3,f[0],f[1],t8,t9,t10,t11,g[0],g[1],g[2],g[3],f[2],f[3]);						// o[3] = 2f0f3+2f1f2+t8t3+t9g2+t10g1+t11g0
	_ff_sum_5_mults_d1(t2,f[0],f[1],t8,t9,t10,g[0],g[1],g[2],f[1],f[2]);								// o[2] = 2f0f2+f1^2+t8g2+t9g1+t10g0
	_ff_sum_3_mults_d1(t1,f[0],t8,t9,g[0],g[1],f[1]);											// o[1] = 2f0f1+t8g1+t9g0
	_ff_sum_2_mults(o[0],f[0],t8,g[0],f[0]);													// o[0] = f0^2+t8g0
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7);
}

// Computes x*f^2 mod g for f degree 7, g=x^8-g6x^6-g5x^5-g4x^4-g3x^3g2x^2-g1x-g0 (note signs)
static inline void ff_poly_square_mult_x_mod_8 (ff_t o[8], ff_t f[8], ff_t g[7])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15;

	_ff_square(t15,f[7]);																	// t15 = f7^2
	_ff_dbl_mult(t14,f[6],f[7]);																// t14 = 2f6f7
	_ff_sum_3_mults_d1 (t13,f[5],f[6],t15,g[6],f[6],f[7]);											// t13 = 2f5f7+f6^2+t15g6
	_ff_sum_4_mults_d2(t12,f[4],f[5],t14,t15,g[5],g[6],f[6],f[7]);										// t12 = 2f4f7+2f5f6+t14g6+t15g5
	_ff_sum_6_mults_d2(t11,f[3],f[4],f[5],t13,t14,t15,g[4],g[5],g[6],f[5],f[6],f[7]);						// t11 = 2f3f7+2f4f6+f5^2+t13g6+t14g5+t15g4
	_ff_sum_7_mults_d3(t10,f[2],f[3],f[4],t12,t13,t14,t15,g[3],g[4],g[5],g[6],f[5],f[6],f[7]);				// t10 = 2f2f7+2f3f6+2f4f5+t12g6+t13g5+t14g4+t15g3
	_ff_sum_9_mults_d3(t9,f[1],f[2],f[3],f[4],t11,t12,t13,t14,t15,g[2],g[3],g[4],g[5],g[6],f[4],f[5],f[6],f[7]);	// t9 = 2f1f7+2f2f6+2f3f5+f4^2+t11g6+t12g5+t13g4+t14g3+t15g2
	_ff_sum_10_mults_d4(t8,f[0],f[1],f[2],f[3],t10,t11,t12,t13,t14,t15,g[1],g[2],g[3],g[4],g[5],g[6],f[4],f[5],f[6],f[7]);// t8 = 2f0f7+2f1f6+2f2f5+2f3f4+t10g6+t12g5+t12g4+t13g3+t14g2+t15g1
	
	_ff_sum_11_mults_d3(t7,f[0],f[1],f[2],f[3],t9,t10,t11,t12,t13,t14,t15,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);	// o[7] = 2f0f6+2f1f5+2f2f4+f3^2+t9g6+t10g5+t11g4+t12g3+t13g2+t14g1+t15g0
	_ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t8,t9,t10,t11,t12,t13,t14,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);	// o[6] = 2f0f5+2f1f4+2f2f3++t8g6+t9g5+t10g4+t11g3+t12g2+t13g1+t14g0
	_ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t8,t9,t10,t11,t12,t13,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);	// o[5] = 2f0f4+2f1f3+f2^2+t8g5+t9g4+t10g3+t11g2+t12g1+t13g0
	_ff_sum_7_mults_d2(t4,f[0],f[1],t8,t9,t10,t11,t12,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);					// o[4] = 2f0f3+2f1f2+t8g4+t9g3+t10g2+t11g1+t12g0
	_ff_sum_6_mults_d1(t3,f[0],f[1],t8,t9,t10,t11,g[0],g[1],g[2],g[3],f[1],f[2]);							// o[3] = 2f0f2+f1^2+t8g3+t9g2+t10g1+t11g0
	_ff_sum_4_mults_d1(t2,f[0],t8,t9,t10,g[0],g[1],g[2],f[1]);										// o[2] = 2f0f1+t8g2+t9g1+t10g0
	_ff_sum_3_mults(t1,f[0],t8,t9,g[0],g[1],f[0]);													// o[1] = f0^2+t8g1+t9g0
	_ff_mult(o[0],t8,g[0]);																	// o[0] = t8g0
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7);
}


static inline void ff_poly_mult_8_8 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_6_mults(t5,f[0],f[1],f[2],f[3],f[4],f[5],g[0],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_7_mults(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6],g[0],g[1],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_8_mults(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_7_mults(o[8],f[1],f[2],f[3],f[4],f[5],f[6],f[7],g[1],g[2],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_6_mults(o[9],f[2],f[3],f[4],f[5],f[6],f[7],g[2],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_5_mults(o[10],f[3],f[4],f[5],f[6],f[7],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_4_mults(o[11],f[4],f[5],f[6],f[7],g[4],g[5],g[6],g[7]);
	_ff_sum_3_mults(o[12],f[5],f[6],f[7],g[5],g[6],g[7]);
	_ff_sum_2_mults(o[13],f[6],f[7],g[6],g[7]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[14],f[7],g[7]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7);
}

static inline void ff_poly_mult_8_9 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4,t5, t6, t7, t8;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_6_mults(t5,f[0],f[1],f[2],f[3],f[4],f[5],g[0],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_7_mults(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6],g[0],g[1],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_8_mults(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_8_mults(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_7_mults(o[9],f[1],f[2],f[3],f[4],f[5],f[6],f[7],g[2],g[3],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_6_mults(o[10],f[2],f[3],f[4],f[5],f[6],f[7],g[3],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_5_mults(o[11],f[3],f[4],f[5],f[6],f[7],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_4_mults(o[12],f[4],f[5],f[6],f[7],g[5],g[6],g[7],g[8]);
	_ff_sum_3_mults(o[13],f[5],f[6],f[7],g[6],g[7],g[8]);
	_ff_sum_2_mults(o[14],f[6],f[7],g[7],g[8]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[15],f[7],g[8]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8);
}

// 36M+71A
static inline void ff_poly_square_9o (ff_t o[], ff_t f[])
{
	ff_t w[8];
	register ff_t t0;

	ff_poly_mult_4_5(w,f+5,f);
	ff_poly_square_4(o+10,f+5);
	ff_poly_square_5(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[5],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[6],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[7],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[8],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[10],t0);
	_ff_add(t0,w[6],w[6]); _ff_addto(o[11],t0);
	_ff_add(t0,w[7],w[7]); _ff_addto(o[12],t0);
	_ff_add(o[9],w[4],w[4]);
}

static inline void ff_poly_square_9 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_8_mults_s8(o[9],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_7_mults_s7(o[10],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_6_mults_s6(o[11],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_5_mults_s5(o[12],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_4_mults_s4(o[13],f[5],f[6],f[7],f[8]);
	_ff_sum_3_mults_s3(o[14],f[6],f[7],f[8]);
	_ff_dbl_mult(o[15],f[7],f[8]);
	_ff_square(o[0],f[0]); _ff_square(o[16],f[8]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8);
}

static inline void ff_poly_mult_9_9 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_6_mults(t5,f[0],f[1],f[2],f[3],f[4],f[5],g[0],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_7_mults(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6],g[0],g[1],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_8_mults(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_9_mults(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_8_mults(o[9],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_7_mults(o[10],f[2],f[3],f[4],f[5],f[6],f[7],f[8],g[2],g[3],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_6_mults(o[11],f[3],f[4],f[5],f[6],f[7],f[8],g[3],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_5_mults(o[12],f[4],f[5],f[6],f[7],f[8],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_4_mults(o[13],f[5],f[6],f[7],f[8],g[5],g[6],g[7],g[8]);
	_ff_sum_3_mults(o[14],f[6],f[7],f[8],g[6],g[7],g[8]);
	_ff_sum_2_mults(o[15],f[7],f[8],g[7],g[8]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[16],f[8],g[8]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8);
}

static inline void ff_poly_mult_9_10 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4,t5, t6, t7, t8, t9;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_6_mults(t5,f[0],f[1],f[2],f[3],f[4],f[5],g[0],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_7_mults(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6],g[0],g[1],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_8_mults(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_9_mults(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_9_mults(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_8_mults(o[10],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_7_mults(o[11],f[2],f[3],f[4],f[5],f[6],f[7],f[8],g[3],g[4],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_6_mults(o[12],f[3],f[4],f[5],f[6],f[7],f[8],g[4],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_5_mults(o[13],f[4],f[5],f[6],f[7],f[8],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_4_mults(o[14],f[5],f[6],f[7],f[8],g[6],g[7],g[8],g[9]);
	_ff_sum_3_mults(o[15],f[6],f[7],f[8],g[7],g[8],g[9]);
	_ff_sum_2_mults(o[16],f[7],f[8],g[8],g[9]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[17],f[8],g[9]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9);
}

static inline void ff_poly_square_mod_9 (ff_t o[9], ff_t f[9], ff_t g[8])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16;

    _ff_square(t16,f[8]);
    _ff_dbl_mult(t15,f[7],f[8]);
    _ff_sum_3_mults_d1(t14,f[6],f[7],t16,g[7],f[7],f[8]);
    _ff_sum_4_mults_d2(t13,f[5],f[6],t15,t16,g[6],g[7],f[7],f[8]);
    _ff_sum_6_mults_d2(t12,f[4],f[5],f[6],t14,t15,t16,g[5],g[6],g[7],f[6],f[7],f[8]);
    _ff_sum_7_mults_d3(t11,f[3],f[4],f[5],t13,t14,t15,t16,g[4],g[5],g[6],g[7],f[6],f[7],f[8]);
    _ff_sum_9_mults_d3(t10,f[2],f[3],f[4],f[5],t12,t13,t14,t15,t16,g[3],g[4],g[5],g[6],g[7],f[5],f[6],f[7],f[8]);
    _ff_sum_10_mults_d4(t9,f[1],f[2],f[3],f[4],t11,t12,t13,t14,t15,t16,g[2],g[3],g[4],g[5],g[6],g[7],f[5],f[6],f[7],f[8]);

    _ff_sum_12_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t10,t11,t12,t13,t14,t15,t16,g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t9,t10,t11,t12,t13,t14,t15,t16,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t9,t10,t11,t12,t13,t14,t15,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
    _ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t9,t10,t11,t12,t13,t14,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
    _ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t9,t10,t11,t12,t13,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
    _ff_sum_6_mults_d2(t3,f[0],f[1],t9,t10,t11,t12,g[0],g[1],g[2],g[3],f[2],f[3]);
    _ff_sum_5_mults_d1(t2,f[0],f[1],t9,t10,t11,g[0],g[1],g[2],f[1],f[2]);
    _ff_sum_3_mults_d1(t1,f[0],t9,t10,g[0],g[1],f[1]);
    _ff_sum_2_mults(o[0],f[0],t9,g[0],f[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8);
}

static inline void ff_poly_square_mult_x_mod_9 (ff_t o[9], ff_t f[9], ff_t g[8])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17;

    _ff_square(t17,f[8]);
    _ff_dbl_mult(t16,f[7],f[8]);
    _ff_sum_3_mults_d1(t15,f[6],f[7],t17,g[7],f[7],f[8]);
    _ff_sum_4_mults_d2(t14,f[5],f[6],t16,t17,g[6],g[7],f[7],f[8]);
    _ff_sum_6_mults_d2(t13,f[4],f[5],f[6],t15,t16,t17,g[5],g[6],g[7],f[6],f[7],f[8]);
    _ff_sum_7_mults_d3(t12,f[3],f[4],f[5],t14,t15,t16,t17,g[4],g[5],g[6],g[7],f[6],f[7],f[8]);
    _ff_sum_9_mults_d3(t11,f[2],f[3],f[4],f[5],t13,t14,t15,t16,t17,g[3],g[4],g[5],g[6],g[7],f[5],f[6],f[7],f[8]);
    _ff_sum_10_mults_d4(t10,f[1],f[2],f[3],f[4],t12,t13,t14,t15,t16,t17,g[2],g[3],g[4],g[5],g[6],g[7],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t11,t12,t13,t14,t15,t16,t17,g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7],f[8]);

    _ff_sum_12_mults_d4(t8,f[0],f[1],f[2],f[3],t10,t11,t12,t13,t14,t15,t16,t17,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t9,t10,t11,t12,t13,t14,t15,t16,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
    _ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t9,t10,t11,t12,t13,t14,t15,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
    _ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t9,t10,t11,t12,t13,t14,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
    _ff_sum_7_mults_d2(t4,f[0],f[1],t9,t10,t11,t12,t13,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
    _ff_sum_6_mults_d1(t3,f[0],f[1],t9,t10,t11,t12,g[0],g[1],g[2],g[3],f[1],f[2]);
    _ff_sum_4_mults_d1(t2,f[0],t9,t10,t11,g[0],g[1],g[2],f[1]);
    _ff_sum_3_mults(t1,f[0],t9,t10,g[0],g[1],f[0]);
    _ff_mult(o[0],t9,g[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8);
}

// 42M+95A
static inline void ff_poly_square_10o (ff_t o[], ff_t f[])
{
	ff_t w[9];
	register ff_t t0;
	
	ff_poly_mult_5_5(w,f+5,f);
	ff_poly_square_5(o+10,f+5);
	ff_poly_square_5(o,f);
	// it is worth unwinding this loop, saves 6% or 7%
	_ff_add(t0,w[0],w[0]); _ff_addto(o[5],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[6],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[7],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[8],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[10],t0);
	_ff_add(t0,w[6],w[6]); _ff_addto(o[11],t0);
	_ff_add(t0,w[7],w[7]); _ff_addto(o[12],t0);
	_ff_add(t0,w[8],w[8]); _ff_addto(o[13],t0);
	_ff_add(o[9],w[4],w[4]);
}

static inline void ff_poly_square_10 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_10_mults_s10(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_9_mults_s9(o[10],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_8_mults_s8(o[11],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_7_mults_s7(o[12],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_6_mults_s6(o[13],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_5_mults_s5(o[14],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_4_mults_s4(o[15],f[6],f[7],f[8],f[9]);
	_ff_sum_3_mults_s3(o[16],f[7],f[8],f[9]);
	_ff_dbl_mult(o[17],f[8],f[9]);
	_ff_square(o[0],f[0]); _ff_square(o[18],f[9]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9);
}

static inline void ff_poly_mult_10_10 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_6_mults(t5,f[0],f[1],f[2],f[3],f[4],f[5],g[0],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_7_mults(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6],g[0],g[1],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_8_mults(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_9_mults(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_10_mults(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_9_mults(o[10],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_8_mults(o[11],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_7_mults(o[12],f[3],f[4],f[5],f[6],f[7],f[8],f[9],g[3],g[4],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_6_mults(o[13],f[4],f[5],f[6],f[7],f[8],f[9],g[4],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_5_mults(o[14],f[5],f[6],f[7],f[8],f[9],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_4_mults(o[15],f[6],f[7],f[8],f[9],g[6],g[7],g[8],g[9]);
	_ff_sum_3_mults(o[16],f[7],f[8],f[9],g[7],g[8],g[9]);
	_ff_sum_2_mults(o[17],f[8],f[9],g[8],g[9]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[18],f[9],g[9]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9);
}

static inline void ff_poly_mult_10_11 (ff_t o[], ff_t f[], ff_t g[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
	
	_ff_sum_2_mults(t1,f[0],f[1],g[0],g[1]);
	_ff_sum_3_mults(t2,f[0],f[1],f[2],g[0],g[1],g[2]);
	_ff_sum_4_mults(t3,f[0],f[1],f[2],f[3],g[0],g[1],g[2],g[3]);
	_ff_sum_5_mults(t4,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3],g[4]);
	_ff_sum_6_mults(t5,f[0],f[1],f[2],f[3],f[4],f[5],g[0],g[1],g[2],g[3],g[4],g[5]);
	_ff_sum_7_mults(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6],g[0],g[1],g[2],g[3],g[4],g[5],g[6]);
	_ff_sum_8_mults(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7]);
	_ff_sum_9_mults(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8]);
	_ff_sum_10_mults(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9]);
	_ff_sum_10_mults(t10,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10]);
	_ff_sum_9_mults(o[11],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10]);
	_ff_sum_8_mults(o[12],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10]);
	_ff_sum_7_mults(o[13],f[3],f[4],f[5],f[6],f[7],f[8],f[9],g[4],g[5],g[6],g[7],g[8],g[9],g[10]);
	_ff_sum_6_mults(o[14],f[4],f[5],f[6],f[7],f[8],f[9],g[5],g[6],g[7],g[8],g[9],g[10]);
	_ff_sum_5_mults(o[15],f[5],f[6],f[7],f[8],f[9],g[6],g[7],g[8],g[9],g[10]);
	_ff_sum_4_mults(o[16],f[6],f[7],f[8],f[9],g[7],g[8],g[9],g[10]);
	_ff_sum_3_mults(o[17],f[7],f[8],f[9],g[8],g[9],g[10]);
	_ff_sum_2_mults(o[18],f[8],f[9],g[9],g[10]);
	_ff_mult(o[0],f[0],g[0]); _ff_mult(o[19],f[9],g[10]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9);  _ff_set(o[10],t10);
}

// computes f^2 mod g  for f degree 9 and g(x) = x^10-g8x^8-g6x^6-g5x^5-g4x^4-g3x^3-g2x^2-g1x-g0 (note signs!)
static inline void ff_poly_square_mod_10 (ff_t o[10], ff_t f[10], ff_t g[9])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18;

	_ff_square(t18,f[9]);																	// t18 = f9^2
	_ff_dbl_mult(t17,f[8],f[9]);																// t17 = 2f8f9
	_ff_sum_3_mults_d1 (t16,f[7],f[8],t18,g[8],f[8],f[9]);											// t16 = 2f7f9+f8^2+t18g8
	_ff_sum_4_mults_d2(t15,f[6],f[7],t17,t18,g[7],g[8],f[8],f[9]);										// t15 = 2f6f9+2f7f8+t17g8+t18g7
	_ff_sum_6_mults_d2(t14,f[5],f[6],f[7],t16,t17,t18,g[6],g[7],g[8],f[7],f[8],f[9]);						// t14 = 2f5f9+2f6f8+f7^2+t16g6+t17g7+t18g6
	_ff_sum_7_mults_d3(t13,f[4],f[5],f[6],t15,t16,t17,t18,g[5],g[6],g[7],g[8],f[7],f[8],f[9]);				// t13 = 2f4f9+2f5f8+2f6f7+t15g8+t16g7+t17g6+t18g5
	_ff_sum_9_mults_d3(t12,f[3],f[4],f[5],f[6],t14,t15,t16,t17,t18,g[4],g[5],g[6],g[7],g[8],f[6],f[7],f[8],f[9]);	// t12 = 2f3f9+2f4f8+2f5f7+f6^2+t14g8+t15g7+t16g6+t17g5+t18g4
	_ff_sum_10_mults_d4(t11,f[2],f[3],f[4],f[5],t13,t14,t15,t16,t17,t18,g[3],g[4],g[5],g[6],g[7],g[8],f[6],f[7],f[8],f[9]);	// t11 = 2f2f9+2f3f8+2f4f7+2f5f6+t13g8+t14g7+t15g6+t16g5+t17g4+t18g3
	_ff_sum_12_mults_d4(t10,f[1],f[2],f[3],f[4],f[5],t12,t13,t14,t15,t16,t17,t18,g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[5],f[6],f[7],f[8],f[9]);	// t10 = 2f1f9+2f2f8+2f3f7+2f4f6+f5^2+t12g8+t13g7+t14g6+t15g5+t16g4+t17g3+t18g2
	
	_ff_sum_13_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t11,t12,t13,t14,t15,t16,t17,t18,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[5],f[6],f[7],f[8],f[9]);	// o[9] = 2f0f9+2f1f8+2f2f7+2f3f6+2f4f5+t11g8+t12g7+t13g6+t14g5+t15g4+t16g3+t17g2+t18g1
	_ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t10,t11,t12,t13,t14,t15,t16,t17,t18,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);	// o[8] = 2f0f8+2f1f7+2f2f6+2f3f5+f4^2+t10g8+t11g7+t12g6+t13g5+t14g4+t15g3+t16g2+t17g1+t18g0
	_ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t10,t11,t12,t13,t14,t15,t16,t17,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);	// o[7] = 2f0f7+2f1f6+2f2f5+2f3f4+t10g7+t11g6+t12g5+t13g4+t14g3+t15g2+t16g1+t17g0
	_ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t10,t11,t12,t13,t14,t15,t16,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);	// o[6] = 2f0f6+2f1f5+2f2f4+f3^2+t10g6+t11g5+t12g4+t13g3+t14g2+t15g1+t16g0
	_ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t10,t11,t12,t13,t14,t15,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);	// o[5] = 2f0f5+2f1f4+2f2f3+t10g5+t11g4+t12g3+t13g2+t14g1+t15g0
	_ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t10,t11,t12,t13,t14,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);			// o[4] = 2f0f4+2f1f3+f2^2+t10g4+t11g3+t12g2+t13g1+t14g0
	_ff_sum_6_mults_d2(t3,f[0],f[1],t10,t11,t12,t13,g[0],g[1],g[2],g[3],f[2],f[3]);						// o[3] = 2f0f3+2f1f2+t10t3+t11g2+t12g1+t13g0
	_ff_sum_5_mults_d1(t2,f[0],f[1],t10,t11,t12,g[0],g[1],g[2],f[1],f[2]);								// o[2] = 2f0f2+f1^2+t10g2+t11g1+t12g0
	_ff_sum_3_mults_d1(t1,f[0],t10,t11,g[0],g[1],f[1]);											// o[1] = 2f0f1+t10g1+t11g0
	_ff_sum_2_mults(o[0],f[0],t10,g[0],f[0]);														// o[0] = f0^2+t10g0
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9);
}

// Computes x*f^2 mod g for f degree 9, g=x^10-g8x^8-g6x^6-g5x^5-g4x^4-g3x^3g2x^2-g1x-g0 (note signs)
static inline void ff_poly_square_mult_x_mod_10 (ff_t o[10], ff_t f[10], ff_t g[9])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19;

	_ff_square(t19,f[9]);																	// t19 = f9^2
	_ff_dbl_mult(t18,f[8],f[9]);																// t18 = 2f8f9
	_ff_sum_3_mults_d1 (t17,f[7],f[8],t19,g[8],f[8],f[9]);											// t17 = 2f7f9+f8^2+t19g8
	_ff_sum_4_mults_d2(t16,f[6],f[7],t18,t19,g[7],g[8],f[8],f[9]);										// t16 = 2f6f9+2f7f8+t18g8+t19g7
	_ff_sum_6_mults_d2(t15,f[5],f[6],f[7],t17,t18,t19,g[6],g[7],g[8],f[7],f[8],f[9]);						// t15 = 2f5f9+2f6f8+f7^2+t17g8+t18g7+t19g6
	_ff_sum_7_mults_d3(t14,f[4],f[5],f[6],t16,t17,t18,t19,g[5],g[6],g[7],g[8],f[7],f[8],f[9]);				// t14 = 2f4f9+2f5f8+2f4f7+t16g8+t17g7+t18g6+t19g5
	_ff_sum_9_mults_d3(t13,f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,g[4],g[5],g[6],g[7],g[8],f[6],f[7],f[8],f[9]);	// t13 = 2f3f9+2f4f8+2f3f7+f6^2+t15g8+t16g7+t17g6+t18g5+t19g4
	_ff_sum_10_mults_d4(t12,f[2],f[3],f[4],f[5],t14,t15,t16,t17,t18,t19,g[3],g[4],g[5],g[6],g[7],g[8],f[6],f[7],f[8],f[9]);// t12 = 2f2f9+2f3f8+2f4f7+2f5f6+t14g8+t15g7+t16g6+t17g5+t18g4+t19g3
	_ff_sum_12_mults_d4(t11,f[1],f[2],f[3],f[4],f[5],t13,t14,t15,t16,t17,t18,t19,g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[5],f[6],f[7],f[8],f[9]);// t11 = 2f1f9+2f2f8+2f3f7+2f4f6+f5^2+t13g8+t14g7+t15g6+t16g5+t17g4+t18g3+t19g2
	_ff_sum_13_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t12,t13,t14,t15,t16,t17,t18,t19,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[5],f[6],f[7],f[8],f[9]);// t10 = 2f0f9+2f1f8+2f2f7+2f3f6+2f4f5+t12g8+t13g7+t14g6+t15g5+t16g4+t17g3+t18g2+t19g1
	
	_ff_sum_14_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t11,t12,t13,t14,t15,t16,t17,t18,t19,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);	// o[7] = 2f0f8+2f1f7+2f2f6+2f3f5+f4^2+t11g8+t12g7+t13g6+t14g5+t15g4+t16g3+t17g2+t18g1+t19g0
	_ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t10,t11,t12,t13,t14,t15,t16,t17,t18,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);	// o[7] = 2f0f7+2f1f6+2f2f5+2f3f4+t10g8+t11g7+t12g6+t13g5+t14g4+t15g3+t16g2+t17g1+t18g0
	_ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t10,t11,t12,t13,t14,t15,t16,t17,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);	// o[7] = 2f0f6+2f1f5+2f2f4+f3^2+t10g7+t11g6+t12g5+t13g4+t14g3+t15g2+t16g1+t17g0
	_ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t10,t11,t12,t13,t14,t15,t16,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);	// o[6] = 2f0f5+2f1f4+2f2f3++t10g6+t11g5+t12g4+t13g3+t14g2+t15g1+t16g0
	_ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t10,t11,t12,t13,t14,t15,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);	// o[5] = 2f0f4+2f1f3+f2^2+t10g5+t11g4+t12g3+t13g2+t14g1+t15g0
	_ff_sum_7_mults_d2(t4,f[0],f[1],t10,t11,t12,t13,t14,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);				// o[4] = 2f0f3+2f1f2+t10g4+t11g3+t12g2+t13g1+t14g0
	_ff_sum_6_mults_d1(t3,f[0],f[1],t10,t11,t12,t13,g[0],g[1],g[2],g[3],f[1],f[2]);						// o[3] = 2f0f2+f1^2+t10g3+t11g2+t12g1+t13g0
	_ff_sum_4_mults_d1(t2,f[0],t10,t11,t12,g[0],g[1],g[2],f[1]);										// o[2] = 2f0f1+t10g2+t11g1+t12g0
	_ff_sum_3_mults(t1,f[0],t10,t11,g[0],g[1],f[0]);												// o[1] = f0^2+t10g1+t11g0
	_ff_mult(o[0],t10,g[0]);																	// o[0] = t10g0
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9);
}

//49M+115A
static inline void ff_poly_square_11o (ff_t o[], ff_t f[])
{
	ff_t w[10];
	register ff_t t0;

	ff_poly_mult_5_6(w,f+6,f);
	ff_poly_square_5(o+12,f+6);
	ff_poly_square_6(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[6],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[7],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[8],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[9],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[10],t0);
	_ff_add(t0,w[6],w[6]); _ff_addto(o[12],t0);
	_ff_add(t0,w[7],w[7]); _ff_addto(o[13],t0);
	_ff_add(t0,w[8],w[8]); _ff_addto(o[14],t0);
	_ff_add(t0,w[9],w[9]); _ff_addto(o[15],t0);
	_ff_add(o[11],w[5],w[5]);
}

static inline void ff_poly_square_11 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_10_mults_s10(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_11_mults_s11(t10,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_10_mults_s10(o[11],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_9_mults_s9(o[12],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_8_mults_s8(o[13],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_7_mults_s7(o[14],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_6_mults_s6(o[15],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_5_mults_s5(o[16],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_4_mults_s4(o[17],f[7],f[8],f[9],f[10]);
	_ff_sum_3_mults_s3(o[18],f[8],f[9],f[10]);
	_ff_dbl_mult(o[19],f[9],f[10]);
	_ff_square(o[0],f[0]); _ff_square(o[20],f[10]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10);
}

// computes f^2 mod g  for f degree 10 and g(x) = x^11-g9x^9-g8x^8-...-g1x-g0 (note signs!)
static inline void ff_poly_square_mod_11 (ff_t o[11], ff_t f[11], ff_t g[10])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20;

	_ff_square(t20,f[10]);																	
	_ff_dbl_mult(t19,f[9],f[10]);																
	_ff_sum_3_mults_d1 (t18,f[8],f[9],t20,g[9],f[9],f[10]);										
	_ff_sum_4_mults_d2(t17,f[7],f[8],t19,t20,g[8],g[9],f[9],f[10]);										
	_ff_sum_6_mults_d2(t16,f[6],f[7],f[8],t18,t19,t20,g[7],g[8],g[9],f[8],f[9],f[10]);					
	_ff_sum_7_mults_d3(t15,f[5],f[6],f[7],t17,t18,t19,t20,g[6],g[7],g[8],g[9],f[8],f[9],f[10]);					
	_ff_sum_9_mults_d3(t14,f[4],f[5],f[6],f[7],t16,t17,t18,t19,t20,g[5],g[6],g[7],g[8],g[9],f[7],f[8],f[9],f[10]);	
	_ff_sum_10_mults_d4(t13,f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,t20,g[4],g[5],g[6],g[7],g[8],g[9],f[7],f[8],f[9],f[10]);
	_ff_sum_12_mults_d4(t12,f[2],f[3],f[4],f[5],f[6],t14,t15,t16,t17,t18,t19,t20,g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[6],f[7],f[8],f[9],f[10]);	
	_ff_sum_13_mults_d5(t11,f[1],f[2],f[3],f[4],f[5],t13,t14,t15,t16,t17,t18,t19,t20,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[6],f[7],f[8],f[9],f[10]);	
	
	_ff_sum_15_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t12,t13,t14,t15,t16,t17,t18,t19,t20,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t11,t12,t13,t14,t15,t16,t17,t18,t19,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t11,t12,t13,t14,t15,t16,t17,t18,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
	_ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t11,t12,t13,t14,t15,t16,t17,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
	_ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t11,t12,t13,t14,t15,t16,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);	
	_ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t11,t12,t13,t14,t15,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);		
	_ff_sum_6_mults_d2(t3,f[0],f[1],t11,t12,t13,t14,g[0],g[1],g[2],g[3],f[2],f[3]);						
	_ff_sum_5_mults_d1(t2,f[0],f[1],t11,t12,t13,g[0],g[1],g[2],f[1],f[2]);							
	_ff_sum_3_mults_d1(t1,f[0],t11,t12,g[0],g[1],f[1]);	
	_ff_sum_2_mults(o[0],f[0],t11,g[0],f[0]);												
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10);
}

// Computes x*f^2 mod g for f degree 10, g(x) = x^11-g9x^9-g8x^8-...-g1x-g0 (note signs!)
static inline void ff_poly_square_mult_x_mod_11 (ff_t o[11], ff_t f[11], ff_t g[10])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21;

	_ff_square(t21,f[10]);																
	_ff_dbl_mult(t20,f[9],f[10]);																
	_ff_sum_3_mults_d1 (t19,f[8],f[9],t21,g[9],f[9],f[10]);											
	_ff_sum_4_mults_d2(t18,f[7],f[8],t20,t21,g[8],g[9],f[9],f[10]);										
	_ff_sum_6_mults_d2(t17,f[6],f[7],f[8],t19,t20,t21,g[7],g[8],g[9],f[8],f[9],f[10]);					
	_ff_sum_7_mults_d3(t16,f[5],f[6],f[7],t18,t19,t20,t21,g[6],g[7],g[8],g[9],f[8],f[9],f[10]);				
	_ff_sum_9_mults_d3(t15,f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,g[5],g[6],g[7],g[8],g[9],f[7],f[8],f[9],f[10]);	
	_ff_sum_10_mults_d4(t14,f[3],f[4],f[5],f[6],t16,t17,t18,t19,t20,t21,g[4],g[5],g[6],g[7],g[8],g[9],f[7],f[8],f[9],f[10]);
	_ff_sum_12_mults_d4(t13,f[2],f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,t20,t21,g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_13_mults_d5(t12,f[1],f[2],f[3],f[4],f[5],t14,t15,t16,t17,t18,t19,t20,t21,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_15_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t13,t14,t15,t16,t17,t18,t19,t20,t21,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9],f[10]);
	
	_ff_sum_15_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t11,t12,t13,t14,t15,t16,t17,t18,t19,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
	_ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t11,t12,t13,t14,t15,t16,t17,t18,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);	
	_ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t11,t12,t13,t14,t15,t16,t17,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
	_ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t11,t12,t13,t14,t15,t16,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
	_ff_sum_7_mults_d2(t4,f[0],f[1],t11,t12,t13,t14,t15,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);					
	_ff_sum_6_mults_d1(t3,f[0],f[1],t11,t12,t13,t14,g[0],g[1],g[2],g[3],f[1],f[2]);							
	_ff_sum_4_mults_d1(t2,f[0],t11,t12,t13,g[0],g[1],g[2],f[1]);										
	_ff_sum_3_mults(t1,f[0],t11,t12,g[0],g[1],f[0]);													
	_ff_mult(o[0],t11,g[0]);																
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10);
}


// 54M+140A 
static inline void ff_poly_square_12o (ff_t o[], ff_t f[])
{
	ff_t w[11];
	register ff_t t0;
	
	ff_poly_mult_6_6(w,f+6,f);
	ff_poly_square_6(o+12,f+6);
	ff_poly_square_6(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[6],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[7],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[8],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[9],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[10],t0);
	_ff_add(t0,w[6],w[6]); _ff_addto(o[12],t0);
	_ff_add(t0,w[7],w[7]); _ff_addto(o[13],t0);
	_ff_add(t0,w[8],w[8]); _ff_addto(o[14],t0);
	_ff_add(t0,w[9],w[9]); _ff_addto(o[15],t0);
	_ff_add(t0,w[10],w[10]); _ff_addto(o[16],t0);
	_ff_add(o[11],w[5],w[5]);
}

static inline void ff_poly_square_12 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_10_mults_s10(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_11_mults_s11(t10,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_12_mults_s12(t11,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_11_mults_s11(o[12],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_10_mults_s10(o[13],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_9_mults_s9(o[14],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_8_mults_s8(o[15],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_7_mults_s7(o[16],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_6_mults_s6(o[17],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_5_mults_s5(o[18],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_4_mults_s4(o[19],f[8],f[9],f[10],f[11]);
	_ff_sum_3_mults_s3(o[20],f[9],f[10],f[11]);
	_ff_dbl_mult(o[21],f[10],f[11]);
	_ff_square(o[0],f[0]); _ff_square(o[22],f[11]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11);
}

static inline void ff_poly_square_mod_12 (ff_t o[12], ff_t f[12], ff_t g[11])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22;

	_ff_square(t22,f[11]);
	_ff_dbl_mult(t21,f[10],f[11]);
	_ff_sum_3_mults_d1(t20,f[9],f[10],t22,g[10],f[10],f[11]);
	_ff_sum_4_mults_d2(t19,f[8],f[9],t21,t22,g[9],g[10],f[10],f[11]);
	_ff_sum_6_mults_d2(t18,f[7],f[8],f[9],t20,t21,t22,g[8],g[9],g[10],f[9],f[10],f[11]);
	_ff_sum_7_mults_d3(t17,f[6],f[7],f[8],t19,t20,t21,t22,g[7],g[8],g[9],g[10],f[9],f[10],f[11]);			
	_ff_sum_9_mults_d3(t16,f[5],f[6],f[7],f[8],t18,t19,t20,t21,t22,g[6],g[7],g[8],g[9],g[10],f[8],f[9],f[10],f[11]);
	_ff_sum_10_mults_d4(t15,f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,g[5],g[6],g[7],g[8],g[9],g[10],f[8],f[9],f[10],f[11]);
	_ff_sum_12_mults_d4(t14,f[3],f[4],f[5],f[6],f[7],t16,t17,t18,t19,t20,t21,t22,g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[7],f[8],f[9],f[10],f[11]);	
	_ff_sum_13_mults_d5(t13,f[2],f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,t20,t21,t22,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_15_mults_d5(t12,f[1],f[2],f[3],f[4],f[5],f[6],t14,t15,t16,t17,t18,t19,t20,t21,t22,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[6],f[7],f[8],f[9],f[10],f[11]);
	
	_ff_sum_16_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t12,t13,t14,t15,t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t12,t13,t14,t15,t16,t17,t18,t19,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
	_ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t12,t13,t14,t15,t16,t17,t18,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
	_ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t12,t13,t14,t15,t16,t17,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
	_ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t12,t13,t14,t15,t16,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
	_ff_sum_6_mults_d2(t3,f[0],f[1],t12,t13,t14,t15,g[0],g[1],g[2],g[3],f[2],f[3]);
	_ff_sum_5_mults_d1(t2,f[0],f[1],t12,t13,t14,g[0],g[1],g[2],f[1],f[2]);
	_ff_sum_3_mults_d1(t1,f[0],t12,t13,g[0],g[1],f[1]);
	_ff_sum_2_mults(o[0],f[0],t12,g[0],f[0]);
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11);
}

static inline void ff_poly_square_mult_x_mod_12 (ff_t o[12], ff_t f[12], ff_t g[11])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23;

	_ff_square(t23,f[11]);
	_ff_dbl_mult(t22,f[10],f[11]);
	_ff_sum_3_mults_d1 (t21,f[9],f[10],t23,g[10],f[10],f[11]);
	_ff_sum_4_mults_d2(t20,f[8],f[9],t22,t23,g[9],g[10],f[10],f[11]);
	_ff_sum_6_mults_d2(t19,f[7],f[8],f[9],t21,t22,t23,g[8],g[9],g[10],f[9],f[10],f[11]);
	_ff_sum_7_mults_d3(t18,f[6],f[7],f[8],t20,t21,t22,t23,g[7],g[8],g[9],g[10],f[9],f[10],f[11]);
	_ff_sum_9_mults_d3(t17,f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,g[6],g[7],g[8],g[9],g[10],f[8],f[9],f[10],f[11]);
	_ff_sum_10_mults_d4(t16,f[4],f[5],f[6],f[7],t18,t19,t20,t21,t22,t23,g[5],g[6],g[7],g[8],g[9],g[10],f[8],f[9],f[10],f[11]);
	_ff_sum_12_mults_d4(t15,f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_13_mults_d5(t14,f[2],f[3],f[4],f[5],f[6],t16,t17,t18,t19,t20,t21,t22,t23,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_15_mults_d5(t13,f[1],f[2],f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,t20,t21,t22,t23,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_16_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[6],f[7],f[8],f[9],f[10],f[11]);
	
	_ff_sum_17_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t12,t13,t14,t15,t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
	_ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t12,t13,t14,t15,t16,t17,t18,t19,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
	_ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t12,t13,t14,t15,t16,t17,t18,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
	_ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t12,t13,t14,t15,t16,t17,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
	_ff_sum_7_mults_d2(t4,f[0],f[1],t12,t13,t14,t15,t16,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
	_ff_sum_6_mults_d1(t3,f[0],f[1],t12,t13,t14,t15,g[0],g[1],g[2],g[3],f[1],f[2]);
	_ff_sum_4_mults_d1(t2,f[0],t12,t13,t14,g[0],g[1],g[2],f[1]);
	_ff_sum_3_mults(t1,f[0],t12,t13,g[0],g[1],f[0]);
	_ff_mult(o[0],t12,g[0]);
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11);
}


// 63M+165A
static inline void ff_poly_square_13o (ff_t o[], ff_t f[])
{
	ff_t w[12];
	register ff_t t0;
	
	ff_poly_mult_6_7(w,f+7,f);
	ff_poly_square_6(o+14,f+7);
	ff_poly_square_7(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[7],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[8],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[9],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[10],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[11],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[12],t0);
	_ff_add(t0,w[7],w[7]); _ff_addto(o[14],t0);
	_ff_add(t0,w[8],w[8]); _ff_addto(o[15],t0);
	_ff_add(t0,w[9],w[9]); _ff_addto(o[16],t0);
	_ff_add(t0,w[10],w[10]); _ff_addto(o[17],t0);
	_ff_add(t0,w[11],w[11]); _ff_addto(o[18],t0);
	_ff_add(o[13],w[6],w[6]);
}

static inline void ff_poly_square_13 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_10_mults_s10(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_11_mults_s11(t10,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_12_mults_s12(t11,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_13_mults_s13(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_12_mults_s12(o[13],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_11_mults_s11(o[14],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_10_mults_s10(o[15],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_9_mults_s9(o[16],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_8_mults_s8(o[17],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_7_mults_s7(o[18],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_6_mults_s6(o[19],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_5_mults_s5(o[20],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_4_mults_s4(o[21],f[9],f[10],f[11],f[12]);
	_ff_sum_3_mults_s3(o[22],f[10],f[11],f[12]);
	_ff_dbl_mult(o[23],f[11],f[12]);
	_ff_square(o[0],f[0]); _ff_square(o[24],f[12]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12);
}

// computes f^2 mod g  for f degree 12 and g(x) = x^13-g11x^11-g10x^10-...-g1x-g0 (note signs!)
static inline void ff_poly_square_mod_13 (ff_t o[13], ff_t f[13], ff_t g[12])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24;
	
	_ff_square(t24,f[12]);
	_ff_dbl_mult(t23,f[11],f[12]);																
	_ff_sum_3_mults_d1 (t22,f[10],f[11],t24,g[11],f[11],f[12]);	
	_ff_sum_4_mults_d2(t21,f[9],f[10],t23,t24,g[10],g[11],f[11],f[12]);	
	_ff_sum_6_mults_d2(t20,f[8],f[9],f[10],t22,t23,t24,g[9],g[10],g[11],f[10],f[11],f[12]);		
	_ff_sum_7_mults_d3(t19,f[7],f[8],f[9],t21,t22,t23,t24,g[8],g[9],g[10],g[11],f[10],f[11],f[12]);					
	_ff_sum_9_mults_d3(t18,f[6],f[7],f[8],f[9],t20,t21,t22,t23,t24,g[7],g[8],g[9],g[10],g[11],f[9],f[10],f[11],f[12]);	
	_ff_sum_10_mults_d4(t17,f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,g[6],g[7],g[8],g[9],g[10],g[11],f[9],f[10],f[11],f[12]);
	_ff_sum_12_mults_d4(t16,f[4],f[5],f[6],f[7],f[8],t18,t19,t20,t21,t22,t23,t24,g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[8],f[9],f[10],f[11],f[12]);	
	_ff_sum_13_mults_d5(t15,f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,t24,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_15_mults_d5(t14,f[2],f[3],f[4],f[5],f[6],f[7],t16,t17,t18,t19,t20,t21,t22,t23,t24,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_16_mults_d6(t13,f[1],f[2],f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[7],f[8],f[9],f[10],f[11],f[12]);
	
	_ff_sum_18_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t13,t14,t15,t16,t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t13,t14,t15,t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
	_ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t13,t14,t15,t16,t17,t18,t19,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
	_ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t13,t14,t15,t16,t17,t18,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);	
	_ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t13,t14,t15,t16,t17,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);		
	_ff_sum_6_mults_d2(t3,f[0],f[1],t13,t14,t15,t16,g[0],g[1],g[2],g[3],f[2],f[3]);						
	_ff_sum_5_mults_d1(t2,f[0],f[1],t13,t14,t15,g[0],g[1],g[2],f[1],f[2]);							
	_ff_sum_3_mults_d1(t1,f[0],t13,t14,g[0],g[1],f[1]);	
	_ff_sum_2_mults(o[0],f[0],t13,g[0],f[0]);												
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12);
}

// Computes x*f^2 mod g for f degree 12, g(x) = x^13-g11x^11-g10x^10-...-g1x-g0 (note signs!)
static inline void ff_poly_square_mult_x_mod_13 (ff_t o[13], ff_t f[13], ff_t g[12])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25;		// we could allocate these as an array, performance is essentially the same

	_ff_square(t25,f[12]);													
	_ff_dbl_mult(t24,f[11],f[12]);																
	_ff_sum_3_mults_d1 (t23,f[10],f[11],t25,g[11],f[11],f[12]);		
	_ff_sum_4_mults_d2(t22,f[9],f[10],t24,t25,g[10],g[11],f[11],f[12]);										
	_ff_sum_6_mults_d2(t21,f[8],f[9],f[10],t23,t24,t25,g[9],g[10],g[11],f[10],f[11],f[12]);					
	_ff_sum_7_mults_d3(t20,f[7],f[8],f[9],t22,t23,t24,t25,g[8],g[9],g[10],g[11],f[10],f[11],f[12]);
	_ff_sum_9_mults_d3(t19,f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,g[7],g[8],g[9],g[10],g[11],f[9],f[10],f[11],f[12]);
	_ff_sum_10_mults_d4(t18,f[5],f[6],f[7],f[8],t20,t21,t22,t23,t24,t25,g[6],g[7],g[8],g[9],g[10],g[11],f[9],f[10],f[11],f[12]);
	_ff_sum_12_mults_d4(t17,f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_13_mults_d5(t16,f[3],f[4],f[5],f[6],f[7],t18,t19,t20,t21,t22,t23,t24,t25,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_15_mults_d5(t15,f[2],f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,t24,t25,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_16_mults_d6(t14,f[1],f[2],f[3],f[4],f[5],f[6],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_18_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	
	_ff_sum_18_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t13,t14,t15,t16,t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
	_ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t13,t14,t15,t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);	
	_ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t13,t14,t15,t16,t17,t18,t19,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
	_ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t13,t14,t15,t16,t17,t18,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
	_ff_sum_7_mults_d2(t4,f[0],f[1],t13,t14,t15,t16,t17,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);					
	_ff_sum_6_mults_d1(t3,f[0],f[1],t13,t14,t15,t16,g[0],g[1],g[2],g[3],f[1],f[2]);							
	_ff_sum_4_mults_d1(t2,f[0],t13,t14,t15,g[0],g[1],g[2],f[1]);										
	_ff_sum_3_mults(t1,f[0],t13,t14,g[0],g[1],f[0]);													
	_ff_mult(o[0],t13,g[0]);																
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12);
}


// 70M+198A
static inline void ff_poly_square_14o (ff_t o[], ff_t f[])
{
	ff_t w[13];
	register ff_t t0;
	
	ff_poly_mult_7_7(w,f+7,f);
	ff_poly_square_7(o+14,f+7);
	ff_poly_square_7(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[7],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[8],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[9],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[10],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[11],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[12],t0);
	_ff_add(t0,w[7],w[7]); _ff_addto(o[14],t0);
	_ff_add(t0,w[8],w[8]); _ff_addto(o[15],t0);
	_ff_add(t0,w[9],w[9]); _ff_addto(o[16],t0);
	_ff_add(t0,w[10],w[10]); _ff_addto(o[17],t0);
	_ff_add(t0,w[11],w[11]); _ff_addto(o[18],t0);
	_ff_add(t0,w[12],w[12]); _ff_addto(o[19],t0);
	_ff_add(o[13],w[6],w[6]);
}

static inline void ff_poly_square_14 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_10_mults_s10(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_11_mults_s11(t10,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_12_mults_s12(t11,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_13_mults_s13(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_14_mults_s14(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_13_mults_s13(o[14],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_12_mults_s12(o[15],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_11_mults_s11(o[16],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_10_mults_s10(o[17],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_9_mults_s9(o[18],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_8_mults_s8(o[19],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_7_mults_s7(o[20],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_6_mults_s6(o[21],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_5_mults_s5(o[22],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_4_mults_s4(o[23],f[10],f[11],f[12],f[13]);
	_ff_sum_3_mults_s3(o[24],f[11],f[12],f[13]);
	_ff_dbl_mult(o[25],f[12],f[13]);
	_ff_square(o[0],f[0]); _ff_square(o[26],f[13]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13);
}

static inline void ff_poly_square_mod_14 (ff_t o[14], ff_t f[14], ff_t g[13])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26;

    _ff_square(t26,f[13]);
    _ff_dbl_mult(t25,f[12],f[13]);
    _ff_sum_3_mults_d1(t24,f[11],f[12],t26,g[12],f[12],f[13]);
    _ff_sum_4_mults_d2(t23,f[10],f[11],t25,t26,g[11],g[12],f[12],f[13]);
    _ff_sum_6_mults_d2(t22,f[9],f[10],f[11],t24,t25,t26,g[10],g[11],g[12],f[11],f[12],f[13]);
    _ff_sum_7_mults_d3(t21,f[8],f[9],f[10],t23,t24,t25,t26,g[9],g[10],g[11],g[12],f[11],f[12],f[13]);
    _ff_sum_9_mults_d3(t20,f[7],f[8],f[9],f[10],t22,t23,t24,t25,t26,g[8],g[9],g[10],g[11],g[12],f[10],f[11],f[12],f[13]);
    _ff_sum_10_mults_d4(t19,f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,g[7],g[8],g[9],g[10],g[11],g[12],f[10],f[11],f[12],f[13]);
    _ff_sum_12_mults_d4(t18,f[5],f[6],f[7],f[8],f[9],t20,t21,t22,t23,t24,t25,t26,g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_13_mults_d5(t17,f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_15_mults_d5(t16,f[3],f[4],f[5],f[6],f[7],f[8],t18,t19,t20,t21,t22,t23,t24,t25,t26,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_16_mults_d6(t15,f[2],f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_18_mults_d6(t14,f[1],f[2],f[3],f[4],f[5],f[6],f[7],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);

    _ff_sum_19_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t14,t15,t16,t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t14,t15,t16,t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t14,t15,t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
    _ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t14,t15,t16,t17,t18,t19,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
    _ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t14,t15,t16,t17,t18,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
    _ff_sum_6_mults_d2(t3,f[0],f[1],t14,t15,t16,t17,g[0],g[1],g[2],g[3],f[2],f[3]);
    _ff_sum_5_mults_d1(t2,f[0],f[1],t14,t15,t16,g[0],g[1],g[2],f[1],f[2]);
    _ff_sum_3_mults_d1(t1,f[0],t14,t15,g[0],g[1],f[1]);
    _ff_sum_2_mults(o[0],f[0],t14,g[0],f[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13);
}

static inline void ff_poly_square_mult_x_mod_14 (ff_t o[14], ff_t f[14], ff_t g[13])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27;

    _ff_square(t27,f[13]);
    _ff_dbl_mult(t26,f[12],f[13]);
    _ff_sum_3_mults_d1(t25,f[11],f[12],t27,g[12],f[12],f[13]);
    _ff_sum_4_mults_d2(t24,f[10],f[11],t26,t27,g[11],g[12],f[12],f[13]);
    _ff_sum_6_mults_d2(t23,f[9],f[10],f[11],t25,t26,t27,g[10],g[11],g[12],f[11],f[12],f[13]);
    _ff_sum_7_mults_d3(t22,f[8],f[9],f[10],t24,t25,t26,t27,g[9],g[10],g[11],g[12],f[11],f[12],f[13]);
    _ff_sum_9_mults_d3(t21,f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,g[8],g[9],g[10],g[11],g[12],f[10],f[11],f[12],f[13]);
    _ff_sum_10_mults_d4(t20,f[6],f[7],f[8],f[9],t22,t23,t24,t25,t26,t27,g[7],g[8],g[9],g[10],g[11],g[12],f[10],f[11],f[12],f[13]);
    _ff_sum_12_mults_d4(t19,f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_13_mults_d5(t18,f[4],f[5],f[6],f[7],f[8],t20,t21,t22,t23,t24,t25,t26,t27,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_15_mults_d5(t17,f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_16_mults_d6(t16,f[2],f[3],f[4],f[5],f[6],f[7],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_18_mults_d6(t15,f[1],f[2],f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_19_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);

    _ff_sum_20_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t14,t15,t16,t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
    _ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t14,t15,t16,t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
    _ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t14,t15,t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
    _ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t14,t15,t16,t17,t18,t19,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
    _ff_sum_7_mults_d2(t4,f[0],f[1],t14,t15,t16,t17,t18,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
    _ff_sum_6_mults_d1(t3,f[0],f[1],t14,t15,t16,t17,g[0],g[1],g[2],g[3],f[1],f[2]);
    _ff_sum_4_mults_d1(t2,f[0],t14,t15,t16,g[0],g[1],g[2],f[1]);
    _ff_sum_3_mults(t1,f[0],t14,t15,g[0],g[1],f[0]);
    _ff_mult(o[0],t14,g[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13);
}

static inline void ff_poly_square_15o (ff_t o[], ff_t f[])
{
	ff_t w[14];
	register ff_t t0;

	ff_poly_mult_7_8(w,f+8,f);
	ff_poly_square_7(o+16,f+8);
	ff_poly_square_8(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[8],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[9],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[10],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[11],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[12],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[13],t0);
	_ff_add(t0,w[6],w[6]); _ff_addto(o[14],t0);
	_ff_add(t0,w[8],w[8]); _ff_addto(o[16],t0);
	_ff_add(t0,w[9],w[9]); _ff_addto(o[17],t0);
	_ff_add(t0,w[10],w[10]); _ff_addto(o[18],t0);
	_ff_add(t0,w[11],w[11]); _ff_addto(o[19],t0);
	_ff_add(t0,w[12],w[12]); _ff_addto(o[20],t0);
	_ff_add(t0,w[13],w[13]); _ff_addto(o[21],t0);
	_ff_add(o[15],w[7],w[7]);
}

static inline void ff_poly_square_15 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_10_mults_s10(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_11_mults_s11(t10,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_12_mults_s12(t11,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_13_mults_s13(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_14_mults_s14(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_15_mults_s15(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_14_mults_s14(o[15],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_13_mults_s13(o[16],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_12_mults_s12(o[17],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_11_mults_s11(o[18],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_10_mults_s10(o[19],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_9_mults_s9(o[20],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_8_mults_s8(o[21],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_7_mults_s7(o[22],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_6_mults_s6(o[23],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_5_mults_s5(o[24],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_4_mults_s4(o[25],f[11],f[12],f[13],f[14]);
	_ff_sum_3_mults_s3(o[26],f[12],f[13],f[14]);
	_ff_dbl_mult(o[27],f[13],f[14]);
	_ff_square(o[0],f[0]); _ff_square(o[28],f[14]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8);
	_ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14);
}

// computes f^2 mod g  for f degree 14 and g(x) = x^15-g13x^13-g12x^142...-g1x-g0 (note signs!)
static inline void ff_poly_square_mod_15 (ff_t o[15], ff_t f[15], ff_t g[14])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28;

	_ff_square(t28,f[14]);
	_ff_dbl_mult(t27,f[13],f[14]);																
	_ff_sum_3_mults_d1 (t26,f[12],f[13],t28,g[13],f[13],f[14]);	
	_ff_sum_4_mults_d2(t25,f[11],f[12],t27,t28,g[12],g[13],f[13],f[14]);	
	_ff_sum_6_mults_d2(t24,f[10],f[11],f[12],t26,t27,t28,g[11],g[12],g[13],f[12],f[13],f[14]);		
	_ff_sum_7_mults_d3(t23,f[9],f[10],f[11],t25,t26,t27,t28,g[10],g[11],g[12],g[13],f[12],f[13],f[14]);					
	_ff_sum_9_mults_d3(t22,f[8],f[9],f[10],f[11],t24,t25,t26,t27,t28,g[9],g[10],g[11],g[12],g[13],f[11],f[12],f[13],f[14]);	
	_ff_sum_10_mults_d4(t21,f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,g[8],g[9],g[10],g[11],g[12],g[13],f[11],f[12],f[13],f[14]);
	_ff_sum_12_mults_d4(t20,f[6],f[7],f[8],f[9],f[10],t22,t23,t24,t25,t26,t27,t28,g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[10],f[11],f[12],f[13],f[14]);	
	_ff_sum_13_mults_d5(t19,f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_15_mults_d5(t18,f[4],f[5],f[6],f[7],f[8],f[9],t20,t21,t22,t23,t24,t25,t26,t27,t28,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_16_mults_d6(t17,f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_18_mults_d6(t16,f[2],f[3],f[4],f[5],f[6],f[7],f[8],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_19_mults_d7(t15,f[1],f[2],f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	
	_ff_sum_21_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t15,t16,t17,t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t15,t16,t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
	_ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t15,t16,t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
	_ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t15,t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);	
	_ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t15,t16,t17,t18,t19,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);		
	_ff_sum_6_mults_d2(t3,f[0],f[1],t15,t16,t17,t18,g[0],g[1],g[2],g[3],f[2],f[3]);						
	_ff_sum_5_mults_d1(t2,f[0],f[1],t15,t16,t17,g[0],g[1],g[2],f[1],f[2]);							
	_ff_sum_3_mults_d1(t1,f[0],t15,t16,g[0],g[1],f[1]);	
	_ff_sum_2_mults(o[0],f[0],t15,g[0],f[0]);												
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10);
	 _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13);  _ff_set(o[14],t14);
}

// Computes x*f^2 mod g for f degree 14, g(x) = x^15-g13x^13-g12x^12-...-g1x-g0 (note signs!)
static inline void ff_poly_square_mult_x_mod_15 (ff_t o[15], ff_t f[15], ff_t g[14])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29;

	_ff_square(t29,f[14]);													
	_ff_dbl_mult(t28,f[13],f[14]);																
	_ff_sum_3_mults_d1 (t27,f[12],f[13],t29,g[13],f[13],f[14]);
	_ff_sum_4_mults_d2(t26,f[11],f[12],t28,t29,g[12],g[13],f[13],f[14]);										
	_ff_sum_6_mults_d2(t25,f[10],f[11],f[12],t27,t28,t29,g[11],g[12],g[13],f[12],f[13],f[14]);					
	_ff_sum_7_mults_d3(t24,f[9],f[10],f[11],t26,t27,t28,t29,g[10],g[11],g[12],g[13],f[12],f[13],f[14]);
	_ff_sum_9_mults_d3(t23,f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,g[9],g[10],g[11],g[12],g[13],f[11],f[12],f[13],f[14]);
	_ff_sum_10_mults_d4(t22,f[7],f[8],f[9],f[10],t24,t25,t26,t27,t28,t29,g[8],g[9],g[10],g[11],g[12],g[13],f[11],f[12],f[13],f[14]);
	_ff_sum_12_mults_d4(t21,f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_13_mults_d5(t20,f[5],f[6],f[7],f[8],f[9],t22,t23,t24,t25,t26,t27,t28,t29,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_15_mults_d5(t19,f[4],f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,t29,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_16_mults_d6(t18,f[3],f[4],f[5],f[6],f[7],f[8],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_18_mults_d6(t17,f[2],f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_19_mults_d7(t16,f[1],f[2],f[3],f[4],f[5],f[6],f[7],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_21_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);

	_ff_sum_21_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t15,t16,t17,t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
	_ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t15,t16,t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);	
	_ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t15,t16,t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
	_ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t15,t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
	_ff_sum_7_mults_d2(t4,f[0],f[1],t15,t16,t17,t18,t19,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);					
	_ff_sum_6_mults_d1(t3,f[0],f[1],t15,t16,t17,t18,g[0],g[1],g[2],g[3],f[1],f[2]);							
	_ff_sum_4_mults_d1(t2,f[0],t15,t16,t17,g[0],g[1],g[2],f[1]);										
	_ff_sum_3_mults(t1,f[0],t15,t16,g[0],g[1],f[0]);													
	_ff_mult(o[0],t15,g[0]);																
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10);
	 _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13);  _ff_set(o[14],t14);
}

static inline void ff_poly_square_16o (ff_t o[], ff_t f[])
{
	ff_t w[15];
	register ff_t t0;
	
	ff_poly_mult_8_8(w,f+8,f);
	ff_poly_square_8(o+16,f+8);
	ff_poly_square_8(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[8],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[9],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[10],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[11],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[12],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[13],t0);
	_ff_add(t0,w[6],w[6]); _ff_addto(o[14],t0);
	_ff_add(t0,w[8],w[8]); _ff_addto(o[16],t0);
	_ff_add(t0,w[9],w[9]); _ff_addto(o[17],t0);
	_ff_add(t0,w[10],w[10]); _ff_addto(o[18],t0);
	_ff_add(t0,w[11],w[11]); _ff_addto(o[19],t0);
	_ff_add(t0,w[12],w[12]); _ff_addto(o[20],t0);
	_ff_add(t0,w[13],w[13]); _ff_addto(o[21],t0);
	_ff_add(t0,w[14],w[14]); _ff_addto(o[22],t0);
	_ff_add(o[15],w[7],w[7]);
}


static inline void ff_poly_square_16 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_10_mults_s10(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_11_mults_s11(t10,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_12_mults_s12(t11,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_13_mults_s13(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_14_mults_s14(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_15_mults_s15(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_16_mults_s16(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_15_mults_s15(o[16],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_14_mults_s14(o[17],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_13_mults_s13(o[18],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_12_mults_s12(o[19],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_11_mults_s11(o[20],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_10_mults_s10(o[21],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_9_mults_s9(o[22],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_8_mults_s8(o[23],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_7_mults_s7(o[24],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_6_mults_s6(o[25],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_5_mults_s5(o[26],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_4_mults_s4(o[27],f[12],f[13],f[14],f[15]);
	_ff_sum_3_mults_s3(o[28],f[13],f[14],f[15]);
	_ff_dbl_mult(o[29],f[14],f[15]);
	_ff_square(o[0],f[0]); _ff_square(o[30],f[15]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8);
	_ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14);  _ff_set(o[15],t15);
}


static inline void ff_poly_square_mod_16 (ff_t o[16], ff_t f[16], ff_t g[15])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30;

    _ff_square(t30,f[15]);
    _ff_dbl_mult(t29,f[14],f[15]);
    _ff_sum_3_mults_d1(t28,f[13],f[14],t30,g[14],f[14],f[15]);
    _ff_sum_4_mults_d2(t27,f[12],f[13],t29,t30,g[13],g[14],f[14],f[15]);
    _ff_sum_6_mults_d2(t26,f[11],f[12],f[13],t28,t29,t30,g[12],g[13],g[14],f[13],f[14],f[15]);
    _ff_sum_7_mults_d3(t25,f[10],f[11],f[12],t27,t28,t29,t30,g[11],g[12],g[13],g[14],f[13],f[14],f[15]);
    _ff_sum_9_mults_d3(t24,f[9],f[10],f[11],f[12],t26,t27,t28,t29,t30,g[10],g[11],g[12],g[13],g[14],f[12],f[13],f[14],f[15]);
    _ff_sum_10_mults_d4(t23,f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,g[9],g[10],g[11],g[12],g[13],g[14],f[12],f[13],f[14],f[15]);
    _ff_sum_12_mults_d4(t22,f[7],f[8],f[9],f[10],f[11],t24,t25,t26,t27,t28,t29,t30,g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_13_mults_d5(t21,f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_15_mults_d5(t20,f[5],f[6],f[7],f[8],f[9],f[10],t22,t23,t24,t25,t26,t27,t28,t29,t30,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_16_mults_d6(t19,f[4],f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_18_mults_d6(t18,f[3],f[4],f[5],f[6],f[7],f[8],f[9],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_19_mults_d7(t17,f[2],f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_21_mults_d7(t16,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);

    _ff_sum_22_mults_d8(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_23_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t16,t17,t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t16,t17,t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t16,t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
    _ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t16,t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
    _ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
    _ff_sum_6_mults_d2(t3,f[0],f[1],t16,t17,t18,t19,g[0],g[1],g[2],g[3],f[2],f[3]);
    _ff_sum_5_mults_d1(t2,f[0],f[1],t16,t17,t18,g[0],g[1],g[2],f[1],f[2]);
    _ff_sum_3_mults_d1(t1,f[0],t16,t17,g[0],g[1],f[1]);
    _ff_sum_2_mults(o[0],f[0],t16,g[0],f[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15);
}

static inline void ff_poly_square_mult_x_mod_16 (ff_t o[16], ff_t f[16], ff_t g[15])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31;

    _ff_square(t31,f[15]);
    _ff_dbl_mult(t30,f[14],f[15]);
    _ff_sum_3_mults_d1(t29,f[13],f[14],t31,g[14],f[14],f[15]);
    _ff_sum_4_mults_d2(t28,f[12],f[13],t30,t31,g[13],g[14],f[14],f[15]);
    _ff_sum_6_mults_d2(t27,f[11],f[12],f[13],t29,t30,t31,g[12],g[13],g[14],f[13],f[14],f[15]);
    _ff_sum_7_mults_d3(t26,f[10],f[11],f[12],t28,t29,t30,t31,g[11],g[12],g[13],g[14],f[13],f[14],f[15]);
    _ff_sum_9_mults_d3(t25,f[9],f[10],f[11],f[12],t27,t28,t29,t30,t31,g[10],g[11],g[12],g[13],g[14],f[12],f[13],f[14],f[15]);
    _ff_sum_10_mults_d4(t24,f[8],f[9],f[10],f[11],t26,t27,t28,t29,t30,t31,g[9],g[10],g[11],g[12],g[13],g[14],f[12],f[13],f[14],f[15]);
    _ff_sum_12_mults_d4(t23,f[7],f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,t31,g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_13_mults_d5(t22,f[6],f[7],f[8],f[9],f[10],t24,t25,t26,t27,t28,t29,t30,t31,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_15_mults_d5(t21,f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_16_mults_d6(t20,f[4],f[5],f[6],f[7],f[8],f[9],t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_18_mults_d6(t19,f[3],f[4],f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_19_mults_d7(t18,f[2],f[3],f[4],f[5],f[6],f[7],f[8],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_21_mults_d7(t17,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_22_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);

    _ff_sum_23_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_22_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t16,t17,t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
    _ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t16,t17,t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
    _ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t16,t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
    _ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t16,t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
    _ff_sum_7_mults_d2(t4,f[0],f[1],t16,t17,t18,t19,t20,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
    _ff_sum_6_mults_d1(t3,f[0],f[1],t16,t17,t18,t19,g[0],g[1],g[2],g[3],f[1],f[2]);
    _ff_sum_4_mults_d1(t2,f[0],t16,t17,t18,g[0],g[1],g[2],f[1]);
    _ff_sum_3_mults(t1,f[0],t16,t17,g[0],g[1],f[0]);
    _ff_mult(o[0],t16,g[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15);
}

static inline void ff_poly_square_17o (ff_t o[], ff_t f[])
{
	ff_t w[16];
	register ff_t t0;
	
	ff_poly_mult_8_9(w,f+9,f);
	ff_poly_square_8(o+18,f+9);
	ff_poly_square_9(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[9],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[10],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[11],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[12],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[13],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[14],t0);
	_ff_add(t0,w[6],w[6]); _ff_addto(o[15],t0);
	_ff_add(t0,w[7],w[7]); _ff_addto(o[16],t0);
	_ff_add(t0,w[9],w[9]); _ff_addto(o[18],t0);
	_ff_add(t0,w[10],w[10]); _ff_addto(o[19],t0);
	_ff_add(t0,w[11],w[11]); _ff_addto(o[20],t0);
	_ff_add(t0,w[12],w[12]); _ff_addto(o[21],t0);
	_ff_add(t0,w[13],w[13]); _ff_addto(o[22],t0);
	_ff_add(t0,w[14],w[14]); _ff_addto(o[23],t0);
	_ff_add(t0,w[15],w[15]); _ff_addto(o[24],t0);
	_ff_add(o[17],w[8],w[8]);
}


static inline void ff_poly_square_17 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_10_mults_s10(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_11_mults_s11(t10,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_12_mults_s12(t11,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_13_mults_s13(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_14_mults_s14(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_15_mults_s15(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_16_mults_s16(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_17_mults_s17(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_16_mults_s16(o[17],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_15_mults_s15(o[18],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_14_mults_s14(o[19],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_13_mults_s13(o[20],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_12_mults_s12(o[21],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_11_mults_s11(o[22],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_10_mults_s10(o[23],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_9_mults_s9(o[24],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_8_mults_s8(o[25],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_7_mults_s7(o[26],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_6_mults_s6(o[27],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_5_mults_s5(o[28],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_4_mults_s4(o[29],f[13],f[14],f[15],f[16]);
	_ff_sum_3_mults_s3(o[30],f[14],f[15],f[16]);
	_ff_dbl_mult(o[31],f[15],f[16]);
	_ff_square(o[0],f[0]); _ff_square(o[32],f[16]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8);
	_ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14);  _ff_set(o[15],t15); _ff_set(o[16],t16);
}

// computes f^2 mod g  for f degree 16 and g(x) = x^17-g15x^15-g14x^14-...-g1x-g0 (note signs!)
static inline void ff_poly_square_mod_17 (ff_t o[17], ff_t f[17], ff_t g[16])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
	
	_ff_square(t32,f[16]);
	_ff_dbl_mult(t31,f[15],f[16]);																
	_ff_sum_3_mults_d1 (t30,f[14],f[15],t32,g[15],f[15],f[16]);	
	_ff_sum_4_mults_d2(t29,f[13],f[14],t31,t32,g[14],g[15],f[15],f[16]);	
	_ff_sum_6_mults_d2(t28,f[12],f[13],f[14],t30,t31,t32,g[13],g[14],g[15],f[14],f[15],f[16]);		
	_ff_sum_7_mults_d3(t27,f[11],f[12],f[13],t29,t30,t31,t32,g[12],g[13],g[14],g[15],f[14],f[15],f[16]);					
	_ff_sum_9_mults_d3(t26,f[10],f[11],f[12],f[13],t28,t29,t30,t31,t32,g[11],g[12],g[13],g[14],g[15],f[13],f[14],f[15],f[16]);	
	_ff_sum_10_mults_d4(t25,f[9],f[10],f[11],f[12],t27,t28,t29,t30,t31,t32,g[10],g[11],g[12],g[13],g[14],g[15],f[13],f[14],f[15],f[16]);
	_ff_sum_12_mults_d4(t24,f[8],f[9],f[10],f[11],f[12],t26,t27,t28,t29,t30,t31,t32,g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[12],f[13],f[14],f[15],f[16]);	
	_ff_sum_13_mults_d5(t23,f[7],f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,t31,t32,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_15_mults_d5(t22,f[6],f[7],f[8],f[9],f[10],f[11],t24,t25,t26,t27,t28,t29,t30,t31,t32,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_16_mults_d6(t21,f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_18_mults_d6(t20,f[4],f[5],f[6],f[7],f[8],f[9],f[10],t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_19_mults_d7(t19,f[3],f[4],f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_21_mults_d7(t18,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_22_mults_d8(t17,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	
	_ff_sum_24_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_24_mults_d8(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_23_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t17,t18,t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t17,t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
	_ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t17,t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
	_ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);	
	_ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);		
	_ff_sum_6_mults_d2(t3,f[0],f[1],t17,t18,t19,t20,g[0],g[1],g[2],g[3],f[2],f[3]);						
	_ff_sum_5_mults_d1(t2,f[0],f[1],t17,t18,t19,g[0],g[1],g[2],f[1],f[2]);							
	_ff_sum_3_mults_d1(t1,f[0],t17,t18,g[0],g[1],f[1]);	
	_ff_sum_2_mults(o[0],f[0],t17,g[0],f[0]);												
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10);
	 _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13);  _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16);
}

// Computes x*f^2 mod g for f degree 16, g(x) = x^17-g15x^15-g14x^14-...-g1x-g0 (note signs!)
static inline void ff_poly_square_mult_x_mod_17 (ff_t o[17], ff_t f[17], ff_t g[16])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33;

	_ff_square(t33,f[16]);													
	_ff_dbl_mult(t32,f[15],f[16]);																
	_ff_sum_3_mults_d1 (t31,f[14],f[15],t33,g[15],f[15],f[16]);		
	_ff_sum_4_mults_d2(t30,f[13],f[14],t32,t33,g[14],g[15],f[15],f[16]);										
	_ff_sum_6_mults_d2(t29,f[12],f[13],f[14],t31,t32,t33,g[13],g[14],g[15],f[14],f[15],f[16]);					
	_ff_sum_7_mults_d3(t28,f[11],f[12],f[13],t30,t31,t32,t33,g[12],g[13],g[14],g[15],f[14],f[15],f[16]);
	_ff_sum_9_mults_d3(t27,f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,g[11],g[12],g[13],g[14],g[15],f[13],f[14],f[15],f[16]);
	_ff_sum_10_mults_d4(t26,f[9],f[10],f[11],f[12],t28,t29,t30,t31,t32,t33,g[10],g[11],g[12],g[13],g[14],g[15],f[13],f[14],f[15],f[16]);
	_ff_sum_12_mults_d4(t25,f[8],f[9],f[10],f[11],f[12],t27,t28,t29,t30,t31,t32,t33,g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_13_mults_d5(t24,f[7],f[8],f[9],f[10],f[11],t26,t27,t28,t29,t30,t31,t32,t33,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_15_mults_d5(t23,f[6],f[7],f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,t31,t32,t33,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_16_mults_d6(t22,f[5],f[6],f[7],f[8],f[9],f[10],t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_18_mults_d6(t21,f[4],f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_19_mults_d7(t20,f[3],f[4],f[5],f[6],f[7],f[8],f[9],t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_21_mults_d7(t19,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_22_mults_d8(t18,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_24_mults_d8(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	
	_ff_sum_24_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_24_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_22_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t17,t18,t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
	_ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t17,t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);	
	_ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t17,t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
	_ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t17,t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
	_ff_sum_7_mults_d2(t4,f[0],f[1],t17,t18,t19,t20,t21,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);					
	_ff_sum_6_mults_d1(t3,f[0],f[1],t17,t18,t19,t20,g[0],g[1],g[2],g[3],f[1],f[2]);							
	_ff_sum_4_mults_d1(t2,f[0],t17,t18,t19,g[0],g[1],g[2],f[1]);										
	_ff_sum_3_mults(t1,f[0],t17,t18,g[0],g[1],f[0]);													
	_ff_mult(o[0],t17,g[0]);																
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10);
	 _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13);  _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16);
}

static inline void ff_poly_square_18o (ff_t o[], ff_t f[])
{
	ff_t w[17];
	register ff_t t0;

	ff_poly_mult_9_9(w,f+9,f);
	ff_poly_square_9(o+18,f+9);
	ff_poly_square_9(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[9],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[10],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[11],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[12],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[13],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[14],t0);
	_ff_add(t0,w[6],w[6]); _ff_addto(o[15],t0);
	_ff_add(t0,w[7],w[7]); _ff_addto(o[16],t0);
	_ff_add(t0,w[9],w[9]); _ff_addto(o[18],t0);
	_ff_add(t0,w[10],w[10]); _ff_addto(o[19],t0);
	_ff_add(t0,w[11],w[11]); _ff_addto(o[20],t0);
	_ff_add(t0,w[12],w[12]); _ff_addto(o[21],t0);
	_ff_add(t0,w[13],w[13]); _ff_addto(o[22],t0);
	_ff_add(t0,w[14],w[14]); _ff_addto(o[23],t0);
	_ff_add(t0,w[15],w[15]); _ff_addto(o[24],t0);
	_ff_add(t0,w[16],w[16]); _ff_addto(o[25],t0);
	_ff_add(o[17],w[8],w[8]);
}


static inline void ff_poly_square_18 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_10_mults_s10(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_11_mults_s11(t10,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_12_mults_s12(t11,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_13_mults_s13(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_14_mults_s14(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_15_mults_s15(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_16_mults_s16(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_17_mults_s17(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_18_mults_s18(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_17_mults_s17(o[18],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_16_mults_s16(o[19],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_15_mults_s15(o[20],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_14_mults_s14(o[21],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_13_mults_s13(o[22],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_12_mults_s12(o[23],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_11_mults_s11(o[24],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_10_mults_s10(o[25],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_9_mults_s9(o[26],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_8_mults_s8(o[27],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_7_mults_s7(o[28],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_6_mults_s6(o[29],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_5_mults_s5(o[30],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_4_mults_s4(o[31],f[14],f[15],f[16],f[17]);
	_ff_sum_3_mults_s3(o[32],f[15],f[16],f[17]);
	_ff_dbl_mult(o[33],f[16],f[17]);
	_ff_square(o[0],f[0]); _ff_square(o[34],f[17]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8);
	_ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14);  _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17);
}

static inline void ff_poly_square_mod_18 (ff_t o[18], ff_t f[18], ff_t g[17])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34;

    _ff_square(t34,f[17]);
    _ff_dbl_mult(t33,f[16],f[17]);
    _ff_sum_3_mults_d1(t32,f[15],f[16],t34,g[16],f[16],f[17]);
    _ff_sum_4_mults_d2(t31,f[14],f[15],t33,t34,g[15],g[16],f[16],f[17]);
    _ff_sum_6_mults_d2(t30,f[13],f[14],f[15],t32,t33,t34,g[14],g[15],g[16],f[15],f[16],f[17]);
    _ff_sum_7_mults_d3(t29,f[12],f[13],f[14],t31,t32,t33,t34,g[13],g[14],g[15],g[16],f[15],f[16],f[17]);
    _ff_sum_9_mults_d3(t28,f[11],f[12],f[13],f[14],t30,t31,t32,t33,t34,g[12],g[13],g[14],g[15],g[16],f[14],f[15],f[16],f[17]);
    _ff_sum_10_mults_d4(t27,f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,t34,g[11],g[12],g[13],g[14],g[15],g[16],f[14],f[15],f[16],f[17]);
    _ff_sum_12_mults_d4(t26,f[9],f[10],f[11],f[12],f[13],t28,t29,t30,t31,t32,t33,t34,g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_13_mults_d5(t25,f[8],f[9],f[10],f[11],f[12],t27,t28,t29,t30,t31,t32,t33,t34,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_15_mults_d5(t24,f[7],f[8],f[9],f[10],f[11],f[12],t26,t27,t28,t29,t30,t31,t32,t33,t34,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_16_mults_d6(t23,f[6],f[7],f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_18_mults_d6(t22,f[5],f[6],f[7],f[8],f[9],f[10],f[11],t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_19_mults_d7(t21,f[4],f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_21_mults_d7(t20,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_22_mults_d8(t19,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_24_mults_d8(t18,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);

    _ff_sum_25_mults_d9(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_26_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_24_mults_d8(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_23_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t18,t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t18,t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
    _ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
    _ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
    _ff_sum_6_mults_d2(t3,f[0],f[1],t18,t19,t20,t21,g[0],g[1],g[2],g[3],f[2],f[3]);
    _ff_sum_5_mults_d1(t2,f[0],f[1],t18,t19,t20,g[0],g[1],g[2],f[1],f[2]);
    _ff_sum_3_mults_d1(t1,f[0],t18,t19,g[0],g[1],f[1]);
    _ff_sum_2_mults(o[0],f[0],t18,g[0],f[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17);
}

static inline void ff_poly_square_mult_x_mod_18 (ff_t o[18], ff_t f[18], ff_t g[17])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35;

    _ff_square(t35,f[17]);
    _ff_dbl_mult(t34,f[16],f[17]);
    _ff_sum_3_mults_d1(t33,f[15],f[16],t35,g[16],f[16],f[17]);
    _ff_sum_4_mults_d2(t32,f[14],f[15],t34,t35,g[15],g[16],f[16],f[17]);
    _ff_sum_6_mults_d2(t31,f[13],f[14],f[15],t33,t34,t35,g[14],g[15],g[16],f[15],f[16],f[17]);
    _ff_sum_7_mults_d3(t30,f[12],f[13],f[14],t32,t33,t34,t35,g[13],g[14],g[15],g[16],f[15],f[16],f[17]);
    _ff_sum_9_mults_d3(t29,f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,g[12],g[13],g[14],g[15],g[16],f[14],f[15],f[16],f[17]);
    _ff_sum_10_mults_d4(t28,f[10],f[11],f[12],f[13],t30,t31,t32,t33,t34,t35,g[11],g[12],g[13],g[14],g[15],g[16],f[14],f[15],f[16],f[17]);
    _ff_sum_12_mults_d4(t27,f[9],f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,t34,t35,g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_13_mults_d5(t26,f[8],f[9],f[10],f[11],f[12],t28,t29,t30,t31,t32,t33,t34,t35,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_15_mults_d5(t25,f[7],f[8],f[9],f[10],f[11],f[12],t27,t28,t29,t30,t31,t32,t33,t34,t35,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_16_mults_d6(t24,f[6],f[7],f[8],f[9],f[10],f[11],t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_18_mults_d6(t23,f[5],f[6],f[7],f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_19_mults_d7(t22,f[4],f[5],f[6],f[7],f[8],f[9],f[10],t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_21_mults_d7(t21,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_22_mults_d8(t20,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_24_mults_d8(t19,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_25_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);

    _ff_sum_26_mults_d8(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_25_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_24_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_22_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t18,t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
    _ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t18,t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
    _ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t18,t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
    _ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t18,t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
    _ff_sum_7_mults_d2(t4,f[0],f[1],t18,t19,t20,t21,t22,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
    _ff_sum_6_mults_d1(t3,f[0],f[1],t18,t19,t20,t21,g[0],g[1],g[2],g[3],f[1],f[2]);
    _ff_sum_4_mults_d1(t2,f[0],t18,t19,t20,g[0],g[1],g[2],f[1]);
    _ff_sum_3_mults(t1,f[0],t18,t19,g[0],g[1],f[0]);
    _ff_mult(o[0],t18,g[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17);
}

static inline void ff_poly_square_19o (ff_t o[], ff_t f[])
{
	ff_t w[18];
	register ff_t t0;

	ff_poly_mult_9_10(w,f+10,f);
	ff_poly_square_9(o+20,f+10);
	ff_poly_square_10(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[10],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[11],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[12],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[13],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[14],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[15],t0);
	_ff_add(t0,w[6],w[6]); _ff_addto(o[16],t0);
	_ff_add(t0,w[7],w[7]); _ff_addto(o[17],t0);
	_ff_add(t0,w[8],w[8]); _ff_addto(o[18],t0);
	_ff_add(t0,w[10],w[10]); _ff_addto(o[20],t0);
	_ff_add(t0,w[11],w[11]); _ff_addto(o[21],t0);
	_ff_add(t0,w[12],w[12]); _ff_addto(o[22],t0);
	_ff_add(t0,w[13],w[13]); _ff_addto(o[23],t0);
	_ff_add(t0,w[14],w[14]); _ff_addto(o[24],t0);
	_ff_add(t0,w[15],w[15]); _ff_addto(o[25],t0);
	_ff_add(t0,w[16],w[16]); _ff_addto(o[26],t0);
	_ff_add(t0,w[17],w[17]); _ff_addto(o[27],t0);
	_ff_add(o[19],w[9],w[9]);
}

static inline void ff_poly_square_19 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_10_mults_s10(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_11_mults_s11(t10,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_12_mults_s12(t11,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_13_mults_s13(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_14_mults_s14(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_15_mults_s15(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_16_mults_s16(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_17_mults_s17(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_18_mults_s18(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_19_mults_s19(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_18_mults_s18(o[19],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_17_mults_s17(o[20],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_16_mults_s16(o[21],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_15_mults_s15(o[22],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_14_mults_s14(o[23],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_13_mults_s13(o[24],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_12_mults_s12(o[25],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_11_mults_s11(o[26],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_10_mults_s10(o[27],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_9_mults_s9(o[28],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_8_mults_s8(o[29],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_7_mults_s7(o[30],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_6_mults_s6(o[31],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_5_mults_s5(o[32],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_4_mults_s4(o[33],f[15],f[16],f[17],f[18]);
	_ff_sum_3_mults_s3(o[34],f[16],f[17],f[18]);
	_ff_dbl_mult(o[35],f[17],f[18]);
	_ff_square(o[0],f[0]); _ff_square(o[36],f[18]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8);  _ff_set(o[9],t9); _ff_set(o[10],t10);
	_ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14);  _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18);
}


// computes f^2 mod g  for f degree 18 and g(x) = x^19-g17x^17-g16x^16-...-g1x-g0 (note signs!)
static inline void ff_poly_square_mod_19 (ff_t o[19], ff_t f[19], ff_t g[18])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36;
	
	_ff_square(t36,f[18]);
	_ff_dbl_mult(t35,f[17],f[18]);																
	_ff_sum_3_mults_d1 (t34,f[16],f[17],t36,g[17],f[17],f[18]);	
	_ff_sum_4_mults_d2(t33,f[15],f[16],t35,t36,g[16],g[17],f[17],f[18]);	
	_ff_sum_6_mults_d2(t32,f[14],f[15],f[16],t34,t35,t36,g[15],g[16],g[17],f[16],f[17],f[18]);		
	_ff_sum_7_mults_d3(t31,f[13],f[14],f[15],t33,t34,t35,t36,g[14],g[15],g[16],g[17],f[16],f[17],f[18]);					
	_ff_sum_9_mults_d3(t30,f[12],f[13],f[14],f[15],t32,t33,t34,t35,t36,g[13],g[14],g[15],g[16],g[17],f[15],f[16],f[17],f[18]);	
	_ff_sum_10_mults_d4(t29,f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,t36,g[12],g[13],g[14],g[15],g[16],g[17],f[15],f[16],f[17],f[18]);
	_ff_sum_12_mults_d4(t28,f[10],f[11],f[12],f[13],f[14],t30,t31,t32,t33,t34,t35,t36,g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[14],f[15],f[16],f[17],f[18]);	
	_ff_sum_13_mults_d5(t27,f[9],f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,t34,t35,t36,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_15_mults_d5(t26,f[8],f[9],f[10],f[11],f[12],f[13],t28,t29,t30,t31,t32,t33,t34,t35,t36,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_16_mults_d6(t25,f[7],f[8],f[9],f[10],f[11],f[12],t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_18_mults_d6(t24,f[6],f[7],f[8],f[9],f[10],f[11],f[12],t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_19_mults_d7(t23,f[5],f[6],f[7],f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_21_mults_d7(t22,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_22_mults_d8(t21,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_24_mults_d8(t20,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_25_mults_d9(t19,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);

	_ff_sum_27_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_27_mults_d9(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_26_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_24_mults_d8(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_23_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t19,t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
	_ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
	_ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);	
	_ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);		
	_ff_sum_6_mults_d2(t3,f[0],f[1],t19,t20,t21,t22,g[0],g[1],g[2],g[3],f[2],f[3]);						
	_ff_sum_5_mults_d1(t2,f[0],f[1],t19,t20,t21,g[0],g[1],g[2],f[1],f[2]);							
	_ff_sum_3_mults_d1(t1,f[0],t19,t20,g[0],g[1],f[1]);	
	_ff_sum_2_mults(o[0],f[0],t19,g[0],f[0]);												
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10);
	 _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13);  _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18);
}

// Computes x*f^2 mod g for f degree 18, g(x) = x^19-g17x^17-g16x^16-...-g1x-g0 (note signs!)
static inline void ff_poly_square_mult_x_mod_19 (ff_t o[19], ff_t f[19], ff_t g[18])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37;

	_ff_square(t37,f[18]);
	_ff_dbl_mult(t36,f[17],f[18]);																
	_ff_sum_3_mults_d1 (t35,f[16],f[17],t37,g[17],f[17],f[18]);		
	_ff_sum_4_mults_d2(t34,f[15],f[16],t36,t37,g[16],g[17],f[17],f[18]);										
	_ff_sum_6_mults_d2(t33,f[14],f[15],f[16],t35,t36,t37,g[15],g[16],g[17],f[16],f[17],f[18]);					
	_ff_sum_7_mults_d3(t32,f[13],f[14],f[15],t34,t35,t36,t37,g[14],g[15],g[16],g[17],f[16],f[17],f[18]);
	_ff_sum_9_mults_d3(t31,f[12],f[13],f[14],f[15],t33,t34,t35,t36,t37,g[13],g[14],g[15],g[16],g[17],f[15],f[16],f[17],f[18]);
	_ff_sum_10_mults_d4(t30,f[11],f[12],f[13],f[14],t32,t33,t34,t35,t36,t37,g[12],g[13],g[14],g[15],g[16],g[17],f[15],f[16],f[17],f[18]);
	_ff_sum_12_mults_d4(t29,f[10],f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,t36,t37,g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_13_mults_d5(t28,f[9],f[10],f[11],f[12],f[13],t30,t31,t32,t33,t34,t35,t36,t37,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_15_mults_d5(t27,f[8],f[9],f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,t34,t35,t36,t37,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_16_mults_d6(t26,f[7],f[8],f[9],f[10],f[11],f[12],t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_18_mults_d6(t25,f[6],f[7],f[8],f[9],f[10],f[11],f[12],t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_19_mults_d7(t24,f[5],f[6],f[7],f[8],f[9],f[10],f[11],t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_21_mults_d7(t23,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_22_mults_d8(t22,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_24_mults_d8(t21,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_25_mults_d9(t20,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_27_mults_d9(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	
	_ff_sum_27_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_27_mults_d8(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_25_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_24_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_22_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t19,t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
	_ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t19,t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);	
	_ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t19,t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
	_ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t19,t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
	_ff_sum_7_mults_d2(t4,f[0],f[1],t19,t20,t21,t22,t23,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);					
	_ff_sum_6_mults_d1(t3,f[0],f[1],t19,t20,t21,t22,g[0],g[1],g[2],g[3],f[1],f[2]);							
	_ff_sum_4_mults_d1(t2,f[0],t19,t20,t21,g[0],g[1],g[2],f[1]);										
	_ff_sum_3_mults(t1,f[0],t19,t20,g[0],g[1],f[0]);													
	_ff_mult(o[0],t19,g[0]);																
	_ff_set(o[1],t1);  _ff_set(o[2],t2);  _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10);
	 _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13);  _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18);
}

static inline void ff_poly_square_20o (ff_t o[], ff_t f[])
{
	ff_t w[19];
	register ff_t t0;
	
	ff_poly_mult_10_10(w,f+10,f);
	ff_poly_square_10(o+20,f+10);
	ff_poly_square_10(o,f);
	// it is worth unwinding this loop
	_ff_add(t0,w[0],w[0]); _ff_addto(o[10],t0);
	_ff_add(t0,w[1],w[1]); _ff_addto(o[11],t0);
	_ff_add(t0,w[2],w[2]); _ff_addto(o[12],t0);
	_ff_add(t0,w[3],w[3]); _ff_addto(o[13],t0);
	_ff_add(t0,w[4],w[4]); _ff_addto(o[14],t0);
	_ff_add(t0,w[5],w[5]); _ff_addto(o[15],t0);
	_ff_add(t0,w[6],w[6]); _ff_addto(o[16],t0);
	_ff_add(t0,w[7],w[7]); _ff_addto(o[17],t0);
	_ff_add(t0,w[8],w[8]); _ff_addto(o[18],t0);
	_ff_add(t0,w[10],w[10]); _ff_addto(o[20],t0);
	_ff_add(t0,w[11],w[11]); _ff_addto(o[21],t0);
	_ff_add(t0,w[12],w[12]); _ff_addto(o[22],t0);
	_ff_add(t0,w[13],w[13]); _ff_addto(o[23],t0);
	_ff_add(t0,w[14],w[14]); _ff_addto(o[24],t0);
	_ff_add(t0,w[15],w[15]); _ff_addto(o[25],t0);
	_ff_add(t0,w[16],w[16]); _ff_addto(o[26],t0);
	_ff_add(t0,w[17],w[17]); _ff_addto(o[27],t0);
	_ff_add(t0,w[18],w[18]); _ff_addto(o[28],t0);
	_ff_add(o[19],w[9],w[9]);
}

static inline void ff_poly_square_20 (ff_t o[], ff_t f[])
{
	register ff_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19;
	
	_ff_dbl_mult(t1,f[0],f[1]);
	_ff_sum_3_mults_s3(t2,f[0],f[1],f[2]);
	_ff_sum_4_mults_s4(t3,f[0],f[1],f[2],f[3]);
	_ff_sum_5_mults_s5(t4,f[0],f[1],f[2],f[3],f[4]);
	_ff_sum_6_mults_s6(t5,f[0],f[1],f[2],f[3],f[4],f[5]);
	_ff_sum_7_mults_s7(t6,f[0],f[1],f[2],f[3],f[4],f[5],f[6]);
	_ff_sum_8_mults_s8(t7,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);
	_ff_sum_9_mults_s9(t8,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_10_mults_s10(t9,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_11_mults_s11(t10,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_12_mults_s12(t11,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_13_mults_s13(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_14_mults_s14(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_15_mults_s15(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_16_mults_s16(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_17_mults_s17(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_18_mults_s18(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_19_mults_s19(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_20_mults_s20(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_19_mults_s19(o[20],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_18_mults_s18(o[21],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_17_mults_s17(o[22],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_16_mults_s16(o[23],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_15_mults_s15(o[24],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_14_mults_s14(o[25],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_13_mults_s13(o[26],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_12_mults_s12(o[27],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_11_mults_s11(o[28],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_10_mults_s10(o[29],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_9_mults_s9(o[30],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_8_mults_s8(o[31],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_7_mults_s7(o[32],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_6_mults_s6(o[33],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_5_mults_s5(o[34],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_4_mults_s4(o[35],f[16],f[17],f[18],f[19]);
	_ff_sum_3_mults_s3(o[36],f[17],f[18],f[19]);
	_ff_dbl_mult(o[37],f[18],f[19]);
	_ff_square(o[0],f[0]); _ff_square(o[38],f[19]);
	_ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8);  _ff_set(o[9],t9); _ff_set(o[10],t10);
	_ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14);  _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19);
}

static inline void ff_poly_square_mod_20 (ff_t o[20], ff_t f[20], ff_t g[19])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38;

    _ff_square(t38,f[19]);
    _ff_dbl_mult(t37,f[18],f[19]);
    _ff_sum_3_mults_d1(t36,f[17],f[18],t38,g[18],f[18],f[19]);
    _ff_sum_4_mults_d2(t35,f[16],f[17],t37,t38,g[17],g[18],f[18],f[19]);
    _ff_sum_6_mults_d2(t34,f[15],f[16],f[17],t36,t37,t38,g[16],g[17],g[18],f[17],f[18],f[19]);
    _ff_sum_7_mults_d3(t33,f[14],f[15],f[16],t35,t36,t37,t38,g[15],g[16],g[17],g[18],f[17],f[18],f[19]);
    _ff_sum_9_mults_d3(t32,f[13],f[14],f[15],f[16],t34,t35,t36,t37,t38,g[14],g[15],g[16],g[17],g[18],f[16],f[17],f[18],f[19]);
    _ff_sum_10_mults_d4(t31,f[12],f[13],f[14],f[15],t33,t34,t35,t36,t37,t38,g[13],g[14],g[15],g[16],g[17],g[18],f[16],f[17],f[18],f[19]);
    _ff_sum_12_mults_d4(t30,f[11],f[12],f[13],f[14],f[15],t32,t33,t34,t35,t36,t37,t38,g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_13_mults_d5(t29,f[10],f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,t36,t37,t38,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_15_mults_d5(t28,f[9],f[10],f[11],f[12],f[13],f[14],t30,t31,t32,t33,t34,t35,t36,t37,t38,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_16_mults_d6(t27,f[8],f[9],f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_18_mults_d6(t26,f[7],f[8],f[9],f[10],f[11],f[12],f[13],t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_19_mults_d7(t25,f[6],f[7],f[8],f[9],f[10],f[11],f[12],t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_21_mults_d7(t24,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_22_mults_d8(t23,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_24_mults_d8(t22,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_25_mults_d9(t21,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_27_mults_d9(t20,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);

    _ff_sum_28_mults_d10(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_29_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_27_mults_d9(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_26_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_24_mults_d8(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_23_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
    _ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
    _ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
    _ff_sum_6_mults_d2(t3,f[0],f[1],t20,t21,t22,t23,g[0],g[1],g[2],g[3],f[2],f[3]);
    _ff_sum_5_mults_d1(t2,f[0],f[1],t20,t21,t22,g[0],g[1],g[2],f[1],f[2]);
    _ff_sum_3_mults_d1(t1,f[0],t20,t21,g[0],g[1],f[1]);
    _ff_sum_2_mults(o[0],f[0],t20,g[0],f[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19);
}

static inline void ff_poly_square_mult_x_mod_20 (ff_t o[20], ff_t f[20], ff_t g[19])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39;

    _ff_square(t39,f[19]);
    _ff_dbl_mult(t38,f[18],f[19]);
    _ff_sum_3_mults_d1(t37,f[17],f[18],t39,g[18],f[18],f[19]);
    _ff_sum_4_mults_d2(t36,f[16],f[17],t38,t39,g[17],g[18],f[18],f[19]);
    _ff_sum_6_mults_d2(t35,f[15],f[16],f[17],t37,t38,t39,g[16],g[17],g[18],f[17],f[18],f[19]);
    _ff_sum_7_mults_d3(t34,f[14],f[15],f[16],t36,t37,t38,t39,g[15],g[16],g[17],g[18],f[17],f[18],f[19]);
    _ff_sum_9_mults_d3(t33,f[13],f[14],f[15],f[16],t35,t36,t37,t38,t39,g[14],g[15],g[16],g[17],g[18],f[16],f[17],f[18],f[19]);
    _ff_sum_10_mults_d4(t32,f[12],f[13],f[14],f[15],t34,t35,t36,t37,t38,t39,g[13],g[14],g[15],g[16],g[17],g[18],f[16],f[17],f[18],f[19]);
    _ff_sum_12_mults_d4(t31,f[11],f[12],f[13],f[14],f[15],t33,t34,t35,t36,t37,t38,t39,g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_13_mults_d5(t30,f[10],f[11],f[12],f[13],f[14],t32,t33,t34,t35,t36,t37,t38,t39,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_15_mults_d5(t29,f[9],f[10],f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,t36,t37,t38,t39,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_16_mults_d6(t28,f[8],f[9],f[10],f[11],f[12],f[13],t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_18_mults_d6(t27,f[7],f[8],f[9],f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_19_mults_d7(t26,f[6],f[7],f[8],f[9],f[10],f[11],f[12],t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_21_mults_d7(t25,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_22_mults_d8(t24,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_24_mults_d8(t23,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_25_mults_d9(t22,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_27_mults_d9(t21,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_28_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);

    _ff_sum_29_mults_d9(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_28_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_27_mults_d8(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_25_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_24_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_22_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t20,t21,t22,t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
    _ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t20,t21,t22,t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
    _ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t20,t21,t22,t23,t24,t25,t26,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
    _ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t20,t21,t22,t23,t24,t25,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
    _ff_sum_7_mults_d2(t4,f[0],f[1],t20,t21,t22,t23,t24,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
    _ff_sum_6_mults_d1(t3,f[0],f[1],t20,t21,t22,t23,g[0],g[1],g[2],g[3],f[1],f[2]);
    _ff_sum_4_mults_d1(t2,f[0],t20,t21,t22,g[0],g[1],g[2],f[1]);
    _ff_sum_3_mults(t1,f[0],t20,t21,g[0],g[1],f[0]);
    _ff_mult(o[0],t20,g[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19);
}

static inline void ff_poly_square_mod_23 (ff_t o[23], ff_t f[23], ff_t g[22])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44;

    _ff_square(t44,f[22]);
    _ff_dbl_mult(t43,f[21],f[22]);
    _ff_sum_3_mults_d1(t42,f[20],f[21],t44,g[21],f[21],f[22]);
    _ff_sum_4_mults_d2(t41,f[19],f[20],t43,t44,g[20],g[21],f[21],f[22]);
    _ff_sum_6_mults_d2(t40,f[18],f[19],f[20],t42,t43,t44,g[19],g[20],g[21],f[20],f[21],f[22]);
    _ff_sum_7_mults_d3(t39,f[17],f[18],f[19],t41,t42,t43,t44,g[18],g[19],g[20],g[21],f[20],f[21],f[22]);
    _ff_sum_9_mults_d3(t38,f[16],f[17],f[18],f[19],t40,t41,t42,t43,t44,g[17],g[18],g[19],g[20],g[21],f[19],f[20],f[21],f[22]);
    _ff_sum_10_mults_d4(t37,f[15],f[16],f[17],f[18],t39,t40,t41,t42,t43,t44,g[16],g[17],g[18],g[19],g[20],g[21],f[19],f[20],f[21],f[22]);
    _ff_sum_12_mults_d4(t36,f[14],f[15],f[16],f[17],f[18],t38,t39,t40,t41,t42,t43,t44,g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_13_mults_d5(t35,f[13],f[14],f[15],f[16],f[17],t37,t38,t39,t40,t41,t42,t43,t44,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_15_mults_d5(t34,f[12],f[13],f[14],f[15],f[16],f[17],t36,t37,t38,t39,t40,t41,t42,t43,t44,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_16_mults_d6(t33,f[11],f[12],f[13],f[14],f[15],f[16],t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_18_mults_d6(t32,f[10],f[11],f[12],f[13],f[14],f[15],f[16],t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_19_mults_d7(t31,f[9],f[10],f[11],f[12],f[13],f[14],f[15],t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_21_mults_d7(t30,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_22_mults_d8(t29,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_24_mults_d8(t28,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_25_mults_d9(t27,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_27_mults_d9(t26,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_28_mults_d10(t25,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_30_mults_d10(t24,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_31_mults_d11(t23,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);

    _ff_sum_33_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_33_mults_d11(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_32_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_30_mults_d10(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_29_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_27_mults_d9(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_26_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_24_mults_d8(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_23_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t23,t24,t25,t26,t27,t28,t29,t30,t31,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
    _ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
    _ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
    _ff_sum_6_mults_d2(t3,f[0],f[1],t23,t24,t25,t26,g[0],g[1],g[2],g[3],f[2],f[3]);
    _ff_sum_5_mults_d1(t2,f[0],f[1],t23,t24,t25,g[0],g[1],g[2],f[1],f[2]);
    _ff_sum_3_mults_d1(t1,f[0],t23,t24,g[0],g[1],f[1]);
    _ff_sum_2_mults(o[0],f[0],t23,g[0],f[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22);
}

static inline void ff_poly_square_mult_x_mod_23 (ff_t o[23], ff_t f[23], ff_t g[22])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45;

    _ff_square(t45,f[22]);
    _ff_dbl_mult(t44,f[21],f[22]);
    _ff_sum_3_mults_d1(t43,f[20],f[21],t45,g[21],f[21],f[22]);
    _ff_sum_4_mults_d2(t42,f[19],f[20],t44,t45,g[20],g[21],f[21],f[22]);
    _ff_sum_6_mults_d2(t41,f[18],f[19],f[20],t43,t44,t45,g[19],g[20],g[21],f[20],f[21],f[22]);
    _ff_sum_7_mults_d3(t40,f[17],f[18],f[19],t42,t43,t44,t45,g[18],g[19],g[20],g[21],f[20],f[21],f[22]);
    _ff_sum_9_mults_d3(t39,f[16],f[17],f[18],f[19],t41,t42,t43,t44,t45,g[17],g[18],g[19],g[20],g[21],f[19],f[20],f[21],f[22]);
    _ff_sum_10_mults_d4(t38,f[15],f[16],f[17],f[18],t40,t41,t42,t43,t44,t45,g[16],g[17],g[18],g[19],g[20],g[21],f[19],f[20],f[21],f[22]);
    _ff_sum_12_mults_d4(t37,f[14],f[15],f[16],f[17],f[18],t39,t40,t41,t42,t43,t44,t45,g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_13_mults_d5(t36,f[13],f[14],f[15],f[16],f[17],t38,t39,t40,t41,t42,t43,t44,t45,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_15_mults_d5(t35,f[12],f[13],f[14],f[15],f[16],f[17],t37,t38,t39,t40,t41,t42,t43,t44,t45,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_16_mults_d6(t34,f[11],f[12],f[13],f[14],f[15],f[16],t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_18_mults_d6(t33,f[10],f[11],f[12],f[13],f[14],f[15],f[16],t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_19_mults_d7(t32,f[9],f[10],f[11],f[12],f[13],f[14],f[15],t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_21_mults_d7(t31,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_22_mults_d8(t30,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_24_mults_d8(t29,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_25_mults_d9(t28,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_27_mults_d9(t27,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_28_mults_d10(t26,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_30_mults_d10(t25,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_31_mults_d11(t24,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_33_mults_d11(t23,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);

    _ff_sum_33_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_33_mults_d10(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_31_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_30_mults_d9(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_28_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_27_mults_d8(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_25_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_24_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_22_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t23,t24,t25,t26,t27,t28,t29,t30,t31,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
    _ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t23,t24,t25,t26,t27,t28,t29,t30,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
    _ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t23,t24,t25,t26,t27,t28,t29,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
    _ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t23,t24,t25,t26,t27,t28,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
    _ff_sum_7_mults_d2(t4,f[0],f[1],t23,t24,t25,t26,t27,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
    _ff_sum_6_mults_d1(t3,f[0],f[1],t23,t24,t25,t26,g[0],g[1],g[2],g[3],f[1],f[2]);
    _ff_sum_4_mults_d1(t2,f[0],t23,t24,t25,g[0],g[1],g[2],f[1]);
    _ff_sum_3_mults(t1,f[0],t23,t24,g[0],g[1],f[0]);
    _ff_mult(o[0],t23,g[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22);
}

static inline void ff_poly_square_mod_29 (ff_t o[29], ff_t f[29], ff_t g[28])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56;

    _ff_square(t56,f[28]);
    _ff_dbl_mult(t55,f[27],f[28]);
    _ff_sum_3_mults_d1(t54,f[26],f[27],t56,g[27],f[27],f[28]);
    _ff_sum_4_mults_d2(t53,f[25],f[26],t55,t56,g[26],g[27],f[27],f[28]);
    _ff_sum_6_mults_d2(t52,f[24],f[25],f[26],t54,t55,t56,g[25],g[26],g[27],f[26],f[27],f[28]);
    _ff_sum_7_mults_d3(t51,f[23],f[24],f[25],t53,t54,t55,t56,g[24],g[25],g[26],g[27],f[26],f[27],f[28]);
    _ff_sum_9_mults_d3(t50,f[22],f[23],f[24],f[25],t52,t53,t54,t55,t56,g[23],g[24],g[25],g[26],g[27],f[25],f[26],f[27],f[28]);
    _ff_sum_10_mults_d4(t49,f[21],f[22],f[23],f[24],t51,t52,t53,t54,t55,t56,g[22],g[23],g[24],g[25],g[26],g[27],f[25],f[26],f[27],f[28]);
    _ff_sum_12_mults_d4(t48,f[20],f[21],f[22],f[23],f[24],t50,t51,t52,t53,t54,t55,t56,g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_13_mults_d5(t47,f[19],f[20],f[21],f[22],f[23],t49,t50,t51,t52,t53,t54,t55,t56,g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_15_mults_d5(t46,f[18],f[19],f[20],f[21],f[22],f[23],t48,t49,t50,t51,t52,t53,t54,t55,t56,g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_16_mults_d6(t45,f[17],f[18],f[19],f[20],f[21],f[22],t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_18_mults_d6(t44,f[16],f[17],f[18],f[19],f[20],f[21],f[22],t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_19_mults_d7(t43,f[15],f[16],f[17],f[18],f[19],f[20],f[21],t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_21_mults_d7(t42,f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_22_mults_d8(t41,f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_24_mults_d8(t40,f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_25_mults_d9(t39,f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_27_mults_d9(t38,f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_28_mults_d10(t37,f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_30_mults_d10(t36,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_31_mults_d11(t35,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_33_mults_d11(t34,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_34_mults_d12(t33,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_36_mults_d12(t32,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_37_mults_d13(t31,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_39_mults_d13(t30,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_40_mults_d14(t29,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);

    _ff_sum_42_mults_d14(t28,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_42_mults_d14(t27,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27]);
    _ff_sum_41_mults_d13(t26,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26]);
    _ff_sum_39_mults_d13(t25,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25]);
    _ff_sum_38_mults_d12(t24,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24]);
    _ff_sum_36_mults_d12(t23,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23]);
    _ff_sum_35_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_33_mults_d11(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_32_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_30_mults_d10(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_29_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_27_mults_d9(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_26_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_24_mults_d8(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_23_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t29,t30,t31,t32,t33,t34,t35,t36,t37,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t29,t30,t31,t32,t33,t34,t35,t36,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t29,t30,t31,t32,t33,t34,t35,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
    _ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t29,t30,t31,t32,t33,t34,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
    _ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t29,t30,t31,t32,t33,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
    _ff_sum_6_mults_d2(t3,f[0],f[1],t29,t30,t31,t32,g[0],g[1],g[2],g[3],f[2],f[3]);
    _ff_sum_5_mults_d1(t2,f[0],f[1],t29,t30,t31,g[0],g[1],g[2],f[1],f[2]);
    _ff_sum_3_mults_d1(t1,f[0],t29,t30,g[0],g[1],f[1]);
    _ff_sum_2_mults(o[0],f[0],t29,g[0],f[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22); _ff_set(o[23],t23); _ff_set(o[24],t24); _ff_set(o[25],t25); _ff_set(o[26],t26); _ff_set(o[27],t27); _ff_set(o[28],t28);
}

static inline void ff_poly_square_mult_x_mod_29 (ff_t o[29], ff_t f[29], ff_t g[28])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57;

    _ff_square(t57,f[28]);
    _ff_dbl_mult(t56,f[27],f[28]);
    _ff_sum_3_mults_d1(t55,f[26],f[27],t57,g[27],f[27],f[28]);
    _ff_sum_4_mults_d2(t54,f[25],f[26],t56,t57,g[26],g[27],f[27],f[28]);
    _ff_sum_6_mults_d2(t53,f[24],f[25],f[26],t55,t56,t57,g[25],g[26],g[27],f[26],f[27],f[28]);
    _ff_sum_7_mults_d3(t52,f[23],f[24],f[25],t54,t55,t56,t57,g[24],g[25],g[26],g[27],f[26],f[27],f[28]);
    _ff_sum_9_mults_d3(t51,f[22],f[23],f[24],f[25],t53,t54,t55,t56,t57,g[23],g[24],g[25],g[26],g[27],f[25],f[26],f[27],f[28]);
    _ff_sum_10_mults_d4(t50,f[21],f[22],f[23],f[24],t52,t53,t54,t55,t56,t57,g[22],g[23],g[24],g[25],g[26],g[27],f[25],f[26],f[27],f[28]);
    _ff_sum_12_mults_d4(t49,f[20],f[21],f[22],f[23],f[24],t51,t52,t53,t54,t55,t56,t57,g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_13_mults_d5(t48,f[19],f[20],f[21],f[22],f[23],t50,t51,t52,t53,t54,t55,t56,t57,g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_15_mults_d5(t47,f[18],f[19],f[20],f[21],f[22],f[23],t49,t50,t51,t52,t53,t54,t55,t56,t57,g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_16_mults_d6(t46,f[17],f[18],f[19],f[20],f[21],f[22],t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_18_mults_d6(t45,f[16],f[17],f[18],f[19],f[20],f[21],f[22],t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_19_mults_d7(t44,f[15],f[16],f[17],f[18],f[19],f[20],f[21],t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_21_mults_d7(t43,f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_22_mults_d8(t42,f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_24_mults_d8(t41,f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_25_mults_d9(t40,f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_27_mults_d9(t39,f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_28_mults_d10(t38,f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_30_mults_d10(t37,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_31_mults_d11(t36,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_33_mults_d11(t35,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_34_mults_d12(t34,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_36_mults_d12(t33,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_37_mults_d13(t32,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_39_mults_d13(t31,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_40_mults_d14(t30,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_42_mults_d14(t29,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);

    _ff_sum_42_mults_d14(t28,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27]);
    _ff_sum_42_mults_d13(t27,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26]);
    _ff_sum_40_mults_d13(t26,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25]);
    _ff_sum_39_mults_d12(t25,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24]);
    _ff_sum_37_mults_d12(t24,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23]);
    _ff_sum_36_mults_d11(t23,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_34_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_33_mults_d10(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_31_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_30_mults_d9(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_28_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_27_mults_d8(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_25_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_24_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_22_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t29,t30,t31,t32,t33,t34,t35,t36,t37,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
    _ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t29,t30,t31,t32,t33,t34,t35,t36,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
    _ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t29,t30,t31,t32,t33,t34,t35,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
    _ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t29,t30,t31,t32,t33,t34,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
    _ff_sum_7_mults_d2(t4,f[0],f[1],t29,t30,t31,t32,t33,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
    _ff_sum_6_mults_d1(t3,f[0],f[1],t29,t30,t31,t32,g[0],g[1],g[2],g[3],f[1],f[2]);
    _ff_sum_4_mults_d1(t2,f[0],t29,t30,t31,g[0],g[1],g[2],f[1]);
    _ff_sum_3_mults(t1,f[0],t29,t30,g[0],g[1],f[0]);
    _ff_mult(o[0],t29,g[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22); _ff_set(o[23],t23); _ff_set(o[24],t24); _ff_set(o[25],t25); _ff_set(o[26],t26); _ff_set(o[27],t27); _ff_set(o[28],t28);
}

static inline void ff_poly_square_mod_31 (ff_t o[31], ff_t f[31], ff_t g[30])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60;

    _ff_square(t60,f[30]);
    _ff_dbl_mult(t59,f[29],f[30]);
    _ff_sum_3_mults_d1(t58,f[28],f[29],t60,g[29],f[29],f[30]);
    _ff_sum_4_mults_d2(t57,f[27],f[28],t59,t60,g[28],g[29],f[29],f[30]);
    _ff_sum_6_mults_d2(t56,f[26],f[27],f[28],t58,t59,t60,g[27],g[28],g[29],f[28],f[29],f[30]);
    _ff_sum_7_mults_d3(t55,f[25],f[26],f[27],t57,t58,t59,t60,g[26],g[27],g[28],g[29],f[28],f[29],f[30]);
    _ff_sum_9_mults_d3(t54,f[24],f[25],f[26],f[27],t56,t57,t58,t59,t60,g[25],g[26],g[27],g[28],g[29],f[27],f[28],f[29],f[30]);
    _ff_sum_10_mults_d4(t53,f[23],f[24],f[25],f[26],t55,t56,t57,t58,t59,t60,g[24],g[25],g[26],g[27],g[28],g[29],f[27],f[28],f[29],f[30]);
    _ff_sum_12_mults_d4(t52,f[22],f[23],f[24],f[25],f[26],t54,t55,t56,t57,t58,t59,t60,g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_13_mults_d5(t51,f[21],f[22],f[23],f[24],f[25],t53,t54,t55,t56,t57,t58,t59,t60,g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_15_mults_d5(t50,f[20],f[21],f[22],f[23],f[24],f[25],t52,t53,t54,t55,t56,t57,t58,t59,t60,g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_16_mults_d6(t49,f[19],f[20],f[21],f[22],f[23],f[24],t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_18_mults_d6(t48,f[18],f[19],f[20],f[21],f[22],f[23],f[24],t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_19_mults_d7(t47,f[17],f[18],f[19],f[20],f[21],f[22],f[23],t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_21_mults_d7(t46,f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_22_mults_d8(t45,f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_24_mults_d8(t44,f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_25_mults_d9(t43,f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_27_mults_d9(t42,f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_28_mults_d10(t41,f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_30_mults_d10(t40,f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_31_mults_d11(t39,f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_33_mults_d11(t38,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_34_mults_d12(t37,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_36_mults_d12(t36,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_37_mults_d13(t35,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_39_mults_d13(t34,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_40_mults_d14(t33,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_42_mults_d14(t32,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_43_mults_d15(t31,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);

    _ff_sum_45_mults_d15(t30,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_45_mults_d15(t29,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29]);
    _ff_sum_44_mults_d14(t28,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_42_mults_d14(t27,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27]);
    _ff_sum_41_mults_d13(t26,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26]);
    _ff_sum_39_mults_d13(t25,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25]);
    _ff_sum_38_mults_d12(t24,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24]);
    _ff_sum_36_mults_d12(t23,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23]);
    _ff_sum_35_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_33_mults_d11(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_32_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_30_mults_d10(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_29_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_27_mults_d9(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_26_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_24_mults_d8(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_23_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t31,t32,t33,t34,t35,t36,t37,t38,t39,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t31,t32,t33,t34,t35,t36,t37,t38,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t31,t32,t33,t34,t35,t36,t37,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
    _ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t31,t32,t33,t34,t35,t36,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
    _ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t31,t32,t33,t34,t35,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
    _ff_sum_6_mults_d2(t3,f[0],f[1],t31,t32,t33,t34,g[0],g[1],g[2],g[3],f[2],f[3]);
    _ff_sum_5_mults_d1(t2,f[0],f[1],t31,t32,t33,g[0],g[1],g[2],f[1],f[2]);
    _ff_sum_3_mults_d1(t1,f[0],t31,t32,g[0],g[1],f[1]);
    _ff_sum_2_mults(o[0],f[0],t31,g[0],f[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22); _ff_set(o[23],t23); _ff_set(o[24],t24); _ff_set(o[25],t25); _ff_set(o[26],t26); _ff_set(o[27],t27); _ff_set(o[28],t28); _ff_set(o[29],t29); _ff_set(o[30],t30);
}

static inline void ff_poly_square_mult_x_mod_31 (ff_t o[31], ff_t f[31], ff_t g[30])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61;

    _ff_square(t61,f[30]);
    _ff_dbl_mult(t60,f[29],f[30]);
    _ff_sum_3_mults_d1(t59,f[28],f[29],t61,g[29],f[29],f[30]);
    _ff_sum_4_mults_d2(t58,f[27],f[28],t60,t61,g[28],g[29],f[29],f[30]);
    _ff_sum_6_mults_d2(t57,f[26],f[27],f[28],t59,t60,t61,g[27],g[28],g[29],f[28],f[29],f[30]);
    _ff_sum_7_mults_d3(t56,f[25],f[26],f[27],t58,t59,t60,t61,g[26],g[27],g[28],g[29],f[28],f[29],f[30]);
    _ff_sum_9_mults_d3(t55,f[24],f[25],f[26],f[27],t57,t58,t59,t60,t61,g[25],g[26],g[27],g[28],g[29],f[27],f[28],f[29],f[30]);
    _ff_sum_10_mults_d4(t54,f[23],f[24],f[25],f[26],t56,t57,t58,t59,t60,t61,g[24],g[25],g[26],g[27],g[28],g[29],f[27],f[28],f[29],f[30]);
    _ff_sum_12_mults_d4(t53,f[22],f[23],f[24],f[25],f[26],t55,t56,t57,t58,t59,t60,t61,g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_13_mults_d5(t52,f[21],f[22],f[23],f[24],f[25],t54,t55,t56,t57,t58,t59,t60,t61,g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_15_mults_d5(t51,f[20],f[21],f[22],f[23],f[24],f[25],t53,t54,t55,t56,t57,t58,t59,t60,t61,g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_16_mults_d6(t50,f[19],f[20],f[21],f[22],f[23],f[24],t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_18_mults_d6(t49,f[18],f[19],f[20],f[21],f[22],f[23],f[24],t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_19_mults_d7(t48,f[17],f[18],f[19],f[20],f[21],f[22],f[23],t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_21_mults_d7(t47,f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_22_mults_d8(t46,f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_24_mults_d8(t45,f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_25_mults_d9(t44,f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_27_mults_d9(t43,f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_28_mults_d10(t42,f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_30_mults_d10(t41,f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_31_mults_d11(t40,f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_33_mults_d11(t39,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_34_mults_d12(t38,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_36_mults_d12(t37,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_37_mults_d13(t36,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_39_mults_d13(t35,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_40_mults_d14(t34,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_42_mults_d14(t33,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_43_mults_d15(t32,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_45_mults_d15(t31,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);

    _ff_sum_45_mults_d15(t30,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29]);
    _ff_sum_45_mults_d14(t29,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_43_mults_d14(t28,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27]);
    _ff_sum_42_mults_d13(t27,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26]);
    _ff_sum_40_mults_d13(t26,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25]);
    _ff_sum_39_mults_d12(t25,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24]);
    _ff_sum_37_mults_d12(t24,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23]);
    _ff_sum_36_mults_d11(t23,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_34_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_33_mults_d10(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_31_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_30_mults_d9(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_28_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_27_mults_d8(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_25_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_24_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_22_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t31,t32,t33,t34,t35,t36,t37,t38,t39,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
    _ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t31,t32,t33,t34,t35,t36,t37,t38,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
    _ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t31,t32,t33,t34,t35,t36,t37,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
    _ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t31,t32,t33,t34,t35,t36,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
    _ff_sum_7_mults_d2(t4,f[0],f[1],t31,t32,t33,t34,t35,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
    _ff_sum_6_mults_d1(t3,f[0],f[1],t31,t32,t33,t34,g[0],g[1],g[2],g[3],f[1],f[2]);
    _ff_sum_4_mults_d1(t2,f[0],t31,t32,t33,g[0],g[1],g[2],f[1]);
    _ff_sum_3_mults(t1,f[0],t31,t32,g[0],g[1],f[0]);
    _ff_mult(o[0],t31,g[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22); _ff_set(o[23],t23); _ff_set(o[24],t24); _ff_set(o[25],t25); _ff_set(o[26],t26); _ff_set(o[27],t27); _ff_set(o[28],t28); _ff_set(o[29],t29); _ff_set(o[30],t30);
}

static inline void ff_poly_square_mod_37 (ff_t o[37], ff_t f[37], ff_t g[36])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72;

    _ff_square(t72,f[36]);
    _ff_dbl_mult(t71,f[35],f[36]);
    _ff_sum_3_mults_d1(t70,f[34],f[35],t72,g[35],f[35],f[36]);
    _ff_sum_4_mults_d2(t69,f[33],f[34],t71,t72,g[34],g[35],f[35],f[36]);
    _ff_sum_6_mults_d2(t68,f[32],f[33],f[34],t70,t71,t72,g[33],g[34],g[35],f[34],f[35],f[36]);
    _ff_sum_7_mults_d3(t67,f[31],f[32],f[33],t69,t70,t71,t72,g[32],g[33],g[34],g[35],f[34],f[35],f[36]);
    _ff_sum_9_mults_d3(t66,f[30],f[31],f[32],f[33],t68,t69,t70,t71,t72,g[31],g[32],g[33],g[34],g[35],f[33],f[34],f[35],f[36]);
    _ff_sum_10_mults_d4(t65,f[29],f[30],f[31],f[32],t67,t68,t69,t70,t71,t72,g[30],g[31],g[32],g[33],g[34],g[35],f[33],f[34],f[35],f[36]);
    _ff_sum_12_mults_d4(t64,f[28],f[29],f[30],f[31],f[32],t66,t67,t68,t69,t70,t71,t72,g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_13_mults_d5(t63,f[27],f[28],f[29],f[30],f[31],t65,t66,t67,t68,t69,t70,t71,t72,g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_15_mults_d5(t62,f[26],f[27],f[28],f[29],f[30],f[31],t64,t65,t66,t67,t68,t69,t70,t71,t72,g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_16_mults_d6(t61,f[25],f[26],f[27],f[28],f[29],f[30],t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_18_mults_d6(t60,f[24],f[25],f[26],f[27],f[28],f[29],f[30],t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_19_mults_d7(t59,f[23],f[24],f[25],f[26],f[27],f[28],f[29],t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_21_mults_d7(t58,f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_22_mults_d8(t57,f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_24_mults_d8(t56,f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_25_mults_d9(t55,f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_27_mults_d9(t54,f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_28_mults_d10(t53,f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_30_mults_d10(t52,f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_31_mults_d11(t51,f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_33_mults_d11(t50,f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_34_mults_d12(t49,f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_36_mults_d12(t48,f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_37_mults_d13(t47,f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_39_mults_d13(t46,f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_40_mults_d14(t45,f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_42_mults_d14(t44,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_43_mults_d15(t43,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_45_mults_d15(t42,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_46_mults_d16(t41,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_48_mults_d16(t40,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_49_mults_d17(t39,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_51_mults_d17(t38,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_52_mults_d18(t37,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);

    _ff_sum_54_mults_d18(t36,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_54_mults_d18(t35,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35]);
    _ff_sum_53_mults_d17(t34,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34]);
    _ff_sum_51_mults_d17(t33,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33]);
    _ff_sum_50_mults_d16(t32,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32]);
    _ff_sum_48_mults_d16(t31,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31]);
    _ff_sum_47_mults_d15(t30,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_45_mults_d15(t29,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29]);
    _ff_sum_44_mults_d14(t28,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_42_mults_d14(t27,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27]);
    _ff_sum_41_mults_d13(t26,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26]);
    _ff_sum_39_mults_d13(t25,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25]);
    _ff_sum_38_mults_d12(t24,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24]);
    _ff_sum_36_mults_d12(t23,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23]);
    _ff_sum_35_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_33_mults_d11(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_32_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_30_mults_d10(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_29_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_27_mults_d9(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_26_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_24_mults_d8(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_23_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t37,t38,t39,t40,t41,t42,t43,t44,t45,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t37,t38,t39,t40,t41,t42,t43,t44,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t37,t38,t39,t40,t41,t42,t43,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
    _ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t37,t38,t39,t40,t41,t42,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
    _ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t37,t38,t39,t40,t41,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
    _ff_sum_6_mults_d2(t3,f[0],f[1],t37,t38,t39,t40,g[0],g[1],g[2],g[3],f[2],f[3]);
    _ff_sum_5_mults_d1(t2,f[0],f[1],t37,t38,t39,g[0],g[1],g[2],f[1],f[2]);
    _ff_sum_3_mults_d1(t1,f[0],t37,t38,g[0],g[1],f[1]);
    _ff_sum_2_mults(o[0],f[0],t37,g[0],f[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22); _ff_set(o[23],t23); _ff_set(o[24],t24); _ff_set(o[25],t25); _ff_set(o[26],t26); _ff_set(o[27],t27); _ff_set(o[28],t28); _ff_set(o[29],t29); _ff_set(o[30],t30); _ff_set(o[31],t31); _ff_set(o[32],t32); _ff_set(o[33],t33); _ff_set(o[34],t34); _ff_set(o[35],t35); _ff_set(o[36],t36);
}

static inline void ff_poly_square_mult_x_mod_37 (ff_t o[37], ff_t f[37], ff_t g[36])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73;

    _ff_square(t73,f[36]);
    _ff_dbl_mult(t72,f[35],f[36]);
    _ff_sum_3_mults_d1(t71,f[34],f[35],t73,g[35],f[35],f[36]);
    _ff_sum_4_mults_d2(t70,f[33],f[34],t72,t73,g[34],g[35],f[35],f[36]);
    _ff_sum_6_mults_d2(t69,f[32],f[33],f[34],t71,t72,t73,g[33],g[34],g[35],f[34],f[35],f[36]);
    _ff_sum_7_mults_d3(t68,f[31],f[32],f[33],t70,t71,t72,t73,g[32],g[33],g[34],g[35],f[34],f[35],f[36]);
    _ff_sum_9_mults_d3(t67,f[30],f[31],f[32],f[33],t69,t70,t71,t72,t73,g[31],g[32],g[33],g[34],g[35],f[33],f[34],f[35],f[36]);
    _ff_sum_10_mults_d4(t66,f[29],f[30],f[31],f[32],t68,t69,t70,t71,t72,t73,g[30],g[31],g[32],g[33],g[34],g[35],f[33],f[34],f[35],f[36]);
    _ff_sum_12_mults_d4(t65,f[28],f[29],f[30],f[31],f[32],t67,t68,t69,t70,t71,t72,t73,g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_13_mults_d5(t64,f[27],f[28],f[29],f[30],f[31],t66,t67,t68,t69,t70,t71,t72,t73,g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_15_mults_d5(t63,f[26],f[27],f[28],f[29],f[30],f[31],t65,t66,t67,t68,t69,t70,t71,t72,t73,g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_16_mults_d6(t62,f[25],f[26],f[27],f[28],f[29],f[30],t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_18_mults_d6(t61,f[24],f[25],f[26],f[27],f[28],f[29],f[30],t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_19_mults_d7(t60,f[23],f[24],f[25],f[26],f[27],f[28],f[29],t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_21_mults_d7(t59,f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_22_mults_d8(t58,f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_24_mults_d8(t57,f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_25_mults_d9(t56,f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_27_mults_d9(t55,f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_28_mults_d10(t54,f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_30_mults_d10(t53,f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_31_mults_d11(t52,f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_33_mults_d11(t51,f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_34_mults_d12(t50,f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_36_mults_d12(t49,f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_37_mults_d13(t48,f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_39_mults_d13(t47,f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_40_mults_d14(t46,f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_42_mults_d14(t45,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_43_mults_d15(t44,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_45_mults_d15(t43,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_46_mults_d16(t42,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_48_mults_d16(t41,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_49_mults_d17(t40,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_51_mults_d17(t39,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_52_mults_d18(t38,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_54_mults_d18(t37,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);

    _ff_sum_54_mults_d18(t36,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35]);
    _ff_sum_54_mults_d17(t35,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34]);
    _ff_sum_52_mults_d17(t34,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33]);
    _ff_sum_51_mults_d16(t33,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32]);
    _ff_sum_49_mults_d16(t32,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31]);
    _ff_sum_48_mults_d15(t31,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_46_mults_d15(t30,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29]);
    _ff_sum_45_mults_d14(t29,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_43_mults_d14(t28,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27]);
    _ff_sum_42_mults_d13(t27,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26]);
    _ff_sum_40_mults_d13(t26,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25]);
    _ff_sum_39_mults_d12(t25,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24]);
    _ff_sum_37_mults_d12(t24,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23]);
    _ff_sum_36_mults_d11(t23,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_34_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_33_mults_d10(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_31_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_30_mults_d9(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_28_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_27_mults_d8(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_25_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_24_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_22_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t37,t38,t39,t40,t41,t42,t43,t44,t45,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
    _ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t37,t38,t39,t40,t41,t42,t43,t44,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
    _ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t37,t38,t39,t40,t41,t42,t43,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
    _ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t37,t38,t39,t40,t41,t42,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
    _ff_sum_7_mults_d2(t4,f[0],f[1],t37,t38,t39,t40,t41,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
    _ff_sum_6_mults_d1(t3,f[0],f[1],t37,t38,t39,t40,g[0],g[1],g[2],g[3],f[1],f[2]);
    _ff_sum_4_mults_d1(t2,f[0],t37,t38,t39,g[0],g[1],g[2],f[1]);
    _ff_sum_3_mults(t1,f[0],t37,t38,g[0],g[1],f[0]);
    _ff_mult(o[0],t37,g[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22); _ff_set(o[23],t23); _ff_set(o[24],t24); _ff_set(o[25],t25); _ff_set(o[26],t26); _ff_set(o[27],t27); _ff_set(o[28],t28); _ff_set(o[29],t29); _ff_set(o[30],t30); _ff_set(o[31],t31); _ff_set(o[32],t32); _ff_set(o[33],t33); _ff_set(o[34],t34); _ff_set(o[35],t35); _ff_set(o[36],t36);
}

static inline void ff_poly_square_mod_41 (ff_t o[41], ff_t f[41], ff_t g[40])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80;

    _ff_square(t80,f[40]);
    _ff_dbl_mult(t79,f[39],f[40]);
    _ff_sum_3_mults_d1(t78,f[38],f[39],t80,g[39],f[39],f[40]);
    _ff_sum_4_mults_d2(t77,f[37],f[38],t79,t80,g[38],g[39],f[39],f[40]);
    _ff_sum_6_mults_d2(t76,f[36],f[37],f[38],t78,t79,t80,g[37],g[38],g[39],f[38],f[39],f[40]);
    _ff_sum_7_mults_d3(t75,f[35],f[36],f[37],t77,t78,t79,t80,g[36],g[37],g[38],g[39],f[38],f[39],f[40]);
    _ff_sum_9_mults_d3(t74,f[34],f[35],f[36],f[37],t76,t77,t78,t79,t80,g[35],g[36],g[37],g[38],g[39],f[37],f[38],f[39],f[40]);
    _ff_sum_10_mults_d4(t73,f[33],f[34],f[35],f[36],t75,t76,t77,t78,t79,t80,g[34],g[35],g[36],g[37],g[38],g[39],f[37],f[38],f[39],f[40]);
    _ff_sum_12_mults_d4(t72,f[32],f[33],f[34],f[35],f[36],t74,t75,t76,t77,t78,t79,t80,g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_13_mults_d5(t71,f[31],f[32],f[33],f[34],f[35],t73,t74,t75,t76,t77,t78,t79,t80,g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_15_mults_d5(t70,f[30],f[31],f[32],f[33],f[34],f[35],t72,t73,t74,t75,t76,t77,t78,t79,t80,g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_16_mults_d6(t69,f[29],f[30],f[31],f[32],f[33],f[34],t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_18_mults_d6(t68,f[28],f[29],f[30],f[31],f[32],f[33],f[34],t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_19_mults_d7(t67,f[27],f[28],f[29],f[30],f[31],f[32],f[33],t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_21_mults_d7(t66,f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_22_mults_d8(t65,f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_24_mults_d8(t64,f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_25_mults_d9(t63,f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_27_mults_d9(t62,f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_28_mults_d10(t61,f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_30_mults_d10(t60,f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_31_mults_d11(t59,f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_33_mults_d11(t58,f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_34_mults_d12(t57,f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_36_mults_d12(t56,f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_37_mults_d13(t55,f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_39_mults_d13(t54,f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_40_mults_d14(t53,f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_42_mults_d14(t52,f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_43_mults_d15(t51,f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_45_mults_d15(t50,f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_46_mults_d16(t49,f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_48_mults_d16(t48,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_49_mults_d17(t47,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_51_mults_d17(t46,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_52_mults_d18(t45,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_54_mults_d18(t44,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_55_mults_d19(t43,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_57_mults_d19(t42,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_58_mults_d20(t41,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);

    _ff_sum_60_mults_d20(t40,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_60_mults_d20(t39,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39]);
    _ff_sum_59_mults_d19(t38,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38]);
    _ff_sum_57_mults_d19(t37,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37]);
    _ff_sum_56_mults_d18(t36,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_54_mults_d18(t35,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35]);
    _ff_sum_53_mults_d17(t34,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34]);
    _ff_sum_51_mults_d17(t33,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33]);
    _ff_sum_50_mults_d16(t32,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32]);
    _ff_sum_48_mults_d16(t31,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31]);
    _ff_sum_47_mults_d15(t30,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_45_mults_d15(t29,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29]);
    _ff_sum_44_mults_d14(t28,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_42_mults_d14(t27,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27]);
    _ff_sum_41_mults_d13(t26,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26]);
    _ff_sum_39_mults_d13(t25,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25]);
    _ff_sum_38_mults_d12(t24,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24]);
    _ff_sum_36_mults_d12(t23,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23]);
    _ff_sum_35_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_33_mults_d11(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_32_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_30_mults_d10(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_29_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_27_mults_d9(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_26_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_24_mults_d8(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_23_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t41,t42,t43,t44,t45,t46,t47,t48,t49,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t41,t42,t43,t44,t45,t46,t47,t48,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t41,t42,t43,t44,t45,t46,t47,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
    _ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t41,t42,t43,t44,t45,t46,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
    _ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t41,t42,t43,t44,t45,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
    _ff_sum_6_mults_d2(t3,f[0],f[1],t41,t42,t43,t44,g[0],g[1],g[2],g[3],f[2],f[3]);
    _ff_sum_5_mults_d1(t2,f[0],f[1],t41,t42,t43,g[0],g[1],g[2],f[1],f[2]);
    _ff_sum_3_mults_d1(t1,f[0],t41,t42,g[0],g[1],f[1]);
    _ff_sum_2_mults(o[0],f[0],t41,g[0],f[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22); _ff_set(o[23],t23); _ff_set(o[24],t24); _ff_set(o[25],t25); _ff_set(o[26],t26); _ff_set(o[27],t27); _ff_set(o[28],t28); _ff_set(o[29],t29); _ff_set(o[30],t30); _ff_set(o[31],t31); _ff_set(o[32],t32); _ff_set(o[33],t33); _ff_set(o[34],t34); _ff_set(o[35],t35); _ff_set(o[36],t36); _ff_set(o[37],t37); _ff_set(o[38],t38); _ff_set(o[39],t39); _ff_set(o[40],t40);
}

static inline void ff_poly_square_mult_x_mod_41 (ff_t o[41], ff_t f[41], ff_t g[40])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81;

    _ff_square(t81,f[40]);
    _ff_dbl_mult(t80,f[39],f[40]);
    _ff_sum_3_mults_d1(t79,f[38],f[39],t81,g[39],f[39],f[40]);
    _ff_sum_4_mults_d2(t78,f[37],f[38],t80,t81,g[38],g[39],f[39],f[40]);
    _ff_sum_6_mults_d2(t77,f[36],f[37],f[38],t79,t80,t81,g[37],g[38],g[39],f[38],f[39],f[40]);
    _ff_sum_7_mults_d3(t76,f[35],f[36],f[37],t78,t79,t80,t81,g[36],g[37],g[38],g[39],f[38],f[39],f[40]);
    _ff_sum_9_mults_d3(t75,f[34],f[35],f[36],f[37],t77,t78,t79,t80,t81,g[35],g[36],g[37],g[38],g[39],f[37],f[38],f[39],f[40]);
    _ff_sum_10_mults_d4(t74,f[33],f[34],f[35],f[36],t76,t77,t78,t79,t80,t81,g[34],g[35],g[36],g[37],g[38],g[39],f[37],f[38],f[39],f[40]);
    _ff_sum_12_mults_d4(t73,f[32],f[33],f[34],f[35],f[36],t75,t76,t77,t78,t79,t80,t81,g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_13_mults_d5(t72,f[31],f[32],f[33],f[34],f[35],t74,t75,t76,t77,t78,t79,t80,t81,g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_15_mults_d5(t71,f[30],f[31],f[32],f[33],f[34],f[35],t73,t74,t75,t76,t77,t78,t79,t80,t81,g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_16_mults_d6(t70,f[29],f[30],f[31],f[32],f[33],f[34],t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_18_mults_d6(t69,f[28],f[29],f[30],f[31],f[32],f[33],f[34],t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_19_mults_d7(t68,f[27],f[28],f[29],f[30],f[31],f[32],f[33],t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_21_mults_d7(t67,f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_22_mults_d8(t66,f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_24_mults_d8(t65,f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_25_mults_d9(t64,f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_27_mults_d9(t63,f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_28_mults_d10(t62,f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_30_mults_d10(t61,f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_31_mults_d11(t60,f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_33_mults_d11(t59,f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_34_mults_d12(t58,f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_36_mults_d12(t57,f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_37_mults_d13(t56,f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_39_mults_d13(t55,f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_40_mults_d14(t54,f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_42_mults_d14(t53,f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_43_mults_d15(t52,f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_45_mults_d15(t51,f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_46_mults_d16(t50,f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_48_mults_d16(t49,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_49_mults_d17(t48,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_51_mults_d17(t47,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_52_mults_d18(t46,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_54_mults_d18(t45,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_55_mults_d19(t44,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_57_mults_d19(t43,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_58_mults_d20(t42,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_60_mults_d20(t41,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);

    _ff_sum_60_mults_d20(t40,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39]);
    _ff_sum_60_mults_d19(t39,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38]);
    _ff_sum_58_mults_d19(t38,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37]);
    _ff_sum_57_mults_d18(t37,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_55_mults_d18(t36,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35]);
    _ff_sum_54_mults_d17(t35,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34]);
    _ff_sum_52_mults_d17(t34,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33]);
    _ff_sum_51_mults_d16(t33,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32]);
    _ff_sum_49_mults_d16(t32,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31]);
    _ff_sum_48_mults_d15(t31,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_46_mults_d15(t30,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29]);
    _ff_sum_45_mults_d14(t29,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_43_mults_d14(t28,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27]);
    _ff_sum_42_mults_d13(t27,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26]);
    _ff_sum_40_mults_d13(t26,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25]);
    _ff_sum_39_mults_d12(t25,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24]);
    _ff_sum_37_mults_d12(t24,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23]);
    _ff_sum_36_mults_d11(t23,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_34_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_33_mults_d10(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_31_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_30_mults_d9(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_28_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_27_mults_d8(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_25_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_24_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_22_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t41,t42,t43,t44,t45,t46,t47,t48,t49,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
    _ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t41,t42,t43,t44,t45,t46,t47,t48,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
    _ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t41,t42,t43,t44,t45,t46,t47,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
    _ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t41,t42,t43,t44,t45,t46,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
    _ff_sum_7_mults_d2(t4,f[0],f[1],t41,t42,t43,t44,t45,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
    _ff_sum_6_mults_d1(t3,f[0],f[1],t41,t42,t43,t44,g[0],g[1],g[2],g[3],f[1],f[2]);
    _ff_sum_4_mults_d1(t2,f[0],t41,t42,t43,g[0],g[1],g[2],f[1]);
    _ff_sum_3_mults(t1,f[0],t41,t42,g[0],g[1],f[0]);
    _ff_mult(o[0],t41,g[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22); _ff_set(o[23],t23); _ff_set(o[24],t24); _ff_set(o[25],t25); _ff_set(o[26],t26); _ff_set(o[27],t27); _ff_set(o[28],t28); _ff_set(o[29],t29); _ff_set(o[30],t30); _ff_set(o[31],t31); _ff_set(o[32],t32); _ff_set(o[33],t33); _ff_set(o[34],t34); _ff_set(o[35],t35); _ff_set(o[36],t36); _ff_set(o[37],t37); _ff_set(o[38],t38); _ff_set(o[39],t39); _ff_set(o[40],t40);
}

static inline void ff_poly_square_mod_43 (ff_t o[43], ff_t f[43], ff_t g[42])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84;

    _ff_square(t84,f[42]);
    _ff_dbl_mult(t83,f[41],f[42]);
    _ff_sum_3_mults_d1(t82,f[40],f[41],t84,g[41],f[41],f[42]);
    _ff_sum_4_mults_d2(t81,f[39],f[40],t83,t84,g[40],g[41],f[41],f[42]);
    _ff_sum_6_mults_d2(t80,f[38],f[39],f[40],t82,t83,t84,g[39],g[40],g[41],f[40],f[41],f[42]);
    _ff_sum_7_mults_d3(t79,f[37],f[38],f[39],t81,t82,t83,t84,g[38],g[39],g[40],g[41],f[40],f[41],f[42]);
    _ff_sum_9_mults_d3(t78,f[36],f[37],f[38],f[39],t80,t81,t82,t83,t84,g[37],g[38],g[39],g[40],g[41],f[39],f[40],f[41],f[42]);
    _ff_sum_10_mults_d4(t77,f[35],f[36],f[37],f[38],t79,t80,t81,t82,t83,t84,g[36],g[37],g[38],g[39],g[40],g[41],f[39],f[40],f[41],f[42]);
    _ff_sum_12_mults_d4(t76,f[34],f[35],f[36],f[37],f[38],t78,t79,t80,t81,t82,t83,t84,g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_13_mults_d5(t75,f[33],f[34],f[35],f[36],f[37],t77,t78,t79,t80,t81,t82,t83,t84,g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_15_mults_d5(t74,f[32],f[33],f[34],f[35],f[36],f[37],t76,t77,t78,t79,t80,t81,t82,t83,t84,g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_16_mults_d6(t73,f[31],f[32],f[33],f[34],f[35],f[36],t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_18_mults_d6(t72,f[30],f[31],f[32],f[33],f[34],f[35],f[36],t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_19_mults_d7(t71,f[29],f[30],f[31],f[32],f[33],f[34],f[35],t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_21_mults_d7(t70,f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_22_mults_d8(t69,f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_24_mults_d8(t68,f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_25_mults_d9(t67,f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_27_mults_d9(t66,f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_28_mults_d10(t65,f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_30_mults_d10(t64,f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_31_mults_d11(t63,f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_33_mults_d11(t62,f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_34_mults_d12(t61,f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_36_mults_d12(t60,f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_37_mults_d13(t59,f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_39_mults_d13(t58,f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_40_mults_d14(t57,f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_42_mults_d14(t56,f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_43_mults_d15(t55,f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_45_mults_d15(t54,f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_46_mults_d16(t53,f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_48_mults_d16(t52,f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_49_mults_d17(t51,f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_51_mults_d17(t50,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_52_mults_d18(t49,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_54_mults_d18(t48,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_55_mults_d19(t47,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_57_mults_d19(t46,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_58_mults_d20(t45,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_60_mults_d20(t44,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_61_mults_d21(t43,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);

    _ff_sum_63_mults_d21(t42,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_63_mults_d21(t41,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41]);
    _ff_sum_62_mults_d20(t40,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_60_mults_d20(t39,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39]);
    _ff_sum_59_mults_d19(t38,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38]);
    _ff_sum_57_mults_d19(t37,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37]);
    _ff_sum_56_mults_d18(t36,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_54_mults_d18(t35,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35]);
    _ff_sum_53_mults_d17(t34,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34]);
    _ff_sum_51_mults_d17(t33,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33]);
    _ff_sum_50_mults_d16(t32,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32]);
    _ff_sum_48_mults_d16(t31,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31]);
    _ff_sum_47_mults_d15(t30,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_45_mults_d15(t29,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29]);
    _ff_sum_44_mults_d14(t28,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_42_mults_d14(t27,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27]);
    _ff_sum_41_mults_d13(t26,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26]);
    _ff_sum_39_mults_d13(t25,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25]);
    _ff_sum_38_mults_d12(t24,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24]);
    _ff_sum_36_mults_d12(t23,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23]);
    _ff_sum_35_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_33_mults_d11(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_32_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_30_mults_d10(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_29_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_27_mults_d9(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_26_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_24_mults_d8(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_23_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_21_mults_d7(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_20_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_18_mults_d6(t11,f[0],f[1],f[2],f[3],f[4],f[5],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_17_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],f[5],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_15_mults_d5(t9,f[0],f[1],f[2],f[3],f[4],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_14_mults_d4(t8,f[0],f[1],f[2],f[3],f[4],t43,t44,t45,t46,t47,t48,t49,t50,t51,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_12_mults_d4(t7,f[0],f[1],f[2],f[3],t43,t44,t45,t46,t47,t48,t49,t50,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
    _ff_sum_11_mults_d3(t6,f[0],f[1],f[2],f[3],t43,t44,t45,t46,t47,t48,t49,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
    _ff_sum_9_mults_d3(t5,f[0],f[1],f[2],t43,t44,t45,t46,t47,t48,g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
    _ff_sum_8_mults_d2(t4,f[0],f[1],f[2],t43,t44,t45,t46,t47,g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
    _ff_sum_6_mults_d2(t3,f[0],f[1],t43,t44,t45,t46,g[0],g[1],g[2],g[3],f[2],f[3]);
    _ff_sum_5_mults_d1(t2,f[0],f[1],t43,t44,t45,g[0],g[1],g[2],f[1],f[2]);
    _ff_sum_3_mults_d1(t1,f[0],t43,t44,g[0],g[1],f[1]);
    _ff_sum_2_mults(o[0],f[0],t43,g[0],f[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22); _ff_set(o[23],t23); _ff_set(o[24],t24); _ff_set(o[25],t25); _ff_set(o[26],t26); _ff_set(o[27],t27); _ff_set(o[28],t28); _ff_set(o[29],t29); _ff_set(o[30],t30); _ff_set(o[31],t31); _ff_set(o[32],t32); _ff_set(o[33],t33); _ff_set(o[34],t34); _ff_set(o[35],t35); _ff_set(o[36],t36); _ff_set(o[37],t37); _ff_set(o[38],t38); _ff_set(o[39],t39); _ff_set(o[40],t40); _ff_set(o[41],t41); _ff_set(o[42],t42);
}

static inline void ff_poly_square_mult_x_mod_43 (ff_t o[43], ff_t f[43], ff_t g[42])
{
    register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85;

    _ff_square(t85,f[42]);
    _ff_dbl_mult(t84,f[41],f[42]);
    _ff_sum_3_mults_d1(t83,f[40],f[41],t85,g[41],f[41],f[42]);
    _ff_sum_4_mults_d2(t82,f[39],f[40],t84,t85,g[40],g[41],f[41],f[42]);
    _ff_sum_6_mults_d2(t81,f[38],f[39],f[40],t83,t84,t85,g[39],g[40],g[41],f[40],f[41],f[42]);
    _ff_sum_7_mults_d3(t80,f[37],f[38],f[39],t82,t83,t84,t85,g[38],g[39],g[40],g[41],f[40],f[41],f[42]);
    _ff_sum_9_mults_d3(t79,f[36],f[37],f[38],f[39],t81,t82,t83,t84,t85,g[37],g[38],g[39],g[40],g[41],f[39],f[40],f[41],f[42]);
    _ff_sum_10_mults_d4(t78,f[35],f[36],f[37],f[38],t80,t81,t82,t83,t84,t85,g[36],g[37],g[38],g[39],g[40],g[41],f[39],f[40],f[41],f[42]);
    _ff_sum_12_mults_d4(t77,f[34],f[35],f[36],f[37],f[38],t79,t80,t81,t82,t83,t84,t85,g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_13_mults_d5(t76,f[33],f[34],f[35],f[36],f[37],t78,t79,t80,t81,t82,t83,t84,t85,g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_15_mults_d5(t75,f[32],f[33],f[34],f[35],f[36],f[37],t77,t78,t79,t80,t81,t82,t83,t84,t85,g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_16_mults_d6(t74,f[31],f[32],f[33],f[34],f[35],f[36],t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_18_mults_d6(t73,f[30],f[31],f[32],f[33],f[34],f[35],f[36],t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_19_mults_d7(t72,f[29],f[30],f[31],f[32],f[33],f[34],f[35],t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_21_mults_d7(t71,f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_22_mults_d8(t70,f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_24_mults_d8(t69,f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_25_mults_d9(t68,f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_27_mults_d9(t67,f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_28_mults_d10(t66,f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_30_mults_d10(t65,f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_31_mults_d11(t64,f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_33_mults_d11(t63,f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_34_mults_d12(t62,f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_36_mults_d12(t61,f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_37_mults_d13(t60,f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_39_mults_d13(t59,f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_40_mults_d14(t58,f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_42_mults_d14(t57,f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_43_mults_d15(t56,f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_45_mults_d15(t55,f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_46_mults_d16(t54,f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_48_mults_d16(t53,f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_49_mults_d17(t52,f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_51_mults_d17(t51,f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_52_mults_d18(t50,f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_54_mults_d18(t49,f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_55_mults_d19(t48,f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_57_mults_d19(t47,f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_58_mults_d20(t46,f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_60_mults_d20(t45,f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_61_mults_d21(t44,f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);
    _ff_sum_63_mults_d21(t43,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42]);

    _ff_sum_63_mults_d21(t42,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41]);
    _ff_sum_63_mults_d20(t41,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],g[41],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40]);
    _ff_sum_61_mults_d20(t40,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],g[40],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39]);
    _ff_sum_60_mults_d19(t39,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],g[39],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38]);
    _ff_sum_58_mults_d19(t38,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],g[38],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37]);
    _ff_sum_57_mults_d18(t37,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],g[37],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36]);
    _ff_sum_55_mults_d18(t36,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],g[36],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35]);
    _ff_sum_54_mults_d17(t35,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],g[35],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34]);
    _ff_sum_52_mults_d17(t34,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],g[34],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33]);
    _ff_sum_51_mults_d16(t33,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],g[33],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32]);
    _ff_sum_49_mults_d16(t32,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],g[32],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31]);
    _ff_sum_48_mults_d15(t31,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],g[31],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30]);
    _ff_sum_46_mults_d15(t30,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],g[30],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29]);
    _ff_sum_45_mults_d14(t29,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],g[29],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28]);
    _ff_sum_43_mults_d14(t28,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],g[28],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27]);
    _ff_sum_42_mults_d13(t27,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],g[27],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26]);
    _ff_sum_40_mults_d13(t26,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],g[26],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25]);
    _ff_sum_39_mults_d12(t25,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],g[25],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24]);
    _ff_sum_37_mults_d12(t24,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],g[24],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23]);
    _ff_sum_36_mults_d11(t23,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],g[23],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22]);
    _ff_sum_34_mults_d11(t22,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],g[22],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21]);
    _ff_sum_33_mults_d10(t21,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],g[21],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20]);
    _ff_sum_31_mults_d10(t20,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],g[20],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
    _ff_sum_30_mults_d9(t19,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
    _ff_sum_28_mults_d9(t18,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
    _ff_sum_27_mults_d8(t17,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
    _ff_sum_25_mults_d8(t16,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
    _ff_sum_24_mults_d7(t15,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
    _ff_sum_22_mults_d7(t14,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
    _ff_sum_21_mults_d6(t13,f[0],f[1],f[2],f[3],f[4],f[5],f[6],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
    _ff_sum_19_mults_d6(t12,f[0],f[1],f[2],f[3],f[4],f[5],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11]);
    _ff_sum_18_mults_d5(t11,f[0],f[1],f[2],f[3],f[4],f[5],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[5],f[6],f[7],f[8],f[9],f[10]);
    _ff_sum_16_mults_d5(t10,f[0],f[1],f[2],f[3],f[4],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9]);
    _ff_sum_15_mults_d4(t9,f[0],f[1],f[2],f[3],f[4],t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[4],f[5],f[6],f[7],f[8]);
    _ff_sum_13_mults_d4(t8,f[0],f[1],f[2],f[3],t43,t44,t45,t46,t47,t48,t49,t50,t51,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7]);
    _ff_sum_12_mults_d3(t7,f[0],f[1],f[2],f[3],t43,t44,t45,t46,t47,t48,t49,t50,g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[3],f[4],f[5],f[6]);
    _ff_sum_10_mults_d3(t6,f[0],f[1],f[2],t43,t44,t45,t46,t47,t48,t49,g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5]);
    _ff_sum_9_mults_d2(t5,f[0],f[1],f[2],t43,t44,t45,t46,t47,t48,g[0],g[1],g[2],g[3],g[4],g[5],f[2],f[3],f[4]);
    _ff_sum_7_mults_d2(t4,f[0],f[1],t43,t44,t45,t46,t47,g[0],g[1],g[2],g[3],g[4],f[2],f[3]);
    _ff_sum_6_mults_d1(t3,f[0],f[1],t43,t44,t45,t46,g[0],g[1],g[2],g[3],f[1],f[2]);
    _ff_sum_4_mults_d1(t2,f[0],t43,t44,t45,g[0],g[1],g[2],f[1]);
    _ff_sum_3_mults(t1,f[0],t43,t44,g[0],g[1],f[0]);
    _ff_mult(o[0],t43,g[0]);
    _ff_set(o[1],t1); _ff_set(o[2],t2); _ff_set(o[3],t3); _ff_set(o[4],t4); _ff_set(o[5],t5); _ff_set(o[6],t6); _ff_set(o[7],t7); _ff_set(o[8],t8); _ff_set(o[9],t9); _ff_set(o[10],t10); _ff_set(o[11],t11); _ff_set(o[12],t12); _ff_set(o[13],t13); _ff_set(o[14],t14); _ff_set(o[15],t15); _ff_set(o[16],t16); _ff_set(o[17],t17); _ff_set(o[18],t18); _ff_set(o[19],t19); _ff_set(o[20],t20); _ff_set(o[21],t21); _ff_set(o[22],t22); _ff_set(o[23],t23); _ff_set(o[24],t24); _ff_set(o[25],t25); _ff_set(o[26],t26); _ff_set(o[27],t27); _ff_set(o[28],t28); _ff_set(o[29],t29); _ff_set(o[30],t30); _ff_set(o[31],t31); _ff_set(o[32],t32); _ff_set(o[33],t33); _ff_set(o[34],t34); _ff_set(o[35],t35); _ff_set(o[36],t36); _ff_set(o[37],t37); _ff_set(o[38],t38); _ff_set(o[39],t39); _ff_set(o[40],t40); _ff_set(o[41],t41); _ff_set(o[42],t42);
}

#ifdef __cplusplus
}
#endif

#endif
