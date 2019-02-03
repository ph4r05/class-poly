#include <assert.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include "gmp.h"
#include "ff.h"
#include "ffpoly.h"
#include "ffpolysmall.h"
#if FF_POLY_BIG
#include "ffpolybig.h"
#endif
#include "polyparse.h"
#include "cstd.h"

/*
    Copyright 2008-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// TODO: remove all dependencies on FF_POLY_MAX_DEGREE.  The code should handle arbitrary length polys in all cases, whether FF_POLY_BIG is defined or not.


// f and W can coincide (W = [a1,a2,a3,a4,a6])
void ff_poly_med_weierstrass (ff_t f[4], ff_t W[5])
{
	register ff_t t1, t2;

//printf ("[%ld,%ld,%ld,%ld,%ld]\n", _ff_get_ui(W[0]), _ff_get_ui(W[1]), _ff_get_ui(W[2]), _ff_get_ui(W[3]), _ff_get_ui(W[4]));	
	_ff_set_one (f[3]);
	_ff_mult (t1, W[0], _ff_half);					// t1 = a1/2
	_ff_square (t2, t1);							// c2 = a1^2/4
	_ff_addto (t2, W[1]);						// f2 = t2 = a2+a1^2/4
	_ff_mult(t1,t1,W[2]);						// c4 = a1a3/2
	_ff_add(f[1],t1,W[3]);						// f1 = a4 + a1a3/2
	_ff_mult(t1,W[2], _ff_half);					// t1 = a3/2
	_ff_square(t1,t1);							// c6 = a3^2/4
	_ff_add(f[0],t1,W[4]);						// f0 = a6 + a3^2/4
	_ff_set (f[2],t2);							// f2 = t2
//printf ("y^2 = x^3 + %ld*x^2 + %ld*x + %ld\n", _ff_get_ui(f[2]), _ff_get_ui(f[1]), _ff_get_ui(f[0]));
	// 3M+2S
}	

// f and W can coincide(W = [a1,a2,a3,a4,a6])
void ff_poly_short_weierstrass (ff_t f[4], ff_t W[5])
{
	register ff_t t1, t2, t3, c2, c4, c6;
	
	_ff_set_one (f[3]);
	_ff_set_zero (f[2]);
	
//printf ("[%ld,%ld,%ld,%ld,%ld]\n", _ff_get_ui(W[0]), _ff_get_ui(W[1]), _ff_get_ui(W[2]), _ff_get_ui(W[3]), _ff_get_ui(W[4]));
	_ff_mult (t1, W[0], _ff_half);					// t1 = a1/2
	_ff_square (c2, t1);							// c2 = a1^2/4
	_ff_addto (c2, W[1]);						// c2 = a2+a1^2/4 is the new a2
	_ff_mult(c4,t1,W[2]);						// c4 = a1a3/2
	_ff_addto(c4,W[3]);							// c4 = a4 + a1a3/2 is the new a4
	_ff_mult(t2,W[2], _ff_half);					// t2 = a3/2
	_ff_square(c6,t2);							// c6 = a3^2/4
	_ff_addto(c6,W[4]);							// c6 = a6 + a3^2/4 is the new a6 
//printf ("[%ld,%ld,%ld]\n", _ff_get_ui(c2), _ff_get_ui(c4), _ff_get_ui(c6));
	_ff_mult(t1,c2,_ff_third);						// t1 = c2/3
	_ff_mult(t2,t1,c2);							// t2 = c2^2/3
	_ff_sub (f[1], c4, t2);						// f[1] = c4 - c2^2/3
	_ff_mult(t2,t2,t1);							// t2 = c2^3/9
	_ff_mult(t3,t2,_ff_third);						// t3 = c2^3/27
	_ff_add(t2,t3,t3);							// t2 = 2/27*c2^3
	_ff_mult(t3,t1,c4);							// t3 = c2c4/3
	_ff_subfrom(t2,t3);							// t2 = 2/27*c2^3 - c2c4/3
	_ff_add (f[0], c6, t2);						// f[0] = c6 + 2/27*c2^3 - c2c4/3
//printf ("x^3+%ld*x+%ld\n", _ff_get_ui(f[1]), _ff_get_ui(f[0]));
	// 8M+2S
}	


int ff_poly_parse (ff_t f[], int maxd, char *expr)
{
	int i, d;
	
	d = ui_poly_parse_mod_p (f, maxd, expr, _ff_p);
	for ( i = 0 ; i <= d ; i++ ) _ff_set_ui (f[i], f[i]);
	for ( ; i <= maxd ; i++ ) _ff_set_zero (f[i]);		// pad out with zeroes 
	return d;
}


void ff_poly_print (ff_t f[], int d_f)
{
	char buf[65536];

	ff_poly_sprint (buf, f, d_f);
	out_printf ("%s\n", buf);
}

int ff_poly_sprint (char *s, ff_t f[], int d_f)
{
	char *t;
	int i;
	
	if ( d_f < 0 ) { strcpy (s, "[zero polynomial]");  return strlen(s); }
	t = s;
	if ( d_f >= 2 ) {
		t += gmp_sprintf (t, "[%lu*x^%d", _ff_get_ui(f[d_f]), d_f);
	} else if ( d_f == 1 ) {
		t += gmp_sprintf (t, "[%lu*x", _ff_get_ui(f[d_f]));
	} else {
		t += gmp_sprintf (t, "[%lu", _ff_get_ui(f[d_f]));
	}
	for ( i = d_f-1 ; i >= 0 ; i-- ) {
		if ( _ff_nonzero(f[i]) ) {
			if ( i >= 2 ) {
				t += gmp_sprintf (t, " + %lu*x^%d", _ff_get_ui(f[i]), i);
			} else if ( i == 1 ) {
				t += gmp_sprintf (t, " + %lu*x", _ff_get_ui(f[i]));
			} else {
				t += gmp_sprintf (t, " + %lu", _ff_get_ui(f[i]));
			}
		}
	}
	*t++ = ']';
	*t= '\0';
	return t-s;
}

/*
	Computes the quadratic twist of y^2=f(x) by a non-residue in F_p
	f and g may coincide
*/
void ff_poly_twist (ff_t g[], ff_t f[], int d)
{
	ff_t a, b;
	int i;

	ff_nonresidue(&b);
	if ( (d&0x1) ) {
		// in the odd degree case we can keep the leading coefficient - means that g is monic if f is.
		_ff_set (g[d], f[d]);					// leading coefficient unchanged - means twist is monic if f is
		_ff_set (a,b);
		for ( i = d-1 ; i >= 0 ; i-- ) {
			ff_mult (g[i], b, f[i]);				// g[i] = a^{d-i}f[i] for i = d-1 down to 0
			ff_mult (b, b, a);
		}
	} else {
		for ( i = 0 ; i <= d ; i++ ) _ff_mult(g[i],b,f[i]);
	}
	return;
}


// Naive Euclidean division of a by b, computes q and r such that a = qb+r with deg(r) < deg(b).  Overlap is ok.  Does not require any inversions if b is monic
// w is workspace for d_a+d_b+2 coeffiicents (enough space to store a and b)  If q and r do not overlap b, da+1 coeffs is enough
// q, pd_q, r, and pd_r are all optional
void _ff_poly_div (ff_t q[], int *pd_q, ff_t r[], int *pd_r, ff_t a[], int d_a, ff_t b[], int d_b, ff_t w[])
{
	ff_t *c;
	register ff_t s, t, x;
	register int i;
	int d_c, d_q, d_w, d_s;
	
	assert ( d_b >= 0 );
	if ( d_a < d_b ) { if ( pd_q ) *pd_q = -1; ff_poly_copy (r, pd_r, a, d_a); return; }	// bugfix to correctly handle null pd_q 08/06/13
	
	c = w + d_a+1;

	// if b overlaps q, make a copy of it
	if ( q == b ) { c = w+d_a+1;  ff_poly_copy (c, &d_c, b, d_b); } else { c = b; d_c = d_b; }
	
	ff_poly_copy (w, &d_w, a, d_a);
	d_q = -1;
	if ( ! _ff_one(c[d_c]) ) { _ff_invert (x, c[d_c]); } else { x = 0; }
	while ( d_w >= d_c ) {
		if ( x ) { _ff_mult (s, x, w[d_w]); } else { _ff_set (s, w[d_w]); }
		d_s = d_w - d_c;
		if ( q ) {
			if ( d_s > d_q ) {
				for ( i = d_q+1 ; i < d_s ; i++ ) _ff_set_zero (q[i]);
				d_q = d_s;
				_ff_set (q[d_q], s);
			} else {
				_ff_addto (q[d_s], s);
			}
		}
		for ( i = 0 ; i < d_c ; i++ ) { _ff_mult (t, s, c[i]); _ff_subfrom (w[i+d_s], t); }
		_ff_set_zero(w[d_w]);
		d_w = ff_poly_degree(w,d_w-1);
	}
	if ( pd_q ) *pd_q = d_q;
	if ( r ) ff_poly_copy (r, pd_r, w, d_w);
}


// computes h=af-bg of degree less than max(d_f,d_g) (f and g must both be non-zero).   note that h will not typically be monic, even if f and g are
// assuming d_f >= d_g,  h will be set to lc(g)*f-lc(f)*x^(d_f - d_g)*g. In this case gcd(f,h) = gcd(f,g), but gcd(g,h) = x^k*gcd(f,g) for some k <= d_f-d_g (typically k=0, but not always!)
static inline void ff_poly_reduce (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g)
{
	register ff_t t0,t1;
	register int i,j;
	
	assert ( d_f >= 0 && d_g >= 0 );
	if ( ! d_f || ! d_g ) { *pd_h = -1; return; }

	// avoid unnecessary copying, but allow overlap of h with f or g
	if ( d_f > d_g ) {
		_ff_neg(t0,f[d_f]);
		_ff_set(t1,g[d_g]);
		j = d_f-d_g;
		for ( i = d_f-1 ; i >= j ; i-- ) _ff_sum_2_mults(h[i],t0,t1,f[i],g[i-j]);
		for ( ; i >= 0 ; i-- ) ff_mult(h[i],f[i],t1);
		for ( i = d_f-1 ; i >= 0 ; i-- ) if ( ! _ff_zero(h[i]) ) break;
		*pd_h = i;
	} else {
		_ff_set(t0,f[d_f]);
		_ff_neg(t1,g[d_g]);		
		j = d_g-d_f;
		for ( i = d_g-1 ; i >= j ; i-- ) _ff_sum_2_mults(h[i],t0,t1,f[i-j],g[i]);
		for ( ; i >= 0 ; i-- ) ff_mult(h[i],g[i],t0);
		for ( i = d_g-1 ; i >= 0 ; i-- ) if ( ! _ff_zero(h[i]) ) break;
		*pd_h = i;
	}
}

// given nonzero polys s and t, reduces one modulo the other until one of them is zero
static inline void _ff_poly_gcd_reduce (ff_t s[], int *pd_s, ff_t t[], int *pd_t)
{
	while ( *pd_s >= 0 && *pd_t >= 0 ) {
		if ( *pd_s < *pd_t ) {
			ff_poly_reduce (t, pd_t, s, *pd_s, t, *pd_t);
		} else {
			ff_poly_reduce (s, pd_s, s, *pd_s, t, *pd_t);
		}
	}
}

// uses poly_reduce to get a multiple of the gcd (typically equal to the gcd), then uses ff_poly_gcd to finish the job if necessary
// overlap is ok
int ff_poly_gcd_reduce (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b)
{
	ff_t *s, *t, *x;
	int d_s, d_t, d_x, k;

	if ( d_b > d_a ) { _swap (a, b, s);  _swap (d_a, d_b, d_s); }
	if ( d_b < 0 ) { assert (d_a >= 0 ); ff_poly_copy(g,pd_g, a,d_a); return d_a; }
	
	// remove any factors of x from the gcd
	for ( k = 0 ; _ff_zero(a[k]) && _ff_zero(b[k]) ; k++ );
	a += k;  d_a -= k;  b += k;  d_b -= k;
	if ( ! d_a ) {_ff_set_one(g[k]);  *pd_g = k;  for ( k-- ; k >= 0 ; k-- ) _ff_set_zero(g[k]);  return *pd_g; }
	
	// reduce into s and t rather than copying
	s = ff_poly_alloc (d_a);  t = ff_poly_alloc (d_a);
	ff_poly_reduce (s, &d_s, a, d_a, b, d_b);
	if ( d_s < 0 ) { x = b; d_x = d_b; goto done; }
	ff_poly_reduce (t, &d_t, b, d_b, s, d_s);
	_ff_poly_gcd_reduce (s, &d_s, t, &d_t);
	if ( d_s < 0 ) { x = t; d_x=d_t; } else { x = s; d_x = d_s; }
done:
	// remove any factors of x that may have been introduced in the process of reduction, we know they are not in the gcd, since we removed all factors of x above
	while ( _ff_zero(x[0]) ) { x++; d_x--; };
	ff_poly_copy (g+k,0,x,d_x);  *pd_g = d_x + k;
	for ( k-- ; k >= 0 ; k-- ) _ff_set_zero (g[k]);
	ff_poly_free (s, d_a);  ff_poly_free (t, d_a);
	if ( ! *pd_g ) _ff_set_one(g[0]);		// make g monic in the trivial case
	return *pd_g;
}

// note that g is not made monic unless it is trivial (note that even if a and b are monic, g might not be)
int ff_poly_gcd (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b)
{
	int d_s, d_t;
	ff_t *s, *t, *w;

	if ( d_b > d_a ) { _swap (a, b, t);  _swap (d_a, d_b, d_t); }
	if ( d_b < 0 ) { assert ( d_a >= 0 );  ff_poly_copy (g, pd_g, a, d_a); return d_a; }
	
	w = mem_alloc ((d_a+3*d_b+1)*sizeof(ff_t));  s = w + d_a+d_b+2;  t = s +d_b;
	_ff_poly_div (0, 0, s, &d_s, a, d_a, b, d_b, w);
	if ( d_s < 0 ) { ff_poly_copy (g, pd_g, b, d_b); goto done; }
	_ff_poly_div (0, 0, t, &d_t, b, d_b, s, d_s, w);
	while ( d_s >= 0 && d_t >= 0 ) {
		if ( d_s > d_t ) {
			_ff_poly_div (0, 0, s, &d_s, s, d_s, t, d_t, w);
		} else {
			_ff_poly_div (0, 0, t, &d_t, t, d_t, s, d_s, w);
		}
	}
	if ( d_s < 0 ) ff_poly_copy (g, pd_g, t, d_t); else ff_poly_copy (g, pd_g, s, d_s);
done:
	mem_free (w);
	if ( ! *pd_g ) _ff_set_one(g[0]);
	return *pd_g;
}


// [CCANT] algorithm 3.2.2, makes output monic.  None of the polynomials can overlap.  g may be null.  degree of gcd is returned in any case
// computes u, v, and g such that g=gcd(a,b) = u*a+v*b
int ff_poly_gcdext (ff_t g[], int *pd_g, ff_t u[], int *pd_u, ff_t v[], int *pd_v, ff_t a[], int d_a, ff_t b[], int d_b)
{
	ff_t *f, *q, *r, *t, *t2, *v1, *v3, *w;
	int d_f, d_q, d_r, d_t, d_t2, d_v1, d_v3;
	int n;

	// GCd(0,b) = b = 0*a+1*b
	if ( d_a < 0 ) { if ( g ) ff_poly_copy (g, pd_g, b, d_b); *pd_u = -1;  _ff_set_one(v[0]); *pd_v = 0; return d_b; }
		
	// GCD(a,0) = a = 1*a+0*b
	if ( d_b < 0 ) { if ( g ) ff_poly_copy (g, pd_g, a, d_a); *pd_u= 0;  _ff_set_one(u[0]); *pd_v = -1; return d_a; }
		
	// GCD(a,c) = 1 = 0*a+1/c*c for c a nonzero constant
	if ( ! d_b ) { if ( g ) { *pd_g = 0;  _ff_set_one(g[0]); } *pd_u = -1; _ff_invert(v[0], b[0]); *pd_v = 0; return 0; }
	
	// GCD(c,b) = 1 = 1/c*c+0*b for c a nonzero constant
	if ( ! d_a ) { if ( g ) { *pd_g = 0;  _ff_set_one(g[0]); } _ff_invert (u[0], a[0]); *pd_u = 0; *pd_v = -1; return 0; }
	
	// for simplicity just set n to the max coeffs needed by any of the polys we are using
	n = ( d_a > d_b ? d_a : d_b );
	w = mem_alloc (13*n*sizeof(ff_t));
	f = w+2*n;  q = f+n;  r = q+n;  t = r+n;  t2 = t+2*n;  v1 = t2+2*n;  v3 = v1+2*n;	// t, t2, v1 and v3 may need space for degreee 2n polys

	*pd_u = 0;  	_ff_set_one(u[0]);				// u = 1
	ff_poly_copy (f, &d_f, a, d_a);					// f = a
	d_v1 = -1;								// v1 = 0
	ff_poly_copy (v3, &d_v3, b, d_b);				// v3 = b	

	while ( d_v3 >= 0 ) {
		_ff_poly_div (q, &d_q, r, &d_r, f, d_f, v3, d_v3, w);
		ff_poly_mult (t2, &d_t2, v1, d_v1, q, d_q);
		ff_poly_sub (t, &d_t, u, *pd_u, t2, d_t2);		// t = u-v1q
		// this is a lot of copying, but it means we only use one multiplication
		ff_poly_copy (u, pd_u,  v1, d_v1);			// u = v1
		ff_poly_copy (v1, &d_v1, t, d_t);			// v1 = t
		ff_poly_copy (f, &d_f, v3, d_v3);			// f = v3
		ff_poly_copy (v3, &d_v3, r, d_r);			// v3 = r
	}
	ff_poly_mult (v3, &d_v3, a, d_a, u, *pd_u);
	ff_poly_sub (v1, &d_v1, f, d_f, v3, d_v3);
	_ff_poly_div (v, pd_v, 0, 0, v1, d_v1, b, d_b,w);	// v = (f-au)/b 		(exact division)
	// make GCD monic
	if ( ! _ff_one (f[d_f]) ) {
		_ff_set (w[1], f[d_f]);
		ff_poly_monic (f, &d_f, f, d_f);
		_ff_invert (w[0], w[1]);
		ff_poly_scalar_mult (u, pd_u, w[0], u, *pd_u);
		ff_poly_scalar_mult (v, pd_v, w[0], v, *pd_v);
	}
	if ( g ) ff_poly_copy (g, pd_g, f, d_f);
	mem_free (w);
	return d_f;
}

int ff_poly_inv_mod (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g)
{
	ff_t *v;
	int d_v;
	int n;
	
	v = ff_poly_alloc (d_f);
	n = ff_poly_gcdext (0, 0, h,  pd_h, v, &d_v, f, d_f, g, d_g);
	ff_poly_free (v, d_f);
	return ( n ? 0 : 1);
}


// g = a^e mod f, requires deg(a) < deg(f)
void ff_poly_pow_mod (ff_t g[], int *pd_g, ff_t a[], int d_a, unsigned long e, ff_t f[], int d_f)
{
	ff_poly_modulus_t mod;
	
	assert (d_a < d_f);
	ff_poly_mod_setup (mod, f, d_f);
	ff_poly_pow_modulus (g, pd_g, a, d_a, e, mod);
	ff_poly_mod_clear(mod);
}

// g = a^e mod f, requires deg(a) < deg(f)
void ff_poly_pow_modulus (ff_t g[], int *pd_g, ff_t a[], int d_a, unsigned long e,  ff_poly_modulus_t mod)
{
	ff_t *h;
	int d_h;
	int t;

	assert (d_a < mod->d);
	if ( d_a < 0 || !e ) { _ff_set_one(g[0]);  *pd_g = 0;  return; }

	h = ff_poly_alloc (2*mod->d-2);
	ff_poly_copy (h, &d_h, a, d_a);
	t = ui_len (e) - 2;
	while ( t >= 0 ) {				// TODO: use a faster exponentiation algorithm here
		ff_poly_square (h, h, d_h);
		ff_poly_mod_reduce (h, &d_h, h, 2*d_h, mod);
		if ( e&(1UL<<t) ) {
			ff_poly_mult (h, &d_h, h, d_h, a, d_a);
			ff_poly_mod_reduce (h, &d_h, h, d_h, mod);
		}
		t--;
	}
	ff_poly_copy (g, pd_g, h, d_h);
	ff_poly_free (h, 2*mod->d-2);
}

// g = a^e mod f, requires deg(a) < deg(f)
void ff_poly_pow_mod_mpz (ff_t g[], int *pd_g, ff_t a[], int d_a, mpz_t e, ff_t f[], int d_f)
{
	ff_poly_modulus_t mod;
	
	assert (d_a < d_f);
	ff_poly_mod_setup (mod, f, d_f);
	ff_poly_pow_modulus_mpz (g, pd_g, a, d_a, e, mod);
	ff_poly_mod_clear(mod);
}
	

// g = a^e mod f, requires deg(a) < deg(f)
void ff_poly_pow_modulus_mpz (ff_t g[], int *pd_g, ff_t a[], int d_a, mpz_t e,  ff_poly_modulus_t mod)
{
	ff_t *h;
	int d_h;
	int t;

	assert (d_a < mod->d);
	if ( d_a < 0 || ! mpz_sgn (e) ) { _ff_set_one(g[0]);  *pd_g = 0;  return; }
	
	h = ff_poly_alloc (2*mod->d-2);
	ff_poly_copy (h, &d_h, a, d_a);
	t = mpz_sizeinbase (e, 2) - 2;
	while ( t >= 0 ) {					// TODO: use a faster exponentiation algorithm here
		ff_poly_square (h, h, d_h);
		ff_poly_mod_reduce (h, &d_h, h, 2*d_h, mod);
		if ( mpz_tstbit (e, t) ) {
			ff_poly_mult (h, &d_h, h, d_h, a, d_a);
			ff_poly_mod_reduce (h, &d_h, h, d_h, mod);
		}
		t--;
	}
	ff_poly_copy (g, pd_g, h, d_h);
	ff_poly_free (h, 2*mod->d-2);
}


// attempts to compute h = (f) mod g using slow Tonelli-Shanks, where g is assumed to be irreducible
// h and f may overlap
int ff_poly_sqrt_modulus (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_poly_modulus_t mod)
{
	ff_t *a, *ai, *g;
	int d_g, d_a, d_ai;
	mpz_t m, n, e, w;
	int i, s;

//printf ("computing sqrt of "); ff_poly_print (f, d_f); printf (" mod "); ff_poly_print (mod->g, mod->d);
	if ( d_f < 0 ) { *pd_h = -1; return 1; }
	g = ff_poly_alloc (2*mod->d-2);
	mpz_init (m);  mpz_init (n);
	mpz_ui_pow_ui (m, _ff_p, mod->d);  mpz_sub_ui (m, m, 1);		// m = p^d-1 is the order of the multiplicative group of Fp[x]/(mod)
	mpz_div_2exp (n, m, 1);									// n = (p^d-1)/2
	
	// Determine whether f has a square root or not by computing f^((p^d-1)/2)
	ff_poly_pow_modulus_mpz (g, &d_g, f, d_f, n, mod);
	if ( ! ff_poly_is_one (g, d_g) ) { if ( h == f ) ff_poly_free (g, mod->d-1); mpz_clear(m);  mpz_clear(n); return 0; }
	
	// Get a random non-residue a
	a = ff_poly_alloc (mod->d-1);  ai = ff_poly_alloc (mod->d-1);
	d_a = mod->d-1;
	do {
		ff_poly_randomize (a, d_a);
		ff_poly_pow_modulus_mpz (g, &d_g, a, d_a, n, mod);
	} while ( ff_poly_is_one (g, d_g) );
	if ( ! ff_poly_inv_modulus (ai, &d_ai, a, d_a, mod) ) { printf ("Unable to compute inverse in ff_poly_mod_sqrt, modulus poly not irreducible!\n");  ff_poly_print (mod->g, mod->d); abort(); }
	
	for ( s = 1 ; !mpz_tstbit (n,0) ; s++ ) mpz_div_2exp (n,n,1);		// m = n*2^s with n odd

	mpz_init (e);  mpz_init (w);
	for ( i = 2 ; i <= s ; i++ ) {								// this loop computes e such that f=a^e mod the 2-Sylow
		ff_poly_pow_modulus_mpz (g, &d_g, ai, d_ai, e, mod);
		ff_poly_mult (g, &d_g, g, d_g, f, d_f);
		ff_poly_mod_reduce (g, &d_g, g, d_g, mod);
		mpz_div_2exp (w, m, i);								// w = (p^d-1)/2^i
		ff_poly_pow_modulus_mpz (g, &d_g, g, d_g, w, mod);
		if ( ! ff_poly_is_one (g, d_g) ) mpz_setbit (e, i-1);
	}
	ff_poly_pow_modulus_mpz (g, &d_g, ai, d_ai, e, mod);
	ff_poly_mult (g, &d_g, g, d_g, f, d_f);						// g = f*a^-e
	ff_poly_mod_reduce (g, &d_g, g, d_g, mod);
	mpz_add_ui (n, n, 1);  mpz_div_2exp (n, n, 1);
	ff_poly_pow_modulus_mpz (g, &d_g, g, d_g, n, mod);			// g = f*a^(-e(n+1)/2) is the square-root of the odd part of f
	mpz_div_2exp (e, e, 1);
	ff_poly_pow_modulus_mpz (ai, &d_ai, a, d_a, e, mod);			// ai = a^(e/2) is the square root of the even part of f
	ff_poly_mult (g, &d_g, g, d_g, ai, d_ai);						// g = f*a^(-e(n+1)/2)*a^(e/2) = f*a^(-en) is the square root of f
	ff_poly_mod_reduce (h, pd_h, g, d_g, mod);
	ff_poly_mult (g, &d_g, h, *pd_h, h, *pd_h);					// verify that h^2 = f
	ff_poly_mod_reduce (g, &d_g, g, d_g, mod);
	assert ( ff_poly_equal (g, d_g, f, d_f) );
	ff_poly_free (g, 2*mod->d-2);  ff_poly_free (a, mod->d-1);  ff_poly_free (ai, mod->d-1);
	mpz_clear (m); mpz_clear (n); mpz_clear (e); mpz_clear (w);
	return 1;
}

int ff_poly_sqrt_mod (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g)
{
	ff_poly_modulus_t mod;
	int ret;
	
	ff_poly_mod_setup (mod, g, d_g);
	ret = ff_poly_sqrt_modulus (h, pd_h, f, d_f, mod);
	ff_poly_mod_clear (mod);
	return ret;
}


// computes x^n mod f for depressed monic f 
void ff_poly_xn_mod (ff_t g[], int *pd_g, unsigned long e, ff_t f[], int d_f)
{
	register int i;
	
	assert ( d_f > 1 && _ff_one(f[d_f]) );
	if ( ! e ) { _ff_set_one(g[0]);  *pd_g = 0;  return; }
	
	if ( e < d_f ) { *pd_g = e;  _ff_set_one(g[e]);  for ( i = e-1 ; i >= 0 ; i-- ) _ff_set_zero(g[i]);  return; }

	if ( d_f <= FF_POLY_SMALL_DEGREE && _ff_zero(f[d_f-1]) ) {
		ff_t nf[FF_POLY_SMALL_DEGREE+1];

		assert (_ff_zero(f[d_f-1]));	// for small degrees we require the poly to be depressed 
		// for small mod poly operations all assume the modulus is of the form x^n  - c_{n-2}x^{n-2} - ... - c_1*x - c_0, so negate the non-leading coefficients (leading one is implicit)
		ff_poly_neg(nf,0,f,d_f-2);
		switch (d_f) {
		case 2: ff_exp_ui (g,nf,e>>1); if ( e&1 ) { _ff_set(g[1],g[0]); _ff_set_zero(g[0]); } else { _ff_set_zero(g[1]); } break;
		case 3: ff_poly_xn_mod_d3 (g,e,nf); break;
		case 4: ff_poly_xn_mod_d4 (g,e,nf); break;
		case 5: ff_poly_xn_mod_d5 (g,e,nf); break;
		case 6: ff_poly_xn_mod_d6 (g,e,nf); break;
		case 7: ff_poly_xn_mod_d7 (g,e,nf); break;
		case 8: ff_poly_xn_mod_d8 (g,e,nf); break;
		case 10: ff_poly_xn_mod_d10 (g,e,nf); break;
		case 11: ff_poly_xn_mod_d11 (g,e,nf); break;
		case 12: ff_poly_xn_mod_d12 (g,e,nf); break;
		case 13: ff_poly_xn_mod_d13 (g,e,nf); break;
		case 15: ff_poly_xn_mod_d15 (g,e,nf); break;
		case 17: ff_poly_xn_mod_d17 (g,e,nf); break;
		case 19: ff_poly_xn_mod_d19 (g,e,nf); break;
		case 23: ff_poly_xn_mod_d23 (g,e,nf); break;
		case 29: ff_poly_xn_mod_d29 (g,e,nf); break;
		case 31: ff_poly_xn_mod_d31 (g,e,nf); break;
#if FF_SUPERFAST
		// the compiler takes forever to deal with these inline, so we conditionally compile them in when asked
		case 37: ff_poly_xn_mod_d37 (g,e,nf); break;
		case 41: ff_poly_xn_mod_d41 (g,e,nf); break;
		case 43: ff_poly_xn_mod_d43 (g,e,nf); break;
#endif
		default: ff_poly_xn_mod_small (g,e,nf,d_f);
		}
		*pd_g = ff_poly_degree(g,d_f-1);
	} else {
		ff_poly_modulus_t mod;
		ff_poly_mod_setup (mod, f, d_f);
		ff_poly_xn_modulus (g, pd_g, e,mod);
		ff_poly_mod_clear(mod);
	}
}

// note, we don't need FF_POLY_BIG defined to use this, ff_poly_square is just slower when it can't call ff_poly_square_big
void ff_poly_xn_modulus  (ff_t g[], int *pd_g, unsigned long e, ff_poly_modulus_t mod)
{
	register ff_t *w;
	register int t, d_w;
	int d_g;
	
	if ( mod->d <= FF_POLY_SMALL_DEGREE && _ff_zero(mod->g[mod->d-1]) ) { ff_poly_xn_mod (g, pd_g, e, mod->g, mod->d); return; }
	if ( e < mod->d ) { *pd_g = e;  _ff_set_one(g[e]);  for ( d_g=e-1 ; d_g >= 0 ; d_g-- ) _ff_set_zero(g[d_g]);  return; }

	_ff_set_one(g[1]);  _ff_set_zero(g[0]);  d_g = 1;  w = mod->w1;
	t = ui_len(e)-2;
	while ( t >= 0 ) {
		if ( e&(1UL<<t) ) {
			ff_poly_square (w+1,g,d_g);   _ff_set_zero(w[0]);
			d_w = ff_poly_degree(w,2*d_g+1);			
		} else {
			ff_poly_square (w,g,d_g);   d_w = ff_poly_degree(w,2*d_g);
		}
		ff_poly_mod_reduce (g,&d_g,w,d_w,mod);
		t--;
	}
	*pd_g = d_g;
}

// note, we don't need FF_POLY_BIG defined to use this, ff_poly_square is just slower when it can't call ff_poly_square_big
void ff_poly_xn_modulus_mpz (ff_t g[], int *pd_g, mpz_t e, ff_poly_modulus_t mod)
{
	register ff_t *w;
	register int t, d_w;
	int d_g;
	
	t = mpz_sizeinbase(e,2);
	if ( t <= 64 && mod->d <= FF_POLY_SMALL_DEGREE ) { ff_poly_xn_mod (g, pd_g, mpz_get_ui(e), mod->g, mod->d); return; }
	_ff_set_one(g[1]);  _ff_set_zero(g[0]);  d_g = 1;  w = mod->w1;
	t -= 2;
	while ( t >= 0 ) {
		if ( mpz_tstbit(e,t) ) {
			ff_poly_square (w+1,g,d_g);   _ff_set_zero(w[0]);
			d_w = ff_poly_degree(w,2*d_g+1);			
		} else {
			ff_poly_square (w,g,d_g);   d_w = ff_poly_degree(w,2*d_g);
		}
		ff_poly_mod_reduce (g,&d_g,w,d_w,mod);
		t--;
	}
	*pd_g = d_g;
}

void ff_poly_xn_mod_mpz (ff_t g[], int *pd_g, mpz_t e, ff_t f[], int d_f)
{
	ff_poly_modulus_t mod;
	register int i;
	
	if ( mpz_cmp_ui (e, d_f) < 0 )  { *pd_g = (int)mpz_get_ui(e);  _ff_set_one(g[*pd_g]);  for ( i = *pd_g-1 ; i >= 0 ; i-- ) _ff_set_zero(g[i]);  return; }
	ff_poly_mod_setup (mod,f, d_f);
	ff_poly_xn_modulus_mpz (g, pd_g, e,mod);
	ff_poly_mod_clear(mod);
}


// computes (x+a)^n mod f for depressed monic f of deg > 1
void ff_poly_xpan_mod (ff_t g[], int *pd_g, ff_t a, unsigned long e, ff_t f[], int d_f)
{
	register int i;

	assert ( d_f > 1 && _ff_one(f[d_f]) && _ff_zero(f[d_f-1]) );
	if ( ! e ) { _ff_set_one(g[0]);  *pd_g = 0;  return; }

	if ( d_f <= FF_POLY_SMALL_DEGREE ) {
		ff_t nf[FF_POLY_SMALL_DEGREE+1];

		// for small mod poly operations all assume the modulus is of the form x^n  - c_{n-2}x^{n-2} - ... - c_1*x - c_0, so negate the non-leading coefficients (leading one is implicit)
		for ( i = 0 ; i < d_f ; i++ ) _ff_neg(nf[i],f[i]);
		switch (d_f) {
		case 2: err_printf ("Degree 2 not supported in ff_poly_xpan_mod_ui, use ff_exp_ui\n"); abort();
		case 3: ff_poly_xpan_mod_d3 (g,a,e,nf); break;
		case 4: ff_poly_xpan_mod_d4 (g,a,e,nf); break;
		case 5: ff_poly_xpan_mod_d5 (g,a,e,nf); break;
		case 6: ff_poly_xpan_mod_d6 (g,a,e,nf); break;
		case 7: ff_poly_xpan_mod_d7 (g,a,e,nf); break;
		case 8: ff_poly_xpan_mod_d8 (g,a,e,nf); break;
		case 10: ff_poly_xpan_mod_d10 (g,a,e,nf); break;
		case 11: ff_poly_xpan_mod_d11 (g,a,e,nf); break;
		case 12: ff_poly_xpan_mod_d12 (g,a,e,nf); break;
		case 13: ff_poly_xpan_mod_d13 (g,a,e,nf); break;
		case 15: ff_poly_xpan_mod_d15 (g,a,e,nf); break;
		case 17: ff_poly_xpan_mod_d17 (g,a,e,nf); break;
		case 19: ff_poly_xpan_mod_d19 (g,a,e,nf); break;
		case 23: ff_poly_xpan_mod_d23 (g,a,e,nf); break;
		case 29: ff_poly_xpan_mod_d29 (g,a,e,nf); break;
		case 31: ff_poly_xpan_mod_d31 (g,a,e,nf); break;
#if FF_SUPERFAST
		case 37: ff_poly_xpan_mod_d37 (g,a,e,nf); break;
		case 41: ff_poly_xpan_mod_d41 (g,a,e,nf); break;
		case 43: ff_poly_xpan_mod_d43 (g,a,e,nf); break;
#endif
		default: ff_poly_xpan_mod_small (g,a,e,nf,d_f); break;
		}
		*pd_g = ff_poly_degree(g,d_f-1);
		return;
	} else {
		ff_poly_modulus_t mod;

		ff_poly_mod_setup (mod, f, d_f);
		ff_poly_xpan_modulus (g, pd_g, a,e,mod);
		ff_poly_mod_clear(mod);
	}
}

void ff_poly_xpan_modulus  (ff_t g[], int *pd_g, ff_t a, unsigned long e, ff_poly_modulus_t mod)
{
	register ff_t *w,t0;
	register int i, t, d_w;
	int d_g;

//printf ("Computing (x+%ld)^%d mod f (p=%ld) ", _ff_get_ui(a), e, _ff_p);  ff_poly_print(f,d_f);
	_ff_set_one(g[1]);  _ff_set(g[0],a); d_g = 1;
	w = mod->w1;
	t = ui_len(e)-2;
	while ( t >= 0 ) {
		if ( e&(1UL<<t) ) {
			ff_poly_square (w+1,g,d_g);   
			_ff_mult(w[0],a,w[1]);
			for ( i = 1 ; i < 2*d_g+1 ; i++ ) { _ff_mult(t0,a,w[i+1]); _ff_addto(w[i],t0); }
			d_w = ff_poly_degree(w,2*d_g+1);
		} else {
			ff_poly_square (w,g,d_g);   d_w = ff_poly_degree(w,2*d_g);
		}
		ff_poly_mod_reduce (g,&d_g,w,d_w,mod);
		t--;
	}
	*pd_g = d_g;
}

// Algorithm IPT from Shoup "A Computational Introduction to Number Theory and Algebra" p. 463
// returns true if f is irreducible, false otherwise.  If pnroots is non-null then *pnroots is set to the number of Fp-roots (whether f is irreducible or not).
int ff_poly_irreducible (ff_t f[], int d_f, int *pnroots)
{
	ff_poly_modulus_t mod;
	ff_t X[2], *h, *w;
	int d_h, d_w;
	int i;

	assert ( d_f > 0 && _ff_one(f[d_f]) );
	if ( d_f == 1 ) {  if ( pnroots ) *pnroots = 1;  return 1;  }
	if ( ! _ff_one(f[d_f]) ) { w = ff_poly_alloc (d_f); ff_poly_monic (w, 0, f, d_f); f = w; } else { w = 0; }	
	
	_ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);		// create the poly -x
	h = ff_poly_alloc (d_f-1);  ff_poly_mod_setup (mod, f, d_f);
	ff_poly_xn_modulus (h, &d_h, _ff_p, mod);  i = 1;			// h = x^p mod f
	for ( ;; ) {
		ff_poly_add (mod->w, &d_w, h, d_h, X, 1);			// w = x^(p^i)-x
		ff_poly_gcd (mod->w, &d_w, mod->w, d_w, f, d_f);		// w = gcd(t,f)
		if ( d_w || ++i > d_f/2 ) break;
		ff_poly_pow_modulus (h, &d_h, h, d_h, _ff_p, mod);		// h = x^(p^i) mod f
	}
	ff_poly_free (h, d_f-1);  ff_poly_mod_clear (mod);
	if ( w ) ff_poly_free (w, d_f);
	if ( i > d_f/2 ) return 1;
	if ( pnroots ) { *pnroots = ( i==1 ? d_w : 0 ); }
	return 0;
}

// computes the factorization pattern of  a (not necessarily monic) polynomial f, sets counts[d] to the number of irreducible factors of degree d, returns total count
// if root is non-null and poly has at least one linear factor, r will be set to a root
int ff_poly_factorization_pattern_and_root (int counts[], ff_t f[], int d_f, ff_t *root)
{
	ff_poly_modulus_t mod;
	ff_t *g, *h, X[2], *w;
	int d_g, d_h, d_w;
	int i, k, n;

	assert ( d_f > 0 );
	if ( d_f == 1 ) { if ( counts ) counts[1] = 1; if ( root ) ff_poly_find_root (root, f, d_f); return 1; }
	if ( ! _ff_one(f[d_f]) ) { w = ff_poly_alloc (d_f); ff_poly_monic (w, 0, f, d_f); f = w; } else { w = 0; }	
	
	if ( counts ) for ( n =0, i = 0 ; i <= d_f ; i++ ) counts[i] = 0;
	_ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);				// create the poly -x
	g = ff_poly_alloc (d_f);  h = ff_poly_alloc (d_f-1);
	ff_poly_monic (g, &d_g, f, d_f);								// g is a working copy of f made monic
	ff_poly_mod_setup (mod, g, d_g);
	i = 1;
	ff_poly_xn_modulus (h, &d_h, _ff_p, mod);						// h = x^p mod f
	for (;;) {
		ff_poly_add (mod->w, &d_w, h, d_h, X, 1);					// w = x^{p^i}-x mod g
		for(;;) {
			ff_poly_gcd_reduce (mod->w, &d_w, mod->w, d_w, g, d_g);
			if ( ! d_w ) break;
			k = d_w/i;  assert (k*i==d_w);
			ff_poly_monic(mod->w,0,mod->w,d_w);
			if ( i == 1 && d_w && root && ! ff_poly_find_root (root, mod->w, d_w) ) { err_printf ("Error, couldn't find root of a totaly split poly in ff_poly_factorization_pattern_with_root\n"); abort(); }
			ff_poly_div (g, &d_g, 0, 0, g, d_g, mod->w, d_w);			// TODO: optimize this exact monic division
			if ( counts ) counts[i] += k;
			n += k;
		}
		if ( ++i > d_g/2 ) break;
		ff_poly_pow_modulus (h, &d_h, h, d_h, _ff_p, mod);				// h = x^(p^i) mod f
	}
	ff_poly_mod_clear (mod);  ff_poly_free (g, d_f);  ff_poly_free (h, d_f-1);
	if ( w ) ff_poly_free (w, d_f);
	if ( d_g ) { if ( counts ) counts[d_g]++; n++; }
	return n;
}


// completely factors a square-free monic polynomial f that is the product of d/k irreducible polynomials of degree k, using Cantor-Zassenhaus (Algs 14.8 and 14.10 in GG)
// r must have space for d=(d/k)*k coefficients:  a concatenated list of d/k implicitly monic polys of degree k will be stored in w, each using just d/k coeffs (monic leading coeff is omitted).
// r and f cannot overlap
void ff_poly_equal_degree_factorization (ff_t r[], ff_t f[], int d, int k)
{
	ff_poly_modulus_t mod;
	ff_t *g;
	mpz_t n;
	int d_g;

	assert ( k > 0 && d > 0 && !(d%k) && _ff_one(f[d]) );
	if ( d == k ) { ff_poly_copy (r, 0, f, d-1); return; }
//printf ("EFD: "); ff_poly_print (f,d);
	g = ff_poly_alloc (d);	
	ff_poly_mod_setup (mod, f, d);
	mpz_init (n); mpz_ui_pow_ui (n, _ff_p, k); mpz_sub_ui(n, n, 1); mpz_div_2exp (n, n, 1);				// n = (p^k-1)/2
	do {		
		ff_poly_randomize (g, d-1);  d_g = ff_poly_degree (g, d-1);								// pick a random poly mod f
		if ( d_g <= 0 ) continue;
		ff_poly_pow_modulus_mpz (g, &d_g, g, d_g, n, mod);									// fg= g^n mod
		if ( d_g < 1 ) continue;
		_ff_dec (g[0]);																	// fg= g^n-1 mod 
		ff_poly_gcd_reduce (g, &d_g, g, d_g, f, d);
	} while ( d_g <= 0 );
	mpz_clear (n);
	ff_poly_mod_clear (mod);
	
	assert ( !(d_g%k) );
	
	ff_poly_monic (g, 0, g, d_g);
	ff_poly_equal_degree_factorization (r, g, d_g, k);
	ff_poly_div (g, &d_g, 0, 0, f, d, g, d_g);
	ff_poly_equal_degree_factorization (r+d-d_g, g, d_g, k);
	ff_poly_free (g, d);
}


// completely factors a monic polynomial f, using Cantor-Zassenhaus (Algs 14.13 in GG)
// r must have space for d_f coefficients:  a concatenated list of implicitly monic polys f_1,...,f_t will be written to r, with n[i] set to deg f_i.
// r and f cannot overlap.  the number of factors t is returned
int ff_poly_factors (ff_t r[], int n[], ff_t f[], int d)
{
	ff_poly_modulus_t mod;
	ff_t *g, *h, *s, *w, X[2];
	int d_g, d_h, d_s;
	int i, j, k, o, t;

	assert ( d > 0 && _ff_one(f[d]) && r != f );
//printf ("Factoring monic poly "); ff_poly_print (f, d);

	// first remove any factors of x
	for ( t = o = 0; _ff_zero(f[0]) ; o++ ) { f++;  d--;  _ff_set_zero(r[o]); n[t++] = 1; }
	if ( !d ) return t;
	if ( d == 1 ) { _ff_set(r[o],f[0]); n[t++] = 1; return t; }

	_ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);				// create the poly -x
	w = mem_alloc (10*(d+1)*sizeof(*w));
	g = w+2*d+2;  h = g+d+1;  s = h+d+1;
	ff_poly_copy (g, &d_g, f, d);									// g is a working copy of f that will be whittled down
	ff_poly_mod_setup (mod, g, d_g);
	ff_poly_xn_modulus (h, &d_h, _ff_p, mod);						// h = x^p mod f
	i = 1;
	for (;;) {
		ff_poly_add (s, &d_s, h, d_h, X, 1);							// w = x^{p^k}-x mod g
		for(;;) {
			ff_poly_gcd_reduce (s, &d_s, s, d_s, g, d_g);
			if ( ! d_s ) break;
			k = d_s/i;  assert (k*i == d_s);
			ff_poly_monic(s,0,s,d_s);
			_ff_poly_div (g, &d_g, 0, 0, g, d_g, s, d_s, w);
//printf ("Remaining = "); ff_poly_print (g, d_g);
			ff_poly_equal_degree_factorization (r+o, s, d_s, i);
			for ( j = 0 ; j < k ; j++ ) n[t++] = i;
			o += d_s;
		}
		if ( ++i > d_g/2 ) break;
		ff_poly_pow_modulus (h, &d_h, h, d_h, _ff_p, mod);				// h = x^(p^k) mod f
	}
	if ( d_g ) { ff_poly_copy (r+o, 0, g, d_g-1);  n[t++] = d_g; }
	ff_poly_mod_clear (mod);  mem_free (w);
//printf ("Found %d factors of degrees: ", t);  for ( i = 0 ; i < t ; i++ ) printf ("%d ", n[i]); puts ("");
	return t;
}


int ff_poly_distinct_roots (ff_t r[], ff_t f[], int d_f)						//  only returns distinct roots, f need not be monic
{
	ff_t *w;
	int k;

	assert (d_f >= 0);
	if ( d_f > 2 && ! _ff_one(f[d_f]) ) { w = ff_poly_alloc (d_f); ff_poly_monic (w, 0, f, d_f); f = w; } else { w = 0; }	

	switch (d_f) {
	case 0: return 0;
	case 1: return ff_poly_roots_d1(r,f);
	case 2: k=ff_poly_roots_d2(r,f,d_f);  if ( k==2 && _ff_equal(r[0],r[1]) ) k=1;  return k;
	case 3:
		k=ff_poly_roots_d3(r,f);
		if ( k==3 ) {
			if ( _ff_equal(r[0],r[1]) ) { if ( _ff_equal(r[1],r[2]) ) { k=1; } else { _ff_set (r[1],r[2]); k =2; } }
			else if ( _ff_equal(r[0],r[2]) ) k = 2;
			else if ( _ff_equal(r[1],r[2]) ) k = 2;
		}
		if ( w ) ff_poly_free (w,d_f);
		return k;
	case 4: k=ff_poly_roots_d4(r,f);
		if ( k==2 && _ff_equal(r[0],r[1]) ) k=1;
		if ( k==4 ) {
			if ( _ff_equal(r[0],r[1]) ) {
				if ( _ff_equal(r[1],r[2]) ) {
					if ( _ff_equal(r[2],r[3]) ) {
						k= 1;
					} else {
						_ff_set(r[1],r[3]);
						k = 2;
					}
				} else {
					if ( _ff_equal(r[1],r[3]) ) {
						_ff_set(r[1],r[2]);  k = 2;
					} else {
						_ff_set(r[1],r[3]); k = (_ff_equal(r[1],r[2])?2:3);
					}
				}
			} else if ( _ff_equal(r[0],r[2]) ) {
				if ( _ff_equal(r[0],r[3]) || _ff_equal(r[1],r[3]) ) {
					k = 2;
				} else {
					_ff_set(r[2],r[3]); k= 3;
				}
			} else if ( _ff_equal(r[0],r[3]) ) {
				k = ( _ff_equal(r[1],r[2]) ? 2 : 3 );
			} else if ( _ff_equal(r[1],r[2]) ) {
				_ff_set(r[2],r[3]);
				k = ( _ff_equal(r[1],r[2]) ? 2 : 3 );
			} else if ( _ff_equal(r[1],r[3]) ) {
				k = 3;
			} else {
				k = ( _ff_equal(r[2],r[3]) ? 3 : 4 );
			}
		}
		if ( w ) ff_poly_free (w,d_f);
		return k;
	}
	k = _ff_poly_distinct_roots(r,f,d_f,0);
	if ( w ) ff_poly_free (w,d_f);
	return k;
}

// this function is not particularly optimized - an area for future work.  requires f monic
// split flag indicates that f is a product of linear polynomials (after the initial call, all recursive calls will have this flag set)
// FF_POLY_ONE_ROOT is used to indicate that one *randomly chosen* root is desired -- in this case r will only have space for one entry
// returns the number of distinct roots (even when FF_POLY_ONE_ROOT is set).
int _ff_poly_distinct_roots  (ff_t *r, ff_t f[], int d_f, int flags)
{
	ff_t *g, *h, *w, X[2], c, c2;
	int d_g, d_h;
	int i, n;

	assert ( d_f >= 0 );
	if ( ! d_f ) return 0;

	if ( d_f < 5 ) {
		if ( (flags&FF_POLY_ONE_ROOT) ) {													// AVS bug fix 03/06/2013, don't go past first entry of r when FF_POLY_ONE_ROOT is set
			ff_t rr[4];
			n = ff_poly_distinct_roots(rr,f,d_f);
			if ( n ) r[0] = rr[ff_randomm_ui(n)];
			return n;
		}
		return ff_poly_distinct_roots(r,f,d_f);
	}
	if ( ! _ff_one(f[d_f]) ) { w = ff_poly_alloc (d_f); ff_poly_monic (w, 0, f, d_f); f = w; } else { w = 0; }
	g = ff_poly_alloc (d_f);  h = ff_poly_alloc (d_f);
	if ( (d_f%_ff_p) ) {
		ff_poly_depress_monic(&c,g,f,d_f);	d_g=d_f;										// depress f  to g(x) = f(x-c) to speed up mod f operations, if possible
	} else {
		ff_poly_copy (g, &d_g, f, d_f);  _ff_set_zero (c);
	}
	if ( ! (flags&FF_POLY_SPLIT) ) {
		ff_poly_xn_mod (h, &d_h, _ff_p, g, d_g);											// h = x^p mod f(x-c)
		_ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);
		ff_poly_addto (h, &d_h, X, 1);													// h = x^p-x mod f(x-c)  (let addto handle the addition of -x to get the degree right)
		ff_poly_gcd_reduce (g, &d_g, h, d_h, g, d_g);										// g=gcd(g,h)
		ff_poly_monic (g, 0, g, d_g);
		if ( d_g < 5 ) { n = _ff_poly_distinct_roots(r,g,d_g,flags|FF_POLY_SPLIT); goto translate_roots; }
		if ( (d_g%_ff_p) ) {
			ff_poly_depress_monic_inplace (&c2,g,d_g);
			_ff_addto(c,c2);
		}
	}
	n = d_g;
	do {
		ff_random(X);																// pick a random a
		ff_poly_xpan_mod (h, &d_h, X[0], _ff_p>>1, g, d_g);								// h = (x+a)^{(p-1)/2} mod g
		if ( d_h < 1 ) continue;
		_ff_dec(h[0]);																// h = (x+a)^{(p-1)/2}-1 mod g
		ff_poly_gcd_reduce (h, &d_h, h, d_h, g, d_g);										// h=gcd(g,h)
	} while ( d_h <= 0 );
	ff_poly_monic(h,0,h,d_h);
	if ( (flags&FF_POLY_ONE_ROOT) ) {
		if ( d_h > d_g-d_h ) ff_poly_div (h, &d_h, 0, 0, g, d_g, h, d_h);
		_ff_poly_distinct_roots(r, h, d_h, flags|FF_POLY_SPLIT);
	} else {
		_ff_poly_distinct_roots(r, h, d_h, flags|FF_POLY_SPLIT);
		ff_poly_div (h, &d_h, 0, 0, g, d_g, h, d_h);
		_ff_poly_distinct_roots(r+d_g-d_h, h, d_h, flags|FF_POLY_SPLIT);
	}
translate_roots:	// if we get here, n is not zero (in fact n is at least 5)
	ff_poly_free (g, d_f);  ff_poly_free (h, d_f);
	if ( w ) ff_poly_free (w, d_f);
	if ( ! _ff_zero(c) ) {
		_ff_subfrom(r[0],c);
		if ( ! (flags&FF_POLY_ONE_ROOT) ) for ( i = 1 ; i < n ; i++ ) _ff_subfrom(r[i],c);
	}
	return n;
}

// returns all Fp-roots of f, with multiplicity.  f need not be monic.  Repeated roots will be adjacent to eachother in the returned list.
int ff_poly_roots (ff_t r[], ff_t f[], int d_f)
{
	ff_t *g, *w, t0;
	int i, k, n, d_g;

	assert ( d_f >= 0 );
	if ( d_f > 2 && ! _ff_one(f[d_f]) ) { w = ff_poly_alloc (d_f); ff_poly_monic (w, 0, f, d_f); f = w; } else { w = 0; }	
	switch (d_f) {
	case 0: return 0;
	case 1: return ff_poly_roots_d1(r,f);
	case 2: return ff_poly_roots_d2(r,f,d_f);
	case 3: return ff_poly_roots_d3(r,f);
	case 4: return ff_poly_roots_d4(r,f);
	}
	k = _ff_poly_distinct_roots(r,f,d_f,0);
	if ( !k || k == d_f ) { if ( w ) ff_poly_free (w,d_f); return k; }
	g = ff_poly_alloc (d_f);  ff_poly_copy (g, &d_g, f, d_f);
	for ( i = 0, n = k ; i < k ; i++ ) {
		ff_poly_remove_root (g, g, d_g, r+i);  d_g--;
		for (;;) {
			ff_poly_eval (&t0, g, d_g, r+i);
			if ( ! _ff_zero(t0) ) break;
			ff_poly_remove_root (g, g, d_g, r+i);  d_g--;
			_ff_set (r[n],r[i]); n++;
		}
	}
	ff_organize (r, n);
	ff_poly_free (g, d_f);
	if ( w ) ff_poly_free (w, d_f);
	return n;
}

// returns the degree of gcd(f,x^p-x) which will be the number of distinct roots.  does not require f to be monic
int ff_poly_count_distinct_roots (ff_t f[], int d_f)
{
	ff_t *g, *h, *w, X[4];
	int d_g, d_h;

	assert ( d_f >= 0 );
	if ( ! _ff_one(f[d_f]) ) { w = ff_poly_alloc (d_f); ff_poly_monic (w, 0, f, d_f); f = w; } else { w = 0; }
	
	if ( d_f < 5 ) { d_g = ff_poly_distinct_roots (X, f, d_f);  if ( w ) ff_poly_free (w,d_f); return d_g; }

	g = ff_poly_alloc (d_f);  h = ff_poly_alloc (d_f-1);
	if ( (d_f%_ff_p) ) { ff_poly_depress_monic (X, g, f, d_f);  d_g = d_f; } else { ff_poly_copy (g, &d_g, f, d_f); }
	ff_poly_xn_mod (h, &d_h, _ff_p, g, d_g);					// h = x^p mod g
	_ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);
	ff_poly_addto (h, &d_h, X, 1);							// h = x^p-x mod g
	ff_poly_gcd_reduce (g, &d_g, g, d_g, h, d_h);				// reduces g and h mod eachother until one is zero and the other is the gcd
	ff_poly_free (g, d_f);  ff_poly_free (h, d_f-1);
	if ( w ) ff_poly_free (w, d_f);
	return d_g;
}

// returns the total number of Fp-roots of f, with multiplicity.
int ff_poly_count_roots (ff_t f[], int d_f)
{
	ff_t *g, *h, X[4], *w;
	int d_g, d_h, n;

	assert ( d_f >= 0 );
	if ( ! _ff_one(f[d_f]) ) { w = ff_poly_alloc (d_f); ff_poly_monic (w, 0, f, d_f); f = w; } else { w = 0; }

	if ( d_f < 5 ) { n = ff_poly_roots (X, f, d_f);  if ( w ) ff_poly_free (w,d_f); return n; }

	g = ff_poly_alloc (d_f);  h = ff_poly_alloc (d_f);
	if ( (d_f%_ff_p) ) { ff_poly_depress_monic (X, g, f, d_f);  d_g = d_f; } else { ff_poly_copy (g, &d_g, f, d_f); }
	ff_poly_xn_mod (h, &d_h, _ff_p, g, d_g);					// h = x^p mod g
	_ff_set_zero (X[0]);  _ff_set_one (X[1]);  ff_negate (X[1]);
	ff_poly_addto (h, &d_h, X, 1);							// h = x^p-x mod g
	ff_poly_gcd_reduce (h, &d_h, g, d_g, h, d_h);				// h = gcd(g,h)
	for ( n = d_h ; d_h ; n += d_h ) {
		ff_poly_div (g, &d_g, 0, 0, g, d_g, h, d_h);			// TODO: optimize this exact monic division
		ff_poly_gcd_reduce (h, &d_h, g, d_g, h, d_h);
	}
	ff_poly_free (g, d_f);  ff_poly_free (h, d_f);
	if ( w ) ff_poly_free (w, d_f);
	return n;
}


// Requires f monic, handles degrees 2, 3, and 4, and also x^d+ax+b cases for d=5 or 7
int _ff_poly_quick_discriminant (ff_t disc[1], ff_t f[], int degree)
{
	ff_t A, B, temp1, temp2, temp3;
	
	if ( ! _ff_one (f[degree]) ) return 0;
	switch (degree) {
	case 1: _ff_set_one (disc[0]);  return 1;
	case 2: _ff_square(temp1,f[1]);  _ff_add(temp2,f[0],f[0]);  _ff_x2(temp2); _ff_sub(disc[0],temp1,temp2);  return 1;
	case 3:
		_ff_square(temp1,f[1]);  ff_mult(temp1,temp1,f[1]);  _ff_x2(temp1);  _ff_x2(temp1);								// temp1 = 4f1^3
		_ff_square(temp2,f[0]); _ff_set_i(temp3,27); ff_mult(temp2,temp2,temp3);										// temp2 = 27f2^2
		_ff_add(disc[0],temp1,temp2);										
		ff_negate(disc[0]);																					// disc = -4f1^3-27f0^2
		if ( ! _ff_zero(f[2]) ) {
			_ff_square(temp2,f[2]);  _ff_mult(temp3,f[2],f[0]);  ff_mult(temp1,temp2,temp3); _ff_x2(temp1);  _ff_x2(temp1);		// temp1 = 4f2^3f0
			_ff_square(temp3,f[1]);  _ff_mult(temp3,temp2,temp3);  _ff_subfrom(temp3,temp1);							// temp3 = f2^2f1^2 - 4f2^3f0
			_ff_set_ui(temp2,18);  _ff_mult (temp1,f[2],temp2);  _ff_mult(temp2,f[1],f[0]); ff_mult(temp1,temp1,temp2);		// temp1 = 18f0f1f2
			_ff_addto(disc[0],temp1);
			_ff_addto(disc[0],temp3);																			// disc = -4f2^3f0+f2^2f1^2+18f2f1f0-4f1^3-27f0^2
		}
		return 1;
	case 4:
		if ( ! _ff_zero(f[3]) ) {
			_ff_square(A,f[3]);  ff_mult(A,A,f[0]);  _ff_square(temp1,f[1]); _ff_addto(A,temp1);								// A = f3^2f0 + f1^2
			_ff_mult(temp1,f[2],f[0]);  _ff_x2(temp1);  _ff_x2(temp1); _ff_subfrom(A,temp1);								// A =  f3^2f0 + f1^2 - 4f2f0
			_ff_mult(B,f[3],f[1]);  _ff_add(temp1,f[0],f[0]);  _ff_x2(temp1);  _ff_subfrom(B,temp1);							// B = f3f1 - 4f0
		} else {
			 _ff_square(A,f[1]);  _ff_mult(temp1,f[2],f[0]);  _ff_x2(temp1);  _ff_x2(temp1); _ff_subfrom(A,temp1);				// A =  f1^2 - 4f2f0
			 _ff_add(B,f[0],f[0]);  _ff_x2(B);  ff_negate(B);															// B = -4f0
		}
		_ff_square(temp3,f[2]);  _ff_mult(temp1,temp3,f[2]);  ff_mult(temp1,temp1,A); _ff_x2(temp1); _ff_x2(temp1);			// temp1 = 4f2^3A
		_ff_square(temp2,B);  ff_mult(disc[0],temp2,temp3); _ff_subfrom (disc[0],temp1);									// disc = -4f2^3A + f2^2B^2
		ff_mult(temp2,temp2,B);  _ff_x2(temp2);  _ff_x2(temp2); _ff_subfrom(disc[0],temp2);								// disc = -4f2^3A + f2^2B^2 - 4*B^3
		_ff_set_ui(temp1,18); _ff_mult(temp3,temp1,f[0]); _ff_mult(temp1,A,B); ff_mult(temp1,temp1,temp3);					// temp1 = 18f0BA
		_ff_set_ui(temp2,27);  _ff_square(temp3,A);  ff_mult(temp2,temp2,temp3); _ff_subfrom(temp1,temp2);					// temp1 = 18f0BA - 27A^2
		_ff_addto(disc[0],temp1);																				 // disc = -4f2^3A + f2^2B^2 + 18f0BA - B^3 - 27A^2
		return 1;
	case 5:
		if ( ! _ff_zero(f[4]) || ! _ff_zero(f[3]) || ! _ff_zero(f[2]) ) return 0;
		_ff_square(temp1,f[1]);
		ff_square(temp1,temp1);
		ff_mult(temp1,temp1,f[1]);
		_ff_set_i(temp2,-256);
		_ff_mult(disc[0],temp1,temp2);
		_ff_square(temp1,f[0]);
		ff_square(temp1,temp1);
		_ff_set_i(temp2,-3125);
		ff_mult(temp1,temp1,temp2);
		_ff_addto(disc[0],temp1);
		return 1;
	case 7:
		if ( ! _ff_zero(f[6]) || ! _ff_zero(f[5]) || ! _ff_zero(f[4]) || ! _ff_zero(f[3]) || ! _ff_zero(f[2]) ) return 0;
		_ff_square(temp1,f[1]);
		ff_square(temp2,temp1);
		ff_mult(temp1,temp1,temp2);
		ff_mult(temp1,temp1,f[1]);
		_ff_set_i(temp2,-46656);
		_ff_mult(disc[0],temp1,temp2);
		_ff_square(temp1,f[0]);
		ff_square(temp2,temp1);
		ff_mult(temp1,temp1,temp2);
		_ff_set_i(temp2,-823543);
		ff_mult(temp1,temp1,temp2);
		_ff_addto(disc[0],temp1);
		return 1;
	}
	return 0;
}

// very fast code to test whether gcd(f,g)=1 or not
int ff_poly_trivial_gcd (ff_t f[], int d_f, ff_t g[], int d_g)
{
	ff_t s[FF_POLY_MAX_DEGREE+1], t[FF_POLY_MAX_DEGREE+1];
	register ff_t a, b;
	register ff_t *sp, *se, *tp, *te, *xp;
	register int i;
	
	se = s+d_f;
	for ( sp = se, tp=f+d_f ; sp >= s ; ) *sp-- = *tp--;
	te = t+d_g;
	for ( tp = te, sp=g+d_g ; tp >= t ; ) *tp-- = *sp--;

	while ( se > s && te > t ) {
		i = (se-s) - (te-t);
		if ( i >= 0 ) {
			a = *te;  b = *se;  xp = s+i;
			ff_negate(b);
			for ( sp = s ; sp < xp ; ) { _ff_mult (*sp,*sp,a); sp++; }
			for ( tp = t ; tp <= te ; tp++ ) { _ff_sum_2_mults (*sp, a, b, *tp, *sp); sp++; }
			for ( se-- ; se >= s && !*se ; se-- );
		} else {
			a = *se;  b = *te;  xp = t-i;
			ff_negate(b);
			for ( tp = t ; tp < xp ; ) { _ff_mult (*tp,*tp,a); tp++; }
			for ( sp = s ; sp <= se ; sp++ ) { _ff_sum_2_mults (*tp, a, b, *sp, *tp); tp++;  }
			for ( te-- ; te >= t && !*te ; te-- );
		}
	}
	if ( se < s || te < t ) return 0;
	return 1;
}

// tests whether discriminat is zero by computing gcd(f,f')
int ff_poly_discriminant_nonzero (ff_t f[], int degree)
{
	ff_t df[FF_POLY_MAX_DEGREE];
	int d_df;

	// check for easy special cases first
	if ( _ff_one(f[degree]) && _ff_poly_quick_discriminant (df, f, degree) ) return ( _ff_zero(df[0]) ? 0 : 1 );
	
	ff_poly_derivative (df, &d_df, f, degree);
	return ff_poly_trivial_gcd (f, degree, df, d_df);
}

/*
	Given monic f(x) of degree d, computes linear translation g(x)=(x-f[d-1]/d) which will be a monic poly with d-1 coefficient zero.
	The parameter s is optional, if specified it will be set to f[d-1]/d so that g(x)=f(x-s).  Replaces f by g.
*/
void ff_poly_depress_monic_inplace (ff_t s[1], ff_t f[], int d_f)
{
	ff_t z;
	register ff_t t,c;
	register int i, j;

	if ( _ff_zero(f[d_f-1]) ) { if ( s ) _ff_set_zero(*s); return; }
	ff_invert_small_int(&z,d_f);
	_ff_mult(c,f[d_f-1],z);
	if ( s ) _ff_set(*s,c);								// g(x) = f(x-s) so that g(x+s)=f(x)	consistent with ff_poly_depress_cubic
	ff_negate(c);									// c = -f[d-1]/d, we will replace f by f(x+c)
	for ( i = d_f ; i ; i-- ) {
		for ( j = i ; j < d_f ; j++ ) { _ff_mult(t,c,f[j]); _ff_addto(f[j-1],t); }
		_ff_addto(f[j-1],c);
	}
}


/*
	For monic f, substitute x <- (x - f_{d-1}/[d*f_d]) to make f_{d-1} term zero if necessary
	Returns 1 if poly modified, 0 if not
*/
int mpz_poly_depress_monic_inplace (mpz_t f[], int d)
{
	mpz_t c, s, t, w;
	int i, j;
	
	if ( d < 1 ) return 0;
	if ( ! mpz_sgn(f[d-1]) ) return 0;
	if ( mpz_cmp_ui (f[d],1) != 0 ) return 0;

	// dynamically allocate mpz's here, we don't expect to be called more than once per curve
	mpz_init (c); mpz_init (s); mpz_init (t); mpz_init (w);
	
	// set s to -d*f[d-1]
	mpz_mul_ui (s, f[d-1], d);  mpz_neg (s, s);

	// compute d^(2d) * f ( (x+sd) / d^2 ) = sum_{i=0}^d sum_{k=i}^d binom(k,i) * d^(2*d-2*j) * s^(j-i) * f[j] x^i
	for ( i = 0 ; i < d-1 ; i++ ) {
		mpz_set_ui (w, 1);
		mpz_set_ui (c, 0);
		for ( j = i ; j <= d ; j++ ) {
			mpz_ui_pow_ui (t, d, 2*(d-j));
			mpz_mul (t, t, f[j]);
			mpz_mul (t, t, w);
			mpz_mul_ui (t, t, ui_binomial (j, i));
			mpz_add (c, c, t);
			mpz_mul (w, w, s);
		}
		mpz_set (f[i], c);
	}
	mpz_set_ui (f[d-1], 0);
	
	mpz_clear (c); mpz_clear (s); mpz_clear (t); mpz_clear (w);

	return 1;
}

/*
	Use hardwired fast schoolbook for multiplying small polys.
	Call zn_poly to handle big polys via ff_poly_mult_big, if FF_POLY_BIG is defined
	We prefer d_f <= d_g but will swap if needed
	Overlap is OK
*/
void ff_poly_mult (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g)
{
	ff_t t[FF_POLY_MAX_DEGREE+1];		// only used if f is small, and if FF_POLY_BIG is defined, only entries up to FF_POLY_MULT_SMALL_DEGREE will be used
	ff_t u[2*FF_POLY_MAX_DEGREE+1];
	ff_t *x, *y;
	register int i, d;
	
	if ( f == g && d_f == d_g ) { ff_poly_square (h, f, d_f); *pd_h = ff_poly_degree(h,2*d_f); return; }
	if ( d_f < 0 || d_g < 0 ) { *pd_h = -1; return; }
	if ( d_f > d_g ) { _swap (d_g, d_f, d); _swap (f, g, x); }
	if ( d_g < 10 ) {
		if ( d_f+1 < d_g ) {
			for ( i = 0 ; i <= d_f ; i++ ) _ff_set(t[i],f[i]);
			while ( i < d_g ) { _ff_set_zero(t[i]); i++; }
			x = t;  y = u;
		} else {
			x = f;  y = h;
		}
		if ( d_f < d_g ) {
			switch ( d_g ) {
			case 1: ff_poly_mult_1_2 (y, x, g); break;
			case 2: ff_poly_mult_2_3 (y, x, g); break;
			case 3: ff_poly_mult_3_4 (y, x, g); break;
			case 4: ff_poly_mult_4_5 (y, x, g); break;
			case 5: ff_poly_mult_5_6 (y, x, g); break;
			case 6: ff_poly_mult_6_7 (y, x, g); break;
			case 7: ff_poly_mult_7_8 (y, x, g); break;
			case 8: ff_poly_mult_8_9 (y, x, g); break;
			case 9: ff_poly_mult_9_10 (y, x, g); break;
			}
		} else {
			switch ( d_g ) {
			case 0: _ff_mult(h[0],f[0],g[0]); break;
			case 1: ff_poly_mult_2_2 (y, x, g); break;
			case 2: ff_poly_mult_3_3 (y, x, g); break;
			case 3: ff_poly_mult_4_4 (y, x, g); break;
			case 4: ff_poly_mult_5_5 (y, x, g); break;
			case 5: ff_poly_mult_6_6 (y, x, g); break;
			case 6: ff_poly_mult_7_7 (y, x, g); break;
			case 7: ff_poly_mult_8_8 (y, x, g); break;
			case 8: ff_poly_mult_9_9 (y, x, g); break;
			case 9: ff_poly_mult_10_10 (y, x, g); break;
			}
		}
#if FF_POLY_BIG
	} else if ( d_g <= FF_POLY_MULT_SMALL_DEGREE ) {
#else
	} else {
		if ( d_g > FF_POLY_MAX_DEGREE ) { err_printf ("Exceeded FF_POLY_MAX_DEGREE = %d in ff_poly_mult\n", FF_POLY_MAX_DEGREE); abort(); }
#endif
		// pad f out to the same size as g if needed
		if ( d_f < d_g ) {
			for ( i = 0 ; i <= d_f ; i++ ) _ff_set(t[i],f[i]);
			while ( i <= d_g ) { _ff_set_zero(t[i]); i++; }
			x = t;  y = u;
		} else {
			x = f;  y = h;
		}
		ff_poly_mult_small(y,x,g,d_g);
#if FF_POLY_BIG
	} else {
		y = h;
		ff_poly_mult_big(h,f,d_f,g,d_g);
#endif
	}
	*pd_h = ff_poly_degree (y,d_f+d_g);
	if ( y != h ) ff_poly_copy (h, 0, y, *pd_h);
}

// we don't inline these because we want ffpoly.c to be the only module that depends on FF_POLY_BIG
void ff_poly_square (ff_t g[], ff_t f[], int d_f)
{
#if FF_POLY_BIG
	if ( d_f > FF_POLY_SQUARE_SMALL_DEGREE ) { ff_poly_square_big(g,f,d_f); return; }
#endif
	ff_poly_square_small(g,f,d_f);
}

void ff_poly_mod (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t f[], int d_f)
{
	if ( d_f > FF_POLY_MOD_SMALL_DEGREE )  { ff_poly_mod_big (g,pd_g,a,d_a,f,d_f); return; }
	// note that we don't use ff_poly_mod_small here because it requires the modulus to be monic and depressed, and negates the signs of the coefficients.
	ff_poly_div (0,0,g,pd_g,a,d_a,f,d_f);
}


// Too avoid alloc/free for each call, use ff_poly_from_roots_big directly
void ff_poly_from_roots (ff_t f[], ff_t r[], int d)
{
	if ( d <= 64 ) { ff_poly_from_roots_small (f, r, d); return; }	
#if FF_POLY_BIG
	{
		ff_t *t, *g;
		t = ff_poly_from_roots_big_workspace (d);
		if ( f == r ) g = mem_alloc ((d+1)*sizeof(*g)); else g = f;
		ff_poly_from_roots_big (g, r, d, t);
		if ( f == r ) { memcpy (f, g, (d+1)*sizeof(*g)); mem_free (g); }
		mem_free (t);
	}
#else
	ff_poly_from_roots_naive (f, r, d);		// this is painfully slow, use FF_POLY_BIG if d gets at all big
#endif
}


/*
	Compute the inverse of f mod x^n using Algorithm 9.3 in von zur Gathen & Gerhard.
	Incorporates optimization to efficiently treat n not a power of 2.
	g and f cannot overlap, work must hold 2*n-1 entries
*/
void  ff_poly_inverse_mod_xn (ff_t *g, int *d_g, ff_t *f, int d_f, int n, ff_t *work)
{
	ff_t *h;
	register ff_t t0;
	register int i, j, k, d, e, r;
	int d_h;
	
	if ( ! n ) { *d_g = -1; return; }
	if ( _ff_zero(f[0]) ) { err_printf ("Uninvertible polynomial in ff_poly_inverse_mod_xn!\n"); ff_poly_print(f,d_f); abort(); }
	if ( _ff_one(f[0]) ) _ff_set_one(g[0]); else _ff_invert(g[0],f[0]);
	if ( n==1 ) { *d_g = 0; return; }
	h = work;
	d = 1; r = ui_len(n);  if ( (1<<(r-1)) == n ) r--;
	for (r--;;r--) {
		e = ui_ceil_ratio (n, (1<<r));
		ff_poly_square(h,g,d-1);
		j = ( 2*d-2 < e ? 2*d-2 : e-1 );
		k = ( d_f < e ? d_f : e-1 );
		ff_poly_mult(h,&d_h,h,j,f,k);
		for ( i = 0 ; i < d ; i++ ) {
			_ff_add(t0,g[i],g[i]);
			_ff_sub(g[i],t0,h[i]);
		}
		if ( e > n ) e = n;
		for ( ; i < e ; i++ ) _ff_neg(g[i],h[i]);
		if ( e == n ) break;
		d = e;
	}
	*d_g = ff_poly_degree(g,n-1);
}

/*
	Compute log(f) = - \sum_{i=1}^{n-1}(1/i)(1-f)^i mod x^n.
	f and g may overlap
*/
void ff_poly_log_mod_xn (ff_t *g, int *pd_g, ff_t *f, int d_f, int n)
{
	ff_t t[FF_POLY_MAX_DEGREE+1];
	ff_t u[FF_POLY_MAX_DEGREE+1];
	ff_t h[2*FF_POLY_MAX_DEGREE+1];
	int d_t, d_u, d_h;
	
	if ( d_f > FF_POLY_MAX_DEGREE || n > FF_POLY_MAX_DEGREE ) { err_printf ("Exceeded FF_POLY_MAX_DEGREE = %d in ff_poly_log_mod_xn\n", FF_POLY_MAX_DEGREE); abort(); }
	if ( d_f < 0 || ! _ff_one(f[0]) ) { err_printf ("Constant term must be one in ff_poly_log_mod_xn\n"); ff_poly_print(f,d_f); abort(); }
	if ( n <= 1 ) { *pd_g = -1; return; }
	ff_poly_mod_xn(f,&d_f,n+1);
	ff_poly_derivative (t, &d_t, f, d_f);
	ff_poly_inverse_mod_xn(u, &d_u, f, d_f, n, h);
	ff_poly_mult(h,&d_h,t,d_t,u,d_u);
	ff_poly_mod_xn(h,&d_h,n-1);
	ff_poly_antiderivative(g,pd_g,h,d_h);
	return;
}


/*
	Compute exp(f) = 1+f+1/2*f^2+1/6*f^3+... mod x^n using Brent's algorithm (see p. 4 of BSMS 2008)
	Could be optimized to better handle n not a power of 2.
	f and g may overlap
*/
void  ff_poly_exp_mod_xn (ff_t *g, int *pd_g, ff_t *f, int d_f, int n)
{
	ff_t t[FF_POLY_MAX_DEGREE+1];
	ff_t v[FF_POLY_MAX_DEGREE+1];
	ff_t h[2*FF_POLY_MAX_DEGREE+1];
	register int i, m;
	int d_h, d_t;
	
	if ( d_f > FF_POLY_MAX_DEGREE || n > FF_POLY_MAX_DEGREE ) { err_printf ("Exceeded FF_POLY_MAX_DEGREE = %d in ff_poly_exp_mod_xn\n", FF_POLY_MAX_DEGREE); abort(); }
	if ( ! n ) { *pd_g = -1; return; }
	if ( d_f < 0 || n == 1 ) { *pd_g = 0; _ff_set_one(g[0]); return; }
	if ( ! _ff_zero(f[0]) ) { err_printf ("Error, polynomial has nonzero constant term in ff_poly_exp_mod_xn\n"); ff_poly_print(f,d_f); abort(); }
//printf("exping "); ff_poly_print(f,d_f);
	_ff_set_one(t[0]); _ff_set(t[1],f[1]);  d_t = m = 2;
	while ( m < n ) {
//printf("exp(m=%d,n=%d) t:",m,n); ff_poly_print(t,d_t);
		m += m;
		ff_poly_log_mod_xn (h,&d_h,t,d_t,m);
//printf("log(t):"); ff_poly_print(h,d_h);
		_ff_set_one(v[0]);  if ( d_h >= 0 ) _ff_subfrom(v[0],h[0]);
		for ( i = 1 ; i <= d_h ; i++ ) _ff_sub(v[i],f[i],h[i]);
		for ( ; i <= d_f && i < m ; i++ ) _ff_set(v[i],f[i]);
		ff_poly_mult(h,&d_h,t,d_t,v,i);
		for ( i = 0 ; i < m ; i++ ) _ff_set(t[i],h[i]);
		d_t = ff_poly_degree(t,i-1);
	}
	ff_poly_mod_xn(t,&d_t,n);
//printf("exp final "); ff_poly_print(t,d_t);
	ff_poly_copy(g,pd_g,t,d_t);
}

// setup up data structure to efficiently reduce polys of degree up to max modulo g, see GG Alg 9.5.  max defaults to 2*d_g-1
void ff_poly_mod_setup_max (ff_poly_modulus_t mod, ff_t *g, int d_g, int max)
{
	register int i;
	int d, dw;
	
	if ( max < 2*d_g-1 ) max = 2*d_g-1;
	mod->d = d_g;  mod->e = max-d_g;
	mod->g = ff_poly_alloc (mod->d);
	mod->rgi = ff_poly_alloc(mod->e);
	mod->w = ff_poly_alloc(max);
	mod->w1 = ff_poly_alloc(max);
	ff_poly_monic (mod->g, 0, g, d_g);
	ff_poly_reverse (mod->w, &dw, mod->g, mod->d);
	ff_poly_inverse_mod_xn (mod->rgi, &d, mod->w, dw, mod->e+1, mod->w1);	// compute the inverse of rev_d(g) mod x^{e+1}
	for ( i = d+1 ; i <= mod->e ; i++ ) _ff_set_zero(mod->rgi[i]);					// zero-pad out to expected degree e
}

void ff_poly_mod_clear (ff_poly_modulus_t mod) 
	{ ff_poly_free (mod->g, mod->d); ff_poly_free (mod->rgi, mod->e); ff_poly_free (mod->w, mod->d+mod->e); ff_poly_free (mod->w1, mod->d+mod->e);  mod->g = 0; }

// h and f may coincide
void ff_poly_mod_reduce (ff_t *h, int *d_h, ff_t *f, int d_f, ff_poly_modulus_t mod)
{
	register ff_t t0;
	register int i, n;
	int d;
	
//printf ("Reducing "); ff_poly_print(f,d_f);
//printf ("Modulo "); ff_poly_print(mod->g,mod->d);
	
	// if f is already reduced, just copy it
	if ( d_f < mod->d ) { *d_h = ff_poly_degree(f,d_f); if ( h != f ) for ( i = *d_h ; i >= 0 ; i-- ) _ff_set(h[i],f[i]);  return; }

	n = mod->d+mod->e;
	assert ( d_f <= n );
	
	// compute rev_n(f), zero-padding f as required
//printf ("f = "); ff_poly_print(f,d_f);
	for ( i = 0 ; i <= n ; i++ ) if ( n-i <= d_f ) _ff_set(mod->w[i],f[n-i]); else _ff_set_zero(mod->w[i]);
//printf ("rev_%d(f) = ", n); ff_poly_print(mod->w,n);
//printf ("(rev_%d(g))^-1 mod x^{%d+1} = ", mod->d, mod->e); ff_poly_print(mod->rgi, mod->e);
	
	// compute q* = rev_n(f) * (rev_d g)^-1 (which we only need mod x^{e+1})
	ff_poly_mult (mod->w,&d,mod->rgi,mod->e,mod->w,mod->e);
//printf ("rev_%d(f) * (rev_%d(g))^-1 mod x^{%d+1} = ", n, mod->d, mod->e); ff_poly_print(mod->w,mod->e);

	// compute the quotient q = rev_e(q*)
	for ( d++ ; d <= mod->e ; d++ ) _ff_set_zero(mod->w[d]);
	for ( i = 0 ; i <= mod->e/2 ; i++ ) { _ff_set(t0, mod->w[i]); _ff_set(mod->w[i],mod->w[mod->e-i]); _ff_set(mod->w[mod->e-i],t0); }
//printf ("q = "); ff_poly_print(mod->w, mod->e);

	// compute the remainder as f - qg, which we know has degree < mod->d
	n = ( mod->e < mod->d ? mod->e : mod->d-1 );
	ff_poly_mult(mod->w,&d,mod->w,n,mod->g,mod->d-1); 
	while ( d < mod->d-1 ) { d++; _ff_set_zero(mod->w[d]); }		// zero pad out to expected degree for subtraction
	for ( i = 0 ; i < mod->d ; i++ ) ff_sub(h[i],f[i],mod->w[i]);
	*d_h = ff_poly_degree(h,mod->d-1);
//printf("f-qg = "); ff_poly_print(h,*d_h);
}

void ff_poly_mod_big (ff_t *h, int *d_h, ff_t *f, int d_f, ff_t *g, int d_g)
{
	ff_poly_modulus_t mod;
	ff_poly_mod_setup_max(mod,g,d_g,d_f);
	ff_poly_mod_reduce(h,d_h,f,d_f,mod);
	ff_poly_mod_clear(mod);
}

// f and g can point to the same place
void ff_poly_antiderivative (ff_t *g, int *pd_g, ff_t *f, int d_f)
{
	int i;
	
	if ( d_f < 0 ) { *pd_g = -1; return; }						// antiderivative of 0 is 0
	if ( d_f >= FF_POLY_MAX_DEGREE ) { err_printf ("Exceeded FF_POLY_MAX_DEGREE=%d in ff_poly_antiderivative\n", FF_POLY_MAX_DEGREE); abort(); }
	if ( d_f >= _ff_p ) { err_printf ("Degree %d must be less than field characteristic %ld in ff_poly_antiderivative\n", d_f, _ff_p); abort(); }
	ff_setup_fractions();
	for ( i = d_f+1 ; i >= 2 ; i-- ) _ff_mult(g[i],_ff_frac[i],f[i-1]);		// work from high degree to low, in case f and g overlap
	_ff_set(g[1],f[0]);  _ff_set_zero(g[0]); 
	*pd_g = d_f+1;
	return;
}


// g(x)=f(x+c), works in place, currently d_f must be <= FF_BINOMIAL_MAX
void ff_poly_translate (ff_t g[], int *pd_g, ff_t f[], int d_f, ff_t c)
{
	ff_t t, b[FF_BINOMIAL_MAX+1], cp[FF_BINOMIAL_MAX+1];
	register int i, j;
	
	if ( d_f > FF_BINOMIAL_MAX ) { err_printf ("Exceeded FF_BINOMIAL_MAX in ff_poly_translate (TODO: remove this restriction)\n");  abort(); }
	ff_binomials (d_f);
	_ff_set(cp[1],c);
	for ( i = 2 ; i <= d_f ; i++ ) _ff_mult(cp[i],cp[i-1],cp[1]);
	for ( i = 0 ; i < d_f ; i++ ) {
		for ( j = i+1 ; j <= d_f ; j++ ) _ff_mult(b[j], _ff_binomial(j,i), cp[j-i]);
		ff_dot_product(&t,f+i+1,b+i+1,d_f-i);
		_ff_add (g[i],f[i],t);
	}
	_ff_set (g[d_f],f[d_f]);
}

/*
	Given polys a, b, and c, solves the first-order linear differential equation af'+bf = c mod x^n,
	using the initial condition f(0)=0.  Based on the algorithm of Brent & Kung in Section 2.3 of BSMS2008
	f can overlap any of the inputs.
*/
void ff_poly_fold_mod_xn (ff_t *f, int *pd_f, ff_t *a, int d_a, ff_t *b, int d_b, ff_t *c, int d_c, int n)
{
	ff_t t[2*FF_POLY_MAX_DEGREE+1], h[2*FF_POLY_MAX_DEGREE+1], u[FF_POLY_MAX_DEGREE+1], w[2*FF_POLY_MAX_DEGREE+1];
	int d_t, d_h, d_u;
	
	
	if ( d_a < 0 || _ff_zero(a[0]) ) { err_printf ("Error, a(0)=0 in ff_poly_fold_mod_xn\n"); ff_poly_print(a,d_a); abort(); }
	if ( n < 2 ) { *pd_f = -1; return; }	// intial condition f(0)=0 implies f = 0 mod x
//printf ("Solving af'+bf=c mod x^%d\n",n);
//printf ("    a = "); ff_poly_print(a,d_a);
//printf ("    b = "); ff_poly_print(b,d_b);
//printf ("    c = "); ff_poly_print(c,d_c);
	
	ff_poly_inverse_mod_xn (t, &d_t, a, d_a, n-1,h);	// compute 1/a mod x^{n-1}
//printf ("1/a = "); ff_poly_print(t,d_t);
	ff_poly_mod_xn(b,&d_b,n-1);					// note that this doesn't modify b, only our local copy of d_b
	ff_poly_mod_xn(c,&d_c,n-1);					// ditto
	ff_poly_mult(h,&d_h,t,d_t,b,d_b);				// compute B = b/a mod x^{n-1}
	ff_poly_mod_xn(h,&d_h,n-1);
//printf ("B=b/a= "); ff_poly_print(h,d_h);
	ff_poly_antiderivative (h, &d_h, h, d_h);			// compute Int(B)
//printf ("Int(B)= "); ff_poly_print(h,d_h);
	ff_poly_exp_mod_xn (u, &d_u, h, d_h, n);			// compute J=exp_n(Int(B))
//printf ("J=Exp(Int(B))= "); ff_poly_print(u,d_u);
	ff_poly_mult(h,&d_h,t,d_t,c,d_c);				// compute C = c/a mod x^{n-1}
	ff_poly_mod_xn(h,&d_h,n-1);
//printf ("C=c/a= "); ff_poly_print(h,d_h);
	ff_poly_mult(t,&d_t,u,d_u,h,d_h);				// compute CJ 
	ff_poly_mod_xn(t,&d_t,n-1);
//printf ("CJ= "); ff_poly_print(t,d_t);
	ff_poly_antiderivative(t,&d_t,t,d_t);				// compute Int(CJ) mod x^n
//printf ("Int(CJ)= "); ff_poly_print(t,d_t);
	ff_poly_inverse_mod_xn(h,&d_h,u,d_u,n,w);		// compute 1/J mod x^n
//printf ("1/J= "); ff_poly_print(u,d_u);
	ff_poly_mult(h,&d_h,h,d_h,t,d_t);				// compute 1/J * Int(CJ) mod x^n (note \alpha = 0)
	ff_poly_mod_xn(h,&d_h,n);
	ff_poly_copy(f,pd_f,h,d_h);					// copy into f
//printf ("f=Int(CJ)/J= "); ff_poly_print(f,*pd_f);
}

/*
	Solves the first-order non-linear differential equation S'^2 = C(1+A*S^4+B*S^6) mod x^n
	using the initial conditions S(0)=0 and S'(0) = 1.  Here C is a poly and A and B are scalars.
	Based on 2.4 of BSMS2008

	There is a bug in this code (it occasionally gets the highest order coefficient wrong).
*/
void ff_poly_S_mod_xn (ff_t *S, int *pd_S, ff_t A, ff_t B, ff_t *C, int d_C, int n)
{
	// make static to avoid worrying about stack overflow - this is a temporary hack
	static ff_t f[FF_POLY_MAX_DEGREE+1], f2[2*FF_POLY_MAX_DEGREE+1], f3[2*FF_POLY_MAX_DEGREE+1];
	static ff_t u[2*FF_POLY_MAX_DEGREE+1], a[2*FF_POLY_MAX_DEGREE+1], b[2*FF_POLY_MAX_DEGREE+1], c[2*FF_POLY_MAX_DEGREE+1], h[2*FF_POLY_MAX_DEGREE+1];
	ff_t A4,B6,t0;
	int d_f, d_f2, d_f3, d_h, d_u, d_a, d_b, d_c, s, m, d, i;

//printf("A=%ld, B=%ld, C= ", _ff_get_ui(A), _ff_get_ui(B)); ff_poly_print(C,d_C);
	
	_ff_set_zero(f[0]); _ff_set_one(f[1]);	d_f=1;						// f = x mod x^2
	_ff_add(A4,A,A); _ff_x2(A4);
	_ff_add(t0,B,B); _ff_triple(B6,t0);
	s = 5;													// we can always skip the first two iterations, we know that f = x+O(x^5)
	while ( s < n ) {
//printf ("f: "); ff_poly_print(f,d_f);
		m = 2*s-2;
		ff_poly_square(f2,f,d_f);  d_f2 = 2*d_f;
		ff_poly_mod_xn(f2,&d_f2,m);								// f2 = f^2
//printf ("f^2: "); ff_poly_print(f2,d_f2);
		ff_poly_mult(f3,&d_f3,f,d_f,f2,d_f2);
		ff_poly_mod_xn(f3,&d_f3,m);								// f3 = f^3
//printf ("f^3: "); ff_poly_print(f3,d_f3);
		for ( i = 0 ; i <= d_f3 ; i++ ) _ff_mult(u[i],A4,f3[i]);  d_u = d_f3;	// u = 4Af^3
//printf ("4*A*f3: "); ff_poly_print(u,d_u);
		ff_poly_mult(h,&d_h,f2,d_f2,f3,d_f3);
		ff_poly_mod_xn(h,&d_h,m);								// h = f^5
//printf("f^5: "); ff_poly_print(h,d_h);
		for ( i = 0 ; i <= d_h ; i++ ) _ff_mult(h[i],B6,h[i]);				// h = 6Bf^5
		ff_poly_addto (u,&d_u,h,d_h);								// u = 4Af^3+6Bf^5
//printf ("4*A*f3+6Bf^5: "); ff_poly_print(u,d_u);
		if ( d_C > m-1 ) d = m-1; else d = d_C;
		ff_poly_mult(b,&d_b,C,d,u,d_u);							// b = C(4f^3+6Bf^5)
		ff_poly_mod_xn(b,&d_b,m);
		for ( i = 0 ; i <= d_b ; i++ ) ff_negate(b[i]);					// b = -C(4f^3+6Bf^5)
//printf("b=-C(4f^3+6Bf^5): "); ff_poly_print(b,d_b);
		ff_poly_square(u,f2,d_f2);  d_u=2*d_f2;						// u = f^4
		ff_poly_mod_xn(u,&d_u,m);
		for  ( i = 0 ; i <= d_u ; i++ ) ff_mult(u[i],A,u[i]);
//printf ("A*f^4: "); ff_poly_print(u,d_u);
		ff_poly_square(h,f3,d_f3); d_h = 2*d_f3;						// h = f^6
		ff_poly_mod_xn(h,&d_h,m);
		for  ( i = 0 ; i <= d_h ; i++ ) ff_mult(h[i],B,h[i]);
//printf ("B*f^6: "); ff_poly_print(h,d_h);
		ff_poly_addto(u,&d_u,h,d_h);								// u = f^4+f^6
//printf("d_u=%d\n", d_u);
		if ( d_u<0 ) {_ff_set_one(u[0]); d_u=0; } else _ff_inc(u[0]);		// u = 1+f^4+f^6
//printf ("1+A*f^4+B*f^6: "); ff_poly_print(u,d_u);
		ff_poly_mult(c,&d_c,C,d_C,u,d_u);
		ff_poly_mod_xn(c,&d_c,m);								// c = C(1+f^4+f^6)
//printf ("C*(1+A*f^4+B*f^6): "); ff_poly_print(c,d_c);
		ff_poly_derivative(a,&d_a,f,d_f);							// a = f'
		ff_poly_square(h,a,d_a);  d_h = 2*d_a;						// h = f'^2
		ff_poly_mod_xn(h,&d_h,m);
		ff_poly_subfrom(c,&d_c,h,d_h);							// c = C(1+f^4+f^2) - f'^2
//printf ("c=C*(1+A*f^4+B*f^6)-f'^2: "); ff_poly_print(c,d_c);
		for ( i = 0 ; i <= d_a ; i++ ) _ff_x2(a[i]);						// a = 2f'
//printf ("a=2f': "); ff_poly_print(a,d_a);
		ff_poly_fold_mod_xn(f2,&d_f2,a,d_a,b,d_b,c,d_c,m);
//printf ("f_2: "); ff_poly_print(f2,d_f2);
		ff_poly_addto(f,&d_f,f2,d_f2);								// f = f_1 + f_2
		s = 2*s-1;
	}
	ff_poly_mod_xn(f,&d_f,n);
//printf ("final f: "); ff_poly_print(f,d_f);
	ff_poly_copy(S,pd_S,f,d_f);
}
