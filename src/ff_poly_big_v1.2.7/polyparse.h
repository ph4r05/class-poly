#ifndef _POLYPARSE_INCLUDE_
#define _POLYPARSE_INCLUDE_

/*
    Copyright 2011-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
	The generic parsing functions below support multiple data types via the helper functions setzero, addto, and iszero, which operate on  a single coefficient f[i].
	f points to an array of maxd+1 coefficients, which may be of any datatype.  Note that addto is passed a rational number, and may return 0 if it can't be converted to the desired datatype.

	Note that the addto function is allowed to trash c if it likes (so it can use it as workspace)

	This isn't superfast, but it saves repeating the same code many times, and parsing poly strings is generally not a limiting performance issue.
	
	Returns the degree of the poly, with degree -1 for the zero polynomial.
	Returns -2 if the poly degree exceeds maxd, and returns -3 if an error is encountered.
*/
// poly_parse expects a univariate polynomial
int poly_parse (void *f, int maxd, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg);
// poly_parse_plane_quartic expects a plane quartic in affine or project coords, f should point to an array of 15 coefficients stored in lexicographic order (x^0y^0z^4, x^0y^1z^3, ..., x^4y^0z^0)
int poly_parse_plane_quartic (void *f, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg);

// returns the maximum degree of poly specified by expr (the actual degree could be less due to zero coefficients (e.g. after reduction mod p)).
int poly_parse_max_degree (char *expr);

// parses polynomial expression with rational coefficients into an array of unsigned longs reduced mod p -- returns 0 if any coefficient contains an univertible denominator -- implemented via poly_parse
int ui_poly_parse_mod_p (unsigned long f[], int maxd, char *expr, unsigned long p);
int ui_poly_parse_plane_quartic_mod_p (unsigned long f[], char *expr, unsigned long p);

// parses polynomial expression with integer coefficients into an array of longs -- returns 0 if any coefficient is non-integral, or in the case of overflow
int i_poly_parse (long f[], int maxd, char *expr);
int i_poly_parse_plane_quartic (long f[], char *expr);

int ui_poly_sprint (char *s, unsigned long f[], int df);
void ui_poly_print (unsigned long f[], int df);
int i_poly_sprint (char *s, long f[], int df);
void i_poly_print (long f[], int df);

// no overflow check
static inline long i_poly_eval (long f[], int df, long x)
{
	register long y, *c;
	
	if ( df < 0 ) return 0;
	c = f+df; y = *c;
	while ( c > f ) y = x*y+*(--c);
	return y;
}

// Provided x is reduced mod p < 2^31 and coeffs of f are bounded by LONG_MAX/2, overflow should not be a problem
static inline long i_poly_eval_mod_p (long f[], int df, long x, long p)
{
	register long y, *c;
	
	if ( df < 0 ) return 0;
	c = f+df; y = *c % p;
	while ( c > f ) y = (x*y+*(--c))%p;
	if ( y < 0 ) y += p;
	return y;
}

// parses hyperelliptic curve spec in the form [f] or [f,h] defining the curve y^2 + h(x)y = f(x).  Returns the genus, or 0 for failure
static inline int i_hyperelliptic_curve_parse (long f[], int *pdf, long h[], int *pdh, int maxd, char *str)		// maxd bounds the degree of f, (maxd+1)/2 bounds the degree of h, returns 0/1 for failure/success
{
	char *s;
	int d;
	
	*pdf = i_poly_parse (f, maxd, str);
	if ( *pdf < 0 ) return 0;	// note that f must by nonzero (o.w. curve would be singular)
	for ( s = str ; *s && *s != ',' ; s++ );
	if ( ! s ) { *pdh = -1; return (*pdf-1)/2; }
	*pdh = i_poly_parse (h, (maxd+1)/2, s+1);
	if ( *pdh < -1 ) return 0;
	d = 2*(*pdh) > (*pdf) ? 2*(*pdh) : *pdf;
	return (d-1)/2;
}

static inline int i_hyperelliptic_curve_normalize (long g[], long f[], int df, long h[], int dh)				// sets g to 4*f+h^2, or just to f if h is zero (g is allowed to alias f but not h)
{
	register int i, j;
	
	if ( dh < 0 ) { if ( g != f ) for ( i = 0 ; i <= df ; i++ ) g[i] = f[i];  return df; }
	for ( i = 0 ; i <= df ; i++ ) g[i] = 4*f[i];
	while ( i <= 2*dh ) g[i++] = 0;
	for ( i = 0 ; i <= dh ; i++ ) for ( j = 0 ; j <= dh ; j++ ) g[i+j] += h[i]*h[j];
	if ( 2*dh > df ) df = 2*dh;
	for ( i = 0 ; i <= df ; i++ ) if ( (f[i]&3) ) break;
	if ( i > df ) for ( i = 0 ; i <= df ; i++ ) f[i] /= 4;
	return df;
}

static inline int i_hyperelliptic_curve_sprint (char *s, long f[], int df, long h[], int dh)
{
	char *t = s + i_poly_sprint (s, f, df) ;
	char *u = t + i_poly_sprint (t, h, dh);
	*(t-1) = ','; *t = ' ';
	return u-s;
}

#ifdef __cplusplus
}
#endif

#endif
