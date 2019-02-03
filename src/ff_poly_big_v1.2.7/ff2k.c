#include <assert.h>
#include <stdio.h>
#include "cstd.h"
#include "ff2k.h"

/*
    Copyright 2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/


// primitive polys taken from Hansen & Mullen Supplement in Math Comp Vol 59, No 200, 1992, pp. S47-S50
static ff2k_t ff2k_xkrtab[FF2K_MAXK+1] =
{ 0, 0, 0x3, 0x3, 0x3, 0x5, 0x3, 0x3, 0x1d, 0x11, 0x9, 0x5, 0x53, 0x1b, 0x2b, 0x3, 0x2d,
  0x9, 0x81, 0x27, 0x9, 0x5, 0x3, 0x23, 0x1b,  0x9, 0x47, 0x27, 0x9, 0x5, 0x53, 0x5 };
ff2k_t _ff2k_xk, _ff2k_xkr, _ff2k_x2km2, _ff2k_x2km2r, _ff2k_m;
int _ff2k_k;

void ff2k_setup (int k)
{
	assert ( k >= 1 );
	assert ( k <= FF2K_MAXK );
	_ff2k_k = k;
	_ff2k_xk = 1UL << k;
	_ff2k_m = _ff2k_xk-1;
	_ff2k_xkr = ff2k_xkrtab[k];
	_ff2k_x2km2 = 1UL << (2*k-2);
	_ff2k_x2km2r = _ff2k_xkr << (k-2);
}

void ff2k_poly_print (ff2k_t f[], int d)
{
	register int i;
	
	printf ("[");
	for ( i = d ; i >= 0 ; i-- ) {
		if ( f[i] ) { 
			printf ("+ "); 
			switch (i) { case 0: ff2k_print(f[i]); break; case 1: ff2k_print(f[i]); printf ("*X"); break; default: ff2k_print (f[i]); printf ("*X^%d", i); }
		}
	}
	puts ("]");
}


// pure brute-force root-finding, only for use with very small k (e.g. <= 3)
int ff2k_poly_distinct_roots (ff2k_t r[], ff2k_t f[], int d)
{
	ff2k_t x, t;
	register int i;
	
	_ff2k_set_zero(x);
	i = 0;
	do {
		ff2k_poly_eval (&t, f, d, &x);
		if ( _ff2k_zero(t) ) { _ff2k_set(r[i],x); i++; }
		_ff2k_next(x);
	} while ( ! _ff2k_zero(x) );
	return i;
}

// determine whether the elliptic curve with Weierstrass coefficients w = [a1,a2,a3,a4,a6] is non-singular
int ff2k_WS_nonsingular (ff2k_t w[5])
{
	register ff2k_t a1, a3, t1, t2, t3;
	
	if ( _ff2k_zero(w[0]) ) return ( _ff2k_zero(w[2]) ? 0 : 1);
	
	// compute D = a1^6*a6 + a1^5*a3*a4 + a1^4*a2*a3^2 + a1^4*a4^2 + a1^3*a3^3 + a3^4
	//                   = a1*(a1*(a1*(a1*(a1*(a1*a6+a3*a4)+a2*a3^2+a4^2)+a3^3))) + a3^4
	_ff2k_set (a1,w[0]); _ff2k_set (a3,w[2]);
	_ff2k_mult (t1,a1,w[4]); _ff2k_mult (t2,a3,w[3]); _ff2k_addto (t1,t2);
	_ff2k_mult (t1,t1,a1); _ff2k_square (t3,a3); _ff2k_mult(t2,w[1],t3); _ff2k_addto (t1,t2); _ff2k_square (t2, w[3]); _ff2k_addto (t1, t2);
	_ff2k_mult (t1,t1,a1); _ff2k_mult (t3,t3,a3); _ff2k_addto (t1,t3);
	_ff2k_square (t2, a1); _ff2k_mult (t2,t2,a1); _ff2k_mult (t1, t1, t2); _ff2k_mult(t3,t3,a3); _ff2k_addto (t1,t3);
	return ( _ff2k_zero (t1) ? 0 : 1 );
}

// very naive point-counting for elliptic curves over F_2^k, intended for small k (e.g. 1 or 2)
long ff2k_WS_pointcount (ff2k_t w[])
{
	register ff2k_t x, y, t1, t2;
	long pts;

	if ( ! ff2k_WS_nonsingular (w) ) return -1;
	pts = 0;
	_ff2k_set_zero(y);
	do {
		_ff2k_set_zero (x);
		do {
			_ff2k_square (t1, y);  _ff2k_mult (t2, w[0], y); _ff2k_mult (t2, t2, x); _ff2k_addto (t1, t2); _ff2k_mult (t2, w[2], y); _ff2k_addto (t1, t2);
			_ff2k_add (t2, x, w[1]); _ff2k_mult (t2, t2, x); _ff2k_addto (t2, w[3]); _ff2k_mult (t2, t2, x); _ff2k_addto (t2, w[4]);
			if ( _ff2k_equal (t1, t2) ) pts++;
			_ff2k_next(x);
		} while ( ! _ff2k_zero (x) );
		_ff2k_next(y);
	} while ( ! _ff2k_zero(y) );
//printf ("pointcount over F_%ld is %ld\n", (1L<<_ff2k_k), pts+1);
	return pts+1;	// always one point at infinity
}

// counts points on y^2+h(x)*y=f(x) over current finite field with 2^k elements
// does not check for good reduction
long ff2k_hyperelliptic_pointcount (ff2k_t f[], int df, ff2k_t h[], int dh)
{
	ff2k_t t1, t2, t3, x, y;
	long pts;
	
	pts = 0;
	_ff2k_set_zero (y);
	do {
		_ff2k_set_zero (x);
		do {
			ff2k_poly_eval (&t1, f, df, &x);
			ff2k_poly_eval (&t2, h, dh, &x);
			_ff2k_mult (t3, t2, y);
			_ff2k_square (t2, y);
			_ff2k_addto (t2, t3);
			if ( _ff2k_equal (t1, t2) ) pts++;
			_ff2k_next (x);
		} while ( ! _ff2k_zero (x) );
		_ff2k_next(y);
	} while ( ! _ff2k_zero (y) );
	return pts+ff2k_hyperelliptic_pointcount_at_infinity (df, dh, _ff2k_k);
}



