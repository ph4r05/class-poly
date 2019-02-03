#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "ffext.h"
#include "ffpoly.h"
#include "ffpolysmall.h"
#include "cstd.h"

/*
    Copyright 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

/*
	This module contains performance-oriented code for operations on polynomials of low degree,
	and support functions for computations on genus 1 and genus 2 curves.
	
	It is assumed here that the ff_t type does not require initialization (no support for GMP's mpz type),
	and in a number of cases it is further assumed that the ff_t type fits in an unsigned long (this assumption is verified when made).
*/

/*
	Generic small poly mult/square, used for small degree
	inputs/outputs may be the same array
*/

void ff_poly_mult_small (ff_t o[], ff_t f[], ff_t g[], int n)			// polys must be same size, n=deg, not # coeff
{
	register int i;
	
	for ( i = n ; i >=0 ; i-- ) ff_conv(o+n+i,f+i,g+i,n-i+1);
	for ( i = n-1 ; i >= 0 ; i-- ) ff_conv(o+i,f,g,i+1);
}

void ff_poly_square_small_old (ff_t o[], ff_t f[], int n)			// n = deg, not # coeff
{
	register int i;

	for ( i = n ; i >=0 ; i-- ) ff_unary_conv(o+n+i,f+i,n-i+1);
	for ( i = n-1 ; i >= 0 ; i-- ) ff_unary_conv(o+i,f,i+1);
}

void ff_poly_square_small_big (ff_t o[], ff_t f[], int n)			// n = deg, not # coeff, n must be at least 20
{
	register ff_t *op, *fp;
	register int i;

	if ( n < 0 ) return;
	if ( n < 20 ) { printf ("Call to ff_poly_square_small_big with n=%d < 20\n", n); abort(); }
	op = o+2*n; fp = f+n;
	_ff_square(*op,*fp);  op--; fp--;
	_ff_sum_2_mults_s2_arr (*op,fp);  op--; fp--;
	_ff_sum_3_mults_s3_arr (*op,fp);  op--; fp--;
	_ff_sum_4_mults_s4_arr (*op,fp);  op--; fp--;
	_ff_sum_5_mults_s5_arr (*op,fp);  op--; fp--;
	_ff_sum_6_mults_s6_arr (*op,fp);  op--; fp--;
	_ff_sum_7_mults_s7_arr (*op,fp);  op--; fp--;
	_ff_sum_8_mults_s8_arr (*op,fp);  op--; fp--;
	_ff_sum_9_mults_s9_arr (*op,fp);  op--; fp--;
	_ff_sum_10_mults_s10_arr (*op,fp);  op--; fp--;
	_ff_sum_11_mults_s11_arr (*op,fp);  op--; fp--;
	_ff_sum_12_mults_s12_arr (*op,fp);  op--; fp--;
	_ff_sum_13_mults_s13_arr (*op,fp);  op--; fp--;
	_ff_sum_14_mults_s14_arr (*op,fp);  op--; fp--;
	_ff_sum_15_mults_s15_arr (*op,fp);  op--; fp--;
	_ff_sum_16_mults_s16_arr (*op,fp);  op--; fp--;
	_ff_sum_17_mults_s17_arr (*op,fp);  op--; fp--;
	_ff_sum_18_mults_s18_arr (*op,fp);  op--; fp--;
	_ff_sum_19_mults_s19_arr (*op,fp);  op--; fp--;
	_ff_sum_20_mults_s20_arr (*op,fp);  op--; fp--;
	for ( i = n-20 ; i >=0 ; i-- ) ff_unary_conv(o+n+i,f+i,n-i+1);
	for ( i = n-1 ; i >= 20 ; i-- ) ff_unary_conv(o+i,f,i+1);
	_ff_sum_20_mults_s20_arr (o[19],f);
	_ff_sum_19_mults_s19_arr (o[18],f);
	_ff_sum_18_mults_s18_arr (o[17],f);
	_ff_sum_17_mults_s17_arr (o[16],f);
	_ff_sum_16_mults_s16_arr (o[15],f);
	_ff_sum_15_mults_s15_arr (o[14],f);
	_ff_sum_14_mults_s14_arr (o[13],f);
	_ff_sum_13_mults_s13_arr (o[12],f);
	_ff_sum_12_mults_s12_arr (o[11],f);
	_ff_sum_11_mults_s11_arr (o[10],f);
	_ff_sum_10_mults_s10_arr (o[9],f);
	_ff_sum_9_mults_s9_arr (o[8],f);
	_ff_sum_8_mults_s8_arr (o[7],f);
	_ff_sum_7_mults_s7_arr (o[6],f);
	_ff_sum_6_mults_s6_arr (o[5],f);
	_ff_sum_5_mults_s5_arr (o[4],f);
	_ff_sum_4_mults_s4_arr (o[3],f);
	_ff_sum_3_mults_s3_arr (o[2],f);
	_ff_sum_2_mults_s2_arr (o[1],f);
	_ff_square(o[0],f[0]);
}

void ff_poly_square_small (ff_t g[], ff_t f[], int d_f)
{
	switch (d_f) {
	case 0: _ff_square(g[0],f[0]); return;
	case 1: ff_poly_square_2(g,f); return;
	case 2: ff_poly_square_3(g,f); return;
	case 3: ff_poly_square_4(g,f); return;
	case 4: ff_poly_square_5(g,f); return;
	case 5: ff_poly_square_6(g,f); return;
	case 6: ff_poly_square_7(g,f); return;
	case 7: ff_poly_square_8(g,f); return;
	case 8: ff_poly_square_9(g,f); return;
	case 9: ff_poly_square_10(g,f); return;
	case 10: ff_poly_square_11(g,f); return;
	case 11: ff_poly_square_12(g,f); return;
	case 12: ff_poly_square_13(g,f); return;
	case 13: ff_poly_square_14(g,f); return;
	case 14: ff_poly_square_15(g,f); return;
	case 15: ff_poly_square_16(g,f); return;
	case 16: ff_poly_square_17(g,f); return;
	case 17: ff_poly_square_18(g,f); return;
	case 18: ff_poly_square_19(g,f); return;
	case 19: ff_poly_square_20(g,f); return;
	}
	ff_poly_square_small_big(g,f,d_f); return;
}

// Replaces f by f mod g, with g monic (not verified), result is zero-filled to degree d_g-1
void ff_poly_mod_small_inplace (ff_t *f, int d_f, ff_t *g, int d_g)
{
	register int i, m;
	register ff_t t0;
	ff_t t;

//	printf ("Computing "); poly_print(f,d_f); printf (" mod "); poly_print(g,d_g); puts("");
	if ( d_f < d_g ) return;
	switch (d_g) {
		case 0: return;
		case 1: _ff_neg(t,g[0]); ff_poly_eval(f,f,d_f,&t); return;
		case 2: ff_poly_mod_d2 (f, d_f, g); return;
		case 3: ff_poly_mod_d3 (f, d_f, g); return;
		case 4: ff_poly_mod_d4 (f, d_f, g); return;
		case 5: ff_poly_mod_d5 (f, d_f, g); return;
		case 6: ff_poly_mod_d6 (f, d_f, g); return;
	}
	switch (d_f-d_g) {
	case 0: for ( i = 0 ; i < d_g ; i++ ) { _ff_mult(t0,f[d_g],g[i]); _ff_subfrom(f[i],t0); } return;
	case 1: ff_poly_mod_delta1 (f, g, d_g); return;
	case 2: ff_poly_mod_delta2 (f, g, d_g); return;
	case 3: ff_poly_mod_delta3 (f, g, d_g); return;
	case 4: ff_poly_mod_delta4 (f, g, d_g); return;
	case 5: ff_poly_mod_delta5 (f, g, d_g); return;
	}
	
	m = d_f-d_g+1;  if ( d_g < m ) m = d_g;
	
	// we know that m>=6, so we unroll the first six interations of the top loop and the last 6 of the bottom loop
	_ff_mult(t0,f[d_f],g[d_g-1]); _ff_subfrom(f[d_f-1],t0);
	_ff_sum_2_mults_arr(t0,f+d_f-1,g+d_g-2); _ff_subfrom(f[d_f-2],t0);
	_ff_sum_3_mults_arr(t0,f+d_f-2,g+d_g-3); _ff_subfrom(f[d_f-3],t0);
	_ff_sum_4_mults_arr(t0,f+d_f-3,g+d_g-4); _ff_subfrom(f[d_f-4],t0);
	_ff_sum_5_mults_arr(t0,f+d_f-4,g+d_g-5); _ff_subfrom(f[d_f-5],t0);
	_ff_sum_6_mults_arr(t0,f+d_f-5,g+d_g-6); _ff_subfrom(f[d_f-6],t0);
	for ( i = d_f-7 ; i >= d_f-m ; i-- ) {/*printf ("1) coeff %d, conv %d goff %d\n", i, d_f-i, d_g-d_f+i);*/ ff_conv(&t,f+i+1,g+d_g-d_f+i,d_f-i); _ff_subfrom(f[i],t); }
	
	// only one of these two loops will be used
	for ( ; i >= d_g ; i-- ) {/* printf ("2) coeff %d, conv %d goff 0\n", i, d_g);*/ ff_conv(&t,f+i+1,g,d_g); _ff_subfrom(f[i],t); }
	for ( ; i >= m ; i-- ) { /*printf ("3) coeff %d, conv %d goff %d\n", i, m,i-m+1);*/  ff_conv(&t,f+d_g,g+i-m+1,m); _ff_subfrom(f[i],t); }
	
	for ( ; i >= 6 ; i-- ) { /*printf("4) coeff %d, conv %d goff 0\n", i, i+1);*/ ff_conv(&t,f+d_g,g,i+1); _ff_subfrom(f[i],t); }
	_ff_sum_6_mults_arr(t0,f+d_g,g); _ff_subfrom(f[5],t0);
	_ff_sum_5_mults_arr(t0,f+d_g,g); _ff_subfrom(f[4],t0);
	_ff_sum_4_mults_arr(t0,f+d_g,g); _ff_subfrom(f[3],t0);
	_ff_sum_3_mults_arr(t0,f+d_g,g); _ff_subfrom(f[2],t0);
	_ff_sum_2_mults_arr(t0,f+d_g,g); _ff_subfrom(f[1],t0);
	_ff_mult(t0,f[d_g],g[0]); _ff_subfrom(f[0],t0);

//printf ("Result "); poly_print(f,d_g-1); puts("");
}

void ff_poly_mod_small (ff_t *h, int *d_h, ff_t *f, int d_f, ff_t *g, int d_g)
{
	ff_t s[FF_POLY_SMALL_DEGREE+1], t[FF_POLY_SMALL_DEGREE+1], *fp, *gp;
	int d;

	if ( d_f < d_g ) { if ( h != f ) ff_poly_copy (h, d_h, f, d_f);  return; }
	
	// prepare to do modular reduction in place, copy f if it needs to be preserved
	if ( h != f ) {
		assert ( d_f <= FF_POLY_SMALL_DEGREE );
		ff_poly_copy (s, &d, f, d_f);
		fp = s;
	} else {
		fp = f;
	}
	
	// if g is not monic, copy it and make it monic
	if ( ! _ff_one(g[d_g]) ) {
		ff_poly_monic (t, &d, g, d_g);
		gp = t;
	} else {
		gp = g;
	}
	ff_poly_mod_small_inplace(fp, d_f, gp, d_g);
	*d_h = ff_poly_degree(fp,d_g-1);
	if ( h != fp ) ff_poly_copy (h, &d, fp, *d_h);
}

// computes f^2 mod g where g=x^n - g[n-2]x^{n-2} - ... - g[0] (note signs!)  and f has degree n-1 (n must be > 0)
void ff_poly_square_mod_small_old (ff_t o[], ff_t f[], ff_t g[], int n) 
{
	ff_t t[FF_POLY_SMALL_DEGREE];
	register int i;
	register ff_t t0, t1;
	
	n--;	// set n to degree of f and output
	_ff_square(t[n],f[n]);
	_ff_dbl_mult(t[n-1],f[n-1],f[n]);
	for ( i=n-2 ; i >= 0 ; i-- ) {
		_ff_sq_coeff (t0, f+i, n-i);
		_ff_sum_mults (t1, t+i+2, g+i+1, n-i-1);
		_ff_add(t[i],t0,t1);
	}
	_ff_set(o[n],t[0]);
	for ( i = n-1 ; i >= 0 ; i-- ) {
		_ff_sq_coeff (t0,f,i);
		_ff_sum_mults(t1,t+1,g,i+1);
		_ff_add(o[i],t0,t1);
	}
}

// computes f^2 mod g where g=x^n - g[n-2]x^{n-2} - ... - g[0] (note signs!)  and f has degree n-1 (n must be > 0)
// assumes n > 20
void ff_poly_square_mod_small_big  (ff_t o[], ff_t f[], ff_t g[], int n) 
{
	ff_t t[FF_POLY_SMALL_DEGREE];
	register int i;
	ff_t t0, t1;
	
	if ( n <= 20 || n > FF_POLY_SMALL_DEGREE ) { printf ("Call to ff_poly_square_mod_small_big with n=%d <= 20 or > %d\n", n, FF_POLY_SMALL_DEGREE); abort(); }
	n--;	// set n to degree of f and output
	_ff_square(t[n],f[n]);
	_ff_dbl_mult(t[n-1],f[n-1],f[n]);
	_ff_sum_3_mults_d1(t[n-2],f[n-2],f[n-1],t[n],g[n-1],f[n-1],f[n]);
	_ff_sum_4_mults_d2(t[n-3],f[n-3],f[n-2],t[n-1],t[n],g[n-2],g[n-1],f[n-1],f[n]);
	_ff_sum_6_mults_d2(t[n-4],f[n-4],f[n-3],f[n-2],t[n-2],t[n-1],t[n],g[n-3],g[n-2],g[n-1],f[n-2],f[n-1],f[n]);
	_ff_sum_7_mults_d3(t[n-5],f[n-5],f[n-4],f[n-3],t[n-3],t[n-2],t[n-1],t[n],g[n-4],g[n-3],g[n-2],g[n-1],f[n-2],f[n-1],f[n]);
	_ff_sum_9_mults_d3(t[n-6],f[n-6],f[n-5],f[n-4],f[n-3],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_10_mults_d4(t[n-7],f[n-7],f[n-6],f[n-5],f[n-4],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_12_mults_d4(t[n-8],f[n-8],f[n-7],f[n-6],f[n-5],f[n-4],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_13_mults_d5(t[n-9],f[n-9],f[n-8],f[n-7],f[n-6],f[n-5],t[n-7],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-8],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_15_mults_d5(t[n-10],f[n-10],f[n-9],f[n-8],f[n-7],f[n-6],f[n-5],t[n-8],t[n-7],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-9],g[n-8],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-5],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_16_mults_d6(t[n-11],f[n-11],f[n-10],f[n-9],f[n-8],f[n-7],f[n-6],t[n-9],t[n-8],t[n-7],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-10],g[n-9],g[n-8],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-5],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_18_mults_d6(t[n-12],f[n-12],f[n-11],f[n-10],f[n-9],f[n-8],f[n-7],f[n-6],t[n-10],t[n-9],t[n-8],t[n-7],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-11],g[n-10],g[n-9],g[n-8],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-6],f[n-5],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_19_mults_d7(t[n-13],f[n-13],f[n-12],f[n-11],f[n-10],f[n-9],f[n-8],f[n-7],t[n-11],t[n-10],t[n-9],t[n-8],t[n-7],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-12],g[n-11],g[n-10],g[n-9],g[n-8],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-6],f[n-5],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_21_mults_d7(t[n-14],f[n-14],f[n-13],f[n-12],f[n-11],f[n-10],f[n-9],f[n-8],f[n-7],t[n-12],t[n-11],t[n-10],t[n-9],t[n-8],t[n-7],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-13],g[n-12],g[n-11],g[n-10],g[n-9],g[n-8],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-7],f[n-6],f[n-5],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_22_mults_d8(t[n-15],f[n-15],f[n-14],f[n-13],f[n-12],f[n-11],f[n-10],f[n-9],f[n-8],t[n-13],t[n-12],t[n-11],t[n-10],t[n-9],t[n-8],t[n-7],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-14],g[n-13],g[n-12],g[n-11],g[n-10],g[n-9],g[n-8],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-7],f[n-6],f[n-5],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_24_mults_d8(t[n-16],f[n-16],f[n-15],f[n-14],f[n-13],f[n-12],f[n-11],f[n-10],f[n-9],f[n-8],t[n-14],t[n-13],t[n-12],t[n-11],t[n-10],t[n-9],t[n-8],t[n-7],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-15],g[n-14],g[n-13],g[n-12],g[n-11],g[n-10],g[n-9],g[n-8],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-8],f[n-7],f[n-6],f[n-5],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_25_mults_d9(t[n-17],f[n-17],f[n-16],f[n-15],f[n-14],f[n-13],f[n-12],f[n-11],f[n-10],f[n-9],t[n-15],t[n-14],t[n-13],t[n-12],t[n-11],t[n-10],t[n-9],t[n-8],t[n-7],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-16],g[n-15],g[n-14],g[n-13],g[n-12],g[n-11],g[n-10],g[n-9],g[n-8],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-8],f[n-7],f[n-6],f[n-5],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_27_mults_d9(t[n-18],f[n-18],f[n-17],f[n-16],f[n-15],f[n-14],f[n-13],f[n-12],f[n-11],f[n-10],f[n-9],t[n-16],t[n-15],t[n-14],t[n-13],t[n-12],t[n-11],t[n-10],t[n-9],t[n-8],t[n-7],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-17],g[n-16],g[n-15],g[n-14],g[n-13],g[n-12],g[n-11],g[n-10],g[n-9],g[n-8],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-9],f[n-8],f[n-7],f[n-6],f[n-5],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	_ff_sum_28_mults_d10(t[n-19],f[n-19],f[n-18],f[n-17],f[n-16],f[n-15],f[n-14],f[n-13],f[n-12],f[n-11],f[n-10],t[n-17],t[n-16],t[n-15],t[n-14],t[n-13],t[n-12],t[n-11],t[n-10],t[n-9],t[n-8],t[n-7],t[n-6],t[n-5],t[n-4],t[n-3],t[n-2],t[n-1],t[n],g[n-18],g[n-17],g[n-16],g[n-15],g[n-14],g[n-13],g[n-12],g[n-11],g[n-10],g[n-9],g[n-8],g[n-7],g[n-6],g[n-5],g[n-4],g[n-3],g[n-2],g[n-1],f[n-9],f[n-8],f[n-7],f[n-6],f[n-5],f[n-4],f[n-3],f[n-2],f[n-1],f[n]);
	for ( i = n-20 ; i >=0 ; i-- ) { ff_unary_conv(&t0,f+i,n-i+1); ff_conv(&t1,t+i+2, g+i+1, n-i-1); _ff_add(t[i],t0,t1); }
	_ff_set(o[n],t[0]);
	for ( i = n-1 ; i > 19 ; i-- ) { ff_unary_conv(&t0,f,i+1); ff_conv(&t1,t+1,g,i+1);  _ff_add(o[i],t0,t1); }
	_ff_sum_30_mults_d10(o[19],f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],t[17],t[18],t[19],t[20],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],g[19],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19]);
	_ff_sum_29_mults_d9(o[18],f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],t[17],t[18],t[19],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],g[18],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
	_ff_sum_27_mults_d9(o[17],f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],t[17],t[18],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],g[17],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17]);
	_ff_sum_26_mults_d8(o[16],f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],t[17],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],g[16],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16]);
	_ff_sum_24_mults_d8(o[15],f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],g[15],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15]);
	_ff_sum_23_mults_d7(o[14],f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],g[14],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14]);
	_ff_sum_21_mults_d7(o[13],f[0],f[1],f[2],f[3],f[4],f[5],f[6],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],g[13],f[7],f[8],f[9],f[10],f[11],f[12],f[13]);
	_ff_sum_20_mults_d6(o[12],f[0],f[1],f[2],f[3],f[4],f[5],f[6],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],g[12],f[6],f[7],f[8],f[9],f[10],f[11],f[12]);
	_ff_sum_18_mults_d6(o[11],f[0],f[1],f[2],f[3],f[4],f[5],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],g[11],f[6],f[7],f[8],f[9],f[10],f[11]);
	_ff_sum_17_mults_d5(o[10],f[0],f[1],f[2],f[3],f[4],f[5],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10],f[5],f[6],f[7],f[8],f[9],f[10]);
	_ff_sum_15_mults_d5(o[9],f[0],f[1],f[2],f[3],f[4],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],f[5],f[6],f[7],f[8],f[9]);
	_ff_sum_14_mults_d4(o[8],f[0],f[1],f[2],f[3],f[4],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],f[4],f[5],f[6],f[7],f[8]);
	_ff_sum_12_mults_d4(o[7],f[0],f[1],f[2],f[3],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],f[4],f[5],f[6],f[7]);
	_ff_sum_11_mults_d3(o[6],f[0],f[1],f[2],f[3],t[1],t[2],t[3],t[4],t[5],t[6],t[7],g[0],g[1],g[2],g[3],g[4],g[5],g[6],f[3],f[4],f[5],f[6]);
	_ff_sum_9_mults_d3(o[5],f[0],f[1],f[2],t[1],t[2],t[3],t[4],t[5],t[6],g[0],g[1],g[2],g[3],g[4],g[5],f[3],f[4],f[5]);
	_ff_sum_8_mults_d2(o[4],f[0],f[1],f[2],t[1],t[2],t[3],t[4],t[5],g[0],g[1],g[2],g[3],g[4],f[2],f[3],f[4]);
	_ff_sum_6_mults_d2(o[3],f[0],f[1],t[1],t[2],t[3],t[4],g[0],g[1],g[2],g[3],f[2],f[3]);
	_ff_sum_5_mults_d1(o[2],f[0],f[1],t[1],t[2],t[3],g[0],g[1],g[2],f[1],f[2]);
	_ff_sum_3_mults_d1(o[1],f[0],t[1],t[2],g[0],g[1],f[1]);
	_ff_sum_2_mults(o[0],f[0],t[1],g[0],f[0]);
}

// computes f^2 mod g where g=x^n - g[n-2]x^{n-2} - ... - g[0] (note signs!!)  and f has degree n-1 (n must be > 0)
void ff_poly_square_mod_small (ff_t o[], ff_t f[], ff_t g[], int n)
{
	switch (n) {
	case 3: ff_poly_square_mod_3 (o,f,g); break;
	case 4: ff_poly_square_mod_4 (o,f,g); break;
	case 5: ff_poly_square_mod_5 (o,f,g); break;
	case 6: ff_poly_square_mod_6 (o,f,g); break;
	case 7: ff_poly_square_mod_7 (o,f,g); break;
	case 8: ff_poly_square_mod_8 (o,f,g); break;
	case 9: ff_poly_square_mod_9 (o,f,g); break;
	case 10: ff_poly_square_mod_10 (o,f,g); break;
	case 11: ff_poly_square_mod_11 (o,f,g); break;
	case 12: ff_poly_square_mod_12 (o,f,g); break;
	case 13: ff_poly_square_mod_13 (o,f,g); break;
	case 14: ff_poly_square_mod_14 (o,f,g); break;
	case 15: ff_poly_square_mod_15 (o,f,g); break;
	case 16: ff_poly_square_mod_16 (o,f,g); break;
	case 17: ff_poly_square_mod_17 (o,f,g); break;
	case 18: ff_poly_square_mod_18 (o,f,g); break;
	case 19: ff_poly_square_mod_19 (o,f,g); break;
	case 20: ff_poly_square_mod_20 (o,f,g); break;
	default:	ff_poly_square_mod_small_big (o,f,g,n);
	}
}

// computes the gcd of f and g where f=x^3+ax+b 
void ff_poly_gcd_x3axb (ff_t h[4], int *d_h, ff_t f[4], ff_t g[4], int d_g)
{
	register ff_t t1,t2;
	ff_t r[4], s[4];
	int d_r, d_s;
	
	if ( d_g == 3 ) {
		_ff_neg(r[2],g[2]);
		_ff_mult(t1,f[1],g[3]);
		_ff_sub(r[1],t1,g[1]);
		_ff_mult(t1,f[0],g[3]);
		_ff_sub(f[0],t1,g[0]);
		d_r = ff_poly_degree(r,2);
	} else {
		ff_poly_copy (r,&d_r,g,d_g);
	}
	if ( d_r < 0 ) { ff_poly_copy (h, d_h, f, 3);  return; }
	if ( ! d_r ) { _ff_set_one(h[0]); *d_h = 0;  return; }
	_ff_neg(s[2],r[d_r-1]);
	_ff_mult(t1,f[1],r[d_r]);
	if ( d_r > 1 ) { _ff_sub(s[1],t1,r[d_r-2]); } else { _ff_set(s[1],t1); }
	_ff_mult(s[0],f[0],r[d_r]);
	d_s = ff_poly_degree(s,2);
	if ( d_s == 2 && d_r == 2) {
		_ff_mult(t1,s[1],r[2]);
		_ff_mult(t2,r[1],s[2]);
		_ff_sub(s[1],t1,t2);
		_ff_mult(t1,s[0],r[2]);
		_ff_mult(t2,r[0],s[2]);
		_ff_sub(s[0],t1,t2);
		d_s = ff_poly_degree(s,1);
	}
	if ( d_s < 0 ) { ff_poly_copy(h,d_h,r,d_r);  return; }
	if ( ! d_s ) { _ff_set_one(h[0]); *d_h = 0;  return; }
	if ( d_r == 2 && d_s == 1 ) {
		_ff_mult(t1,r[1],s[1]);
		_ff_mult(t2,s[0],r[2]);
		_ff_sub(r[1],t1,t2);
		ff_mult(r[0],r[0],s[1]);
		d_r = ff_poly_degree(r,1);
		if ( d_r < 0 ) { ff_poly_copy (h, d_h, s, d_s);  return; }
		if ( ! d_r ) { _ff_set_one(h[0]); *d_h = 0;  return; }
	}
	if ( d_s == 2 && d_r == 1 ) {
		_ff_mult(t1,s[1],r[1]);
		_ff_mult(t2,r[0],s[2]);
		_ff_sub(s[1],t1,t2);
		ff_mult(s[0],s[0],r[1]);
		d_s = ff_poly_degree(s,1);
		if ( d_s < 0 ) { ff_poly_copy (h, d_h, r, d_r);  return; }
		if ( ! d_s ) { _ff_set_one(h[0]); *d_h = 0;  return; }
	}
	// we must have d_r==d_s==1 at this point
	_ff_mult(t1,r[0],s[1]);
	_ff_mult(t2,s[0],r[1]);
	_ff_sub(s[0],t2,t1);
	if ( _ff_zero(s[0]) ) { ff_poly_copy(h,d_h,r,d_r);  return; }
	_ff_set_one(h[0]); *d_h = 0;
	// worst case 13M+7A and some copies
	return;	
}


int ff_poly_dbl_root_d3 (ff_t r[1], ff_t f[3])
{
	ff_t g[2];
	register ff_t t0,t1;
	
	ff_poly_dbl_root_linfac_d3(g,f);
	if ( _ff_zero(g[1]) ) return 0;
	_ff_invert(t0,g[1]);
	_ff_mult(t1,t0,g[0]);
	_ff_neg(r[0],t1);
	return 1;
}

// finds the roots of a monic cubic over F_3 that is not necessarily depressed (this is special case not handled elsewhere)
int _ff_poly_roots_d3_mod3 (ff_t r[3], ff_t f[4])
{
	ff_t y, g[3];

	assert ( _ff_p == 3 );
	if ( _ff_zero(f[0]) ) {
		_ff_set_zero(r[0]);
		return 1 + ff_poly_roots_d2 (r+1, f+1, 2);
	} else {
		_ff_set_one(r[0]); ff_poly_eval (&y, f, 3, r);
		if ( _ff_zero(y) ) {
			ff_poly_remove_root_d3 (g, f, r);
			return 1 + ff_poly_roots_d2 (r+1, g, 2);
		} else {
			_ff_inc(r[0]); ff_poly_eval (&y, f, 3, r);
			if ( _ff_zero (y) ) {
				ff_poly_remove_root_d3 (g, f, r);
				return 1 + ff_poly_roots_d2 (r+1, g, 2);
			}
		}
	}
	return 0;
}

// finds the roots of a monic quartic over F_3 that is not necessarily depressed (this is special case not handled elsewhere)
int _ff_poly_roots_d4_mod3 (ff_t r[3], ff_t f[4])
{
	ff_t y, g[4];

	assert ( _ff_p == 3 );
	if ( _ff_zero(f[0]) ) {
		_ff_set_zero(r[0]);
		return 1 + _ff_poly_roots_d3_mod3 (r+1, f+1);
	} else {
		_ff_set_one(r[0]); ff_poly_eval (&y, f, 4, r);
		if ( _ff_zero(y) ) {
			ff_poly_remove_root_d4 (g, f, r);
			return 1 + _ff_poly_roots_d3_mod3 (r+1, g);
		} else {
			_ff_inc(r[0]); ff_poly_eval (&y, f, 4, r);
			if ( _ff_zero (y) ) {
				ff_poly_remove_root_d4 (g, f, r);
				return 1 + _ff_poly_roots_d3_mod3 (r+1, g);
			}
		}
	}
	return 0;
}


/*
	Finds the roots of a monic depreseed cubic x^3+ax+b over F_p.  Does not verify that poly is monic and depressed (in fact does not look at f[3] or f[2] -- this is important!).
	If roots is NULL only the # of roots is determined, which can be significantly quicker (e.g. when D is not a QR)
	The value pd is optional.  If non-null it must point to the sqrt(-D/3) in F_p (this is useful when computing 3-torsion),
	If flags is set to FF_POLY_EXACTLY_ONE_ROOT, the total # of roots is always returned, but r will be updated only when the return value is 1.

	Assumes _ff_p fits in an unsigned long
*/
int _ff_poly_roots_d3 (ff_t roots[3], ff_t f[4], ff_t *pd, int flags)
{
	ff_t g[3], h[2], w[2], v[2], vk[2], d, t, r, s;
	register ff_t a, t0, t1, t2;
	int i, k, sts;
	
	// handle p=3 separately
	if ( _ff_p == 3 ) {
		if ( _ff_zero(f[0]) ) {
			if ( _ff_zero(f[1]) ) { if ( roots && ! (flags&FF_POLY_EXACTLY_ONE_ROOT) ) { _ff_set_zero(roots[0]); _ff_set_zero(roots[1]); _ff_set_zero(roots[2]); } return 3; }
			if ( _ff_one(f[1]) ) { if ( roots ) { _ff_set_zero(roots[0]); return 1; } }
			if ( roots && ! (flags&FF_POLY_EXACTLY_ONE_ROOT) ) { _ff_set_zero(roots[0]); _ff_set_one(roots[1]); _ff_set(roots[2],_ff_negone); }
			return 3;
		}
		if ( _ff_zero(f[1]) ) { if ( roots && ! (flags&FF_POLY_EXACTLY_ONE_ROOT) ) { _ff_neg(roots[0],f[0]); _ff_neg(roots[1],f[0]); _ff_neg(roots[2],f[0]); } return 3; }
		if ( _ff_one(f[1]) ) { if ( roots ) _ff_set(roots[0],f[0]);  return 1; }
		return 0;
	}
	
	// f=x^3+ax
	if ( _ff_zero(f[0]) ) {
		_ff_neg(t,f[1]);  
		if ( ! roots ) { return (ff_residue(t)?3:1); }
		_ff_set_zero(roots[0]);
		if ( ! ff_sqrt(&d,&t) ) return 1;
		if (! (flags&FF_POLY_EXACTLY_ONE_ROOT) ) { _ff_set(roots[1],d);  _ff_neg(roots[2],d); }
		return 3;
	}
	
	// f = x^3+b
	if ( _ff_zero(f[1]) ) {
		_ff_neg(t,f[0]);
		if ( _ff_p1mod3 ) {
			if ( ! ff_cbrt(&d,&t) ) return 0;
			if ( roots && ! (flags&FF_POLY_EXACTLY_ONE_ROOT) ) { _ff_set(roots[0],d);  _ff_mult(roots[1], roots[0], _ff_cbrt_unity);  _ff_mult(roots[2], roots[1], _ff_cbrt_unity); }
			return 3;
		} else {
			if ( roots ) if ( ! ff_cbrt(roots,&t) ) { printf ("Impossible, ff_cbrt failed on input %ld when p=%ld is 2 mod 3\n", _ff_get_ui(t), _ff_p); abort(); }
			return 1;
		}
	}

	if ( ! pd ) {
		_ff_square(t0,f[0]);  _ff_set_ui(t1,27);  _ff_mult(t2,t0,t1);				// t2 = 27f1^2
		_ff_square(t0,f[1]);  _ff_mult(t1,t0,f[1]);  _ff_x2(t1); _ff_x2(t1);		// t1 = 4f0^3
		_ff_add(t0,t1,t2);  _ff_mult(t,t0,_ff_third);						// t = -D/3 = (27f1^2+4f0^3)/3
		if ( _ff_zero(t) ) {											// t=0 => D=0 => curve is singular
			if ( roots && ! (flags&FF_POLY_EXACTLY_ONE_ROOT) ) {			// distinct root is 3b/a (we know a!=0), and (-3b/(2a) is a double root).
				_ff_invert(t0,f[1]);									// 1/a
				_ff_add(t1,f[0],f[0]);  _ff_addto(t1,f[0]);					// 3b
				_ff_mult(roots[0],t0,t1);								// 3b/a
				_ff_neg(t0,roots[0]);
				_ff_mult(roots[1],t0,_ff_half);							// -3b/(2a)
				_ff_set(roots[2],roots[1]);
			}
			return 3;											
		}
		sts = ff_sqrt_ext(&d,&t);										// use extended sqrt so we get a solution in F_p^2 if necessary
	} else {
		sts = 1;
		_ff_set(d,*pd);
	}
	_ff_mult(t2,_ff_half,_ff_third);									// t2=1/6 (used later)
	_ff_mult(t0, _ff_half, f[0]);
	_ff_neg(a,t0);												// a = -f[0]/2 (used later)
	if ( _ff_p1mod3 ) {
		if ( ! sts ) {											// p=1mod3 => -1/3 is a QR => (t not a QR => D not a QR => f has an even # of factors (by Stickleberger), so 1 root
			if ( ! roots ) return 1;								// if all the caller wants is the number of roots, we are done
			_ff_square(t0,t2);									// 1/36
			_ff_mult(s,t,t0);									// s = -D/108 is a not a QR, we will now work in the field F_p^2 = F_p[z]/(z^2-s) to compute (z+a)^1/3
			// compute (z+t1)^((p+1)*m) in F_p^2, this will be in the 3-Sylow subgroup of F_p^2, necessarily in F_p (here p=3^e*m+1, m not divisible by 3)
			// in the process, we compute g=(z+1)^((m-i)/3) where i is congruent to m mod 3 and h=(z+a)^m
			if ( _ff_p3_m1mod3 ) {
				ff_poly_xpan_mod_d2(g,a,(_ff_p3_m-1)/3,&s);			// g = (z+a)^((m-1)/3)
				ff2_square_s(h,g,s);  ff2_mult_s (h,h,g,s);			// h = g^3 = (z+a)^(m-1)
				ff2_mult_zpa_s(h,h,a,s);							// h = g^3*(z+a) = (z+a)^m
			} else {
				ff_poly_xpan_mod_d2(g,a,(_ff_p3_m-2)/3,&s);			// g = (z+a)^((m-2)/3)
				ff2_square_s (h, g, s);  ff2_mult_s (h, h, g, s);			// h = g^3 = (z+a)^(m-2)
				_ff_square(t0,a); _ff_add(w[1],a,a);  _ff_add(w[0],t0,s);	// w = (z+a)^2
				ff2_mult_s (h,h,w,s);							// h = g^3*(z+a)^2 = (z+a)^m
			}
			ff2_norm_s (&t,h,s);									// norm(h)=((z+a)^(m*(p+1)) is in the 3-Sylow of F_p
			if ( ! ff_3Sylow_invcbrt(&r,&t) ) {						// (z+a) is not a cubic residue - this should be impossible
				printf ("(z+%ld)^%ld = ", _ff_get_ui(a), _ff_p3_m); ff_poly_print(h,1);
				printf ("norm(h)=%ld is not a cube in F_%ld\n", _ff_get_ui(t), _ff_p);
				printf ("z+%ld is not a CR in F_p[z]/(z^2-%ld)\n", _ff_get_ui(a), _ff_get_ui(s));
				ff_poly_print(f,3);
				abort();
			}
			// we now know (z+a)^-((p+1)m)/3
			ff2_norm_s(&t,g,s);									// t = norm(g) = g^(p+1)
			ff2_exp_ui_s(h,h,(_ff_p+1)/(3*_ff_p3_m),s);				// compute h^(3^(e-1))
			if ( !_ff_p3_m1mod3 ) {
				// We need to construct (z+a)^n where n = (2(p+1)m+1)/3 = 2(p+1)(m-2)/3 + 2(2p+1)/3 + 1 = 2(p+1)(m-2)/3 + 2(2*3^(e-1)*m+1) + 1
				// We then have (z+a)^n = norm(g)^2*((h^(3^(e-1))^2*(z+a))^2*(z+a)   (note that we also need to square r to get (z+a)^(-2(p+1)m/3) here)
				ff_square(t,t);  ff_square(r,r);
				ff2_square_s(h,h,s);  ff2_mult_zpa_s (h,h,a,s);  ff2_square_s(h,h,s);
			}// when m is 1 mod 3:
				// We need to construct (z+a)^n where n = ((p+1)m+1)/3 = (p+1)(m-1)/3 + (p-1)/3 + 1 = (p+1)(m-1)/3 + 3^(e-1)*m + 1
				// We then have (z+a)^n = norm(g)*h^(3^(e-1))*(z+a) 
			ff2_scalar_mult(g,t,h);								
			ff2_mult_zpa_s(g,g,a,s);								// g = (z+a)^n
			ff2_scalar_mult(g,r,g);								// g = (z+a)^(1/3)
			// We know that (g-f1/(3g)) is a root of f.  To get others we multiply by cube roots of unity (which are conveniently in F_p).
			// We use this to find the root in F_p (there must be exactly one, since we know f has 2 factors by Stickleberger)
			// The multiple h of g yielding a root in F_p must have the property that -f1/3=norm(h), which we can test without inverting, and we then have tr(h) as our root
			_ff_mult(t0,f[1],_ff_third);
			ff_negate(t0);
			for ( i = 0 ; i < 3 ; i++ ) {  ff2_norm_s(&t,g,s);  if ( _ff_equal(t,t0) ) break;  ff2_scalar_mult(g,_ff_cbrt_unity,g); }
			if ( i == 3 ) {
				printf ("g=%ldz+%ld is a cbrt of (z+%ld) in F_p[z]/(z^2-%ld), but couldn't find an F_%ld root of 2-factor f = ", _ff_get_ui(g[1]), _ff_get_ui(g[0]),_ff_get_ui(a), _ff_get_ui(s), _ff_p);
				ff_poly_print(f,3);
				printf ("t0=%ld, cube root of unity is %ld\n", _ff_get_ui(t0), _ff_get_ui(_ff_cbrt_unity));
				abort();
			}
			ff2_trace(roots,g);
			return 1;
		} else {
			ff_mult(d,d,t2);									// d = sqrt(t)/6 = sqrt(f1^3/27+f0^2/4) (this is the square root in Cardano's formula for depressed cubics)
			_ff_add(t,a,d);										// t = -f0/2 + d (first cube to check)
			if ( ! ff_cbrt_invcbrt(&r,&s,&t) ) return 0;
			if ( ! roots || (flags&FF_POLY_EXACTLY_ONE_ROOT) ) return 3;	// we know the roots exist but we don't need to find them
			_ff_mult(t0,s,_ff_third); _ff_mult(t1,t0,f[1]); _ff_neg(s,t1);
			_ff_add(roots[0],r,s);
			_ff_square(t0,_ff_cbrt_unity);							// could cache this in ff.c
			_ff_mult(t1,r,_ff_cbrt_unity);
			_ff_mult(t2,s,t0);
			_ff_add(roots[1],t1,t2);								// second root is wr+w^2s where w is a primitive cbrt of unity
			_ff_mult(t1,r,t0);
			_ff_mult(t2,s,_ff_cbrt_unity);
			_ff_add(roots[2],t1,t2);								// third root is wr^2+ws
			return 3;
		}
	} else {													// p=2mod3
		if ( sts ) {												// p=2mod3 => -1/3 not a QR => (-D/3 a QR => D not a QR => f has an even # of factors (by Stickleberger)) => f has one root
			if ( ! roots ) return 1;								// if we only need the # of factors, we are done.
			ff_mult(d,d,t2);									// d = sqrt(t)/6 = sqrt(f1^3/27+f0^2/4) (this is the square root in Cardano's formula for depressed cubics)
			_ff_add(t,a,d);										// t = -f0/2 + d (first cube)
			ff_cbrt_invcbrt (&r,&s, &t);							// p is 2 mod 3 so we know there is exactly 1 cube root
//_ff_square(t0,r); _ff_mult(t1,t0,r); if ( ! _ff_equal(t1,t) ) { printf ("cube root failed p=%ld, t=%ld\n", _ff_p, _ff_get_ui(t)); abort(); }
//_ff_mult(t0,r,s); if ( ! _ff_one(t0) ) { printf ("invcube root failed p=%ld, t=%ld\n", _ff_p, _ff_get_ui(t)); abort(); }
//			_ff_sub(t,a,d);										// t = -f0/2 - d (second cube)
//			ff_cbrt(&s,&t);
			// for t = -f0/2 - d we have s = -f1/(3r) as a cube-root of t (see Lenstra Key XTR improvements, Alg. 3.1)
			_ff_mult(t0,s,_ff_third); _ff_mult(t1,t0,f[1]); _ff_neg(s,t1);
			_ff_add(t0,r,s);										// there was no choice of cuberoots (and changing sign of d doesn't change r+s=t), so this has got to be the root
			_ff_set(roots[0],t0);
//_ff_square(t0,roots[0]); _ff_addto(t0,f[1]); _ff_multadd(t1,t0,roots[0],f[0]);
//if ( ! _ff_zero(t1) ) { printf ("failed to find root of "); ff_poly_print(f,3); printf ("p=%ld, d=%ld, t=%ld, r=%ld, s=%ld, root=%ld\n", _ff_p, _ff_get_ui(t), _ff_get_ui(r), _ff_get_ui(s), _ff_get_ui(roots[0])); abort(); }
			return 1;
		} else {												// in this case t = -D/3 is not a QR, D=d^2 is a QR, and we know that f has 1 or 3 factors
			// for p=2mod3, handle roots=NULL separately, since we can be much more efficient
			if ( ! roots || (flags&FF_POLY_EXACTLY_ONE_ROOT) ) {
				_ff_square(t0,t2);								// 1/36
				_ff_mult(s,t,t0);								// s = -D/108 is a not a QR, we will now work in the field F_p^2 = F_p[z]/(z^2-s)
				ff_poly_xpan_mod_d2(g,a,(_ff_p+1)/3,&s);			// (z+a)^((p+1)/3) is in F_p iff (z+a)^((p+1)(p-1)/3)=1 iff (z+a) is a cubic residue
				return ( _ff_zero(g[1]) ? 3 : 0 );
			}

			// in this case we will be slightly less efficient and work in the standard basis for F_p^2=F_p[z](z^2-s) because we
			// want don't want to have to compute a new generator of the Sylow 3-subgroup in F_p^2 every time 
			_ff_set(v[0],a);  _ff_mult(v[1],d,t2);						// v is the element of F_p^2 (in the standard basis) whose cuberoot we need
			// this code is slightly inefficient, but we are trying to compute h=v^((p+1)/3) as quickly as possible (if h is not in F_p, f is irreducible and we are done)
			// while also saving values we will need to compute a cube root of v when f splits
			if ( _ff_p3_m1mod3 ){ k=1; ff2_set(vk,v); }
			else { k=2; ff2_square(vk,v); }						// vk = v^k
			ff2_exp_ui(g,v,(_ff_p3_m-k)/3);							// g = v^((m-k)/3)
			ff2_square(h,g);  ff2_mult(h,h,g);						// h = v^(m-k)
			ff2_mult(h,h,vk);	
			if ( _ff_p3_e > 1 ) {
				ff2_exp_ui(w,h,(_ff_p+1)/(3*_ff_p3_m)-1);			// w = v^((3^(e-1)-1)m) = v^(((p+1)-3m)/3)
				_ff_mult(t1,w[1],h[0]);  _ff_mult(t2,w[0],h[1]);		// compute z coefficient of w*h
				_ff_add(t0,t1,t2);
				if ( ! _ff_zero(t0) ) return 0;						// if w*h=v^((p+1)/3) is not in F_p, then v^((p^2-1)/3)=h^(p-1) is not 1 and v is not a cube, and o.w. it is
			} else {
				if ( ! _ff_zero(h[1]) ) return 0;						// here h=v^m=v^((p+1)/3) since e=1
			}
			// We now know v is a cube and there are three roots
			if ( ! roots ) return 3;
			ff2_norm(&t,g);									// t = N(g)=g^(p+1)=v^((p+1)(m-k)/3) where k=m mod 3
			ff2_scalar_mult(g,t,g);								// g = v^((p+2)(m-k)/3)
			// When e==1, life is pretty simple: we know v is a cube, and there is only one cube in the Sylow 3-subgroup, name 1, and its inverse cube root is 1
			// It follows that we can assume v^((p-1)m/3)=v^((p-2)m/3)=1.
			if ( _ff_p3_e == 1 ) {
				if ( k==2 ) {
					ff2_mult(g,g,h);							// g = v^((p+2)(m-2)/3+m)=v^((pm+2m-2p-4+3m)/3)=v^((pm+2m-2(3m-1)-4+3m)/3)=v^(((p-1)m-2)/3)
				} else {										// note that e=1 and k=1 implies g=v^((p+1)(m-1)/3) = v^(((p-1)m-1)/3)
					ff2_square(g,g);							// g=v^((2(p-1)m-2)/3)
				}
				ff2_mult(g,g,v);								// g = v^(((3-k)(p-1)m+1)/3) = v^(1/3}
				ff2_trace(roots,g);								// Every choice of cube root of v corresponds to a root of f, so we can just take the trace
				ff2_setup_cbrt();								// we still need to get the cube root of unity (this is annoying)
				ff2_mult(g,g,_ff2_cbrt_unity);
				ff2_trace(roots+1,g);
				ff2_mult(g,g,_ff2_cbrt_unity);
				ff2_trace(roots+2,g);
				return 3;
			}
			// now we have e>1
			if ( k==2 ) { ff2_square(w,w);  ff2_mult(w,w,h); }			// w = v^((k*3^(e-1)-1)m)
			ff2_mult(g,g,w);									// g = v^(((p-1)m-k)/3)  (note that (p+2)(m-k)/3 + (k*3^(e-1)-1)m = (pm+2m-kp-2k+k3^em-3m)/3
															//                                                                                                     = (pm-m-kp-2k+k(p+1))/3=((p-1)m-k)/3
			ff2_square(h,g);  ff2_mult(h,h,g);						// h = v^((p-1)m-k)
			ff2_mult(h,h,vk);									// h = v^((p-1)m) is in the 3-Sylow subgroup, and it must be a cubic residue
			if ( ! ff2_3Sylow_invcbrt(h,h) ) {						// h = v^(-(p-1)m)/3)  if e=1 
				printf ("Impossible situation in _ff_poly_roots_d3 with p=%ld,m=%ld:  v^((p-1)m)=(%ldz+%ld) is not a cube in 3-Sylow of F_p^2=F_p[z]/(z^2-%ld) ",
					   _ff_p, _ff_p3_m, _ff_get_ui(h[1]), _ff_get_ui(h[0]), _ff_get_ui(s));
				ff_poly_print(f,3);
				abort();
			}
			if ( k==1 ) {
				ff2_square(h,h);								// h = v^(-2(p-1)m/3)
				ff2_square(g,g);								// g = v^(2((p-1)m-1)/3)
				ff2_mult(g,g,v);								// g = v^(2((p-1)m-1)/3+1) = v^((2(p-1)m+1)/3
			} else {
				ff2_mult(g,g,v);								// g = v^((p-1)m-2)/3+1 = v^((p-1)m+1/3
			}
			ff2_mult(g,g,h);									// g = v^(1/3)
			ff2_trace(roots,g);									// Every choice of cube root of v corresponds to a root of f, so we can just take the trace of each
			ff2_mult(g,g,_ff2_cbrt_unity);
			ff2_trace(roots+1,g);
			ff2_mult(g,g,_ff2_cbrt_unity);
			ff2_trace(roots+2,g);
			return 3;
		}
	}
	puts ("Unreachable code in _ff_poly_roots_d3!");
	abort();
}

// Compute roots of poly of degree 2, or less.  We don't assume leading coefficient is non-zero (or that f is monic))
int ff_poly_roots_d2 (ff_t r[2], ff_t f[3], int d_f)
{
	ff_t t1,t2,D;
	
	switch (d_f) {
	case 0:  if ( ! _ff_zero(f[0]) ) return 0;  _ff_set_zero(r[0]); return 1;
	case 1: 
		if ( _ff_one(f[1]) ) { _ff_set_one(t1); } else { if ( _ff_zero(f[1]) ) return ff_poly_roots_d2(r,f,0); _ff_invert(t1,f[1]); }
		ff_negate(t1); _ff_mult(r[0],t1,f[0]);  return 1;
	case 2:
		_ff_square(t1,f[1]);
		_ff_mult(t2,f[2],f[0]);
		_ff_x2(t2);  _ff_x2(t2);
		_ff_subfrom(t1,t2);
		if ( ! ff_sqrt(&D,&t1) ) return 0;
		if ( _ff_one(f[2]) ) {					
			_ff_set(t1,_ff_half);
		} else {
			if ( _ff_zero(f[2]) ) return ff_poly_roots_d2(r,f,1);		// put zero-test here to optimize for the monic case
			_ff_add(t2,f[2],f[2]);
			_ff_invert(t1,t2);
		}
		_ff_sub(t2,D,f[1]);
		_ff_mult(r[0],t1,t2);
		ff_negate(D);
		_ff_sub(t2,D,f[1]);
		_ff_mult(r[1],t1,t2);
		return 2;
	default:
		printf ("Invalid degree %d in ff_poly_roots_d2\n", d_f); abort();
	}
}

/*
	Computes the roots of f(x)=x^4+ax^2+bx+c, solving by radicals.
	Returns all roots, with multiplicity

	Currently r is required (and must have space for 4 roots, even if it is known that there aren't 4 roots)
	We could easily extend to make it optional and only return a count.
	The parameter pd is optional.  If specified it must point to sqrt(-D/3) (this accelerates the factorization of the cubic resolvent).
	Currently the flags parameter is unused
*/
int _ff_poly_roots_d4 (ff_t r[4], ff_t f[5], ff_t *pd, int flags)
{
	ff_t g[4], h[3], s[3], x, y, t;
	register ff_t t0,t1,t2,w1,w2;
	register int i,k,res;

	if ( _ff_p == 3 ) return _ff_poly_roots_d4_mod3 (r, f);
	
	if ( _ff_zero(f[1]) ) {
		// f = x^4+f2x^2+f0 = g(x^2) where g is monic quadratic
		_ff_set_one(g[2]);  _ff_set(g[1],f[2]);  _ff_set(g[0],f[0]);
		k = ff_poly_roots_d2(s,g,2);
		if ( ! k ) return 0;
		if ( ff_sqrt(r,s) ) {
			_ff_neg(r[1],r[0]);
			// could optimize for double root here, but why bother
			if ( ff_sqrt(r+2,s+1) ) { _ff_neg(r[3],r[2]); return 4; }
			return 2;
		}
		if ( ff_sqrt(r,s+1) ) { _ff_neg(r[1],r[0]); return 2; }
		return 0;
	}
	ff_depressed_cubic_resolvent(&t,g,f);
	k = _ff_poly_roots_d3(r,g,pd,0);				// put the roots of g into r for now, they will get replaced by roots of f below.
	if ( k ) {
		res = 0;
		for ( i = 0 ; i < k ; i++ ) {
//			ff_poly_eval (&y, g, 3, r+i);
//			if ( ! _ff_zero(y) ) { printf ("(%d roots) %ld is not a root of ", k, _ff_get_ui(r[i]));  ff_poly_print(g,3); abort(); }
			_ff_sub(x,t,r[i]);		// translate and negate
			if ( ! ff_sqrt(s+i,&x) ) break;
			if ( ++res == 2 ) break;
		}
		if ( res==2 ) {
			_ff_mult(t1,s[0],s[1]);
			ff_negate(t1);
			_ff_invert(t0,t1);
			_ff_mult(s[2],t0,f[1]);
			_ff_square(t1,s[2]);
			_ff_sub(x,r[2],t);
			ff_negate(x);
			if ( ! _ff_equal(t1,x) ) { printf("%ld: division failed s[0]=%ld, s[1]=%ld, d=%ld\n", _ff_p, _ff_get_ui(s[0]), _ff_get_ui(s[1]), _ff_get_ui(t0));  ff_poly_print(f,4); ff_poly_print(g,3); abort(); }
			_ff_add(t1,s[0],s[1]);  _ff_addto(t1,s[2]); ff_mult(t1,t1,_ff_half);
			_ff_set(r[0],t1);
			ff_negate(s[0]); ff_negate(s[2]);
			_ff_add(t1,s[0],s[1]);  _ff_addto(t1,s[2]); ff_mult(t1,t1,_ff_half);
			_ff_set(r[1],t1);
			ff_negate(s[0]); ff_negate(s[1]);
			_ff_add(t1,s[0],s[1]);  _ff_addto(t1,s[2]); ff_mult(t1,t1,_ff_half);
			_ff_set(r[2],t1);
			ff_negate(s[0]); ff_negate(s[2]);
			_ff_add(t1,s[0],s[1]);  _ff_addto(t1,s[2]); ff_mult(t1,t1,_ff_half);
			_ff_set(r[3],t1);
			return 4;
		} else if ( k==1 && res==1 ) {
			_ff_set_one(h[2]);
			_ff_set(h[1],s[0]);
			_ff_square(t1,s[0]);						// t1 = s^2
			_ff_add(w1,t1,f[2]);
			ff_mult(w1,w1,_ff_half);					// w1 = (s^2+f2)/2
			_ff_square(t0,f[2]);
			_ff_add(t2,f[0], f[0]); _ff_x2(t2);			// t2 = 4f0
			_ff_subfrom(t0,t2);						// t0 = f2^2-4f0
			_ff_add(t2,f[2],f[2]);						// t2 = 2f2
			_ff_add(w2,t1,t2);						// w2 = s^2+2f2
			ff_mult(w2,w2,t1);						// w2 = s^4+2f2s^2
			_ff_addto(w2,t0);
			ff_mult(w2,w2,s[0]);					// w2 = s(s^4+2f2s^2+f2^2-4f0)				
			_ff_invert(t0,f[1]);
			ff_mult(t1,_ff_half,t0);
			ff_mult(w2,w2,t1);
			_ff_sub(h[0],w1,w2);
			if ( ff_poly_roots_d2(r,h,2) != 2 ) {
				ff_negate(h[1]);
				_ff_add(h[0],w1,w2);
				if ( ff_poly_roots_d2(r,h,2) != 2 ) { printf("%ld: fail 2-1-1 split with s=%ld\n", _ff_p, _ff_get_ui(s[0]));  ff_poly_print(f,4);  ff_poly_print(g,3);  ff_poly_print(h,2); abort(); }
			}
			return 2;
		} else {
			// check for a double root in the resolvent, we potentially recompute a square-root here, but this is a rare case so we won't worry about it
			if ( k==3 ) {
				if ( _ff_equal(r[0],r[1]) ) {
					_ff_sub(x,t,r[2]);		// translate and negate
				} else if ( _ff_equal(r[0],r[2]) ) {
					_ff_sub(x,t,r[1]);		// translate and negate
				} else if ( _ff_equal(r[1],r[2]) ) {
					_ff_sub(x,t,r[0]);		// translate and negate
				} else {
					return 0;
				}
				if ( ! ff_sqrt(s,&x) ) return 0;
				_ff_mult(r[0],s[0],_ff_half);
				// we don't know the sign, so check the root and negate if needed
				ff_poly_eval(s,f,4,r);
				if ( ! _ff_zero(s[0]) ) ff_negate(r[0]);
				_ff_set(r[1],r[0]);
				return 2;
			}
			return 0;
		}
	} else {
		// this is the annoying case, we know we have a 1,3 split, but getting the value of the root is expensive
		// we need sqrt(-z) in F_p[z]/h(z) where h(z) is the undepressed resolvent
		// so we want sqrt(t-z) in F_p[z]/g(z) which is sqrt(z+t) in F_p[z]/k(z) where k(z)=-g(-z)=z^3+g1z-g0z.  In terms of z^3-rz-s, r=-g1 and s=g0
		_ff_neg(t1,g[1]);
		if ( ! ff3_trsqrt_zpa_mod_rs(&x,t,t1,g[0]) ) { printf ("%ld: impossible failure of ff3_trsqrt_zpa_mod_rs in ff_poly_roots_d4!", _ff_p); ff_poly_print(f,4); ff_poly_print(g,3); abort(); }
		ff_mult(x,x,_ff_half);
		ff_poly_eval(&y,f,4,&x);
		if ( ! _ff_zero(y) ) ff_negate(x);
		_ff_set(r[0],x);
		return 1;
	}
	puts ("Unreachable code in ff_poly_roots_d4!");
	abort();
}

// computes g(x) = h1(h2(x)) mod f(x) with overlap permittted, assumes f=x^4 - f2x^2 - f1x - f0
void ff_poly_compose_mod_4 (ff_t g[4], ff_t h1[4], ff_t h2[4], ff_t f[3])
{
	ff_t w[4],t0;
	
	_ff_mult(w[3],h1[3],h2[3]); _ff_mult(w[2],h1[3],h2[2]); _ff_mult(w[1],h1[3],h2[1]); _ff_mult(w[0],h1[3],h2[0]); _ff_addto(w[0],h1[2]);			// w = h1[3]*h2 + h1[2]
	ff_poly_mult_mod_4 (w,w,h2,f);																						// w = h2*(h1[3]*h2+h1[2])
	_ff_addto(w[0],h1[1]);  _ff_set(t0,h1[0]); ff_poly_mult_mod_4 (g,w,h2,f); _ff_addto(g[0],t0);											// w = h2*(h2*(h1[3]*h2+h1[2])+h1[1])
}

// computes g(x) = h1(h2(x)) mod f(x) with overlap permittted, assumes f = x^6 - f4x^4 - f3x^3 - ... - f0
void ff_poly_compose_mod_6 (ff_t g[6], ff_t h1[6], ff_t h2[6], ff_t f[5])
{
	ff_t w[6],t0;
	
	_ff_mult(w[5],h1[5],h2[5]); _ff_mult(w[4],h1[5],h2[4]); _ff_mult(w[3],h1[5],h2[3]);
	_ff_mult(w[2],h1[5],h2[2]); _ff_mult(w[1],h1[5],h2[1]); _ff_mult(w[0],h1[5],h2[0]); _ff_addto(w[0],h1[4]);								// w = h1[5]*h2 + h1[4]
	ff_poly_mult_mod_6 (w,w,h2,f);  _ff_addto(w[0],h1[3]);																		// w = h2*(h1[5]*h2+h1[4]) + h1[3]
	ff_poly_mult_mod_6 (w,w,h2,f);  _ff_addto(w[0],h1[2]);																		// w = h2*h2*(h1[5]*h2+h1[4])+h1[3]) + h1[2]
	ff_poly_mult_mod_6 (w,w,h2,f);  _ff_addto(w[0],h1[1]);																		// w = h2*(h2*h2*(h1[5]*h2+h1[4])+h1[3])+h1[2]) + h1[1]
	_ff_set(t0,h1[0]); ff_poly_mult_mod_6 (g,w,h2,f); _ff_addto(g[0],t0);															// w = h2*(h2*(h2*h2*(h1[5]*h2+h1[4])+h1[3])+h1[2])+h1[1]) + h[0]
}


// Computes g = h^n mod f, where f is of the form x^2-f1x-f0 and g is degree 1 (possibly with zero leading coefficients),  h and g may overlap
// standard right-to-left binary exp  with a 2-bit window here
void ff_poly_exp_mod_2 (ff_t g[2], ff_t h[2], unsigned long n, ff_t f[2])
{
	ff_t x[4][2];
	register unsigned long m;
	register int i, j;

	if ( !n) { _ff_set_zero(g[1]); _ff_set_one(g[0]); return; }
	_ff_set(x[1][0],h[0]); _ff_set(x[1][1],h[1]);
	ff_poly_square_mod_2 (x[2],x[1],f);
	ff_poly_mult_mod_2 (x[3],x[2],x[1],f);
	i = _asm_highbit(n);
	if ( i&1 ) i--;
	m = 3UL<<i;
	j = (m&n)>>i;
	_ff_set(g[0],x[j][0]); _ff_set(g[1],x[j][1]);
	for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
		ff_poly_square_mod_2 (g,g,f);  ff_poly_square_mod_2 (g,g,f);
		j = (m&n)>>i;
		if ( j ) ff_poly_mult_mod_2 (g,g,x[j],f);
	}
}


// Computes (x+a)^n modulo x^2-f0 (note the sign).  Assumes n < 2^63
void ff_poly_xpan_mod_d2 (ff_t g[2], ff_t a, unsigned long n, ff_t f[1])
{
	register ff_t t1,g0,g1,f0;
	register unsigned long m;
	int i;

	if ( ! n ) { _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n == 1 ) { _ff_set_one(g[1]);  _ff_set(g[0],a); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	_ff_set(f0,f[0]);
	_ff_add(g1,a,a);  _ff_square(g0,a);  _ff_addto(g0,f0);		// g(x) = g1x+g0 = (x+a)^2 mod f(x) where f(x)=x^2-f0
	for (;;) {
		if ( (n&m) ) {									// compute (x+a)g mod f
			_ff_mult(t1,a,g1); _ff_addto(t1,g0);				// t1 = g1a+g0
			_ff_sum_2_mults(g0,g0,g1,f0,a);				// g0 = g0a+g1f0
			_ff_set(g1,t1);
		}
		m >>= 1;
		if ( ! m ) break;
		// square g mod f
		_ff_square(t1,g1);								// t1 = g1^2
		_ff_dbl_mult(g1,g0,g1);							// g1 = 2g0g1
		_ff_sum_2_mults (g0,g0,t1,f0,g0);					// g0 = g0^2+g1^2a
	}
	_ff_set(g[0],g0);  _ff_set(g[1],g1);
}

// Computes x^n modulo a cubic of the form x^3-f1x-f0 (note signs).  Assumes n < 2^63
void ff_poly_xn_mod_d3 (ff_t g[3], unsigned long n, ff_t f[2])
{
	register unsigned long m;
	register int i;

	if ( !n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return; }
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set_zero(g[0]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	
	if ( n&m ) {
		_ff_set_zero(g[2]); _ff_set(g[1],f[1]); _ff_set(g[0],f[0]);		// start with x^3
	} else {
		_ff_set_one(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);			// start with x^2;
	}
	m >>= 1;
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_3 (g,g,f);
		} else {
			ff_poly_square_mod_3 (g,g,f);
		}
		m >>= 1;
	}
}

// Computes (x+a)^n modulo a cubic of the form x^3-f1x-f0 (note signs).  Assumes n < 2^63
void ff_poly_xpan_mod_d3 (ff_t g[3], ff_t a, unsigned long n, ff_t f[2])
{
	register ff_t t1,t2;
	register unsigned long m;
	int i;
	
	if ( ! n ) { _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	_ff_set_zero(g[2]); _ff_set_one(g[1]);  _ff_set(g[0],a);
	if ( n==1 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	while (m) {
		ff_poly_square_mod_3 (g,g,f);
		if ( (n&m) ) {
			_ff_mult(t2,a,g[2]); _ff_addto(t2,g[1]);					// g2 = ag2+g1
			_ff_sum_2_mults(t1,g[1],g[2],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g2f1+g0
			_ff_sum_2_mults(g[0],g[0],g[2],f[0],a);					// g0 = ag0-g2f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2);
		}
		m >>= 1;
	}
}


// Computes g = h^n mod f, where f is of the form x^4-f2x^2-f1x-f0 and g is degree 3 (possibly with zero leading coefficients),  h and g may overlap
// standard right-to-left binary exp  with a 2-bit window here
void ff_poly_exp_mod_4 (ff_t g[4], ff_t h[4], unsigned long n, ff_t f[3])
{
	ff_t x[4][4];
	register unsigned long m;
	register int i, j;

	if ( !n) { _ff_set_zero(g[3]); _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return; }
	for ( i = 0 ; i <= 3 ; i++ ) _ff_set(x[1][i],h[i]);
	ff_poly_square_mod_4 (x[2],x[1],f);
	ff_poly_mult_mod_4 (x[3],x[2],x[1],f);
	i = _asm_highbit(n);
	if ( i&1 ) i--;
	m = 3UL<<i;
	j = (m&n)>>i;
	_ff_set(g[0],x[j][0]); _ff_set(g[1],x[j][1]); _ff_set(g[2],x[j][2]); _ff_set(g[3],x[j][3]);
	for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
		ff_poly_square_mod_4 (g,g,f);  ff_poly_square_mod_4 (g,g,f);
		j = (m&n)>>i;
		if ( j ) ff_poly_mult_mod_4 (g,g,x[j],f);
	}
}


void ff_poly_xn_mod_d4 (ff_t g[4], unsigned long n, ff_t f[3])
{
	register unsigned long m;
	register int i;

	if ( !n ) { _ff_set_zero(g[3]); _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return; }
	if ( n== 1 ) { _ff_set_zero(g[3]); _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set_zero(g[0]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	_ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n&m ) {
		_ff_set_one(g[3]); _ff_set_zero(g[2]);						// start with x^3
	} else {
		_ff_set_zero(g[3]); _ff_set_one(g[2]);						// start with x^2;
	}
	m >>= 1;
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_4 (g,g,f);
		} else {
			ff_poly_square_mod_4 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d4 (ff_t g[4], ff_t a, unsigned long n, ff_t f[3])
{
	register ff_t t1,t2,t3;
	register unsigned long m;
	register int i;
	
	_ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t3,a,g[3]); _ff_addto(t3,g[2]);					// g3 = ag3+g2
			_ff_sum_2_mults(t2,g[2],g[3],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g3f2+g1
			_ff_sum_2_mults(t1,g[1],g[3],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g3f1+g0
			_ff_sum_2_mults(g[0],g[0],g[3],f[0],a);					// g0 = ag0-g3f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3);
		}
		m >>= 1;
		if ( ! m ) return;
		ff_poly_square_mod_4 (g,g,f);
	}
}


void ff_poly_xn_mod_d5 (ff_t g[5], unsigned long n, ff_t f[4])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[4]);
	if ( !n ) {  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return; }
	if ( n== 1 ) {  _ff_set_zero(g[3]); _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set_zero(g[0]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	_ff_set_zero(g[1]); _ff_set_zero(g[0]);
	// we could save one round of squaring by handling the first 3 bits of n initially
	if ( n&m ) {
		_ff_set_one(g[3]); _ff_set_zero(g[2]);						// start with x^3
	} else {
		_ff_set_zero(g[3]); _ff_set_one(g[2]);						// start with x^2;
	}
	m >>= 1;
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_5 (g,g,f);
		} else {
			ff_poly_square_mod_5 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d5 (ff_t g[5], ff_t a, unsigned long n, ff_t f[4])
{
	register ff_t t1,t2,t3,t4;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[4]); _ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t4,a,g[4]); _ff_addto(t4,g[3]);					// g4 = ag4+g3
			_ff_sum_2_mults(t3,g[3],g[4],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g4f3+g2
			_ff_sum_2_mults(t2,g[2],g[4],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g4f2+g1
			_ff_sum_2_mults(t1,g[1],g[4],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g4f1+g0
			_ff_sum_2_mults(g[0],g[0],g[4],f[0],a);					// g0 = ag0-g4f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_5 (g,g,f);
	}
}


void ff_poly_xn_mod_d6 (ff_t g[6], unsigned long n, ff_t f[5])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[5]); _ff_set_zero(g[4]);
	if ( !n ) {  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return; }
	if ( n== 1 ) {  _ff_set_zero(g[3]); _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set_zero(g[0]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	_ff_set_zero(g[1]); _ff_set_zero(g[0]);
	// we could save one round of squaring by handling the first 3 bits of n initially
	if ( n&m ) {
		_ff_set_one(g[3]); _ff_set_zero(g[2]);						// start with x^3
	} else {
		_ff_set_zero(g[3]); _ff_set_one(g[2]);						// start with x^2;
	}
	m >>= 1;
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_6 (g,g,f);
		} else {
			ff_poly_square_mod_6 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d6 (ff_t g[6], ff_t a, unsigned long n, ff_t f[5])
{
	register ff_t t1,t2,t3,t4,t5;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[5]);_ff_set_zero(g[4]); _ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t5,a,g[5]); _ff_addto(t5,g[4]);					// g5 = ag5+g4
			_ff_sum_2_mults(t4,g[4],g[5],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g5f4+g3
			_ff_sum_2_mults(t3,g[3],g[5],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g5f3+g2
			_ff_sum_2_mults(t2,g[2],g[5],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g5f2+g1
			_ff_sum_2_mults(t1,g[1],g[5],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g5f1+g0
			_ff_sum_2_mults(g[0],g[0],g[5],f[0],a);					// g0 = ag0-g5f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_6 (g,g,f);
	}
}


void ff_poly_xn_mod_d7 (ff_t g[7], unsigned long n, ff_t f[6])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 7 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 4 + ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	if ( i == 7 ) {
		_ff_set(g[5],f[5]); _ff_set(g[4],f[4]); _ff_set(g[3],f[3]);_ff_set(g[2],f[2]);_ff_set(g[1],f[1]);_ff_set(g[0],f[0]);
	} else {
		_ff_set_one(g[i]);
	}
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_7 (g,g,f);
		} else {
			ff_poly_square_mod_7 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d7 (ff_t g[7], ff_t a, unsigned long n, ff_t f[6])
{
	register ff_t t1,t2,t3,t4,t5,t6;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[6]); _ff_set_zero(g[5]);_ff_set_zero(g[4]); _ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t6,a,g[6]); _ff_addto(t6,g[5]);					// g6 = ag6+g5
			_ff_sum_2_mults(t5,g[5],g[6],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g6f5+g4
			_ff_sum_2_mults(t4,g[4],g[6],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g6f4+g3
			_ff_sum_2_mults(t3,g[3],g[6],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g6f3+g2
			_ff_sum_2_mults(t2,g[2],g[6],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g6f2+g1
			_ff_sum_2_mults(t1,g[1],g[6],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g6f1+g0
			_ff_sum_2_mults(g[0],g[0],g[6],f[0],a);					// g0 = ag0-g6f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_7 (g,g,f);
	}
}


void ff_poly_xn_mod_d8 (ff_t g[8], unsigned long n, ff_t f[7])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 8 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 4 + ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_8 (g,g,f);
		} else {
			ff_poly_square_mod_8 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d8 (ff_t g[8], ff_t a, unsigned long n, ff_t f[7])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[7]);  _ff_set_zero(g[6]); _ff_set_zero(g[5]);_ff_set_zero(g[4]); _ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t7,a,g[7]); _ff_addto(t7,g[6]);					// g7 = ag7+g6
			_ff_sum_2_mults(t6,g[6],g[7],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g7f6+g5
			_ff_sum_2_mults(t5,g[5],g[7],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g7f5+g4
			_ff_sum_2_mults(t4,g[4],g[7],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g7f4+g3
			_ff_sum_2_mults(t3,g[3],g[7],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g7f3+g2
			_ff_sum_2_mults(t2,g[2],g[7],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g7f2+g1
			_ff_sum_2_mults(t1,g[1],g[7],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g7f1+g0
			_ff_sum_2_mults(g[0],g[0],g[7],f[0],a);					// g0 = ag0-g7f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_8 (g,g,f);
	}
}


void ff_poly_xn_mod_d10 (ff_t g[10], unsigned long n, ff_t f[9])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 10 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 4 + ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_10 (g,g,f);
		} else {
			ff_poly_square_mod_10 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d10 (ff_t g[10], ff_t a, unsigned long n, ff_t f[9])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[9]);  _ff_set_zero(g[8]); _ff_set_zero(g[7]);  _ff_set_zero(g[6]); _ff_set_zero(g[5]);_ff_set_zero(g[4]); _ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t9,a,g[9]); _ff_addto(t9,g[8]);					// g9 = ag9+g8
			_ff_sum_2_mults(t8,g[8],g[9],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g9f8+g7
			_ff_sum_2_mults(t7,g[7],g[9],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g9f7+g6
			_ff_sum_2_mults(t6,g[6],g[9],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g9f6+g5
			_ff_sum_2_mults(t5,g[5],g[9],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g9f5+g4
			_ff_sum_2_mults(t4,g[4],g[9],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g9f4+g3
			_ff_sum_2_mults(t3,g[3],g[9],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g9f3+g2
			_ff_sum_2_mults(t2,g[2],g[9],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g9f2+g1
			_ff_sum_2_mults(t1,g[1],g[9],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g9f1+g0
			_ff_sum_2_mults(g[0],g[0],g[9],f[0],a);					// g0 = ag0-g9f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_10 (g,g,f);
	}
}


void ff_poly_xn_mod_d11 (ff_t g[11], unsigned long n, ff_t f[10])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 11 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 4 + ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_11 (g,g,f);
		} else {
			ff_poly_square_mod_11 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d11 (ff_t g[11], ff_t a, unsigned long n, ff_t f[10])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]);  _ff_set_zero(g[6]); _ff_set_zero(g[5]);_ff_set_zero(g[4]); _ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t10,a,g[10]); _ff_addto(t10,g[9]);				// g10 = ag10+g9
			_ff_sum_2_mults(t9,g[9],g[10],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g10f9+g8
			_ff_sum_2_mults(t8,g[8],g[10],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g10f8+g7
			_ff_sum_2_mults(t7,g[7],g[10],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g10f7+g6
			_ff_sum_2_mults(t6,g[6],g[10],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g10f6+g5
			_ff_sum_2_mults(t5,g[5],g[10],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g10f5+g4
			_ff_sum_2_mults(t4,g[4],g[10],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g10f4+g3
			_ff_sum_2_mults(t3,g[3],g[10],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g10f3+g2
			_ff_sum_2_mults(t2,g[2],g[10],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g10f2+g1
			_ff_sum_2_mults(t1,g[1],g[10],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g10f1+g0
			_ff_sum_2_mults(g[0],g[0],g[10],f[0],a);					// g0 = ag0-g10f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_11 (g,g,f);
	}
}


void ff_poly_xn_mod_d12 (ff_t g[12], unsigned long n, ff_t f[11])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 12 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 4 + ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_12 (g,g,f);
		} else {
			ff_poly_square_mod_12 (g,g,f);
		}
		m >>= 1;
	}
}

void ff_poly_xpan_mod_d12 (ff_t g[12], ff_t a, unsigned long n, ff_t f[11])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]);  _ff_set_zero(g[6]); _ff_set_zero(g[5]);_ff_set_zero(g[4]); _ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t11,a,g[11]); _ff_addto(t11,g[10]);				// g11 = ag11+g10
			_ff_sum_2_mults(t10,g[10],g[11],f[10],a); _ff_addto(t10,g[9]);// g10 = ag10-g11f10+g9
			_ff_sum_2_mults(t9,g[9],g[11],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g11f9+g8
			_ff_sum_2_mults(t8,g[8],g[11],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g11f8+g7
			_ff_sum_2_mults(t7,g[7],g[11],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g11f7+g6
			_ff_sum_2_mults(t6,g[6],g[11],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g11f6+g5
			_ff_sum_2_mults(t5,g[5],g[11],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g11f5+g4
			_ff_sum_2_mults(t4,g[4],g[11],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g11f4+g3
			_ff_sum_2_mults(t3,g[3],g[11],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g11f3+g2
			_ff_sum_2_mults(t2,g[2],g[11],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g11f2+g1
			_ff_sum_2_mults(t1,g[1],g[11],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g11f1+g0
			_ff_sum_2_mults(g[0],g[0],g[11],f[0],a);				// g0 = ag0-g11f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10); _ff_set(g[11],t11);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_12 (g,g,f);
	}
}


void ff_poly_xn_mod_d13 (ff_t g[13], unsigned long n, ff_t f[12])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 13 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 4 + ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_13 (g,g,f);
		} else {
			ff_poly_square_mod_13 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d13 (ff_t g[13], ff_t a, unsigned long n, ff_t f[12])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]);  _ff_set_zero(g[6]); _ff_set_zero(g[5]);_ff_set_zero(g[4]); _ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t12,a,g[12]); _ff_addto(t12,g[11]);					// g12 = ag12+g11
			_ff_sum_2_mults(t11,g[11],g[12],f[11],a); _ff_addto(t11,g[10]);	// g11 = ag11-g12f11+g10
			_ff_sum_2_mults(t10,g[10],g[12],f[10],a); _ff_addto(t10,g[9]);	// g10 = ag10-g12f10+g9
			_ff_sum_2_mults(t9,g[9],g[12],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g12f9+g8
			_ff_sum_2_mults(t8,g[8],g[12],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g12f8+g7
			_ff_sum_2_mults(t7,g[7],g[12],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g12f7+g6
			_ff_sum_2_mults(t6,g[6],g[12],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g12f6+g5
			_ff_sum_2_mults(t5,g[5],g[12],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g12f5+g4
			_ff_sum_2_mults(t4,g[4],g[12],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g12f4+g3
			_ff_sum_2_mults(t3,g[3],g[12],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g12f3+g2
			_ff_sum_2_mults(t2,g[2],g[12],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g12f2+g1
			_ff_sum_2_mults(t1,g[1],g[12],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g12f1+g0
			_ff_sum_2_mults(g[0],g[0],g[12],f[0],a);					// g0 = ag0-g12f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10); _ff_set(g[11],t11); _ff_set(g[12],t12);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_13 (g,g,f);
	}
}


void ff_poly_xn_mod_d15 (ff_t g[15], unsigned long n, ff_t f[14])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 15 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 4 + ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_15 (g,g,f);
		} else {
			ff_poly_square_mod_15 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d15 (ff_t g[15], ff_t a, unsigned long n, ff_t f[14])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t14,a,g[14]); _ff_addto(t14,g[13]);					// g14 = ag14+g13
			_ff_sum_2_mults(t13,g[13],g[14],f[13],a); _ff_addto(t13,g[12]);	// g13 = ag13-g14f13+g12
			_ff_sum_2_mults(t12,g[12],g[14],f[12],a); _ff_addto(t12,g[11]);	// g12 = ag12-g14f12+g11
			_ff_sum_2_mults(t11,g[11],g[14],f[11],a); _ff_addto(t11,g[10]);	// g11 = ag11-g14f11+g10
			_ff_sum_2_mults(t10,g[10],g[14],f[10],a); _ff_addto(t10,g[9]);	// g10 = ag10-g14f10+g9
			_ff_sum_2_mults(t9,g[9],g[14],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g14f9+g8
			_ff_sum_2_mults(t8,g[8],g[14],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g14f8+g7
			_ff_sum_2_mults(t7,g[7],g[14],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g14f7+g6
			_ff_sum_2_mults(t6,g[6],g[14],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g14f6+g5
			_ff_sum_2_mults(t5,g[5],g[14],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g14f5+g4
			_ff_sum_2_mults(t4,g[4],g[14],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g14f4+g3
			_ff_sum_2_mults(t3,g[3],g[14],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g14f3+g2
			_ff_sum_2_mults(t2,g[2],g[14],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g14f2+g1
			_ff_sum_2_mults(t1,g[1],g[14],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g14f1+g0
			_ff_sum_2_mults(g[0],g[0],g[14],f[0],a);				// g0 = ag0-g14f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10); _ff_set(g[11],t11); _ff_set(g[12],t12); _ff_set(g[13],t13); _ff_set(g[14],t14);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_15 (g,g,f);
	}
}


void ff_poly_xn_mod_d17 (ff_t g[17], unsigned long n, ff_t f[16])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 17 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 8+ ((n&m)?4:0);
	m >>= 1;
	i += ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_17 (g,g,f);
		} else {
			ff_poly_square_mod_17 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d17 (ff_t g[17], ff_t a, unsigned long n, ff_t f[16])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t16,a,g[16]); _ff_addto(t16,g[15]);					// g16 = ag16+g15
			_ff_sum_2_mults(t15,g[15],g[16],f[15],a); _ff_addto(t15,g[14]);	// g15 = ag15-g16f15+g14
			_ff_sum_2_mults(t14,g[14],g[16],f[14],a); _ff_addto(t14,g[13]);	// g14 = ag14-g16f14+g13
			_ff_sum_2_mults(t13,g[13],g[16],f[13],a); _ff_addto(t13,g[12]);	// g13 = ag13-g16f13+g12
			_ff_sum_2_mults(t12,g[12],g[16],f[12],a); _ff_addto(t12,g[11]);	// g12 = ag12-g16f12+g11
			_ff_sum_2_mults(t11,g[11],g[16],f[11],a); _ff_addto(t11,g[10]);	// g11 = ag11-g16f11+g10
			_ff_sum_2_mults(t10,g[10],g[16],f[10],a); _ff_addto(t10,g[9]);	// g10 = ag10-g16f10+g9
			_ff_sum_2_mults(t9,g[9],g[16],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g16f9+g8
			_ff_sum_2_mults(t8,g[8],g[16],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g16f8+g7
			_ff_sum_2_mults(t7,g[7],g[16],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g16f7+g6
			_ff_sum_2_mults(t6,g[6],g[16],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g16f6+g5
			_ff_sum_2_mults(t5,g[5],g[16],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g16f5+g4
			_ff_sum_2_mults(t4,g[4],g[16],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g16f4+g3
			_ff_sum_2_mults(t3,g[3],g[16],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g16f3+g2
			_ff_sum_2_mults(t2,g[2],g[16],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g16f2+g1
			_ff_sum_2_mults(t1,g[1],g[16],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g16f1+g0
			_ff_sum_2_mults(g[0],g[0],g[16],f[0],a);				// g0 = ag0-g16f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10); _ff_set(g[11],t11); _ff_set(g[12],t12); _ff_set(g[13],t13); _ff_set(g[14],t14); _ff_set(g[15],t15); _ff_set(g[16],t16);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_17 (g,g,f);
	}
}


void ff_poly_xn_mod_d19 (ff_t g[19], unsigned long n, ff_t f[18])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 19 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 8+ ((n&m)?4:0);
	m >>= 1;
	i += ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_19 (g,g,f);
		} else {
			ff_poly_square_mod_19 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d19 (ff_t g[19], ff_t a, unsigned long n, ff_t f[18])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t18,a,g[18]); _ff_addto(t18,g[17]);					// g18 = ag18+g17
			_ff_sum_2_mults(t17,g[17],g[18],f[17],a); _ff_addto(t17,g[16]);	// g17 = ag17-g18f17+g16
			_ff_sum_2_mults(t16,g[16],g[18],f[16],a); _ff_addto(t16,g[15]);	// g16 = ag16-g18f16+g15
			_ff_sum_2_mults(t15,g[15],g[18],f[15],a); _ff_addto(t15,g[14]);	// g15 = ag15-g18f15+g14
			_ff_sum_2_mults(t14,g[14],g[18],f[14],a); _ff_addto(t14,g[13]);	// g14 = ag14-g18f14+g13
			_ff_sum_2_mults(t13,g[13],g[18],f[13],a); _ff_addto(t13,g[12]);	// g13 = ag13-g18f13+g12
			_ff_sum_2_mults(t12,g[12],g[18],f[12],a); _ff_addto(t12,g[11]);	// g12 = ag12-g18f12+g11
			_ff_sum_2_mults(t11,g[11],g[18],f[11],a); _ff_addto(t11,g[10]);	// g11 = ag11-g18f11+g10
			_ff_sum_2_mults(t10,g[10],g[18],f[10],a); _ff_addto(t10,g[9]);	// g10 = ag10-g18f10+g9
			_ff_sum_2_mults(t9,g[9],g[18],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g18f9+g8
			_ff_sum_2_mults(t8,g[8],g[18],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g18f8+g7
			_ff_sum_2_mults(t7,g[7],g[18],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g18f7+g6
			_ff_sum_2_mults(t6,g[6],g[18],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g18f6+g5
			_ff_sum_2_mults(t5,g[5],g[18],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g18f5+g4
			_ff_sum_2_mults(t4,g[4],g[18],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g18f4+g3
			_ff_sum_2_mults(t3,g[3],g[18],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g18f3+g2
			_ff_sum_2_mults(t2,g[2],g[18],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g18f2+g1
			_ff_sum_2_mults(t1,g[1],g[18],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g18f1+g0
			_ff_sum_2_mults(g[0],g[0],g[18],f[0],a);				// g0 = ag0-g18f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10); _ff_set(g[11],t11); _ff_set(g[12],t12); _ff_set(g[13],t13); _ff_set(g[14],t14); _ff_set(g[15],t15); _ff_set(g[16],t16); _ff_set(g[17],t17); _ff_set(g[18],t18);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_19 (g,g,f);
	}
}


void ff_poly_xn_mod_d23 (ff_t g[23], unsigned long n, ff_t f[22])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 23 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 8+ ((n&m)?4:0);
	m >>= 1;
	i += ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_23 (g,g,f);
		} else {
			ff_poly_square_mod_23 (g,g,f);
		}
		m >>= 1;
	}
}


void ff_poly_xpan_mod_d23 (ff_t g[23], ff_t a, unsigned long n, ff_t f[22])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t22,a,g[22]); _ff_addto(t22,g[21]);					// g22 = ag22+g21
			_ff_sum_2_mults(t21,g[21],g[22],f[21],a); _ff_addto(t21,g[20]);	// g21 = ag21-g22f21+g20
			_ff_sum_2_mults(t20,g[20],g[22],f[20],a); _ff_addto(t20,g[19]);	// g20 = ag20-g22f20+g19
			_ff_sum_2_mults(t19,g[19],g[22],f[19],a); _ff_addto(t19,g[18]);	// g19 = ag19-g22f19+g18
			_ff_sum_2_mults(t18,g[18],g[22],f[18],a); _ff_addto(t18,g[17]);	// g18 = ag18-g22f18+g17
			_ff_sum_2_mults(t17,g[17],g[22],f[17],a); _ff_addto(t17,g[16]);	// g17 = ag17-g22f17+g16
			_ff_sum_2_mults(t16,g[16],g[22],f[16],a); _ff_addto(t16,g[15]);	// g16 = ag16-g22f16+g15
			_ff_sum_2_mults(t15,g[15],g[22],f[15],a); _ff_addto(t15,g[14]);	// g15 = ag15-g22f15+g14
			_ff_sum_2_mults(t14,g[14],g[22],f[14],a); _ff_addto(t14,g[13]);	// g14 = ag14-g22f14+g13
			_ff_sum_2_mults(t13,g[13],g[22],f[13],a); _ff_addto(t13,g[12]);	// g13 = ag13-g22f13+g12
			_ff_sum_2_mults(t12,g[12],g[22],f[12],a); _ff_addto(t12,g[11]);	// g12 = ag12-g22f12+g11
			_ff_sum_2_mults(t11,g[11],g[22],f[11],a); _ff_addto(t11,g[10]);	// g11 = ag11-g22f11+g10
			_ff_sum_2_mults(t10,g[10],g[22],f[10],a); _ff_addto(t10,g[9]);	// g10 = ag10-g22f10+g9
			_ff_sum_2_mults(t9,g[9],g[22],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g22f9+g8
			_ff_sum_2_mults(t8,g[8],g[22],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g22f8+g7
			_ff_sum_2_mults(t7,g[7],g[22],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g22f7+g6
			_ff_sum_2_mults(t6,g[6],g[22],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g22f6+g5
			_ff_sum_2_mults(t5,g[5],g[22],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g22f5+g4
			_ff_sum_2_mults(t4,g[4],g[22],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g22f4+g3
			_ff_sum_2_mults(t3,g[3],g[22],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g22f3+g2
			_ff_sum_2_mults(t2,g[2],g[22],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g22f2+g1
			_ff_sum_2_mults(t1,g[1],g[22],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g22f1+g0
			_ff_sum_2_mults(g[0],g[0],g[22],f[0],a);				// g0 = ag0-g22f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10); _ff_set(g[11],t11); _ff_set(g[12],t12); _ff_set(g[13],t13); _ff_set(g[14],t14); _ff_set(g[15],t15); _ff_set(g[16],t16); _ff_set(g[17],t17); _ff_set(g[18],t18); _ff_set(g[19],t19); _ff_set(g[20],t20); _ff_set(g[21],t21); _ff_set(g[22],t22);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_23 (g,g,f);
	}
}


void ff_poly_xn_mod_d29 (ff_t g[29], unsigned long n, ff_t f[28])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[28]); _ff_set_zero(g[27]); _ff_set_zero(g[26]); _ff_set_zero(g[25]); _ff_set_zero(g[24]); _ff_set_zero(g[23]); _ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 29 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 8+ ((n&m)?4:0);
	m >>= 1;
	i += ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_29 (g,g,f);
		} else {
			ff_poly_square_mod_29 (g,g,f);
		}
		m >>= 1;
	}
}

void ff_poly_xpan_mod_d29 (ff_t g[29], ff_t a, unsigned long n, ff_t f[28])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[28]); _ff_set_zero(g[27]); _ff_set_zero(g[26]); _ff_set_zero(g[25]); _ff_set_zero(g[24]); _ff_set_zero(g[23]); _ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t28,a,g[28]); _ff_addto(t28,g[27]);					// g28 = ag28+g27
			_ff_sum_2_mults(t27,g[27],g[28],f[27],a); _ff_addto(t27,g[26]);	// g27 = ag27-g28f27+g26
			_ff_sum_2_mults(t26,g[26],g[28],f[26],a); _ff_addto(t26,g[25]);	// g26 = ag26-g28f26+g25
			_ff_sum_2_mults(t25,g[25],g[28],f[25],a); _ff_addto(t25,g[24]);	// g25 = ag25-g28f25+g24
			_ff_sum_2_mults(t24,g[24],g[28],f[24],a); _ff_addto(t24,g[23]);	// g24 = ag24-g28f24+g23
			_ff_sum_2_mults(t23,g[23],g[28],f[23],a); _ff_addto(t23,g[22]);	// g23 = ag23-g28f23+g22
			_ff_sum_2_mults(t22,g[22],g[28],f[22],a); _ff_addto(t22,g[21]);	// g22 = ag22-g28f22+g21
			_ff_sum_2_mults(t21,g[21],g[28],f[21],a); _ff_addto(t21,g[20]);	// g21 = ag21-g28f21+g20
			_ff_sum_2_mults(t20,g[20],g[28],f[20],a); _ff_addto(t20,g[19]);	// g20 = ag20-g28f20+g19
			_ff_sum_2_mults(t19,g[19],g[28],f[19],a); _ff_addto(t19,g[18]);	// g19 = ag19-g28f19+g18
			_ff_sum_2_mults(t18,g[18],g[28],f[18],a); _ff_addto(t18,g[17]);	// g18 = ag18-g28f18+g17
			_ff_sum_2_mults(t17,g[17],g[28],f[17],a); _ff_addto(t17,g[16]);	// g17 = ag17-g28f17+g16
			_ff_sum_2_mults(t16,g[16],g[28],f[16],a); _ff_addto(t16,g[15]);	// g16 = ag16-g28f16+g15
			_ff_sum_2_mults(t15,g[15],g[28],f[15],a); _ff_addto(t15,g[14]);	// g15 = ag15-g28f15+g14
			_ff_sum_2_mults(t14,g[14],g[28],f[14],a); _ff_addto(t14,g[13]);	// g14 = ag14-g28f14+g13
			_ff_sum_2_mults(t13,g[13],g[28],f[13],a); _ff_addto(t13,g[12]);	// g13 = ag13-g28f13+g12
			_ff_sum_2_mults(t12,g[12],g[28],f[12],a); _ff_addto(t12,g[11]);	// g12 = ag12-g28f12+g11
			_ff_sum_2_mults(t11,g[11],g[28],f[11],a); _ff_addto(t11,g[10]);	// g11 = ag11-g28f11+g10
			_ff_sum_2_mults(t10,g[10],g[28],f[10],a); _ff_addto(t10,g[9]);	// g10 = ag10-g28f10+g9
			_ff_sum_2_mults(t9,g[9],g[28],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g28f9+g8
			_ff_sum_2_mults(t8,g[8],g[28],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g28f8+g7
			_ff_sum_2_mults(t7,g[7],g[28],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g28f7+g6
			_ff_sum_2_mults(t6,g[6],g[28],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g28f6+g5
			_ff_sum_2_mults(t5,g[5],g[28],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g28f5+g4
			_ff_sum_2_mults(t4,g[4],g[28],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g28f4+g3
			_ff_sum_2_mults(t3,g[3],g[28],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g28f3+g2
			_ff_sum_2_mults(t2,g[2],g[28],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g28f2+g1
			_ff_sum_2_mults(t1,g[1],g[28],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g28f1+g0
			_ff_sum_2_mults(g[0],g[0],g[28],f[0],a);				// g0 = ag0-g28f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10); _ff_set(g[11],t11); _ff_set(g[12],t12); _ff_set(g[13],t13); _ff_set(g[14],t14); _ff_set(g[15],t15); _ff_set(g[16],t16); _ff_set(g[17],t17); _ff_set(g[18],t18); _ff_set(g[19],t19); _ff_set(g[20],t20); _ff_set(g[21],t21); _ff_set(g[22],t22); _ff_set(g[23],t23); _ff_set(g[24],t24); _ff_set(g[25],t25); _ff_set(g[26],t26); _ff_set(g[27],t27); _ff_set(g[28],t28);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_29 (g,g,f);
	}
}

void ff_poly_xn_mod_d31 (ff_t g[31], unsigned long n, ff_t f[30])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[30]); _ff_set_zero(g[29]);  _ff_set_zero(g[28]);  _ff_set_zero(g[27]); _ff_set_zero(g[26]); _ff_set_zero(g[25]); _ff_set_zero(g[24]); _ff_set_zero(g[23]); _ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 31 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 16+ ((n&m)?8:0);
	m >>= 1;
	i += ((n&m)?4:0);
	m >>= 1;
	i += ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	if ( i == 31 ) {
		for ( i = 0 ; i < 30 ; i++ ) _ff_set(g[i],f[i]);
	} else {
		_ff_set_one(g[i]);
	}
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_31 (g,g,f);
		} else {
			ff_poly_square_mod_31 (g,g,f);
		}
		m >>= 1;
	}
}

void ff_poly_xpan_mod_d31 (ff_t g[31], ff_t a, unsigned long n, ff_t f[30])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[30]); _ff_set_zero(g[29]);  _ff_set_zero(g[28]);  _ff_set_zero(g[27]); _ff_set_zero(g[26]); _ff_set_zero(g[25]); _ff_set_zero(g[24]); _ff_set_zero(g[23]); _ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t30,a,g[30]); _ff_addto(t30,g[29]);					// g30 = ag30+g29
			_ff_sum_2_mults(t29,g[29],g[30],f[29],a); _ff_addto(t29,g[28]);	// g29 = ag29-g30f29+g28
			_ff_sum_2_mults(t28,g[28],g[30],f[28],a); _ff_addto(t28,g[27]);	// g28 = ag28-g30f28+g27
			_ff_sum_2_mults(t27,g[27],g[30],f[27],a); _ff_addto(t27,g[26]);	// g27 = ag27-g30f27+g26
			_ff_sum_2_mults(t26,g[26],g[30],f[26],a); _ff_addto(t26,g[25]);	// g26 = ag26-g30f26+g25
			_ff_sum_2_mults(t25,g[25],g[30],f[25],a); _ff_addto(t25,g[24]);	// g25 = ag25-g30f25+g24
			_ff_sum_2_mults(t24,g[24],g[30],f[24],a); _ff_addto(t24,g[23]);	// g24 = ag24-g30f24+g23
			_ff_sum_2_mults(t23,g[23],g[30],f[23],a); _ff_addto(t23,g[22]);	// g23 = ag23-g30f23+g22
			_ff_sum_2_mults(t22,g[22],g[30],f[22],a); _ff_addto(t22,g[21]);	// g22 = ag22-g30f22+g21
			_ff_sum_2_mults(t21,g[21],g[30],f[21],a); _ff_addto(t21,g[20]);	// g21 = ag21-g30f21+g20
			_ff_sum_2_mults(t20,g[20],g[30],f[20],a); _ff_addto(t20,g[19]);	// g20 = ag20-g30f20+g19
			_ff_sum_2_mults(t19,g[19],g[30],f[19],a); _ff_addto(t19,g[18]);	// g19 = ag19-g30f19+g18
			_ff_sum_2_mults(t18,g[18],g[30],f[18],a); _ff_addto(t18,g[17]);	// g18 = ag18-g30f18+g17
			_ff_sum_2_mults(t17,g[17],g[30],f[17],a); _ff_addto(t17,g[16]);	// g17 = ag17-g30f17+g16
			_ff_sum_2_mults(t16,g[16],g[30],f[16],a); _ff_addto(t16,g[15]);	// g16 = ag16-g30f16+g15
			_ff_sum_2_mults(t15,g[15],g[30],f[15],a); _ff_addto(t15,g[14]);	// g15 = ag15-g30f15+g14
			_ff_sum_2_mults(t14,g[14],g[30],f[14],a); _ff_addto(t14,g[13]);	// g14 = ag14-g30f14+g13
			_ff_sum_2_mults(t13,g[13],g[30],f[13],a); _ff_addto(t13,g[12]);	// g13 = ag13-g30f13+g12
			_ff_sum_2_mults(t12,g[12],g[30],f[12],a); _ff_addto(t12,g[11]);	// g12 = ag12-g30f12+g11
			_ff_sum_2_mults(t11,g[11],g[30],f[11],a); _ff_addto(t11,g[10]);	// g11 = ag11-g30f11+g10
			_ff_sum_2_mults(t10,g[10],g[30],f[10],a); _ff_addto(t10,g[9]);	// g10 = ag10-g30f10+g9
			_ff_sum_2_mults(t9,g[9],g[30],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g30f9+g8
			_ff_sum_2_mults(t8,g[8],g[30],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g30f8+g7
			_ff_sum_2_mults(t7,g[7],g[30],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g30f7+g6
			_ff_sum_2_mults(t6,g[6],g[30],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g30f6+g5
			_ff_sum_2_mults(t5,g[5],g[30],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g30f5+g4
			_ff_sum_2_mults(t4,g[4],g[30],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g30f4+g3
			_ff_sum_2_mults(t3,g[3],g[30],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g30f3+g2
			_ff_sum_2_mults(t2,g[2],g[30],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g30f2+g1
			_ff_sum_2_mults(t1,g[1],g[30],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g30f1+g0
			_ff_sum_2_mults(g[0],g[0],g[30],f[0],a);				// g0 = ag0-g30f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10); _ff_set(g[11],t11); _ff_set(g[12],t12); _ff_set(g[13],t13); _ff_set(g[14],t14); _ff_set(g[15],t15); _ff_set(g[16],t16); _ff_set(g[17],t17); _ff_set(g[18],t18); _ff_set(g[19],t19); _ff_set(g[20],t20); _ff_set(g[21],t21); _ff_set(g[22],t22); _ff_set(g[23],t23); _ff_set(g[24],t24); _ff_set(g[25],t25); _ff_set(g[26],t26); _ff_set(g[27],t27);  _ff_set(g[28],t28); _ff_set(g[29],t29); _ff_set(g[30],t30);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_31 (g,g,f);
	}
}

#if FF_SUPERFAST

void ff_poly_xn_mod_d37 (ff_t g[37], unsigned long n, ff_t f[36])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[36]); _ff_set_zero(g[35]);  _ff_set_zero(g[34]);  _ff_set_zero(g[33]); _ff_set_zero(g[32]); _ff_set_zero(g[31]); _ff_set_zero(g[30]); _ff_set_zero(g[29]);  _ff_set_zero(g[28]);  _ff_set_zero(g[27]); _ff_set_zero(g[26]); _ff_set_zero(g[25]); _ff_set_zero(g[24]); _ff_set_zero(g[23]); _ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 37 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 16 + ((n&m)?8:0);
	m >>= 1;
	i += ((n&m)?4:0);
	m >>= 1;
	i += ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_37 (g,g,f);
		} else {
			ff_poly_square_mod_37 (g,g,f);
		}
		m >>= 1;
	}
}

void ff_poly_xpan_mod_d37 (ff_t g[37], ff_t a, unsigned long n, ff_t f[36])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[36]); _ff_set_zero(g[35]);  _ff_set_zero(g[34]);  _ff_set_zero(g[33]); _ff_set_zero(g[32]); _ff_set_zero(g[31]); _ff_set_zero(g[30]); _ff_set_zero(g[29]);  _ff_set_zero(g[28]);  _ff_set_zero(g[27]); _ff_set_zero(g[26]); _ff_set_zero(g[25]); _ff_set_zero(g[24]); _ff_set_zero(g[23]); _ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t36,a,g[36]); _ff_addto(t36,g[35]);					// g36 = ag36+g35
			_ff_sum_2_mults(t35,g[35],g[36],f[35],a); _ff_addto(t35,g[34]);	// g35 = ag35-g36f35+g34
			_ff_sum_2_mults(t34,g[34],g[36],f[34],a); _ff_addto(t34,g[33]);	// g34 = ag34-g36f34+g33
			_ff_sum_2_mults(t33,g[33],g[36],f[33],a); _ff_addto(t33,g[32]);	// g33 = ag33-g36f33+g32
			_ff_sum_2_mults(t32,g[32],g[36],f[32],a); _ff_addto(t32,g[31]);	// g32 = ag32-g36f32+g31
			_ff_sum_2_mults(t31,g[31],g[36],f[31],a); _ff_addto(t31,g[30]);	// g31 = ag31-g36f31+g30
			_ff_sum_2_mults(t30,g[30],g[36],f[30],a); _ff_addto(t30,g[29]);	// g30 = ag30-g36f30+g29
			_ff_sum_2_mults(t29,g[29],g[36],f[29],a); _ff_addto(t29,g[28]);	// g29 = ag29-g36f29+g28
			_ff_sum_2_mults(t28,g[28],g[36],f[28],a); _ff_addto(t28,g[27]);	// g28 = ag28-g36f28+g27
			_ff_sum_2_mults(t27,g[27],g[36],f[27],a); _ff_addto(t27,g[26]);	// g27 = ag27-g36f27+g26
			_ff_sum_2_mults(t26,g[26],g[36],f[26],a); _ff_addto(t26,g[25]);	// g26 = ag26-g36f26+g25
			_ff_sum_2_mults(t25,g[25],g[36],f[25],a); _ff_addto(t25,g[24]);	// g25 = ag25-g36f25+g24
			_ff_sum_2_mults(t24,g[24],g[36],f[24],a); _ff_addto(t24,g[23]);	// g24 = ag24-g36f24+g23
			_ff_sum_2_mults(t23,g[23],g[36],f[23],a); _ff_addto(t23,g[22]);	// g23 = ag23-g36f23+g22
			_ff_sum_2_mults(t22,g[22],g[36],f[22],a); _ff_addto(t22,g[21]);	// g22 = ag22-g36f22+g21
			_ff_sum_2_mults(t21,g[21],g[36],f[21],a); _ff_addto(t21,g[20]);	// g21 = ag21-g36f21+g20
			_ff_sum_2_mults(t20,g[20],g[36],f[20],a); _ff_addto(t20,g[19]);	// g20 = ag20-g36f20+g19
			_ff_sum_2_mults(t19,g[19],g[36],f[19],a); _ff_addto(t19,g[18]);	// g19 = ag19-g36f19+g18
			_ff_sum_2_mults(t18,g[18],g[36],f[18],a); _ff_addto(t18,g[17]);	// g18 = ag18-g36f18+g17
			_ff_sum_2_mults(t17,g[17],g[36],f[17],a); _ff_addto(t17,g[16]);	// g17 = ag17-g36f17+g16
			_ff_sum_2_mults(t16,g[16],g[36],f[16],a); _ff_addto(t16,g[15]);	// g16 = ag16-g36f16+g15
			_ff_sum_2_mults(t15,g[15],g[36],f[15],a); _ff_addto(t15,g[14]);	// g15 = ag15-g36f15+g14
			_ff_sum_2_mults(t14,g[14],g[36],f[14],a); _ff_addto(t14,g[13]);	// g14 = ag14-g36f14+g13
			_ff_sum_2_mults(t13,g[13],g[36],f[13],a); _ff_addto(t13,g[12]);	// g13 = ag13-g36f13+g12
			_ff_sum_2_mults(t12,g[12],g[36],f[12],a); _ff_addto(t12,g[11]);	// g12 = ag12-g36f12+g11
			_ff_sum_2_mults(t11,g[11],g[36],f[11],a); _ff_addto(t11,g[10]);	// g11 = ag11-g36f11+g10
			_ff_sum_2_mults(t10,g[10],g[36],f[10],a); _ff_addto(t10,g[9]);	// g10 = ag10-g36f10+g9
			_ff_sum_2_mults(t9,g[9],g[36],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g36f9+g8
			_ff_sum_2_mults(t8,g[8],g[36],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g36f8+g7
			_ff_sum_2_mults(t7,g[7],g[36],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g36f7+g6
			_ff_sum_2_mults(t6,g[6],g[36],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g36f6+g5
			_ff_sum_2_mults(t5,g[5],g[36],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g36f5+g4
			_ff_sum_2_mults(t4,g[4],g[36],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g36f4+g3
			_ff_sum_2_mults(t3,g[3],g[36],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g36f3+g2
			_ff_sum_2_mults(t2,g[2],g[36],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g36f2+g1
			_ff_sum_2_mults(t1,g[1],g[36],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g36f1+g0
			_ff_sum_2_mults(g[0],g[0],g[36],f[0],a);				// g0 = ag0-g36f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10); _ff_set(g[11],t11); _ff_set(g[12],t12); _ff_set(g[13],t13); _ff_set(g[14],t14); _ff_set(g[15],t15); _ff_set(g[16],t16); _ff_set(g[17],t17); _ff_set(g[18],t18); _ff_set(g[19],t19); _ff_set(g[20],t20); _ff_set(g[21],t21); _ff_set(g[22],t22); _ff_set(g[23],t23); _ff_set(g[24],t24); _ff_set(g[25],t25); _ff_set(g[26],t26); _ff_set(g[27],t27);  _ff_set(g[28],t28); _ff_set(g[29],t29);  _ff_set(g[30],t30); _ff_set(g[31],t31); _ff_set(g[32],t32); _ff_set(g[33],t33);  _ff_set(g[34],t34); _ff_set(g[35],t35); _ff_set(g[36],t36);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_37 (g,g,f);
	}
}

void ff_poly_xn_mod_d41 (ff_t g[41], unsigned long n, ff_t f[40])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[40]); _ff_set_zero(g[39]);  _ff_set_zero(g[38]);  _ff_set_zero(g[37]); _ff_set_zero(g[36]); _ff_set_zero(g[35]);  _ff_set_zero(g[34]);  _ff_set_zero(g[33]); _ff_set_zero(g[32]); _ff_set_zero(g[31]); _ff_set_zero(g[30]); _ff_set_zero(g[29]);  _ff_set_zero(g[28]);  _ff_set_zero(g[27]); _ff_set_zero(g[26]); _ff_set_zero(g[25]); _ff_set_zero(g[24]); _ff_set_zero(g[23]); _ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 41 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 16 + ((n&m)?8:0);
	m >>= 1;
	i += ((n&m)?4:0);
	m >>= 1;
	i += ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_41 (g,g,f);
		} else {
			ff_poly_square_mod_41 (g,g,f);
		}
		m >>= 1;
	}
}

void ff_poly_xpan_mod_d41 (ff_t g[41], ff_t a, unsigned long n, ff_t f[40])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[40]); _ff_set_zero(g[39]);  _ff_set_zero(g[38]);  _ff_set_zero(g[37]); _ff_set_zero(g[36]); _ff_set_zero(g[35]);  _ff_set_zero(g[34]);  _ff_set_zero(g[33]); _ff_set_zero(g[32]); _ff_set_zero(g[31]); _ff_set_zero(g[30]); _ff_set_zero(g[29]);  _ff_set_zero(g[28]);  _ff_set_zero(g[27]); _ff_set_zero(g[26]); _ff_set_zero(g[25]); _ff_set_zero(g[24]); _ff_set_zero(g[23]); _ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t40,a,g[40]); _ff_addto(t40,g[39]);					// g40 = ag40+g39
			_ff_sum_2_mults(t39,g[39],g[40],f[39],a); _ff_addto(t39,g[38]);	// g39 = ag39-g40f39+g38
			_ff_sum_2_mults(t38,g[38],g[40],f[38],a); _ff_addto(t38,g[37]);	// g38 = ag38-g40f38+g37
			_ff_sum_2_mults(t37,g[37],g[40],f[37],a); _ff_addto(t37,g[36]);	// g37 = ag37-g40f37+g36
			_ff_sum_2_mults(t36,g[36],g[40],f[36],a); _ff_addto(t36,g[35]);	// g36 = ag36-g40f36+g35
			_ff_sum_2_mults(t35,g[35],g[40],f[35],a); _ff_addto(t35,g[34]);	// g35 = ag35-g40f35+g34
			_ff_sum_2_mults(t34,g[34],g[40],f[34],a); _ff_addto(t34,g[33]);	// g34 = ag34-g40f34+g33
			_ff_sum_2_mults(t33,g[33],g[40],f[33],a); _ff_addto(t33,g[32]);	// g33 = ag33-g40f33+g32
			_ff_sum_2_mults(t32,g[32],g[40],f[32],a); _ff_addto(t32,g[31]);	// g32 = ag32-g40f32+g31
			_ff_sum_2_mults(t31,g[31],g[40],f[31],a); _ff_addto(t31,g[30]);	// g31 = ag31-g40f31+g30
			_ff_sum_2_mults(t30,g[30],g[40],f[30],a); _ff_addto(t30,g[29]);	// g30 = ag30-g40f30+g29
			_ff_sum_2_mults(t29,g[29],g[40],f[29],a); _ff_addto(t29,g[28]);	// g29 = ag29-g40f29+g28
			_ff_sum_2_mults(t28,g[28],g[40],f[28],a); _ff_addto(t28,g[27]);	// g28 = ag28-g40f28+g27
			_ff_sum_2_mults(t27,g[27],g[40],f[27],a); _ff_addto(t27,g[26]);	// g27 = ag27-g40f27+g26
			_ff_sum_2_mults(t26,g[26],g[40],f[26],a); _ff_addto(t26,g[25]);	// g26 = ag26-g40f26+g25
			_ff_sum_2_mults(t25,g[25],g[40],f[25],a); _ff_addto(t25,g[24]);	// g25 = ag25-g40f25+g24
			_ff_sum_2_mults(t24,g[24],g[40],f[24],a); _ff_addto(t24,g[23]);	// g24 = ag24-g40f24+g23
			_ff_sum_2_mults(t23,g[23],g[40],f[23],a); _ff_addto(t23,g[22]);	// g23 = ag23-g40f23+g22
			_ff_sum_2_mults(t22,g[22],g[40],f[22],a); _ff_addto(t22,g[21]);	// g22 = ag22-g40f22+g21
			_ff_sum_2_mults(t21,g[21],g[40],f[21],a); _ff_addto(t21,g[20]);	// g21 = ag21-g40f21+g20
			_ff_sum_2_mults(t20,g[20],g[40],f[20],a); _ff_addto(t20,g[19]);	// g20 = ag20-g40f20+g19
			_ff_sum_2_mults(t19,g[19],g[40],f[19],a); _ff_addto(t19,g[18]);	// g19 = ag19-g40f19+g18
			_ff_sum_2_mults(t18,g[18],g[40],f[18],a); _ff_addto(t18,g[17]);	// g18 = ag18-g40f18+g17
			_ff_sum_2_mults(t17,g[17],g[40],f[17],a); _ff_addto(t17,g[16]);	// g17 = ag17-g40f17+g16
			_ff_sum_2_mults(t16,g[16],g[40],f[16],a); _ff_addto(t16,g[15]);	// g16 = ag16-g40f16+g15
			_ff_sum_2_mults(t15,g[15],g[40],f[15],a); _ff_addto(t15,g[14]);	// g15 = ag15-g40f15+g14
			_ff_sum_2_mults(t14,g[14],g[40],f[14],a); _ff_addto(t14,g[13]);	// g14 = ag14-g40f14+g13
			_ff_sum_2_mults(t13,g[13],g[40],f[13],a); _ff_addto(t13,g[12]);	// g13 = ag13-g40f13+g12
			_ff_sum_2_mults(t12,g[12],g[40],f[12],a); _ff_addto(t12,g[11]);	// g12 = ag12-g40f12+g11
			_ff_sum_2_mults(t11,g[11],g[40],f[11],a); _ff_addto(t11,g[10]);	// g11 = ag11-g40f11+g10
			_ff_sum_2_mults(t10,g[10],g[40],f[10],a); _ff_addto(t10,g[9]);	// g10 = ag10-g40f10+g9
			_ff_sum_2_mults(t9,g[9],g[40],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g40f9+g8
			_ff_sum_2_mults(t8,g[8],g[40],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g40f8+g7
			_ff_sum_2_mults(t7,g[7],g[40],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g40f7+g6
			_ff_sum_2_mults(t6,g[6],g[40],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g40f6+g5
			_ff_sum_2_mults(t5,g[5],g[40],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g40f5+g4
			_ff_sum_2_mults(t4,g[4],g[40],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g40f4+g3
			_ff_sum_2_mults(t3,g[3],g[40],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g40f3+g2
			_ff_sum_2_mults(t2,g[2],g[40],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g40f2+g1
			_ff_sum_2_mults(t1,g[1],g[40],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g40f1+g0
			_ff_sum_2_mults(g[0],g[0],g[40],f[0],a);				// g0 = ag0-g40f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10); _ff_set(g[11],t11); _ff_set(g[12],t12); _ff_set(g[13],t13); _ff_set(g[14],t14); _ff_set(g[15],t15); _ff_set(g[16],t16); _ff_set(g[17],t17); _ff_set(g[18],t18); _ff_set(g[19],t19); _ff_set(g[20],t20); _ff_set(g[21],t21); _ff_set(g[22],t22); _ff_set(g[23],t23); _ff_set(g[24],t24); _ff_set(g[25],t25); _ff_set(g[26],t26); _ff_set(g[27],t27);  _ff_set(g[28],t28); _ff_set(g[29],t29);  _ff_set(g[30],t30); _ff_set(g[31],t31); _ff_set(g[32],t32); _ff_set(g[33],t33);  _ff_set(g[34],t34); _ff_set(g[35],t35); _ff_set(g[36],t36); _ff_set(g[37],t37);  _ff_set(g[38],t38); _ff_set(g[39],t39); _ff_set(g[40],t40);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_41 (g,g,f);
	}
}

void ff_poly_xn_mod_d43 (ff_t g[43], unsigned long n, ff_t f[42])
{
	register unsigned long m;
	register int i;

	_ff_set_zero(g[42]); _ff_set_zero(g[41]);  _ff_set_zero(g[40]); _ff_set_zero(g[39]);  _ff_set_zero(g[38]);  _ff_set_zero(g[37]); _ff_set_zero(g[36]); _ff_set_zero(g[35]);  _ff_set_zero(g[34]);  _ff_set_zero(g[33]); _ff_set_zero(g[32]); _ff_set_zero(g[31]); _ff_set_zero(g[30]); _ff_set_zero(g[29]);  _ff_set_zero(g[28]);  _ff_set_zero(g[27]); _ff_set_zero(g[26]); _ff_set_zero(g[25]); _ff_set_zero(g[24]); _ff_set_zero(g[23]); _ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( n < 43 ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 16 + ((n&m)?8:0);
	m >>= 1;
	i += ((n&m)?4:0);
	m >>= 1;
	i += ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_43 (g,g,f);
		} else {
			ff_poly_square_mod_43 (g,g,f);
		}
		m >>= 1;
	}
}

void ff_poly_xpan_mod_d43 (ff_t g[43], ff_t a, unsigned long n, ff_t f[42])
{
	register ff_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42;
	register unsigned long m;
	register int i;

	_ff_set_zero(g[42]); _ff_set_zero(g[41]);  _ff_set_zero(g[40]); _ff_set_zero(g[39]);  _ff_set_zero(g[38]);  _ff_set_zero(g[37]); _ff_set_zero(g[36]); _ff_set_zero(g[35]);  _ff_set_zero(g[34]);  _ff_set_zero(g[33]); _ff_set_zero(g[32]); _ff_set_zero(g[31]); _ff_set_zero(g[30]); _ff_set_zero(g[29]);  _ff_set_zero(g[28]);  _ff_set_zero(g[27]); _ff_set_zero(g[26]); _ff_set_zero(g[25]); _ff_set_zero(g[24]); _ff_set_zero(g[23]); _ff_set_zero(g[22]); _ff_set_zero(g[21]); _ff_set_zero(g[20]); _ff_set_zero(g[19]); _ff_set_zero(g[18]); _ff_set_zero(g[17]); _ff_set_zero(g[16]); _ff_set_zero(g[15]); _ff_set_zero(g[14]); _ff_set_zero(g[13]); _ff_set_zero(g[12]); _ff_set_zero(g[11]); _ff_set_zero(g[10]); _ff_set_zero(g[9]); _ff_set_zero(g[8]); _ff_set_zero(g[7]); _ff_set_zero(g[6]); _ff_set_zero(g[5]); _ff_set_zero(g[4]);  _ff_set_zero(g[3]);_ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_zero(g[0]);
	if ( ! n ) { _ff_set_zero(g[2]); _ff_set_zero(g[1]); _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_zero(g[2]); _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_mult(t42,a,g[42]); _ff_addto(t42,g[41]);					// g42 = ag42+g41
			_ff_sum_2_mults(t41,g[41],g[42],f[41],a); _ff_addto(t41,g[40]);	// g41 = ag41-g42f41+g40
			_ff_sum_2_mults(t40,g[40],g[42],f[40],a); _ff_addto(t40,g[39]);	// g40 = ag40-g42f40+g39
			_ff_sum_2_mults(t39,g[39],g[42],f[39],a); _ff_addto(t39,g[38]);	// g39 = ag39-g42f39+g38
			_ff_sum_2_mults(t38,g[38],g[42],f[38],a); _ff_addto(t38,g[37]);	// g38 = ag38-g42f38+g37
			_ff_sum_2_mults(t37,g[37],g[42],f[37],a); _ff_addto(t37,g[36]);	// g37 = ag37-g42f37+g36
			_ff_sum_2_mults(t36,g[36],g[42],f[36],a); _ff_addto(t36,g[35]);	// g36 = ag36-g42f36+g35
			_ff_sum_2_mults(t35,g[35],g[42],f[35],a); _ff_addto(t35,g[34]);	// g35 = ag35-g42f35+g34
			_ff_sum_2_mults(t34,g[34],g[42],f[34],a); _ff_addto(t34,g[33]);	// g34 = ag34-g42f34+g33
			_ff_sum_2_mults(t33,g[33],g[42],f[33],a); _ff_addto(t33,g[32]);	// g33 = ag33-g42f33+g32
			_ff_sum_2_mults(t32,g[32],g[42],f[32],a); _ff_addto(t32,g[31]);	// g32 = ag32-g42f32+g31
			_ff_sum_2_mults(t31,g[31],g[42],f[31],a); _ff_addto(t31,g[30]);	// g31 = ag31-g42f31+g30
			_ff_sum_2_mults(t30,g[30],g[42],f[30],a); _ff_addto(t30,g[29]);	// g30 = ag30-g42f30+g29
			_ff_sum_2_mults(t29,g[29],g[42],f[29],a); _ff_addto(t29,g[28]);	// g29 = ag29-g42f29+g28
			_ff_sum_2_mults(t28,g[28],g[42],f[28],a); _ff_addto(t28,g[27]);	// g28 = ag28-g42f28+g27
			_ff_sum_2_mults(t27,g[27],g[42],f[27],a); _ff_addto(t27,g[26]);	// g27 = ag27-g42f27+g26
			_ff_sum_2_mults(t26,g[26],g[42],f[26],a); _ff_addto(t26,g[25]);	// g26 = ag26-g42f26+g25
			_ff_sum_2_mults(t25,g[25],g[42],f[25],a); _ff_addto(t25,g[24]);	// g25 = ag25-g42f25+g24
			_ff_sum_2_mults(t24,g[24],g[42],f[24],a); _ff_addto(t24,g[23]);	// g24 = ag24-g42f24+g23
			_ff_sum_2_mults(t23,g[23],g[42],f[23],a); _ff_addto(t23,g[22]);	// g23 = ag23-g42f23+g22
			_ff_sum_2_mults(t22,g[22],g[42],f[22],a); _ff_addto(t22,g[21]);	// g22 = ag22-g42f22+g21
			_ff_sum_2_mults(t21,g[21],g[42],f[21],a); _ff_addto(t21,g[20]);	// g21 = ag21-g42f21+g20
			_ff_sum_2_mults(t20,g[20],g[42],f[20],a); _ff_addto(t20,g[19]);	// g20 = ag20-g42f20+g19
			_ff_sum_2_mults(t19,g[19],g[42],f[19],a); _ff_addto(t19,g[18]);	// g19 = ag19-g42f19+g18
			_ff_sum_2_mults(t18,g[18],g[42],f[18],a); _ff_addto(t18,g[17]);	// g18 = ag18-g42f18+g17
			_ff_sum_2_mults(t17,g[17],g[42],f[17],a); _ff_addto(t17,g[16]);	// g17 = ag17-g42f17+g16
			_ff_sum_2_mults(t16,g[16],g[42],f[16],a); _ff_addto(t16,g[15]);	// g16 = ag16-g42f16+g15
			_ff_sum_2_mults(t15,g[15],g[42],f[15],a); _ff_addto(t15,g[14]);	// g15 = ag15-g42f15+g14
			_ff_sum_2_mults(t14,g[14],g[42],f[14],a); _ff_addto(t14,g[13]);	// g14 = ag14-g42f14+g13
			_ff_sum_2_mults(t13,g[13],g[42],f[13],a); _ff_addto(t13,g[12]);	// g13 = ag13-g42f13+g12
			_ff_sum_2_mults(t12,g[12],g[42],f[12],a); _ff_addto(t12,g[11]);	// g12 = ag12-g42f12+g11
			_ff_sum_2_mults(t11,g[11],g[42],f[11],a); _ff_addto(t11,g[10]);	// g11 = ag11-g42f11+g10
			_ff_sum_2_mults(t10,g[10],g[42],f[10],a); _ff_addto(t10,g[9]);	// g10 = ag10-g42f10+g9
			_ff_sum_2_mults(t9,g[9],g[42],f[9],a); _ff_addto(t9,g[8]);	// g9 = ag9-g42f9+g8
			_ff_sum_2_mults(t8,g[8],g[42],f[8],a); _ff_addto(t8,g[7]);	// g8 = ag8-g42f8+g7
			_ff_sum_2_mults(t7,g[7],g[42],f[7],a); _ff_addto(t7,g[6]);	// g7 = ag7-g42f7+g6
			_ff_sum_2_mults(t6,g[6],g[42],f[6],a); _ff_addto(t6,g[5]);	// g6 = ag6-g42f6+g5
			_ff_sum_2_mults(t5,g[5],g[42],f[5],a); _ff_addto(t5,g[4]);	// g5 = ag5-g42f5+g4
			_ff_sum_2_mults(t4,g[4],g[42],f[4],a); _ff_addto(t4,g[3]);	// g4 = ag4-g42f4+g3
			_ff_sum_2_mults(t3,g[3],g[42],f[3],a); _ff_addto(t3,g[2]);	// g3 = ag3-g42f3+g2
			_ff_sum_2_mults(t2,g[2],g[42],f[2],a); _ff_addto(t2,g[1]);	// g2 = ag2-g42f2+g1
			_ff_sum_2_mults(t1,g[1],g[42],f[1],a); _ff_addto(t1,g[0]);	// g1 = ag1-g42f1+g0
			_ff_sum_2_mults(g[0],g[0],g[42],f[0],a);				// g0 = ag0-g42f0
			_ff_set(g[1],t1);  _ff_set(g[2],t2); _ff_set(g[3],t3); _ff_set(g[4],t4); _ff_set(g[5],t5); _ff_set(g[6],t6); _ff_set(g[7],t7); _ff_set(g[8],t8); _ff_set(g[9],t9); _ff_set(g[10],t10); _ff_set(g[11],t11); _ff_set(g[12],t12); _ff_set(g[13],t13); _ff_set(g[14],t14); _ff_set(g[15],t15); _ff_set(g[16],t16); _ff_set(g[17],t17); _ff_set(g[18],t18); _ff_set(g[19],t19); _ff_set(g[20],t20); _ff_set(g[21],t21); _ff_set(g[22],t22); _ff_set(g[23],t23); _ff_set(g[24],t24); _ff_set(g[25],t25); _ff_set(g[26],t26); _ff_set(g[27],t27);  _ff_set(g[28],t28); _ff_set(g[29],t29);  _ff_set(g[30],t30); _ff_set(g[31],t31); _ff_set(g[32],t32); _ff_set(g[33],t33);  _ff_set(g[34],t34); _ff_set(g[35],t35); _ff_set(g[36],t36); _ff_set(g[37],t37);  _ff_set(g[38],t38); _ff_set(g[39],t39); _ff_set(g[40],t40);_ff_set(g[41],t41); _ff_set(g[42],t42);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_43 (g,g,f);
	}
}

#endif

// Computes x^n modulo a depressed poly (x^d and x^{d-1} coeff are assumed 1 and 0, others are *negated*).  Assumes n < 2^63 and d > 8
void ff_poly_xn_mod_small (ff_t g[], unsigned long n, ff_t f[], int d)
{
	register unsigned long m;
	register int i;

	if ( d < 8 || d > FF_POLY_SMALL_DEGREE ) { printf ("Degree must be at least 8 and less than FF_POLY_SMALL_DEGREE=%d in ff_poly_xn_mod_small\n", FF_POLY_SMALL_DEGREE); abort(); }

	for  ( i =  0 ; i < d ; i++ ) _ff_set_zero(g[i]);
	if ( n < d ) { _ff_set_one(g[n]); return; }
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	i = 4 + ((n&m)?2:0);
	m >>= 1;
	i += ((n&m)?1:0);
	m >>= 1;
	_ff_set_one(g[i]);
	while (m) {
		if ( (n&m)) {
			ff_poly_square_mult_x_mod_small (g,g,f,d);
		} else {
			ff_poly_square_mod_small (g,g,f,d);
		}
		m >>= 1;
	}
}

// Computes (x+a)^n modulo a depressed deg d poly (x^d and x^{d-1} coeff are assumed 1 and 0, others are *negated*).  Assumes n < 2^63 and d > 8
void ff_poly_xpan_mod_small (ff_t g[], ff_t a, unsigned long n, ff_t f[], int d)
{
	register ff_t gtop, t0;
	register unsigned long m;
	register int i;

	if ( d < 8 || d > FF_POLY_SMALL_DEGREE ) { printf ("Degree must be at least 8 and less than FF_POLY_SMALL_DEGREE=%d in ff_poly_xn_mod_small\n", FF_POLY_SMALL_DEGREE); abort(); }
	for  ( i =  0 ; i < d ; i++ ) _ff_set_zero(g[i]);
	if ( ! n ) { _ff_set_one(g[0]); return;}
	if ( n== 1 ) { _ff_set_one(g[1]); _ff_set(g[0],a); return; }
	_ff_set_one(g[2]); _ff_add(g[1],a,a);  _ff_square(g[0],a);			// start with (x+a)^2 (for now)
	if ( n == 2 ) return;
	i = _asm_highbit(n)-1;
	m = 1UL<<i;
	for(;;) {
		if ( (n&m)) {
			_ff_set(gtop,g[d-1]);
			_ff_mult(t0,a,gtop);  _ff_add(g[d-1],t0,g[d-2]);
			for ( i = d-2 ; i > 0 ; i-- )
				{ _ff_sum_2_mults(t0,g[i],gtop,f[i],a); _ff_add(g[i],t0,g[i-1]); }
			_ff_sum_2_mults(g[0],g[0],gtop,f[0],a);
		}
		m >>= 1;
		if ( ! m ) break;
		ff_poly_square_mod_small (g,g,f,d);
	}
}

// assumes d_f >= d_g, replaces f by a linear combination of f and g that has lower degree than f, namely g[d_g]*f - f[d_f]*x^{d_f-d_g}*g
static inline void ff_poly_reduce_inplace (ff_t f[], int *d_f, ff_t g[], int d_g)
{
	register ff_t t0,t1;
	register int i, j;

	j = *d_f-d_g;
	_ff_set(t0,g[d_g]);
	for ( i = 0 ; i < j ; i++ ) ff_mult(f[i],f[i],t0);
	_ff_neg(t1,f[*d_f]);
	for ( ; i < *d_f ; i++ ) _ff_sum_2_mults (f[i],t0,t1,g[i-j],f[i]);
	*d_f = ff_poly_degree(f,i-1);
}

static inline void _swap_polys (ff_t **f, int *d_f, ff_t **g, int *d_g)
{
	register ff_t *h;
	register int d_h;
	h = *f; d_h = *d_f;  *f = *g; *d_f = *d_g; *g = h; *d_g = d_h;
}

// the poly with lesser degree (or both) must be monic, destroys f and g
void ff_poly_gcd_small (ff_t h[], int *d_h, ff_t f[], int d_f, ff_t g[], int d_g)
{
	register ff_t t0, t1, t2, t3;
	register int i;

	if ( d_g < d_f )  _swap_polys(&f,&d_f,&g,&d_g);
	if ( d_g == d_f &&  ! _ff_one(f[d_f]) ) _swap_polys(&f,&d_f,&g,&d_g);
	ff_poly_mod_small_inplace (g, d_g, f, d_f);
	d_g = ff_poly_degree(g,d_f-1);
	for (;;) {
		if ( d_f < d_g )  _swap_polys(&f,&d_f,&g,&d_g);
		if ( d_g < 1 ) break;
		if ( d_f !=	d_g+1 ) { ff_poly_reduce_inplace (f, &d_f, g, d_g); continue; }
		_ff_neg(t0,g[d_g]);  _ff_sum_2_mults(t1,f[d_g],f[d_f],g[d_g-1],t0); 					// t1 is alpha
		_ff_mult(t2,t0,f[d_f]);  _ff_square(t3,t0);
		for ( i = d_g-1 ; i ; i-- ) _ff_sum_3_mults(f[i],f[i],t1,t2,g[i-1],g[i],t3);
		_ff_sum_2_mults(f[0],f[0],t1,g[0],t3);
		d_f = ff_poly_degree(f,d_g-1);
	}
	if ( d_g < 0 ) { ff_poly_copy (h, d_h, f, d_f); return; }
	_ff_set_one(h[0]); *d_h = 0;
// for normal degree sequence d_f >= d_g, uses (d_f-d_g+1)d_g + (3d_g^2+5d_g)/2 mults
}

// call faster inline code for small cases, but check for failure that can arise from leading zero coefficients
// assumes d_f = d_g+1
// This code doesn't actually speed things up!
/*
int ff_poly_gcd_linear_tiny (ff_t h[], ff_t f[], int d_f, ff_t g[])
{
	ff_t t[2];

	switch(d_f) {
	case 7: ff_poly_gcd_linear_7_6_reg(t,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],g[0],g[1],g[2],g[3],g[4],g[5],g[6]); break;
	case 6: ff_poly_gcd_linear_6_5_reg(t,f[0],f[1],f[2],f[3],f[4],f[5],f[6],g[0],g[1],g[2],g[3],g[4],g[5]); break;
	case 5: ff_poly_gcd_linear_5_4_reg(t,f[0],f[1],f[2],f[3],f[4],f[5],g[0],g[1],g[2],g[3],g[4]); break;
	case 4: ff_poly_gcd_linear_4_3_reg(t,f[0],f[1],f[2],f[3],f[4],g[0],g[1],g[2],g[3]); break;
	case 3: ff_poly_gcd_linear_3_2_reg(t,f[0],f[1],f[2],f[3],g[0],g[1],g[2]); break;
	}
	if  ( _ff_zero(t[1]) ) return 0;
	_ff_set(h[1],t[1]); _ff_set(h[0],t[0]);
	return 1;
}
*/

static inline void ff_poly_reduce3m_inplace (ff_t *f, ff_t *g, int d)	// f deg d+1 monic, g deg d >= 4 need not be monic.  Replaces f by linear combo of f and g with deg d-2 (possibly trivial, so caller should check lc)
{
	ff_t f1, f2, g2, g3, g32, s0, s1, t0, t1, t2, w1, w2, w3, w4;		// this is more registers than we need, but handy for debugging
	int i;

//ff_poly_print(f,d+1);
//ff_poly_print(g,d);
	 _ff_set(g3,g[d]);  _ff_square(g32,g3);
	_ff_mult (t0,f[d],g3); _ff_sub(w3,g[d-1],t0); _ff_neg(w4,g3);						// w3 = -f3, w4=-f4g3 = -g3 since f4=1
	_ff_sum_3_mults(f2,f[d-1],w3,w4,g[d-2],g[d-1],g32);  _ff_neg(w2,f2);
	_ff_sum_3_mults(f1,f[d-2],w3,w4,g[d-3],g[d-2],g32);  _ff_neg(w1,f1);
//_ff_neg(f3,w3);
//printf ("f3=%ld, f2=%ld, f1=%ld\n", _ff_get_ui(f3), _ff_get_ui(f2), _ff_get_ui(f1));
	_ff_sum_2_mults(g2,g[d-1],g3,w1,f2);
//printf ("g2=%ld\n", _ff_get_ui(g2));
	 _ff_neg(s0,g2);  _ff_mult(s1,w2,g3);
	_ff_mult(t2,s1,w4);
	_ff_sum_2_mults(t1,s0,s1,w3,w4);
	_ff_sum_2_mults(t0,s0,w2,w2,w3);
	ff_mult(s1,s1,g32);
	ff_mult(s0,s0,g32);
//printf ("t2=%ld, t1=%ld, t0=%ld, s1=%ld, s0=%ld\n", _ff_get_ui(t2), _ff_get_ui(t1), _ff_get_ui(t0), _ff_get_ui(s1), _ff_get_ui(s0));
	// 18M + 12A (11 redc) to compute s and t
	
	// we now need to replace f by (t2 t1 t0)g + (s1 s0)f of degree d-2
	_ff_sum_2_mults(w1,s0,t0,g[0],f[0]);
	_ff_sum_4_mults(w2,s0,s1,t0,t1,g[0],g[1],f[0],f[1]);  _ff_set(f[0],w1);
	for ( i = 1 ; i < d-2 ; i++ ) {
		_ff_sum_5_mults (w1,s0,s1,t0,t1,t2,g[i-1],g[i],g[i+1],f[i],f[i+1]);
		_ff_set(f[i],w2); _ff_set(w2,w1);
	}
	_ff_set(f[i],w2);
//ff_poly_print(f,d-2);
	
	// (5d+9)M (d+8) redc
}

static inline void ff_poly_reduce4_inplace (ff_t *f, ff_t *g, int d)	// f deg d+2, g deg d >= 4, neither need be monic.  Replaces f by linear combo of f and g with deg d-2 (possibly trivial, so caller should check lc)
{
	ff_t f1, f2, f3, f4, g2, g3, g32, g33, s0, s1, t0, t1, t2, t3, w1, w2, w3, w4, w5;		// this is way more registers than we need, but handy for debugging
	int i;

//ff_poly_print(f,d+2);
//ff_poly_print(g,d);
	_ff_neg(w5,f[d+2]);  _ff_set(g3,g[d]);
	_ff_sum_2_mults(f4,f[d+1],w5,g[d-1],g3);
	ff_mult(w5,w5,g3); _ff_square(g32,g3); _ff_neg(w4,f4);
	_ff_sum_3_mults(f3,f[d],w4,w5,g[d-2],g[d-1],g32);
	ff_mult(w5,w5,g3); ff_mult(w4,w4,g3); _ff_neg(w3,f3); _ff_mult(g33,g32,g3);
	_ff_sum_4_mults(f2,f[d-1],w3,w4,w5,g[d-3],g[d-2],g[d-1],g33);  _ff_neg(w2,f2);
	_ff_sum_4_mults(f1,f[d-2],w3,w4,w5,g[d-4],g[d-3],g[d-2],g33);  _ff_neg(w1,f1);
//printf ("f4=%ld, f3=%ld, f2=%ld, f1=%ld\n", _ff_get_ui(f4), _ff_get_ui(f3), _ff_get_ui(f2), _ff_get_ui(f1));
	_ff_sum_2_mults(g2,g[d-1],g3,w1,f2);
//printf ("g2=%ld\n", _ff_get_ui(g2));
	 _ff_neg(s0,g2);  _ff_mult(s1,w2,g3);
	_ff_mult(t3,s1,w5);
	_ff_sum_2_mults(t2,s0,s1,w4,w5);	
	_ff_sum_2_mults(t1,s0,s1,w3,w4);
	_ff_sum_2_mults(t0,s0,w2,w2,w3);
	ff_mult(s1,s1,g33);
	ff_mult(s0,s0,g33);
//printf ("t3=%ld, t2=%ld, t1=%ld, t0=%ld, s1=%ld, s0=%ld\n", _ff_get_ui(t3), _ff_get_ui(t2), _ff_get_ui(t1), _ff_get_ui(t0), _ff_get_ui(s1), _ff_get_ui(s0));
	// 30M + 19A (17 redc) to compute s and t
	
	// we now need to replace f by (t3 t2 t1 t0)g + (s1 s0)f of degree d-2
	_ff_sum_2_mults(w1,s0,t0,g[0],f[0]);
	_ff_sum_4_mults(w2,s0,s1,t0,t1,g[0],g[1],f[0],f[1]);  _ff_set(f[0],w1);
	_ff_sum_5_mults(w1,s0,s1,t0,t1,t2,g[0],g[1],g[2],f[1],f[2]); _ff_set(f[1],w2);
	for ( i = 2 ; i < d-2 ; i++ ) {
		_ff_sum_6_mults (w2,s0,s1,t0,t1,t2,t3,g[i-2],g[i-1],g[i],g[i+1],f[i],f[i+1]);
		_ff_set(f[i],w1); _ff_set(w1,w2);
	}
	_ff_set(f[i],w1);
//ff_poly_print(f,d-2);
//puts("");	
	// (6d+17)M (5d+7)A (d+16) redc
}

int _ff_poly_bgcd (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g, int provisional);
int _ff_poly_bgcd_special (ff_t h[], ff_t f0[], int d_f, ff_t g0[], int d_g);

// both polys must be monic and have a non-trivial gcd (this is assumed and *not* verified)
// the function destroys both f and g
// if the gcd has degree > d_h, zero is returned, otherwise the degree of the gcd is returned
void ff_poly_bgcd (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g)
{
//printf ("clf(%d,%d)\n", d_f, d_g);
	if ( d_f < d_g )  _swap_polys(&f,&d_f,&g,&d_g);

	// send small cases to fast inline code, but check for failure
	if ( d_g < 9 && !pd_h ) {
		switch (d_g) {
		case 1: _ff_set(h[1],g[1]); _ff_set(h[0],g[0]); return;
		case 2: ff_poly_mgcd_linear_n_2 (h, f, d_f, g); break;
		case 3: ff_poly_mgcd_linear_n_3 (h, f, d_f, g); break;
		case 4: ff_poly_mgcd_linear_n_4 (h, f, d_f, g); break;
		case 5: ff_poly_mgcd_linear_n_5 (h, f, d_f, g); break;
		case 6: ff_poly_mgcd_linear_n_6 (h, f, d_f, g); break;
		case 7: ff_poly_mgcd_linear_n_7 (h, f, d_f, g); break;
		case 8: ff_poly_mgcd_linear_n_8 (h, f, d_f, g); break;
		}
		if ( ! _ff_zero(h[1]) ) return;
	}
	ff_poly_mod_small_inplace (f, d_f, g, d_g);
	d_f = ff_poly_degree(f,d_g-1);
	if ( d_f > 24 && _ff_poly_bgcd_special (h, f, d_f, g, d_g) ) { if ( pd_h ) *pd_h = 1; return; }
	_ff_poly_bgcd (h, pd_h, g, d_g, f, d_f, 0);
	return;
}

int _ff_poly_bgcd (ff_t h[], int *pd_h, ff_t f[], int d_f, ff_t g[], int d_g, int provisional)
{
	register ff_t t0, t1, t2, t3;
	register int i, k;

	for (;;) {
		if ( d_f < d_g )  _swap_polys(&f,&d_f,&g,&d_g);
		if ( d_g < 2 ) break;
		if ( d_f !=	d_g+1 ) { ff_poly_reduce_inplace (f, &d_f, g, d_g); continue; }
		_ff_neg(t0,g[d_g]);  _ff_sum_2_mults(t1,f[d_g],f[d_f],g[d_g-1],t0); 					// t1 is alpha
		_ff_mult(t2,t0,f[d_f]);  _ff_square(t3,t0);
		for ( i = d_g-1 ; i ; i-- ) _ff_sum_3_mults(f[i],f[i],t1,t2,g[i-1],g[i],t3);
		_ff_sum_2_mults(f[0],f[0],t1,g[0],t3);
		d_f = ff_poly_degree(f,d_g-1);
	}
	if ( d_g == 1 ) { _ff_set(h[1],g[1]); _ff_set(h[0],g[0]); if ( pd_h ) *pd_h = 1; return 1; }
	if ( provisional ) return 0;
	if ( !d_g ) { err_printf ("polys have trivial gcd in ff_poly_bgcd\n"); abort(); }
	if ( !pd_h || d_f > *pd_h ) {
		ff_poly_monic (f,&d_f,f,d_f);
		k = ff_poly_distinct_roots(g,f,d_f);
		if ( ! k ) { err_printf ("poly gcd has no linear factor in ff_poly_bgcd with p=%ld\n\n", _ff_p); ff_poly_print (f, d_f); abort(); }
		if ( (!pd_h && k > 1) || (pd_h && k > *pd_h) ) { err_printf ("gcd has %d distinct roots, exceeding limit %d in ff_poly_bgcd with p=%ld\n", d_f, *pd_h, _ff_p); ff_poly_print(f, d_f);  abort(); }
		ff_poly_from_roots_naive (h,g,k);	// this is stupid, caller is just going to factor it anyway, but we expect the degree to almost alwyas be 1
		if ( pd_h ) *pd_h = k;
		return k;
	}
	for ( i = 0 ; i <= d_f ; i++ ) _ff_set(h[i],f[i]);
	*pd_h = d_f;
	return d_f;
}

// g must be monic (not checked) and have degree one greater than f (checked)
// this function intentionally fails if anything wierd happens
int _ff_poly_bgcd_special (ff_t h[], ff_t *f, int d_f, ff_t *g, int d_g)
{
	ff_t f0[FF_POLY_MAX_DEGREE+1], g0[FF_POLY_MAX_DEGREE+1];

	if ( d_g != d_f+1 || d_f > FF_POLY_MAX_DEGREE || d_g > FF_POLY_MAX_DEGREE ) return 0;
	memcpy(f0,f,(d_f+1)*sizeof(ff_t));  memcpy(g0,g,(d_g+1)*sizeof(ff_t));
	f = f0; g = g0;
	ff_poly_reduce3m_inplace (g,f,d_f);
	d_g = ff_poly_degree(g,d_f-2);
//ff_poly_print(f,d_f);
//ff_poly_print(g,d_g);
	while ( d_g >= 24 ) {
		if ( d_f != d_g+2 ) return 0; 												// don't mess around with non-normal degree sequence case, let ff_poly_bgcd handle it
		ff_poly_reduce4_inplace (f, g, d_g);
//ff_poly_print(f,d_f);
//ff_poly_print(g,d_g);
		d_f = ff_poly_degree(f,d_g-2);
		_swap_polys(&f,&d_f,&g,&d_g);
	}
//puts ("end loop");
//ff_poly_print(f,d_f);
//ff_poly_print(g,d_g);
	if ( ! _ff_poly_bgcd (h, 0, f, d_f, g, d_g, 1) ) { info_printf ("bgcd_special failed\n"); return 0; }
//printf ("bgcd special succeeded h[1]=%ld, h[0]=%ld\n", _ff_get_ui(h[1]), _ff_get_ui(h[0]));
	return 1;
}
