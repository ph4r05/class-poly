#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ff_poly.h"
#include "ff_poly/ffpolysmall.h"
#include "phi_poly.h"
#include "phi_eval.h"
#include "phi_gcd.h"
#include "cstd.h"

/*
    Copyright 2012 Andrew V. Sutherland

    This file is part of classpoly.

    classpoly is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    classpoly is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with classpoly.  If not, see <http://www.gnu.org/licenses/>.
*/

int  phi_surface_gcd_cycle_7_11 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_7, ff_t *phi_11)
{
	ff_t f[13], g[9], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;

	i1=j2= 0; i2=j1=e-1;
	if ( r != r2 ) { i1=j1+1;  i2 = n-1; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	switch ( phi_sparse_factor(phi_7) ) {
	case 2:
		do {
			phi11_s2_eval_ff (f, phi_11, r2[i1]);
			phi7_s2_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_12_7 (h1, f, g);
			phi11_s2_eval_ff (f, phi_11, r2[i2]);
			phi7_s2_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_12_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 3:
		do {
			phi11_s3_eval_ff (f, phi_11, r2[i1]);
			phi7_s3_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_12_7 (h1, f, g);
			phi11_s3_eval_ff (f, phi_11, r2[i2]);
			phi7_s3_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_12_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 4:
		do {
			phi11_s4_eval_ff (f, phi_11, r2[i1]);
			phi7_s4_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_12_7 (h1, f, g);
			phi11_s4_eval_ff (f, phi_11, r2[i2]);
			phi7_s4_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_12_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 6:
		do {
			phi11_s6_eval_ff (f, phi_11, r2[i1]);
			phi7_s6_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_12_7 (h1, f, g);
			phi11_s6_eval_ff (f, phi_11, r2[i2]);
			phi7_s6_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_12_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 8:
		do {
			phi11_s8_eval_ff (f, phi_11, r2[i1]);
			phi7_s8_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_12_7 (h1, f, g);
			phi11_s8_eval_ff (f, phi_11, r2[i2]);
			phi7_s8_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_12_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 12:
		do {
			phi11_s12_eval_ff (f, phi_11, r2[i1]);
			phi7_s12_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_12_7 (h1, f, g);
			phi11_s12_eval_ff (f, phi_11, r2[i2]);
			phi7_s12_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_12_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 24:
		do {
			phi11_s24_eval_ff (f, phi_11, r2[i1]);
			phi7_s24_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_12_7 (h1, f, g);
			phi11_s24_eval_ff (f, phi_11, r2[i2]);
			phi7_s24_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_12_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	default:
		do {
			phi11_eval_ff (f, phi_11, r2[i1]);
			phi7_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_12_7 (h1, f, g);
			phi11_eval_ff (f, phi_11, r2[i2]);
			phi7_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_12_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	}
	return 1;
}

int phi_surface_gcd_cycle_7_13 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_7, ff_t *phi_13)
{
	ff_t f[15], g[9], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;

	i1=j2= 0; i2=j1=e-1;
	if ( r != r2 ) { i1=j1+1;  i2 = n-1; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	switch ( phi_sparse_factor(phi_7) ) {
	case 2:
		do {
			phi13_s2_eval_ff (f, phi_13, r2[i1]);
			phi7_s2_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_14_7 (h1, f, g);
			phi13_s2_eval_ff (f, phi_13, r2[i2]);
			phi7_s2_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_14_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 3:
		do {
			phi13_s3_eval_ff (f, phi_13, r2[i1]);
			phi7_s3_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_14_7 (h1, f, g);
			phi13_s3_eval_ff (f, phi_13, r2[i2]);
			phi7_s3_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_14_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 4:
		do {
			phi13_s4_eval_ff (f, phi_13, r2[i1]);
			phi7_s4_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_14_7 (h1, f, g);
			phi13_s4_eval_ff (f, phi_13, r2[i2]);
			phi7_s4_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_14_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 6:
		do {
			phi13_s6_eval_ff (f, phi_13, r2[i1]);
			phi7_s6_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_14_7 (h1, f, g);
			phi13_s6_eval_ff (f, phi_13, r2[i2]);
			phi7_s6_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_14_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 8:
		do {
			phi13_s8_eval_ff (f, phi_13, r2[i1]);
			phi7_s8_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_14_7 (h1, f, g);
			phi13_s8_eval_ff (f, phi_13, r2[i2]);
			phi7_s8_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_14_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 12:
		do {
			phi13_s12_eval_ff (f, phi_13, r2[i1]);
			phi7_s12_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_14_7 (h1, f, g);
			phi13_s12_eval_ff (f, phi_13, r2[i2]);
			phi7_s12_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_14_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 24:
		do {
			phi13_s24_eval_ff (f, phi_13, r2[i1]);
			phi7_s24_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_14_7 (h1, f, g);
			phi13_s24_eval_ff (f, phi_13, r2[i2]);
			phi7_s24_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_14_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	default:
		do {
			phi13_eval_ff (f, phi_13, r2[i1]);
			phi7_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_14_7 (h1, f, g);
			phi13_eval_ff (f, phi_13, r2[i2]);
			phi7_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_14_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	}
	return 1;
}

int phi_surface_gcd_cycle_7_17 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_7, ff_t *phi_17)
{
	ff_t f[19], g[9], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;

	i1=j2= 0; i2=j1=e-1;
	if ( r != r2 ) { i1=j1+1;  i2 = n-1; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	switch ( phi_sparse_factor(phi_7) ) {
	case 2:
		do {
			phi17_s2_eval_ff (f, phi_17, r2[i1]);
			phi7_s2_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_18_7 (h1, f, g);
			phi17_s2_eval_ff (f, phi_17, r2[i2]);
			phi7_s2_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_18_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 3:
		do {
			phi17_s3_eval_ff (f, phi_17, r2[i1]);
			phi7_s3_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_18_7 (h1, f, g);
			phi17_s3_eval_ff (f, phi_17, r2[i2]);
			phi7_s3_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_18_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 4:
		do {
			phi17_s4_eval_ff (f, phi_17, r2[i1]);
			phi7_s4_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_18_7 (h1, f, g);
			phi17_s4_eval_ff (f, phi_17, r2[i2]);
			phi7_s4_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_18_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 6:
		do {
			phi17_s6_eval_ff (f, phi_17, r2[i1]);
			phi7_s6_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_18_7 (h1, f, g);
			phi17_s6_eval_ff (f, phi_17, r2[i2]);
			phi7_s6_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_18_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 8:
		do {
			phi17_s8_eval_ff (f, phi_17, r2[i1]);
			phi7_s8_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_18_7 (h1, f, g);
			phi17_s8_eval_ff (f, phi_17, r2[i2]);
			phi7_s8_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_18_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 12:
		do {
			phi17_s12_eval_ff (f, phi_17, r2[i1]);
			phi7_s12_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_18_7 (h1, f, g);
			phi17_s12_eval_ff (f, phi_17, r2[i2]);
			phi7_s12_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_18_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 24:
		do {
			phi17_s24_eval_ff (f, phi_17, r2[i1]);
			phi7_s24_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_18_7 (h1, f, g);
			phi17_s24_eval_ff (f, phi_17, r2[i2]);
			phi7_s24_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_18_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	default:
		do {
			phi17_eval_ff (f, phi_17, r2[i1]);
			phi7_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_18_7 (h1, f, g);
			phi17_eval_ff (f, phi_17, r2[i2]);
			phi7_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_18_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	}
	return 1;
}

int phi_surface_gcd_cycle_7_19 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_7, ff_t *phi_19)
{
	ff_t f[21], g[9], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;

	i1=j2= 0; i2=j1=e-1;
	if ( r != r2 ) { i1=j1+1;  i2 = n-1; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	switch ( phi_sparse_factor(phi_7) ) {
	case 2:
		do {
			phi19_s2_eval_ff (f, phi_19, r2[i1]);
			phi7_s2_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_20_7 (h1, f, g);
			phi19_s2_eval_ff (f, phi_19, r2[i2]);
			phi7_s2_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_20_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 3:
		do {
			phi19_s3_eval_ff (f, phi_19, r2[i1]);
			phi7_s3_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_20_7 (h1, f, g);
			phi19_s3_eval_ff (f, phi_19, r2[i2]);
			phi7_s3_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_20_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 4:
		do {
			phi19_s4_eval_ff (f, phi_19, r2[i1]);
			phi7_s4_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_20_7 (h1, f, g);
			phi19_s4_eval_ff (f, phi_19, r2[i2]);
			phi7_s4_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_20_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 6:
		do {
			phi19_s6_eval_ff (f, phi_19, r2[i1]);
			phi7_s6_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_20_7 (h1, f, g);
			phi19_s6_eval_ff (f, phi_19, r2[i2]);
			phi7_s6_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_20_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 8:
		do {
			phi19_s8_eval_ff (f, phi_19, r2[i1]);
			phi7_s8_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_20_7 (h1, f, g);
			phi19_s8_eval_ff (f, phi_19, r2[i2]);
			phi7_s8_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_20_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 12:
		do {
			phi19_s12_eval_ff (f, phi_19, r2[i1]);
			phi7_s12_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_20_7 (h1, f, g);
			phi19_s12_eval_ff (f, phi_19, r2[i2]);
			phi7_s12_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_20_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 24:
		do {
			phi19_s24_eval_ff (f, phi_19, r2[i1]);
			phi7_s24_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_20_7 (h1, f, g);
			phi19_s24_eval_ff (f, phi_19, r2[i2]);
			phi7_s24_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_20_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	default:
		do {
			phi19_eval_ff (f, phi_19, r2[i1]);
			phi7_eval_ff (g, phi_7, r[j1]);
			ff_poly_remove_root (g, g, 8, r+j1-1);
			ff_poly_mgcd_linear_20_7 (h1, f, g);
			phi19_eval_ff (f, phi_19, r2[i2]);
			phi7_eval_ff (g, phi_7, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 8, r+k);
			ff_poly_mgcd_linear_20_7 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	}
	return 1;
}

#define BATCH_INVERTS		128

// computes s^2*g/(x-r/s) assuming g(r/s)=0 for a cubic g (not assumed monic)
static inline void ff_poly_remove_qroot_3 (ff_t g[4], ff_t r[1], ff_t s[1])
{
	register ff_t t1,t2;

	_ff_square(t2,s[0]);
	_ff_sum_2_mults(t1,g[2],g[3],r[0],s[0]);
	_ff_sum_2_mults(g[0],t1,t2,g[1],r[0]);
	_ff_mult(g[1],s[0],t1);
	_ff_mult(g[2],t2,g[3]);
	// 7M + 2A (5 redc)
}

// computes b^3*d^2*phi2(X,a/b)/(X-c/d) assuming (X-c/d) is a root of phi2(X,a/b)
// combines call to phi2_s3_qeval and ff_poly_remove_qroot_3.  In principle this should be noticeably faster (uses 11M instead of 15M), but with complier optimization it  doesn't appear to be much faster).
// currently unused
static inline void phi2_s3_qeval_rem (ff_t g[3], ff_t phi2[], ff_t a, ff_t b, ff_t c, ff_t d)
{
	register ff_t bc, bd, bd2, t1, t2;
	
	_ff_mult (bc, b, c);  _ff_mult (bd, b, d);  _ff_square (bd2, bd);
	_ff_mult (g[2],b,bd2);																// g[2]=b^3d^2
	_ff_square (t1, a); ff_negate(t1);  _ff_sum_2_mults (t2,b,d,t1,bc);								// t2 = b^2c-a^2d
	_ff_mult (g[1],bd,t2);																// g[1]=b^3cd-a^2bd^2
	_ff_mult (t1,a,phi2[4]);
	_ff_sum_2_mults (g[0],bc,bd2,t1,t2);													// g[0]=phi[4]ab^2d^2+bc(b^2c-a^2d)
	// 11M+2A (9 redc)
}


// computes s^3*g/(x-r/s) assuming g(r/s)=0 for a quartic g (not assumed monic)
static inline void ff_poly_remove_qroot_4 (ff_t g[5], ff_t r[1], ff_t s[1])
{
	register ff_t t1,t2,s2,s3;

	_ff_square(s2,s[0]);  _ff_mult(s3,s2,s[0]);
	_ff_sum_2_mults(t2,g[3],g[4],r[0],s[0]);
	_ff_sum_2_mults(t1,g[2],t2,r[0],s2);
	_ff_sum_2_mults(g[0],g[1],t1,r[0],s3);
	_ff_mult(g[1],t1,s[0]);
	_ff_mult(g[2],t2,s2);
	_ff_mult(g[3],g[4],s3);
	// 11M + 3A (8 redc)
}

// computes s^3*g/(x-r/s) assuming g(r/s)=0 for a quartic g (not assumed monic)
static inline void ff_poly_remove_qroot_6 (ff_t g[7], ff_t r[1], ff_t s[1])
{
	register ff_t t1,t2,t3,t4,s2,s3,s4,s5;

	_ff_square(s2,s[0]);  _ff_mult(s3,s2,s[0]); _ff_square(s4,s2); _ff_mult(s5,s2,s3);
	_ff_sum_2_mults(t4,g[5],g[6],r[0],s[0]);
	_ff_sum_2_mults(t3,g[4],t4,r[0],s2);
	_ff_sum_2_mults(t2,g[3],t3,r[0],s3);
	_ff_sum_2_mults(t1,g[2],t2,r[0],s4);
	_ff_sum_2_mults(g[0],g[1],t1,r[0],s5);
	_ff_mult(g[1],t1,s[0]);
	_ff_mult(g[2],t2,s2);
	_ff_mult(g[3],t3,s3);
	_ff_mult(g[4],t4,s4);
	_ff_mult(g[5],g[6],s5);
	// 19M + 5A (14 redc)
}

// computes s^3*g/(x-r/s) assuming g(r/s)=0 for a quartic g (not assumed monic)
static inline void ff_poly_remove_qroot_8 (ff_t g[9], ff_t r[1], ff_t s[1])
{
	register ff_t t1,t2,t3,t4,t5,t6,s2,s3,s4,s5,s6,s7;

	_ff_square(s2,s[0]);  _ff_mult(s3,s2,s[0]); _ff_square(s2,s2); _ff_mult(s5,s2,s3); _ff_square(s6,s3); _ff_mult(s7,s3,s4);
	_ff_sum_2_mults(t6,g[7],g[8],r[0],s[0]);
	_ff_sum_2_mults(t5,g[6],t6,r[0],s2);
	_ff_sum_2_mults(t4,g[5],t5,r[0],s3);
	_ff_sum_2_mults(t3,g[4],t4,r[0],s4);
	_ff_sum_2_mults(t2,g[3],t3,r[0],s5);
	_ff_sum_2_mults(t1,g[2],t2,r[0],s6);
	_ff_sum_2_mults(g[0],g[1],t1,r[0],s7);
	_ff_mult(g[1],t1,s[0]);
	_ff_mult(g[2],t2,s2);
	_ff_mult(g[3],t3,s3);
	_ff_mult(g[4],t4,s4);
	_ff_mult(g[5],t5,s5);
	_ff_mult(g[6],t6,s6);
	_ff_mult(g[7],g[8],s7);
	// 27M + 7A (20 redc)
}

/*
	Enumerates a surface p1-cycle using gcds with Phi_p2, where alpha_p2=alpha_p1^e and |alpha_p1|=n
	(Here alpha_a denotes the class represented by the pos def reduced primeform <a,b,c>).
	Assumes n > e > 1 and that roots[0],...,roots[e-1] are already present.
*/
int phi_surface_qgcd_cycle_2_3 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_2, ff_t *phi_3)
{
	ff_t f[5], g[4], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m, m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	while ( j < n ) {
		m = ( m0 > n-j ? n-j : m0 );
		j0 = j-2;
		for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
			phi3_eval_ff (f, phi_3, r2[i]);
			phi2_qeval_ff (g, phi_2, r[j-1], d[k-1]);
			ff_poly_remove_qroot_3 (g, r+j-2, d+k-2);
			ff_poly_gcd_linear_4_2 (h, f, g);		// this is about the same (maybe slightly faster) as using ff_poly_gcd_linear_4_3
			_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
		}
		if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
		for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
	}
	return 1;
}

int phi_surface_qgcd_cycle_2_5 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_2, ff_t *phi_5)
{
	ff_t f[7], g[4], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m, m0;

	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case 
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	// duplicate loops to optimize for sparser gamma_2 case
	if ( phi_sparse_factor(phi_2)==3 ) {
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi5_s3_eval_ff (f, phi_5, r2[i]);
//				phi2_s3_qeval_rem (g, phi_2, r[j-1], d[k-1], r[j-2], d[k-2]);		doesn't appear to be noticeably faster
				phi2_s3_qeval_ff (g, phi_2, r[j-1], d[k-1]);
				ff_poly_remove_qroot_3 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_6_2 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	} else {
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi5_eval_ff (f, phi_5, r2[i]);
				phi2_qeval_ff (g, phi_2, r[j-1], d[k-1]);
				ff_poly_remove_qroot_3 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_6_2 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}

int phi_surface_qgcd_cycle_2_7 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_2, ff_t *phi_7)
{
	ff_t f[9], g[4], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	// duplicate loops to optimize for sparser gamma_2 case
	if ( phi_sparse_factor(phi_2)==3 ) {
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_s3_eval_ff (f, phi_7, r2[i]);
				phi2_s3_qeval_ff (g, phi_2, r[j-1], d[k-1]);
				ff_poly_remove_qroot_3 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_8_2 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	} else {
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_eval_ff (f, phi_7, r2[i]);
				phi2_qeval_ff (g, phi_2, r[j-1], d[k-1]);
				ff_poly_remove_qroot_3 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_8_2 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}

int phi_surface_qgcd_cycle_2_11 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_2, ff_t *phi_11)
{
	ff_t f[13], g[4], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	// duplicate loops to optimize for sparser gamma_2 case
	if ( phi_sparse_factor(phi_2)==3 ) {
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_s3_eval_ff (f, phi_11, r2[i]);
				phi2_s3_qeval_ff (g, phi_2, r[j-1], d[k-1]);
				ff_poly_remove_qroot_3 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_12_2 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	} else {
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_eval_ff (f, phi_11, r2[i]);
				phi2_qeval_ff (g, phi_2, r[j-1], d[k-1]);
				ff_poly_remove_qroot_3 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_12_2 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}

int phi_surface_qgcd_cycle_2_13 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_2, ff_t *phi_13)
{
	ff_t f[15], g[4], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	// duplicate loops to optimize for sparser gamma_2 case
	if ( phi_3sparse(phi_2) ) {
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_s3_eval_ff (f, phi_13, r2[i]);
				phi2_s3_qeval_ff (g, phi_2, r[j-1], d[k-1]);
				ff_poly_remove_qroot_3 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_14_2 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	} else {
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_eval_ff (f, phi_13, r2[i]);
				phi2_qeval_ff (g, phi_2, r[j-1], d[k-1]);
				ff_poly_remove_qroot_3 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_14_2 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}

int phi_surface_qgcd_cycle_2_p2 (ff_t r[], ff_t *r2, int n, int p2, int e, ff_t *phi_2, ff_t *phi_p2)
{
	ff_t f[PHI_MAX_M+2], g[4], d[BATCH_INVERTS+2], h[2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	while ( j < n ) {
		m = ( m0 > n-j ? n-j : m0 );
		j0 = j-2;
		for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
			phi_eval_ff (f, phi_p2, p2, r2[i]);
			phi2_qeval_ff (g, phi_2, r[j-1], d[k-1]);
			ff_poly_remove_qroot_3 (g, r+j-2, d+k-2);
			ff_poly_gcd_linear_n_2 (h, f, p2+1, g);
			_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
		}
		if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
		for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
	}
	return 1;
}

int phi_surface_qgcd_cycle_3_5 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_3, ff_t *phi_5)
{
	ff_t f[7], g[5], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	switch ( phi_sparse_factor(phi_3) ) {
	case 2:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi5_s2_eval_ff (f, phi_5, r2[i]);
				phi3_s2_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_gcd_linear_6_4 (h, f, g);				// this is faster than removing a root and calling gcd_linear_6_3
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 4:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi5_s4_eval_ff (f, phi_5, r2[i]);
				phi3_s4_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_gcd_linear_6_4 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 8:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi5_s8_eval_ff (f, phi_5, r2[i]);
				phi3_s8_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_gcd_linear_6_4 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	default:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi5_eval_ff (f, phi_5, r2[i]);
				phi3_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_gcd_linear_6_4 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}


int phi_surface_qgcd_cycle_3_7 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_3, ff_t *phi_7)
{
	ff_t f[9], g[5], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	switch ( phi_sparse_factor(phi_3) ) {
	case 2:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_s2_eval_ff (f, phi_7, r2[i]);
				phi3_s2_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_8_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 4:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_s4_eval_ff (f, phi_7, r2[i]);
				phi3_s4_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_8_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 8:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_s8_eval_ff (f, phi_7, r2[i]);
				phi3_s8_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_8_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	default:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_eval_ff (f, phi_7, r2[i]);
				phi3_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_8_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}

int phi_surface_qgcd_cycle_3_11 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_3, ff_t *phi_11)
{
	ff_t f[13], g[5], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	switch ( phi_sparse_factor(phi_3) ) {
	case 2:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_s2_eval_ff (f, phi_11, r2[i]);
				phi3_s2_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_12_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 4:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_s4_eval_ff (f, phi_11, r2[i]);
				phi3_s4_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_12_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 8:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_s8_eval_ff (f, phi_11, r2[i]);
				phi3_s8_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_12_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	default:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_eval_ff (f, phi_11, r2[i]);
				phi3_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_12_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}


int phi_surface_qgcd_cycle_3_13 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_3, ff_t *phi_13)
{
	ff_t f[15], g[5], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	switch ( phi_sparse_factor(phi_3) ) {
	case 2:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_s2_eval_ff (f, phi_13, r2[i]);
				phi3_s2_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_14_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 4:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_s4_eval_ff (f, phi_13, r2[i]);
				phi3_s4_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_14_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 8:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_s8_eval_ff (f, phi_13, r2[i]);
				phi3_s8_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_14_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	default:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_eval_ff (f, phi_13, r2[i]);
				phi3_qeval_ff (g, phi_3, r[j-1], d[k-1]);
				ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_14_3 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}

int phi_surface_qgcd_cycle_3_p2 (ff_t r[], ff_t *r2, int n, int p2, int e, ff_t *phi_3, ff_t *phi_p2)
{
	ff_t f[PHI_MAX_M+2], g[5], d[BATCH_INVERTS+2], h[2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	while ( j < n ) {
		m = ( m0 > n-j ? n-j : m0 );
		j0 = j-2;
		for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
			phi_eval_ff (f, phi_p2, p2, r2[i]);
			phi3_qeval_ff (g, phi_3, r[j-1], d[k-1]);
			ff_poly_remove_qroot_4 (g, r+j-2, d+k-2);
			ff_poly_gcd_linear_n_3 (h, f, p2+1, g);
			_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
		}
		if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
		for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
	}
	return 1;
}


static inline void ff_poly_s24_gcd_linear_8_6 (ff_t h[2], ff_t f[8], ff_t g[7])	// assumes f = x^8+f7x^7+f4x^4+f1x+f0,  g = g6x^6+g5x^5+g1x+g0
{
	register ff_t t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, g6, g62, ng0, ng1, ng5, w, w1;

	_ff_set(g6,g[6]); _ff_square(g62,g6); _ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng5,g[5]);
	_ff_mult(t1,f[7],g6); _ff_addto(t1,ng5);
	_ff_mult(t0,t1,ng5);
	_ff_mult(t5,t0,ng5);
	ff_mult(w,g6,g62);  _ff_mult(t4,f[4],w);
	_ff_mult(t3,ng1,g62);
	_ff_mult(w1,g6,t1);
	_ff_sum_2_mults(t2,g62,w1,ng1,ng0);
	_ff_sum_3_mults(t1,f[1],w1,t0,ng1,ng0,w);
	_ff_sum_2_mults(t0,f[0],t0,ng0,w);
	// 15M
	_ff_sum_2_mults(s0,ng5,g6,t4,t5);
	_ff_mult(w,t5,g6); _ff_neg(w1,w);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_2_mults(s3,t2,t3,s0,w1);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_square(w,t5);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w);
	_ff_sum_2_mults(s0,g[0],t0,s0,w);
	// 30M+16A (19 redc)
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 66M+37A (37 redc)
}

static inline void ff_poly_s12_gcd_linear_8_6 (ff_t h[2], ff_t f[8], ff_t g[7])	// assumes f = x^8+f7x^7+f6x^6+f4x^4+f2x^2+f1x+f0,  g = g6x^6+g5x^5+g3x^3+g1x+g0 
{
	register ff_t t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, g6, g62, g63, ng0, ng1, ng3, ng5, w0, w1;

	_ff_set(g6,g[6]); _ff_square(g62,g6); _ff_mult(g63,g62,g6); _ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng3,g[3]); _ff_neg(ng5,g[5]);
	_ff_mult(t1,f[7],g6); _ff_addto(t1,ng5);
	_ff_sum_2_mults(t0,f[6],t1,ng5,g62);
	_ff_sum_2_mults(t5,t0,g62,ng3,ng5);
	ff_mult(t1,t1,g6);
	_ff_sum_2_mults(t4,f[4],t1,ng3,g63);
	_ff_sum_2_mults(t3,t0,g62,ng1,ng3);
	_ff_sum_3_mults(t2,f[2],t1,g62,ng0,ng1,g63);
	_ff_sum_3_mults(t1,f[1],t0,t1,ng0,ng1,g63);
	_ff_sum_2_mults(t0,f[0],t0,ng0,g63);
	_ff_sum_2_mults(s0,ng5,g6,t4,t5);
	_ff_mult(w0,t5,g6); _ff_neg(w1,w0); _ff_square(w0,t5);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_3_mults(s3,g[3],t2,t3,s0,w1,w0);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w0);
	_ff_sum_2_mults(s0,g[0],t0,s0,w0);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 72M+48A (37 redc)
}


static inline void ff_poly_s24_gcd_linear_12_6 (ff_t h[2], ff_t f[12], ff_t g[7])	// assumes f = x^12+f11x^11+f9x^9+f7x^7+f5x^5+f3x^3+f1x+f0,  g = g6x^6+g5x^5+g1x+g0
{
	register ff_t t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, w, w0, w1, w5, ng0, ng1, g6, g62, g63, g64, g65, g67;
	
	_ff_set(g6,g[6]); _ff_neg(w5,g[5]); _ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]);
	_ff_mult(t5,f[11],g6); _ff_addto(t5,w5);
	_ff_mult(t4,t5,w5); _ff_square(g62,g6); _ff_mult(g63,g62,g6); _ff_square(g64,g62); _ff_mult(g65,g62,g63); _ff_mult(g67,g65,g62);
	_ff_sum_2_mults(t3,f[9],t4,w5,g63);
	_ff_mult(t2,t3,w5);
	_ff_mult(w1,g64,ng1);  _ff_mult(w0,g65,ng0);
	_ff_sum_2_mults(t1,f[7],t2,w5,g65); _ff_addto(t1,w1);
	_ff_sum_2_mults(t0,t1,t5,w1,w5); _ff_addto(t0,w0);
	_ff_sum_4_mults(t5,f[5],t0,t4,t5,w0,w1,w5,g67);
	_ff_mult(w0,g64,ng0); _ff_mult(w1,g63,ng1);
	_ff_sum_2_mults(t4,t3,t4,w0,w1);
	_ff_mult(w0,g63,ng0); _ff_mult(w1,g62,ng1);
	_ff_sum_3_mults(t3,f[3],t2,t3,w0,w1,g67);
	_ff_mult(w0,g62,ng0); _ff_mult(w1,g6,ng1);
	_ff_sum_2_mults(t2,t1,t2,w0,w1);
	_ff_mult(w0,g6,ng0);
	_ff_sum_3_mults(t1,f[1],t0,t1,w0,ng1,g67);
	_ff_sum_2_mults(t0,f[0],t0,ng0,g67);
	_ff_sum_2_mults(s0,w5,g6,t4,t5);
	_ff_mult(w,t5,g6); _ff_neg(w1,w);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_2_mults(s3,t2,t3,s0,w1);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_square(w,t5);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w);
	_ff_sum_2_mults(s0,g[0],t0,s0,w);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 90M
}

static inline void ff_poly_s12_gcd_linear_12_6 (ff_t h[2], ff_t f[12], ff_t g[7])	// assumes g = g6x^6+g5x^5+g3x^3+g1x+g0
{
	register ff_t ng0, ng1, ng3, ng5, g6, g62, g63, g64, g65, g66, g67, t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, s5, w0, w1;
	
	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng3,g[3]); _ff_neg(ng5,g[5]); _ff_set(g6,g[6]); 
	_ff_square(g62,g6); _ff_mult(g63,g62,g6); _ff_square(g64,g62); _ff_mult(g65,g64,g6); _ff_square(g66,g63); _ff_mult(g67,g66,g6);
	ff_mult(s2,g62,ng3); _ff_mult(s4,g64,ng1); _ff_mult(s5,g65,ng0);
	_ff_mult(t5,f[11],g6); _ff_addto(t5,ng5);
	_ff_sum_2_mults(t4,f[10],t5,ng5,g62);
	_ff_sum_2_mults(t3,f[9],t4,ng5,g63); _ff_addto(t3,s2);
	_ff_sum_3_mults(t2,f[8],t3,t5,s2,ng5,g64);
	_ff_sum_3_mults(t1,f[7],t2,t4,s2,ng5,g65); _ff_addto(t1,s4);
	_ff_sum_4_mults(t0,f[6],t1,t3,t5,s4,s2,ng5,g66); _ff_addto(t0,s5);
	_ff_mult(s1,g6,t1); _ff_mult(s2,g62,t2); _ff_mult(s3,g63,t3); _ff_mult(s4,g64,t4); 
	_ff_sum_5_mults(t5,f[5],t0,s2,s4,s5,t5,ng1,ng3,ng5,g67);
	_ff_sum_4_mults(t4,f[4],s1,s3,s4,ng0,ng1,ng3,g67);
	_ff_sum_4_mults(t3,f[3],t0,s2,s3,ng0,ng1,ng3,g67);
	_ff_sum_3_mults(t2,f[2],s1,s2,ng0,ng1,g67);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,g67);
	_ff_sum_2_mults(t0,f[0],t0,ng0,g67);
	
	_ff_sum_2_mults(s0,ng5,g6,t4,t5);
	_ff_mult(w0,t5,g6); _ff_neg(w1,w0); _ff_square(w0,t5);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_3_mults(s3,g[3],t2,t3,s0,w1,w0);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w0);
	_ff_sum_2_mults(s0,g[0],t0,s0,w0);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 101M
}

static inline void ff_poly_s24_gcd_linear_14_6 (ff_t h[2], ff_t f[14], ff_t g[7])	// assumes f = x^14+f13x^13+f12x^12+f10x^10+f8x^8+f6x^6+f4x^4+f2x^2+f1x+f0,  g = g6x^6+g5x^5+g1x+g0
{
	register ff_t t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, w, w0, w1, w5, ng0, ng1, g6, g62, g63, g64, g65, g69;
	
	_ff_set(g6,g[6]); _ff_neg(w5,g[5]); _ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]);
	_ff_square(g62,g6); _ff_mult(g63,g62,g6); _ff_square(g64,g62); _ff_mult(g65,g62,g63); _ff_mult(g69,g65,g64);
	_ff_mult(t1,f[13],g6); _ff_addto(t1,w5);
	_ff_sum_2_mults(t0,f[12],t1,w5,g62);
	_ff_mult(t5,t0,w5);
	_ff_sum_2_mults(t4,f[10],t5,w5,g64);
	_ff_mult(w1,g64,ng1);  _ff_mult(w0,g65,ng0);
	_ff_mult(t3,t4,w5); _ff_addto(t3,w1);
	_ff_square(w,g63);
	_ff_sum_3_mults(t2,f[8],t3,t1,w1,w5,w); _ff_addto(t2,w0);
	_ff_sum_3_mults(t1,t2,t0,t1,w0,w1,w5);
	_ff_square(w,g64);
	_ff_sum_4_mults(t0,f[6],t1,t5,t0,w0,w1,w5,w);
	_ff_sum_3_mults(t5,t0,t4,t5,w0,w1,w5);
	_ff_mult(w0,g64,ng0); _ff_mult(w1,g63,ng1);
	_ff_sum_3_mults(t4,f[4],t3,t4,w0,w1,g69);
	_ff_mult(w0,g63,ng0); _ff_mult(w1,g62,ng1);
	_ff_sum_2_mults(t3,t2,t3,w0,w1);
	_ff_mult(w0,g62,ng0); _ff_mult(w1,g6,ng1);
	_ff_sum_3_mults(t2,f[2],t1,t2,w0,w1,g69);
	_ff_mult(w0,g6,ng0);
	_ff_sum_3_mults(t1,f[1],t0,t1,w0,ng1,g69);
	_ff_sum_2_mults(t0,f[0],t0,ng0,g69);
	_ff_sum_2_mults(s0,w5,g6,t4,t5);
	_ff_mult(w,t5,g6); _ff_neg(w1,w);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_2_mults(s3,t2,t3,s0,w1);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_square(w,t5);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w);
	_ff_sum_2_mults(s0,g[0],t0,s0,w);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 101M
}

static inline void ff_poly_s12_gcd_linear_14_6 (ff_t h[2], ff_t f[12], ff_t g[7])	// assumes g = g6x^6+g5x^5+g3x^3+g1x+g0
{
	register ff_t ng0, ng1, ng3, ng5, g6, g62, g63, g64, g65, t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, s5, w0, w1;

	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng3,g[3]); _ff_neg(ng5,g[5]); _ff_set(g6,g[6]); 
	_ff_square(g62,g6); _ff_mult(g63,g62,g6); _ff_square(g64,g62); _ff_mult(g65,g64,g6);
	ff_mult(s2,g62,ng3); _ff_mult(s4,g64,ng1); _ff_mult(s5,g65,ng0);
	_ff_mult(t1,f[13],g6); _ff_addto(t1,ng5);
	_ff_sum_2_mults(t0,f[12],t1,ng5,g62);
	_ff_sum_2_mults(t5,f[11],t0,ng5,g63); _ff_addto(t5,s2);
	_ff_sum_3_mults(t4,f[10],t5,t1,s2,ng5,g64);
	_ff_sum_3_mults(t3,f[9],t4,t0,s2,ng5,g65); _ff_addto(t3,s4);
	_ff_mult(w0,g65,g6);
	_ff_sum_4_mults(t2,f[8],t3,t5,t1,s4,s2,ng5,w0); _ff_addto(t2,s5);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t1,f[7],t2,t4,t0,t1,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t0,f[6],t1,t3,t5,t0,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_mult(s1,g6,t1); _ff_mult(s2,g62,t2); _ff_mult(s3,g63,t3); _ff_mult(s4,g64,t4); 
	_ff_sum_5_mults(t5,f[5],t0,s2,s4,s5,t5,ng1,ng3,ng5,w0);
	_ff_sum_4_mults(t4,f[4],s1,s3,s4,ng0,ng1,ng3,w0);
	_ff_sum_4_mults(t3,f[3],t0,s2,s3,ng0,ng1,ng3,w0);
	_ff_sum_3_mults(t2,f[2],s1,s2,ng0,ng1,w0);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,w0);
	_ff_sum_2_mults(t0,f[0],t0,ng0,w0);
	
	_ff_sum_2_mults(s0,ng5,g6,t4,t5);
	_ff_mult(w0,t5,g6); _ff_neg(w1,w0); _ff_square(w0,t5);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_3_mults(s3,g[3],t2,t3,s0,w1,w0);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w0);
	_ff_sum_2_mults(s0,g[0],t0,s0,w0);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 113M
}

static inline void ff_poly_s24_gcd_linear_18_6 (ff_t h[2], ff_t f[18], ff_t g[7])	// assumes f is monic with f4, f7, f11, f14 coeefficients zero (true for Phi_17^f(X,j)
															//  g = g6x^6+g5x^5+g1x+g0 (this will be true for Phi_13^f(X,j) and s^6*Phi_5^f(X,r/s))
{
	register ff_t t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, w, w0, w1, w5, ng0, ng1, g6, g62, g63, g64, g65;
	
	_ff_set(g6,g[6]); _ff_neg(w5,g[5]); _ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]);
	_ff_mult(t5,f[17],g6); _ff_addto(t5,w5);
	_ff_mult(t4,t5,w5); _ff_square(g62,g6); _ff_mult(g63,g62,g6); _ff_square(g64,g62); _ff_mult(g65,g62,g63);
	_ff_sum_2_mults(t4,f[16],t5,w5,g62);
	_ff_sum_2_mults(t3,f[15],t4,w5,g63);
	_ff_mult(t2,t3,w5);
	_ff_mult(w1,g64,ng1);  _ff_mult(w0,g65,ng0);
	_ff_sum_2_mults(t1,f[13],t2,w5,g65); _ff_addto(t1,w1);
	_ff_mult(w,g65,g6);
	_ff_sum_3_mults(t0,f[12],t1,t5,w1,w5,w); _ff_addto(t0,w0);
	_ff_sum_3_mults(t5,t0,t4,t5,w0,w1,w5);
	_ff_mult(w,w,g62);
	_ff_sum_4_mults(t4,f[10],t5,t3,t4,w0,w1,w5,w);
	_ff_mult(w,w,g6);
	_ff_sum_4_mults(t3,f[9],t4,t2,t3,w0,w1,w5,w);
	_ff_mult(w,w,g6);
	_ff_sum_4_mults(t2,f[8],t3,t1,t2,w0,w1,w5,w);
	_ff_sum_3_mults(t1,t2,t0,t1,w0,w1,w5);
	_ff_mult(w,w,g62);
	_ff_sum_4_mults(t0,f[6],t1,t5,t0,w0,w1,w5,w);
	_ff_mult(w,w,g6);
	_ff_sum_4_mults(t5,f[5],t0,t4,t5,w0,w1,w5,w);
	_ff_mult(w0,g64,ng0); _ff_mult(w1,g63,ng1);
	_ff_sum_2_mults(t4,t3,t4,w0,w1);
	_ff_mult(w0,g63,ng0); _ff_mult(w1,g62,ng1);
	_ff_sum_3_mults(t3,f[3],t2,t3,w0,w1,w);
	_ff_mult(w0,g62,ng0); _ff_mult(w1,g6,ng1);
	_ff_sum_3_mults(t2,f[2],t1,t2,w0,w1,w);
	_ff_mult(w0,g6,ng0);
	_ff_sum_3_mults(t1,f[1],t0,t1,w0,ng1,w);
	_ff_sum_2_mults(t0,f[0],t0,ng0,w);
	_ff_sum_2_mults(s0,w5,g6,t4,t5);
	_ff_mult(w,t5,g6); _ff_neg(w1,w);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_2_mults(s3,t2,t3,s0,w1);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_square(w,t5);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w);
	_ff_sum_2_mults(s0,g[0],t0,s0,w);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 101M
}

static inline void ff_poly_s12_gcd_linear_18_6 (ff_t h[2], ff_t f[12], ff_t g[7])	// assumes g = g6x^6+g5x^5+g3x^3+g1x+g0
{
	register ff_t ng0, ng1, ng3, ng5, g6, g62, g63, g64, g65, t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, s5, w0, w1;
	
	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng3,g[3]); _ff_neg(ng5,g[5]); _ff_set(g6,g[6]); 
	_ff_square(g62,g6); _ff_mult(g63,g62,g6); _ff_square(g64,g62); _ff_mult(g65,g64,g6);
	ff_mult(s2,g62,ng3); _ff_mult(s4,g64,ng1); _ff_mult(s5,g65,ng0);
	_ff_mult(t5,f[17],g6); _ff_addto(t5,ng5);
	_ff_sum_2_mults(t4,f[16],t5,ng5,g62);
	_ff_sum_2_mults(t3,f[15],t4,ng5,g63); _ff_addto(t3,s2);
	_ff_sum_3_mults(t2,f[14],t3,t5,s2,ng5,g64);
	_ff_sum_3_mults(t1,f[13],t2,t4,s2,ng5,g65); _ff_addto(t1,s4);
	_ff_mult(w0,g65,g6);
	_ff_sum_4_mults(t0,f[12],t1,t3,t5,s4,s2,ng5,w0); _ff_addto(t0,s5);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t5,f[11],t0,t2,t4,t5,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t4,f[10],t5,t1,t3,t4,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t3,f[9],t4,t0,t2,t3,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t2,f[8],t3,t5,t1,t2,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t1,f[7],t2,t4,t0,t1,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t0,f[6],t1,t3,t5,t0,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_mult(s1,g6,t1); _ff_mult(s2,g62,t2); _ff_mult(s3,g63,t3); _ff_mult(s4,g64,t4); 
	_ff_sum_5_mults(t5,f[5],t0,s2,s4,s5,t5,ng1,ng3,ng5,w0);
	_ff_sum_4_mults(t4,f[4],s1,s3,s4,ng0,ng1,ng3,w0);
	_ff_sum_4_mults(t3,f[3],t0,s2,s3,ng0,ng1,ng3,w0);
	_ff_sum_3_mults(t2,f[2],s1,s2,ng0,ng1,w0);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,w0);
	_ff_sum_2_mults(t0,f[0],t0,ng0,w0);
	
	_ff_sum_2_mults(s0,ng5,g6,t4,t5);
	_ff_mult(w0,t5,g6); _ff_neg(w1,w0); _ff_square(w0,t5);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_3_mults(s3,g[3],t2,t3,s0,w1,w0);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w0);
	_ff_sum_2_mults(s0,g[0],t0,s0,w0);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 137M
}

static inline void ff_poly_s24_gcd_linear_20_6 (ff_t h[2], ff_t f[20], ff_t g[7])	// assumes f is monic with f5, f10, f15 coeefficients zero (true for Phi_19^f(X,j)
															//  g = g6x^6+g5x^5+g1x+g0 (this will be true for Phi_13^f(X,j) and s^6*Phi_5^f(X,r/s))
{
	register ff_t t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, w, w0, w1, w5, ng0, ng1, g6, g62, g63, g64, g65;
	
	_ff_set(g6,g[6]); _ff_neg(w5,g[5]); _ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]);
	_ff_square(g62,g6); _ff_mult(g63,g62,g6); _ff_square(g64,g62); _ff_mult(g65,g62,g63);
	_ff_mult(t1,f[19],g6); _ff_addto(t1,w5);
	_ff_sum_2_mults(t0,f[18],t1,w5,g62);
	_ff_sum_2_mults(t5,f[17],t0,w5,g63);
	_ff_sum_2_mults(t4,f[16],t5,w5,g64);
	_ff_mult(w1,g64,ng1);  _ff_mult(w0,g65,ng0);
	_ff_mult(t3,t4,w5); _ff_addto(t3,w1);
	_ff_mult(w,g65,g6);
	_ff_sum_3_mults(t2,f[14],t3,t1,w1,w5,w); _ff_addto(t2,w0);
	_ff_mult(w,w,g6);
	_ff_sum_4_mults(t1,f[13],t2,t0,t1,w0,w1,w5,w);
	_ff_mult(w,w,g6);
	_ff_sum_4_mults(t0,f[12],t1,t5,t0,w0,w1,w5,w);
	_ff_mult(w,w,g6);
	_ff_sum_4_mults(t5,f[11],t0,t4,t5,w0,w1,w5,w);
	_ff_sum_3_mults(t4,t5,t3,t4,w0,w1,w5);
	_ff_mult(w,w,g62);
	_ff_sum_4_mults(t3,f[9],t4,t2,t3,w0,w1,w5,w);
	_ff_mult(w,w,g6);
	_ff_sum_4_mults(t2,f[8],t3,t1,t2,w0,w1,w5,w);
	_ff_mult(w,w,g6);
	_ff_sum_4_mults(t1,f[7],t2,t0,t1,w0,w1,w5,w);
	_ff_mult(w,w,g6);
	_ff_sum_4_mults(t0,f[6],t1,t5,t0,w0,w1,w5,w);
	_ff_mult(w,w,g6);
	_ff_sum_3_mults(t5,t0,t4,t5,w0,w1,w5);
	_ff_mult(w0,g64,ng0); _ff_mult(w1,g63,ng1);
	_ff_sum_3_mults(t4,f[4],t3,t4,w0,w1,w);
	_ff_mult(w0,g63,ng0); _ff_mult(w1,g62,ng1);
	_ff_sum_3_mults(t3,f[3],t2,t3,w0,w1,w);
	_ff_mult(w0,g62,ng0); _ff_mult(w1,g6,ng1);
	_ff_sum_3_mults(t2,f[2],t1,t2,w0,w1,w);
	_ff_mult(w0,g6,ng0);
	_ff_sum_3_mults(t1,f[1],t0,t1,w0,ng1,w);
	_ff_sum_2_mults(t0,f[0],t0,ng0,w);
	_ff_sum_2_mults(s0,w5,g6,t4,t5);
	_ff_mult(w,t5,g6); _ff_neg(w1,w);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_2_mults(s3,t2,t3,s0,w1);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_square(w,t5);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w);
	_ff_sum_2_mults(s0,g[0],t0,s0,w);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 111M
}

static inline void ff_poly_s12_gcd_linear_20_6 (ff_t h[2], ff_t f[12], ff_t g[7])	// assumes g = g6x^6+g5x^5+g3x^3+g1x+g0
{
	register ff_t ng0, ng1, ng3, ng5, g6, g62, g63, g64, g65, t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, s5, w0, w1;
	
	_ff_neg(ng0,g[0]); _ff_neg(ng1,g[1]); _ff_neg(ng3,g[3]); _ff_neg(ng5,g[5]); _ff_set(g6,g[6]); 
	_ff_square(g62,g6); _ff_mult(g63,g62,g6); _ff_square(g64,g62); _ff_mult(g65,g64,g6);
	ff_mult(s2,g62,ng3); _ff_mult(s4,g64,ng1); _ff_mult(s5,g65,ng0);
	_ff_mult(t1,f[19],g6); _ff_addto(t1,ng5);
	_ff_sum_2_mults(t0,f[18],t1,ng5,g62);
	_ff_sum_2_mults(t5,f[17],t0,ng5,g63); _ff_addto(t5,s2);
	_ff_sum_3_mults(t4,f[16],t5,t1,s2,ng5,g64);
	_ff_sum_3_mults(t3,f[15],t4,t0,s2,ng5,g65); _ff_addto(t3,s4);
	_ff_mult(w0,g65,g6);
	_ff_sum_4_mults(t2,f[14],t3,t5,t1,s4,s2,ng5,w0); _ff_addto(t2,s5);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t1,f[13],t2,t4,t0,t1,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t0,f[12],t1,t3,t5,t0,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t5,f[11],t0,t2,t4,t5,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t4,f[10],t5,t1,t3,t4,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t3,f[9],t4,t0,t2,t3,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t2,f[8],t3,t5,t1,t2,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t1,f[7],t2,t4,t0,t1,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_sum_5_mults(t0,f[6],t1,t3,t5,t0,s5,s4,s2,ng5,w0);
	ff_mult(w0,w0,g6);
	_ff_mult(s1,g6,t1); _ff_mult(s2,g62,t2); _ff_mult(s3,g63,t3); _ff_mult(s4,g64,t4); 
	_ff_sum_5_mults(t5,f[5],t0,s2,s4,s5,t5,ng1,ng3,ng5,w0);
	_ff_sum_4_mults(t4,f[4],s1,s3,s4,ng0,ng1,ng3,w0);
	_ff_sum_4_mults(t3,f[3],t0,s2,s3,ng0,ng1,ng3,w0);
	_ff_sum_3_mults(t2,f[2],s1,s2,ng0,ng1,w0);
	_ff_sum_3_mults(t1,f[1],t0,s1,ng0,ng1,w0);
	_ff_sum_2_mults(t0,f[0],t0,ng0,w0);
	
	_ff_sum_2_mults(s0,ng5,g6,t4,t5);
	_ff_mult(w0,t5,g6); _ff_neg(w1,w0); _ff_square(w0,t5);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_3_mults(s3,g[3],t2,t3,s0,w1,w0);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w0);
	_ff_sum_2_mults(s0,g[0],t0,s0,w0);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 137M
}


int phi_surface_qgcd_cycle_5_7 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_5, ff_t *phi_7)
{
	ff_t f[9], g[7], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;

	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	// duplicate loops to optimize for sparser Weber-f
	switch ( phi_sparse_factor(phi_5) ) {
	case 2:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_s2_eval_ff (f, phi_7, r2[i]);
				phi5_s2_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_8_6 (h, f, g);				// this is about 10% faster than removing a root and using gcd_linear_8_5
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 3:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_s3_eval_ff (f, phi_7, r2[i]);
				phi5_s3_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_8_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 4:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_s4_eval_ff (f, phi_7, r2[i]);
				phi5_s4_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_8_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 6:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_s6_eval_ff (f, phi_7, r2[i]);
				phi5_s6_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_8_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 8:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_s8_eval_ff (f, phi_7, r2[i]);
				phi5_s8_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_8_6 (h, f, g);		// we could save 3 mults by using g3=0, but we don't bother
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 12:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_s12_eval_ff (f, phi_7, r2[i]);
				phi5_s12_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_s12_gcd_linear_8_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 24:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_s24_eval_ff (f, phi_7, r2[i]);
				phi5_s24_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_s24_gcd_linear_8_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	default:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi7_eval_ff (f, phi_7, r2[i]);
				phi5_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_8_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}

int phi_surface_qgcd_cycle_5_11 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_5, ff_t *phi_11)
{
	ff_t f[13], g[7], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	switch ( phi_sparse_factor(phi_5) ) {
	case 2:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_s2_eval_ff (f, phi_11, r2[i]);
				phi5_s2_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_12_6 (h, f, g);			// this is slightly faster than removing a root
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 3:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_s3_eval_ff (f, phi_11, r2[i]);
				phi5_s3_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_12_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 4:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_s4_eval_ff (f, phi_11, r2[i]);
				phi5_s4_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_12_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 6:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_s6_eval_ff (f, phi_11, r2[i]);
				phi5_s6_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_12_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 8:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_s8_eval_ff (f, phi_11, r2[i]);
				phi5_s8_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_12_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 12:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_s12_eval_ff (f, phi_11, r2[i]);
				phi5_s12_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_s12_gcd_linear_12_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 24:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_s24_eval_ff (f, phi_11, r2[i]);
				phi5_s24_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_s24_gcd_linear_12_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	default:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi11_eval_ff (f, phi_11, r2[i]);
				phi5_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_12_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
//		err_printf("Unhandled sparse factor %d in phi_surface_qgcd_cycle_5_11\n", phi_sparse_factor(phi_5)); exit (0);
	return 1;
}

int phi_surface_qgcd_cycle_5_13 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_5, ff_t *phi_13)
{
	ff_t f[15], g[7], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	switch ( phi_sparse_factor(phi_5) ) {
	case 2:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_s2_eval_ff (f, phi_13, r2[i]);
				phi5_s2_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_14_6 (h, f, g);					// it is essentially the same speed whether we remove a root or not, so for simplicity we don't
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 3:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_s3_eval_ff (f, phi_13, r2[i]);
				phi5_s3_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_14_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 4:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_s4_eval_ff (f, phi_13, r2[i]);
				phi5_s4_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_14_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 6:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_s6_eval_ff (f, phi_13, r2[i]);
				phi5_s6_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_14_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 8:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_s8_eval_ff (f, phi_13, r2[i]);
				phi5_s8_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_14_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 12:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_s12_eval_ff (f, phi_13, r2[i]);
				phi5_s12_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_s12_gcd_linear_14_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 24:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_s24_eval_ff (f, phi_13, r2[i]);
				phi5_s24_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_s24_gcd_linear_14_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	default:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi13_eval_ff (f, phi_13, r2[i]);
				phi5_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_gcd_linear_14_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}

int phi_surface_qgcd_cycle_5_17 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_5, ff_t *phi_17)
{
	ff_t f[19], g[7], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;

	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	switch ( phi_sparse_factor(phi_5) ) {
	case 2:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi17_s2_eval_ff (f, phi_17, r2[i]);
				phi5_s2_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 18, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 3:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi17_s3_eval_ff (f, phi_17, r2[i]);
				phi5_s3_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 18, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 4:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi17_s4_eval_ff (f, phi_17, r2[i]);
				phi5_s4_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 18, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 6:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi17_s6_eval_ff (f, phi_17, r2[i]);
				phi5_s6_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 18, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 8:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi17_s8_eval_ff (f, phi_17, r2[i]);
				phi5_s8_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 18, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 12:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi17_s12_eval_ff (f, phi_17, r2[i]);
				phi5_s12_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_s12_gcd_linear_18_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 24:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi17_s24_eval_ff (f, phi_17, r2[i]);
				phi5_s24_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_s24_gcd_linear_18_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	default:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi17_eval_ff (f, phi_17, r2[i]);
				phi5_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 18, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}


int phi_surface_qgcd_cycle_5_19 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_5, ff_t *phi_19)
{
	ff_t f[21], g[7], h[2], d[BATCH_INVERTS+2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	switch ( phi_sparse_factor(phi_5) ) {
	case 2:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi19_s2_eval_ff (f, phi_19, r2[i]);
				phi5_s2_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 20, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 3:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi19_s3_eval_ff (f, phi_19, r2[i]);
				phi5_s3_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 20, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 4:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi19_s4_eval_ff (f, phi_19, r2[i]);
				phi5_s4_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 20, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 6:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi19_s6_eval_ff (f, phi_19, r2[i]);
				phi5_s6_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 20, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 8:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi19_s8_eval_ff (f, phi_19, r2[i]);
				phi5_s8_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 20, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 12:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi19_s12_eval_ff (f, phi_19, r2[i]);
				phi5_s12_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_s12_gcd_linear_20_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	case 24:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi19_s24_eval_ff (f, phi_19, r2[i]);
				phi5_s24_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_s24_gcd_linear_20_6 (h, f, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
		break;
	default:
		while ( j < n ) {
			m = ( m0 > n-j ? n-j : m0 );
			j0 = j-2;
			for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
				phi19_eval_ff (f, phi_19, r2[i]);
				phi5_qeval_ff (g, phi_5, r[j-1], d[k-1]);
				ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
				ff_poly_gcd_linear_n_5 (h, f, 20, g);
				_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
			}
			if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
			for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
		}
	}
	return 1;
}


int phi_surface_qgcd_cycle_5_p2 (ff_t r[], ff_t *r2, int n, int p2, int e, ff_t *phi_5, ff_t *phi_p2)
{
	ff_t f[PHI_MAX_M+2], g[6], d[BATCH_INVERTS+2], h[2];
	register int i,j,j0,k;
	int m,m0;
	
	i = 0; j = e;
	if ( r != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	m0 = ( r==r2 && e < BATCH_INVERTS ? e : BATCH_INVERTS );
	_ff_set_one(d[0]); _ff_set_one(d[1]);
	while ( j < n ) {
		m = ( m0 > n-j ? n-j : m0 );
		j0 = j-2;
		for ( k = 2 ; k <= m+1 ; i++, j++, k++ ) {
			phi_eval_ff (f, phi_p2, p2, r2[i]);
//printf("Phi_%d(X,%ld) mod %ld ", p2, _ff_get_ui(r[i]), _ff_p); ff_poly_print(f,p2+1);
			phi5_qeval_ff (g, phi_5, r[j-1], d[k-1]);
//printf ("Phi_3(X,%ld/%ld) mod %ld ", _ff_get_ui(r[j-1]), _ff_get_ui(d[k-1]), _ff_p);  ff_poly_print(g,4);
			ff_poly_remove_qroot_6 (g, r+j-2, d+k-2);
//printf ("Removed root %ld/%ld ", _ff_get_ui(r[j-2]), _ff_get_ui(d[k-2]));  ff_poly_print(g,3);
			ff_poly_gcd_linear_n_5 (h, f, p2+1, g);
//printf ("GCD "); ff_poly_print(h,1);
			_ff_neg(r[j],h[0]);  _ff_set(d[k],h[1]);
		}
		if ( ! ff_parallel_invert_check(d+2,d+2,m) ) return 0;
		for ( k = 2 ; k <= m+1 ; k++ ) ff_mult(r[j0+k],r[j0+k],d[k]);
	}
	return 1;
}


/*
	Enumerates a surface p1-cycle using gcds with Phi_p2, where alpha_p2=alpha_p1^e and |alpha_p1|=n
	(Here alpha_a denotes the class represented by the pos def reduced primeform <a,b,c>).
	Assumes n > e > 1 and that roots[0],...,roots[e-1] are already present.
*/
void phi_surface_gcd_cycle (ff_t r1[], ff_t r2[], int p1, int n, int p2, int e, ff_t *phi1, ff_t *phi2)
{
	ff_t f[PHI_MAX_M+2], g[PHI_MAX_M+2], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;
//	clock_t start,end;

//dbg_printf ("gcd_cycle(%d,%d) n=%d, e=%d\n", p1, p2, n, e);
	
	if ( p1==2 ) {
		switch (p2) {
		case  3:  if ( phi_surface_qgcd_cycle_2_3 (r1, r2, n, e, phi1, phi2) ) return; break;
		case  5:  if ( phi_surface_qgcd_cycle_2_5 (r1, r2,  n, e, phi1, phi2) ) return; break;
		case  7:  if ( phi_surface_qgcd_cycle_2_7 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 11: if ( phi_surface_qgcd_cycle_2_11 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 13: if ( phi_surface_qgcd_cycle_2_13 (r1, r2, n, e, phi1, phi2) ) return; break;
		default: if ( phi_surface_qgcd_cycle_2_p2 (r1, r2, n, p2, e, phi1, phi2) ) return;	break;
		}
	} else if ( p1 == 3 ) {
		switch (p2) {
		case  5:  if ( phi_surface_qgcd_cycle_3_5 (r1, r2, n, e, phi1, phi2) )  return; break;
		case  7:  if ( phi_surface_qgcd_cycle_3_7 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 11: if ( phi_surface_qgcd_cycle_3_11 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 13: if ( phi_surface_qgcd_cycle_3_13 (r1, r2, n, e, phi1, phi2) ) return; break;
		default: if ( phi_surface_qgcd_cycle_3_p2 (r1, r2, n, p2, e, phi1, phi2) ) return;	break;
		}
	} else if ( p1 == 5 ) {
		switch (p2) {
		case  7:  if ( phi_surface_qgcd_cycle_5_7 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 11: if ( phi_surface_qgcd_cycle_5_11 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 13: if ( phi_surface_qgcd_cycle_5_13 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 17: if ( phi_surface_qgcd_cycle_5_17 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 19: if ( phi_surface_qgcd_cycle_5_19 (r1, r2, n, e, phi1, phi2) ) return; break;
		default: if ( phi_surface_qgcd_cycle_5_p2 (r1, r2, n, p2, e, phi1, phi2) ) return;	break;
		}
	} else if ( p1 == 7 ) {
		if ( p2 == 11 ) { if ( phi_surface_gcd_cycle_7_11 (r1, r2, n, e, phi1, phi2) ) return; }
		if ( p2 == 13 ) { if ( phi_surface_gcd_cycle_7_13 (r1, r2, n, e, phi1, phi2) ) return; }
		if ( p2 == 17 ) { if ( phi_surface_gcd_cycle_7_17 (r1, r2, n, e, phi1, phi2) ) return; }
		if ( p2 == 19 ) { if ( phi_surface_gcd_cycle_7_19 (r1, r2, n, e, phi1, phi2) ) return; }
	}
//start = clock();
	i1=j2= 0; i2=j1=e-1; 
	if ( r1 != r2 ) { i1=j1+1;  i2 = n-1; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	do {
		phi_eval_ff (f, phi2, p2, r2[i1]);
		phi_eval_ff (g, phi1, p1, r1[j1]);
		ff_poly_remove_root (g, g, p1+1, r1+j1-1);
		ff_poly_clf (h1, f, p2+1, g, p1);
		phi_eval_ff (f, phi2, p2, r2[i2]);
		phi_eval_ff (g, phi1, p1, r1[j2]);
		k = j2+1;  if ( k==n ) k = 0;
		ff_poly_remove_root (g, g, p1+1, r1+k);
		ff_poly_clf (h2, f, p2+1, g, p1);
		i1++; i2--; if ( i2 < 0 ) i2 = n-1;
		j1++; j2--; if ( j2 < 0 ) j2 = n-1;
		_ff_mult(t0,h1[1],h2[1]); _ff_invert(t1,t0);
		_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r1[j1],t2);
		_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r1[j2],t2);
	} while ( j2 > j1+1 );
//end = clock();
//printf ("gcd %d-cycle of length %d-%d using aux_ell=%d took %ld msecs\n", p1, n, e, p2, delta_msecs(start,end));
}

void phi_surface_gcd_path (ff_t r1[], ff_t r2[], int p1, int n, int p2, int e, ff_t *phi1, ff_t *phi2)
{
	ff_t f[PHI_MAX_M+2], g[PHI_MAX_M+2], h[2];
	register int i, j;

dbg_printf ("gcd_path(%d,%d) n=%d, e=%d\n", p1, p2, n, e);
	
	if ( p1==2 ) {
		switch (p2) {
		case  3:  if ( phi_surface_qgcd_cycle_2_3 (r1, r2, n, e, phi1, phi2) ) return; break;
		case  5:  if ( phi_surface_qgcd_cycle_2_5 (r1, r2,  n, e, phi1, phi2) ) return; break;
		case  7:  if ( phi_surface_qgcd_cycle_2_7 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 11: if ( phi_surface_qgcd_cycle_2_11 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 13: if ( phi_surface_qgcd_cycle_2_13 (r1, r2, n, e, phi1, phi2) ) return; break;
		default: if ( phi_surface_qgcd_cycle_2_p2 (r1, r2, n, p2, e, phi1, phi2) ) return;	break;	// this shouldn't be used unless the 2-volcano has height greater than 1
		}
	} else if ( p1 == 3 ) {
		switch (p2) {
		case  5:  if ( phi_surface_qgcd_cycle_3_5 (r1, r2, n, e, phi1, phi2) )  return; break;
		case  7:  if ( phi_surface_qgcd_cycle_3_7 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 11: if ( phi_surface_qgcd_cycle_3_11 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 13: if ( phi_surface_qgcd_cycle_3_13 (r1, r2, n, e, phi1, phi2) ) return; break;
		default: if ( phi_surface_qgcd_cycle_3_p2 (r1, r2, n, p2, e, phi1, phi2) ) return;	break;	// this shouldn't be used unless the 3-volcano has height greater than 1
		}
	} else if ( p1 == 5 ) {
		switch (p2) {
		case  7: if ( phi_surface_qgcd_cycle_5_7 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 11: if ( phi_surface_qgcd_cycle_5_11 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 13: if ( phi_surface_qgcd_cycle_5_13 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 17: if ( phi_surface_qgcd_cycle_5_17 (r1, r2, n, e, phi1, phi2) ) return; break;
		case 19: if ( phi_surface_qgcd_cycle_5_19 (r1, r2, n, e, phi1, phi2) ) return; break;
		}
	}
	i = 0; j = e;
	if ( r1 != r2 ) { i=j; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	while ( j < n ) {
		phi_eval_ff (f, phi2, p2, r2[i]);
		phi_eval_ff (g, phi1, p1, r1[j-1]);
		ff_poly_remove_root (g, g, p1+1, r1+j-2);
		ff_poly_clf (h, f, p2+1, g, p1);
		ff_poly_roots_d1 (r1+j, h);
		i++; j++;
	}
}



/*
	These functions have been superseded by qgcd_cycle functions that batch inversions

void phi_surface_gcd_cycle_2_3 (ff_t r[], int n, int e, ff_t *phi_2, ff_t *phi_3)
{
	ff_t f[5], g[4], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;
	
	i1=j2= 0; i2=j1=e-1;
	do {
		phi3_eval_ff (f, phi_3, r[i1]);
		phi2_eval_ff (g, phi_2, r[j1]);
		ff_poly_remove_root_d3 (g, g, r+j1-1);
		ff_poly_mgcd_linear_4_2 (h1, f, g);
		phi3_eval_ff (f, phi_3, r[i2]);
		phi2_eval_ff (g, phi_2, r[j2]);
		k = j2+1;  if ( k==n ) k = 0;
		ff_poly_remove_root_d3 (g, g, r+k);
		ff_poly_mgcd_linear_4_2 (h2, f, g);
		i1++; i2--; if ( i2 < 0 ) i2 = n-1;
		j1++; j2--; if ( j2 < 0 ) j2 = n-1;
		_ff_mult(t0,h1[1],h2[1]); _ff_invert(t1,t0);
		_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
		_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
	} while ( j2 > j1+1 );
}

void phi_surface_gcd_cycle_3_5 (ff_t r[], int n, int e, ff_t *phi_3, ff_t *phi_5)
{
	ff_t f[7], g[5], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;

	i1=j2= 0; i2=j1=e-1;
	do {
		phi5_eval_ff (f, phi_5, r[i1]);
		phi3_eval_ff (g, phi_3, r[j1]);
		ff_poly_remove_root_d4 (g, g, r+j1-1);
		ff_poly_mgcd_linear_6_3 (h1, f, g);
		phi5_eval_ff (f, phi_5, r[i2]);
		phi3_eval_ff (g, phi_3, r[j2]);
		k = j2+1;  if ( k==n ) k = 0;
		ff_poly_remove_root_d4 (g, g, r+k);
		ff_poly_mgcd_linear_6_3 (h2, f, g);
		i1++; i2--; if ( i2 < 0 ) i2 = n-1;
		j1++; j2--; if ( j2 < 0 ) j2 = n-1;
		_ff_mult(t0,h1[1],h2[1]); _ff_invert(t1,t0);
		_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
		_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
	} while ( j2 > j1+1 );
}


static inline void ff_poly_f_mgcd_linear_8_6 (ff_t h[2], ff_t f[8], ff_t g[7])	// assumes f = x^8+f7x^7+f4x^4+f1x+f0,  g = g6x^6+g5x^5+g1x+g0 (this will be true for Phi_7^f(X,j) and s^6*Phi_5^f(X,r/s))
{
	register ff_t g5, g52, s0,s1,s2, t0, t1, t2, t3, t4,t5;
	
	_ff_set(g5,g[5]); _ff_sub(s0,f[7],g5); _ff_mult(s1,g[1],g5); _ff_subfrom(s1,g[0]); _ff_mult(s2,g[0],g5);  _ff_square(g52,g5);
	_ff_mult(t5,g52,s0);  _ff_set(t4,f[4]); _ff_neg(t3,g[1]);
	_ff_mult(t0,g[1],s0); _ff_addto(t0,g[0]); _ff_neg(t2,t0);
	_ff_mult(t0,s0,s1); _ff_add(t1,f[1],t0);
	_ff_mult(s1,s0,s2); _ff_add(t0,f[0],s1);
	ff_poly_gcd_linear_6_5_reg (h,g[0],g[1],g[2],g[3],g[4],g[5],g[6],t0,t1,t2,t3,t4,t5);
	// Total of 61M+38A (33 redc)
}

static inline void ff_poly_f_mgcd_linear_12_6 (ff_t h[2], ff_t f[12], ff_t g[7])	// assumes f = x^12+f11x^11+f9x^9+f7x^7+f5x^5+f3x^3+f1x+f0,  g = x^6+g5x^5+g1x+g0 (this will be true for Phi_11^f(X,j) and Phi_5^f(X,j)
{
	register ff_t t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, w, w0, w1, w5, g6;
	
	_ff_neg(w0,g[0]); _ff_neg(w1,g[1]); _ff_neg(w5,g[5]);
	_ff_add(t5,f[11],w5);
	_ff_mult(t4,t5,w5);
	_ff_mult(t3,t4,w5); _ff_addto(t3,f[9]);
	_ff_mult(t2,t3,w5);
	_ff_mult(t1,t2,w5); _ff_addto(t1,w1); _ff_addto(t1,f[7]);
	_ff_sum_2_mults(t0,t1,t5,w1,w5); _ff_addto(t0,w0);
	_ff_sum_3_mults(t5,t0,t4,t5,w0,w1,w5); _ff_addto(t5,f[5]);
	_ff_sum_2_mults(t4,t3,t4,w0,w1);
	_ff_sum_2_mults(t3,t2,t3,w0,w1); _ff_addto(t3,f[3]);
	_ff_sum_2_mults(t2,t1,t2,w0,w1);
	_ff_sum_2_mults(t1,t0,t1,w0,w1); _ff_addto(t1,f[1]);
	ff_mult(t0,t0,w0); _ff_addto(t0,f[0]);
	_ff_mult(s0,t5,w5); _ff_addto(s0,t4); _ff_neg(w1,t5);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_2_mults(s3,t2,t3,s0,w1);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_square(w,t5);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w);
	_ff_sum_2_mults(s0,g[0],t0,s0,w);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 67M+520A (35 redc)
}


static inline void ff_poly_f_mgcd_linear_14_6 (ff_t h[2], ff_t f[14], ff_t g[7])	// assumes f monic with coeff f3,f5,f7,f9,f11=0,  g = x^6+g5x^5+g1x+g0 (this will be true for Phi_13^f(X,j) and Phi_5^f(X,j)
{
	register ff_t t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, w, w0, w1, w5, g6;
	
	_ff_neg(w0,g[0]); _ff_neg(w1,g[1]); _ff_neg(w5,g[5]);
	_ff_add(t1,f[13],w5);
	_ff_mult(t0,t1,w5); _ff_addto(t0,f[12]);
	_ff_mult(t5,t0,w5);
	_ff_mult(t4,t5,w5); _ff_addto(t4,f[10]);
	_ff_mult(t3,t4,w5); _ff_addto(t3,w1);
	_ff_sum_2_mults(t2,t3,t1,w1,w5);  _ff_addto(t2,w0); _ff_addto(t2,f[8]);
	_ff_sum_3_mults(t1,t2,t0,t1,w0,w1,w5);
	_ff_sum_3_mults(t0,t1,t5,t0,w0,w1,w5); _ff_addto(t0,f[6]);
	_ff_sum_3_mults(t5,t0,t4,t5,w0,w1,w5);
	_ff_sum_2_mults(t4,t3,t4,w0,w1); _ff_addto(t4,f[4]);
	_ff_sum_2_mults(t3,t2,t3,w0,w1);
	_ff_sum_2_mults(t2,t1,t2,w0,w1); _ff_addto(t2,f[2]);
	_ff_sum_2_mults(t1,t0,t1,w0,w1); _ff_addto(t1,f[1]);
	ff_mult(t0,t0,w0); _ff_addto(t0,f[0]);
	_ff_mult(s0,t5,w5); _ff_addto(s0,t4); _ff_neg(w1,t5);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_2_mults(s3,t2,t3,s0,w1);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_square(w,t5);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w);
	_ff_sum_2_mults(s0,g[0],t0,s0,w);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 67M+52A (35 redc)
}


static inline void ff_poly_f_mgcd_linear_18_6 (ff_t h[2], ff_t f[18], ff_t g[7])	// assumes f monic with coeff f4,f7,f11,f14=0,  g = x^6+g5x^5+g1x+g0 (this will be true for Phi_17^f(X,j) and Phi_5^f(X,j)
{
	register ff_t t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, w, w0, w1, w5, g6;
	
	_ff_neg(w0,g[0]); _ff_neg(w1,g[1]); _ff_neg(w5,g[5]);
	_ff_add(t5,f[17],w5);
	_ff_mult(t4,t5,w5); _ff_addto(t4,f[16]);
	_ff_mult(t3,t4,w5); _ff_addto(t3,f[15]);
	_ff_mult(t2,t3,w5);
	_ff_mult(t1,t2,w5); _ff_addto(t1,w1); _ff_addto(t1,f[13]);
	_ff_sum_2_mults(t0,t1,t5,w1,w5);  _ff_addto(t0,w0); _ff_addto(t0,f[12]);
	_ff_sum_3_mults(t5,t0,t4,t5,w0,w1,w5);
	_ff_sum_3_mults(t4,t5,t3,t4,w0,w1,w5); _ff_addto(t4,f[10]);
	_ff_sum_3_mults(t3,t4,t2,t3,w0,w1,w5); _ff_addto(t3,f[9]);
	_ff_sum_3_mults(t2,t3,t1,t2,w0,w1,w5); _ff_addto(t2,f[8]);
	_ff_sum_3_mults(t1,t2,t0,t1,w0,w1,w5);
	_ff_sum_3_mults(t0,t1,t5,t0,w0,w1,w5); _ff_addto(t0,f[6]);
	_ff_sum_3_mults(t5,t0,t4,t5,w0,w1,w5); _ff_addto(t5,f[5]);
	_ff_sum_2_mults(t4,t3,t4,w0,w1);
	_ff_sum_2_mults(t3,t2,t3,w0,w1); _ff_addto(t3,f[3]);
	_ff_sum_2_mults(t2,t1,t2,w0,w1); _ff_addto(t2,f[2]);
	_ff_sum_2_mults(t1,t0,t1,w0,w1); _ff_addto(t1,f[1]);
	ff_mult(t0,t0,w0); _ff_addto(t0,f[0]);
	_ff_mult(s0,t5,w5); _ff_addto(s0,t4); _ff_neg(w1,t5);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_2_mults(s3,t2,t3,s0,w1);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_square(w,t5);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w);
	_ff_sum_2_mults(s0,g[0],t0,s0,w);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 87M+71A (42 redc)
}


static inline void ff_poly_f_mgcd_linear_20_6 (ff_t h[2], ff_t f[18], ff_t g[7])	// assumes f monic with coeff f5,f10,f15=0,  g = x^6+g5x^5+g1x+g0 (this will be true for Phi_19^f(X,j) and Phi_5^f(X,j)
{
	register ff_t t0, t1, t2, t3, t4, t5, s0, s1, s2, s3, s4, w, w0, w1, w5, g6;
	
	_ff_neg(w0,g[0]); _ff_neg(w1,g[1]); _ff_neg(w5,g[5]);
	_ff_add(t1,f[19],w5);
	_ff_mult(t0,t1,w5); _ff_addto(t0,f[18]);
	_ff_mult(t5,t0,w5); _ff_addto(t5,f[17]);
	_ff_mult(t4,t5,w5); _ff_addto(t4,f[16]);
	_ff_mult(t3,t4,w5); _ff_addto(t3,w1);
	_ff_sum_2_mults(t2,t3,t1,w1,w5);  _ff_addto(t2,w0); _ff_addto(t2,f[14]);
	_ff_sum_3_mults(t1,t2,t0,t1,w0,w1,w5); _ff_addto(t1,f[13]);
	_ff_sum_3_mults(t0,t1,t5,t0,w0,w1,w5); _ff_addto(t0,f[12]);
	_ff_sum_3_mults(t5,t0,t4,t5,w0,w1,w5); _ff_addto(t5,f[11]);
	_ff_sum_3_mults(t4,t5,t3,t4,w0,w1,w5);
	_ff_sum_3_mults(t3,t4,t2,t3,w0,w1,w5); _ff_addto(t3,f[9]);
	_ff_sum_3_mults(t2,t3,t1,t2,w0,w1,w5); _ff_addto(t2,f[8]);
	_ff_sum_3_mults(t1,t2,t0,t1,w0,w1,w5); _ff_addto(t1,f[7]);
	_ff_sum_3_mults(t0,t1,t5,t0,w0,w1,w5); _ff_addto(t0,f[6]);
	_ff_sum_3_mults(t5,t0,t4,t5,w0,w1,w5);
	_ff_sum_2_mults(t4,t3,t4,w0,w1); _ff_addto(t4,f[4]);
	_ff_sum_2_mults(t3,t2,t3,w0,w1); _ff_addto(t3,f[3]);
	_ff_sum_2_mults(t2,t1,t2,w0,w1); _ff_addto(t2,f[2]);
	_ff_sum_2_mults(t1,t0,t1,w0,w1); _ff_addto(t1,f[1]);
	ff_mult(t0,t0,w0); _ff_addto(t0,f[0]);
	_ff_mult(s0,t5,w5); _ff_addto(s0,t4); _ff_neg(w1,t5);
	_ff_sum_2_mults(s4,t3,t4,s0,w1);
	_ff_sum_2_mults(s3,t2,t3,s0,w1);
	_ff_sum_2_mults(s2,t1,t2,s0,w1);
	_ff_square(w,t5);
	_ff_sum_3_mults(s1,g[1],t0,t1,s0,w1,w);
	_ff_sum_2_mults(s0,g[0],t0,s0,w);
	ff_poly_gcd_linear_5_4_reg (h,t0,t1,t2,t3,t4,t5,s0,s1,s2,s3,s4);
	// 87M+71A (42 redc)
}

int  phi_surface_gcd_cycle_5_7 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_5, ff_t *phi_7)
{
	ff_t f[9], g[7], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;
	
	i1=j2= 0; i2=j1=e-1;
	if ( r != r2 ) { i1=j1+1;  i2 = n-1; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	// duplicate loops to optimize for sparser Weber-f
	switch ( phi_sparse_factor(phi_5) ) {
	case 2:
		do {
			phi7_s2_eval_ff (f, phi_7, r2[i1]);
			phi5_s2_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_m_5 (h1, f, 8, g);
			phi7_s2_eval_ff (f, phi_7, r2[i2]);
			phi5_s2_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_m_5 (h2, f, 8, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 3:
		do {
			phi7_s3_eval_ff (f, phi_7, r2[i1]);
			phi5_s3_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_m_5 (h1, f, 8, g);
			phi7_s3_eval_ff (f, phi_7, r2[i2]);
			phi5_s3_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_m_5 (h2, f, 8, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 4:
		do {
			phi7_s4_eval_ff (f, phi_7, r2[i1]);
			phi5_s4_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_m_5 (h1, f, 8, g);
			phi7_s4_eval_ff (f, phi_7, r2[i2]);
			phi5_s4_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_m_5 (h2, f, 8, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 6:
		do {
			phi7_s6_eval_ff (f, phi_7, r2[i1]);
			phi5_s6_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_m_5 (h1, f, 8, g);
			phi7_s6_eval_ff (f, phi_7, r2[i2]);
			phi5_s6_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_m_5 (h2, f, 8, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 8:
		do {
			phi7_s8_eval_ff (f, phi_7, r2[i1]);
			phi5_s8_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_m_5 (h1, f, 8, g);
			phi7_s8_eval_ff (f, phi_7, r2[i2]);
			phi5_s8_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_m_5 (h2, f, 8, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 12:
		do {
			phi7_s12_eval_ff (f, phi_7, r2[i1]);
			phi5_s12_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_m_5 (h1, f, 8, g);
			phi7_s12_eval_ff (f, phi_7, r2[i2]);
			phi5_s12_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_m_5 (h2, f, 8, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 24:
		do {
			phi7_s24_eval_ff (f, phi_7, r2[i1]);
			phi5_s24_eval_ff (g, phi_5, r[j1]);
			ff_poly_f_mgcd_linear_8_6 (h1, f, g);
			phi7_s24_eval_ff (f, phi_7, r2[i2]);
			phi5_s24_eval_ff (g, phi_5, r[j2]);
			ff_poly_f_mgcd_linear_8_6 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	default:
		do {
			phi7_eval_ff (f, phi_7, r2[i1]);
			phi5_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_m_5 (h1, f, 8, g);
			phi7_eval_ff (f, phi_7, r2[i2]);
			phi5_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_m_5 (h2, f, 8, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	}
	return 1;
}

int  phi_surface_gcd_cycle_5_11 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_5, ff_t *phi_11)
{
	ff_t f[13], g[7], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;
	
	i1=j2= 0; i2=j1=e-1;
	if ( r != r2 ) { i1=j1+1;  i2 = n-1; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	// duplicate loops to optimize for sparser Weber-f
	switch ( phi_sparse_factor(phi_5) ) {
	case 2:
		do {
			phi11_s2_eval_ff (f, phi_11, r2[i1]);
			phi5_s2_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_n_5 (h1, f, 12, g);
			phi11_s2_eval_ff (f, phi_11, r2[i2]);
			phi5_s2_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_n_5 (h2, f, 12, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 3:
		do {
			phi11_s3_eval_ff (f, phi_11, r2[i1]);
			phi5_s3_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_n_5 (h1, f, 12, g);
			phi11_s3_eval_ff (f, phi_11, r2[i2]);
			phi5_s3_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_n_5 (h2, f, 12, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 4:
		do {
			phi11_s4_eval_ff (f, phi_11, r2[i1]);
			phi5_s4_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_n_5 (h1, f, 12, g);
			phi11_s4_eval_ff (f, phi_11, r2[i2]);
			phi5_s4_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_n_5 (h2, f, 12, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 6:
		do {
			phi11_s6_eval_ff (f, phi_11, r2[i1]);
			phi5_s6_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_n_5 (h1, f, 12, g);
			phi11_s6_eval_ff (f, phi_11, r2[i2]);
			phi5_s6_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_n_5 (h2, f, 12, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 8:
		do {
			phi11_s8_eval_ff (f, phi_11, r2[i1]);
			phi5_s8_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_n_5 (h1, f, 12, g);
			phi11_s8_eval_ff (f, phi_11, r2[i2]);
			phi5_s8_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_n_5 (h2, f, 12, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 12:
		do {
			phi11_s12_eval_ff (f, phi_11, r2[i1]);
			phi5_s12_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_n_5 (h1, f, 12, g);
			phi11_s12_eval_ff (f, phi_11, r2[i2]);
			phi5_s12_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_n_5 (h2, f, 12, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	case 24:
		do {
			phi11_s24_eval_ff (f, phi_11, r2[i1]);
			phi5_s24_eval_ff (g, phi_5, r[j1]);
			ff_poly_f_mgcd_linear_12_6 (h1, f, g);
			phi11_s24_eval_ff (f, phi_11, r2[i2]);
			phi5_s24_eval_ff (g, phi_5, r[j2]);
			ff_poly_f_mgcd_linear_12_6 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
		break;
	default:
		do {
			phi11_eval_ff (f, phi_11, r2[i1]);
			phi5_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_n_5 (h1, f, 12, g);
			phi11_eval_ff (f, phi_11, r2[i2]);
			phi5_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_n_5 (h2, f, 12, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]); 	if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	}
	return 1;
}

int  phi_surface_gcd_cycle_5_13 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_5, ff_t *phi_13)
{
	ff_t f[15], g[7], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;
	
	i1=j2= 0; i2=j1=e-1;
	if ( r != r2 ) { i1=j1+1;  i2 = n-1; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	// duplicate loops to optimize for sparser Weber-f
	if (  phi_sparse_factor(phi_5)==24 ) {
		do {
			phi13_s24_eval_ff (f, phi_13, r2[i1]);
			phi5_s24_eval_ff (g, phi_5, r[j1]);
			ff_poly_f_mgcd_linear_14_6 (h1, f, g);
			phi13_s24_eval_ff (f, phi_13, r2[i2]);
			phi5_s24_eval_ff (g, phi_5, r[j2]);
			ff_poly_f_mgcd_linear_14_6 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);
			_ff_invert(t1,t0); if ( _ff_zero(t0) ) return 0;
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	} else {
		do {
			phi13_eval_ff (f, phi_13, r2[i1]);
			phi5_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_n_5 (h1, f, 14, g);
			phi13_eval_ff (f, phi_13, r2[i2]);
			phi5_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_n_5 (h2, f, 14, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	}
	return 1;
}


int  phi_surface_gcd_cycle_5_17 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_5, ff_t *phi_17)
{
	ff_t f[19], g[7], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;
	
	i1=j2= 0; i2=j1=e-1;
	if ( r != r2 ) { i1=j1+1;  i2 = n-1; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	// duplicate loops to optimize for sparser Weber-f
	if (  phi_sparse_factor(phi_5)==24 ) {
		do {
			phi17_s24_eval_ff (f, phi_17, r2[i1]);
			phi5_s24_eval_ff (g, phi_5, r[j1]);
			ff_poly_f_mgcd_linear_18_6 (h1, f, g);
			phi17_s24_eval_ff (f, phi_17, r2[i2]);
			phi5_s24_eval_ff (g, phi_5, r[j2]);
			ff_poly_f_mgcd_linear_18_6 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	} else {
		do {
			phi_eval_ff (f, phi_17, 17, r2[i1]);
			phi5_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_n_5 (h1, f, 18, g);
			phi_eval_ff (f, phi_17, 17, r2[i2]);
			phi5_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_n_5 (h2, f, 18, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	}
	return 1;
}


int  phi_surface_gcd_cycle_5_19 (ff_t r[], ff_t *r2, int n, int e, ff_t *phi_5, ff_t *phi_19)
{
	ff_t f[21], g[7], h1[2], h2[2];
	register ff_t t0, t1,t2;
	register int i1, i2, j1, j2, k;
	
	i1=j2= 0; i2=j1=e-1;
	if ( r != r2 ) { i1=j1+1;  i2 = n-1; }		// hack, if r1 != r2 we r2 actually points to a p2-isogenous parallel p1-path.  e should be 2 in this case
	if (  phi_sparse_factor(phi_5)==24 ) {
		do {
			phi19_s24_eval_ff (f, phi_19, r2[i1]);
			phi5_s24_eval_ff (g, phi_5, r[j1]);
			ff_poly_f_mgcd_linear_20_6 (h1, f, g);
			phi19_s24_eval_ff (f, phi_19, r2[i2]);
			phi5_s24_eval_ff (g, phi_5, r[j2]);
			ff_poly_f_mgcd_linear_20_6 (h2, f, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	} else {
		do {
			phi_eval_ff (f, phi_19, 19, r2[i1]);
			phi5_eval_ff (g, phi_5, r[j1]);
			ff_poly_remove_root (g, g, 6, r+j1-1);
			ff_poly_mgcd_linear_n_5 (h1, f, 20, g);
			phi_eval_ff (f, phi_19, 19, r2[i2]);
			phi5_eval_ff (g, phi_5, r[j2]);
			k = j2+1;  if ( k==n ) k = 0;
			ff_poly_remove_root (g, g, 6, r+k);
			ff_poly_mgcd_linear_n_5 (h2, f, 20, g);
			i1++; i2--; if ( i2 < 0 ) i2 = n-1;
			j1++; j2--; if ( j2 < 0 ) j2 = n-1;
			_ff_mult(t0,h1[1],h2[1]);  if ( _ff_zero(t0) ) return 0;
			_ff_invert(t1,t0);
			_ff_mult(t0,t1,h2[1]); _ff_mult(t2,t0,h1[0]); _ff_neg(r[j1],t2);
			_ff_mult(t0,t1,h1[1]); _ff_mult(t2,t0,h2[0]);  _ff_neg(r[j2],t2);
		} while ( j2 > j1+1 );
	}
	return 1;
}
*/
