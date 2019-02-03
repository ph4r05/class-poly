#ifndef _PHI_EVAL_INCLUDE_
#define _PHI_EVAL_INCLUDE_

#include <stdio.h>
#include "phi_poly.h"

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

// note that m is assumed to be prime in this module

void phi_eval_ff (ff_t f[], ff_t phi[], int m, ff_t J);

static inline void _phi_eval_ff (ff_t f[], ff_t phi[], int m, ff_t J)
{
	ff_t y[PHI_MAX_M+2];
	register int i;
	
	_ff_set_one(y[0]); _ff_set (y[1],J); for ( i = 2 ; i <= m+1 ; i++ ) ff_mult(y[i],y[i-1],J);		// Set y[i]=j^i for i from 0 to m+1
	for ( i = 0 ; i <= m ; i++ ) ff_dot_product(f+i,y,phi+(m+1)*i,m+1);
	_ff_set_one(f[m+1]); _ff_addto(f[0],y[m+1]); 
}

static inline void _phi_s_eval_ff (ff_t f[], ff_t phi[], int m, ff_t J)
{
	ff_t x[PHI_MAX_M+1], y[PHI_MAX_M+PHI_MAX_SPARSE_FACTOR];
	register ff_t t;
	register int i,j,j0,k,n,s;
	
	s = phi_sparse_factor(phi);
	n = m/s+1;
	_ff_set_one(y[0]); _ff_set (y[n],J); _ff_set(t,J);
	for ( i = 2 ; i <= m ; i++ ) {  _ff_mult(t,t,J);  _ff_set(y[(i%s)*n+i/s],t); }
	for ( i = 0 ; i <= m ; i++ ) {
		j0 = (m*(m+1-i))%s;												// j0 is the least b satisfying a + m*b = m+1 mod s (note that m^2=1 mod s for every valid s|24).
		for ( j=j0, k=0 ; j <= m ; j+= s, k++ ) _ff_set(x[k],phi[(m+1)*i+j]);			// this is independent of J and could be done ahead of time, but actually it appears to be slightly faster to copy anyway
		 ff_dot_product(f+i,y+j0*n,x,k); 
	}
	_ff_mult(t,t,J);  _ff_addto(f[0],t);
	_ff_set_one(f[m+1]);
//printf ("%d-sparse phi_%d(X,%ld) mod %ld = ", s, m, _ff_get_ui(J), _ff_p); ff_poly_print(f,m+1);
}

// We rely throughout on the fact that the X^ell*Y^ell coeff is -1

static inline void phi2_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, t1;
	
	_ff_square(J2,J);  _ff_add(t1,J,phi[2]);
	_ff_sum_2_mults(t1,J,J2,t1,phi[1]);       _ff_add(f[0],t1,phi[0]);
	_ff_sum_2_mults(t1,J,J2,phi[5],phi[4]); _ff_add(f[1],t1,phi[3]);
	_ff_mult(t1,J,phi[7]); _ff_addto(t1,phi[6]); _ff_sub(f[2],t1,J2);
	_ff_set_one(f[3]);
	// 6M+7A (4 redc)
}

static inline void phi2_s3_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t t0,t1;
	
	_ff_square(t0,J);
	_ff_mult(t1,t0,J);  _ff_add(f[0],t1,phi[0]);
	_ff_mult(f[1],J,phi[4]);
	_ff_neg(f[2],t0);
	_ff_set_one(f[3]);
	// 3M+2A (3 redc)
}

static inline void phi3_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, t1;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_add(t1,J,phi[3]);
	_ff_sum_3_mults(t1,J,J2,J3,t1,phi[2],phi[1]);  _ff_add(f[0],t1,phi[0]);
	_ff_sum_3_mults(t1,J,J2,J3,phi[7],phi[6],phi[5]); _ff_add(f[1],t1,phi[4]);
	_ff_sum_3_mults(t1,J,J2,J3,phi[11],phi[10],phi[9]); _ff_add(f[2],t1,phi[8]);
	_ff_sum_2_mults(t1,J,J2,phi[14],phi[13]); _ff_addto(t1,phi[12]); _ff_sub(f[3],t1,J3);
	_ff_set_one(f[4]);
	// 13M+13A (6 redc)
}

static inline void phi3_s2_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, t1;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2); _ff_add(t1,J2,phi[2]);		// in fact, the [0,2] coeff appears to always be zero
	ff_mult(t1,t1,J2);  _ff_add(f[0],t1,phi[0]);
	_ff_sum_2_mults(f[1],J,J3,phi[7],phi[5]);
	_ff_mult(f[2],J2,phi[10]);
	_ff_mult(t1,J,phi[13]);  _ff_sub(f[3],t1,J3);
	_ff_set_one(f[4]);
	// 7M+4A (6 redc)
}

static inline void phi3_s4_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, t1;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2);
	_ff_mult(t1,J3,J); _ff_add(f[0],t1,phi[0]);
	_ff_mult(f[1],J,phi[5]);
	_ff_mult(f[2],J2,phi[10]);
	_ff_neg(f[3],J3);	
	_ff_set_one(f[4]);
	// 5M+2A (5 redc)
}

static inline void phi3_s8_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2);
	_ff_mult(f[0],J3,J);
	_ff_mult(f[1],J,phi[5]);
	_ff_set_zero(f[2]);
	_ff_neg(f[3],J3);	
	_ff_set_one(f[4]);
	// 4M+1A (4 redc)
}

static inline void phi5_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, J4, J5, t1;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_square(J4,J2); _ff_mult(J5,J4,J); _ff_add(t1,J,phi[5]);
	_ff_sum_5_mults(t1,J,J2,J3,J4,J5,t1,phi[4],phi[3],phi[2],phi[1]);  _ff_add(f[0],t1,phi[0]);
	_ff_sum_5_mults(t1,J,J2,J3,J4,J5,phi[11],phi[10],phi[9],phi[8],phi[7]); _ff_add(f[1],t1,phi[6]);
	_ff_sum_5_mults(t1,J,J2,J3,J4,J5,phi[17],phi[16],phi[15],phi[14],phi[13]); _ff_add(f[2],t1,phi[12]);
	_ff_sum_5_mults(t1,J,J2,J3,J4,J5,phi[23],phi[22],phi[21],phi[20],phi[19]); _ff_add(f[3],t1,phi[18]);
	_ff_sum_5_mults(t1,J,J2,J3,J4,J5,phi[29],phi[28],phi[27],phi[26],phi[25]); _ff_add(f[4],t1,phi[24]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi[34],phi[33],phi[32],phi[31]); _ff_addto(t1,phi[30]); _ff_sub(f[5],t1,J5);
	_ff_set_one(f[6]);
	// 33M+31A (10 redc)
}

static inline void phi5_s2_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, J4, J5, t1;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_square(J4,J2); _ff_mult(J5,J4,J); _ff_add(t1,J2,phi[4]);
	_ff_sum_2_mults(t1,J2,J4,t1,phi[2]);  _ff_add(f[0],t1,phi[0]);
	_ff_sum_3_mults(f[1],J,J3,J5,phi[11],phi[9],phi[7]);
	_ff_sum_2_mults(t1,J2,J4,phi[16],phi[14]); _ff_add(f[2],t1,phi[12]);
	_ff_sum_3_mults(f[3],J,J3,J5,phi[23],phi[21],phi[19]);
	_ff_sum_2_mults(t1,J2,J4,phi[28],phi[26]); _ff_add(f[4],t1,phi[24]);
	_ff_sum_2_mults(t1,J,J3,phi[33],phi[31]); _ff_sub(f[5],t1,J5);
	_ff_set_one(f[6]);
	// 19M+13A (10 redc)
}

static inline void phi5_s3_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, J4, J5, t0, t1;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_square(J4,J2); _ff_mult(J5,J4,J);	
	_ff_add(t0,phi[3],J3);  _ff_mult(t1,t0,J3); _ff_add(f[0],t1,phi[0]);
	_ff_sum_2_mults(f[1],J,J4,phi[10],phi[7]);
	_ff_sum_2_mults(f[2],J2,J5,phi[17],phi[14]);
	_ff_mult(t1,J3,phi[21]); _ff_add(f[3],t1,phi[18]);
	_ff_sum_2_mults(f[4],J,J4,phi[28],phi[25]);
	_ff_mult(t1,J2,phi[32]); _ff_sub(f[5],t1,J5);
	_ff_set_one(f[6]);
	// 14M + 7A (10 redc)
}

static inline void phi5_s4_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, J4, J5, t1;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_square(J4,J2); _ff_mult(J5,J4,J);
	_ff_mult(f[0],J,J5);
	_ff_sum_2_mults(f[1],J,J5,phi[11],phi[7]);
	_ff_mult(f[2],J4,phi[16]);
	_ff_mult(f[3],J3,phi[21]);
	_ff_mult(f[4],J2,phi[26]);
	_ff_mult(t1,J,phi[31]); _ff_sub(f[5],t1,J5);	
	_ff_set_one(f[6]);
	// 11M+2A (10 redc)
}

static inline void phi5_s6_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, J4, J5;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2); _ff_square(J4,J2); _ff_mult(J5,J4,J);
	_ff_mult(f[0],J,J5);
	_ff_mult(f[1],J,phi[7]);
	_ff_mult(f[2],J2,phi[14]);
	_ff_mult(f[3],J3,phi[21]);
	_ff_mult(f[4],J4,phi[28]);
	_ff_neg(f[5],J5);		
	_ff_set_one(f[6]);
	// 8M+1A (9 redc)
}

static inline void phi5_s8_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J4, J5;
	
	_ff_square(J2,J);  _ff_square(J4,J2); _ff_mult(J5,J4,J);
	_ff_mult(f[0],J,J5);
	_ff_mult(f[1],J,phi[7]);
	_ff_mult(f[2],J4,phi[16]);
	_ff_set_zero(f[3]);
	_ff_mult(f[4],J2,phi[26]);
	_ff_neg(f[5],J5);
	_ff_set_one(f[6]);
	// 7M+2A (10 redc)
}


static inline void phi5_s12_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, J5;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2); _ff_mult(J5,J3,J2);
	_ff_mult(f[0],J,J5);
	_ff_mult(f[1],J,phi[7]);
	_ff_set_zero(f[2]);
	_ff_mult(f[3],J3,phi[21]);
	_ff_set_zero(f[4]);
	_ff_neg(f[5],J5);		
	_ff_set_one(f[6]);
	// 6M+1A (9 redc)
}

static inline void phi5_s24_eval_ff (ff_t f[], ff_t phi[], ff_t J)	
{
	ff_t J2,J4;
	
	_ff_square(J2,J); _ff_square(J4,J2);
	_ff_mult(f[0],J2,J4); _ff_mult(J2,J,J4); _ff_neg(f[5],J2);
	_ff_set_zero(f[4]); _ff_set_zero(f[3]); _ff_set_zero(f[2]);
	_ff_mult(f[1],J,phi[7]);
	_ff_set_one(f[6]);
	// 5M+1A (5 redc)
}

static inline void phi7_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6); _ff_add(t1,J,phi[7]);
    _ff_sum_7_mults(t1,J,J2,J3,J4,J5,J6,J7,t1,phi[6],phi[5],phi[4],phi[3],phi[2],phi[1]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_7_mults(t1,J,J2,J3,J4,J5,J6,J7,phi[15],phi[14],phi[13],phi[12],phi[11],phi[10],phi[9]);  _ff_add(f[1],t1,phi[8]);
    _ff_sum_7_mults(t1,J,J2,J3,J4,J5,J6,J7,phi[23],phi[22],phi[21],phi[20],phi[19],phi[18],phi[17]);  _ff_add(f[2],t1,phi[16]);
    _ff_sum_7_mults(t1,J,J2,J3,J4,J5,J6,J7,phi[31],phi[30],phi[29],phi[28],phi[27],phi[26],phi[25]);  _ff_add(f[3],t1,phi[24]);
    _ff_sum_7_mults(t1,J,J2,J3,J4,J5,J6,J7,phi[39],phi[38],phi[37],phi[36],phi[35],phi[34],phi[33]);  _ff_add(f[4],t1,phi[32]);
    _ff_sum_7_mults(t1,J,J2,J3,J4,J5,J6,J7,phi[47],phi[46],phi[45],phi[44],phi[43],phi[42],phi[41]);  _ff_add(f[5],t1,phi[40]);
    _ff_sum_7_mults(t1,J,J2,J3,J4,J5,J6,J7,phi[55],phi[54],phi[53],phi[52],phi[51],phi[50],phi[49]);  _ff_add(f[6],t1,phi[48]);
    _ff_sum_6_mults(t1,J,J2,J3,J4,J5,J6,phi[62],phi[61],phi[60],phi[59],phi[58],phi[57]);  _ff_addto(t1,phi[56]); _ff_sub(f[7],t1,J7);
    _ff_set_one(f[8]);
    // 61M+57A (14 redc)
}

static inline void phi7_s2_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6); _ff_add(t1,J2,phi[6]);
    _ff_sum_3_mults(t1,J2,J4,J6,t1,phi[4],phi[2]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_4_mults(f[1],J,J3,J5,J7,phi[15],phi[13],phi[11],phi[9]);
    _ff_sum_3_mults(t1,J2,J4,J6,phi[22],phi[20],phi[18]);  _ff_add(f[2],t1,phi[16]);
    _ff_sum_4_mults(f[3],J,J3,J5,J7,phi[31],phi[29],phi[27],phi[25]);
    _ff_sum_3_mults(t1,J2,J4,J6,phi[38],phi[36],phi[34]);  _ff_add(f[4],t1,phi[32]);
    _ff_sum_4_mults(f[5],J,J3,J5,J7,phi[47],phi[45],phi[43],phi[41]);
    _ff_sum_3_mults(t1,J2,J4,J6,phi[54],phi[52],phi[50]);  _ff_add(f[6],t1,phi[48]);
    _ff_sum_3_mults(t1,J,J3,J5,phi[61],phi[59],phi[57]);  _ff_sub(f[7],t1,J7);
    _ff_set_one(f[8]);
    // 33M+25A (14 redc)
}

static inline void phi7_s3_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6); _ff_add(t1,J3,phi[5]);
    _ff_sum_2_mults(f[0],J2,J5,t1,phi[2]);
    _ff_sum_3_mults(f[1],J,J4,J7,phi[15],phi[12],phi[9]);
    _ff_sum_2_mults(t1,J3,J6,phi[22],phi[19]);  _ff_add(f[2],t1,phi[16]);
    _ff_sum_2_mults(f[3],J2,J5,phi[29],phi[26]);
    _ff_sum_3_mults(f[4],J,J4,J7,phi[39],phi[36],phi[33]);
    _ff_sum_2_mults(t1,J3,J6,phi[46],phi[43]);  _ff_add(f[5],t1,phi[40]);
    _ff_sum_2_mults(f[6],J2,J5,phi[53],phi[50]);
    _ff_sum_2_mults(t1,J,J4,phi[60],phi[57]);  _ff_sub(f[7],t1,J7);
    _ff_set_one(f[8]);
    // 24M+15A (14 redc)
}

static inline void phi7_s4_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6); _ff_add(t1,J4,phi[4]);
    _ff_mult(t1,J4,t1);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_2_mults(f[1],J,J5,phi[13],phi[9]);
    _ff_sum_2_mults(f[2],J2,J6,phi[22],phi[18]);
    _ff_sum_2_mults(f[3],J3,J7,phi[31],phi[27]);
    _ff_mult(t1,J4,phi[36]);  _ff_add(f[4],t1,phi[32]);
    _ff_sum_2_mults(f[5],J,J5,phi[45],phi[41]);
    _ff_sum_2_mults(f[6],J2,J6,phi[54],phi[50]);
    _ff_mult(t1,J3,phi[59]);  _ff_sub(f[7],t1,J7);
    _ff_set_one(f[8]);
    // 19M+9A (14 redc)
}

static inline void phi7_s6_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);
    _ff_mult(f[0],J,J7);
    _ff_sum_2_mults(f[1],J,J7,phi[15],phi[9]);
    _ff_mult(f[2],J6,phi[22]);
    _ff_mult(f[3],J5,phi[29]);
    _ff_mult(f[4],J4,phi[36]);
    _ff_mult(f[5],J3,phi[43]);
    _ff_mult(f[6],J2,phi[50]);
    _ff_mult(t1,J,phi[57]);  _ff_sub(f[7],t1,J7);
    _ff_set_one(f[8]);
    // 15M+2A (14 redc)
}

static inline void phi7_s8_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);
    _ff_mult(t1,J,J7);  _ff_add(f[0],t1,phi[0]);
    _ff_mult(f[1],J,phi[9]);
    _ff_mult(f[2],J2,phi[18]);
    _ff_mult(f[3],J3,phi[27]);
    _ff_mult(f[4],J4,phi[36]);
    _ff_mult(f[5],J5,phi[45]);
    _ff_mult(f[6],J6,phi[54]);
    _ff_neg(f[7],J7);
    _ff_set_one(f[8]);
    // 13M+2A (14 redc)
}

static inline void phi7_s12_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J4, J6, J7;

	_ff_square(J2,J);  _ff_square(J4,J2);   _ff_mult(J6,J2,J4);  _ff_mult(J7,J,J6);
	_ff_mult(f[0],J,J7);
	_ff_mult(f[1],J,phi[9]);
	_ff_mult(f[2],J6,phi[22]);
	_ff_set_zero(f[3]);
	_ff_mult(f[4],J4,phi[36]);
	_ff_set_zero(f[5]);
	_ff_mult(f[6],J2,phi[50]);
	_ff_neg(f[7],J7);
	_ff_set_one(f[8]);
	// 9M+1A (9 redc)
}

static inline void phi7_s24_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J4, J6, J7;

	_ff_square(J2,J);  _ff_square(J4,J2);   _ff_mult(J6,J2,J4);  _ff_mult(J7,J,J6);
	_ff_mult(f[0],J,J7);
	_ff_mult(f[1],J,phi[9]);
	_ff_set_zero(f[2]);
	_ff_set_zero(f[3]);
	_ff_mult(f[4],J4,phi[36]);
	_ff_set_zero(f[5]);
	_ff_set_zero(f[6]);
	_ff_neg(f[7],J7);
	_ff_set_one(f[8]);
	// 7M+1A (7 redc)
}

static inline void phi11_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10); _ff_add(t1,J,phi[11]);
    _ff_sum_11_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,t1,phi[10],phi[9],phi[8],phi[7],phi[6],phi[5],phi[4],phi[3],phi[2],phi[1]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_11_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,phi[23],phi[22],phi[21],phi[20],phi[19],phi[18],phi[17],phi[16],phi[15],phi[14],phi[13]);  _ff_add(f[1],t1,phi[12]);
    _ff_sum_11_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,phi[35],phi[34],phi[33],phi[32],phi[31],phi[30],phi[29],phi[28],phi[27],phi[26],phi[25]);  _ff_add(f[2],t1,phi[24]);
    _ff_sum_11_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,phi[47],phi[46],phi[45],phi[44],phi[43],phi[42],phi[41],phi[40],phi[39],phi[38],phi[37]);  _ff_add(f[3],t1,phi[36]);
    _ff_sum_11_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,phi[59],phi[58],phi[57],phi[56],phi[55],phi[54],phi[53],phi[52],phi[51],phi[50],phi[49]);  _ff_add(f[4],t1,phi[48]);
    _ff_sum_11_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,phi[71],phi[70],phi[69],phi[68],phi[67],phi[66],phi[65],phi[64],phi[63],phi[62],phi[61]);  _ff_add(f[5],t1,phi[60]);
    _ff_sum_11_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,phi[83],phi[82],phi[81],phi[80],phi[79],phi[78],phi[77],phi[76],phi[75],phi[74],phi[73]);  _ff_add(f[6],t1,phi[72]);
    _ff_sum_11_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,phi[95],phi[94],phi[93],phi[92],phi[91],phi[90],phi[89],phi[88],phi[87],phi[86],phi[85]);  _ff_add(f[7],t1,phi[84]);
    _ff_sum_11_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,phi[107],phi[106],phi[105],phi[104],phi[103],phi[102],phi[101],phi[100],phi[99],phi[98],phi[97]);  _ff_add(f[8],t1,phi[96]);
    _ff_sum_11_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,phi[119],phi[118],phi[117],phi[116],phi[115],phi[114],phi[113],phi[112],phi[111],phi[110],phi[109]);  _ff_add(f[9],t1,phi[108]);
    _ff_sum_11_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,phi[131],phi[130],phi[129],phi[128],phi[127],phi[126],phi[125],phi[124],phi[123],phi[122],phi[121]);  _ff_add(f[10],t1,phi[120]);
    _ff_sum_10_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,phi[142],phi[141],phi[140],phi[139],phi[138],phi[137],phi[136],phi[135],phi[134],phi[133]);  _ff_addto(t1,phi[132]); _ff_sub(f[11],t1,J11);
    _ff_set_one(f[12]);
}

static inline void phi11_s2_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10); _ff_add(t1,J2,phi[10]);
    _ff_sum_5_mults(t1,J2,J4,J6,J8,J10,t1,phi[8],phi[6],phi[4],phi[2]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_6_mults(f[1],J,J3,J5,J7,J9,J11,phi[23],phi[21],phi[19],phi[17],phi[15],phi[13]);
    _ff_sum_5_mults(t1,J2,J4,J6,J8,J10,phi[34],phi[32],phi[30],phi[28],phi[26]);  _ff_add(f[2],t1,phi[24]);
    _ff_sum_6_mults(f[3],J,J3,J5,J7,J9,J11,phi[47],phi[45],phi[43],phi[41],phi[39],phi[37]);
    _ff_sum_5_mults(t1,J2,J4,J6,J8,J10,phi[58],phi[56],phi[54],phi[52],phi[50]);  _ff_add(f[4],t1,phi[48]);
    _ff_sum_6_mults(f[5],J,J3,J5,J7,J9,J11,phi[71],phi[69],phi[67],phi[65],phi[63],phi[61]);
    _ff_sum_5_mults(t1,J2,J4,J6,J8,J10,phi[82],phi[80],phi[78],phi[76],phi[74]);  _ff_add(f[6],t1,phi[72]);
    _ff_sum_6_mults(f[7],J,J3,J5,J7,J9,J11,phi[95],phi[93],phi[91],phi[89],phi[87],phi[85]);
    _ff_sum_5_mults(t1,J2,J4,J6,J8,J10,phi[106],phi[104],phi[102],phi[100],phi[98]);  _ff_add(f[8],t1,phi[96]);
    _ff_sum_6_mults(f[9],J,J3,J5,J7,J9,J11,phi[119],phi[117],phi[115],phi[113],phi[111],phi[109]);
    _ff_sum_5_mults(t1,J2,J4,J6,J8,J10,phi[130],phi[128],phi[126],phi[124],phi[122]);  _ff_add(f[10],t1,phi[120]);
    _ff_sum_5_mults(t1,J,J3,J5,J7,J9,phi[141],phi[139],phi[137],phi[135],phi[133]);  _ff_sub(f[11],t1,J11);
    _ff_set_one(f[12]);
}

static inline void phi11_s3_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10); _ff_add(t1,J3,phi[9]);
    _ff_sum_3_mults(t1,J3,J6,J9,t1,phi[6],phi[3]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_4_mults(f[1],J,J4,J7,J10,phi[22],phi[19],phi[16],phi[13]);
    _ff_sum_4_mults(f[2],J2,J5,J8,J11,phi[35],phi[32],phi[29],phi[26]);
    _ff_sum_3_mults(t1,J3,J6,J9,phi[45],phi[42],phi[39]);  _ff_add(f[3],t1,phi[36]);
    _ff_sum_4_mults(f[4],J,J4,J7,J10,phi[58],phi[55],phi[52],phi[49]);
    _ff_sum_4_mults(f[5],J2,J5,J8,J11,phi[71],phi[68],phi[65],phi[62]);
    _ff_sum_3_mults(t1,J3,J6,J9,phi[81],phi[78],phi[75]);  _ff_add(f[6],t1,phi[72]);
    _ff_sum_4_mults(f[7],J,J4,J7,J10,phi[94],phi[91],phi[88],phi[85]);
    _ff_sum_4_mults(f[8],J2,J5,J8,J11,phi[107],phi[104],phi[101],phi[98]);
    _ff_sum_3_mults(t1,J3,J6,J9,phi[117],phi[114],phi[111]);  _ff_add(f[9],t1,phi[108]);
    _ff_sum_4_mults(f[10],J,J4,J7,J10,phi[130],phi[127],phi[124],phi[121]);
    _ff_sum_3_mults(t1,J2,J5,J8,phi[140],phi[137],phi[134]);  _ff_sub(f[11],t1,J11);
    _ff_set_one(f[12]);
}

static inline void phi11_s4_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10); _ff_add(t1,J4,phi[8]);
    _ff_sum_2_mults(t1,J4,J8,t1,phi[4]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_3_mults(f[1],J,J5,J9,phi[21],phi[17],phi[13]);
    _ff_sum_3_mults(f[2],J2,J6,J10,phi[34],phi[30],phi[26]);
    _ff_sum_3_mults(f[3],J3,J7,J11,phi[47],phi[43],phi[39]);
    _ff_sum_2_mults(t1,J4,J8,phi[56],phi[52]);  _ff_add(f[4],t1,phi[48]);
    _ff_sum_3_mults(f[5],J,J5,J9,phi[69],phi[65],phi[61]);
    _ff_sum_3_mults(f[6],J2,J6,J10,phi[82],phi[78],phi[74]);
    _ff_sum_3_mults(f[7],J3,J7,J11,phi[95],phi[91],phi[87]);
    _ff_sum_2_mults(t1,J4,J8,phi[104],phi[100]);  _ff_add(f[8],t1,phi[96]);
    _ff_sum_3_mults(f[9],J,J5,J9,phi[117],phi[113],phi[109]);
    _ff_sum_3_mults(f[10],J2,J6,J10,phi[130],phi[126],phi[122]);
    _ff_sum_2_mults(t1,J3,J7,phi[139],phi[135]);  _ff_sub(f[11],t1,J11);
    _ff_set_one(f[12]);
}

static inline void phi11_s6_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10); _ff_add(t1,J6,phi[6]);
    _ff_mult(t1,J6,t1);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_2_mults(f[1],J,J7,phi[19],phi[13]);
    _ff_sum_2_mults(f[2],J2,J8,phi[32],phi[26]);
    _ff_sum_2_mults(f[3],J3,J9,phi[45],phi[39]);
    _ff_sum_2_mults(f[4],J4,J10,phi[58],phi[52]);
    _ff_sum_2_mults(f[5],J5,J11,phi[71],phi[65]);
    _ff_mult(t1,J6,phi[78]);  _ff_add(f[6],t1,phi[72]);
    _ff_sum_2_mults(f[7],J,J7,phi[91],phi[85]);
    _ff_sum_2_mults(f[8],J2,J8,phi[104],phi[98]);
    _ff_sum_2_mults(f[9],J3,J9,phi[117],phi[111]);
    _ff_sum_2_mults(f[10],J4,J10,phi[130],phi[124]);
    _ff_mult(t1,J5,phi[137]);  _ff_sub(f[11],t1,J11);
    _ff_set_one(f[12]);
}

static inline void phi11_s8_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10); _ff_add(t1,J8,phi[4]);
    _ff_mult(f[0],J4,t1);
    _ff_sum_2_mults(f[1],J,J9,phi[21],phi[13]);
    _ff_mult(f[2],J6,phi[30]);
    _ff_sum_2_mults(f[3],J3,J11,phi[47],phi[39]);
    _ff_mult(t1,J8,phi[56]);  _ff_add(f[4],t1,phi[48]);
    _ff_mult(f[5],J5,phi[65]);
    _ff_sum_2_mults(f[6],J2,J10,phi[82],phi[74]);
    _ff_mult(f[7],J7,phi[91]);
    _ff_mult(f[8],J4,phi[100]);
    _ff_sum_2_mults(f[9],J,J9,phi[117],phi[109]);
    _ff_mult(f[10],J6,phi[126]);
    _ff_mult(t1,J3,phi[135]);  _ff_sub(f[11],t1,J11);
    _ff_set_one(f[12]);
}

static inline void phi11_s12_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11;

	_ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10); _ff_add(t1,J3,phi[9]);
	_ff_mult(f[0],J,J11);
	_ff_mult(f[1],J,phi[13]);
	_ff_mult(f[2],J2,phi[26]);
	_ff_mult(f[3],J3,phi[39]);
	_ff_mult(f[4],J4,phi[52]);
	_ff_mult(f[5],J5,phi[65]);
	_ff_mult(f[6],J6,phi[78]);
	_ff_mult(f[7],J7,phi[91]);
	_ff_mult(f[8],J8,phi[104]);
	_ff_mult(f[9],J9,phi[117]);
	_ff_mult(f[10],J10,phi[130]);
	_ff_neg(f[11],J11);
	_ff_set_one(f[12]);
}

static inline void phi11_s24_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, J5, J7, J9, J11;

	_ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J5,J2,J3);  _ff_mult(J7,J2,J5);  _ff_mult(J9,J2,J7);  _ff_mult(J11,J2,J9);
	_ff_mult(f[0],J,J11);
	_ff_mult(f[1],J,phi[13]); _ff_set_zero(f[2]);
	_ff_mult(f[3],J3,phi[39]); _ff_set_zero(f[4]);
	_ff_mult(f[5],J5,phi[65]); _ff_set_zero(f[6]);
	_ff_mult(f[7],J7,phi[91]); _ff_set_zero(f[8]);
	_ff_mult(f[9],J9,phi[117]); _ff_set_zero(f[10]);
	_ff_neg(f[11],J11);
	_ff_set_one(f[12]);
}

static inline void phi13_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12); _ff_add(t1,J,phi[13]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,t1,phi[12],phi[11],phi[10],phi[9],phi[8],phi[7],phi[6],phi[5],phi[4],phi[3],phi[2],phi[1]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[27],phi[26],phi[25],phi[24],phi[23],phi[22],phi[21],phi[20],phi[19],phi[18],phi[17],phi[16],phi[15]);  _ff_add(f[1],t1,phi[14]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[41],phi[40],phi[39],phi[38],phi[37],phi[36],phi[35],phi[34],phi[33],phi[32],phi[31],phi[30],phi[29]);  _ff_add(f[2],t1,phi[28]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[55],phi[54],phi[53],phi[52],phi[51],phi[50],phi[49],phi[48],phi[47],phi[46],phi[45],phi[44],phi[43]);  _ff_add(f[3],t1,phi[42]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[69],phi[68],phi[67],phi[66],phi[65],phi[64],phi[63],phi[62],phi[61],phi[60],phi[59],phi[58],phi[57]);  _ff_add(f[4],t1,phi[56]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[83],phi[82],phi[81],phi[80],phi[79],phi[78],phi[77],phi[76],phi[75],phi[74],phi[73],phi[72],phi[71]);  _ff_add(f[5],t1,phi[70]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[97],phi[96],phi[95],phi[94],phi[93],phi[92],phi[91],phi[90],phi[89],phi[88],phi[87],phi[86],phi[85]);  _ff_add(f[6],t1,phi[84]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[111],phi[110],phi[109],phi[108],phi[107],phi[106],phi[105],phi[104],phi[103],phi[102],phi[101],phi[100],phi[99]);  _ff_add(f[7],t1,phi[98]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[125],phi[124],phi[123],phi[122],phi[121],phi[120],phi[119],phi[118],phi[117],phi[116],phi[115],phi[114],phi[113]);  _ff_add(f[8],t1,phi[112]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[139],phi[138],phi[137],phi[136],phi[135],phi[134],phi[133],phi[132],phi[131],phi[130],phi[129],phi[128],phi[127]);  _ff_add(f[9],t1,phi[126]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[153],phi[152],phi[151],phi[150],phi[149],phi[148],phi[147],phi[146],phi[145],phi[144],phi[143],phi[142],phi[141]);  _ff_add(f[10],t1,phi[140]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[167],phi[166],phi[165],phi[164],phi[163],phi[162],phi[161],phi[160],phi[159],phi[158],phi[157],phi[156],phi[155]);  _ff_add(f[11],t1,phi[154]);
    _ff_sum_13_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,phi[181],phi[180],phi[179],phi[178],phi[177],phi[176],phi[175],phi[174],phi[173],phi[172],phi[171],phi[170],phi[169]);  _ff_add(f[12],t1,phi[168]);
    _ff_sum_12_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,phi[194],phi[193],phi[192],phi[191],phi[190],phi[189],phi[188],phi[187],phi[186],phi[185],phi[184],phi[183]);  _ff_addto(t1,phi[182]); _ff_sub(f[13],t1,J13);
    _ff_set_one(f[14]);
}

static inline void phi13_s2_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12); _ff_add(t1,J2,phi[12]);
    _ff_sum_6_mults(t1,J2,J4,J6,J8,J10,J12,t1,phi[10],phi[8],phi[6],phi[4],phi[2]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_7_mults(f[1],J,J3,J5,J7,J9,J11,J13,phi[27],phi[25],phi[23],phi[21],phi[19],phi[17],phi[15]);
    _ff_sum_6_mults(t1,J2,J4,J6,J8,J10,J12,phi[40],phi[38],phi[36],phi[34],phi[32],phi[30]);  _ff_add(f[2],t1,phi[28]);
    _ff_sum_7_mults(f[3],J,J3,J5,J7,J9,J11,J13,phi[55],phi[53],phi[51],phi[49],phi[47],phi[45],phi[43]);
    _ff_sum_6_mults(t1,J2,J4,J6,J8,J10,J12,phi[68],phi[66],phi[64],phi[62],phi[60],phi[58]);  _ff_add(f[4],t1,phi[56]);
    _ff_sum_7_mults(f[5],J,J3,J5,J7,J9,J11,J13,phi[83],phi[81],phi[79],phi[77],phi[75],phi[73],phi[71]);
    _ff_sum_6_mults(t1,J2,J4,J6,J8,J10,J12,phi[96],phi[94],phi[92],phi[90],phi[88],phi[86]);  _ff_add(f[6],t1,phi[84]);
    _ff_sum_7_mults(f[7],J,J3,J5,J7,J9,J11,J13,phi[111],phi[109],phi[107],phi[105],phi[103],phi[101],phi[99]);
    _ff_sum_6_mults(t1,J2,J4,J6,J8,J10,J12,phi[124],phi[122],phi[120],phi[118],phi[116],phi[114]);  _ff_add(f[8],t1,phi[112]);
    _ff_sum_7_mults(f[9],J,J3,J5,J7,J9,J11,J13,phi[139],phi[137],phi[135],phi[133],phi[131],phi[129],phi[127]);
    _ff_sum_6_mults(t1,J2,J4,J6,J8,J10,J12,phi[152],phi[150],phi[148],phi[146],phi[144],phi[142]);  _ff_add(f[10],t1,phi[140]);
    _ff_sum_7_mults(f[11],J,J3,J5,J7,J9,J11,J13,phi[167],phi[165],phi[163],phi[161],phi[159],phi[157],phi[155]);
    _ff_sum_6_mults(t1,J2,J4,J6,J8,J10,J12,phi[180],phi[178],phi[176],phi[174],phi[172],phi[170]);  _ff_add(f[12],t1,phi[168]);
    _ff_sum_6_mults(t1,J,J3,J5,J7,J9,J11,phi[193],phi[191],phi[189],phi[187],phi[185],phi[183]);  _ff_sub(f[13],t1,J13);
    _ff_set_one(f[14]);
}

static inline void phi13_s3_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12); _ff_add(t1,J3,phi[11]);
    _ff_sum_4_mults(f[0],J2,J5,J8,J11,t1,phi[8],phi[5],phi[2]);
    _ff_sum_5_mults(f[1],J,J4,J7,J10,J13,phi[27],phi[24],phi[21],phi[18],phi[15]);
    _ff_sum_4_mults(t1,J3,J6,J9,J12,phi[40],phi[37],phi[34],phi[31]);  _ff_add(f[2],t1,phi[28]);
    _ff_sum_4_mults(f[3],J2,J5,J8,J11,phi[53],phi[50],phi[47],phi[44]);
    _ff_sum_5_mults(f[4],J,J4,J7,J10,J13,phi[69],phi[66],phi[63],phi[60],phi[57]);
    _ff_sum_4_mults(t1,J3,J6,J9,J12,phi[82],phi[79],phi[76],phi[73]);  _ff_add(f[5],t1,phi[70]);
    _ff_sum_4_mults(f[6],J2,J5,J8,J11,phi[95],phi[92],phi[89],phi[86]);
    _ff_sum_5_mults(f[7],J,J4,J7,J10,J13,phi[111],phi[108],phi[105],phi[102],phi[99]);
    _ff_sum_4_mults(t1,J3,J6,J9,J12,phi[124],phi[121],phi[118],phi[115]);  _ff_add(f[8],t1,phi[112]);
    _ff_sum_4_mults(f[9],J2,J5,J8,J11,phi[137],phi[134],phi[131],phi[128]);
    _ff_sum_5_mults(f[10],J,J4,J7,J10,J13,phi[153],phi[150],phi[147],phi[144],phi[141]);
    _ff_sum_4_mults(t1,J3,J6,J9,J12,phi[166],phi[163],phi[160],phi[157]);  _ff_add(f[11],t1,phi[154]);
    _ff_sum_4_mults(f[12],J2,J5,J8,J11,phi[179],phi[176],phi[173],phi[170]);
    _ff_sum_4_mults(t1,J,J4,J7,J10,phi[192],phi[189],phi[186],phi[183]);  _ff_sub(f[13],t1,J13);
    _ff_set_one(f[14]);
}

static inline void phi13_s4_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12); _ff_add(t1,J4,phi[10]);
    _ff_sum_3_mults(f[0],J2,J6,J10,t1,phi[6],phi[2]);
    _ff_sum_4_mults(f[1],J,J5,J9,J13,phi[27],phi[23],phi[19],phi[15]);
    _ff_sum_3_mults(t1,J4,J8,J12,phi[40],phi[36],phi[32]);  _ff_add(f[2],t1,phi[28]);
    _ff_sum_3_mults(f[3],J3,J7,J11,phi[53],phi[49],phi[45]);
    _ff_sum_3_mults(f[4],J2,J6,J10,phi[66],phi[62],phi[58]);
    _ff_sum_4_mults(f[5],J,J5,J9,J13,phi[83],phi[79],phi[75],phi[71]);
    _ff_sum_3_mults(t1,J4,J8,J12,phi[96],phi[92],phi[88]);  _ff_add(f[6],t1,phi[84]);
    _ff_sum_3_mults(f[7],J3,J7,J11,phi[109],phi[105],phi[101]);
    _ff_sum_3_mults(f[8],J2,J6,J10,phi[122],phi[118],phi[114]);
    _ff_sum_4_mults(f[9],J,J5,J9,J13,phi[139],phi[135],phi[131],phi[127]);
    _ff_sum_3_mults(t1,J4,J8,J12,phi[152],phi[148],phi[144]);  _ff_add(f[10],t1,phi[140]);
    _ff_sum_3_mults(f[11],J3,J7,J11,phi[165],phi[161],phi[157]);
    _ff_sum_3_mults(f[12],J2,J6,J10,phi[178],phi[174],phi[170]); 
    _ff_sum_3_mults(t1,J,J5,J9,phi[191],phi[187],phi[183]);  _ff_sub(f[13],t1,J13);
    _ff_set_one(f[14]);
}

static inline void phi13_s6_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12); _ff_add(t1,J6,phi[8]);
    _ff_sum_2_mults(f[0],J2,J8,t1,phi[2]);
    _ff_sum_3_mults(f[1],J,J7,J13,phi[27],phi[21],phi[15]);
    _ff_sum_2_mults(t1,J6,J12,phi[40],phi[34]);  _ff_add(f[2],t1,phi[28]);
    _ff_sum_2_mults(f[3],J5,J11,phi[53],phi[47]);
    _ff_sum_2_mults(f[4],J4,J10,phi[66],phi[60]);
    _ff_sum_2_mults(f[5],J3,J9,phi[79],phi[73]);
    _ff_sum_2_mults(f[6],J2,J8,phi[92],phi[86]);
    _ff_sum_3_mults(f[7],J,J7,J13,phi[111],phi[105],phi[99]);
    _ff_sum_2_mults(t1,J6,J12,phi[124],phi[118]);  _ff_add(f[8],t1,phi[112]);
    _ff_sum_2_mults(f[9],J5,J11,phi[137],phi[131]);
    _ff_sum_2_mults(f[10],J4,J10,phi[150],phi[144]);
    _ff_sum_2_mults(f[11],J3,J9,phi[163],phi[157]);
    _ff_sum_2_mults(f[12],J2,J8,phi[176],phi[170]);
    _ff_sum_2_mults(t1,J,J7,phi[189],phi[183]);  _ff_sub(f[13],t1,J13);
    _ff_set_one(f[14]);
}

static inline void phi13_s8_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12); _ff_add(t1,J8,phi[6]);
    _ff_mult(f[0],J6,t1);
    _ff_sum_2_mults(f[1],J,J9,phi[23],phi[15]);
    _ff_sum_2_mults(f[2],J4,J12,phi[40],phi[32]);
    _ff_mult(f[3],J7,phi[49]);
    _ff_sum_2_mults(f[4],J2,J10,phi[66],phi[58]);
    _ff_sum_2_mults(f[5],J5,J13,phi[83],phi[75]);
    _ff_mult(t1,J8,phi[92]);  _ff_add(f[6],t1,phi[84]);
    _ff_sum_2_mults(f[7],J3,J11,phi[109],phi[101]);
    _ff_mult(f[8],J6,phi[118]);
    _ff_sum_2_mults(f[9],J,J9,phi[135],phi[127]);
    _ff_sum_2_mults(f[10],J4,J12,phi[152],phi[144]);
    _ff_mult(f[11],J7,phi[161]);
    _ff_sum_2_mults(f[12],J2,J10,phi[178],phi[170]); 
    _ff_mult(t1,J5,phi[187]);  _ff_sub(f[13],t1,J13);
    _ff_set_one(f[14]);
}

static inline void phi13_s12_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12); _ff_add(t1,J12,phi[2]);
    _ff_mult(f[0],J2,t1);
    _ff_sum_2_mults(f[1],J,J13,phi[27],phi[15]);
    _ff_mult(t1,J12,phi[40]);  _ff_add(f[2],t1,phi[28]);
    _ff_mult(f[3],J11,phi[53]);
    _ff_mult(f[4],J10,phi[66]);
    _ff_mult(f[5],J9,phi[79]);
    _ff_mult(f[6],J8,phi[92]);
    _ff_mult(f[7],J7,phi[105]);
    _ff_mult(f[8],J6,phi[118]);
    _ff_mult(f[9],J5,phi[131]);
    _ff_mult(f[10],J4,phi[144]);
    _ff_mult(f[11],J3,phi[157]);
    _ff_mult(f[12],J2,phi[170]);
    _ff_mult(t1,J,phi[183]);  _ff_sub(f[13],t1,J13);
    _ff_set_one(f[14]);
}

static inline void phi13_s24_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J4, J6, J8, J10, J12;

	_ff_square(J2,J);  _ff_square(J4,J2);  _ff_mult(J6,J2,J4);  _ff_mult(J8,J2,J6);  _ff_mult(J10,J2,J8);  _ff_mult(J12,J2,J10);
	_ff_mult(f[0],J2,J12);
	_ff_mult(f[1],J,phi[15]);
	_ff_mult(f[2],J12,phi[40]); _ff_set_zero(f[3]);
	_ff_mult(f[4],J10,phi[66]); _ff_set_zero(f[5]);
	_ff_mult(f[6],J8,phi[92]); _ff_set_zero(f[7]);
	_ff_mult(f[8],J6,phi[118]); _ff_set_zero(f[9]);
	_ff_mult(f[10],J4,phi[144]); _ff_set_zero(f[11]);
	_ff_mult(f[12],J2,phi[170]);
	_ff_mult(J2,J12,J); _ff_neg(f[13],J2);
	_ff_set_one(f[14]);
}

static inline void phi17_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16);   _ff_add(t1,J,phi[17]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,t1,phi[16],phi[15],phi[14],phi[13],phi[12],phi[11],phi[10],phi[9],phi[8],phi[7],phi[6],phi[5],phi[4],phi[3],phi[2],phi[1]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[35],phi[34],phi[33],phi[32],phi[31],phi[30],phi[29],phi[28],phi[27],phi[26],phi[25],phi[24],phi[23],phi[22],phi[21],phi[20],phi[19]);  _ff_add(f[1],t1,phi[18]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[53],phi[52],phi[51],phi[50],phi[49],phi[48],phi[47],phi[46],phi[45],phi[44],phi[43],phi[42],phi[41],phi[40],phi[39],phi[38],phi[37]);  _ff_add(f[2],t1,phi[36]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[71],phi[70],phi[69],phi[68],phi[67],phi[66],phi[65],phi[64],phi[63],phi[62],phi[61],phi[60],phi[59],phi[58],phi[57],phi[56],phi[55]);  _ff_add(f[3],t1,phi[54]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[89],phi[88],phi[87],phi[86],phi[85],phi[84],phi[83],phi[82],phi[81],phi[80],phi[79],phi[78],phi[77],phi[76],phi[75],phi[74],phi[73]);  _ff_add(f[4],t1,phi[72]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[107],phi[106],phi[105],phi[104],phi[103],phi[102],phi[101],phi[100],phi[99],phi[98],phi[97],phi[96],phi[95],phi[94],phi[93],phi[92],phi[91]);  _ff_add(f[5],t1,phi[90]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[125],phi[124],phi[123],phi[122],phi[121],phi[120],phi[119],phi[118],phi[117],phi[116],phi[115],phi[114],phi[113],phi[112],phi[111],phi[110],phi[109]);  _ff_add(f[6],t1,phi[108]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[143],phi[142],phi[141],phi[140],phi[139],phi[138],phi[137],phi[136],phi[135],phi[134],phi[133],phi[132],phi[131],phi[130],phi[129],phi[128],phi[127]);  _ff_add(f[7],t1,phi[126]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[161],phi[160],phi[159],phi[158],phi[157],phi[156],phi[155],phi[154],phi[153],phi[152],phi[151],phi[150],phi[149],phi[148],phi[147],phi[146],phi[145]);  _ff_add(f[8],t1,phi[144]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[179],phi[178],phi[177],phi[176],phi[175],phi[174],phi[173],phi[172],phi[171],phi[170],phi[169],phi[168],phi[167],phi[166],phi[165],phi[164],phi[163]);  _ff_add(f[9],t1,phi[162]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[197],phi[196],phi[195],phi[194],phi[193],phi[192],phi[191],phi[190],phi[189],phi[188],phi[187],phi[186],phi[185],phi[184],phi[183],phi[182],phi[181]);  _ff_add(f[10],t1,phi[180]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[215],phi[214],phi[213],phi[212],phi[211],phi[210],phi[209],phi[208],phi[207],phi[206],phi[205],phi[204],phi[203],phi[202],phi[201],phi[200],phi[199]);  _ff_add(f[11],t1,phi[198]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[233],phi[232],phi[231],phi[230],phi[229],phi[228],phi[227],phi[226],phi[225],phi[224],phi[223],phi[222],phi[221],phi[220],phi[219],phi[218],phi[217]);  _ff_add(f[12],t1,phi[216]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[251],phi[250],phi[249],phi[248],phi[247],phi[246],phi[245],phi[244],phi[243],phi[242],phi[241],phi[240],phi[239],phi[238],phi[237],phi[236],phi[235]);  _ff_add(f[13],t1,phi[234]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[269],phi[268],phi[267],phi[266],phi[265],phi[264],phi[263],phi[262],phi[261],phi[260],phi[259],phi[258],phi[257],phi[256],phi[255],phi[254],phi[253]);  _ff_add(f[14],t1,phi[252]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[287],phi[286],phi[285],phi[284],phi[283],phi[282],phi[281],phi[280],phi[279],phi[278],phi[277],phi[276],phi[275],phi[274],phi[273],phi[272],phi[271]);  _ff_add(f[15],t1,phi[270]);
    _ff_sum_17_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,phi[305],phi[304],phi[303],phi[302],phi[301],phi[300],phi[299],phi[298],phi[297],phi[296],phi[295],phi[294],phi[293],phi[292],phi[291],phi[290],phi[289]);  _ff_add(f[16],t1,phi[288]);
    _ff_sum_16_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,phi[322],phi[321],phi[320],phi[319],phi[318],phi[317],phi[316],phi[315],phi[314],phi[313],phi[312],phi[311],phi[310],phi[309],phi[308],phi[307]);  _ff_addto(t1,phi[306]); _ff_sub(f[17],t1,J17);
    _ff_set_one(f[18]);
}

static inline void phi17_s2_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16);   _ff_add(t1,J2,phi[16]);
    _ff_sum_8_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,t1,phi[14],phi[12],phi[10],phi[8],phi[6],phi[4],phi[2]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_9_mults(f[1],J,J3,J5,J7,J9,J11,J13,J15,J17,phi[35],phi[33],phi[31],phi[29],phi[27],phi[25],phi[23],phi[21],phi[19]);
    _ff_sum_8_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,phi[52],phi[50],phi[48],phi[46],phi[44],phi[42],phi[40],phi[38]);  _ff_add(f[2],t1,phi[36]);
    _ff_sum_9_mults(f[3],J,J3,J5,J7,J9,J11,J13,J15,J17,phi[71],phi[69],phi[67],phi[65],phi[63],phi[61],phi[59],phi[57],phi[55]);
    _ff_sum_8_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,phi[88],phi[86],phi[84],phi[82],phi[80],phi[78],phi[76],phi[74]);  _ff_add(f[4],t1,phi[72]);
    _ff_sum_9_mults(f[5],J,J3,J5,J7,J9,J11,J13,J15,J17,phi[107],phi[105],phi[103],phi[101],phi[99],phi[97],phi[95],phi[93],phi[91]);
    _ff_sum_8_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,phi[124],phi[122],phi[120],phi[118],phi[116],phi[114],phi[112],phi[110]);  _ff_add(f[6],t1,phi[108]);
    _ff_sum_9_mults(f[7],J,J3,J5,J7,J9,J11,J13,J15,J17,phi[143],phi[141],phi[139],phi[137],phi[135],phi[133],phi[131],phi[129],phi[127]);
    _ff_sum_8_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,phi[160],phi[158],phi[156],phi[154],phi[152],phi[150],phi[148],phi[146]);  _ff_add(f[8],t1,phi[144]);
    _ff_sum_9_mults(f[9],J,J3,J5,J7,J9,J11,J13,J15,J17,phi[179],phi[177],phi[175],phi[173],phi[171],phi[169],phi[167],phi[165],phi[163]); 
    _ff_sum_8_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,phi[196],phi[194],phi[192],phi[190],phi[188],phi[186],phi[184],phi[182]);  _ff_add(f[10],t1,phi[180]);
    _ff_sum_9_mults(f[11],J,J3,J5,J7,J9,J11,J13,J15,J17,phi[215],phi[213],phi[211],phi[209],phi[207],phi[205],phi[203],phi[201],phi[199]);
    _ff_sum_8_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,phi[232],phi[230],phi[228],phi[226],phi[224],phi[222],phi[220],phi[218]);  _ff_add(f[12],t1,phi[216]);
    _ff_sum_9_mults(f[13],J,J3,J5,J7,J9,J11,J13,J15,J17,phi[251],phi[249],phi[247],phi[245],phi[243],phi[241],phi[239],phi[237],phi[235]);
    _ff_sum_8_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,phi[268],phi[266],phi[264],phi[262],phi[260],phi[258],phi[256],phi[254]);  _ff_add(f[14],t1,phi[252]);
    _ff_sum_9_mults(f[15],J,J3,J5,J7,J9,J11,J13,J15,J17,phi[287],phi[285],phi[283],phi[281],phi[279],phi[277],phi[275],phi[273],phi[271]);
    _ff_sum_8_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,phi[304],phi[302],phi[300],phi[298],phi[296],phi[294],phi[292],phi[290]);  _ff_add(f[16],t1,phi[288]);
    _ff_sum_8_mults(t1,J,J3,J5,J7,J9,J11,J13,J15,phi[321],phi[319],phi[317],phi[315],phi[313],phi[311],phi[309],phi[307]); _ff_sub(f[17],t1,J17);
    _ff_set_one(f[18]);
}

static inline void phi17_s3_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16);   _ff_add(t1,J3,phi[15]);
    _ff_sum_5_mults(t1,J3,J6,J9,J12,J15,t1,phi[12],phi[9],phi[6],phi[3]); _ff_add(f[0],t1,phi[0]);
    _ff_sum_6_mults(f[1],J,J4,J7,J10,J13,J16,phi[34],phi[31],phi[28],phi[25],phi[22],phi[19]);
    _ff_sum_6_mults(f[2],J2,J5,J8,J11,J14,J17,phi[53],phi[50],phi[47],phi[44],phi[41],phi[38]); 
    _ff_sum_5_mults(t1,J3,J6,J9,J12,J15,phi[69],phi[66],phi[63],phi[60],phi[57]); _ff_add(f[3],t1,phi[54]);
    _ff_sum_6_mults(f[4],J,J4,J7,J10,J13,J16,phi[88],phi[85],phi[82],phi[79],phi[76],phi[73]); 
    _ff_sum_6_mults(f[5],J2,J5,J8,J11,J14,J17,phi[107],phi[104],phi[101],phi[98],phi[95],phi[92]);
    _ff_sum_5_mults(t1,J3,J6,J9,J12,J15,phi[123],phi[120],phi[117],phi[114],phi[111]);  _ff_add(f[6],t1,phi[108]);
    _ff_sum_6_mults(f[7],J,J4,J7,J10,J13,J16,phi[142],phi[139],phi[136],phi[133],phi[130],phi[127]);
    _ff_sum_6_mults(f[8],J2,J5,J8,J11,J14,J17,phi[161],phi[158],phi[155],phi[152],phi[149],phi[146]);
    _ff_sum_5_mults(t1,J3,J6,J9,J12,J15,phi[177],phi[174],phi[171],phi[168],phi[165]);  _ff_add(f[9],t1,phi[162]);
    _ff_sum_6_mults(f[10],J,J4,J7,J10,J13,J16,phi[196],phi[193],phi[190],phi[187],phi[184],phi[181]);
    _ff_sum_6_mults(f[11],J2,J5,J8,J11,J14,J17,phi[215],phi[212],phi[209],phi[206],phi[203],phi[200]);
    _ff_sum_5_mults(t1,J3,J6,J9,J12,J15,phi[231],phi[228],phi[225],phi[222],phi[219]); _ff_add(f[12],t1,phi[216]);
    _ff_sum_6_mults(f[13],J,J4,J7,J10,J13,J16,phi[250],phi[247],phi[244],phi[241],phi[238],phi[235]);
    _ff_sum_6_mults(f[14],J2,J5,J8,J11,J14,J17,phi[269],phi[266],phi[263],phi[260],phi[257],phi[254]);
    _ff_sum_5_mults(t1,J3,J6,J9,J12,J15,phi[285],phi[282],phi[279],phi[276],phi[273]); _ff_add(f[15],t1,phi[270]);
    _ff_sum_6_mults(f[16],J,J4,J7,J10,J13,J16,phi[304],phi[301],phi[298],phi[295],phi[292],phi[289]);
    _ff_sum_5_mults(t1,J2,J5,J8,J11,J14,phi[320],phi[317],phi[314],phi[311],phi[308]); _ff_sub(f[17],t1,J17);
    _ff_set_one(f[18]);
}

static inline void phi17_s4_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16);   _ff_add(t1,J4,phi[14]);
    _ff_sum_4_mults(f[0],J2,J6,J10,J14,t1,phi[10],phi[6],phi[2]);
    _ff_sum_5_mults(f[1],J,J5,J9,J13,J17,phi[35],phi[31],phi[27],phi[23],phi[19]);
    _ff_sum_4_mults(t1,J4,J8,J12,J16,phi[52],phi[48],phi[44],phi[40]);  _ff_add(f[2],t1,phi[36]);
    _ff_sum_4_mults(f[3],J3,J7,J11,J15,phi[69],phi[65],phi[61],phi[57]);
    _ff_sum_4_mults(f[4],J2,J6,J10,J14,phi[86],phi[82],phi[78],phi[74]); 
    _ff_sum_5_mults(f[5],J,J5,J9,J13,J17,phi[107],phi[103],phi[99],phi[95],phi[91]);
    _ff_sum_4_mults(t1,J4,J8,J12,J16,phi[124],phi[120],phi[116],phi[112]);  _ff_add(f[6],t1,phi[108]);
    _ff_sum_4_mults(f[7],J3,J7,J11,J15,phi[141],phi[137],phi[133],phi[129]);
    _ff_sum_4_mults(f[8],J2,J6,J10,J14,phi[158],phi[154],phi[150],phi[146]);
    _ff_sum_5_mults(f[9],J,J5,J9,J13,J17,phi[179],phi[175],phi[171],phi[167],phi[163]); 
    _ff_sum_4_mults(t1,J4,J8,J12,J16,phi[196],phi[192],phi[188],phi[184]);  _ff_add(f[10],t1,phi[180]);
    _ff_sum_4_mults(f[11],J3,J7,J11,J15,phi[213],phi[209],phi[205],phi[201]);
    _ff_sum_4_mults(f[12],J2,J6,J10,J14,phi[230],phi[226],phi[222],phi[218]);
    _ff_sum_5_mults(f[13],J,J5,J9,J13,J17,phi[251],phi[247],phi[243],phi[239],phi[235]);
    _ff_sum_4_mults(t1,J4,J8,J12,J16,phi[268],phi[264],phi[260],phi[256]);  _ff_add(f[14],t1,phi[252]);
    _ff_sum_4_mults(f[15],J3,J7,J11,J15,phi[285],phi[281],phi[277],phi[273]);
    _ff_sum_4_mults(f[16],J2,J6,J10,J14,phi[302],phi[298],phi[294],phi[290]);
    _ff_sum_4_mults(t1,J,J5,J9,J13,phi[319],phi[315],phi[311],phi[307]); _ff_sub(f[17],t1,J17);
    _ff_set_one(f[18]);
}

static inline void phi17_s6_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16);   _ff_add(t1,J6,phi[12]);
    _ff_sum_2_mults(t1,J6,J12,t1,phi[6]); _ff_add(f[0],t1,phi[0]);
    _ff_sum_3_mults(f[1],J,J7,J13,phi[31],phi[25],phi[19]);
    _ff_sum_3_mults(f[2],J2,J8,J14,phi[50],phi[44],phi[38]); 
    _ff_sum_3_mults(f[3],J3,J9,J15,phi[69],phi[63],phi[57]);
    _ff_sum_3_mults(f[4],J4,J10,J16,phi[88],phi[82],phi[76]); 
    _ff_sum_3_mults(f[5],J5,J11,J17,phi[107],phi[101],phi[95]);
    _ff_sum_2_mults(t1,J6,J12,phi[120],phi[114]);  _ff_add(f[6],t1,phi[108]);
    _ff_sum_3_mults(f[7],J,J7,J13,phi[139],phi[133],phi[127]);
    _ff_sum_3_mults(f[8],J2,J8,J14,phi[158],phi[152],phi[146]);
    _ff_sum_3_mults(f[9],J3,J9,J15,phi[177],phi[171],phi[165]);
    _ff_sum_3_mults(f[10],J4,J10,J16,phi[196],phi[190],phi[184]);
    _ff_sum_3_mults(f[11],J5,J11,J17,phi[215],phi[209],phi[203]);
    _ff_sum_2_mults(t1,J6,J12,phi[228],phi[222]); _ff_add(f[12],t1,phi[216]);
    _ff_sum_3_mults(f[13],J,J7,J13,phi[247],phi[241],phi[235]);
    _ff_sum_3_mults(f[14],J2,J8,J14,phi[266],phi[260],phi[254]);
    _ff_sum_3_mults(f[15],J3,J9,J15,phi[285],phi[279],phi[273]);
    _ff_sum_3_mults(f[16],J4,J10,J16,phi[304],phi[298],phi[292]);
    _ff_sum_2_mults(t1,J5,J11,phi[317],phi[311]); _ff_sub(f[17],t1,J17);
    _ff_set_one(f[18]);
}

static inline void phi17_s8_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16);   _ff_add(t1,J8,phi[10]);
    _ff_sum_2_mults(f[0],J2,J10,t1,phi[2]);
    _ff_sum_3_mults(f[1],J,J9,J17,phi[35],phi[27],phi[19]);
    _ff_sum_2_mults(t1,J8,J16,phi[52],phi[44]);  _ff_add(f[2],t1,phi[36]);
    _ff_sum_2_mults(f[3],J7,J15,phi[69],phi[61]);
    _ff_sum_2_mults(f[4],J6,J14,phi[86],phi[78]); 
    _ff_sum_2_mults(f[5],J5,J13,phi[103],phi[95]);
    _ff_sum_2_mults(f[6],J4,J12,phi[120],phi[112]);
    _ff_sum_2_mults(f[7],J3,J11,phi[137],phi[129]);
    _ff_sum_2_mults(f[8],J2,J10,phi[154],phi[146]);
    _ff_sum_3_mults(f[9],J,J9,J17,phi[179],phi[171],phi[163]); 
    _ff_sum_2_mults(t1,J8,J16,phi[196],phi[188]);  _ff_add(f[10],t1,phi[180]);
    _ff_sum_2_mults(f[11],J7,J15,phi[213],phi[205]);
    _ff_sum_2_mults(f[12],J6,J14,phi[230],phi[222]);
    _ff_sum_2_mults(f[13],J5,J13,phi[247],phi[239]);
    _ff_sum_2_mults(f[14],J4,J12,phi[264],phi[256]);
    _ff_sum_2_mults(f[15],J3,J11,phi[281],phi[273]);
    _ff_sum_2_mults(f[16],J2,J10,phi[298],phi[290]);
    _ff_sum_2_mults(t1,J,J9,phi[315],phi[307]); _ff_sub(f[17],t1,J17);
    _ff_set_one(f[18]);
}

static inline void phi17_s12_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16);   _ff_add(t1,J12,phi[6]);
    _ff_mult(f[0],J6,t1);
    _ff_sum_2_mults(f[1],J,J13,phi[31],phi[19]);
    _ff_mult(f[2],J8,phi[44]); 
    _ff_sum_2_mults(f[3],J3,J15,phi[69],phi[57]);
    _ff_mult(f[4],J10,phi[82]); 
    _ff_sum_2_mults(f[5],J5,J17,phi[107],phi[95]);
    _ff_mult(t1,J12,phi[120]);  _ff_add(f[6],t1,phi[108]);
    _ff_mult(f[7],J7,phi[133]);
    _ff_sum_2_mults(f[8],J2,J14,phi[158],phi[146]);
    _ff_mult(f[9],J9,phi[171]);
    _ff_sum_2_mults(f[10],J4,J16,phi[196],phi[184]);
    _ff_mult(f[11],J11,phi[209]);
    _ff_mult(f[12],J6,phi[222]);
    _ff_sum_2_mults(f[13],J,J13,phi[247],phi[235]);
    _ff_mult(f[14],J8,phi[260]);
    _ff_sum_2_mults(f[15],J3,J15,phi[285],phi[273]);
    _ff_mult(f[16],J10,phi[298]);
    _ff_mult(t1,J5,phi[311]); _ff_sub(f[17],t1,J17);
    _ff_set_one(f[18]);
}


static inline void phi17_s24_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, J4, J5, J6, J8, J9, J10, J12, J13, J15, J16;

	_ff_square(J2,J);  _ff_mult(J3,J2,J); _ff_square(J4,J2);  _ff_mult(J5,J4,J); _ff_mult(J6,J5,J); _ff_mult(J8,J2,J6); _ff_mult(J9,J8,J);  _ff_mult(J10,J9,J);  _ff_mult(J12,J2,J10); _ff_mult(J13,J12,J); _ff_mult(J15,J2,J13); _ff_mult(J16,J15,J);
	_ff_mult(f[0],J2,J16);
	_ff_mult(f[1],J,phi[19]);
	_ff_mult(f[2],J8,phi[44]);
	_ff_mult(f[3],J15,phi[69]);  _ff_set_zero(f[4]);
	_ff_mult(f[5],J5,phi[95]);
	_ff_mult(f[6],J12,phi[120]); _ff_set_zero(f[7]);
	_ff_mult(f[8],J2,phi[146]);
	_ff_mult(f[9],J9,phi[171]);
	_ff_mult(f[10],J16,phi[196]); _ff_set_zero(f[11]);
	_ff_mult(f[12],J6,phi[222]);
	_ff_mult(f[13],J13,phi[247]); _ff_set_zero(f[14]);
	_ff_mult(f[15],J3,phi[273]);
	_ff_mult(f[16],J10,phi[298]);
	_ff_mult(J2,J16,J); _ff_neg(f[17],J2);
	_ff_set_one(f[18]);
}

static inline void phi19_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17, J18, J19;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);  _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16);  _ff_mult(J18,J,J17);  _ff_mult(J19,J,J18);  _ff_add(t1,J,phi[19]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,t1,phi[18],phi[17],phi[16],phi[15],phi[14],phi[13],phi[12],phi[11],phi[10],phi[9],phi[8],phi[7],phi[6],phi[5],phi[4],phi[3],phi[2],phi[1]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[39],phi[38],phi[37],phi[36],phi[35],phi[34],phi[33],phi[32],phi[31],phi[30],phi[29],phi[28],phi[27],phi[26],phi[25],phi[24],phi[23],phi[22],phi[21]);  _ff_add(f[1],t1,phi[20]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[59],phi[58],phi[57],phi[56],phi[55],phi[54],phi[53],phi[52],phi[51],phi[50],phi[49],phi[48],phi[47],phi[46],phi[45],phi[44],phi[43],phi[42],phi[41]);  _ff_add(f[2],t1,phi[40]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[79],phi[78],phi[77],phi[76],phi[75],phi[74],phi[73],phi[72],phi[71],phi[70],phi[69],phi[68],phi[67],phi[66],phi[65],phi[64],phi[63],phi[62],phi[61]);  _ff_add(f[3],t1,phi[60]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[99],phi[98],phi[97],phi[96],phi[95],phi[94],phi[93],phi[92],phi[91],phi[90],phi[89],phi[88],phi[87],phi[86],phi[85],phi[84],phi[83],phi[82],phi[81]);  _ff_add(f[4],t1,phi[80]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[119],phi[118],phi[117],phi[116],phi[115],phi[114],phi[113],phi[112],phi[111],phi[110],phi[109],phi[108],phi[107],phi[106],phi[105],phi[104],phi[103],phi[102],phi[101]);  _ff_add(f[5],t1,phi[100]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[139],phi[138],phi[137],phi[136],phi[135],phi[134],phi[133],phi[132],phi[131],phi[130],phi[129],phi[128],phi[127],phi[126],phi[125],phi[124],phi[123],phi[122],phi[121]);  _ff_add(f[6],t1,phi[120]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[159],phi[158],phi[157],phi[156],phi[155],phi[154],phi[153],phi[152],phi[151],phi[150],phi[149],phi[148],phi[147],phi[146],phi[145],phi[144],phi[143],phi[142],phi[141]);  _ff_add(f[7],t1,phi[140]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[179],phi[178],phi[177],phi[176],phi[175],phi[174],phi[173],phi[172],phi[171],phi[170],phi[169],phi[168],phi[167],phi[166],phi[165],phi[164],phi[163],phi[162],phi[161]);  _ff_add(f[8],t1,phi[160]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[199],phi[198],phi[197],phi[196],phi[195],phi[194],phi[193],phi[192],phi[191],phi[190],phi[189],phi[188],phi[187],phi[186],phi[185],phi[184],phi[183],phi[182],phi[181]);  _ff_add(f[9],t1,phi[180]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[219],phi[218],phi[217],phi[216],phi[215],phi[214],phi[213],phi[212],phi[211],phi[210],phi[209],phi[208],phi[207],phi[206],phi[205],phi[204],phi[203],phi[202],phi[201]);  _ff_add(f[10],t1,phi[200]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[239],phi[238],phi[237],phi[236],phi[235],phi[234],phi[233],phi[232],phi[231],phi[230],phi[229],phi[228],phi[227],phi[226],phi[225],phi[224],phi[223],phi[222],phi[221]);  _ff_add(f[11],t1,phi[220]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[259],phi[258],phi[257],phi[256],phi[255],phi[254],phi[253],phi[252],phi[251],phi[250],phi[249],phi[248],phi[247],phi[246],phi[245],phi[244],phi[243],phi[242],phi[241]);  _ff_add(f[12],t1,phi[240]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[279],phi[278],phi[277],phi[276],phi[275],phi[274],phi[273],phi[272],phi[271],phi[270],phi[269],phi[268],phi[267],phi[266],phi[265],phi[264],phi[263],phi[262],phi[261]);  _ff_add(f[13],t1,phi[260]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[299],phi[298],phi[297],phi[296],phi[295],phi[294],phi[293],phi[292],phi[291],phi[290],phi[289],phi[288],phi[287],phi[286],phi[285],phi[284],phi[283],phi[282],phi[281]);  _ff_add(f[14],t1,phi[280]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[319],phi[318],phi[317],phi[316],phi[315],phi[314],phi[313],phi[312],phi[311],phi[310],phi[309],phi[308],phi[307],phi[306],phi[305],phi[304],phi[303],phi[302],phi[301]);  _ff_add(f[15],t1,phi[300]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[339],phi[338],phi[337],phi[336],phi[335],phi[334],phi[333],phi[332],phi[331],phi[330],phi[329],phi[328],phi[327],phi[326],phi[325],phi[324],phi[323],phi[322],phi[321]);  _ff_add(f[16],t1,phi[320]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[359],phi[358],phi[357],phi[356],phi[355],phi[354],phi[353],phi[352],phi[351],phi[350],phi[349],phi[348],phi[347],phi[346],phi[345],phi[344],phi[343],phi[342],phi[341]);  _ff_add(f[17],t1,phi[340]);
    _ff_sum_19_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,phi[379],phi[378],phi[377],phi[376],phi[375],phi[374],phi[373],phi[372],phi[371],phi[370],phi[369],phi[368],phi[367],phi[366],phi[365],phi[364],phi[363],phi[362],phi[361]);  _ff_add(f[18],t1,phi[360]);
    _ff_sum_18_mults(t1,J,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,phi[398],phi[397],phi[396],phi[395],phi[394],phi[393],phi[392],phi[391],phi[390],phi[389],phi[388],phi[387],phi[386],phi[385],phi[384],phi[383],phi[382],phi[381]);  _ff_addto(t1,phi[380]); _ff_sub(f[19],t1,J19);
    _ff_set_one(f[20]);
}

static inline void phi19_s2_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17, J18, J19;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16);  _ff_mult(J18,J,J17);  _ff_mult(J19,J,J18);  _ff_add(t1,J2,phi[18]);
    _ff_sum_9_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,J18,t1,phi[16],phi[14],phi[12],phi[10],phi[8],phi[6],phi[4],phi[2]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_10_mults(f[1],J,J3,J5,J7,J9,J11,J13,J15,J17,J19,phi[39],phi[37],phi[35],phi[33],phi[31],phi[29],phi[27],phi[25],phi[23],phi[21]);
    _ff_sum_9_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,J18,phi[58],phi[56],phi[54],phi[52],phi[50],phi[48],phi[46],phi[44],phi[42]);  _ff_add(f[2],t1,phi[40]);
    _ff_sum_10_mults(f[3],J,J3,J5,J7,J9,J11,J13,J15,J17,J19,phi[79],phi[77],phi[75],phi[73],phi[71],phi[69],phi[67],phi[65],phi[63],phi[61]);
    _ff_sum_9_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,J18,phi[98],phi[96],phi[94],phi[92],phi[90],phi[88],phi[86],phi[84],phi[82]);  _ff_add(f[4],t1,phi[80]);
    _ff_sum_10_mults(f[5],J,J3,J5,J7,J9,J11,J13,J15,J17,J19,phi[119],phi[117],phi[115],phi[113],phi[111],phi[109],phi[107],phi[105],phi[103],phi[101]);
    _ff_sum_9_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,J18,phi[138],phi[136],phi[134],phi[132],phi[130],phi[128],phi[126],phi[124],phi[122]);  _ff_add(f[6],t1,phi[120]);
    _ff_sum_10_mults(f[7],J,J3,J5,J7,J9,J11,J13,J15,J17,J19,phi[159],phi[157],phi[155],phi[153],phi[151],phi[149],phi[147],phi[145],phi[143],phi[141]);
    _ff_sum_9_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,J18,phi[178],phi[176],phi[174],phi[172],phi[170],phi[168],phi[166],phi[164],phi[162]);  _ff_add(f[8],t1,phi[160]);
    _ff_sum_10_mults(f[9],J,J3,J5,J7,J9,J11,J13,J15,J17,J19,phi[199],phi[197],phi[195],phi[193],phi[191],phi[189],phi[187],phi[185],phi[183],phi[181]); 
    _ff_sum_9_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,J18,phi[218],phi[216],phi[214],phi[212],phi[210],phi[208],phi[206],phi[204],phi[202]);  _ff_add(f[10],t1,phi[200]);
    _ff_sum_10_mults(f[11],J,J3,J5,J7,J9,J11,J13,J15,J17,J19,phi[239],phi[237],phi[235],phi[233],phi[231],phi[229],phi[227],phi[225],phi[223],phi[221]);
    _ff_sum_9_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,J18,phi[258],phi[256],phi[254],phi[252],phi[250],phi[248],phi[246],phi[244],phi[242]);  _ff_add(f[12],t1,phi[240]);
    _ff_sum_10_mults(f[13],J,J3,J5,J7,J9,J11,J13,J15,J17,J19,phi[279],phi[277],phi[275],phi[273],phi[271],phi[269],phi[267],phi[265],phi[263],phi[261]);
    _ff_sum_9_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,J18,phi[298],phi[296],phi[294],phi[292],phi[290],phi[288],phi[286],phi[284],phi[282]);  _ff_add(f[14],t1,phi[280]);
    _ff_sum_10_mults(f[15],J,J3,J5,J7,J9,J11,J13,J15,J17,J19,phi[319],phi[317],phi[315],phi[313],phi[311],phi[309],phi[307],phi[305],phi[303],phi[301]);
    _ff_sum_9_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,J18,phi[338],phi[336],phi[334],phi[332],phi[330],phi[328],phi[326],phi[324],phi[322]);  _ff_add(f[16],t1,phi[320]);
    _ff_sum_10_mults(f[17],J,J3,J5,J7,J9,J11,J13,J15,J17,J19,phi[359],phi[357],phi[355],phi[353],phi[351],phi[349],phi[347],phi[345],phi[343],phi[341]);
    _ff_sum_9_mults(t1,J2,J4,J6,J8,J10,J12,J14,J16,J18,phi[378],phi[376],phi[374],phi[372],phi[370],phi[368],phi[366],phi[364],phi[362]);  _ff_add(f[18],t1,phi[360]);
    _ff_sum_9_mults(t1,J,J3,J5,J7,J9,J11,J13,J15,J17,phi[397],phi[395],phi[393],phi[391],phi[389],phi[387],phi[385],phi[383],phi[381]); _ff_sub(f[19],t1,J19);
    _ff_set_one(f[20]);
}

static inline void phi19_s3_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17,J18,J19;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16);  _ff_mult(J18,J,J17);  _ff_mult(J19,J,J18); _ff_add(t1,J3,phi[17]);
    _ff_sum_6_mults(f[0],J2,J5,J8,J11,J14,J17,t1,phi[14],phi[11],phi[8],phi[5],phi[2]);
    _ff_sum_7_mults(f[1],J,J4,J7,J10,J13,J16,J19,phi[39],phi[36],phi[33],phi[30],phi[27],phi[24],phi[21]);
    _ff_sum_6_mults(t1,J3,J6,J9,J12,J15,J18,phi[58],phi[55],phi[52],phi[49],phi[46],phi[43]); _ff_add(f[2],t1,phi[40]);
    _ff_sum_6_mults(f[3],J2,J5,J8,J11,J14,J17,phi[77],phi[74],phi[71],phi[68],phi[65],phi[62]);
    _ff_sum_7_mults(f[4],J,J4,J7,J10,J13,J16,J19,phi[99],phi[96],phi[93],phi[90],phi[87],phi[84],phi[81]);
    _ff_sum_6_mults(t1,J3,J6,J9,J12,J15,J18,phi[118],phi[115],phi[112],phi[109],phi[106],phi[103]); _ff_add(f[5],t1,phi[100]);
    _ff_sum_6_mults(f[6],J2,J5,J8,J11,J14,J17,phi[137],phi[134],phi[131],phi[128],phi[125],phi[122]);
    _ff_sum_7_mults(f[7],J,J4,J7,J10,J13,J16,J19,phi[159],phi[156],phi[153],phi[150],phi[147],phi[144],phi[141]);
    _ff_sum_6_mults(t1,J3,J6,J9,J12,J15,J18,phi[178],phi[175],phi[172],phi[169],phi[166],phi[163]); _ff_add(f[8],t1,phi[160]);
    _ff_sum_6_mults(f[9],J2,J5,J8,J11,J14,J17,phi[197],phi[194],phi[191],phi[188],phi[185],phi[182]);
    _ff_sum_7_mults(f[10],J,J4,J7,J10,J13,J16,J19,phi[219],phi[216],phi[213],phi[210],phi[207],phi[204],phi[201]);
    _ff_sum_6_mults(t1,J3,J6,J9,J12,J15,J18,phi[238],phi[235],phi[232],phi[229],phi[226],phi[223]); _ff_add(f[11],t1,phi[220]);
    _ff_sum_6_mults(f[12],J2,J5,J8,J11,J14,J17,phi[257],phi[254],phi[251],phi[248],phi[245],phi[242]);
    _ff_sum_7_mults(f[13],J,J4,J7,J10,J13,J16,J19,phi[279],phi[276],phi[273],phi[270],phi[267],phi[264],phi[261]);
    _ff_sum_6_mults(t1,J3,J6,J9,J12,J15,J18,phi[298],phi[295],phi[292],phi[289],phi[286],phi[283]); _ff_add(f[14],t1,phi[280]);
    _ff_sum_6_mults(f[15],J2,J5,J8,J11,J14,J17,phi[317],phi[314],phi[311],phi[308],phi[305],phi[302]);
    _ff_sum_7_mults(f[16],J,J4,J7,J10,J13,J16,J19,phi[339],phi[336],phi[333],phi[330],phi[327],phi[324],phi[321]);
    _ff_sum_6_mults(t1,J3,J6,J9,J12,J15,J18,phi[358],phi[355],phi[352],phi[349],phi[346],phi[343]); _ff_add(f[17],t1,phi[340]);
    _ff_sum_6_mults(f[18],J2,J5,J8,J11,J14,J17,phi[377],phi[374],phi[371],phi[368],phi[365],phi[362]);
    _ff_sum_6_mults(t1,J,J4,J7,J10,J13,J16,phi[396],phi[393],phi[390],phi[387],phi[384],phi[381]); _ff_sub(f[19],t1,J19);
    _ff_set_one(f[20]);
}

static inline void phi19_s4_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17, J18, J19;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16); _ff_mult(J18,J,J17);  _ff_mult(J19,J,J18);  _ff_add(t1,J4,phi[16]);
    _ff_sum_4_mults(t1,J4,J8,J12,J16,t1,phi[12],phi[8],phi[4]);  _ff_add(f[0],t1,phi[0]);
    _ff_sum_5_mults(f[1],J,J5,J9,J13,J17,phi[37],phi[33],phi[29],phi[25],phi[21]);
    _ff_sum_5_mults(f[2],J2,J6,J10,J14,J18,phi[58],phi[54],phi[50],phi[46],phi[42]);
    _ff_sum_5_mults(f[3],J3,J7,J11,J15,J19,phi[79],phi[75],phi[71],phi[67],phi[63]);
    _ff_sum_4_mults(t1,J4,J8,J12,J16,phi[96],phi[92],phi[88],phi[84]);  _ff_add(f[4],t1,phi[80]);
    _ff_sum_5_mults(f[5],J,J5,J9,J13,J17,phi[117],phi[113],phi[109],phi[105],phi[101]);
    _ff_sum_5_mults(f[6],J2,J6,J10,J14,J18,phi[138],phi[134],phi[130],phi[126],phi[122]);
    _ff_sum_5_mults(f[7],J3,J7,J11,J15,J19,phi[159],phi[155],phi[151],phi[147],phi[143]);
    _ff_sum_4_mults(t1,J4,J8,J12,J16,phi[176],phi[172],phi[168],phi[164]);  _ff_add(f[8],t1,phi[160]);
    _ff_sum_5_mults(f[9],J,J5,J9,J13,J17,phi[197],phi[193],phi[189],phi[185],phi[181]); 
    _ff_sum_5_mults(f[10],J2,J6,J10,J14,J18,phi[218],phi[214],phi[210],phi[206],phi[202]);
    _ff_sum_5_mults(f[11],J3,J7,J11,J15,J19,phi[239],phi[235],phi[231],phi[227],phi[223]);
    _ff_sum_4_mults(t1,J4,J8,J12,J16,phi[256],phi[252],phi[248],phi[244]);  _ff_add(f[12],t1,phi[240]);
    _ff_sum_5_mults(f[13],J,J5,J9,J13,J17,phi[277],phi[273],phi[269],phi[265],phi[261]);
    _ff_sum_5_mults(f[14],J2,J6,J10,J14,J18,phi[298],phi[294],phi[290],phi[286],phi[282]);
    _ff_sum_5_mults(f[15],J3,J7,J11,J15,J19,phi[319],phi[315],phi[311],phi[307],phi[303]);
    _ff_sum_4_mults(t1,J4,J8,J12,J16,phi[336],phi[332],phi[328],phi[324]);  _ff_add(f[16],t1,phi[320]);
    _ff_sum_5_mults(f[17],J,J5,J9,J13,J17,phi[357],phi[353],phi[349],phi[345],phi[341]);
    _ff_sum_5_mults(f[18],J2,J6,J10,J14,J18,phi[378],phi[374],phi[370],phi[366],phi[362]);
    _ff_sum_4_mults(t1,J3,J7,J11,J15,phi[395],phi[391],phi[387],phi[383]); _ff_sub(f[19],t1,J19);
    _ff_set_one(f[20]);
}

static inline void phi19_s6_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17,J18,J19;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16); _ff_mult(J18,J,J17);  _ff_mult(J19,J,J18);  _ff_add(t1,J6,phi[14]);
    _ff_sum_3_mults(f[0],J2,J8,J14,t1,phi[8],phi[2]);
    _ff_sum_4_mults(f[1],J,J7,J13,J19,phi[39],phi[33],phi[27],phi[21]);
    _ff_sum_3_mults(t1,J6,J12,J18,phi[58],phi[52],phi[46]); _ff_add(f[2],t1,phi[40]);
    _ff_sum_3_mults(f[3],J5,J11,J17,phi[77],phi[71],phi[65]);
    _ff_sum_3_mults(f[4],J4,J10,J16,phi[96],phi[90],phi[84]);
    _ff_sum_3_mults(f[5],J3,J9,J15,phi[115],phi[109],phi[103]);
    _ff_sum_3_mults(f[6],J2,J8,J14,phi[134],phi[128],phi[122]);
    _ff_sum_4_mults(f[7],J,J7,J13,J19,phi[159],phi[153],phi[147],phi[141]);
    _ff_sum_3_mults(t1,J6,J12,J18,phi[178],phi[172],phi[166]); _ff_add(f[8],t1,phi[160]);
    _ff_sum_3_mults(f[9],J5,J11,J17,phi[197],phi[191],phi[185]);
    _ff_sum_3_mults(f[10],J4,J10,J16,phi[216],phi[210],phi[204]);
    _ff_sum_3_mults(f[11],J3,J9,J15,phi[235],phi[229],phi[223]);
    _ff_sum_3_mults(f[12],J2,J8,J14,phi[254],phi[248],phi[242]);
    _ff_sum_4_mults(f[13],J,J7,J13,J19,phi[279],phi[273],phi[267],phi[261]);
    _ff_sum_3_mults(t1,J6,J12,J18,phi[298],phi[292],phi[286]); _ff_add(f[14],t1,phi[280]);
    _ff_sum_3_mults(f[15],J5,J11,J17,phi[317],phi[311],phi[305]);
    _ff_sum_3_mults(f[16],J4,J10,J16,phi[336],phi[330],phi[324]);
    _ff_sum_3_mults(f[17],J3,J9,J15,phi[355],phi[349],phi[343]);
    _ff_sum_3_mults(f[18],J2,J8,J14,phi[374],phi[368],phi[362]);
    _ff_sum_3_mults(t1,J,J7,J13,phi[393],phi[387],phi[381]); _ff_sub(f[19],t1,J19);
    _ff_set_one(f[20]);
}

static inline void phi19_s8_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17, J18, J19;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16);  _ff_mult(J18,J,J17);  _ff_mult(J19,J,J18); _ff_add(t1,J8,phi[12]);
    _ff_sum_2_mults(f[0],J4,J12,t1,phi[4]);
    _ff_sum_3_mults(f[1],J,J9,J17,phi[37],phi[29],phi[21]);
    _ff_sum_2_mults(f[2],J6,J14,phi[54],phi[46]);
    _ff_sum_3_mults(f[3],J3,J11,J19,phi[79],phi[71],phi[63]);
    _ff_sum_2_mults(t1,J8,J16,phi[96],phi[88]);  _ff_add(f[4],t1,phi[80]);
    _ff_sum_2_mults(f[5],J5,J13,phi[113],phi[105]);
    _ff_sum_3_mults(f[6],J2,J10,J18,phi[138],phi[130],phi[122]);
    _ff_sum_2_mults(f[7],J7,J15,phi[155],phi[147]);
    _ff_sum_2_mults(f[8],J4,J12,phi[172],phi[164]);
    _ff_sum_3_mults(f[9],J,J9,J17,phi[197],phi[189],phi[181]); 
    _ff_sum_2_mults(f[10],J6,J14,phi[214],phi[206]);
    _ff_sum_3_mults(f[11],J3,J11,J19,phi[239],phi[231],phi[223]);
    _ff_sum_2_mults(t1,J8,J16,phi[256],phi[248]);  _ff_add(f[12],t1,phi[240]);
    _ff_sum_2_mults(f[13],J5,J13,phi[273],phi[265]);
    _ff_sum_3_mults(f[14],J2,J10,J18,phi[298],phi[290],phi[282]);
    _ff_sum_2_mults(f[15],J7,J15,phi[315],phi[307]);
    _ff_sum_2_mults(f[16],J4,J12,phi[332],phi[324]);
    _ff_sum_3_mults(f[17],J,J9,J17,phi[357],phi[349],phi[341]);
    _ff_sum_2_mults(f[18],J6,J14,phi[374],phi[366]);
    _ff_sum_2_mults(t1,J3,J11,phi[391],phi[383]); _ff_sub(f[19],t1,J19);
    _ff_set_one(f[20]);
}

static inline void phi19_s12_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
    register ff_t t1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14, J15, J16, J17,J18,J19;

    _ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_mult(J4,J,J3);  _ff_mult(J5,J,J4);  _ff_mult(J6,J,J5);  _ff_mult(J7,J,J6);  _ff_mult(J8,J,J7);  _ff_mult(J9,J,J8);  _ff_mult(J10,J,J9);
   _ff_mult(J11,J,J10);  _ff_mult(J12,J,J11);  _ff_mult(J13,J,J12);  _ff_mult(J14,J,J13);  _ff_mult(J15,J,J14);  _ff_mult(J16,J,J15);  _ff_mult(J17,J,J16); _ff_mult(J18,J,J17);  _ff_mult(J19,J,J18);  _ff_add(t1,J12,phi[8]);
    _ff_mult(f[0],J8,t1);
    _ff_sum_2_mults(f[1],J,J13,phi[33],phi[21]);
    _ff_sum_2_mults(f[2],J6,J18,phi[58],phi[46]);
    _ff_mult(f[3],J11,phi[71]);
    _ff_sum_2_mults(f[4],J4,J16,phi[96],phi[84]);
    _ff_mult(f[5],J9,phi[109]);
    _ff_sum_2_mults(f[6],J2,J14,phi[134],phi[122]);
    _ff_sum_2_mults(f[7],J7,J19,phi[159],phi[147]);
    _ff_mult(t1,J12,phi[172]); _ff_add(f[8],t1,phi[160]);
    _ff_sum_2_mults(f[9],J5,J17,phi[197],phi[185]);
    _ff_mult(f[10],J10,phi[210]);
    _ff_sum_2_mults(f[11],J3,J15,phi[235],phi[223]);
    _ff_mult(f[12],J8,phi[248]);
    _ff_sum_2_mults(f[13],J,J13,phi[273],phi[261]);
    _ff_sum_2_mults(f[14],J6,J18,phi[298],phi[286]);
    _ff_mult(f[15],J11,phi[311]);
    _ff_sum_2_mults(f[16],J4,J16,phi[336],phi[324]);
    _ff_mult(f[17],J9,phi[349]);
    _ff_sum_2_mults(f[18],J2,J14,phi[374],phi[362]);
    _ff_mult(t1,J7,phi[387]); _ff_sub(f[19],t1,J19);
    _ff_set_one(f[20]);
}

static inline void phi19_s24_eval_ff (ff_t f[], ff_t phi[], ff_t J)
{
	register ff_t J2, J3, J4, J6, J7, J8, J9, J11, J12, J13, J14, J16, J17, J18;

	_ff_square(J2,J);  _ff_mult(J3,J2,J); _ff_square(J4,J2); _ff_mult(J6,J4,J2); _ff_mult(J7,J6,J); _ff_mult(J8,J2,J6); _ff_mult(J9,J8,J); 
	_ff_mult(J11,J2,J9); _ff_mult(J12,J11,J); _ff_mult(J13,J12,J); _ff_mult(J14,J13,J); _ff_mult(J16,J14,J2); _ff_mult(J17,J16,J); _ff_mult(J18,J17,J);
	_ff_mult(f[0],J2,J18);
	_ff_mult(f[1],J,phi[21]);
	_ff_mult(f[2],J6,phi[46]);
	_ff_mult(f[3],J11,phi[71]);
	_ff_mult(f[4],J16,phi[96]);  _ff_set_zero(f[5]);
	_ff_mult(f[6],J2,phi[122]);
	_ff_mult(f[7],J7,phi[147]);
	_ff_mult(f[8],J12,phi[172]);
	_ff_mult(f[9],J17,phi[197]);  _ff_set_zero(f[10]);
	_ff_mult(f[11],J3,phi[223]);
	_ff_mult(f[12],J8,phi[248]);
	_ff_mult(f[13],J13,phi[273]);
	_ff_mult(f[14],J18,phi[298]);  _ff_set_zero(f[15]);
	_ff_mult(f[16],J4,phi[324]);
	_ff_mult(f[17],J9,phi[349]);
	_ff_mult(f[18],J14,phi[374]);
	_ff_mult(J2,J18,J); _ff_neg(f[19],J2);
	_ff_set_one(f[20]);
}

static inline void phi2_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t J0,J1,J2,t1;
	
	_ff_square(t1,b); _ff_mult(J1,t1,a); _ff_mult(J0,t1,b);
	_ff_square(t1,a); _ff_mult(J2,t1,b);
	_ff_sum_4_mults(f[0],J0,J1,J2,t1,a,phi[2],phi[1],phi[0]);
	_ff_sum_3_mults(f[1],J0,J1,J2,phi[5],phi[4],phi[3]);
	_ff_sum_2_mults(t1,J0,J1,phi[7],phi[6]);  _ff_sub(f[2],t1,J2);
	_ff_set(f[3],J0);
	// 14M+7A (8 redc)
}

static inline void phi2_s3_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t t0,t1,t2,t3;
	
	_ff_square(t0,a); _ff_square(t1,b); _ff_mult(t2,t1,b);
	_ff_sum_2_mults(f[0],t2,t0,a,phi[0]);
	_ff_mult(t3,a,t1); _ff_mult(f[1],t3,phi[4]);
	_ff_mult(t3,b,t0); _ff_neg(f[2],t3);
	_ff_set(f[3],t2);
}

static inline void phi3_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t J0,J1,J2,J3,t0,t1,t2;
	
	_ff_square(t0,a);  _ff_square(t1,b); _ff_mult(t2,a,b);
	_ff_square(J0,t1); _ff_mult(J1,t1,t2); _ff_mult(J2,t0,t1);  _ff_mult(J3,t0,t2);
	_ff_sum_5_mults(f[0],J0,J1,J2,J3,t0,t0,phi[3],phi[2],phi[1],phi[0]);
	_ff_sum_4_mults(f[1],J0,J1,J2,J3,phi[7],phi[6],phi[5],phi[4]);
	_ff_sum_4_mults(f[2],J0,J1,J2,J3,phi[11],phi[10],phi[9],phi[8]);
	_ff_sum_3_mults(t1,J0,J1,J2,phi[14],phi[13],phi[12]); _ff_sub(f[3],t1,J3);
	_ff_set(f[4],J0);
}

static inline void phi3_s2_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t J0,J1,J2,J3,t0,t1,t2;
	
	_ff_square(t0,a);  _ff_square(t1,b); _ff_mult(t2,a,b);
	_ff_square(J0,t1); _ff_mult(J1,t1,t2); _ff_mult(J2,t0,t1);  _ff_mult(J3,t0,t2);
	_ff_sum_3_mults(f[0],J0,J2,t0,t0,phi[2],phi[0]);
	_ff_sum_2_mults(f[1],J1,J3,phi[7],phi[5]);
	_ff_sum_2_mults(f[2],J0,J2,phi[10],phi[8]);
	_ff_mult(t1,J1,phi[13]);  _ff_sub(f[3],t1,J3);
	_ff_set(f[4],J0);
}

static inline void phi3_s4_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t J0,J1,J2,J3,t0,t1,t2;
	
	_ff_square(t0,a);  _ff_square(t1,b); _ff_mult(t2,a,b);
	_ff_square(J0,t1); _ff_mult(J1,t1,t2); _ff_mult(J2,t0,t1);  _ff_mult(J3,t0,t2);
	_ff_sum_2_mults(f[0],J0,t0,t0,phi[0]);
	_ff_mult(f[1],J1,phi[5]);
	_ff_mult(f[2],J2,phi[10]);
	_ff_neg(f[3],J3);
	_ff_set(f[4],J0);
}

static inline void phi3_s8_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t J0,J1,J3,t0,t1,t2;
	
	_ff_square(t0,a);  _ff_square(t1,b); _ff_mult(t2,a,b);
	_ff_square(J0,t1); _ff_mult(J1,t1,t2); _ff_mult(J3,t0,t2);
	_ff_mult(f[0],t0,t0);
	_ff_mult(f[1],J1,phi[5]);
	_ff_set_zero(f[2]);
	_ff_neg(f[3],J3);
	_ff_set(f[4],J0);
}

static inline void phi5_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t J0,J1,J2,J3,J4,J5,t0,t1,t2,t3,t4;
	
	_ff_square(t0,a);  _ff_square(t1,b); _ff_mult(t2,a,b); _ff_square(t3,t0); _ff_square(t4,t1);
	_ff_mult(J0,t1,t4); _ff_mult(J1,t2,t4); _ff_mult(J2,t0,t4); _ff_square(t4,t2); _ff_mult(J3,t4,t2); _ff_mult(J4,t1,t3); _ff_mult(J5,t2,t3);
	_ff_sum_7_mults(f[0],J0,J1,J2,J3,J4,J5,t0,t3,phi[5],phi[4],phi[3],phi[2],phi[1],phi[0]);
	_ff_sum_6_mults(f[1],J0,J1,J2,J3,J4,J5,phi[11],phi[10],phi[9],phi[8],phi[7],phi[6]);
	_ff_sum_6_mults(f[2],J0,J1,J2,J3,J4,J5,phi[17],phi[16],phi[15],phi[14],phi[13],phi[12]);
	_ff_sum_6_mults(f[3],J0,J1,J2,J3,J4,J5,phi[23],phi[22],phi[21],phi[20],phi[19],phi[18]);
	_ff_sum_6_mults(f[4],J0,J1,J2,J3,J4,J5,phi[29],phi[28],phi[27],phi[26],phi[25],phi[24]);
	_ff_sum_5_mults(t1,J0,J1,J2,J3,J4,phi[34],phi[33],phi[32],phi[31],phi[30]); _ff_sub(f[5],t1,J5);
	_ff_set(f[6],J0);
}

static inline void phi5_s2_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t J0,J1,J2,J3,J4,J5,t0,t1,t2,t3,t4;
	
	_ff_square(t0,a);  _ff_square(t1,b); _ff_mult(t2,a,b); _ff_square(t3,t0); _ff_square(t4,t1);
	_ff_mult(J0,t1,t4); _ff_mult(J1,t2,t4); _ff_mult(J2,t0,t4); _ff_square(t4,t2); _ff_mult(J3,t4,t2); _ff_mult(J4,t1,t3); _ff_mult(J5,t2,t3);
	_ff_sum_4_mults(f[0],J0,J2,J4,t0,t3,phi[4],phi[2],phi[0]);
	_ff_sum_3_mults(f[1],J1,J3,J5,phi[11],phi[9],phi[7]);
	_ff_sum_3_mults(f[2],J0,J2,J4,phi[16],phi[14],phi[12]);
	_ff_sum_3_mults(f[3],J1,J3,J5,phi[23],phi[21],phi[19]);
	_ff_sum_3_mults(f[4],J0,J2,J4,phi[28],phi[26],phi[24]);
	_ff_sum_2_mults(t1,J1,J3,phi[33],phi[31]); _ff_sub(f[5],t1,J5);
	_ff_set(f[6],J0);
}

static inline void phi5_s3_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t J0,J1,J2,J3,J4,J5,t0,t1,t2,t3,t4;
	
	_ff_square(t0,a);  _ff_square(t1,b); _ff_mult(t2,a,b); _ff_square(t3,t0); _ff_square(t4,t1);
	_ff_mult(J0,t1,t4); _ff_mult(J1,t2,t4); _ff_mult(J2,t0,t4); _ff_square(t4,t2); _ff_mult(J3,t4,t2); _ff_mult(J4,t1,t3); _ff_mult(J5,t2,t3);
	_ff_sum_3_mults(f[0],J0,J3,t0,t3,phi[3],phi[0]);
	_ff_sum_2_mults(f[1],J1,J4,phi[10],phi[7]);
	_ff_sum_2_mults(f[2],J2,J5,phi[17],phi[14]);
	_ff_sum_2_mults(f[3],J0,J3,phi[21],phi[18]);
	_ff_sum_2_mults(f[4],J1,J4,phi[28],phi[25]);
	_ff_mult(t1,J2,phi[32]); _ff_sub(f[5],t1,J5);
	_ff_set(f[6],J0);
}

static inline void phi5_s4_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t J0,J1,J2,J3,J4,J5,t0,t1,t2,t3,t4;
	
	_ff_square(t0,a);  _ff_square(t1,b); _ff_mult(t2,a,b); _ff_square(t3,t0); _ff_square(t4,t1);
	_ff_mult(J0,t1,t4); _ff_mult(J1,t2,t4); _ff_mult(J2,t0,t4); _ff_square(t4,t2); _ff_mult(J3,t4,t2); _ff_mult(J4,t1,t3); _ff_mult(J5,t2,t3);
	_ff_sum_2_mults(f[0],J2,t0,t3,phi[2]);
	_ff_sum_2_mults(f[1],J1,J5,phi[11],phi[7]);
	_ff_sum_2_mults(f[2],J0,J4,phi[16],phi[12]);
	_ff_mult(f[3],J3,phi[21]);
	_ff_mult(f[4],J2,phi[26]);
	_ff_mult(t1,J1,phi[31]); _ff_sub(f[5],t1,J5);
	_ff_set(f[6],J0);
}

static inline void phi5_s6_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t J0,J1,J2,J3,J4,J5,t0,t1,t2,t3,t4;
	
	_ff_square(t0,a);  _ff_square(t1,b); _ff_mult(t2,a,b); _ff_square(t3,t0); _ff_square(t4,t1);
	_ff_mult(J0,t1,t4); _ff_mult(J1,t2,t4); _ff_mult(J2,t0,t4); _ff_square(t4,t2); _ff_mult(J3,t4,t2); _ff_mult(J4,t1,t3); _ff_mult(J5,t2,t3);
	_ff_sum_2_mults(f[0],J0,t0,t3,phi[0]);
	_ff_mult(f[1],J1,phi[7]);
	_ff_mult(f[2],J2,phi[14]);
	_ff_mult(f[3],J3,phi[21]);
	_ff_mult(f[4],J4,phi[28]);
	_ff_neg(f[5],J5);
	_ff_set(f[6],J0);
}

static inline void phi5_s8_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t J0,J1,J2,J4,J5,t0,t1,t2,t3,t4;
	
	_ff_square(t0,a);  _ff_square(t1,b); _ff_mult(t2,a,b); _ff_square(t3,t0); _ff_square(t4,t1);
	_ff_mult(J0,t1,t4); _ff_mult(J1,t2,t4); _ff_mult(J2,t0,t4); _ff_mult(J4,t1,t3); _ff_mult(J5,t2,t3);
	_ff_mult(f[0],t0,t3);
	_ff_mult(f[1],J1,phi[7]);
	_ff_mult(f[2],J4,phi[16]);
	_ff_set_zero(f[3]);
	_ff_mult(f[4],J2,phi[26]);
	_ff_neg(f[5],J5);
	_ff_set(f[6],J0);
}

static inline void phi5_s12_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t t0,t1,t2,t3,t4,w;
	
	_ff_square(t0,a);  _ff_square(t1,t0); _ff_square(t2,b); _ff_square(t3,t2);  _ff_mult(t4,a,b);
	_ff_mult(f[6],t2,t3); _ff_mult(w,t1,t4); _ff_neg(f[5],w);
	_ff_set_zero(f[4]);
	_ff_mult(w,t4,t0); ff_mult(w,w,t2); _ff_mult(f[3],w,phi[21]);
	_ff_set_zero(f[2]);
	_ff_mult(w,t3,t4); _ff_mult(f[1],w,phi[7]);
	_ff_mult(f[0],t0,t1);
	// 13M + 1A (13redc)
}

static inline void phi5_s24_qeval_ff (ff_t f[], ff_t phi[], ff_t a, ff_t b)
{
	register ff_t t0,t1,t2,t3,t4,w;
	
	_ff_square(t0,a);  _ff_square(t1,t0); _ff_square(t2,b); _ff_square(t3,t2);  _ff_mult(t4,a,b);
	_ff_mult(f[6],t2,t3); _ff_mult(w,t1,t4); _ff_neg(f[5],w);
	_ff_set_zero(f[4]); _ff_set_zero(f[3]); _ff_set_zero(f[2]);
	_ff_mult(w,t3,t4); _ff_mult(f[1],w,phi[7]);
	_ff_mult(f[0],t0,t1);
	// 10M + 1A (9redc)
}

#endif
