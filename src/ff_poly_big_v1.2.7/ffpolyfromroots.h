#ifndef _FFPOLYFROMROOTS_H_
#define _FFPOLYFROMROOTS_H_

/*
    Copyright 2009-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include "ff.h"

#ifdef __cplusplus
extern "C" {
#endif

void ff_poly_from_roots_tiny (ff_t o[], ff_t r[], int d);
void ff_poly_from_roots_small (ff_t o[], ff_t r[], int d);
void ff_poly_from_roots_med (ff_t o[], ff_t r[], int d);
void ff_poly_from_roots_64 (ff_t o[], ff_t r[]);
void ff_poly_from_roots_128 (ff_t o[], ff_t r[]);
void ff_poly_from_roots_256 (ff_t o[], ff_t r[]);

static inline void ff_poly_from_roots_2 (ff_t o[], ff_t r[])
{
	ff_t t0;
	
	_ff_add(t0,r[0],r[1]); _ff_neg(o[1],t0);
	_ff_mult(o[0],r[0],r[1]);
}

static inline void ff_poly_from_roots_3 (ff_t o[], ff_t r[])
{
	register ff_t t0, t1, t2;
	
	_ff_mult(t0,r[0],r[1]); _ff_add(t1,r[0],r[1]); _ff_add(t2,t1,r[2]); _ff_neg(o[2],t2);
	_ff_mult(t2,t1,r[2]); _ff_add(o[1],t2,t0);
	_ff_mult(t2,t0,r[2]); _ff_neg(o[0],t2);
}

static inline void ff_poly_from_roots_4 (ff_t o[], ff_t r[])
{
	register ff_t f0,f1,g0,g1,t0,t1;
	
	_ff_add(t0,r[0],r[1]); _ff_neg(f1,t0);  _ff_mult(f0,r[0],r[1]);
	_ff_add(t0,r[2],r[3]); _ff_neg(g1,t0);  _ff_mult(g0,r[2],r[3]);
	_ff_add(o[3],f1,g1);
	_ff_mult(t0,f1,g1); _ff_add(t1,f0,g0); _ff_add(o[2],t0,t1);
	_ff_sum_2_mults(o[1],f0,f1,g0,g1);
	_ff_mult(o[0],f0,g0);
}

static inline void ff_poly_from_roots_5 (ff_t o[5], ff_t r[5])
{
   ff_t f[3], g[2];
    register ff_t t0, t1;

    ff_poly_from_roots_3(f,r);  ff_poly_from_roots_2(g,r+3);     _ff_add(o[4],f[2],g[1]);
    _ff_mult(t0,f[2],g[1]); _ff_add(t1,f[1],g[0]); _ff_add(o[3],t0,t1);
    _ff_sum_2_mults_arr(t0,f+1,g);  _ff_add(o[2],t0,f[0]);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_6 (ff_t o[6], ff_t r[6])
{
   ff_t f[3], g[3];
    register ff_t t0, t1;

    ff_poly_from_roots_3(f,r);  ff_poly_from_roots_3(g,r+3);     _ff_add(o[5],f[2],g[2]);
    _ff_mult(t0,f[2],g[2]); _ff_add(t1,f[1],g[1]); _ff_add(o[4],t0,t1);
    _ff_sum_2_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[3],t0,t1);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_7 (ff_t o[7], ff_t r[7])
{
   ff_t f[4], g[3];
    register ff_t t0, t1;

    ff_poly_from_roots_4(f,r);  ff_poly_from_roots_3(g,r+4);     _ff_add(o[6],f[3],g[2]);
    _ff_mult(t0,f[3],g[2]); _ff_add(t1,f[2],g[1]); _ff_add(o[5],t0,t1);
    _ff_sum_2_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[4],t0,t1);
    _ff_sum_3_mults_arr(t0,f+1,g);  _ff_add(o[3],t0,f[0]);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_8 (ff_t o[8], ff_t r[8])
{
   ff_t f[4], g[4];
    register ff_t t0, t1;

    ff_poly_from_roots_4(f,r);  ff_poly_from_roots_4(g,r+4);     _ff_add(o[7],f[3],g[3]);
    _ff_mult(t0,f[3],g[3]); _ff_add(t1,f[2],g[2]); _ff_add(o[6],t0,t1);
    _ff_sum_2_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[5],t0,t1);
    _ff_sum_3_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[4],t0,t1);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_9 (ff_t o[9], ff_t r[9])
{
   ff_t f[5], g[4];
    register ff_t t0, t1;

    ff_poly_from_roots_5(f,r);  ff_poly_from_roots_4(g,r+5);     _ff_add(o[8],f[4],g[3]);
    _ff_mult(t0,f[4],g[3]); _ff_add(t1,f[3],g[2]); _ff_add(o[7],t0,t1);
    _ff_sum_2_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[6],t0,t1);
    _ff_sum_3_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[5],t0,t1);
    _ff_sum_4_mults_arr(t0,f+1,g);  _ff_add(o[4],t0,f[0]);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_10 (ff_t o[10], ff_t r[10])
{
   ff_t f[5], g[5];
    register ff_t t0, t1;

    ff_poly_from_roots_5(f,r);  ff_poly_from_roots_5(g,r+5);     _ff_add(o[9],f[4],g[4]);
    _ff_mult(t0,f[4],g[4]); _ff_add(t1,f[3],g[3]); _ff_add(o[8],t0,t1);
    _ff_sum_2_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[7],t0,t1);
    _ff_sum_3_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[6],t0,t1);
    _ff_sum_4_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[5],t0,t1);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_11 (ff_t o[11], ff_t r[11])
{
   ff_t f[6], g[5];
    register ff_t t0, t1;

    ff_poly_from_roots_6(f,r);  ff_poly_from_roots_5(g,r+6);     _ff_add(o[10],f[5],g[4]);
    _ff_mult(t0,f[5],g[4]); _ff_add(t1,f[4],g[3]); _ff_add(o[9],t0,t1);
    _ff_sum_2_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[8],t0,t1);
    _ff_sum_3_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[7],t0,t1);
    _ff_sum_4_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[6],t0,t1);
    _ff_sum_5_mults_arr(t0,f+1,g);  _ff_add(o[5],t0,f[0]);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_12 (ff_t o[12], ff_t r[12])
{
   ff_t f[6], g[6];
    register ff_t t0, t1;

    ff_poly_from_roots_6(f,r);  ff_poly_from_roots_6(g,r+6);     _ff_add(o[11],f[5],g[5]);
    _ff_mult(t0,f[5],g[5]); _ff_add(t1,f[4],g[4]); _ff_add(o[10],t0,t1);
    _ff_sum_2_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[9],t0,t1);
    _ff_sum_3_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[8],t0,t1);
    _ff_sum_4_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[7],t0,t1);
    _ff_sum_5_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[6],t0,t1);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_13 (ff_t o[13], ff_t r[13])
{
   ff_t f[7], g[6];
    register ff_t t0, t1;

    ff_poly_from_roots_7(f,r);  ff_poly_from_roots_6(g,r+7);     _ff_add(o[12],f[6],g[5]);
    _ff_mult(t0,f[6],g[5]); _ff_add(t1,f[5],g[4]); _ff_add(o[11],t0,t1);
    _ff_sum_2_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[10],t0,t1);
    _ff_sum_3_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[9],t0,t1);
    _ff_sum_4_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[8],t0,t1);
    _ff_sum_5_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[7],t0,t1);
    _ff_sum_6_mults_arr(t0,f+1,g);  _ff_add(o[6],t0,f[0]);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_14 (ff_t o[14], ff_t r[14])
{
   ff_t f[7], g[7];
    register ff_t t0, t1;

    ff_poly_from_roots_7(f,r);  ff_poly_from_roots_7(g,r+7);     _ff_add(o[13],f[6],g[6]);
    _ff_mult(t0,f[6],g[6]); _ff_add(t1,f[5],g[5]); _ff_add(o[12],t0,t1);
    _ff_sum_2_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[11],t0,t1);
    _ff_sum_3_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[10],t0,t1);
    _ff_sum_4_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[9],t0,t1);
    _ff_sum_5_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[8],t0,t1);
    _ff_sum_6_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[7],t0,t1);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_15 (ff_t o[15], ff_t r[15])
{
   ff_t f[8], g[7];
    register ff_t t0, t1;

    ff_poly_from_roots_8(f,r);  ff_poly_from_roots_7(g,r+8);     _ff_add(o[14],f[7],g[6]);
    _ff_mult(t0,f[7],g[6]); _ff_add(t1,f[6],g[5]); _ff_add(o[13],t0,t1);
    _ff_sum_2_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[12],t0,t1);
    _ff_sum_3_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[11],t0,t1);
    _ff_sum_4_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[10],t0,t1);
    _ff_sum_5_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[9],t0,t1);
    _ff_sum_6_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[8],t0,t1);
    _ff_sum_7_mults_arr(t0,f+1,g);  _ff_add(o[7],t0,f[0]);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_16 (ff_t o[16], ff_t r[16])
{
   ff_t f[8], g[8];
    register ff_t t0, t1;

    ff_poly_from_roots_8(f,r);  ff_poly_from_roots_8(g,r+8);     _ff_add(o[15],f[7],g[7]);
    _ff_mult(t0,f[7],g[7]); _ff_add(t1,f[6],g[6]); _ff_add(o[14],t0,t1);
    _ff_sum_2_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[13],t0,t1);
    _ff_sum_3_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[12],t0,t1);
    _ff_sum_4_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[11],t0,t1);
    _ff_sum_5_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[10],t0,t1);
    _ff_sum_6_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[9],t0,t1);
    _ff_sum_7_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[8],t0,t1);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_17 (ff_t o[17], ff_t r[17])
{
   ff_t f[9], g[8];
    register ff_t t0, t1;

    ff_poly_from_roots_9(f,r);  ff_poly_from_roots_8(g,r+9);     _ff_add(o[16],f[8],g[7]);
    _ff_mult(t0,f[8],g[7]); _ff_add(t1,f[7],g[6]); _ff_add(o[15],t0,t1);
    _ff_sum_2_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[14],t0,t1);
    _ff_sum_3_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[13],t0,t1);
    _ff_sum_4_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[12],t0,t1);
    _ff_sum_5_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[11],t0,t1);
    _ff_sum_6_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[10],t0,t1);
    _ff_sum_7_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[9],t0,t1);
    _ff_sum_8_mults_arr(t0,f+1,g);  _ff_add(o[8],t0,f[0]);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_18 (ff_t o[18], ff_t r[18])
{
   ff_t f[9], g[9];
    register ff_t t0, t1;

    ff_poly_from_roots_9(f,r);  ff_poly_from_roots_9(g,r+9);     _ff_add(o[17],f[8],g[8]);
    _ff_mult(t0,f[8],g[8]); _ff_add(t1,f[7],g[7]); _ff_add(o[16],t0,t1);
    _ff_sum_2_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[15],t0,t1);
    _ff_sum_3_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[14],t0,t1);
    _ff_sum_4_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[13],t0,t1);
    _ff_sum_5_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[12],t0,t1);
    _ff_sum_6_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[11],t0,t1);
    _ff_sum_7_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[10],t0,t1);
    _ff_sum_8_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[9],t0,t1);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_19 (ff_t o[19], ff_t r[19])
{
   ff_t f[10], g[9];
    register ff_t t0, t1;

    ff_poly_from_roots_10(f,r);  ff_poly_from_roots_9(g,r+10);     _ff_add(o[18],f[9],g[8]);
    _ff_mult(t0,f[9],g[8]); _ff_add(t1,f[8],g[7]); _ff_add(o[17],t0,t1);
    _ff_sum_2_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[16],t0,t1);
    _ff_sum_3_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[15],t0,t1);
    _ff_sum_4_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[14],t0,t1);
    _ff_sum_5_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[13],t0,t1);
    _ff_sum_6_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[12],t0,t1);
    _ff_sum_7_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[11],t0,t1);
    _ff_sum_8_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[10],t0,t1);
    _ff_sum_9_mults_arr(t0,f+1,g);  _ff_add(o[9],t0,f[0]);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_20 (ff_t o[20], ff_t r[20])
{
   ff_t f[10], g[10];
    register ff_t t0, t1;

    ff_poly_from_roots_10(f,r);  ff_poly_from_roots_10(g,r+10);     _ff_add(o[19],f[9],g[9]);
    _ff_mult(t0,f[9],g[9]); _ff_add(t1,f[8],g[8]); _ff_add(o[18],t0,t1);
    _ff_sum_2_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[17],t0,t1);
    _ff_sum_3_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[16],t0,t1);
    _ff_sum_4_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[15],t0,t1);
    _ff_sum_5_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[14],t0,t1);
    _ff_sum_6_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[13],t0,t1);
    _ff_sum_7_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[12],t0,t1);
    _ff_sum_8_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[11],t0,t1);
    _ff_sum_9_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[10],t0,t1);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_21 (ff_t o[21], ff_t r[21])
{
   ff_t f[11], g[10];
    register ff_t t0, t1;

    ff_poly_from_roots_11(f,r);  ff_poly_from_roots_10(g,r+11);     _ff_add(o[20],f[10],g[9]);
    _ff_mult(t0,f[10],g[9]); _ff_add(t1,f[9],g[8]); _ff_add(o[19],t0,t1);
    _ff_sum_2_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[18],t0,t1);
    _ff_sum_3_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[17],t0,t1);
    _ff_sum_4_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[16],t0,t1);
    _ff_sum_5_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[15],t0,t1);
    _ff_sum_6_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[14],t0,t1);
    _ff_sum_7_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[13],t0,t1);
    _ff_sum_8_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[12],t0,t1);
    _ff_sum_9_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[11],t0,t1);
    _ff_sum_10_mults_arr(t0,f+1,g);  _ff_add(o[10],t0,f[0]);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_22 (ff_t o[22], ff_t r[22])
{
   ff_t f[11], g[11];
    register ff_t t0, t1;

    ff_poly_from_roots_11(f,r);  ff_poly_from_roots_11(g,r+11);     _ff_add(o[21],f[10],g[10]);
    _ff_mult(t0,f[10],g[10]); _ff_add(t1,f[9],g[9]); _ff_add(o[20],t0,t1);
    _ff_sum_2_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[19],t0,t1);
    _ff_sum_3_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[18],t0,t1);
    _ff_sum_4_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[17],t0,t1);
    _ff_sum_5_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[16],t0,t1);
    _ff_sum_6_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[15],t0,t1);
    _ff_sum_7_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[14],t0,t1);
    _ff_sum_8_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[13],t0,t1);
    _ff_sum_9_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[12],t0,t1);
    _ff_sum_10_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[11],t0,t1);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_23 (ff_t o[23], ff_t r[23])
{
   ff_t f[12], g[11];
    register ff_t t0, t1;

    ff_poly_from_roots_12(f,r);  ff_poly_from_roots_11(g,r+12);     _ff_add(o[22],f[11],g[10]);
    _ff_mult(t0,f[11],g[10]); _ff_add(t1,f[10],g[9]); _ff_add(o[21],t0,t1);
    _ff_sum_2_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[20],t0,t1);
    _ff_sum_3_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[19],t0,t1);
    _ff_sum_4_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[18],t0,t1);
    _ff_sum_5_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[17],t0,t1);
    _ff_sum_6_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[16],t0,t1);
    _ff_sum_7_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[15],t0,t1);
    _ff_sum_8_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[14],t0,t1);
    _ff_sum_9_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[13],t0,t1);
    _ff_sum_10_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[12],t0,t1);
    _ff_sum_11_mults_arr(t0,f+1,g);  _ff_add(o[11],t0,f[0]);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_24 (ff_t o[24], ff_t r[24])
{
   ff_t f[12], g[12];
    register ff_t t0, t1;

    ff_poly_from_roots_12(f,r);  ff_poly_from_roots_12(g,r+12);     _ff_add(o[23],f[11],g[11]);
    _ff_mult(t0,f[11],g[11]); _ff_add(t1,f[10],g[10]); _ff_add(o[22],t0,t1);
    _ff_sum_2_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[21],t0,t1);
    _ff_sum_3_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[20],t0,t1);
    _ff_sum_4_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[19],t0,t1);
    _ff_sum_5_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[18],t0,t1);
    _ff_sum_6_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[17],t0,t1);
    _ff_sum_7_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[16],t0,t1);
    _ff_sum_8_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[15],t0,t1);
    _ff_sum_9_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[14],t0,t1);
    _ff_sum_10_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[13],t0,t1);
    _ff_sum_11_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[12],t0,t1);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_25 (ff_t o[25], ff_t r[25])
{
   ff_t f[13], g[12];
    register ff_t t0, t1;

    ff_poly_from_roots_13(f,r);  ff_poly_from_roots_12(g,r+13);     _ff_add(o[24],f[12],g[11]);
    _ff_mult(t0,f[12],g[11]); _ff_add(t1,f[11],g[10]); _ff_add(o[23],t0,t1);
    _ff_sum_2_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[22],t0,t1);
    _ff_sum_3_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[21],t0,t1);
    _ff_sum_4_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[20],t0,t1);
    _ff_sum_5_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[19],t0,t1);
    _ff_sum_6_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[18],t0,t1);
    _ff_sum_7_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[17],t0,t1);
    _ff_sum_8_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[16],t0,t1);
    _ff_sum_9_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[15],t0,t1);
    _ff_sum_10_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[14],t0,t1);
    _ff_sum_11_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[13],t0,t1);
    _ff_sum_12_mults_arr(t0,f+1,g);  _ff_add(o[12],t0,f[0]);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_26 (ff_t o[26], ff_t r[26])
{
   ff_t f[13], g[13];
    register ff_t t0, t1;

    ff_poly_from_roots_13(f,r);  ff_poly_from_roots_13(g,r+13);     _ff_add(o[25],f[12],g[12]);
    _ff_mult(t0,f[12],g[12]); _ff_add(t1,f[11],g[11]); _ff_add(o[24],t0,t1);
    _ff_sum_2_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[23],t0,t1);
    _ff_sum_3_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[22],t0,t1);
    _ff_sum_4_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[21],t0,t1);
    _ff_sum_5_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[20],t0,t1);
    _ff_sum_6_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[19],t0,t1);
    _ff_sum_7_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[18],t0,t1);
    _ff_sum_8_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[17],t0,t1);
    _ff_sum_9_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[16],t0,t1);
    _ff_sum_10_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[15],t0,t1);
    _ff_sum_11_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[14],t0,t1);
    _ff_sum_12_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[13],t0,t1);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_27 (ff_t o[27], ff_t r[27])
{
   ff_t f[14], g[13];
    register ff_t t0, t1;

    ff_poly_from_roots_14(f,r);  ff_poly_from_roots_13(g,r+14);     _ff_add(o[26],f[13],g[12]);
    _ff_mult(t0,f[13],g[12]); _ff_add(t1,f[12],g[11]); _ff_add(o[25],t0,t1);
    _ff_sum_2_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[24],t0,t1);
    _ff_sum_3_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[23],t0,t1);
    _ff_sum_4_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[22],t0,t1);
    _ff_sum_5_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[21],t0,t1);
    _ff_sum_6_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[20],t0,t1);
    _ff_sum_7_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[19],t0,t1);
    _ff_sum_8_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[18],t0,t1);
    _ff_sum_9_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[17],t0,t1);
    _ff_sum_10_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[16],t0,t1);
    _ff_sum_11_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[15],t0,t1);
    _ff_sum_12_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[14],t0,t1);
    _ff_sum_13_mults_arr(t0,f+1,g);  _ff_add(o[13],t0,f[0]);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_28 (ff_t o[28], ff_t r[28])
{
   ff_t f[14], g[14];
    register ff_t t0, t1;

    ff_poly_from_roots_14(f,r);  ff_poly_from_roots_14(g,r+14);     _ff_add(o[27],f[13],g[13]);
    _ff_mult(t0,f[13],g[13]); _ff_add(t1,f[12],g[12]); _ff_add(o[26],t0,t1);
    _ff_sum_2_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[25],t0,t1);
    _ff_sum_3_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[24],t0,t1);
    _ff_sum_4_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[23],t0,t1);
    _ff_sum_5_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[22],t0,t1);
    _ff_sum_6_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[21],t0,t1);
    _ff_sum_7_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[20],t0,t1);
    _ff_sum_8_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[19],t0,t1);
    _ff_sum_9_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[18],t0,t1);
    _ff_sum_10_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[17],t0,t1);
    _ff_sum_11_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[16],t0,t1);
    _ff_sum_12_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[15],t0,t1);
    _ff_sum_13_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[14],t0,t1);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_29 (ff_t o[29], ff_t r[29])
{
   ff_t f[15], g[14];
    register ff_t t0, t1;

    ff_poly_from_roots_15(f,r);  ff_poly_from_roots_14(g,r+15);     _ff_add(o[28],f[14],g[13]);
    _ff_mult(t0,f[14],g[13]); _ff_add(t1,f[13],g[12]); _ff_add(o[27],t0,t1);
    _ff_sum_2_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[26],t0,t1);
    _ff_sum_3_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[25],t0,t1);
    _ff_sum_4_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[24],t0,t1);
    _ff_sum_5_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[23],t0,t1);
    _ff_sum_6_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[22],t0,t1);
    _ff_sum_7_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[21],t0,t1);
    _ff_sum_8_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[20],t0,t1);
    _ff_sum_9_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[19],t0,t1);
    _ff_sum_10_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[18],t0,t1);
    _ff_sum_11_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[17],t0,t1);
    _ff_sum_12_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[16],t0,t1);
    _ff_sum_13_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[15],t0,t1);
    _ff_sum_14_mults_arr(t0,f+1,g);  _ff_add(o[14],t0,f[0]);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_30 (ff_t o[30], ff_t r[30])
{
   ff_t f[15], g[15];
    register ff_t t0, t1;

    ff_poly_from_roots_15(f,r);  ff_poly_from_roots_15(g,r+15);     _ff_add(o[29],f[14],g[14]);
    _ff_mult(t0,f[14],g[14]); _ff_add(t1,f[13],g[13]); _ff_add(o[28],t0,t1);
    _ff_sum_2_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[27],t0,t1);
    _ff_sum_3_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[26],t0,t1);
    _ff_sum_4_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[25],t0,t1);
    _ff_sum_5_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[24],t0,t1);
    _ff_sum_6_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[23],t0,t1);
    _ff_sum_7_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[22],t0,t1);
    _ff_sum_8_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[21],t0,t1);
    _ff_sum_9_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[20],t0,t1);
    _ff_sum_10_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[19],t0,t1);
    _ff_sum_11_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[18],t0,t1);
    _ff_sum_12_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[17],t0,t1);
    _ff_sum_13_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[16],t0,t1);
    _ff_sum_14_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[15],t0,t1);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_31 (ff_t o[31], ff_t r[31])
{
   ff_t f[16], g[15];
    register ff_t t0, t1;

    ff_poly_from_roots_16(f,r);  ff_poly_from_roots_15(g,r+16);     _ff_add(o[30],f[15],g[14]);
    _ff_mult(t0,f[15],g[14]); _ff_add(t1,f[14],g[13]); _ff_add(o[29],t0,t1);
    _ff_sum_2_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[28],t0,t1);
    _ff_sum_3_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[27],t0,t1);
    _ff_sum_4_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[26],t0,t1);
    _ff_sum_5_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[25],t0,t1);
    _ff_sum_6_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[24],t0,t1);
    _ff_sum_7_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[23],t0,t1);
    _ff_sum_8_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[22],t0,t1);
    _ff_sum_9_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[21],t0,t1);
    _ff_sum_10_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[20],t0,t1);
    _ff_sum_11_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[19],t0,t1);
    _ff_sum_12_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[18],t0,t1);
    _ff_sum_13_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[17],t0,t1);
    _ff_sum_14_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[16],t0,t1);
    _ff_sum_15_mults_arr(t0,f+1,g);  _ff_add(o[15],t0,f[0]);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_32 (ff_t o[32], ff_t r[32])
{
   ff_t f[16], g[16];
    register ff_t t0, t1;

    ff_poly_from_roots_16(f,r);  ff_poly_from_roots_16(g,r+16);     _ff_add(o[31],f[15],g[15]);
    _ff_mult(t0,f[15],g[15]); _ff_add(t1,f[14],g[14]); _ff_add(o[30],t0,t1);
    _ff_sum_2_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[29],t0,t1);
    _ff_sum_3_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[28],t0,t1);
    _ff_sum_4_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[27],t0,t1);
    _ff_sum_5_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[26],t0,t1);
    _ff_sum_6_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[25],t0,t1);
    _ff_sum_7_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[24],t0,t1);
    _ff_sum_8_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[23],t0,t1);
    _ff_sum_9_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[22],t0,t1);
    _ff_sum_10_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[21],t0,t1);
    _ff_sum_11_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[20],t0,t1);
    _ff_sum_12_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[19],t0,t1);
    _ff_sum_13_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[18],t0,t1);
    _ff_sum_14_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[17],t0,t1);
    _ff_sum_15_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[16],t0,t1);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_33 (ff_t o[33], ff_t r[33])
{
   ff_t f[17], g[16];
    register ff_t t0, t1;

    ff_poly_from_roots_17(f,r);  ff_poly_from_roots_16(g,r+17);     _ff_add(o[32],f[16],g[15]);
    _ff_mult(t0,f[16],g[15]); _ff_add(t1,f[15],g[14]); _ff_add(o[31],t0,t1);
    _ff_sum_2_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[30],t0,t1);
    _ff_sum_3_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[29],t0,t1);
    _ff_sum_4_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[28],t0,t1);
    _ff_sum_5_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[27],t0,t1);
    _ff_sum_6_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[26],t0,t1);
    _ff_sum_7_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[25],t0,t1);
    _ff_sum_8_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[24],t0,t1);
    _ff_sum_9_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[23],t0,t1);
    _ff_sum_10_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[22],t0,t1);
    _ff_sum_11_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[21],t0,t1);
    _ff_sum_12_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[20],t0,t1);
    _ff_sum_13_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[19],t0,t1);
    _ff_sum_14_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[18],t0,t1);
    _ff_sum_15_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[17],t0,t1);
    _ff_sum_16_mults_arr(t0,f+1,g);  _ff_add(o[16],t0,f[0]);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_34 (ff_t o[34], ff_t r[34])
{
   ff_t f[17], g[17];
    register ff_t t0, t1;

    ff_poly_from_roots_17(f,r);  ff_poly_from_roots_17(g,r+17);     _ff_add(o[33],f[16],g[16]);
    _ff_mult(t0,f[16],g[16]); _ff_add(t1,f[15],g[15]); _ff_add(o[32],t0,t1);
    _ff_sum_2_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[31],t0,t1);
    _ff_sum_3_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[30],t0,t1);
    _ff_sum_4_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[29],t0,t1);
    _ff_sum_5_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[28],t0,t1);
    _ff_sum_6_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[27],t0,t1);
    _ff_sum_7_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[26],t0,t1);
    _ff_sum_8_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[25],t0,t1);
    _ff_sum_9_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[24],t0,t1);
    _ff_sum_10_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[23],t0,t1);
    _ff_sum_11_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[22],t0,t1);
    _ff_sum_12_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[21],t0,t1);
    _ff_sum_13_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[20],t0,t1);
    _ff_sum_14_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[19],t0,t1);
    _ff_sum_15_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[18],t0,t1);
    _ff_sum_16_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[17],t0,t1);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_35 (ff_t o[35], ff_t r[35])
{
   ff_t f[18], g[17];
    register ff_t t0, t1;

    ff_poly_from_roots_18(f,r);  ff_poly_from_roots_17(g,r+18);     _ff_add(o[34],f[17],g[16]);
    _ff_mult(t0,f[17],g[16]); _ff_add(t1,f[16],g[15]); _ff_add(o[33],t0,t1);
    _ff_sum_2_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[32],t0,t1);
    _ff_sum_3_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[31],t0,t1);
    _ff_sum_4_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[30],t0,t1);
    _ff_sum_5_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[29],t0,t1);
    _ff_sum_6_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[28],t0,t1);
    _ff_sum_7_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[27],t0,t1);
    _ff_sum_8_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[26],t0,t1);
    _ff_sum_9_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[25],t0,t1);
    _ff_sum_10_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[24],t0,t1);
    _ff_sum_11_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[23],t0,t1);
    _ff_sum_12_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[22],t0,t1);
    _ff_sum_13_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[21],t0,t1);
    _ff_sum_14_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[20],t0,t1);
    _ff_sum_15_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[19],t0,t1);
    _ff_sum_16_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[18],t0,t1);
    _ff_sum_17_mults_arr(t0,f+1,g);  _ff_add(o[17],t0,f[0]);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_36 (ff_t o[36], ff_t r[36])
{
   ff_t f[18], g[18];
    register ff_t t0, t1;

    ff_poly_from_roots_18(f,r);  ff_poly_from_roots_18(g,r+18);     _ff_add(o[35],f[17],g[17]);
    _ff_mult(t0,f[17],g[17]); _ff_add(t1,f[16],g[16]); _ff_add(o[34],t0,t1);
    _ff_sum_2_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[33],t0,t1);
    _ff_sum_3_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[32],t0,t1);
    _ff_sum_4_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[31],t0,t1);
    _ff_sum_5_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[30],t0,t1);
    _ff_sum_6_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[29],t0,t1);
    _ff_sum_7_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[28],t0,t1);
    _ff_sum_8_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[27],t0,t1);
    _ff_sum_9_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[26],t0,t1);
    _ff_sum_10_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[25],t0,t1);
    _ff_sum_11_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[24],t0,t1);
    _ff_sum_12_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[23],t0,t1);
    _ff_sum_13_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[22],t0,t1);
    _ff_sum_14_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[21],t0,t1);
    _ff_sum_15_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[20],t0,t1);
    _ff_sum_16_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[19],t0,t1);
    _ff_sum_17_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[18],t0,t1);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_37 (ff_t o[37], ff_t r[37])
{
   ff_t f[19], g[18];
    register ff_t t0, t1;

    ff_poly_from_roots_19(f,r);  ff_poly_from_roots_18(g,r+19);     _ff_add(o[36],f[18],g[17]);
    _ff_mult(t0,f[18],g[17]); _ff_add(t1,f[17],g[16]); _ff_add(o[35],t0,t1);
    _ff_sum_2_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[34],t0,t1);
    _ff_sum_3_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[33],t0,t1);
    _ff_sum_4_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[32],t0,t1);
    _ff_sum_5_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[31],t0,t1);
    _ff_sum_6_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[30],t0,t1);
    _ff_sum_7_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[29],t0,t1);
    _ff_sum_8_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[28],t0,t1);
    _ff_sum_9_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[27],t0,t1);
    _ff_sum_10_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[26],t0,t1);
    _ff_sum_11_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[25],t0,t1);
    _ff_sum_12_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[24],t0,t1);
    _ff_sum_13_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[23],t0,t1);
    _ff_sum_14_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[22],t0,t1);
    _ff_sum_15_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[21],t0,t1);
    _ff_sum_16_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[20],t0,t1);
    _ff_sum_17_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[19],t0,t1);
    _ff_sum_18_mults_arr(t0,f+1,g);  _ff_add(o[18],t0,f[0]);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_38 (ff_t o[38], ff_t r[38])
{
   ff_t f[19], g[19];
    register ff_t t0, t1;

    ff_poly_from_roots_19(f,r);  ff_poly_from_roots_19(g,r+19);     _ff_add(o[37],f[18],g[18]);
    _ff_mult(t0,f[18],g[18]); _ff_add(t1,f[17],g[17]); _ff_add(o[36],t0,t1);
    _ff_sum_2_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[35],t0,t1);
    _ff_sum_3_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[34],t0,t1);
    _ff_sum_4_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[33],t0,t1);
    _ff_sum_5_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[32],t0,t1);
    _ff_sum_6_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[31],t0,t1);
    _ff_sum_7_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[30],t0,t1);
    _ff_sum_8_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[29],t0,t1);
    _ff_sum_9_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[28],t0,t1);
    _ff_sum_10_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[27],t0,t1);
    _ff_sum_11_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[26],t0,t1);
    _ff_sum_12_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[25],t0,t1);
    _ff_sum_13_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[24],t0,t1);
    _ff_sum_14_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[23],t0,t1);
    _ff_sum_15_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[22],t0,t1);
    _ff_sum_16_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[21],t0,t1);
    _ff_sum_17_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[20],t0,t1);
    _ff_sum_18_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[19],t0,t1);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_39 (ff_t o[39], ff_t r[39])
{
   ff_t f[20], g[19];
    register ff_t t0, t1;

    ff_poly_from_roots_20(f,r);  ff_poly_from_roots_19(g,r+20);     _ff_add(o[38],f[19],g[18]);
    _ff_mult(t0,f[19],g[18]); _ff_add(t1,f[18],g[17]); _ff_add(o[37],t0,t1);
    _ff_sum_2_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[36],t0,t1);
    _ff_sum_3_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[35],t0,t1);
    _ff_sum_4_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[34],t0,t1);
    _ff_sum_5_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[33],t0,t1);
    _ff_sum_6_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[32],t0,t1);
    _ff_sum_7_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[31],t0,t1);
    _ff_sum_8_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[30],t0,t1);
    _ff_sum_9_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[29],t0,t1);
    _ff_sum_10_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[28],t0,t1);
    _ff_sum_11_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[27],t0,t1);
    _ff_sum_12_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[26],t0,t1);
    _ff_sum_13_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[25],t0,t1);
    _ff_sum_14_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[24],t0,t1);
    _ff_sum_15_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[23],t0,t1);
    _ff_sum_16_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[22],t0,t1);
    _ff_sum_17_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[21],t0,t1);
    _ff_sum_18_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[20],t0,t1);
    _ff_sum_19_mults_arr(t0,f+1,g);  _ff_add(o[19],t0,f[0]);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_40 (ff_t o[40], ff_t r[40])
{
   ff_t f[20], g[20];
    register ff_t t0, t1;

    ff_poly_from_roots_20(f,r);  ff_poly_from_roots_20(g,r+20);     _ff_add(o[39],f[19],g[19]);
    _ff_mult(t0,f[19],g[19]); _ff_add(t1,f[18],g[18]); _ff_add(o[38],t0,t1);
    _ff_sum_2_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[37],t0,t1);
    _ff_sum_3_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[36],t0,t1);
    _ff_sum_4_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[35],t0,t1);
    _ff_sum_5_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[34],t0,t1);
    _ff_sum_6_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[33],t0,t1);
    _ff_sum_7_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[32],t0,t1);
    _ff_sum_8_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[31],t0,t1);
    _ff_sum_9_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[30],t0,t1);
    _ff_sum_10_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[29],t0,t1);
    _ff_sum_11_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[28],t0,t1);
    _ff_sum_12_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[27],t0,t1);
    _ff_sum_13_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[26],t0,t1);
    _ff_sum_14_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[25],t0,t1);
    _ff_sum_15_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[24],t0,t1);
    _ff_sum_16_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[23],t0,t1);
    _ff_sum_17_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[22],t0,t1);
    _ff_sum_18_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[21],t0,t1);
    _ff_sum_19_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[20],t0,t1);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_41 (ff_t o[41], ff_t r[41])
{
   ff_t f[21], g[20];
    register ff_t t0, t1;

    ff_poly_from_roots_21(f,r);  ff_poly_from_roots_20(g,r+21);     _ff_add(o[40],f[20],g[19]);
    _ff_mult(t0,f[20],g[19]); _ff_add(t1,f[19],g[18]); _ff_add(o[39],t0,t1);
    _ff_sum_2_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[38],t0,t1);
    _ff_sum_3_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[37],t0,t1);
    _ff_sum_4_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[36],t0,t1);
    _ff_sum_5_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[35],t0,t1);
    _ff_sum_6_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[34],t0,t1);
    _ff_sum_7_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[33],t0,t1);
    _ff_sum_8_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[32],t0,t1);
    _ff_sum_9_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[31],t0,t1);
    _ff_sum_10_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[30],t0,t1);
    _ff_sum_11_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[29],t0,t1);
    _ff_sum_12_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[28],t0,t1);
    _ff_sum_13_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[27],t0,t1);
    _ff_sum_14_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[26],t0,t1);
    _ff_sum_15_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[25],t0,t1);
    _ff_sum_16_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[24],t0,t1);
    _ff_sum_17_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[23],t0,t1);
    _ff_sum_18_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[22],t0,t1);
    _ff_sum_19_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[21],t0,t1);
    _ff_sum_20_mults_arr(t0,f+1,g);  _ff_add(o[20],t0,f[0]);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_42 (ff_t o[42], ff_t r[42])
{
   ff_t f[21], g[21];
    register ff_t t0, t1;

    ff_poly_from_roots_21(f,r);  ff_poly_from_roots_21(g,r+21);     _ff_add(o[41],f[20],g[20]);
    _ff_mult(t0,f[20],g[20]); _ff_add(t1,f[19],g[19]); _ff_add(o[40],t0,t1);
    _ff_sum_2_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[39],t0,t1);
    _ff_sum_3_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[38],t0,t1);
    _ff_sum_4_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[37],t0,t1);
    _ff_sum_5_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[36],t0,t1);
    _ff_sum_6_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[35],t0,t1);
    _ff_sum_7_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[34],t0,t1);
    _ff_sum_8_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[33],t0,t1);
    _ff_sum_9_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[32],t0,t1);
    _ff_sum_10_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[31],t0,t1);
    _ff_sum_11_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[30],t0,t1);
    _ff_sum_12_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[29],t0,t1);
    _ff_sum_13_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[28],t0,t1);
    _ff_sum_14_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[27],t0,t1);
    _ff_sum_15_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[26],t0,t1);
    _ff_sum_16_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[25],t0,t1);
    _ff_sum_17_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[24],t0,t1);
    _ff_sum_18_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[23],t0,t1);
    _ff_sum_19_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[22],t0,t1);
    _ff_sum_20_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[21],t0,t1);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_43 (ff_t o[43], ff_t r[43])
{
   ff_t f[22], g[21];
    register ff_t t0, t1;

    ff_poly_from_roots_22(f,r);  ff_poly_from_roots_21(g,r+22);     _ff_add(o[42],f[21],g[20]);
    _ff_mult(t0,f[21],g[20]); _ff_add(t1,f[20],g[19]); _ff_add(o[41],t0,t1);
    _ff_sum_2_mults_arr(t0,f+20,g+19); _ff_add(t1,f[19],g[18]);  _ff_add(o[40],t0,t1);
    _ff_sum_3_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[39],t0,t1);
    _ff_sum_4_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[38],t0,t1);
    _ff_sum_5_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[37],t0,t1);
    _ff_sum_6_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[36],t0,t1);
    _ff_sum_7_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[35],t0,t1);
    _ff_sum_8_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[34],t0,t1);
    _ff_sum_9_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[33],t0,t1);
    _ff_sum_10_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[32],t0,t1);
    _ff_sum_11_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[31],t0,t1);
    _ff_sum_12_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[30],t0,t1);
    _ff_sum_13_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[29],t0,t1);
    _ff_sum_14_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[28],t0,t1);
    _ff_sum_15_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[27],t0,t1);
    _ff_sum_16_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[26],t0,t1);
    _ff_sum_17_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[25],t0,t1);
    _ff_sum_18_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[24],t0,t1);
    _ff_sum_19_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[23],t0,t1);
    _ff_sum_20_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[22],t0,t1);
    _ff_sum_21_mults_arr(t0,f+1,g);  _ff_add(o[21],t0,f[0]);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_44 (ff_t o[44], ff_t r[44])
{
   ff_t f[22], g[22];
    register ff_t t0, t1;

    ff_poly_from_roots_22(f,r);  ff_poly_from_roots_22(g,r+22);     _ff_add(o[43],f[21],g[21]);
    _ff_mult(t0,f[21],g[21]); _ff_add(t1,f[20],g[20]); _ff_add(o[42],t0,t1);
    _ff_sum_2_mults_arr(t0,f+20,g+20); _ff_add(t1,f[19],g[19]);  _ff_add(o[41],t0,t1);
    _ff_sum_3_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[40],t0,t1);
    _ff_sum_4_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[39],t0,t1);
    _ff_sum_5_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[38],t0,t1);
    _ff_sum_6_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[37],t0,t1);
    _ff_sum_7_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[36],t0,t1);
    _ff_sum_8_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[35],t0,t1);
    _ff_sum_9_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[34],t0,t1);
    _ff_sum_10_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[33],t0,t1);
    _ff_sum_11_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[32],t0,t1);
    _ff_sum_12_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[31],t0,t1);
    _ff_sum_13_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[30],t0,t1);
    _ff_sum_14_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[29],t0,t1);
    _ff_sum_15_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[28],t0,t1);
    _ff_sum_16_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[27],t0,t1);
    _ff_sum_17_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[26],t0,t1);
    _ff_sum_18_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[25],t0,t1);
    _ff_sum_19_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[24],t0,t1);
    _ff_sum_20_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[23],t0,t1);
    _ff_sum_21_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[22],t0,t1);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_45 (ff_t o[45], ff_t r[45])
{
   ff_t f[23], g[22];
    register ff_t t0, t1;

    ff_poly_from_roots_23(f,r);  ff_poly_from_roots_22(g,r+23);     _ff_add(o[44],f[22],g[21]);
    _ff_mult(t0,f[22],g[21]); _ff_add(t1,f[21],g[20]); _ff_add(o[43],t0,t1);
    _ff_sum_2_mults_arr(t0,f+21,g+20); _ff_add(t1,f[20],g[19]);  _ff_add(o[42],t0,t1);
    _ff_sum_3_mults_arr(t0,f+20,g+19); _ff_add(t1,f[19],g[18]);  _ff_add(o[41],t0,t1);
    _ff_sum_4_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[40],t0,t1);
    _ff_sum_5_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[39],t0,t1);
    _ff_sum_6_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[38],t0,t1);
    _ff_sum_7_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[37],t0,t1);
    _ff_sum_8_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[36],t0,t1);
    _ff_sum_9_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[35],t0,t1);
    _ff_sum_10_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[34],t0,t1);
    _ff_sum_11_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[33],t0,t1);
    _ff_sum_12_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[32],t0,t1);
    _ff_sum_13_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[31],t0,t1);
    _ff_sum_14_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[30],t0,t1);
    _ff_sum_15_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[29],t0,t1);
    _ff_sum_16_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[28],t0,t1);
    _ff_sum_17_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[27],t0,t1);
    _ff_sum_18_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[26],t0,t1);
    _ff_sum_19_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[25],t0,t1);
    _ff_sum_20_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[24],t0,t1);
    _ff_sum_21_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[23],t0,t1);
    _ff_sum_22_mults_arr(t0,f+1,g);  _ff_add(o[22],t0,f[0]);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_46 (ff_t o[46], ff_t r[46])
{
   ff_t f[23], g[23];
    register ff_t t0, t1;

    ff_poly_from_roots_23(f,r);  ff_poly_from_roots_23(g,r+23);     _ff_add(o[45],f[22],g[22]);
    _ff_mult(t0,f[22],g[22]); _ff_add(t1,f[21],g[21]); _ff_add(o[44],t0,t1);
    _ff_sum_2_mults_arr(t0,f+21,g+21); _ff_add(t1,f[20],g[20]);  _ff_add(o[43],t0,t1);
    _ff_sum_3_mults_arr(t0,f+20,g+20); _ff_add(t1,f[19],g[19]);  _ff_add(o[42],t0,t1);
    _ff_sum_4_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[41],t0,t1);
    _ff_sum_5_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[40],t0,t1);
    _ff_sum_6_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[39],t0,t1);
    _ff_sum_7_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[38],t0,t1);
    _ff_sum_8_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[37],t0,t1);
    _ff_sum_9_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[36],t0,t1);
    _ff_sum_10_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[35],t0,t1);
    _ff_sum_11_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[34],t0,t1);
    _ff_sum_12_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[33],t0,t1);
    _ff_sum_13_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[32],t0,t1);
    _ff_sum_14_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[31],t0,t1);
    _ff_sum_15_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[30],t0,t1);
    _ff_sum_16_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[29],t0,t1);
    _ff_sum_17_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[28],t0,t1);
    _ff_sum_18_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[27],t0,t1);
    _ff_sum_19_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[26],t0,t1);
    _ff_sum_20_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[25],t0,t1);
    _ff_sum_21_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[24],t0,t1);
    _ff_sum_22_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[23],t0,t1);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_47 (ff_t o[47], ff_t r[47])
{
   ff_t f[24], g[23];
    register ff_t t0, t1;

    ff_poly_from_roots_24(f,r);  ff_poly_from_roots_23(g,r+24);     _ff_add(o[46],f[23],g[22]);
    _ff_mult(t0,f[23],g[22]); _ff_add(t1,f[22],g[21]); _ff_add(o[45],t0,t1);
    _ff_sum_2_mults_arr(t0,f+22,g+21); _ff_add(t1,f[21],g[20]);  _ff_add(o[44],t0,t1);
    _ff_sum_3_mults_arr(t0,f+21,g+20); _ff_add(t1,f[20],g[19]);  _ff_add(o[43],t0,t1);
    _ff_sum_4_mults_arr(t0,f+20,g+19); _ff_add(t1,f[19],g[18]);  _ff_add(o[42],t0,t1);
    _ff_sum_5_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[41],t0,t1);
    _ff_sum_6_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[40],t0,t1);
    _ff_sum_7_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[39],t0,t1);
    _ff_sum_8_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[38],t0,t1);
    _ff_sum_9_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[37],t0,t1);
    _ff_sum_10_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[36],t0,t1);
    _ff_sum_11_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[35],t0,t1);
    _ff_sum_12_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[34],t0,t1);
    _ff_sum_13_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[33],t0,t1);
    _ff_sum_14_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[32],t0,t1);
    _ff_sum_15_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[31],t0,t1);
    _ff_sum_16_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[30],t0,t1);
    _ff_sum_17_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[29],t0,t1);
    _ff_sum_18_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[28],t0,t1);
    _ff_sum_19_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[27],t0,t1);
    _ff_sum_20_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[26],t0,t1);
    _ff_sum_21_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[25],t0,t1);
    _ff_sum_22_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[24],t0,t1);
    _ff_sum_23_mults_arr(t0,f+1,g);  _ff_add(o[23],t0,f[0]);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_48 (ff_t o[48], ff_t r[48])
{
   ff_t f[24], g[24];
    register ff_t t0, t1;

    ff_poly_from_roots_24(f,r);  ff_poly_from_roots_24(g,r+24);     _ff_add(o[47],f[23],g[23]);
    _ff_mult(t0,f[23],g[23]); _ff_add(t1,f[22],g[22]); _ff_add(o[46],t0,t1);
    _ff_sum_2_mults_arr(t0,f+22,g+22); _ff_add(t1,f[21],g[21]);  _ff_add(o[45],t0,t1);
    _ff_sum_3_mults_arr(t0,f+21,g+21); _ff_add(t1,f[20],g[20]);  _ff_add(o[44],t0,t1);
    _ff_sum_4_mults_arr(t0,f+20,g+20); _ff_add(t1,f[19],g[19]);  _ff_add(o[43],t0,t1);
    _ff_sum_5_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[42],t0,t1);
    _ff_sum_6_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[41],t0,t1);
    _ff_sum_7_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[40],t0,t1);
    _ff_sum_8_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[39],t0,t1);
    _ff_sum_9_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[38],t0,t1);
    _ff_sum_10_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[37],t0,t1);
    _ff_sum_11_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[36],t0,t1);
    _ff_sum_12_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[35],t0,t1);
    _ff_sum_13_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[34],t0,t1);
    _ff_sum_14_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[33],t0,t1);
    _ff_sum_15_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[32],t0,t1);
    _ff_sum_16_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[31],t0,t1);
    _ff_sum_17_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[30],t0,t1);
    _ff_sum_18_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[29],t0,t1);
    _ff_sum_19_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[28],t0,t1);
    _ff_sum_20_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[27],t0,t1);
    _ff_sum_21_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[26],t0,t1);
    _ff_sum_22_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[25],t0,t1);
    _ff_sum_23_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[24],t0,t1);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_49 (ff_t o[49], ff_t r[49])
{
   ff_t f[25], g[24];
    register ff_t t0, t1;

    ff_poly_from_roots_25(f,r);  ff_poly_from_roots_24(g,r+25);     _ff_add(o[48],f[24],g[23]);
    _ff_mult(t0,f[24],g[23]); _ff_add(t1,f[23],g[22]); _ff_add(o[47],t0,t1);
    _ff_sum_2_mults_arr(t0,f+23,g+22); _ff_add(t1,f[22],g[21]);  _ff_add(o[46],t0,t1);
    _ff_sum_3_mults_arr(t0,f+22,g+21); _ff_add(t1,f[21],g[20]);  _ff_add(o[45],t0,t1);
    _ff_sum_4_mults_arr(t0,f+21,g+20); _ff_add(t1,f[20],g[19]);  _ff_add(o[44],t0,t1);
    _ff_sum_5_mults_arr(t0,f+20,g+19); _ff_add(t1,f[19],g[18]);  _ff_add(o[43],t0,t1);
    _ff_sum_6_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[42],t0,t1);
    _ff_sum_7_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[41],t0,t1);
    _ff_sum_8_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[40],t0,t1);
    _ff_sum_9_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[39],t0,t1);
    _ff_sum_10_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[38],t0,t1);
    _ff_sum_11_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[37],t0,t1);
    _ff_sum_12_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[36],t0,t1);
    _ff_sum_13_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[35],t0,t1);
    _ff_sum_14_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[34],t0,t1);
    _ff_sum_15_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[33],t0,t1);
    _ff_sum_16_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[32],t0,t1);
    _ff_sum_17_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[31],t0,t1);
    _ff_sum_18_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[30],t0,t1);
    _ff_sum_19_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[29],t0,t1);
    _ff_sum_20_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[28],t0,t1);
    _ff_sum_21_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[27],t0,t1);
    _ff_sum_22_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[26],t0,t1);
    _ff_sum_23_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[25],t0,t1);
    _ff_sum_24_mults_arr(t0,f+1,g);  _ff_add(o[24],t0,f[0]);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_50 (ff_t o[50], ff_t r[50])
{
   ff_t f[25], g[25];
    register ff_t t0, t1;

    ff_poly_from_roots_25(f,r);  ff_poly_from_roots_25(g,r+25);     _ff_add(o[49],f[24],g[24]);
    _ff_mult(t0,f[24],g[24]); _ff_add(t1,f[23],g[23]); _ff_add(o[48],t0,t1);
    _ff_sum_2_mults_arr(t0,f+23,g+23); _ff_add(t1,f[22],g[22]);  _ff_add(o[47],t0,t1);
    _ff_sum_3_mults_arr(t0,f+22,g+22); _ff_add(t1,f[21],g[21]);  _ff_add(o[46],t0,t1);
    _ff_sum_4_mults_arr(t0,f+21,g+21); _ff_add(t1,f[20],g[20]);  _ff_add(o[45],t0,t1);
    _ff_sum_5_mults_arr(t0,f+20,g+20); _ff_add(t1,f[19],g[19]);  _ff_add(o[44],t0,t1);
    _ff_sum_6_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[43],t0,t1);
    _ff_sum_7_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[42],t0,t1);
    _ff_sum_8_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[41],t0,t1);
    _ff_sum_9_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[40],t0,t1);
    _ff_sum_10_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[39],t0,t1);
    _ff_sum_11_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[38],t0,t1);
    _ff_sum_12_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[37],t0,t1);
    _ff_sum_13_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[36],t0,t1);
    _ff_sum_14_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[35],t0,t1);
    _ff_sum_15_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[34],t0,t1);
    _ff_sum_16_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[33],t0,t1);
    _ff_sum_17_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[32],t0,t1);
    _ff_sum_18_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[31],t0,t1);
    _ff_sum_19_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[30],t0,t1);
    _ff_sum_20_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[29],t0,t1);
    _ff_sum_21_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[28],t0,t1);
    _ff_sum_22_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[27],t0,t1);
    _ff_sum_23_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[26],t0,t1);
    _ff_sum_24_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[25],t0,t1);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_51 (ff_t o[51], ff_t r[51])
{
   ff_t f[26], g[25];
    register ff_t t0, t1;

    ff_poly_from_roots_26(f,r);  ff_poly_from_roots_25(g,r+26);     _ff_add(o[50],f[25],g[24]);
    _ff_mult(t0,f[25],g[24]); _ff_add(t1,f[24],g[23]); _ff_add(o[49],t0,t1);
    _ff_sum_2_mults_arr(t0,f+24,g+23); _ff_add(t1,f[23],g[22]);  _ff_add(o[48],t0,t1);
    _ff_sum_3_mults_arr(t0,f+23,g+22); _ff_add(t1,f[22],g[21]);  _ff_add(o[47],t0,t1);
    _ff_sum_4_mults_arr(t0,f+22,g+21); _ff_add(t1,f[21],g[20]);  _ff_add(o[46],t0,t1);
    _ff_sum_5_mults_arr(t0,f+21,g+20); _ff_add(t1,f[20],g[19]);  _ff_add(o[45],t0,t1);
    _ff_sum_6_mults_arr(t0,f+20,g+19); _ff_add(t1,f[19],g[18]);  _ff_add(o[44],t0,t1);
    _ff_sum_7_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[43],t0,t1);
    _ff_sum_8_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[42],t0,t1);
    _ff_sum_9_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[41],t0,t1);
    _ff_sum_10_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[40],t0,t1);
    _ff_sum_11_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[39],t0,t1);
    _ff_sum_12_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[38],t0,t1);
    _ff_sum_13_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[37],t0,t1);
    _ff_sum_14_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[36],t0,t1);
    _ff_sum_15_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[35],t0,t1);
    _ff_sum_16_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[34],t0,t1);
    _ff_sum_17_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[33],t0,t1);
    _ff_sum_18_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[32],t0,t1);
    _ff_sum_19_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[31],t0,t1);
    _ff_sum_20_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[30],t0,t1);
    _ff_sum_21_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[29],t0,t1);
    _ff_sum_22_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[28],t0,t1);
    _ff_sum_23_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[27],t0,t1);
    _ff_sum_24_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[26],t0,t1);
    _ff_sum_25_mults_arr(t0,f+1,g);  _ff_add(o[25],t0,f[0]);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_52 (ff_t o[52], ff_t r[52])
{
   ff_t f[26], g[26];
    register ff_t t0, t1;

    ff_poly_from_roots_26(f,r);  ff_poly_from_roots_26(g,r+26);     _ff_add(o[51],f[25],g[25]);
    _ff_mult(t0,f[25],g[25]); _ff_add(t1,f[24],g[24]); _ff_add(o[50],t0,t1);
    _ff_sum_2_mults_arr(t0,f+24,g+24); _ff_add(t1,f[23],g[23]);  _ff_add(o[49],t0,t1);
    _ff_sum_3_mults_arr(t0,f+23,g+23); _ff_add(t1,f[22],g[22]);  _ff_add(o[48],t0,t1);
    _ff_sum_4_mults_arr(t0,f+22,g+22); _ff_add(t1,f[21],g[21]);  _ff_add(o[47],t0,t1);
    _ff_sum_5_mults_arr(t0,f+21,g+21); _ff_add(t1,f[20],g[20]);  _ff_add(o[46],t0,t1);
    _ff_sum_6_mults_arr(t0,f+20,g+20); _ff_add(t1,f[19],g[19]);  _ff_add(o[45],t0,t1);
    _ff_sum_7_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[44],t0,t1);
    _ff_sum_8_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[43],t0,t1);
    _ff_sum_9_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[42],t0,t1);
    _ff_sum_10_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[41],t0,t1);
    _ff_sum_11_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[40],t0,t1);
    _ff_sum_12_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[39],t0,t1);
    _ff_sum_13_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[38],t0,t1);
    _ff_sum_14_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[37],t0,t1);
    _ff_sum_15_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[36],t0,t1);
    _ff_sum_16_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[35],t0,t1);
    _ff_sum_17_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[34],t0,t1);
    _ff_sum_18_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[33],t0,t1);
    _ff_sum_19_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[32],t0,t1);
    _ff_sum_20_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[31],t0,t1);
    _ff_sum_21_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[30],t0,t1);
    _ff_sum_22_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[29],t0,t1);
    _ff_sum_23_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[28],t0,t1);
    _ff_sum_24_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[27],t0,t1);
    _ff_sum_25_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[26],t0,t1);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_53 (ff_t o[53], ff_t r[53])
{
   ff_t f[27], g[26];
    register ff_t t0, t1;

    ff_poly_from_roots_27(f,r);  ff_poly_from_roots_26(g,r+27);     _ff_add(o[52],f[26],g[25]);
    _ff_mult(t0,f[26],g[25]); _ff_add(t1,f[25],g[24]); _ff_add(o[51],t0,t1);
    _ff_sum_2_mults_arr(t0,f+25,g+24); _ff_add(t1,f[24],g[23]);  _ff_add(o[50],t0,t1);
    _ff_sum_3_mults_arr(t0,f+24,g+23); _ff_add(t1,f[23],g[22]);  _ff_add(o[49],t0,t1);
    _ff_sum_4_mults_arr(t0,f+23,g+22); _ff_add(t1,f[22],g[21]);  _ff_add(o[48],t0,t1);
    _ff_sum_5_mults_arr(t0,f+22,g+21); _ff_add(t1,f[21],g[20]);  _ff_add(o[47],t0,t1);
    _ff_sum_6_mults_arr(t0,f+21,g+20); _ff_add(t1,f[20],g[19]);  _ff_add(o[46],t0,t1);
    _ff_sum_7_mults_arr(t0,f+20,g+19); _ff_add(t1,f[19],g[18]);  _ff_add(o[45],t0,t1);
    _ff_sum_8_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[44],t0,t1);
    _ff_sum_9_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[43],t0,t1);
    _ff_sum_10_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[42],t0,t1);
    _ff_sum_11_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[41],t0,t1);
    _ff_sum_12_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[40],t0,t1);
    _ff_sum_13_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[39],t0,t1);
    _ff_sum_14_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[38],t0,t1);
    _ff_sum_15_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[37],t0,t1);
    _ff_sum_16_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[36],t0,t1);
    _ff_sum_17_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[35],t0,t1);
    _ff_sum_18_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[34],t0,t1);
    _ff_sum_19_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[33],t0,t1);
    _ff_sum_20_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[32],t0,t1);
    _ff_sum_21_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[31],t0,t1);
    _ff_sum_22_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[30],t0,t1);
    _ff_sum_23_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[29],t0,t1);
    _ff_sum_24_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[28],t0,t1);
    _ff_sum_25_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[27],t0,t1);
    _ff_sum_26_mults_arr(t0,f+1,g);  _ff_add(o[26],t0,f[0]);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_54 (ff_t o[54], ff_t r[54])
{
   ff_t f[27], g[27];
    register ff_t t0, t1;

    ff_poly_from_roots_27(f,r);  ff_poly_from_roots_27(g,r+27);     _ff_add(o[53],f[26],g[26]);
    _ff_mult(t0,f[26],g[26]); _ff_add(t1,f[25],g[25]); _ff_add(o[52],t0,t1);
    _ff_sum_2_mults_arr(t0,f+25,g+25); _ff_add(t1,f[24],g[24]);  _ff_add(o[51],t0,t1);
    _ff_sum_3_mults_arr(t0,f+24,g+24); _ff_add(t1,f[23],g[23]);  _ff_add(o[50],t0,t1);
    _ff_sum_4_mults_arr(t0,f+23,g+23); _ff_add(t1,f[22],g[22]);  _ff_add(o[49],t0,t1);
    _ff_sum_5_mults_arr(t0,f+22,g+22); _ff_add(t1,f[21],g[21]);  _ff_add(o[48],t0,t1);
    _ff_sum_6_mults_arr(t0,f+21,g+21); _ff_add(t1,f[20],g[20]);  _ff_add(o[47],t0,t1);
    _ff_sum_7_mults_arr(t0,f+20,g+20); _ff_add(t1,f[19],g[19]);  _ff_add(o[46],t0,t1);
    _ff_sum_8_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[45],t0,t1);
    _ff_sum_9_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[44],t0,t1);
    _ff_sum_10_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[43],t0,t1);
    _ff_sum_11_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[42],t0,t1);
    _ff_sum_12_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[41],t0,t1);
    _ff_sum_13_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[40],t0,t1);
    _ff_sum_14_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[39],t0,t1);
    _ff_sum_15_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[38],t0,t1);
    _ff_sum_16_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[37],t0,t1);
    _ff_sum_17_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[36],t0,t1);
    _ff_sum_18_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[35],t0,t1);
    _ff_sum_19_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[34],t0,t1);
    _ff_sum_20_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[33],t0,t1);
    _ff_sum_21_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[32],t0,t1);
    _ff_sum_22_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[31],t0,t1);
    _ff_sum_23_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[30],t0,t1);
    _ff_sum_24_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[29],t0,t1);
    _ff_sum_25_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[28],t0,t1);
    _ff_sum_26_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[27],t0,t1);
    _ff_sum_27_mults_arr(o[26],f,g);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_55 (ff_t o[55], ff_t r[55])
{
   ff_t f[28], g[27];
    register ff_t t0, t1;

    ff_poly_from_roots_28(f,r);  ff_poly_from_roots_27(g,r+28);     _ff_add(o[54],f[27],g[26]);
    _ff_mult(t0,f[27],g[26]); _ff_add(t1,f[26],g[25]); _ff_add(o[53],t0,t1);
    _ff_sum_2_mults_arr(t0,f+26,g+25); _ff_add(t1,f[25],g[24]);  _ff_add(o[52],t0,t1);
    _ff_sum_3_mults_arr(t0,f+25,g+24); _ff_add(t1,f[24],g[23]);  _ff_add(o[51],t0,t1);
    _ff_sum_4_mults_arr(t0,f+24,g+23); _ff_add(t1,f[23],g[22]);  _ff_add(o[50],t0,t1);
    _ff_sum_5_mults_arr(t0,f+23,g+22); _ff_add(t1,f[22],g[21]);  _ff_add(o[49],t0,t1);
    _ff_sum_6_mults_arr(t0,f+22,g+21); _ff_add(t1,f[21],g[20]);  _ff_add(o[48],t0,t1);
    _ff_sum_7_mults_arr(t0,f+21,g+20); _ff_add(t1,f[20],g[19]);  _ff_add(o[47],t0,t1);
    _ff_sum_8_mults_arr(t0,f+20,g+19); _ff_add(t1,f[19],g[18]);  _ff_add(o[46],t0,t1);
    _ff_sum_9_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[45],t0,t1);
    _ff_sum_10_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[44],t0,t1);
    _ff_sum_11_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[43],t0,t1);
    _ff_sum_12_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[42],t0,t1);
    _ff_sum_13_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[41],t0,t1);
    _ff_sum_14_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[40],t0,t1);
    _ff_sum_15_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[39],t0,t1);
    _ff_sum_16_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[38],t0,t1);
    _ff_sum_17_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[37],t0,t1);
    _ff_sum_18_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[36],t0,t1);
    _ff_sum_19_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[35],t0,t1);
    _ff_sum_20_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[34],t0,t1);
    _ff_sum_21_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[33],t0,t1);
    _ff_sum_22_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[32],t0,t1);
    _ff_sum_23_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[31],t0,t1);
    _ff_sum_24_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[30],t0,t1);
    _ff_sum_25_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[29],t0,t1);
    _ff_sum_26_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[28],t0,t1);
    _ff_sum_27_mults_arr(t0,f+1,g);  _ff_add(o[27],t0,f[0]);
    _ff_sum_27_mults_arr(o[26],f,g);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_56 (ff_t o[56], ff_t r[56])
{
   ff_t f[28], g[28];
    register ff_t t0, t1;

    ff_poly_from_roots_28(f,r);  ff_poly_from_roots_28(g,r+28);     _ff_add(o[55],f[27],g[27]);
    _ff_mult(t0,f[27],g[27]); _ff_add(t1,f[26],g[26]); _ff_add(o[54],t0,t1);
    _ff_sum_2_mults_arr(t0,f+26,g+26); _ff_add(t1,f[25],g[25]);  _ff_add(o[53],t0,t1);
    _ff_sum_3_mults_arr(t0,f+25,g+25); _ff_add(t1,f[24],g[24]);  _ff_add(o[52],t0,t1);
    _ff_sum_4_mults_arr(t0,f+24,g+24); _ff_add(t1,f[23],g[23]);  _ff_add(o[51],t0,t1);
    _ff_sum_5_mults_arr(t0,f+23,g+23); _ff_add(t1,f[22],g[22]);  _ff_add(o[50],t0,t1);
    _ff_sum_6_mults_arr(t0,f+22,g+22); _ff_add(t1,f[21],g[21]);  _ff_add(o[49],t0,t1);
    _ff_sum_7_mults_arr(t0,f+21,g+21); _ff_add(t1,f[20],g[20]);  _ff_add(o[48],t0,t1);
    _ff_sum_8_mults_arr(t0,f+20,g+20); _ff_add(t1,f[19],g[19]);  _ff_add(o[47],t0,t1);
    _ff_sum_9_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[46],t0,t1);
    _ff_sum_10_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[45],t0,t1);
    _ff_sum_11_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[44],t0,t1);
    _ff_sum_12_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[43],t0,t1);
    _ff_sum_13_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[42],t0,t1);
    _ff_sum_14_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[41],t0,t1);
    _ff_sum_15_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[40],t0,t1);
    _ff_sum_16_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[39],t0,t1);
    _ff_sum_17_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[38],t0,t1);
    _ff_sum_18_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[37],t0,t1);
    _ff_sum_19_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[36],t0,t1);
    _ff_sum_20_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[35],t0,t1);
    _ff_sum_21_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[34],t0,t1);
    _ff_sum_22_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[33],t0,t1);
    _ff_sum_23_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[32],t0,t1);
    _ff_sum_24_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[31],t0,t1);
    _ff_sum_25_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[30],t0,t1);
    _ff_sum_26_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[29],t0,t1);
    _ff_sum_27_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[28],t0,t1);
    _ff_sum_28_mults_arr(o[27],f,g);
    _ff_sum_27_mults_arr(o[26],f,g);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_57 (ff_t o[57], ff_t r[57])
{
   ff_t f[29], g[28];
    register ff_t t0, t1;

    ff_poly_from_roots_29(f,r);  ff_poly_from_roots_28(g,r+29);     _ff_add(o[56],f[28],g[27]);
    _ff_mult(t0,f[28],g[27]); _ff_add(t1,f[27],g[26]); _ff_add(o[55],t0,t1);
    _ff_sum_2_mults_arr(t0,f+27,g+26); _ff_add(t1,f[26],g[25]);  _ff_add(o[54],t0,t1);
    _ff_sum_3_mults_arr(t0,f+26,g+25); _ff_add(t1,f[25],g[24]);  _ff_add(o[53],t0,t1);
    _ff_sum_4_mults_arr(t0,f+25,g+24); _ff_add(t1,f[24],g[23]);  _ff_add(o[52],t0,t1);
    _ff_sum_5_mults_arr(t0,f+24,g+23); _ff_add(t1,f[23],g[22]);  _ff_add(o[51],t0,t1);
    _ff_sum_6_mults_arr(t0,f+23,g+22); _ff_add(t1,f[22],g[21]);  _ff_add(o[50],t0,t1);
    _ff_sum_7_mults_arr(t0,f+22,g+21); _ff_add(t1,f[21],g[20]);  _ff_add(o[49],t0,t1);
    _ff_sum_8_mults_arr(t0,f+21,g+20); _ff_add(t1,f[20],g[19]);  _ff_add(o[48],t0,t1);
    _ff_sum_9_mults_arr(t0,f+20,g+19); _ff_add(t1,f[19],g[18]);  _ff_add(o[47],t0,t1);
    _ff_sum_10_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[46],t0,t1);
    _ff_sum_11_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[45],t0,t1);
    _ff_sum_12_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[44],t0,t1);
    _ff_sum_13_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[43],t0,t1);
    _ff_sum_14_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[42],t0,t1);
    _ff_sum_15_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[41],t0,t1);
    _ff_sum_16_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[40],t0,t1);
    _ff_sum_17_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[39],t0,t1);
    _ff_sum_18_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[38],t0,t1);
    _ff_sum_19_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[37],t0,t1);
    _ff_sum_20_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[36],t0,t1);
    _ff_sum_21_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[35],t0,t1);
    _ff_sum_22_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[34],t0,t1);
    _ff_sum_23_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[33],t0,t1);
    _ff_sum_24_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[32],t0,t1);
    _ff_sum_25_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[31],t0,t1);
    _ff_sum_26_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[30],t0,t1);
    _ff_sum_27_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[29],t0,t1);
    _ff_sum_28_mults_arr(t0,f+1,g);  _ff_add(o[28],t0,f[0]);
    _ff_sum_28_mults_arr(o[27],f,g);
    _ff_sum_27_mults_arr(o[26],f,g);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_58 (ff_t o[58], ff_t r[58])
{
   ff_t f[29], g[29];
    register ff_t t0, t1;

    ff_poly_from_roots_29(f,r);  ff_poly_from_roots_29(g,r+29);     _ff_add(o[57],f[28],g[28]);
    _ff_mult(t0,f[28],g[28]); _ff_add(t1,f[27],g[27]); _ff_add(o[56],t0,t1);
    _ff_sum_2_mults_arr(t0,f+27,g+27); _ff_add(t1,f[26],g[26]);  _ff_add(o[55],t0,t1);
    _ff_sum_3_mults_arr(t0,f+26,g+26); _ff_add(t1,f[25],g[25]);  _ff_add(o[54],t0,t1);
    _ff_sum_4_mults_arr(t0,f+25,g+25); _ff_add(t1,f[24],g[24]);  _ff_add(o[53],t0,t1);
    _ff_sum_5_mults_arr(t0,f+24,g+24); _ff_add(t1,f[23],g[23]);  _ff_add(o[52],t0,t1);
    _ff_sum_6_mults_arr(t0,f+23,g+23); _ff_add(t1,f[22],g[22]);  _ff_add(o[51],t0,t1);
    _ff_sum_7_mults_arr(t0,f+22,g+22); _ff_add(t1,f[21],g[21]);  _ff_add(o[50],t0,t1);
    _ff_sum_8_mults_arr(t0,f+21,g+21); _ff_add(t1,f[20],g[20]);  _ff_add(o[49],t0,t1);
    _ff_sum_9_mults_arr(t0,f+20,g+20); _ff_add(t1,f[19],g[19]);  _ff_add(o[48],t0,t1);
    _ff_sum_10_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[47],t0,t1);
    _ff_sum_11_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[46],t0,t1);
    _ff_sum_12_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[45],t0,t1);
    _ff_sum_13_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[44],t0,t1);
    _ff_sum_14_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[43],t0,t1);
    _ff_sum_15_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[42],t0,t1);
    _ff_sum_16_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[41],t0,t1);
    _ff_sum_17_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[40],t0,t1);
    _ff_sum_18_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[39],t0,t1);
    _ff_sum_19_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[38],t0,t1);
    _ff_sum_20_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[37],t0,t1);
    _ff_sum_21_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[36],t0,t1);
    _ff_sum_22_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[35],t0,t1);
    _ff_sum_23_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[34],t0,t1);
    _ff_sum_24_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[33],t0,t1);
    _ff_sum_25_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[32],t0,t1);
    _ff_sum_26_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[31],t0,t1);
    _ff_sum_27_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[30],t0,t1);
    _ff_sum_28_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[29],t0,t1);
    _ff_sum_29_mults_arr(o[28],f,g);
    _ff_sum_28_mults_arr(o[27],f,g);
    _ff_sum_27_mults_arr(o[26],f,g);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_59 (ff_t o[59], ff_t r[59])
{
   ff_t f[30], g[29];
    register ff_t t0, t1;

    ff_poly_from_roots_30(f,r);  ff_poly_from_roots_29(g,r+30);     _ff_add(o[58],f[29],g[28]);
    _ff_mult(t0,f[29],g[28]); _ff_add(t1,f[28],g[27]); _ff_add(o[57],t0,t1);
    _ff_sum_2_mults_arr(t0,f+28,g+27); _ff_add(t1,f[27],g[26]);  _ff_add(o[56],t0,t1);
    _ff_sum_3_mults_arr(t0,f+27,g+26); _ff_add(t1,f[26],g[25]);  _ff_add(o[55],t0,t1);
    _ff_sum_4_mults_arr(t0,f+26,g+25); _ff_add(t1,f[25],g[24]);  _ff_add(o[54],t0,t1);
    _ff_sum_5_mults_arr(t0,f+25,g+24); _ff_add(t1,f[24],g[23]);  _ff_add(o[53],t0,t1);
    _ff_sum_6_mults_arr(t0,f+24,g+23); _ff_add(t1,f[23],g[22]);  _ff_add(o[52],t0,t1);
    _ff_sum_7_mults_arr(t0,f+23,g+22); _ff_add(t1,f[22],g[21]);  _ff_add(o[51],t0,t1);
    _ff_sum_8_mults_arr(t0,f+22,g+21); _ff_add(t1,f[21],g[20]);  _ff_add(o[50],t0,t1);
    _ff_sum_9_mults_arr(t0,f+21,g+20); _ff_add(t1,f[20],g[19]);  _ff_add(o[49],t0,t1);
    _ff_sum_10_mults_arr(t0,f+20,g+19); _ff_add(t1,f[19],g[18]);  _ff_add(o[48],t0,t1);
    _ff_sum_11_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[47],t0,t1);
    _ff_sum_12_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[46],t0,t1);
    _ff_sum_13_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[45],t0,t1);
    _ff_sum_14_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[44],t0,t1);
    _ff_sum_15_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[43],t0,t1);
    _ff_sum_16_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[42],t0,t1);
    _ff_sum_17_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[41],t0,t1);
    _ff_sum_18_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[40],t0,t1);
    _ff_sum_19_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[39],t0,t1);
    _ff_sum_20_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[38],t0,t1);
    _ff_sum_21_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[37],t0,t1);
    _ff_sum_22_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[36],t0,t1);
    _ff_sum_23_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[35],t0,t1);
    _ff_sum_24_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[34],t0,t1);
    _ff_sum_25_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[33],t0,t1);
    _ff_sum_26_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[32],t0,t1);
    _ff_sum_27_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[31],t0,t1);
    _ff_sum_28_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[30],t0,t1);
    _ff_sum_29_mults_arr(t0,f+1,g);  _ff_add(o[29],t0,f[0]);
    _ff_sum_29_mults_arr(o[28],f,g);
    _ff_sum_28_mults_arr(o[27],f,g);
    _ff_sum_27_mults_arr(o[26],f,g);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_60 (ff_t o[60], ff_t r[60])
{
   ff_t f[30], g[30];
    register ff_t t0, t1;

    ff_poly_from_roots_30(f,r);  ff_poly_from_roots_30(g,r+30);     _ff_add(o[59],f[29],g[29]);
    _ff_mult(t0,f[29],g[29]); _ff_add(t1,f[28],g[28]); _ff_add(o[58],t0,t1);
    _ff_sum_2_mults_arr(t0,f+28,g+28); _ff_add(t1,f[27],g[27]);  _ff_add(o[57],t0,t1);
    _ff_sum_3_mults_arr(t0,f+27,g+27); _ff_add(t1,f[26],g[26]);  _ff_add(o[56],t0,t1);
    _ff_sum_4_mults_arr(t0,f+26,g+26); _ff_add(t1,f[25],g[25]);  _ff_add(o[55],t0,t1);
    _ff_sum_5_mults_arr(t0,f+25,g+25); _ff_add(t1,f[24],g[24]);  _ff_add(o[54],t0,t1);
    _ff_sum_6_mults_arr(t0,f+24,g+24); _ff_add(t1,f[23],g[23]);  _ff_add(o[53],t0,t1);
    _ff_sum_7_mults_arr(t0,f+23,g+23); _ff_add(t1,f[22],g[22]);  _ff_add(o[52],t0,t1);
    _ff_sum_8_mults_arr(t0,f+22,g+22); _ff_add(t1,f[21],g[21]);  _ff_add(o[51],t0,t1);
    _ff_sum_9_mults_arr(t0,f+21,g+21); _ff_add(t1,f[20],g[20]);  _ff_add(o[50],t0,t1);
    _ff_sum_10_mults_arr(t0,f+20,g+20); _ff_add(t1,f[19],g[19]);  _ff_add(o[49],t0,t1);
    _ff_sum_11_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[48],t0,t1);
    _ff_sum_12_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[47],t0,t1);
    _ff_sum_13_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[46],t0,t1);
    _ff_sum_14_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[45],t0,t1);
    _ff_sum_15_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[44],t0,t1);
    _ff_sum_16_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[43],t0,t1);
    _ff_sum_17_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[42],t0,t1);
    _ff_sum_18_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[41],t0,t1);
    _ff_sum_19_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[40],t0,t1);
    _ff_sum_20_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[39],t0,t1);
    _ff_sum_21_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[38],t0,t1);
    _ff_sum_22_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[37],t0,t1);
    _ff_sum_23_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[36],t0,t1);
    _ff_sum_24_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[35],t0,t1);
    _ff_sum_25_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[34],t0,t1);
    _ff_sum_26_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[33],t0,t1);
    _ff_sum_27_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[32],t0,t1);
    _ff_sum_28_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[31],t0,t1);
    _ff_sum_29_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[30],t0,t1);
    _ff_sum_30_mults_arr(o[29],f,g);
    _ff_sum_29_mults_arr(o[28],f,g);
    _ff_sum_28_mults_arr(o[27],f,g);
    _ff_sum_27_mults_arr(o[26],f,g);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_61 (ff_t o[61], ff_t r[61])
{
   ff_t f[31], g[30];
    register ff_t t0, t1;

    ff_poly_from_roots_31(f,r);  ff_poly_from_roots_30(g,r+31);     _ff_add(o[60],f[30],g[29]);
    _ff_mult(t0,f[30],g[29]); _ff_add(t1,f[29],g[28]); _ff_add(o[59],t0,t1);
    _ff_sum_2_mults_arr(t0,f+29,g+28); _ff_add(t1,f[28],g[27]);  _ff_add(o[58],t0,t1);
    _ff_sum_3_mults_arr(t0,f+28,g+27); _ff_add(t1,f[27],g[26]);  _ff_add(o[57],t0,t1);
    _ff_sum_4_mults_arr(t0,f+27,g+26); _ff_add(t1,f[26],g[25]);  _ff_add(o[56],t0,t1);
    _ff_sum_5_mults_arr(t0,f+26,g+25); _ff_add(t1,f[25],g[24]);  _ff_add(o[55],t0,t1);
    _ff_sum_6_mults_arr(t0,f+25,g+24); _ff_add(t1,f[24],g[23]);  _ff_add(o[54],t0,t1);
    _ff_sum_7_mults_arr(t0,f+24,g+23); _ff_add(t1,f[23],g[22]);  _ff_add(o[53],t0,t1);
    _ff_sum_8_mults_arr(t0,f+23,g+22); _ff_add(t1,f[22],g[21]);  _ff_add(o[52],t0,t1);
    _ff_sum_9_mults_arr(t0,f+22,g+21); _ff_add(t1,f[21],g[20]);  _ff_add(o[51],t0,t1);
    _ff_sum_10_mults_arr(t0,f+21,g+20); _ff_add(t1,f[20],g[19]);  _ff_add(o[50],t0,t1);
    _ff_sum_11_mults_arr(t0,f+20,g+19); _ff_add(t1,f[19],g[18]);  _ff_add(o[49],t0,t1);
    _ff_sum_12_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[48],t0,t1);
    _ff_sum_13_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[47],t0,t1);
    _ff_sum_14_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[46],t0,t1);
    _ff_sum_15_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[45],t0,t1);
    _ff_sum_16_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[44],t0,t1);
    _ff_sum_17_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[43],t0,t1);
    _ff_sum_18_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[42],t0,t1);
    _ff_sum_19_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[41],t0,t1);
    _ff_sum_20_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[40],t0,t1);
    _ff_sum_21_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[39],t0,t1);
    _ff_sum_22_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[38],t0,t1);
    _ff_sum_23_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[37],t0,t1);
    _ff_sum_24_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[36],t0,t1);
    _ff_sum_25_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[35],t0,t1);
    _ff_sum_26_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[34],t0,t1);
    _ff_sum_27_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[33],t0,t1);
    _ff_sum_28_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[32],t0,t1);
    _ff_sum_29_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[31],t0,t1);
    _ff_sum_30_mults_arr(t0,f+1,g);  _ff_add(o[30],t0,f[0]);
    _ff_sum_30_mults_arr(o[29],f,g);
    _ff_sum_29_mults_arr(o[28],f,g);
    _ff_sum_28_mults_arr(o[27],f,g);
    _ff_sum_27_mults_arr(o[26],f,g);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_62 (ff_t o[62], ff_t r[62])
{
   ff_t f[31], g[31];
    register ff_t t0, t1;

    ff_poly_from_roots_31(f,r);  ff_poly_from_roots_31(g,r+31);     _ff_add(o[61],f[30],g[30]);
    _ff_mult(t0,f[30],g[30]); _ff_add(t1,f[29],g[29]); _ff_add(o[60],t0,t1);
    _ff_sum_2_mults_arr(t0,f+29,g+29); _ff_add(t1,f[28],g[28]);  _ff_add(o[59],t0,t1);
    _ff_sum_3_mults_arr(t0,f+28,g+28); _ff_add(t1,f[27],g[27]);  _ff_add(o[58],t0,t1);
    _ff_sum_4_mults_arr(t0,f+27,g+27); _ff_add(t1,f[26],g[26]);  _ff_add(o[57],t0,t1);
    _ff_sum_5_mults_arr(t0,f+26,g+26); _ff_add(t1,f[25],g[25]);  _ff_add(o[56],t0,t1);
    _ff_sum_6_mults_arr(t0,f+25,g+25); _ff_add(t1,f[24],g[24]);  _ff_add(o[55],t0,t1);
    _ff_sum_7_mults_arr(t0,f+24,g+24); _ff_add(t1,f[23],g[23]);  _ff_add(o[54],t0,t1);
    _ff_sum_8_mults_arr(t0,f+23,g+23); _ff_add(t1,f[22],g[22]);  _ff_add(o[53],t0,t1);
    _ff_sum_9_mults_arr(t0,f+22,g+22); _ff_add(t1,f[21],g[21]);  _ff_add(o[52],t0,t1);
    _ff_sum_10_mults_arr(t0,f+21,g+21); _ff_add(t1,f[20],g[20]);  _ff_add(o[51],t0,t1);
    _ff_sum_11_mults_arr(t0,f+20,g+20); _ff_add(t1,f[19],g[19]);  _ff_add(o[50],t0,t1);
    _ff_sum_12_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[49],t0,t1);
    _ff_sum_13_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[48],t0,t1);
    _ff_sum_14_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[47],t0,t1);
    _ff_sum_15_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[46],t0,t1);
    _ff_sum_16_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[45],t0,t1);
    _ff_sum_17_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[44],t0,t1);
    _ff_sum_18_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[43],t0,t1);
    _ff_sum_19_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[42],t0,t1);
    _ff_sum_20_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[41],t0,t1);
    _ff_sum_21_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[40],t0,t1);
    _ff_sum_22_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[39],t0,t1);
    _ff_sum_23_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[38],t0,t1);
    _ff_sum_24_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[37],t0,t1);
    _ff_sum_25_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[36],t0,t1);
    _ff_sum_26_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[35],t0,t1);
    _ff_sum_27_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[34],t0,t1);
    _ff_sum_28_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[33],t0,t1);
    _ff_sum_29_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[32],t0,t1);
    _ff_sum_30_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[31],t0,t1);
    _ff_sum_31_mults_arr(o[30],f,g);
    _ff_sum_30_mults_arr(o[29],f,g);
    _ff_sum_29_mults_arr(o[28],f,g);
    _ff_sum_28_mults_arr(o[27],f,g);
    _ff_sum_27_mults_arr(o[26],f,g);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

static inline void ff_poly_from_roots_63 (ff_t o[63], ff_t r[63])
{
   ff_t f[32], g[31];
    register ff_t t0, t1;

    ff_poly_from_roots_32(f,r);  ff_poly_from_roots_31(g,r+32);     _ff_add(o[62],f[31],g[30]);
    _ff_mult(t0,f[31],g[30]); _ff_add(t1,f[30],g[29]); _ff_add(o[61],t0,t1);
    _ff_sum_2_mults_arr(t0,f+30,g+29); _ff_add(t1,f[29],g[28]);  _ff_add(o[60],t0,t1);
    _ff_sum_3_mults_arr(t0,f+29,g+28); _ff_add(t1,f[28],g[27]);  _ff_add(o[59],t0,t1);
    _ff_sum_4_mults_arr(t0,f+28,g+27); _ff_add(t1,f[27],g[26]);  _ff_add(o[58],t0,t1);
    _ff_sum_5_mults_arr(t0,f+27,g+26); _ff_add(t1,f[26],g[25]);  _ff_add(o[57],t0,t1);
    _ff_sum_6_mults_arr(t0,f+26,g+25); _ff_add(t1,f[25],g[24]);  _ff_add(o[56],t0,t1);
    _ff_sum_7_mults_arr(t0,f+25,g+24); _ff_add(t1,f[24],g[23]);  _ff_add(o[55],t0,t1);
    _ff_sum_8_mults_arr(t0,f+24,g+23); _ff_add(t1,f[23],g[22]);  _ff_add(o[54],t0,t1);
    _ff_sum_9_mults_arr(t0,f+23,g+22); _ff_add(t1,f[22],g[21]);  _ff_add(o[53],t0,t1);
    _ff_sum_10_mults_arr(t0,f+22,g+21); _ff_add(t1,f[21],g[20]);  _ff_add(o[52],t0,t1);
    _ff_sum_11_mults_arr(t0,f+21,g+20); _ff_add(t1,f[20],g[19]);  _ff_add(o[51],t0,t1);
    _ff_sum_12_mults_arr(t0,f+20,g+19); _ff_add(t1,f[19],g[18]);  _ff_add(o[50],t0,t1);
    _ff_sum_13_mults_arr(t0,f+19,g+18); _ff_add(t1,f[18],g[17]);  _ff_add(o[49],t0,t1);
    _ff_sum_14_mults_arr(t0,f+18,g+17); _ff_add(t1,f[17],g[16]);  _ff_add(o[48],t0,t1);
    _ff_sum_15_mults_arr(t0,f+17,g+16); _ff_add(t1,f[16],g[15]);  _ff_add(o[47],t0,t1);
    _ff_sum_16_mults_arr(t0,f+16,g+15); _ff_add(t1,f[15],g[14]);  _ff_add(o[46],t0,t1);
    _ff_sum_17_mults_arr(t0,f+15,g+14); _ff_add(t1,f[14],g[13]);  _ff_add(o[45],t0,t1);
    _ff_sum_18_mults_arr(t0,f+14,g+13); _ff_add(t1,f[13],g[12]);  _ff_add(o[44],t0,t1);
    _ff_sum_19_mults_arr(t0,f+13,g+12); _ff_add(t1,f[12],g[11]);  _ff_add(o[43],t0,t1);
    _ff_sum_20_mults_arr(t0,f+12,g+11); _ff_add(t1,f[11],g[10]);  _ff_add(o[42],t0,t1);
    _ff_sum_21_mults_arr(t0,f+11,g+10); _ff_add(t1,f[10],g[9]);  _ff_add(o[41],t0,t1);
    _ff_sum_22_mults_arr(t0,f+10,g+9); _ff_add(t1,f[9],g[8]);  _ff_add(o[40],t0,t1);
    _ff_sum_23_mults_arr(t0,f+9,g+8); _ff_add(t1,f[8],g[7]);  _ff_add(o[39],t0,t1);
    _ff_sum_24_mults_arr(t0,f+8,g+7); _ff_add(t1,f[7],g[6]);  _ff_add(o[38],t0,t1);
    _ff_sum_25_mults_arr(t0,f+7,g+6); _ff_add(t1,f[6],g[5]);  _ff_add(o[37],t0,t1);
    _ff_sum_26_mults_arr(t0,f+6,g+5); _ff_add(t1,f[5],g[4]);  _ff_add(o[36],t0,t1);
    _ff_sum_27_mults_arr(t0,f+5,g+4); _ff_add(t1,f[4],g[3]);  _ff_add(o[35],t0,t1);
    _ff_sum_28_mults_arr(t0,f+4,g+3); _ff_add(t1,f[3],g[2]);  _ff_add(o[34],t0,t1);
    _ff_sum_29_mults_arr(t0,f+3,g+2); _ff_add(t1,f[2],g[1]);  _ff_add(o[33],t0,t1);
    _ff_sum_30_mults_arr(t0,f+2,g+1); _ff_add(t1,f[1],g[0]);  _ff_add(o[32],t0,t1);
    _ff_sum_31_mults_arr(t0,f+1,g);  _ff_add(o[31],t0,f[0]);
    _ff_sum_31_mults_arr(o[30],f,g);
    _ff_sum_30_mults_arr(o[29],f,g);
    _ff_sum_29_mults_arr(o[28],f,g);
    _ff_sum_28_mults_arr(o[27],f,g);
    _ff_sum_27_mults_arr(o[26],f,g);
    _ff_sum_26_mults_arr(o[25],f,g);
    _ff_sum_25_mults_arr(o[24],f,g);
    _ff_sum_24_mults_arr(o[23],f,g);
    _ff_sum_23_mults_arr(o[22],f,g);
    _ff_sum_22_mults_arr(o[21],f,g);
    _ff_sum_21_mults_arr(o[20],f,g);
    _ff_sum_20_mults_arr(o[19],f,g);
    _ff_sum_19_mults_arr(o[18],f,g);
    _ff_sum_18_mults_arr(o[17],f,g);
    _ff_sum_17_mults_arr(o[16],f,g);
    _ff_sum_16_mults_arr(o[15],f,g);
    _ff_sum_15_mults_arr(o[14],f,g);
    _ff_sum_14_mults_arr(o[13],f,g);
    _ff_sum_13_mults_arr(o[12],f,g);
    _ff_sum_12_mults_arr(o[11],f,g);
    _ff_sum_11_mults_arr(o[10],f,g);
    _ff_sum_10_mults_arr(o[9],f,g);
    _ff_sum_9_mults_arr(o[8],f,g);
    _ff_sum_8_mults_arr(o[7],f,g);
    _ff_sum_7_mults_arr(o[6],f,g);
    _ff_sum_6_mults_arr(o[5],f,g);
    _ff_sum_5_mults_arr(o[4],f,g);
    _ff_sum_4_mults_arr(o[3],f,g);
    _ff_sum_3_mults_arr(o[2],f,g);
    _ff_sum_2_mults_arr(o[1],f,g);
    _ff_mult(o[0],f[0],g[0]);
}

#ifdef __cplusplus
}
#endif

#endif
