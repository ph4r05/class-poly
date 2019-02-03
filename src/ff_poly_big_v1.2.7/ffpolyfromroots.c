#include <stdio.h>
#include <stdlib.h>
#include "ffpolyfromroots.h"
#include "cstd.h"

/*
    Copyright 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

void ff_poly_from_roots_64 (ff_t o[64], ff_t r[64])
{
   ff_t f[32], g[32];
    register ff_t t0, t1;

    ff_poly_from_roots_32(f,r);  ff_poly_from_roots_32(g,r+32);     _ff_add(o[63],f[31],g[31]);
    _ff_mult(t0,f[31],g[31]); _ff_add(t1,f[30],g[30]); _ff_add(o[62],t0,t1);
    _ff_sum_2_mults_arr(t0,f+30,g+30); _ff_add(t1,f[29],g[29]);  _ff_add(o[61],t0,t1);
    _ff_sum_3_mults_arr(t0,f+29,g+29); _ff_add(t1,f[28],g[28]);  _ff_add(o[60],t0,t1);
    _ff_sum_4_mults_arr(t0,f+28,g+28); _ff_add(t1,f[27],g[27]);  _ff_add(o[59],t0,t1);
    _ff_sum_5_mults_arr(t0,f+27,g+27); _ff_add(t1,f[26],g[26]);  _ff_add(o[58],t0,t1);
    _ff_sum_6_mults_arr(t0,f+26,g+26); _ff_add(t1,f[25],g[25]);  _ff_add(o[57],t0,t1);
    _ff_sum_7_mults_arr(t0,f+25,g+25); _ff_add(t1,f[24],g[24]);  _ff_add(o[56],t0,t1);
    _ff_sum_8_mults_arr(t0,f+24,g+24); _ff_add(t1,f[23],g[23]);  _ff_add(o[55],t0,t1);
    _ff_sum_9_mults_arr(t0,f+23,g+23); _ff_add(t1,f[22],g[22]);  _ff_add(o[54],t0,t1);
    _ff_sum_10_mults_arr(t0,f+22,g+22); _ff_add(t1,f[21],g[21]);  _ff_add(o[53],t0,t1);
    _ff_sum_11_mults_arr(t0,f+21,g+21); _ff_add(t1,f[20],g[20]);  _ff_add(o[52],t0,t1);
    _ff_sum_12_mults_arr(t0,f+20,g+20); _ff_add(t1,f[19],g[19]);  _ff_add(o[51],t0,t1);
    _ff_sum_13_mults_arr(t0,f+19,g+19); _ff_add(t1,f[18],g[18]);  _ff_add(o[50],t0,t1);
    _ff_sum_14_mults_arr(t0,f+18,g+18); _ff_add(t1,f[17],g[17]);  _ff_add(o[49],t0,t1);
    _ff_sum_15_mults_arr(t0,f+17,g+17); _ff_add(t1,f[16],g[16]);  _ff_add(o[48],t0,t1);
    _ff_sum_16_mults_arr(t0,f+16,g+16); _ff_add(t1,f[15],g[15]);  _ff_add(o[47],t0,t1);
    _ff_sum_17_mults_arr(t0,f+15,g+15); _ff_add(t1,f[14],g[14]);  _ff_add(o[46],t0,t1);
    _ff_sum_18_mults_arr(t0,f+14,g+14); _ff_add(t1,f[13],g[13]);  _ff_add(o[45],t0,t1);
    _ff_sum_19_mults_arr(t0,f+13,g+13); _ff_add(t1,f[12],g[12]);  _ff_add(o[44],t0,t1);
    _ff_sum_20_mults_arr(t0,f+12,g+12); _ff_add(t1,f[11],g[11]);  _ff_add(o[43],t0,t1);
    _ff_sum_21_mults_arr(t0,f+11,g+11); _ff_add(t1,f[10],g[10]);  _ff_add(o[42],t0,t1);
    _ff_sum_22_mults_arr(t0,f+10,g+10); _ff_add(t1,f[9],g[9]);  _ff_add(o[41],t0,t1);
    _ff_sum_23_mults_arr(t0,f+9,g+9); _ff_add(t1,f[8],g[8]);  _ff_add(o[40],t0,t1);
    _ff_sum_24_mults_arr(t0,f+8,g+8); _ff_add(t1,f[7],g[7]);  _ff_add(o[39],t0,t1);
    _ff_sum_25_mults_arr(t0,f+7,g+7); _ff_add(t1,f[6],g[6]);  _ff_add(o[38],t0,t1);
    _ff_sum_26_mults_arr(t0,f+6,g+6); _ff_add(t1,f[5],g[5]);  _ff_add(o[37],t0,t1);
    _ff_sum_27_mults_arr(t0,f+5,g+5); _ff_add(t1,f[4],g[4]);  _ff_add(o[36],t0,t1);
    _ff_sum_28_mults_arr(t0,f+4,g+4); _ff_add(t1,f[3],g[3]);  _ff_add(o[35],t0,t1);
    _ff_sum_29_mults_arr(t0,f+3,g+3); _ff_add(t1,f[2],g[2]);  _ff_add(o[34],t0,t1);
    _ff_sum_30_mults_arr(t0,f+2,g+2); _ff_add(t1,f[1],g[1]);  _ff_add(o[33],t0,t1);
    _ff_sum_31_mults_arr(t0,f+1,g+1); _ff_add(t1,f[0],g[0]);  _ff_add(o[32],t0,t1);
    _ff_sum_32_mults_arr(o[31],f,g);
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

// this is slightly slower than using ff_poly_mks_from_roots
void ff_poly_from_roots_128 (ff_t o[], ff_t r[])
{
	ff_t f[64], g[64];
	register ff_t t1;
	register int i;
	
	ff_poly_from_roots_64(f,r);  ff_poly_from_roots_64(g,r+64);
	_ff_add(o[127],f[63],g[63]);
	for ( i = 1 ; i < 64 ; i++ ) { ff_conv(o+127-i,f+64-i,g+64-i,i); _ff_add(t1,f[63-i],g[63-i]); _ff_addto(o[127-i],t1); }
	for ( ; i ; i-- ) ff_conv(o+i-1,f,g,i);
}

// this is slower than using ff_poly_mks_from_roots
void ff_poly_from_roots_256 (ff_t o[], ff_t r[])
{
	ff_t f[128], g[128];
	register ff_t t1;
	register int i;
	
	ff_poly_from_roots_128(f,r);  ff_poly_from_roots_128(g,r+128);
	_ff_add(o[255],f[127],g[127]);
	for ( i = 1 ; i < 128 ; i++ ) { ff_conv(o+255-i,f+128-i,g+128-i,i); _ff_add(t1,f[127-i],g[127-i]); _ff_addto(o[255-i],t1); }
	for ( ; i ; i-- ) ff_conv(o+i-1,f,g,i);
}

void ff_poly_from_roots_tiny (ff_t o[], ff_t r[], int d)
{
	switch (d) {
	case 1:  _ff_neg(o[0],r[0]); break;
	case 2:  ff_poly_from_roots_2(o,r); break;
	case 3:  ff_poly_from_roots_3(o,r); break;
	case 4:  ff_poly_from_roots_4(o,r); break;
	case 5:  ff_poly_from_roots_5(o,r); break;
	case 6:  ff_poly_from_roots_6(o,r); break;
	case 7:  ff_poly_from_roots_7(o,r); break;
	case 8:  ff_poly_from_roots_8(o,r); break;
	case 9:  ff_poly_from_roots_9(o,r); break;
	case 10:  ff_poly_from_roots_10(o,r); break;
	case 11:  ff_poly_from_roots_11(o,r); break;
	case 12:  ff_poly_from_roots_12(o,r); break;
	case 13:  ff_poly_from_roots_13(o,r); break;
	case 14:  ff_poly_from_roots_14(o,r); break;
	case 15:  ff_poly_from_roots_15(o,r); break;
	case 16:  ff_poly_from_roots_16(o,r); break;
	case 17:  ff_poly_from_roots_17(o,r); break;
	case 18:  ff_poly_from_roots_18(o,r); break;
	case 19:  ff_poly_from_roots_19(o,r); break;
	case 20:  ff_poly_from_roots_20(o,r); break;
	case 21:  ff_poly_from_roots_21(o,r); break;
	case 22:  ff_poly_from_roots_22(o,r); break;
	case 23:  ff_poly_from_roots_23(o,r); break;
	case 24:  ff_poly_from_roots_24(o,r); break;
	case 25:  ff_poly_from_roots_25(o,r); break;
	case 26:  ff_poly_from_roots_26(o,r); break;
	case 27:  ff_poly_from_roots_27(o,r); break;
	case 28:  ff_poly_from_roots_28(o,r); break;
	case 29:  ff_poly_from_roots_29(o,r); break;
	case 30:  ff_poly_from_roots_30(o,r); break;
	case 31:  ff_poly_from_roots_31(o,r); break;
	case 32:  ff_poly_from_roots_32(o,r); break;
	default:
		err_printf ("d=%d > 32 in ff_poly_from_roots_tiny\n", d);
		abort();
	}
	_ff_set_one(o[d]);
	return;
}

void ff_poly_from_roots_small (ff_t o[], ff_t r[], int d)
{
	switch (d) {
	case 32:  ff_poly_from_roots_32(o,r); break;
	case 33:  ff_poly_from_roots_33(o,r); break;
	case 34:  ff_poly_from_roots_34(o,r); break;
	case 35:  ff_poly_from_roots_35(o,r); break;
	case 36:  ff_poly_from_roots_36(o,r); break;
	case 37:  ff_poly_from_roots_37(o,r); break;
	case 38:  ff_poly_from_roots_38(o,r); break;
	case 39:  ff_poly_from_roots_39(o,r); break;
	case 40:  ff_poly_from_roots_40(o,r); break;
	case 41:  ff_poly_from_roots_41(o,r); break;
	case 42:  ff_poly_from_roots_42(o,r); break;
	case 43:  ff_poly_from_roots_43(o,r); break;
	case 44:  ff_poly_from_roots_44(o,r); break;
	case 45:  ff_poly_from_roots_45(o,r); break;
	case 46:  ff_poly_from_roots_46(o,r); break;
	case 47:  ff_poly_from_roots_47(o,r); break;
	case 48:  ff_poly_from_roots_48(o,r); break;
	case 49:  ff_poly_from_roots_49(o,r); break;
	case 50:  ff_poly_from_roots_50(o,r); break;
	case 51:  ff_poly_from_roots_51(o,r); break;
	case 52:  ff_poly_from_roots_52(o,r); break;
	case 53:  ff_poly_from_roots_53(o,r); break;
	case 54:  ff_poly_from_roots_54(o,r); break;
	case 55:  ff_poly_from_roots_55(o,r); break;
	case 56:  ff_poly_from_roots_56(o,r); break;
	case 57:  ff_poly_from_roots_57(o,r); break;
	case 58:  ff_poly_from_roots_58(o,r); break;
	case 59:  ff_poly_from_roots_59(o,r); break;
	case 60:  ff_poly_from_roots_60(o,r); break;
	case 61:  ff_poly_from_roots_61(o,r); break;
	case 62:  ff_poly_from_roots_62(o,r); break;
	case 63:  ff_poly_from_roots_63(o,r); break;
	case 64:  ff_poly_from_roots_64(o,r); break;
	default:
		if ( d < 32 ) { ff_poly_from_roots_tiny(o,r,d);  return; }
		err_printf ("d=%d > 64 in ff_poly_from_roots_small\n", d);
		abort();
	}
	_ff_set_one(o[d]);
	return;
}

// this is slightly slower than using ff_poly_mks_from_roots
void ff_poly_from_roots_med (ff_t o[], ff_t r[], int n)
{
	ff_t f[64], g[64];
	register ff_t t1;
	register int i,m;
	
	if ( n <= 64 ) { ff_poly_from_roots_small(o,r,n); return; }
	m = (n+1)/2;
	ff_poly_from_roots_small(f,r,m);  ff_poly_from_roots_small(g,r+m,n-m);
	_ff_set_one(o[n]);
	_ff_add(o[n-1],f[m-1],g[n-m-1]);
	if ( (n&1) ) {
		for ( i = 1 ; i < n-m ; i++ ) { ff_conv(o+n-i-1,f+n-m-i,g+m-i,i); _ff_add(t1,f[n-m-i-1],g[m-i-1]); _ff_addto(o[n-i-1],t1); }
		ff_conv(o+i,f+1,g,i);
		for ( ; i ; i-- ) ff_conv(o+i-1,f,g,i);
	} else {
		for ( i = 1 ; i < m ; i++ ) { ff_conv(o+n-i-1,f+m-i,g+m-i,i); _ff_add(t1,f[m-i-1],g[m-i-1]); _ff_addto(o[n-i-1],t1); }
		for ( ; i ; i-- ) ff_conv(o+i-1,f,g,i);
	}
}
