#ifndef _FF_INCLUDE_
#define _FF_INCLUDE_

/*
    Copyright 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <assert.h>
#include <math.h>
#include <gmp.h>
#include <limits.h>
#include "asm.h"
#include "ntutil.h"
#include "cstd.h"

#ifdef __cplusplus
extern "C" {
#endif

// For a stripped-down, self-contained implementation if 32 or 64 bit single precision
// finite fields using a Montgomery representation.  All you need is ff.h, asm.h, and ff.c

// THIS CODE IS NOT PORTABLE, but it is fast.  It assumes an unsigned is 32 bits and that an
// unsigned long is 64 bits and (most critically) uses  inline assembly directives.  These
// should work on any AMD/Intel compatible instruction set.

// The FF_MONTGOMERY flag should generally be set, the only reason to unset it is for 
// testing or benchmarking.  It can be helpful to find bugs in code that breaks when
// the field representation changes (e.g. when the bit string "000...01" is not the identity element).

#define FF_NO_EXT_SETUP			0				// eliminate need for ffext.c (and ffpoly.c)  means no sqrt/cbrt functions can be called, for testing only!

#ifndef ULONG_BITS
#if ULONG_MAX == 18446744073709551615UL
#define ULONG_BITS 64
#else
#error ff_poly requires an unsigned long to hold (at least) 64 bits
#endif
#endif

#define FF_BIG_P					0				// supports p up to 61 (could be 62) bits, otherwise support is limited to 56 bits (but substantially faster: 20-30% in places)

#define FF_FAST					1				// define to turn off certain error checking
#define FF_SUPERFAST				1				// define to enable inlining that is painfully expensive to compile but generates significantly faster code for root finding with prime degrees between 30 and 50

#define FF_MAX_PARALLEL_INVERTS		32768			// make big enough to not force chunking in hecurve1(and to support modpoly ell up to 65535)
#define FF_FRACTIONS				511
#define FF_WORDS					1				// must be 1
#define FF_HALF_WORDS				2				// 
#define FF_MONTGOMERY				1				// if 0 FF_HALF_WORDS must be 1
#if FF_BIG_P
#define FF_NAIL_BITS				2				// must be at least 2 (need to be able to add three elements without overflow) 
#else
#define FF_NAIL_BITS				7				// can add 63 elements without overflow
#endif
#define FF_BITS					((FF_HALF_WORDS*ULONG_BITS/2) - FF_NAIL_BITS)
#define FF_MONTGOMERY_RBITS		(FF_HALF_WORDS*ULONG_BITS/2)	// note that R=B are identical - only one word is ever used
#define FF_ITAB_SIZE				(3*FF_MONTGOMERY_RBITS+1)

#define FF_MAX_CHAIN_LEN			42
#define FF_BINOMIAL_MAX			16


#if FF_MONTGOMERY && FF_HALF_WORDS == 1
typedef unsigned ff_t;
#else
typedef unsigned long ff_t;
#endif

// temporary value used by macros - could remove by replacing macros with inlines, probably should
extern ff_t _ff_t1;

// The globals below all assume a single modulus environment.  Multiple moduli can be managed by
// copying and resetting these (provided this is reasonably infrequent).  We don't put them into a
// dynamically allocated context to save the cost of dereferencing a pointer every time we want to
// perform a field operation.

extern ff_t _ff_p;
extern ff_t _ff_2p;
extern ff_t _ff_2g;					// generator of the Sylow 2-subgroup, necessarily a quadratic non-residue
extern ff_t _ff_2gi;					// inverse of _ff_2g
extern ff_t _ff_2Sylow_tab[64][2];
extern ff_t _ff_3g;					// generator of the Sylow 3-subgroup, necessarily a cubic non-residue
extern ff_t _ff_half;					// 1/2 = (p+1)/2
extern ff_t _ff_third;				// 1/3 = (p+1)/3 or (2p+1)/3 for p=2mod3 or 1mod3 (resp.)
extern ff_t _ff_fourth;				// 1/2*1/2
extern ff_t _ff_frac[FF_FRACTIONS+1];	// _ff_frac[i] = 1/i for i > 0.  Only set if ff_invert_fractions is called.		
extern ff_t _ff_negone;
extern ff_t _ff_negthree;
extern int _ff_p2_e;
extern unsigned long _ff_p2_m;		// p = 2^e*m+1
extern int _ff_p3_e;
extern unsigned long _ff_p3_m;		// p = 3^e*m+1 when p=1mod3
extern int _ff_p3_m1mod3;
extern ff_t *_ff_mont_itab;
extern ff_t _ff_mont_R;
extern ff_t _ff_mont_R2;
extern unsigned long _ff_mont_pni;
extern int _ff_p1mod3;

extern int _ff_cbrt_setup;
extern ff_t _ff_cbrt_unity;			// MUST CALL ff_cbrt_setup() or ff_cbrt() to initialize.  Set to 1 if p=2mod3

extern ff_t _ff_binomial_tab[((FF_BINOMIAL_MAX+1)*(FF_BINOMIAL_MAX+2))/2];
int _ff_binomial_top;


void ff_setup_ui (unsigned long p);	// note - DOES NOT CHECK PRIMALITY


// WARNING - several of the macros below are decidedly unsafe.  In general, if it starts with an underscore,
// don't use side effects, and be aware that many can't be used in expressions.
// These could be converted to inline functions, but many involve assignments and would require passing
// by reference and then dereferencing a pointer (this is handled by reference types in C++). 
// A smart compiler would produce inline function code that is just as fast as these macros, but we won't rely on this.

#define _ff_swap(a,b)				(_ff_t1=(a),(a)=(b),(b)=_ff_t1)

// Input/Output
#define _ff_sprint(s,x)				sprintf(s,"%lu", _ff_get_ui(x))
#define _ff_set_raw(z,x)				((z) = (x))								// useful for setting random field elements and other tricks.  assumes 0<=x<p
#define _ff_get_raw(z)				(z)
#if FF_MONTGOMERY
#define _ff_set_ui(z,x)				((z) = ff_montgomery1_mult(( (x) < _ff_p ? (x) : ui_mod(x,_ff_p) ), _ff_mont_R2));
#define _ff_rset_ui(z,x)				((z) = ff_montgomery1_mult((x), _ff_mont_R2))// assumes 0<=x<p
#define _ff_get_ui(x)				(ff_montgomery1_mult((x),1UL))
#else
#define _ff_set_ui(z,x)				((z) = (ui_mod(x,_ff_p))
#define _ff_rset_ui(z,x)				((z) = (x))			 					// assumes 0<=x<p
#define _ff_get_ui(x)				(x)
#endif
// signed conversions - you must use these for negative values
#define _ff_get_i(x)					((long)_ff_get_ui(x))
#define _ff_set_i(x,z)				if ( (z) < 0 ) { _ff_set_ui(x,-(z));  ff_negate(x); } else { _ff_set_ui(x,z); }
#define _ff_set_mpz(x,Z)				_ff_set_ui(x, mpz_fdiv_ui(Z,_ff_p))


// Basic operations and 0/1 - IMPORTANT: the binary value 1 is not the field identity in a Montgomery representation (but 0 is zero)
#define _ff_set_null(z)				((z) =_ff_p)		// sets z to a value that does not represent an element of F_p (this works in Montgomery rep also)
#define _ff_null(z)					((z) == _ff_p)
#define _ff_set_zero(z)				((z)=0)
#define _ff_zero(x)					(!(x))
#define _ff_nonzero(x)				(x)
#define _ff_parity(x)					((x)&0x1)
#define _ff_set(z,x)					((z)=(x))
#define _ff_equal(x,y)				((x) == (y))
#define _ff_low_word(x)				(x)
#if FF_MONTGOMERY
#define _ff_set_one(z)				_ff_set(z,_ff_mont_R)
#define _ff_one(z)					_ff_equal(z,_ff_mont_R)
#define _ff_neg_one(z)				 _ff_equal(z,_ff_negone) 
#define _ff_pm_one(z)				( _ff_one(z) || _ff_neg_one(z) )
#else
#define _ff_set_one(z)				((z)=1UL)
#define _ff_one(z)					((z)==1UL)
#define _ff_neg_one(z)				((z)==_ff_p-1UL)
#define _ff_pm_one(z)				( _ff_one(z) || _ff_neg_one(z) )
#endif
// end basic operations and 0/1

// Core arithmetic operations - these may be applied to unreduced values that fit in the word limit
// NOTE: _ff_add(x,y,x) and _ff_sub(x,y,x) are not valid! (but _ff_add(x,x,y) and _ff_sub(x,x,y) are fine)
#define _ff_core_addto(z,x)			((z) += (x))
#define _ff_core_subfrom(z,x)			((z) -= (x))							// assumes z dominates x
#define _ff_core_shiftl(z)				((z) <<= 1)
#define _ff_core_shiftr(z)				((z) >>= 1)
#define _ff_core_ge(x,y)				((x) >= (y))
#define _ff_addto(z,x)				{register ff_t _ff_t; _ff_set(_ff_t,z); _ff_core_addto(_ff_t,x);_ff_core_red(_ff_t); _ff_set(z,_ff_t);}
#define _ff_add(z,x,y)				{_ff_set(z,x);_ff_addto(z,y);}
#define _ff_subfrom(z,x)				{register ff_t _ff_t; _ff_set(_ff_t,z); _ff_core_dom(_ff_t,x); _ff_core_subfrom(_ff_t,x); _ff_set(z,_ff_t);}
#define _ff_sub(z,x,y)				{_ff_set(z,x);_ff_subfrom(z,y);}
#if FF_MONTGOMERY
#define _ff_core_inc(z)				_ff_core_addto(z,_ff_mont_R)
#else
#define _ff_core_inc(z)				((z)++)
#endif

// derived arithmetic operations - note that inputs cannot overlap outputs! use higher level versions if needed
// Note that core_red requires z < 2p
#define _ff_core_red(z)				{ if (_ff_core_ge(z,_ff_p) ) _ff_core_subfrom (z,_ff_p); }
#define _ff_core_dom(z,x)			{ if ( !_ff_core_ge(z,x) ) _ff_core_addto (z,_ff_p); }
#define _ff_neg(z,x)					{ if (_ff_nonzero(x) ) {_ff_set(z,_ff_p); _ff_core_subfrom(z,x); } else { _ff_set_zero(z); } }
#define _ff_x2(z)					_ff_addto(z,z);		// this is much faster than shifting
#define _ff_inc(z)					{_ff_core_inc(z);_ff_core_red(z);}
#if FF_MONTGOMERY
#define _ff_dec(z)					_ff_subfrom((z),_ff_mont_R)
#define _ff_add_one(z,x)				_ff_add(z,x,_ff_mont_R);
#define _ff_sub_one(z,x)				_ff_sub(z,x,_ff_mont_R);
#else
#define _ff_dec(z)					{ if (z) { (z)--; } else { (z) = _ff_p-1; } }
#define _ff_add_one(z,x)				_ff_add(z,x,1);
#define _ff_sub_one(z,x)				_ff_sub(z,x,1);
#endif
// end core arithmetic operations

// Safer versions - overlap ok, but still need to watch side effects and can't use in expressions
#define ff_negate(z)				{ if (_ff_nonzero(z) ) {_ff_set(_ff_t1,_ff_p); _ff_core_subfrom(_ff_t1,z);  _ff_set(z,_ff_t1); } else { _ff_set_zero(z); } }
#define ff_add(z,x,y)				{_ff_set(_ff_t1,x);_ff_addto(_ff_t1,y);_ff_set(z,_ff_t1);}
#define ff_sub(z,x,y)					{_ff_set(_ff_t1,x);_ff_subfrom(_ff_t1,y);_ff_set(z,_ff_t1);}

// higher arithmetic operations
#if FF_MONTGOMERY
#define _ff_mult(z,x,y)				((z) = ff_montgomery1_mult(x,y))
#define _ff_sum_2_mults(z,x0,x1,y0,y1)	((z) = ff_montgomery1_sum_2_mults(x0,x1,y0,y1))
#define _ff_sum_2_mults_d1(z,x0,x1,y0,y1)	((z) = ff_montgomery1_sum_2_mults_d1(x0,x1,y0,y1))
#define _ff_sum_2_mults_s2(z,x,y)			((z) = ff_montgomery1_dbl_mult(x,y))
#define _ff_dbl_mult(z,x,y)				((z) = ff_montgomery1_sum_2_mults_s2(x,y))
#define _ff_sum_3_mults(z,x0,x1,x2,y0,y1,y2)			((z) = ff_montgomery1_sum_3_mults(x0,x1,x2,y0,y1,y2))
#define _ff_sum_3_mults_d1(z,x0,x1,x2,y0,y1,y2)		((z) = ff_montgomery1_sum_3_mults_d1(x0,x1,x2,y0,y1,y2))
#define _ff_sum_3_mults_d2(z,x0,x1,x2,y0,y1,y2)		((z) = ff_montgomery1_sum_3_mults_d1(x0,x1,x2,y0,y1,y2))
#define _ff_sum_3_mults_s3(z,x0,x1,x2)			((z) = ff_montgomery1_sum_3_mults_s3(x0,x1,x2))
#define _ff_sum_4_mults(z,x0,x1,x2,x3,y0,y1,y2,y3)	((z) = ff_montgomery1_sum_4_mults(x0,x1,x2,x3,y0,y1,y2,y3))
#define _ff_sum_4_mults_d1(z,x0,x1,x2,x3,y0,y1,y2,y3)	((z) = ff_montgomery1_sum_4_mults_d1(x0,x1,x2,x3,y0,y1,y2,y3))
#define _ff_sum_4_mults_d2(z,x0,x1,x2,x3,y0,y1,y2,y3)	((z) = ff_montgomery1_sum_4_mults_d2(x0,x1,x2,x3,y0,y1,y2,y3))
#define _ff_sum_4_mults_s4(z,x0,x1,x2,x3)			((z) = ff_montgomery1_sum_4_mults_s4(x0,x1,x2,x3))
#define _ff_sum_5_mults(z,x0,x1,x2,x3,x4,y0,y1,y2,y3,y4)			((z) = ff_montgomery1_sum_5_mults(x0,x1,x2,x3,x4,y0,y1,y2,y3,y4))
#define _ff_sum_5_mults_d1(z,x0,x1,x2,x3,x4,y0,y1,y2,y3,y4)		((z) = ff_montgomery1_sum_5_mults_d1(x0,x1,x2,x3,x4,y0,y1,y2,y3,y4))
#define _ff_sum_5_mults_s5(z,x0,x1,x2,x3,x4)					((z) = ff_montgomery1_sum_5_mults_s5(x0,x1,x2,x3,x4))
#define _ff_sum_6_mults(z,x0,x1,x2,x3,x4,x5,y0,y1,y2,y3,y4,y5)	((z) = ff_montgomery1_sum_6_mults(x0,x1,x2,x3,x4,x5,y0,y1,y2,y3,y4,y5))
#define _ff_sum_6_mults_d2(z,x0,x1,x2,x3,x4,x5,y0,y1,y2,y3,y4,y5)	((z) = ff_montgomery1_sum_6_mults_d2(x0,x1,x2,x3,x4,x5,y0,y1,y2,y3,y4,y5))
#define _ff_sum_6_mults_d1(z,x0,x1,x2,x3,x4,x5,y0,y1,y2,y3,y4,y5)	((z) = ff_montgomery1_sum_6_mults_d1(x0,x1,x2,x3,x4,x5,y0,y1,y2,y3,y4,y5))
#define _ff_sum_6_mults_s6(z,x0,x1,x2,x3,x4,x5)				((z) = ff_montgomery1_sum_6_mults_s6(x0,x1,x2,x3,x4,x5))
#define _ff_sum_7_mults(z,x0,x1,x2,x3,x4,x5,x6,y0,y1,y2,y3,y4,y5,y6)		((z) = ff_montgomery1_sum_7_mults(x0,x1,x2,x3,x4,x5,x6,y0,y1,y2,y3,y4,y5,y6))
#define _ff_sum_7_mults_d3(z,x0,x1,x2,x3,x4,x5,x6,y0,y1,y2,y3,y4,y5,y6)		((z) = ff_montgomery1_sum_7_mults_d3(x0,x1,x2,x3,x4,x5,x6,y0,y1,y2,y3,y4,y5,y6))
#define _ff_sum_7_mults_d2(z,x0,x1,x2,x3,x4,x5,x6,y0,y1,y2,y3,y4,y5,y6)		((z) = ff_montgomery1_sum_7_mults_d2(x0,x1,x2,x3,x4,x5,x6,y0,y1,y2,y3,y4,y5,y6))
#define _ff_sum_7_mults_s7(z,x0,x1,x2,x3,x4,x5,x6)						((z) = ff_montgomery1_sum_7_mults_s7(x0,x1,x2,x3,x4,x5,x6))
#define _ff_sum_8_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7)	((z) = ff_montgomery1_sum_8_mults(x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7))
#define _ff_sum_8_mults_d2(z,x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7)	((z) = ff_montgomery1_sum_8_mults_d2(x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7))
#define _ff_sum_8_mults_s8(z,x0,x1,x2,x3,x4,x5,x6,x7)					((z) = ff_montgomery1_sum_8_mults_s8(x0,x1,x2,x3,x4,x5,x6,x7))
#define _ff_sum_9_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,y0,y1,y2,y3,y4,y5,y6,y7,y8)	((z) = ff_montgomery1_sum_9_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,y0,y1,y2,y3,y4,y5,y6,y7,y8))
#define _ff_sum_9_mults_d2(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,y0,y1,y2,y3,y4,y5,y6,y7,y8)	((z) = ff_montgomery1_sum_9_mults_d2(x0,x1,x2,x3,x4,x5,x6,x7,x8,y0,y1,y2,y3,y4,y5,y6,y7,y8))
#define _ff_sum_9_mults_d3(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,y0,y1,y2,y3,y4,y5,y6,y7,y8)	((z) = ff_montgomery1_sum_9_mults_d3(x0,x1,x2,x3,x4,x5,x6,x7,x8,y0,y1,y2,y3,y4,y5,y6,y7,y8))
#define _ff_sum_9_mults_s9(z,x0,x1,x2,x3,x4,x5,x6,x7,x8)					((z) = ff_montgomery1_sum_9_mults_s9(x0,x1,x2,x3,x4,x5,x6,x7,x8))
#define _ff_sum_10_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)	((z) = ff_montgomery1_sum_10_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9))
#define _ff_sum_10_mults_d3(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)	((z) = ff_montgomery1_sum_10_mults_d3(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9))
#define _ff_sum_10_mults_d4(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)	((z) = ff_montgomery1_sum_10_mults_d4(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9))
#define _ff_sum_10_mults_s10(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9)					((z) = ff_montgomery1_sum_10_mults_s10(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9))
#define _ff_sum_11_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10)	((z) = ff_montgomery1_sum_11_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10))
#define _ff_sum_11_mults_d3(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10)	((z) = ff_montgomery1_sum_11_mults_d3(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10))
#define _ff_sum_11_mults_s11(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)							((z) = ff_montgomery1_sum_11_mults_s11(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10))
#define _ff_sum_12_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)	((z) = ff_montgomery1_sum_12_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11))
#define _ff_sum_12_mults_d3(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)	((z) = ff_montgomery1_sum_12_mults_d3(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11))
#define _ff_sum_12_mults_d4(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)	((z) = ff_montgomery1_sum_12_mults_d4(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11))
#define _ff_sum_12_mults_s12(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)								((z) = ff_montgomery1_sum_12_mults_s12(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11))
#define _ff_sum_13_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12)	((z) = ff_montgomery1_sum_13_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12))
#define _ff_sum_13_mults_d4(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12)	((z) = ff_montgomery1_sum_13_mults_d4(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12))
#define _ff_sum_13_mults_d5(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12)	((z) = ff_montgomery1_sum_13_mults_d5(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12))
#define _ff_sum_13_mults_s13(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)								((z) = ff_montgomery1_sum_13_mults_s13(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12))
#define _ff_sum_14_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13)	((z) = ff_montgomery1_sum_14_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13))
#define _ff_sum_14_mults_d4(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13)	((z) = ff_montgomery1_sum_14_mults_d4(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13))
#define _ff_sum_14_mults_s14(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13)									((z) = ff_montgomery1_sum_14_mults_s14(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13))
#define _ff_sum_15_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14)	((z) = ff_montgomery1_sum_15_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14))
#define _ff_sum_15_mults_d4(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14)	((z) = ff_montgomery1_sum_15_mults_d4(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14))
#define _ff_sum_15_mults_d5(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14)	((z) = ff_montgomery1_sum_15_mults_d5(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14))
#define _ff_sum_15_mults_s15(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14)									((z) = ff_montgomery1_sum_15_mults_s15(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14))
#define _ff_sum_16_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15)	((z) = ff_montgomery1_sum_16_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15))
#define _ff_sum_16_mults_d5(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15)	((z) = ff_montgomery1_sum_16_mults_d5(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15))
#define _ff_sum_16_mults_d6(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15)	((z) = ff_montgomery1_sum_16_mults_d6(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15))
#define _ff_sum_16_mults_s16(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15)											((z) = ff_montgomery1_sum_16_mults_s16(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15))
#define _ff_sum_17_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16)	((z) = ff_montgomery1_sum_17_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16))
#define _ff_sum_17_mults_d5(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16)	((z) = ff_montgomery1_sum_17_mults_d5(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16))
#define _ff_sum_17_mults_s17(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16)												((z) = ff_montgomery1_sum_17_mults_s17(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16))
#define _ff_sum_18_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17)	((z) = ff_montgomery1_sum_18_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17))
#define _ff_sum_18_mults_d5(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17)	((z) = ff_montgomery1_sum_18_mults_d5(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17))
#define _ff_sum_18_mults_d6(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17)	((z) = ff_montgomery1_sum_18_mults_d6(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17))
#define _ff_sum_18_mults_s18(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17)												((z) = ff_montgomery1_sum_18_mults_s18(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17))
#define _ff_sum_19_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18)	((z) = ff_montgomery1_sum_19_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18))
#define _ff_sum_19_mults_d6(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18)	((z) = ff_montgomery1_sum_19_mults_d6(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18))
#define _ff_sum_19_mults_d7(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18)	((z) = ff_montgomery1_sum_19_mults_d7(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18))
#define _ff_sum_19_mults_s19(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18)													((z) = ff_montgomery1_sum_19_mults_s19(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18))
#define _ff_sum_20_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19)	((z) = ff_montgomery1_sum_20_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19))
#define _ff_sum_20_mults_d6(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19)	((z) = ff_montgomery1_sum_20_mults_d6(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19))
#define _ff_sum_20_mults_s20(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19)														((z) = ff_montgomery1_sum_20_mults_s20(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19))
#define _ff_sum_21_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20)		((z) = ff_montgomery1_sum_21_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20))
#define _ff_sum_21_mults_d6(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20)	((z) = ff_montgomery1_sum_21_mults_d6(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20))
#define _ff_sum_21_mults_d7(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20)	((z) = ff_montgomery1_sum_21_mults_d7(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20))
#define _ff_sum_21_mults_s21(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20)															((z) = ff_montgomery1_sum_21_mults_s21(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20))
#define _ff_sum_22_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21)		((z) = ff_montgomery1_sum_22_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21))
#define _ff_sum_22_mults_d7(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21)	((z) = ff_montgomery1_sum_22_mults_d7(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21))
#define _ff_sum_22_mults_d8(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21)	((z) = ff_montgomery1_sum_22_mults_d8(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21))
#define _ff_sum_22_mults_s22(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21)																	((z) = ff_montgomery1_sum_22_mults_s22(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21))
#define _ff_sum_23_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22)		((z) = ff_montgomery1_sum_23_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22))
#define _ff_sum_23_mults_d7(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22)	((z) = ff_montgomery1_sum_23_mults_d7(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22))
#define _ff_sum_23_mults_s23(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22)																	((z) = ff_montgomery1_sum_23_mults_s23(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22))
#define _ff_sum_24_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23)		((z) = ff_montgomery1_sum_24_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23))
#define _ff_sum_24_mults_d7(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23)	((z) = ff_montgomery1_sum_24_mults_d7(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23))
#define _ff_sum_24_mults_d8(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23)	((z) = ff_montgomery1_sum_24_mults_d8(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23))
#define _ff_sum_24_mults_s24(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23)																		((z) = ff_montgomery1_sum_24_mults_s24(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23))
#define _ff_sum_25_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24)	((z) = ff_montgomery1_sum_25_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24))
#define _ff_sum_25_mults_d8(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24)	((z) = ff_montgomery1_sum_25_mults_d8(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24))
#define _ff_sum_25_mults_d9(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24)	((z) = ff_montgomery1_sum_25_mults_d9(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24))
#define _ff_sum_25_mults_s25(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24)																		((z) = ff_montgomery1_sum_25_mults_s25(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24))
#define _ff_sum_26_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25)	((z) = ff_montgomery1_sum_26_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25))
#define _ff_sum_26_mults_d8(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25)	((z) = ff_montgomery1_sum_26_mults_d8(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25))
#define _ff_sum_26_mults_s26(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25)																			((z) = ff_montgomery1_sum_26_mults_s26(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25))
#define _ff_sum_27_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26)		((z) = ff_montgomery1_sum_27_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26))
#define _ff_sum_27_mults_d8(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26)	((z) = ff_montgomery1_sum_27_mults_d8(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26))
#define _ff_sum_27_mults_d9(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26)	((z) = ff_montgomery1_sum_27_mults_d9(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26))
#define _ff_sum_27_mults_s27(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26)																					((z) = ff_montgomery1_sum_27_mults_s27(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26))
#define _ff_sum_28_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27)		((z) = ff_montgomery1_sum_28_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27))
#define _ff_sum_28_mults_d9(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27)		((z) = ff_montgomery1_sum_28_mults_d9(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27))
#define _ff_sum_28_mults_d10(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27)		((z) = ff_montgomery1_sum_28_mults_d10(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27))
#define _ff_sum_28_mults_s28(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27)																							((z) = ff_montgomery1_sum_28_mults_s28(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27))
#define _ff_sum_29_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28)		((z) = ff_montgomery1_sum_29_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28))
#define _ff_sum_29_mults_d9(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28)		((z) = ff_montgomery1_sum_29_mults_d9(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28))
#define _ff_sum_29_mults_s29(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28)																							((z) = ff_montgomery1_sum_29_mults_s29(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28))
#define _ff_sum_30_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29)		((z) = ff_montgomery1_sum_30_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29))
#define _ff_sum_30_mults_d9(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29)		((z) = ff_montgomery1_sum_30_mults_d9(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29))
#define _ff_sum_30_mults_d10(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29)		((z) = ff_montgomery1_sum_30_mults_d10(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29))
#define _ff_sum_30_mults_s30(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29)																								((z) = ff_montgomery1_sum_30_mults_s30(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29))
#define _ff_sum_31_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30)		((z) = ff_montgomery1_sum_31_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30))
#define _ff_sum_31_mults_d10(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30)		((z) = ff_montgomery1_sum_31_mults_d10(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30))
#define _ff_sum_31_mults_d11(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30)		((z) = ff_montgomery1_sum_31_mults_d11(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30))
#define _ff_sum_31_mults_s31(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30)																								((z) = ff_montgomery1_sum_31_mults_s31(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30))
#define _ff_sum_32_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31)		((z) = ff_montgomery1_sum_32_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31))
#define _ff_sum_32_mults_d10(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31)		((z) = ff_montgomery1_sum_32_mults_d10(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31))
#define _ff_sum_33_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32)		((z) = ff_montgomery1_sum_33_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32))
#define _ff_sum_33_mults_d10(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32)		((z) = ff_montgomery1_sum_33_mults_d10(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32))
#define _ff_sum_33_mults_d11(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32)		((z) = ff_montgomery1_sum_33_mults_d11(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32))
#define _ff_sum_34_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33)    ((z) = ff_montgomery1_sum_34_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33))
#define _ff_sum_34_mults_d11(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33)    ((z) = ff_montgomery1_sum_34_mults_d11(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33))
#define _ff_sum_34_mults_d12(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33)    ((z) = ff_montgomery1_sum_34_mults_d12(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33))
#define _ff_sum_35_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34)    ((z) = ff_montgomery1_sum_35_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34))
#define _ff_sum_35_mults_d11(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34)    ((z) = ff_montgomery1_sum_35_mults_d11(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34))
#define _ff_sum_35_mults_d12(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34)    ((z) = ff_montgomery1_sum_35_mults_d12(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34))
#define _ff_sum_36_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35)    ((z) = ff_montgomery1_sum_36_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35))
#define _ff_sum_36_mults_d11(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35)    ((z) = ff_montgomery1_sum_36_mults_d11(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35))
#define _ff_sum_36_mults_d12(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35)    ((z) = ff_montgomery1_sum_36_mults_d12(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35))
#define _ff_sum_37_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36)    ((z) = ff_montgomery1_sum_37_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36))
#define _ff_sum_37_mults_d12(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36)    ((z) = ff_montgomery1_sum_37_mults_d12(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36))
#define _ff_sum_37_mults_d13(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36)    ((z) = ff_montgomery1_sum_37_mults_d13(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36))
#define _ff_sum_38_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37)    ((z) = ff_montgomery1_sum_38_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37))
#define _ff_sum_38_mults_d12(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37)    ((z) = ff_montgomery1_sum_38_mults_d12(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37))
#define _ff_sum_38_mults_d13(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37)    ((z) = ff_montgomery1_sum_38_mults_d13(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37))
#define _ff_sum_39_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38)    ((z) = ff_montgomery1_sum_39_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38))
#define _ff_sum_39_mults_d12(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38)    ((z) = ff_montgomery1_sum_39_mults_d12(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38))
#define _ff_sum_39_mults_d13(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38)    ((z) = ff_montgomery1_sum_39_mults_d13(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38))
#define _ff_sum_40_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39)    ((z) = ff_montgomery1_sum_40_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39))
#define _ff_sum_40_mults_d13(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39)    ((z) = ff_montgomery1_sum_40_mults_d13(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39))
#define _ff_sum_40_mults_d14(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39)    ((z) = ff_montgomery1_sum_40_mults_d14(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39))
#define _ff_sum_41_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40)    ((z) = ff_montgomery1_sum_41_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40))
#define _ff_sum_41_mults_d13(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40)    ((z) = ff_montgomery1_sum_41_mults_d13(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40))
#define _ff_sum_41_mults_d14(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40)    ((z) = ff_montgomery1_sum_41_mults_d14(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40))
#define _ff_sum_42_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41)    ((z) = ff_montgomery1_sum_42_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41))
#define _ff_sum_42_mults_d13(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41)    ((z) = ff_montgomery1_sum_42_mults_d13(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41))
#define _ff_sum_42_mults_d14(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41)    ((z) = ff_montgomery1_sum_42_mults_d14(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41))
#define _ff_sum_43_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42)    ((z) = ff_montgomery1_sum_43_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42))
#define _ff_sum_43_mults_d14(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42)    ((z) = ff_montgomery1_sum_43_mults_d14(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42))
#define _ff_sum_43_mults_d15(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42)    ((z) = ff_montgomery1_sum_43_mults_d15(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42))
#define _ff_sum_44_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43)    ((z) = ff_montgomery1_sum_44_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43))
#define _ff_sum_44_mults_d14(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43)    ((z) = ff_montgomery1_sum_44_mults_d14(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43))
#define _ff_sum_44_mults_d15(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43)    ((z) = ff_montgomery1_sum_44_mults_d15(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43))
#define _ff_sum_45_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44)    ((z) = ff_montgomery1_sum_45_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44))
#define _ff_sum_45_mults_d14(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44)    ((z) = ff_montgomery1_sum_45_mults_d14(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44))
#define _ff_sum_45_mults_d15(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44)    ((z) = ff_montgomery1_sum_45_mults_d15(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44))
#define _ff_sum_46_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45)    ((z) = ff_montgomery1_sum_46_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45))
#define _ff_sum_46_mults_d15(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45)    ((z) = ff_montgomery1_sum_46_mults_d15(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45))
#define _ff_sum_46_mults_d16(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45)    ((z) = ff_montgomery1_sum_46_mults_d16(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45))
#define _ff_sum_47_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46)    ((z) = ff_montgomery1_sum_47_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46))
#define _ff_sum_47_mults_d15(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46)    ((z) = ff_montgomery1_sum_47_mults_d15(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46))
#define _ff_sum_47_mults_d16(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46)    ((z) = ff_montgomery1_sum_47_mults_d16(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46))
#define _ff_sum_48_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47)    ((z) = ff_montgomery1_sum_48_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47))
#define _ff_sum_48_mults_d15(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47)    ((z) = ff_montgomery1_sum_48_mults_d15(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47))
#define _ff_sum_48_mults_d16(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47)    ((z) = ff_montgomery1_sum_48_mults_d16(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47))
#define _ff_sum_49_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48)    ((z) = ff_montgomery1_sum_49_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48))
#define _ff_sum_49_mults_d16(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48)    ((z) = ff_montgomery1_sum_49_mults_d16(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48))
#define _ff_sum_49_mults_d17(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48)    ((z) = ff_montgomery1_sum_49_mults_d17(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48))
#define _ff_sum_50_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49)    ((z) = ff_montgomery1_sum_50_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49))
#define _ff_sum_50_mults_d16(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49)    ((z) = ff_montgomery1_sum_50_mults_d16(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49))
#define _ff_sum_50_mults_d17(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49)    ((z) = ff_montgomery1_sum_50_mults_d17(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49))
#define _ff_sum_51_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50)    ((z) = ff_montgomery1_sum_51_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50))
#define _ff_sum_51_mults_d16(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50)    ((z) = ff_montgomery1_sum_51_mults_d16(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50))
#define _ff_sum_51_mults_d17(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50)    ((z) = ff_montgomery1_sum_51_mults_d17(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50))
#define _ff_sum_52_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51)    ((z) = ff_montgomery1_sum_52_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51))
#define _ff_sum_52_mults_d17(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51)    ((z) = ff_montgomery1_sum_52_mults_d17(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51))
#define _ff_sum_52_mults_d18(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51)    ((z) = ff_montgomery1_sum_52_mults_d18(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51))
#define _ff_sum_53_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52)    ((z) = ff_montgomery1_sum_53_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52))
#define _ff_sum_53_mults_d17(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52)    ((z) = ff_montgomery1_sum_53_mults_d17(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52))
#define _ff_sum_53_mults_d18(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52)    ((z) = ff_montgomery1_sum_53_mults_d18(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52))
#define _ff_sum_54_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53)    ((z) = ff_montgomery1_sum_54_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53))
#define _ff_sum_54_mults_d17(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53)    ((z) = ff_montgomery1_sum_54_mults_d17(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53))
#define _ff_sum_54_mults_d18(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53)    ((z) = ff_montgomery1_sum_54_mults_d18(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53))
#define _ff_sum_55_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54)    ((z) = ff_montgomery1_sum_55_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54))
#define _ff_sum_55_mults_d18(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54)    ((z) = ff_montgomery1_sum_55_mults_d18(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54))
#define _ff_sum_55_mults_d19(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54)    ((z) = ff_montgomery1_sum_55_mults_d19(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54))
#define _ff_sum_56_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55)    ((z) = ff_montgomery1_sum_56_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55))
#define _ff_sum_56_mults_d18(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55)    ((z) = ff_montgomery1_sum_56_mults_d18(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55))
#define _ff_sum_56_mults_d19(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55)    ((z) = ff_montgomery1_sum_56_mults_d19(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55))
#define _ff_sum_57_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56)    ((z) = ff_montgomery1_sum_57_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56))
#define _ff_sum_57_mults_d18(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56)    ((z) = ff_montgomery1_sum_57_mults_d18(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56))
#define _ff_sum_57_mults_d19(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56)    ((z) = ff_montgomery1_sum_57_mults_d19(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56))
#define _ff_sum_58_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57)    ((z) = ff_montgomery1_sum_58_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57))
#define _ff_sum_58_mults_d19(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57)    ((z) = ff_montgomery1_sum_58_mults_d19(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57))
#define _ff_sum_58_mults_d20(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57)    ((z) = ff_montgomery1_sum_58_mults_d20(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57))
#define _ff_sum_59_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58)    ((z) = ff_montgomery1_sum_59_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58))
#define _ff_sum_59_mults_d19(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58)    ((z) = ff_montgomery1_sum_59_mults_d19(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58))
#define _ff_sum_59_mults_d20(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58)    ((z) = ff_montgomery1_sum_59_mults_d20(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58))
#define _ff_sum_60_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59)    ((z) = ff_montgomery1_sum_60_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59))
#define _ff_sum_60_mults_d19(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59)    ((z) = ff_montgomery1_sum_60_mults_d19(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59))
#define _ff_sum_60_mults_d20(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59)    ((z) = ff_montgomery1_sum_60_mults_d20(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59))
#define _ff_sum_61_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60)    ((z) = ff_montgomery1_sum_61_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60))
#define _ff_sum_61_mults_d20(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60)    ((z) = ff_montgomery1_sum_61_mults_d20(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60))
#define _ff_sum_61_mults_d21(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60)    ((z) = ff_montgomery1_sum_61_mults_d21(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60))
#define _ff_sum_62_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61)    ((z) = ff_montgomery1_sum_62_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61))
#define _ff_sum_62_mults_d20(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61)    ((z) = ff_montgomery1_sum_62_mults_d20(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61))
#define _ff_sum_62_mults_d21(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61)    ((z) = ff_montgomery1_sum_62_mults_d21(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61))
#define _ff_sum_63_mults(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62)    ((z) = ff_montgomery1_sum_63_mults(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62))
#define _ff_sum_63_mults_d20(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62)    ((z) = ff_montgomery1_sum_63_mults_d20(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62))
#define _ff_sum_63_mults_d21(z,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62)    ((z) = ff_montgomery1_sum_63_mults_d21(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62))
#define _ff_sum_1_mults_arr(z,x,y)		_ff_mult(z,(x)[0],(y)[0])
#define _ff_sum_2_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_2_mults((x)[0],(x)[1],(y)[0],(y)[1]))
#define _ff_sum_3_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_3_mults((x)[0],(x)[1],(x)[2],(y)[0],(y)[1],(y)[2]))
#define _ff_sum_4_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_4_mults((x)[0],(x)[1],(x)[2],(x)[3],(y)[0],(y)[1],(y)[2],(y)[3]))
#define _ff_sum_5_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_5_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4]))
#define _ff_sum_6_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_6_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5]))
#define _ff_sum_7_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_7_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6]))
#define _ff_sum_8_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_8_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7]))
#define _ff_sum_9_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_9_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8]))
#define _ff_sum_10_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_10_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9]))
#define _ff_sum_11_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_11_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10]))
#define _ff_sum_12_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_12_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11]))
#define _ff_sum_13_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_13_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12]))
#define _ff_sum_14_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_14_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13]))
#define _ff_sum_15_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_15_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14]))
#define _ff_sum_16_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_16_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15]))
#define _ff_sum_17_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_17_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16]))
#define _ff_sum_18_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_18_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17]))
#define _ff_sum_19_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_19_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18]))
#define _ff_sum_20_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_20_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19]))
#define _ff_sum_21_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_21_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20]))
#define _ff_sum_22_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_22_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20],(y)[21]))
#define _ff_sum_23_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_23_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20],(y)[21],(y)[22]))
#define _ff_sum_24_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_24_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20],(y)[21],(y)[22],(y)[23]))
#define _ff_sum_25_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_25_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20],(y)[21],(y)[22],(y)[23],(y)[24]))
#define _ff_sum_26_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_26_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20],(y)[21],(y)[22],(y)[23],(y)[24],(y)[25]))
#define _ff_sum_27_mults_arr(z,x,y)		((z) = ff_montgomery1_sum_27_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(x)[26],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20],(y)[21],(y)[22],(y)[23],(y)[24],(y)[25],(y)[26]))
#define _ff_sum_28_mults_arr(z,x,y) 	((z) = ff_montgomery1_sum_28_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(x)[26],(x)[27],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20],(y)[21],(y)[22],(y)[23],(y)[24],(y)[25],(y)[26],(y)[27]))
#define _ff_sum_29_mults_arr(z,x,y) 	((z) = ff_montgomery1_sum_29_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(x)[26],(x)[27],(x)[28],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20],(y)[21],(y)[22],(y)[23],(y)[24],(y)[25],(y)[26],(y)[27],(y)[28]))
#define _ff_sum_30_mults_arr(z,x,y) 	((z) = ff_montgomery1_sum_30_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(x)[26],(x)[27],(x)[28],(x)[29],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20],(y)[21],(y)[22],(y)[23],(y)[24],(y)[25],(y)[26],(y)[27],(y)[28],(y)[29]))
#define _ff_sum_31_mults_arr(z,x,y) 	((z) = ff_montgomery1_sum_31_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(x)[26],(x)[27],(x)[28],(x)[29],(x)[30],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20],(y)[21],(y)[22],(y)[23],(y)[24],(y)[25],(y)[26],(y)[27],(y)[28],(y)[29],(y)[30]))
#define _ff_sum_32_mults_arr(z,x,y) 	((z) = ff_montgomery1_sum_32_mults((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(x)[26],(x)[27],(x)[28],(x)[29],(x)[30],(x)[31],(y)[0],(y)[1],(y)[2],(y)[3],(y)[4],(y)[5],(y)[6],(y)[7],(y)[8],(y)[9],(y)[10],(y)[11],(y)[12],(y)[13],(y)[14],(y)[15],(y)[16],(y)[17],(y)[18],(y)[19],(y)[20],(y)[21],(y)[22],(y)[23],(y)[24],(y)[25],(y)[26],(y)[27],(y)[28],(y)[29],(y)[30],(y)[31]))
#define _ff_sum_1_mults_s1_arr(z,x)	_ff_square(z,(x)[0])
#define _ff_sum_2_mults_s2_arr(z,x)		((z) = ff_montgomery1_sum_2_mults_s2((x)[0],(x)[1]))
#define _ff_sum_3_mults_s3_arr(z,x)		((z) = ff_montgomery1_sum_3_mults_s3((x)[0],(x)[1],(x)[2]))
#define _ff_sum_4_mults_s4_arr(z,x)		((z) = ff_montgomery1_sum_4_mults_s4((x)[0],(x)[1],(x)[2],(x)[3]))
#define _ff_sum_5_mults_s5_arr(z,x)		((z) = ff_montgomery1_sum_5_mults_s5((x)[0],(x)[1],(x)[2],(x)[3],(x)[4]))
#define _ff_sum_6_mults_s6_arr(z,x)		((z) = ff_montgomery1_sum_6_mults_s6((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5]))
#define _ff_sum_7_mults_s7_arr(z,x)		((z) = ff_montgomery1_sum_7_mults_s7((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6]))
#define _ff_sum_8_mults_s8_arr(z,x)		((z) = ff_montgomery1_sum_8_mults_s8((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7]))
#define _ff_sum_9_mults_s9_arr(z,x)		((z) = ff_montgomery1_sum_9_mults_s9((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8]))
#define _ff_sum_10_mults_s10_arr(z,x)		((z) = ff_montgomery1_sum_10_mults_s10((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9]))
#define _ff_sum_11_mults_s11_arr(z,x)		((z) = ff_montgomery1_sum_11_mults_s11((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10]))
#define _ff_sum_12_mults_s12_arr(z,x)		((z) = ff_montgomery1_sum_12_mults_s12((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11]))
#define _ff_sum_13_mults_s13_arr(z,x)		((z) = ff_montgomery1_sum_13_mults_s13((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12]))
#define _ff_sum_14_mults_s14_arr(z,x)		((z) = ff_montgomery1_sum_14_mults_s14((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13]))
#define _ff_sum_15_mults_s15_arr(z,x)		((z) = ff_montgomery1_sum_15_mults_s15((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14]))
#define _ff_sum_16_mults_s16_arr(z,x)		((z) = ff_montgomery1_sum_16_mults_s16((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15]))
#define _ff_sum_17_mults_s17_arr(z,x)		((z) = ff_montgomery1_sum_17_mults_s17((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16]))
#define _ff_sum_18_mults_s18_arr(z,x)		((z) = ff_montgomery1_sum_18_mults_s18((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17]))
#define _ff_sum_19_mults_s19_arr(z,x)		((z) = ff_montgomery1_sum_19_mults_s19((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18]))
#define _ff_sum_20_mults_s20_arr(z,x)		((z) = ff_montgomery1_sum_20_mults_s20((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19]))
#define _ff_sum_21_mults_s21_arr(z,x)		((z) = ff_montgomery1_sum_21_mults_s21((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20]))
#define _ff_sum_22_mults_s22_arr(z,x)		((z) = ff_montgomery1_sum_22_mults_s22((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21]))
#define _ff_sum_23_mults_s23_arr(z,x)		((z) = ff_montgomery1_sum_23_mults_s23((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22]))
#define _ff_sum_24_mults_s24_arr(z,x)		((z) = ff_montgomery1_sum_24_mults_s24((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23]))
#define _ff_sum_25_mults_s25_arr(z,x)		((z) = ff_montgomery1_sum_25_mults_s25((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24]))
#define _ff_sum_26_mults_s26_arr(z,x)		((z) = ff_montgomery1_sum_26_mults_s26((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25]))
#define _ff_sum_27_mults_s27_arr(z,x)		((z) = ff_montgomery1_sum_27_mults_s27((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(x)[26]))
#define _ff_sum_28_mults_s28_arr(z,x)	 	((z) = ff_montgomery1_sum_28_mults_s28((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(x)[26],(x)[27]))
#define _ff_sum_29_mults_s29_arr(z,x) 		((z) = ff_montgomery1_sum_29_mults_s29((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(x)[26],(x)[27],(x)[28]))
#define _ff_sum_30_mults_s30_arr(z,x) 		((z) = ff_montgomery1_sum_30_mults_s30((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(x)[26],(x)[27],(x)[28],(x)[29]))
#define _ff_sum_31_mults_s31_arr(z,x) 		((z) = ff_montgomery1_sum_31_mults_s31((x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7],(x)[8],(x)[9],(x)[10],(x)[11],(x)[12],(x)[13],(x)[14],(x)[15],(x)[16],(x)[17],(x)[18],(x)[19],(x)[20],(x)[21],(x)[22],(x)[23],(x)[24],(x)[25],(x)[26],(x)[27],(x)[28],(x)[29],(x)[30]))

#define _ff_sum_mults(z,x,y,n)  ((z) = ff_montgomery1_sum_mults (x,y,n))
#define _ff_sq_coeff(z,x,k)	((z) = ff_montgomery1_sq_coeff(x,k))

#define _ff_triple_square(z,x)			((z)=ff_montgomery1_triple_square(x))
#define _ff_square_double_sub(z,x,y)	((z)=ff_montgomery1_square_double_sub(x,y))
#define _ff_diff_mult_sub(z,v,w,x,y)		((z)=ff_montgomery1_diff_mult_sub(v,w,x,y))
#define _ff_diff_mult(z,v,w,x)			((z)=ff_montgomery1_diff_mult(v,w,x))
#define _ff_square_sub_sub(z,w,x,y)	((z)=ff_montgomery1_square_sub_sub(w,x,y))

#define _ff_square(z,x)				_ff_mult(z,x,x)
#define _ff_invert(z,x)				((z)=ff_montgomery1_invert(x))
#else
#define _ff_mult(z,x,y)				((z) = ((x)*(y)) % _ff_p)
#define _ff_square(z,x)				_ff_mult(z,x,x)
#define _ff_invert(z,x)				((z) = ff_ui_inverse(x,_ff_p))
#endif
#define _ff_div2(z,x)					_ff_mult(z,x,_ff_half)

#define _ff_incmult(z,x,w)			{ _ff_set(z,x);_ff_inc(z);_ff_mult(z,z,w); }	// could be optimized
#define _ff_multadd(z,x,y,a)			{_ff_mult(z,x,y); _ff_addto(z,a); }		// no optimization here
#define ff_multadd(z,x,y,a)			_ff_multadd(z,x,y,a)					// overlap of x and z ok but a can't overlap

#define _ff_triple(z,x)				{ _ff_add(z,x,x); _ff_addto(z,x); }		// z and x cannot overlap

// These macros are used for intermediate operations where reduction mod p may not necessary
// With single precision, its generally better to keep things reduced, so we don't attempt to optmize these
#define _ff_qneg(z,x)				_ff_neg(z,x);
#define _ff_qnegsum(z,x,y)			{ _ff_neg(z,x); _ff_subfrom(z,y); }
#define _ff_qsub(z,x,y)				_ff_sub(z,x,y);
#define _ff_qsubfrom(z,x)				_ff_subfrom(z,x);
#define _ff_qadd(z,x,y)				_ff_add(z,x,y);
#define _ff_qaddto(z,x)				_ff_addto(z,x);


unsigned long ff_montgomery1_invert (unsigned long x);
unsigned long ff_ui_inverse (unsigned long a, unsigned long m);
void ff_exp_ui (ff_t o[1], ff_t a[1], unsigned long e);
int ff_precompute_exp_chain (int chain[FF_MAX_CHAIN_LEN], unsigned long e);
void ff_exp_chain (ff_t o[1], ff_t a[1], int chain[], int len);
int ff_invsqrt (ff_t o[1], ff_t a[1], int ext);							// computes 1/sqrt(a), returns 1 if sqrt is rational, 0 if sqrt is in F_p^2. 
void ff_setup_2g (void);	
int ff_cbrt_invcbrt (ff_t o[1], ff_t *oi, ff_t a[1]);							// computes cbrt(a) and optionally invcbrt(a).  returns 1 for success, 0 if no cbrt in F_p
int ff_3Sylow_invcbrt (ff_t o[1], ff_t a[1]);								// computes a^{-1/3} for a in the 3-Sylow subgroup, returns 0 if a is not a residue in the 3-Sylow
static inline int ff_cbrt (ff_t o[1], ff_t a[1]) { return ff_cbrt_invcbrt(o,0,a); }

ff_t *ff_binomials (int n);											// returns pointer to array of binom(n,0), binom(n,1), ..., binom(n,n)
#define _ff_binomial(n,k)				_ff_binomial_tab[((n*n+n)>>1)+k]	// n must be <= _ff_bonimial_top, call ff_binomials(n) to set this

// called automatically by ff_sqrt, but also used in ffext.c - note that for p=3mod4, e=1 and -1 generates the 2-Sylow subgroup
void ff_setup_2g(void);

// ditto for cbrt, only relevant when p=1mod3
void _ff_setup_3g(void);
static inline void ff_setup_3g(void) { if ( ! _ff_cbrt_setup ) _ff_setup_3g(); }
void _ff_setup_fractions(void);
static inline void ff_setup_fractions () { if ( ! _ff_frac[1] ) _ff_setup_fractions(); }

void ff_parallel_invert (ff_t z[], ff_t x[], unsigned n);					// z and x may overlap
int ff_parallel_invert_check (ff_t z[], ff_t x[], unsigned n);				// returns zero if any x is zero, rather than aborting

void ff_complementary_products (ff_t o[], ff_t a[], int n);				// o and a cannot overlap, o must have space for n+ceil(log_2(n)) entries
int ff_cornacchia (long *x, long *y, long d);							// solves p=x^2+dy^2

#define ff_mult(z,x,y)				_ff_mult(z,x,y)					// overlap is ok for these macros (but not side effects)
#define ff_square(z,x)				_ff_square(z,x)
#define ff_invert(z,x)				{ if ( ! _ff_one(x) ) _ff_invert(z,x); else _ff_set(z,x); }
void ff_invert_small_int (ff_t z[], int x);

// end higher arithmetic operations

// random elements generation
extern gmp_randstate_t ff_randstate;
extern int ff_randstate_init;
static inline void ff_randinit (void) { if ( ! ff_randstate_init ) { gmp_randinit_default (ff_randstate);  gmp_randseed_ui (ff_randstate, cstd_seed()); ff_randstate_init = 1; } }
static inline unsigned long ff_randomm_ui (unsigned long m) { ff_randinit(); return gmp_urandomm_ui (ff_randstate, m); }
static inline unsigned long ff_randomb_ui (int b)  { ff_randinit(); return gmp_urandomb_ui (ff_randstate, b); }
#define _ff_random(x)				((x) = ff_randomm_ui(_ff_p))
static inline void ff_random (ff_t *x) { *x = ff_randomm_ui (_ff_p); }		// note, with Montgomery representation, x is *not* set to the integer returned by ff_random_ui, but it is still a random element of F_p
#define _ff_random_nz(x)				((x) = 1+ff_randomm_ui(_ff_p-1))
static inline void ff_random_nz (ff_t *x) { *x = 1+ff_randomm_ui (_ff_p-1); }


void ff_organize (ff_t a[], int n);				// organize sorts the array into an order convenient for searching.  This is *not* the ordering 1 < 1+1 < 1+1+1 < ...
int ff_find (ff_t x, ff_t a[], int n);				// returns the index of the element x in the organized array a[] or -1 if not found.

// Inline Functions
static inline int ff_is_negation (ff_t x, ff_t y)
{
	register ff_t t;
	
	_ff_set(t,x);
	_ff_core_addto(t,y);
	return ( _ff_equal(t,_ff_p) || _ff_zero(t) );
}

#if FF_MONTGOMERY && FF_WORDS == 1 && FF_HALF_WORDS == 1
static  inline unsigned  ff_montgomery1_mult (unsigned  x, unsigned  y)
{
	register unsigned long z;
	register unsigned a;

	z = (unsigned long)x*y;
	a = z*_ff_mont_pni;
	z  += ((unsigned long)a * _ff_p);
	z >>= 32;
	if ( z >= _ff_p ) z -= _ff_p;
	return z;
}
#endif



#if FF_MONTGOMERY && FF_WORDS == 1 && FF_HALF_WORDS == 2

#include "ffmontgomery64.h"

// this is pretty close to optimal  (on an AMD athlon/phenom).  Blocking at 9 or 10 might be slightly faster, but 16 is definitely slower.
static inline void ff_dot_product (ff_t z[1], ff_t x[], ff_t y[], int n)
{
	register int i;
	register ff_t t0,t1;

#if ! FF_FAST
	if ( _ff_p & 0xF800000000000000UL ) { printf ("p=%ld too large in ff_dot_product\n", _ff_p); exit (0); }
#endif

	if ( n < 8 ) {
		switch (n) {
		case 1: _ff_mult(z[0],x[0],y[0]); break;
		case 2: _ff_sum_2_mults(z[0],x[0],x[1],y[1],y[0]); break;
		case 3: _ff_sum_3_mults(z[0],x[0],x[1],x[2],y[2],y[1],y[0]); break;
		case 4: _ff_sum_4_mults(z[0],x[0],x[1],x[2],x[3],y[3],y[2],y[1],y[0]); break;
		case 5: _ff_sum_5_mults(z[0],x[0],x[1],x[2],x[3],x[4],y[4],y[3],y[2],y[1],y[0]); break;
		case 6: _ff_sum_6_mults(z[0],x[0],x[1],x[2],x[3],x[4],x[5],y[5],y[4],y[3],y[2],y[1],y[0]); break;
		case 7: _ff_sum_7_mults(z[0],x[0],x[1],x[2],x[3],x[4],x[5],x[6],y[6],y[5],y[4],y[3],y[2],y[1],y[0]); break;
		default: _ff_set_zero(z[0]);
		}
		return;
	}

	_ff_set_zero(t0);
	for ( i = 0 ; i <= n-16 ; i+=8 ) {
		_ff_sum_8_mults(t1,x[i],x[i+1],x[i+2],x[i+3],x[i+4],x[i+5],x[i+6],x[i+7],y[i+7],y[i+6],y[i+5],y[i+4],y[i+3],y[i+2],y[i+1],y[i]);
		_ff_addto(t0,t1);
	}
	switch ( n-i ) {
	case 8: _ff_sum_8_mults(t1,x[i],x[i+1],x[i+2],x[i+3],x[i+4],x[i+5],x[i+6],x[i+7],y[i+7],y[i+6],y[i+5],y[i+4],y[i+3],y[i+2],y[i+1],y[i]); break;
	case 9: _ff_sum_9_mults(t1,x[i],x[i+1],x[i+2],x[i+3],x[i+4],x[i+5],x[i+6],x[i+7],x[i+8],y[i+8],y[i+7],y[i+6],y[i+5],y[i+4],y[i+3],y[i+2],y[i+1],y[i]); break;
	case 10: _ff_sum_10_mults(t1,x[i],x[i+1],x[i+2],x[i+3],x[i+4],x[i+5],x[i+6],x[i+7],x[i+8],x[i+9],y[i+9],y[i+8],y[i+7],y[i+6],y[i+5],y[i+4],y[i+3],y[i+2],y[i+1],y[i]); break;
	case 11: _ff_sum_11_mults(t1,x[i],x[i+1],x[i+2],x[i+3],x[i+4],x[i+5],x[i+6],x[i+7],x[i+8],x[i+9],x[i+10],y[i+10],y[i+9],y[i+8],y[i+7],y[i+6],y[i+5],y[i+4],y[i+3],y[i+2],y[i+1],y[i]); break;
	case 12: _ff_sum_12_mults(t1,x[i],x[i+1],x[i+2],x[i+3],x[i+4],x[i+5],x[i+6],x[i+7],x[i+8],x[i+9],x[i+10],x[i+11],y[i+11],y[i+10],y[i+9],y[i+8],y[i+7],y[i+6],y[i+5],y[i+4],y[i+3],y[i+2],y[i+1],y[i]); break;
	case 13: _ff_sum_13_mults(t1,x[i],x[i+1],x[i+2],x[i+3],x[i+4],x[i+5],x[i+6],x[i+7],x[i+8],x[i+9],x[i+10],x[i+11],x[i+12],y[i+12],y[i+11],y[i+10],y[i+9],y[i+8],y[i+7],y[i+6],y[i+5],y[i+4],y[i+3],y[i+2],y[i+1],y[i]); break;
	case 14: _ff_sum_14_mults(t1,x[i],x[i+1],x[i+2],x[i+3],x[i+4],x[i+5],x[i+6],x[i+7],x[i+8],x[i+9],x[i+10],x[i+11],x[i+12],x[i+13],y[i+13],y[i+12],y[i+11],y[i+10],y[i+9],y[i+8],y[i+7],y[i+6],y[i+5],y[i+4],y[i+3],y[i+2],y[i+1],y[i]); break;
	case 15: _ff_sum_15_mults(t1,x[i],x[i+1],x[i+2],x[i+3],x[i+4],x[i+5],x[i+6],x[i+7],x[i+8],x[i+9],x[i+10],x[i+11],x[i+12],x[i+13],x[i+14],y[i+14],y[i+13],y[i+12],y[i+11],y[i+10],y[i+9],y[i+8],y[i+7],y[i+6],y[i+5],y[i+4],y[i+3],y[i+2],y[i+1],y[i]); break;
	default: _ff_set_zero(t1);	// unreachable
	}
	_ff_add(z[0],t0,t1);
}


static inline void ff_conv (ff_t z[1], ff_t x[], ff_t y[], int n)
{
	register int i;
	register ff_t t0,t1;

#if ! FF_FAST
	if ( _ff_p & 0xF800000000000000UL ) { printf ("p=%ld too large in ff_dot_product\n", _ff_p); exit (0); }
#endif

	if ( n < 8 ) {
		switch (n) {
		case 0: _ff_set_zero(z[0]); break;
		case 1: _ff_mult(z[0],x[0],y[0]); break;
		case 2: _ff_sum_2_mults_arr(z[0],x,y); break;
		case 3: _ff_sum_3_mults_arr(z[0],x,y); break;
		case 4: _ff_sum_4_mults_arr(z[0],x,y); break;
		case 5: _ff_sum_5_mults_arr(z[0],x,y); break;
		case 6: _ff_sum_6_mults_arr(z[0],x,y); break;
		case 7: _ff_sum_7_mults_arr(z[0],x,y); break;
		}
		return;
	}

	_ff_set_zero(t0);
	_ff_set_zero(t1);		// avoid compiler warning
	for ( i = 0 ; i <= n-16 ; i+=8 ) {
		_ff_sum_8_mults_arr(t1,x+i,y+n-8-i);
		_ff_addto(t0,t1);
	}
	switch ( n-i ) {
	case 8: _ff_sum_8_mults_arr(t1,x+i,y); break;
	case 9: _ff_sum_9_mults_arr(t1,x+i,y); break;
	case 10: _ff_sum_10_mults_arr(t1,x+i,y); break;
	case 11: _ff_sum_11_mults_arr(t1,x+i,y); break;
	case 12: _ff_sum_12_mults_arr(t1,x+i,y); break;
	case 13: _ff_sum_13_mults_arr(t1,x+i,y); break;
	case 14: _ff_sum_14_mults_arr(t1,x+i,y); break;
	case 15: _ff_sum_15_mults_arr(t1,x+i,y); break;
	}
	_ff_add(z[0],t0,t1);
}

// assumes n > 8
static inline void ff_conv_big (ff_t z[1], ff_t x[], ff_t y[], int n)
{
	register int i;
	register ff_t t0,t1;

	_ff_set_zero(t0);
	for ( i = 0 ; i <= n-16 ; i+=8 ) {
		_ff_sum_8_mults_arr(t1,x+i,y+n-8-i);
		_ff_addto(t0,t1);
	}
	switch ( n-i ) {
	case 8: _ff_sum_8_mults_arr(t1,x+i,y); break;
	case 9: _ff_sum_9_mults_arr(t1,x+i,y); break;
	case 10: _ff_sum_10_mults_arr(t1,x+i,y); break;
	case 11: _ff_sum_11_mults_arr(t1,x+i,y); break;
	case 12: _ff_sum_12_mults_arr(t1,x+i,y); break;
	case 13: _ff_sum_13_mults_arr(t1,x+i,y); break;
	case 14: _ff_sum_14_mults_arr(t1,x+i,y); break;
	case 15: _ff_sum_15_mults_arr(t1,x+i,y); break;
	default: _ff_set_zero(t1);
	}
	_ff_add(z[0],t0,t1);
}


// computes x[0]*x[n-1]+x[1]*x[n-2]+...+x[n-2]*x[1]+x[n-1]*x[0]
// this is significantly slower than using the hardwired versions for n < 16
static inline void ff_unary_conv (ff_t z[1], ff_t x[], int n)
{
	register int i;
	register ff_t t0,t1;

#if ! FF_FAST
	if ( _ff_p & 0xF800000000000000UL ) { printf ("p=%ld too large in ff_dot_product\n", _ff_p); exit (0); }
#endif

	if ( n < 16 ) {
		switch (n) {
		case 0: _ff_set_zero(z[0]); break;
		case 1: _ff_square(z[0],x[0]); break;
		case 2: _ff_sum_2_mults_s2_arr(z[0],x); break;
		case 3: _ff_sum_3_mults_s3_arr(z[0],x); break;
		case 4: _ff_sum_4_mults_s4_arr(z[0],x); break;
		case 5: _ff_sum_5_mults_s5_arr(z[0],x); break;
		case 6: _ff_sum_6_mults_s6_arr(z[0],x); break;
		case 7: _ff_sum_7_mults_s7_arr(z[0],x); break;
		case 8: _ff_sum_8_mults_s8_arr(z[0],x); break;
		case 9: _ff_sum_9_mults_s9_arr(z[0],x); break;
		case 10: _ff_sum_10_mults_s10_arr(z[0],x); break;
		case 11: _ff_sum_11_mults_s11_arr(z[0],x); break;
		case 12: _ff_sum_12_mults_s12_arr(z[0],x); break;
		case 13: _ff_sum_13_mults_s13_arr(z[0],x); break;
		case 14: _ff_sum_14_mults_s14_arr(z[0],x); break;
		case 15: _ff_sum_15_mults_s15_arr(z[0],x); break;
		}
		return;
	}

	_ff_set_zero(t0); _ff_set_zero(t1);
	for ( i = 0 ; i <= n/2-16 ; i+= 8 ) {
		_ff_sum_8_mults(t1,x[i],x[i+1],x[i+2],x[i+3],x[i+4],x[i+5],x[i+6],x[i+7],x[n-i-8],x[n-i-7],x[n-i-6],x[n-i-5],x[n-i-4],x[n-i-3],x[n-i-2],x[n-i-1]);
		_ff_addto(t0,t1);  _ff_addto(t0,t1);
	}

	switch (n-2*i) {
	case 16: _ff_sum_16_mults_s16_arr(t1,x+i);  break;
	case 17: _ff_sum_17_mults_s17_arr(t1,x+i);  break;
	case 18: _ff_sum_18_mults_s18_arr(t1,x+i);  break;
	case 19: _ff_sum_19_mults_s19_arr(t1,x+i);  break;
	case 20: _ff_sum_20_mults_s20_arr(t1,x+i);  break;
	case 21: _ff_sum_21_mults_s21_arr(t1,x+i);  break;
	case 22: _ff_sum_22_mults_s22_arr(t1,x+i);  break;
	case 23: _ff_sum_23_mults_s23_arr(t1,x+i);  break;
	case 24: _ff_sum_24_mults_s24_arr(t1,x+i);  break;
	case 25: _ff_sum_25_mults_s25_arr(t1,x+i);  break;
	case 26: _ff_sum_26_mults_s26_arr(t1,x+i);  break;
	case 27: _ff_sum_27_mults_s27_arr(t1,x+i);  break;
	case 28: _ff_sum_28_mults_s28_arr(t1,x+i);  break;
	case 29: _ff_sum_29_mults_s29_arr(t1,x+i);  break;
	case 30: _ff_sum_30_mults_s30_arr(t1,x+i);  break;
	case 31: _ff_sum_31_mults_s31_arr(t1,x+i);  break;
	}
	_ff_add(z[0],t0,t1);
}


// assumes n > 16
static inline void ff_unary_conv_big (ff_t z[1], ff_t x[], int n)
{
	register int i;
	register ff_t t0,t1;

	_ff_set_zero(t0);
	for ( i = 0 ; i <= n/2-16 ; i+= 8 ) {
		_ff_sum_8_mults(t1,x[i],x[i+1],x[i+2],x[i+3],x[i+4],x[i+5],x[i+6],x[i+7],x[n-i-8],x[n-i-7],x[n-i-6],x[n-i-5],x[n-i-4],x[n-i-3],x[n-i-2],x[n-i-1]);
		_ff_addto(t0,t1);  _ff_addto(t0,t1);
	}

	switch (n-2*i) {
	case 16: _ff_sum_16_mults_s16_arr(t1,x+i);  break;
	case 17: _ff_sum_17_mults_s17_arr(t1,x+i);  break;
	case 18: _ff_sum_18_mults_s18_arr(t1,x+i);  break;
	case 19: _ff_sum_19_mults_s19_arr(t1,x+i);  break;
	case 20: _ff_sum_20_mults_s20_arr(t1,x+i);  break;
	case 21: _ff_sum_21_mults_s21_arr(t1,x+i);  break;
	case 22: _ff_sum_22_mults_s22_arr(t1,x+i);  break;
	case 23: _ff_sum_23_mults_s23_arr(t1,x+i);  break;
	case 24: _ff_sum_24_mults_s24_arr(t1,x+i);  break;
	case 25: _ff_sum_25_mults_s25_arr(t1,x+i);  break;
	case 26: _ff_sum_26_mults_s26_arr(t1,x+i);  break;
	case 27: _ff_sum_27_mults_s27_arr(t1,x+i);  break;
	case 28: _ff_sum_28_mults_s28_arr(t1,x+i);  break;
	case 29: _ff_sum_29_mults_s29_arr(t1,x+i);  break;
	case 30: _ff_sum_30_mults_s30_arr(t1,x+i);  break;
	case 31: _ff_sum_31_mults_s31_arr(t1,x+i);  break;
	default: _ff_set_zero(t1);
	}
	_ff_add(z[0],t0,t1);
}

#endif

int _qsort_ff_ui_cmp (const void *a, const void *b);
static inline void ff_sort (ff_t a[], int n)			// sorts according to the ordering 1 < 1+1 < 1+1+1 < ... (this is not very efficient, use ff_organize if possible)
{
	if ( n <= 1 ) return;
	if ( n == 2 ) { if ( _ff_get_ui(a[0]) > _ff_get_ui (a[1]) ) _ff_swap(a[0], a[1]); return; }
	qsort (a, n, sizeof(a[0]), _qsort_ff_ui_cmp);
}

static inline int ff_residue (ff_t z) { ff_t t;  if ( _ff_zero(z) ) return 1;  ff_exp_ui(&t,&z,(_ff_p-1)/2);  return _ff_one(t); } // this is faster than using Legendre unless z is very small
static inline void ff_nonresidue (ff_t o[1]) { if ( ! _ff_2g ) ff_setup_2g(); _ff_set(o[0],_ff_2g); }
static inline void ff_nonresidue_inverse (ff_t o[1]) { if ( ! _ff_2g ) ff_setup_2g(); _ff_set(o[0],_ff_2gi); }
static inline int ff_sqrt(ff_t o[1], ff_t a[1]) { ff_t t;  if ( ! ff_invsqrt(&t,a,0) ) return 0;  _ff_mult(o[0],t,a[0]);  return 1; }
int ff_norm_equation (long *t, long *v, long D);	// solves norm equation 4p=t^2-v^2D for current p (assumes p splits completely in the ring class field for D)

int _ff_fast_fourth_root (ff_t o[1], ff_t x[1]);		// requires p = 3 mod 4	
int _ff_fast_sixth_root (ff_t o[1], ff_t x[1]);		// requires p = 11 mod 12
int _ff_fast_eighth_root (ff_t o[1], ff_t x[1]);		// requires p = 3 mod 4	
int _ff_fast_twelfth_root (ff_t o[1], ff_t x[1]);		// requires p = 3 mod 4	
int _ff_fast_sixteenth_root (ff_t o[1], ff_t x[1]);	// requires p = 3 mod 4	

static inline int ff_fourth_root (ff_t o[1], ff_t x[1]) { if ( (_ff_p&3) == 1 ) { if ( ! ff_sqrt(o,x) ) return 0; return ff_sqrt(o,o); } return _ff_fast_fourth_root (o,x); }
static inline int ff_sixth_root (ff_t o[1], ff_t x[1]) { if ( (_ff_p&3) == 1 || _ff_p1mod3 ) { if ( ! ff_sqrt(o,x) ) return 0; return ff_cbrt(o,o); } else return _ff_fast_sixth_root (o,x); }
static inline int ff_eighth_root (ff_t o[1], ff_t x[1]) { if ( (_ff_p&3) == 1 ) { if ( ! ff_sqrt(o,x) ) return 0; if ( ! ff_sqrt(o,o) ) return 0; return ff_sqrt(o,o); } return _ff_fast_eighth_root (o,x); }
static inline int ff_twelfth_root (ff_t o[1], ff_t x[1]) { if ( (_ff_p&3) == 1 || _ff_p1mod3 ) { if ( ! ff_fourth_root(o,x) ) return 0; return ff_cbrt(o,o); } else return _ff_fast_twelfth_root (o,x); }
static inline int ff_sixteenth_root (ff_t o[1], ff_t x[1]) { if ( (_ff_p&3) == 1 ) { if ( ! ff_sqrt(o,x) ) return 0; if ( ! ff_sqrt(o,o) ) return 0; if ( ! ff_sqrt(o,o) ) return 0;  return ff_sqrt(o,o); } return _ff_fast_sixteenth_root (o,x); }

	
static inline int ff_sqrt_ext(ff_t o[1], ff_t a[1]) { register int sts; ff_t t;  sts=ff_invsqrt(&t,a,1);  _ff_mult(o[0],t,a[0]);  return sts; }
int _ff_2Sylow_invsqrt (ff_t o[1], ff_t a[1], int ext);
static inline int ff_2Sylow_invsqrt (ff_t o[1], ff_t a[1], int ext)			// computes a^{-1/2} for a in the 2-Sylow subgroup, returns 0 if a^{-1/2} is not in F_p and returns a^{-1/2}/z
{
	// handle easy cases inline, otherwise call _ff_2Sylow_invsqrt (if e <= 2 this will never happen)
	ff_setup_2g();
	if ( _ff_one(a[0]) ) { _ff_set_one(o[0]);  return 1; }				// use 1 as the inverse square root of 1
	if ( _ff_equal(a[0],_ff_2g) ) { _ff_set(o[0],_ff_2gi);  return 0; }		// if e=1, this must happen
	if ( _ff_equal(a[0],_ff_negone) ) { _ff_set(o[0], _ff_2Sylow_tab[_ff_p2_e-2][1]); return 1;}
	if ( _ff_equal(a[0],_ff_2gi) ) { _ff_set_one (o[0]);  return 0; }			// if e=2, this must happen
	return _ff_2Sylow_invsqrt (o, a,ext);
}

static inline void ff_least_nonresidue (ff_t z[1])
{
	unsigned long n;

	for ( n = 2 ; ui_legendre (n, _ff_p) >= 0 ; n++ );
	_ff_set_ui (z[0],n);
}
static inline void ff_least_sextic_nonresidue (ff_t z[1])
{
	unsigned long n;
	ff_t t;

	assert (_ff_p1mod3);
	for ( n = 2 ; ; n++ ) {
		if ( ui_legendre (n, _ff_p) >= 0 ) continue;
		_ff_set_ui (z[0], n );
		ff_exp_ui (&t,z,(_ff_p-1)/3);
		if ( ! _ff_one(t) ) return;
	}
}


static inline void ff_exp_small (ff_t o[1], ff_t a[1], int e)
{
	register ff_t t0,t1;
	
	switch (e) {
	case 0: _ff_set_one(o[0]); break;
	case 1: _ff_set(o[0],a[0]); break;
	case 2: _ff_square(o[0],a[0]); break;
	case 3: _ff_square(t0,a[0]); _ff_mult(o[0],t0,a[0]); break;
	case 4: _ff_square(t0,a[0]); _ff_square(o[0],t0); break;
	case 5: _ff_square(t0,a[0]); _ff_square(t1,t0); _ff_mult(o[0],t1,a[0]); break;
	case 6: _ff_square(t0,a[0]); _ff_square(t1,t0); _ff_mult(o[0],t0,t1); break;
	case 8: _ff_square(t0,a[0]); _ff_square(t1,t0); _ff_square(o[0],t1); break;
	default: ff_exp_ui(o,a,e);
	}
}

#ifdef __cplusplus
}
#endif

#endif
