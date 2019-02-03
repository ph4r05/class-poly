#ifndef _ECURVE_INCLUDE_
#define _ECURVE_INCLUDE_

/*
     Copyright (c) 2011-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <gmp.h>
#include "ff_poly.h"
//#include "mpzutil.h"

#ifdef __cplusplus
extern "C" {
#endif

// some counters used for performance diagnostics
extern unsigned long ecurve_steps;
extern unsigned long ecurve_retries;

#define ECURVE_ORDER_RETRIES		20
#define ECURVE_SHORT_INTERVAL		5			// don't bother with full BSGS search for very short intervals
#define ECURVE_4TOR_MINP			(1<<23)		// don't compute 4-torsion for p smaller than this
#define ECURVE_8TOR_MINP			(1<<26)		// don't check any 8-torsion for p smaller than this
#define ECURVE_MOD5_MINP			(1L<<32)		// don't use 5-torsion data for p smaller than this

	
// Jacobian coordinates for elliptic curves, represents the affine point s(x/z,y/z), in Mumford rep: u(t)=t-x, v(t)=y.
struct ecp_j_struct {
	ff_t x, y, z;		
};
typedef struct ecp_j_struct ecp_j_t;
	
// reduced Chudnovsky Jacobian coordinate representation for genus 1
// Represents the affine point (x/z2,y/z3), or in Mumford rep: u(t)=t-x, v(t)=y.
struct ecp_jc_struct {
	ff_t x, y, z2, z3;				// note that we don't maintain z, this saves a field multiplication when composing elements
};
typedef struct ecp_jc_struct ecp_jc_t;


// converts affine coord to reduced Chudnovsky Jacobian coords
static inline void ecurve_A_to_JC (ecp_jc_t *p, ff_t x, ff_t y)
	{ _ff_set(p->x,x); _ff_set(p->y,y); _ff_set_one(p->z2); _ff_set_one(p->z3); }

// converts reduced Chudnovsky Jacobian to affine coords, given the inverse of z3, overlap ok
static inline void ecurve_JC_to_A (ff_t *px, ff_t *py, ecp_jc_t *p, ff_t z3inv)
	{ register ff_t t0,t1;  ff_mult(*py,p->y,z3inv);  _ff_mult(t0,p->z2,z3inv); _ff_square(t1,t0);  ff_mult(*px,p->x,t1); }
void ecurve_JC_to_A_parallel (ff_t x[], ff_t y[], ecp_jc_t p[], int n);

static inline void ecurve_JC_set_A (ecp_jc_t *p, ff_t x, ff_t y) { _ff_set(p->x,x); _ff_set(p->y,y); _ff_set_one(p->z2); _ff_set_one(p->z3); }
static inline int ecurve_JC_id (ecp_jc_t *p) { return _ff_zero(p->z2); }
static inline void ecurve_JC_set_id (ecp_jc_t *p) { _ff_set_zero(p->x);   _ff_set_zero(p->y);   _ff_set_zero(p->z2);   _ff_set_zero(p->z3); }
static inline int ecurve_JC_2tor (ecp_jc_t *p) { return _ff_zero(p->y); }
static inline void ecurve_JC_invert (ecp_jc_t *p) { ff_negate(p->y); }
static inline void ecurve_JCJC (ecp_jc_t *p1, ecp_jc_t *p2, ff_t f1);

// Sets s[i] = s[0]+i*(x0,y0) for i < n.  Returns k<n if s[k] is the identity (in which case no further elements are computed), otherwise returns n.
int ecurve_AJC_steps (ecp_jc_t s[], ff_t x0, ff_t y0, int n, ff_t f1);

// Computes p=n*(x0,y0) where (x0,y0) !=1 is in affine coordinates and p is in reduced Chudnovsky Jacobian coordinates
void ecurve_AJC_exp_ui (ecp_jc_t *p, ff_t x0, ff_t y0, unsigned long n, ff_t f1);

// Computes p = n*a where p and a are both in reduced Chudnovsky Jacobian coordinates
void ecurve_JC_exp_ui (ecp_jc_t *p, ecp_jc_t *a, unsigned long n, ff_t f1);

// There are many more JC and AJC related functions used internally to ecurve.c that we could easily export here (and will, as needed)

// scalar multiplication in affine coords (converts to JC internally)  returns 0 if the result is the identity, 1 ow.
int ecurve_exp_ui (ff_t *x2, ff_t *y2, ff_t x1, ff_t y1, unsigned long n, ff_t f1);
	
long ecurve_order_F3 (long *pd, ff_t f[4]);
long ecurve_order (long *pd, ff_t f[4]);
long ecurve_prime_order (ff_t f[4], int low);												// low=1 searches only for order < p
int ecurve_group_structure (long n[2], long N, long d, ff_t f[4]);
void ecurve_p_basis (ecp_jc_t *b1, long *q1, ecp_jc_t *b2, long *q2, long p, long N, ff_t f[4]);		// computes a basis (b1,b2) for the p-Sylow subgroup, given the group order N
int ecurve_test_exponent (long e, ff_t f[4]);
int ecurve_fast_trace_sign (ff_t f[4], long t);
//int ecurve_test_order (ppf_t n, ff_t f[4], int maxtest);
//int ecurve_test_order2 (ppf_t n[2], ff_t f[4]);

long ecurve_fastorder (ff_t x, ff_t y, long e, ff_t f1);			// compute the order of (x,y) given multiple e of the order (e.g. group order), f1=A


int ecurve_3tor (ff_t f[4]);							// f must be of the form x^3+ax+b
int ecurve_4tor (int *o, ff_t f[4], int flag8);				// f must be of the form x^3+ax+b (flag8 set indicates that Z/8Z should also be checked)
int ecurve_4tor_test (ff_t *tt, ff_t f[4]);					// returns true if y^2=f(x) has a 4-torsion pt, false otherwise.  If specified, tt must point to a root of f
int ecurve_mod3 (ff_t f[4], int *m);						// if m is non-null, the return value is mod *m, where *m is 3 or 3^2
int ecurve_mod4 (ff_t f[4], int odd_flag);					// returns the order of the elliptic curve y^2=f(x) modulo 4, caller may set odd_flag to 1 if it is already known that the group order is odd
int ecurve_mod5 (int a[2], int *m, ff_t f[4]);				// computes info on #E(Fp) mod 5.  a[] is an array of possible values mod *m, where *m is 5 or 5^2.  returns number of possible a's (1 or 2), or 0 if no info is easily available
int ecurve_halve (ff_t *x, ff_t *y, ff_t f[4]);				// Replaces P1=(x,y) with a point P2 such that 2P2 = +/-P1 (returns 0 if no such point exists)
int ecurve_halve_x (ff_t x1[1], ff_t x0, ff_t f[4]);			// Given the x-coord x2 of a point P2 on f(x)=x^3+f1x+f0, computes the (possible) x-coord x1 of a point P1 which, when doubled, yields P2 (does not verify P1 is on curve)
int ecurve_verify_2depth (ff_t x1, int k, int min_flag, ff_t x0, ff_t f[4]);		// Given the *unique* root x0 of f(x)=x^3+Ax+B and the x-coord x1!=x0 of a point P1, verifies that x1 can be halved exactly k times  (but no more).
int ecurve_2Sylow (ff_t f[4]);							// returns size of the 2-Sylow subgroup of the elliptic curve y^2=f(x)=x^3+f1x+f0, provided it is cyclic, returns -1 otherwise.

static inline int ecurve_to_jinv (ff_t j[1], ff_t f[4])
{
	ff_t t1, t2, t3;
	
	_ff_square(t1,f[1]);  ff_mult(t1,t1,f[1]);  _ff_x2(t1);  _ff_x2(t1);								// t1 = 4f1^3
	_ff_square(t2,f[0]); _ff_set_ui(t3,27); ff_mult(t2,t2,t3);									// t2 = 27f2^2
	_ff_add(t3,t1,t2);
	if ( _ff_zero(t2) ) return 0;
	ff_invert(t2,t3);
	_ff_set_ui(t3,1728);
	ff_mult(t1,t1,t3);
	_ff_mult(j[0],t1,t2);																// j = 1728 * 4f1^3 / (4f1^3+27f2^2)
	return 1;
}

// computes 2 j-invariants in parallel, saving an inversion.  Curves are assumed to be non-singular
static inline void ecurve_to_jinv2 (ff_t j1[1], ff_t j2[1], ff_t f1[4], ff_t f2[4])
{
	ff_t t1,t2,t3,num1,num2,den1,den2,c;
	
	_ff_set_ui(c,27); 
	_ff_square(t1,f1[1]);  _ff_mult(num1,t1,f1[1]);  _ff_x2(num1);  _ff_x2(num1);					// num1 = 4f1[1]^3
	_ff_square(t1,f1[0]); _ff_mult(t2,t1,c);  _ff_add(den1,num1,t2);							// den1 = (4f1[1]^3+27f1[0]^2)
	_ff_square(t1,f2[1]);  _ff_mult(num2,t1,f2[1]);  _ff_x2(num2);  _ff_x2(num2);					// num2 = 4f2[1]^3
	_ff_square(t1,f2[0]);  _ff_mult(t2,t1,c); _ff_add(den2,num2,t2);							// den2 = (4f2[1]^3+27f2[0]^2)
	_ff_mult(t2,den1,den2);
	_ff_invert(t1,t2); _ff_set_ui(c,1728);
	_ff_mult(t2,t1,den2); _ff_mult(t3,num1,c); _ff_mult(j1[0],t2,t3);							// j1 = 1728 * 4f1^3 / (4f1^3+27f1^2)
	_ff_mult(t2,t1,den1); _ff_mult(t3,num2,c); _ff_mult(j2[0],t2,t3);							// j2 = 1728 * 4f2^3 / (4f2^3+27f2^2)
}

static inline void ecurve_from_jinv (ff_t f[4], ff_t j[1])
{
	ff_t t0, t1;
	
	_ff_set_one(f[3]);
	_ff_set_zero(f[2]);
	if ( _ff_zero(j[0]) ) { _ff_set_zero(f[1]); _ff_set_one(f[0]); return; }
	_ff_set_ui(t1, 1728);
	if ( _ff_equal(j[0],t1) ) { _ff_set_one(f[1]); _ff_set_zero(f[0]); return; }
	_ff_sub(t0,t1,j[0]);
	ff_invert(t1,t0);
	_ff_mult(t0,t1,j[0]);
	_ff_add(f[0],t0,t0);
	_ff_add(f[1],f[0],t0);
}

static inline void ecurve_from_jinv_trace (ff_t f[4], ff_t J, long t)
{
	int sign;
	
	ecurve_from_jinv(f,&J);
	sign = ecurve_fast_trace_sign (f,t);
	if ( ! sign ) { printf ("Error, fast trace sign failed\n"); exit(0); }
	if ( sign < 0 ) ff_poly_twist(f,f,3);
}

unsigned long ui_randomb (unsigned long b);		// declared in mpzutil.h but we want to avoid including this

// generate a random affine point on y^2=x^3+f1*x+f0
static inline int ecurve_random_point (ff_t px[1], ff_t py[1], ff_t f[2])
{
	ff_t x, y, z;
	
	do {
		_ff_random(x);
		_ff_square (y, x); _ff_addto (y, f[1]); _ff_mult (y, y, x); _ff_addto (y, f[0]);
	} while ( ! ff_sqrt (&z, &y) );
	if ( ui_randomb(1) ) ff_negate (z);							// flip square root 50% of the time
	_ff_set(px[0],x);
	_ff_set(py[0],z);
	return 1;
}


// doubles a point p1 in reduced Chudnovsky Jacobian coords (p3 is the output, may be equal to p1)
// This code requires one more multiplication than doubling in standard Jacobian coordinates (but two fewer additions, which makes it a close call)
// Surprisingly, it is actually slightly faster than doubling in Jacobian coords when tested on an AMD Athlon 64 (YMMV).
static inline void ecurve_2JC (ecp_jc_t *p3, ecp_jc_t *p1, ff_t f1)
{
	register ff_t a, b, c, t0, t1, t2;
	
	_ff_square(t0,p1->x); _ff_add(t1,t0,t0); _ff_addto(t1,t0); 					// t1 = 3x^2
	_ff_square(t2,p1->z2);  _ff_mult(t0,f1,t2); _ff_add(b,t1,t0);					// b = 3x^2+f1*z^4	(note that f1=a4 in 13.2.1.c  of HECHECC p.282)
	_ff_add(c,p1->y,p1->y); _ff_square(t2,c); _ff_mult(a,t2,p1->x);				// a=4xy^2, c=2y, t2=4y^2
	_ff_add(t0,a,a); _ff_square(t1,b); _ff_sub(p3->x,t1,t0);						// x = b^2-2a
	_ff_mult(p3->z2,p1->z2,t2); _ff_mult(t1,c,t2); _ff_mult(p3->z3,p1->z3,t1);		// z2=4y^2z2, z3=8y^3z3, t1=8y^3
	_ff_mult(c,t1,p1->y);												// c = 8y^4
	_ff_sub(t0,a,p3->x); _ff_mult(t2,t0,b); _ff_sub(p3->y,t2,c);					// y = b(a-x)-c   -- note we use the new x value here
	// 11M+8A
}

// sums p1 in reduced Chudnovsky Jacobian coords and p2 in reduced Chudnovsky Jacobian coords  and puts the result in p1.
// Handles identity and will double if needed.  Cost is 13M+7A, one less than with standard Chudnovsky Jacobian coords.
static inline void ecurve_JCJC (ecp_jc_t *p1, ecp_jc_t *p2, ff_t f1)
{
	register ff_t a, b, c, d, e, f, t0, t1, t2;

	// we need to handle identity cases in general (although in some cases we might be able to rule this out)
	if ( _ff_zero(p2->z2) ) return;
	if ( _ff_zero(p1->z2) ) { *p1=*p2; return; }
	_ff_mult(a,p1->x,p2->z2);  _ff_mult(b,p2->x,p1->z2);  _ff_sub(e,b,a);			// a = x1z2^2, b = x2z1^2, e = a-b
	_ff_mult(c,p1->y,p2->z3);  _ff_mult(d,p2->y,p1->z3);  _ff_sub(f,d,c);			// c = y1z2^3, d = y2z1^3, f = d-c
	if ( _ff_zero(e) && _ff_zero(f) ) { ecurve_2JC (p1,p1,f1); return; }				// must double if pts are equal (11M+8A), inverses will end up with z2=0, so we let them through
	_ff_square(t1,e); ff_mult(t0,p1->z2,p2->z2);  _ff_mult(p1->z2,t0,t1);			// z^2 = z1^2z2^2e^2, t1=e^2
	_ff_mult(t2,t1,e); ff_mult(t0,p1->z3,p2->z3);  _ff_mult(p1->z3,t0,t2);			// z^3 = z1^3z2^3e^3, t2=e^3
	ff_mult(t1,t1,a);													// t1=ae^2
	_ff_square(t0,f); _ff_sub(p1->x,t0,t2); _ff_add(t0,t1,t1); _ff_subfrom(p1->x,t0);	// x = f^2-e^3-2ae^2
	_ff_sub(t0,t1,p1->x); _ff_neg(t1,t2); _ff_sum_2_mults(p1->y,f,c,t1,t0);			// y = f(ae^2-x) - ce^3 -- note we use the new x here
	// 13M+7A
}

#ifdef __cplusplus
}
#endif

#endif
