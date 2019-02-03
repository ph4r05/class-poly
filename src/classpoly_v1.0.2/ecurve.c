#define TEST_P 	1  	// only used when debugging code is uncommented, will cause debug info to be printed for primes above TEST_P

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <memory.h>
#include <gmp.h>
#include "ff_poly.h"
#include "mpzutil.h"
#include "ecurve.h"
#include "cstd.h"

/*
    Copyright (c) 2008-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

/*
	This module implements operations for genus 1 curves, both low level
	point operations and functions for exponentiation computing the
	group order, group structure, etc...
	
	These functions are heavily optimized and responsible for a more than a 2x
	improvement in genus 1 performance from  smalljac version 3 over version 2.
	
	We begin with basic point arithmetic functions (mostly inlined), the higher
	level functions appear below.  	At the end of the module are a collection of functions
	that were implement and test but are not used because they were found to be
	sub-optimal in the current application.
*/

long ecurve_JC_pp_order (ecp_jc_t *a, long p, long q, ff_t f1);
int ecurve_dlog (ecp_jc_t *a, ecp_jc_t *b, long o, ff_t f1);
int ecurve_bsgs_search (long *exp, ecp_jc_t *b, long low, long high, int m, int a1, int a2, int pflag, ff_t f1);
long ecurve_fastorder (ff_t x, ff_t y, long k, ff_t f1);
int ecurve_fastorder2 (ppf_t n, ecp_jc_t a[1], ppf_t e, int verify, ff_t f1);
int ecurve_test_exponent (long e, ff_t f[4]);
void ecurve_p_reduce (ecp_jc_t *b1, long *q1, ecp_jc_t *b2, long *q2, long p, ff_t f1);

unsigned long hecurve_expbits;
unsigned long hecurve_steps;
unsigned long hecurve_retries;

/*
	We use a reduced form of the Chudnovsky Jacobian representation (JC) which uses (x,y,z^2,z^3) to represent the affine point (x/z^2,y/z^3), but does not maintain z.
	Not computing z saves 1 field multiplication when compsing an affine point with JC point (usually called mixed addition).  Squaring (doubling) costs an extra multiplication
	but two fewer additions and the performance is comparable (in fact, when tested on an AMD Athlon 64, JC doubling was slighly faster, not clear why).

	In terms of the Mumford representation used elsewhere, the point point (x,y,z^2,z^3) is u[1]=z, u[0]=-x/z^2, v[0]=y/x^3. (notice the sign on x)
        When z is 1 this corresponds to u(t)=(t-x) and v(t)=y, so that x is a root of u and y=v(x).

	Negating y inverts an element, as in the affine rep, and 2-torsion elements have y=0.
	Any element with z=0 (or z2 for reduced Chudnovsky) is considered the identity, regardless of the value of x and y.

	As we use multiplicative notation for group operations elsewhere, we generally speak of the elliptic curve group operation multiplicatively,
        so we multiply, square, and exponentiate points, rather than adding, doubling, or multiplying by a scalar.

	In our operation counts we don't distinguish field multiplication and squaring.  For prime fields where p fits in a machine word, they are effectively the same.
        We do count additions, as these are not negligible--roughly speaking, 1M is about 2.5A.  Field inversions are horribly expensive relative to Montgomery multiplication,
        costing 40M or more (for p~2^30, say).  Whenever possible, we do n field inversions in parallel for a cost of I+3(n-1)M, effectively 3M per inverted element, for n large.
*/

// We trust the compiler to be smart about register allocation and don't hassle with macros (the speed difference appears to be negligible in testing)
// It would be nice to have C++ reference data types here...


static inline int ecurve_JC_cmp (ecp_jc_t *p1, ecp_jc_t *p2)
{
	register ff_t t0,t1;
	
	_ff_mult(t0,p1->x,p2->z2); _ff_mult(t1,p2->x,p1->z2);
	if ( ! _ff_equal(t0,t1) ) return 0;
	_ff_mult(t0,p1->y,p2->z3); _ff_mult(t1,p2->y,p1->z3);
	return ( _ff_equal(t0,t1) ? 1 : -1 );
}

// All points in p[] must be non-trivial, otherwise an attempt to invert 0 will result
void ecurve_JC_to_A_parallel (ff_t x[], ff_t y[], ecp_jc_t p[], int n)
{
	register int i;
	for ( i = 0 ; i < n ; i++ ) _ff_set(y[i],p[i].z3);
	ff_parallel_invert(x,y,n);
	for ( i = 0 ; i < n ; i++ ) ecurve_JC_to_A (x+i, y+i, p+i, x[i]);
}

// squares (doubles) an affine point into reduced Chudnovsky Jacobian coords
static inline void ecurve_2AJC (ecp_jc_t *p, ff_t x, ff_t y, ff_t f1)
{
	register ff_t t0,t1,t2,a,b;

	_ff_add(t0,y,y);  _ff_mult(t1,t0,y); _ff_add(p->z2,t1,t1); _ff_mult(a,x,p->z2);		// a=4xy^2, t1=2y^2, z2=4y^2
	_ff_mult(p->z3,t0,p->z2);											// z3=8y^3
	_ff_square(t1,x); _ff_add(b,t1,t1); _ff_addto(b,t1); _ff_addto(b,f1);				// b = 3x^2+f1*1^4
	_ff_square(t0,b); _ff_add(t2,a,a); _ff_sub(p->x,t0,t2);						// x = b^2-2a
	_ff_sub(t1,a,p->x); _ff_mult(t2,t1,b); _ff_mult(t0, y,p->z3); _ff_sub(p->y,t2,t0);	// y = b(a-x)-8y^4 (note we use the new x here and new z3=8y^3)
	// 7M+9A
}

// squares (doubles) an affine point into reduced Chudnovsky Jacobian coords and also sets az4 to f1*z3^4=16y^4 (cached value for squaring used below)
static inline void ecurve_2AJC_cache (ecp_jc_t *p, ff_t x, ff_t y, ff_t f1, ff_t *az4)
{
	register ff_t t0,t1,t2,a,b;

	_ff_add(t0,y,y);  _ff_mult(t1,t0,y); _ff_add(p->z2,t1,t1); _ff_mult(a,x,p->z2);		// a=4xy^2, t1=2y^2, z2=4y^2
	_ff_mult(p->z3,t0,p->z2);											// z3=8y^3
	_ff_square(t1,x); _ff_add(b,t1,t1); _ff_addto(b,t1); _ff_addto(b,f1);				// b = 3x^2+f1*1^4
	_ff_square(t0,b); _ff_add(t2,a,a); _ff_sub(p->x,t0,t2);						// x = b^2-2a
	_ff_sub(t1,a,p->x); _ff_mult(t2,t1,b); _ff_mult(t0, y,p->z3); _ff_sub(p->y,t2,t0);	// y = b(a-x)-8y^4 (note we use the new x here and new z3=8y^3)
	_ff_add(*az4,t0,t0); 
	// 7M+10A
}


// same as above except the parameter az4 contains a_4*z1^4 = f1*z1^4 and is updated to hold a_4*z3^4 (cost 1M less, 1A more)
static inline void ecurve_2JC_cache (ecp_jc_t *p3, ecp_jc_t *p1, ff_t *az4)
{
	register ff_t a, b, c, t0, t1, t2;
	
	_ff_square(t0,p1->x); _ff_add(t1,t0,t0); _ff_addto(t1,t0);  _ff_add(b,t1,*az4);	// b = 3x^2+f1*z^4
	_ff_add(c,p1->y,p1->y); _ff_square(t2,c); _ff_mult(a,t2,p1->x);				// a=4xy^2, c=2y, t2=4y^2
	_ff_add(t0,a,a); _ff_square(t1,b); _ff_sub(p3->x,t1,t0);						// x = b^2-2a
	_ff_mult(p3->z2,p1->z2,t2); _ff_mult(t1,c,t2); _ff_mult(p3->z3,p1->z3,t1);		// z2=4y^2z2, z3=8y^3z3, t1=8y^3
	_ff_mult(c,t1,p1->y); _ff_add(t0,c,c); ff_mult(*az4,*az4,t0);					// c = 8y^4, az4 = az4*16y^4
	_ff_sub(t0,a,p3->x); _ff_mult(t2,t0,b); _ff_sub(p3->y,t2,c);					// y = b(a-x)-c   -- note we use the new x value here
	// 10M+9A
}

// multiplies (adds) a non-identity point p1 in reduced Chudnovsky Jacobians coords by a (non-identity) affine point and puts result in p
// using the reduced Chudnosky form (no z, only z^2 and z^3) saves one field mult (this is faster than any of the alternatives in table 13.3 of HECHECC on p.284)
static inline void ecurve_AJC (ecp_jc_t *p, ecp_jc_t *p1, ff_t x0, ff_t y0, ff_t f1)
{
	register ff_t a, c, e, f, t0, t1, t2;
	
	if ( _ff_zero(p1->z2) ) { _ff_set(p->x,x0); _ff_set(p->y,y0); _ff_set_one(p->z2); _ff_set_one(p->z3); return; }
	_ff_mult(a,x0,p1->z2);											// a = x0z^2, b = x*1^2=x (z0=1)
	_ff_sub(e,p1->x,a);												// e = a-b
	_ff_mult(c,y0,p1->z3);											// c = y0z^3, d=y*1^3=y (z0=1)
	_ff_sub(f,p1->y,c);												// f = d-c
	if ( _ff_zero(e) && _ff_zero(f) ) { ecurve_2AJC (p, x0, y0, f1); return; }		// must use doubling code here, but at least it's an affine double (7M+9A) (inverses are ok)
	_ff_square(t0,e); _ff_mult(t1,t0,e); _ff_mult(t2,t0,a);					// e2=e^2, t1=e^3, t2=ae^2
	ff_mult(p->z2,p1->z2,t0); ff_mult(p->z3,p1->z3,t1);					// z2 = 1*z2*e^2, z3=1*z3*e^3
	_ff_square(t0,f); _ff_sub(p->x,t0,t1); _ff_add(t0,t2,t2); _ff_subfrom(p->x,t0);	// x = f^2-e^3-2ae^2
	_ff_sub(t0,t2,p->x); _ff_neg(t2,t1); _ff_sum_2_mults(p->y,f,c,t2,t0);		// y = f(ae^2-x) - ce^3
	// 10M+6A (9 redc)
}


/*
	Computes p=(x,y)^n where (x,y) !=1 is in affine coordinates and p is in Jacobian coordinates
	It takes advantage of the fact that all additions are of the form J+A, requiring only 11M, doubling is done in J+J form, using 10M
*/
void ecurve_AJC_exp_ui (ecp_jc_t *p, ff_t x0, ff_t y0, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	register ff_t negy0;
	int i;

	if ( n == 0 ) { _ff_set_zero(p->z2); return; }
	if ( n == 1 ) { _ff_set(p->x,x0); _ff_set(p->y,y0); _ff_set_one(p->z2); _ff_set_one(p->z3);  return; }
	ecurve_2AJC(p,x0,y0,f1);
	if ( n == 2 ) return;
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits)-2;						// we know the top two bits of the NAF are 10
	_ff_neg(negy0,y0);
	m = (1UL<<i);
	for ( ; m ; m >>= 1 ) {
		ecurve_2JC(p,p,f1);	 					// 11M+8A
		if ( m&pbits ) ecurve_AJC(p,p,x0,y0,f1);		// 10M+6A
		if ( m&nbits ) ecurve_AJC(p,p,x0,negy0,f1);	// 10M+6A
	}
}

/*
	Computes p=a^e where a!=1 is in reduced Jacobian Chudnovsky coords, and so is p. overlap is ok.
*/
void ecurve_JC_exp_ui (ecp_jc_t *p, ecp_jc_t *a, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	ecp_jc_t t, ai;
	int i;
	
	if ( n == 0 ) { _ff_set_zero(p->z2); return; }
	if ( n == 1 ) { *p=*a; return; }
	ecurve_2JC(&t,a,f1);							// 11M+8A

	if ( n == 2 ) { *p = t; return; }
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits)-2;
	ai=*a;  ecurve_JC_invert(&ai);
	m = (1UL<<i);
	for ( ; m ; m >>= 1 ) {
		ecurve_2JC(&t,&t,f1);					// 11M+8A
		if ( m&pbits ) ecurve_JCJC(&t,a,f1);			// 13M+7A
		if ( m&nbits ) ecurve_JCJC(&t,&ai,f1);		// 13M+7A
	}
	*p = t;
}

// scalar mult in affine coords, returns 0 if the result is the identity, 1 ow
int ecurve_exp_ui (ff_t *x2, ff_t *y2, ff_t x1, ff_t y1, unsigned long n, ff_t f1)
{
	ecp_jc_t p[1];
	ff_t zinv;
	
	ecurve_AJC_exp_ui (p, x1, y1, n, f1);
	if ( ecurve_JC_id(p) ) return 0;
	ff_invert (zinv, p->z3);
	ecurve_JC_to_A (x2, y2 ,p, zinv);
	return 1;
}



// Combines precomputed values p[i]=p[0]^(2^i) to compute p[0]^n.  Assumes all required powers are present: NAF potentially requires one more bit!!
// CAUTION: if this is false, the resulting behavior is very strange and unpredictable.  You have been warned...
static inline void ecurve_JC_exp_powers(ecp_jc_t *o, ecp_jc_t *p, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	int i;

	// handle small/common n quickly
	switch (n) {
	case 0:	_ff_set_zero(o->z2); return; 
	case 1:	*o=p[0]; return;
	case 2:	*o=p[1]; return;
	case 3:	*o=p[0]; ecurve_JCJC(o,p+1,f1); return;
	case 4:	*o=p[2]; return;
	case 5:	*o=p[2]; ecurve_JCJC(o,p,f1); return;
	case 6:	*o=p[2]; ecurve_JCJC(o,p+1,f1); return;
	case 7:	*o=p[0]; ecurve_JC_invert(o); ecurve_JCJC(o,p+3,f1); return;
	case 8:	*o=p[3]; return;
	case 9:	*o=p[3]; ecurve_JCJC(o,p,f1); return;
	case 12:	*o=p[3]; ecurve_JCJC(o,p+2,f1); return;
	case 36:	*o=p[5]; ecurve_JCJC(o,p+2,f1); return;
	case 60:	*o=p[2]; ecurve_JC_invert(o); ecurve_JCJC(o,p+6,f1); return;
	case 180: *o=p[7]; ecurve_JCJC(o,p+5,f1); ecurve_JCJC(o,p+4,f1); ecurve_JCJC(o,p+2,f1); return;
	}
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits);
	*o = p[i];
	i-=2;
	m = (1UL<<i);
	for ( ; m ; m >>= 1, i-- ) {
		if ( (m&pbits) ) ecurve_JCJC(o,p+i,f1);								// 13M+8A
		if ( (m&nbits) ) {													// 13M+10A (slightly more expensive)
			ecurve_JC_invert(o); ecurve_JCJC(o,p+i,f1); ecurve_JC_invert(o);
		}
	}
}

// Sets p[i] = p[0]^(2^i) for i in [0,k]. Input is an affine pt not equal to the identity, output is in reduced Jacobian Chudnovsky coords.
void ecurve_AJC_powers (ecp_jc_t p[], ff_t x0, ff_t y0, int k, ff_t f1)
{
	register ecp_jc_t *e;
	ff_t az4;
	
	_ff_set(p->x,x0);  _ff_set(p->y,y0);  _ff_set_one(p->z2); _ff_set_one(p->z3);
	ecurve_2AJC_cache (p+1, x0, y0, f1, &az4);						// 7M + 10A
	for ( e=p+k,p+=2 ; p <= e ; p++ ) ecurve_2JC_cache(p,p-1,&az4);	// 10M + 9A
}

// Sets p[i] = p[0]^(2^i) for i in [1,k]. Input is a point not equal to the identity, input and output is in reduced Jacobian Chudnovsky coords.
void ecurve_JC_powers (ecp_jc_t p[], int k, ff_t f1)
{
	register ecp_jc_t *e;

	// don't bother caching
	for ( e=p+k,p++ ; p <= e ; p++ ) ecurve_2JC(p,p-1,f1);		 	// 11M + 8A
}

// Sets s[i] = s[0]*(x0,y0)^i for i < n.  Returns k<n if s[k] is the identity (in which case no further elements are computed), otherwise returns n.
int ecurve_AJC_steps (ecp_jc_t s[], ff_t x0, ff_t y0, int n, ff_t f1)
{
	register ecp_jc_t *end;
	
	if ( ecurve_JC_id(s) ) return 0;
	end = s+n;
	for ( s++ ; s < end ; s++ ) {
		ecurve_AJC (s,s-1,x0,y0,f1);								// 10M + 6A
		if ( ecurve_JC_id(s) ) break;
	}
	hecurve_steps += n-(end-s);
	return n-(end-s);
}

// Sets s[-i] = s[0]*(x0,-y0)^i for i < n.  Returns k<n if s[k] is the identity (in which case no further elements are computed), otherwise returns n.
// Used for downward stepping - NOTE THAT s SHOULD POINT TO THE LAST ENTRY IN THE ARRAY, we're walking backwards here
int ecurve_AJC_dsteps_1 (ecp_jc_t *s, ff_t x0, ff_t y0, int n, ff_t f1)
{
	register ecp_jc_t *end;
	register ff_t negy0;
	
	if ( ecurve_JC_id(s) ) return 0;
	end = s-n;
	_ff_neg(negy0,y0);
	for ( s-- ; s > end ; s-- ) {
		ecurve_AJC (s,s+1,x0,negy0,f1);							// 10M + 6A
		if ( ecurve_JC_id(s) ) break;
	}
	hecurve_steps += n-(s-end);
	return n-(s-end);
}

// sets s[2*i+j] = s[0]*(x0,y0)^(i+j)*(x1,y1)^i for 2*i+j < n with j=0 or 1.  Returns k<n if s[k] is the identity (in which case no further elements are computed), otherwise returns n.
int ecurve_AJC_steps_2 (ecp_jc_t s[], ff_t x0, ff_t y0, ff_t x1, ff_t y1, int n, ff_t f1)
{
	register ecp_jc_t *end;
	register int i;
	
	if ( ecurve_JC_id(s) ) return 0;
	end = s+n;
	for ( i=1,s++ ; s < end ; i++,s++ ) {
		if ( i&1) {
			ecurve_AJC (s,s-1,x0,y0,f1);		// 10M + 6A
		} else {
			ecurve_AJC (s,s-1,x1,y1,f1);		// 10M + 6A
		}
		if ( ecurve_JC_id(s) ) break;
	}
	hecurve_steps += n-(end-s);
	return n-(end-s);
}

static inline int inv_mod3 (int e) { return ((e%3)==1?1:2); }

// handle F_3 as a special case, don't worry about speed.  Assumes f monic, but does check its discriminant (returns -1 for singular curves)
long ecurve_order_F3 (long *pd, ff_t f[4])
{
	ff_t x, y, t;
	long pts;
	
	
	// verify that disc(x^3+f2*x^2+f1*x+f0) = f2^2*(f1^2-f0*f2) - f1^3 != 0 in F_3, return -1 ow
	_ff_square (x, f[1]); _ff_mult (y, f[0], f[2]);  _ff_sub(t, x, y);  _ff_square (y, f[2]);  _ff_mult (t, t, y); _ff_mult (x, x, f[1]); if  ( _ff_equal (t, x) ) return -1;

	pts = 1;
	_ff_set_zero(x);  ff_poly_eval (&y,f,3,&x);  if ( _ff_zero(y) ) pts++; else if ( _ff_one(y) ) pts+=2;
	_ff_inc(x);  ff_poly_eval (&y,f,3,&x);  if ( _ff_zero(y) ) pts++; else if ( _ff_one(y) ) pts+=2;
	_ff_inc(x);  ff_poly_eval (&y,f,3,&x);  if ( _ff_zero(y) ) pts++; else if ( _ff_one(y) ) pts+=2;

	// if pd != 0, set *pd = gcd(m,6) where E(F_p) = Z/mZ x Z/nZ with m dividing n.  An exhaustive search finds that m=1 unless f(x) = x^3 + 2*x
	if ( pd ) *pd = ( _ff_zero(f[2]) && _ff_zero(f[0]) && _ff_equal (f[1], _ff_half) ? 2 : 1 );
	return pts;
}

/*
	Fast group order computation for elliptic curves y^2=f(x) over F_p, supports for 2 < p < 2^40, optimized for p~2^30.
	Assumes f monic and, for p >3, of the form x^3+f1*x+f0 with nonzero discriminant (so the curve is not singular).

	Returns the group order and if pd is non-null, sets *pd to gcd(m,6), where the group structure is Z/mZ x Z/nZ with m dividing n (possibly m=1).
	This information can be used to speed up group structure computations.
*/
long ecurve_order (long *pd, ff_t f[4])
{
	long r, d, e, E, min, max, low, high, exp, M;
	ecp_jc_t t[1];
	ff_t g[4],x,y,*h;
	int a,a1,a2,i,m,twist,s2known,tor3,tor4,flag8,sts;

	// TODO: This code needs to be rewritten to use mod3/mod4/mod5 functions
	
	assert (_ff_one(f[3]));
	if ( _ff_p == 3 ) return ecurve_order_F3 (pd, f);
	assert (_ff_zero(f[2]));
	
	/*
		We compute the 3-torsion subgroup first, because this information may lead us to work in the twist.
		The function ecurve_3tor() computes the 3-torsion subgroup of y^2=f(x), but if it is trivial, it will try to
		determine the 3-torsion subgroup of the twist, using the factorization pattern of the 3-division polynomial (a quartic).
	*/
	h = f;  twist = 0;
	tor3 = ecurve_3tor(f);	
	if ( tor3<0 ) { ff_poly_twist(g,f,3); h = g; twist = 1; tor3=-tor3; }		// negative return value indicates 3-torsion in the twist (but not in the primary)
	d = ( tor3==9 ? 3 : 1);										// d is a known divisor of |G|/lambda(|G|), i.e. it divides the group exponent (lambda(G)) and its square divides |G|.

	if ( _ff_p < ECURVE_4TOR_MINP ) {								// when p is small, only compute 2-torsion but not 4-torsion (the marginal benefit does not justify the cost)
		 i = ff_poly_roots_d3(0,h);
		if ( i == 0 ) {											// no roots means there is no 2-torsion and the group order must be odd
			s2known = 1;  tor4 = 1;								// the flag s2known indicates we know the entire 2-Sylow subgroup
		} else {												// this means that once we factor it out, the remaining exponent is known to be odd
			s2known = 0;  tor4 = 2;								// if there is one root, we have 2-torsion Z/2Z, otherwise there are 3 roots and 2-torsion Z/2Z x Z/2Z
			if ( i > 1 ) { tor4 = 4;  d *= 2; }			
		}		
	} else {
		s2known = 0;
		flag8 = ( _ff_p < ECURVE_8TOR_MINP ? 0 : 1);					// flag8 set if we also want to check whether the 2-Sylow subgroup is a cyclic group containing Z/8Z
		d *= ecurve_4tor(&tor4,h,flag8);							// compute the 4-torsion subgroup (and, optionally, check Z/8Z as well)
		if ( flag8 ) {
			if ( tor4 <= 4 ) s2known = 1;							// In these cases we know the entire 2-Sylow subgoup
		} else {
			if ( tor4 < 4 || (tor4==4 && d==2) ) s2known = 1;			// if the 4-torsion subgroup contains no elements of order 4, it must be equal to the 2-Sylow subgroup.
		}
	}
	e = tor4*tor3/d;											// e is our known divisor of lambda(|G|)
	E = d*e;													// E is our known divisor of |G|
	m = (s2known==1?2:1)*(tor3==1?3:1);							// we know that gcd(|G|/tor4,m) = 1, which we can use to speed up the BSGS search
	
	a = ( tor3==1 && _ff_p1mod3 ? -1 : 0 );							// if neither the curve or its twist has 3-torsion and p=1mod3, we must have a_p=0 and group order 2 mod 3
															// For any divisor x of |G| this means |G|/x must be equal to -1/x mod 3 (and also mod 6 if m=6)
	r = (long)(2.0*sqrt(_ff_p));
	exp = _ff_p+1;
	min = exp-r;
	max = exp+r;
	low = _ui_ceil_ratio(min,E);
	high = max/E;
	if ( low > high ) { printf ("Error, no multiple of torsion derived e=%ld exists in interval [%ld,%ld] for p=%ld\n", e, min, max, _ff_p); abort(); }
	if ( pd ) { *pd = ( twist && tor3==9 ? d/3 : d );  if ( !((*pd)&3) ) *pd/=2; }	// set *pd to reflect 2-torsion and 3-torsion information in y^2=f(x) (but not the twist)

//printf("%lu: tor3=%d, tor4=%d, e=%ld, d=%ld, low=%ld, high=%ld, min=%ld, max=%ld ", _ff_p, tor3, tor4, e, d, low, high, min, max); ff_poly_print(h,3);
	
	// Note that for a non-CM curve the probability that ecurve_bsgs succeeds on the first try is close to 1 and grows with p (say 99.9% for p~2^20)
	sts = 0;
	for ( i = 0 ; i < ECURVE_ORDER_RETRIES ; i++ ) {
		if ( low==high ) {exp=low; break; }
		if ( ! ecurve_random_point(&x,&y,h) ) continue;				// this can fail, e.g.  y^2=x^3+2x+2 over F_3 has no finite rational points
//printf ("%lu: e=%ld, Random point (%ld,%ld) on curve ", _ff_p, e, _ff_get_ui(x), _ff_get_ui(y));  ff_poly_print(h,3);
		ecurve_AJC_exp_ui (t, x, y, e, h[1]);  if ( ecurve_JC_id(t) ) continue;
		// handle order 2 elements here so that bsgs search can assume |t| > 2
		if ( ecurve_JC_2tor(t) ) { 
			exp = 2;
 		//On paper, the code below should speed things up (it uses no inversions and fewer operations for small intervals), but in testing it actually slows things down slightly
		//} else if ( high-low < ecurve_SHORT_INTERVAL ) {
		//	if ( ecurve_bsgs_short (&exp, &t, low, high, h[1]) ) break;
		} else {
			if ( a ) { a = -inv_mod3(E); if ( m==6 && a==-2 ) a = -5; }
			if ( a ) { a1 = m+a; a2 = 0; } else { a1 = 1; a2 = m-1; }
			if ( (sts=ecurve_bsgs_search (&exp, t, low, high, m, a1, a2, 0, h[1])) ) break;
		}
		e *= exp;  E *= exp;
		low = _ui_ceil_ratio(min,E);
		high = max/E;
		// there are lot's of optimizations we could insert here, but none of them apply often enough to be worth doing, so keep it simple.
		hecurve_retries++;
	}
	if ( sts < 0 ) { printf ("Error, BSGS search of entire Hasse-Weil interval failed for point(%ld/%ld,%ld/%ld), low=%ld, high=%ld, m=%d, e=%ld, E=%ld, p=%ld, curve ",
					 _ff_get_ui(t->x), _ff_get_ui(t->z2), _ff_get_ui(t->y), _ff_get_ui(t->z3), low, high, m, e, E, _ff_p); ff_poly_print(h,3);  abort(); }
	if (  i < ECURVE_ORDER_RETRIES || low==high ) {		// AVS 11/12/2008 fixed a bug that could arise when low==high occurs at the same time ecurve_ORDER_RETRIES is reached
		if ( low==high ) exp=low;
		E *= exp;
		if ( twist ) E = 2*(_ff_p+1)-E;
		if ( E < min || E > max ) { printf ("Error, computed order %ld is not in Hasse-Weil interval [%ld,%ld] for p=%ld, h= ", E, min, max, _ff_p); ff_poly_print(h,3); abort(); }
		return E;
	}
	if ( e > 2*r ) { printf("Error, exponent %ld for p=%ld has unique multiple in [%ld,%ld] (low=%ld,high=%ld) that was not detected for curve ", e, _ff_p, min, max, low, high);  ff_poly_print(h,3); abort(); }
	if ( e < sqrt(min) ) { printf ("Error, exponent %ld is impossibly small for p=%ld\n", e, _ff_p); abort(); }
	
	/*
		With very high probability (exponential in ecurve_ORDER_RETRIES), e is now the group exponent.   For all but 21 primes (the largest of which is 547)
		this uniquely determines the group order (this result is due to Cremona and Harley, or see Washington, "Elliptic Curves", Prop. 4.19).
		In fact, taking our knowledge of 2-torsion and 3-torsion into account, the group order is uniquely determined in every case (AVS: to be written up).
		Unfortunately, we can't prove that e is definitely the group exponent without computing the group structure.
	
		Instead, we apply results of Schoof and Mestre (see Washington, Prop. 4.18) which tells us that either the curve or its twist contains an element whose
		order has a unique multiple in the Hasse interval provided p > 229.  If one examines the exceptional cases for p<=229, one finds that if we also
		consider the information provided by 2-torsion, in every case but two (y^2=x^3+2 over F_13 and F_19) we obtain a known divisor of the group order
		with a unique multiple in the Hasse interval, and the two exceptional cases are addressed if we consider 3-torsion. (AVS: to be written up).

		In the two 3-torsion cases, we will have chosen h to be the twist with 3-torsion and should have succeeded above.  If not we mustn't have computed
		the full group exponent and we should try again (by recursively calling ourselves below)
	
		For the 2-torsion cases, the twist will also have 2-torsion, we just need to use this information when we try to find a unique multiple of the exponent
		in the twist below.  If we fail, it must mean we got unlucky and should try again.
	
		In summary, this yields a Las Vegas algorithm with provably correct output, for all odd p (for p=2 the curve y^2=f(x) is singular).
	*/
	
	m = (d%2?1:2);		// m is a known divisor of the twist group order based on the rank of the 2-torsion subgroup, note that we don't keep track of d anymore and reuse it below

	for ( M = e*_ui_ceil_ratio(min,e) ; M <= max ; M += e ) { d = M/e;  if ( !(e%d) && !((_ff_p-1)%d) ) break; }
	if ( M > max ) { printf ("Error: no multiple M of ambiguous e=%ld satisfies M/e | gcd(e,p-1) (p=%ld).\n", e, _ff_p); abort(); }
	if ( ! twist ) { ff_poly_twist(g,f,3); h = g; } else { h = f; }	

	do {
		r = 2*(_ff_p+1)-M;
		// Attempt to prove that r is the order of the twist.  We don't bother trying to be efficient here, we expect to succeed on the first try.
		for ( i = 0 ; i < ECURVE_ORDER_RETRIES ; i++ ) {
			if ( ! ecurve_random_point(&x,&y,h) ) { printf("hecurve_random_point failed!"); abort(); }
			ecurve_AJC_exp_ui (t, x, y, r, h[1]);
			if ( ! ecurve_JC_id(t) ) break;											// can't be the right M
			exp = m*ecurve_fastorder(x,y,r,h[1]);									// compute the order of our random point, times our 2-torsion derived info
			low = exp*_ui_ceil_ratio(min,exp);										// low is the first multiple of exp >= min
			if ( low <= max && low+exp > max ) { return ( twist ? 2*(_ff_p+1)-M : M ); }	// success
		}
		// try another candidate, if there is one (this can happen for small p)
		for ( M+= e; M <= max ; M += e ) { d = M/e;  if ( !(e%d) && !((_ff_p-1)%d) ) break; }			
	} while ( M <= max );
	// Otherwise, try again -- we must have not gotten the full group exponent (this essentially never happens unless ecurve_ORDER_RETRIES is set quite low).
	return ecurve_order(pd,f);
}


/*
	Computes the order of the specified curve, assuming it is prime (which allows several optimizations).
	If low is set, only checks for order less than p
	Returns zero otherwise.
*/
long ecurve_prime_order (ff_t f[4], int flags)
{
	long r, min, max, exp;
	ecp_jc_t t[1];
	ff_t x,y;
	int a[2],a1,a2,i,m,n,mod4,sts,lowhalf;

	lowhalf = flags&1; flags >>= 1;
	if ( _ff_p < 31 ) { exp = ecurve_order (0, f); return ( ui_is_prime(exp) ? exp : 0 ); }		// let ecurve_order handle special cases for small p
	if ( _ff_poly_roots_d3 (0,f,0,0) ) return 0;										// make sure curve has odd order first (even though caller may have already verified that f doesn't split 1,2)
	a1 = ecurve_mod3(f,&m);
	if ( ! a1 ) return 0;
	a2 = 0;
	mod4 = ecurve_mod4(f,1);
	while ( (a1&3) != mod4 ) a1+=m;
	m *= 4;
	if ( _ff_p >= ECURVE_MOD5_MINP ) {
		i = ecurve_mod5(a,&n,f);
		if ( i==1 ) {
			if ( ! a[0] ) return 0;
			while ( (a1%n) != a[0] ) a1 += m;    m *= n;
		} else if ( i==2 ) {
			for ( a2 = a1 ; (a2%n) != a[1] ; a2 += m );
			while ( (a1%n) != a[0] ) a1 += m;    m *= n;
		}
	}
	r = (long)(2.0*sqrt(_ff_p));
	min = _ff_p+1-r;
	max = ( lowhalf ? _ff_p-1 : _ff_p+1+r );
	for ( i = 0 ; ! ecurve_random_point(&x,&y,f) &&  i < ECURVE_ORDER_RETRIES ; i++ );
	if ( i == ECURVE_ORDER_RETRIES ) return 0;						// if we cant' find a  non-trivial point, assume order is 1, which isn't prime
	if ( _ff_zero(y) ) return 0;									// handle order 2 elements here so that bsgs search can assume |t| > 2
	ecurve_A_to_JC(t,x,y);
	// note that since p > 29, there can only be one prime exponent in the Hasse interval
	sts = ecurve_bsgs_search (&exp, t, min, max, m, a1, a2, 1, f[1]);
//if ( _ff_p > TEST_P ) { printf ("bsgs returned sts=%d, exp=%ld\n", sts, exp); printf ("p=%ld, a1=%ld, a2=%ld, m=%ld, min=%ld, max=%ld\n", _ff_p,a1,a2,m,min,max); ff_poly_print(f,3); }
	if ( sts < 0 || exp < min || exp > max ) return 0;					 // exp should actually never be greater than max, but we won't hold bsgs to this
	return ( ui_is_prime(exp) ? exp : 0 );
}


/*
	Given the group order N, compute the isomorphism type of the group structure of an elliptic curve, of the form Z/n1Z x Z/n2Z with n1 dividing n2 (and also p-1), possibly n1=1.
	The parameter d specifies gcd(6,n1), which can be derived from 2-torsion and 3-torsion info by ecurve_order above (specify 0 if not known).
	The rank of the group (1 or 2) is returned, and n[0]=n2 if n1=1, o.w. n[0]=n1 and n[1]=n2.  Note that n2 is the group exponent.

	It is possible to compute the group invariants (n1 and n2) in probabilistic (Las Vegas) polynomial time using the Weil pairing via Miller's algorithm
	(although this doesn't actually determine a basis, i.e. independent elements of order n1 and n2).

	We take a simpler generic approach which is very fast in the typical case.  We may in the future want to invoke (refinements of) Miller's algorithm for
	hard cases (large prime dividing p-1 whose square divides N).  As it stands, the algorithm below has complexity O(q^(1/2)) where q is the largest prime dividing p-1
	whose square divides N.  This is O(N^(1/4)) in the worst case, but this rarely happens and in practice computing the group structure takes about 10% the
	time it takes to compute the group order (for a non-CM curve).  (It might be interesting to do a performance comparison wth Miller's algorithm at some point)
	
	We just return the group invariants here, and optimize for computing these as quickly as possible, but the algorithm below can easily be modified to return a basis
	without changing the complexity (but the constant factors will be slightly worse).  This is a Las Vegas algorithm, i.e. provably correct output and bounded expected running time.
*/
int ecurve_group_structure (long n[2], long N, long d, ff_t f[4])
{
	long n1,n2,e;
	ecp_jc_t b1, b2;
	unsigned long p[MAX_UI_PP_FACTORS], q[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS];
	ff_t x,y;
	long q0, q1, q2, r;
	int i, j, k;
	
	/*
		For each prime q|N, if q does not divide p-1 we know the q-Sylow subgroup is cyclic.  Even when q does divide p-1, if q^2 does not divide N, we
		again know the q-Sylow is cyclic (of prime order in this case).  Applying just these 2 facts will often suffice to determine the group structure without
		using any group operations, particularly if information in d is also applied (e.g. 2-torsion and 3-torsion info).
	*/
	k = ui_factor(p,h,ui_gcd(_ff_p-1,N));
	for ( i = j = 0 ; i < k ; i++ ) {												// make a list of primes p[i] dividing p-1 whose square divides N
		if ( p[i]==2 && d && d%2 ) continue;									// d=gcd(6,n1) not divisible by 2 implies 2-Sylow must be cyclic
		if ( p[i]==3 && d && d%3 ) continue;									// d=gcd(6,n1) not divisible by 3 implies 3-Sylow must be cyclic
		if ( !(N%(p[i]*p[i])) ) p[j++] = p[i];
	}
	n1 = 1;  n2 = N;
	for ( i = 0 ; i < j ; i++ ) {
		for ( q[i] = p[i]*p[i] ; !(N%(q[i]*p[i])) ; q[i] *= p[i] );						// determine power of p[i] dividing N
		n2 /= q[i];														// remove q[i] from n2
	}
//printf("%ld: N=%ld , j=%d, n2=%ld\n", _ff_p, N, j, n2); ff_poly_print(f,3);
	// Now determine the p[i]-Sylow subgroups for p[i] dividing p-1 and p[i]^2 dividing N (if there are none, we have no work to do).
	for ( i = 0 ; i < j ; i++ ) {
		if ( d ) for ( q0 = 1 ; ! (d%(p[i]*q0)) ; q0 *= p[i] );	else q0 = 1;				// if d was specified, use d to compute q0, a known divisor of q1 (and q2), possibly 1
		e = N/(q[i]*n1);													// the p[i]-Sylow subgroup must lie in the image of the e-power map (we can take out the n1 we know)
		for ( q1 = q2 = q0, k=0 ; q1*q2 < q[i] && k < 1000 ; k++ ) {				// bound the number of retries as a safety valve, just in case the group order was invalid
			if ( ! ecurve_random_point(&x,&y,f) ) continue;						// get a random finite point
			ecurve_AJC_exp_ui (&b1, x, y, e, f[1]);								// get an element of the p[i]-Sylow via the e-power map
			r = ecurve_JC_pp_order (&b1,p[i],q[i],f[1]);							// compute its order, which will be a power of p[i]
			if ( r <= q1 ) continue;											// if the order of b1 is not strictly larger than q1, it isn't going to help us
			if ( q2==q0 ) { b2=b1; q2 = r; continue; }							// first time through we will just set q2
			ecurve_p_reduce (&b1, &r, &b2, &q2, p[i], f[1]);						// reduce to a basis for <b1,b2>
			q1 = ( r > q0 ? r : q0);											// note that if r<q0, then q1 is not the order of b1, but we don't care
		}
		if ( k == 1000 ) { printf ("Group structure computation failed on %lu-Sylow of size %lu with group of order %lu on curve over F_%lu: ", p[i], q[i], N, _ff_p); ff_poly_print(f,3); abort(); }
		n1 *= q1;  n2 *= q2;
	}
	if ( n1*n2 != N ) { printf ("bug in ecurve_group_structure, %ld*%ld != %ld,  F_%lu curve ", n1, n2, N, _ff_p);  ff_poly_print(f,3); abort(); }		// sanity check
	if ( n1 == 1 ) { n[0] = n2; return 1; }
	n[0] = n1;  n[1] = n2;
	return 2;
}

int ecurve_fast_trace_sign (ff_t f[4], long t)
{
	ecp_jc_t P[1];
	ff_t x, y;
	long N1, N2;
	int i;
	
	N1 = _ff_p+1-t;
	N2 = _ff_p+1+t;
	for ( i = 0 ; i < 100 ; i++ ) {
		if ( ! ecurve_random_point(&x,&y,f) ) continue;
		ecurve_AJC_exp_ui(P,x,y,N1,f[1]);
		if ( ! ecurve_JC_id(P) ) return -1;
		ecurve_AJC_exp_ui(P,x,y,N2,f[1]);
		if ( ! ecurve_JC_id(P) ) return 1;
	}
	return 0;
}


/*
	Las Vegas algorithm for testing whether the factored integer n is equal to the group order.

	We rely on the caller to select the curve or its twist appropriately, and handle special
	cases for small p.  If maxtest is specified, an error is returned if the order has not been
	unambiguously tested after maxtest random elements.

	It is assumed that n lies in the Hasse interval
*/
int ecurve_test_order (ppf_t n, ff_t f[4], int maxtest)
{
	ecp_jc_t t[1];
	ppf_t d,e;
	ff_t x,y;
	long m, N, o, r, min, max;
	int k;
	
	r = (long)(2.0*sqrt(_ff_p));
	min = _ff_p+1-r;
	max = _ff_p+1+r;
	N = ppf_eval(n);
//printf("N=%ld="); for ( i = 0 ; i< n->w ; i++ ) printf("%d^%d ", n->p[i], n->h[i]); puts("");
	
	ppf_copy (e,n);
	m = 1;
	for ( k = 0 ; ; k++ ) {
		if ( ! ecurve_random_point(&x,&y,f) ) continue;
		ecurve_AJC_exp_ui(t,x,y,m,f[1]);
//printf ("Random point (%ld,%ld)^%d = (%ld,%ld,%ld,%ld)\n", _ff_get_ui(x), _ff_get_ui(y), m, _ff_get_ui(t->x), _ff_get_ui(t->y), _ff_get_ui(t->z2), _ff_get_ui(t->z3));
		if ( ! ecurve_fastorder2 (d, t,e,1,f[1]) ) return 0;
//printf ("Fast order returned "); for ( i = 0 ; i < d->w ; i++ ) printf ("%d^%d ", d->p[i], d->h[i]); puts("");
		ppf_divexact (e, d);
		m *= ppf_eval(d);
		o = m*_ui_ceil_ratio(min,m);								// o is the first multiple of m >= min	
		if ( o > max ) return 0;									// if o > max, n can't be the order
		if ( o+m > max ) return (o==N?1:0);						// if only one multiple in Hasse interval, we can determine things
		if ( maxtest && k >= maxtest ) { printf ("hit maxtest N=%ld, p=%ld, ", N, _ff_p); ff_poly_print(f,3); return 0; }
	}
	return 1;
}

int ecurve_test_order2 (ppf_t  n[2], ff_t f[4])
{
	ecp_jc_t P[1];
	ff_t ft[4],*g[2], *swap;
	ppf_t d[2],e[2], tmp,tmp2;
	ff_t x,y;
	long m[2],r,o,N[2],min,max,offset,m3,k3;
	register int i, k, swapped;
	
	r = (long)(2.0*sqrt(_ff_p));
	min = _ff_p+1-r;
	max = _ff_p+1+r;
	if ( max < _ff_p || max >= (1L << 62) - 2 ) { err_printf ("p=%ld too large in ecurve_test_order2!\n", _ff_p); abort(); }
	ff_poly_twist (ft,f,3);
	N[0] = ppf_eval(n[0]);
	N[1] = ppf_eval(n[1]);
	if ( N[0]+N[1] != 2*_ff_p+2 ) { printf ("invalid trace order %ld+%ld != 2p+2 for p=%ld\n", N[0], N[1], _ff_p); abort(); }
	ppf_copy (e[0],n[0]);
	ppf_copy (e[1],n[1]);
	m[0] = m[1] = 1;
	g[0] = f;  g[1] = ft;
	swapped = 0;
	for(i=0;;i++) {
		
		// this should never happen, but just in case...
		if ( i > 30 ) {
			ff_t j;
			if ( ! ecurve_to_jinv(&j,f) ) { err_printf ("singular curve in ecurve_test_order2!\n"); }
			else  { err_printf ("ecurve_test_order2 got stuck on p=%ld, j=%ld m=%ld mt=%ld ", _ff_p, _ff_get_ui(j), m[0], m[1]); }
			ff_poly_print(f,3);
			return 0;
		}

		k = i&1;
		while ( ! ecurve_random_point(&x,&y,g[k]) );
		ecurve_AJC_exp_ui(P,x,y,m[k],g[k][1]);
		if ( ! ecurve_fastorder2 (d[k], P,e[k],1,g[k][1]) ) {
			if ( swapped ) return 0;
			if ( (N[0]%m[1]) || (N[1]%m[0]) ) return 0;
			swap = g[0]; g[0] = g[1]; g[1] = swap;
			swapped = 1;
			if ( m[1] == 1 && m[0] == 1 ) continue;
			ppf_copy(tmp,n[1]);  ppf_divexact(tmp,e[1]); ppf_copy(tmp2,n[0]); ppf_divexact(tmp2,tmp);
			ppf_copy(tmp,n[0]);  ppf_divexact(tmp,e[0]); ppf_copy(e[1],n[1]); ppf_divexact(e[1],tmp);
			ppf_copy(e[0],tmp2);
			m3 = m[0]; m[0] = m[1]; m[1] = m3;
			continue;
		}
		ppf_divexact (e[k], d[k]);
		m[k] *= ppf_eval(d[k]);

//		M = (m[0]/ui_gcd(m[0],m[1]))*m[1];							// M can easily overflow here, so don't bother with this optimization
//		if ( M <= 2*r/3 ) continue;

		offset = (2*_ff_p+2)%m[1];									// orders permitted by twist lie in an arithmetic squence of the form a*m[1]+offset, with offset < m[1]
		if ( ! i_aseq_intersection (&m3,&k3,m[0],0,m[1],offset) || k3 < 0 )  { printf ("Error! Intersection of possible orders is empty!\n"); abort(); }
		if (  m3 < 0 || m3 > max || m3+k3 > max ) {
			o = k3;  m3 = -1;
		} else {
			o = m3*_ui_ceil_ratio(min,m3)+k3; if ( o-m3 >= min ) o-= m3;	// o is the first multiple of M >= min.  Note that min <= m3*ceil(min/m3) + k3 < min + m3 + k3 <= min + max <= 2(p+1)
		}
		// if there is only one possible order in the Hasse interval, we're done (one way or the other)
		if ( m3 < 0 || o+m3 > max ) return ( (o==N[0] || o == N[1]) ? 1 : 0 );
		// if there are two possible orders and both lie in {N0,N1}, we're done
		if ( o+2*m3 > max && ( (o == N[0] && o+m3 == N[1]) || (o==N[1] && o+m3 == N[0]) ) ) return 1;
	}
}


/*
	Slow Monte Carlo algorithm to test the group exponent, use for debugging/testing only.
*/
int ecurve_test_exponent (long e, ff_t f[4])
{
	ecp_jc_t t[1];
	ff_t x,y;
	long m, n;
	int i;
	
	n = 1;
	for ( i = 0 ; i < 100 ; i++ ) {
		if ( ! ecurve_random_point(&x,&y,f) ) continue;
		ecurve_AJC_exp_ui(t,x,y,e,f[1]);
		if ( ! ecurve_JC_id(t) ) return 0;
		m = ecurve_fastorder (x,y,e,f[1]);
		if ( n%m ) n = m*(n/ui_gcd(m,n));
	}
	return (n==e);
}

/*
	Given an elliptic curve with p-Sylow subgroup of order q, computes a basis (b1,b2) for the p-Sylow subgroup, with |b1|=q1 dividing |b2|=q2.
*/
void ecurve_p_basis (ecp_jc_t *b1, long *q1, ecp_jc_t *b2, long *q2, long p, long N, ff_t f[4])
{
	ff_t x,y;
	long q, r, e;
	int i;
	
	// set (b1,b2) to the trivial subgroup of the p-Sylow and attempt to expand this subgroup until it has order q
	// this "expansion" isn't necessarily monotonic, but the value of q2 is non-decreasing
	// the expected number of iterations is bounded by 2+2/p (we get an element of maximal order with prob 1-1/p on the first try and then get a basis with prob 1-1/p on the second try)

	for ( e = N/p, q = p ; e*q == N ; e /= p, q *= p );
	q /= p;  e = N / q;
	
	ecurve_JC_set_id (b1); *q1=1; ecurve_JC_set_id (b2); *q2 = 1;
	for ( i = 0 ; (*q1)*(*q2) < q && i < 100 ; i++ ) {							// bound the number of iterations, just in case the group order is wrong
		if ( ! ecurve_random_point(&x,&y, f) ) continue;						// get a random non-trivial point
		ecurve_AJC_exp_ui (b1, x, y, e, f[1]);								// get an element of the p-Sylow via the e-power map (overwrite b1, this is ok)
		r = ecurve_JC_pp_order (b1,p,q,f[1]);								// compute its order
		if ( r <= *q1 ) continue;											// if the order of b is not strictly larger than q1, it isn't going to help us so don't waste time on it
		if ( *q2 == 1 ) { *b2=*b1; *q2 = r; continue; } else *q1=r;				// first time through we just set q2, otherwise replace q1
		ecurve_p_reduce (b1, q1, b2, q2, p, f[1]);							// reduce to a basis for <b1,b2> -- note that this could actually reduce the size of our subgroup, but it cannot reduce q2
	}
	if ( i == 100 ) { printf ("p-Sylow structure computation failed on %ld-Sylow of (alleged) size %ld on curve over F_%ld: ", p, q, _ff_p); ff_poly_print(f,3); abort(); }
//printf ("Basis (%ld,%ld,%ld,%ld) of order %ld and (%ld,%ld,%ld,%ld) of order %ld\n", _ff_get_ui(b1->x), _ff_get_ui(b1->y), _ff_get_ui(b1->z2), _ff_get_ui(b1->z3),  *q1, _ff_get_ui(b2->x), _ff_get_ui(b2->y), _ff_get_ui(b2->z2), _ff_get_ui(b2->z3), *q2);
	return;
}


/*
	Given non-trivial elements b1 and b2 of the p-Sylow subgroup with |b1|=q1 and |b2|=q2  powers of p (p here is a divisor of the group order, not the characteristic of the field),
	The function below computes a basis for <b1,b2> and updates b1, b2, q1, and q2 appropriately, with q1 dividing q2 (possibly q1=1 if <b1,b2> is cyclic).

	The complexity is O(log_p(q1)*(sqrt(p)+log(q2)).  This algorithm is based on "Structure computation and discrete logarithms in finite abelian p-groups", Math Comp 80 (2011), pp.477-500.
*/
void ecurve_p_reduce (ecp_jc_t *b1, long *q1, ecp_jc_t *b2, long *q2, long p, ff_t f1)
{
	ecp_jc_t t[1],a1[1],a2[1];
	long k;
	
	if ( *q1 > *q2 ) { t[0] = *b1;  *b1 = *b2; *b2 = t[0];  k = *q1; *q1 = *q2; *q2 = k; }	// swap inputs if needed so that q1 <= q2

	ecurve_JC_exp_ui(a2,b2,*q2/p,f1);											// a2=b2^(q2/p), <a2> is the subgroup of <b2> of order p
	while (*q1>1) {
		ecurve_JC_exp_ui(a1,b1,(*q1)/p,f1);										// a1=b1^(q1/p), <a1> is the subgroup of <b1> of order p
		k = ecurve_dlog (a2,a1,p,f1);											// compute least k>0 s.t. a1=a2^k, if possible
		if ( ! k ) break;													// if we can't, then a1 and a2 are independent, hence so are b1 and b2, and we are done
		ecurve_JC_exp_ui(t,b2,(*q2)/(*q1)*k,f1);  ecurve_JC_invert(t);					// compute t=b2^-(q2/q1*k)
		ecurve_JCJC (b1,t,f1);												// replace b1 with b1*b2^-(q2/q1*k), this reduces the order of b1 to at most q1/p since
																		// (b1*b2^-(q2/q1*k))^(q1/p) = a1*a2^{-k} = 1, and we do not change <b1,b2> by doing this.
		k = ecurve_JC_pp_order(b1,p,*q1,f1);									// recompute q1 (it could be less than q1/p)
		if ( k >= *q1 ) { printf ("%lu: Error in ecurve_p_reduce(%ld,%ld) element order %ld not reduced\n", _ff_p, *q1, *q2, k); abort(); }		// sanity check
		*q1 = k;
	}
}


long ecurve_JC_pp_order (ecp_jc_t *a, long p, long q, ff_t f1)
{
	register long r;
	ecp_jc_t t[1];
	
	if ( ecurve_JC_id(a) ) return 1;
	if ( p==q ) return p;
	ecurve_JC_exp_ui (t,a,p,f1);
	if ( ecurve_JC_id(t) ) return p;
	for ( r = p*p ; r < q ; r *= p ) {
		ecurve_JC_exp_ui (t,t,p,f1);
		if ( ecurve_JC_id(t) ) return r;
	}
	return q;
}

/*
	Fast and simple small hash table lookup used by ecurve_bsgs_search and also by ecurve_dlog.
*/

#if FF_HALF_WORDS == 1
#define BSGS_MAX_STEPS		256					// guaranteed to handle  p<2^31 (single word ff_t), assuming at least 2 torsion is known, enlarge if needed but it's nice to fit in L1 cache
#else
#define BSGS_MAX_STEPS		4096				// will handle BSGS up to 2^44, assuming 2 and 3 torsion are known
#endif
#define BSGS_TABSIZE			BSGS_MAX_STEPS		// don't make this too big, it takes time to initialize it.  A few collisions won't kill us.
#define BSGS_TABMASK			((unsigned long)(BSGS_TABSIZE-1))

static ecp_jc_t babys[BSGS_MAX_STEPS];
static ecp_jc_t giants[BSGS_MAX_STEPS];

static ff_t stepzs[2*BSGS_MAX_STEPS];
#if FF_HALF_WORDS == 1
static unsigned char hashtab[BSGS_TABSIZE];
#else
static unsigned short hashtab[BSGS_TABSIZE];
#endif
static struct tab_entry {
	ff_t x;
	short i;
	short next;
} entries[BSGS_MAX_STEPS+1];
short nexttabentry;

static inline void tab_clear() { memset(hashtab,0,sizeof(hashtab)); nexttabentry = 1; } // don't use entry 0

// we require inserts to have unique x values, return -1 if unique, otherwise return index value for existing entry
static inline int tab_insert(ff_t x, short i)
{
	register struct tab_entry *p;
	register int n,h;

	h = x&BSGS_TABMASK;					// brutally simple hash function, probably not worth making this more sophisticated
	n = hashtab[h];
	if ( n ) {
		p=entries+n;
		for(;;) {
			if ( _ff_equal(p->x,x) ) return p->i;
			if ( ! p->next ) break;
			p=entries+p->next;
		}
	}
	p = entries+nexttabentry;
	_ff_set(p->x,x);
	p->i = i;
	p->next = n;
	hashtab[h] = nexttabentry++;
	return -1;
}

static inline int tab_lookup(ff_t x)
{
	register struct tab_entry *p;
	register int n;

	n = hashtab[x&BSGS_TABMASK];			// brutally simple hash function, probably not worth making this more sophisticated
	if ( ! n ) return -1;
	p=entries+n;
	for(;;) {
		if ( _ff_equal(p->x,x) ) return p->i;
		if ( ! p->next ) return -1;
		p=entries+p->next;
	}
	return -1;
}

/*
	Computes the discrete log of b wrt a given the order of a using a BSGS search.  Optimized for small searches (shares table lookup code with ecurve_bsgs_search).
	Returns 0 if b is not an element of <a> (this is not assumed).

	Note that we are only called for prime values whose square divides the group order, so we only use O(N^(1/4)) steps in the worst case, and even this is very rare
	since the prime must also divide p-1.
*/
int ecurve_dlog (ecp_jc_t *a, ecp_jc_t *b, long o, ff_t f1)
{
	ecp_jc_t c[1];
	long bsteps, gstep, gsteps;
	ff_t zinv[2], baby_x[1], baby_y[1], giant_x[1], giant_y[1];
	ff_t t0, t1;
	register int i, j;
		
	// hard-wire small o cases
	if ( ecurve_JC_id(a) ) return o;
	switch ( o ) {
	case 1:	return 0;
	case 2:	return ( ecurve_JC_cmp (a,b) ? 1 : 0 );
	case 3: 	i = ecurve_JC_cmp (a,b);  return ( i ? (i>0?1:2) : 0 );
	}

//printf ("dlog(%d) a=(%d,%d,%d,%d) b=(%d,%d,%d,%d)\n", o, _ff_get_ui(a->x), _ff_get_ui(a->y), _ff_get_ui(a->z2), _ff_get_ui(a->z3), _ff_get_ui(b->x), _ff_get_ui(b->y), _ff_get_ui(b->z2), _ff_get_ui(b->z3));
	
	// For small searches we will just search the whole range so we only match once (using just one inversion).   The only optimization we use is the usual one for inverses.
	bsteps = (long)sqrt(o/2);
	gstep = 2*bsteps+1;
	gsteps = _ui_ceil_ratio(o,gstep);
	if ( bsteps > BSGS_MAX_STEPS || gsteps > BSGS_MAX_STEPS ) { printf ("%ld: order %ld caused BSGS_MAX_STEPS to be exceeded in ecurve_dlog\n", _ff_p, o); abort(); }
//printf("bsteps=%d, gsteps=%d, giant step=%d\n", bsteps, gsteps, gstep);
	
	// Convert baby and giant step to affine coords for stepping
	ecurve_JC_exp_ui (c,a,gstep,f1);
	ecurve_JC_invert(c);
	_ff_set(zinv[0],a->z3);
	_ff_set(zinv[1],c->z3);
	ff_parallel_invert (zinv, zinv, 2);
	ecurve_JC_to_A (baby_x, baby_y, a, zinv[0]);
	ecurve_JC_to_A (giant_x, giant_y, c, zinv[1]);
//printf ("baby step(a) = (%d,%d)\n", _ff_get_ui(baby_x[0]), _ff_get_ui(baby_y[0]));
//printf ("giant step(a^-%d) = (%d,%d)\n", gstep, _ff_get_ui(giant_x[0]), _ff_get_ui(giant_y[0]));

	// Step
	babys[0] = *a;  giants[0] = *b;
	j = ecurve_AJC_steps (babys, baby_x[0], baby_y[0], bsteps, f1);							// baby steps
	if ( j < bsteps ) { printf ("%ld: baby step hit identity in dlog with o=%ld\n", _ff_p, o); abort(); }
	j = ecurve_AJC_steps (giants, giant_x[0], giant_y[0], gsteps, f1);							// giant steps
	if ( j < gsteps ) return j*gstep;

	// Convert to affine to get unique x coords
	for ( i = j = 0 ; i < bsteps ; i++, j++ ) _ff_set(stepzs[j], babys[i].z2);
	for ( i = 0 ; i < gsteps ; i++, j++ ) _ff_set(stepzs[j], giants[i].z2);
	ff_parallel_invert(stepzs, stepzs, j);																	// we assume j < FF_MAX_PARALLEL_INVERTS here
	for ( i = j = 0 ; i < bsteps ; i++,j++ ) ff_mult(babys[i].x, babys[i].x, stepzs[j]);								// only invert x-coord here, we never need to invert y
	for ( i = 0 ; i < gsteps ; i++,j++ ) ff_mult(giants[i].x, giants[i].x, stepzs[j]);

/*
printf ("affine baby x's:\n", j);
for ( i = 0 ; i < bsteps ; i++ ) printf ("   %d: %ld\n", i+1, _ff_get_ui(babys[i].x));
printf ("affine giant x's:\n", j);
for ( i = 0 ; i < gsteps ; i++ ) printf ("   %ld: %ld\n", i*gstep, _ff_get_ui(giants[i].x));
*/
	// Populate the table with babys.  Inserts should never fail (baby steps can't be inverses provided bsteps < o/2)
	tab_clear();
	for ( i = 0 ; i < bsteps ; i++ ) if ( (j=tab_insert(babys[i].x,i)) >= 0 ) break;
	if ( i < bsteps ) { printf("%ld: baby step insert failed in dlog\n", _ff_p); abort(); }

	// Now match giant steps
	for ( i = 0 ; i < gsteps ; i++ ) {
		if ( (j=tab_lookup(giants[i].x)) >= 0 ) {
			_ff_mult(t0,babys[j].y,giants[i].z3);	_ff_mult(t1,giants[i].y,babys[j].z3);								// normalize y values for comparison
			if ( _ff_equal(t0,t1) ) return i*gstep+j+1;
			return (i?i*gstep-j-1:o-j-1);
		}
	}
	return 0;
}


/*
	sets exp to the unique multiple of k in [low,high] if there is one (and returns 1) or sets exp=o and returns 0
*/
static inline int set_exp (long *exp, long o, long low, long high)
{
	register int i;
	i = ui_multiples_in_range(o,low,high);
	if ( ! i ) return -1;
	if ( i == 1 ) { *exp = (high/o)*o;  return 1; } else { *exp = o; return 0; }
}


/*
	Uses BSGS (and fastorder) to compute the order of the non-trivial element b, given that some multiple of |b| lies in [low,high] subject to (optional)
	modularity constraints specifed by a1, a2, and m.  If m=1 then a1 and a2 are ignored and no constraint is imposed.
	If m is greater than 1, then the group order is known to *not* be a multiple of m (if it is known to be a multiple of m, caller should exponentiate and shrink the interval).
	In this case a1 must be nonzero and if a2 is zero then it means we should search for an integer o in [low,high] that is congruent to a1 mod m and a multiple of |b|.
	If a2 is also nonzero, then o may be congruent to either a1 or a2 mod m.

	If the return value is 1, then *exp=o is the unique multiuple of k in [low,high].
        If the return value is 0, then *exp=|b|.  This usually means there is more than one multiple of k in [low,high], but not necessarily.
	If the return value is -1, then there is no multiple of |b| in [low,high] that satisfies the specified modularity constraint
	
	If pflag is set, then 1 is returned whenever any multiple *exp=o of |b| in [low,high] is found (without worrying about whether it is unique, satisfies constraints, or equal to the order of |b|)
	It is assumed that when pflag is used the caller will test whether exp is prime and if so uniquely determine the group order (for p > 29).
*/

int ecurve_bsgs_search (long *exp, ecp_jc_t *b, long low, long high, int m, int a1, int a2, int pflag, ff_t f1)
{
	ecp_jc_t p[64], bstep[1], gstep[1];
	ff_t zinv[4], b_x, b_y, bstep_x, bstep_y, gstep_x, gstep_y;
	register ff_t t0, t1;
	register int i,j,k,bsteps,gsteps,dsteps,usteps,tot_gsteps;
	long o, o1, o2, gbase, bspan, gspace, gap;
	
// if ( _ff_p > TEST_P ) printf ("%lu: ecurve_bsgs pt (%lu,%lu,%lu,%lu), low=%ld, high=%ld, m=%d, a1=%d, a2=%d, pflag=%d\n", _ff_p, _ff_get_ui(b->x), _ff_get_ui(b->y),_ff_get_ui(b->z2),_ff_get_ui(b->z3), low, high, m, a1, a2,pflag); 
	
	assert ( low <= high);
	if ( m==1 ) { a1 = a2 = 0; } else if ( m==2 ) { a1 = 1; a2 = 0; }							// just in case caller messed up
	if ( a2&& a1 > a2 ) { i = a1; a1 = a2; a2 = i; }										// make sure a1 is the smaller value
	bsteps = (long)sqrt(0.5*(a2?2.0:1.0)*(double)(high-low+1) / (double)m);					// compute the number of baby steps we need, note range is inclusive (and assumed to be non-empty)
	if ( ! bsteps ) bsteps = 1;															// always take at least one step, just to avoid some special cases
	if ( bsteps > BSGS_MAX_STEPS ) {  printf ("Exceeded BSGS_MAX_STEPS=%d! p=%ld, low=%ld, high=%ld, m=%d, a1=%d, a2=%d\n", BSGS_MAX_STEPS, _ff_p, low, high, m, a1, a2);  abort(); }
//	if ( a2 && (bsteps&1) ) bsteps++;													// make sure bsteps is even if we are covering two a-values
	bspan = m * bsteps;															// bspan is the distance between the identity and the last baby step
	gspace= 2*bspan;																// giant step spacing can be made as large as 2*bspan+1 (due to inverses), but we use 2*bspan because we want bspan to divide gspan...
	if ( m==1 ) gspace++;															// ...except when m=1, where we can actually use 2*bspan+1.

	/*
		Pick first giant step to be a multiple of as large a power of 2 as possible (and still lie in the interval [low,high]).
		We then need to adjust to make it congruent to a mod m, but this changes only a few of the low order bits
		(An optimal approach might search for the integer in [low,high] congruent to a mod m with minimal NAF hamming weight, but this would probably take more time than it is worth).
		
		The overall savings is only a few percent (it cuts the cost of computing the first giant step by about 1/6).
		However, it's an easy optimization, so there is no reason not to do it.
		
		This idea was suggested by Dan Bernstein.
	*/
	for ( k = ui_lg_floor(high-low) ;; k++) { o1 = (1L<<(k+1));  o = o1*_ui_ceil_ratio(low,o1);  if ( o > high ) break; }
	o1 = 1L<<k;
	gbase = o1*_ui_ceil_ratio(low,o1);													// gbase is now a multiple of a largish power of 2 in [low,high]
retry:
	if ( m > 1 ) { i = (gbase-a1)%m; gbase -= i; }										// make gbase congruent to a1 mod m
	while ( gbase < low ) gbase += m;													// stay in [low,high]
	if ( gbase > high ) gbase -= m;
	if ( gbase < low ) { if ( !a2 ) { return -1; } else { a1 = a2; a2 = 0; goto retry; } }			// if we can't find a giant step in [low,high], there is no point in continuing.
	for ( dsteps = (gbase-low)/gspace ; gbase - (dsteps-1)*gspace - bspan > low ; dsteps++ );
	while ( gbase-(dsteps-1)*gspace <= 0 ) gbase += m;									// make sure we don't step on zero or negative values
	for ( usteps = (high-gbase)/gspace ; gbase + (usteps-1)*gspace + bspan < high ; usteps++ );
	if ( dsteps+usteps+1 > BSGS_MAX_STEPS )  { printf ("Exceeded BSGS_MAX_STEPS=%d! p=%ld, low=%ld, high=%ld, m=%d, a1=%d, a2=%d\n", BSGS_MAX_STEPS, _ff_p, low, high, m, a1, a2);  abort(); }
	gsteps = dsteps+usteps-1;														// dsteps and usteps both include starting step, so total is one less then the sum, note that if a2>0 we actually take 2 sets of giant steps
	if ( a2 ) gap = a2-a1; else gap = 0;													// if a2 is set, we will take two sets of giant steps, offset by gap
	
// if ( _ff_p > TEST_P ) printf ("%lu: bsteps=%d, gsteps=%d(%d,%d), bspan=%ld, gstep=%ld, first giant = %ld, gap = %ld\n", _ff_p, bsteps, gsteps, dsteps, usteps, bspan, gspace, gbase, (a2?gap:0)); 
	
	// compute 2^k powers of b sufficient to compute b^gspace (note gstep is certainly bigger than gap) and b^gbase
	o = _ui_max(gspace,gbase);
	k = ui_lg_floor(o);
	p[0] = *b;
	ecurve_JC_powers (p, k+1, f1);													// compute binary powers of b for exponentiating, add 1 extra power for NAF
	ecurve_JC_exp_powers(bstep, p, m, f1);												// interval between baby steps
	ecurve_JC_exp_powers(gstep, p, gspace, f1);											// interval between giant steps
	
	// If either step is the identity we won't get very far, so we need to check for this (yes, it happens)
	if ( ecurve_JC_id(bstep) || ecurve_JC_id(gstep) ) {
		o = ( ecurve_JC_id(bstep) ? m : gspace );
		if ( pflag ) { *exp=o; return 1; }												// if pflag is set, just set exp to our known multiple of the order of b and return
		ff_invert(zinv[0],b->z3);  ecurve_JC_to_A (&b_x, &b_y, b, zinv[0]);						// Convert b to affine coords so we can do a fast order computation to get the exact order of b (which may properly divide m)
		o = ecurve_fastorder (b_x, b_y, o, f1);											// do the order computation
		return set_exp(exp,o,low,high);												// return either the order o or the unique multiple of o in [low,high] if there is one
	}
	
	// We need to convert the steps to affine coords so that we can take them more efficiently.  Convert them and the base point b to affine coords using one inversion.
	_ff_set(zinv[0],b->z3); _ff_set(zinv[1],bstep->z3);  _ff_set(zinv[2],gstep->z3);  ff_parallel_invert (zinv, zinv, 3);
	ecurve_JC_to_A(&b_x,&b_y, b, zinv[0]);  ecurve_JC_to_A(&bstep_x,&bstep_y, bstep, zinv[1]);  ecurve_JC_to_A(&gstep_x,&gstep_y, gstep, zinv[2]);

/*if ( _ff_p > TEST_P ) {
printf("affine base point (%ld,%ld)\n", _ff_get_ui(b_x),_ff_get_ui(b_y));
printf("affine baby step (%ld,%ld)\n", _ff_get_ui(bstep_x),_ff_get_ui(bstep_y));
printf("affine giant step (%ld,%ld)\n", _ff_get_ui(gstep_x),_ff_get_ui(gstep_y));
}*/
	
	// Take baby steps
	babys[0] = bstep[0];															// first step is equal to the interval size m
	j = ecurve_AJC_steps (babys, bstep_x, bstep_y, bsteps, f1); 								// take all the baby steps
	if ( j < bsteps ) { 																// if a baby step is the identity, we know a multiple of the element order, but need to get the exact value.
		o = (j+1)*m;	if ( pflag ) { *exp=o; return 1; }									// j-th baby step is (j*m)-th power of b (i.e. scalar multiple)
		o = ecurve_fastorder (b_x, b_y, o, f1);  return set_exp(exp,o,low,high);
	}	
	
	// Take first giant step, which will be somewhere in the middle of the giants array (depending on the dsteps/usteps split)
	ecurve_JC_exp_powers(giants+dsteps-1, p, gbase, f1);									// First step goes in the middle (and is counted as the first downward and the first upward step)
	j = ecurve_AJC_dsteps_1 (giants+dsteps-1, gstep_x, gstep_y, dsteps, f1);					// Take the downward steps first (note that we walk backwards here)
	if ( j < dsteps ) {																// if a giant step hit the identity, we're done
		o = gbase-j*gspace;  if ( pflag ) { *exp=o; return 1; }								// it might actually be faster to continue, but its simpler to catch this know and its rare in any case
		o = ecurve_fastorder (b_x, b_y, o, f1);  return set_exp(exp,o,low,high);
	}
	j = ecurve_AJC_steps (giants+dsteps-1, gstep_x, gstep_y, usteps, f1);						// Now take the upward steps
	if ( j < usteps ) {																// if a giant step hit the identity, we're done
		o = gbase+j*gspace;  if ( pflag ) { *exp=o; return 1; }
		o = ecurve_fastorder (b_x, b_y, o, f1);  return set_exp(exp,o,low,high);
	}
	gbase -= (dsteps-1)*gspace;														// adjust gbase to match giants[0]

	if ( a2 ) {																		// if we have two values of a mod m to check, take a second set of giant steps
		ecp_jc_t gapstep[1];
		ecurve_JC_exp_powers(gapstep, p, gap, f1);										// gap between a1 giant steps and a2 giant steps, reuse gstep variable for this
		if ( ecurve_JC_id(gapstep) ) {
			o = gap;  if ( pflag ) { *exp=o; return 1; }
			o = ecurve_fastorder (b_x, b_y, o, f1);  return set_exp(exp,o,low,high);
		}
		giants[gsteps] = giants[0];													// clone the (real) base giant step, which is congruent to a1 mod m
		ecurve_JCJC (giants+gsteps,gapstep,f1);											// shift it by gap=a2-a1 to make it congruent to a2 mod m
		j = ecurve_AJC_steps (giants+gsteps, gstep_x, gstep_y, gsteps, f1);					// Now take another complete set of steps congruent to a2 mod m all in one pass
		if ( j < gsteps ) {															// if a giant step hit the identity, we're done
			o = gbase+gap+j*gspace;  if ( pflag ) { *exp=o; return 1; }
			o = ecurve_fastorder (b_x, b_y, o, f1);  return set_exp(exp,o,low,high);
		}
	}
	tot_gsteps = (a2?2:1) * gsteps;
	
	// Convert to affine to get unique x coords
	for ( i = j = 0 ; i < bsteps ; i++, j++ ) _ff_set(stepzs[j], babys[i].z2);
	for ( i = 0 ; i < tot_gsteps ; i++, j++ ) _ff_set(stepzs[j], giants[i].z2);
	ff_parallel_invert(stepzs, stepzs, j);													// IMPORTANT: we assume j < FF_MAX_PARALLEL_INVERTS here, but we could batch if needed
	for ( i = j = 0 ; i < bsteps ; i++,j++ ) ff_mult(babys[i].x, babys[i].x, stepzs[j]);				// only convert x-coords, we only need to compare y's when we get a match
	for ( i = 0 ; i < tot_gsteps ; i++,j++ ) ff_mult(giants[i].x, giants[i].x, stepzs[j]);
	
/*if ( _ff_p > TEST_P ) {
printf ("affine baby x's:\n", j);
for ( i = 0 ; i < bsteps ; i++ ) printf ("   %d(%d): %ld\n", i, (i+1)*m, _ff_get_ui(babys[i].x));
printf ("affine giant x's:\n", j);
for ( i = 0 ; i < gsteps ; i++ ) printf ("   %ld: %ld\n", gbase+i*gspace, _ff_get_ui(giants[i].x));
if ( a2 ) for ( i = 0 ; i < gsteps ; i++ ) printf ("   %ld: %ld\n", gbase+gap+i*gspace, _ff_get_ui(giants[gsteps+i].x));
}*/

	// Populate the table with babys.  Insert will fail if we try to insert the same x value twice
	tab_clear(); for ( i = 0 ; i < bsteps ; i++ ) if ( (j=tab_insert(babys[i].x,i)) >= 0 ) break;

	// If we encountered two baby steps with the same affine x-coord, then we know a small multiple of the element order
	if ( i < bsteps ) {																// must must have 0 <= j < i < bsteps
		_ff_mult(t0,babys[j].y,babys[i].z3);	_ff_mult(t1,babys[i].y,babys[j].z3);					// normalize y values for comparison
		if ( _ff_equal(t0,t1) ) o = (j-i)*m; else o = (j+i+2)*m;								// gap is (j-i)*m if babys i and j are identical, o.w. (j+i+2)*m (note baby[0] is b^m, baby[1] is b^2m, etc...)
		if ( pflag ) { *exp = o; return 1; }
		o = ecurve_fastorder (b_x, b_y, o, f1);  return set_exp(exp,o,low,high);
	}
	
	// Now match giant steps by looking them up int the table of baby steps
	o1 = o2 = 0;
	for ( i = 0 ; i < tot_gsteps ; i++ ) {
		if ( (j=tab_lookup(giants[i].x)) >= 0 ) {
			_ff_mult(t0,babys[j].y,giants[i].z3);	_ff_mult(t1,giants[i].y,babys[j].z3);				// normalize y values for comparison
			if ( _ff_zero(t0) ) { o = 2*(j+1)*m;  return set_exp(exp,o,low,high); }				// handle 2-torsion case
			if ( _ff_equal(t0,t1) ) k = -(j+1)*m; else k = (j+1)*m;							// subtract baby index if giant=baby, add it if giant=baby^-1
			o = ( i < gsteps ? gbase+i*gspace+k : gbase+gap+(i-gsteps)*gspace+k );			// compute multiple of |b|, adjust for gap if match is in the second set of giant steps
//if ( _ff_p > TEST_P ) printf ("p=%lu, k=%d, gbase=%ld, gap=%ld, i=%d, gsteps=%d, o=%ld\n", _ff_p, k, gbase, gap, i, gsteps, o); 
			if ( o < low || o > high ) continue;											// this can happen in some edge cases, we ignore this because we are hoping for a unique multiple in [low,high]
			if ( o==o1 ) continue;													// ignore duplicate matches, this can happen!
			if ( o1 )  { o2 = o; break; }
			o1 = o;
		}
	}
	if ( ! o1 ) return -1;
	if ( ! o2 ) { *exp = o1; return 1; }													// we know the o1 is the unique multiple of the element order in [low,high]
	o=i_abs(o2-o1);																// otherwise, we know |o2-o1| is a multiple of the element order and there are at least two multiples in [low,high]
	if ( pflag ) { *exp=o; return 1; }													// we know there are two multiples in the interval, but when pflag is set we always return 1 or -1
	*exp = (m==1 ? o : ecurve_fastorder (b_x, b_y, o, f1) );									// if m is not 1 we don't know that o1 and o2 are adjacent multiples, but a fastorder computation will figure it out
	return 0;
}

/*
	Compute the order of affine point (x,y) given that (x,y)^e=1 using the classical algorithm.

	This algorithm is not particularly fast, but it doesn't need to be, as it is not invoked for most p.
	The average number of bits exponentiated per prime for a non-CM curve (fixed curve, varying p) is less than 1.
*/
long ecurve_fastorder (ff_t x, ff_t y, long e, ff_t f1)
{
	ecp_jc_t b[MAX_UI_PP_FACTORS];
	unsigned long q, p[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS];
	register long o;
	register int i, j, k, w;

	w = ui_factor(p,h,e);				// takes about 2 microseconds for e~2^32 (2.5MHz AMD-64), but e is generally much smaller
	if ( w == 1 && h[0] == 1 ) return e;	// not hard when the exponent is prime
		
	for ( i = 0 ; i < w ; i++ ) {
		for ( q=p[i],j=1 ; j < h[i] ; j++ ) q *= p[i]; 
		o = e/q;
		ecurve_AJC_exp_ui (b+i, x, y, o, f1);
	}
	o = 1;
	for ( i = 0 ; i < w ; i++ ) {
		for ( j = 0 ; j < h[i]-1 ; j++ ) {
			if ( ecurve_JC_id(b+i) ) break;
			ecurve_JC_exp_ui (b+i, b+i, p[i], f1);
		}
		if ( ! ecurve_JC_id (b+i) ) j++;
		for ( k = 0 ; k < j ; k++ ) o *= p[i];
	}
//printf ("fastorder computed order of (%ld,%ld) is %ld from exponent %ld\n", _ff_get_ui(x), _ff_get_ui(y), o, e);
	return o;
}

/*
	Computes n = |a| given that a^e = 1 (but see verify flag below). 

	Both n and e are factored integers.  The assumption is that the factorization is into prime powers,
	but it could be pseudo-primes.  In the latter case, n will be the smallest compatible divisor of e
	which is an exponent of |a| (compatible means that n doesn't split any of the pseudoprimes in e).

	The algorithm used here is optimized for small exponents which we expect to only have a handful
	of prime powers (almost all 64-bit integers will have log(log(2^64)) < 4 distinct prime factors). 
	We compute b1 = a^{e/q} and b2 = a^q, where q is the largest prime power p^h of e (we expect e/q and
	q to be roughly the same size) The results of this computation will be used to determine the
	exact power of p dividing |a| (usually q), and we then are left with the problem of computing |b2|
	using the exponent e/q which is roughly half the size of e.

	In the most likely cases, the total number of group operations is less than 2n where n is the number of bits in
	e, which is substantially better than the standard algorithm which seperately computes a^{e/q} for all the
	prime-powers q dividing e, probably using more than 4n group operations.

	For much larger values of e or, more importantly, e with many prime factors, one should instead use the
	recursive algorithm of Celler and Leedham-Green, or (better) the Snowball algorithm [2,Algorithm 7.4].

	It is assumed that the factored exponent e, in NAF, fits in a long.
	
	The boolean flag verify is used to indicate that is is not actually known that a^e=1, and the computation
	should either return 0 if it is found this is not the case, or verify its result before returning 1
	(if verify is 0, the return value is always 1).
*/

int ecurve_fastorder2 (ppf_t n, ecp_jc_t a[1], ppf_t e, int verify, ff_t f1)
{
	ecp_jc_t b1[1], b2[1], ap[64];
	ppf_t e2;
	register unsigned long q, m;
	register int i, j, k;

//printf ("Fastorder (%ld,%ld,%ld,%ld) e = ", _ff_get_ui(a->x), _ff_get_ui(a->y), _ff_get_ui(a->z2), _ff_get_ui(a->z3)); for ( i = 0 ; i < e->w ; i++ ) printf("%d^%d ", e->p[i], e->h[i]); puts("");
	
	// Catch trivial case quickly - bottoms out recursion
	if ( ecurve_JC_id (a) ) { n->w = 0;  return 1; }
	if ( e->w == 0 ) return 0;
	
	// if the order of a is a prime power p^h, things are simple
	if ( e->w == 1 ) {
		n->p[0] = e->p[0];  n->h[0] = 1;  n->w = 1;						// n = p
		b1[0] = a[0];
		for ( i = 1 ; i < e->h[0] ; i++ ) {
			ecurve_JC_exp_ui (b1, b1, e->p[0], f1);
			if ( ecurve_JC_id (b1) ) break;
			n->h[0]++;											// n*= p
		}
		if ( i == e->h[0] && verify ) {									// verify only if asked and if we didn't hit the identity
			ecurve_JC_exp_ui (b1, b1, e->p[0], f1);
			if ( ! ecurve_JC_id (b1) ) return 0;
		}
		return 1;
	}
	
	k = e->w-1;													// p = e->p[k] is the largest prime in e (note ppf factors are ordered by prime)
	for ( i = 0 ; i < k ; i++ ) { e2->p[i] = e->p[i]; e2->h[i] = e->h[i]; }		// copy first k-1 primes in e into e2
	for ( i = 1, q = e->p[k] ; i < e->h[k] ; i++, q *= e->p[k] );				// compute q = p^h
	e2->w = k;													// e2 is now e / q
	m = ppf_eval (e2);												// compute value of m = e/q
	i = ui_len(m);  j = ui_len(q); if ( j > i ) i = j;
	ap[0] = a[0];
	ecurve_JC_powers (ap, i, f1);
	ecurve_JC_exp_powers (b2, ap, q, f1);							// compute b2 = a^q first (if we only need one of b1=a^m or b2=a^q, its more likely to be b2)
	if ( ecurve_JC_id (b2) ) {										// don't need m
		if ( e->h[k] == 1 ) {										// if q is prime, we're done
			n->p[0] = e->p[k];  n->h[0] = 1;  n->w = 1;
			return 1;
		}
		e2->w = 1;  e2->p[0] = e->p[k];  e2->h[0] = e->h[k];				// e2 = q
		return ecurve_fastorder2 (n, a, e2, 0, f1);					// more convenient to recurse to finish up prime-power case, no need to verify
	}
	ecurve_JC_exp_powers (b1, ap, m, f1);							// compute b1 = a^m
	if ( ecurve_JC_id (b1) )										// dont need q
		return ecurve_fastorder2 (n, a, e2, 0, f1);					// e2 = m is our new (smaller) exponent of a (don't need to verify, we know a^m=1)

	//  We now know that both m and q contain a divisor of |a| (or e isn't an exponent at all)
	
	for ( j = 1 ; j < e->h[k] ; j++ ) {									// Are all the powers of p in q are needed?
		ecurve_JC_exp_ui (b1, b1, e->p[k], f1);
		if ( ecurve_JC_id (b1) ) break;
	}
	if ( j < e->h[k] ) {												// not all the powers of p are needed (rare)
		verify = 0;												// no need to verify if we hit the identity
		for ( i = 1, q = e->p[k] ; i < j ; i++, q *= e->p[k] );				// recompute q
		ecurve_JC_exp_powers (b2, ap, q, f1);						// recompute b2 = a^q
	}
	
	// Now we know the exact power of p in |a| is p^j
	if ( ! ecurve_fastorder2 (n, b2, e2, verify, f1) ) return 0;			// Recurse to compute |b2| which is now not divisible by p
	n->p[n->w] = e->p[k];   n->h[n->w] = j;  n->w++;					// set n = |a| = p^j|b2|
	return 1;
}

/*
	BSGS search hardwired for very short intervals (typically less than 50), which uses only one or two baby steps.
	These intervals can arise when the first order computation does not uniquely determine the exponent of the group (or when p is small).

	We assume b does not have 2 torsion (hence its order is at least 3) and high > low.
	We rely on the existence of an exponent of b in [low,high] and do not always verify this (e.g., if high=low+1 and b^low != id then we conclude b^high == id)
*/
int ecurve_bsgs_short (long *exp, ecp_jc_t *b, long low, long high, ff_t f1)
{
	long range;
	ecp_jc_t b2[1],b5[1],g[1];
	register ff_t t0, t1;
	register long e;
	long o1,o2;
	int i;

	range = high-low+1;
	switch (range) {
	case 2:		// E no computation other than the g exponentiation, cost E
		ecurve_JC_exp_ui(g, b, low, f1);
		*exp = ( ecurve_JC_id(g) ? low : high );
		return 1;
	case 3:		// E+4/3M (on average)
		ecurve_JC_exp_ui(g, b, low+1, f1);
		if ( ecurve_JC_id(g) ) { *exp = low+1;  return 1; }
		// it must be the case that g = +/- b, we don't verify this but simply compare y values to distinguish
		_ff_mult(t0,b->y,g->z3);  _ff_mult(t1,g->y,b->z3);
		*exp =  ( _ff_equal (t0, t1) ? low : high );
		return 1;
	case 4:		// E+1/4D+2M = E+5M+2A (roughly) (D indicates doubling(squaring) a point on the curve, which costs 11M+8A field operations)
		ecurve_JC_exp_ui(g, b, low+1, f1);
		if ( ecurve_JC_id(g) ) { *exp = low+1;  return 1; }
		_ff_mult(t0,b->x,g->z2);  _ff_mult(t1,g->x,b->z2);
		if ( _ff_equal (t0, t1) ) {
			_ff_mult(t0,b->y,g->z3);  _ff_mult(t1,g->y,b->z3);
			if ( ! _ff_equal (t0, t1) ) { *exp = low+2;  return 1; }
			// we know b^low = id, but it could also be that b^high = id if b has order 3, so we need to check this
			ecurve_2JC (b2,b, f1);
			_ff_mult(t0,b->x,b2->z2);  _ff_mult(t1,b2->x,b->z2);					// given b != id, b^3==id iff b and b^2 have the same x coord
			if ( _ff_equal(t0,t1) ) { *exp = 3; return 0; }
			*exp = low;
			return 1;
		}
		*exp = high;
		return 1;
	case 5:		// < E+4/5D+5M (on average)
		ecurve_JC_exp_ui(g, b, low+2, f1);
		if ( ecurve_JC_id(g) ) { *exp = low+2;  return 1; }
		ecurve_2JC (b2,b, f1);
		if ( ecurve_JC_id(b2) ) { *exp = 2; return 0; }
		if ( ecurve_JC_2tor(b2) ) { *exp = 4; return 0; }
		_ff_mult(t0,b->x,b2->z2);  _ff_mult(t1,b2->x,b->z2);
		if ( _ff_equal(t0,t1) ) { *exp = 3; return 0; }
		// we now know |b|>4, so there is a unique multiple in the interval
		_ff_mult(t0,b->x,g->z2);  _ff_mult(t1,g->x,b->z2);
		if ( _ff_equal (t0, t1) ) {
			_ff_mult(t0,b->y,g->z3);  _ff_mult(t1,g->y,b->z3);
			*exp = ( _ff_equal(t0,t1) ? low+1 : low+3 );
			return 1;
		}
		_ff_mult(t0,b2->y,g->z3);  _ff_mult(t1,g->y,b2->z3);
		*exp =  ( _ff_equal (t0, t1) ? low : low+4 );
		return 1;
	}
	ecurve_2JC (b2,b, f1);
	if ( ecurve_JC_id(b2) ) { *exp = 2; return 0; }
	if ( ecurve_JC_2tor (b2) ) { *exp = 4; return 0; }	
	_ff_mult(t0,b->x,b2->z2);  _ff_mult(t1,b2->x,b->z2);
	if ( _ff_equal(t0,t1) ) { *exp = 3; return 0; }
	ecurve_2JC (b5,b2, f1);
	ecurve_JCJC (b5,b, f1);
	if ( ecurve_JC_id(b5) ) { *exp = 5; return 0; }
	// we know |b|>5
	e = low+2;
	ecurve_JC_exp_ui(g, b, e, f1);
	o1 = o2 = 0;
	while ( e < high+3 ) {
		if ( ecurve_JC_id(g) ) { o2 = o1; o1 = e; goto next; }
		_ff_mult(t0,b->x,g->z2);  _ff_mult(t1,g->x,b->z2);
		if ( _ff_equal (t0, t1) ) {
			_ff_mult(t0,b->y,g->z3);  _ff_mult(t1,g->y,b->z3);
			i = ( _ff_equal(t0,t1) ? -1 : 1 );
			o2 = o1;  o1 = e+i; goto next;
		}
		_ff_mult(t0,b2->x,g->z2);  _ff_mult(t1,g->x,b2->z2);
		if ( _ff_equal (t0, t1) ) {
			_ff_mult(t0,b2->y,g->z3);  _ff_mult(t1,g->y,b2->z3);
			i = ( _ff_equal(t0,t1) ? -2 : 2 );
			o2 = o1;  o1 = e+i; goto next;
		}
next:	if ( o2 ) break;
		ecurve_JCJC (g,b5,f1);
		e += 5;
	}
	if ( ! o1 ) { printf ("%ld: No match found in short interval [%ld,%ld]\n", _ff_p, low, high);  abort(); }
	if ( o2 ) { *exp = o1-o2; return 0; }
	*exp = o1;
	// total cost is roughly E + 2D + (1+(range/5-1))*A + (range/5)*4M (D indicates point doubling, A indicates point addition)
	return 1;
}

/*
	The rest of this file contains code that is not currently used because it was found to to be suboptimal in testing (on an AMD Athlon) but
	may be useful on other platforms and/or other applications.
*/
#if 0

/*
	Below are elliptic curve arithmetic functions for Jacobian coordinates.  The reduced Chudnovsky coordinates
        are faster at everything except 2JC versus 2J is (11M+8A) vs (10M+10A), but this difference is negligible (1M~2.5A),
        and in testing it was found to be faster to use JC everywhere.
*/

// converts Jacobian to affine coords, given the inverse of z
static inline void ecurve_J_to_A (ff_t *px, ff_t *py, ecp_j_t *p, ff_t zinv)
{
	register ff_t t0,t1;
	
	_ff_square(t0,zinv);
	_ff_mult(*px,p->x,t0);
	_ff_mult(t1,t0,zinv);
	_ff_mult(*py,p->y,t1);
	// 4M
}

// converts Jacobian to reduced Chudnovsky Jacobian coords
static inline void ecurve_J_to_JC (ecp_jc_t *o, ecp_j_t *p)
{
	_ff_set(o->x,p->x);
	_ff_set(o->y,p->y);
	_ff_square(o->z2,p->z);
	_ff_mult(o->z3,o->z2,p->z);
	// 2M
}

// squares an affine point into Jacobian coords
static inline void ecurve_2AJ (ecp_j_t *p, ff_t x, ff_t y, ff_t f1)
{
	register ff_t t0,t1,t2,t3,a,b;

	_ff_add(p->z,y,y);  _ff_mult(t1,p->z,y); _ff_add(t2,t1,t1); _ff_mult(a,x,t2);		// z=2y, a=4xy^2, t1=2y^2, t2=4y^2
	_ff_mult(t3,t1,t2);													// t3=8y^4
	_ff_square(t1,x); _ff_add(b,t1,t1); _ff_addto(b,t1); _ff_addto(b,f1);				// b = 3x^2+f1*1^4
	_ff_square(t0,b); _ff_add(t2,a,a); _ff_sub(p->x,t0,t2);						// x = b^2-2a
	_ff_sub(t1,a,p->x); _ff_mult(t2,t1,b); _ff_sub(p->y,t2,t3);					// y = b(a-x)-8y^4 (note we use the new x here)
	// 6M+9A
}


// squares a point p1 in Jacobian coords (p3 is the output, may be equal to p1)
static inline void ecurve_2J (ecp_j_t *p3, ecp_j_t *p1, ff_t f1)
{
	register ff_t a, b, c, t0, t1, t2;
	
	_ff_square(t0,p1->x); _ff_add(t1,t0,t0); _ff_addto(t1,t0); 					// t1 = 3x^2
	_ff_square(t0,p1->z); _ff_square(t2,t0);  _ff_mult(t0,f1,t2); _ff_add(b,t1,t0);		// b = 3x^2+f1*z^4	(note that f1=a4 in 13.2.1.c  of HECHECC p.282)
	_ff_add(c,p1->y,p1->y); ff_mult(p3->z,p1->z,c);							// c=2y, z = 2yz
	_ff_mult(t2,c,p1->y); _ff_add(t0,t2,t2);  _ff_mult(a,t0,p1->x);					// a=4xy^2,t2=2y^2
	_ff_add(t0,a,a); _ff_square(t1,b); _ff_sub(p3->x,t1,t0);						// x = b^2-2a
	_ff_square(t0,t2); _ff_add(c,t0,t0);										// c = 8y^4
	_ff_sub(t0,a,p3->x); _ff_mult(t2,t0,b); _ff_sub(p3->y,t2,c);					// y = b(a-x)-c   -- note we use the new x value here
	// 10M+10A
}


// multiplies a point p in Jacobian coords by a (non-identity) affine point (p is an input and an output)
static inline void ecurve_AJ (ecp_j_t *p, ff_t x0, ff_t y0, ff_t f1)
{
	register ff_t a, c, e, f, t0, t1, t2;
	
	if ( _ff_zero(p->z) ) { _ff_set(p->x,x0); _ff_set(p->y,y0); _ff_set_one(p->z); return; }
	_ff_square(t0,p->z);  _ff_mult(a,x0,t0);									// a = x0z^2, b = x*1^2=x (since z0=1), and t0=z^2
	_ff_sub(e,p->x,a);													// e = a-b
	_ff_mult(t1,t0,p->z); _ff_mult(c,y0,t1);									// c = y0z^3, d=y*1^3=y
	_ff_sub(f,p->y,c);													// f = d-c
	if ( _ff_zero(e) && _ff_zero(f) ) { ecurve_2AJ (p, x0,y0,f1); return; }		// must use doubling code here, but at least its an affine double (6M+9A)
	_ff_square(t0,e); _ff_mult(t1,t0,e); _ff_mult(t2,t0,a);						// t1=e^3, t2=ae^2
	_ff_square(t0,f); _ff_sub(p->x,t0,t1); _ff_add(t0,t2,t2); _ff_subfrom(p->x,t0);		// x = f^2-e^3-2ae^2
	_ff_sub(t0,t2,p->x); _ff_mult(t2,f,t0); _ff_mult(t0,t1,c); _ff_sub(p->y,t2,t0);		// y = f(ae^2-x) - ce^3
	ff_mult(p->z,p->z,e);												// z = 1*z*e
	// 11M+6A
}


// multiplies p1 in Jacobian coords by p2 Jacobian coords  and puts the result in p1.  Handles identity and will square (double) if needed
// this is the slowest case and should be avoided when possible, AJ, AJC, and JCJC are all better
static inline void ecurve_JJ (ecp_j_t *p1, ecp_j_t *p2, ff_t f1)
{
	register ff_t a, b, c, d, e, f, t0, t1, t2;
	
	// we need to handle identity cases in general (although in some cases we might be able to rule this out)
	if ( _ff_zero(p2->z) ) return;
	if ( _ff_zero(p1->z) ) { *p1=*p2; return; }
	_ff_square(t0,p2->z);  _ff_mult(a,p1->x,t0);								// a = x1z2^2, and t0=z2^2
	_ff_square(t1,p1->z);  _ff_mult(b,p2->x,t1);								// b = x2z1^2, and t1=z1^2
	_ff_sub(e,b,a);														// e = a-b
	_ff_mult(t2,t0,p2->z); _ff_mult(c,p1->y,t2);								// c = y1z2^3
	_ff_mult(t2,t1,p1->z); _ff_mult(d,p2->y,t2);								// d = y2z1^3
	_ff_sub(f,d,c);														// f = d-c
	if ( _ff_zero(e) && _ff_zero(f) ) { ecurve_2J (p1,p1,f1); return; }			// must double if pts are equal (10M+10A), inverses will end up with z=0, so we let them through
	_ff_square(t0,e); _ff_mult(t1,t0,e); _ff_mult(t2,t0,a);						// t1=e^3, t2=ae^2
	_ff_square(t0,f); _ff_sub(p1->x,t0,t1); _ff_add(t0,t2,t2); _ff_subfrom(p1->x,t0);	// x = f^2-e^3-2ae^2
	_ff_sub(t0,t2,p1->x); _ff_mult(t2,f,t0); _ff_mult(t0,t1,c); _ff_sub(p1->y,t2,t0);	// y = f(ae^2-x) - ce^3	-- note we use the new x here
	ff_mult(t0,p1->z,p2->z);	_ff_mult(p1->z,t0,e);							// z = z1z2e
	// 16M+7A (ouch)
}


// Computes p=(x,y)^n where (x,y) !=1 is in affine coordinates and p is in Jacobian coordinates
// It takes advantage of the fact that all additions are of the form J+A, requiring only 11M, doubling is done in J+J form, using 10M
void ecurve_AJ_exp_ui (ecp_j_t *p, ff_t x0, ff_t y0, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	register ff_t negy0;
	int i;

	if ( n == 0 ) { _ff_set_zero(p->z); return; }
	if ( n == 1 ) { _ff_set(p->x,x0); _ff_set(p->y,y0); _ff_set_one(p->z); return; }
	ecurve_2AJ(p,x0,y0,f1);
	if ( n == 2 ) return;
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits)-2;						// we know the top two bits of the NAF are 10
hecurve_expbits+=i+1;
	_ff_neg(negy0,y0);
	m = (1UL<<i);
	for ( ; m ; m >>= 1 ) {
		ecurve_2J(p,p,f1);	 					// 10M+10A
		if ( m&pbits ) ecurve_AJ(p,x0,y0,f1);		// 11M+6A
		if ( m&nbits ) ecurve_AJ(p,x0,negy0,f1);		// 11M+6A
	}
}


// Computes p=a^e where a!=1 is in Jacobian coords, and so is p.
// This is slow and should be avoided.  Overlap is ok
void ecurve_J_exp_ui (ecp_j_t *p, ecp_j_t *a, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	ecp_j_t t, ai;
	int i;
	
	if ( n == 0 ) { _ff_set_zero(p->z); return; }
	if ( n == 1 ) { *p=*a; return; }
	ecurve_2J(&t,a,f1);							// 10M+10A
	if ( n == 2 ) { *p = t; return; }
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits)-2;
hecurve_expbits+=i+1;
	ai=*a;  ff_negate(ai.y);
	m = (1UL<<i);
	for ( ; m ; m >>= 1 ) {
		ecurve_2J(&t,&t,f1);						// 10M+10A
		if ( m&pbits ) ecurve_JJ(&t,a,f1);			// 16M+7A
		if ( m&nbits ) ecurve_JJ(&t,&ai,f1);			// 16M+7A
	}
	*p = t;
}


// Combines precomputed values p[i]=p[0]^(2^i) to compute p[0]^n.  Assumes all required powers are present: NAF potentially requires one more bit!
// CAUTION: if this is false, the resulting behavior is very strange and unpredictable.  You have been warned.
void ecurve_J_exp_powers(ecp_j_t *o, ecp_j_t *p, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	ecp_j_t t;
	int i;

	// handle small n quickly
	switch (n) {
	case 0:	_ff_set_zero(o->z); return; 
	case 1:	*o=p[0]; return;
	case 2:	*o=p[1]; return;
	case 3:	*o=p[0]; ecurve_JJ(o,p+1,f1); return;
	case 4:	*o=p[2]; return;
	case 5:	*o=p[2]; ecurve_JJ(o,p,f1); return;
	case 6:	*o=p[2]; ecurve_JJ(o,p+1,f1); return;
	case 7:	*o=p[0]; ff_negate(o->y); ecurve_JJ(o,p+3,f1); return;
	case 8:	*o=p[3]; return;
	}
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits);
	*o = p[i];
	i-=2;
hecurve_expbits+=i+1;
	m = (1UL<<i);
	for ( ; m ; m >>= 1, i-- ) {
		if ( (m&pbits) ) ecurve_JJ(o,p+i,f1);						 	// 16M+7A (these are expensive)
		if ( (m&nbits) )
			{ ff_negate(o->y); ecurve_JJ(o,p+i,f1); ff_negate(o->y); }	// 16M+9A (slightly more expensive)
	}
}

// Sets p[i] = p[0]^(2^i) for i in [0,k]. Input is an affine pt not equal to the identity, output is in Jacobian coords.
void ecurve_AJ_powers (ecp_j_t p[], ff_t x0, ff_t y0, int k, ff_t f1)
{
	ff_t a, b, c, x, y, z, t0, t1, t2;
	register ecp_j_t *e;
	
	_ff_set(p->x,x0);  _ff_set(p->y,y0);  _ff_set_one(p->z);
	ecurve_2AJ (p+1, x0, y0, f1);
	for ( e=p+k,p+=2 ; p <= e ; p++ ) ecurve_2J(p,p-1,f1); 			// 10M + 10A
}
#endif
