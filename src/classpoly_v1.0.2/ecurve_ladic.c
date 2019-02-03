#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "ff_poly.h"
#include "ff_poly/ffpolysmall.h"
#include "ecurve.h"

/*
    Copyright (c) 2010-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

static ff_t phi3[10];			// coefficients ab of X^aY^b terms are stored in order as 00,10,11,20,21,22,30,31,32,33, but we know 00 coeff is zero and 33 coeff is -1
static ff_t phi3_p;

void _phi3_reduce (void)
{
	register ff_t t1,t2;
	
	_ff_set_ui(t1,1000000000UL);  _ff_set_ui(t2,1855425871872UL);  _ff_mult(phi3[1],t1,t2);		// the 1,0 coeff 1855425871872000000000 is bigger than 64 bits, so split it into 2 factors
	_ff_set_i(phi3[2],-770845966336000000L);	// 1,1 coeff
	_ff_set_ui (phi3[3],452984832000000UL);	// 2,0 coeff
	_ff_set_ui(phi3[4],8900222976000UL); 		// 2,1 coeff
	_ff_set_ui(phi3[5],2587918086UL); 			// 2,2 coeff
	_ff_set_ui(phi3[6],36864000UL);			// 3,0 coeff
	_ff_set_i(phi3[7],-1069956);				// 3,1 coeff
	_ff_set_ui(phi3[8],2232UL);				// 3,2 coeff
	phi3_p = _ff_p;
}
static inline void phi3_reduce (void) { if ( phi3_p != _ff_p ) _phi3_reduce(); }

static inline void phi3_eval (ff_t f[5], ff_t J)
{
	register ff_t J2, J3, t1;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_add(t1,J,phi3[6]);
	_ff_sum_3_mults(f[0],J,J2,J3,t1,phi3[3],phi3[1]); 								 // note constant coeff is zero for phi3
	_ff_sum_3_mults(t1,J,J2,J3,phi3[7],phi3[4],phi3[2]); _ff_add(f[1],t1,phi3[1]);
	_ff_sum_3_mults(t1,J,J2,J3,phi3[8],phi3[5],phi3[4]); _ff_add(f[2],t1,phi3[3]);
	_ff_sum_2_mults(t1,J,J2,phi3[8],phi3[7]); _ff_addto(t1,phi3[6]); _ff_sub(f[3],t1,J3);		// X^3Y^3 coeff is -1
	_ff_set_one(f[4]);
	// 13M+13A (6 redc)
}

// computes the 3-division poly for the elliptic curve y^2=x^3+f1x+f0, made monic by assuming p > 3.
// if pd is non-null then *pd is set to a square-root of -D/3, where D is the discriminant of the 3-div poly, which is helpful when factoring it
// f_3/3 = x^4 + 2f1x^2 +4f0x - f1^2/3
static inline void ecurve_div3 (ff_t div[5], ff_t *pd, ff_t f[4])
{
	ff_t t0, t1, t2, d;
	
	_ff_square(t2,f[1]); _ff_neg(t1,_ff_third);
	ff_mult(div[0],t1,t2);  _ff_add(div[1],f[0],f[0]); _ff_x2(div[1]);  _ff_add(div[2],f[1],f[1]);  _ff_set_zero(div[3]);  _ff_set_one(div[4]);
	if ( ! pd ) return;	
	_ff_square(t1,f[0]); _ff_set_ui(d,27); _ff_mult(t0,t1,d);							// t0 = 27f0^2
	_ff_square(t2,f[1]); _ff_mult(t1,t2,f[1]); _ff_x2(t1); _ff_x2(t1);						// t1 = 4f1^3, t2 = f1^2
	_ff_add(d, t0, t1);														// d = 4f1^3+27f0^2 = -discriminant of f
	_ff_add(t0,_ff_third, _ff_third); _ff_x2(t0); _ff_square(t1,t0);						// t0 = (4/3)^2 = 16/9
	ff_mult(*pd,d,t0);														// d^2 = 256/81*D_f^2 = -D/3 where D is discriminant of the 3-div poly
}

// computes the 4-division poly for the elliptic curve y^2=x^3+f1x+f0, made monic.
// f_4/2 = x^6 + 5f1x^4+20f0x^3 - 5f1^2x^2-4f0f1x-8f0^2-f1^3
static inline void ecurve_div4 (ff_t div[7], ff_t f[4])
{
	ff_t  t0, t1, t2, a2, a, b, c5;
	
	_ff_set(a,f[1]); _ff_set(b,f[0]); _ff_set_one(div[6]); _ff_set_zero(div[5]);
	_ff_set_ui(c5,5); _ff_mult(div[4],c5,a);
	_ff_mult(t0,c5,b); _ff_x2(t0); _ff_add(div[3],t0,t0);
	_ff_square(a2,a); _ff_mult(t1,c5,a2); _ff_neg(div[2],t1);
	_ff_mult(t1,a,b); _ff_x2(t1); _ff_x2(t1); _ff_neg(div[1],t1);
	_ff_square(t1,b); _ff_x2(t1); _ff_x2(t1); _ff_x2(t1); _ff_mult(t2,a2,a); _ff_addto(t1,t2); _ff_neg(div[0],t1);
}
	
/*
	ecurve_3tor computes the size of the 3-torsion subgroup of y^2=f(x) = x^3+ax+b.
	If the 3-torsion is trivial, it will instead attempt to determine the 3-torsion subgroup of the twist, and if it is non-trivial
	will return a negative value whose absolute value is the size of the 3-torsion subgroup in the twist.

        This means that a return value of 1 indicates that neither the curve nor its twist have non-trivial 3-torsion.
	In this situation, if p=1mod3, then a_p=0mod3 and the group order must be congruent to 2 mod 3.

	If the curve is singular (discriminant of f is zero), 0 is returned.
*/
int ecurve_3tor (ff_t f[4])
{
	ff_t div[5],r[4],g[4],t0,t1,t2,d;
	int k;
	
#if FF_WORDS > 1
	err_printf ("ecurve_3tor only supports single word finite fields\n"); exit (0);
#endif

	if ( _ff_p == 3 ) { if ( _ff_zero(f[1]) ) return 0; else return 1; }					// no curve y^2=x^3+ax+b over F_3 has 3 torsion
	ecurve_div3 (div,&d,f);

	// handle p=2mod3 ourselves since we can easily save some time.  note that in this case, either both the curve and its twist have 3-torsion subgroup Z/3Z, or both have trivial 3-torsion
	// we also know that the 3-division poly must have an odd number of factors, so it splits 1,1,2 or 4
	if ( ! _ff_p1mod3 ) {
		if ( _ff_zero(div[1]) ) { k = _ff_poly_roots_d4(r,div,&d,0); return ( k ? 3 : 1); }						// the biquadratic case is handled quickly by ff_poly_roots_d4, so let it
		ff_depressed_cubic_resolvent(&t0,g,div);
		k = _ff_poly_roots_d3(r,g,&d,FF_POLY_EXACTLY_ONE_ROOT);
		if ( k != 1 ) { printf ("p=%ld is 2mod3 but cubic resolvent of 3div poly has k=%d!=1 roots!\n  f=", _ff_p, k); ff_poly_print(f,3); exit(0); }
		_ff_sub(t1,t0,r[0]);
		return ( ff_residue(t1) ? 3 : 1);
	}
	
	// Now handle p=1mod3 by computing the roots of the 3-division polynomial (we actually only need to know the number of roots and the value of one of them)
	k = _ff_poly_roots_d4(r,div,&d,0);
	if ( ! k ) return 1;
	_ff_square(t1,r[0]);  _ff_addto(t1,f[1]); _ff_mult(t2,t1,r[0]); _ff_addto(t2,f[0]);		// t2 = f(r[0])
//	if ( _ff_zero(t2) ) { printf ("impossible - non-id point (%lu,0) with 2-tor and 3-tor!\n", _ff_get_ui(r[j]));  ff_poly_print(f,3);  ff_poly_print(div,4); exit(0); }
	
	// We now rely on the fact that the 3-division poly of the curve and its twist are intimately related.  They have the same factorization pattern,
	// and either all the roots of the 3-division poly correspond to points on the curve, or none of them do, so we only need to check one.
	// In the latter case, it must be that all the roots of the 3-division poly of the twist correspond to points on the twist.
	if ( ff_residue(t2) ) return ( k==1?3:9); else return (k==1?-3:-9);
}

int ecurve_mod3 (ff_t f[4], int *m)
{
	ff_t div[5], g[5], a1[4], a2[4], b1[4], b2[4], a[4], b[4], r[4], d, t0;
	register ff_t t1,t2,w0,w1,w2,x;
	register int i, k;

	if ( m ) *m = 3;															// default modulus
	if ( _ff_p == 3 ) return 1;													// Every elliptic curve of the form y^2=x^3+Ax+B over F_3 has order 1 mod 3
	ecurve_div3 (div,&d,f);
	if ( _ff_p1mod3 ) {
		// division poly has an even number of roots, so it splits 1,1,1,1 or 1,3 or 2,2
		k = _ff_poly_roots_d4(r,div,&d,0);
		if ( ! k )  return 2;														// if there are no roots we must have a 2,2 split and the trace is 0 mod 3
		_ff_set(x,r[0]);
		_ff_square(t1,x);  _ff_addto(t1,f[1]); _ff_mult(t2,t1,x); _ff_addto(t2,f[0]);				// t2 = f(x)
		if ( k < 4 && m ) {
			// if k is not 4 then k=1 and we must be on the floor of the 3-volcano.  Note that t^2=4q mod 3 is not zero, so we can't have a supersingular curve.
			// Use our point of order 3 to get a 3-isogenous curve and then check if it is above the floor of the 3-volcano (which gives us the trace mod 9)
			// We apply Velu's formula, summing just one point.  s = 6x^2+2A, u= = 4f(x)+sx, A'=A-5s, B'=B-7u, isogenous curve is y^2=x^3+A'x+B'
			// We can then check whether the isogenous curve has a neighbor besides f in the 3-volcano (if so, it is above us and the volcano has height > 0.
			// But we don't want to mess with j=0 or j=1728.
			_ff_set_ui(d,6);
			_ff_square(t1,x); _ff_mult(w1,t1,d);  _ff_addto(w1,f[1]); _ff_addto(w1,f[1]); _ff_add(w0,d,_ff_negone); _ff_mult(w2,w0,w1); _ff_sub(g[1],f[1],w2);
			_ff_add(w2,t2,t2); _ff_x2(w2); _ff_mult(w0,w1,x); _ff_addto(w2,w0); _ff_inc(d); _ff_mult(w1,d,w2); _ff_sub(g[0],f[0],w1);  _ff_set_one(g[3]); _ff_set_zero(g[2]);
			ecurve_to_jinv2(a,b,f,g);												// a[0] is the j-inv of f, b[0] is the j-inv of the 3-isogenous curve g
			if ( ! _ff_zero(a[0]) ) {												// Don't mess with j-invariant 0, j=1728 is OK
				phi3_reduce();
				phi3_eval(g,b[0]);  ff_poly_remove_root_d4 (g,g,a);						// Compute Phi_3(X,b[0]) / (X-a[0])
				ff_depress_cubic(&d,g);
				if ( _ff_poly_roots_d3(0,g,0,0) ) k = 4;								// if g has a root then the isogenous curve is not on the floor and t^2=4q mod 9.  Set k=4, since the situation is effectively the same.
			}
		}
		// we now are in case (i) of Schoof's them and know t^2 = 4q = 1 mod 3
		if ( ff_residue(t2) ) {
			// We know the curve has a point of order 3 and t=-1 mod 3.  If k=4 then the curve order is divisible by 9
			if ( ! m ) return 0;
			if ( k==4 ) { *m = 9;  return 0; }
			return 0;
		} else {
			// We know the twist has a point of order 3 and t=1 mod 3.  If k=4 the the twist order is divisible by 9
			if ( m && k==4 ) *m = 9;
			if ( k < 4 ) return 1;
			switch ( (_ff_p%9) ) {
			case 1: return 4;	// t = 7 mod 9
			case 4: return 1;   // t = 4 mod 9
			case 7: return 7;	// t = 1 mod 9
			}
		}
	} else {
		// division poly has an odd number of roots, so it splits 1,1,2 or 4
		if ( _ff_zero(div[1]) ) {													// the biquadratic case is handled quickly by ff_poly_roots_d4, so let it
			k = _ff_poly_roots_d4(r,div,&d,0);
		} else {																// handle this case ourselves, since we can save a little time
			ff_depressed_cubic_resolvent(&t0,a,div);
			k = _ff_poly_roots_d3(r,a,&d,FF_POLY_EXACTLY_ONE_ROOT);
			if ( k != 1 ) { printf ("p=%ld is 2mod3 but cubic resolvent of 3div poly has k=%d!=1 roots!\n  f=", _ff_p, k); ff_poly_print(f,3); exit(0); }
			_ff_sub(t1,t0,r[0]);
			k = ( ff_residue(t1) ? 2 : 0 );
		}
		if ( k ) return 0;														// if the division poly has roots, they yield points on both the curve and its twist, needn't check that f(x) is a residue
		/*
			Now comes the fun part, we need to use Schoof's algorithm to figure out the sign of the trace (which we know is +/- 1 mod 3).  We also know that p=q is -1 mod 3.
			We have the equation (*)   (x^(q^2), y^(q^2)) - (x,y) = +/- (x^q,y^q) mod div, (y^2-f), and we know that the x-coords match (regardless of the sign), we just need to check the y-coords
		
			We compute a1 = x^q mod div,  y^q = y*b1 where b1 = f^((q-1)/2) mod div,  a2 = x^(q^2) = a1^q mod div, and y^(q^2) = y*b2 where b2 = b1^(q+1) mod div.
		
			To test the equation (*), we set a=a1-a2, b=b1-b2, and check whether b(a^2(2a1+a2)-b^2f) - (b1+1)a^3 is zero mod div, which
			corresponds to computing (x^q,y^q) + (x^(q^2),y^(q^2)) - (x,y)  with some reorganization.  If this is zero then the trace t is -1, otherwise it is +1.
		*/
		_ff_neg(g[0],div[0]); _ff_neg(g[1],div[1]); _ff_neg(g[2],div[2]);					// negate non-leading coeffs of div for small mod poly operations which assume this
		ff_poly_exp_mod_4 (b1,f,(_ff_p-1)/2,g);										// b1 = y^(q-1) = f^((q-1)/2) mod g
		ff_poly_square_mod_4 (b2,b1,g);
		ff_poly_mult_mod_4(b2,b2,f,g);											// temporarily set b2 = b1^2*f = f^q = f(x^q) mod div, so we can derive a1=x^q from it
		
		/*
			As described in Section 4.1 of Gaudry-Morain "Fast algorithms for computing the eiegenvalue in the Schoof-Elkies-Atkin algorithm",
			we will attempt to derive a1=x^q mod g from f(x^q) mod g by computing GCD(f(X)-f(x^q),g(X)) viewed as an element of R[X] where R = Fq[x]/g(x)
		*/
		_ff_add(t2,g[2],f[1]);  _ff_add(t1,g[1],f[0]);  _ff_square(w2,t2); _ff_neg(w0,g[0]);
		_ff_subfrom(b2[0],t1); ff_poly_square_mod_4(b,b2,g); _ff_addto(b2[0],t1);			// temporarily set b = (f(x^q)-t1)^2
		_ff_sum_2_mults(w1,f[1],w0,t2,w2); _ff_addto(b[0],w1);							// b =t2^2f1-t2t0 + (f(x^q)-t1)^2 now holds the denominator of a1
		_ff_add(w1,g[0],w2); _ff_neg(w0,w1);
		_ff_sum_3_mults(a[0],f[0],g[0],w0,b2[0],t1,w2);
		_ff_mult(a[1],w0,b2[1]); _ff_mult(a[2],w0,b2[2]); _ff_mult(a[3],w0,b2[3]);			// a = t2^2f0+t0t1+(t0-t2^2)f(x^q) now holds the numerator of a1
		if ( _ff_zero(b[3]) ) {														// if lc of the denominator is zero, revert to standard computation
			ff_poly_xn_mod_d4 (a1,_ff_p,g);										// a1 = x^q mod g
		} else {
			ff_poly_invert_mod_4(a2,&d,b,g);										// a2 = d/b
			_ff_invert(t0,d); ff_negate(t0);											// we could in theory avoid this, but it's a hassle
			_ff_mult(a1[0],t0,a2[0]); _ff_mult(a1[1],t0,a2[1]);
			_ff_mult(a1[2],t0,a2[2]); _ff_mult(a1[3],t0,a2[3]); 
			ff_poly_mult_mod_4(a1,a1,a,g);										// a1 = a/b
		}
		ff_poly_compose_mod_4 (a2,a1,a1,g);										// a2 = a1(a1(x)) = x^(q^2) mod g
		ff_poly_compose_mod_4 (b2,b1,a1,g);  ff_poly_mult_mod_4(b2,b2,b1,g);			// b2 = y^(q^2-1) = b1^(q+1) = b1*b1^q = b1*b1(x^q) = b1*b1(a1(x)) mod g
		for ( i = 0 ; i < 4 ; i++ ) {
			_ff_sub(a[i],a1[i],a2[i]); _ff_sub(b[i],b1[i],b2[i]);							// a = a1-a2, b=b1-b2
			_ff_x2(a1[i]); _ff_addto(a1[i],a2[i]);										// replace a1 by 2a1+a2 (temporarily)
		}
		_ff_inc(b1[0]);															// replace b1 by b1+1
		ff_poly_square_mod_4 (a2,a,g);											// replace a2 by a^2
		ff_poly_mult_mod_4 (a1,a1,a2,g);											// replace a1 by a^2(2a1+a2) (in terms of the original a1,a2 defined above)
		ff_poly_mult_mod_4 (a2,a2,a,g);											// replace a2 by a^3
		ff_poly_mult_mod_4 (b1,b1,a2,g);											// replace b1 by (b1+1)a^3 (in terms of the original b1 above)
		ff_poly_square_mod_4 (b2,b,g);											// replace b2 by b^2
		ff_poly_mult_mod_4 (b2,b2,f,g);											// replace b2 by b^2f
		for ( i = 0 ; i < 4 ; i++ ) { _ff_subfrom (a1[i],b2[i]); }							// replace a1 by a^2(2a1+a2)-b^2f (in terms of the originals)
		ff_poly_mult_mod_4 (a1,a1,b,g);											// replace a1 with b(a^2(2a1+a2)-b^2f) (in terms of the originals)
		for ( i = 0 ; i < 4 ; i++ ) if ( ! _ff_equal (a1[i],b1[i]) ) break;
		return ( i < 4 ? 2 : 1 );													//if a1=b1 then t=-1 mod 3 and order = p+1-t = 2+1+1 = 1 mod 3, ow t=+1 mod 3 and order is 2 mod 3
	}
	printf ("Unhandled case in ecurve_mod3 for p=%ld\n", _ff_p); ff_poly_print(f,3); exit (0);
}

/*
	Computes #E(Fp) mod 4, where E is a nonsingular elliptic curve y^2=f(x)=x^3+f1*x+f0.
*/
int ecurve_mod4 (ff_t f[4], int odd)
{
	ff_t r,D,g[3],h[3],m[2];
	register ff_t t0,t1,t2,c27;
	int k;
	
	if ( ! odd ) {
		k = _ff_poly_roots_d3 (&r,f,0,FF_POLY_EXACTLY_ONE_ROOT);
		if ( k==3 ) return 0;
		if ( k==1 ) {
			// Apply Prop. 1 of "Constructing elliptic curves with prescribed torsion" with n=1, in which case we have 4-torsion iff f'(r) is a residue
			_ff_square(t0,r); _ff_add(t1,t0,t0); _ff_addto(t1,t0); _ff_addto (t1,f[1]);		// compute f'(r) = 3r^2 + f1
			if ( ff_residue(t1) ) return 0; else return 2;
		}
	}
	// when f is irreducible, apply Nogami & Morikawa ICISC 2004 LCNS 3506 pp. 249-260
	// compute discriminant D=-4f1^3 - 27f0^2 (note this is equal to -108k^3-108k^2 in equation (26) of the paper)
	_ff_square(t0,f[1]);  _ff_mult(t1,t0,f[1]);  _ff_x2(t1); _ff_x2(t1);
	_ff_square(t2,f[0]);  _ff_set_ui(c27,27); _ff_mult(t0,t2,c27);  _ff_addto(t0,t1); _ff_neg(D,t0);
	if ( (_ff_p&3) == 1 ) {
		ff_exp_ui(&r,&D,(_ff_p-1)>>2);
		return ( _ff_one(r) ? 3 : 1 );
	} else {
		// this is the slow case, we need to compute (x^p - x)^3 + 3(x^p-x) + sqrt(D) modulo f  (as in equation (17) of the paper)
		_ff_neg(m[0],f[0]); _ff_neg(m[1],f[1]);			// note that small mod poly functions want the modulus signs negated
		ff_poly_xn_mod_d3 (g,_ff_p,m);  _ff_dec(g[1]);		// compute x^p-x mod f
		ff_poly_square_mod_3 (h,g,m);
		_ff_add(t0,f[1],f[1]); _ff_addto(t0,f[1]);  _ff_addto(h[0],t0);
		ff_poly_mult_mod_3 (h,h,g,m);					// h now holds (x^p-x)^3 + 3(x^p-x), which should be equal to a sqrt of D (this is the constant term of f~ negated)
		if ( ! _ff_zero(h[2]) || ! _ff_zero(h[1]) ) { if ( odd < 0 ) return 0; err_printf ("For p=%ld , (x^p-x)^3 + 3*f1*(x^p-x) = ", _ff_p); ff_poly_print(h,2); printf ("modulo "); ff_poly_print(f,3); puts("Should have degree zero!"); exit (0); }
		_ff_square(t0,h[0]);
		if ( ! _ff_equal(t0,D) ) { if ( odd < 0 ) return 0; err_printf ("With p=%ld, constant term of f~ is not a square-root of the discriminant of f ", _ff_p); ff_poly_print(f,3); exit (0); }
		return ( ff_residue(h[0]) ? 1: 3 );
	}
}

/*
	Given a point (x,y) on y^2=f(x)=x^3+f1x+f0, replaces (x,y) with a point (u,v) s.t. (u,v) composed with itself yields (x,y) or its inverse (we don't distinguish these cases!)
	Returns 1 for success, 0 if no such point exists.
*/
int ecurve_halve (ff_t *x, ff_t *y, ff_t f[4])
{
	ff_t g[5], r[4];
	register ff_t t0, t1, t2;
	int i, k;
	
	// construct g(z) = z^4 - 2(f1+3x^2)z^2 - 8y^2z + f1^2 - 3x(y^2+f1x + 3f0).  If (u,v)*(u,v)=(x,y) then g(u-x)=0 (note translation by x to kill degree 3 coeff of g).
	_ff_set_one(g[4]);  _ff_set_zero(g[3]);
	_ff_square(t0,*x);  _ff_add(t1,t0,t0); _ff_addto(t1,t0);  _ff_addto(t1,f[1]); _ff_x2(t1); _ff_neg(g[2],t1);			// g[2] = -2(f1+3x^2)
	_ff_add(t0,*y,*y); _ff_square(t1,t0); _ff_x2(t1); _ff_neg(g[1],t1);										// g[1] = -8y^2
	_ff_square(t0,*y); _ff_mult(t1,f[1],*x); _ff_addto(t0,t1); _ff_add(t1,f[0],f[0]); _ff_addto(t1,f[0]); _ff_addto(t0,t1);	// t0 = y^2 + f1x + 3f0
	_ff_add(t1,*x,*x); _ff_addto(t1,*x); _ff_mult(t2,t0,t1); _ff_square(t1,f[1]);  _ff_sub(g[0],t1,t2);				// g[0] = f1^2 - 3x(y^2+f1x+3f0)
//printf ("(%d,%d) havling poly: ", _ff_get_ui(*x), _ff_get_ui(*y)); ff_poly_print(g,4);
	
	k  = ff_poly_roots_d4 (r, g);
//printf ("%d roots\n", k);
	for ( i = 0 ; i < k ; i++ ) {
		_ff_add(t2,r[i],*x);
		_ff_square(t0,t2);  _ff_addto(t0,f[1]);  _ff_mult(t1,t0,t2); _ff_add(g[1],t1,f[0]);						// g[1] = f(r[i]+x)
		if ( ff_sqrt(g,g+1) ) break;
	}
	if ( i == k ) return 0;
//printf ("found (%d,%d)^2 = (%d,%d)\n", _ff_get_ui(t2), _ff_get_ui(g[0]), _ff_get_ui(*x), _ff_get_ui(*y));
	_ff_set(*x,t2); _ff_set(*y,*g);
	return 1;
}


/*
	Given the *unique* root x0 of f(x)=x^3+Ax+B and the x-coord x1 of a point P1 on y^2=f(x) (x1!=x0),
	verifies that x1 can be halved exactly k times (i.e. P1=2^k*P for some P, but P1=2^{k+1}Q does not hold for any Q).
	If min_flag is set, only verifies that x1 can be halved at least k times.

	This algorithm is based on Miret-Moreno-Rio-Valls "Determining the 2-Sylow subgroup of an elliptic curve over a finite field" in the
	case that the 2-Sylow is cyclic and non-trivial (equivalently, f has a unique root).  We work here with f(x)=x^3+Ax+b rather
        than x(x^2+alpha*x+beta), but use x0 to translate.

	Expected cost is 5/2*k+1 square roots, or 5/2*(k-1)+1 if min is set
*/
int ecurve_verify_2depth (ff_t x1, int k, int min_flag, ff_t x0, ff_t f[4])
{
	ff_t d, delta;
	register ff_t alpha, beta, fac, t0, b, t2;
	
//printf ("verify 2depth x1=%ld, k=%d, x0=%ld  f(x)=", _ff_get_ui(x1), k, _ff_get_ui(x0)); ff_poly_print(f,3);
	if ( !k && min_flag ) return 1;
	_ff_triple(alpha,x0);
	_ff_mult(beta,alpha,x0); _ff_addto(beta,f[1]);							// beta = 3x0^2+A
//printf ("alpha=%ld, beta=%ld\n", _ff_get_ui(alpha), _ff_get_ui(beta));
	_ff_add(fac,beta,beta); _ff_x2(fac);
	_ff_sub(t2,x1,x0);
	for ( ; k > 0 ; k-- ) {
		_ff_add(t0,t2,alpha); _ff_mult(delta,t0,t2); _ff_addto(delta,beta);		// delta = t2^2+alpha*t2+beta
//printf("t2=%ld, delta=%ld\n", _ff_get_ui(t2), _ff_get_ui(delta));
		if ( ! ff_sqrt(&delta,&delta) ) return 0;
		if ( min_flag && k==1 ) return 1;
//printf("sqrt(delta)=%ld\n", _ff_get_ui(delta));
		_ff_add(b,t2,delta); _ff_x2(b);
		_ff_square(t0,b); _ff_sub(d,t0,fac);
//printf("d=%ld\n", _ff_get_ui(d));
		if ( ! ff_sqrt(&d,&d) ) {
			_ff_sub(b,t2,delta); _ff_x2(b);
			_ff_square(t0,b); _ff_sub(d,t0,fac);
			if ( ! ff_sqrt(&d,&d) ) { printf ("unexpected failure in verify 2depth p=%ld!", _ff_p); ff_poly_print(f,3);  printf("x0=%ld\n", _ff_get_ui(x0));  printf ("k=%d\n", k); exit(0); }							// failure should be impossible here
		}
		_ff_add(t0,b,d);
		_ff_mult(t2,t0,_ff_half);
	}
//printf("final t2=%ld, %s\n", _ff_get_ui(t2), (ff_residue(t2)?"is a residue, failure":"is not a residue, success"));
	return  ( ff_residue(t2) ? 0 : 1 );
}

/*
	Given the x-coord x2 of a point P2 on f(x)=x^3+f1x+f0, computes the (possible) x-coord x1 of a point P1 which, when doubled, yields P2.
	Does not verify that f(x1) is a quadratic residue, but if P2 is not a 2-torsion point then it must be (o.w. P1 is on the twist but 2*P1 isn't)
*/
int ecurve_halve_x (ff_t x2[1], ff_t x1, ff_t f[4])
{
	ff_t g[5], r[4];
	register ff_t t0, t1, t2;
	int k;
	
	// construct g(z) = z^4 - 2(f1+3x^2)z^2 - 8y^2z + f1^2 - 3x(y^2+f1x + 3f0).  If (u,v)*(u,v)=(x,y) then g(u-x)=0 (note translation by x to kill degree 3 coeff of g).
	_ff_set_one(g[4]);  _ff_set_zero(g[3]);
	_ff_square(t0,x1);  _ff_add(t1,t0,t0); _ff_addto(t1,t0);  _ff_addto(t1,f[1]); _ff_x2(t1); _ff_neg(g[2],t1);			// g[2] = -2(f1+3x^2)
	_ff_addto(t0,f[1]); _ff_mult(t2,t0,x1); _ff_addto(t2,f[0]);												// t2=y^2=x^3+f1*x+f0
	_ff_add(t1,t2,t2); _ff_x2(t1); _ff_x2(t1); _ff_neg(g[1],t1);												// g[1] = -8y^2
	_ff_mult(t0,f[1],x1); _ff_addto(t2,t0); _ff_triple(t0,f[0]); _ff_addto(t2,t0);									// t0 = y^2+f1x+3f0
	_ff_triple(t1,x1); _ff_mult(t0,t1,t2); _ff_square(t1,f[1]); _ff_sub(g[0],t1,t0);								// g[0] = f1^2 - 3x(y^2+f1x+3f0)
	
	k  = ff_poly_roots_d4 (r, g);
	if ( ! k ) return 0;
	_ff_add(x2[0],x1,r[0]);
	return 1;
}



/*
	Returns the size of the 2-Sylow subgroup of the elliptic curve y^2=f(x)=x^3+f1x+f0, provided it is cyclic.
	Otherwise the return value is -1, which indicates that Z/2Z x Z/2Z is a subgroup (and the group order is divisible by 4)
        (we could compute the entire 2-Sylow in this case, but it would be much more time consuming)
*/
int ecurve_2Sylow (ff_t f[4])
{
	ff_t r[3];
	ff_t x, y;
	int n;
	
	n = ff_poly_roots_d3(r,f);
	if ( ! n ) return 1;
	if ( n > 1) return 0;
	
	_ff_set(x,r[0]);  _ff_set_zero(y);
	for ( n = 2 ; ecurve_halve (&x, &y, f) ; n<<= 1 );
	return n;
}


/*
	Computes the order and rank of the 4-torsion subgroup of the elliptic curve y^2=f(x)=x^3+f1x+f0. 
	Set o to the order and returns d to indicate the group Z/dZ x Z/eZ where d*e = o, d divides e (and may be 1)

	If flag8 is set, also check for Z/8Z (for a random curve with full Galois image in GL(2,Z/8Z), these occurs with probability 1/8).
	We could also check for Z/2ZxZ/8Z, or even Z/4ZxZ/8Z or Z/8ZxZ/8Z (these occur with prob. 3/64, 3/512 and 1/1536 resp);
	
	The latter two are rare enough not to be worth the trouble, but Z/2xZ/8Z might be worth doing.  Note, however, that Z/2xZ/4Z
	contains 4 points of order 4 (two pairs of inverses) and we would need to compute two of these.
*/
int ecurve_4tor (int *o, ff_t f[4], int flag8)
{
	ff_t r[3];
	ff_t u, v, v2;
	register ff_t t0,t1,t2;
	register int i, n;
	
	_ff_set_zero(t2);	// to avoid compiler warning
	n = ff_poly_roots_d3(r,f);
	for ( i = 0 ; i < n ; i++ ) {
		_ff_square(t0,r[i]); 	_ff_add(t1,t0,t0); _ff_addto(t1,t0); _ff_add(u,t1,f[1]);			// u = 3x^2+a where (x,0) is a 2-torsion point (x=r[i] was a root of f)
		if ( ff_sqrt(&v,&u) ) {													// v is a root of the translated halving poly
			_ff_add(t2,r[i],v);													// t2=x+v is a root of the halving poly
			_ff_square(t0,t2);  _ff_addto(t0,f[1]);  _ff_mult(t1,t0,t2); _ff_add(u,t1,f[0]);		// u = f(t2)
			if ( ff_sqrt(&v,&u) ) break;											// if successful, (t2,v) is a point of order 4
			_ff_sub(t2,r[i],v);													// t2=x-v is also a root of the halving poly
			_ff_square(t0,t2);  _ff_addto(t0,f[1]);  _ff_mult(t1,t0,t2); _ff_add(u,t1,f[0]);		// u = f(t2)
			if ( ff_sqrt(&v,&u) ) break;											// if successful, (t2,v) is a point of order 4
		}
	}
	if ( i == n ) { *o = n+1;  return (n==3?2:1); }									// no pts of order 4, so 4-torsion subgroup = 2-torsion subgroup
//printf ("%ld: found point(%ld,%ld) with order 4 (i=%d,n=%d) on curve y^2 = ", _ff_p, _ff_get_ui(t2), _ff_get_ui(v), i, n); ff_poly_print(f,3);
	if ( n==1 ) {																// rank 1 case, with a point of order 4
		if ( ! flag8 ) { *o = 4;  return 1; }
//printf("%ld: flag8 set (Z/4Z), attempting to halve order 4 point (%ld,%ld) on curve y^2 = ", _ff_p, _ff_get_ui(t2), _ff_get_ui(v)); ff_poly_print(f,3);
		_ff_set(u,t2);
		*o = ( ecurve_halve(&u,&v,f) ? 8 : 4 );										// check for a point of order 8 (this takes us outside the 4-torsion subgroup)
//if ( *o == 8 ) puts ("Succeeded"); else puts ("Failed");
		return 1;
	}
	if ( i > 0 ) { *o = 8; return 2; }	
		
	// We are here if the 2-rank is 2 and we succesfully halved the first point of order 2 that we tried.  We need to check one more to distinguish Z/2Z x Z/4Z from Z/4Z x Z/4Z.
	// Don't change v or t2 from above, we may need them below
	i = 1;
	_ff_square(t0,r[i]); 	_ff_add(t1,t0,t0); _ff_addto(t1,t0); _ff_add(u,t1,f[1]);				// u = 3x^2+a where (x,0) is a 2-torsion point (x=r[0] was a root of f)
	if ( ff_sqrt(&v2,&u) ) {														// v2 is a root of the translated halving poly
		_ff_add(t1,r[i],v2);														// t1is a root of the halving poly
		_ff_square(t0,t1);  _ff_addto(t0,f[1]);  ff_mult(t1,t1,t0); _ff_addto(t1,f[0]);			// t1 = f(t1)
		if ( ff_residue(t1) ) { *o = 16; return 4; }
	}
	*o = 8;
	return 2;
}

/*
	Fast test for 4-torsion based on Proposition 1 of "Constructing elliptic curves with prescribed torsion over finite fields" [AVS 2008].
	If non-null, the parameter tt points to the x-coord of a 2-torsion point (i.e. a root of f).
*/
int ecurve_4tor_test (ff_t *tt, ff_t f[4])
{
	ff_t r[3], g[3];
	register ff_t t0,t1;
	register int x0,n;
	
	if ( tt ) {
		_ff_square(t0,*tt); _ff_add(t1,t0,t0); _ff_addto(t1,t0); _ff_addto(t1,f[1]);
		x0 = ff_residue(t1);
		if ( (_ff_p&3)==3 ) {
			if ( x0 ) return 1;
			ff_poly_remove_root_d3(g,f,tt);
			if ( ! ff_poly_roots_d2(r,g,2) ) return 0;
			_ff_square(t0,r[0]); _ff_add(t1,t0,t0); _ff_addto(t1,t0); _ff_addto(t1,f[1]);
			return ff_residue(t1);
		} else {			
			ff_poly_remove_root_d3(g,f,tt);
			if ( ! ff_poly_roots_d2(r,g,2) ) return x0;
			if ( x0 ) {
				_ff_sub(t0,*tt,r[0]);
				return ff_residue(t0);
			} else {
				_ff_sub(t0,r[0],r[1]);
				return ff_residue(t0);
			}
		}
	}
	
	n = ff_poly_roots_d3(r,f);
	if ( ! n ) return 0;
	_ff_square(t0,r[0]); _ff_add(t1,t0,t0); _ff_addto(t1,t0); _ff_addto(t1,f[1]);
	x0 = ff_residue(t1);
	if ( n == 1 ) return x0;
	if ( (_ff_p&3)==1 ) {
		if ( x0 ) {
			_ff_sub(t0,r[0],r[1]);
			return ff_residue(t0);
		} else {
			_ff_sub(t0,r[1],r[2]);
			return ff_residue(t0);
		}
	} else {
		if ( x0 ) return 1;
		_ff_square(t0,r[1]); _ff_add(t1,t0,t0); _ff_addto(t1,t0); _ff_addto(t1,f[1]);
		return ff_residue(t1);
	}
}

#define PHI5_COEFFS		((5+1)*(5+2)/2)

char *Phi5_str[20] = {
"141359947154721358697753474691071362751004672000",	// 0,0 coeff at 0
"53274330803424425450420160273356509151232000",	    	// 1,0 coeff at 1
"-264073457076620596259715790247978782949376",		// 1,1 coeff at 2
"6692500042627997708487149415015068467200",			// 2,0 coeff at 3
"36554736583949629295706472332656640000",			// 2,1 coeff at 4
"5110941777552418083110765199360000",				// 2,2 coeff at 5
"280244777828439527804321565297868800",				// 3,0 coeff at 6
"-192457934618928299655108231168000",					// 3,1 coeff at 7
"26898488858380731577417728000",						// 3,2 coeff at 8
"-441206965512914835246100",							// 3,3 coeff at 9
"1284733132841424456253440",							// 4,0 coeff at 10
"128541798906828816384000",							// 4,1 coeff at 11
"383083609779811215375",								// 4,2 coeff at 12
"107878928185336800",								// 4,3 coeff at 13
"1665999364600",										// 4,4 coeff at 14
"1963211489280",										// 5,0 coeff at 15
"-246683410950",										// 5,1 coeff at 16
"2028551200",										// 5,2 coeff at 17
"-4550940",											// 5,3 coeff at 18
"3720",												// 5,4 coeff at 19 (and we know 5,5, coeff is -1)
};

// mixed partial derivative of Phi5
char *Phi5xy_str[15] = {
"-264073457076620596259715790247978782949376",	// 0,0 coeff at 0
"73109473167899258591412944665313280000",	    	// 1,0 coeff at 1
"20443767110209672332443060797440000",			// 1,1 coeff at 2
"-577373803856784898965324693504000",				// 2,0 coeff at 3
"161390933150284389464506368000",				// 2,1 coeff at 4
"-3970862689616233517214900",						// 2,2 coeff at 5
"514167195627315265536000",						// 3,0 coeff at 6
"3064668878238489723000",						// 3,1 coeff at 7
"1294547138224041600",							// 3,2 coeff at 8
"26655989833600",								// 3,3 coeff at 9
"-1233417054750",									// 4,0 coeff at 10
"20285512000",									// 4,1 coeff at 11
"-68264100",										// 4,2 coeff at 12
"74400",											// 4,3 coeff at 13
"-25",											// 4,4 coeff at 14
};

static int Phi5_init;
static mpz_t Phi5[20];
static ff_t phi5[20];
static mpz_t Phi5xy[20];
static ff_t phi5xy[20];
static unsigned long phi5_redp;

static void _phi5_reduce (void)
{
	register int i;
	
	if ( ! Phi5_init ) { for ( i = 0 ; i < 20 ; i++ ) mpz_init_set_str(Phi5[i], Phi5_str[i], 0);  for ( i = 0 ; i < 15 ; i++ ) mpz_init_set_str(Phi5xy[i], Phi5xy_str[i], 0);  Phi5_init = 1; }
	for ( i = 0 ; i < 20 ; i++ ) _ff_set_mpz(phi5[i], Phi5[i]);
	for ( i = 0 ; i < 15 ; i++ ) _ff_set_mpz(phi5xy[i], Phi5xy[i]);
	phi5_redp = _ff_p;
}

static inline void phi5_reduce (void) { if ( phi5_redp == _ff_p ) return;  _phi5_reduce(); }

// Evaluate Phi_5(J,Y), (d/dX Phi_5)(J,Y) and (d/dX^2 Phi_5)(J,Y)
static inline void phi5_eval (ff_t f[7], ff_t J)
{
	register ff_t J2, J3, J4, J5, t1;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_square(J4,J2); _ff_mult(J5,J4,J); _ff_add(t1,J,phi5[15]);
	_ff_sum_5_mults(t1,J,J2,J3,J4,J5,t1,phi5[10],phi5[6],phi5[3],phi5[1]);  _ff_add(f[0],t1,phi5[0]);
	_ff_sum_5_mults(t1,J,J2,J3,J4,J5,phi5[16],phi5[11],phi5[7],phi5[4],phi5[2]); _ff_add(f[1],t1,phi5[1]);
	_ff_sum_5_mults(t1,J,J2,J3,J4,J5,phi5[17],phi5[12],phi5[8],phi5[5],phi5[4]); _ff_add(f[2],t1,phi5[3]);
	_ff_sum_5_mults(t1,J,J2,J3,J4,J5,phi5[18],phi5[13],phi5[9],phi5[8],phi5[7]); _ff_add(f[3],t1,phi5[6]);
	_ff_sum_5_mults(t1,J,J2,J3,J4,J5,phi5[19],phi5[14],phi5[13],phi5[12],phi5[11]); _ff_add(f[4],t1,phi5[10]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5[19],phi5[18],phi5[17],phi5[16]); _ff_addto(t1,phi5[15]); _ff_sub(f[5],t1,J5);
	_ff_set_one(f[6]);
	// 33M+31A (10 redc)
}

// Evaluate  (d/dX Phi_5)(J,Y) and (d/dX^2 Phi_5)(J,Y)
static inline void phi5_eval_x_xx (ff_t fx[6], ff_t fxx[6], ff_t J)
{
	register ff_t J2, J3, J4, t0, t1, c5;

	_ff_square (J2,J); _ff_mult(J3,J2,J); _ff_square(J4,J2);		// we could avoid these 3 mults be resuing values from phi5_eval
	
	_ff_add(t1,J3,phi5[6]); _ff_addto(t1,J3);
	_ff_x2(J); _ff_add(t0,J2,J2); _ff_addto(J2,t0);  _ff_x2(J3); _ff_x2(J3);  _ff_set_ui(c5,5); ff_mult(J4,c5,J4);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5[15],phi5[10],t1,phi5[3]);  _ff_add(fx[0],t1,phi5[1]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5[16],phi5[11],phi5[7],phi5[4]); _ff_add(fx[1],t1,phi5[2]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5[17],phi5[12],phi5[8],phi5[5]); _ff_add(fx[2],t1,phi5[4]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5[18],phi5[13],phi5[9],phi5[8]); _ff_add(fx[3],t1,phi5[7]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5[19],phi5[14],phi5[13],phi5[12]); _ff_add(fx[4],t1,phi5[11]);
	_ff_sum_3_mults(t1,J,J2,J3,phi5[19],phi5[18],phi5[17]); _ff_addto(t1,phi5[16]); _ff_sub(fx[5],t1,J4);
	// 24M+31A (7 redc)	

	
	_ff_add(t0,J,J); _ff_addto(J,t0);  _ff_x2(J2); _ff_x2(J2); ff_mult(J3,c5,J3); _ff_inc(c5);
	_ff_sum_4_mults(t1,J,J2,J3,J4,c5,phi5[15],phi5[10],phi5[6]); _ff_addto(t1,phi5[3]); _ff_add(fxx[0],t1,phi5[3]);
	_ff_sum_3_mults(t1,J,J2,J3,phi5[16],phi5[11],phi5[7]);  _ff_addto(t1,phi5[4]); _ff_add(fxx[1],t1,phi5[4]);
	_ff_sum_3_mults(t1,J,J2,J3,phi5[17],phi5[12],phi5[8]);  _ff_addto(t1,phi5[5]); _ff_add(fxx[2],t1,phi5[5]);
	_ff_sum_3_mults(t1,J,J2,J3,phi5[18],phi5[13],phi5[9]);  _ff_addto(t1,phi5[8]); _ff_add(fxx[3],t1,phi5[8]);
	_ff_sum_3_mults(t1,J,J2,J3,phi5[19],phi5[14],phi5[13]);  _ff_addto(t1,phi5[12]); _ff_add(fxx[4],t1,phi5[12]);
	_ff_sum_2_mults(t1,J,J2,phi5[19],phi5[18]); _ff_addto(t1,phi5[17]); _ff_addto(t1,phi5[17]); _ff_sub(fxx[5],t1,J3);
	// 19M+34A (7 redc)
}

static inline void phi5xy_eval (ff_t f[], ff_t J)
{
	register ff_t J2, J3, J4, t1;

	_ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_square(J4,J2);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5xy[10],phi5xy[6],phi5xy[3],phi5xy[1]);  _ff_add(f[0],t1,phi5xy[0]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5xy[11],phi5xy[7],phi5xy[4],phi5xy[2]); _ff_add(f[1],t1,phi5xy[1]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5xy[12],phi5xy[8],phi5xy[5],phi5xy[4]); _ff_add(f[2],t1,phi5xy[3]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5xy[13],phi5xy[9],phi5xy[8],phi5xy[7]); _ff_add(f[3],t1,phi5xy[6]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5xy[14],phi5xy[13],phi5xy[12],phi5xy[11]); _ff_add(f[4],t1,phi5xy[10]);
	// 23M+20A (8 redc)
}

static inline void phi5x_eval (ff_t f[], ff_t J)
{
	register ff_t J2, J3, J4, t0, t1;
	
	_ff_square(J2,J);  _ff_mult(J3,J,J2);  _ff_square(J4,J2); _ff_add(t1,J3,phi5[6]); _ff_addto(t1,J3);
	_ff_x2(J); _ff_add(t0,J2,J2); _ff_addto(J2,t0);  _ff_x2(J3); _ff_x2(J3);  _ff_set_ui(t0,5); ff_mult(J4,t0,J4);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5[15],phi5[10],t1,phi5[3]);  _ff_add(f[0],t1,phi5[1]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5[16],phi5[11],phi5[7],phi5[4]); _ff_add(f[1],t1,phi5[2]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5[17],phi5[12],phi5[8],phi5[5]); _ff_add(f[2],t1,phi5[4]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5[18],phi5[13],phi5[9],phi5[8]); _ff_add(f[3],t1,phi5[7]);
	_ff_sum_4_mults(t1,J,J2,J3,J4,phi5[19],phi5[14],phi5[13],phi5[12]); _ff_add(f[4],t1,phi5[11]);
	_ff_sum_3_mults(t1,J,J2,J3,phi5[19],phi5[18],phi5[17]); _ff_addto(t1,phi5[16]); _ff_sub(f[5],t1,J4);
	// 28M+31A (10 redc)	
}

/*
	Determines the factorization of F(X)=Phi_5(X,j0) (assumes caller has evaluated Phi_5)
	If there are no roots, the return value is -r, where r=2,3,6 when f splits as 2,2,2, 3,3, or 6.
	If there is a root, the return value is r=1,2,4,5 when f splits 1,...,1, 1,1,2,2, 1,1,4, or 1,5, and root[0] is set to one of the roots.
	If a problem occurs (usually due to j=0, j=1728, supersingular curve, or a double root), 0 is returned.
	Even when successful, the value is not guaranteed to always be correct, the caller should
	verify the final result.
*/
int phi5_factor (ff_t roots[2], ff_t F[7])
{
	ff_t f[7], g[7], h[6], a[6], r[7], s;
	int res, k, d_h;
	
	// Depress F into a local poly f that we will destroy
	ff_poly_depress_monic (&s,f,F,6);																	// depress f for faster mod poly operations
	_ff_set_one(g[6]); _ff_set_zero(g[5]);
	_ff_neg(g[4],f[4]); _ff_neg(g[3],f[3]); _ff_neg(g[2],f[2]); _ff_neg(g[1],f[1]); _ff_neg(g[0],f[0]); 					// put f in form g=x^6 - f4x^4 - ... - f0 expected by mod poly functions
	ff_poly_xn_mod_d6(a,_ff_p,g);																		// a = x^p mod f
	_ff_set(h[0],a[0]); _ff_add(h[1],a[1],_ff_negone);
	_ff_set(h[2],a[2]); _ff_set(h[3],a[3]); _ff_set(h[4],a[4]); _ff_set(h[5],a[5]);									// h = x^p-x mod f
	d_h = ff_poly_degree(h,5);
	if ( d_h < 0 ) {																					// if h=0 then f is a product of 6 linear factors
		_ff_poly_distinct_roots (r, f, 6, FF_POLY_ONE_ROOT|FF_POLY_SPLIT);										// use _ff_poly_distinct_roots to get a root (it will probabilistically split them)
		_ff_sub(roots[0],r[0],s);
		return 1;
	}
	ff_poly_monic (h, &d_h, h, d_h);
	ff_poly_gcd_small(h,&d_h,f,6,h,d_h);
	if ( d_h > 2 ) return 0;																			// this can happen in rare cases, typically when j=0,1728
	k = ff_poly_roots_d2 (r, h, d_h);
	res = _ff_p%5;  res = ( (res==1 || res==4) ? 1 : 0 );													// res = 1 if p is a residue mod 5
	if ( k ) { _ff_sub(roots[0],r[0],s);  if ( k==1 ) return (res?5:0); _ff_sub(roots[1],r[1],s); return (res?2:4); }
	 if ( res ) return -3;																				// apply Schoof/Stickelberger to distinguish 3,3 case (even # factors iff (p|5)=1 )
	ff_poly_compose_mod_6 (h,a,a,g);																	// compute (x^p)^p mod f via composition (this works because we are in characteristic p)
	_ff_dec(h[1]);																					// h = x^(p^2)-x mod f
	d_h = ff_poly_degree(h,5);
	if ( d_h < 0 ) return -2; 																			// if f splits 2,2,2 then it divides x^(p^2)-x and the gcd must be trivial
	return -6;																					// otherwise f must be irreducible (but we don't actually verify this)
}

/*
Given 5-isogenous j-invariants j1 and j2 with j1 the j-invariant of E1: y^2=x^3+a1*x+b1,
the following magma code fragment computes a degree 2 factor of the 5-division poly of E1,
assuming there are no exceptional cases (e.g. partials vanish, j-invariants 0,1728, etc...).

	RXY<X,Y>:=PolynomialRing(F,2);
	Rx<x>:=PolynomialRing(F);
	Phi5:=RXY!ClassicalModularPolynomial(5);
	phix:=F!Evaluate(Evaluate(Derivative(Phi5,1),1,j1),2,j2);
	phiy:=F!Evaluate(Evaluate(Derivative(Phi5,2),1,j1),2,j2);
	phixx:=F!Evaluate(Evaluate(Derivative(Derivative(Phi5,1),1),1,j1),2,j2);
	phiyy:=F!Evaluate(Evaluate(Derivative(Derivative(Phi5,2),2),1,j1),2,j2);
	phixy:=F!Evaluate(Evaluate(Derivative(Derivative(Phi5,1),2),1,j1),2,j2);
	K1:=j1/(j1-1728);
	K2:=j2/(j2-1728);
	R1:=18*b1/a1;
	k1:=R1*j1;
	R2:=-k1*phix/(5*j2*phiy);
	k2:=R2*j2;
	A2:=R2^2*K2;
	a2:=-A2/48;
	z:=5*k2*(phiyy/phiy-2*phixy/phix)-k1*phixx/phix;
	p1:=5*(z/2 - R1*(K1/4+1/3) + 5*R2*(K2/4+1/3));
	f1:=-p1/2;
	f2:=p1^2/8 + (5^3*a2-a1)/12 +2/5*a1;
	F5:=x^2+f1*x+f2;
	print Factorization(F5), Factorization(DivisionPolynomial(EllipticCurve(x^3+a1*x+b1),5));
*/

void ff_setup_fifth(void);
extern ff_t _ff_fifth;																					// _ff_half, _ff_third and _ff_fourth are already defined in ff.h

#define ECURVE_MOD5_VERIFY	1																	// This seems to be unnecessary, but it doesn't cost much (1-2 percent)

/*
	ecurve_mod5 attempts to compute the order of the curve y^2=f(x) mod 5 (or 25)
	If 5 is an Elkies prime (and the curve is not supersingular and has j-inv != 0,1728) it will set a[0] to the exact order mod 5,
	and mod 25 if it is easy to determine that the 5-volcano has non-zero height.  The return value is 1 in this case.
	If 5 is an Atkin prime, either it will determine the order is p+1 mod 5 (i.e. t=0), or it will determine 2 possible values for
	the curve order, depending on the sign of t, which will be stored in a[0] and a[1].  The return value is 2 in this case.
	In rare situations (e.g. j-invariant 0 or 1728), no information is obtained and we return 0.
*/
int ecurve_mod5 (int a[2], int *m, ff_t f[4])
{
	ff_t phi1[7], Phix[6], Phiy[6], Phixx[6], Phiyy[5], Phixy[5], phi2[7], g1[2], g2[2], h[3], roots[2];
	ff_t j1, j2,phix,phiy,phixy,phixx,phiyy;
	register ff_t k1, k2, a1,a2, t0, t1, t2, c5, c, twelfth, z, p1, K1,K2,R1,R2,A2;									// way more registers than necessary, but we trust the compiler to deal with this
	int r, p, retries;

	*m = 5;																						// in general we give info mod 5, mod 25 if and when we can
	if ( _ff_zero(f[0]) || _ff_zero(f[1]) ) return 0;															// rule out 0 and 1728 right off the bat
	_ff_set_ui(c5,5); ff_setup_fifth();
	if ( ! ecurve_to_jinv(&j1,f) ) { err_printf ("Singular curve in ecurve_mod5\n");  ff_poly_print(f,3); abort(); }
	phi5_reduce();
	phi5_eval (phi1, j1);
	r = phi5_factor (roots, phi1);
	if ( ! r ) return 0;																				// if a problem occurs, return 0 which means we don't know anything
	p = _ff_p%5;
	if ( r==2 || r==-2 ) { a[0] = (p+1)%5; return 1; }														// By Prop 6.2 of Schoof 1995, if r=2 we must have t=0 (for either a 1,1,2,2 or 2,2,2 split).
	if ( r==-3 ) {
		if ( p == 1 ) { a[0] = 1; a[1] = 3; return 2; } else if ( p==4 ) { a[0] = 2; a[1] = 3; return 2; }				// We have a 3,3 split and t^2 = p (so p=1,4 mod 5)
		err_printf ("_ff_p=%lu=%lu mod 5 for j=%lu with r=%d in ecurve_mod5 is impossible!\n", _ff_p, _ff_p % 5, _ff_get_ui(j1), r); return 0;
	}
	if ( r==-6 ) {
		if ( p == 2 ) { a[0] = 2 ; a[1] = 4; return 2; } else if ( p==3 ) { a[0] = 1; a[1] = 2; return 2; }				// We have a 6 split (i.e. none) and t^2 = 3p  (so p=2,3 mod 5)
		err_printf ("_ff_p=%ld=%lu mod 5 for j=%lu with r=%d in ecurve_mod5 is impossible!\n", _ff_p, _ff_p % 5, _ff_get_ui(j1), r); return 0;
	}
	if ( r < 0 ) { err_printf ("_ff_p=%ld, j=%ld, r=%d unhandled in ecurve_mod5\n", _ff_p, _ff_get_ui(j1), r); return 0; }
	
	/*
		We have an Elkies prime and the 5-isogenous pair of j-invariants (j1,j2)
		We will compute a quadratic factor h of the 5-division polynomial.  See Schoof 1995 or pp. 125-129 of Blake-Seroussi-Smart
		The code below is optimized to do the minimum necessary to compute just the two coefficients of h that we need, which simplifies things
	
		Start by computing things that only depend on j1, in case we want to retry with a different j2 (in order to get the quadratic factor to split)
	*/
	phi5_eval_x_xx(Phix,Phixx,j1);																		// Phix = Phi_x(j1,Y), Phixx = Phi_xx(j1,Y)
	_ff_set(Phiy[0],phi1[1]); _ff_add(Phiy[1],phi1[2],phi1[2]); _ff_add(t1,phi1[3],phi1[3]);
	_ff_add(Phiy[2],t1,phi1[3]); _ff_add(t1,phi1[4],phi1[4]); _ff_add(Phiy[3],t1,t1);
	_ff_mult(Phiy[4],c5,phi1[5]); _ff_set_ui(Phiy[5],6);														// Phiy = d/dy Phi(j1,Y) = Phi_y(j1,Y)
	_ff_set(Phiyy[0],Phiy[1]); _ff_add(Phiyy[1],Phiy[2],Phiy[2]); _ff_add(t1,Phiy[3],Phiy[3]);
	_ff_add(Phiyy[2],t1,Phiy[3]); _ff_add(t1,Phiy[4],Phiy[4]); _ff_add(Phiyy[3],t1,t1);
	_ff_mult(Phiyy[4],c5,Phiy[5]); 																		// Phiyy = d/dy Phi_y(j1,Y) = Phi_yy(j1,Y)
	_ff_set(Phixy[0],Phix[1]); _ff_add(Phixy[1],Phix[2],Phix[2]); _ff_add(t1,Phix[3],Phix[3]);
	_ff_add(Phixy[2],t1,Phix[3]); _ff_add(t1,Phix[4],Phix[4]); _ff_add(Phixy[3],t1,t1);
	_ff_mult(Phixy[4],c5,Phix[5]); 																		// Phixy = d/dy Phi_x(j1,Y) = Phi_xy(j1,Y)
	_ff_set(j2,roots[0]);  retries=0;
retry:
	ff_poly_eval (&phix, Phix, 5, &j2);
	ff_poly_eval (&phiy, Phiy, 5, &j2);
	ff_poly_eval (&phixx, Phixx, 5, &j2);
	ff_poly_eval (&phiyy, Phiyy, 4, &j2);
	ff_poly_eval(&phixy,Phixy,4,&j2);
//printf ("phix=%ld, phiy=%ld, phixx=%ld, phiyy=%ld, phixy=%ld\n", _ff_get_ui(phix), _ff_get_ui(phiy), _ff_get_ui(phixx), _ff_get_ui(phiyy), _ff_get_ui(phixy));
	_ff_set(a1,f[1]); _ff_set_ui(c,1728); _ff_set(phi1[0],a1); _ff_set(phi1[1],j2); _ff_sub(phi1[2],j1,c); _ff_sub(phi1[3],j2,c);
	_ff_set(phi1[4],phix); _ff_set(phi1[5],phiy);
	if ( ! ff_parallel_invert_check (phi1,phi1,6) ) return 0;													// invert a1, j2, j1-1728, j2-1728, phix, phiy in parallel (store in phi1, which we no longer need)
	_ff_mult(K1,j1,phi1[2]); _ff_mult(K2,j2,phi1[3]);														// K1 = j1/(j1-1728), K2=j2/(j2-1728)
	_ff_set_ui(c,18); _ff_mult(t0,c,f[0]); _ff_mult(R1,t0,phi1[0]); _ff_mult(k1,R1,j1);								// R1 = 18*b1/a1 = -B1/A1 = -E6/E4, k1 = j1' = R1*j1
	_ff_mult(t0,phix,phi1[5]);
	_ff_neg(t1,t0);_ff_mult(t0,k1,t1); _ff_mult(t1,_ff_fifth,phi1[1]); _ff_mult(R2,t0,t1);								// R2 = -k1*phix/(5*j2*phiy) = -B2/A2 = -E6(q^l)/E4(q^l)
	_ff_mult(k2,R2,j2);																				// k2 = R2*j2
	_ff_set_one(t0); _ff_add(t1,t0,t0); _ff_neg(t0,t1); _ff_mult(t1,t0,phi1[4]);
	_ff_sum_2_mults(t0,phiyy,phixy,t1,phi1[5]); _ff_mult(t1,c5,k2); _ff_mult(t2,k1,phixx); ff_negate(t2);
	_ff_sum_2_mults(z,t1,t2,phi1[4],t0);																// z = ell*k2*(phiyy/phiy - 2*phixy/phix)-k1*phixx/phix
	_ff_square(t0,R2); _ff_mult(A2,t0,K2);																// A2 = E4(q^l) = R2^2*K2,
	_ff_square(t0,_ff_half); _ff_mult(twelfth,t0,_ff_third); _ff_mult(c,t0,twelfth); _ff_mult(t0,c,A2); _ff_neg(a2,t0);			// a2 = -A2/48
	_ff_mult(t0,K1,_ff_fourth); _ff_addto(t0,_ff_third); _ff_neg(t1,t0); _ff_mult(t0,K2,_ff_fourth); _ff_addto(t0,_ff_third);
	_ff_mult(t2,c5,t0); _ff_sum_3_mults(t0,z,R1,R2,t2,t1,_ff_half); _ff_mult(p1,c5,t0);								// p1:=5*(z/2 - R1*(K1/4+1/3) + 5*R2*(K2/4+1/3))
	_ff_mult(t0,p1,_ff_half); _ff_neg(h[1],t0);																// h[1] = -p1/2
	_ff_square(t1,t0); _ff_mult(t2,t1,_ff_half); _ff_square(t0,c5); _ff_mult(c,t0,c5); _ff_mult(t1,c,a2); _ff_subfrom(t1,a1);
	_ff_mult(t0,t1,twelfth); _ff_addto(t2,t0); _ff_mult(t0,_ff_fifth,f[1]); _ff_addto(t2,t0); _ff_add(h[0],t2,t0);				// h[0] = p1^2/8 + (5^3*a2-a1)/12 + 2/5*a1 (note that we have simplified this using ell=5)
	_ff_set_one(h[2]);
/*
		From the value of p%5 and r, we can determine exactly when the quadratic h splits.
		This will occur whenever the 5-div poly has a root, which occurs precisely when either E or its twist has a non-trivial 5-torsion point.
		This means that if p+1+/-t = 0 mod 5 (with either sign) then h must split and otherwise not.
		We know from Prop 6.2 of Schoof 1995 "Counting points on elliptic curves over finite fields" that (modulo 5):
		
			p=1,r=1,5 implies t=+/-2 and n=0,4, so h splits, but when p=1,r=2 we have t=0 and n=2, so h does not split
			p=2,r=4 implies t=+/-2 and n=0,1, hence h splits.  Otherwise p=2 is not an Elkies prime
			p=3,r=4 implies t=+/-1 and n=3,0, hence h splits.  Otherwise p=3 is not an Elkies prime
			p=4,r=2 implies t=0 and n=0, hence h splits, but when p=4, r=1,5 we have t=+/-1 and n=+/-1 so h does not split
			
		But note that the cases where r=2 have t=0 which was already handled above.
		
		We verify that h actually is a divisor of the 5-division poly, if requested, although in fact this *never* seems to fail (perhaps we can prove this?)
		Recall that f_5 = F^2*f_4 - f_3^3 where F=4f for the curve y^2=f(x)=x^3+ax+b.
		When h splits, we can just evaluate f_5 at a root of h, otherwise we compute f_5 mod h.  In both cases we use the formula above, rather than multiplying out f_5.
	*/

/*
ff_poly_print(f,3);
printf ("p=%ld = %d mod 5, r=%d, 5-division poly factor is: ", _ff_p, p, r);
ff_poly_print(h,2);
*/
	ecurve_div3(phi1,0,f);  ecurve_div4(phi2,f);														// grab the 3 and 4 division polys which we will use to evaluate/reduce the 5-division poly
//printf ("division polys: "); ff_poly_print(phi1,4); ff_poly_print(phi2,6);
	if ( r==4 || p==1 ) {
		if ( ! ff_poly_roots_d2(h,h,2) ) { if ( ! retries ) { _ff_set(j2,roots[1]);  retries++; goto retry; } else { printf ("retry failed for j1=%ld, p=%ld!\n", _ff_get_ui(j1), _ff_p); ff_poly_print(f,3); exit(0); } }
		_ff_set(phix,h[0]); _ff_square(t2,phix); _ff_addto(t2,f[1]); _ff_mult(t1,phix,t2); _ff_add(phiy,t1,f[0]);		// phiy = f(phix) where phix is a root of h
		if ( ECURVE_MOD5_VERIFY ) {
			ff_poly_eval(&phixx,phi2,6,&phix); _ff_x2(phixx);
			ff_poly_eval(&phiyy,phi1,4,&phix); _ff_add(t0,phiyy,phiyy); _ff_addto(phiyy,t0);
			_ff_x2(phiy); _ff_x2(phiy); _ff_square(t0,phiy); _ff_mult(t1,t0,phixx); _ff_square(t0,phiyy); _ff_mult(t2,t0,phiyy);
			if ( ! _ff_equal(t1,t2) ) { err_printf ("Verification failed for p=%lu, j1=%lu, j2=%lu, r=%d, x=%lu is not a root of the 5-division poly!\n", _ff_p, _ff_get_ui(j2), _ff_get_ui(j2), r, _ff_get_ui(phix)); ff_poly_print(f,3); return 0; }
		}
		if ( r==1 ) *m = 25;																			// we could check the height of the 5-volcano when r=5, but tests show it is not worth doing (at least for p < 2^40)
		// OK, we have a root phix of the 5-division poly, we just need to test whether it corresponds to a point on the curve or its twist
		if ( ff_residue(phiy) ) { a[0]=0; return 1; }		// point is on the curve, so we have a point of order 5 (or 25, when p+1-t=0 mod 5 and t^2=4p mod 25 we have p=1 and p+1-t=0 mod 25)
		// twist has a point of order 5, so w have p+1+t=0
		switch ( p ) {
		case 1: if ( *m==5 ) { a[0] = 4; return 1; }		// t=-2, p+1-t = 4 mod 5
			// Since r=1 we know that t=-2=3 mod 5 and t^2=4p mod 25
			switch (_ff_p%25 ) {
			case 1: a[0] =  4; return 1;				// t = 23, p+1-t = 4 mod 25
			case 6: a[0] =  14; return 1;				// t = 18, p+1-t = 14 mod 25
			case 11: a[0] =  24; return 1;			// t = 13, p+1-t = 24 mod 25
			case 16: a[0] =  9; return 1;				// t = 8, p+1-t = 9 mod 25
			case 21: a[0] =  19; return 1;			// t = 3, p+1-t = 19 mod 25
			default:err_printf ("Unhandled case p=%lu with r=%d in ecurve_mod5\n", _ff_p, r); exit(0);
			}
		case 2: a[0] = 1; return 1;					// t=2, p+1-t = 1 mod 5
		case 3: a[0] = 3; return 1;					// t=1, p+1-t = 3 mod 5
		default: err_printf ("Unhandled case p=%lu with r=%d in ecurve_mod5\n", _ff_p, r); exit(0);
		}
	}
	/*
		If we reach this point we must have p=4 and r=1 or 5 with t=+/-1.
		As a sanity check, we may wish to verify that the 5-division poly f_5 is actually zero mod h.
	*/
	if ( p != 4 ) { err_printf ("Unhandled case in ecurve_mod5\n"); exit (0); }
	ff_negate(h[1]); ff_negate(h[0]);																	// negate coefficients for mod poly operations
	ff_poly_mod_2(phi2,6,h);  																		// reduce division poly (psi_4/4y) mod h
	if ( ECURVE_MOD5_VERIFY ) {
		ff_poly_mod_2(phi1,4,h); _ff_triple(g1[0],phi1[0]);  _ff_triple(g1[1],phi1[1]); 							// phi1 = f_3 mod h
		ff_poly_square_mod_2 (g2,g1,h); 
		ff_poly_mult_mod_2(g1,g1,g2,h);																// g1= f_3^3 mod h
		_ff_set(phi1[0],f[0]); _ff_set(phi1[1],f[1]); _ff_set_zero(phi1[2]); _ff_set_one(phi1[3]);						// copy f into phi1 so we can reduce it inplace
		ff_poly_mod_2(phi1,3,h); _ff_x2(phi1[1]); _ff_x2(phi1[1]); _ff_x2(phi1[0]); _ff_x2(phi1[0]); 					// phi2 = F = 4f mod h
		ff_poly_square_mod_2(phi1,phi1,h);																// phi1 = F^2 mod h
		_ff_add(g2[0],phi2[0],phi2[0]); _ff_add(g2[1],phi2[1],phi2[1]);										// g2= 2*(psi_4/4y) = f_4 mod h
		ff_poly_mult_mod_2(g2,g2,phi1,h);																// g2 = F^2*f_4 mod h
		if ( ! _ff_equal(g1[0],g2[0]) || ! _ff_equal(g1[1],g2[1]) )
			{ err_printf ("Verification failed for p=%lu, j1=%lu, j2=%lu, r=%d, h is not a factor of the 5-division poly!\n", _ff_p, _ff_get_ui(j1), _ff_get_ui(j2), r); ff_poly_print(f,3); ff_poly_print(h,2); return 0; }
	}

	// Now we just need to determine whether t is +/- 1, equivalently, whether the eigen value lambda is -/+2
	_ff_set_one(phi1[3]); _ff_set_zero(phi1[2]); _ff_set(phi1[1],f[1]); _ff_set(phi1[0],f[0]);							// copy f so we reduce inplace
	ff_poly_mod_2(phi1,3,h);  ff_poly_exp_mod_2(phi1,phi1,(_ff_p+3)/2,h);	
	_ff_set_ui(t0,8); ff_mult(phi1[0],phi1[0],t0);  ff_mult(phi1[1],phi1[1],t0);									// phi1 = y^q * 8*f*y = 8*f^{(q+3)/2} mod h
	if ( _ff_equal(phi1[0],phi2[0]) && _ff_equal(phi1[1],phi2[1]) ) { a[0]= 1; return 1; }							// lambda = 2, t = -1, p+1-t = 1 mod 5
	_ff_add(t0,phi1[0],phi2[0]); _ff_add(t1,phi1[1],phi2[1]);
	if ( ! _ff_zero(t0) || ! _ff_zero(t1) ) { err_printf ("Neither sign of lambda works in ecurve_mod5 with p=%ld=4 mod 5 ", _ff_p); ff_poly_print(f,3); ff_poly_print(h,2); exit (0); }
	if ( r==1 ) *m = 25;																				// we could check the height of the 5-volcano when r=5, but tests show it is not worth doing (at least for p < 2^40)
	if ( *m == 5 ) { a[0] = 4; return 1; }																// lambda = -2, t = 1, p+1-t = 4 mod 5
	// Since r=1 we know that t=1 mod 5 and that t^2=4p mod 25
	switch ( _ff_p%25 ) {
	case 4: a[0] =  9; return 1;			// t = 21, p+1-t = 9 mod 25
	case 9: a[0] =  4; return 1;			// t = 6, p+1-t = 4 mod 25
	case 14: a[0] =  24; return 1;		// t = 16, p+1-t = 24 mod 25
	case 19: a[0] =  19; return 1;		// t = 1, p+1-t = 19 mod 25
	case 24: a[0] =  14; return 1;		// t = 11, p+1-t = 14 mod 25
	default: err_printf ("Unhandled case p=%lu with r=%d in ecurve_mod5\n", _ff_p, r); exit(0);
	}
}
