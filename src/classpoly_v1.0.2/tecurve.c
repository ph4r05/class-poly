#include <stdlib.h>
#include <stdio.h>
#include "ff_poly.h"
#include "ecurve.h"
#include "tecurve.h"
#include "bigX1.h"
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

#define TECURVE_MAX_TT_N		30

static int genus0_N[TECURVE_MAX_N+1] = { 0,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	

/*
	See Howe "On the group orders of elliptic curves over finite fields",  Theorem 1.1.

	For N<=37 odd, returns approximate density of curves over F_p with order divisible by N,
        and divisible by 3 if N%3 is nonzero and t3>FILTER_3TOR_1, or not divisible by 3 if N%3 is nonzero and t3=FILTER_3TOR_1.
	(if t3==0 then 3-torsion is ignored).
*/


double tecurve_odd_density (long p, int N, int t3)
{
	double r1, r2, r3;
	
	if ( N > 37 || !(N&1) ) { printf ("N=%d must be <= 31 and odd in tecurve_odd_density\n", N); exit (0); }
	switch (N) {
	case 1: r1=1.0; break;
	case 9: r1 = ( (p%3)==1 ? 11.0/72.0 : 1.0/6.0 ); break;
	case 15:
		r1 = ( (p%3)==1 ? 3.0/8.0 : 1.0/2.0 );
		r2 = ( (p%5)==1 ? 5.0/24.0 : 1.0/4.0 );
		r1 *= r2;
		break;
	case 21:
		r1 = ( (p%3)==1 ? 3.0/8.0 : 1.0/2.0 );
		r2 = ( (p%7)==1 ? 7.0/48.0 : 1.0/6.0 );
		r1 *= r2;
		break;
	case 25: r1 = ( (p%5)==1 ? 29.0/600.0 : 1.0/20.0 ); break;
	case 27: r1 = ( (p%9)==1 ? 11.0/216.0 : 1.0/18.0 ); break;
	case 33:
		r1 = ( (p%3)==1 ? 3.0/8.0 : 1.0/2.0 );
		r2 = ( (p%11)==1 ? 11.0/120.0 : 1.0/10.0 );
		r1 *= r2;
		break;
	case 35:
		r1 = ( (p%5)==1 ? 5.0/24.0 : 1.0/4.0 );
		r2 = ( (p%7)==1 ? 7.0/48.0 : 1.0/6.0 );
		r1 *= r2;
		break;
	default:
		r1 = ( (p%N)==1 ? (double)N / (double)(N*N-1) : 1.0/(double)(N-1) );
	}
	r3 = 1.0;
	if ( (N%3) ) {
		if ( t3==FILTER_3TOR_1 ) r3 = ( (p%3)==1 ? 5.0/8.0 : 1.0/2.0 );
		if ( t3>FILTER_3TOR_1 ) r3 = ( (p%3)==1 ? 3.0/8.0 : 1.0/2.0 );
	}
//printf("p=%d, p mod N=%d p mod 3 =%d, N=%d, t3=%d, r1=%f, r3=%f, r1*r3=%f\n", p, p%N, p%3, N, t3, r1, r3, r1*r3);
	return r1*r3;
}


/*
	Returns density of curves with order divisible by N and also 2^k, with t3 as above.
*/
double tecurve_free2_density (long p, int N, int k, int t3, int j_flags)
{
	double r,r2;
	int b,c,j;
	
	for ( j = 0 ; !(N&1) ; j++ ) N>>=1;
	if ( j > k ) k=j;   // if ( j < k ) k=j;	fixed bug 9/2/09
	r = tecurve_odd_density(p,N,t3);
	if  ( ! k ) return r;
	if ( k==1 ) { r *= 2.0/3.0; }
	else if ( k==2 ) { r*= 5.0/12.0; }		// treat k>2 the same as k=2 (not any more 8/29/2009)
	else {
		b = k/2;
		c = ((k&1)?b+1:b);
		if ( (p%(1<<c))==1 ) {
			r2 = (double)((1<<(b+1))+(1<<b) - 1) / (double)((1<<(k+b-1))*3);
		} else {
			r2 = 1.0/(double)(1<<(k-1));
		}
		r *= r2;
	}
	// if ( (j_flags&FILTER_JCUBE) && (p%3)==1 ) r *= 1.0/3.0;			don't use this right now, since we don't account for the cost
	return r;
}
	
/*
	Returns density of curves with order divisible by N and 2^k, but not 2^{k+1}, with t3 as above.
*/
double tecurve_fix2_density (long p, int N, int k, int t3, int j_flags)
{
	double r,r2,r3;
	int b, c, j;
	
	for ( j = 0 ; !(N&1) ; j++ ) N>>=1;
	if ( j > k ) return TECURVE_MAX_COST;
	r = tecurve_odd_density(p,N,t3);
	if ( k==0 ) {
		// distinguish #E=1mod4 vs 3mod4 when p=1mod4
		r2 = ( (p&2) ? 0.5 : 0.25 );		// this is not true but reflects the fact that we only filter out half the cases when k=0, not 2/3 as we could
	} else if ( k == 1 ) {
		r2 = 0.25;		// 2/3 - 5/12
	} else {
		b = k/2;
		c = ((k&1)?b+1:b);
		if ( (p%(1<<c))==1 ) {
			r2 = (double)((1<<(b+1))+(1<<b) - 1) / (double)((1<<(k+b-1))*3);
		} else {
			r2 = 1.0/(double)(1<<(k-1));
		}
		k++;
		b = k/2;
		c = ((k&1)?b+1:b);
		if ( (p%(1<<c))==1 ) {
			r3 = (double)((1<<(b+1))+(1<<b) - 1) / (double)((1<<(k+b-1))*3);
		} else {
			r3 = 1.0/(double)(1<<(k-1));
		}
		r2-=r3;
	}
	// if ( (j_flags&FILTER_JCUBE) && (p%3)==1 ) r *= 1.0/3.0;			don't use this right now, since we don't account for the cost
	return r*r2;
}
	
	
/*
	This module generates curves with prescribed torsion using X_1(N) as described in
	"Constructing elliptic curves with prescribed torsion over finite fields", Sutherland 2008.
	
	We also rely on Atkin and Morain "Finding suitable curves for the elliptic curve factorization method", which provides
	parameterizations that incorporate a base point (x,y) on the curve which is NOT one of the specified torsion points.
*/

static ff_t c_2, c_3, c_4, c_6, c_8, c_9, c_12, c_24, c_27, c_36, c_54, c_81, c_108, c_216, d_4, d_36, d_216, last_p;

void _tecurve_setup_constants()
{
	register ff_t t0;
	
	// precompute small constants using additions to avoid Montgomery reductions.
	_ff_set_one(t0); _ff_add(c_2,t0,t0) ;_ff_add(c_4,c_2,c_2); _ff_add(c_8,c_4,c_4); _ff_add_one(c_9,c_8);
	_ff_add_one(c_3,c_2); _ff_add(c_6,c_3,c_3);  _ff_add(c_12,c_4,c_8);  _ff_add(c_24,c_12,c_12); _ff_add(c_36,c_24,c_12);
	_ff_add(c_27,c_24,c_3);  _ff_add(c_54,c_27,c_27); _ff_add(c_81,c_54,c_27); _ff_add(c_108,c_54,c_54); _ff_add(c_216,c_108,c_108); 
	_ff_mult(t0,_ff_half,_ff_third); _ff_square(d_36,t0); _ff_mult(d_216,t0,d_36);  _ff_square(d_4,_ff_half);
	last_p = _ff_p;
}
static inline void tecurve_setup_constants () { if ( _ff_p != last_p ) _tecurve_setup_constants(); }

// computes discriminant of x^3+f[1]x + f[0], equal to - 4*f[1]^3 - 27*f[0]^2
static inline void tecurve_discriminant_AB (ff_t d[1], ff_t A, ff_t B)
{
	register ff_t t0,t1,t2;
	_ff_square(t0,A);  _ff_mult(t1,t0,A);  _ff_square(t0,B); 
	_ff_mult(t2,c_4,t1); _ff_mult(t1,c_27,t0); _ff_addto(t2,t1); _ff_neg(d[0],t2);
}
static inline void tecurve_discriminant (ff_t d[1], ff_t f[]) { tecurve_discriminant_AB(d,f[1],f[0]); }

static inline void tecurve_AB_eval_df(ff_t y[1], ff_t x, ff_t A)
	{ register ff_t t0,t1;  _ff_square(t0,x); _ff_triple(t1,t0); _ff_add(y[0],t1,A); }

static inline void tecurve_AB_eval_f(ff_t y[1], ff_t x, ff_t A, ff_t B)
	{ register ff_t t0,t1;  _ff_square(t0,x); _ff_addto(t0,A);  _ff_mult(t1,t0,x); _ff_add(y[0],t1,B); }

// converts y^2=4x^3+b2x^2+2b4x+b6 to y^2=x^3+Ax+B using 7M+3A
static inline void tecurve_AB_from_b246 (ff_t A[1], ff_t B[1], ff_t b2, ff_t b4, ff_t b6)
{
	register ff_t t0, t2, c4, c6;
	
	_ff_square(t2,b2);  _ff_mult(c4,c_24,b4); _ff_subfrom(c4,t2);		// compute -c4=24b4-b2^2
	_ff_mult(t0,c_36,b4); _ff_subfrom(t2,t0); _ff_mult(c6,b2,t2);
	_ff_mult(t0,c_216,b6); _ff_addto(c6,t0);						// compute -c6=b2(b2^2-36b4)+216b6
	_ff_mult(A[0],c_27,c4); _ff_mult(B[0],c_54,c6);
}

// converts Tate normal form E(b,c): y^2+(1-c)xy-by=x^3-bx^2 to short Weierstrass form y^2=x^3+Ax+B using 10M+6A
// stores b2 in specified location for reuse
static inline void tecurve_AB_from_bc (ff_t A[1], ff_t B[1], ff_t *b2, ff_t b, ff_t c)
{
	register ff_t t0,b4,b6;
	
	// E(b,c) in Weierstrass notation has a1=1-c, a2=a3=-b and a4=a6=0
	// compute b2=a1^2+4a2 = (c-1)^2-4b,  b4=2a4+a1a3=(c-1)b, b6=a3^2+4a6=b^2
	_ff_sub_one(b6,c); _ff_square(t0,b6); _ff_add(b4,b,b); _ff_x2(b4); _ff_sub(*b2,t0,b4); _ff_mult(b4,b6,b); _ff_square(b6,b);
	tecurve_AB_from_b246 (A, B, *b2, b4, b6);
}

// computes c=s(r-1) and b=rc and then converts E(b,c) to y^2=x^3+Ax+B using 12M+10A
static inline void tecurve_AB_from_rs (ff_t A[1], ff_t B[1], ff_t *b2, ff_t r, ff_t s)
{
	register ff_t b,c;
	_ff_sub_one(b,r); _ff_mult(c,b,s); _ff_mult(b,c,r);
	tecurve_AB_from_bc(A,B,b2,b,c);
}

// translates x-coord on E(b,c) (or b246) to AB curve (only needs b2) using 1M+3A
static inline void tecurve_shift_x_AB_from_bc (ff_t x1[1], ff_t x0, ff_t b2)
	{ register ff_t t0,t1;  _ff_mult(t0,c_12,x0); _ff_addto(t0,b2); _ff_add(t1,t0,t0); _ff_add(x1[0],t0,t1); }

// translates point on b246 curve to AB curve (only needs b2) using 2M+3A
static inline void tecurve_shift_AB_from_b246 (ff_t x1[1], ff_t y1[1], ff_t x0, ff_t y0, ff_t b2)
{
	register ff_t t0,t1;
	_ff_mult(t0,c_12,x0); _ff_addto(t0,b2); _ff_add(t1,t0,t0); _ff_add(x1[0],t0,t1); _ff_mult(y1[0],c_108,y0);
}


int tecurve_genus0_curves (ff_t A[], ff_t B[], ff_t *tt_x, ff_t *tor_x, ff_t *x, ff_t *y, int n, int N);
int tecurve_get_X1_points (ff_t x[], ff_t y[], int n, int N);
int tecurve_map_X1_points (ff_t r[], ff_t s[], ff_t x[], ff_t y[], int n, int N);
int tecurve_get_tt_points (ff_t tt_x[], ff_t r[], ff_t s[], int n, int N);
int tecurve_sws (ff_t A[], ff_t B[], ff_t r[], ff_t s[], ff_t *tt_x, int n, int N);
int tecurve_filter (ff_t A[], ff_t B[], ff_t *tt_x, ff_t *tor_x, ff_t *base_x, ff_t *base_y, int n, int N, int s2_flag, int t3_flag, int j_flag);

int tecurve_random_curves_x (ff_t f1[], ff_t x[], ff_t y[], int n, int N, int s2_flag, int t3_flag, int j_flags)
{
	ff_t f0[TECURVE_MAX_CURVES], r[TECURVE_MAX_CURVES], s[TECURVE_MAX_CURVES], _tt_x[TECURVE_MAX_CURVES], _tor_x[TECURVE_MAX_CURVES];
	ff_t *tt_x, *tor_x, *base_x, *base_y;
	ff_t f[4], b2;
	register int i;
	
	if ( n > TECURVE_MAX_CURVES ) { printf ("Exceeded TECURVE_MAX_CURVES=%d, reduced n\n", TECURVE_MAX_CURVES); n = TECURVE_MAX_CURVES; }
	if ( N > TECURVE_MAX_N ) { printf ("Exceeded TECURVE_MAX_N=%d\n", TECURVE_MAX_N); return 0; }
	
	tt_x = tor_x = base_x = base_y = 0;
	if ( !(N&1) && N <= TECURVE_MAX_TT_N )
		if ( s2_flag == FILTER_2SYLOW_2 || s2_flag == FILTER_2SYLOW_R14 || s2_flag == FILTER_2SYLOW_R2 || s2_flag == FILTER_2SYLOW_D4 || (s2_flag&FILTER_2SYLOW_C_LG) ) tt_x = _tt_x;
	if ( (s2_flag&FILTER_2SYLOW_C_LG) ) tor_x=_tor_x;
	tecurve_setup_constants();
	if ( genus0_N[N] ) {
		base_x  = x; base_y = y; 
		if ( N==2 && s2_flag == FILTER_2SYLOW_R2 ) return 0;											// our 2-torsion method precludes this
		n = tecurve_genus0_curves (f1, f0, tt_x, tor_x, base_x, base_y, n, N);
	} else {
		n = tecurve_get_X1_points (x, y, n, N);
		n = tecurve_map_X1_points (r, s, x, y, n, N);
		if ( tt_x ) n = tecurve_get_tt_points (tt_x, r, s, n, N);
		for ( i = 0 ; i < n ; i++ ) {
			tecurve_AB_from_rs (f1+i,f0+i,&b2,r[i],s[i]);
			if ( tt_x ) tecurve_shift_x_AB_from_bc (tt_x+i,tt_x[i],b2);
			if ( tor_x ) _ff_triple(tor_x[i],b2);
		}
	}
	n = tecurve_filter (f1, f0, tt_x, tor_x, base_x, base_y, n, N, s2_flag, t3_flag, j_flags);							// we could apply j_flags before getting tt points
	if ( ! base_x ) {
		_ff_set_one(f[3]); _ff_set_zero(f[2]);
		for( i = 0 ; i < n ; i++ ) {
			_ff_set(f[1],f1[i]); _ff_set(f[0],f0[i]);
			ecurve_random_point (x+i, y+i, f);
		}
	}
	return n;
}


int tecurve_genus0_curves (ff_t A[], ff_t B[], ff_t *tt_x, ff_t *tor_x, ff_t *x, ff_t *y, int n, int N)
{
	ff_t inv[TECURVE_MAX_CURVES], den[TECURVE_MAX_CURVES];
	ecp_jc_t P[TECURVE_MAX_CURVES];
	ff_t f[4],a[3],a2[3];
	ff_t u,v,b2;
	register ff_t b,c,d,t0,t1,t2,t3,t4,t5,t6;
	register int i,j;
	
	switch (N) {
	case 1:
		ff3_setup();
		_ff_set_one(f[3]); 
		for ( i = 0 ; i < n ; i++ ) {
			_ff_set_one(a[2]);  _ff_random(a[1]); _ff_random(a[0]);					// get random element of F_p^3 which is not in F_p (a[2]==1 guarantees this)
			ff3_norm(y+i,a);												// save its norm
			ff3_square(a,a);												// square it, then compute f(x)=-g(-x) where g(x) is the minpoly of a
			ff3_trace(f+2,a); 												// f[2] = tr(a) = (a+a^p+a^(p^2)) since g[2]= - tr(a)
			ff3_exp_p (a2,a);  ff3_mult (a2,a2,a);								// a2 = a*a^p
			ff3_trace(f+1,a2);												// f[1] = (a*a^p+a^p*a^(p^2)+a^(p^2)*a)
			_ff_square(f[0],y[i]);												// f[0] = N(a) since g[0] = - N(a)
			ff_depress_cubic(x+i,f);											// depress to make f2=0, still irreducible, translates 0 to x
			_ff_set(A[i],f[1]);  _ff_set(B[i],f[0]);
			if ( tor_x ) _ff_set(tor_x[i],x[i]);										// just set torsion point to base point in this case
		}
		break;
	case 2:																// Our construction for 2-torsion guarantees rank 1.  This is not universal, but it is useful.
		ff2_setup();
		_ff_set_one(f[3]); _ff_set_one(a[1]);
		for ( i = 0 ; i< n ; i++ ) {
			_ff_random(a[0]);												// get a random element a of F_p^2 which is not in F_p (a[1]==1 guarantees this)
			ff2_norm(y+i,a);												// construct f(x)=g(x)(x+g(0)) where g(x) is the min poly of a, so (0,g(0)^2) is on y^2=f(x) and (-g(0),0) is 2-torsion pt
			ff2_trace(&u,a);  
			_ff_sub(f[2],y[i],u);
			_ff_set_one(t0);  _ff_subfrom(t0,u); _ff_mult(f[1],t0,y[i]);
			_ff_square(f[0],y[i]);
			ff_depress_cubic(x+i,f);											// depress g to make g2 zero, shifts base point
			_ff_set(A[i],f[1]);  _ff_set(B[i],f[0]);			
			if ( tt_x ) _ff_sub(tt_x[i],x[i],y[i]);									// 2-torsion point, shifted
			if ( tor_x ) _ff_sub(tt_x[i],x[i],y[i]);
		}
		break;
	case 3:																// Our construction for 3-torsion is not universal, since it is a 1-parameter family, but it hits a constant fraction of the possible j-invariants
		for ( i = 0 ; i < n ; i++ ) {	
			_ff_random(t1);												// pick a1 at random to construct y^2+a1xy-a1y=x^3 with 3-torsion pt (0,a1) and base pt (a1,a1)
			_ff_square(b2,t1); _ff_add(t0,b2,c_24); _ff_mult(t3,t0,b2);				// t2=b2=-b4=b6=a1^2, t3=c4=b2^2-24b4=a1^4+24a1^2
			_ff_mult(t0,c_27,t3); _ff_neg(A[i],t0);								// A=-27c4
			_ff_mult(t0,c_12,b2); _ff_addto(t0,t3); _ff_addto(t0,c_216);
			_ff_mult(t3,t0,b2); _ff_mult(B[i],c_54,t3);								// B=-54c6, -c6 = b2(b2^2-36b4)+216b6 = a1^2(c4+12a1^2+216)
			_ff_mult(t0,c_12,t1); _ff_addto(t0,b2); _ff_add(t3,t0,t0);
			_ff_add(x[i],t0,t3); 	_ff_addto(t1,b2); _ff_mult(y[i],c_108,t1);				// base point is (3(12a1+a1^2),108(a1^2+a1))
			if ( tor_x ) _ff_triple(tor_x[i],b2);
		}
		break;
	case 4:	// c=0
		for ( i = 0 ; i < n ; i++ ) {
			_ff_random(d); 
			_ff_square(t0,d); _ff_set_one(t2); _ff_subfrom(t2,t0); _ff_mult (b,t2,_ff_half);	// b = (1-d^2)/2 (random d),  c=0 in Tate normal form.  (d+1,d+1-b) is point of infinite order on E(b,0)
			_ff_add(t0,b,b); _ff_x2(t0); _ff_set_one(b2); _ff_subfrom(b2,t0);			// t2=b2=a1^2+4a2=1-4b, note b4=-b and b6=b^2
			_ff_x2(t0); _ff_add_one(t1,b); _ff_mult(t4,t0,t1); _ff_x2(t4); _ff_inc(t4);		// t4=c4=b2^2-24b4=1+16b(b+1)
			_ff_mult(t1,t4,c_27); _ff_neg(A[i],t1);								// A=-27c4
			_ff_sub(t1,b,c_2); _ff_mult(t3,t0,t1); _ff_addto(t3,b); _ff_subfrom(t3,c_3);	// c6=-b2^3+36b2b4-216b6=8b(8b(b-2)+b-3)-1
			_ff_mult(t1,t0,t3); _ff_dec(t1); _ff_mult(t3,t1,c_54); _ff_neg(B[i],t3);		// B=-54c6
			_ff_inc(d); _ff_mult(t0,c_12,d); _ff_addto(t0,b2); _ff_add(t1,t0,t0);
			_ff_add(x[i],t1,t0); _ff_subfrom(d,b); _ff_add(t0,d,d);					// base point x-coord: 36(d+1)+3b2
			_ff_addto(t0,d); _ff_mult(y[i],c_108,t0);								// base point y-coord: 108(2(d+1-b)+(d+1)-b) = 108*3*(d+1-b)
			if ( tt_x ) tecurve_shift_x_AB_from_bc (tt_x+i, b, b2);
			if ( tor_x ) _ff_triple(tor_x[i],b2);
		}
		break;
	case 5:	// b=c
		for ( i = 0 ; i < n ; i++ ) {
			_ff_random(d); _ff_add(t1,d,d); _ff_addto(t1,d); _ff_inc(t1); _ff_x2(t1);
			_ff_set(A[i],d); _ff_set(B[i],t1);										// pick random d, store it in A, and store 6d+2 in B for inversion
			if ( _ff_zero(t1) ) i--;
		}
		ff_parallel_invert(inv,B,n);												// invert in parallel
		for ( i = 0 ; i < n ; i++ ) {
			_ff_set(d,A[i]); _ff_mult(t0,d,inv[i]);
			_ff_add(t2,d,d); _ff_x2(t2); _ff_inc(t2); _ff_mult(b,t0,t2);					// set b=d(4d+1)/(6d+2), note c=b for 5-torsion
			tecurve_AB_from_bc(A+i,B+i,&b2,b,b);
			_ff_mult(t0,c_12,d); _ff_addto(t0,b2); _ff_add(t3,t0,t0);
			_ff_add(x[i],t0,t3); _ff_inc(d); _ff_mult(t0,b,d);							// base point x-coord: 36d+3b2
			_ff_mult(y[i],c_108,t0);											// base point y-coord: 54b(d+1)    x^2(2x+1)/(3x+1)
			if ( tor_x ) _ff_triple(tor_x[i],b2);
		}
		break;
	case 6:	// s=1, b=c(c+1)
		for ( i = 0 ; i < n ; i++ ) {
			_ff_random(d); _ff_add(t1,d,d); _ff_inc(t1); _ff_square(t2,t1);				// pick random d, compute (d+2)(2d+1)^2 for inversion
			_ff_add(t1,d,c_2); _ff_set(A[i],d); _ff_mult(B[i],t1,t2); 					// store d in A and invesion value in B
			if ( _ff_zero(B[i]) ) i--;
		}
		ff_parallel_invert(inv,B,n);												// invert in parallel
		for ( i = 0 ; i < n ; i++ ) {
			_ff_set(d,A[i]); _ff_add(t0,d,d); _ff_addto(t0,d); _ff_addto(t0,c_2);
			_ff_mult(t1,t0,inv[i]); _ff_x2(t1); _ff_mult(t4,t0,t1);						// t4=2(3d+2)^2/((d+2)(2d+1)^2)=y-coord of base pt on b246 curve y^2=4x^3+b2x^2+2b4x+b6
			_ff_add_one(t2,d); _ff_mult(t3,t1,t2); _ff_sub(c,t3,t4);					// t3=2(d+1)(3d+2)/((d+2)(2d+1)^2)=x-coord of base pt on b246 curve, c = x-y, b = c+c^2
			_ff_square(b,c); _ff_addto(b,c);
			tecurve_AB_from_bc(A+i,B+i,&b2,b,c);								// we could save a mult here by using c^2 but we don't bother
			tecurve_shift_AB_from_b246 (x+i,y+i,t3,t4,b2);
			if ( tt_x ) tecurve_shift_x_AB_from_bc(tt_x+i,c,b2);						// the x-coord of 2-torsion point on E(b,c) is c
			if ( tor_x ) _ff_triple(tor_x[i],b2);
		}
		break;
	case 7:	// r=s
		_ff_set_one(f[3]); _ff_set_zero(f[2]);
		_ff_set_ui(f[1],34965); _ff_set_ui(t0,58266); _ff_neg(f[0],t0);					// setup parametrization curve T^2=S^3 + 34965*S - 58266 derived from Atkins-Morain
		ecurve_random_point(&u,&v,f);
		ecurve_JC_set_A (P,u,v);
		n = ecurve_AJC_steps (P, u, v, n/2, f[1]);									// note that steps will stop if the identity is reached so we only use distinct points (n may be reduced here)
		ecurve_JC_to_A_parallel (A,B,P,n);
		for ( i = 0 ; i < n ; i++ ) { _ff_set(A[n+i],A[i]);  _ff_neg(B[n+i],B[i]); } n *= 2;		// double points by using inverses (which yields non-isomorphic curves)
		_ff_set_ui(t1,3906);													// we could precompute these constants
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_mult(t0,c_6,A[i]); _ff_sub(den[j],t0,t1);							// compute denominator 6S - 3906 we need to invert
			if ( _ff_zero(den[j]) ) {j--; continue; }								// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(A[j],A[i]); _ff_set(B[j],B[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_ui(t5,297);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_mult(t0,c_27,A[i]); _ff_addto(t0,B[i]); _ff_subfrom(t0,t5);				// compute numerater of X: T+27S-297
			_ff_mult(t1,t0,inv[i]);  _ff_square(t4,t1);								// t1 = X = (T+27S-297)/(6S-3906), t4=X^2
			_ff_mult(t2,c_9,t1); _ff_subfrom(t2,t4); _ff_sub(t0,A[i],c_36);
			_ff_subfrom(t0,c_3);  _ff_mult(t3,t0,d_36); _ff_x2(t3);  _ff_addto(t2,t3);		// t2 = Y = 9X - X^2 + (S-39)/18
			_ff_sub(c,t4,t1); _ff_mult(b,c,t1);									// r=s=X, c=s(r-1)=X^2-X, b=cr=cX
			tecurve_AB_from_bc (A+i,B+i,&b2,b,c);
			_ff_mult(t0,t1,t2); _ff_x2(t1); ff_negate(t1);							// t1=x=-2X and t0=y=XY are coords of base pt on b246 curve
			tecurve_shift_AB_from_b246 (x+i, y+i, t1, t0, b2);
			if ( tor_x ) _ff_triple(tor_x[i],b2);
		}
		break;
	case 8: // r=1/X, s=2-X
		_ff_set_one(f[3]); _ff_set_zero(f[2]);
		_ff_neg(f[1],c_27); _ff_set_ui(f[0],8694); 									// setup parametrization curve T^2=S^3 -27*S - 8694 derived a la Atkins-Morain by supposing x=(d-2)/d=-rs
		ecurve_random_point(&u,&v,f);
		ecurve_JC_set_A (P,u,v);
		n = ecurve_AJC_steps (P, u, v, n/2, f[1]);									// note that steps will stop if the identity is reached so we only use distinct points (n may be reduced here)
		ecurve_JC_to_A_parallel (A,B,P,n);
		for ( i = 0 ; i < n ; i++ ) { _ff_set(A[n+i],A[i]);  _ff_neg(B[n+i],B[i]); } n *= 2;		// double points by using inverses (which yields non-isomorphic curves)
		_ff_set_ui(t4,90); _ff_set_ui(t5,243);										// we could precompute these constants but it hardly seems worth doing
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_add(t0,A[i],A[i]); _ff_add(t1,t0,A[i]); _ff_add(t2,t1,t1); _ff_add(t3,t2,t1);	// compute t0=2S, t2=6S and t3=9S
			_ff_addto(t3,B[i]); _ff_sub(A[j],t3,t5); _ff_sub(B[j],t2,t4);					// compute numerator T+9S-243 and denominator 6S-90 of X, we need to save and invert both
			_ff_mult(den[j],A[j],B[j]);											// invert product
			if ( _ff_zero(den[j]) ) {j--; continue; }								// skip zero denominators (should be very rare)
			_ff_set(x[j],t0);													// stash 2S in x array
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_square(t0,A[i]); _ff_square(t1,B[i]);
			_ff_mult(t2,t0,inv[i]); _ff_mult(t3,t1,inv[i]);							// compute t2=X=(T+9S-243)/(6S-90) and t3=1/X=r
			_ff_sub(t4,c_2,t2); _ff_sub_one(t1,t3); _ff_mult(c,t1,t4); _ff_mult(b,t3,c);	// compute t4=s=2-X, c=s(r-1), and b=rc
			tecurve_AB_from_bc (A+i,B+i,&b2,b,c);
			_ff_add(t0,t2,t2); _ff_subfrom(t0,c_3); _ff_square(t1,t0);
			_ff_mult(t0,c_9,t1); _ff_sub(t1,x[i],t0);								// use stashed 2S value in x[i]
			_ff_addto(t1,c_12); _ff_addto(t1,c_3); _ff_mult(t5,d_36,t1);				// compute t5=Y=(2S-9(2X-3)^2+15)/36
			_ff_mult(t0,t3,t4); _ff_neg(t1,t0); _ff_mult(t0,t3,t5); _ff_mult(t2,t0,t1);		// base point is t1=x=-rs and t2=y=rxY
			tecurve_shift_AB_from_b246 (x+i, y+i, t1, t2, b2);
			if ( tt_x ) {
				_ff_sub_one(t1,t3); _ff_mult(t0,t1,t3);							// 2-torsion point 4P has x-coord r(r-1) on E(b,c)
				tecurve_shift_x_AB_from_bc(tt_x+i,t0,b2);
			}
			if ( tor_x ) _ff_triple(tor_x[i],b2);
		}
		break;		
	case 9: //s=X, r=s(s-1)+1
		_ff_set_one(f[3]); _ff_set_zero(f[2]);
		_ff_neg(f[1],c_9); _ff_set(f[0],c_9);										// setup parametrization curve T^2=S^3 - 9*S + 9 from Atkins-Morain
		ecurve_random_point(&u,&v,f);
		ecurve_JC_set_A (P,u,v);
		n = ecurve_AJC_steps (P, u, v, n/2, f[1]);									// note that steps will stop if the identity is reached so we only use distinct points (n may be reduced here)
		ecurve_JC_to_A_parallel (A,B,P,n);
		for ( i = 0 ; i < n ; i++ ) { _ff_set(A[n+i],A[i]);  _ff_neg(B[n+i],B[i]); } n *= 2;	// double points by using inverses (which yields non-isomorphic curves)
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_sub_one(den[j],A[i]);											// need to invert S-1
			if ( _ff_zero(den[j]) ) {j--; continue; }								// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(A[j],A[i]); _ff_set(B[j],B[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_sub_one(t0,B[i]); _ff_mult(t1,t0,inv[i]); _ff_add(t0,t1,c_3);
			_ff_mult(t2,t0,_ff_half); _ff_square(t3,t1); _ff_dec(t3);					// t2=s=X=1/2((T-1)/(S-1)+3)
			_ff_mult(t0,t3,_ff_half); _ff_sub(t3,A[i],t0); _ff_sub_one(t0,t2); 			// t3=Y=S-(((T-1)/(S-1))^2-1)/2
			_ff_mult(t1,t0,t2); _ff_mult(c,t1,t2); _ff_inc(t1); _ff_mult(b,c,t1); 			// t1=r=s(s-1)+1, c=s(r-1), b=cr
			 tecurve_AB_from_bc (A+i,B+i,&b2,b,c);	
			_ff_square(t4,t2); _ff_add(t0,t2,t2); _ff_dec(t0); _ff_mult(t1,t0,t4);			// t1=x0=(2s-1)s^2
			_ff_square(t0,t4); _ff_mult(t2,t0,t3);								// t2=y0=s^4Y/2
			tecurve_shift_AB_from_b246 (x+i, y+i, t1, t2, b2);
			if ( tor_x ) _ff_triple(tor_x[i],b2);
		}
		break;
	case 10: //s=X r=s^2/(s-(s-1)^2)
		_ff_set_one(f[3]); _ff_set_zero(f[2]);
		_ff_set_ui(f[1],864); _ff_set_ui(t0,22896); _ff_neg(f[0],t0);					// setup parametrization curve T^2=S^3 +864S - 22896 derived from Atkins-Morain
		ecurve_random_point(&u,&v,f);
		ecurve_JC_set_A (P,u,v);
		n = ecurve_AJC_steps (P, u, v, n/2, f[1]);							// note that steps will stop if the identity is reached so we only use distinct points (n may be reduced here)
		ecurve_JC_to_A_parallel (A,B,P,n);
		for ( i = 0 ; i < n ; i++ ) { _ff_set(A[n+i],A[i]);  _ff_neg(B[n+i],B[i]); } n *= 2;	// double points by using inverses (which yields non-isomorphic curves)
		_ff_set_ui(t5,105);
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_sub(t0,A[i],t5); _ff_mult(den[j],c_12,t0);							// need to invert 12(S-105)
			if ( _ff_zero(den[j]) ) {j--; continue; }								// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(A[j],A[i]); _ff_set(B[j],B[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_ui(t6,1107); _ff_mult(t4,d_36,d_4); _ff_add_one(t3,d_4);				// setup constants t6=1107, t4=1/144, t3=5/4
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_add(t0,B[i],t6); _ff_mult(t1,t0,inv[i]);	_ff_square(t2,t1);				// compute t1=(T+1107)/(12(S-105))=X-5/4 and t2=(X-5/4)^2
			_ff_addto(t1,t3); _ff_add(t0,A[i],A[i]); _ff_set(A[j],t1);					// compute X=(T+1107)/(12(S-105))+5/4 and save it in A
			_ff_addto(t0,t5); _ff_mult(d,t0,t4); _ff_sub(B[j],d,t2);					// compute Y=(2S+105)/144-(X-5/4)^2 and save it in B
			_ff_sub_one(t0,t1); _ff_square(t2,t0); _ff_sub(den[j],t1,t2);				// compute (s-(s-1)^2) for inversion (needed to compute r)
			if ( _ff_zero(den[j]) ) {j--; continue; }								// skip zero denominators (should be very rare)
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_set(t4,A[i]); _ff_square(t2,t4); _ff_mult(t3,t2,inv[i]); _ff_sub_one(t0,t3);	// set t4=s=X
			_ff_mult(c,t4,t0); _ff_mult(b,c,t3);									// compute t3=r=s^2/(s-(s-1)^2), c=s(r-1), b=cr
			_ff_mult(t0,t3,t4); _ff_neg(t1,t0); _ff_square(t0,t3); _ff_mult(t2,t0,B[i]);		// compute t1=x0=-rs and t2=y0=r^2Y
			 _ff_x2(t2); tecurve_AB_from_bc (A+i,B+i,&b2,b,c);	
			tecurve_shift_AB_from_b246 (x+i, y+i, t1, t2, b2);
			if ( tt_x ) {
				_ff_sub_one(t1,t4); _ff_mult(t2,t1,t4); _ff_mult(t0,t2,t3);				// 2-torsion point 5P has x-coord rs(s-1) on E(b,c)
				tecurve_shift_x_AB_from_bc(tt_x+i,t0,b2);
			}
			if ( tor_x ) _ff_triple(tor_x[i],b2);
		}
		break;
	case 12: // We use a modified form of the Montgomery-Suyama parametrization: y^2=4x^3 - (3a^2+6-1/a^2)*x^2 + 4*a^2x  (i.e., b2= - (3a^2+6-1/a^2), b4=2a^2, and b6=0)
		// note, we can't use inverse points for this parametrization
		_ff_set_one(f[3]); _ff_set_zero(f[2]);
		_ff_neg(f[1],c_12); _ff_set_zero(f[0]);									// use parametrization curve v^2=u^3-12u as suggested on p.243 of [Montgomery 1987]
		ecurve_random_point(&u,&v,f);
		ecurve_JC_set_A (P,u,v);
		n = ecurve_AJC_steps (P, u, v, n, f[1]);								// note that steps will stop if the identity is reached so we only use distinct points (n may be reduced here)
		ecurve_JC_to_A_parallel (A,B,P,n);
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_sub(t0,A[i],c_4); _ff_mult(t1,t0,A[i]); _ff_sub(x[j],t1,c_12);				// numerator of a is u^2-4u-12, store in x
			_ff_add(t0,A[i],c_12); _ff_mult(t1,t0,A[i]); _ff_sub(y[j],t1,c_12);			// denominator of a is u^2+12u-12, store in y
			_ff_mult(den[j],x[j],y[j]);											// we invert the product of the numerator and denominator of a
			if ( _ff_zero(den[j]) ) {j--; continue; }								// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(A[j],A[i]);												// don't bother copying B[i], we don't need to know v
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_square(t0,x[i]); _ff_mult(t6,t0,inv[i]); _ff_square(t2,t6);				// compute t6=a and t2=a^2
			_ff_mult(t5,y[i],inv[i]); _ff_mult(t1,t5,y[i]); _ff_square(t3,t1);				// compute t5=1/(u^2-4u-12) and t3=1/a^2
			_ff_add(t4,t2,t2); _ff_add(t1,t4,t2); _ff_add(t0,t1,c_6); _ff_sub(b2,t3,t0);	// compute b2=1/a^2-3a^2-6 and t4=b4=2a^2, note b6=0
			_ff_add_one(t3,t1); _ff_mult(t1,t3,d_4);								// compute t1=x0=(3a^2+1)/4
			_ff_square(t0,A[i]); _ff_addto(t0,c_12); _ff_dec(t2); _ff_mult(t3,t0,t2);
			_ff_mult(t0,t3,t5); _ff_mult(t2,t0,_ff_half);							// compute t2=y0=(a^2-1)(u^2+12)/(2(u^2-4u-12))
			tecurve_shift_AB_from_b246 (x+i, y+i, t1, t2, b2);
			_ff_mult(t0,c_24,t4); _ff_square(t1,b2); _ff_subfrom(t0,t1);				// compute t0=-c4=24b4-b2^2
			_ff_mult(A[i],c_27,t0); _ff_mult(t0,c_36,t4); _ff_subfrom(t1,t0);			// A=-27c4
			_ff_mult(t0,t1,b2); _ff_mult(B[i],c_54,t0);							// compute t0=-c6=b2(b2^2-36b4) (note b6=0) and B=-54c6
			if ( tt_x ) { _ff_triple(t0,b2); _ff_set(tt_x[i],t0); }						// 2-torsion point on b246 curve is (0,0), so translation is just 3b2
			if ( tor_x ) tecurve_shift_x_AB_from_bc(tor_x+i,t6,b2);					// point of order 4 has x-coord a
		}
		break;
	default:
		printf ("Invalid N=%d in tecurve_genus0_curves\n", N); exit (0);
	}
	return n;
}

int tecurve_get_X1_points (ff_t x[], ff_t y[], int n, int N)
{
	ecp_jc_t P[TECURVE_MAX_CURVES];
	ff_t f[BIGX1_MAX_YDEG+1], r[BIGX1_MAX_YDEG], u, v;
	register ff_t bx,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11;
	register int i,j,k;
	int d;

	// The explicit forms of X_1(N) used below are all taken from [Sutherland 2008] (see Table 3).
	switch (N) {
	case 11: // X_1(11): y^2=x^3-432x+8208 (use short Weierstrass form)
		_ff_set_one(f[3]); _ff_set_zero(f[2]);
		_ff_set_ui(t0,432); _ff_neg(f[1],t0); _ff_set_ui(f[0],8208);														
		ecurve_random_point(&u,&v,f);
		ecurve_JC_set_A (P,u,v);
		n = ecurve_AJC_steps (P, u, v, n/2, f[1]);	
		ecurve_JC_to_A_parallel (x,y,P,n);
		for ( i = 0 ; i < n ; i++ ) { _ff_set(x[n+i],x[i]); _ff_neg(y[n+i],y[i]); } n *= 2;
		break;
	case 13: // X_1(13):   y^2 + (x^3 + x^2 + 1)y - x^2 - x = 0
		// can only use 1 root
		_ff_set_one(f[2]);
		for ( i = 0 ; i < n ; i++ ) {
			do {
				_ff_random(bx);  _ff_add_one(t1,bx); _ff_mult(t0,t1,bx); _ff_neg(f[0],t0);  ff_mult(t2,t0,bx); _ff_add_one(f[1],t2);
			} while ( ! (k=ff_poly_roots_d2(r,f,2)) );
			_ff_set(x[i],bx);  _ff_set(y[i],r[0]);
		}
		break;
	case 14: // X_1(14): y^2 = x^3 - 675x + 13662 (use short Weierstrass form)
		_ff_set_one(f[3]); _ff_set_zero(f[2]);
		_ff_set_ui(t0,675); _ff_neg(f[1],t0); _ff_set_ui(f[0],13662);														
		ecurve_random_point(&u,&v,f);
		ecurve_JC_set_A (P,u,v);
		n = ecurve_AJC_steps (P, u, v, n/2, f[1]);	
		ecurve_JC_to_A_parallel (x,y,P,n);
		for ( i = 0 ; i < n ; i++ ) { _ff_set(x[n+i],x[i]); _ff_neg(y[n+i],y[i]); } n *= 2;
		break;
	case 15: // X_1(15): y^2 = x^3 - 27x + 8694
		_ff_set_one(f[3]); _ff_set_zero(f[2]);
		_ff_neg(f[1],c_27); _ff_set_ui(f[0],8694);			
		ecurve_random_point(&u,&v,f);
		ecurve_JC_set_A (P,u,v);
		n = ecurve_AJC_steps (P, u, v, n/2, f[1]);	
		ecurve_JC_to_A_parallel (x,y,P,n);
		for ( i = 0 ; i < n ; i++ ) { _ff_set(x[n+i],x[i]); _ff_neg(y[n+i],y[i]); } n *= 2;
		break;
	case 16: // X_1(16):  y^2 + (x^3 + x^2 - x + 1)y + x^2
		// can only use 1 root
		_ff_set_one(f[2]);
		for ( i = 0 ; i < n ; i++ ) {
			do {
				_ff_random(bx);  _ff_square(f[0],bx); _ff_add_one(t1,bx); _ff_mult(t2,f[0],t1); _ff_sub_one(t1,bx); _ff_sub(f[1],t2,t1);
			} while ( ! (k=ff_poly_roots_d2(r,f,2)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]);
		}
		break;
	case 17: // X_1(17): y^4 + (x^3 + x^2 - x + 2)y^3 + (x^3 - 3x + 1)y^2 - (x^4 + 2x)y + x^3 + x^2
		_ff_set_one(f[4]);
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx); _ff_square(t0,bx); _ff_mult(t1,t0,bx); _ff_add(f[0],t0,t1); _ff_sub(t3,f[0],bx); _ff_add(f[3],t3,c_2);
				_ff_triple(t3,bx); _ff_sub(t4,t1,t3); _ff_add_one(f[2],t4); _ff_add(t3,t1,c_2); _ff_mult(t4,t3,bx); _ff_neg(f[1],t4);
			} while ( ! (k=ff_poly_roots_d4(r,f)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 18: // X_1(18): y^2 + (x^3 - 2x^2 + 3x + 1)y + 2x
		// can use both roots
		_ff_set_one(f[2]);
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx); _ff_add(f[0],bx,bx); _ff_add(t0,f[0],bx); _ff_square(t1,bx); _ff_sub(t2,bx,c_2); _ff_mult(t3,t1,t2); _ff_addto(t3,t0); _ff_add_one(f[1],t3);
			} while ( ! (k=ff_poly_roots_d2(r,f,2)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;	
	case 19: // X_1(19): y^5 - (x^2 + 2)y^4 - (2x^3 + 2x^2 + 2x - 1)y^3 + (x^5 + 3x^4 + 7x^3 + 6x^2 + 2x)y^2 - (x^5 + 2x^4 + 4x^3 + 3x^2)y + x^3 + x^2
		_ff_set_one(f[5]);
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx); _ff_square(t2,bx); _ff_add(t0,t2,c_2); _ff_neg(f[4],t0); _ff_mult(t3,t2,bx); _ff_add(f[0],t2,t3); _ff_add(t0,f[0],bx);
				_ff_x2(t0); _ff_sub_one(t4,t0); _ff_neg(f[3],t4); _ff_square(t4,t2); _ff_mult(t5,t4,bx); _ff_add(t6,t4,t4); _ff_triple(t0,f[0]);
				_ff_add(t1,t5,t6); _ff_addto(t1,t0); _ff_addto(t1,t3); _ff_neg(f[1],t1); _ff_addto(t1,t4); _ff_addto(t1,t0); _ff_add(t0,bx,bx); _ff_add(f[2],t0,t1);
			} while ( ! (k=ff_poly_distinct_roots(r,f,5)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;	
	case 20: // X_1(20): y^3 + (x^2 + 3)y^2 + (x^3 + 4)y + 2
		_ff_set_one(f[3]); _ff_set(f[0],c_2);
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx); _ff_square(t2,bx); _ff_add(f[2],t2,c_3); _ff_mult(t3,t2,bx); _ff_add(f[1],t3,c_4);
			} while ( ! (k=ff_poly_roots_d3(r,f)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;	
	case 21: // X_1(21): y^4 + (3x^2 + 1)y^3 + (x^5 + x^4 + 2x^2 + 2x)y^2 + (2x^4 + x^3 + x)y + x^3
		_ff_set_one(f[4]);
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx); _ff_square(t2,bx); _ff_add(t1,t2,t2); _ff_add(t0,t1,t2); _ff_add_one(f[3],t0); _ff_mult(f[0],t2,bx); _ff_square(t4,t2); _ff_add(t0,t4,t4);
				_ff_addto(t0,f[0]); _ff_add(f[1],t0,bx); _ff_add(t0,bx,bx); _ff_addto(t1,t0); _ff_addto(t1,t4); _ff_mult(t5,t4,bx); _ff_add(f[2],t1,t5);
			} while ( ! (k=ff_poly_roots_d4(r,f)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;	
	case 22: // X_1(22):  y^4 + (x^3 + 2x^2 + x + 2)y^3 + (x^5 + x^4 + 2x^3 + 2x^2 + 1)y^2 + (x^5 - x^4 - 2x^3 - x^2 - x)y - x^4 - x^3
		_ff_set_one(f[4]);
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx); _ff_square(t2,bx); _ff_mult(t3,t2,bx); _ff_square(t4,t2); _ff_mult(t5,t3,t2); _ff_add(t1,t3,t4); _ff_neg(f[0],t1);
				_ff_addto(t1,t2); _ff_addto(t1,t3); _ff_add(t0,t1,bx);  _ff_sub(f[1],t5,t0); _ff_addto(t1,t2); _ff_inc(t1); _ff_add(f[2],t5,t1);
				_ff_x2(t2); _ff_add(t0,t2,t3); _ff_add(t1,bx,c_2); _ff_add(f[3],t0,t1);
			} while ( ! (k=ff_poly_roots_d4(r,f)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;	
	case 23: // X_1(23): y^7 + (x^5 - x^4 + x^3 + 4x^2 + 3)y^6 + (x^7 + 3x^5 + x^4 + 5x^3 + 7x^2 - 4x + 3)y^5 + (2x^7 + 3x^5 - x^4 - 2x^3 - x^2 - 8x + 1)y^4
		       //                       + (x^7 - 4x^6 - 5x^5 - 6x^4 - 6x^3 - 2x^2 - 3x)y^3 - (3x^6 - 5x^4 - 3x^3 - 3x^2 - 2x)y^2 + (3x^5 + 4x^4 + x)y - (x^4+2x^3+x^2)
		_ff_set_one(f[7]);
		for ( i = 0 ; i < n ; ) {
			do {
				// we could be slightly more efficient here, but this is all negligible compared to finding the root of a degree 7 poly
				// we cache t11=2x, t10=3x^6, t9=2x^2, t8=3x^5 and generally use tn=x^n (but reuse t6 and t11).  t0 and t1 are work variables
				_ff_random(bx); _ff_square(t2,bx); _ff_mult(t3,t2,bx); _ff_square(t4,t2); _ff_mult(t5,t4,bx); _ff_sub(t0,t5,t4); _ff_addto(t0,t3); _ff_add(t9,t2,t2); _ff_add(t6,t9,t9);
				_ff_addto(t0,t6); _ff_add(f[6],t0,c_3); _ff_triple(t8,t5); _ff_mult(t7,t4,t3); _ff_add(t0,t7,t8); _ff_addto(t0,t4); _ff_triple(t1,t3); _ff_addto(t0,t1); _ff_addto(t0,t3);
				_ff_addto(t0,t3); _ff_addto(t0,t6); _ff_addto(t0,t9); _ff_addto(t0,t2); _ff_add(t11,bx,bx); _ff_add(t6,t11,t11); _ff_subfrom(t0,t6); _ff_add(f[5],t0,c_3); _ff_add(t0,t7,t7);
				_ff_addto(t0,t8); _ff_subfrom(t0,t4); _ff_add(t1,t3,t3); _ff_subfrom(t0,t1); _ff_subfrom(t0,t2); _ff_x2(t6); _ff_subfrom(t0,t6); _ff_add_one(f[4],t0); _ff_add(t0,t4,t1);
				_ff_addto(t0,t2); _ff_neg(f[0],t0); _ff_add(t0,t4,t4); _ff_add(t1,t0,t0); _ff_add(t0,t8,t1); _ff_add(f[1],t0,bx); _ff_mult(t6,t2,t4); _ff_triple(t10,t6); _ff_sub(t0,t10,t1);			
			        _ff_addto(t0,f[0]); _ff_subfrom(t0,t3); _ff_subfrom(t0,t9) _ff_subfrom(t0,t11); _ff_neg(f[2],t0); _ff_sub(t0,t7,t10); _ff_subfrom(t0,t6); _ff_subfrom(t0,f[1]);
				_ff_subfrom(t0,t5); _ff_subfrom(t0,t5); _ff_addto(t0,f[0]); _ff_addto(t0,f[0]); _ff_subfrom(t0,t3); _ff_subfrom(t0,t3); _ff_sub(f[3],t0,t11);
			} while ( ! (k=ff_poly_distinct_roots(r,f,7)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 24: // X_1(24): y^5 + (x^4 + 4x^3 + 3x^2 - x - 2)y^4 - (2x^4 + 8x^3 + 7x^2 - 1)y^3 - (2x^5 + 4x^4 - 3x^3 - 5x^2 - x)y^2 + (2x^5 + 5x^4 + 2x^3)y + x^6 + x^5
		_ff_set_one(f[5]);
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx); _ff_square(t2,bx); _ff_mult(t3,t2,bx); _ff_square(t4,t2); _ff_add(t8,t3,t3); _ff_add(t1,t8,t8); _ff_add(t0,t4,t1); _ff_add(t9,t2,t2);
				_ff_add(t10,t9,t2); _ff_addto(t0,t10); _ff_add(t1,t0,t0); _ff_subfrom(t0,bx); _ff_sub(f[4],t0,c_2); _ff_add(t0,t1,t2); _ff_dec(t0); _ff_neg(f[3],t0);
				_ff_mult(t5,t2,t3); _ff_add(t6,t5,t5); _ff_add(t1,t4,t4); _ff_add(t0,t1,t1); _ff_add(t1,t0,t6); _ff_sub(t0,t1,t8); _ff_subfrom(t0,t3); _ff_subfrom(t0,t9);
				_ff_subfrom(t0,t10); _ff_subfrom(t0,bx); _ff_neg(f[2],t0);  _ff_addto(t1,t4); _ff_add(f[1],t1,t8); _ff_mult(t6,t5,bx); _ff_add(f[0],t5,t6);
			} while ( ! (k=ff_poly_distinct_roots(r,f,5)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 25: //X_1(25): y^8 + (4x^2 + 7x - 4)y^7 - (x^5 - x^4 - 14x^3 - 4x^2 + 24x - 6)y^6 - (x^7 + 4x^6 - 3x^5 - 18x^4 + 15x^3 + 33x^2 - 30x + 4)y^5
			//                      -  (x^8 + 2x^7 - 8x^6 - 14x^5 + 24x^4 + 17x^3 - 41x^2 + 16x - 1)y^4 + (x^8 + 6x^7 + 3x^6 - 20x^5 - 3x^4 + 28x^3 - 19x^2 + 3x)y^3
			//                      -  (3x^7 + 9x^6 - 3x^5 - 13x^4 + 11x^3 - 3x^2)y^2 + (3x^6 + 4x^5 - 4x^4 + x^3)y - x^5;
		_ff_set_one(f[8]);
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx); _ff_triple(t3,bx); _ff_add(t9,t3,t3); _ff_add(t7,t9,bx); _ff_add(t0,t7,bx); _ff_neg(t1,t0); _ff_add(t4,t1,t1); _ff_add(t6,t4,t1);		// t9=2x, t0=8x. t1=-8x
							_ff_sub(t0,t6,t9); _ff_neg(t5,t0);																			// x: f7=7, f6=-24, f5=30, f4=-16, f3=3
				_ff_square(t11,bx); _ff_triple(t2,t11); _ff_add(t0,t2,t11); _ff_addto(t6,t0); _ff_addto(t7,t0); _ff_sub(f[7],t7,c_4); _ff_add(t1,t0,t0); 				// t0=4x^2, t1=8x^2
								_ff_add(t9,t1,t1); _ff_add(t10,t9,t2); _ff_subfrom(t3,t10); _ff_add(t0,t9,t9); _ff_add(t9,t0,t11); _ff_subfrom(t5,t9);		// t9=16x^2, t10=19x^2, t9=33x^2
								_ff_addto(t9,t1); _ff_addto(t4,t9);																		// x^2: f7=f6=4, f5=-33, f4=41, f3=-19, f2=3
				_ff_mult(t1,t11,bx); _ff_add(t7,t1,t1); _ff_add(t9,t7,t1); _ff_add(t8,t7,t7); _ff_x2(t8); _ff_add(t0,t8,t9); _ff_subfrom(t2,t0);					// t1=x^3, t7=2x^3, t8=8x^3, t9=3x^3, t0=11x^3, t10=14x^3
							       _ff_add(t10,t0,t9); _ff_addto(t6,t10); _ff_add(t0,t10,t1); _ff_subfrom(t5,t0); _ff_addto(t0,t7); _ff_subfrom(t4,t0);			// x^3: f6=14, f5=-15, f4=-17, f3=28, f2=-11, f1=1
							       _ff_x2(t10); _ff_addto(t3,t10);
				_ff_mult(t11,t1,bx); _ff_addto(t6,t11); _ff_add(t7,t11,t11); _ff_add(t8,t7,t11); _ff_subfrom(t3,t8); _ff_add(t9,t7,t7); _ff_subfrom(t1,t9);		// t7=2x^4, t8=3x^4, t9=4x^4, t10=6x^4, t7=12x^2
							       _ff_add(t10,t8,t8); _ff_add(t7,t10,t10); _ff_add(t0,t7,t11); _ff_addto(t2,t0); _ff_add(t0,t7,t10); _ff_addto(t5,t0);			// x^4: f6=1, f5=18, f4=24, f3=-3, f2=13, f1=-4
							       _ff_addto(t0,t10); _ff_subfrom(t4,t0);
				ff_mult(t11,t11,bx); _ff_neg(f[0],t11); _ff_subfrom(t6,t11); _ff_add(t7,t11,t11); _ff_add(t8,t7,t11); _ff_addto(t2,t8); _ff_addto(t5,t8);			// t7=2x^5, t8=3x^5, t9=6x^5
								_ff_add(t0,t8,t11); _ff_addto(t1,t0); _ff_add(t9,t8,t8); _ff_add(t0,t9,t9); _ff_addto(t0,t7); _ff_addto(t4,t0);
								_ff_addto(t0,t9); _ff_subfrom(t3,t0);																		// x^5: f6=-1, f5=3, f4=14, f3=-20, f2=3, f1=4, f0=-1
				ff_mult(t11,t11,bx); _ff_add(t7,t11,t11); _ff_add(t0,t11,t7); _ff_add(f[1],t1,t0); _ff_addto(t3,t0); _ff_addto(t0,t11); _ff_subfrom(t5,t0);
								_ff_x2(t0); _ff_addto(t4,t0); _ff_addto(t0,t11); _ff_subfrom(t2,t0);												// x^6: f5=-4, f4=8, f3=3, f2=-9, f1=3
				ff_mult(t11,t11,bx); _ff_subfrom(t5,t11); _ff_add(t0,t11,t11); _ff_subfrom(t4,t0); _ff_addto(t0,t11); _ff_sub(f[2],t2,t0); _ff_x2(t0); _ff_addto(t3,t0); // x^7: f5=-1, f4=-2, f3=6, f2=-3
				ff_mult(t11,t11,bx); _ff_subfrom(t4,t11); _ff_add(f[3],t3,t11); _ff_add(f[6],t6,c_6); _ff_sub(f[5],t5,c_4); _ff_add_one(f[4],t4);					// x^8: f4=-1, f3=1, x^0: f4=1, f5=-1, f6=6, f7=-4 (set above)
			} while ( ! (k=ff_poly_distinct_roots(r,f,8)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}	
		break;
	case 26: // X_1(26): y^6 + (3x^2 + 4x - 2)y^5 + (3x^4 + 10x^3 - 9x + 1)y^4 + (x^6 + 7x^5 + 8x^4 - 14x^3 - 11x^2 + 6x)y^3
		       //                       + (x^7 + 4x^6 - x^5 - 13x^4 + 2x^3 + 10x^2 - x)y^2 - (x^6 - 7x^4 - 4x^3 + 2x^2)y - x^4 - x^3;
		_ff_set_one(f[6]); _ff_set_ui(t7,7);
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx);  _ff_neg(t2,bx); _ff_add(t8,bx,bx); _ff_add(t5,t8,t8); _ff_add(t3,t5,t8); _ff_add(t9,t3,t8); _ff_addto(t9,bx); _ff_neg(t4,t9);		// x: f5=4, f4=-9, f3=6, f2=-1, f1=f0=0
				_ff_square(t11,bx); _ff_add(t10,t11,t11); _ff_neg(t1,t10); _ff_add(t8,t10,t11); _ff_addto(t5,t8); _ff_triple(t9,t8); _ff_addto(t9,t11);
							       _ff_addto(t2,t9); _ff_addto(t9,t11); _ff_neg(t8,t9); _ff_addto(t3,t8);												// x^2: f5=3, f4=0, f3=-11, f2=10, f1=-2, f0=0
				ff_mult(t11,t11,bx);  _ff_neg(t0,t11); _ff_add(t9,t11,t11); _ff_addto(t2,t9); _ff_add(t8,t9,t9); _ff_addto(t1,t8); _ff_add(t10,t8,t8);
						                 _ff_addto(t10,t9); _ff_addto(t4,t10); _ff_addto(t10,t8); _ff_subfrom(t3,t10);										// x^3: f5=0, f4=10, f3=-14, f2=2, f1=4, f0=-1
				ff_mult(t11,t11,bx); _ff_sub(f[0],t0,t11); _ff_add(t8,t11,t11); _ff_add(t9,t8,t11); _ff_addto(t4,t9); _ff_add(t10,t8,t9); _ff_add(t0,t8,t10);		// note t0 is now free
								_ff_addto(t1,t0); _ff_addto(t0,t11); _ff_addto(t3,t0); _ff_addto(t0,t10); _ff_subfrom(t2,t0);							// x^4: f5=0, f4=3, f3=8, f2=-13, f1=7, f0=-1
				ff_mult(t11,t11,bx); _ff_subfrom(t2,t11); _ff_mult(t0,t7,t11); _ff_addto(t3,t0);													// x^5: f5=f4=0, f3=7, f2=-1, f1=f0=0
				ff_mult(t11,t11,bx); _ff_add(f[3],t3,t11); _ff_sub(f[1],t1,t11); _ff_add(t0,t11,t11); _ff_x2(t0); _ff_addto(t2,t0);							// x^6: f5=f4=0, f3=1, f2=4, f1=-1, f0=0
				ff_mult(t11,t11,bx); _ff_add(f[2],t2,t11);																					// x^7: f5=f4=f3=0, f2=1, f0=0
				_ff_sub(f[5],t5,c_2); _ff_add_one(f[4],t4);																					// x^0: f5=-2, f4=1, f3=f2=f1=f0=0
			} while ( ! (k=ff_poly_distinct_roots(r,f,6)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 27:
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx);  d=X27_eval(f,bx);
			} while ( ! (k=ff_poly_distinct_roots(r,f,d)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 28:
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx);  d=X28_eval(f,bx);
			} while ( ! (k=ff_poly_distinct_roots(r,f,d)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 29:
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx);  d=X29_eval(f,bx);
			} while ( ! (k=ff_poly_distinct_roots(r,f,d)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 30:
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx);  d=X30_eval(f,bx);
			} while ( ! (k=ff_poly_distinct_roots(r,f,d)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 31:
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx);  d=X31_eval(f,bx);
			} while ( ! (k=ff_poly_distinct_roots(r,f,d)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 32:
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx);  d=X32_eval(f,bx);
			} while ( ! (k=ff_poly_distinct_roots(r,f,d)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 33:
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx);  d=X33_eval(f,bx);
			} while ( ! (k=ff_poly_distinct_roots(r,f,d)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 34:
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx);  d=X34_eval(f,bx);
			} while ( ! (k=ff_poly_distinct_roots(r,f,d)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 35:
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx);  d=X35_eval(f,bx);
			} while ( ! (k=ff_poly_distinct_roots(r,f,d)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	case 36:
		for ( i = 0 ; i < n ; ) {
			do {
				_ff_random(bx);  d=X36_eval(f,bx);
			} while ( ! (k=ff_poly_distinct_roots(r,f,d)) );
			_ff_set(x[i],bx); _ff_set(y[i],r[0]); i++;
			for ( j=1 ; j < k && i < n ; i++, j++) { _ff_set(x[i],bx); _ff_set(y[i],r[j]); }
		}
		break;
	default:
		printf ("Invalid N=%d in tecurve_get_X1_points\n", N); exit (0);
	}
	return n;
}

int tecurve_map_X1_points (ff_t r[], ff_t s[], ff_t x[], ff_t y[], int n, int N)
{
	ff_t inv[TECURVE_MAX_CURVES], den[TECURVE_MAX_CURVES], rd[TECURVE_MAX_CURVES], sd[TECURVE_MAX_CURVES];
	ff_t c0,c1,c2,c3,c4,c5,c6,c7;		// trust the compiler to put these in registers when and if it makes sense
	register ff_t t0,t1,t2,t3,t4;
	register int i, j;
	
	// We map the point (x,y) on the optimized plane model for X_1(N) to a pair (r,s) defining a curve in Tate normal form E(b,c) with c=s(r-1), b=cr
	// These are taken from [Sutherland 2008] (see Table 3 for N=11,14,15 and Table 7 for the rest).
	switch (N) {
	case 11: // r=(y+108)/216,  s=1+(y-108)/(6x+72)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_add(t0,x[i],c_12); _ff_mult(den[j],t0,c_6);																	// denominator we need to invert is 6(x+12)
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(y[j],y[i]);																						// don't bother copying x, we don't need it anymore
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add(t0,y[i],c_108); _ff_mult(r[i],d_216,t0);																	// r=(y+108)/216
			_ff_sub(t0,y[i],c_108); _ff_mult(t1,t0,inv[i]); _ff_add_one(s[i],t1);													// s = 1=(y-108)/(6x+72)
		}
		break;
	case 13: // r = 1 - xy, s = 1 - xy/(y + 1)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_add_one(den[j],y[i]);																					// denominator we need to invert is 6(y+1)
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_mult(t0,x[i],y[i]); _ff_sub(r[i],t1,t0); _ff_mult(t2,t0,inv[i]); _ff_sub(s[i],t1,t2);										// r=1-xy, s=1-xy/(y+1)
		}
		break;		
	case 14: // r = 1 + (108x - 36y + 3564)/(3x^2 - xy - 342x + 75y + 999), s = (6x - 234)/(9x - y - 135)
		_ff_set_ui(c0,342); _ff_set_ui(c1,75); _ff_set_ui(c2,999); _ff_set_ui(c3,135);										
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_add(t0,x[i],x[i]); _ff_addto(t0,x[i]); _ff_subfrom(t0,y[i]); _ff_subfrom(t0,c0); _ff_mult(t1,t0,x[i]); _ff_mult(t0,c1,y[i]);			// rd = (3x^2 - xy - 342x + 75y + 999)
			_ff_addto(t1,t0); _ff_add(rd[j],t1,c2); _ff_mult(t0,c_9,x[i]); _ff_subfrom(t0,y[i]); _ff_sub(sd[j],t0,c3); _ff_mult(den[j],rd[j],sd[j]);	// sd = (9x - y - 135), invert product of rd and sd
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_ui (c0,234); _ff_set_ui(c1,3564);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_mult(t1,rd[i],inv[i]); _ff_mult(t0,c_6,x[i]); _ff_subfrom(t0,c0); _ff_mult(s[i],t0,t1);									// s = (6x-234) / sd
			_ff_mult(t1,sd[i],inv[i]); _ff_add(t0,x[i],x[i]); _ff_addto(t0,x[i]); _ff_subfrom(t0,y[i]); _ff_mult(t2,t0,c_36); _ff_addto(t2,c1);
			_ff_mult(t0,t1,t2); _ff_add_one(r[i],t0);																		// r = 1 + (36(3x-y)+3564) / rd
		}
		break;
	case 15: // r = 1 - t/(x^2 y - 189x^2 + 42xy - 4050x - 3y^2 + 441y - 1701), s = 1 - t/(x^2 y - 81x^2 + 6xy - 3402x - 3y^2 + 981y - 35721), t = (6x - 90)(18x + 6y - 918)
		_ff_set_ui(c0,189); _ff_set_ui(c1,42); _ff_set_ui(c2,4050); _ff_set_ui(c3,147); _ff_set_ui(c4,1701); _ff_set_ui(c5,3402); _ff_set_ui(c6,327); _ff_set_ui(c7,35721);
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub(t0,y[i],c0); _ff_mult(t1,t0,x[i]); _ff_mult(t0,c1,y[i]); _ff_addto(t1,t0); _ff_subfrom(t1,c2); _ff_mult(t2,t1,x[i]);
			_ff_sub(t1,c3,y[i]); _ff_mult(t0,t1,y[i]); _ff_triple(t1,t0); _ff_subfrom(t1,c4); _ff_add(rd[j],t1,t2);							// rd = x(x(y - 189) + 42y - 4050) + 3y(147-y)- 1701
			_ff_sub(t0,y[i],c_81); _ff_mult(t1,t0,x[i]); _ff_mult(t0,c_6,y[i]); _ff_addto(t1,t0); _ff_subfrom(t1,c5); _ff_mult(t2,t1,x[i]);
			_ff_sub(t1,c6,y[i]); _ff_mult(t0,t1,y[i]); _ff_triple(t1,t0); _ff_subfrom(t1,c7); _ff_add(sd[j],t1,t2);   							// sd = x(x(y - 81) + 6y - 3402) + 3y(327-y) - 35721
			_ff_mult(den[j],rd[j],sd[j]);																				// invert product of rd and sd
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_ui(c0,90); _ff_set_ui(c1,918); _ff_set_one(c2);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_mult(t0,c_6,x[i]); _ff_triple(t1,t0); _ff_subfrom(t0,c0); _ff_mult(t2,c_6,y[i]); _ff_addto(t1,t2); _ff_subfrom(t1,c1); _ff_mult(t2,t0,t1);	// t =  (6x - 90)(18x + 6y - 918)
			_ff_mult(t0,sd[i],inv[i]); _ff_mult(t1,t0,t2); _ff_sub(r[i],c2,t1);  _ff_mult(t0,rd[i],inv[i]); _ff_mult(t1,t0,t2); _ff_sub(s[i],c2,t1);		// r = 1 - t /rd, s = 1 -t /sd
		}
		break;
	case 16: // r = (x^2 - xy + y^2 + y)/(x^2 + x - y - 1), s = (x - y)/(x + 1)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_square(t0,x[i]); _ff_addto(t0,x[i]); _ff_subfrom(t0,y[i]); _ff_sub_one(rd[j],t0);  _ff_add_one(sd[j],x[i]);						// rd = x^2x-y-1, sd=x+1
			_ff_mult(den[j],rd[j],sd[j]);																				// invert product of rd and sd
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_mult(t1,inv[i],rd[i]); _ff_sub(t0,x[i],y[i]); _ff_mult(s[i],t0,t1); _ff_mult(t2,x[i],t0); _ff_square(t0,y[i]);						// s = (x-y) / sd
			_ff_addto(t2,t0); _ff_addto(t2,y[i]); _ff_mult(t1,inv[i],sd[i]); _ff_mult(r[i],t1,t2);										// r = (x^2-xy+y^2+y) / rd
		}
		break;
	case 17: // r = (x^2 + x - y)/(x^2 + xy + x - y^2 - y), s =(x + 1)/(x + y + 1)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_add(t0,x[i],y[i]); _ff_add_one(sd[j],t0); _ff_mult(t0,sd[j],x[i]); _ff_square(t1,y[i]); _ff_subfrom(t0,t1); _ff_sub(rd[j],t0,y[i]);		// sd=x+y+1, rd=x*sd-y^2-y
			_ff_mult(den[j],rd[j],sd[j]);																				// invert product of rd and sd
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_mult(t1,inv[i],rd[i]); _ff_add_one(t0,x[i]); _ff_mult(s[i],t0,t1); 													// s = (x+1) / sd
			_ff_mult(t2,x[i],t0); _ff_subfrom(t2,y[i]); _ff_mult(t1,inv[i],sd[i]); _ff_mult(r[i],t1,t2);									// r = (x(x+1)-y) / rd
		}
		break;
	case 18: // r18:=(x^2 - xy - 3x + 1)/((x-1)^2(xy+1)), s= (x^2 - 2x - y)/(x^2 - xy - 3x - y^2 - 2y);
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub_one(t0,x[i]); _ff_square(t1,t0); _ff_mult(t2,x[i],y[i]); _ff_add_one(t3,t2); _ff_mult(rd[j],t1,t3);						// rd = (x-1)^2(xy+1)
			_ff_subfrom(t1,t2); _ff_add_one(t0,y[i]); _ff_square(t2,t0); _ff_subfrom(t1,t2); _ff_sub(sd[j],t1,x[i]);						// sd = (x-1)^2 - xy - (y+1)^2 - x
			_ff_mult(den[j],rd[j],sd[j]);																				// invert product of rd and sd
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add_one(t0,y[i]); _ff_square(t2,t0); _ff_addto(t2,sd[i]); _ff_mult(t0,sd[i],inv[i]); _ff_mult(r[i],t0,t2);						// r = (sd+(y+1)^2) / rd
			_ff_sub(t0,x[i],c_2); _ff_mult(t1,t0,x[i]); _ff_subfrom(t1,y[i]); _ff_mult(t0,rd[i],inv[i]); _ff_mult(s[i],t0,t1);					// s = (x(x-2)-y) / sd
		}
		break;
	case 19: // r = 1+x(x+y)(y-1)/((x+1)(x^2-xy+2x-y^2+y)), s =1+x(y-1)/((x+1)(x-y+1))
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_add_one(t0,x[i]); _ff_sub(t1,t0,y[i]); _ff_mult(sd[j],t0,t1); _ff_sub_one(t1,y[i]); _ff_square(t2,t1); _ff_sub(t1,sd[j],t2);			// sd=(x+1)*(x+1-y), rd=(x+1)(sd-(y+1)^2)
			_ff_mult(rd[j],t0,t1); _ff_mult(den[j],rd[j],sd[j]);																// invert product of rd and sd
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_sub_one(t0,y[i]); _ff_mult(t1,t0,x[i]); _ff_mult(t0,rd[i],inv[i]); _ff_mult(t2,t0,t1); _ff_add_one(s[i],t2);						// s=1+x(y-1) / sd
			_ff_add(t0,x[i],y[i]); _ff_mult(t2,t0,t1); _ff_mult(t0,sd[i],inv[i]); _ff_mult(t1,t0,t2); _ff_add_one(r[i],t1);						// r = 1+(x+y)x(y-1) / rd
		}
		break;
	case 20: // r=1+(x^3+xy+x)/((x-1)^2(x^2-x+y+1)), s=1+(x^2+y+1)/((x-1)(x^2-x+y+2));
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub_one(t1,x[i]); _ff_square(t0,x[i]); _ff_subfrom(t0,x[i]); _ff_addto(t0,y[i]); _ff_inc(t0); _ff_add_one(t2,t0);				// sd=(x-1)(x^2-x+y+2)
			_ff_mult(sd[j],t1,t2); _ff_square(t2,t1); _ff_mult(rd[j],t0,t2);  _ff_mult(den[j],rd[j],sd[j]);									// rd=(x-1)^2(x^2-x+y+1)
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_square(t0,x[i]); _ff_addto(t0,y[i]); _ff_inc(t0); _ff_mult(t1,rd[i],inv[i]); _ff_mult(t2,t0,t1); _ff_add_one(s[i],t2);				// s = 1 + (x^2+y+1) / sd
			_ff_mult(t1,t0,x[i]); _ff_mult(t0,sd[i],inv[i]); _ff_mult(t2,t0,t1); _ff_add_one(r[i],t2);										// r = 1 + x(x^2+y+1) / rd
		}
		break;
	case 21: // r=1+(y^2+y)(xy+y+1)/((xy+1)(xy+1-y^2)), s=1+(y^2+y)/(xy+1)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_mult(t0,x[i],y[i]); _ff_add_one(sd[j],t0);  _ff_square(t1,y[i]); _ff_sub(rd[j],sd[j],t1); _ff_mult(den[j],rd[j],sd[j]);				// sd=xy+1, rd=sd*(sd-y^2) since sd divides rd, we invert rd and actually store rd/sd in rd
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_square(t0,y[i]); _ff_addto(t0,y[i]); _ff_mult(t1,inv[i],rd[i]); _ff_mult(t2,t0,t1); _ff_add_one(s[i],t2);						// s = 1 + (y^2+y) / sd
			_ff_add(t1,sd[i],y[i]); _ff_mult(t2,t0,t1); _ff_mult(t0,t2,inv[i]); _ff_add_one(r[i],t0);										// r = 1 + (y^2+y)(xy+1+y) / rd (where rd is inv and xy+1 is taken from sd)
		}		
		break;
	case 22: // r=(x^2y + x^2 + xy + y)/(x^3 + 2x^2 + y) s=(xy + y)/(x^2 + y);
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_square(t0,x[i]); _ff_add(t1,t0,y[i]); _ff_set(sd[j],t1); _ff_addto(t1,t0); _ff_mult(t2,t0,x[i]); _ff_add(rd[j],t1,t2);				// sd=x^2+y, rd=sd+x^3+x^2
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add_one(t0,x[i]); _ff_mult(t1,t0,y[i]); _ff_mult(t0,inv[i],rd[i]); _ff_mult(s[i],t0,t1);									// s = (xy+y) / sd
			_ff_add_one(t0,y[i]); _ff_square(t2,x[i]); _ff_mult(t3,t0,t2); _ff_addto(t1,t3); _ff_mult(t0,inv[i],sd[i]); _ff_mult(r[i],t0,t1);			// r = (x^2(y+1)+(xy+y)) / rd
		}		
		break;
	case 23: // r = (x^2 + x + y + 1)/(x^2- x*y),  s = (x + y + 1)/x
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub(rd[j],x[i],y[i]); _ff_mult(den[j],rd[j],x[i]);																// sd=x, rd=sd*(x-y), since sd divides rd we invert rd and store rd/sd in rd (don't need sd)
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add(t0,x[i],y[i]); _ff_inc(t0); _ff_mult(t1,inv[i],rd[i]); _ff_mult(s[i],t0,t1);
			_ff_square(t2,x[i]); _ff_addto(t0,t2); _ff_mult(r[i],t0,inv[i]);														// s=(x+y+1)/x, r=(x^2+x+y+1)/inv
		}		
		break;
	case 24: // r=(x^2 + x - y + 1)/(x^2 + x*y - y^2 + y), s=(x + 1)/(x + y);
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_add(sd[j],x[i],y[i]); _ff_mult(t0,x[i],sd[j]); _ff_square(t1,y[i]); _ff_subfrom(t0,t1); _ff_add(rd[j],t0,y[i]);					// sd=x+y, rd=x*sd-y^2+y
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add_one(t1,x[i]); _ff_mult(t0,inv[i],rd[i]); _ff_mult(s[i],t0,t1); _ff_subfrom(t1,y[i]); _ff_square(t0,x[i]);					// s = (x+1) / sd
			_ff_addto(t0,t1); _ff_mult(t1,inv[i],sd[i]); _ff_mult(r[i],t0,t1);														// r = (x^2+(x+1)-y) / rd
		}		
		break;
	case 25: // r = (x^2 + xy + y^2 - y)/(x^2 + x + y - 1);,s = (x + y)/(x + 1)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_add_one(sd[j],x[i]); _ff_mult(t0,sd[j],x[i]); _ff_dec(t0); _ff_add(rd[j],t0,y[i]);										// sd = x+1, rd = x(x+1)+y-1
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add(t0,x[i],y[i]); _ff_mult(t1,rd[i],inv[i]); _ff_mult(s[i],t0,t1);													// s = (x+y) / sd
			_ff_mult(t1,t0,x[i]); _ff_square(t0,y[i]); _ff_addto(t1,t0); _ff_subfrom(t1,y[i]); _ff_mult(t0,sd[i],inv[i]); _ff_mult(r[i],t0,t1);		// r = x(x+y)+y^2-y / rd
		}		
		break;
	case 26: // r=(x^3y + 3x^2y - x^2 + xy^2) /((x+1)(x^2y+x^2+3xy+y^2)), s=(xy - x)/(xy + y);
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_mult(t0,x[i],y[i]); _ff_add(sd[j],t0,y[i]); _ff_triple(t1,t0); _ff_square(t0,y[i]); _ff_addto(t1,t0);							// sd = xy+y
			_ff_square(t0,x[i]); _ff_add_one(t2,y[i]); _ff_mult(t3,t0,t2); _ff_addto(t1,t3); _ff_add_one(t0,x[i]); _ff_mult(rd[j],t0,t1);			// rd=(x+1)(x^2y+x^2+3xy+y^2)
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_sub(t0,sd[i],y[i]); _ff_sub(t2,t0,x[i]); _ff_mult(t1,rd[i],inv[i]); _ff_mult(s[i],t1,t2);									// s = (sd-y-x) / sd
			_ff_square(t1,x[i]); _ff_triple(t2,x[i]); _ff_addto(t2,t1); _ff_addto(t2,y[i]);
			_ff_mult(t3,t0,t2); _ff_subfrom(t3,t1); _ff_mult(t2,sd[i],inv[i]); _ff_mult(r[i],t2,t3);										// r = (xy(x^2+3x+y)-x^2) / rd
		}		
		break;
	case 27: // r=(-x^3 - x^2 - x - y) / ((x^2+x-1)y - x), s= (-x^2 - x - y)/(x(y-1) - y)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub_one(t0,y[i]);  _ff_mult(t1,t0,x[i]); _ff_sub(sd[j],t1,y[i]);														// sd = x(y-1)-y
			_ff_square(t0,x[i]); _ff_addto(t0,x[i]); _ff_dec(t0); _ff_mult(t1,t0,y[i]); _ff_sub(rd[j],t1,x[i]);								// rd = (x^2+x-1)y - x
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_square(t2,x[i]); _ff_add(t1,t2,x[i]); _ff_addto(t1,y[i]); _ff_neg(t0,t1); _ff_mult(t1,rd[i],inv[i]); _ff_mult(s[i],t0,t1);			// s = (-x^2 - x - y) / sd
			_ff_mult(t1,t2,x[i]); _ff_subfrom(t0,t1); _ff_mult(t1,sd[i],inv[i]); _ff_mult(r[i],t0,t1);										// r = (-x^3 - x^2 - x - y) / rd
		}		
		break;
	case 28: // r=1+(x*y+y) / ((y-1)*(x*y-x+2*y-1)), s=1-(x*y+y) / ((y-1)*(x-y+1))
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub_one(t0,y[i]);  _ff_sub(t1,x[i],t0); _ff_mult(sd[j],t0,t1);														// sd = (y-1) * (x-(y-1))
			_ff_mult(t1,t0,x[i]); _ff_addto(t1,t0); _ff_addto(t1,y[j]); _ff_mult(rd[j],t0,t1);											// rd = (y-1) * (x(y-1) + (y-1) + y)
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add_one(t1,x[i]); _ff_mult(t0,t1,y[i]);  _ff_mult(t1,rd[i],inv[i]); _ff_mult(t2,t0,t1); _ff_neg(t1,t2); _ff_add_one(s[i],t1);		// s = 1 - (y(x+1)) / sd
			_ff_mult(t1,sd[i],inv[i]);  _ff_mult(t2,t0,t1);  _ff_add_one(r[i],t2);													// r = 1 + (y(x+1)) / rd
		}		
		break;
	case 29: // r=(-x^3 - x^2 - x - y) / ((x^2+x-1)y - x), s= (-x^2 - x - y)/(x(y-1) - y)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub_one(t0,y[i]);  _ff_mult(t1,t0,x[i]); _ff_sub(sd[j],t1,y[i]);														// sd = x(y-1)-y
			_ff_square(t0,x[i]); _ff_addto(t0,x[i]); _ff_dec(t0); _ff_mult(t1,t0,y[i]); _ff_sub(rd[j],t1,x[i]);								// rd = (x^2+x-1)y - x
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_square(t2,x[i]); _ff_add(t1,t2,x[i]); _ff_addto(t1,y[i]); _ff_neg(t0,t1); _ff_mult(t1,rd[i],inv[i]); _ff_mult(s[i],t0,t1);			// s = (-x^2 - x - y) / sd
			_ff_mult(t1,t2,x[i]); _ff_subfrom(t0,t1); _ff_mult(t1,sd[i],inv[i]); _ff_mult(r[i],t0,t1);										// r = (-x^3 - x^2 - x - y) / rd
		}		
		break;
	case 30: // r=1 + y(x+1)/(x^2*y - x*y + x), s= 1 + y(x + 1)/(x^2*y + x)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_mult(t0,x[i],y[i]); _ff_add_one(t1,t0); _ff_mult(sd[j],x[i],t1);  _ff_sub(rd[j],sd[j],t0);									// sd = x(xy+1), rd = sd - xy
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add_one(t1,x[i]); _ff_mult(t0,t1,y[i]);  _ff_mult(t1,rd[i],inv[i]); _ff_mult(t2,t0,t1); _ff_add_one(s[i],t2);					// s = 1 + y(x+1)/sd
			_ff_mult(t1,sd[i],inv[i]);  _ff_mult(t2,t0,t1);  _ff_add_one(r[i],t2);													// r = 1 + y(x+1)/rd
		}		
		break;
	case 31: // r=(-x^3 - x^2 - x - y) / ((x^2+x-1)y - x), s= (-x^2 - x - y)/(x(y-1) - y)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub_one(t0,y[i]);  _ff_mult(t1,t0,x[i]); _ff_sub(sd[j],t1,y[i]);														// sd = x(y-1)-y
			_ff_square(t0,x[i]); _ff_addto(t0,x[i]); _ff_dec(t0); _ff_mult(t1,t0,y[i]); _ff_sub(rd[j],t1,x[i]);								// rd = (x^2+x-1)y - x
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_square(t2,x[i]); _ff_add(t1,t2,x[i]); _ff_addto(t1,y[i]); _ff_neg(t0,t1); _ff_mult(t1,rd[i],inv[i]); _ff_mult(s[i],t0,t1);			// s = (-x^2 - x - y) / sd
			_ff_mult(t1,t2,x[i]); _ff_subfrom(t0,t1); _ff_mult(t1,sd[i],inv[i]); _ff_mult(r[i],t0,t1);										// r = (-x^3 - x^2 - x - y) / rd
		}		
		break;
	case 32: // r=(-x^3 - x^2 - x - y) / ((x^2+x-1)y - x), s= 1-x(x+y)/(x(y-1) - y)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub_one(t0,y[i]);  _ff_mult(t1,t0,x[i]); _ff_sub(sd[j],t1,y[i]);														// sd = x(y-1)-y
			_ff_square(t0,x[i]); _ff_addto(t0,x[i]); _ff_dec(t0); _ff_mult(t1,t0,y[i]); _ff_sub(rd[j],t1,x[i]);								// rd = (x^2+x-1)y - x
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add(t0,x[i],y[i]); _ff_mult(t1,t0,x[i]); _ff_mult(t2,rd[i],inv[i]); _ff_mult(t0,t1,t2);  _ff_neg(t1,t0); _ff_add_one(s[i],t1);			// s = 1- x(x+y) / sd
			_ff_add_one(t0,x[i]); _ff_square(t2,x[i]); _ff_mult(t1,t0,t2); _ff_addto (t1,x[i]); _ff_addto(t1,y[i]); _ff_neg(t0,t1);
			_ff_mult(t1,sd[i],inv[i]); _ff_mult(r[i],t0,t1);																	// r = (-x^3 - x^2 - x - y) / rd
		}		
		break;
	case 33: // r=1+x(y-1)(y(x-1)+1) / ((x-y+1)(xy(x-1)+2x-y+1)), s = 1+ x(y-1) / (x(x-y+2)-y+1)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub(t2,x[i],y[i]);  _ff_inc(t2);  _ff_add_one(t0,t2); _ff_mult(t1,t0,x[i]); _ff_subfrom(t1,y[i]); _ff_add_one(sd[j],t1);			// sd = x(x-y+2)-y+1
			_ff_add(t1,t2,x[i]); _ff_mult(t0,x[i],y[i]); _ff_sub_one(t3,x[i]); _ff_mult(t4,t0,t3); _ff_addto(t4,t1);  _ff_mult(rd[j],t2,t4);			// rd = (x-y+1)(xy(x-1)+2x-y+1)
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_sub_one(t0,y[i]); _ff_mult(t1,t0,x[i]); _ff_mult(t2,rd[i],inv[i]); _ff_mult(t0,t1,t2);  _ff_add_one(s[i],t0);					// s = 1+x(y-1) / sd
			_ff_sub_one(t0,x[i]); _ff_mult(t2,t0,y[i]); _ff_inc(t2); _ff_mult(t0,t1,t2); _ff_mult(t1,sd[i],inv[i]); _ff_mult(t2,t0,t1);
			_ff_add_one(r[i],t2);																						// r = 1+x(y-1)(y(x-1)+1)/ rd
		}		
		break;
	case 34: // r=(-x^3 - x^2 - x - y) / ((x^2+x-1)y - x), s= 1-x(x+y)/(x(y-1) - y)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub_one(t0,y[i]);  _ff_mult(t1,t0,x[i]); _ff_sub(sd[j],t1,y[i]);														// sd = x(y-1)-y
			_ff_square(t0,x[i]); _ff_addto(t0,x[i]); _ff_dec(t0); _ff_mult(t1,t0,y[i]); _ff_sub(rd[j],t1,x[i]);								// rd = (x^2+x-1)y - x
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add(t0,x[i],y[i]); _ff_mult(t1,t0,x[i]); _ff_mult(t2,rd[i],inv[i]); _ff_mult(t0,t1,t2);  _ff_neg(t1,t0); _ff_add_one(s[i],t1);			// s = 1- x(x+y) / sd
			_ff_add_one(t0,x[i]); _ff_square(t2,x[i]); _ff_mult(t1,t0,t2); _ff_addto (t1,x[i]); _ff_addto(t1,y[i]); _ff_neg(t0,t1);
			_ff_mult(t1,sd[i],inv[i]); _ff_mult(r[i],t0,t1);																	// r = (-x^3 - x^2 - x - y) / rd
		}		
		break;
	case 35: // r=1+x(y+1)(x+y) / ((x-1)(x(x-y-2)-y(y+1))), s = 1+ x(y+1) / (x(x-y-2)+y+1)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub(t0,x[i],y[i]);  _ff_subfrom(t0,c_2); _ff_mult(t2,t0,x[i]); _ff_add_one(t1,y[i]); _ff_add(sd[j],t1,t2);						// sd = x(x-y-2)+y+1
			_ff_mult(t0,t1,y[i]);  _ff_subfrom(t2,t0); _ff_sub_one(t0,x[i]); _ff_mult(rd[j],t0,t2);										// rd = (x-1)(x(x-y-2)-y(y+1))
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add_one(t0,y[i]); _ff_mult(t1,t0,x[i]); _ff_mult(t2,rd[i],inv[i]); _ff_mult(t0,t1,t2);  _ff_add_one(s[i],t0);					// s = 1+x(y+1) / sd
			_ff_add(t0,x[i],y[i]); _ff_mult(t2,t0,t1);  _ff_mult(t1,sd[i],inv[i]); _ff_mult(t0,t1,t2);  _ff_add_one(r[i],t0);						// r = 1+x(y+1)(x+y)/ rd
		}		
		break;
	case 36: // r=(-x^3 - x^2 - x - y) / ((x^2+x-1)y - x), s= -(x(x+1)+y)/(x(y-1) - y)
		for ( i = j=0 ; i < n ; i++,j++ ) {
			_ff_sub_one(t0,y[i]);  _ff_mult(t1,t0,x[i]); _ff_sub(sd[j],t1,y[i]);														// sd = x(y-1)-y
			_ff_square(t0,x[i]); _ff_addto(t0,x[i]); _ff_dec(t0); _ff_mult(t1,t0,y[i]); _ff_sub(rd[j],t1,x[i]);								// rd = (x^2+x-1)y - x
			_ff_mult(den[j],rd[j],sd[j]);
			if ( _ff_zero(den[j]) ) {j--; continue; }																		// skip zero denominators (should be very rare)
			if ( i == j ) continue;
			_ff_set(x[j],x[i]); _ff_set(y[j],y[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		_ff_set_one(t1);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_add_one(t0,x[i]); _ff_mult(t1,t0,x[i]); _ff_add(t2,t1,y[i]); _ff_neg(t0,t2);  _ff_mult(t2,rd[i],inv[i]); _ff_mult(s[i],t0,t2);			// s = - x(x+y) / sd
			_ff_mult(t0,t1,x[i]); _ff_addto(t0,x[i]); _ff_addto(t0,y[i]); _ff_neg(t1,t0); _ff_mult(t2,sd[i],inv[i]); _ff_mult(r[i],t1,t2);			// r = (-x^3 - x^2 - x - y) / rd
		}		
		break;
	default:
		printf ("Invalid N=%d in tecurve_get_X1_points\n", N); exit (0);
	}
	return n;
}


int tecurve_get_tt_points (ff_t tt_x[], ff_t r[], ff_t s[], int n, int N)
{
	ff_t den[TECURVE_MAX_CURVES], inv[TECURVE_MAX_CURVES];
	register ff_t t0,t1,t2,t3,t4;
	register int i, j;
	
	switch (N) {
	case 14: // P7x = rs(r - 1)(s - 1)(rs - 2r + 1) / (r - s)^2
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_sub(t0,r[i],s[i]);
			if ( _ff_zero(t0) ) { j--; continue; }
			_ff_square(den[j],t0);
			if ( i == j ) continue;
			_ff_set(r[j],r[i]); _ff_set(s[j],s[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_sub_one(t0,r[i]); _ff_sub_one(t1,s[i]); _ff_mult(t2,r[i],s[i]); _ff_sub(t3,t2,t0); _ff_subfrom(t3,r[i]);
			_ff_mult(t4,t0,t1); _ff_mult(t0,t2,t3); _ff_mult(t1,inv[i],t4); _ff_mult(tt_x[i],t0,t1);
		}
		break;
	case 16: // P8x =  r(r - 1)(r - s)(r - s^2 + s - 1) / (rs - 2r + 1)^2
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_sub(t0,s[i],c_2); _ff_mult(t1,t0,r[i]); _ff_inc(t1);
			if ( _ff_zero(t1) ) { j--; continue; }
			_ff_square(den[j],t1);
			if ( i == j ) continue;
			_ff_set(r[j],r[i]); _ff_set(s[j],s[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_sub_one(t0,r[i]); _ff_sub(t1,r[i],s[i]); _ff_square(t3,s[i]); _ff_sub(t2,r[i],t3); _ff_addto(t2,s[i]); _ff_dec(t2);
			_ff_mult(t3,r[i],t0); _ff_mult(t0,t1,t2); _ff_mult(t1,inv[i],t3); _ff_mult(tt_x[i],t0,t1);
		}
		break;
	case 18: // P9x = s(r - 1)(rs - 2r + 1)(rs^2 - 3rs + r + s^2 ) / (r - s^2 + s - 1)^2
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_square(t1,s[i]); _ff_sub(t0,r[i],t1); _ff_addto(t0,s[i]); _ff_dec(t0);
			if ( _ff_zero(t0) ) { j--; continue; }
			_ff_square(den[j],t0);	
			if ( i == j ) continue;
			_ff_set(r[j],r[i]); _ff_set(s[j],s[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_square(t2,s[i]); _ff_triple(t0,s[i]); _ff_sub(t1,t2,t0); _ff_inc(t1); _ff_mult(t0,r[i],t1); _ff_addto(t0,t2);
			_ff_sub(t2,s[i],c_2); _ff_mult(t1,t2,r[i]); _ff_inc(t1); _ff_sub_one(t2,r[i]);
			_ff_mult(t3,s[i],t2); _ff_mult(t2,t0,t1); _ff_mult(t1,inv[i],t3); _ff_mult(tt_x[i],t1,t2);
		}
		break;
	case 20: // P10x = rs(r - s^2 + s - 1)(r^2 - rs^3 + 3rs^2 - 4rs + s) / (rs^2 - 3rs + r + s^2 )^2
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_square(t2,s[i]); _ff_triple(t0,s[i]); _ff_sub(t1,t2,t0); _ff_inc(t1); _ff_mult(t0,r[i],t1); _ff_addto(t0,t2);
			if ( _ff_zero(t0) ) { j--; continue; }
			_ff_square(den[j],t0);	
			if ( i == j ) continue;
			_ff_set(r[j],r[i]); _ff_set(s[j],s[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_square(t3,s[i]); _ff_sub(t0,r[i],t3); _ff_addto(t0,s[i]); _ff_dec(t0);
			_ff_triple(t1,s[i]); _ff_subfrom(t1,t3); _ff_subfrom(t1,c_4); _ff_mult(t2,t1,s[i]); _ff_addto(t2,r[i]); _ff_mult(t1,t2,r[i]); _ff_addto(t1,s[i]);
			_ff_mult(t2,t0,t1); _ff_mult(t3,r[i],s[i]); _ff_mult(t0,t2,t3); _ff_mult(tt_x[i],t0,inv[i]);
		}
		break;
	case 22: // P11x = rs(r-1)(s-1)(rs^2-3rs+r+s^2)(r^2s-3r^2+rs+3r-s^2-1) /  (r^2 - rs^3 + 3rs^2 - 4rs + s)^2
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_square(t0,s[i]); _ff_sub(t1,s[i],t0); _ff_sub(t0,r[i],t1); _ff_dec(t0); _ff_addto(t1,s[i]); _ff_addto(t1,s[i]); _ff_subfrom(t1,c_4);
			_ff_mult(t2,t1,s[i]); _ff_addto(t2,r[i]); _ff_mult(t1,t2,r[i]); _ff_addto(t1,s[i]); 
			if ( _ff_zero(t1) ) { j--; continue; }
			_ff_square(den[j],t1);	
			if ( i == j ) continue;
			_ff_set(r[j],r[i]); _ff_set(s[j],s[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_square(t3,s[i]); _ff_triple(t0,s[i]); _ff_sub(t1,t3,t0); _ff_inc(t1); _ff_mult(t0,r[i],t1); _ff_addto(t0,t3);
			_ff_sub(t1,s[i],c_3); _ff_mult(t2,t1,r[i]); _ff_addto(t2,s[i]); _ff_addto(t2,c_3); _ff_mult(t1,t2,r[i]); _ff_subfrom(t1,t3); _ff_dec(t1);
			_ff_sub_one(t2,r[i]); _ff_sub_one(t3,s[i]); _ff_mult(t4,t2,t3); _ff_mult(t2,t0,t1); _ff_mult(t1,r[i],s[i]); _ff_mult(t3,t4,inv[i]);
			_ff_mult(t0,t2,t3); _ff_mult(tt_x[i],t0,t1);
		}
		break;
	case 24: // P12x = (r-1)(r^2 - rs^3 + 3rs^2 - 4rs + s)(r^3 - r^2s^4 + 5r^2s^3 - 9r^2s^2 + 4r^2s - 2r^2 - rs^3 + 6rs^2 - 3rs + r - s^3) / ((s-1)(r^2s - 3r^2 + rs + 3r - s^2 - 1))^2
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_sub(t1,s[i],c_3); _ff_mult(t2,t1,r[i]); _ff_addto(t2,s[i]); _ff_addto(t2,c_3); _ff_mult(t1,t2,r[i]); _ff_square(t2,s[i]);
			_ff_subfrom(t1,t2); _ff_dec(t1); _ff_sub_one(t0,s[i]); _ff_mult(t2,t0,t1);
			if ( _ff_zero(t2) ) { j--; continue; }
			_ff_square(den[j],t2);
			if ( i == j ) continue;
			_ff_set(r[j],r[i]); _ff_set(s[j],s[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_sub(t0,c_4,s[i]); _ff_inc(t0);  _ff_mult(t1,t0,s[i]); _ff_subfrom(t1,c_9); _ff_mult(t0,r[i],t1); _ff_subfrom(t0,s[i]); _ff_addto(t0,c_6); 
			_ff_mult(t1,t0,s[i]); _ff_add(t2,r[i],r[i]); _ff_x2(t2); _ff_addto(t1,t2); _ff_subfrom(t1,c_3); _ff_mult(t0,t1,s[i]); _ff_sub_one(t1,r[i]); 
			_ff_square(t2,t1); _ff_addto(t0,t2); _ff_mult(t1,t0,r[i]); _ff_square(t2,s[i]); _ff_mult(t0,t2,s[i]); _ff_sub(t3,t1,t0);									// t3 = r*(s*(s*(r*(s*(5-s)-9)-s+6)+4*r-3)+(r-1)^2)-s^3
			_ff_square(t0,s[i]); _ff_sub(t1,s[i],t0); _ff_sub(t0,r[i],t1); _ff_dec(t0); _ff_addto(t1,s[i]); _ff_addto(t1,s[i]); _ff_subfrom(t1,c_4);
			_ff_mult(t2,t1,s[i]); _ff_addto(t2,r[i]); _ff_mult(t1,t2,r[i]); _ff_addto(t1,s[i]); _ff_sub_one(t0,r[i]); _ff_mult(t2,t0,t1);
			_ff_mult(t0,t2,t3); _ff_mult(tt_x[i],t0,inv[i]);			
		}
		break;
	case 26: // P13x = rs(r-1)(s-1)(r-s)(r^2s - 3r^2 + rs + 3r - s^2 - 1)(r^2s^3 - 5r^2s^2 + 6r^2s - r^2 + rs^4 - 3rs^3 + 6rs^2 - 7rs + r + s)
			//              /  (r^3 - r^2s^4 + 5r^2s^3 - 9r^2s^2 + 4r^2s - 2r^2 - rs^3 + 6rs^2 - 3rs + r - s^3)^2
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_sub(t0,c_4,s[i]); _ff_inc(t0);  _ff_mult(t1,t0,s[i]); _ff_subfrom(t1,c_9); _ff_mult(t0,r[i],t1); _ff_subfrom(t0,s[i]); _ff_addto(t0,c_6); 
			_ff_mult(t1,t0,s[i]); _ff_add(t2,r[i],r[i]); _ff_x2(t2); _ff_addto(t1,t2); _ff_subfrom(t1,c_3); _ff_mult(t0,t1,s[i]); _ff_sub_one(t1,r[i]); 
			_ff_square(t2,t1); _ff_addto(t0,t2); _ff_mult(t1,t0,r[i]); _ff_square(t2,s[i]); _ff_mult(t0,t2,s[i]); _ff_sub(t2,t1,t0);									// t2 = r*(s*(s*(r*(s*(5-s)-9)-s+6)+4*r-3)+(r-1)^2)-s^3
			if ( _ff_zero(t2) ) { j--; continue; }
			_ff_square(den[j],t2);
			if ( i == j ) continue;
			_ff_set(r[j],r[i]); _ff_set(s[j],s[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_sub(t1,s[i],c_3); _ff_mult(t2,t1,r[i]); _ff_addto(t2,s[i]); _ff_addto(t2,c_3); _ff_mult(t1,t2,r[i]); _ff_square(t2,s[i]); _ff_subfrom(t1,t2); _ff_dec(t1); 
			_ff_sub(t0,s[i],c_4); _ff_dec(t0); _ff_mult(t3,t0,s[i]); _ff_addto(t3,c_6); _ff_mult(t0,t3,r[i]); _ff_mult(t3,t2,s[i]); _ff_addto(t0,t3); _ff_add(t4,s[i],s[i]);
			_ff_subfrom(t4,t2); _ff_triple(t3,t4); _ff_addto(t0,t3); _ff_subfrom(t0,c_6); _ff_dec(t0); _ff_mult(t3,t0,s[i]); _ff_subfrom(t3,r[i]); _ff_inc(t3);
			_ff_mult(t0,t3,r[i]); _ff_addto(t0,s[i]);																							// t0 = r*(s*(r*(s*(s-5)+6)+s^3-3*s^2+6*s-7)-r+1)+s
			_ff_mult(t2,t0,t1); _ff_sub(t0,r[i],s[i]); _ff_sub_one(t1,r[i]); _ff_mult(t3,t0,t1); _ff_mult(t0,t2,t3);
			_ff_sub_one(t1,s[i]); _ff_mult(t2,t0,t1); _ff_mult(t0,r[i],s[i]); _ff_mult(t1,t0,t2);
			_ff_mult(tt_x[i],t1,inv[i]);	
		}
		break;
	case 28: // P14x  = rs(r-1)(r^3 - r^2*s^4 + 5*r^2*s^3 - 9*r^2*s^2 + 4*r^2*s - 2*r^2 - r*s^3 + 6*r*s^2 - 3*r*s + r - s^3)
	               //                  (r^3 - r^2*s^5 + 7*r^2*s^4 - 18*r^2*s^3 + 19*r^2*s^2 - 10*r^2*s - r*s^5 + 4*r*s^4 - 5*r*s^2 + 5*r*s - s^5 + s^4 - s^3 + s^2 - s)
		       //               /  (r - s)^2 (r^2*s^3 - 5*r^2*s^2 + 6*r^2*s - r^2 + r*s^4 - 3*r*s^3 + 6*r*s^2 - 7*r*s + r + s)^2
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_sub(t1,s[i],c_3); _ff_mult(t0,t1,s[i]); _ff_addto(t0,c_6); _ff_mult(t1,t0,s[i]); _ff_subfrom(t1,c_4); _ff_subfrom(t1,c_3); _ff_mult(t0,t1,s[i]); 
			_ff_sub(t1,t0,r[i]);	_ff_inc(t1);																								// t1 = -r+s(s(s(s-3)+6)-7)+1
			_ff_sub(t0,s[i],c_4); _ff_dec(t0); _ff_mult(t3,t0,s[i]); _ff_addto(t3,c_6); _ff_mult(t0,t3,r[i]); _ff_mult(t3,t0,s[i]); 									// t3 = s(r(s(s-5)+6))
			_ff_addto(t3,t1); _ff_mult(t1,t3,r[i]);  _ff_addto(t1,s[i]); _ff_sub(t0,r[i],s[i]); _ff_mult(t2,t0,t1);  												// t1 = r*(s*(r*(s*(s-5)+6))-r+s(s(s(s-3)+6)-7)+1)+s
			if ( _ff_zero(t2) ) { j--; continue; }
			_ff_square(den[j],t2);
			if ( i == j ) continue;
			_ff_set(r[j],r[i]); _ff_set(s[j],s[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_sub(t0,c_4,s[i]); _ff_inc(t0);  _ff_mult(t1,t0,s[i]); _ff_subfrom(t1,c_9); _ff_mult(t0,r[i],t1); _ff_subfrom(t0,s[i]); _ff_addto(t0,c_6); 
			_ff_mult(t1,t0,s[i]); _ff_add(t2,r[i],r[i]); _ff_x2(t2); _ff_addto(t1,t2); _ff_subfrom(t1,c_3); _ff_mult(t0,t1,s[i]); _ff_sub_one(t1,r[i]); 
			_ff_square(t2,t1); _ff_addto(t0,t2); _ff_mult(t1,t0,r[i]); _ff_square(t2,s[i]); _ff_mult(t0,t2,s[i]); _ff_sub(t2,t1,t0);									// t2 = r*(s*(s*(r*(s*(5-s)-9)-s+6)+4*r-3)+(r-1)^2)-s^3
			_ff_sub(t0,s[i],c_6); _ff_dec(t0); _ff_mult(t1,t0,s[i]);  _ff_add(t4,c_9,c_9); _ff_addto(t1,t4);  _ff_mult(t0,t1,s[i]);
			_ff_inc(t4); _ff_subfrom(t0,t4);	_ff_mult(t1,t0,s[i]); _ff_addto(t1,c_8); _ff_addto(t1,c_2);  _ff_mult(t0,t1,s[i]);_ff_sub(t1,r[i],t0); _ff_mult(t3,r[i],t1);
			_ff_sub(t0,s[i],c_4); _ff_mult(t1,t0,s[i]); _ff_mult(t0,t1,s[i]);
			_ff_add_one(t4,c_4); _ff_addto(t0,t4); _ff_mult(t1,t0,s[i]); _ff_subfrom(t1,t4); _ff_mult(t0,t1,s[i]); _ff_subfrom(t3,t0); _ff_mult(t4,t3,r[i]);				// t4 = r(r(r-s(s(s(s(s-7)+18)-19)+10))-s(s(s(s(s-4))+5)-5))
			_ff_sub_one(t0,s[i]); _ff_mult(t1,t0,s[i]); _ff_inc(t1); _ff_mult(t0,t1,s[i]); _ff_dec(t0); _ff_mult(t1,t0,s[i]); _ff_inc(t1); _ff_mult(t0,t1,s[i]);	_ff_subfrom(t4,t0);	// t4 -= s(s(s(s(s-1)+1)-1)+1)
			_ff_sub_one(t0,r[i]); _ff_mult(t1,t0,r[i]); _ff_mult(t0,t1,s[i]); _ff_mult(t1,t0,t2); _ff_mult(t0,t1,t4);
			_ff_mult(tt_x[i],t0,inv[i]);
		}
		break;
	case 30: // P15x  = s(r-s)(rs-2r+1)(r^2*s^3 - 5*r^2*s^2 + 6*r^2*s - r^2 + r*s^4 - 3*r*s^3 + 6*r*s^2 - 7*r*s + r + s)
		       //              (r^3*s^2 - 4*r^3*s + 2*r^3 + 3*r^2*s^2 + 2*r^2*s - 2*r^2 - r*s^5 + 4*r*s^4 - 10*r*s^3 + 6*r*s^2 - 3*r*s + r + s^4)
	               //               /   (r^3 - r^2*s^5 + 7*r^2*s^4 - 18*r^2*s^3 + 19*r^2*s^2 - 10*r^2*s - r*s^5 + 4*r*s^4 - 5*r*s^2 + 5*r*s - s^5 + s^4 - s^3 + s^2 - s)^2
		for ( i = j = 0 ; i < n ; i++,j++ ) {
			_ff_sub(t0,s[i],c_6); _ff_dec(t0); _ff_mult(t1,t0,s[i]);  _ff_add(t4,c_9,c_9); _ff_addto(t1,t4);  _ff_mult(t0,t1,s[i]);
			_ff_inc(t4); _ff_subfrom(t0,t4);	_ff_mult(t1,t0,s[i]); _ff_addto(t1,c_8); _ff_addto(t1,c_2);  _ff_mult(t0,t1,s[i]);_ff_sub(t1,r[i],t0); _ff_mult(t3,r[i],t1);
			_ff_sub(t0,s[i],c_4); _ff_mult(t1,t0,s[i]); _ff_mult(t0,t1,s[i]);
			_ff_add_one(t4,c_4); _ff_addto(t0,t4); _ff_mult(t1,t0,s[i]); _ff_subfrom(t1,t4); _ff_mult(t0,t1,s[i]); _ff_subfrom(t3,t0); _ff_mult(t4,t3,r[i]);				// t4 = r(r(r-s(s(s(s(s-7)+18)-19)+10))-s(s(s(s(s-4))+5)-5))
			_ff_sub_one(t0,s[i]); _ff_mult(t1,t0,s[i]); _ff_inc(t1); _ff_mult(t0,t1,s[i]); _ff_dec(t0); _ff_mult(t1,t0,s[i]); _ff_inc(t1); _ff_mult(t0,t1,s[i]);	_ff_subfrom(t4,t0);	// t4 -= s(s(s(s(s-1)+1)-1)+1)
			if ( _ff_zero(t4) ) { j--; continue; }
			_ff_square(den[j],t4);
			if ( i == j ) continue;
			_ff_set(r[j],r[i]); _ff_set(s[j],s[i]);
		}
		n = j;
		ff_parallel_invert(inv,den,n);
		for ( i = 0 ; i < n ; i++ ) {
			_ff_sub(t0,r[i],s[i]); _ff_mult(t1,t0,s[i]);  _ff_sub(t0,s[i],c_2); _ff_mult(t2,t0,r[i]); _ff_inc(t2); _ff_mult(t4,t1,t2);									// t4 = s(r-s)(rs-2r+1)
			_ff_sub(t1,s[i],c_3); _ff_mult(t0,t1,s[i]); _ff_addto(t0,c_6); _ff_mult(t1,t0,s[i]); _ff_subfrom(t1,c_4); _ff_subfrom(t1,c_3); _ff_mult(t0,t1,s[i]); 
			_ff_sub(t1,t0,r[i]);	_ff_inc(t1);																								// t1 = -r+s(s(s(s-3)+6)-7)+1
			_ff_sub(t0,s[i],c_4); _ff_dec(t0); _ff_mult(t3,t0,s[i]); _ff_addto(t3,c_6); _ff_mult(t0,t3,r[i]); _ff_mult(t3,t0,s[i]); 									// t3 = s(r(s(s-5)+6))
			_ff_addto(t3,t1); _ff_mult(t1,t3,r[i]);  _ff_addto(t1,s[i]); 	_ff_mult(t3,t1,t4);																// t1 = r*(s*(r*(s*(s-5)+6))-r+s(s(s(s-3)+6)-7)+1)+s
			_ff_sub(t0,s[i],c_4) _ff_mult(t2,t0,s[i]); _ff_add(t1,t2,c_8); _ff_addto(t1,c_2); _ff_mult(t0,t1,s[i]); _ff_subfrom(t0,c_6); _ff_mult(t1,t0,s[i]);
			_ff_addto(t1,c_3); _ff_mult(t0,t1,s[i]); _ff_sub_one(t4,t0);
			_ff_triple(t0,s[i]); _ff_addto(t0,c_2); _ff_mult(t1,t0,s[i]); _ff_sub(t2,t1,c_2);
			_ff_sub(t0,s[i],c_4); _ff_mult(t1,t0,s[i]); _ff_addto(t1,c_2); _ff_mult(t0,t1,r[i]); _ff_addto (t0,t2);
			_ff_mult(t1,t0,r[i]);  _ff_subfrom(t1,t4); _ff_mult(t0,t1,r[i]);
			_ff_square(t2,s[i]); _ff_square(t4,t2); _ff_addto(t0,t4); _ff_mult(t1,t0,t3);																	// t0 = r(r(r(s(s-4)+2)+s(3s+2)-2)-(s(s(s(s(s-4)+10)-6)+3)-1)+s^4
			_ff_mult(tt_x[i],t1,inv[i]);
		}
		break;
	default:
		printf ("Invalid N=%d in tecurve_get_tt_points\n", N); exit (0);
	}
	return n;
}

int tecurve_filter (ff_t A[], ff_t B[], ff_t *tt_x, ff_t *tor_x, ff_t *base_x, ff_t *base_y, int n, int N, int s2_flag, int t3_flag, int j_flags)
{
//ff_t J[TECURVE_MAX_CURVES];
	ff_t D, f[4],r[2],x,y;
	register ff_t t0,t1,tt;
	register int i, j, k, m, M;
	int r1, min_flag;
	
//	if ( N > 36 ) { printf ("N=%ld, y^2 = x^3 + %ld*x + %ld\n", N, _ff_get_ui(A[0]), _ff_get_ui(B[0])); if ( tt_x ) printf("two-torsion point (%ld,0)\n", _ff_get_ui(tt_x[0])); }
	
	min_flag = 0;
	if ( (s2_flag&FILTER_2SYLOW_C_LG) ) {
		if ( s2_flag&FILTER_2SYLOW_C_LG_MIN ) min_flag = 1;		// min_flag means we require a point with order divisible by 2^k, but not exactly 2^k
		m = s2_flag&(FILTER_2SYLOW_C_LG-1);
		if ( min_flag ) {
			while ( !(N&((1<<(m+1))-1)) ) m++;
		} else {
			if ( !(N&((1<<(m+1))-1)) ) { printf ("Unsatisfiable FILTER_2SYLOW_C_LG condition N=%d, m=%d\n", N, m); exit (0); }
		}
	} else {
		m = -1;
	}
	if ( t3_flag == FILTER_3TOR_NOT_1 && !(N%3) ) t3_flag = 0;
	if ( s2_flag == FILTER_2SYLOW_NOT_1 && !(N&1) ) s2_flag = 0;
/*
for ( i = 0 ; i < n ; i++ ) {
tecurve_discriminant_AB (&D,A[i],B[i]);  if ( _ff_zero(D) ) { J[i] = _ff_p; continue; }
_ff_set(f[1],A[i]); _ff_set(f[0],B[i]);
ff_poly_g1_to_jinv(J+i, f);
//printf ("%d) j=%ld A=%ld B=%ld\n", i, _ff_get_ui(J[i]), _ff_get_ui(A[i]), _ff_get_ui(B[i]));
}
*/

//	printf("N=%d, s2_flag=%d, t3_flag=%d, m=%d, tt_x=%d, tor_x=%d\n", N, s2_flag, t3_flag, m, (tt_x?1:0), (tor_x?1:0));
	_ff_set_one(f[3]); _ff_set_zero(f[2]);
	for ( i = j = 0 ; i < n ; i++, j++ ) {
		tecurve_discriminant_AB (&D,A[i],B[i]);
		if ( _ff_zero(D) ) goto skip;
/*
_ff_set(f[1],A[i]); _ff_set(f[0],B[i]);
ff_poly_g1_to_jinv(r, f);
for ( k = 0 ; k < n ; k++ ) if ( k != i && _ff_equal(r[0],J[k]) ) break;
if ( k < n  ) printf ("Duplicate j invariant %ld over p=%ld (i=%d,k=%d) A=%ld, B=%ld for N=%d\n", _ff_get_ui(r[0]), _ff_p, i, k, _ff_get_ui(A[i]), _ff_get_ui(B[i]),  N); 
*/
		if ( _ff_p1mod3 && (j_flags & FILTER_JCUBE) ) {
			_ff_mult(x,D,d_4);  if ( ! ff_cbrt(&y,&x) ) goto skip;															// j is a perfect cube iff (-12f1)^3 / j = D/4 is a perfect cube
		}
/*		if ( j_flags & FILTER_WEBER ) {
			_ff_set(f[1],A[i]); _ff_set(f[0],B[i]);
			ff_poly_g1_to_jinv(&x,f);  _ff_set_ui(t0,1728);  _ff_sub(y,x,t0);
printf ("p=%ld, A=%ld, B=%ld, j=%ld, Weber test=%d\n", _ff_p, _ff_get_ui(A[i]), _ff_get_ui(B[i]), _ff_get_ui(x), ff_residue(y));
		}
*/		if ( !min_flag && m==0 ) { if ( ! ff_residue(D) ) goto skip; else goto next; }											// don't completely filter m=0 case, let rank 2 through
		if ( m > 0 ) {
			_ff_set(f[1],A[i]); _ff_set(f[0],B[i]);
			if ( N&1 ) {
				k = ff_poly_roots_d3(r,f);
				if ( k > 1 || ( m&& !k) ) goto skip;
				_ff_set(tt,r[0]);
				tecurve_AB_eval_df(&y,tt,A[i]);
				if ( ! ff_sqrt(&x,&y) ) { if ( m== 1 ) goto next; else goto skip; }
				if ( m==1 ) goto skip;
				_ff_add(t0,tt,x);
				tecurve_AB_eval_f(&y,t0,A[i],B[i]);
				if ( ! ff_residue(y) ) _ff_sub(t0,tt,x);																	// if tt+x is not the x-ccord of 4-torsion point, tt-x must be
				_ff_set(x,t0);
				k = 2;
			} else if ( (N&3)==2 ) {
				r1 = ( N==2 ? 1 : 0 );																			// r1 = 1if we know  that the 2-Sylow has rank 1, 0 otherwise
				if ( (_ff_p&3) == 1 && !r1 ) { if ( ff_residue(D) ) goto skip; else r1=1; }									// when p=1mod4 check rank of 2-Sylow first
				if ( tt_x ) {
					_ff_set(tt,tt_x[i]);
				} else {
					if ( ! ff_poly_roots_d3(r,f) ) { printf ("Fatal error in tecurve, no 2-torsion for N=%d\n", N); exit(0); }
					_ff_set(tt,r[0]);
				}
				tecurve_AB_eval_df(&y,tt,A[i]);
				k = ff_sqrt(&x,&y);
				if ( (m==1 && k) || (m>1 && !k) ) goto skip;
				if ( !r1 ) if ( ff_residue(D) ) goto skip;																// if we haven't confirmed rank of 2-Sylow is 1, do so now
				if ( m==1 ) goto next;
				_ff_add(t0,tt,x);
				tecurve_AB_eval_f(&y,t0,A[i],B[i]);
				if ( ! ff_residue(y) ) _ff_sub(t0,tt,x);																	// if tt+x is not the x-ccord of 4-torsion point, tt-x must be
				_ff_set(x,t0);
				k = 2;
			} else {
				if ( ff_residue(D) ) goto skip;
				for ( k = 2, M=N>>2 ; !(M&1) ; k++, M>>=1 );
				if ( k > m ) goto skip;
				if ( tt_x ) {
					_ff_set(tt,tt_x[i]);
				} else {
					if ( ! ff_poly_roots_d3(r,f) ) { printf ("Fatal error in tecurve, no 2-torsion for N=%d\n", N); exit(0); }
					_ff_set(tt,r[0]);
				}
				_ff_set(x,tor_x[i]);
			}
			// if we get here we know that m >=k>=2 and x is the x-coord of a point of order 2^k
			if ( ! ecurve_verify_2depth (x,m-k,min_flag,tt,f) ) goto skip;
//			for ( ; k < m ; k++ ) if ( ! ff_poly_g1_halve_x (&x,x,f) ) goto skip;
//			if ( ff_poly_g1_halve_x (&x,x,f) ) goto skip;
		} else {
			switch ( s2_flag ) {
			case FILTER_2SYLOW_1:
				if ( ! ff_residue(D) ) goto skip;																	// if discriminant is not a residue, then g(x) has 2 factors (Stickleberger) and y^2=f(x) has 2-torsion
				break;																					// we could still have 2-torsion rank 2 here (g(x) with 3 factors), but we don't filter this
			case FILTER_2SYLOW_1MOD4:
				if ( (_ff_p&3)==3 ) { if ( ! ff_residue(D) ) goto skip; else break; }										// if p is 3 mod 4 treat as FILTER_2SYLOW_1
				ff_exp_ui(r,&D,(_ff_p-1)>>2);																	// if p is 1 mod 4, then apply a simplified version of the test of Nogami-Morikawa (ICISC 2004)
				_ff_neg(t0,r[0]);
				if ( ! _ff_one(t0) ) goto skip;
				break;
			case FILTER_2SYLOW_3MOD4:
				if ( (_ff_p&3)==3 ) { if ( ! ff_residue(D) ) goto skip; else break; }										// if p is 3 mod 4 treat as FILTER_2SYLOW_1
				ff_exp_ui(r,&D,(_ff_p-1)>>2);																	// if p is 1 mod 4, apply apply a simplified version of the test of Nogami-Morikawa (ICISC 2004)
				if ( ! _ff_one(r[0]) ) goto skip;
				break;
			case FILTER_2SYLOW_2:																			// note that tor=2 case already ensures 2-torsion rank is 1
				if ( ! tt_x ) {
					_ff_set(f[1],A[i]); _ff_set(f[0],B[i]);
					if ( ff_poly_roots_d3(r,f) != 1 ) goto skip;														// this is costly, much cheaper if 2-torsion point tt is known by construction
					tecurve_AB_eval_df(&y,r[0],A[i]);
					if ( ff_residue(y) ) { j--; continue;  };
				} else {
					if ( (_ff_p&3) == 1 ) if ( N > 2 && ff_residue(D) ) goto skip;
					tecurve_AB_eval_df(&y,tt_x[i],A[i]);
					if ( ff_residue(y) ) goto skip;
					if ( (_ff_p&3) == 3 ) if ( N > 2 && ff_residue(D) ) goto skip;
				}
				break;
			case FILTER_2SYLOW_NOT_1:																		// currently we force rank 1 to achive non-trivial 2-torsion for odd torsion case
				if ( !(N&1) ) break;
				// fall through
			case FILTER_2SYLOW_R1:
				if ( n > 2 && ff_residue(D) ) goto skip;															// if discriminant is a residue than we have 2-torsion rank 0 or 2 (but not 1)		
				break;
			case FILTER_2SYLOW_R14:
				if ( !(N&3) ) { if ( ff_residue(D) ) goto skip; else break; }											// if 4-torsion is already present, just verify rank
				if ( ! tt_x ) {
					_ff_set(f[1],A[i]); _ff_set(f[0],B[i]);
					if ( ff_poly_roots_d3(r,f) != 1 ) goto skip;														// this is costly, much cheaper if 2-torsion point tt is known by construction
					tecurve_AB_eval_df(&y,r[0],A[i]);
					if ( ! ff_residue(y) ) goto skip;																// if t0 is not a residue, 4-torsion is not possible  (o.w. 4-torsion must happen in rank 1 case)
				} else {
					if ( n > 2 && ff_residue(D) ) goto skip;														// if discriminant is a residue than we have 2-torsion rank 0 or 2 (but not 1)
					tecurve_AB_eval_df(&y,tt_x[i],A[i]);
					if ( ! ff_residue(y) ) goto skip;																// if t0 is not a residue, 4-torsion is not possible  (o.w. 4-torsion must happen in rank 1 case)
				}
				break;
			case FILTER_2SYLOW_R2:
				if ( ! _ff_zero(tt_x) ) {
					_ff_square(t0,tt_x[i]); _ff_add(t1,t0,A[i]); _ff_x2(t1); _ff_x2(t1); _ff_subfrom(t0,t1);						// compute discriminant of g(x)/(x-tt) (we know g(tt)=0 since tt is the x-coord of 2-torsion pt)
					if ( ! ff_residue(t0) ) goto skip;
				} else {
					_ff_set(f[1],A[i]); _ff_set(f[0],B[i]);
					if ( ff_poly_roots_d3(0,f) < 3 ) goto skip;
				}
				break;
			case FILTER_2SYLOW_D4:
				if ( !(N&3) ) break;
				if ( ! tt_x ) {
					_ff_set(f[1],A[i]); _ff_set(f[0],B[i]);
					k = ff_poly_roots_d3(r,f);
					if ( ! k ) goto skip;
					if ( k>1 ) break;
					tecurve_AB_eval_df(&y,r[0],A[i]);
					if ( ! ff_residue(y) ) goto skip;																// if f'(tt) is not a residue, 4-torsion is not possible, o.w. 4-torsion must happen in rank 1 case
				} else {
					tecurve_AB_eval_df(&y,tt_x[i],A[i]);
					if ( ff_residue(y) ) break;																	// if f'(tt) is a residue, either we have 4-torsion or 2-rank 2
					if ( N==2 )goto skip;																	// in tor==2 case we know the 2-rank is 1, so we can't have 4-torsion
					if ( ! ff_residue(D) ) goto skip;																// Given f'(tt) not a residue, we must have 2-rank 2 to get N divisible by 4
				}
				break;
			}
		}
next:
		if ( t3_flag ) {
			_ff_set(f[1],A[i]); _ff_set(f[0],B[i]);
			k = ecurve_3tor(f);
			switch ( t3_flag ) {
			case FILTER_3TOR_1: if ( k> 1 ) { j--; continue;  } else break;
			case FILTER_3TOR_3: if ( k != 3 ) { j--; continue;  } else break;
			case FILTER_3TOR_NOT_1: if ( k <= 1 ) { j--; continue;  } else break;
			case FILTER_3TOR_9: if ( k != 9 ){ j--; continue;  } else break;
			}
		}
		if ( i==j ) continue;
		_ff_set(A[j],A[i]); _ff_set(B[j],B[i]);
		if ( tt_x ) _ff_set(tt_x[j],tt_x[i]);
		if ( base_x ) { _ff_set(base_x[j],base_x[i]); _ff_set(base_y[j],base_y[i]); }
		continue;
skip:	j--;
	}
	return j;
}

/*
	Computes a list of "random" curves of the form y^2=x^3+Ax+B containing an n-torision point and
	(optionally) satisfying constraints on the Sylow 2-subgroup as specified by s2_flag and constraints
	on the 3-torsion subgroup as specified by t3_flag.

	Note that incompatible n-torsion, s2_flag, and mod3 settings are not checked - the caller shouldn't ask for the impossible.
	For example,  FILTER_3TOR_9 should not be set unless p=1mod 3.

	It is also necessary that p be large enough for a curve with the specified parameters to exist.
	For example, for p=1481 there are no curves which have 13-torsion, 3-torsion, and full 2-torsion.
	Even for p=10039, there are no curves with 13-torsion and full 3-torsion and full 2-torsion (in general, asking for full torsion is a bad idea).

	We also compute a "random" point (x,y) on each curve, which with high probability, is *not* an n-torision point.

void tecurve_random_curves (ff_t f1[], ff_t x[], ff_t y[], int n, int tor, int s2_flag, int t3_flag)
{
	register ff_t b,c,d,f,m,r,s,t,u,t0,t1,t2,t3,t4,tt,bx,by;						// more variables than we need
	ff_t g[8],a[7];
	ff_t D,T;
	int i,k,basept,want_tt;
	
	if ( tor==2 && s2_flag == FILTER_2SYLOW_R2 ) s2_flag = 0;				// our 2-torsion method precludes this
	want_tt = (s2_flag == FILTER_2SYLOW_2 || s2_flag == FILTER_2SYLOW_R14 || s2_flag == FILTER_2SYLOW_R2 || s2_flag == FILTER_2SYLOW_D4 );
	tecurve_setup_constants();
for ( k = 0 ; k < n ; k++ ) {
retry:
	_ff_set_zero(tt);
	basept = 0;
	switch (tor) {
	case   1: // this case does not use the E(b,c) parameterization that the other families use
		ff3_random(a);												// get random element of F_p^3 which is not in F_p (ff3_random guarantees this)
		ff3_minpoly(g,a);											// compute its minimum polynomial, irreducible over F_p
		ff_depress_cubic(&T,g);										// depress to get f2=0, still irreducible
		goto check;
	case   2:  // this case does not use the E(b,c) parameterization that the other families use
		ff2_random(g);												// get a random element a of F_p^2 which is not in F_p
		ff2_norm(&D,g);											// construct f(x)=g(x)(x+g(0)) where g(x) is the min poly of a (irred quadratic over F_p)
		ff2_trace(&T,g);  
		_ff_set_one(g[3]);
		_ff_sub(g[2],D,T);
		_ff_set_one(t0);  _ff_subfrom(t0,T); _ff_mult(g[1],t0,D);
		_ff_square(g[0],D);
		ff_depress_cubic(&T,g);										// depress g to make g2 zero
		if ( s2_flag ) _ff_sub(tt,T,D);									// 2-torsion point is (tt,0)
		goto check;
	case   3:	// this case does not use the E(b,c) parameterization that the other families use
		// The curve y^2=x^3+(2bc-b^4/3)x + (2b^6/27 - 2b^3c/3 + c^2) has 3-torsion point (b^2/3,c) (c nonzero)
		_ff_random(b);  _ff_random_nz(c);
		_ff_square(t0,b);  _ff_mult(t1,t0,_ff_third);  _ff_mult(d,t0,t1);			// t1 = b^2/3, d = b^4/3
		_ff_mult(f,b,c); _ff_x2(f);										// f = 2bc
		_ff_sub(g[1],f,d);											// f1 = 2bc - b^4/3
		_ff_mult(t0,t1,d); _ff_x2(t0); _ff_mult(d,t0,_ff_third);				// d = 2/27*b^6
		_ff_mult(t0,t1,f); _ff_subfrom(d,t0); _ff_square(t0,c);				// d = 2/27*b^6 - 2u^3v/3, t0 = c^2
		_ff_add(g[0],d,t0);											// f0 = 2/27*b^6 - 2u^3v/3 + c^2
		_ff_set_one(g[3]);  _ff_set_zero(g[2]);
		goto check;
	case   4:
		_ff_random(d); _ff_square(t0,d); _ff_set_one(t1); _ff_subfrom(t1,t0); _ff_mult(b,t1,_ff_half); _ff_set_zero(c);						// b=(1-d^2)/2 (random d) c=0 (any b works, but this gives us an easy base point)
		_ff_add_one(bx,d);  _ff_sub(by,bx,b); basept = 1;																	// (bx,by)=(d+1,d+1-b) is infinite-order pt on E(b,c) over Q (our base point)
		// don't bother with tt, only relevant for R2 check
		break;
	case   5:
		_ff_random(bx);																							// pick random x-coordinate for our base point
		_ff_add(t0,bx,bx); _ff_add(t1,t0,t0); _ff_inc(t1); _ff_addto(t0,t1); _ff_inc(t0); if ( _ff_zero(t0) ) goto retry;						// invert 6x+2
		_ff_invert(b,t0); _ff_mult(t0,b,t1); _ff_mult(c,t0,bx); 																// c = x(4x+1)/(6x+2) (random x)
		_ff_square(t1,bx); _ff_add(t0,bx,bx); _ff_inc(t0); _ff_mult(by,t0,t1); _ff_mult(t0,by,b); _ff_add(by,t0,t0);							// y = x^2(2x+1)/(3x+1) (note this is in Kubert coords, not Atkins-Morain)
		_ff_set(b,c);  basept = 1;																					// b = c
		break;
	case   6: _ff_random_nz(c); _ff_square(t0,c); _ff_add(b,t0,c); _ff_set(tt,c);													// b=c+c^2 (random c), 2-torsion point has x-coord c
		break;
	case   7: _ff_random_nz(r); _ff_square(t0,r); _ff_mult(t1,t0,r); _ff_sub(b,t1,t0); _ff_sub(c,t0,r); break;								// b=cr c = r(r-1) (random r=s)
	case   8: _ff_random_nz(r); _ff_sub_one(t1,r); _ff_add(t0,r,t1); _ff_mult(b,t0,t1); _ff_invert(t0,r); _ff_mult(c,t0,b); break;				// b=(2r-1)(r-1) c=b/r (random r)
	case   9: _ff_random_nz(s); _ff_sub_one(t1,s); _ff_mult(r,s,t1); _ff_mult(c,s,r); _ff_inc(r); _ff_mult(b,c,r); break;						// b=cr c=s(r-1) r=s(s-1)+1 (random s)
	case 10: _ff_random_nz(s); _ff_sub_one(d,s); _ff_square(t0,d); _ff_sub(t1,s,t0); if ( _ff_zero(t1) ) goto retry;
			_ff_invert(t0,t1); _ff_square(t1,s); _ff_mult(r,t0,t1); _ff_sub_one(t0,r); _ff_mult(c,s,t0); _ff_mult(b,c,r);						// b=cr c=s(r-1) r=s^2/(s-(s-1)^2) (random s)
			if ( want_tt ) { _ff_mult(t0,r,s); _ff_mult(tt,t0,d); }																// 2-torsion point has x-coord 5P_x = rs(s-1)
			break;
	case 11:
		do {																									// X_1(11): V^2+V+U^3+U^2=0
			_ff_random_nz(bx); _ff_square(t0,bx); _ff_mult(t1,t0,bx); _ff_addto(t1,t0);
			_ff_x2(t1); _ff_x2(t1); _ff_dec(t1); _ff_neg(g[0],t1);																// compute discriminant D of quadratic directly
		} while ( ! ff_sqrt(a,g) );
		_ff_add_one(t0,a[0]); _ff_mult(t2,t0,_ff_half);  _ff_neg(r,t2);  _ff_invert(t1,bx); _ff_mult(s,t1,t2); _ff_inc(s);						// V'=V+1 = (sqrt(D)+1)/2, r-1=-V', s = V'/U + 1
		_ff_mult(c,r,s); _ff_inc(r); _ff_mult(b,c,r);																			// c=s(r-1), b=cr
		break;
	case 12:
		_ff_random_nz(bx); _ff_invert(d,bx); _ff_sub_one(t0,bx);  _ff_add(t1,t0,t0); _ff_add(r,d,t1); _ff_add(t1,r,t0); _ff_mult(s,t1,d);			// r = (2*U^2-2*U+1) / U,  s = (3*U^2-3*U+1) / U^2
		_ff_sub_one(t0,r); _ff_mult(c,s,t0); _ff_mult(b,c,r);																	// c = s(r-1), b = cr
		break;
	case 13:
		_ff_set_one(g[2]);
		do {																									// get a random point on X_1(13): V^2+(U^3+U^2+1)V-U^2-U=0
			_ff_random(bx);  _ff_add_one(t1,bx); _ff_mult(t0,t1,bx); _ff_neg(g[0],t0);  ff_mult(t2,t0,bx); _ff_add_one(g[1],t2);
		} while ( ! ff_poly_roots_d2(a,g,2) );
		_ff_mult(t2,t1,a[0]); _ff_add(t3,t1,a[0]);  _ff_addto(t2,t3); _ff_mult(d,t2,t3); if ( _ff_zero(d) ) goto retry; _ff_invert(d,d);			// invert (U'+V)(U'V+U'+V) where U'=U+1
		_ff_mult(t0,t3,d); _ff_mult(t3,t0,a[0]); _ff_neg(r,t3);  _ff_mult(t0,t2,d); _ff_mult(s,t0,t1);									// r-1 = -V(U'V+U'+V), s = U'/(U'+V)
		_ff_mult(c,r,s); _ff_inc(r); _ff_mult(b,c,r);
		break;

	case 14:
		_ff_set_one(g[2]);
		do {																									// get a random point on X_1(14): V^2 + (U+1)V - U^3 + U = 0
			_ff_random(bx); _ff_add_one(g[1],bx); _ff_square(t0,bx); _ff_dec(t0); _ff_mult(t2,t0,bx); _ff_neg(g[0],t2);
		} while ( ! ff_poly_roots_d2(a,g,2) );
		_ff_dec(bx);  _ff_sub(t2,bx,a[0]); _ff_mult(t3,bx,a[0]); _ff_addto(t3,t2); _ff_mult(d,t2,t3); if( _ff_zero(d) ) goto retry; _ff_invert(d,d);	// invert (U'-V)(U'V+U'-V) where U'=U-1
		_ff_mult(t1,t3,d); _ff_mult(s,bx,t1); _ff_mult(t1,t2,d); _ff_mult(r,t1,a[0]);												// s = U'/(U'-V), r-1 = V/(U'V+U'-V)
		_ff_mult(c,r,s); _ff_inc(r); _ff_mult(b,c,r);																			// c = s(r-1), b = cr
		if ( want_tt ) {
			_ff_sub(t0,r,s); _ff_square(d,t0); if ( _ff_zero(d) ) goto retry; _ff_invert(d,d);
			_ff_dec(s); _ff_mult(t0,r,s); _ff_subfrom(t0,r); _ff_inc(t0); _ff_mult(t1,t0,s); _ff_mult(t0,t1,b); _ff_mult(tt,t0,d);				// tt = 7P x-coord = s(s-1)r(r-1)(rs-2r+1)/(r-s)^2 = b(s-1)(r(s-1)-r+1)/(r-s)^2
		}
		break;
	case 15:
		_ff_set_one(g[2]);
		do {																									// get a random point on X_1(15): V^2+(U+1)V-U^3-U^2=0												
			_ff_random(bx); _ff_square(t0,bx); _ff_mult(t1,t0,bx); _ff_add(t,t0,t1); _ff_neg(g[0],t); _ff_add_one(g[1],bx);					// remember bx=U, t0=U^2, t1=U^3
		} while ( ! ff_poly_roots_d2(a,g,2) );
		_ff_square(t1,g[1]); _ff_sub(t2,t1,bx); _ff_mult(t3,t1,t2); _ff_mult(d,t3,a[0]); if ( _ff_zero(d) ) goto retry; _ff_invert(d,d);			// invert (U+1)^2(U^2+U+1)V
		_ff_add(t0,bx,a[0]); _ff_mult(t3,t0,bx); _ff_neg(t0,t3); _ff_mult(t3,d,t1); _ff_mult(s,t0,t3); _ff_mult(t3,d,t2); _ff_mult(r,t0,t3);		// s-1 = -U(U+V)/((U^2+U+1)V), r-1 = -U(U+V)/((U+1)^2*V)
		_ff_inc(s); _ff_mult(c,s,r); _ff_inc(r); _ff_mult(b,c,r);																	// c = s(r-1), b = cr
		break;
	case 16:
		_ff_set_one(g[2]);
		do {																									// get a random point on X_1(16): V^2+(U^3+U^2-U+1)V+U^2
			_ff_random(bx);  _ff_square(g[0],bx); _ff_add_one(t1,bx); _ff_mult(t2,g[0],t1); _ff_sub_one(t1,bx); _ff_sub(g[1],t2,t1);
		} while ( ! ff_poly_roots_d2(a,g,2) );
		_ff_add(t2,bx,a[0]);  _ff_mult(t0,t1,t2); _ff_mult(d,a[0],t0); if ( _ff_zero(d) ) goto retry; _ff_invert(d,d);							// invert V(U+V)(U-1)
		_ff_add_one(t0,a[0]);  _ff_mult(t3,t1,d); _ff_mult(r,t0,t3); _ff_mult(t3,t2,d); _ff_mult(s,t0,t3); _ff_inc(s);							// r-1 = (V+1)/(V(U+V)), s-1 = (V+1)/(V(U-1))
		_ff_mult(c,s,r); _ff_inc(r); _ff_mult(b,c,r);																			// c = s(r-1), b = cr
		break;
	case 17:
		_ff_set_one(g[4]);
		do {																									// X_1(17): V^4+(U^3+U^2-U+2)V^3+(U^3-3U+1)V^2-U(U^3+2)V+U^3+U^2
			_ff_random(bx); _ff_square(t0,bx); _ff_mult(t1,t0,bx); _ff_add(g[0],t0,t1); _ff_sub(t3,g[0],bx); _ff_add(g[3],t3,c_2);
			_ff_add(t3,bx,bx); _ff_addto(t3,bx); _ff_sub(t4,t1,t3); _ff_add_one(g[2],t4);
			_ff_add(t3,t1,c_2); _ff_mult(t4,t3,bx); _ff_neg(g[1],t4);
		} while ( ! ff_poly_roots_d4(a,g) );
		_ff_add(t1,bx,a[0]); _ff_inc(t1); _ff_mult(t3,bx,a[0]); _ff_sub(t2,t3,a[0]); _ff_dec(t2); _ff_mult(t4,t2,a[0]);						// invert V(UV-V-1)(U+V+1)
		 _ff_mult(d,t4,t1); if ( _ff_zero(d) ) goto retry; _ff_invert(d,d);  _ff_addto(t3,bx);											// r-1 = (UV+U)/(V(UV-V-1), s-1 = -(UV+V)/(V(U+V+1)
		_ff_mult(t0,d,t1); _ff_mult(r,t0,t3);  _ff_neg(t1,t3); _ff_mult(t0,d,t2); _ff_mult(s,t0,t1); _ff_inc(s);
		_ff_mult(c,s,r); _ff_inc(r); _ff_mult(b,c,r);																			// c = s(r-1), b = cr
		break;
	case 18:
		_ff_set_one(g[2]);
		do {																									// X_1(18): V^2 + (U^3+2U^2+3U-1)V - 2U = 0
			_ff_random(bx); _ff_square(t0,bx); _ff_mult(t1,t0,bx); _ff_add(t2,bx,bx); _ff_neg(g[0],t2);
			_ff_add(t3,t0,bx); _ff_addto(t2,t3); _ff_addto(t2,t0); _ff_dec(t2); _ff_add(g[1],t1,t2);
		} while ( ! ff_poly_roots_d2(a,g,2) );
		_ff_sub(t2,t1,t0); _ff_mult(t0,bx,a[0]); _ff_subfrom(t2,t0); _ff_subfrom(t2,a[0]); _ff_inc(t2); _ff_mult(d,t2,t1);
		if ( _ff_zero(d) ) goto retry;  _ff_invert(d,d);																		// invert (U^3-U^2-UV-V+1)U^3
		_ff_dec(t3); _ff_addto(t3,a[0]); _ff_mult(t0,d,t1); _ff_mult(r,t0,t3); _ff_mult(t0,d,t2); _ff_mult(s,t0,t3); _ff_inc(s);					// r-1 = (U^2+U+V-1) / (U^3-U^2-UV-V+1), s-1 = (U^2+U+V-1)/ U^3
		_ff_mult(c,s,r); _ff_inc(r); _ff_mult(b,c,r);																			// c = s(r-1), b = cr
		if ( want_tt ) {																								
			_ff_add(t0,r,s); _ff_dec(t0); _ff_square(t3,s);  _ff_subfrom(t0,t3); _ff_square(d,t0); if ( _ff_zero(d) ) goto retry; _ff_invert(d,d);	// invert 9Px denominator = (r-s^2+s-1)^2
			_ff_mult(t1,r,s);  _ff_sub(t0,s,c_3); _ff_mult(t2,t0,t1); _ff_addto(t2,r); _ff_addto(t2,t3);
			_ff_x2(r); _ff_subfrom(t1,r); _ff_inc(t1); _ff_mult(t0,c,t1); _ff_mult(t1,t0,t2); _ff_mult(tt,t1,d);								// c_2-torsion point has x-coord 9Px = c(rs-2r+1)(rs^2-3rs+r+s^2) / (r-s^2+s-1)^2
		}																										
		break;
	case 19:
		_ff_set_one(g[5]);																							// X_1(19): V^5 - dV^4 - 2U(d-U^2+1)V^3 - U(U^4-2U^2+7U-2)V^2 - dUV - U^2
		do {																									//      where d = U^3+2U^2-3U+1
			_ff_random(bx);  _ff_square(t0,bx); _ff_neg(g[0],t0); _ff_mult(t1,t0,bx); _ff_add(t2,t0,t0); _ff_add(d,t1,t2);
			_ff_add(t4,bx,bx); _ff_addto(t4,bx); _ff_subfrom(d,t4); _ff_inc(d); _ff_neg(g[4],d); _ff_mult(t3,d,bx); _ff_neg(g[1],t3);
			_ff_subfrom(d,t0);  _ff_inc(d); _ff_mult(t3,d,bx); _ff_x2(t3); _ff_neg(g[3],t3); _ff_square(t3,t0); _ff_subfrom(t3,t2);
			_ff_x2(t4); _ff_addto(t4,bx); _ff_addto(t3,t4); _ff_subfrom(t3,c_2); _ff_mult(t4,t3,bx); _ff_neg(g[2],t4);
		} while ( ! ff_poly_distinct_roots (a,g,5) );
		_ff_sub_one(t1,a[0]); _ff_add(t2,bx,t1); _ff_add(t3,t2,a[0]); _ff_add(t4,bx,a[0]); _ff_mult(t0,bx,t4);  _ff_addto(t3,t0);
		_ff_mult(t0,t1,t2); _ff_mult(d,t0,t3); if ( _ff_zero(d) ) goto retry;  _ff_invert(d,d);											// invert (V-1)(U+V-1)(U^2+UV+U+2V-1)
		_ff_inc(bx); _ff_mult(t0,bx,t4); _ff_mult(t1,d,t2); _ff_mult(r,t0,t1); _ff_mult(t1,t3,d); _ff_mult(s,t1,t4); _ff_inc(s);					// r-1 = ((U+1)(U+V))/((V-1)(U^2+UV+U+2V-1), s-1 = (U+V)/((V-1)(U+V-1)
		_ff_mult(c,s,r); _ff_inc(r); _ff_mult(b,c,r);																			// c = s(r-1), b = cr
		break;
	case 20:
		_ff_set_one(g[3]);
		_ff_set(g[0],c_2);
		do {																									// X_1(20): V^3+ (U^2+3)V^2 + (U^3+4)V + 2
			_ff_random(bx); _ff_square(t0,bx); _ff_mult(t1,t0,bx); _ff_add(g[2],t0,c_3); _ff_add(g[1],t1,c_4);
		} while ( ! ff_poly_roots_d3(a,g) );
		_ff_add_one(t1,a[0]); _ff_square(t2,t1); _ff_sub(t3,t2,t0); _ff_addto(t3,bx); _ff_sub_one(t0,bx); _ff_sub(t4,t2,t0);
		_ff_mult(t,t2,t3); _ff_mult(d,t,t4); if ( _ff_zero(d) ) goto retry; _ff_invert(d,d); _ff_mult(t,a[0],t0); _ff_addto(t2,t);					// invert (V+1)^2((V+1)^2-U^2+U)((V-1)^2-(U-1))
		_ff_mult(t,t2,t0); _ff_add(t0,bx,a[0]); _ff_mult(t2,t,t0);  ff_mult(t0,t4,d); _ff_mult(r,t0,t2);									// r-1 = (U-1)(U+V)((V+1)^2+V(U-1)) /  (V+1)^2((V+1)^2-U^2+U)
		_ff_mult(t0,d,t1); _ff_mult(t4,t0,t3); _ff_mult(s,t4,t); _ff_inc(s);														// s-1 = (U-1)((V+1)^2+V(U-1))/((V-1)^2-(U-1))
		_ff_mult(c,s,r); _ff_inc(r); _ff_mult(b,c,r);																			// c = s(r-1), b = cr		
		break;
	case 21:
		_ff_set_one(g[4]);
		do {																									// X_1(21): V^4+(3U^2+1)+(U^5+U^4+2U^2+2U)V^2+(2U^4+U^3+U)V+U^3
			_ff_random(bx); _ff_square(t0,bx); _ff_mult(t1,t0,bx); _ff_add(t2,t0,t0); _ff_addto(t2,t0); _ff_add_one(g[3],t2);
			_ff_square(t3,t0); _ff_add(t4,t3,t1); _ff_add_one(t2,bx); _ff_x2(t2); _ff_addto(t4,t2); _ff_mult(g[2],bx,t4);
			_ff_x2(t3); _ff_addto(t3,t1); _ff_add(g[1],t3,bx); _ff_set(g[0],t1);
		} while ( ! ff_poly_roots_d4(a,g) );
		_ff_mult(t1,bx,a[0]); _ff_inc(t1); _ff_square(t0,a[0]); _ff_sub(t2,t1,t0); _ff_mult(d,t1,t2); if ( _ff_zero(d) ) goto retry; _ff_invert(d,d);	// invert (UV+1)(UV+1-V^2)
		_ff_add(t3,t0,a[0]); _ff_mult(t4,t2,d); _ff_mult(s,t3,t4); _ff_inc(s); _ff_addto(t1,a[0]); _ff_mult(t2,t1,t3); _ff_mult(r,t2,d);			// s-1 = (V^2+V)/(UV+1), r-1 = (UV+1+V)(V^2+V)/((UV+1)(UV+1-V^2)
		_ff_mult(c,s,r); _ff_inc(r); _ff_mult(b,c,r);																			// c = s(r-1), b = cr		
		break;
	case 22:
		_ff_set_one(g[4]);																							// X_1(22): V^4 + (U^3+2U^2+U+2)V^3 + (U^5+U^4+2U^3+2U^2+1)V^2
		do {																									//                           (U^5-U^4-2U^3-U^2-U)V - U^4-U^3
			_ff_random(bx); _ff_square(t0,bx); _ff_mult(t1,t0,bx); _ff_square(t2,t0); _ff_add(t3,t1,t2); _ff_neg(g[0],t3);
			_ff_addto(t3,t1); _ff_addto(t3,t0); _ff_mult(t4,t0,t1); _ff_sub(t,t4,t3); _ff_sub(g[1],t,bx); _ff_addto(t4,t3);
			_ff_addto(t4,t0); _ff_add_one(g[2],t4);  _ff_add(t,t0,t0); _ff_add(t4,t1,t); _ff_addto(t4,bx); _ff_add(g[3],t4,c_2);
		} while ( ! ff_poly_roots_d4(a,g) );
		_ff_addto(t1,t); _ff_addto(t1,a[0]); _ff_add(t2,t0,a[0]); _ff_mult(d,t1,t2); if ( _ff_zero(d) ) goto retry; _ff_invert(d,d);				// invert (U^3+2U^2+V)(U^2+V)
		_ff_mult(t3,bx,a[0]); _ff_add(t4,t3,a[0]); _ff_mult(t,t1,d); _ff_mult(s,t,t4); _ff_add_one(t1,a[0]); _ff_mult(t,t0,t1);					// s = (UV+V)/(U^2+V), r = (U^2V+U^2+UV+V)/(U^3+2U^2+V)
		_ff_addto(t,t4); _ff_mult(t4,d,t2); _ff_mult(r,t,t4);
		_ff_sub_one(t0,r); _ff_mult(c,s,t0); _ff_mult(b,c,r);																	// c = s(r-1), b = cr		
		if ( want_tt ) {																								
			_ff_square(t0,r); _ff_square(t1,s); _ff_mult(t2,r,s); _ff_mult(t3,t2,s); _ff_mult(t4,s,t3); _ff_sub(d,t0,t4);
			_ff_add(t,t3,t3); _ff_addto(t,t3); _ff_addto(d,t); _ff_add(t,t2,t2); _ff_add(t4,t,t); _ff_subfrom(d,t4); _ff_addto(d,s);
			 if ( _ff_zero(d) ) break; _ff_square(t4,d); _ff_invert(d,t4); _ff_addto(t,t2); _ff_sub(t4,t3,t); _ff_addto(t4,r); _ff_addto(t4,t1);	// invert p11x denominator = (r^2-rs^3+3rs^2-4rs+s)^2, but don't retry on failure
			_ff_sub_one(t,s); ff_mult(t,t,b); ff_mult(t,t,t4);  _ff_sub(t4,s,c_3); ff_mult(t4,t4,t0); _ff_addto(t4,t2); _ff_addto(t4,r);
			_ff_addto(t4,r); _ff_addto(t4,r); _ff_subfrom(t4,t1); _ff_dec(t4); ff_mult(t,t,t4); _ff_mult(tt,t,d);							// p11x numerator is b(s-1)(rs^2-3rs+r+s^2)(r^2s-3r^2+rs+3r-s^2-1)
		}																										
		break;
	case 23:
		_ff_set_one(g[7]);							 // V^7 + (U^4-4U^2)V^6 - (U^6+U^5-2U^4-4U^3-3U^2+4U)V^5 + (2U^7-3U^5-6U^4+2U^3+5U^2-U)V^4-(U^8+3U^7-9U^6-U^5+7U^4-U)V^3
		do {									//	    + U(U-1)^2(3U^5-5U^3-3U^2+3U-1)V^2 - U^4(U-1)^3(3U-4)V+U^2(U-1)^6
			_ff_random(bx); _ff_square(t0,bx); _ff_square(t2,t0); _ff_add(t1,t0,t0); _ff_x2(t1); _ff_sub(g[6],t2,t1); _ff_mult(t1,t0,bx);
			_ff_mult(t3,t1,t0); _ff_add(t,t3,g[6]); _ff_subfrom(t,t1); _ff_subfrom(t,t1); _ff_add(t4,bx,bx); _ff_addto(t4,bx);
			_ff_subfrom(t4,c_4); _ff_subfrom(t,t4); _ff_mult(d,t,bx); _ff_neg(g[5],d); _ff_sub_one(m,bx); _ff_square(by,m); ff_mult(m,m,by); 
			_ff_mult(t,t2,m); _ff_mult(d,t,t4); _ff_neg(g[1],d); _ff_square(t,m); _ff_mult(g[0],t0,t); _ff_sub(t,t3,t1); _ff_subfrom(t,t0); _ff_addto(t,bx);
			_ff_add(t4,t,t); _ff_addto(t,t4); _ff_add(t4,t1,t1); _ff_subfrom(t,t4); _ff_dec(t); _ff_mult(t4,t,by); _ff_mult(g[2],t4,bx);  
			_ff_mult(t4,t1,c_3); _ff_add(t,t2,t4); _ff_mult(t4,t0,c_8); _ff_subfrom(t,t4); _ff_subfrom(t,t0); _ff_subfrom(t,bx); _ff_addto(t,c_8); _ff_dec(t);
			_ff_mult(d,t,t2); _ff_subfrom(d,bx); _ff_neg(g[3],d);  _ff_square(t,t1); _ff_subfrom(t,t2); _ff_subfrom(t,t1); _ff_addto(t,t0); _ff_addto(t,bx); _ff_x2(t); 
			_ff_subfrom(t,t2); _ff_x2(t1); _ff_x2(t1); _ff_subfrom(t,t1); _ff_add(t4,bx,bx); _ff_addto(t4,bx); _ff_addto(t,t4); _ff_mult(d,t,bx); _ff_sub(g[4],d,bx);
		} while ( ! ff_poly_distinct_roots(a,g,7) );
		_ff_sub(t2,t0,bx); _ff_subfrom(t2,a[0]); _ff_mult(d,t2,bx); if ( _ff_zero(d) ) goto retry; _ff_invert(d,d);							// invert U(U^2-U-V)
		_ff_mult(t3,t2,d); _ff_sub(t1,bx,a[0]); _ff_mult(s,t3,t1);
		_ff_dec(t1); _ff_neg(t0,t1); _ff_mult(t2,t0,a[0]); _ff_mult(t3,bx,d); _ff_mult(r,t2,t3);										// s = (U-V)/U, r-1 = V(V-U+1)/(U^2-U-V)
		_ff_mult(c,r,s); _ff_inc(r); _ff_mult(b,c,r);																			// c = s(r-1), b = cr		
		break;
	case 24:
		_ff_set_one(g[5]);
		do {
			_ff_random(bx); _ff_square(t0,bx); _ff_mult(t1,t0,bx); _ff_square(t2,t0); _ff_add(d,t1,t1); _ff_x2(d);
			_ff_sub(t,t2,d); _ff_add(m,t0,t0); _ff_addto(m,t0); _ff_addto(t,m); _ff_addto(t,bx); _ff_sub(g[4],t,c_2);
			_ff_add(by,t2,t2); _ff_add(t4,d,d); _ff_sub(t,by,t4); _ff_add(t4,m,m); _ff_addto(t4,t0); _ff_addto(t,t4); _ff_dec(t); _ff_neg(g[3],t);
			_ff_sub(t,by,d); _ff_subfrom(t,m); _ff_add(d,bx,bx); _ff_x2(d); _ff_addto(d,bx); _ff_addto(t,d); _ff_dec(t); _ff_mult(g[2],bx,t);
			_ff_add(t,t0,t0); _ff_subfrom(t,d); _ff_addto(t,c_2); _ff_mult(t4,t,t1); _ff_neg(g[1],t4); _ff_sub(t,t0,bx); _ff_mult(g[0],t,t2);
		} while ( ! ff_poly_distinct_roots(a,g,5) );
		_ff_sub(t1,bx,a[0]); _ff_mult(t2,t1,bx); _ff_sub_one(t3,a[0]); _ff_mult(t4,t3,a[0]); _ff_subfrom(t2,t4); _ff_mult(d,t1,t2);
		if ( _ff_zero(d) ) goto retry; _ff_invert(d,d);  _ff_mult(t4,t2,d); _ff_sub_one(t,bx); _ff_mult(s,t,t4); _ff_sub(t4,t0,bx);				// invert (U-V)(U^2-UV-V^2+V)
		_ff_subfrom(t4,t3); _ff_mult(t2,t1,d); _ff_mult(r,t2,t4);																// s= (U-1)/(U-V), r = (U^2-U-V+1)/(U^2-UV-V^2+V), 
		_ff_sub_one(t0,r); _ff_mult(c,s,t0); _ff_mult(b,c,r);																	// c = s(r-1), b = cr		
		break;
	default: printf("Error, unhandled torsion parameter tor=%d\n", tor);  exit (0);
	}
//printf("b=%ld, c=%ld\n", _ff_get_ui(b), _ff_get_ui(c));	
	_ff_sub_one(t1,c);  _ff_square(t2,t1); _ff_add(t0,b,b); _ff_x2(t0); _ff_subfrom(t2,t0);											// t2=b2=a1^2+4a2^2=(1-c)^2-4b  (SAVE t2=b2)
	_ff_mult(t4,t1,b); _ff_square(t3,b);																					// t4=b4=2a4+a1a3=(1-c)b, t3=b6=a3^2+4a6=b^2
//printf("b2=%ld,b4=%ld,b6=%ld\n", _ff_get_ui(t2), _ff_get_ui(t4), _ff_get_ui(t3));
	_ff_mult(t0,c_24,t4); _ff_square(t1,t2); _ff_subfrom(t0,t1); _ff_mult(g[1],c_27,t0);												// A=-27c4=27(b4-b2^2)
	_ff_mult(t0,c_36,t4);  _ff_sub(t4,t1,t0); _ff_mult(t0,t4,t2); _ff_mult(t1,c_216,t3); _ff_addto(t0,t1); _ff_mult(g[0],c_54,t0);				// B = -54c6=54(b2(b2^2-36b4)+216b6)
	_ff_set_one(g[3]);  _ff_set_zero(g[2]);
//printf ("A=%ld, B=%ld\n", _ff_get_ui(g[1]), _ff_get_ui(g[0]));
//printf("tt on E(%ld,%ld) is %ld\n", _ff_get_ui(b), _ff_get_ui(c), _ff_get_ui(tt));
	if ( s2_flag && ! _ff_zero(tt) ) { _ff_mult(t0,tt,c_12); _ff_addto(t0,t2); _ff_add(tt,t0,t0); _ff_addto(tt,t0); }							// shift two torsion point x-coord to 36x+3b2 = 3(12x+b2)
//printf ("shifted to %ld on ", _ff_get_ui(tt));
//ff_poly_print(g,3);
check:
	tecurve_discriminant (&D,g);
	if ( _ff_zero(D) ) goto retry;																						// retry if curve is singular
	switch ( s2_flag ) {
	case FILTER_2SYLOW_1:
		if ( ! ff_residue(D) ) goto retry;																				// if discriminant is not a residue, then g(x) has 2 factors (Stickleberger) and y^2=f(x) has 2-torsion		
		break;																									// we could still have 2-torsion rank 2 here (g(x) with 3 factors), but we don't filter this
	case FILTER_2SYLOW_2:																							// note that tor=2 case already ensures 2-torsion rank is 1
		if ( _ff_zero(tt) ) {										
			if ( ff_poly_roots_d3(a,g) != 1 ) goto retry;																	// this is costly, much cheaper if 2-torsion point tt is known by construction
			_ff_set(tt,a[0]);
			_ff_square(t1,tt); _ff_add(t0,t1,t1); _ff_addto(t0,t1); _ff_addto(t0,g[1]);											// set t0 to 3tt^2+f1
			if ( ff_residue(t0) ) goto retry;
		} else {
			if ( (_ff_p&3) == 1 ) if ( tor > 2 && ff_residue(D) ) goto retry;
			_ff_square(t1,tt); _ff_add(t0,t1,t1); _ff_addto(t0,t1); _ff_addto(t0,g[1]);											// set t0 to 3tt^2+f1
			if ( ff_residue(t0) ) goto retry;
			if ( (_ff_p&3) == 3 ) if ( tor > 2 && ff_residue(D) ) goto retry;
		}
		break;
	case FILTER_2SYLOW_NOT_1:																						// currently we force rank 1 to achive non-trivial 2-torsion for odd torsion case
		if ( !(tor&1) ) break;
		// fall through
	case FILTER_2SYLOW_R1:
		if ( tor > 2 && ff_residue(D) ) goto retry;																		// if discriminant is a residue than we have 2-torsion rank 0 or 2 (but not 1)		
		break;
	case FILTER_2SYLOW_R14:
		if ( !(tor&3) ) { if ( ff_residue(D) ) goto retry; else break; }														// if 4-torsion is already present, just verify rank
		if ( _ff_zero(tt) ) {										
			if ( ff_poly_roots_d3(a,g) != 1 ) goto retry;																	// this is costly, much cheaper if 2-torsion point tt is known by construction
			_ff_set(tt,a[0]);
		} else {
			if ( tor > 2 && ff_residue(D) ) goto retry;																	// if discriminant is a residue than we have 2-torsion rank 0 or 2 (but not 1)
		}
		_ff_square(t1,tt); _ff_add(t0,t1,t1); _ff_addto(t0,t1); _ff_addto(t0,g[1]);												// set t0 to 3tt^2+f1
		if ( ! ff_residue(t0) ) goto retry;																				// if t0 is not a residue, 4-torsion is not possible  (o.w. 4-torsion must happen in rank 1 case)
		break;
	case FILTER_2SYLOW_R2:
		if ( ! _ff_zero(tt) ) {
			_ff_square(t0,tt); _ff_add(t1,t0,g[1]); _ff_x2(t1); _ff_x2(t1); _ff_subfrom(t0,t1);										// compute discriminant of g(x)/(x-tt) (we know g(tt)=0 since tt is the x-coord of 2-torsion pt)
			if ( ! ff_residue(t0) ) goto retry;
		} else {
			if ( ff_poly_roots_d3(0,g) < 3 ) goto retry;
		}
		break;
	case FILTER_2SYLOW_D4:
		if ( !(tor&3) ) break;
		if ( _ff_zero(tt) ) {
			i = ff_poly_roots_d3(a,g);
			if ( ! i ) goto retry;
			if ( i>1 ) break;
			_ff_square(t1,a[0]); _ff_add(t0,t1,t1); _ff_addto(t0,t1); _ff_addto(t0,g[1]);								
			if ( ! ff_residue(t0) ) goto retry;																			// if f'(tt) is not a residue, 4-torsion is not possible, o.w. 4-torsion must happen in rank 1 case
		} else {
			_ff_square(t1,tt); _ff_add(t0,t1,t1); _ff_addto(t0,t1); _ff_addto(t0,g[1]);
			if ( ff_residue(t0) ) break;																				// if f'(tt) is a residue, either we have 4-torsion or 2-rank 2
			if ( tor==2 ) goto retry;																					// in tor==2 case we know the 2-rank is 1, so we can't have 4-torsion
			if ( ! ff_residue(D) ) goto retry;																			// Given f'(tt) not a residue, we must have 2-rank 2 to get N divisible by 4
		}
		break;
	}
	if ( t3_flag ) {
		i = ff_poly_g1_3tor(g);
		switch ( t3_flag ) {
		case FILTER_3TOR_1: if ( i > 1 ) goto retry; break;
		case FILTER_3TOR_3: if ( i != 3 ) goto retry; break;
		case FILTER_3TOR_NOT_1: if ( i <= 1 ) goto retry; break;
		case FILTER_3TOR_9: if ( i != 9 ) goto retry; break;
		}
	}
	_ff_set(f1[k],g[1]);
	if ( basept ) {
		// translate base pt using x = 36bx+3u and y=216by-108((c-1)bx+b) where (bx,by) is the base pt on Kubert curve E(b,c).  Note that u=b2=(1-c)^2-4b is saved in t2 above
		_ff_add(t0,t2,t2); _ff_addto(t0,t2); _ff_mult(t1,c_36,bx); _ff_add(x[k],t0,t1);
		_ff_sub_one(t0,c); _ff_mult(t1,t0,bx); _ff_addto(t1,b); _ff_add(t0,by,by); _ff_subfrom(t0,t1); _ff_mult(y[k],t0,c_108);
//printf ("Base pt (%ld,%ld) shifted to (%ld,%ld)\n", _ff_get_ui(bx), _ff_get_ui(by), _ff_get_ui(x[k]), _ff_get_ui(y[k]));
	} else {
		// generate a random point the hard way
		ecurve_random_point(x+k,y+k,g);
	}
//printf("p=%ld: ", _ff_p); ff_poly_print(g,3);
//printf("Base pt (%ld,%ld)\n", _ff_get_ui(x[k]), _ff_get_ui(y[k]));
}
}
*/
