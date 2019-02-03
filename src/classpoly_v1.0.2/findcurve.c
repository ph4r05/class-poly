#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mpzutil.h"
#include "cstd.h"
#include "ecurve.h"
#include "ff_poly.h"
#include "tecurve.h"
#include "pickprimes.h"
#include "findcurve.h"

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

int ecurve_test_order2 (ppf_t n[2], ff_t f[4]);	// not defined in ecurve.h

#define MAX_CURVES	128
#define MIN_CURVES	32

unsigned long gop_counter;

int scancurves (int *psignm, ff_t x[], ff_t y[], int n, unsigned long pbits0, unsigned long nbits0, unsigned long pbits1, unsigned long nbits1, ff_t f1[]);
int scancurves1 (ff_t x1[], ff_t y1[], int n, int chain[FF_MAX_CHAIN_LEN], int chain_len, ff_t f1[]);
int compute_signed_exp_chain (int chain[FF_MAX_CHAIN_LEN], unsigned long e);

// searches for an elliptic curve over F_p with a_p = +/-t
// returns 0 for failure, otherwise +/- 1 to indicate the sign of t in a_p.
int findcurve (ff_t f[4], long t, int flags, long *cnt)
{
	time_t start, end;
	ff_t x1[MAX_CURVES];
	ff_t y1[MAX_CURVES];
	ff_t f1[MAX_CURVES];
	ff_t J, J1728;
	ff_t g[4], *h;
	ppf_t N[2];
	unsigned long pbits0, nbits0, pbits1, nbits1;
	long xcnt;
	register ff_t t0,t1;
	register int i, k;
	int sign, tsign, fsign, tor, s2_flag, t3_flag, j_flag, twist, nofilter;
	
	if ( ! cnt ) cnt = &xcnt;
	*cnt = 0;
	start = clock();
	ppf_factor(N[0],_ff_p+1-t); ppf_factor(N[1],_ff_p+1+t);
	tsign = 1;

	// It's rather silly to recompute these, since they were previously computed when the primes are selected, but it's a pain to remember them and pass them in
	if ( ! (flags&FINDCURVE_NOFILTER) ) split_prime_rating_new (_ff_p, t, &twist, &tor, &s2_flag, &t3_flag);

	if ( !(flags&FINDCURVE_NOFILTER) && twist<3 ) {
		if ( twist==1 ) { ui_NAF (&pbits0, &nbits0, _ff_p-t+1); /*chain_len =  compute_signed_exp_chain(chain,_ff_p-t+1); */fsign=1; }
		else {  ui_NAF (&pbits0, &nbits0, _ff_p+t+1); /*chain_len = compute_signed_exp_chain(chain, _ff_p+t+1);*/ fsign=-1; }
		ui_NAF (&pbits1, &nbits1, 0);
	} else {
		ui_NAF (&pbits0, &nbits0, _ff_p+1);
		ui_NAF (&pbits1, &nbits1, t);
		fsign = 0;
	}
//out_printf("p=%ld(%ld), t=%ld(%ld), tor=%d, s2_flag=%d, t3_flag=%d, fsign=%d\n", _ff_p, _ff_p%tor, t, t%tor, tor, s2_flag, t3_flag,  fsign);

	j_flag = ( (flags & FINDCURVE_FILTER_JCUBE) ? FILTER_JCUBE : 0 );
	if ( flags & FINDCURVE_FILTER_WEBER ) j_flag |= FILTER_WEBER;
	_ff_set_ui(J1728,1728);
	_ff_set_one(f[3]);  _ff_set_zero(f[2]);
	memset(x1,0,sizeof(ff_t)*MAX_CURVES);
	nofilter = flags&FINDCURVE_NOFILTER;
	for (;;) {
		if ( !nofilter && *cnt > 100*sqrt(_ff_p) ) { if ( _ff_p > 1000 ) info_printf ("Disabling curve filtering in find_curve due to slow progress with tor=%d at p=%ld\n", tor, _ff_p); nofilter = 1; }
		i = k = 0;
		do {
			if ( nofilter ) {
				do {
					_ff_random(f[1]); _ff_random(f[0]);
					ff_poly_x3axb_disc(&J,f);
				} while ( _ff_zero(J) );
				while ( ! ecurve_random_point (x1+k,y1+k,f) );
				_ff_set(f1[k],f[1]);
				k = k+1;
			} else{
				k+=tecurve_random_curves_x (f1+k,x1+k,y1+k,MAX_CURVES-k,tor,s2_flag,t3_flag,j_flag);
				// if tecurve_random_curves fails repeatedly, it is almost certainly due to a pathological problem caused by a small F_p (this happens primarily with tor=12).
				if ( k == 0 && ++i == 3 ) { dbg_printf ("Reducing tor=%d due to tecurve_random_curves failing thrice in a row with p=%ld\n", tor, _ff_p); if ( !(tor&1) ) tor /= 2; else nofilter =1; }
			}
		} while ( k < MIN_CURVES );
		*cnt += k;
		
//		The benefit of using scancurves1 is negligible so for simplicity we don't.
/*		
		if ( fsign ) {
			sign = fsign;
			k = scancurves1 (x1,y1, k, chain, chain_len, f1);
		} else {
			k = scancurves (&sign, x1,y1, k, pbits0, nbits0, pbits1, nbits1, f1);
		}
*/
		k = scancurves (&sign, x1,y1, k, pbits0, nbits0, pbits1, nbits1, f1);
		if ( k >= 0 ) {
			_ff_set (f[1], f1[k]);
			_ff_square(t1,y1[k]);
			_ff_square(t0,x1[k]);
			_ff_addto(t0,f1[k]);
			ff_mult(t0,t0,x1[k]);
			_ff_sub(f[0],t1,t0);
			if ( ! ecurve_to_jinv(&J,f) ) continue;											// make sure J-invariant is not 0 or 1728 to avoid annoying special cases
			if ( _ff_zero(J) || _ff_equal(J,J1728) ) continue;
			if ( fsign ) sign = fsign;
			if ( sign != tsign ) { ff_poly_twist(g,f,3); h = g;  } else { h = f; }
//printf ("Testing orders N=[%ld,%ld] for curve ", _ff_p+1-t, _ff_p+1+t); ff_poly_print (h, 3);
			if ( ! ecurve_test_order2 (N,h) ) { /* puts ("failed."); */ continue; }
//puts ("success!");
			end = clock();
			return sign;
		}
	}
	end=clock();
	err_printf ("p=%ld, t=%ld, tor=%d, failed after %ld curves in %ld msecs\n", _ff_p, t, tor, *cnt, delta_msecs(start,end));
	return 0;
}

static void inline parallel_invert (ff_t z[], ff_t x[], unsigned n)
{
	ff_t c[MAX_CURVES];
	register ff_t u, v;
	register unsigned i;

	if ( ! n ) return;
	_ff_set (c[0], x[0]);
	for ( i = 1 ; i < n ; i++ ) _ff_mult (c[i], c[i-1], x[i]);
	_ff_invert (u, c[n-1]);
	for ( i = n-1 ; i > 0 ; i-- ) {
		_ff_mult (v, c[i-1], u);
		_ff_mult (u, u, x[i]);
		_ff_set (z[i], v);
	}
	_ff_set (z[0], u);
}

// doubles (x1,y1) in place, given inv=1/(2y1) (a=f1)
static void inline affine_double (ff_t *x1, ff_t *y1, ff_t inv, ff_t f1)
{
	register ff_t t1, t2, t3;

	if ( _ff_zero(*y1) ) { _ff_set_zero(*x1); return; }
	_ff_square(t1,*x1);	_ff_add(t2,t1,t1);  _ff_add(t2,t2,t1); _ff_addto(t2,f1);	// t2 = 3x1^2+a
	_ff_mult(t1,t2,inv);												// t1 = lambda = (3x1^2+a)/(2y1)
	_ff_square(t3,t1);   _ff_add(t2,*x1,*x1);   _ff_subfrom(t3,t2);				// t3 = x3 = lamda^2 - 2x1
	_ff_sub(t2,*x1,t3);												// t2 = x1-x3
	_ff_set(*x1,t3);													// x3 = t3
	_ff_mult(t3,t1,t2);  _ff_sub(t1,t3,*y1);								// y3 = (x1-x3)lambda - y1
	_ff_set(*y1,t1);
	// 4M + 6A
}

// doubles (0,y1) into (x3,y3), given inv=1/(2y1) (a=f1)
static void inline affine_double_x0 (ff_t *x3, ff_t *y3, ff_t y1, ff_t inv, ff_t f1)
{
	register ff_t t1, t2, t3;

	_ff_mult(t1,f1,inv);												// t1 = lambda = (3x1^2+a)/(2y1) = a/(2y1)
	_ff_square(t2,t1); 												// t2 = x3 = lamda^2 - 2x1 = lambda^2
	_ff_mult(t3,t1,t2);  _ff_addto(t3,y1);								// t3 = (x3-x2)lambda + y1 = lambda^3 + y1
	_ff_set(*x3,t2);													// x3 = t2
	_ff_neg(*y3,t3);												// y3 = -t3
	// 3M + 2A
}

//  adds (x1,y1) to (x2,y2) (in place) given inv=1/(x2-x1) and assuming (x1,y1) != (x2,y2)
static void inline affine_add (ff_t *x2, ff_t *y2, ff_t x1, ff_t y1, ff_t inv)
{
	register ff_t t1, t2, t3, t4;

	if ( _ff_zero(x1)&& _ff_zero(y1) ) return;									// (x1,y1)=(0,0) is the identity
	if ( _ff_zero(*x2) && _ff_zero(*y2) ) { _ff_set(*x2,x1); _ff_set(*y2,y1); return; }	// (x2,y2)=(0,0) is the identity
	if ( _ff_equal(x1,*x2) ) { _ff_set_zero(*x2); _ff_set_zero(*y2); return; }			// assume pts are opposites
	
	_ff_sub(t2,*y2,y1);  _ff_mult(t1,t2,inv);								// t1 = lambda = (y2-y1)/(x2-x1)
	_ff_square(t3,t1); _ff_subfrom(t3,x1);  _ff_subfrom(t3,*x2);				// t3 = x3 = lambda^2 - x1 - x2
	_ff_sub(t2,x1,t3); _ff_mult(t4,t1,t2); _ff_sub(*y2,t4,y1);					// y3 = (x1-x3)lambda - y1
	_ff_set(*x2,t3);													// x3 = t3
	// 3M + 5A
}

//  subtracts (x1,y1) to (x2,y2) (in place) given inv=1/(x2-x1) and assuming (x1,y1) != (x2,y2)
static void inline affine_sub (ff_t *x2, ff_t *y2, ff_t x1, ff_t y1, ff_t inv)
{
	register ff_t t1, t2, t3, t4;

//printf("subtracting (%ld,%ld) to (%ld,%ld) using inv=%ld\n", _ff_get_ui(x1), _ff_get_ui(y1), _ff_get_ui(*x2), _ff_get_ui(*y2), _ff_get_ui(inv));

	if ( _ff_zero(x1)&& _ff_zero(y1) ) return;									// (x1,y1)=(0,0) is the identity
	if ( _ff_zero(*x2) && _ff_zero(*y2) ) { _ff_set(*x2,x1); _ff_neg(*y2,y1); return; }	// (x2,y2)=(0,0) is the identity
	if ( _ff_equal(x1,*x2) ) { _ff_set_zero(*x2); _ff_set_zero(*y2); return; }			// assume pts are equal
	
	_ff_add(t2,*y2,y1);  _ff_mult(t1,t2,inv);								// t1 = lambda = (y2+y1)/(x2-x1)
	_ff_square(t3,t1); _ff_subfrom(t3,x1);  _ff_subfrom(t3,*x2);				// t3 = x3 = lambda^2 - x1 - x2
	_ff_sub(t2,x1,t3); _ff_mult(t4,t1,t2); _ff_add(*y2,t4,y1);					// y3 = (x1-x3)lambda + y1
	_ff_set(*x2,t3);													// x3 = t3
//printf("result (%ld,%ld)\n", _ff_get_ui(*x2), _ff_get_ui(*y2));
	// 3M + 5A
}


// doubles n affine points (in place) using one inversion
static void inline parallel_double (ff_t x[], ff_t y[], int n, ff_t f1[])
{
	ff_t s[MAX_CURVES];
	ff_t t[MAX_CURVES+1];
	register ff_t t0, t1;
	register int i;

	_ff_set_one(t[0]);
	for ( i = 0 ; i < n ; i++ ) {
		if ( _ff_zero(y[i]) ) { _ff_set_one(s[i]); } else { _ff_add(s[i], y[i], y[i]); }		// just set 2y to a non-zero value (e.g. 1) for 2-torsion elements (the inverse won't get used by affine_double)
		_ff_mult(t[i+1], t[i], s[i]);
	}
	_ff_invert (t1, t[n]);
	for ( i = n-1 ; i > 0 ; i-- ) {
		_ff_mult (t0, t[i], t1);
		_ff_mult (t1, t1, s[i]);
		affine_double (x+i,y+i,t0,f1[i]);
	}
	affine_double(x,y,t1,f1[i]);
	// Total cost is (7M+7A)*n + I - 2M
gop_counter += n;
}

// sets (x[i],y[i]) = (x[i],y[i]) + (x1[i],y1[i]) for i = 0 to n-1
static void inline parallel_add (ff_t x[], ff_t y[], ff_t x1[], ff_t y1[], int n, ff_t f1[])
{
	ff_t s[MAX_CURVES];
	ff_t t[MAX_CURVES+1];
	register ff_t t0, t1;
	register int i;

	_ff_set_one(t[0]);
	for ( i = 0 ; i < n ; i++ ) {
		if ( _ff_equal(x[i],x1[i]) ) {
			if ( ! _ff_zero(y[i]) && _ff_equal(y[i],y1[i]) ) { _ff_add(s[i], y[i], y[i]); } else { _ff_set_one(s[i]); }
		} else {
			_ff_sub(s[i], x[i], x1[i]);
		}
		_ff_mult(t[i+1], t[i], s[i]);
	}
	_ff_invert (t1, t[n]);
	for ( i = n-1 ; i > 0 ; i-- ) {
		_ff_mult (t0, t[i], t1);
		_ff_mult (t1, t1, s[i]);
		// one could play with the order of these tests to try to improve branch prediction, but it's not obvious this would help
		if ( ! _ff_zero(x[i]) ) { affine_add (x+i,y+i,x1[i],y1[i],t0); continue; }
		if ( _ff_zero(y[i]) ) { _ff_set(y[i],y1[i]); continue; }
		if ( ! _ff_equal(y[i],y1[i]) ) { _ff_set_zero(y[i]);  continue; }
		affine_double (x+i,y+i,t0, f1[i]);
	}
gop_counter += n;
	// i==0 here, but oddly enough the compiler likes the generic version better (matches loop code)
	if ( ! _ff_equal(x[i],x1[i]) ) { affine_add (x+i,y+i,x1[i],y1[i],t1); return; }
	if ( _ff_zero(y[i]) ) { _ff_set(y[i],y1[i]); return; }
	if ( ! _ff_equal(y[i],y1[i]) ) { _ff_set_zero(y[i]);  return; }
	affine_double (x+i,y+i, t1, f1[i]);
	// total cost is approximately (6M+6A)*n + I - 2M
}


// sets (x[i],y[i]) = (x[i],y[i])  - (x1[i],y1[i]) for i = 0 to n-1
static void inline parallel_sub (ff_t x[], ff_t y[], ff_t x1[], ff_t y1[], int n, ff_t f1[])
{
	ff_t s[MAX_CURVES];
	ff_t t[MAX_CURVES+1];
	register ff_t t0, t1;
	register int i;

	_ff_set_one(t[0]);
	for ( i = 0 ; i < n ; i++ ) {
		if ( _ff_equal(x[i],x1[i]) ) {
			if ( ! _ff_zero(y[i]) && _ff_equal(y[i],y1[i]) ) { _ff_add(s[i], y[i], y[i]); } else { _ff_set_one(s[i]); }
		} else {
			_ff_sub(s[i], x[i], x1[i]);
		}
		_ff_mult(t[i+1], t[i], s[i]);
	}
	_ff_invert (t1, t[n]);
	for ( i = n-1 ; i > 0 ; i-- ) {
		_ff_mult (t0, t[i], t1);
		_ff_mult (t1, t1, s[i]);
		// one could play with the order of these tests to try to improve branch prediction, but it's not obvious this would help
		if ( ! _ff_zero(x[i]) ) { affine_sub (x+i,y+i,x1[i],y1[i],t0); continue; }
		if ( _ff_zero(y[i]) ) { _ff_neg(y[i],y1[i]); continue; }
		if ( _ff_equal(y[i],y1[i]) ) { _ff_set_zero(y[i]);  continue; }
		affine_double (x+i,y+i,t0, f1[i]);
	}
gop_counter += n;
	// i==0 here, but oddly enough the compiler likes the generic version better (matches loop code)
	if ( ! _ff_equal(x[i],x1[i]) ) { affine_sub (x+i,y+i,x1[i],y1[i],t1); return; }
	if ( _ff_zero(y[i]) ) { _ff_neg(y[i],y1[i]); return; }
	if ( _ff_equal(y[i],y1[i]) ) { _ff_set_zero(y[i]);  return; }
	affine_double (x+i,y+i, t1, f1[i]);
	// total cost is approximately (6M+4A)*n + I - 2M
}

// Tests whether (0,y1[i])^j=(0,y1[i])^{+/-k} on the curve y^2=x^3+f1*x+(y1[i])^2 for i from 0 to n-1 with j=pbits0-nbits0 and k=pbits1-nbits1 (NAF exponents)
// assumes (0,0) is not a point on the curve, that pbits0&nbits0==pbits1&nbits1==0, and that j>0
// sets *psign to +1 or -1 to indicate sign of match
int scancurves (int *psign, ff_t x1[], ff_t y1[], int n, unsigned long pbits0, unsigned long nbits0, unsigned long pbits1, unsigned long nbits1, ff_t f1[])
{
	ff_t x[MAX_CURVES];
	ff_t y[MAX_CURVES];
	ff_t Px[MAX_CURVES], Py[MAX_CURVES], Qx[MAX_CURVES], Qy[MAX_CURVES];
	register int i,len, Pset, Qset;
	register unsigned long m,max;

//_ff_square(t0,y1[0]);
//printf ("\np=%ld, n=%ld, P1=(0,%ld)  f=x^3+%ld*x+%ld, %d=%x(hex)-%x(hex), %d=%x(hex)-%x(hex)\n", _ff_p, n, _ff_get_ui(y1[0]), _ff_get_ui(f1[0]), _ff_get_ui(t0), pbits0-nbits0, pbits0, nbits0, pbits1-nbits1, pbits1, nbits1);
	Pset = Qset = 0;														// flags to indicates P and Q values unset
	len = n*sizeof(ff_t);
	if ( pbits0&1 ) {memcpy (Px,x1,len); memcpy(Py,y1,len); Pset = 1; }
	if ( nbits0&1 ) { memcpy (Px,x1,len);  for ( i = 0 ; i < n ; i++ ) {_ff_neg(Py[i],y1[i]); } Pset = 1; }
	if ( pbits1&1 ) { memcpy (Qx,x1,len); memcpy(Qy,y1,len); Qset = 1; }
	if ( nbits1&1 ) { memcpy (Qx,x1,len);  for ( i = 0 ; i < n ; i++ ) {_ff_neg(Qy[i],y1[i]); } Qset = 1; }
//if ( Pset ) printf ("P=(%ld,%ld)\n", _ff_get_ui(Px[0]),_ff_get_ui(Py[0])); else printf ("P=null\n");
//if ( Qset ) printf ("Q=(%ld,%ld)\n", _ff_get_ui(Qx[0]), _ff_get_ui(Qy[0])); else printf ("Q=null\n");
		
	memcpy (x,x1,len); memcpy (y,y1,len);
	parallel_double(x,y,n,f1);
	m = 2;
	max = ( pbits0 > pbits1 ? pbits0 : pbits1 );
	for (;;) {
//printf ("m=%lx(hex)  (%ld,%ld)\n", m, _ff_get_ui(x[0]), _ff_get_ui(y[0]));
		if ( (pbits0&m) ) { if ( ! Pset ) { memcpy(Px,x,len); memcpy(Py,y,len); Pset = 1; } else { parallel_add(Px,Py,x,y,n,f1); } }
		if ( (nbits0&m) ) { if ( ! Pset ) { memcpy(Px,x,len); for ( i = 0 ; i < n ; i++ ) {_ff_neg(Py[i],y[i]); } Pset = 1; } else { parallel_sub(Px,Py,x,y,n,f1); } }
		if ( (pbits1&m) ) { if ( ! Qset ) { memcpy(Qx,x,len); memcpy(Qy,y,len); Qset = 1; } else { parallel_add(Qx,Qy,x,y,n,f1); } }
		if ( (nbits1&m) ) { if ( ! Qset ) { memcpy(Qx,x,len); for ( i = 0 ; i < n ; i++ ) {_ff_neg(Qy[i],y[i]); } Qset = 1; } else { parallel_sub(Qx,Qy,x,y,n,f1); } }
//if ( Pset ) printf ("P=(%ld,%ld)\n", _ff_get_ui(Px[0]),_ff_get_ui(Py[0])); else printf ("P=null\n");
//if ( Qset ) printf ("Q=(%ld,%ld)\n", _ff_get_ui(Qx[0]), _ff_get_ui(Qy[0])); else printf ("Q=null\n");

		m <<= 1;
		if ( m > max ) break;
		parallel_double(x,y,n,f1);
	}
//printf("%d*(0,%ld)=(%ld,%ld)\n", pbits0-nbits0,_ff_get_ui(y1[0]),_ff_get_ui(Px[0]),_ff_get_ui(Py[0]));
//printf("%d*(0,%ld)=(%ld,%ld)\n", pbits1-nbits1,_ff_get_ui(y1[0]),_ff_get_ui(Qx[0]),_ff_get_ui(Qy[0]));
	if ( ! Qset ) {
		for ( i = 0 ; i < n ; i++ )  if ( _ff_zero(Px[i]) && _ff_zero(Py[i]) ) break;
	} else {
		for ( i = 0 ; i < n ; i++ ) {
			if ( ! _ff_zero(Px[i]) ) {
				if ( _ff_equal(Px[i],Qx[i]) ) break;
			} else {
				// don't count matches involving the identity, it likely just means the base point has small order
				if ( ! _ff_zero(Py[i]) && _ff_zero(Qx[i]) && ! _ff_zero(Qy[i]) )  break;
//printf ("ignored id match in scancurves\n");
			}
		}
	}
	if ( i == n ) return -1;
	*psign = ( _ff_equal(Py[i],Qy[i]) ? 1 : -1 );		// if y-values differ, points must be inverses
	return i;
}

// returns index of curve for which e*(x[i],y[i])=0 where e is represented as a signed chain based on a sliding window NAF with window size 3
// This is only marginally faster than scancurves (2 or 3 percent)
int scancurves1 (ff_t x1[], ff_t y1[], int n, int chain[FF_MAX_CHAIN_LEN], int chain_len, ff_t f1[])
{
	ff_t Px[MAX_CURVES], Py[MAX_CURVES], Qx[3][MAX_CURVES], Qy[3][MAX_CURVES];
	register int i,j,len;
//long e;

//printf("Curve y^2=x^3+%dx+b with base point P=(%ld,%ld) over F_%ld\n", _ff_get_ui(f1[0]), _ff_get_ui(x1[0]), _ff_get_ui(y1[0]), _ff_p);
	// Q[i] = (2i+1)(x,y)
	len = n*sizeof(ff_t);
	memcpy (Qx[0],x1,len); memcpy (Qy[0],y1,len);
	memcpy (Qx[1],x1,len); memcpy (Qy[1],y1,len);
//printf("P=(%ld,%ld)\n", _ff_get_ui(Qx[1][0]), _ff_get_ui(Qy[1][0]));
	parallel_double(Qx[1],Qy[1],n,f1);
//printf("2*P=(%ld,%ld)\n", _ff_get_ui(Qx[1][0]), _ff_get_ui(Qy[1][0]));
	memcpy (Qx[2],Qx[1],len); memcpy (Qy[2],Qy[1],len);
//	memcpy (Qx[3],Qx[1],len); memcpy (Qy[3],Qy[1],len);
	parallel_add(Qx[1],Qy[1],x1,y1,n,f1);
//printf("3*P=(%ld,%ld)\n", _ff_get_ui(Qx[1][0]), _ff_get_ui(Qy[1][0]));
	parallel_add(Qx[2],Qy[2],Qx[1],Qy[1],n,f1);
//printf("5*P=(%ld,%ld)\n", _ff_get_ui(Qx[2][0]), _ff_get_ui(Qy[2][0]));
//	parallel_add(Qx[3],Qy[3],Qx[2],Qy[2],n,f1);
//printf("7*P=(%ld,%ld)\n", _ff_get_ui(Qx[3][0]), _ff_get_ui(Qy[3][0]));

//printf("chain[0]=%d, chain[1]=%d\n", chain[0], chain[1]);
	memcpy (Px,Qx[chain[0]],len);  memcpy (Py,Qy[chain[0]],len);
//e = 2*chain[0]+1;
//printf("%ld*P=(%ld,%ld)\n", e, _ff_get_ui(Px[0]), _ff_get_ui(Py[0]));
	for ( j = 0 ; j < chain[1] ; j++ ) parallel_double(Px,Py,n,f1);
//e <<= chain[1];
//printf("%ld*P=(%ld,%ld)\n", e, _ff_get_ui(Px[0]),_ff_get_ui(Py[0]));
	
	for ( i = 2 ; i < chain_len ; i+=2 ) {
//printf("chain[%d]=%d, chain[%d]=%d\n", i, chain[i], i+1, chain[i+1]);
		if ( chain[i] >= 0 ) {
//puts("add");
			parallel_add(Px,Py,Qx[chain[i]],Qy[chain[i]],n,f1);
		} else {
//puts("sub");
			parallel_sub(Px,Py,Qx[-1-chain[i]],Qy[-1-chain[i]],n,f1);
		}
//e += 2*chain[i]+1;
//printf("%ld*P=(%ld,%ld)\n", e, _ff_get_ui(Px[0]), _ff_get_ui(Py[0]));
		for ( j = 0 ; j < chain[i+1] ; j++ ) parallel_double(Px,Py,n,f1);
//e <<=chain[i+1];
//printf("%ld*P=(%ld,%ld)\n", e, _ff_get_ui(Px[0]),_ff_get_ui(Py[0]));
	}
	for ( i = 0 ; i < n ; i++ ) if ( _ff_zero(Px[i]) && _ff_zero(Py[i]) ) break;
	if ( i == n ) return -1;
//printf ("Succeeded for curve %d\n", i);
	return i;
}

int compute_signed_exp_chain (int chain[FF_MAX_CHAIN_LEN], unsigned long e)
{
	unsigned long pbits, nbits;
	long x;
	int len;
	register int i, j, k;
	
	ui_NAF (&pbits, &nbits, e);
	len = ff_precompute_exp_chain (chain, pbits|nbits);
	for ( i = len-1, j=0 ; i > 0 ; i-= 2 ) {
		j += chain[i];
		k = (2*chain[i-1]+1) & (nbits>>j);
		if ( k&(~5) ) { printf ("bug in signed_exp_chain for e=%ld, k=%d\n", e, k); exit (0); }
		chain[i-1] -= k;
	}
	x = 2*chain[0]+1; x <<= chain[1];
	for ( i = 2 ; i < len ; i+= 2 )  {
		x += 2*chain[i]+1;
		x <<= chain[i+1];
	}
	if ( x != e ) { printf ("verification failed in signed_exp_chain for e=%ld\n", e); exit (0); }
	return len;
}

/*
	Fixing the base point with x-value zero speeds things up slightly (less than 5%), but isn't well suited to multi-exponentiation, so the code below is not used.
*/


//  adds (0,y1) to (x2,y2) given inv=1/x2 and assuming (0,y1) != (x2,y2)
static void inline affine_add_x0 (ff_t *x2, ff_t *y2, ff_t y1, ff_t inv)
{
	register ff_t t1, t2, t3;
	
	if ( _ff_zero(y1) ) return;											// (0,y1) is the identity
	if ( _ff_zero(*x2) ) { if ( _ff_zero(*y2) ) { _ff_set(*y2,y1); } else { _ff_set_zero(*y2); } return; }		// (0,y2) is either the identity, or the opposite of (0,y1)
	
	_ff_sub(t2,*y2,y1);  _ff_mult(t1,t2,inv);								// 1 = lambda = (y2-y1)/(x2-x1)
	_ff_square(t3,t1); _ff_subfrom(t3,*x2);								// t3 = x3 = lambda^2 - x1 - x2
	_ff_mult(t2,t1,t3);  _ff_addto(t2,y1);  _ff_neg(*y2,t2);					// y3 = (x1-x3)lambda - y1
	_ff_set(*x2,t3);													// x3 = t3
	// 3M + 4A
}

//  subtracts (0,y1) from (x2,y2) given inv=1/x2 and assuming (0,y1) != -(x2,y2)
static void inline affine_sub_x0 (ff_t *x2, ff_t *y2, ff_t y1, ff_t inv)
{
	register ff_t t1, t2, t3;
	
	if ( _ff_zero(y1) ) return;											// (0,y1) is the identity
	if ( _ff_zero(*x2) ) { if ( _ff_zero(*y2) ) { _ff_neg(*y2,y1); } else { _ff_set_zero(*y2); } return; }		// (0,y2) is either the identity, or equal to (0,y1)

	_ff_add(t2,*y2,y1);  _ff_mult(t1,t2,inv);								// 1 = lambda = (y2+y1)/(x2-x1)
	_ff_square(t3,t1); _ff_subfrom(t3,*x2);								// t3 = x3 = lambda^2 - x1 - x2
	_ff_mult(t2,t1,t3);  _ff_sub(*y2,y1,t2);								// y3 = (x1-x3)lambda + y1
	_ff_set(*x2,t3);													// x3 = t3
	// 3M +  3A
}

// sets (x[i],y[i]) = (x[i],y[i]) + (0,y1[i]) for i = 0 to n-1
static void inline parallel_add_x0 (ff_t x[], ff_t y[], ff_t y1[], int n, ff_t f1)
{
	ff_t s[MAX_CURVES];
	ff_t t[MAX_CURVES+1];
	register ff_t t0, t1;
	register int i;

	_ff_set_one(t[0]);
	for ( i = 0 ; i < n ; i++ ) {
		if ( _ff_zero(x[i]) ) {
			if ( ! _ff_zero(y[i]) && _ff_equal(y[i],y1[i]) ) { _ff_add(s[i], y[i], y[i]); } else { _ff_set_one(s[i]); }
		} else {
			_ff_set(s[i], x[i]);											// probably could avoid this copy if we were clever
		}
		_ff_mult(t[i+1], t[i], s[i]);
	}
	_ff_invert (t1, t[n]);
	for ( i = n-1 ; i > 0 ; i-- ) {
		_ff_mult (t0, t[i], t1);
		_ff_mult (t1, t1, s[i]);
		// one could play with the order of these tests to try to improve branch prediction, but it's not obvious this would help
		if ( ! _ff_zero(x[i]) ) { affine_add_x0 (x+i,y+i,y1[i],t0); continue; }
		if ( _ff_zero(y[i]) ) { _ff_set(y[i],y1[i]); continue; }
		if ( ! _ff_equal(y[i],y1[i]) ) { _ff_set_zero(y[i]);  continue; }
		affine_double_x0 (x+i,y+i,y[i],t0, f1);
	}
gop_counter += n;
	if ( ! _ff_zero(x[i]) ) { affine_add_x0 (x+i,y+i,y1[i],t1); return; }
	if ( _ff_zero(y[i]) ) { _ff_set(y[i],y1[i]); return; }
	if ( ! _ff_equal(y[i],y1[i]) ) { _ff_set_zero(y[i]);  return; }
	affine_double_x0 (x+i,y+i,y[i],t1, f1);
	// total cost is approximately (6M+4A)*n + I - 2M
}
