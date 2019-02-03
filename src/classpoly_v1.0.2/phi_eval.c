#include "ff_poly.h"
#include "phi_poly.h"
#include "phi_eval.h"

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

void phi_eval_ff (ff_t f[], ff_t phi[], int m, ff_t J)
{
	switch ( phi_sparse_factor(phi) ) {
	case 2:
		switch (m) {
		case 3: phi3_s2_eval_ff(f,phi,J); return;
		case 5: phi5_s2_eval_ff(f,phi,J); return;
		case 7: phi7_s2_eval_ff(f,phi,J); return;
		case 11: phi11_s2_eval_ff(f,phi,J); return;
		case 13: phi13_s2_eval_ff(f,phi,J); return;
		case 17: phi17_s2_eval_ff(f,phi,J); return;
		case 19: phi19_s2_eval_ff(f,phi,J); return;
		}
	case 3:
		switch (m) {
		case 2: phi2_s3_eval_ff(f,phi,J); return;
		case 5: phi5_s3_eval_ff(f,phi,J); return;
		case 7: phi7_s3_eval_ff(f,phi,J); return;
		case 11: phi11_s3_eval_ff(f,phi,J); return;
		case 13: phi13_s3_eval_ff(f,phi,J); return;
		case 17: phi17_s3_eval_ff(f,phi,J); return;
		case 19: phi19_s3_eval_ff(f,phi,J); return;
		}
	case 4:
		switch (m) {
		case 3: phi3_s4_eval_ff(f,phi,J); return;
		case 5: phi5_s4_eval_ff(f,phi,J); return;
		case 7: phi7_s4_eval_ff(f,phi,J); return;
		case 11: phi11_s4_eval_ff(f,phi,J); return;
		case 13: phi13_s4_eval_ff(f,phi,J); return;
		case 17: phi17_s4_eval_ff(f,phi,J); return;
		case 19: phi19_s4_eval_ff(f,phi,J); return;
		}
	case 6:
		switch (m) {
		case 5: phi5_s6_eval_ff(f,phi,J); return;
		case 7: phi7_s6_eval_ff(f,phi,J); return;
		case 11: phi11_s6_eval_ff(f,phi,J); return;
		case 13: phi13_s6_eval_ff(f,phi,J); return;
		case 17: phi17_s6_eval_ff(f,phi,J); return;
		case 19: phi19_s6_eval_ff(f,phi,J); return;
		}
	case 8:
		switch (m) {
		case 3: phi3_s8_eval_ff(f,phi,J); return;
		case 5: phi5_s8_eval_ff(f,phi,J); return;
		case 7: phi7_s8_eval_ff(f,phi,J); return;
		case 11: phi11_s8_eval_ff(f,phi,J); return;
		case 13: phi13_s8_eval_ff(f,phi,J); return;
		case 17: phi17_s8_eval_ff(f,phi,J); return;
		case 19: phi19_s8_eval_ff(f,phi,J); return;
		}
	case 12:
		switch (m) {
		case 5: phi5_s12_eval_ff(f,phi,J); return;
		case 7: phi7_s12_eval_ff(f,phi,J); return;
		case 11: phi11_s12_eval_ff(f,phi,J); return;
		case 13: phi13_s12_eval_ff(f,phi,J); return;
		case 17: phi17_s12_eval_ff(f,phi,J); return;
		case 19: phi19_s12_eval_ff(f,phi,J); return;
		}
	case 24:
		switch (m) {
		case 5: phi5_s24_eval_ff(f,phi,J); return;
		case 7: phi7_s24_eval_ff(f,phi,J); return;
		case 11: phi11_s24_eval_ff(f,phi,J); return;
		case 13: phi13_s24_eval_ff(f,phi,J); return;
		case 17: phi17_s24_eval_ff(f,phi,J); return;
		case 19: phi19_s24_eval_ff(f,phi,J); return;
		}
	default:
		switch (m) {
		case 2: phi2_eval_ff(f,phi,J); return;
		case 3: phi3_eval_ff(f,phi,J); return;
		case 5: phi5_eval_ff(f,phi,J); return;
		case 7: phi7_eval_ff(f,phi,J); return;
		case 11: phi11_eval_ff(f,phi,J); return;
		case 13: phi13_eval_ff(f,phi,J); return;
		case 17: phi17_eval_ff(f,phi,J); return;
		case 19: phi19_eval_ff(f,phi,J); return;
		}
	}
	if ( phi_sparse_factor(phi) > 1 ) _phi_s_eval_ff(f,phi,m,J); else _phi_eval_ff(f,phi,m,J);
}


/*
	Returns f(x)=b^(m+1) * phi_m (x,a/b) as a poly in x of degree m+1 with coefficients in the current finite field.

	This function is not currently used, all the small cases where it is worthwhile are handled in phi_eval.h
*/

void phi_qeval_ff (ff_t f[], ff_t phi[], int m, ff_t a, ff_t b)
{
	ff_t ya[PHI_MAX_M+2];
	ff_t yb[PHI_MAX_M+2];
	register int i;
			
	_ff_set_one(ya[0]); _ff_set (ya[1],a); for ( i = 2 ; i <= m+1 ; i++ ) ff_mult(ya[i],ya[i-1],a);		// Set ya[i]=a^i for i from 0 to m+1
	_ff_set (yb[m],b); for ( i = m-1 ; i>=0 ; i-- ) ff_mult(yb[i],yb[i+1],b);						// Set yb[i]=b^(m+1-i) for i from 0 to m
	for ( i = 0 ; i <= m ; i++ ) ff_mult(ya[i],ya[i],yb[i]);
	for ( i = 0 ; i <= m ; i++ ) ff_dot_product(f+i,ya,phi+(m+1)*i,m+1);
	_ff_set(f[m+1],ya[0]); _ff_addto(f[0],ya[m+1]); 
}
