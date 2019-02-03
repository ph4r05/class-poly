#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "ff_poly.h"
#include "ff_poly/ffpolysmall.h"
#include "class_inv.h"
#include "bipoly.h"
#include "phi_poly.h"
#include "phi_fj_strings.h"

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

#define PHI_FJ_CACHESIZE		2

// ordered from best to worst
int atkin_enum_levels[15] = { 71, 59, 47, 41, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3, 0 };
int single_eta_enum_levels[5] = { 13, 7, 5, 3, 0 };
int double_eta_enum_levels[10] = { 39, 35, 26, 21, 15, 14, 10, 9, 6, 0 };
int double_eta_enum_heights[10] = { 28, 24, 21, 16, 12, 12, 9, 6, 6, 0 };

static struct phi_fj_cache_entry {
	bipoly_mpz_t *Phi_fj;
	bipoly_ff_t *phi_fj, *phi_fji;
	int phi_fj_inv;
	int phi_fj_terms;
	int phi_fj_fdeg;
	int phi_fj_jdeg;
} Phi_fj_cache[PHI_FJ_CACHESIZE];

bipoly_mpz_t *Phi_fj;
static bipoly_ff_t *phi_fj, *phi_fji;
static ff_t phi_fj_redp, phi_fji_redp;
static int phi_fj_srt, phi_fji_srt;
static int Phi_fj_lru;

int phi_fj_inv, phi_fj_terms, phi_fj_fdeg, phi_fj_jdeg;

/*
	Pick a class invariant.  For the moment, restrict to cases where we can enumerate.
	We avoid calling inv_height_factor for non-specific invariants here because we don't want to have to load any phi_fj polys.
*/
int inv_pick_invariant (long D)
{
	register int inv, best_inv, best_hf, hf, i, N;
	
	best_inv = 0;  best_hf = 1;
	if ( inv_good_discriminant(D,INV_F) ) { best_inv = INV_F; best_hf = 72; }
	for ( i = 0 ; (N=atkin_enum_levels[i]) && best_hf < 4*(N+1) ; i++ ) {
		inv = INV_ATKIN+N;
		if ( inv_ramified(D,inv) && inv_good_discriminant(D,inv) ) { best_inv = inv; best_hf = 2*(N+1); }
		if ( best_hf > (N+1)/2 ) continue;
		if ( inv_good_discriminant(D,inv) ) { best_inv = inv;  best_hf = (N+1)/2; }
	}
	for ( inv = 1 ; inv < INV_MAX_SPECIFIC_ENUM ; inv++ ) {
		if ( ! inv_class_invariant(inv) ) continue;
		if ( ! inv_enum(inv) ) continue;
		hf = inv_height_factor(inv);
		if ( inv_ramified(D,inv) && (inv_atkin(inv)||inv_double_eta(inv)) ) hf *= 4;
		if ( hf <= best_hf ) continue;
		if ( inv_good_discriminant(D,inv) ) { best_inv = inv; best_hf = hf; }
	}
	for ( i = 0 ; (N=double_eta_enum_levels[i]) && best_hf < 2*double_eta_enum_levels[i] ; i++ ) {
		hf = double_eta_enum_heights[i];
		inv = INV_DOUBLE_ETA+N;
		if ( inv_ramified(D,inv) ) {
			if ( ! inv_good_discriminant(D,inv) ) continue;
			if ( 4*hf > best_hf ) { best_inv = inv; best_hf = 4*hf; }
			break;
		}
		if ( hf <= best_hf ) continue;
		if ( inv_good_discriminant(D,inv) ) { best_inv = inv;  best_hf = hf; }
	}
	for ( i = 0 ; (N=single_eta_enum_levels[i]) && best_hf < 2*(N+1); i++ ) {
		inv  = INV_SINGLE_ETA+N;
		if ( best_hf >= (N+1) ) break;
		if ( inv_good_discriminant(D,inv) ) { best_inv = inv; best_hf = (N+1); }
	}
	return best_inv;
}


// note we assume here we have already got the right invariant
static inline void inv_reduce_phi_fj (int inv, int srt)
{
	inv_load_phi_fj(inv);
	if ( _ff_p != phi_fj_redp ) { bipoly_reduce_ff (phi_fj, Phi_fj, phi_fj_terms); bipoly_sort_ff (phi_fj, phi_fj_terms, srt);  phi_fj_redp = _ff_p;  phi_fj_srt = srt; }
	if ( srt != phi_fj_srt ) { bipoly_sort_ff (phi_fj, phi_fj_terms, srt);  phi_fj_srt = srt; }
}

// We always involute the index 0 coefficient
static inline void inv_reduce_phi_fji (int inv, int srt)
{
	ff_t Ns[PHI_SINGLE_ETA_MAX_FDEG+1], t0;
	int N, s, i, p1, p2;
	
	inv_reduce_phi_fj (inv, srt);
	if ( _ff_p != phi_fji_redp ) {
		N = inv_degree(&p1,&p2,inv); s = inv_s(N);  _ff_set_ui(t0,N);
		_ff_set_one(Ns[0]); ff_exp_small(Ns+1,&t0,s);
		for ( i = 2 ; i <= PHI_SINGLE_ETA_MAX_FDEG ; i++ ) _ff_mult(Ns[i],Ns[i-1],Ns[1]);
		for ( i = 0 ; i < phi_fj_terms ; i++ ) {
			if ( phi_fj[i].e[0] > 0 ) _ff_mult(phi_fji[i].c,phi_fj[i].c,Ns[phi_fj[i].e[0]-1]); else _ff_set_one(phi_fji[i].c);											// We assume the constant coefficient in X is N^s
			phi_fji[i].e[0] = N+1-phi_fj[i].e[0];
			phi_fji[i].e[1] = phi_fj[i].e[1];
		}
		bipoly_sort_ff (phi_fji, phi_fj_terms, srt);
		phi_fji_srt = srt;
		phi_fji_redp = _ff_p;
	}
	if ( srt != phi_fji_srt ) { bipoly_sort_ff (phi_fji, phi_fj_terms, srt);  phi_fji_srt = srt; }
}


void _ff_inv_from_j (ff_t *x, ff_t *j, int inv)
{
	switch ( inv ) {
	case INV_J: _ff_set(*x,*j); break;
	case INV_G2: ff_inv_gamma2_from_j(x,j); break;
	case INV_F8: ff_inv_f8_from_j(x,j); break;
	case INV_F4: ff_inv_f4_from_j(x,j); break;
	case INV_F2: ff_inv_f2_from_j(x,j); break;
	case INV_F: ff_inv_f_from_j(x,j); break;
	case INV_T6: ff_inv_t6_from_j(x,j); break;
	case INV_T2: ff_inv_t2_from_j(x,j); break;
	case INV_F3: ff_inv_f3_from_j(x,j); break;
	case INV_T: ff_inv_t_from_j(x,j); break;
	case INV_U8: ff_inv_u8_from_j(x,j); break;
	case INV_U2: ff_inv_u2_from_j(x,j); break;
	case INV_U: ff_inv_u_from_j(x,j); break;
	case INV_W3E2: case INV_W5E2: case INV_W7E2: ff_inv_single_eta_from_j(x,j,inv); break;
	case INV_W3W5: case INV_W2W3E2: case INV_W2W5E2: case INV_W2W7E2: case INV_W3W3E2: case INV_W3W3: case INV_W5W5:
	case INV_W2W3: case INV_W2W5: case INV_W2W7: case INV_W2W11E2: case INV_W2W13: case INV_W3W7: case INV_W3W11E2: ff_inv_double_eta_from_j(x,j,inv); break;
	default:
		if ( inv >= INV_ATKIN && inv < INV_ATKIN_END ) { ff_inv_atkin_from_j(x,j,inv); break; }
		if ( inv >= INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END ) { ff_inv_single_eta_from_j(x,j,inv); break; }
		if ( inv >= INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) { ff_inv_double_eta_from_j(x,j,inv); break; }
		printf("Don't know how to handle invariant %d in _ff_inv_from_j\n", inv); abort ();
	}
}

void _ff_j_from_inv (ff_t *j, ff_t *x, int inv)
{
	switch (inv) {
	case INV_J: _ff_set(*j,*x); break;
	case INV_G2: ff_inv_j_from_gamma2 (j,x); break;
	case INV_F: ff_inv_j_from_f (j,x); break;
	case INV_F8: ff_inv_j_from_f8 (j,x); break;
	case INV_F4: ff_inv_j_from_f4 (j,x); break;
	case INV_F2: ff_inv_j_from_f2 (j,x); break;
	case INV_T6: ff_inv_j_from_t6 (j,x); break;
	case INV_T2: ff_inv_j_from_t2 (j,x); break;
	case INV_F3: ff_inv_j_from_f3 (j,x); break;
	case INV_T: ff_inv_j_from_t (j,x); break;
	case INV_U8: ff_inv_j_from_u8 (j,x); break;
	case INV_U2: ff_inv_j_from_u2 (j,x); break;
	case INV_U: ff_inv_j_from_u (j,x); break;
	case INV_W3E2: case INV_W5E2: case INV_W7E2: ff_inv_j_from_single_eta (j,x,inv); break;
	case INV_W3W5: case INV_W2W3E2: case INV_W2W5E2: case INV_W2W7E2: case INV_W3W3E2: case INV_W3W3: case INV_W5W5:
	case INV_W2W3: case INV_W2W5: case INV_W2W7: case INV_W2W11E2: case INV_W2W13:  case INV_W2W17: case INV_W3W7:case INV_W3W11E2:  ff_inv_j_from_double_eta (j,x,inv); break;
	
	default:
		if ( inv >= INV_ATKIN && inv < INV_ATKIN_END ) { ff_inv_j_from_atkin (j,x,inv);break; }
		if ( inv >= INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END ) { ff_inv_j_from_single_eta (j,x,inv);break; }
		if ( inv >= INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) { ff_inv_j_from_double_eta (j,x,inv);break; }
		printf("Don't know how to handle invariant %d in ff_j_from_inv\n", inv); abort ();
	}
}

void ff_2j_from_inv (ff_t j[2], ff_t *x, int inv)
{
	switch (inv) {	
	case INV_W3E2: case INV_W5E2: case INV_W7E2: ff_inv_2j_from_single_eta (j,x,inv); break;
	case INV_W3W5: case INV_W2W3E2: case INV_W2W5E2: case INV_W2W7E2:  case INV_W3W3E2: case INV_W3W3: case INV_W5W5:
	case INV_W2W3: case INV_W2W5: case INV_W2W7: case INV_W2W11E2: case INV_W2W13: case INV_W2W17:  case INV_W3W7: case INV_W3W11E2: ff_inv_2j_from_double_eta(j,x,inv); break;
	default:
		if ( inv >= INV_ATKIN && inv < INV_ATKIN_END ) { ff_inv_2j_from_atkin (j,x,inv);break; }
		if ( inv >= INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END ) { ff_inv_2j_from_single_eta (j,x,inv);break; }
		if ( inv >= INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) { ff_inv_2j_from_double_eta (j,x,inv);break; }
		printf("Don't know how to handle invariant %d in ff_2j_from_inv \n", inv); abort ();
	}
}

int ff_multi_j_from_inv (ff_t j[], ff_t *x, int inv, int n, int desired)
{
	ff_t f[PHI_MAX_JDEG+1];
	int d_f;

	inv_load_phi_fj (inv);
	if ( phi_fj_jdeg > n ) return -phi_fj_jdeg;
	if ( phi_fj_jdeg == 1 ) { _ff_j_from_inv (j, x, inv); return 1; }
	if ( desired == 1 && phi_fj_jdeg == 2 ) { _ff_j_from_inv (j, x, inv); return 1; }	// we assume that for j degree 2 we can always get a good j
	inv_reduce_phi_fj (inv,1);
	bipoly_eval_ff (f, phi_fj_jdeg, phi_fj, phi_fj_terms, 0, *x);
	d_f = ff_poly_degree(f,phi_fj_jdeg);
	ff_poly_monic (f,&d_f,f,d_f);
	return ff_poly_distinct_roots (j, f, d_f);
}

int ff_inv_from_2j (ff_t x[2], ff_t *j1, ff_t *j2, int inv)
{
	switch (inv) {
	case INV_W3E2: case INV_W5E2: case INV_W7E2: return ff_inv_single_eta_from_2j (x,j1,j2,inv);
	case INV_W3W5: case INV_W2W3E2: case INV_W2W5E2: case INV_W2W7E2:  case INV_W3W3E2: case INV_W3W3: case INV_W5W5:
	case INV_W2W3: case INV_W2W5: case INV_W2W7: case INV_W2W11E2: case INV_W2W13: case INV_W2W17:  case INV_W3W7: case INV_W3W11E2:  return ff_inv_double_eta_from_2j(x,j1,j2,inv);
	default:
		if ( inv >= INV_ATKIN && inv < INV_ATKIN_END ) return ff_inv_atkin_from_2j (x,j1,j2,inv);
		if ( inv >= INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END ) return ff_inv_single_eta_from_2j (x,j1,j2,inv);
		if ( inv >= INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) return ff_inv_double_eta_from_2j (x,j1,j2,inv);
	}
	printf("Don't know how to handle invariant %d in ff_inv_from_2j \n", inv); abort ();
}

int ff_j_from_2inv (ff_t *j, ff_t *x1, ff_t *x2, int inv)
{
	switch (inv) {
	case INV_W3W5: case INV_W2W3E2: case INV_W2W5E2: case INV_W2W7E2:  case INV_W3W3E2: case INV_W3W3: case INV_W5W5:
	case INV_W2W3: case INV_W2W5: case INV_W2W7: case INV_W2W11E2: case INV_W2W13:  case INV_W2W17: case INV_W3W7: case INV_W3W11E2:  return ff_inv_j_from_2double_eta (j,x1,x2,inv);
	default:
		if ( inv >= INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) return ff_inv_j_from_2double_eta (j,x1,x2,inv);
	}
	printf("Don't know how to handle invariant %d in ff_inv_from_2j \n", inv); abort ();
}

// requires _ff_p=2 mod 3
void ff_inv_f8_from_j (ff_t *f8, ff_t *J)
{
	ff_t tf[4], r[3];
	ff_t t0;
	register int i, n;
	
	if ( _ff_p1mod3 ) { printf ("p=%ld is not 2 mod 3 in ff_inv_f8_from_j\n", _ff_p); abort (); }
	ff_cbrt(&t0,J);				// we know the cube root (gamma_2) exists and is unique
	_ff_set_one(tf[3]); _ff_set_zero(tf[2]); _ff_neg(tf[1],t0); _ff_set_ui(tf[0],16); ff_negate(tf[0]);		// X^3 - gamma_2*X - 16
	n = ff_poly_distinct_roots(r, tf, 3);			// we know that f^8 is a root of tf
	for ( i = 0 ; i < n ; i++ ) if ( ff_residue(r[i]) ) break;	// if  r[i] is a residue, it has an eighth root
	if ( i == n ) { printf ("Unable to determine f^8 from j=%ld over F_%ld, numroots=%d\n", _ff_get_ui(*J), _ff_p, n); abort (); }
	_ff_set(*f8,r[i]);					
}

// requires _ff_p=2 mod 3 and _ff_p=3 mod 4
void ff_inv_f_from_j (ff_t *f, ff_t *J)
{
	ff_t tf[4], tr[3];
	ff_t t0;
	register int i, n;
	
	if ( _ff_p1mod3 || (_ff_p&3)==1) { printf ("p=%ld is not 2 mod 3 or is not 3 mod 4 in ff_inv_f_from_j\n", _ff_p); abort (); }
	ff_cbrt(&t0,J);							// we know the cube root (gamma_2) exists and is unique
	_ff_set_one(tf[3]); _ff_set_zero(tf[2]); _ff_neg(tf[1],t0); _ff_set_ui(tf[0],16); ff_negate(tf[0]);
	n = ff_poly_distinct_roots(tr, tf, 3);			// we know that f^8 is a root of tf = X^3-gamma_2*X-16

	for ( i = 0 ; i < n ; i++ )  if ( ff_eighth_root(f,tr+i) ) break;
	if ( i == n ) { printf ("Unable to determine f from j=%ld over F_%ld, numroots=%d\n", _ff_get_ui(*J), _ff_p, n); abort (); }
}

// requires _ff_p=2 mod 3 and _ff_p=3 mod 4
void ff_inv_h2_from_j (ff_t *f, ff_t *J)
{
	ff_t tf[4], tr[3];
	ff_t t0;
	register int i, n;
	
	if ( _ff_p1mod3 || (_ff_p&3)==1) { printf ("p=%ld is not 2 mod 3 or is not 3 mod 4 in ff_inv_h2_from_j\n", _ff_p); abort (); }
	ff_cbrt(&t0,J);							// we know the cube root (gamma_2) exists and is unique
	ff_mult (t0, t0, _ff_fourth);  ff_mult (t0, t0, _ff_fourth);  ff_mult (t0, t0, _ff_fourth);  	// set t0 = gamma_2/256
	_ff_set_one(tf[3]); _ff_set_zero(tf[2]); _ff_neg(tf[1],t0); _ff_set_ui(tf[0],_ff_fourth);
	n = ff_poly_distinct_roots(tr, tf, 3);			// we know that h^8 is a root of tf = X^3-gamma_2/256*X+1/4

	for ( i = 0 ; i < n ; i++ )  if ( ff_fourth_root(f,tr+i) ) break;
	if ( i == n ) { printf ("Unable to determine f from j=%ld over F_%ld, numroots=%d\n", _ff_get_ui(*J), _ff_p, n); abort (); }
}

// requires _ff_p=2 mod 3
void ff_inv_t6_from_j (ff_t *t6, ff_t *J)
{
	ff_t tf[3], tr[2];
	ff_t t0, t1;
	register int i, n;
	
	if ( _ff_p1mod3 ) { printf ("p=%ld is not 2 mod 3 In ff_inv_t6_from_j\n", _ff_p); abort (); }
	ff_cbrt(&t0,J);							// we know the cube root (gamma_2) exists and is unique
	_ff_set_ui(t1,6);
	_ff_addto(t0,t1);
	_ff_set_one(tf[2]); _ff_neg(tf[1],t0); _ff_set_ui(tf[0],27);	ff_negate(tf[0]);
	n = ff_poly_distinct_roots(tr, tf, 2);			// we know that t^6 is a root of tf
	for ( i = 0 ; i < n ; i++ ) if ( ff_residue(tr[i]) ) break;	// if it has a sqrt it has a sixth root
	if ( i == n ) { printf ("Unable to determine t^6 from j=%ld over F_%ld, numroots=%d\n", _ff_get_ui(*J), _ff_p, n); ff_poly_print(tf,2); abort (); }
	_ff_set(*t6,tr[i]);			
}

// requires _ff_p=2 mod 3
void ff_inv_t_from_j (ff_t *t, ff_t *J)
{
	ff_t tf[3], tr[2];
	ff_t t0, t1;
	register int i, n;
	
	if ( _ff_p1mod3 ) { printf ("p=%ld is not 2 mod 3 In ff_inv_t6_from_j\n", _ff_p); abort (); }
	ff_cbrt(&t0,J);							// we know the cube root (gamma_2) exists and is unique
	_ff_set_ui(t1,6);
	_ff_addto(t0,t1);
	_ff_set_one(tf[2]); _ff_neg(tf[1],t0); _ff_set_ui(tf[0],27);	ff_negate(tf[0]);
	n = ff_poly_distinct_roots(tr, tf, 2);			// we know that t^6 is a root of tf
	for ( i = 0 ; i < n ; i++ ) if ( ff_sixth_root(&t0,tr+i) ) break;
	// t0 is a solution to x^12-(\gamma_2+6)*x^6-27, a candidate value for t
	if ( i == n ) { printf ("Unable to determine t^2 from j=%ld over F_%ld, numroots=%d\n", _ff_get_ui(*J), _ff_p, n);ff_poly_print(tf,2);  abort (); }
	_ff_set(*t,t0);
}

// Requires p = 3 mod 4.  Always returns a residue (this relies on the implementation of ff_eighth root  via exponentiation for p=3 mod 4).
// Note that phi_fj(X,j) should always have 6 roots, r, 1-r, 1/r, 1-1/r, 1/(1-r),r/(r-1), exactly 4 of which are residues (only one of 1-r and 1-1/r can be)
// Thus the 4 roots are either r, 1-r, 1/r, and 1/(1-r)  or r, 1-1/r, 1/r, r/(r-1).
void ff_inv_u_from_j (ff_t *u, ff_t *j)
{
	ff_t f[7], r[6];
	ff_t t0, t1, t2;
	register int i, n;

	// Borwein+Borwein p. 115:  256(1-lambda+lambda^2)^3 = lambda^2(1-lambda^2)j
	// with X=lambda, j=Y, we get:  X^6 - 3*X^5 + (6 - Y/256)*X^4 + (Y/128- 7)*X^3 - (6 - Y/256)X^2 - 3*X + 1 = 0
	_ff_square(t0,_ff_fourth);  _ff_square(t2,t0); _ff_mult(t0,t2,*j); _ff_set_ui(t1,6); _ff_subfrom(t1,t0);	// t1 = 6 - j/256
	_ff_x2(t0); _ff_set_ui(t2,7);  _ff_subfrom(t0,t2);	// t0 = j/128 - 7
	 _ff_set_ui(t2,3); ff_negate(t2);  // t2 = -3;
	_ff_set_one(f[6]);_ff_set(f[5],t2); _ff_set(f[4],t1); _ff_set(f[3],t0); _ff_set(f[2],t1); _ff_set(f[1],t2);  _ff_set_one(f[0]);
	n = ff_poly_distinct_roots(r,f,6);
	for ( i = 0 ; i < n ; i++ )  if ( ff_eighth_root(u,r+i) ) break;
	if ( i == n ) { printf ("Unable to determine u from j=%ld over F_%ld, numroots=%d\n", _ff_get_ui(*j), _ff_p, n); abort (); }
//printf ("Converted j=%ld to lambda=%ld, u=%ld\n", _ff_get_ui(*j), _ff_get_ui(r[i]), _ff_get_ui(*u));
}


void _inv_free_phi_fj (int i)
{
	struct phi_fj_cache_entry *x;
	
	if ( i < 0 || i > PHI_FJ_CACHESIZE ) { err_printf ("Inavlid phi_fj cache index %d\n", i); abort (); }
	x = Phi_fj_cache+i;
	if ( x->Phi_fj ) bipoly_free_mpz(x->Phi_fj, x->phi_fj_terms); x->Phi_fj = 0;
	if ( x->phi_fj ) free (x->phi_fj); x->phi_fj = 0;
	if ( x->phi_fji ) free (x->phi_fji); x->phi_fji = 0;
	x->phi_fj_inv = 0;  x->phi_fj_terms = x->phi_fj_fdeg = x->phi_fj_jdeg = 0;
}


void _inv_load_phi_fj (int inv)
{
	char filename[256];
	register int i, k;

	if ( !inv || phi_fj_inv == inv ) return;
	phi_fj_redp = phi_fji_redp = 0;
	for ( i = 0 ; i < PHI_FJ_CACHESIZE ; i++ ) if ( Phi_fj_cache[i].phi_fj_inv == inv ) break;
	if ( i < PHI_FJ_CACHESIZE ) {
		phi_fj_inv = Phi_fj_cache[i].phi_fj_inv;
		Phi_fj = Phi_fj_cache[i].Phi_fj;
		phi_fj = Phi_fj_cache[i].phi_fj;
		phi_fji = Phi_fj_cache[i].phi_fji;
		phi_fj_terms = Phi_fj_cache[i].phi_fj_terms;
		phi_fj_fdeg = Phi_fj_cache[i].phi_fj_fdeg;
		phi_fj_jdeg = Phi_fj_cache[i].phi_fj_jdeg;
		Phi_fj_lru = (i+1) % PHI_FJ_CACHESIZE;		// only true for cachsize of 2, but works in any case
//dbg_printf ("Loaded phi_fj (%d,%d)  for inv=%d from cache\n", phi_fj_fdeg, phi_fj_jdeg, inv);
		return;
	}
	for ( i = 0 ; i < PHI_FJ_CACHESIZE ; i++ ) if ( ! Phi_fj_cache[i].phi_fj_inv ) break;
	if ( i == PHI_FJ_CACHESIZE ) { _inv_free_phi_fj (Phi_fj_lru); i = Phi_fj_lru; }
	Phi_fj_cache[i].phi_fj_inv = phi_fj_inv = inv;
	phi_fj_redp = phi_fji_redp = phi_fj_fdeg = 0;
	// handle hard wired cases
	switch ( inv ) {
	case INV_F:	phi_fj_fdeg = 72; phi_fj_jdeg = 1; break;
	case INV_G2:   phi_fj_fdeg = 3; phi_fj_jdeg = 1; break;
	case INV_F8:	phi_fj_fdeg = 9; phi_fj_jdeg = 1; break;
	case INV_F4:   phi_fj_fdeg = 18; phi_fj_jdeg = 1;break;
	case INV_F2:   phi_fj_fdeg = 36; phi_fj_jdeg = 1; break;
	case INV_T6:   phi_fj_fdeg = 6; phi_fj_jdeg = 1; break;
	case INV_T2:   phi_fj_fdeg = 18; phi_fj_jdeg = 1; break;
	case INV_F3:   phi_fj_fdeg = 24; phi_fj_jdeg = 1; break;
	case INV_T:	phi_fj_fdeg = 36; phi_fj_jdeg = 1; break;
	case INV_U8:	phi_fj_fdeg = 6; phi_fj_jdeg = 1; break;
	case INV_U2:	phi_fj_fdeg = 24; phi_fj_jdeg = 1; break;
	case INV_U:	phi_fj_fdeg = 48; phi_fj_jdeg = 1; break;
	}
	if ( phi_fj_fdeg ) { Phi_fj_cache[i].phi_fj_fdeg = phi_fj_fdeg;  Phi_fj_cache[i].phi_fj_jdeg = phi_fj_jdeg; return; }

	sprintf(filename, "%s/phi_%s_j.txt", PHI_DIR, inv_string(inv));
	phi_fj_terms = bipoly_load_mpz (&Phi_fj, filename);
	if ( phi_fj_terms <= 0 ) { err_printf ("Unable to load Phi_fj poly from file %s for invariant %d\n", filename, inv); abort (); }
	bipoly_sort_mpz (Phi_fj, phi_fj_terms, 1);				// keep sorted on j, for evaluation on f
	phi_fj_fdeg = phi_fj_jdeg = 0;
	for ( k = 0 ; k < phi_fj_terms ; k++ ) {
		if ( Phi_fj[k].e[0] > phi_fj_fdeg ) phi_fj_fdeg = Phi_fj[k].e[0];
		if ( Phi_fj[k].e[1] > phi_fj_jdeg ) phi_fj_jdeg = Phi_fj[k].e[1];
	}
	phi_fj = malloc(sizeof(*phi_fj)*phi_fj_terms);	// allocated but do note reduce yet
	if ( inv_single_eta(inv) ) phi_fji = malloc(sizeof(*phi_fj)*phi_fj_terms); else phi_fji = 0;
	Phi_fj_cache[i].Phi_fj = Phi_fj;
	Phi_fj_cache[i].phi_fj = phi_fj;
	Phi_fj_cache[i].phi_fji = phi_fji;
	Phi_fj_cache[i].phi_fj_terms = phi_fj_terms;
	Phi_fj_cache[i].phi_fj_fdeg = phi_fj_fdeg;
	Phi_fj_cache[i].phi_fj_jdeg = phi_fj_jdeg;
dbg_printf ("Loaded phi_fj (%d,%d)  for inv=%d from file %s\n", phi_fj_fdeg, phi_fj_jdeg, inv, filename);
}

static inline void ff_inv_power (ff_t *x, ff_t *w, int inv)
{
	register ff_t t;

	switch (inv) {
	case INV_W2W3:  _ff_square(t,*w); _ff_mult(t,t,*w); ff_square(t,t); _ff_square(*x,t); return;																	// twelfth
	case INV_W3E2: case INV_W2W11E2: case INV_W2W3E2: case INV_W2W5: case INV_W3W3: case INV_W3W11E2:  _ff_square(t,*w); ff_mult(t,t,*w); _ff_square(*x,t); return;	// sixth
	case INV_W2W7: _ff_square(t,*w); _ff_square(*x,t); return;																							// fourth
	case INV_W5E2:  case INV_W3W5: case INV_W2W5E2: case INV_W3W3E2: case INV_W5W5:  case INV_W2W17:  _ff_square(t,*w); ff_mult(*x,t,*w); return;				// cube
	case INV_W7E2: case INV_W2W7E2: case INV_W2W13: case INV_W3W7:  _ff_square(*x,*w); return;															// square
	default:  _ff_set(*x,*w); return;
	}
}

static inline int ff_inv_root (ff_t *x, ff_t *w, int inv)
{
	ff_t t;
	
	// don't bother checking for cube root failure, we assume p=2 mod 3 in these cases
	// We could compute sixth roots more efficiently (especially if p=11 mod 12) but don't worry about this now.
	switch (inv ) {
	case INV_W2W3:  ff_cbrt(&t,w); if ( ! ff_sqrt(&t,&t) ) return 0; if ( ! ff_sqrt(x,&t) ) return 0; else return 1;										// twelfth
	case INV_W3E2: case INV_W2W3E2: case INV_W2W5: case INV_W2W11E2: case INV_W3W3: ff_cbrt(&t,w); if ( ! ff_sqrt(x,&t) ) return 0; else return 1;		// sixth
	case INV_W7E2: case INV_W2W7E2: case INV_W2W13:  case INV_W3W7:  if ( ! ff_sqrt(x,w) ) return 0; else  return 1;								// fourth
	case INV_W5E2: case INV_W3W5: case INV_W2W5E2: case INV_W3W3E2: case INV_W5W5: case INV_W3W11E2:  case INV_W2W17: ff_cbrt(x,w); return 1;	// cube
	case INV_W2W7: if ( ! ff_sqrt(&t,w) ) return 0; if ( ! ff_sqrt(x,&t) ) return 0; else return 1;														// square
	default: _ff_set(*x,*w); return 1;
	}
}


// for testing purposes
int ff_inv_count_for_j (ff_t j, int inv)
{
	ff_t f[PHI_ATKIN_MAX_FDEG+1];
	ff_t r[PHI_ATKIN_MAX_FDEG];

	inv_reduce_phi_fj (inv,0);
	bipoly_eval_ff (f, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, j);
	return ff_poly_distinct_roots(r,f,phi_fj_fdeg);
}


void ff_inv_single_eta_from_j (ff_t *w, ff_t *j, int inv)
{
	ff_t f[PHI_SINGLE_ETA_MAX_FDEG+1];
	ff_t r[PHI_SINGLE_ETA_MAX_FDEG];
	int n;
	
	inv_reduce_phi_fj (inv,0);
	bipoly_eval_ff (f, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, j[0]);
	n = ff_poly_distinct_roots(r,f,phi_fj_fdeg);
	if ( n < 1 || n > 2 ) { printf ("Error, got %d != 1,2 roots for single eta-quotient when converting from j=%ld (p=%ld)\n", n, _ff_get_ui(*j), _ff_p); exit(0); }				// note double root can happen
	if ( ! ff_inv_root (w, r, inv) ) { err_printf ("inv_root failed for j=%ld in ff_inv_single_eta_from_j with p=%ld\n", _ff_get_ui(*j), _ff_p); exit(0); }	
}

void ff_inv_double_eta_from_j (ff_t *w, ff_t *j, int inv)	// may return either w or its conjugate, you must pick one (consistently!)
{
	ff_t f[PHI_DOUBLE_ETA_MAX_FDEG+1];
	ff_t r[PHI_DOUBLE_ETA_MAX_FDEG];
	int i, n;
	
	inv_reduce_phi_fj (inv,0);
	bipoly_eval_ff (f, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, j[0]);
	n = ff_poly_distinct_roots(r,f,phi_fj_fdeg);
	if ( !n || n > 4 ) {
		if ( inv==INV_W3W5 || (inv > INV_MAX_SPECIFIC && ((inv-INV_DOUBLE_ETA)%2)) ) { printf ("Error, got %d != 1,2,4 roots for double eta-quotient when converting from j=%ld (p=%ld)\n", n, _ff_get_ui(*j), _ff_p); exit(0); }
		if ( n > 6 ) { printf ("Error, got %d != 1,2,3,6 roots for double eta-quotient with p1=2 when converting from j=%ld (p=%ld)\n", n, _ff_get_ui(*j), _ff_p); exit(0); }
		for ( i = 0 ; i < n ; i++ ) if ( ff_inv_root (w, r+i, inv) ) return;
		err_printf ("inv_root failed for all %d roots with j=%ld in ff_inv_double_eta_from_j with p=%ld\n", n, _ff_get_ui(*j), _ff_p); exit(0);
	}
	if ( ! ff_inv_root (w, r, inv) ) { err_printf ("inv_root failed  with j=%ld in ff_inv_double_eta_from_j with p=%ld\n", _ff_get_ui(*j), _ff_p); exit(0); }
}

void ff_inv_atkin_from_j (ff_t *a, ff_t *j, int inv)	// may return either w or its conjugate, you must pick one (consistently!)
{
	ff_t f[PHI_ATKIN_MAX_FDEG+1];
	ff_t r[PHI_ATKIN_MAX_FDEG];
	int n;

	inv_reduce_phi_fj (inv,0);
	bipoly_eval_ff (f, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, j[0]);
	n = ff_poly_distinct_roots(r,f,phi_fj_fdeg);
	if ( n < 1 || n > 2 ) { printf ("Error, got %d != 1,2 roots for atkin when converting from j=%ld (p=%ld)\n", n, _ff_get_ui(*j), _ff_p); exit(0); }				// note double root can happen
	_ff_set(*a,r[0]);
}

/*
void ff_inv_2atkin_from_j (ff_t a[2], ff_t *j, int inv)
{
	ff_t f[PHI_ATKIN_MAX_FDEG+1];
	ff_t r[PHI_ATKIN_MAX_FDEG];
	int n;

	inv_reduce_phi_fj (inv,0);
	bipoly_eval_ff (f, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, j[0]);
	n = ff_poly_distinct_roots(r,f,phi_fj_fdeg);
	if ( n < 1 || n > 2 ) { printf ("Error, got %d != 1,2 roots for atkin when converting from j=%ld over F_%ld\n", n, _ff_get_ui(*j), _ff_p); exit(0); }		// note double root can happen
	_ff_set(a[0],r[0]);  _ff_set(a[1],r[1]);
}

void ff_inv_2single_eta_from_j (ff_t a[2], ff_t *j, int inv)
{
	ff_t f[PHI_SINGLE_ETA_MAX_FDEG+1];
	ff_t r[PHI_SINGLE_ETA_MAX_FDEG];
	int n;

	inv_reduce_phi_fj (inv,0);
	bipoly_eval_ff (f, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, j[0]);
	n = ff_poly_distinct_roots(r,f,phi_fj_fdeg);
	if ( n < 1 || n > 2 ) { printf ("Error, got %d != 1,2 roots for single eta when converting from j=%ld over F_%ld\n", n, _ff_get_ui(*j), _ff_p); exit(0); }		// note double root can happen
	_ff_set(a[0],r[0]);  _ff_set(a[1],r[1]);
}

void ff_inv_4double_eta_from_j (ff_t a[4], ff_t *j, int inv)
{
	ff_t f[PHI_DOUBLE_ETA_MAX_FDEG+1];
	ff_t r[PHI_DOUBLE_ETA_MAX_FDEG];
	int n;

	inv_reduce_phi_fj (inv,0);
	bipoly_eval_ff (f, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, j[0]);
	n = ff_poly_distinct_roots(r,f,phi_fj_fdeg);
	if ( n < 1 || n > 4 ) { printf ("Error, got %d != 1,2,3,4 roots for single eta when converting from j=%ld over F_%ld\n", n, _ff_get_ui(*j), _ff_p); exit(0); }		// note double root can happen
	_ff_set(a[0],r[0]);  _ff_set(a[1],r[1]); _ff_set(a[2],r[2]); _ff_set(a[3],r[3]);
}
*/

int ff_inv_atkin_from_2j (ff_t a[2], ff_t *j1, ff_t *j2, int inv)
{
	static ff_t cache_f[PHI_ATKIN_MAX_FDEG+1], cache_j, cache_p;
	ff_t f[PHI_ATKIN_MAX_FDEG+1], g[PHI_ATKIN_MAX_FDEG+1];
	ff_t r[2];
	int n;
	
	if ( _ff_equal(*j1,*j2) ) { err_printf ("j1 == j2 == %ld in ff_inv_atkin_from_2j, p=%ld\n", _ff_get_ui(*j1), _ff_p); abort (); }
	inv_reduce_phi_fj (inv, 0);
//printf ("Converting j=%ld and j=%ld to atkin invariant\n", _ff_get_ui(*j1), _ff_get_ui(*j2));
	// check j1 against cached values, we expect to see the edges of a cycle
	if ( _ff_p == cache_p && _ff_equal(*j1,cache_j) ) {
		ff_poly_copy (f, &n, cache_f, phi_fj_fdeg);
		bipoly_eval_ff (g, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, *j2);
	} else {
		bipoly_eval_ff (f, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, *j1);
		bipoly_eval_ff (g, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, *j2);
	}
	ff_poly_copy (cache_f, &n, g, phi_fj_fdeg);  _ff_set(cache_j,*j2); cache_p = _ff_p;
//ff_poly_print(f,phi_fj_fdeg);
//ff_poly_print(g,phi_fj_fdeg);
	n = ( phi_fj_jdeg > 2 ? 2 : 1 );
	ff_poly_bgcd (r, &n, f, phi_fj_fdeg, g, phi_fj_fdeg);
//ff_poly_print(f,n);
	n = ff_poly_roots_d2(a,r,n);
//if ( n == 2 ) printf ("two atkin invariants %ld and %ld\n", _ff_get_ui(a[0]), _ff_get_ui(a[1])); else printf ("atkin invariant is %ld\n", _ff_get_ui(a[0]));
	return n;
}

int ff_inv_single_eta_from_2j (ff_t w[1], ff_t *j1, ff_t *j2, int inv)
{
	ff_t f[PHI_SINGLE_ETA_MAX_FDEG+1], g[PHI_SINGLE_ETA_MAX_FDEG+1];
	ff_t r[2];
	
	inv_reduce_phi_fj (inv, 0);
	inv_reduce_phi_fji (inv, 0);
	bipoly_eval_ff (f, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, *j1);
	bipoly_eval_ff (g, phi_fj_fdeg, phi_fji, phi_fj_terms, 1, *j2);
//ff_poly_print(f,phi_fj_fdeg);
//ff_poly_print(g,phi_fj_fdeg);
	ff_poly_clf (r, f, phi_fj_fdeg, g, phi_fj_fdeg);
	ff_poly_roots_d1 (f, r);
//printf ("single eta invariant is %ld\n", _ff_get_ui(w[0]));
	ff_inv_root (w, f, inv);
	return 1;	// currently failure is not permitted
}

int ff_inv_double_eta_from_2j (ff_t w[2], ff_t *j1, ff_t *j2, int inv)
{
	ff_t f[PHI_DOUBLE_ETA_MAX_FDEG+1], g[PHI_DOUBLE_ETA_MAX_FDEG+1];
	ff_t r[4];
	int n, d_g;

	if ( _ff_equal(*j1,*j2) ) { err_printf ("j1 == j2 == %ld in ff_inv_double_eta_from2j, p=%ld\n", _ff_get_ui(*j1), _ff_p); abort (); }
////printf ("Converting j=%ld and j=%ld to atkin invariant\n", _ff_get_ui(*j1), _ff_get_ui(*j2));
	inv_reduce_phi_fj (inv, 0);
//bipoly_print_mpz(Phi_atkin, phi_atkin_terms, 'X', 'Y');
//printf ("Reduced mod p=%ld\n", _ff_p); puts("");  bipoly_print_ff(phi_atkin, phi_atkin_terms, 'X','Y'); puts("");

//		ff_poly_derivative (g,&d_g,f,phi_fj_fdeg);
//		ff_poly_gcd_reduce (r,&d_r,f,phi_fj_fdeg,g,d_g);
//		if ( d_r != 1 ) { err_printf ("gcd degree %d not 1 in ff_inv_double_eta_from_2j\n", d_r);  ff_poly_print(r,d_r); printf ("j1=%ld, j2=%ld, p=%ld\n", _ff_get_ui(*j1), _ff_get_ui(*j2), _ff_p); exit(0); }
//	} else {
	bipoly_eval_ff (f, phi_fj_fdeg, phi_fj, phi_fj_terms, 1, *j1);
	bipoly_eval_ff (g, d_g=phi_fj_fdeg, phi_fj, phi_fj_terms, 1, *j2);
	n = ( phi_fj_jdeg > 2 ? 4 : 1 );
	ff_poly_bgcd (r, &n, f, phi_fj_fdeg, g, phi_fj_fdeg);
	if ( n > 2 ) { err_printf ("Got %d > 2 roots of gcd(Psi_f(X,j0),Psi_f(X,j1)) in ff_inv_single_eta_from_2j in ff_inv_double_eta_from_2j\n", n); return -1; }
	n = ff_poly_roots_d2(f,r,n);
	if ( n ) { ff_inv_root(w,f,inv); if ( n > 1 ) ff_inv_root(w+1,f+1,inv); }
	return n;
}

void ff_inv_j_from_t6 (ff_t *j, ff_t *t6)
{
	ff_t g2;
	register ff_t x, y, z;
	
	_ff_invert(y,*t6);
	_ff_square(x,*t6); 
	_ff_set_ui(z,6); ff_mult(z,z,*t6);
	_ff_subfrom(x,z); _ff_set_ui(z,27); _ff_subfrom(x,z);
	_ff_mult(g2,x,y);
	ff_inv_j_from_gamma2(j,&g2);
}

// the formula for lambda=u^8 is taken from Borwein+Borwein p. 115.
void ff_inv_j_from_u8 (ff_t *j, ff_t *u8)
{
	register ff_t t0,t1,t2,t3;

	_ff_square(t2,*u8);
	_ff_set_one(t0);  _ff_subfrom(t0,*u8);  _ff_add(t3,t0,t2); _ff_sub(t0,t3,*u8);  _ff_mult(t1,t2,t0);		// t1 = lambda^2(1-lamda)^2 = lambda^2(1-2lambda+lambda^2),  t3 = 1- lambda + lambda^2
	_ff_invert(t0,t1);
	_ff_square(t1,t3); _ff_mult(t2,t1,t3); _ff_mult(t1,t0,t2); _ff_set_ui(t0,256); _ff_mult(*j,t0,t1);			// j = 1728*J = 256(1-lambda+lambda^2)^3 / (lambda^2*(1-lambda)^2)
}

void ff_inv_j_from_atkin (ff_t *j, ff_t *a, int inv)
{
	ff_t f[PHI_ATKIN_MAX_JDEG+1];
	ff_t r[PHI_ATKIN_MAX_JDEG+1];
	int k, d_f;
	
	inv_reduce_phi_fj (inv,1);
	bipoly_eval_ff (f, phi_fj_jdeg, phi_fj, phi_fj_terms, 0, *a);
	if ( phi_fj_jdeg<=2 ) {
		k = ff_poly_roots_d2 (r, f, phi_fj_jdeg);		// note that roots_d2 will handle leading zero coeff
	} else {
		d_f = ff_poly_degree(f,phi_fj_jdeg);
		ff_poly_monic (f,&d_f,f,d_f);
		k = ff_poly_distinct_roots (r, f, d_f);
	}
	if ( ! k ) { err_printf ("No roots of Phi_fj(X,a) for a=%ld\n", _ff_get_ui(*a)); abort (); } 
	_ff_set(*j,r[0]);
}

void ff_inv_j_from_single_eta (ff_t *j, ff_t *w, int inv)
{
	ff_t f[PHI_SINGLE_ETA_MAX_JDEG+1];
	ff_t r[PHI_SINGLE_ETA_MAX_JDEG+1];
	
#if PHI_SINGLE_ETA_MAX_JDEG != 1
	err_printf ("PHI_SINGLE_ETA_MAX_JDEG=%d is not equal to 1 in ff_inv_j_from_single_eta\n", PHI_SINGLE_ETA_MAX_JDEG); exit(0);
#endif

	inv_reduce_phi_fj (inv,1);
	ff_inv_power(r,w,inv);
	bipoly_eval_ff (f, phi_fj_jdeg, phi_fj, phi_fj_terms, 0, r[0]);
	if ( _ff_zero(f[1]) ) { err_printf ("No roots of Phi_fj(w,Y) for w=%ld\n", _ff_get_ui(r[0])); abort (); } 
	ff_poly_roots_d1 (j,f);
}

void ff_inv_j_from_double_eta (ff_t *j, ff_t *w, int inv)
{
	ff_t f[PHI_DOUBLE_ETA_MAX_JDEG+1];
	ff_t r[PHI_DOUBLE_ETA_MAX_JDEG+1];
	int k, d_f;
	
	inv_reduce_phi_fj (inv,1);
	ff_inv_power(r,w,inv);
	bipoly_eval_ff (f, phi_fj_jdeg, phi_fj, phi_fj_terms, 0, r[0]);
	if ( phi_fj_jdeg<=2 ) {
		k = ff_poly_roots_d2 (r, f, phi_fj_jdeg);		// note that roots_d2 will handle leading zero coeff
	} else {
		d_f = ff_poly_degree(f,phi_fj_jdeg);
		ff_poly_monic (f,&d_f,f,d_f);
		k = ff_poly_distinct_roots (r, f, d_f);
	}
	if ( ! k ) { err_printf ("No roots of Phi_fj(X,a) for w=%ld\n", _ff_get_ui(*w)); abort (); } 
	_ff_set(*j,r[0]);
}


int ff_inv_j_from_2double_eta (ff_t *j, ff_t *w1, ff_t *w2, int inv)
{
	ff_t f[PHI_DOUBLE_ETA_MAX_JDEG+1], g[PHI_DOUBLE_ETA_MAX_JDEG+1], h[PHI_DOUBLE_ETA_MAX_JDEG+1];
	int d_f, d_g, d_h;
	
	inv_reduce_phi_fj (inv,1);
	ff_inv_power(h,w1,inv);  ff_inv_power(h+1,w2,inv);
	bipoly_eval_ff (f, phi_fj_jdeg, phi_fj, phi_fj_terms, 0, h[0]);
	bipoly_eval_ff (g, phi_fj_jdeg, phi_fj, phi_fj_terms, 0, h[1]);
	d_f = ff_poly_degree(f,phi_fj_jdeg);  d_g = ff_poly_degree(g,phi_fj_jdeg);
	ff_poly_monic (f,&d_f,f,d_f);
	ff_poly_gcd_small (h, &d_h, f, d_f, g, d_g);
	if ( d_h > 1 ) { err_printf ("GCD has degree %d > 1 in ff_inv_j_from_2double_eta\n", d_h); abort (); }
	if ( d_h < 1 ) return 0;
	if ( j ) ff_poly_roots_d1 (j, h);
	return 1;
}


void ff_inv_2j_from_atkin (ff_t j[2], ff_t *a, int inv)
{
	ff_t f[PHI_ATKIN_MAX_JDEG+1];
	ff_t r[PHI_ATKIN_MAX_JDEG+1];
	int k, d_f;
	
	inv_reduce_phi_fj (inv,1);
	bipoly_eval_ff (f, phi_fj_jdeg, phi_fj, phi_fj_terms, 0, *a);
	if ( phi_fj_jdeg<=2 ) {
		k = ff_poly_roots_d2 (r, f, phi_fj_jdeg);		// note that roots_d2 will handle leading zero coeff
	} else {
		d_f = ff_poly_degree(f,phi_fj_jdeg);
		ff_poly_monic (f,&d_f,f,d_f);
		k = ff_poly_distinct_roots (r, f, d_f);
	}
	if ( k!=2 ) { err_printf ("Got %d!=2 roots for Phi_fj(a,Y) for a=%ld\n", k, _ff_get_ui(*a)); abort (); } 
	_ff_set(j[0],r[0]); _ff_set(j[1],r[1]);
}

void ff_inv_2j_from_single_eta (ff_t j[2], ff_t *w, int inv)
{
	ff_t f[PHI_SINGLE_ETA_MAX_JDEG+1];
	ff_t t;

#if PHI_SINGLE_ETA_MAX_JDEG != 1
	err_printf ("PHI_SINGLE_ETA_MAX_JDEG=%d is not equal to 1 in ff_inv_2j_from_single_eta\n", PHI_SINGLE_ETA_MAX_JDEG); exit(0);
#endif

	inv_reduce_phi_fj (inv,1);
	ff_inv_power(&t,w,inv);
	bipoly_eval_ff (f, phi_fj_jdeg, phi_fj, phi_fj_terms, 0, t);
	if ( _ff_zero(f[1]) ) { err_printf ("No roots of Phi_fj(w,Y) for w=%ld\n", _ff_get_ui(t)); abort (); } 
	ff_poly_roots_d1 (j,f);
	inv_reduce_phi_fji (inv,1);
	bipoly_eval_ff (f, phi_fj_jdeg, phi_fji, phi_fj_terms, 0, t);
	if ( _ff_zero(f[1]) ) { err_printf ("No roots of Phi_fji(w,Y) for w=%ld\n", _ff_get_ui(t)); abort (); } 
	ff_poly_roots_d1 (j+1, f);
}

void ff_inv_2j_from_double_eta (ff_t j[2], ff_t *w, int inv)
{
	ff_t f[PHI_DOUBLE_ETA_MAX_JDEG+1];
	ff_t r[PHI_DOUBLE_ETA_MAX_JDEG+1];
	int k, d_f;
	
	inv_reduce_phi_fj (inv,1);
	ff_inv_power(r,w,inv);
	bipoly_eval_ff (f, phi_fj_jdeg, phi_fj, phi_fj_terms, 0, r[0]);
	if ( phi_fj_jdeg<=2 ) {
		k = ff_poly_roots_d2 (r, f, phi_fj_jdeg);		// note that roots_d2 will handle leading zero coeff
	} else {
		d_f = ff_poly_degree(f,phi_fj_jdeg);
		ff_poly_monic (f,&d_f,f,d_f);
		k = ff_poly_distinct_roots (r, f, d_f);
	}
	if ( k!=2 ) { err_printf ("Got %d!=2 roots for Phi_fj(w,Y) for w=%ld\n", k, _ff_get_ui(*w)); abort (); } 
	_ff_set(j[0],r[0]); _ff_set(j[1],r[1]);
}
