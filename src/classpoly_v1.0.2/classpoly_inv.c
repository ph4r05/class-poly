#include <stdlib.h>
#include <stdio.h>
#include "ff_poly.h"
#include "phi_poly.h"
#include "class_inv.h"
#include "classpoly.h"
#include "classpoly_crt.h"
#include "classpoly_inv.h"
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


int classpoly_inv_setup (classpoly_inv_t Hinv, classpoly_t H, int inv, classpoly_crt_t crt)
{
	clock_t start, end;
	ff_t t0, t1, j;
	unsigned long z;
	mpz_t M, X, Y, Z, HT;
	register int i, stabcnt;

	memset (Hinv, 0, sizeof(*Hinv));
	Hinv->inv = inv;
	inv_degree (&Hinv->p1, &Hinv->p2, inv);
	Hinv->inverted = inv_inverted_involution (inv);  Hinv->negated = inv_negated_involution (inv);
	if ( inv_double_eta(inv) && kronecker_p(H->D,Hinv->p1) && kronecker_p(H->D,Hinv->p2) ) {  	// if neither p1 or p2 are ramified in w_{p1,p2} we need to distinguish two distinct class polys
		Hinv->ambiguous_H = 1;
		if ( ! qform_n_rep (Hinv->Nrep, Hinv->p1*Hinv->p2, H->D) ) { err_printf ("qform_n_rep failed with inv_p1=%d, inv_p2=%d, D=%ld\n", Hinv->p1, Hinv->p2, H->D); return 0; }
	}
	Hinv->units = inv_units (inv);
	Hinv->ambiguous_sign = inv_sqrt_invariant(inv);
	mpz_init (Hinv->T);
	if ( Hinv->ambiguous_sign && ! ((H->pres->enum_cnt&1) && Hinv->units) ) {
		start = clock();
		mpz_init (M);  mpz_init (X);  mpz_init (Y);  mpz_init (Z);  mpz_init (HT);
		mpz_set_ui(M,1);
		stabcnt = 0;
		while ( stabcnt < 3 && classpoly_crt_next_prime (crt) ) {
			if ( ! ff_classpoly_find_jinv (&j, H) ) { err_printf ("Fatal error, ff_classpoly_find_jinv failed for p=%ld\n", _ff_p); abort(); }
			classpoly_crt_cache_j (crt, j);
			ff_classpoly_setup (H);
			ff_classpoly_inv_from_j (&j, &j, Hinv, H);
			if ( ! ff_classpoly_enum_roots (H, j, inv) ) { err_printf ("Call to ff_classpoly_enum_roots failed\n"); abort(); }
			ff_classpoly_inv_disambiguate (H, Hinv);
			// We assume that the candidate polys are for invariants f and -f, so the sum of the traces is always zero and the product is minus the square of the trace
			_ff_set_zero(t0); for ( i = 0 ; i < H->pres->enum_cnt ; i++ ) _ff_addto(t0, H->roots[i]); 
			_ff_square(t1,t0); z = _ff_get_ui(t1); 
			if ( crt_mpz_ui(X, Y, Z, M, z, _ff_p) ) stabcnt++; else stabcnt = 0;
			mpz_set(Z, X);  mpz_set(M,Y);
		}
		classpoly_crt_restart(crt);
		if ( stabcnt < 3 ) { err_printf ("square of trace did not stabilize\n"); abort(); }
		crt_mpz_mod_to_Z(Z,M);
		end = clock();
		info_printf ("Used %ld msecs and %d primes to compute square of trace\n", delta_msecs(start,end), crt->jcnt);
		if ( ! mpz_sgn(Z) ) { char filename[256]; FILE *fp; err_printf ("Unable to use zero trace to disambiguate class polys\n");
sprintf (filename, "H_%ld.txt", -H->D);  fp = fopen (filename, "w");  fclose (fp);			// temporary hack for testing with hilbertbatch
			return 0; }
		if ( mpz_sgn(Z) < 0 || ! mpz_perfect_square_p(Z) ) { err_printf ("Error, square of trace is not a perfect square!\n"); abort(); }
		mpz_sqrt(Hinv->T,Z);
		if ( dbg_level >= 2 ) gmp_printf("Classpoly trace is %Zd\n\n", Hinv->T);
	}
	return 1;
}

void classpoly_inv_clear (classpoly_inv_t Hinv)
	{ mpz_clear (Hinv->T); }

void _ff_classpoly_inv_from_j (ff_t *F, ff_t *J, classpoly_inv_t Hinv, classpoly_t H)
{
	struct phi_vshape_struct s;
	ff_t t1;
	register int i,j,p1,p2;

	p1 = Hinv->p1;  p2 = Hinv->p2;
	/*
		In situations where there is more than one f-invariant related to a given j-invariant j1, we look for p1*p2-isogenous j2 such that there is a unique f-invariant related to j1 and j2
		i.e., we may compute f as gcd(Psi^f(X,j1),Psi^f(X,j2)) where Psi^f(X,J) is the polynomial relation f(z) and j(z), see Lemma 2 in Enge-Sutherland "Class invariants by the CRT method"
	*/
	phi_vshape (&s, H->v);
	for ( i = 0 ; i < s.k ; i++ ) if ( p1 == s.p[i] ) break;
	j = ( i < s.k ? s.h[i] : 0 );
	if ( ! phi_next_surface_neighbor (&t1, 0, p1, j, *J, 0, 0) ) { err_printf ("Unable to find %d-neighbor of J=%ld, volcano depth=%d\n", p1, _ff_get_ui(*J), j); abort(); }
	dbg_printf ("%d-neighbor of j=%ld is %ld\n", p1, _ff_get_ui(*J), _ff_get_ui(t1));
	if ( p2 > 1 ) {
		if ( p1 == p2 ) {
			if ( ! phi_next_surface_neighbor (&t1, 0, p1, j, t1, J, 0) ) { err_printf ("Unable to find %d-neighbor of J=%ld, volcano depth=%ld\n", p1, _ff_get_ui(t1), t1); abort(); }
			dbg_printf ("%d-neighbor of j=%ld is %ld\n", p1, _ff_get_ui(*J), _ff_get_ui(t1));
		} else {
			for ( i = 0 ; i < s.k ; i++ ) if ( p2 == s.p[i] ) break;
			j = ( i < s.k ? s.h[i] : 0 );
			if ( ! phi_next_surface_neighbor (&t1, 0, p2, j, t1, 0, 0) ) { err_printf ("Unable to find %d-neighbor of J=%ld, volcano depth=%ld\n", p1, _ff_get_ui(t1), t1); abort(); }
			dbg_printf ("%d-neighbor is j=%ld\n", p2, _ff_get_ui(t1));
			if ( _ff_equal (t1,*J) ) { err_printf ("Discriminant %ld not supported for invariant w_{%d,%d}, the primes over p1 and p2 must lie in distinct classes of cl(D)\n", H->D, p1, p2); abort(); }
		}
	}
	ff_inv_from_2j (F, J,&t1,Hinv->inv);
}

int _ff_classpoly_inv_disambiguate (classpoly_t H, classpoly_inv_t Hinv)
{
	torsor_t T;
	long e[IQ_MAX_GENS];
	ff_t t1,t2;
	int negate;
	register int i, j;
	
	if ( Hinv->ambiguous_H ) {
		if ( ! torsor_setup (T, H, 0) ) { err_printf ("call to torsor_setup failed\n"); abort(); }
		torsor_orient_evec (e, Hinv->Nrep, T);  j = torsor_action (0, e, T);
		if ( ! ff_j_from_2inv (&t1,H->roots,H->roots+j,Hinv->inv) ) {
			if ( Hinv->inverted ) ff_parallel_invert (H->roots, H->roots, H->pres->enum_cnt);
			if ( Hinv->negated ) for ( i = 0 ; i < H->pres->enum_cnt ; i++ ) ff_negate(H->roots[i]);
		}
	}
	if ( Hinv->ambiguous_sign ) {
		negate = 0;
		if ( (H->pres->enum_cnt&1) && Hinv->units ) {
			_ff_set_one(t1);  ff_negate(t1);
			for ( i = 0 ; i < H->pres->enum_cnt ; i++ ) ff_mult(t1,t1,H->roots[i]);
			if ( ! _ff_one(t1) ) {
				ff_negate(t1); if ( ! _ff_one(t1) ) { err_printf ("Classpoly with ambiguous sign doesn't have constant term +/-1, this should never happen!\n"); return 0; }
				negate = 1;
			}
		} else {
			_ff_set_mpz(t1,Hinv->T);  if ( _ff_zero(t1) ) return 0;
			_ff_set(t2,0);  for ( i = 0 ; i < H->pres->enum_cnt; i++ ) _ff_addto(t2,H->roots[i]);
			if ( ! _ff_equal (t1, t2) ) {
				ff_negate(t2); if ( !_ff_equal(t1,t2) ) { err_printf ("Trace %ld mod %ld is incorrect, it should be +/-%ld, this should never happen!\n", _ff_get_ui(t2), _ff_p, _ff_get_ui(t1)); return 0; }
				negate = 1;
			}			
		}
		if ( negate ) for ( i = 0 ; i < H->pres->enum_cnt ; i++ ) ff_negate(H->roots[i]);
	}
	return 1;
}
