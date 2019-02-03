#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "ff_poly.h"
#include "classpoly.h"
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

#define BIG_D	1000000

int main (int argc, char *argv[])
{
	long D;
	int inv;
	mpz_t P;
	char *filename;
	int picked;

	if ( argc < 2 ) { fprintf (stderr, "classpoly D [inv P filename verbosity]\n"); fprintf (stderr, "Version %s (ff_poly version %s)\n", CLASSPOLY_VERSION_STRING, FF_POLY_VERSION_STRING); return 0; }

	_cstd_seed = 42;
	
	mpz_util_init();
	
	D = atol(argv[1]);  if ( D > 0 ) D = -D;
	if ( ! discriminant_test(D) ) { err_printf ("%ld is not a valid discriminant.\n", D); return -1; }
	if ( argc > 5 ) { dbg_level = atoi (argv[5]); } else { dbg_level = ( -D > BIG_D ? 1 : 0 ); }
	picked = 0;
	if ( argc > 2 ) {
		inv = atoi(argv[2]);
		if ( inv == - 1 ) { inv = inv_pick_invariant(D); picked = 1; info_printf ("Using invariant %s (%d)\n", inv_string(inv), inv); }
		if ( ! inv_good_invariant(inv) ) { err_printf ("Unknown or invalid invariant %d specified\n", inv); return -1; }
		if ( ! inv_enum(inv) ) { err_printf ("Unsupported invariant %d specified\n", inv); return -1; }
		if ( ! inv_good_discriminant (D, inv) ) { err_printf ("D=%ld is not a good discriminant for invariant %s (%d)\n", D, inv_string(inv), inv); return -1; }
	} else {
		if ( D >= -4 ) inv = 0; else { inv = inv_pick_invariant(D); picked = 1; }
		info_printf ("Using invariant %s (%d)\n", inv_string(inv), inv);
	}
	
	mpz_init (P);
	if ( argc > 3 )  mpz_eval_expr (P, argv[3]);
	filename = ( argc > 4 ? argv[4] : 0 );
	if ( ! compute_classpoly (D, inv, P, filename) ) {
		if ( ! picked || ! inv ) return 0;
		err_printf ("Unable to compute classpoly using invariant %s (%d), reverting to j-invariant\n", inv_string(inv), inv);
		if ( ! compute_classpoly (D, 0, P, filename) ) { err_printf ("Error computing Hilbert class polynomial!\n"); return-1; }
	}
	mpz_clear (P);
	return 0;
}
