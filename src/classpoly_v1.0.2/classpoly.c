#include <stdlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <gmp.h>
#include <errno.h>
#include "mpzutil.h"
#include "evec.h"
#include "qform.h"
#include "ff_poly.h"
#include "class_inv.h"
#include "ecurve.h"
#include "phi_poly.h"
#include "classpoly.h"
#include "findcurve.h"
#include "polycosts.h"
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

int ecurve_test_order2 (ppf_t n[2], ff_t f[4]);	// not defined in ecurve.h

char _H_dir_str[1024];

#define BUFSIZE	65536

static char buf[BUFSIZE];

/*
	Loads a class polynomial for discriminant D from the ascii file H_d.txt.  Reads invariant from file.
	This function is obsolete.
*/
int classpoly_load_mpz (mpz_t H[], int n, int *inv, long D)
{
	FILE *fp;
	char *s;
	char filename[256];
	long d, i;
	
	if ( D > 0 ) D = -D;
	d = -D;
	sprintf(filename,"%s/H_%ld.txt", H_DIR, d);
	info_printf ("Loading %s\n", filename);
	fp = fopen(filename,"r");
	if ( ! fp ) { if ( errno == ENOENT ) { info_printf ("File %s not found\n", filename); return 0; } else { err_printf ("Unexpected error while attempting to open file %s\n", filename); perror("fopen"); abort(); } }
	if ( ! fgets(buf,sizeof(buf),fp) ) { fclose(fp); return 0; }
	if ( buf[0] != 'I' ) { err_printf ("Expected I in first line of %s\n", filename); fclose (fp); return 0; }
	*inv = atoi(buf+2);
	if ( ! fgets(buf,sizeof(buf),fp) ) { fclose(fp); return 0; }
	if ( buf[0] != 'D' ) { err_printf ("Expected D in second line of %s\n", filename); fclose (fp); return 0; }
	i = atoi(buf+2);
	if ( i != D ) { err_printf ("Discriminant mismatch in file %s\n", filename); fclose (fp); return 0; }
	for ( i = 0 ; fgets(buf,sizeof(buf),fp) ; i++ ) {
		if ( i==0 && buf[0] == 'P' ) { err_printf ("Class poly in file %s is reduced mod P, expected unreduced class poly\n", filename); fclose (fp); return 0; }
		if ( i>=n ) { err_printf ("Specified coefficient array of size %d too small for class polynomial with discriminant %ld\n", n, D); fclose (fp); return 0; }
		for ( s = buf ; *s && *s != '*' ; s++ );
		if ( ! *s ) { err_printf ("Unable to process record %ld in file %s, increase buffer size?\n", i, filename); fclose (fp); return 0; }
		*s = '\0';
		mpz_set_str(H[i],buf,10);
	}
	fclose (fp);
	return i-1;
}


int classpoly_init (classpoly_t H, long D, int inv)
{
	unsigned long p[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS];
	int i, k;

	k = discriminant_factor_conductor (p, h, D);
	if ( k < 0 ) { err_printf ("D=-%ld is not a valid discriminant\n", D);  return 0; }
	if ( k > MAX_P32 ) { err_printf ("Conductor of discriminant %ld has too many prime factors\n", D); abort(); }
	if ( ! inv_good_discriminant (D, inv) ) { err_printf ("D=%ld is not a good discriminant for specified invariant %d\n", D, inv); return 0; }
	if ( k ) { dbg_printf ("D=%ld is not a fundamental discriminant, conductor is ",D);  for ( i = 0 ; i < k ; i++ ) dbg_printf ("%ld^%ld ", p[i], h[i]);  dbg_printf("\n"); }
	memset (H, 0, sizeof(*H));
	H->D = D;
	H->inv = inv;
	H->u = mem_alloc (sizeof(struct ifactors_struct));
	H->u->k = k;
	for ( i = 0 ; i < k ; i++ ) {H->u->p[i] = (int)p[i];  H->u->h[i] = (int) h[i]; }
	return 1;
}


void classpoly_setup_find_jinv (classpoly_t H)
{
	long D1;
	int i, j;

	if ( H->H1 ) return;
	if ( ! H->u ) { err_printf ("You must call classpoly_init before calling classpoly_setup_find_jinv\n"); abort(); }
	/*
		Here we modify step 2 of Algorithm 1.2 to handle large primes ell dividing the conductor by using H_{D/ell^2}
		A curve E with End(E) lying in the correct maximal order O will have the correct power of ell dividing the conductor of O iff E is on the floor of its ell-volcano.
		Since ell > PHI_MAX_M, which is the largest prime that can divide v(p) (where 4p=tr(E)^2-v(p)^2D)), E is on the floor iff H_{D/(ell*ell)}(j(E)) != 0.
	
		So for each prime ell > PHI_MAX_M dividing the conductor of D we create a classpoly for H_{D/ell^2} and store this in a list of classpoly's H1_list.
		Note that we assume PHI_MAX_M > log|D|, so this does not increase the computational complexity (and in practice it takes negligible time,
		we compute fewer than log|D| class polys, all of which have discriminant of absolute value smaller than |D| by at least a log^2|D| factor).
	*/
	for ( i = 0 ; i < H->u->k ; i++ ) if ( H->u->p[i] > PHI_MAX_M ) break;
	if ( i < H->u->k ) {
		H->H1_len = H->u->k-i;
		H->H1 = malloc (H->H1_len*sizeof(*H->H1));
		for ( j = 0 ; j < H->H1_len ; j++ ) {
			D1 = H->D/(H->u->p[i+j]*H->u->p[i+j]);
			if ( ! classpoly_create (H->H1+j, D1, 0) ) { err_printf ("error creating classpoly for D1=%ld\n", D1); abort(); }
			info_printf ("Created classpoly for D1=%ld to handle prime %d dividing the conductor of %ld\n", D1, H->u->p[i+j], H->D);
		}
	}
	return;
}


int classpoly_setup_enum_roots (classpoly_t H, long ellfilter, long ell0, int enum_inv)
{
	if ( ! ellfilter ) ellfilter = 1;
	if ( H->pres && H->pres->enum_inv == enum_inv && !(H->pres->ellfilter % ellfilter)  && H->pres->ell0 == ell0) return 1;
	if ( ! H->pres ) H->pres = mem_alloc (sizeof (*H->pres));
	dbg_printf ("Computing presentation for cl(%ld) using initial ellfilter=%ld, ell0=%ld, enum_inv=%d\n", H->D, ellfilter, ell0, enum_inv);
	if ( classgroup_pcp_setup (H->pres, H->D, ellfilter, ell0, enum_inv) < 0 ) return 0;
	// qform_table_free ();
	dbg_printf ("Computed presentation for cl(%ld) with %d generators\n", H->D, H->pres->k);
	// the sanity check below is unnecessary but doesn't hurt
	if ( H->C ) {
		if ( ell0 ) {
			if ( H->F_d*H->G_d != H->pres->h/2 ) { err_printf ("F_d=%d * G_d=%d is not equal to half the class number h=%ld for D=%ld\n", H->F_d, H->G_d, H->pres->h, H->D); abort(); }
		} else {
			if ( H->F_d*H->G_d != H->pres->h ) { err_printf ("F_d=%d * G_d=%d is not equal to the class number h=%ld for D=%ld\n", H->F_d, H->G_d, H->pres->h, H->D); abort(); }
		}
	}
	// reallocate in case enum_cnt has changed
	if ( H->roots ) mem_free (H->roots);
	H->roots = mem_alloc ((H->pres->enum_cnt+PHI_MAX_H+1)*sizeof(ff_t));

	return 1;
}


// find suborder datastructure, create if needed
struct suborder_struct *classpoly_get_suborder (classpoly_t H, long index)
{
	unsigned long p[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS];
	struct suborder_struct *o;
	int i;
	
	for ( i = 0 ; H->suborders[i] ; i++ ) if ( H->suborders[i]->index == index ) break;
	if ( H->suborders[i] ) return H->suborders[i];
	if ( i >= MAX_SUBORDERS ) { err_printf ("Exceeded MAX_SUBORDERS=%d\n", MAX_SUBORDERS);  abort(); }
	if ( LONG_MAX/(index*index) <= -H->D ) { err_printf ("index %ld is too large for discriminant %ld\n", index, H->D); abort(); }
	H->suborders[i] = o = mem_alloc (sizeof(*H->suborders[i]));
	o->D = index*index*H->D;
	o->index = index;
	o->u = malloc (sizeof(*o->u));
	o->u->k = ui_factor (p, h, index);
	for ( i = 0 ; i < o->u->k ; i++ ) {o->u->p[i] = (int)p[i];  o->u->h[i] = (int)h[i]; }
	dbg_printf ("Created suborder of index %ld relative to D=%ld\n", index, H->D);
	return o;
}


int classpoly_suborder_setup_enum_roots (classpoly_t H, long ellfilter, long ell0, long index, int enum_inv)
{
	struct suborder_struct *o;

	if ( ! ellfilter ) ellfilter = 1;
	o = classpoly_get_suborder (H, index);
	if ( ! o ) { err_printf ("Error creating suborder with index %ld\n", index); abort(); }
	if ( o->pres && o->pres->enum_inv == enum_inv && !(o->pres->ellfilter % ellfilter)  && o->pres->ell0 == ell0) return 1;
	if ( ! o->pres ) o->pres = mem_alloc (sizeof (*o->pres));
	if ( ! inv_enum (enum_inv) ) { err_printf ("Invariant %d (%s) is not suitable for enumeration\n", enum_inv, inv_string(enum_inv)); abort(); }
	dbg_printf ("Computing presentation for cl(%ld^2*%ld)\n", index, H->D);
	if ( classgroup_pcp_setup (o->pres, o->index*o->index*H->D, ellfilter, ell0, enum_inv) < 0 ) return 0;
	dbg_printf ("Computed presentation for cl(%ld) with %d generators\n", o->D, o->pres->k);
	o->roots = mem_alloc ((o->pres->enum_cnt+PHI_MAX_H+1)*sizeof(ff_t));
	return 1;
}


int classpoly_create (classpoly_t H, long D, int inv)
{	
	if ( ! classpoly_load (H, D, inv) ) {
		if ( ! compute_classpoly (D, inv, 0, 0) ) { printf ("call to compute_classpoly failed with D=%ld, inv=%d\n", D, inv); return 0; }
		if ( ! classpoly_load (H, D, inv) ) {err_printf ("classpoly_load failed unexpectedly with D=%ld, inv=%d\n", D, inv); return 0; }
	}
	return 1;
}

int classpoly_load (classpoly_t H, long D, int inv)
{
	FILE *fp;
	char *s;
	char filename[256];
	long loadD;
	int i, H_inv, F_d, G_d;

	if ( D > 0 ) D = -D;
	sprintf(filename,"H_%ld.txt", -D);
	fp = fopen(filename,"r");
	if ( ! fp ) {
		sprintf(filename,"%s/H_%ld.txt", H_DIR, -D);
		fp = fopen(filename,"r");
		if ( ! fp ) { if ( errno == ENOENT ) { info_printf ("File %s not found\n", filename); return 0; } else { err_printf ("Unexpected error while attempting to open file %s\n", filename); perror("fopen"); abort(); } }
	}
	info_printf ("Loading %s\n", filename);
	H_inv = 0; F_d = G_d = 0;  loadD = 0;
	while ( fgets(buf,BUFSIZE,fp) ) {
		switch (buf[0]) {
		case 'I': H_inv = atoi (buf+2); break;
		case 'D': loadD = atol (buf+2); break; 
		case 'P': if ( buf[2] != '0' ) { err_printf ("non-zero modulus P in classpoly file %s, expected poly over Z\n", filename); fclose(fp); return 0; }
		case 'F': F_d = atoi(buf+2); break;
		case 'G': G_d = atoi(buf+2); break;
		}
		if ( ! isalpha(buf[0]) ) break;
	}
	if ( ! F_d && ! G_d ) {
		// support legacy file format in which the degree was not specified (scan the file for highest degree term)
		G_d = 1;  F_d = 0;
		do {
			for ( s = buf ; *s && *s != '^' ; s++ );
			if ( ! *s ) break;
			i = atoi(s+1);
			if ( i > F_d ) F_d = i;
		} while ( fgets(buf, BUFSIZE,fp) );
		rewind (fp);
		while ( fgets(buf,BUFSIZE,fp) && isalpha(buf[0]) );
	}
	if ( ! F_d || ! G_d ) { err_printf ("One of F= or G= is missing in classpoly file %s\n", filename);  fclose(fp); return 0; }
	if ( F_d > FF_POLY_MAX_DEGREE || G_d > FF_POLY_MAX_DEGREE ) { err_printf ("F_d=%d or G_d=%d exceeds FF_POLY_MAX_DEGREE=%d\n", F_d, G_d, FF_POLY_MAX_DEGREE); fclose(fp); return 0; }
	if ( loadD != D && loadD != -D ) { err_printf ("Discriminant mismatch in file %s\n", filename);  fclose(fp); return 0; }
	if ( inv >= 0 && inv != H_inv ) { fclose(fp); return 0; }
	if ( ! classpoly_init (H, D, inv) ) { fclose(fp); return 0; }
	H->D = D; H->inv = H_inv;  H->F_d = F_d;  H->G_d = G_d;
	H->C = mem_alloc((F_d*G_d+1)*sizeof(mpz_t));   H->c = mem_alloc ((F_d*G_d+1)*sizeof(ff_t));
	for ( i = 0 ; i < F_d*G_d+1 ; i++ ) {
		if ( i && ! fgets(buf,BUFSIZE,fp) ) { puts ("barf"); err_printf ("Unable to read %d of %d coefficients from classpoly file %s\n", i+1, F_d*G_d+1, filename); fclose(fp); return 0; }
		for ( s = buf ; *s && *s != '*' ; s++ );
		if ( ! *s ) { err_printf ("Unable to process coefficient %d in file %s, increase buffer size?\n", i, filename); fclose(fp); return 0; }
		*s = '\0';
		mpz_init_set_str(H->C[i],buf,10);
		if ( i == F_d && mpz_cmp_ui(H->C[i],1) != 0 ) { err_printf ("Coefficient %d inf classpoly file %s is not 1, expected F to be monic!\n", F_d, filename); fclose(fp); return 0; }
	}
	fclose (fp);
	H->pres = 0;
	H->p = 0;
	info_printf ("Loaded class polynomial for D=%ld inv=%d from file %s\n", D, inv, filename);
	return 1;
}

void classpoly_clear (classpoly_t H)
{
	int i;
	
	if ( H->C ) for ( i = 0 ; i < H->F_d*H->G_d+1 ; i++ ) mpz_clear (H->C[i]); mem_free (H->C);
	if ( H->c ) mem_free (H->c);
	if ( H->pres ) mem_free (H->pres);
	if ( H->roots ) mem_free (H->roots);  
	if ( H->H1 ) for ( i = 0 ; i < H->H1_len ; i++ ) classpoly_clear (H->H1+i); 
	if ( H->u ) mem_free (H->u);
	for ( i = 0 ; H->suborders[i] ; i++ ) {
		if ( H->suborders[i]->u ) mem_free (H->suborders[i]->u);
		if ( H->suborders[i]->pres ) mem_free (H->suborders[i]->pres);
		if ( H->suborders[i]->roots ) mem_free (H->suborders[i]->roots);
	}
	memset (H, 0, sizeof(*H));
}


void ff_classpoly_reduce (classpoly_t H)
{
	register int i;

	if ( ! _ff_p ) { err_printf ("Finite field not set in ff_classpoly_reduce!\n"); abort(); }
	if ( ! H->C ) { err_printf ("Call to ff_classpoly_reduce with no class poly loaded!\n"); abort(); }
	if ( _ff_p != H->redp ) { for ( i = 0 ; i < H->F_d*H->G_d+1 ; i++ ) _ff_set_mpz (H->c[i], H->C[i]);  H->redp = _ff_p; }
}

int ff_classpoly_root (ff_t r[1], classpoly_t H)
{
	ff_t G[FF_POLY_MAX_DEGREE+1];
	register ff_t d, x, y, t;
	ff_t t0;
	register int i, j, k;
	
	if ( ! H->C )  { err_printf ("Call to ff_classpoly_root with no classpoly loaded!\n"); abort(); }
	// reduce mod p
	if ( ! ff_classpoly_setup (H) ) { err_printf ("p=%ld does not split completely in the ring class field for D=%ld\n", _ff_p, H->D); abort(); }
	ff_classpoly_reduce (H);
	if ( dbg_level >= DEBUG_LEVEL ) { printf ("F mod %ld = ", _ff_p);  ff_poly_print(H->c, H->F_d); }
	// Get a root r of F
	if ( ! ff_poly_find_root(r,H->c,H->F_d) ) { err_printf ("Couldn't find a root of F mod p=%ld, in ff_classpoly_root, F =", _ff_p); ff_poly_print(H->c, H->F_d); abort(); }
	dbg_printf("r = %ld is a root of F mod p\n", _ff_get_ui(r[0]));
	if ( H->G_d == 1 ) return 1;

	// Compute F'(r)
	_ff_set_ui(d,H->F_d);  _ff_set(x,r[0]);  _ff_set(y,d);
	for ( i = H->F_d-1 ; i ; i-- ) { ff_mult(y,y,x);  _ff_dec(d); _ff_mult(t,d,H->c[i]); _ff_addto(y,t); }
	if ( _ff_zero(y) ) {		// this will happen if F has a double root mod p, which is possible for certain p 
		out_printf ("F and F' have r=%ld as a common root mod %ld, removing r and trying again\n", _ff_get_ui(r[0]), _ff_p);
		ff_poly_remove_root (G, H->c, H->F_d, r);
		ff_poly_remove_root (G, G, H->F_d-1, r);
		for ( k = H->F_d-2 ; k ; k-- ) {
			if ( ! ff_poly_find_root(r,G,k) ) return 0;	// should never happen
			_ff_set_ui(d,H->F_d);  _ff_set(x,r[0]);  _ff_set(y,d);
			for ( i = H->F_d-1 ; i ; i-- ) { ff_mult(y,y,x);  _ff_dec(d); _ff_mult(t,d,H->c[i]); _ff_addto(y,t); }
			if ( ! _ff_zero(y) ) break;
			ff_poly_remove_root (G, G, k, r);
		}
		if ( k <= 0 ) { err_printf ("Every root of F is also a root of F' mod %ld, giving up!\n", _ff_p); ff_poly_print(H->c,H->F_d); exit (0); }
	}
	dbg_printf("F'(r)=%ld\n", _ff_get_ui(y));
	_ff_invert(d,y);
	dbg_printf("1/F'(r)=%ld\n", _ff_get_ui(d));

	// Compute the coefficients G_i(r)/F'(r) of G
	for ( i = H->F_d+1, j = 0 ; j < H->G_d-1 ; j++, i+= H->F_d )
		{ ff_poly_eval(&t0,H->c+i,H->F_d-1,r); ff_mult(G[j],d,t0); }
	_ff_set(G[H->G_d-1],r[0]);  _ff_set_one(G[H->G_d]);
	 if ( dbg_level >= DEBUG_LEVEL ) { printf ("G = ");ff_poly_print(G,H->G_d); }

	// return a root of G
	if ( ! ff_poly_find_root(r,G,H->G_d) ) { err_printf ("Couldn't find a root of G mod p=%ld, in ff_classpoly_root, G =", _ff_p); ff_poly_print(G, H->G_d); abort(); }
	dbg_printf ("ff_classpoly_root found root %ld of classpoly for D=%ld, inv=%d, mod p=%ld\n", _ff_get_ui(r[0]), H->D, H->inv, _ff_p);
	return 1;
}

int ff_classpoly_jinv (ff_t j[1], classpoly_t H)
{
	ppf_t N[2];
	ff_t r[1], f[4], jlist[PHI_MAX_JDEG];
	register int i,k;

	if ( ! ff_classpoly_root(r,H) ) return 0;
	if ( H->inv == INV_J ) { _ff_set(j[0],r[0]); return 1; }
	dbg_printf("Obtaining j-invariant from classpoly for D=%ld, inv=%d\n", H->D, H->inv);
	k = ff_multi_j_from_inv (jlist, r, H->inv, PHI_MAX_JDEG, 1);
	dbg_printf("%d candidate j-invariants\n", k);
	if ( k > 2 ) {
		dbg_printf ("Found %d j's for root %ld of class poly for inv=%d (%s), testing orders\n", k, _ff_get_ui(r[0]), H->inv, inv_string(H->inv));
		ppf_factor(N[0],_ff_p+1-H->t);  ppf_factor(N[1],_ff_p+1+H->t);
		for ( i = 0 ; i < k ; i++ ) {
			ecurve_from_jinv(f,jlist+i);
			if ( ecurve_test_order2 (N, f) ) break;
		}
		if ( i == k ) { err_printf ("None of the candidate j's for inv=%d(%s) class poly root %ld yield a curve with the correct order over F_%ld!\n", H->inv, inv_string(H->inv), _ff_get_ui(r[0]),_ff_p); abort(); }
		// TODO: verify that the endomorphism ring is correct!
		dbg_printf ("Succeeded after %d tries\n", i+1);
		_ff_set(j[0],jlist[i]);
	} else {
		_ff_set(j[0],jlist[0]);
	}
	dbg_printf ("Obtained j-invariant %ld from root %ld of H_%ld mod %ld (inv=%d)\n", _ff_get_ui(j[0]), _ff_get_ui(r[0]), -H->D, _ff_p, H->inv);
	return 1;
}


int ff_classpoly_suborder_root (ff_t r[1], classpoly_t H, long index)
{
	struct suborder_struct *o;
	struct ifactors_struct *u;
	long q, v;
	int i, k, d, e;
	
	if ( index < 0 ) { err_printf ("Invalid index=%ld\n", index); abort(); }
	if ( ! (o = classpoly_get_suborder (H, index)) ) return 0;
	if ( ! ff_classpoly_root (r, H) ) return 0;
	if ( index <= 1 ) return 1;		// treat index 0 as 1 (default)
	if ( H->v % index ) { err_printf ("The prime p=%ld does not split completely in the ring class field of discriminant %ld^2*%ld\n", _ff_p, index, H->D); return 0; }
	v = H->v / index;
	u = o->u;
	for ( i = 0 ; i < u->k ; i++ ) {
		for ( d = u->h[i], q = u->p[i] ; ! (v%q) ; d++, q *= u->p[i] );					// compute depth d of u->p[i] volcano
		if ( H->u ) {
			for ( k = 0 ; k < H->u->k ; k++ ) if ( H->u->p[k] == u->p[i] ) break;
			if ( k < H->u->k ) k = H->u->h[k]; else k = 0;
		} else {
			k = 0;
		}
		for ( e = k ; e < k+ u->h[i] ; e++ ) phi_descend (r, r[0], e, d, u->p[i], H->inv);
		dbg_printf ("Descended to j=%ld at level %d of %d-volcano of depth %d\n", _ff_get_ui(r[0]), e, u->p[i], d);
	}
	return 1;
}


int ff_classpoly_suborder_jinv (ff_t j[1], classpoly_t H, long index)
{
	struct suborder_struct *o;
	struct ifactors_struct *u;
	long q, v;
	int i, k, d, e;
	
	if ( index < 0 ) { err_printf ("Invalid index=%ld\n", index); abort(); }
	if ( ! (o = classpoly_get_suborder (H, index)) ) return 0;
	if ( ! ff_classpoly_jinv (j, H) ) return 0;
	if ( index <= 1 ) return 1;		// treat index 0 as 1 (default)
	if ( H->v % index ) { err_printf ("The prime p=%ld does not split completely in the ring class field of discriminant %ld^2*%ld\n", _ff_p, index, H->D); return 0; }
	v = H->v / index;
	u = o->u;
	for ( i = 0 ; i < u->k ; i++ ) {
		for ( d = u->h[i], q = u->p[i] ; ! (v%q) ; d++, q *= u->p[i] );					// compute depth d of u->p[i] volcano
		if ( H->u ) {
			for ( k = 0 ; k < H->u->k ; k++ ) if ( H->u->p[k] == u->p[i] ) break;
			if ( k < H->u->k ) k = H->u->h[k]; else k = 0;
		} else {
			k = 0;
		}
		for ( e = k ; e < k+ u->h[i] ; e++ ) phi_descend (j, j[0], e, d, u->p[i], INV_J);
		dbg_printf ("Descended to j=%ld at level %d of %d-volcano of depth %d\n", _ff_get_ui(j[0]), e, u->p[i], d);
	}
	return 1;
}


ff_t *ff_classpoly_jinvs (classpoly_t H)
{
	ff_t r;

	if ( ! H->C ) { err_printf ("Call to ff_classpoly_jinvs with no classpoly loaded!\n"); abort(); }
	if ( ! ff_classpoly_jinv (&r, H) ) { err_printf ("Call to ff_classpoly_jinv failed unexpectedly in ff_classpoly_jinvs\n"); return 0; }
	return ff_classpoly_enum_jinvs (H, r);
}

ff_t *ff_classpoly_suborder_jinvs (classpoly_t H, long u)
{
	ff_t j;
	
	if ( ! ff_classpoly_suborder_jinv (&j, H, u) ) return 0;
	dbg_printf ("ff_classpoly_suborder_jinv(%ld,%ld) mod %ld returned j-invariant %ld\n", H->D, u, _ff_p, _ff_get_ui(j));
	return ff_classpoly_suborder_enum_jinvs (H, u, j);
}

// Note this currently only works for enum invariants (i.e., those for which we can compute modpolys), ff_classpoly_enum_roots will fail otherwise
ff_t *ff_classpoly_roots (classpoly_t H)
{
	ff_t r;

	if ( ! H->C ) { err_printf ("Call to ff_classpoly_jinvs with no classpoly loaded!\n"); abort(); }
	if ( ! ff_classpoly_root (&r, H) ) { err_printf ("Call to ff_classpoly_jinv failed unexpectedly in ff_classpoly_jinvs\n"); return 0; }
	return ff_classpoly_enum_roots (H, r, H->inv);
}

ff_t *ff_classpoly_suborder_roots (classpoly_t H, long u)
{
	ff_t r;
	
	if ( ! ff_classpoly_suborder_root (&r, H, u) ) return 0;
	dbg_printf ("ff_classpoly_suborder_root(%ld,%ld) mod %ld returned j-invariant %ld\n", H->D, u, _ff_p, _ff_get_ui(r));
	return ff_classpoly_suborder_enum_roots (H, u, r, H->inv);
}

int ff_classpoly_isroot (classpoly_t H, ff_t r)
{
	ff_t y;
	
	if ( ! H->F_d ) { err_printf ("Call to ff_classpoly_isroot with no classpoly loaded!\n"); abort(); }
	if ( H->G_d > 1 ) { err_printf ("Call to ff_classpoly_isroot with G_d > 1!\n"); abort(); }
	ff_classpoly_reduce (H);
	ff_poly_eval (&y, H->c, H->F_d, &r);
	return ( _ff_zero(y) ? 1 : 0 );
}


int ff_classpoly_find_jinv (ff_t J[1], classpoly_t H)
{
	struct phi_vshape_struct s;
	ff_t f[4];
	int j_filter;
	int i, j;
	
	if ( ! ff_classpoly_setup (H) ) return 0;
	classpoly_setup_find_jinv (H);
	if ( H->D == -3 ) { _ff_set_zero(J[0]); return 1; }
	if ( H->D == -4 ) { _ff_set_ui(J[0], 1728); return 1; }
	j_filter =  (_ff_p1mod3&&(H->v%3)&&((-H->D)%3) ? FINDCURVE_FILTER_JCUBE : 0 );
	phi_vshape (&s, H->v);
	for (;;) {
		if ( ! findcurve (f, H->t, j_filter, 0) ) { err_printf ("find_curve failed with p=%ld, t=%ld\n", _ff_p, H->t); abort(); }
		if ( _ff_zero (f[0]) || _ff_zero(f[1]) ) continue;
		ecurve_to_jinv(J,f);
		dbg_printf ("Found curve with j-invariant %ld over F_%ld with trace %ld and v=%ld\n", _ff_get_ui(J[0]), _ff_p, H->t, H->v);
		for ( i = j = 0 ; (i < H->u->k && H->u->p[i] <= PHI_MAX_M) || j < s.k ; ) {
			if ((i < H->u->k && H->u->p[i] <= PHI_MAX_M) && (j == s.k || H->u->p[i] < s.p[j]) ) {
				phi_adjust_level (J, J[0], H->u->h[i], H->u->h[i], H->u->p[i], 0);  i++;
			} else if ( (j < s.k && s.p[j] <= PHI_MAX_M) && (i==H->u->k || s.p[j] < H->u->p[i]) ) {
				phi_adjust_level (J, J[0], 0, s.h[j], s.p[j], 0);  j++;
			} else {
				if ( H->u->p[i] != s.p[j] || i >= H->u->k || i > PHI_MAX_M || j >= s.k ) { err_printf ("Bug, factors of u and v don't match up!\n"); exit (0); }
				phi_adjust_level (J, J[0], H->u->h[i], H->u->h[i]+s.h[j], H->u->p[i], 0);  i++; j++;
			}
		}
		dbg_printf ("Curve with j-invariant %ld over F_%ld is at the correct level for ell <= %d\n", _ff_get_ui(J[0]), _ff_p, PHI_MAX_M);
		for ( i = 0 ; i < H->H1_len ; i++ ) {
			if ( ff_classpoly_isroot(H->H1+i, J[0]) ) { dbg_printf ("invariant %ld is a root of H_%ld (wrong endomorphism ring), retrying (without filter)\n", _ff_get_ui(J[0]), H->H1[i].D); j_filter = FINDCURVE_NOFILTER;  break; }
			dbg_printf ("j-invariant %ld is not a root of H_D1 for D1=%ld\n", _ff_get_ui(J[0]), H->H1[i].D);
		}
		if ( i == H->H1_len ) break;
	}
	return 1;
}


ff_t *ff_classpoly_enum_roots (classpoly_t H, ff_t r, int enum_inv)
{
	long cnt;
	int i;
	
	if ( ! ff_classpoly_setup (H) ) return 0;
	if ( ! H->pres || H->pres->enum_inv != enum_inv ) { err_printf ("you must call classpoly_setup_enum_roots prior to calling classpoly_enum_roots (enum_inv=%d)\n", enum_inv); abort(); }
	dbg_printf ("enumerating roots for D=%ld mod %ld starting from %ld using inv %d\n", H->D, _ff_p, _ff_get_ui(r), H->pres->enum_inv);
	i = (H->pres->ell0?1:0);				// note that the 2-d array aux_e is already shifted to handle ell0 > 0, we only shift the linear lists (yes this is a bit ugly...)
	cnt = phi_enum_roots (H->roots, r, H->pres->ell+i, H->pres->n+i, H->pres->o+i, H->pres->aux_ell+i, H->pres->aux_e, H->pres->k-i, H->v, H->pres->enum_inv);
	if ( cnt != H->pres->enum_cnt ) { err_printf ("Only enumerated %ld of %ld expected j-invariants\n", cnt, H->pres->enum_cnt); abort(); }
	H->roots_p = _ff_p;
	return H->roots;
}


ff_t *ff_classpoly_suborder_enum_roots (classpoly_t H, long u, ff_t r, int enum_inv)
{
	struct suborder_struct *o;
	long cnt;

	if ( ! ff_classpoly_setup (H) ) return 0;
	if ( ! (o = classpoly_get_suborder (H, u)) ) return 0;
	if ( ! o->pres || o->pres->enum_inv != enum_inv ) { err_printf ("you must call classpoly_suborder_setup_enum_roots prior to calling classpoly_suborder_setup_enum_roots (enum_inv=%d)\n", enum_inv); abort(); }
	dbg_printf ("(suborder) enumerating roots for D=%ld mod %ld starting from %ld using inv %d\n", o->D, _ff_p, _ff_get_ui(r), enum_inv);
	cnt = phi_enum_roots (o->roots, r, o->pres->ell, o->pres->n, o->pres->o, o->pres->aux_ell, o->pres->aux_e, o->pres->k, H->v, o->pres->enum_inv);
	if ( cnt != o->pres->enum_cnt ) { err_printf ("Only enumerated %ld of %ld expected j-invariants\n", cnt, o->pres->enum_cnt); abort(); }
	o->roots_p = _ff_p;
	return o->roots;
}


int torsor_setup (torsor_t T, classpoly_t H, long u)
{
	struct suborder_struct *o;
	struct classgroup_pcp_struct *pres;
	ff_t *roots;
	long e1[IQ_MAX_GENS];
	register int i, j;
	
	if ( u ) {
		if ( ! (o = classpoly_get_suborder (H, u)) ) return 0;
		if ( o->roots_p != _ff_p ) { err_printf ("You must enumerate roots before attempting to create a torsor (D=%ld, u=%ld, _ff_p=%ld, roots_p=%ld)\n", H->D, u, _ff_p, o->roots_p); abort(); }
		pres = o->pres;  roots = o->roots;
	} else {
		pres = H->pres;  roots = H->roots;
		if ( H->roots_p != _ff_p ) { err_printf ("You must enumerate roots before attempting to create a torsor (D=%ld, _ff_p=%ld, roots_p=%ld)\n", H->D, _ff_p, H->roots_p); abort(); }
	}
	if ( ! pres || ! pres->orientable ) { err_printf ("You must call setup and orient the presentation before attempting to create a torsor\n"); abort(); }
	if ( dbg_level >= DEBUG_LEVEL ) { printf ("torsor enumeration: "); for ( i = 0 ; i < pres->enum_cnt ; i++ ) printf ("%ld  ", _ff_get_ui(roots[i])); puts (""); }
	
	assert ( pres->k <= IQ_MAX_GENS );
	T->pres = pres;  T->roots = roots;
	for ( i = 0 ; i < pres->k ; i++ ) {
		// if generator doesn't require orientation, just copy its power relation and continue
		if ( pres->orient_p[i] <= 0 ) { T->signs[i] = 1; for (  j = 0 ; j < i ; j++ ) evec_ri(T->r,i)[j] = evec_ri(pres->r,i)[j]; continue; }
		// get rep of orientation element and express it in terms of the (partially) oriented presentation
		for ( j = 0 ; j < i ; j++ ) e1[j] = ( T->signs[j] < 0 ? pres->o[j] - pres->orient_reps[i][j] : pres->orient_reps[i][j] );  e1[j] = pres->orient_reps[i][j];  for ( j++ ;  j < pres->k ; j++ ) e1[j] = 0;
		evec_reduce (e1, pres->n, T->r, pres->k);							
		j = evec_to_index (e1,pres->m,pres->k);
		if ( pres->orient_q[i] > 1 ) {
			dbg_printf ("Verifying (%ld,%ld)-path between [%d]=%ld and [%d]=%ld to orient ell=%ld...", pres->orient_p[i], pres->orient_q[i], 0, _ff_get_ui(roots[0]), j, _ff_get_ui(roots[j]), pres->ell[i]);
			T->signs[i] = ( phi_poly_verify_2path (roots[0], roots[j], pres->orient_p[i], pres->orient_q[i], pres->enum_inv) ? 1 : -1 );
		} else {
			dbg_printf ("Verifying %ld-edge between [%d]=%ld and [%d]=%ld to orient %ld...", pres->orient_p[i], 0, _ff_get_ui(roots[0]), j, _ff_get_ui(roots[j]), pres->ell[i]);
			T->signs[i] = (  phi_poly_verify_edge (roots[0], roots[j], pres->orient_p[i], pres->enum_inv) ? 1 : -1 );
		}
		// update power relation (this is necessary even when sign is 1, since earlier elements may have changed sign)
		for ( j = 0 ; j < i ; j++ ) e1[j] = ( T->signs[i]*T->signs[j] < 0 ? pres->o[j] -evec_ri(pres->r,i)[j] : evec_ri(pres->r,i)[j] );  while ( j < pres->k ) e1[j++] = 0;
		evec_reduce (e1, pres->n, T->r, pres->k);
		for ( j = 0 ; j < pres->k ; j++ ) evec_ri(T->r,i)[j] = e1[j];
		dbg_printf ("sign for generator with norm %ld is %d\n", pres->ell[i], T->signs[i]);
			
		// The code below reverifies the oriented presentation as a sanity check, only for testing, this is not necessary and can eventually be removed
		for ( j = 0 ; j <= i ; j++ ) e1[j] = ( T->signs[j] < 0 ? pres->o[j] - pres->orient_reps[i][j] : pres->orient_reps[i][j] );  while ( j < pres->k ) e1[j++] = 0;
		evec_reduce (e1, pres->n, T->r, pres->k);							
		j = evec_to_index (e1,pres->m,pres->k);
		if ( pres->orient_q[i] > 1 ) {
			dbg_printf ("Reverifying (%ld,%ld)-path between [%d]=%ld and [%d]=%ld to orient ell=%ld...", pres->orient_p[i], pres->orient_q[i], 0, _ff_get_ui(roots[0]), j, _ff_get_ui(roots[j]), pres->ell[i]);
			if ( ! phi_poly_verify_2path (roots[0], roots[j], pres->orient_p[i], pres->orient_q[i], pres->enum_inv) ) { err_printf("Orientation failed\n"); abort(); }
		} else {
			dbg_printf ("Reverifying %ld-edge between [%d]=%ld and [%d]=%ld to orient %ld...", pres->orient_p[i], 0, _ff_get_ui(roots[0]), j, _ff_get_ui(roots[j]), pres->ell[i]);
			if ( ! phi_poly_verify_edge (roots[0], roots[j], pres->orient_p[i], pres->enum_inv) ) { err_printf("Orientation failed\n"); abort(); }
		}
		dbg_printf ("Orientation sign for ell=%ld confirmed.\n", pres->ell[i]);
		// end reverification code
	}
if ( dbg_level >= DEBUG_LEVEL ) {
for ( i = 0 ; i < pres->k ; i++ ) {
	printf ("Generator a[%d] with norm %ld has orientation %d\n", i, pres->ell[i], T->signs[i]);
	printf ("   original relation: "); evec_print_relation (pres->n, pres->r, i);
	printf ("   oriented relation: "); evec_print_relation (pres->n, T->r, i);
}
}
	return 1;
}


// We want to avoid situations where alpha_i^{+/-2} = alpha_j^2 (or = alpha_0*alpha_j^2 if ell0 flag is set), with |alpha_i| = |alpha_j| = 4 (or have 4th powers in <ell0> but not 2nd powers in <ell0>) and j < i
// These cases cause problems when enumerating roots via gcds
// returns the index of the first bad generator, or -1 if no bad generators are found
int classgroup_pcp_check_generators (classgroup_pcp_t pres)
{
	long *n, *r;
	long e1[IQ_MAX_GENS];
	register long *ei;
	register int i,i0,j,k,s;
	
	n = pres->n;  r = pres->r;  k = pres->k;  s = ( pres->ell0 ? 1 : 0 );
	for ( i = s+1 ; i < k ; i++ ) {
		if ( n[i] != 2 ) continue;
		ei = evec_ri(r,i);
		for ( j = s ; j < i ; j++ ) if ( ei[j] ) break;
		if ( j == i ) continue;
		for ( i0 = s ; i0 < i ; i0++ ) {
			if ( (4%n[i0]) ) continue;
			evec_clear(e1,k);  e1[i0] = 4;  evec_reduce(e1,n,r,k);
			for ( j = s; j < i ; j++ ) if ( e1[j] ) break;
			if ( j<i ) continue;								// alpha_i0^4 is not trivial or in <alpha_0>
			evec_clear(e1,k);  e1[i0] = 2;  evec_reduce(e1,n,r,k);	// compute alpha_i0^2
			for ( j =s ; j < i ; j++ ) if ( e1[j] != ei[j] ) break;
			if ( j== i ) return i;
			evec_inverse(e1,e1,n,r,k);						// compute alpha_i0^{-2}
			for ( j =s ; j < i ; j++ ) if ( e1[j] != ei[j] ) break;
			if ( j== i ) return i;
		}
	}
	return -1;
}


// This function must be called while qform_table has been filled by a call to qform_generators
void classgroup_pcp_setup_auxiliary_primes (classgroup_pcp_t pres, int enum_inv)		// enum_inv is only used to estimate costs
{
	long e1[IQ_MAX_GENS], e2[IQ_MAX_GENS];
	long p, wt, best_wt, worst_wt;
	long *ell, *n, *o, *aux_ell, ell0, D;
	register int i, j, k;

	ell = pres->ell;  n = pres->n;  o = pres->o;  aux_ell = pres->aux_ell;  ell0 = pres->ell0;  k = pres->k; D = pres->D;
	for ( i = 0 ; i < k ; i++ ) {
		aux_ell[i] = 0;
		if ( o[i] <= 2 ) continue;							// we don't need to compute either auxilliary ell or orienting elements for generators of order 2
		p = ell[i];
		// we don't know the volcano depth or what p will be, so assume trivial volcano (for ell > 2) and p ~ |D|, inv j, which is not always the best thing to do!	
		best_wt = worst_wt = poly_root_enum_cost (ell[i], (ell[i] > 2 ? 0 : 1), n[i], -D, enum_inv);
		info_printf ("root enum cost for ell=%ld with relative order %ld is %ld nsecs\n", p, n[i], best_wt);
		wt = 0;
		while ( wt < worst_wt ) {						// once we exceed worst_wt, we can't get back under it with a larger ell' (gcd step costs increase monotonically)
			p = qform_next_prime(p,D);
			if ( p > PHI_MAX_M ) break;
			if ( !(pres->ellfilter%p) || !(D%(p*p)) ) continue;
			if ( ! qform_prime_rep (e1, p, D) ) { err_printf("Couldn't get rep of auxiliary ell=%ld for D=%ld\n", p, D); exit (0); }
			if ( ! e1[i] ) continue;						// we want a non-trivial index for ell[i]
			for ( j = i+1 ; j < k ; j++ ) if ( e1[j] ) break;
			if ( j < k ) continue;						// we need zero indices for all ell[j] with j > i
			if ( ell0 && e1[0] ) continue;					// we don't want ell0 to appear in the rep
			if ( o[i] > 3 && n[i] > 2 ) {					// don't use auxiliary primes with order 3 elements or one > ell
				// switch to inverse to compute weight for aux if it is closer (but do not use this for orient!)
				if ( e1[i] > n[i]/2 ) {
					if ( ! qform_iprime_rep (e2, p, D) )  { err_printf("Couldn't get rep of inverse of auxiliary ell=%ld\n", p); exit (0); }
					if ( ell0 && e2[0] ) evec_copy (e2, e1, k);// switch back if this forces us to use ell0
				} else {
					evec_copy (e2, e1, k);
				}
				// currently enum_roots won't use auxilliary primes with exponents less than 3, so skip these
				if ( e2[i] > 2 && wt < worst_wt ) {
					wt = poly_gcd_enum_cost (ell[i], p, e2[i], n[i], 0, -D, enum_inv);								 //16*e2[i]*ell[i]*ell[i] + (n[i]-e2[i])*p*p;
					info_printf ("ell=%ld    aux: %ld = ", ell[i], p);
					for ( j = 0 ; j <= i ; j++ ) info_printf ("%ld^%ld ", ell[j], e2[j]); info_printf("(cost %ld nsecs)\n", wt);
					if ( wt < best_wt ) {
						best_wt = wt;
						aux_ell[i] = p;  for (  j = (ell0?1:0) ; j < k ; j++ ) pres->aux_e[(ell0?i-1:i)][(ell0?j-1:j)] = e2[j];		// shift aux reps down 1 to cover ell0 when set
					}
				}
			}
		}
		if ( aux_ell[i] ) { info_printf ("ell[%d]=%ld using aux_ell=%ld (cost %ld nsecs)\n", i, ell[i], aux_ell[i], best_wt); }
		else { info_printf ("ell[%d]=%ld not using aux_ell\n", i, ell[i]);}
	}
}


int classgroup_pcp_setup (classgroup_pcp_t pres, long D, long ellfilter, long ell0, int enum_inv)		// enum_inv is only used to estimate costs
{
	int i, j;

	if ( ! ellfilter ) ellfilter = 1;
	if ( enum_inv ) ellfilter = ui_lcm (ellfilter, inv_level(enum_inv));
	dbg_printf ("classgroup_pcp_setup D=%ld, ellfilter=%ld, ell0=%ld, enum_inv=%d\n", D, ellfilter, ell0, enum_inv);
	memset (pres, 0, sizeof(*pres));
	pres->D = D;
	pres->ellfilter = ellfilter;
	pres->ell0 = ell0;
	pres->enum_inv = enum_inv;
	dbg_printf ("Computing generators for cl(%ld)\n", pres->D);
	pres->k = qform_generators (pres->ell, pres->n, pres->r, pres->D, pres->ellfilter, pres->ell0, 0, 0);
	if ( ! pres->k ) { pres->h = pres->enum_cnt = 1; if ( pres->ell0 ) { err_printf ("ell0=%ld set for D=%ld with odd class number 1!\n", pres->ell0, pres->D); abort(); } return 0; }
	if ( pres->k <= 0 ) return pres->k;
	dbg_printf ("Checking for bad generators\n");
	while  ( (i = classgroup_pcp_check_generators (pres)) >= 0 ) {
		info_printf ("bad generator ell[%d]=%ld from presentation (a^2=b^2 problem), pres_ellfilter=%ld\n", i, pres->ell[i], pres->ellfilter);
		if  ( ui_len(pres->ell[i])+ui_len(pres->ellfilter) > 63 ) { err_printf ("Too many bad generators excluded from presentation, can't create polycyclic presentation for D=%ld\n", D); return -1; }
		pres->ellfilter *= pres->ell[i];
		pres->k = qform_generators (pres->ell, pres->n, pres->r, pres->D, pres->ellfilter, pres->ell0, 0, 0);
	}
	for ( i = 0 ; i < pres->k ; i++ ) if ( pres->ell[i] > PHI_MAX_M ) { err_printf ("Norms in classgroup presentation for D=%ld exceed PHI_MAX_M=%d\n", D, PHI_MAX_M); return -1; }
	evec_n_to_m (pres->m, pres->n, pres->k);
	evec_orders (pres->o, pres->n, pres->r, pres->k);
	pres->h = pres->n[0];  for ( i = 1 ; i < pres->k ; i++ ) pres->h *= pres->n[i];
	if ( ell0 ) { if ( pres->ell[0] != ell0 || pres->o[0] != 2 ) { err_printf ("ell0=%ld does not appear at index 0 with order 2!\n", ell0); return -1; } }
	pres->enum_cnt = ( ell0 ? pres->h/2 : pres->h );
	if ( dbg_level >= 2 ) {
		printf ("Computed polycyclic presentation for cl(%ld) with %d generators: \n", D, pres->k);
		for ( i = 0 ; i < pres->k ; i++ ) { printf ("  %ld^%ld = ", pres->ell[i], pres->n[i]);  for ( j = 0 ; j < i ; j++ ) printf ("%ld^%ld ", pres->ell[j], evec_ri(pres->r, i)[j]);  printf ("\t\t(order %ld)\n", pres->o[i]); }
	}
	dbg_printf ("Computing auxilliary primes for gcd enumeration\n");
	classgroup_pcp_setup_auxiliary_primes (pres, enum_inv);
	dbg_printf ("Computing orienting primes\n");
	if ( ! classgroup_pcp_orient (pres) ) out_printf ("Warning, unable to orient presentation for D=%ld\n", D);
	return pres->k;
}

// Select auxilliary primes (or prime pairs) to support orienting presentations, see Section 4.3 of "Class invariants by the CRT method".
int classgroup_pcp_orient (classgroup_pcp_t pres)
{
	register int i, j;
	long e1[IQ_MAX_GENS], e2[IQ_MAX_GENS], e3[IQ_MAX_GENS], e4[IQ_MAX_GENS], e5[IQ_MAX_GENS], e6[IQ_MAX_GENS];
	long ell, p, q, D, *ps, *qs;

	// mark elements that require orientation by setting orient_ell[i][0] to 0, otherwise set orient_ell[0][i] to -1
	ps = pres->orient_p;  qs = pres->orient_q;
	for ( i = 0 ; i < pres->k ; i++ ) { ps[i] = -1; if ( pres->o[i] > 2 ) break; }		// We only need to orient elments that have order > 2 and occur after the first element with order > 2 in the presentation
	for ( i++ ; i < pres->k ; i++ ) ps[i] = ( pres->o[i] > 2 ? 0 : -1 );

	D = pres->D;
		
	info_printf ("Determining orienting primes (or pairs of primes) for presentation for cl(%ld)\n", pres->D);
	pres->orientable = 0;
	for ( i = 0 ; i < pres->k ; i++ ) {
		if ( ps[i] ) continue;
		ell = pres->ell[i];
		p = ell;
		while ( ! ps[i] ) {
			p = qform_next_prime(p,D);
			if ( p > PHI_MAX_ORIENT_P ) break;
			if ( !(pres->ellfilter%p) || !(D%(p*p)) ) continue;
			if ( ! qform_prime_rep (e1, p, D) ) { err_printf("Couldn't get rep of primeform with norm %ld in cl(%ld)\n", p, D); abort(); }
			if ( ! e1[i] ) continue;						// we want a non-trivial index for ell[i]
			for ( j = i+1 ; j < pres->k ; j++ ) if ( e1[j] ) break;
			if ( j < pres->k ) continue;					// we need zero indices for all ell[j] with j > i
			if ( pres->ell0 && e1[0] ) continue;			// we don't want ell0 to appear in the rep

			// make sure that swapping ell[i] with its inverse in the rep of ells[i][0] yields an element distinct from both ells[i][0] *and* its inverse
			for ( j = 0 ; j < i ; j++ ) e2[j] = e1[j];  e2[j] = pres->o[j]-e1[i];  j++; while ( j < pres->k ) e2[j++] = 0;
			evec_reduce(e2,pres->n,pres->r,pres->k);
			for ( j = 0 ; j <= i ; j++ ) if ( e1[j] != e2[j] ) break;
			if ( j > i ) continue;
			for ( j = 0 ; j <= i ; j++ ) e3[j] = ( e1[j] ? pres->o[j]-e1[j] : 0 ); while ( j < pres->k ) e3[j++] = 0;
			evec_reduce(e3,pres->n,pres->r,pres->k);
			for ( j = 0 ; j <= i ; j++ ) if ( e2[j] != e3[j] ) break;
			if ( j > i ) continue;					
			ps[i] = p;  qs[i] = 1;
			for ( j = 0 ; j <= pres->k ; j++ ) pres->orient_reps[i][j] = e1[j];
			info_printf ("ell=%ld    orienting norm: %ld = ", ell, p);
			for ( j = 0 ; j <= i ; j++ ) info_printf ("%ld^%ld ", pres->ell[j], e1[j]);  info_printf ("\n");
		}
		if ( ps[i] ) continue;
		// if we couldn't orient ell[i] with a single prime, try using a product of primes p*q, with p>q.
		p = qform_next_prime(ell,D);
		while ( ! ps[i] ) {
			p = qform_next_prime(p,D);
			if ( p > PHI_MAX_ORIENT_P ) break;
			if ( !(pres->ellfilter%p) || !(D%(p*p)) ) continue;
			if ( ! qform_prime_rep (e4, p, D) ) { err_printf("Couldn't get rep of primeform with norm %ld in cl(%ld)\n", p, D); abort(); }
			if ( ! qform_iprime_rep (e5, p, D) ) { err_printf("Couldn't get rep of primeform with norm %ld in cl(%ld)\n", p, D); abort(); }
			q = ell;
			while ( ! qs[i] ) {
				q = qform_next_prime(q,D);
				if ( q >= p ) break;
				if ( !(pres->ellfilter%q) || !(D%(q*q)) ) continue;
				if ( ! qform_prime_rep (e3, q, D) ) { err_printf("Couldn't get rep of auxiliary ell=%ld\n", p); exit (0); }
				evec_compose (e1,e3,e4,pres->n,pres->r,pres->k);
				if ( ! e1[i] ) continue;						// we want a non-trivial index for ell[i]
				for ( j = i+1 ; j < pres->k ; j++ ) if ( e1[j] ) break;
				if ( j < pres->k ) continue;					// we need zero indices for all ell[j] with j > i
				if ( pres->ell0 && e1[0] ) continue;			// we don't want ell0 to appear in the rep
				
				// we now need to ensure that swapping ell[i] with its inverse in the rep of p*q yields an element distinct from p*q, p*q^-1, p^-1*q, and p^-1*q^-1
				for ( j = 0 ; j < i ; j++ ) e2[j] = e1[j];  e2[j] = pres->o[j]-e1[i];  j++; while ( j < pres->k ) e2[j++] = 0;
				evec_reduce(e2,pres->n,pres->r,pres->k);
				if ( evec_equal(e1,e2,pres->k) ) continue;
				evec_inverse_o(e6,e1,pres->n,pres->o,pres->r,pres->k);
				if ( evec_equal(e6,e2,pres->k) ) continue;
				evec_compose(e6,e3,e5,pres->n,pres->r,pres->k);
				if ( evec_equal(e6,e2,pres->k) ) continue;
				evec_inverse_o(e6,e6,pres->n,pres->o,pres->r,pres->k);
				if ( evec_equal(e6,e2,pres->k) ) continue;
				ps[i] = p; qs[i] = q;   for ( j = 0 ; j <= pres->k ; j++ ) pres->orient_reps[i][j] = e1[j];
				info_printf ("ell=%ld    orienting norm: %ld*%ld = ", ell, p, q);
				for ( j = 0 ; j <= i ; j++ ) info_printf ("%ld^%ld ", pres->ell[j], e1[j]);  info_printf ("\n");
			}
		}
		if ( ! ps[i] ) { out_printf ("Couldn't find an orienting prime (or pair of primes) for D=%ld, ell=%ld with PHI_MAX_ORIENT_P=%d\n", D, ell, PHI_MAX_ORIENT_P); return 0; }
	}
	pres->orientable = 1;
	return 1;
}


void classgroup_pcp_print (classgroup_pcp_t pres)
{
	long a, b, c;
	register int i, j;

	for ( i = 0 ; i < pres->k ; i++ ) {
		qform_primeform(&a, &b, &c, pres->ell[i], pres->D, 0);
		info_printf ("%ld: (%ld,%ld,%ld) relative order %ld, order %ld, ", pres->ell[i], a, b, c, pres->n[i], pres->o[i]);
		for ( j = 0 ; j< i ; j++ ) {
			if ( ! j ) { info_printf ("%ld^%ld = %ld^%ld ", pres->ell[i], pres->n[i], pres->ell[j], evec_ri(pres->r,i)[j]); }
			else { info_printf ("%ld^%ld ", pres->ell[j], evec_ri(pres->r,i)[j]); }
		}
		info_printf ("\n");
	}
}


long classgroup_pcp_vfilter (classgroup_pcp_t pres)
{
	int i;
	long vfilter;
	
	// we really only care about ell[0], and possibly ell[1] if  ell[0] is the norm of an element with order 2
	vfilter = 1;
	i = ( pres->ell0 ? 1 : 0 );
	if ( pres->k > i ) {
		if ( pres->ell[i] > 3 && (vfilter % pres->ell[i]) ) vfilter *= pres->ell[i];
		if ( pres->k > i+1 ) {
			if ( pres->n[i]==2 && pres->ell[i+1] > 3 && (vfilter % pres->ell[i+1]) ) vfilter *= pres->ell[i+1];
		}
	}
	return vfilter;
}
