#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <gmp.h>
#include "ff_poly.h"
#include "ff_poly/ffpolysmall.h"
#include "mpzutil.h"
#include "phi_eval.h"
#include "phi_gcd.h"
#include "polycosts.h"
#include "phi_poly.h"
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

#define PHI_MAX_COEFF		(((PHI_MAX_M+1)*(PHI_MAX_M+2))/2)

#define PHI_SMALL_M		31

char _phi_dir_str[4096];

/*
	A nonzero table entry m in row ell and column h indicates that if ell* is greater than m than we should just use Phi_ell (i.e. find roots)
	and not compute gcds with Phi_ell and Phi_ell_aux.  The default is the later.

	this is obsolete, superceded by polycosts.h, to be removed
*/
/*
static unsigned char phi_maxgcd[PHI_SMALL_M+1][PHI_MAX_ENUM_H+1] = {
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 13, 13, 17, 19, 19, 23, 23, },
{ 17, 29, 37, 43, 53, 53, 61, 67 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 31, 67, 89, 107, 113, 131, 149, 157 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 43, 101, 137, 163, 191, 211, 229, 241 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 71, 191, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 83, 241, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 103, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 113, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 139, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 181, 0, 0, 0, 0, 0, 0, 0 },
{ 0, 0, 0, 0, 0, 0, 0, 0 },
{ 193, 0, 0, 0, 0, 0, 0, 0 },
};
*/

// currently we don't free any of these polys once they get loaded.  Even though INV_MAX is big (1000), only one or two invariants are typically used in any particular execution
static mpz_t **Phi_mpz[INV_MAX+1];//[PHI_MAX_M+1];	// For prime m<= PHI_MAX_M, the mth entry points to list of coefficients for Phi_m(X,Y)
static ff_t **Phi_ff[INV_MAX+1];//[PHI_MAX_M+1];
static ff_t *Phi_ff_p[INV_MAX+1];//[PHI_MAX_M+1];		// holds last value of p used in reduction, allows reuse 
mpz_t *Phi_mpz_mod[PHI_MAX_M+2];			 		// only classical modular poly (INV_J)

struct phi_vshape_struct phi_vshapes[PHI_MAX_V+1];
int phi_maxv;

int phi_find_surface_pp (ff_t *sJ, ff_t phi[], int m, int h, ff_t J);
int phi_find_surface_pq (ff_t *sJ, ff_t phi1[], int m1, ff_t phi2[], int m2, ff_t J);
int phi_next_surface_neighbor (ff_t nJ[1], ff_t phi[], int m, int h, ff_t J, ff_t *pJ, int inv);


void _phi_vshape (struct phi_vshape_struct *s, int v)
{
	unsigned long p[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS];
	int i, k;

	k = ui_factor (p, h, v);
	for ( i = 0 ; i < k ; i++ ) {
		if ( p[i] > PHI_MAX_M ) break;
		if ( i >= PHI_MAX_V_P ) { err_printf ("v=%d has more than PHI_MAX_V_P=%d prime factors\n", v, PHI_MAX_V_P); abort(); }
		s->p[i] = (unsigned char) p[i];  s->h[i] = (unsigned char) h[i];
	}
	s->k = i;
}


int phi_poly_setup (void)
{
	int v;
	
	for ( v = 2 ; v <= PHI_MAX_V ; v++ ) _phi_vshape (phi_vshapes+v, v);
	phi_maxv = PHI_MAX_V;	
	return 1;
}

void phi_vshape (struct phi_vshape_struct *s, int v)
{
	if ( ! phi_maxv ) phi_poly_setup();
	if ( v > 0 && v <= phi_maxv ) *s = phi_vshapes[v]; else _phi_vshape (s, v);
}


void phi_poly_init (phi_poly_t phi, int m, int inv)
{
	int i, n;
	
	phi->m = m;
	phi->d = phi_degree (m);
	n = (phi->d*(phi->d+1))/2;
	phi->C = mem_alloc (n*sizeof(*phi->C));
	for ( i = 0 ; i < n ; i++ ) mpz_init (phi->C[i]);
	n = (phi->d+1)*(phi->d+1);
	phi->c = mem_alloc (n*sizeof(*phi->c));
	phi->w1 = mem_alloc ((phi->d+1)*sizeof(*phi->w1));
	phi->w2 = mem_alloc ((phi->d+1)*sizeof(*phi->w2));
	phi->p = 0;
	phi->inv = inv;
}


void phi_poly_clear (phi_poly_t phi)
{
	int i;
	
	for ( i = 0 ; i < (phi->d*(phi->d+1))/2 ; i++ ) mpz_clear (phi->C[i]);
	mem_free (phi->C);  mem_free (phi->c);  mem_free (phi->w1);  mem_free (phi->w2);
	memset (phi, 0, sizeof (*phi));
}


int phi_poly_create (phi_poly_t phi, int m, int inv)
{
	char buf[256];
	
	phi_poly_init (phi, m, inv);
	if ( ! phi_load_mpz (phi->C, m, inv) ) {
		if ( ui_is_prime (m) ) {
			sprintf (buf, "modpoly %d %d", m, inv);
		} else {
			sprintf (buf, "modpolycomp %d %d", m, inv);
		}
		if ( system(buf) < 0 ) { err_printf ("Error creating modular polynomial for m=%d, inv=%d\n", m, inv); phi_poly_clear (phi); return 0; }
		if ( ! phi_load_mpz (phi->C, m, inv) ) {  err_printf ("phi_load_mpz failed unexpectedly with m=%d, inv=%d\n", m, inv); phi_poly_clear (phi); return 0; }
	}
	return 1;
}


int phi_poly_load (phi_poly_t phi, int m, int inv)
{
	phi_poly_init (phi, m, inv);
	if ( ! phi_load_mpz (phi->C, m, inv) ) {  err_printf ("Error loading modular polynomial for m=%d, inv=%d\n", m, inv); phi_poly_clear (phi); return 0; }
	return 1;
}

/*
	Loads the coefficients of the modular polynomial phi_m(x,y) for the invariant f from the ascii file phif_m.txt.
	These polynomials relate f(z) to f(mz) for some modular function f which yields class invariants.
	Symmetry is required. 
        Here m is *not* assumed to be prime (so the degree need not be m+1, it is determined by calling phi_degree).
	(but everywhere else in this module it is assumed to be prime so that Phi_m has degree m+1)

	Supports both old and new formats:

	Old format:
		There are (m+1)(m+2)/2 coefficients to load (the terms x^(m+1) and y^(m+1) are implicit, and
		the coefficient for x^j*y^k is the same as that for x^k*y^j).
		The value P[i] will be set to the ith line of phi_m.txt (starting at i=0).  The coefficient for x^j*y^k with
		j >= k will be in P[i] where i=j(j+1)/2+k.
	New format:
		Every coefficient is listed, up to symmetry (which is assumed), in the form [j,k] c, meaning that c is
		the coefficient of x^j*y^k, where j >= k
*/
int phi_load_mpz (mpz_t Phi[], int m, int inv)
{
	clock_t start;
	FILE *fp;
	char filename[256];
	char buf[65536];
	register char *s;
	register int i,j,k,n,d;
	int type, phi_cnt, sts;
	
	// look in local directory first
	sprintf (filename, "phi_%s_%d.txt", inv_string(inv), m);
	fp = fopen(filename, "r");
	if ( ! fp ) {
		sprintf (filename, "%s/phi_%s_%d.txt", PHI_DIR, inv_string(inv), m);
		fp = fopen(filename, "r");
		if ( ! fp ) { err_printf ("Unable to locate file %s for invariant %d (%s)\n", filename, inv, inv_string(inv));  return 0; }
	}
	if ( ! fgets(buf,sizeof(buf),fp) ) { err_printf ("Error reading file %s\n", filename); fclose (fp); return 0; }
	start=0; // shut up compiler
	if ( m > PHI_MAX_M ) { out_printf ("Loading %s\n", filename);  start=clock(); }
	d = phi_degree (m);
	
	sts = 0;
	
	// skip header records
	while ( isalpha(buf[0]) ) if ( ! fgets(buf,sizeof(buf),fp) ) { err_printf ("Error reading file %s\n", filename); goto done; }
	phi_cnt =(d*(d+1))/2;
	type = 0;
	for ( i = 0 ;; i++ ) {
		if ( buf[strlen(buf)-1] != '\n' ) { err_printf ("Record %d in %s too long for buffer in phi_load_mpz, increase buffer size\n", i, filename); goto done; }
		if ( buf[0]=='[' ) {
			if ( i && !type ) { err_printf ("Inconsistent record format at line %d in %s\n", i, filename); goto done; }
			if ( ! i ) type = 1;	// note that we rely on Phi[i] having all its entries initialized to zero for type 1
			j = atoi(buf+1);
			for ( s = buf ; *s && *s != ',' ; s++ );
			if ( *s != ',' ) {err_printf ("Invalid record format at line %d in %s\n", i, filename); goto done; }
			k = atoi(s+1);
			if ( j < 0 || j > d || k < 0 || k > d ) {err_printf ("Invalid record exponent vector [%d,%d] at line %d in %s\n", j, k, i, filename); goto done; }
			for ( ; *s && *s != ']' ; s++ ); s++;
			if ( j==d && k == 0 ) {
				// ignore X^d and Y^d coefficients, which we know must be 1
				if ( atoi(s) != 1 ) { err_printf ("Invalid record format at line %d in %s\n", i, filename); goto done; }
			} else {
				n = phi_index(j,k);
				if ( n < 0 || n >= phi_cnt ) {err_printf ("Invalid record exponent vector [%d,%d] at line %d in %s\n", j, k, i, filename); goto done; }
				mpz_set_str(Phi[n],s,10);
			}
		} else {
			if ( type ) {err_printf ("Inconsistent record format at line %d in %s\n", i, filename); goto done; }
			mpz_set_str(Phi[i],buf,10);
		}
		if ( !  fgets(buf,sizeof(buf),fp) ) break;
	}
	if ( ! type && i+1 != phi_cnt ) { err_printf ("Only found %d lines in phi_%d.txt, expected %d\n", i, m, phi_cnt); goto done; }
	if ( m > PHI_MAX_M ) { out_printf ("Loaded %s in %.3f secs\n", filename, (double)(clock()-start)/CLOCKS_PER_SEC); }
	sts = 1;
done:	
	fclose (fp);
	return sts;
}

// for phi_ff we store a full (d+1) x (d+1) array for fast dot products
// only supports new format
int phi_load_ff (ff_t phi[], int m, int inv)
{
	FILE *fp;
	char filename[256];
	char buf[65536];
	mpz_t X;
	register char *s;
	register int i,j,k,n,d;
	int sts;
	
	// look in local directory first
	sprintf (filename, "phi_%s_%d.txt", inv_string(inv), m);
	fp = fopen(filename, "r");
	if ( ! fp ) {
		sprintf (filename, "%s/phi_%s_%d.txt", PHI_DIR, inv_string(inv), m);
		fp = fopen(filename, "r");
		if ( ! fp ) { err_printf ("Unable to locate file %s for invariant %d (%s)\n", filename, inv, inv_string(inv));  return 0; }
	}
	if ( ! fgets(buf,sizeof(buf),fp) ) { err_printf ("Error reading file %s\n", filename); fclose (fp); return 0; }
	info_printf ("Loading %s\n", filename);
	d = phi_degree (m);
	
	mpz_init (X);
	sts = 0;
	
	n = d+1;
	// skip header records
	while ( isalpha(buf[0]) ) if ( ! fgets(buf,sizeof(buf),fp) ) { err_printf ("Error reading file %s\n", filename); goto done; }
	for ( i = 0 ;; i++ ) {
		if ( buf[strlen(buf)-1] != '\n' ) { err_printf ("Record %d in %s too long for buffer in phi_load_mpz, increase buffer size\n", i, filename); goto done; }
		if ( buf[0]!='[' ) { err_printf ("Unsupported record format at line %d in %s\n", i, filename); goto done; }
		j = atoi(buf+1);
		for ( s = buf ; *s && *s != ',' ; s++ );
		if ( *s != ',' ) {err_printf ("Invalid record format at line %d in %s\n", i, filename); goto done; }
		k = atoi(s+1);
		if ( j < 0 || j > d || k < 0 || k > d ) {err_printf ("Invalid record exponent vector [%d,%d] at line %d in %s\n", j, k, i, filename); goto done; }
		for ( ; *s && *s != ']' ; s++ );  s++;
		mpz_set_str(X,s,10);
		_ff_set_mpz (phi[j*n+k], X);
		_ff_set (phi[k*n+j],phi[j*n+k]);		// set symmetric coefficient also
 		if ( !  fgets(buf,sizeof(buf),fp) ) break;
	}
	sts = 1;
done:	
	mpz_clear (X);
	fclose (fp);
	return sts;
}


// for phi_ff we store a full (m+1) x (m+1) array for fast dot products
void _phi_reduce (int m, int inv)
{
	struct phi_header *hdr;
	register int i, j, k, n;
	
	if ( m <= 0 || m > PHI_MAX_M ) { err_printf ("Invalid m=%d in _phi_reduce\n", m); abort(); }
	if ( ! Phi_ff_p[inv] ) Phi_ff_p[inv] = calloc(PHI_MAX_M+1,sizeof(ff_t));
	if ( ! Phi_ff[inv] ) Phi_ff[inv] = calloc(PHI_MAX_M+1,sizeof(ff_t));
	if ( ! Phi_mpz[inv] ) Phi_mpz[inv] = calloc(PHI_MAX_M+1,sizeof(mpz_t));

	n = phi_count(m);
	if ( ! Phi_mpz[inv][m] ) {
		Phi_mpz[inv][m] = malloc(n * sizeof(mpz_t));
		hdr = malloc((m+1)*(m+1)*sizeof(ff_t)+sizeof(struct phi_header));
		hdr->inv = inv;  hdr->sparse = inv_sparse_factor(inv);
		Phi_ff[inv][m] = (ff_t *) (hdr+1);
		for ( i = 0 ; i < n ; i++ ) mpz_init(Phi_mpz[inv][m][i]);
		if ( ! phi_load_mpz(Phi_mpz[inv][m], m, inv) ) { printf ("Couldn't load required phi_%d coefficients\n", m); abort(); }
	}
	for ( i = 0, k = 0 ; i <= m ; i++ ) {
		for ( j = 0 ; j <= i ; j++, k++ ) {
			_ff_set_mpz(Phi_ff[inv][m][(m+1)*i+j], Phi_mpz[inv][m][k]);
			if( j < i ) _ff_set(Phi_ff[inv][m][(m+1)*j+i], Phi_ff[inv][m][(m+1)*i+j]);
		}
	}
	if ( k != n ) { printf ("error, phi_count mismatch in _phi_reduce\n"); abort(); }
/*
	We use far more memory than is needed when Phi is sparse, but this doesn't seem to impact performance and would be a pain to change (see phi_eval.h)
	
	s = phi_sparse_factor(Phi_ff[inv][m]);
	if ( s > 1 ) {
		n = m/s+1;
		for ( i = 0 ; i <= m ; i++ ) {
			j0 = (m*(m+1-i))%s;												// j0 is the least b satisfying a + m*b = m+1 mod s (note that m^2=1 mod s for every valid s|24).
			for ( j=j0, k=0 ; j <= m ; j+= s, k++ ) _ff_set(phi_s[n*i+k], Phi_ff[inv][m][(m+1)*i+j]);
		}
	}
*/
	Phi_ff_p[inv][m] = _ff_p;
}

// reduces the coefficients of phi_m(x,y) in Phi_mpz mod p into Phi_ff, where p is the characteristic of the current finite field
static inline void phi_reduce (int m, int inv)
	{ if ( m < 0 || m > PHI_MAX_M || ! Phi_ff_p[inv] || _ff_p != Phi_ff_p[inv][m] ) _phi_reduce(m,inv); }

ff_t *phi_poly_ff (int m, int inv) 
	{ phi_reduce(m, inv); return Phi_ff[inv][m]; }

// reduces the coefficients of phi_m(x,y) in Phi_mpz mod P into Phi_mpz_mod (hardwired for j with m prime)
void phi_reduce_mpz_mod (int m, mpz_t P)
{
	register int i, n;
	
	if ( m <= 0 || m > PHI_MAX_M ) { err_printf ("Invalid m=%d in phi_reduce_mpz_mod\n", m); abort(); }
	if ( ! Phi_mpz[INV_J] ) Phi_mpz[INV_J] = calloc (PHI_MAX_M+1,sizeof(mpz_t));
	n = phi_count(m);
	if ( ! Phi_mpz[INV_J][m] ) {
		Phi_mpz[INV_J][m] = malloc(((m+1)*(m+2))/2 * sizeof(mpz_t));
		Phi_mpz_mod[m] = malloc(((m+1)*(m+2))/2 * sizeof(mpz_t));
		for ( i = 0 ; i < n ; i++ ) { mpz_init(Phi_mpz[INV_J][m][i]); mpz_init (Phi_mpz_mod[m][i]); }
		if ( ! phi_load_mpz(Phi_mpz[INV_J][m], m, INV_J) ) { printf ("Couldn't load required phi_%d coefficients\n", m); abort(); }
	}
	for ( i = 0 ; i < n ; i++ ) mpz_mod(Phi_mpz_mod[m][i], Phi_mpz[INV_J][m][i], P);
}

/*
	TODO: optimize the ell=2 case below to avoid computing the constant coefficient of Phi_2(X,J),
	since it is not needed to compute Phi_2(X,J)/(X-xJ) -- this will save 2 mults

	YES, DO THIS!!
*/

static inline int phi_get_neighbor (ff_t S[], ff_t phi[], int m, ff_t J, ff_t *xJ, int inv)
{
	ff_t f[PHI_MAX_M+2];

	if ( ! phi ) { phi_reduce(m,inv); phi=Phi_ff[inv][m]; }
	if ( m==2 && xJ ) {  			// Hardwire this case for speed
		ff_t t1, t2, D;
		phi2_eval_ff (f, phi, J);
		ff_poly_remove_root_d3 (f,f,xJ);
		_ff_square(t1,f[1]);
		_ff_mult(t2,f[2],f[0]);
		_ff_x2(t2);  _ff_x2(t2);
		_ff_subfrom(t1,t2);
		if ( ! ff_sqrt(&D,&t1) ) return 0;
		if ( ui_randomm(2) ) ff_negate(D);		// randomize choice of sqrt!
		_ff_set(t1,_ff_half);
		_ff_sub(t2,D,f[1]);
		_ff_mult(S[0],_ff_half,t2);
		return 1;
	}
	phi_eval_ff(f,phi,m,J);
// COMMENT
//for ( int i=0 ; i <= m ; i++ ) for (int j=0 ; j <= m ; j++ ) printf (" + %ld*X^%d*Y^%d", _ff_get_ui(phi[i*(m+1)+j]), i, j);  puts ("");
//printf ("phi_get_neighbor ell=%d j=%ld inv=%d, xJ=%ld\n", m, _ff_get_ui(J), inv, xJ ? _ff_get_ui(*xJ) : -1L);
//ff_poly_print (f, m+1);
	if ( xJ ) { ff_poly_remove_root(f,f,m+1,xJ); m--; }
//ff_poly_print (f, m+1);
//printf ("find root returned %d\n", ff_poly_find_root(S,f,m+1));
	return ff_poly_find_root(S,f,m+1);
}

static inline int phi_get_neighbors (ff_t S[], ff_t phi[], int m, ff_t J, ff_t *xJ, int inv)
{
	ff_t f[PHI_MAX_M+2];

	if ( ! phi ) { phi_reduce(m,inv); phi=Phi_ff[inv][m]; }
	if ( m==2 && xJ ) {  			// Hardwire this case for speed
		ff_t t1, t2, D;
		phi2_eval_ff (f, phi, J);
		ff_poly_remove_root_d3 (f,f,xJ);
		_ff_square(t1,f[1]);
		_ff_mult(t2,f[2],f[0]);
		_ff_x2(t2);  _ff_x2(t2);
		_ff_subfrom(t1,t2);
		if ( ! ff_sqrt(&D,&t1) ) return 0;
		if ( ui_randomm(2) ) ff_negate(D);		// randomize choice of sqrt!
		_ff_set(t1,_ff_half);
		_ff_sub(t2,D,f[1]);
		_ff_mult(S[0],_ff_half,t2);
		if ( _ff_zero(D) ) return 1;
		ff_negate(D);
		_ff_sub(t2,D,f[1]);
		_ff_mult(S[1],_ff_half,t2);
		return 2;
	}
	phi_eval_ff (f, phi, m, J);
	if ( xJ ) {  ff_poly_remove_root (f,f,m+1,xJ); m--; }
	return ff_poly_distinct_roots(S,f,m+1);
}

int phi_neighbors (ff_t r[],  int m, ff_t J, ff_t *xJ, int inv)
{
	phi_reduce (m, inv);
	return phi_get_neighbors (r,Phi_ff[inv][m],m,J,xJ,inv);
}

int phi_next_surface_neighbor (ff_t nJ[1], ff_t phi[], int m, int h, ff_t J, ff_t *pJ, int inv)
{
	ff_t S[PHI_MAX_M+2], T[PHI_MAX_M+2], P[PHI_MAX_ENUM_H+2];
	int i,j,k;

//printf("next_surface_neigbor(m=%d,J=-%ld,pJ=-%ld)\n", m, _ff_p-_ff_get_ui(J), (pJ?_ff_p-_ff_get_ui(*pJ):-1));	
	if ( h > PHI_MAX_ENUM_H ) { printf ("Exceeded PHI_MAX_ENUM_H=%d on a %d-volcano with j=%ld mod %ld\n", PHI_MAX_ENUM_H, m, _ff_get_ui(J), _ff_p); abort(); }
	k = phi_get_neighbors (S, phi, m, J, pJ, inv);
//printf ("Surface node J=%ld has %d forward neighbors (*pJ=%ld) w.r.t. m=%d\n", _ff_get_ui(J), k, (pJ?_ff_get_ui(*pJ):-1),m);
	if ( ! k ) return 0;  						 // If there is a double root and pJ is set, then k will be zero (by lemma 2.6 of "Isogeny volcanoes...", a double root is the only root on the surface).  
	if ( k == 1 || (!pJ&&k==2) ) { _ff_set(nJ[0], S[0]); /*printf ("Returned neighbor -%ld\n", _ff_p-_ff_get_ui(nJ[0]));*/ return 1; }				// this means J is static wrt m
	// If we get here then J is not static wrt m
	if ( ! h ) { printf ("phi_next_surface_neighbor: h=0 for m=%d with J=%ld over F_%ld, but k=%d?! (*pJ=%ld)\n", m, _ff_get_ui(J), _ff_p, k, (pJ?_ff_get_ui(*pJ):-1));  abort(); }
	// (at least) one of the k neighbors is on the surface, we need to find one
	for ( i = 0 ; i < k-1 ; i++ ) {
		// check an arbitrary path of length h to see if we hit the floor
		_ff_set(P[0],J);  _ff_set(P[1],S[i]);
		for ( j = 1 ; j <= h ; j++ ) {
			if ( ! phi_get_neighbors(T, phi, m, P[j], P+j-1, inv) ) break;
			_ff_set(P[j+1],T[0]);
		}
		// as a sanity check, make sure we didn't hit the floor too soon
		if ( j < h ) { printf ("phi_next_surface_neighbor: hit the floor after j=%d<h=%d steps, m=%d, J=%ld over F_%ld\n", j, h, m, _ff_get_ui(J), _ff_p); abort(); }
		if ( j > h ) break;		// we didn't hit the floor, so we must have taken a step along the surface
	}
	// if i==k-1 then the first k-1 neighbors are down neighbors, so we assume the last one is on the surface
	_ff_set(nJ[0], S[i]);
//printf ("Returning neighbor -%ld\n", _ff_p-_ff_get_ui(nJ[0]));
	return 1;
}

/*
	Returns a path of length n along the surface of an m-volcano of height h starting from surface node J.
	We assume that the Jacobi symbol (D/m) is 1 (always the case when enumerating via the action of prime forms in the class group).
	From Table 1 of Fouqet&Morain, this means that every surface node has two neighbors on the surface, although these may be the same neighbor in a cycle of length 2.
	*** W must hold space for  n+h nodes ***

	For the sake of speed and simplicity, n is assumed to be less than or equal to the cycle length (the caller should know this).
*/
int phi_surface_path (ff_t W[], int n, ff_t phi[], int m, int h, ff_t J, ff_t *nJ, int inv)
{
	ff_t T[PHI_MAX_ENUM_H+2][PHI_MAX_M+2];
	ff_t f[7];
	int i,j,k,w,x;

dbg_printf ("surface path m=%d, h=%d, J=%ld, nJ=%ld, inv=%d\n", m, h, _ff_get_ui(J), nJ ? _ff_get_ui(*nJ) : -1, inv);
	if ( h > PHI_MAX_ENUM_H ) { printf ("Exceeded PHI_MAX_ENUM_H=%d on a %d-volcano with j=%ld mod %ld\n", PHI_MAX_ENUM_H, m, _ff_get_ui(J), _ff_p); abort(); }
	if ( n==1 ) { _ff_set(W[0],J); return 1; }
	if ( nJ ) {
		k = phi_get_neighbors (T[0]+1, phi, m, J, nJ, inv);
		_ff_set(T[0][0], *nJ); k++;		// insert known neighbor first
	} else {
		k = phi_get_neighbors (T[0], phi, m, J, 0, inv);
	}
dbg_printf("get_neighbors_%d(%ld) returned: ", m, _ff_get_ui(J)); for ( i = 0 ; i < k ; i++ ) dbg_printf ("%ld ", _ff_get_ui(T[0][i])); dbg_printf("\n");
	switch (k) {
	case 0:  printf ("phi_surface_path: No %d-neighbors for J=%ld over F_%ld on %d-volcano of height %d?!\n", m, _ff_get_ui(J), _ff_p, m, h); abort();
	case 1:  if ( n != 2 ) { printf ("phi_surface_path: Request for path of length %d on a cycle of size 2 (m=%d, p=%ld, J=%ld, inv=%d)\n", n, m, _ff_p, _ff_get_ui(J),inv); abort(); }
		       if ( h != 0 ) { printf ("phi_surface_path: h=%d but only %d neighbors found for surface node %ld on %d-volcano over F_%ld\n", h, k, _ff_get_ui(J), m, _ff_p); abort(); }
		       _ff_set(W[0],J);  _ff_set(W[1],T[0][0]);  return 2;
	case 2:
		/*
			If m=2 the only way we can have 2 neighbors is if we have a double root.  This can only happen for |D| <= 16 (by Thm 2.2 of Fouquet-Morain),
			and if it does we must have a 2-cycle.  We do need to handle this case, because it happens for D=-15.
		*/
		if ( m==2 ) {
			if ( n != 2 ) { printf ("phi_surface_path: Request for path of length %d on a cycle of size 2 (m=%d, p=%ld, J=%ld, inv=%d)\n", n, m, _ff_p, _ff_get_ui(J),inv); abort(); }
			_ff_set(W[0],J);
			// the double root is the neighbor on the surface, whom will have exactly one neighbor other than J, the other neighbor of J has either 0 or 2 neighbors that are not J
			k = phi_get_neighbors (T[1], phi, m, T[0][0], &J, inv);
			if ( k==1 ) _ff_set(W[1],T[0][0]); else _ff_set(W[1],T[0][1]);
			return 2;
		}
		// We know that J is static wrt m (so m > 2)
		if ( h ) { printf ("phi_surface_path:h=%d but only %d neighbors found for surface node %ld on %d-volcano of height %d over F_%ld\n", h, k, _ff_get_ui(J), m, h, _ff_p); abort(); }
		
		// note that we aren't using the other root here and we could.
		_ff_set(W[0],J); _ff_set(W[1],T[0][0]);
		// optimize code for small m with separate loops	
		if ( m==3 ) {
			for ( w = 2 ; w < n ; w++ ) {
				phi3_eval_ff (f, phi, W[w-1]);
				ff_poly_remove_root_d4 (f,f,W+w-2);
				if ( ! ff_poly_roots_d3(W+w,f) ) break;
			}
		} else if ( m == 5 ) {
			if ( phi_sparse_factor(phi)==24 ) {
				for ( w = 2 ; w < n ; w++ ) {
					phi5_s24_eval_ff (f, phi, W[w-1]);
					ff_poly_remove_root (f,f,6,W+w-2);
					if ( !  ff_poly_find_root(W+w,f,5) ) break;
				}
			} else {
				for ( w = 2 ; w < n ; w++ ) {
					phi5_eval_ff (f, phi, W[w-1]);
					ff_poly_remove_root (f,f,6,W+w-2);
					if ( ! ff_poly_find_root(W+w,f,5) ) break;
				}
			}
		} else {
			for ( w = 2 ; w < n ; w++ ) {
				k = phi_get_neighbors (W+w, phi, m, W[w-1], W+w-2, inv);
				if ( k != 1 ) { err_printf ("phi_surface_path: Found %d non-previous neighbors on a %d-volcano of height 0 for J=%ld over F_%ld\n", k, m, _ff_get_ui(J), _ff_p); abort(); }
			}
		}
		if ( w != n ) { err_printf ("phi_surface_path failed, only achieved %d of %d steps on %d-volcano of height 0 for J=%ld over F_%ld\n", w, n, m, _ff_get_ui(J), _ff_p); abort(); }
		return w;
	}
	if ( ! h ) { printf ("phi_surface_path: h=0 for m=%d with J=%ld over F_%ld, but k=%d?!\n", m, _ff_get_ui(J), _ff_p, k);  abort(); }
// k need not be m+1 for invariants other than j
//	if ( k != m+1 ) { printf ("phi_surface_path: expected m+1 neighbors, not %d for m=%d with J=%ld over F_%ld\n", k, m, _ff_get_ui(J), _ff_p);  abort(); }
	
	// we could optimize for the case n=3 here
	
//printf ("phi_surface_path: h=%d, m=%d, p=%ld, starting at node %ld\n", h, m, _ff_p, _ff_get_ui(J));
	// If we get here then J is not static wrt m.  Each surface node has m+1 distinct neighbors, 2 on the surface
	_ff_set(W[0],J);  w = 1;
	for ( x = 0 ;; x++ ) {
		_ff_set(W[w],T[(w-1)%h][x]);										// get next neighbor of last known surface node to attempt to extend the path
//printf ("extending from [%d]=%ld\n", w, _ff_get_ui(W[w]));
		if ( x==k-1 && w==n-1 ) return n;									// If we have to test the last neighbor, we know it's on the surface, and if we are done there is no need to extend
		for ( j = w ;; ) {												// walk forward until we hit the floor or finish
			if ( ! phi_get_neighbors(T[j%h], phi, m, W[j], W+j-1, inv) ) break;		// we must get 0 or m neighbors here, but we don't verify this
			_ff_set(W[j+1],T[j%h][0]);  j++;
//printf ("added nbr [%d]=%ld\n", j, _ff_get_ui(W[j]));
			if ( j==w+h ) { w++; x= 0;  k=m; }							// if we have our path by h nodes, we know W[w] is on the surface
			if ( w==n ) return w;
		}
	}
}

// if flag is set, requires a unique common neighor.  otherwise it will allow 2 common neighbors and return both
static inline int phi_common_neighbor (ff_t nbr[1], ff_t J1, ff_t phi1[], int m1, ff_t J2, ff_t phi2[], int m2, int flag)
{
	ff_t f[PHI_MAX_M+2], g[PHI_MAX_M+2], h[PHI_MAX_M+2];
	int d_h;
// COMMENT
dbg_printf ("Finding common neighbor of %ld on %d-volcano and %ld on %d-volcano\n", _ff_get_ui(J1), m1, _ff_get_ui(J2), m2);
	phi_eval_ff (g, phi1, m1, J1);
	phi_eval_ff (f, phi2, m2, J2);
	ff_poly_gcd_small (h, &d_h, f, m2+1, g, m1+1);						// We should always have a unique common neighbor, but for debugging purposes don't assume it
	if ( d_h == 1 ) 	{ ff_poly_roots_d1 (nbr, h); return 1; }
	if ( flag || d_h != 2 ) {
		err_printf ("While finding a common neighbor of %ld on %d-volcano and %ld on %d-volcano (p=%ld):\n", _ff_get_ui(J1), m1, _ff_get_ui(J2), m2, _ff_p);
		err_printf ("gcd has degree %d in phi_common_neighbor, aborting\n", d_h);   abort();
	}
	if ( ! ff_poly_roots_d2 (nbr, h, d_h) ) {
		err_printf ("While finding a common neighbor of %ld on %d-volcano and %ld on %d-volcano (p=%ld):\n", _ff_get_ui(J1), m1, _ff_get_ui(J2), m2, _ff_p);
		err_printf("Degree 2 gcd in phi_common_neighor has no roots, aborting\n"); abort();
	}
	return 2;
}

// Requires a predecessor to J1 which is used to ensure there is a unique common neighbor
static inline void phi_common_neighbor_pred (ff_t nbr[1], ff_t J1, ff_t phi1[], int m1, ff_t J2, ff_t phi2[], int m2, ff_t J0)
{
	ff_t f[PHI_MAX_M+2], g[PHI_MAX_M+2];
	ff_t h[2];
//dbg_printf ("Finding common neighbor of %ld on %d-volcano and %ld on %d-volcano, using pred %ld\n", _ff_get_ui(J1), m1, _ff_get_ui(J2), m2, _ff_get_ui(J0));
	phi_eval_ff (g, phi1, m1, J1);
	ff_poly_remove_root (g, g, m1+1, &J0);
	phi_eval_ff (f, phi2, m2, J2);
	ff_poly_clf (h, f, m2+1, g, m1);
	ff_poly_roots_d1 (nbr, h);
//dbg_printf ("Found neighbor %ld\n", _ff_get_ui(nbr[0]));
}

// Used in uncertain situations, verifies result and returns FALSE if no root of phi1 is found.  J0 is predecessor of J1 on m1-volcano, required
static inline int phi_common_neighbor_verify (ff_t nbr[1], ff_t J1, ff_t phi1[], int m1, ff_t J2, ff_t phi2[], int m2, ff_t J0)
{
	ff_t f[PHI_MAX_M+2], g[PHI_MAX_M+2], h[PHI_MAX_M+2];
	int d_h;

//dbg_printf ("Attempting to find a common neighbor of %ld on %d-volcano and %ld on %d-volcano\n", _ff_get_ui(J1), m1, _ff_get_ui(J2), m2);
	phi_eval_ff (g, phi1, m1, J1);
	ff_poly_remove_root (g, g, m1+1, &J0);
	phi_eval_ff (f, phi2, m2, J2);
	ff_poly_gcd_small (h, &d_h, f, m2+1, g, m1);					// use poly_gcd, not poly_clf, since they might not have a linear factor in common (even with root removed)
	if ( d_h < 1 ) return 0;
	if ( d_h > 1 ) {
		printf ("While attempting to find a common neighbor of %ld on %d-volcano and %ld on %d-volcano (p=%ld):\n", _ff_get_ui(J1), m1, _ff_get_ui(J2), m2, _ff_p);
		printf ("gcd had degree %d in phi_common_neighbor_verify, aborting \n", d_h);  abort();
	}
	ff_poly_roots_d1 (nbr, h);
//dbg_printf ("Found neighbor %ld\n", _ff_get_ui(nbr[0]));
	return 1;
}

// finds a common m1-neighbor of J1 and m2-neighbor of J2, given J0 an m2-neighbor of J1 and an m1-neighbor of J2
static inline void phi_common_neighbor_corner (ff_t nbr[1], ff_t J1, ff_t phi1[], int m1, int h1, ff_t J2, ff_t phi2[], int m2, int h2, ff_t J0)
{
	ff_t f[PHI_MAX_M+2];
	ff_t nbrs[2], nJ1, nJ2, r;
	
//dbg_printf ("Finding common neighbor of %ld on %d-volcano and %ld on %d-volcano, using corner %ld\n", _ff_get_ui(J1), m1, _ff_get_ui(J2), m2, _ff_get_ui(J0));
	if ( phi_common_neighbor (nbrs, J1, phi1, m1, J2, phi2, m2,0) == 2 ) {
		if ( ! _ff_equal(nbrs[0],nbrs[1]) ) {
//dbg_printf ("Two possible corners %ld and %ld\n", _ff_get_ui(nbrs[0]), _ff_get_ui(nbrs[1]));
			// note we assume phi is non-null so that inv=0 parameter will be ignored
			if ( ! phi_next_surface_neighbor (&nJ2, phi1, m1, h1, J2, &J0, 0) ) { err_printf ("Unable to resolve ambiguity in phi_common_neighbor_corner, m1=%d has order 2\n", m1); abort(); }
//dbg_printf ("Next %d-neighbor of %ld is %ld, using pred %ld\n", m1, _ff_get_ui(J2), _ff_get_ui(nJ2), _ff_get_ui(J0));
			if ( ! phi_next_surface_neighbor (&nJ1, phi1, m1, h1, nbrs[0], &J1, 0) ) { err_printf ("Unable to resolve ambiguity in phi_common_neighbor_corner, m1=%d has order 2\n", m1); abort(); }
//dbg_printf ("Next %d-neighbor of %ld is %ld, using pred %ld\n", m1, _ff_get_ui(nbrs[0]), _ff_get_ui(nJ1), _ff_get_ui(J1));
			phi_eval_ff (f, phi2, m2, nJ2);
			ff_poly_eval (&r, f, m2+1, &nJ1);
			if ( ! _ff_zero(r) ) {
				_ff_set(nbrs[0],nbrs[1]);	// nbr[1] must be the right choice
			} else {
				if ( ! phi_next_surface_neighbor (&nJ1, phi1, m1, h1, nbrs[1], &J1, 0) ) { err_printf ("Unable to resolve ambiguity in phi_common_neighbor_corner, m1=%d has order 2\n", m1); abort(); }
//dbg_printf ("Next %d-neighbor of %ld is %ld, using pred %ld\n", m1, _ff_get_ui(nbrs[1]), _ff_get_ui(nJ1), _ff_get_ui(J1));	
				phi_eval_ff (f, phi2, m2, nJ2);
				ff_poly_eval (&r, f, m2+1, &nJ1);
				if ( _ff_zero(r) ) { err_printf ("2 distinct choices in phi_common_neighbor_corner both extended successfully, m1=%d, m2=%d, don't know what to do!\n", m1, m2); abort(); }
			}
		}
	}
	_ff_set(nbr[0],nbrs[0]);
}


int phi_poly_verify_edge (ff_t J1, ff_t J2, int m, int inv)
{
	ff_t f[PHI_MAX_M+2];
	ff_t r;
	
	phi_reduce (m, inv);
	phi_eval_ff (f, Phi_ff[inv][m], m, J2);
	ff_poly_eval(&r,f,m+1,&J1);
	return ( _ff_zero(r) ? 1 : 0 );
}


int phi_poly_verify_2path (ff_t J1, ff_t J2, int m1, int m2, int inv)
{
	ff_t f[PHI_MAX_M+2], g[PHI_MAX_M+2], h[PHI_MAX_M+2];
	int d_h;
	phi_reduce (m1, inv);  phi_reduce (m2, inv);
	phi_eval_ff (f, Phi_ff[inv][m1], m1, J1);  phi_eval_ff(g, Phi_ff[inv][m2], m2, J2);
	ff_poly_gcd_small (h, &d_h, f, m1+1, g, m2+1);
	if ( d_h == 1 ) return 1;
	if ( ! d_h ) return 0;
	ff_poly_monic(h,&d_h,h,d_h);
	return ff_poly_count_distinct_roots(h,d_h);
}


/*
	Given a path V of length n on an m1-volcano, and W[0] m2-isogenous to V[0],
	Extends the path W to length n on an m1-volcano, with W[i] m2-isogenous to V[i]

	Uses gcds unless m2 is too large to make it helpful.  Always uses GCD to get V[1] to ensure consistent orientation
*/
int phi_surface_parallel_path (ff_t W[], ff_t V[], int n, ff_t phi1[], int m1, int h1, ff_t phi2[], int m2, int inv, int cycle)
{
	ff_t h[2], W2;

//dbg_printf ("phi_surface_parallel_path(%d,%d) inv=%d, cycle=%d\n", m1, m2, inv, cycle);
//dbg_printf ("vertex W[0]=%ld on %d-volcano is %d-isogenous to V[0]=%ld, and V[1] = %ld\n", _ff_get_ui(W[0]), m1, m2, _ff_get_ui(V[0]), _ff_get_ui(V[1]));
	if ( phi_common_neighbor (h, W[0], phi1, m1, V[1], phi2, m2,0) == 2 ) {
		if ( ! _ff_equal(h[0],h[1]) ) {
			// we have 2 distinct choices, we need to figure out which one is correct
			if ( n > 2 ) {
				if ( ! phi_common_neighbor_verify (&W2, h[0], phi1, m1, V[2], phi2, m2, W[0]) ) {
					_ff_set(h[0],h[1]);	// h[1] must be the right choice
				} else if ( phi_common_neighbor_verify (&W2, h[1], phi1, m1, V[2], phi2, m2, W[0]) ) {
					err_printf ("2 distinct choices in phi_surface_parallel_path both extended successfully, m1=%d, m2=%d, don't know what to do!\n", m1, m2); abort();
				}
			} else {
				err_printf ("2 distinct choices in phi_surface_parallel_path with n=2, m1=%d, m2=%d, don't know what to do!\n", m1, m2); abort();
			}
		}
	}
	_ff_set(W[1],h[0]);
	if ( n==2 ) return 2;

//dbg_printf ("W[1]=%ld\n", _ff_get_ui(W[1]));
	if ( poly_root_step_cost (m1, h1, _ff_p, inv) < poly_gcd_step_cost(m1,m2,_ff_p,inv) ) return phi_surface_path (W, n, phi1, m1, h1, W[0], W+1, inv);
	if ( cycle ) phi_surface_gcd_cycle (W, V, m1, n, m2, 2, phi1, phi2);
	else phi_surface_gcd_path (W, V, m1, n, m2, 2, phi1, phi2);
	return n;
	
/*	for ( i = 2 ; i < n ; i++ ) {
		phi_eval_ff (g, phi1, m1, W[i-1]);
		ff_poly_remove_root (g, g, m1+1, W+i-2);	
		phi_eval_ff (f, phi2, m2, V[i]);
		ff_poly_clf (h, f, m2+1, g, m1);
		ff_poly_roots_d1 (W+i, h);
	}
*/
}

int old_phi_enum_roots (ff_t roots[], ff_t J0, long p[], long n[], int r, int v, int inv)
{
	long e[IQ_MAX_GENS];
	ff_t pJ[IQ_MAX_GENS], J[IQ_MAX_GENS], nJ;
	int h[IQ_MAX_GENS];
	struct phi_vshape_struct s;
	register int i,t;

	if ( ! r ) { _ff_set(roots[0],J0); return 1; }	// class group is trivial (h(D)=1)
	phi_vshape (&s, v);
	
	// set h[i] to the height of p[i] in v
	for ( i = 0 ; i < r ; i++ ) { h[i] = 0; for ( t = 0 ; t < s.k ; t++ ) if ( s.p[t] == p[i] ) h[i] = s.h[t]; }
	for ( i = 0 ; i < r ; i++ ) { _ff_set(J[i],J0);  e[i] = 0; phi_reduce (p[i], inv); }
	t = 0; //_ff_set(roots[t++],J0);
	for (;;) {
		t += phi_surface_path (roots+t, n[0], Phi_ff[inv][p[0]], p[0], h[0], J[0], 0, inv);  e[0] = n[0]-1;
		for ( i = 0 ; i < r && e[i]==n[i]-1 ; i++ );
		if ( i == r ) break;
		if ( ! phi_next_surface_neighbor (&nJ, Phi_ff[inv][p[i]], p[i], h[i], J[i], (e[i]?pJ+i:0),inv) ) { printf ("phi_next_surface_neighbor failed for J=%ld with m=%ld, v=%d\n", _ff_get_ui(J[i]), p[i], v); abort(); }
		/*_ff_set(roots[t++],nJ);*/  _ff_set(pJ[i],J[i]);  _ff_set(J[i],nJ);  e[i]++;
		while ( i-- > 0 ) { e[i] = 0; _ff_set(J[i],J[i+1]); }
	}
	return t;
}

// *** IMPORTANT *** roots must have room for h+PHI_MAX_H+1 entries!!! 
int phi_enum_roots (ff_t roots[], ff_t J0, long p[], long n[], long o[], long aux_p[], long aux_map[IQ_MAX_GENS][IQ_MAX_GENS], int r, int v, int inv)
{
	long e[IQ_MAX_GENS], M[IQ_MAX_GENS];
	int h[IQ_MAX_GENS], off[IQ_MAX_GENS], poff[IQ_MAX_GENS];
	int use_aux_p0[IQ_MAX_GENS];
	int aux_p0, aux_off00, aux_p1, aux_off10, aux_off10_verified;
	struct phi_vshape_struct s;
	register int i,j,t;

	if ( ! r ) { _ff_set(roots[0],J0); return 1; }	// class group is trivial (h(D)=1)
	phi_vshape (&s, v);

	// set h[i] to the height of p[i] in v
	for ( i = 0 ; i < r ; i++ ) {
		if ( p[i] > PHI_MAX_M ) { err_printf ("prime ell=%ld exceeds PHI_MAX_M in phi_enum_roots\n", p[i]); abort(); }
		h[i] = 0; for ( t = 0 ; t < s.k ; t++ ) if ( s.p[t] == p[i] ) h[i] = s.h[t];
	}
	for ( i = 0 ; i < r ; i++ ) { phi_reduce (p[i], inv); e[i] = 0; off[i] = 0; }
	for ( M[0]=1, i = 1 ; i < r ; i++ ) M[i] = M[i-1]*n[i-1];
	
//	if ( aux_p && aux_p[0] && p[0] <= PHI_SMALL_M && phi_maxgcd[p[0]][h[0]] && aux_p[0] > phi_maxgcd[p[0]][h[0]] ) aux_p[0] = 0;			// check gcd cost in lookup table to make sure it is worth using

	// Note that we don't bother maintaining e[0], it implicitly runs from 0 to n[0]-1 each time a path is enumerated

	// For simplicity, we don't optimize the case where aux_p[0] is equivalent to p[0] or p[0]^2.   Even though these are essentially "best cases", they will never occur in large computations (i.e. when p[0]^3 < |D/3|)
	if ( aux_p && aux_p[0] && aux_map[0][0] > 2 ) {
		aux_p0 = aux_p[0];  aux_off00 = aux_map[0][0];
	} else {
		aux_p0 = 0;  aux_off00 = 0;
	}
	
	if ( aux_p0 ) {
// COMMENT
dbg_printf ("auxilliary prime for ell0=%ld is %ld with exponent %d, n0=%ld\n", p[0], aux_p[0], aux_off00, n[0]);
		phi_reduce (aux_p0, inv);
		phi_surface_path (roots, aux_off00, Phi_ff[inv][p[0]], p[0], h[0], J0, 0, inv);
		phi_surface_gcd_cycle (roots, roots, p[0], n[0], aux_p0, aux_off00, Phi_ff[inv][p[0]], Phi_ff[inv][aux_p0]);
	} else {
dbg_printf ("Enumerating surface path of length %ld for ell0=%ld\n", n[0], p[0]);
		phi_surface_path (roots, n[0], Phi_ff[inv][p[0]], p[0], h[0], J0, 0, inv);
	}
	t = n[0];
dbg_printf ("initial %ld-path: ", p[0]); for ( i = 0 ; i < t ; i++ ) dbg_printf ("%ld ", _ff_get_ui(roots[i])); dbg_printf("\n");
	if ( r == 1 ) return t;
	
	// Don't bother with the auxilliary prime for p0 if parallel path with p_i is better
	if ( aux_p0 ) for ( i = 1 ; i < r ; i++ ) use_aux_p0[i] = ( poly_gcd_step_cost (p[0], aux_p0, _ff_p, inv) < poly_gcd_step_cost (p[0], p[i], _ff_p, inv) ? 1 : 0 );
	else for ( i = 1 ; i < r ; i++ ) use_aux_p0[i] = 0;

	// Currently we only use auxilliar primes for ell0 and ell1.  We may want to generalize this in the future
	aux_off10_verified = 0;
	if ( aux_p && aux_p[1] ) {
dbg_printf ("auxilliary prime for ell1=%ld is %ld with exponent %ld\n", p[1], aux_p[1], aux_map[1][1]);
		aux_p1 = aux_p[1]; phi_reduce (aux_p1, inv); aux_off10 = aux_map[1][0];
	} else {
		aux_p1 = 0;  aux_off10 = 0;
	}

	for ( i = 1 ; i < r && e[i]==n[i]-1 ; i++ );
	while ( i < r ) {
		// note i is always > 0
		for ( j = i+1 ; j < r && ! e[j] ; j++ );
		if ( j < r ) {
			if ( e[i] ) {
				phi_common_neighbor_pred (roots+t, roots[off[i]], Phi_ff[inv][p[i]], p[i], roots[t-M[j]], Phi_ff[inv][p[j]], p[j], roots[poff[i]]);
			} else {
				// This call will abort if it doesn't find a unique common neighbor, we need to be sure things are oriented consistently
				//phi_common_neighbor (roots+t, roots[off[i]], Phi_ff[inv][p[i]], p[i], roots[t-M[j]], Phi_ff[inv][p[j]], p[j], 1);				
				phi_common_neighbor_corner (roots+t, roots[off[i]], Phi_ff[inv][p[i]], p[i], h[i], roots[t-M[j]], Phi_ff[inv][p[j]], p[j], h[j], roots[poff[j]]);				
			}
dbg_printf ("next %ld-nbr: [%d]=%ld (via %ld-nbr [%d]=%ld and %ld-nbr [%ld]=%ld)\n", p[i], t, _ff_get_ui(roots[t]), p[i], off[i], _ff_get_ui(roots[off[i]]), p[j], t-M[j], _ff_get_ui(roots[t-M[j]]));
		} else {
			if ( aux_p1 && i==1 && e[1] && e[1]+1 >= aux_map[1][1] ) {
dbg_printf ("Using auxiliary prime ell=%d to obtain %ld-neighbor at position %d (e[1]=%ld, aux_map[1][1]=%ld)\n", aux_p1, p[i], t, e[1], aux_map[1][1]);
				j = (e[1]-aux_map[1][1]+1)*M[1] + aux_off10;
dbg_printf ("Using reference point [%d]=%ld\n", j, _ff_get_ui(roots[j]));
				if ( ! aux_off10_verified ) {
					if ( ! phi_common_neighbor_verify (roots+t, roots[off[i]], Phi_ff[inv][p[1]], p[1], roots[j], Phi_ff[inv][aux_p1], aux_p1, roots[poff[1]]) ) {
						if ( n[0] != o[0] ) goto cant_flip;		// we can't flip the offset if we don't have the whole cycle (this can occur when N=ell0 is ramified)
dbg_printf ("Flipping offset for auxilliary prime ell=%d\n", aux_p1);
						j -= aux_off10;  aux_off10 = n[0]-aux_off10; j+= aux_off10;
						phi_common_neighbor_pred (roots+t, roots[off[i]], Phi_ff[inv][p[1]], p[1], roots[j], Phi_ff[inv][aux_p1], aux_p1, roots[poff[1]]); 		// This should always work.
					}
					aux_off10_verified = 1;
				} else {
					phi_common_neighbor_pred (roots+t, roots[off[i]], Phi_ff[inv][p[1]], p[1], roots[j], Phi_ff[inv][aux_p1], aux_p1, roots[poff[1]]);
				}
			} else {
cant_flip:
				if ( ! phi_next_surface_neighbor (roots+t, Phi_ff[inv][p[i]], p[i], h[i], roots[off[i]], (e[i]?roots+poff[i]:0),inv) ) { printf ("phi_next_surface_neighbor failed for J=%ld with m=%ld, v=%d\n", _ff_get_ui(roots[off[i]]), p[i], v); abort(); }
			}
dbg_printf ("next %ld-nbr: [%d]=%ld (via %ld-nbr [%d]=%ld)\n", p[i], t, _ff_get_ui(roots[t]), p[i], off[i], _ff_get_ui(roots[off[i]]));
		}
		poff[i]=off[i];  off[i] = t;  e[i]++;
		for ( j = i-1 ; j ; j-- ) { e[j] = 0; off[j] = off[j+1]; }
		if ( use_aux_p0[i]  ) {
dbg_printf("Calling parallel path with aux_p0=%d, with aux_off00=%d\n", aux_p0, aux_off00);
			phi_surface_parallel_path (roots+t, roots+poff[i], aux_off00, Phi_ff[inv][p[0]], p[0], h[0], Phi_ff[inv][p[i]], p[i], inv, 0);
			phi_surface_gcd_cycle (roots+t, roots+t, p[0], n[0], aux_p0, aux_off00, Phi_ff[inv][p[0]], Phi_ff[inv][aux_p0]);			
		} else {
dbg_printf("Calling parallel path with ell0=%ld n0=%ld ell[i]=%ld, not using aux_p0\n", p[0], n[0], p[i]);
//			if ( inv || aux_p ) { 	// temporary hack for testing
				phi_surface_parallel_path (roots+t, roots+poff[i], n[0], Phi_ff[inv][p[0]], p[0], h[0], Phi_ff[inv][p[i]], p[i], inv, (n[0]==o[0]?1:0));		// check whether we have a cycle or not
//			} else {	// temporary hack for testing
//				phi_surface_path (roots+t, n[0], Phi_ff[inv][p[0]], p[0], h[0], roots[t], 0, inv);
//			}
		}
		t += n[0];
//dbg_printf ("next %d-path: ", p[0]); for ( i = t-n[0] ; i < t ; i++ ) dbg_printf ("%ld ", _ff_get_ui(roots[i])); dbg_printf("\n");
		for ( i = 1 ; i < r && e[i]==n[i]-1 ; i++ );
	}
	return t;
}

int phi_get_graph (ff_t E1[], ff_t E2[], int n, ff_t phi[], int m, ff_t J)
{
	ff_t cJ;
	int i, j, k, v;
	
	_ff_set(cJ,J);
	k = v = 0;
	for(;;) {
		if ( k+m+1 > n ) { printf ("Too many edges in phi_get_graph, exceeded %d\n", n); abort(); }
		j = phi_get_neighbors (E2+k, phi, m, cJ, 0,0);
		if ( ! j ) { if ( k ) printf ("Error isolated node in non-trivial graph!\n"); abort(); return 0; }
		for ( i = 0 ; i < j ; i++ ) _ff_set(E1[k+i], cJ);
		k += j;
		do {
			for ( i = 0 ; i < k ; i++ ) if ( _ff_equal(E1[i],E2[v]) ) break;
			if ( i == k ) break;
			v++;
		} while ( v < k );
		if ( v == k ) break;
		_ff_set(cJ,E2[v]);
	}
	return k;
}

// Algorithm DESCEND from Hilbert CRT paper, p. 12
void phi_descend (ff_t *J1, ff_t J0, int k, int h, int m, int inv)
{
	ff_t S[PHI_MAX_M+1];
	ff_t P[PHI_MAX_H+1];
	register int i, j, pi, ni;

dbg_printf ("Descending from J0=%ld at level %d on %d-volcano of height %d, inv=%d\n", _ff_get_ui(J0), k, m, h, inv);
	if ( k >= h ) { printf ("Error, attempt to descend from level %d on a %d-volcano of height %d\n", k, m, h); abort(); }
	phi_reduce (m,inv);										// we may be called externally, so make sure we have phi loaded and reduced
	if ( k == 0 ) {
		_ff_set(P[0],J0);
		i = 1;
		if ( ! phi_get_neighbor(P+i, Phi_ff[inv][m], m, J0, 0, inv) ) { printf ("Error, unable to descend from non-floor node %ld on %d-volcano of height %d (p=%ld)\n", _ff_get_ui(J0), m, h, _ff_p); abort(); }
// COMMENT
dbg_printf ("First  vertex in path %ld\n", _ff_get_ui(P[i]));
		for(j=1;;j++) {
			ni=(i<h?i+1:0);  pi = (i?i-1:h);
dbg_printf ("Current vertex in path %ld, Previous vertex in path %ld\n", _ff_get_ui(P[i]), _ff_get_ui(P[pi]));
			if ( ! phi_get_neighbor(P+ni, Phi_ff[inv][m], m, P[i], P+pi, inv) ) { dbg_printf ("Hit the floor."); break; }
dbg_printf ("Next vertex in path %ld\n", _ff_get_ui(P[ni]));
			i = ni;
			if ( j > 100*h ) { printf ("Error, unable to find floor after %d steps, claimed height value h=%d of %d-volcano containing j=%ld is suspect (inv=%d, p=%ld)\n", j, h, m, _ff_get_ui(J0), inv, _ff_p); abort(); }
		}
		if ( j < h ) { printf ("Error, hit the floor after %d steps from the surface on a %d-volcano of height %d\n", j, m, h); abort(); }
		_ff_set(*J1,P[(j-h+1)%(h+1)]);
		return;
	}
	i = phi_get_neighbors (S, Phi_ff[inv][m], m, J0, 0, inv);
	if ( i != m+1 ) { printf ("Error node at level %d on %d-volcano of height %d has degree %d\n", k, m, h, i); exit(0); }
	_ff_set(P[0],J0); _ff_set(P[1],S[0]);
	for(i=1;i < h-k;i++)
		if ( ! phi_get_neighbor(P+i+1, Phi_ff[inv][m], m, P[i], P+i-1, inv) ) { printf ("Error, hit the floor after %d steps from level %d on a %d-volcano of height %d\n", i, k, m, h); abort(); }
	if ( ! phi_get_neighbor(P, Phi_ff[inv][m], m, P[i], P+i-1, inv) ) { _ff_set(*J1,S[0]); } else { _ff_set(*J1,S[1]); }
	return;
}

// Algorithm ASCEND from Hilbert CRT paper p. 13
void phi_ascend (ff_t *J1, ff_t J0, int k, int h, int m, int inv)
{
	ff_t S[PHI_MAX_M+1];
	ff_t pJ,J,nJ;
	register int i, j, n;

dbg_printf ("Ascending from J0=%ld at level %d on %d-volcano of height %d, inv=%d\n", _ff_get_ui(J0), k, m, h, inv);
	if ( k <= 0 ) { printf ("Error, attempt to ascend from level %d on a %d-volcano of height %d\n", k, m, h); abort(); }
	phi_reduce (m,inv);										// we may be called externally, so make sure we have phi loaded and reduced
	n = phi_get_neighbors (S, Phi_ff[inv][m], m, J0, 0, inv);
	if ( n == 1 ) {
		if ( k != h ) { printf ("Error, only one neighber of %ld at level %d > 0 on a height %d on %d-volcano\n", _ff_get_ui(J0), k, h, m); exit(0); }
		_ff_set(*J1,S[0]);
		return;
	}
	if ( n != m+1 || k >= h) { printf ("Error, found %d neighbors of vertex %ld at level %d of a %d-volcano of height %d\n", n, _ff_get_ui(J0), k, m, h); exit(0); }
	// For the moment, don't optimize by assuming the last neighbor must be the parent after eliminating the other m, explicitly check
	for ( i = 0 ; i <= m ; i++ ) {
		_ff_set(pJ,J0); _ff_set(J,S[i]);
		for( j=1 ; j < h-k ; j++ ) {
			if ( ! phi_get_neighbor(&nJ, Phi_ff[inv][m], m, J, &pJ, inv) ) { printf ("Error, unable to walk a path of length %d from %ld at level %d on %d-volcano of height %d\n", h-k, _ff_get_ui(J0), k, m, h); abort(); }
			_ff_set(pJ,J);  _ff_set(J,nJ);
		}
		if ( phi_get_neighbor(&nJ, Phi_ff[inv][m], m, J, &pJ, inv) ) break;
	}
	if ( i > m ) { printf ("Error, unable to find parent node from vertex %ld at level %d of a %d-volcano of height %d\n", _ff_get_ui(J0), k, m, h); abort(); }
	 _ff_set(*J1,S[i]);
	return;
}

// input must be a j-invariant
int phi_find_surface (ff_t *outJ, ff_t inJ, unsigned long up[], unsigned long uh[], int uw, unsigned long vp[], unsigned long vh[], int vw)
{
	register int i, j;
	register ff_t t0;
	ff_t J;

	_ff_set(J,inJ);
	for ( i = j = 0 ; i < uw || j < vw ; ) {
		if ( i < uw && (j==vw || up[i] < vp[j]) ) {
			phi_adjust_level (&J, J, uh[i], uh[i], up[i], 0);  i++;
		} else if ( j < vw && (i==uw || vp[j] < up[i]) ) {
			phi_adjust_level (&J, J, 0, (int)vh[j], (int)vp[j], 0);  j++;
		} else {
			if ( up[i] != vp[j] || i >= uw || j >= vw ) { err_printf ("Bug, factors of u and v don't match up!\n"); abort(); }
			phi_adjust_level (&J, J, uh[i], uh[i]+(int)vh[j], up[i], 0);  i++; j++;
		}
	}
	// make sure we didn't bump into 0 or 1728
	_ff_set_ui(t0,1728);
	if ( _ff_zero(J) || _ff_equal(J,t0) ) return 0;
	_ff_set(*outJ,J);
	return 1;
}


// attempts to find the floor via a path begining J0->J1 of length <= h, for h > 0, where J0 is not on the floor
// returns the length of the path or h+1 if unsuccessful
// S contains the neighbors of J1, not including J0
// input must be a j-invariant
int phi_test_height(ff_t J0, ff_t J1, int h, ff_t phi[], int m, int inv)
{
	ff_t S[PHI_MAX_M+1];
	ff_t J, pJ;
	register int k, len;
	
	k = phi_get_neighbors (S, phi, m, J1, &J0, inv);
	if ( ! k ) return 1;
	_ff_set(J,S[0]);  _ff_set(pJ,J1);
	for ( len = 2 ; len <= h ; len++ ) {
		k = phi_get_neighbors (S, phi, m, J, &pJ, inv);
		if ( ! k ) break;
		_ff_set(pJ,J);  _ff_set(J,S[0]);
	}
	return len;
}

// Algorithm FINDLEVEL from Hilbert CRT paper p. 11.  Follows convention in the paper
// that the surface is at level 0 and the floor is at level h (unlike phi_find_surface_pp) below.
// Note that phi_get_neighbors returns a count of *distinct* neighbors, so we modify the algorithm somewhat
int phi_find_level(ff_t J, int h, int m, int inv)
{
	ff_t S[PHI_MAX_M+1];
	register int i, j, k;
	
	if ( h == 0 ) return 0;
	phi_reduce (m,inv);										// we may be called externally, so make sure we have phi loaded and reduced
	k = phi_get_neighbors (S, Phi_ff[inv][m], m, J, 0, inv);
//dbg_printf ("Finding level of vertex %ld with %d neighbors on %d-volcano of height %d\n", _ff_get_ui(J), k, m, h);
	if ( k == 1 ) return h;									// we must be on the floor
	if ( k < m+1 ) return 0;									// we must be on the surface to get a double root
	// we are at a level k>h and have m+1 neighbors
	// probe two neighbors in search of a descending path
	i = phi_test_height (J,S[0],h,Phi_ff[inv][m],m,inv);
	i = _min(i,h);
	j = phi_test_height (J,S[1],i,Phi_ff[inv][m],m,inv);
	return h-_min(i,j);
}

// Step 2b of Algorithm 1.2 from Hilbert CRT paper p. 12.  Adjusts J0 to J1 at level k on an m-volcano of height h
void phi_adjust_level (ff_t *J1, ff_t J0, int k, int h, int m, int inv)
{
	register int k0;

	k0 = phi_find_level(J0, h, m, inv);
// COMMENT
dbg_printf ("adjusting level of vertex %ld from %d to %d on %d-volcano of height %d\n", _ff_get_ui(J0), k0, k, m, h);
	while ( k0 > k ) { phi_ascend (J1,J0,k0,h,m,inv);  k0--; _ff_set(J0,*J1); }		// note ascending reduces the level (surface is level 0)
	while ( k0 < k ) { phi_descend (J1,J0,k0,h,m,inv);  k0++; _ff_set(J0,*J1); }		// and descending increases the level (floor is level h)
dbg_printf ("adjusted vertex %ld on level %d on %d-volcano of height %d\n", _ff_get_ui(J0), k, m, h);
}

// *** The algorithm below treats the floor as level 0, not level h.  It predates the notatation in the Hilbert CRT paper and is not currently used ***

/*
	Finds the j-invariant of a curve on the surface of the volcano of m-isogenies containing J,
	where J is the j-invariant of a curve with trace t, given that m^h is a maximal prime power
	dividing v, where, 4p = t^2-v^2D and D is the (negative) fundamental discriminant of
	the endomorphism ring of a curve at the surface.

	If J is already at the surface, we always return *sJ=J (this is relied upon above)

	Note that it is possible to get a double root in Phi_m(X,J), but only when J is at the surface,
	and only when D <= 4m^2 (we don't rely on this latter fact here).
	(see Lemma 2.6 of Fouquet & Morain "Isogeny Volcanoes and the SEA algorithm")

	Note that this algorithm assumes j-invariants are being used, but this is always the case in the Hilbert CRT algorithm
	since find_curve returns a j-invariant.
*/
int phi_find_surface_pp (ff_t *sJ,  ff_t phi[], int m, int h, ff_t J)
{
	ff_t S[PHI_MAX_M+1], S1[PHI_MAX_M+1], S2[PHI_MAX_M+1], J1, J2, J3;
	int i,j,k,up,level;

	k = phi_get_neighbors (S, phi, m, J, 0, 0);
	if ( k == m ) { _ff_set(*sJ,J); return 1; }						// must be a double root
	if ( h == 1 ) {
		if ( k == 1 ) _ff_set(*sJ,S[0]); else _ff_set(*sJ,J);
	} else if ( h == 2 ) {
		if ( k == 1 ) {
			// J is at level 0 and its single neighbor J1 is at level 1
			_ff_set(J1,S[0]);
			k = phi_get_neighbors (S, phi, m, J1, &J, 0);
			for ( i = 0 ; i < k ; i++ ) {
				_ff_set(J2, S[i]);
				j = phi_get_neighbors (S, phi , m, J2, 0, 0);
				if ( j == m+1 ) break;
				if ( j == m ) { _ff_set(*sJ,J2); return 1; }		// must be a double root
				if ( j != 1 ) { printf ("Error1, k=%d %d-isogenies found for j-invariant %ld over F_%ld, expected 1\n", k, m, _ff_get_ui(J2), _ff_p); abort(); }
			}
			if ( i < k ) _ff_set(*sJ,J2); else _ff_set(*sJ,S[k]);		// if first k failed, last neighbor must be at the surface
		} else {
			// find first of m+1 neighbors in S which is not at level 0
			for ( i = 0 ; i < m ; i++ ) {
				k = phi_get_neighbors (S2, phi , m, S[i], 0, 0);
				if ( k == m+1 ) break;
				if ( k == m ) { _ff_set(*sJ,S[i]); return 1; }		// must be a double root
				if ( k != 1 ) { printf ("Error2, k=%d %d-isogenies found for j-invariant %ld over F_%ld, expected 1\n", k, m, _ff_get_ui(S[i]), _ff_p); abort(); }
			}
			if ( i > 0 ) { _ff_set(*sJ,S[i]); return 1; }				// if one of our neighbors was on the floor, we are at level 1 and found an up edge or know it must be the last edge
			k = phi_get_neighbors (S2, phi, m, S[1], 0, 0);		// if first neighbor was not on the floor, check second neighbor also
			if ( k == 1 ) { _ff_set(*sJ,S[0]); return 1; }			// we must be at level 1 and first neighbor is at the surface
			_ff_set(*sJ,J);									// we have two neighbors who aren't on the floor, so we must be  at the surface
		}
	} else if ( m==2 && h==3 ) {							// specializing to 2^3 allows us to reduce the cost by nearly a factor of 2 over the general case (~5 vs 9 nbr checks)
		if ( k == 1 ) {
			_ff_set(J1,S[0]);
			phi_get_neighbors (S, phi, m, J1, &J, 0);
			// S now holds neighbors of node J1 at level 1, excluding J (which is at level 0)
			_ff_set(J2, S[0]);
			k = phi_get_neighbors (S2, phi, m, J2, &J1, 0);
			if ( !k ) {
				_ff_set(J2, S[1]);
				k = phi_get_neighbors (S2, phi, m, J2, &J1, 0);
			}
			if ( k != m ) { printf ("Error3, k=%d %d-isogenies found for j-invariant %ld over F_%ld, expected  %d\n", k+1, m, _ff_get_ui(S[1]), _ff_p, m+1); abort(); }
			// S2 now holds (exactly m=2) neighbors of node J2 at level 2, excluding J1 (at level 1), one of which is at the surface
			_ff_set(J3,S2[0]);								// J3 is at level 1 or 3
			k = phi_get_neighbors (S, phi, m, J3, &J2, 0);
			if ( k == m-1 ) { _ff_set(*sJ,J3); return 1; }			// must be a double root
			k = phi_get_neighbors (S, phi, m, S[0], &J3, 0);
			if ( ! k ) _ff_set(*sJ,S2[1]); else _ff_set(*sJ,J3);		// if a neigbhor of J3 is at the floor (no neighbors but J3), then J3 is at level 1 and the other neighbor of J2 is at level 3
														// otherwise, J3 is at level 3.
		} else {
			// we know that J is at level > 0
			_ff_set(J1,S[0]);
			k = phi_get_neighbors (S1, phi, m, J1, &J, 0);
			if ( k == m-1 ) { _ff_set(*sJ,J1); return 1; }			// must be a double root
			if ( ! k ) {
				// we know that J is at level 1 and J1 is at level 0
				_ff_set(J2,S[1]);
				k = phi_get_neighbors (S2, phi, m, J2, &J,0);
				if ( ! k ) {
					_ff_set(J2,S[2]);
					k = phi_get_neighbors (S2, phi, m, J2, &J,0);
					if ( ! k ) { printf ("Error4, k=%d %d-isogenies found for j-invariant %ld over F_%ld, expected >= %d\n", k+1, m, _ff_get_ui(S[1]), _ff_p, m); abort(); }
				}
				// we know that J2 is at level 2 and that S2 holds two of its neighbors (excluding J), one of which is at level  3 (the other is at level 1)
				k = phi_get_neighbors (S, phi, m, S2[0], &J2,0);
				if ( k == m-1 ) { _ff_set(*sJ,S2[0]); return 1; }	// must be a double root
				if ( ! k ) { printf ("Error5, k=%d %d-isogenies found for j-invariant %ld over F_%ld, expected >= %d\n", k+1, m, _ff_get_ui(S[1]), _ff_p, m); abort(); }
				k = phi_get_neighbors (S, phi, m, S[0], S2,0);
				if ( ! k ) _ff_set(*sJ,S2[1]); else _ff_set(*sJ,S2[0]);
			} else {
				_ff_set(J2,S[1]);
				k = phi_get_neighbors (S2, phi, m, J2, &J,0);
				if ( k == m-1 ) { _ff_set(*sJ,J1); return 1; }			// must be a double root
				if ( ! k ) {
					// we know that J is at level 1 and that J1 is at level 2, so one of J1's neighbors in S1 (which excludes J) is at level 3 (and the other is at level 1)
					k = phi_get_neighbors (S, phi, m, S1[0], &J1,0);
					if ( k == m-1 ) { _ff_set(*sJ,S1[0]); return 1; }	// must be a double root
					if ( ! k ) { printf ("Error6, k=%d %d-isogenies found for j-invariant %ld over F_%ld, expected >= %d\n", k+1, m, _ff_get_ui(S[1]), _ff_p, m); abort(); }
					k = phi_get_neighbors (S, phi, m, S[0], S1,0);
					if ( ! k ) _ff_set(*sJ,S1[1]); else _ff_set(*sJ,S1[0]);
				} else {
					// we know that J is at level 2 or 3
					_ff_set(J3,S1[0]);
					k = phi_get_neighbors (S1, phi, m, S1[0], &J1,0);
					if ( k == m-1 ) { _ff_set(*sJ,J3); return 1; }			// must be a double root
					if ( ! k ) {
						// we know that J1 is at level 1 and J is at level 2.  either J2 is at level 3, or J2 is at level 1 and J's other neighbor is at level 3
						_ff_set(J3,S2[0]);
						k = phi_get_neighbors (S1, phi, m, S2[0], &J2,0);
						if ( ! k ) _ff_set(*sJ,S[2]); else _ff_set(*sJ,J2);		// if k == 0 then J2 is at level 1 and o.w. not
					} else {
						_ff_set(J3,S2[0]);
						k = phi_get_neighbors (S1, phi, m, S2[0], &J2,0);
						if ( k == m-1 ) { _ff_set(*sJ,J3); return 1; }		// must be a double root
						if ( ! k ) { _ff_set(*sJ,J1); return 1; } 			// we know that J2 is at level 1, J is at level 2, and J1 cannot be at level 1 so it must be at level 3
						// neither J1 nor J2 are at level 1, so J cannot be at level 2 and must be at level 3
						_ff_set(*sJ,J);
					}
				}
			}
		}
	} else {
		// general case -- performs an expected binom(h,2)*m/2 factorizations of phi_m(x,j)
		up = -1; j = 0;
		if ( k == 1 ) {
			_ff_set(J1,S[0]);
			k = phi_get_neighbors(S,phi,m,J1,&J,0);
			_ff_set(J,J1);
			level = 1;
		} else {
			i = phi_test_height (J,S[0],h,phi,m,0);
			j = phi_test_height (J,S[1],h,phi,m,0);
			level = ( i < j ? i : j );
			if ( level >= h ) { _ff_set(*sJ,J); return 1; }					// if two attempts failed to hit the floor then we must be at the surface
			if ( i < j ) up = 1;
			if ( j < i ) up = 0;
			if ( up < 0 ) j = 2; else j = 0;								// minor optimization: remember that the first two neighbors are down and avoid retesting them
		}
		for(;;) {
			// we now know the level of J, which is greater than 0, and it has k neigbors in S
			if ( up < 0 ) {
				// if we don't know which way is up, search for a neighbor which doesn't lead down
				for ( up = j ; up < k-1 ; up++ ) if ( phi_test_height (J,S[up],level,phi,m,0) > level ) break;
				// note that if the first k-1 neighbors are down, the kth must be up
			}
			if ( ++level == h ) break;
			_ff_set(J1,S[up]);
			k = phi_get_neighbors(S,phi,m,J1,&J,0);
			_ff_set(J,J1);
			up = -1;  j = 0;
		}
		_ff_set(*sJ,S[up]);
	}
	return 1;
}

/*
debug code for ff_surface_path

				printf("Prevroot: %ld\n", _ff_get_ui(W[w-2]));
				printf("Roots: "); for ( i = 0 ; i < k ; i++ ) printf("%ld ", _ff_get_ui(W[w+i])); puts("");
				phi_eval_ff(T[0],phi,m,J); ff_poly_print(T[0],m+1); 
				ff_poly_remove_root (T[0],T[0],m+1,W+w-2);
				ff_poly_print(T[0],m);
				printf ("ff_poly roots returned %d\n", ff_poly_distinct_roots(W+w,T[0],m));
				abort();
*/
