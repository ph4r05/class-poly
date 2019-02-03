#ifndef _PHI_POLY_INCLUDE_
#define _PHI_POLY_INCLUDE_

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

/*
	The module phi_poly contains functions that use precomputed modular polynomials Phi_m(X,Y) for reasonably small values of m < PHI_MAX_M,
	including functions for navigating isogeny volcanoes.  Except in a few cases (e.g. the load functions), m must be prime.
	
	See modpoly.h for functions that support the computation of modular polynomials.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "ff_poly.h"
#include "iqclass.h"				// this is only here to get IQ_MAX_GENS as the array dimension for aux_map in phi_enum_roots -- this dependency could be removed

// directory for modulary polynomial files phi_*.txt
#define PHI_DIR				phi_dir()
#define PHI_DIR_NAME			"phi_files"
extern char _phi_dir_str[4096];
static inline char *phi_dir (void)
{
	char *s;
	if ( _phi_dir_str[0] ) return _phi_dir_str;
	s = getenv("CLASSPOLY_PHI_FILES");
	if ( s ) snprintf(_phi_dir_str, 4096, "%s", s);
	else {
	    s = getenv("HOME");
	    if ( s ) sprintf (_phi_dir_str, "%s/%s", s, PHI_DIR_NAME);
	    else strcpy (_phi_dir_str, PHI_DIR_NAME);
	}    
	return _phi_dir_str;
}

#define PHI_MAX_M				199			// must be at least 127
#define PHI_MAX_V				1023		// bound only used to determine precomputed shapes, v may be much larger than this
#define PHI_MAX_H				31
#define PHI_MAX_ENUM_H		10
#define PHI_MAX_SPARSE_FACTOR	24
#define PHI_MAX_ORIENT_P		199

#define PHI_MAX_V_P			8			// maximum number of distinct primes that may divide v

struct phi_vshape_struct {
	unsigned char p[PHI_MAX_V_P];			// must be big enough to hold PHI_MAX_M
	unsigned char h[PHI_MAX_V_P];
	int k;
};

extern mpz_t *Phi_mpz_mod[PHI_MAX_M+2];

// phi_header is hidden behind the first element of the arrays Phi_ff (for historical reasons).
struct phi_header {
	int32_t inv;
	int32_t sparse;
};

static inline int phi_sparse_factor (ff_t *phi) { register struct phi_header *hdr;  hdr = ((struct phi_header *)phi) - 1; return hdr->sparse; }
static inline int phi_3sparse (ff_t *phi) { return (phi_sparse_factor(phi)%3?0:1); }
static inline int phi_inv (ff_t *phi) { register struct phi_header *hdr;  hdr = ((struct phi_header *)phi) - 1; return hdr->inv; }

int phi_poly_setup(void);
void phi_vshape (struct phi_vshape_struct *s, int v);
int phi_load_mpz (mpz_t Phi[], int m, int inv);
int phi_load_ff (ff_t phi[], int m, int inv);

int phi_find_surface (ff_t *outJ, ff_t inJ, unsigned long up[], unsigned long uh[], int uw, unsigned long vp[], unsigned long vh[], int vw);
// IMPORTANT roots must have room for h+PHI_MAX_H+1 entries!
int phi_enum_roots (ff_t roots[], ff_t J0, long p[], long n[], long o[], long aux_p[], long aux_map[IQ_MAX_GENS][IQ_MAX_GENS], int r, int v, int inv);
int phi_next_surface_neighbor (ff_t nJ[1], ff_t phi[], int m, int h, ff_t J, ff_t *pJ, int inv);		// outside callers should set phi=0
int phi_neighbors (ff_t r[],  int m, ff_t J, ff_t *xJ, int inv);								// same as phi_roots except xJ may point to a neighbor that is to be excluded
static inline int phi_roots (ff_t r[],  int m, ff_t J, int inv)
	{ return phi_neighbors (r,m,J,0,inv); }
ff_t *phi_poly_ff (int m, int inv);													// returns ptr to Phj_m^g mod p, only supports prime m <= PHI_MAX_M
void phi_eval_ff (ff_t f[], ff_t phi[], int m, ff_t J);
int phi_poly_verify_edge (ff_t J1, ff_t J2, int m, int inv);
int phi_poly_verify_2path (ff_t J1, ff_t J2, int m1, int m2, int inv);

int phi_surface_path (ff_t W[], int n, ff_t phi[], int m, int h, ff_t J, ff_t *nJ, int inv);			// *** W must hold space for  n+h nodes ***  phi may be null
void phi_surface_gcd_cycle (ff_t r1[], ff_t *r2, int p1, int n, int p2, int e, ff_t *phi1, ff_t *phi2);	// Assumes n > e > 1 and that roots[0],...,roots[e-1] are already present in r2, which may equal r1 (if null, uses r2=r1)

void phi_adjust_level (ff_t *J1, ff_t J0, int k, int h, int m, int inv);
int phi_find_level(ff_t J, int h, int m, int inv);
void phi_ascend (ff_t *J1, ff_t J0, int k, int h, int m, int inv);
void phi_descend (ff_t *J1, ff_t J0, int k, int h, int m, int inv);

static inline int phi_index(int i, int j)												// returns index of X^iY^j coeff for 0<=i,j<=ell (note that X^(ell+1) and Y^(ell+1) coeff are implicit and are *not* present)
	{ return ( i>=j ? ((i*(i+1))>>1) + j : ((j*(j+1))>>1) + i ); }
	
static inline int phi_degree (int m)
{
	unsigned long p[MAX_UI_PP_FACTORS], e[MAX_UI_PP_FACTORS];
	int i, k, d;
	
	if ( ui_is_prime(m) ) return m+1;		// this is just a table lookup for primes up to 2^16 (or even 2^20)
	k = ui_factor (p, e, m);
	for ( d = m, i = 0 ; i < k ; i++ ) d = d/p[i] * (p[i]+1);
	return d;
}	

static inline int phi_count (int m) { int d;  d = phi_degree(m); return  (d*(d+1))/2; }		// upper bound on the number of distinct coefficients of Phi_m
static inline int phi_offset (long j, long m, long d, long s)								// returns the least i for which X^i*Y^j has a potentially nonzero coeff in Phi_m of degree d with sparseness factor s
	{ return (int)((d+(s-1)*m*j)%s); }												// note that we force conversion to longs to avoid worrying about overflow in the multiplication


static inline void phi_exps(int *i, int *j, int k)
{
	register int m;
	m = sqrt(2*k); if ( m ) m--;
	while ( ((m*(m+1))>>1) + m+1 < k ) m++;
	*i = m;
	*j = k - ((m*(m+1))>>1);
}

// higher level functions for managing reduction and evaluation of modular polys Phi_m
// more general and simpler, but not as fast as the more specialized functions

typedef struct phi_poly_struct {
	int m;		// mod poly Phi_m parameterizing m-isogenies
	int d;		// degree of Phi_m, equal to phi_degree(m)
	mpz_t *C;		// array of (d(d+1))/2 mpz_t coefficients (symmetric coefficients only stored once)  only loaded into memory for m <= PHI_MAX_M
	ff_t *c;		// array of (d+1)^2 ff_t coefficients (symmetric coefficients stored twice)
	ff_t *w1;		// workspace of d+1 coefficients
	ff_t *w2;		// workspace of d+1 coefficients
	long p;		// phi = Phi mod p
	FILE *fp;		// file pointer used to read coefficients on demand when m is large
	int inv;
} phi_poly_t[1];

static inline int phi_poly_degree (phi_poly_t phi) { return phi->d; }

void phi_poly_init (phi_poly_t phi, int m, int inv);
void phi_poly_clear (phi_poly_t phi);
int phi_poly_load (phi_poly_t phi, int m, int inv);
int phi_poly_create (phi_poly_t phi, int m, int inv);


static inline void ff_phi_poly_reduce (phi_poly_t phi)
{
	register ff_t *c;
	register mpz_t *C;
	register int i, j, k, n;
	
	if ( phi->p == _ff_p ) return;
	n = phi->d + 1;
	c = phi->c;  C = phi->C;
	for ( i = 0, k = 0 ; i < phi->d ; i++ ) {
		for ( j = 0 ; j <= i ; j++, k++ ) {
			_ff_set_mpz(c[n*i+j], C[k]);
			_ff_set(c[n*j+i], c[n*i+j]);			// unnecessary when i=j but it costs less to set than to check 
		}
	}
	_ff_set_one(c[phi->d]);  _ff_set_one(c[n*phi->d]);
	for ( i = 1 ; i < phi->d ; i++ ) _ff_set_zero(c[n*i+phi->d]);
	phi->p = _ff_p;
}

// evaluates Phi_m(X,Y) at X=x
static inline void ff_phi_poly_eval (ff_t f[], phi_poly_t phi, ff_t x)
{
	register int i,n;
	
	ff_phi_poly_reduce (phi);
	n = phi->d+1;
	_ff_set_one(phi->w1[0]); _ff_set (phi->w1[1],x); for ( i = 2 ; i < n ; i++ ) ff_mult(phi->w1[i],phi->w1[i-1],x);		// Set w[i]=x^i for i from 0 to d
	for ( i = 0 ; i < n ; i++ ) ff_dot_product(f+i,phi->w1,phi->c+n*i,n);
}

// evaluates d/dX Phi_m(X,Y) at X=x
static inline void ff_phi_poly_eval_derivative (ff_t f[], phi_poly_t phi, ff_t x)
{
	register int i,n;
	
	ff_phi_poly_reduce (phi);
	n = phi->d+1;
	_ff_set_one(phi->w1[0]); _ff_set (phi->w1[1],x); for ( i = 2 ; i < n ; i++ ) ff_mult(phi->w1[i],phi->w1[i-1],x);		// Set w[i]=x^i for i from 0 to d
	for ( i = 0 ; i < n; i++ ) {
		ff_poly_derivative (phi->w2, 0, phi->c+n*i, phi->d);
		ff_dot_product(f+i,phi->w1,phi->w2,phi->d);
	}
}


#endif
