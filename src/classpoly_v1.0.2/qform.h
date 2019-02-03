#ifndef _QFORM_INCLUDE_
#define _QFORM_INCLUDE_

#include <math.h>
#include <stdint.h>
#include "mpzutil.h"
#include "iqclass.h"

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

#define FIB_PRIME				11400714819323198549UL								// FIB_PRIME/2^64 is approximately (sqrt(5)-1)/2

int qform_reduce (long *pu, long *pv, long *pw);
void qform_compose (long *pu3, long *pv3, long *pw3, long u1, long v1, long w1, long u2, long v2, long w2, long L);
void qform_square (long *pu3, long *pv3, long *pw3, long u, long v, long w, long L);
static inline void qform_invert (long *pu3, long *pv3, long *pw3) { if ( *pu3 != *pv3 && *pu3 != *pw3 ) *pv3 = -*pv3; }
static inline int qform_is_order2 (long a, long b, long c) { return ( a > 1 && ( b==0 || a==b || a==c ) ); }
static inline int qform_is_2torsion (long a, long b, long c) { return ( a == 1 || b==0 || a==b || a==c ); }
static inline int qform_is_identity (long a, long b, long c) { return (a==1); }
static inline void qform_set_identity (long *pa, long *pb, long *pc, long D) { if ( !(D&3) ) { *pa = 1; *pb = 0; *pc = -D/4; } else { *pa = 1; *pb = 1; *pc=(1-D)/4; } }
static inline void qform_set (long *pa, long *pb, long *pc, long a, long b, long c) { *pa = a; *pb = b; *pc = c; }
static inline int qform_equal (long a, long b, long c, long u, long v, long w) { return (a==u&&b==v); }
static inline int qform_inverses (long a, long b, long c, long u, long v, long w) { return ( a==u && b==-v ? 1: 0 ); }
void qform_exp (long *pu, long *pv, long *pw, long u, long v, long w, long e, long L);
long qform_order (long a, long b, long c, long L);
double qform_sum_jbounds (long D, int n);
void qform_sum_jbounds_subgroup (double *psummax, double *pmaxdelta, long p[IQ_MAX_GENS], long n[IQ_MAX_GENS], int k, int d0, int d1, long e, long D, double hf);

long qform_class_number (long D);		// returns the order of the classgroup
long qform_class_exponent (long D);	// returns the exponent of the classgroup

// computes a set of small primeform generators, and estimates upper bound on the #bits in the coefficients of the Hilbert class polynomials
int qform_generators (long p[IQ_MAX_GENS], long n[IQ_MAX_GENS], long map[IQ_MAX_RLEN], long D, long ellfilter, long ell0, double *pbits, double hf);
long qform_table_lookup (long a, long b, long c);
int qform_table_entry (long *pa, long *pb, long i);
double qform_hilbert_height_bound (long D, int deg);							// caller must call qform_generators first!
long qform_classpoly_height_bound (long D, double hf, long ell0, long *pmaxc);		// caller must call qform_generators first!

long qform_next_prime (long a, long D);									// returns least prime q > a for which legendre(D,p) != -1
void qform_table_init (long N, long D);
void qform_table_free ();

int qform_small_primeform (long *pa, long *pb, long *pc, long D, long a0);					
int qform_primeform (long *pa, long *pb, long *pc, long p, long D, int rflag);					// computes the unique form with a=p, b >= 0 or returns zero if none exists.  if rflag is 0, it leaves the form unreduced
																				// if rflag is -1 it reduces it, if rflag = 1 it returns zero if the form is not already reduced
int qform_next_primeform (long *a, long *b, long *c, long *ell, long D, long B, long ellfilter);		// When B is nonzero, looks for reduced primitive primeforms up to B (of norm prime to ellfilter)
																				// Otherwise we get the next primitive non-principal primeform, even if not reduced (of norm prime to ellfilter)

int qform_least_representing_sequence (long D, long h, long *ellmax);							// returns least k s.t. every elt of cl(D) is a product of a subset of the first k primeforms (ellmax is the norm of the kth primeform)


static inline int qform_primeform_is_principal (long n, long D)
	{ long a, b, c;  qform_primeform(&a,&b,&c,n,D,-1); return qform_is_identity(a,b,c); }
int qform_nform (long *pa, long *pb, long *pc, long n, long D);								// returns reduced power-product of primeforms corresponding to the prime factorization of n
int qform_ppf_form (long *pa, long *pb, long *pc, ppf_t N, long D);							// returns reduced power-product of primeforms corresponding to the prime factorization of N
static inline int qform_nform_is_principal (long n, long D)
	{ long a, b, c;  qform_nform(&a,&b,&c,n,D); return qform_is_identity(a,b,c); }
static inline int qform_nform_is_2torsion (long n, long D)
	{ long a, b, c;  qform_nform(&a,&b,&c,n,D); return qform_is_2torsion(a,b,c); }

int qform_form_rep (long e[IQ_MAX_GENS], long a, long b, long c, long D)	;					// D must match current table D! (and for the next 3 functions that use this)
static inline int qform_prime_rep (long e[IQ_MAX_GENS], long q, long D)
	{ long a,b,c; if ( ! qform_primeform (&a,&b,&c,q,D,-1) ) return 0; return qform_form_rep (e, a, b, c, D); }
static inline int qform_iprime_rep (long e[IQ_MAX_GENS], long q, long D)
	{ long a,b,c; if ( ! qform_primeform (&a,&b,&c,q,D,-1) ) return 0; qform_invert (&a,&b,&c); return qform_form_rep (e, a, b, c, D); }
static inline int qform_n_rep (long e[IQ_MAX_GENS], long n, long D)
	{ long a,b,c; if ( ! qform_nform (&a,&b,&c,n,D) ) return 0; return qform_form_rep (e, a, b, c, D); }
static inline int qform_ppf_rep (long e[IQ_MAX_GENS], ppf_t N, long D)
	{ long a,b,c; if ( ! qform_ppf_form (&a,&b,&c,N,D) ) return 0; return qform_form_rep (e, a, b, c, D); }

	
int qform_primeform_test_exp (long p, long e[], int k, long D);								// returns 1 if primeform p exponentiated to e[0], e[1], ..., e[k-1] is the identity, 0 if not.
	
long qform_fastorder_ppf (long a, long b, long c, ppf_t M, long L);
static inline long qform_fastorder (long a, long b, long c, long m, long L) { ppf_t M;  ppf_factor(M,m); return qform_fastorder_ppf(a,b,c,M,L); }

long qform_order_naive (long a, long b, long c, long L);

long qform_primeform_order (long p, long D);
long qform_primeform_mod_discrete_log (long p, long q, long e1, long n1, long n2, long D);		// assumes dl is = e1 mod n1, where n1 divides n2, the order of the base primeform p
long qform_discrete_log (long a, long b, long c, long u, long v, long w, long n, long L);
long qform_discrete_log_ph_ppf (long a, long b, long c, long u, long v, long w, ppf_t M, long L);
long qform_discrete_log_naive (long a, long b, long c, long u, long v, long w, long n, long L);
static inline long qform_discrete_log_ph (long a, long b, long c, long u, long v, long w, long m, long L)
	{ ppf_t M;  if ( m < 16 ) return qform_discrete_log_naive(a,b,c,u,v,w,m,L); ppf_factor(M,m); return qform_discrete_log_ph_ppf(a,b,c,u,v,w,M,L); }

	
static inline long qform_primeform_discrete_log (long p, long q, long n, long D)
{
	long a,b,c, u, v, w, L;

	if ( ! qform_primeform (&a,&b,&c,p,D,-1) ) return -1;
	if ( ! qform_primeform (&u,&v,&w,q,D,-1) ) return -1;
	L = ceil(pow(-D,0.25));
	return qform_discrete_log (a,b,c,u,v,w,n,L);
}

static inline long qform_primeform_discrete_log_ph (long p, long q, long n, long D)
{
	long a,b,c, u, v, w, L;

	if ( ! qform_primeform (&a,&b,&c,p,D,-1) ) return -1;
	if ( ! qform_primeform (&u,&v,&w,q,D,-1) ) return -1;
	L = ceil(pow(-D,0.25));
	return qform_discrete_log_ph (a,b,c,u,v,w,n,L);
}


static inline void qform_primeform_exp (long *pa, long *pb, long *pc, long p, long e, long D)
{
	long L;
	
	if ( ! qform_primeform (pa,pb,pc,p,D,-1) ) { qform_set_identity(pa,pb,pc,D); return; }			// this should never fail, but we don't want to try exponentiating if it does
	L = ceil(pow(-D,0.25));
	qform_exp (pa,pb,pc,*pa,*pb,*pc,e,L);
}

static inline long qform_primeform_fastorder (long p, long D, long m)							// m is a multiple of the order
{
	long a,b,c, L;
	ppf_t M;

	if ( ! qform_primeform (&a,&b,&c,p,D,0) ) return 0;
	qform_reduce(&a,&b,&c);
	L = ceil(pow(-D,0.25));
	if ( m < 16 ) return qform_order_naive (a,b,c,L);
	ppf_factor(M,m);
	return qform_fastorder_ppf (a,b,c,M,L);
}

static inline long qform_primeform_exp_fastorder (long p, long e, long D, long m)					// exponentiates by e, then computes the order, using the multiple m
{
	long a,b,c,L;
	ppf_t M;

	if ( ! qform_primeform (&a,&b,&c,p,D,-1) ) return 0;
	L = ceil(pow(-D,0.25));
	qform_exp (&a,&b,&c,a,b,c,e,L);
	if ( m < 16 ) return qform_order_naive (a,b,c,L);
	ppf_factor(M,m);
	return qform_fastorder_ppf (a,b,c,M,L);
}

// we assume an unsigned long holds 64 bits here
static inline uint32_t qform_hash32 (long a, long b) { return (uint32_t) ((((unsigned long)(a^b))*FIB_PRIME)>>32); }

static inline long qform_c(long a, long b, long D) { return (b*b-D)/(a<<2); }

// returns true iff the products p1*p2, p1*p2^-1, p1^-1*p2, and p1^-1*p2^-1 are all distinct in cl(D)
static inline long qform_p1p2_check (long p1, long p2, long D)
{
	long a1,a2,b1,b2,c1,c2,L;
	
	L = ceil(pow(-D,0.25));
	if ( ! qform_primeform (&a1,&b1,&c1,p1,D,-1) ) return 0;
	qform_square (&a1,&b1,&c1,a1,b1,c1,L);
	if ( ! qform_primeform (&a2,&b2,&c2,p2,D,-1) ) return 0;
	qform_square (&a2,&b2,&c2, a2,b2,c2,L);
	if ( a1==a2 && (b1==b2 || b1==-b2) ) return 0;
	return 1;
}

static inline void qform_evec_print_ell (long ell[], long e[], int k)
{
	register int i;
	
	for ( i = 0 ; i < k ; i++ ) { if ( i ) printf (" "); printf ("%ld^%ld", ell[i], e[i]); }
}

#define K1				2114.567							// upper bound on |j(tau)-1/q|, from Enge p. 1094 (derived using Brisebarre + Philibert)
#define LOGK1				7.656606
#define FUDGE			0.0000000000001
#define LOGFUDGE			-29.9336

/*
	Compute an upper bound B on log|j(z)| for z = (-b+sqrt(D))/2a.  This doesn't depend on b.
	The value returned is an upper bound on log(M_k) where M_k is as in Lemma 8 of the HCP paper.
*/
static inline double jbound (long a, long D)
{
	double B;
	
	B = MPZ_PI*sqrt(-D)/(double)a + FUDGE;		// B = log|1/q|
	if ( LOGK1-B < LOGFUDGE ) {
		B += FUDGE;						// if B is huge then log(1+k1q) is smaller than FUDGE
	} else {
		B += log (1.0+K1*exp(-B));
	}
	return B;
}

#endif
