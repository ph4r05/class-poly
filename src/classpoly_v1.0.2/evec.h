#ifndef _EVEC_INCLUDE_
#define _EVEC_INCLUDE_

#include "ntutil.h"

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
	This module contains inlines for fast operations on polycyclic presentations.
	
	A polycyclic presentation (pcp) for a generic group G is specified by a triple (n,r,k),
	where n is a k-element array of integers (relative orders) and r is a binom(k,2)
	array of integers (power relations), 	with the property that
	
		a[i]^n[i] = prod_{j<i} a[j]^r[binom(i,2)+j]    for i in [0..k-1].
		
	Here (a[0],...,a[k-1]) is a polycyclic sequence for G.  The actual group elements a[i]
	are implicit and not used by any of the functions in this module -- given a pcp for G,
	we can effectively build our own black box for G without knowing the a[i].
	
	The integer n[i] is defined to be the least integer s.t. a[i]^n[i] \in <a[0],...,a[i-1]>,
	and is called the relative order of a[i].  For each c = binom(i,2) and j in [0..i-1],
	the integer r[c+j] lies in [0..n[i]-1] (and is uniquely determined).
	
	See Section 5 of "Computing Hilbert class polynomials with the Chinese Remainder Theorem",
	A.V. Sutherland, Math Comp 80 (2011), pp. 501-538 for more details.

	We assume (without verification) that all group orders fit in a 64-bit long
	This implies that we cannot have a polycyclic presentation with more than 63 generators
*/

#define EVEC_MAX_LEN	63

static inline void evec_clear (long e[], int k)
	{ register int i; for ( i = 0 ; i < k ; i++ ) e[i] = 0; }
	
static inline void evec_copy (long e1[], long e2[], int k)
	{ register int i;  if ( e1 != e2 ) for ( i = 0 ; i < k ; i++ ) e1[i] = e2[i]; }
	
static inline int evec_equal (long e1[], long e2[], int k)
	{ register int i; for ( i = 0 ; i < k && e1[i] == e2[i] ; i++ ); return ( i==k ? 1 : 0 ); }
	
static inline void evec_print (long e[], int k)
{
	register int i;
	printf ("["); for ( i = 0 ; i < k ; i++ ) { if ( i ) printf (" "); printf ("%ld", e[i]); } printf ("]");
}

static inline void evec_print_ell (long e[], long ell[], int k)
{
	register int i;
	printf ("["); for ( i = 0 ; i < k ; i++ ) { if ( i ) printf (" "); printf ("%ld^%ld", ell[i], e[i]); } printf ("]");
}

// r is an array of concatenated evecs in which the ith evec has length i (starting with i=0)
// returns a pointer to the ith evec, at offset  binom(i,2)
static inline long *evec_ri (long r[], int i) { return r+((i*(i-1))>>1); }


// given relative orders n[i], computes m[i] = prod_{j<=i} n[i] for i from 1 to k
// m[i] is the order of the subgroup spanned by the first i generators
static inline void evec_n_to_m (long m[], long n[], int k)
{
	register int i;
	
	m[0] = n[0];
	for ( i = 1 ; i < k ; i++ ) m[i] = m[i-1]*n[i];
}

// converts an evec to an integer index corresponding to the multi-radix representation
// of the evec with moduli corresponding to the subgroup orders m[i]
static inline long evec_to_index (long e[], long m[], int k)
{
	register long index;
	register int i;
	
	index = e[0];
	for ( i = 1 ; i < k ; i++ ) index += e[i]*m[i-1];
	return index;
}

// converts an index back to an evec
static inline void index_to_evec (long e[], long index, long m[], int k)
{
	register int i;
	
	for ( i = k-1 ; i > 0 ; i-- ) {
		e[i] = index / m[i-1];
		index -= e[i]*m[i-1];
	}
	e[0] = index;
}

// reduces evec e so that e[i] < n[i] (e[i] are assumed to be nonnegative) using pcp (n,r,k)
// does not check for overflow -- this could be an issue for groups larger than 2^31/k.
// TODO: consider reducing modulo the group order 
static inline void evec_reduce (long e[], long n[], long r[], int k)
{
	register long q, *ri;
	register int i, j;

	if ( ! k ) return;
	for ( i = k-1 ; i > 0 ; i-- ) {
		if ( e[i] >= n[i] ) {
			q = e[i]/n[i];
			ri = evec_ri(r,i);
			for ( j = 0 ; j < i ; j++ ) e[j] += q*ri[j];
			e[i] -= q*n[i];
		}
	}
	e[0] %= n[0];
}

// computes e3 = log(a^e1*a^e2) in terms of the given polycyclic presentation  (here a denotes the implicit vector of generators)
static inline void evec_compose (long e3[], long e1[], long e2[], long n[], long r[], int k)
{
	register int i;
	
	for ( i = 0 ; i < k ; i++ ) e3[i] = e1[i]+e2[i];
	evec_reduce (e3, n, r, k);
}

// computes e2 = log((a^e1)^x) in terms of the given polycyclic presentation (here a denotes the implicit vector of generators)
static inline void evec_exp (long e2[], long e1[], long x, long n[], long r[], int k)
{
	register int i;
	
	for ( i = 0 ; i < k ; i++ ) e2[i] = x*e1[i];
	evec_reduce (e2, n, r, k);
}

// performs group operation using indices, returns a*b
static inline long index_compose (long a, long b, long m[], long n[], long r[], int k)
{
	long e1[EVEC_MAX_LEN], e2[EVEC_MAX_LEN];

	index_to_evec (e1, a, m, k);
	index_to_evec (e2, b, m, k);
	evec_compose (e1, e1, e2, n, r, k);
	return evec_to_index (e1, m, k);
}

// e1 and e2 may overlap
// note that this function is not very efficient because it does not know the orders of the elements in the presentation, only the relative orders
static inline void evec_inverse (long e2[], long e1[], long n[], long r[], int k)
{
	long e3[EVEC_MAX_LEN], e4[EVEC_MAX_LEN];
	register int i;
	
	for ( i = 0 ; i < k ; i++ ) e4[i] = 0;
	evec_copy(e3,e1,k);
	// We have e1 + e4 = e3 which we maintain throughout while making e1 the zero vector
	for ( i = k-1 ;i >= 0 ; i-- ) {
		if ( e3[i] ) {
			e4[i] += n[i]-e3[i];  evec_reduce(e4,n,r,k);
			e3[i] = n[i]; evec_reduce(e3,n,r,k);
		}
	}
	evec_copy(e2,e4,k);
}

// e1 and e2 may overlap
// this is a faster way to compute inverses, if the presentation element orders are known (these are specified in the array o, the array n holds the relative orders)
static inline void evec_inverse_o (long e2[], long e1[], long n[], long o[], long r[], int k)
{
	register int j;
	
	for ( j = 0 ; j < k ; j++ ) e2[j] = ( e1[j] ? o[j]-e1[j] : 0 );
	evec_reduce(e2,n,r,k);
}

// computes the exponent vector of the group element g=g1^-1*g2 (so that g1*g=g2), where i1 is the index of g1 and i2 is the index of g2 in the multi-radix representation
static inline void evec_index_delta (long e[], long i1, long i2, long m[], long n[], long o[], long r[], int k)
{
	long e1[EVEC_MAX_LEN], e2[EVEC_MAX_LEN];

	index_to_evec (e1, i1, m, k);  index_to_evec (e2, i2, m, k);
	evec_inverse_o (e1, e1, n, o, r, k);
	evec_compose (e, e1, e2, n, r, k);
}

// returns 1 if the e is the 0 vector
static inline int evec_trivial (long e[], int k)
{
	register int i;
	
	for ( i = 0 ; i < k ; i++ ) if ( e[i] ) return 0;
	return 1;
}

// Computes the order of the group element a^e using the pcp (n,r,k)
static inline long evec_order (long e[], long n[], long r[], int k)
{
	long f[EVEC_MAX_LEN];
	register long o, m;
	register int i, j;
	
	evec_copy (f, e, k);
	for ( o = 1, i = k-1 ; i >= 0 ; i-- ) {
		if ( f[i] ) {
			m = n[i]/ui_gcd(f[i],n[i]);
			for ( j = 0 ; j < k ; j++ ) f[j] *= m;
			evec_reduce (f, n, r, k);
			o *= m;
		}
	}
	return o;
}

// Computes orders o[] for each generator using relative orders n[] and power relations r[]
static void inline evec_orders (long o[], long n[], long r[], int k)
{
	long e1[EVEC_MAX_LEN];
	register int i;
	
	evec_clear(e1,k);
	for ( i = 0 ; i < k ; i++ ) {
		e1[i] = 1;  if ( i ) e1[i-1] = 0;
		o[i] = evec_order(e1,n,r,k);
	}
}

// computes r for inverse generators (relative orders don't change)
static void inline evec_invert_generators (long ri[], long o[], long n[], long r[], int k)
{
	long e1[EVEC_MAX_LEN], *me;
	register int i,j;
	
	for ( i = 0 ; i < k ; i++ ) {
		me = evec_ri(r,i);
		for ( j = 0 ; j < i ; j++ ) e1[j] = o[j] - me[j];  while ( j < k ) e1[j++] = 0;
		evec_reduce (e1,n,r,k);
		me = evec_ri(ri,i);
		for ( j = 0 ; j < i ; j++ ) me[j] = e1[j];
}
}

// replaces e by e+1, where e is interpreted as a multi-radix integer modulo the relative orders n[i], returns zero if the result is 0, 1 ow
static inline int evec_inc (long e[], long n[], int k)
{
	register int i;
	
	e[0]++;														// increment low digit
	for ( i = 0 ; i < k-1 && e[i] == n[i] ; i++ ) { e[i] = 0; e[i+1]++; }			// propogate carries
	if ( e[k-1]==n[k-1] ) { e[k-1] = 0; return 0; }
	return 1;
}

// replaces e by e+1, where e is interpreted as a multi-radix integer modulo the prefix subgroup specified by d0, d1, and h, where d0 <= d1 < k and h divides n[d1]
static inline int evec_inc_subgroup (long e[], long n[], int d0, int d1, long h)
{
	register int i;

	if ( d1 > d0 ) e[d0]++; else e[d0] += h;								// increment low digit
	for ( i = d0 ; i < d1-1 && e[i] == n[i] ; i++ ) { e[i] = 0; e[i+1]++; }		// propogate carries
	if ( d1 > d0 && e[d1-1] == n[d1-1] ) { e[d1-1] = 0; e[d1] += h; }
	if ( e[d1]==n[d1] ) { e[d1] = 0; return 0; }
	return 1;
}


static inline void evec_print_relation (long n[], long r[], int i)
{
	register long *e;
	register int j;
	
	e = evec_ri (r,i);
	printf ("a[%d]^%ld = ", i, n[i]);
	for ( j = 0 ; j < i ; j++ ) if ( e[j] ) printf ("a[%d]^%ld ", j, e[j]);
	puts ("");
}


// given polycyclic presentation (n,r,k) and a prime power q=p^e that is a multiple of the order of the p-Sylow subgroup,
// computes a polycyclic presentation (sn,sr,sk) for the p-Sylow with basis (a_g[0],...,a_g[sk-1]) and returns sk
int evec_sylow_presentation (int g[], long sn[], long sr[], int k, long n[], long r[], long q);

// given a polycyclic presentation (n,r,k) for a p-group, computes a basis with generators b_i = prod_j a_j^B[i,j] of order o[i],
// where B is an k'-by-k matrix stored in a linear array of length k^2 (k' <= k is the rank, but we initially need space for a k by k matrix)
// returns the number of generators, which is equal to the p-rank.
// NB: the array r is *trashed* during this computation, the caller should save a copy of r if needed
int evec_sylow_basis (long B[], long n[], long r[], int k);

#endif
