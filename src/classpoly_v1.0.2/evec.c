#include <stdlib.h>
#include <stdio.h>
#include "evec.h"
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

// given polycyclic presentation (n,r,k) and a prime power q=p^e that is a multiple of the order of the p-Sylow subgroup,
// computes a polycyclic presentation (sn,sr,sk) for the p-Sylow with basis (a_g[0],...,a_g[sk-1]) and returns sk
int evec_sylow_presentation (int g[], long sn[], long sr[], int k, long n[], long r[], long q)
{
	long gi[EVEC_MAX_LEN];
	register long mi, *e1, *e2;
	register int i, j, t, si, sk;
	
	for ( i = 0 ; i < k ; i++ ) {
		sn[i] = ui_gcd (n[i], q);
		mi = ui_inverse (n[i]/sn[i], q);
		e1 = evec_ri(r,i);  e2 = sr+(e1-r);
		for ( j = 0 ; j < i ; j++ ) e2[j] = (e1[j] * mi) % q;
	}
	sk = 0;
	for ( i = 0 ; i < k ; i++ ) {
		if ( sn[i] == 1 ) {
			// eliminant redundant generator from all subsequent relations (but do not compress or reduce them)
			e1 = evec_ri (sr,i);
			for ( j = i+1 ; j < k ; j++ ) {
				e2 = evec_ri (sr, j);
				if ( e2[i] ) {
					for ( t = 0 ; t < i ; t++ ) e2[t] = (e2[t]+e2[i]*e1[t]) % q;
					e2[i] = 0;
				}
			}
			gi[i] = -1;
		} else {
			gi[i] = sk;
			g[sk++] = i;
		}
	}
	for ( i = si = 0 ; i < k ; i++ ) {
		if ( sn[i] > 1 ) {
			sn[si] = sn[i];
			e1 = evec_ri (sr, i);  e2 = evec_ri (sr, si);
			for ( j = 0 ; j < i ; j++ ) if ( gi[j] >= 0 ) e2[gi[j]] = e1[j];
			evec_reduce (evec_ri (sr,si), sn, sr, si);
			si++;
		}
	}
	return sk;
}

// swaps row i and j in a k x k matrix B stored as a linear array of k^2 elements
// also swaps upstream relations (i.e.  those above i and j, does not worry about any others)
static inline void swap_generators (long B[], long n[], long r[], int k, int i, int j)
{
	register long x, *e, *u, *v;
	register int h;

	if ( i == j ) return;
	u = B+i*k;  v = B+j*k;
	for ( h = 0 ; h < k ; h++ ) { x = u[h]; u[h] = v[h]; v[h] = x; }
	for ( h = (i>j?i:j)+1 ; h < k ; h++ ) { e = evec_ri(r,h); x = e[i]; e[i] = e[j]; e[j] = x; }
	x = n[i]; n[i] = n[j]; n[j] = x;
//printf ("Swapped generators %d and %d\n", i, j);
}

static inline void print_matrix (long B[], int rows, int cols)
{
	register int bi, bj;
	for ( bi = 0 ; bi < rows ; bi++ ) { for ( bj = 0 ; bj < cols ; bj++ ) printf ("%ld ", B[bi*cols+bj]); puts (""); }
}


// given a pcp (n,r,k) for a p-group, computes a basis with generators b_i = prod_j a_j^B[i,j] of order o[i],
// where B is an k'-by-k matrix stored in a linear array of length k^2 (k' <= k is the rank, but we initially need space for a k by k matrix)
// returns the number of generators, which is equal to the p-rank
// destroys r in the process
int evec_sylow_basis (long B[], long n[], long r[], int k)
{
	register long q, s, si, u, v, *e1, *e2, *e3, *b1, *b2, x, N;
	register int h, i, j, nz;
	
	if ( ! k ) return 0;
	
	// set B to the identity matrix
	memset (B, 0, k*k*sizeof(B[0]));
	for ( i = 0 ; i < k ; i++ ) B[i*k+i] = 1;
	
	// compute the group order (all we need is the exponent, but we don't know it yet)
	for ( i = 1, N = n[0] ; i < k ; i++ ) N *= n[0];
	
//printf ("Starting basis with generator a[%d] of order %ld\n", 0, n[0]);
	
	// we know a[0] is a basis for <a[0]> with order n[0].  we will extend this to a basis for all of G (possibly modifying and/or swapping elements as we go)
	for ( i = 1 ; i < k ; i++ ) {

//printf ("Presentation relations at i=%d:\n", i);
//for ( j = 0 ; j < k ; j++ ) if ( j < i ) printf ("a[%d]^%ld = 1 (basis element)\n", j, n[j]); else evec_print_relation (n, r, j);	
//printf ("Generator matrix:\n");
//print_matrix (B, k, k);

		e1 = evec_ri(r,i);  q = n[i];
		if ( q == 1 ) continue;									// skip redundant generators
		for ( j = 0 ; j < i ; j++ ) if ( e1[j] % q ) break;
		if ( j == i ) {
//printf ("Extending basis with (modified) generator a[%d] using relation ", i); evec_print_relation (n, r, i);
			// extend basis with new element a[i]*a^-(e1/q) of order q=n[i]
			b1 = B + i*k;
			for ( j = 0 ; j < i ; j++ ) {
				if ( e1[j] ) {
					b2 = B + j*k;
					x = n[j] - e1[j] / q;
					for ( h = 0 ; h < k; h++ ) b1[h] += x*b2[h];		// note that we *do* need to look at the entire rows b1 and b2, not just up to i
				}
			}
			// update relations impacted by the change to generator i
			for ( j = i+1 ; j < k ; j++ ) {
				e2 = evec_ri(r,j);
				if ( e2[i] ) {
					for ( h = 0 ; h < i ; h++ ) if ( e1[h] ) e2[h] = (e2[h] + e2[i]*e1[h]/q) % n[h];
				}
			}
//printf ("Extended basis to %d generators (some may be redundant)\n", i+1);
		} else {
			u = n[j] / ui_gcd (n[j], e1[j]);
			swap_generators (B, n, r, k, i-1, j);						 // move the basis element we plan to swap to position i-1 so that the basis stays contiguous
			j = i-1;
			v = n[j]/u;		// v = gcd (n[j],e1[j])
			s = e1[j]/v;		// e1[j] = v*s with s prime to p
//printf ("Swapping basis element a[%d] of order %ld with generator a[%d] of order at least %ld, u=%ld, v=%ld, s=%ld\n", j, n[j], i, n[i]*u, u, v, s);
			swap_generators (B, n, r, k, i, j);
			// we want to compute si = 1/s modulo a sufficient large power of p to ensure it works as an inverse in every exponent
			si = ui_inverse (s, N);
			e2 = evec_ri(r,j);									// note that e1 still points to the i-th relation, set e2 to point to the jth relation (where j=i-1)
//printf ("Updating relation "); evec_print_relation (n, r, j);
//printf ("and also relation "); evec_print_relation (n, r, i);
			/*
				Let x and y denote the initial values of a[i] and a[i-1], respectively, prior to the swap above
				We have relations x^q = y^(r*s) * a[i-2]^e[i-2] * ... * a[0]^e[0]  and y^(r*u) = 1
				We want to replace these with relations x^(q*u) = a[i-2]^(u*e[i-2]) * ... * a[0]^(u*e[0])) and y^r = x^(q*si) * a[i-2]^(-si*e[i-2]) * ... * a[0]^(-si*e[0]))
				Note that e2 points to the location of the new relation vector for x, while e1 points to the location of the old relation vector for x, which will be the new relation vector for y
			*/
			nz = 0;
			for ( h = 0 ; h < j ; h++ ) {
				if ( e1[h] ) {
					e2[h] = (u*e1[h]) %n[h];
					e1[h] = n[h] - ((si*e1[h])%n[h]);				// reduction mod n[h] is fine here (and in the line above) because these are basis elements (h < j = i-1)
					if ( e2[h] ) nz = 1;
				} else {
					e2[h] = 0;
				}
			}
			n[j] = q*u;										// replaces n[j] with n[i]*n[j]/v = q*u
			n[i] = v;											// replaces n[i] with v, so the product n[i]*n[j] is unchanged
			e1[j] = (si*q) % n[j];									// we can't simply reduce mod n[j] here because a[j] is not necessarily a basis element!
			x = (si*q) / n[j];
			if ( x ) {
				for ( h = 0 ; h < j ; h++ ) 
					if ( e2[h] ) e1[h] = (e1[h] +x*e2[h])%n[h];		// reduction mod h[h] is fine (now we are dealing with basis elements
			}

			// if the relative order of what is now a[i] decreased, we need to reduce any upstream relations (we don't need to worry about a[j], since its relative order increased)
			// note that if n[i] is now 1, this amounts to removing a[i] from all relations 
			if ( u > 1 ) {
				for ( j = i+1 ; j < k ; j++ ) {
					e3 = evec_ri(r,j);
					if ( e3[i] >= n[i] ) {
						  x = e3[i] / n[i];  e3[i] = e3[i] % n[i];
						for ( h = 0 ; h < i-1 ; h++ )
							if ( e1[h] ) e3[h] = (e3[h] + x*e1[h]) % n[h];		// reduction for basis elements is easy
						e3[h] += x*e1[h];								// handle reduction for a[i-1], which is not necessarily a basis element
						if ( e3[h] >= n[h] ) {
							 x = e3[h] / n[h];  e3[h] = e3[h] % n[h]; 
							for ( h = 0 ; h < i-1 ; h++ )
								if ( e2[h] ) e3[h] = (e3[h] + x*e2[h]) % n[h];	// note e2 points to the (i-1)-th relation
						}
					}
				}
			}	
//printf ("Updated relation "); evec_print_relation (n, r, i-1);
//printf ("and also relation "); evec_print_relation (n, r, i);
			i -= 1+nz;										// decrement i by 1 or 2 (and then increment it) so that we will process a[i], and also a[i-1] if its relation vector is nonzero, again
		}
	}
	// now remove redundant generators
	for ( i = j = 0 ; i < k ; i++ ) {
		if ( n[i] > 1 ) {
			evec_copy (B+j*k, B+i*k, k);
			n[j++] = n[i];
		}
	}
//printf ("Computed basis of rank %d\n", j);
//print_matrix (B, j, k);
	return j;
}
