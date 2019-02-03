#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpzutil.h"
#include "evec.h"
#include "qform.h"
#include "iqclass.h"
#include "table.h"
#include "bitmap.h"
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

#define QFORM_MAX_H		50000000						// we don't expect to deal with big class groups here, if we do, we need to switch to a baby-steps giant-steps algorithm

struct qform_table_entry {
	int32_t a, b;
	double jb;
} *qform_table;
long qform_table_size;
long qform_table_D;
long qform_table_p[IQ_MAX_GENS];
long qform_table_n[IQ_MAX_GENS];
long qform_table_m[IQ_MAX_GENS];
long qform_table_k;
long qform_table_ell0;
int qform_table_ell0_order;

double qform_min_jbound;

static uint32_t qform_table_next = 1;							// don't use entry 0 (implicitly, this is always the identity

int _qform_generators (long p[IQ_MAX_GENS], long n[IQ_MAX_GENS], long map[IQ_MAX_RLEN], double *S, long D, long L, long h, long ellfilter, long ell0);

#define FUDGE	0.0000000000001
#define LOGFUDGE	-29.9336

/*
	Compute an upper bound on log|j(z)| for z = (-b+sqrt(D))/2a.  This doesn't depend on b.
	The value returned is an upper bound on log(M_k) where M_k is as in Lemma 8 of the HCP paper.
*/
static inline double qform_jbound (long a, long D)
{
	double B;
	
	B = jbound (a, D);
	if ( B < qform_min_jbound ) qform_min_jbound = B;
	return B;
}


// there is no need to insert the identity, it is implicitly present 
static inline void qform_insert (long a, long b, long c)
{
//printf ("insert <%ld,%ld,%ld> at %d\n",a,b,c, qform_table_next);
	// we could check for identity, but for the sake of speed we assume the caller avoids this (or at least doesn't do it too often)
	if ( qform_table_next >= qform_table_size ) {
		//printf ("Exceeded qform table size %ld, increasing by 25%\n", qform_table_size);
		qform_table_size = (5*qform_table_size)/4;
		qform_table = realloc(qform_table, qform_table_size*sizeof(*qform_table));
	}
	qform_table[qform_table_next].a = a;
	qform_table[qform_table_next].b = b;
	qform_table[qform_table_next].jb = qform_jbound(a,b*b-4*a*c);	// cache qform_jbound
	table_insert (qform_table_next++, qform_hash32(a,b));			// no check for duplicates
}

// hash inverses to same table index
static inline void qform_insert_hi (long a, long b, long c)
{
	// we could check for identity, but for the sake of speed we assume the caller avoids this (or at least doesn't do it too often)
	if ( qform_table_next >= qform_table_size ) {
		//printf ("Exceeded qform table size %ld, increasing by 25%\n", qform_table_size);
		qform_table_size = (5*qform_table_size)/4;
		qform_table = realloc(qform_table, qform_table_size*sizeof(*qform_table));
	}
	qform_table[qform_table_next].a = a;
	qform_table[qform_table_next].b = b;
	table_insert (qform_table_next++, qform_hash32(a,(b<0?-b:b)));			// no check for duplicates
}

// returns 0 for identity, -1 if not found
static inline long qform_lookup (long a, long b, long c)
{
	uint32_t matches[TABLE_MAX_MATCHES];
	register int i, k;
//printf ("lookup <%ld,%ld,%ld>\n",a,b,c);
	if ( qform_is_identity(a,b,c) ) return 0;
	k = table_lookup(matches, qform_hash32(a,b));
	for ( i = 0 ; i < k ; i++ ) if ( a==qform_table[matches[i]].a && b==qform_table[matches[i]].b ) break;
	return ( i < k ? (long)matches[i] : -1 );
}

// returns 0 for identity, -1 if not found, doesn't distinguish inverses
static inline long qform_lookup_hi (long a, long b, long c)
{
	uint32_t matches[TABLE_MAX_MATCHES];
	register int i, k;

	if ( qform_is_identity(a,b,c) ) return 0;
	k = table_lookup(matches, qform_hash32(a,(b<0?-b:b)));
	for ( i = 0 ; i < k ; i++ ) if ( a==qform_table[matches[i]].a && ((b==qform_table[matches[i]].b) || (-b == qform_table[matches[i]].b)) ) break;
	return ( i < k ? (long)matches[i] : -1 );
}

long qform_table_lookup (long a, long b, long c)
	{ return qform_lookup(a,b,c); }

int qform_table_entry (long *pa, long *pb, long i)
{
	long c;
	
	if ( i < 0 || i >= qform_table_next ) return 0;
	if ( ! i ) { qform_set_identity (pa, pb, &c, qform_table_D); return 1; }
	*pa = qform_table[i].a;
	*pb = qform_table[i].b;
	return 1;
}

void qform_table_init (long N, long D)
{
	dbg_printf ("qform_table_init N=%ld, D=%ld\n", N, D);
	if ( N > qform_table_size ) {
		N=(1L<<ui_len(N))-1;
		if ( qform_table ) { dbg_printf("Increasing qform table size to %ld\n", N);  qform_table_free(); }
		qform_table = malloc((N+1)*sizeof(*qform_table));
		dbg_printf ("qform_table alloced %ld bytes for N=%ld\n", (N+1)*sizeof(*qform_table), N);
		qform_table_size = N+1;
	}
	qform_table_next = 1;
	table_init(ui_len(N));
}

void qform_table_free ()
{
	if ( qform_table ) { table_free (); free (qform_table); qform_table = 0; qform_table_size = 0;  qform_table_D = 0; qform_table_k = 0; qform_table_ell0 = 0; qform_table_ell0_order = 0; }
}

int _qsort_cmp_i32 (const void *a, const void *b)
{
	int32_t *m, *n;
	
	m = (int32_t *)a;  n = (int32_t *)b;
	if ( *m < *n ) return -1;
	if ( *m > *n ) return 1;
	return 0;
}

// give logx and logy, returns (approximately) log(x+y)
static inline double log_add (double logx, double logy)
 { return ( logx > logy ? logx + log(1.0+exp(logy-logx)) : logy + log(1.0+exp(logx-logy)) ); }

// computes the product of two monic polys in log rep, output array will have f_d+g_d entries
static inline void log_poly_mult_monic (double h[], double f[], int f_d, double g[], int g_d)
{
	register int i, j;
	
	for ( i = 0 ; i < g_d ; i++ ) h[i] = 0;
	for ( j = 0; j < f_d ; i++,j++ ) h[i] = f[j];											// initialize h to lt(g) * (f-lt(f))
	for ( i = 0 ; i < f_d ; i++ ) for ( j = 0 ; j < g_d ; j++ ) h[i+j] = log_add(h[i+j],f[i]+g[j]);		// add (f-lt(f)) * (g-lt(g)) to h
	for ( j = 0 ; j < g_d ; i++, j++ ) h[i] = log_add(h[i],g[j]);								// add lt(f) * (g-lt(g)) to h
}


// O(n^2) algorithm to compute prod (X+r[i]), doesn't set monic coeff, f may equal r
void log_buildpoly (double f[], double r[], int n)
{
	register double t;
	register int i,j;

	f[n-1] = r[n-1];
	for( i = n-2 ; i >=0 ; i-- ) {
		t = r[i];
		f[i] = f[i+1]+t;
		for ( j = i+1 ; j < n-1 ; j++ ) f[j] = log_add (f[j],t+f[j+1]);
		f[j] = log_add(f[j],t);
	}
}

/*
	This works and with hf=1 is essentially the same as the value given by Lemma 8, which is returned by qform_generators.
	When hf>1 it gives a more "rigorous" upper bound than simply dividing the height bound by hf, but in fact it is much less accurate.
	It's clear why it substantially overestimates in many cases.  Consider, for example, Weber-f, where we know a priori that the
	constant term is 1, but the value computed here will be Omega(h).  What is not so clear is why the heuristic of dividing the
	height bound for H_D by the height factor hf actually works so well.

	Currently this function is not used.  It is also very slow.
*/
long qform_classpoly_height_bound (long D, double hf, long ell0, long *pmaxc)
{
	double *roots, max;
	register long i, j, h, k, d, maxi;
	
	if ( qform_table_D != D ) { err_printf ("qform table not loaded for D=%ld\n", D); exit (0); }
	if ( ell0 && (qform_table_p[0] != ell0 || qform_table_n[0] != 2 ) ) { err_printf ("specified ell0=%ld either is not the first generator or is not ramified\n", ell0); exit(0); }
	h = qform_table_next;  k = (ell0 ? 2 : 1);  d = h/k;
	roots = malloc(h*sizeof(*roots));
	roots[0] = qform_jbound(1,D)/hf;
	for ( i = k, j=1 ; i < h ; i+=k, j++ ) roots[j] = qform_jbound(qform_table[i].a,D)/hf;
	log_buildpoly(roots,roots,d);
	max = roots[0]; maxi = 0;
	for ( i = 1 ; i < d ; i++ ) if ( roots[i] > max ) { max = roots[i]; maxi = i; }
	if ( pmaxc ) *pmaxc = maxi;
	return (long)ceil(MPZ_LOG2E*max+0.5);
}

/*
	Computes b = log_2 B, where B is the bound in lemma 8 of the CRT paper.
	To facilate the computation of the square-root of class polynomials whose height bounds are scaled from this bound, the caller may optionally specified the degree (e.g. h(D)/2).
	If deg is zero, it will be treated as h(D)
*/
double qform_hilbert_height_bound (long D, int deg)
{
	double b, y, sum;
	register long i;
	
	if ( qform_table_D != D ) { err_printf ("qform table not loaded for D=%ld\n", D); abort(); }
	for ( sum = qform_jbound(1,D), i = 1 ; i < qform_table_next ; i++ ) sum += qform_jbound(qform_table[i].a,D);
	dbg_printf ("sum(log(M_i)) = %.3f bits\n", sum);
	b = exp(qform_min_jbound);
	i = (long)floor((double)(deg+1)/(b+1));
	y = ( i>0 && i < deg ? logfac(deg) - logfac(i) - logfac(deg-i) : 0 );
	dbg_printf ("deg=%d, m=%ld, log(binom(h,m)) <= %.4f\n", deg, i, y);
	b = MPZ_LOG2E*(y - i*qform_min_jbound + sum);
	info_printf ("b = log2(binom(h,m)) - m*log2(M_h) + sum(log2(M_i)) = %.4f - %ld * %.4f + %.4f = %.4f bits\n", MPZ_LOG2E*y, i, MPZ_LOG2E*qform_min_jbound, MPZ_LOG2E*sum, b);
	return b;
}


/*
	sums log|j(z)| for z = (-b+sqrt(D))/2a over the unique subgroup H presented by the prefix of (p,n) defined by d1 and e
	The subgroup H consists of elements of the form p1^e1*...*pk^ek where k=d1, 0<=e_i<n[i], and ek is a multiple of e (necessarily 0 if e=n[d1])
	Assumes table is loaded for D, ignores the presentation below d0 (used to handle ramified ell0, set d0=1)
*/
void qform_sum_jbounds_subgroup (double *psummax, double *pmaxdelta, long p[IQ_MAX_GENS], long n[IQ_MAX_GENS], int k, int d0, int d1, long e, long D, double hf)
{
	long m[IQ_MAX_GENS];
	register long subgroup, cosets, subgroup_i, coset_i, i, j, n0, n1, n2, n3, n4;
	register double x, sum, max, maxsum, summax, maxdelta, totsummax, totmaxdelta;
	
	if ( qform_table_D != D ) { err_printf ("qform table not loaded for D=%ld\n", D); exit (0); }
	if ( d0 < 0 || d0> d1 || d1 >= k || e <= 0 || (n[d1]%e) ) { err_printf ("invalid prefix subgroup in qform_sum_jbounds_subgroup d0=%d, d1=%d, e=%ld\n", d0, d1, e); exit (0); }
	evec_n_to_m (m,n,k);
	if ( m[k-1] != qform_table_next ) { err_printf ("presentation order doesn't match table in qform_sum_jbounds_subgroup (next=%d)\n", qform_table_next); exit (0); }

	if ( d0 ) n0 = m[d0-1]; else n0 = 1;
	if ( d1 ) n1 = m[d1-1]/n0; else n1 = 1;
	n2 = e;  n3 = n[d1]/e;
	n4 = m[k-1]/m[d1];
	subgroup = n1*n3;
	cosets = n2*n4;
//printf ("n0=%d, subgroup=%d, cosets=%d\n", n0, subgroup, cosets);
	totsummax = totmaxdelta = 0.0;
	for ( j = 0 ; j < n0; j++ ) {
		summax = maxsum = maxdelta = 0.0;
		max = sum = ( !j ? qform_jbound(1,D) : qform_jbound(qform_table[j].a,D)) / hf;						// handle the identity element separately
//printf ("b_%d%d=%0.2f\n", 1, 1, MPZ_LOG2E*qform_jbound(1,D));
		for ( coset_i = 0 ; coset_i < cosets ; coset_i++ ) {
			for ( subgroup_i = (coset_i?0:1) ; subgroup_i < subgroup ; subgroup_i++ ) {
				i = n0 * (subgroup_i%n1+(subgroup_i/n1)*(n1*n2)+(coset_i%n2)*n1+(coset_i/n2)*(n1*n2*n3)) + j;
//				x = qform_jbound(qform_table[i].a,D) / hf;
				x = qform_table[i].jb / hf;
//printf ("b_%d%d=%0.2f\n", coset_i+1, subgroup_i+1, MPZ_LOG2E*x*hf);
				if ( x > max ) max = x;
				sum += x;
			}
//printf ("coset %d: max=%0.2f, sum=%0.2f, delta=%0.2f\n", coset_i+1, MPZ_LOG2E*max*hf, MPZ_LOG2E*sum*hf, MPZ_LOG2E*(sum-max)*hf);
			if ( sum-max > maxdelta ) maxdelta = sum-max;
			summax += max;
			max = sum = 0;
		}
		if ( sum-max > maxdelta ) maxdelta = sum-max;
		summax += max;
		totsummax += summax;  totmaxdelta += maxdelta;
//printf ("summax=%.0f, maxdelta=%.0f, \n", MPZ_LOG2E*summax, MPZ_LOG2E*maxdelta);
	}
	*psummax = totsummax/j; *pmaxdelta = totmaxdelta/j;			// heuristically, taking the average works well (not clear why)
}


/*
	sums log|j(z)| for z = (-b+sqrt(D))/2a over the n largest values of z occuring for <a,b,c> in cl(D)
	Assumes table is loaded for D
*/
double qform_sum_jbounds (long D, int n)
{
	int32_t *alist;
	double s;
	register int i;
	
	if ( qform_table_D != D ) { err_printf ("qform table not loaded for D=%ld\n", D); exit (0); }
	if ( n > qform_table_next ) { err_printf ("n=%d exceeds h(D)=%d\n", n, qform_table_next); exit (0); }
	alist = malloc(qform_table_next*sizeof(int32_t));
	alist[0] = 1;
	for ( i = 1 ; i < qform_table_next ; i++ ) alist[i] = qform_table[i].a;
	qsort(alist, qform_table_next, sizeof(int32_t), _qsort_cmp_i32);
	for ( i = 0, s=0.0 ; i < n ; i++ ) s += qform_jbound (alist[i],D);
	free(alist);
	return s;
}

// soft-O(h(D_0)) algorithm, where D_0 is the fundamental discriminant dividing D
long qform_class_number (long D)
{
	unsigned long vp[MAX_UI_PP_FACTORS], vh[MAX_UI_PP_FACTORS];
	long p[IQ_MAX_GENS], n[IQ_MAX_GENS], h, v, D0, L, hB;
	int i, j, k, r;
	
	k = discriminant_factor_conductor(vp,vh,D);
	if ( k < 0  ) { printf ("Invalid discriminant %ld in qform_class_number\n", D); exit (0); }
	if ( k > 0 ) {
		for ( v = 1, i = 0; i < k ; i++ ) for ( j = 0 ; j < vh[i] ; j++ ) v *= vp[i];
		D0 = D/(v*v);
	} else {
		D0 = D;
	}
	
	if ( qform_table_D != D0 ) {
		L = ceil(pow(-D0,0.25));
		hB = sqrt(-D0);													// start at about twice the average and increase if needed
		if ( hB > QFORM_MAX_H ) hB = QFORM_MAX_H;
		qform_table_init(hB,D0);
		r = _qform_generators (p, n, 0, 0, D0, L, 0, 0, 0);
		for ( i = 0, h=1 ; i < r ; i++ ) h *= n[i];
	} else {
		h = qform_table_k ? qform_table_m[qform_table_k-1] : 1;
	}
	
	if ( k > 0 ) {		
		for ( i = 0 ; i < k ; i++ ) {
			j = kronecker_p(D0,vp[i]);
			h *= vp[i] - j;
			for ( j = 1 ; j < vh[i] ; j++ ) h *= vp[i];
		}
		if ( D0 == -3 ) {
			if ( (h%3) ) { printf ("Error, non-integral h(D) computed for D=%ld\n", D); exit (0); }
			h /= 3;
		} else if ( D0 == -4 ) {
			if ( (h&1) ) { printf ("Error, non-integral h(D) computed for D=%ld\n", D); exit (0); }
			h >>= 1;
		}
	}
	dbg_printf ("h(%ld) = %ld\n", D, h);
	return h;
}

// soft-O(h(D)) algorithm to compute the exponent of the classgroup
long qform_class_exponent (long D)
{
	long p[IQ_MAX_GENS], n[IQ_MAX_GENS], o[IQ_MAX_GENS], r[IQ_MAX_RLEN], L, e, hB;
	int i, k;
	
	// compute a polycyclic presentation for the classgroup
	L = ceil(pow(-D,0.25));
	hB = sqrt(-D);													// start at about twice the average, will be increased if needed
	if ( hB > QFORM_MAX_H ) hB = QFORM_MAX_H;
	qform_table_init(hB,D);
	k = _qform_generators (p, n, r, 0, D, L, 0, 0, 0);
	if ( ! k ) return 1; else if ( k == 1 ) return n[0];
	
	// compute the orders of each element in the polycyclic presentation
	evec_orders (o, n, r, k);
	
	// compute the LCM of the orders in the polycyclic presentation
	// In a finite abelian group, the LCM of the orders of any set of generators is equal to the LCM of the orders of all elements in the group, i.e. the group exponent.
	for ( e = o[0], i = 1 ; i < k ; i++ ) e = ui_lcm (e, o[i]);
	return e;
}



/*
	Given a discriminant D, finds a set of generators with small prime norms for the class group.
	Returns the number of generators, their norms p[i], and integers h[i]>=2 with the property
	that every element of the class group can be written uniquely in the form \prod f_i^{e_i}
	where f_i is the primeform with norm p[i] and e_i is a nonnegative integer less than h[i].
	
	Note that h[i] may be less than the order of f_i---the set of primeforms generates the
	classgroup but may not be a basis for it.
	
	The primeform f_i is (p[i],b,c) where 0<=b<=p[i], with b^2 = D mod p[i] and b^2 = D mod 4.
	(this uniquely determines b, and then c=(b^2-D)/4p[i]).

	*pbits is set to an upper bound on the # of bits in the largest coefficient in the Hilbert class polynomial H_D
	if hf is nonzero, the logarithmic height bound on each root is reduced by this factor when computing *pbits

	returns the number of generators
*/
int qform_generators (long p[IQ_MAX_GENS], long n[IQ_MAX_GENS], long map[IQ_MAX_RLEN], long D, long ellfilter, long ell0, double *pbits, double hf)
{
	long hB, L, a, b, c, h, h2, v;
	int i, k, w;
	double x,y;

	v=discriminant_conductor(D);
	if ( ! v ) { err_printf ("Invalid discriminant %ld\n", D); exit(0); }
	h = ( v == 1 ? 	0 : qform_class_number(D) );

	L = ceil(pow(-D,0.25));
	hB = sqrt(-D);														// start at about twice the average and increase if needed
	if ( hB > QFORM_MAX_H ) hB = QFORM_MAX_H;
	qform_table_init(hB,D);
	k = _qform_generators (p, n, map, pbits, D, L, h, ellfilter, ell0);
	// if we didn't compute h and ellfilter > 1 we need to verify that we got all of cl(D)
	// in this case we know that D is fundamental and we only need to check that the reduced prime forms 
	// excluded by ellfilter can be obtained from our set of generators
	if ( !h && ellfilter > 1 ) {
		unsigned long ps[MAX_UI_PP_FACTORS], hs[MAX_UI_PP_FACTORS];
		w = ui_factor(ps,hs,ellfilter);
		for ( i = 0 ; i < w ; i++ ) {
			if ( ! qform_primeform(&a,&b,&c,ps[i],D,1) ) continue;
			if ( qform_lookup(a,b,c) < 0 ) break;		// note that even if k=0, qform_lookup will return 0 for the identity
		}
		if ( i < w ) {
			info_printf ("Couldn't generate cl(D) with ellfilter=%ld, computing h(D) and switching to non-reduced primeforms.\n", ellfilter);
			qform_table_free();
			h = qform_class_number(D);
			qform_table_init(hB, D);
			k = _qform_generators (p, n, map, pbits, D, L, h, ellfilter, ell0);
		}
	}
	
//	qform_table_free();		don't free here, caller may want to use the table

	if ( ! pbits ) return k;
/*
x = 1.0;	
for ( i = 1 ; i < qform_table_next ; i++ ) x += 1.0/qform_table[i].a;
x *= MPZ_PI*sqrt(-D);
printf ("pi*sqrt(D)*sum 1/a = %.2f\n", x);
*/
	
	/*
		*pbits currently holds log(\prod M_k), where the M_k are as in lemma 8 of the HCP paper.
		
		We want to compute lg(B), where B = binom(h,m) * M_h^{-m}\prod M_k and M_h is the least M_k,
		which is stored in qform_min_jbound, and m = floor((h+1)/(M_h+1)).
	
	*/
// 	dbg_printf ("sum(log|j(a)|) = %.4f\n", *pbits);
	dbg_printf ("sum(log(M_i)) = %.4f bits\n", *pbits);
	if ( ! h ) for ( i = 0, h=1 ; i < k ; i++ ) h*=n[i];			// compute class number if we haven't already
	x = exp(qform_min_jbound);
	dbg_printf ("M_h = %.4f\n", x);
	h2 = ( ell0 ? h/2 : h );								// reduce h by factor of 2 if ell0 is set
//	i = ( x > 2 ? floor((double)(h+1)/(x+1)) : h/2 );
	i = floor((double)(h2+1)/(x+1));
	y = ( i>0 && i < h2 ? logfac(h2) - logfac(i) - logfac(h2-i) : 0 );
	dbg_printf ("h=%ld, m=%d, log(binom(h,m)) <= %.4f\n", h2, i, y);
	x = MPZ_LOG2E*(y - i*qform_min_jbound + *pbits);
	info_printf ("b = log2(binom(h,m)) - m*log2(M_h) + sum(log2(M_i)) = %.4f - %d * %.4f + %.4f = %.4f bits\n", MPZ_LOG2E*y, i, MPZ_LOG2E*qform_min_jbound, MPZ_LOG2E*(*pbits), x);
	*pbits = x / hf;										// this is much more accurate then scaling the root bounds (why is this ?!)
	dbg_printf ("height bound= %.0f bits\n", *pbits);
	return k;
}

long qform_next_prime (long p, long D)			// returns least prime q > a for which legendre(D,p) != -1 and q does not divide the conductor of D
{	
	long v;
	
	v = discriminant_conductor(D);
	if ( v <= 0 ) { err_printf ("Invalid discriminant D=%ld in qform_next_prime\n", D); exit (0); }
	do {
		p = ui_next_prime(p);
	} while ( !(v%p) || kronecker_p(D,p) == -1 );
	return p;
}

int qform_form_rep (long e[IQ_MAX_GENS], long a, long b, long c, long D)
{
	long index;
	
	if ( D != qform_table_D ) { err_printf ("D=%ld does not match table_D=%ld in qform_form_rep!\n", D, qform_table_D); abort(); }
	index = qform_lookup(a,b,c);
	if ( index < 0 ) return 0;
	index_to_evec (e, index, qform_table_m, qform_table_k);
	return 1;	
}

/*
	When B is specified, we only look for reduced primitive primeforms up to B (subject to ellfilter)
	Otherwise we get the next primitive non-principal primeform, even if not reduced (again subject to ellfilter)
*/
int qform_next_primeform (long *a, long *b, long *c, long *ell, long D, long B, long ellfilter)
{
	if ( B ) {
		while ( *ell <= B && qform_small_primeform(a,b,c,D,*ell) ) {
			*ell = *a;
			if ( *ell <= B && (ellfilter%(*ell)) ) return 1;
		}
		return 0;
	} else {
		// if no bound is specified, get next primeform, whether its reduced or not
		for ( *ell = ui_next_prime(*ell) ; ; *ell = ui_next_prime(*ell) ) {
			if ( ! (ellfilter%(*ell)) ) continue;
			if ( ! prime_to_conductor(D,*ell) ) continue;
			if ( qform_primeform(a,b,c,*ell,D,0) ) {
				qform_reduce(a,b,c);
				if ( ! qform_is_identity(*a,*b,*c) ) break;
			}
		}
		return 1;
	}
}

// Implements Algorithm 2.1 of Hilbert CRT paper
int _qform_generators (long p[IQ_MAX_GENS], long n[IQ_MAX_GENS], long map[(IQ_MAX_GENS*(IQ_MAX_GENS+1)) / 2], double *S, long D, long L, long h, long ellfilter, long ell0)
{
	long m[IQ_MAX_GENS];
	long a, b, c, ell, u, v, w, x, y, z, index, B, h0;
	register double s;
	register int i, k, next, top;

	dbg_printf ("_qform_generators D=%ld, h=%ld, ell0=%ld\n", D, h, ell0);
	qform_table_D = 0;
	qform_table_k = 0;
	qform_min_jbound = -D;
	if ( h == 1 || (D >= -163 && (D==-3 || D==-4 || D==-7 || D==-8 || D==-11 || D==-12 || D==-16 || D==-19 || D==-27|| D==-28 || D==-43 || D==-67 || D==-163)) )
		{ if ( S ) *S = qform_jbound(1, D); qform_table_D = D; return 0; }
	if ( ! h ) B = 6*log(-D)*log(-D);	 else B =0;											// use ERH Bach bound if class number is not given
	if ( ! ellfilter ) ellfilter = 1;
	ell = 1;
	if ( ell0 ) {
		if ( ui_is_small_prime(ell0) ) {
			if ( ! qform_primeform (&a,&b,&c,ell0,D,-1) ) { err_printf ("couldn't create primeform for ell0=%ld in cl(%ld)\n", ell0, D); exit (0); }
		} else {
			if ( ! qform_nform (&a,&b,&c,ell0,D) ) { err_printf ("couldn't create primeform for ell0=%ld in cl(%ld)\n", ell0, D); exit (0); }
		}
		qform_table_ell0 = ell0;
		ell = ell0;
	}
	if ( ! ell0 || qform_is_identity(a,b,c) ) {
		if ( ! qform_next_primeform (&a,&b,&c,&ell,D,B,ellfilter) ) {
			if ( ! ellfilter ) { err_printf ("couldn't find a generator with B=%ld and ellfilter=%ld\n", B, ellfilter);  exit(0); }
			if ( S ) *S = 1.0;
			return 0;
		}
	}
	//  Populate table with cyclic group generated by the smallest primeform
	//  Note that we don't insert the identity element, but it is implicitly present as entry 0 (qform_lookup checks for the identity)
	s = qform_jbound(1, D);														// a=1 for the identity element
//dbg_printf ("%d: (%ld,%ld,%ld)\n", qform_table_next, a, b, c);
	qform_insert (a,b,c);  s += qform_jbound(a,D); 
	qform_square (&u,&v,&w,a,b,c,L);
	for ( i = 2 ; ! qform_is_identity(u,v,w) ; i++ ) {
//dbg_printf ("%d: (%ld,%ld,%ld)\n", qform_table_next, u, v, w);
		qform_insert (u,v,w);  s += qform_jbound(u,D); 
		qform_compose (&u, &v, &w, u, v, w, a, b, c, L);
	}
//dbg_printf("primeform %ld has order %d\n", ell, i);
	p[0] = ell; n[0] = i; m[0] = i;
	k = 1;  h0 = i;
	if ( ell0 ) { qform_table_ell0_order = h0; ell = 1; }
	// Now compute relative orders for all primeforms of norm <= B, expanding the table to included all spanned elements as we go
	while ( (!h || h0 < h) && qform_next_primeform(&a,&b,&c,&ell,D,B,ellfilter) ) {
//dbg_printf ("pa=%ld, h0=%ld, h=%ld\n", ell, h0, h);
		u = a; v=b; w=c;  i=1;
		top = qform_table_next;
		while ( (index=qform_lookup(u,v,w)) < 0 ) {
//dbg_printf ("%d: (%ld,%ld,%ld)\n", qform_table_next, u, v, w);
			// for each f=(u,v,w)=(a,b,c)^i that is not in the current subgroup H, add the coset fH to the table
			qform_insert(u,v,w);  s += qform_jbound(u,D);									// identity element is not in table
			for ( next = 1 ; next < top ; next++ ) {
				qform_compose (&x, &y, &z, u, v, w, qform_table[next].a, qform_table[next].b, qform_c(qform_table[next].a,qform_table[next].b,D), L);
//dbg_printf ("%d: (%ld,%ld,%ld)\n", qform_table_next, x, y, z);
				qform_insert (x,y,z);  s += qform_jbound(x,D); 
			}
			if ( i++ == 1 ) {
				qform_square (&u, &v, &w, a, b, c, L);
			} else {
				qform_compose (&u, &v, &w, u, v, w, a, b, c, L);
			}
		}
		if ( i > 1 ) {
//dbg_printf ("primeform %ld matched table entry %d\n", a, index);
			if ( k >= IQ_MAX_GENS ) { err_printf ("Exceeded IQ_MAX_GENS=%d in qform_generators\n", IQ_MAX_GENS); exit (0); }
			if ( map ) { index_to_evec (evec_ri(map,k), index, m, k);  }
			m[k] = m[k-1]*i;
			p[k] = ell;  n[k++] = i;
			h0 *= i;
		}
	}
	if ( S ) *S = s;
	if ( h && h0 != h ) { err_printf ("Error, h0=%ld not equal to h=%ld in _qform_generators\n", h0, h); exit(0); }
	qform_table_k = k;
	for ( i = 0 ; i  < k ; i++ ) { qform_table_p[i] = p[i]; qform_table_n[i] = n[i]; qform_table_m[i] = m[i]; }
	evec_n_to_m (qform_table_m, qform_table_n, k);
	qform_table_D = D;
	return k;
}

/*
not currently used

int qform_new_generators (long new_p[IQ_MAX_GENS], long new_n[IQ_MAX_GENS], long ell0, long ellfilter, long n[IQ_MAX_GENS], long map[IQ_MAX_GENS][IQ_MAX_GENS], int k)
{
	long h0, h, D;
	long a, b, c, m[IQ_MAX_GENS];
	long i, j, ord, new_k, ell;
	
	if ( ! qform_table_D ) { err_printf ("qform_new_generators called without qform_table_D set, must call qform_generators first (and don't free the table)\n"); exit (0); }
	D = qform_table_D;
	h = qform_table_next;
	if ( ! qform_primeform (&a,&b,&c,ell0,D,0) ) { err_printf ("Invalid ell0=%ld in qform_new_generators, couldn't construct primeform with D=%ld\n", ell0, D); exit (0); }
	if ( (i = qform_lookup(a,b,c)) < 0 ) { err_printf ("Couldn't find primeform (%ld,%ld,%ld) in table in qform_new_generators\n", a, b, c); exit (0); }
	bitmap_alloc(h);
	evec_n_to_m (m, n, k);
	for ( j = i, ord = 1 ; ! bitmap_test(j) ; ord++ ) { bitmap_set(j);  j = qform_index_compose(j, i, m, n, map, k); }
	new_p[0] = ell0; new_n[0] = ord;  h0 = h;  ell = ell0;  new_k = 1;
	while ( h0 < h ) {
		if ( ! qform_next_primeform(&a,&b,&c,&ell,D,0,ellfilter) ) break;
		if ( (i = qform_lookup(a,b,c)) < 0 ) { err_printf ("Couldn't find primeform (%ld,%ld,%ld) in table in qform_new_generators\n", a, b, c); exit (0); }
		for ( j = i, ord = 1 ; ! bitmap_test(j) ; ord++ ) { bitmap_set(j);  j = qform_index_compose(j, i, m, n, map, k); }
		if ( new_k >= IQ_MAX_GENS )  { err_printf ("Exceeded IQ_MAX_GENS=%d in qform_new_generators\n", IQ_MAX_GENS); exit (0); }
		new_p[new_k] = ell;  new_n[new_k++] = ord;
	}
	bitmap_free();
	if ( h0 != h ) { err_printf ("Error, h0=%ld not equal to h=%ld in qform_new_generators\n", h0, h); exit(0); }
	return new_k;
}
*/


// Determines the least integer k such that every element of cl(D) can be expressed as a product of a subset of the first k primeforms.
// ellmax is set to the norm of the kth prime form.  Used to test the conjecture in "A Pollard-type algorithm for finding short product
// representations in finite groups" by Bisson and Sutherland.
int qform_least_representing_sequence (long D, long h, long *ellmax)
{
	long a, b, c, ell, u, v, w,  n, L;
	register int k, next, top;

	if ( h == 1 ) { if ( ellmax ) *ellmax = 0;  return 0; }
	qform_table_init (h, D);
	L = ceil(pow(-D,0.25));
	ell = 1;  n = 1; 
	for ( k = 1 ;; k++ ) {
		if ( ! qform_next_primeform(&a,&b,&c,&ell,D,0,1) ) { printf ("unbounded call to qform_next_primeform failed with ell=%ld!\n", ell); abort(); }
		top = qform_table_next;
		if ( qform_lookup (a,b,c) < 0 ) { qform_insert (a,b,c); if ( ++n == h ) break; }
		for ( next = 1 ; next < top ; next++ ) {
			qform_compose (&u, &v, &w, a, b, c, qform_table[next].a, qform_table[next].b, qform_c(qform_table[next].a,qform_table[next].b,D), L);
			if ( qform_lookup (u,v,w) < 0 ) { qform_insert (u,v,w); if ( ++n == h ) break; }
		}
		if ( n == h ) break;	
	}
	qform_table_free ();
	if ( ellmax ) *ellmax = ell;
	return k;
}


// primeform needn't be reduced
long qform_primeform_order (long p, long D)
{
	long a,b,c,L;
	
	if ( ! qform_primeform (&a,&b,&c,p,D,-1) ) return 0;
	L = ceil(pow(-D,0.25));
	return qform_order(a,b,c,L);
}

int qform_primeform_test_exp (long p, long e[], int k, long D)
{
	long a,b,c,L;
	register long i;
	
	if ( ! qform_primeform (&a,&b,&c,p,D,-1) ) return -1;
	L = ceil(pow(-D,0.25));
	for ( i = 0 ; i < k ; i++ ) {
		qform_exp (&a,&b,&c,a,b,c,e[i],L);
		if ( qform_is_identity (a,b,c) ) return 1;
	}
	return 0;
}

long qform_nform_order (long n, long D)
{
	long a,b,c,L;
	
	if ( ! qform_nform (&a,&b,&c,n,D) ) return 0;
	L = ceil(pow(-D,0.25));
	return qform_order(a,b,c,L);
}

// assumes h is nonzero, p^h fits in a long
static inline long _qform_pp_fastorder (long a, long b, long c, long p, int h, long L)
{
	long u, v, w, n;
	int i;

	if ( qform_is_identity (a, b, b) ) return 1;
	if ( h==1 ) return p;
	qform_exp (&u, &v, &w, a, b, c, p, L);
	if ( qform_is_identity (u, v, w) ) return p;
//if ( h==1 ) { err_printf ("didn't get expected identity: <%ld,%ld,%ld>\n", u,v,w); exit (0); }
	for ( i = 1, n = p*p ; i < h-1 ; i++, n *= p ) {
		qform_exp (&u, &v, &w, u, v, w, p, L);
		if ( qform_is_identity (u, v, w) ) return n;
	}
	return n;
}

// CLG fastorder, not optimized at all
long qform_fastorder_ppf (long a, long b, long c, ppf_t M, long L)
{
	ppf_t M1, M2;
	long u, v, w, m1, m2;

	if ( M->w == 1 ) return _qform_pp_fastorder (a, b, c, M->p[0], M->h[0], L);
	ppf_split (M1, M2, M);
	m1 = ppf_eval(M1);  m2 = ppf_eval(M2);
	qform_exp(&u,&v,&w,a,b,c,m2,L);
	m2 = qform_fastorder_ppf (u, v, w, M1, L);
	qform_exp(&u,&v,&w,a,b,c,m1,L);
	m1 = qform_fastorder_ppf (u, v, w, M2, L);
	return m1*m2;
}



long qform_order_naive (long a, long b, long c, long L)
{
	long u, v, w;
	register long i;
	
	if ( qform_is_identity(a,b,c) ) return 1;
	if ( ! b ) return 2;
	qform_square(&u,&v,&w,a,b,c,L);
	for ( i = 2 ; ; i++ ) {
		if ( qform_is_identity(u,v,w) ) return i;
		qform_compose(&u,&v,&w,u,v,w,a,b,c,L);
	}
}

// BSGS order computation
long qform_order (long a, long b, long c, long L)
{
	long bu, bv, bw, gu, gv, gw;
	register long i,j;

	if ( L < 8 ) return qform_order_naive(a,b,c,L);
	qform_table_init(12*L*log(L),0);						// L=D^{1/4}, Schur bound is h<=D^{1/2}log(D), Terr BSGS needs at most sqrt(h) = 2sqrt(2)sqrt(h) < 12*L*log*(L). (could use Littlewood/Shanks bound here)
	bu=gu=a; bv=gv=b; bw=gw=c;
	for ( i = 1 ;; i++ ) {
		if ( qform_is_identity(bu,bv,bw) ) { qform_table_free(); return i; }
		if ( (j=qform_lookup(gu,gv,gw)) >= 0 ) { qform_table_free(); return (i*(i+1))/2-j; }
		qform_insert(bu,bv,bw);
		qform_compose(&bu,&bv,&bw,bu,bv,bw,a,b,c,L);
		qform_compose(&gu,&gv,&gw,gu,gv,gw,bu,bv,bw,L);
	}
}

// brute force for testing and small base cases
long qform_discrete_log_naive (long a, long b, long c, long u, long v, long w, long n, long L)
{
	long x,y,z;
	register long i;
	
	if  ( qform_is_identity(u,v,w) ) return 0;
	if ( ! n ) return -1;
	if ( a==u && b==v ) return 1;
	if ( n==1 ) return -1;
	qform_square(&x,&y,&z,a,b,c,L);
	for ( i = 2 ; i < n ; i++ ) {
		if ( x==u ) {
			if ( y==v ) return i;
			if ( y==-v ) return n-i;
		}
		qform_compose(&x,&y,&z,x,y,z,a,b,c,L);
	}
	return -1;
}

long qform_pp_discrete_log (long a, long b, long c, long u, long v, long w, long p,  int h, long L)
{
	long a1,b1,c1,u1,v1,w1;
	register long e1,e2,q1,q2;
	register int i,h1,h2;
	
//printf ("pp DL <%ld,%ld,%ld> <%ld,%ld,%ld> %d^%d\n",a,b,c,u,v,w,p,h);
	if ( h == 1 ) return qform_discrete_log (a,b,c,u,v,w,p,L);
	if ( h < 3 && p < 4 ) return qform_discrete_log_naive (a,b,c,u,v,w,p*p,L);
	h1 = (h+1)/2;  h2 = h-h1;
	for ( q1 = p, i = 1 ; i < h1 ; i++ ) q1 *=p;
	qform_exp(&a1,&b1,&c1,a,b,c,q1,L);  qform_exp(&u1,&v1,&w1,u,v,w,q1,L);
	e1 = qform_pp_discrete_log (a1,b1,c1,u1,v1,w1,p,h2,L);														// (alpha^q1)^e1 = beta^q1
	if ( e1 < 0 ) return e1;
	qform_exp(&u1,&v1,&w1,a,b,c,e1,L);  qform_invert(&u1,&v1,&w1); qform_compose(&u1,&v1,&w1,u1,v1,w1,u,v,w,L);		// beta*alpha^{-e1}
	for ( q2 = p, i = 1 ; i < h2 ; i++ ) q2 *=p;
	qform_exp(&a1,&b1,&c1,a,b,c,q2,L); 
	e2= qform_pp_discrete_log (a1,b1,c1,u1,v1,w1,p,h1,L);														// (alpha^q2)^e2 = beta*alpha^{-e1}
	if ( e2 < 0 ) return e2;
//printf ("pp DL computed %ld*%ld+%ld  = %ld\n", q2,e2,e1,q2*e2+e1);
	return q2*e2+e1;
}

long qform_discrete_log_ph_ppf (long a, long b, long c, long u, long v, long w, ppf_t M, long L)
{
	ppf_t M1, M2;
	long m1,m2,e1,e2,a1,b1,c1,u1,v1,w1,x,y,n;

//printf ("PH_DL <%ld,%ld,%ld> <%ld,%ld,%ld> with %d factor M\n", a,b,c,u,v,w,M->w,M);
	if ( M->w == 1 ) return qform_pp_discrete_log (a,b,c,u,v,w,M->p[0],M->h[0],L);
	ppf_split (M1, M2, M);
	m1 = ppf_eval(M1);  m2 = ppf_eval(M2);
	ui_gcd_ext (m1,m2,&x,&y);
//printf ("%ld*%ld + %ld*%ld = 1\n", x,m1,y,m2);
	qform_exp(&a1,&b1,&c1,a,b,c,m1,L);  qform_exp(&u1,&v1,&w1,u,v,w,m1,L);
	e1 = qform_discrete_log_ph_ppf (a1,b1,c1,u1,v1,w1,M2,L);
	if ( e1 < 0 ) return e1;
//printf ("e1=%ld\n", e1);
	qform_exp(&a1,&b1,&c1,a,b,c,m2,L);  qform_exp(&u1,&v1,&w1,u,v,w,m2,L);
	e2 = qform_discrete_log_ph_ppf (a1,b1,c1,u1,v1,w1,M1,L);
	if ( e2 < 0 ) return e2;
//printf ("e2=%ld\n", e2);
	n = i_mod (e1*x*m1+e2*y*m2, m1*m2);
//printf ("%ld*%ld*%ld + %ld*%ld*%ld mod %ld = %ld\n", e1, x, m1, e2, y, m2, m1*m2, n);
/*
if ( n >= 0 ) { 
	qform_exp(&u1,&v1,&w1,a,b,c,n,L);
	if ( ! qform_equal(u1,v1,w1,u,v,w) ) { err_printf ("Discrete PH log failed <%ld,%ld,%ld>^%d != <%ld,%ld,%ld>, m=%d\n", a, b, c, n, u, v, w, m1*m2); exit (0); }
} else {
	if ( qform_discrete_log(a,b,c,u,v,w,m1*m2,L) >= 0 ) { err_printf ("Discrete log PH failed <%ld,%ld,%ld>^%d = <%ld,%ld,%ld>, m=%d\n", a,b,c,n,u,v,w,m1*m2); exit (0); }
}
*/
	return n;
}


long qform_discrete_log (long a, long b, long c, long u, long v, long w, long n, long L)
{
	long bu, bv, bw, gu, gv, gw;
	register long i,j,k,bsteps,gsteps;

	if ( n < 16 ) return qform_discrete_log_naive (a,b,c,u,v,w,n,L);
	if ( qform_is_identity(u,v,w) ) return 0;
	bsteps = sqrt(n)/2;
	qform_table_init(2*bsteps,0);
	gsteps = ui_ceil_ratio(n,2*bsteps);

//printf ("BSGS start <%ld,%ld,%ld> <%ld,%ld,%ld> n=%ld\n", a,b,c,u,v,w,n);
	qform_insert_hi(a,b,c);
//printf ("baby insert 1: %ld,%ld\n", a, b);
	qform_square(&bu,&bv,&bw,a,b,c,L);
	if ( qform_is_identity(bu,bv,bw) ) { err_printf ("bad n=%ld in qform_discrete_log for <%ld,%ld,%ld> with order 1 or 2\n", n, a, b, c); exit (0); }
	qform_insert_hi(bu,bv,bw);
//printf ("baby insert 2: %ld,%ld\n", bu, bv);
	for ( i = 2 ; i < bsteps ; i++ ) {
		qform_compose(&bu,&bv,&bw,bu,bv,bw,a,b,c,L);
		qform_insert_hi(bu,bv,bw);
//printf ("baby insert %d: %ld,%ld\n", i+1, bu, bv);
	}
	qform_set (&gu,&gv,&gw,u,v,w);
	qform_invert(&gu,&gv,&gw);
	qform_square(&bu,&bv,&bw,bu,bv,bw,L);
	for ( i = 0 ; i < gsteps ; i++ ) {
//printf ("giant lookup %d: %ld,%ld\n", i, gu, gv);
		if ( (j=qform_lookup_hi(gu,gv,gw)) >= 0 ) break;
		qform_compose(&gu,&gv,&gw,gu,gv,gw,bu,bv,bw,L);
	}
//printf ("BSGS end i=%d, j=%d\n", i, j);
	if ( i < gsteps ) {
		if ( qform_table[j].b == gv ) k = ( i ? 2*i*bsteps-j : ( j ? n-j : j ) );
		else k = 2*i*bsteps+j;
	} else {
		k = -1;
	}
	qform_table_free();
/*
if ( k >= 0 ) { 
	qform_exp(&bu,&bv,&bw,a,b,c,k,L);
	if ( ! qform_equal(bu,bv,bw,u,v,w) ) { err_printf ("Discrete log failed <%ld,%ld,%ld>^%d != <%ld,%ld,%ld>, found match with i=%d, j=%d, bsteps=%d, gsteps=%d, n=%d\n", a, b, c, k, u, v, w, i, j, bsteps, gsteps, n); exit (0); }
} else {
	if ( (k=qform_discrete_log_slow(a,b,c,u,v,w,n,L)) >= 0 ) { err_printf ("Discrete log failed to find <%ld,%ld,%ld>^%d = <%ld,%ld,%ld> with bsteps=%d, gsteps=%d, n=%d\n", a,b,c,k,u,v,w,bsteps,gsteps,n); exit (0); }
}
*/
	return k;
}

// n2 is the order of alpha_p, and e1 must be congruent to DL(alpha_p,alpha_q) mod n1, where n1 | n2.
long qform_primeform_mod_discrete_log (long p, long q, long e1, long n1, long n2, long D)
{
	long a,b,c,u,v,w,x,y,z,e2,L,o;
	
	if ( ! qform_primeform (&a,&b,&c,p,D,0) ) return -1;
	qform_reduce(&a,&b,&c);
	if ( ! qform_primeform (&u,&v,&w,q,D,0) ) return -1;
	qform_reduce(&u,&v,&w);
	L = ceil(pow(-D,0.25));
	qform_exp (&x,&y,&z,a,b,c,e1,L); 	qform_invert (&x,&y,&z); 
	qform_compose(&u,&v,&w,u,v,w,x,y,z,L);					// (u,v,w) = beta*alpha^{-e1}
	qform_exp (&a,&b,&c,a,b,c,n1,L);							// (a,b,c) = alpha^n1
	o = n2/n1;
	e2 = qform_discrete_log (a,b,c,u,v,w,o,L);
//printf ("D=%ld, p=%d,q=%d, n1=%d, e1=%d, n2=%d, e2=%d\n", D, p, q, n1, e1, n2, e2);
	if ( e2 < 0 ) return -1;
	return e1+e2*n1;
}


// simple binary exponentiation, overlap ok
void qform_exp (long *pu, long *pv, long *pw, long u, long v, long w, long e, long L)
{
	long a, b, c;
	register long m;
	register int i;

	switch (e) {
	case 0: qform_set_identity (pu, pv, pw, v*v-4*u*w); return;
	case 1: qform_set (pu,pv,pw,u,v,w); return;
	case 2: qform_square (pu, pv, pw, u, v, w, L); return;
	case 3: qform_square (&a,&b,&c,u,v,w,L); qform_compose (pu, pv, pw, a, b, c, u, v, w, L); return;
	case 4: qform_square (&a,&b,&c,u,v,w,L); qform_square (pu, pv, pw, a, b, c, L); return;
	case 5: qform_square (&a,&b,&c,u,v,w,L); qform_square (&a, &b, &c, a, b, c, L); qform_compose (pu, pv, pw, a, b, c, u, v, w, L); return;
	case 6: qform_square (&a,&b,&c,u,v,w,L); qform_compose (&a, &b, &c, a, b, c, u, v, w, L); qform_square (pu, pv, pw, a, b, c, L); return;
	case 7: qform_square (&a,&b,&c,u,v,w,L); qform_compose (&a, &b, &c, a, b, c, u, v, w, L); qform_square (&a, &b, &c, a, b, c, L); qform_compose(pu, pv, pw, a, b, c, u, v, w, L); return;
	case 8: qform_square (&a,&b,&c,u,v,w,L); qform_square (&a, &b, &c, a, b, c, L); qform_square (pu, pv, pw, a, b, c, L); return;
	}
		
	a=u; b=v; c=w;
	i = ui_len(e)-1;
	m = (1L<<(i-1));
	while (m) {
		qform_square (&a, &b, &c, a, b, c, L);
		if ( (e&m) ) qform_compose (&a, &b, &c, a, b, c, u, v, w, L);
		m >>=1;
	}
	*pu=a; *pv=b; *pw=c;
	return;	
}


// finds the reduced form of discriminant D with *pa > a0 prime and minimal, with *pa <= MPZ_MAX_SMALL_PRIME, or returns 0 if no such form exists
// note that if we have to search past MPZ_MAX_TINY_PRIME, things slow down
int qform_small_primeform (long *pa, long *pb, long *pc, long D, long a0)
{
	register long p;
	register int i;

	if ( a0 >= MPZ_MAX_SMALL_PRIME ) return 0;
	if ( D > -3*a0*a0 ) return 0;												// |b|<=a<=c implies D=b^2-4ac <= -3a^2
	for ( i = ui_small_prime_index(a0) + 1 ; i <= MPZ_SMALL_PRIMES ; i++ ) {
		p = ui_small_prime(i);
		if ( D > -3*p*p ) return 0;
		if ( ! prime_to_conductor(D,p) ) continue;
		if ( qform_primeform (pa, pb, pc, p, D, 1) ) return 1;
	}
	return 0;
}


int qform_primeform (long *pa, long *pb, long *pc, long p, long D, int rflag)
{
	long b;
//printf ("computing primeform for p=%ld, D=%ld, with rflag=%d\n", p, D, rflag);
	if ( ! prime_to_conductor(D,p) ) { dbg_printf ("Request for primeform with norm %ld not prime to the conductor of D=%ld\n", p, D); return 0; }
	if ( p==2 ) {
		// need to handle p=2 separately: we want D=b^2-8c where b = 0, 1, or 2 and c >= 2.  Each case determines D mod 8 and we must have D < -8
		switch (D&7) {
		case 0: *pb = 0; *pc = -D/8; break;
		case 1: *pb = 1; *pc = (1-D)/8; break;
		case 4: *pb = 2; *pc = (4-D)/8; break;
		default: return 0;
		}
		*pa = 2;
	} else {
		if ( (b = tiny_sqrt_modprime(D,p)) < -1 ) b = i_sqrt_modprime(D,p);
		if ( b < 0 ) return 0;
		if ( D&3 ) {													// if D is 1 mod 4
			if ( !(b&1) ) b = p-b;										// make b odd
		} else {														// if D is 0 mod 4
			if ( b&1 ) b = p-b;											// make b even
		}
		*pc = (b*b - D)/(4*p);
		*pa = p;
		*pb = b;
	}
//printf ("primeform D=%ld, (%ld,%ld,%ld), rflag=%d\n", D, *pa, *pb, *pc, rflag);

	if ( *pa==*pb && !((*pc)%p) ) return 0;									// form not primitive
	if ( rflag && *pc < p ) {
		if ( rflag>0 ) return 0;											// require that form is reduced for rflag > 0
		qform_reduce (pa,pb,pc);										// reduce it if rflag < 0
	}
	return 1;
}

// computes product of primeforms over primes appearing in the prime factorization of n (including multiplicity)
int qform_nform (long *pa, long *pb, long *pc, long n, long D)
{
	unsigned long p[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS];
	long u,v,w,x,y,z, L;
	int i, j, k;
	
	L = ceil(pow(-D,0.25));
	k  = ui_factor(p,h,n);
	for ( i = 0 ; i < k ; i++ ) {
		if ( ! qform_primeform (&x,&y,&z,p[i],D,-1) ) return 0;
		for ( j = 0 ; j < h[i] ; j++ ) if ( i || j ) qform_compose (&u,&v,&w,u,v,w,x,y,z,L); else { u=x; v=y; w=z; }
	}
	*pa=u; *pb=v; *pc=w;
	return 1;
}


// computes product of primeforms over primes appearing in the prime factorization of n (including multiplicity)
int qform_ppf_form (long *pa, long *pb, long *pc, ppf_t N, long D)
{
	long u,v,w,x,y,z, L;
	int i;
	
	L = ceil(pow(-D,0.25));
	for ( i = 0 ; i < N->w ; i++ ) {
		if ( ! qform_primeform (&x,&y,&z,N->p[i],D,-1) ) { err_printf ("Couldn't get primeform for ell=%ld\n", N->p[i]); return 0; }
		qform_exp (&x,&y,&z,x,y,z,N->h[i],L);
		if ( i ) qform_compose (&u,&v,&w,u,v,w,x,y,z,L); else { u=x; v=y; w=z; }
	}
	*pa=u; *pb=v; *pc=w;
	return 1;
}


// assumes a > 0, but b may be negative
static inline long gcdext(long a, long b, long *x, long *y)
{
	register long d;
	
	if ( b < 0 ) {
		d = ui_gcd_ext ((unsigned long)a, (unsigned long)(-b), x, y);
		if ( y ) *y = -*y;
	} else {
		d = ui_gcd_ext ((unsigned long)a, (unsigned long)b, x, y);
	}
	return d;
}
		

// Algorithm 5.4.2 in Cohen  "A Course in Computational Algebraic Number Theory" p. 243
// This code may go into an infinite loop if the input is bogus (e.g. a and c must be positive)
int qform_reduce (long *pa, long *pb, long *pc)
{
	register long x, q, r, a, b, c;

#if ! QFORM_FAST
	if ( *pa <= 0 || *pc <= 0 ) { printf ("Form must be positive definite! a = %ld, c = %ld\n", *pa, *pc);  exit (0); }
#endif
	a = *pa;  b = *pb;  c = *pc;
	for (;;) {
		if ( a > c ) { _swap(a,c,x); b = -b; }
		if ( b > a || -b >= a ) {
			x = 2*a;
			q = b/x;  r = b-q*x;
			if ( r <= a ) { r += x; q--; }
			if ( r > a ) { r -= x; q++; }
			c -= ((b+r)*q)>>1;
			b = r;
		} else {
			break;
		}
	}
	if ( a==c && b < 0 ) b = -b;
#if ! QFORM_FAST
	if ( a <= 0 || c <= 0 || a >c || b > a || -b >= a ) { printf ("Error in qform_reduce, result (%ld,%ld,%ld) not reduced\n", a, b, c); exit (0); }
#endif
	*pa = a;  *pb = b;  *pc = c;
	return 1;
}

// Algorithm 1 NUCOMP from Jacobson and van der Poorten "Computational Aspects of NUCOMP" in ANTS V (2002) p. 120
// L = |D|^{1/4} is a required precomputed value, overlap is ok
void qform_compose (long *pu3, long *pv3, long *pw3, long u1, long v1, long w1, long u2, long v2, long w2, long L)
{
	register long Bx, By, Cy, Dy, F, G, H, Q1, Q2, Q3, Q4, ax, ay, bx, by, cx, cy, dx, dy, m, q, s, t, ell, z;		// we certainly don't need so many variables, but we'll let the compiler sort this out
	long b, c, x, y;

//printf("compose (%ld,%ld,%ld)*(%ld,%ld,%ld)\n", u1, v1, w1, u2, v2, w2);
//printf("L=%ld\n", L);

//step1:
	if ( w1 < w2 ) { _swap(u1,u2,x); _swap(v1,v2,x); _swap(w1,w2,x); }
	s = (v1+v2)/2;
	m = v2-s;
//step2:
	F = gcdext (u1, u2, &c, &b);
	if ( !(s%F) ) {  G = F;  /*Ax = G;*/ Bx = m*b;  By = u1/G;  Cy = u2/G;  Dy = s/G; goto step5; }	// Ax is never used
//step3:
	G = gcdext (F, s, &x, &y);
	H = F/G;  By = u1/G;  Cy = u2/G;  Dy = s/G;
//step4:
	ell = y*(b*(w1%H)+c*(w2%H)) % H;
	Bx = b*(m/H) + ell*(By/H);
step5:
	bx = Bx%By;  by = By;
//step5a:
	x = 1;  y = 0;  z = 0;
step5b:
	if ( (by > L || -by > L) && bx ) goto step5c;
	if ( z & 0x1 ) { by = -by;  y = -y; }
	ax = G*x;  ay = G*y;
	goto step6;
step5c:
	q = by/bx;  t = by - q*bx;
	by = bx;  bx = t;
	t = y - q*x;  y = x;  x = t;  z++;
	goto step5b;
step6:
	if ( z ) goto step7;
	Q1 = Cy*bx;  cx = (Q1-m)/By;
	dx = (bx*Dy-w2)/By;
	*pu3 = by*Cy;  *pw3 = bx*cx - G*dx;  *pv3 = v2 - 2*Q1;
	goto step8;
step7:
	cx = (Cy*bx-m*x)/By;
	Q1 = by*cx;  Q2 = Q1+m;
	dx = (Dy*bx-w2*x)/By;
	Q3 = y*dx;  Q4 = Q3+Dy;  dy = Q4/x;
	cy = ( bx ? Q2/bx : (cx*dy-w1)/dx );
	*pu3 = by*cy-ay*dy;  *pw3 = bx*cx-ax*dx;  *pv3 = G*(Q3+Q4) - Q1 - Q2;
step8:
	// we don't need to compute gamma, we are in the imaginary quadratic case
//printf("unreduced result (%ld,%ld,%ld)\n", *pu3, *pv3, *pw3);
	qform_reduce (pu3, pv3, pw3);
//printf("reduced result (%ld,%ld,%ld)\n", *pu3, *pv3, *pw3);
}

// Algorithm 2 NUDUPL from Jacobson and van der Poorten "Computational Aspects of NUCOMP" in ANTS V (2002) p. 123
// L = |D|^{1/4} is a required precomputed value, overlap is ok
void qform_square (long *pu3, long *pv3, long *pw3, long u, long v, long w, long L)
{
	register long G, Bx, By, Dy, Q1, ax, bx, by, dx, ay, dy, q, x, t, z;
	long y;

//printf("square (%ld,%ld,%ld)\n", u, v, w);
	
//step1:
	G = gcdext (u, v, 0, &y);
	By = u/G;  Dy = v/G;		// note Ax in Algorithm 2 is never used
//step2:
	Bx = (y*w)%By;
//step3:
	bx = Bx;  by = By;
//step3a:
	x = 1;  y = 0;  z = 0;
step3b:
	if ( (by > L || -by > L) && bx ) goto step3c;
	if ( z&0x1 ) { by=-by;  y=-y; }
	ax = G*x;  ay = G*y;
	goto step4;
step3c:
	q = by/bx;  t = by-q*bx;
	by = bx;  bx = t;
	t = y - q*x;  y = x;  x = t; z++;
	goto step3b;
step4:	
	if ( z ) goto step5;
//if ( ui_len(bx) > 30 && ui_len(-bx) > 30 ) printf("Big bx=%ld\n", bx);
//if ( ui_len(by) > 30 && ui_len(-by) > 30 ) printf("Big bx=%ld\n", bx);
	dx = (bx*Dy-w)/By;
	*pu3 = by*by;  *pw3 = bx*bx - G*dx;  *pv3 = v - 2*bx*by;					// don't bother with optimization to trade mults for squares
	goto step6;
step5:
//if ( ui_len(bx) > 30 && ui_len(-bx) > 30 ) printf("Big bx=%ld\n", bx);
//if ( ui_len(by) > 30 && ui_len(-by) > 30 ) printf("Big bx=%ld\n", bx);
	dx = (bx*Dy-w*x)/By;
	Q1 = dx*y;  dy = Q1+Dy;
	*pu3 = by*by - ay*(dy/x);  *pw3 = bx*bx - ax*dx;  *pv3 = G*(dy+Q1)-2*bx*by;		// ditto
step6:
	qform_reduce (pu3, pv3, pw3);
}
