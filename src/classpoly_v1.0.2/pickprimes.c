#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include "iqclass.h"
#include "mpzutil.h"
#include "torcosts.h"
#include "tecurve.h"
#include "findcurve.h"
#include "pickprimes.h"
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

#define MAX_V		 	1200		// this really should be consistent with PHI_MAX_V
#define MAX_VK		5

static int v_primorials[MAX_VK+1] = { 1, 2, 6, 30, 210, 2310 };

static int32_t *sieve_primes;
static int32_t *sieve_sqrts;
static int sieve_pcnt, sieve_pmax;
static long minL = 3;


double split_prime_rating (long p, long t, int *ptwist, int *ptor, int *ps2, int *pt3);

void setup_sieve_primes (long maxL, long D)
{
	prime_enum_ctx_t *ctx;
	double x;
	int n;
	int32_t a;
	register unsigned p;
	
	if ( maxL <= minL ) { err_printf ("Invalid maxL=%ld in setup_sieve_primes, minL=%ld\n", maxL, minL); exit (0); }
	x = maxL;  x = x/(log(x)) * (1.0 + 3.0/(2.0*log(x)));						// pi(x) < x/log(x) * (1+3/(2log(x))) for x >= 59 (Shoup p.91)
	n = (int)x;  if ( n < 17 ) n = 17;
	if ( ! sieve_primes ) {
		sieve_pmax = n/2+2*sqrt(n)+100;								// assume about half won't yield a square root
		sieve_primes = malloc(sieve_pmax*sizeof(*sieve_primes));
		sieve_sqrts = malloc(sieve_pmax*sizeof(*sieve_sqrts));
		dbg_printf ("sieve_primes malloc %ld bytes\n", sieve_pmax*(sizeof(*sieve_primes)+sizeof(*sieve_sqrts)));
	}
	bitmap_alloc(maxL);
	ctx = fast_prime_enum_start (minL,maxL,0);
	for ( ; (p = fast_prime_enum(ctx)) ; ) {
		if ( sieve_pcnt > n ) { err_printf ("Fatal error, pi(x) upper bound failed %d > %d!\n", sieve_pcnt, n); exit (0); }
		a = (int32_t)i_sqrt_modprime (D,p);
		if ( a < 0 && p > MAX_V ) continue;
		if ( sieve_pcnt >= sieve_pmax ) {
			sieve_pmax = (3*sieve_pmax)/2;
			sieve_primes = realloc(sieve_primes, sieve_pmax*sizeof(*sieve_primes));
			sieve_sqrts = realloc(sieve_sqrts, sieve_pmax*sizeof(*sieve_sqrts));
			dbg_printf ("sieve_primes realloc %ld bytes\n", sieve_pmax*(sizeof(*sieve_primes)+sizeof(*sieve_sqrts)));
		}
		sieve_primes[sieve_pcnt] = (int32_t)p;
		sieve_sqrts[sieve_pcnt] = a;
		sieve_pcnt++;
	}
	fast_prime_enum_end(ctx);
	minL = maxL;
	dbg_printf ("sieve prime count %d, max=%d, maxL=%ld\n", sieve_pcnt, sieve_primes[sieve_pcnt-1], maxL);
}

void free_sieve_primes (void)
{
	if ( sieve_primes ) {
		info_printf ("sieve_primes used %ld bytes\n", sieve_pmax*(sizeof(*sieve_primes)+sizeof(*sieve_sqrts)));
		free(sieve_primes); free(sieve_sqrts); sieve_primes = 0; sieve_sqrts = 0; sieve_pcnt = sieve_pmax = 0;  minL = 3;
	}
}


// This is ridiculously slow and is implemented only to prove this fact, and as a sanity check
int pick_primes_cornacchia (long *plist, int *tlist, int *vlist, int maxn, long D, long h, double bits)
{
	static mpz_t p, t, v, d;
	static int init;
	register double pbits;
	register int n;
	
	if ( ! init ) { mpz_init(p); mpz_init(t); mpz_init(v); mpz_init(d); init = 1; }
	mpz_set_ui(d,-D);
	mpz_set_ui(p,-D/4-1);
	
	for ( n = 0, pbits = 0.0 ; pbits < bits ; ) {
		if ( n >= maxn ) { err_printf ("maxn=%d is too small in pick_primes, D=%ld, bits=%f!\n", maxn, D, bits); exit (0); }
		mpz_nextprime(p,p);
		if ( mpz_cornacchia4 (t, v, d, p) ) {
			if ( ! mpz_sgn(t) ) continue;
			plist[n] = mpz_get_ui(p);  tlist[n] = mpz_get_ui(t);  vlist[n] = mpz_get_ui(v);
			pbits += log2(plist[n]);
			n++;
		}
	}
	return n;
}

struct split_prime_struct {
	long p;
	int32_t t,v;
	double r;
// we could save these values for use by findcurve, but currently we don't bother
	char tor;
	char twist;
	char t3_flag;
	char s2_flag;
};

int split_prime_cmp (const void *a, const void *b)
{
	struct split_prime_struct *p1, *p2;
	
	p1 = (struct split_prime_struct *)a;  p2 = (struct split_prime_struct *)b;
	if ( p1->r > p2->r ) return 1;
	if ( p1->r < p2->r ) return -1;
	if ( p1->v > p2->v ) return 1;		// should use smaller p, but for consistency use v as tiebreaker
	if ( p1->v < p2->v ) return -1;
	return 0;
}

#define SZ_FACTOR		2

int pick_primes (long **pplist, int **ptlist, int **pvlist, long D, long h, double bits, int pfilter, long ellfilter)
{
	struct split_prime_struct *split_primes;
	time_t start, end;
	long split_pcnt;
	double pbits,r,Dbits;
	long t[MAX_V+1], m[MAX_V+1], wts[MAX_V+1], maxL, a0;
	double vxb[MAX_VK+1];
	int torcnts[40];
	int twist,tor,s2,t3;
	long maxv, maxp;
	register long a, L, C, N, p, q;
	register long i,n,v;
	
	// compute bounds per lemma in the appendix of "Computing Hilbert Class Polynomials with the CRT"
	for ( vxb[0] = 1.0, i = 1 ; i <= MAX_VK ; i++ ) vxb[i] = vxb[i-1]*(1.0+2.0/ui_small_prime(i));	
	
	/*
		We assume here that D is fundamental, but this code works reasonably well even when this is not true.
	*/
	
	// Compute wts[v] = H(-v^2D) = sum_{u|v} h(v^2D) = sum_{u|v} h * prod_{p|u}(1-(D/p))  where (D/p) is the Kronecker symbol
	for ( v = 1 ; v <= MAX_V ; v++ ) wts[v] = h;
	// Here we take advantage of the fact that H(-v^2D) is multiplicative
	for ( p = 2 ; p <= MAX_V ; p = ui_next_prime(p) ) {
		n = p-kronecker_p(D,p);
		for ( q = p ; q <= MAX_V ; q *= p ) {
			for ( v = q ; v <= MAX_V ; v += q ) {
				if ( ! (v%(p*q)) ) continue;
				wts[v] *= (q -1)/(p-1)*n+1;
			}
		}
	}
	// Adjust for unit count in cl(D) when D_0=-3 or -4.  Round up to get correct count of #curves with a given trace (note that H(3*v^2) and H(4*v^2) are not integers).
	// Here we do need need to handle the case that D is not fundamental to get the correct value of D_0
	v = discriminant_test(D);
	if ( ! v ) { err_printf ("Invalid disciminant %ld in pick_primes_sieve!\n", D); exit (0); }
	if ( D/(v*v) == -3 ) {
		for ( i = 1 ; i <= MAX_V ; i++ ) wts[i] = ui_ceil_ratio(wts[i],3);
	} else if ( D/(v*v) == -4 ) {
		for ( i = 1 ; i <= MAX_V ; i++ ) wts[i] = ui_ceil_ratio(wts[i],2);		
	}
	
	maxL = (long)sqrt(-D*log(-D)) + 1;
	if ( maxL < 100 ) maxL = 97;
	setup_sieve_primes(maxL,D);
	bitmap_alloc(maxL);
	split_pcnt = (long) (SZ_FACTOR*bits/log2(-D) + 1);
	split_primes = malloc(split_pcnt*sizeof(*split_primes));
	dbg_printf ("split_pcnt=%ld, split_primes malloced %ld bytes\n", split_pcnt, split_pcnt*sizeof(*split_primes));
	pbits = 0.0;  n = 0;
	C = SZ_FACTOR*(-D/h)/4+2;											// set initial curve count bound (twice the min for the least prime)
	Dbits = log2(-D);
	for ( v = 1 ; v < MAX_V ; v++ ) {
		if ( Dbits + 2*log2(v) > 60 ) break;
		switch ( (v*v*D)&0x7 ) {
		case 0: if ( (pfilter&IQ_FILTER_1MOD4) ) { t[v]  = 0; } else { t[v] = 2; m[v] = 4; } break;
		case 5: t[v] = 1; m[v] = 2; break;
		case 4: t[v] = 4; m[v] = 4; break;
		default: t[v] = 0;
		}
		if ( ellfilter>1 && ui_gcd(v,ellfilter) > 1 ) t[v]=0;
		if ( (pfilter&IQ_FILTER_1MOD3) && !(v%3) ) t[v] = 0;
	}
	while ( v < MAX_V ) t[v++] = 0;
	maxv=maxp=0;
	start = clock();
	for(;;) {
		if ( ui_len(C) > 60 ) { err_printf("Found %ld primes with %f bits up to C=%ld for D=%ld, giving up\n", n, pbits, C, D); exit(0); }
		for ( v = 1 ; v <= MAX_V ; v++ ) {
			if ( ! t[v] ) continue;
			for ( i = 1 ; v_primorials[i] <= v; i++ );
			i--;
			if ( (-v*D)/(4*vxb[i]*h) > C ) continue;
			if ( v==MAX_V ) { if ( pbits < bits ) { err_printf ("Need to increase MAX_V, pfilter=%d, v=%ld, D=%ld, C=%ld, vxb=%f\n", pfilter, v, D, C, vxb[i]); exit (0); } else { break; } }
			if ( ui_len(C)+ui_len(wts[v]) > 60 ) continue;				// avoid overflow here
			N = C*wts[v];
			if ( 2*ui_len(v)+ui_len(-D) > 62 ) continue;				// watch out for overflow here
			if ( -v*v*D >= 4*N ) continue;
			L = sqrt(4*N)+1;
			while( L >= sieve_primes[sieve_pcnt-1] ) { maxL=(3*maxL)/2; setup_sieve_primes (maxL,D); }
			for ( i = 0 ; ; i++ ) {
				p = sieve_primes[i];
				if ( p > L ) break;
				a = sieve_sqrts[i];
				if (! (v%p) ) a = 0;
				if ( a < 0 ) continue;
				a = ui_mod(v*a,p);
				a0 = a;
				if ( a < t[v] ) a += p*((t[v]-a)/p);
				if ( a*a-v*v*D != p && a < L ) bitmap_set(a);
				for ( a += p ; a < L ; a+= p ) bitmap_set(a);
				if ( ! a0 ) continue;
				a = p-a0;
				if ( a < t[v] ) a += p*((t[v]-a)/p);
				if ( a*a-v*v*D != p && a < L ) bitmap_set(a);
				for ( a += p ; a < L ; a += p ) bitmap_set(a);
			}
			/*
				we could speed up the loop below by using, say, a 48 wheel that takes into account the mod4 and mod3 filters,
				but it probably won't make much difference.
			*/
			for ( ; t[v] < L ; t[v] += m[v] ) {
				p = (t[v]*t[v]-v*v*D)>>2;
				if ( p < MIN_CRT_PRIME ) continue;
				if ( p > N ) break;
				if ( (pfilter&IQ_FILTER_1MOD4) && (p&0x3)==1 ) continue;
				if ( (pfilter&IQ_FILTER_1MOD3) && (p%3)==1 ) continue;
				if ( bitmap_test(t[v]) ) continue;
//				if ( !(p&1) ) continue;
//				if ( ui_gcd(614889782588491410UL, p) > 1 ) continue;
//		if ( ! ui_is_prime (p) ) { printf ("Sieve error, p=%ld for t=%d, D=%ld is not prime!\n", p, t[v], D); exit (0); }
//		if ( 4*p != t[v]*t[v] - v*v*D ) { printf ("Sieve error, 4*%ld != %d^2 - %ld \n", p, t[v], D); exit (0); }
				split_primes[n].p = p;
				split_primes[n].t = t[v];
				split_primes[n].v = v;
				r = split_prime_rating_new(p,t[v],&twist,&tor,&s2,&t3);
				split_primes[n].tor = tor;
				split_primes[n].twist = twist;
				split_primes[n].s2_flag = s2;
				split_primes[n].t3_flag = t3;
				split_primes[n].r = (r*p)/wts[v];
/*
{
	ff_t f[4];
	long cnt,totcnt,msecs;
	int j;
puts("");
printf ("Relative rating for p=%ld t=%ld is r=%f, using tor=%d, s2=%x, t3=%d, twist=%d\n", p, t[v], r, tor, s2, t3, twist);
printf ("r/rho = %.0f v=%d, wts[v]=%d\n", split_primes[n].r, v, wts[v]);
for ( j = 0;j < 3 ; j++ ) {
ff_setup_ui(p);
start = clock();
totcnt = 0;
for ( i = 0 ; i < 10 ; i++ ) {
	findcurve (f, t[v], FINDCURVE_NOFILTER, &cnt);  totcnt += cnt;
}
end = clock();
printf ("NO_FILTER: ave cnt = %7ld, ave time %5.1f msecs\n", totcnt/10, (double)delta_msecs(start,end)/10.0);
msecs = delta_msecs(start,end);
start = clock();
totcnt = 0;
for ( i = 0 ; i < 10 ; i++ ) {
	findcurve (f, t[v], 0, &cnt);  totcnt += cnt;
}
end = clock();
printf ("   FILTER: ave cnt = %7ld, ave time %5.1f msecs, ratio %f (vs %f)\n", totcnt/10,  (double)delta_msecs(start,end)/10.0, (double)delta_msecs(start,end)/msecs, r);
}
}
*/				pbits += log2(p);
				if ( p > maxp ) maxp=p;
				if ( v > maxv ) maxv = v;
				n++;
				if ( n >= split_pcnt ) { split_pcnt = (3*split_pcnt)/2;  dbg_printf ("split_primes reallocing %ld bytes\n", split_pcnt*sizeof(*split_primes));  split_primes = realloc(split_primes,split_pcnt*sizeof(*split_primes)); dbg_printf ("split_primes realloc %ld bytes\n", split_pcnt*sizeof(*split_primes)); }
			}
			bitmap_clear(L);
			if ( n >= split_pcnt ) break;
		}
		if ( pbits >= SZ_FACTOR*bits ) break;
		C = (3*C)/2;
	}
	free_sieve_primes();
	bitmap_free();
	end = clock();
	info_printf ("sieved %ld split primes up to z=%ld obtaining a total of %f > %f bits, maxp=%ld, maxv=%ld %ld msecs\n", n, C, pbits, bits, maxp, maxv, delta_msecs(start,end));
	qsort(split_primes,n,sizeof(*split_primes),split_prime_cmp);
	for ( i = 0, pbits = 0.0 ; i < n && pbits < bits ; i++ ) pbits += log2(split_primes[i].p);
	if ( pbits < bits ) { printf ("error, didn't get the expected number of bits\n"); exit (0); }
	n = i;
	*pplist = malloc(n*sizeof(**pplist));
	if ( ptlist ) *ptlist = malloc (n*sizeof(**ptlist));
	if ( pvlist ) *pvlist = malloc(n*sizeof(**pvlist));
	for ( i = 0 ; i < 40 ; i++ ) torcnts[i] = 0;
	for ( i = 0 ; i < n ; i++ ) {
		(*pplist)[i] = split_primes[i].p;
		if ( ptlist ) (*ptlist)[i] = split_primes[i].t;
		if ( pvlist ) (*pvlist)[i] = split_primes[i].v;
		torcnts[(int)split_primes[i].tor]++;

//split_prime_rating_new(split_primes[i].p,split_primes[i].t,&twist,&tor,&s2,&t3);
//out_printf("%d) %ld, t=%d, v=%d, r=%.3f, 1/rho=%ld ", i+1, split_primes[i].p,split_primes[i].t, split_primes[i].v, split_primes[i].r, split_primes[i].p/wts[split_primes[i].v]);
//out_printf(" (twist=%d, tor=%d, s2=%x, t3=%d)\n", twist, tor, s2, t3);

	}
	mem_free(split_primes);
	if ( dbg_level > 0 ) { info_printf("Torsion:\n");   for ( i = 0 ; i < 40 ; i++ ) if ( torcnts[i] ) info_printf ("    %ld (%.3f)", i, (double)torcnts[i]/n); info_printf("\n"); }
//dbg_printf ("%d primes, %f bits, maxp = %ld, maxv = %ld\n", i, pbits, p, v);
	return n;
}

#define TWIST_DOUBLE_RATIO		(9.0/16.0)

double split_prime_rating (long p, long t, int *ptwist, int *ptor, int *ps2, int *pt3)
{
	double minr, r;
	int j_flags;
	long n1, n2;
	int m1, m2, k1, k2, t1, t2, fix;
	char tormap[TECURVE_MAX_N+1];
	register int i, j, k;
	
	n1 = p+1-t;
	n2 = p+1+t;
	if ( (p%3)==1 && (4*p-t*t)%3 ) j_flags = FILTER_JCUBE; else j_flags = 0;
	t1 = ( (n1%3) ? FILTER_3TOR_1 : FILTER_3TOR_NOT_1);
	t2 = ( (n2%3) ? FILTER_3TOR_1 : FILTER_3TOR_NOT_1);
	for ( m1 = n1, k1 = 0 ; !(m1&1) ; m1>>=1, k1++ );
	for ( m2 = n2, k2 = 0 ; !(m2&1) ; m2>>=1, k2++ );
	tormap[1] = ( (n1&1) ? 3 : 0 );
	for ( i = 2 ; i <= TECURVE_MAX_N ; i++ ) {
		tormap[i] = 0;
		if ( !(n1%i) ) tormap[i] |= 1;
		if ( !(n2%i) ) tormap[i] |= 2;
	}
	minr = TECURVE_MAX_COST;
	fix = 0;
	for ( j = 0 ; j < 2 ; j++ ) {
		for ( i = 1 ; i <= TECURVE_MAX_N ; i++ ) {
			if ( tormap[i]&1 ) {
				for ( k = 0 ; k <= k1 && k <= 2 ; k++ ) {
					r = tecurve_free2_density(p,i,k,j*t1,j_flags)*tecurve_free2_cost(p,i,k,j*t1,j_flags);
//printf("Free rating tor=%d, k=%d, r=%f\n",i,k,r);
					if ( r && r < minr ) { minr=r;  *ptwist=1; fix=0; *ptor=i; *ps2=k; *pt3=j*t1; }
				}
				r = tecurve_fix2_density(p,i,k1,j*t1,j_flags)*tecurve_fix2_cost(p,i,k1,j*t1,j_flags);
//printf("Fix rating tor=%d, k=%d, r=%f\n",i,k1,r);
				if ( r && r < minr ) { minr=r;  *ptwist=1; fix=1; *ptor=i; *ps2=k1; *pt3=j*t1; }
			}
			if ( tormap[i]&2 ) {
				for ( k = 0 ; k <= k2 && k <= 2 ; k++ ) {
					r = tecurve_free2_density(p,i,k,j*t2,j_flags)*tecurve_free2_cost(p,i,k,j*t2,j_flags);
//printf("Free rating tor=%d, k=%d, r=%f\n",i,k,r);
					if ( r && r < minr ) { minr=r;  *ptwist=2; fix=0; *ptor=i; *ps2=k; *pt3=j*t2; }
				}
				r = tecurve_fix2_density(p,i,k2,j*t2,j_flags)*tecurve_fix2_cost(p,i,k2,j*t2,j_flags);
//printf("Fix rating tor=%d, k=%d, r=%f\n",i,k1,r);
				if ( r < minr ) { minr=r;  *ptwist=2; fix=1; *ptor=i; *ps2=k2; *pt3=j*t2; }
			}
			if ( tormap[i]==3 && j*t1==j*t2 ) {
				for ( k = 0 ; k <= ui_min(k1,k2) && k <= 2 ; k++ ) {
					r = tecurve_free2_density(p,i,k,j*t1,j_flags)*tecurve_free2_cost(p,i,k,j*t1,j_flags)*TWIST_DOUBLE_RATIO;
//printf("*Free rating tor=%d, k=%d, r=%f\n",i,k,r);
					if ( r && r < minr ) {minr=r;  *ptwist=3; fix=0; *ptor=i; *ps2=k; *pt3=j*t1; }
				}
				if ( k1==k2 ) {
					r = tecurve_fix2_density(p,i,k1,j*t1,j_flags)*tecurve_fix2_cost(p,i,k1,j*t1,j_flags)*TWIST_DOUBLE_RATIO;
//printf("Fix rating tor=%d, k=%d, r=%f\n",i,k1,r);
					if ( r && r < minr ) { minr=r;  *ptwist=3; fix=1; *ptor=i; *ps2=k1; *pt3=j*t1; }
				}
			}
		}
	}
	if ( fix ) {
		*ps2 |= FILTER_2SYLOW_C_LG;
	} else {
		if ( *ps2 ==1 ) *ps2 = FILTER_2SYLOW_NOT_1;
		else if ( *ps2>=2 ) *ps2 = FILTER_2SYLOW_D4;
		else *ps2 = 0;
	}

//printf ("p=%ld, tor = %d, r=%f\n", p, *ptor, minr);
	return minr;
}


double split_prime_rating_new (long p, long t, int *ptwist, int *ptor, int *ps2, int *pt3)
{
	struct torctab_rec *torctab;
	double tormod[32];
	int TORCTAB_SIZE;
	long n1, n2;
	int b1, b2, b12;
	int i;

	if ( (p%3)==2 ) {
		if ( (p&3)==3 ) {
			torctab = torctab1;
			TORCTAB_SIZE = TORCTAB1_SIZE;
		} else {
			torctab = torctab2;
			TORCTAB_SIZE = TORCTAB2_SIZE;
		}
	} else {
		if ( (p&3)==3 ) {
			torctab = torctab3;
			TORCTAB_SIZE = TORCTAB3_SIZE;
		} else {
			torctab = torctab4;
			TORCTAB_SIZE = TORCTAB4_SIZE;
		}
	}
	for ( i = 0 ; i < 32 ; i++ ) tormod[i] = 1.0;
	if ( (p%5)==1 ) tormod[5] = tormod[10] = tormod[15] = 6.0/5.0;
	if ( (p%7)==1 ) tormod[7] = tormod[14] = 8.0/7.0;
	if ( (p%11)== 1 ) tormod[11] = 12.0/11.0;
	if ( (p%13)==1 ) tormod[13] = 14.0/13.0;
	if ( (p%17)==1 ) tormod[17] = 18.0/17.0;
	if ( (p%19)==1 ) tormod[19] = 20.0/19.0;
	if ( (p%23)==1 ) tormod[23] = 24.0/23.0;
	if ( (p%29)==1 ) tormod[29] = 30.0/29.0;
	if ( (p%31)==1 ) tormod[31] = 32.0/31.0;
	
	n1 = p+1-t;
	n2 = p+1+t;
	b1 = b2 = b12 = -1;
	for ( i = 0 ; i < TORCTAB_SIZE ; i++ )
		if ( ! (n1%torctab[i].m) && ( !torctab[i].fix2 || (n1%(2*torctab[i].m)) ) && ( ! torctab[i].fix3 || (n1%(3*torctab[i].m)) ) )
			if ( b1 < 0 || torctab[i].rating*tormod[torctab[i].N] < torctab[b1].rating*tormod[torctab[b1].N] ) b1= i;
	if ( b1 < 0) { printf ("Error, no suitable torsion found in torctab for n1 = %ld!\n", n1); exit (0); }
	for ( i = 0 ; i < TORCTAB_SIZE ; i++ )
		if ( ! (n2%torctab[i].m) && ( !torctab[i].fix2 || (n2%(2*torctab[i].m)) ) && ( ! torctab[i].fix3 || (n2%(3*torctab[i].m)) ) )
			if ( b2 < 0 || torctab[i].rating*tormod[torctab[i].N] < torctab[b2].rating*tormod[torctab[b2].N]  ) b2= i;
	if ( b2 < 0 ) { printf ("Error, no suitable torsion found in torctab for n2 = %ld!\n", n2); exit (0); }
	for ( i = 0 ; i < TORCTAB_SIZE ; i++ )
		if ( ! (n1%torctab[i].m) && ( !torctab[i].fix2 || (n1%(2*torctab[i].m)) ) && ( ! torctab[i].fix3 || (n1%(3*torctab[i].m)) ) )
			if ( ! (n2%torctab[i].m) && ( !torctab[i].fix2 || (n2%(2*torctab[i].m)) ) && ( ! torctab[i].fix3 || (n2%(3*torctab[i].m)) ) )
				if ( b12 < 0 || torctab[i].rating*tormod[torctab[i].N] < torctab[b12].rating*tormod[torctab[b12].N]  ) b12 = i;
	if ( b12 < 0 ) { printf ("Error, no suitable torsion found in torctab for n1 = %ld, n2 = %ld!\n", n1, n2); exit (0); }
	if ( b1 > b2 ) {
		if ( torctab[b2].rating / TWIST_DOUBLE_RATIO > torctab[b12].rating ) {
			*ptwist = 3; *ptor = torctab[b12].N; *ps2 = torctab[b12].s2_flag; *pt3 = torctab[b12].t3_flag;  return torctab[b12].rating*tormod[torctab[b12].N] ;
		} else {
			*ptwist = 2; *ptor = torctab[b2].N; *ps2 = torctab[b2].s2_flag; *pt3 = torctab[b2].t3_flag;  return torctab[b2].rating*tormod[torctab[b2].N] /TWIST_DOUBLE_RATIO;
		}
	} else {
		if ( torctab[b1].rating / TWIST_DOUBLE_RATIO > torctab[b12].rating ) {
			*ptwist = 3; *ptor = torctab[b12].N; *ps2 = torctab[b12].s2_flag; *pt3 = torctab[b12].t3_flag;  return torctab[b12].rating*tormod[torctab[b12].N] ;
		} else {
			*ptwist = 1; *ptor = torctab[b1].N; *ps2 = torctab[b1].s2_flag; *pt3 = torctab[b1].t3_flag;  return torctab[b1].rating*tormod[torctab[b1].N] /TWIST_DOUBLE_RATIO;
		}
	}
}

