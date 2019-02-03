#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "zn_poly/zn_poly.h"
#include "zn_poly/zn_poly_internal.h"
#include "ff.h"
#include "ffpoly.h"
#include "ffpolysmall.h"
#include "ffpolybig.h"
#include "cstd.h"

/*
    Copyright 2008-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

int FFPI_MAX_BASE_N = 96;			// crossover tested on an AMD Phenom II 3.0 GHz with FF_BIG_PRIME=0 and zn_poly version 0.9 (11/3/2009)

void _fft_poly_conv_mult (ff_t h[], ff_t f[], ff_t g[], ff_t w[], ff_t wi[], ff_t half[], int n, int k);

#define _max(a,b)  ((a)>(b)?(a):(b))

static zn_mod_t _ff_poly_zn_mod_ctx;
static unsigned long _ff_poly_zn_mod_p;
static int _ff_poly_zn_mod_init;
static ff_t *_ff_poly_big_a;
static ff_t *_ff_poly_big_b;
static int _ff_poly_alloc_degree;

void ff_poly_zn_mod_setup (unsigned long p, int d)
{
	if ( p != _ff_poly_zn_mod_p ) {
		if ( _ff_poly_zn_mod_init ) zn_mod_clear(_ff_poly_zn_mod_ctx);
		zn_mod_init(_ff_poly_zn_mod_ctx, _ff_p);
		_ff_poly_zn_mod_p = _ff_p;
		_ff_poly_zn_mod_init = 1;
	}
	if ( d > _ff_poly_alloc_degree ) {
		_ff_poly_alloc_degree = (1<<ui_len(d));
		_ff_poly_big_a = realloc(_ff_poly_big_a,_ff_poly_alloc_degree*sizeof(ff_t));
		_ff_poly_big_b = realloc(_ff_poly_big_b,_ff_poly_alloc_degree*sizeof(ff_t));
	}
}

struct product_tree_struct {
	int *nodes;
	int *level[FF_PRODUCT_TREE_MAX_LEVELS];
	int level_mind[FF_PRODUCT_TREE_MAX_LEVELS];
	int h;
	int d;
};
static struct product_tree_struct trees[FF_PRODUCT_TREE_CACHE_SIZE];
static int next_tree;

void ff_poly_build_tree (struct product_tree_struct *tree, int deg)
{
	register int i, j;

	if ( deg > FF_PRODUCT_TREE_MAX_DEGREE ) { err_printf ("d=%d exceeds FF_PRODUCT_TREE_MAX_DEGREE=%d in ff_poly_build_tree\n", deg, FF_PRODUCT_TREE_MAX_DEGREE); abort(); }
	if ( ! tree->nodes ) {
		tree->nodes = malloc ((1<<FF_PRODUCT_TREE_MAX_LEVELS)*sizeof(*tree->nodes));
		for ( tree->level[0] = tree->nodes, i = 0 ; tree->level[i] + (1<<(i+2)) + (1<<(i+3)) <= tree->nodes+(1<<FF_PRODUCT_TREE_MAX_LEVELS) ;  i++ ) tree->level[i+1] = tree->level[i] + (1<<(i+2));
	}
	tree->level[0][0] = (deg+1)/2;  tree->level[0][1] = deg - tree->level[0][0];
	tree->level_mind[0] = tree->level[0][1];
	for ( i = 1 ; tree->level[i-1][0] > FF_PRODUCT_TREE_BASE_CUTOFF ; i++ ) {
		tree->level_mind[i] = deg;
		for ( j = 0 ; j < (1<<i) ; j++ ) {
			tree->level[i][2*j] = (tree->level[i-1][j]+1)/2;  tree->level[i][2*j+1] = tree->level[i-1][j] - tree->level[i][2*j];
			if ( tree->level[i][2*j+1] < tree->level_mind[i] ) tree->level_mind[i] = tree->level[i][2*j+1];
		}
	}
	tree->h = i;
	tree->d = deg;
	info_printf ("Built product tree for d=%d, h=%d, d1=%d\n", tree->d, tree->h, tree->level[tree->h-1][0]);
}


static inline struct product_tree_struct *ff_poly_get_product_tree (int d)
{
	struct product_tree_struct *tree;
	register int i;
	
	for ( i = 0 ; i < FF_PRODUCT_TREE_CACHE_SIZE ; i++ )  if ( trees[i].d == d ) return trees+i;
	tree = trees+next_tree++;
	if ( next_tree == FF_PRODUCT_TREE_CACHE_SIZE ) next_tree = 0;
	ff_poly_build_tree (tree, d);
	return tree;
}

void ff_poly_from_roots_big (ff_t f[], ff_t r[], int d, ff_t t[])
{
	struct product_tree_struct *tree;
	ff_t *s;
	register ff_t t0, t1;
	register unsigned long fudge;
	register int i, j, k, m, n, d1, d2;

	if ( d <= FF_PRODUCT_TREE_BASE_CUTOFF ) { ff_poly_from_roots_small (f, r, d); return; }
	if ( r == f ) { err_printf ("f and r cannot overlap in ff_poly_from_roots_big\n"); abort(); }
	tree = ff_poly_get_product_tree (d);
	
	s = t + d +1;
	i = tree->h-1;
	for ( j = 0, k = 0 ; j < (1<<(i+1)) ; k += n, j++ )  {
		n = tree->level[i][j];
		ff_poly_from_roots_small (s+k+j,r+k,n);
//printf ("(%d) ", n); ff_poly_print(s+k+j,n);
	}
	for ( m = 0 ; m < d+j ; m++ ) s[m]=_ff_get_ui(s[m]);										// convert out of Montgomery rep for calls to zn_poly
	ff_poly_zn_mod_setup(_ff_p, d+1);
	_ff_set_one(t1);
	for ( ; i >= 0 ; i-- ) {
		ff_square(t1,t1);
		n = tree->level_mind[i];
		fudge = _zn_array_mul_fudge(n+2,n+1, 0, _ff_poly_zn_mod_ctx);
		if ( fudge == _zn_array_mul_fudge(n+2, n+2, 0, _ff_poly_zn_mod_ctx) && fudge == _zn_array_mul_fudge(n+1,n+1,0,_ff_poly_zn_mod_ctx) ) {
			_ff_set_ui(t0,fudge);  ff_mult(t1,t1,t0);
			for ( j = 0, k=0 ; j < (1<<i) ; k += d1+d2, j++ ) {
				d1 = tree->level[i][2*j];  d2 = tree->level[i][2*j+1];
				_zn_array_mul (t, s+k+2*j, d1+1, s+k+2*j+d1+1, d2+1, 1, _ff_poly_zn_mod_ctx);
				memcpy (s+k+j,t,sizeof(*t)*(d1+d2+1));
//printf ("(%d,%d) ", d1, d2);  for ( m = d1+d2 ; m ; m-- ) printf ("%d*X^%d + ", s[k+j+m], m); printf ("%d\n", s[k+j]);
//printf ("fudge factor is %ld\n", _ff_get_ui(t1));
			}
		} else {
			// with the current zn_poly implementation this never seems to happen, but in theory it could
			info_printf ("inconsistent fudge factor for n=%d\n", n);
			for ( j = 0, k=0 ; j < (1<<i) ; k += d1+d2, j++ ) {
				d1 = tree->level[i][2*j];  d2 = tree->level[i][2*j+1];
				zn_array_mul (t, s+k+2*j, d1+1, s+k+2*j+d1+1, d2+1, _ff_poly_zn_mod_ctx);	
				memcpy (s+k+j,t,sizeof(*t)*(d1+d2+1));
//printf ("(%d,%d) ", d1, d2);  for ( m = d1+d2 ; m ; m-- ) printf ("%d*X^%d + ", f[k+j+m], m); printf ("%d\n", f[k+j]);
//printf ("fudge factor is %ld\n", _ff_get_ui(t1));
			}
		}
	}
	for ( i = 0 ; i < d+1 ; i++ )  { _ff_set_ui (t0,s[i]); _ff_mult(f[i],t0,t1); }							// back into Montgomery rep and adjust for fudge factor
	
}

void ff_poly_square_big (ff_t h[], ff_t f[], int d)
{
	register int i;

	ff_poly_zn_mod_setup(_ff_p, d);
	for ( i = 0 ; i <= d ; i++ ) _ff_poly_big_a[i] = _ff_get_ui(f[i]);
	zn_array_mul (h, _ff_poly_big_a, d+1, _ff_poly_big_a, d+1, _ff_poly_zn_mod_ctx);
	for ( i = 0 ; i <= 2*d ; i++ ) _ff_set_ui(h[i],h[i]);
}

void ff_poly_mult_big (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b)
{
	register int i;

	ff_poly_zn_mod_setup(_ff_p, _max(d_a,d_b));
	// this conversion is wasteful
	for ( i = 0 ; i <= d_a ; i++ ) _ff_poly_big_a[i] = _ff_get_ui(a[i]);
	for ( i = 0 ; i <= d_b ; i++ ) _ff_poly_big_b[i] = _ff_get_ui(b[i]);
	if ( d_a < d_b ) {
		zn_array_mul (c, _ff_poly_big_b, d_b+1, _ff_poly_big_a, d_a+1, _ff_poly_zn_mod_ctx);
	} else {
		zn_array_mul (c, _ff_poly_big_a, d_a+1, _ff_poly_big_b, d_b+1, _ff_poly_zn_mod_ctx);
	}
	for ( i = 0 ; i <= d_a+d_b ; i++ ) _ff_set_ui(c[i],c[i]);
}

void ff_poly_mult_bigx (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b)
{
	ff_poly_zn_mod_setup(_ff_p,0);
	if ( d_a < d_b ) {
		zn_array_mul (c, b, d_b+1, a, d_a+1, _ff_poly_zn_mod_ctx);
	} else {
		zn_array_mul (c, a, d_a+1, b, d_b+1, _ff_poly_zn_mod_ctx);
	}
}


// no conversion - used for performance tests
void ff_poly_square_bigx (ff_t h[], ff_t f[], int d)
{
	ff_poly_zn_mod_setup(_ff_p, 0);
	zn_array_mul (h, f, d+1, f, d+1, _ff_poly_zn_mod_ctx);
}

void ff_poly_print_raw(ff_t f[], int d)
{
	register int i;
	
	if ( d == 0 ) { printf ("%ld;\n", f[0]); return; }
	if ( d == 1 ) { printf ("%ld*x + %ld;\n", f[1], f[0]); return; }
	printf ("%ld*x^%d", f[d], d);
	for ( i = d-1 ; i > 1 ; i-- ) if ( f[i] ) printf (" +%ld*x^%d", f[i], i);
	printf (" + %ld*x + %ld;\n", f[1], f[0]);
}

void ffpi_alloc_ctx (ffpi_ctx_t ctx, int n)
{
	int *cntp;
	register int i, j;

	// allocate one extra entry in each work variable, just in case we want to store a poly of degree n
	ctx->c = malloc((n+1)*sizeof(ff_t));
	ctx->r = malloc((n+1)*sizeof(ff_t));
	ctx->s = malloc((n+1)*sizeof(ff_t));
	ctx->t = malloc((2*n+4+(n+1)/FF_PRODUCT_TREE_BASE_CUTOFF)*sizeof(ff_t));
	ctx->mbase = malloc (n*FFPI_MAX_BASE_N*sizeof(ff_t));
	ctx->tree_cnts = malloc(n*sizeof(int));
	ctx->n = n;
	cntp = ctx->tree_cnts;
	for ( i = 0, j = n ; j > FFPI_MAX_BASE_N ; i++, j = j-j/2 );
	if ( i > FFPI_MAX_LEVELS ) { err_printf ("Exceeded FFPI_MAX_LEVELS=%d with n=%d in ffpi_alloc_ctx\n", FFPI_MAX_LEVELS, n); abort(); }
	ctx->levels = i;  ctx->level_cnts[i] = cntp++;  ctx->level_lens[i] = 1;  ctx->level_cnts[i][0] = n;
	for ( i-- ; i >= 0 ; i-- ) {
		ctx->level_cnts[i] = cntp;  ctx->level_lens[i] = 2*ctx->level_lens[i+1]; cntp += ctx->level_lens[i];
		for ( j = 0 ; j < ctx->level_lens[i+1] ; j++ ) {
			ctx->level_cnts[i][2*j] = ctx->level_cnts[i+1][j]/2;
			ctx->level_cnts[i][2*j+1] = ctx->level_cnts[i+1][j] - ctx->level_cnts[i+1][j]/2;
		}
	}
	ctx->mtree[0] = malloc ((n+ctx->level_lens[0])*sizeof(ff_t));
	for ( i = 1 ; i < ctx->levels ; i++ ) ctx->mtree[i] = malloc ((n+ctx->level_lens[i])*sizeof(ff_t));
}

void ffpi_free_ctx (ffpi_ctx_t ctx)
{
	free (ctx->mbase);
	free (ctx->s);
	free (ctx->t);
}

// computes s[i] = 1 / prod_{j!=ii} (x[i]-x[j]}, where the x[i] must be distinct (not verified)  s and x cannot overlap!
// currently this is implemented naively - we could do this in O(M(n)log(n)) time
void ffpi_compute_si_naive (ff_t s[], ff_t x[], int n)
{
	register ff_t t0, t1;
	register int i, j;
	
	for ( i = 0 ; i < n ; i++ ) {
		_ff_set_one(t1);
		for ( j = 0 ; j < i ; j++ ) { _ff_sub(t0,x[i],x[j]);  ff_mult(t1,t1,t0); }
		for ( j++ ; j < n ; j++ ) { _ff_sub(t0,x[i],x[j]);  ff_mult(t1,t1,t0); }
		_ff_set(s[i],t1);
	}
	ff_parallel_invert(s,s,n);
}

/*
	This is a simple non-recursive version of Steps 1 and 2 of Alg. 10.11 in von zur Gathen and Gerhard.
	For the moment we just use a naive O(n^2) approach, since we expect ffpi_interpolate to be called much more often (e.g. O(n) times)
	We may want to change this in the future if it will be used in places where only one or a few interpolations use the same abcissa
*/
void ffpi_setup_ctx (ffpi_ctx_t ctx, ff_t x[], int n, int combo_only)
{
	register ff_t *b, *m, *xx;
	register int i, j, k, u, v;

	if ( n != ctx->n ) { printf ("Error, inconsistent n=%d != %d in ffpi_ctx\n", n, ctx->n); abort(); }
//printf ("x= "); for ( i = 0 ; i < n ; i++ ) printf ("%ld ", _ff_get_ui(x[i])); puts("");
	if ( ! combo_only ) ffpi_compute_si_naive (ctx->s, x, n);
	xx = x;  m = ctx->mtree[0]; b = ctx->mbase;
	for ( u = 0 ;u < ctx->level_lens[0] ; u++ ) {
		k = ctx->level_cnts[0][u];
		// compute m = prod_i (x-x[i])
		ff_poly_from_roots_big (m,xx,k,ctx->t);
//printf("m="); ff_poly_print(m,k);
		/*
			Compute polys m_i=m/(x-[i]) and transpose coefficients so that the first n entries in mbase are the constant coefficients, the next n entries the x coefficients, etc...
			Store temporarily in ctx->c.  This will make forming linear combinations of the m_i via dot-products quick and cache-friendly
		*/
		for ( i = 0 ; i < k ; i++ ) {
			ff_poly_remove_root(ctx->c,m,k,xx++);
//printf("m/(x-x[%d]) = ", i); ff_poly_print(ctx->c,k-1);
			for ( j = 0 ; j < k ; j++ ) _ff_set(b[j*k+i],ctx->c[j]);
		}
		m+=(k+1); b += k*k;
	}
	if ( ! ctx->levels ) return;
	
	// convert level 0 m polys out of Montgomery rep for calls to ff_poly_mult_bigx
	m = ctx->mtree[0];
	for ( u = 0 ; u < ctx->level_lens[0] ; u++ ) {
		k = ctx->level_cnts[0][u];
		for ( i = 0 ; i <= k ; i++ ) m[i] = _ff_get_ui(m[i]);
//printf("m="); ff_poly_print_raw(m,k);
		m += (k+1);
	}
	for ( v = 1 ; v < ctx->levels ; v++ ) {
		m = ctx->mtree[v];  b = ctx->mtree[v-1];
		for ( u = 0 ; u < ctx->level_lens[v] ; u++ ) {
			j = ctx->level_cnts[v-1][2*u];  k = ctx->level_cnts[v-1][2*u+1];
			if ( ctx->level_cnts[v][u] != j+k ) { printf("Inconsistent level cnts in ffpi_ctx %d != %d!\n", ctx->level_cnts[v][u], j+k); abort(); }
			ff_poly_mult_bigx (m, b, j, b+j+1, k);
//printf("m="); ff_poly_print_raw(m,k);
			b += j+k+2;
			m += j+k+1;
		}
	}
}

static inline int imin(int a, int b) { return (a < b ? a : b); }

// o specifies the number of output coefficients, which may be less than the number of input coordinates (yielding a slight speedup)
void _ffpi_interpolate (ff_t f[], ff_t y[], int n, int o, int combo, ffpi_ctx_t ctx)
{
	ff_t *c, *r, *b, *m;
	register int i,j,k,u,v;
	
	if ( n != ctx->n ) { printf ("Error, inconsistent n in ffpi_ctx\n"); abort(); }
	if ( combo ) {
		memcpy (ctx->c, y, n*sizeof(y[0]));
	} else {
		for ( i = 0 ; i < n ; i++ ) _ff_mult(ctx->c[i],ctx->s[i],y[i]);
	}

	c = ctx->c;  r = ctx->r;  b = ctx->mbase;
//printf ("c = "); for ( i = 0 ; i < n ; i++ ) printf ("%ld ", _ff_get_ui(c[i])); puts ("");
	for ( u = 0 ; u < ctx->level_lens[0] ; u++ ) {
		k = ctx->level_cnts[0][u];
//printf("k[%d]=%d\n",u, k);
//printf("transposed m[%d]\n",u); for ( i = 0 ; i < k ; i++ ) { for ( j = 0 ; j < k ; j++ ) printf ("%ld ",_ff_get_ui(b[k*i+j])); puts(""); }
//printf ("c[%d] = ", u); for ( i = 0 ; i < k ; i++ ) printf ("%ld ", _ff_get_ui(c[i])); puts ("");
		for ( i = 0 ; i < k ; i++, b+=k ) ff_dot_product(r+i,b,c,k);
//printf("r[%d]=",u); ff_poly_print(r,k-1);
		c += k; r+= k;
	}
	if ( ctx->levels == 0 ) {
		for ( i = 0 ; i < o ; i++ ) _ff_set(f[i],ctx->r[i]);
		return;
	} else {
		// switch out of Montgomery rep for calls to ff_poly_mult_bigx
		for ( i = 0 ; i < n ; i++ ) ctx->r[i] = _ff_get_ui(ctx->r[i]);
	}
	for ( v = 1 ; v <= ctx->levels ; v++ ) {
		m = ctx->mtree[v-1];  r = ctx->r;
		for ( u = 0 ; u < ctx->level_lens[v] ; u++ ) {
			j = ctx->level_cnts[v-1][2*u];  k = ctx->level_cnts[v-1][2*u+1];
//printf ("m1:= "); ff_poly_print_raw(m+j+1,k);
//printf ("r0:= "); ff_poly_print_raw(r,j-1);
			ff_poly_mult_bigx (ctx->t, m+j+1, imin(k,o), r, imin(j-1,o));
//printf ("t0:= "); ff_poly_print_raw (ctx->t, j+k-1);
//printf ("m0:= "); ff_poly_print_raw(m,j);
//printf ("r1:= "); ff_poly_print_raw(r+j,k-1);
			ff_poly_mult_bigx (ctx->t+j+k,m, imin(j,o), r+j, imin(k-1,o));
//printf ("t1:= "); ff_poly_print_raw (ctx->t+j+k, j+k-1);
			for ( i = 0 ; i < j+k ; i++ ) _ff_add(r[i],ctx->t[i],ctx->t[j+k+i]);
//printf("r:= "); ff_poly_print_raw(r,j+k-1);
			m += j+k+2;
			r += j+k;
		}
	}
	for ( i = 0 ; i < o ; i++ ) _ff_set_ui(f[i],ctx->r[i]);
}

/*
	This recursive version of Alg 10.9 is quite slow, use ffpi_combo instead.
*/

// t should have space for 4n+2lg(n) entries
void ff_poly_linear_combo (ff_t f[], ff_t x[], ff_t c[], ff_t *m, ff_t t[], int n)
{
	ff_t *r0, *r1, *m0, *m1;
	int n0, n1;
	register int i;

	if ( n == 1 ) {
		_ff_set(f[0],c[0]);
		if ( m ) { _ff_neg(m[0],x[0]);  m[1]=1; /* _ff_set_one(m[1]);*/ }		// note we use non-montgomery rep here
//if ( m ) printf ("m[%d] = ", n); ff_poly_print_raw(m,n);
		return;
	}
	n0 = n/2;  n1 = n-n0;
	r0 = t; r1 = t+n0;  m0 = t+n; m1 = t+n+n0+1;
	ff_poly_linear_combo (r0, x, c, m0, t+2*n+2, n0);
	ff_poly_linear_combo (r1, x+n0, c+n0, m1, t+2*n+2, n1);
	if ( m ) ff_poly_mult_bigx(m,m0,n0,m1,n1);
//if ( m ) { printf ("m[%d] = ", n); ff_poly_print_raw(m,n); }
//printf ("r0 = "); ff_poly_print_raw(r0,n0-1);
//printf ("m1 = "); ff_poly_print_raw(m1,n1);
	ff_poly_mult_bigx(f,r0,n0-1,m1,n1);
//printf ("result = "); ff_poly_print_raw(f,n-1);
//printf ("r1 = "); ff_poly_print_raw(r1,n1-1);
//printf ("m0 = "); ff_poly_print_raw(m0,n0);
	ff_poly_mult_bigx(m1,r1,n1-1,m0,n0);
//printf ("result = "); ff_poly_print_raw(m1,n-1);
	for ( i = 0 ; i < n ; i++ ) _ff_addto(f[i],m1[i]);
//printf("f[%d] = ", n); ff_poly_print_raw(f,n-1);
}

/*
	This is slow, use ffpi_interpolate
*/

//  If flags is 1, the values s[i]=1/prod_{j \ne i}x[j] are assume to have been precomputed, t should have space for 5n+2lg(n) entries
void ff_poly_recursive_interpolate (ff_t f[], ff_t x[], ff_t y[], ff_t s[], ff_t t[], int n, int flags)
{
	ff_t *xx;
	register ff_t t0, t1;
	register int i,j;

//printf ("Interpolating from points: "); for ( i = 0 ; i < n ; i++ ) printf ("(%ld,%ld)  ", _ff_get_ui(x[i]), _ff_get_ui(y[i])); puts("");
	if ( ! flags ) {
		for ( i = 0 ; i < n ; i++ ) {
			_ff_set_one(t1);
			for ( j = 0 ; j < n ; j++ ) {
				if ( j==i ) continue;
				_ff_sub(t0,x[i],x[j]);
				ff_mult(t1,t1,t0);
			}
			_ff_set(t[i],t1);
		}
		ff_parallel_invert(s,t,n);
//printf ("s[] = "); for ( i = 0 ; i < n ; i++ ) printf ("%ld ", _ff_get_ui(s[i])); puts("");
	}	
	
	for ( i = 0 ; i < n ; i++ ) { _ff_mult(t0,s[i],y[i]); t[i]=_ff_get_ui(t0); }
	xx = t+n;
	for ( i = 0 ; i < n ; i++ ) xx[i] = _ff_get_ui(x[i]);
//printf ("c[] = "); for ( i = 0 ; i < n ; i++ ) printf ("%ld ", t[i]); puts("");
	ff_poly_linear_combo (f,xx,t,0,t+2*n, n);
	for ( i = 0 ; i < n ; i++ ) _ff_set_ui(f[i],f[i]);
	
//printf ("f = "); ff_poly_print(f,n-1);
	for ( i = 0 ; i < n ; i++ ) {
		ff_poly_eval(t,f,n-1,x+i);
		if ( ! _ff_equal(t[0],y[i]) ) { printf ("Interpolation failed, f(%ld)=%ld != %ld\n", _ff_get_ui(x[i]), _ff_get_ui(t[0]), _ff_get_ui(y[i])); abort(); }
	}

}


/*
	Historical poly_from_roots code starts here.
*/


int FFT_POLYROOTS_MKS_BASE = 6;		// hardwired below, you must change ff_poly_from_roots call also
int FFT_POLYROOTS_MKS_CUTOFF =4;		
int FFT_POLYROOTS_KC_CUTOFF = 110;
int FFT_POLYROOTS_FK_CUTOFF1 = 257;
int FFT_POLYROOTS_FK_CUTOFF2 = 257;//2049;
int FFT_POLYROOTS_KURATSUBA_BASE = 4;
int FFT_POLYROOTS_FFT_BASE = 8;

static inline void zn_array_mul_x (unsigned long *res, unsigned long *op1, size_t n1, unsigned long *op2, size_t n2, zn_mod_t mod)
{
	register size_t z1, z2, m1, m2;
	
	for ( z2 = 0 ; z2 < n2 && ! op2[z2] ; z2++ );  m2 = n2-z2;
	if ( ! m2 ) { memset(res,0,sizeof(unsigned long)*(n1+n2-1)); return; }
	for ( z1 = 0 ; z1 < n1 && ! op1[z1] ; z1++ );  m1 = n1-z1;  m2 = n2-z2;
	if ( ! m1 ) { memset(res,0,sizeof(unsigned long)*(n1+n2-1)); return; }
	if ( m1 < m2 ) {
		zn_array_mul (res+z1+z2, op2+z2, m2, op1+z1, m1, mod);
	} else {
		zn_array_mul (res+z1+z2, op1+z1, m1, op2+z2, m2, mod);
	}
	for ( m1 = 0 ; m1 < z1+z2 ; m1++ ) res[m1] = 0;
}

/*
	We pad to a power of 2 but rely on zn_array_mul_x to be smart about stripping the padding away.
	In theory, this algorithm would be faster if it did a balanced recursion (the current approach can lead to very
	unbalanced multiplications for certain d), but in practical tests it is not clear that this would give much
	improvement (if any).  For simplicity we stick with our power of 2 approach using zn_array_mul_x
	(note that stripping the padding *is* a big improvment, 30-40% in some cases).
*/

void ff_poly_mks_from_roots_n (ff_t f[], ff_t r[], int d, ff_t t[])
{
	register ff_t *s, *u;
	register int i, j, k, m1, m2, N;

	if ( d <= FFT_POLYROOTS_MKS_CUTOFF ) { ff_poly_from_roots_naive (f, r, d); return; }

	k = ui_len(d-1);
	N = (1<<k);
	s = t+N;  u = s+N;
	m1 = 1<<FFT_POLYROOTS_MKS_BASE;  m2 = m1<<1;
	if ( N <= m1 ) {
		for ( i = d ; i < N ; i++ ) _ff_set_zero(r[i]);
		switch (k) {
		case 3: ff_poly_from_roots_8(f,r); break;
		case 4: ff_poly_from_roots_16(f,r); break;
		case 5: ff_poly_from_roots_32(f,r); break;
		case 6: ff_poly_from_roots_64(f,r); break;
		}
		for ( i = 0 ; i < d ; i++ ) _ff_set(f[i], f[i+N-d]);
		_ff_set_one(f[d]);
		return;
	}
	
	k = (d%m1? m1*((d/m1)+1) : d);
	for ( i = d ; i < k ; i++ ) _ff_set_zero (r[i]);							// zero pad list of roots
	for ( i = 0 ; i < k ; i+= m1 ) ff_poly_from_roots_64 (f+i,r+i);				// setup monic base polys (monic coeff is implicit and overwritten)
	for ( ; i < N ; i++ ) _ff_set_zero(f[i]);								// zero pad list of polys
	for ( i = 0 ; i < k ; i++ ) f[i]=_ff_get_ui(f[i]);							// convert out of Montgomery (note addition still works in std rep)
	ff_poly_zn_mod_setup(_ff_p, d);
	for ( k=FFT_POLYROOTS_MKS_BASE+1 ; m1<N ; m1 <<= 1, m2 <<= 1, k++ ) {
		for ( i = 0 ; i < d ; i+= m2 ) {
			for ( j = 0 ; j < m1 ; j++ ) _ff_add(u[j],f[i+j],f[i+m1+j]);			// add pair of polys (minus monic coefficient) into u before multiplying	
			zn_array_mul_x (s, f+i, m1, f+i+m1, m1, _ff_poly_zn_mod_ctx);	// call Harvey's multipoint Kronecker substitution code for the heavy lifting
			_ff_set_zero(s[m2-1]);									// zn_array_mul returns deg 2m1-2 poly, we need to clear the top coefficient
			for ( j = 0 ; j < m1 ; j++ ) {								// copy/add product + shifted sum back into f:  (x^m1+ f)(x^m1+g) = x^m2 + (f+g)x^m1 + fg
				_ff_set(f[i+j],s[j]);
				_ff_add(f[i+m1+j],s[m1+j],u[j]);
			}
		}
	}
	k = N-d;
	for ( i = 0 ; i < d ; i++ )  _ff_set_ui (f[i],f[i+k]);						// shift and convert back into Montgomery rep
	_ff_set_one (f[d]);
}


void ff_poly_mks_from_roots_o (ff_t f[], ff_t r[], int d, ff_t t[])
{
	register ff_t *s, *u;
	register int i, j, k, m1, m2, N;

	if ( d <= FFT_POLYROOTS_MKS_CUTOFF ) { ff_poly_from_roots_naive (f, r, d); return; }

	k = ui_len(d-1);
	N = (1<<k);
	s = t+N;  u = s+N;
	m1 = 1<<FFT_POLYROOTS_MKS_BASE;  m2 = m1<<1;
	if ( N <= m1 ) {
		for ( i = d ; i < N ; i++ ) _ff_set_zero(r[i]);
		switch (k) {
		case 3: ff_poly_from_roots_8(f,r); break;
		case 4: ff_poly_from_roots_16(f,r); break;
		case 5: ff_poly_from_roots_32(f,r); break;
		case 6: ff_poly_from_roots_64(f,r); break;
		}
		for ( i = 0 ; i < d ; i++ ) _ff_set(f[i], f[i+N-d]);
		_ff_set_one(f[d]);
		return;
	}
	
	k = (d%m1? m1*((d/m1)+1) : d);
	for ( i = d ; i < k ; i++ ) _ff_set_zero (r[i]);							// zero pad list of roots
	for ( i = 0 ; i < k ; i+= m1 ) ff_poly_from_roots_64 (f+i,r+i);				// setup monic base polys (monic coeff is implicit and overwritten)
	for ( ; i < N ; i++ ) _ff_set_zero(f[i]);								// zero pad list of polys
	for ( i = 0 ; i < k ; i++ ) f[i]=_ff_get_ui(f[i]);							// convert out of Montgomery (note addition still works in std rep)
	ff_poly_zn_mod_setup(_ff_p, d);
	for ( k=FFT_POLYROOTS_MKS_BASE+1 ; m1<N ; m1 <<= 1, m2 <<= 1, k++ ) {
		for ( i = 0 ; i < d ; i+= m2 ) {
			for ( j = 0 ; j < m1 ; j++ ) _ff_add(u[j],f[i+j],f[i+m1+j]);			// add pair of polys (minus monic coefficient) into u before multiplying
			zn_array_mul (s, f+i, m1, f+i+m1, m1, _ff_poly_zn_mod_ctx);		// call Harvey's multipoint Kronecker substitution code for the heavy lifting
			_ff_set_zero(s[m2-1]);									// zn_array_mul returns deg 2m1-2 poly, we need to clear the top coefficient
			for ( j = 0 ; j < m1 ; j++ ) {								// copy/add product + shifted sum back into f:  (x^m1+ f)(x^m1+g) = x^m2 + (f+g)x^m1 + fg
				_ff_set(f[i+j],s[j]);
				_ff_add(f[i+m1+j],s[m1+j],u[j]);
			}
		}
	}
	k = N-d;
	for ( i = 0 ; i < d ; i++ )  _ff_set_ui (f[i],f[i+k]);						// shift and convert back into Montgomery rep
	_ff_set_one (f[d]);
}


/*
	The rest of the code in this module is for testing and benchmarking, it is not currently used and not fully debugged.
*/

int fft_precompute (fft_ctx_t ctx, int n, unsigned long maxint)
{
	ff_t *s, t;
	unsigned long bytes;
	register int i, j, N;
	
	if ( n > FFT_MAX_N || maxint >= FFT_DOUBLE_MAXINT ) return 0;		// ensures epsilon < 1/4 in modular CRT
	N = (1<<n);
	ctx->n = n;
	ctx->primes = ( maxint > FFT_SINGLE_MAXINT ? 2 : 1 );
	bytes = (ctx->primes+3)*N*sizeof(ff_t);
	ctx->t = malloc(bytes);
	if ( ! ctx->t ) { printf ("memory allocation of %ld bytes failed in ff_poly_fft_precomputes\n", bytes); abort(); }
	ctx->p[0] = ( ctx->primes > 1 ? FFT_PRIME_1 : FFT_PRIME_0 );
	ctx->p[1] = ( ctx->primes > 1 ? FFT_PRIME_2 : 0 );
	s = ctx->t + 3*N;
	for ( i = 0 ; i < ctx->primes ; i++ ) {
		ctx->w[i] = s;  s+= N/2; ctx->wi[i] = s;  s+= N/2;
		ff_setup_ui (ctx->p[i]);
		if ( ctx->primes == 2 ) {
			_ff_set_ui(t,ctx->p[1-i]); ff_invert (t,t);
			ctx->a[i] = _ff_get_ui(t); ctx->y[i] = (double)ctx->a[i]/(double)ctx->p[i]; 
			if ( i == 0 ) { _ff_set_ui(t,FFT_ROOT_1); } else { _ff_set_ui(t,FFT_ROOT_2); }
		} else {
			_ff_set_ui(t,FFT_ROOT_0);
		}
		_ff_set_one(ctx->w[i][0]);  _ff_set_one(ctx->wi[i][0]);
		for ( j = 0 ; j < (ctx->primes>1?FFT_H_2:FFT_H_1) - n ; j++ ) ff_square(t,t);
		_ff_set(ctx->w[i][1],t);
		for ( j = 2 ; j < N/2 ; j++ ) _ff_mult(ctx->w[i][j], ctx->w[i][j-1],ctx->w[i][1]);
		for ( j = 1 ; j < N/2 ; j++ ) _ff_neg(ctx->wi[i][j], ctx->w[i][N/2-j]);
		_ff_set_one(ctx->half[i][0]);  _ff_set(ctx->half[i][1],_ff_half);
		for ( j = 1 ; j <= n ; j++ ) _ff_mult(ctx->half[i][j],ctx->half[i][j-1],_ff_half);
	}
	return 1;
}

void fft_poly_mult_ff (ff_t h[], ff_t f[], ff_t g[], int n, fft_ctx_t ctx)
{
	unsigned long p;
	register double r;
	register int i, m, N;

	N = (1<<n);
	m = N>>1;
	if ( _ff_p < FFT_SINGLE_MAXINT ) {
		p = _ff_p;
		for ( i = 0 ; i < m ; i++ ) { f[i] = _ff_get_ui(f[i]); g[i] = _ff_get_ui(g[i]); }					// pull current field elements into Z, we don't bother with zeroes
//printf ("f/Z = "); for ( i = (1<<n)-1 ; i >= 0 ; i-- ) printf ("%ld*X^%d  ", f[i], i); puts("");
//printf ("g/Z = "); for ( i = (1<<n)-1 ; i >= 0 ; i-- ) printf ("%ld*X^%d  ", g[i], i); puts("");
		fft_poly_mult_z (h, f, g, n, 0, ctx);												// multiply to get results in Z (reduced mod fft_p[0])
//printf ("h/Z = "); for ( i = (1<<n)-1 ; i >= 0 ; i-- ) printf ("%ld*X^%d  ", h[i], i); puts("");
		ff_setup_ui(p);
		for ( i = 0 ; i < N ; i++ ) _ff_set_ui(h[i],h[i]);										// reduce back into orginal finite field
	} else {
		register ff_t *s, *u, *v;
		register ff_t t0, t1, t2, t3, t4, t5;
		
		/*
			We use the modified CRT to compute the product mod ff_p by computing it mod fft_p[0] and fft_p[1] as described in section 5.1 of
			Agashe, Lauter, and Venkatesan "Constructing elliptic curves with a known number of points over a prime field" (2003).
			(this algorithm is originally due to Montgomery/Silverman "And FFT extension to the P-1 factoring algorithm", we use K=2 here)
		*/
		p = _ff_p;
		s = ctx->t;  u = s+N;  v = u+N;
		for ( i = 0 ; i < m ; i++ ) { u[i] = f[i] = _ff_get_ui(f[i]); v[i] = g[i] = _ff_get_ui(g[i]); }		// pull current field elements into Z and save a copy , but don't bother with zeroes
		for ( ; i < N ; i++ ) { _ff_set_zero(ctx->t[i]); _ff_set_zero(ctx->t[N+i]); }				// zero pad copy
		fft_poly_mult_z (h, f, g, n, 0, ctx);												// multiply to get results in Z (reduced mod fft_p[0])
//printf ("h/%ld = ", ctx->p[0]); for ( i = (1<<n)-1 ; i >= 0 ; i-- ) printf ("%ld*x^%d  ", h[i], i); puts("");
		fft_poly_mult_z (s, u, v, n, 1, ctx);												// multiply again to get results in Z (reduced mod fft_p[1])
//printf ("h/%ld = ", ctx->p[1]); for ( i = (1<<n)-1 ; i >= 0 ; i-- ) printf ("%ld*x^%d  ", s[i], i); puts("");
		ff_setup_ui(p);
		_ff_set_ui(t0,ctx->p[1]); _ff_set_ui(t1,ctx->p[0]); _ff_mult(t2,t0,t1);					// t2=M=p[0]*p[1] mod ff_p
		_ff_set_ui(t3,ctx->a[0]);  ff_mult(t0,t0,t3);										// t1=a_0*M_0 (note M_0=p[1])
		_ff_set_ui(t3,ctx->a[1]);  ff_mult(t1,t1,t3);										// t1=a_1*M_1 (note M_1=p[])			
		for ( i = 0 ; i < N ; i++ ) {
			_ff_set_ui(t3,h[i]);  _ff_mult(t4,t3,t0);
			_ff_set_ui(t3,s[i]);  _ff_mult(t5,t3,t1);  _ff_addto(t4,t5);							// t4 = a0x0M0+a1x1M1 mod ff_p
			r = ctx->y[0]*h[i]+ctx->y[1]*s[i];											// compute r = sum a_i*x_i/m_i 
			_ff_set_ui(t3,(unsigned long)(r+0.5));										// round r to nearest integer and pull into ff_p
			_ff_mult(t5,t2,t3);														// t5 = rM mod ff_p
//printf ("x0=%ld a0=%ld m0=%ld M0=%ld   ", h[i], ctx->a[0], ctx->p[0], ctx->p[1]);
//printf ("x1=%ld a1=%ld m1=%ld M1=%ld   ", s[i], ctx->a[1], ctx->p[1], ctx->p[0]);
//printf ("r=%.6f\n", r);
			_ff_sub (h[i],t4,t5);														// place final result in h (in ff_p)
		}
	}
}


// NOTE: the inputs and outputs are treated as integers stored in type ff_t, not field elements, and must be less than the value of maxint specified in ff_poly_fft_precompute
// output values are positive integer residues mod fft_p[i]
void fft_poly_mult_z (ff_t h[], ff_t f[], ff_t g[], int n, int i, fft_ctx_t ctx)
{
	register int j;
	
	if ( _ff_p != ctx->p[i] ) ff_setup_ui(ctx->p[i]);
	for ( j = 0 ; j < (1<<(n-1)) ; j++ ) { _ff_rset_ui(f[j],f[j]);  _ff_rset_ui(g[j],g[j]); }		// (trivially) reduce into finite field
	_fft_poly_conv_mult (h, f, g, ctx->w[i], ctx->wi[i], ctx->half[i], n, ctx->n-n);
	for ( j = 0 ; j < (1<<n) ; j++ ) h[j] = _ff_get_ui(h[j]);							// lift mod fft_p result back to Z
}

/*
	Computes the product h=fg, where f and g have degree < 2^{n-1}, destroying f and g in the process.
	w[i]=w^i where w is a 2^n-th root of unity, for i from 0 to 2^{n-1}-1.  wi are the inverses of w
	n must be at least 3

	There is a subtle bug here somewhere, haven't bothered to find it since switching to zn_poly
*/
void _fft_poly_conv_mult (ff_t h[], ff_t f[], ff_t g[], ff_t w[], ff_t wi[], ff_t half[], int n, int k)
{
	register ff_t t1, t2, t3, t4;
	register int i, j, k1, k2, k3, m, m2, s, N;

//printf ("fft_poly_conv_mult n=%d, k=%d\n", n, k);  printf ("    "); ff_poly_print(f,(1<<(n-1))-1);  printf("    "); ff_poly_print(g,(1<<(n-1))-1);	
	N = 1<<n; 
	// descend
	m2 = N>>1;  m = m2>>1;  s = k+1;
	for ( i = 0 ; i < m2 ; i++ ) { _ff_mult(f[m2+i],f[i],w[i<<k]);  _ff_mult(g[m2+i],g[i],w[i<<k]); }		// unwind the first iteration of the loop, since top half of f and g are zero (for mult)
	while ( m > 2 ) {
		for ( i = 0 ; i < N ; i += m2 ) {
			for ( j = 0 ; j < m ; j++ ) {
				k1 = i+j;  k2 = k1+m; k3 = j<<s;
				_ff_set(t2,f[k2]); _ff_sub(t1,f[k1],t2); _ff_addto(f[k1],t2);  _ff_mult(f[k2],t1,w[k3]);
				_ff_set(t2,g[k2]); _ff_sub(t1,g[k1],t2); _ff_addto(g[k1],t2);  _ff_mult(g[k2],t1,w[k3]);
			}
		}
		m >>= 1;  m2 >>= 1;  s++;
	}
	// base - this could be optimized further, but keep it simple for now.
	for ( i = 0 ; i < N ; i+=4 ) {
		_ff_mult(t1,f[i],g[i]);  _ff_mult(t2,f[i+1],g[i+3]);  _ff_mult(t3,f[i+2],g[i+2]);  _ff_mult(t4,f[i+3],g[i+1]);  _ff_addto(t1,t2); _ff_addto(t3,t4);  _ff_add(h[i],t1,t3);
		_ff_mult(t1,f[i],g[i+1]);  _ff_mult(t2,f[i+1],g[i]);  _ff_mult(t3,f[i+2],g[i+3]);  _ff_mult(t4,f[i+3],g[i+2]);  _ff_addto(t1,t2); _ff_addto(t3,t4);  _ff_add(h[i+1],t1,t3);
		_ff_mult(t1,f[i],g[i+2]);  _ff_mult(t2,f[i+1],g[i+1]);  _ff_mult(t3,f[i+2],g[i]);  _ff_mult(t4,f[i+3],g[i+3]);  _ff_addto(t1,t2); _ff_addto(t3,t4);  _ff_add(h[i+2],t1,t3);
		_ff_mult(t1,f[i],g[i+3]);  _ff_mult(t2,f[i+1],g[i+2]);  _ff_mult(t3,f[i+2],g[i+1]);  _ff_mult(t4,f[i+3],g[i]);  _ff_addto(t1,t2); _ff_addto(t3,t4);  _ff_add(h[i+3],t1,t3);
	}
	// ascend
	m = 4; m2 = 8; s = k+n-3;
	while ( m < N ) {
		for ( i = 0 ; i < N ; i += m2 ) {
			for ( j = 0 ; j < m ; j++ ) {
				k1 = i+j;  k2 = k1+m; k3 = j<<s;
				_ff_mult(t2,h[k2],wi[k3]); _ff_sub (h[k2],h[k1],t2);  _ff_addto(h[k1],t2);
			}
		}
		m <<= 1;  m2 <<= 1;  s--;
	}
	for ( i = 0 ; i < (1<<n) ; i++ ) ff_mult(h[i],h[i],half[n-2]);								// adjust half power by 2 to reflect the fact that we didn't double in two base case layers
//printf("    "); ff_poly_print(h,(1<<n)-1);	
}

/*
	Multiplies f and g of degree < 2^{n-1} into h of degree <2^n-1.  f is destroyed in the process, g is not modified.
	The arrays f, g are size 2^{n-1}, while h and t are size 2^n (t is workspace, contents ignored), and can't overlap.
	High coefficients of f and g must be zero-filled below 2^{n-1} as needed.

	n must be at least 2.
*/
void ff_poly_Karatsuba_mult (ff_t h[], ff_t f[], ff_t g[], ff_t t[], int n)
{
	register ff_t t0, t1;
	register int i, m2, m1, m3;
	
	// use classical for base case, 4M+1A is slighltly better than 3M+4A
	if ( n==2 ) { _ff_set_zero(h[3]); _ff_mult(h[2],f[1],g[1]); _ff_mult(h[0],f[0],g[0]); _ff_mult(t0,f[0],g[1]); _ff_mult(t1,f[1],g[0]); _ff_add(h[1],t0,t1);  return; }
	
	m2 = (1<<(n-1));  m1 = m2>>1;

	// compute (f0+f1) and (g0+g1) and save them in the first half of t (we need to do this because f and g will get trashed in the recursion)
	for ( i = 0 ; i < m1 ; i++ ) { _ff_add(t[i],f[i],f[m1+i]);  _ff_add(t[m1+i],g[i],g[m1+i]); }
	ff_poly_Karatsuba_mult (h, f, g, t+m2, n-1);						// compute h0 = f0g0 (trashing f0 and g0 in the process)
	ff_poly_Karatsuba_mult (h+m2, f+m1, g+m1, t+m2, n-1);			// compute h1= f1g1 (trashing f1 and g1 in the process)
	ff_poly_Karatsuba_mult (f, t, t+m1, t+m2, n-1);					// compute h2= (f0+f1)(g0+g1) and store in f
	m3 = m2+m1;
	for ( i = 0 ; i < m1 ; i++ ) {									// add (h2-h0-h1)x^m2 into h
		_ff_sub(t0,f[i],h[i]);  _ff_subfrom(t0,h[m2+i]);					// t0 = (h2-h0-h1)[i]
		_ff_sub(t1,f[m1+i],h[m1+i]);  _ff_subfrom(t1,h[m3+i]);			// t1 = (h2-h0-h1)[m1+i]
		_ff_addto(h[m1+i],t0);  _ff_addto(h[m2+i],t1);				// update h
	}
}

/*
	Given d roots in f, replaces f with poly of degree d+1 having these roots.  f should have space for min(N,d+1) coefficients, where N is the least power of 2 >= d.
	t is workspace for 5N/2 coefficients in F_p
*/
void ff_poly_Karatsuba_from_roots (ff_t f[], ff_t r[], ff_t t[], int d)
{
	register ff_t *s, *u;
	register int i, j, k, m1, m2, N;

	if ( d < FFT_POLYROOTS_KC_CUTOFF || (d>128&&d<190) ) { ff_poly_from_roots_naive (f, r, d); return; }

	N = (1<<ui_len(d-1));  s = t+N;  u = s+N;
	m1 = 1<<FFT_POLYROOTS_KURATSUBA_BASE;  m2 = m1<<1;
	k = (d%m1? m1*((d/m1)+1) : d);
	for ( i = d ; i < k ; i++ ) _ff_set_zero (r[i]);							// zero pad list of roots
	// setup monic base polys (the monic coefficient is implicit and overwritten by the low coefficient of the next poly)
	for ( i = 0 ; i < k ; i+= m1 ) ff_poly_from_roots_naive (f+i,r+i,m1);
	for ( ; i < N ; i++ ) _ff_set_zero(f[i]);								// zero pad list of polys
	for ( k=FFT_POLYROOTS_KURATSUBA_BASE+1 ; m1<N ; m1 <<= 1, m2 <<= 1, k++ ) {
		for ( i = 0 ; i < d ; i+= m2 ) {
			for ( j = 0 ; j < m1 ; j++ ) _ff_add(u[j],f[i+j],f[i+m1+j]);			// add pair of polys (minus monic coefficient) into u before multiplying
			ff_poly_Karatsuba_mult (s, f+i, f+i+m1, t, k);					// multiply pair of polys (minus monic coefficient) into s
			for ( j = 0 ; j < m1 ; j++ ) {								// copy/add product + shifted sum back into f:  (x^m1+ f)(x^m1+g) = x^m2 + (f+g)x^m1 + fg
				_ff_set(f[i+j],s[j]);
				_ff_add(f[i+m1+j],s[m1+j],u[j]);
			}
		}
	}
	if ( d < N ) {													// shift back poly to get degree d result (eliminates padded zero roots)
		k = N-d;
//		for ( i = 0 ; i < k ; i++ ) if  ( ! _ff_zero (f[i]) ) { printf ("Nonzero coefficient %d when zero expected\n", i); abort(); } 
		for ( i = 0 ; i <= d ; i++ )  _ff_set (f[i],f[i+k]);
	}
	_ff_set_one (f[d]);
}


/*
	Given d roots in f, replaces f with poly of degree d+1 having these roots.  f should have space for min(N,d+1) coefficients, where N is the least power of 2 >= d.
	t is workspace for 3N coefficients in F_p
*/
void ff_poly_conv_from_roots (ff_t f[], ff_t r[], ff_t t[], int d, fft_ctx_t ctx)
{
	register ff_t t0, *f1, *f2;
	register int i, j, k, m1, m2, N;
	
	if ( ctx->primes==1 ) {
		if ( d < FFT_POLYROOTS_FK_CUTOFF1 ) { ff_poly_Karatsuba_from_roots(f,r,t,d); return; }
		if ( d < 210 ) { ff_poly_from_roots_naive (f, r, d); return; }
	} else {
		if ( d < FFT_POLYROOTS_FK_CUTOFF2 ) { ff_poly_Karatsuba_from_roots(f,r,t,d); return; }
		if ( ctx->primes==2 && d<310 ) { ff_poly_from_roots_naive (f, r, d); return; }
	}
	N = (1<<ui_len(d-1));  f1 = t+N;  f2 = f1+N;
	m1 = 1<<FFT_POLYROOTS_FFT_BASE;  m2 = m1<<1;
	k = (d%m1? m1*((d/m1)+1) : d);
	for ( i = d ; i < k ; i++ ) _ff_set_zero (r[i]);							// zero pad list of roots to multiple of m1
	// setup monic base polys (the monic coefficient is implicit and overwritten by the low coefficient of the next poly)
	for ( i = 0 ; i < k ; i+= m1 ) ff_poly_from_roots_naive (f+i,r+i,m1);
	for ( ; i < N ; i++ ) _ff_set_zero(f[i]);								// zero pad list of polys
	for ( k=FFT_POLYROOTS_FFT_BASE+1 ; m1<N ; m1 <<= 1, m2 <<= 1, k++ ) {
		for ( i = 0 ; i < d ; i+= m2 ) {
			for ( j = 0 ; j < m1 ; j++ ) { _ff_set(f1[j],f[i+j]);  _ff_set(f2[j],f[i+m1+j]); }
			for ( ; j < m2 ; j++ ) { _ff_set_zero(f1[j]); _ff_set_zero(f2[j]); }	// zero pad for FFT mult
//printf ("f1: X^%d ", m1);  for ( x = m1-1 ; x >= 0 ; x-- ) printf ("+ %ld*X^%d ", _ff_get_ui(f1[x]),x); puts("");
//printf ("f2: X^%d ", m1);  for ( x = m1-1 ; x >= 0 ; x-- ) printf ("+ %ld*X^%d ", _ff_get_ui(f2[x]),x); puts("");
			fft_poly_mult_ff (t, f1, f2, k, ctx);							// multiply pair of polys (minus monic coefficient) into t
//printf ("t: "); ff_poly_print(t,m2);
			for ( j = 0 ; j < m1 ; j++ ) {								// copy/add product + shifted sum back into f:  (x^m1+ f)(x^m1+g) = x^m2 + (f+g)x^m1 + fg
				_ff_add(t0,f[i+j],f[i+m1+j]);
				_ff_set(f[i+j],t[j]);  _ff_add(f[i+m1+j],t[m1+j],t0);
			}
		}
	}
	if ( d < N ) {													// shift back poly to get degree d result (eliminates padded zero roots)
		k = N-d;
//		for ( i = 0 ; i < k ; i++ ) if  ( ! _ff_zero (f[i]) ) { printf ("Nonzero coefficient %d when zero expected\n", i); abort(); } 
		for ( i = 0 ; i <= d ; i++ )  _ff_set (f[i],f[i+k]);
	}
	_ff_set_one (f[d]);
}


/*
	ff_poly_fft computes the DFT of poly of degree < 2^n (coefficients up to 2^n-1 must be set to zero if absent),
	given w[i] = w^i where w is a 2^n-th root of unity in the current field.  Trashes the input in the process.

	If rev=0, then the input poly f is in a and it is assumed to have degree < 2^{n-1}.  The output is put in a with a[i] = f(w^{revbit(i)})
	If rev!=0, then the input poly f is in b, and the output is put in b with b[i] = rev(f)(w^i)
	(where rev(f) is the poly obtained by swapping coefficient i with revbit(i))

	This algorithm is adapted from von zur Gathen and Gerhard, Algorithm 8.16, made non-recursive (see picture on p.232).
	n must be at least 3.  We have unwound the first level and last two levels of the main loop.
*/
void ff_poly_fft (ff_t a[], ff_t w[], int n, int nomult)
{
	register ff_t t,t2;
	register int i, j, k, k2, m, m2, s, N;

	N = (1<<n);
	if ( nomult ) {
		m2 = N;  s = 0;
	} else {
		m2 = N>>1;  s = 1;				// in the multipication case, we assume top half of the input coefficients are zero, so we only need to multiply by w[i]
		for ( i = 0 ; i < m2 ; i++ ) _ff_mult(a[m2+i],a[i],w[i]);
	}
	m = m2>>1;
	while ( m>2 ) {					// we unwind the last two loops below
		for ( i = 0 ; i < N ; i+=m2 )
			for ( j = 0 ; j < m ; j++ ) { k = i+j;  k2 = k+m; _ff_set(t2,a[k2]); _ff_sub(t,a[k],t2); _ff_addto(a[k],t2);  _ff_mult(a[k2],t,w[j<<s]); }
		m >>= 1;  m2 >>= 1;  s++;
	}
	k = (1<<s);
	for ( i = 0 ; i < N ; i+= 4 ) {			// we save 3N mults by handling the last two levels, doesn't make all that much difference though (memory access dominates)
		_ff_sub(t,a[i],a[i+2]); _ff_addto(a[i],a[i+2]);  _ff_set(a[i+2],t);
		_ff_sub(t,a[i+1],a[i+3]); _ff_addto(a[i+1],a[i+3]);  _ff_mult(a[i+3],t,w[k]);
		_ff_sub(t,a[i],a[i+1]); _ff_addto(a[i],a[i+1]);  _ff_set(a[i+1],t);
		_ff_sub(t,a[i+2],a[i+3]); _ff_addto(a[i+2],a[i+3]);  _ff_set(a[i+3],t);
	}
}

/*
	Given a[i]=f(w[2i]) for i from 0 to 2^{n-1}-1, where w[1] is a 2^nth root of unity returns a[i]=f(w[i]) for i from 0 to 2^n.
	Note that f is implied by (and recovered form) a[i] and need not be specified.  b and c are auxillary arrays that need to be of size 2^n
*/
/*
void ff_poly_fft_extend (ff_t a[], ff_t b[], ff_t c[], ff_t w[], ff_t w2[], int n, ff_t n_inv)
{
	register ff_t t;
	register int i, j, k, N;
	
	N = 1<<n;
	for ( i = 0 ; N ; i++ ) _ff_set(a[i+i],b[i]);			// we are given f(w[2i]), we just need to spread them out (and save them from being trashed below
	ff_poly_fft (b, c, w2, n, 0);					// recover coefficients of f
	ff_poly_fft_rev (c, b, n);						// we need to reverse indices (wasn't done by ff_poly_fft)
	for ( i = 0 ; i < N ; i++ ) {					// complete recovery of coefficients of f and multiply by w[i] to form r_1*
		ff_mult(t,c[i],n_inv);	
		ff_mult(c[i],t,w[i]);
	}
	ff_poly_fft (b, c, w2, n, 0);					// compute r_1*(w2[i])
	ff_poly_fft_rev (c, b, n);						// reverse indices
	for ( i = 0 ; i < N ; i++ ) _ff_set(a[i+i+1],c[i]);	// insert in odd positions of output array
}*/

/*
	Multiplies f and g of degree < 2^{n-1} into h of degree <2^n, destroying f and g in the process.
	The arrays f, g, and h are all size 2^n, while w and wi are size 2^{n-1}.  n_inv is 1/n mod _ff_p.  None of the arrays can overlap.
	High coefficients of f and g must be zero-filled as necessary.
*/
void ff_poly_fft_mult (ff_t h[], ff_t f[], ff_t g[], ff_t w[], ff_t wi[], int n, ff_t N_inv)
{
	register int i;
	
	ff_poly_fft (f,w,n,0);
	ff_poly_fft (g,w,n,0);
	for ( i = 0 ; i < (1<<n) ; i++ ) _ff_mult(h[i],f[i],g[i]);
	ff_poly_fft_rev (f,h,n);
	ff_poly_fft (f,wi,n,1);
	ff_poly_fft_rev (h,f,n);
	for ( i = 0 ; i < (1<<n) ; i++ ) ff_mult(h[i],h[i],N_inv);
}

