#ifndef _FFPOLYBIG_INCLUDE_
#define _FFPOLYBIG_INCLUDE_

/*
    Copyright 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include "zn_poly/zn_poly.h"										// this include isn't necessary, but it helps to flag a link-time dependency for ffpolybig.c earlier
#include "ff.h"
#include "ffext.h"
#include "ffpoly.h"
#include "cstd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FF_PRODUCT_TREE_MAX_LEVELS				21
#define FF_PRODUCT_TREE_BASE_CUTOFF			64				// cannot exceed the max supported by ff_poly_from_roots_small, which is currently 64
#define FF_PRODUCT_TREE_MAX_DEGREE			(1<<25)			
#define FF_PRODUCT_TREE_CACHE_SIZE				4

static inline ff_t *ff_poly_from_roots_big_workspace (int d)
	{ return mem_alloc ((2*d+2+(2*d)/FF_PRODUCT_TREE_BASE_CUTOFF)*sizeof(ff_t)); }
void ff_poly_from_roots_big (ff_t f[], ff_t r[], int d, ff_t w[]);				// w needs to hold space for 2d+2+floor(2d/BASE_CUTOFF) entries, w=f  is ok, but w and f cannot overlap r

void ff_poly_square_big (ff_t h[], ff_t f[], int d);
void ff_poly_mult_big (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b);
void ff_poly_mult_bigx (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b);			// assumes integer format, no conversion out of montgomery
void ff_poly_square_bigx (ff_t h[], ff_t f[], int d);						// assumes integer format, no conversion out of montgomery


#define FFPI_MAX_LEVELS	FF_PRODUCT_TREE_MAX_LEVELS
extern int FFPI_MAX_BASE_N;

typedef struct {
	ff_t *c, *r, *s, *t;
	int n;
	int levels;
	int level_lens[FFPI_MAX_LEVELS+1];
	int *level_cnts[FFPI_MAX_LEVELS+1];
	int *tree_cnts;
	ff_t *mtree[FFPI_MAX_LEVELS+1];
	ff_t *mbase;
} ffpi_ctx_t[1];

// strictly speaking, there is no reason the interpolation code needs to be in this module, it could be made independent of zn_poly

void ffpi_compute_si_naive (ff_t s[], ff_t x[], int n);							// computes s[i] = 1 / prod_{j!=i} (x[i]-x[j]), requires x[i] distinct and s and x cannot overlap
void ffpi_alloc_ctx (ffpi_ctx_t ctx, int n);
void ffpi_free_ctx (ffpi_ctx_t ctx);
void ffpi_setup_ctx (ffpi_ctx_t ctx, ff_t x[], int n, int combo_only);				// if combo_only flag is set, only ffpi_combo may be called subsequently (no interpolation)
void _ffpi_interpolate (ff_t f[], ff_t y[], int n, int o, int combo, ffpi_ctx_t ctx);		// if combo is set, computes f(X) = sum y_i g(X)/(X-x_i), otherwise interpolates f(X) s.t. f(x_i)=y_i, computing the first o coefficients
static inline void ffpi_combo (ff_t f[], ff_t c[], int n, ffpi_ctx_t ctx)				// computes f(X) = sum c_i g(X)/(X-x_i)  where g(X) = prod (X-x_i) (the x_i are specified in ffpi_setup_ctx)
	{ _ffpi_interpolate (f,c,n,n,1,ctx); }
static inline void ffpi_interpolate(ff_t f[], ff_t y[], int n, ffpi_ctx_t ctx)			// computes f(X) of degree less than n such that f(x_i)=y_i, overlap is ok
	{ _ffpi_interpolate (f,y,n,n,0,ctx); }
	
/*
	The FFT code is not fully debugged and is not currently used, since zn_poly is faster.
*/


#define FFT_PRIME_0		(197*(1UL<<55)+1)	// we want this as close to 2^63 as we can get, but no bigger (ff.h needs the top bit)
#define FFT_ROOT_0		607
#define FFT_H_1			55
#define FFT_SINGLE_MAXINT	266414583UL			// floor(sqrt(FFT_PRIME_0))/10

/*
	Note that when we use a pair of primes, they can't be too big, to be able to represent p in a double to an accuracy of less than 1/4. 
	Assuming 56 bits of precision, we want p to be comfortably less than 54 bits
*/

#define FFT_PRIME_1		(19*(1UL<<46)+1)
#define FFT_ROOT_1		45
#define FFT_PRIME_2		(27*(1UL<<46)+1)
#define FFT_ROOT_2		67
#define FFT_H_2			46
#define FFT_MAX_N			28
#define FFT_DOUBLE_MAXINT	1266637395197952UL	// FFT_PRIME_1/3

/*
small values for testing

#define FFT_PRIME_1		(15*(1UL<<10)+1)
#define FFT_ROOT_1		84
#define FFT_PRIME_2		(13*(1UL<<10)+1)
#define FFT_ROOT_2		7
#define FFT_H				10
#define FFT_SINGLE_MAXINT	100			// floor(sqrt(FFT_PRIME_1))/3
#define FFT_MAX_N			10
*/


typedef struct {
	ff_t half[2][FFT_MAX_N];
	ff_t *w[2];
	ff_t *wi[2];
	ff_t *t;					// workspace
	unsigned long p[2];
	unsigned long a[2];		// a[i]=1/p[1-i] mod p[i]
	double y[2];				// y[i] = (double)a[i]/p[i]
	int n;
	int primes;
} fft_ctx_t[1];

int fft_precompute (fft_ctx_t ctx, int n, unsigned long maxint);

void fft_poly_mult_z (ff_t h[], ff_t f[], ff_t g[], int n, int i, fft_ctx_t ctx);
void fft_poly_mult_ff (ff_t h[], ff_t f[], ff_t g[], int n, fft_ctx_t ctx);

void ff_poly_fft (ff_t a[], ff_t w[], int n, int mult);
void ff_poly_fft_mult (ff_t h[], ff_t f[], ff_t g[], ff_t w[], ff_t wi[], int n, ff_t N_inv);

/*
// this is obsolete test code, use ff_poly_from_roots above
void ff_poly_conv_mult (ff_t h[], ff_t f[], ff_t g[], int n, fft_ctx_t *ctx);
void ff_poly_Karatsuba_mult (ff_t h[], ff_t f[], ff_t g[], ff_t t[], int n);
void ff_poly_Karatsuba_from_roots (ff_t f[], ff_t r[], ff_t t[], int d);
void ff_poly_conv_from_roots (ff_t f[], ff_t r[], ff_t t[], int d, fft_ctx_t ctx);
*/

// a and b cannot overlap
static inline void ff_poly_fft_rev (ff_t a[], ff_t b[], int n)
{
	register uint32_t i,j;
	
	for ( i = 0 ; i < (1<<n) ; i++ ) { j = ui32_revbit(i)>>(32-n); _ff_set(a[i],b[j]); }
}

#ifdef __cplusplus
}
#endif

#endif
