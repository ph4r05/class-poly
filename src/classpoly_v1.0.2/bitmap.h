#ifndef _INCLUDE_BITMAP_
#define _INCLUDE_BITMAP_

/*
	Copyright (c) 2009, 2012, 2014 Andrew V. Sutherland
	See LICENSE file for license details.
*/


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <memory.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef uint64_t *bitmap_t;

static inline long bm_w64 (long n) { return (n>>6)+((n&0x3F)?1:0); }

// fast inline bitmap implementation -- define BITMAP_NOCHECK to disable range checking (only makes a slight difference)

static uint64_t bitmap_mask_tab[64] = { 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80,
	0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000,
	0x10000, 0x20000, 0x40000, 0x80000, 0x100000, 0x200000, 0x400000, 0x800000,
	0x1000000, 0x2000000, 0x4000000, 0x8000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000,
	0x100000000, 0x200000000, 0x400000000, 0x800000000, 0x1000000000, 0x2000000000, 0x4000000000, 0x8000000000,
	0x10000000000, 0x20000000000, 0x40000000000, 0x80000000000, 0x100000000000, 0x200000000000, 0x400000000000, 0x800000000000,
	0x1000000000000, 0x2000000000000, 0x4000000000000, 0x8000000000000, 0x10000000000000, 0x20000000000000, 0x40000000000000, 0x80000000000000,
	0x100000000000000, 0x200000000000000, 0x400000000000000, 0x800000000000000, 0x1000000000000000, 0x2000000000000000, 0x4000000000000000, 0x8000000000000000 };

static char bitmap_byte_weights[256] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
								  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                                                                  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                                                                  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};
	
// it is convenient to allocate an extra zero longword past the end of our bitmaps (this is used by bm_residue_set)
static inline uint64_t *_bm_alloc (long n) {n = bm_w64(n)+1; uint64_t *bm = malloc(n*sizeof(uint64_t)); bm[n-1] = 0;  bm[n-2] = 0; return bm; }	// important to clear trailing bits even though rest is uninitialized
static inline uint64_t *bm_alloc(long n) { n = bm_w64(n)+1; uint64_t *bm = malloc(n*sizeof(uint64_t)); for ( register long i = 0 ; i < n ; i++ ) bm[i] = 0;  return bm; }
static inline void bm_free(uint64_t *bm) { free(bm); }
static inline void bm_clear(uint64_t *bm, long n) { register long i; for ( i = 0 ; i < (n>>6) ; i++ ) bm[i] = 0;  if ( (n&0x3F) ) bm[i] &= ~(bitmap_mask_tab[n&0x3f]-1); }
static inline void bm_setall(uint64_t *bm, long n) { register long i; for ( i = 0 ; i < (n>>6) ; i++ ) bm[i] = 0xFFFFFFFFFFFFFFFFUL; if ( (n&0x3F) ) bm[i] |= bitmap_mask_tab[n&0x3f]-1; }
static inline uint64_t *bm_alloc_set(long n) { uint64_t *bm = _bm_alloc (n); bm_setall (bm, n); return bm; }
static inline void bm_copy(uint64_t *bm2, uint64_t *bm1, long n) { register long i; n = bm_w64(n)+1;  for ( i = 0 ; i < n ; i++ ) bm2[i] = bm1[i]; }
static inline void bm_set(uint64_t *bm, long i) { bm[i>>6] |= bitmap_mask_tab[i&0x3f]; }
static inline void bm_unset(uint64_t *bm, long i) { bm[i>>6] &= ~bitmap_mask_tab[i&0x3f]; }
static inline int bm_test(uint64_t *bm, long i) { return ((bm[i>>6] & bitmap_mask_tab[i&0x3f])?1:0); }
static inline int bm_test_and_set(uint64_t *bm, long i) { register uint64_t x = bm[i>>6];  register uint64_t y = x | bitmap_mask_tab[i&0x3f];  bm[i>>6] = y;  return y==x; }
static inline void bm_complement (uint64_t *bm, long n) { register long i, j; j=bm_w64(n); for ( i = 0 ; i < j ; i++ ) bm[i] = ~bm[i];  if ( (n&0x3f) ) bm[n>>6] &= bitmap_mask_tab[n&0x3f]-1; }
static inline int bm_allset(uint64_t *bm, long n)	// returns 1 if bits 0..n-1 are all set, 0 otherwise
{
	register long i;
	uint64_t m;
	
	for ( i = 0 ; i < (n>>6) ; i++ ) if ( bm[i] != 0xFFFFFFFFFFFFFFFFUL ) return 0;
	m = bitmap_mask_tab[n&0x3f]-1;
	if ( (bm[i]&m) != m ) return 0;
	return 1;
}

static inline long bm_min_set(uint64_t *bm, long n)	// returns the index of the least set bit with index < n, or n if none
{
	register long i, j, k;
	
	for ( i = 0 ; i < (n>>6) && !bm[i] ; i++ );
	k = ( i < (n>>6 ) ? 64 : (n&0x3f) );
	for ( j = 0 ; j < k ; j++ ) if ( (bm[i] & bitmap_mask_tab[j]) ) break;
	return (i<<6)+j;
}

static inline long bm_min_unset(uint64_t *bm, long n)	// returns the index of the least unset bit with index < n, or n if none
{
	register long i, j, k;
	
	for ( i = 0 ; i < (n>>6) && bm[i] == 0xFFFFFFFFFFFFFFFFUL ; i++ );
	k = ( i < (n>>6 ) ? 64 : (n&0x3f) );
	for ( j = 0 ; j < k ; j++ ) if ( !(bm[i] & bitmap_mask_tab[j]) ) break;
	return (i<<6)+j;
}

static inline long bm_next_set(uint64_t *bm, long m, long n)	// returns the index of the least set bit with index in [m,n), or n if none
{
	register long i, j, k, imax;
	
	if ( m >= n ) return n;
	i = m>>6;  imax = n>>6;
	k = ( i < imax ? 64 : (n&0x3f) );
	for ( j=(m&0x3f) ; j < k ; j++ ) if ( (bm[i]&bitmap_mask_tab[j]) ) break;
	if ( j < k ) return (i<<6)+j;
	if ( i == imax ) return n;
	for ( i++ ; i < imax && !bm[i]  ; i++ );
	k = ( i < imax ? 64 : (n&0x3f) );
	for ( j = 0 ; j < k ; j++ ) if ( (bm[i]&bitmap_mask_tab[j]) ) break;
	return (i<<6)+j;
}
static inline long bm_next_unset(uint64_t *bm, long m, long n)	// returns the index of the least unset bit with index in [m,n), or n if none
{
	register long i, j, k, imax;
	
	if ( m >= n ) return n;
	i = m>>6;  imax = n>>6;
	k = ( i < imax ? 64 : (n&0x3f) );
	for ( j=(m&0x3f) ; j < k ; j++ ) if ( ! (bm[i]&bitmap_mask_tab[j]) ) break;
	if ( j < k ) return (i<<6)+j;
	if ( i == imax ) return n;
	for ( i++ ; i < imax && bm[i] == 0xFFFFFFFFFFFFFFFFUL ; i++ );
	k = ( i < imax ? 64 : (n&0x3f) );
	for ( j = 0 ; j < k ; j++ ) if ( ! (bm[i]&bitmap_mask_tab[j]) ) break;
	return (i<<6)+j;
}

// xor 64*w bits from bm+off to cm
static inline void bm_xor_bits (uint64_t *cm, uint64_t *bm, long off, long w)
{
	register int a,aa;
	register long i,j;
	
	a = off&0x3f;
	j = off>>6;
	if ( !a ) {
		for ( i = 0 ; i < w ; i++, j++ ) cm[i] |= bm[j];
	} else {
		aa = 64-a;
		for ( i = 0 ; i < w ; i++, j++ ) cm[i] |= (bm[j]>>a) | (bm[j+1]<<aa);
	}
}

// sets ith bit in rm to 1 if there is a bit k=i mod m set in bm.  relies on bm being padded with 0 bits (plus a zero longword)
static inline void bm_residue_set (uint64_t *rm, long m, uint64_t *bm, long n)
{
	uint64_t mask;
	register long off, w;

	w = bm_w64(m);
	for ( off = 0 ; off < w ; off++ ) rm[off] = 0;
	for ( off = 0 ; off < n-m ; off += m ) bm_xor_bits (rm, bm, off, w);
	bm_xor_bits (rm, bm, off, bm_w64(n-off));		// note that here we potentially rely on there being an extra zero longword padding bm
	mask = (1UL<<(m&0x3f)) - 1UL;
	if ( mask ) rm[w-1] &= mask;					// clear high bits past m in last longword
}

// sets ith bit in rm to 1 if there is a bit k=i mod m set in bm with i <= 64*w < m.  relies on bm being padded with 0 bits (plus a zero longword)
static inline void bm_partial_residue_set (uint64_t *rm, long w, long m, uint64_t *bm, long n)
{
	long k;
	register long off;

	for ( off = 0 ; off < w ; off++ ) rm[off] = 0;
	for ( off = 0 ; off < n-m ; off += m ) bm_xor_bits (rm, bm, off, w);
	k = bm_w64(n-off);  if ( k < w ) w = k;
	bm_xor_bits (rm, bm, off, w);					// note that here we potentially rely on there being an extra zero longword padding bm
}

// returns the number of bits in [0..n-1] that are set
static inline long bm_weight (uint64_t *bm, long n)
{
	unsigned char mask;
	unsigned char *s, *e;
	long w;

	s = (unsigned char *)bm;
	e = s + (n>>3);
	for ( w = 0 ; s < e ; w += bitmap_byte_weights[*s++] );
	mask = (1<<(n&0x7)) - 1;
	w += bitmap_byte_weights[*s&mask];
	return w;
}

// returns true if the set defined by the set bits of A is contined in the set defined by the set bits of B
static inline int bm_subset (uint64_t *A, uint64_t *B, long n)
{
	long w = bm_w64(n);
	for ( long i = 0 ; i < w ; i++ ) if ( A[i] & ~B[i] ) return 0;
	return 1;
}

// The default bitmap array is not shared, each module gets its own static copy that it can reuse (use bm routines if you want finer control)
static uint64_t *bitmap;
static long bitmap_words;

#ifdef BITMAP_NOCHECK
static inline void bitmap_check (long i) {}
#else
static inline void bitmap_check (long i) { if ( i < 0 || (i>>6) >= bitmap_words ) { printf ("Range check failed for index %ld with bitmap of length %ld, program aborting...\n", i, bitmap_words<<6); abort(); } }
#endif
// pad bitmap with extra zero longword
static inline void bitmap_alloc(long n) { register long i;  bitmap_words = bm_w64(n); bitmap = realloc(bitmap,(bitmap_words+1)*sizeof(uint64_t)); for ( i = 0 ; i <= bitmap_words ; i++ ) bitmap[i] = 0; }
static inline void bitmap_free(void) { bm_free(bitmap); bitmap = 0; bitmap_words = 0; }
static inline void bitmap_clear(long n) { if ( n) bitmap_check(n-1); bm_clear(bitmap, n); }
static inline void bitmap_set(long i) { bitmap_check(i); bm_set(bitmap,i); }
static inline void bitmap_unset(long i) { bitmap_check(i); bm_unset(bitmap,i); }
static inline int bitmap_test(long i) { bitmap_check(i); return bm_test(bitmap,i); }
static inline int bitmap_test_and_set(long i) { bitmap_check(i); return bm_test_and_set(bitmap,i); }
static inline int bitmap_allset(long n) { if ( n) bitmap_check(n-1); return bm_allset(bitmap, n); }
static inline long bitmap_min_unset(long n) { if ( n) bitmap_check(n-1); return bm_min_unset(bitmap, n); }
static inline long bitmap_next_unset(long m, long n) { if ( n) bitmap_check(n-1); return bm_next_unset(bitmap, m, n); }

#ifdef __cplusplus
}
#endif

#endif
