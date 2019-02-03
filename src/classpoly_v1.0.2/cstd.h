#ifndef _CSTD_INCLUDE_			// handles circular inclusions
#define _CSTD_INCLUDE_

#ifdef __cplusplus
extern "C" {
#endif

/*
	Copyright (c) 2007-2014 Andrew V. Sutherland
	See LICENSE file for license details.
*/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <memory.h>
#include <math.h>

#define LONG_BITS		((int)sizeof(long)*CHAR_BIT)

/*
	Header file for general purpose definitions and library extensions
	that we want to include everywhere.
*/

#define _swap(a,b,w)		(w=a,a=b,b=w)

#define DEBUG_LEVEL		2
#define INFO_LEVEL			1
#define WARN_LEVEL		-1
#define ERROR_LEVEL		-2
#define OUTPUT_LEVEL		0

int dbg_level;		// global shared variable to set debug/error reporting level

/*
	Output macros for error and debug reporting.
	These are potentially dangerous, should switch to inlines using variable length argument lists!
*/

#define dbg_setlevel(i)		(dbg_level = (i))
#define dbg_printf			if ( dbg_level >= DEBUG_LEVEL ) printf 
#define info_printf			if ( dbg_level >= INFO_LEVEL ) printf
#define warn_printf			if ( dbg_level >= WARN_LEVEL ) printf
#define err_printf			if ( dbg_level >= ERROR_LEVEL ) printf
#define out_printf			if ( dbg_level >= OUTPUT_LEVEL ) printf

#define gmp_dbg_printf		if ( dbg_level >= DEBUG_LEVEL ) gmp_printf
#define gmp_warn_printf		if ( dbg_level >= WARN_LEVEL ) gmp_printf
#define gmp_info_printf		if ( dbg_level >= INFO_LEVEL ) gmp_printf
#define gmp_err_printf		if ( dbg_level >= ERROR_LEVEL ) gmp_printf
#define gmp_out_printf		if ( dbg_level >= OUTPUT_LEVEL ) gmp_printf

#define delta_msecs(s,t)		(1000UL*(t-s)/CLOCKS_PER_SEC)
#define delta_nsecs(s,t)		(1000000000UL*(t-s)/CLOCKS_PER_SEC)			// assumes 64 bit UL
#define delta_wall_msecs(w1,w2)	\
	(1000*((w2)->tv_sec - (w1)->tv_sec) + ((w2)->tv_usec - (w1)->tv_usec)/1000)

// this function should be used by anyone wanting a seed for a random number generator, so that when debugging you can fix the seed to get deterministic behavior
unsigned long _cstd_seed;		// you can set this to any fixed non-zero value (in your main function, not here)
static inline unsigned long cstd_seed(void)
{
	if ( ! _cstd_seed ) {
		_cstd_seed = (((unsigned long)gethostid())<<32) + getpid();		// make sure seed is different in different processes/hosts, even if seeded at (approximately) the same time
		_cstd_seed *= time(0);
	}
	return _cstd_seed;
}

#define LOG2			0.69314718055994530941723212145818
//#define log2(x)		(log(x)/LOG2)			no longer necessary with gcc 4.4.3

#define _min(a,b)		((a)<(b)?(a):(b))		// generic min, not safe with side-effects

// Some handy utility functions
static inline unsigned long _asm_highbit (unsigned long x) { asm ("bsrq %0, %0" : "=r" (x) : "0" (x)); return x; }
static inline unsigned long _asm_lowbit (unsigned long x) { asm ("bsfq %0, %0" : "=r" (x) : "0" (x)); return x; }
static inline unsigned long ui_ceil_ratio (unsigned long a, unsigned long b)
	{ register unsigned long x = a/b;  return ( x*b<a ? x+1 : x ); }
static inline int ui_len (unsigned long x) { return LONG_BITS-__builtin_clzl(x); }
static inline unsigned long ui_min (unsigned long a, unsigned long b)
	{ return ( a < b ? a : b ); }
static inline unsigned long ui_max (unsigned long a, unsigned long b)
	{ return ( a > b ? a : b ); }

// this is about twice as fast as x%m on an AMD Athlon 64, but it is not valid for large x (unless you make mi a long double, but then things slow down)  -- use with caution!
static inline unsigned long ui_mod_i (unsigned long x, unsigned long m, double mi)
	{ register unsigned long t, z;  t = mi * x - 0.5;  z = x-t*m;  if ( z >= m ) z-= m; return z; }

static inline unsigned long ui_mod (unsigned long x, unsigned long m) { return x%m; }	// play it safe
static inline long i_mod(long x, long m) { register long t;  if ( x >= 0 ) return ui_mod(x,m);  t = - (long)ui_mod((unsigned long)(-x),m);  return (t<0?t+m:t); }

static inline uint32_t ui32_revbit(uint32_t x)
{ 
	register uint32_t n = x;
	n = ((n >> 1) & 0x55555555) | ((n << 1) & 0xaaaaaaaa);
	n = ((n >> 2) & 0x33333333) | ((n << 2) & 0xcccccccc);
	n = ((n >> 4) & 0x0f0f0f0f) | ((n << 4) & 0xf0f0f0f0);
	n = ((n >> 8) & 0x00ff00ff) | ((n << 8) & 0xff00ff00);
	n = ((n >> 16) & 0x0000ffff) | ((n << 16) & 0xffff0000);
	return n;
}

static inline uint64_t ui64_revbit(uint64_t x)
{ 
	register uint64_t n = x;
	n = ((n >> 1) & 0x5555555555555555) | ((n << 1) & 0xaaaaaaaaaaaaaaaa);
	n = ((n >> 2) & 0x3333333333333333) | ((n << 2) & 0xcccccccccccccccc);
	n = ((n >> 4) & 0x0f0f0f0f0f0f0f0f) | ((n << 4) & 0xf0f0f0f0f0f0f0f0);
	n = ((n >> 8) & 0x00ff00ff00ff00ff) | ((n << 8) & 0xff00ff00ff00ff00);
	n = ((n >> 16) & 0x0000ffff0000ffff) | ((n << 16) & 0xffff0000ffff0000);
	n = ((n >> 32) & 0x0000ffffffff) | ((n << 32) & 0xffffffff00000000);
	return n;
}

// mem_alloc guarantees success - it will abort if out of memory, zero fills memory
static inline void *mem_alloc (unsigned long bytes)
{
	void *ptr;

	if ( bytes < 4 ) fprintf (stderr, "Warning mem_alloc size %lu < 4 - ptr error?!\n", bytes);
	ptr = malloc (bytes);
	if ( ! ptr ) { fprintf (stderr, "Fatal error, attempted memory allocation of %lu bytes failed.\n", bytes);  abort(); }
	dbg_printf ("Allocated %lu bytes at %lx\n", bytes, (unsigned long)ptr);
	memset (ptr, 0, bytes);
	return ptr;
}

static inline void mem_free (void *ptr)
{
	dbg_printf ("freed %lx\n", (unsigned long)ptr); 
	free (ptr);
}

// handle simple exponentials in string to integer conversions
static inline long atol_exp (char *str)
{
	register char *s;
	long i,b,n,N;
	
	for ( s = str ; *s && *s != 'e' && *s != 'E' ; s++ );
	if ( *s ) { s++; N = b = atol (str); n = atol(s); for ( i = 1; i < n ; i++ )  N *= b; } else N = atol(str);
	return N;
}

// Functions to enumerate all t-tuples of integers in [0..n-1] in lexicographic order c[t],c[t-1],...,c[1], using Knuth Alg 7.2.1.3T
// c[0] holds an auxiliary variable, combination is c[1],c[2],...,c[t], where 0 < t < n, space for 2 sentinels is also required, so c must have t+3 entries allocated
static inline void lex_combo_first (int c[], int n, int t)
	{ register int j;  for ( j = 1 ; j <= t ; j++ ) c[j] = j-1;  c[t+1] = n;  c[t+2] = 0;  c[0] = t; }
		
static inline int lex_combo_next (int c[], int n, int t)
{
	register int j;

	if ( c[0] ) { c[c[0]] = c[0]; c[0]--; return 1; }
	if ( c[1] +1 < c[2] ) { c[1]++; return 1; }
	j = 2;  c[1] = 0;
	while ( c[j]+1 == c[j+1] ) { j++;  c[j-1] = j-2; }
	if ( j > t ) return 0;
	c[j]++;  c[0] = j-1;
	return 1;
}

#ifdef __cplusplus
}
#endif

#endif
