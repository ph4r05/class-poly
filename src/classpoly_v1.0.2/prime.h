#ifndef _PRIME_INCLUDE_
#define _PRIME_INCLUDE_

/*
	Copyright (c) 2007-2014 Andrew V. Sutherland
	See LICENSE file for license details
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PRIME_BIG			1			// make nonzero to support prime enumeration to 2^43, otherwise default is 2^39

#if PRIME_BIG

#define PRIME_SMALL_BITS				22
#define PRIME_SMALL_PRIMES				295947				// pi(2^PRIME_SMALL_BITS)
#define PRIME_MAX_SMALL_PRIME			4194301				// largest prime less than 2^PRIME_SMALL_BITS
#define PRIME_MAX_PRIME				17592186044399L		// largest prime less than 2^(2*PRIME_SMALL_BITS)

#else

#define PRIME_SMALL_BITS				20
#define PRIME_SMALL_PRIMES				82025
#define PRIME_MAX_SMALL_PRIME			1048573
#define PRIME_MAX_PRIME				1099511627689L

#endif

#define PRIME_BITS						(2*PRIME_SMALL_BITS)
#define PRIME_MAX_SMALL_INTEGER		((1<<PRIME_SMALL_BITS)-1)
#define PRIME_MAX_ENUM				(1L<<PRIME_BITS)

#define PRIME_MAX_PRIMORIAL_W			15							// largest w such that P_w = 2*3*...*p_w fits in a (64-bit) long
#define PRIME_SMALL_PRIMORIAL_W		6							// P_6 = 2*3*5*7*11*13 = 30030
#define PRIME_SMALL_PRIMORIAL			(2*3*5*7*11*13)				// must be less than MPZ_MAX_SMALL_INTEGER
#define PRIME_MAX_WHEEL_W				9							// needs 200 mb map to construct a wheel this large

struct wheel_struct {
	unsigned long n;
	unsigned long phi;
	unsigned char *gaps;
	unsigned maxgap;
};
typedef struct wheel_struct wheel_t;

struct prime_enum_ctx_struct {
	long L, end;											// usually L = end, but for windowed enumeration L may need to be greater than end
	long s[PRIME_BITS];										// for i > 1, s[i] is the index of the greatest prime q st q^i < p
	long e[PRIME_BITS];									// for i > 1, e[i] is the index of the least prime q st q^i >= L
	long q[PRIME_BITS];									// for i > 1 q[i] = s^i
	long h;
	long w;
	long *wi;
	long *wv;
	long p, bufp;
	long pi;
	long j;
	long start;
	wheel_t *wheel;
	unsigned char *map;
	unsigned char *mp;
	int powers;
	int k;												// if powers > 1, then q[k] is the minimum q[i] for 2 <= i <= powers
	int d;												// if powers > 1, then the most recently enumerated prime power is b^d
	long b;
	// entries below are for lookahead in prime_enum_w
	unsigned char *la_gaps;
	int la_i, la_j, la_n, la_window;
	long la_p, la_ep, la_end;
};
typedef struct prime_enum_ctx_struct prime_enum_ctx_t;

// prime_setup will be automatically called if/when needed
extern int _prime_inited;
void _prime_setup();
static inline void prime_setup() { if ( ! _prime_inited ) _prime_setup(); }
int *prime_small_primes ();
void prime_cleanup();

prime_enum_ctx_t *prime_enum_start_w (long start, long end, long window);	// enumerates primes in [start,end],  and optionally logs primes in window
																// use window < 0 for trailing window, window >0 centered, window=0 for no window
prime_enum_ctx_t *prime_enum_start (long start, long end, int powers);		// enumerates primes (or prime powers p^k with k <= powers) in [start,end]
long _prime_enum (prime_enum_ctx_t *ctx);
long _prime_enum_w (prime_enum_ctx_t *ctx);
static inline long prime_enum (prime_enum_ctx_t *ctx) { return ctx->la_window ? _prime_enum_w(ctx) : _prime_enum(ctx); }
long prime_enum_powers (prime_enum_ctx_t *ctx);
static inline int prime_enum_exp (prime_enum_ctx_t *ctx)					// if powers is set, this returns the exponent k of the most recently enumerated prime power p^k
	{ return ctx->d; }
static inline int prime_enum_base (prime_enum_ctx_t *ctx)					// if powers is set, this returns the base p of the most recently enumerated prime power p^k if it was not prime (0 ow)
	{ return ctx->d > 1 ? ctx->b : 0; }
long prime_power_enum (prime_enum_ctx_t *ctx);
void prime_enum_end (prime_enum_ctx_t *ctx);

void prime_start_logging (long start, int window);							// explictly turn on logging to keep a record of primes enumerated within a trailing window
int prime_check_log (long n);
void prime_stop_logging (void);
	
long prime_pi (long x);												// enumerates primes up to x, use prime_table_pi if you will be making many calls to prime_pi

int is_small_prime (int p);
int nth_small_prime (int n);											// returns the nth prime for 0 < n <= PRIME_SMALL_PRIMES (2 is the 1st prime)
int small_prime_pi (int p);											// returns pi(n) for n <= PRIME_MAX_SMALL_INTEGER
long pi_bound (long n);

long primorial (int w);
long primorial_phi (int w);

long next_prime (long p);											// for large primes just calls GMP, create a prime table if you care about speed

#define PRIME_TABLE_MAX_LEVELS			16
#define PRIME_TABLE_LG_PRIMES_PER_NODE	3
#define PRIME_TABLE_PRIMES_PER_NODE		(1<<PRIME_TABLE_LG_PRIMES_PER_NODE)
#define PRIME_TABLE_LG_GAPS_PER_LEAF	6
#define PRIME_TABLE_GAPS_PER_LEAF		(1<<PRIME_TABLE_LG_GAPS_PER_LEAF)
#define PRIME_TABLE_MAX_PRIME			304599508537L				// gaps up to PRIME_TABLE_MAX_PRIME are all <= 500
#define PRIME_TABLE_MAX_PI				11992433550L
	
struct gaptab_struct {
	long si, ti;
	long sp, tp;
	unsigned char *gaps;
	long *levels[PRIME_TABLE_MAX_LEVELS];
	int top;
};
typedef struct gaptab_struct gap_table_t[1], prime_table_t[1];

void gap_table_init (prime_table_t tab, long off, unsigned char *gaps, long n);	// gaps must be 0-terminated and are automatically expanded by a factor of 2 (so largest gap possible is 510)
void prime_table_init (prime_table_t tab, long start_pi, long end_pi);
void gap_table_clear (gap_table_t tab);
static inline void prime_table_clear (prime_table_t tab) { gap_table_clear (tab); }

// returns the least integer p >= x that lies in the table (0 if none) and sets *pgap to point to gap between p and the next entry in the table (0 if last entry)
static inline long gap_table_next_gap (gap_table_t tab, unsigned char **pgap, long x)
{
	register long *pp, p, j;
	register int i;
	
	if ( x > tab->tp ) { *pgap = 0; return 0; }
	if ( x <= tab->sp ) { *pgap = tab->gaps; return tab->sp; }
	pp = tab->levels[i=tab->top];
	for(;;) {
		while ( *(pp+1) <= x ) pp++;
		if ( ! i ) break;
		pp = tab->levels[i-1] + ((pp-tab->levels[i]) << PRIME_TABLE_LG_PRIMES_PER_NODE);
		i--;		
	}
	p = *pp;
	j = (pp-tab->levels[0]) << PRIME_TABLE_LG_GAPS_PER_LEAF;
	while ( p < x ) p += 2*tab->gaps[j++];
	*pgap = tab->gaps+j;
	return p;
}
static inline long prime_table_next_gap (prime_table_t tab, unsigned char **pgap, long x)
	{ return gap_table_next_gap (tab, pgap, x); }

static inline long prime_table_gap_pi (prime_table_t tab, unsigned char *pgap)
	{ long n = tab->si+pgap-tab->gaps;  if ( n < tab->si || n > tab->ti ) return 0; else return n; }

static inline long gap_table_nth_number (gap_table_t tab, long n)
{
	register long i, p;
	
	if ( n < tab->si || n > tab->ti ) return 0;
	n -= tab->si;
	i = n>>PRIME_TABLE_LG_GAPS_PER_LEAF;
	p = tab->levels[0][i];
	for ( i <<= PRIME_TABLE_LG_GAPS_PER_LEAF ; i < n ; i++ ) p += 2*tab->gaps[i];
	return p;
}
static inline long prime_table_nth_prime (prime_table_t tab, long n)
	{ return gap_table_nth_number (tab, n); }
	
// returns the nth entry p and sets *pgap to point to gap between p and the entry after p (0 if last entry)
static inline long gap_table_nth_gap (gap_table_t tab, unsigned char **pgap, long n)
{
	register long i, p;
	
	if ( n < tab->si || n > tab->ti ) { *pgap = 0;  return 0; }
	n -= tab->si;
	i = n>>PRIME_TABLE_LG_GAPS_PER_LEAF;
	p = tab->levels[0][i];
	for ( i <<= PRIME_TABLE_LG_GAPS_PER_LEAF ; i < n ; i++ ) p += 2*tab->gaps[i];
	*pgap = tab->gaps+n;
	return p;
}
static inline long prime_table_nth_gap (prime_table_t tab, unsigned char **pgap, long n)
	{ return gap_table_nth_gap (tab, pgap, n); }
	
// returns a count of the number of integers in the gap talbe that are <= x
static inline long gap_table_count (gap_table_t tab, long x)
{
	unsigned char *gap;
	unsigned long p;
	
	p = gap_table_next_gap (tab, &gap, x);
	if ( ! p ) return 0;
	return tab->si + gap-tab->gaps - ( p > x ? 1 : 0 );
}
static inline long prime_table_pi (prime_table_t tab, long x)
	{ return gap_table_count (tab, x); }
		
#ifdef __cplusplus
}
#endif
	
#endif
