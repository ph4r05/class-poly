#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <stdint.h>
#include <gmp.h>
#include "prime.h"
#include "bitmap.h"

/*
	Copyright 2007-2014 Andrew V. Sutherland
	See LICENSE file for details
*/

#define LONG_BITS		((int)sizeof(long)*CHAR_BIT)

// Some handy inlines
static inline unsigned long ui_ceil_ratio (unsigned long a, unsigned long b)
	{ register unsigned long x = a/b;  return ( x*b<a ? x+1 : x ); }
static inline int ui_len (unsigned long x) { return LONG_BITS-__builtin_clzl(x); }
static inline unsigned long ui_min (unsigned long a, unsigned long b)
	{ return ( a < b ? a : b ); }
static long ipow ( long a,  long e)
	{ register long b, i;  for ( b=1, i=0 ; i < e ; b*=a, i++ ); return b; }
static inline void *_calloc (size_t size) { return calloc (size,1); }
	
int _prime_inited;

static int small_primes[PRIME_SMALL_PRIMES+1];
static int small_prime_index[PRIME_MAX_SMALL_INTEGER+1];
static short small_primorial_map[PRIME_SMALL_PRIMORIAL];		// i-th entry is 0 if gcd(i,PRIME_SMALL_PRIMORIAL)>1, otherwise it indicates that i is the k-th integer in F_P*, where 1 is the 1st.

static long _primorials[PRIME_MAX_PRIMORIAL_W+1];
static long _primorial_phis[PRIME_MAX_PRIMORIAL_W+1];
wheel_t _small_wheels[PRIME_SMALL_PRIMORIAL_W+1];

static unsigned char _wheel_gaps0[1] = { 1 };		// The trivial wheel - every number is relatively prime to 1
static unsigned char _wheel_gaps1[2] = { 2 };		// The first wheel - odd numbers

static long prime_logsize, prime_logstart, prime_logindex, prime_logvalue;

void prime_cleanup()
{
	register int i;
	
	if ( ! _prime_inited ) return;
	for ( i = 2 ; i <= PRIME_SMALL_PRIMORIAL_W ; i++ ) free (_small_wheels[i].gaps);
	_prime_inited = 0;
}

int *prime_small_primes ()
{
	prime_setup();
	return small_primes;
}

void _prime_setup ()
{
	register unsigned char *mp, *np;
	register unsigned char *_small_map;
	register int  i, j, k, p, maxp, w, gap, maxgap;

	if ( _prime_inited ) return;
	
	assert (LONG_BITS >= 64);									// life is too short to support 32-bit code

	_primorials[0] = 1;
	_primorial_phis[0] = 1;
	_primorials[1] = 2;
	_primorial_phis[1] = 1;
	_small_wheels[0].n = 1;
	_small_wheels[0].phi = 1;
	_small_wheels[0].gaps = _wheel_gaps0;
	_small_wheels[1].n = 2;
	_small_wheels[1].phi = 1;
	_small_wheels[1].gaps = _wheel_gaps1;
	_small_map = (unsigned char *) _calloc (PRIME_MAX_SMALL_INTEGER+1);
	*(_small_map+PRIME_MAX_SMALL_INTEGER) = 1;					// mark the end of the map to avoid rolling past it
	_small_map[1] = 1;											// roll the first prime to get things started
	_small_map[3] = 1;
	small_primes[1] = 2;
	for ( i = 2 ; i <= PRIME_SMALL_PRIMORIAL_W ; i++ ) {
		for ( mp = _small_map+2 ; ! *mp ; mp++ );					// find first marked entry > 1 in previous map
		*mp = 0;												// must be the next prime.  clear it.
		p = mp-_small_map;
		small_primes[i] = p;
		_small_wheels[i].n = _primorials[i] = _primorials[i-1]*p;
		_small_wheels[i].phi = _primorial_phis[i] = _primorial_phis[i-1]*(p-1);
		mp = _small_map+_primorials[i-1] +1;
		for ( k = 1 ; k < p ; k++ ) {								// roll the previous wheel (p-1) more times
			for ( j = 0 ; j < _primorial_phis[i-1] ; j++ ) { *mp = 1;  mp += _small_wheels[i-1].gaps[j]; }
		}
		if ( mp != _small_map+_primorials[i]+1 ) { fprintf (stderr, "mpz_util_init: primorial wheel alignment error\n");  abort(); }
		*mp = 1;												// mark primorial+1 entry to bound last gap
		np = _small_map+_primorials[i];
		for ( mp = _small_map+p ; mp < np ; mp += p ) *mp = 0;		// clear the multiples of p
		_small_wheels[i].gaps = (unsigned char *) _calloc(_primorial_phis[i]+4);	// make alloc big enough to avoid small alloc warning
		maxgap = 0;
		for ( mp = _small_map+1, j = 0 ; j < _primorial_phis[i] ; j++ ) {
			for ( np = mp+1 ; ! *np ; np++ );
			gap = np-mp;
			if ( gap > maxgap ) {
				maxgap = gap;
				if ( gap > UCHAR_MAX ) { fprintf (stderr, "gap %d too large to store in primorial wheel %d\n", gap, i);  abort(); }
			}
			_small_wheels[i].gaps[j] = gap;
			mp = np;
		}
		if ( np != _small_map+_primorials[i]+1 ) { fprintf (stderr, "mpz_util_init: primorial wheel alignment error\n");  abort(); }
		_small_wheels[i].maxgap = maxgap;
		//printf ("Wheel %d, maxgap %d\n", i, maxgap);
	}
	w = --i;
	// roll last wheel across rest of the small map
	mp = _small_map+_primorials[w]+1;
	np = _small_map+PRIME_MAX_SMALL_INTEGER;
	while ( mp < np ) {
		for ( j = 0 ; j < _primorial_phis[w] ; j++ ) {
			mp += _small_wheels[w].gaps[j];
			if ( mp >= np ) break;
			*mp = 1;
		}
	}
	// now sieve for primes up to PRIME_MAX_SMALL_INTEGER
	maxp = (unsigned) floor(sqrt((double)PRIME_MAX_SMALL_INTEGER));
	mp = _small_map + small_primes[i];
	for (;;) {
		mp++;
		while ( ! *mp ) mp++;
		p = mp-_small_map;
		if ( p > maxp ) break;
		small_primes[++i] = p;
		for ( j = 0, k=1 ; ; j++ ) {
			if ( j == _small_wheels[w].phi ) j = 0;
			k += _small_wheels[w].gaps[j];
			if ( k*p > PRIME_MAX_SMALL_INTEGER ) break;
			_small_map[k*p] = 0;
		}
	}
	mp = _small_map+small_primes[i]+1;
	for (;;) {
		while ( ! *mp ) mp++;
		if ( mp >= np ) break;
		small_primes[++i] = mp-_small_map;
		mp++;
	}
	if ( i != PRIME_SMALL_PRIMES || small_primes[i] != PRIME_MAX_SMALL_PRIME ) { fprintf (stderr, "Error sieving small primes\n");  abort(); }
	free (_small_map);
	// index primes
	for ( i = 1 ; i <= PRIME_SMALL_PRIMES ; i++ ) small_prime_index[small_primes[i]] = i;
	// spread prime indexes to get pi(n) table
	for ( i = 0, j = 0 ; i <= PRIME_MAX_SMALL_INTEGER ; i++ ) {
		if ( small_prime_index[i] ) {
			j = small_prime_index[i];
		} else {
			small_prime_index[i] = j;
		}
	}

	// Compute various primorials and prime products
	for ( i = PRIME_SMALL_PRIMORIAL_W+1 ; i <= PRIME_MAX_PRIMORIAL_W ; i++ ) {
		_primorials[i] = _primorials[i-1]*small_primes[i];
		_primorial_phis[i] = _primorial_phis[i-1]*(small_primes[i]-1);
	}

	// roll small primorial wheel to create map of integers mod PRIME_MAX_SMALL_PRIMORIAL
	w = PRIME_SMALL_PRIMORIAL_W;
	for ( j = 0, k = 1; j < _small_wheels[w].phi ; j++ ) {
		small_primorial_map[k] = j+1;
		k += _small_wheels[w].gaps[j];
	}
	
	_prime_inited = 1;
}

int is_small_prime (int n)
{
	if ( n < 2 || n > PRIME_MAX_SMALL_PRIME ) return 0;
	prime_setup();
	return ( small_primes[small_prime_index[n]] == n ? 1 : 0 );
}

int nth_small_prime (int n)
{
	if ( n <= 0 ) return 0;
	if ( ! n || n > PRIME_SMALL_PRIMES ) { fprintf (stderr, "Invalid call to nth_small_prime with n = %d > PRIME_SMALL_PRIMES = %d\n", n, PRIME_SMALL_PRIMES);  abort(); }
	prime_setup();
	return small_primes[n];
}

int small_prime_pi (int x)				// returns pi(n) for n <= PRIME_MAX_SMALL_INTEGER
{
	if ( x < 2 ) return 0;
	if ( x > PRIME_MAX_SMALL_INTEGER ) { fprintf (stderr, "Invalid call to small_prime_pi with x = %d > PRIME_MAX_SMALL_INTEGER = %d\n", x, PRIME_MAX_SMALL_INTEGER);  abort(); }
	prime_setup();
	return small_prime_index[x];
}

long primorial (int w)
{
	if ( w > PRIME_MAX_PRIMORIAL_W ) { fprintf (stderr, "Requested primorial P_%d is too large\n", w);  abort(); }
	prime_setup();
	return _primorials[w];
}

long primorial_phi (int w)
{
	if ( w > PRIME_MAX_PRIMORIAL_W ) { fprintf (stderr, "Requested primorial P_%d is too large\n", w);  abort(); }
	prime_setup();
	return _primorial_phis[w];
}

// returns explicit upper bound on pi(n) = (x/log x)(1+3/(2log x))  based on Shoup p.91 (from Rosser and Schoenfeld)
long pi_bound (long n)
{
	double x,y;
	
	if ( n < 59 ) return 18;
	y = log((double)n);
	x = (double) n / y;
	x *= (1.0 + 3.0/(2.0*y));
	return ((long)ceil(x));
}


static wheel_t *_primorial_wheel_alloc (int w)
{
	wheel_t *wheel;
	char *map, *mp, *np;
	unsigned char *small_gaps, gap;
	long i, j, p, small_phi;
	
	prime_setup ();
	if ( w <= PRIME_SMALL_PRIMORIAL_W ) return _small_wheels+w;
	if ( w > PRIME_MAX_WHEEL_W ) {fprintf (stderr, "Requested wheel exceeds PRIME_MAX_WHEEL_W %d > %d\n", w, PRIME_MAX_WHEEL_W);  abort(); }
	wheel = (wheel_t*) _calloc (sizeof(*wheel));
	wheel->n = _primorials[w];
	wheel->phi = _primorial_phis[w];
	wheel->gaps = (unsigned char *) _calloc (_primorial_phis[w]);
	map = (char *) _calloc (wheel->n+1);
	np = map + wheel->n + 1;
	*np = 1;		// end marker to bound the last gap;
	small_phi = _small_wheels[PRIME_SMALL_PRIMORIAL_W].phi;
	small_gaps = _small_wheels[PRIME_SMALL_PRIMORIAL_W].gaps;
	for ( i = PRIME_SMALL_PRIMORIAL_W+1 ; i <= w ; i++ ) {
		p = small_primes[i];
		mp = map+p;
		j = 0;
		while ( mp < np ) {		// wheel the prime over the map
			*mp = 1;
			if ( j == small_phi ) j = 0;
			mp += p*small_gaps[j++];
		}
	}
	// compute the gaps by rolling the small wheel and skipping marked entries
	gap = 0;
	i = j = 0;
	mp = map+1;
	while ( mp < np ) {
		if ( j == small_phi  ) j = 0;
		mp += small_gaps[j];
		gap += small_gaps[j];
		if ( ! *mp ) {
			wheel->gaps[i++] = gap;
			if ( gap > wheel->maxgap ) wheel->maxgap = gap;
			gap = 0;
		}
		j++;
	}
	if ( mp != np) { fprintf (stderr, "wheel alignment error creating wheel for w = %d, %ld != %ld\n", w, mp-map, np-map);  abort(); }
	if ( i != wheel->phi-1 ) { fprintf (stderr, "gap count error creating wheel for w = %d, %ld != %ld\n", w, i, wheel->phi);  abort(); }
	wheel->gaps[i] = 2;	// the last gap is always 2
	return wheel;		
}


static void _primorial_wheel_free (wheel_t *wheel)
{
	if ( wheel->n <= PRIME_SMALL_PRIMORIAL ) return;
	free (wheel->gaps);
	free (wheel);
}

/*
	Fast prime enumeration for p <= L <= 2^PRIME_BITS based on a wheeled sieve.
*/
prime_enum_ctx_t *prime_enum_start (long start, long L, int powers)
{
	prime_enum_ctx_t *ctx;
	register long k, m;
	register long i, n, p;

	if ( start > L ) { fprintf (stderr, "Invalid prime enumeration, start must be less than or equal to L\n");  return 0; }
	prime_setup ();
	if ( L > (1L<<PRIME_BITS) ) { fprintf (stderr, "Attempted prime enumeration too large %ld > 2^%d\n", L, PRIME_BITS);  abort(); }
	prime_stop_logging();								// avoid any possible confusion with earlier log
	
	ctx = (prime_enum_ctx_t *) _calloc (sizeof(*ctx));
	memset (ctx, 0, sizeof(*ctx));
	p = small_primes[PRIME_SMALL_PRIMORIAL_W+1];
	ctx->L = ctx->end = L;
	ctx->h = ui_len(L);										// 2^h <= L <=2^{h+1}
	ctx->powers = ui_min(powers,ctx->h);						// indicates enumeration of prime powers p^h with h <= powers
	m = start;  if ( m < 2 ) m = 2;
	for ( i = 2 ; i <= powers || i <= 2 ; i++ ) {
		n = (long) floor(pow((double)L,1.0/(double)i));			// n^i <= L < (n+1)^i
		ctx->e[i] = small_prime_index[n];					// e(i) = pi(L^{1/i})			// we actually don't currently use anything other than e[2]
		n = (long) floor(pow((double)m,1.0/(double)i));			// n^i <= m < (n+1)^i
		// set q[i] to the least ith power of a prime (identified by its index s[i]) that is greater than or equal to m (be careful to handle the case that m is itself a prime power)
		ctx->s[i] = small_prime_index[n];
		ctx->q[i] = ipow (small_primes[ctx->s[i]], i);
		if ( ctx->q[i] < m ) {
			ctx->s[i]++;
			ctx->q[i] = ipow (small_primes[ctx->s[i]], i);
		}
	}
	ctx->k = 2;
	for ( i = 2 ; i <= ctx->powers ; i++ )
		if ( ctx->q[i] < ctx->q[ctx->k] ) ctx->k = i;				// q[k] is the least ith power of a prime greater than m, for any i in [2,powers]

	ctx->wheel = _primorial_wheel_alloc (PRIME_SMALL_PRIMORIAL_W);
	ctx->w = ctx->e[2];									// p_w <= L^{1/2} < p_{w+1} so it suffices to sieve by the first w primes
	ctx->wv = (long *) _calloc ((ctx->w+1)*sizeof(*ctx->wv));		// the entry wv[i] stores the value of the i-th prime modulo PRIME_SMALL_PRIMORIAL = ctx->wheel->n
	ctx->wi = (long *) _calloc ((ctx->w+1)*sizeof(*ctx->wi));		// the entry wi[i] stores the index of the next wheel gap for the i-th prime - see below
	for ( i =1; i <= ctx->w ; i++ ) ctx->wv[i] = small_primes[i];		// these values aren't necessarily reduced mod PRIME_SMALL_PRIMORIAL, but that's ok, the code below doesn't assume this
	ctx->map = (unsigned char*) _calloc (ctx->wheel->n);
	/*
		The only relevant entries of our map are ones relatively prime to PRIME_SMALL_PRIMORIAL, since we enumerate the map via the wheel.
		Thus we only need to sieve multiples of primes whose index lies in [PRIME_SMALL_PRIMORIAL_W+1,ctx->w] and we only need to worry
	        about those multiples which are relatively prime to PRIME_SMALL_PRIMORIAL.  Thus for each prime, we effectively use the same wheel
	        to enumerate these multiples.  The value wi[i] holds the index into the wheel for the i-th prime.
	
		If the starting point is far from 0, we can skip ahead by starting at the first multiple of PRIME_SMALL_PRIMORIAL below start, call it n.
	        For simplicity, we avoid special code for the first ctx->w primes by always forcing enumeration of the first w primes
		For each prime p <= ctx->w we need to compute the least multiple kp > start where k is relatively prime to PRIME_SMALL_PRIMORIAL
	        and then set wv[i] = kp - n and set wi[i] = index of gap between k and next integer relatively prime to PRIME_SMALL_PRIMORIAL in wheel.
	
		To keep this simple, we set m to the first multiple of p*PRIME_SMALL_PRIMORIAL below start, set wv[i] = p, wi[i] = 0 and roll forward from there.
	*/
	if ( start > 2 ) {
		mpz_t P;
		// This is a bit awkward but removes the need for any conditional code in prime_enum to handle startup.
		// We want to enumerate just up to the greatest prime below start.  To facilitate this we backup start as required
		mpz_init (P);
		n = ui_len(start);
		for ( k = 2 ;; k += n ) {
			if ( k > start ) k = start;
			mpz_set_ui (P, start-k);
			mpz_nextprime (P, P);
			if ( mpz_cmp_ui(P,start) < 0 ) break;
		}
		do {
			m = mpz_get_ui(P);
			mpz_nextprime (P, P);
		} while ( mpz_cmp_ui(P,start) < 0 );
		mpz_clear (P);
		start = m;
	} else {
		start = 0;
	}
	if ( start > PRIME_SMALL_PRIMORIAL && start > small_primes[ctx->w] ) {
		n = (start/PRIME_SMALL_PRIMORIAL) * PRIME_SMALL_PRIMORIAL;
		for ( i = PRIME_SMALL_PRIMORIAL_W+1 ; i <= ctx->w ; i++ ) {
			p = small_primes[i];
			ctx->wv[i] = p;
			m = ui_ceil_ratio(n,p)*(long)p;										// m is the least multiple of p >= n - could be > PRIME_MAX_PRIME
			m /= p;
			k = m % PRIME_SMALL_PRIMORIAL;
			while ( ! small_primorial_map[k] ) k++;								// find least k >= m%PRIME_SMALL_PRIMORIAL relatively prime to PRIME_SMALL_PRIMORIAL
			m = (long)p*((m/PRIME_SMALL_PRIMORIAL)*PRIME_SMALL_PRIMORIAL+k) - n;
			ctx->wv[i] = (long) m;
			ctx->wi[i] = small_primorial_map[k]-1;								// gap index i gives gap between (i+1)st and (i+2)nd integers relatively prime to PRIME_SMALL_PRIMORIAL
		}
		ctx->pi = PRIME_SMALL_PRIMES;		// just need to make it bigger than ctx->w and PRIME_SMALL_PRIMORIAL_W
		ctx->p = n-1;
	} else {
		ctx->pi = 1;
		ctx->p = 0;
	}
	ctx->j = ctx->wheel->phi-1;

	if ( start ) {
		while ( (p=prime_enum(ctx)) < start );
		if ( p != start ) { fprintf (stderr, "Unable to find start prime, p = %lu, start = %lu, L = %lu\n", p, start, L);  abort(); }
	}
	return ctx;
}

/*
	Once called, this function will keep a log that allows fast primality testing for all primes in [p-window,p]
	that have been enumerated ufsing prime_enum since logging was enabled (p is the most recently logged prime).
	
	*** IMPORTANT: logging can be active with only one enumeration, the log is statically allocated!!! ***

	TODO: change to allocate private bitmap

	Logging typically adds 5-10% to the cost of prime enumeration.

	Note, that the log only includes primes that have been explicitly enumerated after logging was started (not including start)
	Tests of primes that were not logged will fail, even if they are in [p-window,p]
*/

void prime_start_logging (long start, int window)
{
	if ( window < 4096 ) window = 4096;								// make sure window is comfortably bigger than the largest prime gap we expect to ever see (largest known gap is 1476 at around 60-bits)
	if ( (window&0xFFF) ) window += 4096 - (window&0xFFF);				// round up to a multiple of 4096
	window += 4096;												// add a buffer so we can clear efficiently 
	prime_logsize = window/2;										// we only log odd values, so window size is cut in half
	bitmap_alloc(prime_logsize);
	prime_logindex = 0;  prime_logvalue = start;
	if ( !(prime_logvalue&1) ) prime_logvalue++;					// make sure logvalue is odd, we won't log 2
	prime_logstart = start;
}
void prime_stop_logging () { prime_logvalue = 0; }

static inline void prime_log (long p)
{
	register unsigned long *ptr, m;
	register long newindex;
	
	newindex = prime_logindex + ((p-prime_logvalue)>>1);
	if ( (newindex&0x7FF) < (prime_logindex&0x7FF) ) {
		if ( newindex >= prime_logsize ) newindex -= prime_logsize;
		m = (newindex>>6) & 0x7FFFFFFFFFFFFFE0L;						// offset of unsigned long at the start of the 2048-bit block containing newindex
		ptr = bitmap+m;
		// clear the next 32*64=2048 bits corresponding to 4096 odd entries, enough to cover the next prime gap
		*ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;
		*ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;
		*ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;
		*ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;  *ptr++ = 0;
	}
	prime_logindex = newindex;
	prime_logvalue = p;
	bitmap_set(prime_logindex);									// if p is 2 then this logs 3 (which will get logged again when p=3).  That's ok.
}

/*
	Returns 1 (or 0) if n is in the log window and marked (or not).  Returns -1 if n is odd and not in the log window.  Returns 1 if n=2 and 0 if n>2 is even.
*/
int prime_check_log (long n)
{
	register long offset;
	
	if ( ! (n&1) ) return (n==2?1:0);									// handle even n here
	offset = (prime_logvalue - n)>>1;
	if ( n <= prime_logstart || n > prime_logvalue || offset >= prime_logsize ) return -1;
	offset = prime_logindex - ((prime_logvalue - n)>>1);
	if ( offset < 0 ) offset += prime_logsize;
	return bitmap_test(offset);
}

long prime_enum_powers (prime_enum_ctx_t *ctx)
{
	register int i;
	register long p;
	
	p = ( ctx->bufp ? ctx->bufp : _prime_enum (ctx) );
	if ( ! p || ctx->q[ctx->k] < p ) {
		i = ctx->k;
		ctx->bufp = p;
		p = ctx->q[i];
		ctx->d = ctx->k;
		ctx->b = small_primes[ctx->s[i]];
		ctx->s[i]++;
		ctx->q[i] = ipow (small_primes[ctx->s[i]], i);
		// find new min prime power
		for ( i = 2 ; i <= ctx->powers ; i++ ) if ( ctx->q[i] < ctx->q[ctx->k] ) ctx->k = i;
	} else {
		ctx->d = 1;
//		ctx->b = p;
		ctx->bufp = 0;
	}
	if ( p > ctx->end ) return 0;
	return p;
}


long _prime_enum (prime_enum_ctx_t *ctx)
{
	register long gap;
	register unsigned char *mp, *np;
	register long i, p;

	// The first ctx->w primes need to be enumerated seperately, since they have all been sieved out of the map
	// Note that if ctx->w > PRIME_SMALL_PRIMORIAL this means the first map is effectively empty, but that's ok.
	// Prime powers are only relevant for the first ctx->w primes since ctx->w was chosen to ensure pi(ctx->w)^2 > L
	if ( ctx->pi <= ctx->w || ctx->pi <= PRIME_SMALL_PRIMORIAL_W ) {
		p = small_primes[ctx->pi++];
		if ( p > ctx->end ) return 0;
		if ( prime_logvalue ) prime_log (p);
		return p;
	}
	for (;;) {
		if ( ctx->j == ctx->wheel->phi-1 ) {
			memset (ctx->map, 0, ctx->wheel->n);
			np = ctx->map + ctx->wheel->n;
			for ( i = PRIME_SMALL_PRIMORIAL_W+1 ; i <= ctx->w ; i++ ) {
				p = small_primes[i];
				mp = ctx->map+ctx->wv[i];
				while ( mp < np ) {
					*mp = 1;
					gap = p*ctx->wheel->gaps[ctx->wi[i]++];
					if ( ctx->wi[i] == ctx->wheel->phi ) ctx->wi[i] = 0;
					ctx->wv[i] += gap;
					mp += gap;
				}
				ctx->wv[i] -= ctx->wheel->n;
			}
			ctx->mp = ctx->map+1;
			if ( ctx->p ) {
				ctx->p += 2;									// last gap is always 2
				ctx->j = 0;
			} else {
				ctx->p = 1+ctx->wheel->gaps[0]; 					// special case - skip 1
				ctx->mp += ctx->wheel->gaps[0]; 
				ctx->j = 1;
			}
		} else {
			gap = ctx->wheel->gaps[ctx->j++];
			ctx->mp += gap;
			ctx->p += gap;
		}
		if ( ctx->p > ctx->end ) return 0;
		if ( ! *ctx->mp ) {
			if ( prime_logvalue ) prime_log (ctx->p);
			return ctx->p;
		}
	}
}

void prime_enum_end (prime_enum_ctx_t *ctx)
{
	if ( ctx->la_gaps ) free (ctx->la_gaps);
	free (ctx->map);
	free (ctx->wv);
	free (ctx->wi);
	_primorial_wheel_free (ctx->wheel);
	free (ctx);
}

/*
	Windowed version of fast prime enumeration.  Supports both trailing window (specify window < 0)
	and a centered window (window/2 on both sides);
	The external functionality is identical to the standard fast prime enum calls (without support for prime powers),
	but it will (greatly) accelerate calls to ui_is_prime for integers in the window.
*/
prime_enum_ctx_t *prime_enum_start_w (long start, long end, long window)
{
	prime_enum_ctx_t *ctx;
	register long size, gap, ep, np, mp;

	if ( window < 0 ) {
		if ( ! (ctx = prime_enum_start (start,end,0)) ) return 0;
		prime_start_logging(start,-window);
		return ctx;
	} else if ( window == 0 ) {
		return prime_enum_start(start,end,0);
	} else {
		if ( (window&1) ) window++;
		window += 4096;								// pad to cover prime gaps (this is surely overkill, the largest known gap is 1476, which covers all p < 2^60).
		if ( ! (ctx = prime_enum_start (start-window/2,end+window/2,0)) ) return 0;
	}
	ctx->la_end = end;
	size = window/8+256;
	if ( size&0xFF ) size += 256 - (size&0xFF);
	ctx->la_gaps = _calloc (size);						// allocate enough memory to be sure we have at least 1 byte per prime in any interval [N,N+window/2]
	ctx->la_n = size;
	ctx->la_i = ctx->la_j = 0;
	ctx->la_window = window/2;
	prime_start_logging(start-window/2,window);
	for ( np = _prime_enum(ctx) ; np && np < start ; np = _prime_enum(ctx) );
	assert ( np );
	ctx->la_p = ep = start;
	mp = ctx->la_p+ctx->la_window;
	while ( ep <=  mp ) {
		if ( ! np ) np = _prime_enum (ctx);
		if ( ! np ) { ep = 0; break; }
		gap = (np-ep);									// we could slice all but the gap between 2 and 3 in half, but there isn't much advantage to doing so
		ep = np;  np = 0;
		while ( gap &0xFF00 ) {							// special code to handle gaps >= 256, rarely used
			ctx->la_gaps[ctx->la_j++] = 255;				// note that a gap of 255 is impossible, so this value necessarily indicates a segmented gap
			gap -=255;
		}
		ctx->la_gaps[ctx->la_j++] = gap;
	}
	ctx->la_ep = ep;
	return ctx;
}

long _prime_enum_w (prime_enum_ctx_t *ctx)
{
	register long p, ep, mp, np, gap;

	p = ctx->la_p;
	do {											// compute next p using gap list, handling segmentation as needed
		if ( ctx->la_i == ctx->la_n ) ctx->la_i = 0;			// wrap as needed
		gap = ctx->la_gaps[ctx->la_i++];
		p += gap; 
	} while ( gap == 255 );
	if ( p > ctx->la_end ) return 0;
	ctx->la_p = p;
	ep = ctx->la_ep;
	if ( ! ep ) return p;
	mp = ctx->la_p+ctx->la_window;
	while ( ep <= mp ) {
		np = _prime_enum (ctx);
		if ( ! np ) { ep = 0; break; }
		gap = (np-ep);									// we could slice all but the gap between 2 and 3 in half, but there isn't much advantage to doing so
		ep = np;
		while ( gap &0xFF00 ) {							// special code to handle gaps >= 256, rarely used
			if ( ctx->la_j == ctx->la_n ) ctx->la_j = 0;		// wrap as needed
			ctx->la_gaps[ctx->la_j++] = 255;				// note that a gap of 255 is impossible, so this value necessarily indicates a segmented gap
			gap -=255;
		}
		if ( ctx->la_j == ctx->la_n ) ctx->la_j = 0;			// wrap as needed
		ctx->la_gaps[ctx->la_j++] = gap;
	}
	ctx->la_ep = ep;
	return p;
}

long next_prime (long p)
{
	mpz_t P;
	
	if ( p < 2 ) return 2;
	if ( p < PRIME_MAX_SMALL_PRIME ) {
		prime_setup();
		return nth_small_prime(small_prime_pi(p)+1);
	}
	mpz_init_set_ui (P, p);
	mpz_nextprime(P,P);
	if ( ! mpz_fits_slong_p (P) ) { mpz_clear (P); return 0; }
	p = mpz_get_si (P);
	mpz_clear (P);
	return p;
}

long prime_pi (long x)
{
	prime_enum_ctx_t *ctx;
	register long n;
	
	if ( x <= PRIME_MAX_SMALL_INTEGER ) return small_prime_pi ((int) x);
	ctx = prime_enum_start (0,x,0);
	for ( n = 0 ; _prime_enum(ctx) ; n++ );
	prime_enum_end (ctx);
	return n;
}


static inline long *_build_leaves (long *tp, long off, unsigned char *gaps, long numgaps)
{
	long *leaves;
	register unsigned char *gap;
	register long i, j, n;
	
	leaves = malloc ((numgaps/PRIME_TABLE_GAPS_PER_LEAF+2)*sizeof(*leaves));
	for ( n = off, gap = gaps, j = 0, i = PRIME_TABLE_GAPS_PER_LEAF ; *gap ; i++ ) {
		if ( i == PRIME_TABLE_GAPS_PER_LEAF ) { leaves[j++] = n; i = 0; }
		n += 2*(*gap++);
	}
	assert (gap-gaps == numgaps);
	*tp = n;
	leaves[j] = LONG_MAX;
	return leaves;
}

static inline long *_summarize_level (long *detail, long len)
{
	long *summary;
	register long i, j;
	
	summary = malloc ((len/PRIME_TABLE_PRIMES_PER_NODE+2)*sizeof(*summary));
	for ( i = j = 0 ; i < len ; i += PRIME_TABLE_PRIMES_PER_NODE ) { summary[j++] = detail[i]; }
	summary[j] = LONG_MAX;
	return summary;
}

void gap_table_init (prime_table_t tab, long off, unsigned char *gaps, long n)
{
	assert ( !gaps[n] );
	tab->si = 0;  tab->ti = n;
	tab->sp = off;
	tab->gaps = gaps;
	tab->levels[tab->top=0] = _build_leaves (&tab->tp, tab->sp, tab->gaps, n);
	n /= PRIME_TABLE_GAPS_PER_LEAF;
	while ( n ) {
		assert ( tab->top+1 < PRIME_TABLE_MAX_LEVELS );
		tab->levels[tab->top+1] = _summarize_level (tab->levels[tab->top], n);  n /= PRIME_TABLE_PRIMES_PER_NODE;
		tab->top++;
	}
	if ( tab->top ) {
		free (tab->levels[tab->top]);
		tab->levels[tab->top--] = 0;
	}
}

void prime_table_init (prime_table_t tab, long start_pi, long end_pi)
{
	prime_enum_ctx_t *ctx;
	unsigned char *gaps;
	register long n, p, firstp, lastp;

	assert (start_pi > 1 && end_pi > start_pi && end_pi <= PRIME_TABLE_MAX_PI );
	lastp = end_pi*ceil(log(end_pi)+log(log(end_pi))-0.5);			// upper bound on p_end given by Thm 5.21.ii on p. 91 of Shoup, valid for end_pi >= 20
	if ( lastp < 100 ) lastp = 100;
	gaps = malloc(end_pi-start_pi+1);
	ctx = prime_enum_start (0, lastp, 0);
	for ( n = 0 ; n < start_pi ; n++ ) p=prime_enum(ctx);
	firstp = lastp = p;
	for ( n = 0 ; n < end_pi-start_pi ; n++ ) {
		p = prime_enum(ctx);
		gaps[n] = (p-lastp)/2;
		lastp = p;
	}
	gaps[n] = 0;
	prime_enum_end (ctx);
	gap_table_init (tab, firstp, gaps, n);
	tab->si = start_pi;  tab->ti = end_pi;
}

void gap_table_clear (gap_table_t tab)
{
	tab->si = tab->ti = 0;
	tab->sp = tab->tp = 0;
	free (tab->gaps);  tab->gaps = 0;
	for ( int i = 0 ; i <= tab->top ; i++ ) { free (tab->levels[i]); tab->levels[i] = 0; }
	tab->top = 0;
}

