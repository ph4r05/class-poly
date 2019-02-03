#ifndef _NTUTIL_INCLUDE_
#define _NTUTIL_INCLUDE_

/*
	Copyright (c) 2007-2014 Andrew V. Sutherland
	See LICENSE file for license details.
*/

#include <assert.h>
#include <math.h>

/*
	Elementary number-theoretic functions implemented as inlines
*/

	
// Algorithm 1.4.10, p. 29 of [CANT] trimmed down to handle a >= 0, b > 0 odd (for the sake of space and speed)
static inline int ui_legendre (unsigned long a, unsigned long b)
{
	register unsigned long r;
	register int k, v;
	
	if ( a > b ) a %= b;
	k = 1;
	while ( a ) {
		for ( v = 0 ; ! (a&0x1) ; v++ ) a >>= 1;
		if ( v&0x1 )  if ( (b&0x7) ==3 || (b&0x7) ==5 ) k = -k;
		if ( (a&0x3) == 3 && (b&0x3) == 3 ) k = -k;
		r = a;   a = (b < (r+r) ? b-r : b%r);  b = r;
	}
	return ( b == 1 ? k : 0 );
}

// computes the legendre symbol for any integer a and odd prime b (in fact b can be any positive odd number)
static inline int legendre (long a, long b)
{
	register int s;
	
	s = 1;
	if ( a < 0 ) { s = ((b&3)==1 ? 1 : -1); a = -a; }
	return s*ui_legendre (a,b);
}


// Algorithm 1.4.10, p. 29 of [CANT] for a,b >= 0
static inline int ui_kronecker (unsigned long a, unsigned long b)
{
	register unsigned long r;
	register int k, v;
	
	if ( ! b ) return ( a==1 ? 1 : 0 );
	k = 1;
	if ( !(b&1) ) {
		if ( !(a&1) ) return 0;
		for ( v = 0 ; ! (b&0x1) ; v++ ) b >>= 1;
		if ( (v&1) && ( (a&0x7) == 3 || (a&0x7) == 5 ) ) k = -1;
	}
	if ( a > b ) a %= b;
	while ( a ) {
		for ( v = 0 ; ! (a&0x1) ; v++ ) a >>= 1;
		if ( v&0x1 )  if ( (b&0x7) ==3 || (b&0x7) ==5 ) k = -k;
		if ( (a&0x3) == 3 && (b&0x3) == 3 ) k = -k;
		r = a;   a = (b < (r+r) ? b-r : b%r);  b = r;
	}
	return ( b == 1 ? k : 0 );
}

// computes the kronecker symbol for arbitrary integers a and b
static inline int kronecker (long a, long b)
{
	int s;
	
	if ( b < 0 ) { s = (a<0?-1:1); b = -b; } else s = 1;
	if ( a < 0 ) { long c; if ( ! b ) return (a==-1?1:0); for ( c = b ; !(c&1) ; c >>= 1 ); s *= ((c&3)==3 ? -1 : 1); a = -a; }
	return s*ui_kronecker (a,b);
}

// uses standard Euclidean algorithm to invert a mod m.
static inline unsigned long ui_inverse (unsigned long a, unsigned long m)
{
	register unsigned long q, r, r0, r1;		// amazingly in this day and age, this register declaration gives a 10% improvement
	register long t, t0, t1;

	if ( a >= m ) a %= m;
	if ( a == 0 ) return 0;

	t1 = 1;  t0 = 0;  r0 = m;  r1 = a;
	while ( r1 > 0 ) {
		q = r0/r1;
		r = r0 - q*r1;
		r0 = r1;  r1 = r;
		t = t0 - q*t1;
		t0 = t1;  t1 = t;
	}
	if ( r0 > 1 ) return 0;
	if ( t0 < 0 ) return m - ((unsigned long)(-t0));
	return (unsigned long)t0;
}

// simple Euclidean gcd, not extended
static inline unsigned long ui_gcd (unsigned long a, unsigned long b)
{
	register unsigned long q, r, r0, r1;
	
	if ( a < b ) { r0 = b;  r1 = a; } else { r0 = a; r1 = b; }
	while ( r1 > 0 ) {
		q = r0/r1;  r = r0 - q*r1;
		r0 = r1;  r1 = r;
	}
	return r0;
}
static inline unsigned long ui_lcm (unsigned long a, unsigned long b) { return (a/ui_gcd(a,b))*b; }

static inline unsigned long ui_binomial (int n, int k)
{
	unsigned long a, b;
	int i;
	
	assert ( n <= 28 );	// to avoid 64-bit overflow, require n <= 28
	if ( k > n/2 ) k = n-k;
	if ( ! k ) return 1;
	for ( a = n, b=1, i = 1 ; i < k ; i++ ) { a *=(n-i);  b *= (i+1); }
	return a/b;
}


/*
static int _mod64res[64] = { 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, };
static int _mod63res[63] = { 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 };
static int _mod65res[65] = { 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, };
static int _mod11res[11] = { 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0 };
*/

#define _mod64res	0x202021202030213
#define _mod63res	0x402483012450293
#define _mod65res	0x218a019866014613
#define _mod11res	0x23b


// return sqrt(n) if n is a perfect nonzero square, 0 ow.
static inline long i_sqrt(long n)
{
	register long x;
	register int i;
	
	if ( n <= 0 ) return 0;
	if ( ! (_mod64res&(1UL<<(n&0x3F))) ) return 0;
	if ( ! (_mod63res&(1UL<<(n%63))) ) return 0;
	if ( (i=(n%65)) < 64 && ! (_mod65res&(1UL<<i)) ) return 0;
	if ( ! (_mod11res&(1UL<<(n%11))) ) return 0;
	x = (long)(sqrt(n)+0.1);		// avoid potential rounding problems
	return (x*x == n ? x : 0);
}

#endif
