#ifndef _FF2K_INCLUDE_
#define _FF2K_INCLUDE_

/*
    Copyright 2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

// very rudimentary implementation of arithmetic in F_{2^k}
// this is not meant to be fast, just functional

#define FF2K_MAXK			31		// it would not be hard to increase this to 63

typedef unsigned long ff2k_t;

extern int _ff2k_k;
extern ff2k_t _ff2k_xk, _ff2k_xkr, _ff2k_x2km2, _ff2k_x2km2r, _ff2km, _ff2k_m;

void ff2k_setup (int k);

#define _ff2k_set(o,a)		((o) = (a))
#define _ff2k_equal(a,b)	((a) == (b))
#define _ff2k_add(o,a,b)	((o) = (a)^(b))
#define _ff2k_addto(o,a)	((o) ^= (a))
#define _ff2k_sub(o,a,b)	_ff2k_add(o,a,b)
#define _ff2k_subfrom(o,a)	_ff2k_addto(o,a)
#define _ff2k_inc(o)		((o) ^= 1UL)
#define _ff2k_zero(a)		(!(a))
#define _ff2k_set_zero(o)	((o)=0)
#define _ff2k_one(a)		((a) == 1UL)
#define _ff2k_set_one(o)	((o)=1UL)
#define _ff2k_negate(o)
#define _ff2k_next(o)		(++(o) == _ff2k_xk ? (o)=0 : (o))

static inline ff2k_t ff2k_multiply (ff2k_t a, ff2k_t b)
{
	register ff2k_t m, r, o;

	for ( o = 0 ; b ; b >>= 1, a <<= 1 ) if ( b&1 ) o ^= a;
	for ( m = _ff2k_x2km2, r = _ff2k_x2km2r ; m >= _ff2k_xk ; m>>=1, r>>=1 ) if ( o&m ) o ^= r;
	o&=_ff2k_m;
	return o;
}

#define _ff2k_mult(o,a,b)		((o)=ff2k_multiply(a,b))
#define _ff2k_square(o,a)		((o)=ff2k_multiply(a,a))	// currently squaring is not optimized

static inline void ff2k_primitive_root (ff2k_t o[1]) { if ( _ff2k_k == 1 ) o[0]=1; else o[0]=2; }

static inline int ff2k_poly_set_mpz (ff2k_t f[], mpz_t F[], int d)
	{ register int i; for ( i = 0 ; i <= d ; i++ ) if ( mpz_tstbit(F[i],0) ) _ff2k_set_one(f[i]); else _ff2k_set_zero(f[i]);  for ( i-- ; i >= 0 && _ff2k_zero(f[i]) ; i-- ); return i; }

// computes o = f(x)
static inline void ff2k_poly_eval (ff2k_t o[1], ff2k_t f[], int d, ff2k_t x[1])
{
	register ff2k_t t,y;
	register int i;
	
	if ( d < 0 ) { _ff2k_set_zero(o[0]); return; }
	_ff2k_set (y,f[d]);
	_ff2k_set (t,*x);
	for ( i = d-1 ; i >= 0 ; i-- ) { _ff2k_mult(y, y, t); _ff2k_addto (y,f[i]); }
	_ff2k_set(o[0],y);
}

// brute-force root-finding
int ff2k_poly_distinct_roots (ff2k_t r[], ff2k_t f[], int d);

// naive point-counting for elliptic curves in WS form over F_2^k, returns -1 if curve is singular
long ff2k_WS_pointcount (ff2k_t w[]);

// naive point-counting for hyperelliptic curves over F_2^k, assumes curve is not singular but does not check
long ff2k_hyperelliptic_pointcount (ff2k_t f[], int df, ff2k_t h[], int dh);

// returns # of points at infinity (assuming good reduction), 0, 1, or 2.
static inline int ff2k_hyperelliptic_pointcount_at_infinity (int df, int dh, int k)
	{ int g = ((df>2*dh?df:2*dh)-1)/2;  if ( dh <= g) return 1;  if ( df == 2*g+2 && (k&1) ) return 0;  return 2; }

static inline void ff2k_print (ff2k_t x)
	{ register int i; if ( _ff2k_zero(x) ) { printf ("(0)"); return; }printf ("("); for ( i = _ff2k_k-1 ; i >= 0 ; i-- ) { if ( x&(1<<i) ) switch (i) { case 0: printf("+1"); break; case 1: printf ("+z"); break; default: printf ("+z^%d", i); } } printf (")"); }
		
void ff2k_poly_print (ff2k_t f[], int d);

#ifdef __cplusplus
}
#endif

#endif
