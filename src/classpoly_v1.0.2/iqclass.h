#ifndef _IQCLASS_INCLUDE_
#define _IQCLASS_INCLUDE_

#include <gmp.h>
#include "mpzutil.h"
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

#define IQ_MAX_GENS			16								// maximum number of generators for class group (this should be plenty for small h)
#define IQ_MAX_RLEN			((IQ_MAX_GENS*(IQ_MAX_GENS-1))/2)	// maximum length of power relations in polycyclic presentation for the class group

#define IQ_FILTER_1MOD3		0x1L							// forces split primes to not be congruent to 1 mod 3
#define IQ_FILTER_1MOD4		0x2L							// forces split primes to not be congruent to 1 mod 3
#define IQ_FILTER_1MOD8		0x4L							// forces split primes to not be congruent to 1 mod 8
#define IQ_REQUIRE_1MODM		0x8L							// forces split primes to be congruent to 1 mod m, where m is specified in the upper 32 bits

// computes kronecker symbol (a/p) for prime p (including 2)
static inline int kronecker_p (long a, unsigned long p)
{
	register int k;
	
	if ( p == 2 ) return (a&1) * ( (((a&7)-3) * ((a&7)-5)) ?1:-1 );
	if ( a < 0 ) { k = ((p&3)==3?-1:1); a = -a; } else { k = 1; }
	return k*ui_legendre (a, p);
}

// Finds the least t > u s.t. 4p=t^2+dv^2 with p prime (note D=-d), subject to filter.  Set *p to the prime and *t to t
int next_split_prime (long *p, long *t, long d, long u, long v, long filter);
int mpz_next_split_prime (mpz_t P, mpz_t T, long d, long v);

// Given a fundamental discriminant D < -4 and v, computes h(v^2D)/h(D) (does not verify that D is fundamental)
long relative_h (long D, long v);

// Given a fundamental discriminant D < -4 and v, computes the sum of h(u^2D)/h(D) over all the divisors u of v
long sum_relative_h (long D, long v);

// Factors the conductor of the imaginary quadratic discriminant D.
// Returns the number of prime-power factors, 0 if D is fundamental, and -1 if D is not a valid discriminant
static inline int discriminant_factor_conductor (unsigned long p[], unsigned long h[], long D)
{
	register long i, j, k, D0;
	
	if ( D >= 0 ) return -1;
	i = (-D)&3;			// don't assume 2's complement is used for negative integers
	if ( i != 0 && i != 3 ) return -1;
	k = ui_factor(p,h,-D);
	D0 = -1;
	for ( i = 0 ; i < k ; i++ ) {
		if ( h[i]&1 ) D0 *= p[i];
		h[i] >>= 1;
	}
	if ( (D0&3) != 1 ) { 
		if ( p[0] != 2 || h[0] < 1 ) return -1;
		h[0]--;
	}
	for ( i = j = 0 ; i < k ; i++ ) if ( h[i] ) { p[j] = p[i]; h[j] = h[i]; j++; }
	return j;
}


// Tests whether D is the discriminant of an imaginary quadratic number field (D is negative!)
// Returns 0 if not, and otherwise the value v such that D/v^2 is a fundamental discriminant
#define discriminant_test	discriminant_conductor
#define discriminant_is_fundamental(D)	(discriminant_conductor(D)==1)

static inline long discriminant_conductor (long D)
{
	unsigned long p[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS];
	register long i, j, k, v;
	
	k = discriminant_factor_conductor (p, h, D);
	if ( k < 0 ) return 0;
	if ( k == 0 ) return 1;
	for ( v = 1, i = 0 ; i < k ; i++ ) for ( j = 0 ; j < h[i] ; j++ ) v *= p[i];
	return v;
}

/*
	Returns true if the prime p does not divide the conductor of D (assumes p is prime).
*/
static inline int prime_to_conductor(long D, long p)
{
	register long b;
	
	if ( p >2 ) return (D%(p*p));
	b = D&0xF;					// if 2 divides the conductor of D then D=0 mod 16 (when D_0 is 0 mod 4) or D=4 mod 16 (when D_0 is 1 mod 4)
	return ( b && b!=4 );
}

#endif
