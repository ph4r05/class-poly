#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "mpzutil.h"
#include "iqclass.h"
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


// Finds the least t > u s.t. 4p=t^2+dv^2 with p an odd prime (note D=-d).  Set *p to the prime and *t to t
int next_split_prime (long *p, long *t, long d, long u, long v, long filter)
{
	register long m,n;

	if ( (d&7)==7 && (v&1) ) { *p=0; *t = u; return 0; }
	if ( filter&IQ_FILTER_1MOD3 ) {
		if ( !(d%3) ) { printf ("D cannot be divisible by 3 when filtering 1mod3 is set!\n"); exit (0); }
		if ( !(v%3) ) { *p = 0; *t = u; return 0; }					// v=0 mod 3 forces p=1mod 3
	}
	if ( filter&IQ_FILTER_1MOD8 && !(v&7) ) { *p = 0; *t = u; return 0; }	// v=0 mod 8 forces p=1 mod 8
	if ( filter&IQ_FILTER_1MOD4 && !(v&3) ) { *p = 0; *t = u; return 0; }	// v=0 mod 4 forces p=1 mod 4
	if ( ui_len(d)+ui_len(v*v) > 62 ) { *p = 0; *t = u; return 0; }			// avoid overflow
	m = 0;													// shut up compiler
	if ( (filter&IQ_REQUIRE_1MODM) ) m = filter>>32;
	for ( u++ ; ; u++ ) {
		n = u*u+d*v*v;
		if ( n < 0 ) { *p = 0; *t=u; return 0; }						// also watch for overflow here
		if ( n&3 ) continue;
		n>>=2;
		if ( (filter&IQ_FILTER_1MOD3) && (n%3) != 2 ) continue;		// n must be 2 mod 3 if filtering 1 mod 3
		if ( (filter&IQ_FILTER_1MOD4) && (n&3) != 3 ) continue;		// n must be 3 mod 4 if filtering 1 mod 4
		if ( (filter&IQ_FILTER_1MOD8) && (n&7) == 1 ) continue;		// require n not be 1 mod 8 
		if ( (filter&IQ_REQUIRE_1MODM) && (n%m) != 1 ) continue;		// require n to be 1 mod m
		if ( ui_is_prime(n) ) break;
	}
	*p = n;
	*t = u;
	return 1;
}



// Finds the least t > T s.t. 4p=t^2+dv^2 with p prime (note D=-d).  Set P to the prime and T to t
int mpz_next_split_prime (mpz_t P, mpz_t T, long d, long v)
{
	if ( (d&7)==7 && (v&1) ) return 0;
	for ( mpz_add_ui (T,T,1) ; ; mpz_add_ui(T,T,1) ) {
		mpz_mul(P,T,T); mpz_add_ui(P,P,d*v*v);
		if ( mpz_tstbit(P,0) || mpz_tstbit(P,1) ) continue;
		mpz_div_2exp (P,P,2);
		if ( mpz_probab_prime_p(P, 5) ) break;
	}
	return 1;
}



// Given a fundamental discriminant D < -4 and v, computes h(v^2D)/h(D)
long relative_h (long D, long v)
{
	unsigned long p[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS];
	register long m, n;
	register int i, w;
	
	if ( D >= -4 || v <= 0 ) { printf ("D must be < -4 and v > 0 in relative_h (D=%ld, v=%ld)\n", D, v); exit (0); }
	w = ui_factor (p,h,v);
	m = n = 1;
	for ( i = 0 ; i < w ; i++ ) {
		m *= p[i] - kronecker_p(D,p[i]);
		n *= p[i];
	}
	return (m*v)/n;
}


// Given a fundamental discriminant D < -4 and v, computes the sum of h(u^2D)/h(D) over all the divisors u of v
long sum_relative_h (long D, long v)
{
	unsigned long p[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS],  m[MAX_UI_PP_FACTORS], e[MAX_UI_PP_FACTORS+1];
	register long sum, x;
	register int i, j, w;

	if ( D >= -4 || v <= 0 ) { printf ("D must be < -4 and v > 0 in sum_relative_h (D=%ld, v=%ld)\n", D, v); exit (0); }
	
	w = ui_factor (p, h, v);
	for ( i = 0 ; i < w ; i++ ) { m[i] = (p[i] - kronecker_p(D,p[i])); e[i] = 0; }
	sum = 0;
	do {
		x = 1;
		for ( i = 0 ; i < w ; i++ ) {
			if ( e[i] ) {
				x *= m[i];
				for ( j = 1 ; j < e[i] ; j++ ) x *= p[i];
			}
		}
		sum += x;
		e[0]++;
		for ( i = 0 ; i < w && e[i] > h[i] ; i++ ) { e[i] = 0; e[i+1]++; }
	} while ( i < w ); 
	return sum;
}
