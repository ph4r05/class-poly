#ifndef _POLYCOSTS_INCLUDE_
#define _POLYCOSTS_INCLUDE_

#include <limits.h>
#include <stdio.h>
#include "mpzutil.h"
#include "class_inv.h"
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

#define POLY_COST_PRIMES				31
#define POLY_COST_MAX_PRIME			127
#define POLY_COST_ENUM_ENTRIES			496
#define POLY_COST_MAX_DEGREE			511
#define POLY_COST_MAX_LOG_ROOTS		25		// for build costs, degree=roots can be up to 2^26

extern int poly_root_enum_costs_32[POLY_COST_PRIMES];
extern int poly_root_enum_costs_40[POLY_COST_PRIMES];
extern int poly_root_enum_costs_48[POLY_COST_PRIMES];
extern int poly_gcd_cycle_costs[(POLY_COST_PRIMES-1)*(POLY_COST_PRIMES-1)];
extern int poly_gcd_cycle_costs_s2[(POLY_COST_PRIMES-1)*(POLY_COST_PRIMES-1)];
extern int poly_gcd_cycle_costs_s3[(POLY_COST_PRIMES-1)*(POLY_COST_PRIMES-1)];
extern int poly_gcd_cycle_costs_s4[(POLY_COST_PRIMES-1)*(POLY_COST_PRIMES-1)];
extern int poly_gcd_cycle_costs_s6[(POLY_COST_PRIMES-1)*(POLY_COST_PRIMES-1)];
extern int poly_gcd_cycle_costs_s8[(POLY_COST_PRIMES-1)*(POLY_COST_PRIMES-1)];
extern int poly_gcd_cycle_costs_s12[(POLY_COST_PRIMES-1)*(POLY_COST_PRIMES-1)];
extern int poly_gcd_cycle_costs_s24[(POLY_COST_PRIMES-1)*(POLY_COST_PRIMES-1)];
extern int poly_one_root_costs_48[POLY_COST_MAX_DEGREE+1];
extern int poly_split_root_costs_48[POLY_COST_MAX_DEGREE+1];
extern long poly_build_costs_32[10*POLY_COST_MAX_LOG_ROOTS+1];
extern long poly_build_costs_40[10*POLY_COST_MAX_LOG_ROOTS+1];
extern long poly_build_costs_48[10*POLY_COST_MAX_LOG_ROOTS+1];

static inline long poly_root_step_cost (int ell, int d, long p, int inv)
{
	long x;
	int i, n;

	// currently we ignore the invariant here -- sparseness helps instantiation, but this cost is negligible compared to the root-finding cost
	n = ui_len(p);
	if ( d > 1 || ( d > 0 && ell > 2 ) ) {
		if ( ell > POLY_COST_MAX_DEGREE ) return LONG_MAX;
		x = (n*poly_split_root_costs_48[ell])/48 + 16*ell*ell;							// a bit pessimistic
		return x * (1 + (ell-1)*d/2);											// pessimistic heuristic, based on Prop. 4 of HCP paper
	} else {
		if ( ell > POLY_COST_MAX_PRIME ) return  (11*n)/4 * (long)ell * (long)ell;		// pessimistic heuristic
		i = ui_small_prime_index(ell)-1;
		if ( n < 36 ) return (long)(n*poly_root_enum_costs_32[i])/32;
		if ( n < 44 ) return (long)(n*poly_root_enum_costs_40[i])/40;
		return (long)(n*poly_root_enum_costs_48[i])/48;
	}
}

static inline long poly_root_enum_cost (int ell, int d, long n, long p, int inv) { return (n-1) * poly_root_step_cost(ell,d,p,inv); }

// p2 > p1, we assume we are using gcd_cycle (performs 2 inversions in parallel), which is a bit optimistic...
static inline long poly_gcd_step_cost (int ell1, int ell2, long p, int inv)
{
	int i,j,d;
	
	d=inv_level(inv);
	if ( d > 1 && ( !(ell1%d) || !(ell2%d) ) ) { printf("ell1=%d or ell2=%d divides level %d in poly_gcd_step_cost!\n", ell1, ell2, d); exit (0); }
	if ( ell1 > POLY_COST_MAX_PRIME ) return 5 * (long)ell1 * (long)ell1 + 2 * (long) ell2 * (long)(ell2);	// heuristic for largish ell1, ell2
	i = ui_small_prime_index(ell1)-1;
	if ( ell2 > POLY_COST_MAX_PRIME ) {
		switch(inv_sparse_factor(inv)){
		case 2: return poly_gcd_cycle_costs_s2[i*(POLY_COST_PRIMES-1)] +  2 * (long) (ell2+ell1) * (long)(ell2-ell1);
		case 3: return poly_gcd_cycle_costs_s3[i*(POLY_COST_PRIMES-1)] +  2 * (long) (ell2+ell1) * (long)(ell2-ell1);
		case 6: return poly_gcd_cycle_costs_s6[i*(POLY_COST_PRIMES-1)] +  2 * (long) (ell2+ell1) * (long)(ell2-ell1);
		case 8: return poly_gcd_cycle_costs_s8[i*(POLY_COST_PRIMES-1)] +  2 * (long) (ell2+ell1) * (long)(ell2-ell1);
		case 12: return poly_gcd_cycle_costs_s12[i*(POLY_COST_PRIMES-1)] +  2 * (long) (ell2+ell1) * (long)(ell2-ell1);
		case 24: return poly_gcd_cycle_costs_s24[i*(POLY_COST_PRIMES-1)] +  2 * (long) (ell2+ell1) * (long)(ell2-ell1);
		default:	return poly_gcd_cycle_costs[i*(POLY_COST_PRIMES-1)] +  2 * (long) (ell2+ell1) * (long)(ell2-ell1);
		}
	}
	j = ui_small_prime_index(ell2)-ui_small_prime_index(ell1)-1;
	switch(inv_sparse_factor(inv)){
	case 2:  return poly_gcd_cycle_costs_s2[i*(POLY_COST_PRIMES-1)+j];
	case 3:  return poly_gcd_cycle_costs_s3[i*(POLY_COST_PRIMES-1)+j];
	case 6:  return poly_gcd_cycle_costs_s6[i*(POLY_COST_PRIMES-1)+j];
	case 8:  return poly_gcd_cycle_costs_s8[i*(POLY_COST_PRIMES-1)+j];
	case 12:  return poly_gcd_cycle_costs_s12[i*(POLY_COST_PRIMES-1)+j];
	case 24:  return poly_gcd_cycle_costs_s24[i*(POLY_COST_PRIMES-1)+j];
	default: return poly_gcd_cycle_costs[i*(POLY_COST_PRIMES-1)+j];
	}
}

static inline long poly_gcd_enum_cost (int ell1, int ell2, int e, int n, int d, long p, int inv)
{
//printf ("cost(%d,%d,%d,%d,%d):  %d*%ld + %d*%ld = %ld\n", ell1, ell2, e, n, d, (e-1), poly_root_step_cost(ell1, d, p, inv), (n-e), poly_gcd_step_cost(ell1,ell2,p,inv), e*poly_root_step_cost(ell1, d, p, inv) + (n-e)*poly_gcd_step_cost(ell1,ell2,p,inv));
	return (e-1)*poly_root_step_cost(ell1, d, p, inv) + (n-e)*poly_gcd_step_cost(ell1,ell2,p,inv);
}

static inline long poly_build_cost (int deg, long p)
{
	long n;
	int i;
	
	n = ui_len(p);
	i = ceil(log2(deg)*10);
	if ( i > 10*POLY_COST_MAX_LOG_ROOTS ) {
		// if i is too big, scale linearly
		if ( n < 36 ) return (n * poly_build_costs_32[10*POLY_COST_MAX_LOG_ROOTS] * i) / (32*(10*POLY_COST_MAX_LOG_ROOTS));
		if ( n < 44 ) return (n * poly_build_costs_40[10*POLY_COST_MAX_LOG_ROOTS] * i) / (40*(10*POLY_COST_MAX_LOG_ROOTS));
		return (n * poly_build_costs_48[10*POLY_COST_MAX_LOG_ROOTS] * i) / (48*(10*POLY_COST_MAX_LOG_ROOTS));
	} else {
		if ( n < 36 ) return (n*poly_build_costs_32[i])/32;
		if ( n < 44 ) return (n*poly_build_costs_40[i])/40;
		return (n*poly_build_costs_48[i])/48;
	}
}

#endif
