#ifndef _TECURVE_INCLUDE_
#define _TECURVE_INCLUDE_

#include "ff_poly.h"

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

#define TECURVE_MAX_N			36
#define TECURVE_MAX_CURVES	1024

#define FILTER_2SYLOW_NONE		0			// no constraints on 2-Sylow
#define FILTER_2SYLOW_1		1			// 2-Sylow has order 1
#define FILTER_2SYLOW_NOT_1	2			// 2-Sylow has order > 1
#define FILTER_2SYLOW_2		3			// 2-Sylow has order 2
#define FILTER_2SYLOW_D4		4			// 2-Sylow has order at least 4
#define FILTER_2SYLOW_R1		5			// 2-Sylow has rank 1
#define FILTER_2SYLOW_R14		6			// 2-Sylow has rank 1 and order at least 4
#define FILTER_2SYLOW_R2		7			// 2-Sylow has rank 2
#define FILTER_2SYLOW_1MOD4	8			// group order is 1 mod 4 (2-Sylow is trvial)
#define FILTER_2SYLOW_3MOD4	9			// group order s 3 mod 4 (2-Sylow is trivial)
#define FILTER_2SYLOW_MAX		9

#define FILTER_2SYLOW_C_LG		0x80		// 2-Sylow is cyclic of order 2^k, where k is the low 7 bits.
#define FILTER_2SYLOW_C_LG_MIN	0x100		// used in conjuction with FILTER_2SYLOW_C_LG.  If set, only a lower bound is enforced.

#define FILTER_3TOR_NONE		0			// no constraints on 3-tor
#define FILTER_3TOR_1			1			// 3-torsion subgroup is trivial
#define FILTER_3TOR_NOT_1		2			// 3-torsion subgroup is non-trivial
#define FILTER_3TOR_3			3			// 3-torsion subgroup is Z/3Z
#define FILTER_3TOR_9			4			// 3-torsion subgroup is Z/3Z x Z/3Z
#define FILTER_3TOR_MAX		4

#define FILTER_JCUBE			1			// if set, ensures j-invariant is a cube in F_p
#define FILTER_WEBER			2			// if set, j-1728 must be a nonresidue in F_p

void tecurve_random_curves (ff_t f1[], ff_t x[], ff_t y[], int n, int tor, int s2_flag, int t3_flag);
int tecurve_random_curves_x (ff_t f1[], ff_t x[], ff_t y[], int n, int tor, int s2_flag, int t3_flag, int j_flags);

#define TECURVE_FREE_K		6
#define TECURVE_FIX_K			6
#define TECURVE_MAX_COST		1000000.0						// returned for unsupported settings

#include "tecurvecosts.h"

double tecurve_free2_density (long p, int N, int k, int t3, int jflags);
double tecurve_fix2_density (long p, int N, int k, int t3, int jflags);

// For fix2 costs, curve order must be divisible by 2^k and not by 2^{k+1} (and in fact is forced to have cyclic 2-Sylow)
static inline double tecurve_fix2_cost (long p, int N, int k, int t3, int jflags)
{
	int i;
	
	N--;
	i = (N<<2) + ((p%3)==1?2:0) + ((p&2)?0:1);
	if ( k >= TECURVE_FIX_K ) return TECURVE_MAX_COST;
	return tecurve_fix2_costs[i][k+(t3?TECURVE_FIX_K:0)];
}

// For free2 costs, order must be divisible by 2^k (for k=0 this means no constraint on 2-Sylow)
static inline double tecurve_free2_cost (long p, int N, int k, int t3, int jflags)
{
	int i;
	
	N--;
	if ( k >= TECURVE_FREE_K ) k=TECURVE_FREE_K-1;
	i = (N<<2) + ((p%3)==1?2:0) + ((p&2)?0:1);
	return tecurve_free2_costs[i][k+(t3?TECURVE_FREE_K:0)];
}

#endif
