#ifndef _PICKPRIMES_INCLUDE_
#define _PICKPRIMES_INCLUDE_

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

#define MIN_CRT_PRIME		10//20000		// avoid primes that are so small that we may find no curves satisfying torsion constraints

int pick_primes_cornacchia (long *plist, int *tlist, int *vlist, int maxn, long D, long h, double bits);
int pick_primes (long **pplist, int **ptlist, int **pvlist, long D, long h, double bits, int pfilter, long ellfilter);
int pick_primes_old (long *plist, int *tlist, int *vlist, int n, long D, long h, double bits, int filter);
double split_prime_rating (long p, long t, int *ptwist, int *ptor, int *ps2, int *pt3);
double split_prime_rating_new (long p, long t, int *ptwist, int *ptor, int *ps2, int *pt3);

#endif
