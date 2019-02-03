#ifndef _CLASSPOLY_CRT_INCLUDE_
#define _CLASSPOLY_CRT_INCLUDE_

#include <gmp.h>
#include "ff_poly.h"
#include "crt.h"

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

struct classpoly_crt_struct {
	unsigned long *p;
	ff_t *jcache;
	int pcnt;
	int jcnt;
	int ccnt;
	long D;
	int inv;
	mpz_t P;
	crt_file_t *fp;
	ecrt_context_t ecrt;
	int index;
};
typedef struct classpoly_crt_struct classpoly_crt_t[1];

void classpoly_crt_start (classpoly_crt_t crt, long D, int inv, unsigned long *primes, int pcnt, int ccnt, mpz_t P);
static inline void classpoly_crt_restart (classpoly_crt_t crt) { crt->index = -1; }
void classpoly_crt_end (classpoly_crt_t crt);	// use end to terminate early and cleanup, finish to terminate normally
int classpoly_crt_finish (classpoly_crt_t crt, long *totbits, long *maxcbits, FILE *fp);

static inline void classpoly_crt_cache_j (classpoly_crt_t crt, ff_t j)
	{ if ( ! crt->jcache ) crt->jcache = mem_alloc (crt->pcnt*sizeof(ff_t)); crt->jcache[crt->index] = j; crt->jcnt++; }

static inline int classpoly_crt_cached_j (ff_t *j, classpoly_crt_t crt)
	{ if ( crt->index >= crt->jcnt ) return 0;  *j = crt->jcache[crt->index];  return 1; }
		
static inline int classpoly_crt_next_prime (classpoly_crt_t crt)
	{ if ( ++crt->index >= crt->pcnt ) return 0; ff_setup_ui(crt->p[crt->index]); return 1; }
	
static inline void classpoly_crt_skip_prime (classpoly_crt_t crt)
{
	register int i;
	
	crt->pcnt--;
	if ( crt->jcnt >0 ) crt->jcnt--;
	for ( i = crt->index ; i < crt->pcnt ; i++ ) {
		crt->p[i] = crt->p[i+1];
		if ( i < crt->jcnt ) crt->jcache[i] = crt->jcache[i+1];
	}
	crt->index--;
}	

static inline void classpoly_crt_update (classpoly_crt_t crt, unsigned long *coeffs)
{
	if ( mpz_sgn(crt->P) ) {
		ecrt_update (crt->ecrt, crt->index, coeffs, crt->ccnt);
	} else {
		crt_write_files (crt->fp, crt->p[crt->index], coeffs, crt->ccnt);
	}
}

#endif
