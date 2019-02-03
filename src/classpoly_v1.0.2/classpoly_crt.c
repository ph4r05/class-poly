#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "crt.h"
#include "classpoly_crt.h"
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

void classpoly_crt_start (classpoly_crt_t crt, long D, int inv, unsigned long *primes, int pcnt, int ccnt, mpz_t P)
{
	clock_t start, end;
	char prefix[32];
	int fcnt;
	
	memset (crt, 0, sizeof(crt[0]));
	crt->D = D;
	crt->inv = inv;
	crt->p = primes;
	crt->pcnt = pcnt;
	crt->ccnt = ccnt;
	if ( P ) mpz_init_set (crt->P, P); else mpz_init (crt->P);
	if ( mpz_sgn(crt->P) ) {
		start = clock();
		sprintf (prefix, "H_%ld", -D);
		ecrt_init (crt->ecrt, primes, pcnt, ccnt, crt->P, 0, 0, prefix);
		end = clock();
		info_printf ("Explicit CRT precomputation completed in %ld msecs\n", delta_msecs(start,end));
	} else {
		if ( ! (fcnt=crt_create_files (&crt->fp, -D, 0, ccnt)) ) { err_printf("Call to crt_create_files failed\n"); abort(); }
		info_printf ("Created %d CRT output files\n", fcnt);				
	}
	crt->index = -1;
}

void classpoly_crt_end (classpoly_crt_t crt)
{
	if ( crt->fp ) crt_close_files (crt->fp, crt->ccnt);
	mpz_clear (crt->P);
	if ( crt->jcache ) mem_free (crt->jcache);
	mem_free (crt->p);
}

int classpoly_crt_output_coefficient (mpz_t c, int i, void *ctx)
	{ gmp_fprintf ((FILE*)ctx, "%Zd*X^%d + \n", c, i);  return 1; }

int classpoly_crt_finish (classpoly_crt_t crt, long *totbits, long *maxcbits, FILE *fp)
{
	mpz_t X;
	clock_t start, end;
	register int i;

	assert (crt->index == crt->pcnt);
	if ( totbits ) *totbits = 0;
	if ( maxcbits ) *maxcbits = 0;
	start = clock();
	if ( crt->fp ) {
		crt_close_files(crt->fp, crt->ccnt);
		info_printf ("Wrote %ld bytes of data to CRT files\n", crt->ccnt*crt->pcnt*sizeof(unsigned long));
	}
	gmp_fprintf (fp, "I=%d\n", crt->inv);
	gmp_fprintf (fp, "D=%ld\n", crt->D);
	if ( mpz_sgn(crt->P) ) gmp_fprintf (fp, "P=%Zd\n", crt->P);

	if ( mpz_sgn(crt->P) ) {		
		// all ecrt cases treated as "small" in the sense of section 6 of "Computing Hilbert class polynomials with the Chinese Remainder Theorem")
		ecrt_finalize (crt->ecrt);
		mpz_init (X);
		for ( i = 0 ; ecrt_next_coeff(X, crt->ecrt) ; i++ ) gmp_fprintf (fp, "%Zd*X^%d + \n", X, i);
		mpz_clear (X);
		if ( i != crt->ccnt ) { err_printf ("ECRT postcomputation returned only %d of %d expected coefficients\n", i, crt->ccnt); mpz_clear (crt->P); return 0; }
		if ( totbits ) *totbits = crt->ccnt * mpz_sizeinbase(crt->P,2);
		if ( maxcbits ) *maxcbits = mpz_sizeinbase(crt->P,2);
		mpz_clear (crt->P);
	} else {
		mpz_clear (crt->P);
		if ( ! crt_process_files (totbits, maxcbits, -crt->D, 0, crt->ccnt, crt->p, crt->pcnt, classpoly_crt_output_coefficient, fp) ) { err_printf ("Error processing CRT output files\n");  return 0; }
	}
	gmp_fprintf (fp, "1*X^%d\n", crt->ccnt);
	end = clock();
	if ( crt->jcache ) mem_free (crt->jcache);
	mem_free (crt->p);
	info_printf ("Completed CRT postcomputation in %ld msecs\n", delta_msecs(start,end));
	return 1;
}
