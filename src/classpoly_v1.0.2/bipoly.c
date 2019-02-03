#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <memory.h>
#include <string.h>
#include <gmp.h>
#include "ff_poly.h"
#include "bipoly.h"
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


static inline int max (int a, int b) { return ( a > b ? a : b ); }

int bipoly_load (bipoly_t phi, char *filename, int sort_v)
{
	register int i, k;

	phi->num_terms = bipoly_load_mpz (&phi->mpz_terms, filename);
	if ( phi->num_terms <= 0 ) return 0;
	bipoly_sort_mpz (phi->mpz_terms, phi->num_terms, sort_v);  phi->mpz_sort_v = sort_v;
	phi->deg[0] = phi->deg[1] = -1;
	for ( k = 0 ; k < phi->num_terms ; k++ ) for ( i = 0 ; i < 2 ; i++ ) if ( phi->mpz_terms[k].e[i] > phi->deg[i] ) phi->deg[i] = phi->mpz_terms[k].e[i];
	phi->ff_terms = malloc (sizeof(*phi->ff_terms)*phi->num_terms);
	phi->phi_z = malloc (max(phi->deg[0],phi->deg[1])*phi->num_terms*sizeof(*phi->phi_z));
	phi->p = 0;
	dbg_printf ("Loaded bipoly of bi-degree (%d,%d) with %d terms  for from file %s\n", phi->deg[0], phi->deg[1], phi->num_terms, filename);
	return phi->num_terms;
}

void bipoly_clear (bipoly_t phi)
{
	int k;
	
	for ( k = 0 ; k < phi->num_terms ; k++ )  mpz_clear (phi->mpz_terms[k].c);
	mem_free (phi->mpz_terms);
	mem_free (phi->ff_terms);
	mem_free (phi->phi_z);
}

void bipoly_reduce (bipoly_t phi)
{
	if ( phi->p == _ff_p ) return;
	bipoly_reduce_ff (phi->ff_terms, phi->mpz_terms, phi->num_terms);
	phi->ff_sort_v = phi->mpz_sort_v;
}

int bipoly_eval (bipoly_t phi, int v, ff_t z)
{
	if ( phi->p != _ff_p ) bipoly_reduce (phi);
	v = ( v ? 1 : 0 );			// treat non-zero as 1
	if ( v == phi->v && _ff_equal (phi->z, z) ) return phi->phi_z_d;
	phi->v = v;
	if ( phi->v == phi->ff_sort_v )	// sort on variable *not* being instantiated -- this should happen very rarely if mpz terms are sorted properly
		{ phi->ff_sort_v = 1-phi->v; bipoly_sort_ff (phi->ff_terms, phi->num_terms, phi->ff_sort_v); }
	_ff_set (phi->z, z);
	bipoly_eval_ff (phi->phi_z, phi->deg[1-v], phi->ff_terms, phi->num_terms, phi->v, phi->z);
	phi->phi_z_d = ff_poly_degree (phi->phi_z, phi->deg[1-v]);
	if ( dbg_level >= DEBUG_LEVEL ) {
		if ( v ) printf ("Phi_f(x,%ld) = ", _ff_get_ui(z)); else  printf ("Phi_f(%ld,x) = ", _ff_get_ui(z));
		ff_poly_print (phi->phi_z, phi->phi_z_d);
	}
	return phi->phi_z_d;
}

// assumes r is big enough (max bi-degree suffices)
int bipoly_roots (ff_t r[], bipoly_t phi, int v, ff_t z)
{
	bipoly_eval (phi, v, z);
	ff_poly_monic (phi->phi_z, 0, phi->phi_z, phi->phi_z_d);
	return ff_poly_roots (r, phi->phi_z, phi->phi_z_d);
}


#define issign(c)		((c)=='+'||(c)=='-')
#define isnumeric(c)	(isdigit(c)||issign(c))
#define iscaret(c)		((c)=='^')
#define ismult(c)		((c)=='*')

char buf[65536];
char buf2[65536];

int bipoly_load_mpz (bipoly_mpz_t **pphi, char *file)
{
	static mpz_t c;
	static int init;
	bipoly_mpz_t *phi;
	FILE *fp;
	char v1, v2;
	char *s, *se;
	int i, j, k, m, sign;

	if ( ! init ) { mpz_init(c); init = 1; }
	fp = fopen (file, "r");
	if ( ! fp ) { printf ("Error loading bipoly from file %s\n", file); return -1; }
	info_printf ("Bipoly loading %s\n", file);
	v1 = v2 = 0;
	m = 0;
	while ( fgets(buf, sizeof(buf), fp) ) {
		s = buf;
		for (;;) {
			while ( isspace(*s) ) s++;
			if ( ! *s ) break;
			if ( issign(*s) ) { sign=1; s++; } else { sign = 0; }
			while ( isspace(*s) ) s++;
			if ( ! *s && sign ) { if ( ! fgets(buf,sizeof(buf),fp) ) { printf ("Expected another line at the end of bipoly file %s\n", file); abort(); } for ( s = buf ; isspace(*s) ; s++ ); }
			while ( isdigit(*s) ) s++;
			while ( isspace(*s) ) s++;
			for ( k = 0 ; ismult(*s) || isalpha(*s) ; k++ ) {
				if ( k > 1 ) { printf ("too many mults in bipoly load for file %s\n", file);  fclose(fp); return -1; }
				if ( ismult(*s) ) { s++;  while ( isspace(*s) ) s++; }
				if ( ! isalpha(*s) ) { printf ("Unexpected character %s in bipoly file %s\n", s, file); fclose(fp); return -1; }
				if ( ! v1 ) {
					v1 = *s;
				} else if ( *s != v1 ) {
					if ( ! v2 ) {
						v2 = *s;
					} else {
						if ( *s != v2 ) { printf ("Unexpected character %s in bipoly file %s\n", s, file); fclose(fp); return -1; }
					}
				}
				s++;
				while ( isspace(*s) ) s++;
				if ( iscaret(*s) ) { s++; while (isdigit(*s) ) s++; }
				while ( isspace(*s) ) s++;
			}
			m++;
		}
	}
	rewind(fp);
	phi = malloc(m*sizeof(*phi));
	i = 0;
	while ( fgets(buf, sizeof(buf), fp) ) {
		s = buf;
		for (;;) {
			while ( isspace(*s) ) s++;
			if ( ! *s ) break;
			if ( i >= m ) { printf ("miscounted terms in bipoly file %s\n", file); abort(); }
			mpz_init(phi[i].c);
			if ( issign(*s) ) { sign = ( *s++ == '-' ? -1 : 1 ); } else { sign = 0; }
			while ( isspace(*s) ) s++;
			if ( ! *s && sign ) { if ( ! fgets(buf,sizeof(buf),fp) ) { printf ("Expected another line at the end of bipoly file %s\n", file); abort(); } for ( s = buf ; isspace(*s) ; s++ ); }
			if ( isdigit(*s) ) {
				for ( se = s ; isdigit(*se) ; se++ );
				memcpy (buf2,s,se-s); buf2[se-s]='\0';
				mpz_set_str(phi[i].c, buf2, 10);
				s = se;
			} else {
				mpz_set_ui(phi[i].c,1);
			}
			if ( sign < 0 ) mpz_neg(phi[i].c,phi[i].c);
			while ( isspace(*s) ) s++;
			phi[i].e[0] = phi[i].e[1] = 0;
			for ( k = 0 ; ismult(*s)|| isalpha(*s) ; k++ ) {
				if ( k > 1 ) { printf ("too many mults in bipoly load for file %s\n", file); abort(); }
				if ( ismult(*s) ) { s++;  while ( isspace(*s) ) s++; }
				if ( *s == v1 ) { j = 0; } else if ( *s== v2 ) { j = 1; } else { printf ("Unexpected character %c in bipoly file %s\n", *s, file); exit(0); }
				phi[i].e[j] = 1;
				s++;
				while ( isspace(*s) ) s++;
				if ( iscaret(*s) ) {
					s++;
					for ( se = s ; isdigit(*se) ; se++ );
					memcpy(buf2,s,se-s); buf2[se-s]='\0';
					phi[i].e[j] = atoi(buf2);
					s = se;
				}
				while ( isspace(*s) ) s++;
			}
			if ( dbg_level >= DEBUG_LEVEL ) gmp_printf ("Loaded X^%dY^%d coefficient %Zd\n", phi[i].e[0], phi[i].e[1],phi[i].c);
			i++;
		}
	}
	fclose(fp);
	if ( i != m ) { err_printf ("miscounted terms in bipoly file %s i=%d, t=%d\n", file, i, m); abort(); }
	bipoly_sort_mpz (phi, m, 0);	// sort on first variable encountered by default
	*pphi = phi;
	dbg_printf ("Loaded %d terms from %s\n", m, file);
	return m;
}

int bipoly_str_to_mpz (bipoly_mpz_t **pphi, bipoly_str_t Phi[], int m)
{
	bipoly_mpz_t *phi;
	register int i;
	
	phi = malloc(m*sizeof(*phi));
	for ( i = 0 ; i < m ; i++ ) {
		phi[i].e[0] = Phi[i].e[0];
		phi[i].e[1] = Phi[i].e[1];
		mpz_init_set_str (phi[i].c, Phi[i].c, 0);
	}
	*pphi = phi;
	return m;
}

void bipoly_free_mpz (bipoly_mpz_t *phi, int m)
{
	int i;
	
	for ( i = 0 ; i < m ; i++ ) mpz_clear (phi[i].c);
	free(phi);
}

void bipoly_print_mpz (bipoly_mpz_t phi[], int m, char v1, char v2)
{
	static mpz_t c;
	static int init;
	int i, sign;
	
	if ( ! init ) { mpz_init(c); init = 1; }
	for ( i = 0 ; i < m ; i++ ) {
		if ( mpz_sgn(phi[i].c) < 0 ) { sign=-1; mpz_neg(c,phi[i].c); } else { sign=1; mpz_set(c,phi[i].c); }
		if ( i ) printf (" ");
		if ( i && sign > 0 ) printf ("+ ");
		if ( sign < 0 ) printf ("- ");
		if ( mpz_cmp_ui(c,1) != 0 ) {
			gmp_printf("%Zd", c);
			if ( phi[i].e[0] || phi[i].e[1] ) printf ("*");
		} else {
			if ( ! phi[i].e[0] && ! phi[i].e[1] ) printf ("1");
		}
		if ( phi[i].e[0] ) {
			printf("%c",v1);
			if ( phi[i].e[0] > 1 ) printf ("^%d", phi[i].e[0]);
			if ( phi[i].e[1] ) printf("*");
		}
		if ( phi[i].e[1] ) {
			printf("%c",v2);
			if ( phi[i].e[1] > 1 ) printf ("^%d", phi[i].e[1]);
		}
	}
	puts("");
}

/*
	For convenience when evaluating, we sort terms by increasing degree in the sorted variable
	and decreasing degree in the unsorted variable.  Note that evaluation can only be performed
	using the unsorted variable.
*/

int _qsort_bipoly_mpz_cmp_v0 (const void *a, const void *b)
{
	bipoly_mpz_t *A, *B;
	
	A = (bipoly_mpz_t *)a;  B= (bipoly_mpz_t *) b;
	if ( A->e[0] < B->e[0] ) return -1;
	if ( A->e[0] > B->e[0] ) return 1;
	if ( A->e[1] < B->e[1] ) return 1;
	if ( A->e[1] > B->e[1] ) return -1;
	return 0;
}

int _qsort_bipoly_mpz_cmp_v1 (const void *a, const void *b)
{
	bipoly_mpz_t *A, *B;
	
	A = (bipoly_mpz_t *)a;  B= (bipoly_mpz_t *) b;
	if ( A->e[1] < B->e[1] ) return -1;
	if ( A->e[1] > B->e[1] ) return 1;
	if ( A->e[0] < B->e[0] ) return 1;
	if ( A->e[0] > B->e[0] ) return -1;
	return 0;
}

void bipoly_sort_mpz (bipoly_mpz_t phi[], int t, int v)
{
	if ( ! v ) {
		qsort (phi, t, sizeof(phi[0]), _qsort_bipoly_mpz_cmp_v0);
	} else {
		qsort (phi, t, sizeof(phi[0]), _qsort_bipoly_mpz_cmp_v1);
	}
}


/*
	Sort in increasing order on sorted variable, and then in decreasing order on the other variable (so we can apply Horner's method)
*/
int _qsort_bipoly_ff_cmp_v0 (const void *a, const void *b)
{
	bipoly_ff_t *A, *B;
	
	A = (bipoly_ff_t *)a;  B= (bipoly_ff_t *) b;
	if ( A->e[0] < B->e[0] ) return -1;
	if ( A->e[0] > B->e[0] ) return 1;
	if ( A->e[1] < B->e[1] ) return 1;
	if ( A->e[1] > B->e[1] ) return -1;
	return 0;
}

int _qsort_bipoly_ff_cmp_v1 (const void *a, const void *b)
{
	bipoly_ff_t *A, *B;
	
	A = (bipoly_ff_t *)a;  B= (bipoly_ff_t *) b;
	if ( A->e[1] < B->e[1] ) return -1;
	if ( A->e[1] > B->e[1] ) return 1;
	if ( A->e[0] < B->e[0] ) return 1;
	if ( A->e[0] > B->e[0] ) return -1;
	return 0;
}

void bipoly_sort_ff (bipoly_ff_t phi[], int t, int v)
{
	if ( ! v ) {
		qsort (phi, t, sizeof(phi[0]), _qsort_bipoly_ff_cmp_v0);
	} else {
		qsort (phi, t, sizeof(phi[0]), _qsort_bipoly_ff_cmp_v1);
	}
}

void bipoly_reduce_ff (bipoly_ff_t phi[], bipoly_mpz_t PHI[], int t)
{
	register int i;
	
	for ( i = 0 ; i < t ; i++ )
		{ _ff_set_mpz(phi[i].c,PHI[i].c); phi[i].e[0] = PHI[i].e[0]; phi[i].e[1] = PHI[i].e[1];  /*gmp_printf ("Reduced X^%dY^%d coefficient from %Zd to %ld\n", phi[i].e[0], phi[i].e[1], PHI[i].c, _ff_get_ui(phi[i].c));*/ }
}


void bipoly_print_ff (bipoly_ff_t phi[], int m, char v1, char v2)
{
	ff_t c;
	int i, sign;
	
	for ( i = 0 ; i < m ; i++ ) {
		_ff_set (c, phi[i].c);
		sign = 1;
		if ( i ) printf (" ");
		if ( i && sign > 0 ) printf ("+ ");
		if ( sign < 0 ) printf ("- ");
		if ( ! _ff_one(c) ) {
			printf("%ld", _ff_get_ui(c));
			if ( phi[i].e[0] || phi[i].e[1] ) printf ("*");
		} else {
			if ( ! phi[i].e[0] && ! phi[i].e[1] ) printf ("1");
		}
		if ( phi[i].e[0] ) {
			printf("%c",v1);
			if ( phi[i].e[0] > 1 ) printf ("^%d", phi[i].e[0]);
			if ( phi[i].e[1] ) printf("*");
		}
		if ( phi[i].e[1] ) {
			printf("%c",v2);
			if ( phi[i].e[1] > 1 ) printf ("^%d", phi[i].e[1]);
		}
	}
	puts("");
}

int bipoly_eval_ff (ff_t f[], int n, bipoly_ff_t phi[], int m, int v, ff_t x)
{
	ff_t y; 
	register int i, j, d, flag;

	if ( v<0 || v>1 ) { printf ("invalid variable index %d in bipoly_eval!\n", v); abort(); }
	d = phi[m-1].e[1-v];	
	if ( d > n ) { printf ("overflow in bipoly_eval, degree %d exceeds %d (m=%d, v=%d, x=%ld)\n", d, n, m, v, _ff_get_ui(x)); abort(); }
	n = d;
	for ( i = 0 ; i < n ; i++ ) _ff_set_zero(f[i]);
	j = phi[0].e[1-v];  _ff_set_zero(y);
	for ( i = 0 ; i < m ; i++ ) {
//printf ("i=%d,j=%d, %ld*X^%dY^%d\n", i, j, _ff_get_ui(phi[i].c), phi[i].e[0], phi[i].e[1]);
		_ff_addto(y,phi[i].c);
//printf ("y=%ld\n", _ff_get_ui(y));
		flag = ( i==m-1 || phi[i+1].e[1-v] != j );					// flag is set if we are switching degrees in the sorted (unevaluated) variable
//printf ("flag = %d\n", flag);
		if ( flag ) d = phi[i].e[v]; else d = phi[i].e[v] - phi[i+1].e[v];
//printf ("d = %d\n", d);
		while ( d > 0 ) { ff_mult(y,y,x); d--; }
//printf ("y = %ld\n", _ff_get_ui(y));
		if ( flag ) { _ff_set(f[j],y);  /*printf ("Set x^%d coeff in f to %ld\n", j, _ff_get_ui(y));*/ _ff_set_zero(y); if ( i<m-1 ) j = phi[i+1].e[1-v]; }
	}
	return n;
}


int bipoly_eval_mod_mpz (mpz_t f[], int n, bipoly_mpz_t phi[], int m, int v, mpz_t x, mpz_t P)
{
	static mpz_t y; 
	static int init;
	register int i, j, d, flag;

	if ( ! init ) { mpz_init(y); init = 1; }
	if ( v<0 || v>1 ) { printf ("invalid variable index %d in bipoly_eval!\n", v); abort(); }
	d = phi[m-1].e[1-v];
	if ( d > n ) { printf ("overflow in bipoly_eval, degree %d exceeds %d\n", d, n); abort(); }
	n = d;
	for ( i = 0 ; i < d ; i++ ) mpz_set_ui(f[i],0);
	j = phi[0].e[1-v];  mpz_set_ui(y,0);
	for ( i = 0 ; i < m ; i++ ) {
		mpz_add(y,y,phi[i].c);  mpz_mod(y,y,P);
		flag = ( i==m-1 || phi[i+1].e[1-v] != j );					// flag is set if we are switching degrees in the sorted (unevaluated) variable
		if ( flag ) d = phi[i].e[v]; else d = phi[i].e[v] - phi[i+1].e[v];
		while ( d > 0 ) { mpz_mul(y,y,x); mpz_mod(y,y,P); d--; }
		if ( flag ) { mpz_set(f[j],y); mpz_set_ui(y,0); if ( i<m-1 ) j = phi[i+1].e[1-v]; }
	}
	return n;
}
