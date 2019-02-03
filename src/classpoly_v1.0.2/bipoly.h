#ifndef _BIPOLY_INCLUDE_
#define _BIPOLY_INCLUDE_

#include <gmp.h>
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

struct mpz_bipt_struct {
	int e[2];
	mpz_t c;
};

struct ff_bipt_struct {
	int e[2];
	ff_t c;
};

struct str_bipt_struct {
	int e[2];
	char *c;
};

typedef struct mpz_bipt_struct bipoly_mpz_t;
typedef struct ff_bipt_struct bipoly_ff_t;
typedef struct str_bipt_struct bipoly_str_t;


typedef struct bipoly_struct {
	bipoly_mpz_t *mpz_terms;
	bipoly_ff_t *ff_terms;
	int mpz_sort_v;
	int ff_sort_v;
	int num_terms;
	int deg[2];					// bi-degree 
	int v;						// instantiation variable (should always be the opposite of ff_sort_v)
	ff_t z;						// instantiation value
	ff_t *phi_z;					// univariate poly in variable sort_v obtained by instantiating poly with z substituted for the other variable (1-sort_v)
	int phi_z_d;					// degree of phi_z (could be less than deg[v])
	unsigned long p;
} bipoly_t[1];

int bipoly_load (bipoly_t phi, char *filename, int sort_v);		// sort_v should be the opposite of the variable you expect to evaluate most frequently
void bipoly_clear (bipoly_t phi);
void bipoly_reduce (bipoly_t phi);						// reduce into the current finite field (as indicated by _ff_p)
int bipoly_eval (bipoly_t phi, int v, ff_t z);					// substitute z for the variable v and evaluate to obtain a univariate poly phi_z (calls bipoly_reduce first)
int bipoly_roots (ff_t r[], bipoly_t phi, int v, ff_t z);			// find all roots of the instantiaed poly phi_z (calls bipoly_eval first and makes the result monic, if required)

// low level functions

int bipoly_load_mpz (bipoly_mpz_t **pphi, char *file);
int bipoly_str_to_mpz (bipoly_mpz_t **pphi, bipoly_str_t Phi[], int m);
void bipoly_free_mpz (bipoly_mpz_t phi[], int t);
void bipoly_print_mpz (bipoly_mpz_t phi[], int m, char v1, char v2);
void bipoly_sort_mpz (bipoly_mpz_t phi[], int m, int v);
int bipoly_eval_mod_mpz (mpz_t f[], int n, bipoly_mpz_t phi[], int m, int v, mpz_t x, mpz_t P);

void bipoly_sort_ff (bipoly_ff_t phi[], int m, int v);
void bipoly_reduce_ff (bipoly_ff_t phi[], bipoly_mpz_t Phi[], int m);
int bipoly_eval_ff (ff_t f[], int n, bipoly_ff_t phi[], int m, int v, ff_t x);		// Assumes phi is sorted on the non-evaluated variable
void bipoly_print_ff (bipoly_ff_t phi[], int m, char v1, char v2);

#endif
