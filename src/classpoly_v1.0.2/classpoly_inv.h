#ifndef _CLASSPOLY_INV_INCLUDE_
#define _CLASSPOLY_INV_INCLUDE_

#include "ff_poly.h"
#include "class_inv.h"
#include "classpoly.h"

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

typedef struct classpoly_inv_struct {
	int inv;
	int p1, p2;
	int inverted, negated, units;
	int ambiguous_H, ambiguous_sign;
	mpz_t T;
	int use_trace;
	long Nrep[IQ_MAX_GENS];
} classpoly_inv_t[1];

int classpoly_inv_setup (classpoly_inv_t Hinv, classpoly_t H, int inv, classpoly_crt_t crt);
void classpoly_inv_clear (classpoly_inv_t Hinv);

void _ff_classpoly_inv_from_j (ff_t *F, ff_t *J, classpoly_inv_t Hinv, classpoly_t H);
static inline void ff_classpoly_inv_from_j (ff_t *F, ff_t *J, classpoly_inv_t Hinv, classpoly_t H)
	{ if ( Hinv->p1 > 1 ) _ff_classpoly_inv_from_j (F,J,Hinv,H); else ff_inv_from_j (F,J,Hinv->inv); }
	
int _ff_classpoly_inv_disambiguate (classpoly_t H, classpoly_inv_t Hinv);
static inline int ff_classpoly_inv_disambiguate (classpoly_t H, classpoly_inv_t Hinv)
	{ if ( Hinv->ambiguous_H || Hinv->ambiguous_sign ) return _ff_classpoly_inv_disambiguate (H, Hinv); else return 1; }

#endif
