#ifndef _CLASS_INV_H_
#define _CLASS_INV_H_

#include <math.h>
#include <gmp.h>
#include "ff_poly.h"
#include "bipoly.h"
#include "qform.h"
#include "iqclass.h"

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

#define INV_J			0
#define INV_F			1							// Weber-f
#define INV_F2		2
#define INV_F3		3
#define INV_F4		4
#define INV_G2		5							// gamma_2
#define INV_W2W3   	6
#define INV_W7E2		7
#define INV_F8		8
#define INV_W3W3   	9
#define INV_W2W5   	10
#define INV_T			11							// Ramanujan
#define INV_T2		12
#define INV_T6		13
#define INV_W2W7   	14
#define INV_W3W5   	15
#define INV_U			16							// lambda^{1/8}, modpolys are antisymmetric, not currently supported
#define INV_U2		17							// lambda^{1/4}
#define INV_U8		18							// lambda
#define INV_W3E2		19							// 6th root of the canonical power, not a class invariant
#define INV_W5E2		20							//
#define INV_W3W7   	21
#define INV_W2W11E2	22
#define INV_W2W3E2   	23
#define INV_W2W5E2   	24
#define INV_W5W5		25
#define INV_W2W13     	26
#define INV_W2W7E2   	27
#define INV_W3W3E2   	28
#define INV_W3W11E2	33
#define INV_W2W17	34
#define INV_MAX_SPECIFIC_ENUM		28
#define INV_MAX_SPECIFIC			34				// increase as required, just needs to stay below INV_ATKIN
static char inv_specific_double_eta[INV_MAX_SPECIFIC+1] = { 0,0,0,0,0,0,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };

#define INV_ATKIN					100				// general atkin A_N are specified as 100+N where N < 300
#define INV_ATKIN_END				400
#define INV_ATKIN_MAX_N			300
#define INV_SINGLE_ETA				400
#define INV_SINGLE_ETA_END			500
#define INV_SINGLE_ETA_MAX_N		100
#define INV_DOUBLE_ETA				500				// canonical double eta-quotients w_{p1,p2}^s are specified as 500 + p1*p2 where p1*p2 < 500 
#define INV_DOUBLE_ETA_END			1000			// YOU MUST MODIFY inv_composite_level if you increase this
#define INV_DOUBLE_ETA_MAX_P		241
#define INV_DOUBLE_ETA_MAX_DEG	726
#define INV_MAX					1500

#define PHI_MAX_FDEG				799
#define PHI_MAX_JDEG				24

// no longer used
#define QR71_MASK	 0x8B2394B8B58ECBBFUL	// bit i is set iff i-1is a residue mod 71
#define QR47_MASK  0x4351B2753DFUL			// bit i is set iff i is a residue mod 47 (i < 47)
#define QR59_MASK  0x22B62183E7B92BBUL		// bit i is set iff i is a residue mod 59 (i < 59)
#define QR41_MASK	 0x7B382B50737UL			// bit i is set iff i is a residue mod 41 (i < 41)
#define QR13_MASK  0x161BUL				// bit i is set iff i is a residue mod 13 (i < 13)
#define QR7_MASK	 0x17UL					// bit i is set iff i is a residue mod 7 (i < 7)

extern int phi_fj_inv, phi_fj_fdeg, phi_fj_jdeg;

int inv_pick_invariant (long D);

static inline int inv_weber(int inv) { return ( inv==INV_F || inv==INV_F2 || inv==INV_F3 || inv==INV_F4 || inv==INV_F8 ? 1 : 0 ); }
static inline int inv_atkin(int inv) { return ( inv >= INV_ATKIN && inv < INV_ATKIN_END ? 1 : 0 ); }
static inline int inv_single_eta(int inv) { return ( inv==INV_W3E2 || inv==INV_W5E2 || inv==INV_W7E2 || (inv >= INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END) ? 1 : 0 ); }
static inline int inv_double_eta(int inv) { if ( inv > 0 && inv <= INV_MAX_SPECIFIC && inv_specific_double_eta[inv] )  return 1; return ( (inv >= INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END) ? 1 : 0 ); }

// edge invariants have two isogenous j's related to each f (this applies also to single eta quotients under the Atkin-Lehner involution)
static inline int inv_edge_invariant (int inv) { return ( inv_atkin(inv) || inv_single_eta(inv) || inv_double_eta(inv) ? 1 : 0 );  }

// invariants that are square roots (which includes 4th roots, 8th roots, etc...) may require special handling due to the sign ambiguity
// note, we could treat INV_WpWqE2 without ambiguity in cases where INV_WpWq is an invariant (since we could compute this and square it), but we don't bother
static inline int inv_sqrt_invariant (int inv)
{
	switch (inv) {
	case INV_F: case INV_F3: case INV_T: case INV_W3E2: case INV_W7E2: case INV_W2W3E2: case INV_W2W7E2: case INV_W2W3:
	case INV_W2W5: case INV_W2W7: case INV_W3W3: case INV_W2W13: case INV_W3W7: case INV_W3W11E2: case INV_W2W11E2: return 1;
	}
	return 0;
}

static inline int inv_units (int inv)
	{ return ( inv_double_eta(inv) || inv_weber(inv) ? 1 : 0 ); }

static inline int inv_inverted_involution (int inv) { return inv_double_eta(inv); }
// determined by trial and error
static inline int inv_negated_involution (int inv) { return ( inv==INV_F || inv==INV_W3W5 || inv==INV_W3W7|| inv==INV_W3W3 || inv==INV_DOUBLE_ETA+15 || inv==INV_DOUBLE_ETA+35 ? 1 : 0 ); }

static inline int inv_level (int inv)		// actually just returns the square-free part of the level, which is all we care about
{
	switch (inv) {
	case INV_J: return 1;
	case INV_G2: case INV_W3W3E2: return 3;
	case INV_F8: case INV_F4: case INV_F2: case INV_F: case INV_T: case INV_T2: case INV_T6: return 6;
	case INV_F3: case INV_U: case INV_U2: case INV_U8: return 2;
	case INV_W3E2: case INV_W3W3: return 6;
	case INV_W5E2: case INV_W5W5: return 15;
	case INV_W7E2: case INV_W2W7E2: case INV_W2W7: return 14;
	case INV_W3W5: return 15;
	case INV_W2W3E2: case INV_W2W3: return 6;
	case INV_W2W5E2: case INV_W2W5: return 30;
	case INV_W2W11E2: return 66;
	case INV_W2W13: return 26;
	case INV_W2W17: return 102;
	case INV_W3W7: return 42;
	 case INV_W3W11E2: return 66;
	}
	if ( inv > INV_ATKIN && inv < INV_ATKIN_END ) return inv-INV_ATKIN;
	if ( inv > INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END ) return inv-INV_SINGLE_ETA;
	if ( inv > INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) return inv-INV_DOUBLE_ETA;
	return 0;
}

// Where applicable, returns N=p1*p2 (possibly p2=1) s.t. two j's related to the same f are N-isogenous, and 0 otherwise.
// This is often (but not necessarily) equal to the level.
static inline int inv_degree (int *p1, int *p2, int inv)
{
	int p, N;
	
	*p1 = 1; *p2 = 1;
	switch (inv) {
	case INV_W3E2: return (*p1=3);
	case INV_W5E2: return (*p1=5);
	case INV_W7E2: return (*p1=7);
	case INV_W3W5: return (*p1=3)*(*p2=5);
	case INV_W2W3E2: case INV_W2W3: return (*p1=2)*(*p2=3);
	case INV_W2W5E2: case INV_W2W5: return (*p1=2)*(*p2=5);
	case INV_W2W7E2: case INV_W2W7: return (*p1=2)*(*p2=7);
	case INV_W2W11E2: return (*p1=2)*(*p2=11);
	case INV_W2W13: return (*p1=2)*(*p2=13);
	case INV_W2W17: return (*p1=2)*(*p2=17);
	case INV_W3W7: return (*p1=3)*(*p2=7);
	case INV_W3W11E2: return (*p1=3)*(*p2=11);
	case INV_W3W3E2: case INV_W3W3: return (*p1=3)*(*p2=3);
	case INV_W5W5: return (*p1=5)*(*p2=5);
	}
	if ( inv > INV_ATKIN && inv < INV_ATKIN_END ) return (*p1=inv-INV_ATKIN);
	if ( inv > INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END ) return (*p1=inv-INV_SINGLE_ETA);
	if ( inv > INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) {
		N = inv-INV_DOUBLE_ETA;  p = 0;
		if ( !(N%2) ) p = 2;
		else if ( !(N%3) ) p = 3;
		else if ( !(N%5) ) p = 5;
		else if ( !(N%7) ) p = 7;
		else if ( !(N%11) ) p = 11;
		else if ( !(N%13) ) p = 13;
		else if ( !(N%17) ) p = 17;
		else if ( !(N%19) ) p = 19;
		if ( ! p ) return 0;
		*p1 = p;
		*p2 = N/p;
		return N;
	}
	return 0;
}

// Certain invariants require that D not have 2 in it's conductor, but this doesn't apply to every invariant with even level so we handle it separately
static inline int inv_odd_conductor (int inv)
{	switch (inv) { case INV_F: case INV_U: case INV_U2: case INV_U8: case INV_W3W3: case INV_W3E2: case INV_W3W7: case INV_W7E2: return 1; } return 0; }


// returns 12/gcd(12,N-1) for N prime
static inline int inv_s(int N)
{
	switch (N%12) {
	case 1:	return 1;
	case 2:	return 12;
	case 3:	return 6;
	case 5:	return 3;
	case 7:	return 2;
	case 11:	return 6;
	}
	return 0;
}

static inline int inv_atkin_enum(int inv) 	// only true for atkin invariants where j has deg 2
{
	int N;
	
	if ( ! inv_atkin(inv) ) return 0;
	N = inv_level(inv);
	if ( N==41 || N == 47 || N == 59 || N == 71 ) return 1;
	if ( N < 32 ) return 1;
	return 0;
}

static inline int inv_single_eta_enum(int inv) // only true for eta invariants where j has deg 1
{
	int N;
	
	if ( ! inv_single_eta(inv) ) return 0;
	N = inv_level(inv);
	if ( N < 8 || N == 13 ) return 1;
	return 0;
}

static inline int inv_double_eta_enum(int inv) // only true for double eta invariants where j has deg 2
{
	int p2,N;

	if  ( ! inv_double_eta(inv) ) return 0;
	N = inv_degree(&p2,&p2,inv);
	if ( N==26 || N==35 || N==39 ) return 1;
	if ( N < 22 ) return 1;
	return 0;
}


static inline int inv_enum(int inv)
{
	// we don't yet support enumeration for INV_U do to anti-symmetry, but we could
	if ( inv <= INV_MAX_SPECIFIC_ENUM ) return (inv==INV_W5W5||inv==INV_W2W11E2||inv==INV_W3W11E2||inv==INV_W2W17?0:1);
	if ( inv > INV_ATKIN && inv < INV_ATKIN_END ) return inv_atkin_enum(inv);
	if ( inv > INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END ) return inv_single_eta_enum(inv);
	if ( inv > INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) return inv_double_eta_enum(inv);
	return 0;
}

void _inv_load_phi_fj (int inv);
static inline void inv_load_phi_fj (int inv)
	{ if ( inv != phi_fj_inv ) _inv_load_phi_fj(inv); }

static inline double inv_enum_height_factor (int inv)
{
	register int N;
	
	switch (inv) {
	case INV_J:	return 1;
	case INV_G2:   return 3;
	case INV_F:	return 72;
	case INV_F2:   return 36;
	case INV_F3:   return 24;
	case INV_F4:   return 18;
	case INV_F8:	return 9;
	case INV_T:	return 36;
	case INV_T2:   return 18;
	case INV_T6:   return 6;
	case INV_U8:	return 6;
	case INV_U2:	return 24;
	case INV_U:	return 48;
	case INV_W3E2: return 24;
	case INV_W5E2: return 18;
	case INV_W7E2: return 16;
	case INV_W2W3: return 72;
	case INV_W3W3: return 36;
	case INV_W2W5: return 54;
	case INV_W2W7: return 48;
	case INV_W3W5: return 36;
	case INV_W2W11E2: return 21.6;
	case INV_W2W17:	return 40.5;
	case INV_W2W13: return 42;
	case INV_W3W7: return 32;
	case INV_W3W11E2: return 28.5;
	case INV_W2W3E2: return 36;
	case INV_W2W5E2: return 27;
	case INV_W2W7E2: return 24;
	case INV_W3W3E2: return 18;
	case INV_W5W5: return 22.5;
	}
	if ( inv > INV_ATKIN && inv < INV_ATKIN_END ) {
		N = inv-INV_ATKIN;
		if ( N < 32 || N == 41 || N == 47 || N == 59 || N == 71 ) return (N+1)/2;
	} else if ( inv > INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END ) {
		N = inv-INV_SINGLE_ETA;
		return N+1;		// we only suppert enum invariants for single eta, they all have degree 1 in j
	} else if ( inv > INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) {
		N = inv-INV_DOUBLE_ETA;
		switch (N) {
		case 6: return 6;
		case 9: return 6;
		case 10: return 9;
		case 14: return 12;
		case 15: return 12;
		case 21: return 16;
		case 26: return 21;
		case 35: return 24;
		case 39: return 28;		
		}
	}
	return 0;
}

static inline double inv_height_factor (int inv)
{
	double hf;
	
	hf = inv_enum_height_factor(inv);
	if ( hf ) return hf;
	if ( phi_fj_inv && inv == phi_fj_inv ) return (double) phi_fj_fdeg / phi_fj_jdeg;
	printf ("Height factor for non enum inv=%d unknown, phi_fj not loaded\n", inv);
	return 1;
}

static inline int inv_sparse_factor (int inv)
{
	switch (inv) {
	case INV_G2: case INV_F8: case INV_W5E2: case INV_W3W5: case INV_W2W5E2: case INV_W3W3E2: case INV_T2: case INV_W5W5: return 3;
	case INV_F4: case INV_W3E2: case INV_W2W3E2: case INV_W2W5: case INV_W3W3: case INV_T: return 6;
	case INV_W7E2: case INV_W2W7E2: case INV_W2W13: case INV_W3W7: return 2;
	case INV_F: return 24;
	case INV_F2: case INV_W2W3: return 12;
	case INV_F3: case INV_U:return 8;
	case INV_U2: case INV_W2W7: return 4;
	}
	return 1;
}

static const char *inv_strings[INV_MAX_SPECIFIC+1] = { "j", "f", "f2", "f3", "f4", "g2", "w2w3e1", "w7e2", "f8", "w3w3e1", "w2w5e1", "t", "t2", "t6", "w2w7e1", "w3w5e1", "u", "u2", "u8",
								     "w3e2", "w5e2", "w3w7e1", "w2w11e2", "w2w3e2", "w2w5e2", "w5w5e1", "w2w13e1", "w2w7e2", "w3w3e2", "h2", "?", "?", "?", "w3w11e2", "w2w17e1" };
static char inv_string_buf[16];
								     
static inline const char *inv_string (int inv) {
	if ( 0 <= inv  && inv <= INV_MAX_SPECIFIC ) return inv_strings[inv];
	if ( inv > INV_ATKIN && inv < INV_ATKIN_END ) { sprintf(inv_string_buf, "a%d", inv-INV_ATKIN); return inv_string_buf; }
	if ( inv > INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END ) { sprintf(inv_string_buf, "w%d", inv-INV_SINGLE_ETA); return inv_string_buf; }
	if ( inv > INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) {
		int p1, p2;
		if ( inv_degree(&p1,&p2,inv) ) { sprintf(inv_string_buf, "w%dw%d", p1, p2); return inv_string_buf; }
	}
	return (char*)"?";
}

void _ff_inv_from_j (ff_t *x, ff_t *j, int inv);
void _ff_j_from_inv (ff_t *j, ff_t *x, int inv);
static inline void ff_inv_from_j (ff_t *x, ff_t *j, int inv) { if ( inv) _ff_inv_from_j(x,j,inv); else _ff_set(*x,*j); }
static inline void ff_j_from_inv (ff_t *j, ff_t *x, int inv) { if (inv) _ff_j_from_inv(j,x,inv); else _ff_set(*j,*x); }
int ff_inv_from_2j (ff_t x[2], ff_t *j1, ff_t *j2, int inv);		// may return up to 2 invariants
int ff_j_from_2inv (ff_t *j, ff_t *x1, ff_t *x2, int inv);			// this function returns 0 on failure and can be used to test an edge
void ff_2j_from_inv (ff_t j[2], ff_t *x, int inv);
int ff_multi_j_from_inv (ff_t j[], ff_t *x, int inv, int n, int desired);

void mpz_j_from_inv (mpz_t J, mpz_t X, mpz_t P, int inv);	// this function is only available if class_inv_mpz.o and the zp_poly library are linked in

static inline int inv_good_invariant (int inv)
{
	if ( inv >= 0 && inv <= INV_MAX_SPECIFIC ) return 1;
	if ( inv > INV_ATKIN && inv < INV_ATKIN_END ) {
		register int N;
		N = inv-INV_ATKIN;
		if ( N > INV_ATKIN_MAX_N || N==2 || ! ui_is_prime(N) ) return 0;
		return 1;
	}
	if ( inv > INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END ) {
		register int N;
		N = inv-INV_SINGLE_ETA;
		if ( N > INV_SINGLE_ETA_MAX_N || N==2 || ! ui_is_prime(N) ) return 0;
		return 1;
	}
	if ( inv > INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) {
		int p1, p2;
		if ( ! inv_degree (&p1, &p2, inv) ) return 0;
		if ( ! ui_is_prime(p2) ) return 0;		// rely on inv_degree to check p1
		return 1;
	}
	return 0;
}


static inline int inv_class_invariant (int inv)
{
	if ( ! inv_good_invariant(inv) ) return 0;
	if ( inv > INV_MAX_SPECIFIC ) return 1;
	if ( inv == INV_U || inv == INV_U2 || inv == INV_U8 || inv == INV_W7E2 || inv == INV_W3E2 ) return 0;
	return 1;
}

static inline int inv_ramified (long D, int inv)		// assumes D is a good discriminant for inv, which implies that the level is prime to the conductor
{
	int p1,p2,N;
	
	N = inv_degree(&p1,&p2,inv);
	if ( N <= 1 ) return 0;
	return ( !(D%p1) && !(D%p2) ? 1 : 0);
}

// this has gotten a bit big to inline, but its convenient to have it in the header file to avoid forcing code to link with class_inv.o
static inline int inv_good_discriminant (long D, int inv)	// assumes D is a discriminant, returns true if the specified function yields class invariants for this D
{
	long a,b,c;
	int N;

	if ( ! inv ) return 1;
	switch (inv) {
	case INV_G2: 	if ( ! (D%3) ) return 0; else return 1;
	case INV_F: case INV_F2: case INV_F3: case INV_F4: case INV_F8:   if ( ! (D%3) ) return 0;  return ( ((-D)&7)==7 ? 1 : 0 );
	case INV_T: case INV_T2: case INV_T6: return ( ((-D)%24)==11 ? 1 : 0 );
	case INV_W5E2: if ( !(D%3) ) return 0; inv = INV_SINGLE_ETA+5; break;
	case INV_W3W5: if ( !(D%3) ) return 0; inv = INV_DOUBLE_ETA+15; break;
	case INV_W3W3E2: if ( !(D%3) ) return 0; inv = INV_DOUBLE_ETA+9; break;
	case INV_W3W3: if ( !(D&1) || !(D%3) ) return 0; inv = INV_DOUBLE_ETA+9; break;
	case INV_W2W3E2: if ( !(D%3) ) return 0;  inv = INV_DOUBLE_ETA+6; break;
	case INV_W2W3: if ( ((-D)&7)!=7 ) return 0; if ( !(D%3) ) return 0;  inv = INV_DOUBLE_ETA+6; break;
	case INV_W2W5: if ( ((-D)%80)==20 ) return 0; if ( !(D%3) ) return 0;  inv = INV_DOUBLE_ETA+10; break;
	case INV_W2W5E2: if ( !(D%3) ) return 0;  inv = INV_DOUBLE_ETA+10; break;
	case INV_W2W7E2: if ( ((-D)%112)==84 ) return 0; inv = INV_DOUBLE_ETA+14; break;
	case INV_W2W7: if ( ((-D)&7)!=7 ) return 0; inv = INV_DOUBLE_ETA+14; break;
	case INV_W2W11E2:  if ( ((-D)&7)!=7 ) return 0;  if ( !(D%3) ) return 0;  inv = INV_DOUBLE_ETA+22; break;
	case INV_W2W13: if ( ((-D)%208)==52 ) return 0; inv = INV_DOUBLE_ETA+26; break;
	case INV_W2W17: if ( !(D%3) ) return 0;  inv = INV_DOUBLE_ETA+34; break;
	case INV_W3W7: if ( !(D&1) || !((-D)%21) ) return 0; inv = INV_DOUBLE_ETA+21; break;
	case INV_W3W11E2: if ( !(D&1) || !(D%3) ) return 0; inv = INV_DOUBLE_ETA+33; break;
	case INV_W5W5: if ( ! (D%3) ) return 0; inv = INV_DOUBLE_ETA+25; break;
	}
	if ( inv > INV_ATKIN && inv < INV_ATKIN_END ) {
		N = inv - INV_ATKIN;
		if ( ! qform_primeform(&a,&b,&c,N,D,-1) ) return 0;	// this checks that N is non-inert and prime to the conductor
		if ( qform_is_identity(a,b,c) ) return 0;				// don't allow N to be principal
		if ( (D%N) && qform_is_2torsion(a,b,c) ) return 0;	// if N is unramified, require it to have order > 2
		return 1;
	}
	if ( inv > INV_SINGLE_ETA && inv < INV_SINGLE_ETA_END ) {
		N = inv - INV_SINGLE_ETA;
		if ( ! qform_primeform(&a,&b,&c,N,D,-1) ) return 0;	// this checks that N is non-inert and prime to the conductor
		if ( qform_is_identity(a,b,c) ) return 0;				// don't allow N to be principal
		if ( (D%N) ) return 0;							// require N ramified 
		return 1;
	}
	if ( inv > INV_DOUBLE_ETA && inv < INV_DOUBLE_ETA_END ) {
		int i1, i2, p1, p2;

		// By Corollary 3.1 of Enge-Schertz Constructing elliptic curves over finite fields using double eta-quotients,
		// we need p1 != p2 to both be non-inert and prime to the conductor, and if p1=p2=p we want p split and prime to the conductor
		// we exclude the case that p1=p2 divides the conductor, even though this does yield class invariants
		if ( ! (N=inv_degree (&p1,&p2,inv)) ) return 0;
		i1 = kronecker_p(D,p1);   if ( i1 < 0 ) return 0;
		if ( p1 == p2 && !i1 ) return 0;								// exclude ramified case for w_{p,p}
		i2 = kronecker_p(D,p2);   if ( i2 < 0 ) return 0;
		if ( ! qform_primeform(&a,&b,&c,p1,D,-1) ) return 0;				// this also verifies that p1 is prime to the conductor
		if ( qform_is_identity(a,b,c) ) return 0;							// don't allow p1 to be principal
		if ( i1 && qform_is_2torsion(a,b,c) ) return 0;						// if p1 is unramified, require it to have order > 2
		if ( p1==p2 ) {
			qform_square(&a,&b,&c,a,b,c,ceil(pow(-D,0.25)));
			if ( qform_is_2torsion(a,b,c) ) return 0;						// if p1=p2 we need p1*p1 to be distinct from its inverse
			return 1;
		}
		if ( ! qform_primeform(&a,&b,&c,p2,D,-1) ) return 0;				// this also verifies that p2 is prime to the conductor
		if ( qform_is_identity(a,b,c) ) return 0;							// don't allow p2 to be principal
		if ( i2 && qform_is_2torsion(a,b,c) ) return 0;						// if p2 is unramified, require it to have order > 2
		if ( i1 > 0 && i2 > 0 ) if ( ! qform_p1p2_check(p1,p2,D) ) return 0;	// if p1 and p2 are split, we also require p1*p2, p1*p2^-1,p1^-1*p2, and p1^-1*p2^-1 to be distinct
		if ( !i1 && !i2 ) {
			if ( qform_nform_is_principal (N,D) ) return 0;					// if both p1 and p2 are ramified, make sure their product is not principal
		}
		return 1;
	}
	return 0;
}


static inline int inv_good_discriminant_for_modpoly (long D, int inv)
{
	register int i;
	
	// here we handle functions which do not necessarily yield class invariants but for which the modpoly code works anyway
	// also, modpoly doesn't want N ramified, but for single eta-quotients this is necessary to get class polys over Z, so we need different handling
	switch (inv) {
	case INV_W3E2: if ( ((-D)%3) !=2 ) return 0; if ( ((-D)&7)!=7 ) return 0; return 1;			// doesn't yield class invariants
	case INV_W7E2: if ( ((-D)&7)!=7 ) return 0; i = (-D)%7; return ( i==3 || i==5 || i==6 ? 1 : 0 );	// doesn't yield class invariants
	case INV_U8: case INV_U2: case INV_U: return ( ((-D)&7)==7 ? 1 : 0 );					// doesn't yield class invariants
	case INV_W5E2: if ( !(D%3) ) return 0;  i = (-D)%5; return ( i==1 || i==4 ? 1 : 0);			// *does* yield class polys in Z[X] if 5 is ramified, but we don't want 5 ramified for modpoly, we want it split
	}
	return inv_good_discriminant (D, inv);
}

static inline int inv_pfilter (int inv)
{
	switch (inv) {
	case INV_G2: case INV_T: case INV_T2: case INV_T6: case INV_W3E2: case INV_W3W3: case INV_W3W3E2: case INV_W5W5:
	case INV_W5E2: case INV_W3W5: case INV_W2W5: case INV_W2W3E2: case INV_W2W5E2: case INV_W2W17: return IQ_FILTER_1MOD3; 			// ensure unique cuberoots
	case INV_U: case INV_U2: case INV_W2W7: return IQ_FILTER_1MOD4;																// ensure at most two 4th/8th roots
	case INV_F: case INV_F4: case INV_F8: case INV_F2: case INV_F3: case INV_W2W3:case INV_W2W11E2: case INV_W3W11E2:  return IQ_FILTER_1MOD3 |  IQ_FILTER_1MOD4;				// ensure unique cuberoots and  at most two 4th/8th roots
	}
	return 0;
}

static inline int inv_good_prime (long p, int inv)
{
	switch (inv) {
	case INV_G2:  case INV_T6: case INV_T2: case INV_T: case INV_W2W3E2:  case INV_W5E2: case INV_W3W3: case INV_W3W3E2: case INV_W3W5:
	case INV_W2W5E2: case INV_W2W5: case INV_W5W5:  case INV_W2W17: return ( (p%3) == 2 ? 1 : 0 );
	case INV_F4: case INV_F8: if ( (p%3) != 2 ) return 0; return ((p&3)==1?0:1);
	case INV_F3: case INV_F2: case INV_F: case INV_W2W3: case INV_W3E2: case INV_W2W11E2: case INV_W3W11E2: if ( (p%3) != 2 ) return 0;  return ((p&3)==1?0:1);
	case INV_U: case INV_U2: case INV_W2W7: return ((p&3)==1?0:1);
	}
	return 1;
}

// Currently these all require _ff_p = 2 mod 3 so that gamma2 is uniquely determined
// When 2^k-th roots are required, p must also be 3 mod 4.  In this case the returned value is always the one that is a residue (but this might not be the one you want)
static inline void ff_inv_gamma2_from_j (ff_t *gamma2, ff_t *J) { ff_cbrt(gamma2,J); }
void ff_inv_f8_from_j (ff_t *f8, ff_t *J);
void ff_inv_f_from_j (ff_t *f2, ff_t *J);	
void ff_inv_t6_from_j (ff_t *t6, ff_t *j);
void ff_inv_t_from_j (ff_t *t, ff_t *j);
void ff_inv_u_from_j (ff_t *t, ff_t *j);
void ff_inv_h2_from_j (ff_t *t, ff_t *j);

static inline void ff_inv_f4_from_j (ff_t *f4, ff_t *j)
	{ ff_t f;  ff_inv_f_from_j(&f,j); _ff_square(f,f); _ff_square(*f4,f); }
static inline void ff_inv_f3_from_j (ff_t *f3, ff_t *j)
	{ ff_t f;  ff_inv_f_from_j(&f,j); _ff_square(*f3,f); ff_mult(*f3,*f3,f); }
static inline void ff_inv_f2_from_j (ff_t *f2, ff_t *j)
	{ ff_t f;  ff_inv_f_from_j(&f,j); _ff_square(*f2,f); }
static inline void ff_inv_t2_from_j (ff_t *t2, ff_t *j)
	{ ff_t t;  ff_inv_t_from_j(&t,j); _ff_square(*t2,t); }
static inline void ff_inv_u8_from_j (ff_t *u8, ff_t *j)
	{ ff_t u;  ff_inv_u_from_j(&u,j); _ff_square(u,u); _ff_square(u,u); _ff_square(*u8,u); }
static inline void ff_inv_u2_from_j (ff_t *u2, ff_t *j)
	{ ff_t u;  ff_inv_u_from_j(&u,j); _ff_square(*u2,u); }

void ff_inv_single_eta_from_j (ff_t *w, ff_t *j, int inv);
void ff_inv_double_eta_from_j (ff_t *w, ff_t *j, int inv);
void ff_inv_atkin_from_j (ff_t *a, ff_t *j, int inv);
void ff_inv_2atkin_from_j (ff_t a[2], ff_t *j, int inv);
void ff_inv_2single_eta_from_j (ff_t w[2], ff_t *j, int inv);
void ff_inv_2double_eta_from_j (ff_t w[2], ff_t *j, int inv);
int ff_inv_atkin_from_2j (ff_t *a, ff_t *j1, ff_t *j2, int inv);
int ff_inv_single_eta_from_2j (ff_t *w, ff_t *j1, ff_t *j2, int inv);
int ff_inv_double_eta_from_2j (ff_t *w, ff_t *j1, ff_t *j2, int inv);
int ff_inv_j_from_2double_eta (ff_t *j, ff_t *x1, ff_t *x2, int inv);
int ff_inv_count_for_j (ff_t j, int inv);		// counts the roots of Psi_f(X,j) -- for testing.

void ff_inv_j_from_t6 (ff_t *j, ff_t *t6);
void ff_inv_j_from_u8 (ff_t *j, ff_t *u8);
void ff_inv_j_from_atkin (ff_t *j, ff_t *a, int inv);
void ff_inv_j_from_single_eta (ff_t *j, ff_t *w, int inv);
void ff_inv_j_from_double_eta (ff_t *j, ff_t *w, int inv);
void ff_inv_2j_from_atkin (ff_t j[2], ff_t *a, int inv);
void ff_inv_2j_from_single_eta (ff_t j[2], ff_t *w, int inv);
void ff_inv_2j_from_double_eta (ff_t j[2], ff_t *w, int inv);

// the mpz_* functions below are implemented in class_inv_mpz.c
void mpz_j_from_gamma2 (mpz_t J, mpz_t G2, mpz_t P);
void mpz_j_from_f8 (mpz_t J, mpz_t F8, mpz_t P);
void mpz_j_from_f4 (mpz_t J, mpz_t F4, mpz_t P);
void mpz_j_from_f3 (mpz_t J, mpz_t F3, mpz_t P);
void mpz_j_from_f2 (mpz_t J, mpz_t F2, mpz_t P);
void mpz_j_from_f (mpz_t J, mpz_t F, mpz_t P);
void mpz_j_from_t6 (mpz_t J, mpz_t T6, mpz_t P);
void mpz_j_from_t2 (mpz_t J, mpz_t T2, mpz_t P);
void mpz_j_from_t (mpz_t J, mpz_t T2, mpz_t P);
void mpz_j_from_u8 (mpz_t J, mpz_t U8, mpz_t P);
void mpz_j_from_u2 (mpz_t J, mpz_t U8, mpz_t P);
void mpz_j_from_u (mpz_t J, mpz_t U, mpz_t P);
void mpz_j_from_atkin (mpz_t J, mpz_t A, mpz_t P, int inv);
void mpz_j_from_single_eta (mpz_t J, mpz_t A, mpz_t P, int inv);
void mpz_j_from_double_eta (mpz_t J, mpz_t A, mpz_t P, int inv);


static inline void ff_inv_j_from_gamma2 (ff_t *j, ff_t *g2)
	{ register ff_t x;   _ff_square(x,*g2);  ff_mult(*j,*g2,x); }
	
static inline void ff_inv_j_from_f8 (ff_t *j, ff_t *f8)
{
	register ff_t x,y,c16;
	ff_t g2;
	
	// gamma_2 = (f8^3-16) / f8
	_ff_invert(y,*f8);  _ff_set_ui(c16,16); _ff_square(x,*f8); ff_mult(x,x,*f8); _ff_subfrom(x,c16); _ff_mult(g2,x,y);
	ff_inv_j_from_gamma2(j,&g2);
}

static inline void ff_inv_j_from_h8 (ff_t *j, ff_t *h8)
{
	register ff_t x,y;
	ff_t g2;
	
	// gamma_2 = 4*(4*h8^3+1) / h8
	_ff_invert(y,*h8);  _ff_square(x,*h8); ff_mult(x,x,*h8); _ff_x2(x); _ff_x2(x); _ff_inc (x); _ff_x2(x);  _ff_x2(x); _ff_mult(g2,x,y);
	ff_inv_j_from_gamma2(j,&g2);
}

static inline void ff_inv_j_from_f4 (ff_t *j, ff_t *f4)
	{ ff_t f8;  _ff_square(f8,*f4); ff_inv_j_from_f8(j,&f8); }

static inline void ff_inv_j_from_f2 (ff_t *j, ff_t *f2)
	{ ff_t f8;  _ff_square(f8,*f2); ff_square(f8,f8); ff_inv_j_from_f8(j,&f8); }

static inline void ff_inv_j_from_f (ff_t *j, ff_t *f)
	{ ff_t f8; _ff_square(f8,*f); ff_square(f8,f8); ff_square(f8,f8);  ff_inv_j_from_f8(j,&f8); }
	
static inline void ff_inv_j_from_f3 (ff_t *j, ff_t *f3)
	{ ff_t f; ff_cbrt(&f,f3); ff_inv_j_from_f (j,&f); }

static inline void ff_inv_j_from_t2 (ff_t *j, ff_t *t2)
	{ ff_t t6;  _ff_square(t6,*t2); ff_mult(t6,t6,*t2); ff_inv_j_from_t6(j,&t6); }
	
static inline void ff_inv_j_from_t (ff_t *j, ff_t *t)
	{ ff_t t6;  _ff_square(t6,*t); ff_mult(t6,t6,*t); ff_square(t6,t6); ff_inv_j_from_t6(j,&t6); }

static inline void ff_inv_j_from_u2 (ff_t *j, ff_t *u2)
	{ ff_t u8; _ff_square(u8,*u2); ff_square(u8,u8);  ff_inv_j_from_u8(j,&u8); }
	
static inline void ff_inv_j_from_u (ff_t *j, ff_t *u)
	{ ff_t u8; _ff_square(u8,*u); ff_square(u8,u8); ff_square(u8,u8);  ff_inv_j_from_u8(j,&u8); }
	
static inline void ff_inv_j_from_h2 (ff_t *j, ff_t *h2)
	{ ff_t h8;  _ff_square(h8,*h2); ff_square(h8,h8); ff_inv_j_from_f8(j,&h8); }
	
static inline void ff_inv_eval_derivative_g2 (ff_t *o, ff_t *g2)
	{ ff_t t;  _ff_square(t,*g2); _ff_triple(*o,t); }

static inline void ff_inv_eval_derivative2_g2 (ff_t *o, ff_t *g2)
	{ ff_t t;  _ff_add(t,*g2,*g2); _ff_triple(*o,t); }

static inline void ff_inv_eval_derivative_f (ff_t *o, ff_t *f)
{
	register ff_t f24, t0, t1, t2;
	
	// psi(x) = (x^24-16)^3/x^24, so d/dx psi = 48(x^72 - 24*x^48 + 2048) / x^25
	_ff_square(t1,*f); _ff_mult(f24,t1,*f); ff_square(f24,f24); ff_square(f24,f24); ff_square(f24,f24);
	_ff_square(t1,f24); _ff_set_ui(t2,24); _ff_sub(t0,f24,t2); _ff_mult(t2,t0,t1); _ff_set_ui(t1,2048); _ff_addto(t2,t1);
	_ff_mult(t1,f24,*f); ff_invert(t0,t1); _ff_mult(t1,t0,t2); _ff_set_ui(t0,48); _ff_mult(*o,t0,t1);
}

static inline void ff_inv_eval_derivative2_f (ff_t *o, ff_t *f)
{
	register ff_t f24, f2, t0, t1, t2;
	
	// psi(x) = (x^24-16)^3/x^24, so d/dx d/dx psi = 48(47x^72 - 552x^48 - 51200) / x^26
	_ff_square(f2,*f); _ff_mult(f24,f2,*f); ff_square(f24,f24); ff_square(f24,f24); ff_square(f24,f24);
	_ff_square(t1,f24); _ff_set_ui(t2,47); _ff_mult(t0,t2,f24); _ff_set_ui(t2,552); _ff_subfrom(t0,t2); _ff_mult(t2,t0,t1); _ff_set_ui(t0,51200); _ff_subfrom(t2,t0);
	_ff_mult(t0,f24,f2); ff_invert(t1,t0); _ff_mult(t0,t1,t2); _ff_set_ui(t1,48); _ff_mult(*o,t0,t1);
}

static inline void ff_inv_eval_derivative (ff_t *o, ff_t *f, int inv)
{
	switch (inv) {
	case INV_J: _ff_set_one(*o); return;
	case INV_G2: ff_inv_eval_derivative_g2 (o, f); return;
	case INV_F: ff_inv_eval_derivative_f (o, f); return;
	}
	printf ("Unhandled invariant %d in ff_inv_eval_derivative\n", inv); exit (0);
}

static inline void ff_inv_eval_derivative2 (ff_t *o, ff_t *f, int inv)
{
	switch (inv) {
	case INV_J: _ff_set_zero(*o); return;
	case INV_G2: ff_inv_eval_derivative2_g2 (o, f); return;
	case INV_F: ff_inv_eval_derivative2_f (o, f); return;
	}
	printf ("Unhandled invariant %d in ff_inv_eval_derivative\n", inv); exit (0);
}

#endif
