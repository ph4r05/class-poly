#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <ctype.h>
#include "polyparse.h"
#include "ntutil.h"
#include "cstd.h"

/*
    Copyright 2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#define CBUFSIZE					65536		// maximum length of a string constant (integer or rational)
#define UI_POLY_BUFSIZE				65536		// maximum length of a polynomial expression with coefficients that fit in an unsigned long (this needs to fit on the stack!)


static inline int isqdigit(char c) { return (isdigit(c) || c == '/' ); }

int poly_parse (void *f, int maxd, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
{
	mpq_t c;
	char cbuf[CBUFSIZE];
	char v;
	register char *r, *s, *t;
	int i, d, md, sign;

	mpq_init (c);
	for ( i = 0 ; i <= maxd ; i++ ) (*setzero)(f,i);

	// assume first alpha character is the indeterminate
	for ( s = expr ; *s && ! isalpha(*s) && *s != ']' && *s != ',' ; s++ );
	if ( ! isalpha(*s) ) {	// constant poly
		s = expr;
		if ( *s == '[' || *s == '(' ) s++;
		while ( isspace(*s) ) s++;
		if ( *s == '+' || *s == '-' ) { r=s; s++; } else r = 0;					// skip leading sign which gmp doesn't like (for reasons that completely escape me)
		for ( t = cbuf ; isqdigit(*s) ; *t++ = *s++ );
		*t='\0';
		mpq_set_str (c, cbuf, 0);
		if ( r && *r == '-' ) mpq_neg (c,c);
		if ( ! mpq_sgn(c) ) d = -1;
		else { (*addto)(f,0,c,arg); d = 0; }
		goto done;
	}
	v = *s;
	s = expr;

	if ( *s == '[' || *s == '(' ) s++;
	md = -1;
	for(;;) {
		while ( isspace(*s) ) s++;
		// handle signs explicitly to deal with things like '+ -3' which gmp doesn't like and we need to handle
		// allow double signs (but not triple signs)
		sign = 1;
		if ( *s == '-' || *s == '+' ) {
			if ( *s == '-' ) sign = -sign;
			for ( s++ ; isspace(*s) ; s++);
			if ( *s == '-' || *s == '+' ) {
				if ( *s == '-' ) sign = -sign;
				for ( s++ ; isspace(*s) ; s++);
			}
		}
		if ( isdigit(*s) ) {
			for ( t = cbuf ; isqdigit(*s) ; *t++ = *s++ );
			*t = '\0';
			if ( mpq_set_str (c, cbuf, 0) < 0 ) { d = -3; goto done; }
			mpq_canonicalize (c);
			while ( isqdigit(*s) ) s++;
			while ( isspace(*s) ) s++;			
		} else if ( *s == v ) {
			mpq_set_ui (c, 1, 1);
		} else {
			break;	// simply stop parsing when we hit something we don't know how to handle -- this could be a null terminator but need not be.
		}
		if ( sign < 0 ) mpq_neg(c,c);
		if ( *s == '*' ) { s++; while ( isspace(*s) ) s++; }
		if ( *s == v ) {
			s++;
			if ( *s == '^' ) s++;
			if ( isdigit(*s) ) {
				i = atoi (s);
				while ( isdigit(*s) ) s++;
			} else {
				i = 1;
			}
			if ( maxd >= 0 && i > maxd )  { d = -2; goto done; }
		} else {
			i = 0;
		}
		if ( ! (*addto) (f, i, c, arg) ) { d = -3; goto done; }
		if ( i > md ) md = i;
	}
	for ( d = md ; d >= 0 && (*iszero)(f,d) ; d-- );
	if ( ! d && (*iszero)(f,d) ) d = -1;
done:
	mpq_clear (c);
	return d;
}

int poly_parse_plane_quartic (void *f, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
{
	mpq_t c;
	char cbuf[CBUFSIZE];
	int off0[5] = { 0, 5, 9, 12, 14 };
	char v[3];
	register char *s, *t;
	int vd[3];
	int i, vi, d, sign;
	
	// assume first alpha character is the indeterminate
	for ( s = expr ; *s  && *s != ']' && *s != ')' && ! isalpha(*s) ; s++ );
	if ( ! isalpha(*s) ) return -2;
	v[0] = *s;
	for ( s++ ; (*s && *s != ']' && *s != ')' && ! isalpha(*s)) || *s == v[0] ; s++ );
	if ( ! isalpha(*s) ) return -2;
	if ( *s < v[0] ) { v[1] = v[0];  v[0] = *s; } else v[1] = *s;
	for ( s++ ; (*s && *s != ']' && *s != ')' && ! isalpha(*s)) || *s == v[0] || *s == v[1] ; s++ );
	if ( isalpha(*s) ) {
		if ( *s < v[0] ) { v[2] = v[1];  v[1] =  v[0];  v[0] = *s; }
		else if ( *s < v[1] ) { v[2] = v[1]; v[1] = *s; }
		else v[2] = *s;
	} else {
		v[2] = '\0';
	}
	s = expr;
	if ( *s == '[' || *s == '(' ) s++;

	mpq_init (c);
	for ( i = 0 ; i < 15 ; i++ ) (*setzero)(f,i);

	d = -2;
	for(;;) {
		while ( isspace(*s) ) s++;
		if ( ! *s || *s == ']' || *s == ')' ) break;
		// handle signs explicitly to deal with things like '+ -3' which gmp doesn't like and we need to handle
		// allow double signs (but not triple signs)
		sign = 1;
		if ( *s == '-' || *s == '+' ) {
			if ( *s == '-' ) sign = -sign;
			for ( s++ ; isspace(*s) ; s++);
			if ( *s == '-' || *s == '+' ) {
				if ( *s == '-' ) sign = -sign;
				for ( s++ ; isspace(*s) ; s++);
			}
		}
		if ( isdigit(*s) ) {
			for ( t = cbuf ; isqdigit(*s) ; *t++ = *s++ );
			*t = '\0';
			if ( mpq_set_str (c, cbuf, 0) < 0 ) break;
			mpq_canonicalize (c);
			while ( isqdigit(*s) ) s++;
			while ( isspace(*s) ) s++;			
		} else if ( *s == v[0] || *s == v[1] || (v[2] && *s == v[2]) ) {
			mpq_set_ui (c, 1, 1);
		} else {
			break;	// simply stop parsing when we hit something we don't know how to handle -- this could be a null terminator but need not be.
		}
		if ( sign < 0 ) mpq_neg(c,c);
		if ( *s == '*' ) { s++; while ( isspace(*s) ) s++; }
		
		if ( *s == v[0] || *s == v[1] || (v[2] && *s == v[2]) ) {
			vd[0]=vd[1]=vd[2]=0;
			d = 4;
			for (;;) {
				if ( !*s || (*s != v[0] && *s != v[1] && *s != v[2]) ) break;
				if ( *s == v[0] ) { vi = 0; } else if ( *s == v[1] ) { vi = 1; } else { vi = 2; }
				if ( vd[vi] ) { d = -2; break; }
				s++;  while ( isspace(*s) ) s++;
				if ( *s == '^' ) { s++; while ( isspace(*s) ) s++; }
				if ( isdigit(*s) ) {
					vd[vi] = atoi (s);
					while ( isdigit(*s) ) s++;
				} else {
					vd[vi] = 1;
				}
				while ( isspace(*s) ) s++;
				if ( *s == '*' )  { s++; while ( isspace(*s) ) s++; }
			}
			if ( d == -2 ) break;
			if ( v[2] ) {
				if ( vd[0]+vd[1]+vd[2] != 4 ) { d = -2; break; }
			} else {
				if ( vd[0]+vd[1] > 4 ) { d = -2; break; }
			}
			i = off0[vd[0]] + vd[1];
		} else {
			i = 0;
		}
		if ( ! (*addto) (f, i, c, arg) ) { d = -2; break; }
	}
	mpq_clear (c);
	return d;
}

void _nop_setzero (void *f, int i) { return; }
int _nop_addto (void *f, int i, mpq_t c, void *arg) { return 1; }
int _nop_iszero (void *f, int i) { return 0; }

int poly_parse_max_degree (char *expr) { return poly_parse (0, -1, expr, _nop_setzero, _nop_addto, _nop_iszero, 0); }

void _ui_coeff_setzero (void *f, int i) { ((unsigned long*)f)[i] = 0; }
int _ui_mod_p_coeff_addto (void *f, int i, mpq_t c, void *arg)
{
	unsigned long p, x;
	
	p = *((unsigned long *)arg);
	if (  mpz_cmp_ui (mpq_denref(c),1) == 0 ) {
		((unsigned long *)f)[i] = mpz_fdiv_ui (mpq_numref(c), p);
	} else {
		x = ui_inverse (mpz_fdiv_ui(mpq_denref(c),p), p);
		if ( ! x ) return 0;
		// use the numerater of c as an mpz_t that we can use for workspace to avoid worrying about overflow (we are allowed to trash c)
		mpz_set_ui (mpq_numref(c), mpz_fdiv_ui (mpq_numref(c), p));
		mpz_mul_ui (mpq_numref(c), mpq_numref(c), x);
		mpz_set_ui (mpq_numref(c), mpz_fdiv_ui (mpq_numref(c),p));
		mpz_add_ui (mpq_numref(c), mpq_numref(c), ((unsigned long *)f)[i]);
		((unsigned long *)f)[i] = mpz_fdiv_ui (mpq_numref(c), p);
	}
	return 1;
}
int _ui_coeff_iszero (void *f, int i) { return ((unsigned long *)f)[i] ? 0 : 1; }

int ui_poly_parse_mod_p (unsigned long f[], int maxd, char *expr, unsigned long p) { return poly_parse (f, maxd, expr, _ui_coeff_setzero, _ui_mod_p_coeff_addto, _ui_coeff_iszero, &p); }
int ui_poly_parse_plane_quartic_mod_p (unsigned long f[], char *expr, unsigned long p) { return poly_parse_plane_quartic (f, expr, _ui_coeff_setzero, _ui_mod_p_coeff_addto, _ui_coeff_iszero, &p); }

void _i_coeff_setzero (void *f, int i) { ((long*)f)[i] = 0; }
int _i_coeff_addto (void *f, int i, mpq_t c, void *unused)
{ 
	long *pc;
	
	if ( mpz_cmp_ui (mpq_denref(c),1) != 0 ) return 0;
	pc = (long*)f + i;
	if ( *pc > 0 ) mpz_add_ui (mpq_numref(c), mpq_numref(c), *pc);
	else if ( *pc < 0 ) mpz_sub_ui  (mpq_numref(c), mpq_numref(c), -(*pc));
	if ( ! mpz_fits_slong_p (mpq_numref(c)) ) return 0;
	*pc = mpz_get_si(mpq_numref(c));
	return 1;
}
int _i_coeff_iszero (void *f, int i) { return ((long *)f)[i] ? 0 : 1; }

int i_poly_parse (long f[], int maxd, char *expr) { return poly_parse (f, maxd, expr, _i_coeff_setzero, _i_coeff_addto, _i_coeff_iszero, 0); }
int i_poly_parse_plane_quartic (long f[], char *expr) { return poly_parse_plane_quartic (f, expr, _i_coeff_setzero, _i_coeff_addto, _i_coeff_iszero, 0); }


void ui_poly_print (unsigned long f[], int d_f)
{
	char buf[UI_POLY_BUFSIZE];

	ui_poly_sprint (buf, f, d_f);
	puts (buf);
}


int ui_poly_sprint (char *s, unsigned long f[], int d_f)
{
	char *t;
	int i;
	
	if ( d_f < 0 ) { strcpy (s, "[zero polynomial]");  return strlen(s); }
	t = s;
	if ( d_f >= 2 ) {
		if ( f[d_f] != 1 ) t += sprintf (t, "[%lux^%d", f[d_f], d_f); else  t += sprintf (t, "[x^%d", d_f);
	} else if ( d_f == 1 ) {
		if ( f[d_f] != 1 ) t += sprintf (t, "[%lux", f[d_f]); else  t += sprintf (t, "[x");
	} else {
		t += sprintf (t, "[%lu", f[d_f]);
	}
	for ( i = d_f-1 ; i >= 0 ; i-- ) {
		if ( f[i] ) {
			if ( i >= 2 ) {
				t += sprintf (t, " + %lux^%d", f[i], i);
			} else if ( i == 1 ) {
				t += sprintf (t, " + %lux", f[i]);
			} else {
				t += sprintf (t, " + %lu", f[i]);
			}
		}
	}
	*t++ = ']';
	*t= '\0';
	return t-s;
}

void i_poly_print (long f[], int df)
{
	char buf[UI_POLY_BUFSIZE];

	i_poly_sprint (buf, f, df);
	puts (buf);
}

int i_poly_sprint (char *s, long f[], int df)
{
	register char *t;
	register int i;
	
	while ( df >= 0 && ! f[df] ) df--;
	if ( df < 0 ) { strcpy (s, "[0]");  return strlen(s); }
	t = s;
	if ( df >= 2 ) {
		if ( f[df] == 1 ) { t += sprintf (t, "[x^%d", df); } else if ( f[df] == -1 ) { t += sprintf (t, "[-x^%d", df); } else { t += sprintf (t, "[%ld*x^%d", f[df], df); }
	} else if ( df == 1 ) {
		if ( f[df] == 1 ) { t += sprintf (t, "[x"); } else if ( f[df] == -1 ) { t += sprintf (t, "[-x"); } else { t += sprintf (t, "[%ld*x", f[df]); }
	} else {
		t += sprintf (t, "[%ld", f[df]);
	}
	for ( i = df-1 ; i >= 0 ; i-- ) {
		if ( f[i] > 1) {
			t += sprintf (t, " + %ld", f[i]); if ( i ) *t++ = '*';
		} else if ( f[i] == 1 ) {
			t += sprintf (t, " + ");
		} else if ( f[i] == -1 ) {
			t += sprintf (t, " - ");
		} else if ( f[i] < -1 ) {
			t += sprintf (t, " - %ld", -f[i]); if ( i ) *t++ = '*';
		}
		if ( f[i] ) {
			if ( i >= 2 ) {
				t += sprintf (t, "x^%d", i);
			} else if ( i == 1 ) {
				t += sprintf (t, "x");
			} else {
				if ( f[i] == 1 || f[i] == -1 ) *t++ = '1';
			}
		}
	}
	*t++ = ']';
	*t= '\0';
	return t-s;
}
