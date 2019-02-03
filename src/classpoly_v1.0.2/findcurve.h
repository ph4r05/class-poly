#ifndef _FINDCURVE_INCLUDE_
#define _FINDCURVE_INCLUDE_

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

#define FINDCURVE_NOFILTER		2
#define FINDCURVE_FILTER_JCUBE	4
#define FINDCURVE_FILTER_WEBER	8

extern unsigned long gop_counter;
double picktor (int *ptor, int *ps2_flag, int *pt3_flag, long p, long t);
int findcurve (ff_t f[4],  long t, int flags, long *cnt);
int scancurves (int *psignm, ff_t x[], ff_t y[], int n, unsigned long pbits0, unsigned long nbits0, unsigned long pbits1, unsigned long nbits1, ff_t f1[]);

#endif
