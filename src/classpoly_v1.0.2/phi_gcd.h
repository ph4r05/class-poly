#ifndef _PHI_GCD_INCLUDE_
#define _PHI_GCD_INCLUDE_

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

void phi_surface_gcd_path (ff_t r1[], ff_t r2[], int p1, int n, int p2, int e, ff_t *phi1, ff_t *phi2);
void phi_surface_gcd_cycle (ff_t r1[], ff_t r2[], int p1, int n, int p2, int e, ff_t *phi1, ff_t *phi2);

#endif
