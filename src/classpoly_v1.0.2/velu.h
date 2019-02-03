#ifndef _VELU_INCLUDE_
#define _VELU_INCLUDE_

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

#define VELU_MAX_ELL		131072

void velu_isogeny (ff_t *J1, ff_t J0, long t, long ell);				// _ff_p+1-t must be divisible by ell and ell must be prime
int velu_verify_Sylow_cyclic (ff_t J, long t, long ell);				// verifies that elliptic curve with j-invariant J and trace t has cyclic ell-Sylow subgroup

#endif
