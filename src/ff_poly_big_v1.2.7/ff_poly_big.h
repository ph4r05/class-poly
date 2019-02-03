#ifndef _FF_POLY_BIG_INCLUDE_
#define _FF_POLY_BIG_INCLUDE_

/*
    Copyright 2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// main include file for ff_poly library, exposes the standard set of functionality
// users wanting direct access to internal functions (e.g. in ffpolysmall.h) should include them directly

// FF_POLY_BIG is defined so that fast polynomial multiplication will be used for large degree polys
// This requires linking with the zn_poly library

#define FF_POLY_BIG	1

#include "ff_poly/ff.h"
#include "ff_poly/ff2k.h"
#include "ff_poly/ffext.h"
#include "ff_poly/ffpoly.h"
#include "ff_poly/ffpolybig.h"

#endif
