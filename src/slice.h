#ifndef _SLICESAMPER_
#define _SLICESAMPLER_
//
// started August 1 2006
// Peter Beerli Tallahassee
// $Id$
//
#include "migration.h"
extern MYREAL slice (MYREAL *startval, long which, world_fmt * world, MYREAL  (*func) (MYREAL, MYREAL, bayes_fmt *, long));
extern MYREAL expslice (MYREAL *startval, long which, world_fmt * world, MYREAL  (*func) (MYREAL, MYREAL, bayes_fmt *, long));
extern void expallslice (MYREAL *values, MYREAL *rateval, long which, world_fmt * world);
#endif
