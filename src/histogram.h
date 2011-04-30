// histogram.h
#ifndef _BAYESREREAD_
#define _BAYESREREAD_
//
// started October 2007
// (c) Peter Beerli Tallahassee 2007
// $Id:$
//
#include "migration.h"
#ifdef ZNZ
extern void read_from_bayesmdim_minimal_info(znzFile mdimfile, world_fmt *world,option_fmt *options, data_fmt *data);
extern void read_bayes_fromfile(znzFile mdimfile, world_fmt *world, option_fmt *options);
#else
extern void read_from_bayesmdim_minimal_info(FILE *mdimfile, world_fmt *world,option_fmt *options, data_fmt *data);
extern void read_bayes_fromfile(FILE *mdimfile, world_fmt *world, option_fmt *options);
#endif


#endif
