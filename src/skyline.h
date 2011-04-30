#ifndef _SKYLINEUPDATE_
#define _SKYLINEUPDATE_
/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
    variation over time routines   R O U T I N E S

    Peter Beerli 2006, Tallahassee
    beerli@fsu.edu

    Copyright 2006 Peter Beerli, Tallahassee

    This software is distributed free of charge for non-commercial use
    and is copyrighted. Of course, we do not guarantee that the software
    works and are not responsible for any damage you may cause or have.


 $Id$
  */
#include "migration.h"
extern void calculate_expected_values(tetra **eventbins, long *eventbinnum, MYREAL eventinterval, MYREAL interval, MYREAL age, long from, long to, long * lineages, long numpop, world_fmt *world);
extern void setup_expected_events (world_fmt * world, option_fmt * options);
extern void destroy_expected_events (world_fmt * world);
extern void print_expected_values(world_fmt * world, option_fmt *options);
extern void debug_skyline(world_fmt *world, char text[]);
#endif /*_SKYLINEUPDATE_*/
