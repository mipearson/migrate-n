#ifndef _EVENTUPDATE_
#define _EVENTUPDATE_
/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
    gathering migration and coalescence events  routines   R O U T I N E S

    Peter Beerli 2006, Tallahassee
    beerli@fsu.edu

    Copyright 2006 Peter Beerli, Tallahassee

    This software is distributed free of charge for non-commercial use
    and is copyrighted. Of course, we do not guarantee that the software
    works and are not responsible for any damage you may cause or have.


$Id$
  */
#include "migration.h"
extern void increase_mighist (mighistloci_fmt * mighistlocus);
extern void setup_mighist (world_fmt * world, option_fmt * options);
extern void destroy_mighist (world_fmt * world);
extern void calculate_event_values(duo **eventbins, long *eventbinnum, 
				   MYREAL eventinterval, MYREAL interval, 
				   MYREAL age, long from, long to, 
				   long * lineages, long numpop, boolean is_last);
extern void setup_event_events (world_fmt * world, option_fmt * options);
extern void destroy_event_events (world_fmt * world);
extern void print_mighist_output (FILE * out, world_fmt * world, MYREAL *sums, boolean mrca);
extern void print_event_values(world_fmt * world);
extern void store_events (world_fmt * world, timelist_fmt * ltl, long np, long rep);

#endif /*_EVENTUPDATE_*/
