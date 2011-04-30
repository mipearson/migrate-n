#ifndef REPORTER_INCLUDE
#define REPORTER_INCLUDE
/** \file reporter.h */
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 
 R E P O R T E R   R O U T I N E S 
                                                                                                               
 Peter Beerli 1999, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: reporter.h 1458 2008-11-28 02:38:18Z beerli $
-------------------------------------------------------*/

#include "migration.h"

extern void convergence_check (world_fmt * world, boolean progress);
extern void calc_chain_s(MYREAL *cs, MYREAL *cm, world_fmt *world, long replicate);
extern void convergence_check_bayes (world_fmt *world, long maxreplicate);
extern void chain_means (MYREAL *thischainmeans, world_fmt * world);
extern void convergence_progress(FILE *file, world_fmt *world);
extern MYREAL single_chain_var(world_fmt *world, long T, MYREAL *variance, MYREAL *autoc, MYREAL *effsample);
extern boolean max_ess(const MYREAL * ess, const long n, const MYREAL minimum);
extern void print_bayes_ess(FILE *file, world_fmt *world, long numparam, int offset, 
				 MYREAL *autocorr, MYREAL *effsample);
extern void calculate_ess_frombayes(world_fmt *world, long T, MYREAL *params, long locus, MYREAL *autoc, MYREAL *effsample);
extern void collect_ess_values(world_fmt *world);
#endif
