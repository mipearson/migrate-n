#ifndef MCMC_INCLUDE
#define MCMC_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 M C M C   R O U T I N E S 
 
 Markov Monte Carlo stuff: treechange, acceptance
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: mcmc.h 1323 2008-07-25 19:13:48Z beerli $
-------------------------------------------------------*/

extern long tree_update (world_fmt * world, long g);
extern void free_timevector (timelist_fmt * timevector);
extern void count_migrations (node * p, long *count);
extern void traverse_tagnew(node *theNode, node *origin);
#endif
