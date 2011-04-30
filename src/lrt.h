#ifndef _LRT_H_
#define _LRT_H_
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 Likelihood ratio test   R O U T I N E S 
 
 moved out from world.c                                                                                                               
 Peter Beerli 2000, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: lrt.h 1323 2008-07-25 19:13:48Z beerli $
-------------------------------------------------------*/

extern void print_lratio_test (world_fmt * world, long *Gmax);
extern void lrt_connect (char *thisString, char *connect2, world_fmt *world);

#endif
