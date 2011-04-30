
#ifndef FST_INCLUDE
#define FST_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 F S T   R O U T I N E S 
 
 calculates FST
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: fst.h 1323 2008-07-25 19:13:48Z beerli $
-------------------------------------------------------*/
extern void fst_type (int type);
extern void calc_fst (world_fmt * world, data_fmt * data);

#endif /* FST_INCLUDE */
