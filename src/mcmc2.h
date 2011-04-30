

#ifndef MCMC2_INCLUDE
#define MCMC2_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 M C M C 2   R O U T I N E S 
 
 Tree changing routines
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: mcmc2.h 1323 2008-07-25 19:13:48Z beerli $
-------------------------------------------------------*/

extern void coalesce1p (proposal_fmt * proposal);
extern void pretendcoalesce1p (proposal_fmt * proposal);

#endif /*MCMC2_INCLUDE */
