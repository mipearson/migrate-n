#ifndef _AIC_H_
#define _AIC_H_
#define AICTEST
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 aic model  test   R O U T I N E S 
 
 Peter Beerli 2001, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: aic.h 1323 2008-07-25 19:13:48Z beerli $
-------------------------------------------------------*/
#include "migration.h"

extern void akaike_information (world_fmt * world, long *Gmax);
extern long find_paramnum(world_fmt *world, char *connect);

#endif
