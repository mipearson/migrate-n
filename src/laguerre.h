#ifndef _LAGUERRE_H_
#define _LAGUERRE_H_
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 Laguerre integration  R O U T I N E S 
 
 Peter Beerli 2000, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: laguerre.h 1323 2008-07-25 19:13:48Z beerli $
-------------------------------------------------------*/

extern void integrate_laguerre (long categs, MYREAL *rate,
                                    MYREAL *probcat,
                                    MYREAL (*func) (MYREAL, helper_fmt *),
                                    helper_fmt * helper, MYREAL *result,
                                    MYREAL *rmax);
extern void initgammacat (long categs, MYREAL alpha, MYREAL theta1,
                              MYREAL *rate, MYREAL *probcat);


#endif
