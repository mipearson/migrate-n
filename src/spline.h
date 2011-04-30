#ifndef __SPLINEH__
#define __SPLINEH__
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 S P L I N E   R O U T I N E S 
 
 interface part to adaptive spline routines
 from the AMS library
 
 Peter Beerli 1999, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: spline.h 1640 2010-02-06 16:28:21Z beerli $
-------------------------------------------------------*/

extern int dbvssc_ (MYREAL *x, MYREAL *y, long *np, long *n, long *k,
                        long *opt, MYREAL *d0, MYREAL *dnp, MYREAL *d20,
                        MYREAL *d2np, long *constr, MYREAL *eps,
                        MYREAL (*beta) (MYREAL *), MYREAL (*betai) (MYREAL *), MYREAL (*rho) (MYREAL *),
                        MYREAL (*rhoi) (MYREAL *), long *kmax, long *maxstp, long *errc,
                        MYREAL *d, MYREAL *d2, long *diagn, MYREAL *work,
                        long *nwork);

extern int dbvsse_ (MYREAL *x, MYREAL *y, long *np, long *n, long *k,
                        MYREAL *xtab, long *ntab, long *sbopt, long *y0opt,
                        long *y1opt, long *y2opt, long *errc, MYREAL *d,
                        MYREAL *d2, MYREAL *y0tab, MYREAL *y1tab, MYREAL *y2tab,
                        long *erre, MYREAL *work, long *nwork);



#endif
