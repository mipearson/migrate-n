#ifndef __BROYDENH__
#define __BROYDENH__
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 P A R A M E T E R E S T I M A T I O N   R O U T I N E S 
 
 estimates parameter for each locus
 using a Broyden minimization
 
 
 Peter Beerli 1997, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: broyden.h 1323 2008-07-25 19:13:48Z beerli $
 
-------------------------------------------------------*/

extern void broyden (world_fmt * world, MYREAL **covariance, long numg,
                         long chain, char type, char ***plane);

extern MYREAL absmaxvec (MYREAL *v, long n);
extern void create_nr (nr_fmt * nr, world_fmt * world, long G,
                           long profilenum, long thislocus, long repkind,
                           long rep);
extern void reset_hess (MYREAL **hess, long n);
extern void destroy_nr (nr_fmt * nr, world_fmt * world);
extern void alloc_apg (MYREAL ****apg, long repstop, long loci, long G);

extern MYREAL calc_locus_like (nr_fmt * nr, MYREAL *param, MYREAL *lparam,
                                   long locus);
extern void calc_param (nr_fmt * nr, MYREAL *param, MYREAL lamda);
extern void param_all_adjust (MYREAL *param, nr_fmt *nr);//worldoption_fmt * wopt,
// long numpop);
#ifdef LONGSUM
extern void gradient_longsum (MYREAL *d, nr_fmt * nr, long locus);
#else /*LONGSUM*/
extern void gradient (MYREAL *d, nr_fmt * nr, long locus);
#endif /*LONGSUM*/
extern void grad2loggrad (MYREAL *param, long *indeks, MYREAL *d, long nn,
                              long profilenum);
extern MYREAL probG (MYREAL *param, MYREAL *lparam, tarchive_fmt * tl,
                         nr_fmt * nr, long locus);

extern void calc_cov (MYREAL **dd, MYREAL *d, MYREAL *param, long n);

extern void calc_dv (MYREAL *dv, MYREAL **hess, MYREAL *gxv, long n);
extern MYREAL calc_line (helper_fmt * helper, MYREAL a, MYREAL b, MYREAL c,
                             MYREAL (*psi) (MYREAL, helper_fmt *));
extern void calc_hessian (MYREAL **hess, long n, MYREAL *delta, MYREAL *gama);
extern MYREAL psi (MYREAL lamda, helper_fmt * helper, MYREAL *param,
                       MYREAL *lparam);
extern void log_param0 (MYREAL *param, MYREAL *lparam, long nn);
extern void copies2lcopies (timearchive_fmt * atl);
extern void create_apg0 (MYREAL *apg0, nr_fmt * nr, timearchive_fmt * tyme,
                             long locus);
extern void create_multiapg0 (MYREAL *apg0, nr_fmt * nr, long rep,
                                  long locus);
extern void force_symmetric_d (MYREAL *gxv, long model, nr_fmt * nr, long nn);

extern void print_contribution (nr_fmt * nr, timearchive_fmt ** atl, long G);

extern MYREAL ln_copies (long n);
#endif
