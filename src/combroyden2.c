/** \file combroyden.c */
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 M A X I M I Z E R      R O U T I N E S
 
 - calculates single locus ML
 - multilocus ML with constant mutation rate
 - multilocus ML with gamma deviated mutation rate 
 
  using Broyden-Fletcher-Goldfarb-Shanno method
 
 
 Peter Beerli 1996-1998, Seattle
 beerli@fsu.edu
 
 
    Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
    This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
$Id: combroyden2.c 1854 2011-03-24 18:07:30Z beerli $
 
-------------------------------------------------------*/

#include "migration.h"
#include "sighandler.h"
#include "broyden.h"
#include "joint-chains.h"
#include "world.h"
#include "integrate.h"
#include "migrate_mpi.h"
#include "reporter.h"
#include "tools.h"
#include "laguerre.h"
#include "gammalike.h"
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif


//if using slownet then this should be activated
// otherwise I will use a define in defintions.h that that
// replaces all occurrences with CALCLIKE with calc_loci_like
// for speed reasons we don't want any further indirection
// for which we have no use.
#ifdef SLOWNET
MYREAL (*calc_like) (helper_fmt *, MYREAL *, MYREAL *);
void (*calc_gradient) (nr_fmt *, helper_fmt *, MYREAL *);
void (*setup_param0) (world_fmt *, nr_fmt *, long, long, long, long, long,
                      boolean);
#endif

// external definitions -------------------------------
extern void debug_plot (helper_fmt * helper);
extern int myID;

/* prototypes ------------------------------------------- */
/* public function used in main.c */
long estimateParameter (long rep, long G, world_fmt * world,
                        option_fmt * options, MYREAL **dd, long chain,
                        char type, long kind, long repkind);

/* calculates likelihood */
MYREAL calc_loci_like (helper_fmt * helper, MYREAL *param, MYREAL *lparam);

/* first derivatives for ML with constant mutation rate */
void simple_loci_derivatives (MYREAL *d, nr_fmt * nr,
                              timearchive_fmt ** tyme, long locus);
/* progress-reporter */
void print_menu_finalestimate (boolean progress, long locus, long loci);

/* start parameter for multilocus estimator, creates simple average */
void calc_meanparam (world_fmt * world, long n, long repstart, long repstop);

void expected_param(MYREAL *param, world_fmt *world);

/* multilocus first derivatives with or without gamma distrib. mut. rate */
void combine_gradient (nr_fmt * nr, helper_fmt * helper, MYREAL *gxv);
#ifdef SLOWNET
void slow_net_combine_gradient (nr_fmt * nr, helper_fmt * helper,
                                MYREAL *gxv);
#endif
/* this helps to use the simple single locus material */
void copy_and_clear_d (nr_fmt * nr);
void add_back_d (nr_fmt * nr);
/* rescales parameter with theta1 for gamma calculation */
void set_gamma_param (MYREAL *paramn, MYREAL *paramo,
                      MYREAL *lparamn, MYREAL *lparamo,
                      MYREAL theta, nr_fmt * nr);
/* calculates norm, this creates a stopping criteria */
MYREAL norm (MYREAL *d, long size);
/* used in multilocus  line search */
MYREAL psilo (MYREAL lamda, helper_fmt * helper, MYREAL *param,
              MYREAL *lparam);
/* used in integral for gamma distr. mutation rate */
MYREAL phi (MYREAL t, helper_fmt * b);

/* calculates new multilocus parameter values */
void calc_loci_param (nr_fmt * nr, MYREAL *lparam, MYREAL *olparam,
                      MYREAL *dv, MYREAL lamda, long nnn);
/* used to reset the BFSG stuff */
void reset_work (MYREAL **work, long nnn);

/* linesearch of BFGS
   and helper functions for this linesearch */
MYREAL line_search (MYREAL *dv, MYREAL *gxv,
                    MYREAL (*fpsi) (MYREAL lamda, helper_fmt * helper,
                                    MYREAL *param, MYREAL *lparam),
                    helper_fmt * thishelper, MYREAL oldL, MYREAL *tempparam,
                    MYREAL *templparam);

void replace_with (int mode, MYREAL *a, MYREAL *b, MYREAL *c, MYREAL m,
                   MYREAL *la, MYREAL *lb, MYREAL *lc, MYREAL ll);
MYREAL psil (MYREAL (*fpsi) (MYREAL lamda, helper_fmt * helper), MYREAL a,
             MYREAL b, helper_fmt * helper, MYREAL *gxv, MYREAL *dv,
             MYREAL oldL, MYREAL *val1, MYREAL *val2);
MYREAL quadratic_lamda (MYREAL lamda, MYREAL a, MYREAL b, MYREAL c, MYREAL la,
                        MYREAL lb, MYREAL lc);


void print_reset_big_param (FILE * file, MYREAL average, MYREAL maxtheta,
                            long pop);

long do_profiles (world_fmt * world, nr_fmt * nr, MYREAL *likes,
                  MYREAL *normd, long kind, long rep, long repkind);

void correct_values (world_fmt * world, nr_fmt * nr);

void which_calc_like (long repkind);

void set_replicates (world_fmt * world, long repkind, long rep,
                     long *repstart, long *repstop);

void prepare_broyden (long kind, world_fmt * world, boolean * multilocus);

long multilocus_gamma_adjust (boolean multilocus, boolean boolgamma,
                              world_fmt * world, long repstart, long repstop);
void profile_setup (long profilenum, long *profiles, MYREAL *param,
                    MYREAL *value, long *nnn);

void setup_parameter0_standard (world_fmt * world, nr_fmt * nr, long repkind,
                                long repstart, long repstop, long loci,
                                long kind, boolean multilocus);

void report_broyden (MYREAL newL, MYREAL normd, long trials,
                     boolean boolgamma, MYREAL theta1, MYREAL alpha,
                     FILE * logfile);
boolean guard_broyden (MYREAL newL, MYREAL oldL, long trials, MYREAL normd,
                       MYREAL *normd20);

void fill_helper (helper_fmt * helper, MYREAL *param, MYREAL *lparam,
                  world_fmt * world, nr_fmt * nr);

/* maximizer routine */
long maximize (MYREAL **thisparam, world_fmt * world, nr_fmt * nr, 
	       MYREAL **hess,
	       long analystype,
               long repkind);

void set_both_delta (nr_fmt * nr,
                     MYREAL *delta, MYREAL *den, MYREAL *deo,
                     MYREAL *gama, MYREAL *gxv, MYREAL *gxv0, long size);
void set_expxv (MYREAL *expa, MYREAL *a, long size);
void set_expparam (MYREAL *expa, MYREAL *a, long size);
void set_logparam (MYREAL *loga, MYREAL *a, long size);
void finish_broyden (MYREAL newL, MYREAL normd, long trials,
                     world_fmt * world, nr_fmt * nr, MYREAL *expxv, long nn,
                     long analystype, long repkind);
void fill_profile_index (nr_fmt * nr);


/* Functions ++++++++++++++++++++++++++++++++++++++++++++++++*/
/* public functions---------------------------------------------- */
/* driver for maximizing single locus, multilocus functions      */
long
estimateParameter (long rep, long G, world_fmt * world, option_fmt * options,
                   MYREAL **dd, long chain, char type, long kind,
                   long repkind)
//timearchive_fmt * tyme
{
    long trials = 0;
    MYREAL normd = 0.;
    MYREAL likes = 0.;
    MYREAL *param;
#ifndef MPI
    long repstart, repstop;
#endif

    boolean multilocus;
    nr_fmt *nr;
    nr = (nr_fmt *) mycalloc (1, sizeof (nr_fmt));
    if (kind == MULTILOCUS)
    {
        world->locus = world->loci;
        print_menu_finalestimate (world->options->progress, world->loci,
                                  world->loci);
    }

    create_nr (nr, world, G, 0, world->locus, repkind, rep);
    set_replicates (world, world->repkind, world->rep, &nr->repstart,
                    &nr->repstop);
    //    fprintf(stdout,"%i> DEBUG sumfile: world-rep=%li world->repkind=%li nr->repstart=%li nr->repstop=%li\n",myID, world->rep, world->repkind, nr->repstart, nr->repstop);
    prepare_broyden (kind, world, &multilocus);
#ifndef INTEGRATEDLIKE
    //    if(world->options->integrated_like)
      SETUPPARAM0 (world, nr, world->repkind, nr->repstart, nr->repstop,
                 world->loci, kind, multilocus);
#else
      SETUPPARAM0 (world, nr, world->repkind, nr->repstart, nr->repstop,
                 world->loci, kind, multilocus);
#endif
    doublevec1d (&param, nr->partsize + 1);

#ifdef MPI
    if (myID == MASTER)
        memcpy (param, world->param0, sizeof (MYREAL) * nr->numpop2);
#else
    if (kind == MULTILOCUS || repkind != SINGLECHAIN)
    {
        set_replicates (world, world->repkind, world->rep, &repstart, &repstop);
        calc_meanparam (world, world->numpop2, repstart, repstop);
        memcpy (param, world->param0, sizeof (MYREAL) * nr->numpop2);
    }
  //    {
    //    expected_param(param,world);
     // }
#endif
    memcpy (param, world->param0, sizeof (MYREAL) * nr->numpop2);

#ifdef LONGSUM

    memcpy (param+nr->numpop2+
            ((world->options->gamma && kind==MULTILOCUS) ? 1 : 0) , world->flucrates, sizeof (MYREAL) * nr->numpop*3);
#endif

    if (strpbrk (world->options->custm2, "0c"))
    {
        if (kind == SINGLELOCUS
                && world->options->custm2[world->numpop2] == 'c')
            trials = maximize (&param, world, nr, dd, kind, repkind);
        else 
            trials = do_profiles (world, nr, &likes, &normd, kind, rep, repkind);
    }
    else
      {
        trials = maximize (&param, world, nr, dd, kind, repkind);
      }
    if (kind == SINGLELOCUS && world->options->progress)
    {
        print_menu_chain (type, chain, world->atl[rep][world->locus].T,
                          world, options, rep);
        FPRINTF (stdout, "           Maximization cycles needed: %li\n           Norm of first derivatives:  %f\n",
                 trials, nr->normd);
        if (world->options->logfile)
        {
            FPRINTF (world->options->logfile,
                     "           Maximization cycles needed: %li\n           Norm of first derivatives:  %f\n", trials, nr->normd);
        }
    }
#ifdef MPI
    if (myID==MASTER  && world->options->plotnow)
#else
      if (kind == MULTILOCUS && world->options->plotnow && world->loci > 1)
#endif
    {
      if (((world->options->replicate && repkind != SINGLECHAIN) ||
                (!world->options->replicate && repkind == SINGLECHAIN)) && (world->numpop > 1))
#ifdef MPI
	{
	  if(world->loci==1)
	    create_plot (world, world->plane[0], nr, world->loci, MULTILOCUSPLOT); 
	  else
	    create_plot (world, world->plane[world->loci], nr, world->loci, MULTILOCUSPLOT); 
	}
#else
      create_plot (world, world->plane[world->loci], nr, world->loci, MULTILOCUSPLOT); 
#endif
    }
#ifdef LONGSUM
    memcpy (world->flucrates, param+nr->numpop2+
            ((world->options->gamma && kind==MULTILOCUS) ? 1 : 0), sizeof (MYREAL) * nr->numpop*3);
#endif

    world->trials = trials;
    world->normd = nr->normd;
    myfree(param);
    destroy_nr (nr, world);
    return trials;
}


/* calculates the likelihood over all loci for the new parameter set */
MYREAL
calc_loci_like (helper_fmt * helper, MYREAL *param, MYREAL *lparam)
{
#ifndef MPI
    long locus;
    MYREAL *mu_rates = helper->nr->world->options->mu_rates;
#endif

    nr_fmt *nr = helper->nr;
    world_fmt *world = helper->nr->world;
    MYREAL logres = 0;
    if (helper->multilocus)
    {
#ifdef MPI
        /* must be MASTER node!
           mutation rate is constant log(L(loci)) = Sum[log(L(locus))) */
        //      printf("%i> before mpi_likelihood_master\n",myID);
        logres = mpi_likelihood_master (param, lparam, world, nr,
                                        helper, world->who);
#else

        if (!helper->boolgamma)
        {
            for (locus = 0; locus < world->loci; locus++)
            {
                if (nr->skiploci[locus])
                    continue;
                nr->locilikes[locus] =
                    calc_locus_like (nr, param, lparam, locus) + mu_rates[locus];
                logres += nr->locilikes[locus];
            }
        }
        else
        {
            logres = gamma_loci_like (helper, param, lparam, helper->weight);
        }
#endif /* NOT MPI */

    }
    else    /* ->singlelocus */
        logres = calc_locus_like (nr, param, lparam, nr->world->locus);
    nr->llike = logres;
    return nr->llike;
}

#ifdef SLOWNET
/* calculates the likelihood over all loci for the new parameter set */
MYREAL
slow_net_calc_loci_like (helper_fmt * helper, MYREAL *param, MYREAL *lparam)
{
    long locus;
    nr_fmt *nr = helper->nr;
    world_fmt *world = helper->nr->world;
    MYREAL *mu_rates = world->options->mu_rates;
    MYREAL logres = 0;
    if (helper->multilocus)
    {
        if (!helper->boolgamma)
        {
            for (locus = 0; locus < world->loci; locus++)
            {
                if (nr->skiploci[locus])
                    continue;
                nr->locilikes[locus] =
                    calc_locus_like (nr, param, lparam, locus) + mu_rates[locus];
                logres += nr->locilikes[locus];
            }
        }
        else
        {
            logres = gamma_loci_like (helper, param, lparam, helper->weight);
        }
    }
    else    /* ->singlelocus */
        logres = calc_locus_like (nr, param, lparam, nr->world->locus);
    nr->llike = logres;
    return nr->llike;
}
#endif /*SLOWNET*/
/*
MYREAL
phi (MYREAL t, helper_fmt *b)
{
  MYREAL ll, alpha, beta;
  long locus;
  MYREAL weight;
  helper_fmt *helper;
  nr_fmt *nr;
  timearchive_fmt **atl;
  helper = (helper_fmt *) b;
  locus = (long) helper->locus;
  weight = (MYREAL) helper->weight;
  atl = (timearchive_fmt **) helper->atl;
  nr = (nr_fmt *) helper->nr;
 
  alpha = nr->oparam[nr->numpop2];
  beta = nr->oparam[0] / alpha;
  set_gamma_param (nr->param, nr->oparam, nr->lparam, nr->olparam, t, nr);
  ll = calc_locus_like (nr, nr->param, nr->lparam, locus);
//#ifdef LAGUERRE
  ll = -t / beta + (alpha - 1.) * LOG (t) + ll  + weight;
//#else
//  ll = EXP (-t / beta + (alpha - 1.) * LOG (t) + ll  + weight);
//#endif
  return ll;
}
 
*/
/* private functions---------------------------------------------- */
/* derivatives */
void
combine_gradient (nr_fmt * nr, helper_fmt * helper, MYREAL *gxv)
{
#ifndef MPI
    long locus;
#endif

    timearchive_fmt **tyme = NULL;
    //    MYREAL save_a, denom, fla, fha, la, ha;
    long nnn = nr->partsize - nr->profilenum;
    memset (gxv, 0, sizeof (MYREAL) * nnn);
    if (helper->analystype == SINGLELOCUS)
    {
        simple_loci_derivatives (gxv, nr, tyme, nr->world->locus);
    }
    else
    {
#ifdef MPI
        if (myID == MASTER)
        {
            mpi_gradient_master (nr, nr->world, nr->world->who);
            memcpy (gxv, nr->d, sizeof (MYREAL) * nr->partsize);
        }
        else   //SLAVE
        {
            error ("not allowed to execute mpi_gradient_worker()\n");
        }
#else
        if (!helper->boolgamma)
        {
            for (locus = 0; locus < nr->world->loci; locus++)
            {
                if (nr->skiploci[locus])
                    continue;
                memset (nr->d, 0, sizeof (MYREAL) * nnn);
                simple_loci_derivatives (nr->d, nr, tyme, locus);
                add_vector (gxv, nr->d, nnn);
            }
        }
        else
        {
            //#ifdef LAGUERRE
            //if(nr->normd<5.)
            //   gamma_loci_difference(helper);
            //else
            gamma_loci_derivative (helper);
            memcpy (gxv, nr->d, sizeof (MYREAL) * nr->partsize);
            /*   save_a = helper->xv[nr->partsize-1];
            la= save_a - save_a/100.;
            helper->xv[nr->partsize-1]=la;
            helper->expxv[nr->partsize-1]=EXP(la);
            denom = -LGAMMA(helper->expxv[nr->partsize-1]) - (helper->xv[0]-helper->xv[nr->partsize-1]) * helper->xv[nr->partsize-1];
            fla = gamma_loci_like (helper, helper->expxv, helper->xv, denom);
            ha= save_a + save_a/100.;
            helper->xv[nr->partsize-1]=ha;
            helper->expxv[nr->partsize-1]=EXP(ha);
            denom = -LGAMMA(helper->expxv[nr->partsize-1]) - (helper->xv[0]-helper->xv[nr->partsize-1]) * helper->xv[nr->partsize-1];
            fha = gamma_loci_like (helper, helper->expxv, helper->xv, denom);
            nr->d[nr->partsize-1] = (fla - fha)/(la-ha);
            memcpy (gxv, nr->d, sizeof (MYREAL) * nr->partsize);*/
            //#else
            //   dt (gxv, helper);
            //endif  ****LAGUERRE****
        }
#endif /*MPI*/

    }
    force_symmetric_d (gxv, nr->world->options->migration_model, nr, nnn);
    grad2loggrad (helper->expxv, nr->indeks, gxv, nnn, nr->profilenum);
}

#ifdef SLOWNET
void
slow_net_combine_gradient (nr_fmt * nr, helper_fmt * helper, MYREAL *gxv)
{
    long locus;
    timearchive_fmt **tyme = NULL;

    long nnn = nr->partsize - nr->profilenum;
    memset (gxv, 0, sizeof (MYREAL) * nnn);
    if (helper->analystype == SINGLELOCUS)
    {
        simple_loci_derivatives (gxv, nr, tyme, nr->world->locus);
        force_symmetric_d (gxv, nr->world->options->migration_model, nr, nnn);
        grad2loggrad (helper->expxv, nr->indeks, gxv, nnn, nr->profilenum);
    }
    else
    {
        if (!helper->boolgamma)
        {
            for (locus = 0; locus < nr->world->loci; locus++)
            {
                if (nr->skiploci[locus])
                    continue;
                memset (nr->d, 0, sizeof (MYREAL) * nnn);
                simple_loci_derivatives (nr->d, nr, tyme, locus);
                add_vector (gxv, nr->d, nnn);
            }
            force_symmetric_d (gxv, nr->world->options->migration_model, nr,
                               nnn);
            grad2loggrad (helper->expxv, nr->indeks, gxv, nnn, nr->profilenum);
        }
        else
        {
	  //#ifdef LAGUERRE
            gamma_loci_derivative (helper);
            memcpy (gxv, nr->d, sizeof (MYREAL) * nr->partsize);
	    //#else
	    //
            //dt (gxv, helper);
	    //#endif

            force_symmetric_d (gxv, nr->world->options->migration_model,
                               nr, nnn);
            grad2loggrad (helper->expxv, nr->indeks, gxv, nnn, nr->profilenum);
        }
    }
}
#endif

void
simple_loci_derivatives (MYREAL *d, nr_fmt * nr, timearchive_fmt ** tyme,
                         long locus)
{
    long g;
#ifndef LONGSUM

    gradient(d, nr, locus);
#else /*LONGSUM*/

    gradient_longsum (d, nr, locus);
#endif /*LONGSUM*/

    for (g = 0; g < nr->partsize - nr->profilenum; g++)
        d[g] /= nr->PGC[locus];
}

void
set_gamma_param (MYREAL *param, MYREAL *oparam,
                 MYREAL *lparam, MYREAL *olparam, MYREAL theta, nr_fmt * nr)
{
    long i;
    MYREAL logtheta;
    param[0] = theta;
    lparam[0] = logtheta = LOG (theta);
    for (i = 1; i < nr->numpop; i++)
    {
        lparam[i] = olparam[i] + logtheta - olparam[0];
        param[i] = EXP (lparam[i]);
    }
    for (i = nr->numpop; i < nr->numpop2; i++)
    {
        lparam[i] = olparam[i] + olparam[0] - logtheta;
        param[i] = EXP (lparam[i]);
    }
    lparam[i] = olparam[i];
    param[i] = EXP (olparam[i]);
}


void
copy_and_clear_d (nr_fmt * nr)
{
    memcpy (nr->od, nr->d, sizeof (MYREAL) * nr->partsize);
    memset (nr->d, 0, sizeof (MYREAL) * nr->partsize);
}

void
add_back_d (nr_fmt * nr)
{
    long pop;
    for (pop = 0; pop < nr->partsize; pop++)
    {
        nr->d[pop] += nr->od[pop];
    }
}


void
print_menu_finalestimate (boolean progress, long locus, long loci)
{
    static char nowstr[STRSIZE];
    if (progress)
    {
        get_time (nowstr, "%H:%M:%S");
        if (locus < loci)
            FPRINTF (stdout, "%s   Parameter estimation for locus %li\n", nowstr,
                     locus);
        else
            FPRINTF (stdout, "%s   Parameter estimation over all loci\n", nowstr);
    }
}

void
print_reset_big_param (FILE * file, MYREAL average, MYREAL maxtheta, long pop)
{
    static boolean checkparam[MAXPOP];
    static boolean done = FALSE;
    static char nowstr[STRSIZE];
    if (!done)
    {
        done = TRUE;
        memset (&checkparam, 0, sizeof (boolean) * MAXPOP);
    }
    if (!checkparam[pop])
    {
        get_time (nowstr, "%H:%M:%S");
        FPRINTF (file,
                 "%s   Theta for population %li was reset to the average\n",
                 nowstr, pop);
        FPRINTF (file, "          of %f, it exceeded the maximum of %f \n",
                 average, maxtheta);
    }
}

void
calc_meanparam (world_fmt * world, long n, long repstart, long repstop)
{
    long r, pop, locus, loci = world->loci;


    memset (world->param0, 0, sizeof (MYREAL) * n);
    for (r = repstart; r < repstop; r++)
    {
        for (locus = 0; locus < loci; locus++)
        {
            for (pop = 0; pop < n; pop++)
            {
                if (!world->data->skiploci[locus])
                    world->param0[pop] += world->atl[r][locus].param[pop];
            }
        }
    }
    for (pop = 0; pop < n; pop++)
    {
        world->param0[pop] /= (loci - world->skipped) * (repstop - repstart);
    }
}


MYREAL
norm (MYREAL *d, long size)
{
    int i;
    MYREAL summ = 0.;
    for (i = 0; i < (int) size; i++)
    {
        summ += d[i] * d[i];
    }
    return sqrt (summ);
}


MYREAL
psilo (MYREAL lamda, helper_fmt * helper, MYREAL *param, MYREAL *lparam)
{
    MYREAL like;
    nr_fmt *nr = helper->nr;
    calc_loci_param (nr, lparam, helper->xv, helper->dv, lamda, nr->partsize);
    set_expparam (param, lparam, nr->partsize);
    fill_helper (helper, param, lparam, nr->world, nr);
    like = CALCLIKE (helper, param, lparam);
    //printf("%f ",like);
    return -like;
}

void
calc_loci_param (nr_fmt * nr, MYREAL *lparam, MYREAL *olparam, MYREAL *dv,
                 MYREAL lamda, long nnn)
{
    long i, ii, z = 0;
    if (nr->profilenum == 0)
    {
        for (i = 0; i < nnn; i++)
            lparam[i] = (MAX (-30., MIN (olparam[i] - lamda * dv[i], 30.)));
    }
    else
    {
        for (i = 0; i < nr->partsize - nr->profilenum; i++)
        {
            ii = nr->indeks[z++];
            lparam[ii] = (MAX (-30., MIN (olparam[ii] - lamda * dv[i], 30.)));
        }
    }
    param_all_adjust (lparam, nr);
    //#ifdef LAGUERRE

    if (nr->world->locus >= nr->world->loci && nr->world->options->gamma)
    {
        if (lparam[nr->numpop2] > 9.903487553)
        {
            lparam[nr->numpop2] = 9.903487553;
        }
        initgammacat (nr->categs, EXP (lparam[nr->numpop2]),1.0 /*EXP (lparam[0])*/,
                      nr->rate, nr->probcat);
        //}
    }
    //#endif

}


void
reset_work (MYREAL **work, long nnn)
{
    long i, j;
    for (i = nnn; i < nnn + 5; i++)
    {
        for (j = 0; j < nnn + 1; j++)
            work[i][j] = 0.0;
    }
}


MYREAL
calc_d_g (MYREAL *d, MYREAL *g, long n)
{
    long i;
    MYREAL sum = 0.;
    for (i = 0; i < n; i++)
        sum += d[i] * g[i];
    return sum;
}


#define STARTLAMBDA 1
#define M_      0

/*
 returns LAMBDA (how far to jump in the search direction)
   needed in the the multilocus broyden-fletcher-shanno maximizer
 
   DV    search direction
   GXV   first derivatives of function to minimize
   FUNC  the function to minimize
   NR    all additional variables needed to calculate FUNC
   OLDL  -log(likelihood) of FUNC with LAMBDA=0 (previously calculated)
 
   needs functions: psil, replace_with, quadratic_lamda
 */
MYREAL
line_search (MYREAL *dv, MYREAL *gxv,
             MYREAL (*fpsi) (MYREAL lamda, helper_fmt * helper, MYREAL *param,
                             MYREAL *lparam), helper_fmt * thishelper,
             MYREAL oldL, MYREAL *tempparam, MYREAL *templparam)
{
  //#ifdef MYREAL == float
  //const MYREAL eps = FLT_EPSILON ;
  //#else
  const MYREAL eps = DBL_EPSILON ;
  //#endif

    long trials = 0;
    MYREAL locallamda = STARTLAMBDA, ql;
    MYREAL a, b, c;
    MYREAL psiabc;
    MYREAL ll, la, lb, lc;
    MYREAL m = 1.0;
    MYREAL denom;
    MYREAL nom;
    helper_fmt helper;
    helper = *thishelper;
    //  partcopy_helper(thishelper, &helper);
    a = 0.0;
    b = 1.0;

    la = -oldL;
    lb = -oldL;   //(*fpsi) (b, &helper, tempparam, templparam);
    //  ll = (*fpsi) (m, &helper, tempparam, templparam);
    //if(ll > la)
    //  c =  -1.;
    //else
    //  c = -2.;
    c = calc_d_g (dv, gxv, helper.nr->partsize - helper.nr->profilenum);
    lc = (*fpsi) (c, &helper, tempparam, templparam);
    while (trials++ < NTRIALS)
    {
        //            printf ("*%3li> %f %f %f / %f %f %f\n", trials, a, b, c, la, lb, lc);
      denom = (2. * (b * la - c * la - a * lb + c * lb +  a * lc - b * lc));
      if((fabs(denom)-eps)>eps)
	{
	  nom = (c * c * (lb - la) + b * b * (la - lc) +
		 a * a * (lc - lb));
	}
      else
	return 1.; 
      locallamda = nom / denom;
      psiabc = ((la - lb) / (-a + b) + (la - lc) / (a - c)) / (-b + c);
      if ((psiabc <= 0.0) || (locallamda >= m))
        {
            if ((a == m || b == m || c == m))
                return m;
            else
            {
                ll = (*fpsi) (m, &helper, tempparam, templparam);
                //      replace_with(M_,&a,&b,&c,m,&la,&lb,&lc,ll);
                replace_with ((int) m, &a, &b, &c, m, &la, &lb, &lc, ll);
                continue;
            }
        }
        ll = (*fpsi) (locallamda, &helper, tempparam, templparam);
        ql = quadratic_lamda (locallamda, a, b, c, la, lb, lc);
        if ((fabs (ll - MYMIN3 (la, lb, lc)) <= BIGEPSILON)
                || (fabs (ll - ql) <= BIGEPSILON))
            return locallamda;
        else
        {
            /*  if(((a<b<c) || (c<b<a)) && (lb < MIN(la,lc)))
               return locallamda;
               if(((b<a<c) || (c<a<b)) && (la < MIN(lb,lc)))
               return locallamda;
               if(((a<c<b) || (b<c<a)) && (lc < MIN(la,lb)))
               return locallamda; */
            replace_with ((int) STARTLAMBDA, &a, &b, &c, locallamda, &la, &lb, &lc,
                          ll);
            m = MYMAX3 (a, b, c);
        }
    }

    return locallamda;
}

void
replace_with (int mode, MYREAL *a, MYREAL *b, MYREAL *c, MYREAL m, MYREAL *la,
              MYREAL *lb, MYREAL *lc, MYREAL ll)
{
    MYREAL ma, mb, mc;
    if (mode == STARTLAMBDA)
    {
        ma = *la;
        mb = *lb;
        mc = *lc;
    }
    else
    {
        ma = ll - *la;
        mb = ll - *lb;
        mc = ll - *lc;
    }
    if (ma > mb)
    {
        if (ma > mc)
        {
            *a = m;
            *la = ll;
        }
        else
        {
            *c = m;
            *lc = ll;
        }
    }
    else
    {
        if (mb > mc)
        {
            *b = m;
            *lb = ll;
        }
        else
        {
            *c = m;
            *lc = ll;
        }
    }
}

/* not used
MYREAL
psil (MYREAL (*fpsi) (MYREAL lamda, helper_fmt * helper), MYREAL a, MYREAL b,
      helper_fmt * helper, MYREAL *gxv, MYREAL *dv, MYREAL oldL, MYREAL *val1,
      MYREAL *val2)
{
  long i;
  long nn = helper->nr->partsize;
  MYREAL value = 0.0;
 
  if (a == 0.0 && b == 0)
    {
      for (i = 0; i < nn; i++)
 {
   value += dv[i];  // gxv[i]  ;
 } 
      *val1 = -oldL;
      *val2 = -oldL;
      return value;
    }
  else
    {
      if (a == 0)
 value = (((*val2) = (*fpsi) (b, helper)) - ((*val1) = -oldL)) / b;
      else
 {
   if (b == 0)
     value =
       (((*val2) = -oldL) - ((*val1) = (*fpsi) (a, helper))) / (-a);
   else
     value =
       (((*val2) = (*fpsi) (b, helper)) - ((*val1) =
        (*fpsi) (a,
          helper))) / (b -
         a);
 }
      return value;
    }
}
*/

MYREAL
quadratic_lamda (MYREAL lamda, MYREAL a, MYREAL b, MYREAL c, MYREAL la,
                 MYREAL lb, MYREAL lc)
{
    MYREAL alpha, beta, gama;
    MYREAL aa, bb, cc;
    aa = a * a;
    bb = b * b;
    cc = c * c;
    alpha =
        ((c - b) * la + (a - c) * lb +
         (b - a) * lc) / ((b - a) * (c - a) * (c - b));
    beta =
        ((cc - bb) * la + (aa - cc) * lb +
         (bb - aa) * lc) / ((b - a) * (b - c) * (c - a));
    gama =
        (b * cc * la - bb * c * la + aa * c * lb - a * cc * lb - aa * b * lc +
         a * bb * lc) / ((b - a) * (c - a) * (c - b));


    return alpha * lamda * lamda + beta * lamda + gama;

}

#define TAU1  2.6180339887
#define TAU   1.6180339887
#define BSTART  -2.360679775
MYREAL golden_section(MYREAL (*fpsi) (MYREAL lamda, helper_fmt * helper, MYREAL *param,
    MYREAL *lparam), helper_fmt * thishelper,
    MYREAL *tempparam, MYREAL *templparam)
{
    MYREAL a = -10.;
    MYREAL b = BSTART;
    MYREAL c = 10.;
    MYREAL bl;
    MYREAL d=1.;
    MYREAL dl = -MYREAL_MAX;
    MYREAL ll = 0.;
    helper_fmt helper;
    helper = *thishelper;
    bl = (*fpsi) (b, &helper, tempparam, templparam);
    while(fabs(ll - dl)>EPSILON)
      {
        ll = dl;
        if((c-b) > (b-a))
          {
            d = b + (c - b)/ TAU1;
            dl = (*fpsi) (d, &helper, tempparam, templparam);
            if(dl > bl)
                c  = d;
            else
              {
                a  = b;
                b  = d;
                bl = dl;
              }
          }            
        else
          {
            d = a + (b - a)/ TAU1;
            dl = (*fpsi) (d, &helper, tempparam, templparam);
            if(dl > bl)
                a = d;
            else
              {
                c = b;
                b = d;
                bl = dl;
              }
          }        
      }
    return d;
}

long
do_profiles (world_fmt * world, nr_fmt * nr, MYREAL *likes,
             MYREAL *normd, long kind, long rep, long repkind)
{
  MYREAL **hess;
    char *p;
    long z = 0;
    long trials;
    long i = 0;
    p = world->options->custm2;
    while (*p != '\0')
    {
        if (*p == '0' || *p == 'c')
        {
            nr->profiles[z] = i;
            nr->values[z] = world->param0[i];
            z++;
        }
        p++;
        i++;
    }
    nr->profilenum = z;
    if(kind!=SINGLELOCUS)
      world->locus = world->loci;
    // use a private copy of the hessian so not to clobber the
    // version we use to print out the -fisher information
    z = world->numpop2 + 
      (world->locus == world->loci ? 
       world->options->gamma : 0);
      doublevec2d(&hess,z,z);
    trials = maximize (&world->param0, world, nr, hess, kind, repkind);
    myfree(hess[0]);
    myfree(hess);
    return trials;
}

void
correct_values (world_fmt * world, nr_fmt * nr)
{
    long pop, i, ii, z;
    MYREAL ssum;
    long nsum;
    long elem = (world->options->gamma ? nr->numpop2 : nr->partsize);
    for (pop = 0; pop < world->numpop; pop++)
    {
        if (world->param0[pop] >= BIGGEST_THETA)
        {
            ssum = 0;
            nsum = 0;
            for (i = 0; i < world->numpop; i++)
            {
                if (world->param0[i] <= BIGGEST_THETA)
                {
                    ssum += world->param0[i];
                    nsum++;
                }
            }
            if (nsum > 0 && ssum > 0.0)
                world->param0[pop] = ssum / (MYREAL) nsum;
            else
            {
	      if (!(ssum > 0.0))
                    world->param0[pop] = SMALLEST_THETA;
                else
                    world->param0[pop] = BIGGEST_THETA;
            }
            if (world->options->writelog)
                print_reset_big_param (world->options->logfile,
                                       world->param0[pop], BIGGEST_THETA, pop);
            //        print_reset_big_param (stdout, world->param0[pop], BIGGEST_THETA,
            //                               pop);
        }
    }
    z = 0;
    for (i = 0; i < elem - nr->profilenum; i++)
    {
        ii = (nr->profilenum > 0) ? nr->indeks[z++] : i;
        if (world->param0[ii] >= BIGGEST_MIGRATION)
            world->param0[ii] = BIGGEST_MIGRATION;
        if (world->param0[ii] < SMALLEST_MIGRATION)
            world->param0[ii] = SMALLEST_MIGRATION;
    }
}

/// decide which conditional likelihood calculator to use
void
which_calc_like (long kind)
{
    switch (kind)
    {
    case PROFILE:
#ifdef SLOWNET

        calc_like =
            (MYREAL (*)(helper_fmt *, MYREAL *, MYREAL *))
            slow_net_calc_loci_like;
        calc_gradient =
            (void (*)(nr_fmt *, helper_fmt *, MYREAL *))
            slow_net_combine_gradient;
        setup_param0 =
            (void (*)
             (world_fmt *, nr_fmt *, long, long, long, long, long,
              boolean)) setup_parameter0_slowmpi;
#endif

        break;
    default:
#ifdef SLOWNET

        calc_like =
            (MYREAL (*)(helper_fmt *, MYREAL *, MYREAL *)) calc_loci_like;
        calc_gradient =
            (void (*)(nr_fmt *, helper_fmt *, MYREAL *)) combine_gradient;
        setup_param0 =
            (void (*)
             (world_fmt *, nr_fmt *, long, long, long, long, long,
              boolean)) setup_parameter0_mpi;
#endif

        break;
    }
}

/// calculates the start and end of a replication step, sets repstop = repstart+1
void
set_replicates (world_fmt * world, long repkind, long rep, long *repstart,
                long *repstop)
{
    if (world->options->replicate)
    {
        *repstart = (repkind == SINGLECHAIN) ? rep : 0;
        *repstop = (repkind == SINGLECHAIN) ? rep + 1 : world->repstop;
    }
    else
    {
        *repstart = 0;
        *repstop = 1;
    }
}



void
prepare_broyden (long kind, world_fmt * world, boolean * multilocus)
{
    if (kind == SINGLELOCUS || (kind == PROFILE && world->loci == 1))
    {
        if (world->loci == 1)
            world->locus = 0;
        *multilocus = FALSE;
    }
    else
    {
        *multilocus = TRUE;
    }
}


long
multilocus_gamma_adjust (boolean multilocus, boolean boolgamma,
                         world_fmt * world, long repstart, long repstop)
{
    long nn;
    if (multilocus)
    {
        nn = boolgamma ? world->numpop2 + 1 : world->numpop2;
        if (boolgamma)
        {
            world->param0 =
                (MYREAL *) myrealloc (world->param0,
                                    sizeof (MYREAL) * (nn > 0 ? nn : 1));
            world->param0[world->numpop2] = world->options->alphavalue;

        }
    }
    else
        nn = world->numpop2;
    return nn;
}

/* profiler setup--------------------------
   needs the passed in values generated in profiles setup (profile.c)
*/
void
profile_setup (long profilenum, long *profiles, MYREAL *param, MYREAL *value,
               long *nnn)
{
    long i;
    if (profilenum > 0)
    {
        *nnn -= profilenum;
        for (i = 0; i < profilenum; i++)
            param[profiles[i]] = LOG (value[i]);
    }
    if (*nnn == 0)
        *nnn = 1;
}


void
setup_parameter0_standard (world_fmt * world, nr_fmt * nr, long repkind,
                           long repstart, long repstop, long loci, long kind,
                           boolean multilocus)
{
    long locus, r;
    if (multilocus)
    {
        for (locus = 0; locus < loci; locus++)
        {
            if (repkind == SINGLECHAIN)
            {
                for (r = repstart; r < repstop; r++)
                    create_apg0 (nr->apg0[r][locus], nr,
                                 &world->atl[r][locus], locus);
            }
            else
            {
                if (kind != PROFILE)
                {
                    for (r = repstart; r < repstop; r++)
                        create_apg0 (nr->apg0[r][locus], nr,
                                     &world->atl[r][locus], locus);
                    interpolate_like (nr, locus);
                }
                //              else
                //               {
                for (r = repstart; r < repstop; r++)
                    create_multiapg0 (nr->apg0[r][locus], nr, r, locus);
                //                }
            }
        }
    }
    else    //single locus
    {
        if (repkind == SINGLECHAIN)
        {
            for (r = repstart; r < repstop; r++)
                create_apg0 (nr->apg0[r][world->locus], nr,
                             &world->atl[r][world->locus], world->locus);
        }
        else
        {
            if (kind != PROFILE)
            {
                for (r = repstart; r < repstop; r++)
                    create_apg0 (nr->apg0[r][world->locus], nr,
                                 &world->atl[r][world->locus], world->locus);
                interpolate_like (nr, world->locus);
            }
            for (r = repstart; r < repstop; r++)
                create_multiapg0 (nr->apg0[r][world->locus], nr, r, world->locus);
        }
    }
}

void
report_broyden (MYREAL newL, MYREAL normd, long trials,
                boolean boolgamma, MYREAL theta1, MYREAL alpha,
                FILE * logfile)
{
    if (boolgamma)
        FPRINTF (stdout, "%i:%li> Log(L)=%f Norm=%f Alpha=%f Theta_1=%f\n",
                 myID, trials, newL, normd, alpha, theta1);
    else
        FPRINTF (stdout, "%i:%li> Log(L)=%f Norm=%f Theta_1=%f\n",
                 myID, trials, newL, normd, theta1);
    if (logfile != NULL)
    {
        if (boolgamma)
            FPRINTF (logfile, "%i:%li> Log(L)=%f Norm=%f Alpha=%f Theta_1=%f\n",
                     myID, trials, newL, normd, alpha, theta1);
        else
            FPRINTF (logfile, "%i:%li> Log(L)=%f Norm=%f Theta_1=%f\n",
                     myID, trials, newL, normd, theta1);
    }
}

boolean
guard_broyden (MYREAL newL, MYREAL oldL, long trials, MYREAL normd,
               MYREAL *normd20)
{
    if ((trials + 1) % 20 == 0)
    {
        if (fabs (normd - *normd20) < EPSILON)
            return TRUE;
        *normd20 = normd;
    }
    return FALSE;
}

/* for profiling: filling of index array */
void
fill_profile_index (nr_fmt * nr)
{
    long i, ii;
    long profilenum = nr->profilenum;
    if ((profilenum > 0) && profilenum < nr->partsize)
    {
        for (i = 0, ii = 0; i < nr->partsize; i++)
        {
            if (!find (i, nr->profiles, profilenum))
                nr->indeks[ii++] = i;
        }
    }
}

// maximizer
// - MPI

// likelihood
// - MPI
// - Gamma
// - custom migration matrix aware
// - profile aware
// - lrt aware [special case of cumstom migration]

// maximizer sceleton
// analystype = (SINGLELOCUS, MULTILOCUS, PROFILE)
long
maximize (MYREAL **thisparam, world_fmt * world, nr_fmt * nr, MYREAL **hess, long analystype,
          long repkind)
{
  //#ifdef MYREAL == float
  //const MYREAL eps = FLT_EPSILON ;
  //const MYREAL numbermax = FLT_MAX ;
  //#else
  const MYREAL eps = DBL_EPSILON ;
  const MYREAL numbermax = MYREAL_MAX ;
  //#endif
  boolean stop;
  long repstart, repstop;
  //  long Gmax;
  MYREAL two = 2.;
  MYREAL newL;
  MYREAL oldL;
  MYREAL lamda;
  long trials;
  MYREAL normd = 1000000;
  MYREAL normd20 = MYREAL_MAX;
  worldoption_fmt *wopt = world->options;
  long nnn = nr->partsize; // keeps track of profile or standard
  long nn;   // keeps track of the full parameter array
  
  MYREAL *param;
  //    MYREAL **hess;
  MYREAL *delta;
  MYREAL *gama;
  MYREAL *dv = NULL;
  MYREAL *gxv;
  MYREAL *gxv0;
  MYREAL *xv0;
  MYREAL *xv;
  MYREAL *expxv;
  MYREAL *temp;
  //MYREAL *tempparam;  //just pointer into temp
  //MYREAL *templparam;  //just pointer into temp
  helper_fmt helper;
  
  setdoublevec1d (&param, *thisparam, nnn);
  //  long reset_z = 0;
  helper.boolgamma = world->locus == world->loci ? wopt->gamma : FALSE;
  if (helper.boolgamma)
    param[nnn - 1] = world->options->alphavalue;
  //doublevec2d (&hess, nnn, nnn);
  reset_hess (hess, nnn);
  doublevec1d (&delta, nnn);
  doublevec1d (&gama, nnn);
  doublevec1d (&xv, nnn + 1);
  doublevec1d (&temp, 2 * (nnn + 1));
 //xcode tempparam = temp;
 //xcode templparam = temp + nnn + 1;
  helper.xv = xv;
  //  expxv = xv + nnn;
  //memcpy(expxv,param,sizeof(MYREAL)*nnn);
  doublevec1d (&gxv, nnn);
  doublevec1d (&gxv0, nnn);
  setdoublevec1d (&expxv, param, nnn);
  helper.expxv = expxv;
  set_logparam (xv, expxv, nnn);
  setdoublevec1d (&xv0, xv, nnn);
  helper.analystype = analystype;
  set_replicates (world, world->repkind, world->rep, &repstart, &repstop);
  nr->repstart = repstart;
  nr->repstop = repstop;
  nr->normd = normd;
  prepare_broyden (analystype, world, &helper.multilocus);
  nnn = multilocus_gamma_adjust (helper.multilocus, helper.boolgamma,
				 world, repstart, repstop);
#ifdef LONGSUM
  
  nnn=nr->partsize;
#endif
  
  nn = nnn;
    //if(helper.boolgamma)
    //{
    // save_a = nr->world->options->custm[nr->world->numpop2];
    // world->options->custm[nr->world->numpop2]='c';
    // world->options->custm2[nr->world->numpop2]='c';
    //}
    profile_setup (nr->profilenum, nr->profiles, xv, nr->values, &nnn);
    fill_profile_index (nr);
    helper.lamda = 0.;
    helper.dv = dv;
    calc_loci_param (nr, xv, xv0, gxv, 0., nnn);
    set_expparam (expxv, xv, nn);
    //SETUPPARAM0 (world, nr, world->repkind, repstart, repstop,
    //             world->loci, analystype, helper.multilocus);
    //   FPRINTF(stdout,"analystype is %li\n",analystype);

    fill_helper (&helper, expxv, xv, world, nr);
    oldL = CALCLIKE (&helper, expxv, xv);
    //  oldL= -MYREAL_MAX;
    fill_helper (&helper, expxv, xv, world, nr);
    CALCGRADIENT (nr, &helper, gxv);
    memcpy (gxv0, gxv, sizeof (MYREAL) * nnn);
    memcpy (xv0, xv, sizeof (MYREAL) * nn);
    setdoublevec1d (&dv, gxv, nnn);
    helper.dv = dv;
    /* Maximization cycle ---------------------------------------------*/
    for (trials = 0; trials < NTRIALS; trials++)
    {
#ifdef __MWERKS__
        eventloop ();
#endif

        newL = -MYREAL_MAX;
   //     if (helper.boolgamma)
    //        lamda = golden_section (psilo, &helper, tempparam, templparam);
        //line_search (dv, gxv, psilo, &helper, oldL, tempparam, templparam);
     //   else
     //       lamda = golden_section (psilo, &helper, tempparam, templparam);
                //line_search (dv, gxv, psilo, &helper, oldL, tempparam, templparam);
        helper.lamda = lamda = 1.;
        while ((newL < oldL || MYISNAN ((float) newL) || newL <= -numbermax) &&
                fabs (lamda) > 10. * eps)
        {
            calc_loci_param (nr, xv, xv0, dv, lamda, nn);
            set_expparam (expxv, xv, nn);
            fill_helper (&helper, expxv, xv, world, nr);
            newL = CALCLIKE (&helper, expxv, xv);
            lamda /= two;
        }
        nr->normd = normd = norm (gxv, nnn);
	//#ifdef LAGUERRE

        if (nr->world->locus >= nr->world->loci && nr->world->options->gamma)
        {
            //    initgammacat (nr->categs, expxv[nr->numpop2], expxv[0],
            //    nr->rate, nr->probcat);
            if (fabs (oldL - newL) < 1e-9)
                lamda = 0.;
        }
	//#endif
        if (world->options->verbose && newL > oldL)
            {
	      if(helper.boolgamma)
		report_broyden (newL, normd, trials, helper.boolgamma,
				expxv[0], expxv[nr->numpop2],
				world->options->logfile);
	      else
		report_broyden (newL, normd, trials, helper.boolgamma,
				expxv[0], -1.0,
				world->options->logfile);

	    }
        stop = guard_broyden (newL, oldL, trials, normd, &normd20);
        if (normd < EPSILON || stop)
            break;   // stopping criteria
        else if (oldL >= newL)
            lamda = 0;
        else
            oldL = newL;
        /* reset sec deri. matrix if lamda goes to 0 and retry */
        //      if(reset_z > 10)
        //{
        //  lamda=0;
        //  reset_z = 0;
        //}
        if (fabs (lamda) <= eps * 10.)
        {
            reset_hess (hess, nnn);
            two = -two;
            memcpy (dv, gxv, sizeof (MYREAL) * nnn);
            continue;
        }
        //swap(gxv0,gxv);
        memcpy (gxv0, gxv, sizeof (MYREAL) * nnn);
        //calc_loci_param (nr, xv, xv0, dv, lamda * two, nn);
        //set_expparam(expxv,xv,nn);
        fill_helper (&helper, expxv, xv, world, nr);
        CALCGRADIENT (nr, &helper, gxv);
        set_both_delta (nr, delta, xv, xv0, gama, gxv, gxv0, nnn);
        calc_hessian (hess, nnn, delta, gama);
        calc_dv (dv, hess, gxv, nnn);
        memcpy (xv0, xv, sizeof (MYREAL) * nn);
    }
    /* end maximizer cycle ------------------------------------*/
    /*if(helper.boolgamma)
    {
       world->options->custm[nr->world->numpop2]=save_a;
       world->options->custm2[nr->world->numpop2]=save_a;
    }*/

    finish_broyden (newL, normd, trials, world, nr, expxv, nn,
                    analystype, repkind);
    memcpy ((*thisparam), expxv, sizeof (MYREAL) * nn);
    myfree(param);
    //myfree(hess[0]);
    //myfree(hess);
    myfree(delta);
    myfree(gama);
    myfree(dv);
    myfree(gxv);
    myfree(gxv0);
    myfree(xv0);
    myfree(xv);
    myfree(expxv);
    myfree(temp);
    fflush (stderr);
    return trials;
}

void
finish_broyden (MYREAL newL, MYREAL normd, long trials,
                world_fmt * world, nr_fmt * nr, MYREAL *expxv, long nn,
                long analystype, long repkind)
{
    long loci = world->loci;

    timearchive_fmt **tyme = world->atl;
    memcpy (world->param0, expxv, sizeof (MYREAL) * nn);
    //correct_values (world, nr); // resets huge thetas to average

    //if (world->options->verbose)
    //    FPRINTF (stdout, "\n");

    //#ifdef INTEGRATEDLIKE
    //if(analystype==MULTILOCUS && world->loci==1)
    //  {
    //	analystype=SINGLELOCUS;
    //	world->locus=0;
    //   }
    //#endif

    switch (analystype)
    {
    case MULTILOCUS:
        //    tyme[0][loci].param = (MYREAL *) myrealloc (tyme[0][loci].param,
        //                                              sizeof (MYREAL) * nn);
        tyme[0][loci].param_like = newL;
        tyme[0][loci].normd = normd;
        tyme[0][loci].trials = trials;
        memmove (tyme[0][loci].param, expxv, sizeof (MYREAL) * nn);
        nr->normd = normd;
        convergence_check (world, world->options->verbose);
        break;
    case SINGLELOCUS:
        world->param_like = newL;
        if (repkind != SINGLECHAIN)
        {
	  //	  printf("%i> single locus combine replicates into rep=%li\n",myID,world->repstop);
            tyme[world->repstop][world->locus].param_like = newL;
            tyme[world->repstop][world->locus].normd = normd;
            tyme[world->repstop][world->locus].trials = trials;
            memmove (tyme[world->repstop][world->locus].param, expxv,
                     sizeof (MYREAL) * nn);
        }
        else
        {
	  //printf("%i> more loci put replicate into rep=%li\n",myID,world->rep);
            tyme[world->rep][world->locus].param_like = newL;
            tyme[world->rep][world->locus].normd = normd;
            tyme[world->rep][world->locus].trials = trials;
            memmove (tyme[world->rep][world->locus].param, expxv,
                     sizeof (MYREAL) * nn);
        }
	if(world->options->gelman)
	  convergence_check (world, world->options->verbose);
        if (world->loci == 1)
        {
            if ((!world->options->gelman &&
                    world->param_like < world->options->lcepsilon &&
                    world->options->plotnow && !world->options->simulation) ||
                    (world->options->plotnow && world->options->gelman &&
                     world->convergence->gelmanmeanmaxR[1] < GELMAN_MYSTIC_VALUE))
            {
                if ((world->options->replicate && repkind != SINGLECHAIN) ||
                        (!world->options->replicate && repkind == SINGLECHAIN))
                    create_plot (world, world->plane[world->locus], nr, 
                                       nr->atl[world->rep][0].T, SINGLELOCUSPLOT);
            }
        }
        nr->normd = normd;
        break;
    case PROFILE:
        //      memmove (nr->profparam, world->param0, sizeof (MYREAL) * nr->partsize);
        nr->llike = newL;
        nr->normd = normd;
    }

    if (world->options->verbose && world->repkind == SINGLECHAIN)
    {
        if (analystype == SINGLELOCUS)
            print_contribution (nr, nr->atl, nr->atl[world->rep][world->locus].T);
    }
}


void
set_both_delta (nr_fmt * nr,
                MYREAL *delta, MYREAL *den, MYREAL *deo,
                MYREAL *gama, MYREAL *gxv, MYREAL *gxv0, long size)
{
    long i, ii;
    if ((nr->profilenum > 0) && nr->profilenum < nr->partsize)
    {
        for (i = 0; i < size; i++)
        {
            ii = nr->indeks[i];
            delta[i] = den[ii] - deo[ii];
            gama[i] = gxv[i] - gxv0[i];
        }
    }
    else
    {
        for (i = 0; i < size; i++)
        {
            delta[i] = den[i] - deo[i];
            gama[i] = gxv[i] - gxv0[i];
        }
    }
}

/// calculates the EXP of an array \todo move this to tools.c
void
set_expparam (MYREAL *expa, MYREAL *a, long size)
{
    long i;
    for (i = 0; i < size; i++)
        expa[i] = EXP (a[i]);
}

/// calculates the log of an array \todo move this to the tools.c
void
set_logparam (MYREAL *loga, MYREAL *a, long size)
{
    long i;
    for (i = 0; i < size; i++)
        loga[i] = LOG (a[i]);
}

/// fills the helper structure with parameters and pointers to world and nr
void
fill_helper (helper_fmt * helper, MYREAL *param, MYREAL *lparam,
             world_fmt * world, nr_fmt * nr)
{
    MYREAL alpha = 0;
    MYREAL theta1 = param[0];
    MYREAL denom = 0;
    helper->boolgamma = world->locus < world->loci ?
        (helper->analystype==PROFILE ? world->options->gamma : FALSE) : world->options->gamma;
    if (helper->boolgamma)
    {
        alpha = param[nr->numpop2];
        if (alpha <= 0.0)
            error ("no no! alpha=0\n");
        denom = LGAMMA (alpha) + LOG (theta1 / alpha) * alpha;
        //  denom = LGAMMA (alpha) + LOG (1. / alpha) * alpha;
    }
    if (nr->lparam != lparam)
        memcpy (nr->lparam, lparam, sizeof (MYREAL) * nr->partsize);
    if (nr->param != param)
        memcpy (nr->param, param, sizeof (MYREAL) * nr->partsize);

    helper->nr = nr;
    helper->ll = -MYREAL_MAX;
    helper->weight = denom;
    helper->atl = world->atl;
}

void expected_param(MYREAL *param, world_fmt *world)
{
    long i, tree;
    //long j;
    tarchive_fmt *thistree;
    long T = world->atl[world->replicate][world->locus].T;
    memset(param,0,sizeof(MYREAL) * world->numpop2);
    for(tree=0; tree < T; tree++)
      {
        thistree = &world->atl[world->replicate][world->locus].tl[tree];
        for(i=0;i<world->numpop; i++)
          {
            param[i] += thistree->kt[i] / 2.;
//xcode   j or i?         for(j=world->mstart[i];i<world->mend[i]; i++)
             for(i=world->mstart[i];i<world->mend[i]; i++)
                param[i] += 1./thistree->km[i];
          }
      }
    for(i=0;i<world->numpop2; i++)
      {
        param[i] /= T;
      }
}
