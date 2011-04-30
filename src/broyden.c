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
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 $Id: broyden.c 1800 2011-01-29 13:40:00Z beerli $
-------------------------------------------------------*/
/*! \file broyden.c
Basic routines for maximum likelihood parameter estimation.

*/

#include "migration.h"
#include "tools.h"
#include "world.h"
#include "random.h"
#include "combroyden.h"
#include "joint-chains.h"
#include "sighandler.h"
#include "migrate_mpi.h"

#ifdef ALTIVEC
#include "altivec.h"
#endif

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

extern MYREAL norm (MYREAL *d, long size);


/* prototypes ----------------------------------------- */
MYREAL absmaxvec (MYREAL *v, long n);
void create_nr (nr_fmt * nr, world_fmt * world, long G, long profilenum,
                long thislocus, long repkind, long rep);
void reset_hess (MYREAL **hess, long n);
void destroy_nr (nr_fmt * nr, world_fmt * world);
MYREAL calc_locus_like (nr_fmt * nr, MYREAL *param, MYREAL *lparam,
                        long locus);
void param_all_adjust (MYREAL *xv, nr_fmt *nr); //worldoption_fmt * wopt, long numpop);
#ifdef LONGSUM
void gradient_longsum (MYREAL *d, nr_fmt * nr, long locus);
#else /*LONGSUM*/
void gradient (MYREAL *d, nr_fmt * nr, long locus);
#endif /*LONGSUM*/
MYREAL probG (MYREAL *param, MYREAL *lparam, tarchive_fmt * tl, nr_fmt * nr,
              long locus);
MYREAL probG2 (MYREAL *param, MYREAL *lparam, MYREAL *sm, MYREAL *kt,
               MYREAL *km, MYREAL *p, MYREAL *mindex, int *msta, int *me,
               long numpop);
void probG3 (MYREAL *apg, MYREAL *apg0r, timearchive_fmt * tyme, long numpop,
             MYREAL *apgmax, MYREAL *param, MYREAL *lparam, MYREAL *sm);
MYINLINE  MYREAL probG4 (MYREAL *fullparam, MYREAL *data, int stop);
MYREAL sum_mig (MYREAL *param, long msta, long msto);
void calc_cov (MYREAL **dd, MYREAL *d, MYREAL *param, long n);
boolean is_singular (MYREAL **dd, long n);
void print_contribution (nr_fmt * nr, timearchive_fmt ** atl, long G);

void calc_dv (MYREAL *dv, MYREAL **hess, MYREAL *gxv, long n);
MYREAL calc_line (helper_fmt * helper, MYREAL a, MYREAL b, MYREAL c,
                  MYREAL (*psi) (MYREAL lamda, helper_fmt * helper));
void calc_hessian (MYREAL **hess, long n, MYREAL *delta, MYREAL *gama);
MYREAL psi (MYREAL lamda, helper_fmt * helper, MYREAL *param, MYREAL *lparam);
void grad2loggrad (MYREAL *param, long *indeks, MYREAL *d, long nn,
                   long profilenum);
void log_param0 (MYREAL *param, MYREAL *lparam, long nn);
void copies2lcopies (timearchive_fmt * atl);
void create_apg0 (MYREAL *apg0, nr_fmt * nr, timearchive_fmt * tyme,
                  long locus);


void force_sametheta (MYREAL *param, worldoption_fmt * wopt, long numpop);
void force_samemigration (MYREAL *param, worldoption_fmt * wopt, long numpop);
void calc_same_d (MYREAL *grad, nr_fmt * nr, long nn, long start, long stop);
void calc_symmetric_d (MYREAL *grad, nr_fmt * nr, long nn, long start,
                       long stop);
void force_symmetric_d (MYREAL *gxv, long model, nr_fmt * nr, long nn);
void check_symmetric_d (MYREAL *gxv, nr_fmt * nr, long nn);
void check_matrix_arbitrary (MYREAL *param, worldoption_fmt * wopt,
                             long numpop, long which_profile);
void print_menu_contribution (FILE * file, long contribution[]);
void alloc_apg (MYREAL ****apg, long repstop, long loci, long G);


void quadratic_constants (MYREAL *xguess,
                          MYREAL low, MYREAL mid, MYREAL high,
                          MYREAL xlow, MYREAL xmid, MYREAL xhigh);
void symmetric_other (long i, long numpop, long *other, long *otherpop);
MYREAL ln_copies (long n);

MYREAL normal_func_ok(MYREAL *param, MYREAL *param0, long numpop2);
MYREAL normal_func_no(MYREAL *param, MYREAL *param0, long numpop2);
MYREAL normal_func_gradient_ok(MYREAL p1, MYREAL p0);
MYREAL normal_func_gradient_no(MYREAL p1, MYREAL p0);

void unset_penalizer_function(boolean inprofiles);

MYREAL (*normal_func)(MYREAL *, MYREAL *, long);
MYREAL (*normal_func_gradient)(MYREAL, MYREAL);

#ifdef LONGSUM
MYREAL probG_longsum (MYREAL *param, MYREAL *lparam, tarchive_fmt * tl, nr_fmt * nr, long locus);
MYREAL probG4_longsum (MYREAL *fullparam, longsum_fmt *data, long longsumlen, long numpop, long numpop2);
#endif /*LONGSUM*/

/* Functions ++++++++++++++++++++++++++++++++++++++++++++++++*/

///
/// Calculates part of the BFSG maximizer routine.
/// Multiplies the approximate hessian with the first derivative vector
/// and returns the product
void
calc_dv (MYREAL *dv, MYREAL **hess, MYREAL *gxv, long n)
{
    long i, j;
    memset (dv, 0, sizeof (MYREAL) * n);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            dv[i] += hess[i][j] * gxv[j];
        }
    }
}

// line searcher
// finds the maximum in a direction
// this should be replaced with something more efficient.
#define PP 0.61803399
#define QQ 0.38196601
#define MOVE3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)




MYREAL
calc_line_new (helper_fmt * helper, MYREAL a, MYREAL b, MYREAL c,
               MYREAL (*fpsi) (MYREAL lamda, helper_fmt * helper))
{
    MYREAL xhigh = c;
    MYREAL xlow = a;
    MYREAL low, high, mid, xmid;
    MYREAL xguess=0.;
    int panic = 0;
    low = (*fpsi) (xlow, helper);
    high = (*fpsi) (xhigh, helper);
    xmid = (xlow + xhigh) * HALF;
    while (panic++ < 10000 && fabs (xhigh - xlow) > EPSILON)
    {
        mid = (*fpsi) (xmid, helper);
        quadratic_constants (&xguess, low, mid, high, xlow, xmid, xhigh);
        if (MYISNAN ((float) xguess))
            return 1.;
        if (xguess < xmid)
        {
            xhigh = xmid;
            high = mid;
        }
        else
        {
            xlow = xmid;
            low = mid;
        }
        xmid = xguess;
    }
    return xlow;
}


void
quadratic_constants (MYREAL *xguess,
                     MYREAL low, MYREAL mid, MYREAL high,
                     MYREAL xlow, MYREAL xmid, MYREAL xhigh)
{
    MYREAL xhigh2 = xhigh * xhigh;
    MYREAL xlow2 = xlow * xlow;
    MYREAL xmid2 = xmid * xmid;
    MYREAL divisor =        ((xlow - xmid) * (xlow * xmid - xlow * xhigh - xmid * xhigh + xhigh2));
    MYREAL a = -((-(xmid * low) + xhigh * low + xlow * mid - xhigh * mid -
                  xlow * high + xmid * high) /divisor);
    MYREAL b =
        -((xmid2 * low - xhigh2 * low - xlow2 * mid + xhigh2 * mid +
           xlow2 * high - xmid2 * high) / divisor);
    //    MYREAL c = -((-(xmid2*xhigh*low) + xmid*xhigh2*low +
    //    xlow2*xhigh*mid - xlow*xhigh2*mid -
    //    xlow2*xmid*high + xlow*xmid2*high)/
    //  (xlow2*xmid - xlow*xmid2 - xlow2*xhigh +
    //   xmid2*xhigh + xlow*xhigh2 - xmid*xhigh2));
    *xguess = -b / (2. * a);
    // printf("quadr:{%f %f}  {%f %f} {%f %f}\n", xlow, low, *xguess, a * (*xguess * *xguess) + b * *xguess + c, xhigh, high);
}



MYREAL
calc_line (helper_fmt * helper, MYREAL a, MYREAL b, MYREAL c,
           MYREAL (*fpsi) (MYREAL lamda, helper_fmt * helper))
{
    /* a < b < c AND psia > psib < psic */
    MYREAL d, psib, psic;
    d = c;
    if ((fabs (c - b)) > (fabs (b - a)))
    {
        c = b + QQ * (c - b);
    }
    else
    {
        c = b;
        b = b - QQ * (b - a);
    }
    psib = (*fpsi) (b, helper);
    psic = (*fpsi) (c, helper);
    while ((fabs (d - a)) > (EPSILON * (fabs (b) + fabs (c))))
    {
        if (psic < psib)
        {
            MOVE3 (a, b, c, PP * b + QQ * d);
            psib = psic;
            psic = (*fpsi) (c, helper);
        }
        else
        {
            MOVE3 (d, c, b, PP * c + QQ * a);
            psic = psib;
            psib = (*fpsi) (b, helper);
        }
        //      printf("b=%f(%f
        //), c=%f(%f) {%f, %f}\n",b,psib,c,psic, helper->nr->param[0],
        //     helper->nr->param[helper->nr->partsize-1]);
    }
    if (psib < psic)
    {
        return b;
    }
    else
    {
        return c;
    }
}

MYREAL
psi (MYREAL lamda, helper_fmt * helper, MYREAL *param, MYREAL *lparam)
{
    MYREAL like;
    calc_loci_param (helper->nr, helper->nr->lparam, helper->xv,
                     helper->dv, lamda, helper->nr->partsize);
    set_expparam (helper->nr->param, helper->nr->lparam, helper->nr->partsize);
    fill_helper (helper, helper->nr->param, helper->nr->lparam,
                 helper->nr->world, helper->nr);
    like = CALCLIKE (helper, helper->nr->param, helper->nr->lparam);
    return -like;
}


void
calc_hessian (MYREAL **hess, long n, MYREAL *delta, MYREAL *gama)
{
    MYREAL **dd, *temp, t;
    long i, j, k;
    MYREAL dtg;
    temp = (MYREAL *) mycalloc (n, sizeof (MYREAL));
    dd = (MYREAL **) mycalloc (n, sizeof (MYREAL *));
    dd[0] = (MYREAL *) mycalloc (n * n, sizeof (MYREAL));
    dtg = delta[0] * gama[0];
    for (i = 1; i < n; i++)
    {
        dd[i] = dd[0] + n * i;
        dtg += delta[i] * gama[i];
    }
    if (dtg > 0.0 || dtg < 0.0)
        dtg = 1. / dtg;
    else
    {
        reset_hess (hess, n);
        myfree(temp);
        myfree(dd[0]);
        myfree(dd);
        return;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            temp[i] += gama[j] * hess[j][i];
        }
    }
    t = 0.0;
    for (i = 0; i < n; i++)
        t += temp[i] * gama[i];
    t = (1.0 + t * dtg) * dtg;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                dd[i][j] += delta[i] * gama[k] * hess[k][j];
            }
        }
    }
    for (i = 0; i < n; i++)
    {
        temp[i] = 0.0;
        for (j = 0; j < n; j++)
        {
            temp[i] += hess[i][j] * gama[j];
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            dd[i][j] += temp[i] * delta[j];
            dd[i][j] *= dtg;
            hess[i][j] += delta[i] * delta[j] * t - dd[i][j];
        }
    }
    myfree(temp);
    myfree(dd[0]);
    myfree(dd);
}

MYREAL
absmaxvec (MYREAL *v, long n)
{
    long i;
    MYREAL max = fabs (v[0]);
    for (i = 1; i < n; i++)
    {
        if (max < fabs (v[i]))
            max = fabs (v[i]);
    }
    return max;
}

///
/// creates the scratchpad for the maximum likelihood calculation, used during the MLE calculations
/// and also used during the profiles, the structure nr holds many parameters for further use in other functions
/// 
void
create_nr (nr_fmt * nr, world_fmt * world, long G,
           long profilenum, long thislocus, long repkind, long rep)
{
    long i;
    long partsize;
    // setting up local variables
    nr->numpop = world->numpop;
    nr->numpop2 = world->numpop2;
    nr->skiploci = world->data->skiploci;
    nr->world = world;
    nr->mstart = world->mstart;
    nr->mend = world->mend;
    nr->repkind = repkind;
    nr->repstart = (repkind == SINGLECHAIN) ? rep : 0;
    nr->repstop = (repkind == SINGLECHAIN) ? rep + 1 : world->repstop;
    nr->profilenum = profilenum;
    // number of trees stored
    nr->apg_max = (MYREAL *) mycalloc (2 * (world->loci + 1), sizeof (MYREAL));
    nr->PGC = nr->apg_max + world->loci + 1; //memory alocated in apg_max
    if (world->options->gamma)
    {
        nr->categs = GAMMA_INTERVALS;
        nr->rate = (MYREAL *) mycalloc (2 * nr->categs, sizeof (MYREAL));
        nr->probcat = nr->rate + nr->categs; // memory allocated in nr->rate
        nr->partsize = world->numpop2 + 1;
    }
    else
    {
        if (!world->options->gamma && world->locus >= world->loci)
        {
            nr->categs = 0;
            nr->rate = NULL;
            nr->partsize = world->numpop2;
        }
        else
        {
            nr->categs = 0;
            nr->rate = NULL;
            nr->partsize = world->numpop2;
        }
    }
#ifdef LONGSUM
    nr->partsize += world->numpop * 3;
#endif
    partsize = nr->partsize;

    nr->atl = world->atl;

    nr->datalike = (MYREAL *) mymalloc ((G+1) * sizeof (MYREAL));
    
    // parts is a long vector containing all standard parts that have a fixed length
    // this excludes the datalikes because those use G (the number of trees), they might get reallocated during the
    // lifetime of nr and so have their own memory block.
    nr->parts = (MYREAL *) mycalloc (9 * partsize + (world->loci + 1), sizeof (MYREAL));
    nr->d = nr->parts + partsize;           // memory allocated in nr->parts
    nr->od = nr->d + partsize;              // memory allocated in nr->parts
    nr->dv = nr->od + partsize;             // memory allocated in nr->parts
    nr->delta = nr->dv + partsize;          // memory allocated in nr->parts
    nr->gdelta = nr->delta + partsize;      // memory allocated in nr->parts
    nr->param = nr->gdelta + partsize;      // memory allocated in nr->parts
    nr->lparam = nr->param + partsize ;     // memory allocated in nr->parts
    nr->values = nr->lparam  + partsize;    // memory allocated in nr->parts
    nr->locilikes = nr->values + partsize;  // memory allocated in nr->parts

    nr->profiles = (long *) mycalloc (2 * partsize, sizeof (long));//holds vector for profiles and indeks
    nr->indeks = nr->profiles + partsize;  //memory allocated in nr->profiles
    
    // assignment of parameters into scratch pad
    memcpy (nr->param, world->param0, partsize * sizeof (MYREAL));
    for (i = 0; i < nr->partsize; ++i)
    {
        if (nr->param[i] > 0.0)
            nr->lparam[i] = LOG (nr->param[i]);
        else
            nr->lparam[i] = -MYREAL_MAX;
    }
    // allocate earlier in world: alloc_apg (&nr->apg0, world->repstop, world->loci, G);
    // allicate earlier in world:  alloc_apg (&nr->apg, world->repstop, world->loci, G);
    nr->apg = world->apg;
    nr->apg0 = world->apg0;
}

void
alloc_apg (MYREAL ****apg, long repstop, long loci, long G)
{
    long j, i;
    (*apg) = (MYREAL ***) mycalloc (repstop, sizeof (MYREAL **));
    for (j = 0; j < repstop; j++)
    {
        (*apg)[j] = (MYREAL **) mycalloc (loci + 1, sizeof (MYREAL *));
        (*apg)[j][0] =
            (MYREAL *) mycalloc ((1 + G) * (loci + 1), sizeof (MYREAL));
        for (i = 1; i < loci + 1; i++)
        {
            (*apg)[j][i] = (*apg)[j][0] + i * (G + 1);
        }
    }
}

void
reuse_nr (nr_fmt * nr, world_fmt * world, long G, long profilenum,
          long thislocus, long repkind, long rep)
{
    long i;
    nr->repkind = repkind;
    nr->repstart = (repkind == SINGLECHAIN) ? rep : 0;
    nr->repstop = (repkind == SINGLECHAIN) ? rep + 1 : world->repstop;
    nr->profilenum = profilenum;
    if (world->options->gamma && world->locus >= world->loci)
    {
        nr->categs = GAMMA_INTERVALS;
        nr->rate = (MYREAL *) myrealloc (nr->rate, 2 * nr->categs * sizeof (MYREAL));
        nr->probcat = nr->rate + nr->categs;
        nr->partsize = world->numpop2 + 1;
    }
    memset (nr->parts, 0, nr->partsize * sizeof (MYREAL));
    memset (nr->d, 0, nr->partsize * sizeof (MYREAL));
    memset (nr->od, 0, nr->partsize * sizeof (MYREAL));
    memset (nr->dv, 0, nr->partsize * sizeof (MYREAL));
    memset (nr->delta, 0, nr->partsize * sizeof (MYREAL));
    memset (nr->gdelta, 0, nr->partsize * sizeof (MYREAL));
    memcpy (nr->param, world->param0, nr->partsize * sizeof (MYREAL));
    for (i = 0; i < nr->partsize; ++i)
    {
        if (nr->param[i] > 0.0)
            nr->lparam[i] = LOG (nr->param[i]);
        else
            nr->lparam[i] = -MYREAL_MAX;
    }
    nr->datalike = (MYREAL *) myrealloc (nr->datalike, (G+1) * sizeof (MYREAL));
    memset (nr->datalike, 0, G * sizeof (MYREAL));
    memset (nr->locilikes, 0, (world->loci + 1) * sizeof (MYREAL));
    /*
     for (j = world->repstop - 1; j >= 0; j--)
     {
      myfree(nr->apg0[j][0]);
      myfree(nr->apg[j][0]);
      myfree(nr->apg0[j]);
      myfree(nr->apg[j]);
     }
     myfree(nr->apg0);
     myfree(nr->apg);
     alloc_apg (&nr->apg0, world->repstop, world->loci, G);
     alloc_apg (&nr->apg, world->repstop, world->loci, G);*/
}

void
reset_hess (MYREAL **hess, long n)
{
    long pop;
    //  printf("resetting hessian\n");
    memset (hess[0], 0, sizeof (MYREAL) * n * n);
    for (pop = 1; pop < n; pop++)
    {
        hess[pop] = hess[0] + pop * n;
        hess[pop][pop] = 1.0;
    }
    hess[0][0] = 1.0;
}

void
destroy_nr (nr_fmt * nr, world_fmt * world)
{
    // the large holding vectors are freed and all related points will go away, too.
    myfree(nr->parts);
    myfree(nr->datalike);
    myfree(nr->apg_max); 
    myfree(nr->profiles);

    if (nr->categs != 0)
    {
        myfree(nr->rate); // probcat is also freed, was part of rate-array
    }
    myfree(nr);
}

///
/// calculate the likelihood of the parameter estimates
#ifndef ALTIVEC
MYREAL
calc_locus_like (nr_fmt * nr, MYREAL *param, MYREAL *lparam, long locus)
{
    long g, r, j, copies;
    const long numpop = nr->world->numpop;
    const long numpop2 = nr->world->numpop2;
    int msta, msto;
    int stop = (int) (numpop2 + numpop + numpop);
    MYREAL gsum = 0;
    MYREAL ***apg0;
    MYREAL apgmax = -MYREAL_MAX;
    //  MYREAL ***apgall = nr->apg;
    MYREAL *apg;
    MYREAL *apg0r;
    timearchive_fmt tyme;
    tarchive_fmt *tl;
    MYREAL *geo = nr->world->data->geo;
    MYREAL *lgeo = nr->world->data->lgeo;
    MYREAL *locallparam;
    MYREAL *localparam;
    MYREAL *sm;
    MYREAL mu_rate = nr->world->options->mu_rates[locus];
    MYREAL inv_mu_rate = 1./ mu_rate;
    MYREAL lmu_rate = nr->world->options->lmu_rates[locus];
    MYREAL normaldev;
#ifdef LONGSUM

    MYREAL *rates;
    MYREAL *lrates;
    MYREAL *rtimes;
#endif /*LONGSUM*/
#ifndef LONGSUM
    // experiment with unrolling the loop in probG4 by 3, the 3 added here guarantees
    // that enough elements are present so that the loop can gracefully finish
#ifdef FAST_TRY
    stop += 3 - (stop % 3);
#endif
    locallparam = (MYREAL *) mycalloc (stop, sizeof (MYREAL));
#else /*LONGSUM*/

    locallparam = (MYREAL *) mycalloc ((numpop2 + numpop + numpop + 9 * numpop),
                                     sizeof (MYREAL));
    rates = locallparam + numpop2 + numpop + numpop;
    lrates = rates + 3 * numpop;
    rtimes = rates +  6 * numpop;
    memcpy (rates, param+nr->partsize-3*numpop,sizeof(MYREAL)*3*numpop);
    memcpy (lrates, lparam+nr->partsize-3*numpop,sizeof(MYREAL)*3*numpop);
    memcpy (rtimes, nr->world->flucrates+ 3 * numpop ,sizeof(MYREAL)*3*numpop);

#endif /*LONGSUM*/

    localparam = locallparam + numpop2;

    sm = localparam + numpop;
    memcpy (localparam, param, sizeof (MYREAL) * numpop);
    memcpy (locallparam, lparam, sizeof (MYREAL) * numpop2);

    nr->PGC[locus] = 0.0;
    apg0 = nr->apg0;

    for (r = 0; r < nr->numpop; r++)
    {
        locallparam[r] = LOG2 - (locallparam[r] + lmu_rate);
        localparam[r] = -1. / (localparam[r] * mu_rate); // minus, so that we can loop over all in probG4
        msta = nr->mstart[r];
        msto = nr->mend[r];
        for (j = msta; j < msto; j++)
        {
            sm[r] -= geo[j] * param[j] * inv_mu_rate; //minus, so that we can loop over all in probG4
            locallparam[j] += lgeo[j] - lmu_rate;
        }
    }
    gsum = 0;
    for (r = nr->repstart; r < nr->repstop; r++)
    {
        apg = nr->apg[r][locus];
        apg0r = apg0[r][locus];
        tyme = nr->atl[r][locus];
        normaldev =  (*normal_func)(param, tyme.param0, numpop2);
        //first element
        tl = &(tyme.tl[0]);
        copies = tl->copies - 1;
        gsum += copies;
#ifdef INTEGRATEDLIKE
	tl->lcopies = ln_copies(copies);
#endif
        if (copies > 0)
        {
#ifndef LONGSUM
            apg[0] = tl->lcopies + probG4 (locallparam, tl->data, stop) - apg0r[0] + normaldev;
#else /*LONGSUM*/

            apg[0] = tl->lcopies + probG4_longsum(locallparam,tl->longsum,tl->longsumlen, numpop,numpop2)
                     -apg0r[0] + normaldev;
#endif /*LONGSUM*/

        }
        else
            apg[0] = -MYREAL_MAX;
        if (apgmax < apg[0])
            apgmax = apg[0];
        //other elements
        for (g = 1; g < tyme.T; g++)
        {
            tl = &(tyme.tl[g]);
            gsum += tl->copies;
#ifdef INTEGRATEDLIKE
	tl->lcopies = ln_copies(tl->copies);
#endif

#ifndef LONGSUM

            apg[g] = tl->lcopies + probG4 (locallparam, tl->data, stop) - apg0r[g]
                     +normaldev;
#else /*LONGSUM*/

            apg[g] = tl->lcopies + probG4_longsum(locallparam,tl->longsum,tl->longsumlen, numpop,numpop2)
                     -apg0r[g]   + normaldev;
#endif /*LONGSUM*/

            if (apg[g] > apgmax)
                apgmax = apg[g];
        }
    }    // end replicates
    for (r = nr->repstart; r < nr->repstop; r++)
    {
        apg = nr->apg[r][locus];
        //apg0r = apg0[r][locus];
        tyme = nr->atl[r][locus];
        // first element
   //xcode     tl = &(tyme.tl[0]);
        // first element
        apg[0] -= apgmax;
        nr->PGC[locus] += EXP (apg[0]);
        // all other elements
        for (g = 1; g < tyme.T; g++)
        {
            apg[g] -= apgmax;
            if (apg[g] > -40.)
                nr->PGC[locus] += EXP (apg[g]);
        }
    }    // replicates
    nr->apg_max[locus] = apgmax;
    nr->llike = apgmax + LOG (nr->PGC[locus]) - LOG (gsum);
    myfree(locallparam);
    return nr->llike;
}
#else /*ALTIVEC*/

MYREAL
calc_locus_like (nr_fmt * nr, MYREAL *param, MYREAL *lparam, long locus)
{
    long g, r, j, copies;
    const long numpop = nr->world->numpop;
    const long numpop2 = nr->world->numpop2;
    int msta, msto;
    MYREAL gsum = 0;
    MYREAL ***apg0;
    MYREAL apgmax = -MYREAL_MAX;
    //  MYREAL ***apgall = nr->apg;
    MYREAL *apg;
    MYREAL *apg0r;
    timearchive_fmt tyme;
    tarchive_fmt *tl;
    MYREAL *geo = nr->world->data->geo;
    MYREAL *lgeo = nr->world->data->lgeo;
    MYREAL *locallparam;
    MYREAL *localparam;
    MYREAL *sm;
    MYREAL mu_rate = nr->world->options->mu_rates[locus];
    MYREAL inv_mu_rate = 1./ mu_rate;
    MYREAL lmu_rate = nr->world->options->lmu_rates[locus];
    MYREAL normaldev;
    long size4 = numpop2 + numpop + numpop;
    long size;
    FloatVec *vlocal;
    FloatVec *vdata;
    size4 -= (size4 % 4) - 4;
    size = size4 * QUARTER;
    vlocal = (FloatVec *) mycalloc (MAX(16,size), sizeof(FloatVec));
    vdata = (FloatVec *) mycalloc (MAX(16,size), sizeof(FloatVec));
    locallparam = (MYREAL *) mycalloc ((size4), sizeof (MYREAL));
    localparam = locallparam + numpop2;

    sm = localparam + numpop;
    memcpy (localparam, param, sizeof (MYREAL) * numpop);
    memcpy (locallparam, lparam, sizeof (MYREAL) * numpop2);

    nr->PGC[locus] = 0.0;
    apg0 = nr->apg0;

    for (r = 0; r < nr->numpop; r++)
    {
        locallparam[r] = LOG2 - (locallparam[r] + lmu_rate);
        localparam[r] = -1. / (localparam[r] * mu_rate); // minus, so that we can loop over all in probG4
        msta = nr->mstart[r];
        msto = nr->mend[r];
        for (j = msta; j < msto; j++)
        {
            sm[r] -= geo[j] * param[j] * inv_mu_rate; //minus, so that we can loop over all in probG4
            locallparam[j] += lgeo[j] - lmu_rate;
        }
    }
    gsum = 0;
    load_float_floatvec(vlocal, locallparam,numpop2+numpop+numpop);
    for (r = nr->repstart; r < nr->repstop; r++)
    {
        apg = nr->apg[r][locus];
        apg0r = apg0[r][locus];
        tyme = nr->atl[r][locus];
        normaldev =  (*normal_func)(param, tyme.param0, numpop2);
        //first element
        tl = &(tyme.tl[0]);
        copies = tl->copies - 1;
        gsum += copies;
        if (copies > 0)
        {
//            load_float_floatvec(vdata, tl->data,numpop2+numpop+numpop);
            apg[0] = tl->lcopies + vdot_product (vlocal, tl->vdata, size) - apg0r[0] + normaldev;
        }
        else
            apg[0] = -MYREAL_MAX;
        if (apgmax < apg[0])
            apgmax = apg[0];
        //other elements
        for (g = 1; g < tyme.T; g++)
        {
            tl = &(tyme.tl[g]);
            gsum += tl->copies;
            load_float_floatvec(vdata, tl->data,numpop2+numpop+numpop);
            apg[g] = tl->lcopies + vdot_product (vlocal,vdata, size) - apg0r[g] + normaldev;
            if (apg[g] > apgmax)
                apgmax = apg[g];
        }
    }    // end replicates
    for (r = nr->repstart; r < nr->repstop; r++)
    {
        apg = nr->apg[r][locus];
        apg0r = apg0[r][locus];
        tyme = nr->atl[r][locus];
        // first element
        tl = &(tyme.tl[0]);
        // first element
        apg[0] -= apgmax;
        nr->PGC[locus] += EXP (apg[0]);
        // all other elements
        for (g = 1; g < tyme.T; g++)
        {
            apg[g] -= apgmax;
            if (apg[g] > -40.)
                nr->PGC[locus] += EXP (apg[g]);
        }
    }    // replicates
    nr->apg_max[locus] = apgmax;
    nr->llike = apgmax + LOG (nr->PGC[locus]) - LOG (gsum);
    myfree(locallparam);
    myfree(vlocal);
    myfree(vdata);
    return nr->llike;
}
#endif /*ALTIVEC*/



///
/// Lookup table for the calculation of log(copy_number_of_trees)
/// Most often there are only a few copies of a specific tree.
MYREAL
ln_copies (long n)
{
    switch (n)
    {
    case 0:
        return -MYREAL_MAX;
    case 1:
        return 0.;
    case 2:
        return 0.69314718055994530942;
    case 3:
        return 1.0986122886681096914;
    case 4:
        return 1.3862943611198906188;
    case 5:
        return 1.6094379124341003746;
    case 6:
        return 1.7917594692280550008;
    case 7:
        return 1.9459101490553133051;
    case 8:
        return 2.0794415416798359283;
    case 9:
        return 2.1972245773362193828;
    case 10:
        return 2.3025850929940456840;
    case 11:
        return 2.3978952727983705441;
    case 12:
        return 2.4849066497880003102;
    case 13:
        return 2.5649493574615367361;
    case 14:
        return 2.6390573296152586145;
    case 15:
        return 2.7080502011022100660;
    case 16:
        return 2.7725887222397812377;
    case 17:
        return 2.8332133440562160802;
    case 18:
        return 2.8903717578961646922;
    case 19:
        return 2.9444389791664404600;
    case 20:
        return 2.9957322735539909934;
    default:
        return LOG ((MYREAL) n);
    }
}

/// make sure that all 'm' marked population sizes are all the same
void
force_sametheta (MYREAL *param, worldoption_fmt * wopt, long numpop)
{
    long i, n = 0;
    MYREAL summ = 0, logsumm;
	char *custm2 = wopt->custm2;
    for (i = 0; i < numpop; i++)
    {
        if (custm2[i] == 'm')
        {
            summ += EXP (param[i]);
            n++;
        }
    }
    summ /= n;
    logsumm = LOG (summ);
    for (i = 0; i < numpop; i++)
    {
        if (custm2[i] == 'm')
            param[i] = logsumm;
    }
}

/// make sure that all 'm' marked migration rates are the same
void
force_samemigration (MYREAL *param, worldoption_fmt * wopt, long numpop)
{
    long i;
    MYREAL summ = 0, logsumm;
    long numpop2 = numpop * numpop;
    long n = 0;
	char *custm2 = wopt->custm2;

    for (i = numpop; i < numpop2; i++)
    {
        if (custm2[i] == 'm')
        {
            summ += EXP (param[i]);
            n++;
        }
    }
    summ /= n;
    logsumm = LOG (summ);
    for (i = numpop; i < numpop2; i++)
    {
        if (custm2[i] == 'm')
            param[i] = logsumm;
    }
}

/// adjust all parameters according to their settings in custom migration matrix
void
param_all_adjust (MYREAL *xv, nr_fmt *nr)
{
    MYREAL *param = xv;
    worldoption_fmt * wopt= nr->world->options;
    long numpop = nr->world->numpop;
    MYREAL which_profile= -1;
    char *custm = NULL;
    long zeron = 0;
    long constn = 0;
    long symn = 0;
    long sym2n = 0;
    long from, to;
    twin_fmt *syms = NULL;
    quad_fmt *sym2s = NULL;
    long *zeros = NULL;
    long *consts = NULL;
    long z = 0, i;
    char *p;
    MYREAL mm, lmm;

    custm = wopt->custm2;
    zeron = wopt->zeron;
    constn = wopt->constn;
    symn = wopt->symn;
    sym2n = wopt->sym2n;

    if (symn > 0)
        syms = wopt->symparam;
    if (sym2n > 0)
        sym2s = wopt->sym2param;
    if (zeron > 0)
        zeros = wopt->zeroparam;
    if (constn > 0)
        consts = wopt->constparam;

    p = custm;
    z = (long) strcspn(p, "m");
    if (z < numpop)
        force_sametheta (param, wopt, numpop);
    p = custm + numpop;
    z = (long) strcspn(p, "m");
    
	if (z < numpop*(numpop-1))
        force_samemigration (param, wopt, numpop);
    
	for (i = 0; i < zeron; i++)
        param[zeros[i]] = -30.; //-MYREAL_MAX;
    
	for (i = 0; i < constn; i++)
    {
        if (consts[i] < numpop)
            param[consts[i]] = LOG (wopt->thetag[consts[i]]);
        if (consts[i] >= numpop)
        {
            if (consts[i] >= numpop * numpop)
                param[consts[i]] = LOG (wopt->alphavalue);
            else
            {
                m2mm (consts[i], numpop, &from, &to);
                param[consts[i]] = LOG (wopt->mg[to*(numpop-1)+from]);
            }
        }
    }

    // this is so weird because of the profiled parameters which should not change of course
    for (i = 0; i < symn; i++)
    {
        if(nr->profilenum>0)
        {
            for(z=0;z<nr->profilenum;z++)
            {
                which_profile = nr->profiles[z];
                if(which_profile == (long) syms[i][0])
                {
                    param[syms[i][1]] = param[syms[i][0]];
                    break;
                }
                else
                {
                    if(which_profile == (long) syms[i][1])
                    {
                        param[syms[i][0]] = param[syms[i][1]];
                        break;
                    }
                    else
                    {
                        mm = (EXP (param[syms[i][0]]) + EXP (param[syms[i][1]]))  * HALF;
                        param[syms[i][0]] = param[syms[i][1]] = LOG (mm);
                        break;
                    }
                }
            }
        }
        else
        {
            mm = (EXP (param[syms[i][0]]) + EXP (param[syms[i][1]])) * HALF;
            param[syms[i][0]] = param[syms[i][1]] = LOG (mm);
        }
    }
    for (i = 0; i < sym2n; i++)
    {
        if(nr->profilenum>0)
        {
            for(z=0;z<nr->profilenum;z++)
            {
                which_profile = nr->profiles[z];
                if(which_profile == (long) sym2s[i][0])
                {
                    //printf(".1.\n");
                    mm = param[sym2s[i][0]] + param[sym2s[i][2]];
                    param[sym2s[i][1]] = mm - param[sym2s[i][3]];

                    break;
                }
                else
                {
                    if(which_profile == (long) sym2s[i][1])
                    {
                        //printf(".2.\n");
                        mm = param[sym2s[i][1]] + param[sym2s[i][3]];
                        param[sym2s[i][0]] = mm - param[sym2s[i][2]];
                        break;
                    }
                    else
                    {
                        mm = (EXP (param[sym2s[i][0]] + param[sym2s[i][2]]) +
                              EXP (param[sym2s[i][1]] + param[sym2s[i][3]])) * HALF;
                        param[sym2s[i][0]] = (lmm = LOG (mm)) - param[sym2s[i][2]];
                        param[sym2s[i][1]] = lmm - param[sym2s[i][3]];
                        break;
                    }
                }
            }
        }
        else
        {
            mm = (EXP (param[sym2s[i][0]] + param[sym2s[i][2]]) +
                  EXP (param[sym2s[i][1]] + param[sym2s[i][3]])) * HALF;
            param[sym2s[i][0]] = (lmm = LOG (mm)) - param[sym2s[i][2]];
            param[sym2s[i][1]] = lmm - param[sym2s[i][3]];
        }
    }
}



//======================
/*  switch (wopt->migration_model)
{
    case MATRIX:
  break;
    case MATRIX_ARBITRARY:
  check_matrix_arbitrary (xv, wopt, numpop, which_profile);
  break;
    case MATRIX_SAMETHETA:
  force_sametheta (xv, wopt, numpop);
  break;
    case ISLAND:
  force_sametheta (xv, wopt, numpop);
  force_samemigration (xv, wopt, numpop);
  break;
    case ISLAND_VARTHETA:
  force_samemigration (xv, wopt, numpop);
  break;
}
}
*/

///
/// checks and corrects the migration matrix according to custom migration matrix
/// which_profile is necessary to avoid accidental changes in the variables that need profiling
void
check_matrix_arbitrary (MYREAL *param, worldoption_fmt * wopt, long numpop, long which_profile)
{
    char *custm = NULL;
    long zeron = 0;
    long constn = 0;
    long symn = 0;
    long sym2n = 0;
    long from, to;
    //   boolean done = FALSE;
    twin_fmt *syms = NULL;
    quad_fmt *sym2s = NULL;
    long *zeros = NULL;
    long *consts = NULL;
    long z = 0, i;
    long count;
    char *p;
    MYREAL mm, lmm;
    //     long locus = 0;
    custm = wopt->custm2;
    zeron = wopt->zeron;
    constn = wopt->constn;
    symn = wopt->symn;
    sym2n = wopt->sym2n;
    if (symn > 0)
        syms = wopt->symparam;
    if (sym2n > 0)
        sym2s = wopt->sym2param;
    if (zeron > 0)
        zeros = wopt->zeroparam;
    if (constn > 0)
        consts = wopt->constparam;
    p = custm;
    z = 0;
    for (count = 0; count < numpop; count++)
    {
        if (*p == 'm')
            z++;
        p++;
    }
    if (z > 0)
        force_sametheta (param, wopt, numpop);
    p = custm + numpop;
    z = 0;
    for (count = numpop; count < numpop * numpop; count++)
    {
        if (*p == 'm')
            z++;
        p++;
    }
    if (z > 0)
        force_samemigration (param, wopt, numpop);
    for (i = 0; i < zeron; i++)
    {
        param[zeros[i]] = -30.; //-MYREAL_MAX;
    }

    for (i = 0; i < constn; i++)
    {
        if (consts[i] < numpop)
            param[consts[i]] = LOG (wopt->thetag[consts[i]]);
        if (consts[i] >= numpop)
        {
            if (consts[i] >= numpop * numpop)
                param[consts[i]] = LOG (wopt->alphavalue);
            else
            {
                m2mm (consts[i], numpop, &from, &to);
                param[consts[i]] = LOG (wopt->mg[to*(numpop-1)+from]);
            }
        }
    }

    for (i = 0; i < symn; i++)
    {
        if(which_profile == syms[i][0])
            param[syms[i][1]] = param[syms[i][0]];
        else
        {
            if(which_profile == syms[i][1])
                param[syms[i][0]] = param[syms[i][1]];
            else
            {
                mm = (EXP (param[syms[i][0]]) + EXP (param[syms[i][1]])) * HALF;
                param[syms[i][0]] = param[syms[i][1]] = LOG (mm);
            }
        }
    }
    for (i = 0; i < sym2n; i++)
    {
        if(which_profile == sym2s[i][0])
            mm = EXP (param[sym2s[i][0]] + param[sym2s[i][2]]);
        else
        {
            if(which_profile == sym2s[i][1])
                mm = EXP (param[sym2s[i][0]] + param[sym2s[i][2]]);
            else
            {
                mm = (EXP (param[sym2s[i][0]] + param[sym2s[i][2]]) +
                      EXP (param[sym2s[i][1]] + param[sym2s[i][3]])) * HALF;
            }
        }
        param[sym2s[i][0]] = (lmm = LOG (mm)) - param[sym2s[i][2]];
        param[sym2s[i][1]] = lmm - param[sym2s[i][3]];
    }
}


/// calculates the gradient of the likelihood function
void
gradient (MYREAL *d, nr_fmt * nr, long locus)
{
    boolean found=FALSE;
    MYREAL tk1,m1,th1,l1;
    //MYREAL nm1,th2;
    //MYREAL nm2, l2, m2, tk2;
    long z, r;
    long other , otherpop;
    long g, i, ii, pop, msta, msto;
    long numpop = nr->numpop;
    long nn = nr->partsize - nr->profilenum;
    MYREAL expapg, *thetas;
    //MYREAL *m;
    MYREAL *geo = nr->world->data->geo;
    tarchive_fmt *tl;
    MYREAL tet;
    MYREAL rtet;
    MYREAL *normalgrad;
    MYREAL mu_rate = nr->world->options->mu_rates[locus];
    MYREAL inv_mu_rate = 1. / mu_rate;
	
    normalgrad = (MYREAL *) mymalloc(nr->numpop2 * sizeof(MYREAL));
    memset (d, 0, sizeof (MYREAL) * nr->numpop2);
    thetas = nr->param;
  //xcode  m = nr->param + numpop;
    for (r = nr->repstart; r < nr->repstop; r++)
    {
        for (pop = 0; pop < nr->numpop2; pop++)
        {
            normalgrad[pop] = (*normal_func_gradient)(thetas[pop], nr->world->atl[r][locus].param0[pop]);
        }
        for (g = 0; g < nr->world->atl[r][locus].T; g++)
        {
            if (nr->apg[r][locus][g] < -40.)
                continue;
            tl = nr->world->atl[r][locus].tl;
            expapg = EXP (nr->apg[r][locus][g]);
            for (pop = 0; pop < numpop; pop++)
            {
                tet = thetas[pop];
                // the following lines are doing this derivative, but have reduced the division operation (speed)
                // nr->parts[pop] = -tl[g].p[pop] / tet + tl[g].kt[pop] / (tet * tet * mu_rate) - normalgrad[pop];
                rtet = 1./(tet * tet * mu_rate);
                nr->parts[pop] = (- mu_rate * tet * (tl[g].p[pop] + normalgrad[pop] * tet)  + tl[g].kt[pop]) * rtet ;
                msta = nr->mstart[pop];
                msto = nr->mend[pop];
             //xcode   z = 0;
             //xcode   found=FALSE;
                for (i = msta; i < msto; i++)
                {
                    if (nr->param[i] > 0.)
                    {
                        if(nr->world->options->custm2[i]=='S')
                        {
                            tk1 = geo[i]* tl[g].km[pop] * inv_mu_rate;
                            l1 = tl[g].mindex[i];
                            m1 = nr->param[i];
                            th1 = nr->param[pop];
                            found =find (i, nr->profiles, nr->profilenum);
                            if (!found)
                            {
                                nr->parts[i] = ((l1 / m1) - tk1) - normalgrad[i];
                            }
                            symmetric_other(i, numpop, &other, &otherpop);
                        //xcode    th2 = nr->param[otherpop];
                           //xcode nm1 = m1 * th1;
                            // speed opt: nr->parts[pop]  += tk1 * m1/th1 - l1/th1;
                            nr->parts[pop]  += tk1 * (m1 - l1)/th1;
                        }
                        else
                            nr->parts[i] = ((tl[g].mindex[i] / (nr->param[i]))
                                            - geo[i] * tl[g].km[pop] * inv_mu_rate) - normalgrad[i];
                    }
                }
            }
            z = 0;
            for (i = 0; i < nn; i++)
            {
                ii = (nr->profilenum > 0) ? nr->indeks[z++] : i;
                d[i] += expapg * nr->parts[ii];
            }
        }
    }
    myfree(normalgrad);
}

#ifdef LONGSUM

long which_rate(MYREAL *rawtimes, long pop, MYREAL time_end)
{
    long i;
    long end = pop * 3 + 3 -1;
    i = pop * 3;
    while(i<end && time_end > rawtimes[i+1])
    {
        i++;
    }
    return i;
}

void
gradient_longsum (MYREAL *d, nr_fmt * nr, long locus)
{
    long z, r;
    long g, i, j, ii, pop;
    long numpop = nr->numpop;
    MYREAL expapg, *thetas, *m;
    MYREAL *geo = nr->world->data->geo;
    tarchive_fmt *tl;
    longsum_fmt *longsum;
    longsum_fmt *tint;
    MYREAL *treeparts;
    MYREAL tet;
    MYREAL time_end;
    MYREAL rawt;
    // long intervals[3];
    long numpop3 = nr->numpop * 3;
    long nn = nr->partsize - nr->profilenum - numpop3;
    MYREAL *rawtimes = nr->world->flucrates+numpop3;
    MYREAL *rawrates = nr->param+nn+nr->profilenum;
    // MYREAL *rawlrates = nr->lparam+nn;
    MYREAL mu_rate = nr->world->options->mu_rates[locus];
    // MYREAL * rates = mycalloc(numpop3 * 3,sizeof(MYREAL));
    // MYREAL *lrates = rates + numpop3;
    long ptime, ptimeold;
    MYREAL timeinterval;
    //interval[0]=intervals[1]=intervals[2]=0;
    FPRINTF(stdout,"rawrates=%f %f %f @ %f %f %f \n",rawrates[0],rawrates[1],rawrates[2],rawrates[3],rawrates[4],rawrates[5]);
    treeparts = mycalloc(nr->numpop2+numpop3,sizeof(MYREAL));
    memset (d, 0, sizeof (MYREAL) * nr->partsize);
    thetas = nr->param;
    m = nr->param + numpop;
    for (r = nr->repstart; r < nr->repstop; r++)
    {
        memset(nr->parts,0,sizeof(MYREAL)*nr->partsize);
        for (g = 0; g < nr->world->atl[r][locus].T; g++)
        {
            if (nr->apg[r][locus][g] < -40.)
                continue;
            memset(treeparts,0,sizeof(MYREAL)*nr->partsize);
            tl = nr->world->atl[r][locus].tl;
            expapg = EXP (nr->apg[r][locus][g]);
            longsum = nr->world->atl[r][locus].tl[g].longsum;
            time_end=0.;
            // for(pop=0; pop < numpop; pop++)
            // rates[pop]=1.;
            for(j=0; j< nr->world->atl[r][locus].tl[g].longsumlen; j++)
            {
                tint = &longsum[j];
                time_end += tint->interval;
                switch(tint->eventtype)
                {
                case 'c':
                    pop = tint->to;
                    tet = thetas[pop];
                    treeparts[pop] -= 1. / tet  ;

                    ptime = which_rate(rawtimes,pop, time_end);
                    treeparts[ptime+nn] -= 1./rawrates[ptime];
                    break;
                case 'm':
                    treeparts[tint->fromto] += 1./(nr->param[tint->fromto]);
                    break;
                }
                for (pop = 0; pop < nr->numpop; pop++)
                {
                    ptime = which_rate(rawtimes,pop, time_end);
                    //     intervals[ptime] += 1;
                    timeinterval = tint->interval;
                    tet = thetas[pop];
                    rawt = rawtimes[ptime];
                    if(((time_end - tint->interval) < rawt) && (time_end > rawt))
                    {
                        ptimeold = which_rate(rawtimes, pop, time_end - tint->interval);
                        //      intervals[ptimeold] += 1;
                        if(ptime-ptimeold>1)
                        {
                            timeinterval =rawtimes[ptimeold+1]-(time_end - tint->interval);
                            treeparts[ptimeold+nn] += timeinterval * tint->lineages2[pop] / (tet * mu_rate * rawrates[ptimeold] * rawrates[ptimeold]);
                            treeparts[pop] += timeinterval * tint->lineages2[pop] / (tet * tet * mu_rate * rawrates[ptimeold]);

                            timeinterval = rawtimes[ptime] - rawtimes[ptimeold+1];
                            treeparts[ptimeold+1+nn] += timeinterval * tint->lineages2[pop] / (tet * mu_rate * rawrates[ptimeold+1] * rawrates[ptimeold+1]);
                            treeparts[pop] += timeinterval * tint->lineages2[pop] / (tet * tet * mu_rate * rawrates[ptimeold+1]);

                            timeinterval = time_end - rawtimes[ptime];
                            treeparts[ptime+nn] += timeinterval * tint->lineages2[pop] / (tet * mu_rate * rawrates[ptime] * rawrates[ptime]);
                            treeparts[pop] += timeinterval * tint->lineages2[pop] / (tet * tet * mu_rate * rawrates[ptime]);
                        }
                        else
                        {
                            timeinterval =rawtimes[ptime]-(time_end - tint->interval);
                            treeparts[ptimeold+nn] += timeinterval * tint->lineages2[pop] / (tet * mu_rate * rawrates[ptimeold] * rawrates[ptimeold]);
                            treeparts[pop] += timeinterval * tint->lineages2[pop] / (tet * tet * mu_rate * rawrates[ptimeold]);
                            timeinterval = time_end - rawtimes[ptime];
                            treeparts[ptime+nn] += timeinterval * tint->lineages2[pop] / (tet * mu_rate * rawrates[ptime] * rawrates[ptime]);
                            treeparts[pop] += timeinterval * tint->lineages2[pop] / (tet * tet * mu_rate * rawrates[ptime]);
                        }
                    }
                    else
                    {
                        treeparts[ptime+nn] += timeinterval * tint->lineages2[pop] / (tet * mu_rate * rawrates[ptime] * rawrates[ptime]);
                        treeparts[pop] += timeinterval * tint->lineages2[pop] / (tet * tet * mu_rate * rawrates[ptime]);
                    }
                }
                for (pop = nr->numpop; pop < nr->numpop2; pop++)
                {
                    treeparts[pop] -= tint->interval * geo[pop] * tint->lineages[(long)((pop - numpop)/(numpop-1))] / mu_rate;
                }
            }
            for (pop = 0; pop < nr->numpop2; pop++)
            {
                nr->parts[pop] = treeparts[pop] - (*normal_func_gradient)(nr->param[pop],nr->world->atl[r][locus].param0[pop]);
            }
            for (pop = nr->numpop2; pop < nr->partsize; pop++)
            {
                nr->parts[pop] = treeparts[pop];
            }
            //   for (pop = 0; pop < 3; pop++)
            //    {
            //  FPRINTF(stdout,"%li ", intervals[pop]);
            //  }
            // FPRINTF(stdout,"\n");

            z = 0;
            for (i = 0; i < nn; i++)
            {
                ii = (nr->profilenum > 0) ? nr->indeks[z++] : i;
                d[i] += expapg * nr->parts[ii];
            }
            for (i = nn; i < nr->partsize; i++)
            {
                if((i-nn) % 3 == 0) // each population has 3 time values and the first one is always 1
                    d[i] = 0. ;
                else
                    d[i] += expapg * nr->parts[i];//da rates
            }
        }
    }
    myfree(treeparts);
}
#endif /*LONGSUM*/
/* calculates the derivatives for the different
migration models
*/
void
force_symmetric_d (MYREAL *gxv, long model, nr_fmt * nr, long nn)
{
    switch (model)
    {
    case MATRIX:
        break;
    case MATRIX_ARBITRARY:
        check_symmetric_d (gxv, nr, nn);
        break;
    case MATRIX_SAMETHETA:
        calc_symmetric_d (gxv, nr, nn, 0, nr->numpop);
        break;
    case ISLAND:
        calc_symmetric_d (gxv, nr, nn, 0, nr->numpop);
        calc_symmetric_d (gxv, nr, nn, nr->numpop, nr->numpop2);
        break;
    case ISLAND_VARTHETA:
        calc_symmetric_d (gxv, nr, nn, nr->numpop, nr->numpop2);
        break;
    }
}

void
calc_same_d (MYREAL *grad, nr_fmt * nr, long nn, long start, long stop)
{
    long i, ii;
    long z = 0;
    MYREAL dt = 0;
    long nsum = 0;
    char *custm = nr->world->options->custm2;
    if (nr->profilenum == 0)
    {
        for (i = start; i < stop; i++)
        {
            if (custm[i] == 'm')
            {
                dt += grad[i];
                nsum++;
            }
        }
        dt /= nsum;
     //xcode   z = 0;
        for (i = start; i < stop; i++)
        {
            if (custm[i] == 'm')
                grad[i] = dt;
        }
    }
    else
    {
        for (i = 0; i < nn; i++)
        {
            ii = (nr->profilenum > 0) ? nr->indeks[z++] : i;
            if (ii >= start && ii < stop && custm[ii] == 'm')
            {
                dt += grad[i];
                nsum++;
            }
        }
        dt /= nsum;
        z = 0;
        for (i = 0; i < nn; i++)
        {
            ii = (nr->profilenum > 0) ? nr->indeks[z++] : i;
            if (ii >= start && ii < stop && custm[ii] == 'm')
                grad[i] = dt;
        }
    }
}

void
calc_symmetric_d (MYREAL *grad, nr_fmt * nr, long nn, long start, long stop)
{
    long i, ii;
    long z = 0;
    MYREAL dt = 0;
    if (nr->profilenum == 0)
    {
        for (i = start; i < stop; i++)
        {
            dt += grad[i];
        }
        dt /= (stop - start);
     //xcode   z = 0;
        for (i = start; i < stop; i++)
        {
            grad[i] = dt;
        }
    }
    else
    {
        for (i = 0; i < nn; i++)
        {
            ii = (nr->profilenum > 0) ? nr->indeks[z++] : i;
            if (ii >= start && ii < stop)
                dt += grad[i];
        }
        dt /= (stop - start);
        z = 0;
        for (i = 0; i < nn; i++)
        {
            ii = (nr->profilenum > 0) ? nr->indeks[z++] : i;
            if (ii >= start && ii < stop)
                grad[i] = dt;
        }
    }
}

void
check_symmetric_d (MYREAL *grad, nr_fmt * nr, long nn)
{
    char *custm = NULL;
    //long zeron = 0;
    //long constn = 0;
    long symn = 0;
    long sym2n = 0;
    //long mmn = 0;
    long numpop = 0;
    /*static boolean done = FALSE; */
    twin_fmt *syms = NULL;
    quad_fmt *sym2s = NULL;
   //xcode  long *zeros = NULL;
   //xcode  long *consts = NULL;
   //xcode  long *mms = NULL;
    long count;
    long z = 0, zz = 0, i;
    char *p;
    MYREAL mm, sq1, sq2;


    //xcode custm = nr->world->options->custm2;
    //xcode zeron = nr->world->options->zeron;
    //xcode constn = nr->world->options->constn;
    symn = nr->world->options->symn;
    sym2n = nr->world->options->sym2n;
    //xcode mmn = nr->world->options->mmn;
    if (symn > 0)
        syms = nr->world->options->symparam;
    if (sym2n > 0)
        sym2s = nr->world->options->sym2param;
    //xcode if (zeron > 0)
    //xcode    zeros = nr->world->options->zeroparam;
    //xcode if (constn > 0)
    //xcode    consts = nr->world->options->constparam;
    //xcode if (mmn > 0)
    //xcode    mms = nr->world->options->mmparam;
    //xcode p = custm;
    custm = nr->world->options->custm2;
    numpop = nr->numpop;
    p = custm;
    z = 0;
    for (count = 0; count < numpop; count++)
    {
        if (*p == 'm')
            z++;
        p++;
    }
    if (z > 0)
        calc_same_d (grad, nr, nn, 0, nr->numpop);
    p = custm + nr->numpop;
    z = 0;
    for (count = numpop; count < numpop * numpop; count++)
    {
        if (*p == 'm')
            z++;
        p++;
    }
    if (z > 0)
        calc_same_d (grad, nr, nn, nr->numpop, nr->numpop2);
    /* if there are symmetric migration rates M go here*/
    for (i = 0; i < symn; i++)
    {
        if (nr->profilenum > 0)
        {
            z = 0;
            while (nr->indeks[z] != syms[i][0] && z++ < nr->partsize)
                ;
            zz = 0;
            while (nr->indeks[zz] != syms[i][1] && zz++ < nr->partsize)
                ;
            if (z < nr->partsize)
                mm = grad[z];
            else
                mm = 0;
            if (zz < nr->partsize)
                mm += grad[zz];
            //else
            //  mm += 0;
            grad[z] = grad[zz] = mm;
        }
        else
        {
            mm = grad[syms[i][0]] + grad[syms[i][1]];
            grad[syms[i][0]] = grad[syms[i][1]] = mm;
        }
    }
    /* if there are symmetric migration rates 4Nm */
    if (nr->profilenum > 0)
    {
        for (i = 0; i < sym2n; i++)
        {
            sq1 = nr->param[sym2s[i][2]] ;// * nr->param[sym2s[i][2]]);
            sq2 = nr->param[sym2s[i][3]];// * nr->param[sym2s[i][3]]);
            z = 0;
            while (nr->indeks[z] != sym2s[i][0] && z++ < nr->partsize)
                ;
            zz = 0;
            while (nr->indeks[zz] != sym2s[i][1] && zz++ < nr->partsize)
                ;
            if (z < nr->partsize)
                mm = grad[z] / sq1;
            else
                mm = 0;
            if (zz < nr->partsize)
                mm += grad[zz] / sq2;
            //else
            //  mm += 0;
            //  mm /= 2.;
            grad[z] = mm * sq1;
            grad[zz] = mm * sq2;
        }
    }
    else
    {
        for (i = 0; i < sym2n; i++)
        {
            sq1 = nr->param[sym2s[i][2]];// * nr->param[sym2s[i][2]];
            sq2 = nr->param[sym2s[i][3]];// * nr->param[sym2s[i][3]];
            mm = (grad[sym2s[i][0]] / sq1 + grad[sym2s[i][1]] / sq2);
            grad[sym2s[i][0]] = mm * sq1;
            grad[sym2s[i][1]] = mm * sq2;
        }
    }
}

void
symmetric_other (long i, long numpop, long *other, long *otherpop)
{
    long pop;
    m2mm (i, numpop, otherpop, &pop);
    *other = mm2m (pop, *otherpop, numpop);
}


/* derivatives to log derivatives */
void
grad2loggrad (MYREAL *param, long *indeks, MYREAL *d, long nn,
              long profilenum)
{
    long i, ii, z = 0;
    for (i = 0; i < nn; i++)
    {
        ii = (profilenum > 0) ? indeks[z++] : i;
        d[i] = -param[ii] * d[i];
        /*the minus is for minimizing  the function -L
         instead of maximizing L
         */
    }
}

MYINLINE MYREAL
probG (MYREAL *param, MYREAL *lparam, tarchive_fmt * tl, nr_fmt * nr,
       long locus)
{
    const int numpop = (int) nr->numpop;
    int i, j, msta, msto;
    MYREAL result = 0., sm;
    MYREAL *geo = nr->world->data->geo;
    MYREAL *lgeo = nr->world->data->lgeo;
    MYREAL mu_rate = nr->world->options->mu_rates[locus];
    MYREAL inv_mu_rate = 1./ mu_rate;
    MYREAL lmu_rate = nr->world->options->lmu_rates[locus];
    for (i = 0; i < numpop; i++)
    {
        if (lparam[i] <= -MYREAL_MAX)
            return MYREAL_MAX;
        result += tl->p[i] * (LOG2 - (lparam[i] + lmu_rate));
        msta = nr->mstart[i];
        msto = nr->mend[i];
        sm = 0.0;
        for (j = msta; j < msto; j++)
        {
            if (param[j] > 0.0)
            {
                result += tl->mindex[j] * (lgeo[j] + lparam[j] - lmu_rate);
                sm += geo[j] * param[j] * inv_mu_rate;
            }
        }
        result += -tl->km[i] * sm - tl->kt[i] * (1./param[i]) * inv_mu_rate;
    }
    return result;
}

MYINLINE  MYREAL
probG4 (MYREAL *fullparam, MYREAL *data, int stop)
{
    int i;
    //    int i1;
    //    int i2;
    //int stop = (int) (numpop2 + numpop + numpop);
    //MYREAL r1;
    //MYREAL r2;
    //MYREAL r3;
    MYREAL result = 0.;
    // fullparam is a linearized
    //     log(param) [size: numpop * numpop] [Log(2/theta), Log(M)] for example 3pop:  3 log(2/theta), 6 log(M)
    //     param [thetas are 1/param] [size: numpop * numpop] 3
    //     sum(migparam) [size: numpop]    3
    //    total of 12
    // data is linearized
    //     p      [size: numpop]  3
    //     mindex [sizes: numpop * (numpop -1)] 6
    //     kt     [size: numpop]  3
    //     km     [size: numpop] 3
    // total 12
    //
    //  calculation: p * log(2/theta) + mindex * log(M)
    // - kt * 1/theta - km * sm
    for (i = 0; i < stop; i++)
    {
        result += fullparam[i] * data[i];
    }
    // stop is at least 3:
    // population in dataset  stop
    //            1            3
    //            2            8
    //            3           15
    //            4           24
    //for (i = 0, i1= 1, i2 = 2; i2 < stop ; i++, i1++, i2++)
    //{
    //  r1 = fullparam[i] * data[i];
    //  r2 = fullparam[i1] * data[i1];
    //  r3 = fullparam[i2] * data[i2];
    //  result += r1 + r2 + r3;
    //}

    return result;
}

#ifdef LONGSUM
MYREAL
probG_longsum (MYREAL *param, MYREAL *lparam, tarchive_fmt * tl, nr_fmt * nr,
               long locus)
{
    const int numpop = (int) nr->numpop;
    long interval, i,  pop, fromto;
    MYREAL result = 0.;
    MYREAL *geo = nr->world->data->geo;
    MYREAL *lgeo = nr->world->data->lgeo;
    MYREAL mu_rate = nr->world->options->mu_rates[locus];
    MYREAL lmu_rate = nr->world->options->lmu_rates[locus];
    longsum_fmt *tint;
    MYREAL *sm;

    sm = (MYREAL *) mycalloc(nr->numpop, sizeof(MYREAL));
    for(pop=0;pop<nr->numpop; pop++)
    {
        for (i = nr->world->mstart[pop]; i < nr->world->mend[pop]; i++)
            sm[pop] += geo[i] * param[i] / mu_rate;
    }

    for(interval=0; interval < tl->longsumlen; interval++)
    {
        tint = &tl->longsum[interval];
        switch(tint->eventtype)
        {
        case 'c':
            result += (LOG2 - (lparam[tint->to] + lmu_rate));
            break;
        case 'm':
            if (param[tint->fromto] > 0.0)
            {
                fromto = tint->fromto;
                result += (lgeo[fromto] + lparam[fromto] - lmu_rate);
            }
            break;
        }
        for (i = 0; i < numpop; i++)
        {
            pop = tint->to;
            result -= tint->interval*tint->lineages2[pop] / (param[pop] * mu_rate);
            result -= tint->interval * tint->lineages[tint->to] * sm[pop];
        }
    }
    myfree(sm);
    return result;
}


MYREAL
probG4_longsum (MYREAL *fullparam, longsum_fmt *data, long longsumlen, long numpop, long numpop2)
{
    int interval, i;
    MYREAL result = 0.;
    MYREAL timeinterval;
    MYREAL rawt;
    long treelen=longsumlen;
    longsum_fmt *tint;
    MYREAL time_end;
    MYREAL *sm = fullparam + numpop2 + numpop;
    // we allow three rates per population [uhh too many parameters]
    MYREAL *rawrates = fullparam + numpop2 + numpop + numpop;
    MYREAL *rawlrates = rawrates + 3*numpop;
    MYREAL *rawtimes = rawrates + 6*numpop;
    //MYREAL *rates;
    //MYREAL *lrates;
    long ptime, ptimeold;
    // MYREAL *rtimes;
    // long numpop3 = 3 * numpop;
    //rates = mycalloc( numpop * 6, sizeof(MYREAL));
    // for(pop=0; pop < numpop3; pop++)
    // rates[pop] = 1.0;
    // lrates = rates + numpop3;
    //  rtimes = mycalloc( numpop3, sizeof(MYREAL));
    // fullparam is a linearized
    //     log(param) [size: numpop * numpop] [Log(2/theta), Log(M)]
    //     param [thetas are 1/param] [size: numpop ]
    //     sum(migparam) [size: numpop]
    // data is linearized
    //     p      [size: numpop]
    //     mindex [sizes: numpop * (numpop -1)]
    //     kt     [size: numpop]
    //     km     [size: numpop]
    //
    //  calculation: p * log(2/theta) + mindex * log(M)
    //  for(pop=0; pop < numpop3; pop += 3)
    //   rates[pop] = rawrates[pop];
    time_end = 0.;
    for(interval=0;interval<treelen;interval++)
    {
        tint = &data[interval];
        time_end += tint->interval;
        ptime = which_rate(rawtimes, tint->to, time_end);
        switch(tint->eventtype)
        {
        case 'c':
            result +=  fullparam[tint->to] - rawlrates[ptime] ;
            break;
        case 'm':
            if (fullparam[tint->fromto] > -MYREAL_MAX)
                result += fullparam[tint->fromto];
            break;
        }
        // - kt * 1/theta - km * sm
        for (i = 0; i < numpop; i++)
        {
            timeinterval = tint->interval;
            ptime = which_rate(rawtimes, i, time_end);
            rawt = rawtimes[ptime];
            if(((time_end - tint->interval) < rawt) && (time_end > rawt))
            {
                ptimeold = which_rate(rawtimes, i, time_end - tint->interval);
                if(ptime-ptimeold>1)
                {
                    timeinterval =  rawtimes[ptimeold+1] - (time_end - tint->interval);
                    result += timeinterval * (tint->lineages2[i] / rawrates[ptimeold] * fullparam[i+numpop2] + tint->lineages[i] * sm[i]);
                    timeinterval =  rawtimes[ptime] - rawtimes[ptimeold+1];
                    result += timeinterval * (tint->lineages2[i] / rawrates[ptimeold+1] * fullparam[i+numpop2] + tint->lineages[i] * sm[i]);
                    timeinterval = time_end - rawtimes[ptime];
                    result +=  timeinterval * (tint->lineages2[i] / rawrates[ptime] * fullparam[i+numpop2] + tint->lineages[i] * sm[i]);
                }
                else
                {
                    timeinterval =  rawtimes[ptime] - (time_end - tint->interval);
                    result += timeinterval * (tint->lineages2[i] / rawrates[ptimeold] * fullparam[i+numpop2] + tint->lineages[i] * sm[i]);
                    timeinterval = time_end - rawtimes[ptime];
                    result +=  timeinterval * (tint->lineages2[i] / rawrates[ptime] * fullparam[i+numpop2] + tint->lineages[i] * sm[i]);
                    //     FPRINTF(stderr,"on border timeend %f, interval %f, rawrate %f\n",time_end,tint->interval, rawrates[ptime]);
                }
            }
            else
            {
                result +=  timeinterval * (tint->lineages2[i] / rawrates[ptime] * fullparam[i+numpop2] + tint->lineages[i] * sm[i]);
            }
            //   FPRINTF(stdout,"p=%i t=%f, u=%f ptime=%li rt=%f res=%f\n", i, time_end, timeinterval, ptime, rawtimes[ptime], result);
        }
    }
    // FPRINTF(stdout,"\n");
    //myfree(rates);
    //  myfree(rtimes);
    return result;
}
#endif /*LONGSUM*/

void
probG3 (MYREAL *apg, MYREAL *apg0r, timearchive_fmt * tyme,
        long numpop, MYREAL *apgmax,
        MYREAL *param, MYREAL *lparam, MYREAL *sm)
{
    long g, g1, g2, g3, i, j;
    long msta, me;
    MYREAL temp, temp1, temp2;
    MYREAL result1, result2, result3, result4;
    MYREAL tparam;
    tarchive_fmt *tl1;
    tarchive_fmt *tl2;
    tarchive_fmt *tl3;
    tarchive_fmt *tl4;
    long remaind = tyme->T % 4;

    if (remaind > 0)
    {
        for (g = 0; g < remaind; g++)
        {
            tl1 = &(tyme->tl[g]);
            result1 = 0.;
            for (i = 0; i < numpop; i++)
            {
                msta = numpop + i * (numpop - 1);
                me = msta + numpop - 1;
                temp = (LOG2 - lparam[i]);
                result1 += tl1->p[i] * temp;
                for (j = msta; j < me; j++)
                {
                    if (param[j] > 0.0)
                    {
                        result1 += tl1->mindex[j] * lparam[j];
                    }
                }
                result1 += -tl1->km[i] * sm[i] - tl1->kt[i] / param[i];
            }
            apg[g] = result1 - apg0r[g];
            if (apg[g] > *apgmax)
                *apgmax = apg[g];
        }
    }
    for (g = remaind; g < tyme->T; g += 4)
    {
        g1 = g + 1;
        g2 = g1 + 1;
        g3 = g2 + 1;
        tl1 = &(tyme->tl[g]);
        tl2 = &(tyme->tl[g1]);
        tl3 = &(tyme->tl[g2]);
        tl4 = &(tyme->tl[g3]);
        result1 = 0.;
        result2 = 0.;
        result3 = 0.;
        result4 = 0.;
        for (i = 0; i < numpop; i++)
        {
            tparam = 1. / param[i];
            msta = numpop + i * (numpop - 1);
            me = msta + numpop - 1;
            temp = (LOG2 - lparam[i]);
            result1 += tl1->p[i] * temp;
            result2 += tl2->p[i] * temp;
            result3 += tl3->p[i] * temp;
            result4 += tl4->p[i] * temp;
            for (j = msta; j < me; j++)
            {
                if (param[j] > 0.0)
                {
                    result1 += tl1->mindex[j] * lparam[j];
                    result2 += tl2->mindex[j] * lparam[j];
                    result3 += tl3->mindex[j] * lparam[j];
                    result4 += tl4->mindex[j] * lparam[j];
                }
            }
            result1 += -tl1->km[i] * sm[i] - tl1->kt[i] * tparam;
            result2 += -tl2->km[i] * sm[i] - tl2->kt[i] * tparam;
            result3 += -tl3->km[i] * sm[i] - tl3->kt[i] * tparam;
            result4 += -tl4->km[i] * sm[i] - tl4->kt[i] * tparam;
        }
        apg[g] = result1 - apg0r[g];
        apg[g1] = result2 - apg0r[g1];
        apg[g2] = result3 - apg0r[g2];
        apg[g3] = result4 - apg0r[g3];
        temp1 = MAX (apg[g], apg[g1]);
        temp2 = MAX (apg[g2], apg[g3]);
        if ((temp = MAX (temp1, temp2)) > *apgmax)
            *apgmax = temp;
    }
}

MYREAL
probG2 (MYREAL *param, MYREAL *lparam, MYREAL *sm,
        MYREAL *kt, MYREAL *km, MYREAL *p, MYREAL *mindex,
        int *msta, int *me, long numpop)
{
    int i, j;
    MYREAL result = 0.;
    for (i = 0; i < numpop; i++)
    {
        if (lparam[i] <= -MYREAL_MAX)
            return MYREAL_MAX;
        //      result += p[i] * (LOG2 - lparam[i]);
        for (j = msta[i]; j < me[i]; j++)
        {
            if (param[j] > 0.0)
                result += mindex[j] * lparam[j];
        }
        result += p[i] * lparam[i] - km[i] * sm[i] - kt[i] * param[i];
        //      result += -km[i] * sm[i] - kt[i] / param[i];
    }
    return result;
}


boolean
is_singular (MYREAL **dd, long n)
{
    long i, j;
    MYREAL temp;
    boolean singular = FALSE;
    for (i = 0; i < n; i++)
    {
        temp = 0.0;
        for (j = 0; j < n; j++)
        {
            temp += dd[i][j];
        }
        if (!(temp > 0.0 || temp < 0.0))
        {
            singular = TRUE;
            break;
        }
    }
    for (i = 0; i < n; i++)
    {
        temp = 0.0;
        for (j = 0; j < n; j++)
        {
            temp += dd[i][j];
        }
        if (!(temp < 0.0 || temp > 0.0))
        {
            singular = TRUE;
            break;
        }
    }
    return singular;
}


void
calc_cov (MYREAL **dd, MYREAL *d, MYREAL *param, long n)
{
    long i, j;
    if (!is_singular (dd, n))
        invert_matrix (dd, n);
    else
    {
        reset_hess (dd, n);
        return;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            dd[i][j] /= (param[i] * param[j]);
            dd[j][i] = dd[i][j];
        }
        dd[i][i] = (dd[i][i] - param[i] * d[i]) / (param[i] * param[i]);
    }
    if (!is_singular (dd, n))
        invert_matrix (dd, n);
    else
        reset_hess (dd, n);
}


void
print_contribution (nr_fmt * nr, timearchive_fmt ** atl, long G)
{
    long g, i, r;
    MYREAL temp = 0, temp2 = 0, maxtemp;
    FILE *mixfile = nr->world->options->mixfile;
    long events = 0;
    long contribution[11];
    long locus = nr->world->locus;
    tarchive_fmt *tyme;
    for (g = 0; g < 11; g++)
        contribution[g] = 0;
    if (nr->world->options->mixplot)
        FPRINTF (mixfile, "\n\n");
    maxtemp = -MYREAL_MAX;
    for (r = nr->repstart; r < nr->repstop; r++)
    {
        for (g = 0; g < G; g++)
        {
            tyme = atl[r][nr->world->locus].tl;
            temp = nr->apg[r][locus][g] + nr->apg0[r][locus][g];
            temp2 = temp + nr->world->likelihood[g];
            if (temp2 > maxtemp)
                maxtemp = temp2;
            if (nr->world->options->mixplot)
            {
                events = 0;
                for (i = nr->world->numpop; i < nr->world->numpop2; i++)
                    events += (long) tyme[g].mindex[i];
                for (i = 0; i < tyme[g].copies; i++)
                    FPRINTF (mixfile, "%li %f %f ", events, temp,
                             nr->world->likelihood[g]);
            }
            temp2 -= maxtemp;
            if (temp2 > -20)
            {
                contribution[9 - (long) (fabs (temp2) / 2)] +=
                    (g > 0) ? tyme[g].copies : tyme[g].copies - 1;
            }
            contribution[10] += (g > 0) ? tyme[g].copies : tyme[g].copies - 1;
        }
    }
    if (nr->world->options->progress)
    {
        print_menu_contribution (stdout, contribution);
        if (nr->world->options->writelog)
            print_menu_contribution (nr->world->options->logfile, contribution);
    }
}

void
print_menu_contribution (FILE * file, long contribution[])
{
    long g;
    FPRINTF (file, "           log(P(g|Param) * P(data|g)\n");
    FPRINTF (file, "                            -20 to ");
    for (g = -18; g <= 0; g += 2)
    {
        FPRINTF (file, "%4li ", g);
    }
    FPRINTF (file, "  All\n");
    FPRINTF (file, "           Counts                  ");
    for (g = 0; g < 10; g++)
    {
        FPRINTF (file, "%4li ", contribution[g]);
    }
    FPRINTF (file, "%5li\n", contribution[10]);
}

/* calculates log(parameters)*/
void
log_param0 (MYREAL *param, MYREAL *lparam, long nn)
{
    long i;
    for (i = 0; i < nn; i++)
    {
        if (param[i] > 0)
            lparam[i] = LOG (param[i]);
        else
            lparam[i] = -30; //-MYREAL_MAX;
    }
}


void
copies2lcopies (timearchive_fmt * atl)
{
    long g;
    if (atl->tl[0].copies > 1)
        atl->tl[0].lcopies = ln_copies (atl->tl[0].copies - 1);
    else
        atl->tl[0].lcopies = -MYREAL_MAX;
    for (g = 1; g < atl->T; g++)
    {
        atl->tl[g].lcopies = ln_copies (atl->tl[g].copies);
    }
}

void
create_apg0 (MYREAL *apg0, nr_fmt * nr, timearchive_fmt * tyme, long locus)
{
    long g;
    long copies;
    /* Prob(G|Param0) */
    //  FPRINTF(stdout,"first 4 theta0: %f %f %f %f \n",tyme->param0[0], tyme->param0[1], tyme->param0[2], tyme->param0[3]);
    for (g = 0; g < tyme->T; g++)
    {
        if (g > 0)
            copies = tyme->tl[g].copies;
        else
            copies = tyme->tl[g].copies - 1;
        if (copies == 0)
            apg0[g] = -MYREAL_MAX;
        else
        {
#ifdef LONGSUM
            apg0[g] =probG_longsum (tyme->param0, tyme->lparam0, &tyme->tl[g], nr, locus);
#else /*LONGSUM*/
	    //#ifdef INTEGRATEDLIKE
	    //apg0[g] = 0;
	    //#else
            apg0[g] =  probG (tyme->param0, tyme->lparam0, &tyme->tl[g], nr, locus);
	    //#endif /*INTEGRATEDLIKE*/
#endif /*LONGSUM*/

        }
    }
    //   FPRINTF(stdout,"first 4 apg0: %f %f %f %f \n",apg0[0], apg0[1], apg0[2], apg0[3]);

}



MYREAL
sum_mig (MYREAL *param, long msta, long msto)
{
    long i;
    MYREAL result = 0.0;
    for (i = msta; i < msto; i++)
    {
        result += param[i];
    }
    return result;
}

// penalizes if theta is too far from theta0
// std22 = std*std*2
MYREAL normal_func_ok(MYREAL *param, MYREAL *param0, long numpop2)
{
    long i;
    MYREAL p0;
    MYREAL result=0.0;
    MYREAL diff;
    for(i=0; i<numpop2;i++)
    {
        if(param0[i]>0.0)
        {
            diff = param[i]-(p0=param0[i]);
            result += -(diff*diff/(2. * p0 *p0)) - 0.91893853320467274178 - log(p0);
        }
    }
    return result;
}

MYREAL normal_func_no(MYREAL *param, MYREAL *param0, long numpop2)
{
    return 0.0;
}

MYREAL normal_func_gradient_ok(MYREAL p1, MYREAL p0)
{

    return  (p1-p0)/(p0*p0);//STD2;
}

MYREAL normal_func_gradient_no(MYREAL p1, MYREAL p0)
{
    return 0.0;
}

void unset_penalizer_function(boolean inprofiles)
{
    if(inprofiles)
    {
        normal_func = (MYREAL (*)(MYREAL *, MYREAL *, long)) normal_func_no;
        normal_func_gradient = (MYREAL (*)(MYREAL , MYREAL)) normal_func_gradient_no;
    }
    else
      {   ///debug was originaly OK and not NO
        normal_func = (MYREAL (*)(MYREAL *, MYREAL *,  long)) normal_func_no;
        normal_func_gradient = (MYREAL (*)(MYREAL, MYREAL)) normal_func_gradient_no;
    }
}

/*
 ALTIVEC garbage dump
 
 MYREAL
 calc_locus_like (nr_fmt * nr, MYREAL *param, MYREAL *lparam, long locus)
 {
  long g, r, j, copies,  r4;
  const long numpop = nr->world->numpop;
  const long numpop2 = nr->world->numpop2;
  int msta, msto;
  MYREAL gsum = 0;
  MYREAL ***apg0;
  MYREAL apgmax = -MYREAL_MAX;
  //  MYREAL ***apgall = nr->apg;
  MYREAL *apg;
  MYREAL *apg0r;
  timearchive_fmt tyme;
  tarchive_fmt *tl;
  MYREAL *geo = nr->world->data->geo;
  MYREAL *lgeo = nr->world->data->lgeo;
  FloatVec *floatparam;
  FloatVec *floatlparam;
  FloatVec *floatgeo;
  FloatVec *floatlgeo;
  FloatVec *locallparam;
  FloatVec *localparam;
  FloatVec *sm;
  float mu_rate = (float) nr->world->options->mu_rates[locus];
  float invmu_rate = 1./mu_rate;
  float lmu_rate = (float) nr->world->options->lmu_rates[locus];
  MYREAL normaldev;
  long sizepops = numpop - (numpop % 4) + 4;
  long sizeall = numpop2 - (numpop2 % 4) + 4;
  float log2 = LOG2;
  vector float vlog2;
  vector float vlmu;
  vector float vmu;
  vector float vimu;
  vector float vtmp1, vtmp2, vtmp3;
  vector float minuszero = vec_zero();
  long size = sizeall + sizepops + sizepops;
  floatparam = (FloatVec *) mycalloc (size, sizeof (FloatVec));
  floatlparam = (FloatVec *) mycalloc (size, sizeof (FloatVec));
  load_double_floatvec(floatparam, param,numpop2);
  load_double_floatvec(floatlparam, lparam,numpop2);
  floatgeo = (FloatVec *) mycalloc (size, sizeof (FloatVec));
  floatlgeo = (FloatVec *) mycalloc (size, sizeof (FloatVec));
  load_double_floatvec(floatgeo, geo,numpop2);
  load_double_floatvec(floatlgeo, lgeo,numpop2);
  locallparam = (FloatVec *) mycalloc (size, sizeof (FloatVec));
  localparam = locallparam + sizeall;
  sm = localparam + sizepops;
  memcpy (localparam, floatparam, sizeof (FloatVec) * numpop);
  memcpy (locallparam, floatlparam, sizeof (FloatVec) * numpop2);
 
  // generate vlog2={LOG2, LOG2, LOG2, LOG2}
  vlog2 = load_float_splat(&log2);
  // generate vlmu={lmu_rate, ....}
  vlmu = load_float_splat(&lmu_rate);
  // generate vmu={mu_rate, ....}
  vmu = load_float_splat(&mu_rate);
  // generate inverse of vmu vimu={1/mu_rate, ....}
  vimu = load_float_splat(&invmu_rate);
 
  nr->PGC[locus] = 0.0;
  apg0 = nr->apg0;
 
  for (r = 0; r < sizepops/4; r++)
    {
   // y = log2 - logparam + logmurate
   locallparam[r].v = vec_sub(vlog2, vec_add(locallparam[r].v,vlmu));
   // y = 1/(logparam * murate);
   vtmp1 = vec_madd(localparam[r].v,vmu, minuszero);
   vtmp2 = vec_re(vtmp1);
   // Newton-Raphson improvement
   vtmp3 = vec_madd(vtmp2, vec_nmsub(vtmp2, vtmp1, vec_float_one()),vtmp2);
   // changing sign for later streamlining
   // is there a better way to do that?
   localparam[r].v = vec_sub(minuszero,vtmp3);
   msta = nr->mstart[r];
   msto = nr->mend[r];
   zz = 0;
   for (j = msta; j < msto; j++)
     {
    sm[zz].f[ -= geo[j] * param[j] / mu_rate; //minus, so that we can loop over all in probG4
     locallparam[j] += lgeo[j] - lmu_rate;
    }
 
   }
// for(j=0; j< sizeall /4 ; j++)
//  {
//sm[r] -= geo[j] * param[j] / mu_rate; //minus, so that we can loop over all in probG4
//locallparam[j] += lgeo[j] - lmu_rate;
// not correctly done yet
//  floatparam[j].v = vec_madd(floatgeo[j].v, vec_madd(floatparam[j].v,vimu,minuszero),minuszero);
// floatlparam[j].v = vec_sub(floatlgeo[j].v, vlmu);
//}
//memcpy(locallparam+numpop, floatlparam+numpop,sizeof(float)*numpop*(numpop-1));
 
gsum = 0;
for (r = nr->repstart; r < nr->repstop; r++)
{
 apg = nr->apg[r][locus];
 apg0r = apg0[r][locus];
 tyme = nr->atl[r][locus];
 normaldev =  (*normal_func)(param, tyme.param0, numpop2);
 //first element
 tl = &(tyme.tl[0]);
 copies = tl->copies - 1;
 gsum += copies;
 
 if (copies > 0)
   {
  apg[0] =  tl->lcopies + probG4 (locallparam, tl->data, numpop,
                                  numpop2) - apg0r[0]
  + normaldev;
  ;
   }
 else
  apg[0] = -MYREAL_MAX;
 if (apgmax < apg[0])
  apgmax = apg[0];
 //other elements
 for (g = 1; g < tyme.T - tyme.T % 4; g++)
   {
  //   tl = &(tyme.tl[g]);
  gsum +=  tyme.tl[g]->copies+tyme.tl[g+1]->copies+tyme.tl[g+2]->copies+tyme.tl[g+3]->copies;
  //   apg[g] = tl->lcopies + probG4 (locallparam, tl->data, numpop, numpop2) - apg0r[g]
  //    +normaldev;
  //
  temp.v = fourdot_products(locallparam, tyme.tl[g]->data,
       tyme.tl[g+1]->data,,tyme.tl[g+2]->data,,tyme.tl[g+3]->data, minuszero);
  maxtemp.v = vec_max(temp.v, maxtemp.v);
  //
  //if (apg[g] > apgmax)
  // apgmax = apg[g];
  apg[g] = temp.f[0];
  apg[g+1] = temp.f[1];
  apg[g+2] = temp.f[2];
  apg[g+3] = temp.f[3];
  apgmax = maxtemp.f[3];
   }
}    // end replicates
for (r = nr->repstart; r < nr->repstop; r++)
{
 apg = nr->apg[r][locus];
 apg0r = apg0[r][locus];
 tyme = nr->atl[r][locus];
 // first element
 tl = &(tyme.tl[0]);
 // first element
 apg[0] -= apgmax;
 nr->PGC[locus] += EXP (apg[0]);
 // all other elements
 for (g = 1; g < tyme.T; g++)
   {
  apg[g] -= apgmax;
  if (apg[g] > -40.)
   nr->PGC[locus] += EXP (apg[g]);
   }
}    // replicates
nr->apg_max[locus] = apgmax;
nr->llike = apgmax + LOG (nr->PGC[locus]) - LOG (gsum);
myfree(locallparam);
return nr->llike;
}*/
