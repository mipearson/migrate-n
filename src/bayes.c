/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
    Bayesian   R O U T I N E S

    Peter Beerli 2003, Tallahassee
    beerli@fsu.edu

    Copyright 2003-2009 Peter Beerli, Tallahassee

    This software is distributed free of charge for non-commercial use
    and is copyrighted. Of course, we do not guarantee that the software
    works and are not responsible for any damage you may cause or have.


 $Id: bayes.c 1854 2011-03-24 18:07:30Z beerli $
    */
/*! \file bayes.c 

bayes.c contains functions to handle the Bayesian implementation of migrate
it is triggered through the options->bayes_infer=TRUE.

*/
#include <stdlib.h>
#include "definitions.h"
#include "migration.h"
#include "bayes.h"
#include "random.h"
#include "tools.h"
#include "sighandler.h"
#include "world.h"
#include "tree.h"
#include "slice.h"
#include "reporter.h"
#ifdef PRETTY
#include "pretty.h"
#endif
#ifdef MPI
#include "migrate_mpi.h"
#else
extern int myID;
#endif

// counter for bayesallfile on a per locus basis
//extern long * mdimfilecount;
extern float page_height;

extern void reprecalc_world(world_fmt *world, long that);
#ifdef MPI
extern void print_marginal_like(float *temp, long *z, world_fmt * world);
#else
extern void print_marginal_like(char *temp, long *z, world_fmt * world);
#endif
extern void  print_marginal_order(char *buf, long *bufsize, world_fmt *world);

//boolean bayes_accept (MYREAL newval, MYREAL oldval, MYREAL heat, MYREAL hastingsratio);
void bayes_print_accept(FILE * file,  world_fmt *world);
//long bayes_update (world_fmt * world);

// function pointers for the bayes machinery
MYREAL (*log_prior_ratiotheta) (MYREAL, MYREAL, bayes_fmt *, long);
MYREAL (*log_prior_ratiomig) (MYREAL, MYREAL, bayes_fmt *, long);
MYREAL (*log_prior_ratiorate) (MYREAL, MYREAL, bayes_fmt *, long);
MYREAL (*log_prior_ratioall) (MYREAL, MYREAL, bayes_fmt *, long);
MYREAL (*log_prior_theta) (world_fmt *, long);//for heating
MYREAL (*log_prior_mig) (world_fmt *, long);//for heating
MYREAL (*log_prior_rate) (world_fmt *, long);//for heating
// for scaling multiple loci correctly
// 
MYREAL (*log_prior_theta1) (world_fmt *, long, MYREAL);
MYREAL (*log_prior_mig1) (world_fmt *, long, MYREAL);
MYREAL (*log_prior_rate1) (world_fmt *, long, MYREAL);
//
MYREAL (*propose_newtheta) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL *);
MYREAL (*propose_newmig) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL *);
MYREAL (*propose_newrate) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL *);

MYREAL (*hastings_ratiotheta) (MYREAL , MYREAL, MYREAL, MYREAL, bayes_fmt *, long );
MYREAL (*hastings_ratiomig) (MYREAL , MYREAL, MYREAL, MYREAL, bayes_fmt *, long );
MYREAL (*hastings_ratiorate) (MYREAL , MYREAL, MYREAL, MYREAL, bayes_fmt *, long );

MYREAL propose_uni_newparam (MYREAL param, MYREAL mean,
                             MYREAL delta, MYREAL minparam, MYREAL maxparam, MYREAL *r);
MYREAL propose_exp_newparam (MYREAL param,MYREAL mean,
                             MYREAL delta, MYREAL minparam, MYREAL maxparam, MYREAL *r);
MYREAL propose_expb_newparam (MYREAL param,
                              MYREAL mean, MYREAL delta, MYREAL minparam, MYREAL maxparam, MYREAL *r);
MYREAL propose_mult_newparam (MYREAL param,
                              MYREAL mean, MYREAL delta, MYREAL minparam, MYREAL maxparam, MYREAL *r);
MYREAL propose_gamma_newparam (MYREAL param,MYREAL mean,
                             MYREAL delta, MYREAL minparam, MYREAL maxparam, MYREAL *r);

MYREAL      log_prior_ratio_uni  (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
MYREAL      log_prior_ratio_exp  (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
MYREAL      log_prior_ratio_wexp (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
MYREAL      log_prior_ratio_mult (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);
MYREAL      log_prior_ratio_gamma  (MYREAL newparam, MYREAL oldparam, bayes_fmt *bayes, long which);

MYREAL      log_prior_uni  (world_fmt * world, long numparam);//for heating
MYREAL      log_prior_exp  (world_fmt * world, long numparam);//for heating
MYREAL      log_prior_wexp (world_fmt * world, long numparam);//for heating
MYREAL      log_prior_mult (world_fmt * world, long numparam);//for heating
MYREAL      log_prior_gamma  (world_fmt * world, long numparam);//for heating

MYREAL      log_prior_uni1  (world_fmt * world, long numparam, MYREAL);
MYREAL      log_prior_exp1  (world_fmt * world, long numparam, MYREAL);
MYREAL      log_prior_wexp1 (world_fmt * world, long numparam, MYREAL);
MYREAL      log_prior_mult1 (world_fmt * world, long numparam, MYREAL);
MYREAL      log_prior_gamma1  (world_fmt * world, long numparam, MYREAL);

MYREAL hastings_ratio_uni(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);
MYREAL hastings_ratio_exp(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);
MYREAL hastings_ratio_expb(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);
MYREAL hastings_ratio_mult(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);
MYREAL hastings_ratio_gamma(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);


MYINLINE MYREAL probWait(long *lines,MYREAL *locallparam, long numpop, long numpop2);

void calculate_credibility_interval(world_fmt * world, long locus);

void calc_hpd_credibility(bayes_fmt *bayes,long locus, long numpop2, long numparam);

void bayes_combine_loci(world_fmt * world);
void print_locus_histogram_header(FILE *bayesfile, MYREAL *deltas, char *custm2, long numpop, long numparam, boolean usem, boolean mu, world_fmt *world);
void print_locus_histogram(FILE *bayesfile, bayes_fmt * bayes, long locus, long numparam);
void print_loci_histogram(FILE *bayesfile, bayes_fmt * bayes, long locus, long numparam);

void bayes_set_param(MYREAL *param, MYREAL newparam, long which, char *custm2, long numpop);

void adjust_bayes_min_max(world_fmt* world, MYREAL **mini, MYREAL **maxi, MYREAL **adjmini, MYREAL **adjmaxi);
void bayes_progress(world_fmt *world);
#ifdef ZNZ
void	  print_bayes_tofile(znzFile mdimfile, MYREAL *params, bayes_fmt *bayes, world_fmt *world, long numpop2, long mdiminterval, long locus, long replicate, long step);
#else
void	  print_bayes_tofile(FILE *mdimfile, MYREAL *params, bayes_fmt *bayes, world_fmt *world, long numpop2, long mdiminterval, long locus, long replicate, long step);
#endif
#ifdef DEBUG
extern void print_bf_values(world_fmt * world);
#endif

void recalc_timelist (world_fmt * world, MYREAL new_ratio , MYREAL old_ratio);


/// \brief Decide which prior distribution for the THETA parameter to use
/// Decide which prior distribution to use for THETA: the functionpointers propose_newparam_x will hold 
/// either the Exponential prior distribution or a Uniform prior distribution or .. other priors ..
/// each prior distribution has its own specific hastings ratio that will be calculated in the
/// function ptr hastings_ratio
void
which_theta_prior (int kind)
{
    switch (kind)
    {
    case UNIFORMPRIOR:
      log_prior_ratiotheta = (MYREAL (*) (MYREAL, MYREAL, bayes_fmt *, long)) log_prior_ratio_uni;
      log_prior_theta = (MYREAL (*) (world_fmt *, long)) log_prior_uni;
      log_prior_theta1 = (MYREAL (*) (world_fmt *, long, MYREAL)) log_prior_uni1;
      propose_newtheta = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL *)) propose_uni_newparam;
      hastings_ratiotheta = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long)) hastings_ratio_uni;
      break;
    case EXPPRIOR:  
      log_prior_ratiotheta = (MYREAL (*) (MYREAL,  MYREAL, bayes_fmt *, long)) log_prior_ratio_exp;
      log_prior_theta = (MYREAL (*) (world_fmt *, long)) log_prior_exp;
      log_prior_theta1 = (MYREAL (*) (world_fmt *,  long, MYREAL)) log_prior_exp1;
      propose_newtheta = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_exp_newparam;
      hastings_ratiotheta = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long )) hastings_ratio_exp;
      break;
    case WEXPPRIOR:
      log_prior_ratiotheta = (MYREAL (*) (MYREAL,  MYREAL, bayes_fmt *, long)) log_prior_ratio_wexp;
      log_prior_theta = (MYREAL (*) (world_fmt *, long)) log_prior_wexp;
      log_prior_theta1 = (MYREAL (*) (world_fmt *,  long, MYREAL)) log_prior_wexp1;
      propose_newtheta = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_expb_newparam;
      hastings_ratiotheta = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long)) hastings_ratio_expb;
      break;
    case MULTPRIOR:  
      log_prior_ratiotheta = (MYREAL (*) (MYREAL,  MYREAL,bayes_fmt *, long)) log_prior_ratio_mult;
      log_prior_theta = (MYREAL (*) (world_fmt * , long)) log_prior_mult;
      log_prior_theta1 = (MYREAL (*) (world_fmt * ,  long, MYREAL)) log_prior_mult1;
      propose_newtheta = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_mult_newparam;
      hastings_ratiotheta = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long )) hastings_ratio_mult;
      break;
    case GAMMAPRIOR:
      log_prior_ratiotheta = (MYREAL (*) (MYREAL,  MYREAL, bayes_fmt *, long)) log_prior_ratio_gamma;
      log_prior_theta1 = (MYREAL (*) (world_fmt *,  long, MYREAL)) log_prior_gamma1;
      propose_newtheta = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_gamma_newparam;
      hastings_ratiotheta = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long)) hastings_ratio_gamma;
      break;
    default:
      error("Prior distribution hookup for THETA failed");
      break;
    }
}
/// \brief Decide which prior distribution for the MIG parameter to use
/// Decide which prior distribution to use for MIG: the functionpointers propose_newparam_x will hold 
/// either the Exponential prior distribution or a Uniform prior distribution or .. other priors ..
/// each prior distribution has its own specific hastings ratio that will be calculated in the
/// function ptr hastings_ratio
void
which_mig_prior (int kind)
{
    switch (kind)
    {
    case UNIFORMPRIOR:
      log_prior_ratiomig = (MYREAL (*) (MYREAL, MYREAL, bayes_fmt *, long)) log_prior_ratio_uni;
      log_prior_mig = (MYREAL (*) (world_fmt * , long)) log_prior_uni;
      log_prior_mig1 = (MYREAL (*) (world_fmt * , long, MYREAL)) log_prior_uni1;
      propose_newmig = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL *)) propose_uni_newparam;
      hastings_ratiomig = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long)) hastings_ratio_uni;
      break;
    case EXPPRIOR:  
      log_prior_ratiomig = (MYREAL (*) (MYREAL,  MYREAL, bayes_fmt *, long)) log_prior_ratio_exp;
      log_prior_mig = (MYREAL (*) (world_fmt * , long)) log_prior_exp;
      log_prior_mig1 = (MYREAL (*) (world_fmt * , long, MYREAL)) log_prior_exp1;
      propose_newmig = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_exp_newparam;
      hastings_ratiomig = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long )) hastings_ratio_exp;
      break;
    case WEXPPRIOR:
      log_prior_ratiomig = (MYREAL (*) (MYREAL,  MYREAL, bayes_fmt *, long)) log_prior_ratio_wexp;
      log_prior_mig = (MYREAL (*) (world_fmt * , long)) log_prior_wexp;
      log_prior_mig1 = (MYREAL (*) (world_fmt * , long, MYREAL)) log_prior_wexp1;
      propose_newmig = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_expb_newparam;
      hastings_ratiomig = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long)) hastings_ratio_expb;
      break;
    case MULTPRIOR:  
      log_prior_ratiomig = (MYREAL (*) (MYREAL,  MYREAL,bayes_fmt *, long)) log_prior_ratio_mult;
      log_prior_mig = (MYREAL (*) (world_fmt * , long)) log_prior_mult;
      log_prior_mig1 = (MYREAL (*) (world_fmt * , long, MYREAL)) log_prior_mult1;
      propose_newmig = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_mult_newparam;
      hastings_ratiomig = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long )) hastings_ratio_mult;
      break;
    case GAMMAPRIOR:  
      log_prior_ratiomig = (MYREAL (*) (MYREAL,  MYREAL,bayes_fmt *, long)) log_prior_ratio_gamma;
      log_prior_mig = (MYREAL (*) (world_fmt * , long)) log_prior_gamma;
      log_prior_mig1 = (MYREAL (*) (world_fmt * , long, MYREAL)) log_prior_gamma1;
      propose_newmig = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_gamma_newparam;
      hastings_ratiomig = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long )) hastings_ratio_gamma;
      break;
    default:
      error("Prior distribution hookup for Migration parameters failed");
      break;
    }
}
/// \brief Decide which prior distribution for the mutation rate modifier parameter to use
/// Decide which prior distribution to use for mutation rate modifier: 
/// the functionpointers propose_newparam_x will hold 
/// either the Exponential prior distribution or a Uniform prior distribution or .. other priors ..
/// each prior distribution has its own specific hastings ratio that will be calculated in the
/// function ptr hastings_ratio
void
which_rate_prior (int kind)
{
    switch (kind)
    {
    case UNIFORMPRIOR:
      log_prior_ratiorate = (MYREAL (*) (MYREAL, MYREAL, bayes_fmt *, long)) log_prior_ratio_uni;
      log_prior_rate = (MYREAL (*) (world_fmt * , long)) log_prior_uni;
      log_prior_rate1 = (MYREAL (*) (world_fmt * , long, MYREAL)) log_prior_uni1;
      propose_newrate = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL *)) propose_uni_newparam;
      hastings_ratiorate = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long)) hastings_ratio_uni;
      break;
    case EXPPRIOR:  
      log_prior_ratiorate = (MYREAL (*) (MYREAL,  MYREAL, bayes_fmt *, long)) log_prior_ratio_exp;
      log_prior_rate = (MYREAL (*) (world_fmt * , long)) log_prior_exp;
      log_prior_rate1 = (MYREAL (*) (world_fmt * , long, MYREAL)) log_prior_exp1;
      propose_newrate = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_exp_newparam;
      hastings_ratiorate = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long )) hastings_ratio_exp;
      break;
    case WEXPPRIOR:
      log_prior_ratiorate = (MYREAL (*) (MYREAL,  MYREAL, bayes_fmt *, long)) log_prior_ratio_wexp;
      log_prior_rate = (MYREAL (*) (world_fmt * , long)) log_prior_wexp;
      log_prior_rate1 = (MYREAL (*) (world_fmt * , long, MYREAL)) log_prior_wexp1;
      propose_newrate = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_expb_newparam;
      hastings_ratiorate = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long)) hastings_ratio_expb;
      break;
    case MULTPRIOR:  
      log_prior_ratiorate = (MYREAL (*) (MYREAL,  MYREAL,bayes_fmt *, long)) log_prior_ratio_mult;
      log_prior_rate = (MYREAL (*) (world_fmt * , long)) log_prior_mult;
      log_prior_rate1 = (MYREAL (*) (world_fmt * , long, MYREAL)) log_prior_mult1;
      propose_newrate = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_mult_newparam;
      hastings_ratiorate = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long )) hastings_ratio_mult;
      break;
    case GAMMAPRIOR:  
      log_prior_ratiorate = (MYREAL (*) (MYREAL,  MYREAL,bayes_fmt *, long)) log_prior_ratio_gamma;
      log_prior_rate = (MYREAL (*) (world_fmt * , long)) log_prior_gamma;
      log_prior_rate1 = (MYREAL (*) (world_fmt * , long, MYREAL)) log_prior_gamma1;
      propose_newrate = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, MYREAL, MYREAL * )) propose_gamma_newparam;
      hastings_ratiorate = (MYREAL (*) (MYREAL, MYREAL, MYREAL, MYREAL, bayes_fmt *, long )) hastings_ratio_gamma;
      break;
    default:
      error("Prior distribution hookup for mutation rate parameter failed");
      break;
    }
}



/// \brief fixes start parameter that are outside of the bayesian mini/maxi values
///
/// Fixes parameter start values that are outside of the Bayesian minimum and maximum 
/// values. Only called in a Bayesian analysis
void bayes_check_and_fix_param(world_fmt *world, option_fmt *options)
{
    long i;
    long frompop, topop;
    MYREAL theta;
    MYREAL *maxparam = world->bayes->maxparam;
    MYREAL *minparam = world->bayes->minparam;
    
    for (i=0; i< world->numpop; i++)
    {
      if(world->param0[i] > maxparam[i])
	world->param0[i] = maxparam[i];
      else
	if(world->param0[i] < minparam[i])
	  world->param0[i] = minparam[i];
    }
    
    for (i=world->numpop; i< world->numpop2; i++)
    {
        if(options->custm2[i]!='0')
        {
	  if(world->options->usem)
	    {
	      if(world->param0[i] > maxparam[i])
                world->param0[i] = maxparam[i];
	      else
                if(world->param0[i] < minparam[i])
		  world->param0[i] = minparam[i];
	    }
	  else
	    {
	      m2mm (i, world->numpop, &frompop, &topop);
	      theta = world->param0[topop];
	      if(world->param0[i]*theta > maxparam[i])
                world->param0[i] = maxparam[i]/theta;
	      else
                if(world->param0[i]*theta < minparam[i])
		  world->param0[i] = minparam[i]/theta;
	    }
	}
    }
}



/// \brief Calculate Prob(g|param)Prob(D|G) from world->treetimes.
///
/// Calculate Prob(g|param)Prob(D|G) from world->treetimes
// uses as input    vtlist *tl = world->treetimes->tl;
// uses as input    world->treetimes->T
// to allow for custom p(g*|theta) calculation, but currently not needed
//MYREAL probg_treetimesX(world_fmt *world, vtlist *tl, long T);
MYREAL probg_treetimes(world_fmt *world)
//{
//  return probg_treetimesX(world,world->treetimes->tl,world->treetimes->T);
//}
//MYREAL probg_treetimesX(world_fmt *world, vtlist *tl, long T)
{
    const long numpop = world->numpop;
    const long numpop2 = world->numpop2;
    const long locus = world->locus;
    const MYREAL like = world->likelihood[world->G];
    long T = world->treetimes->T;
    vtlist *tl = world->treetimes->tl; // comment this out to
    // making probgtreetimes more general
    // so that it can be used in the proposal structure in MCMC1.c
    vtlist *tli = &(tl[0]);
    vtlist *tli1;
    
    long i;
    long r, j;
    long tlif = tli->from;
    long tlit = tli->to;
    int msta, msto;
    
    MYREAL deltatime = tl[0].age;
    MYREAL sumprob = 0.;
    MYREAL eventprob=0.;
   //xcode MYREAL k;
    
    const MYREAL *geo = world->data->geo;
    const MYREAL *lgeo = world->data->lgeo;
    MYREAL *param0 = world->param0;
    MYREAL *locallparam;
    MYREAL *localparam;
    MYREAL *sm;
    const MYREAL mu_rate = world->options->mu_rates[locus];
    const MYREAL lmu_rate = world->options->lmu_rates[locus];

#ifndef LONGSUM
    MYREAL pw;	
    locallparam = (MYREAL *) mycalloc ((numpop2 + numpop + numpop), sizeof (MYREAL));
#else /*LONGSUM*/	
    MYREAL *rates;
    MYREAL *lrates;
    MYREAL *rtimes;
    long partsize = numpop2 + world->bayes->mu + numpop * 3;
    locallparam = (MYREAL *) mycalloc ((numpop2 + numpop + numpop + 9 * numpop), sizeof (MYREAL));
    rates = locallparam + numpop2 + numpop + numpop;
    lrates = rates + 3 * numpop;
    rtimes = rates +  6 * numpop;
    // rates and lrates=log(rates) hold the multipliers for Theta in the the timeslices
    memcpy (rates, localparam+partsize - 3 * numpop, sizeof(MYREAL) * 3 * numpop);
    memcpy (lrates, locallparam+partsize - 3 * numpop, sizeof(MYREAL) * 3 * numpop);
    // rtimes hold the ages at the bottom of a rate timeslice 
    memcpy (rtimes, world->flucrates + 3 * numpop, sizeof(MYREAL) * 3 * numpop);
#endif /*LONGSUM*/
    
    // pointers into the localparam vector
    localparam = locallparam + numpop2;
    sm = localparam + numpop;
    
    // fill localparam vector with data
    memcpy (localparam, param0, sizeof (MYREAL) * numpop);	
    for (r = 0; r < numpop; r++)
      {
	locallparam[r] = LOG2 - (log(param0[r]) + lmu_rate);
        localparam[r] = -1. / (localparam[r] * mu_rate); // minus, so that we can loop over all in probG4
        msta = world->mstart[r];
        msto = world->mend[r];
        for (j = msta; j < msto; j++)
	  {
            if (param0[j] > 0.0)
	      {
                sm[r] -= geo[j] * param0[j] / mu_rate; //minus, so that we can loop over all in probG4
                locallparam[j] = LOG(param0[j]) + lgeo[j] - lmu_rate;
	      }
	  }
      }

    if(tl[0].eventnode->type == 't')
      {
	eventprob = 0.0;
	deltatime = 0.0; //the first tip might be in the past,
	// this makes sure that if all sample have dual date > 0 that
	// we do not calculate strange probabilities, this is a problem only for the
	// first entry. 
      }
    else
      {
	// k = tli->lineages[tlit];
	eventprob = ((tlif==tlit) ? (locallparam[tlif]) : locallparam[mm2m(tlif,tlit, numpop)]);
      }
    sumprob = deltatime * probWait(tli->lineages, locallparam, numpop, numpop2) + eventprob;
	
    for(i=1; i<T-1;i++)
	{
	  tli = &tl[i];
	  tli1 = &tl[i-1];
	  tlif = tli->from;
	  tlit = tli->to;
	  
	  if(tli->eventnode->type == 't')
	    {
	      eventprob = 0.0;
	    }
	  else
	    {
	     // k = tli->lineages[tlit];
	      eventprob = ((tlif==tlit) ? (locallparam[tlif]) :  locallparam[mm2m(tlif,tlit, numpop)]);
	    }
	  deltatime = (tli->age - tli1->age);
#ifndef LONGSUM        
	  pw = probWait(tli->lineages, locallparam, numpop, numpop2);
	  sumprob += deltatime * pw + eventprob;
#else
	  error("Bayesian method of multiple timeslices is not yet implemented");
#endif
	}
    myfree(locallparam);
    //printf("ln p(d|m)=%f+%f=%f\n",sumprob, like ,sumprob + like);
    // each point probability adds a log(1/mu) to the sumprob
    // assuming mu is constant this will be log(treetimes) + log(1/mu) equivalents
    // assuming mu is 1 we simply divide the result by treetimes
    // as a result the values will be still too high posterior estimates
    // will be off by the same quantity log(mu)
    return /*-LOG(world->treetimes->T-1) +*/ (sumprob + like) * world->heat;
}


///
/// calculates the waiting time part for prob(g|param)
MYINLINE
MYREAL probWait(long *lines, MYREAL *locallparam, long numpop, long numpop2)
{
  MYREAL *invtheta = locallparam + numpop2;
  MYREAL *msum = invtheta + numpop;
  long j;
  const long line0 = lines[0];
  const long line01 = (line0 * line0) - line0;
  long line1;
  long line11;
  long line2;
  long line21;
  long line3;
  long line31;
  MYREAL probm = 0.;
  MYREAL probth = 0.;
  switch(numpop)
  {
  case 1:
    probth = line01 * invtheta[0];
    //should be zero! probm = msum[0] * line0;
    return probth; // + probm;   
  case 2:
    line1 = lines[1];
    line11 = (line1 * line1) - line1;
    probm = msum[0] * line0 + msum[1] * line1;
    probth = invtheta[0] * line01 + invtheta[1] * line11;
    return probth + probm;
  case 3:
    line1 = lines[1];
    line11 = (line1 * line1) - line1;
    line2 = lines[2];
    line21 = (line2 * line2) - line2;
    probm = msum[0] * line0 + msum[1] * line1 + msum[2] * line2;
    probth = invtheta[0] * line01 + invtheta[1] * line11 + invtheta[2] * line21;
    return probth + probm;
  case 4:
    line1 = lines[1];
    line11 = (line1 * line1) - line1;
    line2 = lines[2];
    line21 = (line2 * line2) - line2;
    line3 = lines[3];
    line31 = (line3 * line3) - line3;
    probm = msum[0] * line0 + msum[1] * line1 + msum[2] * line2 + msum[3] * line3;
    probth = invtheta[0] * line01 + invtheta[1] * line11 + invtheta[2] * line21 + invtheta[3] * line31;
    return probth + probm;
  default:
    line1 = lines[1];
    line11 = (line1 * line1) - line1;
    line2 = lines[2];
    line21 = (line2 * line2) - line2;
    line3 = lines[3];
    line31 = (line3 * line3) - line3;
    probm = msum[0] * line0 + msum[1] * line1 + msum[2] * line2 + msum[3] * line3;
    probth = invtheta[0] * line01 + invtheta[1] * line11 + invtheta[2] * line21 + invtheta[3] * line31;
    for(j=4; j < numpop; j++)
      {
	const long linex = lines[j];
        probth += ((linex * linex) - linex) * invtheta[j];
	probm += linex * msum[j];
      }
    return probth + probm;
  }
}
///
/// do we accept parameter update
boolean
bayes_accept (MYREAL newval, MYREAL oldval, MYREAL heat, MYREAL hastingsratio)
{
#ifdef SKYTEST
  printf("in skytest: %s %li\n", __file__,__line__);
  return TRUE;
#else
    MYREAL diff = (newval - oldval)  + hastingsratio;// for thermodynamic integration we should not heat this
    // for the standard heating scheme we perhaps should but there is no visible improvement in the runs
    // when we heat this portion for swapping, if we will find some differences I will need to add a flag
    // that differentiates between thermodynamic integration and standard heating.
    if (diff >= 0.0)
        return TRUE;
    if (LOG (RANDUM ()) < diff)
        return TRUE;
    else
        return FALSE;
#endif
}

///
/// log prior ratios for all parameters
MYREAL log_prior_ratio_all(world_fmt *world, MYREAL *newvals)
{
  long np = world->numpop;
  long np2 = world->numpop2;
  long npp = np2 + (world->bayes->mu * world->loci);
  MYREAL logmax = -MYREAL_MAX; 
  long pop;
  MYREAL ratio = 0;
  MYREAL *vals;
  vals = (MYREAL *) mycalloc(npp,sizeof(MYREAL));
  for(pop = 0; pop < np; pop++)
    {
      vals[pop] = log_prior_ratiotheta(newvals[pop], -1., world->bayes,0);
      if(logmax < vals[pop])
	logmax = vals[pop];
    }
  for(pop = np; pop < np2; pop++)
    {
      vals[pop] = log_prior_ratiomig(newvals[pop], -1., world->bayes,0);
      if(logmax < vals[pop])
	logmax = vals[pop];
    }
  if(world->bayes->mu)
    {
      vals[pop] = log_prior_ratiorate(newvals[pop], -1., world->bayes,0);
      if(logmax < vals[pop])
	logmax = vals[pop];
    }
  ratio = 0;
  for(pop = 0; pop < npp; pop++)
    ratio += EXP(vals[pop]-logmax);
  myfree(vals);
  return log(ratio)+logmax;
}

///
/// Log Prior distribution ratios between old and new parameter:
/// UNIFORM
/// $$p[x] = \frac{1}{b-a} $$
/// in reality I use a window to propose new parameters
/// but the distribution should be still the same $p[x]$
/// the prior distribution for the new parameter and the old are the same 
MYREAL log_prior_ratio_uni(MYREAL newparam, 
			   MYREAL oldparam, 
			   bayes_fmt * bayes, 
			   long which)
{
  if((newparam > bayes->maxparam[which]) || (newparam < bayes->minparam[which]))
    return -HUGE;
  else
    return 0.0 ;
}

///
/// Uniform prior distribution for theta or migration rate used in heating acceptance
MYREAL log_prior_uni(world_fmt *world, long numparam)
{
  long numpop = world->numpop;
  long start = ((numparam < numpop || numpop==1) ? 0 : numpop);
  long stop = ((start == numpop) ? world->numpop2 : numpop);
  long i;
  MYREAL * param0 = world->param0;
  bayes_fmt * bayes = world->bayes;
  MYREAL value=0.0;
  MYREAL mu;
  if(numparam>stop) // rate change
    {
      mu = world->options->mu_rates[world->locus];
      i = world->numpop2 + world->locus;
      if((mu > bayes->maxparam[i]) || (mu < bayes->minparam[i]))
	return -HUGE;
      else
	return  -LOG(bayes->maxparam[i] - bayes->minparam[i]);
    }
  // migration or theta parameters
  for(i = start; i < stop; i++)
    {
      if(!strchr("0c", world->options->custm2[i]))
	{
	  if((param0[i] > bayes->maxparam[i]) || (param0[i] < bayes->minparam[i]))
	    return -HUGE;
	  else
	    value += -LOG(bayes->maxparam[i] - bayes->minparam[i]);
	}
    }
  return value ;
}

///
/// uniform prior distribution, returns log(1/(a-b))
MYREAL log_prior_uni1(world_fmt *world, long numparam, MYREAL val)
{
  long numpop = world->numpop;
  long start = ((numparam < numpop || numpop==1) ? 0 : numpop);
  long stop = ((start == numpop) ? world->numpop2 : numpop);
  long i;
  //MYREAL * param0 = world->param0;
  bayes_fmt * bayes = world->bayes;
  MYREAL value=0.0;
  //MYREAL mu;
  if(numparam>stop) // rate change
    {
      //mu = world->options->mu_rates[world->locus];
      i = world->numpop2 + world->locus;
      if((val > bayes->maxparam[i]) || (val < bayes->minparam[i]))
	return -HUGE;
      else
	return  -LOG(bayes->maxparam[i] - bayes->minparam[i]);
    }
  // migration or theta parameters
  for(i = start; i < stop; i++)
    {
      if(!strchr("0c", world->options->custm2[i]))
	{
	  if((val > bayes->maxparam[i]) || (val < bayes->minparam[i]))
	    return -HUGE;
	  else
	    value += -LOG(bayes->maxparam[i] - bayes->minparam[i]);
	}
    }
  return value ;
}

///
/// Log Prior distribution ratios between old and new parameter:
/// EXPONENTIAL
/// $$p[x] = Integrate[Exp[-u/mean]/mean, {u, a, x}]/-Exp[-b/mean] + Exp[-a/mean] $$
/// the ratio between old and new parameter prior will be then
/// $$
///\frac{e^{-\frac{a}{\text{mean}}}-e^{-\frac{x}{\text{mean}}
///}}{e^{-\frac{a}{\text{mean}}}-e^{-\frac{\text{x0}}{\tex
///						     t{mean}}}}
/// $$
///
/// The log(hastingsratio) for this move is not zero but cancels with the log_prior_ratio_exp
/// so instead of calculating stuff that goes away we set in both places the value to zero 
MYREAL log_prior_ratio_exp(MYREAL newparam, MYREAL oldparam, bayes_fmt * bayes, long which)
{
  if((newparam > bayes->maxparam[which]) || (newparam < bayes->minparam[which]))
    return -HUGE;
  else
    return 0.;
  //  MYREAL imean = 1./bayes->priormean[which];
  //MYREAL precalc = exp(-bayes->minparam[which] *imean);
  //MYREAL xo = EXP(-oldparam * imean);
  //MYREAL xn = EXP(-newparam * imean);
  //MYREAL idenom = 1./(precalc - xo );
  //MYREAL nom = (precalc - xn );
  //return log(nom * idenom);
}

///
/// Exponential prior distribution for theta or migration rate used in heating acceptance
/// Out[12]//TextForm=
///   b/mean                  (a - c)/mean
/// E       (-a + c + (-1 + E            ) mean)
/// --------------------------------------------
///                a/mean    b/mean
///              -E       + E
/// a: lower bound
/// b: upper bound
/// c: parameter value
/// (Power(E,b/m)*(-1 + Power(E,(a - c)/m)))/(Power(E,a/m) - Power(E,b/m))
MYREAL log_prior_exp(world_fmt *world, long numparam)
{
  long numpop = world->numpop;
  long start = ((numparam < numpop || numpop==1) ? 0 : numpop);
  long stop = ((start == numpop) ? world->numpop2 : numpop);
  MYREAL * param0 = world->param0;
  MYREAL minp;
  MYREAL maxp;
  MYREAL mean;
  MYREAL imean;
  bayes_fmt *bayes = world->bayes;
  MYREAL value = 0.0;
  MYREAL mu;
  long i;
  if(numparam>stop) // rate change
    {
      i = world->numpop2 + world->locus;
      mu = world->options->mu_rates[world->locus];
      if((mu > bayes->maxparam[i]) || (mu < bayes->minparam[i]))
	    return -HUGE;
	  else
	    {
	      maxp = bayes->maxparam[i];
	      minp = bayes->minparam[i];
	      mean = bayes->meanparam[i];
	      imean = 1.0 / mean;
	      //(Power(E,b/m)*(-1 + Power(E,(a - c)/m)))/(Power(E,a/m) - Power(E,b/m))
	      //return log((EXP(maxp*imean)*(-1.0 + EXP((minp - param0[i])*imean)))/(EXP(minp*imean) - EXP(maxp*imean)));
	      //	      (Power(E,a/m)*(1 - Power(E,(b - c)/m)))/(Power(E,a/m) - Power(E,b/m))
	      	      return log(1.0-(EXP(minp*imean)*(1.0 - EXP((maxp - param0[i])*imean)))/(EXP(minp*imean) - EXP(maxp*imean)));
	      //	      return log(EXP(maxp * imean) * (-minp + param0[i] + EXP((minp-param0[i]) * imean))/
	      //		   (EXP(maxp * imean) - EXP(minp * imean)));
	    }
    }
  // migration and coalescence rate
  for(i = start; i < stop; i++)
    {
      if(!strchr("0c", world->options->custm2[i]))
	{
	  maxp = bayes->maxparam[i];
	  minp = bayes->minparam[i];
	  mean = bayes->meanparam[i];
	  imean = 1.0 / mean;
	  if((param0[i] > bayes->maxparam[i]) || (param0[i] < bayes->minparam[i]))
	    return -HUGE;
	  else
	    {
	      //(Power(E,b/m)*(-1 + Power(E,(a - c)/m)))/(Power(E,a/m) - Power(E,b/m))
	      //value += log((EXP(maxp*imean)*(-1.0 + EXP((minp - param0[i])*imean)))/(EXP(minp*imean) - EXP(maxp*imean)));
	      	      value += log(1.0-(EXP(minp*imean)*(1.0 - EXP((maxp - param0[i])*imean)))/(EXP(minp*imean) - EXP(maxp*imean)));
	      //log(EXP(maxp * imean) * (-minp + param0[i] + EXP((minp-param0[i]) * imean))/
	      //		   (EXP(maxp * imean) - EXP(minp * imean)));
	    }
	}
    }
  return value;
}

///
/// return probability of prior_exp with value
MYREAL log_prior_exp1(world_fmt *world, long numparam, MYREAL value)
{
  long numpop = world->numpop;
  long start = ((numparam < numpop || numpop==1) ? 0 : numpop);
  long stop = ((start == numpop) ? world->numpop2 : numpop);
  //MYREAL * param0 = world->param0;
  MYREAL minp;
  MYREAL maxp;
  MYREAL mean;
  MYREAL imean;
  bayes_fmt *bayes = world->bayes;
  MYREAL localvalue = 0.0;
  //MYREAL mu;
  long i;
  if(numparam>stop) // rate change
    {
      i = world->numpop2 + world->locus;
      //      mu = world->options->mu_rates[world->locus];
      if((value > bayes->maxparam[i]) || (value < bayes->minparam[i]))
	    return -HUGE;
	  else
	    {
	      maxp = bayes->maxparam[i];
	      minp = bayes->minparam[i];
	      mean = bayes->meanparam[i];
	      imean = 1.0 / mean;
	      //(Power(E,b/m)*(-1 + Power(E,(a - c)/m)))/(Power(E,a/m) - Power(E,b/m))
	      //return log((EXP(maxp*imean)*(-1.0 + EXP((minp - value)*imean)))/(EXP(minp*imean) - EXP(maxp*imean)));
	      	      return log(1.0-(EXP(minp*imean)*(1.0 - EXP((maxp - value)*imean)))/(EXP(minp*imean) - EXP(maxp*imean)));
	      //	      return log(EXP(maxp * imean) * (-minp + value + EXP((value) * imean))/
	      //		   (EXP(maxp * imean) - EXP(minp * imean)));
	    }
    }
  // migration and coalescence rate
  for(i = start; i < stop; i++)
    {
      if(!strchr("0c", world->options->custm2[i]))
	{
	  maxp = bayes->maxparam[i];
	  minp = bayes->minparam[i];
	  mean = bayes->meanparam[i];
	  imean = 1.0 / mean;
	  if((value > bayes->maxparam[i]) || (value < bayes->minparam[i]))
	    return -HUGE;
	  else
	    {
	      //(Power(E,b/m)*(-1 + Power(E,(a - c)/m)))/(Power(E,a/m) - Power(E,b/m))
	      //localvalue += log((EXP(maxp*imean)*(-1.0 + EXP((minp - value)*imean)))/(EXP(minp*imean) - EXP(maxp*imean)));
	      localvalue += log(1.-(EXP(minp*imean)*(1.0 - EXP((maxp - value)*imean)))/(EXP(minp*imean) - EXP(maxp*imean)));
	      //	      localvalue += log(EXP(maxp * imean) * (-minp + value + EXP((minp-value) * imean))/
	      //	   (EXP(maxp * imean) - EXP(minp * imean)));
	    }
	}
    }
  return localvalue;
}


///
/// Log Prior distribution ratios between old and new parameter:
/// EXPONENTIAL
/// $$p[x] = Integrate[Exp[-u/mean]/mean, {u, a, x}]/-Exp[-b/mean] + Exp[-a/mean] $$
/// the ratio between old and new parameter prior will be then
/// $$
///\frac{e^{-\frac{a}{\text{mean}}}-e^{-\frac{x}{\text{mean}}
///}}{e^{-\frac{a}{\text{mean}}}-e^{-\frac{\text{x0}}{\tex
///						     t{mean}}}}
/// $$
/// see under log_prior_ratio_exp(), this needs more careful checking as here we
/// depend on current theta and almost certainly the hastings ratio is different from the
/// one with the exp proposal.
MYREAL log_prior_ratio_wexp(MYREAL newparam, MYREAL oldparam, bayes_fmt * bayes, long which)
{
  if((newparam > bayes->maxparam[which]) || (newparam < bayes->minparam[which]))
    return -HUGE;
  else
    return 0.;
}

MYREAL log_prior_wexp(world_fmt *world, long numparam)
{
  return log_prior_exp(world, numparam);
}

MYREAL log_prior_wexp1(world_fmt *world, long numparam, MYREAL value)
{
  return log_prior_exp1(world, numparam, value);
}
///
/// Log Prior distribution ratios between old and new parameter:
/// MULTIPLIER
/// $$p[x] = 1/(lambda m)
/// the ratio between old and new parameter prior will be then
/// $$
///      
/// $$
MYREAL log_prior_ratio_mult(MYREAL newparam, MYREAL oldparam, bayes_fmt * bayes, long which)
{
  if((newparam > bayes->maxparam[which]) || (newparam < bayes->minparam[which]))
    return -HUGE;
  else
    {
      if(newparam > 0.0)
	return log(newparam/oldparam);
      else
	return -MYREAL_MAX;
    }
}


// incorrect
MYREAL log_prior_mult(world_fmt *world, long numparam)
{
  long numpop = world->numpop;
  long start = ((numparam < numpop || numpop==1) ? 0 : numpop);
  long stop = ((start == numpop) ? world->numpop2 : numpop);
  MYREAL * param0 = world->param0;
  MYREAL minp;
  MYREAL maxp;
  MYREAL mean;
  MYREAL imean;
  bayes_fmt *bayes = world->bayes;
  MYREAL value = 0.0;
  MYREAL mu;
  long i;

  if(numparam>stop) // rate change
    {
      i = world->numpop2 + world->locus;
      mu = world->options->mu_rates[world->locus];
      if((mu > bayes->maxparam[i]) || (mu < bayes->minparam[i]))
	    return -HUGE;
	  else
	    {
	      maxp = bayes->maxparam[i];
	      minp = bayes->minparam[i];
	      mean = bayes->meanparam[i];
	      imean = 1.0 / mean;
	      value += log(EXP(maxp * imean) * (-minp + param0[i] + EXP((minp-param0[i]) * imean))/
			   (EXP(maxp * imean) - EXP(minp * imean)));
	    }

    }
  for(i = start; i < stop; i++)
    {
      if(!strchr("0c", world->options->custm2[i]))
	{
	  maxp = bayes->maxparam[i];
	  minp = bayes->minparam[i];
	  mean = bayes->meanparam[i];
	  imean = 1.0 / mean;
	  if((param0[i] > bayes->maxparam[i]) || (param0[i] < bayes->minparam[i]))
	    return -HUGE;
	  else
	    {
	      value += log(EXP(maxp * imean) * (-minp + param0[i] + EXP((minp-param0[i]) * imean))/
			   (EXP(maxp * imean) - EXP(minp * imean)));
	    }
	}
    }
  return value;
}

///
/// 
MYREAL log_prior_mult1(world_fmt *world, long numparam, MYREAL val)
{
  long numpop = world->numpop;
  long start = ((numparam < numpop || numpop==1) ? 0 : numpop);
  long stop = ((start == numpop) ? world->numpop2 : numpop);
  //MYREAL * param0 = world->param0;
  MYREAL minp;
  MYREAL maxp;
  MYREAL mean;
  MYREAL imean;
  bayes_fmt *bayes = world->bayes;
  MYREAL value = 0.0;
  //MYREAL mu;
  long i;

  if(numparam>stop) // rate change
    {
      i = world->numpop2 + world->locus;
      //      mu = world->options->mu_rates[world->locus];
      if((val > bayes->maxparam[i]) || (val < bayes->minparam[i]))
	    return -HUGE;
	  else
	    {
	      maxp = bayes->maxparam[i];
	      minp = bayes->minparam[i];
	      mean = bayes->meanparam[i];
	      imean = 1.0 / mean;
	      value += log(EXP(maxp * imean) * (-minp + val + EXP((minp-val) * imean))/
			   (EXP(maxp * imean) - EXP(minp * imean)));
	    }

    }
  for(i = start; i < stop; i++)
    {
      if(!strchr("0c", world->options->custm2[i]))
	{
	  maxp = bayes->maxparam[i];
	  minp = bayes->minparam[i];
	  mean = bayes->meanparam[i];
	  imean = 1.0 / mean;
	  if((val > bayes->maxparam[i]) || (val < bayes->minparam[i]))
	    return -HUGE;
	  else
	    {
	      value += log(EXP(maxp * imean) * (-minp + val + EXP((minp-val) * imean))/
			   (EXP(maxp * imean) - EXP(minp * imean)));
	    }
	}
    }
  return value;
}

//incorrect
///
/// Log Prior gamma distribution ratios between old and new parameter:
/// GAMMA
/// $$log[gam[x1]/gam[x2]] = -x1/b + x2/b + (-1 + a) Log[x1] + (1 - a) Log[x2]$$
MYREAL log_prior_ratio_gamma(MYREAL newparam, 
			   MYREAL oldparam, 
			   bayes_fmt * bayes, 
			   long which)
{
  MYREAL a = bayes->alphaparam[which];//meanparam holds alpha
  MYREAL ib = bayes->meanparam[which]/a; //assumes that 
  if((newparam > bayes->maxparam[which]) || (newparam < bayes->minparam[which]))
    return -HUGE;
  else
    return ib * (oldparam - newparam) + (a-1) * log(newparam/oldparam);
}

///
/// Gamma prior distribution for theta or migration rate used in heating acceptance
/// log(p[x1]]= -x1/b - a Log[b] + (-1 + a) Log[x1] - Log[Gamma[a]]
MYREAL log_prior_gamma(world_fmt *world, long numparam)
{
  long numpop = world->numpop;
  long start = ((numparam < numpop || numpop==1) ? 0 : numpop);
  long stop = ((start == numpop) ? world->numpop2 : numpop);
  long i;
  MYREAL * param0 = world->param0;
  bayes_fmt * bayes = world->bayes;
  MYREAL a = bayes->alphaparam[numparam];//meanparam holds alpha
  MYREAL ib = bayes->meanparam[numparam]/a; //assumes that 

  for(i = start; i < stop; i++)
    {
      if(!strchr("0c", world->options->custm2[i]))
	{
	  if((param0[i] > bayes->maxparam[i]) || (param0[i] < bayes->minparam[i]))
	    return -HUGE;
	}
    }
  return -param0[numparam] * ib + a * log(ib) + (a - 1.) * log(param0[numparam]) - LGAMMA(a) ;
}

///
/// 
MYREAL log_prior_gamma1(world_fmt *world, long numparam, MYREAL val)
{
  long numpop = world->numpop;
  long start = ((numparam < numpop || numpop==1) ? 0 : numpop);
  long stop = ((start == numpop) ? world->numpop2 : numpop);
  long i;
  //MYREAL * param0 = world->param0;
  bayes_fmt * bayes = world->bayes;
  MYREAL a = bayes->alphaparam[numparam];//meanparam holds alpha
  MYREAL ib = bayes->meanparam[numparam]/a; //assumes that 

  for(i = start; i < stop; i++)
    {
      if(!strchr("0c", world->options->custm2[i]))
	{
	  if((val > bayes->maxparam[i]) || (val < bayes->minparam[i]))
	    return -HUGE;
	}
    }
  return -val * ib + a * log(ib) + (a - 1.) * log(val) - LGAMMA(a) ;
}

///
/// scaling prior used for multiple loci posterior
MYREAL scaling_prior(world_fmt *world, long numparam, MYREAL val)
{ 
    if(numparam < world->numpop2)
      {	
	if(numparam<world->numpop)
	  return (*log_prior_theta1)(world,numparam,val);
	else
	  return (*log_prior_mig1)(world,numparam,val);
      }
    else
      {
	return (*log_prior_rate1)(world, numparam,val);
      }
    return -HUGE;
}
MYREAL calculate_prior(world_fmt *world)
{
  MYREAL nom1 = (*log_prior_theta)(world, world->numpop);
  MYREAL  nom2 = (*log_prior_mig)(world, world->numpop2);
    if(!world->bayes->mu)
      {	
	return nom1 + nom2;
      }
    else
      {
	return nom1 + nom2 +  (*log_prior_rate)(world, world->numpop2);
      }
    return -HUGE;
}


///
/// Uniform flat prior, coding based on Rasmus Nielsen, veryfied with mathematica
/// the correlation among values is dependent on delta
MYREAL
propose_uni_newparam (MYREAL param, MYREAL mean, MYREAL delta, MYREAL minparam, MYREAL maxparam, MYREAL *r)
{
    MYREAL np;
    MYREAL rr = (*r) - 0.5;
    MYREAL thisdelta = rr * delta;
    np = param + thisdelta;
    if (np > maxparam)
      {
	np =  2. * maxparam - np;
      }
 
    else
    {
        if (np < minparam)
	  {
            np = 2. * minparam - np;
	  }
    }   
   return np;
}



///
/// Hastings ratio calculator for uniform distribution
MYREAL hastings_ratio_uni(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
  //  return  (newparam - oldparam)/bayes->priormean[whichparam];
      return 0.;
}



/// \brief Exponential prior distribution

/// Exponential prior using a uniform distribution centered around the old value
/// to allow only smaller changes to the parameters the range for the uniform is
/// arbitrarliy fixed.
///
/// We use the function
///   p[x] = Integrate[Exp[-u/mean]/mean,{u,a,x}]/K
/// K is a constant to integrate to 1 and it is
/// $$
/// -e^{-b/mean} + e^{-a/mean}
/// $$
/// where a is the lower bound and b is an upper bound,
/// The solution to find theta_new is to solve the above for x.
/// and it is
/// $$
/// - m log(-(e^{-\frac{a}{mean}}(r-1))  + 
///                  e^{-\frac{b}{mean}} r)
/// $$
///
MYREAL
propose_exp_newparam (MYREAL param, MYREAL mean,
					  MYREAL delta, MYREAL minparam, MYREAL maxparam, MYREAL * r)
{
    MYREAL np = 0.;
    MYREAL rr = *r;
    np = - mean * log(-(EXP(-minparam/mean)*(rr-1)) + (EXP(-maxparam/mean)*rr));
    return np;
}

///
/// Hastings ratio calculator for exponential distribution
MYREAL hastings_ratio_exp(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
   return 0.; 
}


///
/// uses an exponential prior with boundaries minparm and maxparm and uses also a window around the old parameter
/// the window is of arbitrary size
/// (see propose_exp_newparam() for details)
MYREAL
propose_expb_newparam (MYREAL param, MYREAL mean, MYREAL delta, MYREAL minparam, MYREAL maxparam, MYREAL *r)
{
    MYREAL a = MAX(param - delta, 0.0);
    MYREAL b = param + delta;
    MYREAL rr = (*r);
    MYREAL np = - mean * log(-(EXP(-a/mean)*(rr-1)) + (EXP(-b/mean)* rr));
    while(np < minparam || np > maxparam)
    {
        if(np < minparam)
        {
            np = 2*minparam - np;
        }
        else
        {
            np = 2*maxparam - np;
        }
    }
    return np;
}

///
/// Hastings ratio calculator for exponential distribution
MYREAL hastings_ratio_expb(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
    return 0. ;
}

///
/// Multiplier move, the parameter is multiplied with the rate 1/(lambda*m) 
/// 
/// Requirements: 
///     a < m < b
///     a = 1/b 
///      
/// 
MYREAL
propose_mult_newparam (MYREAL param, MYREAL mean, MYREAL delta, MYREAL minparam, MYREAL maxparam, MYREAL *r)
{
    
    MYREAL multiplier;  // minmult < multiplier < maxmult
    MYREAL maxmult = delta;
    //MYREAL minmult = 1/maxmult;  
    MYREAL lambda = 2. * log(maxmult); // tuning parameter \lambda = 2 ln(b)
    MYREAL np;
    
    multiplier = EXP(lambda * ((*r) - 0.5));
    
    np = multiplier * param;
    (*r) = multiplier;                 
                     
    while(np > maxparam)
    {
        np = maxparam * maxparam / np;
    }
    return np;
}

///
/// Hastings ratio calculator for exponential distribution with multiplicator
/// using the old parameter as mean value
/// simply the jacobian from x->param, matrix if derivatives [should be seocnd right?]
/// [check this with the green paper] 
/// first derivatives 
MYREAL hastings_ratio_mult(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
    // q(x|x')     1/(lambda(1/m))    | dx'/dx       dx'/du  |
    // ------------      -----------------------------    |                      |   =
    // q(x'|x)     1/(lambda m)       | du'/dx       du'/du  |   
    //
    //          |  dmx/dx       dmx/dm     |    
    // = m^2    |                          | = m^2  |(m (-1/m^2) - x zero| = m^2 m/m^2 = m 
    //          |  d(1/m)/dx    d(1/m)/dm  | 
    
    return log(r); // rate multiplier
}

///
/// Gamma prior, NOT WORKING YET
MYREAL
propose_gamma_newparam (MYREAL param, MYREAL mean, MYREAL delta, MYREAL minparam, MYREAL maxparam, MYREAL *r)
{
  //use of incompletegamma(x/alpha,alpha)
  //where x the mean
  // Integrate[(b^-a E^(-(x/b)) x^(-1 + a))/Gamma[a], {x, 0, y}, 
  //  Assumptions -> {a > 0, b > 0}]/(Integrate[(
  //					     b^-a E^(-(x/b)) x^(-1 + a))/Gamma[a], {x, 0, up}, 
  //					    Assumptions -> {a > 0, b > 0}] - 
  //				  Integrate[(b^-a E^(-(x/b)) x^(-1 + a))/Gamma[a], {x, 0, low}, 
  //					    Assumptions -> {a > 0, b > 0}]) // FullSimplify
  // = (Gamma[a] - Gamma[a, y/b])/(Gamma[a, low/b] - Gamma[a, up/b])
  // will give the window between low and up using a(lpha) and b(eta) 
  return -1;
}



///
/// Hastings ratio calculator for gamma distribution
MYREAL hastings_ratio_gamma(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
      return -1.;
}


///
/// verbose log to stdout: report all changes [too much output, this will slow down run a lot]
void print_bayes_verbose(long which, world_fmt *world, MYREAL newparam, boolean success)
{
  long i;
  long frompop;
  long topop;
  MYREAL theta;

  if(world->options->verbose)
    {
      fprintf(stdout,"%i> <%li> ", myID, which);
      for(i=0;i<world->numpop;i++)
	fprintf(stdout,"%f ", which==i ? newparam : world->param0[i]);
      if(world->options->usem)
	{
	  for(i=world->numpop;i<world->numpop2;i++)
	    fprintf(stdout,"%f ", which==i ? newparam : world->param0[i]);
	}
      else
	{
	  for(i=world->numpop;i<world->numpop2;i++)
	    {
	      m2mm(i,world->numpop,&frompop,&topop);
	      theta = world->param0[topop];
	      fprintf(stdout,"%f ", which==i ? newparam * theta : 
		      world->param0[i] * theta);
		  }
	      }
	    if(world->bayes->mu)
	      fprintf(stdout,"%f ", which==world->numpop2 ? newparam : world->options->mu_rates[world->locus]);
	    fprintf(stdout,"%f %c\n",world->param_like,success ? '*': ' ');
	  }
}

///
/// standard metropolis proposal
MYREAL 
uniform_proposal(long which, world_fmt * world, MYREAL *oldparam, boolean *success)
{
#ifdef RANGEDEBUG
  int rangedebug=0;
#endif
  long        topop;
  long        frompop;

  const long  numpop2 = world->numpop2;

  MYREAL      mean;
  MYREAL      newparam;
  MYREAL      oldval;
  MYREAL      r;
  MYREAL      newval;
  MYREAL      hastingsratio;
  MYREAL      theta;
  MYREAL      nm;
  const MYREAL      murate = world->options->mu_rates[world->locus];
  
  bayes_fmt * bayes = world->bayes;

  const boolean verbose = (world->heat > (1.0  - SMALLEPSILON)) && myID==MASTER && world->options->verbose;

  // calculate the probability for the old parameter set
  // this is done from scratch because the tree might have changed in the last step
  mean = bayes->priormean[which];

  oldval = probg_treetimes(world);
  r = UNIF_RANDUM();
  // draw a new parameter from the prior distribution
  // for migration parameters we need to distinguish whether the 
  // prior is in terms of M or xNm
  if(which < world->numpop)
    {
      newparam = (*propose_newtheta) (world->param0[which],
				      mean,
				      bayes->delta[which],
				      bayes->minparam[which],
				      bayes->maxparam[which], &r);
      hastingsratio = (*hastings_ratiotheta)(oldparam[which], newparam, bayes->delta[which], r, bayes, which);
      // add the prior ratio to the log of the hastings ratios
      hastingsratio += (*log_prior_ratiotheta)(oldparam[which], newparam, bayes, which);
      bayes_set_param(world->param0,newparam,which,world->options->custm2, world->numpop);
      reprecalc_world(world, which);
      // calculate prob(G|params)prob(D|G)
      newval = probg_treetimes(world);
    }
  else
    {
      if(which < numpop2)
	{
	  if(world->options->usem)
	    {
	      newparam = (*propose_newmig) (world->param0[which],
					    mean,
					    bayes->delta[which],
					    bayes->minparam[which],
					    bayes->maxparam[which], &r);
	    }
	  else
	    {  
	      m2mm(which,world->numpop,&frompop,&topop);
	      theta = world->param0[topop];
	      nm = world->param0[which] * theta;
	      if(nm > bayes->maxparam[which])
		nm = bayes->maxparam[which];
	      if(nm < bayes->minparam[which])
		nm = bayes->minparam[which];
	      newparam = (*propose_newmig) (nm,
					    mean,
					    bayes->delta[which],
					    bayes->minparam[which],
					    bayes->maxparam[which], &r);
	      newparam /= theta;
#ifdef RANGEDEBUG
	      if(rangedebug==1)
		{
	      printf("%i> min=%f old=%f new=%f (%f) max=%f\n",myID, bayes->minparam[which], nm, newparam*theta, newparam, bayes->maxparam[which]);
	      if(newparam*theta>bayes->maxparam[which])
		warning("shit");
		}
#endif
	    }
	  hastingsratio = (*hastings_ratiomig)(oldparam[which], newparam, bayes->delta[which], r, bayes, which);
	  // add the prior ratio to the log of the hastings ratios
	  hastingsratio += (*log_prior_ratiomig)(oldparam[which], newparam, bayes, which);
	  bayes_set_param(world->param0,newparam,which,world->options->custm2, world->numpop);
	  reprecalc_world(world, which);
	  // calculate prob(G|params)
	  newval = probg_treetimes(world);
	}
      else
	{
	  newparam = (*propose_newrate) (murate,
					mean,
					bayes->delta[which],
					bayes->minparam[which],
					bayes->maxparam[which], &r);
	  hastingsratio = (*hastings_ratiorate)(murate, newparam, bayes->delta[which], r, bayes, which);
	  // add the prior ratio to the log of the hastings ratios
	  hastingsratio += (*log_prior_ratiorate)(murate, newparam, bayes, which);
	  world->options->mu_rates[world->locus] = newparam;
	  world->options->lmu_rates[world->locus] = log(newparam);
	  // change the tree length because of the rate
	  recalc_timelist(world, world->options->mu_rates[world->locus] , oldparam[which]);
	  // calculate prob(G|params)
	  newval = probg_treetimes(world);
	}
    }
  //Acceptance or rejection of the new value
  *success = bayes_accept(newval, oldval,world->heat, hastingsratio);
  if(*success)
    {
      //      if(world->cold && which==numpop2)
      //	printf("*\n");
      /* if(world->in_burnin)
	{
	  bayes->delta[which]*=0.9;
	  } does not work yet*/
      if(verbose)
	print_bayes_verbose(which,world, newparam, *success);
      return newval;
    }
  else
    {
      //      if(world->cold && which==numpop2)
      //	printf("\n");
      if(world->options->prioralone)
	{
	  *success = TRUE;
	  if(verbose)
	    print_bayes_verbose(which,world, newparam, *success);
	  return newval;
	}
      /*      if(world->in_burnin)
	{
	  bayes->delta[which]*=1.1;
	  }does not work yet*/
      if(verbose)
	print_bayes_verbose(which,world, newparam, *success);
      return oldval;
    }
}

void traverse_adjust(node *theNode, MYREAL new_old_ratio)
{
  MYREAL age = 0.;

  if(theNode == NULL)
    return;

  if(theNode->type != 't')
    {
      //      fprintf(stdout,"(%i%c-%f/%f)",myID,theNode->type,theNode->tyme,theNode->tyme*new_old_ratio);
      traverse_adjust(theNode->next->back, new_old_ratio);
      if(theNode->type != 'm')
	traverse_adjust(theNode->next->next->back, new_old_ratio);
    }
  else
    {
      //fprintf(stdout,"(%i%c-%f/%f)\n",myID,theNode->type,theNode->tyme,theNode->tyme*new_old_ratio);
    }
  age = theNode->tyme * new_old_ratio;
  adjust_time_all (theNode, age);
}


///
/// changes branchlength of all branches and the likelihood needs to be calculated completely
void
recalc_timelist (world_fmt * world, MYREAL new_ratio, MYREAL old_ratio)
{
  MYREAL new_old_ratio = new_ratio / old_ratio;
  if (fabs(new_old_ratio - 1.0) < SMALLEPSILON)
    return;  
  //    printf("\n\n+++++++++++++\n%f\n++++++++++++++++++++++++++++++++++++++++++++++++\n",new_old_ratio);
  traverse_adjust(world->root, new_old_ratio);
  set_v (world->root->next->back);
  first_smooth(world,world->locus);
  construct_tymelist (world, &world->treetimes[0]);
#ifdef RATE_DEBUG
  //if(world->cold)
  printf("%i> HEAT: %f, RATE: newrate=%f oldrate=%f oldlike=%f ",myID, world->heat, new_ratio, old_ratio,world->likelihood[world->G]);
#endif
  world->likelihood[world->G] = treelikelihood(world);
#ifdef RATE_DEBUG
  //  if(world->cold)
    printf("newlike=%f\n",world->likelihood[world->G]);
#endif
}

void
recalc_timelist_1 (world_fmt * world, MYREAL new_old_ratio)
{
    long k;
    MYREAL age;

    node *theNode;
    if(new_old_ratio < SMALLEPSILON)
      {
	new_old_ratio = SMALLEPSILON;
      }
    printf("\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    for (k = 0; k < world->treetimes[0].T; k++)
      {
	world->treetimes[0].tl[k].age *= new_old_ratio;
	age = world->treetimes[0].tl[k].age;
	theNode = world->treetimes[0].tl[k].eventnode;
	adjust_time_all (theNode, age);
      }
    for (k = 0; k < world->treetimes[0].T; k++)
      {
	theNode = world->treetimes[0].tl[k].eventnode;
	if (theNode->top == 0 )
	  error("Node is not top\n");
	if(theNode->type=='i')
	  printf("**p(%li):%f l(%li):%f r(%li):%f |||| %li %li\n",theNode->id, theNode->tyme, theNode->next->back->id, theNode->next->back->tyme,  theNode->next->next->back->id, theNode->next->next->back->tyme, world->treetimes[0].tl[k].lineages[0],world->treetimes[0].tl[k].lineages[1]);
	if(theNode->type=='m')
	  printf("**p(%li):%f l(%li):%f r: ---- |||| %li %li\n",theNode->id, theNode->tyme,  theNode->next->back->id, theNode->next->back->tyme, world->treetimes[0].tl[k].lineages[0],world->treetimes[0].tl[k].lineages[1]);
	if(theNode->type=='t')
	  printf("**p(%li):%f l:---- r: ---- |||| %li %li\n",theNode->id, theNode->tyme, world->treetimes[0].tl[k].lineages[0],world->treetimes[0].tl[k].lineages[1]);
	if(theNode->type !='t' && (theNode->next->back->tyme > theNode->tyme))
	  {
	    warning("***********************Time mess: above=%f this=%f\n,",theNode->next->back->tyme , theNode->tyme);
	  }
      }
    treelikelihood(world);
    //    first_smooth(world,world->locus);
    //construct_tymelist (world, &world->treetimes[0]);
}


///
/// Update for parameters using Bayesian inference
long
bayes_update (world_fmt * world)
{
  static unsigned long count=1;
  long counttemp = world->options->lsteps * (world->options->replicatenum > 0 ? world->options->replicatenum : 1);
  
  const unsigned long countmax = (long) (2 * counttemp / (log10((MYREAL)counttemp)));
  
  const boolean     mu = world->bayes->mu;
  boolean           success=FALSE;
  
  long              ba=0;
  long              which;
  long              w; // w is a function of bayes->map[which]
  const long        numpop = world->numpop;
  const long        numpop2 = world->numpop2;
  const long        numpop2rate = world->numpop2 + (mu ? 1 : 0) ;
  // preparation for parameter selection using RANDINT(0..npx)
  const long        npx = numpop2rate - 1;
  long              type;
  //  MYREAL            dummy=0;//for expallslice
  MYREAL            *oldparam=NULL;
  MYREAL            newval = -HUGE;
  const MYREAL      murate = world->options->mu_rates[world->locus];
  bayes_fmt         *bayes = world->bayes;
  
  const boolean     writelog = (myID==MASTER && 
				world->options->progress && 
				(world->heat > 1.0 - SMALLEPSILON));
  const MYREAL      updateratio = world->options->updateratio;
  
  
  
  // progress reporter nested here so that we do need to ask whether we run in Bayes mode or not
  if(!world->in_burnin && writelog)
    count++;
  if(writelog && (count % countmax) == 0)
    {
      bayes_progress(world);
    }    
  // choose genealogy or parameter update
  if(RANDUM() < updateratio) 
    {
      return -1L; // update the tree
    }
  
  // select parameter
  which /*= world->bayes->paramnum*/ = RANDINT(0,npx);
  while(bayes->map[which][1] == INVALID)
    {
      which /*= world->bayes->paramnum*/ = RANDINT(0,npx);
    }
  // savecopy the old parameters
  // we use a memcopy because the custom migration matrix can be filled with m and M and s and S
  // this will change multiple parameters at the same time
  doublevec1d(&oldparam,numpop2rate);
  memcpy(oldparam, world->param0,sizeof(MYREAL) * world->numpop2);
  if(mu)
    oldparam[numpop2] = murate;
  if(which < numpop)
    {
      w = world->bayes->map[which][1]; 
      type = THETAPRIOR;
      if(world->options->slice_sampling[THETAPRIOR])
	{
	  newval = expslice(&world->param0[w], which, world, log_prior_ratiotheta);
	  //expallslice(world->param0, &dummy, which, world);
	  //	  printf("posterior=%f like=%f prior=%f\n",newval,world->likelihood[world->locus],newval - world->likelihood[world->locus]);
	  //newval = world->param0[which];
	  success = TRUE;
	}
      else
	{
	  newval = uniform_proposal(w, world, oldparam, &success);
	}	
    }
  else
    {
      if(which < numpop2)
	{
	  w = world->bayes->map[which][1];
	  type = MIGPRIOR;
	  if(world->options->slice_sampling[MIGPRIOR])
	    {	  
	      newval = expslice(&world->param0[w], which, world, log_prior_ratiomig);
	      //expallslice(world->param0, &dummy, which, world);
	      //newval = world->param0[which];
	      success = TRUE;
	    }
	  else
	    {
	      newval = uniform_proposal(which, world, oldparam, &success);
	    }
	}
      else
	{
	  type = RATEPRIOR;
	  if(world->options->slice_sampling[RATEPRIOR])
	    {
	      newval = slice(&world->options->mu_rates[world->locus], which, world, log_prior_ratiorate);
	      //expallslice(world->param0, &world->options->mu_rates[world->locus], which, world);
	      //newval = world->options->mu_rates[which];
		success = TRUE;
	    }
	  else
	    {
	      newval = uniform_proposal(which, world, oldparam, &success);
	    }
	}
    }
  if(success)
    {
      world->logprior = calculate_prior(world);
      bayes->oldval = newval;
      world->param_like = newval;
      ba = 1;
      bayes->accept[which] += ba;
#ifdef AUTOTUNE
      if(world->options->burn_in)
	{
	  // formula: delta = (pR - 1) * (min-max)
	  // for pR=0.33 and min=0 and max=0.1 --> delta=0.066
	  // pR = delta/(min-max) + 1
	  // Binom(n,k) pR^k (1-pR)^(n-k)
	  // n=1:k=1: delta/(min-max)+1
	  //      delta/(-10*delta)+1 ==> 1-1/10 = 0.9
	  // n=1:k=0: (1-delta/(min-max)-1)) ==> -1/10 = -0.1
	  if(world->options->has_autotune && world->options->autotune < bayes->accept[which]/(1.0+bayes->trials[which]))
	    {
	      if(bayes->delta[which] < bayes->maxparam[which])
		{
		  bayes->delta[which] *= 1.0101; /*0.99^(r/(r-1))*/
		}
	    }
	  //printf("%10li  %f\n",which , bayes->delta[which]);
	}
#endif
    }
  else
    {
#ifdef AUTOTUNE
      if(world->options->burn_in)
	{
	  if(world->options->has_autotune && world->options->autotune > bayes->accept[which]/(1.0+bayes->trials[which]))
	    {
	      if(bayes->delta[which] < bayes->maxparam[which])
		{
		  bayes->delta[which] *= 0.990;
		}
	    }
	}
#endif
      memcpy(world->param0, oldparam, sizeof(MYREAL) * world->numpop2);
      if(type==RATEPRIOR)
	{
	  recalc_timelist(world, oldparam[which], world->options->mu_rates[world->locus]);
	}
      else
	{
	  reprecalc_world(world, which);
	}
      if(mu)
	{
	  world->options->mu_rates[world->locus] = oldparam[numpop2];
	  world->options->lmu_rates[world->locus] = log(oldparam[numpop2]);
	}
      bayes->oldval = newval; // uniform_proposal delivers either old or new value see success
      ba = 0;
    }
  bayes->trials[which] += 1;
  // put in here the prior distribution recorder
  // record_prior(world->bayes->priorhist, world->locus, which, 
  // 
  myfree(oldparam);
  return ba;
}



///
/// fill bayes record in array for histograms and bayesfile
/// set parameter meaning according to option settings
void bayes_save_parameter(world_fmt *world, long pnum, long step)
{
  long i;
  long frompop      = 0;
  long topop        = 0;
  long n            = world->numpop2 + (world->bayes->mu) ;
  long nn           = 2 + n;
  long numpop       = world->numpop;
  long numpop2      = world->numpop2;
  long nnpnum       = nn * pnum;
  boolean mu        = world->bayes->mu;
  //  MYREAL murate     = world->options->meanmu[world->locus] * world->options->mu_rates[world->locus];
  MYREAL murate     = world->options->mu_rates[world->locus];
  MYREAL inheritance_scalar = world->options->inheritance_scalars[world->locus];
  MYREAL * param0   = world->param0;
  worldoption_fmt *wopt = world->options;
  bayes_fmt * bayes = world->bayes;

  (bayes->params+(nnpnum))[0] = bayes->oldval;
    (bayes->params+(nnpnum))[1] = world->likelihood[world->G];
    
    if(wopt->usem)
    {
      memcpy(bayes->params+(nnpnum+2), param0,sizeof(MYREAL)* numpop2);
      if(inheritance_scalar != 1.0)
	{
	  for(i = 0; i < numpop; i++)
	    {
	      (bayes->params+(nnpnum + 2))[i] = param0[i] * inheritance_scalar;
	    }
	}
    }
    else
    {
        memcpy(bayes->params+(nnpnum+2), param0,sizeof(MYREAL)* numpop);
	if(inheritance_scalar != 1.0)
	  {
	    for(i = 0; i < numpop; i++)
	      {
		(bayes->params+(nnpnum + 2))[i] = param0[i] * inheritance_scalar;
	      }
	}
        for(i = numpop; i < numpop2; i++)
        {
            m2mm(i, numpop,&frompop,&topop);
            (bayes->params+(nnpnum + 2))[i] = param0[i] * param0[topop];
        }
    }
    if(mu)
      {
	// we only write one rate per record because the locus is known
	(bayes->params+(nnpnum+2))[numpop2] = murate;
      }
    
    if(wopt->has_bayesmdimfile)
      {
	print_bayes_tofile(world->bayesmdimfile, bayes->params+(nn*pnum), bayes, world, numpop2, bayes->mdiminterval, world->locus, world->replicate, step);
      }
}

///
/// print all parameters continously to a file. This may produce HUGE files
/// and overrun your disk quota 
#ifdef ZNZ
void	  print_bayes_tofile(znzFile mdimfile, MYREAL *params, bayes_fmt *bayes, world_fmt *world, long numpop2, long mdiminterval, long locus, long replicate, long step)
#else
void	  print_bayes_tofile(FILE *mdimfile, MYREAL *params, bayes_fmt *bayes, world_fmt *world, long numpop2, long mdiminterval, long locus, long replicate, long step)
#endif
{
  //  long i;
  long j, j0;
  long interval = mdiminterval > 0 ? mdiminterval : 1;
  boolean *visited;
  long nn = world->numpop2 + (world->bayes->mu);
#ifdef MPI
  float *temp;
  long z=0;
#else
  int fmt = 5;
  long digits;
  long c=0;
  char *temp;
#endif
  bayes->mdimfilecount[locus] += 1;
  if (bayes->mdimfilecount[locus] % interval == 0)
    {
      visited = (boolean*) mycalloc(nn+2,sizeof(boolean));
#ifdef MPI
      temp = (float*) mycalloc(nn+9+world->options->heated_chains+2, sizeof(float));
      temp[0] = (float) step;
      temp[1] = (float) locus;
      temp[2] = (float) replicate;
      temp[3] = (float) params[0];
      temp[4] = (float) params[1];
      temp[5] = (float) params[0] - params[1];
      temp[6] = (float) world->logprior;
      temp[7] = (float) world->treetimes->T-1;
      temp[8] = (float) world->treelen;
      z=9;
#else
      temp = (char*) mycalloc(bayes->linesize, sizeof(char));
      c = sprintf(temp,"%li\t%li\t%li\t%f\t%f\t%f\t%f\t%li\t%f",step, locus+1, replicate+1, params[0], params[1],params[0] - params[1], world->logprior, world->treetimes->T-1, world->treelen);
#endif
      for (j0=0; j0 < nn ; j0++)// the first value is the posterior Log(P(D|G)P(G|p)
	                             // second is log(P(D|G)
	{
	  if(bayes->map[j0][1] == INVALID)
	    continue;
	  else
	    {
	      j = bayes->map[j0][1] + 2;
	      if(visited[j]==TRUE)
		continue;
	      else
		visited[j]=TRUE;
	    }
#ifdef MPI
	  temp[z] = params[j];
	  z++;
	  //  fprintf (stdout, "%f (%li) ",temp[j+6],z);
#else
	  if(params[j] < SMALLEPSILON)
	    {
	      c += sprintf(temp+c,"\t0");
	    }
	  else
	    {
	      digits = (long) log10(params[j]);
	      switch(digits)
		{
		case -8:
		case -6:
		case -5:
		  fmt = 10;
		  break;
		case -4:
		case -3:
		  fmt = 8;
		  break;
		case -2:
		case -1:
		  fmt= 5;
		  break;
		case 0:
		case 1:
		  fmt = 4;
		  break;
		case 2:
		  fmt = 2;
		  break;
		case 3:
		  fmt = 1;
		  break;
		case 4:
		case 5:
		case 6:
		case 7:
		case 8:
		  fmt = 0;
		  break;
		default:
		  if(digits<-8)
		    fmt=20;
		  else
		    fmt = 5;
		  break;
		}
	      c += sprintf(temp+c,"\t%.*f",fmt, params[j]);
	    }
#endif
	}
#ifdef BFDEBUG
#ifdef MPI
      print_marginal_like(temp, &z, world);
#else
      print_marginal_like(temp, &c, world);
#endif
#endif
#ifdef MPI
#ifdef PARALIO
      mpi_mdim_send(&world->mpi_bayesmdimfile,temp,z); 
#else
      mpi_mdim_send(temp,z); 
#endif
#else
      c += sprintf(temp+c,"\n");
      if(!bayes->has_linesize)
	{
	  bayes->linesize = TWO * c;
	  bayes->has_linesize = TRUE;
#ifdef ZNZ
	  unsigned long bytes = bayes->linesize > ONEMEGABYTE ? bayes->linesize : ONEMEGABYTE;
	  znzbuffer(mdimfile,bytes);
#endif
	}
#ifdef ZNZ
      //znzprintf(mdimfile,"%s", temp);
      znzwrite(temp, c, sizeof(char), mdimfile);
#else
      FPRINTF(mdimfile,"%s\n", temp);
#endif 
#endif
      myfree(temp);
      myfree(visited);
    }
}



///
/// Save the Bayesian results for printout into bayesfile
void bayes_save(world_fmt *world, long step)
{
  long np           = 2 + world->numpop2 + (world->bayes->mu * world->loci);
  long pnum         = world->bayes->numparams;
  long allocparams  = world->bayes->allocparams;
  
  if(world->options->has_bayesmdimfile)
    {
      bayes_save_parameter(world, 0, step);
    }
  else
    {
      bayes_save_parameter(world, pnum, step);
      pnum++;
      if(pnum >= allocparams)
	{
	  allocparams += 1000; 
	  world->bayes->params = (MYREAL *) myrealloc(world->bayes->params,sizeof(MYREAL)*allocparams*np);
	}
      world->bayes->numparams = pnum;
      world->bayes->allocparams = allocparams;
    }
}

long setup_bayes_map(longpair *map, char *custm2, long numpop, long numpop2, long size)
{
  long invalid_count = 0 ;
  long j1, j2;
  long n = strlen(custm2);
  long s = MIN(n,numpop2);
  long t = MIN(s,size);
  long i;
  long frompop;
  long topop;
  long first = 0;

  for (i=0; i < t; i++)
    {
      map[i][0] = i;
      switch(custm2[i])
	{
	case '*':
	  map[i][1] = i;
	  break;
	case '0':
	case 'c':
	  map[i][1] = INVALID;
	  invalid_count++;
	  break;
	case 's':
	case 'S':
	  m2mm(i,numpop,&frompop,&topop);
	  if(frompop != topop)
	    {
	      j1 = numpop + topop * (numpop-1) + ((frompop < topop) ? frompop : frompop-1);
	      j2 = numpop + frompop * (numpop-1) + ((topop < frompop) ? topop : topop-1);
	      //printf("j1=%li, j2=%li\n",j1,j2);
	      map[i][1] = j1 < j2 ? j1 : j2;
	    }
	  else
	    map[i][1] = -2; //solves miscoding of sizes with "symmetric" to "mean"
	  break;
	case 'm':
	case 'M':
	  map[i][1] = -2;
	  break;
	default:
	  error("bayes_map needs to know custom migration matrix setting");
	}
    }
  for (i=0; i < numpop; i++)
    {      
      if(map[i][1] == -2)
	{
	  first = map[i][0];
	  break;
	}
    }
  for (i=0; i < numpop; i++)
    {      
      if(map[i][1] == -2)
	{
	  map[i][1] = first;
	}
    }
  for (i=numpop; i < t; i++)
    {      
      if(map[i][1] == -2)
	{
	  first = map[i][0];
	  break;
	}
    }
  for (i=numpop; i < t; i++)
    {      
      if(map[i][1] == -2)
	{
	  map[i][1] = first;
	}
    }

  for (i=t; i < size; i++)
    {
      map[i][0] = i;
      map[i][1] = i;
    }
  return invalid_count;
}

///
/// Initialize the Bayesian framwork
void bayes_init(bayes_fmt *bayes, world_fmt *world, option_fmt *options)
{
  long size = world->numpop2;
  long invalids;
  bayes->numpop2 = size;
  bayes->mu = options->bayesmurates;
  bayes->linesize = MAXBUFSIZE;
  bayes->has_linesize=FALSE;
  bayes->progresslinesize = world->numpop2 * HUNDRED + 5 * HUNDRED;
  if(bayes->mu)
    {
      size += 1; // this is used with bayes->histograms[locus].params so we only need one rate and not loci# rates.
    }
  bayes->map = (longpair *) mycalloc(size, sizeof(longpair));
  invalids = setup_bayes_map(bayes->map, world->options->custm2, world->numpop, world->numpop2, size);
  if(invalids == size)
    {
      options->updateratio = 1.0;
      world->options->updateratio = 1.0;
    }
  bayes->oldval = -MYREAL_MAX;
  bayes->allocparams = 1;
  bayes->numparams = 0;
  bayes->paramnum = 0;
  // datastore for several variables
  bayes->datastore = (MYREAL *) mycalloc(8 * size + 2,sizeof(MYREAL));
  // pointers into datastore
  bayes->priormean = bayes->datastore;
  bayes->delta = bayes->priormean + size;
  bayes->minparam = bayes->delta + size;
  bayes->maxparam = bayes->minparam + size;
  bayes->meanparam = bayes->maxparam + size;
  bayes->deltahist = bayes->meanparam + size;
  //old  bayes->deltahist = (MYREAL *) mycalloc(npp, sizeof(MYREAL));
  bayes->datastore2 = (long *) mycalloc(2 * (size + 2),sizeof(long));
  bayes->accept = bayes->datastore2;
  bayes->trials = bayes->accept + size + 1;

  // records for all bayes derived values
  bayes->params = (MYREAL *) mycalloc(bayes->allocparams * (size+2),sizeof(MYREAL));
  //  bayes->params[0] = (MYREAL *) mycalloc(size+1,sizeof(MYREAL));

  // set counter for mdim file (typically called bayesallfile)
  if(world->cold && options->has_bayesmdimfile)
    {
      bayes->mdimfilecount = (long *) mycalloc(world->loci,sizeof(long));
    }
}


/// initialize the Bayes histogram structure, adds an additional element for the
/// summary over loci.
void bayes_init_histogram(world_fmt * world, option_fmt * options)
{
  long sumloc = world->loci > 1 ? 1 : 0;
    bayes_fmt *bayes = world->bayes;
    long loc;
    long np = world->numpop2;
    long npp = np + (bayes->mu); // the params are like...params...murate_locus 
    long nppt = np + (bayes->mu);// * world->loci); // this includes all different locus-rates
    long pa;

    bayeshistogram_fmt *hist;
    bayes->scaling_factors = (MYREAL *) mycalloc(npp,sizeof(MYREAL));
    bayes->histtotal = (MYREAL *) mycalloc(world->loci * nppt, sizeof(MYREAL));
    bayes->maxmaxvala = -HUGE;
    bayes->prettyhist = options->bayespretty;
    bayes->mdiminterval = options->bayesmdiminterval;
    bayes->histogram = (bayeshistogram_fmt *) mycalloc(world->loci + 1,sizeof(bayeshistogram_fmt));

    for(loc=0; loc < world->loci + sumloc; loc++)
    {
        hist = &(bayes->histogram[loc]);
        hist->bins = (long *) mycalloc(npp, sizeof(long));
        
        hist->results = NULL; //calloc(hist->binsum, sizeof(MYREAL));   // contains histogram, size is bins*numparam
                              // on a per parameter basis
                              // structure has a data storage vectors and the following are all pointers into it
        hist->numparam = npp;    // number of parameters: thetas + migrates + murate
        hist->datastore = (MYREAL *) mycalloc(3 * npp + 8 * nppt, sizeof(MYREAL)); // data storage, size is numparam*11
                                                                        // pointers into data storage
        hist->minima = hist->datastore;    // contains minimal values for each parameter
        hist->maxima = hist->datastore + npp;    // contains maximal values for each parameter
        hist->adjmaxima = hist->datastore + 2*npp;// holds maxima values used in histogram [are smaller than maxima]
	hist->cred50l  = hist->datastore + 3*npp;    // holds 50%-credibility margins (<all lower values>, 
	hist->cred50u = hist->datastore + 3*npp+nppt;   //<all high values>)
	hist->cred95l = hist->datastore + 3*npp+2*nppt;    // holds 95%-credibility margins (<all lower values>)
	hist->cred95u = hist->datastore + 3*npp+3*nppt;   //<all high values>)
	hist->modes = hist->datastore + 3*npp+4*nppt;    // holds 95%-credibility margins (<all lower values>, <all high values>)
	hist->medians = hist->datastore + 3*npp+5*nppt;
	hist->means = hist->datastore + 3*npp+6*nppt;            
	hist->stds = hist->datastore + 3*npp+7*nppt;            
	
	for(pa=0; pa < world->numpop; pa++)
	  {
	    if(!strchr("c", world->options->custm2[pa]))
	      {
		hist->bins[pa] = options->bayespriortheta->bins;        
		hist->binsum += options->bayespriortheta->bins;
		hist->minima[pa] = HUGE;
	      }
	  }
	for(pa=world->numpop; pa < world->numpop2; pa++)
	  {
	    if(!strchr("0c", world->options->custm2[pa]))
	      {
		hist->bins[pa] = options->bayespriorm->bins;        
		hist->binsum += options->bayespriorm->bins;
		hist->minima[pa] = HUGE;
	      }
	    else
	      {
		hist->bins[pa] = 0;        
		hist->minima[pa] = HUGE;
	      }
	  }
	if(world->bayes->mu)
	  {
		hist->bins[pa] = options->bayespriorrate->bins;        
		hist->binsum += options->bayespriorrate->bins;
		hist->minima[pa] = HUGE;
	  }
    }
}

///
/// selects the specific set of prior parameters according to the prior options setting
/// for each parameter with array_count i
MYINLINE  void select_prior_param(int selector, long i, bayes_fmt *bayes, const prior_fmt *prior)       
{
    switch(selector)
    {
    case MULTPRIOR:
    case WEXPPRIOR:
    case EXPPRIOR:
      bayes->priormean[i] = prior->mean; // fill mean for the call to bayes_epxb_newparam
      bayes->delta[i] =  prior->delta ; //(prior->min + prior->max)/(20.); // 1/10 of the max span
      break;
    case UNIFORMPRIOR:
      bayes->priormean[i] = prior->mean; // fill mean for the call to bayes_epxb_newparam
      bayes->delta[i] =  prior->delta; //(prior->max - prior->min)/(10.); // 1/10 of the max span
      break;
    default:
      error("Problems with the specification of the prior distribution");
      break;
    }
    bayes->minparam[i] = prior->min;
    bayes->maxparam[i] = prior->max;
    bayes->meanparam[i] = prior->mean; 
    bayes->deltahist[i] = (prior->max - prior->min)/ prior->bins;
    //printf("%i> deltahist[pa=%li] = %f\n",myID, i, bayes->deltahist[i]); 
}

/// fill the Bayesian framework with values
void bayes_fill(world_fmt *world, option_fmt *options)
{
    long i;
    long locus;

    bayes_fmt * bayes = world->bayes;
    which_theta_prior(options->bayesprior[THETAPRIOR]);
    which_mig_prior(options->bayesprior[MIGPRIOR]);
    which_rate_prior(options->bayesprior[RATEPRIOR]);
    
    for(i=0; i< world->numpop;i++)
    {
        select_prior_param(options->bayesprior[THETAPRIOR], i, bayes, options->bayespriortheta);
    }

    for(i=world->numpop; i< world->numpop2;i++)
    {        
        select_prior_param(options->bayesprior[MIGPRIOR], i, bayes, options->bayespriorm);
    }
    // the rate parameter is set as an additional parameter as the numpop^2 element out of 0..numpop-1 for Theta
    // numpop .. numpop^2-1 for migration.
    if(bayes->mu)
      {
	// set start mu_rate that is compatible with min/max of prior
	select_prior_param(options->bayesprior[RATEPRIOR], world->numpop2, bayes, options->bayespriorrate);
        for(locus=0; locus < options->muloci;locus++)
	  {
	    world->options->mu_rates[locus] = options->bayespriorrate->mean;
	    world->options->lmu_rates[locus] = log(options->bayespriorrate->mean);
	  }
	for(locus=options->muloci;locus < world->loci; locus++)
	  {
	    world->options->mu_rates[locus] = world->options->mu_rates[options->muloci-1];
	    world->options->lmu_rates[locus] = world->options->lmu_rates[options->muloci-1];
	  }
       
      }
    // attach to the custm migration matrix
    bayes->custm2 = world->options->custm2;
}

/// resetting the bayes storage machinery
void bayes_reset(world_fmt * world)
{
  bayes_fmt *bayes = world->bayes;
  long size = world->numpop2 + 2 + (bayes->mu * world->loci);
  if(world->options->bayes_infer)
    {
      bayes->allocparams = 1;
      bayes->params = (MYREAL *) myrealloc(bayes->params, bayes->allocparams * size * sizeof(MYREAL));
      bayes->numparams = 0; // each locus start a new set overwriting the old, allocparam is not reset
    }
}

/// free the Bayesian framework
void bayes_free(world_fmt *world)
{
    long i;
    long sumloc = world->loci > 1 ? 1 : 0;
    if(world->bayes->histogram != NULL)
      {
	for(i=0;i < world->loci + sumloc; i++)
	  {
	    myfree(world->bayes->histogram[i].bins);
	    myfree(world->bayes->histogram[i].datastore);

	    if(myID==MASTER && !world->data->skiploci[i])
	    {
	      myfree(world->bayes->histogram[i].results);
	      //	      printf("bayes results freed for locus %li\n",i);
	      myfree(world->bayes->histogram[i].set95);
		      //printf("bayes set95 freed for locus %li\n",i);

	    }
	  }
	myfree(world->bayes->histogram);
      }
    myfree(world->bayes->datastore);
    myfree(world->bayes->datastore2);
    myfree(world->bayes->params);
    myfree(world->bayes->map);
    myfree(world->bayes);
}


/// calculate the Bayesian statistics and prints the statistics
void bayes_stat(world_fmt *world)
{
  MYREAL threshold;
  MYREAL mutot=0;
  long n=0;
  double mu;
  long lmu;
  double meanmu;
  bayes_fmt * bayes = world->bayes;
  bayeshistogram_fmt *hist;
  long locus;
  long l;
  long frompop=0;
  long topop=0;
  long j=0, j0;
  long size;
  long numpop2 = world->numpop2;
#ifdef DEBUG
  print_bf_values(world);
#endif
  // for single locus data one is not calculating the overall values
  long lozi = world->loci > 1 ? world->loci : 0;
#ifdef LONGSUM
  long addon = (world->options->fluctuate ? (world->numpop * 3) : 0);
#else
 //xcode  long addon = 0;
#endif
  char st[7];  
  char *stemp; 
  
  stemp = (char *) mycalloc(LINESIZE,sizeof(char));
  
  //xcode addon += (bayes->mu * world->loci);
  //CHECK THIS! size = numpop2 + addon;
  size = numpop2 + bayes->mu;
  if(world->loci>1)
    {
      bayes_combine_loci(world);
    }
  // print raw histogram data into the bayesfile
  if(world->options->has_bayesfile)
    {
      print_locus_histogram_header(world->bayesfile, bayes->deltahist, bayes->custm2, world->numpop, size, world->options->usem, bayes->mu, world);
      for(locus=0; locus < world->loci; locus++)
	{
	  if(world->data->skiploci[locus])
	    continue;
	  
	  print_locus_histogram(world->bayesfile, bayes, locus, size);
	}
      if(world->loci>1)
	{
	  print_loci_histogram(world->bayesfile, bayes, world->loci, size);
	}
    }
#ifdef PRETTY
  pdf_print_bayestable(world);
  page_height = pdf_loci_histogram(world);
#endif
  FPRINTF(world->outfile,"\n\n\nBayesian estimates\n");
  FPRINTF(world->outfile,"==================\n\n");
  FPRINTF(world->outfile,"Locus Parameter        2.5%%      25.0%%    mode     75.0%%   97.5%%     median   mean\n");
  FPRINTF(world->outfile,"-----------------------------------------------------------------------------------\n");
  for(locus=0; locus <= lozi; locus++)
    {
      if(locus<world->loci)
	{
	  if(world->data->skiploci[locus])
	    continue;
	}
      hist = &bayes->histogram[locus];
      if(locus == world->loci)
	strcpy(st,"  All ");
      else
	sprintf(st,"%5li ",locus + 1);
      
      for(j0=0; j0< world->numpop2; j0++)
        {
	  
	  //            if(strchr("0c", world->options->custm2[j]))
	  //    continue;
	  if(bayes->map[j0][1] == INVALID)
	    continue;
	  else
	    j = bayes->map[j0][1];
	  
	  if(j < world->numpop)
            {
	      FPRINTF(world->outfile,"%5s ", st);
	      FPRINTF(world->outfile,"Theta_%-3li      ",j0+1);
	      FPRINTF(world->outfile, "%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
		      hist->cred95l[j], hist->cred50l[j], hist->modes[j], 
		      hist->cred50u[j], hist->cred95u[j],hist->medians[j], hist->means[j]);
            }
	  else
            {
	      m2mm(j0,world->numpop,&frompop, &topop);
	      if(world->options->usem)
		sprintf(stemp,"M_%li->%li", frompop+1, topop+1);
	      else
		sprintf(stemp,"Theta_%li*M_%li->%li", topop+1, frompop+1, topop+1);
	      
	      FPRINTF(world->outfile,"%5s ", st);
	      FPRINTF(world->outfile, "%-15.15s",stemp);
	      FPRINTF(world->outfile,"%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
		      hist->cred95l[j], hist->cred50l[j], hist->modes[j],
		      hist->cred50u[j], hist->cred95u[j],hist->medians[j], hist->means[j]);
	    }
	  // warning block: issuing a warning when the 75% and the95% percentiles are within 10% of the 
	  // upper prior boundary
	  threshold =  0.9 * bayes->maxparam[j];
	  if(hist->cred95u[j] > threshold && hist->cred50u[j] > threshold)
	    {
	      if(locus == lozi && locus>0)
		record_warnings(world,"Param %li (all loci): Upper prior boundary seems too low! ",j+1);
	      else
		record_warnings(world,"Param %li (Locus %i): Upper prior boundary seems too low! ", j+1,locus+1);
	    } 
	}
      if(bayes->mu)
	{
	  FPRINTF(world->outfile,"%5.5s ", st);
	  if(locus==lozi && lozi>1)
	    {
	      meanmu = 0.;
	      for(l=0;l<world->loci;l++)
		{
		  meanmu += world->options->meanmu[l];
		}
	      meanmu /= world->loci;
	      lmu = (long) floor( log10(hist->modes[j0] * meanmu));
	      mu = meanmu * pow(10. , -lmu);
	    }
	  else
	    {
	      lmu = (long) floor( log10(hist->modes[j0] * world->options->meanmu[locus]));
	      mu = world->options->meanmu[locus] * pow(10. , -lmu);
	    }
	  FPRINTF(world->outfile," mu [10^%li]     ",lmu);
	  FPRINTF(world->outfile, "%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
		  mu*hist->cred95l[j0], mu*hist->cred50l[j0] , mu*hist->modes[j0], 
		  mu*hist->cred50u[j0], mu*hist->cred95u[j0] , mu*hist->medians[j0], 
		  mu*hist->means[j0]);
	}
    } // over all +1 loci
  FPRINTF(world->outfile,"-----------------------------------------------------------------------------------\n");
    if(bayes->mu)
    {
      mutot = 0;
      n = 1;
      for(locus=0; locus < world->loci; locus++)
  	{
  	  if(world->data->skiploci[locus])
  	    continue;
  	  mu = (double) world->options->meanmu[locus];
  	  mutot += (mu - mutot) / n;
  	  n++;
  	  FPRINTF(world->outfile,"(*) Mutation rate for locus %li was set to %8.5e", locus, mu);
  	}
      if(world->loci>1)
  	{
  	  FPRINTF(world->outfile,"(*) Average mutation rate for all loci is  %8.5e", mutot);
  	}
    }
    // print out the acceptance ratios for every parameter and the tree
    if(world->options->progress)
      {
	// final acceptance ratios for the Bayesian run
	bayes_print_accept(stdout,world);
      }
    myfree(stemp);
}


///
/// print out the acceptance ratios for all the different Bayesian updates
void
bayes_print_accept(FILE * file,  world_fmt *world)
{
    long j;             //used to loop over all parameters
    long topop    =0;   // indicator into the parameter vector, specifying originating population 
    long frompop  =0;   // receiving population
    char *stemp;       // string variable holding print-string 
    long trials   =0;   //
    long tc = world->numpop2 + (world->bayes->mu);// * world->loci); //position of genealogy accept rates
    bayes_fmt *bayes = world->bayes; 

    FPRINTF(file,"\n\n\nAcceptance ratios for all parameters and the genealogies\n");
    FPRINTF(file,"---------------------------------------------------------------------\n\n");
    // This needs more attention but will need more stuff to safe
    if(world->options->datatype == 'g')
      {
	FPRINTF(file,"\n\n\nnot available with datatype=Genealogy\n\n\n");
	return;
      }

    FPRINTF(file,"Parameter           Accepted changes               Ratio\n");

    // population sizes
    for(j=0; j < world->numpop; j++)
    {
        if(!strchr("0c", bayes->custm2[j]))
        {
            if((trials=bayes->trials[j])>0)
            {
                FPRINTF(file,"Theta_%-3li             %8li/%-8li         %8.5f\n", j+1, bayes->accept[j],
                        trials, (MYREAL) bayes->accept[j]/trials);
            }
        }
    }
    // migration rates    
    stemp = (char *) mycalloc(LINESIZE,sizeof(char));
    for(j=world->numpop; j < world->numpop2; j++)
    {
        if(!strchr("0c", bayes->custm2[j]))
        {
            if((trials=bayes->trials[j])>0)
            {
                m2mm (j, world->numpop, &frompop, &topop);
		if(world->options->usem)
		  {
		    sprintf(stemp, "M_%li->%li", frompop+1, topop+1);
		  }
		else
		  {
		    sprintf(stemp, "xN_%lim_%li->%li", topop+1, frompop+1, topop+1);
		  }
                FPRINTF(file, "%-12.12s          %8li/%-8li         %8.5f\n", stemp, bayes->accept[j], 
                        trials, (MYREAL) bayes->accept[j]/trials);
            }
        }
    }
    // accepted rate of mutation rate changes
    if(bayes->mu)
      {
	FPRINTF(file, "Rate of mutation rate %8li/%-8li         %8.5f\n", bayes->accept[j], 
		bayes->trials[j], (MYREAL) bayes->accept[j]/bayes->trials[j]);
      }      
    // accepted trees
    if((trials=bayes->trials[tc])>0)
    {
        FPRINTF(file, "Genealogies           %8li/%-8li         %8.5f\n", bayes->accept[tc], (long)
            trials, (MYREAL) bayes->accept[tc]/trials);
    }
    myfree(stemp);
}


///
/// print out the acceptance ratios for all the different Bayesian updates
void
bayes_progress(world_fmt *world)
{
    // printing machinery
    char *buffer;
    long bufsize=0;
    const boolean writelog = world->options->writelog;
    const boolean progress = world->options->progress;
    char spacer[]="";
    //
    char nowstr[STRSIZE]; // will hold time of day
    char temp[LINESIZE]; // for locus printing with rates
    //long locus;
    long j;         //used to loop over all parameters
    long topop=0;     // indicator into the parameter vector, specifying originating population 
    long frompop=0;   // receiving population
    char *stemp;   // string variable holding print-string 
    long tc = world->numpop2; //this ignores alpha among multiple loci
    bayes_fmt *bayes = world->bayes; 

    buffer = (char *) mycalloc(bayes->progresslinesize,sizeof(char));
    stemp = (char *) mycalloc(LINESIZE,sizeof(char));
        
    get_time (nowstr, "%H:%M:%S");
    if(world->options->heating)
      bufsize += sprintf(buffer+bufsize,
                         "\n%s   [NODE:%i, Locus: %li, Heating:ON (Swaps: %li)] ",
                         nowstr,myID, world->locus + 1, world->swapped);
    else
      bufsize += sprintf(buffer+bufsize,"\n%s   [NODE:%i, Locus: %li] ", nowstr,myID, world->locus + 1);

    prognose_time (nowstr, world, /*world->options->lincr*/ 1, bayes->trials[tc], spacer, TRUE);
    bufsize += sprintf(buffer+bufsize,"\n           (prognosed end of run is %s [%f done])\n", nowstr, 
		       ((MYREAL) world->treesdone / (MYREAL) world->treestotal));

    bufsize += sprintf(buffer+bufsize,"\n           Parameter        Acceptance-ratio         Current value\n           ---------        ----------------         -------------\n"); 

    for(j=0; j < world->numpop; j++)
    {
        if(!strchr("0c", bayes->custm2[j]))
        {
            bufsize += sprintf(buffer+bufsize, "           Theta_%-3li           % 10.5f              %10.5f\n", 
                    j+1, (MYREAL) bayes->accept[j]/bayes->trials[j],world->param0[j]);
        }
    }
    // migration rates
    for(j=world->numpop; j < world->numpop2; j++)
    {
        if(!strchr("0c", bayes->custm2[j]))
        {
            m2mm (j, world->numpop, &frompop, &topop);
	    if(world->options->usem)
	      {
		sprintf(stemp, "M_%li->%li", frompop+1, topop+1);
		bufsize += sprintf(buffer+bufsize, 
				   "           %-12.12s        % 10.5f              %10.5f\n", 
				   stemp, (MYREAL) bayes->accept[j]/bayes->trials[j],
				   world->param0[j]);
	      }        
	    else
	      {
		sprintf(stemp, "xN_%lim_%li->%li", topop+1, frompop+1, topop+1);
		bufsize += sprintf(buffer+bufsize, 
				   "           %-12.12s        % 10.5f              %10.5f\n", 
				   stemp, (MYREAL) bayes->accept[j]/bayes->trials[j],
				   world->param0[j]*world->param0[topop]);
	      }
        }
    }
    if(bayes->mu)
      {
	//	for(locus=0; locus < world->loci; locus++)
	//  {
	//
	//    if(world->data->skiploci[locus])
	//      continue;
	    
	    sprintf(temp,"%li",world->locus);
	    bufsize += sprintf(buffer+bufsize, "           Rate%10.10s      % 10.5f              %10.5f\n",
			       temp, 
			       (MYREAL) bayes->accept[j]/bayes->trials[j],
			       world->options->mu_rates[world->locus]);
	    // }
      }
    // accepted trees
    //bufsize += 
    sprintf(buffer+bufsize, "           Genealogies         % 10.5f              %10.5f\n",  
		       (MYREAL) bayes->accept[tc] /((bayes->trials[tc])),
		       world->likelihood[world->numlike-1]);

    myfree(stemp);
    if(progress)
    {
        FPRINTF(stdout,"%s",buffer);
    }
    if(writelog)
    {
        FPRINTF(world->options->logfile,"%s",buffer);
    }
    myfree(buffer);
    print_bayes_ess(stdout, world, world->numpop2 + world->bayes->mu * world->loci + 1 , 11, world->autocorrelation, world->effective_sample);
}




/// adjusts the allocations for the histogram bins for the current locus 
void adjust_bayes_bins(world_fmt * world, long locus)
{
  //MYREAL dif = 0.;
  //long pa;
  //long np = world->numpop2 + (world->bayes->mu * world->loci);
  bayeshistogram_fmt *hist = &(world->bayes->histogram[locus]);
#ifdef USE_ADJUST_BAYES_BIN
  hist->binsum = 0;
  for(pa=0; pa < np; pa++)
    {
      // if custom migration matrix is set to zero, zero the bins else set it
      if(pa < world->numpop2)
	{
	  if(strchr("0c",world->options->custm2[pa]))
	    {
	      hist->bins[pa] = 0;
	    }
	  else
	    {
	      dif = (hist->maxima[pa] - hist->minima[pa]);
	      if(dif <= 0)
		hist->bins[pa] = 0;
	      else
		hist->bins[pa] = (long) (dif / world->bayes->deltahist[pa]) ;
#ifdef DEBUG_MPI
	      fprintf(stdout,"%i> LOCUS %li: bins[%li]=%li = (%f - %f) /%f \n",myID,locus, pa,hist->bins[pa],hist->maxima[pa], hist->minima[pa], world->bayes->deltahist[pa]);
#endif
	      hist->binsum += hist->bins[pa];
	    }
	}
      else
	{
	  dif = (hist->maxima[pa] - hist->minima[pa]);
	  if(dif <= 0.0)
	    hist->bins[pa] = 0;
	  else
	    hist->bins[pa] = (long) (dif / world->bayes->deltahist[pa]) ;
	  hist->binsum += hist->bins[pa];
	}
    }
  // allocate the found number of bins for the histogram
#endif
  //  if(locus == world->loci && world->options->has_bayesmdimfile)
  //  {
      if(world->bayes->histogram[locus].results == NULL)
	{
	  world->bayes->histogram[locus].results = (MYREAL *) mycalloc(hist->binsum + world->loci * 2, sizeof(MYREAL));//tracing a valgrind issue dec 20 2007
	  world->bayes->histogram[locus].set95 = (char *) mycalloc((hist->binsum+world->loci*2)* 2 + 2, sizeof(char));
	  world->bayes->histogram[locus].set50 = world->bayes->histogram[locus].set95 + hist->binsum + (world->loci * 2) + 1;
	  memset(hist->results, 0 , sizeof(MYREAL) * (hist->binsum));    
	}
      // }
}
 
void construct_param_hist(world_fmt *world, long locus, long numpop2, long pa, long numbin, MYREAL *mini,
			  MYREAL *maxi, MYREAL **results,
			  long *total, MYREAL *themean, MYREAL *thestd)
{
  bayes_fmt *bayes = world->bayes;
  long      j;
  long      i;
  long      np = numpop2 + (bayes->mu);
  long      npx = np + 2;
  long      bin;
  long      nb;

  MYREAL    delta;
  MYREAL    value;
  MYREAL    *data = bayes->params;

  for(i=0;i < bayes->numparams; i++)
    {
      delta = bayes->deltahist[pa];
      value = (data+(npx*i))[pa+2];
      *themean += value;
      *thestd += value * value;
      if(value < mini[pa])
	{
	  warning("data is smaller than minimum in function construct_locus_histogram()\n");
	}
      if(value > maxi[pa])
	{
	  warning("data is larger than maximum in function construct_locus_histogram()\n");
	  bin = bayes->histogram[locus].bins[pa] - 1;
	}
      else
	{
	  bin = (long) ((value - mini[pa])/delta);
	  if(bin<0)
	    bin=0;
	}
      //printf("%i> bin=%li reconstituted=%f value=%f, mini=%f, maxi=%f, delta=%f, locus=%li\n",myID,bin, bin*delta + mini[pa], value,mini[pa],maxi[pa],delta,locus);
      if((bin) > bayes->histogram[locus].bins[pa])
	{
	  warning("%i> value not counted for histogram: bin=%li > histbins=%li\n", myID,bin, bayes->histogram[locus].bins[pa]);
	  *total +=1;
	  return;
	}
      nb = 0;
      for(j=0;j<pa;j++)
	nb += bayes->histogram[locus].bins[j];

      (*results)[nb + bin] += 1.;
      *total += 1;
    }
  //for(j=0;j<maxi[pa]/delta;j++)
  //  printf("%i> %f ",myID,(*results)[j]);
  //printf("\n");
}

/// construct bayes histogram: make a at least BAYESNUMBIN slices through the min-max range for each locus
/// while calculating the histogram calculate also the mean and standard deviation.
/// This can be done in the same loop, but is somehwat opaque.
void construct_locus_histogram(world_fmt *world, long locus, MYREAL *mini, MYREAL *maxi, MYREAL **results)
{
  //    long bin;
    long pa;
    long pa0;
    long rpa;
    //    long i;
    //    MYREAL delta;
    long total=0;
    long numbin=0 ;
    boolean *visited;
    long numpop2 = world->numpop2;
    long np = world->numpop2 + (world->bayes->mu);
    //long npx = np + 2;
    MYREAL themean=0.0;
    MYREAL thestd=0.0;
    visited = (boolean *) mycalloc(np, sizeof(boolean));
    for(pa0=0; pa0 < np ; pa0++)
    {

      if(world->bayes->map[pa0][1] == INVALID)
	continue;
      else
	pa = world->bayes->map[pa0][1];
      if(pa>=numpop2)
	{
	  pa = numpop2;// + locus;
	  rpa = numpop2;
	}
      else
	{
	  rpa = pa;
	}
      if(!visited[rpa])
	{
	  themean = 0.0;
	  thestd = 0.0;
	  total = 0;
	  construct_param_hist(world,locus,world->numpop2, rpa,numbin, mini, maxi, results, &total,&themean,&thestd);
	  world->bayes->histogram[locus].means[pa] = themean/total;
	  world->bayes->histogram[locus].stds[pa]  = thestd / total;
	  world->bayes->histtotal[locus*np+pa] = total;
	  //	  printf("histtotal=%li locus=%li pa=%li\n", total, locus, pa);
	  visited[rpa] = TRUE;
	}
    }
    myfree(visited);
}

///
/// find the largest interval for all loci looking at the stored minima and maxima of the histogram
/// and reports these into adjmini and adjmaxi
/// this allows combination of different histograms without storing all raw data
void adjust_bayes_min_max(world_fmt* world, MYREAL **mini, MYREAL **maxi, MYREAL **adjmini, MYREAL **adjmaxi)
{
  //MYREAL delta;
    long pa0, pa;
    long np = world->numpop2 + (world->bayes->mu);

    for(pa0=0; pa0 < np; pa0++)
    {
        // if custom migration matrix is set to zero
        // continue
	  if(world->bayes->map[pa0][1] == INVALID)
	    continue;
	  else
	    pa = world->bayes->map[pa0][1];

	  //      if(pa < world->numpop2)
	  //{
	  //if(strchr("0c", world->options->custm2[pa]))
          //  continue;
	  //}
      //delta = world->bayes->deltahist[pa];
      if((*mini)[pa] < (*adjmini)[pa0])
        {
	  (*adjmini)[pa0] = (*mini)[pa];
        }
      if((*maxi)[pa] > (*adjmaxi)[pa0])
        {
	  (*adjmaxi)[pa0] = (*maxi)[pa];
        }
    }
}


/// find minimum and maximum values for each parameters 
void find_bayes_min_max(world_fmt* world, MYREAL **mini, MYREAL **maxi, MYREAL **adjmaxi)
{
    long pa0, pa;
    long i;
    //long changed_min;
    //long changed_max;
    long np = world->numpop2 + (world->bayes->mu);
    //xcode long npx = np + 2;
    //xcode MYREAL value;

    //xcode MYREAL delta;
    MYREAL minivalue;
    MYREAL maxivalue;

    for(pa0=0; pa0 < np; pa0++)
      {
        // if custom migration matrix is set to zero
        // continue
	//if(pa < world->numpop2)
	//{
	//  if(strchr("0c", world->bayes->custm2[pa]))
        //    continue;
        //}
      if(world->bayes->map[pa0][1] == INVALID)
	continue;
      else
	pa = world->bayes->map[pa0][1];


      //xcode minivalue = HUGE;
      //xcode maxivalue = 0.;
      //xcode changed_min = 0;
      //xcode changed_max = 0;
      //#ifdef MPI
     // delta = world->bayes->deltahist[pa];
      //#endif
      //	  fprintf(stdout,"%i> npx i pa param value minivalue maxivalue\n",myID);
      for(i=0; i < world->bayes->numparams; i++)
        {
#if 0
          value = (world->bayes->params+(npx*i))[pa+2];
	  if( minivalue > value)
	    {
	      minivalue = value;
	      changed_min++;
	    }
	  if(maxivalue < value)
	    {
	      maxivalue = value;
	      changed_max++;
	    }
	  //	  fprintf(stdout,"%i> %li %li %li %f %f %f %f\n", myID, npx, i, pa, (world->bayes->params+(npx*i))[pa+1], value, minivalue, maxivalue);
        }
      if(maxivalue-minivalue<EPSILON)
        {
#ifdef MPI
	  maxivalue += BAYESNUMBIN * delta;
	  minivalue = 0.0; //only good for parameter 0..inf
#else
	  maxivalue += BAYESNUMBIN * EPSILON;        
#endif
#endif
        }
      //      if(changed_max==0)

	maxivalue = world->bayes->maxparam[pa];
      
	// if(changed_min==0)
	minivalue = world->bayes->minparam[pa];
      
      
      // adjust the mini and maxi so that deltahist will fit in
	//      (*mini)[pa0] =  delta * (long) (minivalue / delta);
	//(*maxi)[pa0] =  delta * (long)(1. + maxivalue / delta);
	(*mini)[pa0] =  minivalue;
	(*maxi)[pa0] =  maxivalue;
      //        fprintf(stdout,"find_bayes_min_max() @@@@@ %f - %f @@@@@ [global: %f - %f]\n",(*mini)[pa],(*maxi)[pa],
      //                world->bayes->histogram[world->loci].minima[pa],world->bayes->histogram[world->loci].maxima[pa]);
    }
}

///
/// prints order of parameters for header files in bayesfile and bayesallfile
void print_param_order(char **buf, long *bufsize, long *allocbufsize, world_fmt *world, long numparam)
{
  long pa;
  long pa0;
  long numpop = world->numpop;
  long numpop2 = world->numpop2;
  long frompop, topop;
  long frompop2, topop2;
  char *custm2 = world->options->custm2;
  boolean usem = world->options->usem;
  bayes_fmt *bayes = world->bayes;
  *bufsize += sprintf(*buf + *bufsize,"# Order of the parameters:\n");
  *bufsize += sprintf(*buf + *bufsize,"# Parameter-number Parameter\n");
  for(pa0=0;pa0<numparam;pa0++)
    {
      if(*bufsize > (*allocbufsize - 150))
	{
	  *allocbufsize += LINESIZE;
	  *buf = (char *) realloc(*buf, sizeof(char) * *allocbufsize);
	}
      if(bayes->map[pa0][1] != INVALID)
	{
	  if(pa0 == bayes->map[pa0][1])
	    pa = pa0;
	  else
	    pa = bayes->map[pa0][1];
	  
	  if(pa0 < numpop)
	    {
	      if((pa0 == pa) && (custm2[pa0] == '*'))
		*bufsize += sprintf(*buf + *bufsize,"#@ %6li    %s_%li\n", pa0+1, "Theta",pa0+1);  
	      else
		*bufsize += sprintf(*buf + *bufsize,"#@ %6li    %s_%li = %s_%li   [%c]\n", 
			pa0+1, "Theta",pa0+1, "Theta",pa+1, custm2[pa0]);  
	    }
	  else
	    {
	      // do we estimate mutation rate changes?
	      if(pa0 == numpop2)
		{
		  *bufsize += sprintf( *buf + *bufsize,"#@ %6li    %s\n", pa0+1, "Rate");
		}
	      else
		{
		  m2mm(pa0,numpop,&frompop,&topop);
		  if((pa0==pa) && (custm2[pa0]=='*'))
		    {
		      if(usem)
			*bufsize += sprintf( *buf + *bufsize,"#@ %6li    %s_(%li,%li)\n", pa+1, "M", frompop+1, topop+1);
		      else  
			*bufsize += sprintf( *buf + *bufsize,"#@ %6li    %s_(%li,%li) = %s_(%li,%li)*%s_%li\n", pa+1, "xNm", frompop+1, topop+1, "M", frompop+1, topop+1, "Theta", topop+1);
		    }
		  else
		    {
		      m2mm(pa,numpop,&frompop2,&topop2);
		      if(usem)
			*bufsize += sprintf(*buf + *bufsize,"#@ %6li    %s_(%li,%li) =  %s_(%li,%li) [%c]\n", pa+1, "M", frompop+1, topop+1, "M", frompop2+1, topop2+1, custm2[pa0]);
		      else  
			*bufsize += sprintf(*buf + *bufsize,"#@ %6li    %s_(%li,%li) = %s_(%li,%li)*%s_%li  = %s_(%li,%li) = %s_(%li,%li)*%s_%li [%c]\n", 
				pa+1, "xNm", frompop+1, topop+1, "M", frompop+1, topop+1, "Theta", topop+1,
				"xNm", frompop2+1, topop2+1, "M", frompop2+1, topop2+1, "Theta", topop2+1,
				custm2[pa0]);
		    }
		}
	    }
	}
    }
}

///
/// prints a comment header, using shell script comments for the output of the raw histogram data for bayesdata
///
/// # Raw data for the histogram of the posterior probabilities for all parameters and loci\n
/// # produced by the program migrate-n VERSIONNUMBER (popgen.csit.fsu.edu/migrate.hml)\n
/// # written by Peter Beerli 2004, Tallahassee, if you have problems email to beerli@fsu.edu\n
/// #
/// # The HPC values are indicators whether the parametervalue is in the highest-posterior credibility set,\n
/// # a 0 means it is outside and a 1 means the value is inside the credibility set.\n
/// #
/// # --------------------------------------------------------------\n
/// # Locus Parameter 50%HPC  95%HPC parameter-value count frequency\n 
/// # --------------------------------------------------------------\n
/// #
void print_locus_histogram_header(FILE *bayesfile, MYREAL *deltas, char *custm2, long numpop, long numparam, boolean usem, boolean mu, world_fmt *world)
{
    long pa;
    long bufsize = 0;
    long allocbufsize = LONGLINESIZE;
    char *buf = (char *) mycalloc(allocbufsize, sizeof(char));
    fprintf(bayesfile, "# Raw data for the histogram of the posterior probabilities for all parameters and loci\n");
    fprintf(bayesfile, "# produced by the program migrate-n %s (http://popgen.scs.fsu.edu/migrate.hml)\n",MIGRATEVERSION);
    fprintf(bayesfile, "# written by Peter Beerli 2004-2007, Tallahassee, if you have problems email to beerli@fsu.edu\n");
    fprintf(bayesfile, "#\n");
    fprintf(bayesfile, "# The HPC values are indicators whether the parametervalue is in the highest-posterior credibility set,\n");
    fprintf(bayesfile, "# a 0 means it is outside and a 1 means the value is inside the credibility set.\n");
    fprintf(bayesfile, "#\n");
    print_param_order(&buf, &bufsize,&allocbufsize, world,numparam);
    fprintf(bayesfile, "#%s\n",buf);
    fprintf(bayesfile, "# Delta for Theta and M ");
    for(pa=0;pa<numparam-1; pa++)
    {
        fprintf(bayesfile,"%f ", custm2[pa]!='0' ? deltas[pa] : -99);
    }
    if(!mu)
      fprintf(bayesfile,"%f ", custm2[pa]!='0' ? deltas[pa] : -99);

    fprintf(bayesfile, "\n# ---------------------------------------------------------------------------------\n");
    fprintf(bayesfile, "# Locus Parameter 50%%HPC  95%%HPC parameter-value count frequency cummulative_freq\n");
    fprintf(bayesfile, "# ---------------------------------------------------------------------------------\n");
    myfree(buf);
}        

///
/// print the locus data for a histogram to the file bayesfile the data has a header starting with # so that
/// other programs easiliy can remove ot for processing. the function calculates the total of all bins and
/// then is printing locusnumber parameternumber 95%HPC 50%HPC bincounts frequency for all bins and parameters
/// the HPC columns are 0 and 1, where 0 means no in the credibiliity set.
/// The header is printed in print_locus_histogram_header(bayesfile)
void print_locus_histogram(FILE *bayesfile, bayes_fmt * bayes, long locus, long numparam)
{
  //  char indicator;
    long bin;
    long pa0, pa;
    long numbins = 0;
    long numbinsall = 0;
    long counterup;
    long counterlow=0;
    //long j;
    MYREAL delta; 
    MYREAL value; 
    MYREAL total=0. ;
    MYREAL freq=0.;
    MYREAL sumfreq=0.;
    
    long *bins = bayes->histogram[locus].bins;
    MYREAL *results = bayes->histogram[locus].results;
    //    MYREAL *modes = bayes->histogram[locus].modes;
    //    MYREAL *medians = bayes->histogram[locus].medians;
    //    MYREAL *means = bayes->histogram[locus].means;
    MYREAL *mini = bayes->histogram[locus].minima;
    MYREAL *maxi = bayes->histogram[locus].maxima;
    char *set50 = bayes->histogram[locus].set50;
    char *set95 = bayes->histogram[locus].set95;
    if(results == NULL)
      return;
    for(pa0=0; pa0 < numparam; pa0++)
    {
      if(bayes->map[pa0][1] == INVALID)
	continue;
      else
	pa = bayes->map[pa0][1];
        //fprintf(stdout,"########### (* locus=%li, parameter %li: *) bins=%li; min=%f; max=%f\n",locus+1, pa+1, bins[pa], mini[pa], maxi[pa]);
        delta = (maxi[pa] - mini[pa])/bins[pa];
        value = mini[pa] + delta/2.;
        total = 0. ;
	counterup=0;
	counterlow = 0;
	//	numbins = 0;
	//for(j=0; j < pa; j++)
	//  numbins += bins[j];
	if(pa<pa0)
	  continue;
	numbinsall += bins[pa];
	numbins = numbinsall - bins[pa];
        for(bin=0;bin<bins[pa];bin++)
        {
	  total += results[numbins + bin];
       
	  if(set95[numbins+bin]=='1')
	    {
	      counterlow++;
	      break;
	    }
	}
        for(bin=counterlow;bin<bins[pa];bin++)
        {
            total += results[numbins + bin];
	  if(set95[numbins+bin]=='1')
	    counterup++ ; //upper 97.5
        }
        sumfreq = 0.;
        for(bin=0;bin<bins[pa];bin++)
	  {
	    freq = results[numbins + bin]/total;
	    sumfreq += freq;
	    //fprintf(stdout,"pa0=%li pa=%li freq=%f sumfreq=%f\n", pa0, pa, freq,sumfreq);
	    /*	  indicator = '-';
		  if(value == modes[pa])
		  indicator = 'T'; //mode
		  if(value == medians[pa])
		  indicator = 'M'; //median
		  if(value == means[pa])
		  indicator = 'A'; //mean
		  if(bin==counterlow)
		  indicator='L';
		  if(bin==counterup)
		  indicator='U';*/
	    //	  fprintf(bayesfile,"%li %li %c %c %c %f %li %f %f\n", locus+1, pa+1, set50[numbins+bin], set95[numbins+bin], indicator, value, 
	    fprintf(bayesfile,"%li %li %c %c %f %li %f %f\n", locus+1, pa+1, set50[numbins+bin], set95[numbins+bin], value, 
		    (long) results[numbins + bin],freq,sumfreq);
	    value += delta;
	  }
        fprintf(bayesfile,"\n");
    }
}


///
/// print the loci data for a histogram to the file bayesfile the data has a header starting with # so that
/// other programs easiliy can remove ot for processing. the function calculates the total of all bins and
/// then is printing locusnumber parameternumber 95%HPC 50%HPC bincounts frequency for all bins and parameters
/// the HPC columns are 0 and 1, where 0 means no in the credibiliity set.
/// The header is printed in print_locus_histogram_header(bayesfile)
void print_loci_histogram(FILE *bayesfile, bayes_fmt * bayes, long locus, long numparam)
{
    long bin;
    //long j;
    long pa0, pa;
    long numbins = 0;
    long numbinsall = 0;
    MYREAL delta; 
    MYREAL value; 
    MYREAL total=0. ;
    MYREAL freq=0. ;
    MYREAL sumfreq=0. ;
    long *bins = bayes->histogram[locus].bins;
    MYREAL *results = bayes->histogram[locus].results;
    MYREAL *mini = bayes->histogram[locus].minima;
    MYREAL *maxi = bayes->histogram[locus].maxima;
    char *set50 = bayes->histogram[locus].set50;
    char *set95 = bayes->histogram[locus].set95;
    
    for(pa0=0; pa0 < numparam; pa0++)
    {
      if(bayes->map[pa0][1] == INVALID)
	continue;
      else
	pa = bayes->map[pa0][1];
      

      delta = (maxi[pa] - mini[pa])/bins[pa];
      value = mini[pa] + delta/2.;
      sumfreq =0.;

      //numbins = 0;
      //for(j=0; j < pa; j++)
      //	numbins += bins[j];    
      if(pa<pa0)
	continue;
      numbinsall += bins[pa];
      numbins = numbinsall - bins[pa];

      for(bin=0;bin<bins[pa];bin++)
        {
	  freq = results[numbins + bin];
	  sumfreq += freq;
	  total = (MYREAL) bayes->trials[pa];
	  fprintf(bayesfile,"%li %li %c %c %f %li %f %f\n", locus+1, pa+1, set50[numbins+bin], set95[numbins+bin], value, 
		  (long) (results[numbins + bin] * total) , freq, sumfreq);
	  value += delta;
        }
        fprintf(bayesfile,"\n");
    }
}


void print_bayes_credibility(FILE *file, MYREAL *cred, MYREAL *results, long numpop)
{
    long pa;
    long numpop2 = numpop * numpop;
    for(pa=0; pa < numpop; pa++)
        fprintf(file,"#Theta lower=%f upper=%f\n",cred[pa], cred[pa+ numpop2]);
    for(pa=numpop; pa < numpop2; pa++)
        fprintf(file,"#Migration lower=%f upper=%f\n",cred[pa], cred[pa+ numpop2]);
}

///
/// print header for complete (raw) parameter outfile (options->bayesmdimfilename)
/// Header uses shell script escapes
/// # Migrate MIGRATEVERSION (Peter Beerli, (c) 2008)
/// # Raw results from Bayesian inference: these values can be used to generate
/// # joint posterior distribution of any parameter combination
/// # and also used for TRACER results [splitting utility may be needed]
#ifdef ZNZ
void print_bayes_mdimfileheader(znzFile file, long interval, world_fmt* world, data_fmt *data)
#else
void print_bayes_mdimfileheader(FILE *file, long interval, world_fmt* world, data_fmt *data)
#endif
{
#ifdef MPI
#ifdef PARALIO
  MPI_Request request;
#endif
#endif
#ifdef BFDEBUG
  long i;
#endif
  long pa;
  long pa0;
  long pop;
  long frompop;
  long topop;
  long repmax;
  boolean *visited;
  bayes_fmt *bayes = world->bayes;
  long bufsize = 0;
  long allocbufsize = LONGLINESIZE;
  char *buf = (char *) mycalloc(allocbufsize,sizeof(char)); //this should be plenty for the header (currently 10^4 characters)
  bufsize = sprintf(buf,"# Migrate %s (Peter Beerli, (c) 2010)\n", MIGRATEVERSION);
  bufsize += sprintf(buf+bufsize,"# Raw results from Bayesian inference: these values can be used to generate\n");
  bufsize += sprintf(buf+bufsize,"# joint posterior distribution of any parameter combination\n");
  bufsize += sprintf(buf+bufsize,"# Writing information on parameters (Thetas, M or xNm)\n");
  bufsize += sprintf(buf+bufsize,"# every %li parameter-steps\n", interval+1);
  bufsize += sprintf(buf+bufsize,"# \n");
  bufsize += sprintf(buf+bufsize,"# --  %s\n", "Steps");
  bufsize += sprintf(buf+bufsize,"# --  %s\n", "Locus");
  bufsize += sprintf(buf+bufsize,"# --  %s\n", "Replicates");
  bufsize += sprintf(buf+bufsize,"# --  %s\n", "log(Posterior)");
  bufsize += sprintf(buf+bufsize,"# --  %s\n", "log(prob(D|G))");
  bufsize += sprintf(buf+bufsize,"# --  %s\n", "log(prob(G|Model))");
  bufsize += sprintf(buf+bufsize,"# --  %s\n", "log(prob(Model))");
  bufsize += sprintf(buf+bufsize,"# --  %s\n", "Sum of time intervals on G");
  bufsize += sprintf(buf+bufsize,"# --  %s\n", "Total tree length of G");
  print_param_order(&buf, &bufsize, &allocbufsize, world, world->numpop2);
  if(world->bayes->mu)
    bufsize += sprintf(buf+bufsize,"# %li  %s\n",  world->numpop2+1, "Rate");
  bufsize += sprintf(buf+bufsize,"# \n");
#ifdef BFDEBUG
  print_marginal_order(buf, &bufsize, world);
#endif
    if (world->options->replicate)
    {
        if (world->options->replicatenum == 0)
            repmax = world->options->lchains;
        else
            repmax = world->options->replicatenum;
    }
    else
        repmax = 1;
  bufsize += sprintf(buf+bufsize,"# \n");
  bufsize += sprintf(buf+bufsize,"#$ ------------------------------------------------------------------ \n");
  bufsize += sprintf(buf+bufsize,"#$ begin  [do not change this content]\n");
  bufsize += sprintf(buf+bufsize,"#$ Model=%s\n",world->options->custm);
  bufsize += sprintf(buf+bufsize,"#$ Mode2=%s\n",world->options->custm2);
  bufsize += sprintf(buf+bufsize,"#$ %li %li %li %li %li %i\n",world->loci, world->numpop,
	  world->numpop2, (long) world->options->replicate, repmax, (int) world->options->usem);
  for (pop = 0; pop < world->numpop; pop++)
    {
      bufsize += sprintf(buf+bufsize,"#$ %s\n",data->popnames[pop]);
    }
  bufsize += sprintf(buf+bufsize,"#$ end\n");
  bufsize += sprintf(buf+bufsize,"#$ ------------------------------------------------------------------ \n");
  bufsize += sprintf(buf+bufsize,"# \n");
  bufsize += sprintf(buf+bufsize,"# remove the lines above and including @@@@@, this allows to use\n");
  bufsize += sprintf(buf+bufsize,"# Tracer (http://tree.bio.ed.ac.uk/software/tracer/) to inspect\n");
  bufsize += sprintf(buf+bufsize,"# this file. But be aware that the current Tracer program (October 2006)\n");
  bufsize += sprintf(buf+bufsize,"# only works with single-locus, single-replicate files\n");
  bufsize += sprintf(buf+bufsize,"# The migrate contribution folder contains a command line utility written\n");
  bufsize += sprintf(buf+bufsize,"# in PERL to split the file for Tracer, it's name is mba\n");
  bufsize += sprintf(buf+bufsize,"# @@@@@@@@\n");
  bufsize += sprintf(buf+bufsize,"Steps\tLocus\tReplicate\tlnPost\tlnDataL\tlnPrbGParam\tlnPrior\ttreeintervals\ttreelength");
  
  visited = (boolean*) mycalloc(world->numpop2,sizeof(boolean));
  for(pa0=0;pa0<world->numpop2;pa0++)
    {
      if(bayes->map[pa0][1] == INVALID)
	continue;
      else
	{
	  pa = bayes->map[pa0][1];
	  if(visited[pa]==TRUE)
	    continue;
	  else
	    visited[pa]=TRUE;
	}
      if(pa < world->numpop)
	{
	  bufsize += sprintf(buf+bufsize,"\t%s_%li", "Theta",pa+1);  
	}
      else
	{
	  m2mm(pa,world->numpop,&frompop,&topop);
	  bufsize += sprintf(buf+bufsize,"\t%s_%li_%li", 
		  (world->options->usem ? "M" : "xNm"), frompop+1, topop+1);
	}
    }
  if(world->bayes->mu)
    {
      bufsize += sprintf(buf+bufsize,"\t%s", "murate");
    }
#ifdef BFDEBUG
  for(i=0;i<world->options->heated_chains;i++)
    bufsize += sprintf(buf+bufsize,"\tmL_%f", world->options->heat[i]);
  if(world->options->heated_chains>1)
    bufsize += sprintf(buf+bufsize,"\tmL_thermo");
  bufsize += sprintf(buf+bufsize,"\tmL_harmonic");
#endif
  bufsize += sprintf(buf+bufsize,"\n");

#ifdef MPI
#ifdef PARALIO
  MPI_File_iwrite_shared(world->mpi_bayesmdimfile,buf,bufsize,MPI_CHAR, &request);
#else
#ifdef ZNZ
  znzprintf(file,"%s",buf);
#else
  fprintf(file,"%s",buf);
#endif
#endif
#else
#ifdef ZNZ
  znzprintf(file,"%s",buf);
#else
  fprintf(file,"%s",buf);
#endif
#endif

  myfree(visited);
  myfree(buf);
}

///
/// finds the mode of the results histogram
/// fills the bayes->histogram[locus].modes
void     find_bayes_mode(bayes_fmt *bayes, long locus, long numparam)
{
    
    long bin;
    long j;
    long numbin=0;
    long pa0, pa;
    long *bins = bayes->histogram[locus].bins;
    
    MYREAL tmp;
    MYREAL *modes = bayes->histogram[locus].modes;
    MYREAL *mini = bayes->histogram[locus].minima;
    MYREAL *results = bayes->histogram[locus].results;
    
    for(pa0=0; pa0 < numparam; pa0++)
    {
        // if custom migration matrix is set to zero
        // continue but ignore rate 
      //if(pa < bayes->numpop2)
      //	{
      //	  if(strchr("0c",bayes->custm2[pa]))
      //     continue;
      // }
      if(bayes->map[pa0][1] == INVALID)
	continue;
      else
	pa = bayes->map[pa0][1];
      
      for(j=0; j < pa; j++)
	numbin += bins[j];    

      tmp = 0. ;
      for(bin=0;bin<bins[pa]; bin++)
	{
	  if (tmp < results[numbin + bin])
	    {
	      tmp = results[numbin+bin];
	      modes[pa] = (MYREAL) bin * bayes->deltahist[pa] + mini[pa];
	    }
	}
      //      numbin += bins[pa];
    }
}

///
/// moving average smoothing
/// a window of size 2*el + 1 is averaged and replaces the central value
/// the array with the values is prepended with el zeros and also appended with el zeros,
/// and the output replaces the input array,
/// this methods uses an Epanechnikov kernel [k(x) = 3/4 (1-x^2) for -1<=x<=1
void bayes_smooth(MYREAL *x, long xelem, long el, boolean lastfirst)
{
  MYREAL d;
  MYREAL dd;
    MYREAL xsum;
    MYREAL *xx;
    MYREAL *weight;
    MYREAL total;
    MYREAL totalx=0.;
    MYREAL totalxsum = 0.;
    long i, j, jj, w;
    long el2 = 2 * el + 1;
    weight = (MYREAL *) mycalloc(el2, sizeof(MYREAL));
    if(xelem==0)
      return;
    xx = (MYREAL *) mycalloc(el2 + xelem,sizeof(MYREAL));
    
    memcpy(xx+el, x, sizeof(MYREAL) * xelem);
    if(lastfirst)
      {
	for(i=0; i < el; i++)
	  {
	    xx[i] = xx[el];
	    xx[i+xelem] = xx[xelem-1];
	  }
      }
    weight[el] = 0.75;
    total = 1.;
    d =  2./(el2-1.);
    for(j = el+1; j < el2; j++)
      {
	dd = (-1. + d*j);
	dd *= dd;
	weight[j] = 0.75 * (1. - dd);
	total += 2 * weight[j];
      }
    //    printf("weight=%f\n",total);
    weight[el] /= total;

    for(j=el+1, jj = el-1; j < el2; j++, jj--)
      {
	weight[j] /= total;
	weight[jj] = weight[j];
      }

    for(i=el; i < xelem + el; i++)
    {
        xsum = 0.;
        for(j=i - el, w=0; j < i + el + 1; j++, w++)
            xsum += xx[j] * weight[w];
	totalxsum += xsum;
	totalx += x[i-el];
	//fprintf(stdout,"%20.20f %20.20f\n",x[i-el], xsum);
        x[i-el] = xsum;
	//x[i-el] = xsum/el2;
    }
    //    fprintf(stdout,"%20.20f %20.20f\n",totalx, totalxsum);
    if(totalxsum>0.0)
      totalx /= totalxsum;
    else
      {
	myfree(xx);
	myfree(weight);
	return;
      }
    for(i=0; i < xelem; i++)
    {
      x[i] *=  totalx;
      //fprintf(stdout,"%20.20f\n",x[i]);
    }
    myfree(xx);
    myfree(weight);
}

///
/// find the credibility set by using the Highest Posterior Probability Density(HPD) and the standard point 
/// descriptors. the statistics are filled into the statistics part of the bayes structure (datastore)
void calc_hpd_credibility(bayes_fmt *bayes,long locus, long numpop2, long numparam)
{
    const MYREAL alpha95 = 0.95;
    const MYREAL alpha50 = 0.50;
    
    //long j;
    long li;
    long pa0, pa;
    long rpa; // the many locus-rates are using only a single rate prior and recording (per locus)
    // and we need to distinguish here 
    long numbins=0;
    long numbinsall = 0;
    long locmedian;
    
    pair *parts; // pairs of doubles
    long csum = 0;
    boolean *smoothy;

    MYREAL biggestcolumn;
    MYREAL total;
    MYREAL cdf;
    MYREAL cutoff95;
    MYREAL cutoff50;
    MYREAL delta;
    MYREAL tmp;
    
    MYREAL * mini;
    MYREAL * maxi;
    MYREAL *results;
    long *bins;
    MYREAL *modes;
    MYREAL *medians;
    MYREAL *cred50;
    MYREAL *cred95;
    char *set50;
    char *set95;

    bins = bayes->histogram[locus].bins;
    modes = bayes->histogram[locus].modes;
    medians = bayes->histogram[locus].medians;
    mini =  bayes->histogram[locus].minima;
    maxi =  bayes->histogram[locus].maxima;
    cred50 =  bayes->histogram[locus].cred50l;
    cred95 =  bayes->histogram[locus].cred95l;
    set50 =  bayes->histogram[locus].set50;
    set95 =  bayes->histogram[locus].set95;
    results =  bayes->histogram[locus].results;

    smoothy = (boolean *) mycalloc(numparam, sizeof(boolean));

    for(pa0=0; pa0 < numpop2 + ((numparam-numpop2)>=1); pa0++)
    {
      if(bayes->map[pa0][1] == INVALID)
	continue;
      else
	pa = bayes->map[pa0][1];

      // so that we can check whether the bins got smoothed already or not
      smoothy[pa0] = FALSE;

      if(pa >= numpop2)
	rpa = numpop2;
      else
	rpa = pa;

      if(csum<bins[rpa])
	csum = bins[rpa];
    }
    parts = (pair *) mycalloc(csum, sizeof(pair));
    
    for(pa0=0; pa0 < numpop2 + ((numparam-numpop2)>=1); pa0++)
    {
      // if custom migration matrix is set to zero
      // continue
      //      if(pa < bayes->numpop2)
      //	{
      //	  if(strchr("0c",bayes->custm2[pa]))
      //     continue;
      // }
      if(bayes->map[pa0][1] == INVALID)
	continue;
      else
	pa = bayes->map[pa0][1];
      
      //      numbins=0;
      //for(j=0; j < pa; j++)
      //	numbins += bins[j];    
      if(pa<pa0)
	continue;


      if(pa >= numpop2)
	rpa = numpop2;
      else
	rpa = pa;


      numbinsall += bins[rpa];
      numbins = numbinsall - bins[rpa];

      delta = (maxi[rpa] - mini[rpa]) / bins[rpa];
      
      parts = (pair *) memset(parts,0,csum*sizeof(pair));
      
      //
      // calculate the total and then the average
      total = 0.;
      biggestcolumn = 0.;
      locmedian = 0;
      //
      // smooth the results before we calculate means etc
      if(!smoothy[pa])
	{
	  //	  bayes_smooth(results+numbins,bins[pa], bins[pa]/50);
	  bayes_smooth(results+numbins,bins[rpa], MAX(3,bins[rpa]/50),TRUE);
	  smoothy[pa] = TRUE;
	}
      for(li=0;li<bins[rpa]; li++)
        {
	  tmp = results[numbins + li];
	  parts[li][1] = tmp;
	  total += tmp;
	  parts[li][0] = delta/2 + mini[rpa] + li * delta;
	  if(tmp > biggestcolumn)
            {
	      biggestcolumn = tmp;
	      locmedian = li;
            }
        }
      //
      // mode is the value with largest column
      modes[pa] = parts[locmedian][0];
      //
      // median is at 0.5 of the cumulative probability distribution
      li = 0;
      tmp = 0;
      while(tmp < total/2  && li < bins[rpa]-1)
        {
	  tmp += parts[li][1];
	  li++;
        }
      medians[pa] = parts[li][0];
      //
      // sort parts using the histogram y-axis for easy finding of quantiles
      paired_qsort2(parts, bins[rpa]);
      
      // find HPD crediblity intervals for 50% and 95% credibility intervals
      // for 50% intervals, starting from the highest value (mode) and moving down
      // until the cumulative sum / total is >=0.5
      cdf = 0.;        
      li = bins[rpa];
      while(cdf < alpha50 && li>0)
        {
	  li--;
	  cdf += parts[li][1] /total;
        }
      cutoff50 = parts[li][1]; // anything lower than this will be in the 50% credibility set
      // or 95% interval
      while(cdf < alpha95 && li>0)
        {
	  li--;
	  cdf += parts[li][1] /total;
        }
      cutoff95 = parts[li][1]; // anything lower than this will be in the 95% credibility set
      
      // fill the credibility sets (are printed into Bayesfile
      for(li=numbins;li<numbins + bins[rpa]; li++)
        {
	  if(results[li] < cutoff95)
            {
	      set50[li] = '0';
	      set95[li] = '0';
            }
	  else
            {
	      set95[li] = '1';
	      if(results[li] < cutoff50)
		set50[li] = '0';
	      else
		set50[li] = '1';
            }
        }
      // fill the innermost 50% levels, smooth over adjacent bins
      // 
      // start at the modus and go left
      li = numbins + locmedian;
      while (set50[li] == '1' && li > numbins)
	--li ;
      cred50[pa] = mini[rpa] + (li-numbins) * delta;// + delta/2.;
      while (set95[li] == '1' && li > numbins)
	--li ;
      cred95[pa] = mini[rpa] + (li-numbins) * delta;// + delta/2.;
      
      // start at the modus and go right
      li = numbins + locmedian;
      while (set50[li] == '1' && li < numbins + bins[rpa])
	++li ;
      cred50[pa + numparam] = maxi[rpa] - (bins[rpa] - li + numbins) * delta;// - delta/2.;
      while (set95[li] == '1' && li < numbins + bins[rpa])
	++li ;        
      cred95[pa + numparam] = maxi[rpa] - (bins[rpa] - li + numbins) * delta;// - delta/2.;
      
      //numbins += bins[rpa];
    }    
    myfree(parts);
    myfree(smoothy);
}

///
/// combines over loci
void bayes_combine_loci(world_fmt * world)
{
    long    locus;
    long    loci = world->loci;
    long    pa0;
    long    pa=0;
    long    rpa;
    long    sourcebin;
    long    targetbin;
    long    sourcesumbin;
    long    targetsumbin;
    long    np = world->numpop2 + world->bayes->mu;// *world->loci
    long    bin;
    long    *bins;

    MYREAL  count;
    MYREAL  sumprob;
    MYREAL  maxval;
    MYREAL  *maxvala;
    MYREAL  total;
    MYREAL  value;
    bayes_fmt * bayes = world->bayes;
    //
    // records whether a cell in the all-loci histogram was filled or not
    char    *touched;
    MYREAL  *priors;
    MYREAL *s_priors_minval;
    MYREAL  *elements;
    //
    // source is the locus histogram, already filled in by the locus workers or then by the locus loop
    bayeshistogram_fmt * source;
    //
    MYREAL s_prior;
    MYREAL inverse_loci_minus_1 = 1;
    long real_loci = world->loci;

    // target is the summary over all loci
    bayeshistogram_fmt *target = &bayes->histogram[world->loci];
    MYREAL   *results;
    MYREAL   *minima;
    MYREAL   *maxima;
    boolean  *visited;
    MYREAL   *mini;
    MYREAL   *maxi;
    MYREAL   *adjmaxi;
    MYREAL   *sumlocus;
    MYREAL sum =0.0;
    //MYREAL w;
    long     tts;
    //
    // allocation of space for maximal values for parameters and set to some large negative number (logs!)   
    doublevec1d(&maxvala, np);
    visited = (boolean *) mycalloc(np,sizeof(boolean));
    sumlocus = (MYREAL *) mycalloc(world->loci,sizeof(MYREAL));
    for(pa0=0;pa0<np; pa0++)
      {
        maxvala[pa0] = -MYREAL_MAX;
	if(bayes->map[pa0][1] == INVALID)
	  continue;
	else
	  {
	    pa = bayes->map[pa0][1];
	    if(pa0 >= world->numpop2)
	      rpa=world->numpop2;
	    else
	      rpa = pa;
	    bayes->histogram[loci].minima[rpa] = MYREAL_MAX;
	    bayes->histogram[loci].maxima[rpa] = -MYREAL_MAX;
	  }
      }
    //
    // set the all-loci minima and maxima
    for(locus=0;locus<world->loci; locus++)
      {
	if(world->data->skiploci[locus])
	  {
	    real_loci -= 1;
	    continue;
	  }
	mini =  world->bayes->histogram[locus].minima;
	maxi =  world->bayes->histogram[locus].maxima;
	adjmaxi =  world->bayes->histogram[locus].adjmaxima;

	find_bayes_min_max(world, &mini, &maxi, &adjmaxi);
	adjust_bayes_min_max(world,&mini,  
			     &maxi, 
			     &bayes->histogram[world->loci].minima,
			     &bayes->histogram[world->loci].maxima);
	//
	//	printf("%i> BAYES A LOCI: locus=%li, binsum=%li\n",myID, locus, bayes->histogram[locus].binsum);
      }
    // allocation of all-loci-histogram table into results, and for the 50% and 95% credibility sets
    adjust_bayes_bins(world, world->loci);
    //printf("%i> BAYES OVER LOCI: locus=%li, binsum=%li\n",myID, world->loci, bayes->histogram[world->loci].binsum);
    bins = target->bins;
    results = target->results;
    minima = target->minima;
    maxima = target->maxima;
    //
    // records whether a cell in the all-loci histogram was filled or not
    touched = (char *) mycalloc(bayes->histogram[loci].binsum + (world->loci*2)+1, sizeof(char));
    // record prior to correct for: P(theta)^(1/n-1)
    priors = (MYREAL *) mycalloc(bayes->histogram[loci].binsum + (world->loci*2)+1, sizeof(MYREAL));
    s_priors_minval = (MYREAL *) mycalloc(bayes->histogram[loci].binsum + (world->loci*2)+1, sizeof(MYREAL));
    elements = (MYREAL *) mycalloc(bayes->histogram[loci].binsum + (world->loci*2)+1, sizeof(MYREAL));
    //
    // set the power for the correct scaling of the posterior
    inverse_loci_minus_1 = 1.0/real_loci - 1.0;
    // calculate the locus weights
    for(locus=0;locus<world->loci; locus++)
      {
      	if(world->data->skiploci[locus])
      	  continue;
	for(pa0=0; pa0 < np; pa0++)
	  {
	    if(bayes->map[pa0][1] == INVALID)
	      continue;
	    else
	      pa = bayes->map[pa0][1];
	    
	    if(pa0 >= world->numpop2)
	      pa=world->numpop2;
	    if(visited[pa]==TRUE)
	      continue;
	    visited[pa] = TRUE;
	    sumlocus[locus] += (MYREAL) world->bayes->histtotal[locus*np+pa]; 
	  }
	sum += sumlocus[locus]; 
	// reset visited
	memset(visited,0,sizeof(boolean)*np);
      }
    // set all bins in target->results to a low log(value)
    for(bin=0;bin<bayes->histogram[world->loci].binsum; bin++)
    {
            results[bin] = -100.;
	    touched[bin] = ' ';
    }
    //  
    // combine previous results into target->results
    // over all loci
    for(locus=0;locus<world->loci; locus++)
    {

      	if(world->data->skiploci[locus])
      	  continue;

        sourcesumbin = 0;
        targetsumbin = 0;

        source = &bayes->histogram[locus];
	memset(visited,0,sizeof(boolean)*np);
        // over all parameters
        for(pa0=0; pa0 < np; pa0++)
	  {
	  if(bayes->map[pa0][1] == INVALID)
	    continue;
	  else
	    pa = bayes->map[pa0][1];

	  if(pa0 >= world->numpop2)
	    rpa=world->numpop2;
	  else
	    rpa = pa;
	  
	  if(visited[rpa]==TRUE)
	    continue;
	  visited[rpa] = TRUE;
	  // placeholder for the maximum value of a parameter
	  maxval = maxvala[rpa] ;
	  //w = sumlocus[locus] / sum;
	  //	  printf("w=%f sumlocus[%li]=%f sum=%f\n",w,locus,sumlocus[locus],sum);
	  //maxval = 0.0;
	  //
	  // number of bins (is this needed and not already adjusted?
	  bins[rpa] = (long) ((maxima[rpa] - minima[rpa])/ bayes->deltahist[rpa]);
	  //fprintf(stdout,"%i> LOCI: bins[%li]=%li = (%f - %f) / %f \n",myID,pa,target->bins[pa],target->maxima[pa], target->minima[pa], world->bayes->deltahist[pa]);
	  // find the entry point for source into the target array
	  // this should always work fine without reminder, but rounding error might produce errors
	  // [this seems ill advised as the 0.5 will increase the bin:] therefore the 0.5 addition and flooring
	  targetbin = (long) ( (source->minima[rpa] - minima[rpa]) / world->bayes->deltahist[rpa]);
	  if(source->results != NULL)
	    {
	      for(sourcebin=0; sourcebin < source->bins[rpa]; sourcebin++)
		{
		  tts = targetsumbin + targetbin + sourcebin;
		  if(source->results[sourcesumbin + sourcebin] > 0.)
		    {
		      // frequencies are normalized per locus and need adjustment for all loci
		      count = source->results[sourcesumbin + sourcebin];
		      results[tts] += log(count); //was w*count
#ifndef OLDBF
		      value = source->minima[rpa] + sourcebin * world->bayes->deltahist[rpa];
		      s_prior = scaling_prior(world,rpa,value);
		      if(s_prior > -HUGE)
			{
			  elements[tts] += 1.0;
			  // calculates mean on the fly
			  if(s_priors_minval[tts] > s_prior)
			    {
			      MYREAL pp = priors[tts];
			      pp = log(pp);
			      pp += s_priors_minval[tts] - s_prior;
			      priors[tts] = exp(pp);
			      s_priors_minval[tts] = s_prior;
			      s_prior = 0.0;
			    }
			  else
			    {
			      s_prior -= s_priors_minval[tts];
			    }
			  // running mean of prior for bin
			  //printf("tts=%li oldmean[%f]=%f, s_prior=%f, elem=%f, ", tts,  elements[tts], priors[tts], s_prior, elements[tts]);
			  priors[tts] += (exp(s_prior) - priors[tts]) / elements[tts];
			  //printf("newmean=%f\n", priors[tts]);
#endif
			  touched[tts] = '1';
			  if(maxval < results[tts])
			    maxval = results[tts];
#ifndef OLDBF
			}
#endif
		    }
		} 
	      sourcesumbin += source->bins[rpa];
	      targetsumbin += target->bins[rpa];
	      maxvala[rpa] = maxval;
	    } 
	  } 
    }
    //
    // adjust with the prior
    //long binsum = bins[0];
    //long z = 0;

    /*    for(bin=0;bin<bayes->histogram[world->loci].binsum; bin++)
//    {
 //     if(bin>=binsum)
//	{
//	  z++;
//	  binsum += bins[z];
//	  maxvala[z] = -HUGE;
//	}
 //     if(touched[bin] == '1')
//	{
//	  results[bin] += (log(priors[bin])+s_priors_minval[bin])*inverse_loci_minus_1;
//	  if(results[bin] > maxvala[z])
//	     maxvala[z] = results[bin];
//	}
 //   }
    */
    //
    // adjust the total and use the overflow safeguard
    targetsumbin = 0;
    memset(visited,0,sizeof(boolean)*np);
    for(pa0 = 0; pa0 < np; pa0++)
      {
	if(world->bayes->map[pa0][1] == INVALID)
	  continue;
	else
	  pa = world->bayes->map[pa0][1];

	if(pa0 >= world->numpop2)
	  rpa=world->numpop2;
	else
	  rpa = pa;

	if(visited[rpa]==TRUE)
	  continue;
	visited[rpa] = TRUE;

	maxvala[rpa] = -HUGE;
	for(bin=targetsumbin; bin < targetsumbin + bins[rpa]; bin++)
	  {
	    if(touched[bin] == '1')
	      {
		//		printf("prior overload! bin=%li results=%f (instead of +%f)\n",bin, results[bin], results[bin]+(log(priors[bin])+ s_priors_minval[bin]) * inverse_loci_minus_1);
		results[bin] += (log(priors[bin])+s_priors_minval[bin])*inverse_loci_minus_1;
		if(results[bin] > maxvala[rpa])
		  maxvala[rpa] = results[bin];
	      }
	  }

	
	target->means[rpa] = 0.;
	sumprob = 0.;
	//	fprintf(stdout,"== start summary \nprobability   running mean of parameter value\n");
	for(bin=targetsumbin; bin < targetsumbin + bins[rpa]; bin++)
	  {//minima[pa] + (bin-targetsumbin)/bins[pa] * (maxima[pa]-minima[pa]), was for reporting in debugmode
	    //  fprintf(stdout,"hist[%6.6li]=%20.20f ", bin, results[bin]);
	    if(touched[bin]=='1')
	      {
		results[bin] = EXP(results[bin] - maxvala[rpa]);
		sumprob += results[bin];
	      }
	    else
	      {
		results[bin] = 0.;
	      }
	    //fprintf(stdout,"%20.20f prior=%f priorcorr=%f %20.20f %20.20f\n", results[bin], priors[bin],(log(priors[bin])+s_priors_minval[bin])*inverse_loci_minus_1, sumprob, maxvala[rpa]);
	  }
	if(bayes->maxmaxvala < maxvala[rpa])
	  bayes->maxmaxvala = maxvala[rpa];
	if(sumprob > 0.0)
	  bayes->scaling_factors[rpa] = log(sumprob) + maxvala[rpa];
	else
	  bayes->scaling_factors[rpa] = -HUGE;
        total = 0;
        for(bin=targetsumbin; bin < targetsumbin + target->bins[rpa]; bin++)
	  {
            if(touched[bin]=='1')
	      {
                results[bin] /= sumprob;
                total +=  results[bin];
                target->means[pa] += results[bin] * (minima[rpa] +  world->bayes->deltahist[rpa]/2. 
						     + world->bayes->deltahist[rpa] * (bin-targetsumbin));
	      }
            //fprintf(stdout,"%5.5li %20.20f %20.20f\n", bin, results[bin], total);
	  }
	//	fprintf(stdout,"== end summary \n");
        targetsumbin += bins[rpa];
      }

    myfree(touched);
    myfree(priors);
    myfree(s_priors_minval);
    myfree(elements);
    myfree(maxvala);
    
    // calculate the credibility intervals
    calc_hpd_credibility(bayes, loci, world->numpop2, world->numpop2+(world->bayes->mu));
    myfree(visited);
    myfree(sumlocus);
}


/// calculates the HPC credibility set and also calculates the posterior singlelocus-histogram
void calculate_credibility_interval(world_fmt * world, long locus)
{
    MYREAL *mini;
    MYREAL *maxi;
    MYREAL *adjmaxi;
    //    long i;
    long np = world->numpop2 + (world->bayes->mu * world->loci);

    mini =  world->bayes->histogram[locus].minima;
    maxi =  world->bayes->histogram[locus].maxima;
    adjmaxi =  world->bayes->histogram[locus].adjmaxima;

    find_bayes_min_max(world, &mini, &maxi, &adjmaxi);
    adjust_bayes_bins(world, locus);
    
     // construct histogram
    construct_locus_histogram(world, locus, mini, maxi, &world->bayes->histogram[locus].results);

    // calc_credibility
    calc_hpd_credibility(world->bayes, locus, world->numpop2, np);
}

///
/// adjust the parameters so that the new set is consistent with the custom migration matrix
void 
bayes_set_param (MYREAL *param, MYREAL newparam, long which, char *custm2, long numpop)
{
    char    migtype = custm2[which];

    long    frompop = 0;
    long    topop   = 0;
    long    i;
    
    MYREAL  nmig;

    //  check custm matrix and then decide
    switch(migtype)
    {
        case 'C':
        case 'c':
            break;
        case 's':
            m2mm (which, numpop, &frompop, &topop);
            param[mm2m(frompop,topop,numpop)] = newparam; // makes the 
            param[mm2m(topop,frompop,numpop)] = newparam; // two parameter equal
            break;
        case 'S':
            m2mm (which, numpop, &frompop, &topop);
            param[mm2m(frompop,topop,numpop)] = newparam; // makes the 
            param[mm2m(topop,frompop,numpop)] = param[mm2m(frompop,topop,numpop)] * param[topop] / param[frompop]; // two parameter equal
            break;
        case 'm':
            if(which<numpop)
                {
		  for(i=0;i<numpop; ++i)
		    {
		      if(custm2[i]=='m')
			param[i] = newparam;
		    }
		}
	    else
	      {
		  for(i=numpop;i<numpop*numpop; ++i)
		    {
		      if(custm2[i]=='m')
			param[i] = newparam;
		    }
		}
            break;
        case 'M':
            if(which<numpop)
            {
                for(i=0;i<numpop; ++i)
                {
                    if(strchr("Mm",custm2[i]))
                        param[i] = newparam;
                }
            }
            else
            {
                m2mm (which, numpop, &frompop, &topop);
                nmig = param[topop] * newparam;
                for(i=numpop;i<numpop*numpop; ++i)
                {
                    if(custm2[i]=='M')
                    {
                        m2mm (i, numpop, &frompop, &topop);
                        param[i] = nmig / param[topop]; // makes the 
                    }
                }
            }
            break;
     case '*':
     default:
         param[which] = newparam;
         break;
    }
}
