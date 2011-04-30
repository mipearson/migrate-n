/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
    variation over time routines   R O U T I N E S

    Peter Beerli 2006, Tallahassee
    beerli@fsu.edu

    Copyright 2006 Peter Beerli, Tallahassee

    This software is distributed free of charge for non-commercial use
    and is copyrighted. Of course, we do not guarantee that the software
    works and are not responsible for any damage you may cause or have.


 $Id$
    */
/*! \file skyline.c 

this file contains functions that calculate the expected parameters based on 
individual time intervals. Results will be printed as a histogram over time (similar to the 
skyline plots of Rambaut, Strimmer etc) and also will print a table/file with values for each parameter

*/
#include <stdlib.h>

#include "skyline.h"
#include "random.h"
#include "tools.h"
#include "sighandler.h"
#include "world.h"
#ifdef PRETTY
#include "pretty.h"
#endif
#ifdef MPI
#include "migrate_mpi.h"
#else
extern int myID;
#endif

///
/// calculate the expected Theta for a coalescence event and returns the expected value
MYREAL waiting_theta(long pop, 
		     MYREAL *param, 
		     MYREAL * mig0list, 
		     long *lineages, 
		     MYREAL interval, 
		     MYREAL inv_mu_rate,
		     long numpop)
{
  long i;
  long kpop;
  MYREAL expected = 0.0;
  
  for(i=0;i<numpop; i++)
    {
      if(lineages[i]>0)
	expected += (lineages[i] * (lineages[i]-1) * (1.0/ param[i]) +  mig0list[i] * lineages[i]) ;
    }
  kpop = lineages[pop] * ( lineages[pop] -1);
  expected = 1./interval - (expected - kpop * 1.0 / param[pop])*inv_mu_rate;
  // 
  if(expected > EPSILON)
    return kpop / expected;
  else
    return -90.;
}

///
/// calculate the expected M for a coalescence event and returns the expected value
MYREAL waiting_M(long pop,long to, 
		 MYREAL *param, 
		 MYREAL * mig0list, 
		 long *lineages, 
		 MYREAL interval, 
		 MYREAL inv_mu_rate,
		 long numpop)
{
  long i;
  long kpop;
  MYREAL expected = 0.0;
  MYREAL line;
  for(i=0;i<numpop; i++)
    {
      line = lineages[i];
      if(line>0)
	expected += (lineages[i] * (lineages[i]-1) * (1.0/ param[i]) +  mig0list[i] * lineages[i]);
    }
  kpop = lineages[to];
  expected = 1./interval - (expected - kpop * param[pop]) * inv_mu_rate ;
  if(kpop > EPSILON)
    return expected / kpop;
  else
    return -9.;
}

// returns the standard deviation and mean of a autocorrelated sequence
// using the the autocorrelation to adjust n for the variance
// where new n = n * (1-r)/(1+r)
// using one-pass algorithms see also in reporter for my unweighted mathematica
// function to calculate r. 
// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
// WEST (1979)
//n = 0
//foreach x in the data:
//if n=0 then 
//    n = 1
//    mean = x
//    S = 0
//    sumweight = weight
//  else
//    n = n + 1
//    temp = weight + sumweight
//	S = S + sumweight*weight*(x-mean)^2 / temp
//	mean = mean + (x-mean)*weight / temp
//      r = r  + xo xn + mo (n mo  - x1 - xo) - mn ((n + 1) mn - x1 - xn); 
//    sumweight = temp
//end if
//end for
//    Variance = S * n / ((n-1) * sumweight)  // if sample is the population, omit n/(n-1)
void weighted_average(float xn, float wn, float *xo, float *wo, 
		      float *sumweight, float *mn, float *mo, 
		      float *var, float *n, float *r)
{
  float temp;
  float invtemp;
  float delta;
  
  if(*n<=0)
    {
      *n = 1;
      *mn = xn;
      *xo = xn;
      *wo = wn;
      *mo = xn;
      *var = 0.;
      *sumweight = wn;
      *r = 0.;
    }
  else
    {
      *n += 1;
      temp = wn + *sumweight;
      invtemp = 1. / temp;
      delta = xn - *mn;
      *var += *sumweight * wn * delta * delta * invtemp;
      *mn += delta*wn * invtemp;
      *r += *sumweight * 0.5 * (wn + *wo) * (xn - *mo) * (*xo - *mo) * invtemp;
      *sumweight = temp;
      *xo = xn;
      *wo = wn;
      *mo = *mn;
    }
}

//n = 0
//mean = 0
//M2 = 0
//foreach x in data:
//  n = n + 1
//  delta = x - mean
//  mean = mean + delta/n
//    M2 = M2 + delta*(x - mean)      // This expression uses the new value of mean
//end for
//variance_n = M2/n
//    variance = M2/(n - 1)
void average(float x, float *mean, float *var, float *n)
{
  //float temp;
  //float invtemp;
  float delta;
  if(*n<=0)
    {
      *n = 1;
      *mean = x;
      *var = 0.;
    }
  else
    {
      *n += 1;
      delta = x - *mean;
      *mean += delta / (*n);
      *var += delta * (x - *mean);
    }
}

///
/// calculates average for value and weight for the skyline bins
void calc_bin_average(tetra *bin, float val, float weight)
{
  float sumweight;
  // calculates the weighted mean and variance component of the value
  // the sumweight is not kept but reconstituted from the weight average.
  sumweight = (*bin)[4] * (*bin)[1]; // average weight * n
  if((*bin)[0]>0.0 && val > 10. * (*bin)[0]) //Saveguard against outliers, a better method than this?
    return;
//syntax:  weighted_average(val, weight, &oldval, &oldweight, &sumweight, &value_mean, &old_mean, &sumweight_mean), &count, &autocorrelation);
 weighted_average(val, weight, &(*bin)[6], &(*bin)[7], &sumweight, &((*bin)[0]), &(*bin)[8], &((*bin)[2]), &((*bin)[4]), &(*bin)[5]);
  // calculates the average weight for this bin
  (*bin)[4] -= 1; // reset for weight average calculation
  average(weight,  &((*bin)[1]), &((*bin)[3]), &((*bin)[4]));
}

///
/// Calculate the expected parameters given the event
void calculate_expected_values(tetra **eventbins, 
			       long *eventbinnum, 
			       MYREAL eventinterval, 
			       MYREAL interval, 
			       MYREAL age, 
			       long from, 
			       long to, 
			       long * lineages, 
			       long numpop, 
			       world_fmt *world)
{
  //  static MYREAL maxlweight = -HUGE;
  //  MYREAL oldlweight = -HUGE;
//xcode  boolean change;
  //  long p;

//xcode  float lweight;
  //  float unweight;
  long i,j;
  long pop;
  long diff;
  long allocsize=1;

  float val;
  // float se1;
  //  float se2;

  float weight1 = 0.0F ;
  float weight2 = 0.0F ;

  MYREAL mu_rate = world->options->mu_rates[world->locus];
  MYREAL inv_mu_rate = 1./ mu_rate;
  MYREAL inv_eventinterval = 1. / eventinterval;
  long firstbin = (long) ((age - interval) * inv_eventinterval);
  //long firstbin = (long) (age * inv_eventinterval);
  long bin     = (long) (age * inv_eventinterval);
  float binage  = eventinterval * bin; //this is on the left side of the bin
  MYREAL firstbinage        = eventinterval * firstbin; //this is the age at the left side of first bin
  tetra *bins=NULL;
  //xcode change =FALSE;
  //xcode lweight = 1.0;
  // do we need to report problems when bin number is above 10000?
  // 
  //  if(age > 1000.0 || bin > 10000)
  if(age > 1000.0)
    return;

  if(firstbin < 0) 
    firstbin=0;


  if(bin < 0) 
    error("bin is smaller than zero");

  if(from == to)
    {
      diff = bin - firstbin;
      weight1 = interval * inv_eventinterval; // or =lweight=1.0
      pop  = to;
      val = waiting_theta(pop,world->param0,world->mig0list, lineages,interval, inv_mu_rate, numpop); 
      if(val <= 0.0)
	return;
      //      printf("%i> k(k-1)=%li val=%f pop=%li\n",myID,k,val, pop);
      //se1  = val * val;
      //se2  = se1;
    } 
  else
    {
      diff = 0;
      weight1 = 1.0;
      pop = mm2m(from, to, numpop);

      val = waiting_M(pop, to, world->param0,world->mig0list, lineages,interval, inv_mu_rate, numpop); 


      if(val < 0.0)
	val = 0.0;
      //DEBUG TEST	return;
      //se1 = val * val;
      //se2 = 1.0;
    }
  if(bin >= eventbinnum[pop])
    {
      allocsize      = bin+10;
      eventbins[pop] = (tetra *) myrealloc(eventbins[pop], allocsize * sizeof(tetra));
      for(i=eventbinnum[pop];i < allocsize; i++)
	{
	  for(j=0;j<9;j++)
	    eventbins[pop][i][j] = 0.;
	}
      eventbinnum[pop] = allocsize;
    }
  // weights are the contribution of the estimate to the bin
  if(diff < 0)
      error("time difference is negative -- not possible");

  bins = eventbins[pop];

  switch(diff)
    {
      //start and stop bin are the same the whole interval does not extend over bin boundaries
    case 0: 
      // weight1 is calculated earlier to make sure that migration weights are 1.0 and the
      // coalescence weights are the standard weight1
      calc_bin_average(&bins[firstbin], val, weight1);
      break;
    case 1: 
      // the interval crosses a single bin boundary
      // older bin
      //      weight1 = (age - binage) * inv_eventinterval; 
      weight1 = (-age + interval + firstbinage + eventinterval) * inv_eventinterval;
      calc_bin_average(&bins[firstbin], val, weight1);
      // newer bin
      weight2 = (age - binage) * inv_eventinterval;
      calc_bin_average(&bins[bin], val, weight2);
      break;
    default:
      // start and stop are in different bins and more than one boundary apart
      weight1 = (-age + interval + firstbinage + eventinterval) * inv_eventinterval;
      calc_bin_average(&bins[firstbin], val, weight1);
      for(i=firstbin+1; i < bin; i++)
	{
	  calc_bin_average(&bins[i], val, 1.0);
	}
      weight2 = (age - binage) * inv_eventinterval;
    calc_bin_average(&bins[bin], val, weight2);
    }
  //  printf("@%c %li %li :%li:  %f %f %f %f %f %f %f\n", from==to ? 'c' : 'm', from, to, bin, val, age, bins[bin][0],bins[bin][1], lweight, weight1, weight2); 
  //  if(pop==0 && bin==15)
  //  printf("@ %li %li | %f | %f %f %f\n",pop, bin, val, bins[bin][0], bins[bin][2] / ((bins[bin][4]*((1.- bins[bin][5])/(1.+ bins[bin][5]))-1.) * bins[bin][1]),bins[bin][5]/bins[bin][2]);
}


///
/// set up the skyline plot histogram containers, only when also the migration histograms are recorded
/// this function needs to be called AFTER the setup_mighist() function
/// the main vehicle is the evenbins that contain 9 basic bins for each time bin, each parameter
/// 0) Time
/// 1) Parameter value (average)
/// 2) standard deviation  (average standard deviation)
/// 3) standard deviation of weights (average weight)
/// 4) Number of values put into bin (average count)
/// 5) autocorrelation (average autocorrelation coeff)
/// 6)
/// 7)
/// 8)
void
setup_expected_events (world_fmt * world, option_fmt * options)
{
  long locus, i, j,z;
    long allocsize = 1;
    MYREAL  binsize = world->options->eventbinsize; // in mutational units (=time scale)
    if (world->options->mighist && world->options->skyline)
    {
      for (locus = 0; locus < world->loci; locus++)
        {
	  world->mighistloci[locus].eventbinsize = binsize;
	  world->mighistloci[locus].eventbins = (tetra **) mycalloc (world->numpop2, sizeof (tetra *));
	  world->mighistloci[locus].eventbinnum = (long *) mycalloc (world->numpop2, sizeof (long));
	  for (i = 0; i < world->numpop2; i++)
            {
	      world->mighistloci[locus].eventbinnum[i] = allocsize;
	      world->mighistloci[locus].eventbins[i] = (tetra *) mycalloc (allocsize, sizeof (tetra));
	      for(j=0;j<allocsize;j++)
		{
		  for(z=0;z<9;z++)
		    world->mighistloci[locus].eventbins[i][j][z] = 0.;
		}
            }
        }
    }
}

///
/// Destroy the skyline plot histogram container
void
destroy_expected_events (world_fmt * world)
{
    long locus, i;
    if (world->options->mighist && world->options->skyline)
    {
      for (locus = 0; locus < world->loci; locus++)
        {	  
	  for (i = 0; i < world->numpop2; i++)
            {
	      myfree(world->mighistloci[locus].eventbins[i]);
            }
	  myfree(world->mighistloci[locus].eventbins);
	  myfree(world->mighistloci[locus].eventbinnum);
        }
    }
}


void print_expected_values_list(FILE *file, long locus, tetra **eventbins, MYREAL eventbinsize, long *eventbinnum, long numpop)
{
  long i;
  long pop;
  long frompop;
  long topop;
  long numpop2 = numpop * numpop;
  MYREAL age;
  
  for(pop = 0; pop < numpop2; pop++)
    {  
      age = 0.;
      if(pop < numpop)
	{
	  fprintf(file,"\nLocus: %li   Parameter: %s_%li\n", locus+1, "Theta",pop+1);  
	}
      else
	{
	  m2mm(pop,numpop,&frompop,&topop);
	  fprintf(file,"\nLocus: %li   Parameter: %s_(%li,%li)\n", locus+1, "M", frompop+1, topop+1);  
	}
      fprintf(file,"Time        Parameter       Frequency of visit\n");
      fprintf(file,"----------------------------------------------\n");
      for(i = 0; i < eventbinnum[pop]; i++)
	{
	  age += eventbinsize;
	  if(eventbins[pop][i][1] < SMALL_VALUE)
	    continue;
	  if(eventbins[pop][i][0] < 0)
	    error("nono");
	  fprintf(file,"%10.10f  %10.10f     %10.10f\n", age, eventbins[pop][i][0], eventbins[pop][i][1]);
	}
    } 
}

#define SOME_ELEMENTS 100
void read_expected_values_fromfile(FILE *file, world_fmt *world)
{
  long locus;
  long pop;
  long i;
 // MYREAL age;
  char *input;
  char *inptr;
  tetra **eventbins;
  long *eventbinnum;
//xcode  MYREAL eventbinsize;
  long spacer;
  long * allocsize;
  allocsize = (long *) mycalloc(world->numpop2*world->loci,sizeof(long));
  for(i=0;i<world->numpop2*world->loci;i++)
    {
      allocsize[i]=1;
    }
  input = (char *) mycalloc(SUPERLINESIZE, sizeof(char));
  while(FGETS(input,SUPERLINESIZE,file) != EOF)
    {
      // grab the commentlines
      while(input[0]=='#')
	{
	  FGETS(input,LINESIZE,file);
	}
      // read the skylinefile
      if(input !=NULL)
	{
	  inptr = input;
	  locus      = atol(strsep(&inptr,"\t"))-1;
	  spacer = locus * world->numpop2;
	  if(locus == -1)
	    error("help");
	  pop       = atol(strsep(&inptr,"\t"))-1;
	  i          = atol(strsep(&inptr,"\t"))-1;
	  //xcode age        = atof(strsep(&inptr,"\t"));
      (void) strsep(&inptr,"\t");

	  eventbins =  world->mighistloci[locus].eventbins;
	  eventbinnum = world->mighistloci[locus].eventbinnum;
	  //xcode eventbinsize = world->mighistloci[locus].eventbinsize;
	  eventbinnum[pop] = i+1;
	  if(allocsize[spacer + pop] < eventbinnum[pop])
	    {
	      allocsize[spacer + pop] += SOME_ELEMENTS;
	      eventbins[pop] = (tetra *) myrealloc(eventbins[pop], allocsize[spacer + pop] * sizeof(tetra));
	    }

	  //xcode if(i==0)
	  //xcode   {
	  //xcode     eventbinsize = 2 * age;
	  //xcode   }
	  eventbins[pop][i][0] = atof(strsep(&inptr,"\t"));
	  eventbins[pop][i][1] = atof(strsep(&inptr,"\t"));
	  eventbins[pop][i][2] = atof(strsep(&inptr,"\t"));
	  eventbins[pop][i][3] = 0;
	  eventbins[pop][i][4] = atof(strsep(&inptr,"\t"));
	  eventbins[pop][i][5] = atof(strsep(&inptr,"\t"));
	}
    }
  myfree(allocsize);
}

void print_expected_values_tofile(FILE *file,  world_fmt *world)
{
  long i;
  long pop;
  long frompop;
  long topop;
  long locus;
  long sumloc;
  long numpop = world->numpop;
  long numpop2 = world->numpop2;
  MYREAL age;
  tetra **eventbins;
  long *eventbinnum;
  MYREAL eventbinsize;

  fprintf(file,"# Raw record of the skyline histogram for all parameters and all loci\n");  
  fprintf(file,"# The time interval is set to %f\n", world->options->eventbinsize);  
  fprintf(file,"# produced by the program %s (http://popgen.csit.fsu.edu/migrate.hml)\n",
	  MIGRATEVERSION);  
  fprintf(file,"# written by Peter Beerli 2006, Tallahassee,\n");
  fprintf(file,"# if you have problems with this file please email to beerli@fsu.edu\n");  
  fprintf(file,"#\n");
  fprintf(file,"# Order of the parameters:\n");
  fprintf(file,"# Parameter-number Parameter\n");
  for(pop=0;pop<numpop2;pop++)
    {
      if(pop < numpop)
	{
	  fprintf(file,"# %6li    %s_%li\n", pop+1, "Theta",pop+1);  
	}
      else
	{
	  m2mm(pop,numpop,&frompop,&topop);
	  fprintf(file,"# %6li    %s_(%li,%li)\n", pop+1, (world->options->usem ? "M" : "xNm"), frompop+1, topop+1);  
	}
    }
  fprintf(file,"#\n#----------------------------------------------------------------------------\n");
  fprintf(file,"# Locus Parameter-number Bin Age Parameter-value Parameter-Frequency \n");
  fprintf(file,"#        Standard-deviation Counts-per-bin Autocorrelation-per-bin\n");
  fprintf(file,"#----------------------------------------------------------------------------\n");
  fprintf(file,"# (*) values with -1 were NEVER visited\n");
  if(world->loci>1)
    {
      sumloc =1;
      fprintf(file,"# Locus %li is sum over all loci, when there are more than 1 locus\n", world->loci+1);
    }
  else
    {
      sumloc = 0;
    }

  for(locus=0; locus < world->loci + sumloc; locus++)
    {
      if(!world->data->skiploci[locus])
	{
	  eventbins =  world->mighistloci[locus].eventbins;
	  eventbinnum = world->mighistloci[locus].eventbinnum;
	  eventbinsize = world->mighistloci[locus].eventbinsize;
	  for(pop = 0; pop < numpop2; pop++)
	    {  
	      age = eventbinsize / 2.;
	      for(i = 0; i < eventbinnum[pop]; i++)
		{
		  fprintf(file,"%li\t%li\t%li\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n", 
			  locus+1, pop+1, i+1, age, eventbins[pop][i][0], eventbins[pop][i][1],eventbins[pop][i][2],eventbins[pop][i][4],eventbins[pop][i][5]);
		  age += eventbinsize;
		}
	    } 
	}
    }
}

void prepare_expected_values(world_fmt *world)
{
  MYREAL sum;
  long locus;
  long pop;
  long i;
  long numpop2 = world->numpop2;
  tetra **eventbins = NULL;
  tetra **eventbins_all;
  long *eventbinnum;
  long eventbinnum_allmax=0;
  float * suml;
  float *count;
  float rs;

  suml = (float *) mycalloc(world->numpop2,sizeof(float));
  count = (float *) mycalloc(world->numpop2,sizeof(float));

  for(locus=0; locus < world->loci; locus++)
    {
      if(!world->data->skiploci[locus])
	{
	  eventbins =  world->mighistloci[locus].eventbins;
	  eventbinnum = world->mighistloci[locus].eventbinnum;
	  for(pop = 0; pop < numpop2; pop++)
	    {  
	      sum = (MYREAL) 0.0;
	      	      
	      for(i = 0; i < eventbinnum[pop]; i++)
		{
		  if(eventbins[pop][i][1] <= 0.0)
		    {
		      eventbins[pop][i][0] = 0.0;
		      eventbins[pop][i][1] = 0.0;
		      eventbins[pop][i][2] = 0.0;
		      eventbins[pop][i][3] = 0.0;
		      eventbins[pop][i][4] = 0.0;
		      eventbins[pop][i][5] = 0.0;
		      eventbins[pop][i][6] = 0.0;
		      eventbins[pop][i][7] = 0.0;
		      eventbins[pop][i][8] = 0.0;
		    }
		  sum += eventbins[pop][i][1];
		}

	      for(i = 0; i < eventbinnum[pop]; i++)
		{
		  if(eventbins[pop][i][1] <= 0.0)
		    continue;
		  // calculate standard deviations of values (using weighted formula)
		  // and adjusting for correlated samples
		  eventbins[pop][i][5] /= eventbins[pop][i][2];//calculation of autocorrelation coefficient
		  rs = eventbins[pop][i][5];//autocorrelation
		  rs = eventbins[pop][i][4] * (1.-rs)/(1.+rs);// artificial value correcting for correlation 
		  // rs is the number of observation in bin
		  // standard deviation corrected with autocorrelation
		  // rs is number of obersvations the formula below
		  // assumes that eventbins[][][1] is the average weight and not the sumweight
		  // because originally to formula calls for *= rs/((rs-1) sumweight)
		  // here we do rs/((rs-1) rs * weightaverage)
		  eventbins[pop][i][2] *= 1./((rs-1.) * eventbins[pop][i][1]);

		  if (MYISNAN (eventbins[pop][i][5]))
		    {
		      eventbins[pop][i][5] = -999.;
		    }
		  if (eventbins[pop][i][2]<0. || MYISNAN (eventbins[pop][i][2]))
		    {
		      eventbins[pop][i][2]=-999.;
		      //warning("%i> the Variance for parameter %li in the skyline plots was negative -- this should not happen!",myID,i);
		    }
		  else
		    {
		      eventbins[pop][i][2] = sqrt(eventbins[pop][i][2])/sqrt(eventbins[pop][i][4]);
		    }
		  // calculate standard deviations of weights (this were unweighted)
		  // this seems only adjust for sample versus population but this was done before see up
		  eventbins[pop][i][3] *= eventbins[pop][i][4]/(eventbins[pop][i][4]-1.); 
		  eventbins[pop][i][3] = sqrt(eventbins[pop][i][3]);
		  eventbins[pop][i][1] /= sum; // calculate frequency: sum(weight_i)/sum(sum(weight_i)_j)
		}
	      // cutting off the tail of the skyline to avoid too many zeroes.
	      i = eventbinnum[pop]-1;
	      while(i>0 && eventbins[pop][i][2]<=0.0)
		    {
		      i--;
		    }
	      eventbinnum[pop] = i;
	      // calculating the highest bin number
	      if(eventbinnum[pop] > eventbinnum_allmax)
		eventbinnum_allmax = eventbinnum[pop];
	    }
	}
    }
  if(world->loci>1)
    {  
      if(world->mighistloci[world->loci].eventbins == NULL)
	world->mighistloci[world->loci].eventbins = (tetra **) mycalloc(world->numpop2,sizeof(tetra *));
      if(world->mighistloci[world->loci].eventbinnum == NULL)
      world->mighistloci[world->loci].eventbinnum = (long *) mycalloc(world->numpop2,sizeof(long));
      world->mighistloci[world->loci].eventbinsize = world->mighistloci[0].eventbinsize;
      for(pop=0; pop< world->numpop2 ; pop++)    
	{
	  world->mighistloci[world->loci].eventbins[pop] = (tetra *) mycalloc(eventbinnum_allmax,sizeof(tetra));
	  world->mighistloci[world->loci].eventbinnum[pop] = eventbinnum_allmax;  
	}
      eventbins_all = world->mighistloci[world->loci].eventbins;
      for (locus = 0; locus < world->loci; locus++)
	{
	  if(!world->data->skiploci[locus])
	    {
	      eventbinnum = world->mighistloci[locus].eventbinnum;
	      eventbins = world->mighistloci[locus].eventbins;
	      for(pop=0; pop< world->numpop2 ; pop++)
		{
		  for(i=0 ; i < eventbinnum[pop]; i++)
		    {
		      if(eventbins[pop][i][1] > 0.0)
			{
			  eventbins_all[pop][i][0] += eventbins[pop][i][0] * eventbins[pop][i][1];
			  eventbins_all[pop][i][1] += eventbins[pop][i][1];
			  eventbins_all[pop][i][2] += eventbins[pop][i][2] * eventbins[pop][i][1];
			  eventbins_all[pop][i][3] += eventbins[pop][i][3] * eventbins[pop][i][1];
			  eventbins_all[pop][i][4] += eventbins[pop][i][4];
			  eventbins_all[pop][i][5] += eventbins[pop][i][5] * eventbins[pop][i][1];
			}
		      suml[pop] += eventbins[pop][i][1];
		      count[pop] += eventbins[pop][i][4];
		    }
		}
	    }
	}
      for(pop=0; pop< world->numpop2 ; pop++)
	{
	  for(i=0 ; i < eventbinnum_allmax; i++)
	    {
	      if(eventbins_all[pop][i][1] > 0.0)
		{
		  eventbins_all[pop][i][0] /= eventbins_all[pop][i][1];//average
		  eventbins_all[pop][i][2] /= eventbins_all[pop][i][1];//average standard deviation
		  eventbins_all[pop][i][3] /= eventbins_all[pop][i][1];//average weight
		  eventbins_all[pop][i][4] /= world->loci;//average count:count[pop];
		  eventbins_all[pop][i][5] /= eventbins_all[pop][i][1];//average autocorrelation coefficient
		}
	    }
	}  
    } 
  myfree(suml);
  myfree(count);
}
  
void print_expected_values_title(FILE *file, boolean progress)
{
  if(progress)
    {      
      fprintf(file,"\n\nParameter changes over time\n");
      fprintf(file,"---------------------------\nSEE IN PDF FILE AND SKYLINE FILE\n");
    }
}

void print_expected_values(world_fmt * world, option_fmt *options)
{
  long locus;
  long sumloc = (world->loci > 1) ? 1 : 0;
  if(world->options->skyline)
    {
      if(world->options->datatype != 'g')
	{
	  // prepare skyline histogram for printing
	  prepare_expected_values(world);
	  // print skyline to file
	  print_expected_values_tofile(world->skylinefile, world);
	}
      else
	{
	  read_expected_values_fromfile(world->skylinefile,world);
	}
      // print title to screen
      print_expected_values_title(stdout, world->options->progress);
      // print title to ascii-outfile
      print_expected_values_title(world->outfile, TRUE);
      // print content to stdout and to ascii - outfile
      for(locus=0; locus < world->loci+sumloc; locus++)
      	{
	  if(world->options->verbose)
	    print_expected_values_list(stdout, locus, world->mighistloci[locus].eventbins, 
	    			       world->mighistloci[locus].eventbinsize, 
	    			       world->mighistloci[locus].eventbinnum, world->numpop);
	  print_expected_values_list(world->outfile, locus, world->mighistloci[locus].eventbins, 
	    			     world->mighistloci[locus].eventbinsize, 
	    			     world->mighistloci[locus].eventbinnum, world->numpop);
	}
#ifdef PRETTY
	  pdf_skyline_histogram(world->loci, world->numpop2,  world, FALSE);
	  if(world->options->bayes_infer)
	    pdf_skyline_histogram(world->loci, world->numpop2,  world, TRUE);
#endif
    }
}


void debug_skyline(world_fmt *world, char text[])
{
  long i;
  long pop;
  long locus;
  tetra ** eventbins ;
  long *eventbinnum;
  MYREAL eventbinsize;
  MYREAL age;
  FILE *file = stdout;
  long numpop2 = world->numpop2;
  fprintf(file,"#%i -------- %s ----------\n",myID,text);
  fprintf(file,"#\n#----------------------------------------------------------------------------\n");
  fprintf(file,"# Locus Parameter-number Bin Age Parameter-value(*) Parameter-Frequency(*)\n");
  fprintf(file,"#----------------------------------------------------------------------------\n");
  fprintf(file,"# (*) values with -1 were NEVER visited\n");
  for(locus=0; locus < world->loci; locus++)
    {
      if(!world->data->skiploci[locus])
	{
	  eventbins =  world->mighistloci[locus].eventbins;
	  eventbinnum = world->mighistloci[locus].eventbinnum;
	  eventbinsize = world->mighistloci[locus].eventbinsize;
	  for(pop = 0; pop < numpop2; pop++)
	    {  
	      age = eventbinsize / 2.;
	      for(i = 0; i < eventbinnum[pop]; i++)
		{
		  fprintf(file,"%li %li %li %10.10f %10.10f %10.10f %10.10f\n", 
			  locus+1, pop+1, i+1, age, eventbins[pop][i][0], eventbins[pop][i][1],eventbins[pop][i][0]/ eventbins[pop][i][1]);
		  age += eventbinsize;
		}
	    } 
	}
    }
  fprintf(file,"#%i >>>>>>>>>>>>>>> %s <<<<<<<<<<<<<<<<<END\n",myID,text);
}

