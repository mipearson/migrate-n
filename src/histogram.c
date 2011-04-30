// histogrammer for bayes histogram data
// takes the file bayesallfile and reads it into the histogram structure
// to use the calculate_hpd etc and also call the pretty printer functions
//
// (c) Peter Beerli 2007-2010
//
#include "definitions.h"
#include "migration.h"
#include "bayes.h"
#include "tools.h"
#include "sighandler.h"
#include "reporter.h"
#include "pretty.h"

extern int myID;
#ifdef ZNZ
void read_bayes_fromfile(znzFile mdimfile, world_fmt *world,option_fmt *options);
#else
void read_bayes_fromfile(FILE *mdimfile, world_fmt *world,option_fmt *options);
#endif

#ifdef ZNZ
void read_from_bayesmdim_minimal_info(znzFile mdimfile, world_fmt *world,option_fmt *options, data_fmt *data)
#else
void read_from_bayesmdim_minimal_info(FILE *mdimfile, world_fmt *world,option_fmt *options, data_fmt *data)
#endif
{
  char *input;
  //long nrep;
  long pop;
  long tmp;
  boolean done=FALSE;
  boolean recordedusem=TRUE;

  input = (char *) mycalloc(LINESIZE , sizeof(char));
#ifdef ZNZ
  while(done==FALSE && ZNZFGETS(input,LINESIZE,mdimfile) != EOF)
#else
  while(done==FALSE && FGETS(input,LINESIZE,mdimfile) != EOF)
#endif
    {
      if(input[0] == '#' && strstr(input,"begin"))
	{
#ifdef ZNZ
	  ZNZFGETS(input,LINESIZE,mdimfile);
#else
	  FGETS(input,LINESIZE,mdimfile);
#endif
	  printf("%i>>>>>>> read from bayesallfile <<<<<<<<<<<<<<<<<<<<<<\n",myID);
	  options->custm = (char *) myrealloc(options->custm, sizeof(char) * (strlen(input)+1));
	  strcpy(options->custm,input+9);
	  printf("%i> custom       = %s\n",myID, options->custm);
#ifdef ZNZ
	  ZNZFGETS(input,LINESIZE,mdimfile);
#else
	  FGETS(input,LINESIZE,mdimfile);
#endif	  
	  options->custm2 = (char *) myrealloc(options->custm2, sizeof(char) * (strlen(input)+1));
	  strcpy(options->custm2,input+9);
	  printf("%i> custom2      = %s\n",myID, options->custm2);
#ifdef ZNZ
	  ZNZFGETS(input,LINESIZE,mdimfile);
#else
	  FGETS(input,LINESIZE,mdimfile);
#endif
	  sscanf (input+3, "%li %li %li %li %li %i", &world->loci, &world->numpop,
		  &world->numpop2, &tmp, &options->replicatenum,&recordedusem);
	  printf("%i> loci         = %li\n",myID, world->loci);
	  printf("%i> numpop       = %li\n",myID, world->numpop);
	  printf("%i> numpop^2     = %li\n",myID, world->numpop2);
	  printf("%i> replicate    = %li\n",myID, tmp);
	  printf("%i> replicatenum = %li\n",myID, options->replicatenum);
	  printf("%i> use_M        = %li\n",myID, (long) recordedusem);
	  // fill some more...
	  data->numpop = world->numpop;
	  options->newpops_numpop = world->numpop;
	//xcode   nrep = options->replicatenum;
	//xcode   if (nrep == 0)
	//xcode     nrep = 1;
	  options->replicate = (boolean) tmp;
	  if(options->usem != recordedusem)
	    options->recordedusem = recordedusem;
        data->popnames = (char **) mymalloc (sizeof (char *) * world->numpop);
        for (pop = 0; pop < world->numpop; pop++)
	  {
            data->popnames[pop] = (char *) mycalloc (1, sizeof (char) * LINESIZE);
#ifdef ZNZ
	    ZNZFGETS(input,LINESIZE,mdimfile);
#else
	    FGETS(input,LINESIZE,mdimfile);
#endif
	    sscanf (input+3, "%s", data->popnames[pop]);
	    printf("%i> population = %s\n",myID, data->popnames[pop]);
	  }
	done=TRUE;
	}
    }
  options->muloci = data->loci = world->loci;
  data->skiploci =
    (boolean *) myrealloc (data->skiploci,
			   sizeof (boolean) * (data->loci + 1));
  memset (data->skiploci, 0, sizeof (boolean) * (data->loci + 1));
  data->numpop = world->numpop;
  printf("%i>>>>>>> end read from bayesallfile <<<<<<<<<<<<<<<<<<\n",myID);
  myfree(input);
}
			    

long get_fullbinsum(MYREAL *lowerbound, MYREAL *upperbound, world_fmt *world, option_fmt *options, long locus)
{
  long temp=0;
  long i;
  for(i = 0; i < world->numpop; i++)
    {
      temp += options->bayespriortheta->bins;
      lowerbound[i] = options->bayespriortheta->min;
      upperbound[i] = options->bayespriortheta->max;
    }
  for(i = world->numpop; i < world->numpop2; i++)
    {
      temp += options->bayespriorm->bins;
      lowerbound[i] = options->bayespriorm->min;
      upperbound[i] = options->bayespriorm->max;
    }
  for(i = world->numpop2; i < world->numpop2+world->bayes->mu * world->loci; i++)
    {
      temp += options->bayespriorrate->bins;
      lowerbound[i] = options->bayespriorrate->min;
      upperbound[i] = options->bayespriorrate->max;
    }
  return temp;
}


///
/// read bayesallfile from disk and creates all the needed parts to 
/// recreate the output and pdf output
#ifdef ZNZ
void read_bayes_fromfile(znzFile mdimfile, world_fmt *world,option_fmt *options)
#else
void read_bayes_fromfile(FILE *mdimfile, world_fmt *world,option_fmt *options)
#endif
{

  long *n = NULL;
  const long nn = world->numpop2 + world->bayes->mu * world->loci + 1;// One is for Log(Prob(Data|Model)                                    
  long j0;
  long j;
  //  char *input;
  long step;
  long locus;
  long numpop = world->numpop;
  long numpop2 = world->numpop2;
  long frompop;
  long topop;
 // long replicate;
 // long T ; // was world->treetimes->T, but this construct is freed, we do not need to
  // keep this values as it only used for the output into bayesallfile
 // MYREAL treelength;
  MYREAL post;
  MYREAL like;
//  MYREAL probg;
//  MYREAL prior;
  MYREAL *params;
#ifdef BFDEBUG
  long nnn=0;
  long t;
//  MYREAL lsum;
  long hc = world->options->heated_chains; 
#endif
  char *input;
  char *inptr;
  bayes_fmt * bayes = world->bayes;
  bayeshistogram_fmt *hist;
  long bin;
  MYREAL *delta = bayes->deltahist;
  //MYREAL *delta; 
  MYREAL *lowerbound;
  MYREAL *upperbound;
  MYREAL *autocorrelation;
  MYREAL *ess;
  //  long binwidth;
  long numbins = 0;
  long numbinsall = 0;
  //MYREAL m;
  boolean *done;
  boolean recordedusem = options->usem;
  const long np = world->numpop2 + world->bayes->mu;

  done = (boolean *) mycalloc(world->loci, sizeof(boolean));
  // params can be replaced by a single value for param and rate
  params = (MYREAL *) mycalloc(2+ numpop2 + world->bayes->mu, sizeof(MYREAL));
  autocorrelation = (MYREAL *) mycalloc(2 * world->loci * nn, sizeof(MYREAL));
  ess = autocorrelation + world->loci * nn;
  lowerbound = (MYREAL *) mycalloc(nn, sizeof(MYREAL));
  upperbound = (MYREAL *) mycalloc(nn, sizeof(MYREAL));
  //delta = (MYREAL *) mycalloc(nn-1, sizeof(MYREAL));
  n = (long *) mycalloc(numpop2 + world->bayes->mu * world->loci, sizeof(long));
  input = (char *) mycalloc(SUPERLINESIZE, sizeof(char));
#ifdef ZNZ
  unsigned long bytes = SUPERLINESIZE > ONEMEGABYTE ? SUPERLINESIZE : ONEMEGABYTE;
  znzbuffer(mdimfile,bytes);
  while(ZNZFGETS(input,SUPERLINESIZE,mdimfile) != EOF)
#else
  while(FGETS(input,SUPERLINESIZE,mdimfile) != EOF)
#endif
    {
      // grab the commentlines
      while(input[0]=='#' || input[0]=='S')
	{
#ifdef ZNZ
	  ZNZFGETS(input,SUPERLINESIZE,mdimfile);
#else
	  FGETS(input,SUPERLINESIZE,mdimfile);
#endif
	}
      // read the bayesallfile
      
      if(input !=NULL)
	{
	  inptr = input;
	  step       = atol(strsep(&inptr,"\t"));
	  locus      = atol(strsep(&inptr,"\t"))-1;
	  if(locus == -1)
	    error("help");
	//  replicate  = atol(strsep(&inptr,"\t"))-1;
      (void) strsep(&inptr,"\t");
	  post       = atof(strsep(&inptr,"\t"));
	  like       = atof(strsep(&inptr,"\t"));
	//  probg      = atof( strsep(&inptr,"\t"));
       (void) strsep(&inptr,"\t");
	//  prior      = atof( strsep(&inptr,"\t"));
       (void) strsep(&inptr,"\t");
	//  T = atol(strsep(&inptr,"\t"))+1;
       (void) strsep(&inptr,"\t");
	//  treelength = atof(strsep(&inptr,"\t"));
       (void) strsep(&inptr,"\t");
	  params[0] = post;
	  params[1] = like;

	  if(!done[locus])
	    {

	      done[locus] = TRUE;
	      // allocate the number of bins for the histogram
	      bayes->histogram[locus].binsum = get_fullbinsum(lowerbound, upperbound, world, options, locus);
	      bayes->histogram[locus].results = (MYREAL *) mycalloc(bayes->histogram[locus].binsum + 1, sizeof(MYREAL));
	      bayes->histogram[locus].set95 = (char *) mycalloc(bayes->histogram[locus].binsum* 2 + 2, sizeof(char));
	      bayes->histogram[locus].set50 = world->bayes->histogram[locus].set95 + bayes->histogram[locus].binsum + 1;
	      memset(world->bayes->histogram[locus].results, 0 , sizeof(MYREAL) * (world->bayes->histogram[locus].binsum)); 
   
	    }
	  hist = &bayes->histogram[locus];
	  numbinsall = 0;
	  for(j0=0;j0 < numpop2; j0++)
	    {
	      if(bayes->map[j0][1] == INVALID)
		continue;
	      else
		{
		  j = bayes->map[j0][1];
		}
	      if(j < j0)
		{
		  continue;
		}
	      else
		{
		  params[j+2] =  atof(strsep(&inptr,"\t"));
		  if(j>=numpop && options->usem!=recordedusem)
		    {
		      m2mm(j, numpop,&frompop,&topop);
		      if(recordedusem)// recorded bayesallfile was with M, new request is for 4Nm 
			{
			  params[j+2] *= params[topop+2];
			}
		      else //recorded bayesallfile was with xNm, new request is for M
			{
			  params[j+2] /= params[topop+2];
			}
		    }
		  n[j] += 1;
		  hist->means[j] += (params[j+2] - hist->means[j]) / n[j];
		}
	      numbinsall += hist->bins[j];
	      numbins = numbinsall - hist->bins[j];
	      bin = (long) ((params[j+2]-lowerbound[j]) / delta[j]);
	      hist->minima[j0] = lowerbound[j];
	      hist->maxima[j0] = upperbound[j];
	      hist->results[numbins + bin] += 1.;
	      bayes->histtotal[locus * np + j] += 1;
	    }
	  if(bayes->mu)
	    {
	      numbins = numbinsall;
	      params[j0+2] = atof(strsep(&inptr,"\t"));
	      n[j0+locus] += 1;
	      hist->means[j0] += (params[j0+2] - hist->means[j0]) / n[j0+locus];
	      bin = (long) ((params[j0+2]-lowerbound[j0]) / delta[j0]); 
	      hist->minima[j0] = lowerbound[j0];
	      hist->maxima[j0] = upperbound[j0];
	      hist->results[numbins + bin] += 1.;
	      bayes->histtotal[locus * np + j0] += 1;
	    }
#ifdef BFDEBUG
	  if(world->options->datatype == 'g')
	    {
	      for(t=0;t<hc;t++)
		{
		  world->bf[locus * hc + t] = atof(strsep(&inptr,"\t"));
		}
	      // dummy read of thermo sum up to this point
	      //lsum = atof(strsep(&inptr,"\t"));
          (void ) strsep(&inptr,"\t");
	      // harmonic mean: scaler contains the log value, hm contains 1.
	      world->hmscale[locus] = atof(strsep(&inptr,"\t"));
	      world->hm[locus] = 1.; 
	      //	      calculate_ess_frombayes (world, step, params, locus, &autocorrelation[locus*nn], &ess[locus*nn]);
	      calculate_ess_frombayes (world, step, params, locus, autocorrelation, ess);
	    }
#endif
	}
    }
#ifdef BFDEBUG
  if(world->options->datatype == 'g')
    {
      // reset the archiving machinery
      memset(world->auto_archive,0, sizeof(MYREAL) * 2 * world->numpop2+options->bayesmurates * world->loci + 1);
      nnn = 1;
      for(j=0;j<world->loci;j++)
	{
	  
	  for(t=0;t<nn; t++)
	    {
	      // onepass mean of autocorrelation
	      world->auto_archive[t] += (autocorrelation[j*nn + t] - world->auto_archive[t])/nnn;
	      // summing ess values
	      world->ess_archive[t] += ess[j*nn + t];
	      printf("j=%li t=%li %f\n", j, t, world->ess_archive[t]);
	    }
	  nnn++;
	}
    }
#endif /*BFDEBUG*/
  myfree(params);
  myfree(autocorrelation);
  myfree(lowerbound);
  myfree(upperbound);
  myfree(n);
  myfree(done);
  myfree(input);
}
