/* \File options.c */
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm   
 ------------------------------------------------------- 
 O P T I O N S   R O U T I N E S 
 
 creates options structures,
 reads options from parmfile if present
 
 prints options,
 and finally helps to destroy itself.
 
 Peter Beerli 1996-2006
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2009 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
$Id: options.c 1831 2011-03-20 19:03:13Z beerli $
-------------------------------------------------------*/
//#include <stdio.h>
#include <time.h>
#include "migration.h"
#include "sighandler.h"
#include <stdarg.h>
#include "fst.h"
#include "tools.h"
#include "migrate_mpi.h"
#include "bayes.h"
#include "menu.h"
#include "random.h"
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
#include "options.h"

extern char * generator;

/* parmfile parameter specifications and keywords */
//#define LONGLINESIZE 10000
#define NUMBOOL 31
#define BOOLTOKENS {"menu","interleaved","print-data","mixfile",\
  "moving-steps","freqs-from-data","useroldtree", \
  "autocorrelation", "simulation","plot", "weights",\
  "read-summary","write-summary","unused", "include-unknown", "print-fst",\
  "distfile","geo","gelman-convergence", "randomtree",\
      "fast-likelihood", "aic-modeltest", "use-M", "use-Nm","bayes-update", "bayes-allfile", "bayes-file", "divergence-times","tipdate-file","heated-swap","auto-tune"}
#define NUMNUMBER 64
#define NUMBERTOKENS {"ttratio","short-chains",\
 "short-steps","short-inc","long-chains",\
 "long-steps", "long-inc", "theta", \
 "nmlength","random-seed","migration","mutation",\
 "datatype", "categories","rates","prob-rates", \
 "micro-max", "micro-threshold", "delimiter","burn-in",\
 "infile", "outfile", "mathfile", "title", \
 "long-chain-epsilon","print-tree","progress","l-ratio",\
 "fst-type","profile","custom-migration","sumfile","short-sample",\
 "long-sample", "replicate","cpu","logfile", "seqerror-rate","uep", \
 "uep-rates","uep-bases", "mu-rates","heating","fluctuate", "resistance",\
 "bayes-updatefreq", "bayesfile","bayes-prior", "usertree", "bayes-posteriorbins",\
 "mig-histogram", "bayes-posteriormaxtype", "pdf-outfile",\
 "bayes-allfileinterval", "bayes-priors","skyline","rates-gamma", "bayes-proposals", \
      "generation-per-year", "mutationrate-per-year", "inheritance-scalars", "micro-submodel","random-subset", "population-relabel"};

// myID is a definition for the executing node (master or worker)
extern int myID;

fpos_t thePos;   // used to track pposition in the parmfile

/* prototypes ------------------------------------------- */
void set_usem_related (option_fmt * options);
//void create_options (option_fmt ** options);
//void init_options (option_fmt * options);
//void set_param (world_fmt * world, data_fmt * data, option_fmt * options,long locus);
//void set_profile_options (option_fmt * options);
//void print_menu_options (world_fmt * world, option_fmt * options,
//                         data_fmt * data);
//void print_options (FILE * file, world_fmt * world, option_fmt * options,
//                    data_fmt * data);
//void decide_plot (worldoption_fmt * options, long chain, long chains,
//                  char type);
//void destroy_options (option_fmt * options);

//void read_options_master (option_fmt * options);
//void read_options_worker (char **buffer, option_fmt * options);
///* private functions */
boolean booleancheck (option_fmt * options, char *var, char *value);
//long boolcheck (char ch);
boolean numbercheck (option_fmt * options, char *var, char *value);
//void reset_oneline (option_fmt * options, long position);
void reset_oneline (option_fmt * options, fpos_t * position);
void read_theta (option_fmt * options);
void read_mig (option_fmt * options);
#ifdef MPI
void read_theta_worker (char **buffer, option_fmt * options);
void read_mig_worker (char **buffer, option_fmt * options);
char skip_sspace (char **buffer);
void read_custom_migration_worker (char **buffer, option_fmt * options,
                               char *value, long customnumpop);
#endif
char skip_space (option_fmt * options);
//void read_custom_migration (FILE * file, option_fmt * options, char *value,
//                            long customnumpop);
//void synchronize_param (world_fmt * world, option_fmt * options);
//void resynchronize_param (world_fmt * world);
void specify_migration_type (option_fmt * options);
//void print_options_nomcmc (FILE * file, option_fmt * options,
//                               world_fmt * world);
//void destroy_options (option_fmt * options);
//void set_partmean_mig (long **mmparam, MYREAL *param, char *custm2, long migm,
//                       long numpop2);
//long scan_connect (char *custm2, long start, long stop, int check);
//void set_plot (option_fmt * options);
void set_plot_values (MYREAL **values, MYREAL plotrange[], long intervals, int type);
//void set_grid_param (world_fmt * world, long gridpoints);
void print_arbitrary_migration_table (FILE * file, world_fmt * world, option_fmt *options,
                                      data_fmt * data);
void print_distance_table (FILE * file, world_fmt * world,
                           option_fmt * options, data_fmt * data);

void fillup_custm (long len, world_fmt * world, option_fmt * options);
//

//long save_options_buffer (char **buffer, long *allocbufsize, option_fmt * options, data_fmt *data);
//void print_parm_delimiter(long *bufsize, char **buffer, long *allocbufsize);
//void print_parm_br(long *bufsize, char **buffer, long *allocbufsize);
//
void prior_consistency(prior_fmt *prior, int type);
//
///*======================================================*/
/// initialize filename
void init_filename (char **filename, char initstring[])
{
    *filename = (char *) mycalloc(STRSIZE, sizeof(char));
    strncpy(*filename, initstring, STRSIZE-1);
}

void
create_options (option_fmt ** options)
{
    (*options) = (option_fmt *) mycalloc (1, sizeof (option_fmt));
}

/// initialize options
void init_options (option_fmt * options)
{
    unsigned long i;
    unsigned long timeseed;
    /* General options --------------------------------------- */
    //
    // name length
    options->nmlength = DEFAULT_NMLENGTH;
    options->popnmlength = DEFAULT_POPNMLENGTH;
    options->allelenmlength = DEFAULT_ALLELENMLENGTH;
    //
    // custom migration matrix setup
    options->custm = (char *) mycalloc (1, sizeof (char) * 1000); /// \todo needs
    //options->custm2 is allocated later
    //options->symn = 0;
    //options->sym2n = 0;
    //options->zeron = 0;
    //options->constn = 0;
    //options->mmn = 0;
    //options->tmn = 0;
    
    /* input/output options ---------------------------------- */
    options->menu = TRUE;
    options->progress = TRUE;
    options->verbose = FALSE;
    options->writelog = FALSE;
    options->geo = FALSE;
    options->div = FALSE;
#ifdef UEP
    options->uep = FALSE;
    options->ueprate = 1.0;
    options->uepmu = 0.9999;
    options->uepnu = 0.0001;
    options->uepfreq0=0.5;
    options->uepfreq1=0.5;
    //  options->uep_last = FALSE;
#endif
    options->printdata = FALSE;
    options->printcov = FALSE;
    options->usertree = FALSE;
    options->usertreewithmig = FALSE;
    options->randomtree = TRUE;  //default changed May 29 2006 [convergence statistic]
    options->fastlike = FALSE;
    options->treeprint = _NONE;
    options->treeinc = 1;
    options->printfst = FALSE;
    options->fsttype = THETAVARIABLE;
    options->usem = TRUE;
    options->migvar = PLOT4NM;
    fst_type (options->fsttype);
    options->plot = FALSE;
    options->plotmethod = PLOTALL; /* outfile and mathematica file */
    options->plotvar = PLOT4NM;
    options->plotscale = PLOTSCALELOG;
    options->plotrange[0] = PLANESTART; /* start x axis */
    options->plotrange[1] = PLANEEND; /*end x axis */
    options->plotrange[2] = PLANESTART; /*start y axis */
    options->plotrange[3] = PLANEEND; /* end y axis */
    options->plotintervals = PLANEINTERVALS;
    options->simulation = FALSE;
    options->movingsteps = FALSE;
    options->acceptfreq = 0.0;
    options->mighist = FALSE;
    options->mighist_all = FALSE;
    options->mighist_counter = 0;
    options->mighist_increment = 1;
    options->skyline = FALSE;
    options->eventbinsize = 0.001;
    options->mixplot = FALSE;
    init_filename( &options->parmfilename, PARMFILE);
    init_filename( &options->infilename, INFILE);
    init_filename( &options->outfilename, OUTFILE);
#ifdef PRETTY
    init_filename( &options->pdfoutfilename, PDFOUTFILE);
#endif
#ifdef THERMOCHECK
    init_filename( &options->mixfilename, MIXFILE);
#endif
    init_filename( &options->logfilename, LOGFILE);
    init_filename( &options->mathfilename, MATHFILE);
    init_filename( &options->sumfilename, SUMFILE);
    init_filename( &options->treefilename, TREEFILE);
    init_filename( &options->utreefilename, UTREEFILE);
    init_filename( &options->catfilename, CATFILE);
    init_filename( &options->weightfilename, WEIGHTFILE);

    init_filename( &options->mighistfilename, MIGHISTFILE);
    init_filename( &options->skylinefilename, SKYLINEFILE);

    init_filename( &options->distfilename, DISTFILE);
    init_filename( &options->geofilename, GEOFILE);
    init_filename( &options->divfilename, DIVFILE);
    init_filename( &options->bootfilename, BOOTFILE);
    init_filename( &options->aicfilename, AICFILE);
    init_filename( &options->seedfilename, SEEDFILE);    
#ifdef UEP
    init_filename( &options->uepfilename, UEPFILE);
#endif
    init_filename( &options->bayesfilename, BAYESFILE);
    init_filename( &options->bayesmdimfilename, BAYESMDIMFILE);
    options->use_compressed=1;
    init_filename( &options->datefilename, TIPDATEFILE);
    // allocation of prior parameter settings
    options->bayespriortheta = (prior_fmt *) mycalloc (3, sizeof (prior_fmt));
    options->bayespriorm = options->bayespriortheta + 1;
    options->bayespriorrate = options->bayespriortheta + 2;
    // prior theta
    options->bayespriortheta->min = SMALLEST_THETA;
    options->bayespriortheta->mean = DNA_GUESS_THETA;
    options->bayespriortheta->max = 10. * DNA_GUESS_THETA;
    options->bayespriortheta->bins = BAYESNUMBIN;
    options->bayespriortheta->delta = (options->bayespriortheta->max - options->bayespriortheta->min)/10.;
    // prior M
    options->bayespriorm->min = SMALLEST_MIGRATION;
    options->bayespriorm->mean = DNA_GUESS_MIG;
    options->bayespriorm->max = 10. * DNA_GUESS_MIG;
    options->bayespriorm->bins = BAYESNUMBIN;
    options->bayespriorm->delta = (options->bayespriorm->max - options->bayespriorm->min)/10.;
    // prior rate
    options->bayespriorrate->min = 1.0;
    options->bayespriorrate->mean = 1.0;
    options->bayespriorrate->max = 1.0;
    options->bayespriorrate->bins = BAYESNUMBIN;
    options->bayespriorrate->delta = (options->bayespriorrate->max - options->bayespriorrate->min)/10.;
    for(i=0;i<PRIOR_SIZE;i++)
      options->slice_sampling[i] = TRUE;
#ifdef PRETTY
    options->bayespretty = PRETTY_P100;
#endif
    strcpy (options->title, "\0");
    options->lratio = (lratio_fmt *) mycalloc (1, sizeof (lratio_fmt));
    options->lratio->alloccounter = 1;
    options->lratio->data =
        (lr_data_fmt *) mycalloc (1, sizeof (lr_data_fmt) *options->lratio->alloccounter);
    options->lratio->data[0].value1 =
        (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
    options->lratio->data[0].value2 =
        (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
    options->lratio->data[0].connect =
        (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
    options->aic = FALSE;
    options->fast_aic = FALSE;
    strcpy(options->aictype,"0m");
    options->profile = ALL;
    options->profilemethod = 'p';
    options->df = 1;
    options->qdprofile = FALSE;
    options->printprofsummary = TRUE;
    options->printprofile = TRUE;
    if (options->usem)
        options->profileparamtype = options->plotvar = options->migvar = PLOTM;
    else
        options->profileparamtype = options->plotvar = options->migvar = PLOT4NM;
    /* data options ------------------------------------------ */
    options->datatype = 's';
    options->include_unknown = FALSE;
    options->migration_model = MATRIX;
    options->thetag = (MYREAL *) mycalloc (1, sizeof (MYREAL) * NUMPOP);
    options->mg = (MYREAL *) mycalloc (1, sizeof (MYREAL) * NUMPOP);
    options->gamma = FALSE;
    options->alphavalue = START_ALPHA;
    options->murates = FALSE;
    options->bayesmurates = FALSE;
#ifdef LONGSUM
    options->fluctuate=FALSE;
    options->flucrates = (MYREAL *) mycalloc(NUMPOP * 3, sizeof(MYREAL));
#endif
    options->mu_rates = NULL;
    options->inheritance_scalars = NULL;
    options->newpops = NULL;
    options->prioralone=FALSE;
    /* EP data */
    options->dlm = '\0';
    /* microsat data */
    options->micro_threshold = MICRO_THRESHOLD;
    options->micro_stepnum = MAX_MICROSTEPNUM;
    options->msat_option = SINGLESTEP;
    options->msat_tuning[0] = 0.0;
    options->msat_tuning[1] = 0.5;
    /*sequence data */
    options->seqerror = 0.0;
    options->interleaved = FALSE;
    options->ttratio = (MYREAL *) mycalloc (1, sizeof (MYREAL) * 2);
    options->ttratio[0] = 2.0;
    options->freqsfrom = TRUE;
    options->categs = ONECATEG;
    options->rate = (MYREAL *) mycalloc (1, sizeof (MYREAL));
    options->rate[0] = 1.0;
    options->rcategs = 1;
    options->rrate = (MYREAL *) mycalloc (1, sizeof (MYREAL));
    options->probcat = (MYREAL *) mycalloc (1, sizeof (MYREAL));
    options->autocorr = FALSE;
    options->rrate[0] = 1.0;
    options->probcat[0] = 1.0;

    options->probsum = 0.0;

    options->lambda = 1.0;
    options->weights = FALSE;
    /* random number options --------------------------------- */
    options->autoseed = AUTO;
    options->autoseed = AUTO;
    //#ifndef MAC
    timeseed = (unsigned long) time (NULL) / 4;
    //#else
    //timeseed = (unsigned long) clock () / 4;
    //#endif
    options->inseed = (long) timeseed + 1;
    /* mcmc options ------------------------------------------ */
    options->thetaguess = FST;
    options->migrguess = FST;
    for (i = 0; i < NUMPOP; i++)
    {
        options->thetag[i] = 1.0;
        options->mg[i] = 1.0;
        options->custm[i] = '*';
    }
    options->custm[i] = '\0';
    options->custm2 = (char *) mycalloc (LINESIZE, sizeof (char));
    strncpy (options->custm2, options->custm, i);
    options->numthetag = options->nummg = 0;
    options->bayes_infer = FALSE;
#ifdef INTEGRATEDLIKE
    options->integrated_like = TRUE;
#else
    options->integrated_like = FALSE;
#endif
    options->updateratio = HALF;
    options->bayesmdiminterval = 1;
    options->schains = 10;
    options->sincrement = 100;
    options->ssteps = 500;
    options->lchains = 1;
    options->lincrement = 100;
    options->lsteps = 5000;
    options->burn_in = BURNINPERIOD;
    options->burnin_autostop = ' ';
    //
    // heating options
    options->heating = 0;  /* no heating */
    options->adaptiveheat=NOTADAPTIVE;
    options->heating_interval = 0;
    options->heatedswap_off=FALSE;
    options->heat[0] = COLD;
    options->heat[1] = WARM;
    options->heat[2] = HOT;
    options->heat[3] = VERYHOT;
    options->heated_chains = 1;
    //
    // obscure options to accept more different trees
    options->lcepsilon = LONGCHAINEPSILON;
    options->gelman = TRUE;  // changed default May 29 2006
    options->gelmanpairs = FALSE;
    options->pluschain = PLUSCHAIN;
    //
    // replication options
    options->replicate = FALSE;
    options->replicatenum = 0;
    options->gridpoints = 0;
    /* genealogy summary options----------------------------------- */
    options->readsum = FALSE;
    options->writesum = FALSE;
    /*threading over loci */
    options->cpu = 1;
    //
    // fatal attraction to zero resistance
    options->minmigsumstat = MINMIGSUMSTAT;
    //
    options->mutationrate_year = (MYREAL*) mycalloc(1, sizeof(MYREAL));
    options->mutationrate_year_numalloc = 1;
    options->generation_year = 1.0;
    options->randomsubset = 0;
    options->inheritance_scalars = (MYREAL*) mycalloc(1, sizeof(MYREAL));
    options->inheritance_scalars[0] = 1.0;
    options->inheritance_scalars_numalloc  = 1;
    options->newpops = (long*) mycalloc(1, sizeof(long));
    options->newpops[0] = 1;
    options->newpops_numalloc  = 1;

    options->slice_sticksizes = (MYREAL*) mycalloc(1, sizeof(MYREAL));
    options->slice_sticksizes[0] = 1.0;
    options->heatedswap_off=FALSE;
    options->has_autotune = FALSE;
    options->autotune=0.33;
}

void
read_options_master (option_fmt * options)
{
  const unsigned long linecpmax = LINESIZE -1; 

    long counter = 0;
    //long position = 0;
    // to make windows reading unix parmfile happy
    //    fpos_t thePos;

    char varvalue[LONGLINESIZE];
    char parmvar[LONGLINESIZE];
    char input[LONGLINESIZE];
    char *p, *tmp;
  
    options->parmfile = fopen(options->parmfilename, "r");
    if (options->parmfile)
    {
        counter = 0;
        //position = ftell (options->parmfile);
	fgetpos(options->parmfile, &thePos);
	FGETS (input, LINESIZE, options->parmfile);// read first set of ###
	FGETS (input, LINESIZE, options->parmfile);// read "# Parmfile for Migrate"
	if(strncmp(input,"# Parmfile for Migrate", 20))
	  {
	    usererror("This file is not a parameter file (parmfile) for migrate\nSyntax for calling from the commandline is\nmigrate-n parmfile\n and not (!!!!)\nmigrate-n datafile\n");
	    exit(-1);
	  }
	fprintf(stdout,"Reading parmfile \"%s\"....\n",options->parmfilename);
        while (FGETS (input, LINESIZE, options->parmfile) != EOF)
        {
            counter++;
            if ((input[0] == '#') || isspace ((int) input[0])
                    || input[0] == ';')
            {
	      //position = ftell (options->parmfile);
	      fgetpos(options->parmfile, &thePos);
	      continue;
            }
            else
            {
                if (!(isalnum ((int) input[0]) || strchr ("{}", input[0])))
                {
                    usererror ("The parmfile contains an error on line %li\n",
                               counter);
                }
            }
            if ((p = strchr (input, '#')) != NULL)
                *p = '\n';
            if (!strncmp (input, "end", 3))
                break;
            tmp = strtok (input, "=");
            if (tmp != NULL)
                strncpy (parmvar, tmp, linecpmax);
            else
            {
                if(input[0]!='\0')
                    fprintf (stderr,
                             "WARNING: error in parmfile on line %li with %s\n",
                             counter, input);
                continue;
            }

	    // DEBUG check that we really read stuff
#ifdef DEBUG
	    fprintf(stdout,"%s\n", parmvar);
#endif
            if (!strncmp (parmvar, "theta", 5))
            {
	      //reset_oneline (options, position);
	      reset_oneline (options, &thePos);
                read_theta (options);
                //position = ftell (options->parmfile);
		fgetpos(options->parmfile, &thePos);
                continue;
            }
            if (!strncmp (parmvar, "migration", 5))
            {
	      //reset_oneline (options, position);
	      reset_oneline (options, &thePos);
                read_mig (options);
                //position = ftell (options->parmfile);
		fgetpos(options->parmfile, &thePos);
                continue;
            }
            tmp = strtok (NULL, "\n");
            if (tmp != NULL)
                strncpy (varvalue, tmp, LINESIZE - 1);
            if (!booleancheck (options, parmvar, varvalue))
            {
                if (!numbercheck (options, parmvar, varvalue))
                {
                    warning ("Inappropiate entry in parmfile: %s ignored\n",
                             input);
                }
            }
            //position = ftell (options->parmfile);
	    fgetpos(options->parmfile, &thePos);
        }
    }
}

#ifdef MPI
void
read_options_worker (char **buffer, option_fmt * options)
{
  const unsigned long linecpmax = LINESIZE - 1;

    long counter = 0;
    char *position;
    char varvalue[LONGLINESIZE];
    char parmvar[LONGLINESIZE];
    char input[LONGLINESIZE];
    char *p, *tmp;

    options->buffer = buffer;
    options->parmfile = NULL;
    if (strlen (*buffer) > 0)
    {
        counter = 0;
        position = *buffer;
        while (sgets (input, LINESIZE, buffer) != NULL)
        {
            counter++;
            if ((input[0] == '#') || isspace ((int) input[0])
                    || input[0] == ';')
            {
                position = *buffer;
                continue;
            }
            else
            {
                if (!(isalnum ((int) input[0]) || strchr ("{}", input[0])))
                {
                    usererror ("The parmfile contains an error on line %li\n",
                               counter);
                }
            }
            if ((p = strchr (input, '#')) != NULL)
                *p = '\n';
            if (!strncmp (input, "end", 3))
                break;
            tmp = strtok (input, "=");
            if (tmp != NULL)
                strncpy (parmvar, tmp, linecpmax);
            else
            {
                if(input[0]!='\0')
                    fprintf (stderr,
                             "WARNING: error in parmfile on line %li with %s\n",
                             counter, input);
                continue;
            }
            if (!strncmp (parmvar, "theta", 5))
            {
                *buffer = position;
                read_theta_worker (buffer, options);
                position = *buffer;
                continue;
            }
            if (!strncmp (parmvar, "migration", 5))
            {
                *buffer = position;
                read_mig_worker (buffer, options);
                position = *buffer;
                continue;
            }
            tmp = strtok (NULL, "\n");
            if (tmp != NULL)
                strncpy (varvalue, tmp, LINESIZE - 1);
            if (!booleancheck (options, parmvar, varvalue))
            {
                if (!numbercheck (options, parmvar, varvalue))
                {
                    warning ("Inappropiate entry in parmfile: %s ignored\n",
                             input);
                }
            }
            position = *buffer;
        }
    }
    if(options->bayes_infer)
        options->schains=0;
#ifdef DEBUG_MPI
    printf("%i> custm=%s custm2=%s\n",myID,options->custm, options->custm2);
    printf("%i> generation/year=%f mutationrate/year=%f\n",myID, options->generation_year, options->mutationrate_year[0]);
#endif

}
#endif

void
print_menu_options (world_fmt * world, option_fmt * options, data_fmt * data)
{ 
    if (options->numpop > data->numpop)
        usererror ("Inconsistency between your Menu/Parmfile and your datafile\n \
                   Check the number of populations!\n");
    if (options->progress)
    {
        print_options (stdout, world, options, data);
        if (options->writelog)
            print_options (options->logfile, world, options, data);
    }
}

///
/// prints the value of minimum value of the prior distribution
char * show_priormin(char *tmp, prior_fmt *prior, int priortype)
{
  sprintf(tmp, "%f ", prior->min);
  return tmp;
}

///
/// prints the value of mean value of the prior distribution
/// or for the gamma proposal the alpha value, and the ,multipler proposal the multiplicator
char * show_priormean(char *tmp, prior_fmt *prior, int priortype)
{
  sprintf(tmp, "%f ", prior->mean);
  return tmp;
}

/// prints the value of maximum value of the prior distribution
char * show_priormax(char *tmp, prior_fmt *prior, int priortype)
{
  sprintf(tmp, "%f ", prior->max);
  return tmp;
}

///
/// prints the value of delta value of the prior distribution
/// for some distribution this has no meaning
char * show_priordelta(char *tmp, prior_fmt *prior, int priortype)
{
  switch(priortype)
    {
    case EXPPRIOR:
      sprintf(tmp, "-"); break;
    default:
      sprintf(tmp, "%f ", prior->delta);
      break;
    }
  return tmp;
}

///
/// prints the bins used for the prior/posterior distribution
char * show_priorbins(char *tmp, prior_fmt *prior, int priortype)
{
  sprintf(tmp, "%li ", prior->bins);
  return tmp;
}


void
print_options (FILE * file, world_fmt * world, option_fmt * options,
               data_fmt * data)
{
  int count;
    long i, j, tt;
    char mytext[LINESIZE];
    char mytext1[LINESIZE];
    char mytext2[LINESIZE];
    char mytext3[LINESIZE];
    char mytext4[LINESIZE];
    char seedgen[LINESIZE], spacer[LINESIZE];
    char paramtgen[LINESIZE], parammgen[LINESIZE];
    if (options->datatype != 'g')
    {
        switch ((short) options->autoseed)
        {
        case AUTO:
            strcpy (seedgen, "with internal timer");
            strcpy (spacer, "  ");
            break;
        case NOAUTOSELF:
            strcpy (seedgen, "from parmfile");
            strcpy (spacer, "      ");
            break;
        case NOAUTO:
            strcpy (seedgen, "from seedfile");
            strcpy (spacer, "      ");
            break;
        default:
            strcpy (seedgen, "ERROR");
            strcpy (spacer, " ");
            break;
        }
        switch (options->thetaguess)
        {
        case OWN:
            strcpy (paramtgen, "from guessed values");
            break;
        case FST:
            strcpy (paramtgen, "from an FST-calculation");
            break;
            //    case PARAMGRID:
            //      strcpy (paramtgen, "GRID values around a center");
            //      break;
        case NRANDOMESTIMATE:
            strcpy (paramtgen, "RANDOM start value from N(mean,std) or U(min,max)");
            break;
        case URANDOMESTIMATE:
            strcpy (paramtgen, "RANDOM start value from U(min,max)");
            break;
        default:
            strcpy (paramtgen, "ERROR");
            break;
        }
        switch (options->migrguess)
        {
        case OWN:
            strcpy (parammgen, "from guessed values");
            break;
        case FST:
            strcpy (parammgen, "from the FST-calculation");
            break;
            //    case PARAMGRID:
            //      strcpy (parammgen, "GRID values around a center");
            //      break;
        case NRANDOMESTIMATE:
            strcpy (parammgen, "RANDOM start value from N(mean,std)");
            break;
        case URANDOMESTIMATE:
            strcpy (parammgen, "RANDOM start value from U(min,max)");
            break;
        default:
            strcpy (parammgen, "ERROR");
            break;
        }
    }
    fprintf (file, "Options in use:\n");
    fprintf (file, "---------------\n\n");
    if(options->bayes_infer)
    {
        fprintf (file, "Analysis strategy is                                    Bayesian\n\n");
	fprintf (file, "Proposal distribution:\n");
	fprintf (file, "Parameter group          Proposal type\n");
	fprintf (file, "-----------------------  -------------------\n");
	fprintf (file, "Population size (Theta) %20s\n", 
		 is_proposaltype(options->slice_sampling[THETAPRIOR]));
	if(world->numpop > 1)
	  {
	    fprintf (file, "Migration rate  %7.7s %20s\n",
		     options->usem ? "(M)" : "(xNm)",
		     is_proposaltype(options->slice_sampling[MIGPRIOR]));
	  }
	if(world->bayes->mu)
	  {
	    fprintf (file, "Mutation rate modifier  %20s\n",
		     is_proposaltype(options->slice_sampling[RATEPRIOR]));
	  }
	fprintf(file,"\n\n");
	fprintf (file, "Prior distribution:\n");
	fprintf (file, "Parameter group          Prior type   Minimum    Mean(*)    Maximum    Delta\n");
	fprintf (file, "-----------------------  ------------ ---------- ---------- ---------- ----------\n");
	fprintf (file, "Population size (Theta) %12s %10.10s %10.10s %10.10s %10.10s\n",
		 is_priortype(options->bayesprior[THETAPRIOR]),
		 show_priormin(mytext1, options->bayespriortheta,options->bayesprior[THETAPRIOR]),
		 show_priormean(mytext2, options->bayespriortheta,options->bayesprior[THETAPRIOR]),
		 show_priormax(mytext3, options->bayespriortheta,options->bayesprior[THETAPRIOR]),
		 show_priordelta(mytext4, options->bayespriortheta,options->bayesprior[THETAPRIOR]));
	if(world->numpop > 1)
	  {
	    fprintf (file, "Migration rate  %7.7s %12s %10.10s %10.10s %10.10s %10.10s\n",
		     options->usem ? "(M)" : "(xNm)",
		     is_priortype(options->bayesprior[MIGPRIOR]),
		     show_priormin(mytext1,options->bayespriorm,options->bayesprior[MIGPRIOR]),
		     show_priormean(mytext2,options->bayespriorm,options->bayesprior[MIGPRIOR]),
		     show_priormax(mytext3,options->bayespriorm,options->bayesprior[MIGPRIOR]),
		     show_priordelta(mytext4,options->bayespriorm,options->bayesprior[MIGPRIOR]));
	  }
	if(world->bayes->mu)
	  {
	    fprintf (file, "Mutation rate modifier  %12s %10.10s %10.10s %10.10s %10.10s\n",
		     is_priortype(options->bayesprior[RATEPRIOR]),
		     show_priormin(mytext1,options->bayespriorrate,options->bayesprior[RATEPRIOR]),
		     show_priormean(mytext2,options->bayespriorrate,options->bayesprior[RATEPRIOR]),
		     show_priormax(mytext3,options->bayespriorrate,options->bayesprior[RATEPRIOR]),
		     show_priordelta(mytext4,options->bayespriorrate,options->bayesprior[RATEPRIOR]));
	  }
	fprintf(file,"\n\n\n");
    }
    else
    {
        fprintf (file, "Analysis strategy is                          Maximum likelihood\n\n\n\n");
    }

    switch (options->datatype)
    {
    case 'a':
        fprintf (file, "Datatype: Allelic data\n");
        fprintf (file, "Missing data is %s\n",options->include_unknown ? "included" : "not included");
        break;
    case 'b':
        fprintf (file, "Datatype: Microsatellite data [Brownian motion]\n");
        fprintf (file, "Missing data is %s\n",options->include_unknown ? "included" : "not included");
        break;
    case 'm':
      if(options->msat_option==SINGLESTEP)
        fprintf (file, "Datatype: Microsatellite data [Singlestep model]\n");
      else
        fprintf (file, "Datatype: Microsatellite data [Multistep model (Tune=%f, P_increase=%f)]\n",options->msat_tuning[0], options->msat_tuning[1]);
        fprintf (file, "Missing data is %s\n",options->include_unknown ? "included" : "not included");
        break;
    case 's':
        fprintf (file, "Datatype: DNA sequence data\n");
        break;
    case 'n':
        fprintf (file, "Datatype: Single nucleotide polymorphism data\n");
        break;
    case 'h':
        fprintf (file, "Datatype: Hapmap Single nucleotide polymorphism data\n");
        break;
    case 'u':
        fprintf (file,
                 "Datatype: Single nucleotide polymorphism data (PANEL)\n");
        break;
    case 'f':
        fprintf (file, "Datatype: Ancestral state method\n");
        break;
    case 'g':
        fprintf (file, "Datatype: Genealogy summary of an older run\n");
        break;
    }

    count=0;
    
    fprintf (file, "\nInheritance scalers in use for Thetas (specified scalars=%li)\n",options->inheritance_scalars_numalloc);

    for(i=0;i<data->loci;i++)
      {
	if(i<options->inheritance_scalars_numalloc)
	  fprintf (file, "%2.2f ", options->inheritance_scalars[i]);
	else
	  fprintf (file, "%2.2f ", options->inheritance_scalars[options->inheritance_scalars_numalloc-1]);
	if(count++ > 7)
	  {
	    fprintf(file,"\n");
	    count=0;
	  }
	else 
	  {
	    count++;
	  }
      }
    fprintf (file, "\n[Each Theta uses the (true) ineritance scalar of the first locus as a reference]\n\n");

    if(options->randomsubset > 0)
      {
        fprintf (file, "\nData set was subsampled: used a random sample of size: %li\n\n", options->randomsubset);
      }

    if (options->datatype != 'g')
    {
      fprintf (file, "\n%-80s\n", generator);
#ifndef QUASIRANDOM      
        fprintf (file, "Random number seed (%s)%s%20li\n", seedgen, " ",
                 options->saveseed);
#endif
        fprintf (file, "\nStart parameters:\n");
	fprintf(file,  "   First genealogy was started using a %s\n",
		options->usertree ? "user tree" : (options->randomtree ? "random tree" : 
						   (options->dist ? "tree from userdistances" : 
						    "UPGMA-tree")));
	fprintf(file,  "   Theta values were generated ");
        fprintf (file, " %s\n", paramtgen);
        if (options->thetaguess == OWN)
        {
            fprintf (file, "   Theta = ");
            for (i = 0; i < options->numthetag - 1; i++)
            {
                fprintf (file, "%.5f,", options->thetag[i]);
            }
            fprintf (file, "%.5f\n", options->thetag[i]);
        }
        fprintf (file, "   M values were generated %s\n", parammgen);
        if (options->migrguess == OWN)
        {
            tt = 0;
            if (options->usem)
                fprintf (file, "   M-matrix: ");
            else
                fprintf (file, "   4Nm-matrix: ");
            if (options->nummg == 1)
            {
                fprintf (file, "%5.2f [all are the same]\n", options->mg[0]);
            }
            else
            {
                for (i = 0; i < world->numpop; i++)
                {
                    for (j = 0; j < world->numpop; j++)
                    {
                        if (i != j)
			  {
                            fprintf (file, "%5.2f ", options->mg[tt]);
			    if(tt<options->nummg)
			      tt++;
			  }
                        else
                            fprintf (file, "----- ");
                        if (j > 10)
                            fprintf (file, "\n                ");
                    }
                    fprintf (file, "\n               ");
                }
                fprintf (file, "\n");
            }
        }
    }
    fprintf(file,"\n");
    print_arbitrary_migration_table (file, world, options, data);
    print_distance_table (file, world, options, data);
    fprintf(file,"\n");
    // mutation related material
    if (options->gamma)
    {
        fprintf (file, "Mutation rate among loci is Gamma-distributed\n");
        fprintf (file, "Initial scale parameter alpha = %f\n",
                 options->alphavalue);
        if (options->custm[world->numpop2] == 'c')
            fprintf (file, "and is constant [will not be estimated]\n");
    }
    else /*Gamma*/
    {
      if (options->murates && world->loci > 1)
        {
	  fprintf (file, "Mutation rate among loci is varying with\n");
	  fprintf (file, "   Rates per locus: ");
	  for (i = 0; i < world->loci - 1; i++)
            {
	      fprintf (file, "%.3f, ", options->mu_rates[i]);
	      if ((i + 1) % 6 == 0)
		fprintf (file, "\n                    ");
            }
	  fprintf (file, "%.3f\n", options->mu_rates[i]);
	  if(options->murates_fromdata)
	    fprintf (file, "[Estimated from the data using the Watterson estimator (ignoring migration)]");
	  else
	    fprintf (file, "[User defined]");
        }
      else
	fprintf (file, "Mutation rate is constant %s\n", world->loci > 1 ?
		 "for all loci" : "");
    }
#ifdef UEP
    if (options->uep)
      {
        fprintf (file, "0/1 polymorphism analysis, with 0/1 data in file:%s\n",
                 options->uepfilename);
        fprintf (file, "     with forward mutation rate %f*mu\n",
                 options->uepmu);
        fprintf (file, "     with back    mutation rate %f*mu\n",
                 options->uepnu);
        fprintf (file, "     with base frequencies \"0\"=%f and \"1\"=%f\n",
                 options->uepfreq0,options->uepfreq1);
      }
#endif /*UEP*/
    if (options->datatype != 'g')
      {
        fprintf (file, "\nMarkov chain settings:\n");
	if(!options->bayes_infer)
	  {
	    fprintf (file, "   Short chains (short-chains):         %20li\n",
		     options->schains);
	    fprintf (file, "      Trees sampled (short-inc*samples):%20li\n",
		     options->sincrement * options->ssteps);
	    fprintf (file, "      Trees recorded (short-sample):    %20li\n",
		     options->ssteps);
	  }
        fprintf (file, "   Long chains (long-chains):           %20li\n",
                 options->lchains);
	if(options->bayes_infer)
	  {
	    if(options->replicate)
	      {
		fprintf (file, "      Steps sampled (inc*samples*rep):  %20li\n",
			 options->lincrement * options->lsteps * options->replicatenum);
		fprintf (file, "      Steps recorded (sample*rep):      %20li\n",
			 options->lsteps*options->replicatenum);
	      }
	    else
	      {
		fprintf (file, "      Steps sampled (long-inc*samples): %20li\n",
			 options->lincrement * options->lsteps);
		fprintf (file, "      Steps recorded (long-sample):     %20li\n",
			 options->lsteps);
	      }
	  }
	else
	  {
	    fprintf (file, "      Trees sampled (long-inc*samples): %20li\n",
		     options->lincrement * options->lsteps);
	    fprintf (file, "      Trees recorded (long-sample):     %20li\n",
		     options->lsteps);
	  }
        if (options->replicate)
        {
            if (options->replicatenum == 0)
                fprintf (file, "   Averaging over long chains\n");
            else
	      {
		if(!options->bayes_infer)
		  fprintf (file, "   Averaging over replicates:           %20li\n",
			   options->replicatenum);
		else
		  fprintf (file, "   Combining over replicates:           %20li\n",
			   options->replicatenum);
	      }
        }
        if (options->heating > 0)
        {
            fprintf (file,
                     "   %s heating scheme\n      %li chains with %s temperatures\n      ",
                     options->adaptiveheat!=NOTADAPTIVE ? (options->adaptiveheat==STANDARD ? "Adaptive_Standard" : "Bounded_adaptive") : "Static", options->heated_chains, options->adaptiveheat ? "start" : "" );
            for (i = 0; i < options->heated_chains - 1; i++)
                fprintf (file, "%5.2f,", options->heat[i]);
	    if(!options->heatedswap_off)
	      fprintf (file, "%5.2f\n      Swapping interval is %li\n",
		       options->heat[i], options->heating_interval);
	    else
	      fprintf (file, "%5.2f\n      No swapping\n",
		       options->heat[i]);
        }
        if (options->movingsteps)
        {
            fprintf (file, "   Forcing at least this percentage of new genealogies:%6.2f\n",
                     (MYREAL) options->acceptfreq);
        }
        if (options->burn_in > 0)
        {
	  sprintf(mytext,"%c%li", options->burnin_autostop, options->burn_in);
	  fprintf (file, "   Number of discard trees per chain:   %20.20s\n", mytext);
        }
        if (options->lcepsilon < LONGCHAINEPSILON)
        {
            fprintf (file, "   Parameter-likelihood epsilon:        %20.5f\n",
                     options->lcepsilon);
        }
    }
    fprintf (file, "\nPrint options:\n");
    if (options->datatype != 'g')
    {
        fprintf (file, "   Data file: %46.46s\n", options->infilename);
#ifdef PRETTY
        fprintf (file, "   Output file (ASCII text): %31.31s\n", options->outfilename);
        fprintf (file, "   Output file (PDF):        %31.31s\n", options->pdfoutfilename);
#else
        fprintf (file, "   Output file: %44.44s\n", options->outfilename);
#endif
        if(options->bayes_infer)
        {
            fprintf (file,   "   Posterior distribution: %33.33s\n", options->bayesfilename);
	    if(options->has_bayesmdimfile)
	      fprintf (file, "   All values of Post.Dist:%33.33s\n", options->bayesmdimfilename);
        }
        fprintf (file, "   Print data: %45.45s\n",
                 options->printdata ? "Yes" : "No");
        switch (options->treeprint)
        {
        case _NONE:
            fprintf (file, "   Print genealogies: %38.38s\n", "No");
            break;
        case ALL:
            fprintf (file, "   Print genealogies: %38.38s\n", "Yes, all");
            break;
        case LASTCHAIN:
            fprintf (file, "   Print genealogies: %32.32s%li\n",
                     "Yes, only those in last chain, every ", options->treeinc);
            break;
        case BEST:
            fprintf (file, "   Print genealogies: %38.38s\n",
                     "Yes, only the best");
            break;
        }
    }
    if (options->plot)
    {
        switch (options->plotmethod)
        {
        case PLOTALL:
            sprintf (mytext, "Yes, to outfile and %s", options->mathfilename);
            break;
        default:
            strcpy (mytext, "Yes, to outfile");
            break;
        }
        fprintf (file, "   Plot data: %-46.46s\n", mytext);
        fprintf (file,
                 "              Parameter: %s, Scale: %s, Intervals: %li\n",
                 options->plotvar == PLOT4NM ? "{Theta, 4Nm}" : "{Theta, M}",
                 options->plotscale == PLOTSCALELOG ? "Log10" : "Standard",
                 options->plotintervals);
        fprintf (file, "              Ranges: X-%5.5s: %f - %f\n",
                 options->plotvar == PLOT4NM ? "4Nm" : "M",
                 options->plotrange[0], options->plotrange[1]);
        fprintf (file, "              Ranges: Y-%5.5s: %f - %f\n", "Theta",
                 options->plotrange[2], options->plotrange[3]);
    }
    else
    {
        fprintf (file, "   Plot data: %-46.46s\n", "No");
    }

    if (options->mighist)
    {
      if(options->mighist_all)
	sprintf(mytext,"Yes: All events");
      else
	sprintf(mytext,"Yes: Migration events");
      fprintf (file,
	       "   Frequency histogram of events  %26.26s\n", mytext);
      fprintf (file,
	       "   Time of events are saved in file %24.24s\n",
	       options->mighistfilename);
    }
    if (options->mighist && options->skyline)
    {
      sprintf(mytext,"%3s","Yes");
      fprintf (file,
	       "   Histogram of the parameter values through time %10s\n", mytext);
      fprintf (file,
	       "   Parameters through time are saved in file %15.15s\n",
	       options->skylinefilename);
    }
    if(!options->bayes_infer)
      {
	switch (options->profile)
	  {
	  case _NONE:
	    strcpy (mytext, "No");
	    break;
	  case ALL:
	    strcpy (mytext, "Yes, tables and summary");
	    break;
	  case TABLES:
	    strcpy (mytext, "Yes, tables");
	    break;
	  case SUMMARY:
	    strcpy (mytext, "Yes, summary");
	    break;
	  }
	fprintf (file, "   Profile likelihood: %-36.36s\n", mytext);
	if (options->profile != _NONE)
	  {
	    switch (options->profilemethod)
	      {
	      case 'p':
		fprintf (file, "             Percentile method\n");
		break;
	      case 'q':
		fprintf (file, "             Quick method\n");
		break;
	      case 'f':
		fprintf (file, "             Fast method\n");
		break;
	      case 'd':
		fprintf (file, "             Discrete method\n");
		break;
	      case 's':
		fprintf (file, "             Spline method\n");
		break;
	      default:
		fprintf (file, "             UNKOWN method????\n");
		break;
	      }
	    fprintf (file, "             with df=%li and for Theta and %s\n\n\n\n",
		     options->df, options->profileparamtype ? "M=m/mu" : "4Nm");
	  }
      }
}

/// \brief sets the theta parameters for random start parameter setting
/// 
/// Sets the startparameter to a random number derived from a uniform distribution
/// with minimum options->thetag[0] and maximum in options->thetag[1]
void set_theta_urandomstart(world_fmt *world, option_fmt *options)
{
    long i;
    long pos=0;
    char temp[SUPERLINESIZE];
    if(world->options->progress)
      pos = sprintf(temp,"Random start Theta values:");
    for (i = 0; i < world->numpop; i++)
    {
        world->param0[i] = options->thetag[0] + RANDUM() * (options->thetag[1]-options->thetag[0]);
	if(world->options->progress)
	  pos += sprintf(temp+pos," %f",world->param0[i]);
    }
    if(world->options->progress)
      FPRINTF(stdout,"[%3i] %s\n",myID, temp);
    if(world->options->writelog)
      FPRINTF(world->options->logfile,"[%3i] %s\n",myID, temp);
}

/// \brief sets the theta parameters for random start parameter setting
/// 
/// Sets the startparameter to a random number derived from a normal distribution
/// with men options->thetag[0] and standard deviation options-.thetag[1]
void set_theta_nrandomstart(world_fmt *world, option_fmt *options)
{
    long i;
    long pos=0;
    char temp[SUPERLINESIZE];
    if(world->options->progress)
      pos = sprintf(temp,"Random start Theta values:");
    for (i = 0; i < world->numpop; i++)
    {
        world->param0[i] = rannor (options->thetag[0], options->thetag[1]);
        while (world->param0[i] < 0)
	  {
            world->param0[i] =
                rannor (options->thetag[0], options->thetag[1]);
	  }
	if(world->options->progress)
	  pos += sprintf(temp+pos," %f",world->param0[i]);
    }
    if(world->options->progress)
      FPRINTF(stdout,"%i> %s\n",myID, temp);
    if(world->options->writelog)
      FPRINTF(world->options->logfile,"%i> %s\n",myID, temp);
}


/// \brief sets the theta parameters for OWN start parameter setting
/// 
/// Sets the startparameter to a user-defined value
void set_theta_ownstart(world_fmt *world, option_fmt *options)
{
    long i, ii;
    for (i = 0; i < world->numpop; i++)
    {
        if (i < options->numthetag - 1)
            ii = i;
        else
            ii = options->numthetag - 1;
        if (options->thetag[ii] == 0.0)
            world->param0[i] = SMALLEST_THETA;
        else
        {
            world->param0[i] = options->thetag[ii];
        }
    }
}



/// \brief sets the theta parameters to a start parameter using an FST value
/// 
/// Sets the startparameter to a value from the FST calculation
void set_theta_fststart(world_fmt *world, option_fmt *options, long locus)
{
    long i;
    for (i = 0; i < world->numpop; i++)
    {
        if (world->fstparam[locus][i] > SMALLEST_THETA)
        {
            if (world->fstparam[locus][i] > 100)
                world->param0[i] = 1.0;
            else
                world->param0[i] = world->fstparam[locus][i];
        }
        else
        {
            if (strchr (SEQUENCETYPES, options->datatype))
            {
                world->param0[i] = 0.01;
            }
            else
            {
                world->param0[i] = 1.0;
            }
        }
    }
}


/// \brief sets the migration parameters for random start parameter setting
/// 
/// Sets the startparameter to a random number derived from a uniform distribution
/// with minimum options->mg[0] and maximum in options->mg[1]
void set_mig_urandomstart(world_fmt *world, option_fmt *options)
{
    long i;
    long pos=0;
    char temp[SUPERLINESIZE];
    if(world->options->progress)
      pos = sprintf(temp,"Random start M values:");

    for (i = world->numpop; i < world->numpop2; i++)
    {
        world->param0[i] =  options->mg[0] + RANDUM() * (options->mg[1]-options->mg[0]);
	if(world->options->progress)
	  pos += sprintf(temp+pos," %f",world->param0[i]);
    }
    if(world->options->progress)
      FPRINTF(stdout,"%i> %s\n",myID, temp);
    if(world->options->writelog)
      FPRINTF(world->options->logfile,"%i> %s\n",myID, temp);
}

/// \brief sets the migration parameters for random start parameter setting
/// 
/// Sets the startparameter to a random number derived from a normal distribution
/// with men options->mg[0] and standard deviation options->mg[1]
void set_mig_nrandomstart(world_fmt *world, option_fmt *options)
{
    long i;
    long pos=0;
    char temp[SUPERLINESIZE];
    if(world->options->progress)
      pos = sprintf(temp,"Random start M values:");

    for (i = world->numpop; i < world->numpop2; i++)
    {
        world->param0[i] = rannor (options->mg[0], options->mg[1]);
        while (world->param0[i] <= 0.0)
	  {
            world->param0[i] = rannor (options->mg[0], options->mg[1]);
	  }
	if(world->options->progress)
	  pos += sprintf(temp+pos," %f",world->param0[i]);
    }
    if(world->options->progress)
      FPRINTF(stdout,"%i> %s\n",myID, temp);
    if(world->options->writelog)
      FPRINTF(world->options->logfile,"%i> %s\n",myID, temp);
}

/// \brief sets the migration parameters for OWN start parameter setting
/// 
/// Sets the startparameter to a user-defined value
void set_mig_ownstart(world_fmt *world, option_fmt *options, data_fmt *data)
{
    long i;
    long j;
    long ii;
    long iitest;
    long numpop = world->numpop;
    MYREAL tmp;
    for (i = 0; i < numpop; i++)
    {
        for (j = 0; j < numpop; j++)
        {
            if (i == j)
                continue;
            if ((iitest = mm2m (j, i, numpop) - numpop) < options->nummg)
                ii = iitest;
            else
	      {
                ii = options->nummg - 1;
	      }
            if (options->usem)
                tmp = options->mg[ii];
            else
                tmp = options->mg[ii] / world->param0[i];
	    iitest += numpop;
            if (options->geo)
            {
                world->param0[iitest] = data->geo[iitest] * tmp;
            }
            else
            {
                world->param0[iitest] = tmp;
            }
        }
    }
}

/// \brief sets the migration parameters to a start parameter using an FST value
/// 
/// Sets the startparameter to a value from the FST calculation
void set_mig_fststart(world_fmt *world, option_fmt *options, long locus)
{
    long i;
        
    for (i = world->numpop; i < world->numpop2; i++)
    {
        if (world->fstparam[locus][i] > 0)
        {
            if (world->fstparam[locus][i] > 100)
            {
                world->param0[i] =
                1.0 / world->param0[(i - world->numpop) /
                    (world->numpop)];
                if (world->param0[i] > 10000)
                    world->param0[i] = 10000;
            }
            else
                world->param0[i] = world->fstparam[locus][i];
        }
        else
        {
            world->param0[i] =
            1.0 / world->param0[(i - world->numpop) / (world->numpop)];
            if (world->param0[i] > 10000)
                world->param0[i] = 10000;
        }
        
    }
}


/// \brief set the starting parameters in main structure world from options
///
/// Set the starting parameters in main structure world from options
/// - Random starting parameters:\n
///   will use an average and a standard deviation to generate a starting parameter
/// - Own starting parameters
/// - Starting parameters derived from FST
/// - check that values from FST are not ridicoulously off and that Bayes start is in min/max bounds
/// Once all parameters are set they are synchronized using the custom-migration matrix
/// that will set some values to zero or to the same value. The synchronization will
/// override user-error, but does provide little control about the outcome
/// error control is only done later when user can check the values in the logfile
void
set_param (world_fmt * world, data_fmt * data, option_fmt * options,
           long locus)
{
    // set THETA
    switch (options->thetaguess)
    {
    case NRANDOMESTIMATE:
        set_theta_nrandomstart(world, options);
        break;
        
    case URANDOMESTIMATE:
        set_theta_urandomstart(world, options);
        break;
        
    case OWN:
        set_theta_ownstart(world,options);
        break;

    case FST:
    default:
        set_theta_fststart(world,options, locus);
        break;
    }
    // set MIGRATION
    switch (options->migrguess)
    {
    case NRANDOMESTIMATE:
        set_mig_nrandomstart(world,options);
        break;

    case URANDOMESTIMATE:
        set_mig_urandomstart(world,options);
        break;

    case OWN:
        set_mig_ownstart(world,options, data);
        break;
    case SLATKIN:
    case FST:
    default:
        set_mig_fststart(world,options, locus);
        break;
    }
    if(options->bayes_infer)
    {
        bayes_check_and_fix_param(world,options);
    }
    synchronize_param (world, options);
    if (options->gamma)
        world->options->alphavalue = options->alphavalue;
    //perhaps to come
    //if(options->thetaguess==PARAMGRID && options->migrguess==PARAMGRID)
    //set_grid_param(world,options->gridpoints);
}

/// \todo not yet implemented, is this a good feature?
void
set_grid_param (world_fmt * world, long gridpoints)
{
    static long which = 0;
    static long z = 0;

    static MYREAL level = 0.1;
    static MYREAL bottom = -2.302585093; // => level =  0.1
    MYREAL top = 2.302585093; // => level = 10
    MYREAL len = top - bottom;
    MYREAL diff = len / (gridpoints - 1.);

    if (z >= gridpoints)
    {
        z = 0;
        which++;
    }
    level = EXP (bottom + z * diff);
    world->param0[which] = world->param0[which] * level;
    z++;

}

/// \brief synchronize parameters
void
synchronize_param (world_fmt * world, option_fmt * options)
{
    char type;
    boolean found = FALSE;
    long i, z = 0, zz = 0, zzz = 0, len;
    long ii, jj;
    long ns = 0, ss = 0, ss2 = 0, xs = 0, migm = 0;
    boolean allsymmig = FALSE;
    boolean allsamesize = FALSE;
    boolean partmig = FALSE;
    MYREAL summ;
    MYREAL nsum = 0;
    len = (long) strlen (world->options->custm2);
    // not needed!
    //    world->options->custm2 =
    //    (char *) myrealloc (world->options->custm2,
    //                      sizeof (char) * (world->numpop2 + 2));
    if(world->options->bayes_infer)
        world->bayes->custm2 = world->options->custm2;
    if (len < world->numpop2)
    {
        fillup_custm (len, world, options);
    }
    migm = scan_connect (world->options->custm2, 0, world->numpop, 'm');
    world->options->tmn = migm; // are there mean theta values?
    migm = scan_connect (world->options->custm2, world->numpop, world->numpop2, 'm');
    world->options->mmn = migm;
    for (i = 0; i < world->numpop2; i++)
    {
        if (!(allsymmig && allsamesize))
        {
            type = world->options->custm2[i];
            switch (type)
            {
            case '*':
                xs++;
                break;
            case 's':  // M is symmetric
                if (i >= world->numpop)
                {
                    z = i;
                    m2mm (z, world->numpop, &ii, &jj);
                    zz = mm2m (jj, ii, world->numpop);
                    world->options->symparam = (twin_fmt *) myrealloc
                                               (world->options->symparam, sizeof (twin_fmt) * (ss + 2));
                    world->options->symparam[ss][0] = zz;
                    world->options->symparam[ss++][1] = z;
                    world->options->symn = ss;
                    summ = (world->param0[z] + world->param0[zz]) / 2.;
                    world->param0[zz] = world->param0[z] = summ;
                }
                break;
            case 'S':  // 4Nm is symmetric, not completely
                // implemented yet, -> derivatives.c
                if (i >= world->numpop)
                {
                    z = i;
                    m2mm (z, world->numpop, &ii, &jj);
                    zz = mm2m (jj, ii, world->numpop);
                    zzz = 0;
                    found = FALSE;
                    while (zzz < ss2)
                    {
                        if (world->options->sym2param[zzz][1] == zz)
                            found = TRUE;
                        zzz++;
                    }
                    if (found)
                        break;
                    world->options->sym2param = (quad_fmt *)
                                                 myrealloc (world->options->sym2param,
                                                         sizeof (quad_fmt) * (ss2 + 2));
                    world->options->sym2param[ss2][0] = zz;
                    world->options->sym2param[ss2][1] = z;
                    world->options->sym2param[ss2][2] = ii;
                    world->options->sym2param[ss2++][3] = jj;
                    world->options->sym2n = ss2;
                    summ = (world->param0[z] * world->param0[jj] +
                            world->param0[zz] * world->param0[ii]) / 2.;
                    world->param0[z] = summ / world->param0[jj];
                    world->param0[zz] = summ / world->param0[ii];
                }
                break;
            case 'C':
            case 'c':
                world->options->constparam = (long *) myrealloc
                                             (world->options->constparam, sizeof (long) * (ns + 2));
                world->options->constparam[ns++] = i;
                world->options->constn = ns;
                break;
            case '0':
                z = i;
                world->param0[z] = 0;
                world->options->zeroparam = (long *) myrealloc
                                            (world->options->zeroparam, sizeof (long) * (ns + 2));
                world->options->zeroparam[ns++] = i;
                world->options->zeron = ns;
                break;
            case 'm':
                summ = 0;
                if (i < world->numpop) /*theta */
                {
                    if (!allsamesize)
                    {
                        nsum = 0;
                        allsamesize = TRUE;
                        for (z = 0; z < world->numpop; z++)
                        {
                            if (world->options->custm2[z] == 'm')
                            {
                                nsum++;
                                summ += world->param0[z];
                            }
                        }
                        summ /= nsum;
                        for (z = 0; z < world->numpop; z++)
                            if (world->options->custm2[z] == 'm')
                                world->param0[z] = summ;
                    }
                }
                else  /* migration */
                {
                    summ = 0;
                    if (!partmig)
                    {
                        nsum = 0;
                        partmig = TRUE;
                        for (z = world->numpop; z < world->numpop2; z++)
                        {
                            if (world->options->custm2[z] == 'm')
                            {
                                nsum++;
                                summ += world->param0[z];
                            }
                        }
                        summ /= nsum;
                        for (z = world->numpop; z < world->numpop2; z++)
                            if (world->options->custm2[z] == 'm')
                                world->param0[z] = summ;
                    }
                }
                break;
            case 'M':
                summ = 0;
                if (i < world->numpop) /*theta */
                {
                    if (!allsamesize)
                    {
                        nsum = 0;
                        allsamesize = TRUE;
                        for (z = 0; z < world->numpop; z++)
                        {
                            if (world->options->custm2[z] == 'm')
                            {
                                nsum++;
                                summ += world->param0[z];
                            }
                        }
                        summ /= nsum;
                        for (z = 0; z < world->numpop; z++)
                            if (world->options->custm2[z] == 'm')
                                world->param0[z] = summ;
                    }
                }
                else  /* migration */
                {
                    summ = 0;
                    if (!partmig)
                    {
                        nsum = 0;
                        partmig = TRUE;
                        for (z = world->numpop; z < world->numpop2; z++)
                        {
                            if (world->options->custm2[z] == 'M')
                            {
                                nsum++;
				m2mm (z, world->numpop, &ii, &jj);
                                summ += world->param0[z] * world->param0[jj];
                            }
                        }
                        summ /= nsum;
                        for (z = world->numpop; z < world->numpop2; z++)
			  {
                            if (world->options->custm2[z] == 'M')
			      {
				m2mm (z, world->numpop, &ii, &jj);
                                world->param0[z] = summ/world->param0[jj];
			      }
			  }
                    }
                }
                break;
            default:
                warning ("The migration connection matrix is misspecified\nOnly these are allowed\n* s S m M 0 (=zero)\nSupplied value was: %c\ncustom-migration=%s\n", type, world->options->custm);
		error("[Program exits]");
            }
        }
    }
    //--------gamma stuff
    if (world->options->gamma)
    {
        if (world->options->custm2[world->numpop2] == 'c')
        {
            world->options->constparam = (long *) myrealloc
                                         (world->options->constparam, sizeof (long) * (ns + 2));
            world->options->constparam[ns++] = i;
            world->options->constn = ns;
        }
        else
        {
            world->options->custm2[world->numpop2] = '*';
            world->options->custm2[world->numpop2 + 1] = '\0';
        }
    }
}


void
resynchronize_param (world_fmt * world)
{
    char type;
    long i, j = 0, z = 0, zz = 0;
    //long len;
    long ii, jj;
    long ns = 0, ss = 0, ss2 = 0, xs = 0, migm = 0;
    boolean allsymmig = FALSE;
    boolean allsamesize = FALSE;
    boolean partmig = FALSE;
    MYREAL summ;
    MYREAL nsum = 0;
 	//xcode      len = (long) strlen (world->options->custm2);
    world->options->symn = 0;
    world->options->sym2n = 0;
    world->options->zeron = 0;
    world->options->constn = 0;
    world->options->mmn = 0;
    migm = scan_connect (world->options->custm2, 0, world->numpop, 'm');
    world->options->tmn = migm; // are there mean theta values?
    migm = scan_connect (world->options->custm2, world->numpop, world->numpop2, 'm');
    world->options->mmn = migm;
    for (i = 0; i < world->numpop2; i++)
    {
        if (!(allsymmig && allsamesize))
        {
            type = world->options->custm2[i];
            switch (type)
            {
            case '*':
                xs++;
                break;
            case 's':  // M is symmetric
                if (i >= world->numpop)
                {
                    z = i;
                    m2mm (z, world->numpop, &ii, &jj);
                    zz = mm2m (jj, ii, world->numpop);
                    world->options->symparam = (twin_fmt *) myrealloc
                                               (world->options->symparam, sizeof (twin_fmt) * (ss + 2));
                    world->options->symparam[ss][0] = zz;
                    world->options->symparam[ss++][1] = z;
                    world->options->symn = ss;
                    summ = (world->param0[z] + world->param0[zz]) / 2.;
                    world->param0[zz] = world->param0[z] = summ;
                }
                break;
            case 'S':  // 4Nm is symmetric, not completely
                // implemented yet, -> derivatives.c
                if (i >= world->numpop)
                {
                    z = i;
                    m2mm (z, world->numpop, &ii, &jj);
                    zz = mm2m (jj, ii, world->numpop);
                    world->options->sym2param = (quad_fmt *)
                                                 myrealloc (world->options->sym2param,
                                                         sizeof (quad_fmt) * (ss2 + 2));
                    world->options->sym2param[ss2][0] = zz;
                    world->options->sym2param[ss2][1] = z;
                    world->options->sym2param[ss2][2] = i;
                    world->options->sym2param[ss2++][3] = j;
                    world->options->sym2n = ss2;
                    summ = (world->param0[z] * world->param0[i] +
                            world->param0[zz] * world->param0[j]) / 2.;
                    world->param0[z] = summ / world->param0[i];
                    world->param0[zz] = summ / world->param0[j];
                }
                break;
            case 'C':
            case 'c':
                world->options->constparam = (long *) myrealloc
                                             (world->options->constparam, sizeof (long) * (ns + 2));
                world->options->constparam[ns++] = i;
                world->options->constn = ns;
                break;
            case '0':
                z = i;
                world->param0[z] = 0;
                world->options->zeroparam = (long *) myrealloc
                                            (world->options->zeroparam, sizeof (long) * (ns + 2));
                world->options->zeroparam[ns++] = i;
                world->options->zeron = ns;
                break;
            case 'm':
                summ = 0;
                if (i < world->numpop) /*theta */
                {
                    if (!allsamesize)
                    {
                        allsamesize = TRUE;
                        nsum = 0;
                        for (z = 0; z < world->numpop; z++)
                        {
                            if (world->options->custm2[z] == 'm')
                            {
                                nsum++;
                                summ += world->param0[z];
                            }
                        }
                        summ /= nsum;
                        for (z = 0; z < world->numpop; z++)
                            if (world->options->custm2[z] == 'm')
                                world->param0[z] = summ;
                    }
                }
                else  /* migration */
                {
                    summ = 0;
                    if (!partmig)
                    {
                        nsum = 0;
                        partmig = TRUE;
                        for (z = world->numpop; z < world->numpop2; z++)
                        {
                            if (world->options->custm2[z] == 'm')
                            {
                                nsum++;
                                summ += world->param0[z];
                            }
                        }
                        summ /= nsum;
                        for (z = world->numpop; z < world->numpop2; z++)
                            if (world->options->custm2[z] == 'm')
                                world->param0[z] = summ;
                    }
                }
                break;
            default:
                error ("no defaults allowed in resynchronize_param()\n");
            }
        }
    }
    //--------gamma stuff
    if (world->options->gamma)
    {
        if (world->options->custm2[world->numpop2] == 'c')
        {
            world->options->constparam = (long *) myrealloc
                                         (world->options->constparam, sizeof (long) * (ns + 2));
            world->options->constparam[ns++] = i;
            world->options->constn = ns;
        }
        else
        {
            world->options->custm2[world->numpop2] = '*';
            world->options->custm2[world->numpop2 + 1] = '\0';
        }
    }
}


long
scan_connect (char *custm2, long start, long stop, int check)
{
    long i, summ = 0;
    for (i = start; i < stop; i++)
    {
        
        if (check == custm2[i])
            summ++;
    }
    return summ;
}

void
set_partmean_mig (long **mmparam, MYREAL *param, char *custm2, long migm,
                  long numpop2)
{
    long i, z = 0;
    MYREAL summ = 0;
    long start = (long) sqrt ((MYREAL) numpop2);
    (*mmparam) = (long *) myrealloc ((*mmparam), sizeof (long) * (migm + 2));

    for (i = start; i < numpop2; i++)
    {
        if (custm2[i] == 'm')
        {
            summ += param[i];
            (*mmparam)[z++] = i;
        }
    }
    summ /= migm;
    for (i = start; i < numpop2; i++)
    {
        if (custm2[i] == 'm')
            param[i] = summ;
    }
}

void free_options_filenames(option_fmt * options)
{
    myfree( options->parmfilename);
    myfree( options->infilename);
    myfree( options->outfilename);
#ifdef PRETTY
    myfree( options->pdfoutfilename);
#endif
    myfree( options->logfilename);
    myfree( options->mathfilename);
    myfree( options->sumfilename);
    myfree( options->treefilename);
    myfree( options->utreefilename);
    myfree( options->catfilename);
    myfree( options->weightfilename);
    myfree( options->mighistfilename);
    myfree( options->skylinefilename);
    myfree( options->distfilename);
    myfree( options->geofilename);
    myfree( options->divfilename);
    myfree( options->bootfilename);
    myfree( options->aicfilename);
    myfree( options->seedfilename);
    
#ifdef UEP
    myfree( options->uepfilename);
#endif
    myfree( options->bayesfilename);
    myfree( options->bayesmdimfilename);
    myfree( options->datefilename);
}

void
destroy_options (option_fmt * options)
{
    myfree(options->thetag);
    myfree(options->mg);
    myfree(options->mu_rates);
    myfree(options->inheritance_scalars);
    myfree(options->newpops);
    myfree(options->ttratio);
    myfree(options->rate);
    myfree(options->probcat);
    myfree(options->rrate);

    while (--options->lratio->alloccounter >= 0)
    {
        myfree(options->lratio->data[options->lratio->alloccounter].value1);
        myfree(options->lratio->data[options->lratio->alloccounter].value2);
        myfree(options->lratio->data[options->lratio->alloccounter].connect);
    }
    myfree(options->lratio->data);
    myfree(options->lratio);
    myfree(options->custm);
    myfree(options->custm2);
    if (options->plot)
    {
        myfree(options->plotxvalues);
        myfree(options->plotyvalues);
    }
    myfree(options->bayespriortheta);
    myfree(options->mutationrate_year);
    free_options_filenames(options);
    myfree(options);
}

void
decide_plot (worldoption_fmt * options, long chain, long chains, char type)
{
    if (options->plot && (chain >= chains - 1) && (type == 'l'))
        options->plotnow = TRUE;
    else
        options->plotnow = FALSE;
}

void
set_plot (option_fmt * options)
{
    long intervals = options->plotintervals;
    MYREAL prangex[2];
    MYREAL prangey[2];
    if (!options->plot)
        return;
    prangex[0] = options->plotrange[0];
    prangex[1] = options->plotrange[1];
    prangey[0] = options->plotrange[2];
    prangey[1] = options->plotrange[3];


    options->plotxvalues = (MYREAL *) mycalloc (1, sizeof (MYREAL) * intervals);
    options->plotyvalues = (MYREAL *) mycalloc (1, sizeof (MYREAL) * intervals);
    set_plot_values (&options->plotxvalues, prangex, options->plotintervals,
                     options->plotscale);
    set_plot_values (&options->plotyvalues, prangey, options->plotintervals,
                     options->plotscale);
}

void
set_plot_values (MYREAL **values, MYREAL plotrange[], long intervals,
                 int type)
{
    long i;
    MYREAL diff = 0;
    MYREAL logstart = 0;
    (*values)[0] = plotrange[0];
    (*values)[intervals - 1] = plotrange[1];
    if (type == PLOTSCALELOG)
    {
        logstart = log10 ((*values)[0]);
        diff =
            (log10 ((*values)[intervals - 1]) - logstart) / (MYREAL) (intervals -
                    1);
        for (i = 1; i < intervals - 1; i++)
        {
            (*values)[i] = pow (10., (logstart + i * diff));
        }
    }
    else
    {
        diff =
            ((*values)[intervals - 1] - (*values)[0]) / (MYREAL) (intervals - 1);
        for (i = 1; i < intervals - 1; i++)
        {
            (*values)[i] = (*values)[i - 1] + diff;
        }
    }
}

///
/// save all the options into a file (default) called parmfile;
/// uses save_options_buffer()
long save_parmfile (option_fmt * options, world_fmt * world, data_fmt *data)
{
    FILE *fp; 
    long bufsize;
    long allocbufsize = LONGLINESIZE;
    char *sbuffer;
    char *buffer;
    sbuffer = (char *) mycalloc (allocbufsize, sizeof (char));
    fp = options->parmfile;
    if (fp)
    {
        fclose (fp);
        openfile (&fp, options->parmfilename, "w", NULL);
    }
    else
    {
        openfile (&fp, options->parmfilename, "w", NULL);
    }
    options->parmfile = fp;
    bufsize = save_options_buffer (&sbuffer, &allocbufsize, options, data);
    buffer = sbuffer;
    fprintf (fp, "%s", buffer);
    fflush (fp);
    printf ("\n\n+++ Options were written to file %s in current directory +++\n\n", options->parmfilename);
    myfree(sbuffer);
    return bufsize;
}



/// prints delimiters between sections in the parmfile
void print_parm_delimiter(long *bufsize, char **buffer, long *allocbufsize)
{
    add_to_buffer("################################################################################\n",bufsize,buffer, allocbufsize);
}

/// prints delimiters between sections in the parmfile
void print_parm_smalldelimiter(long *bufsize, char **buffer, long *allocbufsize)
{
	add_to_buffer("#-------------------------------------------------------------------------------\n", bufsize,buffer, allocbufsize);
}

/// prints delimiters between sections in the parmfile
void print_parm_br(long *bufsize, char **buffer, long *allocbufsize)
{
	add_to_buffer("#\n", bufsize,buffer, allocbufsize);
}
/// prints parmfile single comment line
void print_parm_comment(long *bufsize, char **buffer, long *allocbufsize, char message[])
{
    char fp[LINESIZE];
	sprintf(fp,"# %s\n",message);
	add_to_buffer(fp, bufsize,buffer, allocbufsize);
}

/// prints parmfile mutable comment line
void print_parm_mutable_comment(long *bufsize, char **buffer, long *allocbufsize, char string[], ...)
{
    char message[LINESIZE];
	char fp[LINESIZE];
	va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);	
	sprintf(fp,"# %s\n",message);
	add_to_buffer(fp, bufsize,buffer, allocbufsize);
}

/// prints parmfile mutable option line
void print_parm_mutable(long *bufsize, char **buffer, long *allocbufsize, char string[], ...)
{
    char message[LINESIZE];
	char fp[LINESIZE];
	va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);	
	sprintf(fp,"%s\n",message);
	add_to_buffer(fp, bufsize,buffer, allocbufsize);
}

/// prints parmfile fixed option line
void print_parm(long *bufsize, char **buffer, long *allocbufsize, char string[])
{
	char fp[LINESIZE];
	sprintf(fp,"%s\n",string);
	add_to_buffer(fp, bufsize,buffer, allocbufsize);
}
/// prints a title of section in the parmfile
void print_parm_title(long *bufsize, char **buffer, long *allocbufsize, char message[])
{
	print_parm_delimiter(bufsize, buffer, allocbufsize);
	print_parm_comment(bufsize, buffer, allocbufsize, message);
	print_parm_delimiter(bufsize, buffer, allocbufsize);
}



/// print the ttratio into the parmfile buffer
void print_parm_ttratio(long *bufsize, char **buffer,  long *allocbufsize, option_fmt *options)
{
    long i=1;
    char fp[LINESIZE];
    sprintf(fp,"ttratio=%f ", options->ttratio[0]);
    add_to_buffer(fp,bufsize,buffer, allocbufsize);
    while (options->ttratio[i] != 0.0)
    {
        sprintf (fp, "%f ", options->ttratio[i++]);
        add_to_buffer(fp,bufsize,buffer, allocbufsize);
    }
    sprintf (fp, "\n");
    add_to_buffer(fp,bufsize,buffer, allocbufsize);
}

/// print base frequency into parmfile buffer
void print_parm_freqfrom(long *bufsize, char **buffer,  long *allocbufsize, option_fmt * options)
{
    if (options->freqsfrom)
        print_parm(bufsize, buffer, allocbufsize, "freqs-from-data=YES");
    else
      print_parm_mutable(bufsize,buffer, allocbufsize, "freqs-from-data=NO:%f,%f, %f, %f\n", options->freqa,
                           options->freqc, options->freqg, options->freqt);
}

/// print whether one uses categories and in what file they are into parmfile buffer
void print_parm_categs(long *bufsize, char **buffer,  long *allocbufsize, option_fmt * options)
{
    if (options->categs>1)
        print_parm_mutable(bufsize, buffer,  allocbufsize, "categories=%li:%s",options->categs, options->catfilename);
    else
        print_parm(bufsize, buffer, allocbufsize,  "categories=1 #no categories file specified");
}

/// print whether one uses weights and in what file they are into parmfile buffer
void print_parm_weights(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options)
{
    if (options->weights)
      print_parm_mutable(bufsize, buffer, allocbufsize, "weights=YES:%s", options->weightfilename);
    else
        print_parm(bufsize, buffer,  allocbufsize, "weights=NO");
}

/// print whether one uses rate categories, autocorrelation and in what file they are into parmfile buffer
void print_parm_rates(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options)
{
    long i;
    char fp[LINESIZE];
    long lbufsize = 0L;
    char *lbuffer;
    long alloclbufsize;
    // rates
    lbuffer = (char *) mycalloc(LINESIZE,sizeof(char));
    alloclbufsize = LINESIZE;
    for (i = 0; i < options->rcategs; i++)
    {
        sprintf (fp, "%f ", options->rrate[i]);
        add_to_buffer(fp,&lbufsize, &lbuffer, &alloclbufsize);
    }    
    print_parm_mutable(bufsize, buffer, allocbufsize, "rates=%li: %s",options->rcategs, lbuffer);
    lbufsize = 0L;
    //reset the lbuffer
    lbuffer[0] = '\0';
    // probablities
    for (i = 0; i < options->rcategs; i++)
    {
        sprintf (fp, "%f ", options->probcat[i]);
        add_to_buffer(fp,&lbufsize,&lbuffer, &alloclbufsize);
    }
    print_parm_mutable(bufsize, buffer, allocbufsize, "prob-rates=%li: %s",options->rcategs, lbuffer);
    myfree(lbuffer);
    // autocorrelation
    if (!options->autocorr)
      print_parm(bufsize, buffer, allocbufsize, "autocorrelation=NO");
    else
      print_parm_mutable(bufsize, buffer, allocbufsize, "autocorrelation=YES:%f", 1. / options->lambda);
}


/// print data-option into parmfile buffer
void print_parm_datatype(long *bufsize, char **buffer, long * allocbufsize, option_fmt *options)
{
    switch (options->datatype)
    {
        case 'a':
	  print_parm(bufsize, buffer, allocbufsize, "datatype=AllelicData");
            print_parm_mutable(bufsize, buffer, allocbufsize,"include-unknown=%s", options->include_unknown ? "YES" : "NO"); 
            break;
        case 'b':
            print_parm(bufsize, buffer, allocbufsize, "datatype=BrownianMicrosatelliteData");
            print_parm_mutable(bufsize, buffer, allocbufsize, "include-unknown=%s", options->include_unknown ? "YES" : "NO"); 
            break;
        case 'm':
            print_parm(bufsize, buffer, allocbufsize, "datatype=MicrosatelliteData");
	    if(options->msat_option == SINGLESTEP)
	      print_parm_mutable(bufsize, buffer, allocbufsize, "micro-submodel=%li", options->msat_option);
	    else
	      print_parm_mutable(bufsize, buffer, allocbufsize, "micro-submodel=%li:{%f, %f}", options->msat_option, options->msat_tuning[0], options->msat_tuning[1]);
            print_parm_mutable(bufsize, buffer, allocbufsize, "micro-threshold=%li", options->micro_threshold);
            print_parm_mutable(bufsize, buffer, allocbufsize, "include-unknown=%s", options->include_unknown ? "YES" : "NO"); 
            break;
        case 's':
            print_parm(bufsize, buffer, allocbufsize, "datatype=SequenceData");
            print_parm_ttratio(bufsize, buffer,  allocbufsize, options);
            print_parm_freqfrom(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate=%f", options->seqerror);
            print_parm_categs(bufsize, buffer, allocbufsize, options);
            print_parm_rates(bufsize, buffer, allocbufsize, options);
            print_parm_weights(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "interleaved=%s", options->interleaved ? "YES" : "NO"); 
            print_parm_mutable(bufsize, buffer, allocbufsize, "fast-likelihood=%s", options->fastlike ? "YES" : "NO"); 
            break;
        case 'h':
	  print_parm(bufsize, buffer, allocbufsize, "datatype=HapmapSNPfrequencydata\n");
            print_parm_ttratio(bufsize, buffer, allocbufsize, options);
            print_parm_freqfrom(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate=%f", options->seqerror);
            print_parm_categs(bufsize, buffer, allocbufsize, options);
            print_parm_rates(bufsize, buffer, allocbufsize, options);
            print_parm_weights(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "interleaved=%s", options->interleaved ? "YES" : "NO"); 
            print_parm_mutable(bufsize, buffer, allocbufsize, "fast-likelihood=%s", options->fastlike ? "YES" : "NO"); 
            break;
        case 'n':
            print_parm(bufsize, buffer, allocbufsize,"datatype=NucleotidePolymorphismData");
            print_parm_ttratio(bufsize, buffer, allocbufsize, options);
            print_parm_freqfrom(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate=%f", options->seqerror);
            print_parm_categs(bufsize, buffer, allocbufsize, options);
            print_parm_rates(bufsize, buffer, allocbufsize, options);
            print_parm_weights(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "interleaved=%s", options->interleaved ? "YES" : "NO"); 
            print_parm_mutable(bufsize, buffer, allocbufsize, "fast-likelihood=%s", options->fastlike ? "YES" : "NO"); 
            break;
        case 'u':
            print_parm(bufsize, buffer, allocbufsize,"datatype=UnlinkedSNPData");
            print_parm_ttratio(bufsize, buffer, allocbufsize, options);
            print_parm_freqfrom(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate=%f", options->seqerror);
            print_parm_categs(bufsize, buffer, allocbufsize, options);
            print_parm_rates(bufsize, buffer, allocbufsize, options);
            print_parm_weights(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "interleaved=%s", options->interleaved ? "YES" : "NO"); 
            print_parm_mutable(bufsize, buffer, allocbufsize, "fast-likelihood=%s", options->fastlike ? "YES" : "NO"); 
            break;
        case 'f':
            print_parm(bufsize, buffer, allocbufsize,"datatype=F-Ancestral state method");
            print_parm_ttratio(bufsize, buffer, allocbufsize, options);
            print_parm_freqfrom(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "seqerror-rate=%f", options->seqerror);
            print_parm_categs(bufsize, buffer, allocbufsize, options);
            print_parm_rates(bufsize, buffer, allocbufsize, options);
            print_parm_weights(bufsize, buffer, allocbufsize, options);
            print_parm_mutable(bufsize, buffer, allocbufsize, "interleaved=%s", options->interleaved ? "YES" : "NO"); 
            print_parm_mutable(bufsize, buffer, allocbufsize, "fast-likelihood=%s", options->fastlike ? "YES" : "NO"); 
            break;
        case 'g':
            print_parm(bufsize, buffer, allocbufsize,"datatype=GenealogySummaryOlderRun");
            break;
        default:
            error ("the parmfile-writer contains an error");
    }
}
 
void print_parm_tipdate(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data)
{
  long pos;
  long locus;
  char *input;
  if(options->has_datefile)
  {
    print_parm_mutable(bufsize, buffer, allocbufsize, "tipdate-file=YES:%s", options->datefilename); 
    print_parm_mutable(bufsize, buffer, allocbufsize, "generation-per-year=%f", options->generation_year); 
    input = (char *) mycalloc(options->mutationrate_year_numalloc * 50, sizeof(char));
    pos = sprintf(input, "{%20.20f", options->mutationrate_year[0]);
    for(locus=1; locus < options->mutationrate_year_numalloc; locus++)
      {
	pos += sprintf(input + pos,", %20.20f", options->mutationrate_year[locus]);
      }
    	//xcode   pos += 
    sprintf(input + pos,"}");
    print_parm_mutable(bufsize, buffer, allocbufsize, "mutationrate-per-year=%s", input);
    myfree(input);
  }
}

///
/// print the parmfile entry for the inheritance scalar
void print_parm_inheritence(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data)
{
  long pos;
  long locus;
  char *input;
    input = (char *) mycalloc(options->inheritance_scalars_numalloc * 50, sizeof(char));
    pos = sprintf(input, "inheritance-scalars={%20.20f", options->inheritance_scalars[0]);
    for(locus=1; locus < options->inheritance_scalars_numalloc; locus++)
      {
	pos += sprintf(input + pos,", %20.20f", options->inheritance_scalars[locus]);
      }
    //xcode pos += 
    sprintf(input + pos,"}");
    print_parm_mutable(bufsize, buffer, allocbufsize, "%s", input);
    myfree(input);
}
///
/// print the parmfile entry for the population relabeling
void print_parm_newpops(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options, data_fmt *data)
{
  long pos;
  long i;
  char *input;
    input = (char *) mycalloc(options->newpops_numalloc * 50, sizeof(char));
    pos = sprintf(input, "population-relabel={%li", options->newpops[0]);
    for(i=1; i < options->newpops_numalloc; i++)
      {
	pos += sprintf(input + pos,", %li", options->newpops[i]);
      }
    	//xcode   pos += 
    sprintf(input + pos,"}");
    print_parm_mutable(bufsize, buffer, allocbufsize, "%s", input);
    myfree(input);
}

void print_parm_randomsubset(long * bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
    if(options->randomsubset>0)
    {
        print_parm_mutable(bufsize, buffer, allocbufsize, "random-subset=%li", options->randomsubset);
        return;
    }
}
   
void print_parm_usertree(long * bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
    if(options->usertree)
    {
        print_parm_mutable(bufsize, buffer, allocbufsize, "usertree=TREE:%s", options->utreefilename);
        return;
    }
    if (options->randomtree)
    {
        print_parm(bufsize, buffer, allocbufsize, "usertree=RANDOMTREE");
        return;
    }
    if(options->dist)
    {
        print_parm_mutable(bufsize, buffer, allocbufsize, "usertree=DISTANCE:%s", options->distfilename);
        return;
    }
    print_parm(bufsize, buffer, allocbufsize, "usertree=AUTOMATIC");    
}    

/// print the theta starting parameters
void print_parm_theta(long *bufsize, char ** buffer, long *allocbufsize, option_fmt * options)
{
    long i;
    char fp[LINESIZE];
    char *lbuffer;
    long lbufsize = 1;
    long alloclbufsize = LINESIZE;
    lbuffer = (char *) mycalloc(alloclbufsize,sizeof(char));
    
    switch (options->numthetag)
      {
      case 0:
        if (strchr ("snupf", options->datatype))
	  {
	    sprintf (fp, "theta=Own:0.01");
            add_to_buffer(fp, &lbufsize, &lbuffer, &alloclbufsize);
	  }
        else
	  { 
	    sprintf (fp, "theta=Own:1.0");
            add_to_buffer(fp, &lbufsize, &lbuffer,  &alloclbufsize);
	  }
	break;
      case 1:
            sprintf (fp, "theta=Own:%f",options->thetag[0]);
            add_to_buffer(fp, &lbufsize, &lbuffer, &alloclbufsize);
	    break;
      default:
	sprintf (fp, "theta=Own:{");
	add_to_buffer(fp,&lbufsize, &lbuffer, &alloclbufsize);
	for (i = 0; i < options->numthetag - 1; i++)
	  {
	    sprintf (fp, "%f ", options->thetag[i]);
	    add_to_buffer(fp,&lbufsize,&lbuffer, &alloclbufsize);
	  }
	sprintf (fp, "%f}", options->thetag[i]);
	add_to_buffer(fp,&lbufsize, &lbuffer, &alloclbufsize);
	break;
      }
    print_parm_mutable(bufsize, buffer, allocbufsize, "%s", lbuffer);
}

/// print M starting parameters to the parmfile buffer
void print_parm_m(long *bufsize, char **buffer, long *allocbufsize, option_fmt * options)
{
    char fp[LINESIZE];
    long i, j, z, num;
    switch (options->nummg)
    {
        case 0:
            sprintf (fp, "migration=Own:1\n");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            break;
        case 1:
            sprintf (fp, "migration=Own:%f\n", options->mg[0]);
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            break;
        default:
            sprintf (fp, "migration=Own:{ ");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
            z = 0;
            num = (long) (1. + sqrt (4. * (MYREAL) options->nummg + 1.) / 2.);
            for (i = 0; i < num; i++)
            {
                for (j = 0; j < num; j++)
                {
                    if (i == j)
                    {
                        sprintf (fp, "- ");
                        add_to_buffer(fp,bufsize,buffer, allocbufsize);
                    }
                    else
                    {
                        sprintf (fp, "%f ", options->mg[z++]);
                        add_to_buffer(fp,bufsize,buffer, allocbufsize);
                    }
                }
            }
                sprintf (fp, "}\n");
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
    }
}

void print_parm_heating(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
    long i;
    char fp[LINESIZE];
    sprintf (fp, "heating=%s", options->heating ? (options->adaptiveheat!=NOTADAPTIVE  ? (options->adaptiveheat==STANDARD ? "ADAPTIVE_standard" : "Bounded_adaptive") : "YES") : "NO\n");
    add_to_buffer(fp,bufsize,buffer, allocbufsize);

    if (options->heating)
    {
        sprintf (fp, ":%li:{%f,", options->heating_interval, options->heat[0]);
        add_to_buffer(fp,bufsize,buffer, allocbufsize);
        
        for (i = 1; i < options->heated_chains - 1; i++)
        {
            sprintf (fp, "%f,", options->heat[i]);
            add_to_buffer(fp,bufsize,buffer, allocbufsize);
        }
        sprintf (fp, "%f}\nheated-swap=%s\n", options->heat[i], (options->heatedswap_off ? "NO" : "YES") );
        add_to_buffer(fp,bufsize,buffer, allocbufsize);
    }
}

///
/// returns a TRUE when the option is set and sets a filename, returns FALSE when the
/// option is not used and then of course does not set the filename
boolean  set_filename(char *value, char comparison[], char ** filename)
{
    unsigned long len=1;
    char *newvalue;
    char *temp;
    temp = (char *) mycalloc(LINESIZE,sizeof(char));
    upper(value, &temp);
    len = strlen(comparison);
    if (strncmp (temp, comparison, len)==0)
    {
        newvalue = strchr(value,':');
        if(newvalue!=NULL)
            strncpy (*filename, newvalue + 1, 255);
        myfree(temp);
        return TRUE;
    }
    myfree(temp);
    return FALSE;
}

/// assumes that an options will be used and takes the filename after the ":"
void  set_filename_only(boolean check, char *value, char ** filename)
{
    char *newvalue;
    if(check)
    {
        newvalue = strchr(value,':');
        if(newvalue != NULL)
            strncpy (*filename, newvalue + 1, 255);
    }
}
    

void   set_parm_prior_values(int priortype, prior_fmt * prior, char * mytext)
{
  char tmp1[LINESIZE];
  char tmp2[LINESIZE];
  char tmp3[LINESIZE];
  char tmp4[LINESIZE];
  mytext[0]='\0';
  show_priormin(tmp1, prior, priortype);
  strcat(mytext,tmp1);
  switch(priortype)
    {
      //    case SLICE:
      //show_priormax(tmp3, prior, priortype);
      //strcat(mytext,tmp3);
      //break;
      //    case MULTPRIOR:  
      //	  show_priormax(tmp3, prior, priortype);
      //	  show_priordelta(tmp4, prior,priortype);
      //	  strcat(mytext,tmp3);
      //	  strcat(mytext,tmp4);
      //	  break; 
	case EXPPRIOR: 
	  show_priormean(tmp2,prior, priortype);
	  show_priormax(tmp3, prior,priortype);
	  strcat(mytext,tmp2);
	  strcat(mytext,tmp3);
	  break;
	case WEXPPRIOR:
	  show_priormean(tmp2, prior, priortype);
	  show_priormax(tmp3, prior, priortype);
	  show_priordelta(tmp4, prior, priortype);
	  strcat(mytext,tmp2);
	  strcat(mytext,tmp3);
	  strcat(mytext,tmp4);
	  break;
	  //	case GAMMAPRIOR:
	  //show_priormean(tmp2, prior, priortype);
	  //show_priormax(tmp3, prior, priortype);
	  //strcat(mytext,tmp2);	  
	  //strcat(mytext,tmp3);
	  //break;
	case UNIFORMPRIOR:
	default:	  
	  show_priormax(tmp3, prior, priortype);
	  show_priordelta(tmp4, prior, priortype);
	  strcat(mytext,tmp3);
	  strcat(mytext,tmp4);
	  break;
    }
}


/// \brief returns proposal type string
/// returns a string that shows what proposal type is set
char * show_proposaltype(boolean priorset)
{
  if(priorset)
    {
      return  "SLICE Sampler" ;
    }
  else
    {
      return "METROPOLIS-HASTINGS Sampler";
    }
}

/// \brief returns priortype sting
/// returns a string that shows what prior distribution is set
char * show_parmpriortype(int priorset)
{
  switch(priorset)
    {
      //    case MULTPRIOR: return  "MULTPRIOR" ; 
    case EXPPRIOR: return "EXPPRIOR" ; 
    case WEXPPRIOR: return "WEXPPRIOR";
      //case GAMMAPRIOR: return "GAMMAPRIOR";
    case UNIFORMPRIOR: return "UNIFORMPRIOR";
    default: return "UNIFORMPRIOR";
    }
}

///
/// print prior values to buffer for parmfile
void print_parm_proposal(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
  print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-proposals= THETA %s",
		     show_proposaltype(options->slice_sampling[THETAPRIOR]));

  print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-proposals= MIG %s",
		     show_proposaltype(options->slice_sampling[MIGPRIOR]));
  if(options->bayesmurates)
    {
      print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-proposals= RATE %s",
			 show_proposaltype(options->slice_sampling[RATEPRIOR]));
    }
}

///
/// print prior values to buffer for parmfile
void print_parm_prior(long *bufsize, char **buffer, long *allocbufsize, option_fmt *options)
{
  char mytext[LINESIZE]="";

  set_parm_prior_values(options->bayesprior[THETAPRIOR], options->bayespriortheta, mytext);
  print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-priors= THETA %s: %s",
		     show_parmpriortype(options->bayesprior[THETAPRIOR]), mytext);

  set_parm_prior_values(options->bayesprior[MIGPRIOR], options->bayespriorm, mytext);
  print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-priors= MIG %s: %s",
		     show_parmpriortype(options->bayesprior[MIGPRIOR]), mytext);
  if(options->bayesmurates)
    {
      set_parm_prior_values(options->bayesprior[RATEPRIOR], options->bayespriorrate, mytext);
      print_parm_mutable(bufsize, buffer, allocbufsize, "bayes-priors= RATE %s: %s",
			 show_parmpriortype(options->bayesprior[RATEPRIOR]), mytext);
    }
}

///copy mu_rates into world->options
void set_meanmu(worldoption_fmt * wopt, option_fmt * options, long loci)
{
  long i;
  long n = 0;
  MYREAL last;

  if(options->mutationrate_year_numalloc == 0)
    {
      last = 1.0;
    }
  else
    {
      last = 0.0;
    }

  wopt->meanmu = (MYREAL *) mycalloc(loci, sizeof(MYREAL));

  for(i=0;i< options->mutationrate_year_numalloc;i++)
    {
      n++;
      wopt->meanmu[i] = options->mutationrate_year[i];
      last += (wopt->meanmu[i] - last ) / n;
    }

  for(i = options->mutationrate_year_numalloc; i < loci;i++)
      {
	wopt->meanmu[i] = last;
      }
}

///save option->mu_rates into a buffer
long save_mu_rates_buffer (char **buffer, long *allocbufsize, option_fmt * options)
{
    long i;
    long bufsize = 0;
    print_parm_mutable(&bufsize, buffer, allocbufsize, "%li ",options->muloci);
    for(i=0; i < options->muloci; i++)
      {
	print_parm_mutable(&bufsize, buffer, allocbufsize, "%f ",options->mu_rates[i]);
      }
    //    fprintf(stdout,"muloci=%li, muratebuffer:>%s<\n",options->muloci, *buffer, allocbufsize);
    return bufsize;
}

///save option->inheritance_scalars into a buffer
long save_inheritance_scalars_buffer (char **buffer, long *allocbufsize, option_fmt * options)
{
    long i;
    long bufsize = 0;
    print_parm_mutable(&bufsize, buffer, allocbufsize, "%li ",options->inheritance_scalars_numalloc);
    for(i=0; i < options->inheritance_scalars_numalloc; i++)
      {
	print_parm_mutable(&bufsize, buffer, allocbufsize, "%f ",options->inheritance_scalars[i]);
      }
    return bufsize;
}
///save option->newpops into a buffer
long save_newpops_buffer (char **buffer, long *allocbufsize, option_fmt * options)
{
    long i;
    long bufsize = 0;
    print_parm_mutable(&bufsize, buffer, allocbufsize, "%li ",options->newpops_numalloc);
    for(i=0; i < options->newpops_numalloc; i++)
      {
	print_parm_mutable(&bufsize, buffer, allocbufsize, "%f ",options->newpops[i]);
      }
    return bufsize;
}

///save options into a buffer
long save_options_buffer (char **buffer, long *allocbufsize, option_fmt * options, data_fmt *data)
{
    long i;
    char nowstr[LINESIZE] = "----";
    char fp[LINESIZE];
    long bufsize = 0;
    get_time (nowstr, "%c");
    if(options->bayes_infer)
        options->schains=0;
	// header for parmfile
	print_parm_delimiter(&bufsize, buffer, allocbufsize);	
    print_parm_mutable_comment(&bufsize, buffer, allocbufsize,  "Parmfile for Migrate %s-%s [do not remove these first TWO lines]", MIGRATEVERSION,MIGRATESUBVERSION);
    print_parm_comment(&bufsize, buffer, allocbufsize, "generated automatically on");
    print_parm_mutable_comment(&bufsize, buffer, allocbufsize, "%s", nowstr);
	print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "please report problems to Peter Beerli");
	print_parm_comment(&bufsize, buffer, allocbufsize, " email: beerli@fsu.edu");
	print_parm_comment(&bufsize, buffer, allocbufsize," http://popgen.csit.fsu.edu/migrate.html");
	print_parm_delimiter(&bufsize, buffer, allocbufsize);
	// general options
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_title(&bufsize, buffer, allocbufsize, "General options");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_comment(&bufsize, buffer, allocbufsize,"Interactive or batch job usage");
	print_parm_comment(&bufsize, buffer, allocbufsize,"  Syntax: menu= < YES | NO > ");
	print_parm_comment(&bufsize, buffer, allocbufsize,"For batch runs it needs to be set to NO");
	print_parm_mutable(&bufsize, buffer, allocbufsize,"menu=%s", options->menu ? "YES" : "NO ");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_comment(&bufsize, buffer, allocbufsize,"Specification of length of names of indiviudals");	
	print_parm_comment(&bufsize, buffer, allocbufsize,"   Syntax: nmlength=<INTEGER between 0 .. 30>");	
	print_parm_mutable(&bufsize, buffer, allocbufsize,"nmlength=%li", options->nmlength);
	print_parm_br(&bufsize, buffer, allocbufsize);

	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_title(&bufsize, buffer, allocbufsize, "Data options");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_comment(&bufsize, buffer, allocbufsize,"Several different main datatypes are possible:");
	print_parm_comment(&bufsize, buffer, allocbufsize,"INFINITE ALLELE: usable for electrophoretic markers,");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                 other markers with unknown mutation model");
	print_parm_comment(&bufsize, buffer, allocbufsize,"STEPWISE MUTATION: usable for microsatellite data or");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                 other markers with stepwise change");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                 from one allele to another");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                 [singlestep versus multistep model, see micro-submodel option]");
	print_parm_comment(&bufsize, buffer, allocbufsize,"FINITE SITES MUTATION: standard DNA/RNA sequence mutation");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                 model, usable for DNA or RNA contiguous");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                 sequences or varialbe sites only (SNP)");
	print_parm_comment(&bufsize, buffer, allocbufsize,"GENEALOGY SUMMARY: reanalyzing an old migrate run");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
	print_parm_comment(&bufsize, buffer, allocbufsize,"INFINITE ALLELE");
	print_parm_comment(&bufsize, buffer, allocbufsize," Syntax: datatype=ALLELICDATA ");
	print_parm_comment(&bufsize, buffer, allocbufsize,"         include-unknown=<YES | NO> with YES unknown alleles");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               are included into analysis, NO is the default");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_comment(&bufsize, buffer, allocbufsize,"STEPWISE MUTATION");
	print_parm_comment(&bufsize, buffer, allocbufsize," Syntax: datatype=<MICROSATELLITEDATA | BROWNIANDATA");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               MICRO specifies the standard stepwise mutation");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               model, the BROWNIAN is an approximation to this");
	print_parm_comment(&bufsize, buffer, allocbufsize,"         micro-submodel=<1|2:{tune,pinc}>");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                1 means singlestep mutation model (this is the default and the standard");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                2 is the Multistep model (see Watkins 2007 TPB, section 4.2) it needs");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                  two parameters: tune specifies how close the model is to a singlestep model");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                  so tune=0 --> singlestep, tune=1 --> infinite allele model;");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                  the second parameter defines the probability that the repeat number");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                  is increasing, this value cannot be larger than 0.666, I suggest 0.5.");
	print_parm_comment(&bufsize, buffer, allocbufsize,"                  Example: micro-submodel=2:{0.5,0.5}");
	print_parm_comment(&bufsize, buffer, allocbufsize,"         micro-threshold=<INTEGER> Default is 10 [MICRO only, NEEDS TO BE EVEN!],");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               smaller values speed up analysis, but might also");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               crash, large values slow down analysis considerably.");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               Change this value only when you suspect that your");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               data has huge gaps in repeat length.");
	print_parm_comment(&bufsize, buffer, allocbufsize,"         include-unknown=<YES | NO> with YES unknown alleles");
	print_parm_comment(&bufsize, buffer, allocbufsize,"               are included into analysis, NO is the default");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_comment(&bufsize, buffer, allocbufsize,"FINITE SITES MUTATION");
	print_parm_comment(&bufsize, buffer, allocbufsize," Syntax: datatype=<SEQUENCEDATA | NUCLEOTIDE | UNLINKEDSNPS | ANCESTRAL");
	print_parm_comment(&bufsize, buffer, allocbufsize,"        SEQENCEDATA: typical linked stretches of DNA, for example mtDNA");
    print_parm_comment(&bufsize, buffer, allocbufsize,"        NUCLEOTIDE: linked DNA stretches, all invariable sites removed");
    print_parm_comment(&bufsize, buffer, allocbufsize,"        UNLINKEDSNPS: each variable site is a locus, DO NOT USE THIS YET");
    print_parm_comment(&bufsize, buffer, allocbufsize,"        ANCESTRAL: instead taking into account all posible states, use");
    print_parm_comment(&bufsize, buffer, allocbufsize,"               use only the most likely state probability, DON'T USE THIS YET");
	print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize,"         freqs-from-data=<YES | NO: freq(A), freq(C), freq(G), freq(T)>");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               calculate the prior base frequencies from the data,");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               or specify the frequencies");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         ttratio=<RATIO1 RATIO2 ....> Default is 2.0,");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               ratio between transitions and transversions.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         seq-error=<VALUE> Default is 0.0, typical values for ABI 3700 ");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               sequencers after base calling are around 0.001 (1/650)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         categories=<VALUE:CATFILE> The categories are integers or letters");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               specified in file called CATFILE, this assumes that all");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               sites belong to known categories, this can be used to");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               weight third positions etc.");
	print_parm_comment(&bufsize, buffer, allocbufsize, "         rates=<VALUE1 VALUE2 ...> the rates are specified arbitrarily or");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               then are from a Gamma distribution with alpha=x, currently");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                the alpha value gets lost and is not recorded in the parmfile");
	print_parm_comment(&bufsize, buffer, allocbufsize, "         prob-rates=<RATE2 RATE1 ... > These rates can be arbitrary or ");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               generated with gamma-deviated rates and then are derived");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               using Laguerre's quadrature, this should get better");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               results than equal probability methods.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         autocorrelation=<NO | YES:VALUE> Default is NO");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               autocorrelation makes only sense with rates,");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               VALUE should be >1.0");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         weights=<NO | YES:WEIGHTFILE> The weights are specified");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               in file called WEIGHTFILE, this assumes that all sites");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               belong to known weights, this can be used to weight");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               portions of the sequence etc.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         interleaved=<YES | NO> Use either an interleaved or ");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               non-interleaved format. Default is NO,");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               interleaved=YES is discouraged");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         fast-likelihood=<YES | NO> Default is YES, use NO when you");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               have many hundred individuals and get strange errors");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               during a run, NO is scaling the conditional likelihood");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               so that very small values are >0.00000");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         inheritance-scalars={values for each locus}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               these values are multiplied with Theta, for example having");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               two autosomal and a locus on X- and one on Y-chromosome we would give ");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               inheritance-scalars={1 1 0.666 0.25}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               [if all loci have the same scalar, just use {1}, even for many loci]]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         population-relabel={assignment for each location in the infile}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               example is population-relabel={1 2 2}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         random-subset=number");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               allows to subset the dataset randomly, if number > sample in population");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               all samples are taken, if number is smaller then the pop sample is shuffled and");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               and the first number samples are taken");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         usertree=<NO | UPGMA | AUTOMATIC | TREE:TREEFILE | DISTANCE:DISTFILE | RANDOM>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               Default is RANDOM, NO delivers a UPGMA tree using the data");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               with TREE and DISTANCE the user needs to ");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               give a usertreefile or a pairwise distance file, with RANDOM");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               a random tree will be the starting tree");
    print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_br(&bufsize, buffer, allocbufsize);
    
    
	print_parm_datatype(&bufsize, buffer, allocbufsize, options);
    print_parm_tipdate(&bufsize, buffer, allocbufsize, options, data);
    print_parm_inheritence(&bufsize, buffer, allocbufsize, options, data);
    print_parm_newpops(&bufsize, buffer, allocbufsize, options, data);
    print_parm_randomsubset(&bufsize, buffer, allocbufsize, options);
    print_parm_usertree(&bufsize, buffer, allocbufsize, options);

#ifdef UEP
    // unique event polymorphisms
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_title(&bufsize, buffer, allocbufsize,  "Unique event polymorphism options");
	print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax uep=<NO | YES:UEPFILE >");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               Default is NO, with YES the user needs to ");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               give a file for each individual (same order as in datafile");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               the value: indiviudal<10 characters> uep-state<0|1>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         uep-rates= mu[0->1] : nu[0->1] mutation rate between the two UEP states");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         uep-bases= FREQ1 : FREQ2   prior frequency of the two alleles");
	print_parm_br(&bufsize, buffer, allocbufsize);

    if (options->uep)
    {
        print_parm_mutable(&bufsize, buffer, allocbufsize, "uep=YES:%s", options->uepfilename);
        print_parm_mutable(&bufsize, buffer, allocbufsize, "uep-rates=%f:%f", options->uepmu, options->uepnu);
        print_parm_mutable(&bufsize, buffer, allocbufsize, "uep-bases=%f:%f", options->uepfreq0, options->uepfreq1);
    }
    else
    {
        print_parm(&bufsize, buffer, allocbufsize, "uep=NO");
    }
#endif
    //input and output options
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_title(&bufsize, buffer, allocbufsize, "Input options");
	print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_comment(&bufsize, buffer, allocbufsize, "input file location");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax infile=FILEPATH");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "infile=%s", options->infilename);
	print_parm_br(&bufsize, buffer, allocbufsize);

    print_parm_comment(&bufsize, buffer, allocbufsize, "Random number seed specification");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax random-seed=<AUTO | OWN:< seedfile | value >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     AUTO           uses computer system clock to generate seed");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     OWN:seedfile   uses file seedfile with random number seed");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     OWN:value      uses number value for seed");
    switch (options->autoseed)
    {
    case NOAUTO:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "random-seed=OWN:%s",options->seedfilename);
        break;
    case AUTO:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "random-seed=AUTO #OWN:%li", options->inseed);
        break;
    case NOAUTOSELF:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "random-seed=OWN:%li", options->inseed);
        break;
    default:
        error ("error in wrinting parmfile start seed method unknown");
    }
	print_parm_br(&bufsize, buffer, allocbufsize);

    print_parm_comment(&bufsize, buffer, allocbufsize, "Specify the title of the run, will be overridden by title in datafile");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: title=title text [up to 80 characters]");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "title=%s", options->title);
	print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_title(&bufsize, buffer, allocbufsize, "Output options");
	print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Progress report to the window where the program was started");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: progress=<NO | YES | VERBOSE>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         NO       nothing is printed to the console");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         YES      some messages about progress are reported [default]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         VERBOSE  more messages are reported to console");    
    print_parm_mutable(&bufsize, buffer, allocbufsize, "progress=%s",
             options->progress ? (options->verbose ? "VERBOSE" : "YES") : "NO ");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
	print_parm_br(&bufsize, buffer, allocbufsize);

    print_parm_comment(&bufsize, buffer, allocbufsize, "Recording messages to screen into logfile");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax logfile=<NO | YES:logfilename>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      NONE     no recording of progress");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      logfilename  path to logfile");
    if(options->writelog)
        print_parm_mutable(&bufsize, buffer, allocbufsize, "logfile=YES:%s",  options->logfilename);
    else
        print_parm(&bufsize, buffer, allocbufsize, "logfile=NO");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
	print_parm_br(&bufsize, buffer, allocbufsize);

    print_parm_comment(&bufsize, buffer, allocbufsize, "Print the data as read into the program");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax print-data=<NO | YES>");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "print-data=%s", options->printdata ? "YES" : "NO");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
	print_parm_br(&bufsize, buffer, allocbufsize);

    print_parm_comment(&bufsize, buffer, allocbufsize, "Print output to file [default is outfile]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax outfile=outfilename");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "outfile=%s", options->outfilename);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);

#ifdef PRETTY
    print_parm_comment(&bufsize, buffer, allocbufsize, "Print output to a PDF file [default is outfile.pdf]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax pdf-outfile=outfilename.pdf");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "pdf-outfile=%s", options->pdfoutfilename);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);
#endif

    print_parm_comment(&bufsize, buffer, allocbufsize, "Report M (=migration rate/mutation rate) instead of 4Nm or 2 Nm or Nm");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax use-M=<NO | YES> Default is YES, the name 4Nm is ambiguous");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     for non-diploid data");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "use-M=%s", options->usem ? "YES" : "NO");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
	print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_comment(&bufsize, buffer, allocbufsize, "Plotting parameters: migration versus population size, such that Theta1 x immigration_.1");
    print_parm_comment(&bufsize, buffer, allocbufsize, "this shows the sum of all imigrations int a population");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax plot=<NO | YES:<BOTH | OUTFILE>:<LOG | STD>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "         {x-start, x-end, y-start, y-end}:<N | M>:interval>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     NO   do not show a plot");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     YES  show plot with following specifications");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          BOTH    print raw coordinates into MATHFILE and plot to OUTFILE");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          OUTFILE plot only to OUTFILE");
    print_parm_comment(&bufsize, buffer, allocbufsize, "            LOG   scaling of both axes");
    print_parm_comment(&bufsize, buffer, allocbufsize, "            STD   non-log scaling");
    print_parm_comment(&bufsize, buffer, allocbufsize, "            {...} plot range of both parameters");
    print_parm_comment(&bufsize, buffer, allocbufsize, "            N     use xNm to plot immigration, x=<1,2,3,4>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                  depending on the inheritance characteristic of the data");
    print_parm_comment(&bufsize, buffer, allocbufsize, "            M     plot migration rate/mutation rate as immigration axis");
    print_parm_comment(&bufsize, buffer, allocbufsize, "            interval the plot range is broken up into interval intervals");
    if (options->plot)
    {
        print_parm_mutable(&bufsize, buffer, allocbufsize, "plot=YES:%s:%s:{%f,%f,%f,%f}:%1.1s%li",
            options->plotmethod == PLOTALL ? "BOTH" : "OUTFILE",
            options->plotscale == PLOTSCALELOG ? "LOG" : "STD",
            options->plotrange[0],options->plotrange[1], options->plotrange[2], options->plotrange[3],
            options->plotvar == PLOT4NM ? "N" : "M", options->plotintervals);
    }
    else
    {
        print_parm(&bufsize, buffer, allocbufsize, "plot=NO");
    }
	print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Print plot data into a file ");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax: mathtfile=mathfile the values are printed in a mathematica readable way");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "mathfile=%s", options->mathfilename);
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
	print_parm_br(&bufsize, buffer, allocbufsize);

    print_parm_comment(&bufsize, buffer, allocbufsize, "Profile likelihood for each estimated parameter");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax profile=<NONE | <ALL | TABLES | SUMMARY>:");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              <PRECISE | DISCRETE | QUICK | FAST>  >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     NONE    do not calculate profile likelihoods");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     ALL     print individual profile tables and summary [default]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     TABLES  show only tables and no summary");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     SUMMARY show only summary");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          PRECISE  evaluate profile likelihood at percentiles [Default]");
    //    print_parm_comment(&bufsize, buffer, allocbufsize, "          SPLINE   use splines to calculate precentile values [DOES NOT WORK!]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          QUICK    assumes that there is no interaction of parameters");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          FAST     same as QUICK except in last calculation cycle assumes interaction");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          DISCRETE uses fixed mutipliers: 0.02,0.1,0.2,0.5,1,2,5,10,50");
    switch (options->profilemethod)
    {
        case 's':
            sprintf (fp, ":SPLINE");
            break;
        case 'd':
            sprintf (fp, ":DISCRETE");
            break;
        case 'q':
            sprintf (fp, ":QUICKANDDIRTY");
            break;
        case 'f':
            sprintf (fp, ":FAST");
            break;
        case 'p':
        default:
            sprintf (fp, ":PRECISE");
            break;
    }
    
    switch (options->profile)
    {
    case ALL:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "profile=ALL%s",fp);
        break;
    case TABLES:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "profile=TABLES%s",fp);
        break;
    case SUMMARY:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "profile=SUMMARY%s",fp);
        break;
    case _NONE:
    default:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "profile=NONE");
        break;
    }
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
	print_parm_br(&bufsize, buffer, allocbufsize);

    
    print_parm_comment(&bufsize, buffer, allocbufsize, "Print tree into treefile");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax print-tree=< NONE | <ALL | BEST | LASTCHAIN:Increment>:treefile >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        NONE no tree printed [Default, and only choice using parallel");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        ALL  print all visited genealogies [careful this will be huge]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        BEST print only the best tree visited");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        LASTCHAIN print all trees in last chain");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        with increment INCREMENT");
    switch (options->treeprint)
    {
    case _NONE:
        print_parm(&bufsize, buffer, allocbufsize, "print-tree=NONE");
        break;
    case ALL:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "print-tree=ALL:%s",options->treefilename);
        break;
    case BEST:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "print-tree=BEST:%s",options->treefilename);
        break;
    case LASTCHAIN:
        print_parm_mutable(&bufsize, buffer, allocbufsize, "print-tree=LASTCHAIN:%li:%s", options->treeinc, options->treefilename);
        break;
    default:
        print_parm(&bufsize, buffer, allocbufsize, "print-tree=NONE");
        break;
    }
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);
    
    
    print_parm_comment(&bufsize, buffer, allocbufsize, "write intermediate minimal statistics into a file for later use");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax write-summary=<NO | YES:SUMFILE >");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               Default is NO, with YES the user needs to ");
	print_parm_comment(&bufsize, buffer, allocbufsize, "               give a file to record the summary statistics");
     if(options->writesum)
         print_parm_mutable(&bufsize, buffer, allocbufsize, "write-summary=YES:%s",options->sumfilename);
     else
         print_parm(&bufsize, buffer, allocbufsize, "write-summary=NO");
    
     print_parm_br(&bufsize, buffer, allocbufsize);
     print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
     print_parm_br(&bufsize, buffer, allocbufsize);
          
     print_parm_comment(&bufsize, buffer, allocbufsize, "Likelihood ratio test");
     print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax l-ratio=<NO | YES:values_to_test>");
     print_parm_comment(&bufsize, buffer, allocbufsize, "      Values_to_test are compared to the values generated in the run");
     print_parm_comment(&bufsize, buffer, allocbufsize, "  values_to_test={ab..bbab..ba ... a}");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       the {} is a square matrix with values for the population sizes");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       on the diagonal and migration rates off-diagonal");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       the values a for the diagonal can be any of these:");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       number  constant, the value is for example 0.002");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       *       free to vary, the default is * for every parameter");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       m       mean of theta, this can be a subgroup of all thetas");
     print_parm_comment(&bufsize, buffer, allocbufsize, "               for example the theta 1-3 are averaged and thetas 4,5 are estimated");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       the values b for the migration rates can be any of these:");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       number  constant, the value is for example 45.0 or 0.0");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       *       free to vary, the default is * for every parameter");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       m       mean of M_ij, this can be a subgroup of migration rates");
     print_parm_comment(&bufsize, buffer, allocbufsize, "               for example the M_1-3i are averaged and M_4,5i are estimated");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       M       means of 4Nm (diploid), 2Nm (haploid), Nm (mtDNA, Y-chromosome)");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       s       symmetric migration rates M");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       S       symmetric migrants 4Nm");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       an example for 5 populations could look like this:");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       l-ratio=YES:{*s00s s*s00 0s*s0 00s*s s00s*");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       this describes a circular stepping stone model with 5 symmetric rates");
     print_parm_comment(&bufsize, buffer, allocbufsize, "        and independent sizes, a very basic stepping stone with 2 parameters would");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       look like this l-ratio=YES:{mm00m mmm00 0mmm0 00mmm m00mm}");
     print_parm_comment(&bufsize, buffer, allocbufsize, "       [The L-RATIO statement can be repeated]");
     print_parm_comment(&bufsize, buffer, allocbufsize, " Default: l-ratio=NO");
    for (i = 0; i < options->lratio->counter; i++)
    {
        if (options->lratio->data[i].type == MLE)
            print_parm_mutable(&bufsize, buffer, allocbufsize, "l-ratio=%s:%s", "YES",
                     options->lratio->data[i].value1);
        else
            print_parm_mutable(&bufsize, buffer, allocbufsize, "l-ratio=%s","NO");
            
        /* this code is unused, but this needs careful checking else where
        else
            print_parm_mutable(&bufsize, buffer, allocbufsize, "l-ratio=%s:%s:%s", "ARBITRARY",
                     options->lratio->data[i].value1,
                     options->lratio->data[i].value2);
        */
    }
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);
        
    print_parm_comment(&bufsize, buffer, allocbufsize, "AIC model selection [do not use yet, will come in Summer 2004]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax aic-modeltest=<NO | YES:<FAST | EXHAUSTIVE>>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      FAST        [do not use yet]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      EXHAUSTIVE  [do not use yet]");
    if (options->aic)
    {
        if (options->fast_aic)
            print_parm(&bufsize, buffer, allocbufsize, "aic-modeltest=YES:FAST");
        else
            print_parm(&bufsize, buffer, allocbufsize, "aic-modeltest=YES");
    }
    else
        print_parm(&bufsize, buffer, allocbufsize, "aic-modeltest=NO");
    print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);	
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Print a histogram of the time of migration events for each M(i->j)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax  mig-histogram=<NO | <ALL | MIGRATIONEVENTSONLY>:binsize:mighistfile >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        NO            do not record any events");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        ALL           record migration and coalescence event");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        MIGRATIONEVENTSONLY record only migration events");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        binsize has to be in mutation units, with an average Theta=0.01 try 0.001");
    print_parm_comment(&bufsize, buffer, allocbufsize, "Print a histogram of the parameters through time (skyline plot)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax  skyline=<NO | YES>:binsize:skylinefile >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        NO            do not calculate parameter estimates through time");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        YES           calculate parameters through time");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        binsize has to be in mutation units, with an average Theta=0.01 try 0.001");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        If the interval is too fine the output will be very noisy");

    if(options->mighist)
      {
	if(options->mighist_all)
	  print_parm_mutable(&bufsize, buffer, allocbufsize, "mig-histogram=ALL:%f:%s", options->eventbinsize,options->mighistfilename);
	else
	  print_parm_mutable(&bufsize, buffer, allocbufsize, "mig-histogram=YES:%f:%s", options->eventbinsize, options->mighistfilename);
	if(options->skyline)
	  print_parm_mutable(&bufsize, buffer, allocbufsize, "skyline=YES:%f:%s  #needs mig-histogram", options->eventbinsize, options->skylinefilename);
      }
    else
      {
        print_parm(&bufsize, buffer, allocbufsize, "mig-histogram=NO");
        print_parm(&bufsize, buffer, allocbufsize, "skyline=NO #needs mig-histogram=ALL:...");
      }
    print_parm_br(&bufsize, buffer, allocbufsize);


    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_title(&bufsize, buffer, allocbufsize, "Parameter start settings");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax: theta=<FST | OWN:<{value} | {value1, value2, ...., valuen} | NRANDOM:{mean std} | URANDOM{min,max}>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     migrationt=<FST | OWN:<{value} | {value1, value2, ...., valuen} | NRANDOM:{mean std} | URANDOM{min,max}>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       FST     starting parameter are derived from");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               an FST-like calculation (Beerli&Felsenstein 1999"); 
    print_parm_comment(&bufsize, buffer, allocbufsize, "       OWN     starting values are supplied by user");          
    print_parm_comment(&bufsize, buffer, allocbufsize, "          {value}   if only one value is supplied then all population");     
    print_parm_comment(&bufsize, buffer, allocbufsize, "                    have the same starting value");   
    print_parm_comment(&bufsize, buffer, allocbufsize, "          {value1, value2, ..., valuen} each population has its");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                    own starting value, if the number of values is");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                    insuffient, then the last value is the template");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                    for the remaining populations");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       NRANDOM  starting parameter is drawn randomely from a Normal distribution");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          {mean std} with mean and standard deviation");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       URANDOM  starting parameter is drawn randomely from a Uniform distribution");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          {min max} with minimum and maximum values");
    switch(options->thetaguess)
    {
    case FST:
      print_parm(&bufsize, buffer, allocbufsize, "theta=FST");
      break;
    case NRANDOMESTIMATE:
      print_parm_mutable(&bufsize, buffer, allocbufsize, "theta=NRANDOM: {%f %f}",options->thetag[0],options->thetag[1]);
      break;
    case URANDOMESTIMATE:
      print_parm_mutable(&bufsize, buffer, allocbufsize, "theta=URANDOM: {%f %f}",options->thetag[0],options->thetag[1]);
      break;

    default:
      print_parm_theta(&bufsize, buffer, allocbufsize, options);
	break;
    }

    switch (options->migrguess)
    {
    case FST:
      print_parm(&bufsize, buffer, allocbufsize, "migration=FST");
      break;
    case NRANDOMESTIMATE:
      print_parm_mutable(&bufsize, buffer, allocbufsize, "migration=NRANDOM: {%f %f}",options->mg[0],options->mg[1]);
      break;
    case URANDOMESTIMATE:
      print_parm_mutable(&bufsize, buffer, allocbufsize, "migration=URANDOM: {%f %f}",options->mg[0],options->mg[1]);
      break;
    default:
      print_parm_m(&bufsize, buffer, allocbufsize, options);
      break;
    }
	print_parm_br(&bufsize, buffer, allocbufsize);

	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Mutation rate modifiers");
	print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "  Syntax: mutation=<NOGAMMA | CONSTANT | ESTIMATE | GAMMA:alpha | OWN:loci: rate1 rate2 ... rate_loci>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     NOGAMMA      all loci have same mutation rate");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     CONSTANT     all loci have same mutation rate");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     ESTIMATE     BAYESIAN estimate: mutation rate is drawn from prior");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     GAMMA:alpha  ML estimate: mutation rate has Gamma distribution with alpha");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     OWN          mutation rate is different for every locus, but fixed");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        :loci: rate1, ...     number of loci, rate of locus 1, locus 2 etc.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     DATA         mutation rate modifier is deducted from loci in the data");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                  using Watterson's Theta and then scaling all rates Theta_locus/mean(Theta_loci");

    if (options->gamma)
    {
        print_parm_mutable(&bufsize, buffer, allocbufsize, "mutation=%s:%f", "GAMMA", options->alphavalue);
    }
    else
    {
        if (options->murates)
        {
	  if(options->murates_fromdata)
	    {
	      sprintf (fp, "mutation=DATA");
	      add_to_buffer(fp,&bufsize,buffer, allocbufsize);	      
	    }
	  else
	    {
	      sprintf (fp, "mutation=OWN:%li: ", options->muloci);
	      add_to_buffer(fp,&bufsize,buffer, allocbufsize);
	      for (i = 0; i < options->muloci; i++)
		{
		  sprintf (fp, "%f ", options->mu_rates[i]);
		  add_to_buffer(fp,&bufsize,buffer, allocbufsize);
		}
	      sprintf (fp, " \n");
	      add_to_buffer(fp,&bufsize,buffer, allocbufsize);
	    }
	}
        else
        {
	  if(options->bayesmurates)
	    {
	      print_parm(&bufsize, buffer, allocbufsize, "mutation=ESTIMATE");
	    }
	  else
	    {
	      print_parm(&bufsize, buffer, allocbufsize, "mutation=CONSTANT");
	    }
        }
    }
    print_parm_br(&bufsize, buffer, allocbufsize);

    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "FST model");
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_mutable(&bufsize, buffer, allocbufsize, "fst-type=%s", options->fsttype ? "THETA" : "MIGRATION");
	print_parm_br(&bufsize, buffer, allocbufsize);
    
	print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Custom migration model");
	print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: custom-migration={ab..bbab..ba ... a}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       the {} is a square matrix with values for the population sizes");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       on the diagonal and migration rates off-diagonal");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       the values a for the diagonal can be any of these:");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       c       constant, the value needs to be defined in the theta option");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       *       free to vary, the default is * for every parameter");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       m       mean of theta, this can be a subgroup of all thetas");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               for example the theta 1-3 are averaged and thetas 4,5 are estimated");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       the values b for the migration rates can be any of these:");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       c       constant, the value needs to be defined in the migration option");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       *       free to vary, the default is * for every parameter");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       m       mean of M_ij, this can be a subgroup of migration rates");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               for example the M_1-3i are averaged and M_4,5i are estimated");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       M       means of 4Nm (diploid), 2Nm (haploid), Nm (mtDNA, Y-chromosome)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       s       symmetric migration rates M");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       S       symmetric migrants 4Nm");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       an example for 5 populations could look like this:");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       custom-migration={*s00s s*s00 0s*s0 00s*s s00s*");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       this describes a circular stepping stone model with 5 symmetric rates");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        and independent sizes, a very basic stepping stone with 2 parameters would");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       look like this custom-migration={mm00m mmm00 0mmm0 00mmm m00mm}");
    printf("%li: %s\n", strlen(options->custm),options->custm);
    print_parm_mutable(&bufsize, buffer, allocbufsize, "custom-migration={%s}", options->custm);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Influence of geography on migration rate");
    print_parm_comment(&bufsize, buffer, allocbufsize, "a distance matrix between populations changes the migration rate matrix so that");
    print_parm_comment(&bufsize, buffer, allocbufsize, "(genetic?) migration rates =  inferred migration rate / distance ~ a dispersion coefficient");
    print_parm_comment(&bufsize, buffer, allocbufsize, "the geofile contains a number of populations, names for populations (10 characters), they");
    print_parm_comment(&bufsize, buffer, allocbufsize, "need to be in order of the dataset. And the distances between the populations, they do not");
    print_parm_comment(&bufsize, buffer, allocbufsize, "need to be symmetric");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: geo:<NO | YES:filename>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "            NO       distances among populations are considered to be 1 [all equal]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "            YES      distances are read from a file");

    if (options->geo)
    {
        print_parm_mutable(&bufsize, buffer, allocbufsize, "geo=YES:%s", options->geofilename);
    }
    else
    {
        print_parm(&bufsize, buffer, allocbufsize, "geo=NO");
    }
	print_parm_br(&bufsize, buffer, allocbufsize);
	print_parm_br(&bufsize, buffer, allocbufsize);
    // SEARCH STRATEGIES
	print_parm_title(&bufsize, buffer, allocbufsize, "Search strategies");
	print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_comment(&bufsize, buffer, allocbufsize, "MCMC Strategy method");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: bayes-update=< NO | YES>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       NO      maximum likelihood method");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       YES     Bayesian method");
    print_parm_comment(&bufsize, buffer, allocbufsize, "Some of the options are only available in one or other mode");
    print_parm_comment(&bufsize, buffer, allocbufsize, "BAYESIAN OPTIONS");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-updatefreq=VALUE ");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           VALUE      is a ratio between 0 and 1");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                      ratio of how many times the genealogy is updated compared to the parameters");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                      If the value is 0.4 in a 2 population scenario and with 1000000 steps");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                      The tree will be evaluated 400000 times, Theta_1, Theta_2, M_21, and M_12");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                       will be each evaluated 125000 times.");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-posteriorbins=VALUE VALUE");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           VALUE      is the number of bins in the psterior distribution histogram for Theta or M");
#ifdef PRETTY
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-posteriormaxtype=< ALL | P99 | MAXP99 | P100 >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           ALL        plots the WHOLE prior-parameter range\n");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           P99        plots from the minimum prior range value to\n");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                      the 99% percentile value of EACH parameter\n");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           MAXP99     sets all axes from minimum to the maximal\n");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                      99% percentile value of ALL parameter\n");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           P100       plots from the minimum prior range value to\n");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                      the 100% percentile value of EACH parameter\n");
#endif
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-file=<YES:FILENAME|NO>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           FILENAME is the name of the file that will contain");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                   the results for the posterior distribution");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-allfile=<YES:INTERVAL:FILENAME|NO>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           FILENAME is the name of the file that will contain");
    print_parm_comment(&bufsize, buffer, allocbufsize, "                   all parameters of the posterior distribution [HUGE]");
    //OLD    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-allfileinterval=INTERVAL");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           INTERVAL is the interval at which all parameters are written to file\n");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       ");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-proposals= THETA < SLICE | METROPOLIS >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-proposals= MIG < SLICE | METROPOLIS >");
    //    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-proposal= RATE < SLICE | METROPOLIS >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              SLICE uses the slice sampler to propose new parameter values");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              METROPOLIS uses the Metropolis-Hastings sampler");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              (this is done for each parameter group: THETA or MIGration)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       ");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-priors= THETA <UNIFORM unipriorvalues | EXP exppriorvalues | WINDOWEXP wexppriorvalues ");//| MULT multpriorvalues | FIXED>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-priors= MIG <UNIFORM unipriorvalues | EXP exppriorvalues | WINDOWEXP wexppriorvalues ");//| MULT multpriorvalues | FIXED>");
    //    print_parm_comment(&bufsize, buffer, allocbufsize, "       bayes-priors= RATE <UNIFORM unipriorvalues | EXP exppriorvalues | WINDOWEXP wexppriorvalues | MULT multpriorvalues | FIXED>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               unipriorvalues: min max delta");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               exppriorvalues: min mean max");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               wexppriorvalues: min mean max delta");
    //    print_parm_comment(&bufsize, buffer, allocbufsize, "               multpriorvalues: min max mult");
    //print_parm_comment(&bufsize, buffer, allocbufsize, "               fixed = no value is changed during the run");
	print_parm_br(&bufsize, buffer, allocbufsize);

    print_parm_comment(&bufsize, buffer, allocbufsize, "Maximum likelihood OPTIONS");    
    print_parm_comment(&bufsize, buffer, allocbufsize, "       short-chains=VALUE   VALUE is 1..n [Default is 10]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       short-inc=VALUE      VALUE is the number of updates that are not recorded");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       short-sample=VALUE   VALUE is the number of sampled updates");
	print_parm_br(&bufsize, buffer, allocbufsize);

    print_parm_comment(&bufsize, buffer, allocbufsize, "Search OPTIONS for both strategies");    
    print_parm_comment(&bufsize, buffer, allocbufsize, "       long-chains=VALUE   VALUE is 1..n [Default is 3]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       long-inc=VALUE      VALUE is the number of updates that are not recorded");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       long-sample=VALUE   VALUE is the number of sampled updates");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       burn-in=VALUE       VALUE is the number of updates to discard at the beginning");
	print_parm_br(&bufsize, buffer, allocbufsize);

    if(options->bayes_infer)
    {
        print_parm(&bufsize, buffer, allocbufsize, "bayes-update=YES");
        print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-updatefreq=%f", options->updateratio);
	if(options->bayesmurates)
	  print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-posteriorbins=%li %li %li", 
			     options->bayespriortheta->bins, 
			     options->bayespriorm->bins,
			     options->bayespriorrate->bins);
	else
	  print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-posteriorbins=%li %li", 
			     options->bayespriortheta->bins, 
			     options->bayespriorm->bins);
        print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-posteriormaxtype=%s", 
			   (options->bayespretty == PRETTY_P99 ? "P99" :
			   (options->bayespretty == PRETTY_MAX ? "ALL" :
			   (options->bayespretty == PRETTY_P100 ? "TOTAL" : "TOTAL"))));

	if(options->has_bayesfile)
	  print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-file=YES:%s", options->bayesfilename);
	else
	  print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-file=NO");

	if(options->has_bayesmdimfile)
	  {
	    print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-allfile=YES:%li:%s",  
			       options->bayesmdiminterval, options->bayesmdimfilename);
	    //OLD	    print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-allfileinterval=%li", options->bayesmdiminterval);
	  }
	else
	  print_parm_mutable(&bufsize, buffer, allocbufsize, "bayes-allfile=NO");
	
	print_parm_proposal(&bufsize, buffer, allocbufsize, options);
	print_parm_prior(&bufsize, buffer, allocbufsize, options);
    }                
    else
      {
	print_parm(&bufsize, buffer, allocbufsize, "bayes-update=NO");
	print_parm_mutable(&bufsize, buffer, allocbufsize, "short-chains=%li", options->schains);
	print_parm_mutable(&bufsize, buffer, allocbufsize, "short-inc=%li", options->sincrement);
	print_parm_mutable(&bufsize, buffer, allocbufsize, "short-sample=%li", options->ssteps);
      }
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_mutable(&bufsize, buffer, allocbufsize, "long-chains=%li", options->lchains);
    print_parm_mutable(&bufsize, buffer, allocbufsize, "long-inc=%li", options->lincrement);
    print_parm_mutable(&bufsize, buffer, allocbufsize, "long-sample=%li", options->lsteps);
    print_parm_mutable(&bufsize, buffer, allocbufsize, "burn-in=%li%c", options->burn_in, options->burnin_autostop);
    
    
    print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Schemes to improve MCMC searching and/or thermodynamic integration");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Heating schemes {MCMCMC = MC cubed}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: heating=< NO | <YES | ADAPTIVE>:SKIP:TEMPERATURES");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       NO    No heating");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       YES   heating using TEMPERATURES");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       ADAPTIVE adaptive heating using start TEMPERATURES [fails sometimes!!]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       SKIP skip that many comparisons, this lengthens the run by SKIP");
    print_parm_comment(&bufsize, buffer, allocbufsize, "           TEMPERATURES    { 1.0, temp1, temp2, temp3 .. tempn}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "    Example: heating=YES:1:{1.0, 1.2, 3.0,6.0}");
    print_parm_comment(&bufsize, buffer, allocbufsize, "Heating:  swapping chains");
    print_parm_comment(&bufsize, buffer, allocbufsize, "    Syntax: heated-swap=< YES | NO >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        YES  swapping of chains enabled [DEFAULT]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "        NO   swapping of chains disabled");
    print_parm_comment(&bufsize, buffer, allocbufsize, "     Example: heated-swap=YES");
    print_parm_heating(&bufsize, buffer, allocbufsize, options);
    
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "Lengthening chain schemes");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: moving-steps=< NO | YES:VALUE>");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      VALUE   frequency is between 0..1");

    if (options->movingsteps)
    {
        print_parm_mutable(&bufsize, buffer, allocbufsize, "moving-steps=YES:%f", options->acceptfreq);
    }
    else
    {
        print_parm(&bufsize, buffer, allocbufsize, "moving-steps=NO");
    }
	print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: long-chain-epsilon=VALUE");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      VALUE    is between 0..INFINITY");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               the VALUE is the likelihood ratio between the old and thew chain");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               the VALUE depends on the number of parameters: with 1 values of 0.5 are great");
    print_parm_comment(&bufsize, buffer, allocbufsize, "               but with many parameters values and bad data >20 is more reasonable");
    if (options->lcepsilon < LONGCHAINEPSILON)
       print_parm_mutable(&bufsize, buffer, allocbufsize, "long-chain-epsilon=%f", options->lcepsilon);
    else
        print_parm_mutable(&bufsize, buffer, allocbufsize, "long-chain-epsilon=INFINITY");

    
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Convergence statistic [Gelman and Rubin]");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: gelman-convergence=< YES:Pairs|Summary | NO >");    
    print_parm_comment(&bufsize, buffer, allocbufsize, "      NO      do not use Gelman's convergence criterium");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      YES     use Gelman's convergence criteria between chain i, and i-1");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              PAIRS reports all replicate pairs");
    print_parm_comment(&bufsize, buffer, allocbufsize, "              SUM   reports only mean and maxima");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "gelman-convergence=%s", options->gelman ? (options->gelmanpairs ? "Yes:Pairs" : "Yes:Sum" ) : "No");
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax: replicate=< NO | YES:<VALUE | LastChains> >");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      NO     no replication of run");
    print_parm_comment(&bufsize, buffer, allocbufsize, "      YES    replicate run");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          VALUE     number between 2 and many, complete replicates");
    print_parm_comment(&bufsize, buffer, allocbufsize, "          LastChains  replications over last chains");
    if(!options->replicate)
    {
        print_parm(&bufsize, buffer, allocbufsize, "replicate=NO");
    }
    else
    {
        if(options->replicatenum==0)
            print_parm(&bufsize, buffer, allocbufsize, "replicate=YES:LastChains");
        else
            print_parm_mutable(&bufsize, buffer, allocbufsize, "replicate=YES:%li", options->replicatenum);
    }
    print_parm_br(&bufsize, buffer, allocbufsize);

/// \todo remove this cpu section
//    if (options->cpu > 1)
//    {
//        sprintf (fp, "# number of cpus available [not really used right now]cpu=%i", options->cpu);
//		add_to_buffer(fp,&bufsize,buffer, allocbufsize);
//    }

    print_parm_comment(&bufsize, buffer, allocbufsize, "Migration rates are attracted to zero (fatal attraction)");
    print_parm_comment(&bufsize, buffer, allocbufsize, "Resistance is the lowest migration value for all but the last chain");
    print_parm_comment(&bufsize, buffer, allocbufsize, "   Syntax resistance=VALUE");
    print_parm_comment(&bufsize, buffer, allocbufsize, "       VALUE is the lowest migration rate value allowed during all but the last chain");
    print_parm_comment(&bufsize, buffer, allocbufsize, "             typical values are 0.01 or _lower_ for data with sequences and 0.0001 or _lower_ for other data");
    print_parm_mutable(&bufsize, buffer, allocbufsize, "resistance=%f", options->minmigsumstat);
    print_parm_br(&bufsize, buffer, allocbufsize);
    
    print_parm_smalldelimiter(&bufsize, buffer, allocbufsize);
    print_parm_br(&bufsize, buffer, allocbufsize);
    print_parm(&bufsize, buffer, allocbufsize, "end");
    return bufsize;
}

/*private functions============================================= */
///
/// fill the prior information of bayes theta
void set_bayes_options(char *value, option_fmt *options)
{
  int ptype;
  //long keys = 0;
  char paramtype[LINESIZE], priortype[LINESIZE];
  prior_fmt *prior=NULL;
  sscanf(value,"%s%s", paramtype, priortype);
  switch(uppercase(paramtype[0]))
    {
    case 'T'/* THETA*/: prior = options->bayespriortheta; ptype = THETAPRIOR; break;
    case 'M'/* MIG  */: prior = options->bayespriorm; ptype = MIGPRIOR; break;
    case 'R'/* RATE */: prior = options->bayespriorrate; ptype = RATEPRIOR; break;
    default: return;
    }
# ifdef USE_MYREAL_FLOAT 
    switch(uppercase(priortype[0]))
      {
      case 'S'/*slice sampler with uniform prior*/:
	keys = sscanf(value,"%s%s%f%f",  paramtype, priortype, &prior->min, &prior->max);
	options->bayesprior[ptype] = SLICE;
	options->slice_sampling[ptype] = TRUE;
	prior_consistency(options->bayespriortheta, ptype);
	break;
      case 'M'/*multprior   */:
	keys = sscanf(value,"%s%s%f%f%f",  paramtype, priortype, &prior->min, &prior->max, &prior->delta);
	options->bayesprior[ptype] = MULTPRIOR;
	prior_consistency(options->bayespriortheta, ptype);
	break;
      case 'E'/*expprior    */:   
	keys = sscanf(value,"%s%s%f%f%f",  paramtype, priortype, &prior->min, &prior->mean, &prior->max);
	options->bayesprior[ptype] = EXPPRIOR;
	prior_consistency(options->bayespriortheta, ptype);
	break;
      case 'W'/*wexpprior   */:   
	keys = sscanf(value,"%s%s%f%f%f%f",  paramtype, priortype, &prior->min, &prior->mean, &prior->max, &prior->delta);
	options->bayesprior[ptype] = WEXPPRIOR;
	prior_consistency(options->bayespriortheta, ptype);
	break;
      case 'G'/*gammaprior  */:   
	keys = sscanf(value,"%s%s%f%f%f",  paramtype, priortype, &prior->min, &prior->mean, &prior->max);
	options->bayesprior[ptype] = GAMMAPRIOR;
	prior_consistency(options->bayespriortheta, ptype);
	break;
      case 'U'/*uniformprior*/:   
	keys = sscanf(value,"%s%s%f%f%f",  paramtype, priortype, &prior->min, &prior->max,&prior->delta);
	options->bayesprior[ptype] = UNIFORMPRIOR;
	prior_consistency(options->bayespriortheta, ptype);
	prior->mean = (prior->max + prior->min) / 2.;
	//prior->delta = (prior->max - prior->min) / 10.;
	break;
      }
#else
    switch(uppercase(priortype[0]))
      {
      case 'S'/* slice sampler with uniform prior */:
	sscanf(value,"%s%s%lf%lf%lf",  paramtype, priortype, &prior->min, &prior->max,&prior->delta);
	options->bayesprior[ptype] = SLICE;
	options->slice_sampling[ptype] = TRUE;
	prior_consistency(options->bayespriortheta, ptype);
	prior->mean = (prior->max + prior->min) / 2.;
	//prior->delta = (prior->max - prior->min) / 10.;
	break;
      case 'M'/*multprior   */:
	sscanf(value,"%s%s%lf%lf%lf",  paramtype, priortype, &prior->min, &prior->max, &prior->delta);
	options->bayesprior[ptype] = MULTPRIOR;
	prior_consistency(options->bayespriortheta, ptype);
	break;
      case 'E'/*expprior    */:   ;
	sscanf(value,"%s%s%lf%lf%lf",  paramtype, priortype, &prior->min, &prior->mean, &prior->max);
	options->bayesprior[ptype] = EXPPRIOR;
	prior_consistency(options->bayespriortheta, ptype);
	break;
      case 'W'/*wexpprior   */:   ;
	sscanf(value,"%s%s%lf%lf%lf%lf",  paramtype, priortype, &prior->min, &prior->mean, &prior->max, &prior->delta);
	options->bayesprior[ptype] = WEXPPRIOR;
	prior_consistency(options->bayespriortheta, ptype);
	break;
      case 'G'/*gammaprior  */:   ;
	sscanf(value,"%s%s%lf%lf%lf",  paramtype, priortype, &prior->min, &prior->mean, &prior->max);
	options->bayesprior[ptype] = GAMMAPRIOR;
	prior_consistency(options->bayespriortheta, ptype);
	break;
      case 'U'/*uniformprior*/:   ;
	sscanf(value,"%s%s%lf%lf%lf",  paramtype, priortype, &prior->min, &prior->max,&prior->delta);
	options->bayesprior[ptype] = UNIFORMPRIOR;
	prior_consistency(options->bayespriortheta, ptype);
	prior->mean = (prior->max + prior->min) / 2.;
	//	prior->delta = (prior->max - prior->min) / 10.;
	break;
      }
#endif
}


///
/// checks prior consistency
void prior_consistency(prior_fmt *prior, int type)
{
  MYREAL tmp=0.0;
  MYREAL minim = 0.;
  MYREAL maxim = HUGE;
  switch(type)
    {
    case THETAPRIOR:
      minim = SMALLEST_THETA;
      maxim = BIGGEST_THETA;
      break;
    case MIGPRIOR:
      minim = SMALLEST_MIGRATION;
      maxim = BIGGEST_MIGRATION;
      break;
    case RATEPRIOR:
      minim = SMALLEST_RATE;
      maxim = BIGGEST_RATE;
      break;
    }
  if(prior->min < minim)
    prior->min = minim;
  if(prior->max > maxim)
    prior->max = maxim;
  if(prior->min > prior->max)
    {
      warning("maximum of prior is smaller (%f) than minimum (%f)! correct this problem!",prior->max, prior->min);
      tmp = prior->max;
      prior->max = prior->min;
      prior->min = tmp;
    }
}

///
/// fills the prior information into the option structure old version for compatibility
void set_bayes_options_oldstyle(char *value, option_fmt *options, int  priortype)
{
    long keys = 0;
    char *tmp;
    char *keeptmp;
    char priort[LINESIZE], part[LINESIZE], parm[LINESIZE];
    
    error("please change the bayes prior description, there is a more modern way to do that!");
    
    tmp = (char *) mycalloc (LINESIZE, sizeof (char));
    keeptmp = tmp;
    
    options->bayesprior[THETAPRIOR] = priortype;
    options->bayesprior[MIGPRIOR] = priortype;
    options->bayesprior[RATEPRIOR] = priortype;
    switch(priortype)
    {
    case SLICE:
# ifdef USE_MYREAL_FLOAT 
            keys = sscanf(value,"%s%s%f%f%s%f%f", priort, part, &options->bayespriortheta->min, &options->bayespriortheta->max,
                  parm,  &options->bayespriorm->min, &options->bayespriorm->max);
#else
            keys = sscanf(value,"%s%s%lf%lf%s%lf%lf", priort, part, &options->bayespriortheta->min, &options->bayespriortheta->max,
                  parm,  &options->bayespriorm->min, &options->bayespriorm->max);
#endif
            if(keys!=7)
                  error("abort due to problem with bayes settings in parmfile: error in SLICE sampler specification\n");
	    prior_consistency(options->bayespriortheta, THETAPRIOR);
	    prior_consistency(options->bayespriorm, MIGPRIOR);
            options->bayespriortheta->mean = (options->bayespriortheta->max + options->bayespriortheta->min)/2.;
            options->bayespriortheta->delta = (options->bayespriortheta->max - options->bayespriortheta->min)/10.;//about a tenth of the possible range
            options->bayespriorm->mean = (options->bayespriorm->max + options->bayespriorm->min)/2.;
            options->bayespriorm->delta = (options->bayespriorm->max - options->bayespriorm->min)/10.;
	    break;
        case UNIFORMPRIOR:
# ifdef USE_MYREAL_FLOAT 
            keys = sscanf(value,"%s%s%f%f%s%f%f", priort, part, &options->bayespriortheta->min, &options->bayespriortheta->max,
                  parm,  &options->bayespriorm->min, &options->bayespriorm->max);
#else
            keys = sscanf(value,"%s%s%lf%lf%s%lf%lf", priort, part, &options->bayespriortheta->min, &options->bayespriortheta->max,
                  parm,  &options->bayespriorm->min, &options->bayespriorm->max);
#endif
            if(keys!=7)
                  error("abort due to problem with bayes settings in parmfile: error in UNIFORM prior specification\n");
	    prior_consistency(options->bayespriortheta, THETAPRIOR);
	    prior_consistency(options->bayespriorm, MIGPRIOR);

            options->bayespriortheta->mean = (options->bayespriortheta->max + options->bayespriortheta->min)/2.;
            options->bayespriortheta->delta = (options->bayespriortheta->max - options->bayespriortheta->min)/10.;//about a tenth of the possible range
            options->bayespriorm->mean = (options->bayespriorm->max + options->bayespriorm->min)/2.;
            options->bayespriorm->delta = (options->bayespriorm->max - options->bayespriorm->min)/10.;
            break;
        case EXPPRIOR:
# ifdef USE_MYREAL_FLOAT 
            keys = sscanf(value,"%s%s%f%f%f%s%f%f%f", priort, part, 
                          &options->bayespriortheta->min, 
                          &options->bayespriortheta->mean, 
                          &options->bayespriortheta->max,
                          parm,  &options->bayespriorm->min,                           
                          &options->bayespriorm->mean,
                          &options->bayespriorm->max);
#else
            keys = sscanf(value,"%s%s%lf%lf%lf%s%lf%lf%lf", priort, part, 
                          &options->bayespriortheta->min, 
                          &options->bayespriortheta->mean, 
                          &options->bayespriortheta->max,
                          parm,  &options->bayespriorm->min,                           
                          &options->bayespriorm->mean,
                          &options->bayespriorm->max);
#endif
            if(keys!=9)
                error("abort due to problem with bayes settings in parmfile: error in EXP prior specification\n");
	    prior_consistency(options->bayespriortheta, THETAPRIOR);
	    prior_consistency(options->bayespriorm, MIGPRIOR);

            options->bayespriortheta->delta = (options->bayespriortheta->max - options->bayespriortheta->min)/10.;//about a tenth of the possible range
            options->bayespriorm->delta = (options->bayespriorm->max - options->bayespriorm->min)/10.;
            break;
        case WEXPPRIOR:
# ifdef USE_MYREAL_FLOAT 
            keys = sscanf(value,"%s%s%f%f%f%f%s%f%f%f%f", priort, part, 
                          &options->bayespriortheta->min, 
                          &options->bayespriortheta->mean, 
                          &options->bayespriortheta->max,
                          &options->bayespriortheta->delta,
                          parm,  &options->bayespriorm->min,                           
                          &options->bayespriorm->mean,
                          &options->bayespriorm->max,
                            &options->bayespriorm->delta);
#else
            keys = sscanf(value,"%s%s%lf%lf%lf%lf%s%lf%lf%lf%lf", priort, part, 
                          &options->bayespriortheta->min, 
                          &options->bayespriortheta->mean, 
                          &options->bayespriortheta->max,
                          &options->bayespriortheta->delta,
                          parm,  &options->bayespriorm->min,                           
                          &options->bayespriorm->mean,
                          &options->bayespriorm->max,
                        &options->bayespriorm->delta);
#endif
            if(keys!=11)
                error("abort due to problem with bayes settings in parmfile: error in WEXP prior specification\n");
	    prior_consistency(options->bayespriortheta, THETAPRIOR);
	    prior_consistency(options->bayespriorm, MIGPRIOR);

            break;
        case MULTPRIOR:
# ifdef USE_MYREAL_FLOAT 
            keys = sscanf(value,"%s%s%f%f%f%f%s%f%f%f%f", priort, part, 
                          &options->bayespriortheta->min, 
                          &options->bayespriortheta->mean, 
                          &options->bayespriortheta->max,
                          &options->bayespriortheta->delta,
                          parm,  &options->bayespriorm->min,                           
                          &options->bayespriorm->mean,
                          &options->bayespriorm->max,
                          &options->bayespriorm->delta);
#else
            keys = sscanf(value,"%s%s%lf%lf%lf%s%lf%lf%lf", priort, part, 
                          &options->bayespriortheta->min, 
                          &options->bayespriortheta->max,
                          &options->bayespriortheta->delta,
                          parm,  &options->bayespriorm->min,                           
                          &options->bayespriorm->max,
                          &options->bayespriorm->delta);
#endif
            if(keys!=9)
                error("abort due to problem with bayes settings in parmfile: error in WEXP prior specification\n");
	    prior_consistency(options->bayespriortheta, THETAPRIOR);
	    prior_consistency(options->bayespriorm, MIGPRIOR);

                break;
    }
    myfree(keeptmp);
}
///
/// checks all elements in the options file (parmfile) that involve boolean yes no type options
/// typical format is option=<NO | YES<: filename or parameters etc>
long
boolcheck (char ch)
{
    char c = uppercase (ch);
    if ((c == 'F') || (c == 'N'))
        return 0;
    else if ((c == 'T') || (c == 'Y'))
        return 1;
    else
        return -1;
}    /* boolcheck */

boolean
booleancheck (option_fmt * options, char *var, char *value)
{
    long i, check;
    char *booltokens[NUMBOOL] = BOOLTOKENS;
    char *tmp;
    long ltemp;
    char *extension;

    check = boolcheck (value[0]);
    if (check == -1)
      {
        return FALSE;
      }
    i = 0;
    while (i < NUMBOOL && strcmp (var, booltokens[i]))
        i++;
    switch ((short) i)
    {
    case 0:   /*menu = <yes | no> */
        options->menu = (boolean) (check);
        break;
    case 1:   /*interleaved =<yes | no> */
        options->interleaved = (boolean) (check);
        break;
    case 2:   /*print-data = <yes | no> */
        options->printdata = (boolean) (check);
        break;
    case 3:   /* mixplot=<yes:mixfilename | no> */
      options->mixplot = set_filename(value, "YES", &options->mixfilename);
        break;
    case 4:   /* moving-steps = <yes | no> */
        options->movingsteps = (boolean) (check);
        if (options->movingsteps)
        {
            strtok (value, ":");
            tmp = strtok (NULL, " ,\n");
            if (tmp != NULL)
                options->acceptfreq = atof ((char *) tmp);
            else
                options->acceptfreq = 0.1;
        }
        break;
    case 5:   /* freqs-from-data =  <yes | no> */
        options->freqsfrom = (boolean) (check);
        if (!options->freqsfrom)
        {
            strtok (value, ":");
            tmp = strtok (NULL, " ,");
            if (tmp != NULL)
                options->freqa = atof ((char *) tmp);
            tmp = strtok (NULL, " ,");
            if (tmp != NULL)
                options->freqc = atof ((char *) tmp);
            tmp = strtok (NULL, " ,");
            if (tmp != NULL)
                options->freqg = atof ((char *) tmp);
            tmp = strtok (NULL, " ,\n");
            if (tmp != NULL)
                options->freqt = atof ((char *) tmp);
        }
        break;
    case 6:   /* useroldtree =  < NO | YES OBSOLETE use usetree*/
        options->usertree = set_filename(value, "YES", &options->utreefilename); //USERTREE
        break;
    case 7:
        {    /* autocorrelation=<YES:value | NO> */
            options->autocorr = (boolean) (check);
            if (options->autocorr)
            {
                strtok (value, ":");
                tmp = strtok (NULL, " ;\n");
                if (tmp != NULL)
                    options->lambda = 1.0 / atof ((char *) tmp);
            }
            break;
        }
    case 8:   /* simulation =  <yes | no> */
        options->simulation = (boolean) check;
        break;
    case 9:   /* plot =  <not | yes:<outfile | both><:std|:log><:{xs,xe,ys,ye}<:N|:M><#of_intervals>>> */
        options->plot = (boolean) check;
        if (options->plot)
        {
            strtok (value, ":");
            if (uppercase (value[0]) == 'Y' || uppercase (value[0]) == 'Y')
            {
                tmp = strtok (NULL, ":;\n");
                if (tmp == NULL)
                    options->plotmethod = PLOTALL;
                else
                {
                    switch (lowercase (tmp[0]))
                    {
                    case 'o':
                        options->plotmethod = PLOTOUTFILE;
                        break;
                    case 'b':
                        options->plotmethod = PLOTALL;
                        break;
                    default:
                        options->plotmethod = PLOTALL;
                        break;
                    }
                    tmp = strtok (NULL, ":;\n");
                    if (tmp != NULL)
                    {
                        switch (lowercase (tmp[0]))
                        {
                        case 'l':
                            options->plotscale = PLOTSCALELOG;
                            break;
                        case 's':
                            options->plotscale = PLOTSCALESTD;
                            break;
                        default:
                            options->plotscale = PLOTSCALELOG;
                        }
                        tmp = strtok (NULL, ":;\n");
#ifdef USE_MYREAL_FLOAT
                        if (4 !=
                                sscanf (tmp, "{%f,%f,%f,%f",
                                        &options->plotrange[0],
                                        &options->plotrange[1],
                                        &options->plotrange[2],
                                        &options->plotrange[3]))
                            sscanf (tmp, "{%f%f%f%f", &options->plotrange[0],
                                    &options->plotrange[1],
                                    &options->plotrange[2],
                                    &options->plotrange[3]);
#else
                        if (4 !=
                            sscanf (tmp, "{%lf,%lf,%lf,%lf",
                                    &options->plotrange[0],
                                    &options->plotrange[1],
                                    &options->plotrange[2],
                                    &options->plotrange[3]))
                            sscanf (tmp, "{%lf%lf%lf%lf", &options->plotrange[0],
                                    &options->plotrange[1],
                                    &options->plotrange[2],
                                    &options->plotrange[3]);
#endif
                        tmp = strtok (NULL, ":;\n");
                        if (tmp != NULL)
                        {
                            switch (lowercase (tmp[0]))
                            {
                            case 'm':
                                options->plotvar = 1;
                                while (!isdigit (*tmp) && *tmp != '\0')
                                    tmp++;
                                if ((ltemp =
                                            strtol (tmp, (char **) NULL, 10)) > 0)
                                    options->plotintervals = ltemp;
                                break;
                            case 'n':
                            default:
                                options->plotvar = PLOT4NM;
                                while (!isdigit (*tmp) && *tmp != '\0')
                                    tmp++;
                                if ((ltemp =
                                            strtol (tmp, (char **) NULL, 10)) > 0)
                                    options->plotintervals = ltemp;

                                break;
                            }
                        }
                    }
                }
            }
        }
        break;
    case 10:   /* weights =  <yes | no> */
        options->weights = (boolean) check;
        set_filename(value, "YES", &options->weightfilename);
        break;
    case 11:   /* read-summary  <yes | no> */
        options->readsum = (boolean) check;
        options->datatype = 'g';
        break;
    case 12:   /* write-summary =  <yes | no> */
        options->writesum = (boolean) check;
        set_filename(value, "YES", &options->sumfilename);
        break;
    case 13: /*unused [was mig-hist before] */
        break;
    case 14:   /* include-unknown=<yes | no> */
        options->include_unknown = (boolean) check;
        break;
        // old case 14(heating) moved to numbercheck
    case 15:   /* print-fst =  <yes | no> */
        options->printfst = (boolean) check;
        break;
    case 16:   /* distfile =  <yes | no> */
        warning("OBSOLETE option distance=YES|NO ===> use usertree=DISTANCE:distfile or usertree=NO\n");
        options->dist = (boolean) check;
        break;
    case 17:   /* geofile =  <yes:geofile | no> */
        options->geo = (boolean) check;
        set_filename(value, "YES", &options->geofilename);
        break;
    case 18:   /* gelman-convergence =  <yes:pairs|sum | no> */
        options->gelman = (boolean) check;
	options->gelmanpairs = FALSE;
        if (options->gelman)
        {
            get_next_word(&value, ":\n", &tmp);
            if (value != NULL)
            {
                if (uppercase (value[0]) == 'P')
		  {
		    options->gelmanpairs = TRUE;
		  }
	    }
	}
        break;
    case 19:   /* randomtree start =  <yes | no> */
        warning("OBSOLETE option randomtree=YES|NO ===> use usertree=RANDOM or usertree=NO\n");
        options->randomtree = (boolean) check;
        break;
    case 20:   /* fast-likelihood calculator =  <yes | no> */
        options->fastlike = (boolean) check;
        break;
    case 21:   /*aic-modeltest=yes/no */
        options->aic = (boolean) check;
        if (options->aic)
        {
            strtok (value, ":\n");
            tmp = strtok (NULL, ":");
            if (tmp != NULL)
            {
                if (uppercase (tmp[0]) == 'F')
                {
                    options->fast_aic = TRUE;
                    tmp = strtok (NULL, ":");
                    if (tmp != NULL)
                    {
                        options->aicmod = atof (tmp);
                    }
                    else
                        options->aicmod = 2.;
                }
                else
                    options->aicmod = atof (tmp);
            }
            else
                options->aicmod = 2.;
        }
        else
            options->aicmod = 2.;
        break;
    case 22:   //use-M=true/false
        options->usem = (boolean) check;
        set_usem_related (options);
        break;
    case 23:   //alternative to above is use-4Nm=false/true
        options->usem = !(boolean) check;
        set_usem_related (options);
        break;
    case 24: //bayes-update
        options->bayes_infer = (boolean) check;
	if(options->bayes_infer)
	  options->lchains = 1;
        break;
    case 25: // bayes-allfile: write all bayes parameter estimates to a file
        options->has_bayesmdimfile = (boolean) check;
	if(options->has_bayesmdimfile)
	  {
	    get_next_word(&value,":",&tmp);// extract the word YES, value is now interval:filename
	    get_next_word(&value,":",&tmp);// extract the interval, value is now filename
	    if(value==NULL)
	      {
		strncpy (options->bayesmdimfilename, tmp, 255);
	      }
	    else
	      {
		options->bayesmdiminterval = atol(tmp);
		strncpy (options->bayesmdimfilename, value, 255);
	      }
	    unpad(options->bayesmdimfilename," ");
	    extension = strrchr(options->bayesmdimfilename,'.');
	    if(extension!=NULL && !strncmp(extension,".gz",3))
	      {
		options->use_compressed = 1;
	      }
	    else
	      {
		options->use_compressed = 0;
	      }  
	  }
        break;
    case 26: // bayes-file: write marginal parameter posterior distribution
        options->has_bayesfile = (boolean) check;
        set_filename(value, "YES", &options->bayesfilename);
        break;
#ifdef NEWVERSION /* divfile datefile */
    case 27:   /* divfile =  <yes:divfile | no> */
        options->div = (boolean) check;
        set_filename(value, "YES", &options->divfilename);
        break;
#endif
    case 28: /* tipdate-file: <yes:datefile | no >  */
        options->has_datefile = (boolean) check;
        set_filename(value, "YES", &options->datefilename);
        break;
    case 29: /* heated-swap: <yes | no >  */
      options->heatedswap_off = !((boolean) check);
      break;
    case 30: // auto-tune=<YES:acceptance-ratio | NO >
      options->has_autotune = (boolean) check;
      if(options->has_autotune)
	{
	  get_next_word(&value,":",&tmp);// extract the word YES, value is now interval:filename
	  get_next_word(&value,":",&tmp);// extract the interval, value is now filename
	  if(value==NULL)
	    {
	      options->autotune = 0.33;
	    }
	  else
	    {
	      options->autotune = atof(value);
	    }
	}
      break;
    default:
        return FALSE;
    }
    return TRUE;
}    /* booleancheck */

void
set_usem_related (option_fmt * options)
{
    if (options->usem)
    {
        options->plotvar = PLOTM;
        options->migvar = PLOTM;
        options->profileparamtype = PLOTM;

    }
    else
    {
        options->plotvar = PLOT4NM;
        options->migvar = PLOT4NM;
        options->profileparamtype = PLOT4NM;
    }

}

boolean
numbercheck (option_fmt * options, char *var, char *value)
{
    MYREAL musum = 0., lastrate = 1.;
    long i = 0, z, cc = 0;
    char *tmp, *temp, *temp2, *keeptmp;
    char *numbertokens[NUMNUMBER] = NUMBERTOKENS;

    tmp = (char *) mycalloc (LINESIZE, sizeof (char));
    keeptmp = tmp;
    while (i < NUMNUMBER && strcmp (var, numbertokens[i]))
        i++;
    switch ((short) i)
    {
    case 0:   /*ttratio = value */
        z = 0;
        temp = strtok (value, " ,;\n\0");
        while (temp != NULL)
        {
            options->ttratio[z++] = atof (temp);
            options->ttratio =
                (MYREAL *) myrealloc (options->ttratio, sizeof (MYREAL) * (z + 1));
            options->ttratio[z] = 0.0;

            temp = strtok (NULL, " ,;\n\0");
        }
        break;
    case 1:   /*short-chains = value */
        options->schains = atol (value);
        break;
    case 2:   /*short-steps = value */
    case 32:   /* short-sample = value */
        options->ssteps = atol (value);
        break;
    case 3:   /*short-increment = value */
        options->sincrement = atol (value);
        break;
    case 4:   /*long-chains = value */
        options->lchains = atol (value);
        break;
    case 5:   /*long-steps = value */
    case 33:   /*long-sample = value */
        options->lsteps = atol (value);
        break;
    case 6:   /*long-increment = value */
        options->lincrement = atol (value);
        break;
    case 7:
        break;   /* theta: already handled in read_theta() */
    case 8:   /*nmlength = value */
        options->nmlength = strtol (value, (char **) NULL, 10); //atoi (value);
        break;
    case 9:   /* seed = <Auto | seedfile | Own:value> */
        switch (value[0])
        {
        case 'A':
        case 'a':
        case '0':
            options->autoseed = AUTO;
            options->inseed = (long) time (0) / 4 + 1;
            break;
        case 'S':
        case 's':
        case '1':
            options->autoseed = NOAUTO;
            openfile(&options->seedfile, options->seedfilename, "r", NULL);
            if (options->seedfile == NULL)
            {
                usererror ("cannot find seedfile\n");
            }
            fscanf (options->seedfile, "%ld%*[^\n]", &options->inseed);
            fclose (options->seedfile);
            break;
        case 'O':
        case 'o':
        case '2':
            options->autoseed = NOAUTOSELF;
            strtok (value, ":");
            tmp = strtok (NULL, " ;\n");
            if (tmp != NULL)
                options->inseed = atol ((char *) tmp);
            if (options->inseed > 0)
                break;
        default:
            options->autoseed = AUTO;
            options->inseed = (long) time (0) / 4 + 1;
            usererror ("Failure to read seed method, should be\n \
                       random-seed=auto or random-seed=seedfile or random-seed=own:value\nwhere value is a positive integer\nUsing AUTOMATIC seed=%li\n", options->inseed);
            break;
        }
        break;
    case 10:
        break;   /*"migration" fake: this is already handled in read_migrate */
    case 11:   /*mutation= <auto=gamma | nogamma | constant | estimate | own | data> */
        switch (value[0])
        {
        case 'A':  /*automatic */
        case 'a':
        case 'E':  /*automatic */
        case 'e':
	    options->bayesmurates=TRUE;
            options->murates = FALSE;
            options->murates_fromdata = FALSE;
	    options->gamma = FALSE;
	    break;
        case 'G':
        case 'g':
            options->gamma = TRUE;
	    options->bayesmurates=TRUE;
            options->murates = FALSE;
            options->murates_fromdata = FALSE;
            (void) strtok (value, " :");
            temp = strtok (NULL, " ,;\n");
            if (temp != NULL)
                options->alphavalue = atof (temp);
            else
                options->alphavalue = START_ALPHA;
            break;
        case 'O':
        case 'o':
            options->murates = TRUE;
            options->gamma = FALSE;
            options->murates_fromdata = FALSE;
	    options->bayesmurates = FALSE;
            (void) strtok (value, " :");
            temp = strtok (NULL, " ,;\n");
            if (temp != NULL)
            {
                options->muloci = atol (temp);
                options->mu_rates =
                    (MYREAL *) mycalloc (options->muloci, sizeof (MYREAL));
                musum = 0.;
                for (i = 0; i < options->muloci; i++)
                {
                    temp = strtok (NULL, " ,;\n");
                    if (temp == NULL)
                    {
                        while (i < options->muloci)
                        {
                            options->mu_rates[i] = lastrate;
                            musum += options->mu_rates[i];
                            i++;
                        }
                    }
                    lastrate = options->mu_rates[i] = atof (temp);
                    musum += options->mu_rates[i];
                }
                // mean must be 1.
                musum /= options->muloci;
                for (i = 0; i < options->muloci; i++)
                {
                    options->mu_rates[i] /= musum;
                }
            }
            break;
	case 'D': /* mutaton is varying and calculated from the data*/ 
	case 'd':
            options->gamma = FALSE;
            options->murates = TRUE;
            options->murates_fromdata = TRUE;
	    options->bayesmurates = FALSE;
	    break;
        case 'N':  /*nogamma, none, all loci have same mu */
        case 'n':
	case 'C':
	case 'c':
        default:
            options->murates = FALSE;
            options->gamma = FALSE;
            options->murates_fromdata = FALSE;
	    options->bayesmurates = FALSE;
            break;
        }
        break;
    case 12:   /*datatype=<allele|microsatellite|brownian|sequence|f-ancestral states|genealogies> */
        switch (value[0])
        {
        case 'a':
        case 'A':
            options->datatype = 'a';
            break;
        case 'm':
        case 'M':
            options->datatype = 'm';
            break;
        case 'b':
        case 'B':
            options->datatype = 'b';
            break;
        case 's':
        case 'S':
            options->datatype = 's';
            break;
        case 'n':
        case 'N':
            options->datatype = 'n';
            break;
        case 'h':
        case 'H':
            options->datatype = 'h';
	    options->fastlike=FALSE;
            break;
        case 'u':
        case 'U':
            options->datatype = 'u';
            break;
        case 'f':
        case 'F':
            options->datatype = 'f';
            break;
        case 'g':
        case 'G':
            options->datatype = 'g';
            options->readsum = TRUE;
            break;
        default:
            options->datatype = 's';
            break;
        }
        break;
    case 13:   /* categories=<None | value> */
        if (uppercase (value[0] == 'N'))
        {
            options->categs = ONECATEG;
            break;
        }
        else
        {
            options->categs = strtol (value, (char **) NULL, 10);
            /* needs to read auxilliary file catfile */
            sprintf(tmp,"%li",options->categs);
            set_filename(value, tmp, &options->catfilename);
        }
        break;
    case 14:   /*create rates=value:list of rates */
        strncpy (tmp, value, strcspn (value, ":"));
        if (strtol (tmp, (char **) NULL, 10) /*;atoi (tmp) */  > 1)
        {   /* rate categories */
            options->rcategs = strtol (tmp, (char **) NULL, 10);
            options->rrate =
                (MYREAL *) myrealloc (options->rrate,
                                    sizeof (MYREAL) * (options->rcategs + 1));
            (void) strtok (value, " :");
            temp = strtok (NULL, " ,;\n");
            z = 0;
            while (temp != NULL)
            {
                if (z > options->rcategs)
                    usererror ("check parmfile-option  rates, missing rate\n");
                options->rrate[z++] = atof (temp);
                temp = strtok (NULL, " ,;\n");
            }
        }
        break;
    case 15:   /* probabilities for each rate category */
        strncpy (tmp, value, strcspn (value, ":"));
        if (strtol (tmp, (char **) NULL, 10) > 1)
        {   /* probabilities for each rate category */
            options->rcategs = strtol (tmp, (char **) NULL, 10);
            options->probcat =
                (MYREAL *) myrealloc (options->probcat,
                                    sizeof (MYREAL) * (options->rcategs + 1));
            (void) strtok (value, " :");
            temp = strtok (NULL, " ,;\n");
            z = 0;
            while (temp != NULL)
            {
                if (z > options->rcategs)
                    usererror
                    ("check parmfile prob-rates, missing rate probability\n");
                options->probcat[z++] = atof (temp);
                temp = strtok (NULL, " ,;\n");
            }
        }
        break;
    case 16:   /*micro-stepmax */
        // options->micro_stepnum = strtol (value, (char **) NULL, 10);

        break;
    case 17:   /*micro-threshold */
      options->micro_threshold = atol(value);
      if(options->micro_threshold % 2 != 0)
	options->micro_threshold += 1;
      break;
    case 18:   /*delimiter */
      options->dlm = value[0];
      break;
    case 19:   /*burn-in */
      options->burn_in = atol (value);
      if(strchr(value,'a')!=NULL)
	options->burnin_autostop = 'a';
      else
	{
	  if(strchr(value,'e')!=NULL)
	    options->burnin_autostop = 'e';
	  else
	    options->burnin_autostop = ' ';
	}
      break;
    case 20:   /*infilename */
      strcpy (options->infilename, value);
      break;
    case 21:   /*outfilename */
      strcpy (options->outfilename, value);
        break;
    case 22:   /*mathfilename */
        strcpy (options->mathfilename, value);
        break;
    case 23:   /*title */
        strncpy (options->title, value, 80);
        break;
    case 24:   /*long-chain-epsilon */
        options->lcepsilon = atof (value);
        if (options->lcepsilon <= 0)
            options->lcepsilon = LONGCHAINEPSILON;
        break;
    case 25:   /* print tree options */
        switch (uppercase (value[0]))
        {
        case 'N':
            options->treeprint = _NONE;
            break;
        case 'A':
            options->treeprint = ALL;
	    options->treeinc = 1;
	    set_filename(value, "ALL", &options->treefilename);
            break;
        case 'B':
            options->treeprint = BEST;
	    set_filename(value, "BEST", &options->treefilename);
            break;
        case 'L':
            options->treeprint = LASTCHAIN;
	    get_next_word(&value,":",&tmp);// the word LASTCHAIN
	    get_next_word(&value,":",&tmp);
	    options->treeinc = atol(tmp);
	    set_filename(value, "LASTCHAIN", &options->treefilename);
            break;
        default:
            options->treeprint = _NONE;
            break;
        }
        break;
    case 26:   /* progress: No, Yes, Verbose */
        switch (uppercase (value[0]))
        {
        case 'F':
        case 'N':
            options->progress = FALSE;
            options->verbose = FALSE;
            break;
        case 'T':
        case 'Y':
            options->progress = TRUE;
            options->verbose = FALSE;
            break;
        case 'V':
            options->progress = TRUE;
            options->verbose = TRUE;
            break;
        }
        break;
    case 27:   /* l-ratio: <NO | YES>:val1,val2,val3,val4,val5 */
        cc = options->lratio->counter;
        switch (uppercase (value[0]))
        {
        case 'Y':
            options->lratio->data[cc].type = MLE;
            break;
        case 'N':
        default:
            myfree(keeptmp);
            return FALSE;
        }
        (void) strtok (value, ":");
        temp = strtok (NULL, "\n");
        if (temp != NULL)
        {
            temp2 = strchr (temp, ':');
            if (temp2 != NULL)
            {
                strcpy (options->lratio->data[cc].value2, temp2);
                *temp2 = '\0';
                strcpy (options->lratio->data[cc].value1, temp);
            }
            else
                strcpy (options->lratio->data[cc].value1, temp);
        }
        if (cc + 1 == options->lratio->alloccounter)
        {
            options->lratio->alloccounter += 2;
            options->lratio->data =
                (lr_data_fmt *) myrealloc (options->lratio->data,
                                         sizeof (lr_data_fmt) *
                                         (options->lratio->alloccounter+1));
            for (i = cc + 1; i < options->lratio->alloccounter; i++)
            {
                options->lratio->data[i].elem = 0;
                options->lratio->data[i].value1 =
                    (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
                options->lratio->data[i].elem = 0;
                options->lratio->data[i].value2 =
                    (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
                options->lratio->data[i].connect =
                    (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);

            }
        }
      	//xcode     cc = ++options->lratio->counter;
        break;
    case 28:   /* fst-type: <Theta | Migration> */
        switch (uppercase (value[0]))
        {
        case 'T':
            options->fsttype = 'T';
            break;
        case 'M':
        default:
            options->fsttype = 'M';
            break;
        }
        fst_type (options->fsttype);
        break;
    case 29:   /*profile=<NO| NONE | YES | ALL | TABLES | SUMMARY>><: <FAST |  */
        switch (uppercase (value[0]))
        {
        case 'S':
            options->profile = SUMMARY;
            break;
        case 'Y':
        case 'A':
            options->profile = ALL;
            break;
        case 'N':
            options->profile = _NONE;
            break;
        case 'T':
            options->profile = TABLES;
            break;
        default:  /*A */
            options->profile = ALL;
            break;
        }
        (void) strtok (value, ":;\n");
        temp = strtok (NULL, ":;\n");
        if (temp != NULL)
        {
            switch (lowercase (temp[0]))
            {
            case 'p':  /*precise percentiles */
                options->profilemethod = 'p';
                break;
            case 'd':  /*discrete steps see at start of file */
                options->profilemethod = 'd';
                break;
                //case 's':
                //options->profilemethod = 's';
                //break;
            case 'x':  /* x-rated */
            case 'u':  /* uncorrelated */
            case 'q':  /* quick and dirty */
                options->profilemethod = 'q';
                break;
            case 'f':  /* quick and exact mixture */
                options->profilemethod = 'f';
                break;
            default:
                options->profilemethod = 'f';
                options->printprofsummary = TRUE;
                break;
            }
            temp = strtok (NULL, ":;\n");
            if (temp != NULL)
            {
                switch (lowercase (temp[0]))
                {
                case 'm':
                    options->profileparamtype = 1;
                    break;
                default:
                    options->profileparamtype = PLOT4NM;
                }
            }
        }
        set_profile_options (options);
        break;
    case 30:   /* custom-migration:<{> migration matrix and theta on
                           diagonal:
                           0 means not estimated,
                           x means estimated, s means symmetrically
                           estimated, m means all are the same  <}> */
        if (myID == MASTER)
            read_custom_migration (options->parmfile, options, value,
                                   options->numpop);
#ifdef MPI

        else
            read_custom_migration_worker (options->buffer, options, value,
                                          options->numpop);
#endif

        break;
    case 31:   /*sumfilename */
        strcpy (options->sumfilename, value);
        break;
        /*case 32 and case 33 are fallthroughs to 2 and 3 */
    case 34:   /*replicate */
        switch (uppercase (value[0]))
        {
        case 'T':
        case 'Y':
            options->replicate = TRUE;
            temp = strtok (value, ":;\n");
            if (temp != NULL)
            {
                temp = strtok (NULL, ":;\n");
                if (uppercase (temp[0]) == 'L')
                    options->replicatenum = 0;
                else
                    options->replicatenum = strtol (temp, (char **) NULL, 10);
            }
            else
                options->replicatenum = 0;
            break;
        default:
            options->replicate = FALSE;
            options->replicatenum = 0;
        }
        break;
    case 35:   /* cpu number */
        options->cpu = (short) ATOI (value);
        break;
    case 36:   /* logfile=<YES:logfile | NO> do we write a logfile or not */
        options->writelog = set_filename(value,"YES",&options->logfilename);
        break;
    case 37:   /* sequencing error */
        options->seqerror = atof (value);
        if (options->seqerror < 0.0)
        {
            warning
            ("Sequencing error was misspecified in parmfile, reset to 0.0\n");
            options->seqerror = 0.0;
        }
        break;
#ifdef UEP

    case 38:   /* do we have a uep file or not, function returns TRUE if uep=YES:filename*/
        options->uep = set_filename(value, "YES", &options->uepfilename);
        break;
    case 39:   /* do we have uep-rates */
        temp = strtok (value, ":;\n");
        if (temp != NULL)
        {
            options->uepmu = atof (temp);
            temp = strtok (NULL, ":;\n");
            if (temp != NULL)
            {
                options->uepnu = atof (temp);
            }
            options->ueprate = options->uepmu;
        }
        break;
    case 40:   /* do we have uep-bases */
        temp = strtok (value, ":; \n");
        if (temp != NULL)
        {
            options->uepfreq1 = atof (temp);
            temp = strtok (NULL, " :;\n");
            if (temp != NULL)
            {
                options->uepfreq0 = atof (temp);
            }
            else
                options->uepfreq0=1. - options->uepfreq1;
        }
        break;
#endif

    case 41: // mu-rates??????
        break;
    case 42: /* heating=<no | <yes | adaptive | bounded>:numintervals:{temperatures}> */
        switch (uppercase (value[0]))
        {
        case 'A': //adaptive heating on
            options->heating = 1;
            options->adaptiveheat = STANDARD;
            break;
        case 'B': //adaptive heating on
            options->heating = 1;
            options->adaptiveheat = BOUNDED;
            break;
        case 'Y':
        case 'P':
            options->heating = 1;
            options->adaptiveheat = NOTADAPTIVE;
            break;
        case 'N':
        default:
            options->heating=0;
            options->adaptiveheat=NOTADAPTIVE;
            break;
        }
        if (options->heating == 1)
        {
            strtok (value, ":\n");
            tmp = strtok (NULL, ": ");
            if (tmp != NULL)
            {
                options->heating_interval = atol (tmp);
                tmp = strtok (NULL, "{, ");
                if (tmp != NULL)
                {
                    z = 0;
                    while (1)
                    {
                        options->heat[z++] = atof (tmp);
                        tmp = strtok (NULL, ", :\n");
                        if (tmp == NULL || z >= 1000)
			  {
			    if(z>=1000)
			      warning("Ignored the all temperature above 1000 heated chains\n");
                            break;
			  }
                    }
                    options->heated_chains = z;
                }
            }
        }
        break;
#ifdef LONGSUM

    case 43: /* fluctuate=<no | <yes>:{rate_1today, time1today,rate_1middle,time1middle, rate_1past, time1past,rate2_today,time2_today,...}> */
        switch (uppercase (value[0]))
        {
        case 'Y':
        case 'P':
            options->fluctuate = TRUE;
            break;
        case 'N':
        default:
            options->fluctuate = FALSE;
            break;
        }
        if (options->fluctuate)
        {
            strtok (value, ":\n");
            tmp = strtok (NULL, ":{,");
            if (tmp != NULL)
            {
                z = 0;
                while (1)
                {
                    options->flucrates = myrealloc(options->flucrates, sizeof(MYREAL) * (z+1));
                    options->flucrates[z++] = atof (tmp);
                    tmp = strtok (NULL, ", {}:\n");
                    if (tmp == NULL)
                        break;
                }
            }
        }
        if(z%3!=0)
            error("Fluctuating rates and times need to be 3 rates and 3 times per population");
        break;
#endif /*LONGSUM*/
    case 44:   /*resistance = value  [to fatal attraction to zero*/
    	//xcode       z = 0;
        temp = strtok (value, " ,;\n\0");
        if (temp != NULL)
            options->minmigsumstat = atof (temp);
        else
            options->minmigsumstat = MINMIGSUMSTAT;
        break;
    case 45: /* bayes-updatefreq=[0..1] */
        temp = strtok (value, " ,;\n\0");
        if (temp != NULL)
            options->updateratio = atof (temp);
        else
            options->updateratio= HALF;
        break;
    case 49: /*bayes-posteriorbins*/    
        temp = strtok (value, " ,;\n\0");
        if (temp != NULL)
        {
            options->bayespriortheta->bins = atol (temp);
            temp = strtok (NULL, " ,;\n\0");
            if (temp != NULL)
                options->bayespriorm->bins = atol (temp);
            else
                options->bayespriorm->bins = BAYESNUMBIN;
        }
        else
        {
            options->bayespriortheta->bins = BAYESNUMBIN;
            options->bayespriorm->bins = BAYESNUMBIN;
        }
        break;
        
    case 46:   /*bayesfile=FILENAME  set bayesfile and bayes analysis? */
      strncpy (options->bayesfilename, value, 255); //this is legacy code and got replaced
      // by bayes-file=<YES | NO>:FILENAME
      options->has_bayesfile = TRUE;
	    break;
		
    case 47:   /*set bayes-prior*/
        switch(toupper(value[0]))
           {
	     //case 'S':
	     //set_bayes_options_oldstyle(value, options,SLICE);
	     //break;
	   case 'M':
                set_bayes_options_oldstyle(value, options,MULTPRIOR);
                break;
            case 'E':
                    set_bayes_options_oldstyle(value, options,EXPPRIOR);
                    break;
                case 'W':
                    set_bayes_options_oldstyle(value, options,WEXPPRIOR);
                    break;
                case 'U':
                default:
                    set_bayes_options_oldstyle(value, options,UNIFORMPRIOR);
                    break;
           }
        break;
    case 48:   /* usertree =  < NO |  UPGMA | AUTOMATIC | TREE:intreefilename | RANDOM | DISTANCE:distfilename */
        options->usertree = set_filename(value, "T", &options->utreefilename); //USERTREE
        options->randomtree = (value[0] =='R') ? TRUE : FALSE ; //RANDOMTREE
        options->dist = set_filename(value, "D", &options->distfilename); //DISTANCE
        break;
	/* case 49 is further back*/
    case 50: /* mig-histogram=<<yes|all>:histogram-binsize:filename | no> */
      get_next_word(&value,":",&tmp); 
      options->mighist_all = ((uppercase(tmp[0]) == 'A') ? TRUE : FALSE);
      options->mighist     = ((uppercase(tmp[0]) == 'Y') ? TRUE : FALSE);
      if(options->mighist_all || options->mighist)
	{
	  options->mighist = TRUE;
	  get_next_word(&value,":",&tmp);
	  if(value!=NULL)
	    {
	      options->eventbinsize = atof(tmp);
	      strncpy (options->mighistfilename, value, 255);
	    }
	  else
	    {
	      strncpy (options->mighistfilename, tmp, 255);
	    }
	}
      break;
    case 51:
        switch(toupper(value[0]))
           {
	   case 'A':
	     options->bayespretty = PRETTY_MAX;
	     break;
	   case 'P':
	     options->bayespretty = PRETTY_P99;
	     break;
	   case 'M':
	     options->bayespretty = PRETTY_P99MAX;
	     break;
	   case 'T':
	   default:
	     options->bayespretty = PRETTY_P100;
	     break;
	   }
	break;
#ifdef PRETTY
    case 52:   /*PDF-outfile name */
        strcpy (options->pdfoutfilename, value);
        break;
#endif
    case 53:   /*bayes interval for writing all parameters to file */
        options->bayesmdiminterval = atol(value);
        break;
    case 54:   /*set bayes-priorS*/
      // new parmfile variable so that the old stype setting still work
      set_bayes_options(value, options);
      break;
    case 55: /*set skyline-histogram*/
      get_next_word(&value,":",&tmp); 
      options->skyline = tmp[0] == 'Y' ? TRUE : FALSE;
      if(options->skyline)
	{
	  get_next_word(&value,":",&tmp);
	  options->eventbinsize = atof(tmp);
	  strncpy (options->skylinefilename, value, 255);
	}
      break;
    case 56: /* rates-gamma */
      get_next_word(&value,":,; ",&tmp);
      options->seqrate_gamma_num = 0;
      while(tmp!=NULL)
	{
	  options->seqrate_gamma[options->seqrate_gamma_num++] = atof(tmp);
	}
      break;
    case 57: /* Bayes-proposals*/
      get_next_word(&value,":,; ",&tmp);
      if(tmp != NULL)
	{
	  switch(uppercase(tmp[0]))
	    {
	    case 'T':       get_next_word(&value,":,; ",&tmp);
	      if(tmp != NULL)
		{
		  if(uppercase(tmp[0])=='S')
		    options->slice_sampling[THETAPRIOR] = TRUE;
		  else
		    options->slice_sampling[THETAPRIOR] = FALSE;
		}
	      break;
	    case 'M':        get_next_word(&value,":,; ",&tmp);
	      if(tmp != NULL)
		{
		  if(uppercase(tmp[0])=='S')
		    options->slice_sampling[MIGPRIOR] = TRUE;
		  else
		    options->slice_sampling[MIGPRIOR] = FALSE;
		}
	      break;
	    case 'R':       get_next_word(&value,":,; ",&tmp);
	      if(tmp != NULL)
		{
		  if(uppercase(tmp[0])=='S')
		    options->slice_sampling[RATEPRIOR] = TRUE;
		  else
		    options->slice_sampling[RATEPRIOR] = FALSE;
		}
	      break;
	    }
	}
      break;
    case 58: /*generation-per-year*/
      get_next_word(&value,":,; ",&tmp);
      options->generation_year = atof(tmp);
      break;
    case 59: /*mutationrate-per-year absolute mutation rate per year for each locus*/
      options->mutationrate_year[0] = 1.0;
      get_next_word(&value,"{}:,; ",&tmp);
      z=0;
      while(tmp != NULL)
	{
	  if(z >= options->mutationrate_year_numalloc)
	    {
	      options->mutationrate_year_numalloc = z+1;
	      options->mutationrate_year = (MYREAL*) myrealloc(options->mutationrate_year, options->mutationrate_year_numalloc * sizeof(MYREAL));
	    }
	  options->mutationrate_year[z] = atof(tmp);
	  printf("%i> mutationrate_year[%li]=%g %s\n",myID,z, options->mutationrate_year[z],tmp);
	  z++;
	  get_next_word(&value,"{}:,; ",&tmp);
	}
      break;
#ifdef NEWVERSION /*inheritance scalar*/
    case 60: /*inheritance-scalars defines the scalars so that all reference*/
      /* is made to 4 Ne mu when the scalar is used, otherwise the scalar is*/
      /* assume to be 1*/
      options->inheritance_scalars[0] = 1.0;
      get_next_word(&value,"{}:,; ",&tmp);
      z=0;
      while(tmp != NULL)
	{
	  if(z >= options->inheritance_scalars_numalloc)
	    {
	      options->inheritance_scalars_numalloc = z+1;
	      options->inheritance_scalars = (MYREAL*) myrealloc(options->inheritance_scalars, options->inheritance_scalars_numalloc * sizeof(MYREAL));
	    }
	  options->inheritance_scalars[z] = atof(tmp);
	  //printf("%i> inheritance scalar[%li]=%g %s\n",myID,z, options->inheritance_scalars[z],tmp);
	  z++;
	  get_next_word(&value,"{}:,; ",&tmp);
	}
      break;
#endif /*NEWVERSION inheritance scalar*/
    case 61: /* micro-submodel */
      get_next_word(&value,":",&tmp); 
      if(tmp==NULL)
	options->msat_option = atol(value);
      else
	{
	  options->msat_option = atol(tmp);
	  if(value==NULL)
	    break;
	    get_next_word(&value,", ",&tmp);
	  if(tmp[0]=='{')
	    options->msat_tuning[0] = atof(tmp+1);
	  else
	    options->msat_tuning[0] = atof(tmp);
	  get_next_word(&value," }",&tmp);
	  options->msat_tuning[1] = atof(tmp);
	}
      break;

    case 62: /*random-subset*/
      get_next_word(&value,":,; ",&tmp);
      options->randomsubset = atol(tmp);
      break;

    case 63: /* pooling of populations an array of numbers indicates the pooling*/
      set_localities(&value, &tmp, options);
      /*      options->newpops[0] = 1;
      get_next_word(&value,"{}:,; ",&tmp);
      z=0;
      while(tmp != NULL)
	{
	  if(z >= options->newpops_numalloc)
	    {
	      options->newpops_numalloc = z+1;
	      options->newpops = (long*) myrealloc(options->newpops, options->newpops_numalloc * sizeof(long));
	    }
	  options->newpops[z] = atol(tmp);
	  //printf("%i> population relabel [%li]=%li (input: |%s|)\n",myID,z, options->newpops[z],tmp);
	  z++;
	  get_next_word(&value,"{}:,; ",&tmp);
	  }*/
      break;
    default:
      myfree(keeptmp);
      return FALSE;      
    }
    myfree(keeptmp);
    
    return TRUE;
}    /* numbercheck */

/*void
reset_oneline (option_fmt * options, long position)
{
    fseek (options->parmfile, position, SEEK_SET);
    }*/

void
reset_oneline (option_fmt * options, fpos_t * position)
{
    fsetpos (options->parmfile, position);
}

void read_random_theta(option_fmt *options, char ** buffer)
{
   char *tempstr;
   char *tmp;
   char *otempstr;
   char *otmp;
   
   tempstr = (char *) mycalloc(LINESIZE,sizeof(char));
   tmp = (char *) mycalloc(LINESIZE,sizeof(char));
   otempstr = tempstr;
   otmp = tmp;

  options->numthetag = 2;
  options->thetag =
    (MYREAL *) myrealloc (options->thetag, sizeof (MYREAL) * 3);
  if( myID == MASTER)
    {
      fgets (tmp, LINESIZE, options->parmfile);
    }
  else
    {
      sgets (tmp, LINESIZE, buffer);
    }
  get_next_word(&tmp,",{} \n",&tempstr);
  while(strchr("{ ", (tempstr)[0]))
    {
      get_next_word(&tmp,",{} \n",&tempstr);
    }
  options->thetag[0] = atof(tempstr);
  get_next_word(&tmp,",{} \n",&tempstr);
  options->thetag[1] = atof(tempstr);
  myfree(otmp);
  myfree(otempstr);
}



void read_random_mig(option_fmt *options, char ** buffer)
{
   char *tempstr;
   char *tmp;
   char *otempstr;
   char *otmp;
   
   tempstr = (char *) mycalloc(LINESIZE,sizeof(char));
   tmp = (char *) mycalloc(LINESIZE,sizeof(char));
   otempstr = tempstr;
   otmp = tmp;
  options->nummg = 2;
  options->mg = (MYREAL *) myrealloc (options->mg, sizeof (MYREAL) * 3);
  if(myID == MASTER)
    {
      fgets (tmp, LINESIZE, options->parmfile);
    }
  else
    {
      sgets (tmp, LINESIZE, buffer);
    }
  get_next_word(&tmp,",{} \n",&tempstr);
  //       	while(strchr("{ ", tempstr[0]))
  //  tempstr = strsep(&tmp,",{} \n");
  options->mg[0] = atof(tempstr);
  get_next_word(&tmp,",{} \n",&tempstr);
  options->mg[1] = atof(tempstr);
  myfree(otmp);
  myfree(otempstr);
}

///
/// read the values for starte theta from the list in the parmfile
/// called through read_options_master()
void read_theta (option_fmt * options)
{
   long i = 0;
   
   char varvalue[LINESIZE];
   char ch;
   
   char *tempstr;
   char *tmp;
   char *otempstr;
   char *otmp;
   
   tempstr = (char *) mycalloc(LINESIZE,sizeof(char));
   tmp = (char *) mycalloc(LINESIZE,sizeof(char));
   otempstr = tempstr;
   otmp = tmp;
//xcode removed ch out of while
    while ((getc (options->parmfile)) != '=')
        ;
    ch = getc (options->parmfile);
    while (!isspace ((int) ch) && ch != ':' && ch != '{')
    {
        varvalue[i++] = ch;
        ch = getc (options->parmfile);
    }
    varvalue[i]='\0';
    //    printf("%s, %c\n",varvalue,ch);
    switch (uppercase (varvalue[0]))
    {
    case 'N':
        options->thetaguess = NRANDOMESTIMATE;
	read_random_theta(options, NULL);
    case 'U':
        options->thetaguess = URANDOMESTIMATE;
	read_random_theta(options, NULL);
        break;
    case 'F':
    case '_':
        options->thetaguess = FST;
        break;
        //case 'G':
    case 'O':
    case '0':
        options->thetaguess = OWN;
        ch = skip_space (options);
        if (ch == '\0')
            return;
        if (ch == '{')
        {
            ch = skip_space (options);
            while (ch != '}')
            {
                i = 0;
                if (ch == '\0')
                    return;
                while (!strchr(" ,\t}", ch))
                {
                    tmp[i++] = ch;
                    ch = getc (options->parmfile);
                }
                tmp[i] = '\0';
                options->thetag[options->numthetag] = atof (tmp);
                options->numthetag += 1;
                options->thetag = (MYREAL *) myrealloc (options->thetag, sizeof (MYREAL) * (1 + options->numthetag));
                ch = skip_space (options);
            }
        }
        else
        {
            i = 0;
            while (!isspace (ch))
            {
                tmp[i++] = ch;
                ch = getc (options->parmfile);
            }
            tmp[i] = '\0';
            options->thetag[options->numthetag] = atof (tmp);
            options->numthetag += 1;
            options->thetag =
                (MYREAL *) myrealloc (options->thetag,
                                    sizeof (MYREAL) * (1 + options->numthetag));
        }
        options->numpop = options->numthetag;
        if (options->numthetag == 0)
        {
            warning ("You forgot to add your guess value:\n");
            warning ("Theta=Own:{pop1,pop2, ...}\n");
            warning ("or Theta=Own:guess_pop (same value for all)\n");
        }
        break;
    default:
#ifdef MPI
        fprintf(stderr,"%i> died in read_theta()\n",myID);
#endif
        usererror
        ("Failure to read start theta method, should be\ntheta=FST or theta=Own:x.x\n or theta=Own:{x.x, x.x , x.x, .....} or theta=NRANDOM:{mean,std} or URANDOM:{min,max}");
    }
    myfree(otmp);
    myfree(otempstr);
}


void
read_mig (option_fmt * options)
{
    //  char parmvar[LINESIZE];
    long test = 0;
    long i = 0;

    char varvalue[LINESIZE];
    char ch;
 
   char *tempstr;
   char *tmp;
   char *otempstr;
   char *otmp;
    int choint;
    char choice;
   varvalue[0]='#';// to keep the static analyzer happy

   tempstr = (char *) mycalloc(LINESIZE,sizeof(char));
   tmp = (char *) mycalloc(LINESIZE,sizeof(char));
   otempstr = tempstr;
   otmp = tmp;

    /* 1st example:  1.0 (n-island model)
       2nd example: {1.0} (migration matrix model, all the same start values  
       3rd example: the dashes on the diagonal are NECESSARY, {} are facultativ
       -  1.0 0.1
       1.0  -  2.0
       0.9 1.2  -
       to specify real 0.0 you need to use the custom-migration settings.
       0.0 in the table will be change to SMALLES_MIGRATION
     */
//xcode removed ch out of while
   while ((getc (options->parmfile)) != '=');
   //    {
        //      parmvar[i++] = ch;
   // }
    i = 0;
    ch = getc (options->parmfile);
    while (!isspace ((int) ch) && ch != ':' && ch != '{')
    {
        varvalue[i++] = ch;
        ch = getc (options->parmfile);
    }
    choint = (int) varvalue[0];
    choice = uppercase (choint);
    switch (choice)
    {
    case 'N':
        options->migrguess = NRANDOMESTIMATE;
	read_random_mig(options,NULL);
        break;
    case 'U':
        options->migrguess = URANDOMESTIMATE;
	read_random_mig(options,NULL);
        break;
    case 'F':
    case '_':
        options->migrguess = FST;
        break;
        //    case 'G':
    case 'O':
    case '0':
        //      if(varvalue[0]=='G')
        //       options->migrguess = PARAMGRID;
        //      else
        options->migrguess = OWN;
        ch = skip_space (options);
        if (ch == '\0')
            return;
        if (ch == '{')
        {
            options->migration_model = MATRIX;
            ch = skip_space (options);
            while (ch != '}')
	      {
                if ((ch == '\0') || (ch == '#'))
		  return;
                i = 0;
                while (!strchr(" ,\t}",ch))
		  {
                    tmp[i++] = ch;
                    ch = getc (options->parmfile);
		    if(ch == '#')
		      {
			while(ch!='\n' || ch != '\0')
			  ch = getc(options->parmfile);
			return;
		      }
		  }
                tmp[i] = '\0';
                if (strcmp (tmp, "-"))
		  {
                    options->mg[options->nummg] = atof (tmp);
                    options->nummg += 1;
                    options->mg =
		      (MYREAL *) myrealloc (options->mg,
                                            sizeof (MYREAL) * (1 +
                                                               options->nummg));
		  }
                else
		  {
                    test++;
		  }
                ch = skip_space (options);
	      }
            options->numpop = test;
        }
        else
        {
            options->migration_model = ISLAND;
            i = 0;
            options->numpop = 1;
            while (!isspace (ch))
            {
                tmp[i++] = ch;
                ch = getc (options->parmfile);
            }
            options->mg[options->nummg] = atof (tmp);
            options->nummg += 1;
        }
        if (options->nummg == 0)
        {
            warning ("You forgot to add your guess value, use either:\n");
            warning ("migration=FST\n");
            warning ("migration=Own:{migration matrix, diagonal is -}\n");
            usererror
            ("migration=Own:migration_value, all matrix elements have the same value\n");
        }
        break;
    default:
        usererror ("Failure to read start migration method\n");
    }
    myfree(otmp);
    myfree(otempstr);
}

#ifdef MPI
void
read_theta_worker (char **buffer, option_fmt * options)
{
    long i = 0;

    char varvalue[LINESIZE];
    char ch;

    //char *tempstr;
    char *tmp;
    //char *otempstr;
    char *otmp;

    //tempstr = (char *) mycalloc(LINESIZE,sizeof(char));
    tmp = (char *) mycalloc(LINESIZE,sizeof(char));
    //otempstr = tempstr;
    otmp = tmp;

    while ((ch = sgetc (buffer)) != '=')
        ;
    ch = sgetc (buffer);
    while (!isspace ((int) ch) && ch != ':' && ch != '{')
    {
        varvalue[i++] = ch;
        ch = sgetc (buffer);
    }
    switch (uppercase (varvalue[0]))
    {
    case 'N':
        options->thetaguess = NRANDOMESTIMATE;
	read_random_theta(options, buffer);
        break;
    case 'U':
        options->thetaguess = URANDOMESTIMATE;
	read_random_theta(options, buffer);
        break;
        break;
    case 'F':
    case '_':
        options->thetaguess = FST;
        break;
        //case 'G':
    case 'O':
    case '0':
        //perhaps to come
        //      if(varvalue[0]=='G')
        //options->thetaguess = PARAMGRID;
        //else
        options->thetaguess = OWN;
        ch = skip_sspace (buffer);
        if (ch == '\0')
            return;
        if (ch == '{')
        {

            while (ch != '}')
            {
                i = 0;
                ch = skip_sspace (buffer);
                if (ch == '\0')
                    return;
                while (ch != ' ' && ch != ',' && ch != '}')
                {
                    tmp[i++] = ch;
                    ch = sgetc (buffer);
                }
                tmp[i] = '\0';
                options->thetag[options->numthetag] = atof (tmp);
                options->numthetag += 1;
                options->thetag =
                    (MYREAL *) myrealloc (options->thetag,
                                        sizeof (MYREAL) * (1 +
                                                           options->numthetag));

            }
        }
        else
        {
            i = 0;
            tmp[i++] = ch;
            while (!isspace (ch))
            {
                tmp[i++] = ch;
                ch = sgetc (buffer);
            }
            tmp[i] = '\0';
            options->thetag[options->numthetag] = atof (tmp);
            options->numthetag += 1;
            options->thetag =
                (MYREAL *) myrealloc (options->thetag,
                                    sizeof (MYREAL) * (1 + options->numthetag));
        }
        options->numpop = options->numthetag;
        if (options->numthetag == 0)
        {
            warning ("You forgot to add your guess value:\n");
            warning ("Theta=Own:{pop1,pop2, ...}\n");
            warning ("or Theta=Own:guess_pop (same value for all)\n");
        }
        break;
    default:
        usererror
        ("Failure to read start theta method, should be\ntheta=FST or theta=Own:x.x\n or theta=Own:{x.x, x.x , x.x, .....}");
    }
       myfree(otmp);
    //myfree(otempstr);
}


void
read_mig_worker (char **buffer, option_fmt * options)
{
    char ch;
    long test = 0;
    long i = 0;

    //    char input[LINESIZE];
    char varvalue[LINESIZE];

    //char *tempstr;
   char *tmp;
   //char *otempstr;
   char *otmp;

   //tempstr = (char *) mycalloc(LINESIZE,sizeof(char));
   tmp = (char *) mycalloc(LINESIZE,sizeof(char));
   //otempstr = tempstr;
   otmp = tmp;

    /* 1st example:  1.0 (n-island model)
       2nd example: {1.0} (migration matrix model, all the same start values  
       3rd example: the dashes on the diagonal are NECESSARY, {} are facultativ
       -  1.0 0.1
       1.0  -  2.0
       0.9 1.2  -
       to specify real 0.0 you need to use the custom-migration settings.
       0.0 in the table will be change to SMALLES_MIGRATION
     */

    while ((ch = sgetc (buffer)) != '=')
    {
        //      parmvar[i++] = ch;
    }
    i = 0;
    ch = sgetc (buffer);
    while (!isspace ((int) ch) && ch != ':' && ch != '{')
    {
        varvalue[i++] = ch;
        ch = sgetc (buffer);
    }
    switch (uppercase (varvalue[0]))
    {
    case 'N':
        options->migrguess = NRANDOMESTIMATE;
	read_random_mig(options, buffer);
    case 'U':
        options->migrguess = URANDOMESTIMATE;
	read_random_mig(options, buffer);
        break;
    case 'F':
    case '_':
        options->migrguess = FST;
        break;
        //    case 'G':
    case 'O':
    case '0':
        //      if(varvalue[0]=='G')
        //       options->migrguess = PARAMGRID;
        //      else
        options->migrguess = OWN;
        ch = skip_sspace (buffer);
        if (ch == '\0')
            return;
        if (ch == '{')
        {
            options->migration_model = MATRIX;
            while (ch != '}')
            {
                ch = skip_sspace (buffer);
                if ((ch == '\0') || (ch == '}'))
                    return;
                i = 0;
                while (!isspace (ch) && ch != ',' && ch != '}')
                {
                    tmp[i++] = ch;
                    ch = sgetc (buffer);
                }
                tmp[i] = '\0';
                if (strcmp (tmp, "-"))
                {
                    options->mg[options->nummg] = atof (tmp);
                    options->nummg += 1;
                    options->mg =
                        (MYREAL *) myrealloc (options->mg,
                                            sizeof (MYREAL) * (1 +
                                                               options->nummg));
                }
                else
                {
                    test++;
                }
            }
            options->numpop = test;
        }
        else
        {
            options->migration_model = ISLAND;
            i = 0;
            options->numpop = 1;
            tmp[i++] = ch;
            while (!isspace (ch))
            {
                tmp[i++] = ch;
                ch = sgetc (buffer);
            }
            options->mg[options->nummg] = atof (tmp);
            options->nummg += 1;
        }
        if (options->nummg == 0)
        {
            warning ("You forgot to add your guess value, use either:\n");
            warning ("migration=FST\n");
            warning ("migration=Own:{migration matrix, diagonal is -}\n");
            usererror
            ("migration=Own:migration_value, all matrix elements have the same value\n");
        }
        break;
    default:
        usererror ("Failure to read start migration method\n");
    }
    myfree(otmp);
    //myfree(otempstr);
}
#endif

char
skip_space (option_fmt * options)
{
    char ch = getc (options->parmfile);
    while (isspace ((int) ch) || ch == ',')
    {
        ch = getc (options->parmfile);
    }
    if (isalpha (ch))
    {
        ungetc (ch, options->parmfile);
        ch = '\0';
    }
    return ch;
}

#ifdef MPI
char
skip_sspace (char **buffer)
{
    char ch = sgetc (buffer);
    while (isspace ((int) ch) || ch == ',')
    {
        ch = sgetc (buffer);
        if (isalpha (ch))
            return ch;
    }
    return ch;
}
#endif

void
set_profile_options (option_fmt * options)
{
    switch (options->profile)
    {
    case _NONE:
        options->printprofsummary = options->printprofile = FALSE;
        break;
    case ALL:
        options->printprofsummary = options->printprofile = TRUE;
        break;
    case TABLES:
        options->printprofsummary = FALSE;
        options->printprofile = TRUE;
        break;
    case SUMMARY:
        options->printprofsummary = TRUE;
        options->printprofile = FALSE;
        break;
    }
    if (options->profilemethod == 'd')
        options->printprofsummary = FALSE;
}



/* custom-migration:<{> migration matrix and theta on
   diagonal:
   0 means not estimated,
   x means estimated, s means symmetrically
   estimated, m means all are the same, 
   c means remains constant at start value <}>
   example: 
   {* * s
    * c *
    s 0 *}
*/
void
read_custom_migration (FILE * file, option_fmt * options, char *value,
                       long customnumpop)
{
  long cc = customnumpop;
    long zz = 0, z = 0;
    char ch = '\0';
    long lc, numpop, i, j, ii;
    //    long position = 0;
    //    fpos_t thePos;
    if (cc == 0)
        cc = 1000000;
    else
        cc *= cc;

    z = 0;
    zz = 0;
    while (ch != '}' && zz < cc + options->gamma)
    {
        ch = value[z];
        switch (ch)
        {
        case '}':
        case '{':
        case ' ':
        case '\t':
            z++;
            break;
        case '\0':
        case '\n':
            z = 0;
            if (file == stdin)
                printf ("Enter the next value or list of values\n");
            FGETS (value, LINESIZE, file);
            break;
        default:
            options->custm =
                (char *) myrealloc (options->custm, sizeof (char) * (zz + 2));
            options->custm2 =
                (char *) myrealloc (options->custm2, sizeof (char) * (zz + 2));
            switch (ch)
            {
            case 'S':
                options->custm[zz++] = ch;
                break;
            case 'M':  //do we have code for this?
                options->custm[zz++] = ch;
                break;
            case 'x':
            case 'X':
                options->custm[zz++] = '*';
                break;
            case 'c':
            case 'C':
                options->custm[zz++] = 'c';
                break;
            default:
                options->custm[zz++] = tolower (ch);
                break;
            }
            z++;
        }
    }
    options->custm[zz] = '\0';
    lc = (long) strlen (options->custm);
    numpop = (long) sqrt ((MYREAL) lc);
    z = numpop;
    for (i = 0; i < numpop; i++)
    {
        for (j = 0; j < numpop; j++)
        {
            ii = i * numpop + j;
            if (i == j)
                options->custm2[i] = options->custm[ii];
            else
                options->custm2[z++] = options->custm[ii];
        }
    }
    options->custm2[z] = '\0';
    specify_migration_type (options);
    if (file != stdin)
      {
	//position = ftell (options->parmfile);
	fgetpos(options->parmfile, &thePos);
      }
    while (file != stdin && !(strstr (value, "end") || strchr (value, '=')))
    {
      // position = ftell (options->parmfile);
	fgetpos(options->parmfile, &thePos);
        FGETS (value, LINESIZE, file);
    }
    if (file != stdin)
      {
	//reset_oneline (options, position);
	reset_oneline (options, &thePos);
      }
    //printf("%li %li %li : %s\n", zz, z, lc, options->custm);
}

#ifdef MPI
void
read_custom_migration_worker (char **buffer, option_fmt * options,
                              char *value, long customnumpop)
{
  long cc = customnumpop;
    long zz = 0, z = 0;
    char ch = '\0';
    long lc, numpop, i, j, ii;
    char *position;

    if (cc == 0)
        cc = 1000000;
    else
        cc *= cc;

    z = 0;
    zz = 0;
    while (ch != '}' && zz < cc)
    {
        ch = value[z];
        switch (ch)
        {
        case '}':
        case '{':
        case ' ':
        case '\t':
            z++;
            break;
        case '\0':
        case '\n':
            z = 0;
            sgets (value, LINESIZE, buffer);
            break;
        default:
            options->custm =
                (char *) myrealloc (options->custm, sizeof (char) * (zz + 2));
            options->custm2 =
                (char *) myrealloc (options->custm2, sizeof (char) * (zz + 2));
            switch (ch)
            {
            case 'S':
                options->custm[zz++] = ch;
                break;
            case 'M':  //do we have code for this?
                options->custm[zz++] = ch;
                break;
            case 'x':
            case 'X':
                options->custm[zz++] = '*';
                break;
            default:
                options->custm[zz++] = tolower (ch);
                break;
            }
            z++;
        }
    }
    options->custm[zz] = '\0';
    lc = (long) strlen (options->custm);
    numpop = (long) sqrt ((MYREAL) lc);
    z = numpop;
    for (i = 0; i < numpop; i++)
    {
        for (j = 0; j < numpop; j++)
        {
            ii = i * numpop + j;
            if (i == j)
                options->custm2[i] = options->custm[ii];
            else
                options->custm2[z++] = options->custm[ii];
        }
    }
    options->custm2[z] = '\0';
    specify_migration_type (options);
    position = *buffer;
    while (!(strstr (value, "end") || strchr (value, '=')))
    {
        position = *buffer;
        sgets (value, LINESIZE, buffer);
    }
    *buffer = position;
    ;
}
#endif /*MPI*/
void
specify_migration_type (option_fmt * options)
{
    long len = (long) strlen (options->custm);
    long ms = 0, xs = 0, ns = 0, ss = 0, len2, i;
    char *p;
    p = options->custm;
    while (*p != '\0')
    {
        switch (*p)
        {
        case 'm':
            ms++;
            break;
        case 'x':
        case '*':
            xs++;
            break;
        case '0':
            ns++;
            break;
        case 'S':
        case 's':
            ss++;
            break;
        case 'c':
            break;
        }
        p++;
    }
    if (ms >= len)
    {
        options->migration_model = ISLAND;
        return;
    }
    if (xs >= len)
    {
        options->migration_model = MATRIX;
        return;
    }
    if (ns >= len)
    {
        usererror ("Custom migration matrix was completely set to zero?!\n");
        return;
    }
    len2 = (long) sqrt ((MYREAL) len);
    if (ms == len2 && xs == len - len2)
    {
        for (i = 0; i < len2; i++)
        {
            if (options->custm[i * len2 + i] != 'm')
            {
                options->migration_model = MATRIX;
                return;
            }
        }
        options->migration_model = MATRIX_SAMETHETA;
        return;
    }
    if (xs == len2 && ms == len - len2)
    {
        for (i = 0; i < len2; i++)
        {
            if (options->custm[i * len2 + i] != '*')
            {
                options->migration_model = MATRIX;
                return;
            }
        }
        options->migration_model = ISLAND_VARTHETA;
        return;
    }
    options->migration_model = MATRIX_ARBITRARY;
}



void
fillup_custm (long len, world_fmt * world, option_fmt * options)
{
    long i, j, ii, z;
    char *tmp;
    len = (long) strlen (options->custm);
    tmp = (char *) mycalloc (1, sizeof (char) * (world->numpop2 + 2));
    options->custm =
        (char *) myrealloc (options->custm, sizeof (char) * (world->numpop2 + 2));

    options->custm2 =
        (char *) myrealloc (options->custm2, sizeof (char) * (world->numpop2 + 2));
    strncpy (tmp, options->custm, world->numpop2 + options->gamma);
    z = world->numpop;
    for (i = 0; i < world->numpop; i++)
    {
        for (j = 0; j < world->numpop; j++)
        {
            ii = i * world->numpop + j;
            if (ii < len)
                options->custm[ii] = tmp[ii];
            else
                options->custm[ii] = '*';
            if (i == j)
                options->custm2[i] = options->custm[ii];
            else
                options->custm2[z++] = options->custm[ii];
        }
    }
    if (options->gamma)
    {
        if (len <= world->numpop2)
        {
            options->custm[world->numpop2] = '*';
            options->custm2[world->numpop2] = '*';
        }
        else
        {
            options->custm[world->numpop2] = tmp[world->numpop2];
            options->custm2[world->numpop2] = tmp[world->numpop2];
        }
        options->custm[world->numpop2 + 1] = '\0';
        options->custm2[world->numpop2 + 1] = '\0';
    }
    else
    {
        options->custm[world->numpop2] = '\0';
        options->custm2[world->numpop2] = '\0';
    }
    strcpy (world->options->custm, options->custm);
    strcpy (world->options->custm2, options->custm2);
    myfree(tmp);
}

void
print_arbitrary_migration_table (FILE * file, world_fmt * world, option_fmt *options,
                                 data_fmt * data)
{
  long i,j,z;
    char mytext[LONGLINESIZE];
    switch (world->options->migration_model)
    {
    case ISLAND:
        strcpy (mytext, "N-Island migration model");
        fprintf (file, "Migration model:\n   %-44.44s\n", mytext);
        break;
    case ISLAND_VARTHETA:
        strcpy (mytext, "N-Island migration model with variable Theta");
        fprintf (file, "Migration model:\n   %-44.44s\n", mytext);
        break;
    case MATRIX:
        strcpy (mytext, "Migration matrix model with variable Theta ");
        fprintf (file, "Migration model:\n   %-44.44s\n", mytext);
        break;
    case MATRIX_SAMETHETA:
        strcpy (mytext, "Migration matrix model with same Theta");
        fprintf (file, "Migration model:\n   %-44.44s\n", mytext);
        break;
    case MATRIX_ARBITRARY:
    default:
        strcpy (mytext, "Arbitrary migration matrix model");
        fprintf (file, "Migration model: %-44.44s\n", mytext);
        fprintf (file,
                 "[Legend: m = average (average over a group of Thetas or M]\n");
        fprintf (file,
                 "[s = symmetric M, S = symmetric 4Nm,\n 0 = zero, and not estimated,   ]\n");
        fprintf (file, "[* = free to vary, Thetas are on diagonal]\n");
        if (world->options->migration_model == MATRIX_ARBITRARY)
        {
            for (i = 0; i < data->numpop; i++)
            {
	      fprintf (file, "%10.10s     ", data->popnames[i]);
	      for (j = 0; j < data->numpop; j++)
		{
		  z = world->numpop * (options->newpops[i]-1) + (options->newpops[j]-1);
		  fprintf (file, "%c ", world->options->custm[z]);
		}
            fprintf (file, "\n");
            }
            fprintf (file, "\n\n");
        }
        break;
    }
}

void
print_distance_table (FILE * file, world_fmt * world, option_fmt * options,
                      data_fmt * data)
{
    long i, j;

    if (options->geo)
    {
        fprintf (file,
                 "   Geographic distance matrix between locations\n      ");
        for (i = 0; i < world->numpop; i++)
        {
            fprintf (file, "%-10.10s     ", data->popnames[i]);
            for (j = 0; j < world->numpop; j++)
            {
                if (i == j)
                    fprintf (file, "   -   ");
                else
                    fprintf (file, "%6.3f ", data->ogeo[i][j]);
            }
            fprintf (file, "\n      ");
        }
        fprintf (file, "\n");
    }
}

void reorder_populations(world_fmt *world, option_fmt *options, data_fmt *data)
{
  //  long newnumpop=0;
  long pop;
  long lim = data->numpop;
  long i;
  for(pop=0;pop<lim;pop++)
    {
      for(i=0; i<world->sumtips;i++)
	{
	  if(world->nodep[i]->actualpop == pop)
	    {
	      world->nodep[i]->actualpop = world->nodep[i]->pop = options->newpops[pop]-1; 
	    }
	}
    }
}

void set_localities(char **value, char **tmp, option_fmt *options)
{
      long z=0;
      options->newpops[0] = 1;
      get_next_word(value,"{}:,; ",tmp);
      while(*tmp != NULL)
	{
	  if(z >= options->newpops_numalloc)
	    {
	      options->newpops_numalloc = z+1;
	      options->newpops = (long*) myrealloc(options->newpops, options->newpops_numalloc * sizeof(long));
	    }
	  options->newpops[z] = atol(*tmp);
	  //printf("%i> population relabel [%li]=%li (input: |%s|)\n",myID,z, options->newpops[z],*tmp);
	  z++;
	  get_next_word(value,"{}:,; ",tmp);
	}
}
