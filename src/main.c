/*! \file main.c */
/*! \mainpage M I G R A T E
*  \section intro Introduction
*  Migrate is a program to estimate population-genetic parameters from genetic data such as
*  electrophoretic marker data, microsatellite data, DNA or RNA sequence data, and single 
*  nucleotide polymorphisms. Migrate uses the concepts of Maximum likelihood and Bayesian
*  inference to infer parameters based on coalescence theory. For more information related
*  to the use of the program or the theory behind check out 
*  http://popgen.scs.fsu.edu/migrate.html
*  
*  \section install Installation
*  to install the program simply do
*  \subsection step1 configure
*  \subsection step2 make
*  \subsection step3 sudo make install
*
*  \section maintainer Maintainer and Copyrights 
*   Peter Beerli,
*   beerli@fsu.edu
*
*   Copyright 1997-2002 Peter Beerli and Joseph Felsenstein, Seattle WA\n
*   Copyright 2003-2008 Peter Beerli, Tallahassee FL\n
*
*   This software is distributed free of charge for non-commercial use
*   and is copyrighted. Of course, I do not guarantee that the software
*   works and am not responsible for any damage you may cause or have.
*
*/
/* $Id: main.c 1857 2011-03-24 18:24:18Z beerli $*/

#include "migration.h"
#include "sighandler.h"
#include "heating.h"
#include "histogram.h"
#include "world.h"
#include "data.h"
#include "options.h"
#include "mcmc.h"
#include "bayes.h"
#include "broyden.h"
#include "combroyden.h"
#include "menu.h"
#include "random.h"
#include "sequence.h"
#include "tree.h"
#include "tools.h"
#include "profile.h"
#include "aic.h"
#include "lrt.h"
#include "migrate_mpi.h"
#include "migevents.h"
#include "skyline.h"
#include "reporter.h"
#ifdef PRETTY
#include "pretty.h"
#endif /*PRETTY*/
#include <stdio.h>
#ifdef BEAGLE
#include "calculator.h"
#include "mutationmodel.h"
#endif
#ifdef UEP
#include "uep.h"
#endif
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

#ifdef SNOWLEOPARD
#include <dispatch/dispatch.h>
#endif

//use this if you have nice() and need to be nice with others
#ifdef HAVE_NICE
#include <unistd.h>
#endif

#ifdef PRETTY
#ifdef WIN32
#include "pretty-win32.h"
#endif
#endif

/* Definitions for MCMCMC -- heated chains
 * definitions of the first four heated chains, they are ordered from cold to hot */
#define EARTH  universe[0]	/*!< cold chain has always a temperature=1            */
#define VENUS  universe[1]	/*!< second coldest chain                             */
#define MERKUR universe[2]	/*!< thrid coldest chain                              */
#define SUN    universe[3]	/*!< fourth coldest chain                             */

/* GLOBAL VARIABLES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int myID;	  /*!< myID=0 for single-cpu program and master in MPI program, worker nodes myID > 0 */
int myRepID;	  /*!< myID=0 for single-cpu program and master in MPI program, worker nodes myID > 0 */
int color;	  /*!< either 1 or MPI_UNDEFINED, in anticipation of more complicated MPI schemes     */
int numcpu;	  /*!< used if MPI is running otherwise == 1                                          */
int locidone;	  /*!< used in MPI workers                                                            */
long unique_id_global;
#ifdef MPI
#ifdef  USE_MYREAL_FLOAT
#ifdef WINDOWS
#define mpisizeof MPI_FLOAT ;
#else
const MPI_Datatype mpisizeof = MPI_FLOAT ;
#endif
#else
#ifdef WINDOWS
//#define mpisizeof MPI_DOUBLE ;
 MPI_Datatype mpisizeof ;
#else
const  MPI_Datatype mpisizeof = MPI_DOUBLE ;
#endif
#endif
#endif

#ifdef CAUTIOUS
boolean cautious; /*!< forces a check whether a file exists and will ask before overwriting           */
#endif
#ifdef MEMDEBUG
#include <sys/time.h>
FILE *memfile;
struct timeval memt_start, memt_finish;
double memelapsed;
long totalsize;
#endif
#ifdef MPI
filedb_fmt filedb[30];
long filenum;

MPI_Comm comm_world;		/*!< parallel environment that contains knowledge about all workers and master */
//MPI_Comm comm_workers;		/*!< parallel environment that contains knowledge about all workers */
MPI_Group worker_group;
MPI_Group world_group;
#endif
#ifdef SLOWNET
int profiledone;		/*!< used in MPI workers in profiles on slow networks    */
#endif
#ifdef PTHREADS
tpool_t heating_pool;		/*!< when compiled with PTHREADS then holds all threads */
#else
#ifndef tpool_t
#define tpool_t char
#endif
tpool_t heating_pool;
#endif

// random generator related global variables
long *seed;			/*!< contains the seed of the random number */
long *newseed;			/*!< contains the new random seed */
char *generator;		/*!< string that shows what random number generator is used */

// for bayesian output counter, initialized in bayes_init()
long *mdimfilecount;

#ifdef PRETTY
// for pretty printing
pdf_doc doc;
pdf_page page;
pdf_contents canvas;
int page_counter;
char pdf_pagetitle[LINESIZE+1];
char pdf_time[LINESIZE+1];
float page_height;
float left_margin;
float page_width;
#endif


int
setup_locus (long locus, world_fmt * world, option_fmt * options,
             data_fmt * data);

void
condense_time (world_fmt * world, long *step, long *j,
               MYREAL * accepted, long *G, long *steps, long oldsteps,
               long rep);

MYREAL sumbezier(long intervals, MYREAL x0, MYREAL y0, MYREAL x1, MYREAL y1, MYREAL x2, MYREAL y2);
void calculate_BF(world_fmt **universe, option_fmt *options);

long set_repkind (option_fmt * options);

void
heating_init (world_fmt ** universe, int usize, data_fmt * data,
              option_fmt * options);

void
heating_prepare (world_fmt ** universe, int usize,
                 option_fmt * options, data_fmt * data, long rep);

void heating_prepare2 (world_fmt ** universe, int usize);

long
replicate_number (option_fmt * options, long chain, char type,
                  long replicate);

void
combine_loci_standard (char type, option_fmt * options, world_fmt * world,
                       long *Gmax);
void
combine_loci_mpi (char type, option_fmt * options, world_fmt * world,
                  long *Gmax);

void set_penalizer (long chain, long chains, char type,
                    world_fmt ** universe);

void
run_sampler (option_fmt * options, data_fmt * data,
             world_fmt ** universe, int usize, long *outfilepos, long *Gmax);

void run_replicate (long locus,
                    long replicate,
                    world_fmt ** universe,
                    option_fmt * options,
                    data_fmt * data,
                    tpool_t * localheating_pool,
                    int usize, long *treefilepos, long *Gmax);


void
run_locus (world_fmt ** universe, int usize, option_fmt * options,
           data_fmt * data, tpool_t * localheating_pool, long maxreplicate,
           long locus, long *treefilepos, long *Gmax);

void
run_loci (world_fmt ** universe, int usize, option_fmt * options,
          data_fmt * data, tpool_t * localheating_pool, long maxreplicate,
          long *treefilepos, long *Gmax);

void run_one_update (world_fmt * world);

void run_updates (world_fmt ** universe,
                  int usize, option_fmt * options,
                  tpool_t * localheating_pool, long inc, long increment,
                  long step, long steps);

void heated_swap (world_fmt ** universe, worldoption_fmt * options);

void
change_chaintype (long locus, char *type, world_fmt * world,
                  long *increment, long *oldsteps, long *chains,
                  option_fmt * options);

void
prepare_next_chain (world_fmt ** universe, worldoption_fmt * options,
                    char type, long chain, long *chains,
                    long *pluschain, long locus, long replicate);

void print_bayesfactor(FILE *file, world_fmt **universe, option_fmt * options);
#ifdef MPI
//void fix_bayesfactor(world_fmt *world, option_fmt * options);
void      print_marginal_like(float *temp, long *z, world_fmt * world);
#else /*not MPI*/
void      print_marginal_like(char *temp, long *c, world_fmt * world);
#endif

void print_burnin_stop(FILE *file, world_fmt **universe, option_fmt * options);
void
print_heating_progress (world_fmt ** universe,
                        worldoption_fmt * options, long stepinc);

long analyze_olddata (world_fmt * world, option_fmt * options, data_fmt * data, long *outfilepos);
boolean analyze_oldbayesdata(world_fmt *world, option_fmt *options, data_fmt *data, long *outfilepos);
void print_theta0(FILE *file, world_fmt *world, long maxreplicate);
void profile_tables (option_fmt * options, world_fmt * world, long *gmaxptr);

void finish_mac (option_fmt * options, data_fmt * data);

int
setup_locus (long locus, world_fmt * world, option_fmt * options,
             data_fmt * data);

long set_repkind (option_fmt * options);

void heating_prepare2 (world_fmt ** universe, int usize);

long replicate_number (option_fmt * options, long chain, char type,
                       long replicate);

void print_heating_progress2 (FILE * file, worldoption_fmt * options,
                              world_fmt ** universe);


void get_bayeshist (world_fmt * world, option_fmt * options);
void get_treedata (world_fmt * world, option_fmt * options);
void get_mighistdata (world_fmt * world, option_fmt * options);

void change_longsum_times (world_fmt * world);

boolean  check_parmfile(long argcount, char **arguments, char *parmfilename);
boolean set_usemenu(boolean usemenu, boolean fromparmfile);
void check_bayes_options(option_fmt *options);
void reset_bayesmdimfile(world_fmt *world, option_fmt *options);

extern void unset_penalizer_function (boolean inprofiles);

///
/// specifies the classes of override the menu as a parameter to the program
enum override_enum 
  {
    OVERRIDE_NO, OVERRIDE_MENU, OVERRIDE_NOMENU
  };



///
/// the program migrate calculates migration rates and population sizes
/// from genetic data, it allows for a wide variety of data types and many options
/// \callgraph
int
main (int argc, char **argv)
{
  char      type       = 's';
  int       usize      = 1;
  int       usemenu    = OVERRIDE_NO;
  long      locus;
  long      Gmax       = 0;
  long      outfilepos = 0;
  data_fmt  *data;
  boolean   restarted_bayes_bool=FALSE;
#ifdef INTEGRATEDLIKE
  long      maxreplicate;
#endif
  option_fmt *options;
  world_fmt **universe;
#ifdef MPI
  int       server;
#ifdef WINDOWS
 mpisizeof = MPI_DOUBLE ;
#endif
#endif
#ifdef HAVE_NICE
  nice (10); //nice value arbitrarily set to 10, 5 is still nice, 0 is standard
#endif
  // puts(argv[0]);
  
  //---------------------------------------------------------------------------------------------
  // windows specific code for pretty printing
#ifdef PRETTY
#ifndef MPI
#ifdef WIN32
    set_haru_handler();
#endif
#endif
#endif

    //---------------------------------------------------------------------------------------------
    // MPI initialisation
#ifdef MPI
    // parallel version
    filenum = 0;			//setting up the filename database so that worker node can write back to master
    MPI_Init (&argc, &argv);
    comm_world = MPI_COMM_WORLD;
    MPI_Comm_size (comm_world, &numcpu);
    MPI_Comm_rank (comm_world, &myID);
    MPI_Comm_group (comm_world, &world_group);
    server = MASTER;		//server ID
    MPI_Group_excl (world_group, 1, &server, &worker_group);

    locidone = 0;

#ifdef SLOWNET
    // slow network parallel version -- this will turn into the standard version
    profiledone = 0;
    which_calc_like (SINGLELOCUS);
#endif

#else /*MPI*/
    //scalar version of migrate
    myID = MASTER;
    numcpu = 1;
#endif /*MPI*/

    // debug code for memory problems
#ifdef MEMDEBUG
    totalsize = 0 ;
    memfile = fopen("memoryfile","w+");
    gettimeofday(&memt_start, NULL);
#endif

    //---------------------------------------------------------------------------------------------
    // initialization of main container and random number parts
    seed = (long *) mymalloc (sizeof (long) * 3);
    newseed = (long *) mymalloc (sizeof (long) * 3);
    generator = (char *) mycalloc (1,sizeof(char) * 80);
    universe = (world_fmt **) mycalloc (HEATED_CHAIN_NUM, sizeof (world_fmt *));
    
    // set flag to decide whether we use caution when writing files
#ifdef CAUTIOUS
    cautious = FALSE;
#endif
    // try to catch and beautify some error messages
    signalhandling (ON);
    unique_id_global=0;
    // create main data structures
    create_data (&data);
    create_world (&(EARTH), 1L);
    create_options (&options);
    
    // parmfile and menu
    init_options (options);

    //---------------------------------------------------------------------------------------------
    // master initializations, this works for both parallel and single-cpu version
    if (myID == MASTER)
      {
	usemenu = check_parmfile(argc,argv,options->parmfilename);
	read_options_master (options);
	options->menu = set_usemenu(usemenu, options->menu);
	print_menu_title (stdout, options);
	get_menu (options, EARTH, data);
	check_bayes_options(options);
#ifdef MPI        
	MYMPIBARRIER (comm_world);
        broadcast_options_master (options, data);
#endif
      }
    else
      {
#ifdef MPI
	MYMPIBARRIER (comm_world);
	broadcast_options_worker (options);
#endif
	options->menu = FALSE;
      }

    //---------------------------------------------------------------------------------------------
    // data initialization and data reading phase
    usize = (options->heating ? ( MAX (1, options->heated_chains)) : 1 );
    universe = (world_fmt **) myrealloc (universe, usize * sizeof (world_fmt *));
    // opens files on all nodes    
    init_files (EARTH, data, options);
    
    if (options->writelog && myID == MASTER)
      {
        print_menu_title (options->logfile, options);
      }
    
    EARTH->repkind = SINGLECHAIN;
    //fprintf(stderr,"myID=%i myRepID=%i\n",myID, myRepID);    
    unset_penalizer_function (FALSE);	// install penalizer for far jumps in MCMC
    
    //---------------------------------------------------------------------------------------------
    // sampling phase
    if (!options->readsum) // all go here except when reading old runs
      {
        run_sampler (options, data, universe, usize, &outfilepos, &Gmax);
      }
    else
      {
	if(options->bayes_infer)
	  {
	    // reanalyze old or broken bayesrun [not ready yet]
	    restarted_bayes_bool = analyze_oldbayesdata(EARTH, options, data, &outfilepos);
	  }
	else
	  {
	    //reanalyze sumfile
	    Gmax = analyze_olddata (EARTH, options, data, &outfilepos);
	  }
      }
    unset_penalizer_function (TRUE);
    
    //---------------------------------------------------------------------------------------------
    // combining phase for multiple loci and replicates
#ifdef MPI
    //with MPI:the workers stay here much longer than the Master
    combine_loci_mpi (type, options, EARTH, &Gmax);
#else
    combine_loci_standard (type, options, EARTH, &Gmax);
#endif
    // bayes inference 
    if (options->bayes_infer)
      {
	// if bayes intermediate data recording is ON then reset the file
	// for reading for printing and combining, else use the material still in RAM
	if(!restarted_bayes_bool)
	  reset_bayesmdimfile(EARTH, options);
      }
    else
      {
        // gather the treedata(1) into sumfile or from worker nodes
        get_treedata (EARTH, options);
      }

    //---------------------------------------------------------------------------------------------
    // printing main results
    if (myID == MASTER)
      {
	if (options->bayes_infer)
	  {
	    // printe bayes inference main table
	    bayes_stat (EARTH);
	  }
	else
	  {
	    // print ML table
	    print_list (universe, options, data);
	    print_alpha_curve (EARTH, EARTH->atl, &Gmax);
	    if(options->printcov)
	      print_cov(EARTH, EARTH->numpop, EARTH->loci, EARTH->cov);
	    if (EARTH->options->lratio->counter > 0)
	      print_lratio_test (EARTH, &Gmax);
	    if (options->aic)
	      akaike_information (EARTH, &Gmax);
	  }
	fflush (EARTH->outfile);
	//---------------------------------------------------------------------------------------------
	// printing additional results
#ifdef MPI
	get_mighistdata (EARTH, options);
#endif
	if(options->skyline && options->bayes_infer)
	  {
	    print_expected_values(EARTH, options);       
	  }
	if(options->mighist || options->skyline)
	  {
	    print_event_values(EARTH);
#ifdef PRETTY
	    pdf_histogram_legend();
#endif
	  }              
#ifdef UEP
	// print UEP probabilities
	if (options->uep)
	  analyze_uep (EARTH);
#endif
        // profile tables        unset_penalizer_function (TRUE);	//now we calculate profile and shut down
	// the penalizing of far jumps in the maximizer function
      } // end of printing main tables and histograms
    
#ifdef SLOWNET
    which_calc_like (PROFILE);
    if (myID == MASTER)
      {
	// release the workers from the mpi_maximize_worker() function.
	// the workers will stop working on locus-likelihoods and advance to profiles
	mpi_send_stop (EARTH);
      }
    if (!options->bayes_infer)
      {
	profile_tables (options, EARTH, &Gmax);
      }
#else  /* SLOWNET*/
    if (myID == MASTER)
      {
#ifdef INTEGRATEDLIKE
	if(options->integrated_like)
	  {
	    maxreplicate = (options->replicate
			    && options->replicatenum >
			    0) ? options->replicatenum : 1;
	    print_theta0(stdout, EARTH, maxreplicate);
	    
	    if (EARTH->options->lratio->counter > 0)
	      print_lratio_test (EARTH, &Gmax);
	    if (options->aic)
	      akaike_information (EARTH, &Gmax);
	    
	    profile_tables (options, EARTH, &Gmax);
	  }
#else
	if (!options->bayes_infer)
	  {
	    profile_tables (options, EARTH, &Gmax);
	  }
#endif /*INTEGRATEDLIKE*/
      }
#endif /* SLOWNET */
    
#ifdef MPI
    if (myID == MASTER)
      {
#endif
	//#ifdef BFDEBUG	
	if(EARTH->options->bayes_infer)
	//#else
	//if(EARTH->options->bayes_infer  && options->datatype!='g')
	  //#endif
	  {
	// print marginal likelihoods
#ifdef MPI
	    //    fix_bayesfactor(EARTH,options);
#endif
	    print_bayesfactor(EARTH->outfile, universe,options);
	// printing of MCMC run characteristics
	    fprintf(EARTH->outfile,"MCMC run characteristics\n");
	    fprintf(EARTH->outfile,"========================\n\n");
	    bayes_print_accept(EARTH->outfile,EARTH);
#ifdef PRETTY
	    pdf_bayes_print_accept(EARTH);
#endif
	    fprintf(EARTH->outfile,"Autocorrelation and Effective sample size\n");
	    fprintf(EARTH->outfile,"-------------------------------------------------------------------\n\n");
	    print_bayes_ess(EARTH->outfile,EARTH,EARTH->numpop2 + EARTH->bayes->mu * EARTH->loci + 1,2,EARTH->auto_archive, EARTH->ess_archive);
#ifdef PRETTY
	    pdf_bayes_print_ess(EARTH);
#endif
	  }
	
	if(strchr("ae",EARTH->options->burnin_autostop))
	  {
	    print_burnin_stop(EARTH->outfile, universe, options);
	  }
	
	if(EARTH->options->progress)
	  {
#ifdef BFDEBUG
	    if(EARTH->options->bayes_infer)
#else
	    if(EARTH->options->bayes_infer && options->datatype!='g')
#endif
	      {
		fprintf(stdout,"\nMCMC run characteristics: Autocorrelation and Effective sample size\n");
		fprintf(stdout,"-------------------------------------------------------------------\n\n");
		print_bayes_ess(stdout,EARTH,EARTH->numpop2+ EARTH->bayes->mu * EARTH->loci + 1,2,
				     EARTH->auto_archive, EARTH->ess_archive);
		
		print_bayesfactor(stdout, universe,options);
		
	      }
	  }
	if(options->treeprint)
	  {
	    if(EARTH->options->treeinmemory)
	      {
		for(locus=0;locus<EARTH->loci; locus++)
		  {
		    fprintf(EARTH->treefile,"%s", EARTH->treespace[locus]);
		  }
	      }
#ifdef NEXUSTREE
	    FPRINTF(EARTH->treefile,"\nend;\n");
#endif
	  }

	// if adaptive heating print a table with the average temperatures
	//
	if(options->heating && options->adaptiveheat!=NOTADAPTIVE)
	  {
	    fprintf(EARTH->outfile,"\n\n\nAverage temperatures during the run using %s\n",(options->adaptiveheat!=NOTADAPTIVE)
		    ==STANDARD ? "standard adaptive heating scheme" : "bounded adaptive heating scheme" );
	    fprintf(EARTH->outfile,"===========================================================================\n\n");
	    fprintf(EARTH->outfile,"Chain Temperature\n");
	    // locus means indicator for chain
	    for(locus = 0 ; locus < options->heated_chains; locus++)
	      {
		fprintf(EARTH->outfile,"%5li %10.5f\n",locus+1,universe[locus]->averageheat);
	      }
	    fprintf(EARTH->outfile,"Adaptive heating often fails, if the average temperatures are very close together\n");
	    fprintf(EARTH->outfile,"try to rerun using static heating! If you want to compare models using marginal\n");
	    fprintf(EARTH->outfile,"likelihoods then you MUST use static heating\n");
	    pdf_print_averageheat(universe,options);
	  }
	//
	// print warnings into the PDF and the outfile
	//
	print_stored_warnings(EARTH);
	pdf_print_stored_warnings(EARTH);
	//
	//
	//
	pdf_print_end_time(&page_height);
	// write out to PDF file
	pdf_write_file(options);
	print_finish (EARTH, outfilepos);
	// closing all files 
#ifdef MAC
	finish_mac (options, data);
#endif
	
	exit_files (EARTH, data, options);
#ifdef MPI
#    ifndef SLOWNET
	mpi_send_stop (EARTH); //stop workers
	// printf("%i> at barrier in slownet-main before finalize\n", myID);        
#    endif
      }
    MYMPIBARRIER (comm_world);
    MPI_Finalize ();
#endif /*MPI*/
    //	myfree(seed);
    //myfree(newseed);
    //myfree(generator);
    //free_universe(universe, usize, options);
    //destroy_data(data);
    //destroy_options(options);
    // for debugging with MallocDebug on macosx
#ifdef MACMALLOCDEBUG
    while(TRUE)
      {
	sleep(1);
      } //end for debugging
#endif
    return 0;
}				/* main end */

///
/// Calculates the maximum likelihood estimates from the trees that were gathered for each locus
/// \callgraph
void
combine_loci_standard (char type, option_fmt * options, world_fmt * world,
                       long *Gmax)
{
    long kind = MULTILOCUS;
    if (options->readsum)
    {
        if (world->loci == 1)
            kind = SINGLELOCUS;
    }
#ifndef INTEGRATEDLIKE
    if ((world->loci - world->skipped > 1) || (options->readsum))
#else
    if ((world->loci - world->skipped > 1) || (options->readsum))
#endif
    {
        if (options->gamma)
            world->param0[world->numpop2] = world->options->alphavalue;
        world->repkind = set_repkind (options);
        decide_plot (world->options, world->options->lchains,
                     world->options->lchains, 'l');
#ifdef LONGSUM
        change_longsum_times (world);	//multilocus
#endif /*LONGSUM*/
#ifdef INTEGRATEDLIKE
	    (void) estimateParameter (options->replicatenum, *Gmax, world, options,\
				      world->cov[world->loci], 0 /* chain */ ,\
				      type, kind, world->repkind);
#else
	if(!options->bayes_infer)
	  {
	    (void) estimateParameter (options->replicatenum, *Gmax, world, options,
				      world->cov[world->loci], 0 /* chain */ ,
				      type, kind, world->repkind);
	  }
#endif
    }
}

///
/// Calculates the maximum likelihood estimates from the trees that were gathered for each locus
/// \callgraph
void
combine_loci_mpi (char type, option_fmt * options, world_fmt * world,
                  long *Gmax)
{
#ifdef MPI
    long kind = MULTILOCUS;
    if (options->readsum)
    {
        if (world->loci == 1)
            kind = SINGLELOCUS;
    }
    if ((world->loci - world->skipped > 0) || (options->readsum))
    {
        if (options->gamma)
            world->param0[world->numpop2] = world->options->alphavalue;
        world->repkind = set_repkind (options);
        if (myID == MASTER)
        {
#ifdef INTEGRATEDLIKE
	      mpi_gmax_master (world, Gmax);
	      mpi_startparam_master (world);
	      decide_plot (world->options, world->options->lchains,
			   world->options->lchains, 'l');
#else
	  if(!options->bayes_infer)
	    {
	      mpi_gmax_master (world, Gmax);
	      mpi_startparam_master (world);
	      decide_plot (world->options, world->options->lchains,
			   world->options->lchains, 'l');
	    }
#endif
#ifdef LONGSUM
            change_longsum_times (world);	//multilocus
#endif	   /*LONGSUM*/
            // broadcast Gmax to all workers, the worker pendant is in mpi_gmax_worker
#ifndef INTEGRATEDLIKE
	    if (!options->bayes_infer)
	      {
		MYMPIBCAST (Gmax, 1, MPI_LONG, MASTER, comm_world);
		(void) estimateParameter (options->replicatenum, *Gmax, world,
					  options, world->cov[world->loci],
					  0 /* chain */ ,
					  type, kind, world->repkind);
	      }
#else
		MYMPIBCAST (Gmax, 1, MPI_LONG, MASTER, comm_world);
		(void) estimateParameter (options->replicatenum, *Gmax, world,
					  options, world->cov[world->loci],
					  0 /* chain */ ,
					  type, kind, world->repkind);
#endif
	}
        else
        {
#ifndef INTEGRATEDLIKE
            //            printf("%i> before gmax worker\n", myID);
	  if(!options->bayes_infer)
	    {
	      mpi_gmax_worker (world);
	      mpi_startparam_worker (world);
	    }
            //            printf("%i> start worker postlike calculations\n", myID);
            mpi_maximize_worker (world, MULTILOCUS, options->replicatenum);	// receive broadcast for Gmax
#else
	    mpi_gmax_worker (world);
	    mpi_startparam_worker (world);
            mpi_maximize_worker (world, MULTILOCUS, options->replicatenum);	// receive broadcast for Gmax
#endif
        }
    }
    //   else
    //   {
    //       if (myID != MASTER)
    //       {
    //           maxreplicate = (options->replicate
    //                           && options->replicatenum >
    //                           0) ? options->replicatenum : 1;
    //           fprintf(stdout,"about to pack result buffer and send too");
    //           mpi_results_worker(MIGMPI_SUMFILE, world, maxreplicate, pack_result_buffer);
    //           //mpi_maximize_worker(world, options->replicatenum);
    //       }
    //   }
#endif
}

///
/// calculate the number populations to analyze and resize the newpop array if necessary
long  calculate_newpop_numpop(option_fmt *options, data_fmt *data)
{
  long pop;
  long i;
  long newnumpop=0;
  long numpop = data->numpop;
  //  char *strsep(char **, const char *);
#ifdef WIN32
    // dev studio does not understand the #warning preprocessor directive, it uses pragma but I am not ready for that yet
    //#warning "Windows version is hardcoded to not have more than 5000 populations"
    long temp[5000];
#else
  long temp[data->numpop];
#endif
  //use order of populations {1,2,3,4,5,....} , do not start with zero!!
  //reset to options->newpops={1,1,2,1,2,3,4,...}
  if(options->newpops_numalloc < numpop)
    {
      options->newpops = (long *) myrealloc(options->newpops,sizeof(long)*numpop);
      for(i=options->newpops_numalloc;i<numpop;i++)
	{
	  options->newpops[i] = i+1;
	}
      options->newpops_numalloc = numpop;
    }
  memcpy(temp,options->newpops,sizeof(long)*options->newpops_numalloc);
  qsort((void *) temp, options->newpops_numalloc, sizeof(long), longcmp);
  pop = temp[0];
  newnumpop = 1;
  for(i=1;i<numpop;i++)
    {
      if(temp[i] != pop)
	{
	  newnumpop++;
	  pop = temp[i];
	}
    }
  options->newpops_numpop = newnumpop;
  return newnumpop;
}


///
/// the slice sampling stick size is allocated here for the options, this is 
/// done in three different places for: standard, genealogy, MPI
void alloc_sticksize(option_fmt *options, data_fmt *data)
{
  if( options->slice_sticksizes==NULL)
    {
      options->slice_sticksizes = (MYREAL *) mycalloc(data->numpop * data->numpop + 1, sizeof(MYREAL));
    }
  else
    {
      options->slice_sticksizes = (MYREAL *) myrealloc(options->slice_sticksizes, (data->numpop * data->numpop + 1)* sizeof(MYREAL));
    }
}

///
/// runs the MCMC sampling 
void
run_sampler (option_fmt * options, data_fmt * data, world_fmt ** universe,
             int usize, long *outfilepos, long *Gmax)
{
  MYREAL var;
  long i;
  long maxreplicate;
  long treefilepos;
#ifdef MPI
  MPI_Request *irequests;	// contains requests generated in MPI_Isend
  MPI_Status *istatus;		// conatins stata from MPI_Wait_any
  long minnodes;
  long *twolongs;		// message sent to workers contains locus and replicate number
  
  if (myID < data->loci + 1)
    color = 1;
  else
    color = MPI_UNDEFINED;
  
  twolongs = (long *) mycalloc (TWO, sizeof (long));
#endif
  getseed(options);
#ifdef MPI
  if (myID == MASTER)
    {
#endif
      get_data (data->infile, data, options, EARTH);
#ifdef MPI
      //        printf("%i > finished data reading\n", myID);
      if (numcpu == 1)
        {
	  error
            ("This program was compiled to use a parallel computer\n and you tried to run it on only a single node.\nThis will not work because it uses a \n\"single_master-many_worker\" architecture \nand needs at least TWO nodes\n");
        }
      broadcast_data_master (data, options);
    } /*end of myID==MASTER loop for MPI*/
  else
    {
      broadcast_data_worker (data, options);
    }
#endif
  set_plot (options);
  EARTH->cold = TRUE; // this is the cold chain when we use heating, the hotter chains have FALSE
  // sticksizes are filled in init_world but initialized here so make sure that 
  // the heated chains do not produce a memory leak by mutiply initializing the sticksizes
  alloc_sticksize(options,data);
  // filling of options->slice_sticksizes
  // deferred after world->bayes initialization in init_world() 
  calculate_newpop_numpop(options,data);

  init_world (EARTH, data, options);
  //  set_meanmu(EARTH,options);
  if(EARTH->cold)
    single_chain_var (EARTH, 0, &var, NULL, NULL);
  
#ifdef NEXUSTREE
  if (EARTH->options->treeprint == LASTCHAIN)
    {
      FPRINTF(EARTH->treefile,"#NEXUS\n\nbegin trees;");
    }
#endif
  
#ifdef PTHREADS
  tpool_init (&heating_pool, usize, usize, 0);
#endif
  
  if (options->heating)
    {
      //first in universe is the cold chain[EARTH]
      heating_init (universe, usize, data, options);
    }

  //#ifdef MPI
  //  if(myID!=MASTER)
  //  {
  //    if(options->printfst || options->thetaguess==FST || options->migrguess==FST)
  //	{
  //	  calc_simple_param (EARTH, data);
  //	}
  // }
  //#endif
  /* report to screen */
  if (myID == MASTER)
    {
      print_menu_options (EARTH, options, data);
      if (options->progress)
	print_data_summary (stdout, EARTH, options, data);
      /* print to outfile */
#ifdef PRETTY
      pdf_master_init(EARTH, options, data);
#endif 
      *outfilepos = print_title (EARTH, options);
      print_options (EARTH->outfile, EARTH, options, data);
      print_data_summary (EARTH->outfile, EARTH, options, data);
      print_data (EARTH, options, data);
      print_spectra (EARTH, options, data);
      
      if(options->bayes_infer && myID == MASTER)
	{
	  if(options->has_bayesmdimfile)
	    print_bayes_mdimfileheader(EARTH->bayesmdimfile,
				       options->bayesmdiminterval, EARTH, data);
	}
    }
  if (options->lchains < 2 && options->replicatenum == 0)
    {
      options->replicate = 0;
    }
  maxreplicate = (options->replicate
		  && options->replicatenum > 0) ? options->replicatenum : 1;
#ifdef MPI
    minnodes = MAX (0, numcpu - data->loci - 1);
    irequests = (MPI_Request *) mycalloc (minnodes + 1, sizeof (MPI_Request));
    istatus = (MPI_Status *) mycalloc (minnodes + 1, sizeof (MPI_Status));
    if (myID == MASTER)
      {
	mpi_runloci_master (data->loci, EARTH->who, EARTH, options->readsum, options->menu);
#ifndef MPIREPLICANT
        //sent = mpi_send_stop_mcmc_worker (numcpu, data->loci, &comm_world, irequests,istatus, myID);
#endif
    }
    else
    {
        mpi_runloci_worker (universe, usize, options, data,
                            &heating_pool, maxreplicate, &treefilepos, Gmax);
#ifdef MPIREPLICANT
        // the first worker will always be here waiting 
        //if (myID == FIRSTWORKER)
        //{
	//  MYMPIRECV (twolongs, TWO, MPI_LONG, MASTER, 0, comm_world, &status);
#ifdef MPI_DEBUG
	//    fprintf(stdout,"%i> FIRSTWORKER kill off now non-loci-worker\n",myID);
#endif
	//    sent = mpi_send_stop_mcmc_replicateworker(numcpu,data->loci);
	    //                mpi_send_stop_mcmc_worker (numcpu, data->loci, &comm_workers,
	    //                             irequests, istatus, myID);
	// }
#endif /*MPIREPLICANT*/
    }
    myfree(istatus);
    myfree(irequests);
#else /*MPI*/
    run_loci (universe, usize, options, data,
              &heating_pool, maxreplicate, &treefilepos, Gmax);

#endif /*MPI*/

#ifdef BEAGLE
    beagle_stop(universe, usize);
#endif

#ifdef PTHREADS
    tpool_destroy (heating_pool, 1);
#endif
#ifdef MPI
    
    if (myID != MASTER)
    {
#endif
      if (options->heating)
	{
	  for (i = 0; i < options->heated_chains; i++)
	    {
	      //free_tree(universe[i]->root, universe[i]);
	      free_timevector (universe[i]->treetimes);
	    }
	}
      else
	{
	  //free_tree(EARTH->root, EARTH);
	  free_timevector (EARTH->treetimes);
	}
#ifdef MPI
    }
    myfree(twolongs);
#endif
}

/// generate samples of genealogies for all loci
void
run_loci (world_fmt ** universe, int usize, option_fmt * options,
          data_fmt * data, tpool_t * localheating_pool, long maxreplicate,
          long *treefilepos, long *Gmax)
{
    long locus;
    //  long            i;
    for (locus = 0; locus < data->loci; locus++)
    {
      if(!data->skiploci[locus])
	{
	  run_locus (universe, usize, options, data,
		     localheating_pool, maxreplicate, locus, treefilepos, Gmax);
	}
      free_datapart(data,options,locus);
    }
}


/// save all genealogy summaries
void
get_treedata (world_fmt * world, option_fmt * options)
{
#ifdef MPI
    long maxreplicate = (options->replicate
                         && options->replicatenum >
                         0) ? options->replicatenum : 1;
#endif
    
    if (myID == MASTER)
    {
        if (options->writesum)
        {
#ifdef MPI
            //get all sumfile data from the workers using unpack_sumfile_buffer()
            mpi_results_master (MIGMPI_SUMFILE, world, maxreplicate,
                                unpack_sumfile_buffer);
#endif
            
            write_savesum (world);
        }
#ifdef SLOWNET
        else
        {
            //      printf("%i> about to receive the sumfile data\n", myID);
            mpi_results_master (MIGMPI_SUMFILE, world, maxreplicate,
                                unpack_sumfile_buffer);
            //      printf("%i> all sumfile data received and unpacked\n", myID);
        }
#endif
        
    }
}

///
/// save all genealogy summaries
void
get_bayeshist (world_fmt * world, option_fmt * options)
{
#ifdef MPI
  //long i;
    long maxreplicate = (options->replicate
                         && options->replicatenum >
                         0) ? options->replicatenum : 1;
    
    if (myID == MASTER)
    {
      //for(i=0;i<world->loci;i++)
      //	  printf("%i> before moving stuff  %f %f\n",myID, world->hm[i],world->am[i]);
        mpi_results_master (MIGMPI_BAYESHIST, world, maxreplicate,
                            unpack_bayes_buffer);
	//for(i=0;i<world->loci;i++)
	//  printf("%i> after moving stuff  %f %f\n",myID, world->hm[i],world->am[i]);
    }
#endif
}

void 	recalc_skyline_values(world_fmt *world, option_fmt * options, long maxreplicate)
{
  long i, j, locus;
  mighistloci_fmt *aa;
  long * eventnum;
  const float invmax = 1./maxreplicate;
  if(world->options->mighist && world->options->skyline)
    {
      for(locus=0; locus < world->loci; locus++)
	{
	  aa = &world->mighistloci[locus];
	  eventnum = world->mighistloci[locus].eventbinnum;
	  for(j=0; j< world->numpop2; j++)
	    {
	      for (i = 0; i < eventnum[j]; i++)
		{
		  aa->eventbins[j][i][0] *= invmax;
		  aa->eventbins[j][i][1] *= invmax;
		  aa->eventbins[j][i][2] *= invmax;
		  aa->eventbins[j][i][3] *= invmax;
		  aa->eventbins[j][i][4] *= invmax;
		  aa->eventbins[j][i][5] *= invmax;
		}
	    }
	}
    }
}


/// get migration event time data
/// and skyline data when present
void
get_mighistdata (world_fmt * world, option_fmt * options)
{
#ifdef MPI
    long maxreplicate = (options->replicate
                         && options->replicatenum >
                         0) ? options->replicatenum : 1;
    
    if (myID == MASTER)
    {
        if (options->mighist)
	  {
            //get all mighist data from the workers using unpack_mighist_buffer()
            mpi_results_master (MIGMPI_MIGHIST, world, maxreplicate,
                                unpack_mighist_buffer);
	    if(options->skyline)
	      {
		//get all mighist data from the workers using unpack_mighist_buffer()
		fprintf(stdout,"%i> before unpack skyline\n",myID);
		mpi_results_master (MIGMPI_SKYLINE, world, maxreplicate,
				    unpack_skyline_buffer);
	      }
	  }
	if(/*world->options->treeprint ==BEST &&*/ world->options->treeinmemory==TRUE)
	  {
            mpi_results_master (MIGMPI_TREESPACE, world, maxreplicate,
                                unpack_treespace_buffer);
	  }
	recalc_skyline_values(world, options, maxreplicate);
    }
#endif
    
}

/// swap the tree pointer between chains with different temperatures
void
heated_swap (world_fmt ** universe, worldoption_fmt * options)
{
    long sm = 0, mv = 0, vw = 0;
    long r;
    //debug for heating traverse issues
    if(options->heatedswap_off || (options->heating_count++ % options->heating_interval) != 0)
      return;
    switch (r = RANDINT (0, options->heated_chains - 2))
    {
        case 2:
            sm = chance_swap_tree (SUN, MERKUR);
            MERKUR->swapped += sm;
            break;
        case 1:
            mv = chance_swap_tree (MERKUR, VENUS);
            VENUS->swapped += mv;
            break;
        case 0:
            vw = chance_swap_tree (VENUS, EARTH);
            EARTH->swapped += vw;
            break;
        default:
            universe[r]->swapped += chance_swap_tree (universe[r + 1], universe[r]);
    }
}


/// change the type of the chain form short to long and reset the treelist
void
change_chaintype (long locus, char *type, world_fmt * world, long *increment,
                  long *oldsteps, long *chains, option_fmt * options)
{
    if (*type == 's')
    {
        *type = 'l';
        create_treetimelist (world, &(world->treetimes), locus);
        set_bounds (increment, oldsteps, chains, options, *type);
        world->increment = *increment;
    }
}

/// prepare the next chain
void
prepare_next_chain (world_fmt ** universe, worldoption_fmt * options,
                    char type, long chain, long *chains, long *pluschain,
                    long locus, long replicate)
{
    long i;
    EARTH->likelihood[0] = EARTH->likelihood[EARTH->G];
    if (options->heating)
    {
        for (i = 1; i < options->heated_chains; i++)
            universe[i]->likelihood[0] = universe[i]->likelihood[universe[i]->G];
    }
    EARTH->treetimes[0].copies = 0;
    if (type == 'l')
    {
        if (options->replicate && options->replicatenum > 0)
            EARTH->chainlikes[locus][replicate] = EARTH->param_like;
        else
            EARTH->chainlikes[locus][chain] = EARTH->param_like;
    }
    EARTH->start = FALSE;
    if (type == 'l')
    {
        if (chain < *chains + *pluschain)
        {
            if (((EARTH->param_like > options->lcepsilon)
		 || (options->gelman && EARTH->convergence->gelmanmaxRall > GELMAN_MYSTIC_VALUE)) 
		&& !(options->replicate && options->replicatenum == 0))
            {
                (*chains)++;
                (*pluschain)--;
            }
        }
    }
}

///
/// resetting of accept and accept_freq for heated chains
void reset_heated_accept(world_fmt **universe, long unum)
{
  long i;
  for(i=0; i < unum; i++)
    {
      universe[i]->swapped = 0;
      universe[i]->accept = 0L;
      universe[i]->accept_freq = 0.;
    }
}


/// print progress about heated chains
void
print_heating_progress (world_fmt ** universe, worldoption_fmt * options,
                        long stepinc)
{
  MYREAL fstepinc = (MYREAL) stepinc;
    if (options->heating)
    {
#ifdef DEBUG_MPI
      printf("%i> stepinc=%li =========\n",myID, stepinc);
#endif
        universe[0]->accept_freq /= fstepinc;
        universe[1]->accept_freq /= fstepinc;
        universe[2]->accept_freq /= fstepinc;
        universe[3]->accept_freq /= fstepinc;
        if (options->progress)
            print_heating_progress2 (stdout, options, universe);
        if (options->writelog)
            print_heating_progress2 (options->logfile, options, universe);
    }
}

///
/// sub-function of print_heating, in MPI mode, the string printing will minimize data transfer.
void
print_heating_progress2 (FILE * file, worldoption_fmt * options,
                         world_fmt ** universe)
{
  char *plog;
  long plogsize = 0L;
  plog = (char *) mycalloc(LONGLINESIZE,sizeof(char));
#ifdef MPI
  plogsize = sprintf(plog, "[%3i] Swapping between %li temperatures.\n", myID, options->heated_chains);
#else
  plogsize = sprintf(plog, "           Swapping between %li temperatures.\n", options->heated_chains);
#endif
  if (options->heated_chains > 4)
    plogsize += sprintf (plog + plogsize, "               Only the 4 coldest chains are shown\n");
  plogsize += sprintf (plog + plogsize,"           Temperature | Accepted |  Swaps between temperatures\n");
  if (options->heated_chains > 4)
    plogsize += sprintf (plog + plogsize, "           %10.4g  |   %4.2f   |    %5li||\n",
			 universe[3]->averageheat, universe[3]->accept_freq,
			 universe[3]->swapped);
  else
    plogsize += sprintf (plog + plogsize, "           %10.4g  |   %4.2f   |          |\n",
			 universe[3]->averageheat, universe[3]->accept_freq);
  plogsize += sprintf (plog + plogsize, "           %10.4f  |   %4.2f   |    %5li ||\n",
		       universe[2]->averageheat, universe[2]->accept_freq,
		       universe[2]->swapped);
  plogsize += sprintf (plog + plogsize,"           %10.4f  |   %4.2f   |    %5li  ||\n",
		       universe[1]->averageheat, universe[1]->accept_freq,
		       universe[1]->swapped);
  //xcode plogsize += 
    sprintf (plog + plogsize,"           %10.4f  |   %4.2f   |    %5li   |\n",
		       universe[0]->averageheat, universe[0]->accept_freq,
		       universe[0]->swapped);
  FPRINTF(file,"%s",plog);
  myfree(plog);
}


/// analyze old sumfile data
long
analyze_olddata (world_fmt * world,
                 option_fmt * options, data_fmt * data, long *outfilepos)
{
#ifndef MPI
    long locus;
#else
    long listofthings[6];
#endif
    long Gmax = 0;
    long pop;
    options->heating = FALSE;
    options->gelman = FALSE;
    if (myID == MASTER)
    {
        unset_penalizer_function (TRUE);
        Gmax = read_savesum (world, options, data);
        
        data->popnames = (char **) mymalloc (sizeof (char *) * world->numpop);
        for (pop = 0; pop < world->numpop; pop++)
        {
            data->popnames[pop] =
            (char *) mycalloc (1, sizeof (char) * LINESIZE);
            MYSNPRINTF (data->popnames[pop], LINESIZE, "Pop_%li ", pop + 1);
        }
        
#ifdef MPI
#ifdef WINDOWS	
        MYMPIBCAST (data->geo, world->numpop2, MPI_DOUBLE, MASTER, comm_world );
        MYMPIBCAST (data->lgeo, world->numpop2, MPI_DOUBLE, MASTER, comm_world );
#else
        MYMPIBCAST (data->geo, world->numpop2, mpisizeof, MASTER, comm_world );
        MYMPIBCAST (data->lgeo, world->numpop2, mpisizeof, MASTER, comm_world );
#endif
    }
    else
    {
        MYMPIBCAST (listofthings, 6, MPI_LONG, MASTER, comm_world );
        //      printf("listofthings received\n");
        //    fflush(stdout);
        world->loci = listofthings[0];
        world->numpop = listofthings[1];
        world->numpop2 = listofthings[2];
        options->replicate = (boolean) listofthings[3];
        options->replicatenum = listofthings[4];
        Gmax = listofthings[5];
        data->numpop = world->numpop;
        data->geo = (MYREAL *) mycalloc (1, sizeof (MYREAL) * world->numpop2);
        data->lgeo = (MYREAL *) mycalloc (1, sizeof (MYREAL) * world->numpop2);
#ifdef WINDOWS
        MYMPIBCAST (data->geo, world->numpop2, MPI_DOUBLE, MASTER, comm_world);
        MYMPIBCAST (data->lgeo, world->numpop2, MPI_DOUBLE, MASTER, comm_world);
#else
        MYMPIBCAST (data->geo, world->numpop2, mpisizeof, MASTER, comm_world);
        MYMPIBCAST (data->lgeo, world->numpop2, mpisizeof, MASTER, comm_world);
#endif
        options->muloci = data->loci = world->loci;
        init_world (world, data, options);
#endif
        
    }
#ifdef MPI
#ifdef SLOWNET
    mpi_broadcast_results (world, world->loci,
                           pack_sumfile_buffer, unpack_sumfile_buffer);
#endif
#endif
    
    synchronize_param (world, options);
    if (options->plot)
    {
        options->plotnow = TRUE;
        world->options->plotnow = TRUE;
    }
    world->locus = 0;
    if (myID == MASTER)
    {
        *outfilepos = print_title (world, options);
        print_options (stdout, world, options, data);
        print_options (world->outfile, world, options, data);
#ifdef PRETTY
	pdf_master_init(world, options, data);
#endif
        
#ifdef MPI
	MYMPIBARRIER(comm_world);
	printf("%i> before mpi_runlocimaster\n",myID);fflush(stdout);
	mpi_runloci_master (world->loci, world->who, world, options->readsum, options->menu);
	
#else
        
        for (locus = 0; locus < world->loci; locus++)
        {
            if (options->replicate && options->replicatenum > 0)
            {
                world->locus = locus;
                world->repkind = MULTIPLERUN;
#ifdef LONGSUM
                change_longsum_times (world);	//multi run replicates
#endif	       /*LONGSUM*/
		if (!options->bayes_infer)
		  {
		    (void) estimateParameter (options->replicatenum, Gmax, world,
					      options, world->cov[locus], 1, 'l',
					      SINGLELOCUS, world->repkind);
		  }
	    }
        }
#endif
        
    }
#ifdef MPI
    else
    {
      MYMPIBARRIER(comm_world);
      printf("%i> before assignloci_worker\n",myID);fflush(stdout);
      assignloci_worker (world, options, &Gmax);
      // some workers might still wait in assignloci_worker(), will be stopped through the master
    }
#endif
    return Gmax;
} /*analyze_olddata()*/


boolean analyze_oldbayesdata(world_fmt *world, option_fmt *options, data_fmt *data, long *outfilepos)
{
#ifdef MPI
  long pop;
#endif
  long locus;
  // read data from bayesfile
  world->cold=TRUE;
  // reads some minimal information form the header of the bayesallfile
  // this only works with files written with migrate 2.5+
#ifdef MPI
  if(myID==MASTER)
    {
      warning("does not work correctly yet");
      read_from_bayesmdim_minimal_info(world->bayesmdimfile, world, options, data);
      read_geofile (data, options, world->numpop);
      alloc_sticksize(options,data);
      options->newpops_numalloc = world->numpop;
      options->newpops = (long*) mycalloc(options->newpops_numalloc, sizeof(long));
      for(pop=0;pop<world->numpop;pop++)
	options->newpops[pop]=pop+1;
      init_world (world, data, options);
      read_bayes_fromfile(world->bayesmdimfile, world, options);
    }
#else
  read_from_bayesmdim_minimal_info(world->bayesmdimfile, world, options, data);
  read_geofile (data, options, world->numpop);
  alloc_sticksize(options,data);
  options->newpops_numalloc = world->numpop;
  options->newpops = (long*) mycalloc(options->newpops_numalloc, sizeof(long));
  long pop;
  for(pop=0;pop<world->numpop;pop++)
    options->newpops[pop]=pop+1;
  init_world (world, data, options);
  read_bayes_fromfile(world->bayesmdimfile, world, options);
#endif
  //set_meanmu(world,options);

#ifdef PRETTY
  pdf_master_init(world, options, data);
#endif
  for(locus=0;locus<world->loci;locus++)
    {
      if(world->data->skiploci[locus])
	continue;
      if(world->bayes->histogram[locus].results == NULL)
	{
	  world->data->skiploci[locus] = TRUE;
	  continue;
	}
      calc_hpd_credibility(world->bayes, locus, world->numpop2, world->numpop2 + world->bayes->mu);
    }
  return TRUE;
}

/// generate profile likelihood tables
void
profile_tables (option_fmt * options, world_fmt * world, long *gmaxptr)
{
    long len;
#ifndef SLOWNET
#ifdef PRETTY
    long location = 0;
#endif
    long i;
#endif
    
    if (options->profile != _NONE)
    {
      //        world->buffer =
      //  (char *) myrealloc (world->buffer, allocbufsize * sizeof (char));
        if (options->printprofsummary)
            allocate_profile_percentiles (world);
#ifdef HAVE_STRFTIME        
        time (&world->starttime);
#endif        
        len = world->numpop2 + (world->options->gamma ? 1 : 0);
#ifdef LONGSUM        
        len += world->numpop * 3;
#endif /*LONGSUM*/
#ifdef SLOWNET
        if (!options->readsum)
	  {
	    mpi_broadcast_results (world, world->loci,
                                   pack_sumfile_buffer, unpack_sumfile_buffer);
	  }
        mpi_broadcast_results (world,
                               world->loci + ((world->loci == 1) ? 0 : 1),
                               pack_result_buffer, unpack_result_buffer);
        //fprintf(stdout,"%i >>>>>>>>> wait for all at barrier in profiles_tables [main.c:1053]",myID);
        MYMPIBARRIER (comm_world);
        if (myID == MASTER)
        {
            print_profile_title (world);
#ifdef PRETTY
	    pdf_print_profile_title(world);
#endif	    
            mpi_profiles_master (world, len, world->profilewho);
        }
        else
            mpi_profiles_worker (world, gmaxptr);
#else /*SLOWNET -not*/
	if(myID==MASTER && (options->profile == TABLES || options->profile == ALL))
	  {
	    print_profile_title (world);
#ifdef PRETTY
	    pdf_print_profile_title(world);
#endif	         
	  }   
        for (i = 0; i < len; i++)
	  {
            boolean success =  print_profile_likelihood_driver (i, world, gmaxptr);
	    if(options->profile == TABLES || options->profile == ALL)
	      {
		LARGEFPRINTF (world->outfile, world->allocbufsize+2, "%s\n\n", world->buffer);
		//ascii_print_profile_table(options->profilemethod, world);
#ifdef PRETTY
		location = strlen(world->buffer);
		if(success)
		  pdf_print_profile_table(55., &page_height, options->profilemethod, world->buffer + location + 1, world);
#endif 
	      }
            world->buffer[0] = '\0';
	  }
#endif /*SLOWNET end of not*/
        if (myID == MASTER && options->printprofsummary)
        {
            world->allocbufsize = print_profile_percentile (world);
            // print profile summary
            LARGEFPRINTF (world->outfile, world->allocbufsize, "%s\n\n", world->buffer);
#ifdef PRETTY
	    pdf_print_profile_percentile (world);
#endif
            world->buffer[0] = '\0';
            destroy_profile_percentiles (world);
        }
    }
    //myfree(world->buffer);
}

/// finish up mac files \todo remove this function
void
finish_mac (option_fmt * options, data_fmt * data)
{
#ifdef MAC
#if 0
    fixmacfile (options->outfilename);
    if (options->plotmethod == PLOTALL)
        fixmacfile (options->mathfilename);
    if (options->treeprint)
        fixmacfile (options->treefilename);
    if (options->parmfile)
        fixmacfile ("parmfile");
    if (data->sumfile)
        fixmacfile (options->sumfilename);
#endif
#endif
}

/// \brief setup a locus
/// 
/// Set up a the structures for a locus. Run-Parameters are copied into the 
/// major structure WORLD from options and data. a first genealogy is build
/// adjusted to the starting parameters (changing times of nodes) and a first
/// conditional likelihood is calculated for all nodes, also the tree-parallele
/// time structure is generated.
/// Several run controllers are checked: replicate, datatype
/// \retvalue 0 { failed to setup structure for locus because there is no data }
/// \retvalue 1 { succeeded to setup locus-structure } 
int
setup_locus (long locus, world_fmt * world, option_fmt * options,
             data_fmt * data)
{
  world->locus = locus;
  set_param (world, data, options, locus);
  if (world->options->progress)
    print_menu_locus (stdout, world, locus);
  if (world->options->writelog)
    print_menu_locus (world->options->logfile, world, locus);
  world->start = TRUE;
  buildtree (world, options, data, locus);
#ifdef BEAGLE
  finish_mutationmodel(world, data, options);
  init_beagle(world,locus);
  set_beagle_instances(world,TRUE); // sets two instances for beagle
  //  set_beagle_instances(world,FALSE); // so that wee can easily MCMC
  fill_beagle_instances(world);
#endif
  /*  if (world->replicate == 0)
      {
      if (strchr (SNPTYPES, options->datatype))
      {
      // unlinked snps -- does this work at all
      if (options->datatype == 'u')
      {
      world->data->seq->endsite *= (data->seq->addon + 1);
      }
      if (options->datatype == 'n')
      {
      world->data->seq->endsite += data->seq->addon;
      }
      }
      }
  */

  //    if (!options->replicate
  //    || (options->replicate && options->replicatenum == 0)
  //    || (options->replicate && world->replicate == options->replicatenum))
  //{
  if (data->skiploci[locus])	/* skip loci with no tips */
    return 0;
  //}
  create_treetimelist (world, &world->treetimes, locus);
  fix_times (world, options);
#ifdef DEBUG
  //  printf(" fixed times \n");
#endif
  first_smooth (world, locus);
  world->likelihood[0] = treelikelihood (world);
  world->allikemax = world->likelihood[0]; //-MYREAL_MAX;	/* for best tree option */
  world->treelen = 0.0;
  calc_treelength (world->root->next->back, &world->treelen);

#ifdef UEP
  if (options->uep)
    {
        world->treelen = 0.0;
        if (options->ueprate > 0.0)
	  calc_treelength (world->root->next->back, &world->treelen);
        //world->ueplikelihood = ueplikelihood(world);
        //world->likelihood[0] += world->ueplikelihood;
    }
#endif
    return 1;
}

/// condense a single genealogy into a sufficient statistic for the maximization phase
void
condense_time (world_fmt * world, long *step, long *j, MYREAL * accepted,
               long *G, long *steps, long oldsteps, long rep)
{
  static long n = 0;
#ifdef INTEGRATEDLIKE
  long i;
  long z;
  long nn = world->numpop2 + world->bayes->mu * world->loci;
  MYREAL x;
  MYREAL delta;
#endif
  world->rep = rep;
  if(world->options->bayes_infer)
    {
      if (world->in_last_chain)
	{
	  if (*step == 0)
	    return;

	  // syncs treeupdate and paramupdate
	  if(world->bayesaccept == -1)
	    {
	      world->bayes->oldval = probg_treetimes(world);
	    }
	  bayes_save (world, *step * world->options->lincr);
	  store_events(world, world->treetimes, world->numpop, rep);
	  *accepted += *j;
#ifdef INTEGRATEDLIKE
	  n += 1;
	  //	  fprintf(stdout,"\naverage MP-values: %li ", n);
	  for(i=0; i < nn; i++)
	    {
	      if(world->bayes->map[i][1] == INVALID)
		continue;
	      else
		{
		  z  = world->bayes->map[i][1];
		}
	      x = world->param0[z];
	      delta = x - world->atl[rep][world->locus].param0[z];
	      world->atl[rep][world->locus].param0[z] += delta/n;
	      //  fprintf(stdout,"%f ",world->atl[rep][world->locus].param0[z]);
	    }
	  //	  fprintf(stdout,"\n");
#else 
	  return;
#endif
	}
      else
	{
	  warning("Bayesian analysis: do not run multiple long chains, run 1 long chain! [nothing will be wrong but it is a waste of time]");
	  return;
	}
    }
  if (*step == 0)
    {
      copy_time (world, world->treetimes, FIRSTSTEP, 0, world->numpop, rep);
      *accepted += *j;
      //for INTEGRATEDLIKE
      n = 0;
    }
  else
    {
      copy_time (world, world->treetimes, *G, *G + (long) (*j > 0),
		 world->numpop, rep);
      if (*step < *steps)
        {
	  //handled in advance_world * G += (long) (*j > 0);
	  *accepted += *j;
        }
    }
  if (*step >= oldsteps - 1 && world->options->movingsteps)
    {
      if (((MYREAL) (*G + 1)) < world->options->acceptfreq * oldsteps)
        {
	  (*steps)++;
        }
    }
}

void print_theta0(FILE *file, world_fmt *world, long maxreplicate)
{
  long i;
  long z;
  long locus;
  long rep;
  long nn = world->numpop2;
  fprintf(file,"Locus Replicate Parameters\n");
  fprintf(file,"-----------------------------------------------------------\n");
  for(locus=0; locus < world->loci; locus++)
    {
      for(rep = 0; rep < maxreplicate; rep++)
	{
	  fprintf(file,"%5li  %5li ",locus, rep+1);
	  for(i=0; i < nn; i++)
	    {
	      if(world->bayes->map[i][1] == INVALID)
		continue;
	      else
		{
		  z  = world->bayes->map[i][1];
		}
	      fprintf(file,"%f ", world->atl[rep][world->locus].param0[z]);
	    }
	  fprintf(file,"\n");
	}
    }
  fprintf(file,"\n\n");
}

/// set type of replication scheme
long
set_repkind (option_fmt * options)
{
    if (options->replicate)
    {
        if (options->replicatenum == 0)
            return MULTIPLECHAIN;
        else
           return MULTIPLERUN;
    }
    return SINGLECHAIN;
}

/// intialize heating scheme
void
heating_init (world_fmt ** universe, int usize, data_fmt * data,
              option_fmt * options)
{
    long chain;
    char tmp[20];
    for (chain = 1; chain < usize; ++chain)
    {
        MYSNPRINTF (tmp, 20, "%10.5f", options->heat[chain]);
        create_world (&universe[chain], data->loci);
	universe[chain]->cold = FALSE;
        init_world (universe[chain], data, options);
	universe[chain]->heatid = chain;
    }
}

/// prepare for heating scheme: Step I
void
heating_prepare (world_fmt ** universe, int usize,
                 option_fmt * options, data_fmt * data, long rep)
{
    long chain;
    EARTH->heat = options->heat[0];
    for (chain = 1; chain < usize; ++chain)
      {
	buildtree (universe[chain], options, data, EARTH->locus);//DEBUG
        klone (EARTH, universe[chain], options, data, options->heat[chain]);
      }
}

/// prepare for heating scheme: Step I [this function failed in the parallel run, unused for now, replaced by above]
void
heating_prepare_old (world_fmt ** universe, int usize,
                     option_fmt * options, data_fmt * data, long rep)
{
    long chain;
    EARTH->heat = options->heat[0];
    if (rep == 0)
    {
        for (chain = 1; chain < usize; ++chain)
            klone (EARTH, universe[chain], options, data, options->heat[chain]);
    }
    else
    {
        for (chain = 1; chain < usize; ++chain)
            klone_part (EARTH, universe[chain], options, data,
                        options->heat[chain]);
    }
}

/// prepare for heating scheme: Step II
void
heating_prepare2 (world_fmt ** universe, int usize)
{
    long chain;
    universe[0]->averageheat = 1.0;
    for (chain = 1; chain < usize; ++chain)
    {
        universe[chain]->G = 0;
        universe[chain]->averageheat = universe[chain]->options->adaptiveheat!=NOTADAPTIVE ?
            0.0 : 1. / universe[chain]->heat;
    }
    for (chain = 1; chain < usize; ++chain)
    {
        clone_polish (EARTH, universe[chain]);
    }
}

/// \brief generates replicate number
///
/// generates replicate number
/// \rtval 0 no replicate
/// \rtval rep number of rpelicates
long
replicate_number (option_fmt * options, long chain, char type, long replicate)
{
    long rep = 0;
    if (options->replicate && options->replicatenum > 0)
        rep = replicate;
    else
    {
        if (type == 'l' && options->replicate)
            rep = chain;
        else
            rep = 0;
    }
    return rep;
}

/// generate genealogies for a single locus
/// \callgraph
void
run_locus (world_fmt ** universe, int usize, option_fmt * options,
           data_fmt * data, tpool_t * localheating_pool, long maxreplicate,
           long locus, long *treefilepos, long *Gmax)
{
  long i;
  long convergence_len = universe[0]->numpop2 + 3 * universe[0]->numpop;
  long replicate;
  // DEBUG Hapmap stuff
  if(EARTH->data->skiploci[locus])
    return;
  if(options->bayes_infer)
    convergence_len = universe[0]->numpop2 + 1;
  memset(EARTH->convergence->chain_means,0,maxreplicate * convergence_len * sizeof(MYREAL));
  memset(EARTH->convergence->chain_s,0,maxreplicate * convergence_len * sizeof(MYREAL));
  memset(EARTH->convergence->gelmanmeanmaxR,0,maxreplicate * maxreplicate * sizeof(MYREAL));
  for (replicate = 0; replicate < maxreplicate; replicate++)
    {
      run_replicate (locus, replicate, universe, options, data,
		     localheating_pool, usize, treefilepos, Gmax);
      if(EARTH->cold && EARTH->options->replicatenum > 0)
	{
	  if(options->gelman)
	    {
	      chain_means(&EARTH->convergence->chain_means[replicate * convergence_len], EARTH);
	      calc_chain_s(EARTH->convergence->chain_s, EARTH->convergence->chain_means, EARTH, 
			   replicate);
	      if (replicate > 0)
		{
		  convergence_check_bayes(EARTH, maxreplicate);
		  convergence_progress(stdout, EARTH);
		}
	    }
	}
    }    
#ifdef UEP
  if (options->uep)
    show_uep_store (EARTH);
#endif
  if (options->replicate && options->replicatenum > 0)
    {
      EARTH->repkind = MULTIPLERUN;
      if (options->bayes_infer)
        {
	  if(!options->has_bayesmdimfile)
	    calculate_credibility_interval (EARTH, locus);
	  //	  return;
        }
#ifdef LONGSUM
      change_longsum_times (EARTH);
#endif /*LONGSUM*/
      // over multiple replicates if present or single locus
#ifdef INTEGRATEDLIKE
	  (void) estimateParameter (options->replicatenum, *Gmax, EARTH,
				    options, EARTH->cov[locus],
				    options->lchains, /*type */ 'l',
				    SINGLELOCUS, EARTH->repkind);
#else
      if (!options->bayes_infer)
	{
	  (void) estimateParameter (options->replicatenum, *Gmax, EARTH,
				    options, EARTH->cov[locus],
				    options->lchains, /*type */ 'l',
				    SINGLELOCUS, EARTH->repkind);
	}
#endif
    }
  else
    {
      if (options->bayes_infer)
        {
	  if(!options->has_bayesmdimfile)
	    {
	      adjust_bayes_bins(EARTH,locus);
	      calculate_credibility_interval (EARTH, locus);
	    }
	}
    }
  if (options->heating)
    {
      for (i = 0; i < options->heated_chains; i++)
        {
	  free_tree (universe[i]->root, universe[i]);
	  //   free_timevector(universe[i]->treetimes);
	  myfree(universe[i]->nodep);
        }
    }
  else
    {
      free_tree (EARTH->root, EARTH);
      myfree(EARTH->nodep);
      //  free_timevector(EARTH->treetimes);
    }
  if(options->bayes_infer)
    bayes_reset (universe[0]); 
}

/// \brief updates the tree and records acceptance
///
/// updates the tree and records acceptance
/// when in bayesian mode change between changing the tree and the parameters
void
run_one_update (world_fmt * world)
{
  long count = 0;
  long np = world->numpop2;
  //fprintf(stderr,"%i> locus=%li heat=%f\n",myID, world->locus, world->heat);
  //traverse_check(crawlback (world->root->next));
  if (world->options->bayes_infer)
    {
      if(world->bayes->mu)
	{
	  np += 1; 
	}
      world->bayesaccept = bayes_update (world);
      if (world->bayesaccept == -1)	//either do tree updates or parameter updates
        {
	  world->accept += (count = tree_update (world, world->G));
	  world->bayes->accept[np] += count;	//permanent recorder for the tree acceptance, 
	  // for some obscure reason the world->accept gets lost, most likely every cycle;
	  world->bayes->trials[np] += 1;	// counter for trees tried 
        }
      world->param_like = world->bayes->oldval;
    }
  else
    {
      world->accept += tree_update (world, world->G);
    }
}

/// \brief updates all trees, controls updates and heating
///
/// controls updates and temperatures (threaded or unthreaded)
/// updates the tree and/or parameters
/// \callgraph
void
run_updates (world_fmt ** universe,
             int usize,
             option_fmt * options,
             tpool_t * localheating_pool,
             long inc, long increment, long step, long steps)
{
#ifndef  PTHREADS
#ifndef SNOWLEOPARD
    long ii;
#endif
#endif
    if (options->heating)
    {
#ifdef PTHREADS			/*using threads and running on an SMP machine */
        fill_tpool (*localheating_pool, universe, usize);
        tpool_synchronize (*localheating_pool, usize);
#else /*heating but not using threads */
#ifdef SNOWLEOPARD
        dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
        dispatch_apply(options->heated_chains, queue,
                       ^(unsigned long ii) {
#ifdef FUNNYHEAT
			 if(ii!=0)
			   {
			     reset_weight(universe[ii]);
			   }
#endif
                         run_one_update (universe[ii]);
                       }
                       );
#else
        for (ii = 0; ii < options->heated_chains; ii++)
        {
#ifdef FUNNYHEAT
	  if(ii!=0)
	    {
	      reset_weight(universe[ii]);
	    }
#endif
	  run_one_update (universe[ii]);
        }
#endif /*SNOWLEOPARD*/
#endif /* end of not-using threads */
	if(options->bayes_infer && RANDUM() > options->updateratio)
	  {
#ifdef UEP
	    if (options->uep && EARTH->in_last_chain)
	      update_uepanc (EARTH);
#endif
	    return;
	  }
	else
	  {
	    heated_swap (universe, EARTH->options);
	    switch (options->adaptiveheat)
	      {
	      case STANDARD:
		adjust_temperatures (universe, options->heated_chains,
				     inc /*rement */  + step * increment,
				     steps * increment);
		break;
	      case BOUNDED:
		adjust_temperatures_bounded (universe, options->heated_chains,
				     inc /*rement */  + step * increment,
				     steps * increment);
		break;
	      case NOTADAPTIVE:
	      default:
		break;
	      }
	  }
    }
    else
      {		/* no heating */
	run_one_update (EARTH);
      }
#ifdef UEP
    if (options->uep && EARTH->in_last_chain)
      update_uepanc (EARTH);
#endif
}

/// generates trees over the interval increments
void
run_increments (world_fmt ** universe,
                long usize,
                option_fmt * options,
                tpool_t * localheating_pool,
                long increment, long step, long steps)
{
  static long treefilepos; /* write position in the treefile */
    long i;
    for (i = 0; i < increment; i++)
    {
        EARTH->actualinc = i;
        run_updates (universe, usize, options,
                     localheating_pool, i, increment, step, steps);
        
        if (EARTH->likelihood[EARTH->G] > EARTH->maxdatallike)
            EARTH->maxdatallike = EARTH->likelihood[EARTH->G];
    }
    if (EARTH->options->treeprint != _NONE)
      print_tree (EARTH, EARTH->G, &treefilepos);
}

/// \brief run a chain
/// run a chain
void
run_steps (world_fmt ** universe,
           long usize,
           option_fmt * options,
           tpool_t * localheating_pool, long increment, long steps, long rep)
{
    long step=0;
    long oldsteps = steps;
    long ii;
    MYREAL var=0.0;
    if(EARTH->cold && EARTH->options->bayes_infer)
      single_chain_var (NULL, step, &var, NULL, NULL);
#ifdef MPI
    // does not work? check_memory_limit();
#endif
    for (step = 0; step < steps + 1; step++)
    {
#ifdef MPI
      // does not work?? check_memory_limit();
#endif
        EARTH->increment = increment;
        run_increments (universe, usize, options, localheating_pool,
                        increment, step, steps);
#ifdef UEP
        store_uep (EARTH);
#endif
	//	printf("step=%li inc=%li locus=%li\n",step,increment,EARTH->locus);

        condense_time (EARTH, &step, &EARTH->accept,
                       &EARTH->accept_freq, &EARTH->G, &steps, oldsteps, rep);
	if(EARTH->options->bayes_infer)
	  {
	    calculate_BF(universe,options);
	    if(EARTH->cold)
	      {
		single_chain_var (EARTH, step, &var, EARTH->autocorrelation, EARTH->effective_sample);
	      }
	  }
        if (step > 0)
        {
            advance_clone_like (EARTH, EARTH->accept, &EARTH->G);
            EARTH->accept = 0L;
        }
        // if we do heating
        if (EARTH->options->heating)
        {
            for (ii = 1; ii < EARTH->options->heated_chains; ++ii)
            {
                universe[ii]->accept_freq += universe[ii]->accept;
                advance_clone_like (universe[ii],
                                    universe[ii]->accept, &(universe[ii])->G);
                universe[ii]->accept = 0L;
            }
        }
    }
}

/// run all chains for short and then for long
void
run_chains (world_fmt ** universe,
            long usize,
            option_fmt * options,
            tpool_t * localheating_pool,
            long replicate,
            long chains,
            char type,
            long increment, long oldsteps, long *treefilepos, long *Gmax)
{
  long i;
  long steps;
  long chain;
  long rep;
  const long locus = EARTH->locus;
  long kind;
  long pluschain = 0;
  const MYREAL treeupdateratio = (options->bayes_infer ? options->updateratio : 1.0);

  for (chain = 0;
       chain < chains || (type == 'l' && chain >= chains
			  && EARTH->param_like > options->lcepsilon); chain++)
    {
      if (options->heating)
	MERKUR->swapped = VENUS->swapped = EARTH->swapped = 0;
      rep = replicate_number (options, chain, type, replicate);
      EARTH->rep = rep;
      precalc_world (EARTH);
      polish_world (EARTH);
      burnin_chain (EARTH);
      EARTH->G = 0;
      EARTH->accept_freq = 0.;
      EARTH->accept = 0;
      steps = oldsteps;
      EARTH->maxdatallike = EARTH->likelihood[0];
      if (options->heating)
	{
	  heating_prepare2 (universe, usize);

	  for(i=1; i<EARTH->options->heated_chains; i++)
	    {
	      burnin_chain (universe[i]);
	    }
	  reset_heated_accept(universe,EARTH->options->heated_chains);
	}
      set_penalizer (chain, chains, type, universe);
      // run steps, run increments, run updaters 
      run_steps (universe, usize, options, localheating_pool,
		 increment, steps, rep);
      
      decide_plot (EARTH->options, chain, chains, type);
      // prepare for parameter estimation, precompute and copy
      memcpy (EARTH->param00, EARTH->param0,
	      sizeof (MYREAL) * EARTH->numpop2);
#ifndef INTEGRATELIKE
	if(!options->bayes_infer)
	  {
	    memcpy (EARTH->atl[rep][locus].param0, EARTH->param00,
		    sizeof (MYREAL) * EARTH->numpop2);
	    log_param0 (EARTH->atl[rep][locus].param0,
			EARTH->atl[rep][locus].lparam0, EARTH->numpop2);
	    copies2lcopies (&EARTH->atl[rep][locus]);
	  }
#else
	//	    memcpy (EARTH->atl[rep][locus].param0, EARTH->param00,
	//	    sizeof (MYREAL) * EARTH->numpop2);
	log_param0 (EARTH->atl[rep][locus].param0,
		    EARTH->atl[rep][locus].lparam0, EARTH->numpop2);
	copies2lcopies (&EARTH->atl[rep][locus]);
#endif
        // start parameter estimation for 
        // single chain, single replicate
        EARTH->repkind = SINGLECHAIN;
        kind = SINGLELOCUS;
        *Gmax = (EARTH->G > *Gmax) ? EARTH->G + 1 : (*Gmax);
#ifdef LONGSUM
        change_longsum_times (EARTH);
#endif /*LONGSUM*/
	if (!options->bayes_infer)
	  {
	    (void) estimateParameter (rep, *Gmax, EARTH, options,
				      EARTH->cov[locus], chain, type, kind,
				      EARTH->repkind);
	  }
        if (EARTH->options->heating)
	  {
	    //            print_heating_progress (universe, EARTH->options, EARTH->options->heating_interval * increment * steps * treeupdateratio);
            print_heating_progress (universe, EARTH->options, (long) increment * steps * treeupdateratio);
	    reset_heated_accept(universe,EARTH->options->heated_chains);
	  }
        else
	  {
            print_progress (EARTH->options, EARTH, rep,
                            (long) increment * steps * treeupdateratio, (long) EARTH->accept_freq);
	    EARTH->accept_freq = 0. ;
	    EARTH->accept = 0L ;
	  }
        // cleanup and prepare next
        prepare_next_chain (universe, EARTH->options, type, chain,
                            &chains, &pluschain, locus, replicate);
    }				//end chains
}

/// set penalizing distance for jumps in the maximizer on or off
void
set_penalizer (long chain, long chains, char type, world_fmt ** universe)
{
    if (type == 'l')
    {
        if (chain >= chains - 1)
        {
            EARTH->in_last_chain = TRUE;
            unset_penalizer_function (TRUE);
        }
    }
    else
    {
        EARTH->in_last_chain = FALSE;
        unset_penalizer_function (FALSE);
    }
}


/// \brief run a replicate
///
/// run a replicate
/// \callgraph
void
run_replicate (long locus,
               long replicate,
               world_fmt ** universe,
               option_fmt * options,
               data_fmt * data,
               tpool_t * localheating_pool,
               int usize, long *treefilepos, long *Gmax)
{
  //  long i,j;
    long chains = 0;
    char type = 's';
    //long pluschain;
    long runs;
    long increment = 1;
    //long steps;
    long oldsteps = 0;
    /* loop over independent replicates----------------------- */
    /* but make sure that we get single chain estimates, too   */
    EARTH->locus = locus;
    EARTH->repkind = SINGLECHAIN;
    EARTH->replicate = replicate;
    if (setup_locus (locus, EARTH, options, data) == 0)
        return;
    if (options->heating)
        heating_prepare (universe, usize, options, data, replicate);
    type = 's';
    runs = 1;
    // rep = 0;
    //xcode pluschain = options->pluschain;
    //see obscure options
    /* short and long chains ----------------------------- */
    set_bounds (&increment, &oldsteps, &chains, options, type);
    //xcode steps = oldsteps;
    EARTH->increment = increment;
    while (runs-- >= 0)
    {
        if (myID == MASTER && type == 's')
        {
            print_menu_chain (type, FIRSTCHAIN, oldsteps, EARTH,
                              options, replicate);
            if (EARTH->options->treeprint == ALL)
                print_tree (EARTH, 0, treefilepos);
        }
        EARTH->chains = chains;
        /* loop over chains--------------------------- */
        run_chains (universe, usize, options, localheating_pool, replicate, chains, type, increment, oldsteps, treefilepos, Gmax);	//PB 020203
        change_chaintype (locus, &type, EARTH, &increment, &oldsteps,
                          &chains, options);

	/* debug
	if(EARTH->cold)
	  {
	    for(j=0; j < EARTH->numpop2; j++)
	      {
		for(i=0;i < (EARTH->mighistloci[EARTH->locus].eventbinnum[j]); i++)
		  {
		    fprintf(stdout,"(%f,%f) ",EARTH->mighistloci[EARTH->locus].eventbins[j][i][0],
			    EARTH->mighistloci[EARTH->locus].eventbins[j][i][1]);
		  }
		fprintf(stdout,"\n%i>######## %li ######\n",myID, EARTH->mighistloci[EARTH->locus].eventbinnum[j]);
	      }
	  }
	*/
        /* evaluate multiple long chain estimates */
        if (runs < 0 && options->replicate && options->replicatenum == 0)
        {
            EARTH->repkind = MULTIPLECHAIN;
#ifdef LONGSUM
            change_longsum_times (EARTH);
#endif	   /*LONGSUM*/
	if (!options->bayes_infer)
	  {
            (void) estimateParameter (replicate, *Gmax, EARTH, options,
                                      EARTH->cov[locus], chains, type,
                                      SINGLELOCUS, EARTH->repkind);
	  }
        }
    }				//end runs
    //#ifndef MPI
    // collect correlation information of a locus and stick it into the archive,
    // the autocorrelation is averaged but the ess is summed up assuming that ll loci are independent
    if(EARTH->cold && EARTH->options->bayes_infer)
      {
	  collect_ess_values(EARTH);
      }
    //#endif
}


void print_burnin_stop(FILE *file, world_fmt **universe, option_fmt * options)
{
  world_fmt *world = universe[0];
  long z;
  long maxreplicate = (options->replicate
		       && options->replicatenum >
		       0) ? options->replicatenum : 1;
  fprintf(file,"\n\n\nStop of burn-in phase due to convergence\n");
  fprintf(file,"----------------------------------------\n");
  switch(options->burnin_autostop)
    {
    case 'a':
      fprintf(file,"[Stopping criteria was: Variance ratio between the last two groups of 1000 steps < %f]\n\n",world->varheat); break;
    case 'e':
      fprintf(file,"[Stopping criteria was: All effective MCMC sample sizes > %f]\n\n",world->essminimum); 
      break;
    case ' ':
      fprintf(file,"\n\n");
      break;
    }

  fprintf(file,"Locus  Replicate  Steps  Variance ratio (new/old variance)\n");
  fprintf(file,"-----  ---------  ------ ---------------------------------\n");
 for(z=0; z < world->loci * maxreplicate; z++)
    {
      if(world->burnin_stops[z].oldvariance > 0.0)
	{
	  fprintf(file,"%5li  %5li  %10li   %10.4f (%f/%f)\n",1 + world->burnin_stops[z].locus,
		  world->burnin_stops[z].replicate,
		  world->burnin_stops[z].stopstep,
		  world->burnin_stops[z].variance/world->burnin_stops[z].oldvariance,
		  world->burnin_stops[z].variance,
		  world->burnin_stops[z].oldvariance);
	}
      //			    world->burnin_stops[z].worker = myID;
    }
  fprintf(file,"\n\n");
#ifdef PRETTY
  if(file==EARTH->outfile)
    pdf_burnin_stops(EARTH, maxreplicate);
#endif
}

MYREAL combine_scaling_factor(world_fmt *world)
{ 
  const long np = world->numpop2;
  long pop;
  long i;
  MYREAL scaling_factor=0.0;
  bayes_fmt * bayes = world->bayes;
  for(i=0;i<np;i++)
    {
      if(bayes->map[i][1] == INVALID)
	continue;
      else
	{
	  pop  = bayes->map[i][1];
	}
      scaling_factor += exp(bayes->scaling_factors[pop] - bayes->maxmaxvala);
#ifdef DEBUG
      printf("%li  %li %g %g\n",i, pop, scaling_factor, bayes->maxmaxvala);
#endif
    }
  scaling_factor = log(scaling_factor) + bayes->maxmaxvala;
  if(world->options->has_bayesfile)
    {
#ifdef DEBUG
      printf("# Scaling factor %20.20f\n",scaling_factor);
#endif
      fprintf(world->bayesfile, "# Scaling factor %20.20f\n",scaling_factor);
    }
  return scaling_factor;
}

void print_bf_values(world_fmt * world)
{
  long locus;
  long t;
  const long hc = world->options->heated_chains;
  for(locus=0;locus<world->loci;locus++)
    {
      if(world->data->skiploci[locus])
	continue;
      printf("=====> locus=%li ",locus );
      for(t=0; t < hc; t++)
	{
	  printf("%f ",world->bf[locus*hc+t]);
	}
      printf("\n");
    }
}

#ifdef MPI_do_not_use
void fix_bayesfactor(world_fmt *world, option_fmt * options)
{
  return;
  long locus;
  long t;
  const long hc = world->options->heated_chains;
  const long maxreplicate = (options->replicate
		       && options->replicatenum >
		       0) ? options->replicatenum : 1;
  if(options->heating)//-----------------------heating
    {
      for(locus=0;locus<world->loci;locus++)
	{
	  if(world->data->skiploci[locus])
	    continue;
	  // the thermodynamic integration is reporting maxreplicate times too high values
	  // warning("this should fix the thermodynamic/mpi/replicate  problem, but needs to be checked for multilocus data");
	  printf("@@@@@@@ locus=%li %li (%f) ",locus,maxreplicate, world->bf[locus * hc + 0] );
	  for(t=0; t < hc; t++)
	    {
	      world->bf[locus * hc + t] /= maxreplicate;
	      printf("%f ",world->bf[locus*hc+t]);
	    }
	  printf("\n");
	}
    }
}
#endif


void print_bayesfactor(FILE *file, world_fmt **universe, option_fmt * options)
{
  //#ifdef MPI
  //static boolean done=FALSE;
  //#endif
  long t;
  const long hc = EARTH->options->heated_chains;
  long locus;
  const long maxreplicate = (options->replicate
		       && options->replicatenum >
		       0) ? options->replicatenum : 1;
  //  long lsteps = options->lsteps;
  //  MYREAL sum = 0.;
  MYREAL heat0 = 1.;
  MYREAL heat1 = 1.;
  MYREAL heat2 = 1.0;
  MYREAL val0  = 0.;
  MYREAL val1  = 0.;
  MYREAL val2  = 0.;
  MYREAL bfsum = 0.;
  MYREAL bfsum2 = 0.;
  MYREAL approxlsum = 0.;
  MYREAL hsum = 0.;
  MYREAL asum = 0.;
  MYREAL lsum;
  MYREAL lsum0;
  MYREAL scaling_factor = 0.0;
  // calculate the harmonic mean score
  for(locus=0;locus < EARTH->loci; locus++)
    {
      if(EARTH->data->skiploci[locus])
	continue;
      hsum +=  EARTH->hmscale[locus] - log (EARTH->hm[locus]);
    }

  fprintf(file,"\n\n\nLog-Probability of the data given the model (marginal likelihood = log(P(D|thisModel))\n");
  fprintf(file,"--------------------------------------------------------------------\n[Use this value for Bayes factor calculations:\nBF = Exp[log(P(D|thisModel) - log(P(D|otherModel)]\nshows the support for thisModel]\n\n");
#ifdef PRETTY
  if(file==EARTH->outfile)
    pdf_bayes_factor_header(EARTH,options);
#endif
  if(EARTH->loci>1)
    {
      fprintf(file,"\n\nLocus      Raw Thermodynamic score(1a)  Bezier approximated score(1b)     Harmonic mean(2)\n");
      fprintf(file,"------------------------------------------------------------------------------------------\n");
      if(file==EARTH->outfile)
	pdf_bayes_factor_rawscores_header(EARTH,options);
    }
  if(options->heating)//-----------------------heating
    {
      for(locus=0;locus<EARTH->loci;locus++)
	{
	  if(EARTH->data->skiploci[locus])
	    continue;
	  
	  lsum = 0.;
	  lsum0 = 0.;
	  // from very hot to cold, the 
	  for(t=1; t < hc; t++)
	    {
	      heat2 = heat0;
	      val2 = val0;
#ifdef MPI
	      heat0 = 1./options->heat[t-1] ;
	      heat1 = 1./options->heat[t];
#else
	      if(options->datatype!='g')
		{
		  if(options->adaptiveheat!=NOTADAPTIVE)
		    {
		      heat0 = 1./ universe[t-1]->averageheat ;
		      heat1 = 1./ universe[t]->averageheat;
		    }
		  else
		    {
		      heat0 = universe[t-1]->heat ;
		      heat1 = universe[t]->heat;
		    }
		}
	      else
		{
		  heat0 = 1./options->heat[t-1];
		  heat1 = 1./options->heat[t];
		}
#endif
	      //Simpson's rule and trapezoidal are the same when I do only have function values at a and b
	      val0 = EARTH->bf[locus * hc + t-1];
	      val1 = EARTH->bf[locus * hc + t];
	      //we keep last element to adjust for Bezier approximation
	      lsum0 = (heat0 - heat1) * (val0 + val1) * 0.5;
	      lsum += lsum0; 
#ifdef DEBUG
	      //	      printf("%i> sum[%li]=%f temp=(%f %f) values=(%f %f)\n",myID,t,lsum,heat0,heat1,EARTH->bf[locus * hc + t],EARTH->bf[locus * hc + t]);
#endif
	    }
	  //	  (x2 y1 - x1 y2)/(x1 - x2)
	  // this last addition to the lsum calculates the chunk between the last temperature and 
	  // the infinitely hot temperature as a linear approximation of the the second hottest temperature
	  // this is certainly rough, but in simulations with 3 populations one can see that with large number
	  // of temperatures this looks reasonable, and one can approximate the integral more accurately with 
	  // with only a few columns.
	  // using Bezier to approximate nice curve between the last two points to mimick the curve that
	  // can be found with 16 or 32 heated chains, handle points are calculated using adhoc decisions
	  // (comparison with 32 heated chains) using 0.8 of the interval for handle_y1 and the intercept
	  // between the first and the second last point to calculate the handle_y2 see sumbezier()
	  // for implementation. Currently this is not tunable.
	  //	  approxlsum = sumbezier(100, heat1, EARTH->bf[locus * hc + t-1], 
	  //			 heat0, EARTH->bf[locus * hc + t-2], 
	  //			 1.0, EARTH->bf[locus * hc]);
	  approxlsum = sumbezier(100, heat1, val1, 
				 heat0, val0, 
				 heat2, val2);
	  //#ifdef DEBUG
	  //	  fprintf(stdout,"Thermo[%li]=(%f) - (%f) + (%f) = (%f)\n",locus, lsum, lsum0, approxlsum, approxlsum + lsum - lsum0);
	  //	  #endif
	  if(EARTH->loci>1)
	    {
	      fprintf(file,"  %5li%22.2f        %22.2f           %12.2f\n", locus + 1, lsum, lsum-lsum0+approxlsum, EARTH->hmscale[locus] - log (EARTH->hm[locus])); 
#ifdef PRETTY
	      if(file==EARTH->outfile)
		pdf_bayes_factor_rawscores(locus, lsum, lsum-lsum0+approxlsum, EARTH->hmscale[locus] - log (EARTH->hm[locus]));
#endif
	    }
	  bfsum2 += approxlsum + lsum - lsum0;  	  
	  bfsum += lsum; //+ EARTH->bfscale[locus];
	}
      if(EARTH->loci>1)
	{
	  scaling_factor = combine_scaling_factor(EARTH);
	  bfsum  += scaling_factor;
	  bfsum2 += scaling_factor;
	  hsum   += scaling_factor;
	  fprintf(file,"---------------------------------------------------------------------------------------\n");
	  fprintf(file,"  All  %22.2f        %22.2f           %12.2f\n[Scaling factor = %f]\n", bfsum, bfsum2, hsum, scaling_factor); 
#ifdef PRETTY
	  if(file==EARTH->outfile)
	    pdf_bayes_factor_rawscores(-1, bfsum, bfsum2, hsum);
#endif
	}
      else 
	{
	  fprintf(file,"\n\n(1a) Thermodynamic integration: log(Prob(D|Model))=%f (%f with Bezier-approximation[1b]) %s\n", bfsum, bfsum2,  
		  options->adaptiveheat!=NOTADAPTIVE ? "(adaptive heating -> value may most likely be wrong)" : "");
	  fprintf(file,"(2) Harmonic mean:             log(Prob(D|Model))=%f\n",hsum);
	  fprintf(file,"(1) and (2) should give a similar result, (2) is considered more\ncrude than (1), but (1) needs heating with several well-spaced chains\n\n");
	}
    }
  else   // -----------------------------------------not heating
    {
      for(locus=0;locus<EARTH->loci;locus++)
	{
	  if(EARTH->data->skiploci[locus])
	    continue;
	  if(EARTH->loci>1)
	    {
	      fprintf(file,"    %5li           ------                           ------            %12.2f\n", locus+1, EARTH->hmscale[locus] - log (EARTH->hm[locus])); 
#ifdef PRETTY
	      if(file==EARTH->outfile)
		pdf_bayes_factor_rawscores_harmo(locus, EARTH->hmscale[locus] - log (EARTH->hm[locus]));
#endif
	    }
	}
      if(EARTH->loci>1)
	{
	  scaling_factor = combine_scaling_factor(EARTH);
	  hsum += scaling_factor;
	  fprintf(file,"---------------------------------------------------------------------------------------\n");
	  fprintf(file,"      All           ------                           ------            %12.2f\n", hsum); 
#ifdef PRETTY
	  if(file==EARTH->outfile)
	    pdf_bayes_factor_rawscores_harmo(-1, hsum);
#endif
	}
      else
	{
	  fprintf(file,"\n\n(1) Thermodynamic integration: log(Prob(D|Model))= UNAVAILABLE because no heating was used\n");
	  fprintf(file,"(2) Harmonic mean:             log(Prob(D|Model))=%f\n",hsum);
	  fprintf(file,"(1) and (2) should give a similar result, (2) is considered more\ncrude than (1), but (1) needs heating with several well-spaced chains\n\n");
	}
    }
#ifdef PRETTY
  if(file==EARTH->outfile)
    pdf_bayes_factor(EARTH, bfsum, bfsum2, hsum, asum, maxreplicate,scaling_factor);
#endif
}

///
/// integrates over a Bezier curve between two points
/// calculates two handle points that are set to adhoc values
/// so that the x values of the handle are the the x value of the lowest point
/// and the y values are set to about 80% of the min to max interval for the left point
/// and a value that is the the y value from ax + b where a is calculated from a 
/// third point to the right and the second point, the third point is not used for the
/// the Bezier curve otherwise
MYREAL sumbezier(long intervals, MYREAL x0, MYREAL y0, MYREAL x1, MYREAL y1, MYREAL x2, MYREAL y2)
{
  const MYREAL inv_interval = 1./intervals;
  const MYREAL sx0 = x0;
  const MYREAL sx1 = x0;
  const MYREAL sy0 = 0.2 * y0 + 0.8 * y2;
  const MYREAL sy1 = (-x2 * y1 + x1 * y2)/(x1 - x2);
  MYREAL t     = 0.0;
  MYREAL t2    = 0.0;
  MYREAL t3    = 0.0;
  MYREAL onet  = 1.0;
  MYREAL onet2 = 1.0;
  MYREAL onet3 = 1.0;
  MYREAL newx;
  MYREAL newy;
  MYREAL oldx;
  MYREAL oldy;
  MYREAL sum = 0.0;
  // integrate over intervals between x0 and x1 and return sum
  // intialize with t=0.0
  oldx = x0;
  oldy = y0;
  //fprintf(stdout,"\n\n%f %f %f %f %f %f\n",x2,y2,sx0,sy0, x0,y0);
  for(t=inv_interval; t <= 1.0; t += inv_interval)
    {
      onet  = 1.0 - t;
      onet2 = onet * onet;
      onet3 = onet2 * onet;
      onet2 *= 3.0 * t;
      t2 = t * t;
      t3 = t2 * t;
      t2 *=  3.0 * onet;
      //      newx = 3sx0 (1-t)^2 t + 3 sx1 (1-t) t^2 + (1-t)^3 x0 + t^3 x1
      newx = sx0 * onet2 + sx1 * t2 + onet3 * x0 + t3 * x1;
      newy = sy0 * onet2 + sy1 * t2 + onet3 * y0 + t3 * y1;
      //      fprintf(stdout,"%f %f\n",newx,newy);
      sum += (newx - oldx) * (newy + oldy)/2.;
      oldx = newx;
      oldy = newy;
    }
  //  fprintf(stdout,"%f %f %f %f\n\n\n",x1,y1,sx1,sy1);
  return sum;
}

  
/// calculate values for the marginal likelihood using thermodynamic integration
/// based on a method by Friel and Pettitt 2005
/// (http://www.stats.gla.ac.uk/research/TechRep2005/05.10.pdf)
/// this is the same method described in Lartillot and Phillippe 2006 Syst Bio
/// integrate over all temperature using a simple trapezoidal rule
/// prob(D|model) is only accurate with intervals for temperatures from 1/0 to 1/1.
/// reports also the harmonic mean
void calculate_BF(world_fmt **universe, option_fmt *options)
{
  long i;
  MYREAL xx, xx2;
  //  MYREAL logprior;
#ifdef THERMOCHECK
  static long counter=0;
#endif
  long locus = universe[0]->locus;
  long hc = options->heated_chains;
  if(universe[0]->data->skiploci[locus])
    return;
  if(EARTH->likelihood[EARTH->G] <= -HUGE)
    return;
  EARTH->am[locus] += 1;
  locus = EARTH->locus;
  xx = EARTH->likelihood[EARTH->G];

#ifdef THERMOCHECK
  counter++;
#endif

  if(xx <= -HUGE)
    return;

  if(xx > EARTH->hmscale[locus])
    {
      //      EARTH->hm[locus] += EXP(EARTH->hmscale[locus] - xx);
      xx2 = EXP(EARTH->hmscale[locus] - xx);
      EARTH->hm[locus] += (xx2 - EARTH->hm[locus])/ (EARTH->am[locus]);
    }
  else
    {
      EARTH->hm[locus] *= EXP(xx - EARTH->hmscale[locus]);
      //      EARTH->hm[locus] *= xx - EARTH->hmscale[locus];
      EARTH->hmscale[locus] = xx;
      //      EARTH->hm[locus] += 1.;
      EARTH->hm[locus] += (1. - EARTH->hm[locus])/ (EARTH->am[locus]);
    }

#ifdef THERMOCHECK
  if(options->mixplot)
    fprintf(universe[0]->options->mixfile,"@BF %li %li %10.10f %10.10f",counter, locus, EARTH->heat, xx);
  //, EARTH->hmscale[locus]);
#endif

  if(options->heating)
    {
      for (i = 0; i < hc; i++)
	{
	  xx = universe[i]->likelihood[universe[i]->G];
#ifdef THERMOCHECK
	  if(options->mixplot)
	    fprintf(universe[0]->options->mixfile," @T %10.10f %10.10f", universe[i]->heat, xx);
#endif
	  EARTH->bf[locus * hc + i] += (xx - EARTH->bf[locus * hc + i])/ (EARTH->am[locus]);
	  //EARTH->bf[locus * hc + i] += (xx - EARTH->bf[locus * hc + i]);
	  //#ifdef THERMOCHECK
	  //fprintf(universe[0]->options->mixfile,"%f\n", 	  EARTH->bf[locus * hc + i]);
	  //#endif
	}
    }
#ifdef THERMOCHECK
  if(options->mixplot)
  fprintf(universe[0]->options->mixfile,"\n");
#endif
}

#ifdef BFDEBUG
void  print_marginal_order(char *buf, long *bufsize, world_fmt *world)
{
  long i;

  for(i=0;i<world->options->heated_chains;i++)
    *bufsize += sprintf(buf+ *bufsize,"# --  %s = %f\n", "Thermodynamic temperature", world->options->heat[i]);
  *bufsize += sprintf(buf+ *bufsize,"# --  %s\n", "Marginal log(likelihood) [Thermodynamic integration]");
  *bufsize += sprintf(buf+ *bufsize,"# --  %s\n", "Marginal log(likelihood) [Harmonic mean]");
}

#ifdef MPI /* */
void      print_marginal_like(float *temp, long *z, world_fmt * world)
{
  long locus = world->locus;
  long t;
  long hc = world->options->heated_chains; 
  MYREAL lsum; 
  MYREAL heat0, heat1;

  if(world->options->heating)
    {
      lsum = 0.;
      for(t=1; t < hc; t++)
	{
	  heat0 = 1./world->options->heat[t-1] ;
	  heat1 = 1./world->options->heat[t];
	  // this ignores adaptive heating for MPI!!!!
	  temp[*z] = (float) world->bf[locus * hc + t-1];
	  *z += 1;
	  lsum += (heat0 - heat1) * ((world->bf[locus * hc + t-1] + world->bf[locus * hc + t]) * 0.5);
	}
      temp[(*z)++] =  (float) world->bf[locus * hc + t-1];
      temp[(*z)++] =  (float) lsum;
    }
  temp[(*z)++] =  (float) world->hmscale[locus] - log(world->hm[locus]);
}
#else /*not MPI*/
void      print_marginal_like(char *temp, long *c, world_fmt * world)
{
  long locus = world->locus;
  long t;
  long hc = world->options->heated_chains;  
  MYREAL lsum;
  MYREAL heat0, heat1;
  if(world->options->heating)
    {
      lsum = 0.;
      for(t=1; t < hc; t++)
	{
	  if(world->options->adaptiveheat!=NOTADAPTIVE)
	    {
	      heat0 = world->options->averageheat[t-1] ;
	      heat1 = world->options->averageheat[t];
	    }
	  else
	    {
	      heat0 = 1./ world->options->heat[t-1] ;
	      heat1 = 1./ world->options->heat[t];
	    }
	  *c += sprintf(temp+ *c,"\t%f", world->bf[locus * hc + t-1]);
	  lsum += (heat0 - heat1) * ((world->bf[locus * hc + t-1] + world->bf[locus * hc + t]) * 0.5);
	}
      *c += sprintf(temp + *c,"\t%f", world->bf[locus * hc + t-1]);
      *c += sprintf(temp + *c,"\t%f", lsum);
    }
  *c += sprintf(temp + *c,"\t%f", world->hmscale[locus] - log(world->hm[locus]));
}
#endif /*not MPI*/
#endif /*BFDEBUG*/

#ifdef LONGSUM
/// put timepoints into the total tree to allow for bottlenecks and growth 
long
find_chaincutoffs (MYREAL * treetime, world_fmt * world, long locus, long r)
{
    
    long j;
    long T;
    MYREAL totaltime = 0.;
    long copies = 0;
    tarchive_fmt *atl;
    atl = world->atl[r][locus].tl;
    for (j = 0; j < world->atl[r][locus].T; j++)
    {
        T = atl[j].longsumlen - 1;
        copies += atl[j].copies;
        totaltime = atl[j].longsum[T].eventtime / 3.;
        treetime[0] += totaltime;
        treetime[1] += totaltime + totaltime;
        treetime[2] += totaltime + totaltime + totaltime;
    }
    return copies;
}

/// find the cutoffs for the locus, averaging over multiple replicates
void
find_loci_cutoffs (MYREAL * treetime, world_fmt * world)
{
    long r;
    long locus;
    long count = 0;
    long repstart = 0;
    long repstop = 1;
    memset (treetime, 0, sizeof (MYREAL) * 3);
    set_replicates (world, world->repkind, world->rep, &repstart, &repstop);
    
    for (locus = 0; locus < world->loci; locus++)
    {
        for (r = repstart; r < repstop; r++)
        {
            count += find_chaincutoffs (treetime, world, locus, r);
        }
    }
    treetime[0] /= (MYREAL) count;
    treetime[1] /= (MYREAL) count;
    treetime[2] /= (MYREAL) count;
}


/// find the cutoffs avaergin over all replicates 
void
find_replicate_cutoffs (MYREAL * treetime, world_fmt * world)
{
    long r;
    long count = 0;
    long repstart = 0;
    long repstop = 1;
    memset (treetime, 0, sizeof (MYREAL) * 3);
    set_replicates (world, world->repkind, world->rep, &repstart, &repstop);
    
    for (r = repstart; r < repstop; r++)
    {
        count += find_chaincutoffs (treetime, world, world->locus, r);
    }
    treetime[0] /= (MYREAL) count;
    treetime[1] /= (MYREAL) count;
    treetime[2] /= (MYREAL) count;
}


/// change the cutoff times
void
change_longsum_times (world_fmt * world)
{
    long i;
    long numpop3 = 3 * world->numpop;
    long numpop6 = 2 * numpop3;
    MYREAL *treetime;
    treetime = (MYREAL *) mycalloc (3, sizeof (MYREAL));	// we have 3 classes
    if (world->locus < world->loci)
        find_replicate_cutoffs (treetime, world);
    else
        find_loci_cutoffs (treetime, world);
    for (i = numpop3; i < numpop6; i += 3)
    {
        world->flucrates[i] = treetime[0];
        world->flucrates[i + 1] = treetime[1];
        world->flucrates[i + 2] = treetime[2];
    }
    myfree(treetime);
}

#endif /*LONGSUM*/

///
/// check_parmfile checks what is called and also whether we should honor the menu=YES in the file
/// user uses another parmfile than the default when there is a parameter associated with the program call
/// this function checks the first two line of the file to find out whether this is a parmfile or potentially
/// a datafile/ 
/// The addition of an option to override the menu command in the parmfile
/// migrate-n -menu ==> use menu independent of the parmfile setting
/// migrate-n -nomenu ==> does use menu even if parmfile says so

boolean  check_parmfile(long argcount, char **arguments, char *parmfilename)
{
  int argument=0;
  boolean  usemenu = OVERRIDE_NO;
  if (argcount > 1)			
    {
      argument = 1;
      while(argument < argcount)
	{
	  //fprintf(stderr,"<|%s|>\n",argv[argument]);
	  //fflush(stderr);
	  if(arguments[argument][0]!='-')
	    {
	      strncpy (parmfilename, arguments[argument], 
		       (long) strlen (arguments[argument]) + 1);
	    }
	  else
	    {
	      if(strcmp(arguments[argument],"-menu")==0)
		{
		  usemenu = OVERRIDE_MENU;
		}
	      else
		{
		  if(strcmp(arguments[argument],"-nomenu")==0)
		    {
		      usemenu = OVERRIDE_NOMENU;
		    }
		}
	    }
	  argument++;
	}
    }
  return usemenu;
}

/// set the menu overrider 
boolean set_usemenu(boolean usemenu, boolean fromparmfile)
{
  switch(usemenu)
    {
    case OVERRIDE_MENU:
      return TRUE;
      break;
    case OVERRIDE_NOMENU:
      return FALSE;
      break;
    default:
      break;
    }
  return fromparmfile;
} 

/// checks the settings of the number of long an short chain for bayes options and resets useless settings
void check_bayes_options(option_fmt *options)
{
  if (options->bayes_infer)
    {
      if(options->schains>0)
	{
	  //		fprintf(stdout,"NOTICE: with Bayesian inference no short chains are allowed, setting them to ZERO\n");
	  options->schains = 0;
	}
      if(options->lchains>1 || options->lchains < 1)
	{
	  fprintf(stdout,"NOTICE: with Bayesian inference, most efficient runs are with 1 long chain,\nsetting the long chain to 1\n");
	  options->lchains = 1;
	};
    }
}

void reset_bayesmdimfile(world_fmt *world, option_fmt *options)
{
  long locus=0;
  if(world->options->has_bayesmdimfile)
    {
      get_bayeshist (world, options); // this is transferring some material but
      // the parameter file is not filled in anymore [-> small data block]
      // save parameters in file and reload
#ifdef MPI
      // close all for bayesallfile
#ifdef PARALIO
      MPI_File_close(&world->mpi_bayesmdimfile);
#endif
#endif
      if(myID==MASTER)
	{
#ifdef MPI
#ifdef PARALIO
	  world->bayesmdimfile = fopen(options->bayesmdimfilename,"r+");
#else
#ifdef ZNZ
	  znzclose(world->bayesmdimfile);
	  world->bayesmdimfile = znzopen(options->bayesmdimfilename,"r", options->use_compressed);
#else
	  fflush(world->bayesmdimfile);
	  world->bayesmdimfile = freopen(options->bayesmdimfilename,"r+", world->bayesmdimfile);
#endif
#endif
#else
#ifdef ZNZ
	  znzclose(world->bayesmdimfile);
	  world->bayesmdimfile = znzopen(options->bayesmdimfilename,"r", options->use_compressed);
#else
	  fflush(world->bayesmdimfile);
	  world->bayesmdimfile = freopen(options->bayesmdimfilename,"r+", world->bayesmdimfile);
#endif
#endif
	  read_bayes_fromfile(world->bayesmdimfile,world, options);
	  for(locus=0;locus<world->loci;locus++)
	    {
	      if(world->data->skiploci[locus])
		continue;
	      if(world->bayes->histogram[locus].results == NULL)
		{
		  world->data->skiploci[locus] = TRUE;
		  continue;
		}
	      calc_hpd_credibility(world->bayes, locus, world->numpop2, world->numpop2 + world->bayes->mu);
	    }
	}
    }
  else
    {
      // keep everything in RAM
      get_bayeshist (world, options);
    }
}
