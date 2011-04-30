/* ! \file menu.c */
/*------------------------------------------------------
  Maximum likelihood estimation
  of migration rate  and effectice population size
  using a Metropolis-Hastings Monte Carlo algorithm
  -------------------------------------------------------
  M E N U   R O U T I N E S

  presents the menu and its submenus.
  Peter Beerli 1996, Seattle
  beerli@fsu.edu

  Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
  Copyright 2003-2007 Peter Beerli, Tallahassee FL

  This software is distributed free of charge for non-commercial use
  and is copyrighted. Of course, we do not guarantee that the software
  works and are not responsible for any damage you may cause or have.

  $Id: menu.c 1828 2011-03-20 18:53:32Z beerli $

  -------------------------------------------------------*/
#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#include "sequence.h"
#include "fst.h"
#include "options.h"
#include "migrate_mpi.h"
#include "menu.h"
#include "pretty.h"

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
/* prototypes ------------------------------------------- */
//void            print_menu_title(FILE * file, option_fmt * options);
//void            print_menu_accratio(long a, long b, world_fmt * world);
//long            print_title(world_fmt * world, option_fmt * options);
//void            get_menu(option_fmt * options);
/* private functions */
void            setup_datatype(char *datatype, option_fmt * options);
void            menuData(option_fmt * options, char datatype[]);
void            menuInput(option_fmt * options);
void            menuParameters(option_fmt * options);
void            menuStrategy(option_fmt * options);
void            menuHeat(option_fmt * options, char *input);
void            display_ml_mcmc(option_fmt * options);
void            display_bayes_mcmc(option_fmt * options);
//void          menuSequences(option_fmt * options);
void            read_custom_menu_migration(option_fmt * options);
void            read_custom_menu_lratio(option_fmt * options);
char           *custom_migration_type(long type);
void            read_heatvalues(option_fmt * options);
void            menuRandom(MYREAL * param, char type);
void            start_tree_method(option_fmt * options);
void            start_data_method(option_fmt * options);
void            how_many_pop(long *numpop);
void            set_menu_localities(option_fmt *options);
void            set_localities_string(char *loc, option_fmt *options);
char            dialog_lrt(long numpop, lratio_fmt * lratio);
char            menuread_lrt_paramvalue(char *value, long *counter, long numpop2);
void            get_plotmenu(option_fmt * options);
boolean         menuStrategy_ml(option_fmt * options);
boolean         menuStrategy_bayes(option_fmt * options);
long            get_prior(char *input);
void            set_prior(char *output, int * prior, boolean without_rate);
void set_proposal(char *output, boolean *proposal, boolean without_rate);
char * submodeltype(int type);

#ifdef POPMODEL
#include "popmodel.h"
#endif
#define MI_INFILE 1
#define MI_RAND 2
#define MI_TITLE  3
#define MI_SUMREAD 4
#define MI_PROGRESS 5
#define MI_PRINTDATA 6
#define MI_OUTFILE 7
#define MI_PLOT 8
#define MI_PROFILE 9
#define MI_LRT 10
#define MI_AIC 11
#define MI_TREES 12
#define MI_PLOTCOORD 13
#define MI_SUMFILE 14
#define MI_LOGFILE 15
#define MI_UEPFILE 16
#define MI_UEPRATE 17
#define MI_UEPFREQ 18
#define MI_MIGHISTOGRAM 19
#define MI_SKYLINE 20

///
///Enumeration of possible population models
enum popmodel_menu 
{
  WRIGHTFISHER, CANNING, MORAN
};

///
///Enumeration of possible choices for ML menu
enum ml_menu 
{
  MLSTRATEGY, MLSHORTCHAINS, MLSHORTSKIP, MLSHORTSAMPLES, MLLONGCHAINS, MLLONGSKIP, MLLONGSAMPLES,
  MLBURNIN, MLREPLICATE, MLHEAT, MLMOVINGSTEPS, MLEPSILON, MLGELMAN
};

///
///Enumeration of possible choices for Datatype menu
enum dtseq_menu 
{
  DTSEQZERO, DTSEQTYPE, DTSEQTRATIO, DTSEQFREQ, DTSEQSITECATEGS, DTSEQRATES, DTSEQCORR, DTSEQWEIGHT, DTSEQINTERLEAVED, DTSEQERROR, DTSEQFAST, DTSEQTREE, DTSEQINHERITANCE, DTSEQRANDOMSUBSET, DTSEQTIPDATE, DTSEQMUTRATE, DTSEQGENERATION 
};

enum dtmsat_menu {
  DTMSATZERO, DTMSATTYPE, DTMSATMISS, DTMSATTRESH, DTMSATTREE, DTMSATINHERITANCE, DTMSATRANDOMSUBSET, DTMSATTIPDATE, DTMSATMUTRATE, DTMSATGENERATION
};

enum dtbrown_menu 
{
  DTBROWNZERO, DTBROWNTYPE, DTBROWNMISS, DTBROWNTREE, DTBROWNINHERITANCE, DTBROWNRANDOMSUBSET, DTBROWNTIPDATE, DTBROWNMUTRATE, DTBROWNGENERATION
};
enum dtep_menu {
  DTEPZERO, DTEPTYPE, DTEPMISS, DTEPTREE, DTEPINHERITANCE, DTEPRANDOMSUBSET, DTEPTIPDATE, DTEPMUTRATE, DTEPGENERATION
};

///
///Enumeration of possible choices for Bayes menu
enum bayesian_menu {
  BAYESSTRATEGY, BAYESOUT, BAYESMDIMOUT, BAYESBINNING, BAYESPRETTY, BAYESFREQ, BAYESPROPOSAL, BAYESPRIOR, BAYESLCHAINS, BAYESSKIP, BAYESSAMPLES, BAYESBURNIN,
  BAYESREPLICATE, BAYESHEAT, BAYESMOVINGSTEPS, BAYESGELMAN, BAYESPRIORALONE 
};

///
///get user supplied filenames or fill the filenames with default values
void
get_filename(char message[], char thedefault[], char *filename)
{
  char            input[LINESIZE];
  printf("%s\n[Default: %s]\n===> ", message, thedefault);
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (input[0] == '\0')
    strcpy(filename, thedefault);
  else
    strcpy(filename, input);
}

///
///sets the plotoptions default is NO PLOT
void
get_plotoptions(option_fmt * options)
{
  char            input[LINESIZE];
  do {
    printf("  Plot Likelihood surface:\n");
    printf
      ("  (B)oth to outfile and mathfile, (O)utfile only, (N)o plot\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
  }
  while (strchr("BON", (int) uppercase(input[0])) == NULL);
  switch (uppercase(input[0])) {
  case 'B':
    options->plotmethod = PLOTALL;
    break;
  case 'O':
    options->plotmethod = PLOTOUTFILE;
    break;
  case 'N':
  default:
    options->plot = FALSE;
    return;
  }
  get_plotmenu(options);
}

///
///allows the setting of the tickmarks and ranges for the plots
///[remember these produce only ASCII plots
void get_plotmenu(option_fmt * options) {
  char            input[LINESIZE];
  char            answer = 'N';
  do {
    printf("   Current plot settings are:\n");
    printf
      ("              Parameter: %s, Scale: %s, Intervals: %li\n",
       options->plotvar == PLOT4NM ? "{Theta, xNm}" : "{Theta, M}",
       options->plotscale == PLOTSCALELOG ? "Log10" : "Standard",
       options->plotintervals);
    printf("              Ranges: X-%5.5s: %f - %f\n",
	   options->plotvar == PLOT4NM ? "xNm" : "M",	   options->plotrange[0], options->plotrange[1]);
    printf("              Ranges: Y-%5.5s: %f - %f\n",
	   "Theta", options->plotrange[2],
	   options->plotrange[3]);
    printf(" Is this OK? [YES or No]\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    answer = uppercase(input[0]);
    if (answer != 'Y') {
      printf("1    Use as X-axis: xNm  Y-axis: Theta\n");
      printf("2    Use as X-axis: M    Y-axis: Theta\n\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (input[0] == '2')
	options->plotvar = PLOTM;
      else
	options->plotvar = PLOT4NM;
      printf("1     Scale is standard \n");
      printf("2     Scale is in Log10\n\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (input[0] == '2')
	options->plotscale = PLOTSCALELOG;
      else
	options->plotscale = PLOTSCALESTD;
      do {
	printf("How many plot positions?\n");
	printf
	  ("[Default is 36, this changes only in the mathfile]\n\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (strlen(input) > 1)
	  options->plotintervals = ATOL(input);
      }
      while (options->plotintervals < 2);
      printf("Change the ranges of the axes\n");
      printf("              Ranges: X-%5.5s: %f - %f\n",
	     options->plotvar == PLOT4NM ? "xNm" : "M",
	     options->plotrange[0],
	     options->plotrange[1]);
      printf("              Ranges: Y-%5.5s: %f - %f\n",
	     "Theta", options->plotrange[2],
	     options->plotrange[3]);
      printf
	("Give lowest and highest value for X axis:\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (strlen(input) > 1) {
#ifdef USE_MYREAL_FLOAT
	sscanf(input, "%f%f", &options->plotrange[0],
	       &options->plotrange[1]);
#else
	sscanf(input, "%lf%lf", &options->plotrange[0],
	       &options->plotrange[1]);
#endif
      }
      printf
	("Give lowest and highest value for Y axis:\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (strlen(input) > 1) {
#ifdef USE_MYREAL_FLOAT
	sscanf(input, "%f%f", &options->plotrange[2],
	       &options->plotrange[3]);
#else
	sscanf(input, "%lf%lf", &options->plotrange[2],
	       &options->plotrange[3]);
#endif
      }
    }
  }
  while (answer != 'Y');
}


void
print_menu_title(FILE * file, option_fmt * options)
{
  char            nowstr[LINESIZE];
  if (options->menu || options->progress) {
    FPRINTF(file, "  =============================================\n");
    FPRINTF(file, "  MIGRATION RATE AND POPULATION SIZE ESTIMATION\n");
    FPRINTF(file, "  using Markov Chain Monte Carlo simulation\n");
    FPRINTF(file, "  =============================================\n");
#ifdef MPI
    FPRINTF(file, "  Compiled for a PARALLEL COMPUTER ARCHITECTURE\n");
    FPRINTF(file, "  One master and %i compute nodes are available.\n", numcpu-1);
#endif
#ifdef THREAD
    FPRINTF(file, "  Compiled for a SYMMETRIC MULTIPROCESSORS\n");
#endif
#ifdef ALTIVEC
    FPRINTF(file, "  ALTIVEC enabled\n");
#endif
#ifdef FAST_EXP
    FPRINTF(file, "  Fast approximation to Exp() and Log() used\n");
#endif
#ifdef PRETTY
    FPRINTF(file, "  PDF output enabled [%s]\n",
#ifdef A4PAPER
	    "A4-size"
#else
	    "Letter-size"
#endif
);
#endif

    FPRINTF(file, "  Version %s ", MIGRATEVERSION);
    if(strlen(MIGRATESUBVERSION)>0)
      FPRINTF(file, "  [%s]\n", MIGRATESUBVERSION);
    else
      FPRINTF(file, "\n");
    get_time(nowstr, "  %c");
    if (nowstr[0] != '\0')
      FPRINTF(file, "  Program started at %s\n", nowstr);
    FPRINTF(file, "\n\n");
  }
}

void
print_menu_accratio(long a, long b, world_fmt * world)
{
  char            buffer[STRSIZE];
  boolean         writelog = world->options->writelog;
  boolean         progress = world->options->progress;

  if (writelog || progress) 
    {
      sprintf(buffer, "           Acceptance-ratio = %li/%li (%f)\n\n", a, b,
	      (MYREAL) a / (MYREAL) b);
      if (progress)
	FPRINTF(stdout, "%s", buffer);
      if (writelog)
	FPRINTF(world->options->logfile, "%s", buffer);
    }
}

long
print_title(world_fmt * world, option_fmt * options)
{
  char            nowstr[LINESIZE];
  long            len = 45, i, filepos = -1;
  if (!world->options->simulation) 
    {
      FPRINTF(world->outfile, "  ");
      if (options->title[0] != '\0') 
	{
	  len = (long) MAX(strlen(options->title), 45);
	  for (i = 0; i < len; i++)
	    FPRINTF(world->outfile, "=");
	  FPRINTF(world->outfile, "\n  %-*.*s\n  ", (int) len, (int) len,options->title);
	}
      for (i = 0; i < len; i++)
	FPRINTF(world->outfile, "=");
      
      FPRINTF(world->outfile,
	      "\n  MIGRATION RATE AND POPULATION SIZE ESTIMATION\n");
      FPRINTF(world->outfile,
	      "  using Markov Chain Monte Carlo simulation\n  ");
      
      for (i = 0; i < len; i++)
	FPRINTF(world->outfile, "=");
      FPRINTF(world->outfile, "\n  Version %s\n", MIGRATEVERSION);
      get_time(nowstr, "%c");
      if (nowstr[0] != '\0') {
	FPRINTF(world->outfile, "\n  Program started at %s\n", nowstr);
	filepos = ftell(world->outfile);
	/*
	 * DO NOT REMOVE the blanks in the following print
	 * statement, they are overwritten at the program end
	 * with the timestamp
	 */
	FPRINTF(world->outfile,
		"                                                   \n");
      }
      FPRINTF(world->outfile, "\n\n");
    }
  return filepos;
}

///
/// main menu display, if counter is bigger than 100 the program will exit, to 
/// make sure that batch jobs are not filling up the available harddisk space
void
get_menu(option_fmt * options, world_fmt *world, data_fmt *data)
{
  long           counter=0L;
  char            input[LINESIZE];
  char           *datatype;
  datatype = (char *) mycalloc(LINESIZE, sizeof(char));
  if (options->menu) {
    setup_datatype(datatype, options);
    do {
      printf("  Settings for this run:\n");
      printf("  D       Data type currently set to: %-30.30s\n", datatype);
      printf("  I       Input/Output formats\n");
      printf("  P       Parameters  [start, migration model]\n");
      printf("  S       Search strategy\n");
      printf("  W       Write a parmfile\n");
      printf("  Q       Quit the program\n");
      printf("\n\n");
      printf("  To change the settings type the letter for the menu to change\n");
      printf("  Start the program with typing Yes or Y\n===> ");fflush(stdout);
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      printf("\n");
      switch (uppercase(input[0])) {
      case 'D':
	menuData(options, datatype);
	break;
      case 'I':
	menuInput(options);
	break;
      case 'P':
	menuParameters(options);
	break;
      case 'S':
	menuStrategy(options);
	break;
      case 'W':
	save_parmfile(options, world, data);
	break;
      case 'Q':
#ifdef MPI
        MPI_Finalize ();
#endif /*MPI*/
	exit(0);
      default:
	break;
      }
      if(counter++ > 100)
	{
	  warning("Batch job uses menu ON option, and gets non-sensical input, set menu=NO in parmfile\nProgram continues but check!!!!!");
	  input[0]='Y';
	  break;
	}
    }
    while (uppercase(input[0]) != 'Y');
    if (options->usertree && !strchr(SEQUENCETYPES, options->datatype))
      options->usertree = FALSE;
    if (options->dist && !strchr(SEQUENCETYPES, options->datatype))
      options->dist = FALSE;
    //prevents incompatability
    // of tree and distance designation, only for
    //sequence data the datalines and tree - tips can match
    // msats or EP data is ambiguos concerning the tree.
      }
  if (options->datatype == 'g')
    //prevents overwriting sumfile
    options->writesum = FALSE;

#ifdef UEP

  if (options->uep)
    options->randomtree = TRUE;
#endif
  myfree(datatype);
}
/*
 * private
 * functions---------------------------------------------------
 */

boolean         setup_categs(option_fmt * options) {
  boolean         didchangecat = FALSE;
  if              (options->categs > 1)
    didchangecat = TRUE;
  else
    didchangecat = FALSE;
  return didchangecat;
}

boolean         setup_rcategs(option_fmt * options) {
  boolean         didchangercat = FALSE;
  if              (options->rcategs > 1)
    didchangercat = TRUE;
  else
    didchangercat = FALSE;

  return didchangercat;
}

void            setup_datatype(char *datatype, option_fmt * options) {
  switch (options->datatype) {
  case 'a':
    strcpy(datatype, "'infinite' allele model");
    break;
  case 'm':
    if(options->msat_option == SINGLESTEP)
      strcpy(datatype, "Singlestep mutation model");
    else
      strcpy(datatype, "Multistep mutation model");
    break;
  case 'b':
    strcpy(datatype, "Brownian motion model");
    break;
  case 's':
    strcpy(datatype, "DNA sequence model");
    break;
  case 'n':
    strcpy(datatype, "SNP model");
    break;
  case 'h':
    strcpy(datatype, "SNP model (Hapmap data)");
    break;
    //case 'u':
    //strcpy(datatype, "Panel-SNP (panel is population 1)");
    //break;
    //case 'f':
    //strcpy(datatype, "ancestral state reconstruction method");
    //break;
  case 'g':
    strcpy(datatype, "Genealogy summary");
    options->readsum = TRUE;
    break;
  default:
    if (options->readsum) {
      options->datatype = 'g';
      strcpy(datatype, "Genealogy summary");
    } else {
      options->datatype = 's';
      strcpy(datatype, "DNA sequence model");
    }
    break;
  }
}

///
/// set up start tree
void setup_starttree(char *starttree, option_fmt * options) {
  if (options->usertree && strchr(SEQUENCETYPES, options->datatype))
    sprintf(starttree, "is supplied in %s", options->utreefilename);
  else {
    if (options->dist && strchr(SEQUENCETYPES, options->datatype))
      sprintf(starttree, "generates using %s", options->distfilename);
    else {
      if (options->randomtree)
	strcpy(starttree, "random start genealogy");
      else
	strcpy(starttree, "UPGMA based start genealogy");
    }
  }
}


///
/// checks whether the data type can use a start tree or not
/// and also allows to change the datatype
void            check_type_tree(char *input, option_fmt * options) {
  switch (options->datatype) 
    {
    case 'a':
      switch (atol(input)) {
      case DTEPTREE:
	start_tree_method(options);
	break;    
      case DTEPTYPE:
        start_data_method(options);
        break;
      default:
	break;
      }
      break;
    case 'm':
      switch (atol(input)) 
	{
	case DTMSATTREE:
	  start_tree_method(options);
	  break;
	case DTMSATTYPE:
	  start_data_method(options);
	  break;        
	default:
	  break;
	}
      break;
    case 'b':
      switch (atol(input)) 
	{
	case DTBROWNTREE:
	  start_tree_method(options);
	  break;
	case DTBROWNTYPE:
	  start_data_method(options);
	  break;
	default:
	  break;
	}
      break;
    case 's':
    case 'n':
    case 'h':
      switch (atol(input)) 
	{
	case DTSEQTREE:
	  start_tree_method(options);
	  break;
	case DTSEQTYPE:
	  start_data_method(options);
	  break;
	default:
	  break;
	}
      break;
      //case 'u':
      //case 'f':
    case 'g':
      switch (atol(input)) {
      case 1:
	start_data_method(options);
	break;
      default:
	break;
      }
      break;
  default:
    break;
  }
}


///
/// changing tipdate treatment and filename
void change_tipdate(char *input, option_fmt * options)
{
	printf(" Samples (Tips) where sampled at different dates? [No]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'Y') 
	  {
	    options->has_datefile = TRUE;
	    get_filename(" What is the filename that contains the sample dates?", TIPDATEFILE, options->datefilename);
	  } 
	else 
	  {
	    options->has_datefile = FALSE;
	  }
	input[0] = 'X';
}

///
/// changing mutation rate associated with tipdate
void change_mutationrate(char *input, option_fmt * options)
{
  char *in;
  char *tmp;
  long count=0;
  printf(" Give the number of loci or\n when all loci have the same mutation rate, the rate\nor a C to cancel\n===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (uppercase(input[0]) == 'C')
    return;
  if (input[0] == '0') 
    {
      options->mutationrate_year_numalloc = 1;
      options->mutationrate_year[0] = atof(input);
    } 
  else 
    {
      options->mutationrate_year_numalloc = atol(input);
      options->mutationrate_year = (MYREAL *) myrealloc(options->mutationrate_year,
							sizeof(MYREAL)*options->mutationrate_year_numalloc);
      while(count < options->mutationrate_year_numalloc)
	{
	  printf(" Give the mutation rate per year for each locus\n[either in groups or singly use a space to separate]\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  in = input;
	  tmp = strsep(&in,";, ");
	  while(tmp != NULL)
	    {
	      options->mutationrate_year[count++] = atof(tmp);
	      if(in==NULL)
		break;
	      tmp = strsep(&in,";, ");
	    }
	}
    }
  input[0] = 'X';
}

///
/// changing generation time associated with tipdates
void change_generationtime(char *input, option_fmt * options)
{
  printf(" Give the number of generation per year\nor a C to cancel\n===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (uppercase(input[0]) == 'C')
    return;
  options->generation_year = atof(input);
  input[0] = 'X';
}

///
/// changing hwo many individuals are used for the run, subsetting the data randomly
void change_randomsubset(char *input, option_fmt * options)
{
  printf(" How may individuals should be used\nGive a number or A for all or a C to cancel\n===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (uppercase(input[0]) == 'C')
    return;
  if (uppercase(input[0]) == 'A')
    {
      options->randomsubset = 0;
    }
  else
    {
      options->randomsubset = atof(input);
    }
  input[0] = 'X';
}

///
/// changing inheritance scalar for each locus
void change_inheritance(char *input, option_fmt * options)
{
  char *in;
  char *tmp;
  long count=0;
  printf(" Inheritance scalars: Give the number of loci or\nor a C to cancel\n===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (uppercase(input[0]) == 'C')
    return;
  
  options->inheritance_scalars_numalloc = atol(input);
  options->inheritance_scalars = (MYREAL *) myrealloc(options->inheritance_scalars,
						    sizeof(MYREAL)*options->inheritance_scalars_numalloc);
  while(count < options->inheritance_scalars_numalloc)
    {
      printf(" Give the inheritance scaler for each locus\n[either in groups or singly use a space to separate]\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      in = input;
      tmp = strsep(&in,";, ");
      while(tmp != NULL)
	{
	  options->inheritance_scalars[count++] = atof(tmp);
	  if(in==NULL)
	    break;
	  tmp = strsep(&in,";, ");
	}
    }
  input[0] = 'X';
}

///
///let the user set all data related options
void
menuData(option_fmt * options, char datatype[]) {
  static boolean  didchangecat, didchangercat;

  long            i;
  long            z = 0;
  long            numchar = 0;

  char            input[LINESIZE];
  char            starttree[LINESIZE];

  char           *ostr = (char *) mycalloc(LINESIZE, sizeof(char));

  didchangecat = setup_categs(options);
  didchangercat = setup_rcategs(options);

  do {
    setup_datatype(datatype, options);
    setup_starttree(starttree, options);
    printf("  DATATYPE AND DATA SPECIFIC OPTIONS\n\n");
    if(options->datatype=='m')
      printf(" %2i   change Datatype, currently set to: Stepwise mutation model (%10.10s)\n", DTSEQTYPE,submodeltype(options->msat_option));
    else
      printf(" %2i   change Datatype, currently set to:%29.29s\n", DTSEQTYPE, datatype);
    switch (options->datatype) {
      // allelic data --------------------------------------------------------------------
    case 'a':
      printf(" %2i   Discard missing data       %37.37s\n", DTEPMISS, options->include_unknown ? " NO" : "YES");
      printf(" %2i   Start genealogy            %37.37s\n", DTEPTREE, starttree);

#ifdef NEWVERSION
      printf(" %2i   Inheritance scalar set                                        %3.3s\n", DTEPINHERITANCE, 
	     options->inheritance_scalars_numalloc>1 ? "YES" : " NO");
#endif
      if(options->randomsubset > 0)
	printf(" %2i   Pick random subset per population of individuals            %4li\n", DTEPRANDOMSUBSET,  options->randomsubset);
      else
	printf(" %2i   Pick random subset per population of individuals               NO\n", DTEPRANDOMSUBSET);

      if(options->has_datefile)
	{
	  printf(" %2i   Tip date file %49.49s\n", DTEPTIPDATE, options->datefilename);
	  if(options->mutationrate_year_numalloc > 1)
	    sprintf(ostr,"multiple rates");
	  else
	    sprintf(ostr,"%.12f",options->mutationrate_year[0]);
	  printf(" %2i   Mutation rate per locus and year %30.30s\n", DTEPMUTRATE, ostr);
	  sprintf(ostr,"%10.4f",options->generation_year);
	  printf(" %2i   How many generations per year  %30.30s\n", DTEPGENERATION, ostr);
	}
      else
	{
	  printf(" %2i   Tip date file                      None, all tips a contemporary\n", DTEPTIPDATE);
	}
      break;
      // brownian motion data --------------------------------------------------------------------
    case 'b':
      printf(" %2i   Discard missing data:      %37.37s\n", DTBROWNMISS, options->include_unknown ? " NO" : "YES");
      printf(" %2i   Start genealogy            %37.37s\n", DTBROWNTREE, starttree);

#ifdef NEWVERSION
      printf(" %2i   Inheritance scalar set                                        %3.3s\n", DTBROWNINHERITANCE, 
	     options->inheritance_scalars_numalloc>1 ? "YES" : " NO");
#endif
      if(options->randomsubset > 0)
	printf(" %2i   Pick random subset per population of individuals             %4li\n", DTBROWNRANDOMSUBSET,  options->randomsubset);
      else
	printf(" %2i   Pick random subset per population of individuals              NO\n", DTBROWNRANDOMSUBSET);


      if(options->has_datefile)
	{
	  printf(" %2i   Tip date file              %37.37s\n", DTBROWNTIPDATE, options->datefilename);
	  if(options->mutationrate_year_numalloc > 1)
	    sprintf(ostr,"multiple rates");
	  else
	    sprintf(ostr,"%.12f",options->mutationrate_year[0]);
	  printf(" %2i   Mutation rate per locus and year %30.30s\n", DTBROWNMUTRATE, ostr);
	  sprintf(ostr,"%10.4f",options->generation_year);
	  printf(" %2i   Number of generations per year   %30.30s\n",DTBROWNGENERATION, ostr);
	}
      else
	{
	  printf(" %2i   Tip date file                      None, all tips a contemporary\n", DTBROWNTIPDATE);
	}
      break;
    case 'm':
      printf(" %2i   Discard missing data:     %37.37s\n", DTMSATMISS, options->include_unknown ? " NO" : "YES");
      printf(" %2i   Threshold value:                                     %10li\n",
	     DTMSATTRESH, options->micro_threshold);
      printf(" %2i   Start genealogy           %37.37s\n", DTMSATTREE, starttree);

#ifdef NEWVERSION
      printf(" %2i   Inheritance scalar set                                       %3.3s\n", DTMSATINHERITANCE, 
	     options->inheritance_scalars_numalloc>1 ? "YES" : " NO");
#endif
      if(options->randomsubset > 0)
	printf(" %2i   Pick random subset per population of individuals        %4.4li\n", DTMSATRANDOMSUBSET,  options->randomsubset);
      else
	printf(" %2i   Pick random subset per population of individuals              NO\n", DTMSATRANDOMSUBSET);


      if(options->has_datefile)
	{
	  printf(" %2i   Tip date file              %37.37s\n", DTMSATTIPDATE, options->datefilename);
	  if(options->mutationrate_year_numalloc > 1)
	    sprintf(ostr,"multiple rates");
	  else
	    sprintf(ostr,"%.12f",options->mutationrate_year[0]);
	  printf(" %2i   Mutation rate per locus and year %30.30s\n", DTMSATMUTRATE, ostr);
	  sprintf(ostr,"%10.4f",options->generation_year);
	  printf(" %2i   Number of generations per year   %30.30s\n", DTMSATGENERATION, ostr);
	}
      else
	{
	  printf(" %2i   Tip date file                      None, all tips a contemorary\n", DTMSATTIPDATE);
	}

      break;
    case 'u':
    case 'n':
    case 'h':
    case 's':
    case 'f':
      z = 0;
      numchar = 0;
      while (options->ttratio[z] > 0.0)
	numchar += sprintf(ostr + numchar, "%8.4f ", options->ttratio[z++]);
      printf(" %2i   Transition/transversion ratio:   %31s\n", DTSEQTRATIO, ostr);
      printf(" %2i   Use empirical base frequencies?  %30s\n", DTSEQFREQ, (options->freqsfrom ? "YES" : "NO"));
      printf(" %2i   Fixed categories for each site?  %30.30s\n", DTSEQSITECATEGS, options->categs == ONECATEG ?
	     "One category" :
	     "Categories supplied in file");
      printf(" %2i   Site rate variation?", DTSEQRATES);
      if (options->rcategs == 1)
	printf("                                        YES\n");
      else {
	printf("                 %4li categories of regions\n", options->rcategs);
	printf(" %2i   Rates at adjacent sites correlated?", DTSEQCORR);
	if (!options->autocorr)
	  printf("    NO, they are independent\n");
	else
	  printf("YES, mean block length =%6.1f\n",
		 1.0 / options->lambda);
      }
      printf(" %2i   Sites weighted?            %36.36s\n", DTSEQWEIGHT,
	     (options->weights ? "YES" : "NO"));
      printf(" %2i   Input sequences interleaved? %34.34s\n", DTSEQINTERLEAVED,
	     (options->interleaved ? "YES" : "NO, sequential"));
      printf(" %2i   Sequencing error rate? [0.0 = no error]                   %4.3f\n", DTSEQERROR,
	     options->seqerror);
      printf(" %2i   Slow but safer Data likelihood calculation %20.20s\n", DTSEQFAST,
	     options->fastlike ? "NO" : "YES");
      printf(" %2i   Start genealogy           %37.37s\n", DTSEQTREE, starttree);

#ifdef NEWVERSION
      printf(" %2i   Inheritance scalar set                                      %-3.3s\n", DTSEQINHERITANCE, 
	     options->inheritance_scalars_numalloc>1 ? "YES" : " NO");
#endif
      if(options->randomsubset > 0)
	printf(" %2i   Pick random subset per population of individuals           %4li\n", DTSEQRANDOMSUBSET,  options->randomsubset);
      else
	printf(" %2i   Pick random subset per population of individuals             NO\n", DTSEQRANDOMSUBSET);


      if(options->has_datefile)
	{
	  printf(" %2i   Tip date file             %37.37s\n", DTSEQTIPDATE, options->datefilename);
	  if(options->mutationrate_year_numalloc > 1)
	    sprintf(ostr,"multiple rates");
	  else
	    sprintf(ostr,"%.12f",options->mutationrate_year[0]);
	  printf(" %2i   Mutation rate per locus and year %30.30s\n", DTSEQMUTRATE, ostr);
	  sprintf(ostr,"%10.4f",options->generation_year);
	  printf(" %2i   Number of generations per year   %30.30s\n", DTSEQGENERATION, ostr);
	}
      else
	{
	  printf(" %2i   Tip date file                      None, all tips a contemorary\n", DTSEQTIPDATE);
	}
      break;

    case 'g':
      printf("       [Reanalyze an old run]\n");
      break;
    }
    printf("\n\n");
    printf("  Are the settings correct?\n");
    printf("  (Type Y or the number of the entry to change)\n===> ");
    input[0] = '\0';
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (strlen(input) == 0)
      continue;
    check_type_tree(input, options);
    //start tree option for all and switch datatype
    if (strchr("ab", options->datatype)) 
      // this only works if brownian and allele datatype have the same number of options
      {
	switch (atol(input)) 
	  {
	  case DTEPMISS:	/* discard unknowns */
	    printf("  Discard missing Alleles (YES|NO)?\n");
	    printf("  [Default is YES, this is the best choice for most situations]\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    if (strchr("Yy", input[0]))
	      options->include_unknown = FALSE;
	    else
	      options->include_unknown = TRUE;
	    input[0] = '\0';
	    break;
	  case DTEPINHERITANCE:
	    change_inheritance(input,options);
	    break;
	  case DTEPRANDOMSUBSET:
	    change_randomsubset(input,options);
	    break;
	  case DTEPTIPDATE: 
	    change_tipdate(input,options);
	    break;
	  case DTEPMUTRATE:
	    change_mutationrate(input,options);
	    break;
	  case DTEPGENERATION:
	    change_generationtime(input,options);
	    break;
	  default:
	    break;
	  }
      }
    if (options->datatype == 'm') 
      {
	switch (atol(input)) {
	case DTMSATMISS:	/* discard unknowns */
	  printf("  Discard missing Alleles?\n");
	  printf("  [Default is YES, this is the best choice for most situations]\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  if (strchr("Yy", input[0]))
	    options->include_unknown = TRUE;
	  else
	    options->include_unknown = FALSE;
	  break;
	case DTMSATTRESH:	/* micro-threshold */
	  printf("  What is the threshold value? [needs to be even!]\n");
	  printf
	    ("  E.g. if your allele is 24 and the threshold is 10\n");
	  printf("  there is some probability that the allele 24 can\n");
	  printf
	    ("  change to allele 14 (or 38), but there is a probability\n");
	  printf("  of 0.0 (ZERO) to go to 13 (39),\n");
	  printf
	    ("  if you choose this too small, than the program will fail\n");
	  printf
	    ("  The default is set 100, if the biggest difference in the data is smaller the value\n");
	  printf("  will be adjust to that maximal difference\n");
	  printf
	    ("  [the bigger the longer the run and the more\n accurate the estimate]\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  options->micro_threshold = atol(input);
	  if(options->micro_threshold % 2 != 0)
	    options->micro_threshold += 1;
	  input[0] = '\0';
	  break;
	case DTMSATTIPDATE: 
	  change_tipdate(input,options);
	  break;
	case DTMSATMUTRATE:
	  change_mutationrate(input,options);
	  break;
	case DTMSATGENERATION:
	  change_generationtime(input,options);
	  break;
	  case DTMSATINHERITANCE:
	    change_inheritance(input,options);
	    break;
	  case DTMSATRANDOMSUBSET:
	    change_randomsubset(input,options);
	    break;
	default:
	  break;
	}
      }
    if (strchr(SEQUENCETYPES, options->datatype)) {
      switch (atol(input)) {
      case DTSEQTRATIO:
	initratio(options);
	break;
      case DTSEQFREQ:
	options->freqsfrom = !options->freqsfrom;
	if (!options->freqsfrom) {
	  initfreqs(&options->freqa, &options->freqc,
		    &options->freqg, &options->freqt);
	}
	break;
      case DTSEQSITECATEGS:
	if (options->categs == ONECATEG) {
	  options->categs = MANYCATEGS;
	  didchangecat = TRUE;
	} else {
	  options->categs = ONECATEG;
	  didchangecat = FALSE;
	}
	break;
      case DTSEQRATES:
	if (didchangercat) {
	  options->autocorr = FALSE;
	  didchangercat = FALSE;
	  options->rcategs = ONECATEG;
	} else {
	  printf("\n  Regional rates:\n");
	  initcatn(&options->rcategs);
	  options->probcat =
	    (MYREAL *) myrealloc(options->probcat,
				 options->rcategs * sizeof(MYREAL));
	  options->rrate =
	    (MYREAL *) myrealloc(options->rrate,
				 options->rcategs * sizeof(MYREAL));
	  didchangercat = TRUE;
	  options->gammarates =
	    initcategs(options->rcategs, options->rrate,
		       options->probcat);
	  if (!options->gammarates)
	    {
	      initprobcat(options->rcategs, &options->probsum,
			  options->probcat);
	      constrain_rates(options->rcategs, options->rrate,
			      options->probcat);
	    }
	    FPRINTF(stdout,
		    "\n\nRegion type     Rate of change    Probability\n");
	  FPRINTF(stdout,
		  "---------------------------------------------\n");
	  for (i = 0; i < options->rcategs; i++)
	    FPRINTF(stdout, "%9ld%16.3f%17.3f\n",
		    i + 1, options->rrate[i], options->probcat[i]);
	  FPRINTF(stdout, "\n");
	}
	break;
      case DTSEQCORR:
	options->autocorr = !options->autocorr;
	if (options->autocorr)
	  initlambda(options);
	break;
      case DTSEQWEIGHT:
	options->weights = !options->weights;
	break;
      case DTSEQINTERLEAVED:
	options->interleaved = !options->interleaved;
	break;
      case DTSEQERROR:
	FPRINTF(stdout,
		"Enter the sequencing error rate per site\n[Good values are 0.0 (=no error), or 0.01 (1%% error):\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	options->seqerror = atof(input);
	input[0] = '\0';
	break;
      case DTSEQFAST:
	options->fastlike = !options->fastlike;
	break;
      case DTSEQTIPDATE: 
	change_tipdate(input,options);
	break;
      case DTSEQMUTRATE:
	change_mutationrate(input,options);
	break;
      case DTSEQGENERATION:
	change_generationtime(input,options);
	break;
      case DTSEQINHERITANCE:
	change_inheritance(input,options);
	break;
      case DTSEQRANDOMSUBSET:
	change_randomsubset(input,options);
	break;
      default:
	break;
      }
      if (!didchangercat) {
	options->probcat =
	  (MYREAL *) myrealloc(options->probcat, sizeof(MYREAL) * 2);
	options->rrate =
	  (MYREAL *) myrealloc(options->rrate, sizeof(MYREAL) * 2);
	options->rrate[0] = 1.0;
	options->probcat[0] = 1.0;
      }
      if (!didchangecat) {
	options->rate =
	  (MYREAL *) myrealloc(options->rate, sizeof(MYREAL) * 2);
	options->rate[0] = 1.0;
      }
    }
    if (options->datatype != 'g')
      options->readsum = FALSE;
  }
  while (uppercase(input[0]) != 'Y');

  myfree(ostr);

}

void
menuInput(option_fmt * options) {
  long            test = 0;
  long            len = (long) strlen(options->title);
  unsigned long   timeseed;
  char            input[LINESIZE];
  char            treeinc[LINESIZE];
  char            outputstring[LINESIZE];
  char           *stringstep, *string2, *string3, *string4;
  char            treestr[4][20] = {"None", "All", "Best", "Last chain"};
  char            progress[3][8] = {"Verbose", "YES", "NO"};
  char           *extension;
  stringstep = (char *) mymalloc(sizeof(char) * 128);
  string2 = (char *) mymalloc(sizeof(char) * 128);
  string3 = (char *) mymalloc(sizeof(char) * 128);
  string4 = (char *) mymalloc(sizeof(char) * 128);

  do {
    printf("  INPUT/OUTPUT FORMATS (for %s)\n  -------------------- [approach can be changed in SEARCH Strategy]\n\n  INPUT:\n",options->bayes_infer ? "Bayesian approach" : "Maximum likelihood approach");
    printf("  %2i   Datafile name is %46.46s\n", MI_INFILE,
	   options->infilename);

    sprintf(treeinc,"[every %li]",options->treeinc);

    switch (options->autoseed) {
    case AUTO:
      sprintf(stringstep, "YES");
      break;
    case NOAUTO:
      sprintf(stringstep, "NO, use seedfile");
      break;
    case NOAUTOSELF:
      sprintf(stringstep, "NO, seed=%li ", options->inseed);
      break;
    default:
      options->autoseed = AUTO;
      sprintf(stringstep, "YES");
      break;
    }

    printf("  %2i   Use automatic seed for randomisation?  %24.24s\n",
	   MI_RAND, stringstep);
    printf("  %2i   Title of the analysis is", MI_TITLE);

    if (len == 0 || (len == 5 && strstr(options->title, "Migration analysis")))
      printf("%39.39s\n", "<no title given>");
    else
      printf("%39.39s\n", options->title);

    if (options->readsum && !options->bayes_infer) {
      printf("  %2i    Summary of genealogies are read from %s\n",
	     MI_SUMREAD, options->sumfilename);
    }
    else
      {
	if (options->readsum && options->bayes_infer) {
	  printf("  %2i    Bayesian output data are read from %s\n",
		 MI_SUMREAD, options->bayesmdimfilename);
	}
      }
    printf("\n  OUTPUT:\n");

    printf("  %2i   Print indications of progress of run? %25.25s\n",
	   MI_PROGRESS, options->progress ? (options->
					     verbose ? progress[0] :
					     progress[1]) : progress[2]);
    printf("  %2i   Print the data?%48.48s\n",
	   MI_PRINTDATA, options->printdata ? "YES" : "NO");

    printf("  %2i   Outputfile name is  %43.43s\n", MI_OUTFILE,
	   options->outfilename);
#ifdef PRETTY
    printf("                           %43.43s\n", options->pdfoutfilename);
#endif
    if(!options->bayes_infer)
      {
	if (options->plot) 
	  {
	    switch (options->plotmethod) 
	      {
	      case PLOTALL:
		strcpy(string2, "YES, to outfile and mathfile");
		break;
	      default:
		strcpy(string2, "YES, to outfile");
		break;
	      }
	  }
	else 
	  {
	    strcpy(string2, "NO");
	  }
	printf("  %2i   Plot likelihood surface? %38.38s\n",
	       MI_PLOT, string2);
	switch (options->profile) 
	  {
	  case _NONE:
	    strcpy(string3, "NO");
	    break;
	  case ALL:
	    strcpy(string3, "YES, tables and summary");
	    break;
	  case TABLES:
	    strcpy(string3, "YES, tables");
	    break;
	  case SUMMARY:
	    strcpy(string3, "YES, summary");
	    break;
	  }
	printf("  %2i   Profile-likelihood?%44.44s\n", MI_PROFILE, string3);
	if (options->profile != _NONE) 
	  {
	    switch (options->profilemethod) 
	      {
	      case 'p':
		strcpy(string4, "[Percentiles using exact Bisection method]");
		break;
	      case 's':
		strcpy(string4,
		       "[Percentiles using not so exact Spline method]");
		break;
	      case 'd':
		strcpy(string4,
		       "[Evaluation at preset fractions of ML estimate]");
		break;
	      case 'q':
		strcpy(string4, "[Assuming that parameter are independent]");
		break;
	      case 'f':
		strcpy(string4, "[Mixture of 'p' and 'q']");
		break;
	      }
	    printf("%70.70s\n", string4);
	  }
	printf("  %2i   Likelihood-Ratio tests?%40.40s\n",
	       MI_LRT, options->lratio->counter > 0 ? " YES" : "NO");
	
	printf("  %2i   AIC model selection?   %40.40s\n",
	       MI_AIC, options->aic ? (options->fast_aic ? "YES:Fast" : "YES") : "NO");
      } /*if !bayes_infer*/

    switch (options->treeprint) 
      {
      case ALL:
	printf("  %2i   Print genealogies?     %40.40s\n",
	       MI_TREES, treestr[1]);
	break;
      case BEST:
	printf("  %2i   Print genealogies?     %40.40s\n",
	       MI_TREES, treestr[2]);
	break;
      case LASTCHAIN:
	printf("  %2i   Print genealogies?  %19.19s %20.20s\n",
	       MI_TREES, treeinc ,treestr[3]);
	break;
      case _NONE:
      default:
	printf("  %2i   Print genealogies?     %40.40s\n",
	       MI_TREES, treestr[0]);
	break;
      }
    if(!options->bayes_infer)
      {
	printf("  %2i   Plot coordinates are saved in %33.33s\n",
	       MI_PLOTCOORD, options->mathfilename);
	if (options->writesum || options->datatype == 'g') 
	  {
	    if (!options->readsum && options->writesum) 
	      {		
		printf("  %2i   Summary of genealogies saved in %31.31s\n",
		       MI_SUMFILE, options->sumfilename);
	      }
	  }
	else 
	  {
	    printf("  %2i   Summary of genealogies  %39.39s\n",
		   MI_SUMFILE, "will not be saved");
	  }
      }
    if (options->writelog)
      printf("  %2i   Save logging information into %33.33s\n", MI_LOGFILE,
	     options->logfilename);
    else
      printf("  %2i   Save logging information?  %36.36s\n", MI_LOGFILE, "NO");
#ifdef UEP

    if (options->uep) {
      printf("  %2i   Read unique polymorphism? From file %30.30s\n",
	     MI_UEPFILE, options->uepfilename);
      printf("  %2i   Mutation rate for UEP is %10.5f %s\n",
	     MI_UEPRATE, options->ueprate, " x mutation rate");
      printf("  %2i   Base frequency for UEP is \"0\"=%f, \"1\"=%f\n",
	     MI_UEPFREQ, options->uepfreq0, options->uepfreq1);
    } else
      printf("  %2i   Read unique polymorphism? %20.20s\n", MI_UEPFILE, "NO");
#endif
    if (options->mighist) 
      {
	sprintf(outputstring,"%s %s", options->mighistfilename,
		options->mighist_all ? "(all events)" : "(migration events)");
	printf("  %2i   Show event statistics%42.42s\n", MI_MIGHISTOGRAM, outputstring);
	if(options->mighist_increment > 1)
	  sprintf(outputstring,"every %li %s",  options->mighist_increment,"sample steps");
	else
	  sprintf(outputstring,"every %s", "sample step");
	printf("       Events are recorded every     %32.32s\n",outputstring);
	sprintf(outputstring,"%f",options->eventbinsize);
	printf("       Histogram bin width            %32.32s\n",outputstring);

      } 
    else 
      {
	sprintf(outputstring,"NO");
	printf("  %2i   Show event statistics          %32.32s\n", MI_MIGHISTOGRAM, outputstring);
    
      }
    if (options->skyline) 
      {
	remove_trailing_blanks(&options->skylinefilename);
	sprintf(outputstring,"%s", options->skylinefilename);
	printf("  %2i   Record parameter change through time?%26.26s\n", MI_SKYLINE, outputstring);
	sprintf(outputstring,"%f",options->eventbinsize);
	printf("       Histogram bin width            %32.32s\n",outputstring);
      } 
    else 
      {
	sprintf(outputstring,"NO");
	printf("  %2i   Record parameter change through time?%26.26s\n", MI_SKYLINE, outputstring);
      }

    printf("\n\n  Are the settings correct?\n");
    printf
      ("  (type Y to go back to the main menu or the letter for the entry to change)\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    test = atoi(input);
    switch (test) 
      {
      case MI_INFILE:
	printf("  What is the datafile name?\n[Default: %s]\n===> ", options->infilename);
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (input[0] != '\0')
	  {
	    strcpy(options->infilename, input);
	  }
	break;
      case MI_RAND:
	do 
	  {
	    printf("  (A)utomatic or (S)eedfile or (O)wn\n");
	    printf("  Start value for Random-generator seed\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    switch (uppercase(input[0])) 
	      {
	      case 'A':
		options->autoseed = AUTO;
#ifndef MAC
		timeseed = (unsigned long) time(NULL) / 4;
#else
		timeseed = (unsigned long) clock() / 4;
#endif
		options->inseed = (long) timeseed + 1;
		break;
	      case 'S':
		get_filename(" What is the filename that contains the random number seed?", SEEDFILE, options->seedfilename);
		openfile(&options->seedfile, options->seedfilename, "r", NULL);
		if (options->seedfile) {
		  options->autoseed = NOAUTO;
		  fscanf(options->seedfile, "%ld%*[^\n]",
			 &options->inseed);
		  fclose(options->seedfile);
		} else
		  printf("\n\n  There is no seedfile present\n");
		break;
	      case 'O':
		options->autoseed = NOAUTOSELF;
		printf("  Random number seed (best values are x/4 +1)?\n");
		scanf("%ld%*[^\n]", &options->inseed);
		break;
	      }
	  }
	while (options->autoseed < AUTO || options->autoseed > NOAUTOSELF);
	break;
      case MI_TITLE:
	printf("  Enter a title? [max 80 Characters are used]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (input[0] == '\0')
	  options->title[0] = '\0';
	else
	  strncpy(options->title, input, 80);
	break;
      case MI_SUMREAD:
	if(!options->bayes_infer)
	  {
	    printf
	      (" What is the filename for the summary of genealogies\n[Default: %s]\n===> ",
	       SUMFILE);
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    if (input[0] == '\0')
	      strcpy(options->sumfilename, SUMFILE);
	    else {
	      strcpy(options->sumfilename, input);
	    }
	  }
	else
	  {
	    printf
	      (" What is the filename of the recorded Bayes data\n[Default: %s]\n===> ",
	       BAYESMDIMFILE);
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    if (input[0] == '\0')
	      strcpy(options->bayesmdimfilename, BAYESMDIMFILE);
	    else {
	      strcpy(options->bayesmdimfilename, input);
	    }
	    unpad(options->bayesmdimfilename, " ");
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
      case MI_PROGRESS:
	printf("  Progress report during the run? <YES | Verbose | NO>\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	switch (tolower(input[0])) {
	case 'n':
	  options->progress = FALSE;
	  options->verbose = FALSE;
	  break;
	case 'v':
	  options->progress = TRUE;
	  options->verbose = TRUE;
	  break;
	case '\0':
	case 'y':
	default:
	  options->progress = TRUE;
	  options->verbose = FALSE;
	  break;
	}
	input[0] = 'X';
	break;
      case MI_PRINTDATA:
	options->printdata = !options->printdata;
	break;
      case MI_OUTFILE:
	get_filename("  What is the output filename?", OUTFILE, options->outfilename);
#ifdef PRETTY
	strcpy(options->pdfoutfilename,options->outfilename);
	strcat(options->pdfoutfilename,".pdf");
#endif
	break;
      case MI_PLOT:
	options->plot = !options->plot;
	if (options->plot) 
	  {
	    get_plotoptions(options);
	    input[0] = '\0';
	  }
	break;
      case MI_PROFILE:
	do 
	  {
	    printf("  Evaluate profile likelihoods:\n");
	    printf
	      ("  (N)o, (A)all [Tables and Summary],\n  (T)ables, (S)ummary\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	  }
	while (strchr("NATS", (int) uppercase(input[0])) == NULL);
	switch (uppercase(input[0])) 
	  {
	  case 'N':
	    options->profile = _NONE;
	    break;
	  case 'A':
	    options->profile = ALL;
	    break;
	  case 'T':
	    options->profile = TABLES;
	    break;
	  case 'S':
	    options->profile = SUMMARY;
	    break;
	  }
	if (options->profile != _NONE) 
	  {
	do 
	  {
	    printf("  Method to evaluate the profiles?\n");
	    printf("  (P)ercentiles: exact evaluation (slow)\n");
	    printf
	      ("  (F)ast       : Assuming parameters are indepenent\n                + one complete profile-maximization.\n");	    
	    //printf
	    // ("  (S)plines    : percentiles using splines\n                (fails sometimes)\n");
	    printf
	      ("  (D)iscrete   : evaluation at (0.02, 0.05, 0.1, .. , 50)*ML estimate\n");
	    printf
	      ("  (Q)uick      : Assuming parameters are independent\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	  }
	while (strchr("pdqf", (int) lowercase(input[0])) == NULL);
	switch (lowercase(input[0]))
	  //in strchr removed s
	  {
	  case 'p':
	    options->profilemethod = 'p';
	    break;
	  case 's':
	    options->profilemethod = 's';
	    break;
	  case 'd':
	    options->profilemethod = 'd';
	    break;
	  case 'q':
	    options->profilemethod = 'q';
	    break;
	  case 'f':
	    options->profilemethod = 'f';
	    break;
	  default:
	    options->profilemethod = 's';
	    break;
	  }
      }
	set_profile_options(options);
	break;
      case MI_LRT:
	if (options->lratio->counter > 0) 
	  {
	    options->lratio->counter = 0;
	  } 
	else 
	  {
	    FPRINTF(stdout, "  Likelihood-Ratio tests:\n");
	    FPRINTF(stdout, "  -----------------------\n");
	    read_custom_menu_lratio(options);
	  }
	break;
      case MI_AIC:
	printf
	  ("  Use AIC to select minimal migration model?\n[Default: NO]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	switch (uppercase(input[0])) 
	  {
	  case 'F':
	    options->fast_aic = TRUE;
	    options->aic = TRUE;
	    break;
	  case 'Y':
	    options->fast_aic = FALSE;
	    options->aic = TRUE;
	    break;
	  case 'N':
	  default:
	    options->aic = FALSE;
	    options->fast_aic = FALSE;
	    break;
	  }
	input[0] = 'x';
	break;
      case MI_TREES:
	do 
	  {
	    printf("  Print genealogies:\n");
	    printf("  (N)one, (A)all [!], (B)est, (L)ast chain\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	  }
	while (strchr("NABL", (int) uppercase(input[0])) == NULL);
	switch (uppercase(input[0])) 
	  {
	  case 'N':
	    options->treeprint = _NONE;
	    break;
	  case 'A':
	    options->treeprint = ALL;
	    get_filename(" What is the filename for storing recorded genealogies?", TREEFILE, options->treefilename);
	    break;
	  case 'B':
	    options->treeprint = BEST;
	    get_filename(" What is the filename for storing recorded genealogies?", TREEFILE, options->treefilename);
	    break;
	  case 'L':
	    options->treeprint = LASTCHAIN;
	    get_filename(" What is the filename for storing recorded genealogies?", TREEFILE, options->treefilename);
	    printf("Give the increment to record genealogies\n===>");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    sscanf(input,"%li",&options->treeinc);
	    break;
	  default:
	    options->treeprint = _NONE;
	    break;
	  }
	break;
      case MI_SUMFILE:
	printf(" Save genealogy summaries? [NO]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'Y') 
	  {
	    options->writesum = TRUE;
	    get_filename(" What is the filename for the summary of genealogies?", SUMFILE, options->sumfilename);
	  } 
	else 
	  {
	    options->writesum = FALSE;
	  }
	input[0] = 'X';
	break;
      case MI_PLOTCOORD:
	get_filename("  What is the plot coordinate filename?", MATHFILE, options->mathfilename);
	break;
      case MI_LOGFILE:
	printf(" Save logging information? [NO]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'Y') 
	  {
	    options->writelog = TRUE;
	    get_filename(" What is the filename for logging?", LOGFILE, options->logfilename);
	  }
	else 
	  {
	    options->writelog = FALSE;
	  }
	input[0] = 'X';
	break;
#ifdef UEP
      case MI_UEPFILE:
	printf(" UEP estimation? [NO]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'Y') 
	  {
	    options->uep = TRUE;
	    get_filename(" What is the filename of the UEP data?", UEPFILE, options->uepfilename);
	  } 
	else 
	  {
	    options->uep = FALSE;
	  }
      input[0] = 'X';
      break;
      case MI_UEPRATE:
      do 
	{
	  FPRINTF(stdout, "What is mutation rate for the UEP locus\n");
	  FPRINTF(stdout,
		  "[expressed as a ratio of the point mutation rate\n");
	  FPRINTF(stdout, "Reasonable values are 0.0 =< x <=1.0]\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  options->ueprate = atof(input);
	}
      while (options->ueprate < 0.0);
      break;
      case MI_UEPFREQ:
      do 
	{
	  FPRINTF(stdout, "What is the base frequency for the \"1\" allele\n");
	  FPRINTF(stdout, "Reasonable values are 0.0 < x <1.0]\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  options->uepfreq1 = atof(input);
	  options->uepfreq0 = 1. - options->uepfreq1;
	}
      while (options->uepfreq1 < 0.0);
      break;
#endif
      case MI_MIGHISTOGRAM:
	printf(" Display event statistic and save events into file?\n [options are: NO, ALL events, MIGration events]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'N') 
	  {
	    options->mighist = FALSE;
	  } 
	else 
	  {
	    options->mighist = TRUE;
	    if (uppercase(input[0]) == 'A')
	      options->mighist_all = TRUE;
	    else
	      options->mighist_all = FALSE;

	    /*	    printf(" Thin the raw data and sample only every x sample steps,\n Enter increment [default is 1]:\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    options->mighist_increment = atol(input);
	    if(options->mighist_increment < 1)
	      options->mighist_increment = 1;
	    */
	    get_filename(" Specify a filename for the raw event data!", MIGHISTFILE, options->mighistfilename);
	    
	    do
	      {
		printf("Events are recorded per time intervals.\nWhat is the time interval size?\n");
		printf("[Current value: %f, try values that are about 10x smaller than the largest Theta]\n===> ",
		       options->eventbinsize);
		fflush(stdout); FGETS(input, LINESIZE, stdin);
		if(input[0]=='\0' && options->eventbinsize > 0.0)
		  break;
		options->eventbinsize = atof(input);
	      } 
	    while (options->eventbinsize <= 0.0);

	  }
	input[0] = 'X';
	break;
      case MI_SKYLINE:
	printf(" Display time dependent parameters and save them into file? [Choices are NO or  YES]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'N') 
	  {
	    options->skyline = FALSE;
	  } 
	else 
	  {
	    options->skyline = TRUE;
	    do
	      {
		printf("What is the time interval size?\n");
		printf("[Current value: %f, try values that are about 10x smaller than the largest Theta]\n===> ",
		       options->eventbinsize);
		fflush(stdout); FGETS(input, LINESIZE, stdin);
		if(input[0]=='\0' && options->eventbinsize > 0.0)
		  break;
		options->eventbinsize = atof(input);
	      } 
	    while (options->eventbinsize <= 0.0);
	    get_filename(" Specify a filename for the raw skyline output!", SKYLINEFILE, options->skylinefilename);
	    if(!options->mighist)
	      {
		printf("The skyline option needs the coalescence and migration event recording\n");
		get_filename(" Specify a filename for the raw event data!", MIGHISTFILE, options->mighistfilename);		
		options->mighist_all = TRUE;
		options->mighist_increment = 0;
	      }

	  }
	input[0] = 'X';
	break;
      default:
	break;
      }
  }
  while (uppercase(input[0]) != 'Y');
  myfree(stringstep);
  myfree(string2);
  myfree(string3);
  myfree(string4);
}

///sets the start parameter options for theta
void            switch_theta(char *ostr, option_fmt * options) {
  long            i;
  long            numchar = 0;
  switch          (options->thetaguess) {
  case FST:
    sprintf(ostr, "Estimate with FST (Fw/Fb) measure");
    break;
  case URANDOMESTIMATE:
    sprintf(ostr, "Random Theta from Uniform distribution (min=%f, max=%f)", options->thetag[0], options->thetag[1]);
    break;
  case NRANDOMESTIMATE:
    sprintf(ostr, "Random Theta from Normal distribution (mean=%f, std=%f)", options->thetag[0], options->thetag[1]);
    break;
  default:
    numchar += sprintf(ostr, "NO, initial Theta = {");
    for (i = 0; i < options->numthetag - 1; i++) {
      numchar += sprintf(ostr + numchar, "%.5f,", options->thetag[i]);
    }
    sprintf(ostr + numchar, "%.5f}", options->thetag[i]);
  }
}

void            switch_mig(char *ostr, option_fmt * options) {
  switch (options->migrguess) {
  case FST:
    sprintf(ostr, "Estimate with FST (Fw/Fb) measure");
    break;
  case URANDOMESTIMATE:
    sprintf(ostr, "Random migration rates from Uniform distribution (min=%f, max=%f)", options->mg[0], options->mg[1]);
    break;
  case NRANDOMESTIMATE:
    sprintf(ostr, "Random migration rates from Normal distribution (mean=%f, std=%f)", options->mg[0], options->mg[1]);
    break;
  default:
    sprintf(ostr, "NO, initial migration rates are given");
  }
}


///
///allows to set all parameter related issues:start parameters, mutation type, migration model
void
menuParameters(option_fmt * options) 
{
#ifdef POPMODEL
  char            tmp[LINESIZE];
#endif
  char            input[LINESIZE];
  char            custmexplain[LINESIZE];
  char            localities[LINESIZE];
  char           *temp;
  char           *outputstring;
  long            count = 0;
  long            i, j, tt, check = 0, numpop = 0;
  MYREAL          sum = 0.;
  double          tmpdouble;

  outputstring = (char *) mycalloc(LINESIZE, sizeof(char));


  do {
    printf
      ("  PARAMETERS\n  ---------------------------\n  Start parameters:\n");
    switch_theta(outputstring, options);
    printf("  1   Use a simple estimate of theta as start?\n      %64s\n", outputstring);
    switch_mig(outputstring, options);
    printf("  2   Use a simple estimate of migration rate as start?\n      %64s\n", outputstring);
    printf
      ("\n Gene flow parameter and Mutation rate variation among loci:\n");
    printf("  3   Use M for the gene flow parameter   %28.28s\n",
	   options->usem ? "YES [M=m/mu]" : "NO [Theta * M]");
    if(options->gamma)
      {
	printf("  4   Mutation rate is  %46.46s%f\n","Estimated: Gamma distribution with alpha=",(float) options->alphavalue);
      }
    else
      {
	if(options->murates)
	  {
	    if(options->murates_fromdata)
	      {
		printf("  4   Mutation rate is  %46.46s\n","Varying (from data)");
	      }
	    else
	      {
		printf("  4   Mutation rate is  %46.46s\n","Varying (from user)");
	      }
	  }
	else
	  {
	    if(options->bayesmurates)
	      {
		printf("  4   Mutation rate is  %46.46s\n","Estimated (see prior menu: Rate)");
	      }
	    else
	      {
		printf("  4   Mutation rate is  %46.46s\n","Constant");
	      }
	  }
      }
    if (options->migrguess == FST || options->thetaguess == FST) {
      printf("\n  FST-Calculation (for START value):\n");
      printf("  5   Method: %56.56s\n",
	     options->fsttype ==
	     'T' ? "Variable Theta, M symmetric" :
	     "Variable M, Theta is the same for all populations");
      printf("  6   Print FST table:%48.48s\n", options->printfst ? "YES" : "NO");
    }
    strcpy(custmexplain, custom_migration_type(options->migration_model));
    set_localities_string(&localities[0],options);
    printf("\n  Migration model and combination of localities:\n");
    printf("  7   Sampling localities %43.43s\n", localities);
    printf("  8   Model is set to %48.48s\n", custmexplain);
    if(options->geo)
      printf("  9   Geographic distance matrix: YES:%34.34s\n", options->geofilename);
    else
      printf("  9   Geographic distance matrix: %36.36s\n", "NO");
#ifdef POPMODEL
    printf(" 10   Population model: %46.46s\n",
	   which_popmodel(options,tmp));
#endif /*POPMODEL*/
    printf("\n\n  Are the settings correct?\n");
    printf
      ("  (Type Y to go back to the main menu or the letter for an entry to change)\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    switch (input[0]) {
    case '1':
      printf("  Which method? (F)st or (O)wn or (N)ormally or ((U)niform distributed random value\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      switch (uppercase(input[0])) {
      case 'N':
	options->thetaguess = NRANDOMESTIMATE;
	options->thetag = (MYREAL *) myrealloc(options->thetag,
					       sizeof(MYREAL) * 2);
	menuRandom(options->thetag,'N');
	options->numthetag = 2;
	break;
      case 'U':
	options->thetaguess = URANDOMESTIMATE;
	options->thetag = (MYREAL *) myrealloc(options->thetag,
					       sizeof(MYREAL) * 2);
	menuRandom(options->thetag,'U');
	options->numthetag = 2;
	break;
      case 'F':
	options->thetaguess = FST;
	break;
      case 'O':
	options->thetaguess = OWN;
	if (numpop == 0) {
	  printf("  I do not know yet how many populations the data set contains? Specify the number of populations > ");fflush(stdout);
	  scanf("%ld", &numpop);
	  scanf("%*[^\n]");
	}
	options->thetag =
	  (MYREAL *) myrealloc(options->thetag,
			       sizeof(MYREAL) * (numpop + 1));
	printf("  Give the initial Theta estimates:\n");
	for (i = 0; i < numpop; i++) {
	  printf("Population %3li> ", i);fflush(stdout);
#ifdef USE_MYREAL_FLOAT
	  scanf("%f", &options->thetag[i]);
#else
	  scanf("%lf", &options->thetag[i]);
#endif
	}
	options->numthetag = i;
	break;
      default:
	options->thetaguess = FST;
	break;
      }

      break;
    case '2':
      printf("Which method?  (F)st or (O)wn or (N)ormally or (U)niformly distributed random value\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      switch (uppercase(input[0])) {
      case 'U':
	options->migrguess = URANDOMESTIMATE;
	options->mg = (MYREAL *) myrealloc(options->mg,
					   sizeof(MYREAL) * 2);
	menuRandom(options->mg,'U');
	options->nummg = 2;
	break;
      case 'N':
	options->migrguess = NRANDOMESTIMATE;
	options->mg = (MYREAL *) myrealloc(options->mg,
					   sizeof(MYREAL) * 2);
	menuRandom(options->mg,'N');
	options->nummg = 2;
	break;
      case 'F':
	options->migrguess = FST;
	break;
      case 'O':
	options->migrguess = OWN;
	if (numpop == 0) {
	  printf("  I do not know yet how many populations the data set contains? Specify the number of populations > ");
	  fflush(stdout);
	  scanf("%ld", &numpop);
	  scanf("%*[^\n]");
	}
	printf
	  ("  Initial migration rate estimate?\n[give 4Nm for diploid data\n 2Nm for haploid data\n Nm for mtDNA data]\n[If you use M (m/mu) instead of Nm then you need to specify M here !]");
	tt = 0;
	options->mg =
	  (MYREAL *) myrealloc(options->mg,
			       sizeof(MYREAL) * (numpop * numpop + 1));
	for (i = 0; i < numpop; i++) {
	  for (j = 0; j < numpop; j++) {
	    if (j == i) {
	      printf("Population %3li              > --\n",
		     i + 1);
	      continue;
	    }
	    printf("From population %-3li to %-3li> ", j + 1,
		   i + 1);fflush(stdout);
#ifdef USE_MYREAL_FLOAT
	    check = scanf("%f", &options->mg[tt++]);
#else
	    check = scanf("%lf", &options->mg[tt++]);
#endif
	    scanf("%*[^\n]");
	    if (check == 0)
	      break;
	  }
	}
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	options->nummg = tt;
	tt = 0;
	printf
	  ("\n     You typed in the following migration matrix\n\n     ");
	for (i = 0; i < numpop; i++) {
	  for (j = 0; j < numpop; j++) {
	    if (i != j)
	      printf("%5.2f ", options->mg[tt++]);
	    else {
	      printf("----- ");
	    }
	  }
	  printf("\n     ");
	}
	printf("\n     [Press <Return> to continue]\n");
	printf("\n\n");
	getchar();
	break;
      default:
	options->migrguess = FST;
	break;
      }
      break;
    case '3':
      options->usem = !options->usem;
      if (options->usem) {
	options->plotvar = PLOTM;
	options->migvar = PLOTM;
	options->profileparamtype = PLOTM;
      } else {
	options->plotvar = PLOT4NM;
	options->migvar = PLOT4NM;
	options->profileparamtype = PLOT4NM;
      }

      break;
    case '4':
      printf("Mutation rate among loci\n");
      printf
	("(C)onstant   All loci have the same mutation rate [default]\n");
      printf
	("(E)stimate  Mutation rate \n");
      if (options->gamma)
	{
	  printf("       GAMMA with Shape parameter alpha=%f [only ML]\n",
		 options->alphavalue);
	}
      printf("(V)arying     Mutation rates are different among loci [user input]\n");
      printf("(R)elative    Mutation rates estimated from data\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      switch (uppercase(input[0])) {
      case 'E':
	if(!options->bayes_infer)
	  {
	    options->gamma = TRUE;
	    printf("Enter the value of the shape parameter alpha to use\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    sscanf(input,"%lf",&tmpdouble);
	    options->alphavalue = (MYREAL) tmpdouble;
	    options->bayesmurates = FALSE;
	  }
	else
	  {
	    options->bayesmurates = TRUE;
	  }
	options->murates = FALSE;
	options->murates_fromdata=FALSE;

	break;
      case 'V':
	printf
	  ("Enter the number of loci and a rate for each of them\n[the rates are normalized to average to 1.0]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	temp = strtok(input, " ;,;\n");
	if (temp == NULL) {
	  options->murates = FALSE;
	  options->gamma = FALSE;
	  options->murates_fromdata=FALSE;
	  options->bayesmurates = FALSE;
	} else {
	  options->gamma = FALSE;
	  options->murates = TRUE;
	  options->murates_fromdata=FALSE;
	  options->muloci = atol(temp);
	  options->bayesmurates = FALSE;
	  options->mu_rates =
	    (MYREAL *) myrealloc(options->mu_rates,
				 sizeof(MYREAL) * options->muloci);
	  //  printf("%i> menu.c: opt realloc murate size %li\n",myID,  options->muloci * sizeof (MYREAL));			
	  sum = 0.;
	  count = 0;
	  while (temp != NULL) {
	    temp = strtok(NULL, " ;,;\n");
	    if (temp == NULL)
	      break;
	    options->mu_rates[count] = atof(temp);
	    sum += options->mu_rates[count];
	    count++;
	  }
	  for (i = count; i < options->muloci; i++) {
	    options->mu_rates[i] = options->mu_rates[count - 1];
	    sum += options->mu_rates[count];
	  }
	  sum /= options->muloci;
	  for (i = 0; i < options->muloci; i++) {
	    options->mu_rates[i] /= sum;
	  }
	}
	break;
      case 'R':
	options->gamma = FALSE;
	options->murates = TRUE;
	options->murates_fromdata=TRUE;
	options->bayesmurates = FALSE;
	break;
      default:
	options->gamma = FALSE;
	options->murates = FALSE;
	options->murates_fromdata=FALSE;
	options->bayesmurates = FALSE;
	break;
      }
      break;
    case '5':
      printf("Which FST calculation method?\n");
      printf("(T)heta can be different for each population\n");
      printf("   and migration rates are symmetric.\n");
      printf("   (Number of populations >= 2)\n");
      printf("(M)igration rate can be asymmetric\n");
      printf("   and Theta is the same for both populations\n");
      printf("   (Number of populations = 2)\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      switch (uppercase(input[0])) {
      case 'M':
	options->fsttype = 'M';
	fst_type('M');
	break;
      case 'T':
	options->fsttype = 'T';
	fst_type('T');
	break;
      default:
	options->fsttype = 'T';
	fst_type('T');
	break;
      }
      break;
    case '6':
      options->printfst = !options->printfst;
      break;
    case '7':
      set_menu_localities(options);
      break;
    case '8':	/* fill in the custom migration
		 * matrix */
      read_custom_menu_migration(options);
      break;
    case '9':
      options->geo = !options->geo;
      if(options->geo)
	{
	  get_filename("  What is the filename for the distance matrix file?",
		     GEOFILE, options->geofilename);
	}
      break;
#ifdef POPMODEL
    case '10':
      set_popmodel(options);
      break;
#endif
    }
  }
  while (uppercase(input[0]) != 'Y');
  myfree(outputstring);
}

void set_menu_localities(option_fmt *options)
{
  char *input;
  char *tmp;
  int i;
  input = (char *) mycalloc(LINESIZE,sizeof(char));
  tmp = (char *) mycalloc(LINESIZE,sizeof(char));
  printf("Associate sampling locations with populations\n");
  printf("---------------------------------------------\n");
  printf("Default: every sampling location is a population\n");
  printf("Combine localities this way, for example there are\n");
  printf("4 locations: 1, 2, 3, 4 can be combined into 2 populations\n");
  printf("by mapping the 4 position to 1, 1, 1, 2\n\n\n");
  printf("How many localities are in the data set? [for default simply press RETURN]\n");
  printf("> ");fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  if(input[0]!='\0')
    {
      options->newpops_numalloc = atoi(input);
      options->newpops = (long*) myrealloc(options->newpops, options->newpops_numalloc * sizeof(long));
      printf("Enter now the remappings (little checking is done with this, enter exactly the number of locations)\n");
      for(i=1;i<=options->newpops_numalloc;i++)
	{
	  switch((int) log10((double) options->newpops_numalloc))
	    {
	    case 0: 
	    case 1:
	      printf("%2i",i); break;
	    case 2:
	      printf("%3i",i); break;
	    default:  
	      printf("%4i",i); break;
	    }
	}
      printf("\n>");fflush(stdout);
      FGETS(input,LINESIZE,stdin);
      set_localities(&input,&tmp,options);
    }
  myfree(input);
  myfree(tmp);
}

///print the standard question "are the settings correct? ....." at the end of the menu
void            print_bottom_menu_part() 
{
  printf("\n\n  Are the settings correct?\n");
  printf("  (Type Y to go back to the main menu or the number for a menu to change)\n===>");
}

///
///displays and makes changes to the options that manage the analysis strategy and run condition
///
void            menuStrategy(option_fmt * options) {
  boolean         done = FALSE;
  do {
    printf("  SEARCH STRATEGY\n\n");
    if (!options->bayes_infer) {
      display_ml_mcmc(options);
      print_bottom_menu_part();
      done = menuStrategy_ml(options);
    } else {
      display_bayes_mcmc(options);
      print_bottom_menu_part();
      done = menuStrategy_bayes(options);
    }
  } while (done == FALSE);
}

//-------------------------------------------------------------------------------------
// prior menu functions

/// \brief Menu for multiplier proposal
/// Menu for multiplier proposal
void set_mult_prior(int paramgroupm, prior_fmt *prior)
{
  char input[LINESIZE];
  long count = 0;
  do {
    printf("Specify the MINIMUM, MAXIMUM and maximal MULTIPLIER\n[Current setting: {%f %f %f}\n===> ",
	   prior->min, prior->max, prior->delta);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (input[0] == '\0')
      break;
#ifdef USE_MYREAL_FLOAT
    count = sscanf(input, "%f%f%f", &prior->min, &prior->max, &prior->delta);
#else
    count = sscanf(input, "%lf%lf%lf", &prior->min, &prior->max, &prior->delta);
#endif
  } while (count != 3 && count > 0);
}

/// \brief Menu for exponential proposal
/// Menu for exponential proposal
void set_exp_prior(int paramgroupm, prior_fmt *prior)
{
  char input[LINESIZE];
  long count = 0;
  do {
    printf("Specify the MINIMUM, MEAN, MAXIMUM\n[Current setting: {%f %f %f}\n===> ",
	   prior->min, prior->mean, prior->max);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (input[0] == '\0')
      break;
#ifdef USE_MYREAL_FLOAT
    count = sscanf(input, "%f%f%f", &prior->min, &prior->mean, &prior->max);
#else
    count = sscanf(input, "%lf%lf%lf", &prior->min, &prior->mean, &prior->max);
#endif
  } while (count != 3 && count > 0);
}

/// \brief Menu for gamma proposal
/// Menu for gamma proposal, mean = alpha * beta = 1 -> 
void set_gamma_prior(int paramgroupm, prior_fmt *prior)
{
  char input[LINESIZE];
  long count = 0;
  do {
    printf("Specify the MINIMUM, ALPHA, MAXIMUM\n[Current setting: {%f %f %f}\n===> ",
	   prior->min, prior->mean, prior->max);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if(input[0]=='\0')
      return;
#ifdef USE_MYREAL_FLOAT
    count = sscanf(input, "%f%f%f", &prior->min, &prior->mean, &prior->max);
#else
    count = sscanf(input, "%lf%lf%lf", &prior->min, &prior->mean, &prior->max);
#endif
  } while (count != 3 && count > 0);
}

/// \brief Menu for exponential with window proposal
/// Menu for exponential with window proposal
void set_wexp_prior(int paramgroupm, prior_fmt *prior)
{
  char input[LINESIZE];
  long count = 0;
  do {
    printf("Specify the MINIMUM, MEAN, MAXIMUM and DELTA\nUse for DELTA about 1/10 of the MIN-MAX range\n[Current setting: {%f %f %f %f}\n===> ",
	   prior->min, prior->mean, prior->max, prior->delta);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if(input[0]=='\0')
      return;
#ifdef USE_MYREAL_FLOAT
    count = sscanf(input, "%f%f%f%f", &prior->min, &prior->mean, &prior->max, &prior->delta);
#else
    count = sscanf(input, "%lf%lf%lf%lf", &prior->min, &prior->mean, &prior->max, &prior->delta);
#endif
  } while (count != 4 && count > 0);
}

/// \brief Menu for multiplier proposal
/// Menu for multiplier proposal
void set_uni_prior(int paramgroupm, prior_fmt *prior)
{
  char input[LINESIZE];
  long count = 0;
  do {
    printf("Specify the MINIMUM, MAXIMUM, DELTA\nUse for DELTA about 1/10 of the MIN-MAX range\n[Current setting: {%f %f %f}\n===> ",
	   prior->min, prior->max, prior->delta);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if(input[0]=='\0')
      return;
#ifdef USE_MYREAL_FLOAT
    count = sscanf(input, "%f%f%f", &prior->min, &prior->max, &prior->delta);
#else
    count = sscanf(input, "%lf%lf%lf", &prior->min, &prior->max, &prior->delta);
#endif
  } while (count != 3 && count > 0);
  prior->mean = (prior->max + prior->min) / 2.;
  //prior->delta = (prior->max - prior->min) / 10.;
}

/// \brief return "set" "-" depending on prior set
/// returns a short string that shows whether this prior distribution is set
char * is_prior(int priortype, int priorset)
{
  return priortype == priorset ? "Set" : "-";
}

/// \brief returns priortype sting
/// returns a string that shows what prior distribution is set
char * is_proposaltype(boolean proposalset)
{
  switch(proposalset)
    {
    case TRUE: return  "Slice sampling" ; 
    default: return "Metropolis sampling";
    }
}


/// \brief Menu to set all prior distributions
/// menu to set proposal distribution, returns TRUE when standard course of action
/// FALSE when a problem occured.
boolean set_proposal_menu(int paramgroup, option_fmt *options)
{
  char            input[LINESIZE];
  boolean proposal;

  switch(paramgroup)
    {
    case THETAPRIOR:
      printf("Proposal setting for all population sizes [Theta]:\n");
      proposal = options->slice_sampling[THETAPRIOR];
      break;
    case MIGPRIOR:
      printf("Proposal setting for all migration rates [%s]:\n",options->usem ? "M" : "Theta*M");
      proposal = options->slice_sampling[MIGPRIOR];
      break;
    case RATEPRIOR:
      printf("Proposal setting for mutation rate modifier [r]:\n");
      proposal = options->slice_sampling[RATEPRIOR];
      break;
    default:
      return FALSE;
    }
    input[0] = '\0';
    printf("Which proposal distribution? [current: %11s] (choices: Slice or Metropolis)\n===> ", 
    	   is_proposaltype(proposal));
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    switch (uppercase(input[0])) 
      {
      case 'M':
	options->slice_sampling[paramgroup] = FALSE;
	break;
      case 'S':
      default:
	options->slice_sampling[paramgroup] = TRUE;
	break;
      }
    input[0]='\0';
    return TRUE;
  }


/// \brief returns priortype sting
/// returns a string that shows what prior distribution is set
char * is_priortype(int priorset)
{
  switch(priorset)
    {
      //    case SLICE: return  "Uniform/Slice" ; 
    case MULTPRIOR: return  "Multiplier" ; 
    case EXPPRIOR: return "Exponential" ; 
    case WEXPPRIOR: return "Exp window";
    case GAMMAPRIOR: return "Gamma";
    case UNIFORMPRIOR: return "Uniform";
    default: return "-";
    }
}


/// \brief Menu to set all prior distributions
/// menu to set all prior distributions, returns TRUE when standard course of action
/// FALSE when a problem occured.
boolean set_prior_menu(int paramgroup, option_fmt *options)
{
  int tmp;
  char            input[LINESIZE];
  prior_fmt *prior;

  switch(paramgroup)
    {
    case THETAPRIOR:
      printf("Prior setting for all population sizes [Theta]:\n");
      prior = options->bayespriortheta;
      break;
    case MIGPRIOR:
      printf("Prior setting for all migration rates [%s]:\n",options->usem ? "M" : "Theta*M");
      prior = options->bayespriorm;
      break;
    case RATEPRIOR:
      printf("Prior setting for mutation rate modifier [r]:\n");
      prior = options->bayespriorrate;
      break;
    default:
      return FALSE;
    }
  do{
    input[0] = '\0';
    //    printf("\nSet the prior distribution and its parameters\n\n");
    //printf("  %i   Set Slice sampler?                             %11s\n", SLICE, 
    //	   is_prior(SLICE,options->bayesprior[paramgroup]));
    printf("  %i   Set Multiplication prior distribution?         %11s\n", MULTPRIOR, 
	   is_prior(MULTPRIOR,options->bayesprior[paramgroup]));
    printf("  %i   Set Exponential prior distribution?            %11s\n", EXPPRIOR, 
	   is_prior(EXPPRIOR,options->bayesprior[paramgroup]));
    printf("  %i   Set Exponential prior with window distribution?%11s\n", WEXPPRIOR,
	   is_prior(WEXPPRIOR,options->bayesprior[paramgroup]));
    //    printf("  %i   Set Gamma prior distribution?                  %11s\n", GAMMAPRIOR, 
    //	   is_prior(GAMMAPRIOR,options->bayesprior[paramgroup]));
    printf("  %i   Set Uniform prior distribution?                %11s\n\n", UNIFORMPRIOR,
	   is_prior(UNIFORMPRIOR,options->bayesprior[paramgroup]));
    printf("  (Type Y to go back to the main menu or the letter for an entry to change)\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if(uppercase(input[0]) == 'Y')
      return TRUE;
    input[1] = '\0';
    tmp = atoi(input);
    switch (tmp) 
      {
      case UNIFORMPRIOR:
	set_uni_prior(paramgroup, prior);
	options->bayesprior[paramgroup] = UNIFORMPRIOR;	
	break;
      case EXPPRIOR:
	set_exp_prior(paramgroup, prior);
	options->bayesprior[paramgroup] = EXPPRIOR;	
	break;
      case WEXPPRIOR:
	set_wexp_prior(paramgroup, prior);
	options->bayesprior[paramgroup] = WEXPPRIOR;	
	break;
      case MULTPRIOR:
	set_mult_prior(paramgroup, prior);
	options->bayesprior[paramgroup] = MULTPRIOR;	
	break;
      case GAMMAPRIOR:
	set_gamma_prior(paramgroup, prior);
	options->bayesprior[paramgroup] = GAMMAPRIOR;	
	break;
      default:
	set_uni_prior(paramgroup, prior);
	options->bayesprior[paramgroup] = UNIFORMPRIOR;	
	break;
      }
  }    while (uppercase(input[0]) != 'Y');
  input[0]='\0';
  return TRUE;
}


///
/// sets proposal distribution
void
menuProposal(option_fmt * options) {
  char            input[LINESIZE];
  do 
    {
      input[0] = '\0';
      printf ("  Proposal distribution setting [You still need to set the PRIOR distribution!]:\n");
      printf("  1   Set proposal distribution for Theta?          %11s\n", 
	     is_proposaltype(options->slice_sampling[THETAPRIOR]));
      printf("  2   Set proposal distribution for migration?      %11s\n",
	     is_proposaltype(options->slice_sampling[MIGPRIOR]));
      if(options->bayesmurates)
	{      
	  printf("  3   Set mutation rate modifier prior distribution?%11s\n",
	     is_proposaltype(options->slice_sampling[RATEPRIOR]));
	}
      printf("\n  (Type Y to go back to the main menu or the letter for a menu to change)\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      switch (input[0]) 
	{
	case '1':
	  set_proposal_menu(THETAPRIOR, options);
	  break;
	case '2':
	  set_proposal_menu(MIGPRIOR, options);
	  break;
	case '3':
	  if(options->bayesmurates)
	    {      
	      set_proposal_menu(RATEPRIOR, options);
	    }
	  break;
	default:
	  break;
	}
    } while (uppercase(input[0]) != 'Y');;
}
///
/// sets prior distribution parameters
void
menuPrior(option_fmt * options) {
  char            input[LINESIZE];
  do 
    {
      input[0] = '\0';
      printf ("  Prior distribution setting:\n");
      printf("  1   Set Theta prior distribution?                 %11s\n", 
	     is_priortype(options->bayesprior[THETAPRIOR]));
      printf("  2   Set Migration prior distribution?             %11s\n",
	     is_priortype(options->bayesprior[MIGPRIOR]));
      if(options->bayesmurates)
	{      
	  printf("  3   Set mutation rate modifier prior distribution?%11s\n",
	     is_priortype(options->bayesprior[RATEPRIOR]));
	}
      printf("\n  (Type Y to go back to the main menu or the letter for a menu to change)\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      switch (input[0]) 
	{
	case '1':
	  set_prior_menu(THETAPRIOR, options);
	  break;
	case '2':
	  set_prior_menu(MIGPRIOR, options);
	  break;
	case '3':
	  if(options->bayesmurates)
	    {      
	      set_prior_menu(RATEPRIOR, options);
	    }
	  break;
	default:
	  break;
	}
    } while (uppercase(input[0]) != 'Y');;
}

///
///reads bayes menu choices and sets options accordingly
boolean menuStrategy_bayes(option_fmt * options) {
  char            input[LINESIZE];
  char           *priorstring;
  long            count = 0;
  char           *extension;
  //read user input
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  //if user types yes or YES or YES exit and tell menuStrategy that we are done
    if (strchr("Yy", input[0]))
      return TRUE;

  priorstring = (char *) mycalloc(LINESIZE, sizeof(char));
  //make changes to options
  switch (atoi(input)) {
  case BAYESSTRATEGY:
    //strategy change

    if (strchr(input,'0'))
      {
	options->bayes_infer = !options->bayes_infer;
	if(options->bayes_infer)
	  options->lchains=1;
      }
    input[0] = '\0';
    break;
  case BAYESOUT:
    printf(" Save posterior distribution (frequency histogram values)? [YES]\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (uppercase(input[0]) != 'N') 
      {
	options->has_bayesfile = TRUE;
	get_filename("  What is the filename for the posterior distribution (frequency histogram)?",
		     BAYESFILE, options->bayesfilename);
      }
    else
      {
	options->has_bayesfile = FALSE;
      }
    break;
  case BAYESMDIMOUT:
    printf(" Save posterior distribution (all raw parameter values)? [NO]\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (uppercase(input[0]) == 'Y') 
      {
	options->has_bayesmdimfile = TRUE;
	get_filename("  What is the filename for the complete posterior distribution (raw parameter values)?\n[if the extension of the file is .gz then the will be compressed]",
		     BAYESMDIMFILE, options->bayesmdimfilename);
	unpad(options->bayesmdimfilename, " ");
	extension = strrchr(options->bayesmdimfilename,'.');
	if(extension!=NULL && !strncmp(extension,".gz",3))
	  {
	    options->use_compressed = 1;
	  }
	else
	  {
	    options->use_compressed = 0;
	  }  	
	do
	  {
	    printf("Sampling interval for the raw parameter values\n");
	    printf("[Small values can result in a HUGE file!!]\n");
	    printf("[Default is the same as the sampling increment]\n");
	    printf("[Examples: 1 --> saving all parameters every sampling increment]\n");
	    printf("[          2 --> saving every second sampling increment]\n");
	    printf("[        100 --> saving only very hundredth sample]\n");
	    printf("[                if long-inc=50 and here we have here 100 then only]\n");
	    printf("[                every 5,000th step is saved to file]\n");
	    printf("[                with many loci and large run this is still a large file]\n");
	    printf("Current setting: %li\n===> ",options->bayesmdiminterval);
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    if(input[0] == '\0' && options->bayesmdiminterval > 0)
	      break;
	    else
	      options->bayesmdiminterval = atol(input);
	  }
	while(options->bayesmdiminterval<1);
      }
    else
      {
	options->has_bayesmdimfile = FALSE;
      }
    break;
  case BAYESBINNING:
    do {
      printf("Specify the number of bins for THETA and %s!\n", options->usem ? "M" : "Theta*M");
      printf("Current setting:  Bins      Interval width\n"); 
      printf("  Theta           %li        %f\n", options->bayespriortheta->bins, 
	     (options->bayespriortheta->max - options->bayespriortheta->min) / (float) options->bayespriortheta->bins);
      printf("%7.7s            %li        %f\n", options->usem ? "M" : "Theta*M", options->bayespriorm->bins,
	     (options->bayespriorm->max - options->bayespriorm->min) / (float) options->bayespriorm->bins);
      if(options->bayesmurates)
	printf("  Rate            %li        %f\n", options->bayespriorrate->bins, 
	       (options->bayespriorrate->max - options->bayespriorrate->min) / (float) options->bayespriorrate->bins);
      printf("\n==> ");	     
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (input[0] == '\0')
	break;
      if(options->bayesmurates)
	count = sscanf(input, "%li%li%li", &options->bayespriortheta->bins, &options->bayespriorm->bins, &options->bayespriorrate->bins);
      else
	count = sscanf(input, "%li%li", &options->bayespriortheta->bins, &options->bayespriorm->bins);
    } while (count != (options->bayesmurates ? 3 : 2) && count > 0);
    break;
  case BAYESPRETTY:
    printf("Specify the method for plotting the posterior distribution!\n");
    printf("Valid options are:\n");
    printf("ALL     use the range of the prior distribution\n");
    printf("P99     do not plot the values higher than the 99%% percentile for each parameter\n");
    printf("MAXP99  plot up to the maximum 99%% percentile of all parameters.\n");
    printf("TOTAL   plot up to the 100%% percentile of all parameters.\n");
    printf("Current setting: %s\n",
	   (options->bayespretty == PRETTY_P99 ? "P99" :
	   (options->bayespretty == PRETTY_MAX ? "ALL" : 
	   (options->bayespretty == PRETTY_P100 ? "TOTAL" : "MAXP99"))));
    printf("\n\n==> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (input[0] == '\0')
      break;
    switch(uppercase(input[0]))
      {
      case 'A':
	options->bayespretty = PRETTY_MAX;
	break;
      case 'P':
	options->bayespretty = PRETTY_P99;
	break;
      case 'T':
	options->bayespretty = PRETTY_P100;
	break;
      case 'M':
      default:
	options->bayespretty = PRETTY_P99MAX;
	break;
      }
    break;
  case BAYESFREQ:
    do {
      printf("What is the frequency of the tree updates vs the parameter updates?\n");
      printf("[A frequncy of 1.0 updates only trees, frequency of 0.0 updates only parameters]\n");
      printf("[Suggestion for small problems: 0.5, actual value: %f]\n===> ", options->updateratio);
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (input[0] != '\0')  
	  options->updateratio = atof(input);
    } while (options->updateratio > 1);
    break;
  case BAYESPROPOSAL:
    menuProposal(options);
    break;
  case BAYESPRIOR:
    menuPrior(options);
    break;
  case BAYESLCHAINS:
    do {
      printf("  How many Long Chains?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->lchains = atoi(input);
      if (options->lchains < 0)
	printf("  Must be non-negative\n");
    }
    while (options->lchains < 0);
    break;
	   case BAYESSKIP:
	   do {
	     printf("  How many steps (tree changes, parameter changes) to skip?\n===> ");
	     fflush(stdout); FGETS(input, LINESIZE, stdin);
	     options->lincrement = atoi(input);
	     if (options->lincrement <= 0)
	       printf("  Must be positive\n");
	   }
	   while (options->lincrement <= 0);
	   break;
	   case BAYESSAMPLES:
	   do {
	     printf("  How many trees and parameter steps to record?\n===> ");
	     fflush(stdout); FGETS(input, LINESIZE, stdin);
	     options->lsteps = atoi(input);
	     if (options->lsteps <= 0)
	       printf("  Must be a positive integer\n");
	   }
	   while (options->lsteps <= 0);
	   break;
	   case BAYESBURNIN:
	   do {
	     printf("  How many steps to discard? (burn-in)\n===> ");
	     fflush(stdout); FGETS(input, LINESIZE, stdin);
	     if(strchr(input,'a'))
	       {
		 options->burnin_autostop = 'a';
	       }
	     else
	       {
		 options->burnin_autostop = ' ';
	       }
	     options->burn_in = atoi(input);
	     if (options->burn_in <= 0)
	       printf("  Must be a positive integer or zero (0)\n");
	   }
	   while (options->burn_in < 0);
	   break;
	   case BAYESREPLICATE:
	   do {
	     printf("  Summarize over multiple chains? [YES | NO] \n===> ");
	     fflush(stdout); FGETS(input, LINESIZE, stdin);
	     if (uppercase(input[0]) == 'Y') 
	       {
		 options->replicate = TRUE;
		 printf("  How many independent runs\n===> ");
		 fflush(stdout); FGETS(input, LINESIZE, stdin);
		 options->replicatenum = ATOL(input);
		 if (options->replicatenum < 1)
		   printf("  Enter a number >= 1\n");
	       } 
	     else 
	       {
		 options->replicate = FALSE;
		 options->replicatenum = 0;
	       }
	   }
	   while (options->replicatenum < 0);
	   break;
	   case BAYESHEAT:
	   printf
	   ("  Heating scheme? < NO | YES | STATIC | ADAPTIVE | BOUNDED_ADAPTIVE >\n===> ");
	   fflush(stdout); FGETS(input, LINESIZE, stdin);
	   switch (tolower(input[0])) {
	   case 'a':
	     options->heating = 1;
	     options->adaptiveheat = STANDARD;
	     options->heating_interval = 1;
	     menuHeat(options, input);
	     break;
	   case 'b':
	     options->heating = 1;
	     options->adaptiveheat = BOUNDED;
	     options->heating_interval = 1;
	     menuHeat(options, input);
	     break;
	   case 'y':
	   case 's':
	     options->heating = 1;
	     options->adaptiveheat = NOTADAPTIVE;
	     options->heating_interval = 1;
	     menuHeat(options, input);
	     break;
	   case '\0':
	   case 'n':
	   default:
	     options->heating = 0;
	     options->adaptiveheat = NOTADAPTIVE;
	     options->heated_chains = 1;
	     break;
	   }
	   break;
  case BAYESMOVINGSTEPS:
    options->movingsteps = !options->movingsteps;
    if (options->movingsteps) {
      do {
	printf
	  ("  How big should the fraction of new genealogies\n");
	printf
	  ("  of the originally proposed number of samples be?\n[Use this option only when acceptance ratio is small (<0.01)\nRuntime can increase tremendeously]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	options->acceptfreq = atof(input);
	if (options->acceptfreq < 0 || options->acceptfreq > 1)
	  printf("  Range should be between 0 - 1, and not %f\n",
		 options->acceptfreq);
      }
      while (options->acceptfreq < 0 || options->acceptfreq > 1);
    }
    break;
  case BAYESGELMAN:
    printf("  Convergence statistic options: Pairs , Sum, NO\n");
    printf("  [Currently set to %15s]\n===> ",
	   options->gelman ? (options->gelmanpairs ? "YES:Pairs" : "YES:Summary" ) : " NO");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (input[0] == '\0')
      break;
    switch(uppercase(input[0]))
      {
      case 'N':
	options->gelman = FALSE;
	options->gelmanpairs = FALSE;
	break;
      case 'P':
	options->gelman = TRUE;
	options->gelmanpairs = TRUE;
	break;
      default:
	options->gelman = TRUE;
	options->gelmanpairs = FALSE;
	break;
      }
    break;
  case BAYESPRIORALONE:
    printf
      ("  Show only the prior distributions? < NO | YES >\nwith YES no data is analyzed!!!\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    switch (tolower(input[0])) {
    case 'y':
      options->prioralone = TRUE;
      break;
    case 'n':
    default:
      options->prioralone = FALSE;
    }
    break;
  default:
    break;
  }
  myfree(priorstring);
  return FALSE;
}


boolean         menuStrategy_ml(option_fmt * options) 
{
  char            input[LINESIZE];
  //read user input
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  //if user types yes or Yes or YES exit and tell menuStrategy that we are done
  if              (strchr("Yy", input[0]))
    return TRUE;
  
  //make changes to options
  switch          (atoi(input)) {
  case MLSTRATEGY://strategy change
    if (strchr(input,'0'))
      options->bayes_infer = !options->bayes_infer;
    input[0] = '\0';
    break;
  case MLSHORTCHAINS:
    do {
      printf("  How many Short Chains?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->schains = atoi(input);
      if (options->schains < 0)
	printf("  Must be non-negative\n");
    }
    while           (options->schains < 0);
    break;
  case MLSHORTSKIP:
    do {
      printf("  How many trees to skip?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->sincrement = atoi(input);
      if (options->sincrement <= 0)
	printf("  Must be positive\n");
    }
    while (options->sincrement <= 0);
    break;
  case MLSHORTSAMPLES:
    do {
      printf("  How many trees to record?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->ssteps = atoi(input);
      if (options->ssteps <= 0)
	printf("  Must be a positive integer\n");
    }
    while (options->ssteps <= 0);
    break;
  case MLLONGCHAINS:
    do {
      printf("  How many Long Chains?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->lchains = atoi(input);
      if (options->lchains < 0)
	printf("  Must be non-negative\n");
    }
    while (options->lchains < 0);
    break;
  case MLLONGSKIP:
    do {
      printf("  How many trees to skip?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->lincrement = atoi(input);
      if (options->lincrement <= 0)
	printf("  Must be positive\n");
    }
    while (options->lincrement <= 0);
    break;
  case MLLONGSAMPLES:
    do {
      printf("  How many trees to record?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->lsteps = atoi(input);
      if (options->lsteps <= 0)
	printf("  Must be a positive integer\n");
    }
    while (options->lsteps <= 0);
    break;
  case MLBURNIN:
    do {
      printf("  How many genealogies to discard?\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      options->burn_in = atoi(input);
      if (options->burn_in <= 0)
	printf("  Must be a positive integer or zero (0)\n");
    }
    while (options->burn_in < 0);
    break;
  case MLREPLICATE:
    options->replicate = !options->replicate;
    if (options->replicate) {
      do {
	printf("  Summarize over (L)ong chains?\n");
	printf("  or (M)ultiple runs ?\n");
	printf("  [Default is (L)]\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	if (uppercase(input[0]) == 'M') {
	  printf("  How many independent runs\n===> ");
	  fflush(stdout); FGETS(input, LINESIZE, stdin);
	  options->replicatenum = ATOL(input);
	  if (options->lcepsilon < 0)
	    printf("  Enter a number >= 1\n");
	} else
	  options->replicatenum = 0;
      }
      while (options->replicatenum < 0);
    } else {
      options->replicatenum = 0;
    }
    break;
  case MLHEAT:
    printf
      ("  Heating scheme? < NO | YES | ADAPTIVE | BOUNDED_ADAPTIVE>\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    switch (tolower(input[0])) {
    case 'a':
      options->heating = 1;
      options->adaptiveheat = STANDARD;
      options->heating_interval = 1;
      menuHeat(options, input);
      break;
    case 'b':
      options->heating = 1;
      options->adaptiveheat = BOUNDED;
      options->heating_interval = 1;
      menuHeat(options, input);
      break;
    case 'y':
      options->heating = 1;
      options->adaptiveheat = NOTADAPTIVE;
      options->heating_interval = 1;
      menuHeat(options, input);
      break;
    case '\0':
    case 'n':
    default:
      options->heating = 0;
      options->adaptiveheat = NOTADAPTIVE;
      break;
    }
    break;
  case MLMOVINGSTEPS:
    options->movingsteps = !options->movingsteps;
    if (options->movingsteps) {
      do {
	printf
	  ("  How big should the fraction of new genealogies\n");
	printf
	  ("  of the originally proposed number of samples be?\n===> ");
	fflush(stdout); FGETS(input, LINESIZE, stdin);
	options->acceptfreq = atof(input);
	if (options->acceptfreq < 0 || options->acceptfreq > 1)
	  printf("  Range should be between 0 - 1, and not %f\n",
		 options->acceptfreq);
      }
      while (options->acceptfreq < 0 || options->acceptfreq > 1);
    }
    break;
  case MLEPSILON:
    do {
      printf("  Parameter likelihood epsilon?\n[INF is Default]\n===> ");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (uppercase(input[0]) == 'I')
	options->lcepsilon = LONGCHAINEPSILON;
      else
	options->lcepsilon = atof(input);
      if (options->lcepsilon <= 0)
	printf
	  ("  Must be a positive value, be warned: too small values will run the program forever\n");
    }
    while (options->lcepsilon <= 0);
    break;
  case MLGELMAN:
      printf("  Convergence statistic options: Pairs , Sum, NO\n");
      printf("  [Currently set to %10s]\n===> ",
	 options->gelman ? (options->gelmanpairs ? "YES:Pairs" : "YES:Summary" ) : " NO");
      fflush(stdout); FGETS(input, LINESIZE, stdin);
      if (input[0] == '\0')
	break;
      switch(uppercase(input[0]))
	{
	case 'N':
	  options->gelman = FALSE;
	  options->gelmanpairs = FALSE;
	  break;
	case 'P':
	  options->gelman = TRUE;
	  options->gelmanpairs = TRUE;
	  break;
	default:
	  options->gelman = TRUE;
	  options->gelmanpairs = FALSE;
	  break;
	}
      break;
  }
  return FALSE;
}

void            menuHeat(option_fmt * options, char *input) 
{
  if (options->heating > 0) {
    printf("Enter the number of different \"heated\" chains.\nMinimum is 4\n===> ");
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    options->heated_chains = atol(input);
    if (options->heated_chains < 4)
      options->heated_chains = HEATED_CHAIN_NUM;
    printf("Enter the interval between swapping trees\nEnter 0 (zero) for NO swapping\n[Current interval is %li]\n===> ",
	    options->heating_interval);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (input[0] == '\0')
      return;
    if (atol(input) > 0) {
      options->heating_interval = atol(input);
      options->heating = 1;
    } else {
      options->heating = 1;
      options->heating_interval = 0;
      options->heatedswap_off=TRUE;
    }
    read_heatvalues(options);
  }
}

void
display_ml_mcmc(option_fmt * options) {
  char            temp[LINESIZE];

  printf("  0   Strategy: %42.42s\n", options->bayes_infer ?
	 "Bayesian Inference" :
	 "Maximum Likelihood");
  printf("  1   Number of short chains to run?                %6ld\n",
	 options->schains);
  if              (options->schains > 0) {
    printf
      ("  2   Short sampling increment?                     %6ld\n",
       options->sincrement);
    printf
      ("  3   Number of recorded genealogies in short chain?%6ld\n",
       options->ssteps);
  }
  printf("  4   Number of long chains to run?                 %6ld\n",
	 options->lchains);
  if (options->lchains > 0) {
    printf
      ("  5   Long sampling increment?                      %6ld\n",
       options->lincrement);
    printf
      ("  6   Number of recorded genealogies in long chain? %6ld\n",
       options->lsteps);
  }
  printf("  7   Number of genealogies to discard at \n");
  printf("      the beginning of each chain? [Burn-in]       %c%6ld\n",
	 options->burnin_autostop,options->burn_in);
  if (!options->replicate)
    printf
      ("  8   Combine chains or runs for estimates?             NO\n");
  else {
    if (options->replicatenum == 0)
      printf
	("  8   Combine chains for estimates?       YES, long chains\n");
    else
      printf
	("  8   Combine chains for estimates?     YES, over %3li runs\n",
	 options->replicatenum);
  }
  if (options->heating == 0)
    printf
      ("  9   Heating:                                          NO\n");
  else {
    if (options->adaptiveheat!=NOTADAPTIVE)
      sprintf(temp, "Adaptive (%s %3li chains, swap interval is %3li)\n",
	      options->adaptiveheat == STANDARD ? "Standard" : "Bounded", 
	      options->heated_chains, options->heating_interval);
    else
      sprintf(temp, "     YES (%3li chains, swap interval is %3li)\n",
	      options->heated_chains, options->heating_interval);
    printf("  9   Heating: %s", temp);
  }
  printf
    ("\n -------------------------------------------------------------\n");
  printf(" Obscure options (consult the documentation on these)\n\n");
  if (options->movingsteps)
    printf
      (" 10   Sample a fraction of %2.2f new genealogies?       YES\n",
       options->acceptfreq);
  else
    printf
      (" 10   Sample at least a fraction of new genealogies?    NO\n");
  if (options->lcepsilon < LONGCHAINEPSILON)
    sprintf(temp, "%17.2f", options->lcepsilon);
  else
    sprintf(temp, "%17.17s", "infinity");
  printf(" 11   Epsilon of parameter likelihood    %s\n", temp);
  printf(" 12   Use Gelman's convergence criterium?    %13s\n",
	 options->gelman ? (options->gelmanpairs ? "YES:Pairs" : "YES:Summary" ) : " NO");
}

void
display_bayes_mcmc(option_fmt * options) {
  char           *outputstring;
  outputstring = (char *) mycalloc(LINESIZE, sizeof(char));

  printf("  %2i   Strategy:%54.54s\n", BAYESSTRATEGY, options->bayes_infer ?
	 "Bayesian Inference" :
	 "Maximum Likelihood");
  if(options->has_bayesfile)
    sprintf(outputstring,"%s%s",  "YES:", options->bayesfilename);
  else
    sprintf(outputstring,"NO");
  printf(  "  %2i   File for recording posterior distribution?%21s\n", BAYESOUT, outputstring);
  if(options->has_bayesmdimfile)
    sprintf(outputstring,"%s%s",  "YES:", options->bayesmdimfilename);
  else
    sprintf(outputstring,"NO");
  printf(  "  %2i   File for recording all parameter values?  %21s\n", BAYESMDIMOUT, outputstring);
  if(options->has_bayesmdimfile)
    {
      sprintf(outputstring,"[save every %li x sample-increment = %li]", options->bayesmdiminterval, 
	      options->bayesmdiminterval*options->lincrement);
      printf(  "            %58.58s\n", outputstring);
    }
  sprintf(outputstring, "%li, %li", options->bayespriortheta->bins, options->bayespriorm->bins);
  printf("  %2i   Number of bins of posterior [Theta,M]?%25s\n", BAYESBINNING, outputstring);
  printf("  %2i   Plotting type of posterior distribution?%23s\n", BAYESPRETTY, 
	 (options->bayespretty == PRETTY_P99 ? "up to 99% percentile" :
	 (options->bayespretty == PRETTY_MAX ? "prior distr. range" :
	 (options->bayespretty == PRETTY_P100 ? "up to ~100% percentile" : "up to maximal 99%"))));
  printf("  %2i   Frequency of tree updates vs. parameter updates?         %6.2f\n",
	 BAYESFREQ, options->updateratio);
  set_proposal(outputstring, options->slice_sampling, !options->bayesmurates);
  printf("  %2i   Proposal distribution?%41.41s\n",
	 BAYESPROPOSAL, outputstring);
  set_prior(outputstring, options->bayesprior, !options->bayesmurates);
  printf("  %2i   Prior distribution?%44.44s\n",
	 BAYESPRIOR, outputstring);
  printf("  %2i   Number of long chains to run?                            %6ld\n",
	 BAYESLCHAINS, options->lchains);
  printf("  %2i   Sampling increment?                                      %6ld\n",
	 BAYESSKIP, options->lincrement);
  printf("  %2i   Number of recorded steps in chain                  %12ld\n",
	 BAYESSAMPLES, options->lsteps);
  printf("  %2i   Number of steps to discard at \n",
	 BAYESBURNIN);
  printf("       the beginning of chain? [Burn-in]                 %c%12ld\n",
	 options->burnin_autostop, options->burn_in);

  if(!options->replicate)
    sprintf(outputstring,"NO");
  else
    sprintf(outputstring,"YES (%li independent chains)",options->replicatenum);
  printf("  %2i   Running multiple replicates:  %33.33s\n", BAYESREPLICATE, outputstring);

  if(options->heating == 0)
    printf("  %2i   Heating:            %43.43s\n", BAYESHEAT, "NO");
  else
    printf("  %2i   Heating:                       %10s (%3li parallel chains)\n",
	   BAYESHEAT, (options->adaptiveheat==STANDARD ? "ADAPTIVE" : (options->adaptiveheat==BOUNDED ? "BOUNDED" : "STATIC")) , options->heated_chains);
  sprintf(outputstring,"%f",options->acceptfreq);
  printf("  %2i   Sampling at least fraction of new genealogies:       %10.10s\n", BAYESMOVINGSTEPS, outputstring);
  sprintf(outputstring,"%s",options->gelman ? (options->gelmanpairs ? "YES:Pairs" : "YES:Summary") : "NO");
  printf("  %2i   Convergence diagnostic for replicates:            %13.13s\n", BAYESGELMAN, outputstring);

  printf("  %2i   Run analysis without data:                        %13.13s\n", BAYESPRIORALONE, options->prioralone ? "YES" : "NO");

  myfree(outputstring);
}

void
display_lratio_menutext(void) {
  printf("  H0: the ML-estimates and your values are the same\n");
  printf("  H1: the ML-estimates and your values are different\n");
  printf("  You need to specify values for ALL parameters\n");
  printf("  Beware! If you restrict the migration model then some test\n");
  printf("  will not work. This test assumes that your hypothesis\n");
  printf("  is nested in the hypothesis you test against.\n");
  printf("  Currently you cannot specify \"<\" or \">\" comparisons AND\n");
  printf("  averages of migration rates are calculated and not\n");
  printf("  averages  for number of migrants\n\n");
}

void
display_lrt_syntax(void) {
  printf("  Specify values  as an {n x n} matrix\n");
  printf("  Theta values are on the diagonal, migration rates are\n");
  printf("  off-diagonal, newlines are allowed, and\n");
  printf("  values must be separated by space or commas or tabs.\n");
  printf("\n  Syntax:\n");
  printf("      * = sames as the MLE, and not estimated again\n");
  printf("      s = averages of MIGRATION RATES from i->j and j->i\n");
  printf("      m = averages of THETAS or MIGRATION RATES [!!!]\n");
  printf("      value = specify your own value\n\n");
  printf("  [If you want to leave this editor, type \"quit\" or \"end\" and <return>\n");
  //To review the Syntax type \ "syntax\"\n");
}

char
ask_lratio_type(lr_data_fmt * data, long counter) {
  char            answer = '\0';
  char            input[LINESIZE];
  printf("  Is the test against the Maximum likelihood estimates or\n");
  printf("   or against arbitrary values? [(M)LE, (A)rbitrary]\n===> ");
  fflush(stdout); FGETS(input, 1024, stdin);
  switch          (uppercase(input[0])) {
  case 'M':
    data[counter].type = MLE;
    answer = 'M';
    break;
  case 'A':
    answer = 'A';
    data[counter].type = ARBITRARY;
    break;
  case 'Q':
  case 'E':
    answer = 'Q';
    break;
  case 'S':
    answer = 'S';
    break;
  default:
    data[counter].type = MLE;
    break;
  }
  return answer;
}

void
read_custom_menu_lratio(option_fmt * options) {
  //static long   numpop = 0, numpop2 = 0;
  //char input[LINESIZE];
  //long i, z = 0, pointer;
  //char *temp;
  //char *valpointer;
  char            answer = '\0';
  //boolean done = FALSE;
  lratio_fmt     *lratio = options->lratio;
  display_lratio_menutext();
  display_lrt_syntax();
  while           (answer != 'Q') {
    //answer = ask_lratio_type(lratio->data, lratio->counter);
    answer = 'M';
    lratio->data[lratio->counter].type = MLE;
    switch (answer) {
    case 'S':
      display_lrt_syntax();
      break;
    case 'Q':
      break;
    case 'A':
    case 'M':
      how_many_pop(&options->numthetag);
      answer = dialog_lrt(options->numthetag, lratio);
      break;
    }
  }
}

void
how_many_pop(long *numpop) {
  char            input[LINESIZE];
  if              (*numpop == 0) {
    do {
      printf("  How many populations?\n===> ");
      fflush(stdout); FGETS(input, 1024, stdin);
      *numpop = atoi(input);
    }
    while           (*numpop <= 0 && *numpop < 100);
  } else
    printf("  The data set contains %li populations\n", *numpop);
}

void
check_lrt_allocation(lratio_fmt * lratio) {
  long            i;
  if              (lratio->counter + 1 == lratio->alloccounter) {
    lratio->alloccounter += 2;
    lratio->data =
      (lr_data_fmt *) myrealloc(lratio->data, sizeof(lr_data_fmt) *
				lratio->alloccounter);
    for (i = lratio->counter + 1; i < lratio->alloccounter; i++) {
      lratio->data[i].elem = 0;
      lratio->data[i].value1 =
	(char *) mycalloc(1, sizeof(char) * LONGLINESIZE);
      lratio->data[i].value2 =
	(char *) mycalloc(1, sizeof(char) * LONGLINESIZE);
      lratio->data[i].connect =
	(char *) mycalloc(1, sizeof(char) * LONGLINESIZE);

    }
  }
}

char
menuread_lrt_paramvalue(char *value, long *counter, long numpop2) {
  char           *valptr = value;
  long            pointer = 0;
  long            z = 0;
  char           *temp;
  char            input[LONGLINESIZE];
  while           (z < numpop2) {
    fflush(stdout); FGETS(input, LONGLINESIZE, stdin);
    temp = strtok(input, ", \n\t");
    if (temp != NULL) {
      if (strchr("EQ", uppercase(temp[0]))) {
	//(*counter)--;
	return 'Q';
      }
    }
    while           (temp != NULL) {
      sprintf(valptr + pointer, " %s,", temp);
      z++;
      pointer += (long) strlen(temp) + 2;
      temp = strtok(NULL, ", \n\t");
    }
  }
  return ' ';
}

char
dialog_lrt(long numpop, lratio_fmt * lratio) {
  //char          input[LINESIZE];
  //long z = 0;
  char            answer = ' ';
  //char *valptr1, *valptr2;
  long            numpop2 = numpop * numpop;
  check_lrt_allocation(lratio);
  printf("  %li. Likelihood ratio test\n", lratio->counter);
  printf("       Enter now the %li values for the FIRST parameter set:\n",
	 numpop2);
  lratio->data[lratio->counter].type = MLE;
  answer = menuread_lrt_paramvalue(lratio->data[lratio->counter].value1,
				   &lratio->counter, numpop2);
  if              (answer == 'Q')
    return answer;
  if              (lratio->data[lratio->counter].type == MLE) {
    lratio->counter++;
    return answer;
  }
  printf("       Enter now the values for the SECOND parameter set:\n");
  answer = menuread_lrt_paramvalue(lratio->data[lratio->counter].value1,
				   &lratio->counter, numpop2);
  if (answer == 'Q')
    return answer;
  lratio->counter++;
  return answer;
}

void
read_custom_menu_migration(option_fmt * options) {
  char            input[LINESIZE];
  long            z = 0, numpop = 0, numpop2;
  printf("  Specify the migration model as an {n x n} matrix\n");
  printf("  Theta values are on the diagonal, migration rates are\n");
  printf("  off-diagonal, spaces (\" \"), \"{\", \"}\", or newlines\n");
  printf("  are allowed, but not necessary.\n");
  printf("\n  Syntax:\n");
  printf("      * = independent parameter\n");
  printf("      0 = (zero) not estimated]\n");
  printf("      c = (constant) not estimated, taken from start-parameter]\n");
  printf("      s = symmetric migration rates (M=m/mu)\n");
  printf("      S = symmetric migration rates (xNm) \n");
  printf("      m = average of each label group [not k, or s]\n"); 
  //  printf("      a,b,c,... = average of each label group [not k, or s]\n"); 
  //  printf("      1,2,3,... = if both migration rates are labeled then combine populations\n"); 
  //printf("      M = average, either ALL thetas and/or ALL migration rates\n");
 
  if (options->numthetag > 0)
    numpop = options->numthetag;
  else {
    if (options->nummg > 0)
      numpop =
	(long) ((1. + sqrt(1. + (MYREAL) options->nummg * 4.)) / 2.);
    else {
      how_many_pop(&options->numthetag);
      numpop = options->numthetag;
    }
  }
  numpop2 = numpop * numpop;
  printf("\n  You must give %li values\n===> ", numpop2);
  while (z < numpop2) {
    fflush(stdout); FGETS(input, 1024, stdin);
    read_custom_migration(stdin, options, input, numpop);
    z = (long) strlen(options->custm);
  }
}

char           *
custom_migration_type(long type) {
  switch (type) {
  case MATRIX:
    return "Full migration matrix model";
    break;
  case MATRIX_SYMMETRIC:
    return "Symmetric migration matrix model";
    break;
  case MATRIX_SAMETHETA:
    return "Full migration matrix model (same Theta)";
    break;
  case MATRIX_ARBITRARY:
    return "User specified migration matrix model";
    break;
  case ISLAND:
    return "N-population island model";
    break;
  case ISLAND_VARTHETA:
    return "N-population island model (variable Theta)";
    break;
  case STEPSTONE:
    return "Stepping stone model";
    break;
  case CONTINUUM:
    return "Continuum model";
    break;
  case NEIGHBOR:
    return "Isolation by distance model";
    break;
  default:
    return "Illegal migration model";
    break;
  }
}


void
read_heatvalues(option_fmt * options) {
  long            i, z;
  char           *tmp;
  char            input[LINESIZE];
  MYREAL          diff = 0.;
  MYREAL          x;
  FPRINTF(stdout, " ");
  printf("Enter %li \"temperatures\"\n", options->heated_chains);
  printf("[The coldest temperature, which is the first, has to be 1\n");
  printf(" For example: 1 1.5 3 6]\n");
 
  FPRINTF(stdout,
	  "OR give a range of values  [linear increase:      1 - 10]\n");
  FPRINTF(stdout,
	  "                           [exponential increase: 1 @ 10]\n");
  FPRINTF(stdout,
	  "or, most lazily, let me suggest a range [simply type a #]\n");
  printf("===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if (strstr(input, "-")) 
    {
      tmp = strstr(input, "-");
      options->heat[options->heated_chains - 1] = fabs(atof(tmp + 1));
      options->heat[0] = 1.0;
      diff = options->heat[options->heated_chains - 1] - options->heat[0];
      diff /= options->heated_chains - 1.;
      for (i = 1; i < options->heated_chains; i++)
	options->heat[i] = options->heat[i - 1] + diff;
    } 
  else 
    {
      if (strstr(input, "#")) 
	{
	  diff = 1./(options->heated_chains - 1);
	  x = 1 - diff;
	  for (i = 1; i < options->heated_chains-1; i++)
	    {
	      options->heat[i] = 1./ x;
	      x -= diff;
	    }
	  options->heat[i] = 1000000.;
	}
      else
	{
	  if (strstr(input, "@")) 
	    {
	      tmp = strstr(input, "@");
	      options->heat[options->heated_chains - 1] = fabs(atof(tmp + 1));
	      options->heat[0] = 1.0;
	      diff = 2. * (options->heat[options->heated_chains - 1] -
			   options->heat[0]);
	      diff /= pow(2.0, (MYREAL) options->heated_chains) - 2.;
	      for (i = 1; i < options->heated_chains; i++)
		options->heat[i] =
		  options->heat[i - 1] + pow(2.0, i - 1.) * diff;
	    } 
	  else 
	    {
	      z = 0;
	      tmp = strtok(input, " ,");
	      if (tmp != NULL) 
		{
		  options->heat[z++] = atof(tmp);
		}
	      for (;;) 
		{
		  tmp = strtok(NULL, " ,");
		  if (tmp != NULL) 
		    {
		      options->heat[z++] = atof(tmp);
		    } 
		  else 
		    {
		      if (z < options->heated_chains) 
			{
			  printf("%li more values are needed\n===> ",
				  options->heated_chains - z);
			  fflush(stdout); FGETS(input, LINESIZE, stdin);
			  tmp = strtok(input, " ,");
			  if (tmp != NULL) 
			    {
			      options->heat[z++] = atof(tmp);
			    }
			} 
		      else
			break;
		    }
		}
	    }
	}
    }
  printf("Chain         Temperature\n");
  printf("-------------------------\n");
  for (i = options->heated_chains - 1; i >= 0; --i)
    printf("%4li %20.5f\n", i+1, options->heat[i]);
  printf("\n\nIs this correct? [YES or NO]\n===> ");
  fflush(stdout); FGETS(input, LINESIZE, stdin);
  if ('y' != tolower(input[0]))
    read_heatvalues(options);
}

void
menuRandom(MYREAL * param, char type) {
  //long            check;
  if(type == 'N')
    {
      printf
	("Specify a MEAN and a STANDARD DEVIATION\nThis will be used to generate a random value from a Normal distribution\n");
    }
  else
    {
      printf
	("Specify a MINIMUM and a MAXIMUM\nThis will be used to generate a random value from a Uniform distribution\n");
    }
#ifdef USE_MYREAL_FLOAT
  scanf("%f%f", &param[0], &param[1]);
#else
 scanf("%lf%lf", &param[0], &param[1]);
#endif
}

void
start_tree_method(option_fmt * options) {
  char            input[LINESIZE];
  char            compstring[LINESIZE];
  printf("  Start genealogy is created with\n");
  printf("      (a)utomatic\n");
  if              (strchr(SEQUENCETYPES, options->datatype)) {
    printf("      (u)sertree\n");
    printf("      (d)istancematrix\n");
  }
  printf("      (r)andom\n\n===> ");
  strcpy(compstring, (strchr(SEQUENCETYPES, options->datatype) ?
		      "daur" : "dar"));
  do {
    fflush(stdout); FGETS(input, LINESIZE, stdin);
  }
  while (strchr(compstring, (int) (lowercase(input[0]))) == NULL);
  switch (lowercase(input[0])) {
  case 'u':
    options->usertree = TRUE;
    options->dist = FALSE;
    break;
  case 'd':
    options->usertree = FALSE;
    options->dist = TRUE;
    options->randomtree = FALSE;
    break;
  case 'a':
    options->usertree = FALSE;
    options->dist = FALSE;
    options->randomtree = FALSE;
    break;
  case 'r':
    options->usertree = FALSE;
    options->dist = FALSE;
    options->randomtree = TRUE;
    break;
  }
}

char * submodeltype(int type)
{
  switch(type)
    {
    case MULTISTEP:
      return "Multistep method";
    case SINGLESTEP:
      return "Singlestep method";
    default:
      return "Singlestep method";
    }
}

/// \brief Menu for exponential proposal
/// Menu for exponential proposal
void menu_submodel(option_fmt *options)
{
  char input[LINESIZE];
  char text[LINESIZE];
  char val;
  double tune = 0. ;
  double pinc = 0.5;
  long count = 0;
  do {
    if(options->msat_option == MULTISTEP)
      {
	sprintf(text,"MULTISTEP (tune=%f, p_increase=%f)",options->msat_tuning[0], options->msat_tuning[1]);
      }
    else
      {
	sprintf(text,"SINGLESTEP");
      }
    printf("\nSingle step or Multi step model? [ENTER either S or M]\n");
    printf("Single Step model is the standard stepwise mutation model\n");
    printf("The multistep model is between the infinite model and and the\n");
    printf("singlestep model using two additional parameters:\n");
    printf("- chance of increasing repeat length\n");
    printf("- tuning between single step (tune=0) and infinite allele (tune=1)\n");
    printf("[Current value: %s]\n===> ",text);
    fflush(stdout); FGETS(input, LINESIZE, stdin);
    if (input[0] == '\0')
      break;
    count = sscanf(input, "%c", &val);
    if(uppercase(val) == 'S')
      options->msat_option = SINGLESTEP;
    else
      {
	options->msat_option = MULTISTEP;
	do
	  {
	    printf("MULTISTEP model\nGive two numbers separated by a space:\nFor the tune (range 0 to 1) and\nfor the probability of repeat number increase (0 to 2/3)\n===> ");
	    fflush(stdout); FGETS(input, LINESIZE, stdin);
	    count = sscanf(input, "%lf%lf", &tune, &pinc);
	    options->msat_tuning[0] = (MYREAL) tune;
	    options->msat_tuning[1] = (MYREAL) pinc;
	  } while (count != 2 && count > 0);
	count = 1;
      }
  } while (count != 1 && count > 0);
}


void
start_data_method(option_fmt * options) {
  char            input[LINESIZE];
  do {
    printf("  (a)llele model\n");
    printf("  (m)icrosatellite model [SLOW! Ladder model: %s]\n", submodeltype(options->msat_option));
    printf("  (b)rownian microsatellite model [Brownian motion model]\n");
    printf("  (s)equence model\n");

    printf("  (n)ucleotide polymorphism (SNP)\n");
    //printf("  (u)nlinked nucleotide polymorphism (SNP) using a PANEL\n");
    printf("  (g)enealogy summaries\n\n===> ");
    fflush(stdout); 
    FGETS(input, LINESIZE, stdin);
  }
  while           (strchr("ambnsguf", (int) (lowercase(input[0]))) == NULL);
  options->datatype = input[0];

  if(options->datatype == 'm')
    {
      menu_submodel(options);
    }
  if (!strchr(SEQUENCETYPES, options->datatype))
    options->usertree = FALSE;
}

long            get_prior(char *input) {
  switch (toupper(input[0])) {
    //case 'S':
    //return SLICE;
    //break;
  case 'M':
    return MULTPRIOR;
    break;
  case 'E':
    return EXPPRIOR;
    break;
  case 'W':
    return WEXPPRIOR;
    break;
  case 'U':
    return UNIFORMPRIOR;
  default:
    return UNIFORMPRIOR;
    break;
  }
}

///
/// prints type of proposal in the menu
void set_proposal(char *output, boolean *proposal, boolean without_rate)
{
  int i;
  int l=0;
  const char types[][LINESIZE]={"Theta","Mig","Rate"};
  
  for(i=THETAPRIOR; i <= RATEPRIOR; i++)
    {
      if(without_rate && i==RATEPRIOR)
	continue;
      switch (proposal[i]) {
      case TRUE:
	l+=sprintf(output+l, " %s:Slice", types[i]);
	break;
      case FALSE:
      default:
	l+=sprintf(output+l, " %s:MH", types[i]);
	break;
	//      default:
	//l+=sprintf(output+l, " %s:?", types[i]);
	//break;
      }
    }
}

///
/// prints type of prior in the menu
void            set_prior(char *output, int *prior, boolean without_rate)
{
  int i;
  int l=0;
  const char types[][LINESIZE]={"Theta","Mig","Rate"};
  for(i=THETAPRIOR; i <= RATEPRIOR; i++)
    {
      if(without_rate && i==RATEPRIOR)
	continue;

      switch (prior[i]) {
	//case SLICE:
	//l+=sprintf(output+l, " Slice");
	//break;
      case MULTPRIOR:
	l+=sprintf(output+l, " %s:Mult.", types[i]);
	break;
      case EXPPRIOR:
	l+=sprintf(output+l, " %s:Exp.", types[i]);
	break;
      case WEXPPRIOR:
	l+=sprintf(output+l, " %s:WExp.", types[i]);
	break;
      case UNIFORMPRIOR:
      default:
	l+=sprintf(output+l, " %s:Unif.", types[i]);
	break;
      }
    }
}

void set_localities_string(char *loc, option_fmt *options)
{
  long z;
  long i;
  if(options->newpops_numalloc>1)
    {
      z=0;
      for(i=0; i<options->newpops_numalloc; i++)
	{
	  z += sprintf(loc+z,"%li,",options->newpops[i]);
	}
    }
  else
    {
      strncpy(loc,"default",7);
    }
}
