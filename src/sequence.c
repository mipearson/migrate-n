/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 S E Q U E N C E S   R O U T I N E S 
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
 Copyright 2002 Peter Beerli and Joseph Felsenstein
 
  This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
 $Id: sequence.c 1860 2011-04-03 15:17:11Z beerli $
 
-------------------------------------------------------*/
/* \file sequence.c

*/
#include "migration.h"
#include "sighandler.h"
#include "data.h"
#include "tools.h"
#include "migrate_mpi.h"
#include "watterson.h"
#include "tree.h"

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
/* prototypes ------------------------------------------- */
void make_sequences (world_fmt * world, option_fmt * options, data_fmt * data,
                     long locus);
void init_sequences (world_fmt * world, option_fmt * options, data_fmt * data,
                     long locus);
void init_sequences2 (world_fmt * world, seqmodel_fmt * seq, long locus);
void initratio (option_fmt * options);
void initfreqs (MYREAL *freqa, MYREAL *freqc, MYREAL *freqg, MYREAL *freqt);
void initcatn (long *categs);
boolean initcategs (long categs, MYREAL *rate, MYREAL *probcat);
void initprobcat (long categs, MYREAL *probsum, MYREAL *probcat);
void init_tbl (world_fmt * world, long locus);
void print_weights (FILE * outfile, world_fmt * world, option_fmt * options,
                    long locus);
void print_tbl (FILE * outfile, world_fmt * world, option_fmt * options,
                long locus);
MYREAL treelike_seq (world_fmt * world, long locus);
MYREAL treelike_snp (world_fmt * world, long locus);
/*private functions */
void getbasefreqs (option_fmt * options, seqmodel_fmt * seq, long locus);
void empiricalfreqs (world_fmt * world, option_fmt * options,
                     seqmodel_fmt * seq, long locus);
void makeweights (world_fmt * world, data_fmt * data, option_fmt *options, long locus);
void makevalues_seq (world_fmt * world, option_fmt * options, data_fmt * data,
                     long locus);
void make_invarsites (world_fmt * world, data_fmt * data, long locus);
void make_invarsites_unlinked (world_fmt * world, data_fmt * data,
                               long locus);
void sitecombine2 (world_fmt * world, data_fmt * data, option_fmt *options, long sites,
                   long locus);
void sitesort2 (world_fmt * world, data_fmt * data, option_fmt *options, long sites, long locus);
void sitescrunch2 (world_fmt * world, long sites, long i, long j, long locus);
void inputoptions (world_fmt * world, option_fmt * options, data_fmt * data,
                   long locus);
void inputweights (world_fmt * world, data_fmt * data, long chars);
void inputcategs (long a, long b, world_fmt * world, option_fmt * options,
                  data_fmt * data);
void initlambda (option_fmt * options);
void printweights (FILE * outfile, world_fmt * world, option_fmt * options,
                   short inc, long chars, short *weight, char *letters);
void print_seqfreqs (FILE * outfile, world_fmt * world, option_fmt * options);

void snp_invariants (contribarr invariants, world_fmt *world, long locus, phenotype x1);

void makevalues_snp (world_fmt * world, option_fmt * options, data_fmt * data,
                     long locus);

void init_sequences_aliases (world_fmt * world, option_fmt * options,
                             data_fmt * data, long locus);


void check_basefreq (option_fmt * options);

extern void swap (void *, void *);



/* allocation things */
void
init_sequences (world_fmt * world, option_fmt * options, data_fmt * data,
                long locus)
{

    long sites = world->data->seq[0]->sites[locus];
    if (world->options->datatype == 'u')
        sites *= 5;
    if (!world->data->seq[0]->done || world->data->seq[0]->alias == NULL)
    {
        world->data->seq[0]->done = TRUE;
        world->data->seq[0]->alias = (long *) mycalloc (1, sizeof (long) * sites);
        world->data->seq[0]->ally = (long *) mycalloc (1, sizeof (long) * sites);
        world->data->seq[0]->savealiasweight = (long *) mycalloc (sites, sizeof (long));
        world->data->seq[0]->aliasweight = (long *) mycalloc (1, sizeof (long) * sites);
        world->data->seq[0]->location = (long *) mycalloc (1, sizeof (long) * sites);
        world->data->seq[0]->category = (long *) mycalloc (1, sizeof (long) * sites);
        world->data->seq[0]->weight = (short *) mycalloc (1, sizeof (short) * sites);
    }
    else
    {
        world->data->seq[0]->alias =
            (long *) myrealloc (world->data->seq[0]->alias, sizeof (long) * sites);
        world->data->seq[0]->ally =
            (long *) myrealloc (world->data->seq[0]->ally, sizeof (long) * sites);
        world->data->seq[0]->savealiasweight =
            (long *) myrealloc (world->data->seq[0]->savealiasweight,
                              sizeof (long) * sites);
        world->data->seq[0]->aliasweight =
            (long *) myrealloc (world->data->seq[0]->aliasweight,
                              sizeof (long) * sites);
        world->data->seq[0]->location =
            (long *) myrealloc (world->data->seq[0]->location, sizeof (long) * sites);
        world->data->seq[0]->category =
            (long *) myrealloc (world->data->seq[0]->category, sizeof (long) * sites);
        world->data->seq[0]->weight =
            (short *) myrealloc (world->data->seq[0]->weight, sizeof (short) * sites);
    }
}

void
init_sequences_aliases (world_fmt * world, option_fmt * options,
                        data_fmt * data, long locus)
{
    inputoptions (world, options, data, locus);
    if (!options->freqsfrom)
        getbasefreqs (options, world->data->seq[0], locus);
    makeweights (world, data, options, locus);
}


/* menu material ----------------------------------------- */
void
initratio (option_fmt * options)
{
    long z = 0;
    char *tmp;
    char input[LINESIZE];
    printf
    ("Transition/transversion ratio?\nEnter a value for each locus, spaced by blanks or commas\n[Minium value > 0.5]");
    FGETS (input, LINESIZE, stdin);
    tmp = strtok (input, " ,\n");
    while (tmp != NULL)
    {
        options->ttratio[z++] = atof (tmp);
        tmp = strtok (NULL, " ,;\n");
        options->ttratio =
            (MYREAL *) myrealloc (options->ttratio, sizeof (MYREAL) * (z + 1));
        options->ttratio[z] = 0.0;
    }
}

void
initfreqs (MYREAL *freqa, MYREAL *freqc, MYREAL *freqg, MYREAL *freqt)
{
    char input[LINESIZE];
    int scanned;
    MYREAL summ = 0;

    printf
    ("Base frequencies for A, C, G, T/U\n (use blanks to separate, if all are equal use a sign [=])?\n");
    for (;;)
    {
        FGETS (input, LINESIZE, stdin);
        if (input[0] == '=')
        {
            scanned = 4;
            *freqa = *freqc = *freqg = *freqt = 0.25;
        }
        else
#ifdef USE_MYREAL_FLOAT
            scanned = sscanf (input, "%f%f%f%f%*[^\n]", freqa, freqc, freqg, freqt);
#else
        scanned = sscanf (input, "%lf%lf%lf%lf%*[^\n]", freqa, freqc, freqg, freqt);
#endif
        if (scanned == 4)
            break;
        else
            printf ("Please enter exactly 4 values.\n");
    };
    // adjust frequencies to a total of 1
    summ = *freqa + *freqc + *freqg + *freqt;
    if (summ != 1.0)
        printf ("Frequency values were adjusted to add up to 1.0.\n");
    *freqa /= summ;
    *freqc /= summ;
    *freqg /= summ;
    *freqt /= summ;
    printf ("Nucleotide frequencies: A=%f, C=%f, G=%f, T=%f\n", *freqa, *freqc,
            *freqg, *freqt);
}


void
initcatn (long *categs)
{    /* initialize category number */

    do
    {
        printf ("Number of categories (1-%d)?\n", MAXCATEGS);
        scanf ("%ld%*[^\n]", categs);
        getchar ();
    }
    while (*categs > MAXCATEGS || *categs < 1);
}


boolean
initcategs (long categs, MYREAL *rate, MYREAL *probcat)
{    /* initialize rate categories */
    long i;
    char input[LINESIZE];
    char rest[LINESIZE];
    int scanned;
    boolean done;

    for (;;)
    {
        printf
	  ("Either enter the Shape parameter alpha for Gamma deviated rates\n*OR* enter the rates for each category (use a space to separate)\n===>");fflush(stdout);
        FGETS (input, LINESIZE, stdin);
        done = TRUE;
        if (count_words (input) == 1)
        {
            gamma_rates (rate, probcat, categs, input);
            return TRUE;
        }
        for (i = 0; i < categs; i++)
        {
#ifdef USE_MYREAL_FLOAT
            scanned = sscanf (input, "%f %[^\n]", &rate[i], rest);
#else
            scanned = sscanf (input, "%lf %[^\n]", &rate[i], rest);
#endif
            if ((scanned < 2 && i < (categs - 1))
                    || (scanned < 1 && i == (categs - 1)))
            {
                printf ("Please enter exactly %ld values.\n", categs);
                done = FALSE;
                break;
            }
            strcpy (input, rest);
        }
        if (done)
            break;
    }
    return FALSE;
}

void
initprobcat (long categs, MYREAL *probsum, MYREAL *probcat)
{
    long i;
    boolean done;
    char input[LINESIZE];
    char rest[LINESIZE];
    int scanned;

    do
    {
        printf ("Probability for each category?");
        printf (" (use a space to separate)\n");
        FGETS (input, LINESIZE, stdin);
        done = TRUE;
        for (i = 0; i < categs; i++)
        {
#ifdef USE_MYREAL_FLOAT
            scanned = sscanf (input, "%f %[^\n]", &probcat[i], rest);
#else
            scanned = sscanf (input, "%lf %[^\n]", &probcat[i], rest);
#endif
            if ((scanned < 2 && i < (categs - 1))
                    || (scanned < 1 && i == (categs - 1)))
            {
                done = FALSE;
                printf ("Please enter exactly %ld values.\n", categs);
                break;
            }
            strcpy (input, rest);
        }
        if (!done)
            continue;
        *probsum = 0.0;
        for (i = 0; i < categs; i++)
            *probsum += probcat[i];

        if (fabs (1.0 - (*probsum)) > 0.001)
        {
            for (i = 0; i < categs; i++)
                probcat[i] /= *probsum;
            printf ("Probabilities were adjusted to add up to one\n");
            for (i = 0; i < categs; i++)
                printf ("  %li> %f\n", i + 1, probcat[i]);
            printf ("\n\n");
        }
    }
    while (!done);
}

///
/// constrains the arbitrary site rate variation to an average of 1.0
void constrain_rates(long categs, MYREAL *rate, MYREAL *probcat)
{
  char input[LINESIZE];
  long i;
  MYREAL mean;
  printf("By default the rates will be constrained to an have an average rate of 1.0\nPress return if that is OK (preferred option) or enter the word RAW\n===>"); fflush(stdout);
  FGETS (input, LINESIZE, stdin);
  if(input[0] == '\0')
    {
      mean = 0.0;
      for (i = 0; i < categs; i++)
	{
	  mean += probcat[i] * rate[i];
	}
      for (i = 0; i < categs; i++)
	{
	  rate[i] /= mean;
	}
    }
}

/*data read material ===================================== */
///
/// read sequence data and linked SNP data 
void
make_sequences (world_fmt * world, option_fmt * options, data_fmt * data,
                long locus)
{
    if (world->sumtips==0)
      {
	data->skiploci[locus] = TRUE;
        world->data->skiploci[locus] = TRUE;
        world->skipped += 1;
      }

  makevalues_seq (world, options, data, locus);
  if (options->freqsfrom)
    {
      empiricalfreqs (world, options, world->data->seq[0], locus);
      getbasefreqs (options, world->data->seq[0], locus);
    } 
#ifdef BEAGLE
  if(world->mutationmodels[world->sublocistarts[locus]].basefreqs==NULL)
    world->mutationmodels[world->sublocistarts[locus]].basefreqs = (double*) calloc(4,sizeof(double));
  world->mutationmodels[world->sublocistarts[locus]].basefreqs[0] =  world->data->seq[0]->basefrequencies[0];
  world->mutationmodels[world->sublocistarts[locus]].basefreqs[1] =  world->data->seq[0]->basefrequencies[1];
  world->mutationmodels[world->sublocistarts[locus]].basefreqs[2] =  world->data->seq[0]->basefrequencies[2];
  world->mutationmodels[world->sublocistarts[locus]].basefreqs[3] =  world->data->seq[0]->basefrequencies[3];
  world->mutationmodels[world->sublocistarts[locus]].numstates=4;
#endif
  
}

/// 
/// read unlinked SNP data (each locus must be only 1 site)
void
make_snp (world_fmt * world, option_fmt * options, data_fmt * data,
          long locus)
{
  makevalues_snp (world, options, data, locus);
  /*if (options->freqsfrom)
    {
    empiricalfreqs (world, world->data->seq, locus); */
  options->freqsfrom = FALSE;
  getbasefreqs (options, world->data->seq[0], locus);
  /* }
   */
}

/* private functions================================== */
void
getbasefreqs (option_fmt * options, seqmodel_fmt * seq, long locus)
{
  long l;

  register MYREAL freqa, freqc, freqg, freqt, freqr, freqy;
  register MYREAL /*freqar, freqcy,*/ freqgr, freqty;

  MYREAL aa, bb;
  if (locus == 0)
    seq->ttratio = options->ttratio[0];
  else
    {
      for (l = 1; l <= locus; l++)
        {
	  if (options->ttratio[l] == 0.0)
            {
	      seq->ttratio = options->ttratio[l - 1];
	      break;
            }
	  seq->ttratio = options->ttratio[l];
        }
      if (l > locus)
	seq->ttratio = options->ttratio[locus];
    }
  check_basefreq (options);
  
  seq->basefrequencies[NUC_A] = options->freqa;
  seq->basefrequencies[NUC_C] = options->freqc;
  seq->basefrequencies[NUC_G] = options->freqg;
  seq->basefrequencies[NUC_T] = options->freqt;    
  freqa = seq->basefrequencies[NUC_A];
  freqc = seq->basefrequencies[NUC_C];
  freqg = seq->basefrequencies[NUC_G];
  freqt = seq->basefrequencies[NUC_T];
  seq->basefrequencies[NUC_R] = freqa + freqg;
  seq->basefrequencies[NUC_Y] = freqc + freqt;
  freqr = seq->basefrequencies[NUC_R];
  freqy = seq->basefrequencies[NUC_Y];
  seq->basefrequencies[NUC_AR]= freqa / freqr;
  seq->basefrequencies[NUC_CY]= freqc / freqy;
  seq->basefrequencies[NUC_GR]= freqg / freqr;
  seq->basefrequencies[NUC_TY]= freqt / freqy;
  //freqar = seq->basefrequencies[NUC_AR];
  //freqcy = seq->basefrequencies[NUC_CY];
  freqgr = seq->basefrequencies[NUC_GR];
  freqty = seq->basefrequencies[NUC_TY];

  aa =
    seq->ttratio * (freqr) * (freqy) - freqa * freqg -
    freqc * freqt;
  bb = freqa * (freqgr) + freqc * (freqty);
  seq->xi = aa / (aa + bb);
  seq->xv = 1.0 - seq->xi;
  if (seq->xi <= 0.0)
    {
      warning ("This transition/transversion ratio (%f)\n",seq->ttratio);
      warning ("is impossible with these base frequencies (%f, %f, %f, %f)!\n",freqa,freqc,freqg,freqt);
      seq->xi = 0.00001; // do not set this to zero because of the 1/(fracchange=xi*(...))
      seq->xv = 0.99999;
      seq->ttratio =
	(freqa * freqg +
	 freqc * freqt) / ((freqr) * (freqy));
      
      warning (" Transition/transversion parameter reset\n");
      warning ("  so transition/transversion ratio is %10.6f\n\n",
	       (seq->ttratio));
    }
  // use 1/frac as precomputation speed up
  seq->fracchange = 1. / (
			  (seq->xi) * (2. * freqa * (freqgr) +
				       2. * freqc * (freqty)) + (seq->xv) * (1.0 -
										      freqa *
										      freqa -
										      freqc *
										      freqc -
										      freqg *
										      freqg -
										      freqt *
										      freqt));
}

/*===================================================*/

void
makeweights (world_fmt * world, data_fmt * data, option_fmt *options, long locus)
{
    /* make up weights vector to avoid duplicate computations */
    long i;
    seqmodel_fmt *seq = world->data->seq[0];
    world->data->seq[0]->endsite = 1;
    for (i = 0; i < seq->sites[locus]; i++)
    {
        seq->alias[i] = i + 1;
        seq->ally[i] = 0;
        seq->aliasweight[i] = seq->weight[i];
        seq->location[i] = 0;
    }
    sitesort2 (world, data, options, seq->sites[locus], locus);
    sitecombine2 (world, data, options, seq->sites[locus], locus);
    sitescrunch2 (world, seq->sites[locus], 1, 2, locus);
    for (i = 1; i <= seq->sites[locus]; i++)
    {
        if (seq->aliasweight[i - 1] > 0)
            seq->endsite = i;
    }
    for (i = 1; i <= seq->endsite; i++)
    {
        seq->location[seq->alias[i - 1] - 1] = i;
        seq->ally[seq->alias[i - 1] - 1] =
            seq->alias[i - 1];
    }
    init_sequences2 (world, seq, locus);
    memcpy(seq->savealiasweight, seq->aliasweight, sizeof(long) * seq->endsite);
}    /* makeweights */

void
init_sequences2 (world_fmt * world, seqmodel_fmt * seq, long locus)
{
    if (world->contribution == NULL)
      {
        world->contribution =
            (contribarr *) mymalloc ((4 + seq->endsite) * sizeof (contribarr));
	//	printf("%i> temp=%f size=%li: world->contribution allocated\n",myID, world->heat,(4 + seq->endsite) * sizeof (contribarr));
      }
    else
      {
        world->contribution =
	  (contribarr *) myrealloc (world->contribution,
                                    (4 + seq->endsite) * sizeof (contribarr));
	//	printf("%i> temp=%f size=%li: world->contribution REallocated\n",myID, world->heat,(4 + seq->endsite) * sizeof (contribarr));
      }
}


void set_nucleotide(MYREAL *treedata, const char nucleotide, const MYREAL seqerr)
{
  long b;
  const MYREAL seqerr3 = seqerr / 3.;
  const MYREAL seqerr23 = 2. * seqerr / 3.;
  const MYREAL oneseqerr = 1.0 - seqerr;
  const MYREAL oneseqerr23 = 1.0 - seqerr23;
  const MYREAL oneseqerr3 = 1.0 - seqerr3;

  for (b = 0; b < 4; b++)
    treedata[b] = 0.0 + seqerr3;
  switch (nucleotide)
    {
    case 'A':
      treedata[NUC_A] = oneseqerr;
      break;
      
    case 'C':
      treedata[NUC_C] = oneseqerr;
      break;
      
    case 'G':
      treedata[NUC_G] = oneseqerr;
      break;
      
    case 'T':
      treedata[NUC_T] = oneseqerr;
      break;
      
    case 'U':
      treedata[NUC_T] = oneseqerr;
      break;
      
    case 'M':
      treedata[NUC_A] = oneseqerr;
      treedata[NUC_C] = oneseqerr;
      break;
      
    case 'R':
      treedata[NUC_A] = oneseqerr23;
      treedata[NUC_G] = oneseqerr23;
      break;
      
    case 'W':
      treedata[NUC_A] = oneseqerr23;
      treedata[NUC_T] = oneseqerr23;
      break;
      
    case 'S':
      treedata[NUC_C] = oneseqerr23;
      treedata[NUC_G] = oneseqerr23;
      break;
      
    case 'Y':
      treedata[NUC_C] = oneseqerr23;
      treedata[NUC_T] = oneseqerr23;
      break;
      
    case 'K':
      treedata[NUC_G] = oneseqerr23;
      treedata[NUC_T] = oneseqerr23;
      break;
      
    case 'B':
      treedata[NUC_A] = seqerr;
      treedata[NUC_C] = oneseqerr3;
      treedata[NUC_G] = oneseqerr3;
      treedata[NUC_T] = oneseqerr3;
      break;
      
    case 'D':
      treedata[NUC_A] = oneseqerr3;
      treedata[NUC_C] = seqerr;
      treedata[NUC_G] = oneseqerr3;
      treedata[NUC_T] = oneseqerr3;
      break;
      
    case 'H':
      treedata[NUC_A] = oneseqerr3;
      treedata[NUC_C] = oneseqerr3;
      treedata[NUC_G] = seqerr;
      treedata[NUC_T] = oneseqerr3;
      break;
      
    case 'V':
      treedata[NUC_A] = oneseqerr3;
      treedata[NUC_C] = oneseqerr3;
      treedata[NUC_G] = oneseqerr3;
      treedata[NUC_T] = seqerr;
      break;
      
    case 'N':
      for (b = 0; b < 4; b++)
	treedata[b] = 1.0;
      break;
      
    case 'X':
      for (b = 0; b < 4; b++)
	treedata[b] = 1.0;
      break;
      
    case '?':
      for (b = 0; b < 4; b++)
	treedata[b] = 1.0;
      break;
      
    case 'O':
      for (b = 0; b < 4; b++)
	treedata[b] = 1.0;
      break;
#ifdef GAP                            
    case '-':
      treedata[4] = 1.0;
      break;
#else
    case '-':
      for (b = 0; b < 4; b++)
	treedata[b] = 1.0;
      break;
#endif
    }
}


void
makevalues_seq (world_fmt * world, option_fmt * options, data_fmt * data,
                long locus)
{
  seqmodel_fmt *seq = world->data->seq[0];
  long ii, ind, j, k, l, pop, top;
  long z=0;
  node **treenode = world->nodep;
  MYREAL seqerr = options->seqerror;
  long zpop;
  // we need to recalculate the sumtips
  world->sumtips = 0;

  for (pop = 0; pop < data->numpop; pop++)
    {
      zpop = 0;
      top = max_shuffled_individuals(options, data, pop, locus);
      for (ii = 0; ii < top; ii++)
	{
	  ind = data->shuffled[pop][locus][ii];
	  zpop++;
	  //	  if(options->include_unknowns)
	  // {
	      if (!options->usertree)
		{
		  allocate_tip (world, options, &(treenode[z]), pop, locus, z, ind,
				data->indnames[pop][ind]);
		  //treenode[z]->tyme = find_tipdate(treenode[z]->nayme, pop,world);
		  world->data->sampledates[pop][locus][ind].id = treenode[z]->id;
		  z++;
		}	    
	      for (k = 0; k < seq->endsite; k++)
		{
		  j = seq->alias[k] - 1;
		  for (l = 0; l < world->options->rcategs; l++)
		    {                    
		      set_nucleotide(treenode[z-1]->x.s[k][l],data->yy[pop][ind][locus][0][j], seqerr);
		    }
		}
        }
      data->numalleles[pop][locus] = zpop;
      world->sumtips += zpop;
    }
}



void
makevalues_snp (world_fmt * world, option_fmt * options, data_fmt * data,
                long locus)
{
  long i, ii, j, k, l, pop;
  long b, f;
  long ind;
  seqmodel_fmt *seq = world->data->seq[0];
  node **treenode = world->nodep;
  f = seq->addon;
  for (k = 0; k < seq->endsite * (f + 1);
       k += (1 + seq->addon))
    {
      j = seq->alias[k / (f + 1)] - 1;
      i = -1;
      for (pop = 0; pop < data->numpop; pop++)
	{
	  for (ii = 0; ii < data->numalleles[pop][locus]; ii++)
	    {
	      ind = data->shuffled[pop][locus][ii];
	      i++;
	      if (k==0 && !options->usertree)
		{
		  strcpy (treenode[i]->nayme, data->indnames[pop][ind][locus]);
		  treenode[i]->tyme = find_tipdate(treenode[i]->nayme, pop,world);
		}
	      for (l = 0; l < world->options->rcategs; l++)
		{
		  if (pop == 0)
		    {
		      //                        treenode[i]->x.s[k][l] = vec_float_one();
		      // reset all values for the snp columns
		      //  treenode[i]->x.s[k][l] = vec_zero();
		      for (b = 1; b <= 4; b++)
			{
			  memset (treenode[i]->x.s[k + b][l], 0,
				  sizeof (MYREAL) * 4);
			  treenode[i]->x.s[k + b][l][b - 1] = 1.0;
			}
		      continue;
		    }
		  else
		      {
                        for (b = 0; b <= 4; b++)
			  memset (treenode[i]->x.s[k + b][l], 0,
				  sizeof (MYREAL) * 4);
		      }
		    
                    switch (data->yy[pop][ind][locus][0][j])
		      {
		      case 'A':
                        treenode[i]->x.s[k][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_A] = 1.0;
                        break;
			
                    case 'C':
                        treenode[i]->x.s[k][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_C] = 1.0;
                        break;

                    case 'G':
                        treenode[i]->x.s[k][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_G] = 1.0;
                        break;
                    case 'T':
                    case 'U':
                        treenode[i]->x.s[k][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_T] = 1.0;
                        break;

                    case 'M':
                        treenode[i]->x.s[k][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_C] = 1.0;
                        break;

                    case 'R':
                        treenode[i]->x.s[k][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_G] = 1.0;
                        break;

                    case 'W':
                        treenode[i]->x.s[k][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_T] = 1.0;
                        break;

                    case 'S':
                        treenode[i]->x.s[k][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_G] = 1.0;
                        break;

                    case 'Y':
                        treenode[i]->x.s[k][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_T] = 1.0;
                        break;

                    case 'K':
                        treenode[i]->x.s[k][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_T] = 1.0;
                        break;

                    case 'B':
                        treenode[i]->x.s[k][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_T] = 1.0;
                        break;

                    case 'D':
                        treenode[i]->x.s[k][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_T] = 1.0;
                        break;

                    case 'H':
                        treenode[i]->x.s[k][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_T] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_T] = 1.0;
                        break;

                    case 'V':
                        treenode[i]->x.s[k][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_A] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_C] = 1.0;
                        treenode[i]->x.s[k + 1][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 2][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 3][l][NUC_G] = 1.0;
                        treenode[i]->x.s[k + 4][l][NUC_G] = 1.0;
                        break;

                    case 'N':
                    case 'X':
                    case '-':
                    case 'O':
                    case '?':
                        for (b = 0; b < 4; b++)
                        {
                            treenode[i]->x.s[k][l][b] = 1.0;
                            treenode[i]->x.s[k + 1][l][b] = 1.0;
                            treenode[i]->x.s[k + 2][l][b] = 1.0;
                            treenode[i]->x.s[k + 3][l][b] = 1.0;
                            treenode[i]->x.s[k + 4][l][b] = 1.0;
                        }
                        break;
                    default:
                        break;

                    }
                }
            }
        }
    }
}


void
make_invarsites (world_fmt * world, data_fmt * data, long locus)
{
    long i, k, z, ii, pop, l;
    long endsite = world->data->seq[0]->endsite;
    node **treenode = world->nodep;
    i = -1;
    for (pop = 0; pop < data->numpop; pop++)
    {
        for (ii = 0; ii < data->numalleles[pop][locus]; ii++)
        {
            i++;
            z = 0;
            for (k = endsite; k < endsite + 4; k++)
            {
                for (l = 0; l < world->options->rcategs; l++)
                {
                    memset (treenode[i]->x.s[k][l], 0, sizeof (MYREAL) * 4);
                    treenode[i]->x.s[k][l][z] = 1.0;
                }
                z++;
            }
        }
    }
}

void
make_invarsites_unlinked (world_fmt * world, data_fmt * data, long locus)
{
    long i, k, z, ii, pop, l;
    long endsite = world->data->seq[0]->endsite;
    node **treenode = world->nodep;
    i = -1;
    for (pop = 0; pop < data->numpop; pop++)
    {
        for (ii = 0; ii < data->numalleles[pop][locus]; ii++)
        {
            i++;
            z = 0;
            if (pop == PANEL)
            {
                for (k = endsite * 5; k < endsite * 5 + 4; k++)
                {
                    for (l = 0; l < world->options->rcategs; l++)
                    {
#ifdef ALTIVEC 
                        memset (treenode[i]->x.s[k][l].f, 0, sizeof (MYREAL) * 4);
                        treenode[i]->x.s[k][l].f[z] = 1.0;
#else
                        memset (treenode[i]->x.s[k][l], 0, sizeof (MYREAL) * 4);
                        treenode[i]->x.s[k][l][z] = 1.0;
#endif
                    }
                    z++;
                }
            }
            else
            {
                for (k = endsite * 5; k < endsite * 5 + 4; k++)
                {
                    for (l = 0; l < world->options->rcategs; l++)
                    {
#ifdef ALTIVEC 
                        treenode[i]->x.s[k][l].f[0] = 1.0; //putting ? for data
                        treenode[i]->x.s[k][l].f[1] = 1.0;
                        treenode[i]->x.s[k][l].f[2] = 1.0;
                        treenode[i]->x.s[k][l].f[3] = 1.0;
#else
                        treenode[i]->x.s[k][l][0] = 1.0; //putting ? for data
                        treenode[i]->x.s[k][l][1] = 1.0;
                        treenode[i]->x.s[k][l][2] = 1.0;
                        treenode[i]->x.s[k][l][3] = 1.0;
#endif
                    }
                }
            }
        }
    }
}

void
sitesort2 (world_fmt * world, data_fmt * data, option_fmt *options, long sites, long locus)
{
  long gap, i, i1, j, jj, jg, k, kk, kkk, itemp, pop, z = 0;
    boolean flip, tied, samewt;
    seqmodel_fmt *seq;
    long *tempsum, *temppop;
    long numind;
    MYREAL a1,a2;
    tempsum = (long *) mycalloc (1, sizeof (long) * data->numpop);
    temppop = (long *) mycalloc (1, sizeof (long) * data->numpop);
    for (i = 0; i < data->numpop; i++)
    {
      if(options->randomsubset > 0)
	{
	  numind = (options->randomsubset < data->numind[i][locus] ? options->randomsubset : data->numind[i][locus]);
	}
      else
	{
	  numind = data->numind[i][locus] ;
	}
      if (numind > 0)
        {
	  temppop[z] = i;
	  if (z == 0)
	    tempsum[z] = numind;
	  else
	    tempsum[z] = tempsum[z - 1] + numind;
            z++;
        }
    }
    seq = world->data->seq[0];
    gap = sites / 2;
    while (gap > 0)
    {
        for (i = gap + 1; i <= sites; i++)
        {
            j = i - gap;
            flip = TRUE;
            while (j > 0 && flip)
            {
                jj = seq->alias[j - 1];
                jg = seq->alias[j + gap - 1];
                samewt = ((seq->weight[jj - 1] != 0)
                          && (seq->weight[jg - 1] != 0))
                         || ((seq->weight[jj - 1] == 0) && (seq->weight[jg - 1] == 0));
                tied = samewt
                       && (seq->category[jj - 1] == seq->category[jg - 1]);
                flip = ((!samewt) && (seq->weight[jj - 1] == 0))
                       || (samewt
                           && (seq->category[jj - 1] > seq->category[jg - 1]));
                k = 0;
                pop = 0;
                kk = -1;
                while (k < world->sumtips && tied)
                {
                    if (k == tempsum[pop])
                    {
                        kk = 0;
                        pop++;
                    }
                    else
                    {
                        kk++;
                    }
		    i1 = temppop[pop];
		    kkk = data->shuffled[i1][locus][kk];
		    a1 = data->yy[i1][kkk][locus][0][jj - 1];
		    a2 = data->yy[i1][kkk][locus][0][jg - 1];
                    flip =
                        (a1 > a2);
                    tied = (tied
                            && a1 == a2);
                    k++;
                }
                if (!flip)
                    break;
                itemp = seq->alias[j - 1];
                seq->alias[j - 1] = seq->alias[j + gap - 1];
                seq->alias[j + gap - 1] = itemp;
                itemp = seq->aliasweight[j - 1];
                seq->aliasweight[j - 1] = seq->aliasweight[j + gap - 1];
                seq->aliasweight[j + gap - 1] = itemp;
                j -= gap;
            }
        }
        gap /= 2;
    }
    myfree(tempsum);
    myfree(temppop);
}    /* sitesort2 */


void
sitecombine2 (world_fmt * world, data_fmt * data, option_fmt *options, long sites, long locus)
{
    long i, j, k, kk, pop, z = 0;
    boolean tied, samewt;
    seqmodel_fmt *seq;
    long *tempsum, *temppop;
    long numind;
    long kkk;
    tempsum = (long *) mycalloc (1, sizeof (long) * data->numpop);
    temppop = (long *) mycalloc (1, sizeof (long) * data->numpop);
    if(options->randomsubset > 0)
      {
	tempsum[0] = (options->randomsubset < data->numind[0][locus] ? options->randomsubset : data->numind[0][locus]);
      }
    else
      {
	tempsum[0] = data->numind[0][locus];
      }
    
    for (i = 0; i < data->numpop; i++)
    {
      numind = ((options->randomsubset > 0) && (options->randomsubset < data->numind[0][locus])) ? options->randomsubset : data->numind[i][locus];
      if (numind > 0)
        {
	  temppop[z] = i;
	  if (z == 0)
	    tempsum[z] = numind;
	  else
                tempsum[z] = tempsum[z - 1] + numind;
            z++;
        }
    }

    seq = world->data->seq[0];
    i = 1;
    while (i < sites)
    {
        j = i + 1;
        tied = TRUE;
        while (j <= sites && tied)
        {
            samewt = ((seq->aliasweight[i - 1] != 0)
                      && (seq->aliasweight[j - 1] != 0))
                     || ((seq->aliasweight[i - 1] == 0)
                         && (seq->aliasweight[j - 1] == 0));
            tied = samewt
                   && (seq->category[seq->alias[i - 1] - 1] ==
                       seq->category[seq->alias[j - 1] - 1]);
            k = 0;
            pop = 0;
            kk = -1;
            while (k < world->sumtips && tied)
            {
                if (k == tempsum[pop])
                {
                    kk = 0;
                    pop++;
                }
                else
                {
                    kk++;
                }
		kkk = data->shuffled[pop][locus][kk];
                tied = (tied
                        && data->yy[temppop[pop]][kkk][locus][0][seq->
                                                                alias[i - 1] -
                                                                1] ==
                        data->yy[temppop[pop]][kkk][locus][0][seq->alias[j - 1] -
                                                             1]);
                k++;
            }
            if (!tied)
                break;
            seq->aliasweight[i - 1] += seq->aliasweight[j - 1];
            seq->aliasweight[j - 1] = 0;
            seq->ally[seq->alias[j - 1] - 1] = seq->alias[i - 1];
            j++;
        }
        i = j;
    }
    myfree(temppop);
    myfree(tempsum);
}    /* sitecombine2 */


void
sitescrunch2 (world_fmt * world, long sites, long i, long j, long locus)
{
    /* move so positively weighted sites come first */
    /* used by dnainvar, dnaml, dnamlk, & restml */
    long itemp;
    boolean done, found;
    seqmodel_fmt *seq;
    seq = world->data->seq[0];
    done = FALSE;
    while (!done)
    {
        //found = FALSE;
        if (seq->aliasweight[i - 1] > 0)
            i++;
        else
        {
            if (j <= i)
                j = i + 1;
            if (j <= sites)
            {
                //found = FALSE;
                do
                {
                    found = (seq->aliasweight[j - 1] > 0);
                    j++;
                }
                while (!(found || j > sites));
                if (found)
                {
                    j--;
                    itemp = seq->alias[i - 1];
                    seq->alias[i - 1] = seq->alias[j - 1];
                    seq->alias[j - 1] = itemp;
                    itemp = seq->aliasweight[i - 1];
                    seq->aliasweight[i - 1] = seq->aliasweight[j - 1];
                    seq->aliasweight[j - 1] = itemp;
                }
                else
                    done = TRUE;
            }
            else
                done = TRUE;
        }
        done = (done || i >= sites);
    }
}    /* sitescrunch2 */

void
inputoptions (world_fmt * world, option_fmt * options, data_fmt * data,
              long locus)
{
    long i;
    long sites = world->data->seq[0]->sites[locus];
    for (i = 0; i < sites; i++)
        world->data->seq[0]->category[i] = 1;
    for (i = 0; i < sites; i++)
        world->data->seq[0]->weight[i] = 1;
    if (options->weights)
        inputweights (world, data, sites);
    world->data->seq[0]->weightsum = 0;
    for (i = 0; i < sites; i++)
        world->data->seq[0]->weightsum += world->data->seq[0]->weight[i];
    if (world->options->categs > 1)
    {
        inputcategs (0, sites, world, options, data);
    }
}    /* inputoptions */

void
inputweights (world_fmt * world, data_fmt * data, long chars)
{
    /* input the character weights, 0-9 and A-Z for weights 0 - 35 */
    char ch;
    long i;
	char input[1024];

    ch = getc (data->weightfile);
    while (ch == '#')
    {
        FGETS(input, LINESIZE,data->weightfile);
        ch = getc (data->weightfile);
    }
    ungetc (ch, data->weightfile);
    for(i=0;i<chars;i++)
      {
	world->data->seq[0]->weight[i] = 1;
        ch = getc (data->weightfile);
	
	if (isdigit ((int) ch))
	  world->data->seq[0]->weight[i] = ch - '0';
	else 
	  {
	    if (isalpha ((int) ch))
	      {
		ch = uppercase (ch);
		world->data->seq[0]->weight[i] = (short) (ch - 'A' + 10);
	      }
	    else
	      {
		if(isspace((int) ch))
		  {
		    i--;
		    continue;
		  }
		else
		  printf ("ERROR: Bad weight character: %c\n", ch);
		exit (EXIT_FAILURE);
	      }
	  }
      }
}    /* inputweights */

void
inputcategs (long a, long b, world_fmt * world, option_fmt * options,
             data_fmt * data)
{
    /* input the categories, 1-9 */
    char ch;
    long i;
	char input[LINESIZE];
    ch = getc (data->catfile);
    while (ch == '#')
    {
        FGETS(input, LINESIZE,data->catfile);
        ch = getc (data->catfile);
    }
    ungetc (ch, data->catfile);
    fscanf (data->catfile, "%ld", &options->categs);
    options->rate =
        (MYREAL *) myrealloc (options->rate, sizeof (MYREAL) * options->categs);
    for (i = 0; i < options->categs; i++)
    {
#ifdef USE_MYREAL_FLOAT
        fscanf (data->catfile, "%f", &options->rate[i]);
#else
        fscanf (data->catfile, "%lf", &options->rate[i]);
#endif
    }

    for (i = a; i < b; i++)
    {
        if ((ch >= '1') && (ch <= ('0' + world->options->categs)))
            world->data->seq[0]->category[i] = ch - '0';
        else
        {
			if(isspace((int) ch))
			  {
				i--;
				continue;
			  }
			else
			  {
				printf
            ("BAD CATEGORY CHARACTER: %c -- CATEGORIES ARE CURRENTLY 1-%ld\n",
             ch, world->options->categs);
            exit (EXIT_FAILURE);
			  }
        }
    }
}    /* inputcategs */




void
empiricalfreqs (world_fmt * world, option_fmt * options, seqmodel_fmt * seq,
                long locus)
{
    /* Get empirical base frequencies from the data */
    long i, j, k;
    MYREAL summ, suma, sumc, sumg, sumt, w;
    long snps = world->options->datatype == 'u' ? 5 : 1;
    options->freqa = 0.25;
    options->freqc = 0.25;
    options->freqg = 0.25;
    options->freqt = 0.25;
    for (k = 1; k <= 8; k++)
    {
        suma = 0.0;
        sumc = 0.0;
        sumg = 0.0;
        sumt = 0.0;
        for (i = 0; i < world->sumtips; i++)
        {
            for (j = 0; j < seq->endsite * snps; j += snps)
            {
                w = (MYREAL) seq->aliasweight[j / snps];
                summ = (options->freqa) * world->nodep[i]->x.s[j][0][NUC_A] + (options->freqc) * world->nodep[i]->x.s[j][0][NUC_C]
                    + (options->freqg) * world->nodep[i]->x.s[j][0][NUC_G]
                + (options->freqt) * world->nodep[i]->x.s[j][0][NUC_T];
                suma += w * (options->freqa) * world->nodep[i]->x.s[j][0][NUC_A] / summ;
                sumc += w * (options->freqc) * world->nodep[i]->x.s[j][0][NUC_C] / summ;
                sumg += w * (options->freqg) * world->nodep[i]->x.s[j][0][NUC_G] / summ;
                sumt += w * (options->freqt) * world->nodep[i]->x.s[j][0][NUC_T] / summ;

            }
        }
        summ = suma + sumc + sumg + sumt;
        options->freqa = suma / summ;
        options->freqc = sumc / summ;
        options->freqg = sumg / summ;
        options->freqt = sumt / summ;
    }
}    /* empiricalfreqs */


void
initlambda (option_fmt * options)
{
    while (options->lambda <= 1.0)
    {
        printf
        ("Mean block length of sites having the same rate\n (needs to be greater than 1)?\n");
#ifdef USE_MYREAL_FLOAT
        scanf ("%f%*[^\n]", &options->lambda);
#else
        scanf ("%lf%*[^\n]", &options->lambda);
#endif
        getchar ();
    }
    options->lambda = 1.0 / options->lambda;
}



void
init_tbl (world_fmt * world, long locus)
{
    /* Define a lookup table. Precompute values and print them out in tables */
    long i, j;
    MYREAL sumrates;
    long categs = world->options->categs;
    long rcategs = world->options->rcategs;
    world->tbl = (valrec ***) mymalloc (rcategs * sizeof (valrec **));
    for (i = 0; i < rcategs; i++)
    {
        world->tbl[i] = (valrec **) mymalloc (categs * sizeof (valrec *));
        for (j = 0; j < categs; j++)
            world->tbl[i][j] = (valrec *) mymalloc (sizeof (valrec));
    }
    for (i = 0; i < rcategs; i++)
    {
        for (j = 0; j < categs; j++)
        {
            world->tbl[i][j]->rat =
                world->options->rrate[i] * world->options->rate[j];
            world->tbl[i][j]->ratxi =
                world->tbl[i][j]->rat * world->data->seq[0]->xi;
            world->tbl[i][j]->ratxv =
                world->tbl[i][j]->rat * world->data->seq[0]->xv;
        }
    }
    sumrates = 0.0;
    for (i = 0; i < world->data->seq[0]->oldsite; i++)
    {
        for (j = 0; j < rcategs; j++)
            sumrates +=
                world->data->seq[0]->aliasweight[i] * world->options->probcat[j] *
                world->tbl[j][world->data->seq[0]->
                              category[world->data->seq[0]->alias[i] - 1] - 1]->rat;
    }
    sumrates /= (MYREAL) world->data->seq[0]->sites[locus];
    for (i = 0; i < rcategs; i++)
        for (j = 0; j < categs; j++)
        {
            world->tbl[i][j]->rat /= sumrates;
            world->tbl[i][j]->ratxi /= sumrates;
            world->tbl[i][j]->ratxv /= sumrates;
        }
}    /* inittable */

void
print_weights (FILE * outfile, world_fmt * world, option_fmt * options,
               long locus)
{
    if (options->weights)
    {
        if ((options->printdata) || (options->progress && outfile == stdout))
        {
            printweights (outfile, world, options, 0,
                          world->data->seq[0]->sites[locus],
                          world->data->seq[0]->weight, "Sites");
        }
    }
}

void
print_tbl (FILE * outfile, world_fmt * world, option_fmt * options,
           long locus)
{
    long i;

    option_fmt *opt;


    opt = options;
    if (opt->rcategs > 1)
    {
        FPRINTF (outfile, "\nRegion type     Rate of change    Probability\n");
        FPRINTF (outfile, "---------------------------------------------\n");
        for (i = 0; i < opt->rcategs; i++)
            FPRINTF (outfile, "%9ld%16.3f%17.3f\n", i + 1, opt->rrate[i],
                     opt->probcat[i]);
        FPRINTF (outfile,"\n");
        if (opt->autocorr)
            FPRINTF (outfile,
                     "Expected length of a patch of sites having the same rate = %8.3f\n",
                     1. / opt->lambda);
        FPRINTF (outfile,"\n");
    }
    if (opt->categs > 1)
    {
        FPRINTF (outfile, "Site category   Rate of change\n");
        FPRINTF (outfile, "------------------------------\n");
        for (i = 0; i < opt->categs; i++)
            FPRINTF (outfile, "%9ld%16.3f\n", i + 1, opt->rate[i]);
    }
    if ((opt->rcategs > 1) || (opt->categs > 1))
        FPRINTF (outfile, "\n");
}


void
printweights (FILE * outfile, world_fmt * world, option_fmt * options,
              short inc, long chars, short *weight, char *letters)
{
    /* print out the weights of sites */
    long i, j;
    FPRINTF (outfile, "\n    %s are weighted as follows:\n", letters);
    for (i = 0; i < chars; i++)
    {
        if (i % 60 == 0)
        {
            FPRINTF (outfile, "\n");
            for (j = 1; j <= options->nmlength + 3; j++)
               FPRINTF (outfile, " ");
        }
        FPRINTF (outfile, "%hd", weight[i + inc]);
        if ((i + 1) % 10 == 0 && (i + 1) % 60 != 0)
            FPRINTF (outfile, " ");
    }
    FPRINTF (outfile, "\n\n");
}    /* printweights */



void
print_seqfreqs (FILE * outfile, world_fmt * world, option_fmt * options)
{
  if(outfile==NULL)
    return;
  
  if (world->locus == 0 || outfile == stdout)
    {
        if (options->freqsfrom)
            FPRINTF (outfile, "\nEmpirical ");
        FPRINTF (outfile, "Base Frequencies\n");
        FPRINTF (outfile,
                 "------------------------------------------------------------\n");
        FPRINTF (outfile,
                 "Locus     Nucleotide                        Transition/\n");
        FPRINTF (outfile,
                 "          ------------------------------  Transversion ratio\n");
        FPRINTF (outfile, "          A       C       G       T(U)\n");
        FPRINTF (outfile,
                 "------------------------------------------------------------\n");
    }
    FPRINTF (outfile, "%4li      %6.4f  %6.4f  %6.4f  %6.4f    %10.5f\n",
             world->locus + 1, options->freqa, options->freqc, options->freqg,
             options->freqt, world->data->seq[0]->ttratio);
    if (outfile == stdout)
        FPRINTF (outfile, "\n");

}

MYREAL
treelike_seq (world_fmt * world, long locus)
{
    const seqmodel_fmt *seq = world->data->seq[0];
    const MYREAL freqa = seq->basefrequencies[NUC_A];
    const MYREAL freqc = seq->basefrequencies[NUC_C];
    const MYREAL freqg = seq->basefrequencies[NUC_G];
    const MYREAL freqt = seq->basefrequencies[NUC_T];

    contribarr tterm;
    contribarr like;
    contribarr nulike;
    contribarr clai;
    //  long size = sizeof(MYREAL) * world->options->rcategs;
    MYREAL summ, sum2, sumc, sumterm, lterm;
    long i, j, k, lai;
    MYREAL scale;
    node *p;
    sitelike *x1;
    worldoption_fmt *opt;
    opt = world->options;
    p = crawlback (world->root->next);
    summ = 0.0;

    if (opt->rcategs == 1)
    {
        for (i = 0; i < seq->endsite; i++)
        {
            x1 = &(p->x.s[i][0]);
            scale = p->scale[i];
            tterm[0] =
                freqa * (*x1)[0] + freqc * (*x1)[1] +
                freqg * (*x1)[2] + freqt * (*x1)[3];
            summ += seq->aliasweight[i] * (LOG (tterm[0]) + scale);
        }
    }
    else
    {
        for (i = 0; i < seq->endsite; i++)
        {
            scale = p->scale[i];
            //k = seq->category[seq->alias[i] - 1] - 1;
            for (j = 0; j < opt->rcategs; j++)
            {
                x1 = &(p->x.s[i][j]);
                tterm[j] =
                    freqa * (*x1)[0] + freqc * (*x1)[1] +
                    freqg * (*x1)[2] + freqt * (*x1)[3];
            }
            sumterm = 0.0;
            for (j = 0; j < opt->rcategs; j++)
                sumterm += opt->probcat[j] * tterm[j];
            lterm = LOG (sumterm) + scale;
            for (j = 0; j < opt->rcategs; j++)
                clai[j] = tterm[j] / sumterm;
            swap (clai, world->contribution[i]);
            //      memcpy (world->contribution[i], clai, size);
            summ += seq->aliasweight[i] * lterm;
            if(MYISNAN((float)summ))
	      {
		// error("summ is not a number, should not happen\n");
		summ = -HUGE;
	      }

        }
        for (j = 0; j < opt->rcategs; j++)
            like[j] = 1.0;
        for (i = 0; i < seq->sites[locus]; i++)
        {
            sumc = 0.0;
            for (k = 0; k < opt->rcategs; k++)
                sumc += opt->probcat[k] * like[k];
            sumc *= opt->lambda;
            if ((seq->ally[i] > 0) && (seq->location[seq->ally[i] - 1] > 0))
            {
                lai = seq->location[seq->ally[i] - 1];
                swap (world->contribution[lai - 1], clai);
                //memcpy (clai, world->contribution[lai - 1], size);
                for (j = 0; j < opt->rcategs; j++)
                    nulike[j] = ((1.0 - opt->lambda) * like[j] + sumc) * clai[j];
            }
            else
            {
                for (j = 0; j < opt->rcategs; j++)
                    nulike[j] = ((1.0 - opt->lambda) * like[j] + sumc);
            }
            swap (nulike, like);
            //memcpy (like, nulike, size);
        }
        sum2 = 0.0;
        for (i = 0; i < opt->rcategs; i++)
            sum2 += opt->probcat[i] * like[i];
        summ += LOG (sum2);
    }
    return summ;
}    /* treelikelihood */

void
snp_invariants_original (contribarr invariants, long endsite, long rcategs,
                seqmodel_fmt * seq, phenotype x1)
{
    long i, j;
    MYREAL *val;
    register MYREAL freqa, freqc, freqg, freqt;
    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];

        memset (invariants, 0, sizeof (contribarr));
    for (j = 0; j < rcategs; j++)
      {
	for (i = endsite; i < endsite + 4; i++)
	  {
	    val = x1[i][j];
            invariants[j] += 
	      (freqa * val[0] + freqc * val[1] +
	       freqg * val[2] + freqt * val[3]);
	  }
      }
    for (j = 0; j < rcategs; j++)
      invariants[j] = 1. - invariants[j];
}

void
snp_invariants (contribarr invariants, world_fmt *world, long locus, phenotype x1)
{
    worldoption_fmt *opt;
    seqmodel_fmt *seq =  world->data->seq[0];
    //MYREAL summ;

    sitelike *val;
    long i, j;
    register MYREAL freqa, freqc, freqg, freqt;
    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];

    opt = world->options;
    //summ = 0.0;
    /* snp invariants */
    for (j = 0; j < opt->rcategs; j++)
      {
	invariants[j] = 0.0;
	for (i = seq->endsite - seq->addon ; i < seq->endsite; i++)
	  {
	    val = &(x1[i][j]);
            invariants[j] +=
	      (freqa * (*val)[0] + freqc * (*val)[1] +
	       freqg * (*val)[2] + freqt * (*val)[3]);
	  }
	invariants[j] = 1 - invariants[j];
      }
    //    printf("invariant0=%f\n",invariants[0]);
}    /* snp_invariants*/

MYREAL
treelike_snp (world_fmt * world, long locus)
{
    worldoption_fmt *opt;
    seqmodel_fmt *seq =  world->data->seq[0];
    MYREAL scale;
    contribarr tterm;
    contribarr invariants;
    contribarr like;
    contribarr nulike;
    MYREAL summ, sum2, sumc, sumterm, lterm;
    long i, j, k, lai;

    node *p;
    sitelike *x1;
    register MYREAL freqa, freqc, freqg, freqt;
    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];

    opt = world->options;
    p = crawlback (world->root->next);
    summ = 0.0;
    /* snp invariants */
    snp_invariants (invariants,world, locus, p->x.s);
    for (j = 0; j < opt->rcategs; j++)
      {
	if(invariants[j] <= 0.0)
	  {
#ifdef DEBUG
	    printf("invariants are funny %f\n", invariants[j]);
#endif
	    return -MYREAL_MAX;
	  }
      }
    for (i = 0; i < seq->endsite - seq->addon; i++)
    {
        scale = p->scale[i];
        for (j = 0; j < opt->rcategs; j++)
        {
            x1 = &(p->x.s[i][j]);
            tterm[j] =
                (freqa * (*x1)[0] + freqc * (*x1)[1] +
                 freqg * (*x1)[2] + freqt * (*x1)[3]) /invariants[j];
        }
        sumterm = 0.0;
        for (j = 0; j < opt->rcategs; j++)
            sumterm += opt->probcat[j] * tterm[j];
        lterm = LOG (sumterm) + scale;
        for (j = 0; j < opt->rcategs; j++)
	  world->contribution[i][j] = tterm[j] / sumterm;

        summ += seq->aliasweight[i] * lterm;
    }    /* over endsite - 4[snp-invariants] */
    for (j = 0; j < opt->rcategs; j++)
        like[j] = 1.0;
    for (i = 0; i < seq->sites[locus]; i++)
    {
        sumc = 0.0;
        for (k = 0; k < opt->rcategs; k++)
            sumc += opt->probcat[k] * like[k];
        sumc *= opt->lambda;
        if ((seq->ally[i] > 0) && (seq->location[seq->ally[i] - 1] > 0))
        {
            lai = seq->location[seq->ally[i] - 1];
            for (j = 0; j < opt->rcategs; j++)
                nulike[j] =
                    ((1.0 - opt->lambda) * like[j] +
                     sumc) * world->contribution[lai - 1][j];
        }
        else
        {
            for (j = 0; j < opt->rcategs; j++)
                nulike[j] = ((1.0 - opt->lambda) * like[j] + sumc);
        }
        swap (nulike, like);
        //memcpy (like, nulike, size);
    }
    sum2 = 0.0;
    for (i = 0; i < opt->rcategs; i++)
        sum2 += opt->probcat[i] * like[i];
    summ += LOG (sum2);
    return summ;
}    /* treelike_snp */


MYREAL
treelike_snp_unlinked (world_fmt * world, long locus)
{
    //worldoption_fmt *opt;
    seqmodel_fmt *seq = world->data->seq[0];
    MYREAL scale;
    contribarr tterm, invariants;
    MYREAL summ, datasum = 0, lterm, result = 0;
    long i, ii;

    node *p;
    sitelike *x1;

    register MYREAL freqa, freqc, freqg, freqt;

    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];

    //opt = world->options;
    //seq = world->data->seq[0];
    p = crawlback (world->root->next);
    summ = 0.0;
    /* snp invariants */
    snp_invariants (invariants, world, locus, p->x.s);
    /* no rate categories used */
    for (i = 0; i < seq->endsite; i++)
    {
        ii = i / 5;
        scale = p->scale[i];
        x1 = &(p->x.s[i][0]);
        tterm[0] =
            (freqa * (*x1)[0] + freqc * (*x1)[1] +
             freqg * (*x1)[2] + freqt * (*x1)[3]);
        if (i % 5 == 0)
        {
            lterm = LOG (tterm[0]) + scale;
            summ = 0;
            datasum = seq->aliasweight[ii] * lterm;
        }
        else
            summ += pow (tterm[0], (MYREAL) seq->aliasweight[ii]);
        if (((i + 1) % 5) == 0 && i != 0)
	  {
	    if(summ > 0.)
	      {
		result +=
		  datasum + LOG ((1 - EXP (LOG (summ) - datasum)) / invariants[0]);
	      }
	  }
    }
    return result;
}    /* treelike_snp_unlinked */

void
check_basefreq (option_fmt * options)
{

    if (options->freqa == 0. || options->freqc == 0. || options->freqt == 0.
            || options->freqg == 0.)
    {
        options->freqa = 0.25;
        options->freqc = 0.25;
        options->freqg = 0.25;
        options->freqt = 0.25;
    }
}


void
copy_seq (world_fmt * original, world_fmt * kopie)
{
    long sites;
    seqmodel_fmt *kseq;
    seqmodel_fmt *oseq;
    kseq = kopie->data->seq[0];
    oseq = original->data->seq[0];
    sites = oseq->sites[original->locus];
    memcpy(kseq->basefrequencies,oseq->basefrequencies,sizeof(MYREAL)*BASEFREQLENGTH);
    kseq->aa = oseq->aa;
    kseq->bb = oseq->bb;
    kseq->endsite = oseq->endsite;
    kseq->xi = oseq->xi;
    kseq->xv = oseq->xv;
    kseq->ttratio = oseq->ttratio;
    kseq->fracchange = oseq->fracchange;
    memcpy (kseq->sites, oseq->sites, sizeof (long) * original->loci);
    memcpy (kseq->alias, oseq->alias, sizeof (long) * sites);
    memcpy (kseq->ally, oseq->ally, sizeof (long) * sites);
    memcpy (kseq->category, oseq->category, sizeof (long) * sites);
    memcpy (kseq->weight, oseq->weight, sizeof (short) * sites);
    kseq->weightsum = oseq->weightsum;
    memcpy (kseq->aliasweight, oseq->aliasweight, sizeof (long) * sites);
    memcpy (kseq->location, oseq->location, sizeof (long) * sites);
    kseq->addon = oseq->addon;
}

void find_rates_fromdata(data_fmt * data, option_fmt * options)
{
  long locus;
  long pop;
  long ind;
  long sites;
  long i;
  long n;
  long s;
  long **v;
  long maxsites = 1;
  long *numind;
  char *reference;
  char *indseq;
  long loci = data->loci;
  long numpop = data->numpop;
  long segreg;
  MYREAL mean=0.;
  MYREAL delta = 0.;

  numind = (long *) mycalloc (data->loci, sizeof (long));

  for (locus=0; locus < loci; locus++)
    {
      maxsites = (data->seq[0]->sites[locus] > maxsites ? data->seq[0]->sites[locus] : maxsites);
      for(pop=0;pop < data->numpop; pop++)
	{
	  numind[locus] += data->numind[pop][locus];
	}
    }

  v = (long **) mycalloc (loci, sizeof (long *));
  v[0] = (long *) mycalloc (loci * maxsites, sizeof (long));
  for (i = 1; i < loci; i++)
    {
      v[i] = v[0] + maxsites * i;
    }
  
  for (locus=0; locus< loci; locus++)
    {
      reference = data->yy[0][0][locus][0];
      for(pop=0;pop < numpop; pop++)
	{
	  for(ind=0; ind < data->numind[pop][locus];ind++)
	    {
	      indseq=data->yy[pop][ind][locus][0];
	      if(strcmp(indseq,reference)!=0)
		{
		  for(sites=0; sites < data->seq[0]->sites[locus]; sites++)
		    {
		      if(v[locus][sites] == 0)
			{
			  s = (long)(indseq[sites] != reference[sites]);
			  v[locus][sites] += s;
			}
		    }
		}
	    }
	}
    }
  //n = 0;
  //mean = 0.;
  if(options->mu_rates == NULL)
    {
      options->mu_rates = (MYREAL *) mycalloc( loci, sizeof(MYREAL));
      //      printf("%i> opt murate size %li\n",myID,loci * sizeof (MYREAL));
    }
  n = 0;
  mean = 0. ;
  options->muloci = loci;
  fprintf(stdout,"Relative mutation rate among loci estimated from the data\n");
  fprintf(stdout,"---------------------------------------------------------\n");
  for (locus=0; locus < loci; locus++)
    {
      segreg = 0;
      for(sites=0; sites < data->seq[0]->sites[locus]; sites++)
	{
	  segreg += (long) (v[locus][sites] > 0); 
	}
      // originally U use just the plain watterson estimates but
      // this leads to strange results with highly unequal sequences
      // now the watterson's theta is adjusted per site
      // and that is used to generate the relative rate from the data
      // PB April 3 2011
      options->mu_rates[locus] = (watterson(segreg,numind[locus]) + 0.0000001)/data->seq[0]->sites[locus];
      n = n + 1;
      delta = options->mu_rates[locus] - mean;
      mean += delta/n;
    }
  if(mean > 0.0)
    {
      for (locus=0; locus < loci; locus++)
	{
	  options->mu_rates[locus] /= mean;
	  if(options->progress)
	    {
	      fprintf(stdout,"Locus %2li: rate=%f (Watterson's overall Theta = %f)\n",locus+1, options->mu_rates[locus],options->mu_rates[locus] * mean);
	    }
	}
    }
  else
    {
      options->murates = FALSE;
      options->murates_fromdata = FALSE;
      warning("Relative mutation rates estimation from data was requested but attempt failed --> reset to all-equal rates");
    }
  fprintf(stdout,"\n\n\n");
  myfree(v[0]);
  myfree(v);
  myfree(numind);
}


void free_seq(seqmodel_fmt **seq, long seqnum)
{
  long i;
  for(i=0;i<seqnum;i++)
    {
      if(seq[i]->links!=NULL)
	{
	  myfree(seq[i]->links);
	}
      myfree(seq[i]->sites);
    }
  myfree(seq[0]);
  myfree(seq);
}
