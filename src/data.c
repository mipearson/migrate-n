/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 D A T A   R O U T I N E S 
 
 creates data structures,
 read data (Electrophoretic loci, sequences, microsats),
 feeds data into tree (?),
 prints data,
 destroys data.
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2007 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
$Id: data.c 1814 2011-03-20 17:21:07Z beerli $

-----------------------------------
add this to the microsatellite data, this will allow to use
a known repeat lengths to read number of sites instead of 
number of repeats, the datafile will need an additional line
perhaps use something like #@ 2 2 2 3  .... on the second line
so that other non microsatellite data models 
can ignore it
y = Map[If[(a = Mod[#, rlen]) != 0, 
      If[Random[] < 0.5, # - (rlen + a), # + rlen - a], #] &, x]/
   rlen  - Min[y] + 1 
-------------------------------------------------------*/
/*! \file data.c

Data manipulation routines

*/


#include <string.h>

#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#include "migrate_mpi.h"
#include "data.h"
#include "sequence.h"
#include "random.h"
#ifdef PRETTY
#include "pretty.h"
#endif
extern long number_genomes (int type);
extern void jumble (long *s, long n);

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
/* prototypes ----------------------------------------- */
//void create_data (data_fmt ** data);
//void get_data (FILE * infile, data_fmt * data, option_fmt * options);
//void print_data (world_fmt * world, option_fmt * options, data_fmt * data);
//long find_missing(data_fmt *data, long pop, long locus);
//void print_data_summary (FILE * file, world_fmt * world, option_fmt * options,
//                         data_fmt * data);
//short findAllele (data_fmt * data, char s[], long locus);
//void free_datapart (data_fmt * data, option_fmt * options, long locus);
/*private functions */
//void init_data_structure1 (data_fmt ** data, option_fmt * options);
void read_header (FILE * infile, data_fmt * data, option_fmt * options);
void read_sites (data_fmt * data, world_fmt *world, option_fmt *options);
//void init_data_structure2 (data_fmt ** data, option_fmt * options, long pop);
//void init_data_structure3 (data_fmt * data);
void read_popheader (FILE * infile, data_fmt * data, option_fmt * options, long pop, long genomes);
void read_indname (FILE * file, data_fmt * data, long pop, long ind,
                   long locus, long nmlength);
void read_popdata (FILE * file, data_fmt * data, long pop,
                   option_fmt * options);
void read_microalleles (FILE * infile, data_fmt * data, option_fmt *options, long pop, long ind);
void read_alleles (FILE * infile, data_fmt * data, long pop, long ind);
long read_ind_seq (FILE * infile, data_fmt * data, option_fmt * options,
                   long locus, long pop, long ind, long baseread);
long read_ind_seq_oneliner (FILE * infile, data_fmt * data, option_fmt * options,
		       long pop, long ind, long baseread);

void read_hapmap (FILE * infile, data_fmt * data, option_fmt * options,
	     long locus, long pop);
void read_distance_fromfile (FILE * dfile, long tips, long nmlength,
                             MYREAL **m);
void finish_read_seq (FILE * infile, data_fmt * data, option_fmt * options,
                      long pop, long baseread);
void print_alleledata (world_fmt * world, data_fmt * data,
                       option_fmt * options);
void print_microdata (world_fmt * world, data_fmt * data,
                       option_fmt * options);
void print_seqdata (world_fmt * world, option_fmt * options, data_fmt * data);

void print_header (FILE * outfile, long pop, world_fmt * world,
                   option_fmt * options, data_fmt * data);
void create_alleles (data_fmt * data, option_fmt *options);
void addAllele (data_fmt * data, char s[], long locus, long *z);
void set_numind (data_fmt * data);
void print_seq_pop (long locus, long pop, world_fmt * world,
                    option_fmt * options, data_fmt * data);
void print_seq_ind (long locus, long pop, long ind, world_fmt * world,
                    option_fmt * options, data_fmt * data);
void print_locus_head (long locus, world_fmt * world, option_fmt * options,
                       data_fmt * data);
void find_delimiter (char *title, char *dlm);
void read_geofile (data_fmt * data, option_fmt * options, long numpop);
void read_datefile (data_fmt * data, option_fmt * options, long numpop);
void read_uep_fromfile (FILE * uepfile, long tips, long nmlength, int **uep,
                        long *uepsites, long datatype);
void read_uepfile (data_fmt * data, option_fmt * options, long numpop);


/*=====================================================*/
/// creates the data structure
void
create_data (data_fmt ** data)
{
    (*data) = (data_fmt *) mycalloc (1, sizeof (data_fmt));
}

/*
void
init_data (data_fmt * data)
{
 
}
*/

///
/// free the data module 
void
destroy_data (data_fmt * data)
{
  long ind;
  long indalloc;
  long locus;
  long pop;
  long loci = data->loci;
  long numpop = data->numpop;

  // free data from init_data_structure3
  for (locus = 0; locus < loci; locus++)
    {
      myfree(data->allele[locus]);
    }
  myfree(data->allele);
  myfree(data->subloci);
  myfree(data->maxalleles);
  myfree(data->skiploci);

  // free data from init_data_structure2
  for(pop=0; pop < numpop ; pop++)
    {
      indalloc = -1;
      for(locus=0; locus < loci; locus++)
	{
	  if(indalloc < data->numind[pop][locus])
	    indalloc = data->numind[pop][locus];
	}
      for (ind = 0; ind < indalloc; ind++)
	{
	  for(locus=0; locus < loci; locus++)
	    {
	      myfree(data->indnames[pop][ind][locus]);
	      myfree(data->yy[pop][ind][locus]);
	    }
	  myfree(data->indnames[pop][ind]);
	  myfree(data->yy[pop][ind]);
	}
      myfree(data->indnames[pop]);
      myfree(data->yy[pop]);
    }
  myfree(data->indnames);
  myfree(data->yy);

  // data->yy were already freed in free_datapart()

  // free data from init_structure_1
  myfree(data->popnames[0]);
  myfree(data->numind[0]);
  myfree(data->numalleles[0]);
  myfree(data->popnames);
  myfree(data->numind);
  myfree(data->numalleles);
  
  if(data->position!=NULL)
    myfree(data->position);
  myfree(data->geo);
  myfree(data->lgeo);
  if(data->ogeo != NULL)
    {
      myfree(data->ogeo[0]);
      myfree(data->ogeo);
    }
  // free sampledates
  for(locus=0;locus<data->loci;locus++)
    {
      for(pop=0;pop < data->numpop; pop++)
	{
	  myfree(data->sampledates[pop][locus]);
	}
    }
  myfree(data->sampledates[0]);
  myfree(data->sampledates);
  myfree(data);
}
 
///
/// shuffles (jumbles) the individuals so that we can take the first y to subsample the dataset
void shuffle(long **shuffled, long n, boolean shuffle_ON)
{
  long i;
  for(i=0;i<n;i++)
    {
      (*shuffled)[i] = i;
    }
  if(shuffle_ON)
    {
        jumble (*shuffled, n);
    }
}

void shuffle_data(data_fmt *data, option_fmt *options)
{
  long pop;
  long locus;
  data->shuffled = (long ***) mycalloc(data->numpop,sizeof(long **));
  for(pop=0; pop < data->numpop; pop++)
    {
      data->shuffled[pop] = (long **) mycalloc(data->loci, sizeof(long *));
      for(locus=0; locus < data->loci; locus++)
	{
	  data->shuffled[pop][locus] = (long *) mycalloc(data->numind[pop][locus], sizeof(long));
	  if(options->randomsubset>0 && options->randomsubset <  data->numind[pop][locus])
	    {	
	      shuffle(&data->shuffled[pop][locus],data->numind[pop][locus],TRUE);
	    }
	  else
	    {
	      shuffle(&data->shuffled[pop][locus],data->numind[pop][locus],FALSE);
	    }
	}
    }
}


long max_shuffled_individuals(option_fmt *options, data_fmt *data, long pop, long locus)
{
  long top;
  if(options->randomsubset>0 && options->randomsubset <  data->numind[pop][locus])
    {	
      top = options->randomsubset;
    }
  else
    {
      top = data->numind[pop][locus];
    }
  return top;
}


//---------------------------------------------------------------
// ..number_pop.number_loci.<delimiter>.<title [no numbers to start]
// locu_specification {sNumber_of_sites}{n{number of sites}{m1}{a1}
// [parentheses () mark link loci
// example:
// (n3 s120 m1)(m1 m1 m1 a1 s100)(n1)(n2)(m1)
// contains 5 linked loci, the first two are compounds
// the system is to have for each individual haplotype a set of 
// loci so that it may look like this
// ind11......ACA ATTAGA 13 15 5 3 2 A ATACG A AT 13
// ind12......ACA ATTCGA 12 15 6 3 2 A ATACG T CT 5
void
get_new_data (FILE * infile, data_fmt * data, option_fmt * options, world_fmt *world)
{
  //long locus;
  long pop;
  long genomes=1;
  data->hasghost = FALSE;
  // read how many populations and loci and delimiter and title
  // if we have the new haplotype input then there should be no delimiter
  // needed
  read_header (infile, data, options);
  init_data_structure1 (&data, options);
  //
  // here I need to add the msat_lengths support
  //
  // this will set up all the locus structure, but may need some additional
  // init_data-structure1a()
  read_sites(data, world, options);
  // snps etc need something like
  //        data->seq[0]->addon = 4;
  // hapmap data is peculiar check that
  // if not sequence the fracchange [for F84] needs to be 1.0
  //        data->seq[0]->fracchange = 1.0;

  if (options->progress)
    fprintf (stdout, "\n\n");
  if (options->writelog)
    fprintf (options->logfile, "\n\n");
  
  for (pop = 0; pop < data->numpop; pop++)
    {
      read_popheader (infile, data, options, pop, genomes);
      if (options->progress)
	fprintf (stdout, "Reading (%li) %s ...\n", pop, data->popnames[pop]);
      if (options->writelog)
	fprintf (options->logfile, "Reading (%li) %s ...\n", pop, data->popnames[pop]);
      init_data_structure2 (&data, options, pop);
      read_popdata (infile, data, pop, options);
    }
  read_geofile (data, options, data->numpop);
#ifdef UEP

    read_uepfile (data, options, data->numpop);
#endif
    read_datefile(data, options, data->numpop);
    if (options->progress)
        fprintf (stdout, "\n\n");
    init_data_structure3 (data);

    switch (options->datatype)
    {
    case 'a':
        create_alleles (data, options);
        break;
    case 'b':
        create_alleles (data, options);
        for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = XBROWN_SIZE;
        break;
    case 'm':
        create_alleles (data, options);
        for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = options->micro_stepnum;
        break;
    default: /*DNA types*/
      for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = 4;
    }

    shuffle_data(data, options);

    if(options->murates_fromdata)
      {
	if (strchr (SEQUENCETYPES, options->datatype))
	  {
	    find_rates_fromdata(data, options);
	  }
      }
}

void
get_data (FILE * infile, data_fmt * data, option_fmt * options, world_fmt *world)
{
  long locus;
    long pop;
    long genomes=1;
    data->hasghost = FALSE;
    read_header (infile, data, options);
    genomes =  number_genomes (options->datatype);
    init_data_structure1 (&data, options);
    data->datatype[0] = options->datatype;
    switch (options->datatype)
    {
    case 's':
    case 'f':// read standard sequence data
      read_sites (data,world, options);
        break;
    case 'n': // read snp data
      read_sites (data,world, options);
        data->seq[0]->addon = 4;
        break;
    case 'h': //read data from hapmap allele frequency files
      // there is no sites line in these files!
      for(locus=0;locus<data->loci;locus++)
	data->seq[0]->sites[locus] = 1;
      data->seq[0]->addon = 4;
        break;
    case 'u': // read unlinked snp
      read_sites (data,world, options);
        data->seq[0]->addon = 4;
        break;
    default:
      read_sites(data,world, options); // checks whether there is a line for microsat repeat numbers
      data->seq[0]->fracchange = 1.0;
      break;
    }
    if (options->progress)
        fprintf (stdout, "\n\n");
    if (options->writelog)
        fprintf (options->logfile, "\n\n");
    for (pop = 0; pop < data->numpop; pop++)
    {
      read_popheader (infile, data, options, pop, genomes);
        if (options->progress)
            fprintf (stdout, "Reading (%li) %s ...\n", pop, data->popnames[pop]);
        if (options->writelog)
            fprintf (options->logfile, "Reading (%li) %s ...\n", pop, data->popnames[pop]);
        init_data_structure2 (&data, options, pop);
        read_popdata (infile, data, pop, options);
    }
    read_geofile (data, options, data->numpop);
#ifdef UEP

    read_uepfile (data, options, data->numpop);
#endif
    read_datefile(data, options, data->numpop);
    if (options->progress)
        fprintf (stdout, "\n\n");
    init_data_structure3 (data);

    switch (options->datatype)
    {
    case 'a':
        create_alleles (data, options);
        break;
    case 'b':
        create_alleles (data, options);
        for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = XBROWN_SIZE;
        break;
    case 'm':
        create_alleles (data, options);
        for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = options->micro_stepnum;
        break;
    default: /*DNA types*/
      for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = 4;
    }

    shuffle_data(data, options);

    if(options->murates_fromdata)
      {
	if (strchr (SEQUENCETYPES, options->datatype))
	  {
	    find_rates_fromdata(data, options);
	  }
      }
}

/* private functions ========================================== */

void
init_data_structure1 (data_fmt ** data, option_fmt * options)
{
    long pop;
    long locus;
    long numpop = (*data)->numpop;
    long loci   = (*data)->loci;
    
    (*data)->ogeo = NULL;
    (*data)->geo = NULL;
    if ((*data)->yy == NULL)
    {
        (*data)->yy = (char *****) mymalloc (sizeof (char ****) * numpop);
        (*data)->seq = (seqmodel_fmt **) mycalloc (1, sizeof (seqmodel_fmt *));
        (*data)->seq[0] = (seqmodel_fmt *) mycalloc (1, sizeof (seqmodel_fmt));
        (*data)->popnames =(char **) mymalloc (sizeof (char *) * numpop);
        (*data)->popnames[0] =(char *) mycalloc (numpop * STRSIZE,sizeof(char));
        (*data)->indnames = (char ****) mymalloc (sizeof (char ***) * numpop);
        (*data)->numind = (long **) mymalloc (sizeof (long *) * numpop);
        (*data)->numind[0] = (long *) mymalloc (sizeof (long) * numpop * loci);
        (*data)->numalleles = (long **) mymalloc (sizeof (long *) * numpop);
        (*data)->numalleles[0] = (long *) mymalloc (sizeof (long) * numpop * loci);
	(*data)->locitypes = (char *) mycalloc(loci, sizeof(char));
	(*data)->locusname = (char **) mycalloc(loci, sizeof(char*));
        for (pop = 1; pop < numpop; pop++)
        {
	  (*data)->popnames[pop] = (*data)->popnames[0] + pop * STRSIZE;
	  (*data)->numind[pop] = (*data)->numind[0] + pop * loci;
	  (*data)->numalleles[pop] =  (*data)->numalleles[0] + pop * loci;
        }
        (*data)->seq[0]->sites = (long *) mycalloc (loci, sizeof (long));
        (*data)->seq[0]->links = (boolean *) mycalloc (loci, sizeof (boolean));
        (*data)->position = (long *) mycalloc (loci, sizeof (long));
	(*data)->datatype = (long *) mycalloc(loci, sizeof(long));
	(*data)->numdatatypealloc = loci;
	(*data)->numsublocialloc = loci;
	(*data)->subloci = (long *) mycalloc(loci,sizeof(long));
	(*data)->numrepeatnumbers = (long *) calloc(loci, sizeof(long));
	(*data)->repeatnumbers = (long **) calloc(loci, sizeof(long *));
	(*data)->repeatlength = (long *) calloc(loci, sizeof(long));
	for(locus=0;locus<loci;locus++)
	  {
	    (*data)->numrepeatnumbers[locus] = 1;
	    (*data)->repeatnumbers[locus] = (long *) calloc(1, sizeof(long));
	  }
    }
    else
    {
        error ("Problem with initialization of data matrix yy\n");
    }
}


void
init_data_structure2 (data_fmt ** data, option_fmt * options, long pop)
{
    long ind, locus;
    long indalloc = -1;
    for(locus=0;locus<(*data)->loci;locus++)
    {
        if(indalloc < (*data)->numind[pop][locus])
            indalloc = (*data)->numind[pop][locus];
    }
    if (indalloc == 0)
        indalloc = 2;

    (*data)->yy[pop] = (char ****) mymalloc (sizeof (char ***) * indalloc);
    (*data)->indnames[pop] = (char ***) mycalloc (1, sizeof (char **) * indalloc);

    for (ind = 0; ind < indalloc; ind++)
    {
        (*data)->indnames[pop][ind] =
            (char **) mymalloc (sizeof (char *) * (*data)->loci);
        (*data)->yy[pop][ind] =
            (char ***) mymalloc (sizeof (char **) * (*data)->loci);
        for (locus = 0; locus < (*data)->loci; locus++)
        {
	  (*data)->indnames[pop][ind][locus] =
            (char *) mycalloc (1, sizeof (char) * (1 + options->nmlength));
            if (!strchr (SEQUENCETYPES, options->datatype))
            {
                (*data)->yy[pop][ind][locus] =
                    (char **) mycalloc (1, sizeof (char *) * 2);
                (*data)->yy[pop][ind][locus][0] =
                    (char *) mycalloc (1, sizeof (char) * (options->allelenmlength+1));
                (*data)->yy[pop][ind][locus][1] =
                    (char *) mycalloc (1, sizeof (char) * (options->allelenmlength+1));
            }
            else
            {
                (*data)->yy[pop][ind][locus] =
                    (char **) mycalloc (1, sizeof (char *));
                (*data)->yy[pop][ind][locus][0] =
                    (char *) mycalloc (1,
                                     sizeof (char) * ((*data)->seq[0]->sites[locus]+1));
            }
        }
    }
}


short
findAllele (data_fmt * data, char s[], long locus)
{
    short found = 0;
    while ((strcmp (s, data->allele[locus][found])
            && data->allele[locus][found][0] != '\0'))
        found++;
    return found;
}

void
free_datapart (data_fmt * data, option_fmt * options, long locus)
{
    long ind, pop;
    //  long genomes = number_genomes (options->datatype);
    for (pop = 0; pop < data->numpop; pop++)
    {
        for (ind = 0; ind < data->numind[pop][locus]; ind++)
        {
            if (!strchr (SEQUENCETYPES, options->datatype))
            {
                myfree(data->yy[pop][ind][locus][0]);
                myfree(data->yy[pop][ind][locus][1]);
                myfree(data->yy[pop][ind][locus]);
            }
            else
            {
                myfree(data->yy[pop][ind][locus][0]);
                myfree(data->yy[pop][ind][locus]);
            }
        }
    }
}


void
init_data_structure3 (data_fmt * data)
{
    long locus, pop, maxi;
    data->allele =
        (allele_fmt **) mycalloc (1, sizeof (allele_fmt *) * data->loci);
    data->subloci = (long *) mycalloc(data->loci,sizeof(long));
    for (locus = 0; locus < data->loci; locus++)
    {
      data->subloci[locus] = 1;
        maxi = 0;
        for (pop = 0; pop < data->numpop; pop++)
            maxi += data->numalleles[pop][locus];
        data->allele[locus] =
            (allele_fmt *) mycalloc (1, sizeof (allele_fmt) * maxi);
    }
    data->maxalleles = (long *) mycalloc (1, sizeof (long) * data->loci);
    data->skiploci =
        (boolean *) mycalloc (1, sizeof (boolean) * (data->loci + 1));
}

///
/// read the first line of the data file
/// \param infile datafilename
/// \param data   data structure that holds all the data
/// \param options structure that contain all option information
/// \reval none
void read_header (FILE * infile, data_fmt * data, option_fmt * options)
{
    char input[LINESIZE], *p;
    char title[LINESIZE];
    strcpy(title,"\0");
    input[0] = '#';
    while(input[0]=='#')
      {
	FGETS (input, sizeof (input), infile);
	if ((p = (char *) strpbrk (input, CRLF)) != NULL)
	  *p = '\0';
      }
    switch (lowercase (input[0]))
    {
    case 'a':
        sscanf (input, "%1s%ld%ld%[^\n]", &options->datatype, &(data->numpop),
                &(data->loci), title);
        find_delimiter (title, &data->dlm);
	if(!(title[0] == '\0'))
	  strcpy(options->title,title);
        break;
    case 'b':
    case 'm':
        sscanf (input, "%1s%ld%ld%1s%[^\n]", &options->datatype,
                &(data->numpop), &(data->loci), &data->dlm, title);
	if(!(title[0] == '\0'))
	  strcpy(options->title,title);
        break;
    case 's':
    case 'n':
    case 'h': //hapmap data
    case 'u':
    case 'f':
        sscanf (input, "%1s%ld%ld%[^\n]", &options->datatype, &(data->numpop),
                &(data->loci), title);
	if(!(title[0] == '\0'))
	  strcpy(options->title,title);
        break;
    case 'g':   /* fall through if a menu change forces to analyze data
                               instead of using the already sampled genealogies */
        if (options->datatype == 'g')
            break;
        else
            memmove (input, input + 1, (strlen (input) - 1) * sizeof (char));
    default:
      if(input[0]== '<')
	{
	  usererror ("This data file may contain an XML or HTML tag,\nand cannot be read properly, check the data formatting section in the manual!");
	  exit(-1);
	}
        switch (options->datatype)
        {
        case 'a':
            sscanf (input, "%ld%ld%[^\n]", &(data->numpop), &(data->loci),
                    title);
            find_delimiter (title, &data->dlm);
	    
	if(!(title[0] == '\0'))
	  strcpy(options->title,title);
            break;
        case 'b':
        case 'm':
            sscanf (input, "%ld%ld%1s%[^\n]", &(data->numpop), &(data->loci),
                    &(data->dlm), title);
	    if(!(title[0] == '\0'))
	      strcpy(options->title,title);
            break;
        case 's':
        case 'n':
	case 'h': // hapmap data
        case 'u':
        case 'f':
            sscanf (input, "%ld%ld%[^\n]", &(data->numpop), &(data->loci),
                    title);
	    if(!(title[0] == '\0'))
	      strcpy(options->title,title);
            break;
        default:
            usererror ("Datatype is wrong, please use a valid data type!");
        }
    }
    options->datatype = lowercase (options->datatype);
}

void
find_delimiter (char *title, char *dlm)
{
    char *p = title;
    long z = 0;
    while (*p == ' ')
    {
        p++;
        z++;
    }
    if (isalnum (*p))
        memmove (title, p, sizeof (char) * (strlen (title) - z));
    else
    {
        *dlm = *p;
        p++;
        while (*p == ' ')
        {
            p++;
            z++;
        }
        memmove (title, p, sizeof (char) * (strlen (title) - z));
    }
}

void set_datatype(char input, long locus, char ** locitypes, long *alloc)
{
  if(locus >= *alloc)
    {
      *alloc += SUBLOCICHUNKS;
      *locitypes = (char *) myrealloc(*locitypes, sizeof(char) * (*alloc));
    }
  (*locitypes)[locus] = input;
}

//==================================================================
/// read the number of sites for each locus in the dataset,
/// does not assume a fixed line length, but assumes that at the end of the line is either
/// a \n or \r or \r\l (similar to the sequence reader) to accommodate windows, mac and
/// unix line ends.
/// The generalized sites reader can also read link status using a parenthesis notation, 
/// and allows to give short names to the
/// loci, example: (ldh=s234 agdh=n3) gpi=s100  str1=b1 (str2=m1 a1)
/// the labels a, s, n etc are the same used for the datatype and will mark each locus
/// for the datatype=h it is assumed that all loci are the same and of type n
/// plain numbers are considered to be sequence loci with no name
void
read_sites (data_fmt * data, world_fmt *world, option_fmt *options)
{
  char * val=NULL;
  boolean linked=TRUE;
  long allocbufsize=LINESIZE;
    long locus;
    char *input;
    char * word;
    long oldi;
    long i;
    long len;
    long z=0; //z is adding subloci and is reset with a closing ')'
    long zz=0; //zz is used for mutationmodels and is not reset at closing ')'
    input = (char *) mycalloc(allocbufsize ,sizeof(char));
    word = (char *) mycalloc(allocbufsize ,sizeof(char));
    // check whether the microsats are not repeats but fragmentlength  
    data->has_repeats = TRUE;
    for (locus = 0; locus < data->loci; locus++)
    {
      read_word(data->infile, input);
      switch(input[0])
	{
	case  '#':
	  if(input[1]=='@')
	    {
	      // line is microsatellite instruction
	      if(input[2]=='M')
		{
		  FGETS2(&input,&allocbufsize,data->infile); // read the line 
		  locus--;
		  oldi=0;
		  // long i;
		  len=0;
		  data->has_repeats = FALSE;
		  for(i = 0 ; i < data->loci; i++)
		    {
		      len += read_word_delim(input+len,word," ;,\t");
		      if(word[0]!='\0')
			{
			  data->repeatlength[i] = atol(word);
			  oldi=i;
			}
		      else
			{
			  data->repeatlength[i] = data->repeatlength[oldi];
			}
		    }
		}
	      else
		{
		  FGETS2(&input,&allocbufsize,data->infile); // read the line 
		  locus--;
		}
	      myfree(input);
	      myfree(word);
	      return;
	    }
	  else
	    {
	      // line is a comment
	      FGETS2(&input,&allocbufsize,data->infile); // read the line 
	      locus--;
	      continue;
	    }
	  break;
	  
	  case '(':
	    // unlinked loci are all on one line!
	    // this is a change from the old format, dna is now treated the same way as other loci
	    // I am not sure whether this really works, but would make things simpler,
	    // still unclear how to deal with haplotype versus diplotype information
	    data->oneliner=TRUE;
	    linked=TRUE;
	    if(z>=data->numsublocialloc)
	      {
		data->numsublocialloc += SUBLOCICHUNKS;
		data->subloci = (long *) myrealloc(data->subloci, data->numsublocialloc * sizeof(long));
	      }
	    data->subloci[z] += 1; //we opened a ( and therefore start a linkage group
	    if(input[1]!='\0' && isalpha(input[1]))
	      {
		read_word_delim(input,word,"=");
		if(input[0]!='\0')
		  {
		    data->locusname[z] = (char *) mycalloc(strlen(word),sizeof(char));
		    strcpy(data->locusname[locus],word+1);
		    // what is this good for?
		    if(!strcmp(word,input))
		      val=strchr((char *) "bmsnuh", input[0]);
		    else
		      val=strchr((char *) "bmsnuh", input[1]);

		    if(val==0)
		      {
                        set_datatype((input[0]=='(' ? input[1] : input[0]),locus, &data->locitypes, &data->numlocitypesalloc);
			zz++;
		      }
		    else
		      {
			// no datatype specified in the model therewe assume there is a unique data type specified in the parmfile or menu
                        set_datatype(data->datatype[0],locus, &data->locitypes, &data->numlocitypesalloc);
		      }
		  }
		else
		  {
		    // no datatype specified in the model therefore we assume there is a unique data type specified in the parmfile or menu
                    set_datatype(data->datatype[0],locus, &data->locitypes, &data->numlocitypesalloc);
		  }
	      }
	    else
	      {
		if(strlen(input)==1 || input[1]==' ')
		  {
		    locus--;
		    continue;
		  }
	      }
	    break;
	  case ')':
	    z=0; // reset subloci for the next locus
	    break;
	  default:
	    if (!strchr (SEQUENCETYPES, options->datatype))
	      {
                unread_word(data->infile, input);
		myfree(input);
		myfree(word);
		return;
	      }
	    if(isalpha(input[0]))
	      {
		read_word_delim(input,word,"=");
		data->locusname[locus] = (char *) mycalloc(strlen(word),sizeof(char));
		strcpy(data->locusname[locus],word);
	      }
	  }
	if(linked)
	  data->seq[0]->links[locus] = TRUE;
	else
	  data->seq[0]->links[locus] = FALSE;
	if(input[0]=='(')
	  {
	    if(strchr("bmsnuh",input[1]))
	      {
		val = input+2;
	      }
	    else
	      {
		val = input+1;
	      }
	  }
	else
	  {
	    if(strchr("bmsnuh",input[0]))
	      {
		val = input+1;
	      }
	    else
	      {
		val = input;
	      }
	  }
        data->seq[0]->sites[locus] = atoi (val);
	val = NULL;
        if (data->seq[0]->sites[locus] == 0)
        {
	  warning ("This does look like sequence data\n");
	  warning ("I just read a number of sites=0\n");
	  warning ("If you use the wrong data type, the program\n");
	  usererror ("will crash anyway, so I stop now\n");
        }
    }
    FGETS2(&input,&allocbufsize,data->infile);
    while(input[0]=='#')
      FGETS2(&input,&allocbufsize,data->infile);
    myfree(input);
    myfree(word);
}

//old version
void
read_old_sites (data_fmt * data)
{
  char * val=NULL;
  boolean linked=TRUE;
    long locus;
    char *input;
    char * word;

    input = (char *) mycalloc(LONGLINESIZE ,sizeof(char));
    word = (char *) mycalloc(LONGLINESIZE ,sizeof(char));
    
    for (locus = 0; locus < data->loci-1; locus++)
    {
        read_word(data->infile, input);
	switch(input[0])
	  {
	  case  '#':
	    FGETS(input,LINESIZE,data->infile);
	    locus--;
	    continue;
	    break;

	  case '(':
	    // unlinked loci are all on one line!
	    // this is a change from the old format, dna is now treated the same way as other loci
	    // I am not sure whether this really works, but would make things simpler,
	    // still unclear how to deal with haplotype versus diplotype information
	    data->oneliner=TRUE;
	    if(input[1]!='\0' && isalpha(input[1]))
	      {
		read_word_delim(input,word,"=");
		if(input[0]!='\0')
		  {
		    data->locusname[locus] = (char *) mycalloc(strlen(word),sizeof(char));
		    strcpy(data->locusname[locus],word+1);
		    if(!strcmp(word,input))
		      val=strchr((char *) "bmsnuh", input[0]);
		    else
		      val=strchr((char *) "bmsnuh", input[1]);


                    if(val==0)
                      {
                        set_datatype((input[0]=='(' ? input[1] : input[0]),locus, &data->locitypes, &data->numlocitypesalloc);
                        //zz++;
                      }
                    else
                      {
                        // no datatype specified in the model therefore we assume there is a                                                  
                        // unique data type specified in the parmfile or menu                                                                 
                        set_datatype(data->datatype[0],locus, &data->locitypes, &data->numlocitypesalloc);
                      }
		  }
		else
		  {
                    // no datatype specified in the model therefore we assume there is a unique data type specified in the parmfile or menu   
                    set_datatype(data->datatype[0],locus, &data->locitypes, &data->numlocitypesalloc);

		  }
	      }
	    else
	      {
		if(strlen(input)==1 || input[1]==' ')
		  {
		    locus--;
		    continue;
		  }
	      }
	    break;
	  case '|':
	    if(linked==FALSE)
	      linked=TRUE;
	    else
	      linked=FALSE;
	    break;
	  default:
	    if(isalpha(input[0]))
	      {
		read_word_delim(input,word,"=");
		data->locusname[locus] = (char *) mycalloc(strlen(word),sizeof(char));
		strcpy(data->locusname[locus],word);
	      }
	  }
	if(linked)
	  data->seq[0]->links[locus] = TRUE;
	else
	  data->seq[0]->links[locus] = FALSE;
	if(input[0]=='(')
	  {
	    if(strchr("bmsnuh",input[1]))
	      {
		val = input+2;
	      }
	    else
	      {
		val = input+1;
	      }
	  }
	else
	  {
	    if(strchr("bmsnuh",input[0]))
	      {
		val = input+1;
	      }
	    else
	      {
		val = input;
	      }
	  }
        data->seq[0]->sites[locus] = atoi (val);
	val = NULL;
        if (data->seq[0]->sites[locus] == 0)
        {
	  warning ("This does look like sequence data\n");
	  warning ("I just read a number of sites=0\n");
	  warning ("If you use the wrong data type, the program\n");
	  usererror ("will crash anyway, so I stop now\n");
        }
    }
    FGETS(input,LINESIZE,data->infile);
    while(input[0]=='#')
      FGETS(input,LINESIZE,data->infile);
    if(linked)
      data->seq[0]->links[locus] = TRUE;
    else
      data->seq[0]->links[locus] = FALSE;
    if(strchr("bmsnuh",input[0]))
      {
	val = input+1;
      }
    else
      {
	val = input;
	while(!strchr("123456789",*val))
	      val++;
      }
    data->seq[0]->sites[locus] = atoi (val);
    
    myfree(input);
    myfree(word);
}


void
read_popheader (FILE * infile, data_fmt * data, option_fmt * options, long pop, long genomes)
{
    boolean havepopname = FALSE;
    long    minlength   = 0;
    long    lo;
    long    locus;
    char   *input;

    input = (char *) mycalloc(LINESIZE,sizeof(char));

    // allows that sequence data can have different numbers of individuals for different loci
    // data syntax changes: #ind1 #ind2 #IND3 .... pop_name
    havepopname=FALSE;
    // with multiple loci we need to separate the nubers of individuals per locus
    // from the population title, and also take care in case the number of 
    // individuals is not specified for all loci.
    if(data->loci>1) 
    {
        read_word(data->infile, input);
	while(input[0]=='#')
	  {
	    FGETS(input,LINESIZE,data->infile);//read rest of line and discard
	    read_word(data->infile, input);    //read first word on next line
	  }
        data->numind[pop][0] = atol(input);//set first numind value, this must be always present
	data->numalleles[pop][0] = data->numind[pop][0] * genomes;
        for(locus=1; locus < data->loci; locus++)
        {
	  read_word(infile, input);// if we encounter a # we treat the rest of the line as comment
	    while(input[0]=='#')
	      {
		FGETS(input,LINESIZE,data->infile);
		read_word(data->infile, input);
	      }
            if(isdigit(input[0]) && havepopname == FALSE )
	      {
                data->numind[pop][locus] = atol(input);
		data->numalleles[pop][locus] = data->numind[pop][locus] * genomes;
	      }
            else
	      {
		// encountered a letter and assume this is the population name
                unread_word(infile, input);
                FGETS(input,LINESIZE,infile);
                havepopname=TRUE;
		if(input!=NULL)
		  {
		    minlength = strlen(input);
		    minlength = MIN(minlength,80);
		    strncpy(data->popnames[pop],input,minlength);
		  }
                break;//this leaves the numind empty of not specified, the must be filled in later
	      }
        }
        if(!havepopname)
	  {
            read_word(infile, input);
	    if(input!=NULL)
	      {
                unread_word(infile, input);
                FGETS(input,LINESIZE,infile);
		minlength = strlen(input);
		minlength = MIN(minlength,80);
		strncpy(data->popnames[pop],input,minlength);
	      }
	  }
	
        // fills numind for additional locus in case the numind was not specified
        for(lo=locus; lo < data->loci; lo++)
	  {
            data->numind[pop][lo] = data->numind[pop][locus-1];
	    data->numalleles[pop][lo] = data->numind[pop][lo] * genomes;
	  }
    }
    else
      {
        // only one locus so we can use old scheme [see below]
        FGETS(input,LINESIZE,infile);
	while(input[0]=='#')
	  {
	    FGETS(input,LINESIZE,data->infile);
	  }
        sscanf (input, "%ld%[^\n]", &(data->numind[pop][0]), data->popnames[pop]);
	data->numalleles[pop][0] = data->numind[pop][0] * genomes;
      }
    
    translate (data->popnames[pop], ' ', '_');
    translate (data->popnames[pop], '\t', '_');
    unpad(data->popnames[pop],"_");
    myfree(input);
}


void
read_indname (FILE * file, data_fmt * data, long pop, long ind, long locus, long nmlength)
{
    long i = 0;
    char ch;
    char input[LINESIZE];
    while (i < nmlength)
    {
        ch = getc (file);
	while(ch =='#')
	  {
	    FGETS(input,LINESIZE,data->infile);
	    ch = getc (file);
	  }
        if(!strchr("\r\n",ch))
            data->indnames[pop][ind][locus][i++] = ch;
        if(strchr("\t",ch))
            break;
    }
    data->indnames[pop][ind][locus][nmlength] = '\0';
}

void
read_popdata (FILE * infile, data_fmt * data, long pop, option_fmt * options)
{
    long ind, baseread = 0;
    long locus = 0;
    if(options->datatype=='h')
      {
	read_hapmap(infile, data, options, locus, pop);
      }
    else
      {
	for (ind = 0; ind < data->numind[pop][0]; ind++)
	  {
	    read_indname (infile, data, pop, ind, locus, options->nmlength);
	    switch (options->datatype)
	      {
	      case 'a':
	      case 'b':
	      case 'm':
		if (data->dlm == '\0')
		  read_alleles (infile, data, pop, ind);
		else
		  read_microalleles (infile, data, options, pop, ind);
		break;
	      case 's':
	      case 'n':
	      case 'u':
	      case 'f':
		if(data->oneliner)
		  {
		    baseread = read_ind_seq_oneliner (infile, data, options, pop, ind, 0);
		  }
		else
		  {
		    baseread = read_ind_seq (infile, data, options, locus, pop, ind, 0);
		  }
		break;
	      default:
		usererror
		  ("Wrong datatype, only the types a, m, s, n\n       (electrophoretic alleles, \n       microsatellite data,\n       sequence data,\n       SNP polymorphism)\n        are allowed.\n");
		break;
	      }
	  }
	if (!strchr (SEQUENCETYPES, options->datatype))
	  return;
	else
	  {
	    if(!data->oneliner)
	      finish_read_seq (infile, data, options, pop, baseread);
	  }
      }
}

long check_list(long la1, long nla1, long locus, data_fmt *data)
{
  if(la1 > data->numrepeatnumbers[locus])
    {
      data->repeatnumbers[locus] = (long *) realloc(data->repeatnumbers[locus], sizeof(long) * (la1+1));
      memset(data->repeatnumbers[locus]+data->numrepeatnumbers[locus], 0, sizeof(long)*(la1+1-data->numrepeatnumbers[locus]));
      data->numrepeatnumbers[locus]= la1+1;
      data->repeatnumbers[locus][la1]=nla1;
    }
  else
    {
      if(data->repeatnumbers[locus][la1]!=0)
	return data->repeatnumbers[locus][la1];
    }
  return la1;
}

void len2repeat(char *a1, long rlen, long locus, data_fmt *data)
{
  long la1 = (long) atol(a1);
  long a=0;
  if((a=la1 % rlen) != 0)
    {
      //      long correction = (((((float) rlen) / 2.) < a) ? (-a) : (rlen/2 == a ? (UNIF_RANDUM()<0.5 ? -a : a) : (rlen - a))); 
      long correction = (((((float) rlen) / 2.) < a) ? (-a) : (rlen/2 == a ? (UNIF_RANDUM()<0.5 ? -a : a) : (rlen - a))); 
      //long nla = check_list(la1,la1+correction, locus, data);
      long nla = (la1+correction)/rlen;
      sprintf(a1,"%li",nla);
    }
  else
    {
      sprintf(a1,"%li",la1/rlen);
    }
}


void
read_microalleles (FILE * infile, data_fmt * data, option_fmt *options, long pop, long ind)
{
    char *input, *isave, dlm[2], ddlm[2], *p, *a, *a1, *a2;
    long locus, i;
    input = (char *) mycalloc (1, sizeof (char) * (SUPERLINESIZE + 1));
    isave = input;
    a = (char *) mycalloc (1, sizeof (char) * LINESIZE);
    a1 = (char *) mycalloc (1, sizeof (char) * LINESIZE);
    a2 = (char *) mycalloc (1, sizeof (char) * LINESIZE);
    dlm[0] = data->dlm, dlm[1] = '\0';
    ddlm[0] = ' ', ddlm[1] = '\0';
    FGETS (input, SUPERLINESIZE, infile);
    if ((p = (char *) strpbrk (input, CRLF)) != NULL)
        *p = '\0';
    for (locus = 0; locus < data->loci; locus++)
    {
        while (isspace ((int) *input))
            input++;
        if (input[0] == '\0')
            FGETS (input, SUPERLINESIZE, infile);
        i = 0;
        while (strchr(" \t",input[i])==NULL && input[i] != dlm[0])
        {
            a1[i] = input[i];
            i++;
        }
        a1[i] = '\0';
        input += i;
        i = 0;
        if (input[i] == dlm[0])
        {
            input++;
            while (strchr(" \t",input[i])==NULL && input[i] != '\0')
            {
                a2[i] = input[i];
                i++;
            }
            a2[i] = '\0';
            if (a2[0] == '\0')
            {
                strcpy (a2, a1);
            }
            input += i;
        }
        else
        {
            strcpy (a2, a1);
        }
        strncpy (data->yy[pop][ind][locus][0], a1, options->allelenmlength);
        strncpy (data->yy[pop][ind][locus][1], a2, options->allelenmlength);
    }
    myfree(a);
    myfree(a1);
    myfree(a2);
    myfree(isave);
}

void
read_alleles (FILE * infile, data_fmt * data, long pop, long ind)
{
    char *input, *isave, *p, *a;
    long locus;
    a = (char *) mycalloc (1, sizeof (char) * LINESIZE);

    input = (char *) mycalloc (1, sizeof (char) * SUPERLINESIZE);
    isave = input;
    FGETS (input, SUPERLINESIZE, infile);
    if ((p = (char *) strpbrk (input, CRLF)) != NULL)
        *p = '\0';
    for (locus = 0; locus < data->loci; locus++)
    {
        while (isspace ((int) *input))
        {
            input++;
        }
        if (sscanf (input, "%s", a) == 1)
        {
            input += (long) strlen (a);
        }

        data->yy[pop][ind][locus][0][0] = a[0];
        data->yy[pop][ind][locus][0][1] = '\0';
        if (a[1] == '\0')
        {
            data->yy[pop][ind][locus][1][0] = a[0];
            data->yy[pop][ind][locus][1][1] = '\0';
        }
        else
        {
            data->yy[pop][ind][locus][1][0] = a[1];
            data->yy[pop][ind][locus][1][1] = '\0';
        }
    }
    myfree(a);
    myfree(isave);
}

void
read_hapmap (FILE * infile, data_fmt * data, option_fmt * options,
              long locus, long pop)
{
  long ind;
    char label1;
    char label2;
    long label1count;
    long label2count;
    long label12count;
    char *input;
    long genomes = number_genomes(options->datatype);
    input = (char *) mycalloc(LINESIZE,sizeof(char));
    for(locus=0;locus<data->loci;locus++)
      {
	read_word(infile,input);
	if(input[0]=='#')
	  {
	    fprintf(stdout,"%s\n",input);
	    FGETS(input,LINESIZE,infile);
	    locus--;
	    continue;
	  }
	data->position[locus] = atol(input);
	if(options->printdata)
	  {
	    fprintf(stdout,"%20li|",data->position[locus]);
	  }
	data->seq[0]->sites[locus]=1;
	read_word(infile,input);
	label1 = input[0];;
	read_word(infile,input);
	label1count = atol(input);
	
	read_word(infile,input);
	label2 = input[0];
	read_word(infile,input);
	label2count = atol(input);
	read_word(infile,input);
	label12count = atol(input);
	if(options->printdata)
	  {
	    fprintf(stdout, " freq[%c]=%f freq[%c]=%f\n",label1,
		    (float) label1count/label12count,label2, (float) label2count/label12count);
	  }
	//was messing with the reading, use fraction data->numind[pop][locus] = label12count;
	//if(options->randomsubset > 0 && options->randomsubset < data->numind[pop][locus])
	//  data->numalleles[pop][locus] = options->randomsubset * genomes;
	//else
	data->numalleles[pop][locus] = data->numind[pop][locus] * genomes;
	label1count = floor( ((label1count/label12count * data->numalleles[pop][locus]))+0.5 );
	label12count = data->numalleles[pop][locus] - label1count;
	for(ind=0;ind < label1count; ind++)
	  {
	    data->yy[pop][ind][locus][0][0] = label1;
	  }
	for(ind=label1count;ind < label12count; ind++)
	  {
	    data->yy[pop][ind][locus][0][0] = label2;
	  }
      }
    myfree(input);
}

long
read_ind_seq (FILE * infile, data_fmt * data, option_fmt * options,
              long locus, long pop, long ind, long baseread)
{
    long j;
    char charstate;
    j = (options->interleaved) ? baseread : 0;
    charstate = getc (infile);
    ungetc ((int) charstate, infile);
    if(options->printdata)
      {
    	fprintf(stdout,"%s|",data->indnames[pop][ind][locus]);
      }
    while (j < data->seq[0]->sites[locus]
            && !(options->interleaved && strchr(CRLF,charstate)))
    {
        charstate = getc (infile);
        if (strchr(CRLF,charstate))
        {
#ifdef INTERLEAVED
            if (options->interleaved)
            {
                while(strchr(CRLF,charstate=getc(infile)))
                    ;
                ungetc ((int) charstate, infile);
                return j;
            }
            else
#endif
                charstate = ' ';
        }
        if (charstate == ' '
                || (charstate >= '0' && charstate <= '9') || charstate == '\\')
            continue;
        charstate = uppercase (charstate);
	if(options->printdata)
	{
	  fprintf(stdout, "%c",charstate);
	}
        if ((strchr ("ABCDGHKMNRSTUVWXY?O-", (int) charstate)) == NULL)
        {
            printf
            ("ERROR: BAD BASE: %c AT POSITION %5ld OF INDIVIDUUM %3li in POPULATION %ld\n",
             charstate, j, ind+1, pop+1);
            printf
            ("Last complete sequences was in population %li, individual %li and locus %li:\n%s",
             pop + 1, ind, locus+1, data->indnames[pop][ind - 1][locus]);
            for (j = 0; j < data->seq[0]->sites[locus]; j++)
                printf ("%c", data->yy[pop][ind - 1][locus][0][j]);
            exit (EXIT_FAILURE);
        }
        data->yy[pop][ind][locus][0][j++] = charstate;
    }
    charstate = getc (infile); /* swallow the \n or \r */
    while(charstate == '\r' || charstate == '\n' || charstate == '\t' || charstate == ' ' || charstate == ';')
      {
	charstate = getc(infile);
      }
    if(charstate!='\n')
      ungetc((int) charstate, infile);

    if(options->printdata)
     {
       fprintf(stdout,"\n");
     }
    return j;
}

///
/// reads sequence data that comes where all loci are (1) on one line and (2) are not interleaved
long
read_ind_seq_oneliner (FILE * infile, data_fmt * data, option_fmt * options,
              long pop, long ind, long baseread)
{
    long j = 0;
    long locus=0;
    char charstate;
    charstate = getc (infile);
    ungetc ((int) charstate, infile);
    if(options->printdata)
      {
    	fprintf(stdout,"%s",data->indnames[pop][ind][locus]);
      }
    for(locus = 0; locus < data->loci; locus++)
      {
	j=0;
	if(options->printdata)
	  {
	    fprintf(stdout,"|");
	  }
	while (j < data->seq[0]->sites[locus])
	  {
	    charstate = getc (infile);
	    if (strchr(CRLF,charstate))
	      {
		charstate = ' ';
	      }
	    if (charstate == ' '
                || (charstate >= '0' && charstate <= '9') || charstate == '\\')
            continue;
	    charstate = uppercase (charstate);
	    if(options->printdata)
	      {
		fprintf(stdout, "%c",charstate);
	      }
	    if ((strchr ("ABCDGHKMNRSTUVWXY?O-", (int) charstate)) == NULL)
	      {
		printf
		  ("ERROR: BAD BASE: %c AT POSITION %5ld OF INDIVIDUUM %3li in POPULATION %ld\n",
		   charstate, j, ind+1, pop+1);
		printf
		  ("Last complete sequences was in population %li, individual %li and locus %li:\n%s",
		   pop + 1, ind, locus+1, data->indnames[pop][ind - 1][locus]);
		for (j = 0; j < data->seq[0]->sites[locus]; j++)
		  printf ("%c", data->yy[pop][ind - 1][locus][0][j]);
		exit (EXIT_FAILURE);
	      }
	    data->yy[pop][ind][locus][0][j++] = charstate;
	  }
      }
    charstate = getc (infile); /* swallow the \n or \r */
    while(charstate == '\n' || charstate == '\t' || charstate == ' ' || charstate == ';')
      {
	charstate = getc(infile);
      }
    if(charstate!='\n')
      ungetc((int) charstate, infile);
    if (charstate == '\r')
    {
        charstate = getc (infile); /* swallow the \n */
        if(charstate!='\n')
            ungetc((int) charstate, infile);
    }
    if(options->printdata)
     {
       fprintf(stdout,"\n");
     }
    return j;
}

void
read_distance_fromfile (FILE * dfile, long tips, long nmlength, MYREAL **m)
{
    char input[SUPERLINESIZE];
    long i, j;

    if (dfile != NULL)
    {
        // assumes that the dfile is in PHYLIP format
        // and that all n x n cells are filled.
        FGETS (input, LONGLINESIZE, dfile); //reads header line with
        for (i = 0; i < tips; i++) // of individuals
        {
            //reads the population label
            FGETS (input, nmlength + 1, dfile);
            for (j = 0; j < tips; j++)
            {
#ifdef USE_MYREAL_FLOAT
                fscanf (dfile, "%f", &m[i][j]);
#else
                fscanf (dfile, "%lf", &m[i][j]);
#endif
		if((i!=j) && (m[i][j] < EPSILON))
		  {
		    warning("Reading geofile: adjusting dist[%li][%li]=%f to %f\n",i,j,m[i][j],EPSILON);
		    m[i][j] = EPSILON;
		  }
            }
            // reads the last \n from the
            // data matrix
            FGETS (input, LONGLINESIZE, dfile);
        }
    }
}

#ifdef UEP
// uep function

void
read_uep_fromfile (FILE * uepfile, long tips, long nmlength, int **uep,
                   long *uepsites, long datatype)
{
    char input[LINESIZE];
    long i, j;
    long thistips;
    if (uepfile != NULL)
    {
        // assumes that the uepfile is in PHYLIP format
        // and that all n cells are filled.
        FGETS (input, LINESIZE, uepfile); //reads header line with
        // of individuals and uep sites
        sscanf (input, "%li%li", &thistips, uepsites);
        if (thistips != tips)
            error ("UEP datafile and infile are inconsistent!");
        if (strchr (SEQUENCETYPES, datatype))
        {
            for (i = 0; i < tips; i++)
            {
                uep[i] = (int *) mycalloc (*uepsites, sizeof (int));
                FGETS (input, nmlength + 1, uepfile); //reads each line
                for (j = 0; j < *uepsites; ++j)
                    fscanf (uepfile, "%i", &uep[i][j]);
                // reads the last \n from the data matrix
                FGETS (input, LINESIZE, uepfile);
            }
        }
        else
        {
            for (i = 0; i < tips; i++)
            {
                uep[i] = (int *) mycalloc (*uepsites, sizeof (int));
                uep[i + tips] = (int *) mycalloc (*uepsites, sizeof (int));
                FGETS (input, nmlength + 1, uepfile); //reads each line
                for (j = 0; j < *uepsites; ++j)
                    fscanf (uepfile, "%i", &uep[i][j]);
                // finished reading first allele, no onto the second
                for (j = 0; j < *uepsites; ++j)
                    fscanf (uepfile, "%i", &uep[i + tips][j]);
                // reads the last \n from the data matrix
                FGETS (input, LINESIZE, uepfile);
            }
        }
    }
}
#endif

void
finish_read_seq (FILE * infile, data_fmt * data, option_fmt * options,
                 long pop, long baseread)
{

    long ind, baseread2 = 0, locus = 0;
    if (options->interleaved)
    {
        while (baseread < data->seq[0]->sites[0])
        {
            for (ind = 0; ind < data->numind[pop][0]; ind++)
            {
                baseread2 =
                    read_ind_seq (infile, data, options, locus, pop, ind,
                                  baseread);
            }
            baseread = baseread2;
        }
    }
    for (locus = 1; locus < data->loci; locus++)
    {
        baseread = 0;
        for (ind = 0; ind < data->numind[pop][locus]; ind++)
        {
	  read_indname (infile, data, pop, ind, locus, options->nmlength);
            baseread = read_ind_seq (infile, data, options, locus, pop, ind, 0);
        }
        if (options->interleaved)
        {
            while (baseread < data->seq[0]->sites[locus])
            {
                for (ind = 0; ind < data->numind[pop][locus]; ind++)
                {
                    baseread2 =
                        read_ind_seq (infile, data, options, locus, pop, ind,
                                      baseread);
                }
                baseread = baseread2;
            }
        }
    }
}

long find_missing(data_fmt *data, long pop, long locus)
{
    long missing = 0;
    long ind;
    for(ind=0; ind < data->numind[pop][locus]; ind++)
    {
        if(data->yy[pop][ind][locus][0][0]=='?')
            missing++;
        if(data->yy[pop][ind][locus][1][0]=='?')
            missing++;
    }
    return missing;
}

void
print_data_summary (FILE * file, world_fmt * world, option_fmt * options,
                    data_fmt * data)
{
  int maxlength = 24; // length of the the total string , for popnames
  int length;
  long locus;
  long pop;
  long numind;
  long nummiss;
  char dstring[LINESIZE];
  long *total;
  long *totalmiss;
  total = (long *) mycalloc(data->loci,sizeof(long));
  totalmiss = (long *) mycalloc(data->loci,sizeof(long));
  fprintf (file, "\nSummary of data:\n");
  fprintf(file, "Title:%54.54s\n",options->title);
  
  switch (options->datatype)
    {
    case 'a':
      strcpy (dstring, "Allelic data");
      break;
    case 'f':
    case 's':
      strcpy (dstring, "Sequence data");
      break;
    case 'b':
    case 'm':
      strcpy (dstring, "Microsatellite data");
      break;
    case 'n':
    case 'u':
      strcpy (dstring, "SNP data");
      break;
    case 'h':
      strcpy (dstring, "SNP data (Hapmap data)");
      break;
    default:
      strcpy (dstring, "Unknown data [ERROR]");
      break;
    }
  fprintf (file, "Data file:  %48.48s\n", options->infilename);
  fprintf (file, "Datatype:   %48.48s\n", dstring);
  if(!data->has_repeats)
    {
      fprintf (file, "[Fragment length is translated to repeats]\n");
    }

  fprintf (file, "Number of loci:                         %20li\n\n",
	   data->loci);
  if(options->has_datefile)
    {
      fprintf (file, "Title:%54.54s\n",options->title);
      fprintf (file, "Sample dates:%46.46s\n", options->datefilename);
      fprintf (file, "Generations per year:                        %5.5f\n", options->generation_year);
      fprintf (file, "Mutationrate per year:  %.10g", options->mutationrate_year[0]);
      for(locus=1; locus < options->mutationrate_year_numalloc; locus++)
	{
	  if(locus % 4 == 0)
	    fprintf(file,"                      ");
	  fprintf(file,", %f", options->mutationrate_year[locus]); 
	}
      fprintf(file,"\n"); 
    }
  for (pop = 0; pop < data->numpop; pop++)
    {
      length = strlen(data->popnames[pop]);
      if(maxlength < length)
	 maxlength = length;
    }
    if(maxlength > 40)
      maxlength = 40;
  
  if (!strchr (SEQUENCETYPES, options->datatype))
    {
      fprintf (file,"%-*.*s     Locus   Gene copies    \n",maxlength,maxlength,"Population");
      fprintf (file,"%-*.*s             ---------------\n",maxlength, maxlength, " ");
      fprintf (file,"%-*.*s             data  (missing)\n",maxlength, maxlength, " ");
    }
  else
    {
      fprintf (file,"%-*.*s     Locus   Gene copies    \n",maxlength,maxlength,"Population");
    }
  fprintf (file,"----%-*.*s------------------------\n",maxlength,maxlength,"------------------------------------------------------------------");
  for (pop = 0; pop < data->numpop; pop++)
    {
        for(locus=0; locus< data->loci; locus++)
        {
	  if (!strchr (SEQUENCETYPES, options->datatype))
	    {
	      nummiss = find_missing(data,pop,locus);
	      numind = data->numalleles[pop][locus] - nummiss;
	      fprintf (file, "%3li %-*.*s %5li %6li (%li)\n", options->newpops[pop], maxlength, maxlength,
		       (locus==0 ? data->popnames[pop] : " " ), locus+1 , numind, nummiss);
	    }
	  else
	    {
	      nummiss = 0;
	      numind = data->numind[pop][locus];
	      fprintf (file, "%3li %-*.*s %5li    %6li\n", options->newpops[pop], maxlength, maxlength,
		       (locus==0 ? data->popnames[pop] : " "), locus+1 , numind);
	    }
	  total[locus] += numind;
	  totalmiss[locus] += nummiss;
        }
    }
    if (!strchr (SEQUENCETYPES, options->datatype))
      {
        for(locus=0; locus< data->loci; locus++)
	  {
	    fprintf (file,"    %-*.*s %5li %6li (%li)\n",maxlength,maxlength, 
		     (locus == 0 ? "Total of all populations" : " "), locus+1, total[locus], totalmiss[locus]);
	  }
      }
    else
      {
        for(locus=0; locus< data->loci; locus++)
	  {
	    fprintf (file,"    %-*.*s %5li    %6li\n",maxlength,maxlength, 
		     (locus == 0 ? "Total of all populations" : " "), locus+1, total[locus]);
	  }
      }    
    fprintf(file,"\n");
    myfree(total);
    myfree(totalmiss);
    fflush (file);
}

void
print_data (world_fmt * world, option_fmt * options, data_fmt * data)
{
    if (options->printdata)
    {
        switch (options->datatype)
        {
        case 'a':
	  //	  if(options->dlm=='\0')
          //  print_alleledata (world, data, options);
	  //else
	  print_microdata (world, data, options);
	  break;
        case 'b':
        case 'm':
            print_microdata (world, data, options);
            break;
        case 's':
        case 'n':
	case 'h':
        case 'u':
        case 'f':
            print_seqdata (world, options, data);
            break;
        }
    }
}

///
/// calculate allele frequency spectrum and print
/// allele population1 .. populationN All
void print_spectra(world_fmt * world, option_fmt * options,data_fmt * data)
{
  long found;
  long locus;
  long a;
  long pop;
  long pop1;
  long ind;
  MYREAL **total;
  MYREAL allfreq;
  long *maxalleles;
  long *maxallelepop;
  MYREAL ***freq;
  char *thisallele;
  char *thatallele;
  FILE *outfile = world->outfile;
  if (strchr (SEQUENCETYPES, options->datatype))
    return; /*we do not calculate allele frequencies for DNA sequences*/
  // find total number of alleles
  maxalleles = (long *) mycalloc(data->loci, sizeof(long));
  maxallelepop = (long *) mycalloc(data->numpop, sizeof(long));
  for (locus = 0; locus < data->loci; locus++)
    {
      maxalleles[locus] = findAllele(data,"\0",locus);
#ifdef BEAGLE
      world->mutationmodels[world->sublocistarts[locus]].numstates = maxalleles[locus];
      world->mutationmodels[world->sublocistarts[locus]].basefreqs = (double *) calloc(world->mutationmodels[world->sublocistarts[locus]].numstates,sizeof(double));
#endif 
    }

  // create bins for each population and all
  freq = (MYREAL ***) mycalloc(data->numpop, sizeof(MYREAL **));
  doublevec2d(&total,data->numpop,data->loci);
  for (pop = 0; pop < data->numpop; pop++)
    {
      freq[pop] = (MYREAL **) mycalloc(data->loci, sizeof(MYREAL *));
      for (locus = 0; locus < data->loci; locus++)
	{
	  freq[pop][locus] = (MYREAL *) mycalloc(1+maxalleles[locus], sizeof(MYREAL));
	}
    }
  
  // calculate spectrum
  for (pop1 = 0; pop1 < data->numpop; pop1++)
    {
      pop = options->newpops[pop1]-1;
      for (ind = 0; ind < data->numind[pop1][0]; ind++)
        {
	  for (locus = 0; locus < data->loci; locus++)
	    {
	      thisallele = data->yy[pop1][ind][locus][0];
	      found = findAllele(data, thisallele, locus);
	      thatallele = data->yy[pop1][ind][locus][1];
	      if(!strcmp(thisallele,thatallele))
		{
		     if(!strchr(thisallele,'?'))
		       {
			 freq[pop][locus][found] += 2.0 ;
			 total[pop][locus] += 2.0 ;
		       }
		}
	      else
		{
		  if(!strchr(thisallele,'?'))
		    {
		      freq[pop][locus][found] += 1.0 ;
		      total[pop][locus] += 1.0 ;
		    }
		  found = findAllele(data, thatallele, locus);
		  if(!strchr(thatallele,'?'))
		    {
		      freq[pop][locus][found] += 1.0 ;
		      total[pop][locus] += 1.0 ;
		    }
		}
	    }
	}
    }
  for (pop = 0; pop < options->newpops_numpop; pop++)
    {
      for (locus = 0; locus < data->loci; locus++)
	{
	  for(a=0; a < maxalleles[locus]; a++)
	    {
	      freq[pop][locus][a] /= (MYREAL) total[pop][locus];
	    }
	}
    }
  // print in ascii
  //fprintf(stdout,"Allele frequency spectra\n");
  fprintf(outfile,"Allele frequency spectra\n");
  fprintf(outfile,"========================\n\n");
  for (locus = 0; locus < data->loci; locus++)
    {
      fprintf(outfile,"Locus %li\n", locus + 1);
      fprintf(outfile,"Allele  ");
      for (pop1 = 0; pop1 < data->numpop; pop1++)
	{
	  maxallelepop[pop1]=0;
	  pop = options->newpops[pop1]-1;
	  fprintf(outfile,"Pop%-2li  ",pop+1);
	}
      fprintf(outfile,"All\n-------");
      for (pop = 0; pop < data->numpop+1; pop++)
	fprintf(outfile,"-------");
      fprintf(outfile, "\n");
      for(a=0; a < maxalleles[locus]; a++)
	{
	  allfreq = 0.0;
	  fprintf(outfile,"%6s ",data->allele[locus][a]);
	  for (pop1 = 0; pop1 < data->numpop; pop1++)
	    {
	      pop = options->newpops[pop1]-1;
	      if(freq[pop][locus][a]>0.0)
		{
		  maxallelepop[pop] += 1;
		  fprintf(outfile," %1.3f ",freq[pop][locus][a]);
		  allfreq += freq[pop][locus][a];
		}
	      else
		{
		  fprintf(outfile,"   -   ");
		}
	    }
	  fprintf(outfile," %1.3f\n", (MYREAL) allfreq/data->numpop);
#ifdef BEAGLE
	  // debug not ready yet for subloci
	  world->mutationmodels[world->sublocistarts[locus]].basefreqs[a] = allfreq/data->numpop;
#endif
	}
      fprintf(outfile,"Total ");
      for (pop1 = 0; pop1 < data->numpop; pop1++)
	{
	  pop = options->newpops[pop1]-1;
	  fprintf(outfile," %5li ",maxallelepop[pop]);
	}
      fprintf(outfile," %5li\n\n", maxalleles[locus]);
    }
  // printd in PDF
#ifdef PRETTY
  pdf_print_spectra(data, freq, maxalleles);
#endif
  fflush(outfile);
  // cleanup
  for (pop = 0; pop < data->numpop; pop++)
    {
      for (locus = 0; locus < data->loci; locus++)
	{
	  myfree(freq[pop][locus]);
	}
      myfree(freq[pop]);
    }
  //printf("finished with printspectra");
  myfree(freq);
  myfree(maxalleles);
  myfree(maxallelepop);
  free_doublevec2d(total);
}

void
print_alleledata (world_fmt * world, data_fmt * data, option_fmt * options)
{
    long i, pop, ind, locus, mult80;
    for (pop = 0; pop < data->numpop; pop++)
    {
        print_header (world->outfile, pop, world, options, data);
        for (ind = 0; ind < data->numind[pop][0]; ind++)
        {
            fprintf (world->outfile, "%-*.*s ", (int) options->nmlength,
                     (int) options->nmlength, data->indnames[pop][ind][0]);
            mult80 = options->nmlength;
            for (locus = 0; locus < data->loci; locus++)
            {
                mult80 +=
                    1 + (long) (strlen (data->yy[pop][ind][locus][0]));
                if (mult80 >= 80)
                {
                    mult80 = 0;
                    fprintf (world->outfile, "\n");
                    for (i = 0; i < options->nmlength; i++)
                        FPRINTF(world->outfile, " ");
                }
                fprintf (world->outfile, " %c.%c",
                         data->yy[pop][ind][locus][0][0],
                         data->yy[pop][ind][locus][0][1]);
            }
            fprintf (world->outfile, "\n");
        }
        fprintf (world->outfile, "\n");
    }
    fprintf (world->outfile, "\n\n");
    fflush (world->outfile);
}

void
print_microdata (world_fmt * world, data_fmt * data, option_fmt * options)
{
    long i, pop, ind, locus, mult80;
    for (pop = 0; pop < data->numpop; pop++)
    {
        print_header (world->outfile, pop, world, options, data);
        for (ind = 0; ind < data->numind[pop][0]; ind++)
        {
            fprintf (world->outfile, "%-*.*s ", (int) options->nmlength,
                     (int) options->nmlength, data->indnames[pop][ind][0]);
            mult80 = options->nmlength;
            for (locus = 0; locus < data->loci; locus++)
            {
                mult80 +=
                    1 + (long) (strlen (data->yy[pop][ind][locus][0]) +
                    strlen (data->yy[pop][ind][locus][1]));
                if (mult80 >= 80)
                {
                    mult80 = 0;
                    fprintf (world->outfile, "\n");
                    for (i = 0; i < options->nmlength; i++)
                        FPRINTF(world->outfile, " ");
                }
                fprintf (world->outfile, " %s.%-s",
                         data->yy[pop][ind][locus][0],
                         data->yy[pop][ind][locus][1]);
            }
            fprintf (world->outfile, "\n");
        }
        fprintf (world->outfile, "\n");
    }
    fprintf (world->outfile, "\n\n");
    fflush (world->outfile);
}

void
print_seqdata (world_fmt * world, option_fmt * options, data_fmt * data)
{
    long pop, locus;
    for (pop = 0; pop < data->numpop; pop++)
    {
        print_header (world->outfile, pop, world, options, data);
        for (locus = 0; locus < data->loci; locus++)
        {
            print_locus_head (locus, world, options, data);
            print_seq_pop (locus, pop, world, options, data);
        }
    }
    fflush (world->outfile);
}

void
print_header (FILE * outfile, long pop, world_fmt * world,
              option_fmt * options, data_fmt * data)
{
    long i;
    long locus, mult80 = 80;
    char input[LINESIZE];
    fprintf (outfile, "\n%-s", data->popnames[pop]);
    for (i = 0; i < (long) (80 - (long) strlen (data->popnames[pop])); i++)
         fprintf(world->outfile, "-");
    fprintf (outfile, "\n\n");
    if (!strchr (SEQUENCETYPES, options->datatype))
    {
        fprintf (outfile, "%-s  ", (data->loci == 1 ? "locus" : "loci "));
        for (i = 0; i < (options->nmlength - 6); i++)
            fprintf(world->outfile, " ");
        for (locus = 0; locus < data->loci; locus++)
        {
            if (locus * 4 + options->nmlength > mult80)
            {
                mult80 += 80;
                fprintf (outfile, "\n");
                for (i = 0; i < options->nmlength; i++)
                    fprintf (outfile, " ");
            }
            fprintf (outfile, "  %2li", locus + 1);
        }
        fprintf (outfile, "\n%-s\n",
                 strncpy (input, "indiv.", options->nmlength));
    }
}


MYREAL findleastsquare(MYREAL *rawdata, long total, long repeatlength, long shift, MYREAL *startvalue)
{
  long i;
  long j;
  long high;
  long low;
  MYREAL lowvalue=MYREAL_MAX;
  MYREAL highvalue=0.;
  MYREAL sum = 0.0;
  MYREAL value;
  MYREAL oldvalue;
  for(i=0;i<total;i++)
    {
      if(rawdata[i] < lowvalue)
	lowvalue = rawdata[i];
      if(rawdata[i]> highvalue)
	highvalue = rawdata[i];
    }
  high = (long) (highvalue+repeatlength+shift);
  low  = (long) (lowvalue -repeatlength+shift);
  if(low<0){
    low = 0;
  }
  *startvalue = low;
  for(i=0;i<total;i++)
    {
      oldvalue = MYREAL_MAX;
      for(j=low ; j < high; j += repeatlength)
	{
	  value = rawdata[i] - j;
	  value *= value;
	  if(value < oldvalue)
	    {
	      oldvalue = value;
	    }
	}
      sum += oldvalue;
    }
  return sum;
}

void find_allele_repeatlength(data_fmt *data, option_fmt *options, long locus)
{
  long z=0;
  long zz=0;
  long pop;
  long ind;
  long r;
  MYREAL *rawdata;
  char *a1;
  char *a2;
  long total=z;
  MYREAL minLS=MYREAL_MAX;
  MYREAL value;
  long intval;
  MYREAL startvalue= 0.0;
  MYREAL keepstartvalue= 0.0;
  MYREAL leastsquare;
  MYREAL nonfloored;
  MYREAL floored;
  MYREAL diff;
  for (pop = 0; pop < data->numpop; pop++)
    {
      zz += data->numind[pop][locus];
    }
  rawdata = (MYREAL *) calloc(2 * zz, sizeof(MYREAL));

  for (pop = 0; pop < data->numpop; pop++)
    {
      for (ind = 0; ind < data->numind[pop][locus]; ind++)
	{
	  a1 = data->yy[pop][ind][locus][0];
	  a2 = data->yy[pop][ind][locus][1];
	  if (a1[0]!='?')
	    {
	      rawdata[z++] = atof(a1);
	    }
	  if (a2[0]!='?')
	    {
	      rawdata[z++] = atof(a2);
	    }
	}
    }

  total=z;
  minLS=MYREAL_MAX;
  startvalue= 0.0;
  keepstartvalue= 0.0;

  for(r=0;r<data->repeatlength[locus];r++)
    {
      //MYREAL endvalue  = 0.0;
      leastsquare=findleastsquare(rawdata,total,data->repeatlength[locus],r,&startvalue);
      if(leastsquare < minLS)
	{
	  minLS = leastsquare;
	 // keepr = r;
	  keepstartvalue = startvalue;
	}
    }
  for (pop = 0; pop < data->numpop; pop++)
    {
      for (ind = 0; ind < data->numind[pop][locus]; ind++)
	{
	  a1 = data->yy[pop][ind][locus][0];
	  a2 = data->yy[pop][ind][locus][1];
	  if (a1[0]!='?')
	    {
	      value = atof(a1);
	      //Floor[(279.47 - 259)/4 + If[(Random[]) < 1 - FractionalPart[(279.47 - 259)/4 ], 0, 1]], {1000}] // Tally
	      nonfloored = (value - keepstartvalue)/data->repeatlength[locus]; 
	      floored = floor(nonfloored);
	      diff = nonfloored - floored; 
	      intval = (long) floored + (RANDUM() < (1.0-diff) ? 0 : 1); 
	      sprintf(data->yy[pop][ind][locus][0],"%li",intval);
	    }
	  if (a2[0]!='?')
	    {
	      value = atof(a2);
	      nonfloored = (value - keepstartvalue)/data->repeatlength[locus]; 
	      floored = floor(nonfloored);
	      diff = nonfloored - floored; 
	      intval = (long) floored + (RANDUM() < (1.0-diff) ? 0 : 1); 
	      sprintf(data->yy[pop][ind][locus][1],"%li",intval);
	    }
	}
    }
}


void
create_alleles (data_fmt * data, option_fmt *options)
{
  long n=0;
  MYREAL  mean=0.;
  MYREAL  delta=0.;
    long locus, pop, ind;
    long z;
    char a1[DEFAULT_ALLELENMLENGTH];
    char a2[DEFAULT_ALLELENMLENGTH];
    for (locus = 0; locus < data->loci; locus++)
    {
      if(data->repeatlength[locus]!=0)
	find_allele_repeatlength(data, options,locus);
      z = 0;
      for (pop = 0; pop < data->numpop; pop++)
        {
            for (ind = 0; ind < data->numind[pop][locus]; ind++)
            {
                strcpy (a1, data->yy[pop][ind][locus][0]);
                strcpy (a2, data->yy[pop][ind][locus][1]);
                if (!strcmp (a1, a2))
                {
                    addAllele (data, a1, locus, &z);
                }
                else
                {
                    addAllele (data, a1, locus, &z);
                    addAllele (data, a2, locus, &z);
                }
            }
        }
	if(z==0)
	  {
	    data->skiploci[locus] = TRUE;
	    continue;
	  }

        data->maxalleles[locus] = z + 1;
        /* + 1: for all the unencountered alleles */
	if(options->murates_fromdata)
	  {
	    if(options->mu_rates==NULL)
	      {
		options->mu_rates = (MYREAL * ) mycalloc(data->loci,sizeof(MYREAL));
		//		printf("%i> data.c: 1557 murate size %li\n",myID,data->loci * sizeof (MYREAL));
	      }
	    options->mu_rates[locus] = z+1;
	    n = n + 1;
	    delta = options->mu_rates[locus] - mean;
	    mean += delta/n;
	  }
    }
    if(options->murates_fromdata)
      {
	options->muloci = data->loci;
	for (locus=0; locus < data->loci; locus++)
	  {
	    if(!data->skiploci[locus])
	      options->mu_rates[locus] /= mean;
	  }
      }
}

void
addAllele (data_fmt * data, char s[], long locus, long *z)
{
    long found = 0;
    if(!strcmp("?",s))
      return;
    while ((data->allele[locus][found++][0] != '\0')
            && (strcmp (s, data->allele[locus][found - 1])))
        ;
    if (found > (*z))
    {
        strcpy (data->allele[locus][*z], s);
        (*z)++;
    }
}

void
set_numind (data_fmt * data)
{
    long locus, pop;
    for (locus = 1; locus < data->loci; locus++)
    {
        for (pop = 0; pop < data->numpop; pop++)
        {
            data->numind[pop][locus] = data->numind[pop][0];
            data->numalleles[pop][locus] = data->numalleles[pop][0];
        }
    }
}


void
print_seq_pop (long locus, long pop, world_fmt * world, option_fmt * options,
               data_fmt * data)
{
    long ind;
    for (ind = 0; ind < data->numalleles[pop][locus]; ind++)
    {
        print_seq_ind (locus, pop, ind, world, options, data);
    }
}

void
print_seq_ind (long locus, long pop, long ind, world_fmt * world,
               option_fmt * options, data_fmt * data)
{
    long site;
    char blank[2] = " ";
    fprintf (world->outfile, "%-*.*s", (int) options->nmlength,
             (int) options->nmlength, data->indnames[pop][ind][0]);
    fprintf (world->outfile, " %c", data->yy[pop][ind][locus][0][0]);
    for (site = 1; site < data->seq[0]->sites[locus]; site++)
    {
        if ((site) % 60 == 0)
        {
            fprintf (world->outfile, "\n%-*.*s %c", (int) options->nmlength,
                     (int) options->nmlength, blank,
                     data->yy[pop][ind][locus][0][site]);
        }
        else
        {
            if ((site) % 10 == 0)
            {
                fprintf (world->outfile, " ");
            }
            fprintf (world->outfile, "%c", data->yy[pop][ind][locus][0][site]);
        }
    }
    fprintf (world->outfile, "\n");
}


void
print_locus_head (long locus, world_fmt * world, option_fmt * options,
                  data_fmt * data)
{
    char *head;
    head = (char *) mycalloc (1, sizeof (char) * MAX (10, options->nmlength));
    sprintf (head, "Locus %li", locus);
    fprintf (world->outfile, "%-*.*s --------10 --------20 --------30",
             (int) options->nmlength, (int) options->nmlength, head);
    fprintf (world->outfile, " --------40 --------50 --------60\n");

    myfree(head);
}

void
read_geofile (data_fmt * data, option_fmt * options, long numpop)
{
    long i, j, pop;
    long numpop2 = numpop * numpop;
    data->geo = (MYREAL *) mycalloc (1, sizeof (MYREAL) * numpop2);
    data->lgeo = (MYREAL *) mycalloc (1, sizeof (MYREAL) * numpop2);
    if (!options->geo)
    {
        for (i = 0; i < numpop2; i++)
            data->geo[i] = 1.0;
    }
    else
    {
        data->ogeo = (MYREAL **) mycalloc (1, sizeof (MYREAL *) * numpop);
        data->ogeo[0] = (MYREAL *) mycalloc (1, sizeof (MYREAL) * numpop2);
        for (pop = 1; pop < numpop; pop++)
            data->ogeo[pop] = data->ogeo[0] + numpop * pop;
        read_distance_fromfile (data->geofile, numpop, options->nmlength,
                                data->ogeo);
        for (i = 0; i < numpop; i++)
        {
            for (j = 0; j < numpop; j++)
            {
                if(i!=j)
                {
                    data->geo[mm2m (i, j, numpop)] =   1. / data->ogeo[i][j];
                    data->lgeo[mm2m (i, j, numpop)] =  data->ogeo[i][j] > 0.0 ?
                                                       LOG (1. / data->ogeo[i][j]) : -MYREAL_MAX;
                }
            }
        }
    }
}

///
/// read the file with the tip dates and returns the oldest date
MYREAL
read_date_fromfile (FILE * datefile, data_fmt *data, option_fmt *options, long nmlength)
{
    char input[LINESIZE];
    long pop;
    long locus;
    long l;
    long ind;
    MYREAL oldest = 0. ;
    MYREAL youngest = DBL_MAX;
    MYREAL temp;
    char *name;
    boolean backward = TRUE;
    name = (char *) mycalloc(LINESIZE,sizeof(char));
    if (datefile != NULL)
    {
        FGETS (input, LINESIZE, datefile); //title line
	while(input[0]=='#')
	  FGETS(input,LINESIZE,datefile);
	if(input[0]=='F') //this checks the direction of time
	  backward=FALSE; 
	fprintf(stdout,"\n%i> Tip dates from file %s\n----------------------------------------\n", myID, options->datefilename);	
	fprintf(stdout,"Generations per year:            %f\n", options->generation_year);
	for(pop=0;pop<data->numpop;pop++)
	  {
	    FGETS (input, LINESIZE, datefile); //first populations 
	    while(input[0]=='#')
	      FGETS(input,LINESIZE,datefile);
	    for (locus=0; locus < data->loci; locus++)
	      {
		l = (locus >= options->mutationrate_year_numalloc ? options->mutationrate_year_numalloc -1 : locus);
		fprintf(stdout,"Mutation rate of Locus %li: %g\n", locus, options->mutationrate_year[l]);
		for(ind=0; ind < data->numind[pop][locus]; ind++)
		  {
		    FGETS (input, LINESIZE, datefile);
		    while(input[0]=='#')
		      FGETS(input,LINESIZE,datefile);
#ifdef USE_MYREAL_FLOAT
		    sscanf (input, "%s",name);
		    sscanf (input+nmlength, "%f",&temp);
#else
		    sscanf (input, "%s",name);
		    sscanf (input+nmlength, "%lf",&temp);
#endif
		    fprintf(stdout,"Locus %li: Tipdate %*.*s %f %f\n", locus, (int) nmlength, (int) nmlength, input, temp, temp * options->mutationrate_year[l] / options->generation_year);
		    data->sampledates[pop][locus][ind].date = temp;
		    data->sampledates[pop][locus][ind].name = (char *) mycalloc(strlen(name)+1,sizeof(char));
		    strcpy(data->sampledates[pop][locus][ind].name, name);
		    unpad(data->sampledates[pop][locus][ind].name," ");
		    translate(data->sampledates[pop][locus][ind].name,' ', '_');
		    if(oldest < temp)
		      oldest = temp;
		    if(youngest > temp)
		      youngest = temp;
		  }
	      }
	  }
	// are the times forward or backward?
	for(pop=0;pop<data->numpop;pop++)
	  {
	    for (locus=0; locus < data->loci; locus++)
	      {
		for(ind=0; ind < data->numind[pop][locus]; ind++)
		  {
		    if(backward)
		      {
			data->sampledates[pop][locus][ind].date =
			  data->sampledates[pop][locus][ind].date - youngest;
		      }
		    else
		      {
			data->sampledates[pop][locus][ind].date = 
			  (oldest - data->sampledates[pop][locus][ind].date) 
			  + (youngest - oldest);
		      }
		  }
	      }
	  }
    }
    free(name);
    return oldest;
}

///
/// read the file with the tip dates
void
read_datefile (data_fmt * data, option_fmt * options, long numpop)
{
  long locus, pop;
  data->sampledates = (tipdate_fmt ***) mycalloc (data->numpop, sizeof (tipdate_fmt **));
  data->sampledates[0] = (tipdate_fmt **) mycalloc (data->numpop * data->loci, sizeof (tipdate_fmt *));
  for (pop = 1; pop < data->numpop; pop++)
    {
      data->sampledates[pop] = data->sampledates[0] + data->loci * pop;
    }
  for(locus=0;locus<data->loci;locus++)
    {
      for(pop=0;pop < data->numpop; pop++)
	{
	  data->sampledates[pop][locus] = (tipdate_fmt*) mycalloc(data->numind[pop][locus],sizeof(tipdate_fmt));
	}
    }

  if (options->has_datefile)
    {
      data->maxsampledate = read_date_fromfile (data->datefile, data, options, options->nmlength);
      //      printf("%i> in data section maxsampledate=%f\n",myID, data->maxsampledate);
    }
  else
    {
      data->maxsampledate=0.0;
    }
}

#ifdef UEP
void
read_uepfile (data_fmt * data, option_fmt * options, long numpop)
{
    long i;
    long sumtips = 0;

    if (!options->uep)
        return;

    for (i = 0; i < numpop; ++i)
        sumtips += data->numind[i][0];   //Assumes that UEP has the same number of individuals as
    // locus 1 (Is this OK? most dataset with UEP will have 1 locus?????)
    data->uep = (int **) mycalloc (number_genomes (options->datatype) * sumtips,
                                 sizeof (int *));
    read_uep_fromfile (data->uepfile, sumtips, options->nmlength, data->uep,
                       &data->uepsites, options->datatype);
}

#endif
