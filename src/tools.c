/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 H E L P E R     R O U T I N E S 
 
 some math stuff and 
 string and file manipulation routines
 
 
 Peter Beerli started 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2006 Peter Beerli, Tallahassee FL
 
  This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
$Id: tools.c 1805 2011-03-10 21:29:04Z beerli $
-------------------------------------------------------*/
/* \file tools.c */

#include <sys/types.h>
#include <sys/stat.h>

#include "migration.h"
#include "sighandler.h"
#include "data.h"
#include "world.h"
#include "tree.h"
#include "random.h"
#include "broyden.h"
#include "migrate_mpi.h"
#include "options.h"
#include "laguerre.h"
#include "pretty.h"

#include <stdio.h>
#include <math.h>
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

/* external globals*/
extern int myID;
#ifdef CAUTIOUS
extern boolean cautious;
#endif
/* prototypes ------------------------------------------- */
MYREAL lengthof (node * p);
MYINLINE node *crawlback (const node * theNode);
/*node *crawl(node * theNode); */
MYINLINE node *showtop (node * theNode);
void adjust_time (node * theNode, MYREAL tyme);
void adjust_time_all (node * theNode, MYREAL tyme);
void insert_migr_node (world_fmt * world, node * up, node * down,
                       migr_table_fmt * migr_table, long *migr_table_counter);
void children (node * mother, node ** brother, node ** sister);
/* math tools */
MYREAL incompletegamma (MYREAL tx, MYREAL talpha);
MYREAL polygamma (long n, MYREAL z);
void invert_matrix (MYREAL **a, long nsize);
boolean nrcheck (MYREAL **m, MYREAL **tm, MYREAL *v, long nrows, MYREAL *r1,
                 MYREAL *r2, boolean do_newton);
MYREAL rannor (MYREAL mean, MYREAL sd);
char lowercase (int c);
char uppercase (int c);
MYREAL calc_sum (MYREAL *vector, long n);

void gamma_rates (MYREAL *rate, MYREAL *probcat, long categs, char *input);
void calc_gamma (MYREAL alpha, MYREAL *gama, long categs);

MYREAL mylgamma (MYREAL z);

MYREAL logfac (long n);

void onepass_mean_std_start(MYREAL *mean, MYREAL *std, long *n);
void onepass_mean_std_calc(MYREAL *mean, MYREAL *std, long *n, MYREAL x);
void onepass_mean_std_end(MYREAL *mean, MYREAL *std, long *n);

//MYINLINE double fast_exp(double y) ;
//MYINLINE float fast_log (float vval);

MYREAL bezier(MYREAL x0, MYREAL y0, MYREAL x1, MYREAL y1);

/* vector initialization */
void doublevec1d (MYREAL **v, long size);
void doublevec2d (MYREAL ***v, long size1, long size2);
void floatvec2d (float ***v, long size1, long size2);
void charvec1d (char **v, long size1);
void charvec2d (char ***v, long size1, long size2);
void intvec2d (long ***v, long size1, long size2);
void free_doublevec2d (MYREAL **v);
void free_floatvec2d (float **v);
void free_charvec2d (char **v);
void free_intvec2d (long **v);
void setdoublevec1d (MYREAL **v, MYREAL *w, long size);
void add_vector (MYREAL *result, MYREAL *v, long size);

/*filemanipulation */
void init_files (world_fmt * world, data_fmt * data, option_fmt * options);
void exit_files (world_fmt * world, data_fmt * data, option_fmt * options);
void openfile (FILE ** fp, char *filename, char *mode, char *perm);
#ifdef ZNZ
void znzopenfile (znzFile * fp, char *filename, char *mode, int use_compressed);
#endif
long read_savesum (world_fmt * world, option_fmt * options, data_fmt * data);
void write_savesum (world_fmt * world);

/* string manipulation */
void translate (char *text, char from, char to);
long count_words (char *text);
void fprintf2(FILE *file, long filesize, const char *fmt, ...);
void print_line (FILE * outfile, char c, long nn, long flag);
void sprint_line (char *buffer, char c, long nn, long flag);
void add_to_buffer(char *fp, long *bufsize, char **buffer, long *allocbufsize);

/* time reporting */
void get_time (char *nowstr, char ts[]);
/*printing aid */
void print_llike (MYREAL llike, char *strllike);

/* searching and finding*/
boolean find (long i, long *list, long listlen);

/* conversion between the parameter schemes*/

long mstart (long pop, long numpop);
long mend (long pop, long numpop);
long mmstart (long pop, long numpop);
long mmend (long pop, long numpop);
long mm2m (long frompop, long topop, long numpop);
void m2mm (long i, long numpop, long *frompop, long *topop);
long m2mml (long i, long numpop);
long m2mml2 (long i, long topop, long numpop);



/* private functions */
MYREAL alnorm (MYREAL x, int up);
void lu_decomp (MYREAL **m, long *indeks, long nrows);
void lu_substitution (MYREAL **m, long *indeks, MYREAL *v, long nrows);
MYREAL d1mach (long i);
long i1mach (long i);
int dpsifn (MYREAL *x, long *n, long kode, long m, MYREAL *ans, long *nz,
            long *ierr);
MYREAL find_chi (long df, MYREAL prob);
MYREAL probchi (long df, double chi);
MYREAL chisquare (long df, MYREAL alpha);

MYREAL chiboundary (long zeros, long nonzeros, MYREAL alpha);


#ifdef MESS
long *check_collection;
long check_collection_count;
extern long unique_id_global;
#endif

/* vector initialization */
void
doublevec1d (MYREAL **v, long size)
{
    *v = (MYREAL *) mycalloc (size, sizeof (MYREAL));
}

/// allocate a 2D array with MYREALS
void
doublevec2d (MYREAL ***v, long size1, long size2)
{
    long i;
    *v = (MYREAL **) mycalloc (size1, sizeof (MYREAL *));
    (*v)[0] = (MYREAL *) mycalloc (size1 * size2, sizeof (MYREAL));
    for (i = 1; i < size1; i++)
    {
        (*v)[i] = (*v)[0] + size2 * i;
    }
}

/// allocate a 2D array with floats
void
floatvec2d (float ***v, long size1, long size2)
{
    long i;
    *v = (float **) mycalloc (size1, sizeof (float *));
    (*v)[0] = (float *) mycalloc (size1 * size2, sizeof (float));
    for (i = 1; i < size1; i++)
    {
        (*v)[i] = (*v)[0] + size2 * i;
    }
}

///
/// allocates columns and rows for a matrix of strings
/// use freecharvec2d() to free this.
void
charvec2d (char ***v, long size1, long size2)
{
    long i;
    *v = (char **) mycalloc (size1, sizeof (char *));
    (*v)[0] = (char *) mycalloc (size1 * size2, sizeof (char));
    for (i = 1; i < size1; i++)
    {
        (*v)[i] = (*v)[0] + size2 * i;
    }
}
///
/// allocates columns and rows for a matrix of strings
/// use free_intvec2d() to free this.
void
intvec2d (long ***v, long size1, long size2)
{
    long i;
    *v = (long **) mycalloc (size1, sizeof (long *));
    (*v)[0] = (long *) mycalloc (size1 * size2, sizeof (long));
    for (i = 1; i < size1; i++)
    {
        (*v)[i] = (*v)[0] + size2 * i;
    }
}

/// free 2D vector of MYREAL
void
free_doublevec2d (MYREAL **v)
{
    myfree(v[0]);
    myfree(v);
}

/// free 2D vector of floats
void
free_floatvec2d (float **v)
{
    myfree(v[0]);
    myfree(v);
}

void
free_charvec2d (char **v)
{
    myfree(v[0]);
    myfree(v);
}
void
free_intvec2d (long **v)
{
    myfree(v[0]);
    myfree(v);
}

void
setdoublevec1d (MYREAL **v, MYREAL *w, long size)
{
    doublevec1d (v, size);
    memcpy (*v, w, sizeof (MYREAL) * size);
}

void
add_vector (MYREAL *result, MYREAL *v, long size)
{
    long i;
    for (i = 0; i < size; i++)
        result[i] += v[i];
}


MYREAL
logfac (long n)
{
    /* log(n!) values were calculated with Mathematica
       with a precision of 30 digits */
    switch (n)
    {
    case 0:
        return 0.;
    case 1:
      return 0.;
    case 2:
      return 0.693147180559945309417232121458;
    case 3:
      return 1.791759469228055000812477358381;
    case 4:
      return 3.1780538303479456196469416013;
    case 5:
      return 4.78749174278204599424770093452;
    case 6:
      return 6.5792512120101009950601782929;
    case 7:
      return 8.52516136106541430016553103635;
    case 8:
      return 10.60460290274525022841722740072;
    case 9:
        return 12.80182748008146961120771787457;
    case 10:
      return 15.10441257307551529522570932925;
    case 11:
      return 17.50230784587388583928765290722;
    case 12:
      return 19.98721449566188614951736238706;
    case 13:
      return 22.5521638531234228855708498286;
    case 14:
      return 25.1912211827386815000934346935;
    case 15:
      return 27.8992713838408915660894392637;
    case 16:
      return 30.6718601060806728037583677495;
    case 17:
      return 33.5050734501368888840079023674;
    case 18:
      return 36.3954452080330535762156249627;
    case 19:
      return 39.3398841871994940362246523946;
    case 20:
      return 42.3356164607534850296598759707;
    case 21:
      return 45.3801388984769080261604739511;
    case 22:
      return 48.4711813518352238796396496505;
    case 23:
      return 51.6066755677643735704464024823;
    case 24:
      return 54.7847293981123191900933440836;
    case 25:
      return 58.0036052229805199392948627501;
    case 26:
      return 61.2617017610020019847655823131;
    case 27:
      return 64.5575386270063310589513180238;
    case 28:
      return 67.8897431371815349828911350102;
    case 29:
      return 71.2570389671680090100744070426;
    case 30:
      return 74.6582363488301643854876437342;
    default:
        return LGAMMA (n + 1.);
    }
}

///
/// initialize the mean and standard deviation calculations 
/// for the one-pass calculations that allows to throw away 
/// the intermediate results
void onepass_mean_std_start(MYREAL *mean, MYREAL *std, long *n)
{
  (*mean) = 0.0;
  (*std)  = 0.0;
  (*n)    = 0;
}

///
/// calculates the mean and standard deviation using a one-pass filter 
/// that allows to throw away the intermediate results
void onepass_mean_std_calc(MYREAL *mean, MYREAL *std, long *n, MYREAL x)
{
  
  long   nn    = (*n);
  MYREAL mm    = (*mean);
  MYREAL delta = x - mm;
  mm           = mm + delta / nn;
  (*std)       = (*std) + delta * (x - mm);
  (*mean)      = mm;
  (*n)         = nn;
}

///
/// finishes one-pass mean and std calculations
void onepass_mean_std_end(MYREAL *mean, MYREAL *std, long *n)
{
  (*std) = sqrt((*std) / ((*n) - 1));
}

/*garbage collection and recycle*/
#ifdef MESS
boolean find_in_list(long id)
{
  long i;
  for(i=0; i < check_collection_count; i++)
    {
      if(id == check_collection[i])
	return TRUE;
    }
  return FALSE;
}

/// tags the lineage that was added the last time the tree was changed.
///
boolean traverse_check_id(node *theNode, long id)
{
  boolean found = FALSE;
  if(theNode->type=='r')
    { 
      if(theNode->id == id || theNode->next->id == id || theNode->next->next->id == id)
	return TRUE;	
      else
	{
	  found = traverse_check_id(theNode->next->back, id);
	  return found;
	}
    }
  if(theNode->type == 'i')
    {
      if(theNode->id == id || theNode->next->id == id || theNode->next->next->id == id)
	return TRUE;	
      else
	{
	  found = traverse_check_id(theNode->next->back, id);
	  if(found == TRUE)
	    return found;
	  found = traverse_check_id(theNode->next->next->back, id);
	  return found;
	}
    }
  if(theNode->type == 'm')
    {
      if(theNode->id == id || theNode->next->id == id)
	return TRUE;	
      else
	{
	  found = traverse_check_id(theNode->next->back, id);
	  return found;
	}
    }
  if(theNode->type == 't')
    {
      if(theNode->id == id)
	return TRUE;	
      else
	return FALSE;
    }
  error("node with no type");
  return FALSE;
}

void traverse_checker(node *root)
{
  long i;
  long x;
  for(i=0;i<check_collection_count; i++)
    {
      x = check_collection[i];
      if(!traverse_check_id(root,x))
	printf ("Node id %li not found\n",x);
    }
}
#endif /*MESS*/

void start_node_collection(world_fmt *world)
{
  long i;
#ifdef MESS
  unique_id_global=0;
#endif
  world->node_collection_allocated = 2 * HUNDRED  + 2 * world->numpop * MIGRATION_LIMIT;
  world->node_collection = (node **) calloc(world->node_collection_allocated,sizeof(node*));
#ifdef MESS
  check_collection = (long *) calloc(world->node_collection_allocated, sizeof(long));
#endif
  for(i = 0; i < world->node_collection_allocated; i++)
    {
      world->node_collection[i] = (node *) calloc(1, sizeof(node));
#ifdef MESS
      (*world->node_collection[i]).id = unique_id_global++;
#endif
    }
  world->node_collection_count = world->node_collection_allocated;
}

void stop_node_collection(world_fmt *world)
{
  long i;
  for(i=0;i<world->node_collection_count;i++)
    {
      if(world->node_collection[i] != NULL)
	{
	  free(world->node_collection[i]);
	}
    }
  myfree(world->node_collection);
}

void collect_nodelet(world_fmt *world, node *p)
{
#ifdef MESS
  long i;
#endif
  long oldalloc;

  if(world->node_collection_count < world->node_collection_allocated)
    {
      world->node_collection[world->node_collection_count] = p;
#ifdef MESS
      for(i=0;i<check_collection_count;i++)
	{
	  if(p->id == check_collection[i])
	    break;
	}
      if(i==check_collection_count)
	{
	  warning("found node not accounted for: %li",p->id);
	  error("");
	}
      else
	{
	  check_collection_count--;
	  check_collection[i] = check_collection[check_collection_count];
	}
#endif      
      world->node_collection_count++;
      // printf("d%li %c\n",p->id, p->type);
    }
  else
    {
#ifdef DEBUG
      warning("in node collection phase we should not need to allocate new space\n");
#endif
      oldalloc = world->node_collection_allocated;
      world->node_collection_allocated += 300;
      world->node_collection = (node **) realloc(world->node_collection, world->node_collection_allocated * sizeof(node*));
      memset(world->node_collection+oldalloc,0, 300 * sizeof(node*));
      world->node_collection[world->node_collection_count] = p;
#ifdef MESS
      for(i=0;i<check_collection_count;i++)
	{
	  if(p->id == check_collection[i])
	    break;
	}
      if(i==check_collection_count)
	{
	  warning("found node not accounted for: %li",p->id);
	  error("");
	}
      else
	{
	  check_collection_count--;
	  check_collection[i] = check_collection[check_collection_count];
	}
#endif
      world->node_collection_count++;
#ifdef DEBUG
      FPRINTF(stdout, "%i> Heat=%f: Collect_nodelet: recycler increased size to: %li/%li\n", myID, 1./world->heat, world->node_collection_count, world->node_collection_allocated);
#endif
    }
  //printf("\n");
}

///
/// dispense memory for nodes from a list of avialable nodes, allocates more memory 
/// no nodes are available in the list if nodes
node *  dispense_nodelet(world_fmt *world)
{
  long i;
  node *p = NULL;
  if(world->node_collection_count > 0)
    {
      world->node_collection_count--;
      p = world->node_collection[world->node_collection_count];
#ifdef MESS
      check_collection[check_collection_count++] = p->id;
#endif
      world->node_collection[world->node_collection_count] = NULL;
    }
  else
    {
      //all node memory is checked out, we
      //increase available nodes in chunks of 300
      world->node_collection_allocated += 300;
      world->node_collection = (node **) myrealloc(world->node_collection,sizeof(node *)* world->node_collection_allocated);
      // set new slots to NULL
      memset(world->node_collection + world->node_collection_allocated - 300, 0,sizeof(node *)*300);
      // we allocate memory of 300 nodes and start at the 0 element because all
      // elements allocated earlier are checked out and the node_collection array
      // should contain only NULLs
      for(i = 0; i < 300; i++)
	{
	  world->node_collection[i] = (node *) mycalloc(1, sizeof(node));
#ifdef MESS
	  world->node_collection[i]->id = unique_id_global++;
#endif
	}
#ifdef MESS
      check_collection = (long *) realloc(check_collection, world->node_collection_allocated * sizeof(long));
#endif
      world->node_collection_count += 300;
      world->node_collection_count--;
      p = world->node_collection[world->node_collection_count]; // dispense
      world->node_collection[world->node_collection_count] = NULL; // remove from list
#ifdef MESS
      check_collection[check_collection_count++] = p->id;
#endif
    }
  return p;
}

///
/// swaps the nodelet recycler between different temperatured chains
/// 
void swap_node_collection(world_fmt * tthis, world_fmt * that)
{
  node **temp;
  //long tempalloc;
  long tempcount;
  long tempallocated;
  temp = tthis->node_collection;
  tempcount = tthis->node_collection_count;
  tempallocated = tthis->node_collection_allocated;

  tthis->node_collection = that->node_collection;
  tthis->node_collection_count = that->node_collection_count;
  tthis->node_collection_allocated = that->node_collection_allocated;

  that->node_collection = temp;
  that->node_collection_count = tempcount;
  that->node_collection_allocated = tempallocated;
}

/*FILEMANIPULATION======================================================= */
void
init_files (world_fmt * world, data_fmt * data, option_fmt * options)
{
	long *tempdb;
#ifdef MPI
	long i,j;
#endif
    long maxfilenum = MYMAXFILENUM; //see definitions.h all files have a number and a filename;
    tempdb = (long *) mycalloc(maxfilenum,sizeof(long));
	
    if (myID == MASTER)
    {
#ifdef MPI
        setup_filehandle_db((void *) stdout, world, options,data);
#endif
		if (!options->readsum)
		{
			openfile (&data->infile, options->infilename, "r+",  NULL);
			if (options->usertree)
				openfile (&data->utreefile, options->utreefilename, "r+",  NULL);
			if (options->weights)
				openfile (&data->weightfile, options->weightfilename, "r+",  NULL);
			if (options->categs > 1)
				openfile (&data->catfile, options->catfilename, "r+",  NULL);
			if (options->dist)
				openfile (&data->distfile, options->distfilename, "r+",  NULL);
			if (options->writesum)
				openfile (&data->sumfile, options->sumfilename, "w+",  NULL);
		}
		else
		  {
		    if(!options->bayes_infer)
		      openfile (&data->sumfile, options->sumfilename, "r+",  NULL);
		  }

		openfile (&world->outfile, options->outfilename, "w+",  NULL);
#ifdef MPI
		setup_filehandle_db((void *) world->outfile, world, options,data);
#endif

		if (options->treeprint > 0 && (!options->readsum))
		  {
		    openfile (&world->treefile, options->treefilename, "w+",  NULL);
#ifdef MPI
		    setup_filehandle_db((void *) world->treefile, world, options,data);
#endif
		  }

		if (options->mighist)
		  {
		    if(options->datatype != 'g')
		      {
			openfile (&world->mighistfile, options->mighistfilename, "w+", NULL);
		      }
		    else
		      {
			openfile (&world->mighistfile, options->mighistfilename, "r+", NULL);
		      }
#ifdef MPI
			setup_filehandle_db((void *) world->mighistfile, world, options, data);
#endif			
			if(options->skyline)
			  {
			    if(options->datatype != 'g')
			      {
				openfile (&world->skylinefile, options->skylinefilename, "w+", NULL);
			      }
			    else
			      {
				openfile (&world->skylinefile, options->skylinefilename, "r+", NULL);
			      }
			  }
#ifdef MPI
			setup_filehandle_db((void *) world->skylinefile, world, options, data);
#endif			
		  }
        
		if (options->writelog)
		  {
			openfile (&options->logfile, options->logfilename, "w+",  NULL);
#ifdef MPI
			setup_filehandle_db((void *) options->logfile, world, options, data);
#endif			
		  }
		
		if (options->mixplot)
		  {
		    openfile (&options->mixfile, options->mixfilename, "w+",  NULL);
#ifdef MPI
		    setup_filehandle_db((void *) options->mixfile, world, options, data);
#endif			
		  }
		if (options->aic)
		  {
		    openfile (&options->aicfile, options->aicfilename, "w+",  NULL);
#ifdef MPI
		    setup_filehandle_db((void *) options->aicfile, world, options, data);
#endif
		  }			
		if (options->plot)
		  {
		    switch (options->plotmethod)
		      {
		      case PLOTALL:
			openfile (&world->mathfile, options->mathfilename, "w+", 
				  NULL);
#ifdef MPI
			setup_filehandle_db((void *) world->mathfile, world, options, data);
#endif
			break;
		      default:  /*e.g. 0 this create just the plots in outfile */
			break;
		      }
		  }
		if (options->geo)
		  {
		    openfile (&data->geofile, options->geofilename, "r+",  NULL);
		  }
		if (options->has_datefile)
		  {
		    openfile (&data->datefile, options->datefilename, "r+",  NULL);
		  }
#ifdef UEP
		
		if (options->uep)
		  openfile (&data->uepfile, options->uepfilename, "r+",  NULL);
#endif
		if (options->bayes_infer)
		  {
		    if(options->has_bayesfile)
		      openfile (&world->bayesfile, options->bayesfilename, "w+",  NULL);

		    if(options->has_bayesmdimfile)
		      {
			if(options->datatype != 'g')
			  {
#ifdef MPI
#ifndef PARALIO
#ifdef ZNZ
			    znzopenfile (&world->bayesmdimfile, options->bayesmdimfilename, "w", options->use_compressed);
#else
			    openfile (&world->bayesmdimfile, options->bayesmdimfilename, "w+",  NULL);
#endif
#endif			  
#else
#ifdef ZNZ
			    znzopenfile (&world->bayesmdimfile, options->bayesmdimfilename, "w",  options->use_compressed);
#else
			    openfile (&world->bayesmdimfile, options->bayesmdimfilename, "w+",  NULL);
#endif
#endif
			  }
			else
			  {
#ifdef ZNZ
			    znzopenfile (&world->bayesmdimfile, options->bayesmdimfilename, "r",  options->use_compressed);
#else
			    openfile(&world->bayesmdimfile, options->bayesmdimfilename,"r+",NULL);
#endif
			  }
		      }
#ifdef MPI
		    if(options->has_bayesfile)
		      setup_filehandle_db((void *) world->bayesfile, world, options,data);
		    if(options->has_bayesmdimfile)
		      {
#ifdef PARALIO
			MPI_File_open(comm_world,options->bayesmdimfilename,MPI_MODE_CREATE | MPI_MODE_WRONLY ,MPI_INFO_NULL, &world->mpi_bayesmdimfile);
#else
			setup_filehandle_db((void *) world->bayesmdimfile, world, options,data);
#endif
		      }
#endif
		  }			
#ifdef MPI
		tempdb[0]=filenum;
		for(i=0, j=1;i<filenum;i++)
		  {
		    tempdb[j++] = (long) filedb[i].file;
		    tempdb[j++] = (long) filedb[i].handle;
		  }
		MYMPIBCAST (tempdb, maxfilenum, MPI_LONG, MASTER, comm_world);
#endif
    }
    else
      {
	// all non myID==MASTER
#ifdef MPI
	MYMPIBCAST (tempdb, maxfilenum, MPI_LONG, MASTER, comm_world);
	filenum = tempdb[0]; 
	for(i=0, j=1;i<filenum;i++)
	  {
	    filedb[i].file = (FILE *) tempdb[j++];
	    filedb[i].handle =  tempdb[j++];
	  }
	//assumess same order as in the master
	j=1; //zero is stdout
	world->outfile = filedb[j++].file;
	if (options->treeprint > 0 && (!options->readsum))
	  world->treefile = filedb[j++].file;
	if (options->mighist)
	  world->mighistfile = filedb[j++].file;
	if (options->skyline)
	  world->skylinefile = filedb[j++].file;
	if (options->writelog)
	  options->logfile = filedb[j++].file;
	if (options->aic)
	  options->aicfile = filedb[j++].file;
	if (options->plot && options->plotmethod == PLOTALL)
	  world->mathfile = filedb[j++].file;
        if (options->bayes_infer)
	  {
	    if(options->has_bayesfile)
	      world->bayesfile = filedb[j++].file;
	    // TEST: experiment with parallele I/O
	    if(options->has_bayesmdimfile)
#ifdef PARALIO
	      MPI_File_open(comm_world,options->bayesmdimfilename,MPI_MODE_CREATE,MPI_INFO_NULL, &world->mpi_bayesmdimfile);
#else
#ifdef ZNZ
	    world->bayesmdimfile = (znzFile) filedb[j++].file;
#else
	      world->bayesmdimfile = filedb[j++].file;
#endif
#endif
	  }
#endif /*MPI*/
      }
    myfree(tempdb);
}

void
exit_files (world_fmt * world, data_fmt * data, option_fmt * options)
{
    if(myID==MASTER)
    {
      if (!options->readsum)
	{
	  FClose (data->infile);
	  //  if (options->treeprint > 0)
	  //  FClose (world->treefile);
	  if(options->usertree)
	    FClose(data->utreefile);
	  if (options->weights)
	    FClose (data->weightfile);
	  if (options->categs > 1)
	    FClose (data->catfile);
	  if (options->dist)
	    FClose (data->distfile);
	}
      
      if (options->writesum || options->readsum)
	FClose (data->sumfile);
      
      FClose (world->outfile);
      
      if (options->writelog)
	FClose (options->logfile);
      if (options->mighist)
	FClose (world->mighistfile);
      if (options->geo)
	FClose (data->geofile);
      if (options->has_datefile)
	FClose (data->datefile);
      if (options->aic)
	FClose (options->aicfile);
      if (options->plot && options->plotmethod == PLOTALL)
	FClose (world->mathfile);
#ifdef UEP
      
      if (options->uep)
	FClose (data->uepfile);
#endif
      if (options->bayes_infer)
	{
	  if(options->has_bayesfile)
	    FClose (world->bayesfile);
	  if(options->has_bayesmdimfile)
#ifdef ZNZ
	    znzClose (world->bayesmdimfile);
#else
	    FClose (world->bayesmdimfile);
#endif
	}
    }
}

/* string manipulation ================================== */
/* Converts any character from to character to in string text */
void
translate (char *text, char from, char to)
{
    int i, j, gap = 0;
    while (text[gap] == from)
        gap++;
    for (i = gap, j = 0; text[i] != '\0'; i++)
    {
        if (text[i] != from)
        {
            text[j++] = text[i];
        }
        else
        {
            if (text[i - 1] != from)
            {
                text[j++] = to;
            }
        }
    }
    text[j] = '\0';
}
/// unpad from the end of the word until there is none of the
/// specified character left
void
unpad (char *text, char removechars[])
{
  int gap = strlen(text)-1;
  while (gap > 0 && strchr(removechars,text[gap]))
        gap--;
    text[gap+1] = '\0';
}

///
/// get next word from list with list of delimiters
/// uses strsep
void
get_next_word(char **instring, char *delimiters, char **nextword)
{
  *nextword = strsep(instring, delimiters);
  while(*nextword != NULL && strchr(delimiters,(*nextword)[0]))
  {
    //    if((*instring)==NULL)
    //  {
    //	(*nextword) = NULL;
    //	return;
    //      }
      *nextword = strsep(instring,delimiters);
  }
}
/*===============================================
  count words in a string delimited by delimiter
*/
long
count_words (char *text)
{
    long counts = 0;
    char *pt = text;
    while (isspace (*pt) && *pt != '\0')
        pt++;
    while (*pt != '\0')
    {
        while (!isspace (*pt) && *pt != '\0')
            pt++;
        while (isspace (*pt) && *pt != '\0')
            pt++;
        counts++;
    }
    return counts;
}


/*===============================================
 timer utility
 
 ts = "%c" -> time + full date (see man strftime)
      = "%H:%M:%S" -> time hours:minutes:seconds */

void
get_time (char *nowstr, char ts[])
{
#ifdef NOTIME_FUNC
    switch (strlen (ts))
    {
    case 2:
        strcpy (nowstr, " ");
        break;
    case 3:
        strcpy (nowstr, "  ");
        break;
    case 8:
        strcpy (nowstr, "        ");
        break;
    default:
        strcpy (nowstr, " ");
        break;
    }
#else
    time_t nowbin;
    struct tm *nowstruct;
    if (time (&nowbin) != (time_t) - 1)
    {
        nowstruct = localtime (&nowbin);
        strftime (nowstr, LINESIZE, ts, nowstruct);
    }
#endif
}

/*===============================================
 printer utility
 */
void
print_llike (MYREAL llike, char *strllike)
{
    if (fabs (llike) > 10e20)
    {
        sprintf (strllike, "%cInfinity ", llike < 0 ? '-' : ' ');
    }
    else
        sprintf (strllike, "%-10.5f", llike);
}

///
/// remove traling whitespace from a *string
void remove_trailing_blanks(char **filename)
{
  int pos;
  pos = (long) strlen(*filename);
  while (pos>0 && isspace( (*filename)[pos-1]))
    pos--;
  (*filename)[pos] = '\0';   
}

///
/// read a single word from a string up to the position of 
/// the breakchar and return that position 
long read_word_delim(char *input, char *word, char * delim)
{
  long i = 1;
  int ch;
  long count = 0;
  ch = input[0];
  if(ch=='\0')
    return 0;
  while (isspace (ch))
    ch = input[i++];
  while(strchr(delim, ch)==NULL)
    {
      word[count++] = ch;
      ch = input[i++];
    }
  word[count]='\0';
  return i;
}
///
/// read a single word from a file
/// words are delimited by withespace 
void read_word(FILE *infile, char *word)
{
    int ch;
    long count = 0;
    ch = getc(infile);
    while (isspace (ch))
        ch = getc(infile);
    while(strchr(" \t\n\r", ch)==NULL)
    {
        word[count++] = ch;
        ch = getc(infile);
    }
    word[count]='\0';
    ungetc(ch,infile);
}

//========================================
// prepends a string to a file
// a blank will be inserted between the string and the
// text in the file
void unread_word(FILE *infile, char *word)
{
  char c = ' ';
  long count = (long) strlen(word);
  ungetc(c, infile);
  while(count > 0)
    {
      count--;
      ungetc(word[count], infile);
    }
}

#ifdef STRANGEDEBUG
extern int errno;
#endif

#ifdef ZNZ
/// opens a zip file for reading or writing
void
znzopenfile (znzFile * fp, char *filename, char *mode, int use_compressed)
{
  //int trials = 0;
    znzFile of;
    char *file;
    char *p;
    file = (char *) mycalloc(LINESIZE,sizeof(char));
    if ((p = strpbrk (filename, CRLF)) != NULL)
        *p = '\0';
    remove_trailing_blanks(&filename);
    strncpy (file, filename, LINESIZE - 1);
    of = znzopen (file, mode,use_compressed);
    if (znz_isnull(of))
      error("Could not open zipped bayesallfile");
    *fp = of;
    myfree(file);
}
#endif


/// opens a file for reading or writing
void
openfile (FILE ** fp, char *filename, char *mode, char *perm)
{
    int trials = 0;
    FILE *of = NULL;
    char *file;
    char *p;
#ifdef CAUTIOUS
    struct stat sb;
#endif
    file = (char *) mycalloc(LINESIZE,sizeof(char));
    if ((p = strpbrk (filename, CRLF)) != NULL)
        *p = '\0';
    remove_trailing_blanks(&filename);
    strncpy (file, filename, LINESIZE - 1);
#ifdef CAUTIOUS
    if(cautious)
      {
	//	sb = (stat *) mycalloc(1,sizeof(stat));
	if(strchr(mode,'w') && (0 == stat(file, &sb)))
	  {
	    file[0] = '\0';
	    while (file[0] == '\0' && trials++ < 10)
	      {
		printf ("File %s exists, do you want to overwrite it? [YES/No]\n===>",filename);
		FGETS (file, LINESIZE, stdin);
		if(file[0]=='Y' || file[0] == 'y')
		  {
		    file[0] = '\0';
		    while (file[0] == '\0' && trials++ < 10)
		      {
			printf ("Please enter a new filename>");
			FGETS (file, LINESIZE, stdin);
		      }
		    
		  }
	      }
	  }
	//	myfree(sb);
      }
#endif
    while (trials++ < 10)
    {
      of = fopen (file, mode);
#ifdef STRANGEDEBUG
      fprintf(stdout,"original filename -->%s<--\n",filename);
      fprintf(stdout,"used filename     -->%s<--\n",file);
#endif
        if (of!=NULL)
            break;
        else
        {
#ifdef STRANGEDEBUG
	  fprintf(stdout,"mode              -->%s<--\n", mode);
	  fprintf(stdout,"errorcode         -->%i<--\n",errno);
#endif
            switch (*mode)
	      {
	      case 'r':
#ifdef MPI

                printf("%i> ",myID);
#endif

                printf ("Cannot read from file \"%s\"\n", file);
                file[0] = '\0';
                while (file[0] == '\0' && trials++ < 10)
                {
#ifdef MPI
                    printf("%i> ",myID);
#endif

                    printf ("Please enter a new filename for reading>");
                    FGETS (file, LINESIZE, stdin);
                }
                break;
            case 'w':
#ifdef MPI

                printf("%i> ",myID);
#endif

                printf ("Cannot write to file %s\n", file);
                file[0] = '\0';
                while (file[0] == '\0' && trials++ < 10)
                {
#ifdef MPI
                    printf("%i> ",myID);
#endif

                    printf ("Please enter a new filename for writing>");
                    FGETS (file, LINESIZE, stdin);
                }
                break;
            }
        }
    }
    if (trials >= 10)
    {
#ifdef MPI
        printf("%i> ",myID);
#endif

        printf ("You cannot find your file either, so I stop\n\n");
        exit (0);
    }
    *fp = of;
    if (perm != NULL)
        strcpy (perm, file);
    strcpy (filename, file);
    myfree(file);
}

long
read_savesum (world_fmt * world, option_fmt * options, data_fmt * data)
{
    long Gmax = 0;
    FILE *sumfile = data->sumfile;
    timearchive_fmt **ta;
    long nrep = 0;
    long l, i, j, r;
    char input[1024];
    long hits;
    long tmp;
    long z, zz, jj;
    MYREAL tmpm;
    boolean newstyle = TRUE;
#ifdef MPI

    long listofthings[6];
#endif

    FGETS (input, 1024, sumfile);
    if (!strncmp ("# begin genealogy-summary file of migrate", input, 41))
    {
        FGETS (input, 1024, sumfile);
        if (input[0] == '#')
        {
            newstyle = TRUE;
            FGETS (input, 1024, sumfile);
            if(input[0]=='M') // added custom migration model 11-25-02
              {
                options->custm = (char *) myrealloc(options->custm, sizeof(char) * (strlen(input)+1));
                strcpy(options->custm,input+6);
                FGETS (input, 1024, sumfile);
                options->custm2 = (char *) myrealloc(options->custm2, sizeof(char) * (strlen(input)+1));
                strcpy(options->custm2,input+6);
                FGETS (input, 1024, sumfile);
              }
        }
        else
            newstyle = FALSE;
        hits =
            sscanf (input, "%li %li %li %li %li", &world->loci, &world->numpop,
                    &world->numpop2, &tmp, &options->replicatenum);
        if (hits != 5)
        {
            nrep = 1;
            //xcode hits =
                sscanf (input, "%li %li %li", &world->loci, &world->numpop,
                        &world->numpop2);
            options->replicate = FALSE;
            options->replicatenum = 0;
        }
        else
        {
            nrep = options->replicatenum;
            if (nrep == 0)
                nrep = 1;
            options->replicate = (boolean) tmp;
        }
        options->muloci = data->loci = world->loci;
        data->skiploci =
            (boolean *) myrealloc (data->skiploci,
                                 sizeof (boolean) * (data->loci + 1));
        memset (data->skiploci, 0, sizeof (boolean) * (data->loci + 1));
        read_geofile (data, options, world->numpop);
        data->numpop = world->numpop;
        set_plot (options);
	world->cold=TRUE;
        init_world (world, data, options);
        ta = world->atl;
        for (l = 0; l < world->loci; l++)
            for (r = 0; r < nrep; r++)
            {
                FGETS (input, 1024, sumfile);
                fscanf (sumfile, "%li %li %li\n", &ta[r][l].T, &ta[r][l].numpop,
                        &ta[r][l].sumtips);
                Gmax = (ta[r][l].T > Gmax) ? ta[r][l].T : Gmax;
                ta[r][l].allocT = 0;
                increase_timearchive (world, l, ta[r][l].T, world->numpop, r);
#ifdef USE_MYREAL_FLOAT
                fscanf (sumfile, "%g %g %g\n", &ta[r][l].param_like,
                        &ta[r][l].thb, &ta[r][l].alpha);
#else
                fscanf (sumfile, "%lg %lg %lg\n", &ta[r][l].param_like,
                        &ta[r][l].thb, &ta[r][l].alpha);
#endif
                world->chainlikes[l][r] = ta[r][l].param_like;

                for (i = 0; i < world->atl[r][l].T; i++)
                {
#ifdef USE_MYREAL_FLOAT
                    fscanf (sumfile, "%li %f\n", &ta[r][l].tl[i].copies, &ta[r][l].tl[i].lcopies);
#else
                    fscanf (sumfile, "%li %lf\n", &ta[r][l].tl[i].copies, &ta[r][l].tl[i].lcopies);
#endif
                    for (j = 0; j < world->numpop; j++)
                    {
#ifdef USE_MYREAL_FLOAT
                        fscanf (sumfile, "%f %f %f\n", &ta[r][l].tl[i].km[j],
                                &ta[r][l].tl[i].kt[j], &ta[r][l].tl[i].p[j]);
#else               
                        fscanf (sumfile, "%lf %lf %lf\n", &ta[r][l].tl[i].km[j],
                                &ta[r][l].tl[i].kt[j], &ta[r][l].tl[i].p[j]);
#endif
                    }
                    if (newstyle)
                    {
                        for (j = 0; j < world->numpop2; j++)
                        {
#ifdef USE_MYREAL_FLOAT
                            fscanf (sumfile, "%f ", &ta[r][l].tl[i].mindex[j]);
#else
                            fscanf (sumfile, "%lf ", &ta[r][l].tl[i].mindex[j]);
#endif
                        }
                    }
                    else
                    {
                        z = 0;
                        zz = world->numpop;
                        for (j = 0; j < world->numpop; j++)
                        {
                            for (jj = 0; jj < world->numpop; jj++)
                            {
#ifdef USE_MYREAL_FLOAT
                                fscanf (sumfile, "%f ", &tmpm);
#else
                                fscanf (sumfile, "%lf ", &tmpm);
#endif
                                if (j == jj)
                                    ta[r][l].tl[i].mindex[z++] = tmpm;
                                else
                                    ta[r][l].tl[i].mindex[zz++] = tmpm;
                            }
                        }
                    }
                }
                //     ta[r][l].tl[0].lcopies = ta[r][l].tl[0].copies > 1 ? ln_copies(ta[r][l].tl[0].copies)-1. : 0.0;
                for (i = 0; i < world->numpop2; i++)
                {
#ifdef USE_MYREAL_FLOAT
                    fscanf (sumfile, "%g %g\n", &ta[r][l].param[i],
                            &ta[r][l].param0[i]);
#else
                    fscanf (sumfile, "%lg %lg\n", &ta[r][l].param[i],
                            &ta[r][l].param0[i]);
#endif
                }
                log_param0 (ta[r][l].param0, ta[r][l].lparam0, world->numpop2);
                //for(i=0;i<ta[l].T;i++)
                //  {
                //    fscanf(sumfile,"%lg ",&ta[l].likelihood[i]);
                //  }
                //      FGETS(input,1024,sumfile);
#ifdef USE_MYREAL_FLOAT
                fscanf (sumfile, "%li %g\n", &ta[r][l].trials, &ta[r][l].normd);
#else
                fscanf (sumfile, "%li %lg\n", &ta[r][l].trials, &ta[r][l].normd);
#endif           
            }
    }
    else
    {
        error ("This is not a genealogy-summary file for MIGRATE\n");
    }
#ifdef MPI
    listofthings[0] = world->loci;
    listofthings[1] = world->numpop;
    listofthings[2] = world->numpop2;
    listofthings[3] = (long) options->replicate;
    listofthings[4] = options->replicatenum;
    listofthings[5] = Gmax;
    //printf("listofthings sent");
    fflush(stdout);
    MYMPIBCAST (listofthings,5, MPI_LONG, MASTER, comm_world);
#endif

    return Gmax;
}

void
write_savesum (world_fmt * world)
{
    long r, repmax;
    FILE *sumfile = world->data->sumfile;
    timearchive_fmt **ta = world->atl;

    long l, i, j;
    if (world->options->replicate)
    {
        if (world->options->replicatenum == 0)
            repmax = world->options->lchains;
        else
            repmax = world->options->replicatenum;
    }
    else
        repmax = 1;
    fprintf (sumfile,
             "# begin genealogy-summary file of migrate %s ------\n#\n",
             MIGRATEVERSION);
    fprintf (sumfile,"Model=%s\n",world->options->custm);
    fprintf (sumfile,"Mode2=%s\n",world->options->custm2);
    fprintf (sumfile, "%li %li %li %li %li\n", world->loci, world->numpop,
             world->numpop2, (long) world->options->replicate, repmax);
    for (l = 0; l < world->loci; l++)
    {
        for (r = 0; r < repmax; r++)
        {
            fprintf (sumfile,
                     "%li %li ####### locus %li, replicate %li ################\n",
                     l, r, l, r);
	    fprintf(sumfile,"%li %li %li\n",ta[r][l].T, ta[r][l].numpop, ta[r][l].sumtips);

            fprintf (sumfile, "%20.20g %20.20g %20.20g\n", ta[r][l].param_like,
                     ta[r][l].thb, ta[r][l].alpha);
            for (i = 0; i < ta[r][l].T; i++)
            {
                fprintf (sumfile, "%li %f\n", ta[r][l].tl[i].copies,ta[r][l].tl[i].lcopies);
                for (j = 0; j < world->numpop; j++)
                {
                    fprintf (sumfile, "%20.20f %20.20f %f\n",
                             ta[r][l].tl[i].km[j], ta[r][l].tl[i].kt[j],
                             ta[r][l].tl[i].p[j]);
                }
                for (j = 0; j < world->numpop2; j++)
                {
                    fprintf (sumfile, "%f ", ta[r][l].tl[i].mindex[j]);
                    //debug               printf (sumfile, "%li ", ta[l].tl[i].l[j]);
                }
                fprintf (sumfile, "\n");
            }
            for (i = 0; i < world->numpop2; i++)
            {
                fprintf (sumfile, "%20.20e %20.20e\n", ta[r][l].param[i],
                         ta[r][l].param0[i]);
            }
            //  for(i=0;i<ta[l].T;i++)
            //        {
            //          printf(sumfile,"%20.20e ",ta[l].likelihood[i]);
            //        }
            //      printf(sumfile,"\n");
            fprintf (sumfile, "%li %20.20e\n", ta[r][l].trials, ta[r][l].normd);
        }
    }
    fprintf (sumfile, "# end genealogy-summary file of migrate %s ------\n",
             MIGRATEVERSION);
    fflush (sumfile);
}


/*=======================================================*/




/*--------------------------------
creates the length value in a node
*/
MYREAL
lengthof (node * p)
{
    if (p->type == 'm')
        error ("A migration node was feed into lengthof");
    return fabs (p->tyme - crawlback (p)->tyme);
}    /* length */


/*------------------------------------------------
Find the next non-migration node starting
with the theNode, returns to backnode which is not 
a migration, does NOT return always a top-node!
*/
MYINLINE node *
crawlback (const node * theNode)
{
    node *tmp = theNode->back;

    while (tmp->type == 'm')
    {
        tmp = tmp->next->back;
    }
    return tmp;
}

/*--------------------------------------------
returns the last migration node in a branch or 
the node if there is no migration node
 
node *crawl(node * theNode)
{
   node *otmp, *tmp = theNode->back;
 
   otmp = theNode;
   if (tmp == NULL)
   return otmp;
   while (tmp->type == 'm') {
   otmp = tmp->next;
   tmp = tmp->next->back;
   if (tmp == NULL)
   return otmp;
   }
   return otmp;
}
*/

 
MYINLINE node *
showtop (node * theNode)
{
    if (theNode == NULL)
        return NULL;
    else
    {
        if (theNode->top)
        {
            return theNode;
        }
        else
        {
            if (theNode->next->top)
            {
                return theNode->next;
            }
            else
            {
                return theNode->next->next;
            }
        }
    }

}

/* adjust the time in a node to time */
void
adjust_time (node * theNode, MYREAL tyme)
{
    switch (theNode->type)
    {
    case 'm':
        theNode->tyme = theNode->next->tyme = tyme;
	set_dirty(theNode);
        break;
    case 'i':
        theNode->tyme = theNode->next->tyme = theNode->next->next->tyme = tyme;
	set_dirty(theNode);
        break;
    case 'r':
      break;
    case 't':
      break;
    default:
      error ("Wrong node type");
      break;
    }
}
/* adjust the time in any node to time */
void
adjust_time_all (node * theNode, MYREAL tyme)
{
    switch (theNode->type)
    {
    case 'm':
        theNode->tyme = theNode->next->tyme = tyme;
	set_dirty(theNode);
        break;
    case 'i':
        theNode->tyme = theNode->next->tyme = theNode->next->next->tyme = tyme;
	set_dirty(theNode);
        break;
    case 'r':
      theNode->tyme = theNode->next->tyme = theNode->next->next->tyme = tyme;
      set_dirty(theNode);
      break;
    case 't':
      theNode->tyme = tyme;
      break;
    default:
      error ("Wrong node type");
      break;
    }
}

void
insert_migr_node (world_fmt * world, node * up, node * down,
                  migr_table_fmt * migr_table, long *migr_table_counter)
{
  long i;
  //, panic;
    //node *tmp, *tmp2, *oldNode, *oldNode2, *theNode;
  node *theNode;
    if (!up->top)
        error ("up has to be a top-node");
    //xcode theNode = showtop (up)->back;
    if (*migr_table_counter > 0 && up->tyme > migr_table[0].time)
      {
	printf("%i> up->type=%c tyme=%f showtop(up)time=%f mig[0]tabletyme=%f\n", myID, up->type, up->tyme, showtop(up)->tyme, migr_table[0].time);
        error
	  ("insert_migr_node: the first migration node has a wrong time for up");
      }
    if (migr_table[(*migr_table_counter) - 1].from != down->actualpop)
    {
        error ("this should never happen -> wrong choice of nodes\n");
        (*migr_table_counter)--;
    }
    if (((*migr_table_counter) > 0)
            && (migr_table[(*migr_table_counter) - 1].from != down->actualpop))
        error ("problem catched in inser_migr_table");
    theNode=up;
    for (i = 0; i < (*migr_table_counter); i++)
    {
      theNode = add_migration(world, theNode,
			      migr_table[i].from, 
			      migr_table[i].to,
			      migr_table[i].time - theNode->tyme);
    }
    if(down->tyme < migr_table[i-1].time)
      {
	error("Problem with migration addition in coalesce1p()");
      }
    //printf("%li migration nodelets added\n",(*migr_table_counter) * 2);
    down->back = theNode;
    theNode->back = down;
      /* the stuff below will be excised once I am complete sure that
	 this recycler above really works -- it seems currently dec 05 07
        tmp = (node *) mycalloc (1, sizeof (node));
        tmp2 = (node *) mycalloc (1, sizeof (node));
        oldNode = up;
        theNode = up->back;
        panic = 0;
	//printf("+");
        while (theNode->tyme < migr_table[i].time ****&& panic++ < 10000****)
        {
            if (theNode->type == 't' )
            {
                oldNode = theNode;
                theNode = theNode->back;
            }
            else
            {
                oldNode = theNode->next;
                theNode = theNode->next->back;
            }
        }
        tmp->back = oldNode;
        oldNode->back = tmp;
        tmp->number = -999;
        tmp->nayme = NULL;
        tmp->tip = FALSE;
        tmp->top = FALSE;
        tmp->dirty = TRUE;
        tmp->id = world->unique_id++;
        tmp->tyme = migr_table[i].time;
        tmp->type = 'm';
        tmp->actualpop = migr_table[i].to;
        tmp->pop = migr_table[i].from;
        tmp2->tyme = migr_table[i].time;
        tmp2->type = 'm';
        tmp2->id = world->unique_id++;
        tmp2->actualpop = migr_table[i].to;
        tmp2->pop = migr_table[i].from;
        tmp->next = tmp2;
        tmp2->next = tmp;
        tmp2->top = 1;
        oldNode2 = down;
        theNode = down->back;
        while (theNode->tyme > migr_table[i].time)
        {
            oldNode2 = theNode->next;
            theNode = theNode->next->back;
        }
        tmp2->back = oldNode2;
        oldNode2->back = tmp2;
}*/
}


void
children (node * mother, node ** brother, node ** sister)
{
    node *m;

    m = showtop (mother);

    if (m->type == 't')
    {
        error ("this is a tip, so there are no more child nodes\n");
    }
    else
    {
        (*brother) = crawlback (m->next);
        (*sister) = crawlback (m->next->next);
    }
}

/*       Uses Lanczos-type approximation to ln(gamma) for z > 0. */
/*       Reference: */
/*            Lanczos, C. 'A precision approximation of the gamma */
/*                    function', J. SIAM Numer. Anal., B, 1, 86-96, 1964. */
/*       Accuracy: About 14 significant digits except for small regions */
/*                 in the vicinity of 1 and 2. */
/*       Programmer: Alan Miller */
/*                   CSIRO Division of Mathematics & Statistics */
/*       Latest revision - 17 April 1988 */
/* translated and modified into C by Peter Beerli 1997 */
MYREAL
mylgamma (MYREAL z)
{
    MYREAL a[9] = { 0.9999999999995183, 676.5203681218835,
                    -1259.139216722289, 771.3234287757674, -176.6150291498386,
                    12.50734324009056, -0.1385710331296526, 9.934937113930748e-6,
                    1.659470187408462e-7
                  };
    MYREAL lnsqrt2pi = 0.9189385332046727;
    MYREAL result;
    long j;
    MYREAL tmp;
    if (z <= 0.)
    {
        return MYREAL_MAX;  /*this will kill the receiving calculation */
    }
    result = 0.;
    tmp = z + 7.;
    for (j = 9; j >= 2; --j)
    {
        result += a[j - 1] / tmp;
        tmp -= 1.;
    }
    result += a[0];
    result = LOG (result) + lnsqrt2pi - (z + 6.5) + (z - 0.5) * LOG (z + 6.5);
    return result;
}    /* lgamma */

/* ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3
   Computation of the Incomplete Gamma Integral 
   Auxiliary functions required: lgamma() = logarithm of the gamma 
   function, and alnorm() = algorithm AS66 
   in Mathematica this is GammaRegularized[a,0,x] === Gamma[a,0,x]/Gamma[a]
 */
MYREAL
incompletegamma (MYREAL tx, MYREAL talpha)
{
  const  MYREAL eps = DBL_EPSILON ;
  double gama, d_1, d_2, d_3;
  /*static */
  double  a, b, c, an, rn;
  /*static */
  double  pn1, pn2, pn3, pn4, pn5, pn6, arg;
  double x = (double) tx;
  double alpha = (double) talpha;

  gama = 0.;
  /*  Check that we have valid values for X and P */
  if (alpha <= 0. || x < 0.)
    error ("failed in imcompletegamma(): wrong alpha or x\n");
  if (fabs (x) < eps)
        return (MYREAL) gama;

    /*  Use a normal approximation if P > PLIMIT */
    if (alpha > 1e3)
    {
        pn1 =
            sqrt (alpha) * 3. * (pow (x / alpha, (1. / 3.)) + 1. / (alpha * 9.) -
                                 1.);
        gama = alnorm (pn1, FALSE);
        return (MYREAL) gama;
    }

   /*  If X is extremely large compared to P then set GAMMAD = 1 */
    if (x > 1e8)
    {
        gama = 1.;
        return (MYREAL) gama;
    }

    if (x <= 1. || x < alpha)
    {
        /*  Use Pearson's series expansion. */
        /*  (Note that P is not large enough to force overflow in lgamma()). */
        arg = alpha * LOG (x) - x - LGAMMA (alpha + 1.);
        c = 1.;
        gama = 1.;
        a = alpha;
        while (c > 1e-14)
        {
            a += 1.;
            c = c * x / a;
            gama += c;
        }
        arg += LOG (gama);
        gama = 0.;
        if (arg >= -88.)
        {
            gama = EXP (arg);
        }

    }
    else
    {
        /*  Use a continued fraction expansion */
        arg = alpha * LOG (x) - x - LGAMMA (alpha);
        a = 1. - alpha;
        b = a + x + 1.;
        c = 0.;
        pn1 = 1.;
        pn2 = x;
        pn3 = x + 1.;
        pn4 = x * b;
        gama = pn3 / pn4;
        for (;;)
        {
            a += 1.;
            b += 2.;
            c += 1.;
            an = a * c;
            pn5 = b * pn3 - an * pn1;
            pn6 = b * pn4 - an * pn2;
            if (fabs (pn6) > 0.)
            {
                rn = pn5 / pn6;
                /* Computing MIN */
                d_2 = 1e-14, d_3 = rn * 1e-14;
                if ((d_1 = gama - rn, fabs (d_1)) <= MIN (d_2, d_3))
                {
                    arg += LOG (gama);
                    gama = 1.;
                    if (arg >= -88.)
                    {
                        gama = 1. - EXP (arg);
                    }
                    return (MYREAL) gama;
                }
                gama = rn;
            }
            pn1 = pn3;
            pn2 = pn4;
            pn3 = pn5;
            pn4 = pn6;
            if (fabs (pn5) >= 1e37)
            {
                /*  Re-scale terms in continued fraction if terms are large */
                pn1 /= 1e37;
                pn2 /= 1e37;
                pn3 /= 1e37;
                pn4 /= 1e37;
            }
        }
    }
    return (MYREAL) gama;
}    /* incompletegamma() */


/* calculation is replaced by the correct function in
   polygamma.c (which is a translation of a fortran program by amos
 
   driver for the polygamma calculation */
MYREAL
polygamma (long n, MYREAL z)
{
    MYREAL ans;
    long nz, ierr;
    dpsifn (&z, &n, 1, 1, &ans, &nz, &ierr);
    if (n == 0)
        return -ans;
    else
        return ans;
}

/*-------------------------------------------------------*/
/* nrcheck subroutine (used in damped newton raphson proc */
/* syntax: nrcheck(matrix,inversematrix,ncols=nrows,returnval1,returnval2) */
/* mai 95 PB                                             */
boolean
nrcheck (MYREAL **m, MYREAL **tm, MYREAL *v, long nrows, MYREAL *r1,
         MYREAL *r2, boolean do_newton)
{
    long i, j, k;
    MYREAL *tmp, *tmp2, tmp3 = 0.0, tmp4 = 0.0;
    tmp = (MYREAL *) mycalloc (1, sizeof (MYREAL) * nrows);
    tmp2 = (MYREAL *) mycalloc (1, sizeof (MYREAL) * nrows);
    /*first evaluate r1 */
    (*r1) = (*r2) = 0.0;
    for (i = 0; i < nrows; i++)
    {
        (*r1) += v[i] * v[i];
    }
    /*                                       T    */
    for (j = 0; j < nrows; j++)
    {    /* g . G */
        for (k = 0; k < nrows; k++)
        {
            tmp[j] += v[k] * m[j][k];
            tmp2[j] += v[k] * tm[j][k];
        }
    }
    /*                                       T        */
    for (i = 0; i < nrows; i++)
    {    /* g . G . g */
        (*r2) += tmp[i] * v[i];
        tmp3 += tmp2[i] * v[i];
    }
    tmp4 = LOG (fabs ((*r1)));
    tmp4 = tmp4 + tmp4 - LOG (fabs ((*r2)));
    tmp4 = ((*r2) < 0 ? -1 : 1) * EXP (tmp4);
    myfree(tmp);
    if (do_newton && (tmp3 > (tmp4 > 0 ? tmp4 : 0)))
    {
        memcpy (v, tmp2, sizeof (MYREAL) * nrows);
        myfree(tmp2);
        return TRUE;
    }
    myfree(tmp2);
    return FALSE;
}


/*-------------------------------------------------------*/
/* Matrix inversion subroutine                           */
/* The passed matrix will be replaced by its inverse!!!!! */
/* Gauss-Jordan reduction -- invert matrix a in place,   */
/* overwriting previous contents of a.  On exit, matrix a */
/* contains the inverse.                                 */
void
invert_matrix (MYREAL **a, long nsize)
{
    long i, j;
    long *indeks;
    MYREAL *column, **result;
    indeks = (long *) mymalloc (sizeof (long) * nsize);
    column = (MYREAL *) mymalloc (sizeof (MYREAL) * nsize);
    result = (MYREAL **) mymalloc (sizeof (MYREAL *) * nsize);
    for (i = 0; i < nsize; i++)
    {
        result[i] = (MYREAL *) mymalloc (sizeof (MYREAL) * nsize);
    }
    lu_decomp (a, indeks, nsize);
    for (j = 0; j < nsize; j++)
    {
        memset (column, 0, sizeof (MYREAL) * nsize);
        column[j] = 1.0;
        lu_substitution (a, indeks, column, nsize);
        for (i = 0; i < nsize; i++)
            result[i][j] = column[i];
    }
    for (i = 0; i < nsize; i++)
    {
        memcpy (a[i], result[i], sizeof (MYREAL) * nsize);
        myfree(result[i]);
    }
    myfree(result);
    myfree(column);
    myfree(indeks);
}

/*=======================================================*/

/*-------------------------------------------------------*/
/* LU decomposition                                      */
/* after Dahlquist et al. 1974 and Press et al. 1988     */
/* the method's uses Crout's procedure and the pivoting  */
/* described in Press et al.                             */
/* Syntax: lu_decomp(matrix, indeks, nrows)               */
/* matrix will be destroyed and filled with the two      */
/* triangular matrices, indeks is the index vector for the */
/* pivoting and the row change in case of 0 pivot values */
/* nrows is the number of rows and columns in matrix     */
/* april 95 PB                                           */
void
lu_decomp (MYREAL **m, long *indeks, long nrows)
{
    long i, j, k, p, kmax = -1;
    MYREAL *max_row_vals, big, summ, pivot, bigt;
    max_row_vals = (MYREAL *) mycalloc (1, sizeof (MYREAL) * nrows);
    for (i = 0; i < nrows; i++)
    {
        big = 0.0;
        for (j = 0; j < nrows; j++)
        {
            if ((bigt = fabs (m[i][j])) > big)
                big = bigt;
        }
        max_row_vals[i] = 1.0 / big;
        if (big == 0.0)
        {
            error ("Singular matrix detected in lu_decomp\n");
        }
    }
    for (i = 0; i < nrows; i++)
    {
        for (k = 0; k < i; k++)
        {   /* upper half of matrix */
            summ = m[k][i];
            for (p = 0; p < k; p++)
                summ -= m[k][p] * m[p][i];
            m[k][i] = summ;
        }
        big = 0.0;
        for (k = i; k < nrows; k++)
        {   /* lower half of matrix */
            summ = m[k][i];
            for (p = 0; p < i; p++)
                summ -= m[k][p] * m[p][i];
            m[k][i] = summ;
            pivot = fabs (summ) /**max_row_vals[k]*/ ;
            /*  printf(stdout,"i=%li,pivot=%f,big=%f\n",i,pivot,big); */
            if (pivot >= big)
            {
                big = pivot;
                kmax = k;
            }
        }
        if (i != kmax)
        {
            for (p = 0; p < nrows; p++)
            {
                pivot = m[kmax][p];
                m[kmax][p] = m[i][p];
                m[i][p] = pivot;
            }
            max_row_vals[kmax] = max_row_vals[i];
        }
        indeks[i] = kmax;
        if (m[i][i] == 0.0)
            m[i][i] = SMALL_VALUE;
        if (i != nrows - 1)
        {
            pivot = 1. / m[i][i];
            for (k = i + 1; k < nrows; k++)
                m[k][i] *= pivot;
        }
    }
    myfree(max_row_vals);
}    /* end of lu_decomp */

/*-------------------------------------------------------*/
/* LU substitution                                       */
/* after Dahlquist et al. 1974 and Press et al. 1988     */
/* needs first the evaluation LU decomposition           */
/* Syntax: lu_substition(matrix, indeks, vector, nrows)   */
/* matrix = LU decomposed matrix, indeks = order of matrix */
/* vector = value vector, nrows = number of rows/columns */
/* april 95 PB                                           */
void
lu_substitution (MYREAL **m, long *indeks, MYREAL *v, long nrows)
{
    long i, j;
    MYREAL summ;
    for (i = 0; i < nrows; i++)
    {
        summ = v[indeks[i]];
        v[indeks[i]] = v[i];
        for (j = 0; j < i; j++)
            summ -= m[i][j] * v[j];
        v[i] = summ;
    }
    for (i = nrows - 1; i >= 0; i--)
    {
        summ = v[i];
        for (j = i + 1; j < nrows; j++)
            summ -= m[i][j] * v[j];
        v[i] = summ / m[i][i];
    }
}


/* Algorithm AS66 Applied Statistics (1973) vol22 no.3
   Evaluates the tail area of the standardised normal curve
   from x to infinity if upper is .true. or
   from minus infinity to x if upper is .false. */
MYREAL
alnorm (MYREAL x, int up)
{
    /* Initialized data */
    /* *** machine dependent constants ????????????? */
    /*static */ MYREAL zero = 0.;
    /*static */
    MYREAL a1 = 5.75885480458;
    /*static */
    MYREAL a2 = 2.62433121679;
    /*static */
    MYREAL a3 = 5.92885724438;
    /*static */
    MYREAL b1 = -29.8213557807;
    /*static */
    MYREAL b2 = 48.6959930692;
    /*static */
    MYREAL c1 = -3.8052e-8;
    /*static */
    MYREAL c2 = 3.98064794e-4;
    /*static */
    MYREAL c3 = -.151679116635;
    /*static */
    MYREAL c4 = 4.8385912808;
    /*static */
    MYREAL c5 = .742380924027;
    /*static */
    MYREAL one = 1.;
    /*static */
    MYREAL c6 = 3.99019417011;
    /*static */
    MYREAL d1 = 1.00000615302;
    /*static */
    MYREAL d2 = 1.98615381364;
    /*static */
    MYREAL d3 = 5.29330324926;
    /*static */
    MYREAL d4 = -15.1508972451;
    /*static */
    MYREAL d5 = 30.789933034;
    /*static */
    MYREAL half = .5;
    /*static */
    MYREAL ltone = 7.;
    /*static */
    MYREAL utzero = 18.66;
    /*static */
    MYREAL con = 1.28;
    /*static */
    MYREAL p = .398942280444;
    /*static */
    MYREAL q = .39990348504;
    /*static */
    MYREAL r = .398942280385;

    /*static */
    MYREAL y, result;

    if (x < zero)
    {
        up = !up;
        x = -x;
    }
    if (x <= ltone || (up && x <= utzero))
    {
        y = half * x * x;
        if (x > con)
        {
            result =
                r * EXP (-y) / (x + c1 +
                                d1 / (x + c2 +
                                      d2 / (x + c3 +
                                            d3 / (x + c4 +
                                                  d4 / (x + c5 +
                                                        d5 / (x + c6))))));
            return ((!up) ? one - result : result);
        }
        result =
            half - x * (p - q * y / (y + a1 + b1 / (y + a2 + b2 / (y + a3))));
        return ((!up) ? one - result : result);
    }
    else
    {
        return ((!up) ? 1.0 : 0.);
    }
    /*fake */ return -99;
}    /* alnorm */

/* dpsifn.c -- translated by f2c (version 19950808).
   and hand-patched by Peter Beerli Seattle, 1996
   SUBROUTINE DPSIFN (X, N, KODE, M, ANS, NZ, IERR)
 
   C***BEGIN PROLOGUE  DPSIFN
   C***PURPOSE  Compute derivatives of the Psi function.
   C***LIBRARY   SLATEC
   C***CATEGORY  C7C
   C***TYPE      MYREAL PRECISION (PSIFN-S, DPSIFN-D)
   C***KEYWORDS  DERIVATIVES OF THE GAMMA FUNCTION, POLYGAMMA FUNCTION,
   C             PSI FUNCTION
   C***AUTHOR  Amos, D. E., (SNLA)
   C***DESCRIPTION
   C
   C         The following definitions are used in DPSIFN:
   C
   C      Definition 1
   C         PSI(X) = d/dx (ln(GAMMA(X)), the first derivative of
   C                  the log GAMMA function.
   C      Definition 2
   C                     K   K
   C         PSI(K,X) = d /dx (PSI(X)), the K-th derivative of PSI(X).
   C   ___________________________________________________________________
   C      DPSIFN computes a sequence of SCALED derivatives of
   C      the PSI function; i.e. for fixed X and M it computes
   C      the M-member sequence
   C
   C                    ((-1)**(K+1)/GAMMA(K+1))*PSI(K,X)
   C                       for K = N,...,N+M-1
   C
   C      where PSI(K,X) is as defined above.   For KODE=1, DPSIFN returns
   C      the scaled derivatives as described.  KODE=2 is operative only
   C      when K=0 and in that case DPSIFN returns -PSI(X) + LN(X).  That
   C      is, the logarithmic behavior for large X is removed when KODE=2
   C      and K=0.  When sums or differences of PSI functions are computed
   C      the logarithmic terms can be combined analytically and computed
   C      separately to help retain significant digits.
   C
   C         Note that CALL DPSIFN(X,0,1,1,ANS) results in
   C                   ANS = -PSI(X)
   C
   C     Input      X is MYREAL PRECISION
   C           X      - Argument, X .gt. 0.0D0
   C           N      - First member of the sequence, 0 .le. N .le. 100
   C                    N=0 gives ANS(1) = -PSI(X)       for KODE=1
   C                                       -PSI(X)+LN(X) for KODE=2
   C           KODE   - Selection parameter
   C                    KODE=1 returns scaled derivatives of the PSI
   C                    function.
   C                    KODE=2 returns scaled derivatives of the PSI
   C                    function EXCEPT when N=0. In this case,
   C                    ANS(1) = -PSI(X) + LN(X) is returned.
   C           M      - Number of members of the sequence, M.ge.1
   C
   C    Output     ANS is MYREAL PRECISION
   C           ANS    - A vector of length at least M whose first M
   C                    components contain the sequence of derivatives
   C                    scaled according to KODE.
   C           NZ     - Underflow flag
   C                    NZ.eq.0, A normal return
   C                    NZ.ne.0, Underflow, last NZ components of ANS are
   C                             set to zero, ANS(M-K+1)=0.0, K=1,...,NZ
   C           IERR   - Error flag
   C                    IERR=0, A normal return, computation completed
   C                    IERR=1, Input error,     no computation
   C                    IERR=2, Overflow,        X too small or N+M-1 too
   C                            large or both
   C                    IERR=3, Error,           N too large. Dimensioned
   C                            array TRMR(NMAX) is not large enough for N
   C
   C         The nominal computational accuracy is the maximum of unit
   C         roundoff (=D1MACH(4)) and 1.0D-18 since critical constants
   C         are given to only 18 digits.
   C
   C         PSIFN is the single precision version of DPSIFN.
   C
   C *Long Description:
   C
   C         The basic method of evaluation is the asymptotic expansion
   C         for large X.ge.XMIN followed by backward recursion on a two
   C         term recursion relation
   C
   C                  W(X+1) + X**(-N-1) = W(X).
   C
   C         This is supplemented by a series
   C
   C                  SUM( (X+K)**(-N-1) , K=0,1,2,... )
   C
   C         which converges rapidly for large N. Both XMIN and the
   C         number of terms of the series are calculated from the unit
   C         roundoff of the machine environment.
   C
   C***REFERENCES  Handbook of Mathematical Functions, National Bureau
   C                 of Standards Applied Mathematics Series 55, edited
   C                 by M. Abramowitz and I. A. Stegun, equations 6.3.5,
   C                 6.3.18, 6.4.6, 6.4.9 and 6.4.10, pp.258-260, 1964.
   C               D. E. Amos, A portable Fortran subroutine for
   C                 derivatives of the Psi function, Algorithm 610, ACM
   C                 Transactions on Mathematical Software 9, 4 (1983),
   C                 pp. 494-502.
   C***ROUTINES CALLED  D1MACH, I1MACH
   C***REVISION HISTORY  (YYMMDD)
   C   820601  DATE WRITTEN
   C   890531  Changed all specific intrinsics to generic.  (WRB)
   C   890911  Removed unnecessary intrinsics.  (WRB)
   C   891006  Cosmetic changes to prologue.  (WRB)
   C   891006  REVISION DATE from Version 3.2
   C   891214  Prologue converted to Version 4.0 format.  (BAB)
   C   920501  Reformatted the REFERENCES section.  (WRB)
   C***END PROLOGUE  DPSIFN
 
 
 */

/*static*/ long fifteen = 15;
/*static*/
long sixteen = 16;
/*static*/
long five = 5;
/*static*/
long four = 4;
/*static*/
long fourteen = 14;

MYREAL
d1mach (long i)
{
  //#ifdef MYREAL == float
  //const  MYREAL eps = FLT_EPSILON ;
  //const  MYREAL numbermax = FLT_MAX;
  //const  MYREAL numbermin = FLT_MIN;
  //#else
  const  MYREAL eps = DBL_EPSILON ;
  const  MYREAL numbermax = MYREAL_MAX;
  const  MYREAL numbermin = DBL_MIN;
  //#endif

    switch (i)
    {
    case 1:
        return numbermin;
    case 2:
        return numbermax;
    case 3:
        return eps / FLT_RADIX;
    case 4:
        return eps;
    case 5:
        return log10 ((MYREAL) FLT_RADIX);
    }
    usererror ("invalid argument: d1mach(%ld)\n", i);
    return 0;   /* for compilers that complain of missing return values */
}

long
i1mach (long i)
{
    switch (i)
    {
    case 1:
        return 5;   /* standard input */
    case 2:
        return 6;   /* standard output */
    case 3:
        return 7;   /* standard punch */
    case 4:
        return 0;   /* standard error */
    case 5:
        return 32;  /* bits per integer */
    case 6:
        return 1;   /* Fortran 77 value */
    case 7:
        return 2;   /* base for integers */
    case 8:
        return 31;  /* digits of integer base */
    case 9:
        return LONG_MAX;
    case 10:
        return FLT_RADIX;
    case 11:
        return FLT_MANT_DIG;
    case 12:
        return FLT_MIN_EXP;
    case 13:
        return FLT_MAX_EXP;
    case 14:
        return DBL_MANT_DIG;
    case 15:
        return DBL_MIN_EXP;
    case 16:
        return DBL_MAX_EXP;
    }
    usererror ("invalid argument: i1mach(%ld)\n", i);
    return 0;   /* for compilers that complain of missing return values */
}

int
dpsifn (MYREAL *x, long *n, long kode, long m, MYREAL *ans, long *nz,
        long *ierr)
{
    /* Initialized data */

    /*static */ long nmax = 100;
    /*static */
    MYREAL b[22] = { 1., -.5, .166666666666666667,
                     -.0333333333333333333, .0238095238095238095, -.0333333333333333333,
                     .0757575757575757576, -.253113553113553114, 1.16666666666666667,
                     -7.09215686274509804, 54.9711779448621554, -529.124242424242424,
                     6192.1231884057971, -86580.2531135531136, 1425517.16666666667,
                     -27298231.067816092, 601580873.900642368, -15116315767.0921569,
                     429614643061.166667, -13711655205088.3328, 488332318973593.167,
                     -19296579341940068.1
                   };

    /* System generated locals */
    long i1, i2;
    MYREAL d1, d2;


    /* Local variables */
    /*static */
    MYREAL elim, xinc, xmin, tols, xdmy, yint, trmr[100], rxsq;
    /*static */
    long i__, j, k;
    /*static */
    MYREAL s, t, slope, xdmln, wdtol;
    /*static */
    MYREAL t1, t2;
    /*static */
    long fn;
    /*static */
    MYREAL ta;
    /*static */
    long mm, nn, np;
    /*static */
    MYREAL fx, tk;
    /*static */
    long mx, nx;
    /*static */
    MYREAL xm, tt, xq, den, arg, fln, r1m4, r1m5, eps, rln, tol,
    xln, trm[22], tss, tst;
    int i;
	for(i=0;i<22;i++)
        trm[i]=0.0;
    for(i=0;i<100;i++)
        trmr[i]=0.0;
    
    /* Parameter adjustments */
    --ans;

    /* Function Body */
    /* ----------------------------------------------------------------------- */
    /*             BERNOULLI NUMBERS */
    /* ----------------------------------------------------------------------- */

    /* ***FIRST EXECUTABLE STATEMENT  DPSIFN */
    *ierr = 0;
    *nz = 0;
    if (*x <= 0.)
    {
        *ierr = 1;
    }
    if (*n < 0)
    {
        *ierr = 1;
    }
    if (kode < 1 || kode > 2)
    {
        *ierr = 1;
    }
    if (m < 1)
    {
        *ierr = 1;
    }
    if (*ierr != 0)
    {
        return 0;
    }
    mm = m;
    /* Computing MIN */
    //xcode i1 = -fifteen;
    //xcode i2 = sixteen;
    nx = MIN (-i1mach (fifteen), i1mach (sixteen));
    r1m5 = d1mach (five);
    r1m4 = d1mach (four) * .5;
    wdtol = MAX (r1m4, 5e-19);
    /* ----------------------------------------------------------------------- */
    /*     ELIM = APPROXIMATE EXPONENTIAL OVER AND UNDERFLOW LIMIT */
    /* ----------------------------------------------------------------------- */
    elim = (nx * r1m5 - 3.) * 2.302;
    xln = LOG (*x);
L41:
    nn = *n + mm - 1;
    fn = nn;
    t = (fn + 1) * xln;
    /* ----------------------------------------------------------------------- */
    /*     OVERFLOW AND UNDERFLOW TEST FOR SMALL AND LARGE X */
    /* ----------------------------------------------------------------------- */
    if (fabs (t) > elim)
    {
        goto L290;
    }
    if (*x < wdtol)
    {
        goto L260;
    }
    /* ----------------------------------------------------------------------- */
    /*     COMPUTE XMIN AND THE NUMBER OF TERMS OF THE SERIES, FLN+1 */
    /* ----------------------------------------------------------------------- */
    rln = r1m5 * i1mach (fourteen);
    rln = MIN (rln, 18.06);
    fln = MAX (rln, 3.) - 3.;
    yint = fln * .4 + 3.5;
    slope = fln * (fln * 6.038e-4 + .008677) + .21;
    xm = yint + slope * fn;
    mx = (long) xm + 1;
    xmin = (MYREAL) mx;
    if (*n == 0)
    {
        goto L50;
    }
    xm = rln * -2.302 - MIN (0., xln);
    arg = xm / *n;
    arg = MIN (0., arg);
    eps = EXP (arg);
    xm = 1. - eps;
    if (fabs (arg) < .001)
    {
        xm = -arg;
    }
    fln = *x * xm / eps;
    xm = xmin - *x;
    if (xm > 7. && fln < 15.)
    {
        goto L200;
    }
L50:
    xdmy = *x;
    xdmln = xln;
    xinc = 0.;
    if (*x >= xmin)
    {
        goto L60;
    }
    nx = (long) (*x);
    xinc = xmin - nx;
    xdmy = *x + xinc;
    xdmln = LOG (xdmy);
L60:
    /* ----------------------------------------------------------------------- */
    /*     GENERATE W(N+MM-1,X) BY THE ASYMPTOTIC EXPANSION */
    /* ----------------------------------------------------------------------- */
    t = fn * xdmln;
    t1 = xdmln + xdmln;
    t2 = t + xdmln;
    /* Computing MAX */
    d1 = fabs (t), d2 = fabs (t1), d1 = MAX (d1, d2), d2 = fabs (t2);
    tk = MAX (d1, d2);
    if (tk > elim)
    {
        goto L380;
    }
    tss = EXP (-t);
    tt = .5 / xdmy;
    t1 = tt;
    tst = wdtol * tt;
    if (nn != 0)
    {
        t1 = tt + 1. / fn;
    }
    rxsq = 1. / (xdmy * xdmy);
    ta = rxsq * .5;
    t = (fn + 1) * ta;
    s = t * b[2];
    if (fabs (s) < tst)
    {
        goto L80;
    }
    tk = 2.;
    for (k = 4; k <= 22; ++k)
    {
        t = t * ((tk + fn + 1) / (tk + 1.)) * ((tk + fn) / (tk + 2.)) * rxsq;
        trm[k - 1] = t * b[k - 1];
        if ((d1 = trm[k - 1], fabs (d1)) < tst)
        {
            goto L80;
        }
        s += trm[k - 1];
        tk += 2.;
        /* L70: */
    }
L80:
    s = (s + t1) * tss;
    if (xinc == 0.)
    {
        goto L100;
    }
    /* ----------------------------------------------------------------------- */
    /*     BACKWARD RECUR FROM XDMY TO X */
    /* ----------------------------------------------------------------------- */
    nx = (long) xinc;
    np = nn + 1;
    if (nx > nmax)
    {
        goto L390;
    }
    if (nn == 0)
    {
        goto L160;
    }
    xm = xinc - 1.;
    fx = *x + xm;
    /* ----------------------------------------------------------------------- */
    /*     THIS LOOP SHOULD NOT BE CHANGED. FX IS ACCURATE WHEN X IS SMALL */
    /* ----------------------------------------------------------------------- */
    i1 = nx;
    for (i__ = 1; i__ <= i1; ++i__)
    {
        i2 = -np;
        trmr[i__ - 1] = pow (fx, (MYREAL) i2);
        s += trmr[i__ - 1];
        xm += -1.;
        fx = *x + xm;
        /* L90: */
    }
L100:
    ans[mm] = s;
    if (fn == 0)
    {
        goto L180;
    }
    /* ----------------------------------------------------------------------- */
    /*     GENERATE LOWER DERIVATIVES, J.LT.N+MM-1 */
    /* ----------------------------------------------------------------------- */
    if (mm == 1)
    {
        return 0;
    }
    i1 = mm;
    for (j = 2; j <= i1; ++j)
    {
        --fn;
        tss *= xdmy;
        t1 = tt;
        if (fn != 0)
        {
            t1 = tt + 1. / fn;
        }
        t = (fn + 1) * ta;
        s = t * b[2];
        if (fabs (s) < tst)
        {
            goto L120;
        }
        tk = (MYREAL) (fn + 4);
        for (k = 4; k <= 22; ++k)
        {
            trm[k - 1] = trm[k - 1] * (fn + 1) / tk;
            if ((d1 = trm[k - 1], fabs (d1)) < tst)
            {
                goto L120;
            }
            s += trm[k - 1];
            tk += 2.;
            /* L110: */
        }
L120:
        s = (s + t1) * tss;
        if (xinc == 0.)
        {
            goto L140;
        }
        if (fn == 0)
        {
            goto L160;
        }
        xm = xinc - 1.;
        fx = *x + xm;
        i2 = nx;
        for (i__ = 1; i__ <= i2; ++i__)
        {
            trmr[i__ - 1] *= fx;
            s += trmr[i__ - 1];
            xm += -1.;
            fx = *x + xm;
            /* L130: */
        }
L140:
        mx = mm - j + 1;
        ans[mx] = s;
        if (fn == 0)
        {
            goto L180;
        }
        /* L150: */
    }
    return 0;
    /* ----------------------------------------------------------------------- */
    /*     RECURSION FOR N = 0 */
    /* ----------------------------------------------------------------------- */
L160:
    i1 = nx;
    for (i__ = 1; i__ <= i1; ++i__)
    {
        s += 1. / (*x + nx - i__);
        /* L170: */
    }
L180:
    if (kode == 2)
    {
        goto L190;
    }
    ans[1] = s - xdmln;
    return 0;
L190:
    if (xdmy == *x)
    {
        return 0;
    }
    xq = xdmy / *x;
    ans[1] = s - LOG (xq);
    return 0;
    /* ----------------------------------------------------------------------- */
    /*     COMPUTE BY SERIES (X+K)**(-(N+1)) , K=0,1,2,... */
    /* ----------------------------------------------------------------------- */
L200:
    nn = (long) fln + 1;
    np = *n + 1;
    t1 = (*n + 1) * xln;
    t = EXP (-t1);
    s = t;
    den = *x;
    i1 = nn;
    for (i__ = 1; i__ <= i1; ++i__)
    {
        den += 1.;
        i2 = -np;
        trm[i__ - 1] = pow (den, (MYREAL) i2);
        s += trm[i__ - 1];
        /* L210: */
    }
    ans[1] = s;
    if (*n != 0)
    {
        goto L220;
    }
    if (kode == 2)
    {
        ans[1] = s + xln;
    }
L220:
    if (mm == 1)
    {
        return 0;
    }
    /* ----------------------------------------------------------------------- */
    /*     GENERATE HIGHER DERIVATIVES, J.GT.N */
    /* ----------------------------------------------------------------------- */
    tol = wdtol / 5.;
    i1 = mm;
    for (j = 2; j <= i1; ++j)
    {
        t /= *x;
        s = t;
        tols = t * tol;
        den = *x;
        i2 = nn;
        for (i__ = 1; i__ <= i2; ++i__)
        {
            den += 1.;
            trm[i__ - 1] /= den;
            s += trm[i__ - 1];
            if (trm[i__ - 1] < tols)
            {
                goto L240;
            }
            /* L230: */
        }
L240:
        ans[j] = s;
        /* L250: */
    }
    return 0;
    /* ----------------------------------------------------------------------- */
    /*     SMALL X.LT.UNIT ROUND OFF */
    /* ----------------------------------------------------------------------- */
L260:
    i1 = -(*n) - 1;
    ans[1] = pow (*x, (MYREAL) i1);
    if (mm == 1)
    {
        goto L280;
    }
    k = 1;
    i1 = mm;
    for (i__ = 2; i__ <= i1; ++i__)
    {
        ans[k + 1] = ans[k] / *x;
        ++k;
        /* L270: */
    }
L280:
    if (*n != 0)
    {
        return 0;
    }
    if (kode == 2)
    {
        ans[1] += xln;
    }
    return 0;
L290:
    if (t > 0.)
    {
        goto L380;
    }
    *nz = 0;
    *ierr = 2;
    return 0;
L380:
    ++(*nz);
    ans[mm] = 0.;
    --mm;
    if (mm == 0)
    {
        return 0;
    }
    goto L41;
L390:
    *nz = 0;
    *ierr = 3;
    return 0;
}    /* dpsifn_ */




MYREAL
rannor (MYREAL mean, MYREAL sd)
{
    MYREAL r1, r2;
    r1 = RANDUM ();
    r2 = RANDUM ();
    return sd * sqrt (-2. * LOG (r1)) * cos (TWOPI * r2) + mean;
}


char
lowercase (int c)
{
    return (char) tolower (c);
}

char
uppercase (int c)
{
    return (char) toupper (c);
}

void upper(char *from, char **to)
{
    long i=0;
    while(from[i] != '\0')
    {
       (*to)[i] = from[i];
        ++i;
    }
    (*to)[++i] = '\0';
}

MYREAL
find_chi (long df, MYREAL prob)
{
    double a, b, m;
    double xb = 200.0;
    double xa = 0.0;
    double xm = 5.;
    double dprob = (double) prob;
    a = probchi (df, xa);
    m = probchi (df, xm);
    b = probchi (df, xb);
    while (fabs (m - dprob) > EPSILON)
    {
        if (m < dprob)
        {
            b = m;
            xb = xm;
        }
        else
        {
            a = m;
            xa = xm;
        }
        xm = (-(b * xa) + prob * xa + a * xb - dprob * xb) / (a - b); //(xa + xb)/2.;

        m = probchi (df, xm);
    }
    return (MYREAL) xm;
}


MYREAL
probchi (long df, double chi)
{
  const  MYREAL eps = DBL_EPSILON ;
    double prob;
    double v = ((MYREAL) df) / 2.;
    
    if (chi > eps && v > eps)
    {
        //lg = EXP (LGAMMA (v));
        prob = 1. - incompletegamma (chi / 2., v);
    }
    else
        prob = 1.0;
    //  printf("prob=%f v=%f chi=%f lg(v/2)=%f  ig(chi/2,v/2)=%f\n",
    //  prob,v,chi,lg, incompletegamma(chi/2.,v/2.));

    return prob;
}

MYREAL
probchiboundary (MYREAL chi, long zeros, long all)
{
    long nonzeros = all - zeros;

//xcode    MYREAL a, b, m;
    MYREAL m;
    MYREAL xb = 1.0;  //1.0-EPSILON/1000.;
    MYREAL xa = 0.0;  //EPSILON/1000.;
    MYREAL xm = 0.51;
    if(all==0)
        return 1.;
    if (zeros == 0)
    {
        return probchi (all, chi);
    }
    //xcode a = chiboundary (zeros, nonzeros, xa);
    m = chiboundary (zeros, nonzeros, xm);
    //xcode b = chiboundary (zeros, nonzeros, xb);
    while (fabs (m - chi) > EPSILON
            && (fabs (xa - xm) > EPSILON && fabs (xb - xm) > EPSILON))
    {
        if (m < chi)
        {
            //xcode b = m;
            xb = xm;
        }
        else
        {
            //xcode a = m;
            xa = xm;
        }
        xm =   /*(-(b * xa) + chi * xa + a * xb - chi * xb) / (a - b);      */
            (xa + xb) / 2.;

        m = chiboundary (zeros, nonzeros, xm);
    }
    return xm;
}

MYREAL
chiboundary (long zeros, long nonzeros, MYREAL alpha)
{
    //  MYREAL prob;
    MYREAL sum = 0.;
    long i;
    long k = zeros;
    MYREAL freq;
    MYREAL summ;
    //printf("z=%li a=%4.2f ", k,alpha);
    for (i = 0; i <= zeros; i++)
    {
        //      sum = zerovec[i]  * (i+nonzeros) == 0 ? 0. : find_chi(i+nonzeros,alpha);
        freq = EXP (logfac (k) - logfac (k - i) - logfac (i) - LOG2 * k);
        summ = (i + nonzeros) == 0 ? 0. : find_chi (i + nonzeros, alpha);
        //      printf(" %.2f(%.4f)",freq*pow(2.,k),summ);
        sum += freq * summ;
    }
    //printf("\n");
    return sum;
}


MYREAL
chisquare (long df, MYREAL alpha)
{
    const MYREAL table05[] =
        {
            3.84146, 5.99147, 7.81473, 9.48773, 11.0705, 12.5916
        };
    const MYREAL table01[] =
        {
            6.63490, 9.21034, 11.3449, 13.2767, 15.0863, 16.8119
        };

    if (alpha == 0.05)
        return table05[df - 1];
    if (alpha == 0.01)
        return table01[df - 1];
    error ("Chi-distribution for any probability alpha is not implemented");
    return -1;
}

MYREAL
calc_sum (MYREAL *vector, long n)
{
    long i;
    MYREAL summ = 0.0;
    for (i = 0; i < n; i++)
        summ += vector[i];
    return summ;
}

//==========================================
// searching and finding

boolean
find (long i, long *list, long listlen)
{
    long j;
    for (j = 0; j < listlen; j++)
    {
        if (i == list[j])
            return TRUE;
    }
    return FALSE;
}

//====================================================
// conversion between the different parameter schemes
// returns the begining of  mig_.i
long
mstart (long pop, long numpop)
{
    return numpop + pop * numpop - pop;
}

// returns the end of  mig_.i
long
mend (long pop, long numpop)
{
    return numpop + pop * numpop - pop + numpop - 1;
}

///
/// return the first element of population pop
long
mmstart (long pop, long numpop)
{
    return pop * (numpop);
}


///
/// return the element past the last element in of population pop
long
mmend (long pop, long numpop)
{
    return pop * numpop + numpop;
}

/// 
/// Returns the location in a full matrix given the abbreviated matrix
/// given i,j coordinates in an abbreviated matrix {d,d,...,d,a,a,a,...,a,b,b,...,b, ...}
/// it returns the position j* in a linearized full matrix with diagonal element "-"
/// the linearized matrix would look like this {-, a,a,a,...,a, b, - , b,  .... ,b, c, c, -, c, ....},
/// the position of i is the same for the abbreviated and the full matrix, and is not returned because
/// the calling function knows this already.
long
mm2m (long frompop, long topop, long numpop)
{
    if (frompop == topop)
        return (frompop);
    if (frompop < topop)
        return numpop + topop * (numpop - 1) + frompop;
    else
        return numpop + topop * (numpop - 1) + (frompop - 1);
}

///
/// Calulates the j and i from a linear abbreviated matrix
/// position z in {d,d,d,...,d,a,a,...,a,b,b,...,b, c, ...} we know how many populations c (columns or rows) 
/// exist and wnat to calculate the i,j coordinates in the c x c matrix {{d,a,a,...,a},{b,d,b,...},...}
/// in the function z above is im and i = frompop, j is topop
/// frompop and topop contain the result
void
m2mm (long i, long numpop, long *frompop, long *topop)
{
    if (i < numpop)
    {
        *frompop = i;
        *topop = i;
        return;
    }
    else
    {
        (*topop) = (long) (i - numpop) / (numpop - 1);
        (*frompop) = i - numpop - (*topop) * (numpop - 1);
        if (*frompop >= *topop)
            *frompop += 1;
    }
}


///
/// position i in linear abbreviated array
long
m2mml (long i, long numpop)
{
    long topop, frompop;

    if (i < numpop)
    {
        return i * numpop + i;
    }
    else
    {
        topop = (long) (i - numpop) / (numpop - 1);
        frompop = i - numpop - (topop) * (numpop - 1);
        if (frompop >= topop)
            frompop += 1;
        return numpop * topop + frompop;
    }
}


long
mml2m (long pos, long numpop)
{
    long topop = 0, frompop = 0, i = 1;
    if (pos == 0)
        return 0;
    while (pos > numpop * (i++))
        topop++;
    frompop = pos - topop * numpop;
    return mm2m (frompop, topop, numpop);
}


long
m2mml2 (long i, long topop, long numpop)
{
    long frompop;

    if (i < numpop)
    {
        return i * numpop + i;
    }
    else
    {
        frompop = i - numpop - (topop) * (numpop - 1);
        if (frompop >= topop)
            frompop += 1;
        return numpop * topop + frompop;
    }
}

void
gamma_rates (MYREAL *rate, MYREAL *probcat, long categs, char *input)
{
    long i;
    MYREAL alpha = MYREAL_MAX;
    MYREAL value;
    while (!isdigit (*input) && *input != '\0')
        input++;
    if ((value = strtod (input, (char **) NULL)) > 0)
        alpha = value;
    initgammacat (categs, alpha, 1., rate, probcat);
    for (i = 0; i < categs; i++)
    {
        probcat[i] = EXP (probcat[i]);
    }
    //  calc_gamma (alpha, rate, categs);
}

/* calculation of rate values following a gamma distribution for
   given probability values */
void
calc_gamma (MYREAL alpha, MYREAL *gama, long categs)
{
    long i, panic;
    MYREAL low, mid, high, xlow, xhigh, tmp, freq = 0, x = 10, elements =
                (MYREAL) categs;
    freq = -(0.5 / elements); /*so we have midpoints instead of endpoints */
    for (i = 0; i < categs; i++)
    {
        low = 0;
        mid = incompletegamma (10., alpha);
        high = 1.;
        freq += 1. / (elements);
        if (freq < mid)
        {
            high = mid;
            xlow = 0;
            xhigh = 10.;
            x = 5.;
        }
        else
        {
            low = mid;
            xhigh = 1e10;
            xlow = 10.;
            x = 1e5;
        }
        panic = 0;
        while (panic++ < 1000 && fabs (low - high) > 0.0001 && x > 0.000000001)
        {
            mid = incompletegamma (x, alpha);
            if (freq < mid)
            {
                high = mid;
                tmp = x;
                x = (x + xlow) / 2.;
                xhigh = tmp;
            }
            else
            {
                low = mid;
                tmp = x;
                x = (x + xhigh) / 2.;
                xlow = tmp;
            }
        }
        gama[i] = x / alpha;
        //Debug
        //      printf (stderr, "  %li> %f\n", i, gama[i]);

        if (x >= 10e10)
        {
            error ("calc_gamma(): x is too big");
        }
    }
}


void
fprintf2(FILE *file, long filesize, const char *fmt, ...)
{
    char *p;
    va_list ap;
    long bufsize;
    long pallocsize = filesize+strlen(fmt)+1;
    p  = (char *) mycalloc(pallocsize,sizeof(char));
    va_start(ap, fmt);
    bufsize = vsprintf(p, fmt, ap);
    if(bufsize>=pallocsize)
      error("failed in printf2()");
    fprintf(file,"%s", p);
    va_end(ap);
    myfree(p);
}


void
print_line (FILE * outfile, char c, long nn, long flag)
{
    long i, start = 0;
    switch (flag)
    {
    case START:
        start = 2;
        FPRINTF (outfile, "=--");
        break;
    case STOP:
        start = 2;
        FPRINTF (outfile, "==-");
        break;
    default:
        start = 0;
    }
    for (i = start; i < nn; i++)
    {
        FPRINTF (outfile, "%c",c);
    }
    FPRINTF(outfile, "\n");
}


void
sprint_line
(char *buffer, char c, long nn, long flag)
{
    char ch[2];
    long i, start = 0;
    char fp[LINESIZE];
    buffer[0] = '\0';
    ch[0] = c;
    ch[1] = '\0';
    switch (flag)
    {
    case START:
        start = 2;
        sprintf (fp, "=%c%c", c, c);
        strcat (buffer, fp);
        break;
    case STOP:
        start = 2;
        sprintf (fp, "==%c", c);
        strcat (buffer, fp);
        break;
    default:
        start = 0;
    }
    for (i = start; i < nn; i++)
    {
        strcat (buffer, ch);
    }
    strcat (buffer, "\n");
}

char
sgetc (char **buffer)
{
    char ch;
    ch = **buffer;
    (*buffer)++;
    return ch;
}

/// line-end transparent string gets command
char *
sgets (char *s, int size, char **stream)
{
    long ch = '\0';
    long counter = 0;
    while (counter < size - 1)
    {
        ch = **stream;
        (*stream)++;
        switch (ch)
        {
        case '\0':
        case '\r':
        case '\n':
            s[counter] = '\0';
            return s;
        default:
            s[counter] = (char) ch;
            break;
        }
        counter++;
    }
    return s;
}


/// sticks a text into the printing buffer
void add_to_buffer(char *fp, long *bufsize, char **buffer, long *allocbufsize)
{
  long fpsize = (long) strlen (fp) + 1;
  if(*allocbufsize <= (*bufsize + fpsize))
    {
      *allocbufsize += 100 * fpsize;
      (*buffer) = (char *) myrealloc (*buffer, (*allocbufsize) * sizeof (char));
    }
  (*bufsize) += sprintf((*buffer) + (*bufsize),"%s",fp);
}


/// sticks a text into the printing buffer
long print_to_buffer(char **buffer, long *maxbufsize, char *tempbuffer, long *pos, const char *fmt, ...)
{
  long mypos=0;
  char *p = tempbuffer;
  va_list ap;
  //p = (char *) mycalloc(1024,sizeof(char));
  va_start(ap, fmt);
  mypos = vsprintf(p, fmt, ap);
  va_end(ap);
  if((*pos + mypos) < (*maxbufsize))
      {
	(*pos) += sprintf((*buffer) + (*pos), "%s",p);
	//	if(*pos + mypos >= *maxbufsize )
	//  {
	//    printf("%i> pos=%li + mypos=%li >  maxbufsize=%li \n",myID, *pos, mypos, *maxbufsize); 
	//  }
      }
    else
      {
	*maxbufsize = *pos + 4 * mypos; // add some extra space
	(*buffer) = (char *) myrealloc ((*buffer), (*maxbufsize) * sizeof (char));
	(*pos) += sprintf((*buffer) + (*pos), "%s",p);
      }
  //  myfree(p);
  return (*pos);
}

/// sticks a text into the warning buffer that is printed at the end of the PDF and the end of the TEXT file
void record_warnings(world_fmt * world, const char *fmt, ...)
{
  long mypos=0;
  char *p = (char *) calloc(LINESIZE,sizeof(char));
  va_list ap;
  va_start(ap, fmt);
  mypos = vsprintf(p, fmt, ap);
  va_end(ap);
  if((world->warningsize + mypos) < world->warningallocsize)
      {
	world->warningsize += sprintf(world->warning + world->warningsize, "%s\n",p);
      }
    else
      {
	world->warningallocsize = world->warningsize + 4 * mypos; // add some extra space
	if(world->warning != NULL)
	    world->warning = (char *) myrealloc (world->warning, world->warningallocsize * sizeof (char));
	else
	  world->warning = (char *) mycalloc (world->warningallocsize, sizeof (char));
	world->warningsize += sprintf(world->warning + world->warningsize, "%s\n",p);
      }
  myfree(p);
}

void print_warning2(FILE *file, world_fmt *world)
{
      fprintf(file,"\nPOTENTIAL PROBLEMS\n");
      fprintf(file,"-------------------------------------------------------------------\n");
      fprintf(file,"This section reports potential problems with your run. Such reporting\n");      
      fprintf(file,"is tricky. When many parameters in  multilocus analysis\n");      
      fprintf(file,"are estimated then it is very common that some parameters for some loci\n");      
      fprintf(file,"will not be very informative, triggering suggestions (for example to\n");
      fprintf(file,"increase the prior range) that are not sensible. This suggestion tool\n");
      fprintf(file,"will improve with time. If some parameters are flagged inspect the tables\n");
      fprintf(file,"carefully and judge wether an action is required. For example if you run\n");
      fprintf(file,"a Bayesian inference with sequence data, for macroscopic species there is\n");
      fprintf(file,"rarely the need to increase the prior for Theta beyond 0.1, but if you use\n");
      fprintf(file,"microsatellites it is rather common that your prior for Theta has a range\n");
      fprintf(file,"from 0.0 to 100 or more. With many populations (>3) it is also very common\n");      
      fprintf(file,"that some migration routes are estimated poorly because there is no data.\n");
      fprintf(file,"Increasing the range will not help in such situations, reducing number of\n");
      fprintf(file,"parameters may help in such situations.\n\n");
      fprintf(file,"-------------------------------------------------------------------\n");
      if (world->warning[0]=='\0')
	fprintf(file,"No warning was recorded during the run\n\n");
      else
	fprintf(file,"%s\n\n", world->warning);
      fprintf(file,"-------------------------------------------------------------------\n");
      
}



void print_stored_warnings(world_fmt *world)
{
  if(world->warningsize>0)
    {  
      print_warning2(world->outfile,world);
      if(world->options->progress)
	{
	  print_warning2(stdout,world);
	}
    }
}


// calculates a rough approximation to log(val)
// max error = 0.0007
// Submitted by Laurent de Soras, posted on 30 March 2001
// http://www.flipcode.com/cgi-bin/msg.cgi?showThread=Tip-Fastlogfunction&forum=totd&id=-1
// Fast log() Function, by Laurent de Soras:
// Here is a code snippet to replace the slow log() function...
// It just performs an approximation, but the maximum error is below 0.007.
// Speed gain is about x5, and probably could be increased by tweaking the assembly code.
// The function is based on floating point coding.
// It's easy to get floor (log2(N)) by isolating exponent part.
// We can refine the approximation by using the mantissa. This function returns log2(N)
/*
MYINLINE float fast_log2 (float vval)
{
    float vv;
    int * const exp_ptr =  (int *) (&vval);
    int            x = *exp_ptr;
    const int      log_2 = ((x >> 23) & 255) - 128;
    x &= ~(255 << 23);
    x += 127 << 23;
    *exp_ptr = x;

    vval = ((-1.0/3) * vval + 2) * vval - 2.0/3;   // (1)
    vv = vval + log_2;
    return (vv);
}

// The line (1) computes 1+log2(m), m ranging from 1 to 2. The proposed
// formula is a 3rd degree polynomial keeping first derivate
// continuity. Higher degree could be used for more accuracy. For faster
// results, one can remove this line, if accuracy is not the matter (it
// gives some linear interpolation between powers of 2).
//Now we got log2(N), we have to multiply it by ln(2) to get the natural log :


MYINLINE float fast_log (float vval)
{
    float v = fast_log2 (vval) * 0.69314718f;
    //fprintf(stdout,"val= %f fast_log=%f log=%f\n", (float) vval, (float) v, log(vval));
    return v;
}

//#define myEXPA (1048576 / 0.693147180559945309417232121458)
#define myEXPA 1512775.39519518569383584038231 
#define myEXPC 60801 

MYINLINE double fast_exp(double y) 
{ 
    union 
    { 
        double d; 
      // weird the implementation says ifdef but that does no work at all
#ifndef LITTLE_ENDIAN 
        struct { int j, i; } n; 
#else 
        struct { int i, j; } n; 
#endif 
    } eco;
    eco.n.i = (int)(myEXPA*(y)) + (1072693248 - myEXPC);
    eco.n.j = 0; 
    return eco.d; 
} 

*/
