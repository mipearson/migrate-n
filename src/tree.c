/* \file tree.c */
/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
	T R E E B U I L D I N G   R O U T I N E S

	Peter Beerli 1996, Seattle
	beerli@fsu.edu

	Copyright 1997-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
	Copyright 2003-2004 Peter Beerli, Tallahassee FL

	some code in this file are successors of code in dnaml in the PHYLIP
        package of Joseph Felsenstein. Several changes were made to the conditionl
        likelihood methods to improve speed, but the original design idea is Joe's.

	A new model that includes treatment of gaps in the alignment is in the works
        and may superceede all other models because it seems that this model
        is capable to handle all other sequence models. <THIS IS NOY WORKING YET>

        This software is distributed free of charge for non-commercial use
	and is copyrighted. Of course, we do not guarantee that the software
	works and are not responsible for any damage you may cause or have.


$Id: tree.c 1833 2011-03-20 19:09:41Z beerli $

-------------------------------------------------------*/

#include "migration.h"
#include "sighandler.h"
#include "random.h"
#include "options.h"
#include "data.h"
#include "sequence.h"
#include "world.h"
#include "tools.h"
#include "migrate_mpi.h"
#include "mcmc.h"
#ifdef UEP
#include "uep.h"
#endif
#include "tree.h"
#ifdef BEAGLE
#include "calculator.h"
#endif
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

#define NOTIPS 0
#define WITHTIPS 1

extern long unique_id_global;

//// prototypes -------------------------------------------
////
////
//// creating trees
//void buildtree (world_fmt * world, option_fmt * options, data_fmt * data, long locus);
//void create_treetimelist (world_fmt * world,  timelist_fmt ** ltl, long locus);
//void fix_times (world_fmt * world, option_fmt * options);
//void first_smooth (world_fmt * world, long locus);
//void set_dirty (node * p);
//void construct_tymelist (world_fmt * world, timelist_fmt * timevector);
//void timeslices (timelist_fmt ** timevector);
//void add_partlineages (long numpop, timelist_fmt ** timevector);
//// likelihood calculation
//MYREAL treelikelihood (world_fmt * world);
//MYREAL pseudotreelikelihood (world_fmt * world, proposal_fmt * proposal);
MYREAL treelike_anc (world_fmt * world, long locus);
//void set_pop (node * theNode, long pop, long actualpop);
//void pseudonuview (proposal_fmt * proposal, xarray_fmt xx1, MYREAL *lx1,
//                   MYREAL v1, xarray_fmt xx2, MYREAL *lx2, MYREAL v2);
//void ltov (node * p);
//void treeout (FILE * treefile, node * joint, node * p, long s);
//void print_tree (world_fmt * world, long g, long *filepos);
//MYREAL find_tipdate(char * id, long pop, world_fmt *world);
//
/* private functions------------------------------------- */

void allocate_tree (world_fmt * world, option_fmt * options, data_fmt * data,
                    long locus);
/* allocations of nodes */
void allocatetips (world_fmt * world, option_fmt * options, data_fmt * data,
                   long locus);
void allocateinterior (world_fmt * world, data_fmt * data, long locus);
void allocatepoproot (world_fmt * world, data_fmt * data, long locus);
void allocate_tip (world_fmt * world, option_fmt * options, node ** p,
                   long pop, long locus, long a, long ind, char **tipnames);
void alloc_seqx (world_fmt * world, node * theNode);
void     allocate_xseq(xarray_fmt *x, long sites, long categs);

/* first tree material (upgma, distance) */
void set_tree (world_fmt * world, option_fmt * options, data_fmt * data,
               long locus);
// sets migration events into a tree read from the user (dna only)
void set_migrations (world_fmt * world, option_fmt * options, data_fmt * data,
					 long locus);
void distance_EP (char **data, long tips, MYREAL **m);
void distance_micro (char **data, long tips, MYREAL **m);
void distance_sequence (data_fmt * data, long locus, long tips, long sites,
                        long nmlength, MYREAL **m);
void distance_allele (world_fmt * world, option_fmt * options, long locus,
                      long tips, MYREAL **distm);
void randomize_distm(MYREAL **distm, data_fmt * data, long locus,
		     long tips, char *custm);
void constrain_distance_zeromig (MYREAL **m, option_fmt *options, data_fmt * data, long locus,
                                 long tips, char *custm);

void makevalues (world_fmt * world, option_fmt * options, data_fmt * data,
                 long locus);
void make_alleles (world_fmt * world, option_fmt * options, data_fmt * data, long locus);
void make_microsatellites (world_fmt * world, option_fmt * options, data_fmt * data, long locus);
void make_microbrownian (world_fmt * world, option_fmt *options, data_fmt * data, long locus);
void upgma (world_fmt * world, MYREAL **x, long tips, node ** nodep);
void set_top (world_fmt * world, node * p, long locus);
void set_v (node * p);
void calc_sancost (MYREAL **cost, long numpop);


void free_treetimes (world_fmt * world, long size);
void traverseNodes (node * theNode, timelist_fmt ** timevector, long *slice, world_fmt *world, long *tips);
void increase_timelist (timelist_fmt ** timevector);
void allocate_lineages (timelist_fmt **timevector, const long offset, const long numpop);
void smooth (const node * root, node * p, world_fmt * world,
             const long locus);
void which_nuview (char datatype, boolean fastlike, boolean use_gaps, int watkins);
void nuview_allele (node * mother, world_fmt * world, const long locus);
void nuview_micro (node * mother, world_fmt * world, const long locus);
void nuview_brownian (node * mother, world_fmt * world, const long locus);
void nuview_sequence (node * mother, world_fmt * world, const long locus);
void nuview_sequence_slow (node * mother, world_fmt * world,
                           const long locus);
void nuview_ancestral (node * mother, world_fmt * world, const long locus);
void adjustroot (node * r);
MYINLINE MYREAL pseudo_tl_seq (phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
                               proposal_fmt * proposal, world_fmt * world);
MYINLINE  MYREAL pseudo_tl_snp (phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
                                proposal_fmt * proposal, world_fmt * world);
MYINLINE  MYREAL pseudo_tl_snp_unlinked (phenotype xx1, phenotype xx2, MYREAL v1,
                                         MYREAL v2, proposal_fmt * proposal,
                                         world_fmt * world);
MYINLINE MYREAL pseudo_tl_anc (phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
                               proposal_fmt * proposal, world_fmt * world);
MYINLINE void pseudonu_allele (proposal_fmt * proposal, MYREAL **xx1, MYREAL *lx1,
                               MYREAL vv1, MYREAL *xx2, MYREAL lx2, MYREAL vv2);
MYINLINE void pseudonu_micro (proposal_fmt * proposal, MYREAL **xx1, MYREAL *lx1,
                              MYREAL v1, MYREAL *xx2, MYREAL lx2, MYREAL v2);
MYINLINE void pseudonu_brownian (proposal_fmt * proposal, MYREAL **xx1, MYREAL *lx1,
                                 MYREAL v1, MYREAL *xx2, MYREAL lx2, MYREAL v2);
MYINLINE void pseudonu_seq (proposal_fmt * proposal, phenotype xxx1, MYREAL v1,
                            phenotype xxx2, MYREAL v2);
MYINLINE void pseudonu_seq_slow (proposal_fmt * proposal, phenotype xxx1, MYREAL *sxx1,
                                 MYREAL v1, phenotype xxx2, MYREAL *sxx2, MYREAL v2);
MYINLINE void pseudonu_anc (proposal_fmt * proposal, phenotype xxx1, MYREAL v1,
                            phenotype xxx2, MYREAL v2);
void calculate_steps (world_fmt * world);
MYREAL logfac (long n);


boolean treereader (world_fmt * world,  option_fmt *options,  data_fmt * data);
void length_to_times (node * p);
boolean treeread (FILE * file, world_fmt * world, option_fmt *options, node ** pp, node * q);
char processlength (FILE * file, node ** p, option_fmt * options);
node *allocate_nodelet (world_fmt *world, long num, char type);
void find_tips (node * p, node ** nodelist, long *z);
node *add_migration (world_fmt *world, node * p, long from, long to, MYREAL utime);
node *create_interior_node (world_fmt * world, node ** q);
node *create_root_node (world_fmt *world, node ** q);
node *create_tip_node (FILE * file, world_fmt * world, option_fmt *options, node ** q, char *ch);
boolean processbracket (FILE * file, world_fmt *world, node ** p, char *ch);
void set_tree_pop (node * p, long *pop);
void allocate_x (node * p, world_fmt * world, char datatype,
                 boolean withtips);
long find_firstpop (node * p);

void sankoff (world_fmt * world);
MYREAL minimum (MYREAL *vec1, MYREAL *vec2, long n);
void santraverse (world_fmt *world, node * theNode, MYREAL **cost, long numpop);
long ranbest (MYREAL *array, long tie, MYREAL best, long n);
void jumble (long *s, long n);
long number_genomes (int datatype);

/* copy whole tree */
void copy_tree (world_fmt * original, world_fmt * kopie);
node *copy_node (world_fmt * original, node * o, world_fmt * kopie,
                 node * last);
void copy_node_content (world_fmt * original, world_fmt * kopie, node * o,
                        node * t);
void swap_tree (world_fmt * tthis, world_fmt * tthat);
void swap (void *a, void *b);

void free_tree (node * p, world_fmt * world);
void free_tipnodelet (node * p, world_fmt * world);
void free_mignodelet (node * p, world_fmt * world);
void free_nodelet (node * p, long num, world_fmt * world);
void free_nodedata (node * p, world_fmt * world);

MYREAL inverse_logprob_noevent (world_fmt * world, long interval);
MYREAL sum_migprob (world_fmt * world, long pop, long interval);

MYREAL prob_micro_watkins (MYREAL t, long diff, world_fmt * world, pair *helper);
MYREAL prob_micro_singlestep (MYREAL t, long diff, world_fmt * world, pair *helper);

void
debugtreeout (FILE * file, node * joint, node * p, long s);


/** global variable NUVIEW points to function nuview_datatype() */
static void (*nuview) (node *, world_fmt *, long);
static double (*prob_micro) (MYREAL , long , world_fmt *, pair *);

/* ======================================================= */
/*
 * Creates a start-genealogy using a coalescence approach
 *
 * - set NUVIEW according to datatype
 * initializes tree structure - fills tree with data - set_tree():
 * upgma-tree, adjust for times, sankoff() for migration events
 */
/// \brief builds the starting tree
///
/// Builds the starting tree. 
/// \callgraph
void
buildtree (world_fmt * world, option_fmt * options, data_fmt * data,
           long locus)
{
    long pop;
    long genomes = number_genomes (options->datatype);
    free_tree(world->root, world);
    world->sumtips = 0;
    world->migration_counts = 0;
    for (pop = 0; pop < data->numpop; pop++)
      {
	if(options->randomsubset > 0 && options->randomsubset < data->numind[pop][locus])
	  world->sumtips += options->randomsubset * genomes; 
	else
	  world->sumtips += data->numalleles[pop][locus];
      }
    which_nuview (options->datatype, options->fastlike, FALSE, options->msat_option);

    switch (options->datatype)
      {
      case 's':
      case 'n':
      case 'h':
      case 'u':
      case 'f':
	init_sequences (world, options, data, locus);
	init_sequences_aliases (world, options, data, locus);
	break;
      case 'b':
	world->data->seq[0]->endsite = 1;
	data->freq = -10000000000000.;
	break;
      case 'm':
	world->data->seq[0]->endsite = 1;
	/* world->data->freq = 1. / world->options->micro_stepnum; */
	break;
      case 'a':
	world->data->seq[0]->endsite = 1;
	world->data->freq = 1. / (data->maxalleles[locus]);
	world->data->freqlast = 1. - world->data->freq;
      }
    
    if (options->usertree)
      {
        options->usertreewithmig = treereader (world, options, data);
        makevalues (world, options, data, locus);//insert values into tips
      }
    else
      {
        allocate_tree (world, options, data, locus); //allocate nodep
        makevalues (world, options, data, locus); //allocate/insert data into tips
        allocateinterior(world,data,locus); //allocate nodep guts
        allocatepoproot (world, data, locus); //allocate bottom parts
        if (world->data->skiploci[locus])
	  return;
      }
    //    reorder_populations(world, options, data);

    if (strchr (SEQUENCETYPES, world->options->datatype))
      {
        init_tbl (world, locus);
        if (world->cold && world->replicate == 0)
	  {
            print_seqfreqs (world->outfile, world, options);
            print_tbl (world->outfile, world, options, locus);
            print_weights (world->outfile, world, options, locus);
	  }
        if (world->options->progress)
	  {
	    if (world->cold &&  world->replicate == 0)
	      {
		print_seqfreqs (stdout, world, options);
		print_tbl (stdout, world, options, locus);
		print_weights (stdout, world, options, locus);
		if (world->options->writelog)
		  {
		    print_seqfreqs (world->options->logfile, world, options);
		    print_tbl (world->options->logfile, world, options, locus);
		    print_weights (world->options->logfile, world, options,
				   locus);
		  }
	      }
	  }
      }
    if(!options->usertree)
      set_tree (world, options, data, locus);
    else
      {
        if(!options->usertreewithmig)
	  set_migrations (world, options, data, locus);
      }
    if (world->options->datatype == 'b')
      world->data->maxalleles[locus] = XBROWN_SIZE;
    
    // calculate the migration counts for the first tree
    world->migration_counts = 0;
    count_migrations (world->root->next->back, &world->migration_counts);
    //traverse_check(crawlback (world->root->next));printf("@");

}

/*
 * creates the timelist which represents all time intervals on a tree. The
 * timelist is an array of pointers and not a linked list.
 *
 * - allocates memory using an arbitrary value this will be later adjusted if a
 * longer list is needed - construct
 */
void
create_treetimelist (world_fmt * world,  timelist_fmt ** ltl, long locus)
{
    if ((*ltl)->tl == NULL)
    {
      (*ltl)->allocT =  TIMELIST_GUESS;
      (*ltl)->tl = (vtlist *) mycalloc ((*ltl)->allocT, sizeof (vtlist));
        allocate_lineages (ltl, 0, world->numpop);
    }
    (*ltl)->copies = 0;
    construct_tymelist (world, (*ltl));
}


void
allocate_lineages (timelist_fmt **timevector, const long offset, const long numpop)
{
    long i;
    const long allocT = (*timevector)->allocT;
    vtlist * tl;
    long *lineages;

    // speedup and simplification for allocating and freeing
    // accessing is still done the old way through tl[i].lineages
    // but instead having memory scattered the allocation is one long string
    // freeing time for lineages should be reduced this way [I hope]
    if((*timevector)->lineages != NULL)
      {
	(*timevector)->lineages = (long *) myrealloc((*timevector)->lineages,sizeof(long)*numpop*allocT);
	if(offset != allocT)
	  {
	    memset((*timevector)->lineages+offset,0,sizeof(long)*numpop*(allocT-offset));
	  }
      }
    else
      (*timevector)->lineages = (long *) mycalloc(numpop*allocT,sizeof(long));
   
    tl = (*timevector)->tl;
    lineages = (*timevector)->lineages;

    for (i = offset; i < allocT; i++)
      {
	tl[i].lineages = lineages + i * numpop;	
      }
}

/*
 * start first pass trhough the tree to calculate the tree-likleihood
 */
void
first_smooth (world_fmt * world, long locus)
{
    smooth (world->root->next, crawlback (world->root->next), world, locus);
}

/*
 * Marks a node, so that TREELIKELIHOOD() will recalulated values in node
 */
void
set_dirty (node * p)
{
    p->dirty = TRUE;
}
///
/// inserts the from and to into the timelist from the tree
void
timeslices (timelist_fmt ** timevector)
{
    long z;
    vtlist * tl = (*timevector)->tl;
    vtlist * tlz;
    node * eventnode;
    const long T = (*timevector)->T;
 
    for (z = 0; z < T; z++)
    {
      tlz        = &(tl[z]);
      eventnode  = tlz->eventnode;
      tlz->from  = eventnode->pop;
      tlz->to    = eventnode->actualpop;
      tlz->slice = z;
    }
}

void
add_partlineages (long numpop, timelist_fmt ** timevector)
{
  // this points to the MRCA
  const long T = (*timevector)->T - 2;

  long i, pop;
  vtlist *tl = (*timevector)->tl;
  vtlist *tli;
  vtlist *tli1;
  long from;
  long to;
  long *lineages;
  char type;
  // this should add a lineages for the MRCA
  memset(tl[T+1].lineages, 0, numpop * sizeof(long));  
  tl[T+1].lineages[tl[T].eventnode->pop] = 1;
  for (i = T; i >= 0; i--)
    {
      tli1 = &tl[i+1];
      tli = &tl[i];
      lineages = tli->lineages;
      from = tli->from;
      to = tli->to;
      type = tli->eventnode->type;
      memset(lineages, 0, numpop * sizeof(long));  
      // if the node is an internode (a coalescent node) then add an additional line
      // if it is a tip reduce one
      if(type == 't')
	{
	  lineages[to] -= 1;
	}
      else
	{
	  if(type != 'r')
	    {
	      if (from == to)
		{
		  lineages[to] += 1;
		}
	      else
		{
		  lineages[to] += 1;
		  lineages[from] -= 1;
		}
	    }
	}
      // this copies the content from the last timeinterval (tli1) to the next (tli)
      //      printf("%i> lineages %5li:", myID, i);
      for (pop = 0; pop < numpop; pop++)
	{
	  lineages[pop] += tli1->lineages[pop];
	  //  printf(" %li",lineages[pop]);
	}
      //      printf(" %f %c %li-->%li %s\n",tli->eventnode->tyme, tli->eventnode->type, tli->from, tli->to, tli->eventnode->nayme);
    }
}

void
add_partlineages_313 (long numpop, timelist_fmt ** timevector)
{
  // this points to the MRCA
  const long T = (*timevector)->T - 2;

  long i, pop;
  vtlist *tl = (*timevector)->tl;
  vtlist *tli;
  vtlist *tli1;
  long from;
  long to;
  long *lineages;
  char type;
  // this should add a lineages for the MRCA
  memset(tl[T+1].lineages, 0, numpop * sizeof(long));  
  tl[T+1].lineages[tl[T].eventnode->pop] = 1;
  for (i = T; i >= 0; i--)
    {
      tli1 = &tl[i+1];
      tli = &tl[i];
      lineages = tli->lineages;
      from = tli1->from;
      to = tli1->to;
      type = tli1->eventnode->type;
      memset(lineages, 0, numpop * sizeof(long));  
      // if the node is an internode (a coalescent node) then add an additional line
      // if it is a tip reduce one
      if(type == 't')
	{
	  lineages[to] -= 1;
	}
      else
	{
	  if(type != 'r')
	    {
	      if (from == to)
		{
		  lineages[to] += 1;
		}
	      else
		{
		  lineages[to] += 1;
		  lineages[from] -= 1;
		}
	    }
	}
      // this copies the content from the last timeinterval (tli1) to the next (tli)
      //      printf("%i> lineages %5li:", myID, i);
      for (pop = 0; pop < numpop; pop++)
	{
	  lineages[pop] += tli1->lineages[pop];
	  //  printf(" %li",lineages[pop]);
	}
      //      printf(" %f %c %li-->%li %s\n",tli->eventnode->tyme, tli->eventnode->type, tli->from, tli->to, tli->eventnode->nayme);
    }
}

void
add_partlineages_312 (long numpop, timelist_fmt ** timevector)
{
  // this points to the MRCA
  const long T = (*timevector)->T - 2;

  long i, pop;
  vtlist *tl = (*timevector)->tl;
  vtlist *tli;
  vtlist *tli1;
  long from;
  long to;
  long *lineages;
  char type;
  // this should add a lineages for the MRCA
  memset(tl[T+1].lineages, 0, numpop * sizeof(long));  
    tl[T].lineages[tl[T].eventnode->pop] = 1;// the tMRCA has then 1 lineages looking to the root
  tl[T+1].lineages[tl[T].eventnode->pop] = 0;// the tMRCA has then 1 lineages looking to the root
  for (i = T; i > 0; i--)
    {
      tli1 = &tl[i-1];
      tli = &tl[i];
      lineages = tli1->lineages;
      from = tli->from;
      to = tli->to;
      type = tli->eventnode->type;
      memset(lineages, 0, numpop * sizeof(long));  
      // if the node is an internode (a coalescent node) then add an additional line
      // if it is a tip reduce one
      if(type == 't')
	{
	    lineages[to] -= 1;
	}
      else
	{
	  if(type != 'r')
	    {
	      if (from == to)
		{
		  lineages[to] += 1;
		}
	      else
		{
		  lineages[to] += 1;
		  lineages[from] -= 1;
		}
	    }
	}
      // this copies the content from the last timeinterval (tli1) to the next (tli)
      //printf("%i> lineages %5li:", myID, i);
      for (pop = 0; pop < numpop; pop++)
	{
	  lineages[pop] += tli->lineages[pop];
	  //printf(" %li",lineages[pop]);
	}
      //      printf(" %f %c %li-->%li %s lin=%li\n",tli->eventnode->tyme, tli->eventnode->type, tli->from, tli->to, tli->eventnode->nayme, lineages[0]);
    }
}


/*
 calculates the tree-likelihood according to datatype a, m, b, s,u
 */
MYREAL
treelikelihood (world_fmt * world)
{
    long a;
    MYREAL term = 0.0;
    //MYREAL term1 = 0.0;
    node *nn = crawlback (world->root->next);
    //#ifdef BEAGLE
    //reset_beagle(world->beagle);
    //#endif
#ifdef BEAGLE
    printf("TreeLike:");
    term = calcLnL(world,FALSE);
    //    printf("Start LnL from Beagle: %f\n",term);
    return term;
#endif
    set_dirty (nn);
    smooth (world->root->next, crawlback(world->root->next), world, world->locus);
    adjustroot (world->root);
    switch (world->options->datatype)
    {
    case 's':
      term = treelike_seq (world, world->locus);
      break;
    case 'n':
    case 'h':
      term = treelike_snp (world, world->locus);
      break;
    case 'u':
      term = treelike_snp_unlinked (world, world->locus);
      break;
    case 'a':
      //		  printf("@@@ ");
      for (a = 0; a < world->data->maxalleles[world->locus] - 1; a++)
	{
	  term += (world->data->freq * nn->x.a[a]);
	  //  printf("%f ",term);
	}
      term += (world->data->freqlast * nn->x.a[a]);
      //		  printf("%f ",term);
      term = (term != 0.0) ? (LOG (term) + nn->scale[0]) : -MYREAL_MAX;
      //		  printf("\n@@treelike=%f scale=%f\n",term, nn->scale[0]);
      break;
    case 'm':
      for (a = 0; a < world->data->maxalleles[world->locus]; a++)
	term += nn->x.a[a];
      term = (term != 0.0) ? (LOG (term) + nn->scale[0]) : -MYREAL_MAX;
      break;
    case 'b':
      term = nn->x.a[2];
      break;
    case 'f':
      term = treelike_anc (world, world->locus);
      break;
    }
#ifdef UEP
    if (world->options->uep)
        ueplikelihood (world);
#endif
#ifdef DEBUG
    //     printf("@L(D|G)=%f\n",term);
#endif
    return term;
}

/*
 * calculates tree-likelihood using only arrays DOES NOT CHANGE ARRAYS IN THE
 * TREE
 */
MYREAL
pseudotreelikelihood (world_fmt * world, proposal_fmt * proposal)
{
  long a, locus = world->locus;
  /* freq is not different between pop */
  MYREAL term = 0.0;
  switch (world->options->datatype)
    {
    case 's':
#ifdef BEAGLE
      printf("PseudoLike: ");
      return calcLnL(world,world->beagle->instance);
#endif
      term = pseudo_tl_seq (proposal->xf.s, proposal->xt.s, proposal->v,
			    proposal->vs, proposal, world);

      break;
    case 'n':
    case 'h':
      term = pseudo_tl_snp (proposal->xf.s, proposal->xt.s, proposal->v,
			    proposal->vs, proposal, world);
      break;
    case 'u':
      term = pseudo_tl_snp_unlinked (proposal->xf.s, proposal->xt.s,
				     proposal->v, proposal->vs, proposal,
				     world);
      break;
    case 'a':
      //		  printf("@@@");
      for (a = 0; a < world->data->maxalleles[locus] - 1; a++)
	{
	  term += (world->data->freq * proposal->xf.a[a]);
	  //printf("%f ",term);
	}
      term += (world->data->freqlast * proposal->xf.a[a]);
      //printf("%f ",term);
      if (term == 0.0)
	term = -MYREAL_MAX;
      else
	term = (LOG (term) + proposal->mf[0]);
      //			printf("\n@@pseudolike=%f scale=%f\n",term, proposal->mf[0]);
      break;
    case 'b':
      term = proposal->xf.a[2];
      break;
    case 'm':
      for (a = 0; a < world->data->maxalleles[locus]; a++)
	{
	  term += proposal->xf.a[a];
	}
      if (term == 0.0)
	term = -MYREAL_MAX;
      else
	term = (LOG (term) + proposal->mf[0]);
      break;
    case 'f':
      term = pseudo_tl_anc (proposal->xf.s, proposal->xt.s, proposal->v,
			    proposal->vs, proposal, world);
      break;
    default:
      warning("no datatype found for treelikelihood() calculation!");
      term = -MYREAL_MAX;
      break;
    }
#ifdef UEP
  if (proposal->world->options->uep)
    term += pseudo_tl_uep (&(proposal->uf), &proposal->ut, proposal->v,
			   proposal->vs, proposal, world);
#endif
  if(MYFINITE(((double) term)) == 0)
    term = -MYREAL_MAX;
#ifdef DEBUG
  //    printf("   proposed L(D|G)=%f\n",term);
#endif
  return term;
}


/*
 * Calculates the sub-likelihoods but does not change the arrays in the tree,
 * it uses the passed arrays and overwrites the xx1 array DOES NOT CHANGE THE
 * TREE
 */
void
pseudonuview (proposal_fmt * proposal, xarray_fmt xx1, MYREAL *lx1, MYREAL v1,
              xarray_fmt xx2, MYREAL *lx2, MYREAL v2)
{
    switch (proposal->datatype)
    {
    case 'a':
      pseudonu_allele (proposal, &xx1.a, &(lx1[0]), v1, xx2.a, lx2[0], v2);
      break;
    case 'b':
      pseudonu_brownian (proposal, &xx1.a, &(lx1[0]), v1, xx2.a, lx2[0], v2);
      break;
    case 'm':
      pseudonu_micro (proposal, &xx1.a, &(lx1[0]), v1, xx2.a, lx2[0], v2);
      break;
    case 'u':
    case 'n':
    case 'h':
      pseudonu_seq (proposal, xx1.s, v1, xx2.s, v2);
      break;
    case 's':
      //#ifdef BEAGLE
      //prepare_beagle_instances_proposal(proposal,xx1.s, v1,xx2.s, v2, world->beagle);
      //#else
      if (proposal->world->options->fastlike)
	pseudonu_seq (proposal, xx1.s, v1, xx2.s, v2);
      else
	pseudonu_seq_slow (proposal, xx1.s, lx1, v1, xx2.s, lx2, v2);
      //#endif
      break;
    case 'f':
      pseudonu_anc (proposal, xx1.s, v1, xx2.s, v2);
      break;
    }
#ifdef  UEP
    if (proposal->world->options->uep)
    {
        pseudonu_twostate (proposal, &proposal->uf, proposal->umf, v1,
                           &proposal->ut, proposal->umt, v2);
    }
#endif
}


MYINLINE void
pseudonu_allele (proposal_fmt * proposal, MYREAL **xx1, MYREAL *lx1,
                 MYREAL vv1, MYREAL *xx2, MYREAL lx2, MYREAL vv2)
{
    
    long   a;
    long   aa;
    long   locus = proposal->world->locus; /* allele counters */
    long   mal   = proposal->world->data->maxalleles[locus]; /* maxalleles */
    MYREAL freq  = proposal->world->data->freq;
    //MYREAL freqlast = proposal->world->data->freqlast; // same as 1 - k/(k+1) = 1/(k+1) = freq
    MYREAL w1 = 0.0; /* time variables */
    MYREAL w2 = 0.0;
    MYREAL v1;
    MYREAL v2;
    MYREAL v1freq;
    MYREAL v2freq;
    MYREAL pija1;/* summary of probabilities */
    MYREAL pija2;  
    MYREAL x3m = -MYREAL_MAX;
    MYREAL inv_x3m;
    MYREAL *xx3;
    xx3 = (MYREAL *) mymalloc(sizeof(MYREAL) * mal);
    v1 = 1.0 - EXP (-vv1);
    v2 = 1.0 - EXP (-vv2);
    if (v1 >= 1.)
    {
        w1 = 0.0;
        v1 = 1.0;
    }
    else
    {
        w1 = 1.0 - v1;
    }
    if (v2 >= 1.)
    {
        w2 = 0.0;
        v2 = 1.0;
    }
    else
    {
        w2 = 1.0 - v2;
    }
    //    printf("@@pseudonu 1:{%f,%f,%f,%f}   2:(%f,%f,%f,%f}\n",(*xx1)[0],(*xx1)[1],*lx1,v1,
    //	   xx2[0], xx2[1],lx2,v2);
  
    v1freq = v1 * freq;
    v2freq = v2 * freq;

    for (aa = 0; aa < mal; aa++)
    {
        pija1 = pija2 = 0.0;
        for (a = 0; a < mal; a++)
        {
	  if(aa==a)
	    {
	      pija1 += (w1 + v1freq) * (*xx1)[a];
	      pija2 += (w2 + v2freq) * xx2[a];
	    }
	  else
	    {
	      pija1 += v1freq * (*xx1)[a];
	      pija2 += v2freq * xx2[a];
	    }
	}
        xx3[aa] = pija1 * pija2;
        if (xx3[aa] > x3m)
            x3m = xx3[aa];
    }
    inv_x3m = 1./ x3m;
    for (aa = 0; aa < mal; aa++)
    {
        (*xx1)[aa] = xx3[aa] * inv_x3m;
    }
    //    printf("@@finish: %f %f\n",(*xx1)[0],(*xx1)[1]);
    *lx1 = LOG (x3m) + lx2 + *lx1;
    myfree(xx3);
}

MYINLINE 
void
pseudonu_micro (proposal_fmt * proposal, MYREAL **xx1, MYREAL *lx1, MYREAL v1,
                MYREAL *xx2, MYREAL lx2, MYREAL v2)
{
    long a, s, diff, locus = proposal->world->locus; /* allele counters */
    long aa1, aa2;
    long smax = proposal->world->data->maxalleles[locus];
    long margin = proposal->world->options->micro_threshold[locus];
    MYREAL pija1s, pija2s, vv1, vv2;
    MYREAL x3m = -MYREAL_MAX;
    world_fmt *world = proposal->world;
    MYREAL inv_x3m;
    MYREAL *xx3;
    MYREAL *pm1;
    MYREAL *pm2;
    pair *helper = &world->options->msat_tuning;
    xx3 = (MYREAL *) mymalloc(sizeof(MYREAL) * smax);
    vv1 = v1;
    vv2 = v2;

    pm1 = (MYREAL *) mymalloc(sizeof(MYREAL) * 2 * margin);
    pm2 = pm1 +  margin;

    for (diff = 0; diff < margin; diff++)
      {
	pm1[diff] = prob_micro (vv1, diff, world, helper);
	pm2[diff] = prob_micro (vv2, diff, world, helper);
      }
    
    for (s = 0; s < smax; s++)
    {
        pija1s = pija2s = 0.0;
	aa1 = MAX (0, s - margin);
	aa2 = MIN(s + margin,smax);
	for (a = aa1; a < aa2; a++)
	  //for (a = 0; a < smax; a++)
        {
            diff = labs (s - a);
	    if(diff >= margin)
	      continue;

            if ((*xx1)[a] > 0)
            {
                pija1s += pm1[diff] * (*xx1)[a];
            }
            if (xx2[a] > 0)
            {
                pija2s += pm2[diff] * xx2[a];
            }
        }
        xx3[s] = pija1s * pija2s;
        if (xx3[s] > x3m)
            x3m = xx3[s];
    }
    inv_x3m = 1./ x3m;
    for (s = 0; s < smax; s++)
    {
        (*xx1)[s]= xx3[s] * inv_x3m;
    }
    *lx1 += LOG (x3m) + lx2;
    myfree(xx3);
    myfree(pm1);
}

//================================================
// brownian motion calculation for data likelihood: ghost-version
// changed divisions to multiply inverse
// simplified check and include of boundary for vtot
MYINLINE
void
pseudonu_brownian (proposal_fmt * proposal, MYREAL **xx1, MYREAL *lx1,
                   MYREAL v1, MYREAL *xx2, MYREAL lx2, MYREAL v2)
{
    MYREAL vtot, rvtot, c12;
    MYREAL mean1, mean2, mean, vv1, vv2, f1, f2, diff;
    mean1 = (*xx1)[0];
    mean2 = xx2[0];
	
    vv1 = v1 + (*xx1)[1];
    vv2 = v2 + xx2[1];
    vtot = vv1 + vv2;
    if (vtot > 0.0)
    {
        rvtot = 1./vtot;
        f1 = vv2 * rvtot;
        f2 = 1.0 - f1;
        mean = f1 * mean1 + f2 * mean2;
        diff = mean1 - mean2;
        c12 = diff * diff * rvtot;
        (*xx1)[2] = (*xx1)[2] + xx2[2] + MIN (0, -0.5 * (LOG (vtot) + c12) + LOG2PIHALF);
        (*xx1)[1] = vv1 * f1;
        (*xx1)[0] = mean;
    }
    else
    {
        //xcode rvtot = HUGE;
        //xcode vtot = SMALL_VALUE;
        f1 = 0.5;
        (*xx1)[2] = HUGE;
        (*xx1)[1] = vv1 * f1;
        (*xx1)[0] = f1 * (mean1 + mean2);
    }
}

/*
 adjust the variables POP and ACTUALPOP in interior nodes
 */
void
set_pop (node * theNode, long pop, long actualpop)
{
    switch (theNode->type)
    {
		case 'm':
			theNode->pop = theNode->next->pop = pop;
			theNode->actualpop = theNode->next->actualpop = actualpop;
			break;
		case 'i':
		case 'r':
			theNode->pop = theNode->next->pop = theNode->next->next->pop = pop;
			theNode->actualpop = theNode->next->actualpop = actualpop;
			theNode->next->next->actualpop = actualpop;
			break;
		case 't':
			if (theNode->pop != pop)
				error ("Population designation scrambled");
			break;
		default:
			error ("Undefined node type?!");
			break;
    }
}



/*
 * ======================================================= local functions
 */

void
allocate_tree(world_fmt * world, option_fmt * options, data_fmt * data,
              long locus)
{
    long nodenum = 0, pop, numpop = data->numpop;
    for (pop = 0; pop < numpop; pop++)
    {
        nodenum += data->numalleles[pop][locus] * 2;
    }
    world->nodep = (node **) mycalloc (nodenum+1, sizeof (node *));
}

void
allocatetips (world_fmt * world, option_fmt * options,
              data_fmt * data, long locus)
{
    long a, pop;
    long zz = 0;
    for (pop=0; pop < data->numpop; pop++)
    {
        for(a=0; a < data->numalleles[pop][locus]; a++)
        {
	  allocate_tip (world, options, &world->nodep[zz], pop, locus, zz, a,
			data->indnames[pop][a]);
            zz++;
        }
    }
}

void
allocateinterior (world_fmt * world, data_fmt * data, long locus)
{
  long sites = world->data->seq[0]->endsite + world->data->seq[0]->addon;
    node *p;
    long i;
    long temp=0;
    long mini = 0, maxi = 0;
    long numpop = data->numpop;
    for (i = 0; i < numpop; i++)
    {
        temp += data->numalleles[i][locus];
    }
    mini = temp;
    maxi = temp * 2 - 1; 
    for (i = mini; i < maxi; i++)
    {
      p = allocate_nodelet (world, 3, 'i');
        p->top = TRUE;
        p->scale =
            (MYREAL *) mycalloc (sites, sizeof (MYREAL));
        p->s = (MYREAL *) mycalloc (1, sizeof (MYREAL) * world->numpop);
	p->id = i;
        world->nodep[i] = p;
    }
}

node *
allocate_nodelet (world_fmt *world, long num, char type)
{
  boolean isfirst = TRUE;
  long j;
  node *p, *q = NULL, *pfirst = NULL;
  for (j = 0; j < num; j++)
    {
      p = dispense_nodelet(world);
      if(p==NULL)
	{
	  warning( "%i> heat=%f Nodelet dispenser failed to supply nodelet\n",myID,1./world->heat);
	  p = (node *) mymalloc (sizeof (node));
	}
      p->tip = FALSE;
      p->visited = FALSE;
      p->number = unique_id_global;
      p->pop = -1;
      p->actualpop = -1;
      p->type = type;
      p->id = -1;
      p->top = FALSE;
      p->dirty = TRUE;
      p->next = q;
      p->scale = NULL;
      p->s = NULL;
      p->x.s = NULL;
      p->x.a = NULL;
#ifdef UEP
      
      p->uep = NULL;
      p->ux.s = NULL;
      p->ux.a = NULL;
#endif
      p->back = NULL;
      p->nayme = NULL;
      p->v = 0.0;
      p->tyme = 0.0;
      p->length = 0.0;
      if (isfirst)
        {
	  isfirst = FALSE;
	  pfirst = p;
        }
      q = p;
    }
    if(pfirst !=NULL)
	  pfirst->next = q;
  return q;
}

void
allocatepoproot (world_fmt * world, data_fmt * data, long locus)
{
  long sites = world->data->seq[0]->endsite + world->data->seq[0]->addon;
  //long i;
  node *p, *q;//, *qq;
    //long nodenum = 0;  /* arbitrarily to the first */
    q = NULL;
    //p = allocate_nodelet (world, 3, 'i');
    //p->top = TRUE;
    //for (i = 0; i < world->numpop; i++)
    //    nodenum += data->numalleles[i][locus];
    //p->x.a = (MYREAL *) mycalloc (nodenum + 1, sizeof (MYREAL));
    //p->top = TRUE;
    //qq = p;
    //q = NULL;
    p = allocate_nodelet (world, 3, 'r');
    p->top = TRUE;
    p->scale = (MYREAL *) mycalloc (sites, sizeof (MYREAL));
    //p->next->back = qq;
    //qq->back = p->next;
    world->root = p;
}


void
allocate_tip (world_fmt * world, option_fmt * options, node ** p, long pop,
              long locus, long a, long ind, char **tipnames)
{
  const long sites = world->data->seq[0]->endsite + world->data->seq[0]->addon;	
  long i;
  char *tipname;
  if(options->usertree)
        return; //do nothing because the tip is alread allocated
    if(tipnames[locus][0]!='\0')
      tipname = tipnames[locus];
    else
      tipname = tipnames[0];
    (*p) = allocate_nodelet (world, 1, 't');
    (*p)->tip = TRUE;
    (*p)->top = TRUE;
    (*p)->id  = a;
    if(options->has_datefile)
      {
	(*p)->tyme = find_tipdate(tipname, pop, world);
      }
    else
      (*p)->tyme = 0.0;
    (*p)->scale = (MYREAL *) mycalloc (sites, sizeof (MYREAL));
    (*p)->pop = (*p)->actualpop = options->newpops[pop]-1;
    (*p)->s = (MYREAL *) mycalloc (1, sizeof (MYREAL) * world->numpop);
    for (i = 0; i < world->numpop; i++)
      (*p)->s[i] = MYREAL_MAX;
    (*p)->s[options->newpops[pop]-1] = 0;
    if (strchr (SEQUENCETYPES, world->options->datatype))
      {	
        (*p)->nayme =
	  (char *) mycalloc (1, sizeof (char) * (options->nmlength + 3));
	strncpy((*p)->nayme,tipname,options->nmlength);
	unpad((*p)->nayme," ");
	translate((*p)->nayme,' ', '_');
        alloc_seqx (world, (*p));
      }
    else
      {
        (*p)->x.a =
	  (MYREAL *) mycalloc (1,
			       world->data->maxalleles[locus] * sizeof (MYREAL));
        (*p)->nayme =
	  (char *) mycalloc (1, sizeof (char) * (options->nmlength + DEFAULT_ALLELENMLENGTH + 5));//name<A|B>:allele
	strncpy((*p)->nayme,tipname,options->nmlength);
	unpad((*p)->nayme," ");
	translate((*p)->nayme,' ', '_');
      }
#ifdef UEP
    if (world->options->uep)
      {
        (*p)->uep = (int *) mycalloc (world->data->uepsites, sizeof (int));
        (*p)->ux.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
      }
#endif
}


void
makevalues (world_fmt * world, option_fmt * options, data_fmt * data,
            long locus)
{
    switch (world->options->datatype)
    {
    case 'a':
      make_alleles (world, options, data, locus);
      break;
    case 'b':
      make_microbrownian (world, options, data, locus);
      break;
    case 'm':
      make_microsatellites (world, options, data, locus);
      break;
    case 's':
    case 'f':
      //  allocatetips (world, options, data, locus);
      make_sequences (world, options, data, locus);
      world->data->seq[0]->oldsite = world->data->seq[0]->endsite;
      break;
    case 'n':
    case 'h':
      //allocatetips (world, options, data, locus);
      make_sequences (world, options, data, locus);
      make_invarsites (world, data, locus);
      world->data->seq[0]->oldsite = world->data->seq[0]->endsite;
      world->data->seq[0]->endsite += 4;
      break;
    case 'u':
      allocatetips (world, options, data, locus);
      make_snp (world, options, data, locus);
      make_invarsites_unlinked (world, data, locus);
      world->data->seq[0]->oldsite = world->data->seq[0]->endsite;
      world->data->seq[0]->endsite *= (data->seq[0]->addon + 1);
      break;
    default:
      usererror ("Oh yes, it would be nice if there were more\n \
possible datatypes than just an\n				\
allele model, microsatellite model,\n sequence model, or	\
single nucleotide polymorphism model.\n				\
But there are currently no others, so the programs stops\n\n");
			break;
    }
#ifdef UEP
    if (world->options->uep)
    {
        make_uep_values (world, data, locus);
    }
#endif
}

///
///
/// creates the branchlength and adds migration nodes to the
/// start tree
/// - creates rough genetic distance for upgma
/// - upgma
/// - find branches where we need to insert migrations
/// - insert migrations
/// - adjust time of all nodes using the coalescent with migration
void
set_tree (world_fmt * world, option_fmt * options, data_fmt * data,
          long locus)
{
    long     tips = world->sumtips;
    MYREAL   **distm;
    node     **topnodes;

    topnodes = (node **) mycalloc (1, sizeof (node *) * tips);
    doublevec2d(&distm,tips, tips);
        
    if (!options->randomtree)
      {
	/* create a crude distance matrix according to the datatype */
	switch (world->options->datatype)
	  {
	  case 'a':
	  case 'b':
	  case 'm':
	    distance_allele (world, options, locus, tips, distm);
	    break;
	  case 's':
	  case 'n':
	  case 'h':
	  case 'u':
	  case 'f':
	    distance_sequence (data, locus, tips,
			       world->data->seq[0]->sites[locus],
			       options->nmlength, distm);
	    break;
	  }
      }
    else
      {
	randomize_distm (distm, data, locus, tips,
			 world->options->custm);
	//printf("\n\nRANDOM TREE START\n\n");
      }
    constrain_distance_zeromig (distm, options, data, locus, tips,
				world->options->custm);
    //    printf("distm[0][1]=%f distm[1][2]=%f\n",distm[0][1],distm[1][2]);
    //fflush(stdout);

#ifdef UEP
	if (options->uep)
	  constrain_distance_uep (data->uep, world->data->uepsites, distm,
				  tips);
#endif
    if(!options->usertree)
      upgma (world, distm, tips, world->nodep);
    myfree(distm[0]);
    myfree(distm);
    //}
    world->root->tyme = world->root->next->tyme =
      world->root->next->next->tyme = world->root->next->back->tyme + 10000.;
    /* orient the tree up-down, set the length and v */
    set_top (world, world->root->next->back, locus);
    set_v (world->root->next->back);
#ifdef BEAGLE
    long bid=0;
    bid = set_branch_index (world->root->next->back, &bid);
#endif
    /*
     * insert migration nodes into the tree using the Slatkin and
     * Maddison approach (Fitch parsimony)
     */
    memcpy (topnodes, world->nodep, sizeof (node *) * tips);
    //zzz = 0;
    if(!options->usertree) //TRIAL
      allocate_x (world->root, world, world->options->datatype, NOTIPS);
#ifdef UEP
    
    if (world->options->uep)
      {
	//      allocate_uep (world->root, world, world->options->datatype, NOTIPS);
	update_uep (world->root->next->back, world);
	check_uep_root (world->root->next->back, world);
      }
#endif
    //    debugtreeout (stdout, crawlback (world->root->next), crawlback (world->root->next),0);
    sankoff (world);
    //    debugtreeout (stdout, crawlback (world->root->next), crawlback (world->root->next),0);
    myfree(topnodes);
}    /* set_tree */

/*
 creates the branchlength and adds migration nodes to the
 start tree
 - creates rough genetic distance for upgma
 - upgma
 - find branches where we need to insert migrations
 - insert migrations
 - adjust time of all nodes using the coalescent with migration
 */
void
set_migrations (world_fmt * world, option_fmt * options, data_fmt * data,
				long locus)
{
    long tips = world->sumtips;
	
    node **topnodes;
    topnodes = (node **) mycalloc (1, sizeof (node *) * tips);
	// can we do this here ???#ifdef UEP
	//            if (options->uep)
	//            constrain_distance_uep (data->uep, world->data->uepsites, distm,  tips);
	//#endif
	world->root->tyme = world->root->next->tyme =
		world->root->next->next->tyme = world->root->next->back->tyme + 10000.;
	/* orient the tree up-down, set the length and v */
	set_top (world, world->root->next->back, locus);
	set_v (world->root->next->back);
	/*
	 * insert migration nodes into the tree using the Slatkin and
	 * Maddison approach (Fitch parsimony)
	 */
	memcpy (topnodes, world->nodep, sizeof (node *) * tips);
#ifdef UEP
	if (world->options->uep)
	{
		//      allocate_uep (world->root, world, world->options->datatype, NOTIPS);
		update_uep (world->root->next->back, world);
		check_uep_root (world->root->next->back, world);
	}
#endif
	sankoff (world);
	myfree(topnodes);
}    /* set_migrations */


void randomize_distm(MYREAL **distm, data_fmt * data, long locus,
                            long tips, char *custm)
{
  long i;
  long j;
  for(i=0;i<tips; i++)
    {
      distm[i][i] = 0.0;
      for(j=0;j<i;j++)
	{
	  distm[i][j] = distm[j][i] = RANDUM();
	  //printf("%f ",distm[i][j]);
	}
      //printf("\n");
    }
}


void
constrain_distance_zeromig (MYREAL **m, option_fmt *options, data_fmt * data, long locus,
                            long tips, char *custm)
{
    long pop, ind, j, i = 0;
    long *pops;
    pops = (long *) mycalloc (tips, sizeof (long));
    for (pop = 0; pop < data->numpop; pop++)
    {
        for (ind = 0; ind < data->numalleles[pop][locus]; ind++)
	  {
            pops[i++] = options->newpops[pop];
	    //OLD	    pops[i++] = pop;
	  }
    }
    for (i = 0; i < tips; i++)
    {
        for (j = 0; j < i; j++)
        {
            if (custm[pops[i] + pops[i] * pops[j]] == '0')
                m[i][j] = m[j][i] = 1000;
        }
    }
    myfree(pops);
}

void
distance_EP (char **data, long tips, MYREAL **m)
{
    long i, j;
    for (i = 0; i < tips; i++)
    {
        for (j = 0; j < i; j++)
        {
            if (!strcmp (data[i], data[j]))
                m[i][j] = m[j][i] = fabs (rannor (1., 0.1));
            else
                m[i][j] = m[j][i] = fabs (rannor (0., 0.1));
        }
    }
}

void
distance_micro (char **data, long tips, MYREAL **m)
{
    long i, j;
    for (i = 0; i < tips; i++)
    {
        for (j = 0; j < i; j++)
        {
            m[i][j] = m[j][i] = pow (atof (data[i]) - atof (data[j]), 2.);
            m[i][j] = m[j][i] = fabs (rannor (m[i][j], 0.1));
        }
    }
}

/* calculate pairwise distances using allele disimilarity */
void
distance_allele (world_fmt * world, option_fmt * options, long locus,
                 long tips, MYREAL **distm)
{
    char **mdata;
    long pop;
	
    mdata = (char **) mycalloc (1, sizeof (char *) * (tips + 1));
    for (pop = 0; pop < tips; pop++)
    {
        mdata[pop] =
		(char *) mycalloc (1, sizeof (char) * options->allelenmlength);
        strcpy (mdata[pop], strrchr(world->nodep[pop]->nayme,'!')+1);
    }
    if (world->options->datatype == 'a')
        distance_EP (mdata, tips, distm);
    else
        distance_micro (mdata, tips, distm);
    for (pop = 0; pop < tips; pop++)
    {
        myfree(mdata[pop]);
    }
    myfree(mdata);
}

/* calculate  pairwise distances using sequence similarity */
void
distance_sequence (data_fmt * data, long locus, long tips, long sites,
                   long nmlength, MYREAL **m)
{
    long i = 0, j, z, pop, ind;
    char **dat;
    if (data->distfile != NULL)
    {
        read_distance_fromfile (data->distfile, tips, nmlength, m);
    }
    else
    {
        dat = (char **) mymalloc (sizeof (char *) * tips);
        for (pop = 0; pop < data->numpop; pop++)
        {
            for (ind = 0; ind < data->numalleles[pop][locus]; ind++)
            {
                dat[i++] = data->yy[pop][ind][locus][0];
            }
        }
        if (i != tips)
        {
            error ("Mistake in distance_sequence() tips is not equal sum(i)\n");
        }
        for (i = 0; i < tips; i++)
        {
			
            for (j = i + 1; j < tips; j++)
            {
                //to come m[i][j] = m[j][i] = make_ml_distance(dat[i], dat[j], i, j);
				
                for (z = 0; z < sites; z++)
                {
					
                    if (dat[i][z] != dat[j][z])
                    {
                        //m[i][j] = m[j][i] += fabs(rannor(1.0, 0.1));
                        m[i][j] = m[j][i] += 1.0;
                    }
                }
            }
        }
        myfree(dat);
    }
}



void
fix_times (world_fmt * world, option_fmt * options)
{
    long k;
    MYREAL tipdate=0.0;
    node *theNode;
    MYREAL age = world->data->maxsampledate 
      * world->options->meanmu[world->locus]
      * world->options->generation_year  * world->options->mu_rates[world->locus];
    if (!options->usertree)
    {
#ifdef TREEDEBUG
      printf("%i> %s\nage= %f\nmaxsamp=%f\nmut=%f\ngen=%f\nmu_rate=%g\n-----------------------------\n",myID, "start timelist ---------------------", age,world->data->maxsampledate,world->options->meanmu[world->locus],1./world->options->generation_year,world->options->mu_rates[world->locus]);
#endif
      for (k = 0; k < world->treetimes[0].T - 1; k++)
	{
	  theNode = world->treetimes[0].tl[k].eventnode;
	  if(theNode->type == 't')
	    {
	      //tipdate = find_tipdate(theNode->nayme, theNode->pop, world);
	      tipdate = theNode->tyme;
#ifdef TREEDEBUG
	      printf("DEBUG TIPDATE %i> %s %g\n",myID, theNode->nayme, tipdate);
#endif
	      //theNode->tyme = tipdate;
	      world->treetimes[0].tl[k].age = tipdate;	
	    }
	  else
	    {
	      //k here is never 0, so this k-1 seems correct
	      age += inverse_logprob_noevent (world, k-1);
	      world->treetimes[0].tl[k].age = age;	
	      adjust_time (theNode, age);
#ifdef TREEDEBUG
	      printf("%i> %10.10s %g\n",myID, " ", theNode->tyme);
#endif
	    }
        }
      world->treetimes[0].tl[k].age = age + 10000.; /* this is the root */
      adjust_time (world->treetimes[0].tl[k].eventnode, age+10000);
#ifdef TREEDEBUG
      printf("%i> %s %g\n",myID, "end timelist ---------------------", age+10000 );
#endif
    }
    set_v (world->root->next->back);
}

/*
 * creates a UPGMA tree: x     = distance matrix which will be destroyed
 * through the process, tips  = # of sequences/alleles, nodep = treenodes
 * have to be allocated for ALL nodes
 *
 * This code is stripped neighbor-joining code out of phylip v3.6. Only the
 * upgma option is present.
 */
void
upgma (world_fmt * world, MYREAL **x, long tips, node ** nodep)
{
    long nc, nextnode, mini = -900, minj = -900, i, j, jj, ia, iaa, ja, jaa; 
    MYREAL zz = (world->data->maxsampledate 
		 * world->options->meanmu[world->locus]
		 * world->options->generation_year  * world->options->mu_rates[world->locus]+ 1);
    MYREAL total, tmin, bi, bj, /* ti, tj, */ da;
    MYREAL *av;
    long *oc;
    node **cluster;
    long *enterorder;
	MYREAL denom=1.;
    /* First initialization */
    enterorder = (long *) mycalloc (tips, sizeof (long));
    for (ia = 0; ia < tips; ia++)
        enterorder[ia] = ia;
    jumble (enterorder, tips);
    nextnode = tips;
    av = (MYREAL *) mycalloc (tips, sizeof (MYREAL));
    oc = (long *) mymalloc (tips * sizeof (long));
    cluster = (node **) mycalloc (tips, sizeof (node *));
    for (i = 0; i < tips; i++)
        oc[i] = 1;
    for (i = 0; i < tips; i++)
        cluster[i] = nodep[i];
    /* Enter the main cycle */
    for (nc = 0; nc < tips - 1; nc++)
    {
        tmin = 99999.0;
        /* Compute sij and minimize */
        for (jaa = 1; jaa < tips; jaa++)
        {
            ja = enterorder[jaa];
            if (cluster[ja] != NULL)
            {
                for (iaa = 0; iaa < jaa; iaa++)
                {
                    ia = enterorder[iaa];
                    if (cluster[ia] != NULL)
                    {
                        total = x[ia][ja];
			//printf("ia=%li ja=%li x[ia][ja]=%f\n",ia,ja,x[ia][ja]);
                        if (total < tmin)
                        {
                            tmin = total;
                            mini = ia;
                            minj = ja;
                        }
                    }
                }
            }
        }   /* compute lengths and print */
        bi = x[mini][minj] / 2.0 - av[mini];
        bj = x[mini][minj] / 2.0 - av[minj];
        av[mini] += bi;
        nodep[nextnode]->next->back = cluster[mini];
      	//xcode
        if(cluster[mini]!=NULL)
        {
            if(nodep[nextnode]->next!=NULL)
            {
                cluster[mini]->back = nodep[nextnode]->next;
            }
            nodep[nextnode]->next->next->back = cluster[minj];
            cluster[minj]->back = nodep[nextnode]->next->next;
            cluster[mini]->back->v = cluster[mini]->v = bi;
            cluster[minj]->back->v = cluster[minj]->v = bj;
            cluster[mini] = nodep[nextnode];
        }
        adjust_time (nodep[nextnode], (MYREAL) zz++);
        cluster[minj] = NULL;
        nextnode++;
        /* re-initialization */
        denom = (oc[mini] + oc[minj]);
        for (jj = 0; jj < tips; jj++)
        {
            if (cluster[jj] != NULL)
            {
                da = (x[mini][jj] * oc[mini] + x[minj][jj] * oc[minj])/denom ;
                x[mini][jj] = da;
                x[jj][mini] = da;
            }
        }
        for (j = 0; j < tips; j++)
        {
            x[minj][j] = x[j][minj] = 0.0;
        }
        oc[mini] += oc[minj];
    }
    /* the last cycle */
    for (i = 0; i < tips; i++)
    {
        if (cluster[i] != NULL)
            break;
    }
    world->root->next->back = cluster[i];
    cluster[i]->back = world->root->next;
    myfree(av);
    myfree(oc);
    myfree(cluster);
    myfree(enterorder);
}

void
set_top (world_fmt * world, node * p, long locus)
{
	
	
    if (p->type == 't')
    {
        p->top = TRUE;
        //p->tyme = 0.0;
        return;
    }
    p->top = TRUE;
    p->next->top = FALSE;
    if (p->type != 'm')
    {
        p->next->next->top = FALSE;
    }
    set_top (world, p->next->back, locus);
    if (p->type != 'm')
    {
        set_top (world, p->next->next->back, locus);
    }
    if (p == crawlback (world->root->next))
    {
        p->back->top = FALSE;
        p->back->tyme = ROOTLENGTH;
    }
}    /* set_top */


void
set_v (node * p)
{
    if (p->type == 't')
    {
        p->v = p->length = lengthof (p);
        return;
    }
    ltov (p);
    set_v (crawlback (p->next));
    set_v (crawlback (p->next->next));
}    /* set_v */

void
ltov (node * p)
{
    p->v = lengthof (p);
}    /* ltov */

/*
 * cost matrix COST is for the ISLAND and MATRIX model all 1 with a 0
 * diagonal, this will be changable through options, but perhaps this is not
 * so important
 */
void
calc_sancost (MYREAL **cost, long numpop)
{
    long i, j;
    for (i = 0; i < numpop; i++)
    {
        for (j = 0; j < numpop; j++)
        {
            if (i != j)
                cost[i][j] = 1.0;
            else
                cost[i][i] = 0.0;
        }
    }
}

void
sankoff (world_fmt * world)
{
    MYREAL **cost;
    long i;
    cost = (MYREAL **) mymalloc (sizeof (MYREAL *) * world->numpop);
    cost[0] =
        (MYREAL *) mymalloc (sizeof (MYREAL) * world->numpop * world->numpop);
    for (i = 1; i < world->numpop; i++)
    {
        cost[i] = cost[0] + world->numpop * i;
    }
    calc_sancost (cost, world->numpop);
    santraverse (world, crawlback (world->root->next), cost, world->numpop);
    myfree(cost[0]);
    myfree(cost);
}

void
jumble (long *s, long n)
{
    long *temp, i, rr, tn = n;
	
    temp = (long *) mycalloc (1, sizeof (long) * n);
    memcpy (temp, s, sizeof (long) * n);
    for (i = 0; i < n && tn > 0; i++)
    {
        s[i] = temp[rr = RANDINT (0, tn - 1)];
        temp[rr] = temp[tn - 1];
        tn--;
    }
    myfree(temp);
}

void
santraverse (world_fmt *world, node * theNode, MYREAL **cost, long numpop)
{
    long i, ii, tie, which = 0;
    node *p = NULL, *q = NULL, *tmp = NULL, *left = NULL, *right = NULL;
    MYREAL best;
    long *poplist;
    poplist = (long *) mycalloc (1, sizeof (long) * numpop);
    if (theNode->type != 't')
    {
        if (RANDUM () > 0.5)
        {
            left = theNode->next;
            right = theNode->next->next;
        }
        else
        {
            left = theNode->next->next;
            right = theNode->next;
        }
        if (left->back != NULL)
        {
          p = crawlback (left);
	  santraverse (world, p, cost, numpop);
        }
        if (right->back != NULL)
        {
          q = crawlback (right);
	  santraverse (world, q, cost, numpop);
        }
        best = MYREAL_MAX;
        tie = 0;
        for (i = 0; i < numpop; i++)
            poplist[i] = i;
        jumble (poplist, numpop);
        for (ii = 0; ii < numpop; ii++)
        {
            i = poplist[ii];
          	if(p!=NULL && q!=NULL)
            {              
            	theNode->s[i] =
                	minimum (cost[i], p->s, numpop) + minimum (cost[i], q->s, numpop);
            }
            if (theNode->s[i] < best)
            {
                best = theNode->s[i];
                which = i;
                tie = 0;
            }
            else
            {
                if (theNode->s[i] == best)
                {
                    tie++;
                    which = i;
                }
            }
        }
        if (tie != 0)
        {
            theNode->pop = theNode->actualpop = which =
			ranbest (theNode->s, tie, best, numpop);
            theNode->s[which] -= SANKOFF_DELTA;
        }
        else
        {
            theNode->pop = theNode->actualpop = which;
        }
        if (p!=NULL && p->actualpop != which)
        {
	  //	  printf("%i> p: theNode->tyme=%f, p->tyme=%f mignode=%f\n",myID, theNode->tyme,p->tyme,p->tyme + (theNode->tyme - p->tyme) / 2.);
            tmp = add_migration (world, p, theNode->actualpop, p->actualpop,
			     (MYREAL) RANDDOUBLE(0.0, theNode->tyme - p->tyme));
            left->back = tmp;
            tmp->back = left;
            //debug FPRINTF(startfile, "%li %li\n", theNode->actualpop, p->actualpop);
        }
        if (q!=NULL && q->actualpop != which)
        {
            tmp = add_migration (world, q, theNode->actualpop, q->actualpop,
			     (MYREAL) RANDDOUBLE(0.0, theNode->tyme - q->tyme));

	    //printf("%i> q: theNode->tyme=%f, q->tyme=%f mignode=%f\n",myID, theNode->tyme,q->tyme,q->tyme + (theNode->tyme - q->tyme) / 2.);
            right->back = tmp;
            tmp->back = right;
            //debug FPRINTF(startfile, "%li %li\n", theNode->actualpop, q->actualpop);
        }
    }
    myfree(poplist);
}

/* returns minimum for sankoff routine */
MYREAL
minimum (MYREAL *vec1, MYREAL *vec2, long n)
{
    long j;
    MYREAL summ, min = MYREAL_MAX;
    for (j = 0; j < n; j++)
    {
        if (vec2[j] < MYREAL_MAX)
        {
            if ((summ = vec1[j] + vec2[j]) < min)
                min = summ;
        }
    }
    return min;
}

long
ranbest (MYREAL *array, long tie, MYREAL best, long n)
{
    long i;
    long which = RANDINT (0, tie);
    for (i = 0; i < n; i++)
    {
        if (fabs (array[i] - best) < EPSILON)
        {
            if (which == 0)
                return i;
            else
                which--;
        }
    }
    return -2;
}

///
/// allocate memory for sequence data
void
alloc_seqx (world_fmt * world, node * theNode)
{
    //    long j,z;
    long endsite = world->data->seq[0]->endsite;
    if (strchr ("u", world->options->datatype))
    {
        endsite = endsite * (world->data->seq[0]->addon + 1) + world->data->seq[0]->addon;
    }
    if (strchr ("nh", world->options->datatype))
    {
        endsite += world->data->seq[0]->addon;
    }
    allocate_xseq(&theNode->x, endsite, world->options->rcategs);
}

void     allocate_xseq(xarray_fmt *x, long sites, long categs)
{
    long j;
    (*x).s = (phenotype) mycalloc (sites, sizeof (ratelike *));
#ifdef VARMUT
    (*x).s[0] = (ratelike) mycalloc (linkedloci * categs, sizeof (MYREAL) * sites[datamodeltype] * sitelikesize[datamodeltype]);
#else
    (*x).s[0] = (ratelike) mycalloc (categs * sites, sizeof (sitelike));
#endif /*VARMUT*/
    for (j = 1; j < sites; j++)
        (*x).s[j] = (*x).s[0] + categs * j;
    
}

///
/// zeroes sequence data vector
void     zero_xseq(xarray_fmt *x, long sites, long categs)
{
    long j;
    //    (*x).s = (phenotype) mycalloc (sites, sizeof (ratelike *));
    memset((*x).s[0], 0, categs * sites * sizeof (sitelike));
    for (j = 1; j < sites; j++)
        (*x).s[j] = (*x).s[0] + categs * j;
    
}

void
make_alleles (world_fmt * world, option_fmt * options, data_fmt * data, long locus)
{
  long ii,top;
    long pop, ind;
    long zz=0;
    long zpop;
    long iu;
    char a1[DEFAULT_ALLELENMLENGTH];
    char a2[DEFAULT_ALLELENMLENGTH];
    node **nodelist = world->nodep;
    for (pop = 0; pop < data->numpop; pop++)
    {
      // preparation of combining populations
      //      if(world->custm[m2mm(pop,world->numpop)])
      //
      zpop=0;
      top = max_shuffled_individuals(options, data, pop, locus);
      for (ii = 0; ii < top; ii++)
        {
	  ind = data->shuffled[pop][locus][ii];
	  strcpy (a1, data->yy[pop][ind][locus][0]);
	  strcpy (a2, data->yy[pop][ind][locus][1]);
	  if (strcmp (a1, "?"))
            {
	      allocate_tip(world,options, &world->nodep[zz],pop,locus,zz,ind,data->indnames[pop][ind]);
	      strcat(nodelist[zz]->nayme,"A!");
	      strcat(nodelist[zz]->nayme,a1);
	      nodelist[zz++]->x.a[findAllele (data, a1, locus)] = 1.0;
	      zpop++;
            }
	  else
	    {
	      if(options->include_unknown)
		{
		  fprintf(stdout,"? found\n");
		  allocate_tip(world,options, &world->nodep[zz],pop,locus,zz,ind,data->indnames[pop][ind]);
		  for(iu=0;iu<data->maxalleles[locus]; iu++)
		    nodelist[zz]->x.a[iu] = 1.0;
		  strcat(nodelist[zz]->nayme,"A!");
		  strcat(nodelist[zz]->nayme,a1);
		  zz++;
		  zpop++;
		}                    
	    }
	  if (strcmp (a2, "?"))
            {
	      allocate_tip(world,options, &world->nodep[zz],pop,locus,zz,ind,data->indnames[pop][ind]);
	      strcat(nodelist[zz]->nayme,"B!");
	      strcat(nodelist[zz]->nayme,a2);
	      nodelist[zz++]->x.a[findAllele (data, a2, locus)] = 1.0;
	      zpop++;
            }
	  else
	    {
	      if(options->include_unknown)
		{
		  allocate_tip(world,options, &world->nodep[zz],pop,locus,zz,ind,data->indnames[pop][ind]);
		  //		  strcpy (nodelist[zz]->nayme, a1);
		  for(iu=0;iu<data->maxalleles[locus]; iu++)
		    nodelist[zz]->x.a[iu] = 1.0;
		  strcat(nodelist[zz]->nayme,"B!");
		  strcat(nodelist[zz]->nayme,a2);
		  zz++;
		  zpop++;
		}
	    }
        }
      data->numalleles[pop][locus] = zpop;
    }
    world->sumtips = zz;
    if (world->sumtips==0)
      {
	data->skiploci[locus] = TRUE;
        world->data->skiploci[locus] = TRUE;
        world->skipped += 1;
      }
}

void
find_minmax_msat_allele (world_fmt * world, data_fmt * data, long locus,
                         long *smallest, long *biggest)
{
    long pop, ind, tmp;
    *biggest = 0;
    *smallest = LONG_MAX;
    for (pop = 0; pop < data->numpop; pop++)
    {
        for (ind = 0; ind < data->numind[pop][locus]; ind++)
        {
            if (data->yy[pop][ind][locus][0][0] != '?')
            {
                if ((tmp = atoi (data->yy[pop][ind][locus][0])) > *biggest)
                    *biggest = tmp;
                if (tmp < *smallest)
                    *smallest = tmp;
            }
            if (data->yy[pop][ind][locus][1][0] != '?')
            {
                if ((tmp = atoi (data->yy[pop][ind][locus][1])) > *biggest)
                    *biggest = tmp;
                if (tmp < *smallest)
                    *smallest = tmp;
            }
        }
    }
}

void
make_microsatellites (world_fmt * world, option_fmt *options, data_fmt * data, long locus)
{
  long top;
  long pop, ind, ii;//, tmp = 0;//, tips = 0, unknownsum = 0;
  long zz=0;
  long zpop;
  long iu;
  long smallest = 0;
  long biggest = 0;
    long smax;
  node **nodelist = world->nodep;
  char a1[DEFAULT_ALLELENMLENGTH];
  char a2[DEFAULT_ALLELENMLENGTH];
  calculate_steps (world);
  find_minmax_msat_allele (world, data, locus, &smallest, &biggest);
  smax = biggest - smallest;
  if(smax == 0)
    smax = 1;
  if(world->options->micro_threshold[locus] > smax)
    world->options->micro_threshold[locus] = smax;
  smax += 2 * world->options->micro_threshold[locus];
  world->data->maxalleles[locus] = data->maxalleles[locus] = smax;
  //we start with the smallest - margin and go up to biggest + margin
  smallest = smallest - 1 - world->options->micro_threshold[locus];
  if (smallest < 0)
    smallest = 0;
  for (pop = 0; pop < data->numpop; pop++)
    {
      zpop=0;
      top = max_shuffled_individuals(options, data, pop, locus);
      for (ii = 0; ii < top; ii++)
        {
	  ind = data->shuffled[pop][locus][ii];
	  strcpy (a1, data->yy[pop][ind][locus][0]);
	  strcpy (a2, data->yy[pop][ind][locus][1]);
	  if (strcmp (a1, "?"))
            {
	      allocate_tip(world,options, &world->nodep[zz],pop,locus,zz, ind,
			   data->indnames[pop][ind]);
	      nodelist[zz]->x.a =
		(MYREAL *) myrealloc (nodelist[zz]->x.a, sizeof (MYREAL) * smax);
	      strcat(nodelist[zz]->nayme,"A!");
	      strcat(nodelist[zz]->nayme,a1);
	      nodelist[zz++]->x.a[atoi (a1) - smallest] = 1.0;
	      zpop++;
            }
	  else
	    {
	      if(options->include_unknown)
		{
		  allocate_tip(world,options, &world->nodep[zz],pop,locus,zz, ind,data->indnames[pop][ind]);
		  //		  strcpy (nodelist[zz]->nayme, a1);
		  for(iu=0;iu<data->maxalleles[locus]; iu++)
		    nodelist[zz]->x.a[iu] = 1.0;
		  strcat(nodelist[zz]->nayme,"A!");
		  strcat(nodelist[zz]->nayme,a1);
		  zz++;
		  zpop++;
		}
	    }	  
	  if (strcmp (a2, "?"))
	    {
	      allocate_tip(world,options, &world->nodep[zz],pop,locus,zz, ind,
			data->indnames[pop][ind]);
	      nodelist[zz]->x.a =
		(MYREAL *) myrealloc (nodelist[zz]->x.a, sizeof (MYREAL) * smax);
	      strcat(nodelist[zz]->nayme,"B!");
	      strcat(nodelist[zz]->nayme,a2);
	      nodelist[zz++]->x.a[atoi (a2) - smallest] = 1.0;
	      zpop++;
	    }
	  else
	    {
	      if(options->include_unknown)
		{
		  allocate_tip(world,options, &world->nodep[zz],pop,locus,zz, ind,data->indnames[pop][ind]);
		  //		  strcpy (nodelist[zz]->nayme, a1);
		  for(iu=0;iu<data->maxalleles[locus]; iu++)
		    nodelist[zz]->x.a[iu] = 1.0;
		  strcat(nodelist[zz]->nayme,"B!");
		  strcat(nodelist[zz]->nayme,a2);
		  zz++;
		  zpop++;
		}
	    }            
        }
      data->numalleles[pop][locus] = zpop;
    }
  world->sumtips = zz;
  if (world->sumtips==0)
    {
      data->skiploci[locus] = TRUE;
      world->data->skiploci[locus] = TRUE;
      world->skipped += 1;
    }
}

void
make_microbrownian (world_fmt * world, option_fmt * options, data_fmt * data, long locus)
{
  long top;
  long pop, ind, ii;
  char a1[DEFAULT_ALLELENMLENGTH];
  char a2[DEFAULT_ALLELENMLENGTH];
  node **nodelist = world->nodep;
  long zz=0;
  long zpop;
  //    long iu;
  long btotal;
  MYREAL bsum;
    for (pop = 0; pop < data->numpop; pop++)
    {
        btotal=0;
        bsum=0.;
	top = max_shuffled_individuals(options, data, pop, locus);
        for (ii = 0; ii < top; ii++)
	  {
	    ind = data->shuffled[pop][locus][ii];
            strcpy (a1, data->yy[pop][ind][locus][0]);
            strcpy (a2, data->yy[pop][ind][locus][1]);
            if (strcmp (a1, "?"))
	      {
                btotal++;
                bsum += atof (a1);
	      }
            if (strcmp (a2, "?"))
	      {
                btotal++;
                bsum += atof (a2);
	      }
	  }
        bsum /= btotal;
	
        zpop = 0;
        for (ii = 0; ii < top; ii++)
	  {
	    ind = data->shuffled[pop][locus][ii];
	    strcpy (a1, data->yy[pop][ind][locus][0]);
	    strcpy (a2, data->yy[pop][ind][locus][1]);
	    if (strcmp (a1, "?"))
	      {
		allocate_tip(world,options, &world->nodep[zz],pop,locus,zz,ind,data->indnames[pop][ind]);
		strcat(nodelist[zz]->nayme,"A!");
		strcat(nodelist[zz]->nayme,a1);
		nodelist[zz++]->x.a[0] = atof (a1);
		zpop++;
	      }
	    else
	      {
		if(options->include_unknown)
		  {
		    allocate_tip(world,options, &world->nodep[zz],pop,locus,zz,ind,data->indnames[pop][ind]);
		    strcat(nodelist[zz]->nayme,"A!");
		    strcat(nodelist[zz]->nayme,a1);
		    nodelist[zz++]->x.a[0] = bsum; //average over other values
		    // this might not be very sensible at all.
		    zpop++;
		  }
	      }
	    if (strcmp (a2, "?"))
	      {
		allocate_tip(world,options, &world->nodep[zz],pop,locus,zz,ind,data->indnames[pop][ind]);
		strcat(nodelist[zz]->nayme,"B!");
		strcat(nodelist[zz]->nayme,a2);
		nodelist[zz++]->x.a[0] = atof (a2);
		zpop++;
	      }
	    else
	      {
		if(options->include_unknown)
		  {
		    allocate_tip(world,options, &world->nodep[zz],pop,locus,zz,ind,data->indnames[pop][ind]);
		    strcat(nodelist[zz]->nayme,"B!");
		    strcat(nodelist[zz]->nayme,a2);
		    nodelist[zz++]->x.a[0] = bsum;
		    zpop++;
		  }
	      }
	  }
	data->numalleles[pop][locus] = zpop;
    }
    world->sumtips = zz;
    if (world->sumtips==0)
      {
	data->skiploci[locus] = TRUE;
        world->data->skiploci[locus] = TRUE;
        world->skipped += 1;
      }
}


/*---------------------------------------------
free_treetimes frees the ptr_array of timeslist
*/
void
free_treetimes (world_fmt * world, long size)
{
  //  long i;
  while (size >= 0)
    {
      //      for(i=world->treetimes[size].allocT-1; i>= 0 ; i--)
      //	myfree(world->treetimes[size].tl[i].lineages);
      myfree(world->treetimes[size].lineages);
      myfree(world->treetimes[size].tl);
      size--;
    }
}

///
/// construct timelist
void
construct_tymelist (world_fmt * world, timelist_fmt * timevector)
{
    long z = 0;
    long tips = 0;
    //long ii;
    //long T;
    //    long pop;
    MYREAL tmp;
    timevector->numpop = world->numpop;
    traverseNodes (world->root, &timevector, &z, world, &tips);
    timevector->T = z;
#ifdef TREEDEBUG
    printf("timevector->T=%li tips=%li\n",timevector->T,tips);
#endif
    qsort ((void *) timevector->tl, timevector->T, sizeof (vtlist), agecmp);
    if ((*timevector).tl[(*timevector).T - 1].eventnode->type != 'r')
    {
        z = 0;
        while ((*timevector).tl[z].eventnode->type != 'r')
            z++;
		
	tmp = 	(*timevector).tl[(*timevector).T - 1].eventnode->tyme + 10000.;
	(*timevector).tl[z].eventnode->tyme = tmp;
	(*timevector).tl[z].eventnode->next->tyme = tmp;
	(*timevector).tl[z].eventnode->next->next->tyme = tmp;
	
        (*timevector).tl[z].age = (*timevector).tl[z].eventnode->tyme;
        qsort ((void *) timevector->tl, timevector->T, sizeof (vtlist), agecmp);
        warning("construct_tymelist root moved: new time = %f\n",(*timevector).tl[z].eventnode->tyme);
    }
    timeslices (&timevector);
    add_partlineages (world->numpop, &timevector);
#ifdef TREEDEBUG
    fprintf(stderr,"///////////////////\n");
    long T = timevector->T-2;
    long ii, pop;
    for (ii = T; ii >= 0; ii--)
      {
	for (pop = 0; pop < world->numpop; pop++)
	  fprintf(stderr,"%li ",timevector->tl[ii].lineages[pop]);
	fprintf(stderr,"%20.10f %20.10f | %f %c %li -> %li (%li)(%li)\n",timevector->tl[ii].age,timevector->tl[ii].eventnode->tyme, timevector->tl[ii].eventnode->tyme, timevector->tl[ii].eventnode->type, timevector->tl[ii].from, timevector->tl[ii].to, timevector->tl[ii].eventnode->id, showtop(timevector->tl[ii].eventnode->back)->id);
      }
    fprintf(stderr,"/////////////////////////////////////////////\n");
#endif
}



/// find a tipdate
MYREAL 
find_tipdate(char * id, long pop, world_fmt *world)
{
  MYREAL date;
  long ind;
  long slen;
  long locus = world->locus;
  tipdate_fmt *sampledates = world->data->sampledates[pop][locus];

  for(ind=0;ind < world->data->numind[pop][locus]; ind++)
    {
      if(sampledates[ind].name != NULL)
	{
	  slen = strlen(sampledates[ind].name);
	  if(!strncmp(sampledates[ind].name,id, slen))
	    {
	      date = (sampledates[ind].date     
		      * world->options->generation_year 
		      * world->options->meanmu[locus]) * world->options->mu_rates[locus];
	      printf("%i> id='%s' name='%s' realdate=%f date=%f locus=%li pop=%li\n",myID, id,sampledates[ind].name, sampledates[ind].date, date, locus, pop);
	      fflush(stdout);
	      return date;
	    }
	}
    }
  return 0.0;
}

///
/// traverse the tree and writes node-information into the real timelist also
/// takes care that the size of timelist is increased accordingly
void
traverseNodes (node * theNode, timelist_fmt ** timevector, long *slice, world_fmt *world, long *tips)
{
    //#ifdef DEBUG
    // MYREAL tipdate = 0.0;
    //#endif
    if (theNode != NULL)
      {  
          if (theNode->type != 't')
            {
              if (theNode->next->back != NULL)
                {
                  traverseNodes (theNode->next->back, timevector, slice, world, tips);
                }
              if (theNode->type != 'm' && theNode->next->next->back != NULL)
                {
                  traverseNodes (theNode->next->next->back, timevector, slice, world, tips);
                }
              if (theNode->top)
                {
                  /*
                   * Here we are on the save side if we increase the
                   * timelist so never a fence-write can happen
                   */
                  if (*slice >= (*timevector)->allocT-1)
                    {
                      increase_timelist (timevector);
                    }
                  /*
                   * (*timevector)->tl[*slice].pop =
                   * theNode->actualpop;
                   */
                  (*timevector)->tl[*slice].age = theNode->tyme;
                  //mark the visited node this will not work with recycled node pointers
                  //does it work at all?
                  //if(theNode != (*timevector)->tl[*slice].eventnode)
                  //  (*timevector)->tl[*slice].visited = TRUE;
                  //else
                  //  (*timevector)->tl[*slice].visited = FALSE;
                  //
                  (*timevector)->tl[*slice].eventnode = theNode;
                  (*timevector)->tl[*slice].slice = *slice;
                  (*timevector)->tl[*slice].from = theNode->pop;
                  (*timevector)->tl[*slice].to = theNode->actualpop;
                  (*slice) += 1;
                }
              else
                {
                  error("traverseNodes expects to look only at TOP nodes, but received another one and died\n");
                }
            }
          else
            {
              //tipdate = find_tipdate(theNode->nayme, theNode->pop, world);
              //theNode->tyme = tipdate;
              if (*slice >= (*timevector)->allocT-1)
                {
                  increase_timelist (timevector);
                }
              (*timevector)->tl[*slice].age = theNode->tyme;
              (*timevector)->tl[*slice].eventnode = theNode;
              (*timevector)->tl[*slice].slice = *slice;
              (*timevector)->tl[*slice].from = theNode->pop;
              (*timevector)->tl[*slice].to = theNode->actualpop;
              (*tips) += 1;
              (*slice) += 1;
            }
      }
    else
        error("no node?????????");
}

void
increase_timelist (timelist_fmt ** timevector)
{
  (*timevector)->oldT = (*timevector)->allocT;
  (*timevector)->allocT += (*timevector)->allocT / 4; /* increase timelist by
						       * 25% */
  (*timevector)->tl = (vtlist *) myrealloc ((*timevector)->tl,
					    sizeof (vtlist) * ((*timevector)->allocT + 1));
  memset ((*timevector)->tl + (*timevector)->oldT, 0,
	  ((*timevector)->allocT - (*timevector)->oldT) * sizeof (vtlist));
  allocate_lineages (timevector, 0, (*timevector)->numpop);
}


void set_all_dirty (const node * root, node * p, world_fmt * world, const long locus)
{
//  node *left;
//  node *right;
  if (p->type == 'm')
    error ("MIGRATION NODE IN SMOOTH FOUND, PLEASE REPORT!!!!\n");
  if (p->type == 'i')
    {
        set_all_dirty (root, /*right=*/crawlback (p->next), world, locus);
        set_all_dirty (root, /*left=*/crawlback (p->next->next), world, locus);
        p->dirty = TRUE;
	printf("(INT %li, %li, %f)-->%li\n",p->id, p->bid, p->v, showtop(crawlback(p))->id);
    }
  else
    {
	printf("(TIP %li, %li, %f)-->%li\n",p->id, p->bid, p->v, showtop(crawlback(p))->id);
    }
}    /* set_all_dirty */

void
smooth (const node * root, node * p, world_fmt * world, const long locus)
{
 // node *left;
 // node *right;
 //   /* static */ long panic;
    /* only changed lineages are considered */
    if (!p->dirty)
        return;
	
//xcode    if (p == (crawlback (root)))
//xcode         panic = 0;
    if (p->type == 'm')
        error ("MIGRATION NODE IN SMOOTH FOUND, PLEASE REPORT!!!!\n");
    if (p->type == 'i')
    {
        smooth (root, /*right=*/crawlback (p->next), world, locus);
        smooth (root, /*left=*/crawlback (p->next->next), world, locus);
#ifdef BEAGLE
	prepare_beagle_instances(p,left, right, world->beagle);
	//	printf("(PINT %li, %li, %f)-->%li\n",p->id, p->bid, p->v, showtop(crawlback(p))->id);
#else
        (*nuview) (p, world, locus);
#endif
#ifdef UEP		
        if(world->options->uep)
            twostate_nuview (p, world, locus);
#endif		
        p->dirty = FALSE;
    }
#ifdef BEAGLE
    //else
    //  {
    //	printf("(PTIP: %li, %li, %f)-->%li\n",p->id, p->bid, p->v, showtop(crawlback(p))->id);
    //  }
#endif
}    /* smooth */


/// sets the nuview machinery so that depending on datatype the correct
/// conditional likelihood calculator is used
void
which_nuview (char datatype, boolean fastlike, boolean use_gaps, int watkins)
{
  switch (datatype)
    {
    case 'a':
      nuview = (void (*)(node *, world_fmt *, long)) nuview_allele;
      break;
    case 'b':
      nuview = (void (*)(node *, world_fmt *, long)) nuview_brownian;
      break;
    case 'm':
      switch(watkins)
	{
	case MULTISTEP:
	  prob_micro = (double (*)(MYREAL, long, world_fmt *, pair *)) prob_micro_watkins;
	  break;
	case SINGLESTEP:  
	default:
	  prob_micro = (double (*)(MYREAL, long, world_fmt *, pair *)) prob_micro_singlestep;
	}
      nuview = (void (*)(node *, world_fmt *, long)) nuview_micro;
      break;
    case 's':
      if (fastlike && !use_gaps)
	nuview = (void (*)(node *, world_fmt *, long)) nuview_sequence;
      else
	{
	  if(fastlike)
	    nuview = (void (*)(node *, world_fmt *, long)) nuview_sequence_slow;
	  else /*use_gaps*/
	    nuview = (void (*)(node *, world_fmt *, long)) nuview_sequence_slow; //nuview_gaps;	
	 }
      break;
    case 'n':
    case 'h':
    case 'u':
      nuview = (void (*)(node *, world_fmt *, long)) nuview_sequence;
      break;
    case 'f':   /* fitch, reconstruction of ancestral state
		 * method */
      nuview = (void (*)(node *, world_fmt *, long)) nuview_ancestral;
    }
}

void
nuview_allele (node * mother, world_fmt * world, const long locus)
{
  long a;
  long aa;
  long mal = world->data->maxalleles[locus];

  MYREAL freq = world->data->freq;
  //  MYREAL freqlast = world->data->freqlast; see pseudonuview
  MYREAL w1;
  MYREAL w2;
  MYREAL v1;
  MYREAL v2;
  MYREAL v1freq;
  MYREAL v2freq;
  MYREAL pija1;
  MYREAL pija2;
  MYREAL lx1;
  MYREAL lx2;

  MYREAL x3m = -MYREAL_MAX;

  MYREAL *xx1, *xx2;
  MYREAL *xx3;

  node *d1 = NULL; 
  node *d2 = NULL;

  //  MYREAL test = 0.0;
  
  children (mother, &d1, &d2);
  xx1 = d1->x.a;
  xx2 = d2->x.a;
  xx3 = mother->x.a;
  lx1 = d1->scale[0];
  lx2 = d2->scale[0];
  v1 = 1 - EXP (-d1->v);
  v2 = 1 - EXP (-d2->v);
    if (v1 >= 1.)
    {
        w1 = 0.0;
        v1 = 1.0;
    }
    else
    {
        w1 = 1.0 - v1;
    }
    if (v2 >= 1.)
    {
        w2 = 0.0;
        v2 = 1.0;
    }
    else
    {
        w2 = 1.0 - v2;
    }
    //    printf("@@nuview 1:{%f,%f, %f,%f}   2:(%f,%f, %f,%f}\n",xx1[0],xx1[1],lx1,v1,xx2[0],xx2[1],lx2,v2);

    v1freq = v1 * freq;
    v2freq = v2 * freq;
    
    for (aa = 0; aa < mal; aa++)
    {
        pija1 = pija2 = 0.0;
        for (a = 0; a < mal; a++)
        {
	  if(aa==a)
	    {
	      pija1 += (w1 + v1freq) * xx1[a];
	      pija2 += (w2 + v2freq) * xx2[a];
	    }
	  else
	    {
	      pija1 += v1 * freq * xx1[a];
	      pija2 += v2 * freq * xx2[a];
	    }
	}
        xx3[aa] = pija1 * pija2;
        //test += xx3[aa];
        if (xx3[aa] > x3m)
            x3m = xx3[aa];
		
    }
    //if (test <= 0.0)
    //    error ("xx3 is 0 or garbage!");
    for (aa = 0; aa < mal; aa++)
    {
        xx3[aa] /= x3m;
    }
    //    printf("@@finish: (%f=%f) (%f=%f)\n",xx3[0],mother->x.a[0],xx3[1],mother->x.a[1]);
    mother->scale[0] = LOG (x3m) + lx2 + lx1;
}

void
nuview_brownian (node * mother, world_fmt * world, const long locus)
{

    node *d1 = NULL, *d2 = NULL;
    MYREAL rvtot = 0.;
    MYREAL xx1, xx2, c12, diff;
    MYREAL mean1, mean2, mean, v1, v2, vtot, f1, f2;
    children (mother, &d1, &d2);
    mean1 = d1->x.a[0];
    mean2 = d2->x.a[0];
    xx1 = d1->x.a[2];
    xx2 = d2->x.a[2];
	
    v1 = d1->v + d1->x.a[1]; /* di->v == length of branch time(n1) -
		* time(n2) */
    v2 = d2->v + d2->x.a[1]; /* x.a[1] contains the deltav */
    vtot = v1 + v2;
    diff = mean1 - mean2;
    if (vtot > 0.0)
      {
	rvtot = 1./ vtot;
        f1 = v2 * rvtot;
	c12 = diff * diff * rvtot;
      }
    else
      {
        f1 = 0.5;
	c12 = HUGE;
	vtot = SMALL_VALUE;
      }
    f2 = 1.0 - f1;
    mean = f1 * mean1 + f2 * mean2;


    mother->x.a[2] =
        xx1 + xx2 + MIN (0.0, -0.5 * (LOG (vtot) + c12) + LOG2PIHALF);
    /*
     * printf("L=%f , L1=%f, L2=%f, log(vtot=%f)=%f,
     * c12=%f\n",mother->x.a[2], xx1, xx2,vtot,log(vtot),c12);
     */
    mother->x.a[1] = v1 * f1;
    mother->x.a[0] = mean;
	
}


void
nuview_micro (node * mother, world_fmt * world, const long locus)
{
    node *d1 = NULL;
    node *d2 = NULL;
    long s;
    long a;
    long aa1, aa2;
    long diff;
    long margin = world->options->micro_threshold[locus];
    MYREAL vv1, vv2, lx1, lx2;
    MYREAL x3m = -MYREAL_MAX;
    MYREAL pija1s, pija2s;
    MYREAL *xx1 = NULL, *xx2 = NULL;

    MYREAL *pm1;
    MYREAL *pm2;

    long smax = world->data->maxalleles[locus];
    MYREAL *xx3 = NULL;
    // needed for watkins' microsatellite model
    pair *helper = &world->options->msat_tuning;

    children (mother, &d1, &d2);
    vv1 = d1->v;
    vv2 = d2->v;
    xx1 = d1->x.a;
    xx2 = d2->x.a;
    xx3 = mother->x.a;
    lx1 = d1->scale[0];
    lx2 = d2->scale[0];

    pm1 = (MYREAL *) mymalloc(sizeof(MYREAL) * 4 * margin);
    pm2 = pm1 +  2 * margin;

    for (diff = 0; diff < margin; diff++)
      {
	pm1[diff] = prob_micro (vv1, diff, world, helper);
	pm2[diff] = prob_micro (vv2, diff, world, helper);
	//fprintf(stdout,"pm[diff=%li]=(%g %g)\n",diff,pm1[diff],pm2[diff]);
      }

    for (s = 0; s < smax; s++)
    {
        pija1s = pija2s = 0.0;
	aa1 = MAX (0, s - margin);
	aa2 = MIN(s + margin,smax);
	for (a = aa1; a < aa2; a++)
	  //for (a = 0; a < smax; a++)
        {
            diff = labs (s - a);
	    if(diff>=margin)
	      continue;
	    //	    fprintf(stdout,"***pm[diff=%li]=(%g %g) xx[a=%li]=(%f %f)\n",diff,pm1[diff],pm2[diff],a,xx1[a],xx2[a]);
            if (xx1[a] > 0)
            {
                pija1s += pm1[diff] * xx1[a];
            }
            if (xx2[a] > 0)
            {
                pija2s += pm2[diff] * xx2[a];
            }
        }
        xx3[s] = pija1s * pija2s;
        if (xx3[s] > x3m)
            x3m = xx3[s];
    }
    if (x3m == 0.0)
    {
        mother->scale[0] = -MYREAL_MAX;
    }
    else
    {
        for (s = 0; s < smax; s++)
        {
            xx3[s] /= x3m;
        }
        mother->scale[0] = LOG (x3m) + lx1 + lx2;
    }
    myfree(pm1);
}


void
calculate_steps (world_fmt * world)
{
    long k, diff;
    long locus = world->locus;
    const long stepnum = world->options->micro_threshold[locus];
    MYREAL ***steps = world->options->steps;
	
    for (diff = 0; diff < stepnum; diff++)
    {
        for (k = 0; k < stepnum; k++)
        {
			steps[locus][diff][k] = logfac(k + diff) + logfac (k);
        }
    }
}


///
/// calculates probability of a mutation of size diff in time t. used for msat calculations.
/// i=diff
/// prob[i_, t_] := Exp[-t] Sum[(t/2)^(i + 2 k)/((i + k)! k!), {k, 0, Infinity}]
/// calculated as Sum[Exp[(-t + ((log(t)-log(2))*i) + (2*(log(t)-log(2)) * k) - steps(i,k))
MYREAL
prob_micro_singlestep(MYREAL t, long diff, world_fmt * world, pair *helper)
{
  // helper is a dummy so that the call is the same as for prob_micro_watkins
  const long locus = world->locus;
  MYREAL *steps;
  long k;
  long k2;
  long stepnum = world->options->micro_threshold[locus];
  MYREAL temp1;
  MYREAL temp2;
  MYREAL temp3;
  MYREAL temp4;
  MYREAL const_part;
  MYREAL summ = 0.0;
  // MYREAL oldsum = 0.0;
  const MYREAL logt = LOG (t) - LOG2;
  const MYREAL logt2 = 2 * logt;
  if (diff >= stepnum)
    return summ;
  // linking to precalculated "denominator" Log[(i+k)! k!]
  // was precalculated in calculate_steps()
  steps = world->options->steps[locus][diff];
  const_part = -t + logt * diff;
  // loop unrolling to remove possible stalls
  // assumes stepnum is even
  for (k = 0; k < stepnum; k += 2)
    {
      k2 = k + 1;
      temp1 = const_part + logt2 * k - steps[k];
      temp2 = const_part + logt2 * k2 - steps[k2];
      temp3 = EXP(temp1);
      temp4 = EXP(temp2);
      summ += temp3 + temp4;
      //      if (fabs (oldsum - summ) < eps)
      //	break;
      //oldsum = summ;
    }
    return summ;
}

///
/// Calculates probability of a change of repeat number in time. Used for msat calculations.
/// based on Joe Watkins' 2007 in theoretical population biology, chapter 4.2
/// prob[diff,t,tune, upchance]
/// diff: is negative or positive depending on direction
/// t: is the time
/// tune:    0.0 with the single step mutation model
///          1.0 with the infinite allele model
///          all values in between allow multiple steps   
/// upchance: is the probability that the repeat number increases (upchance < 2/3!)
/// helper is the container for the additional variables
MYREAL
prob_micro_watkins (MYREAL t, long diff, world_fmt * world, pair *helper)
{
  const MYREAL tune = (*helper)[0];
  const MYREAL upchance = (*helper)[1]; 
  //const long locus = world->locus;
  MYREAL x;
  //  long stepnum = world->options->micro_threshold[locus];
  MYREAL summ = 0.0;
  // MYREAL oldsum = 0.0;
  const MYREAL delta = 2. * PI / 100.;
  const MYREAL expt = EXP(-t);
  MYREAL oneplustunesq;
  MYREAL tunep;
  MYREAL oneminustunet; 
  MYREAL sinx;
  MYREAL cosx;
  MYREAL invdenom;
  MYREAL first;
  MYREAL second;

  //if (diff >= stepnum)
  //  return summ;

  oneplustunesq = 1 + tune * tune ;
  tunep = (1.0 - tune) * (2.0 * upchance - 1.0) * t; 
  oneminustunet = (1.0 - tune) * t;
  for (x = -PI; x < PI; x += delta)
    {
      cosx = cos(x);
      sinx = sin(x);
      invdenom = (oneplustunesq - 2. * tune * cosx);
      if(invdenom < SMALL_VALUE)
	invdenom = 1. / SMALL_VALUE;
      first = (diff * x + tunep * sinx) * invdenom;
      second = oneminustunet * cosx * invdenom;
      summ += cos(first) * EXP(second);
    }
#ifdef DEBUG
  if(summ == 0.0)
    error("underflow in msat prob calc\n");
#endif
  return expt * INV2PI * delta * summ;
}

///
/// test implementation of Watkins' calculations
//MYREAL
//prob_micro (MYREAL t, long diff, world_fmt * world)
//{
//  return prob_micro_watkins (t, diff, 0.2 , 0.5, world);
//}
//MYINLINE MYREAL
//prob_micro (MYREAL t, long diff, world_fmt * world)
//{
//  return prob_micro_singlestep (t, diff, world);
//}


///
/// conditional likelihoods
void
nuview_sequence (node * mother, world_fmt * world, const long locus)
{
  const long endsite = world->data->seq[0]->endsite;
    long i, j, k;
    register    MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
    register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
    register MYREAL lw1, lw2, yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2, vzsumr1,
      vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
		sumy1, sumy2, ww1, ww2, zz1, zz2;
    register MYREAL freqa, freqc, freqg, freqt;//, freqr, freqy;
     register MYREAL freqar, freqcy, freqgr, freqty;
     register node *q, *r;
     register sitelike *xx1, *xx2, *xx3;
     register long rcategs = world->options->rcategs;
     register long categs = world->options->categs;
     register tbl_fmt tbl = world->tbl;
     register seqmodel_fmt *seq;
     register valrec *tbl00;
     register valrec *tblij;
     register valrec *tbljk;
    seq = world->data->seq[0];
	
    q = crawlback (mother->next);
    r = crawlback (mother->next->next);
	
    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];
    //freqr = seq->basefrequencies[NUC_R];
    //freqy = seq->basefrequencies[NUC_Y];
    freqar = seq->basefrequencies[NUC_AR];
    freqcy = seq->basefrequencies[NUC_CY];
    freqgr = seq->basefrequencies[NUC_GR];
    freqty = seq->basefrequencies[NUC_TY];
    
    lw1 = -q->v * seq->fracchange;
    
    if ((rcategs | categs) == 1)
    {
        tbl00 = tbl[0][0];
        ww1 = EXP (tbl00->ratxi * lw1);
        zz1 = EXP (tbl00->ratxv * lw1);
        ww1zz1 = ww1 * zz1;
        vv1zz1 = (1.0 - ww1) * zz1;
        lw2 = -r->v * seq->fracchange;
        ww2 = EXP (tbl00->ratxi * lw2);
        zz2 = EXP (tbl00->ratxv * lw2);
        ww2zz2 = ww2 * zz2;
        vv2zz2 = (1.0 - ww2) * zz2;
        yy1 = 1.0 - zz1;
        yy2 = 1.0 - zz2;
        for (i = 0; i < endsite; i++)
        {
            xx1 = &(q->x.s[i][0]);
            xx2 = &(r->x.s[i][0]);
            xx3 = &(mother->x.s[i][0]);
			
            xx1t0 = (*xx1)[0];
            xx1t1 = (*xx1)[1];
            xx1t2 = (*xx1)[2];
            xx1t3 = (*xx1)[3];
			
            xx2t0 = (*xx2)[0];
            xx2t1 = (*xx2)[1];
            xx2t2 = (*xx2)[2];
            xx2t3 = (*xx2)[3];
			
            sum1 = yy1 * (freqa * xx1t0 + freqc * xx1t1 + freqg * xx1t2 + freqt * xx1t3);
            sum2 = yy2 * (freqa * xx2t0 + freqc * xx2t1 + freqg * xx2t2 + freqt * xx2t3);
            
            sumr1 = freqar * xx1t0 + freqgr * xx1t2;
            sumr2 = freqar * xx2t0 + freqgr * xx2t2;
            sumy1 = freqcy * xx1t1 + freqty * xx1t3;
            sumy2 = freqcy * xx2t1 + freqty * xx2t3;
            
            vzsumr1 = vv1zz1 * sumr1;
            vzsumr2 = vv2zz2 * sumr2;
            vzsumy1 = vv1zz1 * sumy1;
            vzsumy2 = vv2zz2 * sumy2;
            (*xx3)[0] =
                (sum1 + ww1zz1 * xx1t0 + vzsumr1) * (sum2 +
													 ww2zz2 * xx2t0 +
													 vzsumr2);
            (*xx3)[1] =
                (sum1 + ww1zz1 * xx1t1 + vzsumy1) * (sum2 +
													 ww2zz2 * xx2t1 +
													 vzsumy2);
            (*xx3)[2] =
                (sum1 + ww1zz1 * xx1t2 + vzsumr1) * (sum2 +
													 ww2zz2 * xx2t2 +
													 vzsumr2);
            (*xx3)[3] =
                (sum1 + ww1zz1 * xx1t3 + vzsumy1) * (sum2 +
													 ww2zz2 * xx2t3 +
													 vzsumy2);
        }
    }
    else
    {
        for (i = 0; i < rcategs; i++)
            for (j = 0; j < categs; j++)
            {
                tblij = tbl[i][j];
                tblij->ww1 = EXP (tblij->ratxi * lw1);
                tblij->zz1 = EXP (tblij->ratxv * lw1);
                tblij->ww1zz1 = tblij->ww1 * tblij->zz1;
                tblij->vv1zz1 = (1.0 - tblij->ww1) * tblij->zz1;
            }
				lw2 = -r->v * seq->fracchange;
        for (i = 0; i < rcategs; i++)
            for (j = 0; j < categs; j++)
            {
                tblij = tbl[i][j];
                tblij->ww2 = EXP (tblij->ratxi * lw2);
                tblij->zz2 = EXP (tblij->ratxv * lw2);
                tblij->ww2zz2 = tblij->ww2 * tblij->zz2;
                tblij->vv2zz2 = (1.0 - tblij->ww2) * tblij->zz2;
            }
	for (i = 0; i < endsite; i++)
	  {
	    k = seq->category[seq->alias[i] - 1] - 1;
	    for (j = 0; j < rcategs; j++)
	      {
		tbljk = tbl[j][k];
		ww1zz1 = tbljk->ww1zz1;
		vv1zz1 = tbljk->vv1zz1;
		yy1 = 1.0 - tbljk->zz1;
		ww2zz2 = tbljk->ww2zz2;
		vv2zz2 = tbljk->vv2zz2;
		yy2 = 1.0 - tbljk->zz2;
		xx1 = &(q->x.s[i][j]);
		xx2 = &(r->x.s[i][j]);
		xx3 = &(mother->x.s[i][j]);
		
		xx1t0 = (*xx1)[0];
		xx1t1 = (*xx1)[1];
		xx1t2 = (*xx1)[2];
		xx1t3 = (*xx1)[3];
		
		xx2t0 = (*xx2)[0];
		xx2t1 = (*xx2)[1];
		xx2t2 = (*xx2)[2];
		xx2t3 = (*xx2)[3];
		
						
		sum1 =
		  yy1 * (freqa * xx1t0 + freqc * xx1t1 +
			 freqg * xx1t2 + freqt * xx1t3);
		sum2 =
		  yy2 * (freqa * xx2t0 + freqc * xx2t1 +
			 freqg * xx2t2 + freqt * xx2t3);
		sumr1 = freqar * xx1t0 + freqgr * xx1t2;
		sumr2 = freqar * xx2t0 + freqgr * xx2t2;
		sumy1 = freqcy * xx1t1 + freqty * xx1t3;
		sumy2 = freqcy * xx2t1 + freqty * xx2t3;
		vzsumr1 = vv1zz1 * sumr1;
		vzsumr2 = vv2zz2 * sumr2;
		vzsumy1 = vv1zz1 * sumy1;
		vzsumy2 = vv2zz2 * sumy2;
		(*xx3)[0] =
		  (sum1 + ww1zz1 * xx1t0 + vzsumr1) * (sum2 +
						       ww2zz2 * xx2t0 +
						       vzsumr2);
		(*xx3)[1] =
		  (sum1 + ww1zz1 * xx1t1 + vzsumy1) * (sum2 +
						       ww2zz2 * xx2t1 +
						       vzsumy2);
		(*xx3)[2] =
		  (sum1 + ww1zz1 * xx1t2 + vzsumr1) * (sum2 +
						       ww2zz2 * xx2t2 +
						       vzsumr2);
		(*xx3)[3] =
		  (sum1 + ww1zz1 * xx1t3 + vzsumy1) * (sum2 +
						       ww2zz2 * xx2t3 +
						       vzsumy2);
	      }
	  }
    }
}    /* nuview */

///
/// accurate version of conditional likleihood calculator using a scaler
#ifdef GAP
void nuview_sequence_slow (node * mother, world_fmt * world, const long locus)
{
  nuview_f84gap_slow (mother, world, locus);
}
#else
void nuview_sequence_slow (node * mother, world_fmt * world, const long locus)
{
    static long count = 0;
    const long endsite = world->data->seq[0]->endsite;
    long i, j, k;
    
    register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
    register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;

    register MYREAL lw1, lw2, yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2, vzsumr1,
        vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
        sumy1, sumy2, ww1, ww2, zz1, zz2;
    MYREAL *sxx1 = NULL, *sxx2 = NULL;
    MYREAL sxx3m, tempsxx3m, invsxx3m;
    node *q, *r;
    sitelike *xx1, *xx2, *xx3;
    long rcategs = world->options->rcategs;
    long categs = world->options->categs;
    tbl_fmt tbl = world->tbl;
    register MYREAL freqa, freqc, freqg, freqt;//, freqr, freqy;
    register MYREAL freqar, freqcy, freqgr, freqty;
    seqmodel_fmt *seq;
    valrec *tbl00;
    valrec *tblij;
    valrec *tbljk;
    seq = world->data->seq[0];

    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];
    //freqr = seq->basefrequencies[NUC_R];
    //freqy = seq->basefrequencies[NUC_Y];
    freqar = seq->basefrequencies[NUC_AR];
    freqcy = seq->basefrequencies[NUC_CY];
    freqgr = seq->basefrequencies[NUC_GR];
    freqty = seq->basefrequencies[NUC_TY];

    q = crawlback (mother->next);
    r = crawlback (mother->next->next);
    lw1 = -q->v * seq->fracchange;
    sxx1 = q->scale;
    sxx2 = r->scale;
    tbl00 = tbl[0][0];
    if ((rcategs | categs) == 1)
    {
        ww1 = EXP (tbl00->ratxi * lw1);
        zz1 = EXP (tbl00->ratxv * lw1);
        ww1zz1 = ww1 * zz1;
        vv1zz1 = (1.0 - ww1) * zz1;
        lw2 = -r->v * seq->fracchange;
        ww2 = EXP (tbl00->ratxi * lw2);
        zz2 = EXP (tbl00->ratxv * lw2);
        ww2zz2 = ww2 * zz2;
        vv2zz2 = (1.0 - ww2) * zz2;
        yy1 = 1.0 - zz1;
        yy2 = 1.0 - zz2;
        for (i = 0; i < endsite; i++)
        {
            xx1 = &(q->x.s[i][0]);
            xx2 = &(r->x.s[i][0]);
            xx3 = &(mother->x.s[i][0]);
            
            xx1t0 = (*xx1)[0];
            xx1t1 = (*xx1)[1];
            xx1t2 = (*xx1)[2];
            xx1t3 = (*xx1)[3];
            
            xx2t0 = (*xx2)[0];
            xx2t1 = (*xx2)[1];
            xx2t2 = (*xx2)[2];
            xx2t3 = (*xx2)[3];
            
            
            sum1 =
                yy1 * (freqa * xx1t0 + freqc * xx1t1 + freqg * xx1t2 +
                       freqt * xx1t3);
            sum2 =
                yy2 * (freqa * xx2t0 + freqc * xx2t1 + freqg * xx2t2 +
                       freqt * xx2t3);
            sumr1 = freqar * xx1t0 + freqgr * xx1t2;
            sumr2 = freqar * xx2t0 + freqgr * xx2t2;
            sumy1 = freqcy * xx1t1 + freqty * xx1t3;
            sumy2 = freqcy * xx2t1 + freqty * xx2t3;
            vzsumr1 = vv1zz1 * sumr1;
            vzsumr2 = vv2zz2 * sumr2;
            vzsumy1 = vv1zz1 * sumy1;
            vzsumy2 = vv2zz2 * sumy2;
            (*xx3)[0] =
                (sum1 + ww1zz1 * xx1t0 + vzsumr1) * (sum2 +
                                                     ww2zz2 * xx2t0 +
                                                     vzsumr2);
            (*xx3)[1] =
                (sum1 + ww1zz1 * xx1t1 + vzsumy1) * (sum2 +
                                                     ww2zz2 * xx2t1 +
                                                     vzsumy2);
            (*xx3)[2] =
                (sum1 + ww1zz1 * xx1t2 + vzsumr1) * (sum2 +
                                                     ww2zz2 * xx2t2 +
                                                     vzsumr2);
            (*xx3)[3] =
                (sum1 + ww1zz1 * xx1t3 + vzsumy1) * (sum2 +
                                                     ww2zz2 * xx2t3 +
                                                     vzsumy2);
            mother->scale[i] = sxx1[i] + sxx2[i];
#ifdef DEBUG
	    //	printf("[(%f,%f,%f,%f) %f]",(*xx3)[0],(*xx3)[1],(*xx3)[2],(*xx3)[3],mother->scale[i]);
#endif
        }
#ifdef DEBUG
	//   printf("id=%li\n",mother->id);
#endif
        count++;
        if (count == SCALEINTERVAL)
        {
            count = 0;
            for (i = 0; i < endsite; i++)
            {
                xx3 = &(mother->x.s[i][0]);
                sxx3m = MAX ((*xx3)[0], (*xx3)[1]);
                sxx3m = MAX (sxx3m, (*xx3)[2]);
                sxx3m = MAX (sxx3m, (*xx3)[3]);
		invsxx3m = 1.0/sxx3m;
                (*xx3)[0] *= invsxx3m;
		(*xx3)[1] *= invsxx3m;
		(*xx3)[2] *= invsxx3m;
		(*xx3)[3]  *= invsxx3m;
                mother->scale[i] += LOG (sxx3m);
            }
        }
    }
    else
    {
        for (i = 0; i < rcategs; i++)
            for (j = 0; j < categs; j++)
            {
                tblij = tbl[i][j];
                tblij->ww1 = EXP (tblij->ratxi * lw1);
                tblij->zz1 = EXP (tblij->ratxv * lw1);
                tblij->ww1zz1 = tblij->ww1 * tblij->zz1;
                tblij->vv1zz1 = (1.0 - tblij->ww1) * tblij->zz1;
            }
                lw2 = -r->v * seq->fracchange;
        for (i = 0; i < rcategs; i++)
            for (j = 0; j < categs; j++)
            {
                tblij = tbl[i][j];
                tblij->ww2 = EXP (tblij->ratxi * lw2);
                tblij->zz2 = EXP (tblij->ratxv * lw2);
                tblij->ww2zz2 = tblij->ww2 * tblij->zz2;
                tblij->vv2zz2 = (1.0 - tblij->ww2) * tblij->zz2;
            }
                for (i = 0; i < endsite; i++)
                {
                    k = seq->category[seq->alias[i] - 1] - 1;
                    for (j = 0; j < rcategs; j++)
                    {
                        tbljk = tbl[j][k];
                        ww1zz1 = tbljk->ww1zz1;
                        vv1zz1 = tbljk->vv1zz1;
                        yy1 = 1.0 - tbljk->zz1;
                        ww2zz2 = tbljk->ww2zz2;
                        vv2zz2 = tbljk->vv2zz2;
                        yy2 = 1.0 - tbljk->zz2;
                        xx1 = &(q->x.s[i][j]);
                        xx2 = &(r->x.s[i][j]);
                        xx3 = &(mother->x.s[i][j]);
                        
                        
                        xx1t0 = (*xx1)[0];
                        xx1t1 = (*xx1)[1];
                        xx1t2 = (*xx1)[2];
                        xx1t3 = (*xx1)[3];
                        
                        xx2t0 = (*xx2)[0];
                        xx2t1 = (*xx2)[1];
                        xx2t2 = (*xx2)[2];
                        xx2t3 = (*xx2)[3];
                        
                        
                        
                        sum1 =
                            yy1 * (freqa * xx1t0 + freqc * xx1t1 +
                                   freqg * xx1t2 + freqt * xx1t3);
                        sum2 =
                            yy2 * (freqa * xx2t0 + freqc * xx2t1 +
                                   freqg * xx2t2 + freqt * xx2t3);
                        sumr1 = freqar * xx1t0 + freqgr * xx1t2;
                        sumr2 = freqar * xx2t0 + freqgr * xx2t2;
                        sumy1 = freqcy * xx1t1 + freqty * xx1t3;
                        sumy2 = freqcy * xx2t1 + freqty * xx2t3;
                        vzsumr1 = vv1zz1 * sumr1;
                        vzsumr2 = vv2zz2 * sumr2;
                        vzsumy1 = vv1zz1 * sumy1;
                        vzsumy2 = vv2zz2 * sumy2;
                        (*xx3)[0] =
                            (sum1 + ww1zz1 * xx1t0 + vzsumr1) * (sum2 +
                                                                 ww2zz2 * xx2t0 +
                                                                 vzsumr2);
                        (*xx3)[1] =
                            (sum1 + ww1zz1 * xx1t1 + vzsumy1) * (sum2 +
                                                                 ww2zz2 * xx2t1 +
                                                                 vzsumy2);
                        (*xx3)[2] =
                            (sum1 + ww1zz1 * xx1t2 + vzsumr1) * (sum2 +
                                                                 ww2zz2 * xx2t2 +
                                                                 vzsumr2);
                        (*xx3)[3] =
                            (sum1 + ww1zz1 * xx1t3 + vzsumy1) * (sum2 +
                                                                 ww2zz2 * xx2t3 +
                                                                 vzsumy2);
                    }
                    mother->scale[i] = sxx1[i] + sxx2[i];
                }
                count++;
        if (count == SCALEINTERVAL)
        {
            count = 0;
            for (i = 0; i < endsite; i++)
            {
                sxx3m = -MYREAL_MAX;
                for (j = 0; j < rcategs; j++)
                {
                    xx3 = &(mother->x.s[i][j]);
                    tempsxx3m = MAX ((*xx3)[0], (*xx3)[1]);
                    tempsxx3m = MAX (tempsxx3m, (*xx3)[2]);
                    tempsxx3m = MAX (tempsxx3m, (*xx3)[3]);
                    if (tempsxx3m > sxx3m)
                        sxx3m = tempsxx3m;
                }
                for (j = 0; j < rcategs; j++)
                {
                    xx3 = &(mother->x.s[i][j]);
		    invsxx3m = 1.0 / sxx3m;
                    (*xx3)[0] *= invsxx3m, (*xx3)[1] *= invsxx3m, (*xx3)[2] *= invsxx3m,
                        (*xx3)[3] *= invsxx3m;
                }
                mother->scale[i] += LOG (sxx3m);
            }
        }
    }
}    /* nuview */
#endif /*GAP*/

///
/// conditional likleihood calcculator using a simplified scheme called ancestral
/// with two lineages the one with the higher value will win and is the ancestor
/// non-Altivec version
void
nuview_ancestral (node * mother, world_fmt * world, const long locus)
{
    long i;
    const long endsite = world->data->seq[0]->endsite;
    MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
    MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
    
    MYREAL lw1, lw2, ratio1, yy1, yy2, sum1, sum2;
    node *q, *r;
    sitelike *xx1, *xx2;
    register MYREAL freqa, freqc, freqg, freqt;
    //register MYREAL freqar, freqcy, freqgr, freqty;
    seqmodel_fmt *seq;
    seq = world->data->seq[0];

    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];
    //    freqr = seq->basefrequencies[NUC_R];
    //freqy = seq->basefrequencies[NUC_Y];
    //freqar = seq->basefrequencies[NUC_AR];
    //freqcy = seq->basefrequencies[NUC_CY];
    //freqgr = seq->basefrequencies[NUC_GR];
    //freqty = seq->basefrequencies[NUC_TY];

    q = crawlback (mother->next);
    r = crawlback (mother->next->next);
    lw1 = -q->v * seq->fracchange;
    lw2 = -r->v * seq->fracchange;
    ratio1 = lw1 / (lw1 + lw2);
    yy1 = (1. - ratio1);
    yy2 = ratio1;
    //printf("%f ", q->tyme);
    //    for (i = 0; i < endsite; i++)
    //printf("(%f %f %f %f)", q->x.s[i][0][0], q->x.s[i][0][1], q->x.s[i][0][2], q->x.s[i][0][3]);
    //printf("\n%f ", r->tyme);
    //for (i = 0; i < endsite; i++)
    //printf("(%f %f %f %f)", r->x.s[i][0][0], r->x.s[i][0][1], r->x.s[i][0][2], r->x.s[i][0][3]);
    //printf("\n");
    for (i = 0; i < endsite; i++)
      {
	xx1 = &(q->x.s[i][0]);
	xx2 = &(r->x.s[i][0]);
	
        
	xx1t0 = (*xx1)[0];
            xx1t1 = (*xx1)[1];
            xx1t2 = (*xx1)[2];
            xx1t3 = (*xx1)[3];
	    
            xx2t0 = (*xx2)[0];
            xx2t1 = (*xx2)[1];
            xx2t2 = (*xx2)[2];
            xx2t3 = (*xx2)[3];
			
            
            
            sum1 =
                yy1 * (freqa * xx1t0 + freqc * xx1t1 +
                       freqg * xx1t2 + freqt * xx1t3);
            sum2 =
                yy2 * (freqa * xx2t0 + freqc * xx2t1 +
                       freqg * xx2t2 + freqt * xx2t3);
            if (sum1 == sum2)
                sum1 += RANDUM () > 0.5 ? -1. : 1.;
            if (sum1 > sum2)
                memcpy (mother->x.s[i][0], xx1, sizeof (sitelike));
            else
                memcpy (mother->x.s[i][0], xx2, sizeof (sitelike));
        }
}    /* nuview_ancestral */

	
void
adjustroot (node * r)
{
  r->next->tyme = r->tyme;
  r->next->length = r->length;
  r->next->v = r->v;
  r->next->next->tyme = r->tyme;
  r->next->next->length = r->length;
  r->next->next->v = r->v;
}


/// \brief Calculation of conditional likelihood for new genealogy
///
/// Calculation of conditional likelihood for new genealogy
/// most of the code is unrolled for efficiency and may local variables are
/// declared also for speed
/// non-Altivec version
MYINLINE 
void
pseudonu_seq (proposal_fmt * proposal, phenotype xxx1,
	      MYREAL v1, phenotype xxx2, MYREAL v2)
{
  const long endsite = proposal->world->data->seq[0]->endsite;
  
  long i, j, k;
  
  register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
  register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
  register MYREAL ta, tc,  tg, tt;
  register MYREAL tta, ttc,  ttg, ttt;
  register MYREAL lw1, lw2, yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2, vzsumr1,
    vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
    sumy1, sumy2, ww1, ww2, zz1, zz2;
  register MYREAL freqa, freqc, freqg, freqt, freqar, freqgr, freqcy, freqty;
  register seqmodel_fmt *seq = proposal->world->data->seq[0];
  register MYREAL xxsumr1, xxsumr2, xxsumy1, xxsumy2;
  register valrec *tbl00;
  register valrec *tblij;
  register valrec *tbljk;
  register sitelike *xx1, *xx2 , *xx1copy, *xx2copy;
  register long rcategs = proposal->world->options->rcategs;
  register long categs = proposal->world->options->categs;
  register tbl_fmt tbl = proposal->world->tbl;
  freqa = seq->basefrequencies[NUC_A];
  freqc = seq->basefrequencies[NUC_C];
  freqg = seq->basefrequencies[NUC_G];
  freqt = seq->basefrequencies[NUC_T];
  //freqr = seq->basefrequencies[NUC_R];
  //freqy = seq->basefrequencies[NUC_Y];
  freqar = seq->basefrequencies[NUC_AR];
  freqcy = seq->basefrequencies[NUC_CY];
  freqgr = seq->basefrequencies[NUC_GR];
  freqty = seq->basefrequencies[NUC_TY];

  lw1 = -v1 * seq->fracchange;
  // use shortcut for dataset that do not use mutiple categories
  // if (rcategs == 1 && categs == 1)
    if ((rcategs | categs ) == 1)
    {
      tbl00 = tbl[0][0];
      lw2 = -v2 * seq->fracchange;
      
      ww1 = EXP (tbl00->ratxi * lw1);
      zz1 = EXP (tbl00->ratxv * lw1);
      ww2 = EXP (tbl00->ratxi * lw2);
      zz2 = EXP (tbl00->ratxv * lw2);
      
      ww1zz1 = ww1 * zz1;
      vv1zz1 = (1.0 - ww1) * zz1;
      
      ww2zz2 = ww2 * zz2;
      vv2zz2 = (1.0 - ww2) * zz2;
      
      yy1 = 1.0 - zz1;
      yy2 = 1.0 - zz2;
      
      for (i = 0; i < endsite; i++)
	{
	  xx1 = xx1copy = &(xxx1[i][0]);
	  xx2 = xx2copy = &(xxx2[i][0]);
	  
	  xx1t0 = (*xx1)[0];
	  ta = freqa * xx1t0; 
	  xx1t1 = (*xx1copy)[1];
	  xx2t1 = (*xx2)[1];
	  tc = freqc * xx1t1; 
	  xx1t2 = (*xx1)[2];
	  xx2t2 = (*xx2copy)[2];
	  tg = freqg * xx1t2; 
	  xx1t3 = (*xx1copy)[3];
	  xx2t3 = (*xx2)[3];
	  xx2t0 = (*xx2copy)[0];
	  tt = freqt * xx1t3;
	  
	  tta = freqa * xx2t0;
	  ttc = freqc * xx2t1;
	  ttg = freqg * xx2t2;
	  ttt = freqt * xx2t3;
	  
	  //sum1 = freqa * xx1t0 + freqc * xx1t1 + freqg * xx1t2 + freqt * xx1t3;
	  sum1 = ta + tc + tg + tt;
	  sum2 = tta + ttc + ttg + ttt ;
	  sum1 *= yy1;
	  sum2 *= yy2;
	  
	  sumr1 = freqar * xx1t0 + freqgr * xx1t2;
	  sumr2 = freqar * xx2t0 + freqgr * xx2t2;
	  sumy1 = freqcy * xx1t1 + freqty * xx1t3;
	  sumy2 = freqcy * xx2t1 + freqty * xx2t3;
	  
	  vzsumr1 = vv1zz1 * sumr1;
	  vzsumr2 = vv2zz2 * sumr2;
	  vzsumy1 = vv1zz1 * sumy1;
	  vzsumy2 = vv2zz2 * sumy2;
	  
	  xxsumr1 = sum1 + vzsumr1;
	  xxsumr2 = sum2 + vzsumr2;
	  
	  (*xx1)[0] = (xxsumr1 + ww1zz1 * xx1t0) * (xxsumr2 + ww2zz2 * xx2t0);
	  (*xx1)[2] = (xxsumr1 + ww1zz1 * xx1t2) * (xxsumr2 + ww2zz2 * xx2t2);
	  xxsumy1 = sum1 + vzsumy1;
	  xxsumy2 = sum2 + vzsumy2;
	  (*xx1)[1] = (xxsumy1 + ww1zz1 * xx1t1) * (xxsumy2 + ww2zz2 * xx2t1);
	  (*xx1)[3] = (xxsumy1 + ww1zz1 * xx1t3) * (xxsumy2 + ww2zz2 * xx2t3);
	}
    }
  else
    {
      for (i = 0; i < rcategs; i++)
	for (j = 0; j < categs; j++)
	  {
	    tblij = tbl[i][j];
	    tblij->ww1 = EXP (tblij->ratxi * lw1);
	    tblij->zz1 = EXP (tblij->ratxv * lw1);
	    tblij->ww1zz1 = tblij->ww1 * tblij->zz1;
	    tblij->vv1zz1 = (1.0 - tblij->ww1) * tblij->zz1;
	  }
      lw2 = -v2 * seq->fracchange;
      for (i = 0; i < rcategs; i++)
	for (j = 0; j < categs; j++)
	  {
	    tblij = tbl[i][j];
	    tblij->ww2 = EXP (tblij->ratxi * lw2);
	    tblij->zz2 = EXP (tblij->ratxv * lw2);
	    tblij->ww2zz2 = tblij->ww2 * tblij->zz2;
	    tblij->vv2zz2 = (1.0 - tblij->ww2) * tblij->zz2;
	  }
      for (i = 0; i < endsite; i++)
	{
	  k = seq->category[seq->alias[i] - 1] - 1;
	  for (j = 0; j < rcategs; j++)
	    {
	      tbljk = tbl[j][k];
	      ww1zz1 = tbljk->ww1zz1;
	      vv1zz1 = tbljk->vv1zz1;
	      yy1 = 1.0 - tbljk->zz1;
	      ww2zz2 = tbljk->ww2zz2;
	      vv2zz2 = tbljk->vv2zz2;
	      yy2 = 1.0 - tbljk->zz2;
	      xx1 = &(xxx1[i][j]);
	      xx2 = &(xxx2[i][j]);
	      
	      xx1t0 = (*xx1)[0];
	      xx1t1 = (*xx1)[1];
	      xx1t2 = (*xx1)[2];
	      xx1t3 = (*xx1)[3];
	      
	      xx2t0 = (*xx2)[0];
	      xx2t1 = (*xx2)[1];
	      xx2t2 = (*xx2)[2];
	      xx2t3 = (*xx2)[3];
	      
	      sum1 = yy1 * (freqa * xx1t0 + freqc * xx1t1 +
			    freqg * xx1t2 + freqt * xx1t3);
	      sum2 = yy2 * (freqa * xx2t0 + freqc * xx2t1 +
			    freqg * xx2t2 + freqt * xx2t3);
	      sumr1 = freqar * xx1t0 + freqgr * xx1t2;
	      sumr2 = freqar * xx2t0 + freqgr * xx2t2;
	      sumy1 = freqcy * xx1t1 + freqty * xx1t3;
	      sumy2 = freqcy * xx2t1 + freqty * xx2t3;
	      vzsumr1 = vv1zz1 * sumr1;
	      vzsumr2 = vv2zz2 * sumr2;
	      /* xx3[j][0] */
	      (*xx1)[0] =
		(sum1 + ww1zz1 * xx1t0 + vzsumr1) * (sum2 +
						     ww2zz2 * xx2t0 +
						     vzsumr2);
	      /* xx3[j][2] */
	      (*xx1)[2] =
		(sum1 + ww1zz1 * xx1t2 + vzsumr1) * (sum2 +
						     ww2zz2 * xx2t2 +
						     vzsumr2);
	      vzsumy1 = vv1zz1 * sumy1;
	      vzsumy2 = vv2zz2 * sumy2;
	      /* xx3[j][1] */
	      (*xx1)[1] =
		(sum1 + ww1zz1 * xx1t1 + vzsumy1) * (sum2 +
						     ww2zz2 * xx2t1 +
						     vzsumy2);
	      /* xx3[j][3] */
	      (*xx1)[3] =
		(sum1 + ww1zz1 * xx1t3 + vzsumy1) * (sum2 +
						     ww2zz2 * xx2t3 +
						     vzsumy2);
	    }
	}
    }
#ifdef BEAGLE
  printf(".");
#endif
}    /* pseudonu_seq */

			
///
/// calculates conditiional likleihood on a fake tree before acceptance/rejection scheme
/// non-altivec version
MYINLINE
    void
    pseudonu_seq_slow (proposal_fmt * proposal, phenotype xxx1, MYREAL *sxx1,
                       MYREAL v1, phenotype xxx2, MYREAL *sxx2, MYREAL v2)
{

  const long endsite = proposal->world->data->seq[0]->endsite;
  
  long i, j, k;
  static long count = 0;
  
  register MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
  register MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
  register MYREAL ta, tc,  tg, tt;
  register MYREAL tta, ttc,  ttg, ttt;

  register MYREAL lw1, lw2, yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2, vzsumr1,
    vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
    sumy1, sumy2, ww1, ww2, zz1, zz2;
  register MYREAL freqa, freqc, freqg, freqt, freqar, freqgr, freqcy, freqty;
  register seqmodel_fmt *seq = proposal->world->data->seq[0];
  register MYREAL xxsumr1, xxsumr2, xxsumy1, xxsumy2;
  register  MYREAL sxx3m, tempsxx3m, invsxx3m;
  register  sitelike *xx1, *xx2, *xx1copy, *xx2copy;
  register  valrec *tblij;
  register  valrec *tbljk;
  register  valrec *tbl00;
  register  long rcategs = proposal->world->options->rcategs;
  register  long categs = proposal->world->options->categs;
  register  tbl_fmt tbl = proposal->world->tbl;

  freqa = seq->basefrequencies[NUC_A];
  freqc = seq->basefrequencies[NUC_C];
  freqg = seq->basefrequencies[NUC_G];
  freqt = seq->basefrequencies[NUC_T];
  //freqr = seq->basefrequencies[NUC_R];
  //freqy = seq->basefrequencies[NUC_Y];
  freqar = seq->basefrequencies[NUC_AR];
  freqcy = seq->basefrequencies[NUC_CY];
  freqgr = seq->basefrequencies[NUC_GR];
  freqty = seq->basefrequencies[NUC_TY];
  lw1 = -v1 * seq->fracchange;
  if ((rcategs | categs) == 1)
    {
        tbl00 = tbl[0][0];
        lw2 = -v2 * seq->fracchange;

        ww1 = EXP (tbl00->ratxi * lw1);
        zz1 = EXP (tbl00->ratxv * lw1);
        ww2 = EXP (tbl00->ratxi * lw2);
        zz2 = EXP (tbl00->ratxv * lw2);

        ww1zz1 = ww1 * zz1;
        vv1zz1 = (1.0 - ww1) * zz1;

        ww2zz2 = ww2 * zz2;
        vv2zz2 = (1.0 - ww2) * zz2;

        yy1 = 1.0 - zz1;
        yy2 = 1.0 - zz2;

        for (i = 0; i < endsite; i++)
        {
            xx1 = xx1copy = &(xxx1[i][0]);
            xx2 = xx2copy = &(xxx2[i][0]);
 
	    xx1t0 = (*xx1)[0];
	    ta = freqa * xx1t0; 
	    xx1t1 = (*xx1copy)[1];
	    xx2t1 = (*xx2)[1];
	    tc = freqc * xx1t1; 
	    xx1t2 = (*xx1)[2];
	    xx2t2 = (*xx2copy)[2];
	    tg = freqg * xx1t2; 
	    xx1t3 = (*xx1copy)[3];
	    xx2t3 = (*xx2)[3];
	    xx2t0 = (*xx2copy)[0];
	    tt = freqt * xx1t3;
	    
	    tta = freqa * xx2t0;
	    ttc = freqc * xx2t1;
	    ttg = freqg * xx2t2;
	    ttt = freqt * xx2t3;
	  
	    sum1 = ta + tc + tg + tt;
	    sum2 = tta + ttc + ttg + ttt ;
	    sum1 *= yy1;
	    sum2 *= yy2;

            sumr1 = freqar * xx1t0 + freqgr * xx1t2;
            sumr2 = freqar * xx2t0 + freqgr * xx2t2;
            sumy1 = freqcy * xx1t1 + freqty * xx1t3;
            sumy2 = freqcy * xx2t1 + freqty * xx2t3;

            vzsumr1 = vv1zz1 * sumr1;
            vzsumr2 = vv2zz2 * sumr2;
            vzsumy1 = vv1zz1 * sumy1;
            vzsumy2 = vv2zz2 * sumy2;

	    xxsumr1 = sum1 + vzsumr1;
	    xxsumr2 = sum2 + vzsumr2;
	    
	    (*xx1)[0] = (xxsumr1 + ww1zz1 * xx1t0) * (xxsumr2 + ww2zz2 * xx2t0);
	    (*xx1)[2] = (xxsumr1 + ww1zz1 * xx1t2) * (xxsumr2 + ww2zz2 * xx2t2);
	    xxsumy1 = sum1 + vzsumy1;
	    xxsumy2 = sum2 + vzsumy2;
	    (*xx1)[1] = (xxsumy1 + ww1zz1 * xx1t1) * (xxsumy2 + ww2zz2 * xx2t1);
	    (*xx1)[3] = (xxsumy1 + ww1zz1 * xx1t3) * (xxsumy2 + ww2zz2 * xx2t3);
            sxx1[i] += sxx2[i];
        }
        count++;
        if (count == SCALEINTERVAL)
        {
            count = 0;
            for (i = 0; i < endsite; i++)
            {
                xx1 = &(xxx1[i][0]);
                sxx3m = MAX ((*xx1)[0], (*xx1)[1]);
                sxx3m = MAX (sxx3m, (*xx1)[2]);
                sxx3m = MAX (sxx3m, (*xx1)[3]);
		invsxx3m = 1. / sxx3m;
                (*xx1)[0] *= invsxx3m, (*xx1)[1] *= invsxx3m;
                (*xx1)[2] *= invsxx3m, (*xx1)[3] *= invsxx3m;
                sxx1[i] += LOG (sxx3m);
            }
        }
    }
    else
    {
        lw2 = -v2 * seq->fracchange;
        for (i = 0; i < rcategs; i++)
	  {
            for (j = 0; j < categs; j++)
	      {
                tblij = tbl[i][j];
                tblij->ww1 = EXP (tblij->ratxi * lw1);
                tblij->zz1 = EXP (tblij->ratxv * lw1);
                tblij->ww1zz1 = tblij->ww1 * tblij->zz1;
                tblij->vv1zz1 = (1.0 - tblij->ww1) * tblij->zz1;
                tblij->ww2 = EXP (tblij->ratxi * lw2);
                tblij->zz2 = EXP (tblij->ratxv * lw2);
                tblij->ww2zz2 = tblij->ww2 * tblij->zz2;
                tblij->vv2zz2 = (1.0 - tblij->ww2) * tblij->zz2;
	      }
	  }
	for (i = 0; i < endsite; i++)
	  {
	    k = seq->category[seq->alias[i] - 1] - 1;
	    for (j = 0; j < rcategs; j++)
	      {
		tbljk = tbl[j][k];
		ww1zz1 = tbljk->ww1zz1;
		vv1zz1 = tbljk->vv1zz1;
		yy1 = 1.0 - tbljk->zz1;
		ww2zz2 = tbljk->ww2zz2;
		vv2zz2 = tbljk->vv2zz2;
		yy2 = 1.0 - tbljk->zz2;
		xx1 = xx1copy = &(xxx1[i][j]);
		xx2 = xx2copy = &(xxx2[i][j]);                
                
		xx1t0 = (*xx1)[0];
		ta = freqa * xx1t0; 
		xx1t1 = (*xx1copy)[1];
		xx2t1 = (*xx2)[1];
		tc = freqc * xx1t1; 
		xx1t2 = (*xx1)[2];
		xx2t2 = (*xx2copy)[2];
		tg = freqg * xx1t2; 
		xx1t3 = (*xx1copy)[3];
		xx2t3 = (*xx2)[3];
		xx2t0 = (*xx2copy)[0];
		tt = freqt * xx1t3;
		
		tta = freqa * xx2t0;
		ttc = freqc * xx2t1;
		ttg = freqg * xx2t2;
		ttt = freqt * xx2t3;
		
		sum1 = ta + tc + tg + tt;
		sum2 = tta + ttc + ttg + ttt ;
		sum1 *= yy1;
		sum2 *= yy2;

		sumr1 = freqar * xx1t0 + freqgr * xx1t2;
		sumr2 = freqar * xx2t0 + freqgr * xx2t2;
		sumy1 = freqcy * xx1t1 + freqty * xx1t3;
		sumy2 = freqcy * xx2t1 + freqty * xx2t3;
                
		vzsumr1 = vv1zz1 * sumr1;
		vzsumr2 = vv2zz2 * sumr2;
		vzsumy1 = vv1zz1 * sumy1;
		vzsumy2 = vv2zz2 * sumy2;
                
		xxsumr1 = sum1 + vzsumr1;
		xxsumr2 = sum2 + vzsumr2;
		xxsumy1 = sum1 + vzsumy1;
		xxsumy2 = sum2 + vzsumy2;
		
		(*xx1)[0] = (xxsumr1 + ww1zz1 * xx1t0) * (xxsumr2 + ww2zz2 * xx2t0);
		(*xx1)[1] = (xxsumy1 + ww1zz1 * xx1t1) * (xxsumy2 + ww2zz2 * xx2t1);
	        (*xx1)[2] = (xxsumr1 + ww1zz1 * xx1t2) * (xxsumr2 + ww2zz2 * xx2t2);
		(*xx1)[3] = (xxsumy1 + ww1zz1 * xx1t3) * (xxsumy2 + ww2zz2 * xx2t3);
	      }
	    sxx1[i] += sxx2[i];
	  }
	count++;
        if (count == SCALEINTERVAL)
        {
            count = 0;
            for (i = 0; i < endsite; i++)
            {
	      sxx3m = -MYREAL_MAX;
	      for (j = 0; j < rcategs; j++)
                {
		  xx1 = &(xxx1[i][j]);
		  tempsxx3m = MAX ((*xx1)[0], (*xx1)[1]);
		  tempsxx3m = MAX (tempsxx3m, (*xx1)[2]);
		  tempsxx3m = MAX (tempsxx3m, (*xx1)[3]);
		  if (tempsxx3m > sxx3m)
		    sxx3m = tempsxx3m;
                }
	      invsxx3m = 1. / sxx3m;
	      for (j = 0; j < rcategs; j++)
                {
		  xx1 = &(xxx1[i][j]);
		  (*xx1)[0] *= invsxx3m, (*xx1)[1] *= invsxx3m;
		  (*xx1)[2] *= invsxx3m, (*xx1)[3] *= invsxx3m;
                }
	      sxx1[i] += LOG (sxx3m) + sxx2[i];
            }
        }
    }
}    /* pseudonu_seq */
///
/// calculates ancestral likelihood before acc/reject scheme
/// non-atlivec version
MYINLINE
void
pseudonu_anc (proposal_fmt * proposal, phenotype xxx1, MYREAL v1,
              phenotype xxx2, MYREAL v2)
{
    long i;
    const long endsite = proposal->world->data->seq[0]->endsite;
    MYREAL xx1t0, xx1t1,xx1t2, xx1t3;
    MYREAL xx2t0, xx2t1, xx2t2, xx2t3;
	
    MYREAL lw1, lw2, ratio1, yy1, yy2, sum1, sum2;
    MYREAL freqa, freqc, freqg, freqt;
    seqmodel_fmt *seq;
    sitelike *xx1, *xx2;
    seq = proposal->world->data->seq[0];
    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];
    //    freqr = seq->basefrequencies[NUC_R];
    //freqy = seq->basefrequencies[NUC_Y];
    //freqar = seq->basefrequencies[NUC_AR];
    //freqcy = seq->basefrequencies[NUC_CY];
    //freqgr = seq->basefrequencies[NUC_GR];
    //freqty = seq->basefrequencies[NUC_TY];

    lw1 = -v1 * seq->fracchange;
    lw2 = -v2 * seq->fracchange;
    ratio1 = lw1 / (lw1 + lw2);
    yy1 = (1. - ratio1);
    yy2 = ratio1;
    //for (i = 0; i < endsite; i++)
    //printf("(%f %f %f %f)", xxx1[i][0][0], xxx1[i][0][1], xxx1[i][0][2], xxx1[i][0][3]);
    //printf("\n");
    //for (i = 0; i < endsite; i++)
    //printf("(%f %f %f %f)", xxx2[i][0][0], xxx2[i][0][1], xxx2[i][0][2], xxx2[i][0][3]);
    //printf("\n");
    for (i = 0; i < endsite; i++)
    {
        xx1 = &(xxx1[i][0]);
        xx2 = &(xxx2[i][0]);
		
        
		xx1t0 = (*xx1)[0];
		xx1t1 = (*xx1)[1];
		xx1t2 = (*xx1)[2];
		xx1t3 = (*xx1)[3];
		
		xx2t0 = (*xx2)[0];
		xx2t1 = (*xx2)[1];
		xx2t2 = (*xx2)[2];
		xx2t3 = (*xx2)[3];
		
		
        sum1 =
            yy1 * (freqa * xx1t0 + freqc * xx1t1 + freqg * xx1t2 +
                   freqt * xx1t3);
        sum2 =
            yy2 * (freqa * xx2t0 + freqc * xx2t1 + freqg * xx2t2 +
                   freqt * xx2t3);
        if (sum1 == sum2)
            sum1 += RANDUM () > 0.5 ? -1. : 1.;
        if (sum1 > sum2)
            memcpy (xxx1[i][0], *xx1, sizeof (sitelike));
        else
            memcpy (xxx1[i][0], *xx2, sizeof (sitelike));
    }
}    /* pseudo_nu_anc */

///
/// ancestral conditional method, calculates ancestral tree likelihood
MYINLINE MYREAL
pseudo_tl_anc (phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
               proposal_fmt * proposal, world_fmt * world)
{
  const long endsite = proposal->world->data->seq[0]->endsite;
    contribarr tterm;
    MYREAL summ;
    long i;
    sitelike *x1;
	
    register MYREAL freqa, freqc, freqg, freqt;
    seqmodel_fmt *seq;
    seq = proposal->world->data->seq[0];

    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];

    seq = world->data->seq[0];
    summ = 0.0;
    for (i = 0; i < endsite; i++)
    {
        x1 = &(xx1[i][0]);
        tterm[0] =
            freqa * (*x1)[0] + freqc * (*x1)[1] + freqg * (*x1)[2] +
            freqt * (*x1)[3];
        summ += seq->aliasweight[i] * LOG (tterm[0]);
        //printf("pseudo %3li> %f %f \n", i, tterm[0], summ);
    }
    return summ;
} /*anc*/


///
/// calculates the conditional likelihood on a tree
MYINLINE MYREAL
pseudo_tl_seq (phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
               proposal_fmt * proposal, world_fmt * world)
{
  const long endsite = proposal->world->data->seq[0]->endsite;
    contribarr tterm;
    contribarr clai;
    contribarr like;
    contribarr nulike;
    //long          size = sizeof(MYREAL) * world->options->rcategs;
    MYREAL summ    = 0.0;
    MYREAL sum2    = 0.0;
    MYREAL sumc    = 0.0;
    MYREAL sumterm = 0.0;
    MYREAL lterm;
    long i; 
    long j; 
    long k;
    long lai;
    worldoption_fmt *opt;
    register MYREAL freqa, freqc, freqg, freqt;//, freqr, freqy;
    seqmodel_fmt *seq = proposal->world->data->seq[0];
    sitelike *x1;

    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];
    //freqr = seq->basefrequencies[NUC_R];
    //freqy = seq->basefrequencies[NUC_Y];
	
    opt = world->options;
    summ = 0.0;
    if (opt->rcategs == 1 && opt->categs == 1)
    {
        for (i = 0; i < endsite; i++)
        {
            x1 = &(xx1[i][0]);
            tterm[0] = 
                freqa * (*x1)[0] + freqc * (*x1)[1] + freqg * (*x1)[2] +
                freqt * (*x1)[3];
            summ += seq->aliasweight[i] * (LOG (tterm[0]) + proposal->mf[i]);
        }
    }
    else
    {
        for (i = 0; i < endsite; i++)
        {
            for (j = 0; j < opt->rcategs; j++)
            {
                x1 = &(xx1[i][j]);
                tterm[j] = freqa * (*x1)[0] + freqc * (*x1)[1] + freqg * (*x1)[2] + freqt * (*x1)[3];
            }
            sumterm = 0.0;
            for (j = 0; j < opt->rcategs; j++)
                sumterm += opt->probcat[j] * tterm[j];
            lterm = LOG (sumterm) + proposal->mf[i];
            for (j = 0; j < opt->rcategs; j++)
                clai[j] = tterm[j] / sumterm;
            swap (clai, world->contribution[i]);
            //memcpy(world->contribution[i], clai, size);
            summ += seq->aliasweight[i] * lterm;
        }
        for (j = 0; j < opt->rcategs; j++)
            like[j] = 1.0;
        for (i = 0; i < seq->sites[world->locus]; i++)
        {
            sumc = 0.0;
            for (k = 0; k < opt->rcategs; k++)
                sumc += opt->probcat[k] * like[k];
            sumc *= opt->lambda;
            if ((seq->ally[i] > 0) && (seq->location[seq->ally[i] - 1] > 0))
            {
                lai = seq->location[seq->ally[i] - 1];
                swap (world->contribution[lai - 1], clai);
                //memcpy(clai, world->contribution[lai - 1], size);
                for (j = 0; j < opt->rcategs; j++)
                    nulike[j] = ((1.0 - opt->lambda) * like[j] + sumc) * clai[j];
            }
            else
            {
                for (j = 0; j < opt->rcategs; j++)
                    nulike[j] = ((1.0 - opt->lambda) * like[j] + sumc);
            }
            swap (nulike, like);
            //memcpy(like, nulike, size);
        }
        sum2 = 0.0;
        for (i = 0; i < opt->rcategs; i++)
            sum2 += opt->probcat[i] * like[i];
        summ += LOG (sum2);
    }
    //    printf("summ=%f\n",summ);
    return summ;
}

MYINLINE 
MYREAL
pseudo_tl_snp (phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
               proposal_fmt * proposal, world_fmt * world)
{
  const long endsite = world->data->seq[0]->endsite;
    contribarr tterm, invariants;
    contribarr like, nulike, clai;
    //long          size = sizeof(MYREAL) * world->options->rcategs;
    MYREAL summ, sum2, sumc, sumterm, lterm;
    long i, j, k, lai;
    worldoption_fmt *opt = world->options;
    register MYREAL freqa, freqc, freqg, freqt;//, freqr, freqy;
    seqmodel_fmt *seq = world->data->seq[0];
    sitelike *x1;
    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];
    //freqr = seq->basefrequencies[NUC_R];
    //freqy = seq->basefrequencies[NUC_Y];

    summ = 0.0;
    memset (invariants, 0, sizeof (contribarr));
    snp_invariants (invariants, world, world->locus, xx1);
    if ((opt->rcategs | opt->categs) == 1)
    {
        for (i = 0; i < endsite - seq->addon; i++)
        {
            x1 = &(xx1[i][0]);
            tterm[0] =
                (freqa * (*x1)[0] + freqc * (*x1)[1] +
                 freqg * (*x1)[2] + freqt * (*x1)[3]) / invariants[0];
            if (tterm[0] == 0.0)
	      {
		fprintf(stderr,"Freq/condL=(%f,%f | %f,%f | %f,%f | %f,%f) Invars=%f\n", 
			freqa , (*x1)[0] , freqc , (*x1)[1] ,
			freqg , (*x1)[2] , freqt , (*x1)[3] , invariants[0]);
       
                warning ("Tree incompatible with data\n");
		return -HUGE;
	      }
            lterm = LOG (tterm[0]) + proposal->mf[i];
            summ += seq->aliasweight[i] * lterm;
        }
        like[0] = 1.0;
        for (i = 0; i < seq->sites[world->locus]; i++)
        {
            sumc = opt->lambda * like[0];
            nulike[0] = ((1.0 - opt->lambda) * like[0] + sumc);
            //memcpy(like, nulike, size);
            swap (nulike, like);
        }
        summ += LOG (like[0]);
        return summ;
    }
    else
    {
        for (i = 0; i < endsite - 4; i++)
        {
           // k = seq->category[seq->alias[i] - 1] - 1;
            for (j = 0; j < opt->rcategs; j++)
            {
                x1 = &(xx1[i][j]);
                tterm[j] =
                    (freqa * (*x1)[0] + freqc * (*x1)[1] +
                     freqg * (*x1)[2] +
                     freqt * (*x1)[3]) / invariants[j];
            }
            sumterm = 0.0;
            for (j = 0; j < opt->rcategs; j++)
                sumterm += opt->probcat[j] * tterm[j];
            lterm = LOG (sumterm) + proposal->mf[i];
            for (j = 0; j < opt->rcategs; j++)
                clai[j] = tterm[j] / sumterm;
            swap (clai, world->contribution[i]);
            summ += seq->aliasweight[i] * lterm;
        }
        for (j = 0; j < opt->rcategs; j++)
            like[j] = 1.0;
        for (i = 0; i < seq->sites[world->locus]; i++)
        {
            sumc = 0.0;
            for (k = 0; k < opt->rcategs; k++)
                sumc += opt->probcat[k] * like[k];
            sumc *= opt->lambda;
            if ((seq->ally[i] > 0) && (seq->location[seq->ally[i] - 1] > 0))
            {
                lai = seq->location[seq->ally[i] - 1];
                swap (world->contribution[lai - 1], clai);
                for (j = 0; j < opt->rcategs; j++)
                    nulike[j] = ((1.0 - opt->lambda) * like[j] + sumc) * clai[j];
            }
            else
            {
                for (j = 0; j < opt->rcategs; j++)
                    nulike[j] = ((1.0 - opt->lambda) * like[j] + sumc);
            }
            swap (nulike, like);
            //memcpy(like, nulike, size);
        }
        sum2 = 0.0;
        for (i = 0; i < opt->rcategs; i++)
            sum2 += opt->probcat[i] * like[i];
        summ += LOG (sum2);
        return summ;
    }
}

MYINLINE
MYREAL
pseudo_tl_snp_unlinked (phenotype xx1, phenotype xx2, MYREAL v1, MYREAL v2,
                        proposal_fmt * proposal, world_fmt * world)
{
  const long endsite = world->data->seq[0]->endsite;
    contribarr tterm, invariants;
    MYREAL summ, datasum = 0, lterm, result = 0;
    long i;
    //worldoption_fmt *opt;
    MYREAL freqa, freqc, freqg, freqt;//, freqr, freqy;
    seqmodel_fmt *seq = proposal->world->data->seq[0];
    sitelike *x1;
    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];
    //freqr = seq->basefrequencies[NUC_R];
   // freqy = seq->basefrequencies[NUC_Y];
	
    //opt = world->options;
    seq = world->data->seq[0];
    summ = 0.0;
    memset (invariants, 0, sizeof (contribarr));
    snp_invariants (invariants, world, world->locus,  xx1);
    for (i = 0; i < endsite; i++)
    {
        x1 = &(xx1[i][0]);
        tterm[0] =
            (freqa * (*x1)[0] + freqc * (*x1)[1] +
             freqg * (*x1)[2] + freqt * (*x1)[3]);
        if (tterm[0] == 0.0)
            error ("Tree incompatible with data\n");
		
        if (i % 5 == 0)
        {
            lterm = LOG (tterm[0]) + proposal->mf[i];
            summ = 0;
            datasum = seq->aliasweight[i / 5] * lterm;
        }
        else
            summ += pow (tterm[0], (MYREAL) seq->aliasweight[i / 5]);
        if (((i + 1) % 5) == 0 && i != 0)
            result +=
                datasum + LOG ((1 - EXP (LOG (summ) - datasum)) / invariants[0]);
    }
    //EXP (sum) is the prob(xa | g)
    //              EXP (datasum) is prob(? a | g)
    //              panelsum = invariants is prob(x ? |g)
    // (datasum - sum) / invariants
    // ++some small number business
    return result;
}




MYREAL
treelike_anc (world_fmt * world, long locus)
{
  const long endsite = world->data->seq[0]->endsite;
    contribarr tterm;
	
    MYREAL summ;
    long i;
    node *p;
    seqmodel_fmt *seq = world->data->seq[0];
    register MYREAL freqa, freqc,freqg,freqt;
    sitelike *x1;
    freqa = seq->basefrequencies[NUC_A];
    freqc = seq->basefrequencies[NUC_C];
    freqg = seq->basefrequencies[NUC_G];
    freqt = seq->basefrequencies[NUC_T];

    p = crawlback (world->root->next);
    summ = 0.0;
    for (i = 0; i < endsite; i++)
    {
        x1 = &(p->x.s[i][0]);
        tterm[0] =
            freqa * (*x1)[0] + freqc * (*x1)[1] +
            freqg * (*x1)[2] + freqt * (*x1)[3];
        summ += seq->aliasweight[i] * LOG (tterm[0]);
        //printf("real  %3li> %f %f \n", i, tterm[0], summ);
    }
    return summ;
}    /* treelike_anc */


///
/// write a tree to a diskfile, compile setting will allow to do this as a NEXUS file
void
treeout (FILE * file, node * joint, node * p, long s)
{
    /* write out file with representation of final tree */
    MYREAL x;
    char migstring[30];
    if (p->type == 't')
    {
#ifdef UEP	
        if(p->uep!=NULL)
            FPRINTF (file, "%s [%c]", p->nayme, p->uep[0]);
        else
            FPRINTF (file, "%s ", p->nayme);
#else
        FPRINTF (file, "%s ", p->nayme);
#endif
#ifdef TREECOMMENTS
	FPRINTF(file,"[& t:%.10f]",p->tyme);
#endif
    }
    else
    {
        FPRINTF (file, "(");
        treeout (file, joint, crawlback (p->next), s);
        FPRINTF (file, ",");
        treeout (file, joint, crawlback (p->next->next), s);
        FPRINTF (file, ")");
#ifdef TREECOMMENTS
	FPRINTF(file,"[& c:%.10f]",p->tyme);
#endif
    }
    if (p != joint)
    {
        x = crawlback (p)->tyme - p->tyme;
#ifdef UEP	
	if(p->uep!=NULL)
	  FPRINTF (file, ":%.10f [%c]", x, p->uep[0]);
	else
	  FPRINTF (file, ":%.10f ",  x);
#else
	FPRINTF (file, ":%.10f ", x);
#endif
        p = showtop (p->back);
        while (p->type == 'm')
	  {
            sprintf (migstring, " [&M %li %li:%g]", p->pop, p->actualpop,
                     p->tyme - showtop (p->next->back)->tyme);
            FPRINTF (file, "%s", migstring);
            p = showtop (p->back);
	  }
    }
    else
    {
#ifdef NEXUSTREE
	FPRINTF (file, ";\n");
#else
	FPRINTF (file, ":0;\n");
#endif
    }
}    /* treeout */

///
/// write a tree to a diskfile, compile setting will allow to do this as a NEXUS file
void
debugtreeout (FILE * file, node * joint, node * p, long s)
{
    /* write out file with representation of final tree */
  //  MYREAL x;
    char migstring[30];
    if (p->type == 't')
    {
        FPRINTF (file, "%li ", p->id);
    }
    else
    {
        FPRINTF (file, "(");
        debugtreeout (file, joint, crawlback (p->next), s);
        FPRINTF (file, ",");
        debugtreeout (file, joint, crawlback (p->next->next), s);
        FPRINTF (file, ")<%li>",p->id);
    }
    if (p != joint)
    {
        //x = crawlback (p)->tyme - p->tyme;
        //	FPRINTF (file, ":%.10f ", x);
        p = showtop (p->back);
        while (p->type == 'm')
	  {
            sprintf (migstring, " [%li-%li]", p->next->id, p->id);
            FPRINTF (file, "%s", migstring);
            p = showtop (p->back);
	  }
    }
    else
    {
	FPRINTF (file, ":0\n");

    }
}    /* treeout */

///
/// writes tree in newick format to a string
void
treeout_string (char ** file, long *filesize, long *pos, node * joint, node * p, long s)
{
  char *tmp;
  MYREAL x;
  char migstring[30];
  tmp = (char *) mycalloc(1024,sizeof(char));
  /* write out file with representation of final tree */
  //    long w;
  if (p->type == 't')
    {
#ifdef UEP
        if(p->uep!=NULL)
	  print_to_buffer(file, filesize, tmp, pos, "%s [%c]", p->nayme, p->uep[0]);
        else
	  print_to_buffer(file, filesize, tmp, pos, "%s ", p->nayme);
#else
        print_to_buffer(file, filesize, tmp, pos, "%s ", p->nayme);
#endif		
#ifdef TREECOMMENTS
	print_to_buffer(file, filesize, tmp, pos,"[& t:%.10f]",p->tyme);
#endif
    }
    else
    {
      print_to_buffer(file, filesize, tmp, pos, "(");
      treeout_string (file, filesize, pos, joint, crawlback (p->next), s);
      print_to_buffer(file, filesize, tmp, pos, ",");
      treeout_string (file, filesize, pos, joint, crawlback (p->next->next), s);
      print_to_buffer(file, filesize, tmp, pos, ")");
#ifdef TREECOMMENTS
      print_to_buffer(file, filesize, tmp, pos,"[& c:%.10f]",p->tyme);
#endif

    }
    if (p == joint)
    {
        x = 0.0;
    }
    else
    {
        x = crawlback (p)->tyme - p->tyme;
    }
#ifdef UEP
    if(p->uep!=NULL)
      print_to_buffer(file, filesize, tmp, pos, ":%.10f [%c]", x, p->uep[0]);
    else
      print_to_buffer(file, filesize, tmp, pos, ":%.10f ", x);
#else
    print_to_buffer(file, filesize, tmp, pos, ":%.10f ", x);
#endif
    if (p != joint)
    {
        p = showtop (p->back);
        while (p->type == 'm')
        {
            sprintf (migstring, " [&M %li %li:%g]", p->pop, p->actualpop,
                     p->tyme - showtop (p->next->back)->tyme);
            print_to_buffer(file, filesize, tmp, pos, "%s", migstring);
            p = showtop (p->back);
        }
    }
    else
    {
      print_to_buffer(file, filesize, tmp, pos, ";\n\0");
    }
    myfree(tmp);
}    /* treeout_string */


void
print_tree (world_fmt * world, long g, long *filepos)
{
#ifdef NEXUSTREE
  static long counter = 0;
#endif
  static long count   = 1;
  long pos            = 0;
  long allocval       = 0;
  char *tmp;
  tmp = (char *) calloc(1024,sizeof(char));

  switch (world->options->treeprint)
    {
    case BEST:
      if(world->options->treeinmemory)
	{
	  *filepos = 0;
	  pos = 0;
	  allocval = world->treespacealloc[world->locus];
	  //	  if (world->likelihood[g] >= world->besttreelike[world->locus])
	  if (world->param_like >= world->besttreelike[world->locus])
	    {
	      // printf("%i> print_tree(BEST) at %p  locus=%li bestlike=%f like=%f\n",myID, &world->treespace[world->locus], world->locus,world->besttreelike[world->locus], world->likelihood[g]);
	      world->besttreelike[world->locus] = world->param_like;
	      // world->besttreelike[world->locus] = world->likelihood[g];
	      world->treespace[world->locus] = (char *) myrealloc(world->treespace[world->locus],
								  sizeof(char) * LONGLINESIZE);
	      allocval = LONGLINESIZE;
	      //speed problem	      memset(world->treespace[world->locus],0, sizeof(char) * LONGLINESIZE);
	      world->treespace[world->locus][0]='\0';
	      pos = sprintf (world->treespace[world->locus], "\n[& Locus %li, best ln(L) = %f (c=coalescent node, t=tipnode) ]\n",
			     world->locus + 1, world->likelihood[g]);
	      treeout_string (&(world->treespace[world->locus]), 
			      &allocval,&pos, 
			      crawlback (world->root->next),
			      crawlback (world->root->next), 0);
	      world->treespacealloc[world->locus] = allocval;
	      world->treespacenum[world->locus] = pos;
	    }
	  //	  else
	  //  {
	  //    printf("%i> *** print_tree(BEST) locus=%li bestlike=%f like=%f\n", myID, world->locus, world->besttreelike[world->locus], world->likelihood[g]);
	  //  }
	}
      else
	{
	  if (world->likelihood[g] > world->allikemax)
	    {
	      if (world->allikemax == -MYREAL_MAX)
		{
		  *filepos = ftell (world->treefile);
		}
	      else
		{
		  fseek (world->treefile, *filepos, SEEK_SET);
		}
	      world->allikemax = world->likelihood[g];
	      FPRINTF (world->treefile,
			     "\n[& Locus %li, best ln(L) = %f (c=coalescent node, t=tipnode) ]\n",
		       world->locus + 1, world->likelihood[g]);
	      treeout (world->treefile, crawlback (world->root->next),
		       crawlback (world->root->next), 0);
	    }
	}
      break;
    case ALL:
    case LASTCHAIN:
      if((count++ % world->options->treeinc) == 0)
	{
	  if (world->in_last_chain)
	    {
	      if(world->options->treeinmemory)
		{
		  *filepos = 0;
		  pos = world->treespacenum[world->locus];
		  allocval = world->treespacealloc[world->locus];
#ifdef NEXUSTREE
		  pos = print_to_buffer(&(world->treespace[world->locus]), &allocval, tmp, &pos,
					"\n[& Locus %li, best ln(L) = %f (c=coalescent node, t=tipnode) ]\ntree repl.%li = [&R] ",
					world->locus + 1, world->likelihood[g], counter++);
#else
		  pos = print_to_buffer(&(world->treespace[world->locus]), &allocval, tmp, &pos,
					"\n[& Locus %li, best ln(L) = %f (c=coalescent node, t=tipnode) ]\n",
					world->locus + 1, world->likelihood[g]);
#endif
		  treeout_string (&(world->treespace[world->locus]), 
				  &allocval,&pos, 
				  crawlback (world->root->next),
				  crawlback (world->root->next), 0);
		  world->treespacealloc[world->locus] = allocval;
		  world->treespacenum[world->locus] = pos;
		}
	      else
		{
		  // writing direct to file
		  FPRINTF (world->treefile,
			   "\n[& Locus %li, ln(L) = %f (c=coalescent node, t=tipnode) ]\n",
			   world->locus + 1, world->likelihood[g]);
#ifdef NEXUSTREE
		  FPRINTF (world->treefile,"tree repl.%li = [&R] ",counter++);
#endif
		  treeout (world->treefile, crawlback (world->root->next),
			   crawlback (world->root->next), 0);
		}
	    }
	}
      break;
    case _NONE:
      break;
    default:
      break;
    }
  myfree(tmp);
}


boolean
treereader (world_fmt * world,  option_fmt *options, data_fmt * data)
{
    /*
     * read a  tree with or without migration events from the usertree and set up nodes
     * and pointers
     */
    boolean has_migration=FALSE;
    node **nodelist;
    char *nayme;
    char *temp, *temp2;
    long pop, w, zz, z = 0, zzz = 0;
    world->nodep = (node **) mycalloc (world->sumtips+(world->sumtips+1),sizeof (node *));
    temp = (char *) mymalloc (LINESIZE * sizeof (char));
    temp2 = (char *) mymalloc (LINESIZE * sizeof (char));
    has_migration = treeread (data->utreefile, world, options, &(world->root), NULL);
    length_to_times (world->root->next->back);
    nodelist = (node **) mycalloc (1, sizeof (node *) * (world->sumtips + 1));
    pop = -1;
    set_tree_pop (world->root, &pop);
    allocate_x (world->root, world, world->options->datatype, WITHTIPS);
    find_tips (world->root, nodelist, &z);
    for (pop = 0; pop < data->numpop; pop++)
    {
        for (w = 0; w < data->numind[pop][world->locus]; w++)
        {
            strcpy (temp2, data->indnames[pop][w][world->locus]);
	    unpad(temp2," ");
	    //            temp2[strcspn (temp2, " ")] = '\0';
            sprintf (temp, "%li%s", pop, temp2);
            for (zz = 0; zz < z; zz++)
            {
                nayme = nodelist[zz]->nayme;
		unpad(nayme," _");
                if (!strcmp (temp, nayme) || !strcmp (temp2, nayme))
                {
                    world->nodep[zzz++] = nodelist[zz];
                    break;
                }
            }
        }
    }
    myfree(nodelist);
    myfree(temp);
    myfree(temp2);
    return has_migration;
}


char
processlength (FILE * file, node ** p, option_fmt *options)
{
    char ch;
    long digit, ordzero;
    MYREAL valyew, divisor;
    boolean pointread, minusread;
	
    ordzero = '0';
    pointread = FALSE;
    minusread = FALSE;
    valyew = 0.0;
    divisor = 1.0;
    ch = getc (file);
    digit = ch - ordzero;
    while (((unsigned long) digit <= 9) | (ch == '.') || (ch == '-'))
    {
        if (ch == '.')
            pointread = TRUE;
        else if (ch == '-')
            minusread = TRUE;
        else
        {
            valyew = valyew * 10.0 + digit;
            if (pointread)
                divisor *= 10.0;
        }
        ch = getc (file);
        digit = ch - ordzero;
    }
    if (!minusread)
        (*p)->length = valyew / divisor;
    else
        (*p)->length = 0.0;
    return ch;
}

boolean
treeread (FILE * file, world_fmt * world, option_fmt *options, node ** pp, node * q)
{
    node *p=NULL;
    boolean has_migration=FALSE;
    char ch = getc (file);
    while (ch != ';')
    {
        switch (ch)
        {
			case '(':
				p = create_interior_node (world, &q);
				q = p->next;
				ch = getc (file);
				break;
			case ',':
				q = q->next;
				if (q->top)
				{
					usererror ("Multifurcation handling not yet installed");
				}
					ch = getc (file);
				break;
			case ')':
				p = showtop (q);
				q = p->back;
				ch = getc (file);
				break;
			case ' ':
			case '\n':
			case '\t':
				ch = getc (file);
				break;
			case ':':
				ch = processlength (file, &p, options);
				break;
			case '[':
			  has_migration = processbracket (file, world, &p, &ch);
			  if(has_migration>0)
			    {
			      q->back = p;
			      p->back = q;
			    }
			  break;
			default:
				p = create_tip_node (file, world, options, &q, &ch);
				break;
        }
    }
    if(p!=NULL)
      {
	    p->length = 10000.;
    	(*pp) = showtop (p->back);
      }
    fscanf (file, "%*[^\n]");
    getc (file);
    return has_migration;
}

void
length_to_times (node * p)
{
    node *q;
    if (p->type != 't')
    {
        length_to_times ((p)->next->back);
        if ((p)->type == 'i')
            length_to_times ((p)->next->next->back);
    }
    q = showtop ((p)->back);
    q->tyme = q->next->tyme = q->next->next->tyme = (p)->tyme + (p)->length;
}

void
find_tips (node * p, node ** nodelist, long *z)
{
    if (p->type == 't')
    {
        nodelist[(*z)++] = p;
    }
    else
    {
        if (p->next->back != NULL)
            find_tips (crawlback (p->next), nodelist, z);
        if (p->next->next->back != NULL)
            find_tips (crawlback (p->next->next), nodelist, z);
    }
}

long
find_firstpop (node * p)
{
    static boolean found = FALSE;
    static long pop = -1;
    if (p->type == 'm')
    {
        found = TRUE;
        pop = p->pop;
    }
    else
    {
        if (p->next->back != NULL)
        {
            find_firstpop (p->next->back);
            if (found)
                return pop;
        }
        if (p->next->next->back != NULL)
            find_firstpop (p->next->next->back);
    }
    return pop;
}

/* touches only coalescent nodes! migration nodes are already set */
void
set_tree_pop_old (node * p, long *pop)
{
    if (p->type != 'r')
    {
		
        (*pop) =
		(showtop (p->back)->actualpop !=
		 *pop) ? showtop (p->back)->actualpop : *pop;
    }
    p->actualpop = p->pop = *pop;
    if (p->type != 't')
    {
        if (p->next->back != NULL)
        {
            set_tree_pop (crawlback (p->next), pop);
        }
        if (p->type != 'm' && p->next->next->back != NULL)
        {
            set_tree_pop (crawlback (p->next->next), pop);
        }
    }
}

/// sets the actualpop and pop values deducting from the mgiration nodes that
/// are already set
void
set_tree_pop (node * p, long *pop)
{
  static boolean done=FALSE;
  switch(p->type)
    {
    case 'r':
	  set_tree_pop(p->next->back, pop);
	  *pop = p->actualpop = p->pop = p->next->back->pop;
      if(!done)
	{
	  done=TRUE;
	  if(*pop == -1)
	    *pop = 0;
	  set_tree_pop(p,pop);
	}
      break;
    case 't':
      if (*pop != -1)
	p->actualpop = p->pop = *pop;
      break;
    case 'i':
      set_tree_pop(p->next->back,pop);
      if(*pop != -1)
	p->actualpop = p->pop = p->next->back->pop;
      set_tree_pop(p->next->next->back,pop);
      if(*pop != -1)
	{
	  p->actualpop = p->pop = p->next->next->back->pop;
	  *pop = p->pop;
	}
      break;
    case 'm':
      *pop = p->actualpop;
      set_tree_pop(p->next->back, pop);
      *pop = p->pop;
      break;
    }
}


node *
create_interior_node (world_fmt * world, node ** q)
{
    node *p;
    p = allocate_nodelet (world, 3, 'i');
    p->top = TRUE;
    p->scale = (MYREAL *) mycalloc (world->data->seq[0]->endsite+world->data->seq[0]->addon, sizeof (MYREAL));
    p->s = (MYREAL *) mycalloc (1, sizeof (MYREAL) * world->numpop);
    p->back = *q;
    if ((*q) == NULL)
      create_root_node (world, &p);
    else
        (*q)->back = p;
    return p;
}

node *
create_root_node (world_fmt *world, node ** q)
{
    node *p;
    p = allocate_nodelet (world, 3, 'r');
    p->top = TRUE;
    p->next->back = *q;
    (*q)->back = p->next;
    return p;
}


node *
create_tip_node (FILE * file, world_fmt * world, option_fmt *options, node ** q, char *ch)
{
    long pop;
    node *p;
    char c;
    char *nayme;
    long i = 1;
    nayme = (char *) mycalloc (options->nmlength+1, sizeof (char));
    nayme[0] = (*ch);
    while (strchr ("[):;,\t\n\r", (int) (c = getc (file))) == NULL)
        nayme[i++] = c;
    nayme[i] = '\0';
    p = allocate_nodelet (world, 1, 't');
    p->nayme = (char *) mycalloc (options->nmlength+5, sizeof (char));
    p->top = TRUE;
    p->tip = TRUE;
    // PB is this needed with usertrees
    p->scale = (MYREAL *) mycalloc (world->data->seq[0]->endsite+world->data->seq[0]->addon, sizeof (MYREAL));
    p->s = (MYREAL *) mycalloc (world->numpop, sizeof (MYREAL)  );
    for (i = 0; i < world->numpop; i++)
      {
        p->s[i] = MYREAL_MAX;    
      }

    strncpy (p->nayme, nayme,options->nmlength);
    unpad(p->nayme," ");
    translate(p->nayme,' ', '_');
    sscanf(nayme,"%li ",&pop);		       
    p->s[pop] = 0;
    p->back = *q;
    (*q)->back = p;
    myfree(nayme);
    (*ch) = c;
    return p;
}

boolean
processbracket (FILE * file, world_fmt *world, node ** p, char * ch)
{
    boolean migfound=FALSE;
    long pop1, pop2;
    MYREAL utime;
    char c;
    c = getc (file);
    if (c == '&')
    {
        c = getc (file);
        switch (c)
        {
			case 'M':
#ifdef USE_MYREAL_FLOAT
				fscanf (file, "%li %li:%f", &pop1, &pop2, &utime);
#else
				fscanf (file, "%li %li:%lf", &pop1, &pop2, &utime);
#endif
				/*c=*/getc (file);
				(*p) = add_migration (world, *p, pop1, pop2, utime);
				migfound = TRUE;
				break;
			default:
				while (c != ']')
					c = getc (file);
				break;
        }
    }
    else
    {
        while (c != ']')
            c = getc (file);
    }
    *ch = getc (file);
    return migfound;
}


node *
add_migration (world_fmt *world, node * p, long from, long to, MYREAL utime)
{
    node *tmp;
    tmp = allocate_nodelet (world, 2, 'm');
    tmp->top = TRUE;
    tmp->next->back = p;
    p->back = tmp->next;
    tmp->length = p->length - utime;
    p->length = utime;
    tmp->tyme = p->tyme + utime;
    tmp->pop = tmp->next->pop = from;
    tmp->actualpop = tmp->next->actualpop = to;
    return tmp;
}

void
allocate_x (node * p, world_fmt * world, char datatype, boolean withtips)
{
    if (p->type != 't')
    {
        if (p->next->back != NULL)
            allocate_x (crawlback (p->next), world, datatype, withtips);
        if (p->next->next->back != NULL)
            allocate_x (crawlback (p->next->next), world, datatype, withtips);
        if (strchr (SEQUENCETYPES, world->options->datatype))
        {
            alloc_seqx (world, p);
        }
        else
        {
            if (strchr (SEQUENCETYPES, world->options->datatype))
                p->x.a = (MYREAL *) mycalloc (1, world->sumtips * sizeof (MYREAL));
            else
                p->x.a =
                    (MYREAL *) mycalloc (1,
                                         MAX (world->sumtips,
                                              world->data->maxalleles[world->locus]) *
                                         sizeof (MYREAL));
        }
#ifdef UEP
        if(world->options->uep)
        {
            p->uep = (int *) mycalloc (world->data->uepsites, sizeof (int));
            p->ux.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        }
#endif
		
    }
    else
    {
        if (withtips)
        {
            if (strchr (SEQUENCETYPES, world->options->datatype))
            {
                alloc_seqx (world, p);
            }
            else
            {
                p->x.a =
				(MYREAL *) mycalloc (1,
                                     world->data->maxalleles[world->locus] *
                                     sizeof (MYREAL));
            }
        }
        //#ifdef UEP
        //      if(world->options->uep)
        // {
        //   p->uep = (long *) mycalloc (world->data->uepsites, sizeof (long));
        //   p->ux.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        // }
        //#endif
		
    }
}

long
number_genomes (int datatype)
{
    switch (datatype)
    {
    case 'a':
    case 'b':
    case 'm':
      return 2;
    case 's':
    case 'n':
    case 'h':
    case 'u':
    case 'f':
      return 1;
    default:
      error ("Wrong data type");
    }
    return 0;
	
}


void
copy_tree (world_fmt * original, world_fmt * kopie)
{
    kopie->root = copy_node (original, original->root, kopie, NULL);
}

node *
copy_node (world_fmt * original, node * o, world_fmt * kopie, node * last)
{
    static long z = 0;
	
    node *t = NULL, *t2, *t3;
    if (o == NULL)
        return NULL;
    if (!o->top)
        error ("copy_tree messed up");
	
    switch (o->type)
    {
		case 'r':
			z = 0;
		case 'i':
			t = (node *) mycalloc (1, sizeof (node));
			t2 = (node *) mycalloc (1, sizeof (node));
			t3 = (node *) mycalloc (1, sizeof (node));
			t->next = t2;
			t2->next = t3;
			t3->next = t;
			copy_node_content (original, kopie, o, t);
			copy_node_content (original, kopie, o->next, t2);
			copy_node_content (original, kopie, o->next->next, t3);
			if (o->next->back != NULL)
				t2->back = copy_node (original, o->next->back, kopie, t2);
				if (o->next->next->back != NULL)
					t3->back = copy_node (original, o->next->next->back, kopie, t3);
					t->back = last;
			break;
		case 'm':
			t = (node *) mycalloc (1, sizeof (node));
			t2 = (node *) mycalloc (1, sizeof (node));
			t->next = t2;
			t2->next = t;
			copy_node_content (original, kopie, o, t);
			copy_node_content (original, kopie, o->next, t2);
			t2->back = copy_node (original, o->next->back, kopie, t2);
			t->back = last;
			break;
		case 't':
			t = (node *) mycalloc (1, sizeof (node));
			//kopie->nodep[z++] = t;
			t->next = t;
			copy_node_content (original, kopie, o, t);
			t->back = last;
			break;
    }
    return t;
}

///
/// copy the guts of a node in the tree (but to where?)
void
copy_node_content (world_fmt * original, world_fmt * kopie, node * o,
                   node * t)
{
    long i;
//    long j;
    const long endsite = original->data->seq[0]->endsite;
    const long rcategs = original->options->rcategs;
    t->type = o->type;
    t->number = o->number;
    t->pop = o->pop;
    t->actualpop = o->actualpop;
    t->id = o->id;
    t->top = o->top;
    t->dirty = o->dirty;
    t->v = o->v;
    t->tyme = o->tyme;
    t->length = o->length;
	
    if (t->top && t->type != 'm')
    {
        t->scale = (MYREAL *) mycalloc (endsite, sizeof (MYREAL));
        memcpy (t->scale, o->scale, sizeof (MYREAL) * endsite);
#ifdef UEP
		
        if (original->options->uep)
        {
            t->uep = (int *) mycalloc (original->data->uepsites, sizeof (int));
            memcpy (t->uep, o->uep, sizeof (long) * original->data->uepsites);
            t->ux.s = (pair *) mycalloc (original->data->uepsites, sizeof (pair));
            for (i = 0; i < original->data->uepsites; i++)
            {
                memcpy (t->ux.s[i], o->ux.s[i], sizeof (pair));
            }
        }
#endif
        if (strchr (SEQUENCETYPES, original->options->datatype))
        {
            alloc_seqx (kopie, t);
            memcpy (t->scale, o->scale, sizeof (MYREAL) * endsite);
            for (i = 0; i < endsite; i++)
            {
//                for (j = 0; j < rcategs; j++)
//                {
                    memcpy (t->x.s[i], o->x.s[i], rcategs * sizeof (sitelike));

//                }
            }
        }
        else
        {
            t->x.a =
			(MYREAL *) mycalloc (1,
                                 original->data->maxalleles[original->locus] *
                                 sizeof (MYREAL));
            memcpy (t->x.a, o->x.a,
                    sizeof (MYREAL) *
                    original->data->maxalleles[original->locus]);
        }
    }
    //if (o->s != NULL)
    ///{
    //t->s = (MYREAL *) mycalloc(1, sizeof(MYREAL) * original->numpop);
    //memcpy] (t->s, o->s, sizeof(MYREAL) * original->numpop);
    //
    //      }
    //
    else
        t->s = NULL;
    if (o->nayme != NULL)
    {
        t->nayme = (char *) mycalloc (20, sizeof (char));
        strncpy (t->nayme, o->nayme, 10);
    }
    else
        t->nayme = NULL;
}

void
swap_tree (world_fmt * tthis, world_fmt * tthat)
{
    node *tmp;
    node **nodetemps;
    tmp = tthis->root;
    tthis->root = tthat->root;
    tthat->root = tmp;
    nodetemps = tthis->nodep;
    tthis->nodep = tthat->nodep;
    tthat->nodep = nodetemps;

}

void
calc_treelength (node * p, MYREAL *treelen)
{
    node *pn, *pnn;
    switch (p->type)
    {
		case 't':
			break;
		case 'm':
			error ("yelp\n");
			break;
		case 'i':
			pn = crawlback (p->next);
			calc_treelength (pn, treelen);
			pnn = crawlback (p->next->next);
			calc_treelength (pnn, treelen);
			break;
		default:
        error ("default reached");
    }
    pn = showtop (crawlback (p));
    if (pn->type != 'r')
        *treelen += pn->tyme - p->tyme;
}

MYREAL
calc_pseudotreelength (proposal_fmt * proposal, MYREAL treelen)
{
    MYREAL len = 0.0;
    MYREAL ot = proposal->origin->tyme;
    MYREAL obt = proposal->oback->tyme;
    MYREAL tt = proposal->target->tyme;
    MYREAL rt = proposal->world->root->next->back->tyme;
    //target is not root
    if (proposal->target != proposal->world->root)
    {
        //oback is not root
        if (proposal->oback != proposal->world->root)
        {
            len = treelen - (obt - ot) + (proposal->time - ot);
            //printf("pseudo_treelen: ob!=r t!=r %f\n", len);
        }
        else
        {
            //oback is root
            len = treelen - (obt - ot) - (rt - tt) +
			(proposal->time - ot) + (proposal->time - tt);
            //printf("pseudo_treelen: ob=r t!=r %f\n", len);
        }
    }
    else
        //target is root
    {
        //oback is not root
        if (proposal->oback != proposal->world->root)
        {
            len = treelen - (obt - ot) + (proposal->time - ot) - (rt - tt) +
			(proposal->time - tt);
            //printf("pseudo_treelen: ob!=r t=r %f\n", len);
        }
        else
        {
            //oback is root
            len = treelen - (obt - ot) - (obt - tt) + (proposal->time - ot)
			+ (proposal->time - tt);
            //printf("pseudo_treelen: ob=r t=r %f\n", len);
        }
    }
    return len;
}

void
swap (void *a, void *b)
{
    void *t;
    t = a;
    a = b;
    b = t;
}


void
free_tree (node * p, world_fmt * world)
{
    if (p != NULL)
    {
        if (p->type != 't')
        {
            if (p->next->back != NULL)
            {
                free_tree (p->next->back, world);
            }
            if (p->type != 'm' && p->next->next->back != NULL)
            {
                free_tree (p->next->next->back, world);
            }
        }
        switch (p->type)
        {
	case 'm':
	  free_mignodelet (p, world);
	  break;
	case 't':
	  myfree(p->nayme);
	  free_tipnodelet (p, world);
	  break;
	case 'i':
	  free_nodelet (p, 3, world);
	  break;
	case 'r':
	  free_nodelet (p, 3, world);
	  world->root = NULL;
	  break;
	default:
	  error("error in freeing nodes, a node with no type found");
	  break;
        }
    }
}

///
/// delete tipnode-nodelet and its content data
void
free_tipnodelet (node * p, world_fmt * world)
{
  free_nodedata (p, world);
  collect_nodelet(world, p);
  //myfree(p); 
  //p->id = -p->id;
}

///
/// delete migratenode-nodelets
void
free_mignodelet (node * p, world_fmt * world)
{
  node *q = p->next;
  collect_nodelet(world, q);
  collect_nodelet(world, p);
  //  myfree(q);
  //myfree(p);
}

///
/// delete interior node and its data
void
free_nodelet (node * p, long num, world_fmt * world)
{
    long i;
    node *q , *r;
    switch(num)
      {
      case 1:
	if(p->top == TRUE)
	  free_nodedata (p, world);
	collect_nodelet(world, p);
	//myfree(p);
	break;
      case 2:	
	q = p->next;
	if(p->top == TRUE)
	  free_nodedata (p, world);
	collect_nodelet(world, p);
	//	myfree(p);
	if(q!=NULL)
	  {
	    if(q->top == TRUE)
	      free_nodedata (q, world);
	    collect_nodelet(world, q);
	    // myfree(q);
	  }
	break;
      case 3:	
	q = p->next;
	r = p->next->next;
	if(p->top == TRUE)
	  free_nodedata (p, world);
	collect_nodelet(world, p);
	//	myfree(p);
	//p->id = -p->id;
	if(q!=NULL)
	  {
	    if(q->top == TRUE)
	      free_nodedata (q, world);
	    collect_nodelet(world, q);
	    //	    myfree(q);
	    //q->id = -q->id;
	  }
	if(r!=NULL)
	  {
	    if(r->top == TRUE)
	      free_nodedata (r, world);
	    collect_nodelet(world, r);
	    //	    myfree(r); 
	    //r->id = -r->id;
	  }
	break;
      default:
	for (i = 0; i < num; i++)
	  {
	    q = p->next;
	    if(p->top == TRUE)
	      free_nodedata (p, world);
	    //	    myfree(p);
	    collect_nodelet(world, p);
	    p = q;
	  }
      }
}

void
free_nodedata (node * p, world_fmt * world)
{
  //    long endsite;
    if (strchr (SEQUENCETYPES, world->options->datatype))
    {
      //       endsite = world->data->seq->endsite;
        myfree(p->x.s[0]);
	myfree(p->x.s);
    }
    else
    {
        myfree(p->x.a);
    }
    if (p->s != NULL)
      myfree(p->s);
    if(p->scale!=NULL)
      myfree(p->scale);
    
#ifdef UEP
    if (world->options->uep)
        myfree(p->uep);
#endif
	
}

///
/// calculates the 1/probability that no event happens
/// is safeguarded against a problem on two tip trees 
/// that will return a result of zero and so a inverse of INF
/// INF is replaced by MYREAL_MAX
MYREAL
inverse_logprob_noevent (world_fmt * world, long interval)
{
    long pop, k;
    MYREAL result = 0.0;
    for (pop = 0; pop < world->numpop; pop++)
    {
        k = world->treetimes[0].tl[interval].lineages[pop];
        result +=
            (k * (k - 1) / world->param0[pop]) + sum_migprob (world, pop, interval);
    }
    if(result > 0.0)
      return 1./result;
    else
      {
	if(world->treetimes[0].tl[interval].eventnode->type=='t')
	  {
	    return world->treetimes[0].tl[interval].age;
	  }
	else
	  return MYREAL_MAX;
      }
}

/// calculates the sum of all immigrations into a specific population
MYREAL
sum_migprob (world_fmt * world, long pop, long interval)
{
    long i;
    MYREAL result = 0.0;
    long *lineages = world->treetimes[0].tl[interval].lineages;
    long msta = world->mstart[pop];
    long msto = world->mend[pop];
    for (i = msta; i < msto; i++)
    {
        result += world->param0[i];
    }
    return result * lineages[pop];
}

void debugline(node *up)
{
  node *thenode = up;
  node * down = crawlback(showtop(up));
  while (thenode != down)
    {
      printf("<%li:%li %c>\n",thenode->id, showtop(thenode)->id, thenode->type);
      thenode = showtop(thenode)->back;
    }
}
