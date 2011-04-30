/*! \file=tree.h */
#ifndef TREE_INCLUDE
#define TREE_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 T R E E B U I L D I N G   R O U T I N E S 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 updated 2009

Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2009 Peter Beerli, Tallahassee FL
 
$Id: tree.h 1833 2011-03-20 19:09:41Z beerli $
-------------------------------------------------------*/
#include "migration.h"
void buildtree (world_fmt * world, option_fmt * options, data_fmt * data, long locus);
void create_treetimelist (world_fmt * world, timelist_fmt ** ltl, long locus);
void allocate_lineages (timelist_fmt **timevector, const long offset, const long numpop);
void allocate_tip (world_fmt * world, option_fmt * options, node ** p, long pop, long locus, long a, long ind, char **tipnames);
void fix_times (world_fmt * world, option_fmt * options);
void first_smooth (world_fmt * world, long locus);
void smooth (const node * root, node * p, world_fmt * world, const long locus);
void set_all_dirty (const node * root, node * p, world_fmt * world, const long locus);
void set_dirty (node * p);
void construct_tymelist (world_fmt * world, timelist_fmt * timevector);
//void timeslices (timelist_fmt ** timevector);
void add_partlineages (long numpop, timelist_fmt ** timevector);
MYREAL treelikelihood (world_fmt * world);
MYREAL pseudotreelikelihood (world_fmt * world,proposal_fmt * proposal);
void set_pop (node * theNode, long pop, long actualpop);
void pseudonuview (proposal_fmt * proposal, xarray_fmt xx1, MYREAL *lx1, MYREAL v1, xarray_fmt xx2, MYREAL *lx2, MYREAL v2);
void set_v (node * p);
void ltov (node * p);
//void treeout (FILE * treefile, node * joint, node * p, long s);
void print_tree (world_fmt * world, long g, long *filepos);
MYREAL find_tipdate(char * id, long pop, world_fmt *world);
void allocatetips (world_fmt * world, option_fmt * options, data_fmt * data, long locus);
//void allocate_x (node * p, world_fmt * world, char datatype, boolean withtips);
void     allocate_xseq(xarray_fmt *x, long sites, long categs);
//void     zero_xseq(xarray_fmt *x, long sites, long categs);
//void copy_tree (world_fmt * original, world_fmt * kopie);
void swap_tree (world_fmt * tthis, world_fmt * tthat);
//void free_nodelet (node * p, long num, world_fmt * world);
void free_mignodelet (node * p, world_fmt * world);
void calc_treelength (node * p, MYREAL *treelen);
//MYREAL calc_pseudotreelength (proposal_fmt * proposal, MYREAL treelen);
//void swap (void *a, void *b);
//void free_treetimes (world_fmt * world, long size);
void free_tree (node * p, world_fmt * world);
//long number_genomes (int datatype);
node *add_migration (world_fmt *world, node * p, long from, long to, MYREAL utime);
#endif /*TREE_INCLUDE */
