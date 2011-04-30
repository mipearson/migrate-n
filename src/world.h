#ifndef WORLD_INCLUDE
#define WORLD_INCLUDE
/** \file world.h
*/
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 W O R L D   R O U T I N E S 
 
 creates tree structures,
 calculates smple parameter estimates (FST,...)
 reads tree [has to be done]
 
 prints results,
 and finally helps to destroy itself.
                                                                                                               
 Peter Beerli, started 1996
 beerli@fsu.edu
 
Copyright 1996-2003 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2006 Peter Beerli, Tallahassee FL

This software is distributed free of charge for non-commercial use
and is copyrighted. Of course, we do not guarantee that the software
works and are not responsible for any damage you may cause or have.

$Id: world.h 1515 2009-05-07 01:02:31Z beerli $
-------------------------------------------------------*/

#include "migration.h"
extern void create_world (world_fmt ** world, long loci);
extern void init_world (world_fmt * world, data_fmt * data,
                            option_fmt * options);
extern void increase_timearchive (world_fmt * world, long locus, long sample,
                                      long numpop, long rep);

extern void calc_simple_param (world_fmt * world, data_fmt * data);
extern void set_bounds (long *increment, long *steps, long *chains,
                            const option_fmt * options, const char type);
extern void print_menu_locus (FILE *file, world_fmt * world, long locus);
extern void print_menu_chain (char type, long chain, long steps,
                                  world_fmt * world, option_fmt * options,
                                  long rep);
extern void print_menu_coalnodes (FILE * file, world_fmt * world, long G,
                                      long rep);
extern void print_progress(worldoption_fmt * options, world_fmt * world,
               long rep, long visited, long accepted);
extern void burnin_chain (world_fmt * world);
extern void print_finish (world_fmt * world, long filepos);
extern void copy_time (world_fmt * world, timelist_fmt * ltl, long from,
                           long to, long np, long rep);
extern void create_plot (world_fmt * world, char ***plane,
                                  nr_fmt * nr, long Gloci, boolean multilocus);
extern void print_simresults (world_fmt * world);
extern void print_list (world_fmt ** universe, option_fmt * options,
                            data_fmt * data);
extern void free_universe (world_fmt ** worlds, long numworlds, option_fmt *options);
extern void test_locus_like (nr_fmt * nr, MYREAL *param0, MYREAL *param1,
                                 long df, long locus, world_fmt * world,
                                 long *maxwhich, long maxnum, boolean withhead,
                                 char *this_string);
extern void test_loci_like (nr_fmt * nr, MYREAL *param0, MYREAL *param1,
                                long df, long loci, world_fmt * world,
                                long *maxwhich, long maxnum, boolean withhead,
                                char *this_string);


extern void precalc_world (world_fmt * world);
extern void reprecalc_world (world_fmt * world, long that);
extern void klone (world_fmt * original, world_fmt * kopie,
                       option_fmt * options, data_fmt * data, MYREAL temperature);
extern void klone_part (world_fmt * original, world_fmt * kopie,
                            option_fmt * options, data_fmt * data, MYREAL temperature);
extern void clone_polish (world_fmt * original, world_fmt * kopie);
extern long chance_swap_tree (world_fmt * tthis, world_fmt * that);
extern void advance_clone_like (world_fmt * world, long accepted, long *j);
extern void polish_world (world_fmt * world);
extern void print_alpha_curve (world_fmt * world, timearchive_fmt ** atl, long *gmaxptr);
extern void print_cov (world_fmt * world, long numpop, long loci, MYREAL ***cov);
extern void print_mighist (world_fmt * world);
extern void print_gelmanr (MYREAL average, MYREAL biggest);
extern void prognose_time (char *nowstr, world_fmt * world,
                           long options_increment, long steps, char *spacer, boolean tobuffer);
extern boolean is_same(long i, long j);
#endif /*WORLD_INCLUDE */
