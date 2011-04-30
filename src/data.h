#ifndef DATA_INCLUDE
#define DATA_INCLUDE
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
 
 
 Theta(1)=4 N(1)mu, Theta(2)=4 N(2)mu,
 M(1) = m(1)/mu, and M(2)= m(2)/mu
                                                                                                               
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: data.h 1593 2009-12-06 22:42:27Z beerli $
-------------------------------------------------------*/

#include "migration.h"

void create_data (data_fmt ** data);
void get_data (FILE * infile, data_fmt * data, option_fmt * options, world_fmt *world);
void print_data (world_fmt * world, option_fmt * options,
                     data_fmt * data);
void print_spectra(world_fmt * world, option_fmt * options,data_fmt * data);
void print_data_summary (FILE * file, world_fmt * world,
                             option_fmt * options, data_fmt * data);
long find_missing(data_fmt *data, long pop, long locus);
short findAllele (data_fmt * data, char s[], long locus);
void free_datapart (data_fmt * data, option_fmt * options, long locus);
void read_distance_fromfile (FILE * dfile, long tips, long nmlength,
                                 MYREAL **m);
void read_geofile (data_fmt * data, option_fmt * options, long numpop);

void init_data_structure1 (data_fmt ** data, option_fmt * options);
void init_data_structure2 (data_fmt ** data, option_fmt * options,
                               long pop);
void init_data_structure3 (data_fmt * data);
void create_alleles (data_fmt * data, option_fmt * options);
void set_numind (data_fmt * data);

long max_shuffled_individuals(option_fmt *options, data_fmt *data, long pop, long locus);

long number_genomes (int type);
void destroy_data(data_fmt * data);

#endif /*DATA_INCLUDE */
