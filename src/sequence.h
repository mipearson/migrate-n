#ifndef SEQUENCE_INCLUDE
#define SEQUENCE_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 S E Q U E N C E S   R O U T I N E S 
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: sequence.h 1348 2008-09-11 14:42:13Z beerli $
-------------------------------------------------------*/

extern void initratio (option_fmt * options);
extern void initfreqs (MYREAL *freqa, MYREAL *freqc, MYREAL *freqg,
                           MYREAL *freqt);
extern void initcatn (long *categs);
extern boolean initcategs (long categs, MYREAL *rate, MYREAL *probcat);
extern void initprobcat (long categs, MYREAL *probsum, MYREAL *probcat);
extern void constrain_rates(long categs, MYREAL *rate, MYREAL *probcat);
extern void initlambda (option_fmt * options);

extern void init_sequences (world_fmt * world, option_fmt * options,
			    data_fmt * data, long locus);
extern void init_sequences2 (world_fmt * world, seqmodel_fmt * seq,
			     long locus);
extern void init_tbl (world_fmt * world, long locus);
extern void print_weights (FILE * outfile, world_fmt * world,
			   option_fmt * options, long locus);
extern void print_tbl (FILE * outfile, world_fmt * world,
                           option_fmt * options, long locus);
extern void print_seqfreqs (FILE * outfile, world_fmt * world,
			    option_fmt * options);
extern MYREAL treelike_seq (world_fmt * world, long locus);
extern MYREAL treelike_snp (world_fmt * world, long locus);
extern void snp_invariants (contribarr invariants, world_fmt *world, long locus,  phenotype x1);
extern void make_sequences (world_fmt * world, option_fmt * options,
			    data_fmt * data, long locus);
extern void make_invarsites (world_fmt * world, data_fmt * data, long locus);
extern void make_invarsites_unlinked (world_fmt * world, data_fmt * data,
				      long locus);
extern void make_snp (world_fmt * world, option_fmt * options,
		      data_fmt * data, long locus);
extern MYREAL treelike_snp_unlinked (world_fmt * world, long locus);
extern void copy_seq (world_fmt * original, world_fmt * kopie);
extern void init_sequences_aliases (world_fmt * world, option_fmt * options,
                                        data_fmt * data, long locus);
extern void find_rates_fromdata(data_fmt * data, option_fmt * options);
extern void free_seq(seqmodel_fmt **seq, long seqnum);

#endif /*SEQUENCE_INCLUDE */
