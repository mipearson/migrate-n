// Bayes update scheme
#ifndef _BAYESUPDATE_
#define _BAYESUPDATE_
//
// started November 2000
// (c) Peter Beerli Tallahassee 2000-2006
// $Id: bayes.h 1737 2010-10-21 16:10:07Z beerli $
//
#include "migration.h"
extern boolean acceptBayes (MYREAL newval, MYREAL oldval);
extern long bayes_update (world_fmt *world);
extern void bayes_free(world_fmt *world);
extern void bayes_fill(world_fmt *world, option_fmt *options);
extern void bayes_init(bayes_fmt *bayes, world_fmt *world, option_fmt *options);
extern void bayes_save(world_fmt *world, long step);
extern void bayes_stat(world_fmt *world);
extern long setup_bayes_map(longpair *map, char *custm2, long numpop, long numpop2, long size);
extern void bayes_init_histogram(world_fmt * world, option_fmt * options);
extern void adjust_bayes_bins(world_fmt * world, long locus);
extern void calculate_credibility_interval(world_fmt * world, long locus);
extern void bayes_reset(world_fmt * world);
extern void bayes_check_and_fix_param(world_fmt *world, option_fmt *options);
#ifdef ZNZ
extern void print_bayes_mdimfileheader(znzFile file, long interval, world_fmt *world, data_fmt *data);
#else
extern void print_bayes_mdimfileheader(FILE *file, long interval, world_fmt *world, data_fmt *data);
#endif
extern void bayes_print_accept(FILE * outfile, world_fmt * world);
extern MYREAL probg_treetimes(world_fmt *world);
extern void recalc_timelist (world_fmt * world, MYREAL new_ratio, MYREAL old_ratio);
extern void bayes_smooth(MYREAL *x, long xelem, long el, boolean lastfirst);
extern void bayes_set_param(MYREAL *param, MYREAL newparam, long which, char *custm2, long numpop);
extern MYREAL log_prior_ratio_uni(MYREAL newparam, 
			   MYREAL oldparam, 
			   bayes_fmt * bayes, 
				  long which);
extern MYREAL log_prior_ratio_all(world_fmt *world, MYREAL *newvals);
extern MYREAL calculate_prior(world_fmt *world);
extern void calc_hpd_credibility(bayes_fmt *bayes,long locus, long numpop2, long numparam);
#endif
