/*
 *  pretty.h
 *  part of migrate
 *  will be PDF printout system
 *  Created by Peter Beerli on 7/25/05.
 *  Copyright 2005 Peter Beerli. All rights reserved.
 *
 */
#ifdef PRETTY
#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#include "profile.h"

#undef USE_ENCRYPTION
#define HAVE_BOOLEAN 
//#include <libharuc.h>
#include "haru/libharuc.h"
#define SEARCHUP 0
#define SEARCHDOWN 1
// GENERAL
// initialize PDF structures
extern void pdf_master_init(world_fmt *world, option_fmt *options, data_fmt *data);
extern int pdf_init(void);
// print title on every new page
extern int pdf_new_page(char *title);
// print master title on first page 
extern int pdf_master_title(char *title, float *orig_page_height, float *orig_left_margin);

extern void pdf_print_end_time(float *page_height);
// OPTIONS
// print options onto the first page
extern void pdf_print_options(world_fmt * world, option_fmt *options, data_fmt * data, 
			      float *orig_page_height, float *orig_left_margin);

// DATA
// print data summary
extern void pdf_print_data_summary(world_fmt * world, option_fmt *options, data_fmt * data, 
                              float *orig_page_height, float *orig_left_margin);
// print data
extern void pdf_print_data (world_fmt * world, option_fmt * options, data_fmt * data, float *orig_page_height, float *orig_left_margin);
extern void pdf_print_spectra(data_fmt *data, MYREAL ***freq, long * maxalleles);
extern void pdf_print_averageheat(world_fmt **universe, option_fmt *options);

// BAYES OUTPUTS
// print bayes table
extern void pdf_print_bayestable(world_fmt *world);
// print histogram at location lx ly with widht and height
extern void pdf_histogram(MYREAL *binvals, char *set50, char *set95, long bins, float bindelta, float binmin, float binmax, float lx, float ly, float width, float height, boolean nofreq);
extern void pdf_histogram_plus(MYREAL *binvals, MYREAL *std, char *set50, char *set95, long bins, float bindelta, float binmin, float binmax, float lx, float ly, float width, float height, MYREAL scaler, boolean nofreq, world_fmt * world, float *confidence, long topop);
// print histogram for each locus and overall loci
extern float pdf_loci_histogram(world_fmt *world);
// print acceptance ratios
extern void pdf_bayes_print_accept(world_fmt *world);
// print ESS and autocorrelation
extern void pdf_bayes_print_ess(world_fmt *world);
extern void pdf_bayes_factor_header(world_fmt *world, option_fmt *options);
extern void pdf_bayes_factor_rawscores_header(world_fmt *world, option_fmt *options);
extern void pdf_bayes_factor_rawscores(long locus, MYREAL rawtermo, MYREAL beziertermo, MYREAL harmo);
void pdf_bayes_factor_rawscores_harmo(long locus, MYREAL harmo);
extern void pdf_bayes_factor(world_fmt *world, MYREAL tsum, MYREAL tsum2, MYREAL hsum, MYREAL asum, long maxreplicate, MYREAL scaling_factor);
extern void pdf_burnin_stops(world_fmt *world, long maxreplicate);
// MIGHIST output
//extern void 
//pdf_mig_histogram(histogram_fmt ** histogram,
//                  plotfield_fmt ** plotfield, long loci, long nmigs,
//                  long bins, long *sum, MYREAL ***migtable,  boolean precalc, world_fmt *world);
//extern void pdf_print_mighist_table (world_fmt * world, MYREAL ***migtable,
extern void pdf_print_event_table (world_fmt * world, float *migtable);
extern void pdf_print_eventtime_table(world_fmt *world);
extern void pdf_print_time_table(world_fmt *world, 
				 float *meantime, float *mediantime, float *stdtime, float *freq,
				 boolean mrca);
extern void pdf_event_histogram(long loci, long numparams,  world_fmt *world);
// SKYLINE plot
extern void pdf_skyline_histogram(long loci, long numparams,  world_fmt *world, boolean enlarged);

// ML OUPUTS
// print mcmc table
//extern void pdf_print_mcmctable(float *page_height, float page_width, world_fmt *world, data_fmt *data, option_fmt *options);
extern void pdf_print_results (float *page_height, world_fmt ** universe, option_fmt * options, data_fmt * data);
extern void pdf_print_profile_title(world_fmt *world);
extern void pdf_print_profile_table(float left_margin, float * page_height, char method, char *buffer, world_fmt *world);
extern long prepare_profile_table(char method, long which, MYREAL *likes, MYREAL **param, world_fmt *world, long bufsize);
extern void pdf_print_profile_percentile (world_fmt * world);

extern void pdf_histogram_legend(void);
extern void pdf_print_lrt_header(MYREAL aicfull, long aicfullparamnum,  float *page_width,  float *page_height);
extern void pdf_print_LRT_box(world_fmt* world,char **box,long size,
			      MYREAL like0, MYREAL like1,MYREAL lrt, MYREAL chiprob, MYREAL chiprob2, MYREAL aic, long df, 
			      long aicparamnum, float *page_height, float *page_width);

void pdf_print_stored_warnings(world_fmt *world);
// write PDF content to file
extern int pdf_write_file(option_fmt *options);
#ifndef HAVE_STRSEP
char *strsep(char **, const char *);
#endif
#endif
