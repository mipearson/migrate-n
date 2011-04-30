#ifndef MIGRATE_DEFINITION
/*definitions.h*/
#define MIGRATE_DEFINITION
/*-----------------------------------------------------------------
  Maximum likelihood estimation of migration rates 
  using coalescent trees
 
  Peter Beerli
  School of Computational Science and Biological Sciences
  Florida State University
  Tallahassee, FL 32306-4120
  beerli@fsu.edu
 
  With help of Joe Felsenstein (joe@gs.washington.edu),
  Mary Kuhner and Jon Yamato (mkkuhner@gs.washington.edu)
 
  
Copyright 2002 Peter Beerli and Joseph Felsenstein
Copyright 2003-2008 Peter Beerli
Copyright 2009-2010 Peter Beerli and Michal Palczewski
 
$Id: definitions.h 1862 2011-04-07 13:18:47Z beerli $
 
  *----------------------------------------------------------------
  */
#ifndef MIGRATEVERSION
#define MIGRATEVERSION "3.2.8"
#else
#ifdef SGI
#undef MIGRATEVERSION
#define MIGRATEVERSION "3.2.8"
#endif
#endif
#ifndef MIGRATESUBVERSION
#define MIGRATESUBVERSION "x"
#endif

#define MAINTAINER "Peter Beerli <beerli@fsu.edu>"
/*------------------------------------------------
  Adjusting the source to the system configuration
*/
#ifdef use_floats
// use this when using floats 
#define USE_MYREAL_FLOAT 0
#define MYREAL float
#define MYFLOAT float
#define MYREAL_MAX 3.40282347e+38F
#else
#define USE_MYREAL_DOUBLE 1
#define MYREAL double
#define MYFLOAT float
#define MYREAL_MAX 1.7976931348623157e+308
#endif
#ifdef HAVE_CONFIG_H
#include "conf.h"
#endif
/* we don't want to have the sun specific define of "sun"*/
#ifdef sun
#undef sun
#endif

#define FGETS  myfgets
#define ZNZFGETS  myznzfgets
#define FGETS2  myfgetssafe
#ifdef MPI
#include <mpi.h>
#define FPRINTF mpi_fprintf
#define LARGEFPRINTF mpi_fprintf2
#else
#define FPRINTF fprintf
#define LARGEFPRINTF fprintf2
#endif

/* random number system */
#define LCG
#ifdef MERSENNETWISTER
#undef LCG
#endif
#ifdef SPRNG
#undef LCG
#endif

#ifdef __MWERKS__
    #include "conf.h"
    #undef HAVE_LGAMMA
    #define HAVE_STRFTIME 1
    #define HIGHBITS
    #define LAMARC_MALLOC
    #ifdef macintosh
    #define MAC
    #include <ansi_prefix.mac.h> /*fixes time problems */
#else /*macintosh */
    #define DOS
    #undef HAVE_STRINGS_H
    #undef HAVE_LGAMMA
    #define HAVE_STRFTIME 1
    #define HIGHBITS
    #define LAMARC_MALLOC
    #endif /*windows */
#endif /*MWERKS*/

#ifdef WIN32
    #define DOS
    #undef HAVE_STRINGS_H
    #undef HAVE_LGAMMA
    #define HAVE_STRFTIME
    #undef HAVE_STRSEP
    #include "windows_timveval.h"
#endif

/* we have found no strftime() and replace the time
   output with blanks*/
#ifndef HAVE_STRFTIME
    #define NOTIME_FUNC
#endif

#define STDOUTNUM 0
#define INFILE "infile"
#define INFILENUM 1
#define WEIGHTFILE "weightfile"
#define WEIGHTFILENUM 2
#define CATFILE "catfile"
#define CATFILENUM 3
#define SEEDFILE "seedfile"
#define SEEDFILENUM 4
#define PARMFILE "parmfile"
#define PARMFILENUM 5
#define OUTFILE "outfile"
#define OUTFILENUM 6
#define LOGFILE "logfile"
#define LOGFILENUM 7
#define TREEFILE "treefile"
#define TREEFILENUM 8
#define UTREEFILE "usertree"
#define UTREEFILENUM 9
#define INTREE "intree"
#define INTREENUM 10
#define OUTTREE "outtree"
#define OUTTREENUM 11
#define MATHFILE "mathfile"
#define MATHFILENUM 12
#define SUMFILE "sumfile"
#define SUMFILENUM 13
#define MIGHISTFILE "mighistfile"
#define MIGHISTFILENUM 14
#define DISTFILE "distfile"
#define DISTFILENUM 15
#define GEOFILE "geofile"
#define GEOFILENUM 16
#define BOOTFILE "bootstrapfile"
#define BOOTFILENUM 17
#define UEPFILE "uepfile"
#define UEPFILENUM 18
#define AICFILE "aicfile"
#define AICFILENUM 19
#define BAYESFILE "bayesfile"
#define BAYESFILENUM 20
#define BAYESMDIMFILE "bayesallfile.gz"
#define BAYESMDIMFILENUM 22
#define PDFOUTFILE "outfile.pdf"
#define PDFOUTFILENUM 21
#define SKYLINEFILE "skylinefile"
#define SKYLINEFILENUM 23
#define MYMAXFILENUM 24
#define DIVFILE "divfile"
#define DIVFILENUM 24
#define TIPDATEFILE "datefile"
#define TIPDATEFILENUM 25
#define MIXFILE "mixfile"
#define MIXFILENUM 26
#define BAYESNUMBIN 1500
#define PRETTY_MAX 0
#define PRETTY_P99 1
#define PRETTY_P99MAX 3
#define PRETTY_P100 4
#define THETAPRIOR 0
#define MIGPRIOR   1
#define RATEPRIOR  2
#define PRIOR_SIZE 3 /*used because of slice sampler*/


/* includes: */
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#include <string.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#ifdef MAC
    #include "mac_interface.h"
#endif /*MAC*/

#include <time.h>
#include <math.h>
#include <limits.h>

#ifdef WIN32
    #define MYINLINE __inline
    #undef DBL_EPSILON
//    #undef MYREAL_MAX
#define __restrict /*__restrict*/
    #include <float.h>
    #define MYISNAN _isnan
    #define MYFINITE _finite
    #define MYSNPRINTF _snprintf
#else  
    #define MYISNAN  isnan
    #define MYFINITE finite
    #define MYSNPRINTF snprintf
    #ifdef __MWERKS__
		#define MYINLINE /*inline*/
		#undef DBL_EPSILON
		#undef MYREAL_MAX
		#include <float.h>
	#else
		#define MYINLINE inline
	#endif
#endif

/*if we use this then we guard against allocation problems */
/* redefinitions of calloc, malloc, realloc are in sighandler.h*/
#ifdef HAVE_MALLOCWRAP
#define LAMARC_MALLOC
#endif

#ifdef HAVE_DEBUGMPI
/* add print messages to mpi commands */
#define DEBUG_MPI
#endif

#ifdef HAVE_LGAMMA
#define LGAMMA lgamma
#else
#define LGAMMA mylgamma
#endif

//#include "sighandler.h"
#define FClose(file) if (file) fclose(file) ; file=NULL
#define znzClose(file) if (file) znzclose(file) ; file=NULL
#ifdef WIN32
#ifndef WINDOWS
#define WINDOWS
#endif
#include <windows.h>
#else
typedef int boolean;
//typedef unsigned char boolean;
#endif
#ifndef TRUE
#define TRUE    (boolean) 1
#endif

#ifndef FALSE
#define FALSE   (boolean) 0
#endif

#define INVALID -999

#define SETBITS 32
/* for snp data with panel see sequence.c*/
#define FIRST 0
/* the first population is the panel for snps*/
#define PANEL 0
#define SEQUENCETYPES "snufh" /*s=sequence,n=linked snp, u=unlinked snp, f=ancestral reconstruction h=hapmap snps*/
#define SNPTYPES      "nuh"
#define ALLELETYPES   "amb" /*a=infinite allele,m=stepwise, b=brownian */

#define CRLF "\r\n"

#ifdef COMPAQ
#define EXP(a)  (((a)< -100) ? 0.0 : exp ((a)))
#else
#ifdef FAST_EXP
#define EXP(a)  fast_exp ((a))
#else
#define EXP(a)  exp ((a))
#endif
#endif

// I was experimenting with a fast log implementation
// but it is
//not fast enough #define LOG(a) ((MYREAL) fast_log(((float) a)))
// so back to standard log
#ifdef FAST_EXP
#define LOG(a) ((MYREAL) fast_log(((float) a)))
#else
#define LOG(a) log((a))
#endif

/* defines for speedier calc_like when we do NOT use undefined SLOWNET*/
#ifdef SLOWNET
#define CALCLIKE (*calc_like)
#define CALCGRADIENT (*calc_gradient)
#define SETUPPARAM0 (*setup_param0)
#else
#define CALCLIKE calc_loci_like
#define CALCGRADIENT combine_gradient
#ifdef MPI
#define SETUPPARAM0 setup_parameter0_mpi
#else
#define SETUPPARAM0 setup_parameter0_standard
#endif
#endif
/*=================================================================
  SPECIFIC FOR MIGRATE
  
  */
#define MASTER 0
#define FIRSTWORKER 1
#define STOP_REPLICATORS -1234
/* number of k int k-allele model, the number
   MUST really be bigger than the number of 
   observed alleles, much bigger does not harm */
#define MAXALLELES            93L
#define SMALLBUFSIZE    255L
#define ONEMEGABYTE   1024000
#define MAXBUFSIZE  10000000L
//#define LINESIZE        1024L /* setting this smaller can break */
#define LINESIZE      100000L /* setting this smaller can break */
#define LONGLINESIZE  100000L /* setting this smaller can break */
#define SUPERLINESIZE 100000L /*used to read many many loci-sites and indiv numbers*/ 
#define STRSIZE              255L /* setting this smaller can break */
#define DEFAULT_NMLENGTH      10L /* length of individual names */
#define DEFAULT_ALLELENMLENGTH 6L /* length of allele names */
#define DEFAULT_POPNMLENGTH  100L /* length of world names */
#define NUMPOP 2L
#define BURNINPERIOD 10000L
#define SCALEINTERVAL 2L
#define DEFAULTALLOCSIZE 2L
/* don't change below here ------------------------------------
   some other constants */
/* numerical borders and epsilons*/
#ifndef HUGE
#define HUGE 1e200
#endif 
#define SICK_VALUE              -1
#define ONE                     1L
#define TWO                     2L
#define THREE                   3L
#define FOUR                    4L
#define SMALL_VALUE       10e-21
#define EPSILON         0.000001 /* a small number */
#define EPSILON4        0.000099 /* another small number */
#define SMALLEPSILON       1e-15 /* a smaller number */
#define BIGEPSILON         0.001 /* a not so small number */
#define PERCENTILETOLERANCE 0.01 /* difference larger than this are flagged FAILED*/
#define LONGCHAINEPSILON  10e100 /* an unreasonable big number,
so that the given number of
long chains is used */
#define CHAINVARIANCEDELTA 0.01 /*difference of variances of parts of a chain [see burnin_autostop]*/
#define ESSMINIMUM 200  /*effective sample size minimum see burnin_autostop*/
#define GELMAN_MYSTIC_VALUE 1.2
#define PLUSCHAIN 10
#ifndef MAXLONG
#define MAXLONG ((long)0x7fffffff)
#endif
#ifndef DBL_EPSILON
#include <float.h>
#ifndef DBL_EPSILON
#define DBL_MAX ((MYREAL)1.7976931348623157e308)
#define DBL_EPSILON 2.2204460492503131e-16
#endif
#endif
#ifndef FLT_EPSILON
#define FLT_EPSILON 1.19209290e-07F
#endif
/* some math constants */
#define HALF 0.5
#define QUARTER 0.25
#define LOGDBL_EPSILON -36.04365338911715608 /*N[Log[DBL_EPSILON] ,30] */
#define LOG1 0.
#define LOG2 0.693147180559945309417232121458 /*N[Log[2] ,30] */
#define LOG2PIHALF -0.918938533204672741780329736406 /*N[Log[1/Sqrt[2 Pi]] ,30] */
#define TWOPI 6.28318530717958647692528676656 /*N[Log[2 Pi] ,30] */
#define PI 3.14159265358979323846264338328 /*N[Pi,30]*/
#define INV2PI 0.159154943091895335768883763373 /*N[1/(2 pi)]*/
#define ROOTLENGTH         10000
#define FIRSTCHAIN            -1
#define FIRSTSTEP             -1
/* min/max functions */
#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#define MYMIN3(a,b,c) (((a)<(b))&&((a)<(c)) ? (a) : (((b)<(c)) ? (b) : (c)))
#define MYMAX3(a,b,c) (((a)>(b))&&((a)>(c)) ? (a) : (((b)>(c)) ? (b) : (c)))
/*string conversion functions
  these replace all the atoi atol etc and use the strto[lif] function*/
#define ATOL(a) strtol((a),(char **) NULL, 10)
#define ATOI(a) strtol((a),(char **) NULL, 10)
#define ATOF(a) strtod((a),(char **) NULL, 10)

/* prior setting */
#define SLICE        5
#define GAMMAPRIOR   4
#define MULTPRIOR    3
#define WEXPPRIOR    2
#define EXPPRIOR     1
#define UNIFORMPRIOR 0

/* theta (4 Ne mu) related material */
// ONLY good boundaries when used with DNA
#define SMALLEST_THETA    1e-10
#define DNA_GUESS_THETA    0.01
#define BIGGEST_THETA      1e4
#define SMALLEST_MIGRATION 0.0  
#define DNA_GUESS_MIG      100
#define BIGGEST_MIGRATION   10e9
#define SMALLEST_RATE 1e-20
#define BIGGEST_RATE  1e20 
// defines the value that is used as a minimum summary statistic for migration events on all trees
// version 1.2.4 had a value of 0.00000001, it seems that values of 0.1 are too large.
// in version -1.6.9 this is set to zero, in 1.7 trial with 0.01 PB Dec 11 2002  [this was correct for Nm]
// was set to 1 (M) for some time , but this might be biasing upwards, changing back to small value 
#define MINMIGSUMSTAT 0.000001
/* the migration limit is per population, should be plenty
   for moderate sample sizes */
#define MIGRATION_LIMIT     1000
#define SMALLEST_PROB     1e-100
#define SMALLEST_GAMMA      1e-3
#define BIGGEST_GAMMA        1e9
#define FRACTION_ALONG      0.66
#define TIMELIST_GUESS       500
#define SAMPLETREE_GUESS       1
#define START_ALPHA          10.0
/* random material */
#define AUTO                   0 /*+ use time() for seed */
#define NOAUTO                 1 /*+ seed in parmfile   */
#define NOAUTOSELF             2 /*+ seed in seedfile   */
/* definitions for the first theta0 values */
#define OWN                    0 /*+ start values in parmfile */
#define WATTERSON              1 /*- start values are Watterson estimate*/
#define EWENS                  2 /*- start values are Ewens estimate*/
#define FST                    3 /*+ start values are FST estimate */
#define NRANDOMESTIMATE         4 /*+ start values are RANDOM estimate */
#define URANDOMESTIMATE         5 /*+ start values are RANDOM estimate */
#define PARAMGRID              6 /*+ use a range of values + replicate=yes:# */
/* defines for the migration0 values */
/* define OWN                  0   + start values in parmfile */
/* define FST                  3   + start values are FST estimates*/
/* define RANDOMESTIMATE       4   + start values are RANDOM estimate*/
#define SLATKIN                1 /*- start values with slatkin's method*/
/* definitions for the type of FST we use*/
#define THETAVARIABLE         'T'
#define MVARIABLE             'M'
/* definitions for the sankoff procedure */
/* the SANKOFF_DELTA should be smaller than any cost_ij*/
#define SANKOFF_DELTA 0.1
/* Maximizer and gamma things */
#define STD22 1.0 /*0.1*/  /*"standard deviation" for penalizer 2*std*std=2*0.223607^2 */
#define STD2 1.0 /*0.05*/  /* std*std = 0.223607^2 */
#define LOGSTD 0.0 /*-1.49786613677699549672*/ /*LOG(std)*/
/* #define INVTWOSQRTPILOGSTD currently UNUSED*/ /*0.57892760357232275494*/ /* log(1/(std sqrt[2 pi])) */
#define INVTWOSQRTPILOGSTD -0.91893853320467274178
#define NTRIALS             1000
#define LOCI_NORM          0.00001
#define GAMMA_INTERVALS       10
/* definitions of migration model */
#define MATRIX                 1 /*+ */
#define MATRIX_SYMMETRIC      11 /*+ */
#define MATRIX_SAMETHETA      12 /*+ */
#define MATRIX_ARBITRARY      13 /*+ */
#define ISLAND                 2 /*+ */
#define ISLAND_VARTHETA       21 /*+ */
#define STEPSTONE              3 /*- */
#define CONTINUUM              4 /*- */
#define NEIGHBOR               5 /*- */
/* replication stuff */
#define SINGLECHAIN            0 /*+ */
#define MULTIPLECHAIN          1 /*+ */
#define MULTIPLERUN            2 /*+ */
/* heating stuff */
#define HEATED_CHAIN_NUM       4
#define HEATED_THREAD_NUM      4
#define COLD                   1 /*+ */
#define WARM                   4 /*+ */
#define HOT                    7 /*+ */
#define VERYHOT                10 /*+ */
#define NOTADAPTIVE            0 /*static heating scheme*/
#define STANDARD               1 /*adaptive heating scheme*/
#define BOUNDED                2 /*bounded adaptive scheme*/
/* maximization material*/
#define SINGLELOCUS            0 /*+ */
#define MULTILOCUS             1 /*+ */
#define PROFILE                2 /*+ */
/* print and plot material*/
#define MULTILOCUSPLOT         (boolean) 1
#define SINGLELOCUSPLOT		   (boolean) 0
#define MAXPRINTVALLENGTH 16 /*+ maximal number of printpos */
#define PLOTALL                0 /*+ plot to outfile and mathfile */
#define PLOTOUTFILE            1 /*+ plot only to outifle */
#define PLOT4NM                0 /*+ use 4Nm instead of M=m/mu */
#define PLOTM                  1 /*+ use M=m/mu */
#define PLOTSCALELOG           0 /*+ plot scale in log 10 units */
#define PLOTSCALESTD           1 /*+ plot scale standard units */
#define PLANEINTERVALS        36 /*+ intervals for plot */
#define PLANESTART        0.0001 /*+ plot axes start */
#define PLANEEND           100.0 /*+ plot axes end */
#define PAGEFEED               fprintf(outfile,"\n\f\n")
#define PAGEFEEDWORLD               fprintf(world->outfile,"\n\f\n")
#define PAGEFEEDW                   fprintf(universe[0]->outfile,"\n\f\n")
#define START                  0
#define STOP                   1
#define CONT                   3
/* tree-print options */
#define _NONE                   0 /*+ */
#define ALL                    1 /*+ */
#define BEST                   2 /*+ */
#define LASTCHAIN              3 /*+ */
/* likelihood ration stuff*/
#define LRATIO_STRINGS      1000 /*+ */
#define HUNDRED             100
/* profile likelihood stuff*/
/* ALL and NONE are alread defined*/
#define TABLES                 2 /*+ */
#define SUMMARY                3 /*+ */
#define MAX_PROFILE_TRIALS   100 /*+ */
/* likelihood ratio stuff*/
#define MLE 0
#define ARBITRARY 1
#define MAXPOP   200  /*n*(n-1)+n) */
#define MAXPARAM 20000
#define LNMAXPARAM 9.9034875525361280455 
#define HEADER 1
#define NOHEADER 0
/* migration histogram material*/
#define MIGHIST_ELEM 5
#define MIGHIST_YSIZE 30
/* stop while loops after some time*/
#define PANIC_MAX           1000
/* microsatellite stuff*/
#define MICRO_THRESHOLD      200
#define MAX_MICROSTEPNUM     200
#define XBROWN_SIZE            3
#define SINGLESTEP             1
#define MULTISTEP              2
/* sequence stuff*/
#define MAXCATEGS              9
#define MANYCATEGS             2
#define ONECATEG               1

#define NUC_A                  0
#define NUC_C                  1
#define NUC_G                  2
#define NUC_T                  3
#define NUC_GAP                4
#define NUC_R                  5
#define NUC_Y                  6
#define NUC_AR                 7
#define NUC_CY                 8
#define NUC_GR                 9
#define NUC_TY                10
#define BASEFREQLENGTH        11
#ifndef LAMARC_MAC_INTERFACE
#define eventloop() /* eventloop would go here */
#endif

#define SUBLOCICHUNKS 10

#ifdef powerpc
#ifndef IBM
#include <Accelerate/Accelerate.h>
#endif
#endif
#undef MYINLINE 
#define MYINLINE /*inline*/
#undef PRETTY
#define PRETTY
#endif /*definitions */
