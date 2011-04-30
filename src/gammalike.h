#ifndef _GAMMALIKE_H_
#define _GAMMALIKE_H_
// gamma deviated mutation rate among loci
//
//
// code for the calculation of the gamma likelihood
// and for its derivatives
#include "migration.h"

MYREAL gamma_loci_like (helper_fmt * helper, MYREAL *oparam,
                                   MYREAL *olparam, MYREAL denom);
MYREAL gamma_locus_like (nr_fmt * nr, MYREAL *oparam, MYREAL *olparam,
                                    MYREAL denom, long locus);
void gamma_loci_derivative (helper_fmt * helper);
void gamma_locus_derivative (helper_fmt * helper, long locus);
void gamma_loci_difference (helper_fmt * helper);
#endif
