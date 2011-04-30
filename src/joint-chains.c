/* \file joint-chains.c 

Calculates Geyer's reverse logistic regression

like(param) = Sum_g^G( p(g|param)/ denom)

denom = Sum_j^chains (n_j P(g|param0_j) / L(param_j)
                      
                      \- pick param
                      \- calc chainparamlike
                      \- solve paramlike iteratively
                      \- maximize paramlike
                      
                      
*/

/* Joint chain estimator
Peter Beerli January 2000
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
$Id: joint-chains.c 886 2007-08-18 03:09:50Z beerli $
*/

#include "migration.h"
#include "broyden.h"
#include "combroyden.h"
#include "sighandler.h"

#define EPSILON5 0.00001

///
/// Precalculate the prob(G|param0) for a replicate and a locus
void
create_multiapg0 (MYREAL *apg0, nr_fmt * nr, long rep, long locus)
{
    long g, r;
    MYREAL *tmp;
    MYREAL tmpmax; //overflow protector
    MYREAL reps = LOG ((MYREAL) (nr->repstop - nr->repstart));

    tmp = (MYREAL *) mycalloc (nr->repstop, sizeof (MYREAL));
    // for each tree summary calculate the Sum over replicates
    for (g = 0; g < nr->world->atl[rep][locus].T; g++)
    {
        tmpmax = -MYREAL_MAX;
        for (r = nr->repstart; r < nr->repstop; r++)
        {
            tmp[r] = -reps +
                     probG (nr->world->atl[r][locus].param0,
                            nr->world->atl[r][locus].lparam0,
                            &nr->world->atl[rep][locus].tl[g], nr, locus)
                     - (nr->world->chainlikes[locus][r]);

            if (tmp[r] > tmpmax)
                tmpmax = tmp[r];
        }
        apg0[g] = 0.0;
        for (r = nr->repstart; r < nr->repstop; r++)
            apg0[g] += EXP (tmp[r] - tmpmax);
        apg0[g] = LOG (apg0[g]) + tmpmax;
    }
    myfree(tmp);
}

///
/// Solve the Geyer (1994) equation system
/// \f[
///       L(\Theta_0,M_0) = \sum_{chains} \frac{P(G|\Theta_0,M_0)}{\sum_{chains}\frac{P(G|\Theta_0,M_0)}{L(\Theta_0,M_0)}}
/// \f]
void
interpolate_like (nr_fmt * nr, long locus)
{
    boolean alldone;
    boolean diffdone=FALSE;
    long j, z, r;
    long repdiff = nr->repstop - nr->repstart;
    MYREAL *newlike, *lparam;
    MYREAL *diff, *oldiff;
    MYREAL delta = 0.;
    MYREAL *oldlike;
    MYREAL *chainlike = nr->world->chainlikes[locus];
    //#ifdef INTEGRATEDLIKE
    //return;
    //#endif
    oldlike = (MYREAL *) mycalloc (repdiff, sizeof (MYREAL));
    memcpy (oldlike, chainlike, sizeof (MYREAL) * repdiff);

    newlike = (MYREAL *) mycalloc (repdiff, sizeof (MYREAL));
    diff = (MYREAL *) mycalloc (repdiff, sizeof (MYREAL));
    oldiff = (MYREAL *) mycalloc (repdiff, sizeof (MYREAL));
    lparam = (MYREAL *) mycalloc (nr->numpop2, sizeof (MYREAL));
    /* following material is reverse logistic regression
       a la Geyer 1994 */
    z = 0;
    alldone = FALSE;
    // minimize the difference between oldlike and newlike, using oldiff
    memset (oldiff, 0, sizeof (MYREAL) * repdiff);
    while (!alldone && z++ < 10000)
    {
        alldone = TRUE;

        for (r = nr->repstart; r < nr->repstop; r++)
            create_multiapg0 (nr->apg0[r][locus], nr, r, locus);
        
        for (j = nr->repstart; j < nr->repstop; j++)
        {
            newlike[j] = calc_locus_like (nr, nr->atl[j][locus].param0,
                                          nr->atl[j][locus].lparam0, locus);
            diff[j] = fabs (newlike[j] - oldlike[j]);
            if (delta < diff[j])
                delta = diff[j];
            if (delta > EPSILON)
                alldone = FALSE;
        }

        if (nr->world->options->progress)
        {
            if (z % 100 == 0 && nr->world->options->verbose)
                printf ("           Iteration%6li biggest difference = %f\n", z, delta);
        }
        
        memcpy (oldlike, newlike, sizeof (MYREAL) * repdiff);
        diffdone = TRUE;
        for (j = nr->repstart; j < nr->repstop; j++)
        {
            if (fabs (oldiff[j] - diff[j]) > EPSILON)
            {
                diffdone = FALSE;
                break;
            }
        }
        if (diffdone)
            break;
        else
            memcpy (oldiff, diff, sizeof (MYREAL) * repdiff);
        // reset delta to the first log(L) difference, this will enter again in the loop 
        // above
        delta = diff[0];
    }
    if ((nr->world->options->progress || diffdone) && delta > EPSILON)
    {
        printf ("           Reweighting operation converged to\n");
        printf ("           constant multiplier %f in %li cycles\n", delta, z);
    }
    myfree(newlike);
    myfree(oldlike);
    myfree(diff);
    myfree(oldiff);
    myfree(lparam);
}
