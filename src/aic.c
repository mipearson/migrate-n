/*! \file aic.c */
/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
	AIC model test   R O U T I N E S

	Peter Beerli 2001, Seattle
	beerli@fsu.edu

	Copyright 2001-2002 Peter Beerli and Joseph Felsenstein

	This software is distributed free of charge for non-commercial use
	and is copyrighted. Of course, we do not guarantee that the software
	works and are not responsible for any damage you may cause or have.

$Id: aic.c 1702 2010-06-14 14:18:13Z beerli $ 
-------------------------------------------------------*/

#include "migration.h"
#include "tools.h"
#include "broyden.h"
#include "combroyden.h"
#include "options.h"
#include "migrate_mpi.h"
#include "aic.h"
#include "sighandler.h"

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif


static void print_progressheader_aic(boolean progress, FILE *file, MYREAL mleaic, MYREAL numparam);

static void aic_score (aic_fmt ** aicvec, long *aicnum, nr_fmt * nr,
                long zero, long which, char *temppattern, MYREAL *param0,
                int migtype);

static boolean legal_pattern (char *matrix, long numpop);

void fast_aic_score (aic_fmt ** aicvec, long *aicnum, nr_fmt * nr,
                     long zero, long which, char *temppattern,
                     MYREAL *param0, int migtype);

static void add_aiclist (int migtype, long numparam, long freeparam, long remainnum, long zero, long m,
                  MYREAL *param0, char *temppattern, char *custm2,
                  long *aicnum, aic_fmt ** aicvec, nr_fmt * nr);

static void add_aicvec(MYREAL aic, aic_fmt **aicvec, long *aicnum, nr_fmt *nr, long numparam, long remainnum);

static void aic_print_driver(aic_fmt *aicvec, MYREAL mleaic, long aicnum, nr_fmt *nr, world_fmt *world);

static void aic_score_driver(aic_struct *aic, MYREAL mle, MYREAL mleaic, nr_fmt *nr, world_fmt *world, char *temppattern);

//void akaike_information (world_fmt * world, long *Gmax);

static void print_header_aic (nr_fmt * nr, MYREAL mleaic);

static void print_aicfile(aic_fmt **aicvec, long *aicnum, nr_fmt *nr);

static void print_progress_aic(boolean add
						, aic_fmt **aicvec, long *aicnum, nr_fmt *nr, long numparam, long freeparam);

//static long find_paramnum(world_fmt *world, char *connect);

///
/// changes a linear matrix of ddd mm mm mm
/// to a diagonal matrix dmm mdm mmd
/// it destroys pattern
static char * reshuffle (char *pattern, char *origpattern, long numpop)
{
    long space = 0;
    long i, j, z = numpop;
	
    for (i = 0; i < numpop; i++)
    {
        for (j = 0; j < numpop; j++)
        {
            if (i == j)
                pattern[i * numpop + j + space] = 'x';
            else
                pattern[i * numpop + j + space] = origpattern[z++];
        }
        pattern[i * numpop + j + space] = ' ';
        space++;
    }
    return pattern;
}

/// print the progress header for AIC based model selection
void print_progressheader_aic(boolean progress, FILE *file, MYREAL mleaic, MYREAL numparam)
{
    if (progress)
    {
        FPRINTF (file, "\n\n");
        FPRINTF (file, "           Selecting the best migration model for this run,\n");
        FPRINTF (file, "           This may take a while!\n");
        FPRINTF (file, "           Checking all parameter combinations\n");
        FPRINTF (file, "           with  (AIC= -2 Log(L(Param)+n_param)<%f+%f\n",
                 mleaic, numparam);
    }
}

///
/// print the header for the AIC based model selection
/// the routine reports also progress to screen and logfile
void
print_header_aic (nr_fmt * nr, MYREAL mleaic)
{
    world_fmt *world = nr->world;
	
    print_progressheader_aic(world->options->progress, stdout,mleaic, world->options->aicmod*world->numpop2);
    print_progressheader_aic(world->options->writelog, stdout,mleaic, world->options->aicmod*world->numpop2);
	
    FPRINTF (world->outfile, "\n\n\n Akaike's Information Criterion  (AIC)\n");
    FPRINTF (world->outfile, "=========================================\n\n");
    FPRINTF (world->outfile, "[Linearized migration matrix, x=diagonal]\n");
}

///
/// Calculates AIC score
/// \callgraph
void aic_score_driver(aic_struct *aic, MYREAL mle, MYREAL mleaic, nr_fmt *nr, world_fmt *world, char *temppattern)
{
    char * savecustm;
    char *savecustm2;
	
    savecustm = (char *) mycalloc (2 * world->numpop2, sizeof (char));
    savecustm2 =  savecustm + world->numpop2;
    memcpy (savecustm, world->options->custm2, sizeof (char) * nr->numpop2);
    memcpy (savecustm2, world->options->custm, sizeof (char) * nr->numpop2);
	
    aic->aicnum = 1;
    aic->aicvec = (aic_fmt *) mycalloc (aic->aicnum, sizeof (aic_fmt));
    aic->aicvec[0].mle = mle;
    aic->aicvec[0].aic = mleaic;
    aic->aicvec[0].lrt = 0.0;
    aic->aicvec[0].prob = 1.0;
    aic->aicvec[0].probcorr = 1.0;
    aic->aicvec[0].numparam = nr->partsize;
    aic->aicvec[0].pattern = (char *) mycalloc (nr->partsize + 1, sizeof (char));
    memcpy (aic->aicvec[0].pattern, world->options->custm2, sizeof (char) * nr->partsize);
    if (world->options->fast_aic)
    {
        fast_aic_score (&aic->aicvec, &aic->aicnum, nr, 0, world->numpop,
                        temppattern, aic->param0, world->options->aictype[0]);
        memcpy (world->options->custm2, savecustm, sizeof (char) * nr->numpop2);
        if(world->options->aictype[1]!='\0')
            fast_aic_score (&aic->aicvec, &aic->aicnum, nr, 0, world->numpop,
                            temppattern, aic->param0, world->options->aictype[1]);
    }
    else
    {
        //find aic scores in a branch-and-bound fashion
        // with some parameters set to zero, this needs more
        // investigation because of boundary problems
        aic_score (&aic->aicvec, &aic->aicnum, nr, 0, world->numpop,
                   temppattern, aic->param0, world->options->aictype[0]);
        memcpy (world->options->custm2, savecustm, sizeof (char) * nr->numpop2);
        // aic scores based on averaging M (not 4Nm)
        if(world->options->aictype[1]!='\0')
        {
            aic_score (&aic->aicvec, &aic->aicnum, nr, 0, world->numpop, temppattern,
                       aic->param0, world->options->aictype[1]);
            memcpy (world->options->custm2, savecustm, sizeof (char) * nr->numpop2);
        }
    }
    memcpy (world->options->custm2, savecustm, sizeof (char) * nr->numpop2);
    memcpy (world->options->custm, savecustm2, sizeof (char) * nr->numpop2);
    myfree(savecustm); // savecustom2 is also freed, because it's part of savecustm
}

///
/// print aic results, assume they are already ordered by qsort()
/// printing is:
/// pattern {migration matrix}
/// aic value, number of parameters, MLE of pattern, Prob of standard LRT, Prob of weighted chisquare
void aic_print_driver(aic_fmt *aicvec, MYREAL mleaic, long aicnum, nr_fmt *nr, world_fmt *world)
{
  //#ifdef MYREAL == float
  //const MYREAL eps = FLT_EPSILON ;
  //#else
  const MYREAL eps = DBL_EPSILON ;
  //#endif
  long i;
  char *temppattern;
  boolean mldone=FALSE;
  char dashes[101] = "----------------------------------------------------------------------------------------------------";
  temppattern =  (char *) mycalloc (world->numpop2 + 1 + world->numpop, sizeof (char));
  FPRINTF (world->outfile,
             "%-*.*s  AIC     #param   Log(L)   LRT        Prob    Probc\n",
             (int) MAX (18, nr->partsize + nr->numpop),
             (int) (nr->partsize + nr->numpop), "Pattern");
    FPRINTF (world->outfile,
             "%-*.*s ---------- ---- ---------- ---------- ------- -------\n",
             (int) MAX (18, nr->partsize + nr->numpop),
             (int) (nr->partsize + nr->numpop), dashes);
    
    for (i = 0; i < aicnum; i++)
    {
        FPRINTF (world->outfile, "%-*.*s %10.5f %4li % 10.5f % 10.5f %6.5f %6.5f\n",
                 (int) MAX (18, nr->partsize + nr->numpop),
                 (int) (nr->partsize + nr->numpop),
                 reshuffle (temppattern, aicvec[i].pattern, nr->numpop),
                 aicvec[i].aic, aicvec[i].numparam, aicvec[i].mle,
                 aicvec[i].lrt, aicvec[i].prob, aicvec[i].probcorr);
        if ((aicvec[i].aic - mleaic > eps) && !mldone)
        {
            mldone = TRUE;
            FPRINTF (world->outfile, "%-*.*s %53.53s\n",
                     (int) MAX (18, nr->partsize + nr->numpop),
                     (int) (nr->partsize + nr->numpop),dashes,
					 dashes);
        }
        myfree(aicvec[i].pattern);
    }
    FPRINTF (world->outfile, "\n\n\n");
    myfree(temppattern);
}

///
/// this drives the aic calculation and is called from main.c
/// \callgraph
void
akaike_information (world_fmt * world, long *Gmax)
{
    aic_struct aic;
    nr_fmt *nr;
    long kind = world->loci > 1 ? MULTILOCUS : SINGLELOCUS;
    long repstop=0;
    long repstart=0;
	
    MYREAL mleaic;
    MYREAL mle;
    boolean multilocus=FALSE;
	
    char *temppattern;
	
    prepare_broyden (kind, world, &multilocus);
    world->options->migration_model = MATRIX_ARBITRARY;
	
    temppattern =  (char *) mycalloc (world->numpop2 + 1 + world->numpop, sizeof (char));
    aic.param0 = (MYREAL *) mycalloc (world->numpop2 + 1, sizeof (MYREAL));
	
    set_replicates (world, world->repkind, world->rep, &repstart, &repstop);
	
    if (kind == MULTILOCUS)
    {
        mle = world->atl[0][world->loci].param_like;
        mleaic = -2. * mle + 2. * find_paramnum(world,NULL);
    }
    else
    {
        mle = world->atl[repstop == 1 ? 0 : repstop][0].param_like;
        mleaic = -2. * mle + 2. * find_paramnum(world,NULL);
    }
    nr = (nr_fmt *) mycalloc (1, sizeof (nr_fmt));
	
    create_nr (nr, world, *Gmax, 0, world->loci, world->repkind, world->rep);
	
    SETUPPARAM0 (world, nr, world->repkind, repstart, repstop,
                 world->loci, kind, multilocus);
	
    print_header_aic (nr, mleaic);
	
    if (kind == MULTILOCUS)
        memcpy (aic.param0, nr->world->atl[0][nr->world->loci].param, sizeof (MYREAL) * nr->numpop2);
    else
        memcpy (aic.param0, nr->world->atl[repstop ==1 ? 0 : repstop][nr->world->locus].param, sizeof (MYREAL) * nr->numpop2);
	
    aic_score_driver(&aic, mle, mleaic, nr, world, temppattern);
	
    qsort ((void *) aic.aicvec, (unsigned long) aic.aicnum, sizeof (aic_fmt), aiccmp);
	
    aic_print_driver(aic.aicvec, mleaic, aic.aicnum, nr, world);
	
    myfree(aic.aicvec);
    (void) fflush (world->outfile);
    myfree(aic.param0);
    myfree(temppattern);
    destroy_nr (nr, world);
    //myfree(nr);
}

///
/// Does check all enumeration on one level and then picks only the best model to continue to
/// the next level
/// [this is different to the "branch-bound" algorithm]
void fast_aic_score (aic_fmt ** aicvec, long *aicnum, nr_fmt * nr,
                long zero, long which, char *temppattern,
                MYREAL *param0, int migtype)
{
  //#ifdef MYREAL == float
  //MYREAL numbermax = FLT_MAX ;
  //#else
  const MYREAL numbermax = MYREAL_MAX ;
  //#endif

    long m, ii;
    MYREAL likes = 0;
    MYREAL normd = 0;
    MYREAL aic;
    MYREAL borderaic = (*aicvec)[0].aic + nr->world->options->aicmod * nr->numpop2;
    char savecustm2;
    long remainnum = 0;
    //    boolean mylegal;
    char *custm2 = nr->world->options->custm2;
    char *scustm2;
    long numparam;
    long freeparam;
    aic_fmt *best;
	
    scustm2 = (char *) mycalloc (nr->partsize, sizeof (char));
    memcpy (scustm2, custm2, sizeof (char) * nr->partsize);
    if (migtype == 'm')
        remainnum = 1;
    numparam = zero;
    freeparam = (nr->numpop2 - numparam - 1 + remainnum);
    best = (aic_fmt *) mycalloc (nr->partsize, sizeof (aic_fmt));
    for (m = nr->numpop; m < nr->numpop2; m++)
    {
        best[m].aic = numbermax;
        best[m].numparam = m;
        if (scustm2[m] == migtype)
            continue;
        savecustm2 = custm2[m];
        custm2[m] = migtype;
        memcpy (nr->world->param0, param0, sizeof (MYREAL) * nr->numpop2);
        resynchronize_param (nr->world);
        if ((legal_pattern (nr->world->options->custm2, nr->numpop)))
        {
            (void) do_profiles (nr->world, nr, &likes, &normd, PROFILE,
                         nr->world->rep, nr->world->repkind);
            aic = -2. * nr->llike + 2. * freeparam;
            best[m].aic = aic;
            best[m].numparam = m;
            add_aiclist(migtype, numparam, freeparam, remainnum, zero, m,
                        param0, temppattern, custm2, aicnum, aicvec, nr);
        }
        else
        {   // illegal combination of parameters
            if (nr->world->options->progress)
                FPRINTF (stdout, "           F   %s %20s\n",
                         reshuffle (temppattern, custm2, nr->numpop), "-----");
            (void) fflush (stdout);
            if (nr->world->options->writelog)
                FPRINTF (nr->world->options->logfile, "           F   %s %20s\n",
                         reshuffle (temppattern, custm2, nr->numpop), "-----");
        }
        custm2[m] = savecustm2;
    }
    for (ii = nr->numpop; ii < nr->partsize; ii++)
    {
        if (best[ii].aic < borderaic && custm2[best[ii].numparam] != migtype)
        {
            custm2[best[ii].numparam] = migtype;
            fast_aic_score (aicvec, aicnum, nr, zero + 1, best[ii].numparam,
                            temppattern, param0, migtype);
        }
    }
    myfree(best);
    myfree(scustm2);
}

/// add model and AIC score to the list of models for the final printing
void add_aicvec(MYREAL aic, aic_fmt **aicvec, long *aicnum, nr_fmt *nr, long numparam, long remainnum)
{
    MYREAL lrt = -2. * (nr->llike - (*aicvec)[0].mle);
	
    *aicvec = (aic_fmt *)
		 myrealloc (*aicvec, sizeof (aic_fmt) * (*aicnum + 1));
    (*aicvec)[*aicnum].pattern = (char *)
		mycalloc (nr->partsize + 1, sizeof (char));
    (*aicvec)[*aicnum].aic = aic;
    (*aicvec)[*aicnum].mle = nr->llike;
    (*aicvec)[*aicnum].numparam = nr->numpop2 - numparam - 1 + remainnum;
    memcpy ((*aicvec)[*aicnum].pattern, nr->world->options->custm2, sizeof (char) * nr->partsize);
	
    (*aicvec)[*aicnum].lrt = lrt;
    (*aicvec)[*aicnum].prob = probchi (numparam, (*aicvec)[*aicnum].lrt);
    (*aicvec)[*aicnum].probcorr = probchiboundary ((*aicvec)[*aicnum].lrt, numparam, numparam);
}

/// print AIC model selection results
void
print_aicfile(aic_fmt **aicvec, long *aicnum, nr_fmt *nr)
{
    long i;
    if (nr->world->options->aicfile)
    {
        FPRINTF (nr->world->options->aicfile, "%f %f %li %f  %f ",
                 (*aicvec)[*aicnum].aic, (*aicvec)[*aicnum].lrt,
                 (*aicvec)[*aicnum].numparam,
                 (*aicvec)[*aicnum].prob,
                 (*aicvec)[*aicnum].probcorr);
		
        for (i = 0; i < nr->partsize; i++)
            FPRINTF (nr->world->options->aicfile, "%f ", nr->world->param0[i]);
        FPRINTF (nr->world->options->aicfile, "\n");
    }
}

/// print progress report for AIC model selection
void
print_progress_aic(boolean add
				   , aic_fmt **aicvec, long *aicnum, nr_fmt *nr, long numparam, long freeparam)
{
    MYREAL mle, aic, lrt, prob, probcorr;
    char *temppattern;
    char *custm2 = nr->world->options->custm2;
    temppattern =  (char *) mycalloc (nr->world->numpop2 + 1 + nr->world->numpop, sizeof (char));
    if(add
	   )
    {
        if (nr->world->options->progress)
            FPRINTF (stdout,
                     "           +   %s %20.5f %3li %8.4f %8.4f %6.4f %6.4f\n",
                     reshuffle (temppattern, custm2, nr->numpop), (*aicvec)[*aicnum].aic,
                     freeparam, (*aicvec)[*aicnum].mle,
                     (*aicvec)[*aicnum].lrt, (*aicvec)[*aicnum].prob,
                     (*aicvec)[*aicnum].probcorr);
        if (nr->world->options->writelog)
            FPRINTF (nr->world->options->logfile,
                     "           +   %s %20.5f %3li %8.4f %8.4f %6.4f %6.4f\n",
                     reshuffle (temppattern, custm2, nr->numpop), (*aicvec)[*aicnum].aic,
                     freeparam, (*aicvec)[*aicnum].mle,
                     (*aicvec)[*aicnum].lrt, (*aicvec)[*aicnum].prob,
                     (*aicvec)[*aicnum].probcorr);
    }
    else
    {
        mle = nr->llike;
        aic =  -2. * mle + 2. * freeparam;
        lrt =  -2. * (mle - (*aicvec)[0].mle);
        prob =  probchi (numparam, lrt);
        probcorr =  probchiboundary (lrt, numparam, numparam);
        if (nr->world->options->progress)
            FPRINTF (stdout,
                     "           -   %s %20.5f %3li %8.4f %8.4f %6.4f %6.4f\n",
                     reshuffle (temppattern, custm2, nr->numpop), aic,
                     freeparam, mle, lrt,
                     prob,
                     probcorr);
        if (nr->world->options->writelog)
            FPRINTF (nr->world->options->logfile,
                     "           +   %s %20.5f %3li %8.4f %8.4f %6.4f %6.4f\n",
                     reshuffle (temppattern, custm2, nr->numpop), aic,
                     freeparam, mle, lrt,
                     prob,
                     probcorr);
    }
    (void) fflush (stdout);
    myfree(temppattern);
}

/// add the AIC models and scores to the result list
void
add_aiclist (int migtype, long numparam, long freeparam, long remainnum, long zero, long m,
             MYREAL *param0, char *temppattern, char *custm2,
             long *aicnum, aic_fmt ** aicvec, nr_fmt * nr)
{
    MYREAL aic = -2. * nr->llike + 2. * freeparam;
	
    if (aic < (*aicvec)[0].aic + nr->world->options->aicmod * freeparam)
    {
        if (migtype != 'm' || (migtype == 'm' && nr->world->options->mmn > 1))
        {
            add_aicvec(aic, aicvec, aicnum, nr, numparam, remainnum);
            print_aicfile(aicvec,aicnum,nr);
            print_progress_aic(TRUE, aicvec,aicnum,nr, numparam, freeparam);
            (*aicnum)++;
            aic_score (aicvec, aicnum, nr, zero + 1, m + 1,
                       temppattern, param0, migtype);
        }
    }
    else
    {
        print_progress_aic(FALSE, aicvec,aicnum,nr, numparam, freeparam);
    }
}

/// calculates AIC scores
void
aic_score (aic_fmt ** aicvec, long *aicnum, nr_fmt * nr,
           long zero, long which, char *temppattern, MYREAL *param0,
           int migtype)
{
    long m;
    MYREAL likes = 0;
    MYREAL normd = 0;
    char savecustm2;
    long remainnum = 0;
    //    boolean mylegal;
    char *custm2 = nr->world->options->custm2;
    long numparam = 0;
    long freeparam;
	
    switch (migtype)
    {
		case '0':
			numparam = nr->world->options->zeron;
			remainnum = 0;
			break;
		case 'm':
			numparam = nr->world->options->mmn;
			remainnum = 1;
			if (nr->world->options->custm2[which] == 'm')
				return;
				break;
    }
    freeparam = (nr->numpop2 - numparam - 1 + remainnum);
    for (m = which; m < nr->numpop2; m++)
    {
        savecustm2 = custm2[m];
        custm2[m] = migtype;
        memcpy (nr->world->param0, param0, sizeof (MYREAL) * nr->numpop2);
        resynchronize_param (nr->world);
		
        if ((legal_pattern (nr->world->options->custm2, nr->numpop)))
        {
            (void) do_profiles (nr->world, nr, &likes, &normd, PROFILE,
                         nr->world->rep, nr->world->repkind);
			
            add_aiclist(migtype, numparam, freeparam, remainnum, zero, m, param0,
                        temppattern, nr->world->options->custm2, aicnum, aicvec, nr);
        }
        else
        {
            if (nr->world->options->progress)
            {
                FPRINTF (stdout, "           F   %s %20s\n", reshuffle (temppattern, custm2, nr->numpop), "-----");
                fflush (stdout);
            }
            if (nr->world->options->writelog)
                FPRINTF (nr->world->options->logfile, "           F   %s %20s\n", reshuffle (temppattern, custm2, nr->numpop), "-----");
        }
        custm2[m] = savecustm2;
    }
}

/* not used
/// check how many parameter are set to zero or are set to m (average)
static boolean
check_numparam (long which, long migtype, worldoption_fmt * options,
                long *numparam, long *remainnum)
{
    boolean rc = FALSE;
    switch (migtype)
    {
		case '0':
			*numparam = options->zeron;
			*remainnum = 0;
			break;
		case 'm':
			*numparam = options->mmn;
			*remainnum = 1;
			if (options->custm2[which] == 'm')
				rc = TRUE;
				break;
    }
    return rc;
}
end not used */

///
/// check whether the migration pattern in matrix is legal or not
/// a legal migration pattern, needs to connect all populations
/// eahc population needs to be connected to at least one other population
/// in one or both direction of the migration route  
boolean
legal_pattern (char *matrix, long numpop)
{
  //#ifdef MYREAL == float
  //const MYREAL eps = FLT_EPSILON ;
  //#else
  const MYREAL eps = DBL_EPSILON ;
  //#endif
    long from=0, to=0, i;
    MYREAL summ = -1;
    long oldto;
    for (i = 0; i < numpop; i++)
    {
        if (matrix[i] == '0')
            return FALSE;
    }
    oldto = -1;
    for (i = numpop; i < numpop * numpop; i++)
    {
        m2mm (i, numpop, &from, &to);
        if (oldto != to)
        {
	  if (fabs(summ)< eps )
                return FALSE;
            oldto = to;
            summ = 0;
        }
        summ += (MYREAL) (matrix[i] != '0') + (MYREAL) (matrix[mm2m (to, from, numpop)] != '0');
    }
    if (summ < eps)
        return FALSE;
    return TRUE;
}

long find_paramnum(world_fmt *world, char *connect)
{
    long cc;
    long tm;
    long mm;
    //    long mmm;
    long ms;
    long mss;
    long total;
    char *temp, *temp1, *temp2, *t;
    
    // checkout if GAMMA mutation rate
    total = world->numpop * world->numpop + (world->options->gamma ? 1 : 0) ;
    
    if(connect==NULL)
	{
        // checkout zeroes and constants in custom migration matrix
        // removes 1 for each zero and (constant>0)
        total -= world->options->zeron + world->options->constn;
        // checkout mean values
        // removes number of parameters that are part of the average and adds
        // 1 for average migrations and 1 for average thetas
        total -= (world->options->tmn > 0 ? world->options->tmn - 1 : 0) + (world->options->mmn > 0 ? world->options->mmn - 1 : 0); 
        // checkout symmetric values
        // removes 1 parameter for each migration pair
        total -= world->options->symn + world->options->sym2n;
	}
    else
	{
        temp1 = (char *) mycalloc(strlen(connect)+1, sizeof(char));
        //temp2 = (char *) mycalloc(strlen(connect)+1, sizeof(char));
        temp = (char *) mycalloc(strlen(connect)+1, sizeof(char));
        t = temp;
        strcpy(temp1,connect);
        temp2 = strtok(temp1," ;,:\n\t\r");
        while(temp2!= NULL)
		{
            switch(temp2[0])
			{
                case '*': *t++ = '*'; break;
                case 'c': *t++ = 'c'; break;
                case 'm': *t++ = 'm'; break;
                case 'M': *t++ = 'M'; break;
                case 's': *t++ = 's'; break;
                case 'S': *t++ = 'S'; break;
                case ' ' : break;
                default:
                    *t++ = 'c';
			}
            temp2 = strtok(NULL," ;,:\n\t\r");
		}
        cc = scan_connect (temp, 0, world->numpop2, 'c');
        tm = scan_connect (temp, 0, world->numpop, 'm');
        mm = scan_connect (temp, world->numpop, world->numpop2, 'm');
	//        mmm = scan_connect (temp, world->numpop, world->numpop2, 'M');
        ms = scan_connect (temp, world->numpop, world->numpop2, 's')/2;
        mss = scan_connect (temp, world->numpop, world->numpop2, 'S')/2;
        total -= cc + (tm > 0 ? tm-1 : 0);
        total -= (mm > 0 ? mm-1 : 0 );
        total -= ms;
        total -= mss;
        myfree(temp);
        myfree(temp1);
        //myfree(temp2);
	}
    return total;
}
