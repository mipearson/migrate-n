/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
 F S T   R O U T I N E S
 
 calculates FST
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
 Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
$Id: fst.c 1807 2011-03-18 20:16:19Z beerli $
 
-------------------------------------------------------*/
/*! \file fst.c

*/

#include "migration.h"
#include "tools.h"
#include "sighandler.h"
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

/* prototypes ------------------------------------------- */
void fst_type (int type);

void calc_fst (world_fmt * world, data_fmt * data);
/* private functions */
void frequencies (MYREAL ***f, char *****data, long numpop,
                  long **numind, long loci);
void calc_fw (MYREAL ***f, long numpop, long locus, MYREAL *fw);
void calc_fb (MYREAL ***f, long numpop, long locus, MYREAL *fb);
void calc_seq_fw (data_fmt * data, long numpop, long locus, MYREAL *fw);
void calc_seq_fb (data_fmt * data, long numpop, long locus, MYREAL *fb);
void solveNm_vartheta (MYREAL *fw, MYREAL *fb, long numpop, MYREAL *params);
void solveNm_varm (MYREAL *fw, MYREAL *fb, long numpop, MYREAL *params);
void solveNm_sym (MYREAL *fw, MYREAL *fb, long numpop, MYREAL *params);

/* global variable SOLVENM points to function solveNm_varxx() */
static void (*solveNm) (MYREAL *, MYREAL *, long, MYREAL *);

/*=======================================================*/
void
fst_type (int type)
{
    switch (type)
    {
    case 'T':
        solveNm =
            (void (*)(MYREAL *, MYREAL *, long, MYREAL *)) solveNm_vartheta;
        break;
    case 'M':
        solveNm = (void (*)(MYREAL *, MYREAL *, long, MYREAL *)) solveNm_varm;
        break;
    case 'S':
    default:
        solveNm = (void (*)(MYREAL *, MYREAL *, long, MYREAL *)) solveNm_sym;
        break;
    }
}

void
calc_fst (world_fmt * world, data_fmt * data)
{
    long pop, locus;
    long connections = world->numpop * (world->numpop - 1) / 2;
    MYREAL ***fstfreq,  *fw, *fb, *sumfb, *sumfw;
    if (connections == 0)
        return;
    fw = (MYREAL *) mycalloc (1, sizeof (MYREAL) * world->numpop);
    fb = (MYREAL *) mycalloc (1, sizeof (MYREAL) * connections);
    sumfw = (MYREAL *) mycalloc (1, sizeof (MYREAL) * world->numpop);
    sumfb = (MYREAL *) mycalloc (1, sizeof (MYREAL) * connections);
    if (world->numpop > 2)
    {
        fst_type ('S');
    }
 
    fstfreq = (MYREAL ***) mycalloc (world->numpop, sizeof (MYREAL **));
    for (pop = 0; pop < world->numpop; pop++)
    {
        doublevec2d(&fstfreq[pop],world->loci,255);  
        //fstfreq[pop] =
        //     (MYREAL **) mycalloc (world->loci, sizeof (MYREAL *));
        //for (locus = 0; locus < world->loci; locus++)
        //   fstfreq[pop][locus] =
        //      (MYREAL *) mycalloc (255, sizeof (MYREAL));
    }

    if (!strchr (SEQUENCETYPES, world->options->datatype))
    {
         frequencies (fstfreq, data->yy, world->numpop, data->numind,
                     world->loci);
    }
    for (locus = 0; locus < world->loci; locus++)
    {
        if (strchr (SEQUENCETYPES, world->options->datatype))
        {
            calc_seq_fw (data, world->numpop, locus, fw);
            calc_seq_fb (data, world->numpop, locus, fb);
        }
        else
        {
            calc_fw (fstfreq, world->numpop, locus, fw);
            calc_fb (fstfreq, world->numpop, locus, fb);
        }

        (*solveNm) (fw, fb, world->numpop, world->fstparam[locus]);

        for (pop = 0; pop < world->numpop; pop++)
        {
            sumfw[pop] += fw[pop];
        }
        for (pop = 0; pop < connections; pop++)
        {
            sumfb[pop] += fb[pop];
        }
    }
    for (pop = 0; pop < world->numpop; pop++)
    {
        sumfw[pop] /= world->loci;
    }
    for (pop = 0; pop < connections; pop++)
    {
        sumfb[pop] /= world->loci;
    }
    if(world->loci > 1)
      (*solveNm) (sumfw, sumfb, world->numpop, world->fstparam[world->loci]);
    if (!strchr (SEQUENCETYPES, world->options->datatype))
    {
      for (pop = world->numpop; pop >= 0; pop++)
        {
          free_doublevec2d(fstfreq[pop]);
//          for (locus = world->loci-1; locus >= 0; locus--)
  //              myfree(fstfreq[pop][locus]);
    //        myfree(fstfreq[pop]);
        }
    }
    myfree(fstfreq);

    myfree(fw);
    myfree(fb);
    myfree(sumfw);
    myfree(sumfb);
}


/*=======================================================*/

void
frequencies (MYREAL ***f,  char *****data, long numpop,
             long **numind, long loci)
{
    long **buckets;
    long pop, locus, ind, a;
    long total=0;
    buckets = (long **) mymalloc (sizeof (long *) * numpop);
    for (pop = 0; pop < numpop; pop++)
    {
        buckets[pop] = (long *) mycalloc (255, sizeof (long));
        for (locus = 0; locus < loci; locus++)
        {
            memset (buckets[pop], 0, sizeof (long) * 255);
            total = 0;
            for (ind = 0; ind < numind[pop][locus]; ind++)
            {
                if (data[pop][ind][locus][0][0] != '?')
                {
                    buckets[pop][data[pop][ind][locus][0][0] - '!'] += 1;
                    total += 1;
                }
                if (data[pop][ind][locus][1][0] != '?')
                {
                    buckets[pop][data[pop][ind][locus][1][0] - '!'] += 1;
                    total += 1;
                }
            }
            for (a = 0; a < 255; a++)
            {
                if (total > 0)
                    f[pop][locus][a] = (MYREAL) buckets[pop][a] / (MYREAL) total;
            }
        }
        myfree(buckets[pop]);
    }
    myfree(buckets);
}



void
calc_fw (MYREAL ***f, long numpop, long locus, MYREAL *fw)
{
    long pop, i;
    for (pop = 0; pop < numpop; pop++)
    {
        fw[pop] = 0;
        for (i = 0; i < 255; i++)
        {
            fw[pop] += f[pop][locus][i] * f[pop][locus][i];
        }
    }
}
void
calc_fb (MYREAL ***f, long numpop, long locus, MYREAL *fb)
{
    long i, p1, p2, zz = 0;
    for (p1 = 0; p1 < numpop; p1++)
    {
        for (p2 = p1 + 1; p2 < numpop; p2++)
        {
            fb[zz] = 0.0;
            for (i = 0; i < 255; i++)
            {
                fb[zz] += f[p1][locus][i] * f[p2][locus][i];
            }
            zz++;
        }
    }
}

void
calc_seq_fw (data_fmt * data, long numpop, long locus, MYREAL *fw)
{
    long pop, i, k, j;
    MYREAL nn;
    MYREAL diff;
    for (pop = 0; pop < numpop; pop++)
    {
        fw[pop] = 0;
        nn =
            data->seq[0]->sites[locus] * (data->numind[pop][locus] *
                                       data->numind[pop][locus] -
                                       data->numind[pop][locus]) / 2.;
        for (i = 0; i < data->numind[pop][locus]; i++)
        {
            for (k = i + 1; k < data->numind[pop][locus]; k++)
            {
                diff = 0.;
                for (j = 0; j < data->seq[0]->sites[locus]; j++)
                {
                    diff +=
                        (data->yy[pop][i][locus][0][j] !=
                         data->yy[pop][k][locus][0][j]);
                }
                if (nn > 0)
                    fw[pop] += diff / nn;
            }
        }
        fw[pop] = 1. - fw[pop];
    }
}

void
calc_seq_fb (data_fmt * data, long numpop, long locus, MYREAL *fb)
{
    long i, k, j;
    MYREAL nn, temp;
    MYREAL diff;
    long p1, p2, zz = 0;
    for (p1 = 0; p1 < numpop; p1++)
    {
        for (p2 = p1 + 1; p2 < numpop; p2++)
        {
            temp = 0.;
            nn =
                (MYREAL) data->seq[0]->sites[locus] * data->numind[p1][locus] *
                data->numind[p2][locus];
            for (i = 0; i < data->numind[p1][locus]; i++)
            {
                for (k = 0; k < data->numind[p2][locus]; k++)
                {
                    diff = 0.;
                    for (j = 0; j < data->seq[0]->sites[locus]; j++)
                    {
                        diff +=
                            (data->yy[p1][i][locus][0][j] !=
                             data->yy[p2][k][locus][0][j]);
                    }
                    if (nn > 0)
                        temp += diff / nn;
                }
            }
            fb[zz++] = 1. - temp;
        }
    }
}



void
solveNm_varm (MYREAL *fw, MYREAL *fb, long numpop, MYREAL *params)
{
    long i, p1 /*,p2 */ ;
    /*Version 2.0  MYREAL sumfw = sum(fw,numpop);
    MYREAL sumfb = sum(fb,numpop*(numpop-1)/2); */
    MYREAL first, denom;
    long offset2 = numpop + numpop * (numpop - 1);
    long offset = numpop;
    long numfb = 0;
    denom = (2. * fb[0] + fw[0] + fw[1]);
    if (denom == 0.0)
        first = -999.;
    else
        first = (2. - fw[0] - fw[1]) / denom;

    for (p1 = 0; p1 < numpop; p1++)
    {
        numfb += p1;
        params[p1] = first;
        params[offset2 + p1] = fw[p1];
    }
    denom = ((fb[0] - fw[0]) * (-2. + fw[0] + fw[1]));
    if (denom == 0.0)
        params[offset] = -999;
    else
        params[offset] =
            (2. * fb[0] - fw[0] - 2. * fb[0] * fw[0] + fw[1]) / denom;
    denom = ((fb[0] - fw[1]) * (-2. + fw[0] + fw[1]));
    if (denom == 0.0)
        params[offset + 1] = -999;
    else
        params[offset + 1] =
            (2. * fb[0] - fw[1] - 2. * fb[0] * fw[1] + fw[0]) / denom;
    for (i = 0; i < offset2; i++)
    {
        if (params[i] < 0.)
            params[i] = -999;
    }
    for (i = 0; i < numfb; i++)
    {
        params[offset2 + numpop + i] = fb[i];
    }
}

void
solveNm_vartheta (MYREAL *fw, MYREAL *fb, long numpop, MYREAL *params)
{
    long i;
    MYREAL nom;
    MYREAL denom, fulldenom;
    long offset2 = numpop + numpop * (numpop - 1);
    long numfb = 0;
    nom = (-2. * fb[0] + fw[0] + fw[1]);
    for (i = 0; i < numpop; i++)
    {
        numfb += i;
    }
    denom = -2. * fb[0] * fb[0] + fw[0] * fw[1];

    for (i = 0; i < numpop; i++)
    {
        params[offset2 + i] = fw[i];
        fulldenom = denom + fw[i] * fw[i];
        if (fulldenom > 0)
            params[i] = (nom * (1. - fw[i])) / fulldenom;
        else
            params[i] = -999;
    }
    if (nom == 0.0)
        params[numpop] = -999;
    else
        params[numpop] = 2. * fb[0] / nom;

    for (i = 1; i < numpop * (numpop - 1); i++)
    {
        params[numpop + i] = params[numpop];
    }
    for (i = 0; i < offset2; i++)
    {
        if (params[i] < 0.)
            params[i] = -999;
    }
    for (i = 0; i < numfb; i++)
    {
        params[offset2 + numpop + i] = fb[i];
    }
}

void
solveNm_sym (MYREAL *fw, MYREAL *fb, long numpop, MYREAL *params)
{
    long p1, p2;
    MYREAL nom;
    MYREAL denom, fulldenom;
    MYREAL fbji;
    MYREAL **xxx;
    long offset2 = numpop + numpop * (numpop - 1);
    long z = 0;
    xxx = (MYREAL **) mycalloc (1, sizeof (MYREAL *) * numpop);
    xxx[0] = (MYREAL *) mycalloc (1, sizeof (MYREAL) * numpop * numpop);
    for (p1 = 1; p1 < numpop; p1++)
    {
        xxx[p1] = xxx[0] + numpop * p1;
    }

    for (p1 = 0; p1 < numpop; p1++)
    {
        params[offset2 + p1] = fw[p1];
        params[p1] = 0.;
        for (p2 = p1 + 1; p2 < numpop; p2++)
        {
            fbji = fb[z++];
            params[offset2 + numpop + z - 1] = fbji;
            nom = (-2. * fbji + fw[p1] + fw[p2]);
            denom = -2. * fbji * fbji + fw[p1] * fw[p2];
            fulldenom = (denom + fw[p1] * fw[p1]);
            if (fulldenom > 0.0 || fulldenom < 0.0)
                params[p1] += (nom * (1. - fw[p1])) / fulldenom;
            if (!(nom < 0.0 || nom > 0.0))
                xxx[p2][p1] = -999. ;
            else
                xxx[p2][p1] = xxx[p1][p2] = 2 * fbji / nom;
        }
        params[p1] /= numpop;
    }
    z = 0;
    for (p1 = 0; p1 < numpop; p1++)
    {
        for (p2 = 0; p2 < numpop; p2++)
        {
            if (p1 != p2)
                params[numpop + z++] = xxx[p1][p2];
        }
    }
    for (p1 = 0; p1 < offset2; p1++)
    {
        if (params[p1] < ((MYREAL) 0.0) )
            params[p1] = -999. ;
    }
    myfree(xxx[0]);
    myfree(xxx);
}
