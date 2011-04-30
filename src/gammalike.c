/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
 Gammalike and derivatives  R O U T I N E S
 
 Peter Beerli 2001, Seattle
 beerli@fsu.edu
 
 Copyright 2001-2003 Peter Beerli and Joseph Felsenstein
 Copyright 2004-2005 Peter Beerli

 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
$Id: gammalike.c 1800 2011-01-29 13:40:00Z beerli $
*/
/*! \file gammalike.c 

Calculates the gamma distributed rate variation among loci
and derivatives
*/


#include "migration.h"
#include "tools.h"
#include "broyden.h"
#include "combroyden.h"
#include "laguerre.h"
#include "sighandler.h"

MYREAL gamma_loci_like (helper_fmt *helper, MYREAL *oparam,
                        MYREAL *olparam, MYREAL denom);
MYREAL gamma_locus_like (nr_fmt * nr, MYREAL *oparam,
                         MYREAL *olparam, MYREAL denom, long locus);
MYREAL loggf (MYREAL r, MYREAL alpha, MYREAL theta1);
void gamma_loci_derivative (helper_fmt * helper);
void gamma_locus_derivative (helper_fmt * helper, long locus);
//MYREAL gamma_theta1_derivative (MYREAL theta, MYREAL *param, nr_fmt * nr,
//                              long locus, int *sign);
//MYREAL gamma_other_derivative (long which, MYREAL theta, MYREAL *param,
//                             nr_fmt * nr, long locus, int *sign);
//MYREAL gamma_alpha_derivative (MYREAL theta, MYREAL *param, nr_fmt * nr,
//                             long locus, int *sign);
void gamma_parts_gradient (MYREAL x, MYREAL *param, nr_fmt * nr, long locus);
MYREAL dt_locus_gsum (long which, MYREAL x, MYREAL *param, nr_fmt * nr,
                      timearchive_fmt ** atl, long locus);


MYREAL
gamma_loci_like(helper_fmt *helper, MYREAL *oparam, MYREAL *olparam, MYREAL denom)
{
    long locus;
    MYREAL summ = 0.0;
    nr_fmt *nr = helper->nr;
    for (locus = 0; locus < nr->world->loci; locus++)
    {
        if (nr->skiploci[locus])
            continue;
        nr->locilikes[locus] =
            gamma_locus_like (nr, oparam, olparam, denom, locus);
        summ += nr->locilikes[locus];
    }
    return summ;
}


MYREAL find_gamma_summ(MYREAL la, MYREAL theta1, nr_fmt *nr, MYREAL *oparam, MYREAL *olparam,
                       long locus, helper_fmt *helper);


MYREAL
gamma_loci_like_iterative(helper_fmt * helper, MYREAL *oparam, MYREAL *olparam, MYREAL denom)
{
    // hack to get around the problem of not being able to get a clean derivative of alpha
    // uses a simple maxima finder for alpha
    long locus=helper->nr->world->locus;
    MYREAL fa = 0.0;
    MYREAL fla, fha;
    MYREAL la  = -10;
    MYREAL ha = +4;
    MYREAL a;
    MYREAL fm, m;
    MYREAL theta1 = helper->nr->param[0];
    nr_fmt *nr = helper->nr;
    char save_a = nr->world->options->custm[nr->world->numpop2];
    //  if(save_a == 'c')
    // {
    // return  gamma_loci_like_old (helper, oparam, olparam, denom);
    // }
    nr->world->options->custm[nr->world->numpop2]='c';
    nr->world->options->custm2[nr->world->numpop2]='c';
    fla = find_gamma_summ( la,  theta1, nr,  oparam, olparam,  locus, helper);
    fha = find_gamma_summ( ha,  theta1, nr,  oparam, olparam,  locus, helper);
    a = (ha + la)/2.;
    fa = find_gamma_summ( a,  theta1, nr,  oparam, olparam,  locus, helper);
    m = (la + a) / 2.;
    fm = find_gamma_summ( m,  theta1, nr,  oparam, olparam,  locus, helper);
    while(fabs(fla-fha)> EPSILON && fabs(ha-la) > EPSILON)
    {
        if(fm < fa)
        {
            fla = fm;
            la = m;
        }
        else
        {
            fha = fa;
            ha  = a;
            a = m;
           //xcode  fa = fm;
        }
        m = (ha + la)/2.;
        fa = find_gamma_summ( m,  theta1, nr,  oparam, olparam,  locus, helper);

    }
    nr->param[nr->numpop2] = EXP(m);
    nr->lparam[nr->numpop2] = m;
    nr->world->options->custm[nr->world->numpop2]=save_a;
    nr->world->options->custm2[nr->world->numpop2]=save_a;
    return fa;
}


MYREAL
gamma_loci_like_new_notworking(helper_fmt * helper, MYREAL *oparam, MYREAL *olparam, MYREAL denom)
{
    long locus=0;
    long trials=0;
    MYREAL logalpha= -MYREAL_MAX;
    MYREAL psiabc;
    MYREAL theta1 = helper->nr->param[0];
    nr_fmt *nr = helper->nr;
    char save_a = nr->world->options->custm[nr->world->numpop2];

    MYREAL a, b, c,  m=0.5;
    MYREAL la, lb, lc, ll=0, ql;
    //MYREAL lm;
    a= 0.00000000001;
    b= 1.0;
    c=2.0;
    nr->world->options->custm[nr->world->numpop2]='c';
    nr->world->options->custm2[nr->world->numpop2]='c';
    la = find_gamma_summ( a,  theta1, nr,  oparam, olparam,  locus, helper);
    lb = find_gamma_summ( b,  theta1, nr,  oparam, olparam,  locus, helper);
    lc = find_gamma_summ( c,  theta1, nr,  oparam, olparam,  locus, helper);

    while(trials++ < NTRIALS)
    {
        logalpha =
            (c * c * (lb - la) + b * b * (la - lc) +
             a * a * (lc - lb)) / (2. * (b * la - c * la - a * lb + c * lb +
                                         a * lc - b * lc));
        psiabc = ((la - lb) / (-a + b) + (la - lc) / (a - c)) / (-b + c);
        if ((psiabc <= 0.0) || (logalpha >= m))
        {
            if (a == m)
            {
                ll = la;
                logalpha = m;
                break;
            }
            if (b == m)
            {
                ll = lb;
                logalpha = m;
                break;
            }
            if (c == m)
            {
                ll = lc;
                logalpha = m;
                break;
            }
            ll = find_gamma_summ( m,  theta1, nr,  oparam, olparam,  locus, helper);
            replace_with ((int) m, &a, &b, &c, m, &la, &lb, &lc, ll);
        }
        ll = find_gamma_summ( logalpha,  theta1, nr,  oparam, olparam,  locus, helper);
        ql = quadratic_lamda (logalpha, a, b, c, la, lb, lc);
        if ((fabs (ll - MYMIN3 (la, lb, lc)) <= BIGEPSILON)
                || (fabs (ll - ql) <= BIGEPSILON))
        {
            break;
        }
        else
        {
            replace_with (1, &a, &b, &c, logalpha, &la, &lb, &lc,
                          ll);
            m = MYMAX3 (a, b, c);
        }
    }
    nr->param[nr->numpop2] = EXP(logalpha);
    nr->lparam[nr->numpop2] = logalpha;
    nr->world->options->custm[nr->world->numpop2]=save_a;
    nr->world->options->custm2[nr->world->numpop2]=save_a;
    return ll;
}

MYREAL find_gamma_summ(MYREAL lla, MYREAL theta1, nr_fmt *nr, MYREAL *oparam, MYREAL *olparam,
                       long locus, helper_fmt *helper)
{
    MYREAL la = EXP(lla);
    MYREAL denom =  -LGAMMA(la) - la * log(theta1/la);
    MYREAL summ=0.0;
    MYREAL so = oparam[nr->numpop2];
    MYREAL slo = olparam[nr->numpop2];
    oparam[nr->numpop2] = la;
    olparam[nr->numpop2] = lla;
    if (olparam[nr->numpop2] > 9.903487553)
    {
        olparam[nr->numpop2] = 9.903487553;
    }
    initgammacat (nr->categs, oparam[nr->numpop2],1./* EXP (lparam[0])*/,
                  nr->rate, nr->probcat);
    //calc_loci_param (nr, olparam, oparam, helper->dv, helper->lamda, nr->partsize);
    for (locus = 0; locus < nr->world->loci; locus++)
    {
        if (nr->skiploci[locus])
            continue;
        nr->locilikes[locus] = gamma_locus_like (nr, oparam, olparam, denom, locus);
        summ += nr->locilikes[locus];
    }
    oparam[nr->numpop2] = so;
    olparam[nr->numpop2] = slo;
    return summ;
}

MYREAL
gamma_locus_like (nr_fmt * nr, MYREAL *oparam, MYREAL *olparam, MYREAL denom,
                  long locus)
{
    //  nr_fmt *nr = helper->nr;
    long r;
    MYREAL tempmax = -MYREAL_MAX;
    MYREAL *temp;
    MYREAL summ = 0.;
    MYREAL alpha = oparam[nr->numpop2];
    MYREAL theta1 = oparam[0];
    temp = (MYREAL *) mycalloc (GAMMA_INTERVALS, sizeof (MYREAL));
    for (r = 0; r < GAMMA_INTERVALS; r++)
    {
        set_gamma_param (nr->param, oparam,
                         nr->lparam, olparam, nr->rate[r], nr);
        temp[r] = nr->probcat[r] + loggf (nr->rate[r], alpha, theta1) - denom +
                  calc_locus_like (nr, nr->param, nr->lparam, locus);
        //      printf("%f ",temp[r]);
        if (temp[r] > tempmax)
            tempmax = temp[r];
    }
    //  printf("\n");
    for (r = 0; r < GAMMA_INTERVALS; r++)
        summ += EXP (temp[r] - tempmax);
    myfree(temp);
    return LOG (summ) + tempmax;
}

MYREAL
loggf (MYREAL r, MYREAL alpha, MYREAL theta1)
{
    return -r * alpha / theta1 + LOG (r) * (alpha - 1.);
}

void
gamma_loci_difference (helper_fmt * helper)
{
    long i;
    nr_fmt *nr = helper->nr;
    long size = nr->partsize;
    MYREAL mle, high, low, param1;
    MYREAL *lparam;
    MYREAL *param;
    MYREAL *xv = helper->xv;
    MYREAL *expxv = helper->expxv;
    MYREAL epsilon = EPSILON;
    //  MYREAL logepsilon = log(epsilon);
    memset (nr->d, 0, size * sizeof (MYREAL));
    setdoublevec1d (&lparam, helper->xv, size);
    setdoublevec1d (&param, helper->expxv, size);
    mle = gamma_loci_like (helper, param, lparam, helper->weight);
    for (i = 0; i < size; i++)
    {
        memcpy (param, expxv, sizeof (MYREAL) * size);
        memcpy (lparam, xv, sizeof (MYREAL) * size);
        param1 = (param[i] += epsilon);
        lparam[i] = LOG (param[i]);
        high = gamma_loci_like (helper, param, lparam, helper->weight);
        param[i] -= 2. * epsilon;
        if (param[i] < 0.0)
        {
            low = mle;
            param[i] = expxv[i];
        }
        else
        {
            lparam[i] = LOG (param[i]);
            low = gamma_loci_like (helper, param, lparam, helper->weight);
        }
        nr->d[i] = (high - low) / (param1 - param[i]);
    }
    myfree(lparam);
    myfree(param);
}

void
gamma_loci_derivative (helper_fmt * helper)
{
    long locus;
    memset (helper->nr->d, 0, helper->nr->partsize * sizeof (MYREAL));
    for (locus = 0; locus < helper->nr->world->loci; locus++)
    {
        if (helper->nr->skiploci[locus])
            continue;
        //adds up in nr->d
        gamma_locus_derivative (helper, locus);
    }
}

void
gamma_locus_derivative (helper_fmt * helper, long locus)
{
    long r, i, ii, z;
    MYREAL **temp, *tempmax;
    MYREAL ll;
    MYREAL sumi;
    int **stemp;
    //MYREAL *summll;
    //MYREAL summllsum;
    //MYREAL summllmax;
    MYREAL rate;
    MYREAL lograte;
    MYREAL tmp;
    MYREAL denom = helper->weight;
    nr_fmt *nr = helper->nr;
    MYREAL alpha = helper->expxv[nr->numpop2];
    MYREAL theta1 = helper->expxv[0];
    long partsize = nr->partsize;
    //summll = (MYREAL *) mycalloc(GAMMA_INTERVALS, sizeof(MYREAL));
    tempmax = (MYREAL *) mycalloc (partsize, sizeof (MYREAL));
    temp = (MYREAL **) mycalloc (partsize, sizeof (MYREAL *));
    temp[0] = (MYREAL *) mycalloc (partsize * GAMMA_INTERVALS, sizeof (MYREAL));
    stemp = (int **) mycalloc (partsize, sizeof (int *));
    stemp[0] = (int *) mycalloc (partsize * GAMMA_INTERVALS, sizeof (int));
    tempmax[0] = -MYREAL_MAX;
    for (r = 1; r < partsize; r++)
    {
        temp[r] = temp[0] + r * GAMMA_INTERVALS;
        stemp[r] = stemp[0] + r * GAMMA_INTERVALS;
        tempmax[r] = -MYREAL_MAX;
    }
    //summllmax = -MYREAL_MAX;
    //summllsum=0.0;
    for (r = 0; r < GAMMA_INTERVALS; r++)
    {
        rate = nr->rate[r];
        lograte = LOG (rate);
        set_gamma_param (nr->param, helper->expxv,
                         nr->lparam, helper->xv, rate, nr);
        ll = nr->probcat[r] + calc_locus_like (nr, nr->param, nr->lparam, locus);
        ll +=  -rate * alpha / theta1 + lograte * (alpha - 1.);
        //ll +=  -rate * alpha + lograte * (alpha - 1.);
        gamma_parts_gradient (rate, helper->expxv, nr, locus); //changes nr->parts
        //first of the derivatives: theta1
        //     tmp =  nr->parts[0];
        tmp = rate * alpha / (theta1 * theta1) + nr->parts[0];

        stemp[0][r] = tmp < 0. ? -1 : 1;
        temp[0][r] = ll + LOG (fabs (tmp));
        if (tempmax[0] < temp[0][r])
            tempmax[0] = temp[0][r];
        //all others except alpha
        for (i = 1; i < nr->numpop2; i++)
        {
            stemp[i][r] = (nr->parts[i] < 0.) ? -1 : 1;
            temp[i][r] = ll + LOG (fabs (nr->parts[i]));
            if (tempmax[i] < temp[i][r])
                tempmax[i] = temp[i][r];
        }
        // last derivative: alpha
        tmp = -rate/theta1 + lograte;
        //tmp = - theta1*lograte + rate;
        stemp[i][r] = tmp < 0. ? -1 : 1;
        temp[i][r] = ll + LOG (fabs (tmp));
        if (tempmax[i] < temp[i][r])
            tempmax[i] = temp[i][r];
    }
    z = 0;
    for (ii = 0; ii < nr->partsize - nr->profilenum; ii++)
    {
        i = (nr->profilenum > 0) ? nr->indeks[z++] : ii;
        sumi = 0;
        //   denom = -LGAMMA(alpha) - alpha * log(1./alpha);
        denom = -LGAMMA(alpha) - alpha * log(theta1/alpha);
        for (r = 0; r < GAMMA_INTERVALS; r++)
            sumi +=
                ((MYREAL) stemp[i][r]) * EXP (temp[i][r] -
                                              tempmax[i]  +  denom);
        if (i == 0)
        {
            sumi =
                //     sumi * EXP (tempmax[i] - nr->locilikes[locus]);
                sumi * EXP (tempmax[i] - nr->locilikes[locus]) - alpha / theta1;
        }
        else
        {
            if (i == nr->numpop2)
            {
                sumi = sumi *(EXP (tempmax[i] - nr->locilikes[locus]))
                       - polygamma (0, alpha) - LOG (theta1 / alpha) +1.;
            }
            else
                sumi = sumi * EXP (tempmax[i] - nr->locilikes[locus]);
        }
        nr->d[ii] += sumi;
    }
    myfree(temp[0]);
    myfree(temp);
    myfree(stemp[0]);
    myfree(stemp);
    //myfree(summll);
}

void
gamma_parts_gradient (MYREAL x, MYREAL *param, nr_fmt * nr, long locus)
{
    long g;

    //  MYREAL copysum=0;
    long pop, i;
    MYREAL expapg;  //, summ;
    //MYREAL summm = 0.0;
    long msta, msto;
    MYREAL waitmig;
    MYREAL pointmig;
    tarchive_fmt *tl;
    MYREAL *part;   //, df = 0.0;
    //  MYREAL value ;
    MYREAL theta1 = param[0];
    MYREAL *geo = nr->world->data->geo;
    long rep;
    MYREAL theta_ratio = theta1 / x;
    part = (MYREAL *) mymalloc (sizeof (MYREAL) * nr->partsize);
    memset (nr->parts, 0, sizeof (MYREAL) * nr->partsize);
    for (rep = nr->repstart; rep < nr->repstop; rep++)
    {
        for (g = 0; g < nr->atl[rep][locus].T; g++)
        {
            if (nr->apg[rep][locus][g] > -40.)
            {
                tl = nr->atl[rep][locus].tl;
                //copies = (g > 0) ? tl[g].copies : tl[g].copies - 1;
                expapg = /*copies  */ EXP (nr->apg[rep][locus][g]);
                if (expapg == 0.0)
                    continue;
                //      summm += expapg;
                // population 1 is treated differently
                part[0] = 0.0;
                msta = nr->mstart[0];
                msto = nr->mend[0];
                for (i = msta; i < msto; i++)
                {
                    part[0] += -geo[i] * param[i] * tl[g].km[0] / x +
                               tl[g].mindex[i] / theta1;
                    if (param[i] > 0.)
                        part[i] = ((tl[g].mindex[i] / param[i])
                                   - geo[i] * tl[g].km[0] * theta_ratio);
                }
                // all other populations
                for (pop = 1; pop < nr->numpop; pop++)
                {
                    part[pop] = -tl[g].p[pop] / param[pop] +
                                tl[g].kt[pop] * theta1 / (x * param[pop] * param[pop]);

                    msta = nr->mstart[pop];
                    msto = nr->mend[pop];
                    //              z = 0;
                    waitmig = 0;
                    pointmig = 0;
                    for (i = msta; i < msto; i++)
                    {
                        waitmig += geo[i] * param[i];
                        pointmig += tl[g].mindex[i];
                        if (param[i] > 0.)
                            part[i] = ((tl[g].mindex[i] / param[i])
                                       - geo[i] * tl[g].km[pop] * theta_ratio);
                    }
                    part[0] += pointmig / theta1 - waitmig * tl[g].km[pop] / x +
                               tl[g].p[pop] / theta1 - tl[g].kt[pop] / (param[pop] * x);
                }
                for (pop = 0; pop < nr->world->numpop2; pop++)
                    nr->parts[pop] += expapg * part[pop];
            }
        }
    }
    for (pop = 0; pop < nr->world->numpop2; pop++)
        nr->parts[pop] /= nr->PGC[locus];
    // printf("nr->PGC[%li]=%f == %f\n",locus,nr->PGC[locus],summm);
    myfree(part);
}
