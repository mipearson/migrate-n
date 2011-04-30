/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 Likelihood ratio test   R O U T I N E S 
 
 moved out from world.c                                                                                                               
 Peter Beerli 2000, Seattle
 beerli@fsu.edu
 
 Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
$Id: lrt.c 1800 2011-01-29 13:40:00Z beerli $
 
-------------------------------------------------------*/
/* \file lrt.c 
Calculates likelihood ratio tests
*/
#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#include "broyden.h"
#include "combroyden.h"
#include "options.h"
#include "aic.h"
#include "migrate_mpi.h"
#include "pretty.h"

#ifndef LAGUERRE
//#include "derivatives2.h"
#endif
#include "sort.h"

//#ifdef UEP
//#include "uep.h"
//#endif
#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

void print_lratio_test (world_fmt * world, long *Gmax);
void test_loci_like (nr_fmt * nr, MYREAL *param0,
                     MYREAL *param1, long df, long zeros,
                     long loci, world_fmt * world,
                     long *maxwhich, long maxnum,
                     boolean withhead, char *this_string1);
long set_test_param (MYREAL *param, lr_data_fmt  *lrtdata, world_fmt * world,
                     long lrline, long locus, long *maxwhich,
                     long *maxnum, long *zeros);
void print_lrt_box(world_fmt *world, MYREAL *param0, MYREAL *param1,  long zeros, long elem,
                   char * this_string1, MYREAL like0, MYREAL like1, long df);

long parse_h0(char ***box, long *size, char *thisString, MYREAL *param0, MYREAL *param1 , long elem, world_fmt *world);
void parse_h0part(char ***box, long *size, long *tempsize,
                MYREAL *param, long elem, char head[], world_fmt *world);
void parse_h0string(char ***box, long *size, long *tempsize,
                    char *thisString);

void print_box(world_fmt* world,char **box,long size,
               MYREAL like0, MYREAL like1,MYREAL lrt, MYREAL chiprob, MYREAL chiprob2, MYREAL aic, long df, long aicparamnum);

long simplify_lrtvalues(char * in, char** out);

boolean
mywhitespace (char ch)
{
    if (!(ch == ',' || ch == '\0' || ch == '\t' || ch == '\r' || ch == '\n' || ch == ';'))
    {
        return TRUE;
    }
    return FALSE;
}


void
print_lratio_test (world_fmt * world, long *Gmax)
{
    long c;
    long r, locus;
    long df;
    long zeros;
    int header;
    nr_fmt *nr;
    MYREAL *param0;
    MYREAL *param1;
    long *maxwhich;
    long nparam;
    long maxnum = 0;
    long rep = !world->options->replicate ? 0 :
               (world->options->replicatenum == 0 ?
                world->options->lchains : world->options->replicatenum);
    long repstop = world->repstop;
    param0 = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (world->numpop2 + 1));
    param1 = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (world->numpop2 + 1));
    maxwhich = (long *) mycalloc (1, sizeof (long) * (world->numpop2 + 1));


    if (world->options->progress)
        FPRINTF (stdout, "           Printing likelihood ratio tests\n");

    nr = (nr_fmt *) mycalloc (1, sizeof (nr_fmt));
    create_nr (nr, world, *Gmax, 0, world->loci, world->repkind, world->rep);
    for (locus = 0; locus < world->loci; locus++)
    {
        if (world->options->replicate)
            for (r = 0; r < repstop; r++)
                create_multiapg0 (nr->apg0[r][locus], nr, r, locus);
        else
            create_apg0 (nr->apg0[0][locus], nr, &world->atl[0][locus], locus);
    }
    if (world->loci > 1)
    {
        locus = world->loci;
        rep = 0;
    }
    else
    {
        world->locus=0;
        locus = 0;
    }
    PAGEFEEDWORLD;
    nparam = world->options->gamma ? world->numpop2 + 1 : world->numpop2;
    for (c = 0; c < world->options->lratio->counter; c++)
    {
        header = (c == 0) ? HEADER : NOHEADER;
        if (world->options->lratio->data[c].type == MLE)
        {
            memcpy (param1, world->atl[rep][locus].param,
                    sizeof (MYREAL) * nparam);
            df = set_test_param (param0,
                                 &world->options->lratio->data[c],
                                 world, 0, -1, maxwhich, &maxnum, &zeros);
            test_loci_like (nr, param0, param1,
                            df, zeros, world->loci, world, maxwhich,
                            maxnum, header,
                            world->options->lratio->data[c].value1);
        }
    }
    fflush (world->outfile);
    myfree(param0);
    myfree(param1);
    myfree(maxwhich);
    destroy_nr (nr, world);
}


//remember: param0 is the parameterset to test and
// NOT the parameterset from migrate.
#define BOXSIZE 50 /*defines the print width of the left box containing the H0*/
#define BOXSIZE2 40 /*defines the print width of the left box with the legend*/
void
test_loci_like (nr_fmt * nr, MYREAL *param0, MYREAL *param1, long df,
                long zeros, long loci, world_fmt * world, long *maxwhich,
                long maxnum, boolean withhead, char *this_string1)
{
    
 //   char *teststat, temp[LRATIO_STRINGS];
    MYREAL like1, like0;//, testval, chiprob, chiprob2;
 //   int length;
 //  
    long numparam;
    long i, j, g = 0;
    long elem;
    //long zi;
    long z = 0, w = 0;
    // long pop;
    //MYREAL normd = 0.0;
    long *which;
    MYREAL **hess;
    MYREAL /**values,*/ *saveparam0;
    //MYREAL tparam;
//    long spaces = 0;
    char temp[100];
    helper_fmt helper;
    MYREAL *lparam0;
    MYREAL *lparam1;
    MYREAL aicfull;
    long aicfullparamnum;
//    char *message;
    numparam = world->numpop2 + 1;
    doublevec2d(&hess,numparam, numparam);
    lparam0 = (MYREAL *) mycalloc (numparam, sizeof (MYREAL));
    lparam1 = (MYREAL *) mycalloc (numparam, sizeof (MYREAL));
    which = (long *) mycalloc (1, sizeof (long) * numparam);
//    message = (char *) mycalloc (20, sizeof (char));
    //  values = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (world->numpop2 + 1));
    saveparam0 = (MYREAL *) mymalloc (sizeof (MYREAL) * numparam);
    memcpy (saveparam0, param0, sizeof (MYREAL) * numparam);
    for (i = 0; i < loci; i++)
    {
        for (j = 0; j < world->repstop; j++)
        {
            if (g < world->atl[j][i].T)
                g = world->atl[j][i].T;
        }
    }
    elem = world->options->gamma ? numparam : nr->numpop2;
    nr->skiploci = world->data->skiploci;
    helper.multilocus = world->loci == 1 ? FALSE : TRUE;
    if (maxnum > 0)
    {
        for (i = 0; i < elem; i++)
        {
            if (i != maxwhich[z])
            {
                which[w] = i;
                nr->values[w++] = param0[i];
            }
            else
            {
                if (maxnum > z + 1)
                    z++;
            }
        }
        nr->profilenum = w;
        maximize (&param0, world, nr, hess, PROFILE, world->repkind);
        like0 = nr->llike;
        //xcode normd = nr->normd;
        memcpy (param0, world->param0, sizeof (MYREAL) * nr->partsize);
    }
    else
    {
        set_logparam (lparam0, param0, elem);
        fill_helper (&helper, param0, lparam0, world, nr);
        like0 = CALCLIKE (&helper, param0, lparam0);
    }
    set_logparam (lparam1, param1, elem);
    fill_helper (&helper, param1, lparam1, world, nr);
    like1 = CALCLIKE (&helper, param1, lparam1);
//    sprintf (message, " %f ", chiprob2);

    if (withhead)
      {
        FPRINTF (world->outfile,
                 "==============================================================================\n");
        FPRINTF (world->outfile, "Likelihood ratio tests\n");
        FPRINTF (world->outfile,
                 "==============================================================================\n");
        FPRINTF (world->outfile, "Over all loci\n");
        FPRINTF (world->outfile, "Legend for the LRT tables\n");
        print_line(world->outfile,'-',79,CONT);
        sprintf (temp,"Null-Hypothesis: your test model");
        FPRINTF(world->outfile,"%-*.*s | Log(likelihood) of test model\n",BOXSIZE2,BOXSIZE2,temp);
        sprintf (temp,"=same=");
        FPRINTF(world->outfile,"%-*.*s | Log(likelihood) of full model\n",BOXSIZE2,BOXSIZE2, temp);
        sprintf (temp,"full model (the model under which the");
        FPRINTF(world->outfile,"%-*.*s | Likelihood ratio test value\n",BOXSIZE2,BOXSIZE2, temp);
        sprintf (temp,"genealogies were sampled)");
        FPRINTF(world->outfile,"%-*.*s | Degrees of freedom of test\n",BOXSIZE2,BOXSIZE2,temp);
        sprintf (temp,"[Theta values are on the diagonal of the ");
        FPRINTF(world->outfile,"%-*.*s | Probability*\n",BOXSIZE2,BOXSIZE2,temp);
        sprintf (temp,"Migration matrix, migration rates are ");
        FPRINTF(world->outfile,"%-*.*s | Probability**\n",BOXSIZE2,BOXSIZE2,temp);
        sprintf (temp,"specified as %s]", world->options->usem ? "M" : "Theta * M");
        FPRINTF(world->outfile,"%-*.*s | Akaike's Information Criterion***\n",BOXSIZE2,BOXSIZE2,temp);
        sprintf (temp," ");
        FPRINTF(world->outfile,"%-*.*s | Number of parameters used\n",BOXSIZE2,BOXSIZE2,temp);
        print_line(world->outfile,'-',79,CONT);
        FPRINTF(world->outfile,"  *) Probability under the assumption that parameters have range -Inf to Inf\n");
        FPRINTF(world->outfile," **) Probability under the assumption that parameters have range 0 to Inf\n");
        FPRINTF(world->outfile,"***) AIC: the smaller the value the better the model\n");
        aicfullparamnum = find_paramnum(world,NULL);
        aicfull = -2. * like1 + 2. * aicfullparamnum;
        FPRINTF(world->outfile,"          [the full model has AIC=%f, num(param)=%li]\n\n",aicfull,aicfullparamnum);

	pdf_print_lrt_header(aicfull,aicfullparamnum, &world->page_width, &world->page_height);
      }
    print_lrt_box(world,param0, param1, zeros, elem, this_string1, like0, like1, df);
    myfree(lparam0);
    myfree(lparam1);
    myfree(saveparam0);
    myfree(which);    
    myfree(hess[0]);
    myfree(hess);
}

void print_lrt_box(world_fmt *world, MYREAL *param0, MYREAL *param1, long zeros, long elem,
                   char * this_string1, MYREAL like0, MYREAL like1, long df)
{
    // Alternative plot
    //
    //---------------------------------------------------------------------
    // H0: (p1,p2,p3,...,              |   LRT  = -2(val1 - val2) = val
    //      pi,....,pn) == (v1,v2,     |   df   = x
    //      v3,.... vi,...vn)          |   Prob = x.xx
    //                                 |   Probc= x.xx
    //                                 |   AIC  = x.xxx
    //---------------------------------------------------------------------
    //
    // LRT
    MYREAL testval = -2. * (like0 - like1);
    // standard probability assuming indpendence and range of -inf .. +inf 
    MYREAL chiprob = probchi (df, testval);
    // probability assuming independence and range of 0 .. inf
    MYREAL chiprob2 = probchiboundary (testval, zeros, df);
    // AIC value penalizing for number of paramters
    long aicparamnum = find_paramnum(world,this_string1);
    MYREAL aic = -2. * like0 + 2. * aicparamnum;

    char **box;
    long size = HUNDRED; // we need at least 7 lines to print the right side of the table
    /// \todo  reallocation of likelihood ratio boxes needs reevaluation
    long newsize;
    long i;
    
    if(world->options->progress)
        print_line(stdout,'-',79,CONT);
    print_line(world->outfile,'-',79,CONT);
    box = (char **) mycalloc (size, sizeof (char *));
    for(i=0; i<size; i++)
        box[i] = (char *) mycalloc(LRATIO_STRINGS,sizeof(char));
    newsize = parse_h0(&box,&size, this_string1,param0,param1,elem,world);
    print_box(world,box,newsize,like0,like1,testval, chiprob, chiprob2, aic, df, aicparamnum);
    pdf_print_LRT_box(world,box,newsize,like0,like1,testval, chiprob, chiprob2, aic, df, aicparamnum, &world->page_height, & world->page_width);
    if(world->options->progress)
        print_line(stdout,'-',79,CONT);
    print_line(world->outfile,'-',79,CONT);
    for(i=0; i<size; i++)
        myfree(box[i]);
    myfree(box);
}

long parse_h0(char ***box, long *size, char *thisString, MYREAL *param0, MYREAL *param1 , long elem, world_fmt *world)
{
    long tempsize=0;
    parse_h0part(box, size, &tempsize, param0, elem, "H0:", world);
    tempsize++;
    parse_h0part(box, size, &tempsize, param1, elem, " = ", world);
    tempsize++;
    parse_h0string(box,size,&tempsize,thisString);
    return tempsize;
}

void parse_h0part(char ***box, long *size, long *tempsize,
                MYREAL *param, long elem, char head[], world_fmt *world)
{
    int length;
    long i;
    long zi;
    long pop;
    MYREAL tparam;
    char tmp[LRATIO_STRINGS];
    long counter = 0;
   //xcode  length = (int) MAX (0, 5 - sprintf (tmp, "%i", ((int) param[0])));
    counter = sprintf ((*box)[*tempsize],"%3.3s", head);
    for (i = 0; i < elem; i++)
      {
        zi = mml2m (i, world->numpop);
        if (zi >= world->numpop)
          {
            pop = (zi - world->numpop) / (world->numpop - 1);
            if(world->options->usem)
                tparam = param[zi];
            else
                tparam = param[zi] * param[pop];
          }
        else
            tparam = param[zi];
        length = MAX (0, 5 - sprintf (tmp, "%i", ((int) tparam)));
        counter += sprintf ((*box)[*tempsize] + counter, " %.*f", length, tparam);
        if (counter + length > BOXSIZE)
        {
            (*tempsize)++;
            if(*tempsize>= *size-1)
              {
                *box = (char **) myrealloc(*box,sizeof(char*) * (*tempsize+1));
                (*box)[*tempsize] = (char *) mycalloc(LRATIO_STRINGS,sizeof(char));
                *size = *tempsize+1; 
              }
            counter = sprintf((*box)[*tempsize],"%3.3s","   ");
        }
    }
}

void parse_h0string(char ***box, long *size, long *tempsize,
                  char *thisString)
{
    long i;
    long counter = 0;
    long strsize = (long) strlen(thisString);
    counter = sprintf ((*box)[*tempsize],"[");
    for (i = 0; i < strsize; i++)
      {
        counter += sprintf ((*box)[*tempsize] + counter, "%c",thisString[i]);
        if (counter >= BOXSIZE)
          {
            (*tempsize)++;
            if(*tempsize>= *size-1)
              {
                *box = (char **) myrealloc(*box,sizeof(char*) * (*tempsize+1));
                (*box)[*tempsize] = (char *) mycalloc(LRATIO_STRINGS,sizeof(char));
                *size = *tempsize+1;
              }
            counter = sprintf((*box)[*tempsize]," ");
          }
      }
    sprintf ((*box)[*tempsize] + counter, "]");
}

void print_box(world_fmt* world,char **box,long size,
               MYREAL like0, MYREAL like1,MYREAL lrt, MYREAL chiprob, MYREAL chiprob2, MYREAL aic, long df, long aicparamnum)
{
    long i;
    FPRINTF(world->outfile,"%-*.*s | LnL(test) = %f\n",BOXSIZE,BOXSIZE, box[0],like0);
    FPRINTF(world->outfile,"%-*.*s | LnL(full) = %f\n",BOXSIZE,BOXSIZE, box[1],like1);
    FPRINTF(world->outfile,"%-*.*s | LRT       = %f\n",BOXSIZE,BOXSIZE, box[2],lrt);
    FPRINTF(world->outfile,"%-*.*s | df        = %li\n",BOXSIZE,BOXSIZE,box[3],(long) df);
    FPRINTF(world->outfile,"%-*.*s | Prob      = %f\n",BOXSIZE,BOXSIZE,box[4],chiprob);
    FPRINTF(world->outfile,"%-*.*s | Probc     = %f\n",BOXSIZE,BOXSIZE,box[5],chiprob2);
    FPRINTF(world->outfile,"%-*.*s | AIC       = %f\n",BOXSIZE,BOXSIZE,box[6],aic);
    FPRINTF(world->outfile,"%-*.*s | num(param)= %li\n",BOXSIZE,BOXSIZE,box[7],aicparamnum);
    for(i=8; i< size; i++)
        FPRINTF(world->outfile,"%-*.*s |\n",BOXSIZE,BOXSIZE,box[i]);
    if(world->options->progress)
    {
        FPRINTF(stdout,"%-*.*s | LnL(test) = %f\n",BOXSIZE,BOXSIZE, box[0],like0);
        FPRINTF(stdout,"%-*.*s | LnL(full) = %f\n",BOXSIZE,BOXSIZE, box[1],like1);
        FPRINTF(stdout,"%-*.*s | LRT       = %f\n",BOXSIZE,BOXSIZE, box[2],lrt);
        FPRINTF(stdout,"%-*.*s | df        = %li\n",BOXSIZE,BOXSIZE,box[3],(long) df);
        FPRINTF(stdout,"%-*.*s | Prob      = %f\n",BOXSIZE,BOXSIZE,box[4],chiprob);
        FPRINTF(stdout,"%-*.*s | Probc     = %f\n",BOXSIZE,BOXSIZE,box[5],chiprob2);
        FPRINTF(stdout,"%-*.*s | AIC       = %f\n",BOXSIZE,BOXSIZE,box[6],aic);
        FPRINTF(stdout,"%-*.*s | num(param)= %li\n",BOXSIZE,BOXSIZE,box[7],aicparamnum);
        for(i=8; i< size; i++)
            FPRINTF(stdout,"%-*.*s |\n",BOXSIZE,BOXSIZE,box[i]);
    }
}


long
set_test_param (MYREAL *param, lr_data_fmt *lrtdata, world_fmt * world, long lrline,
                long locus, long *maxwhich, long *maxnum, long *zeros)
{
    long elem = world->options->gamma ? world->numpop2 + 1 : world->numpop2;
    long elements;
    long repstop = !world->options->replicate ? 0 :
        (world->options->replicatenum == 0 ?
         world->options->lchains : world->options->replicatenum);
    char *paramtype;
    char *ss;
    char **custm;
    MYREAL *meanparam;
    //long zzz;
    long zi;
    long el;
    long df = 0;
    long count=0;
    MYREAL mean=0.;
    long offset;
    long limit;
    long zz;
    long numpop = world->numpop;
    long pop1, pop2;
    
    paramtype = (char *) mycalloc (1, sizeof (char) * elem);
    ss = (char *) mycalloc (LONGLINESIZE, sizeof (char));
    charvec2d(&custm,world->numpop2 + 1, LONGLINESIZE);
    strcpy (ss, lrtdata->value1);

    elements = simplify_lrtvalues(ss,custm);

    *zeros = 0;
    
    if(elem != elements)
        warning("Not enough elements in the l-ratio specification\n");

   //xcode  zzz=0;
    if (world->loci - world->skipped > 1)
        meanparam = world->atl[0][world->loci].param;
    else
        meanparam = world->atl[repstop][0].param;

    for(el=0; el < elements; el++)
    {
        zi = mml2m (el, numpop);
        m2mm (zi, numpop, &pop1, &pop2);

        switch (custm[el][0])
        {
            case 'x':
                paramtype[zi] = '-';
                lrtdata->connect[zi]='x';
                param[zi] = meanparam[zi];
                maxwhich[(*maxnum)++] = zi;
                df++;
                break;
            case '*':
                paramtype[zi] = '-';
                lrtdata->connect[zi]='*';
                param[zi] = meanparam[zi];
                break;
            case 't':
            case 'm':
                paramtype[zi] = '-';
                lrtdata->connect[zi]='m';
                if(custm[el][0] != world->options->custm2[zi])
                    df++;
                    mean = 0.0;
                    count = 0;
                    offset = (zi >= world->numpop) ? world->numpop : 0;
                    limit = (zi >= world->numpop) ? world->numpop2 : world->numpop;
                    for (zz = offset; zz < limit; zz++)
                        {
                        if(custm[m2mml(zz,numpop)][0] == 'm')
                            {
                            mean += meanparam[zz];
                            count++;
                            }
                        }
                    mean /= count; //limit - offset;
                    param[zi] = mean;
                break;
            case 'M':
                if(custm[el][0] != world->options->custm2[zi])
                    df++;
                mean = 0.0;
                count = 0;
                offset = (zi >= world->numpop) ? world->numpop : 0;
                limit =
                    (zi >= world->numpop) ? world->numpop2 : world->numpop;
                if (offset < numpop)
                    {
                        paramtype[zi] = '-';
                        lrtdata->connect[zi]='m';
                        for (zz = offset; zz < limit; zz++)
                          {
                            if(custm[m2mml(zz,numpop)][0] == 'm')
                            {
                                mean += meanparam[zz];
                                count++;
                            }
                          }
                        mean /= count; //limit - offset;
                        param[zi] = mean;
                    }
                else
                    {
                        paramtype[zi] = '+';
                        lrtdata->connect[zi]='M';
                        for (zz = offset; zz < limit; zz++)
                        {
                            if(custm[m2mml(zz,numpop)][0] == 'm')
                              {
                                m2mm (zz, numpop, &pop1, &pop2);
                                mean += meanparam[zz] * meanparam[pop2];
                                count++;
                              }
                        }
                        mean /= count;
                        param[zi] = mean;
                    }
                break;
            case 's':
                if(custm[el][0] != world->options->custm2[zi])
                    df++;
                if (zi < world->numpop)
                {
                    paramtype[zi] = '-';
                    lrtdata->connect[zi]='*';
                    param[zi] = meanparam[zi];
                }
                    else
                    {
                        paramtype[zi] = '-';
                        lrtdata->connect[zi]='s';
                        param[zi] =
                            (meanparam[zi] +
                            meanparam[mm2m (pop2, pop1, world->numpop)]) / 2.;
                    }
                    break;
            case 'S':
                if(custm[el][0] != world->options->custm2[zi])
                    df++;
                if (zi < world->numpop)
                {
                    paramtype[zi] = '-';
                    lrtdata->connect[zi]='*';
                    param[zi] = meanparam[zi];
                }
                    else
                    {
                        paramtype[zi] = '+';
                        lrtdata->connect[zi]='S';
                        param[zi] =
                            (meanparam[zi] * meanparam[pop2] +
                            meanparam[mm2m (pop2, pop1, world->numpop)] *
                            meanparam[pop1]) / 2.;
                    }
                    break;
            default:
                paramtype[zi] = '+';
                if(custm[el][0] != world->options->custm2[zi])
                    df++;
                if (custm[el][0] == '0' && (long) strlen(custm[el])==1)
                {
                    if(custm[el][0] != world->options->custm2[zi])
                        (*zeros)++;
                    lrtdata->connect[zi]='0';
                }
                else
                    lrtdata->connect[zi]='c';
                param[zi] = MAX (atof (custm[el]), SMALLEST_THETA);
                break;
        }
        
    }
    for (zi = world->numpop; zi < world->numpop2; zi++)
      {
	if (paramtype[zi] == '+')
	  {
	    zz = (zi - world->numpop) / (world->numpop - 1);
	    if(!world->options->usem)
	      param[zi] /= param[zz];
	  }
      }
    myfree(paramtype);
    myfree(ss);
    free_charvec2d(custm);
    return df;
}


long simplify_lrtvalues(char *in, char **out)
{
  char *ptr;
  long elements=0;
  while(in != NULL)
    {
      
      ptr = strsep(&in,"\n\r\t; {,}");
      if(strlen(ptr)>0)
	{
	  if(ptr[0]=='0' && ((atof(ptr) - SMALLEPSILON) <= 0.0))
	    sprintf(out[elements++], "0");
	  else
	    sprintf(out[elements++],"%s", ptr);
	}
    }
  return elements;
}

// in is the string from the l-ratio option, custom is a database array that
// gets allocated while its filled
// maximal event numbers is LONGLINESIZE [search for use of simplify_lrtvalues()]
long simplify_lrtvalues_old(char * in, char*** out)
{
    long elements=0;
    long count = 0;
    boolean white=FALSE;
    boolean element_ready = TRUE;
    //    (*out)[0] = (char *) mycalloc(LONGLINESIZE,sizeof(char));
    while(*in!='\0')
      {
        switch(*in)
          {
            case '\r':
            case '\n':
            case '\t':
            case ' ':
            case ',':
            case ';':
	      white=TRUE;
	      break; // jump over "whitespace", 
	      // if elements are from the list {*, m, M, s, S} read one character and then restart a new element
	  case 't':
	  case 'm':
	  case 'M':
	  case 's':
	  case 'S':
	  case '*':
	  case 'x':
	    white = FALSE;
	    // fill element and reset counter within element
	    count=0;
	    (*out)[elements][0] = *in;
	    // make new element ready
	    elements++;
	    //    (*out)[elements] = (char *) mycalloc(LONGLINESIZE,sizeof(char));
	    element_ready = TRUE;
	    break;
	  default:
	    if(!element_ready)
	      {
		(*out)[elements][count++] = *in;
	      }
	    else
	      {
		if(white)
		  {
		    count=0;
		    elements++;
		    //    (*out)[elements] = (char *) mycalloc(LONGLINESIZE,sizeof(char));
		    element_ready = FALSE;
		    (*out)[elements][count++] = *in;
		  }
		else
		  {
		    (*out)[elements][count++] = *in;
		  }
	      }
	    white = FALSE;
	    break;
          }
        in++;
    }
    //*out = (char**) myrealloc(*out, elements * sizeof(char*));
    return elements;
}

#if 0
long
set_test_param (MYREAL *param, lr_data_fmt *lrtdata, world_fmt * world, long lrline,
                long locus, long *maxwhich, long *maxnum, long *zeros)
{
    long i = 0, z = 0, zi = 0, zz = 0, zzz = 0, df = 0;
    long count;
    long offset = 0, limit = 0, pop, pop1, pop2;
    long numpop = world->numpop;
    char *tmp, *tmp2, *ss, *custm, *paramtype;
    MYREAL *meanparam, mean;
    long elem = world->options->gamma ? world->numpop2 + 1 : world->numpop2;
    long repstop = !world->options->replicate ? 0 :
                   (world->options->replicatenum == 0 ?
                    world->options->lchains : world->options->replicatenum);
    *zeros = 0;
    ss = (char *) mycalloc (1, sizeof (char) * LOGLINESIZE);
    custm = (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
    tmp = (char *) mycalloc (1, sizeof (char) * LONGLINESIZE);
    paramtype = (char *) mycalloc (1, sizeof (char) * elem);

    *maxnum = 0;
    strcpy (ss, lrtdata->value1);
    strcpy (custm, lrtdata->value1);
    tmp2 = custm;
    zzz = 0;
    while(*tmp2!='\0')
      {
        if(strchr("tmMsS*xc",*tmp2))
            custm[zzz++] = *tmp2;
        if(strchr("0",*tmp2))
          {
            if(!strchr(".0123456789",*(tmp2+1)))
                custm[zzz++] = *tmp2;
          }
        tmp2++;
      }
    custm[zzz]="\0";
    zzz=0;
    if (world->loci - world->skipped > 1)
        meanparam = world->atl[0][world->loci].param;
    else
        meanparam = world->atl[repstop][0].param;
    while (ss[zzz] != '\0')
    {
        tmp[i] = ss[zzz++];
        if (mywhitespace (tmp[i]))
        {
            if (tmp[i] != ' ')
                i++;
        }
        else
        {
            tmp[i] = '\0';
            i = 0;
            zi = mml2m (z, numpop);
            m2mm (zi, numpop, &pop1, &pop2);
            switch (tmp[0])
            {
            case 'x':
                paramtype[zi] = '-';
                lrtdata->connect[zi]='x';
                param[zi] = meanparam[zi];
                maxwhich[(*maxnum)++] = zi;
                z++;
                df++;
                break;
            case '*':
                paramtype[zi] = '-';
                lrtdata->connect[zi]='*';
                param[zi] = meanparam[zi];
                z++;
                break;
            case 't':
            case 'm':
                paramtype[zi] = '-';
                lrtdata->connect[zi]='m';
                zz = atol (tmp) - 1;
                tmp2 = tmp;
                if(tmp[0] != world->options->custm2[zi])
                    df++;
                if (zz < 0)
                {
                    mean = 0.0;
                    count = 0;
                    offset = (zi >= world->numpop) ? world->numpop : 0;
                    limit = (zi >= world->numpop) ? world->numpop2 : world->numpop;
                    for (zz = offset; zz < limit; zz++)
                    {
                        if(custm[m2mml(zz,world->numpop)] == 'm')
                          {
                            mean += meanparam[zz];
                            count++;
                          }
                    }
                    mean /= count; //limit - offset;
                    param[zi] = mean;
                }
                else
                {
                    param[zi] = meanparam[zz];
                }
                z++;
                break;
            case 'M':
                zz = atol (tmp) - 1;
                df++;
                if (zz < 0)
                {
                    mean = 0.0;
                    offset = (zi >= world->numpop) ? world->numpop : 0;
                    limit =
                        (zi >= world->numpop) ? world->numpop2 : world->numpop;
                    if (offset < numpop)
                    {
                        paramtype[zi] = '-';
                        lrtdata->connect[zi]='m';
                        for (zz = offset; zz < limit; zz++)
                            mean += meanparam[zz];
                        mean /= limit - offset;
                        param[zi] = mean;
                    }
                    else
                    {
                        paramtype[zi] = '+';
                        lrtdata->connect[zi]='M';
                        for (zz = offset; zz < limit; zz++)
                        {
                            m2mm (zz, numpop, &pop1, &pop2);
                            mean += meanparam[zz] * meanparam[pop2];
                        }
                        mean /= limit - offset;
                        param[zi] = mean;
                    }
                }
                else
                {
                    paramtype[zi] = '-';
                    lrtdata->connect[zi]='M';
                    param[zi] = meanparam[zz];
                }
                z++;
                break;
            case 's':
                df++;
                if (zi < world->numpop)
                {
                    paramtype[zi] = '-';
                    lrtdata->connect[zi]='*';
                    param[zi] = meanparam[zi];
                    z++;
                }
                else
                {
                    paramtype[zi] = '-';
                    lrtdata->connect[zi]='s';
                    param[zi] =
                        (meanparam[zi] +
                         meanparam[mm2m (pop2, pop1, world->numpop)]) / 2.;
                    z++;
                }
                break;
            case 'S':
                df++;
                if (zi < world->numpop)
                {
                    paramtype[zi] = '-';
                    lrtdata->connect[zi]='*';
                    param[zi] = meanparam[zi];
                    z++;
                }
                else
                {
                    paramtype[zi] = '+';
                    lrtdata->connect[zi]='S';
                    param[zi] =
                        (meanparam[zi] * meanparam[pop2] +
                         meanparam[mm2m (pop2, pop1, world->numpop)] *
                         meanparam[pop1]) / 2.;
                    z++;
                }
                break;
            default:
                paramtype[zi] = '+';
                df++;
                if (tmp[0] == '0')
                  {
                    (*zeros)++;
                    lrtdata->connect[zi]='0';
                  }
                else
                    lrtdata->connect[zi]='c';
                param[zi] = MAX (atof (tmp), SMALLEST_THETA);
                z++;
                break;
            }
        }
    }
    for (zi = world->numpop; zi < world->numpop2; zi++)
    {
        if (paramtype[zi] == '+')
        {
            pop = (zi - world->numpop) / (world->numpop - 1);
            if(!world->options->usem)
                param[zi] /= param[pop];
        }
    }
    myfree(paramtype);
    myfree(ss);
    myfree(tmp);
    myfree(custm);
    return df;
}

#endif
