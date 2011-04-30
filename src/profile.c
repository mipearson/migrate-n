/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 P R O F I L E    L I K E L I H O O D    R O U T I N E S 
 
 Peter Beerli 1997, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
$Id: profile.c 1800 2011-01-29 13:40:00Z beerli $
 
-------------------------------------------------------*/
/*! \file profile.c
Calculation of profile likelihoods and printing of profile likelihood tables
*/

#include "migration.h"

#include "world.h"
#include "laguerre.h"
#include "tools.h"
#include "broyden.h"
#include "combroyden.h"
#include "spline.h"
#include "joint-chains.h"
#include "sighandler.h"

#include "migrate_mpi.h"

#include "profile.h"
#ifdef PRETTY
#include "pretty.h"
#endif

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

extern int myID;

long find_profilelike (MYREAL testlike, MYREAL prob, long whichprob, MYREAL minparam,
                       MYREAL maxparam, MYREAL **testparam, MYREAL *llike,
                       long which, world_fmt * world, boolean QD, boolean QD2,
                       MYREAL *mlparam, MYREAL *normd, nr_fmt * nr);
void calc_profile_likelihood (char method, long which, MYREAL *likes,
                              MYREAL **param, world_fmt * world, long nn,
                              nr_fmt * nr);
void print_profile_likelihood (long which, world_fmt * world, long *gmaxptr);
void prepare_header (char *var, long which, world_fmt * world);
void print_profile_title (world_fmt * world);
long print_profile_table (char method, long which, char *var, MYREAL likes[],
                          MYREAL **param, world_fmt * world);
long print_profile_percentile (world_fmt * world);
void sprint_percentile_header (char *outfile, boolean first);
void allocate_profile_percentiles (world_fmt * world);
void destroy_profile_percentiles (world_fmt * world);
void print_menu_profile (long which, long nn);
MYREAL interpolate (MYREAL xl, MYREAL xc, MYREAL xh, MYREAL low,
                    MYREAL center, MYREAL high, MYREAL testlike);

long warp (long ii);

void prognose_profile_end (time_t starttime, long numpop2, long nn);


void profile_max_percentiles (long which, MYREAL *likes, MYREAL **param,
                              world_fmt * world, long nn, nr_fmt * nr);

void profile_max_spline (char method, long which, MYREAL *likes,
                         MYREAL **param, world_fmt * world, nr_fmt * nr,
                         long nn);

int calc_spline (MYREAL *param, MYREAL *like, long nn, long *constr,
                 MYREAL *diff, MYREAL *diff2, long *diagn, MYREAL *work,
                 long nwork);
void setup_spline (spline_fmt * spline);

void destroy_spline (spline_fmt * spline);

void prepare_spline_nodes (long n, long which, MYREAL **rawparam,
                           MYREAL *rawlikes, long *indeks, MYREAL *param,
                           MYREAL *likes);

void prepare_spline_first (MYREAL *param, MYREAL *likes);

void prepare_spline_last (long n, MYREAL *param, MYREAL *likes);

void calc_spline_diff (MYREAL *diff, MYREAL *param, MYREAL *likes, long n);

MYREAL recalc_profilelike (MYREAL testparam, MYREAL *param, MYREAL *like,
                           long nn, MYREAL *diff, MYREAL *diff2, MYREAL *work,
                           long nwork);


long set_df (long which, char *custm2, long numpop, long numpop2);


long set_profile_param (MYREAL *param, long which, long *allwhich,
                        MYREAL xval, MYREAL *xvals, world_fmt * world);
void addinvar (MYREAL *param, long *fixed, long *allwhich, MYREAL *xvals,
               long start, long len);

long sprint_nice_param (MYREAL parameter, MYREAL bottom, MYREAL top, boolean failed,
                        char *file);

#ifdef PRIVATE
void master_gridder(world_fmt *world, long *gmaxptr);
void calc_grid(FILE *alldim, world_fmt * world, MYREAL *testparam, MYREAL *ltestparam, MYREAL *mlparam, nr_fmt *nr, long which);
#endif




boolean
print_profile_likelihood_driver (long which, world_fmt * world, long *gmaxptr)
{
    static boolean onetheta = FALSE;
    static boolean onemig = FALSE;
    long i, j;
    //world->options->verbose = FALSE;
    world->locus = 0;
    if (which == world->numpop2 && world->options->gamma) //gamma deviated mutation rate
        print_profile_likelihood (which, world, gmaxptr);
    else
    {
        if(which>world->numpop2)
            return 0;
        switch (world->options->custm2[which])
        {
        case 'C':
        case 'c':
        case '0':
            return 0;
        case 'M':
        case 'm':
            if (which < world->numpop && onetheta)
                return 0;
            if (which < world->numpop)
            {
                print_profile_likelihood (which, world, gmaxptr);
                onetheta = TRUE;
            }
            if (which >= world->numpop && onemig)
                return 0;
            if (which >= world->numpop)
            {
                print_profile_likelihood (which, world, gmaxptr);
                onemig = TRUE;
            }
            return 1;
        case 'S': // Nm is symmetrical
            j = (which - world->numpop) % (world->numpop - 1);
            i = (which - world->numpop) / (world->numpop - 1);
            if(world->options->usem)
            {
                print_profile_likelihood (which, world, gmaxptr);
                return 1;
            }
            else
            {
                if (i <= j)
                {
                    print_profile_likelihood (which, world, gmaxptr);
                    return 1;
                }
            }
            return 0;
        case 's': // M is symmetrical
            j = (which - world->numpop) % (world->numpop - 1);
            i = (which - world->numpop) / (world->numpop - 1);
            if(!world->options->usem)
            {
                print_profile_likelihood (which, world, gmaxptr);
                return 1;
            }
            else
            {
                if (i <= j)
                {
                    print_profile_likelihood (which, world, gmaxptr);
                    return 1;
                }
            }
            return 0;
        default:
            print_profile_likelihood (which, world, gmaxptr);
        }
    }
    return 1;
}

void
print_profile_likelihood (long which, world_fmt * world, long *gmaxptr)
{
    long i, ii;
    long bufsize = 0;
    nr_fmt *nr;
    char *var;
    char method = world->options->profilemethod;
    MYREAL **param;
    MYREAL *likes;
    long elem = world->options->gamma ? world->numpop * world->numpop + 1 :
                world->numpop * world->numpop;
#ifdef LONGSUM

    elem += world->numpop * 3;
#endif /*LONGSUM*/

    nr = (nr_fmt *) mycalloc (1, sizeof (nr_fmt));
    create_nr (nr, world, *gmaxptr, which, world->loci,
               world->repkind, world->rep);
    
    var = (char *) mycalloc(LINESIZE, sizeof(char));
    doublevec1d(&likes, GRIDSIZE);
    //doublevec2d(&param,GRIDSIZE, elem);
    param = (MYREAL **) mycalloc(GRIDSIZE,sizeof(MYREAL));
    for(i=0;i<GRIDSIZE;i++)
        param[i] = (MYREAL*) mycalloc(elem+1,sizeof(MYREAL));
    
    if (world->options->progress)
    {
        print_menu_profile (which,
                            world->numpop2 + (long) world->options->gamma);
    }
    //#ifndef INTEGRATEDLIKE
    SETUPPARAM0 (world, nr, world->repkind, nr->repstart, nr->repstop,
                 world->loci, PROFILE, world->loci>1 ?  TRUE : FALSE);
    //#endif
    calc_profile_likelihood (method, which, likes, param, world, GRIDSIZE, nr);
    prepare_header (var, which, world);
    if (world->options->printprofile)
    {
      bufsize = print_profile_table (method, which, var, likes, param, world);
      world->allocbufsize = bufsize;
#ifdef PRETTY
      prepare_profile_table(method, which, likes, param, world, bufsize);
#endif
    }
    if (world->options->printprofsummary)
    {
        strcpy (world->quantiles[which].name, var);
        for (ii = 0; ii < GRIDSIZE; ii++)
        {
            i = warp (ii);
            world->quantiles[which].param[i] = param[i][which];
        }
    }
    destroy_nr (nr, world);
    myfree(var);
    myfree(likes);
    //myfree(param[0]);
    myfree(param);
}


/// 
/// print the profile likelihood table
/// depending on the profile option
/// prints a mark to the precentiles that did not succeed to converge to the 
/// proper percentile
/// \return bufsize size of buffer
long
print_profile_table (char method, long which, char *var, MYREAL likes[],
                     MYREAL **param, world_fmt * world)
{
    const MYREAL probabilities[] = SHOWGRID;
    long failed = 0;
    char star;
    char methodstring[LINESIZE];

    long i, ii, j, likei = 0;
    long bufsize = 1;
    long elem =
        world->options->gamma ? world->numpop * world->numpop +
        1 : world->numpop * world->numpop;

    char fp[LONGLINESIZE]; //large local string buffer
    char **buffer = &world->buffer; // buffer for printing whole table
    char *temp, *var2, temp2[LINESIZE], likestr[LINESIZE], *paramstr;
    MYREAL likemax = -MYREAL_MAX;
#ifdef LONGSUM

    elem += world->numpop * 3;
#endif /*LONGSUM*/

    //memset (*buffer, 0, sizeof (char) * (long) strlen (*buffer));
    (*buffer)[0] = '\0';
    var2 = (char *) mycalloc (1, sizeof (char) * MAX (LINESIZE, (long) (elem / 10. * LINESIZE)));
    paramstr = (char *) mycalloc (1, sizeof (char) * MAX (LINESIZE,(long) (elem / 10. * LINESIZE)));
    temp = var2;
    for (j = 0; j < elem; j++)
    {
        prepare_header (temp2, j, world);
        if (which == j)
            star = '*';
        else
            star = ' ';
        if (((j + 1) % 5) == 0)
        {
            sprintf (temp, "\n                            %c%.7s%c     ", star,
                     temp2, star);
            temp += 28;
        }
        else
            sprintf (temp, "%c%.7s%c    ", star, temp2, star);
        temp += 11;
    }
    switch (method)
    {
    case 'p':
        strcpy (methodstring,
                "Parameters are evaluated at percentiles\nusing bisection method (slow, but exact).");
        break;
    case 's':
        strcpy (methodstring,
                "Parameters are evaluated at percentiles\nusing cubic splines of profiled parameter\n(faster, but not so exact).");
        break;
    case 'd':
        strcpy (methodstring,
                "Parameters are evaluated at pre-defined values\n");
        break;
    case 'q':
        strcpy (methodstring,
                "Parameters are evaluated assuming complete independence\n");
        break;
    case 'f':
        strcpy (methodstring,
                "Parameters are evaluated assuming complete independence\n and then once maximized at the found profiled parameter value");
        break;
    }
    sprintf (fp, "\n\nProfile likelihood for parameter %s\n", var);
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    sprintf (fp, "%s\n", methodstring);
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    sprint_line (fp, '-', 79, CONT);
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    if (method == 'd')
        sprintf (fp, "      Ln(L)     %7.7s     %s\n", var, var2);
    else
        sprintf (fp, "Per.  Ln(L)     %7.7s     %s\n", var, var2);
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    sprint_line (fp, '-', 79, START);
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    for (i = 0; i < GRIDSIZE; i++)
    {
        if (likemax < MIN (likes[i], MYREAL_MAX))
        {
            likemax = likes[i];
            likei = i;
        }
    }
    for (ii = 0; ii < GRIDSIZE; ii++)
    {
        i = warp (ii);
        if (likes[i] >= MYREAL_MAX - EPSILON)
            sprintf (likestr, "     -    ");
        else
            sprintf (likestr, "% 6.3f%c", likes[i], likei == i ? '*' : ' ');
        temp = paramstr;
        for (j = 0; j < elem; j++)
        {
            if (((j + 1) % 5) == 0)
            {
                sprintf (temp, "\n                               ");
                temp += 26;
            }
            if (param[i][j] <= 0)
                sprintf (temp, "      -     ");
            else
            {
                if (param[i][j] < 0.000001)
                    sprintf (temp, " % 6.3e ", param[i][j]);
                else
                    sprintf (temp, "% 10.6f ", param[i][j]);
            }

            temp += 11;
        }
        if (param[i][which] <= 0)
        {
            if (method == 'd')
                sprintf (fp, "        -           -       %s\n", paramstr);
            else
                sprintf (fp, "%3.3f     -           -       %s\n",
                         probabilities[i], paramstr);
            bufsize += (long) strlen(fp);
            *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
            strcat (*buffer, fp);
        }
        else
        {
            if (method == 'd')
                sprintf (fp, "    %8.8s % 10.6g %s\n", likestr,
                         param[i][which], paramstr);
            else
            {
                if (probabilities[i] == 0.5)
                    sprintf (fp, "MLE   %8.8s % 10.6g %s\n",
                             likestr, param[i][which], paramstr);
                else
                {
                    if(world->percentile_failed[which][i])
		      {
                        sprintf (fp, "**** %8.8s % 10.6g %s\n",
                                 likestr,
                                 param[i][which], paramstr);
			failed += 1;
		      }
                    else
                        sprintf (fp, "%4.3f %8.8s % 10.6g %s\n",
                                 probabilities[i], likestr,
                                 param[i][which], paramstr);
                }
            }
            bufsize += (long) strlen(fp);
            *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
            strcat (*buffer, fp);
        }
    }
    sprint_line (fp, '-', 79, STOP);
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    if(failed>0)
      {
	sprintf (fp, "If the percentile column is marked with **** then the convergence to the percentile value failed\nThe values are still _correct_ but not at the percentile of the profile parameter\n");
	bufsize += (long) strlen(fp);
	*buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
	strcat (*buffer, fp);
      }
    myfree(paramstr);
    myfree(var2);
    return bufsize;
}


void
print_profile_title (world_fmt * world)
{
    print_line (world->outfile, '=', 79, CONT);
    FPRINTF (world->outfile, "Profile likelihood tables   [Summary is at the end of the file]\n");
    print_line (world->outfile, '=', 79, CONT);
}


///
/// set the profile method for ascii output
void ascii_method_set(FILE *outfile, char method)
{
    switch (method)
    {
    case 'p':
      fprintf(outfile, "Parameters are evaluated at percentiles using bisection method (slow, but exact).\n");
      break;
    case 's':
      fprintf(outfile,"Parameters are evaluated at percentiles\n");
      fprintf(outfile,"using cubic splines of profiled parameter (faster, but not so exact).\n");
      break;
    case 'd':
    fprintf(outfile,"Parameters are evaluated at pre-defined values\n");
    break;
    case 'q':
     fprintf(outfile, "Parameters are evaluated assuming complete independence\n");
      break;
    case 'f':
      fprintf(outfile, "Parameters are evaluated assuming complete independence\n");
      fprintf(outfile,"and then once maximized at the found profiled parameter value\n");
      break;
    }
    fprintf(outfile,"\n");
}

/*++++++++++++++++++++++++++ replacement for ascii table, not finished *************/
#if 0
///
/// generic ascii table generator
/// \param float left_margin left edge of table
/// \param float * page_height page height 
/// \param int cols number of columns in table if there are too many columns then new line
/// \param int rows number of rows in table, if a new page is needed then header is repeated
/// \param char *** elements all table elements, currently no formatting of these
/// \param char ** header   header row
/// \param char *position position of all elements: L=left, R=right, C=center
/// \param int col_overflow when there are too many columns this is the column to restart
long ascii_table(int cols, int rows, char *buffer, char *position, int col_overflow, float separator)
{
  boolean new_page=FALSE;

  int row;
  int col;
  long location=0;

  char pos;

  float *col_widths;
  float left_margin = 0.;
  float page_width = 130.;
  float right_margin = 130.;
  float *col_starts;
  float oldcol = HUGE;

  char **header;
  char ***elements;

  col_widths = mycalloc(cols,sizeof(float));
  col_starts = mycalloc(cols,sizeof(float));

  charvec2d(&header, cols, LINESIZE);
  elements = (char ***) mycalloc(rows, sizeof(char **));
  for(row=0; row < rows; row++)
    {
      charvec2d(&(elements[row]),cols, LINESIZE);
    }
  
  location = translate_buffer_table(cols, rows, buffer, header, elements);
  ascii_find_col_width(cols, rows, elements, header, col_widths);
  define_col_start(cols, col_widths, col_overflow, position, left_margin, right_margin, separator, col_starts);
  
  ascii_print_table_header(page_height, position, cols, col_starts,header);
  sprint_line (fp, '-', 130, CONT);

  for(row=0; row< rows; row++)
    {
      for(col=0;col < cols; col++)
	{
	  if(col_starts[col] < oldcol)
	    {
	      fprintf(outfile,"\n");
	    }
	  if(position==NULL)
	    pos='L';
	  else
	    pos=position[col];
	  
	  fprintf(outfilecol_starts[col],*page_height, pos, "%s",elements[row][col]);
	  oldcol = col_starts[col];
	}
    }
  pdf_advance(page_height);
  myfree(col_widths);
  myfree(col_starts);
  free_charvec2d(header);
  for(row=0; row < rows; row++)
    free_charvec2d(elements[row]);
  myfree(elements);
  return location;
}


///
/// Printing PROFILE table for ASCII
/// \param method       profile method
/// \param world        master parameter structure
void ascii_print_profile_table(char method, world_fmt *world)
{
  char *thisparam;
  char *savedparam;
  char *ptr;
  char *buffer = world->buffer;
  long position = 0;
  FILE *outfile = world->outfile;

  boolean failed = world->percentile_some_failed;

  thisparam = (char *) mycalloc(LINESIZE,sizeof(char));
  savedparam = thisparam;
  while(buffer[position] != '\0')
    {
      strncpy(thisparam, buffer, LINESIZE-1);
      strsep(&thisparam,"&"); // remove first element
      strsep(&thisparam,"&"); // remove second  element
      ptr = strsep(&thisparam,"&"); // keep third 
      
      fprintf(outfile,"\n\n");
      fprintf(outfile, "Profile likelihood table for parameter %s\n", ptr);
      ascii_method_set(outfile, method);
      position += pdf_table(left_margin, page_height,world->numpop2 + 3, 9, buffer + position, 
			    NULL, 4, 6.);
      pdf_table_footnote(left_margin, page_height, failed);
    }
  myfree(savedparam);
}
#endif /*++++++++++++++++++++++++++ replacement for ascii table, not finished *************/


long
warp (long ii)
{
    long i;
    if (ii < 4)
        i = ii;
    else
    {
        if (ii == 4)
        {
            i = GRIDSIZE - 1;
        }
        else
        {
            i = GRIDSIZE - ii + 3;
        }
    }
    return i;
}


///
/// print profile tables for percentiles, returns the number
/// of failed convergences to a specific percentile.
long
print_profile_percentile (world_fmt * world)
{
    const MYREAL probabilities[] = SHOWGRID; //debug GRID2
    long i, j, jj;
    boolean first = TRUE;
    boolean last = FALSE;
    char fp[LINESIZE];
    char **buffer = &world->buffer;
    long bufsize = 1;
    long failed=0;


    (*buffer)[0]='\0';
    
    sprintf (fp, "\n\n");
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    
    sprint_line (fp, '=', 79, CONT);
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    sprintf (fp,
             "Summary of profile likelihood percentiles of all parameters\n");
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    sprint_line (fp, '=', 79, CONT);
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    sprintf (fp, "\n");
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    sprint_percentile_header (fp, first);
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    for (i = 0; i < world->numpop2 + (long) world->options->gamma; i++)
    {
        if (world->quantiles[i].name[0] == '\0')
            continue;  /* this variable was not estimated, so return */
        sprintf (fp, "%-10.10s  ", world->quantiles[i].name);
        bufsize += (long) strlen(fp);
        *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
        strcat (*buffer, fp);
        for (jj = 0; jj < GRIDSIZE; jj++)
        {
            j = warp (jj);
            if (probabilities[j] > 0.5)
                continue;
            if (world->quantiles[i].param[j] < 0)
                sprintf (fp, "      -       ");
            else
                failed += sprint_nice_param (world->quantiles[i].param[j], 0.000001,
                                   999.99999, world->percentile_failed[i][j], fp);
            bufsize += (long) strlen(fp);
            *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
            strcat (*buffer, fp);
        }
        sprintf (fp, "\n");
        bufsize += (long) strlen(fp);
        *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
        strcat (*buffer, fp);
    }
    sprintf (fp, "\n\n");
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    sprint_percentile_header (fp, last);
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    for (i = 0; i < world->numpop2 + (long) world->options->gamma; i++)
    {
        if (world->quantiles[i].name[0] == '\0')
            continue;  /* this variable was not estimated, so return */
        sprintf (fp, "%-10.10s  ", world->quantiles[i].name);
        bufsize += (long) strlen(fp);
        *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
        strcat (*buffer, fp);
        for (jj = 0; jj < GRIDSIZE; jj++)
        {
            j = warp (jj);
            if (probabilities[j] < 0.5)
                continue;
            if (world->quantiles[i].param[j] < 0)
                sprintf (fp, "      -       ");
            else
                failed += sprint_nice_param (world->quantiles[i].param[j], 0.000001,
                                   999.99999, world->percentile_failed[i][j], fp);
            bufsize += (long) strlen(fp);
            *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
            strcat (*buffer, fp);
        }
        sprintf (fp, "\n");
        bufsize += (long) strlen(fp);
        *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
        strcat (*buffer, fp);
    }

    sprintf (fp, "\n");
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    sprint_line (fp, '-', 80, CONT);
    bufsize += (long) strlen(fp);
    *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
    strcat (*buffer, fp);
    if(failed > 0)
    {
        sprintf (fp, "Values with a * are NOT at the percentiles!\nThe convergence to the percentile value FAILED.\n");
        bufsize += (long) strlen(fp);
        *buffer = (char *) myrealloc(*buffer, sizeof(char) * bufsize);
        strcat (*buffer, fp);
    }
    return bufsize;
}

long
sprint_nice_param (MYREAL parameter, MYREAL bottom, MYREAL top, boolean failed, char *file)
{
  long failedcount=0;
    if (parameter > bottom && parameter < top)
    {
        if(failed)
	  {
            failedcount = 1;
            sprintf (file, " %11.6f* ", parameter);
	  }
        else
            sprintf (file, " %11.6f  ", parameter);
    }
    else
    {
        if(failed)
	  {
            failedcount = 1;
	    sprintf (file, " %11.6g* ", parameter);
	  }
        else
            sprintf (file, " %11.6g  ", parameter);
    }
    return failedcount;
}

void
sprint_percentile_header (char *outfile, boolean first)
{
    const MYREAL probabilities[] = GRID2;
    long i;
    char fp[LINESIZE];
    if (first)
        sprintf (fp, "Parameter                          Lower percentiles\n");
    else
        sprintf (fp, "Parameter                          Upper percentiles\n");
    strcat (outfile, fp);
    sprintf (fp, "            ");
    strcat (outfile, fp);
    sprint_line (fp, '-', 68, CONT);
    strcat (outfile, fp);
    sprintf (fp, "            ");
    strcat (outfile, fp);
    for (i = 0; i < GRIDSIZE - 1; i++)
    {
        if (first)
        {
            if (probabilities[i] >= 0.5)
                break;
        }
        else
        {
            if (probabilities[i] < 0.5)
                continue;
        }
        if (probabilities[i] == 0.5)
            sprintf (fp, "    MLE       ");
        else
            sprintf (fp, "    %5.3f     ", probabilities[i]);
        strcat (outfile, fp);
    }
    if (probabilities[i] == 0.5)
        sprintf (fp, "    MLE\n");
    else
        sprintf (fp, "    %5.3f\n", probabilities[i]);
    strcat (outfile, fp);
    sprint_line (fp, '-', 80, CONT);
    strcat (outfile, fp);
}


void
prepare_header (char *var, long which, world_fmt * world)
{
    long i, j, from, to;
#ifdef LONGSUM

    long len = world->numpop2 + world->options->gamma;
    len += world->numpop * 3;
#endif /*LONGSUM*/

    if (which < world->numpop)
    {
        sprintf (var, "Theta_%li", from = which + 1);
    }
    else
    {
        if (world->numpop > 1)
        {
            i = (which - world->numpop) / (long) (world->numpop - 1);
            j = which - world->numpop - i * (world->numpop - 1);
            to = i + 1;
            from = (i <= j) ? j + 2 : j + 1;
            if (world->options->profileparamtype == PLOT4NM)
                sprintf (var, "4Nm_%li%li", from, to);
            else
                sprintf (var, "  M_%li%li", from, to);
        }
        if (which == world->numpop2 && world->options->gamma)
            sprintf (var, " Alpha  ");
#ifdef LONGSUM

        if(which >= len && world->options->fluctuate)
            sprintf(var," Rate_%li%li  ", (which - len) / world->numpop, ((which - len) % world->numpop)+1);
#endif /*LONGSUM*/

    }
}

void
calc_profile_likelihood (char method, long which, MYREAL *likes,
                         MYREAL **param, world_fmt * world, long nn,
                         nr_fmt * nr)
{
    long i, j;
    switch (method)
    {
    case 'p' /*percentiles */ :
    case 'q':   /*quick and dirty: no correlation between parametes */
    case 'f':   /*mixture of quick and precise */
        profile_max_percentiles (which, likes, param, world, nn, nr);
        break;
    case 's':   /* spline */
    case 'd' /*discrete */ :
        profile_max_spline (method, which, likes, param, world, nr, nn);
        break;
    }
    if (world->options->profileparamtype == PLOT4NM)
    {
        for (i = 0; i < GRIDSIZE; i++)
        {
            for (j = world->numpop; j < world->numpop2; j++)
            {
                param[i][j] *=
                    param[i][(long) (j - world->numpop) / (world->numpop - 1)];
            }
        }
    }
    //else
    //  printf ("           Using M (m/mu) instead of 4Nm\n");
}


void
prognose_profile_end (time_t starttime, long numpop2, long nn)
{
#ifndef NOTIME_FUNC
    static long sumtime = 0;
    static long tally = 0;
    char * nowstr;
    time_t endtime;
    time_t nowbin;
    struct tm *nowstruct;
    nowstr = (char *) mycalloc(STRSIZE,sizeof(char));
    
    if (tally % nn == 0)
    {
        tally++;
        time (&endtime);
        sumtime = (long) (endtime - starttime);
#ifdef MPI

        nowbin = starttime + (sumtime / tally * (nn * numpop2));
#else

        nowbin = starttime + (sumtime / tally * (nn * numpop2));
#endif

        nowstruct = localtime (&nowbin);

        strftime (nowstr, STRSIZE, "%H:%M %B %d %Y", nowstruct);
        FPRINTF (stdout, "           Prognosed end of run: %s\n", nowstr);
    }
    else
        tally++;
    myfree(nowstr);
#endif
}


void
profile_max_percentiles (long which, MYREAL *likes, MYREAL **param,
                         world_fmt * world, long nn, nr_fmt * nr)
{
    const MYREAL probabilities[] = GRID;

    long i;
    //, trials = 0;
    long df;
    long len;
    boolean QD = FALSE, QD2 = FALSE;
    MYREAL prob, normd, maxlike, testlike, minparam = 0, maxparam = 0;
    MYREAL *mlparam;
    long rep = world->options->replicate ? (world->loci > 1 ? 0 : world->repstop) : 0;
    df = (world->options->df > 0) ? world->options->df : set_df (which,
                world->options->
                custm2,
                world->numpop,
                world->numpop2 +
                (long) world->
                options->gamma);

    /* QD assumes that there is NO correlation between parameters
       that is certainly wrong but sometimes close and much faster
     */
    if (world->options->profilemethod == 'q')
        QD = TRUE;
    else
    {
        if (world->options->profilemethod == 'f')
        {
            QD = TRUE;
            QD2 = TRUE;
        }
        else
            QD = FALSE;
    }
    /* the minimum */
    minparam = which < world->numpop ? SMALLEST_THETA : SMALLEST_MIGRATION;
    /* the maximum */
    if (world->loci > 1)
    {
        mlparam = world->atl[rep][world->loci].param;
        maxlike = world->atl[rep][world->loci].param_like;
    }
    else
    {
        mlparam = world->atl[rep][world->locus].param;
        maxlike = world->atl[rep][world->locus].param_like;
    }

    len = world->numpop2 +  (long) world->options->gamma;
#ifdef LONGSUM
    len += world->numpop * 3;
#endif /*LONGSUM*/

    memcpy (param[nn - 1], mlparam, sizeof (MYREAL) * len);

    likes[nn - 1] = maxlike;
    maxparam = param[nn - 1][which];

    for (i = 0; i < nn - 1; i++)
    {
        if (probabilities[i] > 0.5)
            prob = 1 - probabilities[i];
        else
            prob = probabilities[i];
        testlike = maxlike - 0.5 * find_chi (df, prob);
        if (i > 0)
        {
            minparam = param[i - 1][which];
            if (minparam < 0)
                minparam = which < world->numpop ?
                           SMALLEST_THETA : SMALLEST_MIGRATION;
        }
#ifdef __MWERKS__
        eventloop ();
#endif

        len = world->numpop2 +  (long) world->options->gamma;
#ifdef LONGSUM

        len += world->numpop * 3;
#endif /*LONGSUM*/

        memcpy (param[i], mlparam, sizeof (MYREAL) * len);

        //xcode trials =
            find_profilelike (testlike, prob, i, minparam, maxparam, &param[i],
                              &likes[i], which, world, QD, QD2, mlparam, &normd,
                              nr);
        world->percentile_failed[which][i] = ((fabs(likes[i] - testlike)) > PERCENTILETOLERANCE);
        if (world->options->progress)
            prognose_profile_end (world->starttime, world->numpop2 +
                                  (long) world->options->gamma, nn);
    }
}

void
profile_max_spline (char method, long which, MYREAL *likes, MYREAL **param,
                    world_fmt * world, nr_fmt * nr, long nn)
{

    MYREAL deviate[] = DEVIATE;
//xcode    long errc;
    long indeks[] = INDEX;
    const MYREAL probabilities[] = GRID;
    long i;
    //, trials = 0;
    long panic, df;
    long len;
    MYREAL prob, normd, maxlike, testlike;
    //, maxparam = 0;
    //minparam = 0, 
    MYREAL *mlparam;
    MYREAL **hess;
    MYREAL tmp, x, xx, xlow, xhigh, mid, low, high, value;
    long rep =
        world->options->replicate ? (world->loci > 1 ? 0 : world->repstop) : 0;
    spline_fmt *spline;
    /* set the hessian for hte maximizer*/
    doublevec2d(&hess, world->numpop2 + (long) world->options->gamma,
		world->numpop2 + (long) world->options->gamma);
    /* set degree of freedom */
    df =
        (world->options->df > 0) ? world->options->df : set_df (which,
                world->options->
                custm2,
                world->numpop,
                world->numpop2 +
                (long) world->
                options->gamma);
    /* the minimum */
    //xcode minparam = which < world->numpop ? SMALLEST_THETA : SMALLEST_MIGRATION;
    /* the maximum */
    if (world->loci > 1)
    {
        mlparam = world->atl[rep][world->loci].param;
        maxlike = world->atl[rep][world->loci].param_like;
    }
    else
    {
        mlparam = world->atl[rep][world->locus].param;
        maxlike = world->atl[rep][world->locus].param_like;
    }
    len = world->numpop2 +  (long) world->options->gamma;
#ifdef LONGSUM

    len += world->numpop * 3;
#endif /*LONGSUM*/

    memcpy (param[nn - 1], mlparam, sizeof (MYREAL) * len);
    likes[nn - 1] = maxlike;
    //xcode maxparam = param[nn - 1][which];
    /* calculate likelihood at GRIDNODES */
    for (i = 0; i < nn - 1; i++)
    {
        if (which >= world->numpop && mlparam[which] < EPSILON)
            value = deviate[i];
        else
            value = deviate[i] * mlparam[which];
        memcpy (param[i], mlparam, sizeof (MYREAL) * len);
        nr->profilenum = set_profile_param (param[i], which, nr->profiles,
                                            value, nr->values, world);
        //do_profiles(nr->world, nr, likes, &normd, PROFILE,
        //          nr->world->rep, nr->world->repkind);
        maximize (&param[i], world, nr, hess, PROFILE, world->repkind);
        likes[i] = nr->llike;
        normd = nr->normd;
        memcpy (param[i], world->param0, sizeof (MYREAL) * nr->partsize);
        //xcode trials = 1;
        if (method != 's')
            prognose_profile_end (world->starttime,
                                  world->numpop2 + (long) world->options->gamma,
                                  nn);
    }
    if (method == 's')
    {
        spline = (spline_fmt *) mycalloc (1, sizeof (spline_fmt));
        setup_spline (spline);
        prepare_spline_nodes (nn, which, param, likes, indeks, spline->param,
                              spline->like);
        prepare_spline_first (spline->param, spline->like);
        prepare_spline_last (nn, spline->param, spline->like);
        calc_spline_diff (spline->diff, spline->param, spline->like, nn);
        calc_spline (spline->param, spline->like, nn, spline->constr,
                    spline->diff, spline->diff2, spline->diagn, spline->work,spline->nwork);
        for (i = 0; i < nn - 1; i++)
        {
            if (probabilities[i] < 0.5)
            {
                xhigh = spline->param[(nn + 2) / 2];
                high = spline->like[(nn + 2) / 2];
                xlow = spline->param[0];
                low = spline->like[0];
                prob = probabilities[i];
            }
            else
            {
                xlow = spline->param[nn + 1];
                low = spline->like[nn + 1];
                xhigh = spline->param[(nn + 2) / 2];
                high = spline->like[(nn + 2) / 2];
                prob = 1. - probabilities[i];
            }
            testlike = maxlike - 0.5 * find_chi (df, prob);
            x = (xhigh + xlow) / 2.;
            panic = 0;
            mid = testlike - 1000.;
            while (panic++ < MAX_PROFILE_TRIALS
                    && fabs (mid - testlike) > BIGEPSILON)
            {
                mid =
                    recalc_profilelike (x, spline->param, spline->like, nn,
                                        spline->diff, spline->diff2, spline->work,
                                        spline->nwork);
                if (mid < low || mid > high)
                    break;
                if (testlike < mid)
                {
                    high = mid;
                    tmp = x;
                    x = (xlow + tmp) / 2;
                    xhigh = tmp;
                }
                else
                {
                    low = mid;
                    tmp = x;
                    x = (tmp + xhigh) / 2;
                    xlow = tmp;
                }
            }
            xx = EXP (x);
            nr->profilenum = set_profile_param (param[i], which, nr->profiles,
                                                xx, nr->values, world);
            do_profiles (nr->world, nr, likes, &normd, PROFILE,
                         nr->world->rep, nr->world->repkind);
            //      maximize(param[i], world, nr, PROFILE, world->repkind);
            likes[i] = nr->llike;
            normd = nr->normd;
            memcpy (param[i], world->param0, sizeof (MYREAL) * nr->partsize);
            //xcode x = LOG (xx);
            prognose_profile_end (world->starttime,
                                  world->numpop2 + (long) world->options->gamma,
                                  nn);
        }
        destroy_spline (spline);
	myfree(hess[0]);
	myfree(hess);
    }
}

// setup spline: allocs all the necessary arrays and fill with
//               values when necessary.
void
setup_spline (spline_fmt * spline)
{
    const long degree = 3;
    long i;
    spline->ntab = 1;
    spline->nwork = (GRIDSIZE + 2) * 9 + 5 + degree * (degree + 11) / 2 + 9;
    spline->param = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (GRIDSIZE + 2));
    spline->like = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (GRIDSIZE + 2));
    spline->diff = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (GRIDSIZE + 2));
    spline->diff2 = (MYREAL *) mycalloc (1, sizeof (MYREAL) * (GRIDSIZE + 2));
    spline->constr = (long *) mycalloc (1, (GRIDSIZE + 2) * sizeof (long));
    spline->diagn = (long *) mycalloc (1, (GRIDSIZE + 2) * sizeof (long));
    spline->work = (MYREAL *) mycalloc (1, spline->nwork * sizeof (MYREAL));

    for (i = 0; i < GRIDSIZE + 2; i++)
        spline->constr[i] = 3;
}

void
destroy_spline (spline_fmt * spline)
{
    myfree(spline->param);
    myfree(spline->like);
    myfree(spline->diff);
    myfree(spline->diff2);
    myfree(spline->constr);
    myfree(spline->diagn);
    myfree(spline->work);
}

/* calc_spline() calculates the spline function derivatives,
   and adds two additional
   points, these two points bracket the real data points:
   so that param[0]=0 and param[n+1] = 1000 * the last real param,
   likes[0] is set -MYREAL_MAX and likes[n+1] is set using a linear interpolation
   using the last two real data points.
 
   Attention the arrays param, likes, and yy need to have a size
   of  n+2 and not n
 */
int
calc_spline (MYREAL *param, MYREAL *like, long nn, long *constr, MYREAL *diff,
             MYREAL *diff2, long *diagn, MYREAL *work, long nwork)
{
    /* call to spline.c routines */
    long degree = 3;
    long smoothness = 1;
    long opt = 414;  //convex, monotone, derivatives constraints
    MYREAL d0 = 0, dnp = 0, d20 = 0, d2np = 0; //not used
    MYREAL eps = 0.0001;
    long kmax = 10;  //not used
    long maxstp = 10;  //not used
    long errc;
    dbvssc_ (param, like, &nn, &degree, &smoothness, &opt, &d0, &dnp, &d20,
             &d2np, constr, &eps, NULL, NULL, NULL, NULL, &kmax, &maxstp, &errc,
             diff, diff2, diagn, work, &nwork);
    /* end call to spline.c routines */
    return (long) errc;
}

void
prepare_spline_nodes (long n, long which, MYREAL **rawparam, MYREAL *rawlikes,
                      long *indeks, MYREAL *param, MYREAL *likes)
{
    long i;
    for (i = 0; i < n; i++)
    {
        param[i + 1] = LOG (rawparam[indeks[i]][which]);
        likes[i + 1] = rawlikes[indeks[i]];
    }
}

void
prepare_spline_first (MYREAL *param, MYREAL *likes)
{
    const long i = 2;
    const long i0 = 1;
    param[0] = LOG (SMALLEST_THETA);

    likes[0] =
        ((likes[i] - likes[i0]) / (param[i] - param[i0]) * param[0]) +
        (likes[i0] + param[i0] * (likes[i] + likes[i0]) / (param[i] - param[i0]));
}

void
prepare_spline_last (long n, MYREAL *param, MYREAL *likes)
{
    long n0 = n - 1;
    long n1 = n + 1;
    param[n1] = 100000. + param[n];
    likes[n1] =
        ((likes[n] - likes[n0]) / (param[n] - param[n0]) * param[n1]) +
        (likes[n0] + param[n0] * (likes[n] + likes[n0]) / (param[n] - param[n0]));
}


void
calc_spline_diff (MYREAL *diff, MYREAL *param, MYREAL *likes, long n)
{
    long i;
    diff[0] = ((likes[1] - likes[0]) / (param[1] - param[0]));
    diff[n] = ((likes[n] - likes[n + 1]) / (param[n] - param[n + 1]));
    for (i = 1; i < n; i++)
    {
        diff[i] =
            ((likes[i + 1] - likes[i - 1]) / (param[i + 1] - param[i - 1]));
    }
    diff[GRIDMIDDLE] = 0.0;
}


/* recalc_profilelike uses the spline derivatives and returns a new
   loglikelihood value,
   testparam is THE LOG(testparam) !
 */
MYREAL
recalc_profilelike (MYREAL testparam, MYREAL *param, MYREAL *like, long nn,
                    MYREAL *diff, MYREAL *diff2, MYREAL *work, long nwork)
{
    long degree = 3;
    long smoothness = 1;
    long errc;
    long sbopt = 2;
    long ntab = 1;
    long yy0 = 1;
    long yy1 = 1;
    long yy2 = 1;
    long erre;
    MYREAL tmp;

    MYREAL *y0tab, *y1tab, *y2tab;
    //MYREAL *gaga;
    //long i,gagatab=10;
    //gaga= (MYREAL*)calloc(1,sizeof(MYREAL)*10);
    y0tab = (MYREAL *) mycalloc (1, sizeof (MYREAL) * 10);
    y1tab = (MYREAL *) mycalloc (1, sizeof (MYREAL) * 10);
    y2tab = (MYREAL *) mycalloc (1, sizeof (MYREAL) * 10);
    //gaga[0]= -5;
    //for(i=1;i<gagatab;i++)
    //gaga[i] =  gaga[i-1] + 0.5;

    dbvsse_ (param, like, &nn, &degree, &smoothness,
             //gaga, &gagatab,
             &testparam, &ntab, &sbopt, &yy0, &yy1, &yy2, &errc, diff, diff2,
             y0tab, y1tab, y2tab, &erre, work, &nwork);
    //for(i=0;i<gagatab;i++)
    //printf("%f %f\n",gaga[i],y0tab[i]);
    //exit(0);

    tmp = y0tab[0];
    myfree(y0tab);
    myfree(y1tab);
    myfree(y2tab);
    return tmp;

}


long
set_profile_param (MYREAL *param, long which, long *allwhich,
                   MYREAL xval, MYREAL *xvals, world_fmt * world)
{
    boolean has_fixed = FALSE;

    long z = 0, i = 0, zz = 0;//zz is the counter for # fixed elements
    //xcode zzz = 0; 
    long numpop = world->numpop;
    long numpop2 = world->numpop2;
    char *p;
    char *custm2 = world->options->custm2;
    long *fixed = NULL;
    long profnum = 0;
    MYREAL thetax, thetay;
    // find all real zero parameters and adds them to
    // the list of parameters not to  maximize.
    // [inefficient, then this will be called several times
    //  and recalculated for nothing]
    if (strpbrk (world->options->custm, "0c"))
    {
        has_fixed = TRUE;
        fixed = (long *) mycalloc (1, sizeof (long) * (numpop2 + (long) world->options->gamma));
        p = world->options->custm2;
        while (*p != '\0')
        {
            if ((*p == '0' || *p == 'c') && i != which)
            {
                fixed[zz]= i;
                zz++;
            }
            p++;
            i++;

        }
    }
    switch (custm2[which])
    {
    case '0':
        profnum = 0;
        break;
    case '*':
        param[which] = xval;
        allwhich[0] = which, xvals[0] = xval;
        if (has_fixed)
        {
            //xcode zzz = 0;
            addinvar (param, fixed, allwhich, xvals, 1, zz)
                ;
            myfree(fixed)
            ;
        }
        profnum = 1 + zz;
        break;
    case 'S':
        z = 0;
        while (z < world->options->sym2n &&
                which != world->options->sym2param[z][0] &&
                which != world->options->sym2param[z][1])
            z++;
        thetax = param[world->options->sym2param[z][2]];
        thetay = param[world->options->sym2param[z][3]];
        if (which == world->options->sym2param[z][0])
        {
            allwhich[0] = which;
            allwhich[1] = world->options->sym2param[z][1];
            param[world->options->sym2param[z][0]] = xval;
            param[world->options->sym2param[z][1]] = xval * thetax / thetay;
            xvals[0] = xval;
            xvals[1] = xval  * thetax / thetay;
        }
        else
        {
            allwhich[1] = world->options->sym2param[z][0];
            allwhich[0] = which;
            param[world->options->sym2param[z][0]] = xval * thetay / thetax;
            param[world->options->sym2param[z][1]] = xval;
            xvals[1] = xval * thetay / thetax;
            xvals[0] = xval;
        }
        if (has_fixed)
        {
            addinvar (param, fixed, allwhich, xvals, 2, zz)
                ;
            myfree(fixed)
            ;
        }
        profnum = 2 + zz;
        //      profnum = 1 + zz;
        break;
    case 's':
        z = 0;
        while (z < world->options->symn &&
                which != world->options->symparam[z][0] &&
                which != world->options->symparam[z][1])
            z++;
        param[world->options->symparam[z][0]] = xval;
        param[world->options->symparam[z][1]] = xval;
        if (which == world->options->symparam[z][0])
        {
            allwhich[0] = which;
            allwhich[1] = world->options->symparam[z][1];
        }
        else
        {
            allwhich[0] = world->options->symparam[z][0];
            allwhich[1] = which;
        }
        xvals[0] = xval;
        xvals[1] = xval;
        if (has_fixed)
        {
            addinvar (param, fixed, allwhich, xvals, 2, zz)
                ;
            myfree(fixed)
            ;
        }
        profnum = 2 + zz;
        break;
    case 'm':
        if (which < numpop)
        {
            for (i = 0; i < numpop; i++)
            {
                param[i] = xval;
                xvals[i] = xval;
                allwhich[i] = i;
            }
            if (has_fixed)
            {
                addinvar (param, fixed, allwhich, xvals, numpop, zz)
                    ;
                myfree(fixed)
                ;
            }
            profnum = numpop + zz;
            break;
        }
        else
        {
            for (i = numpop; i < numpop2 + (long) world->options->gamma; i++)
            {
                param[i] = xval;
                xvals[i - numpop] = xval;
                allwhich[i - numpop] = i;
            }
            profnum = numpop2 + (long) world->options->gamma - numpop;
            break;
        }
    }
    return profnum;
}

// adds the invariants to the profiles list so that we do not
// evaluate them and run unnecessary cycles in the maximizer
// fixed = the index of the fixed parameters
// allwhich = the variables in the profile list
// start = last element filled element in allwhich
// len = length of the fixed array
// PB 06/30/99
void
addinvar (MYREAL *param, long *fixed, long *allwhich, MYREAL *xvals,
          long start, long len)
{
    long i, z = 0;
    for (i = start; i < len + start; i++)
    {
        allwhich[i] = fixed[z]
                          ;
        xvals[i] = param[fixed[z]]
                   ;
        ++z;
    }

}

long
find_profilelike (MYREAL testlike, MYREAL prob, long whichprob, MYREAL minparam,
                  MYREAL maxparam, MYREAL **testparam, MYREAL *llike,
                  long which, world_fmt * world, boolean QD, boolean QD2,
                  MYREAL *mlparam, MYREAL *normd, nr_fmt * nr)
{
    const MYREAL probabilities[] = GRID;
    MYREAL tmp;
    MYREAL xlow, low;
    MYREAL x, mid = -MYREAL_MAX;
    MYREAL xhigh, high;
    long panic = 0;
    helper_fmt helper;
    MYREAL *ltestparam;
    MYREAL **hess;
    long len = world->numpop2 +  (long) world->options->gamma;
    doublevec2d(&hess,len,len);
#ifdef LONGSUM

    len += world->numpop * 3;
#endif /*LONGSUM*/
    world->locus = world->loci; // in case we have used gamma
    // used to set gamma option correctly, overly complicated, need
    // revamping, because the gamma test is all over the place PB022403
    
    ltestparam = (MYREAL *) mycalloc (len, sizeof (MYREAL));
    *normd = -1;
    prepare_broyden (PROFILE, world, &helper.multilocus);
    if (probabilities[whichprob] < 0.5)
    {
        xhigh = maxparam;
        xlow = whichprob > 0 ? minparam : SMALLEST_THETA / 10000.;
    }
    else
    {
        if (probabilities[whichprob] > 0.5)
        {
            xhigh = maxparam;
            if (maxparam < EPSILON)
                xlow = probabilities[whichprob] >= 0.98 ? 1000000.0 : minparam;
            else
                xlow =
                    probabilities[whichprob] >= 0.98 ? maxparam * 10000 : minparam;
        }
        else
        {
            xlow = minparam;
            xhigh = maxparam;
        }
    }


    // QD=TRUE -> calculate only CALCLIKE for parameter
    // QD2=TRUE -> assumes QD and then adds a maximizer step
    // if both are false maximize all parameters
    if (QD)
    {
        memcpy (*testparam, mlparam, sizeof (MYREAL) * len);
        nr->profilenum = set_profile_param (*testparam, which, nr->profiles,
                                            xlow, nr->values, world);
        (*testparam)[which] = xlow;
        set_logparam (ltestparam, *testparam, nr->partsize);
        //SETUPPARAM0 (world, nr, world->repkind, nr->repstart, nr->repstop,
        // world->loci, PROFILE, helper.multilocus);
        if(world->options->gamma)
            initgammacat (nr->categs, (*testparam)[nr->numpop2],1./* EXP (lparam[0])*/,
                          nr->rate, nr->probcat);
        fill_helper (&helper, *testparam, ltestparam, nr->world, nr);
        low = CALCLIKE (&helper, *testparam, ltestparam);
        memcpy (*testparam, mlparam, sizeof (MYREAL) *
                (world->numpop2 + (long) world->options->gamma));
        nr->profilenum = set_profile_param (*testparam, which, nr->profiles,
                                            xhigh, nr->values, world);
        (*testparam)[which] = xhigh;
        set_logparam (ltestparam, *testparam, nr->partsize);
        if(world->options->gamma)
            initgammacat (nr->categs, (*testparam)[nr->numpop2],1./* EXP (lparam[0])*/,
                          nr->rate, nr->probcat);
        fill_helper (&helper, *testparam, ltestparam, world, nr);
        high = CALCLIKE (&helper, *testparam, ltestparam);
        //    testlike = high - 0.5 * find_chi (1, prob);
    }
    else
    {
        memcpy (*testparam, mlparam, sizeof (MYREAL) * len);
        nr->profilenum = set_profile_param (*testparam, which, nr->profiles,
                                            xlow, nr->values, world);
        maximize (testparam, world, nr, hess, PROFILE, world->repkind);
        low = nr->llike;
        memcpy (*testparam, mlparam, sizeof (MYREAL) *len);
        nr->profilenum = set_profile_param (*testparam, which, nr->profiles,
                                            xhigh, nr->values, world);
        maximize (testparam, world, nr, hess, PROFILE, world->repkind);
        high = nr->llike;

        // testlike = high - 0.5 * find_chi (1, prob);

        //if (fabs(high - world->atl[rep][loc].param_like)>EPSILON*10)
        //    printf
        //    ("%i> PROBLEM?: high=%f ml=%f (high=?=ml) xigh=%f mlparam[which]=%f testlike=%f\n",
		//  myID, high, world->atl[rep][loc].param_like, xhigh,
  // mlparam[which], testlike);
    }
    panic = 0;
    x = (xhigh + xlow) / 2.;
    while (panic++ < MAX_PROFILE_TRIALS && fabs (mid - testlike) > BIGEPSILON)
    {
        if (QD)
        {
            memcpy (*testparam, mlparam, sizeof (MYREAL) * len);
            nr->profilenum = set_profile_param (*testparam, which,
                                                nr->profiles, x,
                                                nr->values, world);
            (*testparam)[which] = x;
            if(nr->world->options->gamma)
                initgammacat (nr->categs, (*testparam)[nr->numpop2],1.,nr->rate, nr->probcat);
            set_logparam (ltestparam, *testparam, nr->partsize);
            fill_helper (&helper, *testparam, ltestparam, world, nr);
            mid = CALCLIKE (&helper, *testparam, ltestparam);
        }
        else
        {
            memcpy (*testparam, mlparam, sizeof (MYREAL) * nr->partsize);
            nr->profilenum = set_profile_param (*testparam, which, nr->profiles,
                                                x, nr->values, world);
            maximize (testparam, world, nr, hess, PROFILE, world->repkind);
            mid = *llike = nr->llike;
            *normd = nr->normd;
        }
        if ((testlike < low) || (testlike > high))
            return MAX_PROFILE_TRIALS;
        if (testlike < mid)
        {
            high = mid;
            tmp = x;
            x = (xlow + tmp) / 2;
            xhigh = tmp;
        }
        else
        {
            low = mid;
            tmp = x;
            x = (tmp + xhigh) / 2;
            xlow = tmp;
        }
    }
    if (QD)
    {

        *llike = mid;
        if (QD2)
        {
            x = (*testparam)[which];
            memcpy (*testparam, mlparam,
                    sizeof (MYREAL) * (nr->numpop2 +
                                       (long) world->options->gamma));
            nr->profilenum =
                set_profile_param (*testparam, which, nr->profiles, x, nr->values,
                                   world);
            maximize (testparam, world, nr, hess, PROFILE, world->repkind);
            *llike = nr->llike;
            *normd = nr->normd;
            memcpy (*testparam, world->param0, sizeof (MYREAL) * nr->partsize);
        }
    }
    myfree(ltestparam);
    myfree(hess[0]);
    myfree(hess);
    return panic;
}


MYREAL
interpolate (MYREAL xl, MYREAL xc, MYREAL xh, MYREAL low, MYREAL center,
             MYREAL high, MYREAL testlike)
{
    MYREAL xc2, xl2, xh2, a, b, c, x;

    xc2 = xc * xc;
    xl2 = xl * xl;
    xh2 = xh * xh;

    a =
        (-(high * xc) + low * xc + center * xh - low * xh - center * xl +
         high * xl) / ((-xc + xh) * (xh - xl) * (-xc + xl));

    b =
        (high * xc2 - low * xc2 - center * xh2 + low * xh2 + center * xl2 -
         high * xl2) / ((-xc + xh) * (xc - xl) * (-xh + xl));

    c =
        -((low * xc2 * xh - low * xc * xh2 - high * xc2 * xl + center * xh2 * xl +
           high * xc * xl2 - center * xh * xl2) / ((-xc + xh) * (xc - xl) * (xh -
                                                   xl)));


    x = (-b - sqrt (b * b - 4 * a * (c - testlike))) / (2 * a);

    return x;
}


void
allocate_profile_percentiles (world_fmt * world)
{
    long i;
    long len = world->numpop2 + (world->options->gamma ? 1 : 0);
#ifdef LONGSUM
    len += world->numpop*3;
#endif /*LONGSUM*/

    world->quantiles =
        (quantile_fmt *) mycalloc (1, sizeof (quantile_fmt) *len);
    world->percentile_failed = (boolean **) mycalloc (len, sizeof (boolean *));
    for (i = 0; i < len; i++)
    {
        world->quantiles[i].name = (char *) mycalloc (1, sizeof (char) * 20);
        world->quantiles[i].param =
            (MYREAL *) mycalloc (1, sizeof (MYREAL) * GRIDSIZE);
        world->percentile_failed[i] =
            (boolean *) mycalloc (GRIDSIZE, sizeof (boolean));
    }
}

///
/// free profile related quantitites: quantiles, and indicator for failed
/// convergence to the percentiles
void
destroy_profile_percentiles (world_fmt * world)
{
    long i;
    long len = world->numpop2 + (world->options->gamma ? 1 : 0);
#ifdef LONGSUM

    len += world->numpop*3;
#endif /*LONGSUM*/

    for (i = 0; i < len; i++)
    {
        myfree(world->quantiles[i].name);
        myfree(world->quantiles[i].param);
        myfree(world->percentile_failed[i]);
    }
    myfree(world->quantiles);
    myfree(world->percentile_failed);
}


void
print_menu_profile (long which, long nn)
{
    char * nowstr;
    nowstr = (char*) mycalloc(LINESIZE,sizeof(char));
    
    get_time (nowstr, "%H:%M:%S");
    FPRINTF (stdout,
             "%s   Calculating profile likelihoods for parameter %2li out of %2li\n",
             nowstr, 1 + which, nn);
    myfree(nowstr);
}


long
set_df (long which, char *custm2, long numpop, long numpop2)
{
    char ch;
    long i;
    long zeros = 0;
    long ms = 0;
    long ss = 0;

    for (i = 0; i < numpop2; i++)
    {
        ch = custm2[i];
        switch (ch)
        {
        case '0':
            zeros++;
            break;
        case 'm':
            ms++;
            break;
        case 's':
            ss++;
            break;
        }
    }
    if (numpop2 - zeros == ms)
        return 2;
    else
    {
        if (custm2[0] == 'm')
            return numpop2 - numpop + 1 - zeros - ss / 2;
        else
            return numpop2 - zeros - ss / 2;
    }
    return -1;
}

#ifdef PRIVATE
void master_gridder(world_fmt *world, long *gmaxptr)
{
    FILE *alldim;
    nr_fmt *nr;
    long rep =     world->options->replicate ? (world->loci > 1 ? 0 : world->repstop) : 0;
    MYREAL *mlparam;
    MYREAL *testparam;
    MYREAL *ltestparam;
    nr = (nr_fmt *) mycalloc (1, sizeof (nr_fmt));
    create_nr (nr, world, *gmaxptr, 0, world->loci,
               world->repkind, world->rep);
    testparam = (MYREAL *) mycalloc ( nr->partsize, sizeof (MYREAL));
    ltestparam = (MYREAL *) mycalloc ( nr->partsize, sizeof (MYREAL));

    // specify the ML parameters
    if (world->loci > 1)
    {
        mlparam = world->atl[rep][world->loci].param;
        //maxlike = world->atl[rep][world->loci].param_like;
    }
    else
    {
        mlparam = world->atl[rep][world->locus].param;
        //maxlike = world->atl[rep][world->locus].param_like;
    }
    memcpy (testparam, mlparam, sizeof (MYREAL) *
            (world->numpop2 + (long) world->options->gamma));
    set_logparam (ltestparam, testparam, nr->partsize);
    alldim = fopen("dimensions","w+");
    calc_grid(alldim, world, testparam, ltestparam, mlparam, nr, 0);
    fclose(alldim);
    myfree(testparam);
    myfree(ltestparam);
    destroy_nr(nr, world);
}

void calc_grid(FILE *alldim, world_fmt * world,  MYREAL *testparam, MYREAL *ltestparam, MYREAL *mlparam, nr_fmt *nr, long which)
{
    long i;
    MYREAL deviate[] = DEVIATE2;
    //MYREAL maxlike;
    helper_fmt helper;
    MYREAL saveparam, lsaveparam;
    // for(j=0; j<nr->numpop2; j++)
    //   {
    if(which<world->numpop2)
    {
        saveparam = testparam[which];
        lsaveparam = ltestparam[which];
        for(i=0; i<GRIDSIZE; i++)
        {
            testparam[which] = deviate[i] * mlparam[which];
            ltestparam[which] = log(testparam[which]);
            if(which==nr->numpop2-1)
            {
                if(world->options->gamma)
                    initgammacat (nr->categs, testparam[nr->numpop2],1., nr->rate, nr->probcat);
                fill_helper (&helper, testparam, ltestparam, world, nr);
                //FPRINTF(alldim,"%li %li ",which,i);
                //for(z=0;z<world->numpop2;z++)
                //   FPRINTF(alldim, "%g ", testparam[z]);
                FPRINTF(alldim, "%f ", CALCLIKE (&helper, testparam, ltestparam));
            }
            calc_grid(alldim, world,testparam,ltestparam, mlparam, nr,which+1);

        }
        FPRINTF(alldim,"\n");
        testparam[which] = saveparam;
        ltestparam[which] =lsaveparam;
    }
    //   }
    //      FPRINTF(world->mathfile, "\n");
}
#endif



