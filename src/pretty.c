/*
 *  pretty.c
 *  migrate-n
 *
 *  Created by Peter Beerli on 7/25/05.
 *  Copyright 2005-2010 Peter Beerli. All rights reserved.
 *
 */
#ifdef PRETTY
#include "pretty.h"
#include "data.h"
#include "options.h"
#include "migevents.h"
#include "menu.h"
#include "bayes.h"

#include <stdarg.h>

#define LINESTRETCH 15
#define NOT_PERCENTILES FALSE
#define PERCENTILES TRUE
#define SQUARE 0
#define DIAMOND 1

#ifndef A4PAPER
#define LETTERADJUST 50
#else
#define LETTERADJUST 0
#endif


extern pdf_doc doc;
extern pdf_page page;
extern pdf_contents canvas;
extern int page_counter;
extern char pdf_pagetitle[LINESIZE+1];
extern char pdf_time[LINESIZE+1];
extern float page_height;
extern float left_margin;
extern int numcpu;
pdf_contents firstcanvas;
float timestampy;
void pdf_printf_right(float x, float y, char string[],...);
int pdf_print_header(char *title);
void pdf_print_contents_at(float x, float y, char *title);
void pdf_draw_line(float xs, float ys, float xe, float ye);
void pdf_migrate_logo(float x, float y, float stretch);
void pdf_migrate_logo_lines(float x, float y, float stretch);
void pdf_print_time(float startx, float *pageheight, char text[]);
boolean  pdf_advance(float *page_height);
void pdf_print_seqdata (float *page_height, float margin, world_fmt * world, data_fmt * data, option_fmt * options);
void pdf_printf_next(float x, float *y, char string[], ...);
void pdf_printf_ralign(float rx, float y, char string[], ...);
void pdf_printf(float x, float y, char align, char string[], ...);
void pdf_print_tableline(float *page_height, float *width, char *fmt, ...);
void pdf_print_result_param (float *lx, float *page_height,  MYREAL *param, long numpop, long pop,
                             boolean usem);
void pdf_table(float left_margin, float *page_height, int cols, int rows, char **buffer, char *position, int col_overflow, float separator);

long nice_element(MYREAL param, char *element, MYREAL lower, MYREAL mid, MYREAL upper, 
                  int low_mid_digits, int mid_upper_digits, char delimiter);
void pdf_linedotplot(long n, float *x, float *y, float dotcolor[], float linecolor[], MYREAL linethickness, float * page_height, float width, float height, boolean has_dots, boolean has_line);

void pdf_print_symbol(float lx, float *page_height, long fontsize, char pos, char *symbolstring);


#ifndef HAVE_STRSEP
/*-
* Copyright (c) 1990, 1993
 *      The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *      This product includes software developed by the University of
 *      California, Berkeley and its contributors.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
            * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

// #include "config.h"

// #include <sys/types.h>

// #include <string.h>
// #include <stdio.h>

// #ifndef HAVE_STRSEP
char *strsep(char **, const char *);
// #endif

// #if defined(LIBC_SCCS) && !defined(lint)
// static char sccsid[] = "@(#)strsep.c    8.1 (Berkeley) 6/4/93";
// #endif /* LIBC_SCCS and not lint */
// #ifndef lint
// static const char rcsid[] =
//   "$FreeBSD: src/lib/libc/string/strsep.c,v 1.2.12.1 2001/07/09 23:30:07 obrien Exp $";
// #endif

/*
 * Get next token from string *stringp, where tokens are possibly-empty
 * strings separated by characters from delim.
 *
 * Writes NULs into the string at *stringp to end tokens.
 * delim need not remain constant from call to call.
 * On return, *stringp points past the last NUL written (if there might
                                                         * be further tokens), or is NULL (if there are definitely no more tokens).
 *
 * If *stringp is NULL, strsep returns NULL.
 */
char *
strsep(register char **stringp, register const char *delim)
{
    register char *s;
    register const char *spanp;
    register int c, sc;
    char *tok;
    
    if ((s = *stringp) == NULL)
        return (NULL);
    for (tok = s;;) {
        c = *s++;
        spanp = delim;
        do {
            if ((sc = *spanp++) == c) {
                if (c == 0)
                    s = NULL;
                else
                    s[-1] = 0;
                *stringp = s;
                return (tok);
            }
        } while (sc != 0);
    }
    /* NOTREACHED */
}
#endif

///
/// returns the quantile at prob1 or prob2, if the value at prob2 is smaller than
/// the value at prob1 then the value at prob1 is return otherwise the value at prob2.
/// We assume range1< range2 
float quantiler(MYREAL *values, float prob1, float prob2, long range1, long range2)
{
  MYREAL *temp;
  MYREAL a,b;
  temp = (MYREAL *) calloc(range2+1,sizeof(MYREAL));
  memcpy(temp, values,sizeof(MYREAL)*range2);
  qsort ((void *) temp, range2, sizeof (MYREAL), numcmp);
  a = temp[((long) (range1 * prob1))];
  b = temp[((long) ((range2-1) * prob2))];
  myfree(temp);
  return (a<b ? b : a);
}

void symbol_Theta(float lx, float ly, int size, long subscript)
{
    char *thetatitle = "Q";
    char tempstring[100];
    int subsize =  (int) (0.75 * size);
    float lxdelta = lx + (float) 0.8 * size;
    float lydelta = ly - (float) 0.5 * subsize ;
    if(subscript >= 0)
      sprintf(tempstring,"%li",subscript);
    else
      sprintf(tempstring," ");
    pdf_contents_set_font_and_size(canvas, "Symbol", (float) size);
    pdf_print_contents_at(lx,ly,thetatitle);
    pdf_contents_set_font_and_size(canvas, "Symbol", (float) subsize);
    pdf_print_contents_at(lxdelta,lydelta,tempstring);
    // backt to default
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
}

void symbol_R(float lx, float ly, int size, long subscript)
{
    char *thetatitle = "R";
    char tempstring[100];
    int subsize =  (int) (0.75 * size);
    float lxdelta = lx + (float) 0.8 * size;
    float lydelta = ly - (float) 0.5 * subsize ;
    if(subscript > -1)
        sprintf(tempstring,"%li",subscript);
    else
      {
	if(subscript < -1)
	  sprintf(tempstring,"combined");
	else
	  sprintf(tempstring," ");
      }
    pdf_contents_set_font_and_size(canvas, "Helvetica", size);
    pdf_print_contents_at(lx,ly,thetatitle);
    if(subscript < -1)
      pdf_contents_set_font_and_size(canvas, "Helvetica", subsize);
    else
      pdf_contents_set_font_and_size(canvas, "Symbol", subsize);
    pdf_print_contents_at(lxdelta,lydelta,tempstring);
    // backt to default
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
}

void symbol_M(float lx, float ly, int size, long subscript1, long subscript2, boolean usem)
{
    char *migt1 = "M";
    char *migt2 = "xNm";
    char *migtitle = usem ? migt1 : migt2;
    int msub = usem ?  size/2 : size+ size/2;
    char tempstring[100];
    int subsize =  (int) (0.75 * size);
    float lxdelta = lx + (float) 0.8 * size + msub;
    float lydelta = ly - (float) 0.5 * subsize ;
    sprintf(tempstring,"%li->%li",subscript1,subscript2);
    pdf_contents_set_font_and_size(canvas, "Helvetica", size);
    pdf_print_contents_at(lx,ly,migtitle);
    pdf_contents_set_font_and_size(canvas, "Symbol", subsize);
    pdf_print_contents_at(lxdelta,lydelta,tempstring);
    // backt to default
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
}

///
/// Draw a simple line from x0,y0 to x0/y0, this function does not change
/// thickness or line type
void pdf_draw_line(float xs, float ys, float xe, float ye)
{
    pdf_contents_move_to(canvas, xs, ys);
    pdf_contents_line_to(canvas, xe, ye);
    pdf_contents_stroke(canvas);
}

///
/// draw a rectangle [from haru examples]
void draw_rect(pdf_contents canvas, float x, float y, const char* label)
{
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, x, y - 10.);
    pdf_contents_show_text(canvas, label);
    pdf_contents_end_text(canvas);
    
    pdf_contents_rectangle(canvas, x, y - 40., 220., 25.);
}

///
/// when there are fewer than need_pixels create a new page and 
/// and push page_height down so that the need_pixel object will fit.
float pdf_page_advance_or_not(float *page_height, float need_pixels)
{
    boolean new_page = FALSE;
    if(*page_height < need_pixels)
      {
        pdf_new_page("");
        *page_height = pdf_contents_get_height(canvas);
        *page_height -= 55. + LINESTRETCH; 
        new_page = TRUE;
      }
    return new_page;
}

///
/// advance a line and check whether we are at the end of the page, if yes then add a new page
boolean pdf_advance(float *page_height)
{
    *page_height -= LINESTRETCH;
    return pdf_page_advance_or_not(page_height, 55.);
}
///
/// advance a half-line and check whether we are at the end of the page, if yes then add a new page
boolean pdf_advance_half(float *page_height)
{
    *page_height -= LINESTRETCH/2.;
    return pdf_page_advance_or_not(page_height, 55.);
}

#define HORIZONTAL 0
#define VERTICAL   1 
///
/// plot a single tick with label at location x/y.
/// \param xs  x coordinate in true paper coordinates (points)
/// \param xy  y coordinate in true paper coordinates (points)
/// \param orientation either HORIZONTAL or VERTICAL, on HORIZONTAL=0 is checked
/// \param ticklength length of a the tick to draw
/// \param value this is a float value to put into the label
/// \param digits number of decimal digits
void   pdf_draw_tick(float xs, float ys, int orientation, float ticklength, float value, int digits)
{
    float w;
    float h;
    float tl = ticklength;
    char *title;
    title = (char *) mycalloc(100,sizeof(char));
    sprintf(title,"%.*f",digits, value);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    h = (float) pdf_contents_get_font_size(canvas);
    if(orientation==HORIZONTAL)
      {
        pdf_print_contents_at(xs-w/2, ys-tl-h, title);
        pdf_draw_line(xs,ys, xs , ys - tl);
      }
    else
      {
        pdf_print_contents_at(xs-tl-w-2, ys-h/2+2, title);
        pdf_draw_line(xs,ys, xs-tl , ys);
      }
    myfree(title);
}

///
/// find pretty values for tick labels
/// using heuristics
float prettytick(float mi, float ma, long ticks, 
		 float *nmi, float *nma, float * newdelta, int * digits)
{
  float newvalue = 0;
  float dist = (ma - mi)/ticks;
  float multiplier[5];
  long  lmi = (long) floor(log10(dist));
  float scale  = pow(10.,lmi);
  long i=0;
  multiplier[0] = 1.;
  multiplier[1] = 2.;
  multiplier[2] = 2.5;
  multiplier[3] = 5.;
  multiplier[4] = 10.;
  if(lmi<3)
    *digits = labs(lmi);
  else
    *digits = 0;
  dist /= scale;
  if(dist < 1.)
    dist *= 10.;
  else
    dist /= 10.;
  if(dist >= 10.)
    dist /= 10.;
  else
    dist *= 10.;
  i = 0;
  while(i < 5 && dist > multiplier[i])
    {
      i++;
    }
  *newdelta = scale * multiplier[i];  
  *nmi = ceil(mi/(*newdelta) - 0.05) * (*newdelta);

  return newvalue;
}

///
/// create a vertical or horizontal axis
/// \param up is either HORIZONTAL=1 or VERTICAL
/// \param xmin  real start value
/// \param xmax  real end value
/// \param ticks number of ticks
/// \param lx    x paper coordinate to start axis
/// \param ly    y paper coordinate to start axis
/// \param length  length on paper in paper coordinate system
void  pdf_create_axes(int up, float xmin, float xmax, long ticks, float lx, float ly, float length)
{
    long i;
    int digits = 3;
    float value = 0.0;
    float plotdelta = length/ticks;
    float realdelta = (xmax - xmin)/ticks;
    float newxmin;
    float newxmax;
    float newrealdelta;
    float newlx;
    float newly;
    // pretty printing of labels
    // find pretty tickmarks
    prettytick(xmin,xmax, ticks, &newxmin, &newxmax, &newrealdelta,&digits);
    // 
    if(up!=HORIZONTAL)
      {
	//xcode newlx = lx;
	newly = ly + (newxmin-xmin)/(xmax-xmin) * length;
	plotdelta *= newrealdelta/realdelta; 
        pdf_draw_line(lx,ly, lx ,ly + length);
        for(i=0; i<ticks; i++)
          {
            value = newxmin+i*newrealdelta;
	    if(value>xmax)
	      break;
            pdf_draw_tick(lx, newly+i*plotdelta, VERTICAL, 3, value, digits);
            }
          }
    else
      {
        pdf_draw_line(lx,ly, lx + length, ly);
	newly = ly;
	newlx = lx + (newxmin-xmin)/(xmax-xmin) * length;;
	plotdelta *= newrealdelta/realdelta; 

        for(i=0; i<ticks; i++)
          {
            value = newxmin+i*newrealdelta;
	    if(value>xmax)
	      break;
            pdf_draw_tick(newlx + i*plotdelta, 
                          newly, HORIZONTAL, 3, value, digits);
          }
      }
}

///
/// plot a dot of the shape i (currently i = square or diamond) at position
/// xs and ys in RGB color is a float vector of 3 values
void pdf_print_dot(float xs, float ys, float width, int shape, float color[])
{
    float w = (float) width / 2.;
    float red = color[0];
    float green = color[1];
    float blue = color[2];
    pdf_contents_set_rgb_fill(canvas, PdfRGBColor(red, green, blue));
    switch(shape)
      {
        case SQUARE:
            pdf_contents_move_to(canvas, xs-w, ys-w);
            pdf_contents_line_to(canvas, xs+w, ys-w);
            pdf_contents_line_to(canvas, xs+w, ys+w);
            pdf_contents_line_to(canvas, xs-w, ys+w);
            pdf_contents_line_to(canvas, xs-w, ys-w);
            break;
        case DIAMOND:
        default:
            pdf_contents_move_to(canvas, xs, ys + w);
            pdf_contents_line_to(canvas, xs+w, ys);
            pdf_contents_line_to(canvas, xs, ys-w);
            pdf_contents_line_to(canvas, xs-w, ys);
            pdf_contents_line_to(canvas, xs, ys + w);
            break;
      }
    pdf_contents_fill(canvas);
    pdf_contents_set_rgb_fill(canvas, PdfRGBColor(0, 0, 0));
}

void findmoments(MYREAL *vals, long n, float *p01, long *pos01, float *p99, long *pos99, float *mm, long *posmm)
{
    long i;
    float val;
    float total = 0. ;
    *mm = - (float) (LONG_MAX);
    *p01 = 0.;
    *p99 = 0.;
    *pos01 = 0;
    *pos99 = 0;
    *posmm = 0;
    for(i=0;i<n; i++)
      {
        val = (float) vals[i];
        total += val;
        if( *mm < val)
          {
            *mm = val;
            *posmm = i;
          }
      }
    val = 0.;
    for(i=0;i<n; i++)
      {
        val += (float) vals[i]/total;
        if(val < 0.01)
          {
            //	    if(*p01  < (float) vals[i])
            //  {
            *p01 = (float) vals[i];
            *pos01 = i;
            //  }
          }
        if(val <= 0.99)
          {
            //	    if(*p99 < (float) vals[i])
            //  {
            *p99 = (float) vals[i];
            *pos99 = i;
            //   }
          }
      }
}

///
/// finds outliers in vals and retunrs some specified quantiles and the position of these
/// values in the vals vector.
void findoutliers(MYREAL *vals, long n, 
                  float *p01, long *pos01, 
                  float *p99, long *pos99, 
                  float *mm, long *posmm)
{
    long i;
    MYREAL * temp;
    if(n==0)
      {
        warning("findoutliers(): The number of values to check for the histogram is zero\n");
        return;
      }
    temp = (MYREAL *) mycalloc((n+1),sizeof(MYREAL));
    memcpy(temp,vals,n*sizeof(MYREAL));
    qsort(temp, n, sizeof(MYREAL), numcmp);
    
    *mm = (float) temp[n-1];
    *p01 = (float) temp[(long)(n * 0.01)];
    *p99 = (float)temp[(long)(n * 0.99)];
    *pos01 = 0;
    *pos99 = 0;
    *posmm = 0;
    for(i=0;i<n; i++)
      {
        if(vals[i] <= *p01)
          {
            *pos01 = i;
          }
        if(vals[i] < *p99)
          {
            *pos99 = i;
          }
        
        if(vals[i] < *mm)
          {
            *posmm = i;
          }
        
      }
    myfree(temp);
}

///
/// only for fractions 0..1
void findminmax_freq_only(MYREAL *vals, long n, float *p00, long *pos00, float *p100, long *pos100)
{
    long i;
    float val;
    float total = 0. ;
    *p00 = 0.;
    *p100 = 0.;
    *pos00 = 0;
    *pos100 = 0;
    for(i=0;i<n; i++)
      {
        val = (float) vals[i];
        total += val;
      }
    val = 0.;
    for(i=0;i<n; i++)
      {
        val += (float) vals[i]/total;
        if(val < 0.000001)
          {
            *p00 = (float) vals[i];
            *pos00 = i;
          }
        if(val < 0.999999)
          {
            *p100 = (float) vals[i];
            *pos100 = i;
          }
      }
}

///
/// create a histogram from bincounts at a specific location and width and height
/// If binmax is -9999 then the maximum is reset to the 99% percentile.
/// If binmax is -999 then the maximum is reset to the 100% percentile
/// if nofreq is TRUE then the y axes is treated literally and not scaled
void pdf_histogram(MYREAL *binvals, char *set50, char *set95, long bins, float bindelta, float binmin, float binmax, float lx, float ly, float width, float height, boolean nofreq)
{
    long i;
    float total;
    float binvalsmax;
    long binvalspos;
    float p99;
    long pos99;
    float p100;
    long pos100;
    float p00;
    long pos00;
    float p01;
    long pos01;
    float x = 0.;
    float delta = width / bins;
    long numbins = bins;
    
    float red[3]={0.99,0.,0.};
    
    if(nofreq)
      {
        findoutliers(binvals,bins,&p01,&pos01, &p99,
                     &pos99, &binvalsmax, &binvalspos);
      }
    else
        findmoments(binvals,bins,&p01,&pos01, &p99,
                    &pos99, &binvalsmax, &binvalspos);
    //fprintf(stdout,"p01=%f, pos01=%li, p99=%f, pos99=%li, binvalsmax=%f, binvalspos=%li\n",p01,pos01, p99,
    //	    pos99, binvalsmax, binvalspos);
    total = 0. ;
    for(i=0;i<bins;i++)
      {
        total += (float) binvals[i];
      }
    
    if(binmax < -9000)
      {
        binmax = binmin + pos99 * bindelta;
        delta = width/pos99;
        numbins = pos99;
      }
    else
      {
        if(binmax < -900)
          {
            findminmax_freq_only(binvals,bins,&p00,&pos00, &p100,
                       &pos100);
            binmax = binmin + pos100 * bindelta;
            delta = width/pos100;
            numbins = pos100;
          }
      }
    if(nofreq)
        pdf_create_axes(VERTICAL, 0, p99, 5, lx, ly, height);
    else
        pdf_create_axes(VERTICAL, 0, binvalsmax/total, 5, lx, ly, height);
    pdf_create_axes(HORIZONTAL, binmin , binmax, 5, lx, ly, width);
    pdf_contents_set_line_width(canvas, delta);
    
    x = delta/1.95;
    for(i=0;i<numbins;i++)
      {
        if(set50[i] == '1')
            pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0.1, 0.1, 0.1));
        else
          {
            if(set95[i] == '1')
                pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0.4, 0.4, 0.4));
            else
                pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0.8, 0.8, 0.8));
          }
        if(nofreq)
          {
            if(binvals[i]/p99 > 1.0)
              {
                pdf_draw_line(lx+x,ly+0,lx+x,ly + height);
                pdf_print_dot(lx+x, ly + height - 3.0, 3, SQUARE, red);
              }
            else
                pdf_draw_line(lx+x,ly+0,lx+x,ly + binvals[i]/p99 * height);
          }
        else
            pdf_draw_line(lx+x,ly+0,lx+x,ly + binvals[i]/binvalsmax * height);
        //fprintf(stderr,"%f %f %f\n",lx+x,ly+ binvals[i]/binvalsmax * height, binvals[i]);
        x += delta;
      }
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    pdf_contents_set_line_width(canvas, 1);
}

///
/// create a histogram from bincounts at a specific location and width and height using *std as +-std dev
/// If binmax is -9999 then the maximum y-value printed is reset to the 99% percentile.
/// If binmax is -999 then the maximum y-value printed is reset to the 100% percentile
/// if nofreq is TRUE then the y axes is treated literally and are not scaled
/// world is needed for the skyline plots so that we can plot the dots for
/// the data points onto the graph.
void pdf_histogram_plus(MYREAL *binvals, MYREAL *std,  
                        char *set50, char *set95, 
                        long bins, float bindelta, float binmin, float binmax, 
                        float lx, float ly, float width, float height, 
                        MYREAL valmax, boolean nofreq, world_fmt * world, float * confidence, long topop)
{
    long locus;
    long ind;
    long pop;
    float ratio;
    float xcoord;
    float averagerate=0.;
    long i;
    long nbmin;
    float total;
    float binvalsmax;
    long binvalspos;
    float p99;
    long pos99;
    float p100;
    long pos100;
    float p00;
    long pos00;
    float p01;
    long pos01;
    float x = 0.;
    float delta = width / bins;
    long numbins = bins;
    //long grouping = 0;
    float r;
    float vold=0.;
    float xold;
    MYREAL *upper;
    MYREAL *lower;
    MYREAL l;
    MYREAL v;
    MYREAL u;
    //    long k;
    float blue[3]={0.,0.,0.99};
    float red[3]={0.99,0.,0.};

    upper = (MYREAL*) mycalloc(numbins,sizeof(MYREAL));
    lower = (MYREAL*) mycalloc(numbins,sizeof(MYREAL));
    // setting lower and upper to the mean values;
    memcpy(upper, binvals,numbins * sizeof(MYREAL));
    memcpy(lower, binvals,numbins * sizeof(MYREAL));
    
    for(i=0; i< numbins; i++)
      {
        upper[i] += fabs(std[i]);
        lower[i] -= fabs(std[i]);
        if(lower[i] < 0.0)
            lower[i] = 0.0;
      }
    if(nofreq)
      {	
        findoutliers(binvals,bins,&p01,&pos01, &p99,
                     &pos99, &binvalsmax, &binvalspos);
      }
    else
        findmoments(upper,bins,&p01,&pos01, &p99,
                    &pos99, &binvalsmax, &binvalspos);
    // fprintf(stdout,"p01=%f, pos01=%li, p99=%f, pos99=%li, binvalsmax=%f, binvalspos=%li\n",p01,pos01, p99,
    //	    pos99, binvalsmax, binvalspos);
    total = 0. ;
    for(i=0;i<bins;i++)
      {
        total += (float) upper[i];
      }
    
    if(binmax < -9000)
      {
        //find the 99% percentile of the upper limit to show the whole skyline
        //the closeup should not use this setting because it may be dominated by
        //few spikes
        findoutliers(upper,bins,&p01,&pos01, &p99,
                     &pos99, &binvalsmax, &binvalspos);
        binmax = binmin + pos99 * bindelta;
        delta = width/pos99;
        numbins = pos99;
      }
    else
      {
        if(binmax < -900)
          {
            findminmax_freq_only(upper,bins,&p00,&pos00, &p100,
                       &pos100);
            binmax = binmin + pos100 * bindelta;
            delta = width/pos100;
            numbins = pos100;
          }
      }
    if(valmax < binvalsmax)
      {
	binvalsmax = valmax;
      }
    if(nofreq)
      {
        //	if(p99 > 5000.)
        //  p99 = 5000.;
        if(binvalsmax < p99)
            p99 = binvalsmax;
        
        pdf_create_axes(VERTICAL, 0, p99, 5, lx, ly, height);
      }
    else
      {
        pdf_create_axes(VERTICAL, 0, binvalsmax/total, 5, lx, ly, height);
      }
    pdf_create_axes(HORIZONTAL, binmin , binmax, 6, lx, ly, width);
    if(world->options->has_datefile)
      {
        ratio = width / (binmax - binmin);
        for(locus = 0; locus < world->loci; locus++)
          {
	    if(world->bayes->mu)
	      averagerate += world->options->meanmu[locus] * world->bayes->histogram[locus].modes[world->numpop2+world->locus];
	    else
	      averagerate += world->options->meanmu[locus];
	  }
	averagerate /= world->loci;
	pop = topop;
	//	for(pop=0; pop < world->numpop; pop++)
	//  {
	if(world->data->numind != NULL)
	  {
	    for(ind = 0; ind < world->data->numind[pop][0]; ind++)
	      {
		xcoord = world->data->sampledates[pop][0][ind].date 
		  * world->options->generation_year 
		  * averagerate;
		pdf_print_dot(lx+(xcoord * ratio), ly-3.0, 3.0, DIAMOND, blue);
	      }
	  }
	    //  }
      }

    
    //    x = delta/1.95;
    x = delta/2.0;

    xold = x;
    if(nofreq)
      {
	nbmin = (long) 2;//I need a better algorithm! MAX(2.,((double) numbins)/200.);
	bayes_smooth(lower,numbins,nbmin, TRUE);
	bayes_smooth(binvals,numbins, nbmin, TRUE);
	bayes_smooth(upper,numbins, nbmin, TRUE);
	vold = binvals[0];
      }
    for(i=0;i<numbins;i++)
      {
        if(set50[i] == '1')
            pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0.1, 0.1, 0.1));
        else
          {
            if(set95[i] == '1')
                pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0.4, 0.4, 0.4));
            else
                pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0.8, 0.8, 0.8));
          }
        if(nofreq)
          {
	    l = lower[i];
	    u = upper[i];
	    v = binvals[i];
	    if(confidence !=NULL && confidence[i]>=1.0)
	      {
		xold = x;
		vold = v;
		x += delta;
		continue;
	      }
	    pdf_contents_set_line_width(canvas, delta);
            if(u/p99 > 1.0)
              {
                if(l/p99 < (1.0 + SMALLEPSILON))
                  {
                    pdf_draw_line(lx+x,ly + ((float) l)/p99 * height,
                                  lx+x,ly + height);
                  }
                pdf_print_dot(lx+x, ly + height, 3, SQUARE, red);
              }
            else
              {
                pdf_draw_line(lx+x,ly + ((float) l)/p99 * height,
                              lx+x,ly + ((float) u)/p99 * height);
              }
            pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0.2, 0.2, 0.2));
            
            if(((float) v)/p99 < (1.0 + SMALLEPSILON))
              {
		//                pdf_draw_line(lx+x,ly + ((float) v)/p99 * height-1.0F,
                //              lx+x,ly + ((float) v)/p99 * height+1.0F);
		pdf_contents_set_line_width(canvas, 1.0);
                pdf_draw_line(lx+xold,ly + ((float) vold)/p99 * height,
                              lx+x,ly + ((float) v)/p99 * height);
		vold  = v;
		xold = x;
              }
	    else
	      {
		if(vold<p99)
		  {
		    pdf_contents_set_line_width(canvas, 1.0);
		    pdf_draw_line(lx+xold,ly + ((float) vold)/p99 * height,
				  lx+x,ly + height);
		  }
		vold  = p99;
		xold = x;		
	      }
          }
        else
          {
            pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0.8, 0.8, 0.8));
            pdf_contents_set_line_width(canvas, delta);
            pdf_draw_line(lx+x,ly+0,lx+x,ly + ((float) binvals[i])/binvalsmax * height);
            pdf_contents_set_line_width(canvas, delta/3);
            pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0.5, 0.5, 0.5));
            pdf_draw_line(lx+x,ly+((float) binvals[i])/binvalsmax * height,
                          lx+x,ly + ((float) upper[i])/binvalsmax * height);
            pdf_draw_line(lx+x,ly+((float) lower[i])/binvalsmax * height,
                          lx+x,ly + ((float) binvals[i])/binvalsmax * height);
          }
        //fprintf(stderr,"%f %f %f\n",lx+x,ly+ binvals[i]/binvalsmax * height, binvals[i]);
        if(confidence != NULL)
          {
            r = confidence[i];
	    if(r>1.0 && r < 0.0)
	      continue;
            pdf_contents_set_line_width(canvas, delta);
            pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(r,r,r));
            pdf_draw_line(lx+x,ly + height+5.0F,
                          lx+x,ly + height+10.0F);
            //	    pdf_print_dot(lx+x, ly + height + 5, 3, SQUARE, red);
          }
        x += delta;
      }
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    pdf_contents_set_line_width(canvas, 1);
    myfree(lower);
    myfree(upper);
}

void pdf_master_init(world_fmt *world, option_fmt *options, data_fmt *data)
{
        ////////////////////////////////////  
        pdf_init();
        // add first page
        pdf_new_page(options->title);
        pdf_master_title(options->title,  &page_height, &left_margin);
        pdf_print_options(world, options, data, &page_height, &left_margin);
	if(options->datatype!='g')
	  {
	    pdf_print_data_summary(world, options, data,  &page_height, &left_margin);
	    pdf_print_data (world, options, data, &page_height, &left_margin);
	  }
        ////////////////////////////////////  
}


int pdf_init()
{
    pdf_type1_fontdef font1_def;
    pdf_type1_fontdef font2_def;
    pdf_type1_fontdef font3_def;
    pdf_type1_fontdef font4_def;
    
    /* Create a new PDF document. */
    doc = pdf_doc_new();
    pdf_doc_new_doc(doc);
    
    /* Add Helvetica Font. */
    font1_def = pdf_create_type1_fontdef(PDF_FONT_HELVETICA);
    pdf_doc_add_type1font(doc, font1_def, NULL, NULL);
    /* Add Helvetica-Oblique Font. */
    font2_def = pdf_create_type1_fontdef(PDF_FONT_HELVETICA_OBLIQUE);
    pdf_doc_add_type1font(doc, font2_def, NULL, NULL);
    /* Add Symbol Font. */
    font3_def = pdf_create_type1_fontdef(PDF_FONT_SYMBOL);
    pdf_doc_add_type1font(doc, font3_def, NULL, NULL);
    /* Add Courier Font. */
    font4_def = pdf_create_type1_fontdef(PDF_FONT_COURIRE);
    pdf_doc_add_type1font(doc, font4_def, NULL, NULL);
    return 0;
}

///
/// generate a new PDF page with a border rectangle and with page numbers and
/// title in top right corner and impressum at bottom left
int pdf_new_page(char *title)
{
    char stemp[LINESIZE];
    if(strlen(title)>1)
        strncpy(pdf_pagetitle,title,80);
    /* Add a page to the document. */
    page = pdf_doc_add_page(doc);
    page_counter += 1;
    /* Get the canvas object of the page. */
    canvas = pdf_page_get_canvas(page);
    
    /*Set current font to "Helvetica" and set the size to 10. */
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    
    pdf_contents_set_line_width(canvas, 1);
    /* draw page rectangle 50pt from border*/
    pdf_contents_rectangle(canvas, 50, 50, 
                           pdf_contents_get_width(canvas) - 100, 
                           pdf_contents_get_height(canvas) - 110);
    pdf_contents_stroke(canvas);
    /* print the title of the analysis*/
    pdf_print_header(pdf_pagetitle);
    /* print the impressum at the bottome*/
    sprintf(stemp,"Migrate %s: (http://popgen.sc.fsu.edu) [program run on %s]",MIGRATEVERSION, pdf_time);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 6);
    pdf_print_contents_at(50, 42, stemp);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    return 0;
}

///
/// start a new page and print the section title on top of page
void pdf_print_section_title(float *page_width, float * page_height, char *title)
{
    float w;
    // setup new page and title
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 18);    
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    *page_height = pdf_contents_get_height(canvas) - 75.;
    *page_width = pdf_contents_get_width(canvas);
    *page_height -= 20.;
    pdf_print_contents_at((*page_width - w)/2., *page_height, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    *page_height -= 20.;
    pdf_draw_line(50, *page_height, *page_width-50., *page_height);
    *page_height -= 20.;
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);    
}

///
/// print title on every page with page number
int pdf_print_header(char *title)
{
    float w;
    float page_height;
    float page_width;
    char *fulltitle;
    
    fulltitle = (char*) mycalloc(255,sizeof(char));
    /* Print the title of the page (with positioning center). */
    sprintf(fulltitle,"%s -- %i",title, page_counter);
    //printf("%s\n",fulltitle);
    w = (float) pdf_contents_get_text_width(canvas, fulltitle, NULL, NULL);
    /* Start to print text. */
    pdf_contents_begin_text(canvas);    
    /* Move the position of the text to center */
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);
    pdf_contents_move_text_pos(canvas, (page_width - w - 50),
                               page_height - 50);
    
    /* Print title with pagenumber to rightadjusted */
    pdf_contents_show_text(canvas, fulltitle);    
    /* Finish to print text. */
    pdf_contents_end_text(canvas);
    
    myfree(fulltitle);
    
    return 0;
}

void pdf_print_contents_at(float x, float y, char *title)
{
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, x, y);
    pdf_contents_show_text(canvas, title);
    pdf_contents_end_text(canvas);
    //fprintf(stderr,"(%f, %f) %s\n",x,y,title);
}

///
/// print first title page
int pdf_master_title(char *title, float *orig_page_height, float *left_margin)
{
#ifdef MPI
    int cpus;
#endif
    float w;
    float page_height;
    float page_width;
    char newtitle[LINESIZE];
    
    get_time (pdf_time, "%H:%M:%S");

    /* Move the position of the text to center */
    
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);
    /* Print the title of the page (with positioning center). */
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 24);
    
    MYSNPRINTF(newtitle,(long) (page_width-110)/10, "%s",title);
    w = (float) pdf_contents_get_text_width(canvas, newtitle, NULL, NULL);
    *left_margin = 55;
    /* Start to print text. */
    pdf_contents_begin_text(canvas);    
    /* Print title centered */
    page_height -= 100;
    pdf_print_contents_at((page_width - w)/2, page_height, newtitle);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    page_height -= 26;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 12);
    page_height -= 24;
    pdf_print_contents_at(55, page_height, "MIGRATION RATE AND POPULATION SIZE ESTIMATION");  
    pdf_advance(&page_height);
    pdf_print_contents_at(55, page_height, "using the coalescent and maximum likelihood or Bayesian inference");    
    pdf_advance(&page_height);
    pdf_printf(55, page_height, 'L', "Migrate-n version %s [%s]",MIGRATEVERSION, MIGRATESUBVERSION);
#ifdef BEAGLE
    pdf_advance(&page_height);
    pdf_printf(55, page_height, 'L',"  Compiled for GPU or other acceleration methods\n");
#endif    
#ifdef MPI
    pdf_advance(&page_height);
    pdf_printf(55, page_height, 'L', "  Compiled for a PARALLEL COMPUTER ARCHITECTURE\n");
    pdf_advance(&page_height);
    cpus = numcpu - 1;
    pdf_printf(55, page_height, 'L', "  One master and %i compute nodes are available.\n", cpus);
#endif
#ifdef THREAD
    pdf_advance(&page_height);
    pdf_printf(55, page_height, 'L',"  Compiled for a SYMMETRIC MULTIPROCESSORS\n");
#endif
#ifdef SNOWLEOPARD
    pdf_advance(&page_height);
    pdf_printf(55, page_height, 'L',"  Compiled for a SYMMETRIC MULTIPROCESSORS\n");
#endif
#ifdef ALTIVEC
    pdf_advance(&page_height);
    pdf_printf(55, page_height,  'L',"  ALTIVEC enabled\n");
#endif
#ifdef FAST_EXP
    pdf_advance(&page_height);
    pdf_printf(55, page_height,  'L',"  Fast approximation to Exp() and Log() used\n");
#endif
    pdf_advance(&page_height);
    pdf_print_time(55, &page_height, "Program started at ");
    firstcanvas = canvas;
    timestampy = page_height;
    // here the end time stamp will be printed at the very end the run
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    // print migrate logo
    pdf_contents_set_line_width(canvas, 1);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    pdf_contents_set_rgb_fill(canvas, PdfRGBColor(0.99, 0, 0));
    pdf_migrate_logo(-1410, -1435 - LETTERADJUST, 100.0);
    pdf_contents_set_rgb_fill(canvas, PdfRGBColor(0, 0.1, 0.9));
    pdf_migrate_logo(-1374, -1435  - LETTERADJUST, 100.0);
    pdf_contents_set_line_width(canvas, 4);
    pdf_migrate_logo_lines(-1410, -1435 - LETTERADJUST, 100.0);
    *orig_page_height = page_height - 3*LINESTRETCH;
    return 0;
}

int pdf_write_file(option_fmt *options)
{
    pdf_doc_write_to_file(doc, options->pdfoutfilename);
    return 0;
}

///
/// print elements of Bayesian table for character variables
float pdf_print_line_element(float lx, float ly, float offset, char *title)
{
    float w=0;
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    if(offset>=0)
        pdf_print_contents_at(lx+offset-w, ly, title);
    else
        pdf_print_contents_at(lx-offset, ly, title);
    return w;
}

///
/// print elements of Bayesian table for float variables
float pdf_print_line_element2(float lx, float ly, float offset, float value, int fmt1, int fmt2)
{
    float w=0;
    char title[100];
    sprintf(title,"%*.*f",fmt1,fmt2,value);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    if(offset>=0)
        pdf_print_contents_at(lx+offset-w, ly, title);
    else
        pdf_print_contents_at(lx-offset, ly, title);
    return w;
}


void pdf_print_line_theta(float lx, float ly, float offset, long j)
{
    char tempstring[100];
    char * thetatitle="Q";
    if(j < 0)
        sprintf(tempstring," ");
    else
        sprintf(tempstring,"%li",j+1);
    pdf_contents_set_font_and_size(canvas, "Symbol", 11);
    pdf_print_contents_at(lx-offset,ly,thetatitle);
    pdf_contents_set_font_and_size(canvas, "Symbol", 8);
    pdf_print_contents_at(lx-offset+10,ly-4,tempstring);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
}

/// print rate line
void pdf_print_line_rate(float lx, float ly, float offset, long j, long exponent)
{
    char tempstring[100];
    char * rtitle="m";
    if(j < 0)
        sprintf(tempstring," ");
    else
        sprintf(tempstring,"%li",j+1);
    pdf_contents_set_font_and_size(canvas, "Symbol", 10);
    pdf_print_contents_at(lx-offset,ly,rtitle); // print mu
    pdf_contents_set_font_and_size(canvas, "Symbol", 8);
    pdf_print_contents_at(lx-offset+10,ly-4,tempstring);//print locus number as subscript

    pdf_contents_set_font_and_size(canvas, "Symbol", 10);
    sprintf(tempstring,"[10");
    pdf_print_contents_at(lx-offset+13,ly,tempstring); //print the scale should look like this [x10-5]
    pdf_contents_set_font_and_size(canvas, "Symbol", 8);
    sprintf(tempstring,"%li",exponent); 
    pdf_print_contents_at(lx-offset+26,ly+4,tempstring); // print the exponent as superscript
    pdf_contents_set_font_and_size(canvas, "Symbol", 10);
    sprintf(tempstring,"]");
    pdf_print_contents_at(lx-offset+35,ly,tempstring);

    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
}

///
/// pretty print the M values
void  pdf_print_line_mig(char *migtitle, float lx, float page_height, float offset, long frompop, long topop)
{
    float w;
    char tostring[100];
    char tempstring[100];
    if(topop < 0)
        sprintf(tostring,"+");
    else
        sprintf(tostring,"%li",topop+1);
    sprintf(tempstring,"%li->%s",frompop+1,tostring);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    pdf_print_contents_at(lx-offset, page_height,migtitle);
    w = (float) pdf_contents_get_text_width(canvas, migtitle, NULL, NULL);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 8);
    pdf_print_contents_at(lx-offset+w,page_height-4,tempstring);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
}


///
/// Print Bayesian table header
/// and report the width of the text in the headerline
void pdf_print_bayestableheader(float page_height, float left_margin, float right_margin, float *offset)
{
    float lx = left_margin;
    float ly = page_height;
    
    pdf_print_line_element(lx, ly, offset[0], "Locus");
    pdf_print_line_element(lx, ly, offset[1], "Parameter");
    pdf_print_line_element(lx, ly, offset[2], "2.5%");
    pdf_print_line_element(lx, ly, offset[3], "25.0%");
    pdf_print_line_element(lx, ly, offset[4], "Mode");
    pdf_print_line_element(lx, ly, offset[5], "75.0%");
    pdf_print_line_element(lx, ly, offset[6], "97.5%");
    pdf_print_line_element(lx, ly, offset[7], "Median");
    pdf_print_line_element(lx, ly, offset[8], "Mean");
}

// Parameter        2.5%%      25.0%%   median    75.0%%   97.5%%     mode     mean\n"

///
/// print Bayesian table
void pdf_print_bayestable(world_fmt *world)
{
    int fmt1;
    int fmt2;
    double mu;
    long lmu;
    double meanmu;
    long j0, j;
    long l;
    long size = world->numpop2 + (world->bayes->mu);
    float lx;
    float *offset;
    bayes_fmt * bayes = world->bayes;
    bayeshistogram_fmt * hist;
    
    float page_height = pdf_contents_get_height(canvas);
    float page_width = pdf_contents_get_width(canvas);
    
    float left_margin = 60;
    float right_margin = page_width - 60;
    
    //    float w;
    
    char st[100];
    char title[100] = "Bayesian Analysis: Posterior distribution table";
    long frompop;
    long topop;
    long locus;
    long lozi = world->loci > 1 ? world->loci : 0;
    
    //column to to right-align the table columns
    offset = (float *) mycalloc(9,sizeof(float));
    offset[0] = -1; //left align
    offset[1] = -40;
    offset[2] = 135;
    offset[3] = 191;
    offset[4] = 247;
    offset[5] = 303;
    offset[6] = 359;
    offset[7] = 415;
    offset[8] = 471;
    pdf_print_section_title(&page_width, &page_height, title);
    //page_height -= 126;
    lx = left_margin;
    page_height -= 20;
    pdf_draw_line(left_margin, page_height, right_margin, page_height);
    pdf_advance(&page_height);
    pdf_print_bayestableheader(page_height, left_margin, right_margin, offset);
    pdf_advance(&page_height);
    pdf_draw_line(left_margin, page_height, right_margin, page_height);
    pdf_advance(&page_height);    
    for(locus=0; locus <= lozi; locus++)
      {
	if(!world->data->skiploci[locus])
	  {
	    mu = 1.0; //used to adjust values for rate estimates for others this is 1.
	    hist = &bayes->histogram[locus];
	    if(locus == world->loci)
	      strcpy(st,"  All ");
	    else
	      sprintf(st,"%5li ",locus + 1);
	    
	    for(j0=0; j0< size; j0++)
	      {
		if(j0 < world->numpop2)
		  {
		    //if(strchr("0c", world->options->custm2[j]))
		    //  continue;
		    if(world->bayes->map[j0][1] == INVALID)
		      continue;
		  }
		j = world->bayes->map[j0][1];
		pdf_print_line_element(lx, page_height, offset[0], st);
		if(j < world->numpop)
		  {
		    pdf_print_line_theta(lx, page_height, offset[1], j0);
		    fmt1 = 8;
		    fmt2 = 5;
		  }
		else
		  {
		    if(j >= world->numpop2)
		      {
			// rate modifier used
			if(locus==lozi && lozi>1)
			  {
			    meanmu = 0.;
			    for(l=0;l<world->loci;l++)
			      {
				meanmu += world->options->meanmu[l];
			      }
			    meanmu /= world->loci;
			    lmu = (long) floor( log10(hist->modes[j] * meanmu));
			    mu = meanmu * pow(10. , -lmu);
			  }
			else
			  {
			    lmu = (long) floor( log10(hist->modes[j]*world->options->meanmu[locus]));
			    mu = world->options->meanmu[locus] * pow(10. , -lmu);
			  }
			pdf_print_line_rate(lx, page_height, offset[1], -1, lmu);
			fmt1 = 8;
			fmt2 = 5;
		      }
		    else
		      {
			m2mm(j0,world->numpop,&frompop, &topop);
			if(world->options->usem)
			  {
			    pdf_print_line_mig("M", lx, page_height, offset[1], frompop, topop);
			    fmt1 = 8;
			    if (strchr (SEQUENCETYPES, world->options->datatype))
			      fmt2 = 1;
			    else
			      fmt2 = 3;
			  }
			else
			  {
			    pdf_print_line_mig("xNm", lx, page_height, offset[1], frompop, topop);
			    fmt1 = 8;
			    fmt2 = 5;
			  }
		      }
		  }
		// reset the to helvetica from symbol or small helvetica
		pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
		pdf_print_line_element2(lx, page_height, offset[2], mu * hist->cred95l[j],fmt1,fmt2);
		pdf_print_line_element2(lx, page_height, offset[3], mu * hist->cred50l[j],fmt1,fmt2);
		pdf_print_line_element2(lx, page_height, offset[4], mu * hist->modes[j],fmt1,fmt2);
		pdf_print_line_element2(lx, page_height, offset[5], mu * hist->cred50u[j],fmt1,fmt2);
		pdf_print_line_element2(lx, page_height, offset[6], mu * hist->cred95u[j],fmt1,fmt2);
		pdf_print_line_element2(lx, page_height, offset[7], mu * hist->medians[j],fmt1,fmt2);
		pdf_print_line_element2(lx, page_height, offset[8], mu * hist->means[j],fmt1,fmt2);
		page_height -= LINESTRETCH;
		if(page_height < 55)
		  {
		    pdf_new_page("");
		    page_height = pdf_contents_get_height(canvas);
		    page_width = pdf_contents_get_width(canvas);
		    page_height -= 50;
		    lx = left_margin;
		    page_height -= 20;
		    pdf_draw_line(left_margin, page_height, right_margin, page_height);
		    page_height -= 20;
		    pdf_print_bayestableheader(page_height, left_margin, right_margin, offset);
		    page_height -= 20;
		    pdf_draw_line(left_margin, page_height, right_margin, page_height);
		    pdf_advance(&page_height);
		  }
	      }
	    pdf_draw_line(left_margin, page_height, right_margin, page_height);
	    pdf_advance(&page_height);        
	  } // only when loci contains info
      } // over all loci
    myfree(offset);
}


///
/// print out the acceptance ratios for all the different Bayesian updates
void
pdf_bayes_print_accept(world_fmt *world)
{
    char title[LINESIZE];
    float w;
    float left_margin = 55;
    float page_width;
    long j0,j=0;             //used to loop over all parameters
    long topop    =0;   // indicator into the parameter vector, specifying originating population 
    long frompop  =0;   // receiving population
    char *stempo;       // string variable holding print-string 
    char *stemp;        // pointer to string, seems to be need to don't get MYREAL free warnings
    long trials   =0;   //
    long tc = world->numpop2 + (world->bayes->mu);// * world->loci); //position of genealogy accept rates
    bayes_fmt *bayes = world->bayes; 

    stempo = (char *) mycalloc(LINESIZE,sizeof(char));
    stemp = stempo;
    
    sprintf(title,"Acceptance ratios for all parameters and the genealogies");
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    
    /* Start to print text. */ 
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    
// This needs more attention but will need more stuff to safe
    if(world->options->datatype == 'g')
      {
	    pdf_print_contents_at(left_margin, page_height, "not available with datatype=Genealogy");
	    pdf_advance(&page_height);
	    myfree(stempo);
	    return;
      }

    pdf_print_contents_at(left_margin, page_height, "Parameter");
    pdf_print_contents_at(250, page_height, "Accepted changes");
    pdf_print_contents_at(450, page_height, "Ratio");   
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
    // population sizes
    for(j0=0; j0 < world->numpop; j0++)
      {
        if(!strchr("c", bayes->custm2[j0]))
          {
            j = world->bayes->map[j0][1];
            if((trials=bayes->trials[j])>0)
              {
                symbol_Theta(left_margin, page_height, 12, j0+1);
                pdf_printf(250, page_height, 'L', "%8li/%-8li",bayes->accept[j],trials);
                pdf_printf(450, page_height, 'L', "%8.5f", (MYREAL) bayes->accept[j]/trials); 
                pdf_advance(&page_height);
              }
          }
      }
    // migration rates
    for(j0=world->numpop; j0 < world->numpop2; j0++)
      {
        if(!strchr("0c", bayes->custm2[j0]))
          {
            j = world->bayes->map[j0][1];
            if((trials=bayes->trials[j])>0)
              {
                m2mm (j0, world->numpop, &frompop, &topop);
                symbol_M(left_margin, page_height, 12, frompop+1, topop+1, world->options->usem);
                pdf_printf(250, page_height, 'L', "%8li/%-8li",bayes->accept[j],trials);
                pdf_printf(450, page_height, 'L', "%8.5f", (MYREAL) bayes->accept[j]/trials); 
                pdf_advance(&page_height);
                memset(stemp,0,sizeof(char)*(LINESIZE-1));        
              }
          }
      }
    // accepted rate of mutation rate changes
    if(bayes->mu)
      {
        if((trials=bayes->trials[j0])>0)
          {
	    for(j=world->numpop2; j < world->numpop2 + world->loci;j++)
	      {
		symbol_R(left_margin, page_height, 12, j+1-world->numpop2);
		pdf_printf(250, page_height, 'L', "%8li/%-8li",bayes->accept[j0],trials);
		pdf_printf(450, page_height, 'L', "%8.5f", (MYREAL) bayes->accept[j0]/trials); 
		pdf_advance(&page_height);
	      }      
	  }
      }
    // accepted trees
    if((trials=bayes->trials[tc])>0)
      {
        pdf_printf(left_margin, page_height,'L', "Genealogies");
        pdf_printf(250, page_height, 'L', "%8li/%-8li",bayes->accept[tc], trials);
        pdf_printf(450, page_height, 'L', "%8.5f", (MYREAL) bayes->accept[tc]/ trials); 
        pdf_advance(&page_height);
      }
    myfree(stempo);
}

///
/// print out the autocorrelation in replicates and the effective sample size
void
pdf_bayes_print_ess(world_fmt *world)
{
    char title[LINESIZE];
    float w;
    float left_margin = 55;
    float page_width;
    long j0,j=0;             //used to loop over all parameters
    long topop    =0;   // indicator into the parameter vector, specifying originating population 
    long frompop  =0;   // receiving population
    char *stempo;       // string variable holding print-string 
    char *stemp;        // pointer to string, seems to be need to don't get MYREAL free warnings
    bayes_fmt *bayes = world->bayes; 
    stempo = (char *) mycalloc(LINESIZE,sizeof(char));
    stemp = stempo;
    
    sprintf(title,"MCMC-Autocorrelation and Effective MCMC Sample Size");
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    
    /* Start to print text. */ 
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    
    pdf_print_contents_at(left_margin, page_height, "Parameter");
    pdf_print_contents_at(250, page_height, "Autocorrelation");
    pdf_print_contents_at(450, page_height, "Effective Sampe Size");   
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
    if(world->options->bayes_infer)
      {
	// population sizes
	for(j0=0; j0 < world->numpop; j0++)
	  {
	    if(!strchr("c", world->options->custm2[j0]))
	      {
		j = world->bayes->map[j0][1];
		symbol_Theta(left_margin, page_height, 12, j0+1);
		pdf_printf(250, page_height, 'L', "%5.5f",world->auto_archive[j]);
		pdf_printf(450, page_height, 'L', "%10.2f", world->ess_archive[j]); 
		pdf_advance(&page_height);
	      }
	  }
	// migration rates
	for(j0=world->numpop; j0 < world->numpop2; j0++)
	  {
	    if(!strchr("0c", bayes->custm2[j0]))
	      {
		j = world->bayes->map[j0][1];
		m2mm (j0, world->numpop, &frompop, &topop);
		symbol_M(left_margin, page_height, 12, frompop+1, topop+1, world->options->usem);
		pdf_printf(250, page_height, 'L', "%5.5f",world->auto_archive[j]);
		pdf_printf(450, page_height, 'L', "%10.2f", world->ess_archive[j]); 
		pdf_advance(&page_height);
		memset(stemp,0,sizeof(char)*(LINESIZE-1));        
	      }
	  }
      }
    j = world->numpop2; //adjusting j to fit  
    // accepted rate of mutation rate changes
    if(world->bayes->mu)
      {
	for(j=world->numpop2; j < world->numpop2 + world->loci;j++)
	  {
	    symbol_R(left_margin, page_height, 12, j+1-world->numpop2);
	    pdf_printf(250, page_height, 'L', "%5.5f",world->auto_archive[j]);
	    pdf_printf(450, page_height, 'L', "%10.2f", world->ess_archive[j]); 
	    pdf_advance(&page_height);
	  }      
      }
    // likelihood of trees
    pdf_print_contents_at(left_margin, page_height,"Ln[Prob(D|G)]");
    pdf_printf(250, page_height, 'L', "%5.5f",world->auto_archive[j]);
    pdf_printf(450, page_height, 'L', "%10.2f", world->ess_archive[j]); 
	    pdf_advance(&page_height);

    myfree(stempo);
}


///
/// print out marginal likelihood calculated from thermodynamic integration, harmonic and arithmetic mean.
void pdf_bayes_factor_header(world_fmt *world, option_fmt *options)
{
    char title[LINESIZE];
    float w;
    float left_margin = 55;
    float page_width;
    sprintf(title,"Log-Probability of the data given the model (marginal likelihood)");
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    
    /* Start to print text. */ 
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height, 'L', "Use this value for Bayes factor calculations:");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L', "BF = Exp[ ln(Prob(D | thisModel) - ln( Prob( D | otherModel)");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L', "or as LBF = 2 (ln(Prob(D | thisModel) - ln( Prob( D | otherModel))");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L', "shows the support for thisModel]");
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
}


void pdf_bayes_factor_rawscores_header(world_fmt *world, option_fmt *options)
{
    float left_margin = 55;
    float page_width;
    page_width = pdf_contents_get_width(canvas);
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height, "Locus");
    pdf_print_contents_at(100, page_height, "Raw thermodynamic score(1a)");
    pdf_print_contents_at(250, page_height, "Bezier approximation score(1b)");
    pdf_print_contents_at(460, page_height, "Harmonic mean(2)");
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
}

void pdf_bayes_factor_rawscores(long locus, MYREAL rawtermo, MYREAL beziertermo, MYREAL harmo)
{
  float page_width;
  page_width = pdf_contents_get_width(canvas);
  if(locus<0)
    {
      pdf_draw_line(50, page_height,  page_width-50, page_height);
      pdf_advance(&page_height);
      pdf_printf_right(left_margin+20, page_height,"All  ");
    }
  else
    pdf_printf_right(left_margin+20, page_height,"%5li", locus+1);
  pdf_printf_right(200, page_height,"   %20.2f", rawtermo);
  pdf_printf_right(360, page_height,"   %20.2f", beziertermo);
  pdf_printf_right(520, page_height,"   %20.2f", harmo);
  pdf_advance(&page_height);
}

void pdf_bayes_factor_rawscores_harmo(long locus, MYREAL harmo)
{
  float page_width;
  page_width = pdf_contents_get_width(canvas);
  if(locus < 0)
    {
      pdf_draw_line(50, page_height, page_width-50, page_height);
      pdf_advance(&page_height);
      pdf_printf_right(left_margin+20, page_height,"All  ");
    }
  else
    pdf_printf(left_margin, page_height,'R',"%5li", locus+1);
  pdf_draw_line(140, page_height, 150, page_height);
  pdf_draw_line(280, page_height, 290, page_height);
  pdf_printf_right(520, page_height,"   %20.2f", harmo);
  pdf_advance(&page_height);
}

///
/// print out marginal likelihood calculated from thermodynamic integration, harmonic and arithmetic mean.
void
pdf_bayes_factor(world_fmt *world,  MYREAL tsum, MYREAL tsum2, MYREAL hsum, MYREAL asum, long maxreplicate, MYREAL scaling_factor)
{
    float left_margin = 55;
    float page_width;
    page_width = pdf_contents_get_width(canvas);
    pdf_advance(&page_height);
    if(world->loci<=1)
      {
	pdf_draw_line(50, page_height, page_width-50, page_height);
	pdf_advance(&page_height);
	pdf_print_contents_at(left_margin, page_height, "Method");
	pdf_print_contents_at(250, page_height, "ln(Prob(D|Model))");
	pdf_print_contents_at(450, page_height, "Notes");   
	pdf_advance(&page_height);
	pdf_advance(&page_height);
	pdf_draw_line(50, page_height, page_width-50, page_height);
	pdf_advance(&page_height);
	pdf_print_contents_at(left_margin, page_height, "Thermodynamic integration");
	if(world->options->heating)
	  {
	    pdf_printf(250, page_height,'L',"%f", tsum);
	    pdf_printf(480, page_height,'L',"%s", "(1a)");
	    pdf_print_contents_at(left_margin, page_height, "");
	    pdf_advance(&page_height);
	    pdf_printf(250, page_height,'L',"%f", tsum2);
	    pdf_printf(480, page_height,'L',"%s", "(1b)");
	  }
	else
	  {
	    pdf_printf(250, page_height,'L',"(not estimated [no heating])");
	    pdf_printf(480, page_height,'L',"%s", "(1)");
	  }
	pdf_advance(&page_height);
	pdf_print_contents_at(left_margin, page_height, "Harmonic mean");
	pdf_printf(250, page_height,'L',"%f", hsum);
	pdf_printf(480, page_height,'L',"%s", "(2)");
	pdf_advance(&page_height);
	//pdf_print_contents_at(left_margin, page_height, "Arithmetic mean");
	//pdf_printf(250, page_height,'L',"%f", asum);
	//pdf_printf(480, page_height,'L',"%s", "(3)");
	//pdf_advance(&page_height);
	pdf_draw_line(50, page_height, page_width-50, page_height);
	pdf_advance(&page_height);
      }
    pdf_printf(left_margin, page_height,'L',"(1a, 1b and 2) is an approximation to the marginal likelihood, make sure the program run long enough!");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L',"(1a, 1b) and (2) should give a similar result, (2) is considered more");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L',"crude than (1), but (1) needs heating with several well-spaced chains,");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L',"(1b) is using a Bezier-curve to get better approximations for runs with low number");
    pdf_advance(&page_height);
    pdf_printf(left_margin, page_height,'L',"of heated chains");
    pdf_advance(&page_height);
    switch(world->options->adaptiveheat)
      {
      case STANDARD:
        pdf_printf(left_margin, page_height,'L',"%s","Adaptive heating was ON, therefore the values of (1) may be incorrect),");
        pdf_advance(&page_height);
	break;
      case BOUNDED:
        pdf_printf(left_margin, page_height,'L',"%s","Adaptive heating with bounds was ON, therefore the values of (1) may be incorrect),");
        pdf_advance(&page_height);
	break;
      }
    if(world->loci>1)
      {
	pdf_printf(left_margin, page_height,'L',"[Scaling factor = %f", scaling_factor);
      }
    pdf_advance(&page_height);
    pdf_advance(&page_height);
}

void pdf_burnin_stops(world_fmt *world, long maxreplicate)
{
    long z;
    float w;
    char title[LINESIZE];
    float left_margin = 55;
    float page_width;
    sprintf(title,"Stop of burnin-in phase due to convergence");
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    
    /* Start to print text. */ 
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height, "Locus");
    pdf_print_contents_at(120, page_height, "Replicate");
    pdf_print_contents_at(210, page_height, "Steps");   
    pdf_print_contents_at(310, page_height, "Variance ratio (new/old variance)");   
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_advance(&page_height);
    for(z=0; z < world->loci * maxreplicate; z++)
      {
        pdf_printf(left_margin, page_height,'L',"%5li", world->burnin_stops[z].locus);
        pdf_printf(120, page_height,'L',"%10li",world->burnin_stops[z].replicate);
        pdf_printf(200, page_height,'L',"%10li", world->burnin_stops[z].stopstep);
        pdf_printf(310, page_height,'L',"%f", world->burnin_stops[z].variance/world->burnin_stops[z].oldvariance);
        pdf_printf(390, page_height,'L',"(%f/%f)",world->burnin_stops[z].variance,
                   world->burnin_stops[z].oldvariance);
        pdf_advance(&page_height);
      }
    pdf_advance(&page_height);
    pdf_advance(&page_height);
}

void pdf_print_stored_warnings(world_fmt *world)
{
  float w;
  char title[LINESIZE];
  float left_margin = 55;
  float page_width;
  char *buffer;
  char *b;
  char *tmp;
  if(world->warningsize>0)
    {
      sprintf(title,"Potential Problems");
      pdf_new_page("");
      pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
      w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
      
      /* Start to print text. */ 
      page_height = pdf_contents_get_height(canvas);
      page_width = pdf_contents_get_width(canvas);
      pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
      pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
      page_height -= 126;
      pdf_draw_line(50, page_height, page_width-50, page_height);
      pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
      
      pdf_advance(&page_height);
      pdf_advance(&page_height);
      pdf_printf_next(left_margin, &page_height,"This section reports potential problems with your run. Such reporting");      
      pdf_printf_next(left_margin, &page_height,"is tricky. When many parameters in  multilocus analysis");      
      pdf_printf_next(left_margin, &page_height,"are estimated then it is very common that some parameters for some loci");      
      pdf_printf_next(left_margin, &page_height,"will not be very informative, triggering suggestions (for example to");
      pdf_printf_next(left_margin, &page_height,"increase the prior range) that are not sensible. This suggestion tool");
      pdf_printf_next(left_margin, &page_height,"will improve with time. If some parameters are flagged inspect the tables");
      pdf_printf_next(left_margin, &page_height,"carefully and judge wether an action is required. For example if you run");
      pdf_printf_next(left_margin, &page_height,"a Bayesian inference with sequence data, for macroscopic species there is");
      pdf_printf_next(left_margin, &page_height,"rarely the need to increase the prior for Theta beyond 0.1, but if you use");
      pdf_printf_next(left_margin, &page_height,"microsatellites it is rather common that your prior for Theta has a range");
      pdf_printf_next(left_margin, &page_height,"from 0.0 to 100 or more. With many populations (>3) it is also very common");
      pdf_printf_next(left_margin, &page_height,"that some migration routes are estimated poorly because there is no data.");
      pdf_printf_next(left_margin, &page_height,"Increasing the range will not help in such situations, reducing number of");
      pdf_printf_next(left_margin, &page_height,"parameters may help in such situations.");     
      pdf_advance(&page_height);
      pdf_advance(&page_height);
      if(world->warning[0] == '\0')
	{
	  pdf_printf(left_margin, page_height,'L',"No warning was recorded during the run");
	}
      else
	{
	  buffer = (char *) mycalloc(strlen(world->warning)+1,sizeof(char));
	  sprintf(buffer,"%s",world->warning);
	  b = buffer;
	  tmp = strsep(&buffer,"\n");
	  while(tmp!=NULL)
	    {
	      pdf_printf(left_margin, page_height,'L',"%s",tmp);
	      pdf_advance(&page_height);
	      tmp = strsep(&buffer,"\n");
	    }
	  myfree(b);
	}
    }
}




/* DEBUG TODO ready to remove this function?*/
#if 0
///
/// plot migration time histograms, assumes that the ascii_printer has already filled the 
/// plotfield table when used with boolean PRECALC (=True)
void
pdf_mig_histogram(histogram_fmt ** histogram,
                  plotfield_fmt ** plotfield, long loci, long numparams,
                  long bins, long *sum, MYREAL ***migtable,  boolean precalc, world_fmt *world)
{
    float left_margin = 55;
    //float page_height;
    float page_width;
    float lx;
    float ly;
    float ph;
    
    long loc, i, j, z, zz;    
    long frompop;
    long topop;
    
    float biggest = 0.;
    float *binning;
    float *binvec;
    float tempmin = MYREAL_MAX;
    float tempmax = -MYREAL_MAX;
    float begin = MYREAL_MAX;
    float end = -MYREAL_MAX;
    float delta;
    
    MYREAL mtime;
    char *set50;
    char *set95;
    long weight;
    char title[LINESIZE];
    float w;
    long bin;
    MYREAL *temp;
    temp = (MYREAL *) mymalloc(sizeof(MYREAL) * bins);
    if (loci > 1)
        sprintf(title,"Migration event histogram over all loci");
    else
        sprintf(title,"Migration event histogram");
    
    // add a new page so that we can print at least four histograms
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    
    /* Start to print text. */ 
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    // first histogram position
    page_height -= 150;
    lx = 100;
    ly = page_height;
    
    set50 = (char *) mycalloc(bins+1, sizeof(char));
    set95 = (char *) mycalloc(bins+1, sizeof(char));
    memset (set50, '0', bins * sizeof(char));
    memset (set95, '0', bins * sizeof(char));
    
    binning = (float *) mycalloc (bins+1, sizeof (float));
    binvec = (float *) mycalloc (bins+1, sizeof (float));
    
    for (loc = 0; loc < loci; loc++)
      {
        for (i = (world->options->mighist_all ? 0 : world->numpop); i < numparams; i++)
          {
            if (migtable[loc][i][2] == NOAVERAGE)
              {
                plotfield[loc][i].print = FALSE;
                continue;
              }
            minmax (&histogram[loc][i], &tempmin, &tempmax);
            if (tempmin < begin)
                begin = tempmin;
            if (tempmax > end)
                end = tempmax;
          }
      }
    delta = (end - begin) / bins;
    binning[0] = begin + 0.5 * delta;
    for (i = 1; i < bins; i++)
        binning[i] = delta + binning[i - 1];
    if(!precalc)
      {
        for (loc = 0; loc < loci; loc++)
          {
            
            for (i = (world->options->mighist_all ? 0 : world->numpop); i < numparams; i++)
              {
                if (migtable[loc][i][2] == NOAVERAGE)
                    continue;
                memset (binvec, 0, sizeof (float) * (bins+1));
                for (j = 0; j < histogram[loc][i].count; j++)
                  {
                    mtime = histogram[loc][i].time[j];
                    weight = histogram[loc][i].weight[j];
                    z = 0;
                    while (mtime > binning[z] && z < bins)
                        z++;
                    binvec[z] += weight;
                  }
                biggest = 0.;
                for (j = 0; j < bins; j++)
                  {
                    plotfield[loc][i].y[j] = (long) binvec[j];
                    plotfield[loc][i].yfreq[j] = binvec[j] = binvec[j] / sum[loc];
                    if (biggest < binvec[j])
                        biggest = binvec[j];
                  }
                for (j = 0; j < bins; j++)
                  {
                    for (zz = 0;
                         zz <
                         (long) (binvec[j] * plotfield[loc][i].ysize / biggest);
                         zz++)
                        plotfield[loc][i].data[j][zz] = '+';
                    plotfield[loc][i].data[j][zz] = '\0';
                  }
              }
          }
      }
    // plot only overall loci
    // keep code that would allow to print everything
    for (loc = loci-1; loc < loci; loc++)
      {
        if (loc == (loci - 1))
          {
            
            //            if (loci > 1)
            //     pdf_printf_next(left_margin, &page_height,"Over all loci");
            //else
            //   pdf_printf_next(left_margin, &page_height,"Locus %li\n",loc + 1);
          }
        else
            pdf_printf_next(left_margin, &page_height,"Locus %li\n",loc + 1);
        
        for (i0 = (world->options->mighist_all ? 0 : world->numpop); i0 < numparams; i0++)
          {
            if (plotfield[loc][i0].print)
              {
                i = world->bayes->map[i0][1];
                m2mm(i0, world->numpop, &frompop, &topop);
                if(frompop==topop)
                  {
                    pdf_print_contents_at(lx-30,ly+125,"Freq. for ");
                    symbol_Theta(lx+12, ly+125,12,frompop+1);
                    pdf_print_contents_at(lx+90,ly-25, "Time [scaled by mutation rate / site / generation]");
                  }
                else
                  {
                    pdf_print_contents_at(lx-30,ly+125,"Freq. for ");
                    symbol_M(lx+12, ly+125, 12, frompop+1, topop+1, world->options->usem);
                    pdf_print_contents_at(lx+90,ly-25, "Time [scaled by mutation rate / site / generation]");
                  }
                pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
                // TODO convert all results and helper-functions to float
                // as we certainly do not print out higher precision than that
                // currentl only plotfield is compliant but pdf_histogram is written
                // the MYREAL that in oter parts of the program needs to be double
                // 
                for(bin=0;bin< bins; bin++)
                    temp[bin] = (MYREAL) plotfield[loc][i].yfreq[bin];
                //
                pdf_histogram(temp, set50, set95, bins, delta, 
                              0., -999, lx, ly, page_width - 55 - lx, 116, FALSE);
                page_height -= 160;
                if(i < (numparams-1))
                  {
                    pdf_page_advance_or_not(&page_height, 50);
                    ph = pdf_contents_get_height(canvas) - 55 - LINESTRETCH;
                    if(page_height >= ph)
                        page_height -= 150;
                    ly = page_height;
                  }
              }
          }
      }
    myfree(binning);
    myfree(binvec);
    myfree(set50);
    myfree(set95);
    myfree(temp);
}

///
/// print the table header 
void pdf_print_mighist_table_header(float lx, const float *offset, const float *roffset, 
                                    float *page_height, float left_margin, float right_margin)
{
    float *ly = page_height;
    pdf_draw_line(left_margin, *page_height, right_margin, *page_height);
    pdf_advance(page_height);
    pdf_print_line_element(lx, *ly, offset[0], "Population");
    pdf_print_line_element(lx, *ly, offset[2], "Time");
    pdf_print_line_element(lx, *ly, offset[5], "Frequency");
    pdf_advance(page_height);
    pdf_draw_line(offset[2], *ly, lx+roffset[4], *page_height);
    pdf_advance(page_height);
    pdf_print_line_element(lx, *ly, offset[0], "From");
    pdf_print_line_element(lx, *ly, offset[1], "To");
    pdf_print_line_element(lx, *ly, offset[2], "Average");
    pdf_print_line_element(lx, *ly, offset[3], "Median");
    pdf_print_line_element(lx, *ly, offset[4], "SE");
    pdf_advance(page_height);
}

void
pdf_print_mighist_table (world_fmt * world, MYREAL ***migtable,
                         long *total)
{
    boolean first=TRUE;
    float lx;
    float *ly;
    float page_height;
    float page_width;
    //  float w;
    float left_margin = 55;
    float ph;
    float right_margin = pdf_contents_get_width(canvas) - 55;
    const     float offset[] = {-1, -40, 135, 191, 247, 303, 359, 415, 471};
    const     float roffset[] = {5, 50, 135, 191, 247, 303, 359, 415, 471};
    long loc, p1;
    long loci1 = world->loci == 1 ? 1 : world->loci + 1;
    char title[LINESIZE];
    sprintf (title, "Summary of %s events", 
             world->options->mighist_all ? "coalescence and migration" : "migration");
    pdf_print_section_title(&page_width, &page_height, title);
    //page_height -= 126;
    lx = 55;
    page_height -= 20;
    //    pdf_draw_line(left_margin, page_height, right_margin, page_height);
    pdf_advance(&page_height);
    ly = &page_height;
    for (loc = 0; loc < world->loci; loc++)
        total[world->loci] += total[loc];
    first = TRUE;
    for (loc = 0; loc < loci1; loc++)
      {    /* Each locus + Summary */
        
        if (loc != world->loci)
            sprintf(title, "Locus %li", loc + 1);
        else
            sprintf(title, "Over all loci");
        
        pdf_print_contents_at(lx,page_height,title);
        pdf_advance(&page_height);
        
        // table header
        if(first)
          {
            pdf_print_mighist_table_header(lx, offset, roffset, 
                                           &page_height, left_margin, right_margin);
            first=FALSE;
          }
        else
          {
            ph = pdf_contents_get_height(canvas) - 55 - LINESTRETCH;
            if(page_height >= ph)
              {
                pdf_print_mighist_table_header(lx, offset, roffset, 
                                               &page_height, left_margin, right_margin);
              }
          }
        
        for (p1 = (world->options->mighist_all ? 0 : world->numpop); p1 < world->numpop2; p1++)
          {
            ph = pdf_contents_get_height(canvas) - 55 - LINESTRETCH;
            if(page_height >= ph)
              {
                pdf_print_mighist_table_header(lx, offset, roffset, 
                                               &page_height, left_margin, right_margin);
              }
            if (migtable[loc][p1][2] == NOAVERAGE)
              {
                continue;
                //                pdf_printf(lx, ly, 'L',"%4li %4li    No migration event encountered\n",
                //           (long) migtable[loc][p1][0] + 1,
                //           (long) migtable[loc][p1][1] + 1);
              }
            else
              {
                pdf_printf(lx + roffset[0], *ly,'L', "%4li", (long) migtable[loc][p1][0] + 1);
                pdf_printf(lx + roffset[1], *ly,'L', "%4li", (long) migtable[loc][p1][1] + 1);
                pdf_printf_ralign(lx + roffset[2], *ly,"%3.5f",migtable[loc][p1][2]);
                pdf_printf_ralign(lx + roffset[3], *ly,"%3.5f",migtable[loc][p1][3]);
                pdf_printf_ralign(lx + roffset[4], *ly,"%3.5f",migtable[loc][p1][4]);
                pdf_printf_ralign(lx + roffset[5], *ly,"%3.5f",migtable[loc][p1][5] / total[loc]);
              }
            pdf_advance(&page_height);
          }
        //        pdf_draw_line(left_margin+offset[2], page_height, left_margin+offset[4], page_height);
        pdf_advance(&page_height);            
        pdf_advance(&page_height);
      }
pdf_advance(&page_height);            
}

void
pdf_print_event_table (world_fmt * world, float *migtable)
{
    boolean first=TRUE;
    float lx;
    float *ly;
    float page_height;
    float page_width;
    //  float w;
    float left_margin = 55;
    float ph;
    float right_margin = pdf_contents_get_width(canvas) - 55;
    const     float offset[] = {-1, -40, 135, 191, 247, 303, 359, 415, 471};
    const     float roffset[] = {5, 50, 135, 191, 247, 303, 359, 415, 471};
    long loc, p1;
    long loci1 = world->loci == 1 ? 1 : world->loci + 1;
    char title[LINESIZE];
    long frompop;
    long topop;
    long migwidth = loci1 * world->numpop2;
    sprintf (title, "Summary of %s events", 
             world->options->mighist_all ? "coalescence and migration" : "migration");
    pdf_print_section_title(&page_width, &page_height, title);
    //page_height -= 126;
    lx = 55;
    page_height -= 20;
    //    pdf_draw_line(left_margin, page_height, right_margin, page_height);
    pdf_advance(&page_height);
    ly = &page_height;
    first = TRUE;
    for (loc = 0; loc < loci1; loc++)
      {    /* Each locus + Summary */
        
        if (loc != world->loci)
            sprintf(title, "Locus %li", loc + 1);
        else
            sprintf(title, "Over all loci");
        
        pdf_print_contents_at(lx,page_height,title);
        pdf_advance(&page_height);
        
        // table header
        if(first)
          {
            pdf_print_mighist_table_header(lx, offset, roffset, 
                                           &page_height, left_margin, right_margin);
            first=FALSE;
          }
        else
          {
            ph = pdf_contents_get_height(canvas) - 55 - LINESTRETCH;
            if(page_height >= ph)
              {
                pdf_print_mighist_table_header(lx, offset, roffset, 
                                               &page_height, left_margin, right_margin);
              }
          }
        
        for (p1 = (world->options->mighist_all ? 0 : world->numpop); p1 < world->numpop2; p1++)
          {
            ph = pdf_contents_get_height(canvas) - 55 - LINESTRETCH;
            if(page_height >= ph)
              {
                pdf_print_mighist_table_header(lx, offset, roffset, 
                                               &page_height, left_margin, right_margin);
              }
            m2mm(p1,world->numpop,&frompop,&topop);
            pdf_printf(lx + roffset[0], *ly,'L', "%4li", frompop + 1);
            pdf_printf(lx + roffset[1], *ly,'L', "%4li", topop + 1);
            pdf_printf_ralign(lx + roffset[2], *ly,"%3.5f", migtable[loc * loci1 + p1]);
            pdf_printf_ralign(lx + roffset[3], *ly,"%3.5f",migtable[migwidth + loc * loci1 + p1]);
            pdf_printf_ralign(lx + roffset[4], *ly,"%3.5f",migtable[2*migwidth + loc * loci1 + p1]);
            pdf_printf_ralign(lx + roffset[5], *ly,"%3.5f",migtable[3*migwidth + loc * loci1 + p1]);
            pdf_advance(&page_height);
          }
        //        pdf_draw_line(left_margin+offset[2], page_height, left_margin+offset[4], page_height);
        pdf_advance(&page_height);            
        pdf_advance(&page_height);
      }
    pdf_advance(&page_height);            
}
#endif
/* END DEBUG TODO ready to remove this function?*/


///
/// fill and stroke bezier curve
void pdf_fill_stroke(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4)
{
    pdf_contents_move_to(canvas, x1, y1);
    pdf_contents_curve_to(canvas, x2, y2, x3, y3, x4, y4);
    pdf_contents_fill_stroke(canvas);
}

///
/// plot migrate logo
void pdf_migrate_logo(float x, float y, float stretch)
{
    pdf_contents_move_to(canvas, x + stretch * 18.705982, y + stretch * 21.16717);
    pdf_contents_curve_to(canvas, x + stretch * 18.705982, 
                          y + stretch * 21.234596, x + stretch * 18.779981,y + stretch * 21.289254,
                          x + stretch * 18.871266, y + stretch * 21.289254);
    pdf_contents_curve_to(canvas, x + stretch * 18.962551, 
                          y + stretch * 21.289254, x + stretch * 19.03655, y + stretch * 21.234596, 
                          x + stretch * 19.03655, y + stretch * 21.167171);
    pdf_contents_curve_to(canvas, x + stretch * 19.03655, 
                          y + stretch *  21.099745, x + stretch * 18.962551, y + stretch * 21.045087, 
                          x + stretch * 18.871266, y + stretch * 21.045087);
    pdf_contents_curve_to(canvas, x + stretch * 18.779981, 
                          y + stretch * 21.045087, x + stretch * 18.705982, y + stretch * 21.099745, 
                          x + stretch * 18.705982, y + stretch * 21.167171);
    pdf_contents_fill_stroke(canvas);
}

void pdf_migrate_logo_lines(float x, float y, float stretch)
{
    pdf_contents_move_to(canvas, x + stretch * 18.773599, y + stretch * 21.177612);
    pdf_contents_line_to(canvas, x + stretch * 18.773599, y + stretch * 20.997156);
    pdf_contents_line_to(canvas, x + stretch * 18.882535, y + stretch * 20.997156);
    pdf_contents_line_to(canvas, x + stretch * 18.882861, y + stretch * 21.177612);
    pdf_contents_stroke(canvas);
    pdf_contents_move_to(canvas, x + stretch * 18.979539,  y + stretch * 21.177612);
    pdf_contents_line_to(canvas, x + stretch * 18.980202,  y + stretch * 20.948318);
    pdf_contents_line_to(canvas, x + stretch * 19.104166,  y + stretch * 20.948318);
    pdf_contents_line_to(canvas, x + stretch * 19.10341 ,  y + stretch * 21.177612);
    pdf_contents_stroke(canvas);
    pdf_contents_move_to(canvas, x + stretch * 19.213683,  y + stretch * 21.177612);
    pdf_contents_line_to(canvas, x + stretch * 19.213103,  y + stretch * 20.891974);
    pdf_contents_line_to(canvas, x + stretch * 19.045865,  y + stretch * 20.891974);
    pdf_contents_line_to(canvas, x + stretch * 19.045865,  y + stretch * 20.948318);
    pdf_contents_stroke(canvas);
    pdf_contents_move_to(canvas, x + stretch * 19.318285,   y + stretch * 21.177612);
    pdf_contents_line_to(canvas, x + stretch * 19.318285,   y + stretch * 20.809329);
    pdf_contents_line_to(canvas, x + stretch * 19.132266,   y + stretch * 20.809329);
    pdf_contents_line_to(canvas, x + stretch * 19.132521,   y + stretch * 20.890561);
    pdf_contents_stroke(canvas);
    pdf_contents_move_to(canvas, x + stretch * 19.228543,  y + stretch * 20.645554);
    pdf_contents_line_to(canvas, x + stretch * 19.229199,  y + stretch * 20.808985);
    pdf_contents_stroke(canvas);
    pdf_contents_move_to(canvas, x + stretch * 18.829904,  y + stretch * 20.647076);
    pdf_contents_line_to(canvas, x + stretch * 18.829904,  y + stretch * 20.996422);
    pdf_contents_stroke(canvas);
    
}

///
/// print time stamp in pdf
void pdf_print_time(float startx, float *pageheight, char text[])
{
    char nowstr[LINESIZE];
    get_time(nowstr, "  %c");
    if (nowstr[0] != '\0')
        pdf_printf(startx, *pageheight, 'L',"%s %s", text, nowstr);
    pdf_advance(pageheight);
}

///
/// prints time stamp for end of run needs pretty.c-global variables firstcanvas and timestampy/// to successfully print timestamp, compare to call of pdf_print_time() function. 
void pdf_print_end_time(float *page_height)
{
    pdf_contents thiscanvas = firstcanvas;
    char nowstr[LINESIZE];
    char *title;
    title = (char *) mycalloc(LINESIZE,sizeof(char));
    
    get_time(nowstr, "  %c");
    sprintf(title,"Program finished at %s", nowstr);
    
    if (nowstr[0] != '\0')
      {    
        pdf_contents_set_font_and_size(thiscanvas, "Helvetica", 12);
        pdf_contents_begin_text(thiscanvas);
        pdf_contents_move_text_pos(thiscanvas, 55., timestampy);
        pdf_contents_show_text(thiscanvas, title);
        pdf_contents_end_text(thiscanvas);
      }
    myfree(title);
}

///
/// add a new page and set the master title
void pdf_title(char *title, float *page_height, float page_width)
{
    float w;
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);  
    pdf_print_contents_at((page_width - w)/2, *page_height - 100, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    pdf_draw_line(50, *page_height - 126, page_width-50, *page_height - 126);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    *page_height -= 126.0f;
}

void find_posterior_min_max(float *minval, float *maxval, long start, long stop, bayes_fmt *bayes, long locus)
{
    long j, pa0,pa;
    long numbins  = 0;
    long *bins = bayes->histogram[locus].bins;
    MYREAL *results = bayes->histogram[locus].results;
    MYREAL *mini = bayes->histogram[locus].minima;
    MYREAL *maxi = bayes->histogram[locus].maxima;
    float p01, p99, mm;
    float val;
    long pos01, pos99, posmm;
    // set all parameter print maxima to the same value for theta and M
    *minval  = HUGE;
    *maxval  = 0.;
    for(pa=0;pa < start; pa++)
        numbins += bins[pa];
    for(pa0=start; pa0 < stop; pa0++)
      {
        // if custom migration matrix is set to zero
        // continue
        //if(strchr("0c",bayes->custm2[pa]))
        //continue;
        if(bayes->map[pa0][1] == INVALID)
            continue;
        else
            pa = bayes->map[pa0][1];
        for(j=0;j < pa; j++)
            numbins += bins[j];
        findmoments(&results[numbins],bins[pa],&p01,&pos01,&p99,&pos99,&mm,&posmm);
        if(mini[pa] < *minval)
            *minval = mini[pa];
        val = pos99 * (maxi[pa]-mini[pa])/bins[pa];
        if(val> *maxval)
            *maxval = val;
        //      numbins += bins[pa];
      }
}

///
/// plot posterior distribution, returns page_height 
float pdf_loci_histogram(world_fmt *world)
{
    bayes_fmt * bayes  = world->bayes;
    
    long locus         = world->loci > 1 ? world->loci : 0;
    long numparam      = world->numpop2 + (world->bayes->mu);
    long z             = 0;
    long numbins       = 0;
    long numbinsall    = 0;
    long p99bins       = 0;
    long *bins         = bayes->histogram[locus].bins;
    long frompop;
    long topop;
    long pa0, pa=0;
    long rpa;
    //long j;
    
    float thetamin     = 0; 
    float migmin       = 0;
    float ratemin      = 0;
    float thetamax     = 0;
    float migmax       = 0;
    float ratemax      = 0;
    float page_height  = pdf_contents_get_height(canvas);
    float page_width   = pdf_contents_get_width(canvas);
    float lx;
    float ly;
    
    char *set50        = bayes->histogram[locus].set50;
    char *set95        = bayes->histogram[locus].set95;
    char title[100]    = "Bayesian Analysis: Posterior distribution over all loci";
    
    float delta; 
    MYREAL *results    = bayes->histogram[locus].results;
    MYREAL *mini       = bayes->histogram[locus].minima;
    MYREAL *maxi       = bayes->histogram[locus].maxima;
    
    // set the title of the section
    pdf_title(title, &page_height, page_width);
    
    switch(bayes->prettyhist)
      {
        case PRETTY_MAX:
        case PRETTY_P99:
        case PRETTY_P100:
            break;
        default:
            find_posterior_min_max(&thetamin, &thetamax, 0, world->numpop, bayes, locus);
            find_posterior_min_max(&migmin, &migmax, world->numpop, world->numpop2, bayes, locus);
            if(bayes->mu)
                find_posterior_min_max(&ratemin, &ratemax, world->numpop2, numparam, bayes, locus);
      }
    
    lx = 100;
    ly = page_height - 190;
    numbinsall = 0;
    for(pa0=0; pa0 < numparam; pa0++)
      {
        if(bayes->map[pa0][1] == INVALID)
          {
            continue;
          }
        else
          {
            pa = bayes->map[pa0][1];
          }

	if(pa < pa0)
	  continue; // does not print multiple copies for 'm' and 's'

	if(pa0>=world->numpop2)
	  {
	    rpa=world->numpop2;
	    pa = pa0;
	  }
	else
	  {
	    rpa=pa;
	  }
	numbinsall += bins[rpa];
	numbins = numbinsall - bins[rpa];
        //for(j=0; j < pa; j++)
        //  {
        //    numbins += bins[j];
        //  }
        if(pa0<world->numpop2)
          {
            m2mm(pa0, world->numpop, &frompop, &topop);
            if(frompop==topop)
              {
                pdf_print_contents_at(lx-30,ly+125, "Freq");
                symbol_Theta(lx+80, ly-25, 12, frompop+1);
              }
            else
              {
                pdf_print_contents_at(lx-30,ly+125, "Freq");
                symbol_M(lx+80, ly-25, 12, frompop+1, topop+1, world->options->usem);
              }
          }
        else
          { 
            pdf_print_contents_at(lx-30,ly+125, "Freq");
            symbol_R(lx+80, ly-25, 12, -2);
          }
        delta = (maxi[rpa] - mini[rpa])/bins[rpa];
        switch(bayes->prettyhist)
          {
            case PRETTY_MAX:
                pdf_histogram(&results[numbins],&set50[numbins], &set95[numbins], 
                              bins[pa],(float) delta, (float) mini[rpa], (float) maxi[rpa],lx,ly,187,116, FALSE);
                break;
            case PRETTY_P99:
                pdf_histogram(&results[numbins],&set50[numbins], &set95[numbins], 
                              bins[rpa],(float) delta, (float) mini[rpa], -9999,lx,ly,187,116, FALSE);
                break;
            case PRETTY_P99MAX:
                if(pa < world->numpop)
                  {
                    p99bins = (long)((thetamax-thetamin)/(float)delta);
                    pdf_histogram(&results[numbins],&set50[numbins], &set95[numbins], 
                                  p99bins,(float) delta, thetamin, thetamax,
                                  lx,ly,187,116, FALSE);
                  }
                else
                  {
                    if(pa >= world->numpop2)
                      {
                        p99bins = (long) ((ratemax - ratemin)/(float)delta);
                        pdf_histogram(&results[numbins],&set50[numbins], &set95[numbins], 
                                      p99bins,(float) delta, ratemin, ratemax,
                                      lx,ly,187,116, FALSE);
                      }
                    else
                      {
                        p99bins = (long) ((migmax - migmin)/(float)delta);
                        pdf_histogram(&results[numbins],&set50[numbins], &set95[numbins], 
                                      p99bins,(float) delta, migmin, migmax,
                                      lx,ly,187,116, FALSE);
                      }
                  }
                break;
            case PRETTY_P100:
            default:
                pdf_histogram(&results[numbins],&set50[numbins], &set95[numbins], 
                              bins[rpa],(float) delta, (float) mini[rpa], -999,lx,ly,187,116, FALSE);
                break;
                
          }
        if(z++ % 2 == 0 && world->numpop > 1)
	      {
            lx = 350;
            ly = page_height - 190;
	      }
	    else
	      {
            lx = 100;
            page_height -= 180;
            if(page_height-50 < 180)
              {
                pdf_new_page("");
                page_height = pdf_contents_get_height(canvas) - 10 ;
              }
            ly = page_height - 190;
	      }
        //numbins += bins[pa];
      }
    return page_height;
}

///
/// prints right-adjusted text, based on margin-width x and stay on the line 
void pdf_printf_right(float x, float y, char string[],...)
{
    //    float     page_width = pdf_contents_get_width(canvas);
    float w;
    char message[LINESIZE];
	char fp[LINESIZE];
	va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);	
    sprintf(fp,"%s",message);
    w = (float) pdf_contents_get_text_width(canvas, fp, NULL, NULL);
    pdf_print_contents_at(/*page_width-*/x-w, y, fp);
}


///
/// prints right-adjusted text, based on margin-width x and jumps to next line 
void pdf_printf_right_next(float x, float *y, char string[],...)
{
    float     page_width = pdf_contents_get_width(canvas);
    float w;
    char message[LINESIZE];
	char fp[LINESIZE];
	va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);	
    sprintf(fp,"%s",message);
    w = (float) pdf_contents_get_text_width(canvas, fp, NULL, NULL);
    pdf_print_contents_at(page_width-x-w, *y, fp);
    pdf_advance(y);
}

///
/// prints right-aligned mutable text. Like vprintf
void pdf_printf_ralign(float rx, float y, char string[], ...)
{
    float w;
    char message[LINESIZE];
	char fp[LINESIZE];
	va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);	
	sprintf(fp,"%s",message);
    w = (float) pdf_contents_get_text_width(canvas, fp, NULL, NULL);
    pdf_print_contents_at(rx-w, y, fp);
}



///
/// Prints aligned text mutable text; like vprintf.
/// The align parameter defines the alignment of the coordinate
/// when align = 'L' it is left aligned, 'C' center aligned, and 'R' rightaligned
void pdf_printf(float x, float y, char align, char string[], ...)
{
    
    char message[LINESIZE];
    char fp[LINESIZE];
    float w;
    va_list args;
    
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);
    
    sprintf(fp,"%s",message);
    
    w = (float) pdf_contents_get_text_width(canvas, fp, NULL, NULL);
    switch(align)
      {
        case 'C':
            w /= 2. ;
            break;
        case 'R':
            break;
        case 'L':
        default:
            w = 0.;
            break;
      }
    pdf_print_contents_at(x+w, y, fp);
}


///
/// prints mutable text and advances to next line.
void pdf_printf_next(float x, float *y, char string[], ...)
{
    char message[LINESIZE];
	char fp[LINESIZE];
	va_list args;
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);	
	sprintf(fp,"%s",message);
    pdf_print_contents_at(x, *y, fp);
    pdf_advance(y);
}

///
/// prints mutable text at x,y and changes x and y when reaching end of line
void pdf_printf_cell(float *x, float *y, float width, char string[], ...)
{
    char message[LINESIZE];
    char fp[LINESIZE];
    float w;
    //double ww;
    va_list args;
    
    va_start (args, string);
    vsprintf (message, string, args);
    va_end (args);	
    sprintf(fp,"%s",message);
    //  pdf_contents_get_char_widths(canvas, fp, &ww);
    //w = (float) ww;
    w = (float) pdf_contents_get_text_width(canvas, fp, NULL, NULL);
    if((w + *x) >=  width)
      {
        *x = 55;
        pdf_advance(y);
      }
    pdf_print_contents_at(*x, *y, fp);
    *x += 50;
}

void    pdf_print_connection_table (float * page_height, world_fmt * world, option_fmt *options, data_fmt * data)
{
    const float offset = 20;
    long i;
    long j;
    long z;
    float right_margin = pdf_contents_get_width(canvas) - 2*55 - offset;
    float left_margin = 55;
    float left_margin2 = 60;
    float lx = left_margin + 105;
    pdf_advance(page_height);		
    pdf_printf_next(left_margin, page_height,"Connection type matrix:");
    pdf_printf_next (left_margin2, page_height, "where m = average (average over a group of Thetas or M,");
    pdf_printf_next (left_margin2, page_height, "s = symmetric M, S = symmetric 4Nm,\n 0 = zero, and not estimated,");
    pdf_printf_next (left_margin2, page_height, "* = free to vary, Thetas are on diagonal\n");
    pdf_advance(page_height);		
    pdf_printf (left_margin, *page_height, 'L', "Population");
    for (i = 0; i < data->numpop; i++)
      {
        pdf_printf (lx, *page_height, 'L', "%3li", options->newpops[i]);
        lx += offset;
        if(lx > right_margin)
          {
            lx = left_margin + 110;
            pdf_advance(page_height);		                
          }
      }
    //xcode lx = left_margin + 110;
    for (i = 0; i < data->numpop; i++)
      {
        lx = left_margin + 110;
        pdf_advance(page_height);		
        pdf_printf (left_margin, *page_height, 'L', "%3li %-15.15s",options->newpops[i],  data->popnames[i]);
        for (j = 0; j < data->numpop; j++)
          {
            z = mm2m(options->newpops[j]-1,options->newpops[i]-1,world->numpop);
            pdf_printf(lx, *page_height, 'L', " %c", world->options->custm2[z]);
            lx += offset;
            if(lx > right_margin)
              {
                lx = left_margin + 90;
                pdf_advance(page_height);		                
              }
          }

      }
    pdf_advance(page_height);		
    pdf_advance(page_height);		
}


///
/// prints order of parameters for PDF output file
void pdf_print_param_order(float *page_height, world_fmt *world)
{
  float left_margin = 60;

  long pa;
  long pa0;
  long numpop = world->numpop;
  long numpop2 = world->numpop2;
  long frompop, topop;
  long frompop2, topop2;
  char *custm2 = world->options->custm2;
  boolean usem = world->options->usem;
  bayes_fmt *bayes = world->bayes;
  long numparam = world->numpop2 + world->bayes->mu;
  pdf_advance(page_height);		
  pdf_printf_next(left_margin, page_height,"Order of parameters:");
  //  fprintf(bayesfile,"# Parameter-number Parameter\n");
  for(pa0=0;pa0<numparam;pa0++)
    {
      if(bayes->map[pa0][1] != INVALID)
	{
	  if(pa0 == bayes->map[pa0][1])
	    pa = pa0;
	  else
	    pa = bayes->map[pa0][1];
	  
	  if(pa0 < numpop)
	    {
	      if((pa0 == pa) && (custm2[pa0] == '*'))
		{
		  pdf_printf(left_margin, *page_height,'L',"%4li ",pa0+1);
		  symbol_Theta(left_margin+80, *page_height,12,pa0+1);
		  pdf_printf(left_margin+230, *page_height,'L',"<displayed>");
		}  
	      else
		{
		  pdf_printf(left_margin, *page_height,'L',"%4li ",pa0+1);
		  symbol_Theta(left_margin+80, *page_height,12,pa0+1);
		  pdf_printf(left_margin+120, *page_height,'L'," = ");
		  symbol_Theta(left_margin+150, *page_height,12,pa+1);
		  pdf_printf(left_margin+190, *page_height,'L',"[%c]",custm2[pa0]);
		  if(pa==pa0)
		    pdf_printf(left_margin+230, *page_height,'L',"<displayed>");
		}
	    }
	  else
	    {
	      // do we estimate mutation rate changes?
	      if(pa0 >= numpop2)
		{
		  pdf_printf(left_margin, *page_height,'L',"%4li ",pa0+1);
		  pdf_printf(left_margin+80, *page_height,'L', "Rate");
		  pdf_printf(left_margin+230, *page_height,'L',"<displayed>");
		}
	      else
		{
		  m2mm(pa0,numpop,&frompop,&topop);
		  if((pa0==pa) && (custm2[pa0]=='*'))
		    {
		      if(usem)
			{
			  pdf_printf(left_margin, *page_height,'L',"%4li ",pa0+1);
			  symbol_M(left_margin+80, *page_height,12,frompop+1,topop+1,TRUE);
			  pdf_printf(left_margin+230, *page_height,'L',"<displayed>");
			}
		      else  
			{
			  pdf_printf(left_margin, *page_height,'L',"%4li ",pa0+1);
			  symbol_M(left_margin+80, *page_height,12,frompop+1,topop+1,FALSE);
			  pdf_printf(left_margin+120, *page_height,'L'," = ");
			  symbol_Theta(left_margin+150, *page_height,12,topop+1);
			  symbol_M(left_margin+190, *page_height,12,frompop+1,topop+1,TRUE);
			  pdf_printf(left_margin+230, *page_height,'L',"<displayed>");
			}
		    }
		  else
		    {
		      m2mm(pa,numpop,&frompop2,&topop2);
		      if(usem)
			{
			  pdf_printf(left_margin, *page_height,'L',"%4li ",pa0+1);
			  symbol_M(left_margin+80, *page_height,12,frompop+1,topop+1,TRUE);
			  pdf_printf(left_margin+120, *page_height,'L'," = ");
			  symbol_M(left_margin+150, *page_height,12,frompop2+1,topop2+1,TRUE);
			  pdf_printf(left_margin+190, *page_height,'L',"[%c]",custm2[pa0]);
			  if(pa==pa0)
			    pdf_printf(left_margin+230, *page_height,'L',"<displayed>");
			}	      
		      else  
			{
			  symbol_M(left_margin+80, *page_height,12,frompop+1,topop+1,FALSE);
			  pdf_printf(left_margin+120, *page_height,'L'," = ");
			  symbol_Theta(left_margin+150, *page_height,12,topop+1);
			  symbol_M(left_margin+190, *page_height,12,frompop+1,topop+1,TRUE);
			  pdf_printf(left_margin+240, *page_height,'L'," = ");
			  symbol_M(left_margin+290, *page_height,12,frompop2+1,topop2+1,FALSE);
			  pdf_printf(left_margin+330, *page_height,'L'," = ");
			  symbol_Theta(left_margin+350, *page_height,12,topop2+1);
			  symbol_M(left_margin+390, *page_height,12,frompop2+1,topop2+1,TRUE);
			  pdf_printf(left_margin+430, *page_height,'L',"[%c]",custm2[pa0]);
			  if(pa==pa0)
			    pdf_printf(left_margin+460, *page_height,'L',"<displayed>");

			}
		    }
		}
	    }
	  pdf_advance(page_height);
	}
    }
  pdf_advance(page_height);		
  pdf_advance(page_height);		
}



void    pdf_print_distance_table (float * page_height, world_fmt * world, option_fmt * options, data_fmt * data)
{
    const float offset = 60;
    long i;
    long j;
    float right_margin = pdf_contents_get_width(canvas) - 2*55 - offset;
    float left_margin = 55;
    float lx=left_margin + 80;
    
    if(!options->geo)
        return;
    
    pdf_advance(page_height);		
    pdf_printf_next(left_margin, page_height, "Distance among populations:");
    pdf_advance(page_height);		
    pdf_printf (left_margin, *page_height, 'L', "Population");
    for (i = 0; i < data->numpop; i++)
      {
        pdf_printf (lx, *page_height, 'L', "%3li", options->newpops[i]);
        lx += offset;
        if(lx > right_margin)
          {
            lx = left_margin + 80;
            pdf_advance(page_height);		                
          }
      }
    for (i = 0; i < data->numpop; i++)
      {
        lx = left_margin + 80;
        pdf_advance(page_height);		
        pdf_printf (left_margin, *page_height, 'L', "%3li %-15.15s",options->newpops[i],  data->popnames[i]);
        for (j = 0; j < data->numpop; j++)
          {
            pdf_printf(lx, *page_height, 'L', " %10.4f ", data->ogeo[j][i]);
            lx += offset;
            if(lx > right_margin)
              {
                lx = left_margin + 80;
                pdf_advance(page_height);		                
              }
          }
      }
    pdf_advance(page_height);		
    pdf_advance(page_height);		
}



void pdf_print_options(world_fmt * world, option_fmt *options, data_fmt * data, float *orig_page_height, float *orig_left_margin)
{
    
    float page_width;
    float page_height = *orig_page_height;
    float left_margin = *orig_left_margin;
    float right_margin = 300. ;
    float w;
    float width[7]={0, 100, 160, 240, 320, 400, 480};
    float width2[2]={0, 240};
    char *title = "Options";
    long i, j, tt, ii;
    char mytext[LINESIZE];
    char mytext1[LINESIZE];
    char mytext2[LINESIZE];
    char mytext3[LINESIZE];
    char mytext4[LINESIZE];
    char seedgen[LINESIZE], spacer[LINESIZE];
    char paramtgen[LINESIZE], parammgen[LINESIZE];
    
    //xcode page_width = pdf_contents_get_width(canvas) - 120;
    
    if (options->datatype != 'g')
      {
        switch ((short) options->autoseed)
          {
            case AUTO:
                strcpy (seedgen, "with internal timer");
                strcpy (spacer, "  ");
                break;
            case NOAUTOSELF:
                strcpy (seedgen, "from parmfile");
                strcpy (spacer, "      ");
                break;
            case NOAUTO:
                strcpy (seedgen, "from seedfile");
                strcpy (spacer, "      ");
                break;
            default:
                strcpy (seedgen, "ERROR");
                strcpy (spacer, " ");
                break;
          }
        switch (options->thetaguess)
          {
            case OWN:
                strcpy (paramtgen, "from guessed values");
                break;
            case FST:
                strcpy (paramtgen, "from the FST-calculation");
                break;
            case NRANDOMESTIMATE:
                strcpy (paramtgen, "RANDOM start value from N(mean,std)");
                break;
            case URANDOMESTIMATE:
                strcpy (paramtgen, "RANDOM start value from U(min,msx)");
                break;
            default:
                strcpy (paramtgen, "ERROR");
                break;
          }
        switch (options->migrguess)
          {
            case OWN:
                strcpy (parammgen, "from guessed values");
                break;
            case FST:
                strcpy (parammgen, "from the FST-calculation");
                break;
            case NRANDOMESTIMATE:
                strcpy (parammgen, "RANDOM start value from N(mean,std)");
                break;
            case URANDOMESTIMATE:
                strcpy (parammgen, "RANDOM start value from U(min,max)");
                break;
            default:
                strcpy (parammgen, "ERROR");
                break;
          }
      }
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 18);    
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    page_width = pdf_contents_get_width(canvas);
    right_margin = page_width - 55;
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    pdf_contents_set_rgb_fill(canvas, PdfRGBColor(0, 0, 0));
    pdf_contents_set_line_width(canvas, 1);
    pdf_print_contents_at((page_width - w)/2, page_height, title);
    page_height -= 20;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    page_height -= 20;
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);    
    
    switch (options->datatype)
      {
        case 'a':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Allelic data");
            pdf_print_contents_at(left_margin, page_height,"Missing data:");
            pdf_printf_right_next(left_margin, &page_height,"%s\n",options->include_unknown ? "included" : "not included");
            break;
        case 'b':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Microsatellite data [Brownian motion]\n");
            pdf_print_contents_at(left_margin, page_height,"Missing data:");
            pdf_printf_right_next(left_margin, &page_height,"%s\n",options->include_unknown ? "included" : "not included");
            break;
        case 'm':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
	    if(options->msat_option == SINGLESTEP)
	      pdf_printf_right_next(left_margin, &page_height,"Microsatellite data [Singlestep model]\n");
	    else
	      pdf_printf_right_next(left_margin, &page_height,"Microsatellite data [Multistep model (Tune=%f, P_increase=%f)]\n",
							       options->msat_tuning[0], options->msat_tuning[1]);
            pdf_print_contents_at(left_margin, page_height,"Missing data:");
            pdf_printf_right_next(left_margin, &page_height,"%s\n",options->include_unknown ? "included" : "not included");
            break;
        case 's':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"DNA sequence data\n");
            break;
        case 'n':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Single nucleotide polymorphism data\n");
            break;
        case 'h':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Single nucleotide polymorphism data(Hapmap formatting)\n");
            break;
        case 'u':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Single nucleotide polymorphism data (PANEL)\n");
            break;
        case 'f':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Ancestral state method\n");
            break;
        case 'g':
            pdf_print_contents_at(left_margin, page_height,"Datatype:");
            pdf_printf_right_next(left_margin, &page_height,"Genealogy summary of an older run\n");
            break;
      }

    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height,"Inheritance scalers in use for Thetas:");
    left_margin += pdf_contents_get_text_width(canvas, "Inheritance scalers in use for Thetas: ", NULL, NULL);
    for(i=0;i<data->loci;i++)
      {
	if(i<options->inheritance_scalars_numalloc)
	  pdf_printf_cell(&left_margin, &page_height, page_width-left_margin, "%2.2f", options->inheritance_scalars[i]);
	else
	  pdf_printf_cell(&left_margin, &page_height, page_width-left_margin, "%2.2f", options->inheritance_scalars[options->inheritance_scalars_numalloc-1]);
      }
    pdf_advance(&page_height);
    left_margin = *orig_left_margin ;
    pdf_print_contents_at(left_margin, page_height,"[Each Theta uses the (true) ineritance scalar of the first locus as a reference]");
    pdf_advance(&page_height);
    
 
    if(options->randomsubset > 0)
      {
        pdf_advance(&page_height);
	pdf_print_contents_at(left_margin, page_height,"Data set was subsampled: used a random sample of size:");
	pdf_printf_right_next(left_margin, &page_height,"%5li", options->randomsubset);
        pdf_advance(&page_height);
      }

    if (options->datatype != 'g')
      {
        pdf_print_contents_at(left_margin, page_height,"Random number seed:");
        pdf_printf_right_next(left_margin, &page_height,"(%s)%s%20li", seedgen, " ",
                              options->saveseed);
        pdf_printf_next(left_margin, &page_height,"Start parameters:");
        pdf_advance(&page_height);
        pdf_print_contents_at(left_margin, page_height,"Theta values were generated");
        pdf_printf_right_next(left_margin, &page_height," %s", paramtgen);
        if (options->thetaguess == OWN)
          {
            //xcode ii = 0;
            left_margin += 10;
            pdf_print_contents_at(left_margin, page_height,"Theta = ");
            left_margin += pdf_contents_get_text_width(canvas, "Theta = ", NULL, NULL);
            for (i = 0; i < options->numthetag; i++)
              {
                pdf_printf_cell(&left_margin, &page_height, page_width, "%.5f", options->thetag[i]);
              }
          }
        
        left_margin = *orig_left_margin ;
        pdf_advance(&page_height);
        
        if(options->usem)
            pdf_print_contents_at(left_margin, page_height,"M values were generated");
        else
            pdf_print_contents_at(left_margin, page_height,"xNm values were generated");
        
        pdf_printf_right_next(left_margin, &page_height," %s", parammgen);
        
        left_margin = *orig_left_margin + 10;
        if (options->migrguess == OWN)
          {
            tt = 0;
            if (options->usem)
                pdf_print_contents_at(left_margin, page_height,"M-matrix:");
            else
                pdf_print_contents_at(left_margin, page_height,"xNm-matrix:");
            pdf_advance(&page_height);
            if (options->nummg == 1)
              {
                pdf_printf_cell(&left_margin, &page_height, page_width, 
                                "%5.2f [all are the same]", options->mg[tt]);//xcode replaced tt++ with tt
              }
            else
              {
                for (i = 0; i < world->numpop; i++)
                  {
                   //xcode  ii = 0;
                    for (j = 0; j < world->numpop; j++)
                      {
                        if (i != j)
                          {
                            pdf_printf_cell(&left_margin, &page_height, page_width, 
                                            " %5.*f, ",options->usem ? 1 : 5, options->mg[tt++]);
                          }
                        else
                          {
                            pdf_printf_cell(&left_margin, 
                                            &page_height, page_width, "   -   ");
                          }
                      }
                    left_margin = *orig_left_margin;
                    pdf_advance(&page_height);
                  }
                left_margin = *orig_left_margin;
                pdf_advance(&page_height);
              }
          }
      }
    pdf_print_connection_table (&page_height, world, options, data);
    pdf_print_distance_table (&page_height, world, options, data);
    pdf_print_param_order(&page_height, world);
    left_margin = *orig_left_margin;
    if (options->gamma)
      {
        pdf_print_contents_at(left_margin, page_height,"Mutation rate among loci:");
        pdf_printf_right_next(left_margin, &page_height,"from a Gamma distribution\n");
        pdf_printf_right_next(left_margin, &page_height,"Initial scale parameter alpha = %f\n",
                              options->alphavalue);
        if (options->custm[world->numpop2] == 'c')
          {
            pdf_printf_next(left_margin + 300, &page_height,"and is constant [will not be estimated]\n");
          }
      }
    else
      {
        if(options->bayesmurates)
          {
            pdf_printf_right_next(left_margin, &page_height,"Mutation rate is estimated %s", world->loci > 1 ?
                                  "for all loci" : "");
          }
        else
          {
            pdf_print_contents_at(left_margin, page_height,"Mutation rate among loci:");
            if (options->murates && world->loci > 1)
              {
                if(options->murates_fromdata)
                    pdf_printf_right_next(left_margin, &page_height,"Varying ([crudely] estimated from data)");
                else
                    pdf_printf_right_next(left_margin, &page_height,"Varying (user input)");
                pdf_print_contents_at(left_margin + 10, page_height,"Rates per locus: ");
                ii=0;
                for (i = 0; i < world->loci-1; i++)
                  {
                    sprintf(mytext,"%.5f, ", options->mu_rates[i]);
                    if (i % 6 == 5)
                      {
                        ii=0;
                        pdf_advance(&page_height);
                        pdf_print_contents_at(left_margin + 100 + ii * 60, page_height,mytext);
                      }
                    else
                      {
                        pdf_print_contents_at(left_margin + 100 + ii * 60, page_height,mytext);
                      }
                    ii++;
                  }
                if(i % 6 == 5)
                  {
                    ii=0;
                    pdf_advance(&page_height);
                  }
                pdf_printf_next(left_margin + 100 + ii * 60, &page_height,"%.5f", options->mu_rates[i]);
              }
            else
                pdf_printf_right_next(left_margin, &page_height,"Mutation rate is constant %s", world->loci > 1 ?
                                      "for all loci" : "");
          }
      }
#ifdef UEP
    if (options->uep)
      {
        pdf_printf_next(left_margin, &page_height,"0/1 polymorphism analysis, with 0/1 data in file:");
        pdf_printf_right_next(left_margin, &page_height,"%s\n", options->uepfilename);
        pdf_printf_right_next(left_margin, &page_height,"with forward mutation rate %f*mu\n",
                              options->uepmu);
        pdf_printf_right_next(left_margin, &page_height,"with back mutation rate %f*mu\n",
                              options->uepnu);
        pdf_printf_right_next(left_margin, &page_height,"with base frequencies \"0\"=%f and \"1\"=%f\n",
                              options->uepfreq0,options->uepfreq1);
      }
#endif
    pdf_advance(&page_height);
    
    if(options->bayes_infer)
      {
        pdf_print_contents_at(left_margin, page_height,"Analysis strategy:");
        pdf_printf_right_next(left_margin, &page_height,"Bayesian inference");
        pdf_advance(&page_height);
        
        pdf_print_contents_at(left_margin, page_height,"Proposal distributions for parameter");
        pdf_advance(&page_height);
        pdf_print_tableline(&page_height, width2, "%s %s", "Parameter", "Proposal");
        pdf_advance(&page_height);
        pdf_print_tableline(&page_height, width2, "%s %s", "Theta",
                            is_proposaltype(options->slice_sampling[THETAPRIOR]));
        pdf_advance(&page_height);
        pdf_print_tableline(&page_height, width2, "%s %s", options->usem ? "M" : "xNm",
                            is_proposaltype(options->slice_sampling[MIGPRIOR]));
        pdf_advance(&page_height);
        if(options->bayesmurates)
          {
            pdf_print_tableline(&page_height, width2, "%s %s", "Rate",
                                is_proposaltype(options->slice_sampling[RATEPRIOR]));
            pdf_advance(&page_height);
          }
        pdf_advance(&page_height);
        
        pdf_print_contents_at(left_margin, page_height,"Prior distribution for parameter");
        pdf_advance(&page_height);
        pdf_print_tableline(&page_height, width, "%s %s %s %s %s %s %s", "Parameter", "Prior", "Minimum",  "Mean*",  "Maximum", "Delta", "Bins");
        pdf_advance(&page_height);
        pdf_print_tableline(&page_height, width, "%s %s %s %s %s %s %s", "Theta",
                            is_priortype(options->bayesprior[THETAPRIOR]),
                            show_priormin(mytext1, options->bayespriortheta,options->bayesprior[THETAPRIOR]),
                            show_priormean(mytext2, options->bayespriortheta,options->bayesprior[THETAPRIOR]),
                            show_priormax(mytext3, options->bayespriortheta,options->bayesprior[THETAPRIOR]),
                            show_priordelta(mytext4,options->bayespriortheta,options->bayesprior[THETAPRIOR]),
                            show_priorbins(mytext, options->bayespriortheta,options->bayesprior[THETAPRIOR]));
        pdf_advance(&page_height);
        
        pdf_print_tableline(&page_height, width, "%s %s %s %s %s %s %s", options->usem ? "M" : "xNm",
                            is_priortype(options->bayesprior[MIGPRIOR]),
                            show_priormin(mytext1, options->bayespriorm,options->bayesprior[MIGPRIOR]),
                            show_priormean(mytext2, options->bayespriorm,options->bayesprior[MIGPRIOR]),
                            show_priormax(mytext3, options->bayespriorm,options->bayesprior[MIGPRIOR]),
                            show_priordelta(mytext4, options->bayespriorm,options->bayesprior[MIGPRIOR]),
                            show_priorbins(mytext, options->bayespriorm,options->bayesprior[MIGPRIOR]));
        pdf_advance(&page_height);
        if(options->bayesmurates)
          {
            pdf_print_tableline(&page_height, width, "%s %s %s %s %s %s %s", "Rate modif.",
                                is_priortype(options->bayesprior[RATEPRIOR]),
				show_priormin(mytext1, options->bayespriorrate,options->bayesprior[RATEPRIOR]),
				show_priormean(mytext2, options->bayespriorrate,options->bayesprior[RATEPRIOR]),
				show_priormax(mytext3, options->bayespriorrate,options->bayesprior[RATEPRIOR]),
				show_priordelta(mytext4, options->bayespriorrate,options->bayesprior[RATEPRIOR]),
				show_priorbins(mytext, options->bayespriorrate,options->bayesprior[RATEPRIOR]));
            pdf_advance(&page_height);
          }
        pdf_advance(&page_height);
      }
    else
      {
        pdf_print_contents_at(left_margin, page_height,"Analysis strategy is");
        pdf_printf_right_next(left_margin, &page_height,"Maximum likelihood");
      }
    
    pdf_advance(&page_height);
    
    if (options->datatype != 'g')
      {
        pdf_print_contents_at(left_margin, page_height,"Markov chain settings:");
        if(!options->bayes_infer)
            pdf_print_contents_at(left_margin + 300, page_height,"Short chain");
        pdf_printf_right_next(left_margin, &page_height, "Long chain");
        pdf_print_contents_at(left_margin, page_height,"Number of chains");
        if(!options->bayes_infer)
            pdf_printf_ralign(left_margin + 340, page_height,"%20li", options->schains);
        pdf_printf_right_next(left_margin+10, &page_height,"%20li", options->lchains);
        
        pdf_print_contents_at(left_margin + 10, page_height,"Recorded steps [a]");
        if(!options->bayes_infer)
            pdf_printf_ralign(left_margin + 340, page_height,"%20li", options->ssteps);
        pdf_printf_right_next(left_margin + 10, &page_height,"%20li", options->lsteps);
        
        pdf_print_contents_at(left_margin + 10, page_height,"Increment (record every x step [b]");
        if(!options->bayes_infer)
            pdf_printf_ralign(left_margin + 340, page_height,"%20li", options->sincrement);
        pdf_printf_right_next(left_margin+10, &page_height,"%20li", options->lincrement);
        if(options->bayes_infer)
          {
            pdf_print_contents_at(left_margin + 10, page_height,"Number of concurrent chains (replicates) [c]");
            pdf_printf_ralign(left_margin + 340, page_height," ");
            pdf_printf_right_next(left_margin + 10, &page_height,"%20li", options->replicate ? 
                                  options->replicatenum : 1);
            pdf_print_contents_at(left_margin + 10, page_height,"Visited (sampled) parameter values [a*b*c]");
            // pdf_printf_ralign(left_margin + 340, page_height,"%20li", 
            //                  options->ssteps * options->sincrement*(options->replicate ? options->replicatenum : 1));
            pdf_printf_right_next(left_margin + 10, &page_height,"%20li", options->lsteps * options->lincrement*(options->replicate ? options->replicatenum : 1));
          } 
        else
          {
            pdf_print_contents_at(left_margin + 10, page_height,"Visited (sampled) genealogies [a*b]");
            pdf_printf_ralign(left_margin + 340, page_height,"%20li", options->ssteps * options->sincrement);
            pdf_printf_right_next(left_margin + 10, &page_height,"%20li", options->lsteps * options->lincrement);
            
          }
        if (options->burn_in > 0)
          {
            pdf_print_contents_at(left_margin + 10, page_height,"Number of discard trees per chain (burn-in)");
            if(!options->bayes_infer)
                pdf_printf(left_margin + 320, page_height,'L', "%c%li",  options->burnin_autostop, (long) options->burn_in);
            pdf_printf_right_next(left_margin + 10, &page_height,"%c%li", options->burnin_autostop, (long) options->burn_in);        
          }
        if (options->movingsteps)
          {
            pdf_print_contents_at(left_margin + 10, page_height,"Forcing percentage of new genealogies");
             //    pdf_printf(left_margin + 360, page_height,"%5.2f",   (MYREAL) options->acceptfreq);
            pdf_printf_right_next(left_margin + 10, &page_height,"%5.2f",   (MYREAL) options->acceptfreq);
          }
        if (options->lcepsilon < LONGCHAINEPSILON)
          {
            pdf_print_contents_at(left_margin + 10, page_height,"Forcing parameter-likelihood improvement");
            //    pdf_printf(left_margin + 360, page_height,"%20.5f",   (MYREAL) options->lcepsilon);
            pdf_printf_right_next(left_margin + 10, &page_height,"%20.5f",   (MYREAL) options->lcepsilon);
          }
        
        pdf_advance(&page_height);
        
        if(options->replicate && !options->bayes_infer)
            pdf_printf_next(left_margin, &page_height,"Multiple Markov chains:");
        else
          {
            if(options->heating > 0)
                pdf_printf_next(left_margin, &page_height,"Multiple Markov chains:");
          }
        if (options->replicate && !options->bayes_infer)
          {
            pdf_print_contents_at(left_margin + 10, page_height,"Averaging over replicates");
            if (options->replicatenum == 0)
                pdf_printf_right_next(left_margin, &page_height,"Long chains only");
            else
                pdf_printf_right_next(left_margin, &page_height,"Over indepedent %li replicates", options->replicatenum);
          }
        
        if (options->heating > 0)
          {
            pdf_printf(left_margin+10, page_height,'L', "%s heating scheme", options->adaptiveheat!=NOTADAPTIVE ? ( options->adaptiveheat==STANDARD ? "Adaptive_standard" : "Bounded_adaptive") : "Static");
            pdf_printf_right_next(left_margin, &page_height, "%li chains with %s temperatures",
                                  options->heated_chains, options->adaptiveheat!=NOTADAPTIVE ? "start values" : "" );
            for (i = options->heated_chains - 1; i >= 0; i--)
                pdf_printf_right(right_margin - i * 50, page_height,"%5.2f ", options->heat[i]);
            pdf_advance(&page_height);
            
            pdf_printf_right_next(left_margin, &page_height,"Swapping interval is %li\n",
                                  options->heating_interval);
          }        
      }
    
    pdf_advance(&page_height);
    pdf_printf_next(left_margin, &page_height,"Print options:\n");
    if (options->datatype != 'g')
      {
        pdf_print_contents_at(left_margin + 10, page_height,"Data file:");
        pdf_printf_right_next(left_margin, &page_height,"%s", options->infilename);
        pdf_print_contents_at(left_margin + 10, page_height,"Output file:");
        pdf_printf_right_next(left_margin, &page_height,"%s", options->outfilename);
        if(options->writesum)
          {
            pdf_print_contents_at(left_margin + 10, page_height,"Summary of genealogies for further run:");
            pdf_printf_right_next(left_margin, &page_height,"%s", options->sumfilename);
          }
        if(options->writelog)
          {
            pdf_print_contents_at(left_margin + 10, page_height,"Log file:");
            pdf_printf_right_next(left_margin, &page_height,"%s", options->logfilename);
          }                              
        
        if(options->bayes_infer)
          {
            pdf_print_contents_at(left_margin + 10, page_height,"Posterior distribution raw histogram file:");
            pdf_printf_right_next(left_margin, &page_height,"%s\n", options->bayesfilename);
          }
        
        pdf_print_contents_at(left_margin + 10, page_height,"Print data:");
        pdf_printf_right_next(left_margin, &page_height,"%45.45s", options->printdata ? "Yes" : "No");
        
        pdf_printf(left_margin + 10, page_height,'L', "Print genealogies [only some for some data type]:");
        switch (options->treeprint)
          {
            case _NONE:
                pdf_printf_right_next(left_margin, &page_height,"None");
                break;
            case ALL:
                pdf_printf_right_next(left_margin, &page_height,"Yes, all");
                break;
            case LASTCHAIN:
                pdf_printf_right_next(left_margin, &page_height,"Yes, only those in last chain");
                break;
            case BEST:
                pdf_printf_right_next(left_margin, &page_height,"Yes, only the best");
                break;
          }
        if (options->mighist)
          {
            pdf_print_contents_at(left_margin + 10, page_height,"Histogram of the frequency of migration events");
            pdf_printf_right_next(left_margin, &page_height,"%s", options->mighistfilename);
          }        
      }
    else
      {
        pdf_print_contents_at(left_margin + 10, page_height,"Data file:");
        pdf_printf_right_next(left_margin, &page_height,"%s", options->infilename);
        pdf_print_contents_at(left_margin + 10, page_height,"Output file:");
        pdf_printf_right_next(left_margin, &page_height,"%s", options->outfilename);        
      }
    if(!options->bayes_infer)
      {
        pdf_print_contents_at(left_margin + 10, page_height,"Plot log(likelihood) surface:");
        if (options->plot)
          {
            switch (options->plotmethod)
              {
                case PLOTALL:
                    sprintf (mytext, "Yes, to outfile and %s", options->mathfilename);
                    break;
                default:
                    strcpy (mytext, "Yes, to outfile");
                    break;
              }
            pdf_printf_right_next(left_margin, &page_height,"%s\n", mytext);
            pdf_printf_right_next(left_margin, &page_height,
                                  "Parameter: %s, Scale: %s, Intervals: %li\n",
                                  options->plotvar == PLOT4NM ? "{Theta, 4Nm}" : "{Theta, M}",
                                  options->plotscale == PLOTSCALELOG ? "Log10" : "Standard",
                                  options->plotintervals);
            pdf_printf_right_next(left_margin, &page_height,"Ranges: X-%5.5s: %f - %f\n",
                                  options->plotvar == PLOT4NM ? "4Nm" : "M",
                                  options->plotrange[0], options->plotrange[1]);
            pdf_printf_right_next(left_margin, &page_height,"Ranges: Y-%5.5s: %f - %f\n", "Theta",
                                  options->plotrange[2], options->plotrange[3]);
          }
        else
          {
            pdf_printf_right_next(left_margin, &page_height,"No");
          }
        switch (options->profile)
          {
            case _NONE:
                strcpy (mytext, "No");
                break;
            case ALL:
                strcpy (mytext, "Yes, tables and summary");
                break;
            case TABLES:
                strcpy (mytext, "Yes, tables");
                break;
            case SUMMARY:
                strcpy (mytext, "Yes, summary");
                break;
          }
        pdf_printf(left_margin+10, page_height,'L', "Profile likelihood:");
        pdf_printf_right_next(left_margin, &page_height,"%s", mytext);
        
        if (options->profile != _NONE)
          {
            switch (options->profilemethod)
              {
                case 'p':
                    pdf_printf_right_next(left_margin, &page_height,"Percentile method\n");
                    break;
                case 'q':
                    pdf_printf_right_next(left_margin, &page_height,"Quick method\n");
                    break;
                case 'f':
                    pdf_printf_right_next(left_margin, &page_height,"Fast method\n");
                    break;
                case 'd':
                    pdf_printf_right_next(left_margin, &page_height,"Discrete method\n");
                    break;
                case 's':
                    pdf_printf_right_next(left_margin, &page_height,"Spline method\n");
                    break;
                default:
                    pdf_printf_right_next(left_margin, &page_height,"UNKOWN method????\n");
                    break;
              }
            pdf_printf_right_next(left_margin, &page_height,"with df=%li and for Theta and %s",
                                  options->df, options->profileparamtype ? "M=m/mu" : "4Nm");
          }
      }
    *orig_page_height = page_height;
    *orig_left_margin = left_margin;
}


void pdf_print_data_summary(world_fmt * world, option_fmt *options, data_fmt * data, 
                            float *orig_page_height, float *orig_left_margin)
{
    long locus;
    long pop;
    long numind;
    long nummiss;
    char dstring[LINESIZE];
    long *total;
    long *totalmiss;
    char *title = "Data summary";
    //    float w;
    float page_height = *orig_page_height;
    float left_margin = *orig_left_margin;
    float page_width;
    
    
    float col1 = left_margin + 300;
    float col2 = left_margin + 320;
    float col3 = left_margin + 400;
    
    
    total = (long *) mycalloc(data->loci,sizeof(long));
    totalmiss = (long *) mycalloc(data->loci,sizeof(long));
    // setup new page and title
    pdf_print_section_title(&page_width, &page_height, title);
    switch (options->datatype)
      {
        case 'a':
            strcpy (dstring, "Allelic data");
            break;
        case 'f':
        case 's':
            strcpy (dstring, "Sequence data");
            break;
        case 'b':
        case 'm':
            strcpy (dstring, "Microsatellite data");
            break;
        case 'n':
      case 'h':
        case 'u':
            strcpy (dstring, "SNP data");
            break;
        default:
            strcpy (dstring, "Unknown data [ERROR]");
            break;
      }
    pdf_print_contents_at(left_margin, page_height,"Datatype:");
    pdf_printf_right_next(left_margin, &page_height,"%s", dstring);
  if(!data->has_repeats)
    {
      pdf_printf_right_next(left_margin, &page_height,"[Fragment length is translated to repeats]");
    }
    pdf_print_contents_at(left_margin, page_height,"Number of loci:");
    pdf_printf_right_next(left_margin, &page_height,"%li", data->loci);
    pdf_advance(&page_height);
    if(options->has_datefile)
      {
	//fprintf (file, "Sample dates:          %s\n", world->datefile);
	pdf_print_contents_at(left_margin, page_height,"Sample dates:");
	pdf_printf_right_next(left_margin, &page_height,"%s", options->datefilename);
	
	//fprintf (file, "Generations per year:  %s\n", options->generation_year);
	pdf_print_contents_at(left_margin, page_height,"Generations per year:");
	pdf_printf_right_next(left_margin, &page_height,"%f", (float) options->generation_year);
	
	//fprintf (file, "Mutationrate per year: %s", options->mutationrate_year[0]);
	pdf_print_contents_at(left_margin, page_height,"Mutationrate per year:");
	pdf_printf_right_next(left_margin, &page_height,"%f",  options->mutationrate_year[0]);
	for(locus=1; locus < options->mutationrate_year_numalloc; locus++)
	  {
	    pdf_printf_right_next(left_margin, &page_height,"%f",  options->mutationrate_year[locus]);
	  }
      }
    pdf_advance(&page_height);
    pdf_print_contents_at(left_margin, page_height,"Population");
    pdf_print_contents_at(col1, page_height,"Locus");
    pdf_print_contents_at(col3, page_height, "Gene copies");    
    pdf_advance(&page_height);
    if (!strchr (SEQUENCETYPES, options->datatype))
      {
        pdf_print_contents_at(col3, page_height,"data");
        pdf_printf_right_next(left_margin, &page_height,"(missing)");
      }
    
    for (pop = 0; pop < data->numpop; pop++)
      {
        if (!strchr (SEQUENCETYPES, options->datatype))
          {
            nummiss = find_missing(data,pop,0);
            numind = data->numalleles[pop][0] - nummiss;
            pdf_printf(left_margin, page_height,'L', "%li %s", options->newpops[pop], data->popnames[pop]);
            pdf_printf_ralign(col2, page_height,"1");
            pdf_printf_ralign(col3, page_height,"%li",numind);
            pdf_printf_right_next(left_margin+10, &page_height,"(%li)", nummiss);
          }
        else
          {
            nummiss = 0;
            numind = (options->randomsubset > 0 && options->randomsubset < data->numind[pop][0]) ? options->randomsubset : data->numind[pop][0];
            pdf_printf(left_margin, page_height, 'L', "%li %s",options->newpops[pop], data->popnames[pop]);
            pdf_printf_ralign(col2, page_height,"1");
            pdf_printf_ralign(col3, page_height,"%li", numind);
            pdf_advance(&page_height);
          }
        total[0] += numind;
        totalmiss[0] += nummiss;
        
        for(locus=1; locus< data->loci; locus++)
          {
            if (!strchr (SEQUENCETYPES, options->datatype))
              {
                nummiss = find_missing(data,pop,locus);
                numind = data->numalleles[pop][locus] - nummiss;
                pdf_printf_ralign(col2, page_height,"%li",locus+1);
                pdf_printf_ralign(col3, page_height,"%li",numind);
                pdf_printf_right_next(left_margin+10, &page_height,"(%li)", nummiss);
              }
            else
              {
                nummiss=0;
                numind = data->numind[pop][locus];
                pdf_printf_ralign(col2, page_height,"%li",locus+1);
                pdf_printf_ralign(col3, page_height,"%li",numind);
                pdf_advance(&page_height);
              }
            total[locus] += numind;
            totalmiss[locus] += nummiss;
          }
      }
    pdf_printf(left_margin, page_height, 'L', 
               "Total of all populations");
    pdf_printf_ralign(col2, page_height,"1");
    pdf_printf_ralign(col3, page_height,"%li",total[0]);
    if (!strchr (SEQUENCETYPES, options->datatype))
      {
        pdf_printf_right_next(left_margin+10, &page_height,"(%li)", totalmiss[0]);
        for(locus=1; locus< data->loci; locus++)
          {
            pdf_printf_ralign(col2, page_height,"%li",locus+1);
            pdf_printf_ralign(col3, page_height,"%li",total[locus]);
            pdf_printf_right_next(left_margin+10, &page_height,"(%li)", totalmiss[locus]);
          }
      }
    else
      {
        pdf_advance(&page_height);
        for(locus=1; locus< data->loci; locus++)
          {
            pdf_printf_ralign(col2, page_height,"%li",locus+1);
            pdf_printf_ralign(col3, page_height,"%li",total[locus]);
            pdf_advance(&page_height);
          }
      }    
    myfree(total);
    myfree(totalmiss);
    *orig_page_height = page_height;
    *orig_left_margin = left_margin;
}

///
/// \param[in] *fmt format string identical to a fomrat string in printf()
/// \reval int returns the number of elements by counting the %
int count_elements(char *fmt)
{
    int element;
    int count = 0;
    for(element=0; element < (int) strlen(fmt); element++)
      {
        if(fmt[element]=='%')
            count++;
      }
    return count;
}


///
/// Print  a line in a table in the PDF file at a specific height on the page and with a given width
/// \param *page_height gets modified by this function, page_height specifies the Y coordinate
/// \param *width contains a list of floating points numbers that specify the columns, negative numbers
/// mean left-adjusted and positive numbers are right-adjusted
/// \param *fmt format string identical to printf (number % needs to match the element in width)
/// \param ...  parameters to print
void pdf_print_tableline(float *page_height, float *width, char *fmt, ...)
{
    boolean start=FALSE;
    boolean stop=FALSE;
    char fmtval;
    char this_fmt[LINESIZE];
    long ival;
    va_list ap;
    double dval;
    long fmti=0;
    char cval, *sval;
    float y = *page_height;
    long ii=0;
    boolean left_aligned=FALSE;
    float offset = 55;
    va_start(ap, fmt);
    while (*fmt)
      {
        if(*fmt == '%')
          {
            if(width[ii] <= EPSILON)
              {
                left_aligned = TRUE;
                offset = 55 - width[ii++];
              }
            else
              {
                left_aligned = FALSE;
                offset = 55 + width[ii++];
              }
            
            start = TRUE;
            *fmt++;
            this_fmt[0] = '%';
            fmti = 1;
            continue;
          }
        else
          {   
            switch(*fmt) 
              {
                case 'c':                       /* string */
                    stop = TRUE;
                    fmtval = 'c';
                    cval = va_arg(ap, int);
                    this_fmt[fmti++] = fmtval;
                    this_fmt[fmti] = '\0';
                    if(left_aligned)
                        pdf_printf(offset, y, 'L', this_fmt, cval);
                    else
                        pdf_printf_ralign(offset, y,this_fmt, cval);
                    break;
                case 's':                       /* string */
                    stop = TRUE;
                    fmtval = 's';
                    sval = va_arg(ap, char *);
                    this_fmt[fmti++] = fmtval;
                    this_fmt[fmti] = '\0';              
                    if(left_aligned)
                        pdf_printf(offset, y,'L', this_fmt, sval);
                    else
                        pdf_printf_ralign(offset, y,this_fmt, sval);
                    break;
                case 'i':                       /* int */
                    stop = TRUE;
                    fmtval = 'i';
                    ival = va_arg(ap, int);
                    this_fmt[fmti++] = fmtval;
                    this_fmt[fmti] = '\0';   
                    if(left_aligned)
                        pdf_printf(offset, y,'L', this_fmt, ival);
                    else
                        pdf_printf_ralign(offset, y,this_fmt, ival);
                    break;
                case 'f':                       /* char */
                    stop = TRUE;
                    fmtval = 'f';
                    dval = va_arg(ap, double);
                    this_fmt[fmti++] = fmtval;
                    this_fmt[fmti] = '\0';                     
                    if(left_aligned)
                        pdf_printf(offset, y,'L', this_fmt, dval);
                    else
                        pdf_printf_ralign(offset, y,this_fmt, dval);
                    break;
                case '\n':
                    pdf_advance(&y);
                    break;
              }
            if(start)
              {
                if(!stop)
                  {
                    if(strchr("1234567890.", *fmt)!=NULL)
                        this_fmt[fmti++] = *fmt;
                    if(*fmt == '-')
                        left_aligned = TRUE;
                  }
                else
                  {
                    start = FALSE;
                    stop = FALSE;
                  }
              }
            *fmt++;
          }
      }
    va_end(ap);
}


void pdf_print_allelelegend(float *column_width, float *page_height, long loci)
{
    long locus;
    char stemp[LINESIZE];
    sprintf(stemp, "%-s", (loci == 1 ? "locus" : "loci "));
    
    pdf_printf(column_width[0], *page_height, 'L', "Indiv.");
    for(locus=1; locus < loci+1; locus++)
      {
        pdf_printf(column_width[locus], *page_height, 'L', "%li",locus);        
        if(column_width[locus+1]<column_width[locus])
          {
            pdf_advance(page_height);
          }
      }
}

#define NOTFIRSTDATAPAGE (boolean) 0
#define FIRSTDATAPAGE (boolean) 1
///
/// print data header for allelic or sequence data
void
pdf_print_dataheader (boolean first, char *title, float *page_height, float *column_width, long pop, world_fmt * world,
                      option_fmt * options, data_fmt * data)
{
    
    float w;
    float page_width;
    float left_margin = 55;
    
    // setup new page and title
    if(first)
      {
        pdf_new_page("");
        pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 18);    
        w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
        *page_height = pdf_contents_get_height(canvas) - 75;
        page_width = pdf_contents_get_width(canvas);
        *page_height -= 20;
        pdf_print_contents_at((page_width - w)/2, *page_height, title);
        pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
        *page_height -= 20;
        pdf_draw_line(50, *page_height, page_width-50, *page_height);
        *page_height -= 20;
        pdf_contents_set_font_and_size(canvas, "Helvetica", 10);    
      }
    else
      {
        pdf_advance(page_height);
        page_width = pdf_contents_get_width(canvas);
        pdf_draw_line(50, *page_height, page_width-50, *page_height);
        pdf_advance(page_height);    
      }
    pdf_printf_next(left_margin, page_height,"\n%-s", data->popnames[pop]);
    pdf_draw_line(50, *page_height, page_width-50, *page_height);
    pdf_advance(page_height);   
    if (!strchr (SEQUENCETYPES, options->datatype))
        pdf_print_allelelegend(column_width, page_height, world->loci);
    //    else
    //    pdf_print_seqlegend(left_margin, page_height, world->loci);
    
    pdf_advance(page_height);            
    pdf_draw_line(50, *page_height, page_width-50, *page_height);
    pdf_advance(page_height);   
}



///
/// print the allele typ data in diploid format 
/// \param page_width page width
/// \param page_height is manipulated and is always at the actual y coordinate on the page
/// \param[in] margin is the margin that is used, either left or right margin
/// \param[in] world contains all main data structures
/// \param[in] data  contains all the biological data structures
/// \param[in] options cintains all option structures
void pdf_print_alleledata (float *page_height, float margin, world_fmt * world, data_fmt * data, option_fmt * options)
{
  long top, ii;
    long pop, ind, locus;
    char stemp[LINESIZE];
    float w;
    float * column_width = (float *) mycalloc(world->loci+2, sizeof(float));
    float pixel;
    float total;
    float oldpage_height;
    float page_width = pdf_contents_get_width(canvas);
    
    // calculate column width for indvidual name and each locus
    // we use a cummulative sum to set the x-coord just right without further
    // calculations
    for (locus = 0; locus < data->loci; locus++)
      {
        for (pop = 0; pop < data->numpop; pop++)
          {
	    if(options->randomsubset > 0 && options->randomsubset < data->numind[pop][locus])
	      {
		top = options->randomsubset;
	      }
	    else
	      {
		top = data->numind[pop][locus];
	      }
	    
	    for (ii = 0; ii < top; ii++)
	      {
		ind = data->shuffled[pop][locus][ii];
		
                w = (float)( 2. + (float) (strlen (data->yy[pop][ind][locus][0]) +
                                           strlen (data->yy[pop][ind][locus][1])));
                if(column_width[locus+2] < w)
                    column_width[locus+2] = w;
              }
          }
      }
    // calculate columnwidth for in page units for individual name and locus column
    column_width[1] = 1. + (float) sprintf(stemp, "%-*.*s", (int) options->nmlength,
                                           (int) options->nmlength, data->indnames[0][0][0]);
    w = (float) pdf_contents_get_text_width(canvas, stemp, NULL, NULL);
    pixel = w / options->nmlength;
    column_width[0] = margin;
    column_width[1] = margin + w;
    total = margin + w ;
    for (locus = 0; locus < data->loci; locus++)
      {
        column_width[locus+2] *= pixel;
        total += column_width[locus+2];
        if(total > (page_width - 55))
          {
            total = margin;
            column_width[locus+2] = total;
          }
        else
          {
            column_width[locus+2] = total;
          }
      }
    
    for (pop = 0; pop < data->numpop; pop++)
      {
        if(pop==0)
            pdf_print_dataheader (FIRSTDATAPAGE, "Allelic data", page_height, column_width, pop, world, options, data);
        else
            pdf_print_dataheader (NOTFIRSTDATAPAGE, "Allelic data",page_height, column_width, pop, world, options, data);
        
        for (ind = 0; ind < data->numind[pop][0]; ind++)
          {
            pdf_printf(column_width[0], *page_height, 'L', "%-*.*s", (int) options->nmlength,
                       (int) options->nmlength, data->indnames[pop][ind][0]);
            for (locus = 0; locus < data->loci; locus++)
              {
                pdf_printf(column_width[locus+1], *page_height, 'L', " %s.%-s",
                           data->yy[pop][ind][locus][0],
                           data->yy[pop][ind][locus][1]);
                if(column_width[locus+2]<column_width[locus+1])
                  {
                    pdf_advance(page_height);
                  }
              }
            oldpage_height = *page_height;
            pdf_advance(page_height);
            if(oldpage_height < *page_height)
                pdf_print_dataheader (NOTFIRSTDATAPAGE, "Allelic data", page_height, column_width, pop, world, options, data);
            
          }
      }
    myfree(column_width);
}




void
pdf_print_data (world_fmt * world, option_fmt * options, data_fmt * data, float *orig_page_height, float *orig_left_margin)
{
    if (options->printdata)
      {
        switch (options->datatype)
          {
            case 'a':
            case 'b':
            case 'm':
                pdf_print_alleledata (orig_page_height,*orig_left_margin, world, data, options);
                break;
            case 's':
            case 'n':
	  case 'h':
            case 'u':
            case 'f':
                
                pdf_print_seqdata (orig_page_height, *orig_left_margin, world, data, options);
                break;
          }
      }
}



void pdf_print_sequence(float left_margin, float right_margin, float *page_height, data_fmt *data, long locus, long pop, long ind)
{
    long site;
    float w = 0.;
    float wtot = 0.;
    char stemp[LINESIZE];
    w = 60.;
    for(site=0; site < data->seq[0]->sites[locus]; site+=10)
      {
        sprintf(stemp,"%-10.10s", &data->yy[pop][ind][locus][0][site]);
        pdf_contents_set_font_and_size(canvas, "Courier", 9);
        pdf_print_contents_at(left_margin + wtot, *page_height, stemp);
        wtot += w;
        if((left_margin + wtot + w) > right_margin)
          {
            wtot = 0;
            pdf_advance(page_height);
          }            
      }
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
}

void pdf_print_seqdata (float *page_height, float margin, world_fmt * world, data_fmt * data, option_fmt * options)
{
  long top;
  long ii;
    long pop, ind, locus=0;
    char stemp[LINESIZE];
    float w;
    float * column_width = (float *) mycalloc(world->loci+2, sizeof(float));
 //xcode   float pixel;
    float oldpage_height;
    float page_width = pdf_contents_get_width(canvas);
    
    // calculate columnwidth for in page units for individual name and locus column
    column_width[1] = 1. + (float) sprintf(stemp, "%-*.*s", (int) options->nmlength,
                                           (int) options->nmlength, data->indnames[0][0][0]);
    w = (float) pdf_contents_get_text_width(canvas, stemp, NULL, NULL);
    //xcode pixel = w / options->nmlength;
    column_width[0] = margin;
    column_width[1] = margin + w;
    column_width[2] = margin + w + (float) 10;
    
    for (pop = 0; pop < data->numpop; pop++)
      {
        if(pop==0)
            pdf_print_dataheader (FIRSTDATAPAGE, "Sequence data", page_height, column_width, pop, world, options, data);
        else
            pdf_print_dataheader (NOTFIRSTDATAPAGE,  "Sequence data", page_height, column_width, pop, world, options, data);

	for (locus = 0; locus < data->loci; locus++)
	  {
	    
	    top = max_shuffled_individuals(options, data, pop, locus);
	    
	    for (ii = 0; ii < top; ii++)
	      {
		ind = data->shuffled[pop][locus][ii];
		pdf_printf(column_width[0], *page_height, 'L', "%-*.*s", (int) options->nmlength,
			   (int) options->nmlength, data->indnames[pop][ind][0]);
                pdf_printf(column_width[1], *page_height, 'L', "%li",locus+1);
                pdf_print_sequence(column_width[1]+10, page_width - 55, page_height, data, locus, pop, ind);
                pdf_advance(page_height);
              }
            oldpage_height = *page_height;
            pdf_advance(page_height);
            pdf_advance(page_height);
            if(oldpage_height < *page_height)
                pdf_print_dataheader (NOTFIRSTDATAPAGE,  "Sequence data", page_height, column_width, pop, world, options, data);            
          }
      }
    myfree(column_width);
}

//void format_helper_vec(atl_fmt *atl, long pop, long locus, long maxrep, long numpop, long *fmt_vector)
//{
//  for (i=0; i < atl.params; i++)
//}
    
void format_helper(MYREAL * matrix, long len, int *fmt1, int *fmt2)
{
    long i;
    long count = 0;
    long maxcount = 0;
    MYREAL val;
    
    for(i=0;i<len; i++)
      {
        val = matrix[i];
        if(val>0.)
            count = (long)(log10(val)) + 1; 
        if(count>maxcount)
            maxcount = count;     
      }
    *fmt1 = 7;
    if(maxcount<2)
      {
        *fmt2 = 5;
      }
    else
      {
        *fmt2 = 2;
      }
}


///
/// print matrix line for various uses
/// - print linear shortened migration matrix in n x n matrix with - as diagonals
/// alowing for some pretty print using a format_helper() the columns are right-aligned
/// based on the width.
void pdf_print_matrix_line(long whichline, float *page_height, float left_margin, float width, long cols,long rows,MYREAL * matrix,boolean ismig)
{
    long i=0;
    long col, row = whichline;
    float offset;
    float page_width = pdf_contents_get_width(canvas);
    int fmt1=5, fmt2=3;
    long pad = ((ismig)?cols:0);
    format_helper(matrix+pad, cols * rows - pad ,&fmt1, &fmt2);
    offset = left_margin;
    if(ismig)
        i = cols;
    for(col = 0; col < cols; col++)
      {
        if(offset > page_width - 55)
          {
            offset = left_margin;
            pdf_advance(page_height);
          }
        if(ismig && row==col)
          {
            pdf_printf_ralign(offset+width - 15, *page_height, "-");
          }
        else
          {
            pdf_printf_ralign(offset + width - 5, *page_height, "%*.*f", fmt1, fmt2, matrix[i]);
            i++;
          }
        offset += width;
      }
    pdf_advance(page_height);
}

///
/// print text using a simple formatter that works on a single paramgraph
/// breaking it approx at the end of the line.
void pdf_print_comment(float lx, float *ly, char *this_text)
{
    pdf_printf(lx, *ly,'L', this_text);
    pdf_advance(ly);
}



/*
 * print header
 * ===========================================================================
 * === Titletext
 * ===========================================================================
 * === Population     Loci  Ln(L)   Theta    xNm [xNe mu] xx     xx     xx
 * xx     xx     xx xx.... -------------- ---- -------- --------
 * ----------------------------------------
 */
/// print the MLE table header
void
pdf_print_result_header (float *lx, float *page_height, char *titletext, world_fmt * world)
{
    long p1, zz;
    float *ly = page_height;
    float page_width;
    
    pdf_print_section_title(&page_width, ly, titletext);
    pdf_advance(page_height);
    pdf_printf(lx[0], *ly,'L',"Population [x]");
    pdf_printf(lx[1], *ly,'L',"Loc.");
    pdf_printf(lx[2], *ly,'L',"Ln(L/L0)");
    symbol_Theta(lx[3], *ly,12,-1);
    //    pdf_printf(lx[3], *ly,'L',"Theta");
    if(world->numpop>1)
      {
        if(world->options->usem)
            pdf_printf(lx[4], *ly,'L',"M (m/mu) [+=receiving population");
        else
            pdf_printf(lx[4], *ly,'L',"xNm [+=receiving population");
      }
    pdf_advance(ly);
    pdf_printf(lx[3], *ly,'L',"[x Ne mu]");
    
    zz = 4;
    for (p1 = 0; p1 < world->numpop; p1++)
      {
        if (zz > 8)
          {
            zz = 4;
            pdf_advance(ly);
          }
        pdf_printf(lx[zz], *ly,'L',"%2li,+", p1 + 1);
        zz++;
      }
    pdf_advance(ly);
    pdf_draw_line(lx[0],*ly, page_width - 55, *ly);
    pdf_advance(ly);
} 

/// print the string for a population
void pdf_print_popstring(float *lx, float *page_height, long pop, world_fmt *world, option_fmt *options, data_fmt *data)
{
    float *ly = page_height;
    //    char popstring[LINESIZE];
    if (options->readsum)
      {
        pdf_printf(lx[0], *ly, 'L', "%2li:",pop+1);
      }
    else
      {
        pdf_printf(lx[0],*ly, 'L', "%2li:%10.10s",pop+1, data->popnames[pop]);
      }
}


/// print the replicate number
void pdf_print_replicate(float lx, float *page_height, world_fmt *world, long maxrep, long rep, long locus)
{
    float *ly = page_height;
    char repstring[LINESIZE];
    sprintf (repstring, "%li", rep + 1);
    pdf_printf(lx, *ly, 'L', "%li %s", locus + 1, maxrep > 1 ? (rep == maxrep - 1 ? " A" : repstring) : "  ");
}

/// print the MLE table content for each population
void
pdf_print_result_population (float *lx, float *page_height, long pop, world_fmt * world,
                             option_fmt * options, data_fmt * data)
{
    float *ly = page_height;
    long skipped = 0, locus;
    long maxrep = world->options->replicate ?
        (world->options->replicatenum > 0 ?
         world->options->replicatenum + 1 : world->options->lchains + 1) : 1;
    long rep;
    pdf_print_popstring(lx, page_height, pop, world, options, data);
    for (locus = 0; locus < world->loci; locus++)
      {
        if (world->data->skiploci[locus])
          {
            skipped++;
            continue;
          }
        for (rep = 0; rep < maxrep; rep++)
          {
            pdf_print_replicate(lx[1], page_height, world, maxrep, rep, locus);
            pdf_printf(lx[2], *ly, 'L', "% 8.3f ",
                       world->atl[rep][locus].param_like);
            pdf_print_result_param (lx, page_height,  world->atl[rep][locus].param,
                                    world->numpop, pop, world->options->usem);
          }
      }
    if (world->loci - skipped > 1)
      {
        pdf_printf(lx[1], *page_height, 'L', "All ");
        //locus is exactly world->loci
        // re is always one because we have only replication of single locus chains
        pdf_printf(lx[2],*page_height, 'L', "% 8.3f ", world->atl[0][locus].param_like);
        pdf_print_result_param (lx, page_height,  world->atl[0][locus].param,
                                world->numpop, pop, world->options->usem);
      }
    /* FPRINTF(world->outfile,"%s\n",sline);     */
}


/// print the parameter for the MLE parameter printout
void
pdf_print_result_param (float *lx, float *page_height,  MYREAL *param, long numpop, long pop,
                        boolean usem)
{
    char    temp[LINESIZE];
    long    i;
    long    zz;
    long    msta = mstart (pop, numpop);
    long    msto = mend (pop, numpop);
    float  *ly   = page_height;
    //int   fmtint = 8;
    //int fmtfloat = 5;
    //int digits;
    MYREAL  tmp  = 0;
    
    // population size
    if (param[pop] <= SICK_VALUE)
        pdf_printf (lx[3], *ly, 'L', "-");
    else
      {
        nice_element(param[pop], temp, 0.0001, 10, 100, 4, 2, '\0');
        pdf_printf (lx[3], *ly, 'L', "%s", temp);
        /*        if (param[pop] < 0.0001 && param[pop] > 0)
            pdf_printf (lx[3], *ly, 'L', "%3.2e ", param[pop]);
        else
          {
            digits = (int) log10(param[pop]);
            if(digits>3)
              {
                fmtfloat=0;
                fmtint=digits;
              }
            pdf_printf (lx[3], *ly, 'L', "%*.*f ", fmtint, fmtfloat, param[pop]);
          }
        */
      }
    
    // migration rate
    zz=4;
    for (i = msta; i < msto; i++)
      {
        if (zz > 8)
          {
            zz = 4;
            pdf_advance(ly);
          }
        if (pop == i - msta)
          {
            pdf_printf (lx[zz], *ly, 'L', "-");
            zz++;;
          }
        if ((param[i] <= SICK_VALUE) || (param[pop] <= SICK_VALUE))
            pdf_printf (lx[zz], *ly, 'L', "-");
        else
          {
            if (usem)
              {
                //tmp = param[i];
                nice_element(param[i], temp, 0.001, 100, 1000, 3, 2, '\0');
                //if (tmp < 0.0001)
                pdf_printf (lx[zz], *ly, 'L', "%s",temp);
                //else
                //pdf_printf (lx[zz], *ly, 'L', "%7.4f ", tmp);
              }
            else
              {
                tmp = param[pop] * param[i];
                nice_element(tmp, temp, 0.0001, 10, 100, 4, 2, '\0');
                pdf_printf (lx[zz], *ly, 'L', "%s",temp);
                /*		if (tmp < 0.00001)
                    pdf_printf (lx[zz], *ly, 'L', " 0.0000 ");
                else
                    pdf_printf (lx[zz], *ly, 'L', "%7.4f ", tmp);
                */
              }
          }
        zz++;
      }
    if (pop == numpop - 1)
        pdf_printf (lx[i-msta+4], *ly, 'L',"-");
    pdf_advance(ly);
}



void 
pdf_print_results (float *page_height, world_fmt ** universe, option_fmt * options, data_fmt * data)
{
    float lx[]={55,140,170,220,270,320,370,420,470,520};
    
    long pop;
    float *ly = page_height;
//xcode    FILE *outfile;
    world_fmt *world = universe[0];
    worldoption_fmt *wopt = world->options;
    char sch[10], lch[10], cva[50];
    long rep = world->loci > 1 ? 0 : (wopt->replicate ? world->repstop : 0);
    //xcode outfile = world->outfile;
    if (options->schains == 1)
        strcpy (sch, "chain");
    else
        strcpy (sch, "chains");
    if (options->lchains == 1)
        strcpy (lch, "chain");
    else
        strcpy (lch, "chains");
    pdf_print_result_header (lx, page_height, "Maximum Likelihood estimates", world);
    for (pop = 0; pop < world->numpop; pop++)
      {
        pdf_print_result_population (lx, page_height, pop, world, options, data);
      }
    pdf_advance(ly);
    pdf_printf(lx[0],*ly,'L', "Comments:");
    pdf_advance(ly);
    pdf_printf(lx[0],*ly,'L', "The x is 1, 2, or 4 for mtDNA, haploid, or diploid data, respectively");
    pdf_advance(ly);
    pdf_printf(lx[0],*ly,'L', "There were %li short %s (%li used trees out of sampled %li)",
               options->schains, sch, options->ssteps, options->sincrement * options->ssteps);
    pdf_advance(ly);
    pdf_printf(lx[0],*ly, 'L', "and %li long %s (%li used trees out of sampled %li)\n",
               options->lchains, lch, options->lsteps,
               options->lincrement * options->lsteps);
    pdf_advance(ly);
    if (wopt->heating)
      {
        if(options->adaptiveheat!=NOTADAPTIVE)
          {
	    if(options->adaptiveheat==STANDARD)
	      {
		pdf_printf(lx[0],*ly,'L',
			   "Adaptive heating with %li chains was active",options->heated_chains);
	      }
	    else
	      {
		pdf_printf(lx[0],*ly,'L',
			   "Bounded adaptive heating with %li chains was active",options->heated_chains);
	      }
            pdf_advance(ly);
            pdf_printf(lx[0],*ly,'L',"check Log file (if present) for average temperatures");
            
            pdf_printf(lx[0],*ly,'L',"Average last chain temp: 1.0",
                       options->heated_chains);
            for(pop=1;pop<options->heated_chains;pop++)
                pdf_printf(lx[pop+2],*ly,'L', ", %f", universe[pop]->averageheat);
            pdf_advance(ly);
          }
        else
            pdf_printf(lx[pop+2],*ly,'L', "Static heating with %li chains was active\n",options->heated_chains);
      }
    if (options->gamma)
      {
        if (world->atl[rep][world->loci].param[world->numpop2] < 10e-9)
            strcpy (cva, "0");
        else
            sprintf (cva, "%f",
                     sqrt (1. / world->atl[rep][world->loci].param[world->numpop2]));
        pdf_printf(lx[0],*ly,'L',"With shape parameter Alpha=%g ([1/CV(mu)]^2; CV(mu)=%s)",
                   world->atl[rep][world->loci].param[world->numpop2],
                   cva);
        pdf_advance(ly);
      }
    if (world->options->replicate)
      {
        if (world->repkind == MULTIPLECHAIN)
            pdf_printf(lx[0],*ly, 'L', "  COMBINATION OF ALL LONG CHAINS");
        else
            pdf_printf(lx[0],*ly, 'L', "  COMBINATION OF %li MULTIPLE RUNS",
                       world->options->replicatenum);
        pdf_advance(ly);
      }
    if (world->atl[rep][world->loci].normd > LOCI_NORM)
      {
        pdf_printf(lx[0],*ly,'L',"  [Last maximization needed %li cycles of maximal %i,",
                   world->atl[rep][world->loci].trials,
                   NTRIALS);
        pdf_advance(ly);
        pdf_printf(lx[0],*ly,'L',"   Norm(first derivatives)=%f (Normal stopping criteria is < %f)]", 
                   world->atl[rep][world->loci].normd, LOCI_NORM);
      }
    pdf_advance(ly);
    pdf_advance(ly);
}

void prepare_profile_header(long which, char method, char *header, world_fmt *world)
{
    const long numpop = world->numpop;
    
    long pop;
    long pop2;
    long from;
    long to;
//xcode     long z;
    long position = 0;
    char *mm;
    
    mm = (char *) mycalloc(LINESIZE,sizeof(char));
    if (world->options->profileparamtype == PLOT4NM)
        sprintf(mm,"%s","xNm");
    else
        sprintf(mm,"%s","M");
    
    if(method=='p')
      {
        position += sprintf(header + position, "Per. &");
      }
    else
      {
        position += sprintf(header + position," &");
      }
    position += sprintf(header + position,"Ln(L)&");
    if(which < numpop)
        position += sprintf(header + position,"%s_%li&","Q",which+1);
    else
      {
        m2mm(which,numpop,&from,&to);
        position += sprintf(header + position,"%s_%li->%li&",mm,from+1,to+1);
      }
    for(pop=0; pop < numpop; pop++)
      {
        position += sprintf(header + position,"%s_%li&","Q",pop+1);
      }
    //xcode z=pop+3;
    for(pop=0; pop < numpop; pop++)
      {
        for(pop2=0; pop2 < numpop; pop2++)
          {
            if(pop!=pop2)
              {
                position += sprintf(header + position,"%s_%li->%li&",mm,pop2+1,pop+1);
              }
          }
      }
    //xcode position += 
    sprintf(header + position,"%%&");
    myfree(mm);
}

///
/// pretty print for tables
/// takes a lower and upper bound below and above it uses scientific notation, the range in
/// between prints fixed point notation below mid with low_mid_digits and above mid with mid_upper_digits. 
long nice_element(MYREAL param, char *element, MYREAL lower, MYREAL mid, MYREAL upper, 
                  int low_mid_digits, int mid_upper_digits, char delimiter)
{
    long position = 0;
    if (param < lower)
        position = sprintf(element,"%.2e", param);
    else
      {
        if (param > upper)
            position = sprintf(element,"%.2e", param);
        else
          {
            if(param > mid)
                position = sprintf(element,"%.*f", mid_upper_digits,param);
            else
                position = sprintf(element,"%.*f", low_mid_digits,param);
          }
      }
    if(delimiter != '\0')
        position += sprintf(element + position,"%c",delimiter);
    return position;
}


/// 
/// prepare the table elements for the profile likelihood table
long
prepare_profile_elements (char method, long which,  MYREAL likes[],
                          MYREAL **param, world_fmt * world, char **elements)
{
    const MYREAL probabilities[] = SHOWGRID;
    MYREAL       likemax         = -MYREAL_MAX;  
    long         i;
    long         ii;
    long         j;
    long         elem     = world->options->gamma ? (world->numpop2 + 1) : world->numpop2;
    long         position = 0L;
    long         likei    = 0L;
    long         failed   = 0L;
    char         likestr[LINESIZE];
    
    for (i = 0; i < GRIDSIZE; i++)
      {
        if (likemax < MIN(likes[i], MYREAL_MAX))
          {
            likemax = likes[i];
            likei   = i;
          }
      }
    
    for(ii=0; ii< GRIDSIZE; ii++)
      {
        i = warp (ii);
        //printf("warped ii=%li --> i=%li\n",ii,i);
        position = 0;
        //first column
        if(method != 'd')
          {
            if(world->percentile_failed[which][i])
              {
                position += sprintf(elements[i] + position,"****&");
                failed++;
              }
            else
              {
                if (probabilities[i]==0.5)
                    position += sprintf(elements[i] + position,"MLE&");
                else
                    position += sprintf(elements[i] + position,"%4.3f&", probabilities[i]);
              }
          }
        else
          {
            position += sprintf(elements[i] + position,"&");
          }
        
        //second column
        if (likes[i] >= MYREAL_MAX - EPSILON)
            sprintf (likestr, "-");
        else
            sprintf (likestr, "%.3f%c", likes[i], likei == i ? '*' : ' ');
        position += sprintf(elements[i] + position,"%s&", likestr);
        // third column
        position += sprintf(elements[i] + position,"%.6g&", param[i][which]);
        // fourth and further columns
        for (j = 0; j < world->numpop; j++)
          {
            position += nice_element(param[i][j],elements[i] + position, 0.0001, 10, 100, 4, 2, '&');
          }
        for (j = world->numpop; j < elem; j++)
          {
            position += nice_element(param[i][j],elements[i] + position, 0.001, 100, 1000, 3, 1, '&');
          }
        if(ii<GRIDSIZE-1)
          {            
            //xcode position += 
            sprintf(elements[i] + position,"@&");
          }
        else
          {            
            //xcode position += 
            sprintf(elements[i] + position,"%%&");
          }
      }
    return failed;
}

///
/// Printing PROFILE table title for PDF
/// \param world        master parameter structure
void pdf_print_profile_title(world_fmt *world)
{
    char *title = "Profile likelihood tables and plots";
    float w;
    float page_width;
    
    // setup new page and title
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 18);    
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    page_height = pdf_contents_get_height(canvas) - 75;
    page_width = pdf_contents_get_width(canvas);
    page_height -= 20;
    pdf_print_contents_at((page_width - w)/2, page_height, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    page_height -= 20;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    page_height -= 20;
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);    
}

void method_set(float lx, float *page_height, char method)
{
    switch (method)
      {
        case 'p':
            pdf_printf(lx, *page_height, 'L', "Parameters are evaluated at percentiles using bisection method (slow, but exact).");
            break;
        case 's':
            pdf_printf(lx, *page_height, 'L',"Parameters are evaluated at percentiles");
            pdf_advance(page_height);
            pdf_printf(lx, *page_height, 'L',"using cubic splines of profiled parameter (faster, but not so exact).");
            break;
        case 'd':
            pdf_printf(lx, *page_height, 'L',"Parameters are evaluated at pre-defined values\n");
            break;
        case 'q':
            pdf_printf(lx, *page_height, 'L', "Parameters are evaluated assuming complete independence\n");
            break;
        case 'f':
            pdf_printf(lx, *page_height, 'L', "Parameters are evaluated assuming complete independence");
            pdf_advance(page_height);
            pdf_printf(lx, *page_height, 'L',"and then once maximized at the found profiled parameter value");
            break;
      }
    pdf_advance(page_height);
}

void pdf_table_footnote(float lx, float * page_height, long failed, boolean percentiles)
{
    char * foot1 = "If the percentile column is marked with **** then the convergence to the percentile value failed";
    char * foot2 = "The values are still CORRECT but not at the percentile of the profile parameter";
    char * foot3 = "Values with a * are NOT at the percentiles!";
    char * foot4 = "The convergence to the percentile value FAILED.";
    if(failed>0)
      {
        pdf_printf(lx, *page_height, 'L', percentiles ? foot3 : foot1);
        pdf_advance(page_height);
        pdf_printf(lx, *page_height, 'L', percentiles ? foot4 : foot2);
        pdf_advance(page_height);
      }
}

///
/// Prepare PROFILE table for PDF and ASCII printing
/// \param method       profile method type ('p'=percentiles)
/// \param which        which parameter is currently profiled 
/// \param likes        likelihoods
/// \param param        parameters
/// \param world        master parameter structure
/// \return failed      how many of the profile calculations failed
long prepare_profile_table(char method, long which, MYREAL *likes, MYREAL **param, world_fmt *world, long bufsize)
{
    boolean failed;
    //  float left_margin = 55;
    long i;
    long ii;
    //  long pop;
    char *header;
    char **elements;
    //  long elem = world->options->gamma ? world->numpop2 + 1 : world->numpop2;
    long position = bufsize;
    
    world->buffer = (char *) myrealloc (world->buffer, (position + SUPERLINESIZE + 1) * sizeof (char));
    memset(world->buffer + position, 0 , (SUPERLINESIZE+1)* sizeof(char));
    
    header = (char *) mycalloc(LINESIZE, sizeof(char));
    charvec2d(&(elements), 9, LINESIZE);
    
    //  for(pop=0;pop < world->numpop2; pop++)
    //  {
    prepare_profile_header(which, method, header, world);
    position += sprintf(world->buffer + position,"%s",header);
    failed = prepare_profile_elements(method, which, likes, param, world, elements);
    for(i=0;i<GRIDSIZE;i++)
      {
        ii = warp(i);
        position += sprintf(world->buffer + position,"%s",elements[ii]);
      }
    //}
    myfree(header);
free_charvec2d(elements);
return failed;
}


///
/// translate buffer table (syntax similar to LaTeX) into table header and table elements
void  translate_buffer_table(long cols, long rows, char **thebuffer, char **header, char ***elements)
{
    char *temp;
    char *buffer; 
    char *savebuffer;
    long z = 0;
    long r = 0;
    //  long position=0;
    
    buffer = (char *) mycalloc(strlen(*thebuffer)+1,sizeof(char));
    savebuffer= buffer;
    strcpy(buffer,*thebuffer);
    while((temp=strsep(&buffer,"&"))!=NULL && buffer!=NULL )
      {
        if(temp[0]=='%')
            break;
        else
          {
            sprintf(header[z++],"%s",temp);
          }
      }
    z=0;
    r=0;
    while((temp=strsep(&buffer,"&"))!=NULL && buffer!=NULL)
      {
        if(r==rows-1 && z>=cols)
          {
            //	  position = (long) strlen(strstr(thebuffer,buffer))+1;
            // memcpy(*thebuffer,buffer,(strlen(buffer)+1)*sizeof(char));
            break;
          }
        if(temp[0]=='%')
            break;
        if(temp[0]=='@')
          {
            z=0;
            r++;
          }
        else
          {
            if(z<cols)
              {
                sprintf(elements[r][z],"%s",temp);
                z++;
              }
            else
              {
                if(r<rows)
                    warning("column counter exceeded available columns: z=%li temp=%s\n",z, temp);
                else
                    break;
              }
          }
      }
    myfree(savebuffer);
    //  return position;
}
///
/// extract column out of  buffer table (syntax similar to LaTeX) into table header and table elements
void  extract_column_buffer_table(long col, long cols, long rows, char **thebuffer, float *x)
{
    char *temp;
    char *buffer; 
    char *savebuffer;
    long z = 0;
    long r = 0;
    //  long position=0;
    
    buffer = (char *) mycalloc(strlen(*thebuffer)+1,sizeof(char));
    savebuffer= buffer;
    strcpy(buffer,*thebuffer);
    while((temp=strsep(&buffer,"&"))!=NULL && buffer!=NULL )
      {
        if(temp[0]=='%')
            break;
        else
          {
	    z++;
	    //            sprintf(header[z++],"%s",temp);
          }
      }
    z=0;
    r=0;
    while((temp=strsep(&buffer,"&"))!=NULL && buffer!=NULL)
      {
        if(r==rows-1 && z>=cols)
          {
            //	  position = (long) strlen(strstr(thebuffer,buffer))+1;
            // memcpy(*thebuffer,buffer,(strlen(buffer)+1)*sizeof(char));
            break;
          }
        if(temp[0]=='%')
            break;
        if(temp[0]=='@')
          {
            z=0;
            r++;
          }
        else
          {
	    if(z==col)
	      x[r] = atof(temp);

            if(z<cols)
              {
                //sprintf(elements[r][z],"%s",temp);
                z++;
              }
            else
              {
                if(r<rows)
                    warning("column counter exceeded available columns: z=%li temp=%s\n",z, temp);
                else
                    break;
              }
          }
      }
    myfree(savebuffer);
    //  return position;
}

///
/// Printing PROFILE table for PDF
/// \param left_margin  left start x coord
/// \param *page_height start y coordinate
/// \param page_width   width of paper
/// \param world        master parameter structure
/// \param data         data structures
/// \param options      options
///
void pdf_print_profile_table(float left_margin, float * page_height, char method, char *buffer, world_fmt *world)
{
    char *thisparam;
    char *savedparam;
    char *ptr;
    //  long position = 0;
    
    boolean failed = world->percentile_some_failed;
    float page_width;
    float red[3]={0.99,0.,0.};
    float *x = (float *) calloc(9,sizeof(float));
    float *y = (float *) calloc(9,sizeof(float));

    thisparam = (char *) mycalloc(LINESIZE,sizeof(char));
    savedparam = thisparam;
    strncpy(thisparam, buffer, LINESIZE-1);
    strsep(&thisparam,"&"); // remove first element
    strsep(&thisparam,"&"); // remove second  element
    ptr = strsep(&thisparam,"&"); // keep third 
    page_width = pdf_contents_get_width(canvas);
    pdf_advance(page_height);
    pdf_draw_line(55,*page_height+5.0, page_width - 50, *page_height+5.0);
    pdf_advance(page_height);
    pdf_printf(left_margin, *page_height,'L',"Profile likelihood table and plot for parameter ");
    pdf_print_symbol(left_margin+205, page_height, 12, 'L', ptr);
    pdf_advance(page_height);
    method_set(left_margin, page_height, method);
    
    pdf_table(left_margin, page_height,world->numpop2 + 3, 9, &buffer/* + position*/, 
              NULL, 4, 6.);
    pdf_table_footnote(left_margin, page_height, failed, NOT_PERCENTILES);
   
    extract_column_buffer_table(2, world->numpop2 + 3, 9, &buffer, x);
    extract_column_buffer_table(1, world->numpop2 + 3, 9, &buffer, y);
    pdf_page_advance_or_not(page_height, 100.);
    pdf_linedotplot(9, x, y, red, red, 1, page_height, page_width - 200., 100., TRUE, TRUE); 
    myfree(savedparam);
}

// generic table generator package
// containing functions:
// pdf_table()
// find_col_width()
// define_col_start()
// pdf_print_table_header()

///
/// finds the width of all columns given the elements of the table
void  find_col_width(int cols, int rows, char ***elements, char **header, float *col_widths)
{
    float w;
    float keepw;
    int row;
    int col;
    
    for(col=0;col < cols; col++)
      {
        keepw = (float) pdf_contents_get_text_width(canvas, header[col], NULL, NULL);
        for(row=0;row < rows; row++)
          {
            w = (float) pdf_contents_get_text_width(canvas, elements[row][col], NULL, NULL);
            if(w > keepw)
                keepw = w;
          }
        col_widths[col] = keepw;
      }
}

///
/// align column
float align_column(char position, float cw, float left_margin)
{
    float loc = 0.;
    switch(position)
      {
        case 'R':
            loc = left_margin + cw;
            break;
        case 'C':
            loc = left_margin + cw/2.;
            break;
        case 'L':
        default:
            loc = left_margin;
            break;
      }
    return loc;
}

///
/// sets X-coordinates to start the columns 
void  define_col_start(float cols, float * col_widths, int col_overflow, char *position, float left_margin, float page_width, float separator, float *col_starts)
{
    int col;
    char pos;
    if(position==NULL)
        pos='L';
    else
        pos=position[0];
    col_starts[0] = align_column(pos, col_widths[0], left_margin);
    for(col=1; col < cols; col++)
      {
        if(position==NULL)
            pos='L';
        else
            pos=position[col];
        col_starts[col] =  align_column(pos, col_widths[col], separator + col_starts[col-1] + col_widths[col-1]);
        if((col_starts[col] + col_widths[col]) > page_width) //correct for 'L', but questionable for 'R'
          {
            if(col_overflow==0)
                col_starts[col] = col_starts[0];
            else
              {
                if(col_overflow <= col)
                    col_starts[col] = col_starts[col_overflow];
                else
                    col_starts[col] = col_starts[0];
              }
            //	  fprintf(stdout,">>> %f %f %f\n",col_starts[col],col_widths[col], separator);
          }
      }
}  

///
/// prints the table header at column positions
void  pdf_print_table_header(float *page_height, char * position, int cols, float *col_starts,char **header)
{
    int col;
    char pos;
    long oldcolstarts = -1;
    for(col=0; col < cols; col++)
      {
        if(position==NULL)
            pos='L';
        else
            pos=position[col];
        if(col_starts[col] < oldcolstarts)
          {
            pdf_advance(page_height);
          }
	pdf_print_symbol(col_starts[col], page_height, 10, pos, header[col]);
        oldcolstarts = col_starts[col];
      }
}

void pdf_print_symbol(float lx, float *page_height, long fontsize, char pos, char *symbolstring)
{
  int topop;
  int frompop;
  switch(symbolstring[0])
    {
    case 'T':
      sscanf(symbolstring,"Theta_%i",&topop);
      symbol_Theta(lx, *page_height, fontsize, topop);
      break;
    case 'Q':
      sscanf(symbolstring,"Q_%i",&topop);
      symbol_Theta(lx, *page_height, fontsize, topop);
      break;
    case 'M':
      if(symbolstring[1]!='_')
	pdf_printf(lx, *page_height, pos, "%s", symbolstring);
      else
	{
	  sscanf(symbolstring,"M_%i->%i",&frompop,&topop);
	  symbol_M(lx, *page_height, fontsize, frompop, topop, TRUE);
	}
      break;
    case 'x':
      sscanf(symbolstring,"xNm_%i->%i",&frompop,&topop);
      symbol_M(lx, *page_height, fontsize, frompop, topop, TRUE);
      break;
    default:
      pdf_printf(lx, *page_height, pos, "%s", symbolstring);
    }
}

///
/// generic table generator using header and elements arrays
/// \param float left_margin left edge of table
/// \param float * page_height page height 
/// \param int cols number of columns in table if there are too many columns then new line
/// \param int rows number of rows in table, if a new page is needed then header is repeated
/// \param char *** elements all table elements, currently no formatting of these
/// \param char ** header   header row
/// \param char *position position of all elements: L=left, R=right, C=center
/// \param int col_overflow when there are too many columns this is the column to restart
long pdf_table2(float left_margin, float *page_height, int cols, int rows, char **header, char ***elements, char *position, int col_overflow, float separator)
{
    boolean new_page=FALSE;
    
    int row;
    int col;
    long location=0;
    
    char pos;
    
    float *col_widths;
    float page_width = pdf_contents_get_width(canvas);
    float right_margin = page_width - 55;
    float *col_starts;
    float oldcol = HUGE;
    
    col_widths = (float *) mycalloc(cols,sizeof(float));
    col_starts = (float *) mycalloc(cols,sizeof(float));
    
    find_col_width(cols, rows, elements, header, col_widths);
    define_col_start(cols, col_widths, col_overflow, position, left_margin, right_margin, separator, col_starts);
    
    pdf_print_table_header(page_height, position, cols, col_starts,header);
    new_page = pdf_advance(page_height);
    pdf_draw_line(55,*page_height+5.0,right_margin, *page_height+5.0);
    
    for(row=0; row< rows; row++)
      {
        if(new_page)
          {
            pdf_print_table_header(page_height, position, cols,col_starts,header);
            //xcode new_page = 
            pdf_advance(page_height);
            pdf_draw_line(55.,*page_height,right_margin, *page_height);
            new_page = pdf_advance(page_height);
          }
        for(col=0;col < cols; col++)
          {
            if(col_starts[col] < oldcol)
              {
                new_page = pdf_advance(page_height);
                if(new_page)
                  {
                    pdf_print_table_header(page_height, position, cols,col_starts,header);
                    //xcode new_page = 
                    pdf_advance(page_height);
                    pdf_draw_line(55.,*page_height,right_margin, *page_height);
                    new_page = pdf_advance(page_height);
                  }
              }
            if(position==NULL)
                pos='L';
            else
                pos=position[col];
            pdf_print_symbol(col_starts[col],page_height,10,pos,elements[row][col]);
	    //            pdf_printf(col_starts[col],*page_height, pos, "%s",elements[row][col]);
            oldcol = col_starts[col];
          }
      }
    pdf_advance(page_height);
    myfree(col_widths);
    myfree(col_starts);
    return location;
}

///
/// generic table generator using a buffer that contains header and elements, these
/// are then split into header and elements, this approach allows to use the existing
/// machinery for MPI.
/// \param float left_margin left edge of table
/// \param float * page_height page height 
/// \param int cols number of columns in table if there are too many columns then new line
/// \param int rows number of rows in table, if a new page is needed then header is repeated
/// \param char *** elements all table elements, currently no formatting of these
/// \param char ** header   header row
/// \param char *position position of all elements: L=left, R=right, C=center
/// \param int col_overflow when there are too many columns this is the column to restart
void pdf_table(float left_margin, float *page_height, int cols, int rows, char **buffer, char *position, int col_overflow, float separator)
{
    int row;
    //  long location=0;
    char **header;
    char ***elements;
    
    charvec2d(&header, cols, LINESIZE);
    elements = (char ***) mycalloc(rows, sizeof(char **));
    for(row=0; row < rows; row++)
      {
        charvec2d(&(elements[row]),cols, LINESIZE);
      }
    
    translate_buffer_table(cols, rows, buffer, header, elements);
    
    pdf_table2(left_margin, page_height, cols, rows, header, elements, position, col_overflow, separator);
    
    free_charvec2d(header);
    for(row=0; row < rows; row++)
        free_charvec2d(elements[row]);
    myfree(elements);
    // return location;
}


///
/// pdf print profile tables for percentiles, returns the number
/// of failed convergences to a specific percentile.
void
pdf_print_profile_percentile (world_fmt * world)
{
    long i, j, jj;
    float page_width = pdf_contents_get_width(canvas);
    char **buffer = &world->buffer;
    long failed=0;
    long cols = 10;
    long maxrows = world->numpop2 + (long) world->options->gamma;
    long rows=0;
    long row;
    char **header;
    char ***elements;
    charvec2d(&header, cols, LINESIZE);
    sprintf(header[0]," ");
    sprintf(header[1],"0.005");
    sprintf(header[2],"0.025");
    sprintf(header[3],"0.05");
    sprintf(header[4],"0.25");
    sprintf(header[5],"MLE");
    sprintf(header[6],"0.75");
    sprintf(header[7],"0.95");
    sprintf(header[8],"0.975");
    sprintf(header[9],"0.995");
    
    elements = (char ***) mycalloc(maxrows, sizeof(char **));
    for(row=0; row < maxrows; row++)
      {
        charvec2d(&(elements[row]),cols, LINESIZE);
      }  
    (*buffer)[0]='\0';
    pdf_print_section_title(&page_width, &page_height, "Summary of profile likelihood percentiles of all parameters");
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    pdf_printf(55.,page_height,'L',"Parameter");
    pdf_printf(55.+100.,page_height,'L',"Percentiles");
    pdf_advance(&page_height);
    pdf_draw_line(100., page_height, page_width-55., page_height);
    pdf_advance(&page_height);
    rows=0;
    for (i = 0; i < maxrows; i++)
      {
        if (world->quantiles[i].name[0] == '\0')
            continue;  /* this variable was not estimated, so return */
        
        sprintf (elements[rows][0], "%-10.10s", world->quantiles[i].name);
        for (j = 0; j < GRIDSIZE; j++)
          {
            jj = warp (j);
            if (world->quantiles[i].param[j] < 0)
              {
                sprintf (elements[rows][j], "-");
              }
            else
              {
                nice_element(world->quantiles[i].param[jj],elements[rows][j+1], 0.0001, 100, 1000, 4 , 2, '\0');
              }
          }
        rows++;
      }
    
    pdf_table2(55., &page_height, cols, rows, header, elements, NULL, 1, 10.);
    pdf_table_footnote(55., &page_height, failed, PERCENTILES);
    
    myfree(header);
    for(row=0; row < maxrows; row++)
        free_charvec2d(elements[row]);
    myfree(elements);
    
}




///
/// plot event histograms
void
pdf_event_histogram(long loci, long numparams,  world_fmt *world)
{
    float   page_width;
    float   lx;
    float   ly;
    float   ph;
    float   delta;    
    float   w;
    float * binning;
    
    long    eventbinnum_allmax = 0L;
    long    loc;
    long    i;
    long    i0;
    long    j;
    long    bins;
    long    frompop;
    long    topop;
    long  * eventbinnum;
    
    char  * set50;
    char  * set95;
    char    title[LINESIZE];
    float total;
    
    MYREAL ** eventbins_all;
    
    duo ** eventbins;
    
    
    if (loci > 1)
        sprintf(title,"Summary of events through time over all loci");
    else
        sprintf(title,"Events through time");
    
    // add a new page so that we can print at least four histograms
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    
    /* Start to print text. */ 
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    // first histogram position
    page_height -= 150;
    lx = 100;
    ly = page_height;
    
    //    for (loc = 0; loc < loci; loc++)
    //  {
    loc = loci > 1 ? loci : 0;
    eventbinnum = world->mighistloci[loc].migeventbinnum;
    for(i0=0; i0< world->numpop2 ; i0++)
      {
        if(world->bayes->map[i0][1] == INVALID)
            continue;
        else
            i = world->bayes->map[i0][1];
        
        if(eventbinnum[i] > eventbinnum_allmax)
	  eventbinnum_allmax = eventbinnum[i];
      }
    //}
    doublevec2d(&eventbins_all,world->numpop2,eventbinnum_allmax);
    total = 0.0F;
    for (loc = 0; loc < loci; loc++)
      {
	eventbinnum = world->mighistloci[loc].migeventbinnum;
	eventbins = world->mighistloci[loc].migeventbins;
	for(i0=0; i0< numparams ; i0++)
	  {
	    if(world->bayes->map[i0][1] == INVALID)
	      continue;
	    else
	      i = world->bayes->map[i0][1];
	    
	    for(j=0 ; j < eventbinnum[i]; j++)
	      {
		if(eventbins[i][j][0] > 0.0)
		  {
		    eventbins_all[i][j] += (MYREAL) eventbins[i][j][0];
		    total += (MYREAL) eventbins[i][j][0];
		  }
	      }
	  }
      }
    bins = eventbinnum_allmax;

    set50 = (char *) mycalloc(bins+1, sizeof(char));
    set95 = (char *) mycalloc(bins+1, sizeof(char));
    memset (set50, '0', bins * sizeof(char));
    memset (set95, '0', bins * sizeof(char));
    
    binning = (float *) mycalloc (bins+1, sizeof (float));
    
    delta = world->mighistloci[0].eventbinsize;
    binning[0] = 0.5 * delta;
    for (i = 1; i < bins; i++)
      binning[i] = delta + binning[i - 1];
    
    for (i0 = (world->options->mighist_all ? 0 : world->numpop); i0 < numparams; i0++)
      {
	if(world->bayes->map[i0][1] == INVALID)
	  continue;
	else
	  i = world->bayes->map[i0][1];
	
	m2mm(i, world->numpop, &frompop, &topop);
	if(frompop==topop)
	  {
	    pdf_print_contents_at(lx-30,ly+125,"Freq. for ");
	    symbol_Theta(lx+12, ly+125,12,frompop+1);
	  }
	else
	  {
	    pdf_print_contents_at(lx-30,ly+125,"Freq. for ");
	    symbol_M(lx+12, ly+125, 12, frompop+1, topop+1, world->options->usem);
	  }
	pdf_print_contents_at(lx+90,ly-25, "Time [scaled by mutation rate / site / generation]");
	pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
	pdf_histogram(eventbins_all[i], set50, set95, bins, delta, 
		      0., -999, lx, ly, page_width - 55 - lx, 116,TRUE);
	page_height -= 160;
	if(i < (numparams-1))
	  {
	    pdf_page_advance_or_not(&page_height, 50);
	    ph = pdf_contents_get_height(canvas) - 55 - LINESTRETCH;
	    if(page_height >= ph)
	      page_height -= 150;
	    ly = page_height;
	  }
      }
    myfree(binning);
    myfree(set50);
    myfree(set95);
    free_doublevec2d(eventbins_all);
}


///
/// plot skyline histograms
void
pdf_skyline_histogram(long loci, long numparams,  world_fmt *world, boolean enlarged)
{
    float   page_width;
    float   lx;
    float   ly;
    float   ph;
    float   delta;    
    float   w;
    float * binning;
    float **confidence;
    float   c;
    float   lasttimebin;
    float   upperlimit;
    float maxtime;
    MYREAL *eventbin1max;
    //    float  confidencesum;
    long    clong;
    long    eventbinnum_allmax = 0L;
    long    loc;
    long    i;
    long    i0;
    long    j;
    long    bins;
    long    adjustbins;
    long    pop;
    long    frompop;
    long    topop;
    long  * eventbinnum;
    
    char  * set50;
    char  * set95;
    char    title[LINESIZE];
    
    
    MYREAL ** eventbins_all;
    MYREAL **std;
    
    tetra ** eventbins;
    
    
    if (loci > 1)
      {
        sprintf(title,"Summary of parameter values through %s over all loci",
                enlarged ? "RECENT time" : "time");
      }
    else
      {
        sprintf(title,"Parameter values through %s",enlarged ? "RECENT time" : "time");
      }
    // add a new page so that we can print at least four histograms
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    
    /* Start to print text. */ 
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    // first histogram position
    page_height -= 150;
    lx = 100;
    ly = page_height;
    
    //    for (loc = 0; loc < loci; loc++)
    // {
    loc = (loci > 1) ? loci : 0;
    eventbinnum = world->mighistloci[loc].eventbinnum;
    for(i0=0; i0< world->numpop2 ; i0++)
      {
        if(world->bayes->map[i0][1] == INVALID)
            continue;
        else
            i = world->bayes->map[i0][1];
        
        if(eventbinnum[i] > eventbinnum_allmax)
            eventbinnum_allmax = eventbinnum[i];
      }
    // }
    doublevec2d(&eventbins_all,world->numpop2, eventbinnum_allmax);
    doublevec2d(&std,world->numpop2, eventbinnum_allmax);
    floatvec2d(&confidence,world->numpop2, eventbinnum_allmax);
    eventbin1max = (MYREAL *) mycalloc(world->numpop2,sizeof(MYREAL));
    loc = (loci > 1) ? loci : 0;
    eventbinnum = world->mighistloci[loc].eventbinnum;
    eventbins = world->mighistloci[loc].eventbins;
    for(i0=0; i0< numparams ; i0++)
      {
        if(world->bayes->map[i0][1] == INVALID)
            continue;
        else
            i = world->bayes->map[i0][1];
        
        for(j=0 ; j < eventbinnum[i]; j++)
          {
            if(eventbins[i][j][1] > 0.0)
              {
                eventbins_all[i][j] = eventbins[i][j][0];
                std[i][j] = (MYREAL) eventbins[i][j][2] * 1.96; // standard error * 1.96
                //old e...[2] changed to stderr: std[i][j] = (MYREAL) eventbins[i][j][2]; // standard deviation
		//never used                std[i][j] /= 0.510204 * eventbins[i][j][4]; // standard error[=std/number in bin] * 1.96
                if(eventbins[i][j][1] > eventbin1max[i])
                    eventbin1max[i] = eventbins[i][j][1];
              }
          }
      }
    
    for(i0=0; i0< numparams ; i0++)
      {
        if(world->bayes->map[i0][1] == INVALID)
            continue;
        else
            i = world->bayes->map[i0][1];

        for(j=0 ; j < eventbinnum[i]; j++)
          {
	    //c = MAX(0.0,MIN(1.0,(float) (eventbins[i][j][1])));
	    clong = (long) eventbins[i][j][4];
	    if(clong < 100 )
	      c = 0.01;
	    else if(clong < 400)
	      c = 0.1;
	    else if (clong < 1600)
	      c = 0.5;
	    else if (clong < 3200)
	      c = 0.70;
	    else if (clong < 6400)
	      c = 0.8;
	    else if (clong < 12800)
	      c = 0.9;
	    else
	      c = 1.0;

            confidence[i][j] = 1.0 - c;
          }
      }
    
    bins = eventbinnum_allmax;
    
    set50 = (char *) mycalloc(bins+1, sizeof(char));
    set95 = (char *) mycalloc(bins+1, sizeof(char));
    memset (set50, '0', bins * sizeof(char));
    memset (set95, '0', bins * sizeof(char));
    
    binning = (float *) mycalloc (bins+1, sizeof (float));
    
    delta = world->mighistloci[0].eventbinsize;
    binning[0] = 0.5 * delta;
    for (i = 1; i < bins; i++)
        binning[i] = delta + binning[i - 1];
    
    for (i0 = (world->options->mighist_all ? 0 : world->numpop); i0 < numparams; i0++)
      {
        if(world->bayes->map[i0][1] == INVALID)
            continue;
        else
            i = world->bayes->map[i0][1];
        
        m2mm(i0, world->numpop, &frompop, &topop);
        if(frompop==topop)
          {
            symbol_Theta(lx-40, ly+125, 12, frompop+1);
            pdf_print_contents_at(lx+90,ly-25, "Time [scaled by mutation rate / generation (DNA: per site, other: per locus)]");
          }
        else
          {
            symbol_M(lx-40, ly+125, 12, frompop+1, topop+1, world->options->usem);
            pdf_print_contents_at(lx+90,ly-25, "Time [scaled by mutation rate / generation (DNA: per site, other: per locus)]");
          }
        pdf_contents_set_font_and_size(canvas, "Helvetica", 10);

	// focus on all recorded times or only on 0..mean of largest population
        if(!enlarged) 
          {
	    // all records
	    upperlimit = quantiler(eventbins_all[i], 1.00, 1.00, (long) bins/2,bins);
	    if(upperlimit < 2*world->bayes->histogram[loc].modes[i])
	      upperlimit =  2*world->bayes->histogram[loc].modes[i];
	    pdf_histogram_plus(eventbins_all[i], std[i], set50, set95, bins, delta, 
			       0., -9999, lx, ly, page_width - 55 - lx, 116, upperlimit, TRUE, world, confidence[i], topop);
          }
        else
          {
	    // up to mean of largest population
            loc = (loci > 1) ? loci : 0;
	    //            m2mm(i, world->numpop, &frompop, &topop);
	    maxtime = world->bayes->histogram[loc].means[0];
	    for(pop=1;pop<world->numpop;pop++)
	      {
		if(world->bayes->histogram[loc].means[pop]>maxtime)
		  maxtime=world->bayes->histogram[loc].means[pop];
	      }
	    //xcode j=eventbinnum[i]-1;
	    adjustbins = ceil(maxtime/delta);
	    // why could bin be smaller than adjustbins? 
            bins = bins < adjustbins ? bins : adjustbins;
	    if(bins==0)
	      continue;
	    lasttimebin = delta * bins;
	    // evaluates how much of the y-axis should be shown
	    upperlimit = quantiler(eventbins_all[i], 1.0, 0.95, (long) bins/2,bins);
            pdf_histogram_plus(eventbins_all[i], std[i], set50, set95, bins, delta, 
                               0., lasttimebin,
                               lx, ly, page_width - 55 - lx, 116, upperlimit ,TRUE, world, confidence[i], topop);
          }
        page_height -= 160;
        if(i < (numparams-1))
          {
            pdf_page_advance_or_not(&page_height, 50);
            ph = pdf_contents_get_height(canvas) - 55 - LINESTRETCH;
            if(page_height >= ph)
                page_height -= 150;
            ly = page_height;
          }
      }
    myfree(binning);
    myfree(set50);
    myfree(set95);
    free_doublevec2d(eventbins_all);
    free_doublevec2d(std);
    free_floatvec2d(confidence);
    myfree(eventbin1max);
}

void pdf_print_spectra(data_fmt *data, MYREAL ***freq, long * maxalleles)
{
    char **header;
    char ***elements;
    long mostalleles =0;
    long pop;
    long locus; 
    long a;
    long *maxallelepop;
    //  float position[10] = {55., 100., 150., 200., 250., 300., 350., 400., 450., 500.};
    float page_width;
    MYREAL allfreq;
    maxallelepop = (long *) mycalloc(data->numpop, sizeof(long));
    charvec2d(&header, data->numpop+2,LINESIZE);
    //print title
    page_width = pdf_contents_get_width(canvas);
    pdf_print_section_title(&page_width, &page_height, "Allele frequency spectra");
    // loop over loci
    pdf_advance(&page_height);
    // header
    sprintf(header[0],"Allele");
    for(pop=1; pop < data->numpop+1; pop++)
        sprintf(header[pop],"Pop%-3li",pop);
    sprintf(header[pop],"All");
    for(locus=0; locus < data->loci; locus++)
      {
        if(mostalleles < maxalleles[locus])
            mostalleles = maxalleles[locus];
      }
    elements = (char ***) mycalloc(mostalleles+1, sizeof(char**));
    for(a=0; a < mostalleles+1; a++)
        charvec2d(&elements[a],data->numpop+2, LINESIZE);
    
    
    for(locus=0; locus < data->loci; locus++)
      {
        memset(maxallelepop,'\0',sizeof(long)*data->numpop);
        pdf_printf(55.,page_height, 'L', "Locus %i",locus+1);
        pdf_advance(&page_height);
        for(a=0; a < maxalleles[locus]; a++)
          {
            sprintf(elements[a][0],"%6s ",data->allele[locus][a]);
            allfreq = 0.0;
            for(pop=0; pop < data->numpop; pop++)
              {
                maxallelepop[pop] += 1;
                sprintf(elements[a][pop+1],"%.3f",freq[pop][locus][a]);
                allfreq += freq[pop][locus][a];
              }
            sprintf(elements[a][pop+1],"%.3f",allfreq/data->numpop);
          }
        sprintf(elements[a][0],"Total");
        for(pop=0; pop < data->numpop; pop++)
          {
            sprintf(elements[a][pop+1],"%li",maxallelepop[pop]);
          }
        sprintf(elements[a][pop+1],"%li",maxalleles[locus]);
        pdf_table2(55., &page_height, data->numpop+2,maxalleles[locus], header, elements, NULL, 2, 10.0);
        pdf_advance(&page_height);
      }
    free_charvec2d(header);
    for(a=0; a < mostalleles+1; a++)
        free_charvec2d(elements[a]);
    myfree(elements);
    myfree(maxallelepop);
}

///
/// print average temperatures
void pdf_print_averageheat(world_fmt **universe, option_fmt *options)
{
  long a;
    char **header;
    char ***elements;
    //  float position[10] = {55., 100., 150., 200., 250., 300., 350., 400., 450., 500.};
    float page_width;
    charvec2d(&header, 2,LINESIZE);
    //print title
    page_width = pdf_contents_get_width(canvas);
    pdf_print_section_title(&page_width, &page_height, "Average temperatures during the run");
    pdf_advance(&page_height);
    // header
    sprintf(header[0],"Chain");
    sprintf(header[1],"Temperatures");
    elements = (char ***) mycalloc(options->heated_chains, sizeof(char**));
    for(a=0; a < options->heated_chains; a++)
        charvec2d(&elements[a],2, LINESIZE);
    for(a=0; a < options->heated_chains; a++)
      {
	sprintf(elements[a][0],"%5li ",a+1);
	sprintf(elements[a][1],"%10.5f ",universe[a]->averageheat);
      }
    pdf_table2(55., &page_height, 2, options->heated_chains, header, elements, NULL, 2, 10.0);
    pdf_advance(&page_height);
    pdf_advance(&page_height);
    pdf_printf_next(55., &page_height,"Adaptive heating often fails, if the average temperatures are very close together\n");
    pdf_printf_next(55., &page_height,"try to rerun using static heating! If you want to compare models using marginal\n");
    pdf_printf_next(55., &page_height,"likelihoods then you MUST use static heating\n");
    pdf_advance(&page_height);
    free_charvec2d(header);
    for(a=0; a < options->heated_chains; a++)
        free_charvec2d(elements[a]);
    myfree(elements);
}

///
/// table frequency of events through time for each locus and all loci
void pdf_print_eventtime_table(world_fmt *world)
{
    float   page_width;
    float   w;
    long    i;
    long    j;
    long    j0;
    long    frompop;
    long    topop;
    long  * eventbinnum;
    long locus;
    duo ** eventbins;
    float  age;
    char    title[LINESIZE];
    float left_margin = 55;
    
    sprintf(title,"Distribution of events trough time");
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    
    /* Start to print text. */ 
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    for(locus=0;locus < world->loci + (world->loci > 1 ? 1 : 0); locus++)
      {
	if(locus<world->loci)
	  {
	    if(world->data->skiploci[locus])
	      continue;
	  }
	eventbins = world->mighistloci[locus].migeventbins;
	eventbinnum = world->mighistloci[locus].migeventbinnum;
	// population sizes
	for(j0=0; j0 < world->numpop2; j0++)
	  {
	    if(world->bayes->map[j0][0] != INVALID)
	      {
		j = world->bayes->map[j0][1];
		
                pdf_advance(&page_height);
                if(locus == world->loci)
		  pdf_printf(left_margin, page_height, 'L' , "All loci");
                else
		  pdf_printf(left_margin, page_height, 'L' , "Locus %li:",locus);
                if(j0 < world->numpop)
		  symbol_Theta(left_margin+80, page_height, 12, j0+1); 
                else
                  {
                    m2mm (j0, world->numpop, &frompop, &topop);
                    symbol_M(left_margin+80, page_height, 12, frompop+1, topop+1, world->options->usem);
                  }	    
                pdf_advance(&page_height);
                pdf_printf(left_margin, page_height, 'L', "Time");
                pdf_printf(left_margin + 200, page_height, 'L', "Frequency of visit");
                pdf_printf(left_margin + 400, page_height, 'L', "MRCA frequency");
                pdf_advance_half(&page_height);
                pdf_draw_line(50, page_height, page_width-50, page_height);
                pdf_advance(&page_height);
                age = -world->mighistloci[locus].eventbinsize/2.;
                for(i = 0; i < eventbinnum[j]; i++)
                  {
                    age += world->mighistloci[locus].eventbinsize;
                    pdf_printf(left_margin, page_height, 'L', "%10.10f", age);
                    pdf_printf(left_margin+200, page_height, 'L', "%8.5f", (MYREAL) eventbins[j][i][0]); 
                    pdf_printf(left_margin+400, page_height, 'L', "%8.5f", (MYREAL) eventbins[j][i][1]); 
                    pdf_advance(&page_height);
                  }
                pdf_draw_line(50, page_height, page_width-50, page_height);
                pdf_advance(&page_height);
              }
          }
      }
}
///
/// average and median time for a events
void pdf_print_time_table(world_fmt *world, 
                          float *meantime, float *mediantime, float *stdtime, float *freq,
                          boolean mrca)
{
    float   page_width;
    float   w;
    //    long    j;
    //    long    j0;
    long    pop;
    long    lp;
    long    frompop;
    long    topop;
//    long  * eventbinnum;
    long    locus;
//    duo **  eventbins;
    char    title[LINESIZE];
    float   left_margin = 55;
    if(mrca)
        sprintf(title,"Time and probability of location of most recent common ancestor");
    else
        sprintf(title,"Summary statistics of events through time");
    pdf_new_page("");
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
    w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
    
    /* Start to print text. */ 
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);
    pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    page_height -= 126;
    pdf_draw_line(50, page_height, page_width-50, page_height);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    for(locus=0;locus < world->loci + (world->loci > 1 ? 1 : 0); locus++)
      {
	if(world->data->skiploci[locus])
	  continue;
      
      //xcode eventbins = world->mighistloci[locus].migeventbins;
      //xcode eventbinnum = world->mighistloci[locus].migeventbinnum;
        
        pdf_advance(&page_height);
        if(locus == world->loci)
            pdf_printf(left_margin, page_height, 'L' , "All loci");
        else
            pdf_printf(left_margin, page_height, 'L' , "Locus %li",locus+1);
        
        pdf_advance(&page_height);
        pdf_printf(left_margin, page_height, 'L', "Population");
        pdf_printf(left_margin + 100, page_height, 'L', "Time");
        pdf_printf(left_margin + 400, page_height, 'L', "Frequency");
        pdf_advance_half(&page_height);
        pdf_draw_line(left_margin + 100, page_height, left_margin+390, page_height);
        pdf_advance(&page_height);
        if(!mrca)
          {
            pdf_printf(left_margin, page_height, 'L', "From");
            pdf_printf(left_margin + 50, page_height, 'L', "To");
          }
        pdf_printf(left_margin + 100, page_height, 'L', "Average");
        pdf_printf(left_margin + 200, page_height, 'L', "Median");
        pdf_printf(left_margin + 300, page_height, 'L', "Std");
        pdf_advance_half(&page_height);
        pdf_draw_line(50, page_height, page_width-50, page_height);
        pdf_advance(&page_height);
        for (pop = (world->options->mighist_all ? 0 : world->numpop); 
             pop <  (mrca ? world->numpop : world->numpop2); pop++)
          {
	    if(world->bayes->map[pop][1] == INVALID)
	      continue;
            lp = locus * world->numpop2 + pop;
            m2mm(pop,world->numpop,&frompop,&topop);
            if(freq[lp]>0.0)
              {
                pdf_printf(left_margin,       page_height, 'L', "%li",frompop+1);
                if(!mrca)
                    pdf_printf(left_margin + 50,  page_height, 'L', "%li", topop+1);
                pdf_printf(left_margin + 100, page_height, 'L', "%f", meantime[lp]);
                pdf_printf(left_margin + 200, page_height, 'L', "%f", mediantime[lp]);
                pdf_printf(left_margin + 300, page_height, 'L', "%f", stdtime[lp]);
                pdf_printf(left_margin + 400, page_height, 'L', "%f", freq[lp]);
                pdf_advance(&page_height);
              }
          }
        pdf_advance(&page_height);
      }
}


///
/// plot event histograms
void
pdf_histogram_legend()
{
  float   page_width;
  float   lx;
  float   w;
  char    title[LINESIZE];

  // add a new page for the legend
  pdf_new_page("");
  sprintf(title,"%s","Legend for Skyline and Event plots");
  pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 16);
  w = (float) pdf_contents_get_text_width(canvas, title, NULL, NULL);
  
  /* Start to print text. */ 
  page_height = pdf_contents_get_height(canvas);
  page_width = pdf_contents_get_width(canvas);
  pdf_print_contents_at((page_width - w)/2, page_height - 100, title);
  pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
  page_height -= 126;
  pdf_draw_line(50, page_height, page_width-50, page_height);
  pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
  // first histogram position
  pdf_advance(&page_height);
  lx = 60;
  pdf_printf_next(lx-5, &page_height, "Skyline plots:");
  pdf_printf_next(lx, &page_height, "Skyline plots visualize the changes of population sizes and migration rates through time");
  pdf_printf_next(lx, &page_height, "(today is on the left side and time is measured into the past. The time scale is in units of");
  pdf_printf_next(lx, &page_height, "expected mutations per generation. To calculate the absolute time scale you must supply an");
  pdf_printf_next(lx, &page_height, "mutation rate per year and the duration of a  generation in years in the data option.");
  pdf_printf_next(lx, &page_height, "You can calculate the absolute time by multiplying the scale by generation time times ");
  pdf_printf_next(lx, &page_height, "mutation rate per year (per site for DNA; per locus for all other datatypes).");
  pdf_advance(&page_height);
  pdf_printf_next(lx, &page_height, "With estimated mutation rate only the combined rate modifier is plotted.");
  pdf_printf_next(lx, &page_height, "[this will change to  mutation rate plot].");
  pdf_advance(&page_height);
  pdf_printf_next(lx, &page_height, "The gray bars cover 1.96 * approximate standard error (std in file skylinefile/number");
  pdf_printf_next(lx, &page_height, "of observations in the bin) up and down from the expected value.");
  pdf_printf_next(lx, &page_height, "The bar with different shades of gray on top of each plot indicates the number of values that were used to");
  pdf_printf_next(lx, &page_height, "to calculate the expected value, white means there were very few and black means");
  pdf_printf_next(lx, &page_height, "that there were many thousands of samples per bin.");
  pdf_advance(&page_height);
  pdf_printf_next(lx, &page_height, "On some plots one can see red squares below the grayscale bar, these suggest that either the");
  pdf_printf_next(lx, &page_height, "upper quantile and/or the main value was higher than the visible part of the  axis.");
  pdf_advance(&page_height);
  pdf_advance(&page_height);

  pdf_printf_next(lx-5, &page_height, "Event histograms:");
  pdf_printf_next(lx, &page_height, "All accepted events (migration events, coalescent events) are recorded and their frequency");
  pdf_printf_next(lx, &page_height, "are shown as histograms over time with recent time on the left side. The frequency plots of");
  pdf_printf_next(lx, &page_height, "populations with constant size and constant immigration rates show histograms that are similar");
  pdf_printf_next(lx, &page_height, "to exponential distribution, if the populations come from a divergence model without migration");
  pdf_printf_next(lx, &page_height, "then the frequency of migration events can show a peak in the past.");
}

///
///
void findminmax(float *vals, const long n, float *min, float *max)
{
  long i;
#ifdef WINDOWS
  *min = 1000000;
  *max = -1000000;
#else
  *min = HUGE;
  *max = -HUGE;
#endif
  for(i=0; i < n; i++)
    {
      float v = (float) vals[i];
      if(*min > v)
	*min = v;
      if(*max < v)
	*max = v;
    }
  //  printf("min=%f, max=%f\n",*min,*max);
}

///
/// Line plot assumes that x and y are a series of coordinates, connection with a colored line with given thickness, dots are marked in color
void pdf_linedotplot(long n, float *x, float *y, float dotcolor[], float linecolor[], MYREAL linethickness, float * page_height, float width, float height, boolean has_dots, boolean has_line)
{
  float minx;
  float miny;
  float maxx;
  float maxy;
  //  float page_width = pdf_contents_get_width(canvas);
  float lx = 100;
  float ly = *page_height - height;
  long i;
  float red = linecolor[0];
  float green = linecolor[1];
  float blue = linecolor[2];
  float xrange;
  float yrange;
  *page_height = ly;
  pdf_advance(page_height);
  findminmax(x, n, &minx, &maxx);
  findminmax(y, n, &miny, &maxy);
  pdf_create_axes(VERTICAL, miny, maxy, 5, lx, ly, height);
  pdf_create_axes(HORIZONTAL, minx , maxx, 9, lx, ly, width);
  if(has_line)
    {
      pdf_contents_set_line_width(canvas, linethickness);
      pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(red, green, blue));
      xrange = maxx - minx;
      yrange = maxy - miny;
      width /= xrange;
      height/= yrange;
      for(i=0; i < n-1; i++)
	{
	  pdf_draw_line(lx+(x[i]-minx)*width, ly+(y[i]-miny)*height,lx+(x[i+1]-minx)*width,ly+(y[i+1]-miny)*height);
	}
      pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0.,0.,0.));//reset to black
    }
  if(has_dots)
    {
      for(i=0; i < n; i++)
	{
	  pdf_print_dot(lx+(x[i]-minx)*width, ly+(y[i]-miny)*height, 3, SQUARE, dotcolor);
	}
    } 
}


/*
==============================================================================
Likelihood ratio tests
==============================================================================
Over all loci
Legend for the LRT tables
-------------------------------------------------------------------------------
  Null-Hypothesis: your test model         | Log(likelihood) of test model
       =same=                                   | Log(likelihood) of full model
       full model (the model under which the    | Likelihood ratio test value
		   genealogies were sampled)                | Degrees of freedom of test
[Theta values are on the diagonal of the | Probability*
 Migration matrix, migration rates are    | Probability**
specified as M]                          | Akaike's Information Criterion***
                                         | Number of parameters used
-------------------------------------------------------------------------------
  *) Probability under the assumption that parameters have range -Inf to Inf
 **) Probability under the assumption that parameters have range 0 to Inf
***) AIC: the smaller the value the better the model
          [the full model has AIC=-2.396909, num(param)=4]
*/
#define OFFSET 300
#define OFFSET2 280
#define ADDED 70
void pdf_print_lrt_header(MYREAL aicfull, long aicfullparamnum,  float *page_width,  float *page_height)
{
  float lx=55;
  float tempheight ;
  pdf_print_section_title(page_width, page_height, "Approximate Likelihood Ratio Tests");
  pdf_advance(page_height);
  lx = 60;
  pdf_printf_next(lx, page_height, "Legend for the likelihood ratio tables");
  pdf_draw_line(lx, *page_height, *page_width-lx, *page_height);
  tempheight = *page_height;
  pdf_advance(page_height);
  pdf_printf(lx, *page_height, 'L', "Null-Hypothesis: your test model");
  pdf_printf_next(lx+OFFSET, page_height, "Log(likelihood) of test model");
  //
  pdf_printf(lx, *page_height, 'L', "is equal to");
  pdf_printf_next(lx+OFFSET, page_height, "Log(likelihood) of full model");
  pdf_printf(lx, *page_height, 'L', "full model (the model under which the");
  pdf_printf_next(lx+OFFSET, page_height, "Likelihood ratio test value");
  pdf_printf(lx, *page_height, 'L', "genealogies were sampled)");
  pdf_printf_next(lx+OFFSET, page_height, "Degrees of freedom of test");
  pdf_printf_next(lx+OFFSET, page_height, "[Theta values are on the diagonal of the");
  pdf_printf_next(lx+OFFSET, page_height, "Probability*");
  pdf_printf(lx, *page_height, 'L',"Migration matrix, migration rates are");
  pdf_printf_next(lx+OFFSET, page_height, "Probability**");
  pdf_printf(lx, *page_height, 'L', "specified as M]");
  pdf_printf_next(lx+OFFSET, page_height, "Akaike's Information Criterion***");
  pdf_printf_next(lx+OFFSET, page_height, "Number of parameters used");
  pdf_draw_line(lx, *page_height, *page_width-lx, *page_height);
  pdf_draw_line(OFFSET2+lx, tempheight, OFFSET2+lx, *page_height);
  pdf_advance(page_height);
  pdf_printf_next(lx, page_height, "  *) Probability under the assumption that parameters have range -Inf to Inf");
  pdf_printf_next(lx, page_height, " **) Probability under the assumption that parameters have range 0 to Inf");
  pdf_printf_next(lx, page_height, "***) AIC: the smaller the value the better the model");
  pdf_printf_next(lx, page_height, "     [the full model has AIC=%f, num(param)=%li]",aicfull, aicfullparamnum);
}

/*
-------------------------------------------------------------------------------
H0: 0.0085 0.0000 262.32 0.0262                    | LnL(test) = 4.213529
 =  0.0085 12.065 262.32 0.0262                    | LnL(full) = 5.198454
[ *, 0, *, *,]                                     | LRT       = 1.969850
                                                   | df        = 1
                                                   | Prob      = 0.160464
                                                   | Probc     = 0.047160
                                                   | AIC       = -2.427058
                                                   | num(param)= 3
-------------------------------------------------------------------------------
*/
//remember: param0 is the parameterset to test and
// NOT the parameterset from migrate.
#define BOXSIZE 50 /*defines the print width of the left box containing the H0*/
#define BOXSIZE2 40 /*defines the print width of the left box with the legend*/
void pdf_print_LRT_box(world_fmt* world,char **box,long size,
		       MYREAL like0, MYREAL like1,MYREAL lrt, MYREAL chiprob, MYREAL chiprob2, MYREAL aic, long df, long aicparamnum, float *page_height, float *page_width)
{
  float lx=60;
  float tempheight;
  pdf_page_advance_or_not(page_height, 200.);
  tempheight = *page_height;
  //printf("%f\n",tempheight);
  pdf_draw_line(lx, *page_height, *page_width-lx, *page_height);
  pdf_advance(page_height);
  pdf_printf(lx, *page_height, 'L', "%-*.*s",BOXSIZE,BOXSIZE, box[0]);  
  pdf_printf(lx+OFFSET, *page_height, 'L', "LnL(test)");
  pdf_printf_next(lx+OFFSET+ADDED, page_height, "= %f",like0);
  pdf_printf(lx, *page_height, 'L', "%-*.*s",BOXSIZE,BOXSIZE, box[1]);  
  pdf_printf(lx+OFFSET, *page_height, 'L', "LnL(full)");
  pdf_printf_next(lx+OFFSET+ADDED, page_height, "= %f",like1);
  pdf_printf(lx, *page_height, 'L', "%-*.*s",BOXSIZE,BOXSIZE, box[2]);  
  pdf_printf(lx+OFFSET, *page_height, 'L', "LRT");
  pdf_printf_next(lx+OFFSET+ADDED, page_height, "= %f", lrt);
  pdf_printf(lx, *page_height, 'L', "%-*.*s",BOXSIZE,BOXSIZE, box[3]);  
  pdf_printf(lx+OFFSET, *page_height, 'L', "df");
  pdf_printf_next(lx+OFFSET+ADDED, page_height, "= %li",(long) df);
  pdf_printf(lx, *page_height, 'L', "%-*.*s",BOXSIZE,BOXSIZE, box[4]);  
  pdf_printf(lx+OFFSET, *page_height, 'L', "Prob");
  pdf_printf_next(lx+OFFSET+ADDED, page_height, "= %f",chiprob);
  pdf_printf(lx, *page_height, 'L', "%-*.*s",BOXSIZE,BOXSIZE, box[5]);  
  pdf_printf(lx+OFFSET, *page_height, 'L', "Probc");
  pdf_printf_next(lx+OFFSET+ADDED, page_height, "= %f",chiprob2);
  pdf_printf(lx, *page_height, 'L', "%-*.*s",BOXSIZE,BOXSIZE, box[6]);  
  pdf_printf(lx+OFFSET, *page_height, 'L', "AIC");
  pdf_printf_next(lx+OFFSET+ADDED, page_height, "= %f",aic);
  pdf_printf(lx, *page_height, 'L', "%-*.*s",BOXSIZE,BOXSIZE, box[7]);  
  pdf_printf(lx+OFFSET, *page_height, 'L', "num(param)");
  pdf_printf_next(lx+OFFSET+ADDED, page_height, "= %li",aicparamnum);
  //pdf_printf(lx+OFFSET, *page_height, 'L', "num(param)= %li",aicparamnum);
  pdf_draw_line(lx, *page_height, *page_width-lx, *page_height);
  pdf_draw_line(OFFSET2+lx, tempheight, OFFSET2+lx, *page_height);
  //pdf_advance(page_height);
}
#endif /*end of PRETTY*/













