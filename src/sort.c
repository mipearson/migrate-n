/* ----------------------------------------------------- */
/* sort.c                       */
/* comparison routeins for various qsorts and bsearch's */
/* ----------------------------------------------------- */

/*                              */
/* P. Beerli                        */
/*
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
$Id: sort.c 1722 2010-10-09 15:27:25Z beerli $ */
/* ----------------------------------------------------- */
/*! \file sort.c */


/* include files */
#include "migration.h"
#include "sort.h"
/* private functions */

/* public functions */
int
charcmp (const void *v1, const void *v2)
{
    if (*(char *) v1 < *(char *) v2)
    {
        return -1;
    }
    else
    {
        if (*(char *) v1 > *(char *) v2)
        {
            return 1;
        }
        else
            return 0;
    }
}

int
stringcmp (const void *v1, const void *v2)
{
    if (strcmp ((char *) v1, (char *) v2) < 0)
    {
        return -1;
    }
    else
    {
        if (strcmp ((char *) v1, (char *) v2) > 1)
            return 1;
        else
            return 0;
    }
}

int
numcmp (const void *v1, const void *v2)
{
    if (*(MYREAL *) v1 < *(MYREAL *) v2)
    {
        return -1;
    }
    else
    {
        if (*(MYREAL *) v1 > *(MYREAL *) v2)
        {
            return 1;
        }
        else
            return 0;
    }
}

int
paircmp (const void *v1, const void *v2)
{
    if (((MYREAL *) v1)[1] < ((MYREAL *) v2)[1])
    {
        return -1;
    }
    else
    {
        if (((MYREAL *) v1)[1] > ((MYREAL *) v2)[1])
        {
            return 1;
        }
        else
            return 0;
    }
}

int
paircmp_first (const void *v1, const void *v2)
{
    if (((MYREAL *) v1)[0] < ((MYREAL *) v2)[0])
    {
        return -1;
    }
    else
    {
        if (((MYREAL *) v1)[0] > ((MYREAL *) v2)[0])
        {
            return 1;
        }
        else
            return 0;
    }
}



void paired_qsort2(pair *x, long xelem)
{
    long oldli;
    long li;
    long elements;
    
    // sort using the second element
    qsort(x,xelem,sizeof(pair),paircmp);
    
    // with elements of the same histogram size sort by their first element
    oldli = 0;
    elements = 1;
    for(li=1; li < xelem; li++)
    {
        if(x[li][1] == x[li-1][1])
        {
            elements++;
            continue;
        }
        else
        {
            if(elements > 1)
                qsort(x+oldli,elements,sizeof(pair),paircmp_first);
            elements=1;
            oldli = li;
        }
    }
}


int
floatcmp (const void *v1, const void *v2)
{
    if (*(float *) v1 < *(float *) v2)
    {
        return -1;
    }
    else
    {
        if (*(float *) v1 > *(float *) v2)
        {
            return 1;
        }
        else
            return 0;
    }
}


int
longcmp (const void *v1, const void *v2)
{
    if (*(long *) v1 < *(long *) v2)
    {
        return -1;
    }
    else
    {
        if (*(long *) v1 > *(long *) v2)
        {
            return 1;
        }
        else
            return 0;
    }
}

int
intcmp (const void *v1, const void *v2)
{
    if (*(int *) v1 < *(int *) v2)
    {
        return -1;
    }
    else
    {
        if (*(int *) v1 > *(int *) v2)
        {
            return 1;
        }
        else
            return 0;
    }
}
 
int
agecmp (const void *x, const void *y)
{
  vtlist * xx = (vtlist *) x;
  vtlist * yy = (vtlist *) y;
  MYREAL xa = xx->age;
  MYREAL ya = yy->age;

  if (xa < ya)
    {
        return -1;
    }
    else
    {
        if (xa != ya)
        {
            return 1;
        }
        else
            return 0;
    }
}

int
delcmp (const void *x, const void *y)
{
    if ((*((node **) x))->id < (*((node **) y))->id)
    {
        return -1;
    }
    else
    {
        if ((*((node **) x))->id == (*((node **) y))->id)
        {
            return 0;
        }
        else
            return 1;
    }
}

int
migr_time_cmp (const void *x, const void *y)
{
    if (((migr_table_fmt *) x)->time < ((migr_table_fmt *) y)->time)
    {
        return -1;
    }
    else
    {
        if (((migr_table_fmt *) x)->time == ((migr_table_fmt *) y)->time)
        {
            return 0;
        }
        else
            return 1;
    }
}


int
searchagecmp (const void *x, const void *y)
{
    MYREAL xx = (MYREAL) *((MYREAL *) x);
    MYREAL age = (MYREAL) ((vtlist *) y)->age;

    if (xx < age)
    {
        return -1;
    }
    else
    {
        if (xx == age)
        {
            return 0;
        }
        else
            return 1;
    }

}


int
aiccmp (const void *x, const void *y)
{
    if (((aic_fmt *) x)->aic < ((aic_fmt *) y)->aic)
    {
        return -1;
    }
    else
    {
        if (((aic_fmt *) x)->aic == ((aic_fmt *) y)->aic)
        {
            return 0;
        }
        else
            return 1;
    }
}
