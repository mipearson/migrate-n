/* spline.c -- translated by f2c (version 19950110)
   from src.f from BVSPIS package (algorithm 7.. from
   the AMS library
 
   the original fortran source is described in:
 
   translated and adapted to interface with the
   profile likelihood code in migrate-n by
   Peter Beerli 1999
 
   how to call these functions: see an example in profile.c
   
   
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
$Id: spline.c 1807 2011-03-18 20:16:19Z beerli $
*/

/* \file spline.c

not used yet

*/

#include <math.h>
#include "migration.h"
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#include <string.h>


#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif


#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a))
#endif
#ifndef MIN
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#endif

#define pow_di(a,b) (pow((MYREAL) (*(a)),(MYREAL) (*(b))))
#define pow_ii(a,b) ((long) pow((MYREAL) (*(a)),(MYREAL) (*(b))))

int
dintrs_ (
	 MYREAL *a, MYREAL *b, MYREAL *c, MYREAL *d, MYREAL *p1, MYREAL *p2);
boolean
dtst_ (MYREAL *a, MYREAL *b, MYREAL *c, MYREAL *d, MYREAL *eps);

 int
dal2_ (
       MYREAL *a2,
       long *np,
       MYREAL *info,
       long *comm, long *part,
       MYREAL *d);

int
dal2dp_ (
	 MYREAL *a2,
	 long *np,
	 MYREAL *info,
	 long *comm, long *part,
	 MYREAL *d);

int
dprj2_ (
	MYREAL *a, MYREAL *b, MYREAL *c, MYREAL *d,
	long *i,
	MYREAL *info,
	long *comm, long *part, long *np,
	MYREAL *p1, MYREAL *p2);

MYREAL
dsl_ (MYREAL *a, MYREAL *b, MYREAL *c);

long
dmnind_ (
	 MYREAL *d, MYREAL *part);

int
dbsear_ (
	 MYREAL *x,
	 long *np,
	 MYREAL *xtab,
	 long *ind);

int
dlspis_ (
	 MYREAL *x, MYREAL *y, MYREAL *d, MYREAL *d2,
	 long *np, long *n, long *k, long *ind,
	 MYREAL *l);


int
dfpsvf_ (
	 MYREAL *a1, MYREAL *a2,
	 long *np,
	 MYREAL *info,
	 long *comm, long *part,
	 MYREAL *eps,
	 long *maxstp,
	 MYREAL (*beta) (MYREAL *),
	 long *errc,
	 MYREAL *dstar);

int
dprj0_ (
	long *i,
	MYREAL *info,
	long *comm, long *part, long *np,
	MYREAL *p1, MYREAL *p2);


MYREAL
d_sign (MYREAL *a, MYREAL *b)
{
    MYREAL x;
    x = (*a >= 0 ? *a : -*a);
    return (*b >= 0 ? x : -x);
}

int dstinf_ (
	 long *opt,
	 MYREAL *d0, MYREAL *dnp,
	 long *constr, long *n, long *k,
	 MYREAL *x, MYREAL *y, MYREAL *d,
	 long *np, long *comm, long *part,
	 MYREAL *eps,
	 MYREAL (*beta) (MYREAL *), MYREAL (*betai) (MYREAL *),
	 MYREAL *daux2, MYREAL *daux3, MYREAL *info,
	 long *errc, long *diagn);

int
dscdrc_ (
	 long *n,
	 MYREAL *x, MYREAL *y, MYREAL *d,
	 long *opt, long *np,
	 MYREAL *eps, MYREAL *d20, MYREAL *d2np,
	 MYREAL (*rho) (MYREAL *), MYREAL (*rhoi) (MYREAL *),
	 MYREAL *a1, MYREAL *a2, MYREAL *h, MYREAL *d2,
	 long *errc);

int
dbve_ (
       MYREAL *x, MYREAL *y,
       long *np, long *n, long *k,
       MYREAL *xtab,
       long *ntab, long *sbopt, long *y0opt, long *y1opt, long *y2opt,
       MYREAL *d, MYREAL *d2,
       long *errc,
       MYREAL *tb, MYREAL *l, MYREAL *laux0, MYREAL *laux1,MYREAL *laux2,MYREAL *y0tab,MYREAL *y1tab,MYREAL *y2tab,
       long *erre);

int
dtrmb_ (
	long *n,
	MYREAL *tb);

int
dsqtab_ (
	 MYREAL *x, MYREAL *y,
	 long *np,
	 MYREAL *xtab,
	 long *ntab, long *y0opt, long *y1opt, long *y2opt, long *n, long *k,
	 MYREAL *d, MYREAL *d2, MYREAL *tb, MYREAL *l, MYREAL *laux0, MYREAL *laux1, MYREAL *laux2, MYREAL *y0tab, MYREAL *y1tab, MYREAL *y2tab);

int
dtdc_ (
       long *np,
       MYREAL *x,
       long *comm, long *part,
       MYREAL *info, MYREAL *dd2, MYREAL *dd3);

int
dprj1_ (
	MYREAL *a, MYREAL *b, MYREAL *c, MYREAL *d,
	long *i,
	MYREAL *info,
	long *comm, long *part, long *np,
	MYREAL *p1, MYREAL *p2);


/* Table of constant values */

/*static*/
long c__5 = 5;
/*static*/
long c__2 = 2;
/*static*/
long c__0 = 0;
/*static*/
long c__1 = 1;

/* Subroutine */ int
dalg1_ (
	MYREAL *a1,
	long *np,
	MYREAL *info,
	long *comm, long *part,
	MYREAL *eps, MYREAL *a2,
	long *errc, long *diagn)
{
    /* Initialized data */

    /*static */ MYREAL fl0 = 0.;

    /* System generated locals */
    long i__1;

    /* Local variables */
    /*static */
    long errc1;
    //    extern /* Subroutine */ int dprj1_ (MYREAL *);
    /*static */
    long i;

    /*  DALG1 implements the algorithm A1[B(0)] described in subr. DBVSSC. */

    /*  The input parameters NP,COMM,PART,EPS and the output parameters */
    /*  ERRC, DIAGN are described in DBVSSC. The input parameter INFO is */
    /*  described in DSTINF. */

    /*  Items of possible interest are: */

    /*  A1: floating array, of bounds 1:2, 0:NP, containing the sequence of */
    /*      the sets  B(i), i=0,1,...,np (see the comments in DBVSSC). */
    /*      More precisely,  B(i) = [a1(1,i),a1(2,i)] . */

    /*  A2: floating array, of bounds 1:2, 0:NP, containing the sequence of */
    /*      the sets  A(i), i=0,1,...,np (see the comments in DBVSSC). */
    /*      More precisely, A(i) = [a2(1,i),a2(2,i)] . */
    /* Parameter adjustments */
    --a2;
    --a1;
    --info;

    /* Function Body */
    errc1 = 0;
    /*  Step 1. */
    a2[1] = a1[1];
    a2[2] = a1[2];
    /*  Step 2. */
    i__1 = *np;
    for (i = 1; i <= i__1; ++i)
    {
L10:
        /*  Step 2.1. */
        dprj1_ (&a2[((i - 1) << 1) + 1], &a2[((i - 1) << 1) + 2],
                &a1[(i << 1) + 1], &a1[(i << 1) + 2], &i, &info[1], comm, part,
                np, &a2[(i << 1) + 1], &a2[(i << 1) + 2]);
        /*  Ignore the constraints in  [x(i),x(i+1)]  when A(i) is empty. */
        if (a2[(i << 1) + 1] > a2[(i << 1) + 2] + *eps)
        {
            info[*comm + 1 + (i - 1)] = fl0;
            errc1 = 1;
            diagn[i - 1] = 1;
            goto L10;
        }
        /* L20: */
    }
    if (errc1 == 1 && *errc == 9)
    {
        *errc = 10;
    }
    else if (errc1 == 1)
    {
        *errc = 6;
    }
    return 0;
}    /* dalg1_ */

/*  --------------------------------------------------------------------- */
/* Subroutine */ int
dalg3_ (
	MYREAL *info,
	long *np, long *comm, long *part, long *opt,
	MYREAL *d0, MYREAL *dnp, MYREAL *eps,
	long *kmax, long *maxstp,
	MYREAL (*beta) (MYREAL *), MYREAL (*betai) (MYREAL *),
	MYREAL *a1, MYREAL *a2, MYREAL *d,
	long *errc, long *diagn)
{
    /* Initialized data */

    /*static */ MYREAL fl2 = 2.;

    /* System generated locals */
    long i__1;
    MYREAL d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    //    extern boolean dtst_ ();

    long i, k, p, q;
    MYREAL dstar, p1, p2;
    //extern /* Subroutine */ int dal2dp_ (), dfpsvf_ (), dintrs_ (), dal2_ ();

    /*  DALG3 computes a sequence of slopes ( d(0), d(1), ..., d(np) ) which
    */
    /*  can be used to compute a shape-preserving interpolating spline with */
    /*  or without boundary conditions, as requested by the user. It is an */
    /*  implementation of the algorithm A3 described in subroutine DBVSSC. */

    /*  The input parameters NP,COMM,PART,OPT,EPS,KMAX,MAXSTP,BETA,BETAI,D */
    /*  and the output parameter ERRC are described in subprogram DBVSSC. */
    /*  The input parameter INFO is described in subprogram DSTINF. */
    /* Parameter adjustments */
    --a2;
    --a1;
    --info;

    /* Function Body */
    p = *opt / 100;
    q = *opt % 100 / 10;
    /*  If kmax.LE.0 it is reset to ink. */
    if (*kmax <= 0)
    {
        *kmax = 10;
    }
    /*  If maxstp.LE.0 it is reset to instp. */
    if (*maxstp <= 0)
    {
        *maxstp = 10;
    }
    /*  Start step 1: store the sets  B(i), i=0,1,...,np , into the array A1.
    */
    i__1 = *np;
    for (i = 1; i <= i__1; ++i)
    {
        dprj0_ (&i, &info[1], comm, part, np, &a1[((i - 1) << 1) + 1],
                &a1[((i - 1) << 1) + 2]);
        /* L10: */
    }
    a1[(*np << 1) + 1] = info[4];
    a1[(*np << 1) + 2] = info[5];
    /*  Reset the first and the last interval if separable boundary condtions
    */
    /*  are required */
    if (q == 3)
    {
        a1[1] = *d0;
        a1[2] = *d0;
        a1[(*np << 1) + 1] = *dnp;
        a1[(*np << 1) + 2] = *dnp;
    }
    /*  Start step 2. Call DALG1 to compute the array A2 containing the */
    /*  sets A(i) , i=0,1,...,np. */
    dalg1_ (&a1[1], np, &info[1], comm, part, eps, &a2[1], errc, diagn);
    /*  Start step 3 (steps 3-7 are activated only if boundary conditions are
    */
    /*  required). */
    if (q == 2)
    {
        /*  Compute  betai(phi(A(0).INT.beta(A(0))) . */
        /* Computing MIN */
        d__2 = (*beta) (&a2[1]), d__3 = (*beta) (&a2[2]);
        d__1 = MIN (d__2, d__3);
        /* Computing MAX */
        d__5 = (*beta) (&a2[1]), d__6 = (*beta) (&a2[2]);
        d__4 = MAX (d__5, d__6);
        dintrs_ (&d__1, &d__4, &a2[(*np << 1) + 1], &a2[(*np << 1) + 2], &p1,
                 &p2);
        /* Computing MIN */
        d__1 = (*betai) (&p1), d__2 = (*betai) (&p2);
        a1[1] = MIN (d__1, d__2);
        /* Computing MAX */
        d__1 = (*betai) (&p1), d__2 = (*betai) (&p2);
        a1[2] = MAX (d__1, d__2);
        /*  Start step 4. */
        if (p1 > p2 + *eps)
        {
            *errc = 5;
            return 0;
        }
        /*  Start step 5 : initialization */
        k = 1;
        i__1 = *np;
        for (i = 1; i <= i__1; ++i)
        {
            a1[(i << 1) + 1] = a2[(i << 1) + 1];
            a1[(i << 1) + 2] = a2[(i << 1) + 2];
            /* L20: */
        }
        /*  Iteration. The loop is stopped if a convergence test is satisfied
        */
        /*  or kmax iterations have already been done. */
L30:
        if (dtst_ (&a1[1], &a1[2], &a2[1], &a2[2], eps) || k == *kmax)
        {
            goto L50;
        }
        /*  Step 5.1 . */
        dalg1_ (&a1[1], np, &info[1], comm, part, eps, &a2[1], errc, diagn);
        /*  Step 5.2 . */
        /* Computing MIN */
        d__2 = (*beta) (&a2[1]), d__3 = (*beta) (&a2[2]);
        d__1 = MIN (d__2, d__3);
        /* Computing MAX */
        d__5 = (*beta) (&a2[1]), d__6 = (*beta) (&a2[2]);
        d__4 = MAX (d__5, d__6);
        dintrs_ (&d__1, &d__4, &a2[(*np << 1) + 1], &a2[(*np << 1) + 2], &p1,
                 &p2);
        /* Computing MIN */
        d__1 = (*betai) (&p1), d__2 = (*betai) (&p2);
        a1[1] = MIN (d__1, d__2);
        /* Computing MAX */
        d__1 = (*betai) (&p1), d__2 = (*betai) (&p2);
        a1[2] = MAX (d__1, d__2);
        /*  If  gamma(A(0))  is empty for some k the problem does not have any
         */
        /*  solution. */
        if (p1 > p2 + *eps)
        {
            *errc = 5;
            return 0;
        }
        i__1 = *np;
        for (i = 1; i <= i__1; ++i)
        {
            a1[(i << 1) + 1] = a2[(i << 1) + 1];
            a1[(i << 1) + 2] = a2[(i << 1) + 2];
            /* L40: */
        }
        ++k;
        goto L30;
L50:
        /*  Start step 7. */
        /*  Assign to dstar a suitable value */
        if (p == 1)
        {
            dstar = (a1[1] + a1[2]) / fl2;
        }
        else
        {
            dstar = info[*comm + *part * *np + 1];
        }
        /*  Check if dstar solves the problem, that is,  beta(dstar)  belongs
        to */
        /*  phi(dstar); if it is not the case, another value for dstar */
        /*  is looked for. */
        dfpsvf_ (&a1[1], &a2[1], np, &info[1], comm, part, eps, maxstp, beta,
                 errc, &dstar);
        if (*errc != 0 && *errc != 6 && *errc != 9 && *errc != 10)
        {
            return 0;
        }
        info[*comm + *part * *np + *np + 1] = (*beta) (&dstar);
        a2[(*np << 1) + 1] = (*beta) (&dstar);
        a2[(*np << 1) + 2] = (*beta) (&dstar);
    }
    /*  Start step 8. */
    if (p == 1)
    {
        dal2_ (&a2[1], np, &info[1], comm, part, d);
    }
    else
    {
        dal2dp_ (&a2[1], np, &info[1], comm, part, d);
    }
    return 0;
}    /* dalg3_ */

/* Subroutine */ int
dbvc_ (
       MYREAL *x, MYREAL *y,
       long *np, long *n, long *k, long *opt,
       MYREAL *d0, MYREAL *dnp, MYREAL *d20, MYREAL *d2np,
       long *constr,
       MYREAL *info,
       long *comm, long *part,
       MYREAL *eps,
       long *kmax, long *maxstp,
       MYREAL *a1, MYREAL *a2, MYREAL *daux2, MYREAL *daux3,
       MYREAL (*beta) (MYREAL *), MYREAL (*betai) (MYREAL *), MYREAL (*rho) (MYREAL *), MYREAL (*rhoi) (MYREAL *),
       MYREAL *d, MYREAL *d2,
       long *errc, long *diagn)
{
    /* System generated locals */
    long i__1;

    /* Local variables */
    //  extern  int dalg3_ ();
    long i, p, q, r;
    //  extern  int dscdrc_ (), dstinf_ ();

    /*  DBVC checks input parameters and computes the required spline. */

    /*  The input parameters X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,COMM,PART,
    */
    /*  EPS,KMAX,MAXSTP,BETA,BETAI,RHO,RHOI and the output parameters */
    /*  D,D2,ERRC,DIAGN are described in subroutine DBVSSC. */
    /*  The other parameters are described in the called subprograms. */
    /*  Check the input parameters NP and OPT. */
    /* Parameter adjustments */
    --daux2;
    --a2;
    --a1;
    --info;

    /* Function Body */
    p = *opt / 100;
    q = *opt % 100 / 10;
    r = *opt % 10;
    if (*np < 2 || (p < 1 || p > 4) || (q < 1 || q > 3) || (r < 0 || r > 4)
            || (p == 1 && r == 0))
    {
        *errc = 1;
        return 0;
    }
    /*  Check the array CONSTR. */
    if (*opt % 10 == 4)
    {
        i__1 = *np - 1;
        for (i = 0; i <= i__1; ++i)
        {
            if (constr[i] < 0 || constr[i] > 3)
            {
                *errc = 2;
                return 0;
            }
            /* L10: */
        }
    }
    /*  Check the input parameters N and K. */
    if (*k < 1 || *k > 2 || *n < *k * 3)
    {
        *errc = 3;
        return 0;
    }
    /*  Check the abscissas of the interpolation points. */
    i__1 = *np - 1;
    for (i = 0; i <= i__1; ++i)
    {
        if (x[i + 1] <= x[i])
        {
            *errc = 4;
            return 0;
        }
        /* L20: */
    }
    /*  Call subprogram DSTINF to set the information array INFO. */
    /*  Initialize the array DIAGN. */
    i__1 = *np - 1;
    for (i = 0; i <= i__1; ++i)
    {
        diagn[i] = 0;
        /* L30: */
    }
    dstinf_ (opt, d0, dnp, constr, n, k, x, y, d, np, comm, part, eps, beta,
             betai, &daux2[1], daux3, &info[1], errc, diagn);
    /*  Call subprogram DALG3 to compute the array D containing the first */
    /*  derivative at initial points. */
    dalg3_ (&info[1], np, comm, part, opt, d0, dnp, eps, kmax, maxstp, beta,
            betai, &a1[1], &a2[1], d, errc, diagn);
    if (*k == 2)
    {
        /*  A  C(2) spline is required. Compute the sequence of second derivat
        i- */
        /*  ves d2(i), i=0,...,np , according to the shape constraints and, if
         */
        /*  possible, to boundary conditions. */
        dscdrc_ (n, x, y, d, opt, np, eps, d20, d2np, rho, rhoi, &a1[1], &a2[1],
                 daux3, d2, errc);
    }
    return 0;
}    /* dbvc_ */


/*int dbvssc_ (x, y, np, n, k, opt, d0, dnp, d20, d2np, constr,
  eps, beta, betai, rho, rhoi, kmax, maxstp, errc, d, d2, diagn, work,
  nwork)
     MYREAL *x, *y;
     long *np, *n, *k, *opt;
     MYREAL *d0, *dnp, *d20, *d2np;
     long *constr;
     MYREAL *eps;
     MYREAL (*beta) (MYREAL *), (*betai) (MYREAL *), (*rho) (MYREAL *), (*rhoi) (MYREAL *);
     long *kmax, *maxstp, *errc;
     MYREAL *d, *d2;
     long *diagn;
     MYREAL *work;
     long *nwork;
*/
int
dbvssc_ (MYREAL *x, MYREAL *y, long *np, long *n, long *k, long *opt,
         MYREAL *d0, MYREAL *dnp, MYREAL *d20, MYREAL *d2np, long *constr,
         MYREAL *eps, MYREAL (*beta) (MYREAL *), MYREAL (*betai) (MYREAL *), MYREAL (*rho) (MYREAL *),
         MYREAL (*rhoi) (MYREAL *), long *kmax, long *maxstp, long *errc, MYREAL *d,
         MYREAL *d2, long *diagn, MYREAL *work, long *nwork)
{
  //    extern /* Subroutine */ int dbvc_ ();
    /*static */
    long i1, i2, i3, i4, i5, i6, i7, i8, i9, i10;

    /*  ------------------------------------------------- */
    /*            Lines 49-549 are comment lines. */
    /*            The actual code begins at line 555. */
    /*  ------------------------------------------------- */
    /*  ABSTRACT: */

    /*  DBVSSC is designed to compute the coefficients (first and, if */
    /*  appropriate, second derivatives) of a shape-preserving spline, of */
    /*  continuity class C(k), k=1,2 , which interpolates a set of data */
    /*  points and, if required, satisfies additional boundary conditions. */
    /*  DBVSSC furnishes the input parameters for DBVSSE, which, in turn, */
    /*  provides to evaluate the spline and its derivatives at a set of */
    /*  tabulation points. */

    /*  The user is allowed to use the following options: */

    /*  - to compute a spline subject to: */
    /*        - no constraint, */
    /*        - monotonicity constraints, */
    /*        - convexity constraints, */
    /*        - monotonicity and convexity constraints, */
    /*        - one of the above constraints in each subinterval; */

    /*  - to impose separable or non-separable boundary conditions on the */
    /*    spline, */

    /*  - to assign the first derivatives d(i), i=0,1,...,np , in input or to
    */
    /*    compute them from the constraints only or as the best approximation
    */
    /*    to a set of optimal values. Although the final sequence of */
    /*    derivatives does in any case satisfy the imposed restrictions on */
    /*    the shape, the resulting graphs may exhibit different behaviors. */


    /*  REMARK: */

    /*  In these comments variable and array names will be denoted with */
    /*  capital letters, and their contents with small letters. Moreover: */
    /*  .IN.   means belonging to; */
    /*  .INT.  stands for intersection. */


    /*  The code has the following structure: */

    /*         DBVSSC */
    /*              DBVC */
    /*                   DSTINF */
    /*                        DMSK1 */
    /*                        DMSK2 */
    /*                             DPRJ0 */
    /*                             DPRJ1 */
    /*                        DTDC */
    /*                             DMNMOD */
    /*                             DMDIAN */
    /*                   DALG3 */
    /*                        DPRJ0 */
    /*                        DALG1 */
    /*                             DPRJ1 */
    /*                        DINTRS */
    /*                        DTST */
    /*                        DFPSVF */
    /*                             DSL */
    /*                             DALG1D */
    /*                                  DPRJ1 */
    /*                        DAL2 */
    /*                             DPRJ2 */
    /*                        DAL2DP */
    /*                             DMNIND */
    /*                             DPRJ2 */
    /*                             DSL */
    /*                        DSCDRC */


    /*  CALLING SEQUENCE: */

    /*       CALL DBVSSC (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,EPS,BETA, */
    /*    *               BETAI,RHO,RHOI,KMAX,MAXSTP,ERRC,D,D2,DIAGN, */
    /*    *               WORK,NWORK) */


    /*  INPUT PARAMETERS: */

    /*  X       : floating array, of bounds 0:NP, containing the data */
    /*            abscissas  x(i), i=0,1,...,np. */
    /*            Restriction: x(i).LT.x(i+1), i=0,1,...,np. */
    /*  Y       : floating array, of bounds 0:NP, containing the data */
    /*            ordinates  y(i), i=0,1,...,np. */
    /*  NP      : long variable, defining the number of interpolation */
    /*            points. Restriction: np.GE.2 . */
    /*  N       : long variable, containing the degree of s. */
    /*            Restriction: n.GE.3 . */
    /*  K       : long variable, containing the class of continuity of s.
    */
    /*            Restriction:  k.EQ.1  or  k.EQ.2  and  n.GE.3*k . */
    /*  OPT     : long variable, containing a control parameter. It is */
    /*            a three-digit decimal of the form  pqr  (that is of */
    /*            numerical value  p*100+q*10+r ) where: */
    /*            r  controls the constraints on the shape. */
    /*            q  controls the boundary conditions and */
    /*            p  controls the computation of the derivatives, */
    /*            More specifically: */
    /*            r=0 (opt=pq0) : no constraint on the shape is required; */
    /*            r=1 (opt=pq1) : monotonicity constraints are required; */
    /*            r=2 (opt=pq2) : convexity constraints are required; */
    /*            r=3 (opt=pq3) : monotonicity and convexity constraints are
    */
    /*                            required; */
    /*            r=4 (opt=pq4) : local constraints for any subinterval are */
    /*                            supplied by the user (see the description */
    /*                            of the array CONSTR); */
    /*            q=1 (opt=p1r) : no boundary condition is imposed; */
    /*            q=2 (opt=p2r) : non-separable boundary conditions are */
    /*                            imposed (see the description of BETA, */
    /*                            BETAI, RHO, RHOI); */
    /*            q=3 (opt=p3r) : separable boundary conditions are imposed */
    /*                            (see the description of D0, DNP, D20, */
    /*                             D2NP); */
    /*            p=1 (opt=1qr) : the sequence of first derivatives */
    /*                            d(0),....,d(np)  is computed using the */
    /*                            constraints only using subroutine DAL2; */
    /*            p=2 (opt=2qr) : the sequence is computed as the constrained
    */
    /*                            best approximation to Bessel derivatives */
    /*                            using subroutine DAL2DP; */
    /*            p=3 (opt=3qr) : the sequence is computed as the constrained
    */
    /*                            best approximation to a set of third order
    */
    /*                            accurate derivative estimates produced in */
    /*                            subroutine DTDC using subroutine DAL2DP */
    /*                            (since this estimates are inherently mono-
    */
    /*                            tonicity preserving, it is not recommended
    */
    /*                            to associate this option with the convexity
    */
    /*                            constraints only); */
    /*            p=4 (opt=4qr) : the sequence is computed as the constrained
    */
    /*                            best approximation to a set of values given
    */
    /*                            in input by the user using DAL2DP; note */
    /*                            that opt.EQ.410 will provide the classical
    */
    /*                            C(k) function interpolating the data and */
    /*                            the derivatives. */
    /*         Restriction: ( p.GE.1 .AND. p.LE.4 ) .AND. */
    /*                      ( q.GE.1.AND. q.LE.3 ) .AND. */
    /*                      ( r.GE.0 .AND. r.LE.4 ) .AND. */
    /*                      .NOT. ( r.EQ.0 .AND. p.EQ.1 ) . */
    /*  D0      : floating variable containing the left separable boundary */
    /*            condition for the first derivative (d(0)=d0). */
    /*            D0 is referenced only when  q=3 . */
    /*  DNP     : floating variable containing the right separable boundary */
    /*            condition for the first derivative (d(np)=dnp). */
    /*            DNP is referenced only when  q=3 . */
    /*  D20     : floating variable containing the left separable boundary */
    /*            condition for the second derivative (d2(0)=d20). */
    /*            D20 is referenced only when  q=3  and  k=2 . */
    /*  D2NP    : floating variable containing the right separable boundary */
    /*            condition for the second derivative (d2(np)=d2np). */
    /*            D2NP is referenced only when  q=3  and  k=2 . */
    /*  EPS     : floating variable, containing a value for computing the */
    /*            relative tolerance of the method. It should be set greater
    */
    /*            or equal to the machine precision. However, if eps.LE.0, */
    /*            DBVSSC resets it to 0.0001 which has turned out to be a */
    /*            good choice for graphical applications. */
    /*  BETA    : user supplied function, which represents non-separable */
    /*            boundary conditions for the first derivatives. */
    /*            BETA is referenced only when  q=2 . */
    /*  BETAI   : user supplied function, which is the inverse of BETA. */
    /*            BETAI is referenced only when  q=2 . */
    /*  RHO     : user supplied function, which represents non-separable */
    /*            boundary conditions for the second derivatives. */
    /*            RHO is referenced only when  q=2  and  k=2 . */
    /*  RHOI    : user supplied function, which is the inverse of RHO. */
    /*            RHOI is referenced only when  q=2  and  k=2 . */
    /*  KMAX    : long variable, containing the number of iterations */
    /*            allowed for selecting the minimal set ASTAR described */
    /*            below. If kmax.LE.0, DBVSSC resets it to 10 . */
    /*            KMAX is referenced only when  q=2 . */
    /*  MAXSTP  : long variable, containing the number of steps allowed */
    /*            to find dstar in the set of admissible values. */
    /*            If maxstp.LE.0, DBVSSC resets it to 10 . */
    /*            MAXSTP is referenced only when  q=2 . */


    /*  INPUT / OUTPUT PARAMETERS: */

    /*  CONSTR  : long array, of bounds  0:NP , containing, in input the */
    /*            desired constraints on the shape for each subinterval. */
    /*            constr(i)=kind , kind=0,1,2,3 , means that none, monotoni-
    */
    /*            city, convexity, monotonicity and convexity constraint is */
    /*            imposed on the subinterval [x(i),x(i+1)]. If constr(i) is */
    /*            not compatible with the data it is relaxed according to */
    /*            their shape (see subroutine DMSK1 for details). So, on out-
    */
    /*            put, CONSTR contains the constraints actually imposed. */
    /*            For example, if the data are convex and on input we have */
    /*            constr(i)=3 , the result in output will be  constr(i)=2. */
    /*            Restriction: constr(i).GE.0 .AND. constr(i).LE.3 . */
    /*            CONSTR is referenced only when  r=4 . */
    /*  D       : floating array, of bounds 0:NP, containing the first */
    /*            derivatives at  x(i), i=0,1,...,np . If  p=4 , d(i) is the
    */
    /*            input value to be approximated by the computed derivative,
    */
    /*            which is then stored in the same location. */
    /*            On output, D is computed by the routine DAL2 if  p=1  and */
    /*            is computed by the routine DAL2DP if  p=2  or  p=3 . */


    /*  OUTPUT PARAMETERS */

    /*  ERRC    : long variable, containing an error flag which displays */
    /*            the status of the code. The status is divided into: severe
    */
    /*            error (error on the input data, no computation has been */
    /*            done), error (some computation has been done and some */
    /*            information or suggestions are available), warning (some */
    /*            requirement is not fulfilled, but the spline's parameters */
    /*            have been computed), success. */
    /*            errc=0 : success, normal return of the code; */
    /*            errc=1 : severe error, incorrect assignment for some of */
    /*                     the values nwork, opt, np; */
    /*            errc=2 : severe error, for some i the restriction */
    /*                     0.LE.constr(i) .AND. constr(i).LE.3  is not */
    /*                     fulfilled; */
    /*            errc=3 : severe error, incorrect assignment for some of */
    /*                     the values n,k; */
    /*            errc=4 : severe error, the restriction x(i).LT.x(i+1) is */
    /*                     not fulfilled for some i; */
    /*            errc=5 : error, the problem does not have any solution */
    /*                     because the set */
    /*                     betai ( phi(a(0,k)) .INT. beta(a(0,k)) ) */
    /*                     is empty for some k. In other words the boundary */
    /*                     conditions cannot be satisfied and the output */
    /*                     parameters are meaningless. */
    /*                     The user is suggested to increase the value of n.
    */
    /*            errc=6 : warning; for some i, the constraints on the */
    /*                     interval  [x(i),x(i+1)]  are too strong and they */
    /*                     have not been considered. There is no guarantee */
    /*                     that the spline is shape-preserving within all */
    /*                     the intervals. More accurate diagnostic details */
    /*                     can be found in the array DIAGN. */
    /*                     The user is suggested to increase the value of n.
    */
    /*            errc=7 : error, dstar such that beta(dstar).IN.phi(dstar) */
    /*                     has not been found. The long parameter maxstp */
    /*                     should be increased. */
    /*                     The output parameters are meaningless. */
    /*            errc=8 : error, both situations described in errc=6 and */
    /*                     errc=7  have occurred. */
    /*            errc=9 : warning, one of the separable boundary conditions
    */
    /*                     d(0)=d0  and/or  d(np)=dnp  are not compatible */
    /*                     with the constraints in  [x(0),x(1)]  and/or */
    /*                     [x(np-1),x(np)]  which have consequently been */
    /*                     relaxed. The user is suggested to increase the */
    /*                     value of n. More accurate diagnostic details can */
    /*                     be found in the array DIAGN. */
    /*            errc=10: warning, both situations described for errc=6 and
    */
    /*                     errc=9 have occurred. */
    /*            errc=11: warning, one of the separable boundary conditions
    */
    /*                     d2(0)=d20  and/or  d2(np)=d2np  are not compatible
    */
    /*                     with the constraints in  [x(0),x(1)]  and/or */
    /*                     [x(np-1),x(np)] . The boundary conditions have */
    /*                     consequently been approximated. The user is */
    /*                     suggested to increase the value of n. */
    /*            errc=12: warning, both situations described for errc=6 and
    */
    /*                     errc=11 have occurred. */
    /*            errc=13: warning, both situations described for errc=9 and
    */
    /*                     errc=11 have occurred. */
    /*            errc=14: warning, both situations described for errc=10 and
    */
    /*                     errc=11 have occurred. */
    /*            errc=15: warning, the non-separable boundary conditions */
    /*                     d2(np)=rho(d2(0))  are not compatible with the */
    /*                     constraints. The boundary conditions have */
    /*                     consequently been approximated. The user is */
    /*                     suggested to increase the value of n. */
    /*            errc=16: warning, both situations described for errc=6 and
    */
    /*                     errc=15 have occurred. */
    /*            errc=17: warning, both situations described for errc=9 and
    */
    /*                     errc=15 have occurred. */
    /*            errc=18: warning, both situations described for errc=10 and
    */
    /*                     errc=15 have occurred. */
    /*  D2      : floating array of bounds 0:NP containing the second */
    /*            derivatives at knots. D2 is computed in subroutine DCDERC .
    */
    /*            D2 is referenced only when  k=2 . */
    /*  DIAGN   : long array of bounds 0:NP-1 containing further */
    /*            diagnostic information: */
    /*            diagn(i)=0 : the constraints in the interval [x(i),x(i+1)]
    */
    /*                         have been satisfied; */
    /*            diagn(i)=1 : the constraints in the interval [x(i),x(i+1)]
    */
    /*                         have not been satisfied; */



    /*  OTHER PARAMETERS: */

    /*  WORK    : floating array, of bounds 1:NKORK, which is used as */
    /*            a work area to store intermediate results. */
    /*            The same array can be used to provide workspace for both */
    /*            the main subroutines  DBVSSC and DBVSSE . */
    /*  NWORK   : long variable containing the size of the work area. */
    /*            Restriction: nwork .GE. comm+(part+7)*np+(n*(n+11))/2+9 */
    /*                           that is */
    /*                         nwork .GE. 5+(2+7)*np+(n*(n+11))/2+9 */


    /*  ------------------------------------------------ */

    /*  METHOD: */

    /*  Let the longs n and k, such that k=1,2 ; n >= 3k , and the */
    /*  sequence of points  (x(i),y(i)), i=0,1,...,np , with */
    /*  x(0) < x(1) < ... <x(np) , be given; let us denote with  BS(n;k) */
    /*  the set of the splines s of degree n and continuity k whose second */
    /*  derivative, in the case k=2 , vanishes at the knots. We are */
    /*  interested in the existence and construction, if possible, of a */
    /*  shape-preserving interpolating spline s of BS(n;k) such that */

    /*            s(x(i)) = y(i) , i=0,1,...,np                          (1)
    */

    /*  and optionally subject to constraints on the shape in each interval */
    /*  [x(i),x(i+1)] . */

    /*  In the case k=2 the zero derivatives of the spline  s.IN.BS(n;k) are
    */
    /*  then modified to assume non-zero values which are not in contrast */
    /*  with the shape constraints and, if possible, satisfy the boundary */
    /*  conditions eventually imposed by the user. For details we refer to */
    /*  the comments in subroutine DCSDRC. */

    /*  Since any s.IN.BS(n;k) is determined by its values and slopes at */
    /*  x(i) , i=0,1,...,np , we can reformulate the problem as follows: */
    /*  compute the values  d(i), i=0,1,...,np , such that the spline s, */
    /*  satisfying (1) and */

    /*            Ds(x(i)) := d(i) , i=0,1,...,np                        (2)
    */

    /*  is shape-preserving. */
    /*  Setting  delta(i) := (y(i+1)-y(i))/(x(i+1)-x(i)) , we have that s is
    */
    /*  increasing (I) ( decreasing (D) ) in [x(i),x(i+1)] if and only if */
    /*  (d(i),d(i+1))  belongs to */

    /*    D(i) := { (u,v).IN.RxR : u >= 0, v >= 0, v =< -u+ n/k delta(i) } */
    /*                                                                    (3)
    */
    /*  ( D(i) := { (u,v).IN.RxR : u =< 0, v =< 0, v >= -u+ n/k delta(i) } )
    */

    /*  s is convex (CVX) ( concave (CNC) ) if and only if (d(i),d(i+1)) */
    /*  belongs to */

    /*    D(i) := { (u,v).IN.RxR : v >= - (k/(n-k)) u + (n/(n-k)) delta(i) ,
    */
    /*                             v =< - ((n-k)/k) u + (n/k) delta(i) } */
    /*                                                                    (4)
    */
    /*  ( D(i) := { (u,v).IN.RxR : v =< - (k/(n-k)) u + (n/(n-k)) delta(i) ,
    */
    /*                             v >= - ((n-k)/k) u + (n/k) delta(i) }  ) */

    /*  and that s is I (D) and CVX (CNC) if and only if (d(i),d(i+1)) */
    /*  belongs to */

    /*             D(i) := { (u,v) satisfying (3) and (4) } . */

    /*  So, if we choose the family of sets D(i) , i=0,1,...,np-1 , according
    */
    /*  to the shape of the data, we have to solve: */

    /*  PROBLEM P1. Does a sequence ( d(0), d(1), ..., d(np) ) such that */
    /*              (d(i),d(i+1)) .IN. D(i) , i=0,1,...,np-1 , exist ? */

    /*  PROBLEM P2. If P1 is feasible, how can a (the best) solution be */
    /*              computed ? */

    /*  Let DPRJ1: RxR -> R and DPRJ2: RxR -> R be, respectively, the */
    /*  projection maps from uv-plane onto the u-axis and v-axis and let us */
    /*  denote with  B(i) := DPRJ1(D(i)) : */

    /*      ALGORITHM A1[B0]. */
    /*        1. Set A(0):=B(0); J:=np. */
    /*        2. For i=1,2,...,np */
    /*           2.1. Set A(i):= DPRJ2( D(i-1).INT.{ A(i-1) x B(i) } ) . */
    /*           2.2. If A(i) is empty, set J:=i and stop. */
    /*        3. Stop. */

    /*  We have the following result: */

    /*  THEOREM 1. P1 has a solution if, and only if, J=np, that is A(i) is */
    /*             not empty , i=0,1,...,np . If ( d(0), d(1), ...,d(np) ) */
    /*             is a solution then  d(i).IN.A(i) , i=0,1,...,np . */

    /*  A solution can be computed with the following algorithm: */

    /*      ALGORITHM A2[A(np),B0]. */
    /*        1. Choose d(np).IN.A(np). */
    /*        2. For i=np-1, np-2, ..., 0 */
    /*           2.1. Choose d(i).IN.DPRJ1( D(i).INT.{ A(i) x { d(i+1) }}). */
    /*        3. Stop. */

    /*  For more theoretical details about A1 and A2 see \1\ , and for */
    /*  practical details see subprograms DALG1, DAL2, DAL2DP. In the latter
    */
    /*  a dynamic programming scheme is used to find the best solution in */
    /*  the feasible set. More specifically, it is possible to compute the */
    /*  values  d(i),i=0,..,np which satisfy the constraints and are as close
    */
    /*  as possible to another sequence which does not satisfy the */
    /*  constraints but is, in some sense, optimal. */

    /*  From a theoretical point of view, algs A1 and A2 give a complete */
    /*  answer to problems P1 and P2. However, it could be pointed out that,
    */
    /*  for practical applications, we would like to have the best possible */
    /*  plot, whether or not P1 has solutions. Let us suppose that the */
    /*  problem is solvable from 0 to j and from j to np, but that alg A1 */
    /*  applied to the whole family of sets  D(i), i=0,1,...,np-1  gives */
    /*  J.eq.j.ne.np ; if we reset  D(j-1) := A(j-1) x B(j) , alg A1 applied
    */
    /*  to this new family of sets will produce J=np . However, it must be */
    /*  recalled that, in this way, we do not consider the constraints in the
    */
    /*  interval [x(j-i),x(j)] and so there is no guarantee that the spline */
    /*  is shape-preserving in this interval. Whenever this fact cannot be */
    /*  accepted it is convenient to rerun the code with a larger value for */
    /*  the degree n , because the domains of constraints enlarge as n */
    /*  increases (see (3) and (4)). */

    /*  It is immediate to see that separable boundary conditions of the form
    */

    /*            d(0) := d0 ; d(np) := dnp */

    /*  can be easily inserted with a reduction of the corresponding */
    /*  admissible sets which does not modify the above theory: */

    /*       D(0) := D(0).INT.{d(0)=d0} ; D(np) := D(np).INT.{d(np)=dnp} */

    /*  In the case k=2 the corresponding conditions  d2(0) = d20 , */
    /*  d2(np) = d2np  are imposed only if not in contrast with the shape of
    */
    /*  the data; otherwise the admissible values for  d2(0) and d2(np) */
    /*  respectively closest to d20 and d2np are chosen. */

    /*  Now, let beta be a continuous function from R to R, with continuous */
    /*  inverse betai, we want to solve the following non-separable boundary
    */
    /*  valued problem: */

    /*  PROBLEM P3. Do sequences ( d(0), d(1), ..., d(np) ) , such that */
    /*              (d(i),d(i+1)).IN.D(i), i=0,1,...,np-1    and */
    /*              d(np) = beta ( d(0) ) , exist ? */

    /*  It is obvious that a solution of this new problem, if it exists, can
    */
    /*  be found amongst the solutions of P1. Let A(0), A(1),...,A(np) be the
    */
    /*  sequence of sets given by alg A1 (we assume that A(i) is not empty, */
    /*  i=0,1,...,np , that is P1 is solvable or, if this is not the case, */
    /*  the constraints have been relaxed ), we can assume that */
    /*  A(np) = phi(A(0)) , where  phi: R -> R is a set valued function */
    /*  (see \1\ for details). It can be demonstrated that: */

    /*  THEOREM 2. P1 is solvable if, and only if, there is  dstar.IN.A(0) */
    /*             such that   beta(dstar).IN.phi({dstar}) . */

    /*  It should be noted that if ( d(0), d(1), ..., d(np) ) satisfies P1, */
    /*       d(0) .IN. betai(phi(A(0)).INT.beta(A(0))) =: gamma(A(0)) */
    /*  and, consequently, the set of admissible values is reduced. If we */
    /*  repeat this procedure, we get a gradually diminishing admissible set
    */
    /*  for d(0). We define */
    /*     ASTAR := lim A(0,m)  where */
    /*     A(0,0) := A(0)   and   A(0,m) := gamma(A(0,m-1)) ; */
    /*  ASTAR is the minimal admissible set for dstar. We can now combine the
    */
    /*  various theorems and algorithms and give the general algorithm to */
    /*  solve P3: */

    /*      ALGORITHM A3. */
    /*        1. Set A(0,0) := B0 ; m:=1. */
    /*        2. Use A1[A(0,0)] for computing phi (A(0,0)). */
    /*        3. Set A(0,1) := gamma(A(0,0)) */
    /*                       = betai(phi(A(0,0)).INT.beta(A(0,0))). */
    /*        4. If A(0,1) is empty, stop (P1 is unsolvable). */
    /*        5. While ( convergence test not satisfied ) do */
    /*           5.1. Use A1[A(0,m)] for computing A(np,m) = phi (A(0,m)). */
    /*           5.2. Set A(0,m+1) := gamma(A(0,m)). */
    /*           5.3. Set m:=m+1. */
    /*        6. Set ASTAR := A(0,m). */
    /*        7. Use A1[{d(0)}] to find dstar.IN.ASTAR such that */
    /*           beta(dstar).IN.phi(dstar). */
    /*        8. Use A2[beta(dstar),dstar] for computing a sequence */
    /*           ( d(0), d(1), ..., d(np) )  which solves P1. */

    /*  In the case k=2 the corresponding condition  d2(np) = beta2(d2(0)) */
    /*  is imposed only if not in contrast with the shape of */
    /*  the data; otherwise the admissible values for  d2(0) and d2(np) */
    /*  closest to the boundary condition are chosen. */

    /*  References */

    /*  \1\ P.Costantini: Boundary Valued Shape-Preserving Interpolating */
    /*      Splines, ACM Trans. on Math. Softw., companion paper. */
    /*  \2\ R.Bellman, S.Dreyfus: Applied Dynamic Programming, Princeton */
    /*      University Press, New York, 1962. */
    /*  \3\ H.T.Huynh: Accurate Monotone Cubic Interpolation, SIAM J. Num. */
    /*      Anal., 30 (1993), 57-100. */

    /*  The ideas involved in Algorithm A3 have been implemented in the code
    */
    /*  in a general form. Since Algorithm A3 resembles closely the abstract
    */
    /*  formulation it could, therefore, be used for several practical */
    /*  problems. The particular case actually treated is reflected in the */
    /*  contents of the information array INFO (see its description in */
    /*  subroutine DSTINF) which contains all the data needed for the */
    /*  construction of the operators DPRJ0, DPRJ1 and DPRJ2. */

    /*  As a consequence, the user has the following options: */

    /*  - to compute a Spline subject to: */
    /*        - no constraint; */
    /*        - monotonicity constraints, */
    /*        - convexity constraints, */
    /*        - monotonicity and convexity constraints, */
    /*        - one of the above constraints in each subinterval, as */
    /*          specified in the corresponding array CONSTR; */

    /*  - to impose separable or non-separable boundary conditions on the */
    /*    spline. In the latter case, the external functions BETA, BETAI, */
    /*    RHO and RHOI must be supplied, */

    /*  - to assign the first derivatives d(i), i=0,1,...,np , in input or to
    */
    /*    compute them from the only constraints or as the best approximation
    */
    /*    to a set of optimal values. Although the final sequence of */
    /*    derivatives does in any case satisfy the imposed restrictions on */
    /*    the shape, the resulting graphs may exhibit different behaviors. */

    /*  See the description of the input parameter OPT for more details. */
    /*  ------------------------------------------------ */
    /*            End of comments. */
    /*  ------------------------------------------------ */
    /*  COMM contains the number of global data referring to the initial */
    /*  points  (x(i),y(i)) stored in the array INFO, described in */
    /*  subroutine DSTINF. */
    /*  PART contains the number of particular data referring to each */
    /*  interval  (x(i),x(i+1)) , i=0,1,...,np , stored in the array INFO. */
    /*  Assign the success value to the error flag. */
    /* Parameter adjustments */
    --work;

    /* Function Body */
    *errc = 0;
    /*  Check the size of the work area. */
    if (*nwork < *np * 9 + 5 + *n * (*n + 11) / 2 + 9)
    {
        *errc = 1;
        return 0;
    }
    /*  Compute indices for the splitting of the work array WORK. */
    i1 = 1;
    i2 = i1 + ((*np << 1) + 5 + *np + 1);
    i3 = i2 + ((*np + 1) << 1);
    i4 = i3 + ((*np + 1) << 1);
    i5 = i4 + (*n * (*n + 1) / 2 + *n);
    i6 = i5 + (*n + 1);
    i7 = i6 + (*n + 1);
    i8 = i7 + (*n + 1);
    i9 = i8 + (*n + 1);
    i10 = i9 + *np - 1;
    /*  DBVSSC is essentially an interfacing routine which relieves the */
    /*  user of a longer calling sequence. The structure of the method can */
    /*  be seen in DBVC and in the subroutines called. */
    dbvc_ (x, y, np, n, k, opt, d0, dnp, d20, d2np, constr, &work[i1], &c__5,
           &c__2, eps, kmax, maxstp, &work[i2], &work[i3], &work[i9],
           &work[i10], beta, betai, rho, rhoi, d, d2, errc, diagn);
    return 0;
}    /* dbvssc_ */

/*  --------------------------------------------------------------------- */
/* Subroutine */
/*int
dbvsse_ (x, y, np, n, k, xtab, ntab, sbopt, y0opt, y1opt,
  y2opt, errc, d, d2, y0tab, y1tab, y2tab, erre, work, nwork)
     MYREAL *x, *y;
     long *np, *n, *k;
     MYREAL *xtab;
     long *ntab, *sbopt, *y0opt, *y1opt, *y2opt, *errc;
     MYREAL *d, *d2, *y0tab, *y1tab, *y2tab;
     long *erre;
     MYREAL *work;
     long *nwork;
     */
int
dbvsse_ (MYREAL *x, MYREAL *y, long *np, long *n, long *k, MYREAL *xtab,
         long *ntab, long *sbopt, long *y0opt, long *y1opt, long *y2opt,
         long *errc, MYREAL *d, MYREAL *d2, MYREAL *y0tab, MYREAL *y1tab,
         MYREAL *y2tab, long *erre, MYREAL *work, long *nwork)
{
  //    extern /* Subroutine */ int dbve_ ();
    /*static */
    long i1, i2, i3, i4, i5, i6, i7, i8;//, i9;
    //, i10;

    /*  ------------------------------------------------- */
    /*            Lines 621-754 are comment lines. */
    /*            the spline at the tabulation points xtab(i) , */
    /*            i=0,1,...,ntab when the option  y0opt=1  is activated. */
    /*  Y1TAB   : floating array, of bounds 0:NTAB, containing the values of
    */
    /*            the first derivative of the spline at the tabulation points
    */
    /*            xtab(i) , i=0,1,...ntab , when the option y1opt=1 is */
    /*            activated. */
    /*  Y2TAB   : floating array, of bounds 0:NTAB, containing the values of
    */
    /*            the second derivative of the spline at the tabulation */
    /*            points xtab(i) , i=0,1,...,ntab , when the option y2opt=1 */
    /*            is activated. */
    /*  ERRE    : long variable, containing an error flag which displays */
    /*            the status of the code. DBVSSE has only two levels of error
    */
    /*            (see DBVSSC for comparison): success and severe error, */
    /*            which means that some incorrect assignment for input data */
    /*            have been set. */
    /*            ERRE=0:  success, normal return of the code; */
    /*            ERRE=1:  severe error, the value errc gives a status of */
    /*                     error, which means that the output of DBVSSC is */
    /*                     meaningless. Check the input parameters of DBVSSC.
    */
    /*            ERRE=2:  severe error, incorrect assignment for some of */
    /*                     the values ntab, sbopt, y0opt, y1opt, y2opt , */
    /*                     nwork; */
    /*            ERRE=3:  severe error, the restriction xtab(i).LT.xtab(i+1)
    */
    /*                     is not fulfilled for some i when sequential search
    */
    /*                     is required; */


    /*  OTHER PARAMETERS: */

    /*  WORK    : floating array, of bounds 1:NKORK, which is used as */
    /*            a work area to store intermediate results. */
    /*            The same array can be used to provide workspace for both */
    /*            the main subroutines  DBVSSC and DBVSSE . */
    /*  NWORK   : long variable containing the size of the work area. */
    /*            Restriction: nwork .GE. comm+(part+7)*np+(n*(n+11))/2+9 */
    /*                           that is */
    /*                         nwork .GE. 3+(2+7)*np+(n*(n+11))/2+9 */
    /*  ------------------------------------------------- */
    /*            End of comments. */
    /*  ------------------------------------------------- */
    /*  Assign the success value to the error flag. */
    /* Parameter adjustments */
    --work;

    /* Function Body */
    *erre = 0;
    /*  Check the size of the work area. */
    if (*nwork < *np * 9 + 5 + *n * (*n + 11) / 2 + 9)
    {
        *erre = 2;
        return 0;
    }
    /*  Compute indices for the splitting of the work array WORK. */
    i1 = 1;
    i2 = i1 + ((*np << 1) + 5 + *np + 1);
    i3 = i2 + ((*np + 1) << 1);
    i4 = i3 + ((*np + 1) << 1);
    i5 = i4 + (*n * (*n + 1) / 2 + *n);
    i6 = i5 + (*n + 1);
    i7 = i6 + (*n + 1);
    i8 = i7 + (*n + 1);
    //i9 = i8 + (*n + 1);
   //xcode  i10 = i9 + *np - 1;
    /*  DBVSSE is essentially an interfacing routine which relieves the */
    /*  user of a longer calling sequence. The structure of the method can */
    /*  be seen in DBVE and in the subroutines called. */
    dbve_ (x, y, np, n, k, xtab, ntab, sbopt, y0opt, y1opt, y2opt, d, d2, errc,
           &work[i4], &work[i5], &work[i6], &work[i7], &work[i8], y0tab, y1tab,
           y2tab, erre);
    return 0;
}    /* dbvsse_ */


/*  --------------------------------------------------------------------- */
/* Subroutine */ int
dalg1d_ (
	 MYREAL *dstar, MYREAL *a1,
	 long *np,
	 MYREAL *info,
	 long *comm, long *part,
	 MYREAL *eps, MYREAL *a2,
	 long *errc1)
{
    /* System generated locals */
    long i__1;

    /* Local variables */
    //    extern /* Subroutine */ int dprj1_ ();
    /*static */
    long i;

    /*  DALG1D computes the sequence of sets A(i), i=0,1,...,np, implementing
    */
    /*  the algorithm A1[{dstar}], that is with A(0)={dstar} (see the com- */
    /*  ments in subroutine DBVSSC for details). */

    /*  The input parameters NP,COMM,PART,EPS are described in DBVSSC; the */
    /*  input parameter INFO is described in DSTINF; the input parameters A1
    */
    /*  and A2 are described in subprogram DALG1. */

    /*  Item of possible interest is: */

    /*  ERRC1  : Integer parameter, containing a control variable which is */
    /*           then used in subr. DFPSVF */
    /*           errc1 = 0 - success, normal return of the subprogram; */
    /*           errc1 = 1 - A(i) is empty for some i. */
    /* Parameter adjustments */
    --a2;
    --a1;
    --info;

    /* Function Body */
    *errc1 = 0;
    /*  Step 1. */
    a2[1] = *dstar;
    a2[2] = *dstar;
    /*  Step 2. */
    i__1 = *np;
    for (i = 1; i <= i__1; ++i)
    {
        dprj1_ (&a2[((i - 1) << 1) + 1], &a2[((i - 1) << 1) + 2],
                &a1[(i << 1) + 1], &a1[(i << 1) + 2], &i, &info[1], comm, part,
                np, &a2[(i << 1) + 1], &a2[(i << 1) + 2]);
        if (a2[(i << 1) + 1] > a2[(i << 1) + 2] + *eps)
        {
            *errc1 = 1;
            return 0;
        }
        /* L10: */
    }
    return 0;
}    /* dalg1d_ */


/*  --------------------------------------------------------------------- */
/* Subroutine */ int
dal2_ (
       MYREAL *a2,
       long *np,
       MYREAL *info,
       long *comm, long *part,
       MYREAL *d)
{
    /* Initialized data */

    MYREAL fl1d2 = .5;

    //    extern /* Subroutine */ int dprj2_ ();
    long i;
    MYREAL p1, p2;

    /*  DAL2 computes a sequence of slopes (d(0),d(1),...,d(np)) implementing
    */
    /*  alg. A2  described in subr. DBVSSC. Each d(i),i=0,1,...,np , is */
    /*  chosen as the midpoint of the interval of all feasible values . */

    /*  The input parameters NP,COMM,PART and the output parameter D are */
    /*  described in DBVSSC; the input parameter INFO is described in DSTINF.
    */

    /*  Item of possible interest is: */

    /*  A2   : floating array, of bounds 1:2, 0:NP; [a2(1,i),a2(2,i)] */
    /*         is the feasible interval for d(i) . */
    /* Parameter adjustments */
    --a2;
    --info;

    /* Function Body */
    d[*np] = (a2[(*np << 1) + 1] + a2[(*np << 1) + 2]) * fl1d2;
    for (i = *np; i >= 1; --i)
    {
        dprj2_ (&a2[((i - 1) << 1) + 1], &a2[((i - 1) << 1) + 2], &d[i], &d[i],
                &i, &info[1], comm, part, np, &p1, &p2);
        d[i - 1] = (p1 + p2) * fl1d2;
        /* L10: */
    }
    return 0;
}    /* dal2_ */

/*  --------------------------------------------------------------------- */
/* Subroutine */ int
dal2dp_ (
	 MYREAL *a2,
	 long *np,
	 MYREAL *info,
	 long *comm, long *part,
	 MYREAL *d)
{
    /* Initialized data */

    MYREAL fl0 = 0.;
    MYREAL fle30 = 1e30;

    /* System generated locals */
    long i__1;
    MYREAL d__1, d__2;

    /* Local variables */
    //    extern /* Subroutine */ int dprj2_ ();
    MYREAL part0[22], part1[22];
    long i, j;
    MYREAL d0, d1, h0, h1, p1, p2, psi1mn;
    //extern long dmnind_ ();
    long jd0, ind;
    //    extern MYREAL dsl_ ();
    MYREAL psi0[22], psi1[22];

    /*  DAL2DP links algorithm A2 and a dynamic programming scheme */
    /*  to select, among the set of all feasible solutions, the sequence */
    /*  ( d(0),d(1), ..., d(np) ) which is the best 2-norm approximation to */
    /*  a set of optimal values. More precisely, if (ds(0),ds(1), ...,ds(np))
    */
    /*  is the sequence of optimal values, DAL2DP use the following dynamic */
    /*  programming relations */

    /*    psi(0;d(0)) := (d(0)-ds(0))**2 */
    /*    psi(j;d(j)) := (d(j)-ds(j))**2 + MIN(psi(j-1;d(j-1))) */

    /*  for describing the objective function */

    /*      SUM  ((d(j) - ds(j)) ** 2 */
    /*    j=0,np */

    /*  For a complete comprehension of the algorithm see the book \2\ */
    /*  quoted in the references of subr. DBVSSC */

    /*  The input parameters NP,COMM,PART and the output parameter D are */
    /*  described in subprogram DBVSSC; the input parameter INFO is described
    */
    /*  in subprogram DSTINF and the input parameter A2 is described in DAL2.
    */
    /*  The constant NSUBD defined below is related to the discretization of
    */
    /*  the admissible domain. */
    /* Parameter adjustments */
    --a2;
    --info;

    /* Function Body */
    ind = *comm + *part * *np + 1;
    /* Computing MAX */
    d__1 = fl0, d__2 = (a2[2] - a2[1]) / 20;
    h0 = MAX (d__1, d__2);
    for (j = 0; j <= 20; ++j)
    {
        part0[j] = a2[1] + j * h0;
        /* L5: */
    }
    part0[21] = dsl_ (&a2[1], &a2[2], &info[ind]);
    d[0] = dsl_ (&a2[1], &a2[2], &info[ind]);
    for (j = 0; j <= 21; ++j)
    {
        /* Computing 2nd power */
        d__1 = part0[j] - d[0];
        psi0[j] = d__1 * d__1;
        /* L10: */
    }
    i__1 = *np;
    for (i = 1; i <= i__1; ++i)
    {
        /* Computing MAX */
        d__1 = fl0, d__2 = (a2[(i << 1) + 2] - a2[(i << 1) + 1]) / 20;
        h1 = MAX (d__1, d__2);
        for (j = 0; j <= 20; ++j)
        {
            part1[j] = a2[(i << 1) + 1] + j * h1;
            /* L15: */
        }
        part1[21] = dsl_ (&a2[(i << 1) + 1], &a2[(i << 1) + 2], &info[ind + i]);
        psi1mn = fle30;
        for (j = 0; j <= 21; ++j)
        {
            d1 = part1[j];
            dprj2_ (&a2[((i - 1) << 1) + 1], &a2[((i - 1) << 1) + 2], &d1, &d1,
                    &i, &info[1], comm, part, np, &p1, &p2);
            d0 = dsl_ (&p1, &p2, &d[i - 1]);
            if (h0 > fl0)
            {
                jd0 = dmnind_ (&d0, part0);
            }
            else
            {
                jd0 = 0;
            }
            /* Computing 2nd power */
            d__1 = d1 - info[ind + i];
            psi1[j] = d__1 * d__1 + psi0[jd0];
            if (psi1[j] < psi1mn)
            {
                psi1mn = psi1[j];
                d[i] = d1;
            }
            /* L20: */
        }
        h0 = h1;
        for (j = 0; j <= 21; ++j)
        {
            psi0[j] = psi1[j];
            part0[j] = part1[j];
            /* L30: */
        }
        /* L40: */
    }
    for (i = *np; i >= 1; --i)
    {
        dprj2_ (&a2[((i - 1) << 1) + 1], &a2[((i - 1) << 1) + 2], &d[i], &d[i],
                &i, &info[1], comm, part, np, &p1, &p2);
        d[i - 1] = dsl_ (&p1, &p2, &d[i - 1]);
        /* L50: */
    }
    return 0;
}    /* dal2dp_ */

/*  --------------------------------------------------------------------- */
MYREAL
dbl_ (
      MYREAL *x,
      long *n,
      MYREAL *l, MYREAL *x0, MYREAL *xn, MYREAL *tb,
      long *flag__,
      MYREAL *laux0)
{
    /* Initialized data */

    MYREAL fl1 = 1.;

    /* System generated locals */
    long i__1;
    MYREAL ret_val, d__1;

    /* Builtin functions */


    /* Local variables */
    MYREAL xnmx;
    long i;
    MYREAL aux, xmx0;

    /*  DBL computes the value assumed by the n-degree Bernstein polynomial */
    /*  of a function  l  in the interval  (x0,xn)  at the point  x . */
    /*  The evaluation is made using a Horner scheme, and the instructions */
    /*  which do not depend upon  x  are executed under the control of */
    /*  the variable  FLAG , for avoiding useless computations in */
    /*  subsequent calls. */
    /*  The degree  n  is supposed greater or equal to  3 . */


    /*  INPUT PARAMETERS */

    /*  X     : floating variable, containing the evaluation point. */
    /*  N     : long variable, containing the degree of Bernstein */
    /*          polynomial. */
    /*  L     : floating array, of bounds  0:N , containing the values */
    /*          of the function  l . */
    /*  X0    : floating variable, containing the left extreme of the */
    /*          interval. */
    /*  XN    : floating variable, containing the right extreme of the */
    /*          interval. */
    /*  TB    : floating array, of bounds  0:N , containing the binomial */
    /*          terms used for computing the Bernstein polynomial. */
    /*  FLAG  : long variable, containing a control parameter. */
    /*          In the case  flag=0  DBL  assumes to perform the first */
    /*          evaluation of the polynomial, and computes the values */
    /*          tb(i)*l(i) , i=0,1,...,n . In the case  flag=1  DBL */
    /*          assumes to perform subsequent evaluations, and uses the */
    /*          values previously computed. */


    /*  OTHER PARAMETERS */

    /*  LAUX0 : floating array, of bounds 0:N used as a work area to store */
    /*          intermediate results. */
    if (*flag__ == 0)
    {
        i__1 = *n;
        for (i = 0; i <= i__1; ++i)
        {
            laux0[i] = tb[i] * l[i];
            /* L10: */
        }
    }
    xnmx = *xn - *x;
    xmx0 = *x - *x0;
    aux = fl1;
    ret_val = laux0[*n];
    for (i = *n - 1; i >= 0; --i)
    {
        aux = xnmx * aux;
        ret_val = laux0[i] * aux + xmx0 * ret_val;
        /* L20: */
    }
    d__1 = *xn - *x0;
    ret_val /= pow_di (&d__1, n);
    return ret_val;
}    /* dbl_ */

/*  --------------------------------------------------------------------- */
MYREAL
dbl1_ (
       MYREAL *x,
       long *n,
       MYREAL *l, MYREAL *x0, MYREAL *xn, MYREAL *tb,
       long *flag__,
       MYREAL *laux1)
{
    /* Initialized data */

    MYREAL fl1 = 1.;

    /* System generated locals */
    long i__1;
    MYREAL ret_val, d__1;

    /* Builtin functions */


    /* Local variables */
    MYREAL xnmx;
    long i;
    MYREAL aux, xmx0;

    /*  DBL1 computes the value assumed by the first derivative of an */
    /*  n-degree Bernstein polynomial of a function  l  in the interval */
    /*  (x0,xn)  at the point  x . */
    /*  The evaluation is made using a Horner scheme, and the instructions */
    /*  which do not depend upon  x  are executed under the control of */
    /*  the variable  FLAG , for avoiding useless computations in */
    /*  subsequent calls. */
    /*  The degree  n  is supposed greater or equal to  3 . */

    /*  INPUT PARAMETERS */

    /*  X     : floating variable, containing the evaluation point. */
    /*  N     : long variable, containing the degree of Bernstein */
    /*          polynomial. */
    /*  L     : floating array, of bounds  0:N , containing the values */
    /*          of the function  l . */
    /*  X0    : floating variable, containing the left extreme of the */
    /*          interval. */
    /*  XN    : floating variable, containing the right extreme of the */
    /*          interval. */
    /*  TB    : floating array, of bounds  0:N-1 , containing the binomial */
    /*          terms used for computing the Bernstein polynomial. */
    /*  FLAG  : long variable, containing a control parameter. */
    /*          In the case  flag=0  DBL1  assumes to perform the first */
    /*          evaluation of the polynomial, and computes the values */
    /*          tb(i)*(l(i+1)-l(i)) , i=0,1,...,n-1 . In the case  flag=1 */
    /*          DBL1 assumes to perform subsequent evaluations, and uses */
    /*          the values previously computed. */


    /*  OTHER PARAMETERS */

    /*  LAUX1 : floating array, of bounds 0:N-1 used as a work area to store
    */
    /*          intermediate results. */
    if (*flag__ == 0)
    {
        i__1 = *n - 1;
        for (i = 0; i <= i__1; ++i)
        {
            laux1[i] = tb[i] * (l[i + 1] - l[i]);
            /* L10: */
        }
    }
    xnmx = *xn - *x;
    xmx0 = *x - *x0;
    aux = fl1;
    ret_val = laux1[*n - 1];
    for (i = *n - 2; i >= 0; --i)
    {
        aux = xnmx * aux;
        ret_val = laux1[i] * aux + xmx0 * ret_val;
        /* L20: */
    }
    d__1 = *xn - *x0;
    ret_val = *n * ret_val / pow_di (&d__1, n);
    return ret_val;
}    /* dbl1_ */

/*  --------------------------------------------------------------------- */
MYREAL
dbl2_ (
       MYREAL *x,
       long *n,
       MYREAL *l, MYREAL *x0, MYREAL *xn, MYREAL *tb,
       long *flag__,
       MYREAL *laux2)
{
    /* Initialized data */

    MYREAL fl1 = 1.;

    /* System generated locals */
    long i__1;
    MYREAL ret_val, d__1;

    /* Builtin functions */


    /* Local variables */
    MYREAL xnmx;
    long i;
    MYREAL aux, xmx0;

    /*  DBL2 computes the value assumed by the second derivative of an */
    /*  n-degree Bernstein polynomial of a function  l  in the interval */
    /*  (x0,xn)  at the point  x . */
    /*  The evaluation is made using a Horner scheme, and the instructions */
    /*  which do not depend upon  x  are executed under the control of */
    /*  the variable  FLAG , for avoiding useless computations in */
    /*  subsequent calls. */
    /*  The degree  n  is supposed greater or equal to  3 . */

    /*  INPUT PARAMETERS */

    /*  X     : floating variable, containing the evaluation point. */
    /*  N     : long variable, containing the degree of Bernstein */
    /*          polynomial. */
    /*  L     : floating array, of bounds  0:N , containing the values */
    /*          of the function  l . */
    /*  X0    : floating variable, containing the left extreme of the */
    /*          interval. */
    /*  XN    : floating variable, containing the right extreme of the */
    /*          interval. */
    /*  TB    : floating array, of bounds  0:N-2 , containing the binomial */
    /*          terms used for computing the Bernstein polynomial. */
    /*  FLAG  : long variable, containing a control parameter. */
    /*          In the case  flag=0  DBL2  assumes to perform the first */
    /*          evaluation of the polynomial, and computes the values */
    /*          tb(i)*(l(i+2)-2*l(i+1)+l(i)) , i=0,1,...,n-2 . */
    /*          In the case  flag=1  DBL2 assumes to perform subsequent */
    /*          evaluations, and uses the values previously computed. */


    /*  OTHER PARAMETERS */

    /*  LAUX2 : floating array, of bounds 0:N-2 used as a work area to store
    */
    /*          intermediate results. */
    if (*flag__ == 0)
    {
        i__1 = *n - 2;
        for (i = 0; i <= i__1; ++i)
        {
            laux2[i] = tb[i] * (l[i + 2] - l[i + 1] - l[i + 1] + l[i]);
            /* L10: */
        }
    }
    xnmx = *xn - *x;
    xmx0 = *x - *x0;
    aux = fl1;
    ret_val = laux2[*n - 2];
    for (i = *n - 3; i >= 0; --i)
    {
        aux = xnmx * aux;
        ret_val = laux2[i] * aux + xmx0 * ret_val;
        /* L20: */
    }
    d__1 = *xn - *x0;
    ret_val = *n * (*n - 1) * ret_val / pow_di (&d__1, n);
    return ret_val;
}    /* dbl2_ */

/*  --------------------------------------------------------------------- */
/* Subroutine */ int
dbntab_ (
	 MYREAL *x, MYREAL *y,
	 long *np,
	 MYREAL *xtab,
	 long *ntab, long *y0opt, long *y1opt, long *y2opt, long *n, long *k,
	 MYREAL *d, MYREAL *d2, MYREAL *tb, MYREAL *l, MYREAL *laux0, MYREAL *laux1, MYREAL *laux2, MYREAL *y0tab,MYREAL  *y1tab, MYREAL *y2tab)
{
    /* System generated locals */
    long i__1;

    /* Local variables */
    long j;
    //    extern /* Subroutine */ int dbsear_ (), dlspis_ ();
    //  extern MYREAL dbl_ ();
    long ind;
    //  extern MYREAL dbl1_ (), dbl2_ ();

    /*  DBNTAB evaluates the spline and/or its first derivative and/or its */
    /*  second derivative at the points  xtab(j) , j=0,1,...,ntab  using */
    /*  a binary search for finding the interval  [x(i),x(i+1)] in which */
    /*  the tabulation point falls. The input (X,Y,NP,XTAB,NTAB,Y0OPT, */
    /*  Y1OPT,Y2OPT,N,K,D,D2,TB) and the output (Y0TAB,Y1TAB,Y2TAB) */
    /*  parameters have been explained in subroutine DBVSSE. For the others */
    /*  see subroutines DTRMB, DLSPIS. */
    /* Parameter adjustments */
    --tb;

    /* Function Body */
    i__1 = *ntab;
    for (j = 0; j <= i__1; ++j)
    {
        /*  Call subprogram  DBSEAR  to compute the index  ind  such that */
        /*       x(ind).LE.xtab(j).LT.x(ind+1) . */
        dbsear_ (x, np, &xtab[j], &ind);
        /*  Call subprogram  DLSPIS  to compute the linear shape-preserving */
        /*  interpolating spline  l  at */
        /*      x(ind)+p*(x(ind+1)-x(ind))/n , p=0,1,...,n . */
        dlspis_ (x, y, d, d2, np, n, k, &ind, l);
        if (*y0opt == 1)
        {
            /*  Evaluate the spline at  xtab(j) . */
            y0tab[j] =
                dbl_ (&xtab[j], n, l, &x[ind], &x[ind + 1],
                      &tb[*n * (*n + 1) / 2], &c__0, laux0);
        }
        if (*y1opt == 1)
        {
            /*  Evaluate the first derivative of the spline at  xtab(j) . */
            y1tab[j] =
                dbl1_ (&xtab[j], n, l, &x[ind], &x[ind + 1],
                       &tb[(*n - 1) * *n / 2], &c__0, laux1);
        }
        if (*y2opt == 1)
        {
            /*  Evaluate the second derivative of the spline at  xtab(j) . */
            y2tab[j] =
                dbl2_ (&xtab[j], n, l, &x[ind], &x[ind + 1],
                       &tb[(*n - 2) * (*n - 1) / 2], &c__0, laux2);
        }
        /* L40: */
    }
    return 0;
}    /* dbntab_ */

/*  --------------------------------------------------------------------- */
/* Subroutine */ int
dbsear_ (
	 MYREAL *x,
	 long *np,
	 MYREAL *xtab,
	 long *ind)
{
    long i1, i2, med;

    /*  Given an ordered set of points  (x(i), i=0,1,...,np)  and the */
    /*  point  xtab , DBSEAR finds the index  ind  such that */

    /*              x(ind) .LE. xtab .LT. x(ind+1) */

    /*  using a standard binary search. DBSEAR  sets  ind=0  or  ind=np-1 */
    /*  whenever  xtab.LT.x(0)  or  x(np).LE.xtab . */


    /*  INPUT PARAMETERS */

    /*  X     : floating array, of bounds  0:NP , containing the set of */
    /*          ordered points. */
    /*  XTAB  : floating variable, containing the point to be processed. */
    /*  NP    : long  variable, defining the number of points of the */
    /*          ordered set. */


    /*  OUTPUT PARAMETERS */

    /*  IND   : long variable, whose value selects the interval in */
    /*          which the point  xtab  falls. */
    if (*xtab <= x[0])
    {
        *ind = 0;
        return 0;
    }
    if (*xtab >= x[*np])
    {
        *ind = *np - 1;
        return 0;
    }
    i1 = 0;
    i2 = *np;
L10:
    if (!(i1 != i2 - 1))
    {
        goto L20;
    }
    med = (i1 + i2) / 2;
    if (*xtab < x[med])
    {
        i2 = med;
    }
    else if (*xtab > x[med])
    {
        i1 = med;
    }
    else
    {
        *ind = med;
        return 0;
    }
    goto L10;
L20:
    *ind = i1;
    return 0;
}    /* dbsear_ */

/*  --------------------------------------------------------------------- */

/*  --------------------------------------------------------------------- */
int
dbve_ (
       MYREAL *x, MYREAL *y,
       long *np, long *n, long *k,
       MYREAL *xtab,
       long *ntab, long *sbopt, long *y0opt, long *y1opt, long *y2opt,
       MYREAL *d, MYREAL *d2,
       long *errc,
       MYREAL *tb, MYREAL *l, MYREAL *laux0, MYREAL *laux1,MYREAL *laux2,MYREAL *y0tab,MYREAL *y1tab,MYREAL *y2tab,
       long *erre)
{
    /* System generated locals */
    long i__1;

    /* Local variables */
    long i;
    //  extern  int dtrmb_ (), dbntab_ (), dsqtab_ ();

    /*  DBVE checks input parameters and evaluates the required spline. */

    /*  The input parameters X,Y,NP,N,K,XTAB,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT, */
    /*  D,D2,ERRC and the output parameters Y0TAB,Y1TAB,Y2TAB,ERRE are */
    /*  described in subroutine DBVSSE. The others are used as work areas */
    /*  and will be eventually described in the subsequent routines. */
    /*  Check the correctness of input data, that is if subroutine DBVSSC */
    /*  has correctly run. */
    /* Parameter adjustments */
    --tb;

    /* Function Body */
    if (*errc == 1 || *errc == 2 || *errc == 3 || *errc == 4 || *errc == 5
            || *errc == 7 || *errc == 8 || *errc == 10)
    {
        *erre = 1;
        return 0;
    }
    /*  Check the input parameters NTAB, SBOPT, Y0OPT, Y1OPT, Y2OPT. */
    if (*ntab < 0 || (*sbopt != 1 && *sbopt != 2)
            || (*y0opt != 0 && *y0opt != 1) || (*y1opt != 0 && *y1opt != 1)
            || (*y2opt != 0 && *y2opt != 1))
    {
        *erre = 2;
        return 0;
    }
    if (*sbopt == 1 && *ntab > 0)
    {
        /*  Check the abscissas of the tabulation points when the sequential
        */
        /*  search is required. */
        i__1 = *ntab - 1;
        for (i = 0; i <= i__1; ++i)
        {
            if (xtab[i + 1] <= xtab[i])
            {
                *erre = 3;
                return 0;
            }
            /* L10: */
        }
    }
    /*  Call subprogram DTRMB to compute the binomial terms needed */
    /*  in the expression of Bernstein polynomials. */
    dtrmb_ (n, &tb[1]);
    if (*sbopt == 1)
    {
        /*  sbopt=1:  sequential search is required. */
        dsqtab_ (x, y, np, xtab, ntab, y0opt, y1opt, y2opt, n, k, d, d2, &tb[1],
                 l, laux0, laux1, laux2, y0tab, y1tab, y2tab);
    }
    else
    {
        /*  sbopt=2: binary search is required. */
        dbntab_ (x, y, np, xtab, ntab, y0opt, y1opt, y2opt, n, k, d, d2, &tb[1],
                 l, laux0, laux1, laux2, y0tab, y1tab, y2tab);
    }
    return 0;
}    /* dbve_ */

/*  --------------------------------------------------------------------- */
int
dfpsvf_ (
	 MYREAL *a1, MYREAL *a2,
	 long *np,
	 MYREAL *info,
	 long *comm, long *part,
	 MYREAL *eps,
	 long *maxstp,
	 MYREAL (*beta) (MYREAL *),
	 long *errc,
	 MYREAL *dstar)
{
    /* Initialized data */

    MYREAL fl1d2 = .5;

    /* System generated locals */
    long i__1, i__2;

    /* Builtin functions */


    /* Local variables */
    long step, errc1;
    MYREAL h;
    long i;
    //  extern  int dalg1d_ ();
    long substp;
    MYREAL mid;
    //    extern MYREAL dsl_ ();

    /*  DFPSVF finds, if possible, dstar.IN.[a1(1,0),a1(2,0)] such that */
    /*               beta(dstar) .INT. phi(dstar)                     (1) */

    /*  The input parameters NP,COMM,PART,EPS,MAXSTP,BETA, and the output */
    /*  parameter ERRC are described in DBVSSC. The input parameters A1 and */
    /*  A2 are described in  DALG1. */
    /* Parameter adjustments */
    --a2;
    --a1;
    --info;

    /* Function Body */
    /*  If the optimum input value of dstar does not belong to */
    /*  [a1(1,0),a1(2,0)] , the nearest extreme of this interval */
    /*  replaces the old value of dstar. */
    *dstar = dsl_ (&a1[1], &a1[2], dstar);
    /*  Compute phi(dstar). */
    dalg1d_ (dstar, &a1[1], np, &info[1], comm, part, eps, &a2[1], &errc1);
    /*  If phi(dstar) is not empty and dstar satisfies (1), it is the desired
    */
    /*  value. */
    if (errc1 == 0
            && (a2[(*np << 1) + 1] - *eps <= (*beta) (dstar)
                && (*beta) (dstar) <= a2[(*np << 1) + 2] + *eps))
    {
        return 0;
    }
    /*  If it is not the case, look for another value. First, check */
    /*  if the midpoint of the interval of all possible values satisfies */
    /*  condition (1). */
    mid = (a1[1] + a1[2]) * fl1d2;
    *dstar = mid;
    dalg1d_ (dstar, &a1[1], np, &info[1], comm, part, eps, &a2[1], &errc1);
    if (errc1 == 0
            && (a2[(*np << 1) + 1] - *eps <= (*beta) (dstar)
                && (*beta) (dstar) <= a2[(*np << 1) + 2] + *eps))
    {
        return 0;
    }
    /*  Second, check if any point of a tabulation of interval */
    /*  [a1(1,0),a1(2,0)]  satisfies the condition (1). The tabulation */
    /*  points are given by jumps of decreasing lenghts with alternate */
    /*  direction with respect to the middle of the interval. */
    h = (a1[2] - a1[1]) * fl1d2;
    i__1 = *maxstp;
    for (step = 1; step <= i__1; ++step)
    {
        h *= fl1d2;
        substp = pow_ii (&c__2, &step) - 1;
        i__2 = substp;
        for (i = 1; i <= i__2; i += 2)
        {
            *dstar = mid - i * h;
            dalg1d_ (dstar, &a1[1], np, &info[1], comm, part, eps, &a2[1],
                     &errc1);
            if (errc1 == 0
                    && (a2[(*np << 1) + 1] - *eps <= (*beta) (dstar)
                        && (*beta) (dstar) <= a2[(*np << 1) + 2] + *eps))
            {
                return 0;
            }
            *dstar = mid + i * h;
            dalg1d_ (dstar, &a1[1], np, &info[1], comm, part, eps, &a2[1],
                     &errc1);
            if (errc1 == 0
                    && (a2[(*np << 1) + 1] - *eps <= (*beta) (dstar)
                        && (*beta) (dstar) <= a2[(*np << 1) + 2] + *eps))
            {
                return 0;
            }
            /* L10: */
        }
        /* L20: */
    }
    /*  Finally, check if condition (1) is satisfied by one of the */
    /*  [a1(1,0),a1(2,0)] extremes. */
    *dstar = a1[1];
    dalg1d_ (dstar, &a1[1], np, &info[1], comm, part, eps, &a2[1], &errc1);
    if (errc1 == 0
            && (a2[(*np << 1) + 1] - *eps <= (*beta) (dstar)
                && (*beta) (dstar) <= a2[(*np << 1) + 2] + *eps))
    {
        return 0;
    }
    *dstar = a1[2];
    dalg1d_ (dstar, &a1[1], np, &info[1], comm, part, eps, &a2[1], &errc1);
    if (errc1 == 0
            && (a2[(*np << 1) + 1] - *eps <= (*beta) (dstar)
                && (*beta) (dstar) <= a2[(*np << 1) + 2] + *eps))
    {
        return 0;
    }
    /*  If dstar satisfying (1) has not been found, send a message resetting
    */
    /*  the error flag errc. */
    if (*errc == 6)
    {
        *errc = 8;
    }
    else
    {
        *errc = 7;
    }
    return 0;
}    /* dfpsvf_ */

/*  --------------------------------------------------------------------- */
int
dintrs_ (
	 MYREAL *a, MYREAL *b, MYREAL *c, MYREAL *d, MYREAL *p1, MYREAL *p2)
{
    /*  DINTRS computes the intersection of the two intervals  [a,b] */
    /*  and [c,d]. [p1,p2] is the result. The output  p1.GT.p2 means that */
    /*  the intervals are disjoint. DINTRS assumes  a.LE.b  and  c.LE.d . */
    *p1 = MAX (*a, *c);
    *p2 = MIN (*b, *d);
    return 0;
}    /* dintrs_ */

/*  --------------------------------------------------------------------- */
int
dlspis_ (
	 MYREAL *x, MYREAL *y, MYREAL *d, MYREAL *d2,
	 long *np, long *n, long *k, long *ind,
	 MYREAL *l)
{
    /* Initialized data */

    MYREAL fl1 = 1.;
    MYREAL fl2 = 2.;
    MYREAL fl4 = 4.;

    /* System generated locals */
    long i__1;
    MYREAL d__1;

    /* Local variables */
    MYREAL h;
    long i;
    MYREAL alpha, q1, q2;

    /*  DLSPIS   evaluates the control points of the Bernstein-Bezier net */
    /*  l:=l(x) of the interpolating spline  s:=s(x) , s.IN.BS(n;k) in the */
    /*  interval  [x(ind),x(ind+1)] . For a description of the function  l */
    /*  see the comments in subroutines  DBVSSC and DSCDRC. Here we only */
    /*  recall that the structure of the net is different for k=1 or k=2 . */

    /*  The input parameters  X,Y,D,D2,NP,N,K  are explained in subroutine */
    /*  SPISE. */

    /*  OTHER PARAMETERS */

    /*  IND   : long variable, used to select the knots interval. */
    /*  L     : floating array, of bounds  0:N , containing the ordinates */
    /*          of the control points. */
    h = x[*ind + 1] - x[*ind];
    /*  Compute the net in the case  k=1 . */
    if (*k == 1)
    {
        l[0] = y[*ind];
        l[1] = y[*ind] + h * d[*ind] / *n;
        l[*n] = y[*ind + 1];
        l[*n - 1] = y[*ind + 1] - h * d[*ind + 1] / *n;
        i__1 = *n - 2;
        for (i = 2; i <= i__1; ++i)
        {
            alpha = (i - fl1) / (*n - fl2);
            l[i] = (fl1 - alpha) * l[1] + alpha * l[*n - 1];
            /* L10: */
        }
    }
    else if (*k == 2)
    {
        /*  Compute the net in the case  k=2 . */
        l[0] = y[*ind];
        l[1] = l[0] + h * d[*ind] / *n;
        /* Computing 2nd power */
        d__1 = h;
        l[2] = d__1 * d__1 * d2[*ind] / (*n * (*n - 1)) + fl2 * l[1] - l[0];
        l[*n] = y[*ind + 1];
        l[*n - 1] = l[*n] - h * d[*ind + 1] / *n;
        /* Computing 2nd power */
        d__1 = h;
        l[*n - 2] =
            d__1 * d__1 * d2[*ind + 1] / (*n * (*n - 1)) + fl2 * l[*n - 1] -
            l[*n];
        alpha = (*n - 4) / 2 / (*n - fl4);
        q1 =
            alpha * (y[*ind + 1] - fl2 * h * d[*ind + 1] / *n) + (fl1 -
                    alpha) *
            (y[*ind] + fl2 * h * d[*ind] / *n);
        q2 =
            (fl1 - alpha) * (y[*ind + 1] - fl2 * h * d[*ind + 1] / *n) +
            alpha * (y[*ind] + fl2 * h * d[*ind] / *n);
        i__1 = *n / 2;
        for (i = 3; i <= i__1; ++i)
        {
            alpha = (i - fl2) / (*n / 2 - fl2);
            l[i] = (fl1 - alpha) * l[2] + alpha * q1;
            /* L30: */
        }
        i__1 = *n - 3;
        for (i = *n / 2 + 1; i <= i__1; ++i)
        {
            alpha = (*n - fl2 - i) / (*n / 2 - fl2);
            l[i] = (fl1 - alpha) * l[*n - 2] + alpha * q2;
            /* L40: */
        }
    }
    return 0;
}    /* dlspis_ */

/*  --------------------------------------------------------------------- */
MYREAL
dmdian_ (MYREAL *a, MYREAL *b, MYREAL *c)
{
    /* System generated locals */
    MYREAL ret_val, d__1, d__2;

    /*  Given three numbers a,b,c , median  is the one which lies between */
    /*  the other two. */
    /* Computing MIN */
    d__1 = MAX (*a, *b), d__2 = MAX (*b, *c), d__1 = MIN (d__1, d__2), d__2 =
                                    MAX (*c, *a);
    ret_val = MIN (d__1, d__2);
    return ret_val;
}    /* dmdian_ */

/*  --------------------------------------------------------------------- */
long
dmnind_ (
	 MYREAL *d, MYREAL *part)
{
    /* Initialized data */

    MYREAL fl30 = 1e30;

    /* System generated locals */
    long ret_val = 0;
    MYREAL d__1;

    /* Local variables */
    long j;
    MYREAL mindis, aux[22];

    /*  DMNIND finds the index of the component of the array PART closest */
    /*  to d . */
    for (j = 0; j <= 21; ++j)
    {
        aux[j] = (d__1 = *d - part[j], fabs (d__1));
        /* L10: */
    }
    mindis = fl30;
    for (j = 0; j <= 21; ++j)
    {
        if (aux[j] < mindis)
        {
            ret_val = j;
            mindis = aux[j];
        }
        /* L20: */
    }
    return ret_val;
}    /* dmnind_ */

/*  --------------------------------------------------------------------- */
MYREAL
dmnmod_ (MYREAL *a, MYREAL *b)
{
    /* Initialized data */

    MYREAL fl1 = 1.;
    MYREAL fl2 = 2.;

    /* System generated locals */
    MYREAL ret_val, d__1, d__2;

    /* Builtin functions */
    //  MYREAL d_sign ();

    /*  Given two real numbers a and b, DMNMOD returns the number between */
    /*  a and b which is closest to zero. */
    /* Computing MIN */
    d__1 = fabs (*a), d__2 = fabs (*b);
    ret_val = (d_sign (&fl1, a) + d_sign (&fl1, b)) / fl2 * MIN (d__1, d__2);
    return ret_val;
}    /* dmnmod_ */

/*  --------------------------------------------------------------------- */
int
dmsk1_ (
	MYREAL *info,
	long *constr, long *comm, long *part, long *ind1, long *np)
{
    /* System generated locals */
    long i__1;

    /* Local variables */
    long i;

    /*  DMSK1 compares the constraints required in input by the user and */
    /*  stored in the array CONSTR with the shape of the data, stored by */
    /*  DSTINF in the array INFO. If the required and real shapes do not */
    /*  agree, DMSK1 resets both INFO and CONSTR with the 'intersection' */
    /*  of the shapes. For example, if info(ind1+i)=11 , that is the data */
    /*  are increasing and convex, and constr(i)=2 , that is only convexity */
    /*  is required, then the output values will be  info(ind1+i)=10 */
    /*  (convexity) and constr(i)=2 (unchanged). If  info(ind1+i)=20 */
    /*  (concavity) and  constr(i)=1 (monotonicity) the output will be */
    /*  info(ind1+i)=constr(i)=0 (no constraints). So, the computations made
    */
    /*  in DALG3 will be based on these new values for selecting the domains
    */
    /*  of admissible derivatives, and CONSTR will contain information on */
    /*  the constraints effectively imposed. */
    /*  Further details on the parameters INFO and IND1 can be found in sub-
    */
    /*  routine DSTINF; CONSTR, COMM, PART, NP are explaained in subroutine */
    /*  DBVSSC. */
    /* Parameter adjustments */
    --info;

    /* Function Body */
    i__1 = *np - 1;
    for (i = 0; i <= i__1; ++i)
    {
        if (info[*ind1 + i] == 0. || constr[i] == 0)
        {
            info[*ind1 + i] = 0.;
            constr[i] = 0;
        }
        else if (info[*ind1 + i] == 1.)
        {
            if (constr[i] == 1)
            {
                info[*ind1 + i] = 1.;
                constr[i] = 1;
            }
            else if (constr[i] == 2)
            {
                info[*ind1 + i] = 0.;
                constr[i] = 0;
            }
            else if (constr[i] == 3)
            {
                info[*ind1 + i] = 1.;
                constr[i] = 1;
            }
        }
        else if (info[*ind1 + i] == 2.)
        {
            if (constr[i] == 1)
            {
                info[*ind1 + i] = 2.;
                constr[i] = 1;
            }
            else if (constr[i] == 2)
            {
                info[*ind1 + i] = 0.;
                constr[i] = 0;
            }
            else if (constr[i] == 3)
            {
                info[*ind1 + i] = 2.;
                constr[i] = 1;
            }
        }
        else if (info[*ind1 + i] == 10.)
        {
            if (constr[i] == 1)
            {
                info[*ind1 + i] = 0.;
                constr[i] = 0;
            }
            else if (constr[i] == 2)
            {
                info[*ind1 + i] = 10.;
                constr[i] = 2;
            }
            else if (constr[i] == 3)
            {
                info[*ind1 + i] = 10.;
                constr[i] = 2;
            }
        }
        else if (info[*ind1 + i] == 20.)
        {
            if (constr[i] == 1)
            {
                info[*ind1 + i] = 0.;
                constr[i] = 0;
            }
            else if (constr[i] == 2)
            {
                info[*ind1 + i] = 20.;
                constr[i] = 2;
            }
            else if (constr[i] == 3)
            {
                info[*ind1 + i] = 20.;
                constr[i] = 2;
            }
        }
        else if (info[*ind1 + i] == 11.)
        {
            if (constr[i] == 1)
            {
                info[*ind1 + i] = 1.;
                constr[i] = 1;
            }
            else if (constr[i] == 2)
            {
                info[*ind1 + i] = 10.;
                constr[i] = 2;
            }
            else if (constr[i] == 3)
            {
                info[*ind1 + i] = 11.;
                constr[i] = 3;
            }
        }
        else if (info[*ind1 + i] == 21.)
        {
            if (constr[i] == 1)
            {
                info[*ind1 + i] = 1.;
                constr[i] = 1;
            }
            else if (constr[i] == 2)
            {
                info[*ind1 + i] = 20.;
                constr[i] = 2;
            }
            else if (constr[i] == 3)
            {
                info[*ind1 + i] = 21.;
                constr[i] = 3;
            }
        }
        else if (info[*ind1 + i] == 12.)
        {
            if (constr[i] == 1)
            {
                info[*ind1 + i] = 2.;
                constr[i] = 1;
            }
            else if (constr[i] == 2)
            {
                info[*ind1 + i] = 10.;
                constr[i] = 2;
            }
            else if (constr[i] == 3)
            {
                info[*ind1 + i] = 12.;
                constr[i] = 3;
            }
        }
        else if (info[*ind1 + i] == 22.)
        {
            if (constr[i] == 1)
            {
                info[*ind1 + i] = 2.;
                constr[i] = 1;
            }
            else if (constr[i] == 2)
            {
                info[*ind1 + i] = 20.;
                constr[i] = 2;
            }
            else if (constr[i] == 3)
            {
                info[*ind1 + i] = 22.;
                constr[i] = 3;
            }
        }
        /* L10: */
    }
    return 0;
}    /* dmsk1_ */

/*  --------------------------------------------------------------------- */
int
dmsk2_ (
	MYREAL *info,
	long *comm, long *part, long *ind1, long *np,
	MYREAL *d0, MYREAL *dnp, MYREAL *eps,
	long *errc, long *diagn)
{
  //    extern int dprj0_ (), dprj1_ ();
    MYREAL a, b, p1, p2;

    /*  This routine controls if the separable boundary conditions d(0)=d0 */
    /*  and d(np)=dnp are compatible with the first and the last domain of */
    /*  constraints. The error flag is reset correspondingly. */
    /*  Details on the parameters INFO, IND1 and COMM, PART, NP, D0, DNP, */
    /*  EPS, ERRC, DIAGN can be found in subroutines */
    /*  DSTINF and DBVSSC respectively. */
    /* Parameter adjustments */
    --info;

    /* Function Body */
    dprj0_ (&c__1, &info[1], comm, part, np, &p1, &p2);
    if (!(p1 - *eps <= *d0 && *d0 <= p2 + *eps))
    {
        info[*ind1] = 0.;
        *errc = 9;
        diagn[0] = 1;
    }
    dprj0_ (np, &info[1], comm, part, np, &a, &b);
    dprj1_ (&a, &b, &info[4], &info[5], np, &info[1], comm, part, np, &p1, &p2);
    if (!(p1 - *eps <= *dnp && *dnp <= p2 + *eps))
    {
        info[*ind1 + *np - 1] = 0.;
        *errc = 9;
        diagn[*np - 1] = 1;
    }
    return 0;
}    /* dmsk2_ */

/*  --------------------------------------------------------------------- */
int
dprj0_ (
	long *i,
	MYREAL *info,
	long *comm, long *part, long *np,
	MYREAL *p1, MYREAL *p2)
{
    /* Initialized data */

    MYREAL fl0 = 0.;

    long kind;
    MYREAL k, n, del;

    /*  Given the long i , DPRJ0 computes the set B(i) performing the */
    /*  projection of D(i) (subset of the (i-1)i-plane) onto the (i-1)-axis.
    */

    /*  The input parameters COMM,PART,NP are described in DBVSSC; the input
    */
    /*  parameter INFO is described in subroutine DSTINF. */

    /*  OUTPUT PARAMETERS: */

    /*  P1  : floating variable, containing the left extreme of the */
    /*        resulting interval. */

    /*  P2  : floating variable, containing the right extreme of the */
    /*        resulting interval. */
    /* Parameter adjustments */
    --info;

    /* Function Body */
    n = info[1];
    k = info[2];
    kind = (long) info[*comm + 1 + (*i - 1)];
    del = info[*comm + *np + 1 + (*i - 1)];
    /*  No constraint */
    if (kind == 0)
    {
        *p1 = info[4];
        *p2 = info[5];
        /*  Increase constraints */
    }
    else if (kind == 1)
    {
        *p1 = fl0;
        *p2 = n * del / k;
        /*  Decrease constraints */
    }
    else if (kind == 2)
    {
        *p1 = n * del / k;
        *p2 = fl0;
        /*  Convexity constraints */
    }
    else if (kind == 10)
    {
        *p1 = info[4];
        *p2 = del;
        /*  Concavity constraints */
    }
    else if (kind == 20)
    {
        *p1 = del;
        *p2 = info[5];
        /*  Increase and convexity */
    }
    else if (kind == 11)
    {
        *p1 = fl0;
        *p2 = del;
        /*  Increase and concavity */
    }
    else if (kind == 21)
    {
        *p1 = del;
        *p2 = n * del / k;
        /*  Decrease and convexity */
    }
    else if (kind == 12)
    {
        *p1 = n * del / k;
        *p2 = del;
        /*  Decrease and concavity */
    }
    else if (kind == 22)
    {
        *p1 = del;
        *p2 = fl0;
    }
    return 0;
}    /* dprj0_ */

/*  --------------------------------------------------------------------- */
int
dprj1_ (
	MYREAL *a, MYREAL *b, MYREAL *c, MYREAL *d,
	long *i,
	MYREAL *info,
	long *comm, long *part, long *np,
	MYREAL *p1, MYREAL *p2)
{
    /* Initialized data */

    MYREAL fl0 = 0.;

    /* System generated locals */
    MYREAL d__1, d__2;

    /* Local variables */
    long kind;
    MYREAL k, n, del;

    /*  Given the set S=[a,b]x[c,d] and the long i , DPRJ1 performs the */
    /*  intersection of S with the domain D(i) and the projection of the */
    /*  resulting set (a subset of (i-1)i-plane) onto the i-axis . */

    /*  The input parameters COMM,PART,NP are described in DBVSSC; the input
    */
    /*  parameter INFO is described in DSTINF. */

    /*  OUTPUT PARAMETERS: */

    /*  P1  : floating variable, containing the left extreme of the */
    /*        resulting interval. */

    /*  P2  : floating variable, containing the right extreme of the */
    /*        resulting interval. */
    /* Parameter adjustments */
    --info;

    /* Function Body */
    n = info[1];
    k = info[2];
    kind = (long) info[*comm + 1 + (*i - 1)];
    del = info[*comm + *np + 1 + (*i - 1)];
    /*  No constraint */
    if (kind == 0)
    {
        *p1 = *c;
        *p2 = *d;
        /*  Increase constraints */
    }
    else if (kind == 1)
    {
        *p1 = fl0;
        *p2 = -(*a) + n * del / k;
        /*  Decrease constraints */
    }
    else if (kind == 2)
    {
        *p1 = -(*b) + n * del / k;
        *p2 = fl0;
        /*  Convexity constraints */
    }
    else if (kind == 10)
    {
        *p1 = -k * *b / (n - k) + n * del / (n - k);
        *p2 = -(n - k) * *a / k + n * del / k;
        /*  Concavity constraints */
    }
    else if (kind == 20)
    {
        *p1 = -(n - k) * *b / k + n * del / k;
        *p2 = -k * *a / (n - k) + n * del / (n - k);
        /*  Increase and convexity */
    }
    else if (kind == 11)
    {
        *p1 = -k * *b / (n - k) + n * del / (n - k);
        *p2 = -(n - k) * *a / k + n * del / k;
        /*  Increase and concavity */
    }
    else if (kind == 21)
    {
        /* Computing MAX */
        d__1 = fl0, d__2 = -(n - k) * *b / k + n * del / k;
        *p1 = MAX (d__1, d__2);
        *p2 = -k * *a / (n - k) + n * del / (n - k);
        /*  Decrease and convexity */
    }
    else if (kind == 12)
    {
        *p1 = -k * *b / (n - k) + n * del / (n - k);
        /* Computing MIN */
        d__1 = fl0, d__2 = -(n - k) * *a / k + n * del / k;
        *p2 = MIN (d__1, d__2);
        /*  Decrease and concavity */
    }
    else if (kind == 22)
    {
        *p1 = -(n - k) * *b / k + n * del / k;
        *p2 = -k * *a / (n - k) + n * del / (n - k);
    }
    *p1 = MAX (*p1, *c);
    *p2 = MIN (*p2, *d);
    return 0;
}    /* dprj1_ */

/*  --------------------------------------------------------------------- */
int
dprj2_ (
	MYREAL *a, MYREAL *b, MYREAL *c, MYREAL *d,
	long *i,
	MYREAL *info,
	long *comm, long *part, long *np,
	MYREAL *p1, MYREAL *p2)
{
    /* Initialized data */

    MYREAL fl0 = 0.;

    /* System generated locals */
    MYREAL d__1, d__2;

    /* Local variables */
    long kind;
    MYREAL k, n, del;

    /*  Given the set s=[a,b]x[c,d] and the long i, DPRJ2 performs the */
    /*  intersection of S with the domain D(i) and the projection of the */
    /*  resulting set (subset of (i-1)i-plane) onto the (i-1)-axis . */

    /*  The input parameters COMM,PART,NP are described in DBVSSC; the input
    */
    /*  parameter INFO is described in DSTINF. */

    /*  OUTPUT PARAMETERS: */

    /*  P1  : floating variable, containing the left extreme of the */
    /*        resulting interval. */

    /*  P2  : floating variable, containing the right extreme of the */
    /*        resulting interval. */
    /* Parameter adjustments */
    --info;

    /* Function Body */
    n = info[1];
    k = info[2];
    kind = (long) info[*comm + 1 + (*i - 1)];
    del = info[*comm + *np + 1 + (*i - 1)];
    /*  No constraints */
    if (kind == 0)
    {
        *p1 = *a;
        *p2 = *b;
        /*  Increase constraints */
    }
    else if (kind == 1)
    {
        *p1 = fl0;
        *p2 = -(*c) + n * del / k;
        /*  Decrease constraints */
    }
    else if (kind == 2)
    {
        *p1 = -(*d) + n * del / k;
        *p2 = fl0;
        /*  Convexity constraints */
    }
    else if (kind == 10)
    {
        *p1 = -(n - k) * *d / k + n * del / k;
        *p2 = -k * *c / (n - k) + n * del / (n - k);
        /*  Concavity constraints */
    }
    else if (kind == 20)
    {
        *p1 = -k * *d / (n - k) + n * del / (n - k);
        *p2 = -(n - k) * *c / k + n * del / k;
        /*  Increase and convexity */
    }
    else if (kind == 11)
    {
        /* Computing MAX */
        d__1 = fl0, d__2 = -(n - k) * *d / k + n * del / k;
        *p1 = MAX (d__1, d__2);
        *p2 = -k * *c / (n - k) + n * del / (n - k);
        /*  Increase and concavity */
    }
    else if (kind == 21)
    {
        *p1 = -k * *d / (n - k) + n * del / (n - k);
        *p2 = -(n - k) * *c / k + n * del / k;
        /*  Decrease and convexity */
    }
    else if (kind == 12)
    {
        *p1 = -(n - k) * *d / k + n * del / k;
        *p2 = -k * *c / (n - k) + n * del / (n - k);
        /*  Decrease and concavity */
    }
    else if (kind == 22)
    {
        *p1 = -k * *d / (n - k) + n * del / (n - k);
        /* Computing MIN */
        d__1 = fl0, d__2 = -(n - k) * *c / k + n * del / k;
        *p2 = MIN (d__1, d__2);
    }
    *p1 = MAX (*p1, *a);
    *p2 = MIN (*p2, *b);
    return 0;
}    /* dprj2_ */

/*  --------------------------------------------------------------------- */
int
dscdrc_ (
	 long *n,
	 MYREAL *x, MYREAL *y, MYREAL *d,
	 long *opt, long *np,
	 MYREAL *eps, MYREAL *d20, MYREAL *d2np,
	 MYREAL (*rho) (MYREAL *), MYREAL (*rhoi) (MYREAL *),
	 MYREAL *a1, MYREAL *a2, MYREAL *h, MYREAL *d2,
	 long *errc)
{
    /* Initialized data */

    MYREAL fl0 = 0.;
    MYREAL fl1 = 1.;
    MYREAL fl2 = 2.;
    MYREAL fl4 = 4.;

    /* System generated locals */
    long i__1;
    MYREAL d__1, d__2;

    /* Local variables */
    MYREAL diff2, a, b, c, e, f, g;
    long i, q;
    MYREAL gama, alpha, p1, p2, q1, q2, dd, hh;
    //  extern  int dintrs_ ();
    //    extern MYREAL dsl_ ();

    /*  DSCDRC computes the sequence  d2(i) , i=0,1,...,np , of second */
    /*  derivatives at the knots. The basic idea is that the vanishing second
    */
    /*  derivatives (which are admissible by virtue of the theory involved in
    */
    /*  the routines called previously) can be locally changed to non-zero */
    /*  values without modifying the monotonicity and/or convexity. */
    /*  Let us consider the restriction to the i-th subinterval of the */
    /*  Bernstein-Bezier net for the C(2) spline with zero derivatives given
    */
    /*  by subroutine DALG3. Let A, B and C be the second, third and */
    /*  (int(n/2))-th point of the net, and let E, F, and G be given by a */
    /*  symmetric construction. */

    /*            B_______________C___G_______________F */
    /*           /           .             .           \ */
    /*          /      .                         .      \ */
    /*         /  .D                                 H.  \ */
    /*        A                                           E */
    /*       /                                             \ */
    /*      /                                               \ */
    /*     /                                                 \ */

    /*  Then the 'intermediate net' obtained inserting the straight lines */
    /*  trough A-C and E-F is shape-preserving and we can take as the 'final
    */
    /*  net' the union of two convex combination of A-B-C , A-D-C and H-F-G ,
    */
    /*  E-H-D respectively. Expressing the net in term of the second */
    /*  derivatives, the points D, B and H,F lead to restriction like */
    /*  d2(i).IN.[a1(1,i),a1(2,i)] , d2(i+1).IN.[a2(1,i),a2(2,i)] */
    /*  This construction must be repeated for all the subintervals and so */
    /*  d2(i) .IN. [a2(1,i-1),a2(2,i-1)].INT.[a1(1,i),a1(2,i)] . */

    /*  The input parameters N,X,Y,D,OPT,NP,EPS,D20,D2NP,RHO,RHOI and the */
    /*  input ones D2,ERRC are documented in subroutine DBVSSC. */
    /* Parameter adjustments */
    --a2;
    --a1;

    /* Function Body */
    i__1 = *np - 1;
    for (i = 0; i <= i__1; ++i)
    {
        h[i] = x[i + 1] - x[i];
        /* L10: */
    }
    i__1 = *np - 1;
    for (i = 0; i <= i__1; ++i)
    {
        /*  Compute the points of the 'original' and 'intermediate' net. */
        a = y[i] + h[i] * d[i] / *n;
        b = y[i] + fl2 * h[i] * d[i] / *n;
        e = y[i + 1] - h[i] * d[i + 1] / *n;
        f = y[i + 1] - fl2 * h[i] * d[i + 1] / *n;
        alpha = (*n - 4) / 2 / (*n - fl4);
        c = alpha * f + (fl1 - alpha) * b;
        g = (fl1 - alpha) * f + alpha * b;
        gama = fl1 / ((*n - 4) / 2 + fl1);
        dd = gama * c + (fl1 - gama) * a;
        hh = gama * g + (fl1 - gama) * e;
        /*  Define the left and the right restriction for the second finite */
        /*  difference of the net. */
        /* Computing MIN */
        d__1 = fl0, d__2 = dd - fl2 * a + y[i];
        a1[(i << 1) + 1] = MIN (d__1, d__2);
        /* Computing MAX */
        d__1 = fl0, d__2 = dd - fl2 * a + y[i];
        a1[(i << 1) + 2] = MAX (d__1, d__2);
        /* Computing MIN */
        d__1 = fl0, d__2 = y[i + 1] - fl2 * e + hh;
        a2[((i + 1) << 1) + 1] = MIN (d__1, d__2);
        /* Computing MAX */
        d__1 = fl0, d__2 = y[i + 1] - fl2 * e + hh;
        a2[((i + 1) << 1) + 2] = MAX (d__1, d__2);
        /* L20: */
    }
    /*  Take the intersection of the left and right restrictions for the */
    /*  same second differences and translate it in terms of the second */
    /*  derivatives. */
    /* Computing 2nd power */
    d__1 = h[0];
    a1[1] = a1[1] * *n * (*n - 1) / (d__1 * d__1);
    /* Computing 2nd power */
    d__1 = h[0];
    a1[2] = a1[2] * *n * (*n - 1) / (d__1 * d__1);
    i__1 = *np - 1;
    for (i = 1; i <= i__1; ++i)
    {
        dintrs_ (&a1[(i << 1) + 1], &a1[(i << 1) + 2], &a2[(i << 1) + 1],
                 &a2[(i << 1) + 2], &p1, &p2);
        /* Computing 2nd power */
        d__1 = h[i];
        a1[(i << 1) + 1] = p1 * *n * (*n - 1) / (d__1 * d__1);
        /* Computing 2nd power */
        d__1 = h[i];
        a1[(i << 1) + 2] = p2 * *n * (*n - 1) / (d__1 * d__1);
        /* L30: */
    }
    /* Computing 2nd power */
    d__1 = h[*np - 1];
    a1[(*np << 1) + 1] = a2[(*np << 1) + 1] * *n * (*n - 1) / (d__1 * d__1);
    /* Computing 2nd power */
    d__1 = h[*np - 1];
    a1[(*np << 1) + 2] = a2[(*np << 1) + 2] * *n * (*n - 1) / (d__1 * d__1);
    /*  The internal derivatives are defined as the admissible value closest
    */
    /*  to the central second divided difference of the data. */
    i__1 = *np - 1;
    for (i = 1; i <= i__1; ++i)
    {
        diff2 =
            ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]) / (h[i] +
                    h[i -
                      1]);
        d2[i] = dsl_ (&a1[(i << 1) + 1], &a1[(i << 1) + 2], &diff2);
        /* L40: */
    }
    q = *opt % 100 / 10;
    if (q == 1)
    {
        /*  No boundary condition is required. Take the first and last */
        /*  derivative as the middle of admissible values. */
        d2[0] = (a1[1] + a1[2]) / fl2;
        d2[*np] = (a1[(*np << 1) + 1] + a1[(*np << 1) + 2]) / fl2;
    }
    else if (q == 2)
    {
        /*  Non-separable boundary conditions are required. Check if these can
         be */
        /*  satisfied by admissible derivatives. */
        /* Computing MIN */
        d__1 = (*rho) (&a1[1]), d__2 = (*rho) (&a1[2]);
        q1 = MIN (d__1, d__2);
        /* Computing MAX */
        d__1 = (*rho) (&a1[1]), d__2 = (*rho) (&a1[2]);
        q2 = MAX (d__1, d__2);
        dintrs_ (&a1[(*np << 1) + 1], &a1[(*np << 1) + 2], &q1, &q2, &p1, &p2);
        if (p1 > p2 + *eps)
        {
            /*  The boundary conditions cannot be satisfied. Set the error fla
            g and */
            /*  define the first and the last derivative as the nearest point
            to the */
            /*  admissible and the boundary interval. */
            if (*errc == 0)
            {
                *errc = 15;
            }
            else if (*errc == 6)
            {
                *errc = 16;
            }
            else if (*errc == 9)
            {
                *errc = 17;
            }
            else if (*errc == 10)
            {
                *errc = 18;
            }
            d__1 = (p1 + p2) / fl2;
            d2[*np] = dsl_ (&a1[(*np << 1) + 1], &a1[(*np << 1) + 2], &d__1);
            d__1 = (*rhoi) (&d2[*np]);
            d2[0] = dsl_ (&a1[1], &a1[2], &d__1);
        }
        else
        {
            /*  It is possible to satisfy the boundary conditions. */
            d2[*np] = (p1 + p2) / fl2;
            d2[0] = (*rhoi) (&d2[*np]);
        }
    }
    else if (q == 3)
    {
        /*  Separable boundary conditions are required. Check if they are */
        /*  compatible with the admissible intervals and, if not, set the */
        /*  error flag and take the admissible points nearest to the boundary
        */
        /*  values. Otherwise take simply the boundary values. */
        if (*d20 < a1[1] - *eps || *d20 > a1[2] + *eps
                || *d2np < a1[(*np << 1) + 1] - *eps
                || *d2np > a1[(*np << 1) + 2] + *eps)
        {
            if (*errc == 0)
            {
                *errc = 11;
            }
            else if (*errc == 6)
            {
                *errc = 12;
            }
            else if (*errc == 9)
            {
                *errc = 13;
            }
            else if (*errc == 10)
            {
                *errc = 14;
            }
        }
        d2[0] = dsl_ (&a1[1], &a1[2], d20);
        d2[*np] = dsl_ (&a1[(*np << 1) + 1], &a1[(*np << 1) + 2], d2np);
    }
    return 0;
}    /* dscdrc_ */

/*  --------------------------------------------------------------------- */
MYREAL
dsl_ (MYREAL *a, MYREAL *b, MYREAL *c)
{
    /* System generated locals */
    MYREAL ret_val;

    /*  Given the interval [a,b] and the number c, dsl is c if c belongs */
    /*  to [a,b], otherwise, it is the nearest extreme to c. */
    if (*c <= *a)
    {
        ret_val = *a;
    }
    else if (*c >= *b)
    {
        ret_val = *b;
    }
    else
    {
        ret_val = *c;
    }
    return ret_val;
}    /* dsl_ */

/*  --------------------------------------------------------------------- */
int
dsqtab_ (
	 MYREAL *x, MYREAL *y,
	 long *np,
	 MYREAL *xtab,
	 long *ntab, long *y0opt, long *y1opt, long *y2opt, long *n, long *k,
	 MYREAL *d, MYREAL *d2, MYREAL *tb, MYREAL *l, MYREAL *laux0, MYREAL *laux1, MYREAL *laux2, MYREAL *y0tab, MYREAL *y1tab, MYREAL *y2tab)
{
    /* System generated locals */
    long i__1, i__2;

    /* Local variables */
    long i, j;
    //  extern  int dlspis_ ();
    //  extern MYREAL dbl_ ();
    long ind;
    //  extern MYREAL dbl1_ (), dbl2_ ();
    long ind1;

    /*  DSQTAB evaluates the spline and/or its first derivative and/or its */
    /*  second derivative at the points  xtab(j) , j=0,1,...,ntab  using */
    /*  a sequential search for finding the interval  [x(i),x(i+1)] in which
    */
    /*  the tabulation point falls. The input (X,Y,NP,XTAB,NTAB,Y0OPT, */
    /*  Y1OPT,Y2OPT,N,K,D,D2,TB) and the output (Y0TAB,Y1TAB,Y2TAB) */
    /*  parameters have been explained in subroutine DBVSSE. For the others */
    /*  see subroutines DTRMB, DLSPIS. */
    /* Parameter adjustments */
    --tb;

    /* Function Body */
    ind = 0;
    ind1 = 1;
    i__1 = *ntab;
    for (j = 0; j <= i__1; ++j)
    {
        /*  Compute the index  ind  such that  x(ind).LE.xtab(j).LT.x(ind+1) .
         */
        if (x[0] <= xtab[j])
        {
            i__2 = *np - 1;
            for (i = ind1; i <= i__2; ++i)
            {
                if (x[i] <= xtab[j])
                {
                    ind = i;
                }
                /* L20: */
            }
        }
        /*  Check if  ind  selects a new subinterval. */
        if (ind != ind1 || j == 0)
        {
            /*  Call subprogram  DLSPIS  to compute the linear shape-preservin
            g */
            /*  interpolating spline  l:=l(x)  at */
            /*      x(ind)+p*(x(ind+1)-x(ind))/n , p=0,1,...,n . */
            dlspis_ (x, y, d, d2, np, n, k, &ind, l);
            if (*y0opt == 1)
            {
                /*  Evaluate the spline at  xtab(j)  using new values of  l .
                */
                y0tab[j] =
                    dbl_ (&xtab[j], n, l, &x[ind], &x[ind + 1],
                          &tb[*n * (*n + 1) / 2], &c__0, laux0);
            }
            if (*y1opt == 1)
            {
                /*  Evaluate the first derivative of the spline at  xtab(j)  u
                sing new */
                /*  values of  l . */

                y1tab[j] =
                    dbl1_ (&xtab[j], n, l, &x[ind], &x[ind + 1],
                           &tb[(*n - 1) * *n / 2], &c__0, laux1);
            }
            if (*y2opt == 1)
            {
                /*  Evaluate the second derivative of the spline at  xtab(j)
                using new */
                /*  values of  l . */
                y2tab[j] =
                    dbl2_ (&xtab[j], n, l, &x[ind], &x[ind + 1],
                           &tb[(*n - 2) * (*n - 1) / 2], &c__0, laux2);
            }
        }
        else
        {
            if (*y0opt == 1)
            {
                /*  Evaluate the spline at  xtab(j)  using old values of  l .
                */
                y0tab[j] =
                    dbl_ (&xtab[j], n, l, &x[ind], &x[ind + 1],
                          &tb[*n * (*n + 1) / 2], &c__1, laux0);
            }
            if (*y1opt == 1)
            {
                /*  Evaluate the first derivative of the spline at  xtab(j)  u
                sing old */
                /*  values of  l . */
                y1tab[j] =
                    dbl1_ (&xtab[j], n, l, &x[ind], &x[ind + 1],
                           &tb[(*n - 1) * *n / 2], &c__1, laux1);
            }
            if (*y2opt == 1)
            {
                /*  Evaluate the second derivative of the spline at  xtab(j)
                using old */
                /*  values of  l . */
                y2tab[j] =
                    dbl2_ (&xtab[j], n, l, &x[ind], &x[ind + 1],
                           &tb[(*n - 2) * (*n - 1) / 2], &c__1, laux2);
            }
        }
        ind1 = ind;
        ind = 0;
        /* L30: */
    }
    return 0;
}    /* dsqtab_ */

/*  --------------------------------------------------------------------- */
int
dstinf_ (
	 long *opt,
	 MYREAL *d0, MYREAL *dnp,
	 long *constr, long *n, long *k,
	 MYREAL *x, MYREAL *y, MYREAL *d,
	 long *np, long *comm, long *part,
	 MYREAL *eps,
	 MYREAL (*beta) (MYREAL *), MYREAL (*betai) (MYREAL *),
	 MYREAL *daux2, MYREAL *daux3, MYREAL *info,
	 long *errc, long *diagn)
{
    /* Initialized data */

    MYREAL fl0 = 0.;
    MYREAL flem4 = 1e-4;
    MYREAL fl2 = 2.;

    /* System generated locals */
    long i__1;
    MYREAL d__1, d__2, d__3;

    /* Local variables */
    MYREAL d2im1;
    //  extern  int dtdc_ (), dmsk1_ (), dmsk2_ ();
    MYREAL iaux0;
    long i, p, q, r;
    MYREAL d2i, iauxnp;
    long ind1, ind2, ind3;

    /*  DSTINF computes the information needed in the other parts of the */
    /*  program using the data-dependent input parameters and stores it in */
    /*  the output array INFO. */

    /*  The parameters OPT,N,K,X,Y,D,NP,COMM,PART,EPS,BETA,BETAI,ERRC,DIAGN */
    /*  are described in subroutine DBVSSC . */

    /*  Items of possible interest are: */

    /*  INFO  : floating array, of bounds 1:COMM+PART*NP+NP+1. It is composed
    */
    /*          of four parts: the first, of bounds 1:comm, contains the */
    /*          global information n, k , the maximum of the first divided */
    /*          differences of initial points and the lower and upper bounds
    */
    /*          for the derivatives, bounds which are used when no */
    /*          constraints are imposed (see the parameter OPT described in */
    /*          DBVSSC) or when the constraints must be relaxed; the second,
    */
    /*          of bounds  comm+1:comm+np, contains information about */
    /*          constraints in the interval (x(i),x(i+1)) , i=0,1,...,np-1 ;
    */
    /*          if: */
    /*          info((comm+1)+i)= 0 - no attribute; */
    /*          info((comm+1)+i)= 1 - the data are increasing; */
    /*          info((comm+1)+i)= 2 - the data are decreasing; */
    /*          info((comm+1)+i)=10 - the data are convex; */
    /*          info((comm+1)+i)=11 - the data are increasing and convex; */
    /*          info((comm+1)+i)=12 - the data are decreasing and convex; */
    /*          info((comm+1)+i)=20 - the data are concave; */
    /*          info((comm+1)+i)=21 - the data are increasing and concave; */
    /*          info((comm+1)+i)=22 - the data are decreasing and concave. */
    /*          The third part, of bounds comm+np+1:comm+part*np, contains */
    /*          the first divided differences of initial points */
    /*              ( y(i+1)-y(i) ) / ( x(i+1)-x(i) ) ,  i=0,1,...,np-1 . */
    /*          The fourth, of bounds comm+part*np+1:comm+part*np+np+1, */
    /*          contains, eventually, the initial estimates of the first */
    /*          derivatives which are then used to compute the constrained */
    /*          best approximation (see the description of the input */
    /*          parameter OPT  and of the array D in subr. DBVSSC). More */
    /*          precisely, having defined  p := opt/100 , if p=2 it contains
    */
    /*          the Bessel estimates, if p=3 it contains a set of third order
    */
    /*          accurate estimates giving a co-monotone cubic Hermite */
    /*          interpolant (see subr. DTDC described later), if p=4 it */
    /*          contains a set of values given by the user; if p=1 this part
    */
    /*          of INFO is not referenced. */
    /* Parameter adjustments */
    --daux2;
    --info;

    /* Function Body */
    r = *opt % 10;
    q = *opt % 100 / 10;
    p = *opt / 100;
    ind1 = *comm + 1;
    ind2 = *comm + *np + 1;
    ind3 = *comm + (*np << 1) + 1;
    /*  Set the first and the second components of INFO to n and k */
    /*  respectively. */
    info[1] = (MYREAL) (*n);
    info[2] = (MYREAL) (*k);
    /*  Compute the first divided differences of the initial points and set */
    /*  info(3) to their maximum. */
    info[3] = fl0;
    i__1 = *np - 1;
    for (i = 0; i <= i__1; ++i)
    {
        info[ind2 + i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        /* Computing MAX */
        d__2 = info[3], d__3 = (d__1 = info[ind2 + i], fabs (d__1));
        info[3] = MAX (d__2, d__3);
        /* L10: */
    }
    /*  Compute the lower and upper bounds for derivatives */
    info[4] = -(*n) * info[3] / *k;
    info[5] = *n * info[3] / *k;
    /*  If eps.LE.0 it is reset to flem4. */
    if (*eps <= fl0)
    {
        *eps = flem4;
    }
    /*  Compute the relative tollerance of the method. */
    *eps *= info[3];
    /*  Set the second part of INFO. Firstly, all the components are */
    /*  initialized with 0. */
    i__1 = *np - 1;
    for (i = 0; i <= i__1; ++i)
    {
        info[ind1 + i] = 0.;
        /* L20: */
    }
    /*  Monotonicity is required: check if the initial points are increasing
    */
    /*  or decreasing in each interval ( x(i), x(i+1) ) , i=0,1,...,np-1 . */
    if (r == 1 || r == 3 || r == 4)
    {
        i__1 = *np - 1;
        for (i = 0; i <= i__1; ++i)
        {
            if (info[ind2 + i] >= fl0)
            {
                ++info[ind1 + i];
            }
            else
            {
                info[ind1 + i] += 2;
            }
            /* L30: */
        }
    }
    /*  Convexity is required: check if the initial points are concave or */
    /*  convex in each interval ( x(i), x(i+1) ) , i=1,...,np-2 . */
    if (r == 2 || r == 3 || r == 4)
    {
        i__1 = *np - 2;
        for (i = 1; i <= i__1; ++i)
        {
            d2im1 = info[ind2 + i] - info[ind2 + (i - 1)];
            d2i = info[ind2 + (i + 1)] - info[ind2 + i];
            if ((d2im1 >= *eps && d2i >= -(*eps))
                    || (d2im1 >= -(*eps) && d2i >= *eps) || (fabs (d2im1) <= *eps
                            && fabs (d2i) <= *eps))
            {
                info[ind1 + i] += 10;
            }
            else if ((d2im1 <= -(*eps) && d2i <= *eps)
                     || (d2im1 <= *eps && d2i <= -(*eps)))
            {
                info[ind1 + i] += 20;
            }
            /* L40: */
        }
        /*  The convexity in the first and in the last interval is defined as
        the */
        /*  second and the second to last, respectively. */
        info[ind1] += (long) info[ind1 + 1] / 10 * 10;
        info[ind1 + (*np - 1)] += (long) info[ind1 + (*np - 2)] / 10 * 10;
    }
    /*  In the case  r=4 , that is when the constraint selection */
    /*  is made on any interval, we compare the kind given by the data with */
    /*  those given by the array CONSTR */
    if (r == 4)
    {
        dmsk1_ (&info[1], constr, comm, part, &ind1, np);
    }
    /*  In the case q=3, the kind in the first and last subinterval */
    /*  is compared with the boundary conditions */
    if (q == 3)
    {
        dmsk2_ (&info[1], comm, part, &ind1, np, d0, dnp, eps, errc, diagn);
    }
    /*  If p=2 the Bessel derivatives are stored in the fourth */
    /*  part of INFO. */
    if (p == 2)
    {
        i__1 = *np - 1;
        for (i = 1; i <= i__1; ++i)
        {
            info[ind3 + i] =
                ((x[i + 1] - x[i]) * info[ind2 + (i - 1)] +
                 (x[i] - x[i - 1]) * info[ind2 + i]) / (x[i + 1] - x[i - 1]);
            /* L50: */
        }
        if (q == 1)
        {
            /*  If no boundary condition is imposed, set the first and last */
            /*  derivatives using the standard formulas for Bessel interpolati
            on. */
            info[ind3] =
                ((x[1] - x[0]) * (fl2 * info[ind2] - info[ind2 + 1]) +
                 (x[2] - x[1]) * info[ind2]) / (x[2] - x[0]);
            info[ind3 + *np] =
                ((x[*np] - x[*np - 1]) * (fl2 * info[ind2 + (*np - 1)] -
                                          info[ind2 + (*np - 2)]) + (x[*np - 1] -
                                                                     x[*np -
                                                                       2]) *
                 info[ind2 + (*np - 1)]) / (x[*np] - x[*np - 2]);
        }
        else
        {
            /*  Compute the first and last derivatives taking into account bot
            h the */
            /*  slopes of the data and the restriction imposed by the boundary
             */
            /*  conditions */
            iaux0 =
                ((x[1] - x[0]) * (fl2 * info[ind2] - info[ind2 + 1]) +
                 (x[2] - x[1]) * info[ind2]) / (x[2] - x[0]);
            iauxnp =
                ((x[*np] - x[*np - 1]) * (fl2 * info[ind2 + (*np - 1)] -
                                          info[ind2 + (*np - 2)]) + (x[*np - 1] -
                                                                     x[*np -
                                                                       2]) *
                 info[ind2 + (*np - 1)]) / (x[*np] - x[*np - 2]);
            info[ind3] =
                ((x[*np] - x[*np - 1]) * iaux0 +
                 (x[1] - x[0]) * (*betai) (&iauxnp)) / (x[1] - x[0] + (x[*np] -
                                                        x[*np -
                                                          1]));
            info[ind3 + *np] =
                ((x[*np] - x[*np - 1]) * (*beta) (&iaux0) +
                 (x[1] - x[0]) * iauxnp) / (x[1] - x[0] + (x[*np] - x[*np - 1]));
        }
        /*  If p=3 then the set of third order accurate estimates, computed by
         */
        /*  subr. DTDC, is stored in the fourth part of INFO. */
    }
    else if (p == 3)
    {
        dtdc_ (np, x, comm, part, &info[1], &daux2[1], daux3);
        /*  If p=4 then the set of values given by the user is stored in */
        /*  the fourth part of INFO. */
    }
    else if (p == 4)
    {
        i__1 = *np;
        for (i = 0; i <= i__1; ++i)
        {
            info[ind3 + i] = d[i];
            /* L60: */
        }
    }
    return 0;
}    /* dstinf_ */

/*  --------------------------------------------------------------------- */
int
dtdc_ (
       long *np,
       MYREAL *x,
       long *comm, long *part,
       MYREAL *info, MYREAL *dd2, MYREAL *dd3)
{
    /* Initialized data */

    MYREAL fl0 = 0.;
    MYREAL fl2 = 2.;
    MYREAL fl3 = 3.;

    /* System generated locals */
    long i__1;
    MYREAL d__1, d__2;

    /* Local variables */
    MYREAL tmin, tmax, f;
    long i;
    MYREAL d1, e1, e2, e3, d2, p1, q1, q2, p2, q12, q32, si, ti;
    //  extern MYREAL dmdian_ (), dmnmod_ ();
    long ind;
    MYREAL tit;
    long ind1;

    /*  Given the initial points ( x(i), y(i) ) , i=0,1,...,np , DTDC */
    /*  computes a sequence  ds(0),ds(1),...,ds(np)  of estimates of the */
    /*  function's derivatives which are third order accurate and furnish a */
    /*  cubic Hermite interpolant preserving monotonicity. */
    /*  The method is composed by the two following fundamental steps: */
    /*  1 - compute an initial estimate of ds(i), i=0,1,...,np , which is */
    /*      third or higher order accurate; */
    /*  2 - refine it imposing the monotonicity constraint. */
    /*  The computation of ds(i) needs the points x(i+j), j = -3,-2,...,3 , */
    /*  consequently, the boundary values ds(0), ds(1), ds(2) and ds(np-2), */
    /*  ds(np-1), ds(np) are computed in an approximate way. Although they */
    /*  are still third order accurate, may not preserve the monotonicity. */
    /*  For more details see \3\ . */

    /*  The input parameter NP,X,COMM,PART are described in subr. DBVSSC; the
    */
    /*  parameter INFO is described in subr. DSTINF . */

    /*  The computed values are stored in the last part of INFO. */
    /* Parameter adjustments */
    --dd2;
    --info;

    /* Function Body */
    ind = *comm + *np + 1;
    ind1 = *comm + *part * *np + 1;
    /*  Compute the second divided differences of the initial points. */
    i__1 = *np - 1;
    for (i = 1; i <= i__1; ++i)
    {
        dd2[i] = (info[ind + i] - info[ind + i - 1]) / (x[i + 1] - x[i - 1]);
        /* L10: */
    }
    /*  Compute the third divided differences of the initial points */
    i__1 = *np - 2;
    for (i = 1; i <= i__1; ++i)
    {
        dd3[i] = (dd2[i + 1] - dd2[i]) / (x[i + 2] - x[i - 1]);
        /* L20: */
    }
    /*  Compute approximate values for  f[x(-1),x(0),x(1),x(2)]  and */
    /*  f[x(np-2),x(np-1),x(np),x(np+1)] ; they are needed for the */
    /*  computation of ds(2) and ds(np-2). */
    dd3[0] =
        dd3[1] + (dd3[2] - dd3[1]) / (x[4] - x[0]) * (x[0] + x[1] - x[2] - x[3]);
    dd3[*np - 1] =
        dd3[*np - 2] + (dd3[*np - 3] - dd3[*np - 2]) / (x[*np - 4] -
                x[*np]) * (x[*np] +
                           x[*np - 1] -
                           x[*np - 2] -
                           x[*np - 3]);
    i__1 = *np - 2;
    for (i = 2; i <= i__1; ++i)
    {
        /*  ds(i) : initialization */
        e1 = dmnmod_ (&dd3[i - 2], &dd3[i - 1]);
        e2 = dmnmod_ (&dd3[i - 1], &dd3[i]);
        e3 = dmnmod_ (&dd3[i], &dd3[i + 1]);
        d__1 = dd2[i - 1] + e1 * (x[i] - x[i - 2]);
        d__2 = dd2[i] + e2 * (x[i] - x[i + 1]);
        q1 = info[ind + i - 1] + (x[i] - x[i - 1]) * dmnmod_ (&d__1, &d__2);
        d__1 = dd2[i] + e2 * (x[i] - x[i - 1]);
        d__2 = dd2[i + 1] + e3 * (x[i] - x[i + 2]);
        q2 = info[ind + i] - (x[i + 1] - x[i]) * dmnmod_ (&d__1, &d__2);
        f = (q1 + q2) / fl2;
        /*  Refinement */
        tit = dmnmod_ (&q1, &q2);
        d1 = dmnmod_ (&dd2[i - 1], &dd2[i]);
        d2 = dmnmod_ (&dd2[i], &dd2[i + 1]);
        p1 = info[ind + i - 1] + d1 * (x[i] - x[i - 1]);
        p2 = info[ind + i] + d2 * (x[i] - x[i + 1]);
        ti = dmnmod_ (&p1, &p2);
        si = dmnmod_ (&info[ind + i - 1], &info[ind + i]);
        /* Computing MIN */
        d__1 = fl0, d__2 = fl3 * si, d__1 = MIN (d__1, d__2), d__2 =
                                                fl3 * ti / fl2, d__1 = MIN (d__1, d__2);
        tmin = MIN (d__1, tit);
        /* Computing MAX */
        d__1 = fl0, d__2 = fl3 * si, d__1 = MAX (d__1, d__2), d__2 =
                                                fl3 * ti / fl2, d__1 = MAX (d__1, d__2);
        tmax = MAX (d__1, tit);
        d__1 = tmin - f;
        d__2 = tmax - f;
        info[ind1 + i] = f + dmnmod_ (&d__1, &d__2);
        /* L50: */
    }
    /*  ds(1): initialization */
    q12 =
        info[ind] + dd2[1] * (x[1] - x[0]) + dd3[0] * (x[1] - x[0]) * (x[1] -
                x[2]);
    q32 =
        info[ind] + dd2[1] * (x[1] - x[0]) + dd3[1] * (x[1] - x[0]) * (x[1] -
                x[2]);
    e1 = dmnmod_ (dd3, &dd3[1]);
    e2 = dmnmod_ (&dd3[1], &dd3[2]);
    q1 = dmdian_ (&info[ind], &q12, &q32);
    d__1 = dd2[1] + e1 * (x[1] - x[0]);
    d__2 = dd2[2] + e2 * (x[1] - x[3]);
    q2 = info[ind + 1] - (x[2] - x[1]) * dmnmod_ (&d__1, &d__2);
    f = (q1 + q2) / fl2;
    /*  refinement */
    tit = dmnmod_ (&q1, &q2);
    d2 = dmnmod_ (&dd2[1], &dd2[2]);
    p1 = info[ind] + dd2[1] * (x[1] - x[0]);
    p2 = info[ind + 1] + d2 * (x[1] - x[2]);
    ti = dmnmod_ (&p1, &p2);
    si = dmnmod_ (&info[ind], &info[ind + 1]);
    /* Computing MIN */
    d__1 = fl0, d__2 = fl3 * si, d__1 = MIN (d__1, d__2), d__2 =
                                            fl3 * ti / fl2, d__1 = MIN (d__1, d__2);
    tmin = MIN (d__1, tit);
    /* Computing MAX */
    d__1 = fl0, d__2 = fl3 * si, d__1 = MAX (d__1, d__2), d__2 =
                                            fl3 * ti / fl2, d__1 = MAX (d__1, d__2);
    tmax = MAX (d__1, tit);
    d__1 = tmin - f;
    d__2 = tmax - f;
    info[ind1 + 1] = f + dmnmod_ (&d__1, &d__2);
    /*  ds(np-1): initialization */
    e1 = dmnmod_ (&dd3[*np - 3], &dd3[*np - 2]);
    e2 = dmnmod_ (&dd3[*np - 2], &dd3[*np - 1]);
    d__1 = dd2[*np - 2] + e1 * (x[*np - 1] - x[*np - 3]);
    d__2 = dd2[*np - 1] + e2 * (x[*np - 1] - x[*np]);
    q1 =
        info[ind + *np - 2] + (x[*np - 1] - x[*np - 2]) * dmnmod_ (&d__1, &d__2);
    q12 =
        info[ind + *np - 2] + dd2[*np - 1] * (x[*np - 1] - x[*np - 2]) + dd3[*np -
                1] *
        (x[*np - 1] - x[*np - 2]) * (x[*np - 1] - x[*np]);
    q32 =
        info[ind + *np - 2] + dd2[*np - 1] * (x[*np - 1] - x[*np - 2]) + dd3[*np -
                2] *
        (x[*np - 1] - x[*np - 2]) * (x[*np - 1] - x[*np]);
    q2 = dmdian_ (&info[ind + *np - 1], &q12, &q32);
    f = (q1 + q2) / fl2;
    /*  Refinement */
    tit = dmnmod_ (&q1, &q2);
    d1 = dmnmod_ (&dd2[*np - 2], &dd2[*np - 1]);
    p1 = info[ind + *np - 2] + d1 * (x[*np - 1] - x[*np - 2]);
    p2 = info[ind + *np - 1] + dd2[*np - 1] * (x[*np - 1] - x[*np]);
    ti = dmnmod_ (&p1, &p2);
    si = dmnmod_ (&info[ind + *np - 2], &info[ind + *np - 1]);
    /* Computing MIN */
    d__1 = fl0, d__2 = fl3 * si, d__1 = MIN (d__1, d__2), d__2 =
                                            fl3 * ti / fl2, d__1 = MIN (d__1, d__2);
    tmin = MIN (d__1, tit);
    /* Computing MAX */
    d__1 = fl0, d__2 = fl3 * si, d__1 = MAX (d__1, d__2), d__2 =
                                            fl3 * ti / fl2, d__1 = MAX (d__1, d__2);
    tmax = MAX (d__1, tit);
    d__1 = tmin - f;
    d__2 = tmax - f;
    info[ind1 + *np - 1] = f + dmnmod_ (&d__1, &d__2);
    /*  ds(0): */
    q1 =
        info[ind] + dd2[1] * (x[0] - x[1]) + dd3[0] * (x[0] - x[1]) * (x[0] -
                x[2]);
    q2 =
        info[ind] + dd2[1] * (x[0] - x[1]) + dd3[1] * (x[0] - x[1]) * (x[0] -
                x[2]);
    info[ind1] = dmdian_ (&info[ind], &q1, &q2);
    /*  ds(np): */
    q1 =
        info[ind + *np - 1] + dd2[*np - 1] * (x[*np] - x[*np - 1]) + dd3[*np -
                1] *
        (x[*np] - x[*np - 2]) * (x[*np] - x[*np - 1]);
    q2 =
        info[ind + *np - 1] + dd2[*np - 1] * (x[*np] - x[*np - 1]) + dd3[*np -
                2] *
        (x[*np] - x[*np - 2]) * (x[*np] - x[*np - 1]);
    info[ind1 + *np] = dmdian_ (&info[ind + *np - 1], &q1, &q2);
    return 0;
}    /* dtdc_ */

/*  --------------------------------------------------------------------- */
int
dtrmb_ (
	long *n,
	MYREAL *tb)
{
    /* Initialized data */

    MYREAL fl1 = 1.;

    /* System generated locals */
    long i__1, i__2;

    /* Local variables */
    long i, k, ind, ind1;

    /*  DTRMB    computes the binomial terms */
    /*      i!/(k!*(i-k)!) , i=1,2,...,n , k=0,1,...,i . */

    /*  INPUT PARAMETERS */

    /*  N     : long variable, containing the largest binomial term */
    /*          needed. */

    /*  OUTPUT PARAMETERS */

    /*  TB    : floating array, of bounds  1:N*(N+1)/2+N  , containing */
    /*          the values   i!/(k!*(i-k)!)  , k=0,1,...,i , in the */
    /*          elements   TB(i*(i+1)/2),...,TB((i*(i+1)/2)+i) . */
    /* Parameter adjustments */
    --tb;

    /* Function Body */
    tb[1] = fl1;
    tb[2] = fl1;
    ind = 1;
    i__1 = *n;
    for (i = 2; i <= i__1; ++i)
    {
        ind1 = ind;
        ind = i * (i + 1) / 2;
        tb[ind] = fl1;
        tb[ind + i] = fl1;
        i__2 = i - 1;
        for (k = 1; k <= i__2; ++k)
        {
            tb[ind + k] = tb[ind1 + k] + tb[ind1 + k - 1];
            /* L10: */
        }
        /* L20: */
    }
    return 0;
}    /* dtrmb_ */

/*  --------------------------------------------------------------------- */
boolean
dtst_ (MYREAL *a, MYREAL *b, MYREAL *c, MYREAL *d, MYREAL *eps)
{
    /* System generated locals */
    MYREAL d__1, d__2;
    boolean ret_val;

    /*  DTST checks if two intervals [a,b] and [c,d] differ less than eps. */
    /*  DTST assumes  a.LE.b  and  c.LE.d . */
    ret_val = (d__1 = *a - *c, fabs (d__1)) <= *eps
              && (d__2 = *b - *d, fabs (d__2)) <= *eps;
    return ret_val;
}    /* dtst_ */
