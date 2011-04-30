/*     ALGORITHM 614 COLLECTED ALGORITHMS FROM ACM. */
/*     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.10, NO. 2, */
/*     JUN., 1984, P. 152-160. */
/* Subroutine */
/*! \file inthp.c */

#include "definitions.h"

#include <stdlib.h>
#include <math.h>


#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif



#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif /* MIN */
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif /* MAX */
// allready defined in definitions.h #define FOUR 4
#define ZERO 0.
// already defined in definitions.h  #define ONE 1


int inthp (MYREAL *a, MYREAL *b, MYREAL *d__, MYREAL (*f) (MYREAL, void *),
           long *m, MYREAL *p, MYREAL *eps, long *inf, void *helper,
           MYREAL *quadr);



int
inthp (MYREAL *a, MYREAL *b, MYREAL *d__, MYREAL (*f) (MYREAL, void *),
       long *m, MYREAL *p, MYREAL *eps, long *inf, void *helper,
       MYREAL *quadr)
{
    /* System generated locals */
    MYREAL r__1;
    MYREAL zero = ZERO;

    /* Local variables */
    MYREAL alfa = 0, exph = 0, exph0 = 0, c__ = 0, h__ = 0;
    long i__ = 0, k = 0, l = 0, n = 0;
    MYREAL s = 0, t = 0, u = 0, v = 0, w = 0, c0 = 0, e1 = 0, h0 = 0, h1 = 0;
    long i1 = 0, l1 = 0, m1 = 0, m2 = 0, n1 = 0;
    MYREAL s1 = 0, v0 = 0, v1 = 0, v2 = 0, w1 = 0, w2 = 0, w3 = 0, w4 = 0, ba =
                                            0, pi = 0, sr = 0, sq2 = 0, cor = 0, sum = 0;
    boolean inf1 = FALSE, inf2 = FALSE;
    MYREAL eps3 = 0, sum1 = 0, sum2 = 0;
    /*        THIS SUBROUTINE COMPUTES INTEGRAL OF FUNCTIONS WHICH */
    /*     MAY HAVE SINGULARITIES AT ONE OR BOTH END-POINTS OF AN */
    /*     INTERVAL (A,B), SEE [1, 2]. IT CONTAINS FOUR DIFFERENT */
    /*     QUADRATURE ROUTINES: ONE OVER A FINITE INTERVAL (A,B), */
    /*     TWO OVER (A,+INFINITY), AND ONE OVER (-INFINITY,+INFINITY). */
    /*     OF THE TWO FORMULAS OVER (A,+INFINITY), THE FIRST (INF=2 */
    /*     BELOW) IS MORE SUITED TO NON-OSCILLATORY INTEGRANDS, WHILE */
    /*     THE SECOND (INF=3) IS MORE SUITED TO OSCILLATORY INTEGRANDS. */
    /*        THE USER SUPPLIES THE INTEGRAND FUNCTION, HE SPECIFIES THE */
    /*     INTERVAL, AS WELL AS THE RELATIVE ERROR TO WHICH THE INTE- */
    /*     GRAL IS TO BE EVALUATED. */
    /*        THE FORMULAS ARE OPTIMAL IN CERTAIN HARDY SPACES H(P,DD), */
    /*     SEE [1, 2]. HERE DD IS AN OPEN DOMAIN IN THE COMPLEX PLANE, */
    /*     A AND B BELONG TO THE BOUNDARY OF DD AND H(P,DD), P.GT.1, IS */
    /*     THE SET OF ALL ANALYTIC FUNCTONS IN DD WHOSE P-TH NORM DEFI- */
    /*     NED AS IN [2] IS FINITE. */
    /*        IF THE USER IS UNABLE TO SPECIFY THE PARAMETERS P AND D */
    /*     OF THE SPACE H(P,DD) TO WHICH HIS INTEGRAND BELONGS, THE */
    /*     ALGORITHM TERMINATES ACCORDING TO A HEURISTIC CRITERION, SEE */
    /*     [2] AND COMMENTS TO EPS. */
    /*        IF THE USER CAN SPECIFY THE PARAMETERS P AND D OF THE */
    /*     SPACE H(P,DD) TO WHICH HIS INTEGRAND BELONGS, THE ALGORITHM */
    /*     TERMINATES WITH AN ANSWER HAVING A GUARANTEED ACCURACY ( DE- */
    /*     TEMINISTIC CRITERION, SEE [1, 2] AND COMMENTS TO EPS). */
    /*     INPUT PARAMETERS */
    /*     A = LOWER LIMIT OF INTEGRATION (SEE COMMENTS TO INF). */
    /*     B = UPPER LIMIT OF INTEGRATION (SEE COMMENTS TO INF). */
    /*     D = A PARAMETER OF THE CLASS H(P,DD) (SEE COMMENTS TO */
    /*         INF). */
    /*         USER SETS D: */
    /*         HEURISTIC TERMINATION */
    /*       = ANY MYREAL NUMBER */
    /*         DETERMINISTIC TERMINATION */
    /*       = A NUMBER IN THE RANGE 0.LT.D.LE.PI/2. */
    /*     F = A NAME OF AN EXTERNAL INTEGRAND FUNCTION TO BE */
    /*         SUPPLIED BY THE USER. F(X) COMPUTES THE VALUE OF */
    /*         A FUNCTION F AT A POINT X. THE STATEMENT */
    /*         ...EXTERNAL F... MUST APPEAR IN THE MAIN PROGRAM. */
    /*     M = MAXIMAL NUMBER OF FUNCTION EVALUATIONS ALLOWED IN */
    /*         THE COMPUTATIONS, M.GE.3.( ALTERED ON EXIT ). */
    /*     P = 0, 1, .GT.1  A PARAMETER OF THE CLASS H(P,DD). */
    /*         USER SETS P: */
    /*       = 0 - HEURISTIC TERMINATION. */
    /*       = 1 - DETERMINISTIC TERMINATION WITH THE INFINITY */
    /*             NORM. */
    /*      .GT.1 -DETERMINISTIC TERMINATION WITH THE P-TH NORM. */
    /*   EPS = A MYREAL NUMBER - THE RELATIVE ERROR BOUND - SEE */
    /*         REMARKS BELOW. ( ALTERED ON EXIT ). */
    /*   INF = 1, 2, 3, 4 - INFORMATION PARAMETER. ( ALTERED ON EXIT ). */
    /*       = 1 SIGNIFIES AN INFINITE INTERVAL (A,B)=MYREAL LINE, */
    /*           A AND B ANY NUMBERS. */
    /*           DETERMINISTIC TERMINATION - */
    /*           DD=STRIP(Z:ABS(IM(Z)).LT.D). */
    /*       = 2 SIGNIFIES A SEMI-INFINITE INTERVAL (A, +INFINITY) */
    /*           USER SUPPLIES A, B ANY NUMBER. */
    /*           QUADRATURE SUITED TO NON-OSCILLATORY INTEGRANDS. */
    /*           DETERMINISTIC TERMINATION - */
    /*           DD=SECTOR(Z:ABS(ARG(Z-A)).LT.D). */
    /*       = 3 SIGNIFIES A SEMI INFINITE INTERVAL (A,+INFINITY) */
    /*           USER SUPPLIES A, B ANY NUMBER. */
    /*           QUADRATURE SUITED TO OSCILLATORY INTEGRANDS. */
    /*           DETERMINISTIC TERMINATION - */
    /*           DD=REGION(Z:ABS(ARG(SINH(Z-A))).LT.D). */
    /*       = 4 SIGNIFIES A FINITE INTERVAL (A,B). */
    /*           USER SUPPLIES A AND B. */
    /*           DETERMINISTIC TERMINATION - */
    /*           DD=LENS REGION(Z:ABS(ARG((Z-A)/(B-Z))).LT.D). */




    /*     OUTPUT PARAMETERS */
    /*     M = THE NUMBER OF FUNCTION EVALUATIONS USED IN THE */
    /*         QUADRATURE. */
    /*   EPS = THE RELATIVE ERROR BOUND (SEE REMARKS BELOW). */
    /*         DETERMINISTIC TERMINATION */
    /*       = THE RELATIVE ERROR REXA BOUND, I.E., */
    /*                 REXA(F,M(OUTPUT)) .LE. EPS. */
    /*         HEURISTIC TERMINATION */
    /*       = MAX(EPS(INPUT),MACHEP). */
    /*   INF = 0, 1 - DETERMINISTIC TERMINATION */
    /*       = 0 COMPUTED QUADRATURE QCOM(F,M(EPS)), SEE REMARKS */
    /*           BELOW. */
    /*       = 1 COMPUTED QUADRATURE QCOM(F,M1), SEE REMARKS */
    /*           BELOW. */
    /*   INF = 2, 3, 4 - HEURISTIC TERMINATION. */
    /*       = 2 INTEGRATION COMPLETED WITH EPS=MAX(EPS(INPUT), */
    /*           MACHEP). WE CAN EXPECT THE RELATIVE ERROR */
    /*           REXA TO BE OF THE ORDER OF EPS (FOR SOME P.GE.1). */
    /*       = 3 INTEGRATION NOT COMPLETED. ATTEMPT TO EXCEED THE */
    /*           MAXIMAL ALLOWED NUMBER OF FUNCTION EVALUATIONS M. */
    /*           TRUNCATION CONDITIONS (SEE [2]) SATISFIED. QUADR */
    /*           SET TO BE EQUAL TO THE LAST TRAPEZOIDAL APPRO- */
    /*           XIMATION. IT IS LIKELY THAT QUADR APPROXIMATES THE */
    /*           INTEGRAL IF M IS LARGE. */
    /*       = 4 INTEGRATION NOT COMPLETED. ATTEMPT TO EXCEED THE */
    /*           MAXIMAL ALLOWED NUMBER OF FUNCTION EVALUATIONS M. */
    /*           TRUNCATION CONDITIONS (SEE [2]) NOT SATISFIED. */
    /*           QUADR SET TO BE EQUAL TO THE COMPUTED TRAPEZOIDAL */
    /*           APPROXIMATION. IT IS UNLIKELY THAT QUADR APPROXIMATES */
    /*           THE INTEGRAL. */
    /*   INF = 10, 11, 12, 13 - INCORRECT INPUT */
    /*       = 10  M.LT.3. */
    /*       = 11  P DOES NOT SATISFY P=0, P=1 OR P.GT.1 OR IN THE */
    /*             CASE OF DETERMINISTIC TERMINATION D DOES NOT */
    /*             SATISFY 0.LT.D.LE.PI/2. */
    /*       = 12  A.GE.B IN CASE OF A FINITE INTERVAL. */
    /*       = 13  INF NOT EQUAL TO 1, 2, 3, OR 4. */
    /*   QUADR = THE COMPUTED VALUE OF QUADRATURE. */




    /*     REMARKS: */
    /*         LET  QEXA(F,M)  ( QCOM(F,M) ) BE THE EXACT (COMPUTED) */
    /*         VALUE OF THE QUADRATURE WITH M FUNCTION EVALUATIONS, */
    /*         AND LET  REXA(F,M) ( RCOM(F,M) ) BE THE RELATIVE ERROR */
    /*         OF QEXA (QCOM) ,I.E., */
    /*            REXA(F,M)=ABS(INTEGRAL(F)-QEXA(F,M))/NORM(F), */
    /*            RCOM(F,M)=ABS(INTEGRAL(F)-QCOM(F,M))/NORM(F), */
    /*         WITH THE NOTATION 0/0=0. */
    /*             DUE TO THE ROUNDOFF ONE CANNOT EXPECT THE ERROR */
    /*         RCOM TO BE LESS THAN THE RELATIVE MACHINE PRECISION */
    /*         MACHEP. THEREFORE THE INPUT VALUE OF EPS IS CHANGED */
    /*         ACCORDING TO THE FORMULA */
    /*                   EPS=MAX(EPS,MACHEP). */
    /*         DETERMINISTIC TERMINATON CASE */
    /*             THE NUMBER OF FUNCTON EVALUATIONS M(EPS) IS COMPUTED */
    /*         SO THAT THE ERROR REXA IS NO GREATER THAN EPS,I.E., */
    /*         (*)     REXA(F,M(EPS)) .LE. EPS . */
    /*         IF M(EPS).LE.M THEN THE QUADRATURE QCOM(F,M(EPS)) IS COM- */
    /*         PUTED. OTHERWISE, WHICH MEANS THAT EPS IS TOO SMALL WITH */
    /*         RESPECT TO M, THE QUADRATURE QCOM(F,M1) IS COMPUTED, WHERE */
    /*         M1=2*INT((M-1)/2)+1. IN THIS CASE EPS IS CHANGED TO THE */
    /*         SMALLEST NUMBER FOR WHICH THE ESTIMATE (*) HOLDS WITH */
    /*         M(EPS)=M1 FUNCTION EVALUATIONS. */
    /*         HEURISTIC TERMINATION CASE */
    /*             WE CAN EXPECT THE RELATIVE ERROR REXA TO BE OF THE */
    /*         ORDER OF EPS, SEE [2]. IF EPS IS TOO SMALL WITH RESPECT */
    /*         TO M THEN THE QUADRATURE QCOM(F,M) IS COMPUTED. */
    /*         ROUNDOFF ERRORS */
    /*             IN BOTH DETERMINISTIC AND HEURISTIC CASES THE ROUND- */
    /*         OFF ERROR */
    /*                    ROFF=ABS(QEXA(F,M)-QCOM(F,M)) */
    /*         CAN BE ESTIMATED BY */
    /*         (**)       ROFF .LE. 3*C1*R*MACHEP, */
    /*         WHERE  R=QCOM(ABS(F),M)+(1+2*C2)/3*SUM(W(I),I=1,2,...M) */
    /*         AND C1 IS OF THE ORDER OF UNITY, C1=1/(1-3*MACHEP), W(I) */
    /*         ARE THE WEIGHTS OF THE QUADRATURE, SEE [2], AND C2 IS */
    /*         A CONSTANT ESTIMATING THE ACCURACY OF COMPUTING FUNCTION */
    /*         VALUES, I.E., */
    /*               ABS(EXACT(F(X))-COMPUTED(F(X))).LE.C2*MACHEP. */
    /*         IF THE INTEGRAND VALUES ARE COMPUTED INACCURATELY, I.E., */
    /*         C2 IS LARGE, THEN THE ESTIMATE (**) IS LARGE AND ONE CAN */
    /*         EXPECT THE ACTUAL ERROR ROFF TO BE LARGE. NUMERICAL TESTS */
    /*         INDICATE THAT THIS HAPPENS ESPECIALLY WHEN THE INTEGRAND */
    /*         IS EVALUATED INACCURATELY NEAR A SINGULARITY. THE WAYS OF */
    /*         CIRCUMVENTING SUCH PITFALLS ARE EXPLAINED IN [2]. */
    /*     REFERENCES: */
    /*     [1] SIKORSKI,K., OPTIMAL QUADRATURE ALGORITHMS IN HP */
    /*            SPACES, NUM. MATH., 39, 405-410 (1982). */
    /*     [2] SIKORSKI,K., STENGER,F., OPTIMAL QUADRATURES IN */
    /*            HP SPACES, ACM TOMS. */

    pi = atan (1.) * 4.;

    /*     CHECK THE INPUT DATA */

    if (*inf < 1 || *inf > 4)
    {
        *inf = 13;
        return 0;
    }
    if (*m < 3)
    {
        *inf = 10;
        return 0;
    }
    if (*p < 1. && *p != 0.)
    {
        *inf = 11;
        return 0;
    }
    if (*p >= 1. && (*d__ <= 0. || *d__ > pi / 2.))
    {
        *inf = 11;
        return 0;
    }
    if (*inf == 4 && *a >= *b)
    {
        *inf = 12;
        return 0;
    }

    sq2 = sqrt (2.);
    i1 = *inf - 2;
    ba = *b - *a;
    n1 = 0;

    /*     COMPUTE THE RELATIVE MACHINE PRECISION AND CHECK */
    /*     THE VALUE OF EPS.  CAUTION...THIS LOOP MAY NOT WORK ON A */
    /*     MACHINE THAT HAS AN ACCURATED ARITHMETIC PROCESS COMPARED */
    /*     TO THE STORAGE PRECISION.  THE VALUE OF U MAY NEED TO BE */
    /*     SIMPLY DEFINED AS THE RELATIVE ACCURACY OF STORAGE PRECISION. */

    u = 1.;
    do
    {
        u /= 10.;
        t = u + 1.;
    }
    while (1. != t);
    u *= 10.;
    if (*eps < u)
    {
        *eps = u;
    }

    if (*p == 0.)
    {
        /*     SET UP DATA FOR THE HEURISTIC TERMINATION */
        h__ = 1.;
        h0 = 1.;
        eps3 = *eps / 3.;
        sr = sqrt (*eps);
        v1 = *eps * 10.;
        v2 = v1;
        m1 = *m - 1;
        n = (long) ((MYREAL) (m1 / 2));
        m2 = n;
        l1 = 0;
        inf1 = TRUE;
        inf2 = FALSE;
    }
    else
    {
        /*     SET UP DATA FOR THE DETERMINISTIC TERMINATION */
        if (*p == 1.)
        {
            alfa = 1.;
        }
        if (*p > 1.)
        {
            alfa = (*p - 1.) / *p;
        }
        c__ =
            pi * 2. / (1. - 1. / EXP (pi * sqrt (alfa))) + pow ((MYREAL) FOUR,
                    alfa) / alfa;
        w = log (c__ / *eps);
        w1 = 1. / (pi * pi * alfa) * w * w;
        n = (long) w1;
        if (w1 > (MYREAL) n)
        {
            ++n;
        }
        if (w1 == 0.)
        {
            n = 1;
        }
        n1 = (n << 1) + 1;
        sr = sqrt (alfa * (MYREAL) n);
        if (n1 <= *m)
        {
            *m = n1;
            n1 = 0;
        }
        else
        {
            /*     EPS TOO SMALL WITH RESPECT TO M. COMPUTE THE NEW EPS */
            /*     GUARANTEED BY THE VALUE OF M. */
            n1 = 1;
            n = (long) ((MYREAL) ((*m - 1) / 2));
            sr = sqrt (alfa * (MYREAL) n);
            *m = (n << 1) + 1;
            *eps = c__ / EXP (pi * sr);
        }
        h__ = *d__ * 2. / sr;
        sum2 = 0.;
        l1 = n;
        k = n;
        inf1 = FALSE;
        inf2 = FALSE;
        h0 = h__;
    }

    /*     INITIALIZE THE QUADRATURE */
    i__ = 0;
    if (*inf == 1)
    {
        sum = (*f) (zero, helper);
    }
    if (*inf == 2)
    {
        r__1 = *a + 1.;
        sum = (*f) (r__1, helper);
    }
    if (*inf == 3)
    {
        r__1 = *a + log (sq2 + 1.);
        sum = (*f) (r__1, helper) / sq2;
    }
    if (*inf == 4)
    {
        r__1 = (*a + *b) / 2.;
        sum = (*f) (r__1, helper) / 4. * ba;
    }

    /*     COMPUTE WEIGHTS, NODES AND FUNCTION VALUES */

L60:
    exph = EXP (h__);
    exph0 = EXP (h0);
    h1 = h0;
    e1 = exph0;
    u = 0.;
    cor = 0.;

L70:
    if (i1 < 0)
    {
        goto L80;
    }
    else if (i1 == 0)
    {
        goto L90;
    }
    else
    {
        goto L100;
    }

L80:
    v = (*f) (h1, helper);
    h1 += h__;
    goto L150;

L90:
    r__1 = *a + e1;
    v = e1 * (*f) (r__1, helper);
    e1 *= exph;
    goto L150;

L100:
    if (*inf == 4)
    {
        goto L140;
    }
    w1 = sqrt (e1 + 1. / e1);
    w2 = sqrt (e1);
    if (e1 < .1)
    {
        goto L110;
    }
    s = log (e1 + w1 * w2);
    goto L130;
L110:
    w3 = e1;
    w4 = e1 * e1;
    c0 = 1.;
    s = e1;
    s1 = e1;
    t = 0.;
L120:
    c0 = -(MYREAL) c0 *(t + .5) * (t * 2. + 1.) / (t * 2. + 3.) / (t + 1.);
    t += 1.;
    w3 *= w4;
    s += c0 * w3;
    if (s == s1)
    {
        goto L130;
    }
    s1 = s;
    goto L120;
L130:
    r__1 = *a + s;
    v = w2 / w1 * (*f) (r__1, helper);
    e1 *= exph;
    goto L150;

L140:
    w1 = e1 + 1.;
    r__1 = (*a + *b * e1) / w1;
    v = e1 / w1 / w1 * (*f) (r__1, helper) * ba;
    e1 *= exph;

    /*     SUMMATION ALGORITHM */

L150:
    ++i__;
    sum1 = u + v;
    if (fabs (u) < fabs (v))
    {
        goto L160;
    }
    cor = v - (sum1 - u) + cor;
    goto L170;
L160:
    cor = u - (sum1 - v) + cor;
L170:
    u = sum1;
    if (i__ < l1)
    {
        goto L70;
    }

    /*     SWITCH TO CHECK TRUNCATION CONDITION ( HEURISTIC */
    /*     TERMINATION) */

    if (inf1)
    {
        goto L190;
    }

    /*     SWITCH TO COMPUTE THE MIDORDINATE APPROXIMATION */
    /*     ( HEURISTIC TERMINATION ) OR TO STOP ( DETERMINIS- */
    /*     TIC TERMINATION) */

    if (inf2)
    {
        goto L210;
    }

    /*     SET UP PARAMETERS TO CONTINUE SUMMATION */

    l1 = k;
L180:
    inf2 = TRUE;
    i__ = 0;
    exph = 1. / exph;
    h0 = -(MYREAL) h0;
    e1 = 1. / exph0;
    h1 = h0;
    h__ = -(MYREAL) h__;
    goto L70;

    /*     TRUNCATION CONDITION */

L190:
    v0 = v1;
    v1 = v2;
    v2 = fabs (v);
    if (v0 + v1 + v2 <= eps3)
    {
        goto L200;
    }
    if (i__ < m2)
    {
        goto L70;
    }
    n1 = 5;
L200:
    if (inf2)
    {
        k = i__;
    }
    if (!inf2)
    {
        l = i__;
    }
    v1 = *eps * 10.;
    v2 = v1;
    m2 = m1 - l;
    if (!inf2)
    {
        goto L180;
    }

    /*     N1=5 - TRUNCATION CONDITION NOT SATISFIED */

    if (n1 == 5)
    {
        goto L260;
    }

    /*     TRUNCATION CONDITION SATISFIED, SUM2=TRAPEZOIDAL */
    /*     APPROXIMATION */

    sum2 = sum1 + cor + sum;
    m2 = (k + l) << 1;

    /*     CHECK THE NUMBER OF FUNCTION EVALUATIONS */

    if (m2 > m1)
    {
        goto L240;
    }

    /*     INITIALIZE ITERATION */

    inf1 = FALSE;
    inf2 = FALSE;
    l1 = l;
    i__ = 0;
    h__ = -(MYREAL) h__;
    h0 = h__ / 2.;
    goto L60;

    /*     P.GE.1 = DETERMINISTIC TERMINATION */

L210:
    if (*p >= 1.)
    {
        goto L220;
    }

    /*     COMPUTE THE MIDORDINATE APPROXIMATION SUM1 */

    h__ = -(MYREAL) h__;
    sum1 = (sum1 + cor) * h__;
    w1 = (sum1 + sum2) / 2.;

    /*     TERMINATION CONDITION */

    if ((r__1 = sum1 - sum2, fabs (r__1)) <= sr)
    {
        goto L230;
    }

    /*     SET UP DATA FOR THE NEXT ITERATION */

    m2 <<= 1;
    if (m2 > m1)
    {
        goto L250;
    }
    i__ = 0;
    k <<= 1;
    l <<= 1;
    l1 = l;
    h__ /= 2.;
    h0 = h__ / 2.;
    sum2 = w1;
    inf2 = FALSE;
    goto L60;

    /*     FINAL RESULTS */

L220:
    *quadr = -(MYREAL) h__ *(sum1 + cor + sum);
    *inf = n1;
    return 0;

L230:
    *quadr = w1;
    *inf = 2;
    *m = m2 + 1;
    return 0;

L240:
    *quadr = sum2;
    *inf = 3;
    *m = k + l + 1;
    return 0;

L250:
    *quadr = w1;
    *inf = 3;
    *m = m2 / 2 + 1;
    return 0;

L260:
    *quadr = u + cor + sum;
    *inf = 4;
    *m = k + l + 1;
    return 0;





}
