#ifndef __INTEGRATE_H
#define __INTEGRATE_H
extern int inthp (MYREAL *a, MYREAL *b, MYREAL *d__,
                      MYREAL (*f) (MYREAL, void *), long *m, MYREAL *p,
                      MYREAL *eps, long *inf, MYREAL *quadr, void *helper);

#endif
