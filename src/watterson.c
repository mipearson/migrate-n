/*
 *  watterson.c
 *  migrate-n
 *
 *  Created by Peter Beerli on June 26 2006.
 *  Copyright 2006 Peter Beerli. All rights reserved.
 *
 */

#include "watterson.h"



double watterson_a(long n)
{
  long k;
  double sum=0.;
  for(k=1; k<n; k++)
    {
      sum += 1. / k;
    }
  return sum;
}
double watterson_b(long n)
{
  long k;
  double sum=0.;
  for(k=1; k<n; k++)
    {
      sum += 1. / (k*k);
    }
  return sum;
}

MYREAL watterson(long segreg, long n)
{
  return (MYREAL) (segreg / watterson_a(n));
}

//wattersonvar[s_, n_] :=
//((an = (Sum[1./k, {k, 1,
// n}]))(thetaw = watterson[s, n]) +  (bn = (Sum[
// 1./(k^2), {k, 1, n}]))thetaw^2 )/(an^2)
MYREAL wattersonvar(double  thetaw, long n)
{
  double an = watterson_a(n);
  double bn = watterson_b(n);
  return (MYREAL) ((an * thetaw + bn * thetaw * thetaw) / (an*an));
}

MYREAL wattersonstd(double  thetaw, long n)
{
  return (MYREAL)(sqrt(wattersonvar(thetaw,n)));
}

