/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 R A N D O M   G E N E R A T O R   R O U T I N E S 
 
 creates options structures,
 reads options from parmfile if present
 
 prints options,
 and finally helps to destroy itself.
                                                                                                               
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2007 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
$Id: random.c 1850 2011-03-21 20:08:31Z beerli $
-------------------------------------------------------*/
/* \file random.c */
#include "sighandler.h"
#include "migration.h"
#include "heating.h"
#include "random.h"
#include "migrate_mpi.h"
#ifdef WIN32
#include <sys/timeb.h>
#else /* WIN32 */
#include <sys/time.h>
#endif /* WIN32 */


#ifdef MERSENNE_TWISTER
#include "SFMT.c"
#endif

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif


#ifdef SNOWLEOPARD
#include <dispatch/dispatch.h>
#endif


#ifdef PTHREADS
extern tpool_t heating_pool;
#endif

/* prototypes ----------------------------------------- */
void getseed (option_fmt * options);
void swap_ptr (long **ptr, long **newptr);

MYREAL randum (void);
#ifdef PTHREADS
MYREAL randum_thread (void);
#endif
long random_integer(long low, long high);

#ifdef QUASIRANDOM
#define MYINDEX 0

#ifdef CUD
int  nextn=1, kk=1,ii=0,jj=0;
void setup_cud(int thiskk)
{
  strncpy(generator,"Quasi-random number generator: Completely Uniform Distributed numbers", 79);
  kk = thiskk;
}

MYINLINE MYREAL get_cud()
{
  MYREAL rr;
  MYREAL value;
  const MYREAL  magicnumbers[]={2.,2.,4.,15.,65.,315.,1586.,8036.,40448.,200401.,972536.,
	       4609846.,21310545.,96017492., 421606654.,1804551131. };
  //double ps[]={2.,3.,5.,7.,11.,13.,17.,19.,23.,29.,31.,37.,41.,47.,53.,59.,61.,67.,71.,73.,79.,83.,89.,97.,101.};
  const MYREAL lps[] =
  {0.69314718055994530942, 1.0986122886681096914, 
   1.6094379124341003746, 1.9459101490553133051, 
   2.3978952727983705441, 2.5649493574615367361, 
   2.8332133440562160802, 2.9444389791664404600, 
   3.1354942159291496908, 3.3672958299864740272, 
   3.4339872044851462459, 3.6109179126442244444, 
   3.7135720667043078039, 3.8501476017100585868, 
   3.9702919135521218341, 4.0775374439057194506, 
   4.1108738641733112488, 4.2046926193909660597, 
   4.2626798770413154213, 4.2904594411483911291, 
   4.3694478524670214942, 4.4188406077965979235, 
   4.4886363697321398383, 4.5747109785033828221, 
   4.6151205168412594509} ; 
  //rr=kk*log(ps[jj]);
#ifdef PTHREADS

    if ((pthread_mutex_lock (&(heating_pool->random_lock))) != 0)
        error ("pthread_mutex_lock failed in get_cud()");
#endif

  rr = kk * lps[jj];
  
  if (jj >= ii)
    { 
      kk++; 
      jj=0;
    }
  else 
    jj++;

  if(kk >= magicnumbers[ii])
    {
      kk = 1; 
      ii++;
    }
  value = rr - floor(rr);
#ifdef PTHREADS

    if ((pthread_mutex_unlock (&(heating_pool->random_lock))) != 0)
        error ("pthread_mutex_unlock failed in get_cud()");
#endif

    return   value;
  //printf("%f\n",rr);
}

MYINLINE MYREAL get_quasi()
{
    return get_cud();
}

#endif

#ifdef KOROBOV

unsigned int x0=137,ii=0; nextn=1;

void setup_korobov()
{
  strncpy(generator,"Quasi-random number generator: Korobov sequence", 79);
}

MYREAL get_korobov()
{  
unsigned long  a =17364, m=65521;
  unsigned long xn;
  double x;

  xn=(a*x0)%m;
   x=(xn+0.0)/(m+0.0);
   x0=xn; nextn++;
   if (nextn>m){ ii++; nextn=2+ii;}
 return x;
}


MYINLINE MYREAL get_quasi()
{
  return get_korobov();
}

#endif

#ifdef WEYL
static int s, nextn=8;

void setup_weyl()
{
  strncpy(generator,"Quasi-random number generator: Weyl sequence", 79);
}


MYINLINE MYREAL get_weyl()
{
    double r=sqrt(7.0),rr;
    int next;

    next=nextn*nextn;
    rr=r*next;

//    rr= (r*nextn-floor(r*nextn))*nextn;
    rr=rr-floor(rr);

    nextn++;
    return rr;
}

MYINLINE MYREAL get_quasi()
{
  return get_weyl();
}
#endif /*WEYL*/

#ifdef HALTON

#define MAX_D 500 
static int s, nextn=8;

double e, quasi[100];
static int prime[MAX_D];
static double iprime[MAX_D];
static int  primroots[][10]={{1, 2, 3, 3, 8, 11,12,14, 7,18},
    {12,13,17,18,29,14,18,43,41,44},
    {40,30,47,65,71,28,40,60,79,89},
    {56,50,52,61,108,56,66,63,60,66},
    {104,76,111,142,71,154,118,84,127,142},
    {84,105,186,178,188,152,165,159,103,205}, 
    {166,173,188,181,91,233,210,217,153,212},
};

static int warnockOpt[]={1,  2,  2,  5,  3,  7,  3,  10,  18, 11, 
    17, 5, 17,  26, 40, 14, 40, 44, 12, 31,
    45, 70,8,   38, 82, 8,  12, 38, 47, 70,
    29, 57, 97, 110,32, 48, 84, 124,155,26,
    69, 83, 157,171, 8, 22, 112,205, 15, 31,
    61, 105,127,212,12, 57, 109,133,179,210,
    231,34, 161,199,222,255,59, 120,218,237,
    278,341,54, 110,176,218,280,369,17, 97, 
    193,221,331,350,419,21, 85, 173,221,243,
    288,424,45, 78, 173,213,288,426,455,138,
}; 

int gohalt(double *,double *, double *);
double get_halton();
int primes();
int power(int, int, int);
int inhalt(int dimen, int atmost, double tiny, double *quasi);
int power(int a, int b, int m)
{ int i,c=1;
    for(i=0;i<b;i++)
        c=(c*a)%m;
    return c;
} 

void setup_halton()
{
  strncpy(generator,"Quasi-random number generator: Halton sequence", 79);
}

double get_halton()
{ double wq[100],dq[100];
    gohalt(quasi,dq,wq);
    return wq[MYINDEX];
} 

int inhalt(int dimen, int atmost, double tiny, double *quasi)
{
    double delta, f;
    int i,m;
    
    // check dimen
    primes();
    
    s=dimen;
    if (s<1||s>1000)
        return(-1);
    
    // compute and check tolerance
    
    e=0.9*(1.0/(atmost*prime[s-1])-10.0*tiny);
    delta=100*tiny*(double)(atmost+1)*log10((double)atmost);
    if (delta>=0.09*(e-10.0*tiny))
        return(-2);
    
    // now compute first vector
    
    m=1;
    for (i=0;i<s;i++)
    {       
        iprime[i]=1.0/iprime[i];
        quasi[i]=iprime[i]; 
        m=i*prime[i];
    }
    
    printf("largest prime=%d, %f \n",prime[s-1], quasi[1]);
    
    nextn=2;
    
    return 0;
}

int gohalt(double *quasi, double *dq, double *wq)
{
    int i, j, k, ytemp[40],xtemp[40], ztemp, ktemp, ltemp, mtemp;
    double r;
    double t,f,g,h;
    
    // generate quasi one compoment at a time using radix prime[k] for 
    // component k
    
    
    for (i=0;i<s;i++)
    {
        t=iprime[i];
        f=1.0-quasi[i];
        g=1.0;
        h=t;
        while ((f-h)<e)
            // this checks whether q+h>1-e
        {
            g=h;
            h*=t;
        }
        quasi[i]=g+h-f;
    }
    
    for(i=0;i<s;i++)
    {	  
        k=0; mtemp=nextn; 
        ltemp=prime[i]; 
        
        while(mtemp!=0){
            ytemp[k]=mtemp%ltemp;
            mtemp=mtemp/ltemp;
            k++; 
        }
        
        //generating Optimal primitive root 
        for(j=0;j<k;j++)
        {
            // xtemp[j] = (ytemp[j]*power(primroots[i/10][i%10], nextn%ltemp, ltemp))%ltemp;
            if(j>=1) 
                xtemp[j] =(warnockOpt[i]*power(primroots[i/10][i%10], ytemp[j], 
                                               ltemp)+ytemp[j-1])%ltemp;
            else xtemp[j] =(warnockOpt[i]*power(primroots[i/10][i%10], ytemp[j],
                                                ltemp))%ltemp;
            xtemp[j] -= ytemp[j];
        }  
        
        dq[i]=0;t=iprime[i]; 
        for(j=0;j<k; j++)
        { 
            dq[i] += xtemp[j]*t;
            t *= iprime[i];
        }
        
        dq[i] += quasi[i];
        
        
        // generating Warnock Optimal sequences
        for(j=0;j<k;j++)	   
        {      
            if(j>=1)
                xtemp[j]= (ytemp[j]*power(warnockOpt[i],i+1,ltemp)+ytemp[j-1])%ltemp;
            else
                xtemp[j]= (ytemp[j]*power(warnockOpt[i],i+1,ltemp))%ltemp;
            
            xtemp[j] -= ytemp[j];
        }
        
        wq[i]=0;t=iprime[i];
        for(j=0;j<k; j++)
        {
            wq[i] += xtemp[j]*t;
            t *= iprime[i];
        }
        
        wq[i] += quasi[i];
    }
    
    nextn++;
    return(0);
}

int primes()
{
    int i, j, a[MAX_D+1];
    for (a[1] = 0, i = 2; i <= MAX_D; i++)
        a[i] = 1;
    for (i = 2; i <= MAX_D/2; i++)
        for (j = 2; j <= MAX_D/i; j++)
            a[i*j] = 0;
    for (i = 1, j = 0; i <= MAX_D; i++){
        if (a[i]){
            prime[j] =i; 
			iprime[j]=i;
            j++;
        }
    }   
    return j;
}       

MYINLINE MYREAL get_quasi()
{
  return get_halton();
}

#endif /*HALTON */


#endif

#ifdef MERSENNE_TWISTER
void setup_mersennetwister()
{
  strncpy(generator,"Pseudo-random number generator: Mersenne twister", 79);
}
#endif

/// random integer
//=============================================
MYINLINE long random_integer(long low, long high)
{
    //Math.floor(Math.random()*(N-M+1))%(N-M+1)+M
    long r;
    long tt = high-low+1;
    r =  low + (long) (RANDUM() * tt);
    if(r<low || r>high)
      {
        // warning("Random %li integer is out of bounds (%li, %li)\n",r, low, high);
	if(r>high)
	  return high;
	else
	  return low;
	//	error("Stop because this should never happen");
      }
    return r;
}


///
/// better random_seed function taken from a discussion on 
/// http://sourceware.org/ml/gsl-discuss/2004-q1/msg00071.html
/// March 30, 2007, Robert G. Brown (http://www.phy.duke.edu/~rgb/)
/// at Duke University Dept. of Physics, Box 90305  Durham, N.C. 27708-030
/// suggested the code 
unsigned long int random_seed()
{
  
  unsigned int seed;
  struct timeval tv;
  FILE *devrandom;
  
  if ((devrandom = fopen("/dev/random","r")) == NULL) 
    {
      gettimeofday(&tv,0);
      seed = tv.tv_sec + tv.tv_usec;
    }
  else 
    {
      fread(&seed,sizeof(seed),1,devrandom);
      fclose(devrandom);
    }
  return(seed);
}

void set_seed(long autoseed, long * inseed)
{
  unsigned long timeseed = 0;
   switch (autoseed)
    {
    case AUTO:
      timeseed = random_seed();
#ifdef LCG
      switch (timeseed % 4)
        {
        case 0:
	  ++timeseed; break;
        case 2:
            ++timeseed;
            break;
        case 1:
        case 3:
            break;
        }
#endif
      *inseed = abs((long) timeseed);
      break;
    case NOAUTO:
    case NOAUTOSELF:
        break;
    default:
        error ("Error: Seed value not defined");
        break;
    }
}

void init_lcg(option_fmt *options)
{
  long i;
  for (i = 0; i <= 2; i++)
    seed[i] = 0;
  i = 0;
  do
    {
      seed[i] = options->inseed & 2047;
      printf("RANDOM SEEDING: i=%li seed[i]=%li inseed=%li\n",i,seed[i], options->inseed);
      options->inseed /= 2048;
      i++;
    }
  while (options->inseed != 0);
  if(i>3)
    error("lcg random number generator is broken! random.c:507");
}

void getseed_lcg (option_fmt *options)
{
    strncpy(generator,"Pseudo-random number generator: Least Congruental Generator", 79);
    set_seed(options->autoseed, &options->inseed);
    options->saveseed = options->inseed;
    init_lcg(options);
}

void getseed_mt (option_fmt *options)
{
    strncpy(generator,"Pseudo-random number generator: Mersenne-Twister", 79);
    set_seed(options->autoseed, &options->inseed);
    options->saveseed = options->inseed;
    init_gen_rand((unsigned long)options->inseed);
}

///
/// set up the quasi random number material
/// this function is dependent on running the standard MT or LCG random number generator first.
void getseed_quasi (option_fmt *options)
{
#ifdef CUD
  setup_cud((int) RANDINT (1, MAXLONG-1));
#endif
#ifdef WEYL
  setup_weyl();
#endif
#ifdef KOROBOV
  setup_korobov();
#endif
#ifdef HALTON
  setup_halton();
#endif

}


#ifdef MPI
/// Random number seed function for MPI usage, the main seed gets distributed to 
/// all nodes and then the nodes pick a random number from the master seed for
/// their own seed. given that the nodes typically analyze different loci
/// it will be unlikely that even when picked by the off chance the same seed on two
/// nodes we get the same stream of events. the node ID gets the randomnumber[#ID] from the stream.
void getseed(option_fmt *options)
{
  long test = -1;
  long i;
#ifdef MERSENNE_TWISTER
  getseed_mt(options);
#else
  getseed_lcg(options);
#endif
  MYMPIBCAST (&options->inseed, 1, MPI_LONG, MASTER, comm_world);
#ifdef DEBUG_MPI
  printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
  printf("WARNING all random number seeds on the nodes are identical to the master     WARNING\n");
  printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
  return;
#endif
  if(myID != MASTER)
    {
      for (i = 0; i <= myID; i++)
	{
	  test = RANDINT (1, MAXLONG-1); 
	}
      options->inseed = test;
      //      printf("RANDOM NUMBER SEED: node=%i seed=%li\n",myID,test);
#ifdef MERSENNE_TWISTER
      getseed_mt(options);
#else
      getseed_lcd(options);
#endif
#ifdef QUASIRANDOM
      getseed_quasi(options);
#endif
    }
  else
    {
      printf("MASTER RANDOM NUMBER SEED: node=%i seed=%li\n",myID,options->inseed);
    }
}
#else
void getseed(option_fmt *options)
{
#ifdef MERSENNE_TWISTER
  getseed_mt(options);
#else
  getseed_lcg(options);
#endif
#ifdef QUASIRANDOM
  getseed_quasi(options);
#endif
}
#endif


void
swap_ptr (long **ptr, long **newptr)
{
    long *temp;
    temp = *ptr;
    *ptr = *newptr;
    *newptr = temp;
}

#ifdef PTHREADS
MYREAL
randum_thread (void)
/* thread save random number generator */
{
    MYREAL value;
    if ((pthread_mutex_lock (&(heating_pool->random_lock))) != 0)
        error ("pthread_mutex_lock failed in random_thread()");
#ifdef MERSENNE_TWISTER
    value = genrand_res53();
#else
    newseed[0] = 1549 * seed[0];
    newseed[1] = newseed[0] / 2048;
    newseed[0] &= 2047;
    newseed[2] = newseed[1] / 2048;
    newseed[1] &= 2047;
    newseed[1] += 1549 * seed[1] + 812 * seed[0];
    newseed[2] += newseed[1] / 2048;
    newseed[1] &= 2047;
    newseed[2] += 1549 * seed[2] + 812 * seed[1];
    swap_ptr (&newseed, &seed);
    seed[2] &= 1023;
    value = (((seed[0] / 2048.0 + seed[1]) / 2048.0 + seed[2]) / 1024.0);
#endif
    if ((pthread_mutex_unlock (&(heating_pool->random_lock))) != 0)
        error ("pthread_mutex_unlock failed in randum_thread()");
    return value;
}    /* threaded randum */
#else
MYREAL
randum (void)
/* Non thread-safe random number generator (faster) */
{
#ifdef MERSENNE_TWISTER
#ifdef MESS
    MYREAL val = genrand_res53();
    fprintf(stdout,"%i> R:%f\n",myID,val);
    return val;
#else
#ifdef SNOWLEOPARD
    //    MYREAL r;
    // dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
    //dispatch_async(queue,
    //		   ^{     
		     MYREAL r = genrand_res53();
		     //		   } 
		     //);
    return r;
#else
    MYREAL r = genrand_res53();
    return r;
#endif
#endif
#else
MYREAL value;
newseed[0] = 1549 * seed[0];
newseed[1] = newseed[0] / 2048;
newseed[0] &= 2047;
newseed[2] = newseed[1] / 2048;
newseed[1] &= 2047;
newseed[1] += 1549 * seed[1] + 812 * seed[0];
newseed[2] += newseed[1] / 2048;
newseed[1] &= 2047;
newseed[2] += 1549 * seed[2] + 812 * seed[1];
swap_ptr (&newseed, &seed);
seed[2] &= 1023;
value = (((seed[0] / 2048.0 + seed[1]) / 2048.0 + seed[2]) / 1024.0);
return value;
#endif /*mersenne_twister*/
}    /* randum */
#endif /*phtreads*/
