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
 
$Id: profile.h 1795 2011-01-16 03:21:11Z beerli $
-------------------------------------------------------*/
extern boolean print_profile_likelihood_driver (long which, world_fmt * world,
            long *gmaxptr);
extern void print_profile_likelihood (long which, world_fmt * world, long *gmaxptr);
extern long print_profile_percentile (world_fmt * world);
extern void allocate_profile_percentiles (world_fmt * world);
extern void destroy_profile_percentiles (world_fmt * world);
extern void print_profile_title (world_fmt * world);
extern long warp (long ii);

#define GRIDSIZE 9
#define GRIDMIDDLE 5
#define GRID    {0.01,0.05,0.10,0.50,0.99,0.95,0.90,0.50,1.0}
#define SHOWGRID {0.005,0.025,0.05,0.25,0.995,0.975,0.95,0.75,0.50}
#define GRID2   {0.005,0.025,0.05,0.25,0.5, 0.75,0.95,0.975,0.995}
#define INDEX {0,1,2,3,8,7,6,5,4}
#define DEVIATE {0.02,0.10,0.20, 0.5, 50.,10., 5., 2., 1.}
#define DEVIATE2 {0.02,0.10,0.20, 0.5, 1., 2., 5., 10., 50. }
//#define DEVIATE {0.002,0.010,0.020, 0.05, 5000.,1000., 500., 200., 1.}
//#define ABSOLUTE {1e-100,1e100}
