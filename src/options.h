/*! \file options.h */
#ifndef OPTIONS_INCLUDE
#define OPTIONS_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 O P T I O N S   R O U T I N E S 
 
 creates options structures,
 reads options from parmfile if present
 set first parameter for mcmc run,
 prints options,
 and finally helps to destroy itself.
                                                                                                               
 Peter Beerli 1996, Seattle
    updated   2009, Tallahassee
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2009 Peter Beerli, Tallahassee FL
 
$Id: options.h 1576 2009-09-09 18:22:38Z beerli $
-------------------------------------------------------*/

#include "migration.h"
///
/// initializes the memory for the option structure
/// \param **options {pointer to the optionstructure pointer}
/// \returns None
void create_options (option_fmt ** options);
///
/// populates the content of the option structure with default values 
/// \param options  {is the option structure that holds all general option data}
/// \returns None
void init_options (option_fmt * options);

///
/// reads the options from the parmfile, the parmfile name default is "parmfile"
/// but this can be overidden using the arguments to the main(). The parmfile is already
/// open and has a pointer to options->parmfile.
/// \param options  {is the option structure that holds all general option data}
/// \returns None
void read_options_master (option_fmt * options);
///
/// reads the options from the buffer that was filled by the master in save_options_buffer().
/// \param buffer { pointer to the buffer string that contains all parameter values 
///                 this buffer is filled by the function save_options_buffer() } 
/// \param options  {is the option structure that holds all general option data}
/// \returns None
void read_options_worker (char **buffer, option_fmt * options);
///
/// reads the custom migration matrix from the parmfile
/// \param file { pointer to the parmfile}  
/// \param options  {is the option structure that holds all general option data}
/// \param value {value string holding }
/// \param customnumpop {number of populations specified in the custom-migration string read from the parmfile}
/// \returns None
void read_custom_migration (FILE * file, option_fmt * options, char *value, long customnumpop);

long scan_connect (char *custm2, long start, long stop, int check);
void fillup_custm (long len, world_fmt * world, option_fmt * options);

void set_param (world_fmt * world, data_fmt * data, option_fmt * options, long locus);
void set_profile_options (option_fmt * options);
void decide_plot (worldoption_fmt * options, long chain, long chains, char type);
void set_plot (option_fmt * options);

long save_options_buffer (char **buffer, long *allocbufsize, option_fmt * options, data_fmt *data);
long save_mu_rates_buffer (char **buffer, long *allocbufsize, option_fmt * options);
long save_parmfile (option_fmt * options, world_fmt * world, data_fmt *data);

void print_menu_options (world_fmt * world, option_fmt * options, data_fmt * data);
void print_options (FILE * file, world_fmt * world, option_fmt * options, data_fmt * data);

///
/// prints the value of minimum value of the prior distribution
/// \param tmp {string to hold the minimum value of the prior distribution}
/// \param prior {pointer to the prior structure}
/// \param priortype {indicator of the priortype}
/// \returns {pointer to tmp}
char * show_priormin(char *tmp, prior_fmt *prior, int priortype);
///
/// prints the value of the mean of the prior distribution
/// \param tmp {string to hold the mean value of the prior distribution}
/// \param prior {pointer to the prior structure}/// \param priortype {indicator of the priortype}
/// \returns {pointer to tmp}
char * show_priormean(char *tmp, prior_fmt *prior, int priortype);
///
/// prints the value of maximum value of the prior distribution
/// \param tmp {string to hold the maximum value of the prior distribution}
/// \param prior {pointer to the prior structure}
/// \param priortype {indicator of the priortype}
/// \returns {pointer to tmp}
char * show_priormax(char *tmp, prior_fmt *prior, int priortype);
///
/// prints the value of delta used for the prior distribution
/// \param tmp {string to hold the delta value used for the prior distribution}
/// \param prior {pointer to the prior structure}
/// \param priortype {indicator of the priortype}
/// \returns {pointer to tmp}
char * show_priordelta(char *tmp, prior_fmt *prior, int priortype);
///
/// prints the value of the number of bins of the posterior and prior distribution
/// \param tmp {string to hold the number of bins for recording of the posterior}
/// \param prior {pointer to the prior structure}
/// \param priortype {indicator of the priortype}
/// \returns {pointer to tmp}
char * show_priorbins(char *tmp, prior_fmt *prior, int priortype);
///
/// adjust the parameter values according to the custom migration matrix
/// to set the parameters to specific values
/// that allow for parameter reduction by using only one average population size
/// or one average migration rate, or symmetric migration rates [symmetric m/mu
/// or symmetric Nm values] or constant parameter values that do not change during the run.
/// \param world {holds all runtime related material including the data and run-options structures}
/// \param options  {is the option structure that holds all general option data}
/// \returns {None}
void synchronize_param (world_fmt * world, option_fmt * options);
///
/// resynchronizes the parameters with the custm migraiton matrix but only recalculates the material that was 
/// changed, is supposed to do less work than the synchronize_param() function.
/// \param world {holds all runtime related material including the data and run-options structures}
/// \returns {None}
void resynchronize_param (world_fmt * world);
///
/// Reorders and relabels the locations(populations), so that the datafile may contain n number of locations that
/// can be combined into m number of new locations that are then used for the anlysis.
/// Examples:\n
/// \verbatim
/// (1,2,3,4,5) --> (1,1,1,2,2)
/// (1,2,3,4,5) --> (2,3,1,1,1)
/// \endverbatim
/// The reordering parameters are held in options->newposts and its size is data->numpop 
/// \param world {holds all runtime related material including the data and run-options structures}
/// \param options  {is the option structure that holds all general option data}
/// \param data  {holds all data related parts}
/// \returns {None}
void reorder_populations(world_fmt *world, option_fmt *options, data_fmt *data);
/// 
/// fills the option->newposts array with localities to reorder 
void set_localities(char **value, char **tmp, option_fmt *options);
///
/// sets the average of the mutation rates among loci and updates option structures
/// in world->options and options.
/// \param wopt  {is the lean local copy of options in world}
/// \param options  {is the option structure that holds all general option data}
/// \param loci {number of loci}
/// \returns {None}
void set_meanmu(worldoption_fmt * wopt, option_fmt * options, long loci);

#endif /*OPTIONS_INCLUDE */
