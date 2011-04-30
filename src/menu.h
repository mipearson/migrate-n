#ifndef MENU_INCLUDE
#define MENU_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 M E N U   R O U T I N E S 
 
 presents the menu and its submenus.                                                                                                               
 Peter Beerli 1996, Seattle
      updated 2009, Tallahassee
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2009 Peter Beerli, Tallahassee FL
 
$Id: menu.h 1546 2009-06-01 16:01:52Z beerli $
-------------------------------------------------------*/

#include "migration.h"
///
/// print the title of the analysis to the log and the screen
/// \param file {stream pointer to which the title string is printed to}
/// \param options  {is the option structure that holds all general option data}
/// \returns {None}
void print_menu_title (FILE * file, option_fmt * options);
///
/// print the acceptance ratio to the screen
/// \param a {likelihood of the newly proposed MCMC state}
/// \param b {likelihood of the old MCMC state
/// \param world {holds all runtime related material including the data and run-options structures}
/// \returns {None}
void print_menu_accratio (long a, long b, world_fmt *world);
///
/// print the title to the outfile
/// \param world {holds all runtime related material including the data and run-options structures}
/// \param options  {is the option structure that holds all general option data}
/// \returns {the position in the output file}
long print_title (world_fmt * world, option_fmt * options);
///
/// presents the menu to the user
/// \param options {is the option structure that holds all general option data}
/// \param world {holds all runtime related material including the data and run-options structures, 
///               but this is not really filled yet because all data related material is missing}
/// \param data  {holds all data related parts}
/// \returns {None}
void get_menu (option_fmt * options, world_fmt *world, data_fmt *data);
///
/// presents title of the prior method
/// \param priorset {is an int that specifies the prior type}
/// \returns {a string with the name of the proposal method}
char * is_priortype(int priorset);
///
/// presents title of the posterior method
/// \param proposalset {is a boolean to specify whether this is SPLICE or MH sampling scheme}
/// \returns {a string with the name of the proposal method}
char * is_proposaltype(boolean proposalset);
#endif /*MENU_INCLUDE */
