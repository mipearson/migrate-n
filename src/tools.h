#ifndef TOOLS_INCLUDE
#define TOOLS_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                       
 
 H E L P E R     R O U T I N E S 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2006 Peter Beerli, Tallahassee FL
 
$Id: tools.h 1796 2011-01-25 08:33:50Z beerli $
-------------------------------------------------------*/

#include "migration.h"

extern MYREAL lengthof (node * p);
extern MYINLINE node *crawlback (const node * theNode);
extern node *crawl (node * theNode);
extern MYINLINE node *showtop (node * theNode);
extern void adjust_time (node * theNode, MYREAL tyme);
extern void adjust_time_all (node * theNode, MYREAL tyme);
extern void insert_migr_node (world_fmt * world, node * up, node * down,
                                  migr_table_fmt * migr_table,
                                  long *migr_table_counter);
extern void children (node * mother, node ** brother, node ** sister);
/* math tools */
extern MYREAL incompletegamma (MYREAL tx, MYREAL talpha);
extern MYREAL polygamma (long n, MYREAL z);
extern void invert_matrix (MYREAL **a, long nsize);
extern boolean nrcheck (MYREAL **m, MYREAL **tm, MYREAL *v, long nrows,
                            MYREAL *r1, MYREAL *r2, boolean do_newton);
extern MYREAL rannor (MYREAL mean, MYREAL sd);
extern MYREAL find_chi (long df, MYREAL prob);
extern MYREAL probchiboundary (MYREAL chi, long zeros, long all);
extern MYREAL chiboundary (long zeros, long nonzeros, MYREAL alpha);
extern MYREAL probchi (long df, MYREAL chi);
extern MYREAL chisquare (long df, MYREAL alpha);
extern void gamma_rates (MYREAL *rate, MYREAL *probcat, long categs,
                             char *input);
#ifndef HAVE_LGAMMA
extern MYREAL mylgamma (MYREAL z);
#endif
extern MYREAL calc_sum (MYREAL *vector, long n);
extern MYREAL logfac (long n);

extern void onepass_mean_std_start(MYREAL *mean, MYREAL *std, long *n);
extern void onepass_mean_std_calc(MYREAL *mean, MYREAL *std, long *n, MYREAL x);
extern void onepass_mean_std_end(MYREAL *mean, MYREAL *std, long *n);

extern boolean traverse_check_id(node *theNode, long id);
extern  void start_node_collection(world_fmt *world);
extern void stop_node_collection(world_fmt *world);
extern void collect_nodelet(world_fmt *world, node *p);
extern node *  dispense_nodelet(world_fmt *world);
extern void swap_node_collection(world_fmt * tthis, world_fmt * that);


extern float fast_log(float vval);
extern MYINLINE double fast_exp(double y); 


extern char lowercase (int c);
extern char uppercase (int c);
extern void upper(char *from, char **to);
extern long read_word_delim(char *input, char *word, char * delim);

/* vector initialization and calc*/
extern void doublevec1d (MYREAL **v, long size);
extern void doublevec2d (MYREAL ***v, long size1, long size2);
extern void floatvec2d (float ***v, long size1, long size2);
extern void charvec2d (char ***v, long size1, long size2);
extern void intvec2d (long ***v, long size1, long size2);
extern void free_doublevec2d (MYREAL **v);
extern void free_floatvec2d (float **v);
extern void free_charvec2d (char **v);
extern void free_intvec2d (long **v);
extern void setdoublevec1d (MYREAL **v, MYREAL *w, long size);
extern void add_vector (MYREAL *result, MYREAL *v, long size);

/*filemanipulation */
extern void init_files (world_fmt * world, data_fmt * data,
                            option_fmt * options);
extern void exit_files (world_fmt * world, data_fmt * data,
                            option_fmt * options);
extern void openfile (FILE ** fp, char *filename, char *mode, char *perm);
#ifdef ZNZ
extern void znzopenfile (znzFile * fp, char *filename, char *mode, int use_compressed);
#endif
extern long read_savesum (world_fmt * world, option_fmt * options,
                              data_fmt * data);
extern void write_savesum (world_fmt * world);

/* string manipulation */
extern void translate (char *text, char from, char to);
extern void unpad (char *text, char removechars[]);
extern void get_next_word(char **instring, char *delimiters, char **nextword);
extern long count_words (char *text);
extern void read_word(FILE *infile, char *word);
extern void unread_word(FILE *infile, char *word);
extern void add_to_buffer(char *fp, long *bufsize, char **buffer, long *allocbufsize);
extern long print_to_buffer(char **buffer, long *maxbufsize, char *tempbuffer, long *pos, const char *fmt, ...);
extern void record_warnings(world_fmt * world, const char *fmt, ...);
extern void print_stored_warnings(world_fmt *world);
extern void remove_trailing_blanks(char **filename);

/* time reporting */
extern void get_time (char *nowstr, char ts[]);
/* printing aid */
extern void print_llike (MYREAL llike, char *strllike);
/* searching and finding*/
extern boolean find (long i, long *list, long listlen);
/* conversion between the different parameter schemes*/
extern long mstart (long pop, long numpop);
extern long mend (long pop, long numpop);
extern long mm2m (long frompop, long topop, long numpop);
extern void m2mm (long i, long numpop, long *frompop, long *topop);
extern long m2mml (long i, long numpop);
extern long m2mml2 (long i, long topop, long numpop);
extern long mmstart (long pop, long numpop);
extern long mmend (long pop, long numpop);
extern long mml2m (long pos, long numpop);
extern void print_line (FILE * outfile, char c, long nn, long flag);
extern void sprint_line (char *buffer, char c, long nn, long flag);
extern void fprintf2(FILE *file, long filesize, const char *fmt, ...);

/* reading from char * buffer */
char sgetc (char **buffer);
char *sgets (char *s, int size, char **stream);
#ifdef CAUTIOUS
extern boolean cautious;
#endif

#endif /*TOOLS_INCLUDE */
