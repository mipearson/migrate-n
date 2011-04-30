/*! \File migrate_mpi.c */
/* MPI parts for migrate
   started November 2000, Seattle
   Peter Beerli beerli@fsu.edu
 
   
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2007 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
$Id: migrate_mpi.c 1861 2011-04-06 21:14:54Z beerli $
*/
#ifdef MPI
#include "migration.h"
#include "tools.h"
#include "sighandler.h"
#include "migrate_mpi.h"
#include "broyden.h"
#include "combroyden.h"
#include "gammalike.h"
#include "profile.h"
#include "pretty.h"
#include "options.h"
#include "tree.h"
#include "world.h"
#include "joint-chains.h"
#include "data.h"
#include "laguerre.h"
#include "random.h"
#ifdef UEP
#include "uep.h"
#endif
#include "bayes.h"
#ifndef WINDOWS
#include <unistd.h>
#endif
/*should go into profile.h*/
#define GRIDSIZE 9

#ifdef PRETTY
extern float page_height;
#endif

extern const MPI_Datatype mpisizeof;

extern void run_replicate(long locus,
                          long replicate,
                          world_fmt **universe,
                          option_fmt *options,
                          data_fmt *data, 
                          tpool_t * heating_pool,
                          int usize,
                          long *treefilepos,
                          long *Gmax);
extern void run_locus (world_fmt ** universe, int usize,
                       option_fmt * options, data_fmt * data,
                       tpool_t * heating_pool, long maxreplicate,
                       long locus, long *treefilepos, long *Gmax);

void mpi_run_locus(world_fmt ** universe, int usize, option_fmt * options,
                   data_fmt * data, tpool_t * heating_pool, long maxreplicate,
                   long locus, long *treefilepos, long *Gmax);
void mpi_runreplicates_worker (world_fmt ** universe, int usize,
                               option_fmt * options, data_fmt * data,
                               tpool_t * heating_pool,
                               long *treefilepos, long *Gmax);
long pack_databuffer (char **buffer, data_fmt * data, option_fmt * options);
void unpack_databuffer (char *buffer, data_fmt * data, option_fmt * options);
void pack_allele_data (char **buffer, long *bufsize, data_fmt * data,
                       long pop, long ind);
void pack_sequence_data (char **buffer, long *bufsize, data_fmt * data,
                         long pop, long ind, long locus);
void mpi_gradient_master (nr_fmt * nr, world_fmt * world, int *who);
void mpi_resultsmaster (long sendtype, world_fmt * world,
                         long maxrep,
                         void (*unpack) (char *buffer, world_fmt * world,
                                         long locus, long maxrep,
                                         long numpop));

void mpi_results_worker (long bufs, world_fmt * world,
                         long maxrep,
                         long (*pack) (MYREAL **buffer, world_fmt * world,
                                       long locus, long maxrep, long numpop));
void assignloci_worker (world_fmt * world, option_fmt *options, long *Gmax);
void assign_worker_cleanup (void);
void swap_atl (long from, long to, world_fmt * world);
long pack_quantile (char **buffer, quantile_fmt quant, long n);
void unpack_quantile (char *buffer, quantile_fmt quant, long n);
long pack_failed_percentiles (char **buffer, boolean *failed, long n);
void unpack_failed_percentiles (char *buffer, boolean *failed, long n);

void handle_message(char *rawmessage,int sender, world_fmt *world);
void handle_mdim(float *values,long n, int sender, world_fmt * world);
void handle_burnin_message(char *rawmessage,int sender, world_fmt * world);

void set_filehandle(char *message, world_fmt *world,
                    void **file, long *msgstart);

void mpi_receive_replicate( int sender, int tag, long locus, long replicate, world_fmt * world);

long unpack_single_bayes_buffer(MYREAL *buffer, bayes_fmt * bayes, world_fmt * world,long locus);
long pack_single_bayes_buffer(MYREAL **buffer, bayes_fmt *bayes, world_fmt *world,long locus);
long pack_single_bayes_buffer_part(MYREAL **buffer, bayes_fmt *bayes, world_fmt *world,long locus);

void unpack_hist_bayes_buffer(MYREAL *buffer, bayes_fmt *bayes, world_fmt *world, long locus);
long pack_hist_bayes_buffer(MYREAL **buffer, bayeshistogram_fmt *hist, world_fmt * world, long startposition);

long unpack_BF_buffer(MYREAL *buffer, long start, long locus, world_fmt * world);
long unpack_ess_buffer(MYREAL *buffer, long start, world_fmt *world);
long pack_BF_buffer(MYREAL **buffer, long start, long locus, world_fmt * world);
long pack_ess_buffer(MYREAL **buffer, long start, world_fmt *world);

void unpack_sumfile_buffer (MYREAL *buffer, world_fmt * world,
                            long locus, long maxrep, long numpop);
void unpack_single_sumfile_buffer (MYREAL *buffer, timearchive_fmt **ta, world_fmt *world,
                                   long locus, long replicate, long numpop, long *startz);
long
pack_sumfile_buffer (MYREAL **buffer, world_fmt * world,
                     long locus, long maxrep, long numpop);

long pack_single_sumfile_buffer(MYREAL **buffer, long z, world_fmt * world,
                                long locus, long replicate, long numpop);

void mpi_send_replicate(int sender, long locus, long replicate, world_fmt * world);
long  mpi_send_stop_mcmc_lociworker(long numcpu, long loci);
long  mpi_send_stop_mcmc_replicateworker(long numcpu, long loci);
void mpi_send_stop_tag (int worker, world_fmt * world);
long  mpi_send_stop_mcmc_worker(long numcpu, long loci, MPI_Comm *comm, MPI_Request *irequests, MPI_Status *istatus, long id);
void mpi_send_stop_tag (int worker, world_fmt * world);
void send_receive_bayes_params(world_fmt *world, long locus);
void handle_replication(int sender,int tag,char *tempstr, world_fmt *world);
boolean in_mpistack(int sender, world_fmt *world);

#ifdef MPICHECK
#include <sys/resource.h>
void set_memory_limit(rlim_t softsize,rlim_t maxsize);
void check_memory_limit();
#endif
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

///
/// Controls all loci in the MPI implementation, uses a load balancing scheme and 
/// distributes the work on the waiting nodes
void
mpi_runloci_master (long loci, int *who, world_fmt *world, boolean options_readsum, boolean menu)
{
    boolean done = FALSE;
    long alldone = 0;
    int tag;
    int sender       = 0;

    long locus;
    long *twolongs;
    long locusdone   = -1;
    long numsent     = 0;
    long tempstrsize = LINESIZE;
    long nbase       = loci + 1;
    long minnodes    = MIN((long) numcpu-1, (long) nbase-1);

    char *tempstr; 
    char *leadstr;
    float *temp;   
    long nn = world->numpop2 + 9 + world->bayes->mu + world->options->heated_chains + 2;

    MPI_Status status;
    MPI_Status *istatus;
    MPI_Request *irequests;

    int ll;
#ifdef IPROBE
    int notwaiting=0;
#endif
    irequests  = (MPI_Request *) mycalloc(minnodes, sizeof(MPI_Request));
    istatus    = (MPI_Status *) mycalloc(minnodes, sizeof(MPI_Status));
    twolongs   = (long *) mycalloc(TWO, sizeof(long));
    tempstr    = (char *) mycalloc(tempstrsize, sizeof(char));
    leadstr    = (char *) mycalloc(SMALLBUFSIZE, sizeof(char));
    temp       = (float *) mycalloc(nn, sizeof(float));
    twolongs[1]= 0;
    //    	printf("%i> MASTER: minnodes = %li\n+++++++++++++++++++++\n",myID, minnodes);
    for (locus = 0; locus < minnodes; locus++)
      {
	twolongs[0] = locus;
	ll = (int) (locus + 1);
	MYMPIISEND(twolongs, TWO, MPI_LONG, ll, ll, comm_world, &irequests[numsent]);
	numsent++;
	//	printf("%i> MASTER: Locus %li sent to worker n%i twolong[0]=%li\n",myID, locus, (MYINT) locus + 1, twolongs[0]);
      }
    // waits until all minodes nodes received their message
    MYMPIWAITALL(minnodes, irequests, istatus);
    //printf("%i> MASTER: after initial send of loci\n",myID);

    locus = 0;   
    world->mpistacknum = 0;
    while((alldone < world->loci) || (world->mpistacknum < numcpu-1))
      {
        done=FALSE;
        while(!done)
	  {
	    memset(tempstr,0,sizeof(char)*tempstrsize);
#ifdef IPROBE
	    if(menu)
	      {
		MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm_world, &notwaiting, &status);
		if(!notwaiting)
		  {
		    sleep(1);
		    get_time (tempstr, "%H:%M:%S");
		    fprintf(stdout,"%i> master waits for workers -- %s\n",myID, tempstr);
		    continue;
		  }
	      }
#endif
	    MYMPIRECV (leadstr, SMALLBUFSIZE, MPI_CHAR, (MYINT) MPI_ANY_SOURCE, 
		       (MYINT) MPI_ANY_TAG, comm_world, &status);
	    sender = status.MPI_SOURCE;
	    tag = status.MPI_TAG;
	    //if(tempstr[0] != 'M')
	      //printff("%i> MASTER: received from n%i with tag %i: %s\n",myID,sender,tag,tempstr);
	    switch(leadstr[0])
	      {
	      case 'M':
		tempstrsize = atol(leadstr+1);
		tempstr = (char*) realloc(tempstr,sizeof(char)*(tempstrsize+1));
		memset(tempstr,0,sizeof(char)*(tempstrsize+1));
		MYMPIRECV (tempstr, tempstrsize, MPI_CHAR, sender, tag,
			   comm_world, &status);
		//		printf("%i> MASTER: message from n%i: %s\n",myID,sender,tempstr);
		handle_message(tempstr,sender, world);
		break;
	      case 'Z': // reading/writing of the raw bayes posterior data , also used to guide
		// multiple replicates on when to start sampling
		MYMPIRECV (temp, nn, MPI_FLOAT, sender, tag,
			   comm_world, &status);
		handle_mdim(temp,nn,sender, world);
		break;
	      case 'B': //burnin stopping rule
		tempstrsize = atol(leadstr+1);
		tempstr = (char*) realloc(tempstr,sizeof(char)*(tempstrsize+1));
		memset(tempstr,0,sizeof(char)*(tempstrsize+1));
		MYMPIRECV (tempstr, tempstrsize, MPI_CHAR, sender, tag,
			   comm_world, &status);
		//printf("%i> MASTER: burn-in message from n%i: %s\n",myID,sender,tempstr);
		handle_burnin_message(tempstr,sender-BURNTAG, world);
		break;
	      case 'R':
		//ignore first character and translate into locusnumber
                locusdone = atol(leadstr+1);
		//if negative this means locus had no data at all -> ignored 
		if(locusdone<0)
		  {
		    locusdone = -locusdone;
		    world->data->skiploci[locusdone]=TRUE;
		    //printf("MASTER received a skiplocus %li\n",locusdone);
		  }
		done=TRUE;
		++alldone;
		//printff("%i> MASTER: locus %li, %li of %li finished **************************\n",myID, locusdone,alldone,world->loci);
                break;
	      case 'N': // need a replicator, this message is sent from a
		// locus-node to the master for distribution among nodes that
		// are waiting for work, the master send will delegate a 
		// replicate and tell also the replicator 
		// where it needs to send the final result.
#ifdef DEBUG_MPI
	        printf("%i> MASTER: received replicator request from n%i using string %s\n",myID, sender, tempstr);
#endif
		// add the mpistack_request list
		handle_replication(sender,tag,leadstr,world);
		break;
	      case 'G':
		// node has free time and could do replicator work
#ifdef DEBUG_MPI
	        printf("%i> MASTER: accepts replicator n%i at stack position %li\n",myID, sender,world->mpistacknum);
#endif
		world->mpistack[world->mpistacknum] = sender;
		world->mpistacknum += 1;
		if(alldone>=world->loci && (world->mpistacknum >= (numcpu-1)) && world->mpistack_requestnum==0)
		  {
#ifdef DEBUG_MPI
		    printf("%i> MASTER: reached mpistacknum=%li\n",myID,world->mpistacknum);
#endif
		    done=TRUE;
		    continue;
		  }
		while(world->mpistacknum > 0 && world->mpistack_requestnum > 0)
		  {
		    world->mpistack_requestnum -= 1;
		    sender = world->mpistack_request[world->mpistack_requestnum].sender;
		    tag = world->mpistack_request[world->mpistack_requestnum].tag;
		    handle_replication(sender,tag,
		      world->mpistack_request[world->mpistack_requestnum].tempstr,world);
		  }
		break;
	      default: /* never go here under normal run condition */
#ifdef DEBUG_MPI
                fprintf(stderr,"%i> message=@%s@\n@%i@> sender=@%i@ tag=@%i@\n",
			myID,leadstr, myID,status.MPI_SOURCE,status.MPI_TAG);
		done=TRUE;
#else
                MPI_Finalize();
                error("DIED because of wrong message from worker");
#endif
                break;
	      }
	  }
        who[locusdone] = sender;
	// if not done with loci send another locus-work to sender (a node 1..numcpu)
        if (numsent < loci)
	  {
            twolongs[0]=numsent;
            MYMPISEND (twolongs, TWO, MPI_LONG, (MYINT) sender, (MYINT) numsent + 1, comm_world);
            numsent++;
	//printff("%i> MASTER: ANOTHER locus %li sent to worker n%i %i\n",myID, numsent, (MYINT) numsent, (MYINT) numsent + 1);
	  }
        else
	  {
            twolongs[0] = 0;
	    //tell workers to stop waiting for new loci
	    if(!in_mpistack(sender, world))
	      {
		MYMPISEND (twolongs, TWO, MPI_LONG, (MYINT) sender, (MYINT) 0, comm_world); 
	      }
	  }
	//printff("%03i> available workers: %li\n     requests pending:      %li\n     total workers=%i\n     loci done=%li\n",myID, world->mpistacknum,world->mpistack_requestnum, numcpu-1, locus+1);
	locus++;
      }
//printf("%03i> available workers: %li\n     requests pending:      %li\n     total workers=%i\n     loci done=%li\n   alldone=%li\n",myID, world->mpistacknum,world->mpistack_requestnum, numcpu-1, locus+1,alldone);
    // stop loci and/or replicate worker that had never the chance to work on a locus 
    // or replicate, but are still listening
    for(locus=0;locus < world->mpistacknum;locus++)
      {
#ifdef DEBUG_MPI        
	fprintf(stdout,"%i> sent kill to node %i\n",myID, world->mpistack[locus]);
#endif
	mpi_send_stop_tag(world->mpistack[locus], world);
      }
    myfree(twolongs);
    myfree(istatus);
    myfree(irequests);
    myfree(tempstr);
    myfree(leadstr);
    myfree(temp);
    //printf("%i> before barrier\n",myID);
    //MYMPIBARRIER(comm_world);
    //printf("%i> leaving mpi_runloci_function()\n", myID);
}


///
/// worker nodes execute this function and serve the master
void mpi_runloci_worker (world_fmt ** universe, int usize,
                    option_fmt * options, data_fmt * data,
                    tpool_t * heating_pool, long maxreplicate,
                    long *treefilepos, long *Gmax)
{
    boolean done = FALSE;

    char *rawmessage;

    long locus;
    long *twolongs;
    long rawmsgsize    = 0;
#ifdef MPIREPLICANT    
    long nbase         = data->loci+1;
#endif

    MPI_Status status;    

    twolongs   = (long *) mycalloc(TWO, sizeof(long));
    rawmessage = (char *) mycalloc(SMALLBUFSIZE,sizeof(char));
#ifdef MPIREPLICANT
    if(myID < nbase)
      {
#endif
        while (!done)
	  {
            MYMPIRECV (twolongs, TWO, MPI_LONG, MASTER, MPI_ANY_TAG,
		       comm_world, &status);
            locus = twolongs[0];
            if (status.MPI_TAG != 0) //stop condition
	      {
		//printf("%i> WORKER: received locus %li (twolongs[0]=%li) from MASTER\n",myID, locus, twolongs[0]);
#ifdef MPIREPLICANT
                mpi_run_locus(universe, usize, options, data, 
                              heating_pool, maxreplicate, locus, treefilepos, Gmax);  
#else
                run_locus (universe, usize, options, data,
                           heating_pool, maxreplicate, locus, treefilepos, Gmax);
#endif
		if(universe[0]->data->skiploci[locus]==FALSE)
		  rawmsgsize = 1 + sprintf(rawmessage,"R%li",locus);
		else
		  rawmsgsize = 1 + sprintf(rawmessage,"R-%li",locus);
                MYMPISEND (rawmessage, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, 
			   (MYINT) (locus + ONE), comm_world);
                /* we want to know what locus we worked for
                   - to control the work sent by master
                   - to use in setup_parameter0() [combroyden2.c] */
                universe[0]->who[locidone++] = locus;
	      }
            else
	      {
                done = TRUE;
		//printf("%i> LOCUSWORKER finished\n",myID);
		mpi_runreplicates_worker (universe, usize, options,  data, heating_pool, treefilepos, Gmax);
		//		printf("%i> LOCUSWORKER and then REPLICANT: finished\n",myID);
	      }
	  }
#ifdef MPIREPLICANT
      }
    else
      {
        mpi_runreplicates_worker (universe, usize, options,  data, heating_pool, treefilepos, Gmax);
      }        
#endif
    myfree(twolongs);
    myfree(rawmessage);
    //printf("%i> WORKER done\n",myID);
    //   printf("%i> before barrier\n",myID);
    //MYMPIBARRIER(comm_world);
}

#ifdef MPIREPLICANT
///
/// attempt to improve MPIruns with replicates
/// each locus is responsible for replication farming and and reporting back to master
/// master <- locus-master <- locus-replicate-worker
/// generate genealogies for a single locus
/// \callgraph
void
mpi_run_locus(world_fmt ** universe, int usize, option_fmt * options,
          data_fmt * data, tpool_t * heating_pool, long maxreplicate,
          long locus, long *treefilepos, long *Gmax)
{
    
    boolean done=FALSE;
    
    char *tempstr;
    //    const long tempstrsize = SMALLBUFSIZE;
    int tag;
    int *who;
    int minnodes;                         // number of replicate-worker nodes [changed meaning!]
    int sender           = 0;
    int nbase            = data->loci;    //number of loci-worker nodes
    int senderlocus      = -1;
    long replicate;
    long i;
    long *temp;
    long numsent         = 0;
    long senderreplicate = 0;

    //    MPI_Request irequest;
    MPI_Request *irequests = NULL;

    MPI_Status status;
    MPI_Status *istatus    = NULL;
    
    temp    = (long *) mycalloc(TWO, sizeof(long));
    tempstr = (char *) mycalloc(SMALLBUFSIZE,sizeof(char));
    who     = (int *) mycalloc(maxreplicate,sizeof(int));

    
    if(maxreplicate>1)
      {
	// number of nodes available for replicate workers,
	// numcpu - nbase - 1 [for masternode] and maxreplicate-1 
	// [locus-worker is doing one replicate itself]  are limiting
	//minnodes =  maxreplicate; //MAX(0,MIN(maxreplicate-1,numcpu-nbase-1))+1;
	minnodes =  MAX(0,MIN(numcpu-nbase-1,maxreplicate-1))+1;
	irequests = (MPI_Request *) mycalloc(minnodes+1,sizeof(MPI_Request));
	istatus   = (MPI_Status *) mycalloc(minnodes+1,sizeof(MPI_Status));
	temp[0] = locus;
	// replicate 1 to maxreplicate should be worked on by other nodes
	// minnodes is the number of nodes free for work on replicates alone
	// so with 3 replicates to work on and 3 total nodes and a single locus we need
	// 1 master, 1 locus-worker and have 1 replicate-worker, the locus-worker and 
	// the replicate worker will both work on replicates, the locus-worker sends 
	// off 1 request (for the replicate worker), because we start the locus worker 
	// with replicate 0, the loop for the replicate workers starts at one
	// and to address this shift the MIN(..,minnodes+1) addresses that
	for (replicate = 1; replicate < minnodes; replicate++)
	  {
            temp[1] = replicate;
	    sprintf(tempstr,"N%i %li %li", myID, locus, replicate);
	    //fprintf(stdout,"%i> request a node from n%i for locus %li and replicate %li\n",myID, MASTER, locus, replicate);
            MYMPIISEND (tempstr, SMALLBUFSIZE, MPI_CHAR, 
			(MYINT) MASTER, (MYINT) (locus + 1 + REPTAG), comm_world, &irequests[numsent]);
            numsent++;   // counter of how many replicates are sent off-node
	  }
        // this should be asynchronous so that the myRepID-master can do some work, too.
        run_replicate(locus, 0, universe, options, data, 
                      heating_pool, usize,
                      treefilepos, Gmax);
        who[0] = myID;
        MYMPIWAITALL(numsent,irequests, istatus); // wait for all replicators to finish
        // set replicate counter to 1 because locus-worker itself finished first replicate
        replicate=1;
        numsent++;
        
        done = FALSE;
        while(!done)
	  {
            // done=TRUE means that
            // no replicator worker is available yet
            // the loci-worker has to do all the work if numsent==1
            if(numsent==1)
	      {
                run_replicate(locus, replicate, universe, options, data, 
                              heating_pool, usize,
                              treefilepos, Gmax);
                who[replicate] = myID;
                replicate++;
                if(replicate >= maxreplicate)
		  done=TRUE;
	      }            
            else
	      {
                memset(irequests,0,sizeof(int)*minnodes);
                memset(istatus,0,sizeof(int)*minnodes);
                MYMPIRECV (tempstr, SMALLBUFSIZE, MPI_CHAR, MPI_ANY_SOURCE, 
			   (MYINT)(locus+1+ REPTAG), comm_world, &status);
                sender = status.MPI_SOURCE;  // repID of the replicator that did the work
                tag = status.MPI_TAG;        // tag is working locus + replicator tag 
                senderlocus = tag-REPTAG;    // locus that was worked on by sender  
                // test so that we can be sure that we got things from a valid replicator worker
                if(sender == myID)
		  {
                    fprintf(stdout,"%i, %i> DIE\nDIE\nDIE\nDIE\nDIE\nDIE\nDIE\nDIE\nDIE\n",myID, myRepID);
                    error("tried to send a replicate to myself using MPI -- this is not allowed\n");
		  }
                // test whether the locus is the same between locus-worker and replicator
                if(senderlocus-1 != locus)
		  warning("%i> !!!!!!!!!!!!! got wrong locus from worker myRepID=%i (my locus %i != its locus %i )\n",
                          myID, sender, locus, senderlocus-1);
                // receive only messages that are prefixed with 'R' and exit on all others
                if(tempstr[0]=='R')
		  {
                    //ignore first character and translate into repnumber
                    senderreplicate = atol(tempstr+1);
		  }
                else
		  {
		    fprintf(stderr,"%i> message=%s\n%i> sender=%i tag=%i\n",myID,tempstr, myID,
			    status.MPI_SOURCE,status.MPI_TAG);
                    error("DIED because of wrong message from worker");
		  }
		// record sender , this record should be filled at the end of this function
                who[senderreplicate] = sender;
		//fprintf(stdout,"%i> senderreplicate=%li\n%i> sender=%i tag=%i\n",myID,senderreplicate, myID,
		//	    status.MPI_SOURCE,status.MPI_TAG);
                mpi_receive_replicate(sender, tag, locus, senderreplicate, universe[0]); 
                replicate++;   
                if(replicate >= maxreplicate)
		  {
                    done=TRUE;
		  }
                temp[0] = locus;
                if (numsent < maxreplicate) //at least one set was worked by the locus-worker
		  {
                    temp[1] = numsent;
		    sprintf(tempstr,"N%i %li %li", myID, locus, numsent);
		    //fprintf(stdout,"%i> request a node from n%i for locus %li and replicate %li\n",myID, MASTER, locus, numsent);
		    MYMPIISEND (tempstr, SMALLBUFSIZE, MPI_CHAR, 
				(MYINT) MASTER, (MYINT) (locus + 1 + REPTAG), comm_world, &irequest);
                    numsent++;
                    MYMPIWAITALL(ONE, &irequest, &status);
		  }
	      }
	  }
      }
    else
      { /* no replicates */
        run_replicate(locus, 0, universe, options, data, 
                      heating_pool, usize,
                      treefilepos, Gmax);
      }    
    myfree(temp);
    myfree(tempstr);
    myfree(who);
    if(maxreplicate>1 && (numcpu-1) > universe[0]->loci)
      {
	const long hc = universe[0]->options->heated_chains;
	long t;
	if(universe[0]->options->heating)
	  {
	    printf("CORRECTING\n");
	    for(t=0; t < hc; t++)
	      {
		universe[0]->bf[locus * hc + t] /= maxreplicate;
	      }
	  }
	myfree(irequests);
	myfree(istatus);
      }
#ifdef UEP
    if (options->uep)
      show_uep_store(universe[0]);
#endif
    if (options->bayes_infer)
      {
	if (options->replicate && options->replicatenum > 0)
	  {
	    (universe[0])->repkind = MULTIPLERUN;
	  }
	// @@@@@@@@@@CHECK@@@@@@@@@@@@@
	if(!universe[0]->options->has_bayesmdimfile)
	  calculate_credibility_interval(universe[0], locus);
      }
    else
      {
	if (options->replicate && options->replicatenum > 0)
	  {
	    (universe[0])->repkind = MULTIPLERUN;
#ifdef LONGSUM
	    change_longsum_times(EARTH);
#endif /*LONGSUM*/        
	    // over multiple replicates if present or single locus
	    //printf("ML over multiple replicates or single locus\n");
	    (void) estimateParameter(options->replicatenum, *Gmax, universe[0], options,
				     (universe[0])->cov[locus], options->lchains, /*type*/ 'l',
				     SINGLELOCUS, (universe[0])->repkind);
	  }
      }

    // cleanup
    if (options->heating)
    {
        for (i = 0; i < options->heated_chains; i++)
        {
	  //free_tree(universe[i]->root, universe[i]);
        }
    }
    else
    {
      //      free_tree(universe[0]->root, universe[0]);        
    }
    bayes_reset(universe[0]);
}


///
/// run replicates on replicate-worker nodes
/// this function is called in main() and will work on any preset locus and any preset replicate
/// the calling function is responsible for the correct assignment to locus and replicate.
void
mpi_runreplicates_worker (world_fmt ** universe, int usize,
                    option_fmt * options, data_fmt * data,
                    tpool_t * heating_pool, 
                    long *treefilepos, long *Gmax)
{
    boolean done = FALSE;

    char *rawmessage;

    int sender;

    long *temp;
    long replicate;
    long locus;
    long rawmsgsize = 0;
    char *ready;
    MPI_Status status;
    ready = mycalloc(SMALLBUFSIZE,sizeof(char));
    temp = (long *) mycalloc(3, sizeof(long));
    rawmessage = (char *) mycalloc(STRSIZE,sizeof(char));
    sprintf(ready,"G%i",myID);
    while (!done)
      {
	// 	printf("%i> REPLICANT:  send ready message \"%s\" to n%i\n", myID, ready, MASTER);
	MYMPISEND (ready, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, myID, comm_world);
        MYMPIRECV (temp, 3, MPI_LONG, (MYINT) MASTER, MPI_ANY_TAG, comm_world, &status);
        sender = (int) temp[0];
        locus = temp[1];
        replicate = temp[2];
	//printf("%i> REPLICANT:  received locus %li and replicate %li from n%i via n%i \n", myID, locus+1, replicate+1, sender, status.MPI_SOURCE);
        if (status.MPI_TAG != 0) //stop condition
          {
            run_replicate(locus, replicate, universe, options, data, heating_pool, usize,treefilepos, Gmax);
	    //printf("%i> after runreplicate in\n",myID);
            rawmsgsize = 1 + sprintf(rawmessage,"R%li ",replicate);
            MYMPISEND (rawmessage, rawmsgsize, MPI_CHAR, (MYINT) sender, 
		       (MYINT) (locus+1+ REPTAG), comm_world);
            mpi_send_replicate(sender, locus, replicate, universe[0]);
	    //printf("%i> after mpi_send_replicate()\n",myID);
	    if(universe[0]->options->bayes_infer)
	      {
		universe[0]->bayes->numparams = 0;
	      }
          }
        else
          {
            done = TRUE;
	    //printf("%i> REPLICANT received KILL signal.\n",myID);
          }
      }
    myfree(ready);
    myfree(temp);
    myfree(rawmessage);
    //    printf("%i> REPLICANT: all replication work is done\n",myID);
}

//---------end replication in MPI
void
assign_worker_cleanup (void)
{
    boolean done = FALSE;

    char *rawmessage;

    int sender;

    long *temp;
    char *ready;
    MPI_Status status;
    ready = mycalloc(SMALLBUFSIZE,sizeof(char));
    temp = (long *) mycalloc(3, sizeof(long));
    rawmessage = (char *) mycalloc(STRSIZE,sizeof(char));
    sprintf(ready,"G%i",myID);
    while (!done)
      {
 	//printf("%i> REPLICANT:  send ready message \"%s\" to n%i\n", myID, ready, MASTER);
	MYMPISEND (ready, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, myID, comm_world);
        MYMPIRECV (temp, 3, MPI_LONG, (MYINT) MASTER, MPI_ANY_TAG, comm_world, &status);
        sender = (int) temp[0];
        if (status.MPI_TAG != 0) //stop condition
          {
	    printf("%i> received real message but do not know what to do with this\n",myID);
          }
        else
          {
            done = TRUE;
          }
      }
    myfree(ready);
    myfree(temp);
    myfree(rawmessage);
}

//---------end replication in MPI
#endif /*MPIREPLICANT*/


///
/// orchestrates the likelihood calculation 
MYREAL
mpi_likelihood_master (MYREAL *param, MYREAL *lparam,
                       world_fmt * world, nr_fmt * nr,
                       helper_fmt * helper, int *who)
{
  int tag;

  long locus;
  long worker;
  long sender;
  long addon    = 1;
  long numelem  = world->numpop2 + (world->options->gamma ? 1 : 0);
  long numelem2 = numelem * 2;
  
  MYREAL logres = 0.0;
  MYREAL *temp;
  MYREAL *tmp;

  MPI_Status status;

    
  doublevec1d(&tmp,world->loci);
  doublevec1d(&temp,numelem2+2);
  
  temp[0] = MIGMPI_LIKE;
  memcpy (temp + 1, param, numelem * sizeof (MYREAL));
  memcpy (temp + 1 + numelem, lparam, numelem * sizeof (MYREAL));
  
  memset (nr->locilikes, 0, sizeof (MYREAL) * world->loci);
    
  //  if(world->loci==1)
  //  addon=1;
  //else
  //addon=0;
   
  for (worker = 1; worker < MIN (world->loci + addon, numcpu); worker++)
    {
      MYMPISEND (temp, (MYINT) numelem2+2, mpisizeof, (MYINT) worker, (MYINT) worker, comm_world);
    }
  for (worker = 1; worker < MIN (world->loci + addon, numcpu); worker++)
    {
      MYMPIRECV (tmp, (MYINT) world->loci, mpisizeof, MPI_ANY_SOURCE,
		 MPI_ANY_TAG, comm_world, &status);
      sender = status.MPI_SOURCE;
      tag = status.MPI_TAG;
      // the worker send a vector of values
      // e.g. (0,0,0,-12,0,-34,0,0,0), these are all loci
      // of which most of them were not evaluated
      // the loop updates the master copy of locilikes
      for (locus = 0; locus < world->loci; locus++)
        {
	  nr->locilikes[locus] += tmp[locus];
        }
    }
  for (locus = 0; locus < world->loci; locus++)
    {
      logres += nr->locilikes[locus];
    }
  nr->llike = logres;
  myfree(temp);
  myfree(tmp);
  return logres;
}


///
/// workers calculate likelihood and wait for master to receive work orders
void
mpi_likelihood_worker (world_fmt * world, helper_fmt * helper, long rep)
{
    long locus;
    long ww;
    
    nr_fmt *nr = helper->nr;
    
    MYREAL *param = helper->expxv;
    MYREAL *lparam = helper->xv;
    MYREAL *mu_rates = world->options->mu_rates;
    
    memset (nr->locilikes, 0, sizeof (MYREAL) * world->loci);
    
    if (world->options->gamma)
      {
        if (lparam[nr->numpop2] > LNMAXPARAM)
	  {
            lparam[nr->numpop2] = LNMAXPARAM;
	  }
        initgammacat (nr->categs, EXP (lparam[nr->numpop2]),1./* EXP (lparam[0])*/,
                      nr->rate, nr->probcat);
      }
    
    for (ww = 0; ww < locidone; ww++)
      {
        locus = nr->world->who[ww];
        if (!world->options->gamma)
	  {
            nr->locilikes[locus] =
	      calc_locus_like (nr, param, lparam, locus) + mu_rates[locus];
	  }
        else
	  {
            helper->locus = locus;
            nr->locilikes[locus] = gamma_locus_like (nr,param,lparam,helper->weight,locus);
	  }
      }
}

///
/// setup start parameters 
void
mpi_startparam_master(world_fmt * world)
{
    int tag;
    int numreceived = 0;

    long i;
    long sender;
    long workerloci = 0;

    MYREAL  *tmp;

    MPI_Status status;

    tmp = (MYREAL*) mycalloc(world->numpop2+1,sizeof(MYREAL));

    while (numreceived < world->loci)
      {
        MYMPIRECV (tmp, world->numpop2+1, mpisizeof, MPI_ANY_SOURCE,
                   MPI_ANY_TAG, comm_world, &status);
        sender = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        workerloci=tmp[0];
        for(i=0; i<world->numpop2; i++)
	  {
	    world->param0[i] += tmp[i+1];
	  }
        numreceived+=workerloci;
      }
    for(i=0; i<world->numpop2; i++)
      world->param0[i] /= world->loci;
    
    myfree(tmp);
}

///
/// set start parameters for workers
void
mpi_startparam_worker (world_fmt * world)
{
    long ww;
    long repstart;
    long repstop;
    long r;
    long i;
    long locus;

    MYREAL *tmp;

    if(locidone>0)
      {
	tmp = (MYREAL*) mycalloc(world->numpop2+1,sizeof(MYREAL));
        set_replicates (world, world->repkind, world->options->replicatenum,
                        &repstart, &repstop);
        tmp[0]=(MYREAL)locidone;
        for (ww = 0; ww < locidone; ww++)
	  {
            locus = world->who[ww];
            for (r = repstart; r < repstop; r++)
	      {
                for(i=0; i < world->numpop2; i++)
		  tmp[i+1] += world->atl[r][locus].param[i];
	      }
	  }
        for(i=1; i < world->numpop2+1; i++)
	  {
	    tmp[i] /= locidone * (repstop-repstart);
	  }
	MYMPISEND (tmp, world->numpop2+1, mpisizeof, MASTER, myID, comm_world);
	myfree(tmp);
      }
}


///
/// orchestrates max(gmax) over all nodes
void
mpi_gmax_master (world_fmt * world, long *Gmax)
{
    int tag;
    int sender;
    int numreceived = 0;

    long tmp;

    MPI_Status status;

    *Gmax = 0.;

    MYMPIBCAST (Gmax, ONE, MPI_LONG, MASTER, comm_world);

    while (numreceived < MIN(world->loci, numcpu - 1))
      {
	MYMPIRECV (&tmp, ONE, MPI_LONG, MPI_ANY_SOURCE,
		   MPI_ANY_TAG, comm_world, &status);
	sender = status.MPI_SOURCE;
	tag = status.MPI_TAG;
        if (*Gmax < tmp)
	  {
            *Gmax = tmp;
	  }
        numreceived++;
      }
    //  do we need this barrier really?
    MYMPIBARRIER(comm_world);
}

///
/// returns the gmax values to the master to evaluate max(gmax)
void
mpi_gmax_worker (world_fmt * world)
{
    long ww;
    long repstart;
    long repstop;
    long r;
    long locus;
    long Gmax = 1;

    MYMPIBCAST (&Gmax, ONE, MPI_LONG, MASTER, comm_world);
#ifdef DEBUG_MPI
    printf("%i> locidone=%i\n",myID,locidone);
#endif
    if(locidone>0)
      {
	set_replicates (world, world->repkind, world->options->replicatenum,
			&repstart, &repstop);
	
	for (ww = 0; ww < locidone; ww++)
	  {
	    locus = world->who[ww];
	    for (r = repstart; r < repstop; r++)
	      {
		if (Gmax < world->atl[r][locus].T)
		  Gmax = world->atl[r][locus].T;
	      }
	  }
	MYMPISEND (&Gmax, ONE, MPI_LONG, MASTER, myID, comm_world);
      }
    //  do we need this barrier really?
    MYMPIBARRIER(comm_world);
}

///
/// first worker (myID=1, myRepId=0) will send stop-message to replication nodes
/// the comm_worker group's master is the first worker who has ID=0 in this group
/// as results we send messages to id=1..x in the comm_worker group, do not mix this 
/// with the id in comm_world that include the master (id=0 there).
long  mpi_send_stop_mcmc_worker(long numcpu, long loci, MPI_Comm *comm, 
				     MPI_Request *irequests, MPI_Status *istatus, long id)
{
    long twolongs[2];
    long *temp;
    long receiver;
    long sent      = 0;
    long xx        = (id==0) ? 0 : 1;
 
    temp = (long *) mycalloc(TWO, sizeof(long));
    twolongs[0]=0;
    twolongs[1]=0;

    for(receiver=loci+1-xx; receiver< numcpu-xx; receiver++)
      {
	MYMPIISEND (temp, TWO, MPI_LONG, (MYINT) receiver, 0, *comm, &irequests[sent]);
        sent++;
      }
    if(sent>0)
      {
	MYMPIWAITALL(sent,irequests, istatus); // wait for all replicators to finish
      }
    myfree(temp);
    return sent;
}

///
/// first worker (myID=1, myRepId=0) will send stop-message to replication nodes
/// the comm_worker group's master is the first worker who has ID=0 in this group
/// as results we send messages to id=1..x in the comm_worker group, do not mix this 
/// with the id in comm_world that include the master (id=0 there).
long  mpi_send_stop_mcmc_lociworker(long numcpu, long loci)
{

  long receiver;
  long *temp;
  long xx       = (myID==0) ? 0 : 1;
  long sent     = 0;
  long minnodes = MIN(numcpu,loci -1);
  
  
  MPI_Request *irequests;
  MPI_Status *istatus;
  irequests = (MPI_Request *) mycalloc(minnodes+1,sizeof(MPI_Request));
  istatus = (MPI_Status *) mycalloc(minnodes+1,sizeof(MPI_Status));
  
  temp = (long *) mycalloc(TWO, sizeof(long));
  
  for(receiver=loci+1-xx; receiver< numcpu-xx; receiver++)
    {
           error("disable because testing why openmpi breaks");
      //      MYMPIISEND (temp, TWO, MPI_LONG, (MYINT) receiver, 0, comm_workers, &irequests[sent]);
      sent++;
    }
  if(sent>0)
    MYMPIWAITALL(sent,irequests, istatus); // wait for all replicators to finish
  myfree(temp);
  myfree(irequests);
  myfree(istatus);
  return sent;
}

///
/// stops all replicate workers
long  mpi_send_stop_mcmc_replicateworker(long numcpu, long loci)
{

  long *temp;
  long receiver;
  long sent      = 0;
  long xx        = (myID==0) ? 0 : 1;
  long minnodes  = labs(numcpu - loci -1);

  MPI_Request *irequests;
  MPI_Status *istatus;
  
  irequests = (MPI_Request *) mycalloc(minnodes+1,sizeof(MPI_Request));
  istatus   = (MPI_Status *) mycalloc(minnodes+1,sizeof(MPI_Status));
  temp      = (long *) mycalloc(TWO, sizeof(long));

  for(receiver=loci+1-xx; receiver< numcpu-xx; receiver++)
    {
      error("disable because testing why openmpi breaks");
      //     MYMPIISEND (temp, 2, MPI_LONG, (MYINT) receiver, 0, comm_workers, &irequests[sent]);
      sent++;
    }
  if(sent>0)
    {
      MYMPIWAITALL(sent,irequests, istatus); // wait for all replicators to finish
    }
  myfree(temp);
  myfree(irequests);
  myfree(istatus);
  return sent;
}

///
/// sends a stop signal to all loci-worker
void
mpi_send_stop (world_fmt * world)
{
  long worker;
  long numelem  = world->numpop2 + (world->options->gamma ? 1 : 0);
  long numelem2 = 2 * numelem;
  
  MYREAL *temp;

  temp = (MYREAL *) mycalloc (numelem2+2, sizeof (MYREAL));
  temp[0] = MIGMPI_END;

  for (worker = 1; worker < numcpu; worker++)
    {
      MYMPISEND (temp, (MYINT) numelem2+2, mpisizeof, (MYINT) worker, (MYINT) 0, comm_world); //end of loci
    }
  myfree(temp);
}

///
/// sends a stop signal to a specific worker used in for assignloci_worker
void mpi_send_stop_tag (int worker, world_fmt * world)
{
  long *temp;

  temp = (long *) mycalloc (TWO, sizeof (MYREAL));

  temp[0] = 0;
  temp[1] = 0;
  MYMPISEND (temp, (MYINT) TWO, MPI_LONG, (MYINT) worker, (MYINT) 0, comm_world);

  myfree(temp);

}

///
/// sends a continue with further analysis to the workers
void
mpi_results_stop (void)
{
    long worker;
    long dummy = 0;

    for (worker = 1; worker < numcpu; worker++)
      {
        MYMPISEND (&dummy, ONE, MPI_LONG, (MYINT) worker, (MYINT) 0, comm_world);
      }
}

///
/// sends out requests and receives all the locus-gradients from the workers
void
mpi_gradient_master (nr_fmt * nr, world_fmt * world, int *who)
{
    int tag;

    long locus;
    long sender;
    long *tempindex;
    long addon    = 1;//(world->loci == 1)? 0 : 1;
    long numelem  = nr->partsize;
    long numelem2 = 2 * numelem;

    MYREAL *temp;

    MPI_Status status;

    temp      = (MYREAL *) mycalloc (numelem2+2, sizeof (MYREAL));
    tempindex = (long *) mycalloc (numelem, sizeof (long));
    temp[0]   = MIGMPI_GRADIENT;

    memcpy (temp + 1, nr->param, numelem * sizeof (MYREAL));
    memcpy (tempindex, nr->indeks, nr->partsize * sizeof (long));
    memcpy (temp + 1 + numelem, nr->lparam, numelem * sizeof (MYREAL));
    temp[numelem2+1] = nr->profilenum;

    for (locus = 1; locus < MIN (world->loci + addon, numcpu); locus++)
      {
        MYMPISEND (temp, (MYINT) numelem2+2, mpisizeof, (MYINT) locus, (MYINT) locus, comm_world);
        MYMPISEND (tempindex, (MYINT) numelem, MPI_LONG, (MYINT) locus, (MYINT) locus, comm_world);
      }
    
    memset (nr->d, 0, sizeof (MYREAL) * numelem);
    for (locus = 1; locus < MIN (world->loci + addon, numcpu); locus++)
      {
	copy_and_clear_d (nr);
	MYMPIRECV (nr->d, (MYINT) numelem, mpisizeof, MPI_ANY_SOURCE,
		   MPI_ANY_TAG, comm_world, &status);
	add_back_d (nr);
	sender = status.MPI_SOURCE;
	tag = status.MPI_TAG;
      }
    myfree(temp);
    myfree(tempindex);
}

///
/// worker returns a gradient to master
void
mpi_gradient_worker (helper_fmt * helper, nr_fmt * nr,
                     timearchive_fmt ** tyme)
{
  long ww, locus;
  
  memset (nr->d, 0, sizeof (MYREAL) * nr->partsize);
  for (ww = 0; ww < locidone; ww++)
    {
      locus = nr->world->who[ww];
      if(!nr->world->options->gamma)
        {
	  copy_and_clear_d (nr);
	  simple_loci_derivatives (nr->d, nr, tyme, locus);
	  add_back_d (nr);
        }
      else
        {
	  gamma_locus_derivative (helper, locus);
        }
    }
}

///
/// provides calculations for the master. calculates likelihoods, gradients, but also
/// delivers results, migration histories, bayesian results and more.
/// workers stay a long time in this function and anwer requests from the master
void
mpi_maximize_worker (world_fmt * world, long kind, long rep)
{
  boolean done = FALSE;
  long locus;
  long repstart;
  long repstop;
  long Gmax;
  long numelem  =  world->numpop2 + (world->options->gamma ?  1 : 0) ;
  long numelem2 = numelem * 2;
  MYREAL *temp;
  MPI_Status status;
  nr_fmt *nr;
  helper_fmt helper;

  temp = (MYREAL *) mycalloc (numelem2 + 2, sizeof (MYREAL));
  helper.xv = (MYREAL *) mycalloc (numelem2, sizeof (MYREAL));
  helper.expxv = (MYREAL *) mycalloc (numelem2, sizeof (MYREAL));
  nr = (nr_fmt *) mycalloc (1, sizeof (nr_fmt));
  set_replicates (world, world->repkind, rep, &repstart, &repstop);
  which_calc_like (world->repkind);
  if(!world->options->bayes_infer)
    {
      MYMPIBCAST (&Gmax, 1, MPI_LONG, MASTER, comm_world);
      create_nr (nr, world, Gmax, 0, world->loci, world->repkind, repstart);
      SETUPPARAM0 (world, nr, world->repkind,
		   repstart, repstop, world->loci, kind, TRUE);
    }
  else
    {
      Gmax=1;
      create_nr (nr, world, Gmax, 0, world->loci, world->repkind, repstart);
    }
  while (!done)
    {
      MYMPIRECV (temp, numelem2+2, mpisizeof, MASTER, MPI_ANY_TAG,
		 comm_world, &status);
      locus = world->locus = status.MPI_TAG - 1;
      switch ((long) temp[0])
        {
        case MIGMPI_LIKE: // returns a log(likelihood) to the master
	  memset (nr->locilikes, 0, sizeof (MYREAL) * (world->loci+1));
	  memcpy (helper.expxv, temp + 1, sizeof (MYREAL) * numelem);
	  memcpy (helper.xv, temp + 1 + numelem,
		  sizeof (MYREAL) * numelem);
	  fill_helper (&helper, helper.expxv, helper.xv, world, nr);
	  mpi_likelihood_worker (world, &helper, rep);
	  MYMPISEND (nr->locilikes, (MYINT) world->loci, mpisizeof, (MYINT) MASTER,
		     (MYINT) (locus + 1), comm_world);
	  break;
        case MIGMPI_GRADIENT: // returns a gradient to the master
	  memcpy (nr->param, temp + 1, sizeof (MYREAL) * (numelem - 1));
	  memcpy (nr->lparam, temp + 1 + numelem,
		  sizeof (MYREAL) * numelem);
	  fill_helper (&helper, nr->param, nr->lparam, world, nr);
	  nr->profilenum = temp[numelem2 + 1];
	  MYMPIRECV (nr->indeks, (MYINT) numelem, MPI_LONG, MASTER, (MYINT) (locus+1), 
		     comm_world, &status);
	  mpi_gradient_worker (&helper, nr, world->atl);
	  MYMPISEND (nr->d, (MYINT) numelem, mpisizeof, (MYINT) MASTER, (MYINT) (locus + 1),
		     comm_world);
	  break;
        case MIGMPI_RESULT: // returns the results (stats for the loci worked on
	  mpi_results_worker ((long) temp[0], world, repstop, pack_result_buffer);
	  break;
	  //        case MIGMPI_PLOTPLANE: // returns the results (stats for the loci worked on
	  //mpi_results_worker ((long) temp[0], world, repstop, pack_plotplane_buffer);
	  //break;
        case MIGMPI_SUMFILE: // returns the tree summary statistic to the master
	  mpi_results_worker ((long) temp[0], world, repstop, pack_sumfile_buffer);
	  break;
        case MIGMPI_TREESPACE: // returns the treefile to the master
	  //	  printf("%i> print treespace to master\n",myID); 
	  mpi_results_worker ((long) temp[0], world, repstop, pack_treespace_buffer);
	  break;
        case MIGMPI_MIGHIST: // returns the histogram of events (top parts + obsolete parts)
	  mpi_results_worker ((long) temp[0], world, repstop, pack_mighist_buffer);
	  break;
        case MIGMPI_SKYLINE: // returns the histogram of events
	  //fprintf(stdout,"%i-------------------\n",myID);
	  //	  debug_skyline(world,"mpiresultsworker-skyline");
	  mpi_results_worker ((long) temp[0], world, repstop, pack_skyline_buffer);
	  //fprintf(stdout,"%i-------------------\n",myID);
	  break;
        case MIGMPI_BAYESHIST: // returns the histogram of the parameters
	  //	  printf("%i> send bayeshist\n",myID);
	  mpi_results_worker ((long) temp[0], world, repstop, pack_bayes_buffer);
	  break;
        case MIGMPI_END:
	  done = TRUE;
	  break;
        default:
	  fprintf (stdout, "%i> Do not understand task --> exit \n", myID);
	  MPI_Finalize();
	  exit (0);
        }
    }
  myfree(temp);
  myfree(helper.xv);
  myfree(helper.expxv);
  destroy_nr (nr, world);
}

///
/// sends out the options to all workers
void
broadcast_options_master (option_fmt * options, data_fmt *data)
{
  long allocbufsize = LONGLINESIZE;
  long bufsize;
  char *buffer;
  //  char *buffermu;
  buffer = (char *) mycalloc (allocbufsize, sizeof (char));
  bufsize = save_options_buffer (&buffer, &allocbufsize, options, data);
  MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
  MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
  myfree(buffer);
}

///
/// receives the options 
void
broadcast_options_worker (option_fmt * options)
{
  //  long i;
  //  char *temp;
  long bufsize;
  char *buffer, *sbuffer;
  MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
  buffer = (char *) mycalloc (bufsize + 1, sizeof (char));
  sbuffer = buffer;
  MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
  read_options_worker (&buffer, options);
  myfree(sbuffer);  
}


///
/// broadcasts the data to all nodes
void
broadcast_data_master (data_fmt * data, option_fmt * options)
{
  long bufsize;
  long allocbuffermusize = 1;
  char *buffer;
  char *buffermu;

  buffer = (char *) mycalloc (allocbuffermusize, sizeof (char));
  bufsize = pack_databuffer (&buffer, data, options);
  bufsize = (long) strlen(buffer)+1;
  MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
  MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
  myfree(buffer);
  if(options->murates)
    {
      buffermu = (char *) mycalloc (1, sizeof (char));
      bufsize = save_mu_rates_buffer(&buffermu,&allocbuffermusize, options);
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      MYMPIBCAST (buffermu, bufsize, MPI_CHAR, MASTER, comm_world);
      myfree(buffermu);
    }
}

///
/// receives the data from the master
void
broadcast_data_worker (data_fmt * data, option_fmt * options)
{
  long bufsize;
  char *buffer;
  char *buffermu, *sbuffermu;
  char *temp;
  long i;
  MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
  buffer = (char *) mycalloc (bufsize, sizeof (char));
  MYMPIBCAST (buffer, bufsize, MPI_CHAR, MASTER, comm_world);
  unpack_databuffer (buffer, data, options);
  myfree(buffer);
  if(options->murates)
    {
      MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
      buffermu = (char *) mycalloc (bufsize + 1, sizeof (char));
      sbuffermu = buffermu;
      MYMPIBCAST (buffermu, bufsize, MPI_CHAR, MASTER, comm_world);
      //fprintf(stdout,"%i> received muratebuffer:>%s<\n",myID, buffermu);
      temp=strsep(&buffermu,"\n");
      options->muloci = atoi(temp);
      if(options->mu_rates==NULL)
	options->mu_rates = (MYREAL *) mycalloc(1+options->muloci, sizeof(MYREAL));
      for(i=0; i< options->muloci; i++)
	{
	  temp = strsep(&buffermu,"\n");
	  options->mu_rates[i] = atof(temp);
	}
      myfree(sbuffermu);
    }
}

///
/// pack the data for shipment to the workers
long
pack_databuffer (char **buffer, data_fmt * data, option_fmt * options)
{
  long locus, pop, ind;
  long bufsize = 0;
  long biggest;
#ifdef UEP
  long sumtips, i;
#endif
  char fp[LONGLINESIZE];
  
  bufsize += LONGLINESIZE;
  *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
  bufsize += 1 + sprintf (fp, "%c %li %li %li\n", options->datatype, (long) data->hasghost,
			  data->numpop, data->loci);
  *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
  strcat (*buffer, fp);
  for (locus = 0; locus < data->loci; locus++)
    {
      bufsize += 1 + sprintf (fp, "%li\n", data->seq[0]->sites[locus]);
      *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
      strcat (*buffer, fp);
    }
  bufsize += 1 + sprintf (fp, "%li %f\n", data->seq[0]->addon, data->seq[0]->fracchange);
  *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
  strcat (*buffer, fp);
  // population data
  for (pop = 0; pop < data->numpop; pop++)
    {
      bufsize += 1 + sprintf (fp, "%s\n", data->popnames[pop]);
      *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
      strcat (*buffer, fp);
      biggest = 0;
      for(locus=0; locus<data->loci; locus++)
        {
	  bufsize += 1 + sprintf (fp, "%li %li\n", data->numind[pop][locus],data->numalleles[pop][locus]);
	  if(biggest < data->numind[pop][locus])
	    biggest = data->numind[pop][locus];
	  *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
	  strcat (*buffer, fp);
        }
      if (!strchr (SEQUENCETYPES, options->datatype))
        {
	  for (ind = 0; ind < biggest; ind++)
            {
	      bufsize += 1 + sprintf (fp, "%*.*s\n", (int) options->nmlength,
				      (int) options->nmlength, data->indnames[pop][ind][0]);
	      *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
	      strcat (*buffer, fp);
	      pack_allele_data (buffer, &bufsize, data, pop, ind);
            }
        }
      else
        {
	  for(locus=0;locus<data->loci; ++locus)
            {
	      for (ind = 0; ind < data->numind[pop][locus]; ind++)
                {
		  bufsize += 1 + sprintf (fp, "%*.*s\n", (int) options->nmlength, (int) options->nmlength,
					  data->indnames[pop][ind][locus]);
		  *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
		  strcat (*buffer, fp);
		  pack_sequence_data (buffer, &bufsize, data, pop, ind, locus);
                }
            }
        }
    }
  // geofile
  if (options->geo)
    {
      for (pop = 0; pop < data->numpop * data->numpop; pop++)
        {
	  bufsize += 1 + sprintf (fp, "%f %f\n", data->geo[pop], data->lgeo[pop]);
	  *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
	  strcat (*buffer, fp);
        }
    }
  // randomsubset
  for (pop = 0; pop < data->numpop; pop++)
    {
      for(locus=0;locus<data->loci;locus++)
	{
	      for(ind=0; ind < data->numind[pop][locus]; ind++)
		{
		  bufsize += 1 + sprintf (fp, "%li\n", data->shuffled[pop][locus][ind]);
		  *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
		  strcat (*buffer, fp);
		}
	}
    }

  // tipdate file
  if (options->has_datefile)
    {
      for(locus=0;locus<data->loci;locus++)
	{
	  for (pop = 0; pop < data->numpop; pop++)
	    {
	      for(ind=0; ind < data->numind[pop][locus]; ind++)
		{
		  bufsize += 1 + sprintf (fp, "%s %20.20f\n", data->sampledates[pop][locus][ind].name, data->sampledates[pop][locus][ind].date);
		  *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
		  strcat (*buffer, fp);
		}
	    }
	}
      bufsize += 1 + sprintf (fp, "%20.20f\n", data->maxsampledate);
      //      printf("%i> pack_data maxsampledate=%f %s\n",myID, data->maxsampledate, fp);
      *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
      strcat (*buffer, fp);
      //printf("%i> pack_data:%s\n""""""""""""""""""""""""\n\n\n",myID, *buffer);
    }

  // uepfile
#ifdef UEP
  if (options->uep)
    {
      sumtips = 0;
      for (pop = 0; pop < data->numpop; ++pop)
	sumtips += data->numind[pop][0];//Assumes UEP is matched by locus 1
      bufsize += 1 + sprintf (fp, "%li %li\n", sumtips, data->uepsites);
      *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
      strcat (*buffer, fp);
      if (strchr (SEQUENCETYPES, options->datatype))
        {
	  for (pop = 0; sumtips; pop++)
            {
	      for (i = 0; i < data->uepsites; i++)
                {
		  bufsize += 1 + sprintf (fp, "%i\n", data->uep[pop][i]);
		  *buffer =
		    (char *) myrealloc (*buffer, sizeof (char) * bufsize);
		  strcat (*buffer, fp);
                }
            }
        }
      else
        {
	  for (pop = 0; sumtips; pop++)
            {
	      for (i = 0; i < data->uepsites; i++)
                {
		  bufsize += 1 + sprintf (fp, "%i %i\n", data->uep[pop][i],
					  data->uep[pop + sumtips][i]);
		  *buffer =
		    (char *) myrealloc (*buffer, sizeof (char) * bufsize);
		  strcat (*buffer, fp);
                }
            }
        }
    }
  
#endif
  return bufsize;
}


void
pack_allele_data (char **buffer, long *bufsize, data_fmt * data, long pop,
                  long ind)
{
    char fp[LONGLINESIZE];
    long locus;
    for (locus = 0; locus < data->loci; locus++)
    {
        *bufsize += 1 + sprintf (fp, "%s %s\n", data->yy[pop][ind][locus][0],
                                 data->yy[pop][ind][locus][1]);
        *buffer = (char *) myrealloc (*buffer, sizeof (char) * *bufsize);
        strcat (*buffer, fp);
    }
}

void
pack_sequence_data (char **buffer, long *bufsize, data_fmt * data, long pop,
                    long ind, long locus)
{
    char *fp;
    //  long locus;
    //  fp = mycalloc (1, sizeof (char));
    // for (locus = 0; locus < data->loci; locus++)
    //   {
    fp = (char *) mycalloc ((2 + data->seq[0]->sites[locus]), sizeof (char));
    sprintf (fp, "%s\n", data->yy[pop][ind][locus][0]);
    *bufsize += 2 + data->seq[0]->sites[locus];
    *buffer = (char *) myrealloc (*buffer, sizeof (char) * *bufsize);
    strcat (*buffer, fp);
    //   }
    myfree(fp);
}

// this function and get_data() do not mix well!
void
unpack_databuffer (char *buffer, data_fmt * data, option_fmt * options)
{
    long locus, pop, ind, i=0;
    long biggest;
#ifdef UEP

    long sumtips;
#endif

    char *buf = buffer;
    char *input;
    char *name;
    long hasghost;
    input = (char *) mycalloc (LONGLINESIZE, sizeof (char));
    name = (char *) mycalloc (LONGLINESIZE, sizeof (char));
    sgets (input, LONGLINESIZE, &buf);
    sscanf (input, "%c%li%li%li", &options->datatype, &hasghost, &data->numpop,
            &data->loci);
    data->hasghost = (boolean) hasghost;
    init_data_structure1 (&data, options);
    for (locus = 0; locus < data->loci; locus++)
    {
        sgets (input, LONGLINESIZE, &buf);
        sscanf (input, "%li", &data->seq[0]->sites[locus]);
    }
    sgets (input, LONGLINESIZE, &buf);
#ifdef USE_MYREAL_FLOAT
    sscanf (input, "%li%f", &data->seq[0]->addon, &data->seq[0]->fracchange);
#else
    sscanf (input, "%li%lf", &data->seq[0]->addon, &data->seq[0]->fracchange);
#endif
    // population data
    for (pop = 0; pop < data->numpop; pop++)
    {
        sgets (input, LONGLINESIZE, &buf);
        sscanf (input, "%s", data->popnames[pop]);
        biggest=0;
        for(locus=0; locus<data->loci; locus++)
        {
            sgets (input, LONGLINESIZE, &buf);
            sscanf (input, "%li %li", &data->numind[pop][locus],&data->numalleles[pop][locus]);
            if(biggest<data->numind[pop][locus])
                biggest = data->numind[pop][locus];
        }
        init_data_structure2 (&data, options, pop);
        if (!strchr (SEQUENCETYPES, options->datatype))
        {
            for (ind = 0; ind < biggest; ind++)
            {
                sgets (input, LONGLINESIZE, &buf);
                sscanf (input, "%s", data->indnames[pop][ind][0]);
                for (locus = 0; locus < data->loci; locus++)
                {
                    sgets (input, LONGLINESIZE, &buf);
                    sscanf (input, "%s %s", data->yy[pop][ind][locus][0],
                            data->yy[pop][ind][locus][1]);
                }
            }
        }
        else
        {
            for (locus = 0; locus < data->loci; locus++)
            {
                for (ind = 0; ind < data->numind[pop][locus]; ind++)
                {
                    sgets (input, LONGLINESIZE, &buf);
                    strncpy(data->indnames[pop][ind][locus],input, options->nmlength);
                    input =(char *) myrealloc (input, sizeof (char) * (1 + data->seq[0]->sites[locus]));
                    sgets (input, 10 + data->seq[0]->sites[locus], &buf);
                    strncpy(data->yy[pop][ind][locus][0],input,1+data->seq[0]->sites[locus]);
                }
            }
        }
    }
    // geofile
    data->geo =
        (MYREAL *) mycalloc (1, sizeof (MYREAL) * data->numpop * data->numpop);
    data->lgeo =
        (MYREAL *) mycalloc (1, sizeof (MYREAL) * data->numpop * data->numpop);
    if (!options->geo)
    {
        for (i = 0; i < data->numpop * data->numpop; i++)
            data->geo[i] = 1.0;
    }
    else
    {
        for (pop = 0; pop < data->numpop * data->numpop; pop++)
        {
            sgets (input, LONGLINESIZE, &buf);
#ifdef USE_MYREAL_FLOAT
	    sscanf (input, "%f%f", &data->geo[pop], &data->lgeo[pop]);
#else
	    sscanf (input, "%lf%lf", &data->geo[pop], &data->lgeo[pop]);
#endif
        }
    }

  // randomsubset
    data->shuffled = (long ***) mycalloc(data->numpop,sizeof(long **));
  for (pop = 0; pop < data->numpop; pop++)
    {
      data->shuffled[pop] = (long **) mycalloc(data->loci,sizeof(long *));
      for(locus=0;locus<data->loci;locus++)
	{
	  data->shuffled[pop][locus] = (long *) mycalloc(data->numind[pop][locus],sizeof(long));
	  for(ind=0; ind < data->numind[pop][locus]; ind++)
	    {
	      sgets (input, LONGLINESIZE, &buf);
	      sscanf (input, "%li", &data->shuffled[pop][locus][ind]);
	    }
	}
    }
  

    data->sampledates = (tipdate_fmt ***) mycalloc (data->numpop, sizeof (tipdate_fmt **));
    data->sampledates[0] = (tipdate_fmt **) mycalloc (data->numpop * data->loci, sizeof (tipdate_fmt *));
    for (pop = 1; pop < data->numpop; pop++)
      {
	data->sampledates[pop] = data->sampledates[0] + data->loci * pop;
      }
    for(locus=0;locus<data->loci;locus++)
      {
	for(pop=0;pop < data->numpop; pop++)
	  {
	    data->sampledates[pop][locus] = (tipdate_fmt*) mycalloc(data->numind[pop][locus],sizeof(tipdate_fmt));
	  }
      }
    data->maxsampledate = 0.0;
    if (options->has_datefile)
      {
	for(locus=0;locus<data->loci;locus++)
	  {
	    for (pop = 0; pop < data->numpop; pop++)
	      {
		for(ind=0; ind < data->numind[pop][locus]; ind++)
		  {
		    sgets (input, LONGLINESIZE, &buf);
		    //printf("%i> %s\n",myID, input);
#ifdef USE_MYREAL_FLOAT
		    sscanf (input, "%s %f", name, &data->sampledates[pop][locus][ind].date);
#else
		    sscanf (input, "%s %lf", name, &data->sampledates[pop][locus][ind].date);

		    data->sampledates[pop][locus][ind].name = (char *) mycalloc(strlen(name)+1,sizeof(char));
		    strcpy(data->sampledates[pop][locus][ind].name, name);
		  }
	      }
	  }
	sgets (input, LONGLINESIZE, &buf);
#ifdef USE_MYREAL_FLOAT
	sscanf(input,"%f",&data->maxsampledate);
#else
	sscanf(input,"%lf",&data->maxsampledate);
#endif

	//	printf("%i> maxsampledate=%f %s\n''''''\n%s\n''''''''''''\n",myID, data->maxsampledate, input, buffer);
      }
#endif
    // uepfile
#ifdef UEP
    if (options->uep)
    {
        sgets (input, LONGLINESIZE, &buf);
        sscanf (input, "%li%li", &sumtips, &data->uepsites);
        data->uep =
            (int **) mycalloc (number_genomes (options->datatype) * sumtips,
                             sizeof (int *));
        if (strchr (SEQUENCETYPES, options->datatype))
        {
            for (pop = 0; sumtips; pop++)
            {
                data->uep[i] = (int *) mycalloc (data->uepsites, sizeof (int));
                for (i = 0; i < data->uepsites; i++)
                {
                    sgets (input, LONGLINESIZE, &buf);
                    sscanf (input, "%i", &data->uep[pop][i]);
                }
            }
        }
        else
        {
            for (pop = 0; sumtips; pop++)
            {
                data->uep[i] = (int *) mycalloc (data->uepsites, sizeof (int));
                data->uep[i + sumtips] =
                    (int *) mycalloc (data->uepsites, sizeof (int));
                for (i = 0; i < data->uepsites; i++)
                {
                    sgets (input, LONGLINESIZE, &buf);
                    sscanf (input, "%i%i", &data->uep[pop][i],
                            &data->uep[pop + sumtips][i]);
                }
            }
        }

    }
#endif
    init_data_structure3 (data);

    switch (options->datatype)
    {
    case 'a':
        create_alleles (data, options);
        break;
    case 'b':
        for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = XBROWN_SIZE;
        break;
    case 'm':
        create_alleles (data, options);
        for (pop = 0; pop < data->loci; pop++)
            data->maxalleles[pop] = options->micro_stepnum;
        break;
    }
    myfree(input);
    myfree(name);
}

#if 0
///
/// unpacks results from buffer to fill data structures in master for final printout
void
unpack_plotplane_buffer (MYREAL *buffer, world_fmt * world,
                      long locus, long maxrep, long numpop)
{

}
///
/// pack plotplane into buffer to ship to master
long
pack_plotplane_buffer (MYREAL **buffer, world_fmt * world,
                    long locus, long maxrep, long numpop)
{
}
#endif


///
/// unpacks results from buffer to fill data structures in master for final printout
void
unpack_result_buffer (MYREAL *buffer, world_fmt * world,
                      long locus, long maxrep, long numpop)
{
  long z=0;
  long rep;
  long pop;
  long addon=0;
  long numpop2 = world->numpop2;
  timearchive_fmt **atl = world->atl;
  MYREAL ***apg0 = world->apg0;

  if (maxrep > 1)
    addon = 1;

  for (rep = 0; rep < maxrep + addon; rep++)
    {
        atl[rep][locus].param_like = buffer[z++];
        for (pop = 0; pop < 4 * numpop2; pop++)
        {
            atl[rep][locus].parameters[pop] = buffer[z++];
        }
	//fprintf(stderr,"%i> unpacked result locus=%li replicate %li\n",myID,locus, rep);
    }
    // apg0
    for (rep = 0; rep < maxrep; rep++)
    {
        for (pop = 0; pop < world->options->lsteps; pop++)
        {
            apg0[rep][locus][pop] = buffer[z++];
        }
    }

    // BF material
    z = unpack_BF_buffer(buffer, z, locus, world);
    // ESS material
    z = unpack_ess_buffer(buffer, z, world);
}

///
/// pack results into buffer to ship to master
long
pack_result_buffer (MYREAL **buffer, world_fmt * world,
                    long locus, long maxrep, long numpop)
{
  long rep;
  long pop;
  long bufsize;
  long z = 0;
  long addon = 0;
  long numpop2 = world->numpop2;
  long nn = numpop2 + (world->bayes->mu * world->loci);
  timearchive_fmt **atl = world->atl;
  MYREAL ***apg0 = world->apg0;
  
  if (maxrep > 1)
    addon = 1;
  
  bufsize = (maxrep+addon) + (maxrep+addon) * 4 * numpop2 + maxrep * world->options->lsteps + \
    world->options->heated_chains * world->loci + 5 * world->loci + 1 + nn * 2;
  (*buffer) = (MYREAL *) myrealloc (*buffer, sizeof (MYREAL) * bufsize);
  memset (*buffer, 0, sizeof (MYREAL) * bufsize);
  
  for (rep = 0; rep < maxrep + addon; rep++)
    {
      (*buffer)[z++] = atl[rep][locus].param_like;
      for (pop = 0; pop < 4 * numpop2; pop++)
	{
	  (*buffer)[z++] =  atl[rep][locus].parameters[pop];
	}
      //      fprintf(stderr,"%i> packed result locus=%li replicate %li\n",myID,locus, rep);
    }
  // apg0
  for (rep = 0; rep < maxrep; rep++)
    {
      for (pop = 0; pop < world->options->lsteps; pop++)
        {
	  (*buffer)[z++] =  apg0[rep][locus][pop];
        }
    }
  // BF material
  if(!world->data->skiploci[locus])
    {
      //fprintf(stderr,"%i> packed result locus=%li replicate %li\n",myID,locus, rep);
      z = pack_BF_buffer(buffer, z, locus, world);
    }
  // ESS material
  z = pack_ess_buffer(buffer, z, world);

#ifdef DEBUG_MPI
  fprintf(stdout,"DEBUG: %i> z=%li, bufsize=%li\n", myID, z, bufsize);
#endif
  if(bufsize >= z)
    {
      bufsize = z;
    }
  else
    {
      fprintf(stderr,"%i> bufsize=%li < z=%li\n",myID,bufsize,z);
      error("pack_results tried to stuff to much into buffer in pack_result_buffer()\n");
    }
  return bufsize;
}




///
/// unpacks replicate samples of migration events, adds numbers and part-vectors
/// to the final array per locus, this function is only used with replicates over
/// multiple loci.
void
unpack_mighist_replicate_buffer_old (MYREAL *buffer, world_fmt * world,
                       long locus, long numpop)
{
  long i; 
  long j;
  long z = 0;
  long nummighist;
  long nummighistold;
  mighistloci_fmt *aa;
  aa = &world->mighistloci[locus];
  nummighist = (long) buffer[z++];
  nummighistold = aa->mighistnum;
  aa->mighistnum += nummighist;
  if(aa->allocsize <= aa->mighistnum)
    {
      aa->mighist = (mighist_fmt *) myrealloc (aa->mighist, sizeof (mighist_fmt) *(aa->mighistnum+1));
      for(j=aa->allocsize; j<=aa->mighistnum; j++)
        {
	  aa->mighist[j].allocsize=1;
	  aa->mighist[j].migeventsize=0;
	  aa->mighist[j].migevents =
	    (migevent_fmt *) mycalloc (1,  sizeof (migevent_fmt) );
	  //printf("%i> first events: alloc=%li size=%li\n", myID, aa->mighist[j].allocsize , aa->mighist[j].migeventsize);
        }
      aa->allocsize = aa->mighistnum+1;
    }
  for (j = nummighistold; j < aa->mighistnum; j++)
    {
      aa->mighist[j].copies = (long) buffer[z++];
      aa->mighist[j].weight = (long) buffer[z++];
      aa->mighist[j].migeventsize = (long) buffer[z++];
      //printf("%i> events: alloc=%li size=%li\n", myID, aa->mighist[j].allocsize , aa->mighist[j].migeventsize);
      aa->mighist[j].allocsize = aa->mighist[j].migeventsize + 1;
      aa->mighist[j].migevents = (migevent_fmt *) myrealloc (aa->mighist[j].migevents,
							     sizeof (migevent_fmt) *
							     aa->mighist[j].allocsize);
      for (i = 0; i < aa->mighist[j].migeventsize; i++)
        {
	  aa->mighist[j].migevents[i].age = buffer[z++];
	  aa->mighist[j].migevents[i].from = (long) buffer[z++];
	  aa->mighist[j].migevents[i].to = (long) buffer[z++];
	  aa->mighist[j].migevents[i].sumlines = (long) buffer[z++];
        }
    }
}
///
/// unpacks replicate samples of migration events, adds numbers and part-vectors
/// to the final array per locus, this function is only used with replicates over
/// multiple loci.
void
unpack_mighist_replicate_buffer(MYREAL *buffer, world_fmt * world,
                       long locus, long numpop)
{
  long i;
  long pop;
  long z         = 0;
  long          oldeventbinnum;
  long         * eventbinnum = NULL;
  duo         ** eventbins;
  mighistloci_fmt *aa;
  
  aa = &world->mighistloci[locus];
  eventbins = aa->migeventbins;
  eventbinnum = aa->migeventbinnum;
  

  for (pop = (world->options->mighist_all ? 0 : world->numpop); 
       pop <  world->numpop2; pop++)
    {
      oldeventbinnum = eventbinnum[pop];
      eventbinnum[pop] = (long) buffer[z++];
      if(oldeventbinnum < eventbinnum[pop])
	{
	  aa->migeventbins[pop] = (duo *) myrealloc(aa->migeventbins[pop], sizeof(duo) * eventbinnum[pop]);
	  memset(aa->migeventbins[pop][oldeventbinnum], 0, sizeof(duo) * (eventbinnum[pop] - oldeventbinnum));
	}
      for (i = 0; i < eventbinnum[pop]; i++)
	{
	  eventbins[pop][i][0] += buffer[z++];
	  eventbins[pop][i][1] += buffer[z++];
	}
    }
}


///
/// receive buffer with treespace content and unpack it.
/// Messing with format assignment to fit into the result_worker/master scheme
void
unpack_treespace_buffer (MYREAL *buffer, world_fmt * world,
                       long locus, long maxrep, long numpop)
{
  MYREAL like = -HUGE;
    char *input;
    char *sbuf;
    char *buf ;
    long size  = strlen((char *) buffer)+1;
    buf = NULL;
    //    printf("%i> in unpack_treespace_buffer() buffer has size %li\n%s\n\n",myID,size, (char*) buffer);
    if(size<=1)
      {
	return;
      } 
    //input = (char*) mycalloc(LINESIZE,sizeof(char));
    //    printf("%i> in unpack_treespace_buffer() after first if\n",myID);
    buf = (char *) mycalloc(size+1,sizeof(char));
    sbuf = buf;
    memcpy(buf, buffer,sizeof(char)*(size));
    //printf("%i> in unpack_treespace_buffer() after memcpy\n",myID);
    //    printf("%i> in unpack_treespace_buffer()  %s\n%s\n",myID,input,buf);
    //    printf("%i> Unpack++++++++++++\nsize=%li\n%s\n++++++++++++\n%s\n++++++++++++++++++++++++++++++++++++++++++++++++++\n",myID, size, buf, (char*) buffer);fflush(stdout);
    if(world->options->treeprint==BEST)
      {
	input = strsep(&buf,"@");
	//printf("%i> in unpack_treespace_buffer() after strsep %s\n",myID,input);fflush(stdout);
	like = atof(input);
	if(like >= world->besttreelike[locus])
	  {
	    world->besttreelike[locus] = like;
	    //	if(world->treespacenum[locus] < size)
	    //  {
	    //printf("%i> in unpack_treespace_buffer() before treespace realloc  %s\n%s\n",myID,input,buf);
	    size = strlen(buf) + 1;
	    world->treespacealloc[locus] = size;
	    world->treespacenum[locus] = size;
	    world->treespace[locus] = (char *) myrealloc(world->treespace[locus],
							 sizeof(char)*(size+1));
	    strcpy(world->treespace[locus],buf);
	  }
      }
    else
      {
	size = strlen(buf) + 1;
	world->treespacealloc[locus] = size;
	world->treespacenum[locus] = size;
	world->treespace[locus] = (char *) myrealloc(world->treespace[locus],
						     sizeof(char)*(size+1));
	strcpy(world->treespace[locus],buf);
      }
    myfree(sbuf);
}

///
/// pack treespace into buffer
/// ship to master; messing with format assignment because the standard transport buffer
/// is a double 
long pack_treespace_buffer (MYREAL **buffer, world_fmt * world,
			    long locus, long maxrep, long numpop)
{
  long thissize = strlen(world->treespace[locus])+1;
  long thisrealsize;
  char *input;
  char *ptr2, *ptr3;
  long pos = 0;
  MYREAL like = -HUGE;
  if(world->treespace[locus]==NULL)
    return 0;
  ptr2 = strstr(world->treespace[locus],"=");
  ptr3 = strstr(world->treespace[locus]," ]\n");
  pos  = ptr3-ptr2;
  
  input = mycalloc(LINESIZE,sizeof(char));
  if(world->options->treeprint == BEST)
    {
      strncpy(input, ptr2+1,pos);
      like = atof(input);
      pos   = sprintf(input,"%f @", like);
      //  fprintf(stdout,"%i> in pack_treespace_buffer() filling %li size\n",myID, thissize);
      // the master routine for this expects a MYREAL *buffer, but
      // the tree is a string
      thissize += pos;
    }
  thisrealsize = (long) (1. + (thissize * sizeof(char) / sizeof(MYREAL)));
  (*buffer) = (MYREAL *) myrealloc(*buffer, (1+thisrealsize) * sizeof(MYREAL));
  memset(*buffer, 0, sizeof(MYREAL) * (1+thisrealsize));
  // the sizeof(char) is NO MISTAKE!
  sprintf((char*)(*buffer), "%s%s", input, world->treespace[locus]);
  //  if(myID==1)
  //  {
  //    printf("%i> pack++++++++++++++++\nsize=%li, realsize=%li [char=%li real=%li]\n%s\n++++++++++++++++++++++++++++++++\n",myID, thissize,thisrealsize, sizeof(char), sizeof(MYREAL),(char *)(*buffer));
  //  }
  return thisrealsize;
}


void
unpack_mighist_buffer (MYREAL *buffer, world_fmt * world,
                       long locus, long maxrep, long numpop)
{
  long i;
  long pop;
  long z         = 0;
  long           oldeventbinnum;
  long         * eventbinnum = NULL;
  duo         ** eventbins;
  mighistloci_fmt *aa;
  
  aa = &world->mighistloci[locus];
  eventbins = aa->migeventbins;
  eventbinnum = aa->migeventbinnum;
  

  for (pop = (world->options->mighist_all ? 0 : world->numpop); 
       pop <  world->numpop2; pop++)
    {
      oldeventbinnum = eventbinnum[pop];
      eventbinnum[pop] = (long) buffer[z++];
      if(oldeventbinnum < eventbinnum[pop])
	{
	  aa->migeventbins[pop] = (duo *) myrealloc(aa->migeventbins[pop], sizeof(duo) * eventbinnum[pop]);
	  memset(aa->migeventbins[pop][oldeventbinnum], 0, sizeof(duo) * (eventbinnum[pop] - oldeventbinnum));
	}

      for (i = 0; i < eventbinnum[pop]; i++)
	{
	  eventbins[pop][i][0] = buffer[z++];
	  eventbins[pop][i][1] = buffer[z++];
	}
    }
}


///
/// pack migration events and coalescence events into  buffer to 
/// ship to master
long pack_mighist_buffer_old (MYREAL **buffer, world_fmt * world,
			  long locus, long maxrep, long numpop)
{
  long i;
  long j;
  long bufsize = 1;
  long z = 0;
  mighistloci_fmt *aa;
  long thin = 1;
  MYREAL *tmp;
  aa = &world->mighistloci[locus];
  for (j = 0; j < aa->mighistnum; j++)
    {
      //aa->mighist[j].migeventsize = 1;
      bufsize += aa->mighist[j].migeventsize;
    }
  // trial code to avoid a crash due to not enough memory available for buffer
  if((tmp = (MYREAL *) myrealloc ((*buffer), sizeof (MYREAL) * (bufsize+1))) == NULL)
    {
      do
	{
	  thin *= 10;
	  warning("Recovery mode for the migration events: thinning the events list by a factor of %li",thin);
	  bufsize = 0;
	  for (j = 0; j < aa->mighistnum; j+=thin)
	    {
	      if(j<aa->mighistnum)
		bufsize += aa->mighist[j].migeventsize;
	    }
	  tmp = (MYREAL *) myrealloc ((*buffer), sizeof (MYREAL) * (bufsize+1));
	} 
      while(tmp == NULL || thin < 10000);
    }
  (*buffer) = tmp;
  memset (*buffer, 0, sizeof (MYREAL) * (bufsize+1));
  
  (*buffer)[z++] = (MYREAL)(aa->mighistnum/thin);
  for (j = 0; j < aa->mighistnum; j+=thin)
    {
      if(j < aa->mighistnum)
	{
#ifdef DEBUG_MPI
	  printf ("%i> packmighistbuffer: %li %li %li\n", myID, aa->mighist[j].copies,
		  aa->mighist[j].weight, aa->mighist[j].migeventsize);
#endif
	  (*buffer)[z++] = (MYREAL) aa->mighist[j].copies;
	  (*buffer)[z++] = (MYREAL) aa->mighist[j].weight;
	  (*buffer)[z++] = (MYREAL) aa->mighist[j].migeventsize;
	  
	  for (i = 0; i < aa->mighist[j].migeventsize; i++)
	    {
	      (*buffer)[z++] = aa->mighist[j].migevents[i].age;
	      (*buffer)[z++] = (MYREAL) aa->mighist[j].migevents[i].from;
	      (*buffer)[z++] = (MYREAL) aa->mighist[j].migevents[i].to;
	      (*buffer)[z++] = (MYREAL) aa->mighist[j].migevents[i].sumlines;
	    }
	}
    }
  if(bufsize < z)
    {
      fprintf(stderr,"%i> bufsize=%li < z=%li\n",myID,bufsize,z);
      error("pack_mighist_buffer() failed\n");
    }
  return z;
}
///
/// pack migration events and coalescence events into  buffer to 
/// ship to master
long pack_mighist_buffer (MYREAL **buffer, world_fmt * world,
			  long locus, long maxrep, long numpop)
{
  long i;
  long pop;
  long z         = 0;
  long bufsize   = 0;
  long         * eventbinnum = NULL;
  duo         ** eventbins;
  mighistloci_fmt *aa;

  aa = &world->mighistloci[locus];
  eventbins = aa->migeventbins;
  eventbinnum = aa->migeventbinnum;
  for (pop = (world->options->mighist_all ? 0 : world->numpop); 
       pop <  world->numpop2; pop++)
    {
      bufsize += 1 + 2 * eventbinnum[pop];
    }
  (*buffer) = (MYREAL *) myrealloc ((*buffer), sizeof (MYREAL) * (bufsize+1));

  for (pop = (world->options->mighist_all ? 0 : world->numpop); 
       pop <  world->numpop2; pop++)
    {
      (*buffer)[z++] = (MYREAL)(eventbinnum[pop]);
      for (i = 0; i < eventbinnum[pop]; i++)
	{
	  (*buffer)[z++] = (MYREAL) eventbins[pop][i][0];
	  (*buffer)[z++] = (MYREAL) eventbins[pop][i][1];
	}
    }
  if(bufsize < z)
    error("pack_mighist_buffer() has a memory problem");
  return z;
}

///
/// receive buffer with skyline content and unpack it
void
unpack_skyline_buffer (MYREAL *buffer, world_fmt * world,
                       long locus, long maxrep, long numpop)
{
  //const   MYREAL invmaxrep = 1./((world->options->replicate
  //                       && world->options->replicatenum >
  //				  0) ? world->options->replicatenum : 1);
  const long np2 = numpop * numpop;
  long i;
  long j;
  long templ;
  long allocsize;
  long z   = 0;
  long *receive_eventbinnum;
  MYREAL temp;
  mighistloci_fmt *aa;

  receive_eventbinnum = (long *) mycalloc(np2, sizeof(long));

  aa = &world->mighistloci[locus];
  // read the number of bins for each parameter
  for(j=0; j< world->numpop2; j++)
    {
      templ = (long) buffer[z++];
      receive_eventbinnum[j] = templ;
      if(templ >= aa->eventbinnum[j])
	{
	  allocsize = templ + 1;
	  aa->eventbins[j] = (tetra *) myrealloc (aa->eventbins[j], sizeof (tetra) * allocsize);
	  //memset(aa->eventbins[j] + aa->eventbinnum[j],0,sizeof(tetra)*(templ-aa->eventbinnum[j]));
	  for(i=aa->eventbinnum[j];i < allocsize; i++)
	    {
	      aa->eventbins[j][i][0] = 0.F;
	      aa->eventbins[j][i][1] = 0.F;
	      aa->eventbins[j][i][2] = 0.F;
	      aa->eventbins[j][i][3] = 0.F;
	      aa->eventbins[j][i][4] = 0.F;
	      aa->eventbins[j][i][5] = 0.F;
	      aa->eventbins[j][i][6] = 0.F;
	      aa->eventbins[j][i][7] = 0.F;
	      aa->eventbins[j][i][8] = 0.F;
	    }
	  aa->eventbinnum[j] = allocsize;
	  //debug_skyline(world,"increased eventbins in unpack-skyline");
	}
    }
  // read time width of bins
  temp = (MYREAL) buffer[z++];
  if(temp - world->options->eventbinsize > EPSILON)
    error("problem with bins size transmission in unpack_skyline...");
  for(j=0; j< world->numpop2; j++)
    {
      for (i = 0; i < receive_eventbinnum[j]; i++)
	{
	  aa->eventbins[j][i][0] += (float) buffer[z++];// * invmaxrep;
	  aa->eventbins[j][i][1] += (float) buffer[z++];// * invmaxrep;
	  aa->eventbins[j][i][2] += (float) buffer[z++];// * invmaxrep;
	  aa->eventbins[j][i][3] += (float) buffer[z++];// * invmaxrep;
	  aa->eventbins[j][i][4] += (float) buffer[z++];// * invmaxrep;
	  aa->eventbins[j][i][5] += (float) buffer[z++];// * invmaxrep;
	}
    }
    //debug_skyline(world,"after unpack_skyline_buffer");
    myfree(receive_eventbinnum);
}

///
/// pack skyline into  buffer to 
/// ship to master
long pack_skyline_buffer (MYREAL **buffer, world_fmt * world,
                     long locus, long maxrep, long numpop)
{
  long j, i;
  long bufsize = world->numpop2 + 1;
  long z = 0L;
  mighistloci_fmt *aa;
  
  aa = &world->mighistloci[locus];
  //printf("Buffer in pack_skyline()=%f (%p)(%p)\n",(*buffer)[0], (*buffer), buffer);
  for (j = 0; j < world->numpop2; j++)
    {
      //      printf("%i> bufsize=%li eventbinnum=%li\n",myID, bufsize, aa->eventbinnum[j]);
      bufsize += 6 * aa->eventbinnum[j];
    }
  //  fprintf(stdout,"%i> bufsize in pack-skyline-buffer %li for locus %li\n", myID, bufsize, locus);
  //myfree(*buffer);
  (*buffer) = (MYREAL *) myrealloc ((*buffer), sizeof (MYREAL) * (bufsize));
  //(*buffer) = (MYREAL *) mycalloc(bufsize,sizeof (MYREAL));
  //  memset (*buffer, 0, sizeof (MYREAL) * (bufsize));
  // record how many bins in for each parameter
  for (j = 0; j < world->numpop2; j++)
    {
      (*buffer)[z++] = (MYREAL) aa->eventbinnum[j];
    }
  // time width of bins 
  (*buffer)[z++] = (MYREAL) world->options->eventbinsize;
  for (j = 0; j < world->numpop2; j++)
    {
      for (i = 0; i < aa->eventbinnum[j]; i++)
	{
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][0];
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][1];
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][2];
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][3];
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][4];
	  (*buffer)[z++] =  (MYREAL) aa->eventbins[j][i][5];
	}
    }
  if(bufsize < z)
    {
      fprintf(stderr,"error: bufsize=%li, z=%li\n", bufsize, z);
      error("pack_skyline_buffer: bufsize is too small");
    }
    return z;
}


///
/// unpack bayes parameters to fit into mpi_results_master()
void
unpack_bayes_buffer (MYREAL *buffer, world_fmt * world,
                       long locus, long maxrep, long numpop)
{
  // this combines the single_bayes and hist_bayes buffer unpacking
  unpack_hist_bayes_buffer(buffer, world->bayes, world, locus); 
}

///
/// pack bayes parameters to fit into mpi_results_worker()
long 
pack_bayes_buffer (MYREAL **buffer, world_fmt * world,
		   long locus, long maxrep, long numpop)
{
  long i;
  long bufsize;
  long z; 
  long sizec       = sizeof(MYREAL);
  //  long numbins     = 0;
  bayes_fmt *bayes = world->bayes;
  long np2         = world->numpop2;
  long npx         = np2 + (long) bayes->mu;
  long npp         = np2 + bayes->mu * world->loci;
  bayeshistogram_fmt  *hist;
  // memory for pack_single_bayes_buffer_part()
  // locus and numparams + accept and trial of genealogy
  bufsize =  4; 
  for(i=0; i < npx; i++)
    {
      if(bayes->map[i][1] != INVALID)
	  {
	    // per parameter
	    //for each parameter acceptance and trial +2 // for each parameter
	    bufsize += 2; 
	  }
    }
  
  // hist_bayes_buffer:
  // max buffer memory needed is (npp + 11*npp + (3 * npp * hist->bins[i])
  hist = &(world->bayes->histogram[locus]);
  for(i=0; i < npx; i++)
    {
      if(bayes->map[i][1] != INVALID &&	!world->options->has_bayesmdimfile)
	  {
	    //for each parameter
	    //                   bins number          +1
	    //                   datastore           +11
	    //                   bins*3              +3*bins
	    bufsize += 12; 
	    bufsize += 3 * hist->bins[i]; //set50, set95, result per bin 
	  }
    }
  // pack_BF_buffer
  //   hmscale, hm    +2
  //   bf: number of heated chains +world->options->heated_chains
  bufsize += 2 + world->options->heated_chains;
  // Autoarchive, ESS buffer:     parameters           +2*(numpop2+loci)
  //                              genealogy            +2
  bufsize += 2*npp + 2;
  //test printf("%i> bufsize=%li (npp=%li, numbins=%li)\n",myID, bufsize, npp, numbins);
  (*buffer) = (MYREAL *) myrealloc(*buffer, sizeof(MYREAL) * bufsize);
  memset (*buffer, 0,sizec * bufsize);
  z = pack_single_bayes_buffer_part(buffer,world->bayes,world,locus);
  //printf("%i> z=%li (%li)\n",myID,z, 2 + 2 * (np2+1));
  if(!world->options->has_bayesmdimfile)
    z = pack_hist_bayes_buffer(buffer, hist, world, z);
  // BF material
  if(!world->data->skiploci[locus])
    {
      // fprintf(stderr,"%i> pack_result_buffer(): packed BF result locus=%li replicate %li\n",myID,locus, -1);
      z = pack_BF_buffer(buffer, z, locus, world);
      if(z > bufsize)
	{
	  fprintf(stderr,"%i> ERROR: allocated bufsize=%li is smaller than used bufsize=%li\n",myID, bufsize, z);
	  error("buffer allocation overflowed");
	}
    }
  // ESS material
  z = pack_ess_buffer(buffer, z, world);
  if(z > bufsize)
    {
      fprintf(stderr,"%i> ERROR: allocated bufsize=%li is smaller than used bufsize=%li\n",myID, bufsize, z);
      error("buffer allocation overflowed");
    }
  return z;
}

///
/// unpack the bayes histogram
void unpack_hist_bayes_buffer(MYREAL *buffer, bayes_fmt *bayes, world_fmt *world, long locus)
{
    long                j, i;
    long                pa;
    long                z = 0;
    long                numbins = 0;
    long                pnum;
    long                tmp1, tmp2;
    long                tmplocus;
    bayeshistogram_fmt  *hist;
    long                total = 0;
    long                np2 = world->numpop2; 
    long                npp = np2 + bayes->mu; 
    long                npp11 = 11 * npp;
    //
    // begin unpack_single_bayes_buffer_part
    tmplocus = (long) buffer[z++]; 
    pnum = (long) buffer[z++];
#ifdef DEBUG_MPI
    fprintf (stdout, "%i> unpack_hist_bayes_buffer() received pnum=%li\n", myID, pnum);fflush(stdout);
#endif
    if(tmplocus!=locus)
    {
        bayes->numparams=0;
        locus = tmplocus;
    }
    for (j = 0; j < np2; ++j)
    {
      if(bayes->map[j][1] != INVALID)
	  {
	    bayes->accept[j] += (tmp1 = (long) buffer[z++]);
	    bayes->trials[j] += (tmp2 = (long) buffer[z++]);
#ifdef DEBUG_MPI
	    fprintf (stdout, "%i> received (acc %li) (trial %li) => (sumacc %li) (sumtrial %li)\n", myID, tmp1, tmp2, bayes->accept[j], bayes->trials[j]);
#endif        
	  }
    }
    // genealogy
    bayes->accept[j] += (tmp1 = (long) buffer[z++]);
    bayes->trials[j] += (tmp2 = (long) buffer[z++]);
    // end unpack_single_bayes_buffer
    //
    if(!world->options->has_bayesmdimfile)
      {
	bayes->histogram[locus].datastore = (MYREAL *) myrealloc(bayes->histogram[locus].datastore, sizeof(MYREAL) * (1+npp11));    
	hist = &(bayes->histogram[locus]);
	hist->minima = hist->datastore;    // contains minimal values for each parameter
        hist->maxima = hist->datastore + npp;    // contains maximal values for each parameter
        hist->adjmaxima = hist->datastore + 2*npp;// holds maxima values used in histogram [are smaller than maxima]
	hist->cred50l  = hist->datastore + 3*npp;    // holds 50%-credibility margins (<all lower values>, 
	hist->cred50u = hist->datastore + 4*npp;   //<all high values>)
	hist->cred95l = hist->datastore + 5*npp;    // holds 95%-credibility margins (<all lower values>)
	hist->cred95u = hist->datastore + 6*npp;   //<all high values>)
	hist->modes = hist->datastore + 7*npp;    // holds 95%-credibility margins (<all lower values>, <all high values>)
	hist->medians = hist->datastore + 8*npp;
	hist->means = hist->datastore + 9*npp;            
	hist->stds = hist->datastore + 10*npp;            


	for(i = 0; i < npp; ++i)
	  {
	    if(bayes->map[i][1] != INVALID)
	      {
		hist->bins[i] = (long) buffer[z++];
		total += hist->bins[i];
	      }
	  } 
	doublevec1d (&hist->results, total * npp);
	hist->set95 = (char *) mycalloc(total * npp * 2 + 2, sizeof(char));
	hist->set50 = hist->set95 + (total * npp + 1);
	
	for(i = 0; i < npp; ++i)
	  {
	    if(bayes->map[i][1] != INVALID)
	      {
		for(j=0;j<11;j++)
		  hist->datastore[11*i+j] = buffer[z++];
	      }
	  }
	numbins = 0;
	for(pa=0; pa < npp; pa++)
	  {
	    for(i=0;i<hist->bins[pa];i++)
	      {
		hist->set50[numbins + i] = (buffer[z++] < 1.0 ? '0' : '1'); 
		hist->set95[numbins + i] = (buffer[z++] < 1.0 ? '0' : '1'); 
		hist->results[numbins + i] = buffer[z++];
	      }
	    numbins += hist->bins[pa];
	    //
	    // CHECK
	    world->bayes->histtotal[locus * npp + pa] = hist->bins[pa];
	    //
	  }
      }
    // BF material
    z = unpack_BF_buffer(buffer, z, locus, world);
    // ESS material
    z = unpack_ess_buffer(buffer, z, world);
}
///
/// Bayes factor material buffer unpacker
long unpack_BF_buffer(MYREAL *buffer, long start, long locus, world_fmt * world)
{
  long i;
  long z = start;
  long hc = world->options->heated_chains;
  MYREAL              temp;
  MYREAL              *ttemp;
  MYREAL              *htemp;
  //  MYREAL              *atemp;
  
  ttemp = calloc(3 * world->loci + world->options->heated_chains + 2, sizeof(MYREAL));
  htemp = ttemp + world->loci;
  //atemp = htemp + world->loci;
  //
  // harmonic mean calculation
  temp = buffer[z++]; //hmscale
  if(temp < world->hmscale[locus])
    {      
      if(world->hm[locus]>0.0)
	world->hm[locus] *= EXP(temp - world->hmscale[locus]);
      world->hmscale[locus] = temp;
      //      printf("%i> locus=%li hmscale=%f hm=%f temp=%f\n",myID,locus, world->hmscale[locus], world->hm[locus], temp);
    }
  htemp[locus] = temp;//hmscale store
  world->hm[locus] += EXP(-htemp[locus]+world->hmscale[locus]) * buffer[z++];

  // thermodynamic integration
  for(i=0;i < world->options->heated_chains; i++)
    { 
      world->bf[locus * hc + i] +=  buffer[z++];
      //printf("%i> ****bf[%li + hc + %li]=%f\n",myID, locus, i, world->bf[locus*hc+i]);
    }

  myfree(ttemp);
  return z;
}

long unpack_ess_buffer(MYREAL *buffer, long start, world_fmt *world)
{
  long pa;
  long z = start;
  static long n=1;
  long npp = world->numpop2 + ((long) world->bayes->mu) * world->loci; 
  //printf("\n%i> npp =%li\n",myID,npp);
  //printf("%i> mu  =%i\n",myID,(int) world->bayes->mu);
  //printf("%i> loci=%li\n",myID,world->loci);
  //printf("%i> unpack_ess_buffer ***unpacking BF[0]=%f hm=%f npp=%li in unpack_hist_bayes_buffer()\n", myID, world->bf[0], world->hm[0], npp); 
  // unpacking autocorrelation and ess buffer
  for(pa=0; pa < npp; pa++)
    {
      if(world->bayes->map[pa][1] != INVALID)
	world->auto_archive[pa] += (buffer[z++] - world->auto_archive[pa])/n;
    }
  //genealogy
  world->auto_archive[pa] += (buffer[z++] - world->auto_archive[pa])/n;
  n++;
  //    printf("%i>>>>>> autoarchive %f\n",myID, world->auto_archive[0]);
  for(pa=0; pa < npp; pa++)
    {
      if(world->bayes->map[pa][1] != INVALID)
	world->ess_archive[pa] += buffer[z++];
    }
  //genealogy
  world->ess_archive[pa] += buffer[z++];
  return z;
}

///
/// Bayes factor material buffer packer
long pack_BF_buffer(MYREAL **buffer, long start, long locus, world_fmt * world)
{
  // buffer memory needs are (2 + #heatedchains)
  long i;
  long z  = start;
  const long hc = world->options->heated_chains;
  (*buffer)[z++] = world->hmscale[locus];
  (*buffer)[z++] = world->hm[locus];
  for(i=0; i < hc; i++)
    {
      (*buffer)[z++] = world->bf[locus * hc + i];
    }
  //  printf("%i> locus=%li send hmscale=%f hm=%f\n",myID, locus,world->hmscale[locus],world->hm[locus]);
  return z;
}

 /// packing autocorrelation and ess buffer
long pack_ess_buffer(MYREAL **buffer, long start, world_fmt *world)
{
  // buffer memory needs are (npp + 1) + (npp+1) (numpop2 + mu*loci)
  long i;
  long z = start;
  const long np2 = world->numpop2;
  const long npp = np2 + ((long) world->bayes->mu * world->loci); 
  //printf("\n%i> npp =%li\n",myID,npp);
  //printf("%i> mu  =%li\n",myID,(long)world->bayes->mu);
  //printf("%i> loci=%li\n",myID,world->loci);
  //printf("%i>> pack_essbuffer: npp=%li autoarchive %f\n",myID, npp, world->auto_archive[0]);
  for(i=0;i<np2;i++)
    {
      if(world->bayes->map[i][1] != INVALID)
	(*buffer)[z++] = world->auto_archive[i];
      //printf("%i> auto: %f\n",myID,world->auto_archive[i]);
    }
  //rates for each locus
  if(world->bayes->mu)
    {
      for(i=np2;i<npp;i++)
	(*buffer)[z++] = world->auto_archive[i];
    }
  // genealogy
  (*buffer)[z++] = world->auto_archive[i];
  //printf("%i> auto G: %f\n",myID,world->auto_archive[i]);
  for(i=0;i<np2;i++)
    {
      if(world->bayes->map[i][1] != INVALID)
	(*buffer)[z++] = world->ess_archive[i];
      //printf("%i> ess: %f\n",myID,world->ess_archive[i]);
    }
  if(world->bayes->mu)
    {
      for(i=np2;i<npp;i++)
	(*buffer)[z++] = world->auto_archive[i];
    }
  //genealogy
  (*buffer)[z++] = world->ess_archive[i];
  //printf("%i> ess G: %f\n",myID,world->ess_archive[i]);
  return z;
}

///
/// pack the bayes histogram
long pack_hist_bayes_buffer(MYREAL **buffer, bayeshistogram_fmt *hist, world_fmt * world, long startposition)
{
  // buffer memory needed is (npp + 11*npp + (3 * npp * hist->bins[i])
  long  j;
  long  i;
  long  npp     = world->numpop2 + world->bayes->mu;
  long  numbins = 0;
  long  z       = startposition;
  bayes_fmt *bayes = world->bayes;
#ifdef DEBUG_MPI
    printf("%i> pack_hist_bayes_buffer: position=%li last value = %f numparams=%li npp=%li\n",myID, startposition, (z > 0) ? (*buffer)[startposition] : -9999., hist->numparam, npp);
#endif
    
    for(i = 0; i < npp; ++i)
      {
	if(bayes->map[i][1] != INVALID)
          {
	    (*buffer)[z++] = (MYREAL) hist->bins[i];
	  }
      }
    //    printf("%i> npp=%li, z=%li\n",myID,npp,z);
    // parameter and mu datastore
    for(i = 0; i < npp; ++i)
      {
	if(bayes->map[i][1] != INVALID)
          {
	    for(j=0;j<11;j++)
	      (*buffer)[z++] = (MYREAL) hist->datastore[11*i+j];
	  }
      }
    numbins = 0;
    //printf("%i> npp2+=%li, z=%li\n",myID,npp+npp11,z);
    for(i=0; i < npp; i++)
      {
	// bins are zero for "c" and "0" parameters
	//	if(bayes->map[i][1] != INVALID)
        //  {
	    for(j=0;j<hist->bins[i];j++)
	      {
		(*buffer)[z++] = (MYREAL) (hist->set50[numbins + j]=='1' ? 1.0 : 0.0);
		(*buffer)[z++] = (MYREAL) (hist->set95[numbins + j]=='1' ? 1.0 : 0.0);
		(*buffer)[z++] = hist->results[numbins + j];
	      }
	    //  }
	numbins += hist->bins[i];
      }
#ifdef DEBUG_MPI
    printf("%i> pack_hist_bayes_buffer: position=%li numbins=%li, last value = %f\n",myID, z, numbins, (*buffer)[z-1]);
#endif
    return z;
}


///
/// unpack bayes parameter buffer, sent from replicant nodes AND lociworker to the master
/// the values will be simply added to the bayes->params, no records of replicates will be done.
long unpack_single_bayes_buffer(MYREAL *buffer,bayes_fmt * bayes, world_fmt * world,long locus)
{
  long i, j;
  long z = 0 ;
  long pnum;
  long tmp1, tmp2;
  long tmplocus;
  long allocparams = world->bayes->allocparams;
  //    long oldallocparams = world->bayes->allocparams;
  long repstart;
  long repstop;
  long nn = 2+world->numpop2 + (world->bayes->mu * world->loci) ;
  set_replicates (world, world->repkind, world->options->replicatenum,
		  &repstart, &repstop);
  tmplocus = (long) buffer[z++];
  pnum = (long) buffer[z++];
  //fprintf (stdout, "%i> received locus=%li (oldlocus=%li) pnum=%li\n", myID, tmplocus, locus, pnum);
  if(tmplocus!=locus)
    world->bayes->numparams=0;
  
  pnum += world->bayes->numparams;
  if(pnum >=world->bayes->allocparams)
    {
      allocparams = pnum + 1;
      world->bayes->params = (MYREAL *) myrealloc(world->bayes->params,sizeof(MYREAL)*allocparams*nn);
    }
  world->bayes->allocparams = allocparams;
  for(i = world->bayes->numparams; i < pnum; ++i)
    {
      // the first element is log(p(d|g)p(g|param))
      (world->bayes->params+(nn*i))[0] = buffer[z++];
      // the second element is log(p(d|g))
      (world->bayes->params+(nn*i))[1] = buffer[z++];
      //fprintf (stdout, "%i> receive params line %li ", myID, i);
      for (j = 2; j < world->numpop2+2; ++j) 
	  {
	    (world->bayes->params+(nn*i))[j] = buffer[z++];
	    //fprintf (stdout, "%f ", (world->bayes->params+(nn*i + 1))[j]);
	  }
      //fprintf (stdout, "\n");
    }
  world->bayes->numparams = pnum;
  // acceptance ratios are added to the ones we have already
  // parameter acceptances
  for (j = 0; j < world->numpop2; ++j)
    {
      if(bayes->map[j][1] != INVALID)
	{
	  tmp1 = (long) buffer[z++];
	  tmp2 = (long) buffer[z++];
	  world->bayes->accept[j] += tmp1;
	  world->bayes->trials[j] += tmp2;
	}
    }
  // the last acceptance is the one for the genealogies
  tmp1 = (long) buffer[z++];
  tmp2 = (long) buffer[z++];
  world->bayes->accept[j] += tmp1;
  world->bayes->trials[j] += tmp2;
  //printf("%i> unpack_single_bayes_buffer(): buffer counter is at %li",myID, z); 
  return z;
}

///
/// Pack bayes parameter buffer, sent from replicant nodes AND lociworker to the master
/// the values will be simply added to the bayes->params, no records of specific replicates are kept.
long pack_single_bayes_buffer(MYREAL **buffer, bayes_fmt *bayes, world_fmt *world,long locus)
{
    long i, j;
    long bufsize;
    long z = 0;
    //    long nump = world->numpop2 + 1;
    long nn = 2 + world->numpop2 + (world->bayes->mu * world->loci) ;
    bufsize = 2 * (world->numpop2 + 1); //acceptance ratio: params + tree
    bufsize += 2; // loci + numparams
    bufsize += world->bayes->numparams * nn;
    //printf("%i> bufsize in pack_single_bayes_buffer()=%li\n",myID,bufsize);
    (*buffer) = (MYREAL *) myrealloc(*buffer, bufsize * sizeof(MYREAL));
    memset (*buffer, 0, bufsize * sizeof(MYREAL));
#ifdef DEBUG_MPI
    fprintf(stdout, "%i>>>>>\n  buffersize=%li, numparams=%li\n>>>>\n", 
    	    myID,
	    bufsize,
	    world->bayes->numparams
	    );
    fflush(stdout);
#endif

    (*buffer)[z++] = (MYREAL) locus;
    (*buffer)[z++] = (MYREAL) bayes->numparams;

    for(i = 0; i < world->bayes->numparams; ++i)
      {
	//the first and second elements are logprob                                                                                   
	(*buffer)[z++] = (bayes->params+(i*nn))[0];
	(*buffer)[z++] = (bayes->params+(i*nn))[1];
	for (j = 0; j < world->numpop2; ++j) //the first and second elements are logprob                                              
	  {
	    if(bayes->map[j][1] != INVALID)
	      (*buffer)[z++] = (bayes->params+(i*nn))[j+2];
	  }
      }
    // for the parameters                                                                                                           
    for (j = 0; j < world->numpop2; ++j)
      {
	if(bayes->map[j][1] != INVALID)
	  {
	    (*buffer)[z++] =  (MYREAL) bayes->accept[j];
	    (*buffer)[z++] =  (MYREAL) bayes->trials[j];
	  }
      }
    // for the genealogy                                                                                                            
    (*buffer)[z++] =  (MYREAL) bayes->accept[j];
    (*buffer)[z++] =  (MYREAL) bayes->trials[j];

    if(bufsize < z)
      error("buffer is too small in pack_single_bayes_buffer()\n");
    return z;
}
///
/// Pack bayes parameter buffer, sent from replicant nodes AND lociworker to the master
/// the acceptance and trial values will be simply added.

long pack_single_bayes_buffer_part(MYREAL **buffer, bayes_fmt *bayes, world_fmt *world,long locus)
{
    long j;
    long z = 0;
    (*buffer)[z++] = (MYREAL) locus;
    (*buffer)[z++] = (MYREAL) bayes->numparams;

    // parameters
    for (j = 0; j < world->numpop2; ++j)
    {
      if(bayes->map[j][1] != INVALID)
	{
	  (*buffer)[z++] = (MYREAL) bayes->accept[j]; 
	  (*buffer)[z++] = (MYREAL)  bayes->trials[j];
	}
    }
    // genealogy
    (*buffer)[z++] = (MYREAL) bayes->accept[j]; 
    (*buffer)[z++] = (MYREAL)  bayes->trials[j];
    return z;
}

///
/// unpack minimal statistic trees
/// \todo there are differences between unpack_sumfile() and and read_savesum() this needs reconciliation
void
unpack_sumfile_buffer (MYREAL *buffer, world_fmt * world,
                       long locus, long maxrep, long numpop)
{
    long replicate;
    timearchive_fmt **ta = world->atl;
    long z=0;
    for (replicate = 0; replicate < maxrep; replicate++)
    {
      unpack_single_sumfile_buffer (buffer, ta, world, locus, replicate, numpop,&z);
    }
}



///
/// unpack minimal statistic trees for a single replicate
void
unpack_single_sumfile_buffer (MYREAL *buffer, timearchive_fmt **ta, world_fmt *world,
                              long locus, long replicate, long numpop, long *startz)
{
    long i, j;
    long z = *startz;
    long numpop2plus = numpop * numpop + 2 * numpop;
    ta[replicate][locus].T = (long) buffer[z++];
    ta[replicate][locus].numpop = (long) buffer[z++];
    ta[replicate][locus].sumtips = (long) buffer[z++];
    ta[replicate][locus].param_like = buffer[z++];
    //fprintf(stderr,"%i> unpack_single_sumfile: start_z=%li, replicate=%li, locus=%li ", myID, *startz, replicate, locus);
    world->chainlikes[locus][replicate] = ta[replicate][locus].param_like;
    //fprintf(stdout,"%i> got sumfile locus %li and replicate %li (%li)\n",myID,locus,replicate, ta[replicate][locus].allocT);
    increase_timearchive (world, locus, ta[replicate][locus].T, world->numpop, replicate);
    for (i = 0; i < ta[replicate][locus].T; i++)
    {
      ta[replicate][locus].tl[i].copies = buffer[z++];
      ta[replicate][locus].tl[i].lcopies =  buffer[z++];
      for (j = 0; j < numpop2plus; j++)
        {
	  ta[replicate][locus].tl[i].data[j] = buffer[z++];
        }
    }
    for (i = 0; i < world->numpop2; i++)
    {
      ta[replicate][locus].param[i] = buffer[z++];
      ta[replicate][locus].param0[i] = buffer[z++]; 
    }
    log_param0 (ta[replicate][locus].param0, ta[replicate][locus].lparam0, world->numpop2);
    ta[replicate][locus].trials = buffer[z++];
    ta[replicate][locus].normd = buffer[z++];
    *startz = z;
    //fprintf(stderr,"<%i> end_z=%li\n", myID, *startz);
}


long pack_single_sumfile_buffer(MYREAL **buffer, long z, world_fmt * world,
                                long locus, long replicate, long numpop)
{
    long i, j; 
    long numpop2 = numpop * numpop;
    long numpop2plus = numpop2 + 2 * numpop;
    timearchive_fmt **ta = world->atl;

    (*buffer)[z++] = (MYREAL) ta[replicate][locus].T;
    (*buffer)[z++] = (MYREAL) ta[replicate][locus].numpop;
    (*buffer)[z++] = (MYREAL) ta[replicate][locus].sumtips;
    (*buffer)[z++] = ta[replicate][locus].param_like;
    
    for (i = 0; i < ta[replicate][locus].T; i++)
      {
	(*buffer)[z++] = (MYREAL) ta[replicate][locus].tl[i].copies;
	(*buffer)[z++] = ta[replicate][locus].tl[i].lcopies;
	for (j = 0; j < numpop2plus; j++)
	  {
	    (*buffer)[z++] = ta[replicate][locus].tl[i].data[j];
	  }
    }
    for (i = 0; i < numpop2; i++)
      {
	(*buffer)[z++] = ta[replicate][locus].param[i];
	(*buffer)[z++] = ta[replicate][locus].param0[i];
      }
    (*buffer)[z++] = (MYREAL) ta[replicate][locus].trials;
    (*buffer)[z++] = ta[replicate][locus].normd;
    return z;    
}

long
pack_sumfile_buffer (MYREAL **buffer, world_fmt * world,
                     long locus, long maxrep, long numpop)
{
  long replicate;
  long bufsize = 1;
  long allocbufsize = 0;
  long z = 0;
  long numpop2 = numpop * numpop;
  long numpop2plus = numpop2 + 2 * numpop;
  timearchive_fmt **ta = world->atl;
  
  for (replicate = 0; replicate < maxrep; replicate++)
    {
      bufsize = ta[replicate][locus].T; 
      bufsize *= numpop2plus;
      bufsize += 2 * ta[replicate][locus].T; 
      bufsize += 2 * numpop2;
      bufsize += 6;
      allocbufsize += bufsize;
      //fprintf(stderr,"%i> bufsize=%li allocbufsize=%li\n", myID, bufsize, allocbufsize);
    }
  (*buffer) = (MYREAL *) myrealloc ((*buffer), allocbufsize * sizeof (MYREAL));
  memset(*buffer, 0, sizeof(MYREAL) * allocbufsize);
  for (replicate = 0; replicate < maxrep; replicate++)
    {
      z = pack_single_sumfile_buffer(buffer, z, world, locus, replicate, numpop);
      //fprintf(stderr,"%i> PACKED SUMFILE FOR REPLICATE %li z=%li locus=%li allocbufsize=%li\n", myID, replicate, z, locus, allocbufsize);
    }
    bufsize = z;
    if(bufsize > allocbufsize)
      {
	fprintf(stderr,"bufsize=%li allocbufsize=%li\n", bufsize, allocbufsize);
        error("allocation exceeded in pack_sumfile_buffer");
      }
    return bufsize;
}


///
/// gather results (sumfiles, results, migrate-histogram, ..) from workers
void
mpi_results_master (long sendtype, world_fmt * world, long maxreplicate,
                    void (*unpack) (MYREAL *buffer, world_fmt * world,
                                    long locus, long maxrep, long numpop))
{
#ifdef DEBUG_MPI
  long ii;
#endif
    long numpop = world->numpop;
    long bufsize = 1;
    // maxreplicate > 1 ---> add 1 [this all needs careful checking]
    // MIGMPI_SUMFILE -----> 0 
    // MIGMPI_HIST    -----> 0
    // still strange? long addon = (maxreplicate>1) ? 1 : ((sendtype == MIGMPI_SUMFILE) ||  (sendtype == MIGMPI_MIGHIST) )? 0 : ((world->loci == 1) ? 0 : 1) ;
    //long addon = (maxreplicate>1) ? 1 : ((world->loci > 1) ? 1 : 0) ;
    long addon = 1;
    //    boolean done = FALSE;
    MYREAL *buffer;
    MYREAL *temp;
    int worker;
    long z, tag, sender;
    MPI_Status status;
    long numelem = world->numpop2 + (world->options->gamma ? 1 : 0);
    long numelem2 = 2 * numelem;

    temp = (MYREAL *) mycalloc (numelem2 + 2, sizeof (MYREAL));
    buffer = (MYREAL *) mycalloc (bufsize+1, sizeof (MYREAL));
    temp[0] = (MYREAL)sendtype;
    temp[1] = (MYREAL) bufsize;
    for (worker = 1; worker < MIN (world->loci + addon, numcpu); worker++)
    {
      //printf("%i> MASTER requests information from n%i for locus %i\n",myID, worker, worker-1);
        MYMPISEND (temp, numelem2 + 2, mpisizeof, worker, worker, comm_world);
    }
    z = 0;
    while (z < world->loci)
    {
        MYMPIRECV (&bufsize, ONE, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
                  comm_world, &status);
        buffer = (MYREAL *) myrealloc (buffer, sizeof (MYREAL) * (bufsize + 1));
        memset (buffer, 0, sizeof (MYREAL) * (bufsize + 1));
        sender = status.MPI_SOURCE;
        tag = status.MPI_TAG;
#ifdef DEBUG_MPI
	fprintf(stdout, "%i> z=%li worker=%li bufsize=%li -------------------------------------\n",myID, z, sender, bufsize);
#endif
        MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
#ifdef DEBUG_MPI
	fprintf(stdout,"%i>------------------------------------\nbuffer=",myID);
	for(ii=0;ii<bufsize;ii++)
	  fprintf(stdout," %f", buffer[ii]);
	fprintf(stdout,"\n%i>-------------------------------------\n",myID);
#endif
        (*unpack) (buffer, world, tag - 1, maxreplicate, numpop);
	//fprintf(stdout,"%i> unpacked bufsize=%li from node %li\n",myID,bufsize,sender); fflush(stdout);
        z++;
    }
    myfree(buffer);
    myfree(temp);
}

void
mpi_results_worker (long bufs, world_fmt * world, long maxrep,
                    long (*pack) (MYREAL **buffer, world_fmt * world,
                                  long locus, long maxrep, long numpop))
{
    long numpop = world->numpop;
    long ww, locus;
    MYREAL *allbuffer;
    long bufsize = 1;
    allbuffer = (MYREAL *) mycalloc(1, sizeof(MYREAL));
    //fprintf(stdout,"%i> locidone=%i\n",myID, locidone); fflush(stdout);
    for (ww = 0; ww < locidone; ww++)
    {
      locus = world->who[ww];
      bufsize = (*pack) (&allbuffer, world, locus, maxrep, numpop);
#ifdef DEBUG_MPI
      fprintf(stdout,"%i> locus=%li after pack bufsize=%li\n",myID, locus, bufsize); fflush(stdout);
#endif
      MYMPISEND (&bufsize, ONE, MPI_LONG, MASTER, (int) (locus + 1), comm_world);
      //fprintf(stdout,"%i> sending results from locus=%li using bufsize=%li to master \n",myID, locus, bufsize); fflush(stdout);
      MYMPISEND (allbuffer, bufsize, mpisizeof, MASTER, (int) (locus + 1), comm_world);
    }
    myfree(allbuffer);
}

void
mpi_broadcast_results (world_fmt * world, long loci,
                       long (*pack) (MYREAL **buffer, world_fmt * world,
                                     long locus, long maxrep, long numpop),
                       void (*unpack) (MYREAL *buffer, world_fmt * world,
                                       long locus, long maxrep, long numpop))
{
    long locus;
    // long addon = (world->loci == 1) 0 : 1;
    long bufsize=1;
#ifdef DEBUG_MPI
    char nowstr[STRSIZE];
#endif
    MYREAL *allbuffer = NULL;// = &world->buffer;

    long maxreplicate = (world->options->replicate
                         && world->options->replicatenum >
                         0) ? world->options->replicatenum : 1;
    //    allbuffer = (char *) mycalloc (10, sizeof (char));
#ifdef DEBUG_MPI
    get_time (nowstr, "%H:%M:%S");
    if(world->options->progress)
      fprintf(stdout, "%i> Redistributing the data\nResult parts [Time is %s]\n",myID, nowstr);
#endif
    for (locus = 0; locus < loci; locus++)
    {
        if (myID == MASTER)
        {
            bufsize =(*pack) (&allbuffer, world, locus, maxreplicate,
                              world->numpop);
            MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
            MYMPIBCAST (allbuffer, bufsize, mpisizeof, MASTER, comm_world);
#ifdef DEBUG_MPI
            printf("%i> Locus %li results sent\n",myID, locus);
#endif
        }
        else
        {
            MYMPIBCAST (&bufsize, 1, MPI_LONG, MASTER, comm_world);
            allbuffer = (MYREAL *) myrealloc (allbuffer, sizeof (MYREAL) * bufsize + 1);
            MYMPIBCAST (allbuffer, bufsize, mpisizeof, MASTER, comm_world);
            (*unpack)(allbuffer, world, locus, maxreplicate,
                      world->numpop);
#ifdef DEBUG_MPI
            printf("%i> Locus %li results received\n",myID, locus);
#endif
        }
	myfree(allbuffer);
	allbuffer=NULL;
	//        memset (allbuffer, 0, sizeof (char) * bufsize);
    }
    MYMPIBARRIER(comm_world);
    myfree(allbuffer);
}

// slownet profiler
//
#ifdef SLOWNET

///
/// orchestrates and prints the profile likelihoods using MPI
/// when the slownet option for MPI work is chosen
void
mpi_profiles_master (world_fmt * world, long nparam, int *profilewho)
{

  boolean done;
  int     tag;
  int     sender   = 0;

#ifdef PRETTY
  long    location;
#endif
  long    pnum;
  long    pdone;
  long    bufsize;
  long    quantsize;
  long    failuresize;
  long    tempstrsize=STRSIZE;
  long    temp[4];
  long    minnodes = MIN (nparam, numcpu - 1);
  long    numsent  = 0;

  char   *tempstr;
  char   *leadstr;
  char  **buffer = &world->buffer;

  FILE   *outfile = world->outfile;

  MPI_Request *irequests;
  MPI_Status  status;
  MPI_Status  *istatus;

  irequests = (MPI_Request *) mycalloc(minnodes,sizeof(MPI_Request));
  istatus = (MPI_Status *) mycalloc(minnodes,sizeof(MPI_Status));
  
  tempstr = (char*) mycalloc(tempstrsize,sizeof(char));
  leadstr = (char*) mycalloc(STRSIZE,sizeof(char));
  
  for (pnum = 0; pnum < minnodes; pnum++)
    {
      MYMPIISEND (&pnum, 1, MPI_LONG, (MYINT) (pnum + 1), (MYINT) (pnum + 1), comm_world,&irequests[numsent]);
      //fprintf(stdout,"%i>>>>> sent parameter number %li to node %li with tag %li\n",myID,pnum,pnum+1,pnum+1);
        numsent++;
    }
    MYMPIWAITALL(minnodes,irequests, istatus);

    for (pnum = 0; pnum < nparam; pnum++)
    {
        done = FALSE;
        while(!done)
        {
            MYMPIRECV (leadstr, STRSIZE, MPI_CHAR, (MYINT) MPI_ANY_SOURCE, (MYINT) MPI_ANY_TAG, comm_world, &status);
            sender = status.MPI_SOURCE;
            tag = status.MPI_TAG;
            switch(leadstr[0])
            {
            case 'M':
	      tempstrsize = atol(leadstr+1);
	      tempstr = (char*) realloc(tempstr,sizeof(char)*(tempstrsize+1));
	      memset(tempstr,0,sizeof(char)*(tempstrsize+1));
	      //	      fprintf(stdout,"%i> ready to receive %li chars\n",myID,tempstrsize);
	      MYMPIRECV (tempstr, tempstrsize, MPI_CHAR, (MYINT) sender, (MYINT) tag,
			 comm_world, &status);                
	      handle_message(tempstr,sender, world);
	      break;
            case 'P':
	      //fprintf(stdout,"%i> ready to receive 4 longs\n",myID);
                MYMPIRECV (temp, FOUR, MPI_LONG, (MYINT) sender, (MYINT) MPI_ANY_TAG,
                           comm_world, &status);
                pdone = temp[0];
                bufsize = temp[1];
                quantsize = temp[2];
                failuresize = temp[3];
                //fprintf(stdout,"%i> ++++++++++++ bufsize=%li quantsize=%li failuresize=%li from sender %i\n",
                //        myID,bufsize,quantsize,failuresize,sender);
                profilewho[pdone] = sender;
                *buffer =
                    (char *) myrealloc (*buffer, sizeof (char) * (bufsize + quantsize + failuresize + 1));
                //memset (*buffer, 0, sizeof (char) * (bufsize + quantsize + failuresize + 1));
                MYMPIRECV (*buffer, bufsize + quantsize + failuresize, MPI_CHAR, (MYINT) sender, (MYINT) tag,
                          comm_world, &status);
                //fprintf(stdout,"@%s@\n\n@%s@\n\n\n",*buffer+bufsize,*buffer+bufsize+quantsize);
                //fprintf(stdout,"########################\n# %i> buf %li from %i \n########################\n",
                //myID, (long) (long) strlen(*buffer)+1, sender);
                if(world->options->printprofsummary)
                {
                    unpack_quantile ((*buffer) + bufsize, world->quantiles[pdone],
                                     GRIDSIZE);
                    unpack_failed_percentiles ((*buffer) + bufsize + quantsize, world->percentile_failed[pdone],
                                     GRIDSIZE);
//                    fprintf(stdout,"%i>failed=%i %i %i %i %i %i %i\n",myID,world->percentile_failed[pdone][0]
  //                          ,world->percentile_failed[pdone][1]
    //                        ,world->percentile_failed[pdone][2]
      //                      ,world->percentile_failed[pdone][3]
        //                    ,world->percentile_failed[pdone][4]
          //                  ,world->percentile_failed[pdone][5]
            //                ,world->percentile_failed[pdone][6]);
                    memset ((*buffer) + bufsize , 0, sizeof (char) * (quantsize + failuresize));
                }
                // print profile table that is in the first part of the buffer
                fprintf (outfile, "%s\n\n", *buffer);
#ifdef PRETTY
		location = strlen(*buffer);
		pdf_print_profile_table(55.0f, &page_height, world->options->profilemethod, *buffer + location + 1, world);
#endif
                done=TRUE;
                break;
            default:
	      fprintf(stderr,"%i> message=%s\n%i> sender=%i tag=%i\n",myID,tempstr, myID,status.MPI_SOURCE,status.MPI_TAG);
                MPI_Finalize();
                error("DIED because of wrong message from worker");
                break;
            }
        }
        if (numsent < nparam)
        {
            MYMPISEND (&numsent, ONE, MPI_LONG, (MYINT) sender, (MYINT) (numsent + 1), comm_world);
            numsent++;
        }
        else
        {
            // stop worker because there is nothing to do anymore
            MYMPISEND (&nparam, ONE, MPI_LONG, sender, 0, comm_world); //end of parameter list
        }
    }
    // stop workers that did nothing for profiles
    for (sender = MIN (nparam, numcpu - 1) + 1; sender < numcpu ; sender++)
    {
        // stop all nodes to wait for profiles
        MYMPISEND (&nparam, ONE, MPI_LONG, sender, 0, comm_world); 
    }
    myfree(tempstr);
    myfree(leadstr);
    myfree(istatus);
    myfree(irequests);
}

void
mpi_profiles_worker (world_fmt * world, long *gmaxptr)
{
    boolean done = FALSE;
    long pnum;
    long temp[4];
    char *tempstr;
    char *quantilebuffer;
    char *failedbuffer;
    MPI_Status status;
    quantilebuffer = (char *) mycalloc (ONE, sizeof (char));
    failedbuffer = (char *) mycalloc (ONE, sizeof (char));
    tempstr = (char *) mycalloc (STRSIZE, sizeof (char));
    while (!done)
    {
        //fprintf(stdout,"%i> before receive of parameter number\n",myID);
        MYMPIRECV (&pnum, 1, MPI_LONG, (MYINT) MASTER, (MYINT) MPI_ANY_TAG, comm_world, &status);
        //fprintf(stdout,"%i> RECEIVED parameter number %li from %i with tag %i\n",myID,pnum,status.MPI_SOURCE,status.MPI_TAG);
                
        if (status.MPI_TAG != 0) //stop condition
        {
            // fills world->buffer with profile information
            print_profile_likelihood_driver (pnum, world, gmaxptr);
            temp[0] = pnum;
            temp[1] = (long) strlen (world->buffer);
#ifdef PRETTY
	    // the PDF table is also sent in addition to the formatted ascii table
	    temp[1] = (long) strlen(world->buffer + temp[1] + 1) + temp[1] + 1;
#endif
            if(world->options->printprofsummary)
            {
                temp[2] = pack_quantile (&quantilebuffer, world->quantiles[pnum], GRIDSIZE);
                world->buffer =     (char *) myrealloc (world->buffer,
                                                        sizeof (char) * (temp[1] + temp[2] + 1));
                sprintf(world->buffer + temp[1], "%s", quantilebuffer);
                temp[3] = pack_failed_percentiles (&failedbuffer, world->percentile_failed[pnum], GRIDSIZE);
                world->buffer =     (char *) myrealloc (world->buffer,
                                                        sizeof (char) * (temp[1] + temp[2] + temp[3] + 1));
                sprintf(world->buffer+temp[1]+temp[2], "%s", failedbuffer);
            }
            else
            {
                temp[2] = 0;
                temp[3] = 0;
            }
            sprintf(tempstr,"P%li", temp[1]);
            //tempstr[0]='P';
            MYMPISEND (tempstr, STRSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) pnum + 1, comm_world);
            MYMPISEND (temp, FOUR, MPI_LONG, (MYINT) MASTER, (MYINT) pnum + 1, comm_world);
            MYMPISEND (world->buffer, temp[1] + temp[2] + temp[3], MPI_CHAR, (MYINT) MASTER, (MYINT) 
                      pnum + 1, comm_world);
            world->profilewho[profiledone++] = pnum;
        }
        else
        {
            done = TRUE;
        }
    }
    myfree(tempstr);
    myfree(quantilebuffer);
    myfree(failedbuffer);
}

long
pack_quantile (char **buffer, quantile_fmt quant, long n)
{
    long i;
    char fp[LONGLINESIZE];
    long bufsize = LINESIZE;
    *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
    sprintf (*buffer, "QUANTILEBUFFER:\n %s\n", quant.name);
    for (i = 0; i < n; i++)
    {
        bufsize += 1 + sprintf (fp, "%20.20f\n", quant.param[i]);
        *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
        strcat (*buffer, fp);
    }
    bufsize = (long) strlen(*buffer);
    return bufsize;
}

void
unpack_quantile (char *buffer, quantile_fmt quant, long n)
{
    long i;
    char *input;
    char *buf = buffer;
    input = (char*) mycalloc(LONGLINESIZE,sizeof(char));
    sgets (input, LONGLINESIZE, &buf);
    sgets (input, LONGLINESIZE, &buf);
    strcpy (quant.name, input);
    for (i = 0; i < n; i++)
    {
        sgets (input, LONGLINESIZE, &buf);
        quant.param[i] = atof (input);
    }
    myfree(input);
}

/// 
/// pack notice of failure of convergence to the profile likelihood precentiles
/// this assume that n is never bigger than LONGLINESIZE, a safe assumption
/// n is the number of grid points in the profile calculation, currently set to 9
/// (May 19 2004), changing this number will cause large ripple effects. but see
/// under profile_max_precentile()
long
pack_failed_percentiles (char **buffer, boolean *failed, long n)
{
    long i;
    char fp[LONGLINESIZE];
    long bufsize = n + ONE;
    *buffer = (char *) myrealloc (*buffer, sizeof (char) * bufsize);
    memset(*buffer,0,sizeof(char)*bufsize);
    for (i = 0; i < n; i++)
        fp[i] =  failed[i] ? '1' : '0' ;
    fp[i]='\0';
    strcat (*buffer, fp);
    return bufsize;
}

/// 
/// unpack notice of failure of convergence to the profile likelihood precentiles
/// this assume that n is never bigger than LONGLINESIZE, a safe assumption
void
unpack_failed_percentiles (char *buffer, boolean *failed, long n)
{
    long i;
    char *input;
    char *buf = buffer;
    input = (char*) mycalloc(LONGLINESIZE,sizeof(char));
    sgets (input, LONGLINESIZE, &buf);
    //fprintf(stdout,"@%s@\n",input);
    for (i = 0; i < n; i++)
    {
        failed[i] = (input[i] == '1');
    //    fprintf(stdout,"@%i\n",(int) failed[i]);
    }
    myfree(input);
}

#endif

/*
// send the data over all loci/replicates to all nodes
// including the master node, so that all nodes can then 
// start calculating profiles [see calc_profiles()]
//
void distribute_locidata(world_fmt *world)
{
  char *buffer;
  pack_loci_data(world, &buffer);
  MPI_allgather(buffer);
  unpack_loci_data(buffer, world);
  myfree(buffer);
}
 
void pack_loci_data(world_fmt *world, char **buffer)
{
  long replicates = world->options->repl
  *buffer = myrealloc(*buffer,LONGLINESIZE);
  hits = sscanf (input, "%li %li %li %li %li", &world->loci, &world->numpop, &world->numpop2, &tmp, &replicates);  
}
*/

// necessary for analyzing old sumfiles using MPI
//
// master is reusing  mpi_runloci_master()
void
assignloci_worker (world_fmt * world, option_fmt *options, long * Gmax)
{
    boolean done = FALSE;
    long locus;
    MPI_Status status;
    long * twolongs;
    char *locusstring;
    twolongs = (long *) mycalloc(TWO,sizeof(long));
    locusstring = (char *) mycalloc(SMALLBUFSIZE,sizeof(char));

    world->options->progress=FALSE;
    options->progress=FALSE;
    
    while (!done)
    {
        MYMPIRECV (twolongs, TWO, MPI_LONG, (MYINT) MASTER, (MYINT) MPI_ANY_TAG,
		   comm_world, &status); //from mpi_runloci_master() around line migrate_mpi.c:163
        if (status.MPI_TAG != 0) //stop condition
        {
	  locus = twolongs[0];
	  sprintf(locusstring,"R%li",locus);
#ifdef DEBUG_MPI
	  printf("%i>>>>>> received locus %li in assignloci_worker{}\n",myID,locus);
	  swap_atl (locus, locidone, world);
	  printf("%i>>>>>> will send locus %li (%s) in assignloci_worker{}\n",myID,locus,locusstring);
#else
	  swap_atl (locus, locidone, world);
#endif
	  MYMPISEND (locusstring, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) locus + 1, comm_world);
	  /* we want to know what locus we worked for
	     - to control the work sent by master
	     - to use in setup_parameter0() [combroyden2.c] */
	  world->who[locidone++] = locus;
	  //
	  if (options->replicate && options->replicatenum > 0)
            {
	      world->locus = locus;
	      world->repkind = MULTIPLERUN;
#ifdef LONGSUM
	      change_longsum_times (world);   //multi run replicates
#endif         /*LONGSUM*/
	      if (!options->bayes_infer)
		{
		  (void) estimateParameter (options->replicatenum, *Gmax, world,
					    options, world->cov[locus], 1, 'l',
					    SINGLELOCUS, world->repkind);
		}
            }

	  //
        }
        else
	  {
            done = TRUE;
	    assign_worker_cleanup();
#ifdef DEBUG_MPI
	    fprintf(stdout,"%i> STOP: received stop from %i\n",myID, status.MPI_SOURCE);
#endif
	  }
    }
    myfree(locusstring);
}

void
swap_atl (long from, long to, world_fmt * world)
{
    long r;
    timearchive_fmt *tmp;
    for (r = 0; r < world->options->replicatenum; r++)
    {
        tmp = &world->atl[r][to];
        world->atl[r][to] = world->atl[r][from];
        world->atl[r][from] = *tmp;
    }
}


#ifdef SLOWNET
void
setup_parameter0_slowmpi (world_fmt * world, nr_fmt * nr, long repkind,
                          long repstart, long repstop, long loci, long kind,
                          boolean multilocus)
{
    long locus, r;
    if (myID != MASTER)
    {
        if (multilocus)
        {
            for (locus = 0; locus < loci; locus++)
            {
                if (repkind == SINGLECHAIN)
                {
                    for (r = repstart; r < repstop; r++)
                        create_apg0 (nr->apg0[r][locus], nr,
                                     &world->atl[r][locus], locus);
                }
                else
                {
//                    if (kind != PROFILE)
//                    {
                        for (r = repstart; r < repstop; r++)
                            create_apg0 (nr->apg0[r][locus], nr,
                                         &world->atl[r][locus], locus);
                        interpolate_like (nr, locus);
//                    }
//                    else
//                    {
                        for (r = repstart; r < repstop; r++)
                            create_multiapg0 (nr->apg0[r][locus], nr, r, locus);
//                    }
                }
            }
        }
        else   //single locus
        {
            if (repkind == SINGLECHAIN)
            {
                for (r = repstart; r < repstop; r++)
                    create_apg0 (nr->apg0[r][world->locus], nr,
                                 &world->atl[r][world->locus], world->locus);
            }
            else
            {
//                if (kind != PROFILE)
//                {
                    for (r = repstart; r < repstop; r++)
                        create_apg0 (nr->apg0[r][world->locus], nr,
                                     &world->atl[r][world->locus], world->locus);
                    interpolate_like (nr, world->locus);
//                }
                for (r = repstart; r < repstop; r++)
                    create_multiapg0 (nr->apg0[r][world->locus], nr, r,
                                      world->locus);
            }
        }
    }
}
#endif

void
setup_parameter0_mpi (world_fmt * world, nr_fmt * nr, long repkind,
                      long repstart, long repstop, long loci, long kind,
                      boolean multilocus)
{
    long locus, r;
    long ll;
    if (myID != MASTER)
    {
        if (multilocus)
        {
            for (ll = 0; ll < locidone; ll++)
            {
                locus = world->locus = world->who[ll];
                if (repkind == SINGLECHAIN)
                {
                    for (r = repstart; r < repstop; r++)
                        create_apg0 (nr->apg0[r][locus], nr,
                                     &world->atl[r][locus], locus);
                }
                else
                {
//                    if (kind != PROFILE)
//                    {
                        for (r = repstart; r < repstop; r++)
                            create_apg0 (nr->apg0[r][locus], nr,
                                         &world->atl[r][locus], locus);
                        interpolate_like (nr, locus);
//                    }
//                    else
//                    {
                        for (r = repstart; r < repstop; r++)
                            create_multiapg0 (nr->apg0[r][locus], nr, r, locus);
//                    }
                }
            }
        }
        else   //single locus
        {
            if (repkind == SINGLECHAIN)
            {
                for (r = repstart; r < repstop; r++)
                    create_apg0 (nr->apg0[r][world->locus], nr,
                                 &world->atl[r][world->locus], world->locus);
            }
            else
            {
//                if (kind != PROFILE)
//                {
                    for (r = repstart; r < repstop; r++)
                        create_apg0 (nr->apg0[r][world->locus], nr,
                                     &world->atl[r][world->locus], world->locus);
                    interpolate_like (nr, world->locus);
//                }
                for (r = repstart; r < repstop; r++)
                    create_multiapg0 (nr->apg0[r][world->locus], nr, r,
                                      world->locus);
            }
        }
    }
}

///
/// checks whether a node is already in the mpistack list
boolean in_mpistack(int sender, world_fmt *world)
{
  long i;
  for (i=0;i < world->mpistacknum; i++)
    {
      if((world->mpistack[i]) == sender)
	return TRUE;
    }
  return FALSE;
}

void handle_replication(int sender,int tag,char *tempstr, world_fmt *world)
{
  //long  pos=0;
  long i;
  long * temp;
  int realsender;
  long locus;
  long replicate;
  int replicator;
  //  boolean from_locus_sender = (sender <= world->loci);
  //  boolean from_replicator = (sender > world-> loci);
  sscanf(tempstr+1,"%i%li%li", &realsender, &locus, &replicate);
  //if(from_replicator)
  //  warning("these guys should not send to here");
  temp = (long *) mycalloc(3, sizeof(long));
  temp[0] = sender;
  temp[1] = locus;
  temp[2] = replicate;
  if(world->mpistacknum > 0)
    {
      world->mpistacknum -= 1;
      replicator = world->mpistack[world->mpistacknum];
      //printf("%i> checking out replicator:[%4li] %i\n",myID, world->mpistacknum,replicator);
      MYMPISEND(temp, 3, MPI_LONG, replicator, (MYINT) tag, comm_world);
    }
  else
    {
      if(world->mpistack_requestnum>=world->mpistack_request_numalloc)
	{

	  world->mpistack_request = (mpirequest_fmt *) myrealloc(world->mpistack_request,
					      sizeof(mpirequest_fmt)*(world->mpistack_requestnum+10));
	  for(i=world->mpistack_request_numalloc; i < world->mpistack_requestnum + 10; i++)
	    {
	      	world->mpistack_request[i].tempstr = (char *) mycalloc(SMALLBUFSIZE,sizeof(char));
	    }
	  world->mpistack_request_numalloc = world->mpistack_requestnum + 10;
	}
      strcpy(world->mpistack_request[world->mpistack_requestnum].tempstr, tempstr);
      world->mpistack_request[world->mpistack_requestnum].sender = sender;
      world->mpistack_request[world->mpistack_requestnum].tag = tag;
      world->mpistack_requestnum += 1;
      //      printf("%i> MASTER: added request from n%i to request-stack position %li\n",myID, sender, world->mpistack_requestnum);
    }
  myfree(temp);
}

void handle_mdim(float *values,long n, int sender, world_fmt * world)
{
    register long i;
    register long z;

    int digits;
#ifdef ZNZ
    znzFile file = world->bayesmdimfile;
#else
    FILE *file = world->bayesmdimfile;
#endif
    register char *temp;

    temp = (char *) mycalloc(n*20,sizeof(char)); 
    z = sprintf(temp,"%li\t%li\t%li\t%f\t%f\t%f\t%f\t%li\t%f",
	   (long) values[0], (long) values[1]+1,(long) values[2]+1,
	    values[3], values[4], values[5], values[6],(long) values[7], values[8]);
    for(i=9;i<n; i++)
      {
	//fprintf (stdout, "%f (%li) ",values[i], i);
#ifdef BFDEBUG
	if(i > (world->numpop2 +  world->bayes->mu))
	  {
	    z += sprintf(temp+z,"\t%f",values[i]);
	    continue;
	  }
#endif
	if(values[i] < SMALLEPSILON)
	  {
	    z += sprintf(temp+z,"\t0");
	  }
	else
	  {
	    digits = (long) log10(values[i]);
	    switch(digits)
	      {
	      case -8:
	      case -6:
	      case -5:
		//		fmt = 10;
		z += sprintf(temp + z,"\t%.10f", values[i]);
		break;
	      case -4:
	      case -3:
		//fmt = 8;
		z += sprintf(temp + z,"\t%.8f", values[i]);
		break;
	      case -2:
	      case -1:
		//		fmt= 5;
		z += sprintf(temp + z,"\t%.5f", values[i]);
		break;
	      case 0:
	      case 1:
		//		fmt = 4;
		z += sprintf(temp + z,"\t%.4f", values[i]);
		break;
	      case 2:
		//		fmt = 2;
		z += sprintf(temp + z,"\t%.2f", values[i]);
		break;
	      case 3:
		//		fmt = 1;
		z += sprintf(temp + z,"\t%.1f", values[i]);
		break;
	      case 4:
	      case 5:
	      case 6:
	      case 7:
	      case 8:
		//		fmt = 0;
		z += sprintf(temp + z,"\t%.0f", values[i]);
		break;
	      default:
		if(digits<-8)
		  {
		    //		    fmt=20;
		    z += sprintf(temp + z,"\t%.20f", values[i]);
		  }
		else
		  {		   
		    //		    fmt = 5;
		    z += sprintf(temp + z,"\t%.5f", values[i]);
		  }
		break;
	      }
	    //	    fprintf(stdout,"\t%f", values[i]);
	  }
      }
    //    fprintf (stdout, " [%i] \n", myID);
#ifdef ZNZ
    znzprintf(file,"%s\n",temp);
#else
    fprintf(file,"%s\n",temp);
    fflush(file);
#endif
    myfree(temp);
    //fprintf(stdout,"\n");
    // calculate_parallel_convergence(world, values, size);
}

void handle_message(char *rawmessage,int sender, world_fmt * world)
{
    char *rawptr;
    long  pos=0;
    void *file = (void *) stdout;
    rawptr = rawmessage;
    //fprintf(stderr,"%i> handle_message: %s\n",myID, rawmessage);
    set_filehandle(rawmessage, world, &file, &pos);
    //fprintf(stderr,"%i> handle_message: pos=%li\n",myID,pos);
    //fprintf(stderr,"%i> handle_message: %s\n",myID, rawmessage + pos);
    fprintf((FILE *) file,"%s", rawptr + pos);
    fflush((FILE *) file);
}

void handle_burnin_message(char *rawmessage,int sender, world_fmt * world)
{
  static long       z = 0;
  static boolean done = FALSE;
  static long    *replicates;
  long locus;
  long step;
  MYREAL var;
  MYREAL oldvar;
  if(!done)
    {
      replicates = mycalloc(world->loci,sizeof(long));
      done = TRUE;
    }
#ifdef USE_MYREAL_FLOAT
  sscanf(rawmessage,"%li%f%f%li",&locus, &var, &oldvar, &step);
#else
  sscanf(rawmessage,"%li%lf%lf%li",&locus, &var, &oldvar, &step);
#endif
  replicates[locus] += 1;
  if(world->options->verbose)
    {
      fprintf(stdout,"%i> Burn-in on node %i (locus=%li, repl=%li) stopped at step %li with variance-ratio=%.2f/%.2f=%.3f\n", myID, sender, locus, replicates[locus], step, var, oldvar, var/oldvar);
    }
  world->burnin_stops[z].locus        =  locus; 
  world->burnin_stops[z].replicate =  replicates[locus]; 
  world->burnin_stops[z].stopstep = step; 
  world->burnin_stops[z].variance = var; 
  world->burnin_stops[z].oldvariance = oldvar;
  world->burnin_stops[z].worker = sender;
  z++;
}

///
/// sets up a file database so that the master and worker worker end up writing to the same file
/// the workers send the stuff to the master and the master then figures out (using the db)
/// what file pointer the stuff was intended for.
/// needs globals filedb, and filenum
void setup_filehandle_db(FILE *file, world_fmt *world, option_fmt *options, data_fmt *data)
{
    long filehandle = get_filehandle(file, world, options, data);
    filedb[filenum].file = file;
    filedb[filenum++].handle = filehandle;
#ifdef DEBUG_MPI
    fprintf(stdout,"filedb %li: %p %li\n",filenum, file,filehandle);
#endif
}

long retrieve_filehandle(FILE *file)
{
    long i=0;
    long filehandle = 0;
    while(filedb[i].file != file && i<filenum)
        i++;
    if(i!=filenum)
        filehandle = filedb[i].handle;
    return filehandle;
}

long get_filehandle(void *vfile, world_fmt *world, option_fmt *options, data_fmt *data)
{
#ifdef ZNZ
  if(((znzFile) vfile) == world->bayesmdimfile)
        return BAYESMDIMFILENUM;
    FILE * file = (FILE *) vfile;
#else
    FILE * file = (FILE *) vfile;
    if(file == world->bayesmdimfile)
        return BAYESMDIMFILENUM;
#endif
    if(file == stdout)
        return STDOUTNUM;
    if(file == options->logfile)
        return LOGFILENUM;
    if(file == world->outfile)
        return OUTFILENUM;
    if(file == options->aicfile)
        return AICFILENUM;
    if(file == world->mathfile)
        return MATHFILENUM;
    if(file == world->mighistfile)
        return MIGHISTFILENUM;
    if(file == world->skylinefile)
        return SKYLINEFILENUM;
    if(file == world->bayesfile)
        return BAYESFILENUM;
    if(file == world->pdfoutfile)
        return PDFOUTFILENUM;
    if(file == world->treefile)
        return TREEFILENUM;
    return STDOUTNUM;
}

long get_filehandle2(void *vfile, world_fmt *world)
{
  FILE *file ;
#ifdef ZNZ
  if(((znzFile) vfile) == world->bayesmdimfile)
        return BAYESMDIMFILENUM;
#else
    if(((FILE *) vfile) == world->bayesmdimfile)
        return BAYESMDIMFILENUM;
#endif
    file = (FILE*) vfile;
    if(file == stdout)
        return STDOUTNUM;
    if(file == world->options->logfile)
        return LOGFILENUM;
    if(file == world->outfile)
        return OUTFILENUM;
    if(file == world->mighistfile)
        return MIGHISTFILENUM;
    if(file == world->skylinefile)
        return SKYLINEFILENUM;
    if(file == world->options->aicfile)
        return AICFILENUM;
    if(file == world->mathfile)
        return MATHFILENUM;
    if(file == world->bayesfile)
        return BAYESFILENUM;
    if(file == world->treefile)
        return TREEFILENUM;
 //   fprintf(stdout,"@@@@@@@@@@@@@@@wrong wrong wrong@@@@@@@@@@@@@@@@\n");
    return STDOUTNUM;
}

void set_filehandle(char *message, world_fmt *world,
                    void **file, long *msgstart)
{
    char *temp;
    long filenum;
    long i = 1;
    temp = (char *) mycalloc(10,sizeof(char));
    if(message[0] == '\0')
      {
	warning("%i> set_filehandle() the message was NULL",myID);
      }
    temp[0] = message[i];
    while(temp[i-1]!=':' && i < 9 && message[i]!='\0')
      {
#ifdef DEBUG_MPI
	fprintf(stderr,"%i>>>>>> temp     =%s\n",myID,temp);
	fprintf(stderr,"%i>>>>>> temp[%li]=%c\n",myID,i,temp[i]);
	fprintf(stderr,"%i>>>>>> temp     =%s\n",myID,message);
#endif
	i++;
	temp[i-1] = message[i];
      }
    *msgstart = i+1;
    filenum = atol(temp);
    //    fprintf(stdout,"\n@@@@@@@@@@@@@@@@@%li@%s@%li@\n",filenum,temp,i);
    myfree(temp);
    switch(filenum)
      {
      case STDOUTNUM:
	{
	  //		fprintf(stdout,"\n");
	  *file = stdout;
	  return;
	}
      case LOGFILENUM:
	{
	  //	fprintf(stdout," logfile\n");
	  *file = world->options->logfile;
	  return;
	}
      case OUTFILENUM:
	{
	  *file = world->outfile;
	  return ;
	}
      case AICFILENUM:
	{
	  *file = world->options->aicfile;
	  return ;
	}
      case MATHFILENUM:
	{
	  *file = world->mathfile;
	  return ;
	}
      case MIGHISTFILENUM:
	{
	  *file = world->mighistfile;
	  return ;
	}
      case SKYLINEFILENUM:
	{
	  *file = world->skylinefile;
	  return ;
	}
      case BAYESFILENUM:
	{
	  *file = world->bayesfile;
	  return ;
	}
      case BAYESMDIMFILENUM:
	{
	  *file = world->bayesmdimfile;
	  return ;
	}
      case TREEFILENUM:
	{
	  *file = world->treefile;
	  return ;
	}
      case PDFOUTFILENUM:
	{
	  *file = world->pdfoutfile;
	  return ;
	}
      }
    *file = stdout;
    return;
}


void
mpi_fprintf(FILE *file, const char *fmt, ...)
{
    char *p1;
    char *p;
    va_list ap;
    long filehandle = 0;
    long bufsize = 0;
    
    long pallocsize = LINESIZE+strlen(fmt)+1;
    p  = (char *) mycalloc(pallocsize,sizeof(char));
    p1 = (char *) mycalloc(STRSIZE,sizeof(char));
    if(myID!=MASTER)
    {
        filehandle = retrieve_filehandle(file);
        bufsize += sprintf(p, "%c%li:",'M',filehandle);
    }
    va_start(ap, fmt);
    bufsize += vsprintf(p+bufsize, fmt, ap);
    if(bufsize>=pallocsize)
      error("failed in mpi_printf(): problem with buffer size!");
    if(myID!=MASTER)
    {
        sprintf(p1,"M%li",bufsize);
        MYMPISEND (p1, STRSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
        MYMPISEND (p, bufsize, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
    }
    else
        fprintf(file,"%s", p);
    va_end(ap);
    myfree(p);
    myfree(p1);
}
void
mpi_fprintf2(FILE *file, long filesize, const char *fmt, ...)
{
    char *p1;
    char *p;
    va_list ap;
    long filehandle = 0;
    long bufsize = 0;
    
    long pallocsize = filesize+strlen(fmt)+10;//leave room for "M:number"
    
    p  = (char *) mycalloc(pallocsize,sizeof(char));
    p1 = (char *) mycalloc(STRSIZE,sizeof(char));
    if(myID!=MASTER)
    {
        filehandle = retrieve_filehandle(file);
        bufsize = sprintf(p, "%c%li:",'M',filehandle);
    }
    va_start(ap, fmt);
    bufsize += vsprintf(p+bufsize, fmt, ap);
    if(bufsize>=pallocsize)
      {
	warning("Failing because bufsize=%li >= allocsize=%li\n",bufsize,pallocsize);
	error("failed in mpi_printf2()");
      }
    if(myID!=MASTER)
    {
        sprintf(p1,"M%li",bufsize);
        MYMPISEND (p1, STRSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
        MYMPISEND (p, bufsize, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
    }
    else
        fprintf(file,"%s", p);
    va_end(ap);
    myfree(p);
    myfree(p1);
}

///
/// sends raw bayesian parameters to master, using label 'Z' to match on the master side
#ifdef PARALIO
void mpi_mdim_send(MPI_File *file, float *values, long size)
#else
void mpi_mdim_send(float *values, long size)
#endif
{
#ifdef PARALIO
  MPI_Request request;
#endif
    char *p1;
    p1 = (char *) mycalloc(SMALLBUFSIZE,sizeof(char));
    if(myID!=MASTER)
    {
#ifdef PARALIO
      // TEST parallele I/O
      MPI_File_iwrite_shared(*file,values,size,MPI_FLOAT, &request);
#else
      sprintf(p1,"Z%li",size);
      MYMPISEND (p1, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
#ifdef DEBUG_MPI
	fprintf(stdout,"%i> mdimlast=%f\n",myID,values[size-1]);
#endif
	MYMPISEND (values, size, MPI_FLOAT, (MYINT) MASTER, (MYINT) myID+PRINTTAG, comm_world);
#endif
    }
    else
        error("master sends itself bayesallfile stuff");
    myfree(p1);
}


///
/// assembles the data from a replicant (receives materials from mpi_send_replicant()
/// 
void 
mpi_receive_replicate(int sender, int tag, long locus, long replicate, world_fmt * world)
{
    MYREAL *buffer;
    long  bufsize=1;    
    MPI_Status status;
    long z;
    MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
#ifdef DEBUG_MPI
    fprintf(stdout,"%i> WORKER: mpi_receive_replicate received bufsize=%li from sender=%i with tag=%i\n",myID, bufsize, status.MPI_SOURCE, status.MPI_TAG);    
    sender = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    fprintf(stdout,"%i> mpi_receive_replicate received bufsize=%li from sender=%i with tag=%i\n",myID, bufsize, sender, tag);    
#endif
    buffer = (MYREAL *) mycalloc (bufsize, sizeof (MYREAL));
    MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
    //fprintf(stdout,"%i> received bufsize is really %li and bufsize=%li\n",myID,(long) (long) strlen(buffer),bufsize);
    if(world->options->bayes_infer)
    {
        unpack_single_bayes_buffer(buffer,world->bayes,world,locus);
    }
    else
    {
      z=0;
      unpack_single_sumfile_buffer(buffer, world->atl, world, locus, replicate, world->numpop, &z);
      //      fprintf(stdout,"%i> replicate=%li locus=%li\n",myID, replicate, locus);
    }
    if(world->options->mighist)
    {
      MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
#ifdef DEBUG_MPI
      fprintf(stdout,"%i> mpi_receive_replicate received mighistogram bufsize=%li from sender=%i with tag=%i\n",
	      myID, bufsize, status.MPI_SOURCE, status.MPI_TAG);    
#endif
      buffer = (MYREAL *) myrealloc (buffer, bufsize * sizeof (MYREAL));
      MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
      unpack_mighist_replicate_buffer(buffer, world, locus, world->numpop);
    }
    if(world->options->mighist && world->options->skyline)
    {
      MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
      buffer = (MYREAL *) myrealloc (buffer, bufsize * sizeof (MYREAL));
      MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
      unpack_skyline_buffer(buffer, world, locus, -1, world->numpop);
    }
  // send best tree if available
    //  if(world->options->treeprint == BEST && world->options->treeinmemory == TRUE)
  if(world->options->treeinmemory == TRUE)
    {
      MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
      buffer = (MYREAL *) myrealloc (buffer, bufsize * sizeof (MYREAL));
      //      printf("%i> mpi_receive_buffer() received treespace buffersize: %li %p\n", myID, bufsize,  buffer);
      MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
      unpack_treespace_buffer(buffer, world, locus, -1, world->numpop);      
    }

    // BF material
  MYMPIRECV (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) tag, comm_world, &status);
  buffer = (MYREAL *) myrealloc (buffer, bufsize * sizeof (MYREAL));
  MYMPIRECV (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) tag, comm_world, &status);
  if(world->options->bayes_infer)
    {
      if(!world->data->skiploci[locus])
	{
	  //printf("%i> mpi_receive_buffer() received BF buffersize: %li %p\n", myID, bufsize,  buffer);
	  bufsize = unpack_BF_buffer(buffer, 0, locus, world);
	}
    }
  // ESS material
  bufsize = unpack_ess_buffer(buffer, bufsize, world);


    myfree(buffer);
}

///
/// replicant sends data to sub-master
/// 
void 
mpi_send_replicate(int sender, long locus,  long replicate, world_fmt * world)
{
  long    allocbufsize   = 1;
  long    bufsize        = 0;
  long    numpop         = world->numpop;
  long    numpop2        = numpop * numpop;
  long    numpop2plus    = numpop2 + 2 * numpop;
  long    npp            = numpop2 + world->bayes->mu * world->loci;
  MYREAL  *buffer        = NULL;
  timearchive_fmt **ta   = world->atl;
  
  if(world->options->bayes_infer)
    {
	  buffer = (MYREAL *) mycalloc (ONE, sizeof (MYREAL));
	  bufsize = pack_single_bayes_buffer(&buffer,world->bayes,world,locus);
    }
  else
    {
      bufsize = ta[replicate][locus].T;
      bufsize *= numpop2plus;
      bufsize += 2 * ta[replicate][locus].T;
      bufsize += 2 * numpop2;
      bufsize += 6;
      allocbufsize = bufsize;
      buffer = (MYREAL *) mycalloc (allocbufsize, sizeof (MYREAL));    
      bufsize = pack_single_sumfile_buffer(&buffer, 0, world, locus, replicate, world->numpop);
    }
  MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
  MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
  if(world->options->mighist)
    {
      bufsize = pack_mighist_buffer(&buffer, world, locus, -1, numpop);
      MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
      MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
    }

  if(world->options->mighist && world->options->skyline)
    {
      buffer = (MYREAL *) myrealloc (buffer, allocbufsize *  sizeof (MYREAL));    
      bufsize = pack_skyline_buffer(&buffer, world, locus, -1, numpop);
      MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
      MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
    }
  // send best tree if available
  if(/*world->options->treeprint == BEST &&*/ world->options->treeinmemory == TRUE)
    {
      //      printf("%i> send treespace buffer to %li", myID,locus+1+REPTAG);
      bufsize = pack_treespace_buffer(&buffer, world, locus, -1, numpop);
      MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
      MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
    }

    // BF material
  if(world->options->bayes_infer)
    {
      bufsize = world->options->heated_chains * world->loci + 5 * world->loci + 20 * (npp+1);
      buffer = (MYREAL *) myrealloc (buffer, bufsize *  sizeof (MYREAL));
      if(!world->data->skiploci[locus])
	{
	  //fprintf(stderr,"%i> REPLICANT: packed result locus=%li replicate %li\n",myID,locus, replicate);
	  bufsize = pack_BF_buffer(&buffer, 0, locus, world);
	}
    }
  // ESS material
  bufsize = pack_ess_buffer(&buffer, bufsize, world);
  // send BF and ESS material
  MYMPISEND (&bufsize, ONE, MPI_LONG, (MYINT) sender, (MYINT) (locus+1+ REPTAG), comm_world);
  MYMPISEND (buffer, bufsize, mpisizeof, (MYINT) sender, (MYINT) (locus+1 + REPTAG), comm_world);
  
  myfree(buffer);
}

#ifdef MPICHECK
void set_memory_limit(rlim_t softsize,rlim_t maxsize)
{
  struct rlimit r;
  r.rlim_cur=softsize;
  r.rlim_max=maxsize;
  setrlimit(RLIMIT_AS, &r);
}

void check_memory_limit()
{
  struct rlimit r;
  getrlimit(RLIMIT_AS, &r);
  fprintf(stdout, "%i> current memory/data usage: %f\n",myID, (double) r.rlim_cur);
}
#endif
#endif /* MPI */
