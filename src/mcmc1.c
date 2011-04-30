/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effective population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 M C M C   R O U T I N E S 
 
 Markov Monte Carlo stuff: treechange, acceptance
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
 Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
 Copyright 2003-2006 Peter Beerli, Tallahassee FL
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 $Id: mcmc1.c 1807 2011-03-18 20:16:19Z beerli $
 -------------------------------------------------------*/
/* \file mcmc1.c
 
 Tree changer and acceptance rejection scheme
 
 */
#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#include "random.h"
#include "tree.h"
#include "mcmc2.h"
#ifdef UEP
#include "uep.h"
#endif

#include "bayes.h"

#ifdef BEAGLE
#include "calculator.h"
#endif

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

#define MIGRATION_AIR (boolean) 1
#define MIGRATION_IN_TREE (boolean) 0
#define NO_MIGR_NODES 0
#define WITH_MIGR_NODES 1

boolean debugtip;

extern int myID;
extern MYREAL probg_treetimesX(world_fmt *world, vtlist *tl, long T);

extern void     zero_xseq(xarray_fmt *x, long sites, long categs);

/* prototypes ------------------------------------------- */
// metropolize over trees
long tree_update (world_fmt * world, long g);
/* private functions */
void new_localtimelist (timelist_fmt ** ntl, timelist_fmt * otl, long numpop);
void new_proposal (proposal_fmt ** proposal, timelist_fmt * tl,
                   world_fmt * world);
void set_new_proposal (proposal_fmt ** proposal, timelist_fmt * tl, world_fmt * world);
void chooseOrigin (proposal_fmt * proposal);
void construct_localtimelist (timelist_fmt * timevector,
                              proposal_fmt * proposal);
void traverseAllNodes (node * theNode, node *** nodelist, long *node_elem,
                       long *oldnode_elem, int include_migration);
void chooseTarget (proposal_fmt * proposal, timelist_fmt * timevector,
                   node ** bordernodes, long *bordernum);
void findbordernodes (node * theNode, proposal_fmt * proposal, long pop,
                      node *** bordernodes, long *allocsize, long *bordernum, vtlist ** tyme,
                      long gte);
void free_proposal (proposal_fmt * proposal);
void free_timevector (timelist_fmt * timevector);
void prune_timelist (timelist_fmt * oldtv, timelist_fmt * newtv, register node ** __restrict ptr, proposal_fmt * proposal, long numpop);

long xor (node ** ptrl1, node ** ptrl2);
long rmigrcount (proposal_fmt * proposal);

int migrate_old (proposal_fmt * proposal, node * up,
                 long *old_migr_table_counter, boolean air);
int migrate (proposal_fmt * proposal, node * up);
int migrateb (proposal_fmt * proposal, node * up);

int pre_population (proposal_fmt * proposal, vtlist * ltime, long gte,
                    long *slider);

boolean acceptlike (world_fmt * world, proposal_fmt * proposal, long g,
                    timelist_fmt * tyme);
MYREAL eventtime (proposal_fmt * proposal, long pop, vtlist * tentry,
                  char *event);
node *showsister (node * theNode);
void count_migrations (node * p, long *count);
long migration_from (long to, proposal_fmt * proposal);

MYREAL prob_tree (world_fmt * world, timelist_fmt * tyme);
void traverse_check(node *theNode);
void reset_proposal (proposal_fmt ** proposal, world_fmt *world);


/* Functions ++++++++++++++++++++++++++++++++++++++++++++++++*/
void
set_tree_dirty (node * p)
{
    switch (p->type)
    {
        case 'm':
            set_dirty (p);
            set_tree_dirty (p->next->back);
            break;
        case 't':
            break;
        case 'i':
            set_dirty (p);
            set_tree_dirty (p->next->back);
            set_tree_dirty (p->next->next->back);
            break;
        case 'r':
            set_dirty (p);
            set_tree_dirty (p->next->back);
            break;
    }
}

/*=======================================================*/
long
tree_update (world_fmt * world, long g)
{    /*return 1 if tree was accepted, 0 otherwise */
    //    static long treefilepos; /* write position in the treefile */
#ifdef BEAGLE
    int beagle_control=1;
#endif
    //short show=0;
    boolean coalesced;
    //  boolean test;
    char event;
    long slider;
    long bordernum;
    long actualpop = -99, zz;
    MYREAL endtime, nexttime, age;
#ifdef UEP
    
    boolean uepsuccess = FALSE;
#endif
#ifdef TESTING2
    proposal_fmt *proposal = world->proposal; 
#else
    proposal_fmt *proposal=NULL; 
#endif
    /*scratchpad  in which all involved nodes
     are recorded and help arrays, e.g. migration arrays,
     are stored */
    timelist_fmt *timevector = NULL; /* local timelist */
    vtlist *tentry = NULL; /*pointer into timeslice */
    
    /* ---------------------------------------------------------
     initialize local timelist and construct residual timelist 
     find the start node and snip of the branch and node below */
#ifdef __MWERKS__
    
    eventloop ();
#endif
    new_localtimelist (&timevector, &world->treetimes[0], world->numpop);
#ifdef TESTING2
    if(world->has_proposal_details)
        reset_proposal (&proposal, world);//, &world->treetimes[0], world);
    else
    {
    	set_new_proposal (&proposal, &world->treetimes[0], world);
    	world->has_proposal_details = TRUE;
    }
#else
    new_proposal (&proposal, &world->treetimes[0], world);
#endif
    chooseOrigin (proposal);
    //
    if(world->options->bayes_infer)
    {
        world->bayes->starttime = -proposal->origin->tyme; // we start with negative value to
        // indicate that we do not know yet whether we will accept this move or not
        // if the move is accepted then the time will be negated (->positiv), we record
        // all start and stop times one could easily imagine that some parameter values
        // never gets accepted by the genealogy, and when we record this for bayesallfile we would
        // like to know this. see further doen at the acceptance for the proposals
    }
    //
    construct_localtimelist (timevector, proposal);
    tentry = &(*timevector).tl[0];
    age = proposal->origin->tyme;
    zz = 0;
    while ((tentry->age < age || tentry->age - age < SMALLEPSILON)&& zz < (*timevector).T)
    {
        tentry = &(*timevector).tl[zz];
        zz++;
    }
    zz--;
    nexttime = tentry->age;
    if ((*timevector).T > 1)
        endtime = (*timevector).tl[(*timevector).T - 2].age;
    else
        endtime = 0.0;
    proposal->time = 0.0;
#ifdef TREEDEBUG
    fprintf(stdout,"%i> begin zz=%li nextime=%f (pt=%f+age=%f)=%f,pot=%f\n",myID, zz, nexttime, proposal->time, age, age+proposal->time,proposal->origin->tyme);
#endif
    proposal->time = age;
    coalesced = FALSE;
    /*------------------------------------
     main loop: sliding down the tree  */
    slider = 0;
#ifdef DEBUG
    //printf("\nTime    Pop Type Lineages [mutrate=%f]\n",proposal->world->options->mu_rates[world->locus]);
#endif
    while (nexttime <= endtime)
    {
        actualpop =
        (proposal->migr_table_counter >
         0) ? proposal->migr_table[proposal->migr_table_counter -
                                   1].from : proposal->origin->pop;
#ifdef TREEDEBUG
        if(age == 0.0 && proposal->origin->tyme>0.0)
        {
            printf("\n\n%f\n\n",age);
        }
#endif
        proposal->time = eventtime (proposal, actualpop, tentry, &event);
#ifdef TREEDEBUG
        fprintf(stdout,"%i> zz=%li (pt=%f+age=%f)=%f,pot=%f\n",myID, zz, proposal->time, age, age+proposal->time,proposal->origin->tyme);
#endif
        proposal->time += age;
        //	age = proposal->time;
        
        if(proposal->time < proposal->origin->tyme)
        {
            fprintf(stdout,"%i> Proposal failed because of unordered entry in time list, abort this sample\nProposed time=%f origin time=%f nexttime=%f\n", myID, proposal->time, proposal->origin->tyme, nexttime);
            // we end up here when the migration events exceed the upper limit
            free_proposal (proposal);
            free_timevector (timevector);
            return 0;
            //	    error("Proposal of time failed\n");
        }
        if (proposal->time < nexttime)
        {
            if (event == 'm')
            {
                if (!migrate (proposal, proposal->origin))
                {
                    // we end up here when the migration events exceed the upper limit
                    free_proposal (proposal);
                    free_timevector (timevector);
                    //warning("Event was migration event but could not match a lineage to it\n");
                    return 0;
                }
                age = proposal->time;
                continue;
            }
            else
            {   /*coalesce */
                chooseTarget (proposal, timevector, proposal->bordernodes,
                              &bordernum);
                if(bordernum == 0)
                {
                    if (!migrate (proposal, proposal->origin))
                    {
                        // we end up here when the migration events exceed the upper limit
                        free_proposal (proposal);
                        free_timevector (timevector);
                        return 0;
                    }
                    age = proposal->time;
                    continue;
                    //		    warning("chooseTarget failed at age %f\n", proposal->time);
                    //		    free_proposal (proposal);
                    //		    free_timevector (timevector);
                    //		    return 0;
                }
                pretendcoalesce1p (proposal);
#ifdef UEP
                uepsuccess = is_success_pseudo_uep (proposal);
#endif
                
                coalesced = TRUE;
                break;
            }
        }   /*end if proposal->time < nextime */
        age = nexttime;
#ifdef TREEDEBUG
        if(age < proposal->origin->tyme)
        {
            printf("darn: nextime=%f\n",nexttime);
        }
#endif
        zz++;
        if(zz >= timevector->T)
        {
            break;
        }
        tentry = &(*timevector).tl[zz]; /*next entry in timelist */
        nexttime = tentry->age;
    }
    if (!coalesced)
    {
        if (!pre_population(proposal, (*timevector).tl, (*timevector).T - 1, &slider))
        {
            free_proposal (proposal);
            free_timevector (timevector);
            return 0;
        }
        //	fprintf(stderr,"%i> last time=%f, end time = %f (%c)\n", myID, age,(*timevector).tl[(*timevector).T-1].eventnode->tyme,(*timevector).tl[(*timevector).T-1].eventnode->type );
        pretendcoalesce1p (proposal);
#ifdef UEP
        
        uepsuccess = is_success_pseudo_uep (proposal);
#endif
        
    }
    if (
#ifdef UEP
        ((!world->options->uep && !uepsuccess)
         || (world->options->uep && uepsuccess)) &&
#endif
        (acceptlike (world, proposal, g, timevector)))
    {
#ifdef BEAGLE
        printf("Accepted\n");
        change_beagle(world->root->next->back,world->beagle,world->sumtips);
#endif
        if (proposal->time > world->root->tyme)
        {   /*saveguard */
            world->root->tyme += proposal->time;
        }
        coalesce1p (proposal);
        //DEBUG
        //traverse_check(crawlback (proposal->root->next));printf("*");
#ifdef DEBUG
        //traverse_adjust(world->root->next->back, 1.0);
#endif
        // record the time interval that was used for the lineage.
        if(world->options->bayes_infer)
        {
            world->bayes->starttime = -world->bayes->starttime;
            world->bayes->stoptime = proposal->time;
            world->treelen = 0.0;
            calc_treelength (world->root->next->back, &world->treelen);
        }
        //
#ifdef UEP
        world->likelihood[g] = treelikelihood (world);
        if (world->options->uep)
        {
            update_uep (world->root->next->back, world);
            check_uep_root (world->root->next->back, world);
            world->treelen = 0.0;
            calc_treelength (world->root->next->back, &world->treelen);
            world->ueplikelihood = ueplikelihood (world);
            world->likelihood[g] = /*world->ueplikelihood +*/ world->likelihood[g];
        }
#else
#ifdef BEAGLE
        //	set_beagle_dirty(proposal->origin,proposal->target,showtop(world->root->next->back));
        //reset_beagle(world->beagle);
        //smooth (world->root->next, crawlback(world->root->next), world, world->locus);
        //if(world->beagle->numbranches>0)
        //  {
        world->likelihood[g] = proposal->likelihood; //treelikelihood (world);
        if(beagle_control==1)
        {
            printf("Recalculate LnL using treelikelihood()\n");
            world->likelihood[g] = force_beagle_recalculate(world,world->locus);
            printf("%i> %f = newLnL=%f + recalcLnL=%f\n",myID, proposal->likelihood - world->likelihood[g], proposal->likelihood,world->likelihood[g]);
        }
        //  }
#else
        world->likelihood[g] = treelikelihood (world);
#endif /*BEAGLE*/
#endif /*UEP*/
        /* create a new timelist */
        construct_tymelist (world, &world->treetimes[0]);
        //        if (world->options->treeprint != _NONE && world->cold)
        //  print_tree (world, g, &treefilepos);
        world->migration_counts = 0;
        /* report the number of migration on the tree */
        count_migrations (world->root->next->back, &world->migration_counts);
        free_proposal (proposal);
        free_timevector (timevector);
        return 1;   /* new tree accepted */
    }
    else
    {
        // record the time interval that was used for the lineage.
        if(world->options->bayes_infer)
        {
            world->bayes->stoptime = -proposal->time;
        }
        //traverse_check(crawlback (proposal->root->next));printf(".");
        //
    }
    free_proposal (proposal);
    free_timevector (timevector);
    return 0;   /* not accepted */   
}

#ifdef SLATKIN_IMPORTANCE
long slatkin_importance(world_mt *world, long g)
{
    long slice=0;
    // create timelist with times drawn using the parameters
    new_localtimelist (&timevector, &world->treetimes[0], world->numpop);
    while (slice < sumlineages)
    {
        actualpop =
        (proposal->migr_table_counter >
         0) ? proposal->migr_table[proposal->migr_table_counter -
                                   1].from : proposal->origin->pop;
        age = (*timevector).tl[slice].age;
        (*timevector).tl[slice].age = age + eventtime (proposal, actualpop, tentry, &event);
    }
    // calculate which two datapoints to join give the distance
}
#endif

/*=======================================================*/
///
/// allocate the timelist, contains  times and lineages at that time
/// lineages have a large storage device that is accessed from the 
/// tl[i].lineages.
void
new_localtimelist (timelist_fmt ** ntl, timelist_fmt * otl, long numpop)
{
    //    long i;
    (*ntl) = (timelist_fmt *) mycalloc (1, sizeof (timelist_fmt));
    (*ntl)->tl = (vtlist *) mymalloc ((*otl).allocT * sizeof (vtlist));
    (*ntl)->allocT = otl->allocT;
    (*ntl)->T = otl->T;
    (*ntl)->oldT = otl->oldT;
    //memcpy ((*ntl)->tl, otl->tl, otl->allocT * sizeof (vtlist));
    allocate_lineages (ntl, 0, numpop);
    //memcpy ((*ntl)->lineages, otl->lineages, (*ntl)->allocT * numpop * sizeof (long));
}

///
/// reuse a global timelist to store values, this should save 
/// malloc/free calls
void
new_localtimelist_new (timelist_fmt ** ntl, timelist_fmt * otl, long numpop)
{
    // UNFINISHED
    //    long i;
    //    (*ntl) = (timelist_fmt *) mycalloc (1, sizeof (timelist_fmt));
    if((*ntl)->allocT < otl->allocT)
    {
        (*ntl)->tl = (vtlist *) myrealloc ((*ntl)->tl,(*otl).allocT * sizeof (vtlist));
        //      memset((*ntl)->tl+((*ntl)->allocT),0,sizeof(vtlist)*((*otl)->allocT-(*ntl)->allocT));
        (*ntl)->allocT = otl->allocT;
    }
    (*ntl)->T = otl->T;
    memcpy ((*ntl)->tl, otl->tl, otl->allocT * sizeof (vtlist));
    allocate_lineages (ntl, 0, numpop);
    memcpy ((*ntl)->lineages, otl->lineages, (*ntl)->allocT * numpop * sizeof (long));
}



void
new_proposal (proposal_fmt ** proposal, timelist_fmt * tl, world_fmt * world)
{
#ifdef UEP    
    long j;
#endif
    long mal = world->data->maxalleles[world->locus];
    //long listsize = ((*tl).allocT + (*tl).T + 5);
    long listsize = 2*(world->sumtips + 2);
    long sumtips = world->sumtips;
    
    // allocate the scratchpad (contains a shadow of the tree)
    (*proposal) = (proposal_fmt *) mycalloc (1, sizeof (proposal_fmt));
    (*proposal)->likelihood = -HUGE;
    // pointers and values to outside structures
    (*proposal)->world = world;
    (*proposal)->datatype = world->options->datatype;
    (*proposal)->sumtips = world->sumtips;
    (*proposal)->numpop = world->numpop;
    (*proposal)->endsite = world->data->seq[0]->endsite;
    (*proposal)->fracchange = world->data->seq[0]->fracchange;
    (*proposal)->param0 = world->param0;
    (*proposal)->root = world->root;
    (*proposal)->migration_model = world->options->migration_model;
    // precalculated values
    (*proposal)->mig0list = world->mig0list;
    (*proposal)->design0list = world->design0list;
    
    (*proposal)->listsize =  listsize;
    
    // nodes above the picked node + migration nodes
    (*proposal)->aboveorigin =
    (node **) mycalloc (listsize, sizeof (node *));
    
    // node data holding vector for bordernodes, line_f, line_t
    (*proposal)->nodedata =
    (node **) mycalloc (3 * listsize, sizeof (node *));
    // adjacent nodes of the picked nodes
    (*proposal)->bordernodes = (*proposal)->nodedata;
    // line_f nodes
    (*proposal)->line_f =   (*proposal)->bordernodes + listsize;
    (*proposal)->line_t =  (*proposal)->line_f + listsize;
    // mf holds also mt array  
    (*proposal)->mf = (MYREAL *) mycalloc (2* (*proposal)->endsite, sizeof (MYREAL));
    (*proposal)->mt = (*proposal)->mf + (*proposal)->endsite;
    
    if (strchr (SEQUENCETYPES, (*proposal)->datatype))
    {
        allocate_xseq(&(*proposal)->xf, world->data->seq[0]->endsite, world->options->rcategs);
        allocate_xseq(&(*proposal)->xt, world->data->seq[0]->endsite, world->options->rcategs);
    }
    else
    {
        (*proposal)->xf.a = (MYREAL *) mycalloc (1, sizeof (MYREAL) * mal);
        (*proposal)->xt.a = (MYREAL *) mycalloc (1, sizeof (MYREAL) * mal);
    }
    (*proposal)->old_migr_table_counter = 4 * sumtips /* 100 */ ;
    (*proposal)->old_migr_table_counter2 = 4 * sumtips /* 100 */ ;
    (*proposal)->migr_table =
    (migr_table_fmt *) mycalloc (1,
                                 sizeof (migr_table_fmt) *
                                 (*proposal)->old_migr_table_counter);
    (*proposal)->migr_table2 =
    (migr_table_fmt *) mycalloc ((*proposal)->old_migr_table_counter2,
                                 sizeof (migr_table_fmt));
    (*proposal)->migr_table_counter = 0;
    (*proposal)->migr_table_counter2 = 0;
#ifdef UEP
    if (world->options->uep)
    {
        (*proposal)->ueplike =
        (MYREAL **) mycalloc (world->data->uepsites, sizeof (MYREAL *));
        (*proposal)->ueplike[0] =
        (MYREAL *) mycalloc (world->numpop * world->data->uepsites,
                             sizeof (MYREAL));
        for (j = 1; j < world->data->uepsites; ++j)
            (*proposal)->ueplike[j] = (*proposal)->ueplike[0] + j * world->numpop;
        
        (*proposal)->ut.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        (*proposal)->uf.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        (*proposal)->umt = (MYREAL *) mycalloc (world->data->uepsites, sizeof (MYREAL));
        (*proposal)->umf = (MYREAL *) mycalloc (world->data->uepsites, sizeof (MYREAL));
    }
#endif
    
#ifdef BEAGLE
    (*proposal)->leftid   = 0;
    (*proposal)->rightid  = 0;
    (*proposal)->parentid = 0;
    reset_beagle(world->beagle);
#endif
}

///
/// sets new proposal structure using template gproposal, this function will be called most of the time instead of new_proposal()
void
set_new_proposal (proposal_fmt ** proposal, timelist_fmt * tl, world_fmt * world)
{
    long mal = world->data->maxalleles[world->locus];
    //long oldsize=0;
    long newsize =0;
    long sumtips = world->sumtips;
#ifdef UEP    
    long j;
#endif
    long listsize = 2*(world->sumtips + 2);
    (*proposal)->likelihood = -HUGE;
    (*proposal)->sumtips = sumtips;
    // pointers and values to outside structures
    (*proposal)->world = world;
    (*proposal)->datatype = world->options->datatype;
    (*proposal)->numpop = world->numpop;
    (*proposal)->endsite = world->data->seq[0]->endsite;
    (*proposal)->fracchange = world->data->seq[0]->fracchange;
    (*proposal)->param0 = world->param0;
    (*proposal)->root = world->root;
    (*proposal)->migration_model = world->options->migration_model;
    // precalculated values
    (*proposal)->mig0list = world->mig0list;
    (*proposal)->design0list = world->design0list;
    newsize = 4 * listsize;
    if((*proposal)->nodedata == NULL)
        (*proposal)->nodedata = (node **) calloc (newsize, sizeof (node *));
    else
        (*proposal)->nodedata = (node **) realloc ((*proposal)->nodedata, newsize * sizeof (node *));
    //xcode  oldsize = (*proposal)->listsize;
    memset((*proposal)->nodedata, 0, newsize * sizeof (node *));
    (*proposal)->aboveorigin = (*proposal)->nodedata;
    (*proposal)->bordernodes = (*proposal)->nodedata + listsize;
    (*proposal)->line_f =   (*proposal)->bordernodes + listsize;
    (*proposal)->line_t =  (*proposal)->line_f + listsize;
    (*proposal)->listsize =  listsize;
    // mf holds also mt array  
    (*proposal)->mf = (MYREAL *) mycalloc (2* (*proposal)->endsite, sizeof (MYREAL));
    (*proposal)->mt = (*proposal)->mf + (*proposal)->endsite;
    if (strchr (SEQUENCETYPES, (*proposal)->datatype))
    {
        allocate_xseq(&(*proposal)->xf, world->data->seq[0]->endsite, world->options->rcategs);
        allocate_xseq(&(*proposal)->xt, world->data->seq[0]->endsite, world->options->rcategs);
    }
    else
    {
        (*proposal)->xf.a = (MYREAL *) mycalloc (1, sizeof (MYREAL) * mal);
        (*proposal)->xt.a = (MYREAL *) mycalloc (1, sizeof (MYREAL) * mal);
    }
    (*proposal)->old_migr_table_counter = 4 * sumtips /* 100 */ ;
    (*proposal)->old_migr_table_counter2 = 4 * sumtips /* 100 */ ;
    (*proposal)->migr_table =
    (migr_table_fmt *) mycalloc (1,
                                 sizeof (migr_table_fmt) *
                                 (*proposal)->old_migr_table_counter);
    (*proposal)->migr_table2 =
    (migr_table_fmt *) mycalloc ((*proposal)->old_migr_table_counter2,
                                 sizeof (migr_table_fmt));
    (*proposal)->migr_table_counter = 0;
    (*proposal)->migr_table_counter2 = 0;
#ifdef UEP
    if (world->options->uep)
    {
        (*proposal)->ueplike =
        (MYREAL **) mycalloc (world->data->uepsites, sizeof (MYREAL *));
        (*proposal)->ueplike[0] =
        (MYREAL *) mycalloc (world->numpop * world->data->uepsites,
                             sizeof (MYREAL));
        for (j = 1; j < world->data->uepsites; ++j)
            (*proposal)->ueplike[j] = (*proposal)->ueplike[0] + j * world->numpop;
        
        (*proposal)->ut.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        (*proposal)->uf.s = (pair *) mycalloc (world->data->uepsites, sizeof (pair));
        (*proposal)->umt = (MYREAL *) mycalloc (world->data->uepsites, sizeof (MYREAL));
        (*proposal)->umf = (MYREAL *) mycalloc (world->data->uepsites, sizeof (MYREAL));
    }
#endif
    
#ifdef BEAGLE
    (*proposal)->leftid   = 0;
    (*proposal)->rightid  = 0;
    (*proposal)->parentid = 0;
    reset_beagle(world->beagle);
#endif
}


void
jumblenodes (node ** s, long n)
{
    node **temp;
    
    long i, rr, tn = n;
    
    temp = (node **) mycalloc (1, sizeof (node *) * n);
    memcpy (temp, s, sizeof (node *) * n);
    for (i = 0; i < n && tn > 0; i++)
    {
        s[i] = temp[rr = RANDINT (0, tn - 1)];
        temp[rr] = temp[tn - 1];
        tn--;
    }
    myfree(temp);
}

void traverse_check(node *theNode)
{
    if (theNode != NULL)
    {
        if (theNode->type != 't')
        {
            if (theNode->next->back != NULL)
            {
                if(theNode->tyme < showtop(theNode->next->back)->tyme)
                {
                    printf("Problem in traverse_check: id=%li time=%f type=%c time-up-next=%f type_up=%c\n",
                           theNode->id, theNode->tyme, theNode->type, theNode->next->back->tyme, theNode->next->back->type);
                    error("time conflict");
                }
                traverse_check (theNode->next->back);
            }
            if (theNode->type != 'm' && theNode->next->next->back != NULL)
            {
                if(theNode->tyme < theNode->next->next->back->tyme)
                {
                    printf("Problem in traverse_check: id=%li time=%f type=%c time-up-next-next=%f\n",
                           theNode->id, theNode->tyme, theNode->type, theNode->next->next->back->tyme);
                    error("time conflict");
                }
                
                traverse_check (theNode->next->next->back);
            }
        }
    }
    else
    {
        error("Node is NULL????");
    }
}

///
/// tags the lineage that was added the last time the tree was changed.
void traverse_tagnew(node *theNode, node *origin)
{
    if (theNode != NULL)
    {
        theNode->visited = FALSE;
        if (theNode->type != 't'  && theNode->next->next->back != NULL)
        {
            if(theNode->top != 1)
                theNode = showtop(theNode);
            //error("");
            
            traverse_tagnew (theNode->next->back, origin);       
            if (theNode->type != 'm' && theNode->next->next->back != NULL)
            {
                traverse_tagnew (theNode->next->next->back, origin);
            }
        }
        if(theNode==origin)
        {
            while(theNode != NULL && theNode->type != 'r')
            {
                theNode = showtop(showtop(theNode)->back);
                theNode->visited = TRUE;
            }
        }
    }
    else
    {
        error("a missing node pointer encountered, aborted");
    }
}

///
/// pick a node at random from the available list of nodes, this code will even 
/// work with datasets where there are only two nodes.
void
chooseOrigin (proposal_fmt * proposal)
{
    long elem = 0;
    // oldelem = (proposal->sumtips * 2.);
    node *tmp=NULL, **goal;
    //032110 goal = (node **) mycalloc (oldelem, sizeof (node *));
    //#ifdef DEBUG
    //traverse_check(crawlback (proposal->root->next));
    //#endif
    //032110 traverseAllNodes (crawlback (proposal->root->next), &goal, &elem, &oldelem,
    //032110                  NO_MIGR_NODES);
    goal = proposal->world->nodep;
    elem = proposal->world->sumtips * 2;
    switch(elem)
    {
        case 0:
            error("problem with choosing a node for the tree change [choosOrigin()]");
            break;
        case 1:
        case 2:
            tmp = goal[0];
            break;
        default:
            tmp = goal[RANDINT (0, elem - 2)];
            while(tmp->back->type =='r')
                tmp = goal[RANDINT (0, elem - 2)];
            break;
    }
    //032110 myfree(goal);
    proposal->origin = tmp;
    if (proposal->origin != showtop (crawlback (proposal->root->next)))
    {
        proposal->oback = showtop (crawlback (proposal->origin));
        proposal->osister = showsister (proposal->origin);
        if (proposal->oback != showtop (crawlback (proposal->root->next)))
        {
            proposal->ocousin = showsister (proposal->oback);
        }
        else
        {
            proposal->ocousin = NULL;
        }
    }
    if (proposal->origin == NULL)
        error ("Designation of origin for branch removal failed");
}


///
/// construct a list of ordered times with pointers to the tree
/// using the origin to prune the old time list
void
construct_localtimelist (timelist_fmt * timevector, proposal_fmt * proposal)
{
#ifdef TREEDEBUG
    double dif;
    long ii;
#endif
    long z = 0;
    long oz = proposal->listsize;
    //    long numpop = timevector->numpop = proposal->numpop;
    traverseAllNodes (crawlback (proposal->origin)->back,
                      &proposal->aboveorigin, &z, &oz, WITH_MIGR_NODES);
    proposal->aboveorigin[z++] = proposal->oback;
    proposal->aboveorigin[z] = NULL;
    prune_timelist(&proposal->world->treetimes[0],timevector, proposal->aboveorigin, proposal, proposal->world->numpop);
    add_partlineages(proposal->numpop, &timevector);
#ifdef TREEDEBUG
    for(ii=0;ii<timevector->T-2;ii++)
    {
        printf("%i> ii=%li %0.5f  (%li->%li) %c %3li",myID, ii, timevector->tl[ii].age, 
               timevector->tl[ii].from, timevector->tl[ii].to, 
               timevector->tl[ii].eventnode->type,timevector->tl[ii].lineages[0]);
        long jj;
        for(jj = 1;jj < proposal->numpop; jj++)
            printf(" %3li",timevector->tl[ii].lineages[jj]);
        printf(" backtyme:%f dif:%f\n",showtop(timevector->tl[ii].eventnode->back)->tyme,dif=showtop(timevector->tl[ii].eventnode->back)->tyme-timevector->tl[ii].eventnode->tyme);
        if(dif<0.0)
        {
            error("shit");
        }
    }
#endif
}

/*----------------------------------------------------------------------------
 finds all nodes in a tree starting at the root node and crawling up 
 to the tips in a recursive fashion, writing nodeptrs in the nodelist vector
 the flag include_migration is 1 if we want to touch the migration nodes too,
 otherwise =0 -> jump over the migration nodes. for convenience we define the 
 the macros NO_MIGR_NODES=0 and WITH_MIGR_NODES=1 in the treesetup.h file
 PB 1995
 */
void
traverseAllNodes (node * theNode, node *** nodelist, long *node_elem,
                  long *oldnode_elem, int include_migration)
{
    static long counter=0;
    long elem;
    counter++;
    if(theNode->type != 'r' && theNode->tyme > showtop(theNode->back)->tyme)
    {
        printf("%i> TTRAVERSETRAVERSETRAVERSETRAVERSE\nfunction was called with %li times: with node (%c) time %f and node (%c) time %f\n",
               myID, counter, theNode->type,  theNode->tyme, showtop(theNode->back)->type, showtop(theNode->back)->tyme);
        error("traverseAllNodes() failed");
    }
    
    if (include_migration == NO_MIGR_NODES)
    {
        // nodelist needs to be at least twice as long as sumtips
        // otherwise it will run out of reserved memory 
        if (theNode->type != 't')
        {
            if (crawlback (theNode->next) != NULL)
                traverseAllNodes (crawlback (theNode->next), nodelist, node_elem,
                                  oldnode_elem, NO_MIGR_NODES);
            if (theNode->type != 'm' && crawlback (theNode->next->next) != NULL)
                traverseAllNodes (crawlback (theNode->next->next), nodelist,
                                  node_elem, oldnode_elem, NO_MIGR_NODES);
            (*nodelist)[(*node_elem)] = theNode;
            (*node_elem) += 1;
            if (theNode->type == 'm')
            {
                error ("Migration node encountered?! and died!");
            }
            if(*node_elem > *oldnode_elem)
                error("die on the spot, this should not happen. Failed in traverseAllNodes without migration nodes");
        }
        else
        {
            (*nodelist)[(*node_elem)] = theNode;
            (*node_elem) += 1;
            if(*node_elem > *oldnode_elem)
                error("die on the spot, this should not happen. Failed in traverseAllNodes without migration nodes at tip");
        }
    }
    else
    {
        if (theNode->type != 't')
        {
            if (theNode->next->back != NULL)
                traverseAllNodes (theNode->next->back, nodelist, node_elem,
                                  oldnode_elem, WITH_MIGR_NODES);
            if (theNode->type != 'm' && theNode->next->next->back != NULL)
                traverseAllNodes (theNode->next->next->back, nodelist, node_elem,
                                  oldnode_elem, WITH_MIGR_NODES);
            if ((*node_elem) == (*oldnode_elem - 2))
            {
                elem = (*oldnode_elem) + (*oldnode_elem);
                (*nodelist) =
                (node **) myrealloc ((*nodelist), sizeof (node *) * elem);
                memset ((*nodelist) + (*oldnode_elem), 0,
                        sizeof (node *) * (elem - (*oldnode_elem)));
                *oldnode_elem = elem;
            }
            (*nodelist)[(*node_elem)++] = theNode;
        }
        else
        {
            if ((*node_elem) == (*oldnode_elem - 2))
            {
                elem = (*oldnode_elem) + (*oldnode_elem);
                (*nodelist) =
                (node **) myrealloc ((*nodelist), sizeof (node *) * elem);
                memset ((*nodelist) + (*oldnode_elem), 0,
                        sizeof (node *) * (elem - (*oldnode_elem)));
                *oldnode_elem = elem;
            }
            (*nodelist)[(*node_elem)++] = theNode;
        }
    }
}

///
/// copies elements from the old timelist (oldtv)
/// into the new local timeslist (newtv) checking whether they
/// are used after the removal of the residual tree that is above the origin (ptr)
void
prune_timelist (timelist_fmt * oldtv, timelist_fmt * newtv,  register node ** __restrict ptr, proposal_fmt * proposal, long numpop)
{
    register long i    = 0;
    register long j    = 0;
    register node * thenode;
    long          slot = 0;
    vtlist        * tls;
    
    slot=0;
    // go through all slices in old time list and
    // extract node pointers
    for (i = 0; i < (*oldtv).T; i++)
    {
        j = 0;
        thenode = (*oldtv).tl[i].eventnode;
        while ((thenode != ptr[j]) && (ptr[j] != NULL))
        {
            j++;
        }
        // if the comparison reaches the end of the ptr list that holds an NULL element at the end
        // then the node needs to be present in the new time list.
        if (ptr[j] == NULL)
        {
            tls = &((*newtv).tl[slot]);
            tls->eventnode = thenode ;
            tls->age       = (*oldtv).tl[i].age ;
            tls->interval  = (*oldtv).tl[i].interval ;
            if(thenode!=NULL)
            {
                tls->from      = thenode->pop;
                tls->to        = thenode->actualpop;
            }
            tls->slice     = slot;
            slot++;
        }
    }
    (*newtv).T = slot;
    //   if ((*newtv).tl[(*newtv).T - 1].eventnode->type != 'r')
    //{
    //    error ("Root not at the end of local timelist\n");
    //}
    // sort seems to be necessary with dated samples
    qsort ((void *) (*newtv).tl, (*newtv).T, sizeof (vtlist), agecmp);
}

/* replaces nodepointers in list 1 with NULL if they are present in list 2
 returns the first NULL slot in the array.
 */
long
xor (node ** ptrl1, node ** ptrl2)
{
    long i = 0, j = 0, slot = -1;
    /* assumes that there is an NULL element at the end */
    for (i = 0; ptrl1[i] != NULL; j = 0, i++)
    {
        while ((ptrl1[i] != ptrl2[j]) && (ptrl2[j] != NULL))
            j++;
        if (ptrl2[j] != NULL)
        {
            if (slot == -1)
                slot = i;
            ptrl1[i] = NULL;
        }
    }
    return slot;
}

/* migrate() fills the PROPOSAL->MIGRATION_TABLE in the tree or
 PROPOSAL->MIGRATION_TABLE2 when at the bottom of the tree
 
 PROPOSAL proposal-scratchpad
 UP       node above (younger)
 MIGR_TABLE_COUNTER migration array counter, 
 will increase by one during execution
 AIR      if true standard execution, if false updating the last 
 lineage in the residual tree.
 */
int
migrate (proposal_fmt * proposal, node * up)
{
    long tmp;
    long numpop = proposal->numpop;
    long i = proposal->migr_table_counter;
    migr_table_fmt *array = proposal->migr_table;
    if (i > MIGRATION_LIMIT * numpop || numpop < 2)
    {
        return 0;
    }
    if (i > 0)
        array[i].to = array[i - 1].from;
    else
        array[i].to = up->pop;
    tmp = migration_from (array[i].to, proposal);
    
    if(tmp < numpop)
        array[i].from = tmp;
    else
        return 0;
    
    //DEBUG
    // printf("i: %3li %li %li\n",i,proposal->migr_table[i].from,proposal->migr_table[i].to);
    array[i++].time = proposal->time;
    if(proposal->time < proposal->origin->tyme)
    {
        error("in migrate() wrong time found");
    }
    if (i >= proposal->old_migr_table_counter)
    {
        proposal->old_migr_table_counter += 10;
        proposal->migr_table =
        (migr_table_fmt *) myrealloc (proposal->migr_table,
                                      sizeof (migr_table_fmt) *
                                      (proposal->old_migr_table_counter));
    }
    proposal->migr_table_counter = i;
    return 1;
}

int
migrateb (proposal_fmt * proposal, node * up)
{
    long i = proposal->migr_table_counter2;
    migr_table_fmt *array = proposal->migr_table2;
    if (i > MIGRATION_LIMIT * proposal->numpop)
    {
        return 0;
    }
    if (i > 0)
        array[i].to = array[i - 1].from;
    else
        array[i].to = up->pop;
    array[i].from = migration_from (array[i].to, proposal);
    //  printf("t: %3li %li %li\n",i,proposal->migr_table2[i].from,proposal->migr_table2[i].to);
    array[i].time = proposal->time;
    i++;
    if (i >= proposal->old_migr_table_counter2)
    {
        proposal->old_migr_table_counter2 += 10;
        proposal->migr_table2 =
        (migr_table_fmt *) myrealloc (proposal->migr_table2,
                                      sizeof (migr_table_fmt) *
                                      (proposal->old_migr_table_counter2));
    }
    proposal->migr_table_counter2 = i;
    return 1;
}

int
migrate_old (proposal_fmt * proposal, node * up, long *old_migr_table_counter,
             boolean air)
{
    migr_table_fmt *array;
    long i;
    if (air)
    {
        array = proposal->migr_table;
        i = proposal->migr_table_counter;
    }
    else
    {
        array = proposal->migr_table2;
        i = proposal->migr_table_counter2;
    }
    if (i > MIGRATION_LIMIT * proposal->numpop)
    {
        //      FPRINTF (stdout, "migration limit reached\n");
        return 0;
    }
    switch (proposal->migration_model)
    {
        case ISLAND:
        case ISLAND_VARTHETA:
        case MATRIX:
        case MATRIX_SAMETHETA:
        case MATRIX_ARBITRARY:
            if (i > 0)
                array[i].to = array[i - 1].from;
            else
                array[i].to = up->pop;
            array[i].from = migration_from (array[i].to, proposal);
            //      printf("i: %3li %li %li\n",i,proposal->migr_table[i].from,proposal->migr_table[i].to);
            //      printf("b: %3li %li %li\n",i,proposal->migr_table2[i].from,proposal->migr_table2[i].to);
            break;
        case CONTINUUM:
        case STEPSTONE:
            error ("not yet implemented\n");
            break;
        default:
            break;
    }
    array[i++].time = proposal->time;
    if (i >= (*old_migr_table_counter))
    {
        (*old_migr_table_counter) += 10;
        if (air)
        {
            proposal->migr_table =
            (migr_table_fmt *) myrealloc (proposal->migr_table,
                                          sizeof (migr_table_fmt) *
                                          (*old_migr_table_counter));
            //xcode  array = proposal->migr_table;
        }
        else
        {
            proposal->migr_table2 =
            (migr_table_fmt *) myrealloc (proposal->migr_table2,
                                          sizeof (migr_table_fmt) *
                                          (*old_migr_table_counter));
            //xcode array = proposal->migr_table;
        }
    }
    if (air)
    {
        proposal->migr_table_counter = i;
    }
    else
    {
        proposal->migr_table_counter2 = i;
    }
    return 1;
}

/* migration_from() returns the FROM population when there was a migration
 TO        population to migrate to
 PROPOSAL  proposal-scratchpad
 */
long
migration_from_old (long to, proposal_fmt * proposal)
{
    long j, ii, msta, msto;
    MYREAL *geo = proposal->world->data->geo;
    MYREAL *r, rr = UNIF_RANDUM ();
    r = (MYREAL *) mycalloc (1, sizeof (MYREAL) * proposal->numpop);
    msta = mstart (to, proposal->numpop);
    msto = mend (to, proposal->numpop);
    r[0] = proposal->param0[msta] * geo[msta];
    for (j = 1, ii = msta + 1; ii < msto; j++, ii++)
    {
        r[j] = r[j - 1] + geo[ii] * proposal->param0[ii];
    }
    ii = 0;
    while (rr > r[ii] / r[j - 1])
    {
        ii++;
    }
    myfree(r);
    if (ii < to)
        return ii;
    else
        return ++ii;
}

long
migration_from (long to, proposal_fmt * proposal)
{
    long ii = 0;
    MYREAL *r = proposal->world->migproblist[to];
    MYREAL rr = UNIF_RANDUM ();
    while (rr >  r[ii] && ii < proposal->world->numpop-1)
    {
        ii++;
    }
#ifdef TREEDEBUG
    if(r[ii] == 0.0)
    {
        warning("DEBUG: migproblist does not contain a migration probability\n");
    }
#endif
    if (ii < to)
        return ii;
    else
        return ++ii;
}

void
chooseTarget (proposal_fmt * proposal, timelist_fmt * timevector,
              node ** bordernodes, long *bordernum)
{
    long actualpop = -99;
    node *rb = crawlback (proposal->root->next);
    *bordernum = 0;
    proposal->target = NULL;
    proposal->realtarget = NULL;
    if (proposal->migr_table_counter == 0)
        actualpop = proposal->origin->pop;
    else
        actualpop = proposal->migr_table[proposal->migr_table_counter - 1].from;
    if (rb->tyme < proposal->time)
    {
        error ("Wrong Time for action in chooseTarget()\n");
    }
    findbordernodes (rb, proposal, actualpop, &bordernodes, &proposal->listsize, bordernum,
                     &(*timevector).tl, (*timevector).T);
    if (*bordernum > 0)
    {
        // found elegible linages to coalesce to
        proposal->target = bordernodes[RANDINT (0, (*bordernum) - 1)];
        if (proposal->target != rb)
        {
            proposal->tsister = showsister (proposal->target);
            proposal->realtsister = crawlback (proposal->tsister)->back;
        }
        else
            proposal->tsister = NULL;
        proposal->realtarget = proposal->target;
        if (proposal->target->type == 'm')
            proposal->target = crawlback (showtop (proposal->target)->next);
    }
    else
    {
        // no lineage to coalesce to was found.
        proposal->target = NULL;
        proposal->tsister = NULL;
        proposal->realtsister = NULL;
        proposal->realtarget = NULL;
    }
}

void
findbordernodes (node * theNode, proposal_fmt * proposal, long pop,
                 node *** bordernodes, long *allocsize, long *bordernum, vtlist ** tyme,
                 long gte)
{
    node *tmp, *back;
    // we search on the old tree and reaching the oback node that was excised from
    // from the residual tree we jump over it and go up the sister branch and back is
    // going further down than the oback node
    if (theNode == proposal->oback)
    {
        tmp = showtop (crawlback (proposal->osister)->back);
        back = showtop (proposal->oback->back);
    }
    else
    {
        tmp = showtop (theNode);
        back = showtop (theNode->back);
    }
    if (pop == tmp->pop && pop == back->actualpop && tmp->tyme < proposal->time
        && back->tyme > proposal->time)
    {
        // lineage end points bracket the the new proposed time,
        // therefore this branch needs to be included
        (*bordernodes)[(*bordernum)++] = tmp;
        return;
    }
    else
    {
        // the lineage does not fit the populations (it may fit the time, so)
        // the if here checks this
        if (back->tyme < proposal->time)
        {
            // this makes sure that we looked at all lineages that
            // may be electable to receive a coalescent have been looked
            // at, if we reach here the time of the proposal is older than
            // all available lineages on this branch and we safely can stop
            // following this branch.
            return;
        }
        // if we are still working on interior nodes and the proposal time is still 
        // younger than then lineage bracketing nodes we need to continue our search
        // up in the tree following right and left branches, but need to stop when 
        // we encounter a tip node (this is important with dated tips because these can
        // be interspersed in the timelist.
        if (tmp->type != 't')
        {
            if (tmp->next->back != NULL)
                findbordernodes (tmp->next->back, proposal, pop, bordernodes, allocsize,
                                 bordernum, tyme, gte);
            if (tmp->type != 'm' && tmp->next->next->back != NULL)
                findbordernodes (tmp->next->next->back, proposal, pop,
                                 bordernodes, allocsize, bordernum, tyme, gte);
        }
    }
}

/*
 boolean
 same_pop(node * up, MYREAL tyme, long pop)
 {
 node *oldnn = showtop(up->back);
 node *nn = up;
 while (nn->tyme < tyme) {
 oldnn = nn;
 nn = showtop(nn->back);
 }
 if (oldnn->pop == pop && nn->actualpop == pop)
 return TRUE;
 else
 return FALSE;
 }
 */


/* -----------------------------------------------------------------------
 simulates two lineages at once, if we are extending below the last node */
int
pre_population (proposal_fmt * proposal, vtlist * ltime, long gte,
                long *slider)
{
    boolean coalesced = FALSE;
    boolean choice = FALSE;
    long pop1 = -99, pop2 = -98;
    //  long msta1=0, msto1=0;
    //  long msta2=0, msto2=0;
    MYREAL ux;
    MYREAL  rate = proposal->world->options->mu_rates[proposal->world->locus];
    MYREAL  invrate = 1./rate;
    MYREAL age1, denom, rr, r0, r1, horizon, mm, mm2;
    //    MYREAL denom_invrate;
    if (gte > 0)
        proposal->realtarget = ltime[gte - 1].eventnode;
    else
        proposal->realtarget = ltime[0].eventnode->next->back; //?????gte
    if (proposal->realtarget == proposal->oback)
    {
        proposal->realtarget = crawlback (proposal->osister)->back;
    }
    if (proposal->realtarget->type == 'm')
    {
        proposal->target = crawlback (proposal->realtarget->next);
        if (proposal->target == proposal->oback)
        {
            proposal->target = proposal->osister;
        }
    }
    else
    {
        proposal->target = proposal->realtarget;
    }
    proposal->tsister = NULL;
    pop2 = proposal->realtarget->pop;
    pop1 =
    proposal->migr_table_counter >
    0 ? proposal->migr_table[proposal->migr_table_counter -
                             1].from : proposal->origin->pop;
    age1 =
    MAX (proposal->realtarget->tyme,
         proposal->migr_table_counter >
         0 ? proposal->migr_table[proposal->migr_table_counter -
                                  1].time : proposal->origin->tyme);
    horizon = MAX (proposal->oback->tyme, age1);
    while (age1 < horizon)
    {
        mm = proposal->mig0list[pop1] * invrate;
        //mm = 0.0;
        //msta1 = mstart(pop1,proposal->numpop);
        //msto1 = mend(pop1,proposal->numpop);
        //      for (i = msta1; i < msto1; i++)
        //{
        //mm += proposal->param0[i];
        //}
        if (pop1 == pop2)
        {
            denom = mm + (2. / (proposal->param0[pop1] * rate));
            //	    denom_invrate = denom * invrate;
            ux = UNIF_RANDUM();
            //            proposal->time = age1 - LOG(ux) / denom_invrate;
            proposal->time = age1 - LOG(ux) / denom;
            age1 = proposal->time;
            if (age1 < horizon)
            {
                rr = UNIF_RANDUM ();
                r0 = (2. / (proposal->param0[pop1] * rate)) / denom;
                if (rr < r0)
                {
                    return 1;
                }
            }
        }
        else
        {
            denom = mm;
            proposal->time = age1 - LOG (UNIF_RANDUM ()) / denom;
            age1 = proposal->time;
        }
        if (age1 < horizon)
        {
            if (!migrate (proposal, proposal->origin))
                //                      &proposal->old_migr_table_counter, MIGRATION_AIR))
            {
                return 0;
            }
            pop1 =
            proposal->migr_table_counter >
            0 ? proposal->migr_table[proposal->migr_table_counter -
                                     1].from : proposal->origin->pop;
        }
    }
    age1 = horizon;
    while (!coalesced)
    {
        //mm = mm2 = 0;
        //msta1 = mstart(pop1,proposal->numpop);
        //msto1 = mend(pop1,proposal->numpop);
        //msta2 = mstart(pop2,proposal->numpop);
        //msto2 = mend(pop2,proposal->numpop);
        //for (i = msta1; i < msto1; i++)
        //{
        //mm += proposal->param0[i];
        //}
        
        //limits treesize
        if(proposal->time > 1000000)
            return 0;
        
        mm = proposal->mig0list[pop1] * invrate;
        mm2 = proposal->mig0list[pop2] * invrate;
        
        if (pop1 == pop2)
        {
            denom = 2. * mm + (2. / proposal->param0[pop1] * rate);
            //	    denom_invrate = denom * invrate;
            proposal->time = age1 - LOG (UNIF_RANDUM ()) / denom;
            age1 = proposal->time;
            rr = UNIF_RANDUM ();
            r0 = ((2. / proposal->param0[pop1] * rate) / denom);
            r1 = r0 + mm / denom;
            if (rr < r0)
            {
                return 1;
            }
            else
            {
                if (rr < r1)
                {
                    choice = TRUE;
                }
                else
                {
                    choice = FALSE;
                }
            }
        }
        else
        {   /*pop1 not equal pop2 */
            //for (i = msta2; i < msta2; i++)
            //{
            //mm2 += proposal->param0[i];
            //}
            denom = mm + mm2;
            //	    denom_invrate = denom * invrate;
            proposal->time = age1 - LOG (UNIF_RANDUM ()) / denom;
            age1 = proposal->time;
            if (RANDUM () < (mm / denom))
            {
                choice = TRUE;
            }
            else
            {
                choice = FALSE;
            }
        }
        if (choice)
        {
            if (!migrate (proposal, proposal->origin))
                //                      &proposal->old_migr_table_counter, MIGRATION_AIR))
            {
                return 0;  /* migration limit reached */
            }
            pop1 =
            proposal->migr_table_counter >
            0 ? proposal->migr_table[proposal->migr_table_counter -
                                     1].from : proposal->origin->pop;
        }
        else
        {
            if (!migrateb (proposal, proposal->realtarget))
                //                      &proposal->old_migr_table_counter2,
                //                      MIGRATION_IN_TREE))
            {
                return 0;  /* migration limit reached */
            }
            pop2 =
            proposal->migr_table_counter2 >
            0 ? proposal->migr_table2[proposal->migr_table_counter2 -
                                      1].from : proposal->realtarget->pop;
        }
    }
    error ("Reached the end of function without coalescing");
    return -1;   /*makes the compiler happy */
}

#ifdef TESTING2
///
/// freeing the proposal structure, this is the replacement of the real free_proposal method
/// and does not allocate, but reuses old memory
void free_proposal(proposal_fmt *proposal)
{
    // do nothing
}
void reset_simple_proposal_variables(proposal_fmt **proposal)
{
    (*proposal)->mig_removed = FALSE;
    (*proposal)->rr = 0.0;
    (*proposal)->origin = NULL;
    (*proposal)->target = NULL;
    (*proposal)->realtarget = NULL;
    (*proposal)->tsister = NULL;
    (*proposal)->realtsister = NULL;
    (*proposal)->osister = NULL;
    (*proposal)->realosister = NULL;
    (*proposal)->ocousin = NULL;
    (*proposal)->realocousin = NULL;
    (*proposal)->oback = NULL;
    (*proposal)->realoback = NULL;
    //
    (*proposal)->connect = NULL;
    (*proposal)->likelihood = 0.0;
    (*proposal)->time = 0.0;
    (*proposal)->v = 0.0;
    (*proposal)->vs = 0.0;
    //
#ifdef UEP
    (*proposal)->ueplikelihood = 0.0;
#endif
    (*proposal)->migr_table_counter=0;
    (*proposal)->migr_table_counter2=0;
    (*proposal)->timeslice = 0;
    //
    (*proposal)->treelen = 0.0;
#ifdef BEAGLE
    (*proposal)->parentid = 0;
    (*proposal)->leftid = 0;
    (*proposal)->rightid = 0;
#endif
}

void
reset_proposal (proposal_fmt ** proposal, world_fmt *world)
{
    const long listsize = 2 * (world->sumtips + 2);
    const long newsize = 4 * listsize;
    const long mal = (*proposal)->world->data->maxalleles[(*proposal)->world->locus];
    (*proposal)->likelihood = -HUGE;
    // pointers and values to outside structures
    (*proposal)->world = world;
    (*proposal)->datatype = world->options->datatype;
    (*proposal)->sumtips = world->sumtips;
    (*proposal)->numpop = world->numpop;
    (*proposal)->endsite = world->data->seq[0]->endsite;
    (*proposal)->fracchange = world->data->seq[0]->fracchange;
    (*proposal)->param0 = world->param0;
    (*proposal)->root = world->root;
    (*proposal)->migration_model = world->options->migration_model;
    // precalculated values
    (*proposal)->mig0list = world->mig0list;
    (*proposal)->design0list = world->design0list;
    (*proposal)->listsize =  listsize;
    // ..line_t are also reset with this
    long i;
    for(i=0;i<newsize;i++)
        (*proposal)->nodedata[i]=NULL;
    //  memset((*proposal)->nodedata,0,sizeof(node *) * newsize); 
    memset((*proposal)->mf,0, sizeof(MYREAL)*(2 * (*proposal)->endsite)); // (*proposal)->mt is also freed with this
    // resetting pointers
    reset_simple_proposal_variables(proposal);
    
    if (strchr (SEQUENCETYPES, (*proposal)->datatype))
    {
        zero_xseq(&(*proposal)->xf,(*proposal)->world->data->seq[0]->endsite, (*proposal)->world->options->rcategs);
        zero_xseq(&(*proposal)->xt,(*proposal)->world->data->seq[0]->endsite, (*proposal)->world->options->rcategs);
    }
    else
    {
        memset((*proposal)->xf.a, 0, mal);
        memset((*proposal)->xt.a, 0, mal);
    }
    memset((*proposal)->migr_table, 0, sizeof(migr_table_fmt) * (*proposal)->old_migr_table_counter);
    memset((*proposal)->migr_table2, 0, sizeof(migr_table_fmt) * (*proposal)->old_migr_table_counter2);
#ifdef UEP
    
    if ((*proposal)->world->options->uep)
    {
        memset((*proposal)->ueplike, 0, sizeof(MYREAL) * world->numpop * world->data->uepsites);
        for (j = 1; j < world->data->uepsites; ++j)
            (*(*proposal))->ueplike[j] = (*(*proposal))->ueplike[0] + j * world->numpop;
        
        memset((*proposal)->uf.s, 0, world->data->uepsites * sizeof(pair));
        memset((*proposal)->ut.s, 0, world->data->uepsites * sizeof(pair));
        memset((*proposal)->umf, 0, world->data->uepsites * sizeof(MYREAL));
        memset((*proposal)->umt, 0, world->data->uepsites * sizeof(MYREAL));
    }
#endif
}
#else
///
/// freeing the proposal structure, this will be replaced by a function that resets permanently
/// allocated structure [reset_proposal()]
void
free_proposal (proposal_fmt * proposal)
{
    //static long count=0;
    myfree(proposal->aboveorigin);
    // myfree(proposal->bordernodes);
    // myfree(proposal->line_f);
    // myfree(proposal->line_t);
    //  printf("%li ",count++);fflush(stdout);
    myfree(proposal->nodedata); // ..line_t are also freed with this
    myfree(proposal->mf); // proposal->mt is also freed with this
    if (strchr (SEQUENCETYPES, proposal->datatype))
    {
        myfree(proposal->xf.s[0]);
        myfree(proposal->xt.s[0]);
        myfree(proposal->xf.s);
        myfree(proposal->xt.s);
    }
    else
    {
        myfree(proposal->xf.a);
        myfree(proposal->xt.a);
    }
    myfree(proposal->migr_table);
    myfree(proposal->migr_table2);
#ifdef UEP
    
    if (proposal->world->options->uep)
    {
        myfree(proposal->ueplike[0]);
        myfree(proposal->ueplike);
        myfree(proposal->uf.s);
        myfree(proposal->ut.s);
        myfree(proposal->umf);
        myfree(proposal->umt);
        
    }
#endif
    myfree(proposal);
}
#endif /*TESTING*/


void
free_timevector (timelist_fmt * timevector)
{
    //    long i;
    //    for (i = 0; i < timevector->allocT; i++)
    //    {
    //        myfree(timevector->tl[i].lineages);
    //    }
    myfree(timevector->lineages);
    myfree(timevector->tl);
    myfree(timevector);
}

///
/// using a global timevector to save malloc/free pair calls and 
/// so to save time
void
free_timevector_new (timelist_fmt * timevector)
{
    //    memset(timevector->lineages, 0, sizeof(long) * );
    // (timevector->tl);
}

/*----------------------------------------------------------*
 * rejection/acceptance of the new tree according to the likelihood
 * and an acceptance ratio which is higher the better the
 * likelihood values are (-> Metropolis)
 */
boolean
acceptlike (world_fmt * world, proposal_fmt * proposal, long g,
            timelist_fmt * tyme)
{
    //proposal changes when chain is stuck, to accept some new trees without
    //considering data to remove out of sticky area [hopefully]
    static long not_accepted = 0;
    
    const long limit = MIGRATION_LIMIT * world->numpop;
    
    long rm  = 0L;
    long rmc = rmigrcount (proposal);
    
    MYREAL rr;
    MYREAL expo;
    
#ifdef UEP
    node *first;
#endif
    
    rm =  proposal->migr_table_counter + proposal->migr_table_counter2 
    + world->migration_counts - rmc;
    
    if (rm > limit)
    {
        // warning might disturb parallel runs and might be not really useful at all
        //        if (report || g != oldg)
        //        {
        //            warning ("migration limit (%li) exceeded: %li\n",
        //                     MIGRATION_LIMIT * world->numpop, rm);
        //            warning
        //            ("results may be underestimating migration rates for this chain\n");
        //            report = FALSE;
        //            oldg = g;
        //        }
#ifdef BEAGLE
        printf("Reset to: %f\n",world->likelihood[g]);
        //set_beagle_dirty(proposal->origin,proposal->target,showtop(world->root->next->back));
        //calcLnL(world, world->beagle->instance);
#endif
        return FALSE;
    }
#ifdef UEP
    if (world->options->uep)
    {
        first = first_uep2 (proposal->world->root->next->back,
                            proposal->world->root->next->back,
                            proposal->world->data->uepsites);
        proposal->firstuep = first_uep (first, proposal->world->root,
                                        proposal->world->data->uepsites);
        proposal->ueplikelihood = pseudo_ueplikelihood (world, proposal);
        proposal->likelihood = pseudotreelikelihood (world, proposal);
    }
    else
    {
        proposal->likelihood = pseudotreelikelihood (world, proposal);
    }
#else
    //    printf("DEBUG: world->likellihood[%li]=%f heat=%f\n",g,world->likelihood[g],world->heat);
    proposal->likelihood = pseudotreelikelihood (world, proposal);
    //printf("DEBUG: proposal->likellihood=%f heat=%f\n",proposal->likelihood, world->heat);
#endif
    
    //if(proposal->likelihood <= -HUGE &&  world->likelihood[g] <= -HUGE)
    //  {
	//	if(!world->in_burnin)
	//  {
    //long burnin = world->options->burn_in;
    //world->options->burn_in *= 10;
    //fprintf(stdout,"likelihood is -inf, insert a burn-in phase\n");
    //	    burnin_chain(world);
    //world->options->burn_in = burnin;
    //  }
    //	return TRUE;
    // }
    
    if (world->likelihood[g] < proposal->likelihood)
    {
        not_accepted = 0;
        return TRUE;
    }
    if (!world->options->heating)
    {
        expo = proposal->likelihood - world->likelihood[g];
        rr = LOG (RANDUM ());
        if (rr < expo)
        {
            not_accepted = 0;
            return TRUE;
        }
    }
    else
    {
        expo = (proposal->likelihood - world->likelihood[g]) * world->heat;
        rr = LOG (RANDUM ());
        if (rr < expo)
        {
            not_accepted = 0;
            return TRUE;
        }
    }
    if(world->options->prioralone)
    {
        return TRUE;
    }
    return FALSE;
}

long
rmigrcount (proposal_fmt * proposal)
{
    node *p;
    long count = 0;
    for (p = proposal->origin; p != proposal->oback; p = showtop (p->back))
    {
        if (p->type == 'm')
            count++;
    }
    return count;
}


///
/// generates the time interval to the next event (migration or coalescence)
/// and also decides what event happens.
/// This funcion is modified by the rate of the mutation rate that can be
/// estimated in a Bayesian context, and be fixed in a ML context
/// this rate is only influencing the time and not the eventtype
/// because the rate cancels out of the ratio that chooses the event.
MYREAL
eventtime (proposal_fmt * proposal, long pop, vtlist * tentry, char *event)
{
    //    static boolean mig_force = TRUE;
    MYREAL  interval;
    MYREAL  lines;
    MYREAL  denom;
    MYREAL  invdenom;
    MYREAL  rate = proposal->world->options->mu_rates[proposal->world->locus];
    MYREAL  inheritance = proposal->world->options->inheritance_scalars[proposal->world->locus];
    MYREAL  invrate = 1./ rate;
    MYREAL  mm         = proposal->mig0list[pop] * invrate ;
    //long    design0    = proposal->design0list[pop];
    //could be used to calculate p(gn|go): MYREAL nu;
    lines    = 2.0 * tentry->lineages[pop];
    denom    = mm + (lines * (1.0 / (inheritance*proposal->param0[pop]*rate)));
    invdenom = 1.0 / denom;
    interval =  (-(LOG (/*nu = */UNIF_RANDUM ())) * invdenom) ;
    //could be used to calculate p(gn|go): proposal->nu += nu; //adds the probabilities (=random numbers);
    //[overcautious]
    //if (interval < 0.0)
    //    error("interval time is negative");
    if (lines > 0.0)
    {
        if ((UNIF_RANDUM ()) < (mm * invdenom))
        {
            *event = 'm';
            return interval;
        }
        else
        {
            *event = 'c';
            return interval;
        }
    }
    else
    {
        //      printf("mcmc1.c 1653 pop = %li, lines = %li\n", pop, tentry->lineages[pop]);
        *event = 'm';
        return interval;
    }
}

/*--------------------------------------------------------*
 * showsister() 
 * find the sisternode, by going down the branch and up on 
 * the other side again, neglecting the migration nodes.
 */
node *
showsister (node * theNode)
{
    node *tmp = crawlback (theNode);
    
    if (tmp->next->top)
    {
        return crawlback (tmp->next->next);
    }
    else
    {
        if (tmp->next->next->top)
        {
            return crawlback (tmp->next);
        }
        else
        {
            error ("error in treestructure, cannot find sisternode\n");
        }
    }
    return NULL;
}

void
count_migrations (node * p, long *count)
{
    if (p->type != 't')
    {
        if (p->type == 'm')
        {
            *count += 1;
            count_migrations (p->next->back, count);
        }
        else
        {
            count_migrations (p->next->back, count);
            count_migrations (p->next->next->back, count);
        }
    }
}

MYREAL
prob_tree (world_fmt * world, timelist_fmt * tyme)
{
    long j, pop;
    MYREAL mm, cc, ss = 0;
    
    tyme->tl[0].interval = tyme->tl[0].age;
    for (j = 1; j < tyme->T; j++)
    {
        tyme->tl[j].interval = tyme->tl[j].age - tyme->tl[j - 1].age;
    }
    for (j = 0; j < tyme->T - 1; j++)
    {
        mm = cc = 0.0;
        for (pop = 0; pop < world->numpop; pop++)
        {
            mm += tyme->tl[j].lineages[pop] * world->mig0list[pop];
            cc +=
            tyme->tl[j].lineages[pop] * (tyme->tl[j].lineages[pop] -
                                         1) / world->param0[pop];
        }
        ss += -(tyme->tl[j].interval) * (mm + cc);
        if (tyme->tl[j].from == tyme->tl[j].to)
            ss += LOG2 - LOG (world->param0[tyme->tl[j].from]);
        else
            ss += LOG (world->param0[tyme->tl[j].from]);
    }
    return ss;
}
