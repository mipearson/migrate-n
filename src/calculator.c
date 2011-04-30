#ifdef BEAGLE
#include "beagle.h"
#include "migration.h"
#include "tree.h"
#include "tools.h"
#include "sighandler.h"

#define BEAGLE_PARTIALS 7

extern long myID;

void debug_beagle(beagle_fmt *beagle);
long new_id(long id, long sumtips);
/*-----------------------------------------------------------------------------
|	This function sets up the beagle library and initializes all data members.
*/
void init_beagle(world_fmt *world, long locus)
{

  beagle_fmt *beagle = world->beagle;
  beagle->instance_handle    = (int *) mycalloc(world->nummutationmodels[locus],sizeof(int));
  beagle->numallocpartials    = 4 * world->sumtips*world->mutationmodels[0].numstates*world->mutationmodels[0].numsites;
  beagle->partials            = (double *) mycalloc(beagle->numallocpartials, sizeof(double));

  beagle->numallocbranches   = 200;
  beagle->branch_lengths     = (double *) mycalloc(beagle->numallocbranches, sizeof(double));
  beagle->numallocoperations = 200;
  beagle->operations         = (int *) mycalloc(beagle->numallocoperations * BEAGLE_PARTIALS,sizeof(int));//parent,parentIDscalingfactor, lchild,lchildtrans,rchild,rtrans
  beagle->branch_indices     = (int *) mycalloc(beagle->numallocoperations * 2, sizeof(int));

  beagle->scalingfactorsindices= (int *) mycalloc(2*(2 * world->sumtips-1), sizeof(int));
  beagle->scalingfactorscount  = 0; //actual count of which scale factors are updated? this was world->sumtips-1;

  beagle->numallocallyweights = world->mutationmodels[0].numsites;
  beagle->allyweights = (int *) mycalloc(beagle->numallocallyweights,sizeof(int));
  // this is already set and definitely should be set here, left as reminder
  // that the current setting in makevalues() is a hack
  //  world->mutationmodels[world->sublocistarts[locus]].basefreqs = (double *) calloc(world->mutationmodels[world->sublocistarts[locus]].numstates,sizeof(double));
} 

void  set_beagle_instances(world_fmt *world, boolean instance)
{
  int resource=1;//cpu=0, gpu=1
  long i;
  long ii;
  long locus = world->locus;
  beagle_fmt *beagle = world->beagle;
  for(i=0;i<world->nummutationmodels[locus];i++)
    {
      ii = world->sublocistarts[locus] + i;
      beagle->instance_handle[i] = beagleCreateInstance(world->sumtips,       // number of tips
							2 * (2 * world->sumtips-1), // total buffers (= total nodes in tree)
							// twice to accomodate rejections.
						  0,			// number of compact state representation [funny tips]
						  world->mutationmodels[ii].numstates, // number of states (nucleotides etc)
						  world->mutationmodels[ii].numsites,  // number of site patterns
						  1,			// number of rate matrices eigen-decomp buffers
							2*(2 * world->sumtips - 2),// number of rate matrix buffers (=number of branches)
						  world->mutationmodels[ii].numsiterates,// categoryCount
							0, //2*(2 * world->sumtips - 1),   number of scaling buffers
                                                  &resource,		        // List of resources (?)
						  1,			// Number of resources
						  0L,			// preferenceFlags (see BeagleFlags)
						  0L			// requirementFlags (see BeagleFlags) 
						  );

      if (beagle->instance_handle[i] < 0)
	usererror("createInstance returned a negative instance handle (and that's not good)");
      // initialize the instance
      int code = beagleInitializeInstance(beagle->instance_handle[i], NULL);
      
      if (code < 0) 
	{
	  usererror("initializeInstance returned a negative error code (and that is bad)\n\n");
	}
    }
}


long
set_branch_index (node * p,  long *bid)
{
  long bbid = *bid;
    if (p->type != 't')
    {
      bbid = set_branch_index (crawlback (p->next), &bbid);
      bbid = set_branch_index (crawlback (p->next->next), &bbid);
    }
    p->bid = bbid;
    //printf("%li -- %li\n",p->id, bbid);
    return bbid+1;
}    /* set_branch_index */


void adjust_beagle(beagle_fmt *beagle)
{
  beagle->numallocoperations = beagle->numoperations + 10;
  beagle->operations = (int *) myrealloc(beagle->operations, sizeof(int) * BEAGLE_PARTIALS * beagle->numallocoperations);
  beagle->branch_lengths   = (double *) myrealloc(beagle->branch_lengths, sizeof(double) * 2 * beagle->numallocoperations);
  beagle->branch_indices = (int *) myrealloc(beagle->branch_indices, beagle->numallocoperations * 2 * sizeof(int));
  beagle->numallocbranches = beagle->numallocoperations * 2;
}


void prepare_beagle_instances(node *theNode, node * left, node *right, beagle_fmt *beagle)
{
  long parentid = theNode->id;
  long rightid  = right->id;
  long leftid   = left->id;
 
  long rightbid  = right->bid;
  long leftbid   = left->bid;

  double leftbranch = left->v;
  double rightbranch = right->v;
  long ii = beagle->numoperations;
  long i = ii * BEAGLE_PARTIALS ;
  long j = ii * 2;
  if(beagle->numoperations  >=  beagle->numallocoperations)
    {
      adjust_beagle(beagle);
    }
  beagle->branch_indices[j] = leftbid;
  beagle->branch_indices[j+1] = rightbid;

  beagle->operations[i]   = parentid;
  beagle->operations[i+1]   = parentid; //BEAGLE_OP_NONE; //parentid;
  beagle->operations[i+2]   = parentid; //BEAGLE_OP_NONE; //parentid;
  beagle->operations[i+3] = leftid;
  beagle->operations[i+4] = leftbid;
  beagle->operations[i+5] = rightid;
  beagle->operations[i+6] = rightbid;

  beagle->branch_lengths[j]     = leftbranch; 
  beagle->branch_lengths[j+1]   = rightbranch; 

  beagle->scalingfactorsindices[ii] = parentid;
  beagle->scalingfactorscount += 1;
  beagle->numbranches += 2;
  beagle->numoperations++;
  //#ifdef BEAGLEDEBUG
    printf("%li> TREE: %i {%i %i %i %i %i %i %i} {%i %i} {%f %f}\n",myID, beagle->numoperations,
	 beagle->operations[i],   
	 beagle->operations[i+1] ,
	 beagle->operations[i+2] ,
	 beagle->operations[i+3] ,
	 beagle->operations[i+4] ,
	 beagle->operations[i+5] ,
	 beagle->operations[i+6] ,
	 beagle->branch_indices[j],
	 beagle->branch_indices[j+1],                 
	 beagle->branch_lengths[j],
	 beagle->branch_lengths[j+1]);
  //#endif
  //  printf("c");
}

void prepare_beagle_instances_proposal(proposal_fmt *proposal, long trueparentid, long leftid, long leftbid, double leftbranch, 
				       long rightid, long rightbid, double rightbranch, beagle_fmt *beagle)
{
  long i = beagle->numoperations * BEAGLE_PARTIALS;
  long j = beagle->numoperations * 2;
  long nodep_boundary = proposal->world->sumtips * 2 - 1;
  long parentid = trueparentid > nodep_boundary ? trueparentid - nodep_boundary : trueparentid + nodep_boundary; 
  if(beagle->numoperations  >=  beagle->numallocoperations)
    {
      adjust_beagle(beagle);
    }
  printf("PINT: (%li, %li) l:%li,%f, r:%li,%f\n",trueparentid,parentid,leftid,leftbranch,rightid,rightbranch);
  beagle->branch_indices[j] = leftbid;
  beagle->branch_indices[j+1] = rightbid;
  beagle->operations[i]   = parentid;
  beagle->operations[i+1]   = parentid; //BEAGLE_OP_NONE;//parentid;
  beagle->operations[i+2]   = parentid; //BEAGLE_OP_NONE; //parentid;
  beagle->operations[i+3] = leftid;
  beagle->operations[i+4] = leftbid;
  beagle->operations[i+5] = rightid;
  beagle->operations[i+6] = rightbid;

  beagle->branch_lengths[j]     = leftbranch; 
  beagle->branch_lengths[j+1]   = rightbranch; 
  beagle->numbranches += 2;
  beagle->numoperations++;
  //#ifdef BEAGLEDEBUG
  printf("%li> PROP: %i {%i %i %i %i %i %i %i} {%i %i} {%f %f}\n",myID,beagle->numoperations,
	 beagle->operations[i],   
	 beagle->operations[i+1] ,
	 beagle->operations[i+2] ,
	 beagle->operations[i+3] ,
	 beagle->operations[i+4] ,
	 beagle->operations[i+5] ,
	 beagle->operations[i+6] ,
	 beagle->branch_indices[j],
	 beagle->branch_indices[j+1],                 
	 beagle->branch_lengths[j],
	 beagle->branch_lengths[j+1]);
  //#endif
  //printf("p");
}

void
evaluate_beagle_instances_proposal (proposal_fmt * proposal,
				    node * mother,  
				    node * newdaughter, long newdaughter_id, long newdaughter_bid, 
				    MYREAL v)
{
    node *nn = NULL, *d1 = NULL, *d2 = NULL, *oldnn = NULL;

    if (mother->type != 'r')
    {
        children (mother, &d1, &d2);
        if (d1 == newdaughter)
	  {
	    d1 = d2;
	  }
	long id0, id1, id2;
	long bid1, bid2;
	double v1, v2;
	id0  = mother->id; 
	id1  = new_id(newdaughter_id, proposal->world->sumtips);
	bid1 = newdaughter_bid;
	v1 = v;
	id2  = d1->id;
	bid2 = d1->bid;
	v2 = d1->v;
	prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);

        oldnn = mother;
        nn = showtop (crawlback (mother));
        while (nn->type != 'r')
        {
            children (nn, &d1, &d2);
            if (d1 == oldnn)
	      {
		d1 = d2;
		d2 = oldnn;
	      }
	    if(d2->id != oldnn->id)
	      {
		warning("One of the children is not what it should be\n");
	      }
	    id0  = nn->id; 
	    id1  = new_id(d2->id,proposal->world->sumtips);;
	    bid1 = d2->bid;
	    v1 = d2->v;
	    id2  = d1->id;
	    bid2 = d1->bid;
	    v2 = d1->v;
	    prepare_beagle_instances_proposal(proposal, id0, id1, bid1, v1, id2, bid2, v2, proposal->world->beagle);

            oldnn = nn;
            nn = showtop (crawlback (nn));
        }
    }
}




void fill_beagle_instances(world_fmt *world)
{
  beagle_fmt *beagle = world->beagle;
  unsigned long i;
  unsigned long ii;
  unsigned long j;
  //unsigned long zz;
  unsigned long site;
  unsigned int numsites;
  unsigned int numstates;
  unsigned int numweights=0;
  int code;
  long locus = world->locus;
  // set partials for all tipnodes
  // use z to advance through the partial array
  unsigned long z = 0L;

  for(i=0; i<world->nummutationmodels[locus]; i++)
    {
      ii = world->sublocistarts[locus] + i;
      numstates = world->mutationmodels[ii].numstates;
      numsites = world->mutationmodels[ii].numsites;
      // will currently fail with more complex data
      if(numsites + numweights < beagle->numallocallyweights)
	{
	  beagle->numallocallyweights += numsites * 2;
	  beagle->allyweights = (int *) myrealloc(beagle->allyweights,sizeof(int) * beagle->numallocallyweights);
	}
      memcpy(beagle->allyweights+numweights, world->data->seq[0]->aliasweight, sizeof(int) * numsites);
      numweights += numsites;
      for (j = 0; j < world->sumtips; ++j)
	{
	  for(site=0;site<numsites;site++)
	    {
	      if(z > beagle->numallocpartials)
		{
		  beagle->numallocpartials += numsites * numstates;
		  beagle->partials = (double *) myrealloc(beagle->partials,sizeof(double) * beagle->numallocpartials);
		}
	      memcpy(&beagle->partials[z], 
		     &(world->nodep[j]->x.s[site][0][0]),
		     sizeof(double) * numstates);
	      beagle->numpartials++;
#if 0
	      fprintf(stdout,"Site %li Partial: %li (%f ", site, j, 
		      beagle->partials[z]);
	      for(zz=1;zz<world->mutationmodels[ii].numstates;zz++)
		{
		  fprintf(stdout,"%f ",beagle->partials[z+zz]);
		}
	      fprintf(stdout,")\n")
#endif;
	      // advance z
	      z += numstates;
	    }
	  code = beagleSetTipPartials(
			     beagle->instance_handle[i],			// instance
			     j,							// indicator for tips
			     &beagle->partials[z - numsites*numstates]);// inPartials
	  if (code != 0)
	    usererror("setTipPartials encountered a problem");
	}
    }

  
  for(i=0; i<world->nummutationmodels[locus]; i++)
    {
      ii = world->sublocistarts[locus] + i;

      code = beagleSetCategoryRates(beagle->instance_handle[i], world->mutationmodels[ii].siterates);
      if (code != 0)
	usererror("setCategoryRates encountered a problem");

      code = beagleSetEigenDecomposition(beagle->instance_handle[i],		  // instance
				   0,					  // eigenIndex,
				   (const double *)world->mutationmodels[ii].eigenvectormatrix,	// inEigenVectors,
				   (const double *)world->mutationmodels[ii].inverseeigenvectormatrix,// inInverseEigenVectors,
				   world->mutationmodels[i].eigenvalues); // inEigenValues
      
      if (code != 0)
	usererror("setEigenDecomposition encountered a problem");
    }
}
/*-----------------------------------------------------------------------------
|	Calculates the log likelihood by calling the beagle functions
|	updateTransitionMatrices, updatePartials and calculateEdgeLogLikelihoods.
*/
double calcLnL(world_fmt *world, boolean instance)
{
  beagle_fmt *beagle = world->beagle;
  double logL = 0.0;
  unsigned long i;
  unsigned long j;
  //unsigned long z;
  long locus = world->locus;
  long ii;
  double *patternloglike = (double *) calloc(world->maxnumpattern[locus],sizeof(double));
   double *outlike = (double *) calloc(10*world->maxnumpattern[locus],sizeof(double));
  int rootIndex = beagle->operations[BEAGLE_PARTIALS * (beagle->numoperations-1)];//world->root->next->back->id;
  for(i=0; i<world->nummutationmodels[locus]; i++)
    {
      ii = world->sublocistarts[locus] + i;
      int code = beagleUpdateTransitionMatrices(beagle->instance_handle[i],		// instance,
					  0,					// eigenIndex,
					  (const int *) beagle->branch_indices,	// indicators transitionrates for each branch,
					  NULL, 			        // firstDerivativeIndices,
					  NULL,					// secondDervativeIndices,
					  beagle->branch_lengths,		// edgeLengths,
					  beagle->numbranches);			// number branches to update, count

	if (code != 0)
		usererror("updateTransitionMatrices encountered a problem");

	int cumulativeScalingFactorIndex = 0; //BEAGLE_OP_NONE; //this would be the index of the root scaling location 
	
	beagleResetScaleFactors(beagle->instance_handle[i],
			  cumulativeScalingFactorIndex);

	beagleAccumulateScaleFactors(beagle->instance_handle[i],
			       beagle->scalingfactorsindices,
			       beagle->scalingfactorscount,
			       cumulativeScalingFactorIndex);
	

	code = beagleUpdatePartials((const int *) &beagle->instance_handle[i],	// instance
			      1,					        // instanceCount
			      beagle->operations,		                // operations
			      beagle->numoperations,				// operationCount
				    cumulativeScalingFactorIndex);//BEAGLE_OP_NONE);					        // connected to accumulate....
#ifdef BEAGLEDEBUG
	for(j=0;j<2*(world->sumtips * 2 - 1); j++)
	  {
	    beagleGetPartials(beagle->instance_handle[i],j,BEAGLE_OP_NONE, outlike);
	    if(j==world->sumtips * 2 - 1)
	      printf("-----------------------\n");
	    printf("%li: {%f, %f, %f, %f}\n",j, outlike[0],outlike[1],outlike[2],outlike[3]);
	  }
#endif
	if (code != 0)
	  usererror("updatePartials encountered a problem");



	if(beagle->weights==NULL)
	  beagle->weights = (double *) mycalloc(1,sizeof(double));
	//	else
	//  beagle->weights = (double *) myrealloc(beagle->weights,beagle->numoperations * sizeof(double));

	beagle->weights[0]= 1.0;

	// calculate the site likelihoods at the root node
	code = beagleCalculateRootLogLikelihoods(beagle->instance_handle[i],         // instance
					   (const int *) &rootIndex, // bufferIndices
					   (const double *) world->mutationmodels[ii].siteprobs,   // weights
					   (const double *) world->mutationmodels[ii].basefreqs,// stateFrequencies
					   &cumulativeScalingFactorIndex, //scalingfactors index,
					   1,              // count is this correct
	                            patternloglike);         // outLogLikelihoods
	//trash from function above					   &beagle->scalingfactorscount,//size of the scaling factor index
	if (code != 0)
		usererror("calculateRootLogLikelihoods encountered a problem");

	for (j = 0; j < world->mutationmodels[ii].numsites; j++) 
	  {
		logL += beagle->allyweights[j] * patternloglike[j];
		//		printf("%.1f ",patternloglike[j]);
	}
    }
  printf("Log LnL=%f (instance=%li)\n",logL,(long) instance);
#ifdef BEAGLEDEBUG
  debug_beagle(beagle);
#endif
  myfree(patternloglike);
  myfree(outlike);
  return logL; 
}

double force_beagle_recalculate(world_fmt *world, long locus)
{
  reset_beagle(world->beagle);
  set_all_dirty(world->root->next, crawlback (world->root->next), world, locus);
  smooth (world->root->next, crawlback (world->root->next), world, locus);
  return treelikelihood(world);
}

///
/// start out at the root node and traverse up and match the ids with the new ids
void change_beagle(node *theNode, beagle_fmt *beagle, long sumtips)
{
  // move node ids;
  // this works only just after the pseudolikelihood step + the acceptlike with TRUE outcome
  long i;
  long parentid;
  long oldparentid;
  long nodep_boundary = sumtips * 2 - 1;
  if (theNode->type != 't')
    {
      change_beagle(crawlback(theNode->next), beagle,sumtips);
      change_beagle(crawlback(theNode->next->next),beagle, sumtips);
    }
  for(i=0;i<beagle->numoperations;i++)
    {
      parentid = beagle->operations[i*BEAGLE_PARTIALS];
      oldparentid = parentid < nodep_boundary ? parentid + nodep_boundary : parentid - nodep_boundary; 
      if(theNode->id == oldparentid)
	{
	  theNode->id = parentid;
	}
    }
}

long new_id(long id, long sumtips)
{
  const long nodep_boundary = 2 * sumtips - 1;
  return (id < nodep_boundary ? id + nodep_boundary : id - nodep_boundary); 
}

void reset_beagle(beagle_fmt *beagle)
{
  beagle->numoperations = 0;
  beagle->numbranches   = 0;
  beagle->scalingfactorscount = 0;
}

void beagle_stop(world_fmt **universe, long usize)
{
  long u;
  for(u=0;u< usize; u++)
    {
      world_fmt *world = universe[u];
      beagle_fmt *beagle = world->beagle;
      unsigned long i;
      long locus = world->locus;
      for(i=0; i<world->nummutationmodels[locus]; i++)
	{
	  beagleFinalizeInstance(beagle->instance_handle[i]);
	}
    }
}


void set_beagle_dirty(node *origin, node *target, node *mrca)
{
  node *nn = origin;
  while((nn=showtop(crawlback(nn)))!=mrca)
    {
      set_dirty(nn);
    }
  nn = target;
  while((nn=showtop(crawlback(nn)))!=mrca)
    {
      set_dirty(nn);
    }
}

void debug_beagle(beagle_fmt *beagle)
{
  printf("----beagle content---------------------\n");
  printf("operations    : %p\n",beagle->operations);
  printf("     alloc    : %i\n",beagle->numallocoperations);
  printf("       num    : %i\n",beagle->numoperations);
  printf("branch indices: %p\n",beagle->branch_indices);
  printf("     alloc    : %i\n",beagle->numallocbranches);
  printf("       num    : %i\n",beagle->numbranches);
  printf("----beagle content end-----------------\n\n");
}



#endif /*BEAGLE*/
