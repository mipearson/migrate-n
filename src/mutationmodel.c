//
// mutation models
//
//
//
#include "migration.h"
#include "sighandler.h"

///
/// returns the number of states for a specific mutation model
int get_states(char datatype, world_fmt *world)
{
      return world->data->maxalleles[world->locus];
}
///
/// calculates (if necessary) and sets the base frequencies
void set_subloci_basefrequencies(mutationmodel_fmt *mumod, world_fmt *world, data_fmt *data, long sublocus)
{
  printf("sequence hack\n");
  mumod->basefreqs[0] = 0.25;
  mumod->basefreqs[1] = 0.25;
  mumod->basefreqs[2] = 0.25;
  mumod->basefreqs[3] = 0.25;
}

#define JC69 0
#define K2P  1

void set_siterates(long z, world_fmt *world, data_fmt *data, option_fmt *options)
{
  long i;
  mutationmodel_fmt *mumod = &world->mutationmodels[z];
  if(mumod->numsiterates==0)
    {
      mumod->numsiterates = options->rcategs;
      mumod->siterates = (double *) mycalloc(mumod->numsiterates,sizeof(double));
      mumod->siteprobs = (double *) mycalloc(mumod->numsiterates,sizeof(double));
    }
  for(i=0; i<mumod->numsiterates;i++)
    {
      mumod->siterates[i] = options->rrate[i];
      mumod->siteprobs[i] = options->probcat[i];
    }
}

void set_mutationmodel_eigenmaterial(long z, world_fmt *world, data_fmt *data, option_fmt *options) //eigenvectormatrix, inverse eigenvector matrix, eigenvalues
{
  mutationmodel_fmt *mumod = &world->mutationmodels[z];

  //beagle->eigenvectormatrix = world->mutationmodel->eigenvectormatrix;
  // JC69 model eigenvector matrix
  double jc69evec[4 * 4] = {
    1.0,  2.0,  0.0,  0.5,
    1.0,  -2.0,  0.5,  0.0,
    1.0,  2.0, 0.0,  -0.5,
    1.0,  -2.0,  -0.5,  0.0
  };
  double k2pevec[4 * 4] = {-1., 1., 0., -1., 1., 1., -1., 0., -1., 1., 0., 1., 1., 1., 1., 0};

  //beagle->inverseeigenvectormatrix = world->mutationmodel->inverseeigenvectormatrix;
  // JC69 model inverse eigenvector matrix
    double jc69ivec[4 * 4] = {
      0.25,  0.25,  0.25,  0.25,
      0.125,  -0.125,  0.125,  -0.125,
      0.0,  1.0,  0.0,  -1.0,
      1.0,  0.0,  -1.0,  0.0
    };
    double k2pivec[4 * 4] = {-0.25, 0.25, 0., -0.5, 0.25, 0.25, -0.25, 0., -0.25, 0.25, 0., 0.5, 0.25, 0.25, 0.5, 0. };
    
    // JC69 model eigenvalues
    double jc69eval[4] = { 0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333 };
    //K2P needs to be calculated because of kappa
    double kappa = 2.0;

    int named_model = JC69;
    const int numstates = mumod->numstates;
    const int numstates2 = numstates * numstates;
    mumod->eigenvalues = (double*) mycalloc(numstates,sizeof(double));
    mumod->eigenvectormatrix = (double*) mycalloc(numstates2,sizeof(double));
    mumod->inverseeigenvectormatrix = (double*) mycalloc(numstates2,sizeof(double));

    switch ((char) mumod->datatype)
      {
	// Models per datatype
	// Alleles
      case 'A':
      case 'a':
	break;
	// Msats
	// Msats Brownian
	// Sequence
      case 'S':
      case 's':      
	switch(named_model)//this should be through option!
	  {
	  case JC69:
	    memcpy(mumod->eigenvalues,jc69eval,sizeof(double) * 4);
	    memcpy(mumod->eigenvectormatrix,jc69evec,sizeof(double) * 16);
	    memcpy(mumod->inverseeigenvectormatrix,jc69ivec,sizeof(double) * 16);
	    break;
	  case K2P:
	    kappa = 2.0; //mumod->parameters[0];
	    mumod->eigenvalues[0] = -1.0;
	    mumod->eigenvalues[1] = 0.0;
	    mumod->eigenvalues[2] = mumod->eigenvalues[3] = -0.5 * (kappa + 1.0);
	    memcpy(mumod->eigenvectormatrix,k2pevec,sizeof(double) * 16);
	    memcpy(mumod->inverseeigenvectormatrix,k2pivec,sizeof(double) * 16);
	    break;
	    //case F84:
	  default:
	    break;
	}
      break;
      }
  // Sinlge nucleotide is the same as sequence but needs to calculate values for constant site patterns

}


///
/// initialize the mutation model structure
void init_mutationmodel(world_fmt *world, data_fmt *data, option_fmt *options)
{
  long totalmodels = 0L;
  long nummodels;
  const long numloci   = data->loci;
  unsigned long i;
  unsigned long j;
  unsigned long z;
  world->sublocistarts = (long *) mycalloc(data->loci, sizeof(long));
  world->nummutationmodels = (long *) mycalloc(data->loci, sizeof(long));
  world->maxnumpattern = (long *) mycalloc(data->loci, sizeof(long));
  for(i=0;i<numloci;i++)
    {
      nummodels = data->subloci[i];
      totalmodels += nummodels;
      if(i<numloci-1)
	world->sublocistarts[i+1] = totalmodels;
      world->nummutationmodels[i] = nummodels;
    } 
  world->mutationmodels = (mutationmodel_fmt *) mycalloc(totalmodels,sizeof(mutationmodel_fmt));

  // for each locus:
  z=0;
  for(i=0;i<numloci; i++)
    {
      // for each sublocus
      for(j=0; j< data->subloci[i]; j++)
	{
	  mutationmodel_fmt *mumod = &world->mutationmodels[z];
	  mumod->numpatterns = 0L; // number unique patterns
	  mumod->numsites    = 0L; // number of sites
	  mumod->numstates   = 0L; // number of states in model: DNA=4, DNA+gap=5, msat>2
	  mumod->datatype        = data->datatype[z]; //specifices the model
	  z++;
	}
    }
}

///
/// finish the mutation model structure
void finish_mutationmodel(world_fmt *world, data_fmt *data, option_fmt *options)
{
  //unsigned long totalmodels = 0L;
  //unsigned long nummodels;
  // const unsigned long numloci   = world->loci;
  const unsigned long locus = world->locus;
  //unsigned long i;
  unsigned long j;
  unsigned long z = 0L;

  z = world->sublocistarts[locus];
  // for each sublocus
  for(j=0; j< data->subloci[locus]; j++)
    {
	  mutationmodel_fmt *mumod = &world->mutationmodels[z];
	  mumod->numpatterns = world->data->seq[0]->endsite; // number unique patterns
	  mumod->numsites    = world->data->seq[0]->endsite; // number of sites
	  if(world->maxnumpattern[locus] < mumod->numpatterns)
	    {
	      world->maxnumpattern[locus] = mumod->numpatterns;
	    }
	  mumod->datatype        = data->datatype[z]; //specifices the model
	  //	  mumod->numstates   = get_states(mumod->datatype, world); // number of states in model: DNA=4, DNA+gap=5, msat>2
	  //set_subloci_basefrequencies(mumod,world, data,j);  
	  set_mutationmodel_eigenmaterial(z,world, data, options); //eigenvectormatrix, inverse eigenvector matrix, eigenvalues
	  set_siterates(z,world,data,options);
	  z++;
    }
}
