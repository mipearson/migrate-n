#ifndef __MUTATIONMODEL
#define __MUTATIONMODEL
#include "migration.h"
#include "sighandler.h"

///
/// initialize the mutation model structure
void init_mutationmodel(world_fmt *world, data_fmt *data, option_fmt *options);
void finish_mutationmodel(world_fmt *world, data_fmt *data, option_fmt *options);
#endif
