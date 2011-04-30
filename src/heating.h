#ifndef __HEATING_H__
#define __HEATING_H__
/* heating scheme using a thread pool if available
   
Peter Beerli 2000
 
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
$Id: heating.h 1476 2008-12-30 04:38:59Z beerli $
 
*/

#ifdef PTHREADS
#include <stdio.h>
#include <pthread.h>

typedef struct tpool_work
{
    void (*routine) ();
    void *arg;
    struct tpool_work *next;
}
tpool_work_t;

typedef struct _tpool_t
{
    int num_threads;
    int max_queue_size;
    int do_not_block_when_full;
    pthread_t *threads;
    int cur_queue_size;
    tpool_work_t *queue_head;
    tpool_work_t *queue_tail;
    pthread_mutex_t random_lock;
    pthread_mutex_t queue_lock;
    pthread_cond_t queue_done;
    pthread_cond_t queue_not_empty;
    pthread_cond_t queue_not_full;
    pthread_cond_t queue_empty;
    int queue_closed;
    int shutdown;
    int done;
}
_tpool_t;

typedef _tpool_t *tpool_t;

extern void tpool_init (tpool_t * tpoolp, int num_worker_threads,
                            int max_queue_size, int do_not_block_when_full);
extern void fill_tpool (tpool_t tpool, world_fmt ** universe,
                            int universe_size);
extern void wait_tpool (tpool_t tpoolp, int usize);
extern int tpool_destroy (tpool_t tpool, int finish);
extern int tpool_synchronize (tpool_t tpool, int finish);
extern int tpool_add_work (tpool_t tpool, void *routine, void *arg, long z);
#else



#endif /*PTHREADS*/


extern void adjust_temperatures(world_fmt ** universe, long hchains, long step, long steps);
extern void adjust_temperatures_bounded(world_fmt ** universe, long hchains, long step, long steps);
#endif
