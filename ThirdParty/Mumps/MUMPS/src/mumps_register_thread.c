/*
 *
 *  This file is part of MUMPS 5.9.1, released
 *  on Mon Jul 20 09:00:43 UTC 2026
 *
 *
 *  Copyright 1991-2026 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
#if defined(USE_GPU)
#include "mumps_compat.h"
#if defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) || defined(USE_XKBLAS)
/* Async. pinning not available for Windows */
void mumps_register_thread_return()
{
}
#else
#include <pthread.h>
#include <errno.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include "mumps_register_thread.h"
/*
*/
typedef enum mumps_register_request_op {
  MUMPS_REGISTER_REQUEST,
  MUMPS_UNREGISTER_REQUEST
} mumps_register_request_op_t;
/*
*/
typedef struct mumps_request {
  mumps_register_request_op_t op;
  void* ptr;
  size_t size;
} mumps_request_t;
typedef struct mumps_request_queue {
#define MUMPS_SIZE_INQUEUE 32
  pthread_t         thread;
  volatile uint64_t pos_write;
  volatile uint64_t pos_read;
  volatile uint64_t pos_signal;
  pthread_mutex_t   lock;
  pthread_cond_t    cond;
  pthread_cond_t    cond_wait;
  volatile int      is_term;
  mumps_request_t req[MUMPS_SIZE_INQUEUE];
} mumps_request_queue_t;
/* array of lists of requests */
static int all_queue_size = 0;
static mumps_request_queue_t* all_queue = 0;
/*
*/
static void mumps_request_queue_init(mumps_request_queue_t* queue )
{
#if defined(DEBUG_GPU)
  printf("=> mumps_request_queue_init()");
#endif
  queue->pos_write = 0;
  queue->pos_read = 0;
  queue->pos_signal  = 0;
  queue->is_term = 0;
  pthread_mutex_init(&queue->lock, 0);
  pthread_cond_init(&queue->cond, 0);
  pthread_cond_init(&queue->cond_wait, 0);
/*
  assert(0 == pthread_mutex_init(&queue->lock, 0));
  assert(0 == pthread_cond_init(&queue->cond, 0));
  assert(0 == pthread_cond_init(&queue->cond_wait, 0));
*/
  memset(queue->req, 0, sizeof(queue->req));
#if defined(DEBUG_GPU)
  printf("<= mumps_request_queue_init()");
#endif
}
/*
*/
static void mumps_request_queue_destroy( mumps_request_queue_t* queue )
{
#if defined(DEBUG_GPU)
  printf("=> mumps_request_queue_destroy()\n");
#endif
  pthread_mutex_lock(&queue->lock);
  queue->is_term = 1;
  pthread_cond_signal(&queue->cond);
  pthread_mutex_unlock(&queue->lock);
  void* dummy;
  /* error code sur pthread_join ??? */
  pthread_join( queue->thread, &dummy );
  pthread_mutex_destroy(&queue->lock);
  pthread_cond_destroy(&queue->cond);
  pthread_cond_destroy(&queue->cond_wait);
/*
  assert(0 == pthread_mutex_destroy(&queue->lock));
  assert(0 == pthread_cond_destroy(&queue->cond));
  assert(0 == pthread_cond_destroy(&queue->cond_wait));
*/
#if defined(DEBUG_GPU)
  printf("<= mumps_request_queue_destroy()\n");
#endif
}
/* main thread entrypoint */
void* mumps_register_thread( void* dummy )
{
#if defined(DEBUG_GPU)
  printf("=> mumps_register_thread \n");
#endif
  mumps_request_queue_t* queue = (mumps_request_queue_t*)dummy;
  pthread_mutex_lock(&queue->lock);
  while (1)
  {
    while (queue->pos_write == queue->pos_read && !queue->is_term)
      pthread_cond_wait(&queue->cond, &queue->lock);
    assert(queue->is_term || (queue->pos_write > queue->pos_read) );
    if (queue->pos_write > queue->pos_read)
    {
      uint64_t index = queue->pos_read % MUMPS_SIZE_INQUEUE;
      mumps_register_request_op_t op;
      void* ptr;
      size_t size;
      ptr  = queue->req[index].ptr;
      assert( ptr != 0 );
      op   = queue->req[index].op;
      size = queue->req[index].size;
      ++queue->pos_read;
      /* request copied back, unlock queue */
      pthread_mutex_unlock(&queue->lock);
      cudaError_t err;
      if (op == MUMPS_REGISTER_REQUEST)
      {
        err = cudaHostRegister(ptr, size, cudaHostRegisterPortable);
        if (!( (cudaSuccess == err) || (cudaErrorHostMemoryAlreadyRegistered == err)))
        {
          printf("***[%s]: cudaHostRegister ptr=%p, size=%li, error: %i\n", __func__, ptr, size, err);
        }
      }
      else if (op == MUMPS_UNREGISTER_REQUEST)
      {
        err = cudaHostUnregister(ptr);
        if (!( (cudaSuccess == err) || (cudaErrorHostMemoryNotRegistered == err)))
        {
          printf("***[%s]: cudaHostUnregister ptr=%p, size=%li, error: %i\n", __func__, ptr, size, err);
        }
      }
      /* signal ok: */
      pthread_mutex_lock(&queue->lock);
      ++queue->pos_signal;
      /* broadcast to all waiting threads */
      pthread_cond_broadcast(&queue->cond_wait);
    }
    else if (queue->is_term) break;
  }
  pthread_mutex_unlock(&queue->lock);
  return 0;
}
/*
*/
void mumps_init_request_pool(int size_pool)
{
#if defined(DEBUG_GPU)
  printf("=> mumps_init_request_pool()\n");
#endif
  int i, err;
  all_queue_size = size_pool;
  all_queue = (mumps_request_queue_t*)malloc(all_queue_size* sizeof(mumps_request_queue_t));
  for (i=0; i<all_queue_size; ++i)
  {
    mumps_request_queue_init(&all_queue[i]);
    err = pthread_create(&all_queue[i].thread, 0, mumps_register_thread, &all_queue[i] );
    assert(err ==0);
  }
#if defined(DEBUG_GPU)
  printf("<= mumps_init_request_pool()\n");
#endif
}
/*
*/
void mumps_finalize_request_pool()
{
#if defined(DEBUG_GPU)
  printf("=> mumps_finalize_request_pool()\n");
#endif
  int i;
  void* tmp;
  for (i=0; i<all_queue_size; ++i)
  {
    pthread_mutex_lock(&all_queue[i].lock);
    all_queue[i].is_term = 1;
    pthread_cond_signal( &all_queue[i].cond );
    pthread_mutex_unlock(&all_queue[i].lock);
    assert(0==pthread_join( all_queue[i].thread, &tmp ));
  }
  free( all_queue );
  all_queue_size = 0;
  all_queue = 0;
#if defined(DEBUG_GPU)
  printf("<= mumps_finalize_request_pool()\n");
#endif
}
#define HASH_FUNC(pth,sz) (((uint64_t)(pth)) % (sz))
/*
*/
uint64_t mumps_post_request(
    mumps_register_request_op_t op,
    void* ptr, size_t size
)
{
#if defined(DEBUG_GPU)
  printf("=> mumps_post_request\n");
#endif
  /* ne pas oublier d'appeller xkblas_ si configuré avec */
  if (ptr ==0) return (uint64_t)-1;
  /* Hash function from the thread id to request register list */
  int hash = HASH_FUNC( pthread_self(), all_queue_size);
  mumps_request_queue_t* queue = &all_queue[hash];
  pthread_mutex_lock( &queue->lock );
  while (queue->pos_write - queue->pos_read >= MUMPS_SIZE_INQUEUE)
    pthread_cond_wait( &queue->cond_wait, &queue->lock );
  uint64_t index = queue->pos_write++;
  int idx = index % MUMPS_SIZE_INQUEUE;
  queue->req[idx].op   = op;
  queue->req[idx].size = size;
  queue->req[idx].ptr  = ptr;
  pthread_cond_signal( &queue->cond );
  pthread_mutex_unlock( &queue->lock );
#if defined(DEBUG_GPU)
  printf("<= mumps_post_request\n");
#endif
  return index;
}
/*
*/
uint64_t mumps_host_register_async(void *ptr, size_t size)
{
#if defined(DEBUG_GPU)
  printf("RT: mumps_host_register_async \n");
#endif
  return mumps_post_request( MUMPS_REGISTER_REQUEST, ptr, size );
}
/*
*/
uint64_t mumps_host_unregister_async(void *ptr, size_t size)
{
#if defined(DEBUG_GPU)
  printf("RT: mumps_host_unregister_async \n");
#endif
  return mumps_post_request( MUMPS_UNREGISTER_REQUEST, ptr, size );
}
/*
*/
int mumps_request_testwait(
    uint64_t index,
    int flag
)
{
  /* Hash function from the thread id to request register list */
  int hash = HASH_FUNC( pthread_self(), all_queue_size);
  mumps_request_queue_t* queue = &all_queue[hash];
  if ((flag != 2) && ((index == (uint64_t)-1) || (index >= queue->pos_write)))
    return EINVAL;
  switch (flag)
  {
    case 0: /* test: no lock needed, read is atomic, write on pos_signal after lock -> memory is ok */
      if (index <= queue->pos_signal) return 0;
      return EINPROGRESS;
    case 1: /* wait */
      pthread_mutex_lock(&queue->lock);
      while (index >= queue->pos_signal)
        pthread_cond_wait( &queue->cond_wait, &queue->lock );
      pthread_mutex_unlock(&queue->lock);
      return 0;
    case 2: /* wait all */
    {
      pthread_mutex_lock(&queue->lock);
      uint64_t capture_pos_write = queue->pos_write;
      while (capture_pos_write > queue->pos_signal)
        pthread_cond_wait( &queue->cond_wait, &queue->lock );
      pthread_mutex_unlock(&queue->lock);
      return 0;
    }
    default:
      return EINVAL;
  }
  return 0;
}
/* Return 0 if request is completed, else EINPROGRESS or EINVAL if invalid arg.
*/
int mumps_request_test(
    uint64_t index
)
{
  return mumps_request_testwait(
    index,
    0
  );
}
/*
*/
int mumps_request_wait(
    uint64_t index
)
{
  int err = mumps_request_testwait(
    index,
    1
  );
  assert( err != EINPROGRESS );
  return (err == 0 ? 0 : EINVAL);
}
/*
*/
int mumps_request_wait_all(
/*    uint64_t index */
)
{
  int err = mumps_request_testwait(
    0,
    2
  );
  assert( err != EINPROGRESS );
  return (err == 0 ? 0 : EINVAL);
}
#endif
#else
void mumps_register_thread_return()
{
  /*
   * Registering tools will be available in the future.
   */
}
#endif
