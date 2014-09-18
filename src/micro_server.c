/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "micro_server.h"
#include "microscale_information.h"

#include <stdlib.h>
#include <string.h>

struct work_queue{
  pgf_FE2_job *jobs;
  size_t current;
  size_t size;
  size_t max_size;
};
typedef struct work_queue work_queue;

/**
 * Initialize a work queue object.
 */
void work_queue_init(work_queue *queue)
{
  queue->jobs = NULL;
  queue->current = 0;
  queue->size = 0;
  queue->max_size = 0;
}

/**
 * Destroy a work queue object.
 */
void work_queue_destroy(work_queue *queue)
{
  free(queue->jobs);
}

/**
 * Add a job to the end of the work queue.
 */
void work_queue_push(work_queue *queue,
		     const pgf_FE2_job *job)
{
  static size_t increment = 4;
  if(queue->size >= queue->max_size){
    queue->max_size += increment;
    queue->jobs = realloc(queue->jobs,queue->max_size*sizeof(*(queue->jobs)));
  }
  memcpy((queue->jobs) + queue->size,job,sizeof(*job));
  queue->size++;
}

/**
 * Return true if queue is empty.
 */
int work_queue_is_empty(const work_queue *queue)
{
  return (queue->size <= 0);
}

/**
 * Get pointer to current job. Returns NULL if queue is empty.
 */
pgf_FE2_job* work_queue_get_current_job(work_queue *queue)
{
  if(!work_queue_is_empty(queue)){
    return (queue->jobs) + queue->current;
  } else {
    return NULL;
  }
}

/**
 * Make first job current and return pointer. Returns NULL if queue is
 * empty.
 */
pgf_FE2_job* work_queue_first(work_queue *queue)
{
  queue->current = 0;
  return work_queue_get_current_job(queue);
}

/**
 * Make last job current and return pointer. Returns NULL if queue is
 * empty.
 */
pgf_FE2_job* work_queue_last(work_queue *queue)
{
  queue->current = queue->size - 1;
  return work_queue_get_current_job(queue);
}

/**
 * Make next job in queue current. Returns true if queue is not empty.
 */
int work_queue_next(work_queue *queue)
{
  queue->current++;
  if(queue->current >= queue->size) queue->current = 0;
  return (!work_queue_is_empty(queue));
}

/**
 * Remove the current job from the queue.
 */
void work_queue_pop(work_queue *queue)
{
  /* move valid entries after current position to current position. */
  memmove((queue->jobs) + queue->current,
	  (queue->jobs) + queue->current + 1,
	  (queue->size - (queue->current + 1))*sizeof(*(queue->jobs)));

  /* decrement the size -- number of valid jobs in queue */
  queue->size--;
}

/**
 * Sort the work queue jobs according to the sort function.
 */
void work_queue_sort(work_queue *queue,
		     int (*sort)(const void *, const void *))
{
  if(!work_queue_is_empty(queue)){
    qsort(queue->jobs,queue->size,sizeof(*(queue->jobs)),sort);
  }
}

int pgf_new_micro_server_master(const PGFEM_mpi_comm *comm)
{
  int err = 0;

  /* look for partition/server info on pgf_mpi_comm->mm_inter */
  {
    int len = 0;
    char *buffer = NULL;
    {
      int flag = 0;
      MPI_Status rebal_stat;
      while(!flag){
	err += MPI_Iprobe(MPI_ANY_SOURCE,-1,comm->mm_inter,&flag,&rebal_stat);
      }

      /* post appropriate receive for partition information from macroscale */
      err += MPI_Get_count(&rebal_stat,MPI_CHAR,&len);
      buffer = calloc(len,sizeof(*buffer));
      MPI_Request rebal_req;
      err += MPI_Irecv(buffer,len,MPI_CHAR,rebal_stat.MPI_SOURCE,
		       rebal_stat.MPI_TAG,comm->mm_inter,&rebal_req);
      err += MPI_Wait(&rebal_req,MPI_STATUS_IGNORE);
    }

    /* build work_queue for current step */

    /* spawn rebalancing on pgf_mpi_comm->micro */

    /* rebalance master(s) on pgf_mpi_comm->worker_inter */

    /* Complete rebalancing -- mark as NEED_INFO */

    free(buffer);
  }

  /* post receives for jobs to compute (mm_inter)*/

  while(1 /* !work_queue_is_empty(queue) */){

    /* MPI_Waitsome for information. */

    /* Mark jobs in work_queue as READY */

    /* Compute READY jobs in work_queue */

    /* Post response and mark REPLY */

    /* MPI_Testsome for finished replies. Mark DONE */

  }

  return err;
}

enum {REBAL_N_KEEP=0,
      REBAL_N_SEND,
      REBAL_N_RECV,
      REBAL_KEEP_OFF,
      REBAL_SEND_OFF,
      REBAL_RECV_OFF,
      REBAL_N_META};
void pgf_FE2_server_rebalance_build(pgf_FE2_server_rebalance *t,
				    const size_t n_keep,
				    const size_t n_send,
				    const size_t n_recv)
{
  /* currently, pgf_FE2_server_rebalance is simply a typedef for a
     pointer to an integer, but we will have a special layout of the
     data described in the code below. */

  /* allocate buffer */
  *t = calloc(REBAL_N_META + n_keep + 2*(n_send + n_recv),sizeof(**t));

  /* encode separate buffer sizes */
  (*t)[REBAL_N_KEEP] = n_keep;
  (*t)[REBAL_N_SEND] = n_send;
  (*t)[REBAL_N_RECV] = n_recv;

  /* encode offsets */
  (*t)[REBAL_KEEP_OFF] = REBAL_N_META; /* keep_offset */
  (*t)[REBAL_SEND_OFF] = (*t)[3] + n_keep; /* send_offset */
  (*t)[REBAL_RECV_OFF] = (*t)[4] + 2*n_send; /* recv_offset */
}

void pgf_FE2_server_rebalance_destroy(pgf_FE2_server_rebalance *t)
{
  free(*t);
}

int pgf_FE2_server_rebalance_n_keep(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t)[REBAL_N_KEEP];
  else return -1;
}

int* pgf_FE2_server_rebalance_keep_buff(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return *t + (*t)[REBAL_KEEP_OFF];
  else return NULL;
}

int pgf_FE2_server_rebalance_n_send(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t)[REBAL_N_SEND];
  else return -1;
}

int* pgf_FE2_server_rebalance_send_buff(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t) + (*t)[REBAL_SEND_OFF];
  else return NULL;
}

int* pgf_FE2_server_rebalance_send_dest(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t) + (*t)[REBAL_SEND_OFF] + (*t)[REBAL_N_SEND];
  else return NULL;
}

int pgf_FE2_server_rebalance_n_recv(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t)[REBAL_N_RECV];
  else return -1;
}

int* pgf_FE2_server_rebalance_recv_buff(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t) + (*t)[REBAL_RECV_OFF];
  else return NULL;
}

int* pgf_FE2_server_rebalance_recv_src(const pgf_FE2_server_rebalance *t)
{
  if(*t != NULL) return (*t) + (*t)[REBAL_RECV_OFF] + (*t)[REBAL_N_RECV];
  else return NULL;
}
