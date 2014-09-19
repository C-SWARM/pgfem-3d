/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "pgf_fe2_job.h"
#include <limits.h>
#include <assert.h>

/*** FE2_job_comm_buf ***/
void pgf_FE2_job_comm_buf_init(pgf_FE2_job_comm_buf *buf)
{
  buf->buffer_len = 0;
  buf->buffer = NULL;
}

void pgf_FE2_job_comm_buf_build(pgf_FE2_job_comm_buf *buf,
				const size_t buffer_len)
{
  buf->buffer_len = buffer_len;
  buf->buffer = malloc(buffer_len*sizeof(*(buf->buffer)));
}

void pgf_FE2_job_comm_buf_destroy(pgf_FE2_job_comm_buf *buf)
{
  free(buf->buffer);
}

/*** FE2_job ***/

static inline int pgf_FE2_job_private_compare_int(const int *a,
						  const int *b)
{
  return (*a < *b) - (*a > *b);
}

static inline int pgf_FE2_job_private_compare_size_t(const size_t *a,
						     const size_t *b)
{
  return (int) ((*a < *b) - (*a > *b));
}

int pgf_FE2_job_compare_state(const void *a,
			      const void *b)
{
  const pgf_FE2_job *A = a;
  const pgf_FE2_job *B = b;
  return pgf_FE2_job_private_compare_int(&(A->state),&(B->state));
}

int pgf_FE2_job_compare_time(const void *a,
			     const void *b)
{
  const pgf_FE2_job *A = a;
  const pgf_FE2_job *B = b;
  return pgf_FE2_job_private_compare_size_t(&(A->time),&(B->time));
}


int pgf_FE2_job_compare_id(const void *a,
			   const void *b)
{
  const pgf_FE2_job *A = a;
  const pgf_FE2_job *B = b;
  return pgf_FE2_job_private_compare_size_t(&(A->id),&(B->id));
}

void pgf_FE2_job_init(pgf_FE2_job *job,
		      const int id,
		      const int state)
{
  /* poison values */
  job->id = id;
  job->time = 0;
  job->state = state;
  job->comm_buf = malloc(sizeof(*(job->comm_buf)));
  pgf_FE2_job_comm_buf_init(job->comm_buf);
}

void pgf_FE2_job_destroy(pgf_FE2_job *job)
{
  pgf_FE2_job_comm_buf_destroy(job->comm_buf);
  free(job->comm_buf);
  job->comm_buf = NULL;
}

static const int encode_proc_offset = 1e6;
static const int encode_elem_offset = 1e2;
int pgf_FE2_job_compute_encoded_tag(const size_t proc_id,
				    const size_t elem_id,
				    const size_t int_pt)
{
  assert(proc_id <= INT_MAX/encode_proc_offset);
  assert(elem_id < encode_proc_offset);
  assert(int_pt < encode_elem_offset);
 
  return (proc_id * encode_proc_offset
	  + elem_id * encode_elem_offset
	  + int_pt);
}

void pgf_FE2_job_decode_id(const int id,
			   size_t *proc_id,
			   size_t *elem_id,
			   size_t *int_pt)
{
  assert(id >= 0);
  *proc_id = pgf_FE2_job_decode_id_proc(id);
  *elem_id = id % encode_proc_offset;
  *int_pt = *elem_id % encode_elem_offset;
  *elem_id /= encode_elem_offset;
}

int pgf_FE2_job_decode_id_proc(const int id)
{
  assert(id >= 0);
  return id / encode_proc_offset;
}

/**
 * Check the status level and attempt to send or receive information
 * from macroscale if needed. Returns job state on exit.
 */
int pgf_FE2_job_state_check_info(pgf_FE2_job *job,
				 PGFEM_mpi_comm *mpi_comm)
{
  /* Already have the info to compute, return state */
  if(job->state != FE2_STATE_NEED_INFO_REBALANCE 
     && job->state != FE2_STATE_NEED_INFO) return job->state;

  /* probe for message based on job id */
  size_t proc = 0;
  size_t elem = 0;
  size_t int_pt = 0;
  MPI_Status status;
  int msg_waiting = 0;
  pgf_FE2_job_decode_id(job->id,&proc,&elem,&int_pt);
  MPI_Iprobe(proc,job->id,mpi_comm->mm_inter,&msg_waiting,&status);

  /* if there is a message, allocate the comm_buf and post the
     matching receive. */
  if(msg_waiting){
    int len = 0;
    MPI_Get_count(&status,MPI_CHAR,&len);
    pgf_FE2_job_comm_buf_build(job->comm_buf,len);

    /* post blocking receive since we can't move on until we get it
       and the message is for sure there. */
    MPI_Recv(job->comm_buf->buffer,job->comm_buf->buffer_len,
	     MPI_CHAR,proc,job->id,mpi_comm->mm_inter,&status);

    /* update the job state */
    switch(job->state){
    case FE2_STATE_NEED_INFO_REBALANCE:
      job->state = FE2_STATE_HAVE_INFO_REBALANCE;
      break;
    case FE2_STATE_NEED_INFO:
      job->state = FE2_STATE_COMPUTE_READY;
      break;
    default:  /* should never get here */ assert(0); break;
    }
  }
  return job->state;
}
