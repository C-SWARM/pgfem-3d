/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "pgf_fe2_job.h"
#include "macro_micro_functions.h"
#include <time.h>
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
  buf->buffer = NULL;
  buf->buffer_len = 0;
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

void pgf_FE2_job_set_state(pgf_FE2_job *job,
			   const int state)
{
  job->state = state;
}


void pgf_FE2_job_destroy(pgf_FE2_job *job)
{
  if(job->comm_buf != NULL){
    pgf_FE2_job_comm_buf_destroy(job->comm_buf);
  }
  free(job->comm_buf);
  job->comm_buf = NULL;
}

static const int encode_proc_offset = 1e6;
static const int encode_elem_offset = 1e2;
int pgf_FE2_job_compute_encoded_id(const size_t proc_id,
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

int pgf_FE2_job_get_info(pgf_FE2_job *job,
			 const PGFEM_mpi_comm *mpi_comm)
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

void pgf_FE2_job_compute_worker(const size_t job_id,
				const size_t buffer_len,
				MICROSCALE *micro)
{
  /* allocate and receive job buffer */
  char *buf = malloc(buffer_len*sizeof(*buf));
  MPI_Bcast(buf,buffer_len,MPI_CHAR,0,micro->common->mpi_comm);

  /* determine solution index from job id and compute */
  int idx = sol_idx_map_id_get_idx(&(micro->idx_map),job_id);
  int exit_server = 0; /* unused in this implementation */
  microscale_compute_job(idx,buffer_len,buf,micro,&exit_server);
  assert(exit_server == 0);
  free(buf);
}

int pgf_FE2_job_compute(pgf_FE2_job *job,
			MICROSCALE *micro,
			const PGFEM_mpi_comm *mpi_comm)
{
  /* return immediately if not ready to compute */
  if(job->state != FE2_STATE_COMPUTE_READY) return job->state;

  /* record time of job */
  time_t start, finish;
  time(&start);

  /* broadcast information to the microscale */
  assert(mpi_comm->rank_micro == 0);
  static const int n_meta = 2;
  size_t meta[n_meta];
  meta[0] = job->id;
  meta[1] = job->comm_buf->buffer_len;
  MPI_Bcast(meta,n_meta,MPI_LONG,0,micro->common->mpi_comm);
  MPI_Bcast(job->comm_buf->buffer,meta[1],MPI_CHAR,0,
	    micro->common->mpi_comm);

  /* determine solution index from job id and compute */
  int idx = sol_idx_map_id_get_idx(&(micro->idx_map),meta[0]);
  assert(idx >= 0);
  int exit_server = 0; /* unused in this implementation */
  microscale_compute_job(idx,meta[1],job->comm_buf->buffer,
			 micro,&exit_server);
  assert(exit_server == 0);

  /* increment state and post communication to macroscale */
  job->state = FE2_STATE_REPLY_READY;

  /* time to complete job */
  time(&finish);
  job->time = difftime(finish,start) + 1;

  return pgf_FE2_job_reply(job,mpi_comm);
}

int pgf_FE2_job_reply(pgf_FE2_job *job,
		      const PGFEM_mpi_comm *mpi_comm)
{
  /* return immediately if not ready to reply */
  if(job->state != FE2_STATE_REPLY_READY) return job->state;

  /* post non-blocking send of information to the macroscale */
  size_t proc = 0, elem = 0, int_pt = 0;
  pgf_FE2_job_decode_id(job->id,&proc,&elem,&int_pt);
  MPI_Isend(job->comm_buf->buffer,job->comm_buf->buffer_len,
	    MPI_CHAR,proc,job->id,mpi_comm->mm_inter,
	    &(job->comm_buf->request));

  job->state = FE2_STATE_REPLY;
  return job->state;
}

int pgf_FE2_job_complete(pgf_FE2_job *job)
{
  /* return immediately if reply not in progress */
  if(job->state != FE2_STATE_REPLY) return job->state;

  int finished = 0;
  MPI_Test(&(job->comm_buf->request),&finished,MPI_STATUS_IGNORE);
  if(finished){
    /* destroy the comm_buffer rendering it unusable for erroneous
       further communication. */
    pgf_FE2_job_comm_buf_destroy(job->comm_buf);
    job->state = FE2_STATE_DONE;
  }

  return job->state;
}
