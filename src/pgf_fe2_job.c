/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "pgf_fe2_job.h"
#include <limits.h>
#include <assert.h>

void pgf_FE2_job_id_init(pgf_FE2_job_id *id)
{
  /* poison values */
  id->tag = INT_MIN;
  id->micro_sol_idx = -1;
}

void pgf_FE2_job_id_destroy(pgf_FE2_job_id *id)
{
  /* do nothing */
}

static const int encode_proc_offset = 1e6;
static const int encode_elem_offset = 1e2;

void pgf_FE2_job_id_set(pgf_FE2_job_id *id,
			const size_t proc_id,
			const size_t elem_id,
			const size_t int_pt)
{
  id->tag = compute_pgf_FE2_encoded_tag(proc_id,elem_id,int_pt);
}

void pgf_FE2_job_id_get_info(const pgf_FE2_job_id *id,
			     size_t *proc_id,
			     size_t *elem_id,
			     size_t *int_pt)
{
  assert(id->tag >= 0);
  *proc_id = id->tag / encode_proc_offset;
  *elem_id = id->tag % encode_proc_offset;
  *int_pt = *elem_id % encode_elem_offset;
  *elem_id /= encode_elem_offset;
}

int compute_pgf_FE2_encoded_tag(const size_t proc_id,
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

enum pgf_FE2_job_status_state{
  FE2_STATUS_UNDEFINED=-1,
  FE2_STATUS_READY,
  FE2_STATUS_REBALANCE,
  FE2_STATUS_NEED_INFO,
  FE2_STATUS_REPLY,
  FE2_STATUS_DONE};

void pgf_FE2_job_status_init(pgf_FE2_job_status *status)
{
  status->time = 0.0;
  status->have_info = 0;
  status->state = FE2_STATUS_UNDEFINED;
}

void pgf_FE2_job_status_destroy(pgf_FE2_job_status *status)
{
  /* do nothing */
}

int pgf_FE2_job_status_compare_state(const void *a,
				     const void *b)
{
  const pgf_FE2_job_status *A = a;
  const pgf_FE2_job_status *B = b;
  return A->state - B->state;
}

int pgf_FE2_job_status_compare_time(const void *a,
				    const void *b)
{
  const pgf_FE2_job_status *A = a;
  const pgf_FE2_job_status *B = b;
  return (A->time > B->time) - (A->time < B->time);
}

void pgf_FE2_job_data_init(pgf_FE2_job_data *data)
{
  data->len = 0;
  data->max_len = 0;
  data->buffer = NULL;
}

void pgf_FE2_job_data_resize(pgf_FE2_job_data *data,
			     const size_t len)
{
  if(len > data->max_len){
    data->buffer = realloc(data->buffer,len*sizeof(*(data->buffer)));
    data->max_len = len;
  } else {
    data->len = len;
  }
}

void pgf_FE2_job_data_destroy(pgf_FE2_job_data *data)
{
  free(data->buffer);
}

void pgf_FE2_job_init(pgf_FE2_job *job)
{
  pgf_FE2_job_id_init(&(job->id));
  pgf_FE2_job_status_init(&(job->status));
  pgf_FE2_job_data_init(&(job->data));
}

void pgf_FE2_job_destroy(pgf_FE2_job *job)
{
  pgf_FE2_job_id_destroy(&(job->id));
  pgf_FE2_job_status_destroy(&(job->status));
  pgf_FE2_job_data_destroy(&(job->data));
}
