/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */
#pragma once
#ifndef PGF_FE2_JOB_H
#define PGF_FE2_JOB_H

#include <stdlib.h>
#include "PGFEM_mpi.h"

struct pgf_FE2_job_comm_buf{
  size_t buffer_len;
  char *buffer;
  MPI_Request request;
  MPI_Status status;
};
#ifndef TYPEDEF_pgf_FE2_job_comm_buf
#define TYPEDEF_pgf_FE2_job_comm_buf
typedef struct pgf_FE2_job_comm_buf pgf_FE2_job_comm_buf;
#endif

void pgf_FE2_job_comm_buf_init(pgf_FE2_job_comm_buf *buf);
void pgf_FE2_job_comm_buf_build(pgf_FE2_job_comm_buf *buf,
				const size_t buffer_len);
void pgf_FE2_job_comm_buf_destroy(pgf_FE2_job_comm_buf *buf);

enum pgf_FE2_job_state{
  FE2_STATE_UNDEFINED=-1,
  FE2_STATE_NEED_INFO_REBALANCE,  /*< job is rebalancing */
  FE2_STATE_NEED_INFO,            /*< job needs info to compute (rebal done) */
  FE2_STATE_HAVE_INFO_REBALANCE,  /*< have info to compute when rebal done */
  FE2_STATE_COMPUTE_READY,        /*< job is ready to compute */
  FE2_STATE_REPLY_READY,          /*< job computed and ready to reply */
  FE2_STATE_REPLY,                /*< reply comm posted */
  FE2_STATE_DONE                  /*< reply completed, all done */
};               

/**
 * Structure containing an FE2 microscale job metadata and
 * communication buffer(s).
 *
 * The job contains an encoded id, a state counter for tracking
 * progress of the job, a time variable for the length of time spent
 * comupting the job with second resolution.
 */
struct pgf_FE2_job{
  size_t id;
  size_t time;
  pgf_FE2_job_comm_buf *comm_buf;
  int state;
};
#ifndef TYPEDEF_pgf_FE2_job
#define TYPEDEF_pgf_FE2_job
typedef struct pgf_FE2_job pgf_FE2_job;
#endif

void pgf_FE2_job_init(pgf_FE2_job *job,
		      const int id,
		      const int state);
void pgf_FE2_job_destroy(pgf_FE2_job *job);

/**
 * Compute an encoded integer id from the macroscale proc_id, elem_id
 * and element int_pt.
 *
 * Restrictions: proc_id <= INT_MAX/1e6, elem_id < 1e6, int_pt < 1e2
 * \return encoded integer id.
 */
int pgf_FE2_job_compute_encoded_id(const size_t proc_id,
				   const size_t elem_id,
				   const size_t int_pt);

/**
 * Decode an integer id and return the macroscale information.
 */
void pgf_FE2_job_decode_id(const int id,
			   size_t *proc_id,
			   size_t *elem_id,
			   size_t *int_pt);

/**
 * Extract only the macroscale proc_id from an encoded integer id.
 */
int pgf_FE2_job_decode_id_proc(const int id);

int pgf_FE2_job_compare_state(const void *a,
			      const void *b);

int pgf_FE2_job_compare_time(const void *a,
			     const void *b);

int pgf_FE2_job_compare_id(const void *a,
			   const void *b);


#endif
