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
#include "microscale_information.h"

/**
 * Structure for maintaining communication information related to a
 * particular job.
 */
struct pgf_FE2_job_comm_buf{
  size_t buffer_len;                       /**< length of the buffer in bytes */
  char *buffer;
  MPI_Request request;
  MPI_Status status;
};

/**
 * Set a valid initial state for a pgf_FE2_job_comm_buf object.
 *
 * Performs no allocation.
 */
void pgf_FE2_job_comm_buf_init(pgf_FE2_job_comm_buf *buf);

/**
 * Allocate internal structure for a pgf_FE2_job_comm_buf object.
 *
 * Assumes that internal structures are not already allocated. Calling
 * pgf_FE2_job_comm_buf_build on multiple times on an object without
 * intermediat calls to pgf_FE2_job_comm_buf_destroy results in a
 * memory leak.
 */
void pgf_FE2_job_comm_buf_build(pgf_FE2_job_comm_buf *buf,
                const size_t buffer_len);

/**
 * Destroys pgf_FE2_job_comm_buf object and leaves it in a valid
 * (empty) state.
 */
void pgf_FE2_job_comm_buf_destroy(pgf_FE2_job_comm_buf *buf);

/**
 * Enumeration for the progression of job states.
 */
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

/**
 * Sets/allocates an initial valid state for a pgf_FE2_job object.
 */
void pgf_FE2_job_init(pgf_FE2_job *job,
              const int id,
              const int state);

/**
 * Set the job progression state for a pgf_FE2_job object.
 */
void pgf_FE2_job_set_state(pgf_FE2_job *job,
               const int state);

/**
 * Destroy a pgf_FE2_job object.
 */
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

/**
 * Comparator function for 'state' of two pgf_FE2_job
 * objects. Suitable for use in sorting functions etc.
 */
int pgf_FE2_job_compare_state(const void *a,
                  const void *b);

/**
 * Comparator function for 'time' of two pgf_FE2_job objects. Suitable
 * for use in sorting functions etc.
 */
int pgf_FE2_job_compare_time(const void *a,
                 const void *b);

/**
 * Comparator function for 'id' of two pgf_FE2_job objects. Suitable
 * for use in sorting functions etc.
 */
int pgf_FE2_job_compare_id(const void *a,
               const void *b);

/**
 * Check the job state and attempt to get information from macroscale
 * if needed. Returns job state on exit.
 */
int pgf_FE2_job_get_info(pgf_FE2_job *job,
             const PGFEM_mpi_comm *mpi_comm);

/**
 * Check the job state and compute if possible. Returns job state on
 * exit.
 */
int pgf_FE2_job_compute(pgf_FE2_job *job,
            MICROSCALE *micro,
            const PGFEM_mpi_comm *mpi_comm,
            const int mp_id);

/**
 * Initiates computing a job on a worker process (non-master).
 */
void pgf_FE2_job_compute_worker(const size_t job_id,
                const size_t buffer_len,
                MICROSCALE *micro,
                const int mp_id);

/**
 * Check the job state and reply to the macroscale if
 * possible. Returns job state on exit.*/
int pgf_FE2_job_reply(pgf_FE2_job *job,
              const PGFEM_mpi_comm *mpi_comm);

/**
 * Check the job state and complete if possible. Returns job state on
 * exit.
 */
int pgf_FE2_job_complete(pgf_FE2_job *job);

#endif
