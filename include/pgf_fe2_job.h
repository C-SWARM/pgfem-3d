/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */
#pragma once
#ifndef PGF_FE2_JOB_H
#define PGF_FE2_JOB_H

#include <stdlib.h>

/**
 * Structure to maintain the job status in terms of progress and time
 * to completion.
 */
struct pgf_FE2_job_status{
  double time;
  int have_info;
  int state;
};
#ifndef TYPEDEF_pgf_FE2_job_status
#define TYPEDEF_pgf_FE2_job_status
typedef struct pgf_FE2_job_status pgf_FE2_job_status;
#endif

void pgf_FE2_job_status_init(pgf_FE2_job_status *status);
void pgf_FE2_job_status_destroy(pgf_FE2_job_status *status);

int pgf_FE2_job_status_compare_state(const void *a,
				     const void *b);

int pgf_FE2_job_status_compare_time(const void *a,
				    const void *b);

/**
 * Structure to manage data related to a particular job.
 */
struct pgf_FE2_job_data{
  size_t len;
  size_t max_len;
  char *buffer;
};
#ifndef TYPEDEF_pgf_FE2_job_data
#define TYPEDEF_pgf_FE2_job_data
typedef struct pgf_FE2_job_data pgf_FE2_job_data;
#endif

/**
 * Initialize the object.
 */
void pgf_FE2_job_data_init(pgf_FE2_job_data *data);

/**
 * Reset the working size for the data. Calls realloc as needed. Does
 * not ever shrink the amount of memory allocated to the buffer.
 */
void pgf_FE2_job_data_resize(pgf_FE2_job_data *data,
			     const size_t len);

/**
 * Destroy the data object.
 */
void pgf_FE2_job_data_destroy(pgf_FE2_job_data *data);

/**
 * structure to maintain a FE2 microscale job.
 */
struct pgf_FE2_job{
  size_t id;
  size_t micro_sol_idx;
  pgf_FE2_job_status status;
  pgf_FE2_job_data data;
};
#ifndef TYPEDEF_pgf_FE2_job
#define TYPEDEF_pgf_FE2_job
typedef struct pgf_FE2_job pgf_FE2_job;
#endif

void pgf_FE2_job_init(pgf_FE2_job *job);
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

#endif
