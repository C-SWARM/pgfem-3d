/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#pragma once
#ifndef MICRO_SERVER_H
#define MICRO_SERVER_H

#include <stdlib.h>
#include "mpi.h"
//#include "PGFEM_mpi.h" /* include after prelim testing */

/**
 * Uniquely identify a microscale job.
 */
struct pgf_job_id{
  size_t proc_id;
  size_t el_id;
  size_t int_pt;
};
#ifndef TYPEDEF_pgf_job_id
#define TYPEDEF_pgf_job_id
typedef struct pgf_job_id pgf_job_id;
#endif

/**
 * Micro server job object.
 */ 
struct pgf_micro_server_job{
  pgf_job_id job_id;
  size_t have_info;
  int status;
  char *buf;
};
#ifndef TYPEDEF_pgf_micro_server_job
#define TYPEDEF_pgf_micro_server_job
typedef struct pgf_micro_server_job pgf_micro_server_job;
#endif

/**
 * Object to manage list of jobs.
 */
struct pgf_server_job_pool{
  pgf_micro_server_job *jobs;
  size_t n_jobs;
  size_t max_n_jobs;
};
#ifndef TYPEDEF_pgf_server_job_pool
#define TYPEDEF_pgf_server_job_pool
typedef struct pgf_server_job_pool pgf_server_job_pool;
#endif

/**
 * Structure to identify and monitor the transfer of data among
 * microscale servers.
 */
struct pgf_worker_transfer_status{
  size_t comm_type; /** send or receive */
  MPI_Request *request;
  MPI_Status *status;
};

#endif
