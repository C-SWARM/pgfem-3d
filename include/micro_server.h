/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#pragma once
#ifndef MICRO_SERVER_H
#define MICRO_SERVER_H

#include <stdlib.h>
#include "PGFEM_mpi.h"
#include "pgf_fe2_job.h"

struct pgf_FE2_server_stats{
  double avg;
  double std;
  double min;
  double max;
};
typedef struct pgf_FE2_server_stats pgf_FE2_server_stats;

struct pgf_FE2_server_job_list{
  size_t n_jobs;
  pgf_FE2_job *jobs;
  pgf_FE2_server_stats stats;
};
typedef struct pgf_FE2_server_job_list pgf_FE2_server_job_list;

#endif
