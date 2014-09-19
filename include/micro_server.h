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
#include "pgf_fe2_server_rebalance.h"

struct pgf_FE2_micro_server_stats{
  double total;
  double avg;
  double std;
  double min;
  double max;
};
typedef struct pgf_FE2_micro_server_stats pgf_FE2_micro_server_stats;

struct pgf_FE2_micro_server{
  size_t n_jobs;
  pgf_FE2_job *jobs;
  pgf_FE2_micro_server_stats stats;
};
typedef struct pgf_FE2_micro_server pgf_FE2_micro_server;

#endif
