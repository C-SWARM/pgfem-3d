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

enum pgf_FE2_micro_server_event_tags{
  FE2_MICRO_SERVER_EXIT = -5,       /**< Signal the micro server to
				       cleanup and exit */

  FE2_MICRO_SERVER_REBALANCE = -1   /**< Signal the micro server to
				       prepare for a new list of work
				       and possible rebalancing
				       step. */
};

/**
 * Structure for maintaining overall statistics on the server's
 * current work load.
 */
struct pgf_FE2_micro_server_stats{
  double total;
  double avg;
  double std;
  double min;
  double max;
};
typedef struct pgf_FE2_micro_server_stats pgf_FE2_micro_server_stats;

/**
 * Sturucture to organize a server's work.
 */
struct pgf_FE2_micro_server{
  size_t n_jobs;
  pgf_FE2_job *jobs;
  pgf_FE2_micro_server_stats *stats;
};
typedef struct pgf_FE2_micro_server pgf_FE2_micro_server;

#endif
