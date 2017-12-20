/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#pragma once
#ifndef MICRO_SERVER_H
#define MICRO_SERVER_H

#include <stdlib.h>
#include "microscale_information.h"
#include "pgf_fe2_job.h"
#include "pgf_fe2_server_rebalance.h"

/**
 * Specify the minimum tag upper bound prescribed in the MPI
 * standard. It seems that querying MPI_TAG_UB may not be accurate for
 * the implementation... (for example mvapich2-2.0: MPI_TAG_UB - 1
 * gives MPI_TAG_ERR)
 */
#define PGFEM_TAG_UB 32767
/**
 * Enumeration to use for special macro-micro communication tags.
 */
enum pgf_FE2_micro_server_event_tags{
  /** Signal the micro server to cleanup and exit */
  FE2_MICRO_SERVER_EXIT = PGFEM_TAG_UB - 1,

  /** Signal the micro server to prepare for a new list of work and
   * possible rebalancing step.
   */
  FE2_MICRO_SERVER_REBALANCE = PGFEM_TAG_UB -2
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

/**
 * Sturucture to organize a server's work.
 */
struct pgf_FE2_micro_server{
  size_t n_jobs;
  pgf_FE2_job *jobs;
  pgf_FE2_micro_server_stats *stats;
};

/**
 * Initialize a handle to a pgf_FE2_micro_server object. (allow for
 * full encapsulation)
 */
void pgf_FE2_micro_server_init(pgf_FE2_micro_server **server);

/**
 * Allocate and populate a pgf_FE2_micro_server object from data in
 * the pgf_FE2_server_rebalance object.
 */
void pgf_FE2_micro_server_build(pgf_FE2_micro_server *server,
                const pgf_FE2_server_rebalance *rebal);

/**
 * Destroy the pgf_FE2_micro_server object. Points the handle to NULL.
 */
void pgf_FE2_micro_server_destroy(pgf_FE2_micro_server *server);

/**
 * This is the main function that starts the master/worker server
 * processes. This function must be called by ALL micro-scale
 * processes.
 */
int pgf_FE2_micro_server_START(const pgfem3d::MultiscaleComm *mscom,
			       MICROSCALE *micro,
			       const int mp_id);

/**
 * Unpack a server summary (shallow copy) from a buffer and return a
 * handle to the resulting constructed pgf_FE2_micro_server object.
 *
 * Note that the shallow copy omits some internal data, i.e., the
 * communication buffers for each job. These structures are in a valid
 * state such that they may be built (if desired) or cleanly
 * destroyed.
 */
void pgf_FE2_micro_server_unpack_summary(pgf_FE2_micro_server **server,
                     const char *buf);

/**
 * Print the stats structure nicely so it can be grepped from the logs.
 *
 * The printed line(s) begins with "SERVER STATS"
 */
void pgf_FE2_micro_server_stats_print(const pgf_FE2_micro_server *server);

#endif
