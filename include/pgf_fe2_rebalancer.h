/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */
#pragma once
#ifndef PGF_FE2_REBALANCER_H
#define PGF_FE2_REBALANCER_H

#include <stdlib.h>
#include "PGFEM_mpi.h"
#include "pgf_fe2_server_rebalance.h"

/**
 * Allows for switching of different rebalancing heuristics. May
 * incorporate into API in the future.
 */
enum pgf_FE2_rebalancer_heuristic {
  FE2_REBALANCE_NONE,
  FE2_REBALANCE_GREEDY
};

/**
 * Load balancer called on each macroscale prcocess to rebalance the
 * servers. Nominally collective communication on mpi_comm->mm_inter.
 */
pgf_FE2_server_rebalance** pgf_FE2_rebalancer(const PGFEM_mpi_comm *mpi_comm,
					     const size_t total_n_jobs,
					     const size_t max_n_jobs,
					     const int heuristic);

/* The following is REALLY ugly, but I don't want to mess with it
right now, I just want it to work and be externally callable. Will fix
up later (maybe) */

/**
 * Very ugly initialization of data structure I don't want to keep around...
 */
void new_partition_build_set_keep(void **np,
				  const size_t max_n_job,
				  const size_t n_job,
				  const int *job_id,
				  const int *job_time,
				  const int *proc_id);

void new_partitions_void(void **parts,
			 const size_t n_parts,
			 const size_t n_max_job);

/**
 * destroy via opaque handle
 */
void new_partition_destroy_void(void *np);

/**
 * Destroy a collection via opaque handle
 */
void new_partitions_destroy_void(void *np,
				 const size_t n_parts);

/**
 * Greedy algorithm for partitioning the jobs.
 */
void rebalance_partitions_greedy(const size_t n_parts,
				 void *all_parts, /* one object containing all jobs */
				 void *parts); /* pointer to multiple objects */

/**
 * Data format conversion to sturcture I actually plan to keep
 * around/use externally.
 */
void new_partitions_void_to_pgf_FE2_server_rebalance(const int n_parts,
						     const void *np,
						     pgf_FE2_server_rebalance **rb);
#endif
