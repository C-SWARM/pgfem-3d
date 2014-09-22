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
void pgf_FE2_rebalancer(const PGFEM_mpi_comm *mpi_comm,
			const size_t total_n_jobs,
			const size_t max_n_jobs);
#endif
