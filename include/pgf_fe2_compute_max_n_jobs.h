/* HEADER */
/**
 * AUTHORS:
 *  Matthew Mosby, University of Notre Dame <mmosby1[at]nd.edu>
 */
#pragma once
#ifndef PGF_FE2_COMPUTE_MAX_N_JOBS_H
#define PGF_FE2_COMPUTE_MAX_N_JOBS_H

#include "microscale_information.h"
#include "PGFEM_mpi.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */
/**  
 * Compute the maximum number of jobs to allow on a single server.
 * This function must be called by ALL processes simultaneously.
 *
 * The maximum number of jobs is `max_n_jobs = 1.1 * min_n_jobs`,
 * where `min_n_jobs = total_macro_nce / n_servers + rem` and 
 * `rem = 1` if `total_macro_nce % n_servers > 0`. The maximum
 * number of jobs may also be specified at runtime through the
 * command line (see -h) so long as it is greater than min_n_jobs.
 *
 * Communication on mpi_comm->world.
 * \param[in] macro Macroscale domain (=NULL for servers).
 * \param[in] mpi_comm FE2 MPI communicators
 * \param[out] max_n_jobs
 * \return non-zero on internal error
 */
int pgf_FE2_compute_max_n_jobs(const MACROSCALE *macro,
			       const PGFEM_mpi_comm *mpi_comm,
			       int *max_n_jobs);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */


#endif /* #ifndef  */
