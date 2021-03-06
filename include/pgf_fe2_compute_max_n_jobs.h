/* HEADER */
/**
 * AUTHORS:
 *  Matthew Mosby, University of Notre Dame <mmosby1[at]nd.edu>
 */
#pragma once
#ifndef PGF_FE2_COMPUTE_MAX_N_JOBS_H
#define PGF_FE2_COMPUTE_MAX_N_JOBS_H

#include "pgfem3d/MultiscaleCommon.hpp"

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
 * \param[in] micro Microscale domain (=NULL for macro)
 * \param[in] mscom Multiscale communicator structure
 * \param[out] max_n_jobs
 * \return non-zero on internal error
 */
int pgf_FE2_compute_max_n_jobs(const pgfem3d::Macroscale *macro,
			       const pgfem3d::Microscale *micro,
			       const multiscale::MultiscaleCommunicator *mscom,
			       int *max_n_jobs);

#endif /* #ifndef  */
