/**
 * @file This hearder declares the functions called at the macroscale
 * to compute the microscale contributions from a list of jobs. These
 * functions serve as the primary interface to the microscale from the
 * macroscale.
 */
#ifndef PGFEM3D_MS_COHE_JOB_LIST_H
#define PGFEM3D_MS_COHE_JOB_LIST_H

#include "PGFEM_mpi.h"
#include "cohesive_element.h"
#include "microscale_information.h"
#include "ms_cohe_job_info.h"
#include "node.h"
#include "pgfem3d/Solver.hpp"

/** Create the list of jobs to be performed on the communicator
    'ms_comm'. Returns the number of jobs on 'ms_comm' (Gnjobs), the
    number of jobs for each process on 'ms_comm' (n_job_dom), and
    the allocated list of jobs. */
int create_group_ms_cohe_job_list(const long nce,
                                  const COEL *coel,
                                  const Node *node,
                                  const MPI_Comm macro_mpi_comm,
                                  const MPI_Comm ms_comm,
                                  const int group_id,
                                  long *Gn_jobs,
                                  long **n_job_dom,
                                  MS_COHE_JOB_INFO **job_list,
                                  const int mp_id);

/** Update the displacement jump for each job in the list. This is
    the ONLY thing that is updated */
int update_group_ms_cohe_job_list(const long nce,
                                  const COEL *coel,
                                  const Node *node,
                                  const SUPP sup,
                                  const double *sol,
                                  MPI_Comm ms_comm,
                                  MS_COHE_JOB_INFO *job_list);

/** Compute the microscale contributions from the list of jobs. The
    contributions to the macroscale tangent are assembled on the
    processes that own the job. Contains communication at macro and
    micro scales */
int compute_ms_cohe_tan_res(const int compute_micro_eq,
                            const COMMUN macro_pgfem_comm,
                            const MPI_Comm macro_mpi_comm,
                            MS_COHE_JOB_INFO *job_list,
                            pgfem3d::solvers::SparseSystem *macro_solver,
                            MICROSCALE *microscale,
                            const int mp_id);

/** Assemble the residuals computed from the last call to *insert
    function name*. No communication */
int assemble_ms_cohe_res(const MICROSCALE *micro,
                         const MS_COHE_JOB_INFO *jobs,
                         const MPI_Comm macro_mpi_comm,
                         double *macro_loc_res,
                         const int mp_id);

/** destroy the list of jobs */
void destroy_ms_cohe_job_list(const long Gn_job,
                              MS_COHE_JOB_INFO *job_list);

#endif /* #define PGFEM3D_MS_COHE_JOB_LIST_H  */
