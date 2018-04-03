/* HEADER */
#pragma once
#ifndef COMPUTE_MS_COHE_JOB_H
#define COMPUTE_MS_COHE_JOB_H

#include "microscale_information.h"
#include "ms_cohe_job_info.h"

/** Compute the microscale solution, tangents, etc. for a single
    microscale cohesive job. */
int compute_ms_cohe_job(const int job_id,
                        MS_COHE_JOB_INFO *p_job,
                        MICROSCALE *microscale,
                        const int mp_id,
                        int &EXA_metric);

/** Assemble the macroscale cohesive residual to the local vector on
    the owning domain. NO COMMUNICATION */
int assemble_ms_cohe_job_res(const int job_id,
                             const MS_COHE_JOB_INFO *p_job,
                             const MPI_Comm micro_comm,
                             const MPI_Comm macro_comm,
                             double *loc_res,
                             const int mp_id);

#endif
