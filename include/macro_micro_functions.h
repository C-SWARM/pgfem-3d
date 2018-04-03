/* HEADER */
/**
 * AUTHORS:
 *  Matthew Mosby, University of Notre Dame
 */
#pragma once
#ifndef MACRO_MICRO_FUNCTIONS_H
#define MACRO_MICRO_FUNCTIONS_H

#include "microscale_information.h"
#include "ms_cohe_job_info.h"
#include "ms_job_intercomm.h"

struct pgf_FE2_macro_client;
struct PGFEM_mpi_comm;

/* container for passing through Newton Raphson */
struct MS_SERVER_CTX{
  struct pgf_FE2_macro_client *client;
  struct PGFEM_mpi_comm *mpi_comm;
  MACROSCALE *macro;
};

/** */
int compute_n_job_and_job_sizes(const COMMON_MACROSCALE *c,
                                int *n_jobs,
                                int **job_buff_sizes);

/** Start a microscale server process. Non-blocking collective
    communication on mpi_comm->mm_inter and mpi_comm->micro */
int start_microscale_server(const PGFEM_mpi_comm *mpi_comm,
                            const PGFEM_ms_job_intercomm *ic,
                            MICROSCALE *microscale,const int mp_id,int &EXA_metric);

/** compute a microscale job on the master of a microscale work
    group. Collective communication on mpi_comm->micro. */
int micro_job_master(const PGFEM_mpi_comm *mpi_comm,
                     const int idx,
                     const int buff_size,
                     char *in_buffer,
                     char *out_buffer,
                     MICROSCALE *micro,
                     int *exit_server,
                     const int mp_id,
                     int &EXA_metric);

/** compute a microscale job on the slaves of a microscale work
    group. Collective communication on mpi_comm->micro */
int micro_job_slave(const PGFEM_mpi_comm *mpi_comm,
                    MICROSCALE *micro,const int mp_id,int &EXA_metric);

/** compute a microscale job. Collective communication on
    micro->common->mpi_comm. */
int microscale_compute_job(const int idx,      /**< sol/job id */
                           const int buff_len, /**< n_elem in buffer */
                           char *buffer,       /**< buffer of job info */
                           MICROSCALE *micro,
                           int *exit_server,
                           const int mp_id,
                           int &EXA_metric);


/** start computing struff for the macroscale --> initializes and
    communicates job information with waiting microscale
    server. Initializes non-blocing communication on ic->comm */
int start_macroscale_compute_jobs(const PGFEM_ms_job_intercomm *ic,
                                  const MACROSCALE *macro,
                                  const int job_type,
                                  const double *loc_sol,
                                  MS_COHE_JOB_INFO *job_list,
                                  PGFEM_server_ctx *send,
                                  PGFEM_server_ctx *recv);

/** finish the job macroscale job (typically involves assembly to
    tangent/residual). Finalizes non-blocking communication on
    ic->comm. Collective communication on mpi_comm->macro */
int finish_macroscale_compute_jobs(MS_COHE_JOB_INFO *job_list,
                                   MACROSCALE *macro,
                                   PGFEM_server_ctx *send,
                                   PGFEM_server_ctx *recv);

/** See function name. No communication */
int macroscale_update_job_info(const MACROSCALE *macro,
                               const int job_type,
                               const double *loc_sol,
                               MS_COHE_JOB_INFO *job);

/** updates the cohesive element associated with job */
int macroscale_update_coel(const MS_COHE_JOB_INFO *job,
                           MACROSCALE *macro);

#endif /* #ifndef  */
