/* HEADER */
#pragma once
#ifndef MS_JOB_INTERCOMM_H
#define MS_JOB_INTERCOMM_H

#include "PGFEM_mpi.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Job intercommunicator structure */
  typedef struct PGFEM_ms_job_intercomm{
    MPI_Comm comm;
    PGFEM_comm_info *send_info;
    PGFEM_comm_info *recv_info;
  } PGFEM_ms_job_intercomm;

  /** create a job intercommunicator structure */
  int create_PGFEM_ms_job_intercomm(const int nproc_macro,
				    const PGFEM_mpi_comm *comm,
				    const int n_jobs,
				    const int *job_buff_sizes,
				    PGFEM_ms_job_intercomm **intercomm);

  /** destroy a job intercommunicator structure */
  int destroy_PGFEM_ms_job_intercomm(PGFEM_ms_job_intercomm *intercomm);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
