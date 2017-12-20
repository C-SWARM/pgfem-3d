/* HEADER */
#pragma once
#ifndef MS_JOB_INTERCOMM_H
#define MS_JOB_INTERCOMM_H

#include "pgfem3d/MultiscaleCommunication.hpp"

/** Job intercommunicator structure */
typedef struct PGFEM_ms_job_intercomm{
  pgfem3d::net::PGFem3D_Comm comm;
  pgfem3d::MultiscaleComm *send_info;
  pgfem3d::MultiscaleComm *recv_info;
} PGFEM_ms_job_intercomm;

/** create a job intercommunicator structure */
int create_PGFEM_ms_job_intercomm(const int nproc_macro,
				  const pgfem3d::MultiscaleComm *mscom,
				  const int n_jobs,
				  const int *job_buff_sizes,
				  PGFEM_ms_job_intercomm **intercomm);

/** destroy a job intercommunicator structure */
int destroy_PGFEM_ms_job_intercomm(PGFEM_ms_job_intercomm *intercomm);

#endif /* #ifndef  */
