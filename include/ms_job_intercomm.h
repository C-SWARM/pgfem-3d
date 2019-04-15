/* HEADER */
#pragma once
#ifndef MS_JOB_INTERCOMM_H
#define MS_JOB_INTERCOMM_H

#include "pgfem3d/Communication.hpp"

/** Job intercommunicator structure */
typedef struct PGFEM_ms_job_intercomm{
  multiscale::net::MSNET_Comm comm;
  multiscale::MultiscaleCommInfo *send_info;
  multiscale::MultiscaleCommInfo *recv_info;
} PGFEM_ms_job_intercomm;

/** create a job intercommunicator structure */
int create_PGFEM_ms_job_intercomm(const int nproc_macro,
				  const multiscale::MultiscaleCommunicator *mscom,
				  const int n_jobs,
				  const int *job_buff_sizes,
				  PGFEM_ms_job_intercomm **intercomm);

/** destroy a job intercommunicator structure */
int destroy_PGFEM_ms_job_intercomm(PGFEM_ms_job_intercomm *intercomm);

#endif /* #ifndef  */
