/* HEADER */
#pragma once
#ifndef MS_COHE_JOB_INFO_H
#define MS_COHE_JOB_INFO_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#include <stdlib.h> /* for size_t */
#include "PGFEM_io.h"

  /* Job types */
  enum{JOB_NO_COMPUTE_EQUILIBRIUM, /**< do not compute micro equilibrium */
       JOB_COMPUTE_EQUILIBRIUM, /**< compute micro equilibrium */
       JOB_UPDATE, /**< update microscale state and return.*/
       JOB_PRINT, /**< Print microscale if print_flag == true */
       JOB_EXIT}; /**< no more jobs will be computed, close server */

  /** structure containing information required for multiscale
      coheisve modeling at a single integration point. All information
      is from the macroscale */
  typedef struct MS_COHE_JOB_INFO{
    int nnode; /**< # nodes on associated macro element */
    int ndofe; /**< # dofs on associated macro element (nnode*ndim) */
    int elem_id; /**< macro element id */
    int proc_id; /**< processor ID (in macro communicator) that owns
		    the macro element */
    int int_pt;
    int job_type; /**< flag to guide microscale computation */
    int print_flag; /**< flag to print during update */

    double int_wt; /**< macro integration weight including the
		      transformation */

    double *jump; /**< macro jump across interface [ndim]*/
    double *jump_n; /**< macro jump across interface at time (n)
		       [ndim]. Modified by microscale only */
    double *normal; /**< macro normal to interface [ndim]*/
    double *traction; /**< current traction computed at the microscale
			 [ndim]. Modified by microscale only */

    double *traction_n; /**< traction computed at time n. Modified by
			   microscale only. */

    double max_traction; /**< state variable. Maintained by microscale. */
    double max_jump;  /**< state variable. Maintained by microscale. */

    double *shape; /**< shape function for each macro node [nnode]*/

    double *traction_res; /**< macro traction residual [ndofe] */

    double *K_00_contrib; /* element matrix [ndofe*ndofe]*/

    long *loc_dof_ids; /**< macro local dof ids [ndofe] */  
    long *g_dof_ids; /**< macro global dof ids [ndofe] */

    /* solution procedure */
    int tim; /**< time step at macroscale */
    double *times; /**< macro time at (tim-1) (tim) and (tim+1) [3]*/
    int n_step;

  }MS_COHE_JOB_INFO;

  size_t compute_MS_COHE_JOB_INFO_size(const MS_COHE_JOB_INFO *info);   

  int build_MS_COHE_JOB_INFO(MS_COHE_JOB_INFO *info,
			     const int nnode);

  int build_MS_COHE_JOB_INFO_buffer(const size_t buff_len,
				    const char *buffer,
				    MS_COHE_JOB_INFO *info);

  int set_MS_COHE_JOB_INFO(MS_COHE_JOB_INFO *info,
			   const double *normal,
			   const double *jump,
			   const double *shape,
			   const long *loc_dof_ids,
			   const long *g_dof_ids);

  void destroy_MS_COHE_JOB_INFO(MS_COHE_JOB_INFO *info);

  /** pack the job info into a contiguous buffer */
  int pack_MS_COHE_JOB_INFO(const MS_COHE_JOB_INFO *info,
			    const size_t buffer_len,
			    char *buffer);

  /** unpack the buffer into a job. NOTE: the job must be the same
      size (nnodes) */
  int unpack_MS_COHE_JOB_INFO(MS_COHE_JOB_INFO *info,
			      const size_t buffer_len,
			      const char *buffer);

  /** print the job information to an output stream  */
  int print_MS_COHE_JOB_INFO(FILE *out,
			     const MS_COHE_JOB_INFO *info);
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
