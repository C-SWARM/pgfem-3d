#ifndef H__H__PMTL_LINEAR_SOLVER_INTERFACE__H__H
#define H__H__PMTL_LINEAR_SOLVER_INTERFACE__H__H

#include "mpi.h"
#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

int construct_solver(void **m);
int destruct_solver(void **m);

int initialize_linear_system(void *m, int N, int nne);
int update_linear_system_A_IJ(void *m, int *I, int IN,
                                       int *J, int JN, int N, double *values, MPI_Comm comm);
int set_linear_system_b(void *m, double *values, int N);
int set_solver_pc(void *m, int type);
int solve_linear_system(void *m, int N);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
