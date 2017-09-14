/* HEADER */
#pragma once
#ifndef PGFEM3D_SUBDIVISION_H
#define PGFEM3D_SUBDIVISION_H

#include "PGFEM_mpi.h"
#include "crpl.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "sig.h"
#include "supp.h"

typedef struct {
  int step_size;              /// step size
  int step_id;                /// stepping id with subdivision
  int decellerate;            /// decelleration parameter
  int accellerate;            /// accelleration parameter
  double dt_0;                /// saved time step size to accellerate or decellerate
  int reset_variables;        /// if 1, reset variables
  int is_subdivided;          /// 1 if subdivided
  int need_to_update_loading; /// if 1, new time step is updated. Need to update loads too,
                              /// in order to increase loads based on the new time step size
  double loading_factor;      /// loading increments factor
} SUBDIVISION_PARAM;

/// subdevide time step size
///
/// \param[in] was_NR_ok if 1, previous iteration was successful
/// \param[in, out] sp container of subdivision parameters
/// \param[in, out] dt time step size, new value will be updated
/// \param[in, out] times time at t(tim-1), t(tim), and t(tim+1), times[tim] will be updated
/// \param[in] tim time step id
/// \param[in] iter number of iteration taken in the iterative solver
/// \param[in] maximum number of iteration defined in the iterative solver
/// \param[in] alpha physics based evolution parameters
/// \return non-zero on internal error
int subdivision_scheme(int was_NR_ok,
                       SUBDIVISION_PARAM *sp,
                       double *dt,
                       double *times,
                       long tim,
                       int iter,
                       int max_iter,
                       double alpha,
                       MPI_Comm mpi_comm);

double subdiv_arc (long INFO,
                   double *dt,
                   double dt0,
                   long *STEP,
                   long *DIV,
                   long tim,
                   double *times,
                   long *ST,
                   long ne,
                   long ndofd,
                   long npres,
                   Element *elem,
                   CRPL *crpl,
                   EPS *eps,
                   SIG *sig,
                   SUPP sup,
                   double *sup_defl,
                   double *rr,
                   double *d_r,
                   double *D_R,
                   double *f_defl,
                   double *f,
                   long *GAMA,
                   double *DT,
                   long *OME,
                   double stab,
                   double dAL0,
                   double dAL,
                   double dALMAX,
                   double nor_min,
                   double dlm0,
                   long *ITT,
                   long iter,
                   long iter_max,
                   long TYPE,
                   MPI_Comm mpi_comm,
                   const int analysis);

#endif /* #define PGFEM3D_SUBDIVISION_H */
