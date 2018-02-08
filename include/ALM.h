#ifndef PGFEM3D_ALM_H
#define PGFEM3D_ALM_H

#include "PGFEM_mpi.h"
#include "PGFem3D_data_structure.h"
#include "PGFem3D_options.h"
#include "bounding_element.h"
#include "cohesive_element.h"
#include "crpl.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "matgeom.h"
#include "node.h"
#include "pgfem_comm.h"
#include "sig.h"
#include "supp.h"

/** */
double D_lam_ALM (long ndofd,
                  double *BS_rr,
                  double *BS_d_r,
                  double *BS_D_R,
                  double *BS_R,
                  double *BS_DK,
                  double dlm,
                  double dAL,
                  long *DomDof,
                  MPI_Comm mpi_comm);

/** Returns 0. */
double d_ALM2 (long ndofd,
               double *rr,
               double *R,
               double *DK,
               double d_lm);

/** Returns 0. */
double d_lam_ALM2 (long ndofd,
                   double *rr,
                   double *R,
                   double *DK,
                   double dAL,
                   double DET,
                   double DET0,
                   double dlm0,
                   double nor_min,
                   double *dR);

/// D_lam_ALM2_MP
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv array of field variable object
/// \param[in] sol object array for solution scheme
/// \param[in] load object for loading
/// \param[in] com object array for communications
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dlm Arc_length parameter
/// \param[in] dAL Arc_length parameter
/// \param[in] dt times step size
/// \return computed DLM
double D_lam_ALM2_MP(Grid *grid,
                     MaterialProperty *mat,
                     FieldVariables *fv,
                     pgfem3d::Solver *sol,
                     LoadingSteps *load,
                     CommunicationStructure *com,
                     CRPL *crpl,
                     MPI_Comm mpi_comm,
                     const PGFem3D_opt *opts,
                     const Multiphysics& mp,
                     int mp_id,
                     double dlm,
                     double dAL,
                     double dt);

double d_ALM4 (long ndofd,
               double *BS_rr,
               double *BS_DK,
               double dlm,
               long *DomDof,
               MPI_Comm mpi_comm);

double d_lam_ALM4 (long ndofd,
                   double *BS_rr,
                   double *BS_DK,
                   double *BS_dR,
                   double dAL,
                   long *DomDof,
                   MPI_Comm mpi_comm);

double D_lam_ALM4 (long ndofd,
                   double *BS_rr,
                   double *BS_d_r,
                   double *BS_D_R,
                   double *BS_DK,
                   double dlm,
                   double dAL,
                   long *DomDof,
                   MPI_Comm mpi_comm);

#endif /* #define PGFEM3D_ALM_H */
