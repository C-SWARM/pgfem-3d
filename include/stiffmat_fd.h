/* HEADER */
#pragma once
#ifndef STIFFMAT_FD_H
#define STIFFMAT_FD_H

#include "data_structure.h"
#include "bounding_element.h"
#include "cohesive_element.h"
#include "crpl.h"
#include "element.h"
#include "eps.h"
#include "femlib.h"
#include "hommat.h"
#include "macro_micro_functions.h"
#include "matgeom.h"
#include "node.h"
#include "pgfem_comm.h"
#include "PGFEM_mpi.h"
#include "PGFem3D_data_structure.h"
#include "PGFem3D_options.h"
#include "sig.h"
#include "supp.h"
#include "pgfem3d/Solver.hpp"

/// compute element stiffness matrix
///
/// \param[in] fe finite element object
/// \param[out] lk computed local(element level) stiffness matrix
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step
/// \param[in] lm Load multiplier level in Arc Length scheme
/// \param[in] be tangential load vector when periodic and solution scheme is Arc Length
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \return non-zero on internal error
int el_compute_stiffmat_MP(FEMLIB *fe,
                           double *lk,
                           GRID *grid,
                           MATERIAL_PROPERTY *mat,
                           FIELD_VARIABLES *fv,
                           SOLVER_OPTIONS *sol,
                           LOADING_STEPS *load,
                           CRPL *crpl,
                           const PGFem3D_opt *opts,
                           MULTIPHYSICS *mp,
                           int mp_id,
                           double dt,
                           double lm,
                           double *be,
                           double *r_e);

/// Compute stiffnes
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com communication object
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dt time step
/// \param[in] iter number of Newton Raphson interataions
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int stiffmat_fd_MP(GRID *grid,
                   MATERIAL_PROPERTY *mat,
                   FIELD_VARIABLES *fv,
                   SOLVER_OPTIONS *sol,
                   LOADING_STEPS *load,
                   COMMUNICATION_STRUCTURE *com,
                   CRPL *crpl,
                   MPI_Comm mpi_comm,
                   const PGFem3D_opt *opts,
                   MULTIPHYSICS *mp,
                   int mp_id,
                   double dt,
                   long iter,
                   int myrank);

/// Multiscale simulation interface to compute stiffness matrix
///
/// \param[in] c structure of macroscale information
/// \param[in,out] s contains the information for the history-dependent solution
/// \param[in] opts structure PGFem3D option
/// \param[in] iter number of Newton Raphson interataions
/// \param[in] nor_min nonlinear convergence tolerance
/// \param[in] FNR if 1: Full Newton-Raphson
///                   0: only compute stiffnes at the 1st iteration
/// \param[in] myrank current process rank
/// \param[in] nproc   number of total process
/// \return non-zero on internal error
int stiffmat_fd_multiscale(COMMON_MACROSCALE *c,
                           MACROSCALE_SOLUTION *s,
                           const PGFem3D_opt *opts,
                           long iter,
                           double nor_min,
                           long FNR,
                           int myrank,
                           int nproc);

/** Assemble non-local parts as they arrive */
int assemble_nonlocal_stiffmat(const COMMUN pgfem_comm,
                               MPI_Status *sta_r,
                               MPI_Request *req_r,
                               pgfem3d::solvers::SparseSystem *system,
                               double **recv);

#endif /* #ifndef STIFFMAT_FD_H */
