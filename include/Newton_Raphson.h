/* HEADER */
/**
 *  Newton-Raphson solution algorithm.
 *
 * This is the Newton-Raphson nonlinear solution algorithm. It
 * includes a globally convergent line-search algorithm and
 * subdivision procedure. If at any point the linear solution does
 * not converge to the specified tolerance (+ an additional
 * tolerance for "wiggle room"), the load is subdivided and we
 * attempt to solve the prescribed deformation in multiple
 * increments. The prescribed load is likewise subdivided if the
 * solution is not obtained within the specified number of nonlinear
 * iterations or if the line-search algorithm fails to converge. The
 * relative residual is checked for convergence. The error norm is
 * also checked for convergence to ensure that the solution
 * progresses in the case that the relative residual is not less
 * than the tolerance but the total residual is very small (err^2).
 *
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 *  Sangmin Lee, University of Notre Dame, <slee43 [at] nd.edu>
 *  Karel Matous, University of Notre Dame, <kmatous [at] nd.edu>
 */
#pragma once
#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

#include "PGFem3D_data_structure.h"
#include "PGFEM_mpi.h"
#include "element.h"
#include "matgeom.h"
#include "hommat.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "cohesive_element.h"
#include "bounding_element.h"
#include "PGFem3D_options.h"
#include "pgfem_comm.h"

#include "hypre_global.h"
#include "macro_micro_functions.h"
#include "solver_file.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

/// Compute residuals for Newton Raphson iteration
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dts time step sizes a n, and n+1
/// \return non-zero on internal error
long compute_residuals_for_NR(GRID *grid,
                              MATERIAL_PROPERTY *mat,
                              FIELD_VARIABLES *fv,
                              SOLVER_OPTIONS *sol,
                              LOADING_STEPS *load,
                              CRPL *crpl,
                              MPI_Comm mpi_comm,
                              const PGFem3D_opt *opts,
                              MULTIPHYSICS *mp,
                              int mp_id,
                              double t,
                              double *dts,
                              int updated_deformation);

/// Perform Newton Staggered Newton Raphson
///
/// \param[in] print_level print level for a summary of the entire function call
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] VVolume original volume of the domain
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \return time spent for this routine
double Multiphysics_Newton_Raphson(const int print_level,
                                   GRID *grid,
                                   MATERIAL_PROPERTY *mat,
                                   FIELD_VARIABLES *FV,
                                   SOLVER_OPTIONS *SOL,
                                   LOADING_STEPS *load,
                                   COMMUNICATION_STRUCTURE *COM,
                                   PGFem3D_TIME_STEPPING *time_steps,
                                   CRPL *crpl,
                                   MPI_Comm mpi_comm,
                                   const double VVolume,
                                   const PGFem3D_opt *opts,
                                   MULTIPHYSICS *mp);

/// Multiscale simulation interface to perform Newton Raphson iteration
///
/// \param[in] print_level print level for a summary of the entire function call
/// \param[in] c structure of macroscale information
/// \param[in,out] s contains the information for the history-dependent solution
/// \param[in] solver_file structure for storing/updating the data
/// \param[in] ctx container for passing through Newton Raphson
/// \param[in] opts structure PGFem3D option
/// \param[in] sup_defl Prescribed deflection
/// \param[out] pores opening volume of failed cohesive interfaces
/// \param[out] the number of nonlinear steps taken to solve the given increment
/// \return time spent in linear solver (seconds).
double Newton_Raphson_multiscale(const int print_level,
                                 COMMON_MACROSCALE *c,
                                 MACROSCALE_SOLUTION *s,
                                 SOLVER_FILE *solver_file,
                                 MS_SERVER_CTX *ctx,
                                 const PGFem3D_opt *opts,
                                 double *sup_defl,
                                 double *pores,
                                 int *n_step);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef NEWTON_RAPHSON_H */
