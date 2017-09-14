/**
 * AUTHORS:
 * Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 * Karel Matous, University of Notre Dame, kmatous [at] nd.
 */
#ifndef PGFEM3D_FD_RESIDUALS_H
#define PGFEM3D_FD_RESIDUALS_H

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
#include "macro_micro_functions.h"
#include "matgeom.h"
#include "node.h"
#include "sig.h"
#include "solver_file.h"
#include "supp.h"

/// Compute residuals
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
long fd_residuals_MP(Grid *grid,
                     MaterialProperty *mat,
                     FieldVariables *fv,
                     pgfem3d::Solver *sol,
                     LoadingSteps *load,
                     CRPL *crpl,
                     MPI_Comm mpi_comm,
                     const PGFem3D_opt *opts,
                     Multiphysics *mp,
                     int mp_id,
                     double t,
                     double *dts,
                     int updated_deformation);

/// compute the reaction force for multiphysics mode
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
int fd_res_compute_reactions_MP(Grid *grid,
                                MaterialProperty *mat,
                                FieldVariables *fv,
                                pgfem3d::Solver *sol,
                                LoadingSteps *load,
                                CRPL *crpl,
                                MPI_Comm mpi_comm,
                                const PGFem3D_opt *opts,
                                Multiphysics *mp,
                                int mp_id,
                                double t,
                                double *dts);

/// Multiscale simulation interface computing reaction forces
///
/// \param[in] c structure of macroscale information
/// \param[in] s contains the information for the history-dependent solution
/// \param[in] solver_file structure for storing/updating the data
/// \param[in] opts structure PGFem3D option
/// \param[in] dts time step sizes at n, and n+1
/// \return non-zero on internal error
int fd_res_compute_reactions_multiscale(COMMON_MACROSCALE *c,
                                        MACROSCALE_SOLUTION *s,
                                        SOLVER_FILE *solver_file,
                                        const PGFem3D_opt *opts,
                                        double *dts);

#endif /* #define PGFEM3D_FD_RESIDUALS_H */
