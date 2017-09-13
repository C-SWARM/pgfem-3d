#ifndef PGFEM3D_READ_INPUT_FILE_H
#define PGFEM3D_READ_INPUT_FILE_H

#include "Arc_length.h"
#include "PGFEM_mpi.h"
#include "PGFem3D_data_structure.h"
#include "PGFem3D_options.h"
#include "element.h"
#include "material.h"
#include "matgeom.h"
#include "mesh_load.h"
#include "node.h"
#include "supp.h"

/** Function for reading the entire input file. All required space
    is allocated within this function */
int read_input_file(const PGFem3D_opt *opts,
                    MPI_Comm comm,
                    long *nn,
                    long *Gnn,
                    long *ndofn,
                    long *ne,
                    long *lin_maxit,
                    double *lin_err,
                    double *lim_zero,
                    long *nmat,
                    long *n_concentrations,
                    long *n_orient,
                    Node **node,
                    Element **elem,
                    Material **material,
                    MATGEOM *matgeom,
                    SUPP *sup,
                    long *nln,
                    ZATNODE **znod,
                    long *nel_s,
                    ZATELEM **zelem_s,
                    long *nel_v,
                    ZATELEM **zelem_v,
                    const int *fv_ndofn,
                    const int phyicsno,
                    const int *ndim,
                    char **physicsnames);

/// Read mesh info, boundary conditions, and material properties.
/// from main input files (*.in)
///
/// \param[out] grid a mesh object
/// \param[out] mat a material object
/// \param[out] FV array of field variable object
/// \param[out] SOL array of solution scheme object
/// \param[out] load object for loading
/// \param[in] mp multiphysics object
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int read_mesh_file(Grid *grid,
                   MaterialProperty *mat,
                   FieldVariables *FV,
                   pgfem3d::Solver *SOL,
                   LoadingSteps *load,
                   Multiphysics *mp,
                   MPI_Comm mpi_comm,
                   const PGFem3D_opt *opts);

/// Read solver file for time stepping.
///
/// \param[out] time_steps object for time stepping
/// \param[out] mat a material object
/// \param[out] FV array of field variable object
/// \param[out] SOL array of solution scheme object
/// \param[out] load object for loading
/// \param[out] crpl object for lagcy crystal plasticity
/// \param[in] mp multiphysics object
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_solver_file(TimeStepping *ts,
                     MaterialProperty *mat,
                     FieldVariables *FV,
                     pgfem3d::Solver *SOL,
                     LoadingSteps *load,
                     CRPL *crpl,
                     Multiphysics *mp,
                     const PGFem3D_opt *opts,
                     int myrank);


/// Read initial conditions.
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[out] SOL array of solution scheme object
/// \param[out] load object for loading
/// \param[out] ts object for time stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[out] tnm1 if restart, read time step info from the previous run
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_initial_values(Grid *grid,
                        MaterialProperty *mat,
                        FieldVariables *FV,
                        pgfem3d::Solver *SOL,
                        LoadingSteps *load,
                        TimeStepping *ts,
                        PGFem3D_opt *opts,
                        Multiphysics *mp,
                        double *tnm1,
                        int myrank);

/// Read loads increments.
///
/// \param[in] grid a mesh object
/// \param[in] variables object for field variables
/// \param[out] load object for loading
/// \param[in] mp multiphysics object
/// \param[in] tim time step ID
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_and_apply_load_increments(Grid *grid,
                                   FieldVariables *variables,
                                   LoadingSteps *load,
                                   Multiphysics *mp,
                                   long tim,
                                   MPI_Comm mpi_comm,
                                   int myrank);

/// Read read cohesive elements.
///
/// \param[out] grid a mesh object
/// \param[out] mat a material object
/// \param[in] opts structure PGFem3D option
/// \param[in] ensight object
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_cohesive_elements(Grid *grid,
                           MaterialProperty *mat,
                           const PGFem3D_opt *opts,
                           Ensight *ensight,
                           MPI_Comm mpi_comm,
                           int myrank);

#endif /* #define PGFEM3D_READ_INPUT_FILE_H */
