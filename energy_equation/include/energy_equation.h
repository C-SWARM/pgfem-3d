/// Declare energy equation: function for computing stiffness matrix and residual vector
///
/// Authors:
///  Sangmin Lee, [1], <slee43@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN
#ifndef PGFEM3D_ENERGY_EQUATION_H
#define PGFEM3D_ENERGY_EQUATION_H

#include "PGFem3D_data_structure.h"

/// compute residuals for heat conduction problem
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv field variable object
/// \param[in] load object for loading
/// \param[in] mp_id mutiphysics id
/// \param[in] use_updated if use_updated=1, compute residuals updated temperature
///                           use_updated=0, compute residuals using temporal temperature
/// \param[in] dt time step size
/// \return non-zero on internal error
int energy_equation_compute_residuals(Grid *grid,
                                      MaterialProperty *mat,
                                      FieldVariables *fv,
                                      LoadingSteps *load,
                                      const int mp_id,
                                      int use_updated,
                                      double dt);

/// compute stiffness for heat conduction problem
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in] fv field variable object
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com object for communications
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int energy_equation_compute_stiffness(Grid *grid,
                                      MaterialProperty *mat,
                                      FieldVariables *fv,
                                      pgfem3d::Solver *sol,
                                      LoadingSteps *load,
                                      CommunicationStructure *com,
                                      MPI_Comm mpi_comm,
                                      int myrank,
                                      const PGFem3D_opt *opts,
                                      const int mp_id,
                                      double dt);

/// compute flux due to Dirichlet BCs
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv field variable object
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int energy_equation_compute_load4pBCs(Grid *grid,
                                      MaterialProperty *mat,
                                      FieldVariables *fv,
                                      pgfem3d::Solver *sol,
                                      LoadingSteps *load,
                                      int myrank,
                                      const PGFem3D_opt *opts,
                                      const int mp_id,
                                      double dt);

/// update for for print
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv field variable object
/// \param[in] dt time step size
/// \return non-zero on internal error
int update_thermal_flux4print(Grid *grid,
                              MaterialProperty *mat,
                              FieldVariables *fv,
                              double dt);

#endif // #define PGFEM3D_ENERGY_EQUATION_H
