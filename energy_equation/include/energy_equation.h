#pragma once
#ifndef H__PGFEM3D_ENERGY_EQUATION__H
#define H__PGFEM3D_ENERGY_EQUATION__H

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
int energy_equation_compute_residuals(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      LOADING_STEPS *load,
                                      const int mp_id,
                                      int use_updated,
                                      double dt);

/// compute stiffness for heat conduction problem
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in] fv field variable object 
/// \param[in] sol object for solution scheme
/// \param[in] com object for communications
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int energy_equation_compute_stiffness(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      SOLVER_OPTIONS *sol,
                                      COMMUNICATION_STRUCTURE *com,
                                      MPI_Comm mpi_comm,
                                      int myrank,
                                      const PGFem3D_opt *opts,
                                      const int mp_id,
                                      double dt);

/// compute flux due to Dirichlet BCs
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in] fv field variable object 
/// \param[in] sol object for solution scheme
/// \param[in] com object for communications
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error                                      
int energy_equation_compute_load4pBCs(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      SOLVER_OPTIONS *sol,
                                      LOADING_STEPS *load,
                                      int myrank,
                                      const PGFem3D_opt *opts,
                                      const int mp_id,
                                      double dt);
#endif
