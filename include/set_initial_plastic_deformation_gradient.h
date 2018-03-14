#ifndef H__H__SET_INITIAL_PLASTIC_DEFORMATION_GRADIENT__H__H
#define H__H__SET_INITIAL_PLASTIC_DEFORMATION_GRADIENT__H__H

#include "PGFem3D_data_structure.h"
#include "pgfem3d/Solver.hpp"

/// compute new geometry by given intial plastic deformation, and
/// set initial plastic deformation gradient in each element.
int set_initial_plastic_deformation_gradient(Grid *grid,
                                             FieldVariables *fv,
                                             MaterialProperty *mat,                                             
                                             pgfem3d::Solver *sol,
                                             LoadingSteps *load,
                                             const CommunicationStructure *com,
                                             const MPI_Comm mpi_comm,
                                             const PGFem3D_opt *opts,
                                             const Multiphysics& mp,
                                             const int mp_id,
                                             const int myrank);
#endif
