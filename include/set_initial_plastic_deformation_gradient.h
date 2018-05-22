#ifndef H__H__SET_INITIAL_PLASTIC_DEFORMATION_GRADIENT__H__H
#define H__H__SET_INITIAL_PLASTIC_DEFORMATION_GRADIENT__H__H

#include "PGFem3D_data_structure.h"

/// compute new geometry by given intial plastic deformation, and
/// set initial plastic deformation gradient in each element.
int set_initial_plastic_deformation_gradient(Grid *grid,
                                             FieldVariables *fv,
                                             MaterialProperty *mat,
                                             pgfem3d::Solver *sol,
                                             LoadingSteps *load,
                                             const pgfem3d::CommunicationStructure *com,
                                             const PGFem3D_opt *opts,
                                             const Multiphysics& mp,
                                             const int mp_id);
#endif
