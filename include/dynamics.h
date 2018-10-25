#ifndef PGFEM3D_DYNAMICS_H
#define PGFEM3D_DYNAMICS_H

#include "PGFem3D_data_structure.h"
#include "PGFem3D_options.h"
#include "crpl.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "node.h"
#include "sig.h"
#include "supp.h"
#include "hyperelasticity.h"

const constexpr double MIN_DENSITY = 1.0e-16;
const constexpr int         DT_NP1 = 0;
const constexpr int           DT_N = 1;

void DISP_resid_body_force_el(double *f,
                              const int ii,
                              const int ndofn,
                              const int nne,
                              const double *x,
                              const double *y,
                              const double *z,
                              const Element *elem,
                              const HOMMAT *hommat,
                              const Node *node, 
                              double dt, 
                              double t,
                              HyperElasticity *elast,
                              bool is4cm);

struct FEMLIB;

/// compute element residual vector in transient
///
/// \param[in] fe finite element helper object
/// \param[out] be computed element residual vector
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dts time step size at t(n), t(n+1); dts[DT_N] = t(n) - t(n-1)
///                                                dts[DT_NP1] = t(n+1) - t(n)
/// \param[in] t current time
/// \return non-zero on internal error
int residual_with_inertia(FEMLIB *fe,
                          double *be,
                          double *r_e,
                          Grid *grid,
                          MaterialProperty *mat,
                          FieldVariables *fv,
                          pgfem3d::Solver *sol,
                          LoadingSteps *load,
                          CRPL *crpl,
                          const PGFem3D_opt *opts,
                          const Multiphysics& mp,
                          int mp_id,
                          double *dts,
                          double t,
			  int myrank);

/// compute element stiffness matrix in transient
///
/// \param[in] fe finite element helper object
/// \param[out] Ke computed computed element stiffness matrix
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int stiffness_with_inertia(FEMLIB *fe,
                           double *Ks,
                           double *r_e,
                           Grid *grid,
                           MaterialProperty *mat,
                           FieldVariables *fv,
                           pgfem3d::Solver *sol,
                           LoadingSteps *load,
                           CRPL *crpl,
                           const PGFem3D_opt *opts,
                           const Multiphysics& mp,
                           int mp_id,
                           double dt);

#endif /* #define PGFEM3D_DYNAMICS_H */
