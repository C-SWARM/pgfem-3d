#ifndef PGFEM3D_THREE_FIELD_ELEMENT_H
#define PGFEM3D_THREE_FIELD_ELEMENT_H


#include "PGFem3D_data_structure.h"
#include "crpl.h"
#include "elem3d.h"
#include "eps.h"
#include "femlib.h"
#include "new_potentials.h"
#include "sig.h"
#include "tensors.h"
#define N_VOL_TREE_FIELD 1

using pgfem3d::Solver;


/// compute element stiffness matrix in quasi steady state
///
/// \param[in]  fe    finite element helper object
/// \param[out] lk    computed element stiffness matrix
/// \param[in]  r_e   nodal variabls(displacements) on the current element
/// \param[in]  grid  a mesh object
/// \param[in]  mat   a material object
/// \param[in]  fv    object for field variables
/// \param[in]  alpha mid point alpha
/// \param[in]  dt    time step size
/// \return non-zero on internal error
int stiffmat_3f_el(FEMLIB *fe,
                   double *lk,
                   double *r_e,
                   Grid *grid,
                   MaterialProperty *mat,
                   FieldVariables *fv,
                   double alpha,
                   double dt);

/// compute element residual vector in quasi steady state.
///
/// \param[in]  fe    finite element helper object
/// \param[out] lk    computed element stiffness matrix
/// \param[in]  r_e   nodal variabls(displacements) on the current element
/// \param[in]  grid  a mesh object
/// \param[in]  mat   a material object
/// \param[in]  fv    object for field variables
/// \return non-zero on internal error
int residuals_3f_el(FEMLIB *fe,
                    double *f,
                    double *r_e,
                    Grid *grid,
                    MaterialProperty *mat,
                    FieldVariables *fv);

/// compute element residual vector in transient .
/// Total Lagrangian based three-field mixed method for Hyperelasticity.
///
/// \param[in]  fe    finite element helper object
/// \param[out] f     computed element stiffness matrix
/// \param[in]  r_e   nodal variabls(displacements) on the current element
/// \param[in]  grid  a mesh object
/// \param[in]  mat   a material object
/// \param[in]  fv    object for field variables
/// \param[in]  alpha mid point alpha
/// \param[in]  dts   time step size at t(n), t(n+1); dts[DT_N]   = t(n)   - t(n-1)
///                                                   dts[DT_NP1] = t(n+1) - t(n)
/// \return non-zero on internal error                    
int residuals_3f_w_inertia_el(FEMLIB *fe,
                              double *f,
                              double *r_e,
                              Grid *grid,
                              MaterialProperty *mat,
                              FieldVariables *fv,
                              Solver *sol,
                              const double *dts,
                              double *u_nm1,
                              double *u_npa);

/// Update variables during Newton Raphson iterations
///
/// \param[in] grid  a mesh object
/// \param[in] mat   a material object
/// \param[in] fv    object for field variables
/// \param[in] load  object for loading
/// \param[in] opts  structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] dts   time step size at t(n), t(n+1); dts[DT_N]   = t(n)   - t(n-1)
///                                                  dts[DT_NP1] = t(n+1) - t(n)
/// \param[in] alpha mid point alpha
/// \return non-zero on internal error                     
int update_3f_NR(Grid *grid,
                 MaterialProperty *mat,
                 FieldVariables *fv,
                 LoadingSteps *load,
                 const PGFem3D_opt *opts,
                 int mp_id,
                 const double *dts,
                 double alpha);

/// Update variables during Newton Raphson iterations
///
/// \param[in] grid  a mesh object
/// \param[in] mat   a material object
/// \param[in] fv    object for field variables
/// \param[in] load  object for loading
/// \param[in] mp_id mutiphysics id
/// \param[in] dt    time step size
/// \param[in] t     time
/// \param[in] mpi_comm MPI_COMM_WORLD
void update_3f_output_variables(Grid *grid,
                                MaterialProperty *mat,
                                FieldVariables *fv,
                                LoadingSteps *load,
                                const int mp_id,
                                const double dt,
                                const double t,
                                MPI_Comm mpi_comm);
                                
/// compute and set initial conditions for three field mixed method
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv array of field variable object
/// \return non-zero on internal error
void compute_3f_initial_conditions(Grid *grid,
                                   MaterialProperty *mat,
                                   FieldVariables *fv);
                                   
void compute_stress(double *GS, Element *elem, HOMMAT *hommat, long ne, int npres, Node *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm,int analysis);

#endif // #define PGFEM3D_THREE_FIELD_ELEMENT_H
