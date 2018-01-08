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
/// \param[in]  alpha mid point alpha
/// \param[in]  dt    time step size
/// \return non-zero on internal error
int residuals_3f_el(FEMLIB *fe,
                    double *f,
                    double *r_e,
                    Grid *grid,
                    MaterialProperty *mat,
                    FieldVariables *fv);

void residuals_3f_w_inertia_el(double *f,const int ii,
                               const int ndofn,const int nne,const int npres,const int nVol,const int nsd,
                               const double *x,const double *y,const double *z,
                               const Element *elem,const HOMMAT *hommat,const Node *node,
                               const double *dts,SIG *sig,EPS *eps,double alpha, double *r_n_a, double *r_n_1_a);

int update_3f(Grid *grid,
              MaterialProperty *mat,
              FieldVariables *fv,
              LoadingSteps *load,
              const PGFem3D_opt *opts,
              Multiphysics *mp,
              int mp_id,
              const double dt,
              double alpha);              

void update_3f_state_variables(long ne, long ndofn, long npres, double *d_r, double *r,
                               Node *node, Element *elem, HOMMAT *hommat, SUPP sup, EPS *eps, SIG *sig, double dt, double t,
                               MPI_Comm mpi_comm, const int mp_id);

void compute_stress(double *GS, Element *elem, HOMMAT *hommat, long ne, int npres, Node *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm,int analysis);

#endif // #define PGFEM3D_THREE_FIELD_ELEMENT_H
