/**
 * Function definitions of some common placeholder functions to
 * improve code reuse.
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#include "cm_placeholder_functions.h"
#include "utils.h"
#include <ttl/ttl.h>

static const int dim = 3;
static const int tensor = 9;

//ttl declarations
namespace {

  template<int R, int D = 3, class S = double>
  using Tensor = ttl::Tensor<R, D, S>;
      
  template<int R, int D = 3, class S = double *>
  using TensorA = ttl::Tensor<R, D, S>;
    
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
}

int cm_get_var_zero(const Constitutive_model *m,
                    double *var) { *var = 0.0; return 0; }

int cm_get_F_zero(const Constitutive_model *m,
                  double *F) 
{
  memset(F, 0, tensor*sizeof(*F));
  return 0; 
}

int cm_get_F_eye(const Constitutive_model *m,
                 double *F_in)
{
  TensorA<2> F(F_in);
  F = ttl::identity(i,j);
  return 0; 
}

int cm_get_lam_p(const Constitutive_model *m,
                 double *lam_p)
{
  int err = 0;
  Tensor<2> Fp, Cp = {};

  /* get the plastic deformation gradient at (n + 1) */
  m->param->get_pF(m, Fp.data,1);

  /* compute lam_p = sqrt( tr(Cp) / 3 ) */
  /* tr(Cp) = Fp : Fp */
  Cp = Fp(i,k)*Fp(j,k);
  *lam_p = Cp(i,j)*Cp(i,j);
  
  *lam_p = sqrt( *lam_p / 3.0 );
  return err;
}

int cm_compute_null_dMdu(const Constitutive_model *m,
                         CM_Ctx &cm_ctx,
                         const double *Grad_op,
                         const int nne,
                         const int ndofn,
                         double *dM_du)
{
  memset(dM_du, 0, tensor * nne *dim * sizeof(*dM_du));
  return 0;
}

int cm_no_subdiv(const Constitutive_model *m,
                 double *subdiv_param,
                 const double dt) { *subdiv_param = 0.0; return 0; }
                 
int cm_get_eq_plastic_strain(const Constitutive_model *m,
                             double *eq)
{
  return 0.0;
}