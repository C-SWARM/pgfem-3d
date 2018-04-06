/**
 * This is a driver program for the J2 Plasticity + Damage model.
 *
 * REFERENCES:
 *
 * Simo, J. C., and J. W. Ju. "On continuum damage-elastoplasticity at
 * finite strains." Computational Mechanics 5.5 (1989): 375-400.
 *
 *
 * AUTHORS:
 *   Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "allocation.h"
#include "cm_j2_plasticity.h"
#include "constitutive_model.h"
#include "data_structure_c.h"
#include "hommat.h"
#include "index_macros.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>

static void hommat_assign_values(HOMMAT *p_hmat)
{
  /* parameter values from Mosby and Matous, MSMSE (2015) */
  p_hmat->m01 = 0.0;
  p_hmat->E = 0.0;
  p_hmat->G = 6.0;
  p_hmat->m10 = 0.0;
  p_hmat->nu = 0.49;
  p_hmat->devPotFlag = -1;
  p_hmat->volPotFlag = 2;

  /* unused */
  p_hmat->M = NULL;
  p_hmat->L = NULL;
  p_hmat->density = 0.0;
  p_hmat->e1 = 0;
  p_hmat->e2 = 0;
  p_hmat->e3 = 0;
  p_hmat->e4 = 0;
}

static void param_assign_values(double *param,
                                const HOMMAT *p_hmat)
{
  /* elastic props */
  param[0] = p_hmat->G; /* G - reset from HOMMAT */
  param[1] = p_hmat->nu; /* nu - reset from HOMMAT */

  /* plastic props */
  param[2] = 0.0; /* beta */
  param[3] = 3.0e2; /* hp */
  param[4] = 0.5; /* k0 */

  /* damage props */
  param[5] = 0.0; /* mu */
  param[6] = 0.0; /* p1 */
  param[7] = 0.0; /* p2 */
  param[8] = 0.0; /* Yin */

}

static int compute_stress(double * __restrict sig,
                          const Constitutive_model *m,
                          const void *ctx)
{
  int err = 0;
  double p = 0;
  const double kappa = hommat_get_kappa(m->param->p_hmat);
  err += m->param->compute_dudj(m,ctx,&p);
  p *= kappa;
  Matrix<double> S(3,3);
  err += m->param->compute_dev_stress(m,ctx,S.m_pdata);

  /* push stress forward */
  const double *Fe = m->vars_list[0][m->model_id].Fs[0].m_pdata;
  const double J = det3x3(Fe);
  memset(sig,0,9*sizeof(*sig));
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          sig[idx_2(i,j)] += Fe[idx_2(i,k)] * S.m_pdata[idx_2(k,l)] * Fe[idx_2(j,l)] / J;
        }
      }
    }
  }

  sig[0] += p;
  sig[4] += p;
  sig[8] += p;

  return err;
}

static void get_F(const double t,
                  const double nu,
                  double *F)
{
  const double rate = 0.0001;
  memset(F, 0, 9*sizeof(*F));
  /* F[0] = F[4] = 1 - nu * rate * t; */
  F[0] = F[4] = 1;
  F[8] = 1 + rate * t;
}

static int write_data_point(FILE *f,
                            const void *ctx,
                            const Constitutive_model *m,
                            const double t)
{
  int err = 0;
  double sig[9] = {};
  err += compute_stress(sig,m,ctx);
  double J = det3x3(m->vars_list[0][m->model_id].Fs[0].m_pdata);
  fprintf(f,"%e\t%e\t%e\n", t, sig[8], J);
  return err;
}

int main(int argc, char **argv)
{
  int err = 0;
  
  FILE *out = fopen("j2d.dat","w");
  HOMMAT mat; 
  CM_J2P_PARAM p;
  Constitutive_model m;
  m.model_id = 0;
  State_variables *sv = new State_variables;
  m.vars_list = &sv;

   // initialize values   
  hommat_assign_values(&mat);  
  err += p.initialization(&mat, BPA_PLASTICITY);
  param_assign_values(p.model_param, &mat);
  err += m.initialization(&p);
  
  // get the model info
  Model_var_info info;
  err +=  m.param->get_var_info(info);
  err += info.print_variable_info(stdout);
    
  const double dt = 5.0;
  double t = dt;
  double F[9] = {0};
  void *ctx = NULL;
  const int n_step = 20;
  for (int i = 0; i < n_step; i++) {
    printf("STEP [%d]=================================\n",i);
    get_F(t,mat.nu, F);
    err += j2d_plasticity_model_ctx_build(&ctx, F, dt);
    err += m.param->integration_algorithm(&m, ctx);
    err += m.param->update_state_vars(&m);
    //err += write_data_point(stdout,ctx,&m,t);
    err += write_data_point(out,ctx,&m,t);
    err += m.param->destroy_ctx(&ctx);
    t += dt;
    /* printf("\n"); */
  }

  delete sv;
  fclose(out);
  return err;
}
