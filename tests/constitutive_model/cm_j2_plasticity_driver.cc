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
#include "config.h"
#endif

#include "allocation.h"
#include "cm_j2_plasticity.h"
#include "constitutive_model.h"
#include "hommat.h"
#include "utils.h"
#include "index_macros.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>

enum {PARAM_G, 
      PARAM_nu, 
      PARAM_beta, 
      PARAM_hp, 
      PARAM_k0,
      PARAM_mu, 
      PARAM_w_max, 
      PARAM_P1, 
      PARAM_P2, 
      PARAM_Yin,
      PARAM_da,
      PARAM_db,
      PARAM_va,
      PARAM_vb,
      PARAM_NO};

static void hommat_assign_values(HOMMAT *p_hmat)
{
  // parameter values from Mosby and Matous, MSMSE (2015) 
  p_hmat->m01 = 0.0;
  p_hmat->E = 0.0;
  p_hmat->G = 6.0;
  p_hmat->m10 = 0.0;
  p_hmat->nu = 0.49;
  p_hmat->devPotFlag = -1;
  p_hmat->volPotFlag = 2;

  // unused 
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
  // elastic props 
  param[PARAM_G]  = p_hmat->G; // G - reset from HOMMAT 
  param[PARAM_nu] = p_hmat->nu; // nu - reset from HOMMAT 

  // plastic props 
  param[PARAM_beta] = 0.0; // beta 
  param[PARAM_hp]   = 3.0e2; // hp 
  param[PARAM_k0]   = 0.5; // k0 

  // damage props 
  param[PARAM_mu]    = 0.0; // mu 
  param[PARAM_w_max] = 1.0;
  param[PARAM_P1]    = 0.0; // p1 
  param[PARAM_P2]    = 0.0; // p2 
  param[PARAM_Yin]   = 0.0; // Yin  

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
                            CM_Ctx &cm_ctx,
                            const Constitutive_model *m,
                            const double t)
{
  int err = 0;
  double sig[9] = {};
  double J = det3x3(m->vars_list[0][m->model_id].Fs[1].m_pdata);
  err += m->param->update_elasticity(m, cm_ctx, NULL, sig, 0);  
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
  err += p.initialization(&mat, J2_PLASTICITY_DAMAGE);

  param_assign_values(p.model_param, &mat);    
  MATERIAL_CONTINUUM_DAMAGE mat_d;        
  set_damage_parameters(&mat_d, p.model_param[PARAM_P1],  p.model_param[PARAM_P2],
                                p.model_param[PARAM_Yin], p.model_param[PARAM_mu], 
                                p.model_param[PARAM_w_max]);

  MATERIAL_J2_PLASTICITY mat_J2p;
  set_J2_plasticity_parameters(&mat_J2p, 
                               p.model_param[PARAM_hp], 
                               p.model_param[PARAM_beta], 
                               p.model_param[PARAM_k0]);                                  

  p.cm_mat->mat_d = &mat_d;
  p.cm_mat->mat_J2p = &mat_J2p;
  
  err += m.initialization(&p);
  
  // get the model info
  Model_var_info info;
  err +=  m.param->get_var_info(info);
  err += info.print_variable_info(stdout);
    
  const double dt = 5.0;
  double t = dt;
  double F[9] = {0};
  const int n_step = 20;
  for (int i = 0; i < n_step; i++) {
    printf("STEP [%d]=================================\n",i);
    get_F(t,mat.nu, F);

    CM_Ctx cm_ctx;
    cm_ctx.set_tensors_ss(F);
    cm_ctx.set_time_steps_ss(dt);  
              
    err += m.param->integration_algorithm(&m, cm_ctx);
    err += m.param->update_state_vars(&m);
    //err += write_data_point(stdout,cm_ctx,&m,t);
    err += write_data_point(out,cm_ctx,&m,t);
    t += dt;
    /* printf("\n"); */
  }

  delete sv;
  fclose(out);
  return err;
}
