/**
 * This is a driver program for the BPA plasticity model.
 *
 * AUTHORS:
 *   Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "plasticity_model_BPA.h"
#include "_plasticity_model_BPA.h"
#include "allocation.h"
#include "constitutive_model.h"
#include "data_structure.h"
#include "hommat.h"
#include "index_macros.h"
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

using namespace gcm;

void set_mat_values(HOMMAT *m)
{
  /* parameter values from S. Holopanien, Mech. of Mat. (2013) */
  m->m01 = 0.0;
  m->E = 2.3e3;
  m->G = 0.8846e3;
  m->m10 = 0.5 * m->G;
  m->nu = 0.3;
  m->devPotFlag = 2;
  m->volPotFlag = 99;

  /* unused */
  m->M = NULL;
  m->L = NULL;
  m->density = 0.0;
  m->e1 = 0;
  m->e2 = 0;
  m->e3 = 0;
  m->e4 = 0;
}

void get_F(const double t,
           const double nu,
           double *F)
{
  const double rate = 0.001;
  memset(F, 0, 9*sizeof(*F));
  /* compression */
  F[0] = F[4] = 1 + nu * rate * t;
  F[8] = 1 - rate * t;
}

int compute_stress(double * __restrict sig,
                   const Constitutive_model *m,
                   const void *ctx)
{
  int err = 0;
  const HOMMAT *mat = m->param->p_hmat;
  const double kappa = (2* mat->G * (1 + mat->nu)) / (3 * (1 - 2 * mat->nu));
  double p = 0;
  m->param->compute_dudj(m,ctx,&p);
  p *= kappa;
  Matrix<double> S(3,3,0.0);
  m->param->compute_dev_stress(m,ctx,S.m_pdata);

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

int write_data_point(FILE *f,
                     const void *ctx,
                     const Constitutive_model *m,
                     const double t)
{
  int err = 0;
  double sig[9] = {};
  err += compute_stress(sig,m,ctx);
  fprintf(f,"%e\t%e\t%e\n",t,-sig[8],m->vars_list[0][m->model_id].state_vars->m_pdata[0]);
  return err;
}

int main(int argc, char **argv)
{
  int err = 0;
  FILE *in = fopen("bpa.props","r");
  
  int temp1, temp2;
  err += scan_for_valid_line(in);
  CHECK_SCANF(in, "%d", &temp1);
  
  err += scan_for_valid_line(in);
  CHECK_SCANF(in, "%d %d", &temp1, &temp2);
  err += scan_for_valid_line(in);
  
  int brace = fgetc(in);
  assert(brace == '{' && "Expect opening brace as next valid entry");
  
  FILE *out = fopen("bpa.dat","w");
  HOMMAT mat;
  BPA_PARAM p;
  Constitutive_model m;
  m.model_id = 0;
  State_variables *sv = new State_variables;
  m.vars_list = &sv;

  // initialize values
  set_mat_values(&mat);
  err += p.initialization(&mat, BPA_PLASTICITY);
  err += p.read_param(in);
  err += m.initialization(&p);

  // get the model info
  Model_var_info info;
  err +=  m.param->get_var_info(info);
  err += info.print_variable_info(stdout);

  double dt = 1.0;
  const int nstep = 100;
  double t = dt;
  double F[9] = {};
  void *ctx = NULL;
  int conv = 0;

  for (int i = 0; i < nstep; i++) {
    printf("STEP [%d]=================================\n",i);
    get_F(t,mat.nu,F);
    err += plasticity_model_BPA_ctx_build(&ctx,F,dt);
    conv = m.param->integration_algorithm(&m,ctx);
    assert(conv >= 0);
     if(conv < 0) { 
       printf("\t DID NOT CONVERGE, RESTARTING\n\n");
       dt /= 2;
       t -= dt;
       err += m.param->update_state_vars(&m);
       err += m.param->destroy_ctx(&ctx);
       i--;
       continue;
     } else err += conv;
    err += m.param->update_state_vars(&m);
    err += write_data_point(out,ctx,&m,t);
    err += m.param->destroy_ctx(&ctx);
    printf("\n");
    t += dt;
  }

  delete sv;
  fclose(out);
  fclose(in);
  return err;
}
