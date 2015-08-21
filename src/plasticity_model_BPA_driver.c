/**
 * This is a driver program for the BPA plasticity model.
 *
 * AUTHORS:
 *   Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "plasticity_model_BPA.h"
#include "_plasticity_model_BPA.h"
#include "constitutive_model.h"
#include "hommat.h"
#include "utils.h"
#include "data_structure_c.h"
#include "index_macros.h"

Define_Matrix(double);

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

int compute_stress(double * restrict sig,
                   const Constitutive_model *m,
                   const void *ctx)
{
  int err = 0;
  const HOMMAT *mat = m->param->p_hmat;
  const double kappa = (2* mat->G * (1 + mat->nu)) / (3 * (1 - 2 * mat->nu));
  double p = 0;
  m->param->compute_dudj(m,ctx,&p);
  p *= kappa / 2.0 ;
  Matrix_double S;
  Matrix_construct_redim(double,S,3,3);
  m->param->compute_dev_stress(m,ctx,&S);

  /* push stress forward */
  const double *Fe = m->vars.Fs[0].m_pdata;
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

void write_data_point(FILE *f,
                      const double *sig,
                      const double t)
{
  fprintf(f,"%e\t%e\n",t,-sig[8]);
}

int main(int argc, char **argv)
{
  int err = 0;
  FILE *out = fopen("bpa.dat","w");
  HOMMAT *mat = malloc(sizeof(*mat));
  Model_parameters *p = malloc(sizeof(*p));
  Constitutive_model *m = malloc(sizeof(*m));

  /* initialize values */
  set_mat_values(mat);
  err += model_parameters_construct(p);
  err += model_parameters_initialize(p,NULL,NULL,mat,BPA_PLASTICITY);
  err += constitutive_model_construct(m);
  err += constitutive_model_initialize(m,p);
  err += plasticity_model_BPA_set_initial_values(m);

  /* get the model info */
  Model_var_info *info = NULL;
  err +=  m->param->get_var_info(&info);
  err += model_var_info_print(stdout,info);
  err += model_var_info_destroy(&info);

  const double dt = 1.0;
  const int nstep = 100;
  double t = dt;
  double F[9] = {};
  double sig[9] = {};
  void *ctx = NULL;

  write_data_point(out,sig,0);
  for (int i = 0; i < nstep; i++) {
    printf("STEP [%d]=================================\n",i);
    get_F(t,mat->nu,F);
    err += plasticity_model_BPA_ctx_build(&ctx,F,dt);
    err += m->param->integration_algorithm(m,ctx);
    err += m->param->update_state_vars(m);
    err += compute_stress(sig,m,ctx);
    err += plasticity_model_BPA_ctx_destroy(&ctx);
    printf("\n");
    write_data_point(out,sig,t);
    t += dt;
  }

  /* call destructors */
  err += model_parameters_destroy(p);
  err += constitutive_model_destroy(m);

  /* free pointers */
  free(mat);
  free(p);
  free(m);
  fclose(out);
  return err;
}
