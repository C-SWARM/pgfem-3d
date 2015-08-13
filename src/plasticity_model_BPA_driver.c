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

void set_mat_values(HOMMAT *m)
{
  m->m01 = 0.0;
  m->E = 15.0e3;
  m->G = 6.0e3;
  m->m10 = 0.5 * m->G;
  m->nu = 0.25;
  m->devPotFlag = 2;
  m->devPotFlag = 99;

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
  memset(F, 0, 9*sizeof(*F));
  /* compression */
  F[0] = F[4] = 1 + nu * t;
  F[8] = 1 + t;
}

int main(int argc, char **argv)
{
  int err = 0;
  HOMMAT *mat = malloc(sizeof(*mat));
  Model_parameters *p = malloc(sizeof(*p));
  Constitutive_model *m = malloc(sizeof(*m));

  /* initialize values */
  set_mat_values(mat);
  err += model_parameters_construct(p);
  err += model_parameters_initialize(p,NULL,NULL,mat,BPA_PLASTICITY);
  err += constitutive_model_construct(m);
  err += constitutive_model_initialize(m,p);

  /* get the model info */
  Model_var_info *info = NULL;
  err +=  m->param->get_var_info(&info);
  err += model_var_info_print(stdout,info);
  err += model_var_info_destroy(&info);

  const double dt = 0.001;
  const int nstep = 0;
  double t = dt;
  double F[9];
  for (int i = 0; i < nstep; i++) {
    get_F(t,mat->nu,F);

    t += dt;
  }

  /* call destructors */
  err += model_parameters_destroy(p);
  err += constitutive_model_destroy(m);

  /* free pointers */
  free(mat);
  free(p);
  free(m);
  return err;
}
