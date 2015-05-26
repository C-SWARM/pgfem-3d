/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Adetokunbo Adedoyin, [1], <aadedoyi@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */

#include "state_variables.h"
#include "data_structure_c.h"

Define_Matrix(double);

int state_variables_build(State_variables *s)
{
  int err = 0;
  s->Fs = NULL;
  s->state_vars = NULL;
  s->n_Fs = 0;
  return err;
}

int state_variables_destroy(State_variables *s)
{
  int err = 0;
  for (size_t i = 0, e = s->n_Fs; i < e; i++){
    Matrix_cleanup((s->Fs[i]));
  }
  free(s->Fs);
  s->Fs = NULL;
  Matrix_cleanup((*(s->state_vars)));
  free(s->state_vars);
  s->state_vars = NULL;
  s->n_Fs = 0;
  return err;
}

int state_variables_initialize(State_variables *s,
                               const size_t n_Fs,
                               const size_t n_vars)
{
  int err = 0;
  s->n_Fs = n_Fs;
  s->Fs = malloc(n_Fs*sizeof(*(s->Fs)));
  for (size_t i = 0; i < n_Fs; i++) {
    Matrix_construct_redim(double,(s->Fs[i]),3,3);
  }
  /* state variables is column vector */
  s->state_vars = malloc(sizeof(*(s->state_vars)));
  Matrix_construct_redim(double,(*(s->state_vars)),3,1);
  return err;
}

