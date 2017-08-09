/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */

#include "state_variables.h"
#include "data_structure_c.h"
#include "utils.h"

int state_variables_build(State_variables *s)
{
  int err = 0;
  s->Fs = NULL;
  s->state_vars = NULL;
  s->flags = NULL;
  s->n_Fs = 0;
  s->n_flags = 0;
  return err;
}

int state_variables_destroy(State_variables *s)
{
  int err = 0;
  for (size_t i = 0, e = s->n_Fs; i < e; i++){
    Matrix_cleanup(s->Fs[i]);
  }
  free(s->Fs);
  s->Fs = NULL;
  s->n_Fs = 0;

  if (s->state_vars) {
    Matrix_cleanup(s->state_vars[0]);
  }
  free(s->state_vars);
  s->state_vars = NULL;

  free(s->flags);
  s->flags = NULL;
  s->n_flags = 0;

  return err;
}

int state_variables_initialize(State_variables *s,
                               const size_t n_Fs,
                               const size_t n_vars,
                               const size_t n_flags)
{
  int err = 0;
  s->n_Fs = n_Fs;
  s->Fs = malloc(n_Fs*sizeof(*(s->Fs)));
  for(size_t i = 0; i < n_Fs; i++)
  {
    Matrix_construct(double,s->Fs[i]);
    Matrix_eye(s->Fs[i],3);
  }

  /* state variables is column vector */
  s->state_vars = malloc(sizeof(*(s->state_vars)));

  Matrix_construct_init(double, s->state_vars[0], n_vars, 1, 0);

  s->n_flags = n_flags;
  s->flags = calloc(n_flags, sizeof(*(s->flags)));

  return err;
}

size_t state_variables_get_packed_size(const State_variables *s)
{
  /* get the length of all Fs. Note that we are computing the length
     of each tensor just in case. */
  size_t len_Fs = 0;
  for (int i = 0, e = s->n_Fs; i < e; i++) {
    len_Fs += s->Fs[i].m_col * s->Fs[i].m_row;
  }
  const size_t len_vars = s->state_vars->m_col * s->state_vars->m_row;
  return ( (len_Fs + len_vars) * sizeof(double) + s->n_flags * sizeof(int) );
}

int state_variables_pack(const State_variables *s,
                         char *buffer,
                         size_t *pos)
{
  int err = 0;

  /* pack the deformation gradients */
  for (size_t i = 0, e = s->n_Fs; i < e; i++) {
    const size_t len = s->Fs[i].m_row * s->Fs[i].m_col;
    const double *data = s->Fs[i].m_pdata;
    pack_data(data, buffer, pos, len, sizeof(*data));
  }

  /* pack the state variables */
  {
    const size_t len = s->state_vars->m_row * s->state_vars->m_col;
    const double *data = s->state_vars->m_pdata;
    pack_data(data, buffer, pos, len, sizeof(*data));
  }

  /* pack the flags */
  pack_data(s->flags, buffer, pos, s->n_flags, sizeof(*(s->flags)));

  return err;
}

int state_variables_unpack(State_variables *s,
                           const char *buffer,
                           size_t *pos)
{
  int err = 0;

  /* pack the deformation gradients */
  for (size_t i = 0, e = s->n_Fs; i < e; i++) {
    const size_t len = s->Fs[i].m_row * s->Fs[i].m_col;
    double *data = s->Fs[i].m_pdata;
    unpack_data(buffer, data, pos, len, sizeof(*data));
  }

  /* pack the state variables */
  {
    const size_t len = s->state_vars->m_row * s->state_vars->m_col;
    double *data = s->state_vars->m_pdata;
    unpack_data(buffer, data, pos, len, sizeof(*data));
  }

  /* pack the flags */
  unpack_data(buffer, s->flags, pos, s->n_flags, sizeof(*(s->flags)));

  return err;
}
