/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */

#include "data_structure.h"
#include "state_variables.h"
#include "utils.h"


int State_variables::initialization(const size_t n_Fs,
                                    const size_t n_vars,
                                    const size_t n_flags)
{
  int err = 0;
  this->n_Fs = n_Fs;
  this->Fs = new Matrix<double>[n_Fs];

  for(size_t i = 0; i < n_Fs; i++)
  {
    this->Fs[i].initialization(3,3);
    this->Fs[i](1,1) = this->Fs[i](2,2) = this->Fs[i](3,3) = 1.0;
  }

  /* state variables is column vector */
  this->state_vars = new Matrix<double>[1];
  this->state_vars[0].initialization(n_vars,1);

  this->n_flags = n_flags;
  this->flags = new int[n_flags]();

  return err;
}

size_t State_variables::state_variables_get_packed_size(void)
{
  /* get the length of all Fs. Note that we are computing the length
     of each tensor just in case. */
  size_t len_Fs = 0;
  for (int i = 0, e = this->n_Fs; i < e; i++) {
    len_Fs += this->Fs[i].m_col * this->Fs[i].m_row;
  }
  const size_t len_vars = this->state_vars->m_col * this->state_vars->m_row;
  return ( (len_Fs + len_vars) * sizeof(double) + this->n_flags * sizeof(int) );
}

int State_variables::state_variables_pack(char *buffer,
                                          size_t *pos)
{
  int err = 0;

  /* pack the deformation gradients */
  for (size_t i = 0, e = this->n_Fs; i < e; i++) {
    const size_t len = this->Fs[i].m_row * this->Fs[i].m_col;
    const double *data = this->Fs[i].m_pdata;
    pack_data(data, buffer, pos, len, sizeof(*data));
  }

  /* pack the state variables */
  {
    const size_t len = this->state_vars->m_row * this->state_vars->m_col;
    const double *data = this->state_vars->m_pdata;
    pack_data(data, buffer, pos, len, sizeof(*data));
  }

  /* pack the flags */
  pack_data(this->flags, buffer, pos, this->n_flags, sizeof(*(this->flags)));

  return err;
}

int State_variables::state_variables_unpack(const char *buffer,
                                            size_t *pos)
{
  int err = 0;

  /* pack the deformation gradients */
  for (size_t i = 0, e = this->n_Fs; i < e; i++) {
    const size_t len = this->Fs[i].m_row * this->Fs[i].m_col;
    double *data = this->Fs[i].m_pdata;
    unpack_data(buffer, data, pos, len, sizeof(*data));
  }

  /* pack the state variables */
  {
    const size_t len = this->state_vars->m_row * this->state_vars->m_col;
    double *data = this->state_vars->m_pdata;
    unpack_data(buffer, data, pos, len, sizeof(*data));
  }

  /* pack the flags */
  unpack_data(buffer, this->flags, pos, this->n_flags, sizeof(*(this->flags)));

  return err;
}
