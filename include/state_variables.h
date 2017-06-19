/**
 * Interface for the State_variables object.
 *
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */
#pragma once
#ifndef STATE_VARIABLES_H
#define STATE_VARIABLES_H

#include <stdlib.h>
#include "data_structure.h"

using namespace gcm;

/**
 * Object for storing state variables at an integration point.
 */

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

class State_variables
{
  public:
  Matrix<double> *Fs;
  /** Array of handles to deformation gradients, e.g., Fp, Ft, Fe,... */

  /** Handle to vector of state variables */
  Matrix<double> *state_vars;
  int *flags;

  size_t n_Fs;
  size_t n_flags;
  
  State_variables()
  {
    Fs = NULL;
    state_vars = NULL;
    flags = NULL;
    n_Fs = 0;
    n_flags = 0;
  };

  ~State_variables()
  {
    delete this->Fs;
    this->Fs = NULL;
    this->n_Fs = 0;  

    delete this->state_vars;
    this->state_vars = NULL;

    delete this->flags;
    this->flags = NULL;
    this->n_flags = 0;
  }  

  void initialization(const size_t n_Fs,
                      const size_t n_vars,
                      const size_t n_flags);
                      
  size_t state_variables_get_packed_size(void);                     

  int state_variables_pack(char *buffer,
                           size_t *pos);
  int state_variables_unpack(const char *buffer,
                             size_t *pos);
};

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */


#endif
