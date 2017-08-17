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
#include "data_structure_c.h"

/**
 * Object for storing state variables at an integration point.
 */

struct State_variables {
  /** Array of handles to deformation gradients, e.g., Fp, Ft, Fe,... */
  Matrix<double> *Fs;

  /** Handle to vector of state variables */
  Vector<double> *state_vars;
  int *flags;

  size_t n_Fs;
  size_t n_flags;
};

int state_variables_build(State_variables *s);
int state_variables_destroy(State_variables *s);
int state_variables_initialize(State_variables *s, const size_t n_Fs,
    const size_t n_vars, const size_t n_flags);
size_t state_variables_get_packed_size(const State_variables *s);
int state_variables_pack(const State_variables *s, char *buffer, size_t *pos);
int state_variables_unpack(State_variables *s, const char *buffer, size_t *pos);

#endif
