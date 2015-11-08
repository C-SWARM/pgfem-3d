/**
 * Interface for the State_variables object.
 *
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Adetokunbo Adedoyin, [1], <aadedoyi@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */
#pragma once
#ifndef STATE_VARIABLES_H
#define STATE_VARIABLES_H

#include <stdlib.h>

struct Matrix_double;
#ifndef TYPE_MATRIX_DOUBLE
#define TYPE_MATRIX_DOUBLE
typedef struct Matrix_double Matrix_double;
#endif

#ifndef TYPE_VECTOR_DOUBLE
#define TYPE_VECTOR_DOUBLE
typedef struct Matrix_double Vector_double;
#endif
/**
 * Object for storing state variables at an integration point.
 */
struct State_variables {
  /** Array of handles to deformation gradients, e.g., Fp, Ft, Fe,... */
  Matrix_double *Fs;

  /** Handle to vector of state variables */
  Vector_double *state_vars;
  int *flags;

  size_t n_Fs;
  size_t n_flags;
};

#ifndef TYPE_STATE_VARIABLES
#define TYPE_STATE_VARIABLES
typedef struct State_variables State_variables;
#endif

int state_variables_build(State_variables *s);
int state_variables_destroy(State_variables *s);
int state_variables_initialize(State_variables *s,
                               const size_t n_Fs,
                               const size_t n_vars,
                               const size_t n_flags);

#endif
