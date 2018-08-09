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

class State_variables
{
  public:
  Matrix<double> *Fs;
  /** Array of handles to deformation gradients, e.g., Fp, Ft, Fe,... */

  /** Handle to vector of state variables */
  Matrix<double> *state_vars;
  int *flags;
  
  Matrix<double> pressure;
  Matrix<double> tFr;
  
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
    if(Fs)
    {
      delete[] Fs;  
      Fs = NULL;
    }
    
    n_Fs = 0;  
  
    if(state_vars)
    {  
      delete[] state_vars;
      state_vars = NULL;
    }
  
    if(flags)
    {
      delete[] flags;
      flags = NULL;
    }  
    n_flags = 0;    
  }  

  int initialization(const size_t n_Fs,
                      const size_t n_vars,
                      const size_t n_flags);
                      
  size_t state_variables_get_packed_size(void);                     

  int state_variables_pack(char *buffer,
                           size_t *pos);
  int state_variables_unpack(const char *buffer,
                             size_t *pos);
};

/// Object for querying/describing the state variables.
class Model_var_info
{
 public:

  char **F_names;
  char **var_names;
  char **flag_names;
  size_t n_Fs;
  size_t n_vars;
  size_t n_flags;

  /// construct a Model_var_info object.
  Model_var_info()
  {
    F_names    = NULL;
    var_names  = NULL;
    flag_names = NULL;
    n_Fs    = 0;
    n_vars  = 0;
    n_flags = 0;
  };

  /// destroy a Model_var_info object. Assumes full control of all
  /// internal pointers.
  ~Model_var_info();

  /// Print the object to the specified file.
  int print_variable_info(FILE *f);
};

/// set constitutive model information such as variable names.
typedef struct {
  int id;
  char name[1024];
} CMVariableNames;
 
int constitutive_model_info(Model_var_info &info,
                            const int Vars_no, const CMVariableNames *Vars,
                            const int Tens_no, const CMVariableNames *Tens,
                            const int Flag_no, const CMVariableNames *Flag);

#endif
