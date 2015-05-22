/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Adetokunbo Adedoyin, [1], <aadedoyi@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */

#include "constitutive_model.h"

#include "material.h"
#include "hommat.h"
#include "PGFEM_io.h"
#include "data_structure_c.h"

/* model headers */
#include "plasticity_model_none.h"

int constitutive_model_construct(Constitutive_model *m)
{
  int err = 0;
  m->model = NULL;
  err += state_variables_build(&(m->vars));
  return err;
}

int constitutive_model_initialize(Constitutive_model *m,
                                  const Model_parameters *param)
{
  int err = 0;
  if (param == NULL){
    err++;
  } else {
    m->model = param;
    Model_var_info *info = NULL;
    m->model->get_var_info(&info);
    err += state_variables_initialize(&(m->vars),info->n_Fs,info->n_vars);
    err += model_var_info_destroy(&info);
  }
  return err;
}

int constitutive_model_destroy(Constitutive_model *m)
{
  int err = 0;
  /* drop pointer to Model_parameters object (deallocated elsewhere) */
  m->model = NULL;
  err += state_variables_destroy(&(m->vars));
  return err;
}

int model_var_info_destroy(Model_var_info **info)
{
  int err = 0;
  for (size_t i=0, e=(**info).n_Fs; i < e; i++){
    free(F_names[i]);
  }
  free(F_names);

  for (size_t i=0, e=(**info).n_vars; i < e; i++){
    free(var_names[i]);
  }
  free(var_names);

  /* destroy and invalidate pointer */
  free(*info);
  *info = NULL;
  return err;
}

int model_parameters_construct(Model_parameters *p)
{
  int err = 0;
  /* poison all values */
  p->p_mat = NULL;
  p->p_hmat = NULL;
  p->integration_algorithm = NULL;
  p->compute_dev_stress = NULL;
  p->compute_dudj = NULL;
  p->compute_dev_tangent = NULL;
  p->compute_d2udj2 = NULL;
  p->update_state_vars = NULL;
  p->reset_state_vars = NULL;
  p->get_var_info = NULL;
  p->type = -1;
  return err;
}

int model_parameters_initialize(Model_parameters *p,
                          const MATERIAL *p_mat,
                          const HOMMAT *p_hmat,
                          const size_t type)
{
  int err = 0;
  p->p_mat = p_mat;
  p->p_hmat = p_hmat;
  p->type = type;
  switch(type) {
  case NONE:
    err += plasticity_model_none_initialize(p);
    break;
  case CRYSTAL_PLASTICITY:
  case BPA_PLASTICITY:
  default:
    PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",type);
    err++;
    break;
  }
  return err;
}

int model_parameters_destroy(Model_parameters *p)
{
  int err = 0;
  /* drop pointer to material (material free'd elsewhere) */
  p->p_mat = NULL;
  p->p_hmat = NULL;

  /* drop function pointers */
  p->integration_algorithm = NULL;
  p->compute_dev_stress = NULL;
  p->compute_dudj = NULL;
  p->compute_dev_tangent = NULL;
  p->compute_d2udj2 = NULL;
  p->update_state_vars = NULL;
  p->reset_state_vars = NULL;
  p->get_var_info = NULL;

  /* reset counters/flags */
  p->type = -1;

  return err;
}
