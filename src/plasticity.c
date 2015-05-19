/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Adetokunbo Adedoyin, [1], <aadedoyi@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */

#include "plasticity.h"

#include "material.h"
#include "hommat.h"
#include "PGFEM_io.h"
/* #include "data_structure_c.h" */

/* model headers */
#include "plasticity_model_none.h"

int plasticity_construct(Plasticity *p)
{
  int err = 0;
  /* poison all values */
  p->Fs = NULL;
  /* p->state_vars = NULL; (is this possible?) */
  p->p_mat = NULL;
  p->p_hmat = NULL;
  p->integration_algorithm = NULL;
  p->compute_dev_stress = NULL;
  p->compute_pressure = NULL;
  p->compute_dev_tangent = NULL;
  p->compute_pressure_tangent = NULL;
  p->update_state_vars = NULL;
  p->reset_state_vars = NULL;
  p->n_Fs = 0;
  p->type = -1;
  return err;
}

int plasticity_initialize(Plasticity *p,
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

int plasticity_destroy(Plasticity *p)
{
  int err = 0;
  for ( size_t i = 0, e = p->n_Fs; i < e; i++) {
    /* destroy matrix objects */
  }
  free(p->Fs);
  p->Fs = NULL;

  /* destroy state variables vector */
  /* Matrix_destroy(p->state_vars); */

  /* drop pointer to material (material free'd elsewhere) */
  p->p_mat = NULL;
  p->p_hmat = NULL;

  /* drop function pointers */
  p->integration_algorithm = NULL;
  p->compute_dev_stress = NULL;
  p->compute_pressure = NULL;
  p->compute_dev_tangent = NULL;
  p->compute_pressure_tangent = NULL;
  p->update_state_vars = NULL;
  p->reset_state_vars = NULL;

  /* reset counters/flags */
  p->n_Fs = 0;
  p->type = -1;

  return err;
}
