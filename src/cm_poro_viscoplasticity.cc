
/// Define poro-visco plasticity model functions for the constitutive model interface
/// 
/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
/// Alberto Salvadori, [1], <asalvad2@nd.edu>
/// 
/// [1] University of Notre Dame, Notre Dame, IN

//#include <ttl/ttl.h>
#include "cm_poro_viscoplasticity.h"

#include "plasticity_model.h"
#include "constitutive_model.h"
#include "cm_placeholder_functions.h"
#include "new_potentials.h"
#include "data_structure_c.h"
#include "elem3d.h"
#include "gen_path.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "material_properties.h"
#include "hyperelasticity.h"
#include "crystal_plasticity_integration.h"
#include "flowlaw.h"

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

#define MAX_D_GAMMA 0.005

static const int FLAG_end = 0;

enum variable_names {
  VAR_L_n,
  VAR_L_np1,
  VAR_g_n,
  VAR_g_np1,
  VAR_L_nm1,
  VAR_g_nm1,
  VAR_end
};

enum tensor_names {
  TENSOR_Fn,
  TENSOR_pFn,
  TENSOR_Fnp1,
  TENSOR_pFnp1,
  TENSOR_Fnm1,
  TENSOR_pFnm1,
  TENSOR_end
};

enum param_names {
  PARAM_pressure,
  PARAM_tol_hardening,
  PARAM_tol_M,
  PARAM_computer_zero,
  PARAM_NO
};

enum param_index_names {
  PARAM_max_itr_stag,
  PARAM_max_itr_hardening,
  PARAM_max_itr_M,
  PARAM_max_subdivision,
  PARAM_INX_NO
};

Define_Matrix(double);
Define_Matrix(int);

static size_t cm_pvp_get_size(const Constitutive_model *m)
{
  return state_variables_get_packed_size(m->vars_list[0]+m->model_id);
}

static int cm_pvp_pack(const Constitutive_model *m,
                    char *buffer,
                    size_t *pos)
{
  return state_variables_pack(m->vars_list[0]+m->model_id, buffer, pos);
}

static int cm_pvp_unpack(Constitutive_model *m,
                     const char *buffer,
                     size_t *pos)
{
  return state_variables_unpack(m->vars_list[0]+m->model_id, buffer, pos);
}

/// Private structure for use exclusively with this model and
// associated functions.
typedef struct poro_viscoplasticity_ctx {
  double *F;
  double dt;    // time increment
  double alpha; // mid point alpha
  double *eFnpa;
  int is_coulpled_with_thermal;
  double *hFn;
  double *hFnp1;  
} poro_viscoplasticity_ctx;

static int poro_viscoplasticity_int_alg(Constitutive_model *m,
                                        const void *ctx);
                                        
int poro_viscoplasticity_model_initialize(Model_parameters *p)
{
  int err = 0;

  /* set functions */
  p->integration_algorithm = NULL;//poro_viscoplasticity_int_alg;
  p->compute_dev_stress    = NULL;//poro_viscoplasticity_dev_stress;
  p->compute_dudj          = NULL;//poro_viscoplasticity_dudj;
  p->compute_dev_tangent   = NULL;//poro_viscoplasticity_dev_tangent;
  p->update_elasticity     = NULL;//poro_viscoplasticity_model_update_elasticity;
  p->compute_d2udj2        = NULL;//poro_viscoplasticity_d2udj2;
  p->update_state_vars     = NULL;//poro_viscoplasticity_update;
  p->reset_state_vars      = NULL;//poro_viscoplasticity_reset;
  p->get_subdiv_param      = NULL;//poro_viscoplasticity_model_get_subdiv_param;
  p->get_var_info          = NULL;//poro_viscoplasticity_info;
  p->get_F                 = NULL;//poro_viscoplasticity_get_F;
  p->get_Fn                = NULL;//poro_viscoplasticity_get_Fn;
  p->get_Fnm1              = NULL;//poro_viscoplasticity_get_Fnm1;  
  p->get_pF                = NULL;//poro_viscoplasticity_get_pF;
  p->get_pFn               = NULL;//poro_viscoplasticity_get_pFn;
  p->get_pFnm1             = NULL;//poro_viscoplasticity_get_pFnm1;    
  p->get_eF                = NULL;//poro_viscoplasticity_get_eF;
  p->get_eFn               = NULL;//poro_viscoplasticity_get_eFn;
  p->get_eFnm1             = NULL;//poro_viscoplasticity_get_eFnm1;
  p->get_eF_of_hF          = NULL;//poro_viscoplasticity_get_eF_with_thermal;

  p->reset_state_vars_using_temporal   = NULL;//poro_viscoplasticity_reset_using_temporal;
  p->update_np1_state_vars_to_temporal = NULL;//poro_viscoplasticity_update_np1_to_temporal;
  p->save_state_vars_to_temporal       = NULL;//poro_viscoplasticity_save_to_temporal;
    
  p->get_hardening         = NULL; //poro_viscoplasticity_get_hardening;
  p->get_hardening_nm1     = NULL; //poro_viscoplasticity_get_hardening_nm1;  
  p->get_plast_strain_var  = NULL; //poro_viscoplasticity_get_lam_p;
  p->write_restart         = NULL; //poro_viscoplasticity_write_restart;
  p->read_restart          = NULL; //poro_viscoplasticity_read_restart;
  p->destroy_ctx           = NULL; //poro_viscoplasticity_model_ctx_destroy;
  p->compute_dMdu          = NULL; //poro_viscoplasticity_compute_dMdu;
  p->set_init_vals         = NULL; //poro_viscoplastivity_set_init_vals;
  p->read_param            = NULL; //poro_viscoplastivity_read;
  p->get_size              = NULL; //poro_viscoplastivity_get_size;
  p->pack                  = NULL; //poro_viscoplastivity_pack;
  p->unpack                = NULL; //poro_viscoplastivity_unpack;
  p->type                  = CRYSTAL_PLASTICITY;

  p->n_param           = PARAM_NO;
  p->model_param       = calloc(PARAM_NO, sizeof(*(p->model_param)));
  p->n_param_index     = PARAM_INX_NO;
  p->model_param_index = calloc(PARAM_INX_NO, sizeof(*(p->model_param_index)));
  return err;
}