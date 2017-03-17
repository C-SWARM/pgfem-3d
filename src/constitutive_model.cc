/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */

#include "constitutive_model.h"
#include "cm_placeholder_functions.h"
#include "plasticity_model_none.h"
#include "plasticity_model.h"
#include "plasticity_model_BPA.h"
#include "cm_iso_viscous_damage.h"
#include "cm_j2_plasticity.h"
#include "cm_uqcm.h"

#include "hommat.h"
#include "PGFEM_io.h"
#include "PGFEM_mpi.h"
#include "data_structure_c.h"
#include "supp.h"
#include "elem3d.h"
#include "femlib.h"
#include "index_macros.h"
#include "material_properties.h" // <= constitutive model material properties
#include "hyperelasticity.h"     // <= constitutive model elasticity
#include <string.h>
#include "dynamics.h"
#include "PGFem3D_data_structure.h"
#include "get_dof_ids_on_elem.h"

#ifndef _Matrix_double
Define_Matrix(double);
#define _Matrix_double 1
#endif

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

#define MAX(a, b) ((a) >= (b)? (a) : (b))

/// compute nodal temperature in an element
///
/// if time step is subdivided, temperature needs to be subdivided
/// dT   = T(t(n+1)) - T(t(n))
/// Tn_e = T(t(n)) + dT*factor_n
/// Tn_e = T(t(n)) + dT*factor_n
///
/// T(t=n)---------Tn_e----Tnp1_e-----T(t=np1)
/// |<--factor_n-->|       |
/// |<-------factor_np1--->| 
///
/// \param[in] fe container of finite element resources for an element
/// \param[in] grid an object containing all mesh info
/// \param[in] fv_h field variable object for thermal 
/// \param[in] load object for loading
/// \param[in] mp_id mutiphysics id
/// \param[out] Tnp1_e computed nodal temperature for current element at t(n+1)
/// \param[out] Tn_e computed nodal temperature for current element at t(n+1)
/// \param[out] Tnm1_e computed nodal temperature for current element at t(n+1)
/// \param[in] factor_np1 factor of computing temperature for t(n+1) 
/// \param[in] factor_n factor of computing temperature for t(n)
/// \return non-zero on internal error
int get_nodal_temperatures(FEMLIB *fe,
                           GRID *grid,
                           FIELD_VARIABLES *fv_h,
                           LOADING_STEPS *load,
                           int mp_id,
                           double *Tnp1_e,
                           double *Tn_e,
                           double *Tnm1_e,
                           double factor_np1,
                           double factor_n)
{
  int err = 0;

  int nne = fe->nne;
  long *cn = aloc1l(nne);  
  long *nod = (fe->node_id).m_pdata;
  int ndofn = fv_h->ndofn;
  get_dof_ids_on_elem_nodes(0,fe->nne,ndofn,nod,grid->node,cn,mp_id);  

  double T0    = fv_h->u0;
  double *T    = fv_h->u_np1;
  double *dT   = fv_h->d_u;
  double *Tn   = fv_h->temporal->u_n;
  double *Tnm1 = fv_h->temporal->u_nm1;  
  
  SUPP sup = load->sups[mp_id];

  for(int ia=0; ia<nne; ia++)
  {
    Tnm1_e[ia] = Tnm1[nod[ia]];
    const int id = cn[ia];
    const int aid = abs(id) - 1;
    
    double T_n = Tn[nod[ia]];
    double T_np1 = 0.0;

    if (id == 0)
      T_np1 = T0;
    else if(id > 0)
      T_np1 = T[aid] + dT[aid];
    else
      T_np1 = T0 + sup->defl[aid] + sup->defl_d[aid];
   
      Tn_e[ia] = T_n + (T_np1 - T_n)*factor_n;
    Tnp1_e[ia] = T_n + (T_np1 - T_n)*factor_np1;  
  }
  
  free(cn);
  return err;
}

/// compute temperature dependent variables at the integration point
/// 
/// \param[in] fe container of finite element resources for an element
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in] T0 initial temperature
/// \param[in] Tnp1 nodal temperature array at t(n+1)
/// \param[in] Tn nodal temperature array at t(n)
/// \param[in] Tnm1 nodal temperature array at t(n-1)
/// \param[out] hFnp1 computed deformation gradient at t(n+1)
/// \param[out] hFn computed deformation gradient at t(n)
/// \param[out] hFnm1 computed deformation gradient at t(n-1)
/// \return non-zero on internal error 
int compute_temperature_at_ip(FEMLIB *fe,
                              const GRID *grid,
                              const MATERIAL_PROPERTY *mat,
                              const double T0,
                              double *Tnp1,
                              double *Tn,
                              double *Tnm1,
                              double *hFnp1,
                              double *hFn,
                              double *hFnm1)
{
  int err = 0;
  double T     = 0.0;
  double dTn   = 0.0;
  double dTnp1 = 0.0;
  double dTnm1 = 0.0;
  double *N = (fe->N).m_pdata;
  
  for(int ia=0; ia<fe->nne; ia++)
  {
    T     += N[ia]*Tnp1[ia];
    dTnp1 += N[ia]*(Tnp1[ia] - T0);
    dTn   += N[ia]*(Tn[ia]   - T0);
    dTnm1 += N[ia]*(Tn[ia]   - T0);    
  }
  
  const int eid = fe->curt_elem_id;
  const int hmat_id = (grid->element[eid]).mat[2];
  const int mat_id = (mat->hommat[hmat_id]).mat_id;
            
  const double ax = mat->mater[mat_id].ax;
  const double ay = mat->mater[mat_id].ay;
  const double az = mat->mater[mat_id].az;
  
  hFnp1[1] = hFnp1[2] = hFnp1[3] = hFnp1[5] = hFnp1[6] = hFnp1[7] = 0.0;
    hFn[1] =   hFn[2] =   hFn[3] =   hFn[5] =   hFn[6] =   hFn[7] = 0.0;
  hFnm1[1] = hFnm1[2] = hFnm1[3] = hFnm1[5] = hFnm1[6] = hFnm1[7] = 0.0;
  
  hFnp1[0] = 1.0 + ax*dTnp1;
  hFnp1[4] = 1.0 + ay*dTnp1;
  hFnp1[8] = 1.0 + az*dTnp1;
  
  hFn[0] = 1.0 + ax*dTn;
  hFn[4] = 1.0 + ay*dTn;
  hFn[8] = 1.0 + az*dTn;
  
  hFnm1[0] = 1.0 + ax*dTnm1;
  hFnm1[4] = 1.0 + ay*dTnm1;
  hFnm1[8] = 1.0 + az*dTnm1;  
  
  return err;
}

/* add the macroscopic deformation gradient to the _TOTAL_ deformation
   gradient */
static void cm_add_macro_F(const SUPP sup,
                           double * restrict F)
{
  const double * restrict F0 = sup->F0;
  for (int i = 0; i < DIM_3x3; i++) F[i] += F0[i];
}

/// construct Matrix(double) array
///
/// \param[out] F_out Matrix array to be created
/// \param[in] row number of rows of each Matrix
/// \param[in] col number of columns of each Matrix
/// \param[in] num size of array
/// \param[in] initialize if 1, initialize each matrix to zero
/// \return non-zero on internal error
int construct_matrix_array(Matrix(double) **F_out, 
                           const int row, 
                           const int col, 
                           int num, 
                           int initialize)
{
  int err = 0;

  Matrix(double) *F = NULL;
  F = malloc(num*sizeof(Matrix(double)));
  
  if(F==NULL)
  {  
    printf("Memory allocation error [%s:%s:%d]\n", __func__, __FILE__,__LINE__);
    err++;
    return err; 
  }  
  for(int ia = 0; ia < num; ia++)
  {
    if(initialize)
      Matrix_construct_init(double, F[ia],row,col,0.0);
    else
      Matrix_construct_redim(double, F[ia],row,col);      
  }
  
  *F_out = F; 
    
  return err;
}


/// cleanup Matrix(double) array
///
/// \param[out] F_in Matrix array to be deallocated
/// \param[in] num size of array
/// \return non-zero on internal error
int cleanup_matrix_array(Matrix(double) **F_in, int num)
{
  int err = 0;

  Matrix(double) *F = *F_in;
  for(int ia = 0; ia < num; ia++)
    Matrix_cleanup(F[ia]);

  free(F);
  
  *F_in = NULL; 
    
  return err;
}

/* this is a wrapper function for the switch that was copy/pasted
   everywhere. It is no big deal to keep adding to this private
   function's argument list. Just put everything any model might need
   and the switch will handle what is actually used. */
int construct_model_context(void **ctx,
                            const int type,
                            double *F,
                            const double dt,
                            const double alpha,
                            double *eFnpa)
{
  int err = 0;
  switch(type) {
  case TESTING:
  case HYPER_ELASTICITY:
    err += plasticity_model_none_ctx_build(ctx, F, eFnpa, NULL, NULL, 0);
    break;
  case CRYSTAL_PLASTICITY:
    err += plasticity_model_ctx_build(ctx, F, dt,alpha, eFnpa, NULL, NULL, 0);
    break;
  case BPA_PLASTICITY:
    err += plasticity_model_BPA_ctx_build(ctx, F, dt);
    break;
  case ISO_VISCOUS_DAMAGE:
    err += iso_viscous_damage_model_ctx_build(ctx, F, dt);
    break;
  case J2_PLASTICITY_DAMAGE:
    err += j2d_plasticity_model_ctx_build(ctx, F, dt);
    break;
  default:
    PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n", type);
    err++;
    break;
  }
  assert (err == 0);
  return err;
}

/// constructor of constitutive model context
///
/// this is a wrapper function for the switch that was copy/pasted
/// everywhere for the coupled problem. Keep add function's argument list to this private.
/// In adding more constitutive model, just put everything any
/// and the switch will handle what is actually used.
/// 
/// \param[out] ctx constructed context based on the model type
/// \param[in] type constitutive model type
/// \param[in] F total deformation gradient
/// \param[in] dt time step size
/// \param[in] alpha mid point alpha
/// \param[in] eFnpa mid point elastic part of deformation gradient
/// \param[in] hFn thermal part of deformation gradient at t(n)
/// \param[in] hFnp1 thermal part of deformation gradient at t(n+1)
/// \return non-zeoro on internal error
int construct_model_context_with_thermal(void **ctx,
                                         const int type,
                                         double *F,
                                         const double dt,
                                         const double alpha,
                                         double *eFnpa,
                                         double *hFn,
                                         double *hFnp1)
{
  int err = 0;
  switch(type) {
  case TESTING:
  case HYPER_ELASTICITY:
    err += plasticity_model_none_ctx_build(ctx, F, eFnpa, hFn, hFnp1, 1);
    break;  	  		
  case CRYSTAL_PLASTICITY:
    err += plasticity_model_ctx_build(ctx, F, dt,alpha, eFnpa, hFn, hFnp1, 1);
    break;
  default:
    PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n", type);
    err++;
    break;
  }
  assert (err == 0);
  return err;
}

int constitutive_model_construct(Constitutive_model *m)
{
  int err = 0;
  m->param = NULL;
  err += state_variables_build(m->vars_list[0] + m->model_id);
  return err;
}

int constitutive_model_initialize(Constitutive_model *m,
                                  const Model_parameters *p)
{
  int err = 0;
  if (p == NULL){
    err++;
  } else {
    m->param = p;
    Model_var_info *info = NULL;
    m->param->get_var_info(&info);
    err += state_variables_initialize(m->vars_list[0] + m->model_id, info->n_Fs,
                                      info->n_vars, info->n_flags);
    m->param->set_init_vals(m);
    err += model_var_info_destroy(&info);
  }
  return err;
}

int constitutive_model_destroy(Constitutive_model *m)
{
  int err = 0;
  /* drop pointer to Model_parameters object (deallocated elsewhere) */
  m->param = NULL;
  err += state_variables_destroy(m->vars_list[0] + m->model_id);
  return err;
}

int model_var_info_print(FILE *f,
                         const Model_var_info * info)
{
  int err = 0;
  fprintf(f,"F names: ");
  for(int i = 0, e = info->n_Fs; i < e; i++) fprintf(f,"%s ",info->F_names[i]);
  fprintf(f,"\nVar names: ");
  for(int i = 0, e = info->n_vars; i < e; i++) fprintf(f,"%s ",info->var_names[i]);
  fprintf(f,"\nFlag names: ");
  for(int i = 0, e = info->n_flags; i < e; i++) fprintf(f,"%s ",info->flag_names[i]);
  fprintf(f,"\n");
  return err;
}

int model_var_info_destroy(Model_var_info **info)
{
  int err = 0;
  Model_var_info *t_info = *info;
  /* invalidate pointer */
  *info = NULL;

  /* deallocate internal memory */
  for (size_t i=0, e=t_info->n_Fs; i < e; i++){
    if(t_info->F_names[i])
      free(t_info->F_names[i]);
  }
  if(t_info->F_names)
    free(t_info->F_names);

  for (size_t i=0, e=t_info->n_vars; i < e; i++){
    if(t_info->var_names[i])
      free(t_info->var_names[i]);
  }
  if(t_info->var_names)
    free(t_info->var_names);

  for (size_t i=0, e=t_info->n_flags; i < e; i++){
    if(t_info->flag_names[i])
      free(t_info->flag_names[i]);
  }
  if(t_info->flag_names)
    free(t_info->flag_names);

  /* destroy memory for structure */
  free(t_info);

  return err;
}


int model_parameters_construct(Model_parameters *p)
{
  int err = 0;
  /* poison all values */
  p->p_hmat = NULL;
  p->cm_mat = NULL;
  p->cm_elast = NULL;
  p->integration_algorithm = NULL;
  p->compute_dev_stress = NULL;
  p->compute_dudj = NULL;
  p->compute_dev_tangent = NULL;
  p->compute_d2udj2 = NULL;
  p->compute_AST = NULL;
  p->update_elasticity = NULL;
  p->update_state_vars = NULL;
  p->reset_state_vars = NULL;
  p->reset_state_vars_using_temporal = NULL;
  p->update_np1_state_vars_to_temporal = NULL;
  p->save_state_vars_to_temporal = NULL;  

  p->get_var_info = NULL;
  p->get_F        = NULL;
  p->get_Fn       = NULL;
  p->get_Fnm1     = NULL;
  p->get_pF       = NULL;
  p->get_pFn      = NULL;
  p->get_pFnm1    = NULL;
  p->get_eF       = NULL;
  p->get_eFn      = NULL;
  p->get_eFnm1    = NULL;
  p->get_eF_of_hF = NULL;

  p->get_hardening        = NULL;
  p->get_hardening_nm1    = NULL;
  p->get_plast_strain_var = NULL;
  p->get_subdiv_param = cm_no_subdiv;

  p->write_restart = NULL;
  p->read_restart  = NULL;

  p->destroy_ctx  = NULL;
  p->compute_dMdu = NULL;

  p->set_init_vals = NULL;
  p->read_param = NULL;

  p->get_size = NULL;
  p->pack = NULL;
  p->unpack = NULL;

  p->type              = -1;
  p->n_param           = -1;
  p->model_param       = NULL;
  p->n_param_index     = -1;
  p->model_param_index = NULL;
  return err;
}

int model_parameters_initialize(Model_parameters *p,
                                const HOMMAT *p_hmat,
                                const size_t type)
{
  int err = 0;
  p->p_hmat = p_hmat;
  p->type = type;
  
  MATERIAL_CONSTITUTIVE_MODEL *cm_mat = malloc(sizeof(MATERIAL_CONSTITUTIVE_MODEL));
  MATERIAL_ELASTICITY          *mat_e = malloc(sizeof(MATERIAL_ELASTICITY));
  ELASTICITY *elast = malloc(sizeof(ELASTICITY));
    
  set_properties_using_E_and_nu(mat_e,p_hmat->E,p_hmat->nu);
  mat_e->m01 = p_hmat->m01;
  mat_e->m10 = p_hmat->m10;
  mat_e->G   = p_hmat->G;
  mat_e->kappa = hommat_get_kappa(p_hmat);
  mat_e->devPotFlag = p_hmat->devPotFlag;
  mat_e->volPotFlag = p_hmat->volPotFlag;  

  set_properties_constitutive_model(cm_mat,mat_e,NULL);      
  construct_elasticity(elast, mat_e, 1);
  
  p->cm_mat   = cm_mat;
  p->cm_elast = elast;  
  
  switch(type) {
  case TESTING:
  case HYPER_ELASTICITY:
    err += plasticity_model_none_initialize(p);
    break;
  case CRYSTAL_PLASTICITY:
  {
    err += plasticity_model_initialize(p);
    break;
  }
  case BPA_PLASTICITY:
    err += plasticity_model_BPA_initialize(p);
    break;
  case ISO_VISCOUS_DAMAGE:
    err += iso_viscous_damage_model_initialize(p);
    break;
  case J2_PLASTICITY_DAMAGE:
    err += j2d_plasticity_model_initialize(p);
    break;
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

  switch(p->type) {
  case TESTING:
    break;
  case HYPER_ELASTICITY:    
    break;
  case CRYSTAL_PLASTICITY:
    err += plasticity_model_destory(p);
    break;
  case BPA_PLASTICITY:    
    break;
  case ISO_VISCOUS_DAMAGE:
    break;
  case J2_PLASTICITY_DAMAGE:
    break;
  default:
    PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",p->type);
    err++;
    break;
  }

  /* drop pointer to material (material free'd elsewhere) */
  p->p_hmat = NULL;    
  free((p->cm_mat)->mat_e);
  free(p->cm_mat);  
  destruct_elasticity(p->cm_elast);
  free(p->cm_elast);
  
  /* drop function pointers */
  p->integration_algorithm = NULL;
  p->compute_dev_stress = NULL;
  p->compute_dudj = NULL;
  p->compute_dev_tangent = NULL;
  p->compute_d2udj2 = NULL;
  p->compute_AST = NULL;
  p->update_state_vars = NULL;
  p->reset_state_vars = NULL;
  p->reset_state_vars_using_temporal = NULL;
  p->update_np1_state_vars_to_temporal = NULL;
  p->save_state_vars_to_temporal = NULL;
  p->get_var_info = NULL;
  p->get_F        = NULL;  
  p->get_Fn       = NULL;
  p->get_Fnm1     = NULL;
  p->get_pF       = NULL;
  p->get_pFn      = NULL;
  p->get_pFnm1    = NULL;
  p->get_eF       = NULL;
  p->get_eFn      = NULL;
  p->get_eFnm1    = NULL;
  p->get_eF_of_hF = NULL;

  p->get_hardening = NULL;
  p->get_hardening_nm1 = NULL;
  p->get_plast_strain_var = NULL;
  p->get_subdiv_param = NULL;

  p->write_restart = NULL;
  p->read_restart = NULL;

  p->destroy_ctx = NULL;
  p->compute_dMdu = NULL;

  p->set_init_vals = NULL;
  p->read_param = NULL;

  p->get_size = NULL;
  p->pack = NULL;
  p->unpack = NULL;

  /* reset counters/flags */
  p->type = -1;

  /* free model constants */
  p->n_param = -1;
  free(p->model_param);
  p->model_param = NULL;
  p->n_param_index = -1;
  p->model_param_index = NULL;

  return err;
}

// compute stiffness tensor
int constitutive_model_default_update_elasticity(const Constitutive_model *m,
                                                 Matrix_double *eF,
                                                 Matrix_double *L,
                                                 Matrix_double *S,
                                                 const int compute_stiffness)
{
  int err = 0;
 
  ELASTICITY *elast = (m->param)->cm_elast; // get elasticity handle
  double *tempS = elast->S; // temporal pointer to update *L, and *S using elast
  double *tempL = elast->L;
  elast->S = S->m_pdata;
  
  if(compute_stiffness)
    elast->L = L->m_pdata;
  else
    elast->L = NULL;
  
  elast->update_elasticity(elast,eF->m_pdata,compute_stiffness);

  elast->S = tempS;
  elast->L = tempL;
  return err;
}

static int compare_mat_id(const void *a, const void *b)
{
  return ((const HOMMAT *)a)->mat_id - ((const HOMMAT *)b)->mat_id;
}

int read_model_parameters_list(Model_parameters **param_list,
                               const int n_mat,
                               const HOMMAT *hmat_list,
                               FILE *in)
{
  /* see issue #28 */
  /* File format:
     ------
     num_entries
     material_id model_type
     { # begin model info
     ... # specified by model
     } # end model info
     ...
     ------
     caveats:
     - Comments can only go at the end of a line after all of the data
       on that line has been specified
     - No support for sub-braces in model sections (yet)
     - Undefined behavior for duplicate entires (will try to
       re-initialize the object)
   */
  int err = 0;
  if (n_mat <= 0) return 1;
  (*param_list) = malloc(n_mat*sizeof(**param_list));
  int *is_set = calloc(n_mat, sizeof(*is_set));

  int num_entries = -1;
  HOMMAT *key = calloc(1, sizeof(*key));
  err += scan_for_valid_line(in);
  fscanf(in, "%d", &num_entries);

  int i = 0;
  for (i = 0; i < num_entries; i++) {
    int model_type = -1;
    HOMMAT *p_hmat = NULL;
    err += scan_for_valid_line(in);
    if (feof(in)) break;

    fscanf(in, "%d %d", &(key->mat_id), &model_type);
    err += scan_for_valid_line(in);

    int brace = fgetc(in);
    assert(brace == '{' && "Expect opening brace as next valid entry");

    /*
     * NOTE: The material ID in the input files is not necessarily the
     * index of the HOMMAT material. The hmat_list is the reduced set
     * of material properties that are actually used on the
     * domain. Therefore, we need to search for the matching
     * mat_id. Futheremore, no warning is issued if no match is found,
     * but we perform checks to ensire that all materials are
     * specified.
     */

    /* search for matching pointer in hmat_list (assume unique) */
    p_hmat = bsearch(key,hmat_list,n_mat,sizeof(*hmat_list),compare_mat_id);

    /* check for match */
    if (p_hmat != NULL) {
      int idx = p_hmat - hmat_list;
      if (!is_set[idx]) {
        is_set[idx] = 1;
        /* construct and initialize this object */
        
        err += model_parameters_construct(&((*param_list)[idx]) );
        if(model_type==CM_UQCM)
        {
          ((*param_list) + idx)->uqcm = 1;
          err += scan_for_valid_line(in);
          fscanf(in, "%d", &model_type);
        }
        else
          ((*param_list) + idx)->uqcm = 0;
          
        err += model_parameters_initialize(&((*param_list)[idx]),
                                           hmat_list + idx,
                                           model_type);
                                                   
        err += ((*param_list)[idx]).read_param(&((*param_list)[idx]),in);
        ((*param_list)[idx]).mat_id = key->mat_id;
      }
    }

    /* scan to closing brace and continue on to the next entry */
    while(fgetc(in) != '}' && !feof(in)){}
    if (feof(in)) break;
  }

  if (feof(in) && i != num_entries) {
    err++;
    assert(0 && "Prematurely reached EOF");
  }

  int sum = 0;
  for (int i = 0; i < n_mat; i++){
    sum += is_set[i];
  }
  if (sum != n_mat) err++;
  assert(sum == n_mat && "require that all model params are set");

  free(key);
  free(is_set);
  return err;
}

int destroy_model_parameters_list(const int n_mat,
                                  Model_parameters *param_list)
{
  int err = 0;
  if (param_list == NULL) return 0;

  for (int i = 0; i < n_mat; i++) {
    err += model_parameters_destroy(param_list + i);
  }
  free(param_list);
  return 0;
}                                

int init_all_constitutive_model(EPS *eps,
                                const int ne,
                                const ELEMENT *elem,
                                const int n_mat,
                                const Model_parameters *param_list)
{
  int err = 0;
  if (ne <= 0) return 1;
  if(param_list==NULL)
    return err;
  
  for(int i = 0; i < ne; i++) {
    /* aliases */
    EPS *p_eps = eps + i;
    const ELEMENT *p_el = elem + i;
    const Model_parameters *p_param = param_list + (p_el->mat[2]);

    long n_ip = 0;
    int_point(p_el->toe,&n_ip);
    for (int j = 0; j < n_ip; j++)
      err += constitutive_model_initialize((p_eps->model) + j, p_param);
  }
  
  plasticity_model_set_orientations(eps, ne, elem, n_mat, param_list); // nothing will happen if there is no use of the crystal plasticity model
  return err;
}

int constitutive_model_reset_state(EPS *eps,
                                   const int ne,
                                   const ELEMENT *elem)
{
  int err = 0;
  if (ne <= 0) return 1;

  for (int i = 0; i < ne; i++) {
    const ELEMENT *p_el = elem + i;
    long n_ip = 0;
    int_point(p_el->toe,&n_ip);
    for (int j = 0; j < n_ip; j++) {
      Constitutive_model *m = &(eps[i].model[j]);
      m->param->reset_state_vars(m);
    }
  }
  return err;
}

/// save state variables 
///
/// save state variables(t(n-1) and t(n)) to temporal variables(t(n-1) and t(n)) 
/// in order to use when solution step is failed and requires go to initial step
/// 
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \return non-zero on internal error
int constitutive_model_save_state_vars_to_temporal(FIELD_VARIABLES *fv,
                                                   GRID *grid)
{
  int err = 0;
  State_variables *var = fv->temporal->var;
  
  int intg_order = 0;
  const ELEMENT *elem = grid->element;
  
  for(int eid=0; eid<grid->ne; eid++)
  {
    int nne = elem[eid].toe;
    long nint = FEMLIB_determine_integration_type(nne, intg_order);
    for (int ip = 0; ip < nint; ip++) 
    {
      Constitutive_model *m = &(fv->eps[eid].model[ip]);
      m->param->save_state_vars_to_temporal(m, var+m->model_id);
    }
  }
  return err;
}

/// update state variables(t(n+1)) 
///
/// temporary save state variables(t(n+1)) 
/// when coupled physics calls dependent physics, 
/// original soultions at t(n-1) and t(n) and updated solutions t(n+1) are needed.
/// 
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \return non-zero on internal error
int constitutive_model_update_np1_state_vars_to_temporal(FIELD_VARIABLES *fv,
                                                         GRID *grid)
{
  int err = 0;
  State_variables *var = fv->temporal->var;
  
  int intg_order = 0;
  const ELEMENT *elem = grid->element;
  
  for(int eid=0; eid<grid->ne; eid++)
  {
    int nne = elem[eid].toe;
    long nint = FEMLIB_determine_integration_type(nne, intg_order);
    for (int ip = 0; ip < nint; ip++) 
    {
      Constitutive_model *m = &(fv->eps[eid].model[ip]);
      m->param->update_np1_state_vars_to_temporal(m, var+m->model_id);
    }
  }
  return err;
}

/// reset state variables using priori stored values
///
/// reset state variables(t(n-1) and t(n)) using temporal variables(t(n-1) and t(n))
/// when solution step is failed and requires go to initial step
/// 
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \return non-zero on internal error
int constitutive_model_reset_state_using_temporal(FIELD_VARIABLES *fv,
                                                  GRID *grid)
{
  int err = 0;
  State_variables *var = fv->temporal->var;
  
  int intg_order = 0;
  const ELEMENT *elem = grid->element;
  
  for(int eid=0; eid<grid->ne; eid++)
  {
    int nne = elem[eid].toe;
    long nint = FEMLIB_determine_integration_type(nne, intg_order);
    for (int ip = 0; ip < nint; ip++) 
    {
      Constitutive_model *m = &(fv->eps[eid].model[ip]);
      m->param->reset_state_vars_using_temporal(m, var+m->model_id);
    }
  }
  return err;
}

int constitutive_model_update_time_steps(const ELEMENT *elem,
                                         NODE *node,
                                         EPS *eps,
                                         const int ne,
                                         const int nn,
                                         const int ndofn,
                                         const double* r,
                                         const double dt,
                                         const int total_Lagrangian,
                                         const int mp_id)
{
  int nsd = 3;
  int err = 0;
  if (ne <= 0) return 1;
    
  for (int i = 0; i < ne; i++) 
  {
    const ELEMENT *p_el = elem + i;
    long n_ip = 0;
    int_point(p_el->toe,&n_ip);
    for (int j = 0; j < n_ip; j++)     
    {
      Constitutive_model *m = &(eps[i].model[j]);
      m->param->update_state_vars(m);
    }
  }    
  
  /*********************/
  /* Coordinate update */
  /*********************/
  if(total_Lagrangian) {
    for(int n = 0;n<nn; n++) {
      for(int a=0;a<nsd;a++) {
        int II = node[n].id_map[mp_id].id[a];
        if (II != 0){
          if (a == 0)      node[n].x1 = node[n].x1_fd + r[n*ndofn + a];
          else if (a == 1) node[n].x2 = node[n].x2_fd + r[n*ndofn + a];
          else if (a == 2) node[n].x3 = node[n].x3_fd + r[n*ndofn + a];
        }
      }
    }/* end n < nn */
  } else {
    for(int n = 0;n<nn; n++) {
      for(int a=0;a<nsd;a++) {
        int II = node[n].id_map[mp_id].id[a];
        if (II != 0){
          if (a == 0)      node[n].x1 += r[n*ndofn + a];
          else if (a == 1) node[n].x2 += r[n*ndofn + a];
          else if (a == 2) node[n].x3 += r[n*ndofn + a];
        }
      }
    }/* end n < nn */
  }
  
  return err;  
}

int constitutive_model_test(const HOMMAT *hmat, Matrix_double *L_in, int Print_results)
{
  int err = 0;//plasticity_model_test(hmat, L_in, Print_results);
  test_crystal_plasticity_single_crystal();
  return err;
}

/// Common part of computing stiffness Matrix
///
/// \param[out] lk computed element stiffness matrix
/// \param[in] fe finite element helper object
/// \param[in] Fr 2nd order tensor Fr  
/// \param[in] eFnMT 2nd order tensor (eFn*M)'
/// \param[in] eFn 2nd order tensor eFn
/// \param[in] M 2nd order tensor M
/// \param[in] FrTFr 2nd order tensor Fr'*Fr
/// \param[in] eFnM 2nd order tensor eFn*M
/// \param[in] S 2nd order tensor S
/// \param[in] L 4th order elasticity tensor L
/// \param[in] dMdu_all list of 2nd order dMdu tensors
/// \param[in] Jn det(Fn)
/// \return non-zero on internal error 
int compute_stiffness_matrix(double *lk,
                                    const FEMLIB *fe,
                                    const Matrix(double) *Fr,
                                    const Matrix(double) *eFnMT,
                                    const Matrix(double) *eFn,
                                    const Matrix(double) *M,
                                    const Matrix(double) *FrTFr,
                                    const Matrix(double) *eFnM,
                                    const Matrix(double) *S,
                                    const Matrix(double) *L,
                                    double *dMdu_all,
                                    const double Jn)
{
  int err = 0;
  const int nne = fe->nne;
  const int nsd = fe->nsd;
  
  Matrix(double) ST_ab, ST_wg, dMdu; // no memory is created
  ST_ab.m_row = ST_ab.m_col = DIM_3;
  ST_wg.m_row = ST_wg.m_col = DIM_3;
   dMdu.m_row =  dMdu.m_col = DIM_3;
    
  enum {MTeFnT_sAA_eFn,
        MTeFnT_sAA_eFnM,
        sBB,sCC,
        dCdu,
        temp_F_1,
        temp_F_2,
        Fend};
  
  // list of second-order tensors
  Matrix(double) *F2;
  err + construct_matrix_array(&F2,DIM_3,DIM_3,Fend,0); // if not initialized some compiler sets non-numbers and causes
                                                        // problems when do Matrix_Tns2_AxBxC
  for(int ia=0; ia<DIM_3x3; ia++)
  {
     F2[temp_F_1].m_pdata[ia] = 0.0;
     F2[temp_F_2].m_pdata[ia] = 0.0;
     F2[MTeFnT_sAA_eFn].m_pdata[ia] = 0.0;
     F2[MTeFnT_sAA_eFnM].m_pdata[ia] = 0.0;     
  }                                                          
  
  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const int id_ab = idx_4_gen(a,b,0,0,nne,nsd,nsd,nsd);
      ST_ab.m_pdata = (fe->ST)+id_ab;

      int AA  = temp_F_1; 
      int sAA = temp_F_2;
      Matrix_AxB(F2[AA],1.0,0.0,*Fr,1,ST_ab,0);
      Matrix_symmetric(F2[AA],F2[sAA]);
      
      Matrix_Tns2_AxBxC(F2[MTeFnT_sAA_eFn],1.0,0.0,*eFnMT,F2[sAA],*eFn);        
      Matrix_AxB(F2[MTeFnT_sAA_eFnM],1.0,0.0,F2[MTeFnT_sAA_eFn],0,*M,0);

      for(int w=0; w<nne; w++)
      {
        for(int g=0; g<nsd; g++)
        { 
          const int id_wg = idx_4_gen(w,g,0,0,nne,nsd,nsd,nsd);
          ST_wg.m_pdata = (fe->ST)+id_wg;

          dMdu.m_pdata = dMdu_all + id_wg;

          int BB = temp_F_1;
          Matrix_AxB(F2[BB],1.0,0.0,*Fr,1,ST_wg,0); 
          Matrix_symmetric(F2[BB],F2[sBB]);

          int CC = temp_F_1;
          Matrix_AxB(F2[CC], 1.0,0.0,ST_ab,1,ST_wg,0);
          Matrix_symmetric(F2[CC],F2[sCC]);

          // compute F2[dCdu]
          int MTeFnT_FrTFreFn     = temp_F_1;
          int MTeFnT_FrTFreFndMdu = temp_F_2;
          Matrix_Tns2_AxBxC(F2[MTeFnT_FrTFreFn],1.0,0.0,*eFnMT,*FrTFr,*eFn);            
          Matrix_AxB(F2[MTeFnT_FrTFreFndMdu],1.0,0.0,F2[MTeFnT_FrTFreFn],0,dMdu,0);                                    
          Matrix_symmetric(F2[MTeFnT_FrTFreFndMdu],F2[dCdu]);
                      
          Matrix_Tns2_AxBxC(F2[dCdu],1.0,1.0,*eFnMT,F2[sBB],*eFnM);            
          
          // compute F2[MTeFnT_sAA_eFnM]:L:F2[dCdu]
          int L_dCdu = temp_F_1;
          Matrix_Tns4_dd_Tns2(F2[L_dCdu],*L,F2[dCdu]);
          double MTeFnT_sAA_eFnM_L_dCdu = 0.0;
          Matrix_ddot(F2[MTeFnT_sAA_eFnM],F2[L_dCdu],MTeFnT_sAA_eFnM_L_dCdu);
          
          // compute F2[MTeFnT_sCC_eFnM]
          int MTeFnT_sCC_eFnM = temp_F_1;
          Matrix_Tns2_AxBxC(F2[MTeFnT_sCC_eFnM],1.0,0.0,*eFnMT,F2[sCC],*eFnM);
                      
          // compute F2[MTeFnT_sCC_eFnM]:F2[S]
          double MTeFnT_sCC_eFnM_S = 0.0;
          Matrix_ddot(F2[MTeFnT_sCC_eFnM],*S,MTeFnT_sCC_eFnM_S);
          
          // compute F2[MTeFnT_sAA_eFndMdu]
          int MTeFnT_sAA_eFndMdu  = temp_F_1;
          int sMTeFnT_sAA_eFndMdu = temp_F_2;
          Matrix_AxB(F2[MTeFnT_sAA_eFndMdu],1.0,0.0,F2[MTeFnT_sAA_eFn],0,dMdu,0);    
          Matrix_symmetric(F2[MTeFnT_sAA_eFndMdu], F2[sMTeFnT_sAA_eFndMdu]);        

          // compute F2[MTeFnT_sAA_eFndMdu]:F2[S]
          double sMTeFnT_sAA_eFndMdu_S = 0.0;            
          Matrix_ddot(F2[sMTeFnT_sAA_eFndMdu],*S,sMTeFnT_sAA_eFndMdu_S);
          
          const int lk_idx = idx_K(a,b,w,g,nne,nsd);  
                    
          lk[lk_idx] += 1.0/Jn*fe->detJxW*(MTeFnT_sAA_eFnM_L_dCdu + 2.0*sMTeFnT_sAA_eFndMdu_S + MTeFnT_sCC_eFnM_S);            
        }
      }
    }
  }
  
  // destroy second-order tensors
  err += cleanup_matrix_array(&F2, Fend);
  
  return err;
}

/// Common part of computing residual vector
///
/// \param[out] f computed residual vector
/// \param[in] fe finite element helper object
/// \param[in] Fr 2nd order tensor Fr  
/// \param[in] eFnMT 2nd order tensor (eFn*M)'
/// \param[in] eFnM 2nd order tensor eFn*M
/// \param[in] S 2nd order tensor S
/// \param[in] Jn det(Fn)
/// \return non-zero on internal error 
int compute_residual_vector(double *f,
                            const FEMLIB *fe,
                            const Matrix(double) *Fr,
                            const Matrix(double) *eFnMT,
                            const Matrix(double) *eFnM,
                            const Matrix(double) *S,
                            const double Jn)
{
  int err = 0;
  const int nne = fe->nne;
  const int nsd = fe->nsd;
  
  Matrix(double) ST_ab; // no memory is created
  ST_ab.m_row = ST_ab.m_col = DIM_3;

  enum {temp_F_1, temp_F_2, Fend};
  
  // list of second-order tensors
  Matrix(double) *F2;
  err + construct_matrix_array(&F2,DIM_3,DIM_3,Fend,1); 

  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const int id_ab = idx_4_gen(a,b,0,0,nne,nsd,nsd,nsd);
      ST_ab.m_pdata = (fe->ST)+id_ab;
      
      int AA  = temp_F_1; 
      int sAA = temp_F_2;
      Matrix_AxB(F2[AA],1.0,0.0,*Fr,1,ST_ab,0);
      Matrix_symmetric(F2[AA],F2[sAA]);
      
      int MTeFnT_sAA_eFnM = temp_F_1;
      Matrix_Tns2_AxBxC(F2[MTeFnT_sAA_eFnM],1.0,0.0,*eFnMT,F2[sAA],*eFnM);
      
      double MTeFnT_sAA_eFnM_S = 0.0; 
      Matrix_ddot(F2[MTeFnT_sAA_eFnM],*S,MTeFnT_sAA_eFnM_S);       
            
      int fe_id = a*nsd + b;              
      f[fe_id] += 1.0/Jn*fe->detJxW*MTeFnT_sAA_eFnM_S;              
    }
  }
  
  // destroy second-order tensors
  err += cleanup_matrix_array(&F2, Fend);
       
  return err;
}

/// compute ouput variables e.g. effective stress and strain
///
/// Visit each element and compute output variables according to the element model type.
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[in] load object for loading
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \param[in] alpha mid point rule alpha
/// \return non-zero on internal error
int constitutive_model_update_output_variables(GRID *grid,
                                               MATERIAL_PROPERTY *mat,
                                               FIELD_VARIABLES *FV,
                                               LOADING_STEPS *load,
                                               PGFem3D_opt *opts,
                                               MULTIPHYSICS *mp,
                                               int mp_id,
                                               const double dt,
                                               double alpha)
{
  int err = 0;
 
  int total_Lagrangian = 1; 
  if(opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;
 
  FIELD_VARIABLES *fv = FV+mp_id;   
  SIG *sig = fv->sig;
  EPS *eps = fv->eps;
  NODE *node = grid->node;
  ELEMENT *elem = grid->element;  
    
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  
  for(int ia=0; ia<fv->n_coupled; ia++)
  { 
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;    
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }     
    
  //static const double eye[DIM_3x3] = {[0] = 1.0, [4] = 1.0, [8] = 1.0};
  static const double eye[DIM_3x3] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  /* deformation gradient */
  Matrix_double F, eF, pF, S;
  Matrix_construct_redim(double, F,  DIM_3, DIM_3);
  Matrix_construct_init( double, eF,  DIM_3, DIM_3,0.0);
  Matrix_construct_redim(double, pF, DIM_3, DIM_3);
  Matrix_construct_redim(double, S,  DIM_3, DIM_3);
  
  MATERIAL_ELASTICITY mat_e_new;
  MATERIAL_ELASTICITY *mat_e_in;

  for (int i = 0; i < grid->ne; i++)
  {
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe,i,elem,node,0,total_Lagrangian);

    memset(sig[i].el.o,0,6*sizeof(double));
    memset(eps[i].el.o,0,6*sizeof(double));
    double volume = 0.0;


    Matrix(double) Tnp1, Tn, Tnm1;    
    FIELD_VARIABLES *fv_h = NULL;
    
    if(is_it_couple_w_thermal >= 0)
    {
      fv_h = fv->fvs[is_it_couple_w_thermal]; 
      Matrix_construct_init(double, Tnm1,fe.nne,1,0.0);    
      Matrix_construct_init(double, Tnp1,fe.nne,1,0.0); 
      Matrix_construct_init(double, Tn,  fe.nne,1,0.0);
    
      // compute temperature for this element for each nodal point
      int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
      err += get_nodal_temperatures(&fe, grid, fv_h, load, mp_cp_id,
                                    Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,fv->subdivision_factor_np1,fv->subdivision_factor_n);
    } 
    
    for(int ip=0; ip<fe.nint; ip++)
    {
      FEMLIB_elem_basis_V(&fe, ip+1); 
      Constitutive_model *m = &(eps[i].model[ip]);
      const Model_parameters *func = m->param;

      err += func->get_Fn(m, &F);
      err += func->get_pFn(m, &pF);
 
      void *ctx = NULL; 
      if(is_it_couple_w_thermal>=0)
      {
        // compute temperature at the integration point
        double T0 = fv_h->u0;
        double hFnm1[9],hFn[9],hFnp1[9];      
        err += compute_temperature_at_ip(&fe,grid,mat,T0,
                                         Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                         hFnp1,hFn,hFnm1);
        Matrix(double) pFI,hFI;
        Matrix_construct_redim(double, pFI,DIM_3,DIM_3);
        Matrix_construct_redim(double, hFI,DIM_3,DIM_3);
        Matrix_inv(pF,pFI);
        inv3x3(hFnp1,hFI.m_pdata);
        Matrix_Tns2_AxBxC(eF,1.0,0.0,F,hFI,pFI);
        Matrix_cleanup(pFI);
        Matrix_cleanup(hFI);
      }
      else
        err += func->get_eFn(m, &eF);

      err += construct_model_context(&ctx, m->param->type, F.m_pdata,dt,alpha,eF.m_pdata);
          
      if((m->param)->uqcm)
      {       
        double *x_ip = (fe.x_ip).m_pdata;            
      
        ELASTICITY *elast = (m->param)->cm_elast;
        mat_e_in = elast->mat;
        err += material_properties_elasticity_at_ip(mat_e_in, &mat_e_new, x_ip[0], x_ip[1], x_ip[2]);
        elast->mat = &mat_e_new; // should be replaced by original mat_e_in after computation  
        err += (m->param)->update_elasticity(m,ctx,NULL,&S,0);
      
        elast->mat = mat_e_in;
      }
      else
        err += (m->param)->update_elasticity(m,ctx,NULL,&S,0);

      err += (m->param)->destroy_ctx(&ctx);
      // <-- update elasticity part

      /* get aliases to Matrix data for simpler access */
      const double *Sd = S.m_pdata;
      const double *eFd = eF.m_pdata;
      const double *pFd = pF.m_pdata;
      const double eJ = det3x3(eFd);      

      /* store symmetric part of S (PK2) */
      sig[i].il[ip].o[0] = Sd[idx_2(0,0)]; /* XX */
      sig[i].il[ip].o[1] = Sd[idx_2(1,1)]; /* YY */
      sig[i].il[ip].o[2] = Sd[idx_2(2,2)]; /* ZZ */
      sig[i].il[ip].o[3] = Sd[idx_2(1,2)]; /* YZ */
      sig[i].il[ip].o[4] = Sd[idx_2(0,2)]; /* XZ */
      sig[i].il[ip].o[5] = Sd[idx_2(0,1)]; /* XY */

      /* store elastic deformation */
      memcpy(eps[i].il[ip].F, eFd, DIM_3x3 * sizeof(*eFd));

      /* store the hardening parameter */
      err += func->get_hardening(m, &eps[i].dam[ip].wn);

      /* compute/store the plastic strain variable */
      err += func->get_plast_strain_var(m, &eps[i].dam[ip].Xn);

      /* Compute the Cauchy Stress sigma = 1/eJ eF S eF' */
      double sigma[DIM_3x3] = {};
      double temp[DIM_3x3] = {};
      double temp_I[DIM_3x3] = {};
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                  DIM_3,DIM_3,DIM_3, 1.0 / eJ, eFd, DIM_3, Sd, DIM_3,
                  0.0, temp,DIM_3);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                  DIM_3, DIM_3, DIM_3, 1.0, temp, DIM_3, eFd, DIM_3,
                  0.0, sigma, DIM_3);

      volume += fe.detJxW;
      
      /* store symmetric part */
      sig[i].el.o[0] += fe.detJxW*sigma[idx_2(0,0)]; /* XX */
      sig[i].el.o[1] += fe.detJxW*sigma[idx_2(1,1)]; /* YY */
      sig[i].el.o[2] += fe.detJxW*sigma[idx_2(2,2)]; /* ZZ */
      sig[i].el.o[3] += fe.detJxW*sigma[idx_2(1,2)]; /* YZ */
      sig[i].el.o[4] += fe.detJxW*sigma[idx_2(0,2)]; /* XZ */
      sig[i].el.o[5] += fe.detJxW*sigma[idx_2(0,1)]; /* XY */

      /* Compute the logarithmic strain e = 1/2(I - inv(FF'))*/
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                  DIM_3, DIM_3, DIM_3, 1.0, eFd, DIM_3, eFd, DIM_3,
                  0.0, temp, DIM_3);
      err += inv3x3(temp, temp_I);
      /* e <-- temp is the Euler strain */
      for(int ia = 0; ia < DIM_3x3; ia++)
        temp[ia] = 0.5*(eye[ia]-temp_I[ia]);

      /* store symmetric part (also Eng. strain) */
      eps[i].el.o[0] += fe.detJxW*temp[idx_2(0,0)];
      eps[i].el.o[1] += fe.detJxW*temp[idx_2(1,1)];
      eps[i].el.o[2] += fe.detJxW*temp[idx_2(2,2)];

      eps[i].el.o[3] += fe.detJxW*2.0*temp[idx_2(1,2)];
      eps[i].el.o[4] += fe.detJxW*2.0*temp[idx_2(0,2)];
      eps[i].el.o[5] += fe.detJxW*2.0*temp[idx_2(0,1)];

    }
    for(int ia=0; ia<6; ia++)
    {
      sig[i].el.o[ia] = sig[i].el.o[ia]/volume;
      eps[i].el.o[ia] = eps[i].el.o[ia]/volume;
    }

    if(is_it_couple_w_thermal >=0)
    {
      Matrix_cleanup(Tnm1);    
      Matrix_cleanup(Tnp1);
      Matrix_cleanup(Tn);
    } 

    FEMLIB_destruct(&fe);
  }

  Matrix_cleanup(F);
  Matrix_cleanup(eF);
  Matrix_cleanup(pF);
  Matrix_cleanup(S);

  return err;
}

/// compute element residual for mid point rule
///
/// compute residual(n+alpha)
///
/// \param[out] f computed element residual
/// \param[in] m constitutive model object
/// \param[in] ii element id
/// \param[in] ndofn number of degree freedom on node
/// \param[in] pFnp1 2nd order tenosr pF(n+1)
/// \param[in] pFn   2nd order tenosr pF(n)                                     
/// \param[in] Fnp1  2nd order tenosr F(n+1)
/// \param[in] Fn    2nd order tenosr F(n)
/// \param[in] hFnp1 2nd order tenosr hF(n+1)
/// \param[in] hFn   2nd order tenosr hF(n)
/// \param[in] is_it_couple_w_chemical flag for coupling with thermal
/// \param[in] alpha mid point rule alpha
/// \param[in] dt_alpha_1_minus_alpha -dt*(1-alpha) for mid point btw t(n) t(n+1) 
///                                   -dt*alpha     for mid point btw t(n-1) t(n)
/// \param[in] fe finite element helper object
/// \return non-zero on internal error
int residuals_el_constitutive_model_n_plus_alpha(double *f,
                                                 const Constitutive_model *m,
                                                 const int ii,
                                                 const int ndofn,
                                                 const Matrix(double) *pFnp1,
                                                 const Matrix(double) *pFn,                                                
                                                 const Matrix(double) *Fnp1,
                                                 const Matrix(double) *Fn,
                                                 const Matrix(double) *hFnp1,
                                                 const Matrix(double) *hFn,
                                                 const int is_it_couple_w_thermal,                                                 
                                                 const double alpha,
                                                 const double dt_alpha_1_minus_alpha,
                                                 FEMLIB *fe)
{
  // Total Lagrangian based
  int err = 0;
  const int nne = fe->nne;
  const int nsd = fe->nsd;  
  
  enum {M,eFnpa,pFnpa,pFnpa_I,hFnpa,Fnpa,S,MT,Fend}; 
  
  // list of second-order tensors
  Matrix(double) *F2;
  err += construct_matrix_array(&F2,DIM_3,DIM_3,Fend,0); 
  
  for(int ia=0; ia<DIM_3x3; ia++)
  {
        F2[M].m_pdata[ia] = 0.0;
    F2[eFnpa].m_pdata[ia] = 0.0;
  }   

  int compute_stiffness = 0;
   
  mid_point_rule(F2[pFnpa].m_pdata, pFn->m_pdata, pFnp1->m_pdata, alpha, nsd*nsd);
  Matrix_AeqB(F2[hFnpa], 1.0,*hFnp1);
//  mid_point_rule(F2[hFnpa].m_pdata, hFn->m_pdata, hFnp1->m_pdata, alpha, nsd*nsd);
  mid_point_rule( F2[Fnpa].m_pdata,  Fn->m_pdata,  Fnp1->m_pdata, alpha, nsd*nsd);

  if(is_it_couple_w_thermal>0)
  {
    Matrix(double) hFnpa_I;
    Matrix_construct_redim(double,hFnpa_I,DIM_3,DIM_3);
        
    err += inv3x3(F2[hFnpa].m_pdata, hFnpa_I.m_pdata);
    err += inv3x3(F2[pFnpa].m_pdata, F2[pFnpa_I].m_pdata);
    Matrix_AxB(F2[M],1.0,0.0,hFnpa_I,0,F2[pFnpa_I],0);
    Matrix_cleanup(hFnpa_I);    
  }
  else
    err += inv3x3(F2[pFnpa].m_pdata, F2[M].m_pdata);
   
  Matrix_AxB(F2[eFnpa], 1.0,0.0,F2[Fnpa],0,F2[M],0);
        
  Matrix_AeqBT(F2[MT],1.0,F2[M]);  
  Matrix_init(F2[S],0.0);    

  {
    // check that deformation is invertible -> J > 0
    int terr = 0;
    double tJ = 0.0;
    tJ = getJacobian(F2[Fnpa].m_pdata, ii, &terr);
    err += terr;
  }
  
  void *ctx;

  err += construct_model_context(&ctx, m->param->type, Fnp1->m_pdata,0.0,alpha,F2[eFnpa].m_pdata);
  err += (m->param)->update_elasticity(m,ctx,NULL,&F2[S],compute_stiffness);
  err += m->param->destroy_ctx(&ctx);

  if(err==0)
  { 
    double Jn = 1.0; 
    err += compute_residual_vector(f,fe,F2+Fnpa,F2+MT,F2+M,F2+S,Jn/dt_alpha_1_minus_alpha);
  }          
  
  // destroy second-order tensors
  err += cleanup_matrix_array(&F2, Fend);
  return err;
}

int cm_get_subdivision_parameter(double *subdiv_param,
                                 const int ne,
                                 const ELEMENT *elem,
                                 const EPS *eps,
                                 const double dt)
{
  int err = 0;
  *subdiv_param = 0.0;
  double cur_val = 0.0;
  double max_val = 0.0;
  for (int i = 0; i < ne; i++) {
    long n_ip = 0;
    int_point(elem[i].toe, &n_ip);
    for (int ip = 0; ip < n_ip; ip++) {
      err += eps[i].model[ip].param->get_subdiv_param(&(eps[i].model[ip]), &cur_val, dt);
      max_val = MAX(max_val, cur_val);
    }
  }

  *subdiv_param = max_val;
  return err;
}

/// compute element stiffness matrix in transient
///
/// Compute element stiffness matrix based on the mid point rule. When thermal
/// is couled, temperature is assumed constant. Currenlty total Lagrangian is 
/// active (Updated lagrangian is implemented, but accelleration term is not 
/// fully implemented for updated lagrangian)
///
/// \param[in] fe finite element helper object
/// \param[out] lk computed element stiffness matrix
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int stiffness_el_constitutive_model_w_inertia(FEMLIB *fe,
                                              double *lk,
                                              double *r_e,
                                              GRID *grid,
                                              MATERIAL_PROPERTY *mat,
                                              FIELD_VARIABLES *fv,
                                              SOLVER_OPTIONS *sol,
                                              LOADING_STEPS *load,
                                              CRPL *crpl,
                                              const PGFem3D_opt *opts,
                                              MULTIPHYSICS *mp,
                                              int mp_id,
                                              double dt)
{
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int err = 0;
  double alpha = sol->alpha;
  int total_Lagrangian = 1;  
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  
  for(int ia=0; ia<fv->n_coupled; ia++)
  { 
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;    
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }      

// when updated Lagrangian is used should be un-commented  
//  if(opts->cm != 0)
//    total_Lagrangian = 1;

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  SUPP sup = load->sups[mp_id];

  double *u = malloc(sizeof(*u)*nne*nsd);
  double *dMdu_all = malloc(sizeof(*dMdu_all)*DIM_3x3*nne*nsd);
  memset(dMdu_all,0,DIM_3x3*nne*nsd*sizeof(double));

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }
  
  enum {Fn,Fr,Fnp1,pFn,pFnp1,S,
        eFnpa,pFnpa,pFnpa_I,
        eFn,M,eFnM,eFnMT,FrTFr,
        hFnp1,hFn,hFnpa,hFnpa_I,Fend};
  
  Matrix(double) L;  
  Matrix_construct_redim(double,L ,DIM_3x3x3x3,1);  

  Matrix(double) Tnp1, Tn, Tnm1; 
  FIELD_VARIABLES *fv_h = NULL;
  
  if(is_it_couple_w_thermal >= 0)
  {
    fv_h = fv->fvs[is_it_couple_w_thermal]; 
    Matrix_construct_init(double, Tnm1,fe->nne,1,0.0);
    Matrix_construct_init(double, Tnp1,fe->nne,1,0.0); 
    Matrix_construct_init(double, Tn,  fe->nne,1,0.0);
    
    // compute temperature for this element for each nodal point

    int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
    err += get_nodal_temperatures(fe, grid, fv_h, load, mp_cp_id,
                                  Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                  fv->subdivision_factor_np1,fv->subdivision_factor_n);
  }

  // list of second-order tensors
  Matrix(double) *F2;
  err + construct_matrix_array(&F2,DIM_3,DIM_3,Fend,0);
  
  for(int ia=0; ia<DIM_3x3; ia++)
  {
        F2[M].m_pdata[ia] = 0.0;
    F2[eFnpa].m_pdata[ia] = 0.0;
    F2[FrTFr].m_pdata[ia] = 0.0;
        F2[S].m_pdata[ia] = 0.0;
  } 

  int compute_stiffness = 1;      

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    FEMLIB_elem_basis_V(fe, ip);
    FEMLIB_update_shape_tensor(fe);  
    FEMLIB_update_deformation_gradient(fe,ndofn,u,F2+Fnp1);
    
    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);
    
    // get a shortened pointer for simplified CM function calls 
    const Model_parameters *func = m->param;
    
    Matrix_eye(F2[hFnpa_I],DIM_3);
    if(is_it_couple_w_thermal >= 0)
    {
      double T0 = fv_h->u0;
      double hFnm1[9];
      err += compute_temperature_at_ip(fe,grid,mat,T0,
                                       Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                       F2[hFnp1].m_pdata,F2[hFn].m_pdata,hFnm1);
      
      Matrix_AeqB(F2[hFnpa], 1.0, F2[hFnp1]);                                 
      //mid_point_rule(F2[hFnpa].m_pdata, F2[hFn].m_pdata, F2[hFnp1].m_pdata, alpha, DIM_3x3);                                  
      Matrix_inv(F2[hFnpa],F2[hFnpa_I]);                                        
    }
        
    err += func->get_pF(m,&F2[pFnp1]);
    err += func->get_pFn(m,&F2[pFn]);
    err += func->get_Fn(m,&F2[Fn]);

               
    mid_point_rule(F2[pFnpa].m_pdata, F2[pFn].m_pdata, F2[pFnp1].m_pdata, alpha, DIM_3x3);
    mid_point_rule(   F2[Fr].m_pdata,  F2[Fn].m_pdata,  F2[Fnp1].m_pdata, alpha, DIM_3x3);

    err += inv3x3(F2[pFnpa].m_pdata, F2[pFnpa_I].m_pdata);
    
    Matrix_AxB(F2[M],1.0,0.0,F2[hFnpa_I],0,F2[pFnpa_I],0);
    Matrix_AxB(F2[eFnpa], 1.0,0.0,F2[Fr],0,F2[M],0);
    Matrix_AxB(F2[FrTFr],1.0,0.0,F2[Fr],1,F2[Fr],0);  

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, F2[Fnp1].m_pdata,dt,alpha, F2[eFnpa].m_pdata,
                                                  F2[hFn].m_pdata,F2[hFnp1].m_pdata);      
    else
      err += construct_model_context(&ctx, m->param->type, F2[Fnp1].m_pdata,dt,alpha, F2[eFnpa].m_pdata);
        
    err += m->param->compute_dMdu(m, ctx, fe->ST, nne, ndofn, dMdu_all);
    
    // --> update elasticity part
    Matrix_init(L,0.0);
    Matrix_init(F2[S],0.0);    

    err += (m->param)->update_elasticity(m,ctx,&L,&F2[S],compute_stiffness);
    // <-- update elasticity part
    err += m->param->destroy_ctx(&ctx);

    if(err!=0)
      break;
      
    // --> start computing tagent
    //Matrix_AxB(F2[eFnM],1.0,0.0,F2[eFn],0,F2[M],0);
    Matrix_eye(F2[eFn],DIM_3);
    Matrix_AeqB(F2[eFnM],1.0,F2[M]); // eFn = 1 for total_Lagrangian
    Matrix_AeqBT(F2[eFnMT],1.0,F2[eFnM]);
        
    double Jn; // Matrix_det(F2[Fn], Jn);
    Jn = 1.0;

    err += compute_stiffness_matrix(lk,fe,
                                    F2+Fr,F2+eFnMT,F2+eFn,F2+M,F2+FrTFr,F2+eFnM,F2+S,
                                    &L,dMdu_all,Jn);
  }
  
  // free memory for thermal part
  if(is_it_couple_w_thermal>=0)
  {
    Matrix_cleanup(Tnm1);    
    Matrix_cleanup(Tnp1);
    Matrix_cleanup(Tn);
  } 
    
  free(u);
  
  Matrix_cleanup(L);

  // destroy second-order tensors
  err += cleanup_matrix_array(&F2, Fend);

  free(dMdu_all);
 
  return err;
}

/// compute element stiffness matrix in quasi steady state
///
/// Updated Lagrangian and total Lagrangian based. When thermal
/// is couled, temperature is assumed constant.
///
/// \param[in] fe finite element helper object
/// \param[out] lk computed element stiffness matrix
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int stiffness_el_constitutive_model(FEMLIB *fe,
                                    double *lk,
                                    double *r_e,
                                    GRID *grid,
                                    MATERIAL_PROPERTY *mat,
                                    FIELD_VARIABLES *fv,
                                    SOLVER_OPTIONS *sol,
                                    LOADING_STEPS *load,
                                    CRPL *crpl,
                                    const PGFem3D_opt *opts,
                                    MULTIPHYSICS *mp,
                                    int mp_id,
                                    double dt)
{
  int err = 0;
  double alpha = -1.0; // if alpha < 0, no inertia
  int total_Lagrangian = 0;  
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  
  for(int ia=0; ia<fv->n_coupled; ia++)
  { 
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;    
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }      
  
  if(opts->cm != 0)
    total_Lagrangian = 1;

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  SUPP sup = load->sups[mp_id];
  
  double *u = malloc(sizeof(*u)*nne*nsd);
  double *dMdu_all = malloc(sizeof(*dMdu_all)*DIM_3x3*nne*nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }

  enum {Fr,Fnp1,pFnp1,S,
        eFn,M,eFnM,eFnMT,FrTFr,
        hFn,hFnp1,hFnp1_I,temp_F_1,temp_F_2,Fend};
          
  Matrix(double) L;  
  Matrix_construct_redim(double,L ,DIM_3x3x3x3,1);
    
  Matrix(double) Tnp1, Tn, Tnm1;    
  FIELD_VARIABLES *fv_h = NULL;
  
  if(is_it_couple_w_thermal >= 0)
  {
    fv_h = fv->fvs[is_it_couple_w_thermal]; 
    Matrix_construct_init(double, Tnm1,fe->nne,1,0.0);
    Matrix_construct_init(double, Tnp1,fe->nne,1,0.0); 
    Matrix_construct_init(double, Tn,  fe->nne,1,0.0);
    
    // compute temperature for this element for each nodal point
    int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
    err += get_nodal_temperatures(fe, grid, fv_h, load, mp_cp_id, 
                                  Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                  fv->subdivision_factor_np1,fv->subdivision_factor_n);
  }  

  // list of second-order tensors
  Matrix(double) *F2;
  err + construct_matrix_array(&F2,DIM_3,DIM_3,Fend,0);

  for(int ia=0; ia<DIM_3x3; ia++)
  {
     F2[Fnp1].m_pdata[ia] = 0.0;
        F2[M].m_pdata[ia] = 0.0;
      F2[eFn].m_pdata[ia] = 0.0;   
     F2[eFnM].m_pdata[ia] = 0.0;
    F2[FrTFr].m_pdata[ia] = 0.0;
  }   

  int compute_stiffness = 1;      

  MATERIAL_ELASTICITY mat_e_new;
  MATERIAL_ELASTICITY *mat_e_in;

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double Jn    = 1.0;
        
    FEMLIB_elem_basis_V(fe, ip);
    FEMLIB_update_shape_tensor(fe);  
    FEMLIB_update_deformation_gradient(fe,ndofn,u,F2+Fr);
    
    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);
    
    // get a shortened pointer for simplified CM function calls
    const Model_parameters *func = m->param;

    // --> update deformations due to coupled physics
    if(is_it_couple_w_thermal >= 0)
    { 
      // compute temperature at the integration point
      double T0 = fv_h->u0;

      int hFnm1 = temp_F_1;
      err += compute_temperature_at_ip(fe,grid,mat,T0,
                                       Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                       F2[hFnp1].m_pdata,F2[hFn].m_pdata,F2[hFnm1].m_pdata);
      Matrix_inv(F2[hFnp1],F2[hFnp1_I]);
    }
    
    // --> update plasticity part
    if(total_Lagrangian)
    {
      if(sup->multi_scale)
        cm_add_macro_F(sup,F2[Fr].m_pdata);

      // Total Lagrangian formulation Fn = 1, Fnp1 = Fr      
      Matrix_eye(F2[eFn],DIM_3);
      Matrix_copy(F2[Fnp1], F2[Fr]);
      Matrix_eye(F2[hFn], DIM_3);
    }
    else
    {
      if(sup->multi_scale)
      {
        PGFEM_printerr("Multi-scale formulation does not support UL!\n");
        PGFEM_Abort();
      }
 
      int Fn = temp_F_1;
      err += m->param->get_Fn(m, F2+Fn);
      Matrix_AxB(F2[Fnp1],1.0,0.0,F2[Fr],0,F2[Fn],0); // F2[Fn+1] = Fr*Fn
      Matrix_det(F2[Fn], Jn);
    }
        
    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, F2[Fnp1].m_pdata,dt,alpha, NULL,
                                                  F2[hFn].m_pdata,F2[hFnp1].m_pdata);
    else
      err += construct_model_context(&ctx, m->param->type, F2[Fnp1].m_pdata,dt,alpha, NULL);
  
    err += func->compute_dMdu(m, ctx, fe->ST, nne, ndofn, dMdu_all);
    err += func->get_pF(m,F2+pFnp1);
        
    if(total_Lagrangian) // Total Lagrangian formulation, all xFn = 1
    { 
      if(is_it_couple_w_thermal >= 0)
      { 
        int pFnp1_I = temp_F_1; 
        Matrix_inv(F2[pFnp1], F2[pFnp1_I]);        
        Matrix_AxB(F2[M],1.0,0.0,F2[hFnp1_I],0,F2[pFnp1_I],0);        
      }
      else
        Matrix_inv(F2[pFnp1], F2[M]);

      Matrix_AeqB(F2[eFnM],1.0,F2[M]);
      Matrix_AeqBT(F2[eFnMT],1.0,F2[M]);                      
    }
    else
    {
      if(is_it_couple_w_thermal>=0)
      {
        int pFnp1_I = temp_F_1;
        Matrix_inv(F2[pFnp1], F2[pFnp1_I]);
        Matrix_Tns2_AxBxC(F2[temp_F_2],1.0,0.0,F2[hFn],F2[hFnp1_I],F2[pFnp1_I]);          
        int pFn     = temp_F_1;
        err += func->get_pFn(m,F2+pFn);        
        Matrix_AxB(F2[M],1.0,0.0,F2[pFn],0,F2[temp_F_2],0);
         
        int hFn_I = temp_F_2;
        int stepno = 1; // 0 = time step = n-1
                        // 1 = time step = n
                        // 2 = time step = n+1
        Matrix_inv(F2[hFn], F2[hFn_I]);                
        err += m->param->get_eF_of_hF(m,F2+eFn,F2+hFn_I,stepno); 
      }
      else
      {
        int pFnp1_I = temp_F_1;
        Matrix_inv(F2[pFnp1], F2[pFnp1_I]);
        int pFn     = temp_F_2;        
        err += func->get_pFn(m,F2+pFn);        
        Matrix_AxB(F2[M],1.0,0.0,F2[pFn],0,F2[pFnp1_I],0);
        err += m->param->get_eFn(m,F2+eFn);
      }
      Matrix_AxB(F2[eFnM],1.0,0.0,F2[eFn],0,F2[M],0);
      Matrix_AeqBT(F2[eFnMT],1.0,F2[eFnM]);   
    }   
    // <-- update plasticity part 

    Matrix_AxB(F2[FrTFr],1.0,0.0,F2[Fr],1,F2[Fr],0); 

    // --> update elasticity part
    Matrix_init(L,0.0);
    Matrix_init(F2[S],0.0);    

    if((m->param)->uqcm)
    { 
      double *x_ip = (fe->x_ip).m_pdata;      
      ELASTICITY *elast = (m->param)->cm_elast;
      mat_e_in = elast->mat;
      err += material_properties_elasticity_at_ip(mat_e_in, &mat_e_new, x_ip[0], x_ip[1], x_ip[2]);
      elast->mat = &mat_e_new; // should be replaced by original mat_e_in after computation
      err += (m->param)->update_elasticity(m,ctx,&L,&F2[S],compute_stiffness);
      elast->mat = mat_e_in;
    }
    else
      err += (m->param)->update_elasticity(m,ctx,&L,&F2[S],compute_stiffness);
    // <-- update elasticity part   
    
    err += func->destroy_ctx(&ctx);
    if(err!=0)
      break;    

    // start computing tagent
    err += compute_stiffness_matrix(lk,fe,
                                    F2+Fr,F2+eFnMT,F2+eFn,F2+M,F2+FrTFr,F2+eFnM,F2+S,
                                    &L,dMdu_all,Jn);
  }
  free(u);

  /* check diagonal for zeros/nans */
  for (int a = 0; a < nne; a++) {
    for (int b = 0; b < nsd; b++) {
      if ( !isnormal(lk[idx_K(a,b,a,b,nne,nsd)]) ) err++;
    }
  }
  
  Matrix_cleanup(L);

  // free memory for thermal part
  if(is_it_couple_w_thermal>=0)
  {
    Matrix_cleanup(Tnm1);    
    Matrix_cleanup(Tnp1);
    Matrix_cleanup(Tn);
  }

  // destroy second-order tensors
  err += cleanup_matrix_array(&F2, Fend);

  free(dMdu_all);
 
  return err;
}

/// compute element residual vector in transient
///
/// redual = residual(n+alpha) + residual(n-1+alpha)
/// If residual is computed during the iterative solution scheme (Newton iteration),
/// integration algorithm is performed. However, in the case of just checking residual,
/// no integration algorithm will be executed. The switch of running integration algorithm
/// is sol->run_integration_algorithm. 
///
/// \param[in] fe finite element helper object
/// \param[out] f computed element residual vector
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] dts time step size at t(n), t(n+1); dts[DT_N] = t(n) - t(n-1)
///                                                dts[DT_NP1] = t(n+1) - t(n)
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int residuals_el_constitutive_model_w_inertia(FEMLIB *fe,
                                              double *f,
                                              double *r_e,
                                              GRID *grid,
                                              MATERIAL_PROPERTY *mat,
                                              FIELD_VARIABLES *fv,
                                              SOLVER_OPTIONS *sol,
                                              LOADING_STEPS *load,
                                              CRPL *crpl,
                                              const PGFem3D_opt *opts,
                                              MULTIPHYSICS *mp,
                                              const double *dts,
                                              int mp_id,
                                              double dt)
{
  int err = 0;
  double alpha = sol->alpha;
  int total_Lagrangian = 1;  
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  
  for(int ia=0; ia<fv->n_coupled; ia++)
  { 
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;    
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

// when updated Lagrangian is used should be un-commented  
//  if(opts->cm != 0)
//    total_Lagrangian = 1;

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  SUPP sup = load->sups[mp_id];
      
  double *u       = (double *) malloc(sizeof(double)*nne*nsd);
  double *f_npa   = (double *) malloc(sizeof(double)*nne*nsd);
  double *f_nm1pa = (double *) malloc(sizeof(double)*nne*nsd);    
    
  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }
  
  enum {Fnp1,Fn,Fnm1,pFnp1,pFn,pFnm1,
        hFnp1,hFn,hFnm1,Fend}; 
  
  // list of second-order tensors
  Matrix(double) *F2;
  err + construct_matrix_array(&F2,DIM_3,DIM_3,Fend,0); 
  
  Matrix(double) Tnp1, Tn, Tnm1;    
  FIELD_VARIABLES *fv_h = NULL;
  
  if(is_it_couple_w_thermal >= 0)
  {
    fv_h = fv->fvs[is_it_couple_w_thermal]; 
    Matrix_construct_init(double, Tnm1,fe->nne,1,0.0);
    Matrix_construct_init(double, Tnp1,fe->nne,1,0.0); 
    Matrix_construct_init(double, Tn,  fe->nne,1,0.0);
    
    // compute temperature for this element for each nodal point
    int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
    err += get_nodal_temperatures(fe, grid, fv_h, load, mp_cp_id,
                                  Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                  fv->subdivision_factor_np1,fv->subdivision_factor_n);
  }  

  int compute_stiffness = 0;
  
  memset(f_npa, 0, sizeof(double)*nne*ndofn);   
  memset(f_nm1pa, 0, sizeof(double)*nne*ndofn);   
    
  for(int ip = 1; ip<=fe->nint; ip++)
  {
    FEMLIB_elem_basis_V(fe, ip);
    FEMLIB_update_shape_tensor(fe);  
    FEMLIB_update_deformation_gradient(fe,ndofn,u,F2+Fnp1);

    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);
    
    // get a shortened pointer for simplified CM function calls 
    const Model_parameters *func = m->param;
    
    if(is_it_couple_w_thermal >= 0)
    {
      double T0 = fv_h->u0;
      err += compute_temperature_at_ip(fe,grid,mat,T0,
                                       Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                       F2[hFnp1].m_pdata,F2[hFn].m_pdata,F2[hFnm1].m_pdata);                                      
    }
    
    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, F2[Fnp1].m_pdata,dts[DT_NP1],alpha, NULL,
                                                  F2[hFn].m_pdata,F2[hFnp1].m_pdata);      
    else
      err += construct_model_context(&ctx, m->param->type, F2[Fnp1].m_pdata,dts[DT_NP1],alpha, NULL);

    if(sol->run_integration_algorithm)
      err += m->param->integration_algorithm(m,ctx); // perform integration algorithm

    err += m->param->destroy_ctx(&ctx);

    if(err!=0)
      break;    
    
    err += m->param->get_pF(m,    F2+pFnp1);
    err += m->param->get_pFn(m,   F2+pFn);
    err += m->param->get_pFnm1(m, F2+pFnm1); 
    err += m->param->get_Fn(m,    F2+Fn);
    err += m->param->get_Fnm1(m,  F2+Fnm1);
    
    double dt_1_minus_alpha = -dts[DT_NP1]*(1.0-alpha);
    err += residuals_el_constitutive_model_n_plus_alpha(f_npa,m,eid,ndofn,
                                                        F2+pFnp1,F2+pFn,F2+Fnp1,F2+Fn,
                                                        F2+hFnp1,F2+hFn,
                                                        is_it_couple_w_thermal,
                                                        alpha, dt_1_minus_alpha,fe);
                                
    double dt_alpha = -dts[DT_N]*alpha;

    err += residuals_el_constitutive_model_n_plus_alpha(f_nm1pa,m,eid,ndofn,
                                                        F2+pFn,F2+pFnm1,F2+Fn,F2+Fnm1,
                                                        F2+hFnp1,F2+hFn,
                                                        is_it_couple_w_thermal,
                                                        alpha, dt_alpha,fe);
  }
  
  if(err==0)
  {
    for(int a=0; a<nne*nsd; a++)
      f[a] += f_npa[a] + f_nm1pa[a];
  }
  
  // free memory for thermal part
  if(is_it_couple_w_thermal>=0)
  {
    Matrix_cleanup(Tnm1);    
    Matrix_cleanup(Tnp1);
    Matrix_cleanup(Tn);
  }   
  
  free(u);
  free(f_npa);
  free(f_nm1pa);
  
  // destroy second-order tensors
  err += cleanup_matrix_array(&F2, Fend);

  return err;
}

/// compute element residual vector in quasi steady state
///
/// If residual is computed during the iterative solution scheme (Newton iteration),
/// integration algorithm is performed. However, in the case of just checking residual,
/// no integration algorithm will be executed. The switch of running integration algorithm
/// is sol->run_integration_algorithm. 
///
/// \param[in] fe finite element helper object
/// \param[out] f computed element residual vector
/// \param[in] r_e nodal variabls(displacements) on the current element
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int residuals_el_constitutive_model(FEMLIB *fe,
                                    double *f,
                                    double *r_e,
                                    GRID *grid,
                                    MATERIAL_PROPERTY *mat,
                                    FIELD_VARIABLES *fv,
                                    SOLVER_OPTIONS *sol,
                                    LOADING_STEPS *load,
                                    CRPL *crpl,
                                    const PGFem3D_opt *opts,
                                    MULTIPHYSICS *mp,
                                    int mp_id,
                                    double dt)
{
  int err = 0;
  double alpha = -1.0; // if alpha < 0, no inertia
  int total_Lagrangian = 0;
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  
  for(int ia=0; ia<fv->n_coupled; ia++)
  { 
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;    
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }      

  if(opts->cm != 0)
    total_Lagrangian = 1;
      
  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  SUPP sup = load->sups[mp_id];
      
  double *u = (double *) malloc(sizeof(double)*nne*nsd);
  
  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];  
  }
  
  enum {Fr,Fnp1,pFnp1,
        L,S,M,eFnM,eFnMT,
        hFn,hFnp1,hFnp1_I,
        temp_F_1,
        temp_F_2,Fend};


  Matrix(double) Tnp1, Tn, Tnm1;    
  FIELD_VARIABLES *fv_h = NULL;
  
  if(is_it_couple_w_thermal >= 0)
  {
    fv_h = fv->fvs[is_it_couple_w_thermal]; 
    Matrix_construct_init(double, Tnm1,fe->nne,1,0.0);    
    Matrix_construct_init(double, Tnp1,fe->nne,1,0.0); 
    Matrix_construct_init(double, Tn,  fe->nne,1,0.0);

    // compute temperature for this element for each nodal point
    int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
    err += get_nodal_temperatures(fe, grid, fv_h, load, mp_cp_id, 
                                  Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                  fv->subdivision_factor_np1,fv->subdivision_factor_n);
  }
  
  // list of second-order tensors
  Matrix(double) *F2;
  err + construct_matrix_array(&F2,DIM_3,DIM_3,Fend,0);
  
  for(int ia=0; ia<DIM_3x3; ia++)
  {
     F2[Fnp1].m_pdata[ia] = 0.0;
        F2[M].m_pdata[ia] = 0.0;
     F2[eFnM].m_pdata[ia] = 0.0;
     F2[temp_F_1].m_pdata[ia] = 0.0;
     F2[temp_F_2].m_pdata[ia] = 0.0;     
  }    

  int compute_stiffness = 0;

  MATERIAL_ELASTICITY mat_e_new;
  MATERIAL_ELASTICITY *mat_e_in;
     
  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double Jn = 1.0; // if upated Lagrangian, Jn = det(Fn), later updated

    FEMLIB_elem_basis_V(fe, ip);
    FEMLIB_update_shape_tensor(fe);  
    FEMLIB_update_deformation_gradient(fe,ndofn,u,F2+Fr);

    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);
    
    // get a shortened pointer for simplified CM function calls
    const Model_parameters *func = m->param;
    
    // --> update deformations due to coupled physics
    if(is_it_couple_w_thermal >= 0)
    { 
      // compute temperature at the integration point
      double T0 = fv_h->u0;
      double hFnm1[9];      
      err += compute_temperature_at_ip(fe,grid,mat,T0,
                                       Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                       F2[hFnp1].m_pdata,F2[hFn].m_pdata,hFnm1);
      Matrix_inv(F2[hFnp1], F2[hFnp1_I]);                                       
      if(total_Lagrangian)
        Matrix_eye(F2[hFn], DIM_3);
    }        

    // --> update plasticity part
    if(total_Lagrangian)
    {
      if (sup->multi_scale) {
        cm_add_macro_F(sup,F2[Fr].m_pdata);
      }

      // TOTAL LAGRANGIAN FORMULATION Fn = 1, Fnp1 = Fr
      Matrix_copy(F2[Fnp1], F2[Fr]);      
    }
    else
    {
      if (sup->multi_scale) {
        PGFEM_printerr("Multi-scale formulation does not support UL!\n");
        PGFEM_Abort();
      }
 
      int Fn = temp_F_1;
      err += m->param->get_Fn(m, F2+Fn);      
      Matrix_AxB(F2[Fnp1],1.0,0.0,F2[Fr],0,F2[Fn],0);  // compute F2[Fnp1] = Fr*Fn
      Matrix_det(F2[Fn], Jn);   
    }      

    {
      /* check that deformation is invertible -> J > 0 */
      int terr = 0;
      double tJ = 0.0;
      tJ = getJacobian(F2[Fnp1].m_pdata, eid, &terr);
      err += terr;
    }
    
    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, F2[Fnp1].m_pdata,dt,alpha, NULL,
                                                  F2[hFn].m_pdata,F2[hFnp1].m_pdata);
    else
      err += construct_model_context(&ctx, m->param->type, F2[Fnp1].m_pdata,dt,alpha, NULL);      
    
    if(sol->run_integration_algorithm)
      err += m->param->integration_algorithm(m,ctx); // perform integration algorithm

    if(err>0)
    	return err;

    err += func->get_pF(m,F2+pFnp1);

    if(total_Lagrangian)
    {
      if(is_it_couple_w_thermal>=0)
      {
        int pFnp1_I = temp_F_1;
        Matrix_inv(F2[pFnp1], F2[pFnp1_I]);        
        Matrix_AxB(F2[M],1.0,0.0,F2[hFnp1_I],0,F2[pFnp1_I],0);
      }
      else
        Matrix_inv(F2[pFnp1], F2[M]);            

      Matrix_AeqB(F2[eFnM],1.0,F2[M]);
      Matrix_AeqBT(F2[eFnMT],1.0,F2[M]);      

    }
    else
    {
      int eFn = temp_F_1;
      if(is_it_couple_w_thermal>=0)
      {
        int pFnp1_I = temp_F_1;
        Matrix_inv(F2[pFnp1], F2[pFnp1_I]);
        Matrix_Tns2_AxBxC(F2[temp_F_2],1.0,0.0,F2[hFn],F2[hFnp1_I],F2[pFnp1_I]);          
        int pFn     = temp_F_1;
        err += func->get_pFn(m,F2+pFn);        
        Matrix_AxB(F2[M],1.0,0.0,F2[pFn],0,F2[temp_F_2],0);
         
        int hFn_I = temp_F_2;
        int stepno = 1; // 0 = time step = n-1
                        // 1 = time step = n
                        // 2 = time step = n+1
        Matrix_inv(F2[hFn], F2[hFn_I]);                
        err += m->param->get_eF_of_hF(m,F2+eFn,F2+hFn_I,stepno);        
        
      }
      else
      {
        int pFnp1_I = temp_F_1;
        Matrix_inv(F2[pFnp1], F2[pFnp1_I]);
        int pFn     = temp_F_2;        
        err += func->get_pFn(m,F2+pFn);        
        Matrix_AxB(F2[M],1.0,0.0,F2[pFn],0,F2[pFnp1_I],0);
        err += m->param->get_eFn(m,F2+eFn);
      }
      Matrix_AxB(F2[eFnM],1.0,0.0,F2[eFn],0,F2[M],0);
      Matrix_AeqBT(F2[eFnMT],1.0,F2[eFnM]);      
    }
      
    // <-- update plasticity part
      
    // --> update elasticity part
    Matrix_init(F2[S],0.0);    
      
    if((m->param)->uqcm)
    { 
      double *x_ip = (fe->x_ip).m_pdata;      
      ELASTICITY *elast = (m->param)->cm_elast;
      mat_e_in = elast->mat;
      err += material_properties_elasticity_at_ip(mat_e_in, &mat_e_new, x_ip[0], x_ip[1], x_ip[2]);
      elast->mat = &mat_e_new; // should be replaced by original mat_e_in after computation
      err += (m->param)->update_elasticity(m,ctx,NULL,F2+S,compute_stiffness);
      elast->mat = mat_e_in;
    }
    else
      err += (m->param)->update_elasticity(m,ctx,NULL,F2+S,compute_stiffness);
    // <-- update elasticity part
    err += m->param->destroy_ctx(&ctx);

    if(err!=0)
      break;     

    err += compute_residual_vector(f,fe,F2+Fr,F2+eFnMT,F2+eFnM,F2+S,Jn);
  }       
  
  free(u);
  
  // free memory for thermal part
  if(is_it_couple_w_thermal>=0)
  {
    Matrix_cleanup(Tnm1);    
    Matrix_cleanup(Tnp1);
    Matrix_cleanup(Tn);
  }  
  
  // destroy second-order tensors
  err += cleanup_matrix_array(&F2, Fend);
  return err;
}                                   
