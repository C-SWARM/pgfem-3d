/**
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  Aaron Howell, [1], <ahowell3@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "PGFEM_io.h"
#include "PGFEM_mpi.h"
#include "PGFem3D_data_structure.h"
#include "allocation.h"
#include "constitutive_model.h"
#include "constitutive_model_3f.h"
#include "cm_placeholder_functions.h"
#include "cm_iso_viscous_damage.h"
#include "cm_j2_plasticity.h"
#include "cm_uqcm.h"
#include "cm_poro_viscoplasticity.h"
#include "dynamics.h"
#include "hommat.h"
#include "supp.h"
#include "elem3d.h"
#include "femlib.h"
#include "get_dof_ids_on_elem.h"
#include "hommat.h"
#include "hyperelasticity.h"                 // <= constitutive model elasticity
#include "index_macros.h"
#include "material_properties.h"    // <= constitutive model material properties
#include "plasticity_model_none.h"
#include "plasticity_model.h"
#include "plasticity_model_BPA.h"
#include "supp.h"
#include "utils.h"
#include <ttl/ttl.h>
#include <algorithm>
#include "condense.h"
#include "new_potentials.h"

using pgfem3d::Solver;

int
Model_var_info::print_variable_info(FILE *f)
{
  int err = 0;
  fprintf(f,"F names: ");
  for(int i = 0, e = this->n_Fs; i < e; i++) fprintf(f,"%s ",this->F_names[i]);
  fprintf(f,"\nVar names: ");
  for(int i = 0, e = this->n_vars; i < e; i++) fprintf(f,"%s ",this->var_names[i]);
  fprintf(f,"\nFlag names: ");
  for(int i = 0, e = this->n_flags; i < e; i++) fprintf(f,"%s ",this->flag_names[i]);
  fprintf(f,"\n");
  return err;
}

Model_var_info::~Model_var_info()
{
  /// deallocate internal memory
  for (size_t ia=0; ia<this->n_Fs; ia++)
  {
    if(this->F_names[ia])
      delete this->F_names[ia];
  }
  if(this->F_names)
    delete this->F_names;

  for(size_t ia=0; ia<this->n_vars; ia++)
  {
    if(this->var_names[ia])
      delete this->var_names[ia];
  }
  if(this->var_names)
    delete this->var_names;

  for(size_t ia=0; ia<this->n_flags; ia++)
  {
    if(this->flag_names[ia])
      delete this->flag_names[ia];
  }
  if(this->flag_names)
    delete this->flag_names;
}

//ttl declarations
namespace {
using namespace ttl;

constexpr Index<'i'> i;
constexpr Index<'j'> j;
constexpr Index<'k'> k;
constexpr Index<'l'> l;
constexpr Index<'o'> o;
constexpr Index<'p'> p;

constexpr int       DIM_3 = 3;
constexpr int     DIM_3x3 = 9;
constexpr int   DIM_3x3x3 = 27;
constexpr int DIM_3x3x3x3 = 81;

template<class T1, class T2>
inline int inv(T1 &A, T2 &AI) {
  int err = inv3x3(A.data, AI.data);
  return err;
}
}

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
                           Grid *grid,
                           FieldVariables *fv_h,
                           LoadingSteps *load,
                           int mp_id,
                           double *Tnp1_e,
                           double *Tn_e,
                           double *Tnm1_e,
                           double factor_np1,
                           double factor_n)
{
  int err = 0;

  int nne = fe->nne;
  long *cn = (long *) aloc1l(nne);
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
    {
      T_np1 = T0 + sup->defl[aid] + sup->defl_d[aid];
    }

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
                              const Grid *grid,
                              const MaterialProperty *mat,
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
    dTnm1 += N[ia]*(Tnm1[ia] - T0);    
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
                           double * __restrict F)
{
  const double * __restrict F0 = sup->F0;
  for (int i = 0; i < DIM_3x3; i++) F[i] += F0[i];
}

/// this is a wrapper function for the switch that was copy/pasted
/// everywhere. It is no big deal to keep adding to this private
/// function's argument list. Just put everything any model might need
/// and the switch will handle what is actually used.
int construct_model_context(void **ctx,
                            const int type,
                            double *F,
                            const double dt,
                            const double alpha,
                            double *eFnpa,
                            int npa)
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
  case POROVISCO_PLASTICITY:
    err += poro_viscoplasticity_model_ctx_build(ctx, F, dt,alpha, eFnpa, NULL, NULL, 0,npa);
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
                                         double *hFnp1,
                                         int npa)
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
  case POROVISCO_PLASTICITY:
    err += poro_viscoplasticity_model_ctx_build(ctx, F, dt,alpha, eFnpa, hFn, hFnp1, 1,npa);
    break;
  default:
    PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n", type);
    err++;
    break;
  }
  assert (err == 0);
  return err;
}

int
Constitutive_model::initialization(const Model_parameters *p)
{
  int err = 0;
  if (p == NULL)
    err++;
  else
  {
    this->param = p;
    Model_var_info info;
    this->param->get_var_info(info);

    err += this->vars_list[0][this->model_id].initialization(info.n_Fs,
                                                             info.n_vars,
                                                             info.n_flags);
    this->param->set_init_vals(this);
  }
  return err;
}

/// Construct a Model_parameters. Model type will be determined at runtime.
///
/// \param[out] **p        a Model_parameters, p[model_id] will be constructed
/// \param[in]  model_id   index of p
/// \param[in]  model_type constitutive model type
/// \return non-zero on error.
int construct_Model_parameters(Model_parameters **p, int model_id, int model_type)
{
  int err = 0;

  if(p[model_id] != NULL)
  {  
    delete p[model_id];
    p[model_id] = NULL;
  }
  
  switch(model_type)
  {
    case TESTING:
    case HYPER_ELASTICITY:
      p[model_id] = new HE_PARAM;
      break;
    case CRYSTAL_PLASTICITY:
      p[model_id] = new CP_PARAM;
      break;
    case BPA_PLASTICITY:
      p[model_id] = new BPA_PARAM;
      break;
    case ISO_VISCOUS_DAMAGE:
      p[model_id] = new CM_IVD_PARAM;
      break;
    case J2_PLASTICITY_DAMAGE:
      p[model_id] = new CM_J2P_PARAM;
      break;
    case POROVISCO_PLASTICITY:
      p[model_id] = new CM_PVP_PARAM;
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",model_type);
      err++;
      return err;
  }

  if(p[model_id] != NULL)
    p[model_id]->set_nulls();
    
  return err;
}
/// Initialize the Model_parameters object. The object may be used
/// after calling this function. Calling this function on an already
/// initialized object is undefined.
///
/// \param[in] p_hmat material property object
/// \param[in] type   constitutive model type
/// \return non-zero on error.
int
Model_parameters::initialization(const HOMMAT *p_hmat, const size_t type)
{
  int err = 0;

  switch(type)
  {
    case TESTING:
    case HYPER_ELASTICITY:
    case CRYSTAL_PLASTICITY:
    case BPA_PLASTICITY:
    case ISO_VISCOUS_DAMAGE:
    case J2_PLASTICITY_DAMAGE:
    case POROVISCO_PLASTICITY:
      break; // no action
    default:
      PGFEM_printerr("ERROR: Unrecognized model type! (%zd)\n",type);
      err++;
      return err;
  }

  this->p_hmat = p_hmat;
  this->type = type;

  auto* cm_mat = PGFEM_malloc<MATERIAL_CONSTITUTIVE_MODEL>();
  auto* mat_e = PGFEM_malloc<MATERIAL_ELASTICITY>();
  auto* elast = PGFEM_malloc<ELASTICITY>();

  set_properties_using_E_and_nu(mat_e,p_hmat->E,p_hmat->nu);
  mat_e->m01 = p_hmat->m01;
  mat_e->m10 = p_hmat->m10;
  mat_e->G   = p_hmat->G;
  mat_e->kappa = p_hmat->E/(3.0*(1.0-2.0*p_hmat->nu));
  mat_e->devPotFlag = p_hmat->devPotFlag;
  mat_e->volPotFlag = p_hmat->volPotFlag;

  set_properties_constitutive_model(cm_mat,mat_e,NULL);
  construct_elasticity(elast, mat_e, 1);

  this->cm_mat   = cm_mat;
  this->cm_elast = elast;
  err += this->model_dependent_initialization();

  return err;
}

int
Model_parameters::finalization()
{
  /* drop pointer to material (material free'd elsewhere) */
  this->p_hmat = NULL;
  delete (this->cm_mat)->mat_e;
  delete this->cm_mat;
  destruct_elasticity(this->cm_elast);
  delete this->cm_elast;

  /* reset counters/flags */
  this->type = -1;

  /* free model constants */
  this->n_param = -1;
  delete this->model_param;
  this->model_param = NULL;
  this->n_param_index = -1;
  if(this->model_param_index !=NULL)
    delete this->model_param_index;
  this->model_param_index = NULL;
  return 0;
}

/// User defined function that returns the size of the data to be
/// packed/unpacked.
/// Does not modify the CM object or any of the data it holds.
///
/// \param[in] m, CM object with internal data set from the buffer
/// \return size in bytes of the pack/unpack data
int
Constitutive_model::get_size()
{
  return this->vars_list[0][this->model_id].state_variables_get_packed_size();
}

/// User defined function to pack the CM data into a buffer (see pack_data).
/// Does not modify the CM object or any of the data it holds.
///
/// \param[in,out] buffer, a buffer to insert data to
///
/// \param[in,out] pos,    insert position in the buffer. Upon exit - next
///                        insertion position.
/// \return non-zero on error.
int
Constitutive_model::pack(char *buffer, size_t *pos)
{
  return this->vars_list[0][this->model_id].state_variables_pack(buffer, pos);
}

/// User defined function to unpack CM data from a buffer (see also
/// usr_pack, unpack_data).
///
/// \param[in]     buffer, the buffer to read data from
/// \param[in,out] pos,    the position in buffer to begin reading from.
///                        Upon exit - position for next read.
/// \return        non-zero on error.
int
Constitutive_model::unpack(const char *buffer, size_t *pos)
{
  return this->vars_list[0][this->model_id].state_variables_unpack(buffer, pos);
}


// compute stiffness tensor
int constitutive_model_default_update_elasticity(const Constitutive_model *m,
                                                 const double *eF,
                                                 double *L,
                                                 double *S,
                                                 const int compute_stiffness)
{
  int err = 0;

  ELASTICITY *elast = (m->param)->cm_elast; // get elasticity handle
  double *tempS = elast->S; // temporal pointer to update *L, and *S using elast
  double *tempL = elast->L;
  elast->S = S;

  double F[DIM_3x3];
  memcpy(F,eF,sizeof(double)*DIM_3x3);
  if(compute_stiffness)
    elast->L = L;
  else
    elast->L = NULL;

  elast->update_elasticity(elast,F,compute_stiffness);

  elast->S = tempS;
  elast->L = tempL;
  return err;
}

static int compare_mat_id(const void *a, const void *b)
{
  return ((const HOMMAT *)a)->mat_id - ((const HOMMAT *)b)->mat_id;
}

int read_model_parameters_list(const int n_mat,
                               HOMMAT *hmat_list,
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
  double *pF = NULL;
  
  if (n_mat <= 0) return 1;
  int *is_set = (int *) calloc(n_mat, sizeof(int));

  int num_entries = -1;
  HOMMAT *key = (HOMMAT *) calloc(1, sizeof(*key));
  err += scan_for_valid_line(in);
  CHECK_SCANF(in, "%d", &num_entries);

  int i = 0;
  for (i = 0; i < num_entries; i++) {
    int model_type = -1;
    int model_type_uq = 0;

    HOMMAT *p_hmat = NULL;
    err += scan_for_valid_line(in);
    if (feof(in)) break;

    CHECK_SCANF(in, "%d %d", &(key->mat_id), &model_type);
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

    if(key->mat_id>=0)
    {
      /* search for matching pointer in hmat_list (assume unique) */
      p_hmat = static_cast<HOMMAT*>(bsearch(key, hmat_list, n_mat,
                                            sizeof(*hmat_list), compare_mat_id));
  
      /* check for match */
      if (p_hmat != NULL) {
        int idx = p_hmat - hmat_list;
        if (!is_set[idx]) {
          is_set[idx] = 1;
          // construct and initialize this object        
          if(model_type==CM_UQCM)
          {
            model_type_uq = 1;
            err += scan_for_valid_line(in);
            CHECK_SCANF(in, "%d", &model_type);
          }
  
          err += construct_Model_parameters(&(hmat_list[idx].param), 0, model_type);
          hmat_list[idx].param->uqcm = model_type_uq;

          err += hmat_list[idx].param->initialization(hmat_list + idx, model_type);
          err += hmat_list[idx].param->read_param(in);
          hmat_list[idx].param->mat_id = key->mat_id;
        }
      }
    }
    else
    {
      pF = (double *) malloc(sizeof(double)*9);
      err += scan_for_valid_line(in);
      int match = fscanf(in, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", pF+0, pF+1, pF+2, 
                                                                    pF+3, pF+4, pF+5,
                                                                    pF+6, pF+7, pF+8);

      if(match != 9)
      {   
        err++;
        assert(match == 9 && "Did not read expected number of parameters"); 
      }

      err += scan_for_valid_line(in); 
      // not expecting EOF, check and return error if encountered
      if (feof(in)) err ++;
       assert(!feof(in) && "EOF reached prematurely"); 
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

  for (int ia = 0; ia < n_mat; ia++){
    hmat_list[ia].param->pF = pF;
  }

  free(key);
  free(is_set);
  return err;
}

int init_all_constitutive_model(EPS *eps,
                                const int ne,
                                const Element *elem,
                                const int n_mat,
                                const HOMMAT *hmat_list)
{
  int err = 0;
  if (ne <= 0) return 1;

  for(int i = 0; i < ne; i++) {
    /* aliases */
    EPS *p_eps = eps + i;
    const Element *p_el = elem + i;
    const Model_parameters *p_param = hmat_list[p_el->mat[2]].param;

    long n_ip = 0;
    int_point(p_el->toe,&n_ip);
    for (int j = 0; j < n_ip; j++)
      p_eps->model[j].initialization(p_param);
  }

  plasticity_model_set_orientations(eps, ne, elem, n_mat, hmat_list); // nothing will happen if there is no use of the crystal plasticity model
  return err;
}

int constitutive_model_reset_state(EPS *eps,
                                   const int ne,
                                   const Element *elem)
{
  int err = 0;
  if (ne <= 0) return 1;

  for (int i = 0; i < ne; i++) {
    const Element *p_el = elem + i;
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
int constitutive_model_save_state_vars_to_temporal(FieldVariables *fv,
                                                   Grid *grid)
{
  int err = 0;
  State_variables *var = fv->temporal->var;
  const Element *elem = grid->element;
  for(int eid=0; eid<grid->ne; eid++)
  {
    long nint = 1;
    int_point(elem[eid].toe,&nint);
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
int constitutive_model_update_np1_state_vars_to_temporal(FieldVariables *fv,
                                                         Grid *grid)
{
  int err = 0;
  State_variables *var = fv->temporal->var;

  const Element *elem = grid->element;

  for(int eid=0; eid<grid->ne; eid++)
  {
    long nint = 1;
    int_point(elem[eid].toe,&nint);
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
int constitutive_model_reset_state_using_temporal(FieldVariables *fv,
                                                  Grid *grid)
{
  int err = 0;
  State_variables *var = fv->temporal->var;

  const Element *elem = grid->element;

  for(int eid=0; eid<grid->ne; eid++)
  {
    long nint = 1;
    int_point(elem[eid].toe,&nint);
    for (int ip = 0; ip < nint; ip++)
    {
      Constitutive_model *m = &(fv->eps[eid].model[ip]);
      m->param->reset_state_vars_using_temporal(m, var+m->model_id);
    }
  }
  return err;
}

int constitutive_model_update_time_steps(const Element *elem,
                                         Node *node,
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
    const Element *p_el = elem + i;
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

int constitutive_model_test(const HOMMAT *hmat, double *L_in, int Print_results)
{
  int err = 0;
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
                             double *Fr_in,
                             double *eFnMT_in,
                             double *eFn_in,
                             double *M_in,
                             double *FrTFr_in,
                             double *eFnM_in,
                             double *S_in,
                             double *L_in,
                             double *dMdu_all,
                             const double Jn)
{
  int err = 0;
  const int nne = fe->nne;
  const int nsd = fe->nsd;

  Tensor<2, 3, double*> Fr(Fr_in);
  Tensor<2, 3, double*> eFnMT(eFnMT_in);
  Tensor<2, 3, double*> eFn(eFn_in);
  Tensor<2, 3, double*> M(M_in);
  Tensor<2, 3, double*> FrTFr(FrTFr_in);
  Tensor<2, 3, double*> eFnM(eFnM_in);
  Tensor<2, 3, double*> S(S_in);
  Tensor<4, 3, double*> L(L_in);

  Tensor<2,3,double> MTeFnT_sAA_eFn, MTeFnT_sAA_eFnM,
            sBB,sCC, dCdu;

  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const int id_ab = idx_4_gen(a,b,0,0,nne,nsd,nsd,nsd);
      Tensor<2, 3, double*> ST_ab((fe->ST)+id_ab);

      Tensor<2,3,double> AA =  Fr(k,i).to(i,k)*ST_ab(k,j);
      Tensor<2,3,double> sAA = 0.5*(AA(i,j) + AA(j,i).to(i,j));

      MTeFnT_sAA_eFn(i,j) = eFnMT(i,k) * sAA(k,l) * eFn(l,j);
      MTeFnT_sAA_eFnM(i,j) = MTeFnT_sAA_eFn(i,k) * M(k,j);

      for(int w=0; w<nne; w++)
      {
        for(int g=0; g<nsd; g++)
        {
          const int id_wg = idx_4_gen(w,g,0,0,nne,nsd,nsd,nsd);
          Tensor<2, 3, double*> ST_wg((fe->ST)+id_wg);

          Tensor<2, 3, double*> dMdu(dMdu_all + id_wg);

          Tensor<2,3,double> BB =  Fr(k,i).to(i,k) * ST_wg(k,j);
          Tensor<2,3,double> sBB = 0.5 * (BB(i,j) + BB(j,i).to(i,j));

          Tensor<2,3,double> CC =  ST_ab(k,i).to(i,k) * ST_wg(k,j);
          sCC(i,j) = 0.5 * (CC(i,j) + CC(j,i).to(i,j));

          // compute dCdu
          Tensor<2,3,double> MTeFnT_FrTFreFndMdu = eFnMT(i,k)*FrTFr(k,l)*eFn(l,o)*dMdu(o,j);
          dCdu(i,j) = 0.5 * (MTeFnT_FrTFreFndMdu(i,j)
                                    + MTeFnT_FrTFreFndMdu(j,i).to(i,j)) + eFnMT(i,k)*sBB(k,l)*eFnM(l,j);


          // compute MTeFnT_sAA_eFnM:L:dCdu
          Tensor<2,3,double> L_dCdu = L(i,j,k,l) * dCdu(k,l);
          double MTeFnT_sAA_eFnM_L_dCdu = MTeFnT_sAA_eFnM(i,j) * L_dCdu(i,j);

          // compute MTeFnT_sCC_eFnM
          Tensor<2,3,double> MTeFnT_sCC_eFnM = eFnMT(i,k) * sCC(k,l) * eFnM(l,j);

          // compute MTeFnT_sCC_eFnM:S
          double MTeFnT_sCC_eFnM_S = MTeFnT_sCC_eFnM(i,j) * S(i,j);

          // compute MTeFnT_sAA_eFndMdu
          Tensor<2,3,double> MTeFnT_sAA_eFndMdu = MTeFnT_sAA_eFn(i,k) * dMdu(k,j);
          // compute MTeFnT_sAA_eFndMdu:S
            double sMTeFnT_sAA_eFndMdu_S = 0.5*(MTeFnT_sAA_eFndMdu(i,j)
                                                       + MTeFnT_sAA_eFndMdu(j,i).to(i,j))* S(i,j);

          const int lk_idx = idx_K(a,b,w,g,nne,nsd);

          lk[lk_idx] += 1.0/Jn*fe->detJxW*(MTeFnT_sAA_eFnM_L_dCdu + 2.0*sMTeFnT_sAA_eFndMdu_S + MTeFnT_sCC_eFnM_S);
        }
      }
    }
  }

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
                            double *Fr_in,
                            double *eFnMT_in,
                            double *eFnM_in,
                            double *S_in,
                            const double Jn)
{
  int err = 0;
  const int nne = fe->nne;
  const int nsd = fe->nsd;

  Tensor<2, 3, double*> Fr(Fr_in);
  Tensor<2, 3, double*> eFnMT(eFnMT_in);
  Tensor<2, 3, double*> eFnM(eFnM_in);
  Tensor<2, 3, double*> S(S_in);

  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const int id_ab = idx_4_gen(a,b,0,0,nne,nsd,nsd,nsd);
      Tensor<2, 3, double*> ST_ab((fe->ST)+id_ab);

      Tensor<2,3,double> AA = Fr(k,i)*ST_ab(k,j);
      Tensor<2,3,double> sAA = 0.5*(AA(i,j)+AA(j,i));
      double MTeFnT_sAA_eFnM_S = eFnMT(i,k)*sAA(k,l)*eFnM(l,j)*S(i,j);

      int fe_id = a*nsd + b;
      f[fe_id] += 1.0/Jn*fe->detJxW*MTeFnT_sAA_eFnM_S;
    }
  }

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
int constitutive_model_update_output_variables(Grid *grid,
                                               MaterialProperty *mat,
                                               FieldVariables *FV,
                                               LoadingSteps *load,
                                               PGFem3D_opt *opts,
                                               Multiphysics *mp,
                                               int mp_id,
                                               const double dt,
                                               double alpha)
{
  int err = 0;
  
  int total_Lagrangian = 1;
  if(opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;

  FieldVariables *fv = FV+mp_id;
  SIG *sig = fv->sig;
  EPS *eps = fv->eps;
  Node *node = grid->node;
  Element *elem = grid->element;

  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  static const double eye[DIM_3x3] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  /* deformation gradient */
  Tensor<2,3,double> F, eF, pF, S;

  MATERIAL_ELASTICITY mat_e_new;
  MATERIAL_ELASTICITY *mat_e_in;

  for (int ie = 0; ie < grid->ne; ie++)
  {
    FEMLIB fe(ie,elem,node,0,total_Lagrangian);

    memset(sig[ie].el.o,0,6*sizeof(double));
    memset(eps[ie].el.o,0,6*sizeof(double));
    double volume = 0.0;

    Matrix<double> Tnp1, Tn, Tnm1;
    FieldVariables *fv_h = NULL;

    if(is_it_couple_w_thermal >= 0)
    {
      fv_h = fv->fvs[is_it_couple_w_thermal];
      Tnm1.initialization(fe.nne,1);
      Tnp1.initialization(fe.nne,1);
      Tn.initialization(fe.nne,1);

      // compute temperature for this element for each nodal point
      int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
      err += get_nodal_temperatures(&fe, grid, fv_h, load, mp_cp_id,
                                    Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,fv->subdivision_factor_np1,fv->subdivision_factor_n);
    }

    if(is_it_couple_w_chemical >= 0){}

    for(int ip=0; ip<fe.nint; ip++)
    {
      fe.elem_basis_V(ip+1);

      Constitutive_model *m = &(eps[ie].model[ip]);
      const Model_parameters *func = m->param;

      err += func->get_F(m,  F.data,1);
      err += func->get_pF(m,pF.data,1);

      void *ctx = NULL;
      if(is_it_couple_w_thermal>=0)
      {
        // compute temperature at the integration point
        double T0 = fv_h->u0;
        double hFnm1[9],hFn[9],hFnp1[9];
        err += compute_temperature_at_ip(&fe,grid,mat,T0,
                                         Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                         hFnp1,hFn,hFnm1);
        Tensor<2,3,double> pFI,hFI;
        inv(pF,pFI);
        inv3x3(hFnp1,hFI.data);
        eF(i,j) = F(i,k)*hFI(k,l)*pFI(l,j);
      }
      else
        err += func->get_eF(m, eF.data,1);

      err += construct_model_context(&ctx, m->param->type, F.data,dt,alpha,eF.data,-1);
          
      if((m->param)->uqcm)
      {
        double *x_ip = (fe.x_ip).m_pdata;

        ELASTICITY *elast = (m->param)->cm_elast;
        mat_e_in = elast->mat;
        err += material_properties_elasticity_at_ip(mat_e_in, &mat_e_new, x_ip[0], x_ip[1], x_ip[2]);
        elast->mat = &mat_e_new; // should be replaced by original mat_e_in after computation
        err += (m->param)->update_elasticity(m,ctx,NULL,S.data,0);

        elast->mat = mat_e_in;
      }
      else
        err += (m->param)->update_elasticity(m,ctx,NULL,S.data,0);

      err += (m->param)->destroy_ctx(&ctx);
      // <-- update elasticity part

      /* get aliases to Matrix data for simpler access */
      const double *Sd = S.data;
      const double *eFd = eF.data;
      const double eJ = det3x3(eFd);

      /* store symmetric part of S (PK2) */
      sig[ie].il[ip].o[0] = Sd[idx_2(0,0)]; /* XX */
      sig[ie].il[ip].o[1] = Sd[idx_2(1,1)]; /* YY */
      sig[ie].il[ip].o[2] = Sd[idx_2(2,2)]; /* ZZ */
      sig[ie].il[ip].o[3] = Sd[idx_2(1,2)]; /* YZ */
      sig[ie].il[ip].o[4] = Sd[idx_2(0,2)]; /* XZ */
      sig[ie].il[ip].o[5] = Sd[idx_2(0,1)]; /* XY */

      /* store total deformation */
      memcpy(eps[ie].il[ip].F, F.data, DIM_3x3 * sizeof(double));

      /* store the hardening parameter */
      err += func->get_hardening(m, &eps[ie].dam[ip].wn,2);

      /* compute/store the plastic strain variable */
      err += func->get_plast_strain_var(m, &eps[ie].dam[ip].Xn);

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
      sig[ie].el.o[0] += fe.detJxW*sigma[idx_2(0,0)]; /* XX */
      sig[ie].el.o[1] += fe.detJxW*sigma[idx_2(1,1)]; /* YY */
      sig[ie].el.o[2] += fe.detJxW*sigma[idx_2(2,2)]; /* ZZ */
      sig[ie].el.o[3] += fe.detJxW*sigma[idx_2(1,2)]; /* YZ */
      sig[ie].el.o[4] += fe.detJxW*sigma[idx_2(0,2)]; /* XZ */
      sig[ie].el.o[5] += fe.detJxW*sigma[idx_2(0,1)]; /* XY */

      /* Compute the logarithmic strain e = 1/2(I - inv(FF'))*/
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                  DIM_3, DIM_3, DIM_3, 1.0, eFd, DIM_3, eFd, DIM_3,
                  0.0, temp, DIM_3);
      err += inv3x3(temp, temp_I);
      /* e <-- temp is the Euler strain */
      for(int ia = 0; ia < DIM_3x3; ia++)
        temp[ia] = 0.5*(eye[ia]-temp_I[ia]);

      /* store symmetric part (also Eng. strain) */
      eps[ie].el.o[0] += fe.detJxW*temp[idx_2(0,0)];
      eps[ie].el.o[1] += fe.detJxW*temp[idx_2(1,1)];
      eps[ie].el.o[2] += fe.detJxW*temp[idx_2(2,2)];

      eps[ie].el.o[3] += fe.detJxW*2.0*temp[idx_2(1,2)];
      eps[ie].el.o[4] += fe.detJxW*2.0*temp[idx_2(0,2)];
      eps[ie].el.o[5] += fe.detJxW*2.0*temp[idx_2(0,1)];

    }
    for(int ia=0; ia<6; ia++)
    {
      sig[ie].el.o[ia] = sig[ie].el.o[ia]/volume;
      eps[ie].el.o[ia] = eps[ie].el.o[ia]/volume;
    }
  }

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
                                                 double *pFnp1_in,
                                                 double *pFn_in,
                                                 double *Fnp1_in,
                                                 double *Fn_in,
                                                 double *hFnp1_in,
                                                 double *hFn_in,
                                                 const int is_it_couple_w_thermal,
                                                 const double alpha,
                                                 const double dt_alpha_1_minus_alpha,
                                                 FEMLIB *fe,
                                                 const int npa)
{
  // Total Lagrangian based
  int err = 0;
  const int nsd = fe->nsd;

  Tensor<2,3,double> M = {},eFnpa = {},pFnpa,pFnpa_I,hFnpa,Fnpa,S = {},MT;
  Tensor<2, 3, double*> pFnp1(pFnp1_in),pFn(pFn_in),
                        Fnp1(Fnp1_in),Fn(Fn_in),
                        hFnp1(hFnp1_in),hFn(hFn_in);


  int compute_stiffness = 0;

  mid_point_rule(pFnpa.data, pFn.data, pFnp1.data, alpha, nsd*nsd);
  hFnpa = hFnp1(i,j);
//  mid_point_rule(hFnpa.data, hFn.data, hFnp1.data, alpha, nsd*nsd);
  mid_point_rule( Fnpa.data,  Fn.data,  Fnp1.data, alpha, nsd*nsd);

  if(is_it_couple_w_thermal>=0)
  {
    Tensor<2,3,double> hFnpa_I;

    err += inv(hFnpa, hFnpa_I);
    err += inv(pFnpa, pFnpa_I);
    M = hFnpa_I(i,k)*pFnpa_I(k,j);
  }
  else
    err += inv(pFnpa, M);

  eFnpa = Fnpa(i,k)*M(k,j);

  MT = M(i,j).to(j,i);

  {
    // check that deformation is invertible -> J > 0
    int terr = 0;
    getJacobian(Fnpa.data, ii, &terr);
    err += terr;
  }

  void *ctx;

  err += construct_model_context(&ctx, m->param->type, Fnp1.data,0.0,alpha,eFnpa.data,npa);
  err += (m->param)->update_elasticity(m,ctx,NULL,S.data,compute_stiffness);
  err += m->param->destroy_ctx(&ctx);

  if(err==0)
  {
    double Jn = 1.0;
    err += compute_residual_vector(f,fe,Fnpa.data,MT.data,M.data,S.data,Jn/dt_alpha_1_minus_alpha);
  }

  return err;
}

int cm_get_subdivision_parameter(double *subdiv_param,
                                 const int ne,
                                 const Element *elem,
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
      max_val = std::max(max_val, cur_val);
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
                                              Grid *grid,
                                              MaterialProperty *mat,
                                              FieldVariables *fv,
                                              Solver *sol,
                                              LoadingSteps *load,
                                              CRPL *crpl,
                                              const PGFem3D_opt *opts,
                                              Multiphysics *mp,
                                              int mp_id,
                                              double dt)
{
  int err = 0;
  double alpha = sol->alpha;
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

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

  double *u = (double *) malloc(sizeof(*u)*nne*nsd);
  double *dMdu_all = (double *) malloc(sizeof(*dMdu_all)*DIM_3x3*nne*nsd);
  memset(dMdu_all,0,DIM_3x3*nne*nsd*sizeof(double));

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];
  }

  Tensor<2,3,double> Fn,Fr,Fnp1,pFn,pFnp1,S = {},eFnpa = {},pFnpa,pFnpa_I,
            eFn,M = {},eFnM,eFnMT,FrTFr = {},hFnp1,hFn,hFnpa,hFnpa_I;

  Tensor<4,3,double> L;

  Matrix<double> Tnp1, Tn, Tnm1;
  FieldVariables *fv_h = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    fv_h = fv->fvs[is_it_couple_w_thermal];
    Tnm1.initialization(fe->nne,1);
    Tnp1.initialization(fe->nne,1);
    Tn.initialization(fe->nne,1);

    // compute temperature for this element for each nodal point

    int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
    err += get_nodal_temperatures(fe, grid, fv_h, load, mp_cp_id,
                                  Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                  fv->subdivision_factor_np1,fv->subdivision_factor_n);
  }

  if(is_it_couple_w_chemical >= 0)
  {}

  int compute_stiffness = 1;

  for(int ip = 1; ip<=fe->nint; ip++)
  {

    double Jn = 1.0; // det(F2[Fn], Jn);
    double hJ = 1.0;
    double pJ = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    fe->update_deformation_gradient(ndofn,u,Fnp1.data);

    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);

    // get a shortened pointer for simplified CM function calls
    const Model_parameters *func = m->param;

    hFnpa_I = identity(i,j);
    if(is_it_couple_w_thermal >= 0)
    {
      double T0 = fv_h->u0;
      double hFnm1[9];
      err += compute_temperature_at_ip(fe,grid,mat,T0,
                                       Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1);

      hFnpa = hFnp1(i,j);
      hJ = det(hFnp1);
      //mid_point_rule(hFnpa.data, hFn.data, hFnp1.data, alpha, DIM_3x3);
      inv(hFnpa,hFnpa_I);
    }

    err += func->get_pF(m,pFnp1.data,2);
    pJ = det(pFnp1);

    err += func->get_pF(m,pFn.data,1);
    err += func->get_F(m,  Fn.data,1);


    mid_point_rule(pFnpa.data, pFn.data, pFnp1.data, alpha, DIM_3x3);
    mid_point_rule(   Fr.data,  Fn.data,  Fnp1.data, alpha, DIM_3x3);

    err += inv(pFnpa, pFnpa_I);

    M = hFnpa_I(i,k)*pFnpa_I(k,j);
    eFnpa = Fr(i,k)*M(k,j);
    FrTFr = Fr(k,i)*Fr(k,j);

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, Fnp1.data,dt,alpha, eFnpa.data,
                                                  hFn.data,hFnp1.data,1);      
    else
      err += construct_model_context(&ctx, m->param->type, Fnp1.data,dt,alpha, eFnpa.data,1);
        
    err += m->param->compute_dMdu(m, ctx, fe->ST, nne, ndofn, dMdu_all);

    // --> update elasticity part
    memset(L.data, 0, sizeof(double)*DIM_3x3x3x3);
    memset(S.data, 0, sizeof(double)*DIM_3x3);

    err += (m->param)->update_elasticity(m,ctx,L.data,S.data,compute_stiffness);
    // <-- update elasticity part
    err += m->param->destroy_ctx(&ctx);

    if(err!=0)
      break;

    // --> start computing tagent
    // total Lagrangian
    eFn = identity(i,j);
    eFnM = M(i,j);
    eFnMT(j,i) = M(i,j).to(j,i);

    Jn = Jn/pJ/hJ;
    err += compute_stiffness_matrix(lk,fe,
                                    Fr.data,eFnMT.data,eFn.data,M.data,FrTFr.data,eFnM.data,S.data,
                                    L.data,dMdu_all,Jn);
  }

  free(u);
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
int stiffness_el_constitutive_model_1f(FEMLIB *fe,
                                       double *lk,
                                       double *r_e,
                                       Grid *grid,
                                       MaterialProperty *mat,
                                       FieldVariables *fv,
                                       Solver *sol,
                                       LoadingSteps *load,
                                       CRPL *crpl,
                                       const PGFem3D_opt *opts,
                                       Multiphysics *mp,
                                       int mp_id,
                                       double dt)
{
  int err = 0;
  double alpha = -1.0; // if alpha < 0, no inertia
  int total_Lagrangian = 0;
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

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

  double *u = (double *) malloc(sizeof(*u)*nne*nsd);
  double *dMdu_all = (double *) malloc(sizeof(*dMdu_all)*DIM_3x3*nne*nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];
  }

  Tensor<2,3,double> Fr, Fnp1 = {}, pFnp1, S,
            eFn = {}, M = {}, eFnM = {},eFnMT, FrTFr = {},
            hFn, hFnp1, hFnp1_I;

  Tensor<4,3,double> L;

  Matrix<double> Tnp1, Tn, Tnm1;
  FieldVariables *fv_h = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    fv_h = fv->fvs[is_it_couple_w_thermal];
    Tnm1.initialization(fe->nne,1);
    Tnp1.initialization(fe->nne,1);
    Tn.initialization(fe->nne,1);

    // compute temperature for this element for each nodal point
    int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
    err += get_nodal_temperatures(fe, grid, fv_h, load, mp_cp_id,
                                  Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                  fv->subdivision_factor_np1,fv->subdivision_factor_n);
  }

  if(is_it_couple_w_chemical >=0)
  {}

  int compute_stiffness = 1;

  MATERIAL_ELASTICITY mat_e_new;
  MATERIAL_ELASTICITY *mat_e_in;

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double Jn = 1.0;
    double hJ = 1.0;
    double pJ = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    fe->update_deformation_gradient(ndofn,u,Fr.data);

    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);

    // get a shortened pointer for simplified CM function calls
    const Model_parameters *func = m->param;

    // --> update deformations due to coupled physics
    if(is_it_couple_w_thermal >= 0)
    {
      // compute temperature at the integration point
      double T0 = fv_h->u0;

      Tensor<2,3,double> hFnm1;
      err += compute_temperature_at_ip(fe,grid,mat,T0,
                                       Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1.data);
      hJ = det(hFnp1);
      inv(hFnp1,hFnp1_I);
    }

    // --> update plasticity part
    if(total_Lagrangian)
    {
      if(sup->multi_scale)
        cm_add_macro_F(sup,Fr.data);

      // Total Lagrangian formulation Fn = 1, Fnp1 = Fr
      eFn = identity(i,j);
      Fnp1 = Fr(i,j);
      hFn = identity(i,j);
    }
    else
    {
      if(sup->multi_scale)
      {
        PGFEM_printerr("Multi-scale formulation does not support UL!\n");
        PGFEM_Abort();
      }

      Tensor<2,3,double> Fn;
      err += m->param->get_F(m, Fn.data,1);
      Fnp1 = Fr(i,k)*Fn(k,j); // Fn+1 = Fr*Fn
      Jn = det(Fn);
    }

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, Fnp1.data,dt,alpha, NULL,
                                                  hFn.data,hFnp1.data,-1);
    else
      err += construct_model_context(&ctx, m->param->type, Fnp1.data,dt,alpha, NULL,-1);
  
    err += func->compute_dMdu(m, ctx, fe->ST, nne, ndofn, dMdu_all);
    err += func->get_pF(m,pFnp1.data,2);

    pJ = det(pFnp1);


    if(total_Lagrangian) // Total Lagrangian formulation, all xFn = 1
    {
      if(is_it_couple_w_thermal >= 0)
      {
        Tensor<2,3,double> pFnp1_I;
        inv(pFnp1, pFnp1_I);
        M = hFnp1_I(i,k)*pFnp1_I(k,j);
      }
      else
        inv(pFnp1, M);

      eFnM = M(i,j);
      eFnMT(j,i) = M(i,j).to(j,i);
    }
    else
    {
      if(is_it_couple_w_thermal>=0)
      {
        Tensor<2,3,double> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += func->get_pF(m,pFn.data,1);
        M = pFn(i,k)*hFn(k,l)*hFnp1_I(l,o)*pFnp1_I(o,j);

        Tensor<2,3,double> hFn_I;
        int stepno = 1; // 0 = time step = n-1
                        // 1 = time step = n
                        // 2 = time step = n+1
        inv(hFn, hFn_I);
        err += m->param->get_eF_of_hF(m,eFn.data,hFn_I.data,stepno);
      }
      else
      {
        Tensor<2,3,double> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += func->get_pF(m,pFn.data,1);
        M = pFn(i,k)*pFnp1_I(k,j);
        err += m->param->get_eF(m,eFn.data,1);
      }
      eFnM = eFn(i,k)*M(k,j);
      eFnMT(j,i) = eFnM(i,j).to(j,i);
    }
    // <-- update plasticity part

    FrTFr = Fr(k,i)*Fr(k,j);

    // --> update elasticity part
    memset(L.data, 0, sizeof(double)*DIM_3x3x3x3);
    memset(S.data, 0, sizeof(double)*DIM_3x3);

    if((m->param)->uqcm)
    {
      double *x_ip = (fe->x_ip).m_pdata;
      ELASTICITY *elast = (m->param)->cm_elast;
      mat_e_in = elast->mat;
      err += material_properties_elasticity_at_ip(mat_e_in, &mat_e_new, x_ip[0], x_ip[1], x_ip[2]);
      elast->mat = &mat_e_new; // should be replaced by original mat_e_in after computation
      err += (m->param)->update_elasticity(m,ctx,L.data,S.data,compute_stiffness);
      elast->mat = mat_e_in;
    }
    else
      err += (m->param)->update_elasticity(m,ctx,L.data,S.data,compute_stiffness);
    // <-- update elasticity part

    err += func->destroy_ctx(&ctx);
    if(err!=0)
      break;

    // start computing tagent
    Jn = Jn/hJ/pJ;
    err += compute_stiffness_matrix(lk,fe,
                                    Fr.data,eFnMT.data,eFn.data,M.data,FrTFr.data,eFnM.data,S.data,
                                    L.data,dMdu_all,Jn);
  }
  free(u);

  /* check diagonal for zeros/nans */
  for (int a = 0; a < nne; a++) {
    for (int b = 0; b < nsd; b++) {
      if ( !isnormal(lk[idx_K(a,b,a,b,nne,nsd)]) ) err++;
    }
  }

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
int stiffness_el_constitutive_model_3f(FEMLIB *fe,
                                       double *lk,
                                       double *r_e,
                                       Grid *grid,
                                       MaterialProperty *mat,
                                       FieldVariables *fv,
                                       Solver *sol,
                                       LoadingSteps *load,
                                       CRPL *crpl,
                                       const PGFem3D_opt *opts,
                                       Multiphysics *mp,
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
  int Pno   = fv->npres;
  int Vno   = fv->nVol;
  SUPP sup = load->sups[mp_id];

  Matrix<double> u(nne*nsd, 1);
  Matrix<double> dMdu(DIM_3x3*nne*nsd,1);
  Matrix<double> dMdt(DIM_3x3*Vno,1);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u.m_pdata[a*nsd+b] = r_e[a*ndofn+b];
  }

  Tensor<2,3,double> Fr, Fnp1 = {},pFnp1,
            eSd = {}, M = {}, eFn = {},
            hFn,hFnp1,hFnp1_I;

  Tensor<4,3,double> Ld={};

  Matrix<double> Tnp1, Tn, Tnm1;
  FieldVariables *fv_h = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    fv_h = fv->fvs[is_it_couple_w_thermal];
    Tnm1.initialization(fe->nne,1);
    Tnp1.initialization(fe->nne,1);
    Tn.initialization(fe->nne,1);

    // compute temperature for this element for each nodal point
    int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
    err += get_nodal_temperatures(fe, grid, fv_h, load, mp_cp_id,
                                  Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                  fv->subdivision_factor_np1,fv->subdivision_factor_n);
  }

  if(is_it_couple_w_chemical >=0)
  {}

  Matrix<double> Kuu(nne*nsd, nne*nsd, 0.0), Ktu(Vno, nne*nsd, 0.0), Kpu(Pno, nne*nsd, 0.0);
  Matrix<double> Kut(nne*nsd, Vno), Ktt(Vno, Vno, 0.0), Kpt(Pno, Vno, 0.0), Kup(nne*nsd, Pno, 0.0), Ktp(Vno, Pno, 0.0);

  Matrix<double> Nt(Vno,1), Np(Pno,1);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double Jn  = 1.0; // if upated Lagrangian, Jn = det(Fn), later updated
    double hJ  = 1.0;
    double pJ  = 1.0;
    double tJn = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    fe->update_deformation_gradient(ndofn,u.m_pdata,Fr.data);

    fe->elem_shape_function(ip,Pno, Np.m_pdata);
    fe->elem_shape_function(ip,Vno, Nt.m_pdata);

    double theta_r = 0.0;
    double theta_n = 0.0;
    double Pnp1    = 0.0;

    for(int ia=1; ia<=Pno; ia++)
      Pnp1 += (fv->tf.P_np1(eid+1, ia) + fv->tf.dP(eid+1, ia))*Np(ia);

    for(int ia=1; ia<=Vno; ia++)
    {
      theta_r += (fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia))*Nt(ia);
      theta_n += fv->tf.V_n(  eid+1, ia)*Nt(ia);
    }

    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);

    // get a shortened pointer for simplified CM function calls
    const Model_parameters *mp = m->param;
    ELASTICITY *elasticity = mp->cm_elast;

    // --> update deformations due to coupled physics
    if(is_it_couple_w_thermal >= 0)
    {
      // compute temperature at the integration point
      double T0 = fv_h->u0;

      Tensor<2,3,double> hFnm1;
      err += compute_temperature_at_ip(fe,grid,mat,T0,
                                       Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1.data);
      hJ = det(hFnp1);
      inv(hFnp1,hFnp1_I);
    }

    // --> update plasticity part
    if(total_Lagrangian)
    {
      if(sup->multi_scale)
        cm_add_macro_F(sup,Fr.data);

      // Total Lagrangian formulation Fn = 1, Fnp1 = Fr
      eFn = identity(i,j);
      Fnp1 = Fr(i,j);
      hFn = identity(i,j);
    }
    else
    {
      if(sup->multi_scale)
      {
        PGFEM_printerr("Multi-scale formulation does not support UL!\n");
        PGFEM_Abort();
      }

      Tensor<2,3,double> Fn;
      err += m->param->get_F(m, Fn.data,1);
      Fnp1 = Fr(i,k)*Fn(k,j); // Fn+1 = Fr*Fn
      Jn = det(Fn);
    }

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, Fnp1.data,dt,alpha, NULL,
                                                  hFn.data,hFnp1.data,-1);
    else
      err += construct_model_context(&ctx, m->param->type, Fnp1.data,dt,alpha, NULL,-1);
  
    err += mp->get_pF(m,pFnp1.data,2);

    pJ = det(pFnp1);


    if(total_Lagrangian) // Total Lagrangian formulation, all xFn = 1
    {
      if(is_it_couple_w_thermal >= 0)
      {
        Tensor<2,3,double> pFnp1_I;
        inv(pFnp1, pFnp1_I);
        M = hFnp1_I(i,k)*pFnp1_I(k,j);
      }
      else
        inv(pFnp1, M);

      eFn = identity(i,j);
    }
    else
    {
      if(is_it_couple_w_thermal>=0)
      {
        Tensor<2,3,double> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += mp->get_pF(m,pFn.data,1);
        M = pFn(i,k)*hFn(k,l)*hFnp1_I(l,o)*pFnp1_I(o,j);

        Tensor<2,3,double> hFn_I;
        int stepno = 1; // 0 = time step = n-1
                        // 1 = time step = n
                        // 2 = time step = n+1
        inv(hFn, hFn_I);
        err += m->param->get_eF_of_hF(m,eFn.data,hFn_I.data,stepno);
      }
      else
      {
        Tensor<2,3,double> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += mp->get_pF(m,pFn.data,1);
        M = pFn(i,k)*pFnp1_I(k,j);
        err += m->param->get_eF(m,eFn.data,1);
      }
    }

    err += mp->compute_dMdu(m, ctx, fe->ST, nne, ndofn, dMdu.m_pdata);
    err += mp->compute_dMdt(m, ctx, fe->ST, Vno, dMdt.m_pdata);

    // <-- update plasticity part
    err += mp->compute_dev_stress(m,ctx,eSd.data);
    err += mp->compute_dev_tangent(m,ctx,Ld.data);
    double kappa = elasticity->mat->kappa;              
    double dUd_theta, d2Ud_theta2;
    
    double eJ = ttl::det(eFn);
    double MJ = ttl::det(M);
    double theta_e = theta_r*eJ*MJ;

    elasticity->compute_dudj(&dUd_theta, theta_e);
    elasticity->compute_d2udj2(&d2Ud_theta2, theta_e);
    dUd_theta = dUd_theta*kappa;
    d2Ud_theta2 = d2Ud_theta2*kappa;
    // <-- update elasticity part
    
    err += mp->destroy_ctx(&ctx);

    if(err!=0)
      break;

    tJn = Jn;
    Jn = Jn/pJ/hJ;

    Var_Carrier vc;

    vc.set_tenosrs(Fr.data, eFn.data, M.data, pFnp1.data, eSd.data, Ld.data);
    vc.set_scalars(theta_r, theta_n, tJn, Jn, Pnp1, dUd_theta, d2Ud_theta2);

    err += compute_Kuu(fe, Kuu.m_pdata, dMdu.m_pdata, vc);
    err += compute_Kut(fe, Kut.m_pdata, dMdt.m_pdata, vc, Vno, Nt.m_pdata);
    err += compute_Kup(fe, Kup.m_pdata, vc, Pno, Np.m_pdata);
    err += compute_Ktu(fe, Ktu.m_pdata, dMdu.m_pdata, vc, Vno, Nt.m_pdata);
    err += compute_Ktp(fe, Ktp.m_pdata, vc, Vno, Nt.m_pdata, Pno, Np.m_pdata);
    err += compute_Ktt(fe, Ktt.m_pdata, dMdt.m_pdata, vc, Vno, Nt.m_pdata);
  }

  for(int ia=1; ia<=Vno; ia++)
    for(int ib=1; ib<=Pno; ib++)
      Kpt(ib,ia) = Ktp(ia,ib);

  for(int ia=1; ia<=nne*nsd; ia++)
    for(int ib=1; ib<=Pno; ib++)
      Kpu(ib,ia) = Kup(ia,ib);
 
  err += condense_K_3F_to_1F(lk, nne, nsd, Pno, Vno,
                             Kuu.m_pdata, Kut.m_pdata, Kup.m_pdata,
                             Ktu.m_pdata, Ktt.m_pdata, Ktp.m_pdata,
                             Kpu.m_pdata, Kpt.m_pdata, NULL);

  // check diagonal for zeros/nans
  for (int a = 0; a < nne; a++) {
    for (int b = 0; b < nsd; b++) {
      if ( !isnormal(lk[idx_K(a,b,a,b,nne,nsd)]) ) err++;
    }
  }

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
                                    Grid *grid,
                                    MaterialProperty *mat,
                                    FieldVariables *fv,
                                    Solver *sol,
                                    LoadingSteps *load,
                                    CRPL *crpl,
                                    const PGFem3D_opt *opts,
                                    Multiphysics *mp,
                                    int mp_id,
                                    double dt)
{
  int err = 0;
  if(opts->analysis_type==CM)
    err += stiffness_el_constitutive_model_1f(fe,lk,r_e,grid,mat,fv,sol,load,crpl,opts,mp,mp_id,dt);

  if(opts->analysis_type==CM3F)
    err += stiffness_el_constitutive_model_3f(fe,lk,r_e,grid,mat,fv,sol,load,crpl,opts,mp,mp_id,dt);

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
                                              Grid *grid,
                                              MaterialProperty *mat,
                                              FieldVariables *fv,
                                              Solver *sol,
                                              LoadingSteps *load,
                                              CRPL *crpl,
                                              const PGFem3D_opt *opts,
                                              Multiphysics *mp,
                                              const double *dts,
                                              int mp_id,
                                              double dt)
{
  int err = 0;
  double alpha = sol->alpha;
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

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

  double *u       = (double *) malloc(sizeof(double)*nne*nsd);
  double *f_npa   = (double *) malloc(sizeof(double)*nne*nsd);
  double *f_nm1pa = (double *) malloc(sizeof(double)*nne*nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];
  }

  Tensor<2,3,double> Fnp1,Fn,Fnm1,pFnp1,pFn,pFnm1,
        hFnp1,hFn,hFnm1;

  Matrix<double> Tnp1, Tn, Tnm1;
  FieldVariables *fv_h = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    fv_h = fv->fvs[is_it_couple_w_thermal];
    Tnm1.initialization(fe->nne,1);
    Tnp1.initialization(fe->nne,1);
    Tn.initialization(fe->nne,1);

    // compute temperature for this element for each nodal point
    int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
    err += get_nodal_temperatures(fe, grid, fv_h, load, mp_cp_id,
                                  Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                  fv->subdivision_factor_np1,fv->subdivision_factor_n);
  }

  if(is_it_couple_w_chemical >=0)
  {}

  memset(f_npa, 0, sizeof(double)*nne*ndofn);
  memset(f_nm1pa, 0, sizeof(double)*nne*ndofn);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double hJ = 1.0;
    double pJ = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    fe->update_deformation_gradient(ndofn,u,Fnp1.data);

    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);

    if(is_it_couple_w_thermal >= 0)
    {
      double T0 = fv_h->u0;
      err += compute_temperature_at_ip(fe,grid,mat,T0,
                                       Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1.data);
      hJ = det(hFnp1);
    }

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, Fnp1.data,dts[DT_NP1],alpha, NULL,
                                                  hFn.data,hFnp1.data, -1);      
    else
      err += construct_model_context(&ctx, m->param->type, Fnp1.data,dts[DT_NP1],alpha, NULL,-1);

    if(sol->run_integration_algorithm)
      err += m->param->integration_algorithm(m,ctx); // perform integration algorithm

    err += m->param->destroy_ctx(&ctx);

    if(err!=0)
      break;

    err += m->param->get_pF(m, pFnp1.data,2);
    err += m->param->get_pF(m,   pFn.data,1);
    err += m->param->get_pF(m, pFnm1.data,0);
    err += m->param->get_F(m,     Fn.data,1);
    err += m->param->get_F(m,   Fnm1.data,0);

    pJ = det(pFnp1);

    double dt_1_minus_alpha = -dts[DT_NP1]*(1.0-alpha)*pJ*hJ;
    err += residuals_el_constitutive_model_n_plus_alpha(f_npa,m,eid,ndofn,
                                                        pFnp1.data,pFn.data,Fnp1.data,Fn.data,
                                                        hFnp1.data,hFn.data,
                                                        is_it_couple_w_thermal,
                                                        alpha, dt_1_minus_alpha,fe,1);

    double dt_alpha = -dts[DT_N]*alpha*pJ*hJ;

    err += residuals_el_constitutive_model_n_plus_alpha(f_nm1pa,m,eid,ndofn,
                                                        pFn.data,pFnm1.data,Fn.data,Fnm1.data,
                                                        hFnp1.data,hFn.data,
                                                        is_it_couple_w_thermal,
                                                        alpha, dt_alpha,fe,0);
  }

  if(err==0)
  {
    for(int a=0; a<nne*nsd; a++)
      f[a] += f_npa[a] + f_nm1pa[a];
  }

  // free memory for thermal part

  free(u);
  free(f_npa);
  free(f_nm1pa);

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
int residuals_el_constitutive_model_1f(FEMLIB *fe,
                                       double *f,
                                       double *r_e,
                                       Grid *grid,
                                       MaterialProperty *mat,
                                       FieldVariables *fv,
                                       Solver *sol,
                                       LoadingSteps *load,
                                       CRPL *crpl,
                                       const PGFem3D_opt *opts,
                                       Multiphysics *mp,
                                       int mp_id,
                                       double dt)
{
  int err = 0;
  double alpha = -1.0; // if alpha < 0, no inertia
  int total_Lagrangian = 0;
  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;
  // @todo prevent warnings about unused variables, remove once it becomes used.
  (void)is_it_couple_w_chemical;

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

  double *u = PGFEM_malloc<double>(nne*nsd);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];
  }

  Tensor<2,3,double> Fr, Fnp1 = {},pFnp1,
            S = {}, M = {}, eFnM = {},eFnMT,
            hFn,hFnp1,hFnp1_I;

  Matrix<double> Tnp1, Tn, Tnm1;
  FieldVariables *fv_h = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    fv_h = fv->fvs[is_it_couple_w_thermal];
    Tnm1.initialization(fe->nne,1);
    Tnp1.initialization(fe->nne,1);
    Tn.initialization(fe->nne,1);

    // compute temperature for this element for each nodal point
    int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
    err += get_nodal_temperatures(fe, grid, fv_h, load, mp_cp_id,
                                  Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                  fv->subdivision_factor_np1,fv->subdivision_factor_n);
  }
  if(is_it_couple_w_chemical >=0)
  {}

  int compute_stiffness = 0;

  MATERIAL_ELASTICITY mat_e_new;
  MATERIAL_ELASTICITY *mat_e_in;

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double Jn = 1.0; // if upated Lagrangian, Jn = det(Fn), later updated
    double hJ = 1.0;
    double pJ = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    fe->update_deformation_gradient(ndofn,u,Fr.data);

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
                                       hFnp1.data,hFn.data,hFnm1);
      hJ = det(hFnp1);
      inv(hFnp1, hFnp1_I);
      if(total_Lagrangian)
        hFn = identity(i,j);
    }

    // --> update plasticity part
    if(total_Lagrangian)
    {
      if (sup->multi_scale) {
        cm_add_macro_F(sup,Fr.data);
      }

      // TOTAL LAGRANGIAN FORMULATION Fn = 1, Fnp1 = Fr
      Fnp1 = Fr(i,j);
    }
    else
    {
      if (sup->multi_scale) {
        PGFEM_printerr("Multi-scale formulation does not support UL!\n");
        PGFEM_Abort();
      }

      Tensor<2,3,double> Fn;
      err += m->param->get_F(m, Fn.data,1);
      Fnp1 = Fr(i,k)*Fn(k,j);  // compute Fnp1 = Fr*Fn
      Jn = det(Fn);
    }

    {
      /* check that deformation is invertible -> J > 0 */
      int terr = 0;
      getJacobian(Fnp1.data, eid, &terr);
      err += terr;
    }

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, m->param->type, Fnp1.data,dt,alpha, NULL,
                                                  hFn.data,hFnp1.data,-1);
    else
      err += construct_model_context(&ctx, m->param->type, Fnp1.data,dt,alpha, NULL,-1);      
    
    if(sol->run_integration_algorithm)
      err += m->param->integration_algorithm(m,ctx); // perform integration algorithm

    if(err>0)
        return err;

    err += func->get_pF(m,pFnp1.data,2);
    pJ = det(pFnp1);

    if(total_Lagrangian)
    {
      if(is_it_couple_w_thermal>=0)
      {
        Tensor<2,3,double> pFnp1_I;
        inv(pFnp1, pFnp1_I);
        M = hFnp1_I(i,k)*pFnp1_I(k,j);
      }
      else
        inv(pFnp1, M);

      eFnM = M(i,j);
      eFnMT(j,i) = M(i,j).to(j,i);

    }
    else
    {
      Tensor<2,3,double> eFn;
      if(is_it_couple_w_thermal>=0)
      {
        Tensor<2,3,double> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += func->get_pF(m,pFn.data,1);
        M = pFn(i,k)*hFn(k,l)*hFnp1_I(l,o)*pFnp1_I(o,j);

        Tensor<2,3,double> hFn_I;
        int stepno = 1; // 0 = time step = n-1
                        // 1 = time step = n
                        // 2 = time step = n+1
        inv(hFn, hFn_I);
        err += m->param->get_eF_of_hF(m,eFn.data,hFn_I.data,stepno);

      }
      else
      {
        Tensor<2,3,double> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += func->get_pF(m,pFn.data,1);
        M = pFn(i,k)*pFnp1_I(k,j);
        err += m->param->get_eF(m,eFn.data,1);
      }
      eFnM = eFn(i,k)*M(k,j);
      eFnMT(j,i) = eFnM(i,j).to(j,i);
    }

    // <-- update plasticity part

    // --> update elasticity part
    if((m->param)->uqcm)
    {
      double *x_ip = (fe->x_ip).m_pdata;
      ELASTICITY *elast = (m->param)->cm_elast;
      mat_e_in = elast->mat;
      err += material_properties_elasticity_at_ip(mat_e_in, &mat_e_new, x_ip[0], x_ip[1], x_ip[2]);
      elast->mat = &mat_e_new; // should be replaced by original mat_e_in after computation
      err += (m->param)->update_elasticity(m,ctx,NULL,S.data,compute_stiffness);
      elast->mat = mat_e_in;
    }
    else
      err += (m->param)->update_elasticity(m,ctx,NULL,S.data,compute_stiffness);
    // <-- update elasticity part
    err += m->param->destroy_ctx(&ctx);

    if(err!=0)
      break;


    Jn = Jn/pJ/hJ;
    err += compute_residual_vector(f,fe,Fr.data,eFnMT.data,eFnM.data,S.data,Jn);
  }

  free(u);

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
int residuals_el_constitutive_model_3f(FEMLIB *fe,
                                       double *f,
                                       double *r_e,
                                       Grid *grid,
                                       MaterialProperty *mat,
                                       FieldVariables *fv,
                                       Solver *sol,
                                       LoadingSteps *load,
                                       CRPL *crpl,
                                       const PGFem3D_opt *opts,
                                       Multiphysics *mp,
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
  int Pno   = fv->npres;
  int Vno   = fv->nVol;
  SUPP sup = load->sups[mp_id];

  Matrix<double> u(nne*nsd, 1);
  Matrix<double> dMdt(DIM_3x3*Vno,1);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u.m_pdata[a*nsd+b] = r_e[a*ndofn+b]; 
  }

  Tensor<2,3,double> Fr, Fnp1 = {},pFnp1,
            eSd = {}, M = {}, eFn = {},
            hFn,hFnp1,hFnp1_I;

  Tensor<4,3,double> Ld={};

  Matrix<double> Tnp1, Tn, Tnm1;
  FieldVariables *fv_h = NULL;

  if(is_it_couple_w_thermal >= 0)
  {
    fv_h = fv->fvs[is_it_couple_w_thermal];
    Tnm1.initialization(fe->nne,1);
    Tnp1.initialization(fe->nne,1);
    Tn.initialization(fe->nne,1);

    // compute temperature for this element for each nodal point
    int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
    err += get_nodal_temperatures(fe, grid, fv_h, load, mp_cp_id,
                                  Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                  fv->subdivision_factor_np1,fv->subdivision_factor_n);
  }
  if(is_it_couple_w_chemical >=0)
  {}

  Matrix<double> Ru(nne*nsd,1,0.0), Rp(Pno,1,0.0), Rt(Vno,1,0.0);
  Matrix<double> Kut(nne*nsd,Vno,0.0), Kup(nne*nsd,Pno,0.0), Ktp(Vno,Pno,0.0), Ktt(Vno,Vno,0.0),Kpt(Pno,Vno,0.0);

  Matrix<double> Nt(Vno,1), Np(Pno,1);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double Jn  = 1.0; // if upated Lagrangian, Jn = det(Fn), later updated
    double hJ  = 1.0;
    double pJ  = 1.0;
    double tJn = 1.0;

    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    fe->update_deformation_gradient(ndofn,u.m_pdata,Fr.data);

    fe->elem_shape_function(ip,Pno, Np.m_pdata);
    fe->elem_shape_function(ip,Vno, Nt.m_pdata);

    double theta_r = 0.0;
    double theta_n = 0.0;
    double Pnp1    = 0.0;

    for(int ia=1; ia<=Pno; ia++)
      Pnp1 += (fv->tf.P_np1(eid+1, ia) + fv->tf.dP(eid+1, ia))*Np(ia);

    for(int ia=1; ia<=Vno; ia++)
    {
      theta_r += (fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia))*Nt(ia);
      theta_n += fv->tf.V_n(  eid+1, ia)*Nt(ia);
    }

    Constitutive_model *m = &(fv->eps[eid].model[ip-1]);

    // get a shortened pointer for simplified CM function calls
    const Model_parameters *mp = m->param;
    ELASTICITY *elasticity = mp->cm_elast;

    // --> update deformations due to coupled physics
    if(is_it_couple_w_thermal >= 0)
    {
      // compute temperature at the integration point
      double T0 = fv_h->u0;
      double hFnm1[9];
      err += compute_temperature_at_ip(fe,grid,mat,T0,
                                       Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                       hFnp1.data,hFn.data,hFnm1);
      hJ = det(hFnp1);
      inv(hFnp1, hFnp1_I);
      if(total_Lagrangian)
        hFn = identity(i,j);
    }

    // --> update plasticity part
    if(total_Lagrangian)
    {
      if (sup->multi_scale) {
        cm_add_macro_F(sup,Fr.data);
      }

      // TOTAL LAGRANGIAN FORMULATION Fn = 1, Fnp1 = Fr
      Fnp1 = Fr(i,j);
    }
    else
    {
      if (sup->multi_scale) {
        PGFEM_printerr("Multi-scale formulation does not support UL!\n");
        PGFEM_Abort();
      }

      Tensor<2,3,double> Fn;
      err += mp->get_F(m, Fn.data,1);
      Fnp1 = Fr(i,k)*Fn(k,j);  // compute Fnp1 = Fr*Fn
      Jn = det(Fn);
    }

    {
      /* check that deformation is invertible -> J > 0 */
      int terr = 0;
      getJacobian(Fnp1.data, eid, &terr);
      err += terr;
    }

    void *ctx = NULL;
    if(is_it_couple_w_thermal>=0)
      err += construct_model_context_with_thermal(&ctx, mp->type, Fnp1.data,dt,alpha, NULL,
                                                  hFn.data,hFnp1.data,-1);
    else
      err += construct_model_context(&ctx, mp->type, Fnp1.data,dt,alpha, NULL,-1);      
    
    if(sol->run_integration_algorithm)
      err += mp->integration_algorithm(m,ctx); // perform integration algorithm

    if(err>0)
        return err;

    err += mp->get_pF(m,pFnp1.data,2);
    pJ = det(pFnp1);

    if(total_Lagrangian)
    {
      if(is_it_couple_w_thermal>=0)
      {
        Tensor<2,3,double> pFnp1_I;
        inv(pFnp1, pFnp1_I);
        M = hFnp1_I(i,k)*pFnp1_I(k,j);
      }
      else
        inv(pFnp1, M);

      eFn = identity(i,j);
    }
    else
    {
      if(is_it_couple_w_thermal>=0)
      {
        Tensor<2,3,double> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += mp->get_pF(m,pFn.data,1);
        M = pFn(i,k)*hFn(k,l)*hFnp1_I(l,o)*pFnp1_I(o,j);

        Tensor<2,3,double> hFn_I;
        int stepno = 1; // 0 = time step = n-1
                        // 1 = time step = n
                        // 2 = time step = n+1
        inv(hFn, hFn_I);
        err += mp->get_eF_of_hF(m,eFn.data,hFn_I.data,stepno);

      }
      else
      {
        Tensor<2,3,double> pFnp1_I, pFn;
        inv(pFnp1, pFnp1_I);
        err += mp->get_pF(m,pFn.data,1);
        M = pFn(i,k)*pFnp1_I(k,j);
        err += mp->get_eF(m,eFn.data,1);
      }
    }

    err += mp->compute_dMdt(m, ctx, fe->ST, Vno, dMdt.m_pdata);

    // <-- update plasticity part
    err += mp->compute_dev_stress(m,ctx,eSd.data);
    err += mp->compute_dev_tangent(m,ctx,Ld.data);
    double kappa = elasticity->mat->kappa;              
    double dUd_theta, d2Ud_theta2;
    
    double eJ = ttl::det(eFn);
    double MJ = ttl::det(M);
    double theta_e = theta_r*eJ*MJ;

    elasticity->compute_dudj(&dUd_theta, theta_e);
    elasticity->compute_d2udj2(&d2Ud_theta2, theta_e);
    dUd_theta = dUd_theta*kappa;
    d2Ud_theta2 = d2Ud_theta2*kappa;
    // <-- update elasticity part

    err += mp->destroy_ctx(&ctx);

    if(err!=0)
      break;

    tJn = Jn;
    Jn = Jn/pJ/hJ;

    Var_Carrier vc;

    vc.set_tenosrs(Fr.data, eFn.data, M.data, pFnp1.data, eSd.data, Ld.data);
    vc.set_scalars(theta_r, theta_n, tJn, Jn, Pnp1, dUd_theta, d2Ud_theta2);
      
    err += compute_Ru(fe, Ru.m_pdata, vc);
    err += compute_Rp(fe, Rp.m_pdata, Pno, Np.m_pdata, vc);
    err += compute_Rt(fe, Rt.m_pdata, Vno, Nt.m_pdata, vc);

    err += compute_Kut(fe, Kut.m_pdata, dMdt.m_pdata, vc, Vno, Nt.m_pdata);
    err += compute_Kup(fe, Kup.m_pdata, vc, Pno, Np.m_pdata);
    err += compute_Ktp(fe, Ktp.m_pdata, vc, Vno, Nt.m_pdata, Pno, Np.m_pdata);
    err += compute_Ktt(fe, Ktt.m_pdata, dMdt.m_pdata, vc, Vno, Nt.m_pdata);
  }

  Kpt.trans(Ktp);

  err += condense_F_3F_to_1F(f, nne, nsd, Pno, Vno,
                             Ru.m_pdata, Rt.m_pdata, Rp.m_pdata,
                             Kut.m_pdata, Kup.m_pdata, Ktp.m_pdata, Ktt.m_pdata, Kpt.m_pdata);    
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
                                    Grid *grid,
                                    MaterialProperty *mat,
                                    FieldVariables *fv,
                                    Solver *sol,
                                    LoadingSteps *load,
                                    CRPL *crpl,
                                    const PGFem3D_opt *opts,
                                    Multiphysics *mp,
                                    int mp_id,
                                    double dt)
{
  int err = 0;

  if(opts->analysis_type==CM)
    err += residuals_el_constitutive_model_1f(fe,f,r_e,grid,mat,fv,sol,load,crpl,opts,mp,mp_id,dt);

  if(opts->analysis_type==CM3F)
    err += residuals_el_constitutive_model_3f(fe,f,r_e,grid,mat,fv,sol,load,crpl,opts,mp,mp_id,dt);

  return err;
}


/// compute ouput variables e.g. effective stress and strain
///
/// Visit each element and compute output variables according to the element model type.
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv array of field variable object
/// \param[in] load object for loading
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \param[in] alpha mid point rule alpha
/// \return non-zero on internal error
int constitutive_model_update_NR(Grid *grid,
                                 MaterialProperty *mat,
                                 FieldVariables *fv,
                                 LoadingSteps *load,
                                 const PGFem3D_opt *opts,
                                 Multiphysics *mp,
                                 int mp_id,
                                 const double dt,
                                 double alpha)
{
  int err = 0;

  int total_Lagrangian = 1;
  if(opts->cm==UPDATED_LAGRANGIAN)
    total_Lagrangian = 0;

  int is_it_couple_w_thermal  = -1;
  int is_it_couple_w_chemical = -1;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_THERMAL)
      is_it_couple_w_thermal = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  EPS *eps = fv->eps;
  Node *node = grid->node;
  Element *elem = grid->element;

  int ndofn = fv->ndofn;
  int Pno   = fv->npres;
  int Vno   = fv->nVol;
  SUPP sup = load->sups[mp_id];


  for (int eid = 0; eid < grid->ne; eid++)
  {
    FEMLIB fe(eid,elem,node,0,total_Lagrangian);
    int nsd   = fe.nsd;
    int nne   = fe.nne;
    int ndofe = nne*ndofn;

    Matrix<long> cn(ndofe,1);
    long *nod = fe.node_id.m_pdata;

    Matrix<double> dr_e(ndofe,1), r_e(ndofe, 1), du(nne*nsd,1), u(nne*nsd, 1);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,grid->node,cn.m_pdata,mp_id);
        
    // get the deformation on the element
    def_elem_total(cn.m_pdata,ndofe,fv->u_np1,fv->d_u,elem,node,sup,r_e.m_pdata);
    def_elem(cn.m_pdata,ndofe,fv->dd_u,elem,node,dr_e.m_pdata,sup,2);

    for(int a=0;a<nne;a++)
    {
      for(int b=0; b<nsd;b++)
      {
        du.m_pdata[a*nsd+b] = dr_e.m_pdata[a*ndofn+b];
        u.m_pdata[a*nsd+b] = r_e.m_pdata[a*ndofn+b];
      }
    }

    Matrix<double> Ru(nne*nsd,1,0.0), Rp(Pno,1,0.0), Rt(Vno,1,0.0);
    Matrix<double> Ktu(Vno,nne*nsd,0.0), Kpu(nne*nsd,Pno,0.0), Ktp(Vno,Pno,0.0), Ktt(Vno,Vno,0.0),Kpt(Pno,Vno,0.0);
    Matrix<double> dMdu(DIM_3x3*nne*nsd,1);
    Matrix<double> dMdt(DIM_3x3*Vno,1);

    Matrix<double> Tnp1, Tn, Tnm1;
    FieldVariables *fv_h = NULL;

    if(is_it_couple_w_thermal >= 0)
    {
      fv_h = fv->fvs[is_it_couple_w_thermal];
      Tnm1.initialization(fe.nne,1);
      Tnp1.initialization(fe.nne,1);
      Tn.initialization(fe.nne,1);

      // compute temperature for this element for each nodal point
      int mp_cp_id = mp->coupled_ids[mp_id][is_it_couple_w_thermal+1];
      err += get_nodal_temperatures(&fe, grid, fv_h, load, mp_cp_id,
                                    Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,fv->subdivision_factor_np1,fv->subdivision_factor_n);
    }

    if(is_it_couple_w_chemical >= 0){}

    Tensor<2,3,double> tFr, Fr, Fnp1 = {},pFnp1,
              eSd = {}, M = {}, eFn = {},
              hFn,hFnp1,hFnp1_I;

    Tensor<4,3,double> Ld={};

    Matrix<double> Nt(Vno,1), Np(Pno,1);

    for(int ip=1; ip<=fe.nint; ip++)
    {
      double Jn  = 1.0; // if upated Lagrangian, Jn = det(Fn), later updated
      double hJ  = 1.0;
      double pJ  = 1.0;
      double tJn = 1.0;

      fe.elem_basis_V(ip);
      fe.elem_shape_function(ip,Pno, Np.m_pdata);
      fe.elem_shape_function(ip,Vno, Nt.m_pdata);      
      fe.update_shape_tensor();
      fe.update_deformation_gradient(ndofn,u.m_pdata,Fr.data);
      
      double theta_r = 0.0;
      double theta_n = 0.0;
      double Pnp1    = 0.0;

      for(int ia=1; ia<=Pno; ia++)
        Pnp1 += (fv->tf.P_np1(eid+1, ia) + fv->tf.dP(eid+1, ia))*Np(ia);

      for(int ia=1; ia<=Vno; ia++)
      {
        theta_r += (fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia))*Nt(ia);
        theta_n += fv->tf.V_n(  eid+1, ia)*Nt(ia);
      }

      Constitutive_model *m = &(eps[eid].model[ip-1]);

      // get a shortened pointer for simplified CM function calls
      const Model_parameters *mp = m->param;
      ELASTICITY *elasticity = mp->cm_elast;

      err += mp->get_F(m,Fnp1.data,  2);
      err += mp->get_pF(m,pFnp1.data,2);

      // --> update deformations due to coupled physics
      if(is_it_couple_w_thermal >= 0)
      {
        // compute temperature at the integration point
        double T0 = fv_h->u0;
        double hFnm1[9];
        err += compute_temperature_at_ip(&fe,grid,mat,T0,
                                         Tnp1.m_pdata,Tn.m_pdata,Tnm1.m_pdata,
                                         hFnp1.data,hFn.data,hFnm1);
        hJ = det(hFnp1);
        inv(hFnp1, hFnp1_I);
        if(total_Lagrangian)
          hFn = identity(i,j);
      }

      // --> update plasticity part
      if(total_Lagrangian)
      {
        // TOTAL LAGRANGIAN FORMULATION Fn = 1, Fnp1 = Fr
        Fnp1 = Fr(i,j);
      }
      else
      {
        Tensor<2,3,double> Fn;
        err += mp->get_F(m, Fn.data,1);
        Fnp1 = Fr(i,k)*Fn(k,j); // compute Fnp1 = Fr*Fn
        Jn = det(Fn);
      }

      void *ctx = NULL;
      if(is_it_couple_w_thermal>=0)
        err += construct_model_context_with_thermal(&ctx, mp->type, Fnp1.data,dt,alpha, NULL,
                                                    hFn.data,hFnp1.data,-1);
      else
        err += construct_model_context(&ctx, mp->type, Fnp1.data,dt,alpha, NULL,-1);
    
      err += mp->get_pF(m,pFnp1.data,2);
      pJ = det(pFnp1);

      if(total_Lagrangian)
      {
        if(is_it_couple_w_thermal>=0)
        {
          Tensor<2,3,double> pFnp1_I;
          inv(pFnp1, pFnp1_I);
          M = hFnp1_I(i,k)*pFnp1_I(k,j);
        }
        else
          inv(pFnp1, M);

        eFn = identity(i,j);
      }
      else
      {
        if(is_it_couple_w_thermal>=0)
        {
          Tensor<2,3,double> pFnp1_I, pFn;
          inv(pFnp1, pFnp1_I);
          err += mp->get_pF(m,pFn.data,1);
          M = pFn(i,k)*hFn(k,l)*hFnp1_I(l,o)*pFnp1_I(o,j);

          Tensor<2,3,double> hFn_I;
          int stepno = 1; // 0 = time step = n-1
                          // 1 = time step = n
                          // 2 = time step = n+1
          inv(hFn, hFn_I);
          err += mp->get_eF_of_hF(m,eFn.data,hFn_I.data,stepno);

        }
        else
        {
          Tensor<2,3,double> pFnp1_I, pFn;
          inv(pFnp1, pFnp1_I);
          err += mp->get_pF(m,pFn.data,1);
          M = pFn(i,k)*pFnp1_I(k,j);
          err += mp->get_eF(m,eFn.data,1);
        }
      }

      err += mp->compute_dMdu(m, ctx, fe.ST, nne, ndofn, dMdu.m_pdata);
      err += mp->compute_dMdt(m, ctx, fe.ST, Vno, dMdt.m_pdata);

      // <-- update plasticity part
      err += mp->compute_dev_stress(m,ctx,eSd.data);
      err += mp->compute_dev_tangent(m,ctx,Ld.data);
      double kappa = elasticity->mat->kappa;              
      double dUd_theta, d2Ud_theta2;
      
      double eJ = ttl::det(eFn);
      double MJ = ttl::det(M);
      double theta_e = theta_r*eJ*MJ;
  
      elasticity->compute_dudj(&dUd_theta, theta_e);
      elasticity->compute_d2udj2(&d2Ud_theta2, theta_e);
      dUd_theta = dUd_theta*kappa;
      d2Ud_theta2 = d2Ud_theta2*kappa;
      // <-- update elasticity part

      err += (m->param)->destroy_ctx(&ctx);

      tJn = Jn;
      Jn = Jn/pJ/hJ;

      Var_Carrier vc;

      vc.set_tenosrs(Fr.data, eFn.data, M.data, pFnp1.data, eSd.data, Ld.data);
      vc.set_scalars(theta_r, theta_n, tJn, Jn, Pnp1, dUd_theta, d2Ud_theta2);

      err += compute_Ru(&fe, Ru.m_pdata, vc);
      err += compute_Rp(&fe, Rp.m_pdata, Pno, Np.m_pdata, vc);
      err += compute_Rt(&fe, Rt.m_pdata, Vno, Nt.m_pdata, vc);

      err += compute_Ktu(&fe, Ktu.m_pdata, dMdu.m_pdata, vc, Vno, Nt.m_pdata);            
      err += compute_Kup(&fe, Kpu.m_pdata, vc, Pno, Np.m_pdata);
      err += compute_Ktp(&fe, Ktp.m_pdata, vc, Vno, Nt.m_pdata, Pno, Np.m_pdata);
      err += compute_Ktt(&fe, Ktt.m_pdata, dMdt.m_pdata, vc, Vno, Nt.m_pdata);      
    }

    Kpu.trans();
    Kpt.trans(Ktp);  

    Matrix<double> d_theta(Vno, 1), dP(Pno, 1);
    err += compute_d_theta_dP(d_theta.m_pdata, dP.m_pdata, du.m_pdata,
                              nne, nsd, Pno, Vno,
                              Ru.m_pdata, Rt.m_pdata, Rp.m_pdata,
                              Kpu.m_pdata, Ktu.m_pdata, Ktp.m_pdata, Ktt.m_pdata, Kpt.m_pdata);

    for(int ia=1; ia<=Pno; ia++)
      fv->tf.ddP(eid+1,ia) = dP(ia);


    for(int ia=1; ia<=Vno; ia++)
      fv->tf.ddV(eid+1,ia) = d_theta(ia);      
  }

  return err;
}

int cm_write_tensor_restart(FILE *fp, const double *tensor)
{
  int err = 0;
  fprintf(fp, "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
          tensor[0], tensor[1], tensor[2],
          tensor[3], tensor[4], tensor[5],
          tensor[6], tensor[7], tensor[8]);
  return err;
}

int cm_read_tensor_restart(FILE *fp, double *tensor)
{
  int err = 0;
  fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
         &tensor[0], &tensor[1], &tensor[2],
         &tensor[3], &tensor[4], &tensor[5],
         &tensor[6], &tensor[7], &tensor[8]);
  return err;
}
