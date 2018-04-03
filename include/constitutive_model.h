/**
 * Define the interface to generalized constitutive models with
 * their associated data structure(s). The user is responsible for
 * defining and linking the functions associated with a particular
 * model. NEED MORE USAGE DETAILS.
 *
 * Authors:
 *  Matt Mosby, [1], <mmosby1@nd.edu>
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  Alberto Salvadori, [1], <asalvad2@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */
#ifndef PGFEM3D_CONSTITUTIVE_MODEL_H
#define PGFEM3D_CONSTITUTIVE_MODEL_H

#include "PGFem3D_data_structure.h"
#include "PGFem3D_options.h"
#include "crpl.h"
#include "sig.h"
#include "state_variables.h"
#include "supp.h"
#include <stdlib.h>
#include <stdio.h>
#include "material_properties.h"
#include "hyperelasticity.h"

struct EPS;
struct Element;
struct Node;

/// Pre-declare HOMMAT structure
struct HOMMAT;
#ifndef TYPE_HOMMAT
#define TYPE_HOMMAT
typedef struct HOMMAT FEMLIB;
#endif

#ifndef PGFEM3D_DEV_TEST
#define PGFEM3D_DEV_TEST 1
#endif

/**
 * Enumeration for the model type
 */
enum model_type {
  HYPER_ELASTICITY,
  CRYSTAL_PLASTICITY,
  BPA_PLASTICITY,
  ISO_VISCOUS_DAMAGE,
  J2_PLASTICITY_DAMAGE,
  POROVISCO_PLASTICITY,  
  NUM_MODELS,
  ISO_VISCOUS_SPLIT_DAMAGE = 97,
  J2_PLASTICITY_SPLIT_DAMAGE = 98,   
  TESTING=99,
  CM_UQCM=100
};

/// Enumeration for the frame
enum integration_frame {
  UPDATED_LAGRANGIAN,
  TOTAL_LAGRANGIAN,
  MIXED_ANALYSIS_MODE
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

/// Pre-declare the Model_parameters structure
class Model_parameters
{
 public:
  /** Pointer to isotropic material props */
  double *pF;
  const HOMMAT *p_hmat;
  int mat_id; // Global material id, mat_id may not be the same as the hommat id
  int uqcm;   // UQ study through constitutive model 0: no, or yes
  bool cm3f;

  MATERIAL_CONSTITUTIVE_MODEL *cm_mat;
  ELASTICITY *cm_elast;
  
  // set members variables to initial values (zeros and NULLs)   
  void set_nulls(void)
  {
    p_hmat            = NULL;
    mat_id            = -1;
    uqcm              = 0; 
    cm3f              = false;
    cm_mat            = NULL;
    cm_elast          = NULL;    
    type              = -1;
    n_param           = 0;
    model_param       = NULL;
    n_param_index     = 0;
    model_param_index = NULL;
    pF                = NULL;    
  };

  /// Construct a Model_parameters object. The object is left in an
  /// invalid state until model_parameters_initialize is called.
  Model_parameters()
  {
    set_nulls();
  };

  /// Destroy a Model_parameters object.
  virtual ~Model_parameters()
  {
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
        return;
    }    
    model_dependent_finalization();
    finalization();
  };


  /// Initialize the Model_parameters object. The object may be used
  /// after calling this function. Calling this function on an already
  /// initialized object is undefined.
  ///
  /// \param[in] p_hmat material property object
  /// \param[in] type   constitutive model type
  /// \return non-zero on error.
  int initialization(const HOMMAT *p_hmat,
                     const size_t type);

  int finalization(void);
  virtual int model_dependent_initialization(void)
  { return 0; };

  virtual int model_dependent_finalization(void)
  { return 0; };

  /// User defined function for the Constitutive_model integration
  /// algorithm. This function shall be implemented such that it modifies
  /// the internal state to contain the updated values upon exit, i.e.,
  /// subsequent calls to get_F functions will return the correct
  /// values without re-integrating the constitutive model.
  ///
  /// \param[in,out] m,   pointer to a Constitutive_model object.
  /// \param[in,out] ctx, handle to a user defined structure that
  ///                     controls execution of the user-defined integration algorithm. This
  ///                     is the mechanism for passing model-specific information into the
  ///                     general function interface.
  /// \param[out] EXA_metric, exascale metric counter for total number of integration iterations
  /// \return non-zero on internal error that should be handled by the calling function.
  virtual int integration_algorithm(Constitutive_model *m,
                                    const void *usr_ctx,
                                    int &EXA_metric) const { return 0; };

  /// User defined function to compute deviatroic stress tenosr.
  ///
  /// \param[in]     m,   pointer to Constitutive_model object.
  /// \param[in,out] ctx, handle to a user defined structure that
  ///                     controls execution of the user-defined function. This is the
  ///                     mechanism for passing model-specific information into the general
  ///                     function interface.
  /// \param[out]    S,   computed deviatroic 2nd order stress tenosr.
  ///                     that is populated with values from the constitutive model.
  /// \return non-zero value on internal error that should be handled by the calling function.
  virtual int compute_dev_stress(const Constitutive_model *m,
                                 const void *ctx,
                                 double *S) const { return 0; };

  /// User defined function to compute deviatroic stress tenosr.
  ///
  /// \param[in]     m,   pointer to Constitutive_model object.
  /// \param[in]     eC,  eC = eF'*eF
  /// \param[out]    S,   computed deviatroic 2nd order stress tenosr.
  ///                     that is populated with values from the constitutive model.
  /// \return non-zero value on internal error that should be handled by the calling function.
  virtual int compute_dev_stress(const Constitutive_model *m,
                                 double *eC,
                                 double *S) 
  const 
  {    
    cm_elast->compute_PK2_dev(eC, cm_elast->mat, S); 
    return 0; 
  };
                                 
  /// User defined function to compute volumetric stress contributions
  ///
  /// \param[in]     m,     pointer to Constitutive_model object.
  /// \param[in,out] ctx,   handle to a user defined structure that
  ///                       controls execution of the user-defined function. This is the
  ///                       mechanism for passing model-specific information into the general
  ///                       function interface.
  /// \param[outƒ]    value, computed value (passed by reference)
  /// \return non-zero value on internal error that should be handled by the calling function.
  virtual int compute_dudj(const Constitutive_model *m,
                           const void *ctx,
                           double *value) const { return 0; };

  /// User defined function to compute volumetric stress contributions
  ///
  /// \param[in]     m,       pointer to Constitutive_model object.
  /// \param[in]     theta_e, theta_e = J = det(eF)
  /// \return computed volumetric stress contribution
  virtual double compute_dudj(const Constitutive_model *m,
                              const double theta_e,
                              const int npa,
                              const double alpha) 
  const 
  {
    double dudj = 0.0;
    cm_elast->compute_dudj(&dudj, theta_e);  
    return dudj*(cm_elast->mat->kappa);
  };

  ///User defined function to compute deviatroic stiffness tenosr.
  ///
  /// \param[in]     m,   pointer to Constitutive_model object.
  /// \param[in,out] ctx, handle to a user defined structure that
  ///                     controls execution of the user-defined function. This is the
  ///                     mechanism for passing model-specific information into the general
  ///                     function interface.
  /// \param[out]    L,   computed deviatroic 4th order stiffness tenosr.
  ///                     that is populated with values from the constitutive model.
  /// \return non-zero value on internal error that should be handled by the calling function.
  virtual int compute_dev_tangent(const Constitutive_model *m,
                                  const void *ctx,
                                  double *L) const { return 0; };

  ///User defined function to compute deviatroic stiffness tenosr.
  ///
  /// \param[in]     m,   pointer to Constitutive_model object.
  /// \param[in]     eC,  eC = eF'*eF
  /// \param[out]    L,   computed deviatroic 4th order stiffness tenosr.
  ///                     that is populated with values from the constitutive model.
  /// \return non-zero value on internal error that should be handled by the calling function.
  virtual int compute_dev_tangent(const Constitutive_model *m,
                                  double *eC,
                                  double *L)
  const
  { 
    cm_elast->compute_tangent_dev(eC, cm_elast->mat, L);
    return 0; 
  };
                                  
  /// User defined function to compute volumetric elasticity contributions
  ///
  /// \param[in]     m,     pointer to Constitutive_model object.
  /// \param[in,out] ctx,   handle to a user defined structure that
  ///                       controls execution of the user-defined function. This is the
  ///                       mechanism for passing model-specific information into the general
  ///                       function interface.
  /// \param[out]    value, computed value (passed by reference)
  /// \return non-zero value on internal error that should be handled by the calling function.
  virtual int compute_d2udj2(const Constitutive_model *m,
                             const void *ctx,
                             double *value) const { return 0; };

  /// User defined function to compute volumetric elasticity contributions
  ///
  /// \param[in]     m,     pointer to Constitutive_model object.
  /// \param[in]     theta_e, theta_e = J = det(eF)
  /// \return computed volumetric elasticity contribution
  virtual double compute_d2udj2(const Constitutive_model *m,
                                double theta_e,
                                const int npa,
                                const double alpha)
  const 
  {
    double ddudj = 0.0;
    this->cm_elast->compute_d2udj2(&ddudj, theta_e);  
    return ddudj*(cm_elast->mat->kappa);    
  };
                             
  /// User defined function to compute the elastic algorithmic stiffness tangent
  ///
  /// \param[in]  m,                 pointer to Constitutive_model object.
  /// \param[in]  ctx,               a context for the model
  /// \param[out] L,                 4th order elasticity tensor
  /// \param[out] S,                 2nd order PKII tensor
  /// \param[in]  compute_stiffness, if 1 compute elasticity tensor
  ///                                   0 no compute elasticity tensor
  /// \return non-zero on internal error that should be handled by the
  virtual int update_elasticity(const Constitutive_model *m,
                                const void *ctx,
                                double *L,
                                double *S,
                                const int compute_stiffness) const { return 0; };
                                
  /// User defined function to compute the deviatroic part of elastic stiffness tangent
  ///
  /// \param[in]  m,                 pointer to Constitutive_model object.
  /// \param[in]  eF,                elastic deformation gradient
  /// \param[out] L,                 4th order elasticity tensor
  /// \param[out] S,                 2nd order PKII tensor
  /// \param[in]  compute_stiffness, if 1 compute elasticity tensor
  ///                                   0 no compute elasticity tensor
  /// \return non-zero on internal error that should be handled by the
  virtual int update_elasticity_dev(const Constitutive_model *m,
                                    double *eF,
                                    double *L,
                                    double *S,
                                    const int npa,
                                    const double alpha,
                                    const int compute_stiffness = 0) 
  const
  {
    int err = 0;
    double eC[9];
    for(int ia=0; ia<3; ia++){
      for(int ib=0; ib<3; ib++){
         eC[ia*3 + ib] = 0.0;
         for(int ic=0; ic<3; ic++)
           eC[ia*3 + ib] += eF[ic*3 + ia]*eF[ic*3 + ib];
       }
     }
     err += compute_dev_stress(m,eC,S);
     if(compute_stiffness)
       err += compute_dev_tangent(m,eC,L);
             
     return err; 
  };

  /// User defined function that to update the internally
  /// advance: (n + 1) -> (n) the internal state variables.
  ///
  /// \param[in,out] m, pointer to Constitutive_model object.
  /// \return non-zero on internal error that should be handled by the
  /// calling function.
  virtual int update_state_vars(Constitutive_model *m) const { return 0; };

  /// User defined function that to reset the internally
  /// reset: (n + 1) <- (n) the internal state variables.
  ///
  /// \param[in,out] m, pointer to Constitutive_model object.
  /// \return non-zero on internal error that should be handled by the
  /// calling function.
  virtual int reset_state_vars(Constitutive_model *m) const { return 0; };

  /// User defined function to reset state variable from temporal space
  ///
  /// \param[in]     m,   constant reference to a Constitiutive_model object.
  /// \param[in,out] var, contains state variables to store or get
  /// \return non-zero on internal error
  virtual int reset_state_vars_using_temporal(const Constitutive_model *m,
                                              State_variables *var) const { return 0; };

  /// User defined function to upated state variable at n+1 to temporal space
  ///
  /// \param[in]     m,   constant reference to a Constitiutive_model object.
  /// \param[in,out] var, contains state variables to store or get
  /// \return non-zero on internal error
  virtual int update_np1_state_vars_to_temporal(const Constitutive_model *m,
                                                State_variables *var) const { return 0; };

  /// User defined function to store state variable to temporal space
  ///
  /// \param[in]     m,   constant reference to a Constitiutive_model object.
  /// \param[in,out] var, contains state variables to store or get
  /// \return non-zero on internal error
  virtual int save_state_vars_to_temporal(const Constitutive_model *m,
                                          State_variables *var) const { return 0; };

  /// A user described function that allocates a Model_var_info object
  /// and populates the internal structure. This function should fully
  /// allocate the internals of Model_var_info and *copy* data rather
  /// than hold references.
  /// \param[in, out] info, reference to object containing model variable info
  /// \return non-zero on error.
  virtual int get_var_info(Model_var_info &info) const { return 0; };

  /// User defined function to return the total deformation.
  ///
  /// \param[in]  m,      constant reference to a Constitiutive_model object.
  /// \param[out] F,      reference to Matrix object that contains the
  ///                     deformation gradient upon exit.
  /// \param[in]  stepno, time step id 0: n-1
  ///                                 1: n
  ///                                 2: n+1
  /// \return non-zero on internal error
  virtual int get_F(const Constitutive_model *m,
                    double *F,
                    const int stepno) const { return 0; };

  /// User defined function to return the plastic part deformation.
  ///
  /// \param[in]  m,      constant reference to a Constitiutive_model object.
  /// \param[out] F,      reference to Matrix object that contains the
  ///                     deformation gradient upon exit.
  /// \param[in]  stepno, time step id 0: n-1
  ///                                 1: n
  ///                                 2: n+1
  /// \return non-zero on internal error
  virtual int get_pF(const Constitutive_model *m,
                     double *F,
                     const int stepno) const { return 0; };

  /// User defined function to return the elastic deformation.
  ///
  /// \param[in]  m,      constant reference to a Constitiutive_model object.
  /// \param[out] F,      reference to Matrix object that contains the
  ///                     deformation gradient upon exit.
  /// \param[in]  stepno, time step id 0: n-1
  ///                                 1: n
  ///                                 2: n+1
  /// \return non-zero on internal error
  virtual int get_eF(const Constitutive_model *m,
                     double *F,
                     const int stepno) const { return 0; };


  /// User defined function to return the elastic deformation gradient with thermal
  /// expansitions
  ///
  /// \param[in]  m,      constant reference to a Constitiutive_model object.
  /// \param[out] F,      reference to Matrix object that contains the
  ///                     deformation gradient upon exit.
  /// \param[in]  hFI,    inverse of thermal part of the deformation gradient
  /// \param[in]  stepno, time step id 0: n-1
  ///                                 1: n
  ///                                 2: n+1
  /// \return non-zero on internal error
  virtual int get_eF_of_hF(const Constitutive_model *m,
                           double *F,
                           double *hFI,
                           const int stepno) const { return 0; };


  /// User defined function to return a hardening variables.
  ///
  /// \param[in]  m,      constant reference to a Constitiutive_model object.
  /// \param[out] var,    contains the value of the state variable upon exit.
  /// \param[in]  stepno, time step id 0: n-1
  ///                                 1: n
  ///                                 2: n+1
  /// \return non-zero on internal error
  virtual int get_hardening(const Constitutive_model *m,
                            double *var,
                            const int stepno)
  const
  {
    *var = 0.0;
    return 0;
  };

  virtual int get_plast_strain_var(const Constitutive_model *m,
                                   double *lam_p)
  const
  {
    *lam_p = 0.0;
    return 0;
  }

  /// objtaion subdivision parameters
  /// It doese not modify the internal state variables.
  ///
  /// \param[in]  m,   constant reference to a Constitiutive_model object.
  /// \param[in]  t,   constant t (time).
  /// \param[out] var, contains the value of the state variable upon exit.
  /// \return non-zero on internal error
  virtual int get_subdiv_param(const Constitutive_model *m,
                               double *var,
                               const double t)
  const
  {
    *var = 0.0;
    return 0;
  };

  /// A user described function that writes restart file at
  /// integration point
  ///
  /// \param[in] fp, file pointer for writing restart file
  /// \param[in] m,  Constitutive model object
  /// \return non-zero on error.
  virtual int write_restart(FILE *fp,
                            const Constitutive_model *m) const { return 0; };

  /// A user described function that reads restart file at
  /// integration point
  ///
  /// \param[in] fp, file pointer for reading restart file
  /// \param[in] m,  Constitutive model object
  /// \return non-zero on error.
  virtual int read_restart(FILE *fp,
                           Constitutive_model *m) const { return 0; };

  /// User defined function to destroy a context for the model. This
  /// function shall destroy any internally allocated data and
  /// invalidate the handle, i.e., *ctx = NULL.
  ///
  /// \param[in] ctx a context for the model
  /// \return non-zero on error.
  virtual int destroy_ctx(void **ctx) const { return 0; };

  /// User defined function to compute the linearization of the plastic
  /// deformation w.r.t. the displacement variable.
  ///
  /// \param[in]  m,       Constitutive model object
  /// \param[in]  ctx,     model specific user context
  /// \param[in]  Grad_op, 4-index FE gradient operator where the 1st two indices are the node followed by the DOF.
  /// \param[in]  nne,     number of nodes on the element (length of idx 1 in Grad_op)
  /// \param[in]  ndofn,   number of dofs on each node (length of idx 2 in Grad_op)
  /// \param[out] dM_du,   4-index FE linearization of M, same dimensions as Grad_op
  /// \return non-zero on internal error
  virtual int compute_dMdu(const Constitutive_model *m,
                           const void *ctx,
                           const double *Grad_op,
                           const int nne,
                           const int ndofn,
                           double *dM_du)
  const
  {
    // there is no plastic deformation in this formulation, return zeros
    // in dM_du
    memset(dM_du, 0, nne*ndofn*9*sizeof(double));
    return 0;
  };


  /// User defined function to compute the linearization of the plastic
  /// deformation w.r.t. the volume variable.
  ///
  /// \param[in]  m,       Constitutive model object
  /// \param[in]  ctx,     model specific user context
  /// \param[in]  Grad_op, 4-index FE gradient operator where the 1st two indices are the node followed by the DOF.
  /// \param[in]  nne,     number of nodes on the element (length of idx 1 in Grad_op)
  /// \param[in]  ndofn,   number of dofs on each node (length of idx 2 in Grad_op)
  /// \param[out] dM_du,   4-index FE linearization of M, same dimensions as Grad_op
  /// \return non-zero on internal error
  virtual int compute_dMdt(const Constitutive_model *m,
                           const void *ctx,
                           const double *Grad_op,
                           const int Vno,
                           double *dM_dt)
  const
  {
    // there is no plastic deformation in this formulation, return zeros
    // in dM_du
    memset(dM_dt, 0, Vno*9*sizeof(double));
    return 0;
  };

  /// User defined function to set the initial values of the state
  /// variables for the particular model.
  ///
  /// \param[in] m, CM object with internal data set from the buffer
  /// \return non-zero on error.
  virtual int set_init_vals(Constitutive_model *m)
  const
  { return 0;}


  /// User defined function for reading in material parameters from a file.
  ///
  /// \param[in] in, file pointer for reading material parameters
  /// \return non-zero on error.
  virtual int read_param(FILE *in) const { return 0; };

  /*
 /// User defined function that returns the size of the data to be
 /// packed/unpacked.
 /// Does not modify the CM object or any of the data it holds.
 ///
 /// \param[in] m, CM object with internal data set from the buffer
 /// \return size in bytes of the pack/unpack data
 int get_size(const Constitutive_model *m);


 /// User defined function to pack the CM data into a buffer (see pack_data).
 /// Does not modify the CM object or any of the data it holds.
 ///
 /// \param[in,out] buffer, a buffer to insert data to
 ///
 /// \param[in,out] pos,    insert position in the buffer. Upon exit - next
 ///                        insertion position.
 /// \return non-zero on error.
 int pack(const Constitutive_model *m,
 char *buffer,
 size_t *pos);

 /// User defined function to unpack CM data from a buffer (see also
 /// usr_pack, unpack_data).
 ///
 /// \param[out]    m,      CM object with internal data set from the buffer
 /// \param[in]     buffer, the buffer to read data from
 /// \param[in,out] pos,    the position in buffer to begin reading from.
 ///                        Upon exit - position for next read.
 /// \return        non-zero on error.
 int unpack(Constitutive_model *m,
 const char *buffer,
 size_t *pos);
  */
 public:
  size_t type;          /// Model type, see enumeration @model_type
  size_t n_param;       /// array for storing the model constants.

  double *model_param;
  size_t n_param_index; /// array for storing the model integer type constant
  int *model_param_index;
};

/// Allocate a Model_parameters
int construct_Model_parameters(Model_parameters **p, int model_id, int model_type);

/// Allocate and populate a list of Model_parameters given the number
/// of materials.
int read_model_parameters_list(const int n_mat,
                               HOMMAT *hmat_list,
                               FILE *in);

/// General interface to a constitutive model.
///
/// Holds a Model_parameters object.
/// Has a State_variables object.
class Constitutive_model
{
 public:
  const Model_parameters *param;
  int model_id;
  State_variables **vars_list;

  /// Construct a Constitutive_model object. The object is left in an
  /// invalid state until constitutive_model initialize is called.
  Constitutive_model()
  {
    param = NULL;
  }

  /// destructor a Constitutive_model object.
  ~Constitutive_model()
  {
    // drop pointer to Model_parameters object
    this->param = NULL;
    // drop pointer to state variables
    this->vars_list = NULL;
  }

  /// Initialize the Constitutive_model object given the material type
  /// and properties. The object may be used after calling this function.
  ///
  /// \return non-zero on error.
  int initialization(const Model_parameters *param);

  /// User defined function that returns the size of the data to be
  /// packed/unpacked.
  /// Does not modify the CM object or any of the data it holds.
  ///
  /// \return size in bytes of the pack/unpack data
  int get_size(void);


  /// User defined function to pack the CM data into a buffer (see pack_data).
  /// Does not modify the CM object or any of the data it holds.
  ///
  /// \param[in,out] buffer, a buffer to insert data to
  ///
  /// \param[in,out] pos,    insert position in the buffer. Upon exit - next
  ///                        insertion position.
  /// \return non-zero on error.
  int pack(char *buffer,
           size_t *pos);

  /// User defined function to unpack CM data from a buffer (see also
  /// usr_pack, unpack_data).
  ///
  /// \param[in]     buffer, the buffer to read data from
  /// \param[in,out] pos,    the position in buffer to begin reading from.
  ///                        Upon exit - position for next read.
  /// \return        non-zero on error.
  int unpack(const char *buffer,
             size_t *pos);

  /// This function is running integration algorithm such that it modifies
  /// the internal state to contain the updated values upon exit, i.e.,
  /// subsequent calls to get_F functions will return the correct
  /// values without re-integrating the constitutive model.
  ///
  /// \param[in] *Fnp1  deformation gradient at t(n+1)
  /// \param[in] *hFn   thermal part of deformation gradient at t(n)
  /// \param[in] *hFnp1 thermal part of deformation gradient at t(n+1)
  /// \param[in] dt                     time step size
  /// \param[in] alpha                  mid point rule alpha
  /// \param[out] EXA_metric exascale metric counter for total number of integration iterations
  /// \param[in] is_it_couple_w_thermal checking coupling with thermal
  ///                                   if true: apply thermal expansitions
  /// \param[in] tf_factor              (theta/J)^(1/3) used for computing true Fnp1  
  /// \return non-zero on internal error that should be handled by the calling function.
  int run_integration_algorithm(double *Fnp1, 
                                double *hFn,
                                double *hFnp1,
                                double dt,
                                double alpha,
                                int &EXA_metric,
                                int is_it_couple_w_thermal = 0,
                                double tf_factor = 1.0);
};

///
/// Initialize the constitutive model object at each integration point.
///
/// \return non-zero on error.
///
int init_all_constitutive_model(EPS *eps,
                                const int ne,
                                const Element *elem,
                                const int n_mat,
                                const HOMMAT *hmat_list);
/// save state variables
///
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \return non-zero on internal error
int constitutive_model_save_state_vars_to_temporal(FieldVariables *fv,
                                                   Grid *grid);
/// update state variables
///
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \return non-zero on internal error
int constitutive_model_update_np1_state_vars_to_temporal(FieldVariables *fv,
                                                         Grid *grid);

/// reset state variables using priori stored values
///
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \return non-zero on internal error
int constitutive_model_reset_state_using_temporal(FieldVariables *fv,
                                                  Grid *grid);

/// Reset the cinstitutive model at each integration point in the
/// domain by calling its respective reset function.
///
/// \param[in, out] eps,  structure of strains
/// \param[in]      ne,   number of elements
/// \param[in]      elem, element object
/// \return non-zero on error.
int constitutive_model_reset_state(EPS *eps,
                                   const int ne,
                                   const Element *elem);


/// Compute the elastic stress and tangent tensors. Note that the total
/// stress and tangent (not only the deviatoric parts) are returned.
///
/// \param[in] F, *total* deformation gradient
/// \param[in] compute_stiffness, flag to comptue the stiffness tensor
/// (L) if non-zero. Otherwise L is unchanged.
///
/// \return non-zero on internal error.
///
int constitutive_model_default_update_elasticity(const Constitutive_model *m,
                                                 const double *eF,
                                                 double *L,
                                                 double *S,
                                                 const int compute_stiffness);

/// update values for next time step: variables[tn] = variables[tn+1]
/// \return non-zero on error.
int constitutive_model_update_time_steps(const Element *elem,
                                         Node *node,
                                         EPS *eps,
                                         const int ne,
                                         const int nn,
                                         const int ndofn,
                                         const double* r,
                                         const double dt,
                                         const int total_Lagrangian,
                                         const int mp_id);

int constitutive_model_test(const HOMMAT *hmat,
                            double *L_in,
                            int Print_results);

/// compute ouput variables e.g. effective stress and strain
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
                                               const Multiphysics& mp,
                                               int mp_id,
                                               const double dt,
                                               double alpha);

/// Compute the physics-based subdivision paramter for all integration
/// points on the domain.
int cm_get_subdivision_parameter(double *subdiv_param,
                                 const int ne,
                                 const Element *elem,
                                 const EPS *eps,
                                 const double dt);

/// Construct the model context for any model.
int construct_model_context(void **ctx,
                            const int type,
                            double *F,
                            const double dt,
                            const double alpha,
                            double *eFnpa,
                            const int npa);

struct FEMLIB;

/// compute element stiffness matrix in transient
///
/// \param[in] fe finite element helper object
/// \param[out] lk computed element stiffness matrix
/// \param[in] re_np1 nodal variables at t(n+1) in the current element
/// \param[in] re_npa nodal variables at (1-alpha)r(n) + alpha*r(n+1) in the current element
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
/// \param[out] EXA_metric exascale metric counter for total number of integration iterations
/// \return non-zero on internal error
int stiffness_el_constitutive_model_w_inertia(FEMLIB *fe,
                                              double *lk,
                                              double *re_np1,
                                              double *re_npa,
                                              Grid *grid,
                                              MaterialProperty *mat,
                                              FieldVariables *fv,
                                              pgfem3d::Solver *sol,
                                              LoadingSteps *load,
                                              CRPL *crpl,
                                              const PGFem3D_opt *opts,
                                              const Multiphysics& mp,
                                              int mp_id,
                                              double dt,
                                              int &EXA_metric);

/// compute element stiffness matrix in quasi steady state
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
/// \param[out] EXA_metric exascale metric counter for total number of integration iterations
/// \return non-zero on internal error
int stiffness_el_constitutive_model(FEMLIB *fe,
                                    double *lk,
                                    double *r_e,
                                    Grid *grid,
                                    MaterialProperty *mat,
                                    FieldVariables *fv,
                                    pgfem3d::Solver *sol,
                                    LoadingSteps *load,
                                    CRPL *crpl,
                                    const PGFem3D_opt *opts,
                                    const Multiphysics& mp,
                                    int mp_id,
                                    double dt,
                                    int &EXA_metric);

/// compute element residual vector in transient
///
/// \param[in] fe finite element helper object
/// \param[out] f computed element residual vector
/// \param[in] re_np1 nodal variables at t(n+1) in the current element
/// \param[in] re_npa nodal variables at (1-alpha)r(n)   + alpha*r(n+1) in the current element
/// \param[in] re_nma nodal variables at (1-alpha)r(n-1) + alpha*r(n)   in the current element
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
/// \param[out] EXA_metric exascale metric counter for total number of integration iterations
/// \return non-zero on internal error
int residuals_el_constitutive_model_w_inertia(FEMLIB *fe,
                                              double *f,
                                              double *re_np1,
                                              double *re_npa,
                                              double *re_nma,
                                              Grid *grid,
                                              MaterialProperty *mat,
                                              FieldVariables *fv,
                                              pgfem3d::Solver *sol,
                                              LoadingSteps *load,
                                              CRPL *crpl,
                                              const PGFem3D_opt *opts,
                                              const Multiphysics& mp,
                                              const double *dts,
                                              int mp_id,
                                              double dt,
                                              int &EXA_metric);

/// compute element residual vector in quasi steady state
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
/// \param[out] EXA_metric exascale metric counter for total number of integration iterations
/// \return non-zero on internal error
int residuals_el_constitutive_model(FEMLIB *fe,
                                    double *f,
                                    double *r_e,
                                    Grid *grid,
                                    MaterialProperty *mat,
                                    FieldVariables *fv,
                                    pgfem3d::Solver *sol,
                                    LoadingSteps *load,
                                    CRPL *crpl,
                                    const PGFem3D_opt *opts,
                                    const Multiphysics& mp,
                                    int mp_id,
                                    double dt,
                                    int &EXA_metric); 
                                    
int cm_write_tensor_restart(FILE *fp, const double *tensor);

int cm_read_tensor_restart(FILE *fp, double *tensor);

int constitutive_model_update_NR(Grid *grid,
                                 MaterialProperty *mat,
                                 FieldVariables *fv,
                                 LoadingSteps *load,
                                 const PGFem3D_opt *opts,
                                 const Multiphysics& mp,
                                 int mp_id,
                                 const double *dts,
                                 double alpha,
                                 int &EXA_metric);


#endif // #define PGFEM3D_CONSTITUTIVE_MODEL_H
