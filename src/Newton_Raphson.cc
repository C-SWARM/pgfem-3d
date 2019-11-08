#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "LINE.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "Newton_Raphson.h"
#include "PGFEM_io.h"
#include "PGFem3D_data_structure.h"
#include "bounding_element_utils.h"
#include "compute_reactions.h"
#include "constitutive_model.h"
#include "constitutive_model_handle.h"
#include "displacement_based_element.h"
#include "dynamics.h"
#include "energy_equation.h"
#include "enumerations.h"
#include "fd_increment.h"
#include "fd_residuals.h"
#include "get_dof_ids_on_elem.h"
#include "incl.h"
#include "integration.h"
#include "interface_macro.h"
#include "load.h"
#include "macro_micro_functions.h"
#include "matice.h"
#include "matrix_printing.h"
#include "ms_cohe_job_list.h"
#include "pgf_fe2_macro_client.h"
#include "pgf_fe2_micro_server.h"
#include "press_theta.h"
#include "res_fini_def.h"
#include "solver_file.h"
#include "stabilized.h"
#include "stiffmat_fd.h"
#include "subdivision.h"
#include "three_field_element.h"
#include "utils.h"
#include "vol_damage_int_alg.h"
#include "vtk_output.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <cassert>
#include <vector>

#ifndef NR_UPDATE
#define NR_UPDATE 0
#endif

#ifndef NR_PRINT_INTERMEDIATE
#define NR_PRINT_INTERMEDIATE 0
#endif

#ifndef NR_COMPUTE_REACTIONS
#define NR_COMPUTE_REACTIONS 0
#endif

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

#ifndef PFEM_DEBUG_ALL
#define PFEM_DEBUG_ALL 0
#endif

#ifndef DEBUG_MULTISCALE_SERVER
#define DEBUG_MULTISCALE_SERVER 1
#endif

vector<long> ODE_EXA_metric;
vector<long> dof_EXA_metric;

using namespace pgfem3d;
using namespace multiscale::net;

namespace {
using pgfem3d::Solver;
using pgfem3d::solvers::SparseSystem;
  
/* MINIMAL_OUTPUT prints a summary of the entire function call. For
 * any print_level > MINIMAL_OUTPUT, normal output is used. */
enum {MINIMAL_OUTPUT,NORMAL_OUTPUT,VERBOSE_OUTPUT};

const constexpr int periodic = 0;

// max micro substep criteria
const constexpr int    max_n_micro_substep = 2;
const constexpr double alpha_restart_ms    = 2.0;
const constexpr double max_damage_per_step = 0.05;
const constexpr double alpha_restart       = 1.25;

// printing info lines
constexpr const char line_level_0[] = "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::";
constexpr const char line_level_1[] = "======================================================================";
constexpr const char line_level_2[] = "----------------------------------------------------------------------";

struct NR_summary {
  int total_iter;
  double final_rel_res;
  double final_res;
};

struct NR_time_steps {
  double times[3];
  double dt[2];
  int tim;
 };
}


/// check element evolutions
///
/// compute element volumes at t(n) and t(n+1) and apply limit of
/// element volume changes. If this function is used in checking
/// physics based evolution alpha, additional subdivision might be
/// activated to evolve displacements within the range of volume
/// evolution is allowed.
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] com object for communication
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dts time step sizes a n, and n+1
/// \return non-zero on internal error
int is_displacement_acceptable(double *alpha,
                               Grid *grid,
                               FieldVariables *fv,
                               LoadingSteps *load,
                               const PGFem3D_opt *opts,
                               int mp_id)
{
  int err = 0;
  *alpha = 0.0;
  return err;
}

// @todo This function returns. I'm not sure what all the dead code below is
  //       for. This needs to be reviewed by @cp and eliminated if
  //       appropriate. LD
  //
  // SUPP sup = load->sups[mp_id];

  // // @todo This variable is never used, so I commented out this code. It
  // //       should be reviewed by an @cp person and eliminated if
  // //       appropriate. LD
  // //
  // // int is_total_lagrangian = 0;

  // // if(!sup->multi_scale)
  // //   is_total_lagrangian = 1;
  // // else
  // // {
  // //   switch(opts->analysis_type)
  // //   {
  // //     case DISP:
  // //     case TF:
  // //       is_total_lagrangian = 1;
  // //       break;
  // //     case CM:
  // //       if(opts->cm != UPDATED_LAGRANGIAN)
  // //         is_total_lagrangian = 1;
  // //       break;
  // //   }
  // // }

  // double alpha_V = 0.0;

  // for(int e=0; e<grid->ne; e++)
  // {
  //   int nne   = grid->element[e].toe;
  //   int ndofn = fv->ndofn;
  //   int ndofe = nne*ndofn;

  //   long *nod = (long *) malloc(sizeof(long)*nne);
  //   long *cn = aloc1l (ndofe);
  //   double *u_e = aloc1(ndofe);

  //   elemnodes(e,nne,nod,grid->element);
  //   get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,grid->node,cn,mp_id);
  //   def_elem_total(cn,ndofe,fv->u_np1,fv->d_u,grid->element,grid->node,load->sups[mp_id],u_e);

  //   double *X = aloc1(nne);
  //   double *Y = aloc1(nne);
  //   double *Z = aloc1(nne);

  //   double *x_n = aloc1(nne);
  //   double *y_n = aloc1(nne);
  //   double *z_n = aloc1(nne);

  //   double *x_np1 = aloc1(nne);
  //   double *y_np1 = aloc1(nne);
  //   double *z_np1 = aloc1(nne);

  //   for(int n=0; n<nne; n++)
  //   {
  //     long nid = nod[n];

  //     X[n]   = grid->node[nid].x1_fd;
  //     Y[n]   = grid->node[nid].x2_fd;
  //     Z[n]   = grid->node[nid].x3_fd;

  //     x_n[n]   = grid->node[nid].x1;
  //     y_n[n]   = grid->node[nid].x2;
  //     z_n[n]   = grid->node[nid].x3;

  //     x_np1[n] = grid->node[nid].x1_fd + u_e[n*ndofn + 0];
  //     y_np1[n] = grid->node[nid].x2_fd + u_e[n*ndofn + 1];
  //     z_np1[n] = grid->node[nid].x3_fd + u_e[n*ndofn + 2];
  //   }

  //   double ratio = 0.0;
  //   // reference configuration volume
  //   double V_0   = compute_volumes_from_coordinates(X,  Y,  Z,nne);
  //   // current configuration volume at t(n)
  //   double V_n   = compute_volumes_from_coordinates(x_n,  y_n,  z_n,nne);
  //   // current configuration volume at t(n+1)
  //   double V_np1 = compute_volumes_from_coordinates(x_np1,y_np1,z_np1,nne);

  //   double dV_max = 0.2; // limit volume changes within 20%
  //   ratio = fabs(V_np1-V_n)/V_n/dV_max;
  //   alpha_V = (alpha_V > ratio)? alpha_V : ratio;

  //   free(nod);
  //   free(cn);
  //   free(u_e);

  //   free(X);
  //   free(Y);
  //   free(Z);

  //   free(x_n);
  //   free(y_n);
  //   free(z_n);

  //   free(x_np1);
  //   free(y_np1);
  //   free(z_np1);
  // }

  // *alpha = alpha_V;
  // return err;


/// reset variables for Newton Raphson iteration
/// If Newton Raphson is restarted, reset variables to inital
///
/// \param[in] grid a mesh object
/// \param[in,out] fv object for field variables
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] physics_id, MULTIPHYSICS_MECHANICAL, MULTIPHYSICS_THERMAR, ...
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int reset_variables_for_NR(Grid *grid,
                           FieldVariables *fv,
                           CRPL *crpl,
                           int physics_id,
                           const PGFem3D_opt *opts)
{
  int err = 0;

  if(physics_id==MULTIPHYSICS_MECHANICAL)
  {
    switch(opts->analysis_type){
     case FS_CRPL: case FINITE_STRAIN:
      res_fini_def(grid->ne,fv->npres,grid->element,fv->eps,fv->sig,
                   crpl,opts->analysis_type);
      break;
     case STABILIZED:
      res_stab_def(grid->ne,fv->npres,grid->element,fv->eps,fv->sig,opts->stab);
      break;
     case MINI:
      MINI_reset(grid->element,grid->ne,fv->npres,fv->sig);
      break;
     case MINI_3F:
      MINI_3f_reset(grid->element,grid->ne,fv->npres,4,fv->sig,fv->eps);
      break;
     case CM: // intented flow
     case CM3F: 
      constitutive_model_reset_state(fv->eps, grid->ne, grid->element);
      break;
     default: break;
    }
    if(opts->analysis_type == CM3F || opts->analysis_type == TF)
      fv->tf.reset();
      
    /* reset damage variables */
    for(int ia=0; ia<grid->ne; ia++)
    {
      long n_ip;
      int_point(grid->element[ia].toe,&n_ip);
      for(int ja=0; ja<n_ip; ja++)
        reset_damage(&(fv->eps[ia].dam[ja]));

      if(opts->analysis_type == STABILIZED) // get pressure terms too
      {
        int_point(10,&n_ip);
        for(int ja=0; ja<n_ip; ja++){
          reset_damage(&(fv->eps[ia].st[ja].dam));
        }
      }
    }
  }

  for(int ia=0;ia<fv->ndofd;ia++)
    fv->dd_u[ia] = fv->d_u[ia] = fv->f_defl[ia] = fv->f[ia] = 0.0;

  return err;
}


/// When time step is subdivided, loading should be subdivided too.
/// this function updates loading values when it is needed according to
/// the subdivision history. SubdivisionScheme contains information
/// when the load needs to be updated.
///
/// \param[in] sp container of parameters for subdivision
/// \param[in, out] sup_defl displacement increments
/// \param[in] npd number of supported nodes (number of prescribed field variables)
/// \param[in, out] RRn nodal force at t(n)
/// \param[in, out] R nodal force
/// \param[in] ndofd number of degree of freedom in the domain
/// \return non-zero on internal error
int update_load_increments_for_subdivision(const SubdivisionScheme *sp,
                                           double *sup_defl,
                                           const int npd,
                                           double *RRn,
                                           double *R,
                                           const int ndofd)
{
  int err = 0;

  if(sp->need_to_update_loading)
  {
    double factor = sp->loading_factor;
    for(int ia=0; ia<npd; ia++){
      assert(sup_defl != NULL && "sup_defl can't be null");
      sup_defl[ia] = sup_defl[ia] - sup_defl[ia]*factor;
    }

    for(int ia=0; ia<ndofd; ia++)
    {
      RRn[ia] = RRn[ia] + R[ia]*factor;
      R[ia]   =   R[ia] - R[ia]*factor;
    }
  }

  return err;
}

/// Overloaded aboved function which uses vector<double> types in argument
/// instead of pointers
int update_load_increments_for_subdivision(const SubdivisionScheme *sp,
                                           vector<double>& sup_defl,
                                           const int npd,
                                           vector<double>& RRn,
                                           vector<double>& R,
                                           const int ndofd)
{
  int err = 0;

  if(sp->need_to_update_loading)
  {
    double factor = sp->loading_factor;
    
    for(int ia=0; ia<npd; ia++){
      sup_defl[ia] = sup_defl[ia] - sup_defl[ia]*factor;
    }

    for(int ia=0; ia<ndofd; ia++)
    {
      RRn[ia] = RRn[ia] + R[ia]*factor;
      R[ia]   =   R[ia] - R[ia]*factor;
    }
  }

  return err;
}


/// Compute residuals for Newton Raphson iteration
///
/// \param[in,out] non-zero on internal error
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] com object for communication
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dts time step sizes a n, and n+1
/// \param[in] updated_deformation if 1, compute resiual on updated deformation
/// \return compute time taken by this function 
double compute_residuals_for_NR(long *INFO, 
                                Grid *grid,
                                MaterialProperty *mat,
                                FieldVariables *fv,
                                Solver *sol,
                                LoadingSteps *load,
                                CRPL *crpl,
                                const CommunicationStructure *com,
                                const PGFem3D_opt *opts,
                                const Multiphysics& mp,
                                int mp_id,
                                double t,
                                double *dts,
                                int updated_deformation)
{
  double func_time = -CLOCK();

  switch(mp.physics_ids[mp_id])
  {
   case MULTIPHYSICS_MECHANICAL:
    *INFO = fd_residuals_MP(grid,mat,fv,sol,load,crpl,com,opts,mp,mp_id,t,dts,updated_deformation);
    break;
   case MULTIPHYSICS_THERMAL:
    *INFO = energy_equation_compute_residuals(grid,mat,fv,load,mp_id,updated_deformation,dts[DT_NP1]);
    break;
   case MULTIPHYSICS_CHEMICAL: //intented flow, not yet been implemented
   default:
    PGFEM_printf("%s is not defined (compute_residuals_for_NR)\n", mp.physicsname[mp_id]);

    // @todo This abort() was added because otherwise we're returning an
    //       uninitialized value. If this branch is a valid path then @cp
    //       should add a useful value for INFO here. LD
    abort();
  }

  func_time += CLOCK();
  return func_time;
}

/// Compute stiffnes
///
/// \param[in,out] non-zero on internal error
/// \param[in] max_substep microscale maximum subdivision number
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com communication object
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step
/// \param[in] iter number of Newton Raphson interataions
/// \param[in] myrank current process rank
/// \return compute time taken by this function 
double compute_stiffness_for_NR(long *INFO, 
                                int *max_substep,
                                Grid *grid,
                                MaterialProperty *mat,
                                FieldVariables *fv,
                                Solver *sol,
                                LoadingSteps *load,
                                const CommunicationStructure *com,
                                CRPL *crpl,
                                const PGFem3D_opt *opts,
                                const Multiphysics& mp,
                                int mp_id,
                                double dt,
                                long iter)
{
  double func_time = -CLOCK();
  long GInfo;
  int myrank = com->rank;
  
  if(sol->microscale == NULL) // Null the matrix (if not doing multiscale)
    sol->system->zero();

  switch(mp.physics_ids[mp_id])
  {
   case MULTIPHYSICS_MECHANICAL:
     *INFO = stiffmat_fd_MP(grid,mat,fv,sol,load,com,crpl,opts,mp,mp_id,dt,iter);
    break;
   case MULTIPHYSICS_THERMAL:
    *INFO = energy_equation_compute_stiffness(grid,mat,fv,sol,load,com,opts,mp_id,dt);
    break;
   case MULTIPHYSICS_CHEMICAL: //intented flow, not yet been implemented
   default:
    PGFEM_printf("%s is not defined (compute_stiffness_for_NR)\n", mp.physicsname[mp_id]);
  }

  // if INFO value greater than 0, the previous computation has an error
  com->net->allreduce(INFO,&GInfo,1,NET_DT_LONG,NET_OP_MAX,com->comm);

  if(GInfo > 0)
  {
    if(myrank == 0)
    {
      PGFEM_printf("Error detected (stiffmat_fd) %s:%s:%ld.\n"
                   "Subdividing load.\n", __func__, __FILE__, __LINE__);
    }
    *INFO = 1;
  }

  // turn off line search for server-style multiscale
  if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL)
  {
    // complete any jobs before assembly
    MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
    pgf_FE2_macro_client_recv_jobs(ctx->client,ctx->macro,max_substep);
  }

  // Matrix assmbly
  sol->system->assemble();

  func_time += CLOCK();
  return func_time;
}

/// During subdivided Newton iteration, shift solutions from t(n+1) -> t(n) -> t(n-1)
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] com object for communication
/// \param[in] VVolume original volume of the domain
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] NR_t container of time stepping info
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int update_values_for_next_NR(Grid *grid,
                              MaterialProperty *mat,
                              FieldVariables *fv,
                              Solver *sol,
                              LoadingSteps *load,
                              CRPL *crpl,
                              const CommunicationStructure *com,
                              const double VVolume,
                              const PGFem3D_opt *opts,
                              const Multiphysics& mp,
                              NR_time_steps *NR_t,
                              int mp_id)
{
  int err = 0;
  double t = NR_t->times[NR_t->tim+1];
  double dt = NR_t->dt[DT_NP1];

  switch(mp.physics_ids[mp_id])
  {
   case MULTIPHYSICS_MECHANICAL:
     {
       /* increment coheisve elements */
       if(opts->cohesive)
         increment_cohesive_elements(grid->nce,grid->coel,&(fv->pores),grid->node,load->sups[mp_id],fv->d_u,mp_id);

       /* Finite deformations increment */
       switch(opts->analysis_type){
        case FS_CRPL:
        case FINITE_STRAIN:
         fd_increment (grid->ne,grid->nn,fv->ndofn,fv->npres,mat->matgeom,mat->hommat,
                       grid->element,grid->node,load->sups[mp_id],fv->eps,fv->sig,fv->d_u,fv->u_np1,
                       sol->nor_min,crpl,dt,grid->nce,grid->coel,&(fv->pores),com,
                       VVolume,opts, mp_id);
         break;
        case STABILIZED:
         st_increment (grid->ne,grid->nn,fv->ndofn,fv->ndofd,mat->matgeom,mat->hommat,
                       grid->element,grid->node,load->sups[mp_id],fv->eps,fv->sig,fv->d_u,fv->u_np1,
                       sol->nor_min,opts->stab,dt,grid->nce,grid->coel,&(fv->pores),com,
                       opts->cohesive,mp_id);
         break;
        case MINI:
         MINI_increment(grid->element,grid->ne,grid->node,grid->nn,fv->ndofn,
                        load->sups[mp_id],fv->eps,fv->sig,mat->hommat,fv->d_u,com,mp_id);
         break;
        case MINI_3F:
         MINI_3f_increment(grid->element,grid->ne,grid->node,grid->nn,fv->ndofn,
                           load->sups[mp_id],fv->eps,fv->sig,mat->hommat,fv->d_u,com,mp_id);
         break;
        case DISP:
         DISP_increment(grid->element,grid->ne,grid->node,grid->nn,fv->ndofn,load->sups[mp_id],fv->eps,
                        fv->sig,mat->hommat,fv->d_u,fv->u_np1,com,mp_id);
         break;
        case TF:
         update_3f_output_variables(grid,mat,fv,load,mp_id,dt,t);
         break;
        case CM:
        case CM3F:
        {
          switch(opts->cm){
            case HYPER_ELASTICITY: case DISP:
              DISP_increment(grid->element,grid->ne,grid->node,grid->nn,fv->ndofn,load->sups[mp_id],fv->eps,
                             fv->sig,mat->hommat,fv->d_u,fv->u_np1,com,mp_id);
              break;
            case CRYSTAL_PLASTICITY: case BPA_PLASTICITY: case TESTING:
             /* updated later... */
              break;
            default: assert(0 && "undefined CM type"); break;
          }
          break;
        }
        default: break;
       }
       break;
     }
   case MULTIPHYSICS_THERMAL:
    update_thermal_flux4print(grid,mat,fv,dt);

    break;
   case MULTIPHYSICS_CHEMICAL:
    break;
  }
  /* Add deformation increment into displacement vector */
  vvplus(fv->u_np1,fv->d_u,fv->ndofd);

  // update previous time step values, u_n and u_nm1 from current
  //  time step values u_np1
  // For FE2, these vectors are not allocated and are passed as NULL
  if ( fv->u_nm1 != NULL && fv->u_n != NULL )
  {
    for(long a = 0; a<grid->nn; a++)
    {
      for(long b = 0; b<fv->ndofn; b++)
      {
        fv->u_nm1[a*fv->ndofn + b] = fv->u_n[a*fv->ndofn + b];
        long id = grid->node[a].id_map[mp_id].id[b];
        if((opts->analysis_type==CM || opts->analysis_type==CM3F) && opts->cm==UPDATED_LAGRANGIAN)
          // Updated Lagrangian
        {
          if(id>0)
            fv->u_n[a*fv->ndofn + b] = fv->d_u[id-1];
          else
          {
            if(id==0)
              fv->u_n[a*fv->ndofn + b] = fv->u0;
            else
              fv->u_n[a*fv->ndofn + b] = fv->u0 + (load->sups[mp_id])->defl_d[abs(id)-1];
          }
        }
        else
          // Total Lagrangian: DISP, TF
        {
          if(id>0)
            fv->u_n[a*fv->ndofn + b] = fv->u_np1[id-1];
          else
          {
            if(id==0)
              fv->u_n[a*fv->ndofn + b] = fv->u0;
            else
              fv->u_n[a*fv->ndofn + b] = fv->u0 + (load->sups[mp_id])->defl[abs(id)-1] + (load->sups[mp_id])->defl_d[abs(id)-1];
          }
        }
      }
    }
  }

  if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL)
  {
    /* update of internal fv and nodal coordinates */
    if(opts->analysis_type==CM || opts->analysis_type==CM3F)
    {
      switch(opts->cm){
       case UPDATED_LAGRANGIAN:
       case TOTAL_LAGRANGIAN:
        constitutive_model_update_time_steps(grid->element,grid->node,fv->eps,grid->ne,grid->nn,
                                             fv->ndofn,fv->u_n,dt,opts->cm,mp_id);
        break;
       case MIXED_ANALYSIS_MODE:
        constitutive_model_update_time_steps(grid->element,grid->node,fv->eps,grid->ne,grid->nn,
                                             fv->ndofn,fv->u_n,dt,1 /* TL */,mp_id);
        break;
       default: break;
      }
    }

    if(opts->analysis_type==TF || opts->analysis_type==CM3F)
    {
      fv->tf.update_np1_from_increments();
      if(opts->cm == UPDATED_LAGRANGIAN)
        fv->tf.update_for_next_time_step(true);
      else
        fv->tf.update_for_next_time_step();
    }
  }
  // Update time steps
  NR_t->dt[DT_N] = NR_t->dt[DT_NP1];

  /* Null prescribed increment deformation */
  for (int i=0;i<(load->sups[mp_id])->npd;i++){
    (load->sups[mp_id])->defl[i] += (load->sups[mp_id])->defl_d[i];
    (load->sups[mp_id])->defl_d[i] = 0.0;
  }

  for (int i=0;i<fv->ndofd;i++) {
    fv->d_u[i] = 0.0;
    fv->dd_u[i] = 0.0;
    fv->f_defl[i] = 0.0;
    fv->f[i] = 0.0;
  }
  return err;
}

/// check convergence on energy norm
///
/// Very small loading cases, residuals are very small, too, and difficult to converge to the tolarance with
/// a relative norm of the residual. Instead, PGFem3D checks convergence on energy norm which is computed as:
/// ||E|| = ||R*du||
///
/// \param[out] ENORM computed energy norm
/// \param[in] nor absolute residual
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] com container of communications info
/// \param[in] opts structure PGFem3D option
/// \param[in] tim time step id
/// \param[in] iter NR iteration id
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int check_energy_norm(double *ENORM,
                      double nor,
                      FieldVariables *fv,
                      Solver *sol,
                      const CommunicationStructure *com,
                      const PGFem3D_opt *opts,
                      int tim,
                      int iter)
{
  int is_converged = 0;
  int myrank = com->rank;
  double enorm, Genorm;
  // Compute the energy norm E_norm = abs(R ddu/(R0 ddu0))
  LToG(fv->dd_u,fv->BS_x,fv->ndofd,com);
  enorm = ss(fv->BS_f,fv->BS_x,com->DomDof[myrank]);
  com->net->allreduce(&enorm,&Genorm,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);

  if((tim == 0 && iter == 0))
    *ENORM = Genorm;

  if(*ENORM <= sol->computer_zero)
    *ENORM = 1.0;

  enorm = fabs(Genorm/(*ENORM));

  if (myrank == 0)
  {
    PGFEM_printf("INFO: Energy norm = |R'du|/|R0'du0| = [%1.12e] || [%1.12e]\n",enorm,fabs(Genorm));
    fflush(PGFEM_stdout);
  }

  // Check energy norm for convergence
  if(opts->solution_scheme_opt[CVG_CHECK_ON_ENERGY_NORM])
  {
    double err_err = (sol->nor_min)*(sol->nor_min);
    if((enorm < err_err) && (nor < 50*err_err) && iter > 0)
    {
      if(myrank == 0)
        PGFEM_printf("Converged on energy norm\n");
      is_converged = 1;
    }
  }

  return is_converged;
}

/// Set initial velocity to u(n-1)
///
/// \parma[in] u_nm1 saved u(n-1) values
/// \param[in] grid  a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] load object for loading
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dts time step size at t(n), t(n+1); dts[DT_N] = t(n) - t(n-1), dts[DT_NP1] = t(n+1) - t(n),
void apply_initial_velocity_to_nm1(double *u_nm1,
                                   Grid *grid,
                                   MaterialProperty *mat,
                                   FieldVariables *fv,
                                   LoadingSteps *load,
                                   const PGFem3D_opt *opts,
                                   int mp_id,
                                   double *dts)
{
  return;
  // this function is not verified and disabled, but will be used
  // after quick verification test after full implemention of the CM3F with inertia
  for(int ia=0; ia<grid->nn; ia++){
    for(long ib = 0; ib<grid->nsd; ib++){
      long uid = ia*fv->ndofn + ib;
      long id = grid->node[ia].id_map[mp_id].id[ib];
      double u_np1;
      if(id>0)
        u_np1 = fv->u_np1[id-1] + fv->d_u[id-1];
      else
      {
        if(id==0)
          u_np1 = fv->u0;
        else
          u_np1 = fv->u0 + (load->sups[mp_id])->defl[abs(id)-1] + (load->sups[mp_id])->defl_d[abs(id)-1];
      }
      double v0 = (fv->u_n[uid] - u_nm1[uid])/dts[DT_N];
      double a = (u_np1 -v0*dts[DT_NP1] - fv->u_n[uid])/dts[DT_NP1]/dts[DT_NP1];
      fv->u_nm1[uid] = a*dts[DT_N]*dts[DT_N] - v0*dts[DT_N] + fv->u_n[uid]; 
    }        
  }
  if(opts->analysis_type == TF || opts->analysis_type == CM3F)    
    compute_3f_initial_conditions(grid, mat, fv);
   
}

/// Actual Newton Raphson scheme with line search
///
/// \parma[out] solve_time time measure of spent for linear solver
/// \parma[out] stiffmat_loc_time - compute time for constructing A in Ax=b  
/// \parma[out] residuals_loc_time - compute time for constructing b in Ax=b
/// \parma[out] alpha physics based evolution parameter (maximum value from each physics)
/// \param[out] NOR Normalize norm of the residuals
/// \param[out] gam flag for denoting Line search that decreases increment of displacements
///                 in order to make solultion to converge
/// \param[out] ART if 1 Line search was activated
/// \param[out] NR_itr_no number of NR iterations taken in this NR step
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com container of communications info
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] times array of times
///                  t(n-1) = times[tim-1], t(n) = times[tim], t(n+1) = times[tim+1]
/// \param[in] dts time step size at t(n), t(n+1); dts[DT_N] = t(n) - t(n-1), dts[DT_NP1] = t(n+1) - t(n),
/// \param[in] tim time step id
/// \param[in] STEP number of subdivided steps
/// \param[in] DIV subdivision step id
/// \param[in, out] usage struct to get resource usage
/// \return non-zero on internal error
long Newton_Raphson_with_LS(double *solve_time,
                            double *stiffmat_loc_time,
                            double *residuals_loc_time,
                            double *alpha,
                            double *NOR,
                            long *gam,
                            int *ART,
                            int *NR_itr_no,
                            Grid *grid,
                            MaterialProperty *mat,
                            FieldVariables *fv,
                            Solver *sol,
                            LoadingSteps *load,
                            const CommunicationStructure *com,
                            CRPL *crpl,
                            const PGFem3D_opt *opts,
                            const Multiphysics& mp,
                            int mp_id,
                            double *times,
                            double *dts,
                            int tim,
                            long STEP,
                            long DIV,
                            rusage *usage,
                            long master_tim)
{
  long INFO = 0;
  *alpha = 0.0;
  
  double dt = dts[DT_NP1];
  double t = times[tim+1];
  
  // do Taylor updates
  if(sol->FNR == 5 && mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL)
  {
    double F[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
    double *dF  = load->sups[mp_id]->defl_d;
    for(int ia=0; ia<9; ia++)
      F[ia] = F[ia] + load->sups[mp_id]->defl[ia] + load->sups[mp_id]->defl_d[ia];
    
    double X[3] = {};
    double u[3] = {};
    for(int n = 0; n<grid->nn; ++n) {
      X[0] = grid->node[n].x1;
      X[1] = grid->node[n].x2;
      X[2] = grid->node[n].x3;
  
      u[0] = dF[0]*X[0] + dF[1]*X[1] + dF[2]*X[2];
      u[1] = dF[3]*X[0] + dF[4]*X[1] + dF[5]*X[2];
      u[2] = dF[6]*X[0] + dF[7]*X[1] + dF[8]*X[2];
      
      for(int ia=0; ia<fv->ndofn; ++ia){
        int II = grid->node[n].id_map[0].id[ia];
        if (II > 0)
          fv->f[II-1] = fv->dd_u[II-1] = fv->d_u[II-1] = u[ia];
      }
    }
    
    if(opts->analysis_type == CM || opts->analysis_type == CM3F)
    {
      for(long eid=0; eid<grid->ne; eid++)
      {
        long nne = grid->element[eid].toe;
        long nint = 0;
        int_point(nne, &nint);
        for(int ip=0; ip<nint; ++ip){
          Constitutive_model *m = &(fv->eps[eid].model[ip]);
          CM_Ctx cm_ctx;
          cm_ctx.set_tensors_ss(F,NULL,NULL,false);
          cm_ctx.set_time_steps_ss(dts[DT_NP1],dts[DT_N]);
          m->param->integration_algorithm(m, cm_ctx);
        }
      }
      // check from constitutive mode
      cm_get_subdivision_parameter(alpha, grid->ne, grid->element, fv->eps, dt);
    }  

    return INFO;
  }    

  int iter = 0;
  int myrank = com->rank;
  
  double nor, nor2, GNOR;
  nor = nor2 = GNOR = 10.0;

  long GInfo;
  const int NR_REBALANCE = (opts->no_migrate)? FE2_REBALANCE_NONE : FE2_REBALANCE_ADAPTIVE;

  if(sol->FNR == 1 || sol->FNR == 0)
    assert(opts->solverpackage == HYPRE || opts->solverpackage == TRILINOS);

  double alpha_ms = 0.0;
  double ENORM = 1.0;

  double max_damage = 0.0;
  double dissipation = 0.0;
  
  Matrix<double> u_nm1;

  if(master_tim==0 && tim==0 && fv->apply_initial_velocity)
  {
    u_nm1.initialization(fv->ndofd,1,0.0);
    
    for(int ia=0; ia<fv->ndofd; ia++)
      u_nm1.m_pdata[ia] = fv->u_nm1[ia];
        
    apply_initial_velocity_to_nm1(u_nm1.m_pdata, grid, mat, fv, load, opts, mp_id, dts);
  }
      

  while(1) // do until converge
  {
    int max_substep = 0;
    INFO = 0;
    if(sol->FNR == 1 || (sol->FNR == 0 && iter == 0))
    {
      //compute stiffness matrix
      *stiffmat_loc_time += compute_stiffness_for_NR(&INFO,&max_substep,grid,mat,fv,sol,load,com,
                                                     crpl,opts,mp,mp_id,dt,iter);
      if(INFO > 0)
      {
        *ART = 1;
        break; // goto rest
      }

      // turn off line search for server-style multiscale
      if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL)
        *ART = 1;
    }

    /*=== Solve the system of equations ===*/
    SOLVER_INFO s_info;

    *solve_time += sol->system->solveSystem(opts, fv->BS_f, fv->BS_x, tim, iter,
                                            com->DomDof, &s_info);

    // sol->system->printWithRHS("test_", fv->BS_f, com->DomDof[myrank],
    //                                myrank);
    // exit(0);

    if(myrank == 0)
      solve_system_check_error(PGFEM_stdout,s_info);

    double BS_nor = s_info.res_norm;
    int BS_iter   = s_info.n_iter;

    // Check for correct solution
    if (!isfinite(BS_nor))
    {
      INFO = 1;
      *ART = 1;
      if (myrank == 0)
        PGFEM_printf("ERROR in the solver: nor = %f\n",BS_nor);
      break; // goto rest
    }

    if (BS_nor > 500.*(sol->err) || BS_iter < 0)
    {
      INFO = 1;
      *ART = 1;
      if(myrank == 0)
        PGFEM_printf("ERROR in the solver: nor = %8.8e || iter = %d\n",BS_nor,BS_iter);
      break; // goto rest
    }

    /* Transform GLOBAL displacement vector to LOCAL */
    GToL (fv->BS_x, fv->dd_u, fv->ndofd, com);

    /* LINE SEARCH */
    double Gss_temp = 0.0;
    double tmp  = ss(fv->BS_f,fv->BS_f,com->DomDof[myrank]);
    com->net->allreduce(&tmp,&Gss_temp,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
    double LS1 = 1./2.*Gss_temp;
    
    if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL){
    // Pressure and volume change THETA

      switch(opts->analysis_type){
       case FS_CRPL:
       case FINITE_STRAIN:
        press_theta (grid->ne,fv->ndofn,fv->npres,grid->element,grid->node,fv->d_u,fv->dd_u,load->sups[mp_id],mat->matgeom,
                     mat->hommat,fv->eps,fv->sig,iter,sol->nor_min,dt,crpl,opts,mp_id);
        break;
       case MINI:
        MINI_update_bubble(grid->element,grid->ne,grid->node,fv->ndofn,load->sups[mp_id],
                           fv->eps,fv->sig,mat->hommat,fv->d_u,fv->dd_u,iter,mp_id);
        break;
       case MINI_3F:
        MINI_3f_update_bubble(grid->element,grid->ne,grid->node,fv->ndofn,load->sups[mp_id],
                              fv->eps,fv->sig,mat->hommat,fv->d_u,fv->dd_u,iter,mp_id);
        break;
       case TF:
          update_3f_NR(grid,mat,fv,load,opts,mp_id,dts,sol->alpha);
        break;
       case CM3F:
        constitutive_model_update_NR(grid, mat, fv, load, opts, mp, mp_id, dts, sol->alpha);
       default:
        break;
      }
    }
    
    /*************************/
    /* INTEGRATION ALGORITHM */
    /*************************/
    if (opts->analysis_type == FS_CRPL)
    {
      INFO = integration_alg (grid->ne,fv->ndofn,fv->ndofd,fv->npres,crpl,grid->element,
                              grid->node,fv->d_u,fv->dd_u,load->sups[mp_id],mat->matgeom,mat->hommat,
                              fv->eps,fv->sig,tim,iter,dt,sol->nor_min,STEP,0,opts,mp_id);

      /* Gather INFO from all domains */
      // if INFO value greater than 0, the previous computation has an error
      com->net->allreduce(&INFO,&GInfo,1,NET_DT_LONG,NET_OP_MAX,com->comm);
      if (GInfo > 0) // if not 0, an error is detected
      {
        INFO = 1;
        *ART = 1;
        break; // goto rest
      }
    }

    /* Update deformations */
    for (int i=0;i<fv->ndofd;i++) {
      fv->f[i] = fv->d_u[i] + fv->dd_u[i];
      fv->f_u[i] = 0.0;
    }

    INFO = vol_damage_int_alg(grid->ne,fv->ndofn,fv->f,fv->u_np1,grid->element,grid->node,
                              mat->hommat,load->sups[mp_id],dt,iter,com,
                              fv->eps,fv->sig,&max_damage,&dissipation,
                              opts->analysis_type,mp_id);

    bounding_element_communicate_damage(grid->n_be,grid->b_elems,grid->ne,fv->eps,com);
    /* this is not verified and currently active, but fully implemented.
       if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL)
       {
       double element_volume_evolution = 0.0;

       INFO += is_displacement_acceptable(&element_volume_evolution,grid,fv,load,opts,mp_id);
       alpha = 10.0;
       }*/
    // if INFO value greater than 0, the previous computation has an error
    com->net->allreduce(&INFO,&GInfo,1,NET_DT_LONG,NET_OP_MAX,com->comm);
    if(GInfo > 0) // if not 0, an error is detected
    {
      INFO = 1;
      *ART = 1;
      if(myrank == 0)
        PGFEM_printf("Inverted element detected (vol_damage_int_alg).\n");
      break; // goto rest
    }

    /* server-style multiscale */
    if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL)
    {
      /* zero the macroscale tangent */
      sol->system->zero();

      /* start the microscale jobs */
      MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
      if ( iter == 0 ) {
        /* == Do not rebalance if this is the first iteration. ==
         * This is because most often it will be after an update or
         * print operation from the previous step. These operations
         * are approx. equal for all cells and thus the timing
         * information is not valid for rebalancing purposes. */
        pgf_FE2_macro_client_rebalance_servers(ctx->client,ctx->mscom,
                                               FE2_REBALANCE_NONE);
      } else {
        pgf_FE2_macro_client_rebalance_servers(ctx->client,ctx->mscom,
                                               NR_REBALANCE);
      }
    
      unsigned int i;
      double temp_defl[9];
      double newDt;
      //dt0 is the size of a full (nondivided) step
      double dt0 = times[tim+1] - times[0];      

      //if it is the start of a subdivision group, make sure to record that first number
      //t(n+1) = t(n)' + (T (n+1) -t'(n)) * (step + 1)/NS
      //we need tn' to stay constant at whatever it was at step = 0
      if (DIV == 0) {
        if (tim > 0) {
         ctx->macro->sol->times[ctx->macro->sol->tim - 1]  = ctx->macro->sol->times[ctx->macro->sol->tim];
        }
      }
      //if there is no subdivision, then dt is straightforward. Otherwise, it needs to be calculated
      if (STEP > 0) {
        double dtFactor =  (times[tim + 1] - ctx->macro->sol->times[ctx->macro->sol->tim - 1]) * (((double) DIV + 1.0)/((double) STEP));
        ctx->macro->sol->times[ctx->macro->sol->tim + 1] = ctx->macro->sol->times[ctx->macro->sol->tim - 1] + dtFactor;
        newDt = ctx->macro->sol->times[ctx->macro->sol->tim - 1] + dtFactor - ctx->macro->sol->times[ctx->macro->sol->tim];
      } else {
        ctx->macro->sol->times[ctx->macro->sol->tim + 1] = ctx->macro->sol->times[ctx->macro->sol->tim] + dt;
        newDt = dt;
      }
      //once we know the size of the particular step, we can use it to scale grad u
      for (i = 0; i < sizeof(load->sup_defl[mp_id]) ; i++) {
        temp_defl[i] = load->sup_defl[mp_id][i];
        load->sup_defl[mp_id][i] = (newDt/(dt0))*load->sup_defl[mp_id][i];
      }
      pgf_FE2_macro_client_send_jobs(ctx->client,ctx->mscom,ctx->macro,
                                     JOB_COMPUTE_EQUILIBRIUM);
      //gradU needs to be put back afterwards
      for (i = 0; i < sizeof(load->sup_defl[mp_id]) ; i++) {
       load->sup_defl[mp_id][i] = temp_defl[i];
      }


    }

    /* Residuals */
    INFO = 0;
    *residuals_loc_time += compute_residuals_for_NR(&INFO,grid,mat,fv,sol,load,crpl,
                                                    com,opts,mp,mp_id,t,dts, 1);

    // if INFO value greater than 0, the previous computation has an error
    com->net->allreduce(&INFO,&GInfo,1,NET_DT_LONG,NET_OP_MAX,com->comm);
    if (GInfo > 0)// if not 0, an error is detected
    {
      INFO = 1;
      *ART = 1;
      if(myrank == 0)
        PGFEM_printf("Error detected (fd_residuals)\n");

      break; // goto rest
    }

    if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL){
      /* print_array_d(PGFEM_stdout,f_u,ndofd,1,ndofd); */
      MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
      pgf_FE2_macro_client_recv_jobs(ctx->client,ctx->macro,&max_substep);

      /* determine substep factor */
      alpha_ms = ((double) max_substep) / max_n_micro_substep;
      if(alpha_ms > alpha_restart_ms)
      {
        INFO = 1;
        *ART = 1;
        *alpha = alpha_ms;
        if(myrank == 0)
          PGFEM_printf("Too many subdvisions at microscale (alpha_ms = %f).\n",alpha_ms);

        break; // goto rest
      }
    }

    /* Transform LOCAL load vector to GLOBAL */
    LToG (fv->f_u, fv->BS_f_u, fv->ndofd, com);

    /* Compute Euclidian norm */
    for (int i=0; i<com->DomDof[myrank]; i++)
      fv->BS_f[i] = fv->BS_RR[i] - fv->BS_f_u[i];

    nor  = ss (fv->BS_f,fv->BS_f,com->DomDof[myrank]);
    com->net->allreduce(&nor,&GNOR,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
    nor2 = nor = sqrt(GNOR);

    // Reset *NORM if less than convergence tolerance MM 6/27/2012*/
    // take maximum
    if ((tim == 0 && iter == 0) || (iter == 0 && fv->NORM < sol->nor_min))
    {
      if(nor > fv->NORM)
        fv->NORM = nor;
      if(fv->NORM<sol->computer_zero)
        fv->NORM = sol->computer_zero;
    }

    nor /= fv->NORM; // Normalize norm

    if (!isfinite(nor))
    {
      INFO = 1;
      *ART = 1;
      if (myrank == 0)
        PGFEM_printf("ERROR in the algorithm : nor = %f\n",nor);
      break; // goto rest
    }

    // My Line search
    if (*ART == 0)
    {
      INFO = LINE_S3_MP(residuals_loc_time,grid,mat,fv,sol,load,com,crpl,opts,mp,
                        dts,t,mp_id,&nor,&nor2,fv->NORM,LS1,iter,&max_damage,&dissipation,
                        tim,STEP);

      // Gather infos
      // if INFO value greater than 0, the previous computation has an error
      com->net->allreduce(&INFO,&GInfo,1,NET_DT_LONG,NET_OP_MAX,com->comm);
      // ERROR in line search
      if (GInfo > 0)
      {
        if (myrank == 0)
          PGFEM_printf("Error in the line search algorithm\n");

        if(*NOR < 5.0*(sol->nor_min)) //MM 10/2/2012 soft convergence criteria
        {
          if(myrank == 0)
            PGFEM_printf("I will take last solution.\n");

          nor = *NOR;
          sol->gama = 0.0;
          INFO = 0;
          break; // take previous as converged value

        }
        else
        {
          INFO = 1; // MM 10/2/2012 not converged, subdivide
          break; //goto rest
        }
      }

      if (sol->gama != 1.0) *gam = 1;

    }
    else // no line search
      sol->gama = 1.0;

    // Total deformation increment
    *NOR = nor;
    for (int i=0;i<fv->ndofd;i++)
    {
      fv->d_u[i] += (sol->gama)*fv->dd_u[i];
      fv->dd_u[i] *= sol->gama;
    }
    
    if(master_tim==0 && tim==0 && fv->apply_initial_velocity){
      apply_initial_velocity_to_nm1(u_nm1.m_pdata, grid, mat, fv, load, opts, mp_id, dts);
    }
    
    
    
    if((mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL) && 
       (opts->analysis_type == TF || opts->analysis_type == CM3F))
       fv->tf.update_increments_from_NR(sol->gama);

    if (myrank == 0){
      getrusage (RUSAGE_SELF,usage);
      PGFEM_printf("(%ld) IT = %d : R = %8.8e :: ||f||/||f0||",
                   iter,BS_iter,BS_nor);
      PGFEM_printf(" = [%8.8e] || [%8.8e] :: S %ld.%ld, U %ld.%ld\n",
                   nor, nor2, usage->ru_stime.tv_sec,usage->ru_stime.tv_usec,
                   usage->ru_utime.tv_sec, usage->ru_utime.tv_usec);
      fflush(PGFEM_stdout);
    }

    if(NR_UPDATE || PFEM_DEBUG || PFEM_DEBUG_ALL){
      if(opts->analysis_type == MINI){
        MINI_check_resid(fv->ndofn,grid->ne,grid->element,grid->node,mat->hommat,fv->eps,
                         fv->sig,fv->d_u,load->sups[mp_id],fv->RR,com->DomDof,fv->ndofd,
                         com,mp_id);
      }
      if(opts->analysis_type == MINI_3F){
        MINI_3f_check_resid(fv->ndofn,grid->ne,grid->element,grid->node,mat->hommat,fv->eps,
                            fv->sig,fv->d_u,load->sups[mp_id],fv->RR,com->DomDof,fv->ndofd,
                            com,mp_id);
      }
    }

    //EXA metric computation
    if (opts->print_EXA){
      long global_ODE_EXA_metric = 0;
      com->net->reduce(&perIter_ODE_EXA_metric,&global_ODE_EXA_metric,1,NET_DT_LONG,NET_OP_SUM,0,com->comm);
      if (myrank == 0){
        if (opts->print_EXA == 2){
          PGFEM_printf ("ODE operations: (%ld)\n", global_ODE_EXA_metric);
          PGFEM_printf ("DoFs (fv->Gndof * # of ITs): %d * %d = %ld\n", fv->Gndof, BS_iter+1, fv->Gndof * (BS_iter+1));
        }
        ODE_EXA_metric[mp_id] += global_ODE_EXA_metric;          //accumulate exascale metrics per NR iteration
        dof_EXA_metric[mp_id] += fv->Gndof * (BS_iter+1);
      }
    }
    perIter_ODE_EXA_metric = 0;                                 //reset ODE operation counter
    
    if(check_energy_norm(&ENORM,nor2,fv,sol,com,opts,tim,iter))
    {
      INFO = 0;
      break;
    }

    /* Max number of iterations restart */
    if (iter > (sol->iter_max-1) && nor > sol->nor_min)
    {
      if (nor > 20.0*(sol->nor_min) || iter > (sol->iter_max + 2))
      {
        if (nor < 5.0*(sol->nor_min))
        {
          if (myrank == 0)
            PGFEM_printf ("I will take it\n");
          nor = 0.0;
        }
        else
        {
          if (myrank == 0)
            PGFEM_printf ("Error in the iteration : iter > iter_max (%d, %d)\n", sol->iter_max, iter);
          INFO = 1;
          if(*gam == 0)
            *ART = 1;
          break; //goto rest
        }
      }
    } /* end iter > iter_max */
    
    if(nor <= sol->nor_min || nor2 <= sol->computer_zero)
      break;

    iter++;
    *NR_itr_no = iter;


  }/* end while nor > nor_min */

  if(INFO==0)
  {
    // before increment after convergence, check max damage */
    if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL)
    {
      // check from constitutive mode
      if(opts->analysis_type == CM || opts->analysis_type == CM3F)
        cm_get_subdivision_parameter(alpha, grid->ne, grid->element, fv->eps, dt);
      else
        *alpha = max_damage/max_damage_per_step;
      // check displacement increment
      double element_volume_evolution = 0.0;
      is_displacement_acceptable(&element_volume_evolution,grid,fv,load,opts,mp_id);
      *alpha = (*alpha > element_volume_evolution)? *alpha : element_volume_evolution;
    }

    com->net->allreduce(NET_IN_PLACE,alpha,1,NET_DT_DOUBLE,NET_OP_MAX,com->comm);
    if(myrank == 0)
    {
      PGFEM_printf("Physics based evolution thresh (e.g. damage): %f (wmax: %f)\n"
                   "Microscale subdivision alpha_ms: %f (max_substep: %d)\n",
                   *alpha,max_damage_per_step,alpha_ms,max_n_micro_substep);
    }
    //after the solution is found, that time should be recorded
    if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL){
      MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
      ctx->macro->sol->times[ctx->macro->sol->tim] = ctx->macro->sol->times[ctx->macro->sol->tim + 1];
    }
    if(*alpha > alpha_restart || alpha_ms > alpha_restart_ms)
    {
      if(myrank == 0)
        PGFEM_printf("Subdividing to maintain accuracy of the material response.\n");

      *alpha = (*alpha > alpha_ms)? *alpha : alpha_ms;
      INFO = 1;
      *ART = 1;
    }
  }

  // increment converged, output the volume weighted dissipation */
  if(INFO==0)
  {
    // always adapt time step based on largest adaption parameter
    *alpha = (*alpha > alpha_ms)? *alpha : alpha_ms;

    if(myrank == 0)
    {
      com->net->reduce(NET_IN_PLACE,&dissipation,1,NET_DT_DOUBLE,NET_OP_SUM,0,com->comm);
      PGFEM_printf("Dissipation (vol weighted) [current] ||"
                   " [integrated] || (dt): %2.8e %2.8e %2.8e\n",
                   dissipation/dt, dissipation,dt);
    }
    else
    {
      double tmp = 0;
      com->net->reduce(&dissipation,&tmp,1,NET_DT_DOUBLE,NET_OP_SUM,0,com->comm);
    }
    sol->last_residual = nor2;
  }
  if(INFO>0 && master_tim==0 && tim==0 && fv->apply_initial_velocity)
  {
    for(int ia=0; ia<fv->ndofd; ia++)
      fv->u_nm1[ia] = u_nm1.m_pdata[ia];    
  }

  return INFO;
}

/// Perform Newton Raphson iteration. If not converge, do subdivision
///
/// \parma[out] solve_time time measure of spent for linear solver
/// \parma[out] stiffmat_loc_time - compute time for constructing A in Ax=b  
/// \parma[out] residuals_loc_time - compute time for constructing b in Ax=b
/// \param[in] print_level print level for a summary of the entire function call
/// \param[out] is_NR_converged if 1, NR was successful
/// \param[out] alpha physics based evolution rate
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] VVolume original volume of the domain
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] NR_t container of time stepping info
/// \param[in] mp_id mutiphysics id
void perform_Newton_Raphson_with_subdivision(double *solve_time,
                                             double *stiffmat_loc_time, 
                                             double *residuals_loc_time, 
                                             const int print_level,
                                             int *is_NR_converged,
                                             double *alpha,
                                             Grid *grid,
                                             MaterialProperty *mat,
                                             FieldVariables *FV,
                                             Solver *SOL,
                                             LoadingSteps *load,
                                             CommunicationStructure *COM,
                                             TimeStepping *time_steps,
                                             CRPL *crpl,
                                             const double VVolume,
                                             const PGFem3D_opt *opts,
                                             const Multiphysics& mp,
                                             NR_time_steps *NR_t,
                                             int mp_id)
{
  *is_NR_converged = 1;
  *alpha = 0.0;
  int iter = 0;
  
  // use pointers for physics[mp_id]
  Solver                 *sol = SOL + mp_id;
  FieldVariables          *fv = FV  + mp_id;
  CommunicationStructure *com = COM + mp_id;

  sol->gama = 0.0;

  int myrank = com->rank;
  long tim = NR_t->tim;

  if(myrank==0)
    PGFEM_printf(":: t(n) = %e, dts = (%e %e)\n", NR_t->times[tim], NR_t->dt[DT_N], NR_t->dt[DT_NP1]);

  double *dts = NR_t->dt; // short hand of dts
  double dt   = NR_t->dt[DT_NP1];
  double t    = NR_t->times[tim+1];
  double tn_0 = NR_t->times[tim];
  double dt0  = dt;
  fv->subdivision_factor_n   = 0.0;
  fv->subdivision_factor_np1 = 1.0;

  long i, j, INFO, gam = 0;
  int ART;
  double NOR=10.0;
  struct rusage usage;

  double max_damage  = 0.0; // damage substep criteria
  double dissipation = 0.0; // damage dissipation

  /* option '-no-migrate' */
  // const int NR_REBALANCE = (opts->no_migrate)? FE2_REBALANCE_NONE : FE2_REBALANCE_ADAPTIVE;

  switch(opts->analysis_type)
  {
   case STABILIZED:
   case MINI:
   case MINI_3F:
    fv->ndofn = 4;
    break;
   default:
    break;
  }

  sol->n_step = 0;

  /* SUBDIVISION */
  SubdivisionScheme sp(opts->solution_scheme_opt[NO_SUBDIVISION_LIMITS]);  
  INFO = ART = 0;

  fflush(PGFEM_stdout);

  while(1)
  {
    if(INFO==1 && opts->solution_scheme_opt[LINE_SEARCH]==0)
    {
      ART = 1;
      if(myrank==0)
        PGFEM_printf("Imposed to use NO Line search [INFO = %ld, ART = %d]\n", INFO, ART);
    }

    if (INFO == 1 && ART == 0)
    {
      reset_variables_for_NR(grid, fv, crpl, mp.physics_ids[mp_id], opts);
      if(myrank == 0 ) PGFEM_printf("\n** Try without LINE SEARCH **\n\n");

      ART = 1;
    }
    else
    {
      // run subdivision
      sp.do_subdivision(INFO,&dt,NR_t->times,tim,iter,sol->iter_max,*alpha,com);

      // update variables according to the subdivision
      if(sp.reset_variables)
        reset_variables_for_NR(grid, fv, crpl, mp.physics_ids[mp_id], opts);

      update_load_increments_for_subdivision(&sp,load->sup_defl[mp_id],(load->sups[mp_id])->npd,
                                             fv->RRn,fv->R,fv->ndofd);

      if(sp.is_subdivided)
        sol->is_subdivided = sp.step_size;

      if(sp.step_size>sol->max_subdivision && sol->max_subdivision > 0)
      {
        if(myrank==0)
          PGFEM_printf("maximum subdivision no = %d, requested = %d\n", sol->max_subdivision, sp.step_size);

        *is_NR_converged = 0;
        break;
      }

      dts[DT_NP1] = dt; // update dt_np1
      // dt_n is updated when Newton Raphson is completed without error

      gam = ART = 0;

    }

    /* recompute the microscale tangent if restart due to error. */
    if(INFO == 1 && DEBUG_MULTISCALE_SERVER && sol->microscale != NULL){
      /* zero the macroscale tangent */
      sol->system->zero();

      /* start the microscale jobs. Do not compute equilibrium. Use
         d_r (no displacement increments) for displacement dof
         vector */
      MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
      pgf_FE2_macro_client_rebalance_servers(ctx->client,ctx->mscom,
                                             FE2_REBALANCE_NONE);

      unsigned int i = 0;
      double temp_defl[9];
      double newDt;
      double dt0 = NR_t->times[tim+1] - NR_t->times[0];
      int DIV = sp.step_id;
      int STEP = sp.step_size;
      if (DIV == 0) {
        if (tim > 0) {
         ctx->macro->sol->times[ctx->macro->sol->tim - 1]  = ctx->macro->sol->times[ctx->macro->sol->tim];
        }
      }
      if (STEP > 0) {
        double dtFactor =  (NR_t->times[tim + 1] - ctx->macro->sol->times[ctx->macro->sol->tim - 1]) * (((double) DIV + 1.0)/((double) STEP));
        ctx->macro->sol->times[ctx->macro->sol->tim + 1] = ctx->macro->sol->times[ctx->macro->sol->tim - 1] + dtFactor;
        newDt = ctx->macro->sol->times[ctx->macro->sol->tim - 1] + dtFactor - ctx->macro->sol->times[ctx->macro->sol->tim];
      } else {
        ctx->macro->sol->times[ctx->macro->sol->tim + 1] = ctx->macro->sol->times[ctx->macro->sol->tim] + dt;
        newDt = dt;
      }

      for (i = 0; i < sizeof(load->sup_defl[mp_id]) ; i++) {
        temp_defl[i] = load->sup_defl[mp_id][i];
        load->sup_defl[mp_id][i] = (newDt/(dt0))*load->sup_defl[mp_id][i];
      }
      ctx->macro->sol->times[ctx->macro->sol->tim+1] = NR_t->times[tim] + dt;
      pgf_FE2_macro_client_send_jobs(ctx->client,ctx->mscom,ctx->macro,
                                     JOB_NO_COMPUTE_EQUILIBRIUM);

      for (i = 0; i < sizeof(load->sup_defl[mp_id]) ; i++) {
       load->sup_defl[mp_id][i] = temp_defl[i];
      }
    }

    /* reset the error flag */
    INFO = 0;
    int is_sub_cvg = 1;

    while (sp.step_size > sp.step_id)
    {
      fv->subdivision_factor_np1 = (NR_t->times[tim] + dt*(sp.step_id+1) - tn_0)/dt0;
      if(myrank==-1)
        PGFEM_printf(":: subdivision fators : t(n) = %e , t(n+1) = %e\n", fv->subdivision_factor_n, fv->subdivision_factor_np1);

      if (sp.is_subdivided && myrank == 0 )
      {
        PGFEM_printf("\nSTEP = %ld :: NS =  %ld || Time %e | dt = %e\n",
                     sp.step_id,sp.step_size,NR_t->times[tim+1],dt);
      }
      bounding_element_communicate_damage(grid->n_be,grid->b_elems,grid->ne,fv->eps,com);
      if (periodic == 1)
      {
        long nt = time_steps->nt;
        double factor = (NR_t->times[tim] + (sp.step_id+1)*dt)/NR_t->times[nt] * fv->eps[0].load;
        /* Plane strain deviatoric tension */
        if (fv->eps[0].type == 1){
          fv->eps[0].F[0][0] = NR_t->times[nt]/(NR_t->times[nt]
                                                - (NR_t->times[tim] + (sp.step_id+1)*dt)*fv->eps[0].load);
          fv->eps[0].F[1][1] = (1. - factor);
          fv->eps[0].F[2][2] = 1.;
        }
        /* Simple shear */
        if (fv->eps[0].type == 2){
          fv->eps[0].F[0][0] = 1.;
          fv->eps[0].F[1][1] = 1.;
          fv->eps[0].F[2][2] = 1.;
          fv->eps[0].F[0][1] = factor;
        }
        /* Deviatoric tension */
        if (fv->eps[0].type == 3){
          fv->eps[0].F[0][0] = 1./((1. - factor)*(1. - factor));
          fv->eps[0].F[1][1] = (1. - (NR_t->times[tim] + (sp.step_id+1)*dt)/NR_t->times[nt] * fv->eps[0].load1);
          fv->eps[0].F[2][2] = (1. - (NR_t->times[tim] + (sp.step_id+1)*dt)/NR_t->times[nt] * fv->eps[0].load1);
        }

        if (myrank == 0){
          PGFEM_printf("The deformation gradient F\n");
          for (i=0;i<3;i++){
            for (j=0;j<3;j++){
              PGFEM_printf("%12.12f  ",fv->eps[0].F[i][j]);
            }
            PGFEM_printf("\n");
          }
        }

        /* Residuals */
        nulld (fv->f_u,fv->ndofd);
        INFO = 0;
        *residuals_loc_time += compute_residuals_for_NR(&INFO,grid,mat,fv,sol,load,crpl,com,
                                                        opts,mp,mp_id,t,dts,1);

        for (i=0;i<fv->ndofd;i++){
          fv->f[i] = - fv->f_u[i];
          fv->R[i] = fv->RR[i] = 0.0;
        }

      }/* end periodic */
      else
      {

        for (i=0;i<(load->sups[mp_id])->npd;i++)
          (load->sups[mp_id])->defl_d[i] = load->sup_defl[mp_id][i]/sp.step_size;

        /* Compute macro interface deformation gradient */
        if((load->sups[mp_id])->multi_scale)
        {
          INFO = compute_macro_grad_u((load->sups[mp_id])->F0,load->sups[mp_id],opts->analysis_type);

          if(INFO != 0){
            PGFEM_printerr("[%d] ERROR: not enough prescribed displacements"
                           " for interface multiscale modeling!\n"
                           "Must have at least six (6) prescribed displacements.\n"
                           "Check input and try again.\n",myrank);
            PGFEM_Comm_code_abort(com, 0);
          }
          nulld (fv->f_u,fv->ndofd);
          vol_damage_int_alg(grid->ne,fv->ndofn,fv->d_u,fv->u_np1,grid->element,grid->node,
                             mat->hommat,load->sups[mp_id],dt,iter,com,
                             fv->eps,fv->sig,&max_damage,&dissipation,
                             opts->analysis_type,mp_id);

          INFO = 0;
          *residuals_loc_time += compute_residuals_for_NR(&INFO,grid,mat,fv,sol,load,crpl,com,
                                                          opts,mp,mp_id,t,dts,0);
        } else {
          nulld (fv->f_u,fv->ndofd);
        }

        /*  NODE (PRESCRIBED DEFLECTION)
            - SUPPORT COORDINATES generation of the load vector  */
        nulld (fv->f_defl,fv->ndofd);

        compute_load_vector_for_prescribed_BC(grid,mat,fv,sol,load,dt,crpl,
                                              opts,mp,mp_id,com);
	
        /* Generate the load and vectors */
        for (i=0;i<fv->ndofd;i++)  {
          fv->f[i] = fv->R[i]/sp.step_size - fv->f_defl[i] - fv->f_u[i];
          //sum += f[i];
          fv->RR[i] = fv->RRn[i] + fv->R[i]/sp.step_size*(sp.step_id+1);
        }
      }

      /* Transform LOCAL load vector to GLOBAL */
      LToG (fv->f, fv->BS_f, fv->ndofd, com);

      /* Transform LOCAL load vector to GLOBAL */
      LToG (fv->RR, fv->BS_RR, fv->ndofd, com);

      iter = 0;
      double nor, nor2;
      nor = nor2 = 10.0;

      // Newton Raphson iteration with Line search
      INFO = Newton_Raphson_with_LS(solve_time,stiffmat_loc_time,residuals_loc_time,alpha,
                                    &NOR,&gam,&ART,&iter,grid,mat,fv,sol,load,com,
                                    crpl,opts,mp,mp_id,NR_t->times,dts,
                                    tim,sp.step_size,sp.step_id,&usage,time_steps->tim);

      if(INFO!=0)
      {
        is_sub_cvg = 0;
        break; //goto rest;
      }

      sp.decellerate = gam = ART = 0;

      /* increment the step counter */
      (sol->n_step)++;

      /* /\* turn off line search *\/ */
      /* ART = 1; */

      /* microscale update. Overlay microscale update (using f = r+d_r)
         with macroscale update */
      if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL){
        /* start the microscale jobs */
        MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
        pgf_FE2_macro_client_rebalance_servers(ctx->client,ctx->mscom,
                                               FE2_REBALANCE_NONE);
        unsigned int i;
        double temp_defl[9];
        double newDt;
        double dt0 = NR_t->times[tim+1] - NR_t->times[0];
        int DIV = sp.step_id;
        int STEP = sp.step_size;
        if (DIV == 0) {
          if (tim > 0) {
            ctx->macro->sol->times[ctx->macro->sol->tim - 1]  = ctx->macro->sol->times[ctx->macro->sol->tim];
          }
        }

        if (STEP > 0) {
	  PGFEM_printf("\nIn microscale jobs: NR_t->times[tim+1] - NR_t->times[0] = %e - %e = %e", NR_t->times[tim+1], NR_t->times[0], dt0);
          //const double dtFactor = 1.0e-05;// test ok!
	  double dtFactor = (NR_t->times[tim + 1] - NR_t->times[0]) * (((double) DIV + 1.0)/((double) STEP)); // modified by TP to avoid Segmentation Fault
	  //double dtFactor = (NR_t->times[tim + 1] - ctx->macro->sol->times[ctx->macro->sol->tim - 1]) * (((double) DIV + 1.0)/((double) STEP)); // old: --> Segmentation Fault
          ctx->macro->sol->times[ctx->macro->sol->tim + 1] = dtFactor + NR_t->times[0]; //ctx->macro->sol->times[ctx->macro->sol->tim - 1];
          newDt = ctx->macro->sol->times[ctx->macro->sol->tim - 1] + dtFactor - ctx->macro->sol->times[ctx->macro->sol->tim];
	} 
	else {
          ctx->macro->sol->times[ctx->macro->sol->tim + 1] = ctx->macro->sol->times[ctx->macro->sol->tim] + dt;
          newDt = dt;
        }
        for (i = 0; i < sizeof(load->sup_defl[mp_id]) ; i++) {
          temp_defl[i] = load->sup_defl[mp_id][i];
          load->sup_defl[mp_id][i] = (newDt/(dt0))*load->sup_defl[mp_id][i];
        }
        ctx->macro->sol->times[ctx->macro->sol->tim+1] = NR_t->times[tim] + dt;

        pgf_FE2_macro_client_send_jobs(ctx->client,ctx->mscom,ctx->macro,
                                       JOB_UPDATE);
        for (i = 0; i < sizeof(load->sup_defl[mp_id]) ; i++) {
          load->sup_defl[mp_id][i] = temp_defl[i];
        }
      }

      // update converged values and apply increments
      // for next Newton Raphson step while subdividing

      if(sp.step_size>sp.step_id+1)
      {
        update_values_for_next_NR(grid,mat,fv,sol,load,crpl,com,VVolume,opts,mp,NR_t,mp_id);
        fv->subdivision_factor_n   = fv->subdivision_factor_np1;
      }

      /* finish microscale update */
      if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL)
      {
        /* start the microscale jobs */
        int max_substep = 0;
        MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
        pgf_FE2_macro_client_recv_jobs(ctx->client,ctx->macro,&max_substep);
      }

      /************* TEST THE UPDATE FROM N TO N+1  *************/
      if(NR_UPDATE || PFEM_DEBUG || PFEM_DEBUG_ALL){
        for (i=0;i<fv->ndofd;i++) {fv->f_u[i] = 0.0; fv->d_u[i] = 0.0;}

        if(opts->analysis_type == MINI){
          MINI_check_resid(fv->ndofn,grid->ne,grid->element,grid->node,mat->hommat,fv->eps,
                           fv->sig,fv->d_u,load->sups[mp_id],fv->RR,com->DomDof,fv->ndofd,
                           com,mp_id);
        }
        if(opts->analysis_type == MINI_3F){
          MINI_3f_check_resid(fv->ndofn,grid->ne,grid->element,grid->node,mat->hommat,fv->eps,
                              fv->sig,fv->d_u,load->sups[mp_id],fv->RR,com->DomDof,fv->ndofd,
                              com,mp_id);
        }
        INFO = 0;
        *residuals_loc_time += compute_residuals_for_NR(&INFO,grid,mat,fv,sol,load,crpl,com,
                                                        opts,mp,mp_id,t,dts,0);

        for (i=0;i<fv->ndofd;i++) fv->f[i] = fv->RR[i] - fv->f_u[i];
        /* print_array_d(stdout,RR,ndofd,1,ndofd); */

        LToG(fv->f, fv->BS_f, fv->ndofd, com);
        nor = ss(fv->BS_f,fv->BS_f,com->DomDof[myrank]);
        double tmp;
        com->net->allreduce(&nor,&tmp,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
        nor = sqrt (tmp);

        if (myrank == 0) PGFEM_printf("NORM NORM = %12.12e\n",nor);

      }
      /************* TEST THE UPDATE FROM N TO N+1  *************/

      sp.step_id++;
      if(NR_PRINT_INTERMEDIATE){
        if(sp.step_size > sp.step_id){
          /* converged intermediate step, but not last step print output */
          char fname[100];
          sprintf(fname,"%s_%ld",opts->ofname,tim);
          if(myrank == 0){
            VTK_print_master(opts->opath,fname,sol->n_step,com->nproc,opts);
          }
          VTK_print_vtu(opts->opath,fname,sol->n_step,myrank,grid->ne,grid->nn,grid->node,
                        grid->element,load->sups[mp_id],fv->u_np1,fv->tf.P_np1.m_pdata,fv->tf.V_np1.m_pdata,fv->sig,fv->eps,opts,mp_id);
        }
      }

      double *res_trac = aloc1(3);
      bounding_element_compute_resulting_traction(grid->n_be,grid->b_elems,grid->element,grid->node,
                                                  fv->eps,fv->sig,fv->ndofd,com->DomDof,
                                                  com,
                                                  opts->analysis_type,
                                                  res_trac);
      PGFEM_free(res_trac);

      if(sp.step_size > 2 && sp.step_id == 2)
      {
        is_sub_cvg = 0;
        break; //goto rest;
      }

    } //end SUBDIVISION
    if(is_sub_cvg)
    {
      sp.accellerate = 0;

      if(NR_COMPUTE_REACTIONS && !(load->sups[mp_id])->multi_scale){
        compute_reactions(grid->ne,fv->ndofn,fv->npres,fv->u_np1,grid->node,grid->element,mat->matgeom,
                          mat->hommat,load->sups[mp_id],fv->eps,fv->sig,sol->nor_min,crpl,
                          dt,opts->stab,com,opts->analysis_type,mp_id);
      }
      break;
    }
  }
}


/// Computing residuals for dependent physics when multiphysics is staggered
/// needs to use temporal field variables. This function serves to compute
/// residulas with temporal field variables.
///
/// This will not perform any integration algorithm.
///
/// \parma[out] residuals_loc_time - compute time for constructing b in Ax=b
/// \param[in] nor norm of residuals
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] t time at t(n+1)
/// \param[in] dts time step sizes a n, and n+1
/// \param[in] mp_id mutiphysics id
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int compute_coupled_physics_residual_norm(double *residuals_loc_time,
                                          double *nor,
                                          Grid *grid,
                                          MaterialProperty *mat,
                                          FieldVariables *FV,
                                          Solver *SOL,
                                          LoadingSteps *load,
                                          CommunicationStructure *COM,
                                          TimeStepping *time_steps,
                                          CRPL *crpl,
                                          const PGFem3D_opt *opts,
                                          const Multiphysics& mp,
                                          double t,
                                          double *dts,
                                          int mp_id)
{
  int err = 0;
  // use pointers for physics[mp_id]
  Solver                 *sol = SOL + mp_id;
  FieldVariables          *fv = FV  + mp_id;
  CommunicationStructure *com = COM + mp_id;

  // temporal
  // double *u_n   = fv->u_n;
  // double *u_nm1 = fv->u_nm1;
  // State_variables *statv_list = fv->statv_list;

  // double dt = dts[DT_NP1];

  int myrank = com->rank;
  
  for(int ia=0; ia<fv->ndofd; ia++)
    fv->f_u[ia] = 0.0;

  sol->run_integration_algorithm = 0; // turn off running integration algorithm

  long INFO = 0;
  *residuals_loc_time += compute_residuals_for_NR(&INFO,grid,mat,fv,sol,load,crpl,com,
                                                  opts,mp, mp_id,t,dts, 1);
  
  sol->run_integration_algorithm = 1; // reset integration algorithm to be active

  // Transform LOCAL load vector to GLOBAL
  LToG(fv->f_u, fv->BS_f_u, fv->ndofd, com);

  // Compute Euclidian norm
  for(int i=0;i<com->DomDof[myrank]; i++)
    fv->BS_f[i] = fv->BS_RR[i] - fv->BS_f_u[i];

  double LNOR = 0.0;
  double GNOR = 0.0;

  LNOR = ss(fv->BS_f,fv->BS_f,com->DomDof[myrank]);
  com->net->allreduce(&LNOR,&GNOR,1,NET_DT_DOUBLE,NET_OP_SUM,com->comm);
  *nor = sqrt(GNOR);

  return err;
}

/// save variables(t(n) and t(n-1)) to temporal_variables(t(n) and t(n-1))
///
/// \param[in] grid a mesh object
/// \param[in,out] FV array of field variable object
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int save_field_variables_to_temporal(Grid *grid,
                                     FieldVariables *FV,
                                     const PGFem3D_opt *opts,
                                     const Multiphysics& mp,
                                     int mp_id)
{
  FieldVariables *fv = FV + mp_id;

  int err = 0;
  for(int ia=0; ia<(grid->nn)*(fv->ndofn); ia++)
  {
    fv->temporal->u_nm1[ia] = fv->u_nm1[ia];
    fv->temporal->u_n[ia]   = fv->u_n[ia];
  }

  if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL && (opts->analysis_type==CM || opts->analysis_type==CM3F))
    err += constitutive_model_save_state_vars_to_temporal(FV + mp_id, grid);
    
  if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL && opts->analysis_type==CM3F)
  {
    for(int ia=0; ia<(grid->ne)*(fv->nVol); ia++)
    {
      fv->temporal->tf.V_nm1.m_pdata[ia] = fv->tf.V_nm1.m_pdata[ia];
      fv->temporal->tf.V_n.m_pdata[ia]   = fv->tf.V_n.m_pdata[ia];      
    }
    
    for(int ia=0; ia<(grid->ne)*(fv->npres); ia++)
    {
      fv->temporal->tf.P_nm1.m_pdata[ia] = fv->tf.P_nm1.m_pdata[ia];
      fv->temporal->tf.P_n.m_pdata[ia]   = fv->tf.P_n.m_pdata[ia];
    }    
  }    

  return err;
}

/// update variables(t(n+1)) to temporal_variables(t(n+1))
///
/// \param[in] grid a mesh object
/// \param[in,out] FV array of field variable object
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int update_temporal_field_variables_np1(Grid *grid,
                                        FieldVariables *FV,
                                        const PGFem3D_opt *opts,
                                        const Multiphysics& mp,
                                        int mp_id)
{
  int err = 0;
  if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL && (opts->analysis_type==CM || opts->analysis_type==CM3F))
    err += constitutive_model_update_np1_state_vars_to_temporal(FV + mp_id, grid);

  return err;
}

/// reset variables(t(n) and t(n-1)) from temporal_variables(t(n) and t(n-1))
///
/// \param[in] grid a mesh object
/// \param[in,out] FV array of field variable object
/// \param[in] opts structure of PGFem3D options
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int reset_field_variables_using_temporal(Grid *grid,
                                         FieldVariables *FV,
                                         const PGFem3D_opt *opts,
                                         const Multiphysics& mp,
                                         int mp_id)
{
  FieldVariables *fv = FV + mp_id;

  int err = 0;
  for(int ia=0; ia<(grid->nn)*(fv->ndofn); ia++)
  {
    fv->u_nm1[ia] = fv->temporal->u_nm1[ia];
    fv->u_n[ia]   = fv->temporal->u_n[ia];
  }

  if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL && 
    (opts->analysis_type==CM || opts->analysis_type==CM3F))
    err += constitutive_model_reset_state_using_temporal(FV + mp_id, grid);
    
  if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL && opts->analysis_type==CM3F)
  {
    for(int ia=0; ia<(grid->ne)*(fv->nVol); ia++)
    {
      fv->tf.V_nm1.m_pdata[ia] = fv->temporal->tf.V_nm1.m_pdata[ia];
      fv->tf.V_n.m_pdata[ia] = fv->temporal->tf.V_n.m_pdata[ia];      
    }
    
    for(int ia=0; ia<(grid->ne)*(fv->npres); ia++)
    {
      fv->tf.P_nm1.m_pdata[ia] = fv->temporal->tf.P_nm1.m_pdata[ia];
      fv->tf.P_n.m_pdata[ia] = fv->temporal->tf.P_n.m_pdata[ia];      
    }    
  }

  return err;
}

/// reset mesh coordinate (t(n)) from using saved displacement (t(n))
///
/// \param[in] grid a mesh object
/// \param[in,out] FV array of field variable object
/// \param[in] opts structure of PGFem3D options
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int reset_mesh_coordinates(Grid *grid,
                           FieldVariables *FV,
                           const PGFem3D_opt *opts,
                           const Multiphysics& mp,
                           int mp_id)
{
  FieldVariables *fv = FV + mp_id;


  int err = 0;
  for(int ia=0; ia<grid->nn; ia++)
  {
    int id = ia*(fv->ndofn);
    grid->node[ia].x1 = grid->node[ia].x1_fd + fv->temporal->u_n[id+0];
    grid->node[ia].x2 = grid->node[ia].x2_fd + fv->temporal->u_n[id+1];
    grid->node[ia].x3 = grid->node[ia].x3_fd + fv->temporal->u_n[id+2];
  }
  return err;
}

/// reset variables(t(n) and t(n-1)) from temporal_variables(t(n) and t(n-1))
///
/// \param[in] grid a mesh object
/// \param[in,out] FV array of field variable object
/// \param[in] SOL object array for solution scheme
/// \param[in] opts structure of PGFem3D options
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int reset_state_field_variables(Grid *grid,
                                FieldVariables *FV,
                                Solver *SOL,
                                const PGFem3D_opt *opts,
                                const Multiphysics& mp,
                                int mp_id)
{
  int err = 0;

  if(mp.physics_ids[mp_id]==MULTIPHYSICS_MECHANICAL)
    err += reset_mesh_coordinates(grid,FV,opts,mp,mp_id);

  if(SOL[mp_id].is_subdivided)
    err += reset_field_variables_using_temporal(grid,FV,opts,mp,mp_id);

  return err;
}

/// set times for Netwon Raphson iterations
///
/// in Newton_Raphson iteration tim = 1 is the current time step id
/// where times[0] = t(n-1)
///       times[1] = t(n)
///       times[2] = t(n+1)
/// s.t dt(n)   = times[1] - times[0]
///     dt(n+1) = times[2] - times[1]
/// HOWEVER when actual time step is at 0 (time_steps->tim == 0)
/// tim = 0, because there is no t(n-1) values
///
/// \param[in] ts object for time stepping
/// \param[out] NR_t->times store times, NR_t->times[0]   = t(n-1)
///                                      NR_t->times[1]   = t(n)
///                                      NR_t->times[2]   = t(n+1)
///             NR_t->dt store dt,       NR_t->dt[DT_N]   = t(n)   - t(n-1)
///                                      NR_t->dt[DT_NP1] = t(n+1) - t(n)
/// \return non-zero on internal error
int set_time_step_info_for_NR(TimeStepping *ts,
                              NR_time_steps *NR_t)
{
  int err = 0;

  long t_step_id = ts->tim;

  if(ts->tim==0)
  {
    NR_t->tim = 0;
    NR_t->times[0] = ts->times[t_step_id];
    NR_t->times[1] = ts->times[t_step_id+1];
    NR_t->times[2] = 2.0*NR_t->times[1] + NR_t->times[0]; // this is dummy
  }
  else
  {
    NR_t->tim = 1;
    NR_t->times[0] = ts->times[t_step_id-1];
    NR_t->times[1] = ts->times[t_step_id];
    NR_t->times[2] = ts->times[t_step_id+1];
  }
  NR_t->dt[DT_N]   = NR_t->times[1] - NR_t->times[0];
  NR_t->dt[DT_NP1] = NR_t->times[2] - NR_t->times[1];

  return err;
}

/// After complete Newton Raphson for physics X, check residuals of physics Y, and Z
/// that are dependent on physics X.
///
/// \parma[out] residuals_loc_time - compute time for constructing b in Ax=b
/// \param[out] is_cnvged if any residuls of the dependent physics is greater than tolerance,
///                       it is updated to zero
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in] FV array of field variable object
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] NR_time save time step info during previous NR iteration
/// \param[in] mp_id mutiphysics id
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int check_convergence_of_NR_staggering(double *residuals_loc_time,
                                       int *is_cnvged,
                                       Grid *grid,
                                       MaterialProperty *mat,
                                       FieldVariables *FV,
                                       Solver *SOL,
                                       LoadingSteps *load,
				       CommunicationStructure *COM,
                                       TimeStepping *time_steps,
                                       CRPL *crpl,
                                       const PGFem3D_opt *opts,
                                       const Multiphysics& mp,
                                       NR_time_steps *NR_time,
                                       int mp_id)
{
  int err = 0;

  CommunicationStructure *com = COM + mp_id;
  int myrank = com->rank;
  
  // @todo This variable is set but not used ever. I commented it out as dead
  //       code. This should be reviewed by @cp and eliminated if
  //       appropriate. LD
  //
  // int tim = 1;
  // if( time_steps->tim==0)
  //   tim = 0;

  for(int ib = 0; ib<FV[mp_id].n_coupled; ib++)
  {
    int cpled_mp_id = mp.coupled_ids[mp_id][ib+1];
    int is_it_dependent_on_physics_mp_id = 0;
    for(int ic=0; ic<FV[cpled_mp_id].n_coupled; ic++)
    {
      if(mp.coupled_ids[cpled_mp_id][ic+1]==mp_id)
      {
        // phsics[cpled_mp_id] is dependent on phsics[mp_id]
        is_it_dependent_on_physics_mp_id = 1;
        break;
      }
    }

    if(is_it_dependent_on_physics_mp_id==0)
      continue;

    double nor = 0.0;

    compute_coupled_physics_residual_norm(residuals_loc_time,&nor,grid,mat,FV,SOL,load,
                                          COM,time_steps,crpl,opts,mp,
                                          NR_time[cpled_mp_id].times[2],
                                          NR_time[cpled_mp_id].dt,cpled_mp_id);

    double Rn_R = fabs(SOL[cpled_mp_id].last_residual - nor)/FV[cpled_mp_id].NORM;

    long global_ODE_EXA_metric = 0;
    com->net->reduce(&perIter_ODE_EXA_metric,&global_ODE_EXA_metric,1,NET_DT_LONG,NET_OP_SUM,0,com->comm);
    if(myrank==0)
    {
      PGFEM_printf(":: R(%s): ", mp.physicsname[cpled_mp_id]);
      PGFEM_printf("|R| = %e |R|/|R0| = %e, Rn = %e, Rn-R = %e\n",
             nor, nor/FV[cpled_mp_id].NORM, SOL[cpled_mp_id].last_residual,Rn_R);
             
      //EXA metric computation
      if (opts->print_EXA == 2)
        PGFEM_printf ("Convergence ODE operations: (%ld)\n", global_ODE_EXA_metric);
      ODE_EXA_metric[mp_id] += global_ODE_EXA_metric;          //accumulate exascale metrics for the convergence check
    }
    perIter_ODE_EXA_metric = 0;
    
    
    if(nor/FV[cpled_mp_id].NORM>SOL[cpled_mp_id].nor_min && Rn_R > SOL[cpled_mp_id].nor_min)
      *is_cnvged = 0;
  }
  return err;
}

/// set first residual for very stiff problem if needed
///
/// Pre-computing residual helps to converge the solution when
/// iteration step is very stiff. This function computes residual and set
/// as first NORM at the very first time of the non-linear iteration step by
/// perturbing field variables. (delta u is read from solver file and used to
/// perturbate the initial values (fv->u0)) such that fv->u0 + sol->du is used to
/// compute the residual.
///
/// This will not perform any integration algorithm.
///
/// \parma[out] residuals_loc_time - compute time for constructing b in Ax=b
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object. FV.NORM is updateed.
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] com->comm COM->COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int set_0th_residual(std::vector<double> &residuals_time,
                     Grid *grid,
                     MaterialProperty *mat,
                     FieldVariables *FV,
                     Solver *SOL,
                     LoadingSteps *load,
		     CommunicationStructure *COM,
                     TimeStepping *time_steps,
                     CRPL *crpl,
                     const PGFem3D_opt *opts,
                     const Multiphysics& mp)
{
  int err = 0;
  int myrank = COM[0].rank;
  
  if(time_steps->tim>0)
    return err;

  for(int mp_id=0; mp_id<mp.physicsno; mp_id++)
  {
    if(SOL[mp_id].set_initial_residual==0)
      continue;

    for(int ia=0; ia<FV[mp_id].ndofd; ia++)
      FV[mp_id].u_np1[ia] = FV[mp_id].u0 + SOL[mp_id].du;

    double nor = 0.0;
    double residuals_loc_time = 0.0;
    NR_time_steps NR_t;
    set_time_step_info_for_NR(time_steps,&NR_t);
    compute_coupled_physics_residual_norm(&residuals_loc_time,&nor,grid,mat,FV,SOL,load,COM,
                                          time_steps,crpl,opts,mp,NR_t.times[2],
                                          NR_t.dt,mp_id);

    residuals_time[mp_id] += residuals_loc_time; 
    FV[mp_id].NORM = nor; // set first residual

    for(int ia=0; ia<FV[mp_id].ndofd; ia++)
      FV[mp_id].u_np1[ia] = FV[mp_id].u0;

    if(myrank==0)
    {
      PGFEM_printf("INFO. The first residual for the physics, %s, is computed as %e\n", mp.physicsname[mp_id], nor);
      PGFEM_printf("by perturbing the disp. with %e\n", SOL[mp_id].du);
    }
  }
  return err;
}

/// Staggered Newton Raphson iterative solver for multiphysics problems
///
/// \parma[out] solve_time time measure of spent for linear solver
/// \parma[out] stiffmat_loc_time - compute time for constructing A in Ax=b  
/// \parma[out] residuals_loc_time - compute time for constructing b in Ax=b
/// \param[out] iterno number of staggered iterations taken in this routine
/// \param[out] is_SNR_converged if 1, the SNR converged
/// \param[out] alpha_out physics based evolution rate
/// \param[out] NR_t_in container of time stepping info
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] VVolume original volume of the domain
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
void Multiphysics_Newton_Raphson_sub(std::vector<double> &hypre_time,
                                     std::vector<double> &stiffmat_time,
                                     std::vector<double> &residuals_time,
                                     int *iterno,
                                     int *is_SNR_converged,
                                     double *alpha_out,
                                     NR_time_steps *NR_t_in,
                                     Grid *grid,
                                     MaterialProperty *mat,
                                     FieldVariables *FV,
                                     Solver *SOL,
                                     LoadingSteps *load,
                                     CommunicationStructure *COM,
                                     TimeStepping *time_steps,
                                     CRPL *crpl,
                                     const double VVolume,
                                     const PGFem3D_opt *opts,
                                     const Multiphysics& mp)
{
  
  const int print_level = 1;
  *iterno = 0;
  *is_SNR_converged = 1;
  *alpha_out = 0.0;
  double alpha_single = 0.0; // measure of physics based evolution rate for non-coupled physics

  int max_itr = SOL[0].max_NR_staggering; // max. number staggering
  int myrank = COM[0].rank;

  NR_time_steps *NR_time = (NR_time_steps *) malloc(sizeof(NR_time_steps)*mp.physicsno);

  int coupled_physics_no = 0;
  int printed_before = 0;

  for(int mp_id=0; mp_id<mp.physicsno; mp_id++)
  {
    save_field_variables_to_temporal(grid,FV,opts,mp,mp_id);

    NR_time[mp_id].times[0] = NR_t_in->times[0];
    NR_time[mp_id].times[1] = NR_t_in->times[1];
    NR_time[mp_id].times[2] = NR_t_in->times[2];

    NR_time[mp_id].dt[DT_N]   = NR_t_in->dt[DT_N];
    NR_time[mp_id].dt[DT_NP1] = NR_t_in->dt[DT_NP1];
    NR_time[mp_id].tim        = NR_t_in->tim;

    if(FV[mp_id].n_coupled>0)
    {
      coupled_physics_no++;
      continue;
    }
    //print current physics name

    if(!printed_before)
    {
      if(myrank==0)
      {
        PGFEM_printf("%s\n", line_level_1);
        PGFEM_printf(":: Do Newton Raphson for independent physics first\n");
      }
      printed_before = 1;
    }

    if(myrank==0)
    {
      PGFEM_printf("%s\n", line_level_2);
      PGFEM_printf(":: Newton Raphson iteration for : %s\n", mp.physicsname[mp_id]);
    }

    double alpha = 0.0;
    double solve_time = 0.0;
    double stiffmat_loc_time = 0.0;
    double residuals_loc_time = 0.0;
    perform_Newton_Raphson_with_subdivision(&solve_time,&stiffmat_loc_time,&residuals_loc_time,
                                            print_level,is_SNR_converged,&alpha,grid,mat,FV,
                                            SOL,load,COM,time_steps,crpl,VVolume,opts,
                                            mp,NR_time+mp_id,mp_id);

    hypre_time[mp_id] += solve_time;
    stiffmat_time[mp_id] += stiffmat_loc_time;
    residuals_time[mp_id] += residuals_loc_time;

    alpha_single = (alpha_single > alpha)? alpha_single: alpha;

    if(SOL[mp_id].is_subdivided)
    {
      reset_state_field_variables(grid,FV,SOL,opts,mp,mp_id);
      if(myrank==0)
        PGFEM_printf("NR is subdivied, reset t(n)\n");
    }

    if(*is_SNR_converged==0)
      break;

    update_temporal_field_variables_np1(grid,FV,opts,mp,mp_id);
  }

  // if there is no physics coupled with others, Newton_Raphson is done here.
  if(coupled_physics_no>0 && *is_SNR_converged)
  {
    //print Staggered Newton Raphson iteration
    if(myrank==0)
    {
      PGFEM_printf("%s\n", line_level_1);
      PGFEM_printf(":: Staggered Newton Raphson iteration starts for dependent physics \n");
    }

    // Prepare memory for saving load data if subdivided
    // Load will be divided, and applied subsequently while subdivision is stepping.
    double **sup_defl = (double **) malloc(sizeof(double *)*mp.physicsno);
    double **defl     = (double **) malloc(sizeof(double *)*mp.physicsno);
    double **R        = (double **) malloc(sizeof(double *)*mp.physicsno);
    double **RRn      = (double **) malloc(sizeof(double *)*mp.physicsno);

    // start save data
    for(int ia=0; ia<mp.physicsno; ia++)
    {
      int npd = (load->sups[ia])->npd;
      if(npd>0)
      {
        sup_defl[ia] = (double *) malloc(sizeof(double)*npd);
        defl[ia] = (double *) malloc(sizeof(double)*npd);

        for(int ib=0;ib<npd;ib++)
        {
          sup_defl[ia][ib] = load->sup_defl[ia][ib];
          defl[ia][ib] = (load->sups[ia])->defl[ib];
        }
      }
      else
      {
	sup_defl[ia] = NULL;
	defl[ia] = NULL;
      }

      R[ia] = (double *) malloc(sizeof(double)*FV[ia].ndofd);
      RRn[ia] = (double *) malloc(sizeof(double)*FV[ia].ndofd);
      for(int ib=0; ib<FV[ia].ndofd; ib++)
      {
        R[ia][ib] = FV[ia].R[ib];
        RRn[ia][ib] = FV[ia].RRn[ib];
      }
    }

    while(*iterno<max_itr)
    {
      int is_cnvged = 1;

      double alpha_cpled = 0.0; // measure of physics based evolution rate for non-coupled physics
      for(int mp_id=0; mp_id<mp.physicsno; mp_id++)
      {
        if(FV[mp_id].n_coupled >0)
        {
          //print current physics name
          if(myrank==0)
          {
            PGFEM_printf("%s\n", line_level_2);
            PGFEM_printf(":: (%d/%d) Newton Raphson iteration for : %s\n", *iterno, max_itr, mp.physicsname[mp_id]);
          }

          if(*iterno>0)
          {
            // reset time step size
            NR_time[mp_id].times[0] = NR_t_in->times[0];
            NR_time[mp_id].times[1] = NR_t_in->times[1];
            NR_time[mp_id].times[2] = NR_t_in->times[2];

            NR_time[mp_id].dt[DT_N]   = NR_t_in->dt[DT_N];
            NR_time[mp_id].dt[DT_NP1] = NR_t_in->dt[DT_NP1];
            NR_time[mp_id].tim        = NR_t_in->tim;

            int npd = (load->sups[mp_id])->npd;
            if(npd>0)
            {
              for(int ib=0;ib<npd;ib++)
              {
                load->sup_defl[mp_id][ib] = load->sups[mp_id]->defl_d[ib] = sup_defl[mp_id][ib];
                (load->sups[mp_id])->defl[ib] = defl[mp_id][ib];
              }
            }

            for(int ib=0; ib<FV[mp_id].ndofd; ib++)
            {
              FV[mp_id].R[ib]   =   R[mp_id][ib];
              FV[mp_id].RRn[ib] = RRn[mp_id][ib];

              FV[mp_id].d_u[ib]    = 0.0;
              FV[mp_id].dd_u[ib]   = 0.0;
              FV[mp_id].f_defl[ib] = 0.0;
              FV[mp_id].f[ib]      = 0.0;
            }
            
            if((mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL) && 
              (opts->analysis_type == TF || opts->analysis_type == CM3F))
              FV[mp_id].tf.reset();
          }

          SOL[mp_id].is_subdivided = 0;

          double alpha = 0.0;
          double solve_time = 0.0;
          double stiffmat_loc_time = 0.0;
          double residuals_loc_time = 0.0;
          perform_Newton_Raphson_with_subdivision(&solve_time,&stiffmat_loc_time,&residuals_loc_time,
                                                  print_level,is_SNR_converged,&alpha,grid,mat,FV,
                                                  SOL,load,COM,time_steps,crpl,VVolume,
                                                  opts,mp,NR_time+mp_id,mp_id);
	  
          hypre_time[mp_id] += solve_time;
          stiffmat_time[mp_id] += stiffmat_loc_time;
          residuals_time[mp_id] += residuals_loc_time;

          alpha_cpled = (alpha_cpled>alpha)? alpha_cpled: alpha;

          if(*is_SNR_converged==0)
            break;

          update_temporal_field_variables_np1(grid,FV,opts,mp,mp_id);

          residuals_loc_time = 0.0;
          check_convergence_of_NR_staggering(&residuals_loc_time,&is_cnvged,grid,mat,FV,SOL,load,
                                             COM,time_steps,crpl,opts,mp,NR_time,mp_id);

          residuals_time[mp_id] += residuals_loc_time;
        }
      }

      *alpha_out = (alpha_cpled>alpha_single)? alpha_cpled: alpha_single;

      // allways go back to t(n) for next NR
      for(int mp_id=0; mp_id<mp.physicsno; mp_id++)
      {
        if(FV[mp_id].n_coupled == 0)
          continue;

        if(SOL[mp_id].is_subdivided)
        {
          reset_state_field_variables(grid,FV,SOL,opts,mp,mp_id);
          if(myrank==0)
            PGFEM_printf("NR for %s is subdivied, reset t(n)\n", mp.physicsname[mp_id]);
        }
      }

      if(*is_SNR_converged)
      {
        if(is_cnvged==0 && (*iterno) == (max_itr-1)) // last step, but not converged
        {
          *is_SNR_converged = 0;
          break;
        }

        if(is_cnvged)
        {
          if(myrank==0)
          {
            PGFEM_printf("%s\n", line_level_1);
            PGFEM_printf(":: Staggered Newton Raphson is converged. Done.\n");
          }
          break;
        }
      }
      else
        break;

      (*iterno)++;
    }

    // if not converged, reset accumulated prescribed boundary values which have been updated during the NR iterations.
    // free allocated memory
    for(int mp_id=0; mp_id<mp.physicsno; mp_id++)
    {
      int npd = (load->sups[mp_id])->npd;
      if(npd>0)
      {
        if(*is_SNR_converged == 0)
        {
          for(int ib=0;ib<npd;ib++){
	    assert(defl[mp_id] != NULL && "defl[mp_id] can't be null");
            (load->sups[mp_id])->defl[ib] = defl[mp_id][ib];
	  }
        }

        free(sup_defl[mp_id]);
        free(defl[mp_id]);
      }

      free(R[mp_id]);
      free(RRn[mp_id]);
    }
    free(sup_defl);
    free(defl);
    free(R);
    free(RRn);
  }

  if(*is_SNR_converged)
  {
    // update final results
    for(int mp_id=0; mp_id<mp.physicsno; mp_id++)
    {
      int tim = NR_time[mp_id].tim;
      time_steps->tns[mp_id] = NR_time[mp_id].times[tim+1] - NR_time[mp_id].dt[DT_NP1];
      update_values_for_next_NR(grid,mat,FV+mp_id,SOL+mp_id,load,
                                crpl,COM+mp_id,VVolume,opts,mp,NR_time+mp_id,mp_id);
    }
  }

  free(NR_time);
}

/// Staggered Newton Raphson iterative solver for multiphysics problems
///
/// If staggered NR step diverges, time step will be subdivided.
///
/// \parma[out] solve_time time measure of spent for linear solver
/// \parma[out] stiffmat_loc_time - compute time for constructing A in Ax=b  
/// \parma[out] residuals_loc_time - compute time for constructing b in Ax=b
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] VVolume original volume of the domain
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
void Multiphysics_Newton_Raphson(std::vector<double> &hypre_time,
                                 std::vector<double> &stiffmat_time,
                                 std::vector<double> &residuals_time,
                                 Grid *grid,
                                 MaterialProperty *mat,
                                 FieldVariables *FV,
                                 Solver *SOL,
                                 LoadingSteps *load,
                                 CommunicationStructure *COM,
                                 TimeStepping *time_steps,
                                 CRPL *crpl,
                                 const double VVolume,
                                 const PGFem3D_opt *opts,
                                 const Multiphysics& mp)
{
  const int print_level = 1;
  double alpha = 0.0;
  int iterno = 0;

  int myrank = COM[0].rank;

  // if too siff to converge, try pre compute resiudal by perturbing displacement
  // slightly
  if(time_steps->tim==0)
  {
    set_0th_residual(residuals_time,grid,mat,FV,SOL,load,COM,time_steps,crpl,
                     opts,mp);
  }
  perIter_ODE_EXA_metric = 0;

  NR_time_steps NR_t;
  NR_time_steps NR_t_sub;
  set_time_step_info_for_NR(time_steps,&NR_t);
  set_time_step_info_for_NR(time_steps,&NR_t_sub);

  assert(mp.physicsno > 0 && "Number of physics must be positive");

  // if single physics
  //------------------------------------------------------------------------------------------
  if(mp.physicsno==1)
  {

    int is_NR_cvg = 1;
    double solve_time = 0.0;
    double stiffmat_loc_time = 0.0;
    double residuals_loc_time = 0.0;
    perform_Newton_Raphson_with_subdivision(&solve_time,&stiffmat_loc_time,&residuals_loc_time,
                                            print_level,&is_NR_cvg,&alpha,grid,mat,FV,SOL,load,
                                            COM,time_steps,crpl,VVolume,opts,mp,&NR_t,0);
    
    hypre_time[0] += solve_time;
    stiffmat_time[0] += stiffmat_loc_time;
    residuals_time[0] += residuals_loc_time;

    update_values_for_next_NR(grid,mat,FV,SOL,load,crpl,COM,VVolume,
                              opts,mp,&NR_t,0);
    time_steps->times[time_steps->tim] = time_steps->times[time_steps->tim+1] - NR_t.dt[DT_NP1];
  }
  // single physics is done


  // if multiphysics
  //------------------------------------------------------------------------------------------
  else
  {
    // Prepare memory for saving load data if subdivided
    // Load will be divided, and applied subsequently while subdivision is stepping.
    vector<vector<double>> sup_defl(mp.physicsno), R(mp.physicsno), RRn(mp.physicsno);
  
    gcm::Matrix<bool> apply_V0;
    if(time_steps->tim==0)    
      apply_V0.initialization(mp.physicsno, 1, false);
    
    // start save data
    for(int ia=0; ia<mp.physicsno; ia++)
    {
      if(time_steps->tim==0) 
        apply_V0(ia) = FV[ia].apply_initial_velocity;
        
      int npd = (load->sups[ia])->npd;
      if(npd>0)
      {
        sup_defl[ia].resize(npd);
  
        for(int ib=0;ib<npd;ib++)
          sup_defl[ia][ib] = load->sup_defl[ia][ib];
      }
      if (FV[ia].ndofd > 0){
        R[ia].resize(FV[ia].ndofd);
        RRn[ia].resize(FV[ia].ndofd);
        for(int ib=0; ib<FV[ia].ndofd; ib++)
        {
  	  R[ia][ib] = FV[ia].R[ib];
  	  RRn[ia][ib] = FV[ia].RRn[ib];
        }
      }
    } //end mp.physicsno loop
  
    SubdivisionScheme sp(opts->solution_scheme_opt[NO_SUBDIVISION_LIMITS]);  
    int INFO = 0; // if INF0 == 0, no error detected
                  // if INFO > 0 , error detected
  
  
    int iter_max = SOL[0].max_NR_staggering;
    double tn = NR_t.times[NR_t.tim];
  
    while(1)
    {
      // run subdivision
      sp.do_subdivision(INFO,NR_t.dt + DT_NP1,NR_t.times,NR_t.tim,iterno,iter_max,alpha,COM);
  
      for(int ia=0; ia<mp.physicsno; ia++){
        update_load_increments_for_subdivision(&sp,sup_defl[ia],(load->sups[ia])->npd,RRn[ia],R[ia],FV[ia].ndofd);
      }
  
      INFO = 0;
      int is_sub_cvg = 1;
  
      while(sp.step_size > sp.step_id)
      {
        for(int ia=0; ia<mp.physicsno; ia++)
        {
          for (int ib=0;ib<(load->sups[ia])->npd;ib++)
          {
            load->sup_defl[ia][ib]       = sup_defl[ia][ib]/sp.step_size;
            (load->sups[ia])->defl_d[ib] = sup_defl[ia][ib]/sp.step_size;
          }
  	  
          for(int ib=0; ib<FV[ia].ndofd; ib++)
          {
            FV[ia].R[ib]   = R[ia][ib]/sp.step_size;
            FV[ia].RRn[ib] = RRn[ia][ib] + R[ia][ib]/sp.step_size*(sp.step_id+1);
          }
        }
  
        int is_sub_converged = 0;
  
        // 1st level of the subdivision has own time stepping
        // and 2nd level can have its subdivision to. So, prepare time stepping for sub level (2nd level)
        if(NR_t_sub.tim==0)
        {
          NR_t_sub.times[1] = tn + NR_t.dt[DT_NP1];
          NR_t_sub.times[2] = tn + 2.0*NR_t.dt[DT_NP1];
        }
        else
        {
          NR_t_sub.times[0] = tn - NR_t.dt[DT_N];
          NR_t_sub.times[1] = tn;
          NR_t_sub.times[2] = tn + NR_t.dt[DT_NP1];
        }
  
        NR_t_sub.dt[DT_N]   = NR_t.dt[DT_N];
        NR_t_sub.dt[DT_NP1] = NR_t.dt[DT_NP1];
  
        // print status
        if(sp.is_subdivided && myrank == 0)
        {
          if(myrank==0)
          {
            PGFEM_printf("%s\n", line_level_0);
            PGFEM_printf(":: (%ld:%d/%d) Muliphysics Newton Raphson \n", time_steps->tim,sp.step_id, sp.step_size);
            PGFEM_printf(":: time = (%e, %e, %e)\n", NR_t_sub.times[0], NR_t_sub.times[1], NR_t_sub.times[2]);
          }
        }
      
        // perform sub level (Staggered Newton iteration)
        alpha = 0.0; // always reset evolution threshold
        
        if(time_steps->tim==0){          
          for(int ia=0; ia<mp.physicsno; ia++)
            FV[ia].apply_initial_velocity = apply_V0(ia);
        }
       
        Multiphysics_Newton_Raphson_sub(hypre_time,stiffmat_time,residuals_time,&iterno,
                                        &is_sub_converged,&alpha,&NR_t_sub,
                                        grid,mat,FV,SOL,load,COM,time_steps,crpl,
                                        VVolume,opts,mp);
   
        if(myrank==0)
          PGFEM_printf(":: Maximum physics based evolution threshold = %f\n", alpha);
  
        // check convergence
        if(!is_sub_converged)
        {
          INFO = 1;
          is_sub_cvg = 0;
          break;
        }
        sp.decellerate = 0;
        sp.step_id++;
  
        tn += NR_t.dt[DT_NP1];
        NR_t.dt[DT_N] = NR_t.dt[DT_NP1];
  
        if(NR_t_sub.tim==0)
        {
          NR_t_sub.tim=1;
          if(time_steps->tim==0){          
            for(int ia=0; ia<mp.physicsno; ia++)
              apply_V0(ia) = false;
          }
        }
        if(sp.step_size > 2 && sp.step_id == 2)
        {
          is_sub_cvg = 0;
          break; //goto rest;
        }
      } //end while(sp.step_size > sp.step_id)
      if(is_sub_cvg)
      {
        NR_t.dt[DT_N] = NR_t.dt[DT_NP1];
        sp.accellerate = 0;
        break;
      }
  
    } //end while(1)
  
    // update final times achieved overall sudivision
    time_steps->times[time_steps->tim]   = NR_t.times[NR_t.tim+1] - NR_t.dt[DT_NP1];
    time_steps->times[time_steps->tim+1] = NR_t.times[NR_t.tim+1];
  }
}

/// Multiscale simulation interface to perform Newton Raphson iteration
///
/// data structures have defined for multiscale simulation differently.
/// In order to use multiphysics data sturcture for multiscale simulation,
/// data should be reformed from data structure from multiscale to data structure for multiphysics
///
/// \param[in] print_level print level for a summary of the entire function call
/// \param[in] c structure of macroscale information
/// \param[in,out] s contains the information for the history-dependent solution
/// \param[in] solver_file structure for storing/updating the data
/// \param[in] ctx container for passing through Newton Raphson
/// \param[in] opts structure PGFem3D option
/// \param[in] sup_defl Prescribed deflection
/// \param[out] pores opening volume of failed cohesive interfaces
/// \param[out] n_step the number of nonlinear steps taken to solve the given increment
/// \return time spent in linear solver (seconds).
double Newton_Raphson_multiscale(const int print_level,
                                 MultiscaleCommon *c,
                                 MULTISCALE_SOLUTION *s,
                                 SOLVER_FILE *solver_file,
                                 MS_SERVER_CTX *ctx,
                                 const PGFem3D_opt *opts,
                                 double *sup_defl,
                                 double *pores,
                                 int *n_step, int tag_TP)
{
  int mp_id = 0;

  // initialize and define multiphysics
  Multiphysics mp;
  int id = MULTIPHYSICS_MECHANICAL;
  int ndim = c->ndofn;
  int write_no = MECHANICAL_Var_NO;

  mp.write_ids.resize(MECHANICAL_Var_NO);
  vector<int> coupled_ids;
  char *physicsname = (char *) malloc(sizeof(char)*1024);
  {
    for(int ia=0; ia<MECHANICAL_Var_NO; ia++)
      mp.write_ids[ia].push_back(ia); //sets ia as the fist element of each vector

    coupled_ids.push_back(0);
    sprintf(physicsname, "Mechanical");

    mp.physicsno      = 1;
    mp.physicsname    = &physicsname;
    mp.physics_ids    = &id;
    mp.ndim           = &ndim;
    mp.write_no       = &write_no;
    mp.coupled_ids.push_back(coupled_ids);
    mp.total_write_no = MECHANICAL_Var_NO;
  }
  ode_time.resize(mp.physicsno);

  // initialize and define mesh object
  Grid grid;
  grid_initialization(&grid);
  {
    grid.ne          = c->ne;
    grid.nn          = c->nn;
    grid.element     = c->elem;
    grid.b_elems     = NULL;
    grid.node        = c->node;
    grid.nce         = c->nce;
    grid.coel        = c->coel;
  }

  // initialize and define field variables
  FieldVariables fv;
  {
    field_varialbe_initialization(&fv);
    fv.ndofn  = c->ndofn;
    fv.ndofd  = c->ndofd;
    fv.npres  = c->npres;
    fv.sig    = s->sig_e;
    fv.eps    = s->eps;
    fv.u_np1  = s->r;
    fv.f      = s->f;
    fv.d_u    = s->d_r;
    fv.dd_u   = s->rr;
    fv.R      = s->R;
    fv.f_defl = s->f_defl;
    fv.RR     = s->RR;
    fv.f_u    = s->f_u;
    fv.RRn    = s->RRn;
    fv.pores  = *pores;
    fv.BS_x   = s->BS_x;
    fv.BS_f   = s->BS_f;
    fv.BS_RR  = s->BS_RR;
    fv.BS_f_u = s->BS_f_u;
    fv.NORM   = s->NORM;
  }

  /// initialize and define iterative solver object
  Solver sol{};
  {
    if(solver_file==NULL && tag_TP==0) // taylor micro
    {
      sol.nor_min  = c->lin_err;
      sol.FNR      = 5;
      sol.iter_max = c->maxit_nl;
      PGFEM_printf("\n NR-multiscale (micro msmID=%ld) ---> Taylor FNR=%ld\n", tag_TP, sol.FNR);
    }
    else if(solver_file==NULL && tag_TP==1) // pde micro
    {
      sol.nor_min  = c->lin_err;
      sol.FNR      = 1;
      sol.iter_max = c->maxit_nl;
      PGFEM_printf("\n NR-multiscale (micro msmID=%ld) ---> pde FNR=%ld\n", tag_TP, sol.FNR);
    }
    else // macro pde
    {
      sol.nor_min  = solver_file->nonlin_tol;
      sol.FNR      = solver_file->nonlin_method;
      sol.iter_max = solver_file->max_nonlin_iter;
      PGFEM_printf("\n NR-multiscale default macro ---> PDE FNR=%ld\n", sol.FNR);
    }
    sol.n_step     = *n_step;
    sol.system     = c->SOLVER;
    sol.err        = c->lin_err;
    sol.microscale = ctx;
  }

  // initialize and define loading steps object
  LoadingSteps load;
  {
    loading_steps_initialization(&load);
    load.sups     = &(c->supports);
    load.sup_defl = &sup_defl;
  }

  // initialize and define material properties
  MaterialProperty mat;
  {
    material_initialization(&mat);
    mat.hommat  = c->hommat;
    mat.matgeom = c->matgeom;
  }

  /// initialize and define communication structures
  CommunicationStructure com;
  {
    communication_structure_initialization(&com);
    com.Ap     = c->Ap;
    com.Ai     = c->Ai;
    com.DomDof = c->DomDof;
    com.GDof   = c->GDof;
    com.nbndel = c->nbndel;
    com.bndel  = c->bndel;
    com.boot   = c->boot;
    com.net    = c->net;
    com.comm   = c->comm;
    com.rank   = c->rank;
    com.nproc  = c->nproc;
    com.spc    = c->spc;
  }

  /// initialize and define time stepping variable
  TimeStepping ts;
  {
    time_stepping_initialization(&ts);
    if(solver_file==NULL)
      ts.nt = 1;
    else
      ts.nt  = solver_file->n_step;

    ts.tim = s->tim;
    ts.times = s->times;
    ts.dt_n   = 0.0;
    ts.dt_np1 = s->dt;
    ts.print  = NULL;
    ts.tns    = NULL;
  }

  NR_time_steps NR_t;
  set_time_step_info_for_NR(&ts,&NR_t);
  
  int is_NR_cvg = 1;
  double alpha = 0;
  double solve_time = 0.0; 
  double stiffmat_loc_time = 0.0; 
  double residuals_loc_time = 0.0; 
  //record the start of the full timestep
  NR_t.times[0] =  ts.times[ts.tim];
  perform_Newton_Raphson_with_subdivision(&solve_time,&stiffmat_loc_time,&residuals_loc_time,
                                          print_level,&is_NR_cvg,&alpha,&grid,&mat,&fv,&sol,
                                          &load,&com,&ts,s->crpl,c->VVolume,
                                          opts,mp,&NR_t,mp_id);

  update_values_for_next_NR(&grid,&mat,&fv,&sol,&load,s->crpl,&com,c->VVolume,
                            opts,mp,&NR_t,mp_id);
  //put everything back
  ts.times[ts.tim]     = NR_t.times[NR_t.tim + 1] - NR_t.dt[DT_NP1];
  ts.times[ts.tim + 1] = NR_t.times[NR_t.tim + 1];

  *n_step = sol.n_step;
  s->NORM = fv.NORM;
  *pores  = fv.pores;

  free(physicsname);

  return solve_time;
}
