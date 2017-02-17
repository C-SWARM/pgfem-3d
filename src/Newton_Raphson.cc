#include "Newton_Raphson.h"

#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>


#include "PGFem3D_data_structure.h"
#include "PGFEM_io.h"
#include "enumerations.h"
#include "fd_increment.h"
#include "fd_residuals.h"
#include "integration.h"
#include "vol_damage_int_alg.h"
#include "LINE.h"
#include "load.h"
#include "matice.h"
#include "press_theta.h"
#include "res_fini_def.h"
#include "stabilized.h"
#include "stiffmat_fd.h"
#include "subdivision.h"
#include "utils.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "displacement_based_element.h"
#include "matrix_printing.h"
#include "interface_macro.h"
#include "compute_reactions.h"
#include "bounding_element_utils.h"
#include "solve_system.h"
#include "vtk_output.h"
#include "ms_cohe_job_list.h"
#include "macro_micro_functions.h"
#include "three_field_element.h"

#include "pgf_fe2_macro_client.h"
#include "pgf_fe2_micro_server.h"

#include "constitutive_model.h"
#include "dynamics.h"
#include "energy_equation.h"
#include "PGFem3D_data_structure.h"
#include "get_dof_ids_on_elem.h"
#include "solver_file.h"

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

/* MINIMAL_OUTPUT prints a summary of the entire function call. For
 * any print_level > MINIMAL_OUTPUT, normal output is used. */
enum{MINIMAL_OUTPUT,NORMAL_OUTPUT,VERBOSE_OUTPUT} PRINT_LEVEL;

static const int periodic = 0;

typedef struct NR_summary{
  int total_iter;
  double final_rel_res;
  double final_res;
} NR_summary;

typedef struct {
  double times[3];
  double dt[2];

} NR_time_steps;

/** Set the time vector to contain current dt for microscale */
static void set_time_micro(const int tim,
                           double *times,
                           const double dt,
                           const long DIV,
                           double *time_np1)
{
  *time_np1 = times[tim+1];
  times[tim + 1] = times[tim] + dt;
}

/** Reset the time vector to normal state */
static void set_time_macro(const int tim,
                           double *times,
                           const double time_np1)
{
  times[tim + 1] = time_np1;
}

/// check velocity increment 
/// Compute residuals for Newton Raphson iteration
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dts time step sizes a n, and n+1
/// \return non-zero on internal error
int is_displacement_acceptable(double *alpha,
                               GRID *grid,
                               FIELD_VARIABLES *fv,
                               LOADING_STEPS *load,
                               const PGFem3D_opt *opts,
                               int mp_id)
{
  int err = 0;
  *alpha = 0.0;
  return err;

  SUPP sup = load->sups[mp_id];
  int is_total_lagrangian = 0;

  if(sup->multi_scale)
    is_total_lagrangian = 1;
  else
  {
    switch(opts->analysis_type)
    {
      case DISP:
      case TF:
        is_total_lagrangian = 1;
        break;
      case CM:
        if(opts->cm != UPDATED_LAGRANGIAN)
          is_total_lagrangian = 1;
        break;
    }
  }

  double alpha_V = 0.0;
  
  for(int e=0; e<grid->ne; e++)
  {
    int nne   = grid->element[e].toe;
    int ndofn = fv->ndofn;
    int ndofe = nne*ndofn;

    long *nod = (long *) malloc(sizeof(long)*nne);
    long *cn = aloc1l (ndofe);
    double *u_e = aloc1(ndofe);

    elemnodes(e,nne,nod,grid->element);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,grid->node,cn,mp_id);    
    def_elem_total(cn,ndofe,fv->u_np1,fv->d_u,grid->element,grid->node,load->sups[mp_id],u_e);

    double *X = aloc1(nne);
    double *Y = aloc1(nne);
    double *Z = aloc1(nne);

    double *x_n = aloc1(nne);
    double *y_n = aloc1(nne);
    double *z_n = aloc1(nne);

    double *x_np1 = aloc1(nne);
    double *y_np1 = aloc1(nne);
    double *z_np1 = aloc1(nne);
  
    double max_disp[3];
    max_disp[0]  = max_disp[1] = max_disp[2] == 0.0;
 
    for(int n=0; n<nne; n++)
    {
      long nid = nod[n];
 
      X[n]   = grid->node[nid].x1_fd;
      Y[n]   = grid->node[nid].x2_fd;
      Z[n]   = grid->node[nid].x3_fd;     

      x_n[n]   = grid->node[nid].x1;
      y_n[n]   = grid->node[nid].x2;
      z_n[n]   = grid->node[nid].x3;

      x_np1[n] = grid->node[nid].x1_fd + u_e[n*ndofn + 0];
      y_np1[n] = grid->node[nid].x2_fd + u_e[n*ndofn + 1];
      z_np1[n] = grid->node[nid].x3_fd + u_e[n*ndofn + 2];
      
      max_disp[0] = (fabs(max_disp[0]) > fabs(x_np1[n] - x_n[n]))? max_disp[0] : x_np1[n] - x_n[n];
      max_disp[1] = (fabs(max_disp[1]) > fabs(y_np1[n] - y_n[n]))? max_disp[1] : y_np1[n] - y_n[n];
      max_disp[2] = (fabs(max_disp[2]) > fabs(z_np1[n] - z_n[n]))? max_disp[2] : z_np1[n] - z_n[n];
    }
    double max_x[3], min_x[3], dx[3];
    max_x[0] = max_x[1] = max_x[2] = -1.0e+15;
    min_x[0] = min_x[1] = min_x[2] =  1.0e+15;
   
    for(int n = 0; n<nne; n++)
    {
       max_x[0] = (max_x[0] > x_n[n])? max_x[0]: x_n[n];
       max_x[1] = (max_x[1] > y_n[n])? max_x[1]: y_n[n];
       max_x[2] = (max_x[2] > z_n[n])? max_x[2]: z_n[n];

       min_x[0] = (min_x[0] < x_n[n])? min_x[0]: x_n[n];
       min_x[1] = (min_x[1] < y_n[n])? min_x[1]: y_n[n];
       min_x[2] = (min_x[2] < z_n[n])? min_x[2]: z_n[n];
    }

    dx[0] = max_x[0] - min_x[0];
    dx[1] = max_x[1] - min_x[1];
    dx[2] = max_x[2] - min_x[2];
    
    double ratio = 0.0;
    for(int ia=0; ia<3; ia++)
    {
      double ratio_ia = fabs(max_disp[ia])/dx[ia]/0.07;
      ratio = (ratio > ratio_ia) ? ratio : ratio_ia;
    }

    double V_0   = compute_volumes_from_coordinates(X,  Y,  Z,nne);
    double V_n   = compute_volumes_from_coordinates(x_n,  y_n,  z_n,nne);
    double V_np1 = compute_volumes_from_coordinates(x_np1,y_np1,z_np1,nne);
   
    int myrank = 0;
    MPI_Comm mpi_comm = MPI_COMM_WORLD;
    MPI_Comm_rank (mpi_comm,&myrank);
    if(0)
    {
      printf("%d %d: %e -> %e: %e\n", myrank,e,V_n, V_np1, fabs((V_np1-V_n)/V_n));
      printf("dx = (%e %e %e), disp = (%e %e %e)\n", dx[0],dx[1],dx[2], max_disp[0],max_disp[1],max_disp[2]);
    }
    
    double dV_max = 0.2;
    ratio = fabs(V_np1-V_n)/V_n/dV_max;
    alpha_V = (alpha_V > ratio)? alpha_V : ratio;   

    if(0)
    {
      printf("%d: %e -> %e: %e\n", e,V_n, V_np1, fabs((V_np1-V_n)/V_n));
      printf("dx = (%e %e %e), disp = (%e %e %e)\n", dx[0],dx[1],dx[2], max_disp[0],max_disp[1],max_disp[2]);
    }

    free(nod);
    free(cn);
    free(u_e);

    free(X);
    free(Y);
    free(Z); 
   
    free(x_n);
    free(y_n);
    free(z_n);

    free(x_np1);
    free(y_np1);
    free(z_np1);    

//    if(V_np1<0)
//    {  
//      err = 1; 
//      break;
//    }    
  }
  
  *alpha = alpha_V;
  return err;
}                               


/// reset variables for Newton Raphson iteration
/// If Newton Raphson is restarted, reset variables to inital
///
/// \param[in] grid a mesh object
/// \param[in,out] variables object for field variables
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int reset_variables_for_NR(GRID *grid,
                           FIELD_VARIABLES *fv,
                           CRPL *crpl,
                           const PGFem3D_opt *opts)
{
  int err = 0;
  switch(opts->analysis_type){
    case FS_CRPL: case FINITE_STRAIN:
      res_fini_def (grid->ne,fv->npres,grid->element,fv->eps,fv->sig,
              crpl,opts->analysis_type);
      break;
    case STABILIZED:
      res_stab_def (grid->ne,fv->npres,grid->element,fv->eps,fv->sig,opts->stab);
      break;
    case MINI:
      MINI_reset(grid->element,grid->ne,fv->npres,fv->sig);
      break;
    case MINI_3F:
      MINI_3f_reset(grid->element,grid->ne,fv->npres,4,fv->sig,fv->eps);
      break;
    case CM:
      constitutive_model_reset_state(fv->eps, grid->ne, grid->element);
      break;
    default: break;
  }  
  return err;
}

/// Compute residuals for Newton Raphson iteration
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dts time step sizes a n, and n+1
/// \return non-zero on internal error
long compute_residuals_for_NR(GRID *grid,
                              MATERIAL_PROPERTY *mat,
                              FIELD_VARIABLES *fv,
                              SOLVER_OPTIONS *sol,
                              LOADING_STEPS *load,
                              CRPL *crpl,
                              MPI_Comm mpi_comm,
                              const PGFem3D_opt *opts,
                              MULTIPHYSICS *mp,
                              int mp_id,
                              double t,
                              double *dts,
                              int updated_deformation)
{
  long INFO;
  switch(mp->physics_ids[mp_id])
  {
    case MULTIPHYSICS_MECHANICAL:
      INFO = fd_residuals_MP(grid,mat,fv,sol,load,crpl,mpi_comm,opts,mp,mp_id,t,dts,updated_deformation);
      break;
    case MULTIPHYSICS_THERMAL:
        INFO = energy_equation_compute_residuals(grid,mat,fv,load,mp_id,updated_deformation,dts[DT_NP1]);
      break;
    case MULTIPHYSICS_CHEMICAL: //intented flow, not yet been implemented
    default: 
      printf("%s is not defined (compute_residuals_for_NR)\n", mp->physicsname[mp_id]);
  }
     
  return INFO;
}

/// Compute stiffnes
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] variables object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com communication object
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] t time
/// \param[in] dt time step
/// \param[in] iter number of Newton Raphson interataions
/// \param[in] myrank current process rank
/// \return non-zero on internal error
long compute_stiffness_for_NR(GRID *grid,
                              MATERIAL_PROPERTY *mat,
                              FIELD_VARIABLES *fv,
                              SOLVER_OPTIONS *sol,
                              LOADING_STEPS *load,
                              COMMUNICATION_STRUCTURE *com,
                              CRPL *crpl,
                              MPI_Comm mpi_comm,
                              const PGFem3D_opt *opts,
                              MULTIPHYSICS *mp,
                              int mp_id,
                              double dt,
                              long iter,
                              int myrank)
{ 
  long INFO = 0;
  switch(mp->physics_ids[mp_id])
  {
    case MULTIPHYSICS_MECHANICAL:
      INFO = stiffmat_fd_MP(grid,mat,fv,sol,load,com,crpl,mpi_comm,opts,mp,mp_id,dt,iter,myrank);
      break;
    case MULTIPHYSICS_THERMAL:
      INFO = energy_equation_compute_stiffness(grid,mat,fv,sol,load,com,mpi_comm,myrank,opts,mp_id,dt);
      break;
    case MULTIPHYSICS_CHEMICAL: //intented flow, not yet been implemented
    default: 
      printf("%s is not defined (compute_stiffness_for_NR)\n", mp->physicsname[mp_id]);
  }

  return INFO;  
}
int update_values_for_next_NR(GRID *grid,
                               MATERIAL_PROPERTY *mat,
                               FIELD_VARIABLES *fv,
                               SOLVER_OPTIONS *sol,
                               LOADING_STEPS *load,
                               CRPL *crpl,
                               MPI_Comm mpi_comm,
                               const double VVolume,
                               const PGFem3D_opt *opts,
                               MULTIPHYSICS *mp,
                               int mp_id,
                               double t,
                               double *dts)
{
  int err = 0;
  double dt = dts[DT_NP1];
  
  switch(mp->physics_ids[mp_id])
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
                        sol->nor_min,crpl,dt,grid->nce,grid->coel,&(fv->pores),mpi_comm,
                        VVolume,opts, mp_id);
          break;
        case STABILIZED:
          st_increment (grid->ne,grid->nn,fv->ndofn,fv->ndofd,mat->matgeom,mat->hommat,
                        grid->element,grid->node,load->sups[mp_id],fv->eps,fv->sig,fv->d_u,fv->u_np1,
                        sol->nor_min,opts->stab,dt,grid->nce,grid->coel,&(fv->pores),mpi_comm,
                        opts->cohesive,mp_id);
          break;
        case MINI:
          MINI_increment(grid->element,grid->ne,grid->node,grid->nn,fv->ndofn,
                         load->sups[mp_id],fv->eps,fv->sig,mat->hommat,fv->d_u,mpi_comm,mp_id);
          break;
        case MINI_3F:
          MINI_3f_increment(grid->element,grid->ne,grid->node,grid->nn,fv->ndofn,
                            load->sups[mp_id],fv->eps,fv->sig,mat->hommat,fv->d_u,mpi_comm,mp_id);
          break;
        case DISP:
          DISP_increment(grid->element,grid->ne,grid->node,grid->nn,fv->ndofn,load->sups[mp_id],fv->eps,
                         fv->sig,mat->hommat,fv->d_u,fv->u_np1,mpi_comm,mp_id);
          break;
        case TF:
          update_3f_state_variables(grid->ne,fv->ndofn,fv->npres,fv->d_u,fv->u_np1,grid->node,grid->element,mat->hommat,load->sups[mp_id],fv->eps,fv->sig,
                                    dt,t,mpi_comm,mp_id);
          break;
        case CM:
        {
          switch(opts->cm)
          {
            case HYPER_ELASTICITY: case DISP:
              DISP_increment(grid->element,grid->ne,grid->node,grid->nn,fv->ndofn,load->sups[mp_id],fv->eps,
                             fv->sig,mat->hommat,fv->d_u,fv->u_np1,mpi_comm,mp_id);
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
        if(opts->analysis_type==CM && opts->cm==UPDATED_LAGRANGIAN)
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
  
  if(mp->physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL)
  {  
    /* update of internal fv and nodal coordinates */
    if(opts->analysis_type==CM)
    {
      switch(opts->cm){
        case UPDATED_LAGRANGIAN:
        case TOTAL_LAGRANGIAN:
          constitutive_model_update_time_steps_test(grid->element,grid->node,fv->eps,grid->ne,grid->nn,
                                                    fv->ndofn,fv->u_n,dt,opts->cm,mp_id);
          break;
        case MIXED_ANALYSIS_MODE:
          constitutive_model_update_time_steps_test(grid->element,grid->node,fv->eps,grid->ne,grid->nn,
                                                    fv->ndofn,fv->u_n,dt,1 /* TL */,mp_id);
          break;
        default: break;
      }
    }
    
    if(opts->analysis_type==TF)
    {
      int nVol = 1;
      for (int e=0;e<grid->ne;e++)
      {
        if(fv->npres==1)
        {
          fv->eps[e].d_T[2] = fv->eps[e].d_T[1];
          fv->eps[e].d_T[1] = fv->eps[e].d_T[0];
          
        }
        for(int a=0; a<nVol; a++)
        {
          fv->eps[e].T[a*3+2] = fv->eps[e].T[a*3+1];
          fv->eps[e].T[a*3+1] = fv->eps[e].T[a*3+0];
        }
      }
    }
  }
  // Update time steps
  dts[DT_N] = dts[DT_NP1];
  
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

/// Newton Raphson iterative solver
///
/// \param[in] print_level print level for a summary of the entire function call
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] VVolume original volume of the domain
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] times time info for this interation
/// \param[in, out] dts time step sizes from t(n-1) to t(n), and t(n) to t(n+1)
///                 updateed by subdivisions
/// \param[in] mp_id mutiphysics id
/// \param[in] myrank current process rank
/// \return time spent for this routine
double Newton_Raphson_test(const int print_level,
                           GRID *grid,
                           MATERIAL_PROPERTY *mat,
                           FIELD_VARIABLES *FV,
                           SOLVER_OPTIONS *SOL,
                           LOADING_STEPS *load,
                           COMMUNICATION_STRUCTURE *COM,
                           PGFem3D_TIME_STEPPING *time_steps,
                           CRPL *crpl,
                           MPI_Comm mpi_comm,
                           const double VVolume,
                           const PGFem3D_opt *opts,
                           MULTIPHYSICS *mp,
                           double *times,
                           double *dts,
                           int mp_id,
                           int myrank)                                        
{
  // use pointers for physics[mp_id]
  SOLVER_OPTIONS          *sol = SOL + mp_id;
  FIELD_VARIABLES         *fv =  FV  + mp_id;
  COMMUNICATION_STRUCTURE *com = COM + mp_id;
  
  double GNOR; // global norm of residual
  double nor1; //
  sol->gama = 0.0; 

  long tim = 1;  
  if(time_steps->tim == 0)
    tim = 0;
    
  if(myrank==0)
    PGFEM_printf("tim = %ld, time = (%e %e %e), dts = (%e %e)\n", tim, times[0],times[1], times[2],dts[0],dts[1]);  

  double dt = dts[DT_NP1];  
  double t = times[tim+1];
  
  long DIV, ST, GAMA, OME, i, j, N, M, INFO, iter, STEP, ART, GInfo, gam;
  double DT, NOR=10.0, ERROR, LS1, tmp, Gss_temp, nor2, nor;
  char str1[500];
  struct rusage usage;
  
  double enorm, Genorm;
  static double ENORM = 1.0;
  
  double solve_time = 0.0;
  double zero_tol = 1e-15;
  nor2 = 1.0;
  
  /* interface multiscale_modeling */
  double *macro_jump_u = aloc1(3);
  
  /* damage substep criteria */
  const double max_damage_per_step = 0.05;
  const double alpha_restart = 1.25;
  double max_damage = 0.0;
  double alpha = 0.0;
  
  /* max micro substep criteria */
  const int max_n_micro_substep = 2;
  const double alpha_restart_ms = 2.0;
  int max_substep = 0;
  double alpha_ms = 0.0;
  
  /* damage dissipation */
  double dissipation = 0.0;
  
  double BS_nor=0.0;
  int BS_iter;
  
  /* option '-no-migrate' */
  const int NR_REBALANCE = (opts->no_migrate)? FE2_REBALANCE_NONE : FE2_REBALANCE_ADAPTIVE;
  
  switch(opts->analysis_type){
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
  DIV = ST = GAMA = OME = INFO = ART = 0;
  STEP = 1;
  DT = 0.0;
  ERROR = sol->nor_min;
  iter = 0;
  
  /* GOTO REST */
  rest:
    fflush(PGFEM_stdout);
    
    if(INFO==1 && opts->solution_scheme_opt[LINE_SEARCH]==0)
    {
      ART = 1;
      if(myrank==0)
        printf("Imposed to use NO Line search [INFO = %ld, ART = %ld]\n", INFO, ART);        
    }
    
    if (INFO == 1 && ART == 0)
    {
      // reset variables
      if(mp->physics_ids[mp_id]==MULTIPHYSICS_MECHANICAL)
        reset_variables_for_NR(grid, fv, crpl, opts);
      
      for (i=0;i<fv->ndofd;i++) {
        fv->dd_u[i] = fv->d_u[i] = fv->f_defl[i] = fv->f[i] = 0.0;
      }
      
      if (myrank == 0 ) PGFEM_printf("\n** Try without LINE SEARCH **\n\n");
      
      ART = 1;
    } 
    else 
    {      
      if(mp->physics_ids[mp_id]==MULTIPHYSICS_MECHANICAL)
      {  
        subdivision(INFO,&dt,&STEP,&DIV,tim,times,&ST,grid->ne,
                    fv->ndofn,fv->ndofd,fv->npres,grid->element,crpl,fv->eps,fv->sig,
                    load->sups[mp_id],load->sup_defl[mp_id],fv->dd_u,fv->d_u,fv->f_defl,fv->f,fv->RRn,fv->R,&GAMA,
                    &DT,&OME,opts->stab,iter,sol->iter_max,alpha,mpi_comm,
                    opts->analysis_type);
      }
      else
      {
        subdivision(INFO,&dt,&STEP,&DIV,tim,times,&ST,grid->ne,
                    fv->ndofn,fv->ndofd,fv->npres,grid->element,crpl,fv->eps,fv->sig,
                    load->sups[mp_id],load->sup_defl[mp_id],fv->dd_u,fv->d_u,fv->f_defl,fv->f,fv->RRn,fv->R,&GAMA,
                    &DT,&OME,opts->stab,iter,sol->iter_max,alpha,mpi_comm,
                    -1);
      } 
                   
      if(STEP>1)
        sol->is_subdivided = 1;             
    
      dts[DT_NP1] = dt; // update dt_np1
                        // dt_n is updated when Newton Raphson is completed without error

      gam = ART = 0;
      
    }
    
    /* recompute the microscale tangent if restart due to error. */
    if(INFO == 1 && DEBUG_MULTISCALE_SERVER && sol->microscale != NULL){
      /* zero the macroscale tangent */
      ZeroHypreK(sol->PGFEM_hypre,com->Ai,com->DomDof[myrank]);
      
      /* start the microscale jobs. Do not compute equilibrium. Use
       d_r (no displacement increments) for displacement dof
       vector */
      MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
      pgf_FE2_macro_client_rebalance_servers(ctx->client,ctx->mpi_comm,
              FE2_REBALANCE_NONE);
      double tnp1 = 0;
      set_time_micro(tim,times,dt,DIV,&tnp1);
      pgf_FE2_macro_client_send_jobs(ctx->client,ctx->mpi_comm,ctx->macro,
                                     JOB_NO_COMPUTE_EQUILIBRIUM);
      set_time_macro(tim,times,tnp1);
    }
    
    /* reset the error flag */
    INFO = 0;
    while (STEP > DIV){
      if ((STEP > 1 || ST == 1) && myrank == 0 ){
        PGFEM_printf ("\nSTEP = %ld :: NS =  %ld || Time %e | dt = %e\n",
                DIV,STEP,times[tim+1],dt);
      }
      bounding_element_communicate_damage(grid->n_be,grid->b_elems,grid->ne,fv->eps,mpi_comm);
      if (periodic == 1){
        long nt = time_steps->nt;
        double factor = (times[tim] + (DIV+1)*dt)/times[nt] * fv->eps[0].load;
        /* Plane strain deviatoric tension */
        if (fv->eps[0].type == 1){
          fv->eps[0].F[0][0] = times[nt]/(times[nt]
                  - (times[tim] + (DIV+1)*dt)*fv->eps[0].load);
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
          fv->eps[0].F[1][1] = (1. - (times[tim] + (DIV+1)*dt)/times[nt] * fv->eps[0].load1);
          fv->eps[0].F[2][2] = (1. - (times[tim] + (DIV+1)*dt)/times[nt] * fv->eps[0].load1);
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
        INFO = compute_residuals_for_NR(grid,mat,fv,sol,load,crpl,mpi_comm,opts,mp,
                                        mp_id,t,dts, 1);
                
        for (i=0;i<fv->ndofd;i++){
          fv->f[i] = - fv->f_u[i];
          fv->R[i] = fv->RR[i] = 0.0;
        }
        
      }/* end periodic */
      else{
        
        for (i=0;i<(load->sups[mp_id])->npd;i++)
          (load->sups[mp_id])->defl_d[i] = load->sup_defl[mp_id][i]/STEP;
        
        /* Compute macro interface deformation gradient */
        if((load->sups[mp_id])->multi_scale){
          /* INFO = compute_interface_macro_jump_u(macro_jump_u,sup); */
          /* compute_interface_macro_grad_u(sup->F0,sup->lc,macro_jump_u, */
          /* 			       sup->N0); */
          
          INFO = compute_macro_grad_u((load->sups[mp_id])->F0,load->sups[mp_id],opts->analysis_type);
          
          if(INFO != 0){
            PGFEM_printerr("[%d] ERROR: not enough prescribed displacements"
                    " for interface multiscale modeling!\n"
                    "Must have at least six (6) prescribed displacements.\n"
                    "Check input and try again.\n",myrank);
            PGFEM_Comm_code_abort(mpi_comm,0);
          }
          nulld (fv->f_u,fv->ndofd);
          vol_damage_int_alg(grid->ne,fv->ndofn,fv->d_u,fv->u_np1,grid->element,grid->node,
                             mat->hommat,load->sups[mp_id],dt,iter,mpi_comm,
                             fv->eps,fv->sig,&max_damage,&dissipation,
                             opts->analysis_type,mp_id);
          
          compute_residuals_for_NR(grid,mat,fv,sol,load,crpl,mpi_comm,opts,mp,
                                   mp_id,t,dts, 0);
        } else {
          nulld (fv->f_u,fv->ndofd);
        }
        
        /*  NODE (PRESCRIBED DEFLECTION)
            - SUPPORT COORDINATES generation of the load vector  */
        nulld (fv->f_defl,fv->ndofd);

        compute_load_vector_for_prescribed_BC(grid,mat,fv,sol,load,dt,crpl,
                                              opts,mp,mp_id,myrank); 
        
        /* Generate the load and vectors */
        for (i=0;i<fv->ndofd;i++)  {
          fv->f[i] = fv->R[i]/STEP - fv->f_defl[i] - fv->f_u[i];
          //sum += f[i];
          fv->RR[i] = fv->RRn[i] + fv->R[i]/STEP*(DIV+1);
        }
      }
      
      /* Transform LOCAL load vector to GLOBAL */
      LToG (fv->f,fv->BS_f,myrank,com->nproc,fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);
      
      /* Transform LOCAL load vector to GLOBAL */
      LToG (fv->RR,fv->BS_RR,myrank,com->nproc,fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);
      
      iter = 0;
      nor = nor2 = GNOR = 10.0;

      while (nor > ERROR
              && nor2 > zero_tol
              /* && nor2 > ERROR*ERROR */ /* norm of the residual < error^2
               * MM 9/26/2012*/
              ){
        
        if (sol->FNR == 1 || (sol->FNR == 0 && iter == 0)){
          
          assert(opts->solverpackage == HYPRE);
          
          /* Null the matrix (if not doing multiscale)*/
          if(sol->microscale == NULL){
            ZeroHypreK(sol->PGFEM_hypre,com->Ai,com->DomDof[myrank]);
          }
          
          //compute stiffness matrix
          INFO =  compute_stiffness_for_NR(grid,mat,fv,sol,load,com,crpl,mpi_comm,
                                           opts,mp,mp_id,dt,iter,myrank);          
          
          // if INFO value greater than 0, the previous computation has an error
          MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_MAX,mpi_comm); 
          if (GInfo > 0) {
            if(myrank == 0){
              PGFEM_printf("Error detected (stiffmat_fd) %s:%s:%ld.\n"
                      "Subdividing load.\n", __func__, __FILE__, __LINE__);
            }
            INFO = 1;
            ART = 1;
            goto rest;
          }
          
          /* turn off line search for server-style multiscale */
          if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL){
            ART = 1;
            /* complete any jobs before assembly */
            MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
            pgf_FE2_macro_client_recv_jobs(ctx->client,ctx->macro,&max_substep);
          }
          
          /* Matrix assmbly */
          INFO = HYPRE_IJMatrixAssemble((sol->PGFEM_hypre)->hypre_k);
          
        }
        
        /*=== Solve the system of equations ===*/
        SOLVER_INFO s_info;

        solve_time += solve_system(opts,fv->BS_f,fv->BS_x,tim,iter,com->DomDof,&s_info,
                                   sol->PGFEM_hypre,mpi_comm);

        if(mp_id==-1)
        { 
          PGFEM_HYPRE_solve_info *PGFEM_hypre = sol->PGFEM_hypre;
          HYPRE_IJMatrixPrint(PGFEM_hypre->hypre_k,"test_k.txt");
          char fn_f[1024];
          sprintf(fn_f, "test_f.txt.%.5d", myrank);
          
          FILE *fp_f = fopen(fn_f, "w");
          fprintf(fp_f, "%d %d\n", PGFEM_hypre->ilower,
                                    PGFEM_hypre->iupper);
          for(int ia=0; ia<com->DomDof[myrank]; ia++)
            fprintf(fp_f, "%d %e\n", PGFEM_hypre->ilower +ia, fv->BS_f[ia]);
            
          fclose(fp_f);  
          
          exit(0);
        }
       
        if(myrank == 0){
          solve_system_check_error(PGFEM_stdout,s_info);
        }
        BS_nor = s_info.res_norm;
        BS_iter = s_info.n_iter;
        
        /* Check for correct solution */
        if (BS_nor > 500.*(sol->err) || BS_iter < 0) {
          INFO = 1;
          ART = 1;
          if (myrank == 0)
            PGFEM_printf("ERROR in the solver: nor = %8.8e || iter = %d\n",
                    BS_nor,BS_iter);
          goto rest;
        }
        /* Clear hypre errors */
        hypre__global_error = 0;
        
        if (!isfinite(BS_nor)) {
          if (myrank == 0)
            PGFEM_printf("ERROR in the solver: nor = %f\n",BS_nor);
          INFO = 1;
          ART = 1;
          goto rest;
        }
        
        /* Transform GLOBAL displacement vector to LOCAL */
        GToL (fv->BS_x,fv->dd_u,myrank,com->nproc,fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);
        
        /* LINE SEARCH */
        tmp  = ss (fv->BS_f,fv->BS_f,com->DomDof[myrank]);
        MPI_Allreduce(&tmp,&Gss_temp,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
        LS1 = 1./2.*Gss_temp;
        
        /* Pressure and volume change THETA */
        
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
            update_3f(grid->ne,fv->ndofn,fv->npres,fv->d_u,fv->u_np1,fv->dd_u,grid->node,grid->element,mat->hommat,load->sups[mp_id],fv->eps,fv->sig,
                      dt,t,mpi_comm,opts,sol->alpha,fv->u_n,fv->u_nm1,mp_id);
            
            break;
          default:
            break;
        }
        
        /*************************/
        /* INTEGRATION ALGORITHM */
        /*************************/
        if (opts->analysis_type == FS_CRPL) {
          INFO = integration_alg (grid->ne,fv->ndofn,fv->ndofd,fv->npres,crpl,grid->element,
                                  grid->node,fv->d_u,fv->dd_u,load->sups[mp_id],mat->matgeom,mat->hommat,
                                  fv->eps,fv->sig,tim,iter,dt,sol->nor_min,STEP,0,opts,mp_id);
          
          /* Gather INFO from all domains */
          // if INFO value greater than 0, the previous computation has an error
          MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_MAX,mpi_comm); 
          if (GInfo > 0) { // if not 0, an error is detected
            ART = 1;
            INFO = 1;
            goto rest;
          }
        }
        
        /* Update deformations */
        for (i=0;i<fv->ndofd;i++) {
          fv->f[i] = fv->d_u[i] + fv->dd_u[i];
          fv->f_u[i] = 0.0;
        }
        
        INFO = vol_damage_int_alg(grid->ne,fv->ndofn,fv->f,fv->u_np1,grid->element,grid->node,
                                  mat->hommat,load->sups[mp_id],dt,iter,mpi_comm,
                                  fv->eps,fv->sig,&max_damage,&dissipation,
                                  opts->analysis_type,mp_id);
        
        bounding_element_communicate_damage(grid->n_be,grid->b_elems,grid->ne,fv->eps,mpi_comm);
/*        
        if(mp->physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL)
        {        
          double element_volume_evolution = 0.0;
        
          INFO += is_displacement_acceptable(&element_volume_evolution,grid,fv,load,opts,mp_id);
          alpha = 10.0;
        }*/
        // if INFO value greater than 0, the previous computation has an error
        MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_MAX,mpi_comm);
        if (GInfo > 0) { // if not 0, an error is detected
          if(myrank == 0){
            PGFEM_printf("Inverted element detected (vol_damage_int_alg).\n"
                    "Subdividing load.\n");
          }
          INFO = 1;
          ART = 1;
          goto rest;
        }
        
        /* server-style multiscale */
        if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL){
          /* zero the macroscale tangent */
          ZeroHypreK(sol->PGFEM_hypre,com->Ai,com->DomDof[myrank]);
          
          /* start the microscale jobs */
          MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
          if ( iter == 0 ) {
            /* == Do not rebalance if this is the first iteration. ==
             * This is because most often it will be after an update or
             * print operation from the previous step. These operations
             * are approx. equal for all cells and thus the timing
             * information is not valid for rebalancing purposes. */
            pgf_FE2_macro_client_rebalance_servers(ctx->client,ctx->mpi_comm,
                                                   FE2_REBALANCE_NONE);
          } else {
            pgf_FE2_macro_client_rebalance_servers(ctx->client,ctx->mpi_comm,
                                                   NR_REBALANCE);
          }
          double tnp1 = 0;
          set_time_micro(tim,times,dt,DIV,&tnp1);
          pgf_FE2_macro_client_send_jobs(ctx->client,ctx->mpi_comm,ctx->macro,
                                         JOB_COMPUTE_EQUILIBRIUM);
          set_time_macro(tim,times,tnp1);
        }
        
        /* Residuals */
        INFO = compute_residuals_for_NR(grid,mat,fv,sol,load,crpl,mpi_comm,opts,mp,
                                        mp_id,t,dts, 1);
        
        // if INFO value greater than 0, the previous computation has an error
        MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_MAX,mpi_comm);
        if (GInfo > 0) { // if not 0, an error is detected
          if(myrank == 0){
            PGFEM_printf("Error detected (fd_residuals) %s:%s:%ld.\n"
                    "Subdividing load.\n", __func__, __FILE__, __LINE__);
          }
          INFO = 1;
          ART = 1;
          goto rest;
        }
        
        if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL){
          /* print_array_d(PGFEM_stdout,f_u,ndofd,1,ndofd); */
          MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
          pgf_FE2_macro_client_recv_jobs(ctx->client,ctx->macro,&max_substep);
          
          /* determine substep factor */
          alpha_ms = ((double) max_substep) / max_n_micro_substep;
          if(alpha_ms > alpha_restart_ms){
            if(myrank == 0){
              PGFEM_printf("Too many subdvisions at microscale (alpha_ms = %f).\n"
                      "Subdividing load.\n",alpha_ms);
            }
            alpha = alpha_ms;
            INFO = 1;
            ART = 1;
            goto rest;
          }
        }
        
        /* Transform LOCAL load vector to GLOBAL */
        LToG (fv->f_u,fv->BS_f_u,myrank,com->nproc,fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);
        
        /* Compute Euclidian norm */
        for (i=0;i<com->DomDof[myrank];i++)
          fv->BS_f[i] = fv->BS_RR[i] - fv->BS_f_u[i];
        
        nor  = ss (fv->BS_f,fv->BS_f,com->DomDof[myrank]);
        MPI_Allreduce(&nor,&GNOR,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
        nor2 = nor = sqrt(GNOR);
        
        if ((tim == 0 && iter == 0)
        || (iter == 0 && fv->NORM < ERROR)){ /* Reset *NORM if
         * less than convergence tolerance
         * MM 6/27/2012*/
          /* take maximum */
          if(nor > fv->NORM)	fv->NORM = nor;
          if(fv->NORM<sol->computer_zero)
            fv->NORM = sol->computer_zero;
        }
        
        
        /* THIS IS GLOBAL-LOCAL TOLERANCE */
        nor1 = fv->NORM;
        
        /* Normalize norm */
        nor /= nor1;
        
        if (!isfinite(nor)) {
          if (myrank == 0)
            PGFEM_printf("ERROR in the algorithm : nor = %f\n",nor);
          INFO = 1;
          ART = 1;
          goto rest;
        }
        
        /* My Line search */
        if (ART == 0) {
          INFO = LINE_S3_MP(grid,mat,FV,SOL,load,COM,time_steps,crpl,mpi_comm,opts,mp,
                            dts,t,mp_id,&nor,&nor2,nor1,NOR,LS1,iter,&max_damage,&dissipation,
                            tim,STEP);
          
          /* Gather infos */
          // if INFO value greater than 0, the previous computation has an error
          MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_MAX,mpi_comm);
          if (GInfo > 0) { /* ERROR in line search */
            if (myrank == 0){
              PGFEM_printf("Error in the line search algorithm\n");
            }
            if(NOR < 5. * ERROR){ /* MM 10/2/2012 soft convergence criteria */
              if(myrank == 0){
                PGFEM_printf("I will take last solution.\n");
              }
              nor = NOR;
              sol->gama = 0.0;
              INFO = 0;
              break; /* take previous as converged value */
              
            } else { /* MM 10/2/2012 not converged, subdivide */
              INFO = 1;
              goto rest;
            }
          } /* end ERROR */
          
          if (sol->gama != 1.0) gam = 1;
          
        } else {       /* no line search */
          sol->gama = 1.0;
        }
        
        /* Total deformation increment */
        NOR = nor;
        for (i=0;i<fv->ndofd;i++){
          fv->d_u[i] += (sol->gama)*fv->dd_u[i];
          fv->dd_u[i] *= sol->gama;
        }
        
        
        /* Compute the energy norm E_norm = abs(R ddu/(R0 ddu0)) */
        LToG(fv->dd_u,fv->BS_x,myrank,com->nproc,fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);
        enorm = ss(fv->BS_f,fv->BS_x,com->DomDof[myrank]);
        MPI_Allreduce(&enorm,&Genorm,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
        if((tim == 0 && iter == 0)) ENORM = Genorm;
        if(ENORM <= zero_tol) ENORM = 1.0;
        enorm = fabs(Genorm/ENORM);
        
        if (myrank == 0){
          getrusage (RUSAGE_SELF,&usage);
          PGFEM_printf("(%ld) IT = %d : R = %8.8e :: ||f||/||f0||",
                       iter,BS_iter,BS_nor);
          PGFEM_printf(" = [%8.8e] || [%8.8e] :: S %ld.%ld, U %ld.%ld\n",
                       nor, nor2, usage.ru_stime.tv_sec,usage.ru_stime.tv_usec,
                       usage.ru_utime.tv_sec, usage.ru_utime.tv_usec);
          PGFEM_printf("|R'u|/|R0'u0| = [%1.12e] || [%1.12e]\n",enorm,fabs(Genorm));
          fflush(PGFEM_stdout);
        }
        
        if(NR_UPDATE || PFEM_DEBUG || PFEM_DEBUG_ALL){
          if(opts->analysis_type == MINI){
            MINI_check_resid(fv->ndofn,grid->ne,grid->element,grid->node,mat->hommat,fv->eps,
                             fv->sig,fv->d_u,load->sups[mp_id],fv->RR,com->DomDof,fv->ndofd,
                             com->GDof,com->comm,mpi_comm,mp_id);
          }
          if(opts->analysis_type == MINI_3F){
            MINI_3f_check_resid(fv->ndofn,grid->ne,grid->element,grid->node,mat->hommat,fv->eps,
                                fv->sig,fv->d_u,load->sups[mp_id],fv->RR,com->DomDof,fv->ndofd,
                                com->GDof,com->comm,mpi_comm,mp_id);
          }
        }

        /* Check energy norm for convergence */
        if(opts->solution_scheme_opt[CVG_CHECK_ON_ENERGY_NORM])
        {  
          if((enorm < ERROR*ERROR) && (nor2 < 50*ERROR*ERROR) && iter > 0){
            if(myrank == 0){
              PGFEM_printf("Converged on energy norm\n");
            }
            INFO = 0;
            break;
          }
        }
        
        /* Max number of iterations restart */
        if (iter > (sol->iter_max-1) && nor > ERROR) {
          if (nor > 20.*ERROR || iter > (sol->iter_max + 2)) {
            if (nor < 5.*ERROR) {
              if (myrank == 0)
                PGFEM_printf ("I will take it\n");
              nor = 0.0;
            } else {
              if (myrank == 0)
                PGFEM_printf ("Error in the iteration : iter > iter_max (%d, %d)\n", sol->iter_max, iter);
              INFO = 1;
              if (gam == 0)
                ART = 1;
              goto rest;
            }
          }
        } /* end iter > iter_max */
        iter++; BS_nor = 0.0;
        
      }/* end while nor > nor_min */
      
      /* before increment after convergence, check max damage */
      if(mp->physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL)
      {
        // check from constitutive mode  
        if (opts->analysis_type == CM) {
          cm_get_subdivision_parameter(&alpha, grid->ne, grid->element, fv->eps, dt);
        } else {
          alpha = max_damage/max_damage_per_step;
        }
        // check displacement increment
        double element_volume_evolution = 0.0;
        is_displacement_acceptable(&element_volume_evolution,grid,fv,load,opts,mp_id);
        alpha = (alpha > element_volume_evolution)? alpha : element_volume_evolution;
      }
      MPI_Allreduce(MPI_IN_PLACE,&alpha,1,MPI_DOUBLE,MPI_MAX,mpi_comm);
      if(myrank == 0){
        PGFEM_printf("Physics based evolution thresh (e.g. damage): %f (wmax: %f)\n"
                     "Microscale subdivision alpha_ms: %f (max_substep: %d)\n",
                      alpha,max_damage_per_step,alpha_ms,max_n_micro_substep);
      }
      if(alpha > alpha_restart || alpha_ms > alpha_restart_ms){
        if(myrank == 0){
          PGFEM_printf("Subdividing to maintain accuracy of the material response.\n");
        }
        
        /* adapt time step based on largest adaption parameter */
        alpha = (alpha > alpha_ms)? alpha : alpha_ms;
        
        INFO = 1;
        ART = 1;
        goto rest;
      }
      
      /* always adapt time step based on largest adaption parameter */
      alpha = (alpha > alpha_ms)? alpha : alpha_ms;
      
      /* increment converged, output the volume weighted dissipation */
      if(myrank == 0){
        MPI_Reduce(MPI_IN_PLACE,&dissipation,1,
                   MPI_DOUBLE,MPI_SUM,0,mpi_comm);
        PGFEM_printf("Dissipation (vol weighted) [current] ||"
                     " [integrated] || (dt): %2.8e %2.8e %2.8e\n",
                     dissipation/dt, dissipation,dt);
      } else {
        double tmp = 0;
        MPI_Reduce(&dissipation,&tmp,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
      }
      
      ST = GAMA = gam = ART = 0;
      
      /* increment the step counter */
      (sol->n_step)++;
      
      /* /\* turn off line search *\/ */
      /* ART = 1; */
      
      /* microscale update. Overlay microscale update (using f = r+d_r)
       with macroscale update */
      if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL){
        /* start the microscale jobs */
        MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
        pgf_FE2_macro_client_rebalance_servers(ctx->client,ctx->mpi_comm,
                                               FE2_REBALANCE_NONE);
        
        double tnp1 = 0;
        set_time_micro(tim,times,dt,DIV,&tnp1);
        pgf_FE2_macro_client_send_jobs(ctx->client,ctx->mpi_comm,ctx->macro,
                                       JOB_UPDATE);
        set_time_macro(tim,times,tnp1);
      }
      
      // update converged values and apply increments 
      // for next Newton Raphson step while subdividing

      sol->last_residual = nor2;
      if(STEP>DIV+1)
        update_values_for_next_NR(grid,mat,fv,sol,load,crpl,mpi_comm,VVolume,opts,mp,mp_id,t,dts);

      /* finish microscale update */
      if(DEBUG_MULTISCALE_SERVER && sol->microscale != NULL){
        /* start the microscale jobs */
        MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) sol->microscale;
        pgf_FE2_macro_client_recv_jobs(ctx->client,ctx->macro,&max_substep);
      }
      
      /************* TEST THE UPDATE FROM N TO N+1  *************/
      if(NR_UPDATE || PFEM_DEBUG || PFEM_DEBUG_ALL){
        for (i=0;i<fv->ndofd;i++) {fv->f_u[i] = 0.0; fv->d_u[i] = 0.0;}
        
        if(opts->analysis_type == MINI){
          MINI_check_resid(fv->ndofn,grid->ne,grid->element,grid->node,mat->hommat,fv->eps,
                           fv->sig,fv->d_u,load->sups[mp_id],fv->RR,com->DomDof,fv->ndofd,
                           com->GDof,com->comm,mpi_comm,mp_id);
        }
        if(opts->analysis_type == MINI_3F){
          MINI_3f_check_resid(fv->ndofn,grid->ne,grid->element,grid->node,mat->hommat,fv->eps,
                              fv->sig,fv->d_u,load->sups[mp_id],fv->RR,com->DomDof,fv->ndofd,
                              com->GDof,com->comm,mpi_comm,mp_id);
        }
        compute_residuals_for_NR(grid,mat,fv,sol,load,crpl,mpi_comm,opts,mp,
                                 mp_id,t,dts, 0);
        
        for (i=0;i<fv->ndofd;i++) fv->f[i] = fv->RR[i] - fv->f_u[i];
        /* print_array_d(stdout,RR,ndofd,1,ndofd); */
        
        LToG(fv->f,fv->BS_f,myrank,com->nproc,fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);
        nor = ss(fv->BS_f,fv->BS_f,com->DomDof[myrank]);
        MPI_Allreduce(&nor,&tmp,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
        nor = sqrt (tmp);
        
        if (myrank == 0) PGFEM_printf("NORM NORM = %12.12e\n",nor);
        
      }
      /************* TEST THE UPDATE FROM N TO N+1  *************/
      
      DIV++;
      if(NR_PRINT_INTERMEDIATE){
        if(STEP > DIV){
          /* converged intermediate step, but not last step print output */
          char fname[100];
          sprintf(fname,"%s_%ld",opts->ofname,tim);
          if(myrank == 0){
            VTK_print_master(opts->opath,fname,sol->n_step,com->nproc,opts);
          }
          VTK_print_vtu(opts->opath,fname,sol->n_step,myrank,grid->ne,grid->nn,grid->node,
                        grid->element,load->sups[mp_id],fv->u_np1,fv->sig,fv->eps,opts,mp_id);
        }
      }
      
      double *res_trac = aloc1(3);
      bounding_element_compute_resulting_traction(grid->n_be,grid->b_elems,grid->element,grid->node,
                                                  fv->eps,fv->sig,fv->ndofd,com->DomDof,
                                                  com->GDof,com->comm,mpi_comm,
                                                  opts->analysis_type,
                                                  res_trac);
      free(res_trac);
      
      if (STEP > 2 && DIV == 2 /*&& iter < iter_max*/){
        goto rest;
      }
      
    }/*end SUBDIVISION */
    INFO = 0;
    OME = 0;
    free(macro_jump_u);
    
    if(NR_COMPUTE_REACTIONS && !(load->sups[mp_id])->multi_scale){
      compute_reactions(grid->ne,fv->ndofn,fv->npres,fv->u_np1,grid->node,grid->element,mat->matgeom,
                        mat->hommat,load->sups[mp_id],fv->eps,fv->sig,sol->nor_min,crpl,
                        dt,opts->stab,mpi_comm,opts->analysis_type,mp_id);
    }
    
    return (solve_time);    
}

/// compute residuals 
///
/// After converge, check residuals
/// This will not perform any integration algorithm.
///
/// \param[in] nor norm of residuals
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] t time at t(n+1)
/// \param[in] dts time step sizes a n, and n+1
/// \param[in] mp_id mutiphysics id
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int compute_coupled_physics_residual_norm(double *nor,
                                          GRID *grid,
                                          MATERIAL_PROPERTY *mat,
                                          FIELD_VARIABLES *FV,
                                          SOLVER_OPTIONS *SOL,
                                          LOADING_STEPS *load,
                                          COMMUNICATION_STRUCTURE *COM,
                                          PGFem3D_TIME_STEPPING *time_steps,
                                          CRPL *crpl,
                                          MPI_Comm mpi_comm,
                                          const PGFem3D_opt *opts,
                                          MULTIPHYSICS *mp,
                                          double t,
                                          double *dts,
                                          int mp_id,
                                          int myrank)
{
  int err = 0;
  // use pointers for physics[mp_id]
  SOLVER_OPTIONS          *sol = SOL + mp_id;
  FIELD_VARIABLES         *fv =  FV  + mp_id;
  COMMUNICATION_STRUCTURE *com = COM + mp_id;
  
  // temporal  
  double *u_n   = fv->u_n;
  double *u_nm1 = fv->u_nm1;
  State_variables *statv_list = fv->statv_list;
    
  double dt = dts[DT_NP1];  
  
  for(int ia=0; ia<fv->ndofd; ia++)
    fv->f_u[ia] = 0.0;  
  
  sol->run_integration_algorithm = 0; // turn off running integration algorithm
  long INFO = compute_residuals_for_NR(grid,mat,fv,sol,load,crpl,mpi_comm,opts,mp,
                                       mp_id,t,dts, 1);
  sol->run_integration_algorithm = 1; // reset integration algorithm to be active

  // Transform LOCAL load vector to GLOBAL
  LToG(fv->f_u,fv->BS_f_u,myrank,com->nproc,fv->ndofd,com->DomDof,com->GDof,com->comm,mpi_comm);

  // Compute Euclidian norm
  for(int i=0;i<com->DomDof[myrank]; i++)
    fv->BS_f[i] = fv->BS_RR[i] - fv->BS_f_u[i];
  
  double LNOR = 0.0;
  double GNOR = 0.0;
        
  LNOR = ss(fv->BS_f,fv->BS_f,com->DomDof[myrank]);
  MPI_Allreduce(&LNOR,&GNOR,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
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
/// \return time spent for this routine
int save_field_variables_to_temporal(GRID *grid,
                                     FIELD_VARIABLES *FV,
                                     const PGFem3D_opt *opts,
                                     MULTIPHYSICS *mp,
                                     int mp_id)
{
  FIELD_VARIABLES *fv = FV + mp_id;
  
  int err = 0;  
  for(int ia=0; ia<(grid->nn)*(fv->ndofn); ia++)
  {
    fv->temporal->u_nm1[ia] = fv->u_nm1[ia];
    fv->temporal->u_n[ia]   = fv->u_n[ia];
  }
  
  if(mp->physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL && opts->analysis_type==CM)
    err += constitutive_model_save_state_vars_to_temporal(FV + mp_id, grid);
    
  return err;
}

/// update variables(t(n+1)) to temporal_variables(t(n+1)) 
///
/// \param[in] grid a mesh object
/// \param[in,out] FV array of field variable object
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \return time spent for this routine
int update_temporal_field_variables_np1(GRID *grid,
                                        FIELD_VARIABLES *FV,
                                        const PGFem3D_opt *opts,
                                        MULTIPHYSICS *mp,
                                        int mp_id)
{
  FIELD_VARIABLES *fv = FV + mp_id;
  
  int err = 0;
  if(mp->physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL && opts->analysis_type==CM)
    err += constitutive_model_update_np1_state_vars_to_temporal(FV + mp_id, grid);
    
  return err;
}

/// reset variables(t(n) and t(n-1)) from temporal_variables(t(n) and t(n-1))  
///
/// \param[in] grid a mesh object
/// \param[in,out] FV array of field variable object
/// \param[in] mp mutiphysics object
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \return time spent for this routine
int reset_field_variables_using_temporal(GRID *grid,
                                         FIELD_VARIABLES *FV,
                                         const PGFem3D_opt *opts,
                                         MULTIPHYSICS *mp,
                                         int mp_id)
{                                         
  FIELD_VARIABLES *fv = FV + mp_id;
  
  int err = 0;  
  for(int ia=0; ia<(grid->nn)*(fv->ndofn); ia++)
  {
    fv->u_nm1[ia] = fv->temporal->u_nm1[ia];
    fv->u_n[ia]   = fv->temporal->u_n[ia];    
  }
  
  if(mp->physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL && opts->analysis_type==CM)
    err += constitutive_model_reset_state_using_temporal(FV + mp_id, grid);
    
  return err;
}

/// reset mesh coordinate (t(n)) from using saved displacement (t(n))
///
/// \param[in] grid a mesh object
/// \param[in,out] FV array of field variable object
/// \param[in] mp mutiphysics object
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \return time spent for this routine
int reset_mesh_coordinates(GRID *grid,
                           FIELD_VARIABLES *FV,
                           const PGFem3D_opt *opts,
                           MULTIPHYSICS *mp,
                           int mp_id)
{                                         
  FIELD_VARIABLES *fv = FV + mp_id;
  
  
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
/// \param[in] mp mutiphysics object
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] mp_id mutiphysics id
/// \return time spent for this routine
int reset_state_field_variables(GRID *grid,
                                FIELD_VARIABLES *FV,
                                SOLVER_OPTIONS *SOL,
                                const PGFem3D_opt *opts,
                                MULTIPHYSICS *mp,
                                int mp_id)
{
  int err = 0;

  if(mp->physics_ids[mp_id]==MULTIPHYSICS_MECHANICAL)
    err += reset_mesh_coordinates(grid,FV,opts,mp,mp_id);
    
  if(SOL[mp_id].is_subdivided)
    err += reset_field_variables_using_temporal(grid,FV,opts,mp,mp_id);
  
  return err;
}                                

/// set times for Netwon Raphson iterations
///
/// \param[in] ts object for time stepping
/// \param[out] t store times, t[0] = t(n-1)
///                            t[1] = t(n)
///                            t[2] = t(n+1)
/// \param[out] dts store dt , dts[DT_N]   = t(n)   - t(n-1)
///                          , dts[DT_NP1] = t(n+1) - t(n)
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int set_time_step_info_for_NR(PGFem3D_TIME_STEPPING *ts,
                              double *t,
                              double *dts,
                              int mp_id)
{
  int err = 0;

  long t_step_id = ts->tim;
    
//  t[0] = ts->tns[mp_id];
  if(ts->tim==0)
  {
    t[0] = ts->times[t_step_id];
    t[1] = ts->times[t_step_id+1];
    t[2] = 2.0*t[1] + t[0]; // this is dummy
  }
  else
  {
    t[0] = ts->times[t_step_id-1];
    t[1] = ts->times[t_step_id];
    t[2] = ts->times[t_step_id+1];                
  }
  dts[DT_N]   = t[1] - t[0];  
  dts[DT_NP1] = t[2] - t[1];  
      
  return err;
}

/// After complete Newton Raphson for physics X, check residuals of physics Y, and Z 
/// that are dependent on physics X.
///
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
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] NR_time save time step info during previous NR iteration
/// \param[in] mp_id mutiphysics id
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int check_convergence_of_NR_staggering(int *is_cnvged,
                                       GRID *grid,
                                       MATERIAL_PROPERTY *mat,
                                       FIELD_VARIABLES *FV,
                                       SOLVER_OPTIONS *SOL,
                                       LOADING_STEPS *load,
                                       COMMUNICATION_STRUCTURE *COM,
                                       PGFem3D_TIME_STEPPING *time_steps,
                                       CRPL *crpl,
                                       MPI_Comm mpi_comm,
                                       const PGFem3D_opt *opts,
                                       MULTIPHYSICS *mp,
                                       NR_time_steps *NR_time,
                                       int mp_id,
                                       int myrank)
{
  int err = 0;

  int tim = 1;
  if( time_steps->tim==0)
    tim = 0;    
  
  for(int ib = 0; ib<FV[mp_id].n_coupled; ib++)
  {
    int cpled_mp_id = mp->coupled_ids[mp_id][ib+1];
    int is_it_dependent_on_physics_mp_id = 0;
    for(int ic=0; ic<FV[cpled_mp_id].n_coupled; ic++)
    {
      if(mp->coupled_ids[cpled_mp_id][ic+1]==mp_id)
      { 
        // phsics[cpled_mp_id] is dependent on phsics[mp_id] 
        is_it_dependent_on_physics_mp_id = 1;
        break; 
      }
    }
    
    if(is_it_dependent_on_physics_mp_id==0)
      continue;
      
    double nor = 0.0;           

    compute_coupled_physics_residual_norm(&nor, grid,mat,FV,SOL,load,COM,time_steps,
                                          crpl,mpi_comm,opts,mp,
                                          NR_time[cpled_mp_id].times[2],
                                          NR_time[cpled_mp_id].dt,
                                          cpled_mp_id,myrank);

    double Rn_R = fabs(SOL[cpled_mp_id].last_residual - nor)/FV[cpled_mp_id].NORM;

    if(myrank==0)
    {
      printf("Check residuals for %s (tol = %e)\n", mp->physicsname[cpled_mp_id],SOL[cpled_mp_id].nor_min);
      printf("||Residual|| = %e, ||Residual(t=0)|| = %e, ||Residual||/||Residual(t=0)|| = %e, Last residual = %e, Rn-R = %e\n", 
             nor, FV[cpled_mp_id].NORM, 
             nor/FV[cpled_mp_id].NORM, 
             SOL[cpled_mp_id].last_residual,
             Rn_R);
    }
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
/// \param[in] nor norm of residuals
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object. FV.NORM is updateed.
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int set_0th_residual(GRID *grid,
                     MATERIAL_PROPERTY *mat,
                     FIELD_VARIABLES *FV,
                     SOLVER_OPTIONS *SOL,
                     LOADING_STEPS *load,
                     COMMUNICATION_STRUCTURE *COM,
                     PGFem3D_TIME_STEPPING *time_steps,
                     CRPL *crpl,
                     MPI_Comm mpi_comm,
                     const PGFem3D_opt *opts,
                     MULTIPHYSICS *mp,
                     int myrank)
{
  int err = 0;
  
  if(time_steps->tim>0)
    return err;
    
  for(int mp_id=0; mp_id<mp->physicsno; mp_id++)
  {
    if(SOL[mp_id].set_initial_residual==0)
      continue;
          
    for(int ia=0; ia<FV[mp_id].ndofd; ia++)
      FV[mp_id].u_np1[ia] = FV[mp_id].u0 + SOL[mp_id].du;

    double nor = 0.0;
    double times[3];
    double dts[2];
    set_time_step_info_for_NR(time_steps,times,dts,mp_id);
    compute_coupled_physics_residual_norm(&nor, grid,mat,FV,SOL,load,COM,time_steps,
                                          crpl,mpi_comm,opts,mp,times[2],dts,
                                          mp_id,myrank);

    FV[mp_id].NORM = nor; // set first residual
    
    for(int ia=0; ia<FV[mp_id].ndofd; ia++)
      FV[mp_id].u_np1[ia] = FV[mp_id].u0;

    if(myrank==0)
    {  
      printf("INFO. The first residual for the physics, %s, is computed as %e\n", mp->physicsname[mp_id], nor);
      printf("by pertubing the disp. with %e\n", SOL[mp_id].du);
    }
  }
  return err;    
}
      
/// Staggered Newton Raphson iterative solver for multiphysics problems
///
/// \param[in] print_level print level for a summary of the entire function call
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[in] SOL object array for solution scheme
/// \param[in] load object for loading
/// \param[in] COM object array for communications
/// \param[in] time_steps object for time stepping
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] VVolume original volume of the domain
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \return time spent for this routine
double Multiphysics_Newton_Raphson(const int print_level,
                                   GRID *grid,
                                   MATERIAL_PROPERTY *mat,
                                   FIELD_VARIABLES *FV,
                                   SOLVER_OPTIONS *SOL,
                                   LOADING_STEPS *load,
                                   COMMUNICATION_STRUCTURE *COM,
                                   PGFem3D_TIME_STEPPING *time_steps,
                                   CRPL *crpl,
                                   MPI_Comm mpi_comm,
                                   const double VVolume,
                                   const PGFem3D_opt *opts,
                                   MULTIPHYSICS *mp)
{
  double solve_time = 0.0;
  int max_itr = 10; // temporal setup for max. iteration
  //MPI rank
  int myrank;
  MPI_Comm_rank(mpi_comm,&myrank);

  // if too siff to converge, try pre compute resiudal by perturbing displacement 
  // slightly
  set_0th_residual(grid,mat,FV,SOL,load,COM,time_steps,crpl,mpi_comm,opts,mp,myrank);
  
  NR_time_steps *NR_time = (NR_time_steps *) malloc(sizeof(NR_time_steps)*mp->physicsno);

  //print start
  if(myrank==0)
  {  
    printf("======================================================================\n");
    printf(":: Do Newton Raphson for independent physics first\n");
    printf("======================================================================\n");    
  }
  
  // in Newton_Raphson iteration tim = 1 is the current time step id 
  // where times[0] = t(n-1)
  //       times[1] = t(n)
  //       times[2] = t(n+1)
  // s.t dt(n)   = times[1] - times[0]
  //     dt(n+1) = times[2] - times[1]
  // HOWEVER when actual time step is at 0 (time_steps->tim == 0)
  // tim = 0, because there is no t(n-1) values  
  int tim = 1;
  if( time_steps->tim==0)
    tim = 0;  

  int coupled_physics_no = 0;
  for(int mp_id=0; mp_id<mp->physicsno; mp_id++)
  {
    save_field_variables_to_temporal(grid,FV,opts,mp,mp_id);

    // obtain time steps for the current physics
    set_time_step_info_for_NR(time_steps,NR_time[mp_id].times,NR_time[mp_id].dt,mp_id);

    if(FV[mp_id].n_coupled>0)
    {
      coupled_physics_no++; 
      continue;
    } 
    //print current physics name
    if(myrank==0)
    {
      printf("----------------------------------------------------------------------\n");    
      printf("Newton Raphson iteration for : %s\n", mp->physicsname[mp_id]);  
      printf("----------------------------------------------------------------------\n");    
    }
    
    solve_time += Newton_Raphson_test(print_level,
                                      grid,mat,FV,SOL,load,COM,time_steps,
                                      crpl,mpi_comm,VVolume,opts,mp,NR_time[mp_id].times,NR_time[mp_id].dt,mp_id,myrank);
    
    reset_state_field_variables(grid,FV,SOL,opts,mp,mp_id);  
    update_temporal_field_variables_np1(grid,FV,opts,mp,mp_id);                                      
  }

  // if there is no physics coupled with others, Newton_Raphson is done here.
  if(coupled_physics_no>0)
  {
    //print Staggered Newton Raphson iteration
    if(myrank==0)
    {
      printf("======================================================================\n");    
      printf(":: Staggered Newton Raphson iteration starts for dependent physics \n");
      printf("======================================================================\n");    
    }

    for(int ia=0; ia<max_itr; ia++)
    {
      int is_cnvged = 1;
      
      for(int mp_id=0; mp_id<mp->physicsno; mp_id++)
      {
        if(FV[mp_id].n_coupled >0)
        {        
          //print current physics name
          if(myrank==0)
          {
            printf("----------------------------------------------------------------------\n");    
            printf("(%d) Newton Raphson iteration for : %s\n", ia, mp->physicsname[mp_id]);  
            printf("----------------------------------------------------------------------\n");    
          }
                    
          set_time_step_info_for_NR(time_steps,NR_time[mp_id].times,NR_time[mp_id].dt,mp_id);

          SOL[mp_id].is_subdivided = 0;
          solve_time += Newton_Raphson_test(print_level,
                                            grid,mat,FV,SOL,load,COM,time_steps,
                                            crpl,mpi_comm,VVolume,opts,mp,NR_time[mp_id].times,NR_time[mp_id].dt,mp_id,myrank);

          update_temporal_field_variables_np1(grid,FV,opts,mp,mp_id);         
          check_convergence_of_NR_staggering(&is_cnvged,grid,mat,FV,SOL,load,COM,time_steps,
                                             crpl,mpi_comm,opts,mp,NR_time,mp_id,myrank);                                                    
        }        
      }
      
      for(int mp_id=0; mp_id<mp->physicsno; mp_id++)
      {
        if(FV[mp_id].n_coupled == 0)
          continue;
        reset_state_field_variables(grid,FV,SOL,opts,mp,mp_id);      
      }
      
      if(is_cnvged)
      {
        if(myrank==0)
        {
          printf("======================================================================\n");    
          printf(":: Staggered Newton Raphson is converged. Done.\n");
          printf("======================================================================\n");     
        }
        break;
      } 
    }
  }
  
  // update final results
  for(int mp_id=0; mp_id<mp->physicsno; mp_id++)
  {
    time_steps->tns[mp_id] = NR_time[mp_id].times[tim+1] - NR_time[mp_id].dt[DT_NP1];
    update_values_for_next_NR(grid,mat,FV+mp_id,SOL+mp_id,load,crpl,mpi_comm,VVolume,opts,mp,mp_id,NR_time[mp_id].times[tim+1],NR_time[mp_id].dt); 
  }
  
  free(NR_time);

  return solve_time;  
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
/// \param[out] the number of nonlinear steps taken to solve the given increment
/// \return time spent in linear solver (seconds).
double Newton_Raphson_multiscale(const int print_level,
                                 COMMON_MACROSCALE *c,
                                 MACROSCALE_SOLUTION *s,
                                 SOLVER_FILE *solver_file,
                                 MS_SERVER_CTX *ctx,
                                 const PGFem3D_opt *opts,
                                 double *sup_defl,
                                 double *pores,
                                 int *n_step)
{
  // MPI stuff
  int nproc,myrank;
  MPI_Comm_size(c->mpi_comm,&nproc);
  MPI_Comm_rank(c->mpi_comm,&myrank);
  
  int mp_id = 0;
  
  // initialize and define multiphysics
  MULTIPHYSICS mp;  
  int id = MULTIPHYSICS_MECHANICAL;
  int ndim = c->ndofn;
  int write_no = MECHANICAL_Var_NO;

  int *write_ids = (int *) malloc(sizeof(int)*MECHANICAL_Var_NO);
  int *coupled_ids = (int *) malloc(sizeof(int));
  char *physicsname = (char *) malloc(sizeof(char)*1024);
  {
    for(int ia=0; ia<MECHANICAL_Var_NO; ia++)
      write_ids[ia] = ia;
    
    coupled_ids[0] = 0;
    sprintf(physicsname, "Mechanical");
 
    mp.physicsno      = 1;
    mp.physicsname    = &physicsname;
    mp.physics_ids    = &id;
    mp.ndim           = &ndim;
    mp.write_no       = &write_no;
    mp.write_ids      = &write_ids;
    mp.coupled_ids    = &coupled_ids;
    mp.total_write_no = MECHANICAL_Var_NO;
  }
  
  // initialize and define mesh object  
  GRID grid;
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
  FIELD_VARIABLES fv;
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
  SOLVER_OPTIONS sol;
  {
    solution_scheme_initialization(&sol);
    if(solver_file==NULL)
    {
      sol.nor_min  = c->lin_err;
      sol.FNR      = 1; // full NR
      sol.iter_max = c->maxit_nl;      
    }
    else
    {
      sol.nor_min  = solver_file->nonlin_tol;
      sol.FNR      = solver_file->nonlin_method;
      sol.iter_max = solver_file->max_nonlin_iter;
    }
    sol.PGFEM_hypre  = c->SOLVER;
    sol.err          = c->lin_err;
    sol.microscale   = ctx;                     
  }
  
  // initialize and define loading steps object
  LOADING_STEPS load;    
  {
    loading_steps_initialization(&load);
    load.sups     = &(c->supports);
    load.sup_defl = &sup_defl;
  }
  
  // initialize and define material properties
  MATERIAL_PROPERTY mat;
  {
    material_initialization(&mat);
    mat.hommat  = c->hommat;
    mat.matgeom = c->matgeom;
  }
  
  /// initialize and define communication structures
  COMMUNICATION_STRUCTURE com;
  {
    communication_structure_initialization(&com); 
    com.nproc  = nproc;
    com.Ap     = c->Ap;
    com.Ai     = c->Ai;
    com.DomDof = c->DomDof;
    com.comm   = c->pgfem_comm;
    com.GDof   = c->GDof;
    com.nbndel = c->nbndel;
    com.bndel  = c->bndel;
  }
  
  /// initialize and define time stepping variable
  PGFem3D_TIME_STEPPING ts;
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

  int tim = 1;
  if(ts.tim==0)
    tim = 0; 
      
  double times[3];
  double dts[2];
  
  set_time_step_info_for_NR(&ts,times,dts,0);  
  	  
  double solve_time = Newton_Raphson_test(print_level,&grid,&mat,&fv,&sol,&load,&com,&ts,
                                          s->crpl,c->mpi_comm,c->VVolume,
                                          opts,&mp,times,dts,mp_id,myrank);

  update_values_for_next_NR(&grid,&mat,&fv,&sol,&load,s->crpl,c->mpi_comm,c->VVolume,
                            opts,&mp,mp_id,times[tim+1],dts); 
                                          
  *n_step = sol.n_step;
  s->NORM = fv.NORM;
  *pores  = fv.pores; 
  
  free(write_ids);
  free(coupled_ids);
  free(physicsname);
  
  return solve_time;
}                                        
