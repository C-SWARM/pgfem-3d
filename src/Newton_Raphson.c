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

double Newton_Raphson_test(const int print_level,
                           GRID *grid,
                           MATERIAL_PROPERTY *mat,
                           FIELD_VARIABLES *variables,
                           SOLVER_OPTIONS *sol,
                           LOADING_STEPS *load,
                           COMMUNICATION_STRUCTURE *com,
                           CRPL *crpl,             /**< Crystal plasticity stuff */
                           double GNOR,            /**< should be local variable. */
                           double nor1,            /**< should be local variable. */
                           long nt,                /**< _DEPRECATED_ */
                           MPI_Comm mpi_comm,
                           const double VVolume,   /**< original volume of the domain */
                           const PGFem3D_opt *opts /**< structure of options */)
{
  return Newton_Raphson(print_level,sol->n_step,
                        grid->ne,
                        grid->n_be,
                        grid->nn,
                        variables->ndofn,
                        variables->ndofd,
                        variables->npres,
                        sol->tim,
                        sol->times,
                        sol->nor_min,
                        sol->dt_np1,
                        grid->element,
                        grid->b_elems,
                        grid->node,
                        load->sup,
                        load->sup_defl,
                        mat->hommat,
                        mat->matgeom,
                        variables->sig,
                        variables->eps,
                        com->Ap,
                        com->Ai,
                        variables->u_np1,
                        variables->f,
                        variables->d_u,
                        variables->dd_u,
                        variables->R,
                        variables->f_defl,
                        variables->RR,
                        variables->f_u,
                        variables->RRn,
                        crpl,
                        sol->stab,
                        grid->nce,
                        grid->coel,
                        sol->FNR,
                        variables->pores,
                        sol->PGFEM_hypre,
                        variables->BS_x,
                        variables->BS_f,
                        variables->BS_RR,
                        sol->gama,
                        GNOR,
                        nor1,
                        sol->err,
                        variables->BS_f_u,
                        com->DomDof,
                        com->comm,
                        com->GDof,
                        nt,
                        sol->iter_max,
                        variables->NORM,
                        com->nbndel,
                        com->bndel,
                        mpi_comm,
                        VVolume,
                        opts,
                        sol->microscale,
                        sol->alpha,
                        variables->u_n,
                        variables->u_nm1);
}
double Newton_Raphson (const int print_level,
                       int *n_step,
                       long ne,
                       int n_be,
                       long nn,
                       long ndofn,
                       long ndofd,
                       long npres,
                       long tim,
                       double *times,
                       double nor_min,
                       double dt,
                       ELEMENT *elem,
                       BOUNDING_ELEMENT *b_elems,
                       NODE *node,
                       SUPP sup,
                       double *sup_defl,
                       HOMMAT *hommat,
                       MATGEOM matgeom,
                       SIG *sig_e,
                       EPS *eps,
                       int *Ap,
                       int *Ai,
                       double *r, /**< total displacement, i.e. for TL */
                       double *f,
                       double *d_r,
                       double *rr,
                       double *R,
                       double *f_defl,
                       double *RR,
                       double *f_u,
                       double *RRn,
                       CRPL *crpl,
                       double stab,
                       long nce,
                       COEL *coel,
                       long FNR,
                       double *pores,
                       PGFEM_HYPRE_solve_info *PGFEM_hypre,
                       double *BS_x,
                       double *BS_f,
                       double *BS_RR,
                       double gama,
                       double GNOR,
                       double nor1,
                       double err,
                       double *BS_f_u,
                       long *DomDof,
                       COMMUN comm,
                       int GDof,
                       long nt,
                       long iter_max,
                       double *NORM,
                       long nbndel,
                       long *bndel,
                       MPI_Comm mpi_comm,
                       const double VVolume,
                       const PGFem3D_opt *opts,
                       void *microscale,
                       double alpha_alpha,
                       double *r_n,
                       double *r_n_1)
{
  double t = times[tim+1];
  double dts[2];

  if(tim==0)
    dts[DT_N] = times[tim+1] - times[tim];
  else  
    dts[DT_N] = times[tim] - times[tim-1];
    
  dts[DT_NP1] = times[tim+1] - times[tim];  
  
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
  
  /* MPI stuff */
  int nproc,myrank;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);
  
  switch(opts->analysis_type){
    case STABILIZED:
    case MINI:
    case MINI_3F:
      ndofn = 4;
      break;
    default:
      break;
  }
  
  *n_step = 0;
  
  /* SUBDIVISION */
  DIV = ST = GAMA = OME = INFO = ART = 0;
  STEP = 1;
  DT = 0.0;
  ERROR = nor_min;
  iter = 0;
  
  /* introduce imperfection in the displacements (useful for
   * homogeneous deformations in homogeneous materials, i.e. when you
   * should be doing it by hand...) */
  /* if(iter == 0 && tim == 0){ */
  /*   for(int i=0; i< ndofd; i++){ */
  /*     d_r[i] += 0.0003; */
  /*   } */
  /* } */
  
  /* GOTO REST */
  rest:
    fflush(PGFEM_stdout);
    
    if(INFO==1 && opts->solution_scheme_opt[LINE_SEARCH]==0)
    {
      ART = 1;
      if(myrank==0)
        printf("Imposed to use NO Line search [INFO = %ld, ART = %ld]\n", INFO, ART);
    }
    
    if (INFO == 1 && ART == 0){
      
      /* Reset variables */
      switch(opts->analysis_type){
        case FS_CRPL: case FINITE_STRAIN:
          res_fini_def (ne,npres,elem,eps,sig_e,
                  crpl,opts->analysis_type);
          break;
        case STABILIZED:
          res_stab_def (ne,npres,elem,eps,sig_e,stab);
          break;
        case MINI:
          MINI_reset(elem,ne,npres,sig_e);
          break;
        case MINI_3F:
          MINI_3f_reset(elem,ne,npres,4,sig_e,eps);
          break;
        case CM:
          constitutive_model_reset_state(eps, ne, elem);
          break;
        default: break;
      }
      
      for (i=0;i<ndofd;i++) {
        rr[i] = d_r[i] = f_defl[i] = f[i] = 0.0;
      }
      
      if (myrank == 0 ) PGFEM_printf("\n** Try without LINE SEARCH **\n\n");
      
      ART = 1;
    } else {
      
      subdivision (INFO,&dt,&STEP,&DIV,tim,times,&ST,ne,
                   ndofn,ndofd,npres,elem,crpl,eps,sig_e,
                   sup,sup_defl,rr,d_r,f_defl,f,RRn,R,&GAMA,
                   &DT,&OME,stab,iter,iter_max,alpha,mpi_comm,
                   opts->analysis_type);
    
      dts[DT_NP1] = dt; // update dt_np1
                        // dt_n is updated when Newton Raphson is completed without error

      gam = ART = 0;
      
    }
    
    /* recompute the microscale tangent if restart due to error. */
    if(INFO == 1 && DEBUG_MULTISCALE_SERVER && microscale != NULL){
      /* zero the macroscale tangent */
      ZeroHypreK(PGFEM_hypre,Ai,DomDof[myrank]);
      
      /* start the microscale jobs. Do not compute equilibrium. Use
       sol->d_r (no displacement increments) for displacement dof
       vector */
      MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) microscale;
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
      bounding_element_communicate_damage(n_be,b_elems,ne,eps,mpi_comm);
      if (periodic == 1){
        double factor = (times[tim] + (DIV+1)*dt)/times[nt] * eps[0].load;
        /* Plane strain deviatoric tension */
        if (eps[0].type == 1){
          eps[0].F[0][0] = times[nt]/(times[nt]
                  - (times[tim] + (DIV+1)*dt)*eps[0].load);
          eps[0].F[1][1] = (1. - factor);
          eps[0].F[2][2] = 1.;
        }
        /* Simple shear */
        if (eps[0].type == 2){
          eps[0].F[0][0] = 1.;
          eps[0].F[1][1] = 1.;
          eps[0].F[2][2] = 1.;
          eps[0].F[0][1] = factor;
        }
        /* Deviatoric tension */
        if (eps[0].type == 3){
          eps[0].F[0][0] = 1./((1. - factor)*(1. - factor));
          eps[0].F[1][1] = (1. - (times[tim] + (DIV+1)*dt)/times[nt] * eps[0].load1);
          eps[0].F[2][2] = (1. - (times[tim] + (DIV+1)*dt)/times[nt] * eps[0].load1);
        }
        
        if (myrank == 0){
          PGFEM_printf("The deformation gradient F\n");
          for (i=0;i<3;i++){
            for (j=0;j<3;j++){
              PGFEM_printf("%12.12f  ",eps[0].F[i][j]);
            }
            PGFEM_printf("\n");
          }
        }
        
        /* Residuals */
        nulld (f_u,ndofd);
        INFO = fd_residuals (f_u,ne,n_be,ndofn,npres,d_r,r,node,elem,b_elems,
                matgeom,hommat,sup,eps,sig_e,
                nor_min,crpl,dts,t,stab,nce,coel,mpi_comm,opts,alpha_alpha,r_n,r_n_1);
        for (i=0;i<ndofd;i++){
          f[i] = - f_u[i];
          R[i] = RR[i] = 0.0;
        }
        
      }/* end periodic */
      else{
        
        for (i=0;i<sup->npd;i++)
          sup->defl_d[i] = sup_defl[i]/STEP;
        
        /* Compute macro interface deformation gradient */
        if(sup->multi_scale){
          /* INFO = compute_interface_macro_jump_u(macro_jump_u,sup); */
          /* compute_interface_macro_grad_u(sup->F0,sup->lc,macro_jump_u, */
          /* 			       sup->N0); */
          
          INFO = compute_macro_grad_u(sup->F0,sup,opts->analysis_type);
          
          if(INFO != 0){
            PGFEM_printerr("[%d] ERROR: not enough prescribed displacements"
                    " for interface multiscale modeling!\n"
                    "Must have at least six (6) prescribed displacements.\n"
                    "Check input and try again.\n",myrank);
            PGFEM_Comm_code_abort(mpi_comm,0);
          }
          nulld (f_u,ndofd);
          vol_damage_int_alg(ne,ndofn,d_r,r,elem,node,
                             hommat,sup,dt,iter,mpi_comm,
                             eps,sig_e,&max_damage,&dissipation,
                             opts->analysis_type);
          
          fd_residuals (f_u,ne,n_be,ndofn,npres,d_r,r,node,elem,b_elems,
                        matgeom,hommat,sup,eps,sig_e,
                        nor_min,crpl,dts,t,stab,nce,coel,mpi_comm,opts,alpha_alpha,r_n,r_n_1);
        } else {
          nulld (f_u,ndofd);
        }
        
        /*  NODE (PRESCRIBED DEFLECTION)
    - SUPPORT COORDINATES generation of the load vector  */
        nulld (f_defl,ndofd);
        load_vec_node_defl (f_defl,ne,ndofn,elem,b_elems,node,hommat,
                            matgeom,sup,npres,nor_min,
                            sig_e,eps,dt,crpl,stab,r,r_n,opts,alpha_alpha);
        
        /* Generate the load and vectors */
        for (i=0;i<ndofd;i++)  {
          f[i] = R[i]/STEP - f_defl[i] - f_u[i];
          //sum += f[i];
          RR[i] = RRn[i] + R[i]/STEP*(DIV+1);
        }
      }
      
      /* Transform LOCAL load vector to GLOBAL */
      LToG (f,BS_f,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
      
      /* Transform LOCAL load vector to GLOBAL */
      LToG (RR,BS_RR,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
      
      iter = 0;
      nor = nor2 = GNOR = 10.0;

      while (nor > ERROR
              && nor2 > zero_tol
              /* && nor2 > ERROR*ERROR */ /* norm of the residual < error^2
               * MM 9/26/2012*/
              ){
        
        if (FNR == 1 || (FNR == 0 && iter == 0)){
          
          assert(opts->solverpackage == HYPRE);
          
          /* Null the matrix (if not doing multiscale)*/
          if(microscale == NULL){
            ZeroHypreK(PGFEM_hypre,Ai,DomDof[myrank]);
          }
          
          INFO = stiffmat_fd (Ap,Ai,ne,n_be,ndofn,elem,b_elems,nbndel,bndel,
                              node,hommat,matgeom,sig_e,eps,d_r,r,npres,sup,
                              iter,nor_min,dt,crpl,stab,nce,coel,0,0.0,f_u,
                              myrank,nproc,DomDof,GDof,
                              comm,mpi_comm,PGFEM_hypre,opts,alpha_alpha,r_n,r_n_1);
          
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
          if(DEBUG_MULTISCALE_SERVER && microscale != NULL){
            ART = 1;
            /* complete any jobs before assembly */
            MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) microscale;
            pgf_FE2_macro_client_recv_jobs(ctx->client,ctx->macro,&max_substep);
          }
          
          /* Matrix assmbly */
          INFO = HYPRE_IJMatrixAssemble(PGFEM_hypre->hypre_k);
          
        }
        
        /*=== Solve the system of equations ===*/
        SOLVER_INFO s_info;
        solve_time += solve_system(opts,BS_f,BS_x,tim,iter,DomDof,&s_info,
                                   PGFEM_hypre,mpi_comm);
        if(myrank == 0){
          solve_system_check_error(PGFEM_stdout,s_info);
        }
        BS_nor = s_info.res_norm;
        BS_iter = s_info.n_iter;
        
        /* Check for correct solution */
        if (BS_nor > 500.*err || BS_iter < 0) {
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
        GToL (BS_x,rr,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
        
        /* LINE SEARCH */
        tmp  = ss (BS_f,BS_f,DomDof[myrank]);
        MPI_Allreduce(&tmp,&Gss_temp,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
        LS1 = 1./2.*Gss_temp;
        
        /* Pressure and volume change THETA */
        
        switch(opts->analysis_type){
          case FS_CRPL:
          case FINITE_STRAIN:
            press_theta (ne,ndofn,npres,elem,node,d_r,rr,sup,matgeom,
                         hommat,eps,sig_e,iter,nor_min,dt,crpl,opts);
            break;
          case MINI:
            MINI_update_bubble(elem,ne,node,ndofn,sup,
                               eps,sig_e,hommat,d_r,rr,iter);
            break;
          case MINI_3F:
            MINI_3f_update_bubble(elem,ne,node,ndofn,sup,
                                  eps,sig_e,hommat,d_r,rr,iter);
            break;
          case TF:
            update_3f(ne,ndofn,npres,d_r,r,rr,node,elem,hommat,sup,eps,sig_e,
                      dt,t,mpi_comm,opts,alpha_alpha,r_n,r_n_1);
            
            break;
          default:
            break;
        }
        
        /*************************/
        /* INTEGRATION ALGORITHM */
        /*************************/
        if (opts->analysis_type == FS_CRPL) {
          INFO = integration_alg (ne,ndofn,ndofd,npres,crpl,elem,
                                  node,d_r,rr,sup,matgeom,hommat,
                                  eps,sig_e,tim,iter,dt,nor_min,STEP,0,opts);
          
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
        for (i=0;i<ndofd;i++) {
          f[i] = d_r[i] + rr[i];
          f_u[i] = 0.0;
        }
        
        INFO = vol_damage_int_alg(ne,ndofn,f,r,elem,node,
                                  hommat,sup,dt,iter,mpi_comm,
                                  eps,sig_e,&max_damage,&dissipation,
                                  opts->analysis_type);
        
        bounding_element_communicate_damage(n_be,b_elems,ne,eps,mpi_comm);
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
        if(DEBUG_MULTISCALE_SERVER && microscale != NULL){
          /* zero the macroscale tangent */
          ZeroHypreK(PGFEM_hypre,Ai,DomDof[myrank]);
          
          /* start the microscale jobs */
          MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) microscale;
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
        INFO = fd_residuals (f_u,ne,n_be,ndofn,npres,f,r,node,elem,
                             b_elems,matgeom,hommat,sup,
                             eps,sig_e,nor_min,crpl,dts,t,stab,
                             nce,coel,mpi_comm,opts,alpha_alpha,r_n,r_n_1);
        
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
        
        if(DEBUG_MULTISCALE_SERVER && microscale != NULL){
          /* print_array_d(PGFEM_stdout,f_u,ndofd,1,ndofd); */
          MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) microscale;
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
        LToG (f_u,BS_f_u,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
        
        /* Compute Euclidian norm */
        for (i=0;i<DomDof[myrank];i++)
          BS_f[i] = BS_RR[i] - BS_f_u[i];
        
        nor  = ss (BS_f,BS_f,DomDof[myrank]);
        MPI_Allreduce(&nor,&GNOR,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
        nor2 = nor = sqrt(GNOR);
        
        if ((tim == 0 && iter == 0)
        || (iter == 0 && *NORM < ERROR)){ /* Reset *NORM if
         * less than convergence tolerance
         * MM 6/27/2012*/
          /* take maximum */
          if(nor > *NORM)	*NORM = nor;
        }
        
        
        /* THIS IS GLOBAL-LOCAL TOLERANCE */
        nor1 = *NORM;
        
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
          INFO = LINE_S3 (&nor,&nor2,&gama,nor1,NOR,LS1,iter,f_u,
                          ne,n_be,ndofd,ndofn,npres,d_r,r,node,elem,b_elems,
                          matgeom,hommat,sup,eps,sig_e,nor_min,crpl,dts,t,
                          stab,nce,coel,f,rr,RR,tim,
                          /*GNOD *gnod,GEEL *geel,*/
                          BS_f,BS_RR,BS_f_u,DomDof,comm,GDof,STEP,mpi_comm,
                          &max_damage, &dissipation, opts,alpha_alpha,r_n,r_n_1);
          
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
              gama = 0.0;
              INFO = 0;
              break; /* take previous as converged value */
              
            } else { /* MM 10/2/2012 not converged, subdivide */
              INFO = 1;
              goto rest;
            }
          } /* end ERROR */
          
          if (gama != 1.0) gam = 1;
          
        } else {       /* no line search */
          gama = 1.0;
        }
        
        /* Total deformation increment */
        NOR = nor;
        for (i=0;i<ndofd;i++){
          d_r[i] += gama*rr[i];
          rr[i] *= gama;
        }
        
        
        /* Compute the energy norm E_norm = abs(R ddu/(R0 ddu0)) */
        LToG(rr,BS_x,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
        enorm = ss(BS_f,BS_x,DomDof[myrank]);
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
            MINI_check_resid(ndofn,ne,elem,node,hommat,eps,
                             sig_e,d_r,sup,RR,DomDof,ndofd,
                             GDof,comm,mpi_comm);
          }
          if(opts->analysis_type == MINI_3F){
            MINI_3f_check_resid(ndofn,ne,elem,node,hommat,eps,
                                sig_e,d_r,sup,RR,DomDof,ndofd,
                                GDof,comm,mpi_comm);
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
        if (iter > (iter_max-1) && nor > ERROR) {
          if (nor > 20.*ERROR || iter > (iter_max + 2)) {
            if (nor < 5.*ERROR) {
              if (myrank == 0)
                PGFEM_printf ("I will take it\n");
              nor = 0.0;
            } else {
              if (myrank == 0)
                PGFEM_printf ("Error in the iteration : iter > iter_max\n");
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
      if (opts->analysis_type == CM) {
        cm_get_subdivision_parameter(&alpha, ne, elem, eps, dt);
      } else {
        alpha = max_damage/max_damage_per_step;
      }
      MPI_Allreduce(MPI_IN_PLACE,&alpha,1,MPI_DOUBLE,MPI_MAX,mpi_comm);
      if(myrank == 0){
        PGFEM_printf("Damage thresh alpha: %f (wmax: %f)\n"
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
      (*n_step) ++;
      
      /* /\* turn off line search *\/ */
      /* ART = 1; */
      
      /* microscale update. Overlay microscale update (using f = r+d_r)
       with macroscale update */
      if(DEBUG_MULTISCALE_SERVER && microscale != NULL){
        /* start the microscale jobs */
        MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) microscale;
        pgf_FE2_macro_client_rebalance_servers(ctx->client,ctx->mpi_comm,
                                               FE2_REBALANCE_NONE);
        
        double tnp1 = 0;
        set_time_micro(tim,times,dt,DIV,&tnp1);
        pgf_FE2_macro_client_send_jobs(ctx->client,ctx->mpi_comm,ctx->macro,
                                       JOB_UPDATE);
        set_time_macro(tim,times,tnp1);
      }
      
      /* increment coheisve elements */
      if(opts->cohesive){
        increment_cohesive_elements(nce,coel,pores,node,sup,d_r);
      }
      
      /* Finite deformations increment */
      switch(opts->analysis_type){
        case FS_CRPL:
        case FINITE_STRAIN:
          fd_increment (ne,nn,ndofn,npres,matgeom,hommat,
                        elem,node,sup,eps,sig_e,d_r,r,
                        nor_min,crpl,dt,nce,coel,pores,mpi_comm,
                        VVolume,opts);
          break;
        case STABILIZED:
          st_increment (ne,nn,ndofn,ndofd,matgeom,hommat,
                        elem,node,sup,eps,sig_e,d_r,r,
                        nor_min,stab,dt,nce,coel,pores,mpi_comm,
                        opts->cohesive);
          break;
        case MINI:
          MINI_increment(elem,ne,node,nn,ndofn,
                         sup,eps,sig_e,hommat,d_r,mpi_comm);
          break;
        case MINI_3F:
          MINI_3f_increment(elem,ne,node,nn,ndofn,
                            sup,eps,sig_e,hommat,d_r,mpi_comm);
          break;
        case DISP:
          DISP_increment(elem,ne,node,nn,ndofn,sup,eps,
                         sig_e,hommat,d_r,r,mpi_comm);
          break;
        case TF:
          update_3f_state_variables(ne,ndofn,npres,d_r,r,node,elem,hommat,sup,eps,sig_e,
                                    dt,t,mpi_comm);
          break;
        case CM:
        {
          switch(opts->cm)
          {
            case HYPER_ELASTICITY: case DISP:
              DISP_increment(elem,ne,node,nn,ndofn,sup,eps,
                             sig_e,hommat,d_r,r,mpi_comm);
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
      
      /* Add deformation increment into displacement vector */
      vvplus (r,d_r,ndofd);
      
      /* update previous time step values, r_n and r_n_1 from current
       time step values r*/
      /* For FE2, these vectors are not allocated and are passed as NULL */
      if ( r_n_1 != NULL && r_n != NULL )
      {
        for(long a = 0; a<nn; a++)
        {
          for(long b = 0; b<ndofn; b++)
          {
            r_n_1[a*ndofn + b] = r_n[a*ndofn + b];
            long id = node[a].id[b];
            if(opts->analysis_type==CM && opts->cm==UPDATED_LAGRANGIAN)
              // Updated Lagrangian
            {
              if(id>0)
                r_n[a*ndofn + b] = d_r[id-1];
              else
              {
                if(id==0)
                  r_n[a*ndofn + b] = 0.0;
                else
                  r_n[a*ndofn + b] = sup->defl_d[abs(id)-1];
              }
            }
            else
              // Total Lagrangian: DISP, TF
            {
              if(id>0)
                r_n[a*ndofn + b] = r[id-1];
              else
              {
                if(id==0)
                  r_n[a*ndofn + b] = 0.0;
                else
                  r_n[a*ndofn + b] = sup->defl[abs(id)-1] + sup->defl_d[abs(id)-1];
              }
            }
          }
        }
      }
      
      /* update of internal variables and nodal coordinates */
      if(opts->analysis_type==CM)
      {
        switch(opts->cm){
          case UPDATED_LAGRANGIAN:
          case TOTAL_LAGRANGIAN:
            constitutive_model_update_time_steps_test(elem,node,eps,ne,nn,
                                                      ndofn,r_n,dt,opts->cm);
            break;
          case MIXED_ANALYSIS_MODE:
            constitutive_model_update_time_steps_test(elem,node,eps,ne,nn,
                                                      ndofn,r_n,dt,1 /* TL */);
            break;
          default: break;
        }
      }
      
      if(opts->analysis_type==TF)
      {
        int nVol = 1;
        for (int e=0;e<ne;e++)
        {
          if(npres==1)
          {
            eps[e].d_T[2] = eps[e].d_T[1];
            eps[e].d_T[1] = eps[e].d_T[0];
            
          }
          for(int a=0; a<nVol; a++)
          {
            eps[e].T[a*3+2] = eps[e].T[a*3+1];
            eps[e].T[a*3+1] = eps[e].T[a*3+0];
          }
        }
      }
      
      // Update time steps
      dts[DT_N] = dts[DT_NP1];
      
      /* Null prescribed increment deformation */
      for (i=0;i<sup->npd;i++){
        sup->defl[i] += sup->defl_d[i];
        sup->defl_d[i] = 0.0;
      }
      
      for (i=0;i<ndofd;i++) {
        d_r[i] = 0.0;
        rr[i] = 0.0;
        f_defl[i] = 0.0;
        f[i] = 0.0;
      }
      
      /* finish microscale update */
      if(DEBUG_MULTISCALE_SERVER && microscale != NULL){
        /* start the microscale jobs */
        MS_SERVER_CTX *ctx = (MS_SERVER_CTX *) microscale;
        pgf_FE2_macro_client_recv_jobs(ctx->client,ctx->macro,&max_substep);
      }
      
      /************* TEST THE UPDATE FROM N TO N+1  *************/
      if(NR_UPDATE || PFEM_DEBUG || PFEM_DEBUG_ALL){
        for (i=0;i<ndofd;i++) {f_u[i] = 0.0; d_r[i] = 0.0;}
        
        if(opts->analysis_type == MINI){
          MINI_check_resid(ndofn,ne,elem,node,hommat,eps,
                           sig_e,d_r,sup,RR,DomDof,ndofd,
                           GDof,comm,mpi_comm);
        }
        if(opts->analysis_type == MINI_3F){
          MINI_3f_check_resid(ndofn,ne,elem,node,hommat,eps,
                              sig_e,d_r,sup,RR,DomDof,ndofd,
                              GDof,comm,mpi_comm);
        }
        fd_residuals (f_u,ne,n_be,ndofn,npres,d_r,r,node,elem,b_elems,
                      matgeom,hommat,sup,eps,sig_e,
                      nor_min,crpl,dts,t,stab,nce,coel,mpi_comm,opts,alpha_alpha,r_n,r_n_1);
        
        for (i=0;i<ndofd;i++) f[i] = RR[i] - f_u[i];
        /* print_array_d(stdout,RR,ndofd,1,ndofd); */
        
        LToG(f,BS_f,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
        nor = ss(BS_f,BS_f,DomDof[myrank]);
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
            VTK_print_master(opts->opath,fname,*n_step,nproc,opts);
          }
          VTK_print_vtu(opts->opath,fname,*n_step,myrank,ne,nn,node,
                        elem,sup,r,sig_e,eps,opts);
        }
      }
      
      double *res_trac = aloc1(3);
      bounding_element_compute_resulting_traction(n_be,b_elems,elem,node,
                                                  eps,sig_e,ndofd,DomDof,
                                                  GDof,comm,mpi_comm,
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
    
    if(NR_COMPUTE_REACTIONS && !sup->multi_scale){
      compute_reactions(ne,ndofn,npres,r,node,elem,matgeom,
                        hommat,sup,eps,sig_e,nor_min,crpl,
                        dt,stab,mpi_comm,opts->analysis_type);
    }
    
    times[tim] = times[tim+1] - dts[DT_NP1];
    return (solve_time);
    
}
