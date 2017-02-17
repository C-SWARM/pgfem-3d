#include "Arc_length.h"
#include <sys/time.h> 
#include <sys/resource.h>
#include <assert.h>

#include "PGFEM_io.h"
#include "enumerations.h"
#include "solve_system.h"
#include "ALM.h"
#include "fd_increment.h"
#include "fd_residuals.h"
#include "incl.h"
#include "integration.h"
#include "vol_damage_int_alg.h"
#include "LINE.h"
#include "macroscopic_load_AL.h"
#include "matice.h"
#include "out.h"
#include "press_theta.h"
#include "res_fini_def.h"
#include "stabilized.h"
#include "stiffmat_fd.h"
#include "subdivision.h"
#include "utils.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "displacement_based_element.h"
#include "dynamics.h"
#include "PGFem3D_data_structure.h"


#ifndef ARC_DEBUG
#define ARC_DEBUG 0
#endif

#ifndef ARC_UPDATE
#define ARC_UPDATE 0
#endif

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

#ifndef PFEM_DEBUG_ALL
#define PFEM_DEBUG_ALL 0
#endif

static const int periodic = 0;

/// initialize arc length variable object
/// assign defaults (zoro for single member varialbes and NULL for member arrays and structs
///                  except ARC=1)
/// 
/// \param[in, out] arc an object for arc length analysis containing additional variables
/// \return non-zero on internal error
int arc_length_variable_initialization(ARC_LENGTH_VARIABLES *arc)
{
  int err = 0;
  arc->dt0    = 0.0;
  arc->D_R    = NULL;
  arc->U      = NULL;
  arc->DK     = NULL;
  arc->dR     = NULL;
  arc->BS_d_r = NULL;
  arc->BS_D_R = NULL;
  arc->BS_rr  = NULL;
  arc->BS_R   = NULL;
  arc->BS_U   = NULL;
  arc->BS_DK  = NULL;
  arc->BS_dR  = NULL;
  arc->lm     = 0.0;
  arc->dAL0   = 0.0;
  arc->DET0   = 0.0;
  arc->DLM0   = 0.0;
  arc->DLM    = 0.0;
  arc->AT     = 0;  
  arc->ARC    = 1; // ARC LENGTH METHOD
                   // ARC == 0 :: Crisfield
                   // ARC == 1 :: Simo ONLY ARC LENG METHOD FOR PARALELL
  arc->dALMAX = 0.0;
  arc->ITT    = 0;
  arc->DAL    = 0.0;
  return err;
}

/// construct arc length variable object
/// create memory space for member arrays and structs
/// 
/// \param[in, out] arc an object for arc length analysis containing additional variables
/// \param[in] fv an object containing all field variables
/// \param[in] com an object for communication
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int construct_arc_length_variable(ARC_LENGTH_VARIABLES *arc,
                                  FIELD_VARIABLES *fv,
                                  COMMUNICATION_STRUCTURE *com,
                                  int myrank)
{
  int err = 0;
  arc->D_R    = aloc1(fv->ndofd);
  arc->U      = aloc1(fv->ndofd);
  arc->DK     = aloc1(fv->ndofd);
  arc->dR     = aloc1(fv->ndofd);
  arc->BS_d_r = aloc1(com->DomDof[myrank]);
  arc->BS_D_R = aloc1(com->DomDof[myrank]);
  arc->BS_rr  = aloc1(com->DomDof[myrank]);
  arc->BS_R   = aloc1(com->DomDof[myrank]);
  arc->BS_U   = aloc1(com->DomDof[myrank]);
  arc->BS_DK  = aloc1(com->DomDof[myrank]);
  arc->BS_dR  = aloc1(com->DomDof[myrank]);
  return err;
} 

/// destruct arc length variable object
/// free memory spaces for member arrays and structs
/// 
/// \param[in, out] arc an object for arc length analysis containing additional variables
/// \return non-zero on internal error
int destruct_arc_length_variable(ARC_LENGTH_VARIABLES *arc)
{
  int err = 0;
  if(NULL != arc->D_R)    free(arc->D_R);
  if(NULL != arc->U)      free(arc->U);
  if(NULL != arc->DK)     free(arc->DK);
  if(NULL != arc->dR)     free(arc->dR);
  if(NULL != arc->BS_d_r) free(arc->BS_d_r);
  if(NULL != arc->BS_D_R) free(arc->BS_D_R);
  if(NULL != arc->BS_rr)  free(arc->BS_rr);
  if(NULL != arc->BS_R)   free(arc->BS_R);
  if(NULL != arc->BS_U)   free(arc->BS_U);
  if(NULL != arc->BS_DK)  free(arc->BS_DK);
  if(NULL != arc->BS_dR)  free(arc->BS_dR);
  err += arc_length_variable_initialization(arc);
  return err;
} 

double Arc_length(long ne,
		   int n_be,
		   long nn,
		   long ndofn,
		   long ndofd,
		   long npres,
		   long nt,
		   long tim,
		   double *times,
		   double nor_min,
		   long iter_max,
		   double dt,
		   double dt0,
		   ELEMENT *elem,
		   BOUNDING_ELEMENT *b_elems,
		   long nbndel,
		   long *bndel,
		   NODE *node,
		   SUPP sup,
		   double *sup_defl,
		   HOMMAT *hommat,
		   MATGEOM matgeom,
		   SIG *sig_e,
		   EPS *eps,
		   int *Ap,
		   int *Ai,
		   PGFEM_HYPRE_solve_info *PGFEM_hypre,
		   double *RRn,
		   double *f_defl,
		   CRPL *crpl,
		   double stab,
		   long nce,
		   COEL *coel,
		   double *r,
		   double *f,
		   double *d_r,
		   double *D_R,
		   double *rr,
		   double *R,
		   double *RR,
		   double *f_u,
		   double *U,
		   double *DK,
		   double *dR,
		   double *BS_f,
		   double *BS_d_r,
		   double *BS_D_R,
		   double *BS_rr,
		   double *BS_R,
		   double *BS_RR,
		   double *BS_f_u,
		   double *BS_U,
		   double *BS_DK,
		   double *BS_dR,
		   long FNR,
		   double lm,
		   double dAL0,
		   double *DET0,
		   double *DLM0,
		   double *dlmdlm,
		   long gr2,
		   long gr4,
		   SIG *sig_n,
		   char *out_dat,
		   long *print,
		   long *AT,
		   long ARC,
		   double dALMAX,
		   long *ITT,
		   double *DAL,
		   double *pores,
		   /*long nge,
		     GEEL *geel,
		     long ngn,
		     GNOD *gnod*/
		   long *DomDof,
		   int GDof,
		   COMMUN comm,
		   double err,
		   double *NORM,
		   MPI_Comm mpi_comm,
		   const double VVolume,
		   const PGFem3D_opt *opts,
		   const int mp_id)
/*
  ARC == 0 :: Crisfield
  ARC == 1 :: Simo - ONLY THIS CAN BE USED WITH PARALLEL FRAMEWORK
  
  TYPE = 0 :: Time subdivision scheme
  TYPE = 1 :: Arc-length subdivision procedure

*/
{
  double dts[2];
  dts[0] = dts[1] = dt;
  double t = 0.0;
  double *r_n = NULL;
  double *r_n_1 = NULL;
  double alpha_alpha = 0.0;
  
  double nor,nor1,nor2,dlm,dlm0,DLM,DET=0.0,dAL;
  double DT,DDLM,ddlm,ERROR,LS1,gama,pdt,tmp,nor3;
  long iter,INFO,i,j,STEP,N,M,DIV,ST,GAMA,OME,FI,ART,gam,TYPE,GInfo;
  char *error[]={"inf","-inf","nan"},str1[500],jmeno[50];
  FILE *out;
  struct rusage usage;

  int nproc,myrank;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);

  /* damage substep criteria */
  const double max_damage_per_step = 0.05;
  double max_damage = 0.0;
  double alpha = 0.0;

  /* dissipation */
  double dissipation = 0.0;

  /* BlockSolve 95 */
  double BS_nor=0.0;
  int BS_iter;

  /* TEMPORARY bounding element compile testing */
  /* int n_be = 0; */
  /* BOUNDING_ELEMENT *b_elems = NULL; */

  TYPE = 1;
  ERROR = nor_min;
  pdt = dt;
 
  switch(opts->analysis_type){
  case STABILIZED:
  case MINI:
  case MINI_3F:
    ndofn = 4;
    break;
  default:
    break;
  }

  /* SUBDIVISION */
  DT = DDLM = ddlm = dlm = 0.0;
  DIV = ST = GAMA = OME = INFO = ART = gam = FI = iter = 0;
  STEP = 1;
 
 rest:
  /* /\* Force load stepping for ALL restarts *\/ */
  /* ART = 0; */
  /* if (INFO == 1 && ART == 0){ */
    
  /*   /\* Reset variables *\/ */
  /*   if (opts->analysis_type == 2 || opts->analysis_type == 3) */
  /*     res_fini_def (ne,npres,elem,eps,sig_e,crpl); */

  /*   if (opts->analysis_type == 4) */
  /*     res_stab_def (ne,npres,elem,eps,sig_e,stab); */
    
  /*   if(opts->analysis_type == 5) /\* P1+B/P1 *\/ */
  /*     MINI_reset(elem,ne,npres,sig_e); */

  /*   if(opts->analysis_type == 6) /\* P1+B/P1/V1 *\/ */
  /*     MINI_3f_reset(elem,ne,npres,1,sig_e,eps); */

  /*   for (i=0;i<ndofd;i++) { */
  /*     rr[i] = d_r[i] = D_R[i] = f_defl[i] = f[i] = 0.0; */
  /*   } */
    
  /*   if (myrank == 0) PGFEM_printf ("\n** Trying new root selection **\n\n"); */
    
  /*   ART = 1; */
  /* } */
  /* else { */
  /* if (myrank == 0) PGFEM_printf ("\n** Subdividing **\n\n"); */
  dlm0 = subdiv_arc (INFO,&dt,dt0,&STEP,&DIV,tim,times,&ST,
		     ne,ndofd,npres,elem,crpl,eps,sig_e,sup,
		     sup_defl,rr,d_r,D_R,f_defl,f,&GAMA,&DT,
		     &OME,stab,dAL0,*DAL,dALMAX,nor_min,dlm0,
		     ITT,iter,iter_max,TYPE,mpi_comm,opts->analysis_type);
	
	dts[DT_NP1] = dt;	     
    
  if (periodic == 1) ART = 1;
  else               ART = 0;
    
  gam = 0;
  if (TYPE == 1) dAL = dlm0;
  /* } */
  
  dlm0 *= *DLM0/fabs(*DLM0);

  INFO = 0;
  while (STEP > DIV){
    if ( (STEP > 1 || ST == 1) && myrank == 0){
      PGFEM_printf ("\nSTEP = %ld :: NS =  %ld || Time %f | dt = %10.10f\n",
	      DIV,STEP,times[tim+1],dt);
    }

    if(myrank == 0){
      PGFEM_printf("dlm0 = %12.12e\n",dlm0);
    }

    iter = 0; 
    if (periodic == 1){
      if (TYPE == 1) {
	nor = dlm0;
	dlm0 = 0.0;
      }
      /* macroscopic load */
      macroscopic_load_AL (eps[0].type,lm+dlm0,eps);
    }/* end periodic */
    
    if (periodic == 1) nulld (R,ndofd);
    
    assert(opts->solverpackage == HYPRE);
    /* Null the matrix */
    ZeroHypreK(PGFEM_hypre,Ai,DomDof[myrank]);
    stiffmat_fd (Ap,Ai,ne,n_be,ndofn,elem,b_elems,nbndel,bndel,
		 node,hommat,matgeom,sig_e,eps,d_r,r,npres,sup,iter,
		 nor_min,dt,crpl,stab,nce,coel,FNR,lm+dlm0,R,myrank,
		 nproc,DomDof,GDof,comm,mpi_comm,PGFEM_hypre,opts,alpha_alpha,r_n,r_n_1,
		 mp_id);

    /* Assemble the matrix */
    HYPRE_IJMatrixAssemble(PGFEM_hypre->hypre_k);

    
    if (periodic == 1 && TYPE == 1) dlm0 = nor;
    /* Transform LOCAL load vector to GLOBAL */
    if (periodic == 1){
      LToG (R,BS_R,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
    }

    /*** SOLVE THE SYSTEM FOR DIRECTION ***/
    {
      SOLVER_INFO s_info;
      solve_system(opts,BS_R,BS_rr,tim,iter,DomDof,
		   &s_info,PGFEM_hypre,mpi_comm);
      if(myrank == 0){
	solve_system_check_error(stdout,s_info);
      }
      BS_nor = s_info.res_norm;
      BS_iter = s_info.n_iter;
    }

    if(ARC_DEBUG && myrank == 0){
      PGFEM_printf("Completed solve I!\n");
    }

    if (BS_nor > 100.*err || BS_iter < 0) {
      INFO = 1; 
      ART = 1;
      if (myrank == 0){
	PGFEM_printf("ERROR in the BSpar_solve : nor = %8.8e || iter = %d\n",
	       BS_nor,BS_iter);
      }
      goto rest;
    }
     
    sprintf (str1,"%f",BS_nor);
    for (N=0;N<3;N++){
      M = 10;
      M = strcmp(error[N],str1);
      if (M == 0) {
	if (myrank == 0)
	  PGFEM_printf("ERROR in the BSpar_solve : nor = %s\n",error[N]);
	INFO = 1;
	ART = 1;
	goto rest;
      }
    }

    if (ARC_DEBUG && myrank == 0) {
      PGFEM_printf("Starting GToL... ");
    }    

    /* Transform GLOBAL displacement vector to LOCAL */
    GToL (BS_rr,rr,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);

    /* Null R for periodic */
    if (periodic == 1) {nulld (R,ndofd); nulld (BS_R,DomDof[myrank]);}
    
    /* Scaling parameter DK */
    if (tim == 0 && iter == 0){
      if (ARC == 0){ /* diag_KK (ndofd,k,Ap,Ai,DK); */
	for (i=0;i<DomDof[myrank];i++) BS_DK[i] = 1.0;
      }
      if (ARC == 1){ /* diag_KK (ndofd,k,Ap,Ai,DK); */
	for (i=0;i<DomDof[myrank];i++) BS_DK[i] = 1.0;
      }
    }

    /* First load multiplier */
    if (ARC == 0) {
      if (TYPE == 0)
	dAL = d_ALM2 (ndofd,BS_rr,BS_R,BS_DK,dlm0);
      dlm = d_lam_ALM2 (ndofd,BS_rr,BS_R,BS_DK,dAL,DET,*DET0,
			*DLM0,nor_min,dR);
    } 
    if (ARC == 1) {
      if (TYPE == 0)
	dAL = d_ALM4 (ndofd,BS_rr,BS_DK,dlm0,DomDof,mpi_comm);
      dlm = d_lam_ALM4 (ndofd,BS_rr,BS_DK,BS_dR,dAL,DomDof,mpi_comm);
    }

    /* First load multiplier */
    dlm0 = dlm;
    if (myrank == 0)
      PGFEM_printf("dAL = %12.15e :: dALmax = %12.15e || dlm0  = %12.15e\n\n",
	     dAL,dALMAX,dlm0);
    
    /* Limiting the arc length */
    if (dAL > dALMAX){
      if (myrank == 0)
	PGFEM_printf ("*** Arc length too large: Restart with smaller arc ***\n");
      INFO = 1; ART = 1;
      goto rest;
    }
    
    /* macroscopic load */
    if (periodic == 1) macroscopic_load_AL (eps[0].type,lm+dlm,eps);
    
    /* First displacement increment dr_0 */
    for (i=0;i<ndofd;i++) d_r[i] = dlm*rr[i];
    
    /* INCREMENT: L -> G */
    LToG(d_r,BS_d_r,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
    
    /* Pressure and volume change THETA */
    nulld (f_u,ndofd);
    switch(opts->analysis_type){
    case FS_CRPL:
    case FINITE_STRAIN:
      press_theta (ne,ndofn,npres,elem,node,f_u,d_r,sup,
		   matgeom,hommat,eps,sig_e,iter,nor_min,dt,crpl,opts,mp_id);
      break;
    case MINI:
      MINI_update_bubble(elem,ne,node,ndofn,sup,
			 eps,sig_e,hommat,f_u,d_r,iter,mp_id);
      break;
    case MINI_3F:
      MINI_3f_update_bubble(elem,ne,node,ndofn,sup,
			 eps,sig_e,hommat,f_u,d_r,iter,mp_id);
      break;
    default:
      break;
    }
    
    /*************************/
    /* INTEGRATION ALGORITHM */
    /*************************/
    if (opts->analysis_type == FS_CRPL) {
      INFO = integration_alg (ne,ndofn,ndofd,npres,crpl,elem,
			      node,f_u,d_r,sup,matgeom,hommat,
			      eps,sig_e,tim,iter,dt,nor_min,STEP,0,opts,mp_id); 

      /* Gather INFO from all domains */
      MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
      if (GInfo == 1) {
	INFO = 1;
	goto rest;
      }
    }

    vol_damage_int_alg(ne,ndofn,d_r,r,elem,node,
		       hommat,sup,dt,iter,mpi_comm,
		       eps,sig_e,&max_damage,&dissipation,
		       opts->analysis_type,mp_id);
    
    /* Residuals */
    fd_residuals(f_u,ne,n_be,ndofn,npres,d_r,r,node,elem,b_elems,matgeom,
		  hommat,sup,eps,sig_e,nor_min,crpl,dts,t,stab,
		  nce,coel /*,gnod,geel*/,mpi_comm,opts,alpha_alpha,r_n,r_n_1,
		  mp_id);
    
    /* Compute Euclidian norm */
    for (i=0;i<ndofd;i++) f[i] = (lm + dlm)*R[i] - f_u[i];  
    
    /* BS_f : L->G */
    LToG(f,BS_f,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);

    nor = ss(BS_f,BS_f,DomDof[myrank]);

    /* Gather nor from each domain */
    MPI_Allreduce(&nor,&nor3,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    nor = nor2 = sqrt (nor3);
    nor1 = nor;

    sprintf (str1,"%f",nor);
    for (N=0;N<3;N++){
      M = 10;
      M = strcmp(error[N],str1); 
      if (M == 0) {
	if (myrank == 0)
	  PGFEM_printf("ERROR in the algorithm : nor = %s\n",error[N]);
	INFO = 1;
	goto rest;
      }
    }

    if (tim == 0 && iter == 0) *NORM = nor1;
    
    /* THIS IS THE GLOBAL-LOCAL TOLERANCE */
    nor1 = *NORM;
    
    /* Normalize norm */
    if (nor < 1.e-10) {
      if (myrank == 0)
	PGFEM_printf ("Scale your units : Close to the computer precision\n");
      nor = 0.0;
    } else nor /= nor1;

    /* For very bad convergence problems */
    if (fabs(dlm) < nor_min) {
      nor = nor2;
      ERROR = 10.0*nor_min;
      if (myrank == 0)
	PGFEM_printf ("Reducing and changing error || ERROR = %12.12e\n",ERROR);
    } else ERROR = nor_min;
    
    if (myrank == 0) {
      getrusage (RUSAGE_SELF,&usage);
      PGFEM_printf("(%ld) IT = %d : R = %8.8e :: ||f||/||f0|| ="
	     " [%8.8e] || [%8.8e] :: S %ld.%ld, U %ld.%ld\n",
	     iter,BS_iter,BS_nor,nor,nor2,usage.ru_stime.tv_sec,
	     usage.ru_stime.tv_usec,usage.ru_utime.tv_sec,
	     usage.ru_utime.tv_usec);
    }

    /*************/
    /* iter =  1 */
    /*************/    

    iter = 1; DLM = 0.0;
    while (nor > ERROR){

      /* Null the residual vector */
      nulld (f_u,ndofd);
      
      assert(opts->solverpackage == HYPRE);
      /* Null the matrix */
      ZeroHypreK(PGFEM_hypre,Ai,DomDof[myrank]);
      stiffmat_fd (Ap,Ai,ne,n_be,ndofn,elem,b_elems,nbndel,bndel,
		   node,hommat,matgeom,sig_e,eps,d_r,r,npres,sup,iter,
		   nor_min,dt,crpl,stab,nce,coel,FNR,lm+dlm,f_u,myrank,
		   nproc,DomDof,GDof,comm,mpi_comm,PGFEM_hypre,opts,alpha_alpha,r_n,r_n_1,
		   mp_id);
	
      /* Assemble the matrix */
      HYPRE_IJMatrixAssemble(PGFEM_hypre->hypre_k);
      
      /* Transform unequibriated force: L -> G */
      if (periodic != 1){
	for (i=0;i<DomDof[myrank];i++)
	  BS_f_u[i] = BS_R[i];
      } else {
	LToG(f_u,BS_f_u,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
      }

      /*** SOLVE THE SYSTEM EQUATIONS ***/
      {
	SOLVER_INFO s_info;
	solve_system(opts,BS_f_u,BS_rr,tim,iter,DomDof,
		     &s_info,PGFEM_hypre,mpi_comm);
	if(myrank == 0){
	  solve_system_check_error(stdout,s_info);
	}
	BS_nor = s_info.res_norm;
	BS_iter = s_info.n_iter;
      }

      /* Check for correct solution */
      if (BS_nor > 100.*err || BS_iter < 0) {
	INFO = 1;
	ART = 1;
	if (myrank == 0)
	  PGFEM_printf("ERROR in the BSpar_solve : nor = %8.8e || iter = %d\n",
		 BS_nor,BS_iter);
	goto rest;
      }

      sprintf (str1,"%f",BS_nor);
      for (N=0;N<3;N++){
	M = 10;
	M = strcmp(error[N],str1);
	if (M == 0) {
	  if (myrank == 0)
	    PGFEM_printf("ERROR in the BSpar_solve : nor = %s\n",error[N]);
	  INFO = 1;
	  ART = 1;
	  goto rest;
	}
      }
      

      if (myrank == 0) {
	getrusage (RUSAGE_SELF,&usage);
	PGFEM_printf("(%ld.5) Intermediate Solve Info: IT = %d : R = %8.8e ::"
	       " S %ld.%ld, U %ld.%ld\n",(iter-1),BS_iter,BS_nor,
	       usage.ru_stime.tv_sec,usage.ru_stime.tv_usec,
	       usage.ru_utime.tv_sec,usage.ru_utime.tv_usec);
      }

      /* initial disp. G-> L */
      GToL(BS_rr,rr,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
      
      /* Solve i-th incremnt */
      {
	SOLVER_INFO s_info;
	solve_system(opts,BS_f,BS_U,tim,iter,DomDof,
		     &s_info,PGFEM_hypre,mpi_comm);
	if(myrank == 0){
	  solve_system_check_error(stdout,s_info);
	}
	BS_nor = s_info.res_norm;
	BS_iter = s_info.n_iter;
      }

      /* Check for correct solution */
      if (BS_nor > 100.*err || BS_iter < 0) {
	INFO = 1;
	ART = 1;
	if (myrank == 0)
	  PGFEM_printf("ERROR in the BSpar_solve : nor = %8.8e || iter = %d\n",
		 BS_nor,BS_iter);
	goto rest;
      }

      sprintf (str1,"%f",BS_nor);
      for (N=0;N<3;N++){
	M = 10;
	M = strcmp(error[N],str1);
	if (M == 0) {
	  if (myrank == 0)
	    PGFEM_printf("ERROR in the BSpar_solve : nor = %s\n",error[N]);
	  INFO = 1;
	  ART = 1;
	  goto rest;
	}
      }
      
      /* increment disp. G-> L */
      GToL (BS_U,U,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
      
      /* dlam */
      if (ARC == 0){
	if (*AT == 1 || ART == 0){
	  DLM = D_lam_ALM (ndofd,BS_rr,BS_d_r,BS_U,BS_R,BS_DK,dlm,
			   dAL,DomDof,mpi_comm);
	} else {
	  DLM = D_lam_ALM2 (BS_rr,BS_U,BS_R,BS_DK,dlm,lm,dAL,ne,n_be,ndofd,
			    npres,BS_d_r,r,node,elem,b_elems,matgeom,hommat,
			    sup,eps,sig_e,nor_min,crpl,dt,stab,nce,coel,
			    DomDof,GDof,comm,mpi_comm,opts,mp_id);
	}
      }
      if (ARC == 1)
	DLM = D_lam_ALM4 (ndofd,BS_rr,BS_d_r,BS_U,BS_DK,dlm,dAL,
			  DomDof,mpi_comm);
      
      sprintf (str1,"%f",DLM);
      for (N=0;N<3;N++){
	M = 10;
	M = strcmp(error[N],str1);
	if (M == 0) {
	  if (myrank == 0)
	    PGFEM_printf("Complex root in ARC-LENGHT method\n");
	  INFO = 1;
	  ART = 1;
	  goto rest;
	}
      }
      
      /* Periodic loading */
      if (periodic == 1)
	macroscopic_load_AL (eps[0].type,lm+dlm+DLM,eps);
      
      /* Update deformation */
      for (i=0;i<ndofd;i++)
	D_R[i] = U[i] + DLM*rr[i];
      
      /* LINE SEARCH */
      tmp = ss(BS_f,BS_f,DomDof[myrank]);
      MPI_Allreduce(&tmp,&LS1,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
      LS1 *= 1./2;
      
      /* Pressure and volume change THETA */
      switch(opts->analysis_type){
      case FS_CRPL:
      case FINITE_STRAIN:
	press_theta (ne,ndofn,npres,elem,node,d_r,D_R,sup,matgeom,
		     hommat,eps,sig_e,iter,nor_min,dt,crpl,opts,mp_id);
	break;
      case MINI:
	MINI_update_bubble(elem,ne,node,ndofn,sup,
			   eps,sig_e,hommat,d_r,D_R,iter,mp_id);
	break;
      case MINI_3F:
	MINI_3f_update_bubble(elem,ne,node,ndofn,sup,
			      eps,sig_e,hommat,d_r,D_R,iter,mp_id);
	break;
      default:
	break;
      }

      /*************************/
      /* INTEGRATION ALGORITHM */
      /*************************/
      if (opts->analysis_type == FS_CRPL) {
	INFO = integration_alg (ne,ndofn,ndofd,npres,crpl,elem,
				node,d_r,D_R,sup,matgeom,hommat,
				eps,sig_e,tim,iter,dt,nor_min,STEP,0,opts,mp_id);
      
	/* Gather INFO from all domains */
	MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
	if (GInfo == 1) {
	  INFO = 1;
	  goto rest;
	}
      }
      
      /* Update deformations */
      for (i=0;i<ndofd;i++) {
	f[i] = d_r[i] + D_R[i];
	f_u[i] = 0.0;
      } 

      vol_damage_int_alg(ne,ndofn,f,r,elem,node,
			 hommat,sup,dt,iter,mpi_comm,
			 eps,sig_e,&max_damage,&dissipation,
			 opts->analysis_type,mp_id);
      
      /* Residuals */
      fd_residuals (f_u,ne,n_be,ndofn,npres,f,r,node,elem,b_elems,matgeom,
		    hommat,sup,eps,sig_e,nor_min,crpl,dts,t,stab,
		    nce,coel/*,gnod,geel*/,mpi_comm,opts,alpha_alpha,r_n,r_n_1,
		    mp_id);

      /* Compute Euclidean norm */
      for (i=0;i<ndofd;i++)
	f[i] = (lm + dlm + DLM)*R[i] - f_u[i];
      
      /* Residuals L -> G */
      LToG(f,BS_f,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
      
      nor = ss(BS_f,BS_f,DomDof[myrank]);
      MPI_Allreduce(&nor,&tmp,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
      nor = nor2 = sqrt (tmp);
      nor /= nor1;
      
      sprintf (str1,"%f",nor);
      for (N=0;N<3;N++){
	M = 10;
	M = strcmp(error[N],str1);
	if (M == 0){
	  if (myrank == 0)
	    PGFEM_printf("ERROR in the algorithm : nor = %s\n",error[N]);
	  INFO = 1;
	  goto rest;
	}
      }
      
      /* LINE SEARCH */
      if (ART == 0) {
	INFO = ALINE_S3 (ARC,&DLM,&nor,&nor2,&gama,nor1,LS1,iter,ne,n_be,
			 ndofd,ndofn,npres,tim,nor_min,dts,stab,nce,
			 dlm,lm,dAL,d_r,r,D_R,node,elem,b_elems,matgeom,
			 hommat,sup,eps,sig_e,crpl,coel,f_u,f,R,BS_f,
			 BS_R,BS_D_R,BS_d_r,BS_DK,BS_U,BS_rr,
			 DomDof,GDof,comm,STEP,mpi_comm,&max_damage,
			 &dissipation,opts,mp_id);

	/* Gather INFO from all domains */
	MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
	if (GInfo == 1) {
	  if (myrank == 0)
	    PGFEM_printf ("Error in the Line Search 3 algorithm\n");
	  INFO = 1;
	  goto rest;
	}
	if (gama != 1.0)
	  gam = 1;
      }
      
      /* For very bad convergence problems */
      if (fabs(dlm0) < nor_min)
	nor = nor2;
      
      /* Load multiplied update */
      dlm += DLM;
      
      /* Deformation update */
      for (i=0;i<ndofd;i++)
	d_r[i] += D_R[i];
      
      /* Total increment: L -> G */
      LToG(d_r,BS_d_r,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
      
      if (myrank == 0) {
	getrusage (RUSAGE_SELF,&usage);
	PGFEM_printf("(%ld) IT = %d : R = %8.8e :: "
	       "||f||/||f0|| = [%8.8e] || [%8.8e] ::"
	       " S %ld.%ld, U %ld.%ld\n",
	       iter,BS_iter,BS_nor,nor,nor2,usage.ru_stime.tv_sec,
	       usage.ru_stime.tv_usec,usage.ru_utime.tv_sec,
	       usage.ru_utime.tv_usec);
      }
      
      /* Max number of iterations restart */
      if (iter > (iter_max-1) && nor > ERROR) {
	if (nor > 20.*ERROR || iter > (iter_max + 2)) {
	  if (nor < 5.*ERROR) {
	    if (myrank == 0)
	      PGFEM_printf ("I will take it\n");
	    nor = 0.0;
	  }
	  else {
	    if (myrank == 0 )
	      PGFEM_printf ("Error in the iteration : iter > iter_max\n"); 
	    INFO = 1;
	    if (gam == 0 && ARC == 1)
	      ART = 1;
	    goto rest;
	  }
	}
      } /* end iter > iter_max */
      iter++;
      BS_nor = 0.0;
    }/* while nor > nor_min */

    /* before increment after convergence, check max damage */
    alpha = max_damage/max_damage_per_step;
    MPI_Allreduce(MPI_IN_PLACE,&alpha,1,MPI_DOUBLE,MPI_MAX,mpi_comm);
    if(alpha > 1.0){
      if(myrank == 0){
	PGFEM_printf("Damage value (%f) is greater than max. damage/step (%f).\n"
	       "Subdividing to maintain accuracy of damage law.\n",
	       alpha*max_damage_per_step,max_damage_per_step);
      }
      INFO = 1;
      goto rest;
    }

    ST = GAMA = gam = 0;
    if (periodic == 1)
      ART = 1;
    else
      ART = 0;

    /* macroscopic load */
    if (periodic == 1)
      macroscopic_load_AL (eps[0].type,lm+dlm,eps);
    
    /* increment coheisve elements */
    if(opts->cohesive){
      increment_cohesive_elements(nce,coel,pores,node,sup,d_r,mp_id);
    }

    /* Finite deformations increment */
    switch(opts->analysis_type){
    case FS_CRPL:
    case FINITE_STRAIN:
      fd_increment (ne,nn,ndofn,npres,matgeom,hommat,elem,node,sup,
		    eps,sig_e,d_r,r,nor_min,crpl,dt,nce,coel,
		    pores,mpi_comm,VVolume,opts, mp_id);
      break;
    case STABILIZED:
      st_increment (ne,nn,ndofn,ndofd,matgeom,hommat,elem,node,sup,
		    eps,sig_e,d_r,r,nor_min,stab,dt,nce,coel,
		    pores,mpi_comm,opts->cohesive,mp_id);
      break;
    case MINI:
      MINI_increment(elem,ne,node,nn,ndofn,sup,eps,
		     sig_e,hommat,d_r,mpi_comm,mp_id);
      break;
    case MINI_3F:
      MINI_3f_increment(elem,ne,node,nn,ndofn,sup,eps,
			sig_e,hommat,d_r,mpi_comm,mp_id);
      break;
    case DISP:
      DISP_increment(elem,ne,node,nn,ndofn,sup,eps,
		     sig_e,hommat,d_r,r,mpi_comm,mp_id);
      break;
    default:
      break;
    }

    /* Add deformation increment into displacement vector */
    vvplus (r,d_r,ndofd);
    
    /* Null prescribed increment deformation */
    for (i=0;i<sup->npd;i++)
      sup->defl_d[i] = 0.0;

    for (i=0;i<ndofd;i++) {
      dR[i] = d_r[i];
      d_r[i] = D_R[i] = rr[i] = f_defl[i] = f[i] = U[i] = 0.0;
    }
    
    /* Tranform dR : L -> G */
    LToG(dR,BS_dR,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
    
    /* Prevent circling */
    if (tim == 0)
      nor = 1.0;
    else
      nor = (sqrt((dlm - fabs(*dlmdlm))*(dlm - fabs(*dlmdlm)))
	     /sqrt(*dlmdlm**dlmdlm));

    if (nor < 1.e-2 && *dlmdlm < 0.0 && dlm > 0.0)
      *AT = 1;
    else
      *AT = 0;
    
    /* Load multiplicator and Determinant */
    DDLM += dlm;
    lm += dlm;
    *dlmdlm = dlm; 
    ddlm += fabs (dlm);
    *DLM0 = dlm0;
    *DET0 = DET;
    *ITT = iter;
    *DAL = dAL;
    
    if (myrank == 0)
      PGFEM_printf ("Inelastic step [%ld] :: STEP = %ld :: NS =  %ld"
	      " || lm = %12.12f : DDLM = %12.12f : dlm = %12.12f\n\n",
	      tim,DIV,STEP,lm,DDLM,dlm); 

    dlm = 0.0;
    if (periodic == 1 && myrank == 0) {
      PGFEM_printf("\nThe total F\n");
      for (i=0;i<3;i++){
	for (j=0;j<3;j++){
	  PGFEM_printf("%12.12f  ",eps[0].F[i][j]);
	}
	PGFEM_printf("\n");
      }
      PGFEM_printf("\n");
    }
    if (*AT == 1 && myrank == 0)
      PGFEM_printf ("*** Trying to prevent cycling ***\n\n");
    
    /************* TEST THE UPDATE FROM N TO N+1  *************/
    if(PFEM_DEBUG || PFEM_DEBUG_ALL || ARC_UPDATE){
      for (i=0;i<ndofd;i++) {f_u[i] = 0.0; d_r[i] = 0.0;}
      fd_residuals (f_u,ne,n_be,ndofn,npres,d_r,r,node,elem,
		    b_elems,matgeom,hommat,sup,eps,sig_e,
		    nor_min,crpl,dts,t,stab,nce,coel,mpi_comm,opts,alpha_alpha,r_n,r_n_1,
		    mp_id);
      for (i=0;i<ndofd;i++) f[i] = lm*R[i] - f_u[i];
      
      LToG(f,BS_f,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
      nor = ss(BS_f,BS_f,DomDof[myrank]);
      MPI_Allreduce(&nor,&tmp,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
      nor = sqrt (tmp);
      
      if (myrank == 0) PGFEM_printf("NORM NORM = %12.12e || %12.12f\n",nor,lm);
    } /* End debug */
    /************* TEST THE UPDATE FROM N TO N+1  *************/
    
    /* Calculating equvivalent Mises stresses and strains vectors */
    Mises (ne,sig_e,eps,opts->analysis_type);
 
    /********************************************************************/
    /***  PRINT INTERMIDIATE OUTPUT FILES FOR CONSTANT LOAD LEVEL     ***/
    /********************************************************************/
    
    if (ddlm > pdt*dAL0/dt0 && STEP > (DIV+1)){
      
      if (print[tim] == 1 && gr2 == VIS_ELIXIR) {
	
	sprintf (jmeno,"%s_%d.out%ld_%ld",out_dat,myrank,tim,FI);
	
	if ((out = fopen(jmeno,"w")) == NULL ){
	  PGFEM_printf("Output file is not possible to open\n");
	  PGFEM_printf("Check the output file and run program again\n");
	  return (0);
	} 
	
	/*logo (out);*/
	
	PGFEM_fprintf (out,"\n");
	PGFEM_fprintf (out,"FINITE DEFORMA. + CRYSTAL PLASTICITY ANALYSIS  : Step - %ld, Time - %f\n",tim,times[tim+1]);
	PGFEM_fprintf (out,"Number of nodes                                : %ld\n",nn);
	PGFEM_fprintf (out,"Number of elements                             : %ld\n",ne);
	PGFEM_fprintf (out,"Number of equations                            : %ld\n",ndofd);
	PGFEM_fprintf (out,"Number of elements in the matrix - SPARSE      : %d\n",Ap[ndofd]);
	PGFEM_fprintf (out,"Load multiplier level                          : %f\n",lm);
	PGFEM_fprintf (out,"\n");
	getrusage (RUSAGE_SELF,&usage);
	PGFEM_fprintf (out,"Time of solution of the system of equations  -  System %ld.%ld, User %ld.%ld\n",
		 usage.ru_stime.tv_sec,usage.ru_stime.tv_usec,
		 usage.ru_utime.tv_sec,usage.ru_utime.tv_usec);

	/* Print MACRO fileds */
	if (periodic == 1) macro_fields_out (out,eps,opts);
	
	/* Print deformations to output file */
	deform (out,node,elem,nn,ne,ndofd,sup,r);
	/* Print stress to output file */
	stress_out (out,ne,nn,elem,sig_e,sig_n,gr4);
	/* Print strain to output file */
	strain_out (out,ne,elem,eps,opts);
	/* Print Fe */
	deform_grad_out (out,ne,elem,eps);
	/* Print cohesive elements */
	if (opts->cohesive == 1) cohesive_out (out,nce,coel);
	
	fclose (out); 
	
	sprintf (jmeno,"%s_.elx%ld_%ld",out_dat,tim,FI);
	FI++;
	
	/* Print to elix file */
	elixir (jmeno,nn,ne,ndofn,node,elem,sup,r,sig_e,
		sig_n,eps,gr4,nce,coel/*,nge,geel,ngn,gnod*/,opts);
      }
      ddlm = 0.0;
    }/* end ddlm > dAL0 */
    
    DIV++;
    if (STEP > 2 && DIV == 2)
      goto rest;
  }/*end SUBDIVISION */

  return (DDLM);
}

/// Perform arc length analysis.
/// It calls the lagacy Arc_length function which requires many function arguments.
/// Later, this function needs to be combined with Arc_length function. 
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv object for field variables
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] time_steps object for time stepping
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in, out] arc an object for Arc length scheme, cotains variables related to Arc length
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] VVolume original volume of the domain
/// \param[in] opts structure PGFem3D option
/// \return load multiplier
double Arc_length_test(GRID *grid,
                       MATERIAL_PROPERTY *mat,
                       FIELD_VARIABLES *fv,
                       SOLVER_OPTIONS *sol,
                       LOADING_STEPS *load,
                       COMMUNICATION_STRUCTURE *com,
                       PGFem3D_TIME_STEPPING *time_steps, 
                       CRPL *crpl,
                       MPI_Comm mpi_comm,
                       const double VVolume,
                       const PGFem3D_opt *opts,
                       const int mp_id)
{
  ARC_LENGTH_VARIABLES *arc = sol->arc;
  double dALMAX = (time_steps->dt_np1)/(arc->dt0)*(arc->dALMAX);
  char out_dat[500];
  sprintf(out_dat,"%s/%s",opts->opath,opts->ofname);
  return Arc_length(grid->ne,
                    grid->n_be,
                    grid->nn,
                    fv->ndofn,
                    fv->ndofd,
                    fv->npres,
                    time_steps->nt,
                    time_steps->tim,
                    time_steps->times,
                    sol->nor_min,
                    sol->iter_max,
                    time_steps->dt_np1,
                    arc->dt0,
                    grid->element,
                    grid->b_elems,
                    com->nbndel,
                    com->bndel,
                    grid->node,
                    load->sups[mp_id],
                    load->sup_defl[mp_id],
                    mat->hommat,
                    mat->matgeom,
                    fv->sig,
                    fv->eps,
                    com->Ap,
                    com->Ai,
                    sol->PGFEM_hypre,
                    fv->RRn,
                    fv->f_defl,
                    crpl,
                    opts->stab,
                    grid->nce,
                    grid->coel,
                    fv->u_np1,
                    fv->f,
                    fv->d_u,
                    arc->D_R,
                    fv->dd_u,
                    fv->R,
                    fv->RR,
                    fv->f_u,
                    arc->U,
                    arc->DK,
                    arc->dR,
                    fv->BS_f,
                    arc->BS_d_r,
                    arc->BS_D_R,
                    arc->BS_rr,
                    arc->BS_R,
                    fv->BS_RR,
                    fv->BS_f_u,
                    arc->BS_U,
                    arc->BS_DK,
                    arc->BS_dR,
                    sol->FNR,
                    arc->lm,
                    arc->dAL0,
                    &(arc->DET0),
                    &(arc->DLM0),
                    &(arc->DLM),
                    opts->vis_format,
                    opts->smoothing,
                    fv->sig_n,
                    out_dat,
                    time_steps->print,
                    &(arc->AT),
                    arc->ARC,
                    dALMAX,
                    &(arc->ITT),
                    &(arc->DAL),
                    &(fv->pores),
                    com->DomDof,
                    com->GDof,
                    com->comm,
                    sol->err,
                    &(fv->NORM),
                    mpi_comm,
                    VVolume,
                    opts,
                    mp_id);
}
