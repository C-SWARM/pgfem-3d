#include "subdivision.h"
#include <math.h>

#include "PGFEM_io.h"
#include "enumerations.h"
#include "matice.h"
#include "stabilized.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "res_fini_def.h"
#include "elem3d.h"
#include "constitutive_model.h"

static const int MAX_STEP = 10000;
static const double MIN_D_TIME = 1.0e-10;
static const int periodic = 0;

void subdivision (long INFO,
		  double *dt,
		  long *STEP,
		  long *DIV,
		  long tim,
		  double *times,
		  long *ST,
		  long ne,
		  long ndofn,
		  long ndofd,
		  long npres,
		  ELEMENT *elem,
		  CRPL *crpl,
		  EPS *eps,
		  SIG *sig,
		  SUPP sup,
		  double *sup_defl,
		  double *rr,
		  double *d_r,
		  double *f_defl,
		  double *f,
		  double *RRn,
		  double *R,
		  long *GAMA,
		  double *DT,
		  long *OME,
		  double stab,
		  int iter,
		  int max_iter,
		  double alpha,
		  MPI_Comm mpi_comm,
		  const int analysis)
{
  long i;
  double nor,gama;
  
  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);

  if (INFO == 0) { /* Step converged */
    *GAMA = 0; /* reset decelleration flag */
    if (*DIV == 0){ /* 1st subdivision */
      if (tim == 0) *STEP = 1;
      else {/* not 1st step, compute subdivision based on previous
	       step */

	/* compute STEP based on previous convergence properties
	   (times[tim-1] is modified by subdivision procedure) */
	*STEP +=((times[tim+1] - times[tim])
		 /(times[tim] - times[tim-1]) - 1);

	if (*STEP <= 0) *STEP = 1;
	*dt = (times[tim+1] - times[tim])/(*STEP);
      }
    }/* end DIV == 0 */
    else { /* *DIV != 0 */
      /*===== ACCELERATE =====*/
      /* compute the acceleration scaling parameter. NOTE: smaller
	 gama --> larger acceleration */
      if(alpha <= 0.5){ /* physics allows large acceleration */
	if(iter > max_iter){/* max_iter reached, accel. slowly*/
	  *DT = *dt;
	  gama = 0.91; /* slow acceleration */
	  *OME = 0; /* do not over accelerate if catch next step */
	} else if (*OME == 0) { /* first acceleration */
	  *DT = *dt;
	  gama = 0.75;
	  *OME = 1;
	} else {
	  gama = 1./exp(*OME*1./2.);
	  *OME += 1;
	}
      } else if(alpha <= 0.8){ /* physics nearing optimum */
	*DT = *dt;
	gama = 0.8;
	(*OME) = 0;
      } else if(alpha > 0.8){ /* physics very near optimum,
				 asymptotically accelerate */
	*DT = *dt;
	gama = alpha; /* note alpha > 1 automatically decellerates! */
	(*OME) = 0;
      }
  
      *dt = (times[tim+1] - times[tim])/(*STEP)**DIV;
      times[tim] += *dt;
      
      /* compute new load increment */
      if (periodic != 1){
	for (i=0;i<sup->npd;i++) {
	  nor = sup_defl[i]/(*STEP)**DIV;
	  sup_defl[i] -= nor;
	}
	for (i=0;i<ndofd;i++)  {
	  RRn[i] += R[i]/(*STEP)**DIV;
	  nor = R[i]/(*STEP)**DIV;
	  R[i] -= nor;
	}
      }
      
      /* compute new number of steps */
      *dt = *DT/gama;
      nor = (times[tim+1] - times[tim])/(*dt);
      *STEP = round1 (nor);
      if (*STEP <= 0)
	*STEP = 1;
      *dt = (times[tim+1] - times[tim])/(*STEP);
      if (*STEP == 1) 
	*ST = 1;
      *DIV = 0;
    }
  }/* end INFO = 0 i.e., step/sub-step DID NOT converge */
  else{
    *OME = 0; /* reset acceleration flag */
    if (*DIV == 0){/* 1st subdivision */
      if (*GAMA == 0) { /* 1st load stepping */
	if(alpha > 4./3.){/* decellerate based on physics. If would
			     decellerate more based on physics, then
			     do it, otherwise use exponential */
	  *DT = *dt;
	  gama = 1.0/alpha;
	  *GAMA = 0;
	} else {
	  *DT = *dt;
	  gama = 0.75;
	  *GAMA = 1;
	}
      } else {
	gama = 1./exp(*GAMA*1./2.);
	*GAMA += 1;
      }
 
      /* compute number of steps */     
      *dt = *DT*gama;
      nor = (times[tim+1] - times[tim])/(*dt);
      i = round1 (nor);
      if (i <= *STEP)
	*STEP += 1;
      else
	*STEP = i;
      if (*STEP <= 0)
	*STEP = 1;
      *dt = (times[tim+1] - times[tim])/(*STEP);
      *OME = 0;
    }/* end DIV = 0 */
    else{ /* load is subdivided and successfully completed at least
	     one step. Update time, compute remaining load and
	     subdivide again */
      *DT = *dt;
      gama = 0.75;
      *GAMA = 1;
      *dt = (times[tim+1] - times[tim])/(*STEP)**DIV;
      times[tim] += *dt;
      
      /* compute new load increment */
      if (periodic != 1){
	for (i=0;i<sup->npd;i++) {
	  nor = sup_defl[i]/(*STEP)**DIV;
	  sup_defl[i] -= nor;
	}
	for (i=0;i<ndofd;i++)  {
	  RRn[i] += R[i]/(*STEP)**DIV;
	  nor = R[i]/(*STEP)**DIV;
	  R[i] -= nor;
	}
      }
      
      /* compute number of steps */
      *dt = *DT*gama;
      nor = (times[tim+1] - times[tim])/(*dt);
      i = round1 (nor);
      if (i <= (*STEP-*DIV))
	*STEP += 1;
      if (*STEP <= 0)
	*STEP = 1;
      *dt = (times[tim+1] - times[tim])/(*STEP);
      *DIV = 0;
      *OME = 0;
    }/* end DIV != 0 */
    
    if (myrank == 0){
      PGFEM_printf ("\n[%ld] Sub. steps = %ld :: gama = %f ||"
	      " Time %f | dt = %f\n",
	      tim,*STEP,gama,times[tim+1],*dt);
    }
    
    /* Reset variables */
    switch(analysis){
    case FS_CRPL:
    case FINITE_STRAIN:
      res_fini_def (ne,npres,elem,eps,sig,crpl,analysis);
      break;
    case STABILIZED:
      res_stab_def (ne,npres,elem,eps,sig,stab);
      break;
    case MINI:
      MINI_reset(elem,ne,npres,sig);    
      break;
    case MINI_3F:
      MINI_3f_reset(elem,ne,npres,4,sig,eps); 
      break;
    case CM:
      constitutive_model_reset_state(eps, ne, elem);
      break;
    default: break;
    }

    /* reset damage variables */
    for(i=0; i<ne; i++){
      long n_ip;
      int_point(elem[i].toe,&n_ip);
      for(int j=0; j<n_ip; j++){
	reset_damage(&eps[i].dam[j]);
      }
      if(analysis == STABILIZED){/* get pressure terms too */
	int_point(10,&n_ip);
	for(int j=0; j<n_ip; j++){
	  reset_damage(&eps[i].st[j].dam);
	}
      }
    }

    for (i=0;i<ndofd;i++) {rr[i] = d_r[i] = f_defl[i] = f[i] = 0.0;}
    
    if (gama < 1.0 / MAX_STEP || *STEP > MAX_STEP || *dt  < MIN_D_TIME) {
      if (myrank == 0)
	PGFEM_printf ("Error in Subdivision routine\n");
      PGFEM_Comm_code_abort (mpi_comm,i);
    }
  }/* end INFO != 0 */
}


double subdiv_arc (long INFO,
		   double *dt,
		   double dt0,
		   long *STEP,
		   long *DIV,
		   long tim,
		   double *times,
		   long *ST,
		   long ne,
		   long ndofd,
		   long npres,
		   ELEMENT *elem,
		   CRPL *crpl,
		   EPS *eps,
		   SIG *sig,
		   SUPP sup,
		   double *sup_defl,
		   double *rr,
		   double *d_r,
		   double *D_R,
		   double *f_defl,
		   double *f,
		   long *GAMA,
		   double *DT,
		   long *OME,
		   double stab,
		   double dAL0,
		   double dAL,
		   double dALMAX,
		   double nor_min,
		   double dlm0,
		   long *ITT,
		   long iter,
		   long iter_max,
		   long TYPE,
		   MPI_Comm mpi_comm,
		   const int analysis)
{
  long i;
  double nor,gama,DLM0,I1,I2;

  static int last_it  = 1;

  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);
  
  DLM0 = 0.0; 

  if (tim == 0 && iter == 0 && INFO == 0)
    return ((*dt/dt0)*dAL0);

  if (TYPE == 0){/* Time subdivision scheme */
    if (INFO == 0) {
      if (*DIV == 0){
	nor = (times[tim+1] - times[tim])/(times[tim] - times[tim-1]);
	*STEP =  round1 (nor) - 1;
	if (*STEP <= 0)
	  *STEP = 1; 
	*dt = (times[tim+1] - times[tim])/(*STEP);
	
	/*
	  DDT = (times[tim+1] - times[tim])/(times[tim] - times[tim-1]);
	  I1 = *ITT-1;
	  I2 = (iter_max+2)/2;
	  if (I1 <= 0)
	    I1 = 1;
	  *dt = sqrt (I2/I1)*DDT;
	  nor = (times[tim+1] - times[tim])/(*dt);
	  *STEP = round1 (nor);
	  if (*STEP <= 0)
	    *STEP = 1;
	  *dt = (times[tim+1] - times[tim])/(*STEP);
	  PGFEM_printf ("NS =  %ld || I1 = %ld || I2 = %ld ::"
	          " Time %f | dt = %10.10f\n",*STEP,*ITT-1,
		  (iter_max+2)/2,times[tim+1],*dt);
	*/
      }/* end DIV == 0 */
      else{
	if (*OME == 0) {
	  *DT = *dt;
	  gama = 0.75;
	  *OME = 1;
	}
	else {
	  gama = 1./exp(*OME*1./2.);
	  *OME += 1;
	}
	*dt = (times[tim+1] - times[tim])/(*STEP)**DIV;
	times[tim] += *dt;
	
	*dt = *DT/gama;
	nor = (times[tim+1] - times[tim])/(*dt);
	*STEP = round1 (nor);
	if (*STEP <= 0)
	  *STEP = 1;
	
	/* Restart  */
	if (fabs(dlm0) < nor_min) {
	  if (myrank == 0) PGFEM_printf ("Trying increment length restart\n");
	  *STEP = 1;
	}
	*dt = (times[tim+1] - times[tim])/(*STEP);
	if (*STEP == 1)
	  *ST = 1;
	*DIV = 0;
	
	/* DDT = *dt; */
	/* *dt = (times[tim+1] - times[tim])/(*STEP)**DIV; */
	/* times[tim] += *dt; */
	/* I1 = *ITT-1; */
	/* I2 = (iter_max+2)/2; */
	/* if (I1 <= 0) */
	/*   I1 = 1; */
	/* *dt = sqrt (I2/I1)*DDT; */
	  
	/* nor = (times[tim+1] - times[tim])/(*dt); */
	/* *STEP = round1 (nor); */
	/* if (*STEP <= 0) */
	/*   *STEP = 1; */
	  
	/* /\* Restart *\/ */
	/* if (fabs(dlm0) < nor_min) { */
	/*   PGFEM_printf ("Trying increment length restart\n"); */
	/*   *STEP = 1; */
	/*   } */
	  
	/* *dt = (times[tim+1] - times[tim])/(*STEP); */
	/* if (*STEP == 1) */
	/*   *ST = 1; */
	/* *DIV = 0; */
	/* PGFEM_printf ("STEP = %ld :: NS =  %ld || I1 = %ld || I2 = %ld" */
	/*         " :: Time %f | dt = %10.10f\n",*DIV,*STEP,*ITT-1, */
	/* 	  (iter_max+2)/2,times[tim+1],*dt); */

      }
    }/* end INFO = 0 */
    else{
      if (*DIV == 0){
	if (*GAMA == 0) {
	  *DT = *dt;
	  gama = 0.75;
	  *GAMA = 2;
	} else {
	  gama = 1./exp(*GAMA*1./2.);
	  *GAMA += 2;
	}
	
	*dt = *DT*gama;
	nor = (times[tim+1] - times[tim])/(*dt);
	i = round1 (nor);
	if (i <= *STEP)
	  *STEP += 1;
	else
	  *STEP = i;
	if (*STEP <= 0)
	  *STEP = 1;
	
	/* Limit long overflow */
	if (*STEP >= 200000000) {
	  *STEP = 1;
	  *DIV = 0;
	  *GAMA = 0;
	  gama = 0.0;
	}
	
	*dt = (times[tim+1] - times[tim])/(*STEP);
	*OME = 0;
      }/* end DIV = 0 */
      else{
	*DT = *dt;
	gama = 0.75;
	*GAMA = 2;
	
	*dt = (times[tim+1] - times[tim])/(*STEP)**DIV;
	times[tim] += *dt;
	*dt = *DT*gama;
	nor = (times[tim+1] - times[tim])/(*dt);
	i = round1 (nor);
	if (i <= (*STEP-*DIV))
	  *STEP += 1;
	if (*STEP <= 0)
	  *STEP = 1;
	
	/* Limit long overflow */
	if (*STEP >= 200000000) {
	  *STEP = 1;
	  *DIV = 0;
	  *GAMA = 0;
	  gama = 0.0;
	}
	
	*dt = (times[tim+1] - times[tim])/(*STEP);
	*DIV = 0;
	*OME = 0;
      }/* DIV != 0 */
      
      if (myrank == 0)
	PGFEM_printf ("\n[%ld] Sub. steps = %ld :: gama = %8.8e ||"
		" Time %f | dt = %10.10f\n",tim,*STEP,
		gama,times[tim+1],*dt);
      
      /* Reset variables */
      switch(analysis){
      case FS_CRPL:
      case FINITE_STRAIN:
	res_fini_def (ne,npres,elem,eps,sig,crpl,analysis);
	break;
      case STABILIZED:
	res_stab_def (ne,npres,elem,eps,sig,stab);
	break;
      case MINI:
	MINI_reset(elem,ne,npres,sig);  
	break;
      case MINI_3F:
	MINI_3f_reset(elem,ne,npres,4,sig,eps);
	break;
      default: break;
      }

      /* reset damage variables */
      for(i=0; i<ne; i++){
	long n_ip;
	int_point(elem[i].toe,&n_ip);
	for(int j=0; j<n_ip; j++){
	  reset_damage(&eps[i].dam[j]);
	}
	if(analysis == STABILIZED){/* get pressure terms too */
	  int_point(10,&n_ip);
	  for(int j=0; j<n_ip; j++){
	    reset_damage(&eps[i].st[j].dam);
	  }
	}
      }

      for (i=0;i<ndofd;i++) {
	rr[i] = d_r[i] = D_R[i] = f_defl[i] = f[i] = 0.0;
      }
      
      if (gama < nor_min && (*dt/dt0)*dAL0 < nor_min*nor_min) {
	if (myrank == 0)
	  PGFEM_printf ("Error in Subdivision rutine\n"); 
	PGFEM_Comm_code_abort(mpi_comm,0); }
    }/* end INFO = 1 */
    
    DLM0 = (*dt/dt0)*dAL0;
  }/* end TYPE =  0 */
  /************************************************************/
  /************************************************************/
  if (TYPE == 1){/* Arc-length subdivision procedure */

    I1 = (double) *ITT;
    I2 = (double) iter_max; //round1((iter_max+2)/2.);
        
    if (INFO == 0) {
      if (tim == 0){
	DLM0 = dAL0;
      } else {
	if (I1 <= 0) I1 = 1;
	if ( I1 <= last_it && last_it < iter_max){
	  I1 = (I1<last_it)?I1:(last_it - 1);
	}
	DLM0 = sqrt (I2/I1)*dAL;
	if (DLM0 > dALMAX) DLM0 = dALMAX;
      }
    }/* end INFO = 0 */
    else{
      if (*GAMA == 0) {
	gama = 0.75;
	*GAMA = 2;
      } else {
	gama = 1./exp(*GAMA*1./2.);
	*GAMA += 2;
      }
      
      DLM0 = dAL*gama;
      
      if (myrank == 0)
	PGFEM_printf ("\n[%ld] gama = %8.8e || Time %f | dt = %10.10f\n",
		tim,gama,times[tim+1],*dt);
      
      /* Reset variables */
      switch(analysis){
      case FS_CRPL:
      case FINITE_STRAIN:
	res_fini_def (ne,npres,elem,eps,sig,crpl,analysis);
	break;
      case STABILIZED:
	res_stab_def (ne,npres,elem,eps,sig,stab);
	break;
      case MINI:
	MINI_reset(elem,ne,npres,sig);
	break;
      case MINI_3F:
	MINI_3f_reset(elem,ne,npres,4,sig,eps);
	break;
      default: break;
      }

      /* reset damage variables */
      for(i=0; i<ne; i++){
	long n_ip;
	int_point(elem[i].toe,&n_ip);
	for(int j=0; j<n_ip; j++){
	  reset_damage(&eps[i].dam[j]);
	}
	if(analysis == STABILIZED){/* get pressure terms too */
	  int_point(10,&n_ip);
	  for(int j=0; j<n_ip; j++){
	    reset_damage(&eps[i].st[j].dam);
	  }
	}
      }


      for (i=0;i<ndofd;i++) {
	rr[i] = d_r[i] = D_R[i] = f_defl[i] = f[i] = 0.0;
      }
      
      if (DLM0 < nor_min*nor_min) {
	if (myrank == 0)
	  PGFEM_printf ("Error in Subdivision rutine\n");
	PGFEM_Comm_code_abort(mpi_comm,0);
      }
    }/* end INFO = 1 */

    last_it = *ITT; // Save last number of iterations in static variale

    *STEP = 1; *DIV = 0;
    
    if (myrank == 0)
      PGFEM_printf ("I1 = %f || I2 = %f :: dAL = %12.12f\n",I1,I2,dAL);
  }/* end TYPE = 1 */
  
  return (DLM0);
}
