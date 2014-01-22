#include "LINE.h"
#include <math.h>
#include <string.h>

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef INTEGRATION_H
#include "integration.h"
#endif

#ifndef VOL_DAMAGE_INT_ALG_H
#include "vol_damage_int_alg.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef FD_RESIDUALS_H
#include "fd_residuals.h"
#endif

#ifndef MATICE_H
#include "matice.h"
#endif

#ifndef ALM_H
#include "ALM.h"
#endif

#ifndef MACROSCOPIC_LOAD_AL
#include "macroscopic_load_AL.h"
#endif

#ifndef MINI_ELEMENT_H
#include "MINI_element.h"
#endif

#ifndef MINI_3F_ELEMENT_H
#include "MINI_3f_element.h"
#endif

#ifndef DISP_BASED_ELEM_H
#include "bounding_element_utils.h"
#endif

static const int periodic = 0;

long LINE_S1 (double *nor,
	      double *gama,
	      double nor1,
	      double NOR,
	      long iter,
	      double *f_u,
	      long ne,
	      int n_be,
	      long ndofd,
	      long ndofn,
	      long npres,
	      double *d_r,
	      double *r,
	      NODE *node,
	      ELEMENT *elem,
	      BOUNDING_ELEMENT *b_elems,
	      MATGEOM matgeom,
	      HOMMAT *hommat,
	      SUPP sup,
	      EPS *eps,
	      SIG *sig_e,
	      double nor_min,
	      CRPL *crpl,
	      double dt,
	      double stab,
	      long nce,
	      COEL *coel,
	      double *f,
	      double *rr,
	      double *RR,
	      long tim,
	      /*GNOD *gnod,
		GEEL *geel,
	      */double *BS_f,
	      double *BS_RR,
	      double *BS_f_u,
	      long *DomDof,
	      long STEP,
	      COMMUN comm,
	      int GDof,
	      MPI_Comm mpi_comm,
	      double *max_damage,
	      double *dissipation,
	      const PGFem3D_opt *opts)
{
  double GNOR;
  long i,j,N,M,INFO,GInfo;
  char  *error[]={"inf","-inf","nan"},str1[500];

  int myrank,nproc;
  MPI_Comm_rank(mpi_comm,&myrank);
  MPI_Comm_size(mpi_comm,&nproc);
  
  *gama = 1.0; j = 1; INFO = 0;
  while (*nor > NOR && iter != 0){
    if (*nor < 5*nor_min){*nor = 0.0; break;}
    if (*gama == 1.0) *gama = 0.75;
    
    if (myrank == 0) PGFEM_printf ("Gama = %2.3f || nor = %8.8e\n",*gama,*nor);
    
    /* Decrease increment of displacements */
    for (i=0;i<ndofd;i++)  f[i] = *gama*rr[i];
    
    /*************************/
    /* INTEGRATION ALGORITHM */
    /*************************/
    if (opts->analysis_type == FS_CRPL ) {
      INFO = integration_alg (ne,ndofn,ndofd,npres,crpl,elem,node,
			      d_r,f,sup,matgeom,hommat,eps,sig_e,
			      tim,iter,dt,nor_min,STEP,1,opts);
	    
      /* Gather INFO from all domains */
      MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
      if (GInfo == 1) {
	INFO = 1;
	return 1;
      }
    }
    
    /* Update deformations */
    for (i=0;i<ndofd;i++) {
      f[i] += d_r[i];
      f_u[i] = 0.0;
    }
    
    INFO = vol_damage_int_alg(ne,ndofn,f,r,elem,node,
			      hommat,sup,dt,iter,mpi_comm,
			      eps,sig_e,max_damage,dissipation,
			      opts->analysis_type);
      
    MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
    if (GInfo == 1) {
      if(myrank == 0){
	PGFEM_printf("Inverted element detected (vol_damage_int_alg).\n"
	       "Subdividing load.\n");
      }
      INFO = 1;
      return 1;
    }

    /* Residuals */
    fd_residuals (f_u,ne,n_be,ndofn,npres,f,r,node,elem,b_elems,matgeom,
		  hommat,sup,eps,sig_e,nor_min,crpl,dt,stab,
		  nce,coel,mpi_comm,opts);
    
    /* Transform LOCAL load vector to GLOBAL */
    LToG (f_u,BS_f_u,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
    
    /* Compute Euclidian norm */
    for (i=0;i<DomDof[myrank];i++)
      BS_f[i] = BS_RR[i] - BS_f_u[i];

    *nor  = ss (BS_f,BS_f,DomDof[myrank]);
    
    /* Gather residual nor from Domains */
    MPI_Allreduce (nor,&GNOR,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    *nor = sqrt(GNOR);
    *nor /= nor1;
    
    sprintf (str1,"%f",*nor);
    for (N=0;N<3;N++){
      M = 10;
      M = strcmp(error[N],str1);
      
      if (M == 0){
	if (myrank == 0)
	  PGFEM_printf("ERROR in the algorithm : nor = %s\n",error[N]);
	return (1);
      }
    }
    if (*nor < NOR) break;
    
    *gama = 1./exp(j*1./2.);
    j++;
    if (*gama < 0.05){
      if (myrank == 0)
	PGFEM_printf ("Error in the iteration : GAMA = 0.0\n");
      return (1);
    }
  }/* end nor > NOR || GAMA */
  return (INFO);
}

long LINE_S3 (double *nor,
	      double *nor2,
	      double *gama,
	      double nor1,
	      double NOR,
	      double LS1,
	      long iter,
	      double *f_u,
	      long ne,
	      int n_be,
	      long ndofd,
	      long ndofn,
	      long npres,
	      double *d_r,
	      double *r,
	      NODE *node,
	      ELEMENT *elem,
	      BOUNDING_ELEMENT *b_elems,
	      MATGEOM matgeom,
	      HOMMAT *hommat,
	      SUPP sup,
	      EPS *eps,
	      SIG *sig_e,
	      double nor_min,
	      CRPL *crpl,
	      double dt,
	      double stab,
	      long nce,
	      COEL *coel,
	      double *f,
	      double *rr,
	      double *RR,
	      long tim,
	      double *BS_f,
	      double *BS_RR,
	      double *BS_f_u,
	      long *DomDof,
	      COMMUN comm,
	      int GDof,
	      long STEP,
	      MPI_Comm mpi_comm,
	      double *max_damage,
		double *dissipation,
		const PGFem3D_opt *opts)
{
  long i,j,N,M,INFO,GInfo;
  double LS2,slope,tmplam,rhs1,rhs2,AL,a,b,f2,disc,scale,nor3;
  char  *error[]={"inf","-inf","nan"},str1[500];
  
  int myrank,nproc;
  MPI_Comm_rank(mpi_comm,&myrank);
  MPI_Comm_size(mpi_comm,&nproc);

  scale = LS1; LS1 /= scale;
  
  slope = -2.*LS1; 
  
  /* Compute norm LS2 = 1./2. * ss (f,f,ndofd)/scale; */
  nor3 = ss (BS_f,BS_f,DomDof[myrank]);
  /* Gather residual nor from Domains */
  MPI_Allreduce (&nor3,&LS2,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  LS2 *= 1./(2.*scale);

  *gama = 1.0;
  j = 0;
  INFO = 0;
  if (iter == 0)
    return (0);

  if (myrank == 0) 
    PGFEM_printf("nor = %12.8e :: slope = %12.8e || LS1 = %12.8e"
	   "|| LS2 = %12.8e || HU = %12.8e\n",
	   *nor,slope,LS1,LS2,LS1 + 1.e-4**gama*slope); 
  
  while (LS2 > LS1 + 1.e-4**gama*slope){
    
    if (*gama == 1.0){
      tmplam = -slope/(2.*(LS2-LS1-slope));
    }
    else{
      rhs1 = LS2 - LS1 - *gama*slope;
      rhs2 = f2 - LS1 - AL*slope;
      a = (rhs1/(*gama**gama)-rhs2/(AL*AL))/(*gama-AL);
      b = (-AL*rhs1/(*gama**gama)+*gama*rhs2/(AL*AL))/(*gama-AL);
      if (a == 0.0) tmplam = -slope/(2.*b);
      else{
	disc = b*b - 3.*a*slope;
	if (disc < 0.0) tmplam = 0.5**gama;
	else
	  if (b <= 0.0) tmplam = (-b+sqrt(disc))/(3.0*a);
	  else tmplam = -slope/(b+sqrt(disc));
      }/* else a == 0 */
      if (tmplam > 0.5**gama) tmplam = 0.5**gama;
    }/* end *gama == 1*/
    
    /* Gama */
    AL = *gama;
    if (tmplam >= 0.1**gama)
      *gama = tmplam;
    else
      *gama *= 0.1;

    if (*gama < 1.e-1*nor_min)
      return (1);
    
    /* Decrease increment of displacements */
    for (i=0;i<ndofd;i++)
      f[i] = *gama*rr[i];
    
    /*************************/
    /* INTEGRATION ALGORITHM */
    /*************************/
    if (opts->analysis_type == FS_CRPL) {
      INFO = integration_alg (ne,ndofn,ndofd,npres,crpl,elem,
			      node,d_r,f,sup,matgeom,hommat,eps,
			      sig_e,tim,iter,dt,nor_min,STEP,1,opts);
    
      /* Gather INFO from all domains */
      MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
      if (GInfo == 1) {
	INFO = 1;
	return (1);
      }
    }
    
    /* Update deformations */
    for (i=0;i<ndofd;i++) {
      f[i] += d_r[i];
      f_u[i] = 0.0;
    }
    
    INFO = vol_damage_int_alg(ne,ndofn,f,r,elem,node,
			      hommat,sup,dt,iter,mpi_comm,
			      eps,sig_e,max_damage,dissipation,
			      opts->analysis_type);

    bounding_element_communicate_damage(n_be,b_elems,ne,eps,mpi_comm);

    MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
    if (GInfo == 1) {
      if(myrank == 0){
	PGFEM_printf("Inverted element detected (vol_damage_int_alg).\n"
	       "Subdividing load.\n");
      }
      INFO = 1;
      return 1;
    }

    /* Residuals */
    fd_residuals (f_u,ne,n_be,ndofn,npres,f,r,node,elem,b_elems,matgeom,
		  hommat,sup,eps,sig_e,nor_min,crpl,dt,stab,
		  nce,coel/*,gnod,geel*/,mpi_comm,opts);
	
    /* Compute Euclidian norm */
    for (i=0;i<ndofd;i++)
      f[i] = RR[i] - f_u[i];

    /* Local to global transformation */
    LToG(f,BS_f,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
    
    nor3 = ss (BS_f,BS_f,DomDof[myrank]);
    
    /* Gather residual nor from Domains */
    MPI_Allreduce(&nor3,nor,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    nor3 = *nor;
    *nor = *nor2 = sqrt(nor3);
    *nor /= nor1;
    
    sprintf (str1,"%f",*nor);
    for (N=0;N<3;N++){
      M = 10;
      M = strcmp(error[N],str1);

      if (M == 0){
	if (myrank == 0)
	  PGFEM_printf("ERROR in the algorithm LS2 : nor = %s\n",error[N]);
	return (1);
      }
    }
    
    /* New LS and gama */
    f2 = LS2;
    LS2 = 1./2. * nor3/scale;
    
    if (myrank ==0)
      PGFEM_printf ("gama = %12.12f || tmplam = %12.8f ||"
	      "NOR = %12.8e  || LS2 = %12.8e\n",
	      *gama,tmplam,*nor,LS2);
    
  }/* end LINE SEARCH */
  
  return (INFO);
}

long ALINE_S3 (long ARC,
	       double *DLM,
	       double *nor,
	       double *nor2,
	       double *gama,
	       double nor1,
	       double LS1,
	       long iter,
	       long ne,
	       int n_be,
	       long ndofd,
	       long ndofn,
	       long npres,
	       long tim,
	       double nor_min,
	       double dt,
	       double stab,
	       long nce,
	       double dlm,
	       double lm,
	       double dAL,
	       double *d_r,
	       double *r,
	       double *D_R,
	       NODE *node,
	       ELEMENT *elem,
	       BOUNDING_ELEMENT *b_elems,
	       MATGEOM matgeom,
	       HOMMAT *hommat,
	       SUPP sup,
	       EPS *eps,
	       SIG *sig_e,
	       CRPL *crpl,
	       COEL *coel,
	       double *f_u,
	       double *f,
	       double *R/*,
			  GNOD *gnod,
			  GEEL *geel*/,
 	       double *BS_f,
	       double *BS_R,
	       double *BS_D_R,
	       double *BS_d_r,
	       double *BS_DK,
	       double *BS_U,
	       double *BS_rr,
	       long *DomDof,
 	       int GDof,
	       COMMUN comm,
	       long STEP,
	       MPI_Comm mpi_comm,
	      double *max_damage,
		double *dissipation,
		const PGFem3D_opt *opts )
{
  long i,j,N,M,INFO,GInfo;
  double LS2,slope,tmplam,rhs1,rhs2,AL,a,b,f2,disc,scale,tmp;
  char  *error[]={"inf","-inf","nan"},str1[500];

  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);

  scale = LS1;
  LS1 /= scale;
  slope = -2.*LS1; 
  tmp = ss(BS_f,BS_f,DomDof[myrank]);
  MPI_Allreduce(&tmp,&LS2,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
  LS2 *= 1./(2.*scale);

  *gama = 1.0;
  j = 0;
  INFO = 0;
  if (iter == 0) return (0);

  if (myrank == 0)
    PGFEM_printf("nor = %12.8e :: slope = %12.8e || LS1 = %12.8e "
	   "|| LS2 = %12.8e || HU = %12.8e\n",
	   *nor,slope,LS1,LS2,LS1 + 1.e-4**gama*slope);
  
  while (LS2 > LS1 + 1.e-4**gama*slope){
    
    if (*gama == 1.0){
      tmplam = -slope/(2.*(LS2-LS1-slope));
    }
    else{
      rhs1 = LS2 - LS1 - *gama*slope;
      rhs2 = f2 - LS1 - AL*slope;
      a = (rhs1/(*gama**gama)-rhs2/(AL*AL))/(*gama-AL);
      b = (-AL*rhs1/(*gama**gama)+*gama*rhs2/(AL*AL))/(*gama-AL);
      if (a == 0.0) tmplam = -slope/(2.*b);
      else{
	disc = b*b - 3.*a*slope;
	if (disc < 0.0) tmplam = 0.5**gama;
	else
	  if (b <= 0.0) tmplam = (-b+sqrt(disc))/(3.0*a);
	  else tmplam = -slope/(b+sqrt(disc));
      }/* else a == 0 */
      if (tmplam > 0.5**gama) tmplam = 0.5**gama;
    }/* end *gama == 1*/
        
    /* Gama */
    AL = *gama;
    if (tmplam >= 0.1**gama)
      *gama = tmplam;
    else
      *gama *= 0.1;
    if (*gama < 1.e-1*nor_min)
      return (1);
    
    /* Decrease increment of displacements */
    for (i=0;i<DomDof[myrank];i++) {
      BS_D_R[i] = *gama*BS_U[i];
      BS_f[i] = *gama*BS_rr[i];
    }
    
    /* dlam */
    if (ARC == 0){
      *DLM = D_lam_ALM (ndofd,BS_f,BS_d_r,BS_D_R,BS_R,BS_DK,
			dlm,dAL,DomDof,mpi_comm); 
      /* *DLM = D_lam_ALM2 (f,D_R,R,DK,dlm,lm,dAL,ne,ndofd,
	                    ndofn,npres,d_r,r,node,elem,matgeom,
			    hommat,sup,eps,sig_e,nor_min,crpl,dt,
			    stab,nce,coel); */
    }
    if (ARC == 1)
      *DLM = D_lam_ALM4 (ndofd,BS_f,BS_d_r,BS_D_R,BS_DK,dlm,
			 dAL,DomDof,mpi_comm);
    
    sprintf (str1,"%f",*DLM);
    for (N=0;N<3;N++){
      M = 10;
      M = strcmp(error[N],str1);
      
      if (M == 0) {
	if (myrank == 0) PGFEM_printf("Complex root in ARC-LENGTH method\n");
	return (1);
      }
    }
    
    /* Displacement update */
    for (i=0;i<DomDof[myrank];i++)
      BS_D_R[i] = *gama * (BS_U[i] + *DLM*BS_rr[i]);
    
    /* macroscopic load */
    if (periodic == 1) macroscopic_load_AL (eps[0].type,lm+dlm+*DLM,eps);

    /* Global to local transformation */
    GToL(BS_D_R,D_R,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);

    /*************************/
    /* INTEGRATION ALGORITHM */
    /*************************/
    if (opts->analysis_type == FS_CRPL) {
      INFO=integration_alg (ne,ndofn,ndofd,npres,crpl,elem,node,
			    d_r,D_R,sup,matgeom,hommat,eps,sig_e,
			    tim,iter,dt,nor_min,STEP,0,opts); 
      
      /* Gather INFO from all domains */
      MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
      if (GInfo == 1) {
	INFO = 1;
	return (1);
      }
    }
    
    /* Update deformations */
    for (i=0;i<ndofd;i++) {
      f[i] = d_r[i] + D_R[i];
      f_u[i] = 0.0;
    } 

    INFO = vol_damage_int_alg(ne,ndofn,f,r,elem,node,
			      hommat,sup,dt,iter,mpi_comm,
			      eps,sig_e,max_damage,dissipation,
			      opts->analysis_type);
      
    MPI_Allreduce (&INFO,&GInfo,1,MPI_LONG,MPI_BOR,mpi_comm);
    if (GInfo == 1) {
      if(myrank == 0){
	PGFEM_printf("Inverted element detected (vol_damage_int_alg).\n"
	       "Subdividing load.\n");
      }
      INFO = 1;
      return 1;
    }
    
    /* Residuals */
    fd_residuals (f_u,ne,n_be,ndofn,npres,f,r,node,elem,b_elems,matgeom,
		  hommat,sup,eps,sig_e,nor_min,crpl,dt,stab,
		  nce,coel/*,gnod,geel*/,mpi_comm,opts );
    
    /* Compute Euclidean norm */
    for (i=0;i<ndofd;i++)
      f[i] = (lm + dlm + *DLM)*R[i] - f_u[i];

    /* L -> G : f */
    LToG(f,BS_f,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);

    *nor = ss(BS_f,BS_f,DomDof[myrank]);
    MPI_Allreduce(nor,&tmp,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    
    *nor = *nor2 = sqrt (tmp);
    *nor /= nor1;
    sprintf (str1,"%f",*nor);
    for (N=0;N<3;N++){
      M = 10;
      M = strcmp(error[N],str1);
      if (M == 0){
      if (myrank == 0)
	PGFEM_printf("ERROR in the algorithm LS2 : nor = %s\n",error[N]);
      return (1);
      }
    }
    
    /* New LS and gama */
    f2 = LS2;
    LS2 = 1./2. * tmp/scale; 
    
    if (myrank == 0)
      PGFEM_printf ("gama = %12.12f || tmplam = %12.8f || NOR = %12.8e "
	      " || LS2 = %12.8e\n",
	      *gama,tmplam,*nor,LS2);
    
  }/* end LINE SEARCH */
  
  return (INFO);
}
