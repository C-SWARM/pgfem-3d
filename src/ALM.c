#include "ALM.h"
#include <math.h>

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef FD_RESIDUALS_H
#include "fd_residuals.h"
#endif

#ifndef MATICE_H
#include "matice.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef ALM_DEBUG
#define ALM_DEBUG 0
#endif

static const double PII = 3.141592653589793238462643;

double D_lam_ALM (long ndofd,
		  double *BS_rr,
		  double *BS_d_r,
		  double *BS_D_R,
		  double *BS_R,
		  double *BS_DK,
		  double dlm,
		  double dAL,
		  long *DomDof,
		  MPI_Comm mpi_comm)
{
  long i;
  double /* s1,s2, */b,a1,a2,a3,an1,an2,DLM,x1,x2,*p1,*p2,*p3/* ,tmp */;
  int myrank;

  MPI_Comm_rank(mpi_comm,&myrank);

  DLM = 0.0;
  b = 0.0;
  
  p1 = (double *)aloc1 (DomDof[myrank]);
  p2 = (double *)aloc1 (DomDof[myrank]);
  p3 = (double *)aloc1 (DomDof[myrank]);

  for (i=0;i<DomDof[myrank];i++) {
    p1[i] = BS_rr[i]/BS_DK[i];
    p2[i] = (BS_d_r[i] + BS_D_R[i])/BS_DK[i];
    p3[i] = BS_d_r[i];
  }
  
  double send[4], rec[4];
  double p1p1, RR, p1p2, p2p2, p1p3, p2p3, p3p3;

  /*** Can improve performance by using a matrix-matrix product and
       using the appropriate components of the reslting matrix ***/

  /* Compute dot products  and pack */
  send[0] = p1p1 = ss (p1,p1,DomDof[myrank]);
  send[1] = RR	 = ss (BS_R,BS_R,DomDof[myrank]);
  send[2] = p1p2 = ss (p1,p2,DomDof[myrank]);
  send[3] = p2p2 = ss (p2,p2,DomDof[myrank]);

  /* MPI_Allreduce */
  MPI_Allreduce(send,rec,4,MPI_DOUBLE,MPI_SUM,mpi_comm);

  /* unpack recieve container */
  p1p1 = rec[0];
  RR   = rec[1];
  p1p2 = rec[2];
  p2p2 = rec[3];

  /* A1 */ 
  a1 = p1p1 + b*b*RR;
	   
  /* A2 */ 
  a2 = 2.*(b*b*dlm*RR + p1p2);	   

  /* A3 */ 
  a3 = p2p2 - dAL*dAL + b*b*dlm*dlm*RR;

  /* /\******\/ */
  /* /\* A1 *\/ */
  /* /\******\/ */
  /* tmp = ss (p1,p1,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s1,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */
  /* tmp = ss (BS_R,BS_R,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s2,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */

  /* a1 = s1 + b*b*s2; */
  /* /\* PGFEM_printf("s1 = %12.12f || a1 = %12.12f\n",s1,a1); *\/ */
  
  /* /\******\/ */
  /* /\* A2 *\/ */
  /* /\******\/ */
  /* tmp = ss (p1,p2,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s1,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */
  /* a2 = 2.*(b*b*dlm*s2 + s1); */
  /* /\* PGFEM_printf("s2 = %12.12f || a2 = %12.12f\n",s1,a2); *\/ */
  
  /* /\******\/ */
  /* /\* A3 *\/ */
  /* /\******\/ */
  /* tmp = ss (p2,p2,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s1,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */
  /* a3 = s1 - dAL*dAL + b*b*dlm*dlm*s2; */
  /* /\* PGFEM_printf("s3 = %12.12f || a3 = %12.12f\n",s1,a3); *\/ */
  
  if ((a2*a2 - 4.*a1*a3) < 0.0 || 2.*a1 == 0.0){
    if (myrank == 0)
      PGFEM_printf ("(a2*a2 - 4.*a1*a3) = %12.12f || 2*a1 = %12.12f\n",
	      (a2*a2 - 4.*a1*a3),2.*a1); 
    DLM = 1./0.0;
    free(p1);
    free(p2);
    free(p3);
    return (DLM); /* ??? return inf? */
  }
  
  /**** x1,2 = (-b +- sqrt(b^2 - 4*a*c))/(2*a) ****/
  x1 = (-1.*a2 + sqrt (a2*a2 - 4.*a1*a3))/(2.*a1);
  x2 = (-1.*a2 - sqrt (a2*a2 - 4.*a1*a3))/(2.*a1);
  
  p1p1 = p2p2 = p3p3 = 0.0;
  for (i=0;i<DomDof[myrank];i++){
    p1[i] = BS_d_r[i] + BS_D_R[i] + x1*BS_rr[i];
    p2[i] = BS_d_r[i] + BS_D_R[i] + x2*BS_rr[i];

    /* added to reduce # MPI_Allreduce */
    p1p1 += p1[i]*p1[i];
    p2p2 += p2[i]*p2[i];
    p3p3 += p3[i]*p3[i];
  }

  /* compute magnitude of vectors */  
  send[0] = p1p1;
  send[1] = p2p2;
  send[2] = p3p3;
  MPI_Allreduce(send,rec,3,MPI_DOUBLE,MPI_SUM,mpi_comm);
  p1p1 = rec[0];
  p2p2 = rec[1];
  p3p3 = rec[2];

  p1p1 = sqrt(p1p1);
  p2p2 = sqrt(p2p2);
  p3p3 = sqrt(p3p3);

  p1p3 = 0.0;
  p2p3 = 0.0;
  for(i=0; i<DomDof[myrank]; i++){
    p1p3 += p1[i]*p3[i];
    p2p3 += p2[i]*p3[i];
  }
  p1p3 /=(p1p1*p3p3);
  p2p3 /=(p2p2*p3p3);

  send[0] = p1p3;
  send[1] = p2p3;
  MPI_Allreduce(send,rec,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
  p1p3 = rec[0];
  p2p3 = rec[1];

  /* Cos angle */
  an1 = acos (p1p3);
  an2 = acos (p2p3);

  /* nor_vec (p1,DomDof[myrank],mpi_comm); */
  /* nor_vec (p2,DomDof[myrank],mpi_comm); */
  /* nor_vec (p3,DomDof[myrank],mpi_comm); */
  
  /* tmp = ss (p1,p3,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s1,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */
  /* tmp = ss (p2,p3,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s1,1,MPI_DOUBLE,MPI_SUM,mpi_comm); /\* <=== BUG??? s2? *\/ */

  /* /\* Cos angle *\/ */
  /* an1 = acos (s1); */
  /* an2 = acos (s2); */
  
  /* PGFEM_printf ("s1 = %12.12f : s2 = %12.12f || a1 = %5.5f :"
             " a2 = %5.5f\n",s1,s2,an1,an2); */
  
  if (an1 <= an2) DLM = x1;
  else            DLM = x2;
    
  /*
    if (s1 < 0.0 && s2 > 0.0) DLM = x2;
    if (s2 < 0.0 && s1 > 0.0) DLM = x1;
    
    LIN = -1.*a3/a2;
    
    if (s1 > 0.0 && s2 > 0.0){
    
    s1 = sqrt ((x1 - LIN)*(x1 - LIN));
    s2 = sqrt ((x2 - LIN)*(x2 - LIN));
    
    if (s1 <= s2)  DLM = x1;
    else           DLM = x2;
    }
  */
  
  /* if (s1 == 0.0 && s2 == 0.0)  DLM = 0.0;*/
  if (p1p3 == 0.0 && p2p3 == 0.0)  DLM = 0.0;
  
  dealoc1 (p1);
  dealoc1 (p2);
  dealoc1 (p3);
  
  return (DLM);
}

/***********************************************************************/
/***********************************************************************/

double d_ALM2 (long ndofd,
	       double *rr,
	       double *R,
	       double *DK,
	       double d_lm)
{
  double DLM=0;
  return (DLM);
}

double d_lam_ALM2 (long ndofd,
		   double *rr,
		   double *R,
		   double *DK,
		   double dAL,
		   double DET,
		   double DET0,
		   double dlm0,
		   double nor_min,
		   double *dR)
{
  double DLM=0;
  return (DLM);
}

double D_lam_ALM2 (double *BS_rr,
		   double *BS_D_R,
		   double *BS_R,
		   double *BS_DK,
		   double dlm,
		   double lm,
		   double dAL,
		   long ne,
		   int n_be,
		   long ndofd,
		   long npres,
		   double *BS_d_r,
		   double *r,
		   NODE *node,
		   ELEMENT *elem,
		   BOUNDING_ELEMENT *b_elems,
		   MATGEOM matgeom,
		   HOMMAT *hommat,
		   SUPP sup,
		   EPS *eps,
		   SIG *sig,
		   double nor_min,
		   CRPL *crpl,
		   double dt,
		   double stab,
		   long nce,
		   COEL *coel,
		   long *DomDof,
		   int GDof,
		   COMMUN comm 
		   /*,
		     GNOD *gnod,
		     GEEL *geel*/  ,
		   MPI_Comm mpi_comm,
		   const PGFem3D_opt *opts)
{
  double t = 0.0;
  double alpha_alpha = 0.0;
  double *r_n = NULL;
  double *r_n_1 = NULL;
  
  long i;
  double /* s1,s2, */b,a1,a2,a3,DLM,x1,x2,*p1,*p2,
    *p3,an1,an2,nor1,nor2,*BS_f_u,*f_u,/* tmp, */*p12;
  int myrank,nproc;

  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);

  DLM = 0.0;
  b = 0.0;
  
  p1 = (double *)aloc1 (DomDof[myrank]);
  p2 = (double *)aloc1 (DomDof[myrank]);
  p3 = (double *)aloc1 (DomDof[myrank]); 
  f_u = (double *)aloc1 (ndofd);
  BS_f_u = (double *)aloc1 (DomDof[myrank]);
  p12 = (double*)aloc1(ndofd);
  
  for (i=0;i<DomDof[myrank];i++) {
    p1[i] = BS_rr[i]/BS_DK[i];
    p2[i] = (BS_d_r[i] + BS_D_R[i])/BS_DK[i];
    p3[i] = BS_d_r[i]/BS_DK[i];
  }

  double send[4], rec[4];
  double p1p1, RR, p1p2, p2p2, p1p3, p2p3, p3p3;

  /*** Can improve performance by using a matrix-matrix product and
       using the appropriate components of the reslting matrix ***/

  /* Compute dot products  and pack */
  send[0] = p1p1 = ss (p1,p1,DomDof[myrank]);
  send[1] = RR	 = ss (BS_R,BS_R,DomDof[myrank]);
  send[2] = p1p2 = ss (p1,p2,DomDof[myrank]);
  send[3] = p2p2 = ss (p2,p2,DomDof[myrank]);

  /* MPI_Allreduce */
  MPI_Allreduce(send,rec,4,MPI_DOUBLE,MPI_SUM,mpi_comm);

  /* unpack recieve container */
  p1p1 = rec[0];
  RR   = rec[1];
  p1p2 = rec[2];
  p2p2 = rec[3];

  /* A1 */
  a1 = p1p1 + b*b*RR;

  /* A2 */
  a2 = 2.*(b*b*dlm*RR + p1p2);

  /* A3 */
  a3 = p2p2 - dAL*dAL + b*b*dlm*dlm*RR;

  /* /\******\/ */
  /* /\* A1 *\/ */
  /* /\******\/ */
  /* tmp = ss (p1,p1,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s1,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */
  /* tmp = ss (BS_R,BS_R,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s2,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */
  
  /* a1 = s1 + b*b*s2; */
  
  /* /\******\/ */
  /* /\* A2 *\/ */
  /* /\******\/ */
  /* tmp = ss (p1,p2,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s1,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */
  /* a2 = 2.*(b*b*dlm*s2 + s1); */
  
  /* /\******\/ */
  /* /\* A3 *\/ */
  /* /\******\/ */
  /* tmp = ss (p2,p2,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s1,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */
  /* a3 = s1 - dAL*dAL + b*b*dlm*dlm*s2; */
  
  if ((a2*a2 - 4.*a1*a3) < 0.0 || 2.*a1 == 0.0){
    if (myrank == 0)
      PGFEM_printf ("(a2*a2 - 4.*a1*a3) = %12.12f || 2*a1 = %12.12f\n",
	      (a2*a2 - 4.*a1*a3),2.*a1); 
    DLM = 1./0.0;

    dealoc1 (p1);
    dealoc1 (p2);
    dealoc1 (p3);
    dealoc1 (p12);
    dealoc1 (f_u);
    dealoc1 (BS_f_u);

    return (DLM); /* ??? return inf ???*/
  }
  
  /**** x1,2 = (-b +- sqrt(b^2 - 4*a*c))/(2*a) ****/
  x1 = (-1.*a2 + sqrt (a2*a2 - 4.*a1*a3))/(2.*a1);
  x2 = (-1.*a2 - sqrt (a2*a2 - 4.*a1*a3))/(2.*a1);
  
  p1p1 = p2p2 = p3p3 = 0.0;
  for (i=0;i<DomDof[myrank];i++){
    p1[i] = BS_d_r[i] + BS_D_R[i] + x1*BS_rr[i];
    p2[i] = BS_d_r[i] + BS_D_R[i] + x2*BS_rr[i];

    /* added to reduce # MPI_Allreduce */
    p1p1 += p1[i]*p1[i];
    p2p2 += p2[i]*p2[i];
    p3p3 += p3[i]*p3[i];
  } 

  /* compute magnitude of vectors */  
  send[0] = p1p1;
  send[1] = p2p2;
  send[2] = p3p3;
  MPI_Allreduce(send,rec,3,MPI_DOUBLE,MPI_SUM,mpi_comm);
  p1p1 = rec[0];
  p2p2 = rec[1];
  p3p3 = rec[2];

  p1p1 = sqrt(p1p1);
  p2p2 = sqrt(p2p2);
  p3p3 = sqrt(p3p3);

  p1p3 = 0.0;
  p2p3 = 0.0;
  for(i=0; i<DomDof[myrank]; i++){
    p1p3 += p1[i]*p3[i];
    p2p3 += p2[i]*p3[i];
  }
  p1p3 /=(p1p1*p3p3);
  p2p3 /=(p2p2*p3p3);

  send[0] = p1p3;
  send[1] = p2p3;
  MPI_Allreduce(send,rec,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
  p1p3 = rec[0];
  p2p3 = rec[1];

  /* Cos angle */
  an1 = acos(p1p3);
  an2 = acos(p2p3);
  an1 *= 180.0/PII;
  an2 *= 180.0/PII;

  /* nor_vec (p1,DomDof[myrank],mpi_comm); */
  /* nor_vec (p2,DomDof[myrank],mpi_comm); */
  /* nor_vec (p3,DomDof[myrank],mpi_comm); */
  
  /* /\* Cos angle *\/ */
  /* tmp = ss (p1,p3,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s1,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */
  /* tmp = ss (p2,p3,DomDof[myrank]); */
  /* MPI_Allreduce(&tmp,&s2,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */

  /* an1 = acos(s1); */
  /* an2 = acos(s2); */
  /* an1 *= 180.0/PII; */
  /* an2 *= 180.0/PII; */
  
  /* Residuals */
  for (i=0;i<DomDof[myrank];i++){
    p1[i] = BS_d_r[i] + BS_D_R[i] + x1*BS_rr[i];
    p2[i] = BS_d_r[i] + BS_D_R[i] + x2*BS_rr[i];
  }

  for (i=0;i<ndofd;i++ ) f_u[i] = 0.0;
  GToL(p1,p12,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
  fd_residuals (f_u,ne,n_be,ndofd,npres,p12,r,node,elem,b_elems,matgeom,
		hommat,sup,eps,sig,nor_min,crpl,dt,t,stab,nce,
		coel /*,gnod,geel*/,mpi_comm,opts,alpha_alpha,r_n,r_n_1);

  LToG(f_u,BS_f_u,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);
  for (i=0;i<DomDof[myrank];i++) {
    p1[i] = (lm + dlm + x1)*BS_R[i] - BS_f_u[i];
  }
  
  /* nor1 = ss (p1,p1,DomDof[myrank]); */
  /* MPI_Allreduce(&nor1,&tmp,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */
  /* nor1 = sqrt(tmp); */

  for (i=0;i<ndofd;i++){
    f_u[i] = 0.0;
  }

  GToL(p2,p12,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);

  fd_residuals (f_u,ne,n_be,ndofd,npres,p12,r,node,elem,b_elems,matgeom,
		hommat,sup,eps,sig,nor_min,crpl,dt,t,stab,nce,
		coel /*,gnod,geel*/,mpi_comm,opts,alpha_alpha,r_n,r_n_1);

  LToG(f_u,BS_f_u,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);

  for (i=0;i<DomDof[myrank];i++){
    p2[i] = (lm + dlm + x2)*BS_R[i] - BS_f_u[i]; 
  }

  /* nor2 = ss (p2,p2,DomDof[myrank]); */
  /* MPI_Allreduce(&nor2,&tmp,1,MPI_DOUBLE,MPI_SUM,mpi_comm); */
  /* nor2 = sqrt(tmp); */

  send[0] = ss (p1,p1,DomDof[myrank]);
  send[1] = ss (p2,p2,DomDof[myrank]);
  MPI_Allreduce(send,rec,2,MPI_DOUBLE,MPI_SUM,mpi_comm); 
  nor1 = sqrt(rec[0]);
  nor2 = sqrt(rec[1]);
  
  if (myrank == 0)
    PGFEM_printf ("\ns1 = %12.12f : s2 = %12.12f || a1 = %5.5f :"
	    " a2 = %5.5f || n1 = %12.12f : n2 = %12.12f\n",
	    /* s1,s2 */p1p3,p2p3,an1,an2,nor1,nor2);
  
  if (an1 <= an2){
    if (nor1 < nor2) DLM = x1;
    else             DLM = x2;
  }
  else{
    if (nor1 > nor2) DLM = x2;
    else             DLM = x1;
  }

  if (myrank == 0)
    PGFEM_printf("x1 = %12.12f : x2= %12.12f || DLM = %12.12f\n\n",x1,x2,DLM); 
  
  dealoc1 (p1);
  dealoc1 (p2);
  dealoc1 (p3);
  dealoc1 (p12);
  dealoc1 (f_u);
  dealoc1 (BS_f_u);
  
  return (DLM);
}

/*********************************************************************/
/*********************************************************************/

double d_ALM4 (long ndofd,
	       double *BS_rr,
	       double *BS_DK,
	       double dlm,
	       long *DomDof,
	       MPI_Comm mpi_comm)
/*
  SIMO
*/
{
  long i;
  double s1,*p,DLM,tmp;
  int myrank;

  MPI_Comm_rank(mpi_comm,&myrank);

  DLM = 0.0;
  
  p = (double *)aloc1 (DomDof[myrank]);
  for (i=0;i<DomDof[myrank];i++)
    p[i] = BS_rr[i]/BS_DK[i];
  
  tmp = ss (p,p,DomDof[myrank]);
  MPI_Allreduce(&tmp,&s1,1,MPI_DOUBLE,MPI_SUM,mpi_comm);

  DLM = sqrt (dlm*dlm*s1 + dlm*dlm);
  
  dealoc1 (p);
  
  return (DLM);
}

double d_lam_ALM4 (long ndofd,
		   double *BS_rr,
		   double *BS_DK,
		   double *BS_dR,
		   double dAL,
		   long *DomDof,
		   MPI_Comm mpi_comm)
/*
  SIMO
*/
{
  long i;
  double s,s1,DLM,*p1,*p2,*p3,znam,an1,an2,tmp1;

  double p1p1,p2p2,p3p3,p1p3,p2p3;
  double rec[3], send[3];

  int myrank;
  MPI_Comm_rank(mpi_comm,&myrank);

  /* DLM = s = s1 = s2 = tmp1 = 0.; */
  tmp1 = 0.0;
  znam = 1.;
  
  p1 = (double *)aloc1 (DomDof[myrank]);
  p2 = (double *)aloc1 (DomDof[myrank]);
  p3 = (double *)aloc1 (DomDof[myrank]); 
  
  for (i=0;i<DomDof[myrank];i++) {
    p1[i] = BS_rr[i]/BS_DK[i];
    p3[i] = BS_dR[i];
    tmp1 += BS_dR[i];
  }
  
  send[0] = ss (p1,p1,DomDof[myrank]); 
  send[1] = tmp1;
  MPI_Allreduce(send,rec,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
  s = rec[0];
  s1 = rec[1];

  /* For dR == 0.0 */
  if (s1 == 0.0){ /* <<=== BUG ??? (compare to 0.0) */
    dealoc1 (p1);
    dealoc1 (p2);
    dealoc1 (p3);
    return (sqrt(dAL*dAL/(s + 1.)));
  }
  
  p1p1 = p2p2 = p3p3 = 0.0;
  for (i=0;i<DomDof[myrank];i++){
    p1[i] =  BS_rr[i]; /* + */
    p2[i] = -BS_rr[i]; /* - */

    /* added to reduce # MPI_Allreduce */
    p1p1 += p1[i]*p1[i];
    p2p2 += p2[i]*p2[i];
    p3p3 += p3[i]*p3[i];
  }

  /* compute magnitude of vectors */  
  send[0] = p1p1;
  send[1] = p2p2;
  send[2] = p3p3;
  MPI_Allreduce(send,rec,3,MPI_DOUBLE,MPI_SUM,mpi_comm);
  p1p1 = sqrt(rec[0]);
  p2p2 = sqrt(rec[1]);
  p3p3 = sqrt(rec[2]);

  /* compute dot product */
  p1p3 = 0.0;
  p2p3 = 0.0;
  for(i=0; i<DomDof[myrank]; i++){
    p1p3 += p1[i]*p3[i];
    p2p3 += p2[i]*p3[i];
  }

  /* normalize */
  p1p3 /=(p1p1*p3p3);
  p2p3 /=(p2p2*p3p3);

  send[0] = p1p3;
  send[1] = p2p3;
  MPI_Allreduce(send,rec,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
  p1p3 = rec[0];
  p2p3 = rec[1];

  /* Cos angle */
  an1 = acos (p1p3);
  an2 = acos (p2p3);

  if (an1 <= an2) znam = +1.;
  else            znam = -1.;
  
  DLM = znam*sqrt(dAL*dAL/(s + 1.));
  
  dealoc1 (p1);
  dealoc1 (p2);
  dealoc1 (p3);
  
  return (DLM);
}

double D_lam_ALM4 (long ndofd,
		   double *BS_rr,
		   double *BS_d_r,
		   double *BS_D_R,
		   double *BS_DK,
		   double dlm,
		   double dAL,
		   long *DomDof,
		   MPI_Comm mpi_comm)
/*
  SIMO  
*/
{
  long i;
  double s1,s2,DLM,*p1,*p2,*p3;
  int myrank;

  MPI_Comm_rank(mpi_comm,&myrank);

  DLM = 0.0;
  
  p1 = (double *)aloc1 (DomDof[myrank]);
  p2 = (double *)aloc1 (DomDof[myrank]);
  p3 = (double *)aloc1 (DomDof[myrank]);
  
  /* Scale deformation vectors */
  for (i=0;i<DomDof[myrank];i++) {
    p1[i] = BS_rr[i]/BS_DK[i];
    p2[i] = BS_d_r[i]/BS_DK[i];
    p3[i] = BS_D_R[i]/BS_DK[i];
  }
  
  double send[2], rec[2];
  send[0] = ss (p2,p3,DomDof[myrank]);
  send[1] = ss (p1,p2,DomDof[myrank]);
  MPI_Allreduce(send,rec,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
  s1 = rec[0];
  s2 = rec[1];

  DLM = -s1/(s2 + dlm);
  
  dealoc1 (p1);
  dealoc1 (p2);
  dealoc1 (p3);
  
  return (DLM);
}
