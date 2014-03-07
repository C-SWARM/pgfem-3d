/* HEADER */
#include "new_potentials.h"
#include <math.h>
#include <string.h>
#include "mkl_cblas.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef INDEX_MACROS_H
#include "index_macros.h"
#endif

inline static double del(const int i, const int j)
{
  return ( (i==j)? 1.0:0.0 );
}

/*==== Get functions pointers ====*/

devPotentialFuncPtr getDevPotentialFunc(const int flag,
					HOMMAT const * mat)
{
  switch (mat->devPotFlag){
  case 0:
    if(mat->nu > 0.45){
      return devPotential_Mooney_Rivlin;
    } else {
      return devPotential_Linear;
    }

  case 1:
    return devPotential_Mooney_Rivlin;

  case 2:
    return devPotential_Linear;

  default:
    PGFEM_printf("ERROR: Unrecognized deviatoric potential flag (%d)\n",
	   (mat->devPotFlag));
    abort();
  }
}

devStressFuncPtr getDevStressFunc(const int flag,
				  HOMMAT const *mat)
{
  /* The mat->devPotFlag is used but defers to the poisson ratio
     in the case that it is not set. The integer flag is unused, but
     available for future use. */

  switch (mat->devPotFlag) {
  case 0:
    if(mat->nu > 0.45) {
      return devStress_Mooney_Rivlin;
    } else {
      return devStress_Linear;
    }

  case 1:
    return devStress_Mooney_Rivlin;

  case 2:
    return devStress_Linear;

  default:
    PGFEM_printf("ERROR: Unrecognized deviatoric potential flag (%d)\n",
	   (mat->devPotFlag));
    abort();
  }
}

matStiffFuncPtr getMatStiffFunc(const int flag,
				HOMMAT const *mat)
{
   /* The mat->devPotFlag is used but defers to the poisson ratio
     in the case that it is not set. The integer flag is unused, but
     available for future use. */
  switch (mat->devPotFlag) {
  case 0:
    if(mat->nu > 0.45) {
      return matStiff_Mooney_Rivlin;
    } else {
      return matStiff_Linear;
    }

  case 1:
    return matStiff_Mooney_Rivlin;

  case 2:
    return matStiff_Linear;

  default:
    PGFEM_printf("ERROR: Unrecognized deviatoric potential flag (%d)\n",
	   (mat->devPotFlag));
    abort();
  }
}

matSensFuncPtr getMatSensFunc(const int flag,
				HOMMAT const *mat)
{
   /* The mat->devPotFlag is used but defers to the poisson ratio
     in the case that it is not set. The integer flag is unused, but
     available for future use. */
  switch (mat->devPotFlag) {
  case 0:
    if(mat->nu > 0.45) {
      return matSens_Mooney_Rivlin;
    } else {
      return matSens_Linear;
    }

  case 1:
    return matSens_Mooney_Rivlin;

  case 2:
    return matSens_Linear;

  default:
    PGFEM_printf("ERROR: Unrecognized deviatoric potential flag (%d)\n",
	   (mat->devPotFlag));
    abort();
  }
}

UFuncPtr getUFunc(const int flag,
		  HOMMAT const *mat)
{
  /* The mat->volPotFlag is used but defers to the hard-coded
     flag in the case that it is not set. */
  int type;
  if (mat->volPotFlag > 0) {
    type = mat->volPotFlag;
  } else {
    type = flag;
  }

  switch(type){
  case 1:
    return U_Common;

  case 2:
    return U_Doll_Schweizerhof_7;

  case 3:
    return U_Doll_Schweizerhof_8;

  default:
     PGFEM_printf("ERROR: Unrecognized volumetric potential flag (%d)\n",(flag));
    abort();
  }
}

dUdJFuncPtr getDUdJFunc(const int flag,
			HOMMAT const *mat)
{
  /* The mat->volPotFlag is used but defers to the hard-coded
     flag in the case that it is not set. */
  int type;
  if (mat->volPotFlag > 0) {
    type = mat->volPotFlag;
  } else {
    type = flag;
  }

  switch(type){
  case 1:
    return dUdJ_Common;

  case 2:
    return dUdJ_Doll_Schweizerhof_7;

  case 3:
    return dUdJ_Doll_Schweizerhof_8;

  default:
    PGFEM_printf("ERROR: Unrecognized volumetric potential flag (%d)\n",(flag));
    abort();
  }
}

d2UdJ2FuncPtr getD2UdJ2Func(const int flag,
			    HOMMAT const *mat)
{
  /* The mat->volPotFlag is used but defers to the hard-coded
     flag in the case that it is not set. */
  int type;
  if (mat->volPotFlag > 0) {
    type = mat->volPotFlag;
  } else {
    type = flag;
  }

  switch(type){
  case 1:
    return d2UdJ2_Common_new;

  case 2:
    return d2UdJ2_Doll_Schweizerhof_7;

  case 3:
    return d2UdJ2_Doll_Schweizerhof_8;

  default:
    PGFEM_printf("ERROR: Unrecognized volumetric potential flag (%d)\n",(flag));
    abort();
  }
}

d2UdJ2FuncPtr getD3UdJ3Func(const int flag,
			    HOMMAT const *mat)
{
  /* The mat->volPotFlag is used but defers to the hard-coded
     flag in the case that it is not set. */
  int type;
  if (mat->volPotFlag > 0) {
    type = mat->volPotFlag;
  } else {
    type = flag;
  }

  switch(type){
  case 1:
    return d3UdJ3_Common_new;

  case 2:
    return d3UdJ3_Doll_Schweizerhof_7;

  case 3:
    return d3UdJ3_Doll_Schweizerhof_8;

  default:
    PGFEM_printf("ERROR: Unrecognized volumetric potential flag (%d)\n",(flag));
    abort();
  }
}

/*==== Deviatoric potential functions ====*/
void devPotential_Mooney_Rivlin(double const *C,
				HOMMAT const *mat,
				double *W)
{
  /* W_hat(C) = W(C_hat) */
  double *C_hat = aloc1(9);
  memset(C_hat,0,9*sizeof(double));
  cblas_daxpy(9,pow(det3x3(C),-1./3.),C,1,C_hat,1);

  const double trCh = C_hat[0] + C_hat[4] + C_hat[8];
  const double ChCh = cblas_ddot(9,C_hat,1,C_hat,1);

  *W = mat->m10*(trCh-3.) + 0.5*mat->m01*(trCh*trCh - ChCh-6);

  free(C_hat);
}

void devPotential_Linear(double const *C,
			 HOMMAT const *mat,
			 double *W)
{
  /* W_hat = int(S_hat)dC || W_hat = 0|C=I */
  /* S_hat = 2G E || E = 1/2 (C-1) */

 const double CC = cblas_ddot(9,C,1,C,1);
 const double trC = C[0] + C[4] + C[8];

 *W = 0.25*mat->G*(CC - 2*trC + 3);
}

/*==== Deviatoric stress functions ====*/
void devStress_Mooney_Rivlin(double const *C,
			     HOMMAT const *mat,
			     double *S)
{
  double *invC, *ident;
  invC = aloc1(9);
  ident = aloc1(9);

  inverse(C,3,invC);
  ident[0] = ident[4] = ident[8] = 1.0;

  double detC, trC, CC;
  detC = det3x3(C);
  trC = C[0] + C[4] + C[8];
  CC = cblas_ddot(9,C,1,C,1);

  for(int i=0; i<9; i++){
    S[i] =( 
	   (2*mat->m10*pow(detC,-1./3.)
	    *(ident[i] - 1./3.*trC*invC[i]))
	    
	   +(2*mat->m01*pow(detC,-2./3.)
	     *(trC*ident[i] - C[i]
	       + 1./3.*(CC-trC*trC)*invC[i]))
	    );
  }

 free(invC);
 free(ident);
}

void devStress_Linear(double const *C,
		      HOMMAT const *mat,
		      double *S)
{
  double *ident;
  ident = aloc1(9);
  ident[0] = ident[4] = ident[8] = 1.0;

  for(int i=0; i<9; i++){
    S[i] = mat->G*(C[i]-ident[i]);
  }

  free(ident);
}

/*==== Material stiffness functions ====*/
void matStiff_Mooney_Rivlin(double const *C,
			    HOMMAT const *mat,
			    double *L)
{
  double *invC, *ident;
  invC = aloc1(9);
  ident = aloc1(9);

  inverse(C,3,invC);
  ident[0] = ident[4] = ident[8] = 1.0;

  double detC1,detC2, trC, CC;
  detC1 = pow((det3x3(C)),-1./3.);
  detC2 = detC1*detC1;
  trC = C[0] + C[4] + C[8];
  CC = cblas_ddot(9,C,1,C,1);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	for(int l=0; l<3; l++){
	  L[idx_4(i,j,k,l)] = (4*mat->m10*detC1
			       *(
			      	 (1./9.*trC*invC[idx_2(i,j)]
			      	  - 1./3.*ident[idx_2(i,j)])*invC[idx_2(k,l)]
				 
			      	 + 1./3.*(trC*invC[idx_2(i,k)]*invC[idx_2(l,j)]
			      		  -ident[idx_2(k,l)]*invC[idx_2(i,j)])
			      	 )

			       + (4*mat->m01*detC2
			      	  *(
			      	    ident[idx_2(i,j)]*ident[idx_2(k,l)]
			      	    -ident[idx_2(i,k)]*ident[idx_2(l,j)]
				    
			      	    + 2./3.*((C[idx_2(i,j)]
			      		      -trC*ident[idx_2(i,j)])
			      		     *invC[idx_2(k,l)])

			      	    - ((CC-trC*trC)
			      	       *(2./9.*invC[idx_2(i,j)]*invC[idx_2(k,l)]
			      		 + 1./3.*invC[idx_2(i,k)]*invC[idx_2(l,j)]))

			      	    + 2./3.*((C[idx_2(k,l)]
			      		      -trC*ident[idx_2(k,l)])
			      		     *invC[idx_2(i,j)])
			      	    )
			      	  )
			       );
	}
      }
    }
  }

  free(invC);
  free(ident);
}

void matStiff_Linear(double const *C,
		     HOMMAT const *mat,
		     double *L)
{
  double *ident;
  ident = aloc1(9);
  ident[0] = ident[4] = ident[8] = 1.0;
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	for(int l=0; l<3; l++){
	  L[idx_4(i,j,k,l)] = (2*mat->G*ident[idx_2(i,k)]
			       *ident[idx_2(j,l)]);
	}
      }
    }
  }

  free(ident);
}

/*==== Material sensitivity functions ====*/
void matSens_Mooney_Rivlin(double const *C,
			   HOMMAT const *mat,
			   double *K)
{
  double *c_inv = aloc1(9);
  inverse(C,3,c_inv);
  const double *C_I = c_inv; /* get constatnt pointer */
  const double J23 = pow(det3x3(C),-1./3.);
  const double J43 = J23*J23;
  const double trC = C[0]+C[4]+C[8];
  const double CC = cblas_ddot(9,C,1,C,1);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      const int ij = idx_2(i,j);
      for(int k=0; k<3; k++){
	const int ik = idx_2(i,k);
	for(int l=0; l<3; l++){
	  const int kl = idx_2(k,l);
	  const int lj = idx_2(l,j);
	  for(int r=0; r<3; r++){
	    const int ir = idx_2(i,r);
	    const int kr = idx_2(k,r);
	    const int lr = idx_2(l,r);
	    for(int s=0; s<3; s++){
	      const int rs = idx_2(r,s);
	      const int sj = idx_2(s,j);
	      const int sk = idx_2(s,k);
	      const int sl = idx_2(s,l);
	      const int ijklrs = idx_6(i,j,k,l,r,s);

	      /*============ mu_10 term ==================*/
	      double M_10_term = 0.0;
	      M_10_term += 1./9.*(del(r,s)*C_I[ij] - trC*C[ir]*C_I[sj])*C_I[kl];

	      M_10_term -= (1./9.*trC*C_I[ij]-1./3.*del(i,j))*C_I[kr]*C_I[sl];

	      M_10_term += 1./3.*(del(r,s)*C_I[ik]*C_I[lj]
				  - trC*C_I[ir]*C_I[sk]*C_I[lj]
				  - trC*C_I[ik]*C_I[lr]*C_I[sj]
				  + del(k,l)*C_I[ir]*C_I[sj]);

	      M_10_term -= (1./27.*trC*C_I[ij]-1./9.*del(i,j))*C_I[kl]*C_I[rs];

	      M_10_term -= 1./9.*(trC*C_I[ik]*C_I[lj]-del(k,l)*C_I[ij])*C_I[rs];

	      M_10_term *= 8.*mat->m10*J23;

	      /*============ mu_01 term ==================*/
	      double M_01_term = 0.0;
	      M_01_term += 2./3.*(del(i,r)*del(s,j) - del(r,s)*del(i,j))*C_I[kl];

	      M_01_term -= 2./3.*(C[ij]-trC*del(i,j))*C_I[kr]*C_I[sl];

	      M_01_term -= 2.*(C[rs]-trC*del(r,s))*(2./9.*C_I[ij]*C_I[kl]
						    + 1./3.*C_I[ik]*C_I[lj]);

	      M_01_term += (CC-trC*trC)*(2./9.*C_I[ir]*C_I[sj]*C_I[kl]
					 + 2./9.*C_I[ij]*C_I[kr]*C_I[sl]
					 + 1./3.*C_I[ir]*C_I[sk]*C_I[lj]
					 + 1./3.*C_I[ik]*C_I[lr]*C_I[sj]);

	      M_01_term += 2./3.*(del(k,r)*del(s,l) - del(r,s)*del(k,l))*C_I[ij];

	      M_01_term -= 2./3.*(C[kl]-trC*del(k,l))*C_I[ir]*C_I[sj];

	      M_01_term -= 2./3.*(del(i,j)*del(k,l)-del(i,k)*del(l,j))*C_I[rs];

	      M_01_term -= 4./9.*(C[ij]-trC*del(i,j))*C_I[kl]*C_I[rs];

	      M_01_term += 2./3.*(CC-trC*trC)*(2./9.*C_I[ij]*C_I[kl]
					       + 1./3.*C_I[ik]*C_I[lj])*C_I[rs];

	      M_01_term -= 4./9.*(C[kl]-trC*del(k,l))*C_I[ij]*C_I[rs];

	      M_01_term *= 8.*mat->m01*J43;

	      /*============ K ==================*/
	      K[ijklrs] = M_10_term + M_01_term;
	    }
	  }
	}
      }
    }
  }

  free(c_inv);
}

void matSens_Linear(double const *C,
		    HOMMAT const *mat,
		    double *K)
{
  memset(K,0,729*sizeof(double));
}

/*==== Volumetric Potetntial Functions ====*/
void U_Common(double const J,
	      HOMMAT const *mat,
	      double *U)
{
  *U = (0.5*(J-1.0)*(J-1.0));
}

void U_Doll_Schweizerhof_7(double const J,
			   HOMMAT const *mat,
			   double *U)
{
  *U = (0.5*(exp(J-1.0)-log(J)-1.0));
}

void U_Doll_Schweizerhof_8(double const J,
			   HOMMAT const *mat,
			   double *U)
{
  *U = (0.5*(J-1.0)*log(J));
}

/*==== dUdJ functions ====*/
void dUdJ_Common(double const J,
		 HOMMAT const *mat,
		 double *dUdJ)
{
  (*dUdJ) = J-1;
}

void dUdJ_Doll_Schweizerhof_7(double const J,
			      HOMMAT const *mat,
			      double *dUdJ)
{
  (*dUdJ) = (exp(-1 + J) - 1./J)*.5;
}

void dUdJ_Doll_Schweizerhof_8(double const J,
			      HOMMAT const *mat,
			      double *dUdJ)
{
  (*dUdJ) = (-1 + J)/(2.*J) + log(J)*.5;
}

/*==== d2UdJ2 functions ===*/
void d2UdJ2_Common_new(double const J,
		       HOMMAT const *mat,
		       double *d2UdJ2)
{
  (*d2UdJ2) = 1.0;
}

void d2UdJ2_Doll_Schweizerhof_7(double const J,
				HOMMAT const *mat,
				double *d2UdJ2)
{
  (*d2UdJ2) = (exp(-1 + J) + 1./(J*J))*.5;
}

void d2UdJ2_Doll_Schweizerhof_8(double const J,
				HOMMAT const *mat,
				double *d2UdJ2)
{
  (*d2UdJ2) = 1./J - (-1 + J)/(2.*J*J);
}

/*==== d3UdJ3 functions ===*/
void d3UdJ3_Common_new(double const J,
		       HOMMAT const *mat,
		       double *d3UdJ3)
{
  *d3UdJ3 = 0.0;
}

void d3UdJ3_Doll_Schweizerhof_7(double const J,
				HOMMAT const *mat,
				double *d3UdJ3)
{
  *d3UdJ3 = (exp(-1. + J) - 2./(J*J*J))/2.;
}

void d3UdJ3_Doll_Schweizerhof_8(double const J,
				HOMMAT const *mat,
				double *d3UdJ3)
{
  *d3UdJ3 = (-1. + J)/(J*J*J) - 3./(2.*J*J);
}
