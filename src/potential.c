#include "potential.h"
#include <math.h>

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef DEF_GRAD_H
#include "def_grad.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#define DELTA(X,Y) ((double)((X)==(Y)))

deviatoricStressFunctionPtr getDeviatoricStressFunction(const int flag,
							const HOMMAT *material)
{
  /* The material->devPotFlag is used but defers to the poisson ratio
     in the case that it is not set. The integer flag is unused, but
     available for future use. */

  switch (material->devPotFlag) {
  case 0:
    if(material->nu > 0.45) {
      return deviatoricStress_MooneyRivlin;
    } else {
      return deviatoricStress_Linear;
    }

  case 1:
    return deviatoricStress_MooneyRivlin;

  case 2:
    return deviatoricStress_Linear;

  default:
    PGFEM_printf("ERROR: Unrecognized deviatoric potential flag (%d)\n",
	   (material->devPotFlag));
    abort();
  }
}

materialStiffnessFunctionPtr getMaterialStiffnessFunction(const int flag,
							  const HOMMAT *material)
{
   /* The material->devPotFlag is used but defers to the poisson ratio
     in the case that it is not set. The integer flag is unused, but
     available for future use. */
  switch (material->devPotFlag) {
  case 0:
    if(material->nu > 0.45) {
      return materialStiffness_MooneyRivlin;
    } else {
      return materialStiffness_Linear;
    }

  case 1:
    return materialStiffness_MooneyRivlin;

  case 2:
    return materialStiffness_Linear;

  default:
    PGFEM_printf("ERROR: Unrecognized deviatoric potential flag (%d)\n",
	   (material->devPotFlag));
    abort();
  }
}

volumetricPressureFunctionPtr getVolumetricPressureFunction(const int flag,
							    const HOMMAT *material)
{
  /* The material->volPotFlag is used but defers to the hard-coded
     flag in the case that it is not set. */
  int type;
  if (material->volPotFlag > 0) {
    type = material->volPotFlag;
  } else {
    type = flag;
  }

  switch(type){
  case 1:
    return volumetricPressure_Common;

  case 2:
    return volumetricPressure_DollSchweizerhof_7;

  case 3:
    return volumetricPressure_DollSchweizerhof_8;

  default:
     PGFEM_printf("ERROR: Unrecognized volumetric potential flag (%d)\n",(flag));
    abort();
  }
}

linearizedPressureFunctionPtr getLinearizedPressureFunction(const int flag,
							    const HOMMAT *material)
{
  /* The material->volPotFlag is used but defers to the hard-coded
     flag in the case that it is not set. */
  int type;
  if (material->volPotFlag > 0) {
    type = material->volPotFlag;
  } else {
    type = flag;
  }

  switch(type){
  case 1:
    return linearizedPressure_Common;

  case 2:
    return linearizedPressure_DollSchweizerhof_7;

  case 3:
    return linearizedPressure_DollSchweizerhof_8;

  default:
    PGFEM_printf("ERROR: Unrecognized volumetric potential flag (%d)\n",(flag));
    abort();
  }
}

d2UdJ2FunctionPtr getd2UdJ2Function(const int flag,
				    const HOMMAT *material)
{
  /* The material->volPotFlag is used but defers to the hard-coded
     flag in the case that it is not set. */
  int type;
  if (material->volPotFlag > 0) {
    type = material->volPotFlag;
  } else {
    type = flag;
  }

  switch(type){
  case 1:
    return d2UdJ2_Common;

  case 2:
    return d2UdJ2_DollSchweizerhof_7;
    
  case 3:
    return d2UdJ2_DollSchweizerhof_8;

  default:
    PGFEM_printf("ERROR: Unrecognized volumetric potential flag (%d)\n",(flag));
    abort();
  }
}

totalUpFunctPtr getTotalUpFunction(const int flag,
				   const HOMMAT *material)
{
  /* The material->volPotFlag is used but defers to the hard-coded
     flag in the case that it is not set. */
  int type;
  if (material->volPotFlag > 0) {
    type = material->volPotFlag;
  } else {
    type = flag;
  }

  switch(type){
  case 1:
    return totalUp_Common;

  case 2:
    return totalUp_DS7;
    
  case 3:
    return  totalUp_DS8;

  default:
    PGFEM_printf("ERROR: Unrecognized volumetric potential flag (%d)\n",(flag));
    abort();
  }
}

totalUppFunctPtr getTotalUppFunction(const int flag,
				     const HOMMAT *material)
{
  /* The material->volPotFlag is used but defers to the hard-coded
     flag in the case that it is not set. */
  int type;
  if (material->volPotFlag > 0) {
    type = material->volPotFlag;
  } else {
    type = flag;
  }

  switch(type){
  case 1:
    return totalUpp_Common;

  case 2:
    return totalUpp_DS7;
    
  case 3:
    return  totalUpp_DS8;

  default:
    PGFEM_printf("ERROR: Unrecognized volumetric potential flag (%d)\n",(flag));
    abort();
  }
}

/************************************************************************/
/******************** Deviatoric Stress Functions ***********************/
/************************************************************************/
void deviatoricStress_MooneyRivlin(const double * const *C,
				   const HOMMAT *material,
				   double **S)
{
  /* S= 2m10 det(C)^-1/3 (I-1/3tr(C)C^-T)+ 2m01 det(C)^-2/3 (tr(C)I-C+1/3(C:C-tr(C)^2)C^-T) */

  double detC;
  detC = def_grad_det(C);

  double trC;
  trC = C[0][0] + C[1][1] + C[2][2];

  // CC = C:C
  double CC = 0.0;
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      CC += C[i][j]*C[i][j];
    }
  }

  /* Get C_I = C^-1 */
  double **C_I;
  C_I = aloc2(3,3);
  def_grad_inv(C,C_I);

  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      S[i][j] = (2*material->m10*pow(detC,-1./3.)
		 *(DELTA(i,j) - 1./3.*trC*C_I[j][i]) + 
		 2*material->m01*pow(detC,-2./3.)*
		 (trC*DELTA(i,j)-C[i][j] 
		  + 1./3.*(CC - trC*trC)*C_I[j][i]));
    }
  }
  dealoc2(C_I,3);
}

void deviatoricStress_Linear(const double * const *C,
			     const HOMMAT *material,
			     double **S)
{
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) { 
      S[i][j] = material->G*(C[i][j] - DELTA(i,j));
    }
  }
}

/************************************************************************/
/******************* Material Stiffness Functions ***********************/
/************************************************************************/
void materialStiffness_MooneyRivlin(const double * const *C,
				    const HOMMAT *material,
				    const double * const *Fn,
				    const double * const *Fr,
				    double ****L)
{
  // Material stiffness L[I][J][k][K] = dSdC[I][J][M][N] * dCdFr[M][N][k][K]

  // Compute needed values/tensors
  double trC;
  trC = C[0][0] + C[1][1] + C[2][2];

  double CC; // C:C
  CC = 0.0;
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      CC += C[i][j] * C[i][j];
    }
  }

  double detC;
  detC = def_grad_det(C);

  double  **C_I; // C^-1
  C_I = aloc2(3,3);
  def_grad_inv(C,C_I);

  // Compute dSdC and dCdFr
  double ****dSdC;
  double ****dCdFr;
  dSdC = aloc4(3,3,3,3);
  dCdFr = aloc4(3,3,3,3);

  for(int I=0; I<3; I++) {
    for(int J=0; J<3; J++) {
      for(int M=0; M<3; M++) {
	for(int N=0; N<3; N++) {

	  // dSdC
	  dSdC[I][J][M][N] = 2*material->m10*pow(detC,-1./3.)*
	    ((1./9.*trC*C_I[J][I] - 1./3.*DELTA(I,J))*C_I[N][M]
	     + 1./3.*(trC*C_I[J][M]*C_I[N][I] - DELTA(M,N)*C_I[J][I]))
	    + 2*material->m01*pow(detC,-2./3.)*
	    (DELTA(M,N)*DELTA(I,J) - DELTA(I,M)*DELTA(J,N)
	     + 2./3.*(C[I][J] - trC*DELTA(I,J))*C_I[N][M]
	     - (CC - trC*trC)*(2./9.*C_I[J][I]*C_I[N][M]
			       + 1./3.*C_I[J][M]*C_I[N][I])
	     + 2./3.*(C[M][N] - trC*DELTA(M,N))*C_I[J][I]);

	  // dCdFr
	  dCdFr[I][J][M][N] = 0.0;
	  for(int k=0; k<3; k++) {
	    for(int K=0; K<3; K++) {
	      dCdFr[I][J][M][N] += Fn[k][I]*(DELTA(N,k)*Fr[M][K] +
					     Fr[M][k]*DELTA(K,N))*Fn[K][J];
	    }
	  }
	}
      }
    }
  }

  // Compute L = dSdC*dCdFr
  for(int I=0; I<3; I++) {
    for(int J=0; J<3; J++) {
      for(int k=0; k<3; k++) {
	for(int K=0; K<3; K++) {

	  L[I][J][k][K] = 0.0;

	  for(int M=0; M<3; M++) {
	    for(int N=0; N<3; N++) {
	      L[I][J][k][K] += dSdC[I][J][M][N] * dCdFr[M][N][k][K];
	    }
	  }
	}
      }
    }
  }


  dealoc2(C_I,3);
  dealoc4(dSdC,3,3,3);
  dealoc4(dCdFr,3,3,3);
}

void materialStiffness_Linear(const double * const *C,
			      const HOMMAT *material,
			      const double * const *Fn,
			      const double * const *Fr,
			      double ****L)
{
  double ****dSdC;
  double ****dCdFr;
  dSdC = aloc4(3,3,3,3);
  dCdFr = aloc4(3,3,3,3);

  for(int I=0; I<3; I++) {
    for(int J=0; J<3; J++) {
      for(int M=0; M<3; M++) {
	for(int N=0; N<3; N++) {

	  dSdC[I][J][M][N] = material->G*DELTA(I,M)*DELTA(J,N);

	  dCdFr[I][J][M][N] = 0.0;
	  for(int k=0; k<3; k++) {
	    for(int K=0; K<3; K++) {
	     dCdFr[I][J][M][N] += Fn[k][I]*(DELTA(N,k)*Fr[M][K] +
					     Fr[M][k]*DELTA(K,N))*Fn[K][J];
	    }
	  }
	}
      }
    }
  }

  for(int I=0; I<3; I++) {
    for(int J=0; J<3; J++) {
      for(int k=0; k<3; k++) {
	for(int K=0; K<3; K++) {

	  L[I][J][k][K] = 0.0;

	  for(int M=0; M<3; M++) {
	    for(int N=0; N<3; N++) {
	      L[I][J][k][K] += dSdC[I][J][M][N] * dCdFr[M][N][k][K];
	    }
	  }
	}
      }
    }
  }

  dealoc4(dSdC,3,3,3);
  dealoc4(dCdFr,3,3,3);
}

/************************************************************************/
/******************* Volumetric Pressure Functions **********************/
/************************************************************************/
void volumetricPressure_Common(const double Jn,
			       const double Jr,
			       const HOMMAT *material,
			       double *pres)
{
  // Evaluate updated pressure
  (*pres) = Jn*(Jr-1);
}

void totalUp_Common(const double Tn,
		    const double Tr,
		    const HOMMAT *material,
		    double *Up)
{
  (*Up) = Tr*Tn-1;
}

void volumetricPressure_DollSchweizerhof_7(const double Jn,
					   const double Jr,
					   const HOMMAT *material,
					   double *pres)
{
  // Evaluate updated pressure
  // (*pres) = ((Jr-1) + Jr*Jn*(exp(Jn*Jr-1) - exp(Jn-1)))/(2*Jr*Jn);
  // (*pres) = (-exp(-1 + Jn) + 1./Jn +exp(-1 + Jn*Jr)*Jn - 1./Jr)/2.;
  (*pres) = (-exp(-1 + Jn) + 1./Jn)/2. + (exp(-1 + Jn*Jr) - 1/(Jn*Jr))/2.;
}

void totalUp_DS7(const double Tn,
		 const double Tr,
		 const HOMMAT *material,
		 double *Up)
{
  (*Up) = (exp(-1 + Tn*Tr) - 1./(Tn*Tr))*.5;
}

void volumetricPressure_DollSchweizerhof_8(const double Jn,
					   const double Jr,
					   const HOMMAT *material,
					   double *pres)
{
  // Evaluate updated pressure
  // (*pres) = (Jr+Jn*Jr*log(Jr)-1)/(2*Jr*Jn);
  // (*pres) = (-1 + 1./Jn + Jn - 1./Jr - log(Jn) + Jn*log(Jn*Jr))/2.;
  (*pres) = (-1 + Jr - Jn*Jr*log(Jn) + Jn*Jr*log(Jn*Jr))/(2.*Jn*Jr);
}

void totalUp_DS8(const double Tn,
		 const double Tr,
		 const HOMMAT *material,
		 double *Up)
{
  (*Up) = (-1 + Tn*Tr)/(2.*Tn*Tr) + log(Tn*Tr)*.5;
}

/************************************************************************/
/******************* Linearized Pressure Functions **********************/
/************************************************************************/
void linearizedPressure_Common(const double Jn,
			       const double Jr,
			       const double *const *Fn,
			       const double *const *Fr,
			       const HOMMAT *material,
			       double **linPres)
{
  double **Fr_I;
  Fr_I = aloc2(3,3);
  def_grad_inv(Fr,Fr_I);

  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
       linPres[i][j] = Jr*Jn*Fr_I[j][i];
    }
  }
  dealoc2(Fr_I,3);
}

void linearizedPressure_DollSchweizerhof_7(const double Jn,
					   const double Jr,
					   const double *const *Fn,
					   const double *const *Fr,
					   const HOMMAT *material,
					   double **linPres)
{
  double **Fr_I;
  Fr_I = aloc2(3,3);
  def_grad_inv(Fr,Fr_I);

  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      /* linPres[i][j] = 0.5*(Jn*Jr*Fr_I[j][i]*exp(Jn*Jr-1)  */
      /* 			    + Fr_I[j][i]/(Jn*Jr)); */
      /* linPres[i][j] = (0.5*(exp(-1 + Jn*Jr)*Jn*Jn + 1./(Jr*Jr)) */
      /* 		       *Jr*Fr_I[j][i]); */
      linPres[i][j] = ((exp(-1 + Jn*Jr)*Jn + 1/(Jn*Jr*Jr))/2.
		       *Jr*Fr_I[j][i]);
    }
  }
  dealoc2(Fr_I,3);
}

void linearizedPressure_DollSchweizerhof_8(const double Jn,
					   const double Jr,
					   const double *const *Fn,
					   const double *const *Fr,
					   const HOMMAT *material,
					   double **linPres)
{
  double **Fr_I;
  Fr_I = aloc2(3,3);
  def_grad_inv(Fr,Fr_I);

  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      /* linPres[i][j] = Fr_I[j][i]*(1+1/(2*Jn*Jr)); */
      /* linPres[i][j] = 0.5*(1./(Jr*Jr) + Jn/Jr)*Jr*Fr_I[j][i]; */
      linPres[i][j] = ((1 + Jn*Jr)/(2.*Jn*Jr*Jr)
		       *Jr*Fr_I[j][i]);
    }
  }
  dealoc2(Fr_I,3);
}

/************************************************************************/
/*******************         d^2 U/d J^2           **********************/
/************************************************************************/
 void d2UdJ2_Common(const double Jn,
		     const double Jr,
		     const HOMMAT *material,
		     double *d2UdJ2)
{
  (*d2UdJ2) = Jn;
}

void totalUpp_Common(const double Tn,
		     const double Tr,
		     const HOMMAT *material,
		     double *Upp)
{
  (*Upp) = 1.0;
}

void d2UdJ2_DollSchweizerhof_7(const double Jn,
			       const double Jr,
			       const HOMMAT *material,
			       double *d2UdJ2)
{
  //(*d2UdJ2) = (exp(-1 + Jn*Jr)*Jn*Jn + 1./(Jr*Jr))/(2.*Jn);
  //(*d2UdJ2) = 0.5*(exp(-1 + Jn*Jr)*Jn*Jn + 1./(Jr*Jr));
  (*d2UdJ2) = (exp(-1 + Jn*Jr)*Jn + 1/(Jn*Jr*Jr))/2.;
}

void totalUpp_DS7(const double Tn,
		  const double Tr,
		  const HOMMAT *material,
		  double *Upp)
{
  (*Upp) = (exp(-1 + Tn*Tr) + 1./(Tn*Tn*Tr*Tr))*.5;
}

void d2UdJ2_DollSchweizerhof_8(const double Jn,
			       const double Jr,
			       const HOMMAT *material,
			       double *d2UdJ2)
{
  //(*d2UdJ2) = (1 + Jn*Jr)/(2.*Jn*Jr*Jr);
  //(*d2UdJ2) = 0.5*(1./(Jr*Jr) + Jn/Jr);
  (*d2UdJ2) = (1 + Jn*Jr)/(2.*Jn*Jr*Jr);
}

void totalUpp_DS8(const double Tn,
		  const double Tr,
		  const HOMMAT *material,
		  double *Upp)
{
  (*Upp) = 1./(Tn*Tr) - (-1 + Tn*Tr)/(2.*Tn*Tn*Tr*Tr);
}
