#include "def_grad.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

static const int periodic = 0;

void def_grad_get (const long nne,
		   const long ndofn,
		   const double ****AA,
		   const double *r_e,
		   double **F)
{
  long i,j,o,g,ndn;
  double **u;
  double d_ij;
  
  /* 3D */ ndn = 3;
  
  u = aloc2 (ndn,nne);
  
  /* evaluate u */
  for (i=0;i<nne;i++){
    u[0][i] = r_e[i*ndofn+0];
    u[1][i] = r_e[i*ndofn+1];
    u[2][i] = r_e[i*ndofn+2];
  }
  
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){

      if (i == j)
	d_ij = 1.0; 
      else
	d_ij = 0.0;

      F[i][j] = d_ij;
      if (periodic == 1 /*|| gem == 1*/)
	F[i][j] = 0.0;
      for (g=0;g<ndn;g++){
	for (o=0;o<nne;o++){
	  F[i][j] += AA[i][j][g][o]*u[g][o];
	}
      }
    }
  }
  
  dealoc2 (u,ndn);
}

void def_grad_inv (const double * const *F,
		   double **F_I)
{
  double DET,A11,A12,A13,A21,A22,A23,A31,A32,A33;
  
  DET = def_grad_det(F);
  
  A11 = F[1][1]*F[2][2] - F[1][2]*F[2][1];
  A12 = F[0][2]*F[2][1] - F[0][1]*F[2][2];
  A13 = F[0][1]*F[1][2] - F[0][2]*F[1][1];
  
  A21 = F[1][2]*F[2][0] - F[1][0]*F[2][2];
  A22 = F[0][0]*F[2][2] - F[0][2]*F[2][0];
  A23 = F[0][2]*F[1][0] - F[0][0]*F[1][2];
  
  A31 = F[1][0]*F[2][1] - F[1][1]*F[2][0];
  A32 = F[0][1]*F[2][0] - F[0][0]*F[2][1];
  A33 = F[0][0]*F[1][1] - F[0][1]*F[1][0];
  
  
  F_I[0][0] = 1./DET*A11; F_I[0][1] = 1./DET*A12; F_I[0][2] = 1./DET*A13;
  F_I[1][0] = 1./DET*A21; F_I[1][1] = 1./DET*A22; F_I[1][2] = 1./DET*A23;
  F_I[2][0] = 1./DET*A31; F_I[2][1] = 1./DET*A32; F_I[2][2] = 1./DET*A33;
  
  /*
    for (i=0;i<3;i++){
    for (j=0;j<3;j++){
    for (k=0;k<3;k++){
    I[i][k] += F[i][j]*F_I[j][k];
    }
    }
    }
    
    for (i=0;i<3;i++){
    for (j=0;j<3;j++){
    if (i == j && (I[i][j] > 1.0000001 || I[i][j] < +0.999999)) {PGFEM_printf ("Error in def_grad_inv\n"); abort();}
    if (i != j && (I[i][j] > 0.0000001 || I[i][j] < -0.000001)) {PGFEM_printf ("Error in def_grad_inv\n"); abort();}
    }
    }
    
    dealoc2 (I,3);
  */
}

double def_grad_det (const double * const *F)
{
  double J;
  
  J = ((F[0][0]*F[1][1]*F[2][2] 
	+ F[0][1]*F[1][2]*F[2][0] 
	+ F[0][2]*F[1][0]*F[2][1]) - 
       (F[0][2]*F[1][1]*F[2][0] 
	+ F[0][0]*F[1][2]*F[2][1] 
	+ F[0][1]*F[1][0]*F[2][2]));
  
  return (J);
}
