#include "TA_GA.h"
#include <math.h>

#ifndef INCL_H
#include "incl.h"
#endif

void TA_GA (long nss,long mat,CRPL *crpl,double Har,double **eFn,double **S,double *TA,double *GA)
{
  long k,M,N,P,Q;
  double **PP,**AA,pom;
 
  PP = aloc2 (3,3); AA = aloc2 (3,3);
 
  for (k=0;k<nss;k++){
    
    PP[0][0] = crpl[mat].P[k][3]*crpl[mat].P[k][0];  PP[0][1] = crpl[mat].P[k][3]*crpl[mat].P[k][1];  PP[0][2] = crpl[mat].P[k][3]*crpl[mat].P[k][2];
    PP[1][0] = crpl[mat].P[k][4]*crpl[mat].P[k][0];  PP[1][1] = crpl[mat].P[k][4]*crpl[mat].P[k][1];  PP[1][2] = crpl[mat].P[k][4]*crpl[mat].P[k][2];
    PP[2][0] = crpl[mat].P[k][5]*crpl[mat].P[k][0];  PP[2][1] = crpl[mat].P[k][5]*crpl[mat].P[k][1];  PP[2][2] = crpl[mat].P[k][5]*crpl[mat].P[k][2];
    
    TA[k] = 0.0;
    for (M=0;M<3;M++){
      for (N=0;N<3;N++){
	AA[M][N] = 0.0;
	for (P=0;P<3;P++){
	  for (Q=0;Q<3;Q++){
	    AA[M][N] += eFn[P][M]*eFn[P][Q]*S[Q][N];
	  }
	}
	
	TA[k] += AA[M][N]*PP[M][N];
      }
    }
    
    pom = fabs (TA[k]/Har);
    GA[k] = crpl[mat].a*(TA[k]/Har)*pow(pom,1./crpl[mat].m - 1.);
    
  }/* end k < nss */
  
  dealoc2 (PP,3); dealoc2 (AA,3);
}
