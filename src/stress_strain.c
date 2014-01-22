#include "stress_strain.h"
#include <math.h>

void get_GL_strain (double **Fn,double **Fr,double Jr,double Tr,double **E)
{
  long M,N,X,P,Q;
  double d_MN;
  
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      if (M == N) d_MN = 1.0; else d_MN = 0.0;
      E[M][N] = -1./2.*d_MN;
      for (X=0;X<3;X++){
	for (P=0;P<3;P++){
	  for (Q=0;Q<3;Q++){
	    E[M][N] += (1./2.*pow(Tr,+2./3.)*pow(Jr,-2./3.)
			* Fn[X][M]*Fr[P][X]*Fr[P][Q]*Fn[Q][N]);
	  }
	}
      }
    }
  }
}

void get_SPK_stress (double L[3][3][3][3],double **E,double **S)
{
  long M,N,X,P;
  
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      S[M][N] = 0.0;
      for (X=0;X<3;X++){
	for (P=0;P<3;P++){
	  if (L[M][N][X][P] == 0.0) continue;
	  S[M][N] += L[M][N][X][P]*E[X][P];
	}
      }
    }
  }
}
