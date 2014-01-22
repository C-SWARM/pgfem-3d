#include "Re1_Re2_Re3.h"
#include <math.h>

static const int periodic = 0;

/* VVolume is only used in periodic branch that is depricated. Set to
   0 -> poisoned value will result in NaN and prompt investigation. */
static const double VVolume = 0.0;

void Re1_Re2_Re3 (long nne,
		  long ndn,
		  long npres,
		  double ai,
		  double aj,
		  double ak,
		  double J,
		  double *Psi,
		  double Jr,
		  double Jn,
		  double Tr,
		  double Tn,
		  double **Fr_I,
		  double ****ST,
		  double ****FF,
		  double **S,
		  double **f,
		  double pp,
		  double **Re1,
		  double **re1,
		  double *Re2,
		  double *re2,
		  double *Re3,
		  double *re3)
{
  long M,N,U,V,P,R;
  double help1;
  
  /* Residual for displacements */
  for (M=0;M<ndn;M++){
    for (N=0;N<nne;N++){
      
      help1 = 0.0;
      for (U=0;U<3;U++){
	for (V=0;V<3;V++){
	  help1 += Fr_I[V][U]*ST[U][V][M][N];
	}
      }
      
      Re1[M][N] = pp*Jr*help1;
      
      for (P=0;P<3;P++){
	for (R=0;R<3;R++){
	  Re1[M][N] += 1./Jn*pow(Tr,2./3.)*FF[P][R][M][N]*S[P][R];
	}
      }
      
      if (periodic == 1)  re1[M][N] += 1./(VVolume) *ai*aj*ak*J* Re1[M][N];
      else                re1[M][N] +=               ai*aj*ak*J* Re1[M][N];
    }
  }
  
  /* Residuals for pressure and volume change */
  for (M=0;M<npres;M++){
    Re2[M] = 1./Jn*Psi[M]*(Jr*Jn - Tr*Tn);
    Re3[M] = -1./Jn*Tn*pp*Psi[M];
    for (N=0;N<3;N++){
      for (P=0;P<3;P++){
	Re3[M] += 1./Jn*1./3.*pow(Tr,-1./3)*Psi[M]*f[N][P]*S[N][P];
      }
    }
    if (periodic == 1){
      re2[M] += 1./(VVolume) *ai*aj*ak*J* Re2[M];
      re3[M] += 1./(VVolume) *ai*aj*ak*J* Re3[M];
    }
    else{
      re2[M] += ai*aj*ak*J* Re2[M];
      re3[M] += ai*aj*ak*J* Re3[M];
    }
  }  
}
