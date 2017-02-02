#include "tensors.h"
#include <math.h>

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef INCL_H
#include "incl.h"
#endif

#ifndef NULL_H
#include "null.h"
#endif

static const int periodic = 0;

/* VVolume only appears in deprecated 'periodic' branch. Poison value
   to return NaN, prompting investigation of how entered deprecated
   branch of code. */
static const double VVolume = 0.0;

void shape_tensor (long nne,
		   long ndofn,
		   double *N_x,
		   double *N_y,
		   double *N_z,
		   double ****AA)
{
  long i,m,g,o,ndn;
  double d_mg;
  
  /* 3D */ ndn = 3;
  
  for (o=0;o<nne;o++){
    for (g=0;g<ndn;g++){
      for (m=0;m<ndn;m++){
	if (g == m) d_mg = 1.0;  else d_mg = 0.0;
	for (i=0;i<ndn;i++){
	  if (i == 0)  AA[m][i][g][o] = N_x[o]*d_mg;
	  if (i == 1)  AA[m][i][g][o] = N_y[o]*d_mg;
	  if (i == 2)  AA[m][i][g][o] = N_z[o]*d_mg;
	}
      }
    }
  }
}

void matrix_tensor_3D (long mat,
		       HOMMAT *hommat,
		       double A[3][3][3][3])
{
  null_4d (A);

  A[0][0][0][0] = hommat[mat].L[0];
  A[0][0][1][1] = hommat[mat].L[1];
  A[0][0][2][2] = hommat[mat].L[2];

  A[1][1][0][0] = hommat[mat].L[1];
  A[1][1][1][1] = hommat[mat].L[3];
  A[1][1][2][2] = hommat[mat].L[4];
  
  A[2][2][0][0] = hommat[mat].L[2];
  A[2][2][1][1] = hommat[mat].L[4];
  A[2][2][2][2] = hommat[mat].L[5];
  
  A[1][2][1][2] = A[2][1][2][1] 
    = A[1][2][2][1] = A[2][1][1][2] = hommat[mat].L[6];

  A[0][2][0][2] = A[2][0][2][0]
    = A[0][2][2][0] = A[2][0][0][2] = hommat[mat].L[7];

  A[0][1][0][1] = A[1][0][1][0]
    = A[0][1][1][0] = A[1][0][0][1] = hommat[mat].L[8];
}

void tensors_FF_f (long nne,
		   long ndn,
		   double **Fn,
		   double **Fr,
		   double **Fr_I,
		   double Jr,
		   double ****ST,
		   double ****FF,
		   double **f)
{
  long M,N,X,P,V,W,Q,R,U;
  double help1;
  
  /* Tensors FF */
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      for (X=0;X<ndn;X++){
	for (P=0;P<nne;P++){
	  FF[M][N][X][P] = 0.0;
	  
	  help1 = 0.0;
	  for (V=0;V<3;V++){
	    for (W=0;W<3;W++){
	      if (ST[V][W][X][P] == 0.0) continue;
	      help1 += Fr_I[W][V]*ST[V][W][X][P];
	    }
	  }
	  
	  for (Q=0;Q<3;Q++){
	    for (R=0;R<3;R++){
	      for (U=0;U<3;U++){
		FF[M][N][X][P] += (pow(Jr,-2./3.)*Fn[Q][M] 
				   * (1./2.*(ST[R][Q][X][P]*Fr[R][U]
					     + Fr[R][Q]*ST[R][U][X][P])
				      - 1./3.*help1*Fr[R][Q]*Fr[R][U])
				   * Fn[U][N]);
	      }
	    }
	  }
	  
	}
      }
    }
  }
  
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      f[M][N] = 0.0;
      
      for (X=0;X<3;X++){
	for (P=0;P<3;P++){
	  for (Q=0;Q<3;Q++){
	    f[M][N] += pow(Jr,-2./3.) * Fn[X][M]*Fr[P][X]*Fr[P][Q]*Fn[Q][N];
	  }
	}
      }

    }
  }
}

void tensors_FFla (long ii,
		   long ip,
		   long nne,
		   long ndn,
		   EPS *eps,
		   double **Fn,
		   double **Fr,
		   double **Fr_I,
		   double Jr,
		   double **f,
		   double **UU,
		   double **Flam,
		   double **FFlam,
		   double **pGplam,
		   const int analysis)
{
  long M,N,V,W,Q,R,U,X;
  double help1;
  
  /* Tensors FF */
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      FFlam[M][N] = 0.0;
      if (analysis == FS_CRPL)
	pGplam[M][N] = 0.0;
      
      help1 = 0.0;
      for (V=0;V<3;V++){
	for (W=0;W<3;W++){
	  help1 += Fr_I[W][V]*Flam[V][W];
	}
      }
      
      for (Q=0;Q<3;Q++){
	for (R=0;R<3;R++){
	  for (U=0;U<3;U++){
	    FFlam[M][N] += (pow(Jr,-2./3.)*Fn[Q][M] 
			    * (1./2.*(Flam[R][Q]*Fr[R][U]
				      + Fr[R][Q]*Flam[R][U])
			       - 1./3.*help1*Fr[R][Q]*Fr[R][U])
			    * Fn[U][N]);
	   
	     if (analysis == FS_CRPL){
	       for (X=0;X<3;X++){
		if (f[Q][R] != 0.0
		    && eps[ii].il[ip].dUU_Fr[U][X][Q][M] != 0.0)
		  pGplam[M][N] += (f[Q][R]*eps[ii].il[ip].dUU_Fr[U][X][Q][M]
				   *UU[R][N]*Flam[U][X]);

		if (f[Q][R] != 0.0
		    && eps[ii].il[ip].dUU_Fr[U][X][R][N] != 0.0)
		  pGplam[M][N] += (f[Q][R]*UU[Q][M]
				   *eps[ii].il[ip].dUU_Fr[U][X][R][N]
				   *Flam[U][X]);
	      }
	    }

	  }
	}
      }
      
    }
  }/* end M < 3 */  
}

void tensors_aa_bb_dd_mm (long nne,
			  long ndn,
			  long npres,
			  double ai,
			  double aj,
			  double ak,
			  double J,
			  double *Psi,
			  double L[3][3][3][3],
			  double **Fn,
			  double **Fr,
			  double **Fr_I,
			  double Jn,
			  double Jr,
			  double Tn,
			  double Tr,
			  double ****ST,
			  double ****FF,
			  double **f,
			  double ***AA,
			  double ***aa,
			  double ***BB,
			  double ***bb,
			  double **DD,
			  double **dd,
			  double **MM,
			  double **mm)
{
  long M,N,X,P,V,W,R,Q;
  double d_VW,d_PR,help1;
  
  /* Tensors AA and BB*/
  for (Q=0;Q<ndn;Q++){
    for (M=0;M<npres;M++){
      for (N=0;N<nne;N++){
	AA[Q][M][N] = 0.0;
	BB[Q][M][N] = 0.0;
	
	help1 = 0.0;
	for (X=0;X<3;X++){
	  for (P=0;P<3;P++){
	    help1 += Fr_I[P][X]*ST[X][P][Q][N];
	    
	    for (V=0;V<3;V++){
	      for (W=0;W<3;W++){
		if (V == W)
		  d_VW = 1.0;
		else
		  d_VW = 0.0;
		BB[Q][M][N] += (1./Jn*4./3.*pow(Tr,-1./3.)
				*Psi[M]*FF[X][P][Q][N]*L[X][P][V][W]
				*0.5*(pow(Tr,2./3.)*f[V][W] - 0.5*d_VW));
	      }
	    }
	  }
	}
	
	AA[Q][M][N] = Jr*Psi[M]*help1;
	
	if (periodic == 1){
	  aa[Q][M][N] += 1./(VVolume) *ai*aj*ak*J * AA[Q][M][N];
	  bb[Q][M][N] += 1./(VVolume) *ai*aj*ak*J * BB[Q][M][N];
	}
	else{
	  aa[Q][M][N] += ai*aj*ak*J * AA[Q][M][N];
	  bb[Q][M][N] += ai*aj*ak*J * BB[Q][M][N];
	}
      }
    }
  }
  
  /* Tensors DD and MM */
  for (Q=0;Q<npres;Q++){
    for (M=0;M<npres;M++){
      DD[Q][M] = 0.0;  MM[Q][M] = 0.0;
      
      for (N=0;N<3;N++){
	for (X=0;X<3;X++){
	  for (P=0;P<3;P++){
	    for (R=0;R<3;R++){
	      if (P == R) d_PR = 1.0; else d_PR = 0.0;
	      MM[Q][M] += (1./Jn*1./9.*pow(Tr,-4./3.)
			   *Psi[Q]*Psi[M] * f[N][X]*L[N][X][P][R]
			   * 0.5*(pow(Tr,2./3.)*f[P][R] + d_PR));
	    }
	  }
	}
      }
      
      DD[Q][M]  = Tn/Jn*Psi[Q]*Psi[M];
      if (periodic == 1){
	dd[Q][M] += 1./(VVolume) *ai*aj*ak*J * DD[Q][M];
	mm[Q][M] += 1./(VVolume) *ai*aj*ak*J * MM[Q][M];
      }
      else{
	dd[Q][M] += ai*aj*ak*J * DD[Q][M];
	mm[Q][M] += ai*aj*ak*J * MM[Q][M];
      }
    }
  }
}

void tensors_GG (long ii,
		 long ip,
		 long nne,
		 long ndn,
		 double ****ST,
		 double ****FF,
		 double **f,
		 double **UU,
		 EPS *eps,
		 double ****GT,
		 double gT[3][3],
		 double ****gg)
{
  long M,N,P,Q,U,V,W,X;
  
  /* Tensors GT, gg and gT */

  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      
      gT[M][N] = 0.0;
      for (U=0;U<3;U++){
	for (V=0;V<3;V++){
	  if (f[U][V] == 0.0)
	    continue;
	  if (eps[ii].il[ip].dUU_Tr[U][M] != 0.0)
	    gT[M][N] += f[U][V]*eps[ii].il[ip].dUU_Tr[U][M]*UU[V][N];

	  if (eps[ii].il[ip].dUU_Tr[V][N] != 0.0)
	    gT[M][N] += f[U][V]*UU[U][M]*eps[ii].il[ip].dUU_Tr[V][N];
	}
      }
      
      for (P=0;P<ndn;P++){
	for (Q=0;Q<nne;Q++){
	  GT[M][N][P][Q] = 0.0;
	  gg[M][N][P][Q] = 0.0;

	  for (U=0;U<3;U++){
	    for (V=0;V<3;V++){
	      if (FF[U][V][P][Q] != 0.0
		  && eps[ii].il[ip].dUU_Tr[U][M] != 0.0)
		GT[M][N][P][Q] += (FF[U][V][P][Q]
				   *eps[ii].il[ip].dUU_Tr[U][M]*UU[V][N]);

	      if (FF[U][V][P][Q] != 0.0
		  && eps[ii].il[ip].dUU_Tr[V][N] != 0.0)
		GT[M][N][P][Q] += (FF[U][V][P][Q]*UU[U][M]
				   *eps[ii].il[ip].dUU_Tr[V][N]);
	      
	      for (W=0;W<3;W++){
		for (X=0;X<3;X++){
		  if (ST[W][X][P][Q] == 0.0 || f[U][V] == 0.0)
		    continue;
		  if (eps[ii].il[ip].dUU_Fr[W][X][U][M] != 0.0)
		    gg[M][N][P][Q] += (f[U][V]
				       *eps[ii].il[ip].dUU_Fr[W][X][U][M]
				       *UU[V][N]*ST[W][X][P][Q]);

		  if (eps[ii].il[ip].dUU_Fr[W][X][V][N] != 0.0)
		    gg[M][N][P][Q] += (f[U][V]*UU[U][M]
				       *eps[ii].il[ip].dUU_Fr[W][X][V][N]
				       *ST[W][X][P][Q]);
		}
	      }
	    }
	  }
	  
	}
      }/* end P */
    }
  }/* end M */
}

void tensors_aa_bb_dd_mm_plast (long ii,
				long ip,
				long nne,
				long ndn,
				long npres,
				double ai,
				double aj,
				double ak,
				double J,
				double *Psi,
				double L[3][3][3][3],
				double **Fn,
				double **Fr,
				double **Fr_I,
				double Jn,
				double Jr,
				double Tn,
				double Tr,
				double ****ST,
				double ****FF,
				double **f,
				double ***AA,
				double ***aa,
				double ***BB,
				double ***bb,
				double **DD,
				double **dd,
				double **MM,
				double **mm,
				double **S,
				double ****pFFp,
				double **pfp,
				double **UU,
				EPS *eps,
				double ***BB1,
				double ***bb1,
				double ***BB2,
				double ***bb2,
				double **MMp,
				double **mmp,
				double ****pGp)
{
  long Q,M,N,X,P,V,W,R;
  double ****GT,gT[3][3],help1;
  
  GT = aloc4 (3,3,ndn,nne);
  null_2d (gT);
  
  /* Tensors GT, gg and gT */
  tensors_GG (ii,ip,nne,ndn,ST,FF,f,UU,eps,GT,gT,pGp);
  
  /* Tensors BB1 and BB2 */
  for (Q=0;Q<ndn;Q++){
    for (M=0;M<npres;M++){
      for (N=0;N<nne;N++){
	AA[Q][M][N] = 0.0;
	BB[Q][M][N] = 0.0;
	BB1[Q][M][N] = 0.0;
	BB2[Q][M][N] = 0.0;
	
	help1 = 0.0;
	for (X=0;X<3;X++){
	  for (P=0;P<3;P++){
	    if (ST[X][P][Q][N] != 0.0)
	      help1 += Fr_I[P][X]*ST[X][P][Q][N];
	    
	    /* elastic */
	    if (pFFp[X][P][Q][N] != 0.0)
	      BB[Q][M][N] += (1./Jn*2./3.*pow(Tr,-1./3.)
			      *Psi[M]*pFFp[X][P][Q][N]*S[X][P]);
	    /* plastic */
	    if (GT[X][P][Q][N] != 0.0)
	      BB1[Q][M][N] += (1./Jn*pow(Tr,2./3.)
			       *Psi[M]*GT[X][P][Q][N]*S[X][P]);

	    if (pGp[X][P][Q][N] != 0.0)
	      BB2[Q][M][N] += (1./Jn*1./3.*pow(Tr,-1./3.)
			       *Psi[M]*pGp[X][P][Q][N]*S[X][P]);
	    
	    for (V=0;V<3;V++){
	      for (W=0;W<3;W++){
		if (L[X][P][V][W] == 0.0)
		  continue;

		/* elastic */
		if (pFFp[X][P][Q][N] != 0.0 && pfp[V][W] != 0.0)
		  BB[Q][M][N] += (1./Jn*1./3.*pow(Tr,1./3.)
				  *Psi[M]*pFFp[X][P][Q][N]
				  *L[X][P][V][W]*pfp[V][W]);

		/* plastic */
		if (gT[X][P] != 0.0 && pFFp[V][W][Q][N] != 0.0)
		  BB1[Q][M][N] += (1./Jn*1./2.*pow(Tr,4./3.)
				   *Psi[M]*gT[X][P]*L[X][P][V][W]
				   *pFFp[V][W][Q][N]);

		if (pGp[X][P][Q][N] != 0.0 && pfp[V][W] != 0.0)
		  BB2[Q][M][N] += (1./Jn*1./6.*pow(Tr,1./3.)
				   *Psi[M]*pGp[X][P][Q][N]
				   *L[X][P][V][W]*pfp[V][W]);
	      }
	    }
	    
	  }/* end P */
	}/* end X */
	
	AA[Q][M][N] = Jr*Psi[M]*help1;
	
	if (periodic == 1){
	  aa[Q][M][N]  += 1./(VVolume) *ai*aj*ak*J * AA[Q][M][N];
	  bb[Q][M][N]  += 1./(VVolume) *ai*aj*ak*J * BB[Q][M][N];
	  bb1[Q][M][N] += 1./(VVolume) *ai*aj*ak*J * BB1[Q][M][N];
	  bb2[Q][M][N] += 1./(VVolume) *ai*aj*ak*J * BB2[Q][M][N];
	}
	else{
	  aa[Q][M][N]  += ai*aj*ak*J * AA[Q][M][N];
	  bb[Q][M][N]  += ai*aj*ak*J * BB[Q][M][N];
	  bb1[Q][M][N] += ai*aj*ak*J * BB1[Q][M][N];
	  bb2[Q][M][N] += ai*aj*ak*J * BB2[Q][M][N];
	}
      }
    }
  }
  
  /* Tensor MM */
  for (Q=0;Q<npres;Q++){
    for (M=0;M<npres;M++){
      DD[Q][M] = 0.0;
      MM[Q][M] = 0.0;
      MMp[Q][M] = 0.0;
      
      for (N=0;N<3;N++){
	for (X=0;X<3;X++){
	  
	  /* elastic */
	  if (pfp[N][X] != 0.0)
	    MM[Q][M]  += (-1./Jn*1./9.*pow(Tr,-4./3.)
			  *Psi[Q]*Psi[M]*pfp[N][X]*S[N][X]);

	  /* plastic */
	  if (gT[N][X] != 0.0)
	    MMp[Q][M] += (1./Jn*1./3.*pow(Tr,-1./3.)
			  *Psi[Q]*Psi[M]*gT[N][X]*S[N][X]);
	  
	  for (P=0;P<3;P++){
	    for (R=0;R<3;R++){
	      if (L[N][X][P][R] == 0.0)
		continue;

	      /* elastic */
	      if (pfp[N][X] != 0.0 && pfp[P][R] != 0.0)
		MM[Q][M] += (1./Jn*1./9.*pow(Tr,-2./3.)
			     *Psi[Q]*Psi[M] * pfp[N][X]
			     *L[N][X][P][R]*pfp[P][R]);

	      /* plastic */
	      if (pfp[N][X] != 0.0 && gT[P][R] != 0.0)
		MMp[Q][M] += (1./Jn*1./6.*pow(Tr,1./3.)
			      *Psi[Q]*Psi[M] * pfp[N][X]
			      *L[N][X][P][R]*gT[P][R]);
	    }
	  }
	}
      }
      
      DD[Q][M]  = Tn/Jn*Psi[Q]*Psi[M];
      
      if (periodic == 1){
	dd[Q][M]  += 1./(VVolume) *ai*aj*ak*J * DD[Q][M];
	mm[Q][M]  += 1./(VVolume) *ai*aj*ak*J * MM[Q][M];
	mmp[Q][M] += 1./(VVolume) *ai*aj*ak*J * MMp[Q][M];
      }
      else{
	dd[Q][M]  += ai*aj*ak*J * DD[Q][M];
	mm[Q][M]  += ai*aj*ak*J * MM[Q][M];
	mmp[Q][M] += ai*aj*ak*J * MMp[Q][M];
      }
    }
  }
  dealoc4 (GT,3,3,ndn);
}
