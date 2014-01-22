#include "resid_on_elem.h"
#include <math.h>

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef CAST_MACROS_H
#include "cast_macros.h"
#endif

#ifndef DEF_GRAD_H
#include "def_grad.h"
#endif

#ifndef ELEM3D_H
#include "elem3d.h"
#endif

#ifndef INCL_H
#include "incl.h"
#endif

#ifndef INTEGRATION_H
#include "integration.h"
#endif

#ifndef PRESSU_SHAPE_H
#include "pressu_shape.h"
#endif

#ifndef RE1_RE2_RE3_H
#include "Re1_Re2_Re3.h"
#endif

#ifndef STRESS_STRAIN_H
#include "stress_strain.h"
#endif

#ifndef TENSORS_H
#include "tensors.h"
#endif

static const int periodic = 0;

int resid_on_elem (long ii,
		   long ndofn,
		   long nne,
		   long *nod,
		   ELEMENT *elem,
		   NODE *node,
		   MATGEOM matgeom,
		   HOMMAT *hommat,
		   double *x,
		   double *y,
		   double *z,
		   EPS *eps,
		   SIG *sig,
		   double *r_e,
		   long npres,
		   double nor_min,
		   double *fe,
		   CRPL *crpl,
		   double dt,
		   const int analysis)
{
  int err = 0;
  long i,j,k,II,JJ,KK,ip,M,N,P,R,U,W,X,ndn,mat;
  double ksi,eta,zet,*gk,*ge,*gz,*w,ai,aj,ak,J;
  double *N_x,*N_y,*N_z,*Psi,**Fn,**Fr,**Fr_I;
  double ****ST,Jn,Jr,Tn,Tr,pp,L[3][3][3][3],**E;
  double **Re1,*Re2,*Re3,****FF,**f,**re1,*re2;
  double *re3,***AA,***aa,***BB,***bb,**DD,***BB1;
  double ***bb1,***BB2,***bb2,**MMp,**mmp,**FnB;
  double ****pGp,dij,**dd,**MM,**mm,**S,**UU;
  double ****pFFp,**pfp;
  
  /* 3D */
  ndn = 3;

  mat = elem[ii].mat[2];
  
  gk = aloc1 (5);
  ge = aloc1 (5);
  gz = aloc1 (5);
  w = aloc1 (5);
  N_x = aloc1 (nne);
  N_y = aloc1 (nne);
  N_z = aloc1 (nne);
  FnB = aloc2 (3,3);
  Psi = aloc1 (npres);
  Fn = aloc2 (3,3);
  Fr = aloc2 (3,3);
  Fr_I = aloc2 (3,3);
  FF = aloc4 (3,3,ndn,nne);
  f = aloc2 (3,3);
  AA = aloc3 (ndn,npres,nne);
  aa = aloc3 (ndn,npres,nne);
  BB = aloc3 (ndn,npres,nne);
  bb = aloc3 (ndn,npres,nne);
  DD = aloc2 (npres,npres);
  dd = aloc2 (npres,npres);
  MM = aloc2 (npres,npres);
  mm = aloc2 (npres,npres);
  ST = aloc4 (3,3,ndn,nne);
  E = aloc2 (3,3);
  Re1 = aloc2 (ndn,nne);
  re1 = aloc2 (ndn,nne);
  Re2 = aloc1(npres);
  re2 = aloc1 (npres);
  Re3 = aloc1(npres);
  re3 = aloc1 (npres);
  S = aloc2 (3,3);

  
  if (analysis == FS_CRPL){
    UU = aloc2 (3,3);
    pFFp = aloc4 (3,3,ndn,nne);
    pfp = aloc2 (3,3);
    MMp = aloc2 (npres,npres);
    mmp = aloc2 (npres,npres);
    BB1 = aloc3 (ndn,npres,nne);
    bb1 = aloc3 (ndn,npres,nne);
    BB2 = aloc3 (ndn,npres,nne);
    bb2 = aloc3 (ndn,npres,nne);
    pGp = aloc4 (3,3,ndn,nne);
  }
  
  /* Integration */
  integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
  
  ip = 0;
  for (i=0;i<II;i++){
    for (j=0;j<JJ;j++){
      for (k=0;k<KK;k++){
	
	if (nne == 4)  {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	if (nne == 10) {
	  ksi = *(gk+k);
	  eta = *(ge+k);
	  zet = *(gz+k);
	  ai = *(w+k);
	  aj = 1.0;
	  ak = 1.0;
	}
	if (nne == 8)  {
	  ksi = *(gk+i);
	  eta = *(gk+j);
	  zet = *(gk+k);
	  ai = *(w+i);
	  aj = *(w+j);
	  ak = *(w+k);
	}
	
	/* Derivatives of shape functions */
	J = deriv (ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
	
	/****** FOR CRYSTAL PLASTICITY
		Fn denotes eFn || for multi-scale FLUCTUATION Fn *****/
	Fn[0][0] = eps[ii].il[ip].F[0];
	Fn[0][1] = eps[ii].il[ip].F[1];
	Fn[0][2] = eps[ii].il[ip].F[2];

	Fn[1][0] = eps[ii].il[ip].F[3];
	Fn[1][1] = eps[ii].il[ip].F[4];
	Fn[1][2] = eps[ii].il[ip].F[5];

	Fn[2][0] = eps[ii].il[ip].F[6];
	Fn[2][1] = eps[ii].il[ip].F[7];
	Fn[2][2] = eps[ii].il[ip].F[8];

	
	shape_tensor (nne,ndofn,N_x,N_y,N_z,ST);
	def_grad_get (nne,ndofn,CONST_4(double) ST,r_e,Fr); 
	def_grad_inv (Fr,Fr_I);
	Jr = def_grad_det (Fr);
	Jn = def_grad_det (Fn);
	
	/* Check determinants and cleanly exit with error if inverted
	   element */
	if(Jn <= 0.0 || Jr <= 0.0){
	  err += 1;
	  goto exit_function;
	}

	/* Pressure shape functions */
	pressu_shape (npres,ksi,eta,zet,Psi);
	
	Tn = 0.0;
	Tr = 0.0;
	pp = 0.0;
	for (M=0;M<npres;M++){
	  Tn += Psi[M]*eps[ii].T[M];
	  Tr += Psi[M]*eps[ii].d_T[M];
	  pp += Psi[M]*(sig[ii].p[M] + sig[ii].d_p[M]);
	}
	
	/* Set Fn-BAR */
	for (N=0;N<3;N++){
	  for (P=0;P<3;P++){
	    FnB[N][P] = pow(Tn,1./3.)*pow(Jn,-1./3.)*Fn[N][P];
	  }
	}
	
	/*** UNIT CELL APPROACH ***/
	if (periodic == 1) {
	  if (analysis == FS_CRPL){
	    S[0][0] = eps[ii].il[ip].Fp[0];
	    S[0][1] = eps[ii].il[ip].Fp[1];
	    S[0][2] = eps[ii].il[ip].Fp[2];

	    S[1][0] = eps[ii].il[ip].Fp[3];
	    S[1][1] = eps[ii].il[ip].Fp[4];
	    S[1][2] = eps[ii].il[ip].Fp[5];

	    S[2][0] = eps[ii].il[ip].Fp[6];
	    S[2][1] = eps[ii].il[ip].Fp[7];
	    S[2][2] = eps[ii].il[ip].Fp[8];
	    
	    def_grad_inv (S,FnB);
	  }/* analysis == FS_CRPL */
	  else{
	    for (N=0;N<3;N++){
	      for (P=0;P<3;P++){
		if (N == P) dij = 1.0;
		else dij = 0.0;
		FnB[N][P] = dij;
	      }
	    }
	  }
	  
	  for (N=0;N<3;N++){
	    for (P=0;P<3;P++){
	      Fr[N][P] += eps[0].F[N][P] + Fn[N][P];
	    }
	  }
	  
	  Jr = def_grad_det (Fr);
	  def_grad_inv (Fr,Fr_I);
	  Jn = Tn = 1.;
	  
	}/* end PERIODIC */
	
	/* Material stiffness matrix */
	matrix_tensor_3D (elem[ii].mat[2],hommat,L);
	
	/*** CRYSTAL PLASTICITY ***/
	if (analysis == FS_CRPL){
	  UU[0][0] = eps[ii].il[ip].UU[0];
	  UU[0][1] = eps[ii].il[ip].UU[1];
	  UU[0][2] = eps[ii].il[ip].UU[2];

	  UU[1][0] = eps[ii].il[ip].UU[3];
	  UU[1][1] = eps[ii].il[ip].UU[4];
	  UU[1][2] = eps[ii].il[ip].UU[5];

	  UU[2][0] = eps[ii].il[ip].UU[6];
	  UU[2][1] = eps[ii].il[ip].UU[7];
	  UU[2][2] = eps[ii].il[ip].UU[8];
	  
	  for (M=0;M<3;M++){
	    for (N=0;N<3;N++){
	      S[M][N] = 0.0;
	      for (P=0;P<3;P++){
		S[M][N] += FnB[M][P]*UU[P][N];
	      }
	    }
	  }
	  
	  /* Strain */
	  get_GL_strain (S,Fr,Jr,Tr,E);
	}/* analysis == FS_CRPL */
	else
	  /* Strain */
	  get_GL_strain (FnB,Fr,Jr,Tr,E);
	
	/* Stress */
	get_SPK_stress (L,E,S);
	
	/* Tensors used in analysis */
	tensors_FF_f (nne,ndn,FnB,Fr,Fr_I,Jr,ST,FF,f);
	
	/*** CRYSTAL PLASTICITY ***/
	if (analysis == FS_CRPL){
	  
	  for (P=0;P<3;P++){
	    for (R=0;R<3;R++){
	      /* eFn+1 and F* */
	      pfp[P][R] = 0.0;
	      
	      for (U=0;U<3;U++){
		for (W=0;W<3;W++){
		  pfp[P][R] += UU[U][P] * f[U][W] * UU[W][R];
		}
	      }
	      
	      for (M=0;M<ndn;M++){
		for (N=0;N<nne;N++){
		  pFFp[P][R][M][N] = 0.0;
		  
		  for (U=0;U<3;U++){
		    for (W=0;W<3;W++){
		      pFFp[P][R][M][N] += (UU[U][P] * FF[U][W][M][N]
					   * UU[W][R]);
		    }
		  }
		}
	      }
	    }/* end R < 3 */
	  }
	  
	  /* Elasto-Plastic tensors */
	  tensors_aa_bb_dd_mm_plast (ii,ip,nne,ndn,npres,ai,aj,ak,
				     J,Psi,L,Fn,Fr,Fr_I,Jn,Jr,Tn,Tr,
				     ST,FF,f,AA,aa,BB,bb,DD,dd,MM,
				     mm,S,pFFp,pfp,UU,eps,BB1,bb1,
				     BB2,bb2,MMp,mmp,pGp);
	}/* analysis == FS_CRPL */
	else
	  tensors_aa_bb_dd_mm (nne,ndn,npres,ai,aj,ak,J,Psi,L,FnB,
			       Fr,Fr_I,Jn,Jr,Tn,Tr,ST,FF,f,AA,aa,
			       BB,bb,DD,dd,MM,mm);
	
	/* Residuals */
	if (analysis == FS_CRPL){
	  Re1_Re2_Re3 (nne,ndn,npres,ai,aj,ak,J,Psi,Jr,Jn,Tr,Tn,
		       Fr_I,ST,pFFp,S,pfp,pp,Re1,re1,Re2,re2,Re3,re3);
	} else {
	  Re1_Re2_Re3 (nne,ndn,npres,ai,aj,ak,J,Psi,Jr,Jn,Tr,Tn,
		       Fr_I,ST,FF,S,f,pp,Re1,re1,Re2,re2,Re3,re3);
	}
	
	ip++;
      }/* end k */
    }/* end i */
  }/* end j */
  
  /****************************************************************/
  /****************************************************************/
  
  /*** CRYSTAL PLASTICITY ***/
  if (analysis == FS_CRPL){
    /* Tensors BB and bb */
    for (P=0;P<ndn;P++){
      for (M=0;M<npres;M++){
	for (N=0;N<nne;N++){
	  bb[P][M][N] += bb1[P][M][N];
	}
      }
    }
    for (P=0;P<npres;P++){
      for (M=0;M<npres;M++){
	mm[P][M] += mmp[P][M];
      }
    }
  }/* end analysis == FS_CRPL */
  
  for (X=0;X<nne;X++){
    for (N=0;N<ndn;N++){
      fe[X*ndofn+N] = re1[N][X];
      
      for (P=0;P<npres;P++){
	for (R=0;R<npres;R++){
	  
	  fe[X*ndofn+N] += (aa[N][P][X]*1./dd[P][R]*re3[R]
			    + bb[N][P][X]*1./dd[P][R]*re2[R]);
	  
	  for (U=0;U<npres;U++){
	    for (W=0;W<npres;W++){
	      fe[X*ndofn+N] += (aa[N][P][X]*1./dd[P][R]*mm[R][U]
				*1./dd[U][W]*re2[W]);
	    }
	  }
	}
      }
    }
  }/* end X */
  
 exit_function:
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
  dealoc1 (N_x);
  dealoc1 (N_y);
  dealoc1 (N_z);
  dealoc2 (FnB,3);
  dealoc1 (Psi);
  dealoc2 (Fn,3);
  dealoc2 (Fr,3);
  dealoc2 (Fr_I,3);
  dealoc4 (FF,3,3,ndn);
  dealoc2 (f,3);
  dealoc3 (AA,ndn,npres);
  dealoc3 (aa,ndn,npres);
  dealoc3 (BB,ndn,npres);
  dealoc3 (bb,ndn,npres);
  dealoc2 (DD,npres);
  dealoc2 (dd,npres);
  dealoc2 (MM,npres);
  dealoc2 (mm,npres);
  dealoc4 (ST,3,3,ndn);
  dealoc2 (E,3);
  dealoc2 (Re1,ndn);
  dealoc2 (re1,ndn);
  dealoc1 (Re2);
  dealoc1 (re2);
  dealoc1 (Re3);
  dealoc1 (re3);
  dealoc2 (S,3);
 
  if (analysis == FS_CRPL){
    dealoc2 (UU,3);
    dealoc4 (pFFp,3,3,ndn);
    dealoc2 (pfp,3);
    dealoc2 (MMp,npres);
    dealoc2 (mmp,npres);
    dealoc3 (BB1,ndn,npres);
    dealoc3 (bb1,ndn,npres);
    dealoc3 (BB2,ndn,npres);
    dealoc3 (bb2,ndn,npres);
    dealoc4 (pGp,3,3,ndn);
  }

  return err;
}
