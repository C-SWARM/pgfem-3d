#include "stiffmatel_fd.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
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

#ifndef PRESSU_SHAPE_H
#include "pressu_shape.h"
#endif

#ifndef STRESS_STRAIN_H
#include "stress_strain.h"
#endif

#ifndef TENSORS_H
#include "tensors.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

static const int periodic = 0;

/* VVolume only appears in deprecated 'periodic' branch. Poison value
   to return NaN, prompting investigation of how entered deprecated
   branch of code. */
static const double VVolume = 0.0;

int stiffmatel_fd (long ii,
		   long ndofn,
		   long nne,
		   long *nod,
		   double *x,
		   double *y,
		   double *z,
		   ELEMENT *elem,
		   MATGEOM matgeom,
		   HOMMAT *hommat,
		   NODE *node,
		   SIG *sig,
		   EPS *eps,
		   double *r_e,
		   long npres,
		   double nor_min,
		   double *Ks,
		   double dt,
		   CRPL *crpl,
		   long FNR,
		   double lm,
		   double *fe,
		   const int analysis)
{
  int err = 0;
  long  i,j,k,II,JJ,kK,ip,ndn,M,N,X,P,R,Q,U,V,W,Y,Z,mat,ndofe;
  double **B,**S,*N_x,*N_y,*N_z,J,ksi,eta,zet,*w,*gk,*ge,*gz,ai,
    aj,ak,****ST,**Fr,**Fr_I,**Fn,***AA,***aa,Jn,Jr,*Psi,**DD,
    **dd,L[3][3][3][3],**f,**MM,**mm;
  double ***BB,***bb,****FF,help1,help2,help3,****KK,****kk,Tn,
    Tr,pp,**E,****EE,****ee,****GG,****gg,****HH,****hh,help4,
    **UU,****pFFp,**pfp,***BB1,***bb1;
  double ***BB2,***bb2,**MMp,**mmp,****pKKp,****pkkp,****pGGp,
    ****pggp,****pGp,**FnB,dij,**Flam,**FFlam,**pFFplam,**KKlam,
    **EElam,**GGlam,**HHlam;
  double **pKKplam,**pGGplam,**pGplam,**kklam,**eelam,**gglam,
    **hhlam,**pkkplam,**pggplam,*Re2,*re2,*Re3,*re3,*Re2lam,
    *re2lam,*Re3lam,*re3lam;
  
  /*
    FILE *slon;
    if ((slon = fopen("../out/koza.dat","w")) == NULL ){
    PGFEM_printf("Output file is not possible to open\n");
    PGFEM_printf("Check the output file and run program again\n");
    }
  */
  
  /* 3D */
  ndn = 3;
  mat = elem[ii].mat[2];
  ndofe = ndofn*nne;
  
  gk = aloc1(5);
  ge = aloc1(5);
  gz = aloc1(5);
  w = aloc1 (5);
  N_x = aloc1 (nne);
  N_y = aloc1 (nne);
  N_z = aloc1 (nne);
  B = aloc2 (3,3);
  S = aloc2 (3,3);
  E = aloc2 (3,3);
  ST = aloc4 (3,3,ndn,nne);
  Fr = aloc2 (3,3);
  Fr_I = aloc2 (3,3);
  Fn = aloc2 (3,3);
  AA = aloc3 (ndn,npres,nne);
  aa = aloc3 (ndn,npres,nne);
  Psi = aloc1 (npres);
  DD = aloc2 (npres,npres);
  dd = aloc2 (npres,npres);
  f = aloc2 (3,3);
  MM = aloc2 (npres,npres);
  mm = aloc2 (npres,npres);
  BB = aloc3 (ndn,npres,nne);
  bb = aloc3 (ndn,npres,nne);
  FF = aloc4 (3,3,ndn,nne);
  KK = aloc4 (ndn,ndn,nne,nne);
  kk = aloc4 (ndn,ndn,nne,nne);
  EE = aloc4 (ndn,ndn,nne,nne);
  ee = aloc4 (ndn,ndn,nne,nne);
  GG = aloc4 (ndn,ndn,nne,nne);
  gg = aloc4 (ndn,ndn,nne,nne);
  HH = aloc4 (ndn,ndn,nne,nne);
  hh = aloc4 (ndn,ndn,nne,nne);
  FnB = aloc2 (3,3);
  
  if (periodic == 1 && (FNR == 2 || FNR == 3)){
    Flam = aloc2 (3,3);
    FFlam = aloc2 (3,3);
    KKlam = aloc2 (ndn,nne);
    EElam = aloc2 (ndn,nne);
    GGlam = aloc2 (ndn,nne);
    HHlam = aloc2 (ndn,nne);
    kklam = aloc2 (ndn,nne);
    eelam = aloc2 (ndn,nne);
    gglam = aloc2 (ndn,nne);
    hhlam = aloc2 (ndn,nne);
    Re2 = aloc1(npres);
    re2 = aloc1(npres);
    Re3 = aloc1(npres);
    re3 = aloc1(npres);
    Re2lam = aloc1(npres);
    re2lam = aloc1(npres);
    Re3lam = aloc1(npres);
    re3lam = aloc1(npres);
    
    if (analysis == FS_CRPL){
      pFFplam = aloc2 (3,3);
      pKKplam = aloc2 (ndn,nne);
      pGGplam = aloc2 (ndn,nne);
      pGplam = aloc2 (3,3);
      pkkplam = aloc2 (ndn,nne);
      pggplam = aloc2 (ndn,nne);
    }
  }
  
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
    pKKp = aloc4 (ndn,ndn,nne,nne);
    pkkp = aloc4 (ndn,ndn,nne,nne);
    pGGp = aloc4 (ndn,ndn,nne,nne);
    pggp = aloc4 (ndn,ndn,nne,nne);
    pGp = aloc4 (3,3,ndn,nne);
  }
  
  /* Integration */
  integrate (nne,&II,&JJ,&kK,gk,ge,gz,w);
  
  /********************************/
  /* FINITE DEFORMATION STIFFNESS */
  /********************************/
	
  ip = 0;
  for(i=0;i<II;i++){
    for(j=0;j<JJ;j++){
      for (k=0;k<kK;k++){
	
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
	
	/* Derivatives of shape functions and Jacobian of integration */
	J = deriv (ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
	
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

	/* Check determinants and exit cleanly with error if
	   negative */
	if(Jn <= 0.0 || Jr <= 0.0){
	  err += 1;
	  goto exit_function;
	}
	
	/* Pressure shape functions */
	pressu_shape (npres,ksi,eta,zet,Psi);
	
	Tn = 0.0; Tr = 0.0; pp = 0.0;
	for (M=0;M<npres;M++){
	  Tn += Psi[M]*eps[ii].T[M];
	  Tr += Psi[M]*eps[ii].d_T[M];
	  pp += Psi[M]*(sig[ii].p[M] + sig[ii].d_p[M]);
	}
	
	/* Set Fn-BAR */
	for (Q=0;Q<3;Q++){
	  for (R=0;R<3;R++){
	    FnB[Q][R] = pow(Tn,1./3.)*pow(Jn,-1./3.)*Fn[Q][R];
	  }
	}
	
	/********************** UNIT CELL APPROACH ************************/
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
	  }/* end analysis == FS_CRPL */
	  else{
	    for (N=0;N<3;N++){
	      for (P=0;P<3;P++){
		if (N == P)
		  dij = 1.0;
		else
		  dij = 0.0;
		FnB[N][P] = dij;
	      }
	    }
	  }
	  
	  for (N=0;N<3;N++){
	    for (P=0;P<3;P++)
	      Fr[N][P] += eps[0].F[N][P] + Fn[N][P];
	  }
	  
	  Jr = def_grad_det (Fr);
	  def_grad_inv (Fr,Fr_I);
	  Jn = Tn = 1.;
	  
	  if (FNR == 2 || FNR == 3){
	    /* Compression */
	    if (eps[0].type == 1){
	      Flam[0][0] = eps[0].load/((1.-lm*eps[0].load)*(1.-lm*eps[0].load));
	      Flam[1][1] = -eps[0].load;
	      Flam[2][2] = 0.0;
	    }
	  }

	}/* end PERIODIC */
	
	/* Material stiffness matrix */
	matrix_tensor_3D (elem[ii].mat[2],hommat,L);
	
	/*********************** CRYSTAL PLASTICITY *********************************/
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
	}/* end analysis == FS_CRPL */
	else
	  /* Strain */
	  get_GL_strain (FnB,Fr,Jr,Tr,E);
	
	/* Stress */
	get_SPK_stress (L,E,S);
	
	/* Tensors used in analysis */
	tensors_FF_f (nne,ndn,FnB,Fr,Fr_I,Jr,ST,FF,f);
	if (periodic == 1 && (FNR == 2 || FNR == 3))
	  tensors_FFla (ii,ip,nne,ndn,eps,FnB,Fr,Fr_I,Jr,f,
			UU,Flam,FFlam,pGplam,analysis);
	
	/*********************** CRYSTAL PLASTICITY ********************/
	if (analysis == FS_CRPL){
	  
	  for (P=0;P<3;P++){
	    for (R=0;R<3;R++){
	      /* eFn+1 and F* */
	      pfp[P][R] = 0.0;
	      if (periodic == 1 && (FNR == 2 || FNR == 3))
		pFFplam[P][R] = 0.0;
	      
	      for (U=0;U<3;U++){
		for (W=0;W<3;W++){
		  pfp[P][R] += UU[U][P] * f[U][W] * UU[W][R];
		  if (periodic == 1 && (FNR == 2 || FNR == 3))
		    pFFplam[P][R] += UU[U][P] * FFlam[U][W] * UU[W][R];
		}
	      }
	      
	      for (M=0;M<ndn;M++){
		for (N=0;N<nne;N++){
		  pFFp[P][R][M][N] = 0.0;
		  
		  for (U=0;U<3;U++){
		    for (W=0;W<3;W++){
		      if (FF[U][W][M][N] == 0.0)
			continue;
		      pFFp[P][R][M][N] += UU[U][P] * FF[U][W][M][N] * UU[W][R];
		    }
		  }
		}
	      }
	    }/* end R < 3 */
	  }/* end P < 3 */
	  
	  /* Elasto-Plastic tensors */
	  tensors_aa_bb_dd_mm_plast (ii,ip,nne,ndn,npres,ai,aj,ak,J,Psi,L,Fn,Fr,
				     Fr_I,Jn,Jr,Tn,Tr,ST,FF,f,AA,aa,BB,bb,DD,dd,
				     MM,mm,S,pFFp,pfp,UU,eps,BB1,bb1,BB2,bb2,MMp,mmp,pGp);
	}/* end analysis == FS_CRPL */
	else
	  tensors_aa_bb_dd_mm (nne,ndn,npres,ai,aj,ak,J,Psi,L,FnB,Fr,Fr_I,Jn,Jr,
			       Tn,Tr,ST,FF,f,AA,aa,BB,bb,DD,dd,MM,mm);
	
	/***************************************************************/
	/***************************************************************/
	/***************************************************************/
	
	/* Tensor EE */
	for (R=0;R<3;R++){
	  for (U=0;U<3;U++){
	    B[R][U] = 0.0;
	    if (analysis == FS_CRPL){
	      for (V=0;V<3;V++){
		for (W=0;W<3;W++){
		  
		  for (Y=0;Y<3;Y++){
		    for (Z=0;Z<3;Z++){
		      B[R][U] += 1./Jn*pow(Tr,2./3.)*pow(Jr,-2./3.) *
			FnB[R][V]*UU[V][W]*S[W][Y]*UU[Z][Y]*FnB[U][Z];
		    }
		  }
		}
	      }
	    }
	    else{
	      for (V=0;V<3;V++){
		for (W=0;W<3;W++){
		  B[R][U] += 1./Jn*pow(Tr,2./3.)*pow(Jr,-2./3.) *
		    FnB[R][V]*S[V][W]*FnB[U][W];
		}
	      }
	    }/* end else */
	  }/* end U */
	}/* end R */
	
	for (M=0;M<ndn;M++){
	  for (N=0;N<ndn;N++){
	    for (X=0;X<nne;X++){
	      for (P=0;P<nne;P++){
		/* null matrices */ 
		KK[M][N][X][P] = EE[M][N][X][P] = GG[M][N][X][P] = HH[M][N][X][P] = 0.0;
		if (analysis == FS_CRPL) { /* gboa */
		  pKKp[M][N][X][P] = pGGp[M][N][X][P] = 0.0;
		}  /* gboa */
		
		help1 = 0.0;
		help2 = 0.0;
		help3 = 0.0;
		help4 = 0.0;
		for (Q=0;Q<3;Q++){
		  for (R=0;R<3;R++){
		    for (U=0;U<3;U++){
		      for (V=0;V<3;V++){
			/* MATERIAL STIFFNESS MATRIX - Tensor KK */
			if (analysis == FS_CRPL) {
			  if (L[Q][R][U][V] != 0.0 &&
			      pFFp[Q][R][M][X] != 0.0 &&
			      pFFp[U][V][N][P] != 0.0)
			    KK[M][N][X][P] += 1./Jn*pow(Tr,4./3.)*pFFp[Q][R][M][X]*
			      L[Q][R][U][V]*pFFp[U][V][N][P];
			}
			else {
			  if (L[Q][R][U][V] != 0.0 &&
			      FF[Q][R][M][X] != 0.0 &&
			      FF[U][V][N][P] != 0.0)
			    KK[M][N][X][P] += 1./Jn*pow(Tr,4./3.)*FF[Q][R][M][X]*
			      L[Q][R][U][V]*FF[U][V][N][P];
			}
			
			if (analysis == FS_CRPL){/* PLASTIC MATERIAL + GEOMETRIC STIFFNESS */
			  if (L[Q][R][U][V] != 0.0 &&
			      pFFp[Q][R][M][X] != 0.0 &&
			      pGp[U][V][N][P] != 0.0)
			    pKKp[M][N][X][P] += 1./Jn*1./2.*pow(Tr,4./3.)* pFFp[Q][R][M][X]*
			      L[Q][R][U][V]*pGp[U][V][N][P];
			  
			  for (Y=0;Y<3;Y++){
			    if (ST[Q][R][N][P] == 0.0 || FF[U][V][M][X] == 0.0)
			      break;
			    for (Z=0;Z<3;Z++){
			      if (eps[ii].il[ip].dUU_Fr[Q][R][U][Y] != 0.0)
				pGGp[M][N][X][P] += 1./Jn*pow(Tr,2./3.)*ST[Q][R][N][P]*
				  FF[U][V][M][X]*eps[ii].il[ip].dUU_Fr[Q][R][U][Y]*UU[V][Z] *
				  S[Y][Z];
			      if (eps[ii].il[ip].dUU_Fr[Q][R][V][Z] != 0.0)
				pGGp[M][N][X][P] += 1./Jn*pow(Tr,2./3.)*ST[Q][R][N][P]*
				  FF[U][V][M][X]*UU[U][Y]*eps[ii].il[ip].dUU_Fr[Q][R][V][Z] *
				  S[Y][Z];
			    }
			  }/* end Y */
			}/* analysis == FS_CRPL */
			
			if (ST[Q][R][M][X] != 0.0 && ST[U][V][N][P] != 0.0){/* gboa */
			  help1 += (2./9.*Fr_I[R][Q]*Fr_I[V][U] + 1./3.*Fr_I[V][Q]*Fr_I[R][U])
			    *ST[Q][R][M][X]*ST[U][V][N][P];
			  help3 += (Fr_I[R][Q]*Fr_I[V][U] - Fr_I[V][Q]*Fr_I[R][U])*
			    ST[Q][R][M][X]*ST[U][V][N][P];
			}
			
		      }/* end V */
		    }/* end U */
		    if (ST[Q][R][M][X] != 0.0)
		      help2 += Fr_I[R][Q]*ST[Q][R][M][X];
		    if (ST[Q][R][N][P] != 0.0)
		      help4 += Fr_I[R][Q]*ST[Q][R][N][P];
		  }/* end R */
		}/* end Q */
		
		/* Geometric Stiffness */
		EE[M][N][X][P] = Jr*pp*help3;
		
		for (Q=0;Q<3;Q++){
		  for (R=0;R<3;R++){
		    for (U=0;U<3;U++){
		      GG[M][N][X][P] += B[Q][R] * 1./2.*
			(ST[U][Q][M][X]*ST[U][R][N][P] + ST[U][Q][N][P]*ST[U][R][M][X]);
		      HH[M][N][X][P] += B[Q][R] * 
			(help1*Fr[U][Q]*Fr[U][R] - 2./3.*help2 * 1./2.*
			 (ST[U][Q][N][P]*Fr[U][R] + Fr[U][Q]*ST[U][R][N][P])
			 - 2./3.*help4 * 1./2.*
			 (ST[U][Q][M][X]*Fr[U][R] + Fr[U][Q]*ST[U][R][M][X]));
		    }
		  }
		}/* end Q */
		if (periodic == 1){
		  ee[M][N][X][P] += 1./(VVolume) *ai*aj*ak*J * EE[M][N][X][P];
		  gg[M][N][X][P] += 1./(VVolume) *ai*aj*ak*J * GG[M][N][X][P];
		  hh[M][N][X][P] += 1./(VVolume) *ai*aj*ak*J * HH[M][N][X][P];
		  kk[M][N][X][P] += 1./(VVolume) *ai*aj*ak*J * KK[M][N][X][P];
		}
		else{
		  ee[M][N][X][P] += ai*aj*ak*J * EE[M][N][X][P];
		  gg[M][N][X][P] += ai*aj*ak*J * GG[M][N][X][P];
		  hh[M][N][X][P] += ai*aj*ak*J * HH[M][N][X][P];
		  kk[M][N][X][P] += ai*aj*ak*J * KK[M][N][X][P];  
		}
		if (analysis == FS_CRPL){
		  if (periodic == 1){
		    pkkp[M][N][X][P] += 1./(VVolume) *ai*aj*ak*J * pKKp[M][N][X][P];
		    pggp[M][N][X][P] += 1./(VVolume) *ai*aj*ak*J * pGGp[M][N][X][P];
		  }
		  else{
		    pkkp[M][N][X][P] += ai*aj*ak*J * pKKp[M][N][X][P];
		    pggp[M][N][X][P] += ai*aj*ak*J * pGGp[M][N][X][P];
		  }
		}

	      }/* end P */
	    }
	  }
	}/* end  M */
	
	/**************************/
	/* Tangential load vector */
	/**************************/
	
	if (periodic == 1 && (FNR == 2 || FNR == 3)){
	  for (M=0;M<ndn;M++){
	    for (X=0;X<nne;X++){
	      /* null matricies */
	      KKlam[M][X] = EElam[M][X] = GGlam[M][X] = HHlam[M][X] = 0.0;  /* gb */
	      if (analysis == FS_CRPL) {
		pKKplam[M][X] = pGGplam[M][X] = 0.0;
	      }  /* gb */
	      
	      help1 = 0.0;
	      help2 = 0.0;
	      help3 = 0.0;
	      help4 = 0.0;
	      for (Q=0;Q<3;Q++){
		for (R=0;R<3;R++){
		  for (U=0;U<3;U++){
		    for (V=0;V<3;V++){
		      if (analysis == FS_CRPL){/* plastic */
			if (L[Q][R][U][V] != 0.0 &&
			    pFFp[Q][R][M][X] != 0.0 &&
			    pFFplam[U][V] != 0.0){
			  KKlam[M][X] += 1./Jn*pow(Tr,4./3.)*pFFp[Q][R][M][X]*
			    L[Q][R][U][V]*pFFplam[U][V];
			}
		      }
		      else{/* elastic */
			if (L[Q][R][U][V] != 0.0 &&
			    FF[Q][R][M][X] != 0.0 &&
			    FFlam[U][V] != 0.0){
			  KKlam[M][X] += 1./Jn*pow(Tr,4./3.)*FF[Q][R][M][X]*
			    L[Q][R][U][V]*FFlam[U][V];
			}
		      }
		      
		      if (analysis == FS_CRPL){/* plastic */
			if (L[Q][R][U][V] != 0.0 &&
			    pFFp[Q][R][M][X] != 0.0 &&
			    pGplam[U][V] != 0.0){
			  pKKplam[M][X] += 1./Jn*1./2.*pow(Tr,4./3.)* pFFp[Q][R][M][X]*
			    L[Q][R][U][V]*pGplam[U][V];
			}
			
			for (Y=0;Y<3;Y++){
			  if (FF[U][V][M][X] == 0.0) break;
			  for (Z=0;Z<3;Z++){
			    if (eps[ii].il[ip].dUU_Fr[Q][R][U][Y] != 0.0)
			      pGGplam[M][X] += 1./Jn*pow(Tr,2./3.)*Flam[Q][R]*FF[U][V][M][X]*
				eps[ii].il[ip].dUU_Fr[Q][R][U][Y]*UU[V][Z] * S[Y][Z];
			    if (eps[ii].il[ip].dUU_Fr[Q][R][V][Z] != 0.0)
			      pGGplam[M][X] += 1./Jn*pow(Tr,2./3.)*Flam[Q][R]*FF[U][V][M][X]*
				UU[U][Y]*eps[ii].il[ip].dUU_Fr[Q][R][V][Z] * S[Y][Z];
			  }
			}/* end Y */
		      }/* analysis == FS_CRPL */
		      
		      if (ST[Q][R][M][X] != 0.0){
			help1 += (2./9.*Fr_I[R][Q]*Fr_I[V][U] + 1./3.*Fr_I[V][Q]*Fr_I[R][U])*
			  ST[Q][R][M][X]*Flam[U][V];
			help3 += (Fr_I[R][Q]*Fr_I[V][U] - Fr_I[V][Q]*Fr_I[R][U])*
			  ST[Q][R][M][X]*Flam[U][V];
		      }
		      
		    }/* end V */
		  }/* end U */
		  
		  if (ST[Q][R][M][X] != 0.0)
		    help2 += Fr_I[R][Q]*ST[Q][R][M][X];
		  help4 += Fr_I[R][Q]*Flam[Q][R];
		  
		}/* end R */
	      }/* end Q */
	      
	      /* Geometric Stiffness */
	      EElam[M][X] = Jr*pp*help3;
	      
	      for (Q=0;Q<3;Q++){
		for (R=0;R<3;R++){
		  for (U=0;U<3;U++){
		    GGlam[M][X] += B[Q][R] * 1./2.*(ST[U][Q][M][X]*Flam[U][R] + 
						    Flam[U][Q]*ST[U][R][M][X]);
		    HHlam[M][X] += B[Q][R] *
		      (help1*Fr[U][Q]*Fr[U][R] - 2./3.*help2 * 1./2.*
		       (Flam[U][Q]*Fr[U][R] + Fr[U][Q]*Flam[U][R])
		       - 2./3.*help4 * 1./2.*
		       (ST[U][Q][M][X]*Fr[U][R] + Fr[U][Q]*ST[U][R][M][X]));
		  }
		}
	      }/* end Q */
	      
	      eelam[M][X] += 1./(VVolume) *ai*aj*ak*J * EElam[M][X];
	      gglam[M][X] += 1./(VVolume) *ai*aj*ak*J * GGlam[M][X];
	      hhlam[M][X] += 1./(VVolume) *ai*aj*ak*J * HHlam[M][X];
	      kklam[M][X] += 1./(VVolume) *ai*aj*ak*J * KKlam[M][X];
	      
	      if (analysis == FS_CRPL){
		pkkplam[M][X] += 1./(VVolume) *ai*aj*ak*J * pKKplam[M][X];
		pggplam[M][X] += 1./(VVolume) *ai*aj*ak*J * pGGplam[M][X];
	      }

	    }/* end X */
	  }/* end  M */
	  
	  for (M=0;M<npres;M++){
	    Re2lam[M] = Re3lam[M] = 0.0; 
	    
	    Re2[M] = 1./Jn*Psi[M]*(Jr*Jn - Tr*Tn);
	    Re3[M] = -1./Jn*Tn*pp*Psi[M];
	    
	    for (X=0;X<3;X++){
	      for (P=0;P<3;P++){
		
		if (analysis == FS_CRPL)
		  Re3[M] += 1./Jn*1./3.*pow(Tr,-1./3)*Psi[M]*pfp[X][P]*S[X][P];
		else               
		  Re3[M] += 1./Jn*1./3.*pow(Tr,-1./3)*Psi[M]*f[X][P]*S[X][P];
		
		Re2lam[M] += Jr*Psi[M]*Fr_I[P][X]*Flam[X][P];
		
		if (analysis == FS_CRPL){/* plastic */
		  if (pFFplam[X][P] != 0.0) 
		    Re3lam[M] += 1./Jn*2./3.*pow(Tr,-1./3.)*Psi[M]*pFFplam[X][P]*S[X][P];
		  if (pGplam[X][P]  != 0.0) 
		    Re3lam[M] += 1./Jn*1./3.*pow(Tr,-1./3.)*Psi[M]*pGplam[X][P]*S[X][P];
		}
		else{/* elastic */
		  if (FFlam[X][P] != 0.0)
		    Re3lam[M] += 1./Jn*2./3.*pow(Tr,-1./3.)*Psi[M]*FFlam[X][P]*S[X][P];
		}
		
		for (V=0;V<3;V++){
		  for (W=0;W<3;W++){
		    if (L[X][P][V][W] == 0.0) continue;
		    if (analysis == FS_CRPL){/* plastic */
		      if (pFFplam[X][P] != 0.0 && pfp[V][W] != 0.0)
			Re3lam[M] += 1./Jn*1./3.*pow(Tr,1./3.)*Psi[M]*pFFplam[X][P]*
			  L[X][P][V][W]*pfp[V][W];
		      if (pGplam[X][P] != 0.0 && pfp[V][W] != 0.0) 
			Re3lam[M] += 1./Jn*1./6.*pow(Tr,1./3.)*Psi[M]*pGplam[X][P]*
			  L[X][P][V][W]*pfp[V][W];  
		    }
		    else{/* elastic */
		      if (FFlam[X][P] != 0.0 && f[V][W] != 0.0)
			Re3lam[M] += 1./Jn*1./3.*pow(Tr,1./3.)*Psi[M]*FFlam[X][P]*
			  L[X][P][V][W]*f[V][W];
		    }
		  }
		}/* end V */
		
	      }/* end P */
	    }/* end X */
	    
	    re2[M] += 1./(VVolume) *ai*aj*ak*J* Re2[M];
	    re3[M] += 1./(VVolume) *ai*aj*ak*J* Re3[M];
	    
	    re2lam[M] += 1./(VVolume) *ai*aj*ak*J* Re2lam[M];
	    re3lam[M] += 1./(VVolume) *ai*aj*ak*J* Re3lam[M];
	    
	  }/* end M */
	}/* end periodic */
	
	ip++;
      }/* end k */
    }/* end j */
  }/* end i */
  
  /*************************************/
  /* Element Stiffness Matrix Assembly */
  /*************************************/
  
  /******************** CRYSTAL PLASTICITY *********************/
  if (analysis == FS_CRPL){
    /* Tensors BB and bb */
    for (P=0;P<ndn;P++){
      for (M=0;M<npres;M++){
	for (N=0;N<nne;N++){
	  bb1[P][M][N] += bb[P][M][N];
	  bb2[P][M][N] += bb[P][M][N];
	}
      }
    }
    for (P=0;P<npres;P++){
      for (M=0;M<npres;M++){
	mm[P][M] += mmp[P][M];
      }
    }
  }

  /* Stiffness matrix is -- gboa -- */
  for (M=0;M<ndn;M++){
    for (N=0;N<ndn;N++){
      for (X=0;X<nne;X++){
	for (P=0;P<nne;P++){
	  
	  if (analysis == FS_CRPL)
	    KK[M][N][X][P] = kk[M][N][X][P] + gg[M][N][X][P] + hh[M][N][X][P] +
	      ee[M][N][X][P] + pkkp[M][N][X][P] + pggp[M][N][X][P];
	  else              
	    KK[M][N][X][P] = kk[M][N][X][P] + gg[M][N][X][P] +
	      hh[M][N][X][P] + ee[M][N][X][P];
	  
	  for (R=0;R<npres;R++){
	    for (U=0;U<npres;U++){

	      if (analysis == FS_CRPL)
		KK[M][N][X][P] += aa[M][R][X]*(1./dd[R][U])*bb2[N][U][P] + bb1[M][R][X]*
		  (1./dd[R][U])*aa[N][U][P];
	      else      
		KK[M][N][X][P] += aa[M][R][X]*(1./dd[R][U])*bb[N][U][P] + bb[M][R][X]*
		  (1./dd[R][U])*aa[N][U][P];
	      
	      for (V=0;V<npres;V++){
		for (W=0;W<npres;W++){
		  KK[M][N][X][P] += aa[M][R][X]*(1./dd[R][U])*mm[U][V]*
		    (1./dd[V][W])*aa[N][W][P];
		}
	      }
	    }
	  }
	  
	  /* Composition */
	  Ks[(X*ndn+M)*ndofe + P*ndn+N] = KK[M][N][X][P];
	  
	}/* end P */
      }/* end X */
    }/* end N */
  }/* end M */
  
  /**************************/
  /* Tangential load vector */
  /**************************/
  
  if (periodic == 1 && (FNR == 2 || FNR == 3)){/* fext(l) = -dfint/dl */
    for (X=0;X<nne;X++){
      for (M=0;M<ndn;M++){
	if (analysis == FS_CRPL)
	  fe[X*ndofn+M] = -1.*(kklam[M][X] + gglam[M][X] + hhlam[M][X] +
			       eelam[M][X] + pkkplam[M][X] + pggplam[M][X]);
	else              
	  fe[X*ndofn+M] = -1.*(kklam[M][X] + gglam[M][X] + hhlam[M][X] + eelam[M][X]);
	
	for (P=0;P<npres;P++){
	  for (R=0;R<npres;R++){
	    
	    if (analysis == FS_CRPL)
	      fe[X*ndofn+M] += -1.*(aa[M][P][X]*1./dd[P][R]*re3lam[R] +
				    bb1[M][P][X]*1./dd[P][R]*re2lam[R]);
	    else      
	      fe[X*ndofn+M] += -1.*(aa[M][P][X]*1./dd[P][R]*re3lam[R] +
				    bb[M][P][X]*1./dd[P][R]*re2lam[R]);

	    for (U=0;U<npres;U++){
	      for (W=0;W<npres;W++){
		fe[X*ndofn+M] += -1.*(aa[M][P][X]*1./dd[P][R]*mm[R][U]*1./dd[U][W]*re2lam[W]);
	      }
	    }
	  }
	}/* end P */
	
      }
    }/* end X */
  }/* end periodic */
  
  /*********************************************************/
  /* PRINT OUT OF MATERIAL AND GEOMETRIC STIFFNESS MATRICES */
  /*********************************************************/

  /*
    PGFEM_fprintf (slon,"\n");
    PGFEM_fprintf (slon,"\n");
    
    PGFEM_fprintf (slon,"Stiffness from Hu-Washizu\n");
    PGFEM_fprintf (slon,"\n");
    
    for (M=0;M<ndofe;M++){
    for (N=0;N<ndofe;N++){
    PGFEM_fprintf (slon,"%4.5f  ",Ku[M][N]);
    }
    PGFEM_fprintf (slon,"\n");
    }
    PGFEM_fprintf (slon,"\n");
    PGFEM_fprintf (slon,"\n");
    
    PGFEM_fprintf (slon,"Material Stiffness\n");
    PGFEM_fprintf (slon,"\n");
    
    for (X=0;X<nne;X++){
    for (M=0;M<ndn;M++){
    for (P=0;P<nne;P++){
    for (N=0;N<ndn;N++){
    PGFEM_fprintf (slon,"%4.5f  ",kk[M][N][X][P]);
    }
    }
    PGFEM_fprintf (slon,"\n");
    }
    }
    
    PGFEM_fprintf (slon,"\n");
    PGFEM_fprintf (slon,"\n");
    
    PGFEM_fprintf (slon,"Geometric Stiffness\n");
    PGFEM_fprintf (slon,"\n");
    
    for (X=0;X<nne;X++){
    for (M=0;M<ndn;M++){
    for (P=0;P<nne;P++){
    for (N=0;N<ndn;N++){
    PGFEM_fprintf (slon,"%4.5f  ",gg[M][N][X][P]);
    }
    }
    PGFEM_fprintf (slon,"\n");
    }
    }
    
    i = fflush (slon);
    
    fclose (slon);
  */
  
 exit_function:
  dealoc1(gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1(w);
  dealoc1 (N_x);
  dealoc1 (N_y);
  dealoc1 (N_z);
  dealoc2 (B,3);
  dealoc2 (S,3);
  dealoc2 (E,3);
  dealoc4 (ST,3,3,ndn);
  dealoc2 (Fr,3);
  dealoc2 (Fr_I,3);
  dealoc2 (Fn,3);
  dealoc3 (AA,ndn,npres);
  dealoc3 (aa,ndn,npres);
  dealoc1 (Psi);
  dealoc2 (DD,npres);
  dealoc2 (dd,npres);
  dealoc2 (f,3);
  dealoc2 (MM,npres);
  dealoc2 (mm,npres);
  dealoc3 (BB,ndn,npres);
  dealoc3 (bb,ndn,npres);
  dealoc4 (FF,3,3,ndn);
  dealoc4 (KK,ndn,ndn,nne);
  dealoc4 (kk,ndn,ndn,nne);
  dealoc4 (EE,ndn,ndn,nne);
  dealoc4 (ee,ndn,ndn,nne);
  dealoc4 (GG,ndn,ndn,nne);
  dealoc4 (gg,ndn,ndn,nne);
  dealoc4 (HH,ndn,ndn,nne);
  dealoc4 (hh,ndn,ndn,nne);
  dealoc2 (FnB,3);
  
  if (periodic == 1 && (FNR == 2 || FNR == 3)){
    dealoc2 (Flam,3);
    dealoc2 (FFlam,3);
    dealoc2 (KKlam,ndn);
    dealoc2 (EElam,ndn);
    dealoc2 (GGlam,ndn);
    dealoc2 (HHlam,ndn);
    dealoc2 (kklam,ndn);
    dealoc2 (eelam,ndn);
    dealoc2 (gglam,ndn);
    dealoc2 (hhlam,ndn);
    dealoc1 (Re2);
    dealoc1 (re2);
    dealoc1 (Re3);
    dealoc1 (re3);
    dealoc1 (Re2lam);
    dealoc1 (re2lam);
    dealoc1 (Re3lam);
    dealoc1 (re3lam);

    if (analysis == FS_CRPL){
      dealoc2 (pFFplam,3);
      dealoc2 (pKKplam,ndn);
      dealoc2 (pGGplam,ndn);
      dealoc2 (pGplam,3);
      dealoc2 (pkkplam,ndn);
      dealoc2 (pggplam,ndn);
    }
  }/* end periodic */
  
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
    dealoc4 (pKKp,ndn,ndn,nne);
    dealoc4 (pkkp,ndn,ndn,nne);
    dealoc4 (pGGp,ndn,ndn,nne);
    dealoc4 (pggp,ndn,ndn,nne);
    dealoc4 (pGp,3,3,ndn);
  }

  return err;
}
