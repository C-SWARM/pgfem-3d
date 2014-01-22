#include "fd_increment.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef INCL_H
#include "incl.h"
#endif

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef GET_DOF_IDS_ON_ELEM_H
#include "get_dof_ids_on_elem.h"
#endif

#ifndef GET_NDOF_ON_ELEM_H
#include "get_ndof_on_elem.h"
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

#ifndef TA_GA_H
#include "TA_GA.h"
#endif

#ifndef TENSORS_H
#include "tensors.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

static const int periodic = 0;

void fd_increment (long ne,
		     long nn,
		     long ndofn,
		     long npres,
		     MATGEOM matgeom,
		     HOMMAT *hommat,
		     ELEMENT *elem,
		     NODE *node,
		     SUPP sup,
		     EPS *eps,
		     SIG *sig,
		     double *d_r,
		     double *r,
		     double nor_min,
		     CRPL *crpl,
		     double dt,
		     long nce,
		     COEL *coel,
		     double *pores,
		     MPI_Comm mpi_comm,
		     const double VVolume,
		     const PGFem3D_opt *opts)
{
  long ii,i,j,k,ip,II,JJ,KK,nne,*nod;
  long M,N,P,Q,R,ndn,mat,nss,U,W,ndofe,ndofc,*cn;
  double *x,*y,*z,*gk,*ge,*gz,*w,ksi,eta,zet,ai,aj,ak;
  double **Fn,**Fr,*N_x,*N_y,*N_z,Jn,Jr,Tr,Tn,****ST,**E;
  double J,L[3][3][3][3],*r_e,*Psi,**S,**UU,**UU_I,**Fp,*TA;
  double *GA,dij,**Fn1,**FnB,PL,**AA,**Fn1B,**BB,*X,*Y;
  double pom,volume,EL_e,**CC,**FoN,**FlN,**DD,**Fee,**Fpp;
  double GEL_e,GVol,*GPf,LPf[63],*GDp,LDp[9];
  double Gpores;

  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);
  
  /* 3D */ ndn = 3;
  
  gk = aloc1 (5);
  ge = aloc1 (5);
  gz = aloc1 (5);
  w = aloc1 (5);
  Psi = aloc1 (npres);
  Fn = aloc2 (3,3);
  Fr = aloc2 (3,3);
  S = aloc2 (3,3);
  E = aloc2 (3,3);
  Fn1 = aloc2 (3,3);
  UU = aloc2 (3,3);
  UU_I = aloc2 (3,3);
  Fp = aloc2 (3,3);
  FnB = aloc2 (3,3);
  AA = aloc2 (3,3);
  BB = aloc2 (3,3);
  Fn1B = aloc2 (3,3);
  X = aloc1 (3);
  Y = aloc1 (3);
  CC = aloc2 (3,3);
  FoN = aloc2 (3,3);
  FlN = aloc2 (3,3);
  DD = aloc2 (3,3);
  Fee = aloc2 (3,3);
  Fpp = aloc2 (3,3);

  
  /***** UNIT CELL APPROACH *****/
  if (periodic == 1){
    EL_e = 0.0;
    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
	FoN[P][R] = eps[0].FB[P][R];
	eps[0].P[P][R] = eps[0].FB[P][R] = eps[0].Fe[P][R] = 0.0;

	if (opts->analysis_type == FS_CRPL){
	  eps[0].Fp[P][R] = 0.0;
	}
      }
    }
  }/* end PERIODIC */
  
  /* Deformation rate */
  for (P=0;P<3;P++){ 
    for (R=0;R<3;R++){ 
      eps[0].Dp[P][R] = 0.0;
    }
  }
  
  for (ii=0;ii<ne;ii++){
    
    mat = elem[ii].mat[2];
    if (opts->analysis_type == FS_CRPL) 
      nss = crpl[mat].nss;
    else
      nss = 1;
    
    /* Number of element nodes */
    nne = elem[ii].toe;
    /* Nodes on element */
    nod = aloc1l (nne);
    elemnodes (ii,nne,nod,elem);
    /* Element Dof */
    ndofe = get_ndof_on_elem_nodes(nne,nod,node);
    
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);
    N_x = aloc1 (nne);
    N_y = aloc1 (nne);
    N_z = aloc1 (nne);
    ST = aloc4 (3,3,ndn,nne);
    r_e = aloc1 (ndofe);
    TA = aloc1 (nss);
    GA = aloc1 (nss);
    cn = aloc1l (ndofe);
    
    /* Coordinates of element */
    switch(opts->analysis_type){
    case DISP:
      nodecoord_total (nne,nod,node,x,y,z);
      break;
    default:
      nodecoord_updated (nne,nod,node,x,y,z);
      break;
    }
    
    /* code numbers on element */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);
    
    /* deformation on element */
    def_elem (cn,ndofe,d_r,elem,node,r_e,sup,0);
    
    /* Volume of element */
    if (nne == 4)   volume = Tetra_V (x,y,z);
    if (nne == 8)   volume = Hexa_V (x,y,z);
    if (nne == 10)  volume = Tetra_qv_V (nne,ndofn,x,y,z);
    
    /* Integration */
    integrate (nne,&II,&JJ,&KK,gk,ge,gz,w);
    
    /* Null fields */
    for (k=0;k<6;k++){
      eps[ii].el.o[k] = 0.0;
      sig[ii].el.o[k] = 0.0;
      if (opts->analysis_type == FS_CRPL) {
	eps[ii].pl.o[k] = 0.0;
	eps[ii].pl.eq[0] = 0.0;
      }
    }
    eps[ii].GD = 0.0;
    
    ip = 0;
    for (i=0;i<II;i++){
      for (j=0;j<JJ;j++){
	for (k=0;k<KK;k++){
	  
	  if (nne == 4)  {ksi = *(gk+k);
	    eta = *(ge+k);
	    zet = *(gz+k);
	    ai = *(w+k);
	    aj = 1.0;
	    ak = 1.0;
	  }
	  if (nne == 10) {ksi = *(gk+k);
	    eta = *(ge+k);
	    zet = *(gz+k);
	    ai = *(w+k);
	    aj = 1.0;
	    ak = 1.0;
	  }
	  if (nne == 8)  {ksi = *(gk+i);
	    eta = *(gk+j);
	    zet = *(gk+k);
	    ai = *(w+i);
	    aj = *(w+j);
	    ak = *(w+k);
	  }
	  
	  /* Derivatives of shape functions and Jacobian of integration */
	  J = deriv (ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
	  
	  /* eFn -> denotes only elastic part of deformation fpr plasticity */
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
	  Jr = def_grad_det (Fr);
	  Jn = def_grad_det (Fn);
	  
	  /* Pressure shape functions */
	  pressu_shape (npres,ksi,eta,zet,Psi);
	  
	  Tr = 0.0;
	  Tn = 0,0; 
	  for (M=0;M<npres;M++){
	    Tn += Psi[M]*eps[ii].T[M];
	    Tr += Psi[M]*eps[ii].d_T[M];
	  }
	  
	  /* Fn+1 for elasticity and F* = Fr*eFn for plasticity || SET Fn-BAR */
	  for (M=0;M<3;M++){
	    for (N=0;N<3;N++){
	      Fn1[M][N] = 0.0;
	      FnB[M][N] = pow(Tn,1./3.)*pow(Jn,-1./3.)*Fn[M][N];

	      for (P=0;P<3;P++){
		Fn1[M][N] += Fr[M][P]*Fn[P][N];
	      }
	    }
	  }
	  
	  /**** UNIT CELL APPROACH ****/
	  if (periodic == 1) {
	    if (opts->analysis_type == FS_CRPL){
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
	    }
	    else{
	      for (N=0;N<3;N++){
		for (P=0;P<3;P++){
		  if (N == P) dij = 1.0; else dij = 0.0;
		  FnB[N][P] = dij;
		}
	      }
	    }/* end elastic */
	    
	    /* Fn1 : total fluctuation def gradient || Fr total micro
	       deformation gradient */
	    for (N=0;N<3;N++){
	      for (P=0;P<3;P++){
		Fn1[N][P] = Fr[N][P] + Fn[N][P];
		Fr[N][P] += eps[0].F[N][P] + Fn[N][P];
	      }
	    }
	    
	    Jr = def_grad_det (Fr); 
	    Jn = Tn = 1.;
	    
	  }/* end PERIODIC */
	  
	  /**** CRYSTAL PLASTICITY ****/
	  if (opts->analysis_type == FS_CRPL){
	    S[0][0] = eps[ii].il[ip].Fp[0];
	    S[0][1] = eps[ii].il[ip].Fp[1];
	    S[0][2] = eps[ii].il[ip].Fp[2];

	    S[1][0] = eps[ii].il[ip].Fp[3];
	    S[1][1] = eps[ii].il[ip].Fp[4];
	    S[1][2] = eps[ii].il[ip].Fp[5];

	    S[2][0] = eps[ii].il[ip].Fp[6];
	    S[2][1] = eps[ii].il[ip].Fp[7];
	    S[2][2] = eps[ii].il[ip].Fp[8];

	    
	    UU[0][0] = eps[ii].il[ip].UU[0];
	    UU[0][1] = eps[ii].il[ip].UU[1];
	    UU[0][2] = eps[ii].il[ip].UU[2];

	    UU[1][0] = eps[ii].il[ip].UU[3];
	    UU[1][1] = eps[ii].il[ip].UU[4];
	    UU[1][2] = eps[ii].il[ip].UU[5];

	    UU[2][0] = eps[ii].il[ip].UU[6];
	    UU[2][1] = eps[ii].il[ip].UU[7];
	    UU[2][2] = eps[ii].il[ip].UU[8];

	    
	    eps[ii].il[ip].UU[0] = 1.0;
	    eps[ii].il[ip].UU[1] = 0.0;
	    eps[ii].il[ip].UU[2] = 0.0;

	    eps[ii].il[ip].UU[3] = 0.0;
	    eps[ii].il[ip].UU[4] = 1.0;
	    eps[ii].il[ip].UU[5] = 0.0;

	    eps[ii].il[ip].UU[6] = 0.0;
	    eps[ii].il[ip].UU[7] = 0.0;
	    eps[ii].il[ip].UU[8] = 1.0;
	    
	    def_grad_inv (UU,UU_I);
	    
	    /********************************************************/
	    /* EFFECTIVE PLASTIC STRAIN */
	    
	    for (M=0;M<3;M++){
	      for (N=0;N<3;N++){
		if (M == N)
		  dij = 1.0;
		else
		  dij = 0.0;
		BB[M][N] = 1./dt*(dij - UU[M][N]);
	      }
	    }
	    
	    for (M=0;M<3;M++){
	      for (N=0;N<3;N++){
		AA[M][N] = 1./2.*(BB[M][N] + BB[N][M]);
		/* Macroscopic Dp */
		eps[0].Dp[M][N]  += 1./VVolume *ai*aj*ak*J* AA[M][N];
	      }
	    }
	    
	    pom = 0;
	    for (M=0;M<3;M++){
	      for (N=0;N<3;N++){
		pom += AA[M][N]*AA[M][N];
	      }
	    }
	    
	    eps[ii].il[ip].eff += dt*sqrt(2./3.*pom);
	    eps[ii].pl.eq[0] += ai*aj*ak*J/volume*eps[ii].il[ip].eff;
	    /********************************************************/	    
	    
	    for (M=0;M<3;M++){
	      for (N=0;N<3;N++){
		Fp[M][N] = 0.0;
		for (P=0;P<3;P++){
		  Fp[M][N] += UU_I[M][P]*S[P][N];
		}
	      }
	    }

	    /* pFn */
	    AA[0][0] = eps[ii].il[ip].Fp[0];
	    AA[0][1] = eps[ii].il[ip].Fp[1];
	    AA[0][2] = eps[ii].il[ip].Fp[2];

	    AA[1][0] = eps[ii].il[ip].Fp[3];
	    AA[1][1] = eps[ii].il[ip].Fp[4];
	    AA[1][2] = eps[ii].il[ip].Fp[5];

	    AA[2][0] = eps[ii].il[ip].Fp[6];
	    AA[2][1] = eps[ii].il[ip].Fp[7];
	    AA[2][2] = eps[ii].il[ip].Fp[8];

	    
	    /* pFn+1 || Avaluate Plastic deformation gradient at time
	       t+1 */
	    eps[ii].il[ip].Fp[0] = Fp[0][0];
	    eps[ii].il[ip].Fp[1] = Fp[0][1];
	    eps[ii].il[ip].Fp[2] = Fp[0][2];

	    eps[ii].il[ip].Fp[3] = Fp[1][0];
	    eps[ii].il[ip].Fp[4] = Fp[1][1];
	    eps[ii].il[ip].Fp[5] = Fp[1][2];

	    eps[ii].il[ip].Fp[6] = Fp[2][0];
	    eps[ii].il[ip].Fp[7] = Fp[2][1];
	    eps[ii].il[ip].Fp[8] = Fp[2][2];

	    
	    /* eFn+1 || Avaluate Elastic deformation gradient at time t+1 */
	    if (periodic != 1){
	      for (M=0;M<3;M++){
		for (N=0;N<3;N++){
		  S[M][N] = 0.0;
		  for (P=0;P<3;P++){
		    S[M][N] += Fn1[M][P]*UU[P][N];
		  }
		}
	      }
	      for (M=0;M<3;M++){
		for (N=0;N<3;N++){
		  Fn1[M][N] = S[M][N];
		}
	      }
	    }
	  }/* end analysis == FS_CRPL */
	  
	  /* Deformation gradient || Fn+1 or eFn+1 for crystal
	     plasticity || F(fluctuation) for multi-scale modeling */
	  eps[ii].il[ip].F[0] = Fn1[0][0];
	  eps[ii].il[ip].F[1] = Fn1[0][1];
	  eps[ii].il[ip].F[2] = Fn1[0][2];

	  eps[ii].il[ip].F[3] = Fn1[1][0];
	  eps[ii].il[ip].F[4] = Fn1[1][1];
	  eps[ii].il[ip].F[5] = Fn1[1][2];

	  eps[ii].il[ip].F[6] = Fn1[2][0];
	  eps[ii].il[ip].F[7] = Fn1[2][1];
	  eps[ii].il[ip].F[8] = Fn1[2][2];

	  
	  /**** UNIT CELL APPROACH ****/
	  
	  /* eFn+1-BAR */
	  if (opts->analysis_type == FS_CRPL){
	    for (M=0;M<3;M++){
	      for (N=0;N<3;N++){
		Fn1B[M][N] = 0.0;
		for (P=0;P<3;P++){
		  for (Q=0;Q<3;Q++){
		    Fn1B[M][N] +=  (pow(Tr,1./3.)*pow(Jr,-1./3.)
				    *Fr[M][P]*FnB[P][Q]*UU[Q][N]);
		  }
		}
	      }
	    }
	  }
	  else{
	    for (M=0;M<3;M++){
	      for (N=0;N<3;N++){
		Fn1B[M][N] = 0.0;
		for (P=0;P<3;P++){
		  Fn1B[M][N] +=  (pow(Tr,1./3.)*pow(Jr,-1./3.)
				  *Fr[M][P]*FnB[P][N]);
		}
	      }
	    }
	  }
	  
	  if (periodic == 1){
	    
	    /* Elastic : eFn at time n */
	    BB[0][0] = eps[ii].il[ip].Fe[0];
	    BB[0][1] = eps[ii].il[ip].Fe[1];
	    BB[0][2] = eps[ii].il[ip].Fe[2];

	    BB[1][0] = eps[ii].il[ip].Fe[3];
	    BB[1][1] = eps[ii].il[ip].Fe[4];
	    BB[1][2] = eps[ii].il[ip].Fe[5];

	    BB[2][0] = eps[ii].il[ip].Fe[6];
	    BB[2][1] = eps[ii].il[ip].Fe[7];
	    BB[2][2] = eps[ii].il[ip].Fe[8];

	    
	    if (opts->analysis_type == FS_CRPL){
	      /* Total Fn */
	      for (P=0;P<3;P++){
		for (R=0;R<3;R++){
		  FlN[P][R] = 0.0;
		  for (U=0;U<3;U++){
		    FlN[P][R] += BB[P][U]*AA[U][R];
		  }
		}
	      }
	    }
	    else{
	      
	      /* Inverse total def. grad. */
	      def_grad_inv (Fn1B,AA);
	      
	      for (P=0;P<3;P++){
		for (R=0;R<3;R++){
		  FlN[P][R] = BB[P][R]; E[P][R] = 0.0;
		  for (U=0;U<3;U++){
		    E[P][R] += BB[P][U]*AA[U][R];
		  }
		}
	      }
	      
	      /* EFFECTIVE STRAIN */
	      for (M=0;M<3;M++){
		for (N=0;N<3;N++){
		  if (M == N) dij = 1.0; else dij = 0.0;
		  BB[M][N] = 1./dt*(dij - E[M][N]);
		}
	      }
	      
	      for (M=0;M<3;M++){
		for (N=0;N<3;N++){
		  /* Macroscopic Dp */
		  eps[0].Dp[M][N]  += (1./VVolume *ai*aj*ak*J
				       * 1./2.*(BB[M][N] + BB[N][M]));
		}
	      }
	    }/* end elastic */
	    
	    /* Elastic or total deformation gradient */
	    eps[ii].il[ip].Fe[0] = eps[ii].il[ip].Fe1[0] = Fn1B[0][0];
	    eps[ii].il[ip].Fe[1] = eps[ii].il[ip].Fe1[1] = Fn1B[0][1];
	    eps[ii].il[ip].Fe[2] = eps[ii].il[ip].Fe1[2] = Fn1B[0][2];

	    eps[ii].il[ip].Fe[3] = eps[ii].il[ip].Fe1[3] = Fn1B[1][0];
	    eps[ii].il[ip].Fe[4] = eps[ii].il[ip].Fe1[4] = Fn1B[1][1];
	    eps[ii].il[ip].Fe[5] = eps[ii].il[ip].Fe1[5] = Fn1B[1][2];

	    eps[ii].il[ip].Fe[6] = eps[ii].il[ip].Fe1[6] = Fn1B[2][0];
	    eps[ii].il[ip].Fe[7] = eps[ii].il[ip].Fe1[7] = Fn1B[2][1];
	    eps[ii].il[ip].Fe[8] = eps[ii].il[ip].Fe1[8] = Fn1B[2][2];

	  }/* end periodic */
	  
	  /* Material stiffness matrix */
	  matrix_tensor_3D (elem[ii].mat[2],hommat,L);
	  
	  /**** CRYSTAL PLASTICITY ****/
	  if (opts->analysis_type == FS_CRPL){
	    /* Strain */
	    for (M=0;M<3;M++){
	      for (N=0;N<3;N++){
		S[M][N] = 0.0;
		for (P=0;P<3;P++){
		  S[M][N] += FnB[M][P]*UU[P][N];
		}
	      }
	    }
	    get_GL_strain (S,Fr,Jr,Tr,E);
	  }
	  else
	    /* Strain */
	    get_GL_strain (FnB,Fr,Jr,Tr,E);
	  
	  /* Stress */
	  get_SPK_stress (L,E,S);
	  /* PGFEM_printf ("[S] pressure =
	     %12.12f\n",(S[0][0]+S[1][1]+S[2][2])/3.); */
	  
	  /**** CRYSTAL PLASTICITY ****/
	  if (opts->analysis_type == FS_CRPL){
	    
	    /* Get stress and plastic rate on slip systems */
	    TA_GA (nss,mat,crpl,sig[ii].il[ip].Har1,Fn1B,S,TA,GA);
	    
	    /* Update hardenenig and shear stress */
	    pom = 0.0;
	    for (M=0;M<nss;M++) {
	      eps[ii].il[ip].GA[M] = GA[M];
	      sig[ii].il[ip].Tau[M] = TA[M];
	      pom += fabs(GA[M]);
	    }
	    
	    /* Acumulative plastic slip : PLC */
	    eps[ii].il[ip].GAMA += dt*pom;
	    eps[ii].GD += ai*aj*ak*J/volume*eps[ii].il[ip].GAMA;
	    
	    /* Hardening stress update */
	    sig[ii].il[ip].Har = sig[ii].il[ip].Har1;
	  }
	  
	  /* Second Piola Kirchoff stress */
	  sig[ii].il[ip].o[0] = S[0][0];
	  sig[ii].il[ip].o[1] = S[1][1];
	  sig[ii].il[ip].o[2] = S[2][2];

	  sig[ii].il[ip].o[3] = S[1][2];
	  sig[ii].il[ip].o[4] = S[0][2];
	  sig[ii].il[ip].o[5] = S[0][1];
	  
	  /* Elastic Green Lagrange strain */
	  eps[ii].il[ip].o[0] = E[0][0];
	  eps[ii].il[ip].o[1] = E[1][1];
	  eps[ii].il[ip].o[2] = E[2][2];

	  eps[ii].il[ip].o[3] = 2.*E[1][2];
	  eps[ii].il[ip].o[4] = 2.*E[0][2];
	  eps[ii].il[ip].o[5] = 2.*E[0][1];
	  
	  /**** UNIT CELL APPROACH ****/
	  if (periodic == 1){
	    if (opts->analysis_type == FS_CRPL){
	      
	      /* total plastic and elastic def. gradients */
	      def_grad_inv (Fp,BB); def_grad_inv (Fn1B,DD);
	      
	      /* First P-K stress */
	      for (P=0;P<3;P++){
		for (R=0;R<3;R++){
		  AA[P][R] = 0.0;
		  CC[P][R] = 0.0;
		  UU_I[P][R] = 0.0;
		  for (U=0;U<3;U++){
		    CC[P][R] += Fn1[P][U]*BB[U][R];
		    UU_I[P][R] += DD[P][U]*Fn1[U][R];
		    for (Q=0;Q<3;Q++){
		      AA[P][R] += Fn1B[P][U]*S[U][Q]*BB[R][Q];
		    }
		  }
		}
	      }
	      
	      for (P=0;P<3;P++){
		for (R=0;R<3;R++){
		  /* First P-K stress */
		  eps[0].P[P][R]  += 1./VVolume*ai*aj*ak*J * AA[P][R];
		  /* Total def. gradient */
		  eps[0].FB[P][R] += (1./VVolume*ai*aj*ak*J 
				      * pow(Tr,1./3.)*pow(Jr,-1./3.)*Fr[P][R]);
		  /* Elastic def. gradient */
		  eps[0].Fe[P][R] += 1./VVolume*ai*aj*ak*J * DD[P][R];
		  /* Plastic def. gradient */
		  eps[0].Fp[P][R] += (1./VVolume*ai*aj*ak*J 
				      * pow(Tr,1./3.)*pow(Jr,-1./3.)*BB[P][R]);
		  Fee[P][R] += 1./VVolume*ai*aj*ak*J * UU_I[P][R];
		  Fpp[P][R] += (1./VVolume*ai*aj*ak*J 
				* pow(Tr,1./3.)*pow(Jr,-1./3.)*CC[P][R]);
		}
	      }
	    }/* opts->analysis_type == FS_CRPL */
	    else{
	      for (P=0;P<3;P++){
		for (R=0;R<3;R++){
		  AA[P][R] = 0.0;
		  for (U=0;U<3;U++){
		    AA[P][R] += Fn1B[P][U]*S[U][R];
		  }
		}
	      } 
	      for (P=0;P<3;P++){
		for (R=0;R<3;R++){
		  eps[0].P[P][R]  += 1./VVolume*ai*aj*ak*J* AA[P][R];
		  eps[0].FB[P][R] += (1./VVolume*ai*aj*ak*J* pow(Tr,1./3.)
				      *pow(Jr,-1./3.)*Fr[P][R]);
		}
	      }
	    }/* elastic */
	    
	    /* COMPUTE ENERGY || S:dE */
	    for (P=0;P<3;P++){
	      for (R=0;R<3;R++){
		AA[P][R] = (1./dt * (pow(Tr,1./3.)*pow(Jr,-1./3.)
				     *Fr[P][R] - FlN[P][R]));
	      }
	    }

	    for (P=0;P<3;P++){
	      for (R=0;R<3;R++){
		BB[P][R] = 0.0;
		for (U=0;U<3;U++){
		  BB[P][R] += (1./2. *(AA[U][P]*pow(Tr,1./3.)
				       *pow(Jr,-1./3.)*Fr[U][R]
				       + pow(Tr,1./3.)*pow(Jr,-1./3.)
				       *Fr[U][P]*AA[U][R]));
		}
	      }
	    }
	    for (P=0;P<3;P++){
	      for (R=0;R<3;R++){
		EL_e += 1./VVolume *ai*aj*ak*J * S[P][R]*BB[P][R];
	      }
	    }
	  }/* end PERIODIC */
	  
	  /* Solve for Cauchy stress and elastic Almansi tensor */
	  def_grad_inv (Fn1B,AA);
	  pom = def_grad_det (Fn1B);
	  for (M=0;M<3;M++){
	    for (N=0;N<3;N++){
	      UU_I[M][N] = 0.0;
	      BB[M][N] = 0.0;
	      for (P=0;P<3;P++){
		for (Q=0;Q<3;Q++){
		  UU_I[M][N] += 1./pom*Fn1B[M][P]*S[P][Q]*Fn1B[N][Q];
		  BB[M][N] += AA[P][M]*E[P][Q]*AA[Q][N];
		}
	      }
	    }
	  }
	  for (M=0;M<3;M++){
	    for (N=0;N<3;N++){
	      S[M][N] = UU_I[M][N];
	      E[M][N] = AA[M][N] = BB[M][N];
	    }
	  }

	  /* PGFEM_printf ("[sig] S0 = %12.12f : S1 = %12.12f : S2 = %12.12f
	     || pressure = %12.12f\n",S[0][0],S[1][1],S[2][2],
	     (S[0][0]+S[1][1]+S[2][2])/3.); */

	  /**** CRYSTAL PLASTICITY ****/
	  if (opts->analysis_type == FS_CRPL){
	    
	    /* En+1 || Total strain */
	    for (M=0;M<3;M++){
	      for (N=0;N<3;N++){
		if (M == N)
		  dij = 1.;
		else
		  dij = 0;
		UU[M][N] = -1./2.*dij;
		
		AA[M][N] = 0.0;
		for (P=0;P<3;P++){
		  AA[M][N] += Fn1B[M][P]*Fp[P][N];
		  
		  for (Q=0;Q<3;Q++){
		    for (R=0;R<3;R++){
		      UU[M][N] += (1./2.*Fp[P][M]*Fn1B[Q][P]
				   *Fn1B[Q][R]*Fp[R][N]);
		    }
		  }
		}
	      }
	    }
	    
	    def_grad_inv (AA,UU_I);
	    
	    for (M=0;M<3;M++){
	      for (N=0;N<3;N++){
		AA[M][N] = 0.0;
		for (P=0;P<3;P++){
		  for (Q=0;Q<3;Q++){
		    AA[M][N] += UU_I[P][M]*UU[P][Q]*UU_I[Q][N];
		  }
		}
	      }
	    }
	    
	    /* Almansi Plastic Strain on element */
	    eps[ii].pl.o[0] += ai*aj*ak*J/volume*(AA[0][0] - E[0][0]);
	    eps[ii].pl.o[1] += ai*aj*ak*J/volume*(AA[1][1] - E[1][1]);
	    eps[ii].pl.o[2] += ai*aj*ak*J/volume*(AA[2][2] - E[2][2]);
	    
	    eps[ii].pl.o[3] += 2.*ai*aj*ak*J/volume*(AA[1][2] - E[1][2]); 
	    eps[ii].pl.o[4] += 2.*ai*aj*ak*J/volume*(AA[0][2] - E[0][2]); 
	    eps[ii].pl.o[5] += 2.*ai*aj*ak*J/volume*(AA[0][1] - E[0][1]);
	    
	    for (P=0;P<3;P++){
	      for (R=0;R<3;R++){
		eps[ii].il[ip].dUU_Tr_n[P][R] = eps[ii].il[ip].dUU_Tr[P][R];
		for (U=0;U<3;U++){
		  for (W=0;W<3;W++){
		    eps[ii].il[ip].dUU_Fr_n[P][R][U][W] =
		      eps[ii].il[ip].dUU_Fr[P][R][U][W];
		  }
		}/* end U */
	      }
	    }/* end  P < 3 */
	    
	  }/* end analysis == FS_CRPL */

	  /* Compute Logarithmic strain */
	  /* Logarithmic_strain (Fn1B,AA); */

	  /* Almansi or Logarithmic strain */
	  eps[ii].el.o[0] += ai*aj*ak*J/volume*AA[0][0];
	  eps[ii].el.o[1] += ai*aj*ak*J/volume*AA[1][1] ;
	  eps[ii].el.o[2] += ai*aj*ak*J/volume*AA[2][2];

	  eps[ii].el.o[3] += 2.*ai*aj*ak*J/volume*AA[1][2];
	  eps[ii].el.o[4] += 2.*ai*aj*ak*J/volume*AA[0][2];
	  eps[ii].el.o[5] += 2.*ai*aj*ak*J/volume*AA[0][1];
	  
	  /* Cauchy Stress */
	  sig[ii].el.o[0] += ai*aj*ak*J/volume*S[0][0];
	  sig[ii].el.o[1] += ai*aj*ak*J/volume*S[1][1] ;
	  sig[ii].el.o[2] += ai*aj*ak*J/volume*S[2][2];

	  sig[ii].el.o[3] += ai*aj*ak*J/volume*S[1][2];
	  sig[ii].el.o[4] += ai*aj*ak*J/volume*S[0][2];
	  sig[ii].el.o[5] += ai*aj*ak*J/volume*S[0][1];

	  
	  ip++;
	}/* end k < KK */
      }/* end j < JJ */
    }/* end i < II */
    
    /* PGFEM_printf("Jn = %12.20f : Jr = %12.20f || Tn = %12.20f : Tr =
       %12.20f\n",Jn,Jr,Tn,Tr); */
    
    /* PRESSURE AND VOLUME CHANGES */
    for (M=0;M<npres;M++){
      if (periodic == 1){
	eps[ii].T[M]  = eps[ii].d_T[M];
	sig[ii].p[M] += sig[ii].d_p[M];
	
	sig[ii].d_p[M] = 0.0;
      }
      else{
	eps[ii].T[M] *= eps[ii].d_T[M];
	sig[ii].p[M] += sig[ii].d_p[M];
	
	eps[ii].d_T[M] = 1.0;
	sig[ii].d_p[M] = 0.0;
      }
    }

    /* PGFEM_printf("JnJr = %12.20f || TnTr = %12.20f || pn+1 =
       %12.12f\n",Jn*Jr,eps[ii].T[0],sig[0].p[0]); */

    dealoc1l (nod);
    dealoc1 (x);
    dealoc1 (y);
    dealoc1 (z);
    dealoc1 (N_x);
    dealoc1 (N_y);
    dealoc1 (N_z);
    dealoc4 (ST,3,3,ndn);
    dealoc1 (r_e);
    dealoc1 (TA);
    dealoc1 (GA);
    dealoc1l (cn);

  }/* end ii < ne */
  
  if (periodic == 1){
    if (opts->analysis_type == FS_CRPL){
      /* Averaging over the domains */
      GPf = aloc1 (nproc*54);
      
      U = 0;
      for (P=0;P<3;P++){
	for (R=0;R<3;R++){
	  LPf[U+0]  = eps[0].P[P][R];
	  LPf[U+9]  = eps[0].FB[P][R];
	  LPf[U+18] = eps[0].Fe[P][R];
	  LPf[U+27] = eps[0].Fp[P][R];
	  LPf[U+36] = Fee[P][R];
	  LPf[U+45] = Fpp[P][R];
	  U++;
	}
      }
      
      MPI_Allgather (LPf,54,MPI_DOUBLE,GPf,54,MPI_DOUBLE,mpi_comm);
      
      for (M=0;M<nproc;M++){
	if (M == myrank) continue;
	U = 0;
	for (P=0;P<3;P++){
	  for (R=0;R<3;R++){
	    eps[0].P[P][R] += GPf[M*54+U+0];
	    eps[0].FB[P][R] += GPf[M*54+U+9];
	    eps[0].Fe[P][R] += GPf[M*54+U+18];
	    eps[0].Fp[P][R] += GPf[M*54+U+27];
	    Fee[P][R] += GPf[M*54+U+36];
	    Fpp[P][R] += GPf[M*54+U+45];
	    U++;
	  }
	}
      }
      
      dealoc1 (GPf);
      
      pom = 0.0;
      for (P=0;P<3;P++){
	for (R=0;R<3;R++){
	  AA[P][R] = Fpp[P][R];
	  CC[P][R] = Fee[P][R];
	  for (U=0;U<3;U++){
	    AA[P][R] += eps[0].FB[P][U]*eps[0].Fp[U][R];
	    CC[P][R] += eps[0].Fe[P][U]*eps[0].FB[U][R];
	  }
	}
      }
      
      for (P=0;P<3;P++){
	for (R=0;R<3;R++){
	  eps[0].Fe[P][R] = AA[P][R];
	  eps[0].Fp[P][R] = CC[P][R];
	}
      }
      
      for (P=0;P<3;P++){
	for (R=0;R<3;R++){
	  BB[P][R] = 0.0;
	  for (U=0;U<3;U++){
	    BB[P][R] += AA[P][U]*CC[U][R];
	  }
	}
      }
      
      if (myrank == 0){
	PGFEM_printf("FB=Fe*Fp\n");
	for (P=0;P<3;P++){
	  for (R=0;R<3;R++){
	    PGFEM_printf("%12.12f  ",BB[P][R]);
	  }
	  PGFEM_printf("\n");
	}
      }
    }/* end analysis == FS_CRPL */
    else {
      /* Averaging over the domains */
      GPf = aloc1 (nproc*27);
      
      U = 0;
      for (P=0;P<3;P++){
	for (R=0;R<3;R++){
	  LPf[U+0]  = eps[0].P[P][R];
	  LPf[U+9]  = eps[0].FB[P][R];
	  LPf[U+18] = eps[0].Fn[P][R];
	  U++;
	}
      }
      
      MPI_Allgather (LPf,27,MPI_DOUBLE,GPf,27,MPI_DOUBLE,mpi_comm);
      
      for (M=0;M<nproc;M++){
	if (M == myrank) continue;
	U = 0;
	for (P=0;P<3;P++){
	  for (R=0;R<3;R++){
	    eps[0].P[P][R] += GPf[M*27+U+0];
	    eps[0].FB[P][R] += GPf[M*27+U+9];
	    eps[0].Fn[P][R] += GPf[M*27+U+18];
	    U++;
	  }
	}
      }

      dealoc1 (GPf);
    }

    /* MACRO ENERGY || P:dF */
    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
	AA[P][R] = 1./dt * (eps[0].F[P][R] - FoN[P][R]);
      }
    }
    
    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
	BB[P][R] = 0.0;
	for (U=0;U<3;U++){
	  BB[P][R] += (1./2. *(AA[U][P]*eps[0].FB[U][R]
			       + eps[0].FB[U][P]*AA[U][R]));
	}
      }
    }
    
    pom = 0.0;
    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
	pom  += eps[0].P[P][R] * AA[P][R];
	eps[0].Fn[P][R] = eps[0].F[P][R];
      }
    }
    
    MPI_Allreduce(&EL_e, &GEL_e, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);  
    if (myrank == 0)
      PGFEM_printf("P:dF = %12.12f || Aver. of Mic. S:dE = %12.12f\n",pom,GEL_e);
  }/* end periodic == 1 */
  
  /* Averaging over the domains */
  GDp = aloc1 (nproc*9);
  
  U = 0;
  for (P=0;P<3;P++){
    for (R=0;R<3;R++){
      LDp[U]  = eps[0].Dp[P][R];
      U++;
    }
  }	
  
  MPI_Allgather (LDp,9,MPI_DOUBLE,GDp,9,MPI_DOUBLE,mpi_comm);
  
  for (M=0;M<nproc;M++){
    if (M == myrank) continue;
    U = 0;
    for (P=0;P<3;P++){
      for (R=0;R<3;R++){
	eps[0].Dp[P][R] += GDp[M*9+U];
	U++;
      }
    }	
  }
  
  dealoc1 (GDp);
  
  pom = 0.0;
  for (P=0;P<3;P++){
    for (R=0;R<3;R++){
      pom += eps[0].Dp[P][R]*eps[0].Dp[P][R];
    }
  }
  /* Macro-scale effective plastic or elastic strain */
  eps[0].eff += dt*sqrt(2./3.*pom);
  
  /* Coordinates update */
  if (periodic == 1){
    for (ii=0;ii<nn;ii++){
      
      X[0] = node[ii].x1_fd;
      X[1] = node[ii].x2_fd;
      X[2] = node[ii].x3_fd;

      for (P=0;P<3;P++){
	Y[P] = 0.0;
	for (R=0;R<3;R++){
	  Y[P] += eps[0].F[P][R]*X[R];
	}
      }
      
      for (i=0;i<ndn;i++){
	II = node[ii].id[i];
	
	if (i == 0) node[ii].x1 = Y[i];
	if (i == 1) node[ii].x2 = Y[i];
	if (i == 2) node[ii].x3 = Y[i];
	
	if (II > 0){
	  if (i == 0) node[ii].x1 += r[II-1] + d_r[II-1];
	  if (i == 1) node[ii].x2 += r[II-1] + d_r[II-1];
	  if (i == 2) node[ii].x3 += r[II-1] + d_r[II-1];
	}
      }
    }
  }/* end periodic */
  else{
    for (ii=0;ii<nn;ii++){
      for (i=0;i<ndn;i++){
	II = node[ii].id[i];
	if (II > 0){
	  if (i == 0) node[ii].x1 += d_r[II-1];
	  if (i == 1) node[ii].x2 += d_r[II-1];
	  if (i == 2) node[ii].x3 += d_r[II-1];
	}
	if (II < 0){
	  if (i == 0) node[ii].x1 += sup->defl_d[abs(II)-1];
	  if (i == 1) node[ii].x2 += sup->defl_d[abs(II)-1];
	  if (i == 2) node[ii].x3 += sup->defl_d[abs(II)-1];
	}
      }
    }/* end ii < nn */
  }
  
  PL = T_VOLUME (ne,ndofn,elem,node);
  /* Gather Volume from all domains */
  MPI_Reduce (&PL,&GVol,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if (opts->cohesive == 1) {
    MPI_Reduce (pores,&Gpores,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    *pores = Gpores;
  }

  if (myrank == 0) {
    if (opts->cohesive == 0)
      PGFEM_printf ("AFTER DEF - VOLUME = %12.12f\n",GVol);
    else
      PGFEM_printf ("AFTER DEF - VOLUME = %12.12f ||"
	      " Volume of voids %12.12f ||"
	      " Vv/V = %12.12f\n",
	      GVol,Gpores,Gpores/VVolume);
  }
  
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
  dealoc1 (Psi);
  dealoc2 (Fn,3);
  dealoc2 (Fr,3);
  dealoc2 (S,3);
  dealoc2 (E,3);
  dealoc2 (Fn1,3);
  dealoc2 (UU,3);
  dealoc2 (UU_I,3);
  dealoc2 (Fp,3);
  dealoc2 (FnB,3);
  dealoc2 (AA,3);
  dealoc2 (BB,3);
  dealoc2 (Fn1B,3);
  dealoc1 (X);
  dealoc1 (Y);
  dealoc2 (CC,3);
  dealoc2 (FoN,3);
  dealoc2 (FlN,3);
  dealoc2 (DD,3);
  dealoc2 (Fee,3);
  dealoc2 (Fpp,3);
}
