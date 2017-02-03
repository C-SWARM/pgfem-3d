#include "integration.h"
#include <math.h>
#include <string.h>
#include "umfpack.h"
#include "enumerations.h"
#include "incl.h"
#include "elem3d.h"
#include "get_dof_ids_on_elem.h"
#include "cast_macros.h"
#include "def_grad.h"
#include "initial_guess.h"
#include "differentiation_of_UU.h"
#include "null.h"
#include "pressu_shape.h"
#include "stress_strain.h"
#include "TA_GA.h"
#include "tensors.h"
#include "utils.h"

static const int periodic = 0;

double int_har_explicit (long ii,
			 long ip,
			 long mat,
			 SIG *sig,
			 CRPL *crpl,
			 double dt)
{
  long i,nss;
  double Tau,Ga,A,GA,Ts,B;
  
  nss = crpl[mat].nss;
  
  GA = 0.0;
  for (i=0;i<nss;i++){
    A = fabs (sig[ii].il[ip].Tau[i]/sig[ii].il[ip].Har);
    Ga = crpl[mat].a*(sig[ii].il[ip].Tau[i]/sig[ii].il[ip].Har)*pow(A,1./crpl[mat].m - 1.);
    GA += fabs (Ga);
  }
  
  B = fabs (GA/crpl[mat].Go);

  Ts = crpl[mat].Tso*pow(B,crpl[mat].mm);
  
  Tau = sig[ii].il[ip].Har + dt*crpl[mat].TH*((Ts - sig[ii].il[ip].Har)/(Ts - crpl[mat].To));
  
  return (Tau);
}

double int_alg_res (double lam,
		    double dt,
		    long mat,
		    long nss,
		    double *TA,
		    double *GA,
		    CRPL *crpl,
		    double **UU,
		    double **Fp,
		    double **BB)
{
  long M,N,P,k;
  double **AA,**E,**PP,pom,dij,DD;
  
  DD = 0.0;
  
  AA = aloc2 (3,3); E = aloc2 (3,3); PP = aloc2 (3,3);
  
  def_grad_inv (CCONST_2(double) UU,AA);
  for (M=0;M<3;M++){for (N=0;N<3;N++){E[M][N] = 0.0; for (P=0;P<3;P++){E[M][N] += AA[M][P]*Fp[P][N];}}}
  pom = def_grad_det (CCONST_2(double) E);
  
  nulld2 (E,3,3);
  for (k=0;k<nss;k++){
    
    PP[0][0] = crpl[mat].P[k][3]*crpl[mat].P[k][0];  PP[0][1] = crpl[mat].P[k][3]*crpl[mat].P[k][1];  PP[0][2] = crpl[mat].P[k][3]*crpl[mat].P[k][2];
    PP[1][0] = crpl[mat].P[k][4]*crpl[mat].P[k][0];  PP[1][1] = crpl[mat].P[k][4]*crpl[mat].P[k][1];  PP[1][2] = crpl[mat].P[k][4]*crpl[mat].P[k][2];
    PP[2][0] = crpl[mat].P[k][5]*crpl[mat].P[k][0];  PP[2][1] = crpl[mat].P[k][5]*crpl[mat].P[k][1];  PP[2][2] = crpl[mat].P[k][5]*crpl[mat].P[k][2];
    
    for (M=0;M<3;M++){
      for (N=0;N<3;N++){
	E[M][N] += GA[k]*PP[M][N];
      }
    }
  }/* end k < nss */
  
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      if (M == N) dij = 1.0; else dij = 0.0;
      BB[M][N] = E[M][N] - 1./dt*(dij - UU[M][N]) - lam*pom*AA[N][M];
    }
  }
  
  DD = pom - 1.;
  
  dealoc2 (AA,3); dealoc2 (E,3); dealoc2 (PP,3);
  
  return (DD);
}

long int_UU (long ii,
	     long ip,
	     double NOR_MIN,
	     double dt,
	     long mat,
	     double *NOR,
	     long nss,
	     CRPL *crpl,
	     double **Fr,
	     double **FnB,
	     double **Fp,
	     double **BB,
	     double *DD,
	     double **S,
	     double L[3][3][3][3],
	     double *TA,
	     double *GA,
	     double Tr,
	     double Jr,
	     double Har,
	     double *LAM,
	     double **UU,
	     long GAMA,
	     long iter)
{
  long INFO = 0,M,N,P,Q,U,V,W,X,it,k;
  double **eFnB,**FstB,nor,dij,EE[3][3][3][3],FF[3][3][3][3],pom,gama,**PP,**CC,dGdT,*f,**KK,*K,*r,**E,dkl,lam=0.0;
  char  *err[]={"inf","-inf","nan"},str[500];
  
  double *Control = (double *) NULL, *Info = (double *) NULL;
  long *Ai,*Ap,NEL,status = 0;
  void *Symbolic, *Numeric;

  null_4d (EE); null_4d (FF);
  
  eFnB = aloc2 (3,3); FstB = aloc2 (3,3); PP = aloc2 (3,3); CC = aloc2 (3,3); f = aloc1 (10); KK = aloc2 (10,10); r = aloc1 (10); E = aloc2 (3,3);

  /* FstB */
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      FstB[M][N] = 0.0;
      for (P=0;P<3;P++){
	FstB[M][N] += pow(Tr,1./3.)*pow(Jr,-1./3.)*Fr[M][P]*FnB[P][N];
      }
    }
  }
  
  it = 1;
  while (*NOR > NOR_MIN){
    
    /* eFnB */
    for (M=0;M<3;M++){
      for (N=0;N<3;N++){
	eFnB[M][N] = 0.0;
	for (P=0;P<3;P++){
	  eFnB[M][N] += FstB[M][P]*UU[P][N];
	}
      }
    }
    
    /******************************************************************************************************************************************************/
    /************************************************* JACOBI MATRIX **************************************************************************************/
    /******************************************************************************************************************************************************/
    
    for (M=0;M<3;M++){
      for (N=0;N<3;N++){
	for (P=0;P<3;P++){
	  if (P == N) dij = 1.; else dij = 0.;
	  for (Q=0;Q<3;Q++){
	    
	    EE[M][N][P][Q] = 0.0; /* Jkqs */
	    FF[M][N][P][Q] = 0.0; /* Jkqs */
	    
	    for (U=0;U<3;U++){
	      EE[M][N][P][Q] += eFnB[U][P]*FstB[U][M]*S[N][Q];
	      for (V=0;V<3;V++){
		EE[M][N][P][Q] += FstB[U][M]*dij*eFnB[U][V]*S[V][Q];
		
		for (W=0;W<3;W++){
		  for (X=0;X<3;X++){
		    if (L[V][Q][N][W] == 0.0) continue;
		    EE[M][N][P][Q] += eFnB[U][P]*eFnB[U][V]*L[V][Q][N][W]*FstB[X][M]*eFnB[X][W];
		  }
		}
	      }
	    }
	    
	  }/* end Q */
	}/* end P */
      }/* end N */
    }/* end M */
    
    /* SLIP SYSTEMS */
    for (k=0;k<nss;k++){
      pom = fabs (TA[k]/Har);
      dGdT = crpl[mat].a/(crpl[mat].m*Har)*pow(pom,1./crpl[mat].m - 1.);
      
      PP[0][0] = crpl[mat].P[k][3]*crpl[mat].P[k][0];  PP[0][1] = crpl[mat].P[k][3]*crpl[mat].P[k][1];  PP[0][2] = crpl[mat].P[k][3]*crpl[mat].P[k][2];
      PP[1][0] = crpl[mat].P[k][4]*crpl[mat].P[k][0];  PP[1][1] = crpl[mat].P[k][4]*crpl[mat].P[k][1];  PP[1][2] = crpl[mat].P[k][4]*crpl[mat].P[k][2];
      PP[2][0] = crpl[mat].P[k][5]*crpl[mat].P[k][0];  PP[2][1] = crpl[mat].P[k][5]*crpl[mat].P[k][1];  PP[2][2] = crpl[mat].P[k][5]*crpl[mat].P[k][2];
      
      
      for (M=0;M<3;M++){
	for (N=0;N<3;N++){
	  CC[M][N] = 0.0; /* Jk */
	  for (U=0;U<3;U++){
	    for (V=0;V<3;V++){
	      CC[M][N] +=  EE[M][N][U][V]*PP[U][V];
	    }
	  }
	}
      }
      
      for (M=0;M<3;M++){
	for (N=0;N<3;N++){
	  for (P=0;P<3;P++){
	    for (Q=0;Q<3;Q++){
	      /* Plus is for sum through slip systems */
	      FF[M][N][P][Q] += dGdT * PP[M][N]*CC[P][Q]; /* UuJk */
	    }
	  }
	}
      }
      
    }/* end k < nss */
    
    /* Inverse of UU */
    def_grad_inv (CCONST_2(double) UU,PP);
    for (M=0;M<3;M++){for (N=0;N<3;N++){ CC[M][N] = 0.0; for (U=0;U<3;U++){ CC[M][N] += PP[M][U]*Fp[U][N];}}}
    
    /* Det. of pFn+1 */
    pom = def_grad_det (CCONST_2(double) CC);
    
    for (M=0;M<3;M++){
      for (N=0;N<3;N++){
	for (P=0;P<3;P++){
	  if (M == P) dij = 1.; else dij = 0.0;
	  for (Q=0;Q<3;Q++){
	    if (N == Q) dkl = 1.; else dkl = 0.0;
	    
	    FF[M][N][P][Q] += 1./dt*dij*dkl + *LAM*pom*(PP[Q][P]*PP[N][M] + PP[N][P]*PP[Q][M]); /* UuJk */
	  }
	}
	
	CC[M][N] = -pom*PP[N][M];
      }
    }
    
    /******************************************************************************************************************************************************/
    /******************************************************************************************************************************************************/
    /******************************************************************************************************************************************************/
    
    tensor_9x9 (KK,FF,0);
    
    KK[0][9]=KK[9][0] = CC[0][0]; KK[1][9]=KK[9][1] = CC[1][1]; KK[2][9]=KK[9][2] = CC[2][2];
    KK[3][9]=KK[9][3] = CC[1][2]; KK[4][9]=KK[9][4] = CC[0][2]; KK[5][9]=KK[9][5] = CC[0][1];
    KK[6][9]=KK[9][6] = CC[2][1]; KK[7][9]=KK[9][7] = CC[2][0]; KK[8][9]=KK[9][8] = CC[1][0]; KK[9][9] = 0.0;
    
    NEL = 0; for (M=0;M<10;M++) for (N=0;N<10;N++) if (KK[N][M] != 0) NEL ++;
    Ap = aloc1l (11); Ai = aloc1l (NEL); K = aloc1 (NEL);
    
    NEL = 0;
    for (M=0;M<10;M++){
      P = 0;
      for (N=0;N<10;N++){
	if (KK[N][M] != 0){
	  
	  Ai[NEL] = N;
	  K[NEL] = KK[N][M];
	  NEL ++; P++;
	}
	Ap[M+1] = Ap[M] + P;
      }
    }
    
    f[0] = -1.*BB[0][0]; f[1] = -1.*BB[1][1]; f[2] = -1.*BB[2][2]; f[3] = -1.*BB[1][2]; f[4] = -1.*BB[0][2]; f[5] = -1.*BB[0][1]; f[6] = -1.*BB[2][1];
    f[7] = -1.*BB[2][0]; f[8] = -1.*BB[1][0]; f[9] = -1.**DD;
    
    /* Solve for dUU */
    
    umfpack_dl_symbolic (10,10,Ap,Ai,K,&Symbolic,Control,Info);
    status = umfpack_dl_numeric (Ap,Ai,K,Symbolic,&Numeric,Control,Info);
    if (status != 0) {INFO = 1; PGFEM_printf("%ld %ld || Singular system of equations in Int. Algorithm \n",ii,ip); break;}
    (void) umfpack_dl_solve (UMFPACK_A,Ap,Ai,K,r,f,Numeric,Control,Info);
    
    umfpack_dl_free_symbolic (&Symbolic);
    umfpack_dl_free_numeric (&Numeric);
    dealoc1 (K); dealoc1l (Ap); dealoc1l (Ai); 
    
    PP[0][0] = r[0]; PP[0][1] = r[5]; PP[0][2] = r[4];
    PP[1][0] = r[8]; PP[1][1] = r[1]; PP[1][2] = r[3];
    PP[2][0] = r[7]; PP[2][1] = r[6]; PP[2][2] = r[2];    
    lam = r[9];
    
    for (M=0;M<3;M++){for (N=0;N<3;N++){CC[M][N] = UU[M][N] + PP[M][N];}}
    
    /* Strain */
    for (M=0;M<3;M++){for (N=0;N<3;N++){S[M][N] = 0.0; BB[M][N] = 0.0; for (P=0;P<3;P++){S[M][N] += FnB[M][P]*CC[P][N]; BB[M][N] += FstB[M][P]*CC[P][N];}}}
    
    get_GL_strain (S,Fr,Jr,Tr,E);
    /* Stress */
    get_SPK_stress (L,E,S);
    
    /*************/
    /* RESIDUALS */
    /*************/
    
    /* Get stress and plastic rate on slip systems */
    TA_GA (nss,mat,crpl,Har,BB,S,TA,GA);
    
    /* Solve for residuals */ pom = *LAM + lam;
    *DD = int_alg_res (pom,dt,mat,nss,TA,GA,crpl,CC,Fp,BB);
    
    nor = 0.0; for (M=0;M<3;M++){for (N=0;N<3;N++){nor += BB[M][N]*BB[M][N];}}
    
    pom = sqrt (nor + *DD**DD); 
    
    sprintf (str,"%f",pom);
    for (N=0;N<3;N++){M = 10; M = strcmp(err[N],str);
      if (M == 0) {if(iter == 0 || GAMA == 2) PGFEM_printf("%ld %ld || Error in the integration algorithm : nor = %s\n",ii,ip,err[N]); INFO = 1; break;}}
    if (INFO == 1) break;
    
    gama = 1.0; V = 1;
    while (pom > *NOR*100){
      if (gama == 1.0) gama = 0.75;
      
      /* PGFEM_printf ("%ld %ld <%ld> INT :: Gama = %2.3f || nor = %12.10e || NOR = %12.10e\n",ii,ip,it,gama,pom,*NOR); */
      
      for (M=0;M<3;M++){for (N=0;N<3;N++){CC[M][N] = UU[M][N] + gama*PP[M][N];}}
      
      /* Strain */
      for (M=0;M<3;M++){for (N=0;N<3;N++){S[M][N] = 0.0; BB[M][N] = 0.0; for (P=0;P<3;P++){S[M][N] += FnB[M][P]*CC[P][N]; BB[M][N] += FstB[M][P]*CC[P][N];}}}
      
      get_GL_strain (S,Fr,Jr,Tr,E);
      
      /* Stress */
      get_SPK_stress (L,E,S);
      
      /* Stress and Strain on slip systems */
      TA_GA (nss,mat,crpl,Har,BB,S,TA,GA); pom = *LAM + gama*lam;
      
      /* Residuals */
      *DD = int_alg_res (pom,dt,mat,nss,TA,GA,crpl,CC,Fp,BB);
      
      nor = 0.0; for (M=0;M<3;M++){for (N=0;N<3;N++){nor += BB[M][N]*BB[M][N];}}
      
      pom = sqrt (nor + *DD**DD); 
      
      sprintf (str,"%f",pom);
      for (N=0;N<3;N++){M = 10; M = strcmp(err[N],str);
	if (M == 0) {if(iter == 0 || GAMA == 2) PGFEM_printf("%ld %ld || Error in the integration algorithm : nor = %s\n",ii,ip,err[N]); INFO = 1; break;}}
      if (INFO == 1) break;
      
      if (pom < *NOR*100) break;
      
      gama = 1./exp(V*1./2.); V++;
      
      if (gama < 0.01) {if(iter == 0 || GAMA == 2) PGFEM_printf ("%ld %ld || Error in the integration algorithm : gama = 0.0\n",ii,ip); INFO = 1; break;}
    }
    if (INFO == 1) break;
    
    *LAM += gama*lam; *NOR = pom; for (M=0;M<3;M++){for (N=0;N<3;N++){UU[M][N] += gama*PP[M][N];}}
    
    /* PGFEM_printf ("[%ld][%ld] | (%ld) nor = %12.10e : lam = %2.10e\n",ii,ip,it,*NOR,*LAM); */
    
    if (it > 50 && *NOR > NOR_MIN) {if(iter == 0 || GAMA == 2)PGFEM_printf("%ld %ld || Error in the integration algorithm : it > 50 | N = %12.12e\n",ii,ip,*NOR); 
      INFO = 1; break;}
    it++;
  }/* end nor > nor_min */
  
  dealoc2 (eFnB,3); dealoc2 (FstB,3); dealoc2 (PP,3); dealoc2 (CC,3); dealoc1 (f); dealoc2 (KK,10); dealoc1 (r); dealoc2 (E,3);

  return (INFO);
}

long int_HA (long ii,
	     long ip,
	     long nss,
	     long mat,
	     double dt,
	     CRPL *crpl,
	     double *GA,
	     double Har,
	     double *HAR,
	     EPS *eps,
	     const PGFem3D_opt *opts)
{
  long INFO = 0,k,M,N;
  double ga,DD,pom,Gama,Rof,Rom,Xi;
  char  *err[]={"inf","-inf","nan"},str[500];

  ga = 0.0; for (k=0;k<nss;k++) ga += fabs(GA[k]);
  
  if (opts->plc == 0){
    DD = fabs (ga/crpl[mat].Go);
    pom = crpl[mat].Tso*pow(DD,crpl[mat].mm);
    
    /* Compute new hardening */
    *HAR = ((pom - crpl[mat].To)*Har + dt*crpl[mat].TH*pom*ga)/(pom - crpl[mat].To + dt*crpl[mat].TH*ga);
  }
  /* PLC */
  if (opts->plc == 1){
    
    Gama = eps[ii].il[ip].GAMA + dt*ga;
    Rof = (crpl[mat].l1/crpl[mat].l2 + crpl[mat].c1*exp(-1./2.*crpl[mat].l2*Gama))*(crpl[mat].l1/crpl[mat].l2 + crpl[mat].c1*exp(-1./2.*crpl[mat].l2*Gama));
    Rom = 2.0*crpl[mat].c1*crpl[mat].c2/(crpl[mat].b*crpl[mat].l2)*(exp(-1./2.*crpl[mat].l2*Gama) - crpl[mat].c3);
    Xi = pow ((crpl[mat].b*Rom/(sqrt(Rof)*crpl[mat].to*ga)),1./3.);
    
    /* Compute new hardening */
    *HAR = crpl[mat].To + 1./3.*crpl[mat].nu*crpl[mat].b*sqrt(Rof) + crpl[mat].fo*(1.0 - exp(-1.*Xi));
  }
  
  sprintf (str,"%f",*HAR);
  for (N=0;N<3;N++){ M = 10; M = strcmp(err[N],str);
    if (M == 0) {PGFEM_printf("%ld %ld || Error in the integration algorithm (Hardening) : nor = %s\n",ii,ip,err[N]); return (1);}}
  
  /* PGFEM_printf ("[%ld][%ld] :: Har = %12.10f || ga = %12.12e : GAMA = %5.5e\n",ii,ip,*HAR,ga,eps[ii].il[ip].GAMA); */
  
  return (INFO);
}

long integration_alg (long ne,
		      long ndofn,
		      long ndofd,
		      long npres,
		      CRPL *crpl,
		      ELEMENT *elem,
		      NODE *node,
		      double *d_r,
		      double *rr,
		      SUPP sup,
		      MATGEOM matgeom,
		      HOMMAT *hommat,
		      EPS *eps,
		      SIG *sig,
		      long tim,
		      long iter,
		      double dt,
		      double nor_min,
		      long STEP,
		      long GAMA,
		      const PGFem3D_opt *opts,
		      const int mp_id)
{
  
  long INFO = 0, ii,i,j,k,II,JJ,KK,ip,nne,ndofe,*nod,M,N,P,R,U,W,Q,mat,ndn,
    nss,IT,GUESS,IT_MAX,*cn;
  double *r_e,*x,*y,*z,Har,**Fn,**Fr,Jr,Jn,Tn,Tr,*Psi,**E,**S,L[3][3][3][3],
    ksi,eta,zet,*gk,*ge,*gz,aj,ai,ak,****ST,*N_x,*N_y,*N_z,J,*w,**UU,**Fp,*GA,*TA;
  double nor,**BB,DD,nor1,**FnB,*rrr,*r_r,tol,**eFnB,**FstB,**AA,lam,HAR_n,TOLER;
  
  char  *err[]={"inf","-inf","nan"},str[500];

  if (opts->plc == 0) {
    TOLER = nor_min;
    nor_min *= 0.1;
    IT_MAX = 20;
  }
  
  if (opts->plc == 1) {
    TOLER = 0.01*nor_min;
    nor_min *= 0.00001;
    IT_MAX = 50;
  }
  
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
  Fp = aloc2 (3,3);
  UU = aloc2 (3,3);
  BB = aloc2 (3,3);
  FnB = aloc2 (3,3);
  rrr = aloc1 (ndofd);
  AA = aloc2 (3,3);
  eFnB = aloc2 (3,3);
  FstB = aloc2 (3,3);
  
  for (ii=0;ii<ne;ii++){
    
    nne = elem[ii].toe;
    ndofe = nne*ndofn;
    mat = elem[ii].mat[2];
    nss = crpl[mat].nss;
    
    nod = aloc1l (nne);
    r_e = aloc1 (ndofe);
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);
    N_x = aloc1 (nne);
    N_y = aloc1 (nne);
    N_z = aloc1 (nne);
    ST = aloc4 (3,3,ndn,nne);
    GA = aloc1 (nss);
    TA = aloc1 (nss);
    r_r = aloc1 (ndofe);
    cn = aloc1l (ndofe);
    
    /* Nodes on element */
    elemnodes (ii,nne,nod,elem);
    
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
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    /* displacement on element */
    def_elem (cn,ndofe,d_r,elem,node,r_r,sup,0);
    
    /* deformation on element */
    for (M=0;M<ndofd;M++)
      rrr[M] = d_r[M] + rr[M];

    def_elem (cn,ndofe,rrr,elem,node,r_e,sup,0);
    
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
	  
	  /* Derivatives of shape functions and Jacobian of integration */
	  J = deriv (ksi,eta,zet,nne,x,y,z,N_x,N_y,N_z);
	  
	  /* Material stiffness matrix */
	  matrix_tensor_3D (elem[ii].mat[2],hommat,L);
	  
	  /* For CRYSTAL PLASTICITY Fn denotes eFn */
	  Fn[0][0] = eps[ii].il[ip].F[0];
	  Fn[0][1] = eps[ii].il[ip].F[1];
	  Fn[0][2] = eps[ii].il[ip].F[2];
	  
	  Fn[1][0] = eps[ii].il[ip].F[3];
	  Fn[1][1] = eps[ii].il[ip].F[4];
	  Fn[1][2] = eps[ii].il[ip].F[5];

	  Fn[2][0] = eps[ii].il[ip].F[6];
	  Fn[2][1] = eps[ii].il[ip].F[7];
	  Fn[2][2] = eps[ii].il[ip].F[8];
	  
	  Fp[0][0] = eps[ii].il[ip].Fp[0];
	  Fp[0][1] = eps[ii].il[ip].Fp[1];
	  Fp[0][2] = eps[ii].il[ip].Fp[2];

	  Fp[1][0] = eps[ii].il[ip].Fp[3];
	  Fp[1][1] = eps[ii].il[ip].Fp[4];
	  Fp[1][2] = eps[ii].il[ip].Fp[5];

	  Fp[2][0] = eps[ii].il[ip].Fp[6];
	  Fp[2][1] = eps[ii].il[ip].Fp[7];
	  Fp[2][2] = eps[ii].il[ip].Fp[8];
	  
	  shape_tensor (nne,ndofn,N_x,N_y,N_z,ST);
	  def_grad_get (nne,ndofn,CONST_4(double) ST,r_e,Fr);
	  Jr = def_grad_det (CCONST_2(double) Fr);
	  Jn = def_grad_det (CCONST_2(double) Fn);
	  
	  /* Pressure shape functions */
	  pressu_shape (npres,ksi,eta,zet,Psi);
	  
	  Tn = 0.0; 
	  Tr = 0.0;
	  
	  for (M=0;M<npres;M++){
	    Tn += Psi[M]*eps[ii].T[M];
	    Tr += Psi[M]*eps[ii].d_T[M];
	  }
	  
	  /* Set eFn-BAR */
	  for (N=0;N<3;N++){
	    for (P=0;P<3;P++){
	      FnB[N][P] = pow(Tn,1./3.)*pow(Jn,-1./3.)*Fn[N][P];
	    }
	  }
	  
	  /*************************** UNIT CELL APPROACH ******************************/
	  if (periodic == 1) {
	    
	    /* Inverse of plastic deformation gradient at time n */
	    def_grad_inv (CCONST_2(double) Fp,FnB);
	    
	    for (N=0;N<3;N++){
	      for (P=0;P<3;P++)
		Fr[N][P] += eps[0].F[N][P] + Fn[N][P];
	    }
	    
	    Jr = def_grad_det (CCONST_2(double) Fr);
	    Jn = Tn = 1.;
	    
	  }/* end PERIODIC */
	  
	  /* Initial Guess */
	  GUESS = 0;
	guessN:
	  INFO = initial_guess (ii,ip,tim,iter,nne,ndofn,mat,nss,crpl,STEP,GAMA,r_r,
				ST,Fr,Fn,FnB,Jr,Tr,L,sig,eps,UU,&Har,&lam,dt);
	  
	  /*
	    PGFEM_printf ("Guess UU\n");
	    PGFEM_printf ("%12.19e %12.19e %12.19e\n",UU[0][0],UU[0][1],UU[0][2]);
	    PGFEM_printf ("%12.19e %12.19e %12.19e\n",UU[1][0],UU[1][1],UU[1][2]);
	    PGFEM_printf ("%12.19e %12.19e %12.19e\n",UU[2][0],UU[2][1],UU[2][2]);
	    PGFEM_printf ("\n");
	    PGFEM_printf ("Guess :: Har = %12.12f and lam = %12.12f\n",Har,lam);
	  */
	  
	  /********************************************************************************/
	  /********************************************************************************/
	  /********************************************************************************/
	  
	  /**************************  NONLINEAR ITERATION FOR UU  ************************/
	  IT = 0;
	  tol = 10.0;
	  HAR_n = Har; 
	  /* if(ii == 0 && ip == 0)
	     PGFEM_printf ("NORMA = %le | nor_min = %le | TOLER = %le\n",NORMA,nor_min,TOLER); */
	  while (tol > TOLER){
	    
	    /* Strain */
	    for (M=0;M<3;M++){
	      for (N=0;N<3;N++){
		S[M][N] = 0.0;
		BB[M][N] = 0.0;
		for (P=0;P<3;P++){
		  S[M][N] += FnB[M][P]*UU[P][N];
		  for (Q=0;Q<3;Q++)
		    BB[M][N] +=  pow(Tr,1./3.)*pow(Jr,-1./3.)*Fr[M][P]*FnB[P][Q]*UU[Q][N];
		}
	      }
	    }
	    
	    get_GL_strain (S,Fr,Jr,Tr,E);
	    /* Stress */
	    get_SPK_stress (L,E,S);
	    
	    /* Shear stress on slip systems */
	    TA_GA (nss,mat,crpl,Har,BB,S,TA,GA);

	    /* Solve for residuals */
	    DD = int_alg_res (lam,dt,mat,nss,TA,GA,crpl,UU,Fp,BB);
	    
	    nor1 = 0.0;
	    for (M=0;M<3;M++){
	      for (N=0;N<3;N++)
		nor1 += BB[M][N]*BB[M][N];
	    }
	    
	    nor = sqrt (nor1 + DD*DD);
	    
	    sprintf (str,"%f",nor);

	    for (N=0;N<3;N++){
	      M = 10;
	      M = strcmp(err[N],str);
	      if (M == 0) {
		if(iter == 0 || GAMA == 2)
		  PGFEM_printf("%ld %ld || ERROR in the integration algorithm : nor = %s\n",
			 ii,ip,err[N]);
		INFO = 1;
		break;
	      }
	    }

	    if (INFO == 1) break;
	    
	    /* PLC */
	    /* if (opts->plc == 1 && nor < nor_min && IT > 0) break; */
	    
	    /* PGFEM_printf ("[%ld][%ld] |::| [%ld] nor = %2.8e : pom = %2.5f : lam = %2.8e\n",
	       ii,ip,IT,tol,DD+1.,lam); */
	    
	    /* Integration of UU */
	    INFO = int_UU (ii,ip,nor_min,dt,mat,&nor,nss,crpl,Fr,FnB,Fp,BB,&DD,S,L,
			   TA,GA,Tr,Jr,Har,&lam,UU,GAMA,iter);

	    if (INFO == 1)
	      break;
	    
	    /* Intergration of Har */
	    INFO = int_HA (ii,ip,nss,mat,dt,crpl,GA,sig[ii].il[ip].Har,&Har,eps,opts);
	    if (INFO == 1)
	      break;
	    
	    tol = sqrt((Har - HAR_n)*(Har - HAR_n))/sig[ii].il[ip].Har;
	    HAR_n = Har;
	    
	    if (IT > IT_MAX && tol > nor_min){
	      if(iter == 0 || GAMA == 2){
		PGFEM_printf("%ld %ld || Error in the integration algorithm",ii,ip);
		PGFEM_printf(" : IT > %ld :: tol = %5.5e\n",IT_MAX,TOLER);
	      }
	      INFO = 1;
	      break;
	    }

	    IT++;
	  }  /* end tol > nor_min */
	  if (INFO == 1 && iter > 0 && GUESS == 0) {
	    GAMA = 2;
	    GUESS = 1;
	    goto guessN;
	  }

	  if (INFO == 1) break;
	  if (GAMA == 2) GAMA = 0;
	  
	  /*
	    for (M=0;M<3;M++){for (N=0;N<3;N++){S[M][N] = 0.0; BB[M][N] = 0.0;
	    for (P=0;P<3;P++){S[M][N] += FnB[M][P]*UU[P][N];
	    for (Q=0;Q<3;Q++) BB[M][N] +=  Fr[M][P]*FnB[P][Q]*UU[Q][N];}}}
	    get_GL_strain (S,Fr,Jr,Tr,E);
	    get_SPK_stress (L,E,S);
	    
	    PGFEM_printf ("Jr=%12.12f : Tr=%12.12f : Jn=%12.12f : Tn=%12.12f\n",Jr,Tr,Jn,Tn);
	    
	    PGFEM_printf ("Fr\n");
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",Fr[0][0],Fr[0][1],Fr[0][2]);
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",Fr[1][0],Fr[1][1],Fr[1][2]);
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",Fr[2][0],Fr[2][1],Fr[2][2]);
	    PGFEM_printf ("\n");
	    
	    PGFEM_printf ("Fe\n");
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",BB[0][0],BB[0][1],BB[0][2]);
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",BB[1][0],BB[1][1],BB[1][2]);
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",BB[2][0],BB[2][1],BB[2][2]);
	    PGFEM_printf ("\n");
	    
	    def_grad_inv (UU,AA);
	    
	    PGFEM_printf ("rFp\n");
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",AA[0][0],AA[0][1],AA[0][2]);
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",AA[1][0],AA[1][1],AA[1][2]);
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",AA[2][0],AA[2][1],AA[2][2]);
	    PGFEM_printf ("\n");
	    
	    PGFEM_printf ("E\n");
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",E[0][0],E[0][1],E[0][2]);
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",E[1][0],E[1][1],E[1][2]);
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",E[2][0],E[2][1],E[2][2]);
	    PGFEM_printf ("\n");
	    
	    PGFEM_printf ("S\n");
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",S[0][0],S[0][1],S[0][2]);
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",S[1][0],S[1][1],S[1][2]);
	    PGFEM_printf ("%12.19f %12.19f %12.19f\n",S[2][0],S[2][1],S[2][2]);
	    PGFEM_printf ("\n");
	    
	    PGFEM_printf ("UU\n");
	    PGFEM_printf ("%12.19e %12.19e %12.19e\n",UU[0][0],UU[0][1],UU[0][2]);
	    PGFEM_printf ("%12.19e %12.19e %12.19e\n",UU[1][0],UU[1][1],UU[1][2]);
	    PGFEM_printf ("%12.19e %12.19e %12.19e\n",UU[2][0],UU[2][1],UU[2][2]);
	    PGFEM_printf ("\n");
	    
	    scanf ("%ld",&M);
	  */
	  
	  /* Resulting increment of Plastic def. Gradient */
	  eps[ii].il[ip].UU[0] = UU[0][0];
	  eps[ii].il[ip].UU[1] = UU[0][1];
	  eps[ii].il[ip].UU[2] = UU[0][2];

	  eps[ii].il[ip].UU[3] = UU[1][0];
	  eps[ii].il[ip].UU[4] = UU[1][1];
	  eps[ii].il[ip].UU[5] = UU[1][2];

	  eps[ii].il[ip].UU[6] = UU[2][0];
	  eps[ii].il[ip].UU[7] = UU[2][1];
	  eps[ii].il[ip].UU[8] = UU[2][2];
	  
	  /* Resulting hardness and lagrange multipliers */
	  eps[ii].il[ip].lam = lam;
	  sig[ii].il[ip].Har1 = Har;
	  
	  for (P=0;P<3;P++){
	    for (R=0;R<3;R++){
	      eFnB[P][R] = 0.0;
	      FstB[P][R] = 0.0;
	      S[P][R] = 0.0;
	      
	      for (U=0;U<3;U++){
		FstB[P][R] += pow(Tr,1./3.)*pow(Jr,-1./3.)*Fr[P][U]*FnB[U][R];
		S[P][R] += FnB[P][U]*UU[U][R];
		for (W=0;W<3;W++){
		  eFnB[P][R] += pow(Tr,1./3.)*pow(Jr,-1./3.)*Fr[P][U]*FnB[U][W]*UU[W][R];
		}
	      }
	    }
	  }

	  /* Elastic deformation gradient */
	  if (periodic == 1){
	    eps[ii].il[ip].Fe1[0] = eFnB[0][0];
	    eps[ii].il[ip].Fe1[1] = eFnB[0][1];
	    eps[ii].il[ip].Fe1[2] = eFnB[0][2];
	    
	    eps[ii].il[ip].Fe1[3] = eFnB[1][0];
	    eps[ii].il[ip].Fe1[4] = eFnB[1][1];
	    eps[ii].il[ip].Fe1[5] = eFnB[1][2];

	    eps[ii].il[ip].Fe1[6] = eFnB[2][0];
	    eps[ii].il[ip].Fe1[7] = eFnB[2][1];
	    eps[ii].il[ip].Fe1[8] = eFnB[2][2];
	  }

	  /* Strain */
	  get_GL_strain (S,Fr,Jr,Tr,E);
	  /* Stress */
	  get_SPK_stress (L,E,S);
	  
	  /* Get stress and plastic rate on slip systems */
	  TA_GA (nss,mat,crpl,Har,eFnB,S,TA,GA);
	  for (P=0;P<nss;P++)
	    eps[ii].il[ip].GA1[P] = GA[P];

	  /* Differentiation of UU with respect to Fr and Tr */
	  differentiation_of_UU (ii,ip,mat,crpl,Tr,Jr,Fr,FnB,eFnB,FstB,S,L,UU,sig,eps,dt,opts);
	  
	  ip++;
	}/* end k < KK */
	if (INFO == 1) break;
      }/* end j < JJ */
      if (INFO == 1) break;
    }/* end i < II */
    
    dealoc1l (nod);
    dealoc1 (r_e);
    dealoc1 (x);
    dealoc1 (y);
    dealoc1 (z);
    dealoc1 (N_x);
    dealoc1 (N_y);
    dealoc1 (N_z);
    dealoc4 (ST,3,3,ndn);
    dealoc1 (GA);
    dealoc1 (TA);
    dealoc1 (r_r);
    dealoc1l (cn);
    
    if (INFO == 1)
      break;
  }/* end ii < ne */
  
  dealoc1 (gk);
  dealoc1 (ge);
  dealoc1 (gz);
  dealoc1 (w);
  dealoc1 (Psi);
  dealoc2 (Fn,3);
  dealoc2 (Fr,3);
  dealoc2 (S,3);
  dealoc2 (E,3);
  dealoc2 (Fp,3);
  dealoc2 (UU,3);
  dealoc2 (BB,3);
  dealoc2 (FnB,3);
  dealoc1 (rrr);
  dealoc2 (AA,3);
  dealoc2 (eFnB,3);
  dealoc2 (FstB,3);
    
  return (INFO);
}
