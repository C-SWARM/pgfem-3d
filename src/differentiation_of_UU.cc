/* HEADER */

#include "differentiation_of_UU.h"
#include <math.h>
#include "umfpack.h"
#include "def_grad.h"
#include "allocation.h"
#include "null.h"
#include "PLC.h"
#include "utils.h"
#include "TA_GA.h"
#include "tensors.h"

void differentiation_of_UU (long ii,
			    long ip,
			    long mat,
			    CRPL *crpl,
			    double Tr,
			    double Jr,
			    double **Fr,
			    double **FnB,
			    double **eFnB,
			    double **FstB,
			    double **S,
			    double L[3][3][3][3],
			    double **UU,
			    SIG *sig,
			    EPS *eps,
			    double dt,
			    const PGFem3D_opt *opts)
{
  long k,nss,M,N,P,Q,U,V,W,X,Z;
  double pom,dGdT,Har,*TA,*GA,CC[3][3][3][3],DD[3][3][3][3],
    dij,dkl,**FrI,**A,**B,**EE,**PP,AA[3][3][3][3],
    BB[3][3][3][3],**bb,**K,**KK,*inv,*f,**Ab,**Bb,ga;
  double Ts,POM,dGdt,N5,H_n,N3,N4; /*H,N1,N2;*/
  
  double *Control = NULL, *Info = NULL,*XX;
  long *Ai,*Ap,NEL,status;
  void *Symbolic, *Numeric;
  
  null_4d (CC);
  null_4d (DD);
  null_4d (AA);
  null_4d (BB);
  
  nss = crpl[mat].nss;
  
  TA = aloc1 (nss);
  GA = aloc1 (nss);
  FrI = aloc2 (3,3);
  A = aloc2 (3,3);
  B = aloc2 (3,3);
  PP = aloc2 (3,3);
  EE = aloc2 (3,3);
  bb = aloc2 (3,3);
  K = aloc2 (9,9);

  KK = aloc2 (9,9);
  Ab = aloc2 (3,3);
  Bb = aloc2 (3,3);
  f = aloc1 (9);
  XX = aloc1 (9);
  
  /* Invert Fr */
  def_grad_inv ((const double *const *)Fr,FrI);
  
  /* Implicit integration of hardening */
  Har = sig[ii].il[ip].Har1;
  H_n = sig[ii].il[ip].Har;
  
  /* Get stress and plastic rate on slip systems */
  TA_GA (nss,mat,crpl,Har,eFnB,S,TA,GA);
  
  ga = 0.0; dij = 0.0;
  for (k=0;k<nss;k++){
    pom = fabs (TA[k]/Har);
    if (fabs(GA[k]) > 0.0)
      dij += (-1.*GA[k]/fabs(GA[k])*crpl[mat].a
	      *TA[k]/(Har*Har*crpl[mat].m)
	      *pow(pom,1./crpl[mat].m - 1.));
    ga += fabs(GA[k]);
  }
  
  if (opts->plc == 0){
    pom = fabs (ga/crpl[mat].Go);
    Ts = crpl[mat].Tso*pow(pom,crpl[mat].mm);
    
    /*
      H = crpl[mat].TH*((Ts - Har)/(Ts - crpl[mat].To));
      if (pom == 0.0)
        N1 = 0.0;
      else
        N1 = (crpl[mat].TH*(Har - crpl[mat].To)
	      /((Ts - crpl[mat].To)*(Ts - crpl[mat].To))
	      * crpl[mat].Tso*crpl[mat].mm/crpl[mat].Go * 
	      pow(pom,crpl[mat].mm - 1.));

      N2 = crpl[mat].TH*dt/(Ts - crpl[mat].To);
      N3 = dt*(N1*ga + H)/(1. + N2*ga - dt*(N1*ga + H)*dij);
    */
    
    N5 = (dt*crpl[mat].TH*((crpl[mat].Tso*crpl[mat].mm/crpl[mat].Go
			    * pow(pom,crpl[mat].mm - 1.))
			   *ga*(H_n-crpl[mat].To+dt*crpl[mat].TH*ga)
			   + H_n*crpl[mat].To + Ts*(Ts-crpl[mat].To-H_n))
	  /((Ts-crpl[mat].To+dt*crpl[mat].TH*ga)
	    *(Ts-crpl[mat].To+dt*crpl[mat].TH*ga))
	  );
    
    N3 = N5/(1. - N5*dij);    
  }
  /* PLC */
  if (opts->plc == 1){
    N5 = PLC_diff_g_ga (ii,ip,nss,mat,dt,crpl,GA,eps);
    /* if (N5 < 0.0)
         PGFEM_printf("%ld %ld || POZOR N5 = %f < 0.0 :: ga = %f\n",ii,ip,N5,ga); */
    N3 = N5/(1. - N5*dij);
  }
  
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      Ab[M][N] = 0.0; Bb[M][N] = 0.0;
      for (P=0;P<3;P++){
	Ab[M][N] += FstB[P][M]*FstB[P][N];
	Bb[M][N] += eFnB[P][M]*eFnB[P][N];
      }
    }
  }
  
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      
      EE[M][N] = 0.0; /* qs */
      
      for (P=0;P<3;P++){
	if (P == N) dij = 1.; else dij = 0.;
	for (Q=0;Q<3;Q++){
	  
	  CC[M][N][P][Q] = 0.0; /* Ttqs */
	  DD[M][N][P][Q] = 0.0; /* jKqs */
	  
	  for (U=0;U<3;U++){
	    
	    EE[M][N] += (1./(3.*Tr)*(UU[P][M]*eFnB[Q][U]
				     + eFnB[Q][M]*UU[P][U])
			 * FstB[Q][P]*S[U][N]);
	    
	    for (V=0;V<3;V++){
	      if (V == N) dkl = 1.;
	      else dkl = 0.;
	      
	      CC[M][N][P][Q] += ((dij*eFnB[U][V] + eFnB[U][P]*dkl)
				 * FstB[U][M]*S[V][Q]);

	      DD[M][N][P][Q] += (pow(Tr,1./3.)*pow(Jr,-1./3.)
				 *(UU[U][P]*eFnB[M][V] + eFnB[M][P]*UU[U][V])
				 * FnB[N][U]*S[V][Q]);
	      
	      for (W=0;W<3;W++){
		
		if (L[U][Q][N][W] != 0.0)
		  CC[M][N][P][Q] += Bb[P][U]*Ab[M][V]*UU[V][W]*L[U][Q][N][W];
		
		DD[M][N][P][Q] += (-1./3.
				   *(UU[U][P]*eFnB[V][W] + eFnB[V][P]*UU[U][W])
				   *FrI[N][M]*FstB[V][U]*S[W][Q]);
		
		if (L[P][N][Q][U] != 0.0)
		  EE[M][N] += (1./(3.*Tr)*L[P][N][Q][U]*Bb[M][P]
			       *Ab[V][W]*UU[V][Q]*UU[W][U]);
		
		for (X=0;X<3;X++){
		  for (Z=0;Z<3;Z++){
		    if (L[U][Q][V][W] == 0.0) continue;
		    DD[M][N][P][Q] += (pow(Tr,1./3.)*pow(Jr,-1./3.)
				       *L[U][Q][V][W]*Bb[P][U]
				       *UU[X][V]*UU[Z][W]*FstB[M][Z]*FnB[N][X]
				       -1./3.*L[U][Q][V][W]*Bb[P][U]
				       *UU[X][V]*UU[Z][W]*FrI[N][M]*Ab[X][Z]);
		  }
		}
	      }
	    }
	  }
	  
	}/* end Q */
      }/* end P */
    }
  }
  
  POM = 0.0;
  for (M=0;M<3;M++){
    for (N=0;N<3;N++){
      Ab[M][N] = 0.0;
      Bb[M][N] = 0.0;
    }
  }
  
  for (k=0;k<nss;k++){
    pom = fabs (TA[k]/Har);
    
    if (fabs(GA[k]) > 0.0)
      N4 = (GA[k]/fabs(GA[k])*crpl[mat].a/(crpl[mat].m*Har)
	    *pow(pom,1./crpl[mat].m - 1.));
    else
      N4 = 0.0;
    
    PP[0][0] = crpl[mat].P[k][3]*crpl[mat].P[k][0];
    PP[0][1] = crpl[mat].P[k][3]*crpl[mat].P[k][1];
    PP[0][2] = crpl[mat].P[k][3]*crpl[mat].P[k][2];
    
    PP[1][0] = crpl[mat].P[k][4]*crpl[mat].P[k][0];
    PP[1][1] = crpl[mat].P[k][4]*crpl[mat].P[k][1];
    PP[1][2] = crpl[mat].P[k][4]*crpl[mat].P[k][2];

    PP[2][0] = crpl[mat].P[k][5]*crpl[mat].P[k][0];
    PP[2][1] = crpl[mat].P[k][5]*crpl[mat].P[k][1];
    PP[2][2] = crpl[mat].P[k][5]*crpl[mat].P[k][2];
    
    /* Plus is for sum through slip systems */
    for (M=0;M<3;M++){
      for (N=0;N<3;N++){
	
	POM += N4 * EE[M][N]*PP[M][N];
	
	for (U=0;U<3;U++){
	  for (V=0;V<3;V++){
	    Ab[M][N] += N4 * CC[M][N][U][V]*PP[U][V]; /* Tt */
	    Bb[M][N] += N4 * DD[M][N][U][V]*PP[U][V];
	  }
	}
      }
    }
    
  }/* end k < nss */
  
  for (k=0;k<nss;k++){
    pom = fabs (TA[k]/Har);
    dGdT = crpl[mat].a/(crpl[mat].m*Har)*pow(pom,1./crpl[mat].m - 1.);
    dGdt = (-1.*crpl[mat].a*TA[k]/(Har*Har*crpl[mat].m)
	    *pow(pom,1./crpl[mat].m - 1.));

    PP[0][0] = crpl[mat].P[k][3]*crpl[mat].P[k][0];
    PP[0][1] = crpl[mat].P[k][3]*crpl[mat].P[k][1];
    PP[0][2] = crpl[mat].P[k][3]*crpl[mat].P[k][2];
    
    PP[1][0] = crpl[mat].P[k][4]*crpl[mat].P[k][0];
    PP[1][1] = crpl[mat].P[k][4]*crpl[mat].P[k][1];
    PP[1][2] = crpl[mat].P[k][4]*crpl[mat].P[k][2];

    PP[2][0] = crpl[mat].P[k][5]*crpl[mat].P[k][0];
    PP[2][1] = crpl[mat].P[k][5]*crpl[mat].P[k][1];
    PP[2][2] = crpl[mat].P[k][5]*crpl[mat].P[k][2];

    pom = 0.0;
    for (M=0;M<3;M++){
      for (N=0;N<3;N++){
	A[M][N] = 0.0; B[M][N] = 0.0;
	
	pom +=  EE[M][N]*PP[M][N];
	
	for (U=0;U<3;U++){
	  for (V=0;V<3;V++){
	    A[M][N] += CC[M][N][U][V]*PP[U][V];
	    B[M][N] += DD[M][N][U][V]*PP[U][V];
	  }
	}
      }
    }
    
    /* Plus is for sum through slip systems */
    for (M=0;M<3;M++){
      for (N=0;N<3;N++){
	/* Ii */
	bb[M][N] += -1.*dt*dGdT * pom*PP[M][N] - 1.*dt*dGdt*N3 * POM*PP[M][N];

	for (P=0;P<3;P++){
	  for (Q=0;Q<3;Q++){
	    AA[M][N][P][Q] += (+1.*dt*dGdT * A[M][N]*PP[P][Q]
			       + 1.*dt*dGdt*N3* Ab[M][N]*PP[P][Q]); /* TtIi */

	    BB[M][N][P][Q] += (-1.*dt*dGdT * B[M][N]*PP[P][Q]
			       - 1.*dt*dGdt*N3* Bb[M][N]*PP[P][Q]); /* jKIi */
	  }
	}
      }
    }
  }/*end k < nss */
  
  /* TtIi */
  AA[0][0][0][0] += 1.;
  AA[1][1][1][1] += 1.;
  AA[2][2][2][2] += 1.;
  AA[1][2][1][2] += 1.;
  AA[0][2][0][2] += 1.;
  AA[0][1][0][1] += 1.;
  AA[2][1][2][1] += 1.;
  AA[2][0][2][0] += 1.;
  AA[1][0][1][0] += 1.;
  
  /* Invert AA tensor */
  tensor_9x9 (K,AA,0);
  
  NEL = 0;
  for (M=0;M<9;M++)
    for (N=0;N<9;N++)
      if (K[N][M] != 0)
	NEL ++;
  Ap = aloc1l (10);
  Ai = aloc1l (NEL);
  inv = aloc1 (NEL);
  
  NEL = 0;
  for (M=0;M<9;M++){
    P = 0;
    for (N=0;N<9;N++){
      if (K[N][M] != 0){
	
	Ai[NEL] = N;
	inv[NEL] = K[N][M];
	NEL ++;
	P++;
      }
      Ap[M+1] = Ap[M] + P;
    }
  }
  
  /* Solve for dUU */
  (void) umfpack_dl_symbolic (9,9,Ap,Ai,inv,&Symbolic,Control,Info);
  status = umfpack_dl_numeric (Ap,Ai,inv,Symbolic,&Numeric,Control,Info);
  
  if (status == 0){
    for (M=0;M<9;M++){
      for (N=0;N<9;N++)
	if (N == M) f[N] = 1.;
	else f[N] = 0.0;
      
      status = umfpack_dl_solve (UMFPACK_A,Ap,Ai,inv,XX,f,Numeric,Control,Info);
      
      for (N=0;N<9;N++) KK[N][M] = XX[N];
    }
  }
  else{
    for (M=0;M<9;M++)
      for (N=0;N<9;N++)
	KK[N][M] = 0.0;
  }
  
  umfpack_dl_free_symbolic (&Symbolic);
  umfpack_dl_free_numeric (&Numeric);
  dealoc1l (Ap);
  dealoc1l (Ai);
  dealoc1 (inv);
  
  tensor_9x9 (KK,CC,1);
  
  if (status == 0){
    for (M=0;M<3;M++){
      for (N=0;N<3;N++){
	
	eps[ii].il[ip].dUU_Tr[M][N] = 0.0;
	
	for (P=0;P<3;P++){
	  for (Q=0;Q<3;Q++){
	    
	    eps[ii].il[ip].dUU_Tr[M][N] += bb[P][Q]*CC[P][Q][M][N]; /* Tt */
	    
	    eps[ii].il[ip].dUU_Fr[M][N][P][Q] = 0.0;
	    
	    for (U=0;U<3;U++){
	      for (V=0;V<3;V++){
		eps[ii].il[ip].dUU_Fr[M][N][P][Q] += 
		  BB[M][N][U][V]*CC[U][V][P][Q]; /* jkTt */
	      }
	    }
	  }
	}/* end P */
      }
    }
  }/* end status == 0 */
  else{
    for (M=0;M<3;M++){
      for (N=0;N<3;N++){
	eps[ii].il[ip].dUU_Tr[M][N] = eps[ii].il[ip].dUU_Tr_n[M][N];
	for (P=0;P<3;P++){
	  for (Q=0;Q<3;Q++){
	    eps[ii].il[ip].dUU_Fr[M][N][P][Q] =
	      eps[ii].il[ip].dUU_Fr_n[M][N][P][Q];
	  }
	}
      }
    }
  }/* end status != 0 */
  
  dealoc1 (TA);
  dealoc1 (GA);
  dealoc2 (FrI,3);
  dealoc2 (A,3);
  dealoc2 (B,3);
  dealoc2 (PP,3);
  dealoc2 (EE,3);
  dealoc2 (bb,3);
  dealoc2 (K,9);
  dealoc2 (KK,9);

  dealoc2 (Ab,3);
  dealoc2 (Bb,3);
  dealoc1 (f);
  dealoc1 (XX);
}
