#include "initial_guess.h"
#include <math.h>
#include "allocation.h"
#include "def_grad.h"
#include "cast_macros.h"

static const int periodic = 0;

long initial_guess (long ii,
		    long ip,
		    long tim,
		    long iter,
		    long nne,
		    long ndofn,
		    long mat,
		    long nss,
		    CRPL *crpl,
		    long STEP,
		    long GAMA,
		    double *r_r,
		    double ****ST,
		    double **Fr,
		    double **Fn,
		    double **FnB,
		    double Jr,
		    double Tr,
		    double L[3][3][3][3],
		    SIG *sig,
		    EPS *eps,
		    double **UU,
		    double *HAR,
		    double *LAM,
		    double dt)
{
  long INFO = 0,M,N,P,Q;
  double **E,**PP,**BB;
  
  E = aloc2 (3,3); PP = aloc2 (3,3); BB = aloc2 (3,3);
  
  if (periodic == 1){
    if (iter == 0 || GAMA != 0){
      def_grad_inv (CCONST_2(double) Fr,E);
      def_grad_inv (CCONST_2(double) FnB,PP);
      
      BB[0][0] = eps[ii].il[ip].Fe[0];  BB[0][1] = eps[ii].il[ip].Fe[1]; BB[0][2] = eps[ii].il[ip].Fe[2];
      BB[1][0] = eps[ii].il[ip].Fe[3];  BB[1][1] = eps[ii].il[ip].Fe[4]; BB[1][2] = eps[ii].il[ip].Fe[5];
      BB[2][0] = eps[ii].il[ip].Fe[6];  BB[2][1] = eps[ii].il[ip].Fe[7]; BB[2][2] = eps[ii].il[ip].Fe[8];
      
      for (M=0;M<3;M++){for (N=0;N<3;N++){UU[M][N] = 0.0; for (P=0;P<3;P++){
	    for (Q=0;Q<3;Q++){UU[M][N] += pow(Tr,-1./3.)*pow(Jr,1./3.)*PP[M][P]*E[P][Q]*BB[Q][N];}}}}
     
      /* Deviatoric part */
      /* J = def_grad_det (UU); for (M=0;M<3;M++) for (N=0;N<3;N++) UU[M][N] *= pow (J,-1./3.); */
      
      *LAM = 0.0; *HAR = sig[ii].il[ip].Har;
      /* INFO = int_HA (ii,ip,nss,mat,dt,crpl,eps[ii].il[ip].GA,sig[ii].il[ip].Har,HAR); */
    }/* iter == 0 */
    else{
      UU[0][0] = eps[ii].il[ip].UU[0];  UU[0][1] = eps[ii].il[ip].UU[1]; UU[0][2] = eps[ii].il[ip].UU[2];
      UU[1][0] = eps[ii].il[ip].UU[3];  UU[1][1] = eps[ii].il[ip].UU[4]; UU[1][2] = eps[ii].il[ip].UU[5];
      UU[2][0] = eps[ii].il[ip].UU[6];  UU[2][1] = eps[ii].il[ip].UU[7]; UU[2][2] = eps[ii].il[ip].UU[8];
      
      /* deformation tensor */
      
      def_grad_get (nne,ndofn,CONST_4(double) ST,r_r,PP);
      
      for (M=0;M<3;M++){for (N=0;N<3;N++){BB[M][N] = Fr[M][N] - (eps[0].F[M][N] + Fn[M][N] + PP[M][N]);}}
      
      for (M=0;M<3;M++){
	for (N=0;N<3;N++){
	  
	  PP[M][N] = 0.0;
	  for (P=0;P<3;P++){
	    for (Q=0;Q<3;Q++){
	      PP[M][N] += BB[P][Q] * eps[ii].il[ip].dUU_Fr[P][Q][M][N]; /* Tt */
	    }
	  }
	  UU[M][N] += PP[M][N];
	}
      }
      
      /* Deviatoric part */
      /* J = def_grad_det (UU); for (M=0;M<3;M++) for (N=0;N<3;N++) UU[M][N] *= pow (J,-1./3.); */
      
      *LAM = eps[ii].il[ip].lam;  *HAR = sig[ii].il[ip].Har1;
      /* INFO = int_HA (ii,ip,nss,mat,dt,crpl,eps[ii].il[ip].GA1,sig[ii].il[ip].Har1,HAR,eps); */
    }
  }/* end periodic == 1 */
  else{
    if (iter == 0 || GAMA != 0){
      /*
	if (tim == 0 && STEP == 1){for (P=0;P<3;P++){for (Q=0;Q<3;Q++){if(P == Q) dij = 1.0; else dij = 0.0; UU[P][Q] = dij;}}}
	else{
      */
      
      def_grad_inv (CCONST_2(double) Fr,E);
      def_grad_inv (CCONST_2(double) FnB,PP);
      for (M=0;M<3;M++){for (N=0;N<3;N++){UU[M][N] = 0.0; for (P=0;P<3;P++){
	    for (Q=0;Q<3;Q++){UU[M][N] += pow(Tr,-1./3.)*pow(Jr,1./3.)*PP[M][P]*E[P][Q]*FnB[Q][N];}}}}
      
      /* Hardening from step n-1 */
      *LAM = 0.0; *HAR = sig[ii].il[ip].Har;
      /* INFO = int_HA (ii,ip,nss,mat,dt,crpl,eps[ii].il[ip].GA,sig[ii].il[ip].Har,HAR); */
    }
    else{
      UU[0][0] = eps[ii].il[ip].UU[0];  UU[0][1] = eps[ii].il[ip].UU[1]; UU[0][2] = eps[ii].il[ip].UU[2];
      UU[1][0] = eps[ii].il[ip].UU[3];  UU[1][1] = eps[ii].il[ip].UU[4]; UU[1][2] = eps[ii].il[ip].UU[5];
      UU[2][0] = eps[ii].il[ip].UU[6];  UU[2][1] = eps[ii].il[ip].UU[7]; UU[2][2] = eps[ii].il[ip].UU[8];
      
      /* deformation tensor */
      def_grad_get (nne,ndofn,CONST_4(double) ST,r_r,PP);
      
      for (M=0;M<3;M++){for (N=0;N<3;N++){BB[M][N] = Fr[M][N] - PP[M][N];}}
      
      for (M=0;M<3;M++){
	for (N=0;N<3;N++){
	  
	  PP[M][N] = 0.0;
	  for (P=0;P<3;P++){
	    for (Q=0;Q<3;Q++){
	      PP[M][N] += BB[P][Q] * eps[ii].il[ip].dUU_Fr[P][Q][M][N];
	    }
	  }
	  UU[M][N] += PP[M][N];
	}
      }
      
      *LAM = eps[ii].il[ip].lam;  *HAR = sig[ii].il[ip].Har1;
      /* INFO = int_HA (ii,ip,nss,mat,dt,crpl,eps[ii].il[ip].GA1,sig[ii].il[ip].Har1,HAR,eps); */
    }
  }
  
  dealoc2 (E,3); dealoc2 (PP,3); dealoc2 (BB,3);
  
  return (INFO);
}
