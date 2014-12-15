#include "condense.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mkl_cblas.h"
#include "index_macros.h"
#include "cast_macros.h"
#include "allocation.h"
#include "utils.h"


//*
// KttI:       [nVol]*[nVol]
// Kut:     [nne*nsd]*[nVol]
// KutKttI: [nne*nsd]*[nVol]
// Kpt:       [npres]*[nVol]
// KptKttI:   [npres]*[nVol]
	
// C[m,n] = a*A[m,k] x B[k,n] + b*C[m,n]
// cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,a,A,k,B,n,b,C,n);
//*/
void condense_Fupt_to_Fup(double *fe, int nne, int nsd, int npres, int nVol,
                   double *fu, double *ft, double *fp, double *Kut, double *Ktt, double *Kpt)
{ 
  int ndofn = nsd + 1; 
  double *KttI, *KutKttI, *KptKttI;
  KttI = aloc1(nVol*nVol);
  KutKttI = aloc1(nne*nsd*nVol);
  KptKttI = aloc1(npres*nVol);
  
	inverse(Ktt,nVol,KttI);
     
  double *_fu = aloc1(nne*nsd);
  double *_fp = aloc1(npres);
    
     
	// KttI:       [nVol]*[nVol]
	// Kut:     [nne*nsd]*[nVol]
	// KutKttI: [nne*nsd]*[nVol]
	// Kpt:       [npres]*[nVol]
	// KptKttI:   [npres]*[nVol]
        
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,nVol,nVol,1.0,Kut,nVol,KttI,nVol,0.0,KutKttI,nVol);
       
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              npres,nVol,nVol,1.0,Kpt,nVol,KttI,nVol,0.0,KptKttI,nVol);        
              
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,1,nVol,1.0,KutKttI,nVol,ft,1,0.0,_fu,1);
        
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              npres,1,nVol,1.0,KptKttI,nVol,ft,1,0.0,_fp,1);
                
	for(long a = 0; a<nne; a++)
	{
	  for(long b=0; b<ndofn; b++)
	  {
       if(b<nsd)
         fe[a*ndofn + b] = (fu[a*nsd + b] - _fu[a*nsd + b]);
       else
         fe[a*ndofn + b] = (fp[a] - _fp[a]);
	  }
	}  
  
  free(_fu);
  free(_fp);
  
  free(KttI);
  free(KutKttI);
  free(KptKttI);  
  return;
}
/*
void condense_Fupt_to_Fup(double *fe, int nne, int nsd, int npres, int nVol,
                   double *fu, double *ft, double *fp, double *Kut, double *Ktt, double *Kpt)
{ 
  int ndofn = nsd + 1; 
  Matrix(double) KttI, KutKttI, KptKttI, _fu, _fp;
  
	Matrix_construct(double, KttI);
	Matrix_construct(double, KutKttI);
	Matrix_construct(double, KptKttI);
	Matrix_construct(double, _fu;
	Matrix_construct(double, _fp);  

  Matrix_redim(KttI, nVol, nVol);
  Matrix_redim(KutKttI, nne*nsd, nVol);
  Matrix_redim(KptKttI, npres, nVol);
  
	Matrix_inv(Ktt.m_pdata,KttI.m_pdata);
     
  Matrix_redim(_fu, nne*nsd, 1);
  Matrix_redim(_fp, npres, 1);
    
  // Matrix_AxB(C, a, A, AT, B, BT) <-- C = aAxB
  Matrix_AxB(KutKttI, 1.0,     Kut, 0, KttI, 0);
  Matrix_AxB(KptKttI, 1.0,     Kpt, 0, KttI, 0);  
  Matrix_AxB(    _fu, 1.0, KutKttI, 0,   ft, 0);
  Matrix_AxB(    _fp, 1.0, KptKttI, 0,   ft, 0);  
                        
	for(long a = 0; a<nne; a++)
	{
	  for(long b=0; b<ndofn; b++)
	  {
       if(b<nsd)
         fe[a*ndofn + b] = fu[a*nsd + b] - Mat_v(_fu, a*nsd + b+1, 1);
       else
         fe[a*ndofn + b] = fp[a] - Mat_v(_fp, a, 1);
	  }
	}  
  
  Matrix_cleanup(_fu);
  Matrix_cleanup(_fp);
  
  Matrix_cleanup(KttI);
  Matrix_cleanup(KutKttI);
  Matrix_cleanup(KptKttI);  
  return;
}

*/
void condense_Fupt_to_Fu(double *fe, int nne, int nsd, int npres, int nVol,
                   double *fu, double *ft, double *fp, double *Kup, double *Ktp, double *Ktt,double *Kpt)                               
{
  
	int m = nVol;
	
  double *KptI = aloc1(npres*m);
  double *KtpI = aloc1(npres*m);
  double *_fu  = aloc1(nne*nsd);
  
  memset(_fu, 0, sizeof(double)*nne*nsd);

  inverse(Ktp,m,KtpI);  
  inverse(Kpt,m,KptI);  

  double *KupKtpI, *KupKtpIKtt, *KptIFp;         
  KupKtpI    = aloc1(nne*nsd*m);
  KupKtpIKtt = aloc1(nne*nsd*m);
  KptIFp    = aloc1(m);

	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,m,m,1.0,Kup,m,KtpI,m,0.0,KupKtpI,m);
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,m,m,1.0,KupKtpI,m,Ktt,m,0.0,KupKtpIKtt,m);
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              m,1,m,1.0,KptI,m,fp,1,0.0,KptIFp,1);              

	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,1,m,1.0,KupKtpI,m,ft,1,0.0,_fu,1);     
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,1,m,-1.0,KupKtpIKtt,m,KptIFp,1,1.0,_fu,1);  
                               
  for(int a=0; a<nne*nsd; a++)
    fe[a] = fu[a] - _fu[a]; 

  free(KptI);
  free(KtpI);
  free(KupKtpI);
  free(KupKtpIKtt);
  free(KptIFp); 
  free(_fu);  
}

void condense_K2_to_K1(double *K11, int nne, int nsd, int npres,
                   double *Kuu, double *Kup,
                   double *Kpu, double *Kpp)                               
{
  double *KppI, *KupKppI, *_Kuu;
  KppI = aloc1(npres*npres);
  KupKppI = aloc1(nne*nsd*npres);
  _Kuu = aloc1(nne*nsd*nne*nsd);

  inverse(Kpp,npres,KppI);
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,npres,npres,1.0,Kup,npres,KppI,npres,0.0,KupKppI,npres);
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,nne*nsd,npres,1.0,KupKppI,npres,Kpu,nne*nsd,0.0,_Kuu,nne*nsd);
  for(long a=0;a<nne*nsd*nne*nsd; a++)
    K11[a] = Kuu[a] - _Kuu[a];
    
  free(KppI);    
  free(KupKppI);
  free(_Kuu);          
  
}
void condense_K3_to_K2(double *K11, double *K12, double *K21, double *K22, 
                   int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp)                               
{
	double *KttI, *KutKttI, *KptKttI;
  KttI = aloc1(nVol*nVol);
  KutKttI = aloc1(nne*nsd*nVol);
  KptKttI = aloc1(npres*nVol);
   
  inverse(Ktt,nVol,KttI);  
     
  double *_Kuu = aloc1(nne*nsd*nne*nsd);
  double *_Kup = aloc1(nne*nsd*npres);
  
  double *_Kpu = aloc1(npres*nsd*nne);
  double *_Kpp = aloc1(npres*npres);
   
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,nVol,nVol,1.0,Kut,nVol,KttI,nVol,0.0,KutKttI,nVol);
       
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              npres,nVol,nVol,1.0,Kpt,nVol,KttI,nVol,0.0,KptKttI,nVol);        
              
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,nne*nsd,nVol,1.0,KutKttI,nVol,Ktu,nne*nsd,0.0,_Kuu,nne*nsd);
       
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,npres,nVol,1.0,KutKttI,nVol,Ktp,npres,0.0,_Kup,npres);
        
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              npres,nne*nsd,nVol,1.0,KptKttI,nVol,Ktu,nne*nsd,0.0,_Kpu,nne*nsd);
       
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              npres,npres,nVol,1.0,KptKttI,nVol,Ktp,npres,0.0,_Kpp,npres);
  
	for(long a=0; a<nne*nsd*nne*nsd; a++)
	  K11[a] = Kuu[a] - _Kuu[a];
	
	for(long a=0; a<nne*nsd*npres; a++)
	{
	  K12[a] = Kup[a] - _Kup[a];
	  K21[a] = Kpu[a] - _Kpu[a];
	}
	  
	for(long a=0; a<npres*npres; a++)
	  K22[a] = Kpp[a] - _Kpp[a];	  

/*
// DEBUG 
  for(int a = 0; a<nVol; a++)
  {
    for(int b = 0; b<nVol; b++)
      printf("%e ", KttI[a*nVol+b]);
    printf("\n");  
  }
  */
  free(_Kuu);
  free(_Kup);
  free(KttI);
  
  free(_Kpu);
  free(_Kpp);
  
  free(KutKttI);
  free(KptKttI);  
  
  
} 

void condense_Kupt_to_Ku(double *Ks, int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp)                               
{
  /*
  double *K11 = aloc1(nne*nsd*nne*nsd);
  double *K12 = aloc1(nne*nsd*npres);
  
  double *K21 = aloc1(npres*nsd*nne);
  double *K22 = aloc1(npres*npres);

  condense_3_to_2(K11,K12,K21,K22,nne,nsd,npres,nVol,
                  Kuu,Kut,Kup,Ktu,Ktt,Ktp,Kpu,Kpt,Kpp); 
                  
// \/ CBD ---------------------------------------------------------------------
  int d1 = npres;
  int d2 = npres;
  for(int a = 0; a<d1; a++)
  {
    for(int b = 0; b<d2; b++)
      printf("%e ", K12[a*d2 + b]);
    printf("\n");
  }  
  printf("\n\n\n"); 
  
// /\ CBD ---------------------------------------------------------------------                   
                  
  condense_2_to_1(Ks,nne,nsd,npres,
                  K11,K12,K21,K22);                   

  free(K11);
  free(K12);
  free(K21);
  free(K22); */
  
  
	double *KptI, *KtpI;
	int m = nVol;
	
  KptI = aloc1(npres*m);
  KtpI = aloc1(npres*m);

  inverse(Ktp,m,KtpI);  
  inverse(Kpt,m,KptI);  

  double *KupKtpI, *KupKtpIKtt, *KptIKpu;         
  KupKtpI    = aloc1(nne*nsd*m);
  KupKtpIKtt = aloc1(nne*nsd*m);
  KptIKpu    = aloc1(nne*nsd*m);

	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,m,m,1.0,Kup,m,KtpI,m,0.0,KupKtpI,m);
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,m,m,1.0,KupKtpI,m,Ktt,m,0.0,KupKtpIKtt,m);
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              m,nne*nsd,m,1.0,KptI,m,Kpu,nne*nsd,0.0,KptIKpu,nne*nsd);              

	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nne*nsd,nne*nsd,m,1.0,KupKtpIKtt,m,KptIKpu,nne*nsd,1.0,Kuu,nne*nsd);                            
  
  for(int a=0; a<nne*nsd*nne*nsd; a++)
    Ks[a] = Kuu[a];
    
  free(KptI);
  free(KtpI);
  free(KupKtpI);
  free(KupKtpIKtt);
  free(KptIKpu);                  
}
  
void condense_Kupt_to_Kup(double *Ks, int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp)                               
{
  int ndofn = nsd + 1;
  double *K11 = aloc1(nne*nsd*nne*nsd);
  double *K12 = aloc1(nne*nsd*npres);
  
  double *K21 = aloc1(npres*nne*nsd);
  double *K22 = aloc1(npres*npres);
        
  condense_K3_to_K2(K11,K12,K21,K22,nne,nsd,npres,nVol,
                  Kuu,Kut,Kup,Ktu,Ktt,Ktp,Kpu,Kpt,Kpp);
                  
	for(long a = 0; a<nne; a++)
	{
	  for(long b=0; b<ndofn; b++)
	  {
	    for(long c=0; c<nne; c++)
	    {
	      for(long d = 0; d<=ndofn; d++)
	      {
	        int idx = idx_K(a,b,c,d,nne,ndofn);
	        if(b<nsd && d<nsd)
	        {
	          int idx_uu = idx_K(a,b,c,d,nne,nsd);  
	          Ks[idx] = K11[idx_uu];
	        }
	        else if(b<nsd && d==nsd)
	        {
	          int idx_up = idx_K_gen(a,b,c,0,nne,nsd,npres,1);
	          Ks[idx] = K12[idx_up];
	        }
	        else if(b==nsd && d<nsd)
	        {
	          int idx_pu = idx_K_gen(a,0,c,d,npres,1,nne,nsd);
	          Ks[idx] = K21[idx_pu];
	        }
	        else if(b==nsd && d==nsd)
	        {
	          int idx_pp = idx_K_gen(a,0,c,0,npres,1,npres,1);
	          Ks[idx] = K22[idx_pp];
	        }
	      }
	    }
	  }
	} 
  free(K11);
  free(K12);
  free(K21);
  free(K22);
} 

void condense_K_out(double *Ks, int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp)                               
{
  if(npres==4)
  {
    condense_Kupt_to_Kup(Ks,nne,nsd,npres,nVol,
                     Kuu,Kut,Kup,Ktu,Ktt,Ktp,Kpu,Kpt,Kpp);
    return;
  }
                    
  if(npres==1)
  {  
     condense_Kupt_to_Ku(Ks,nne,nsd,npres,nVol,
                     Kuu,Kut,Kup,Ktu,Ktt,Ktp,Kpu,Kpt,Kpp);
    return;
  }
  printf("codensing with pressure # %d is not supported\n", npres);
  exit(0);                                             
}
