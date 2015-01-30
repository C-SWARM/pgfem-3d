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

#include "data_structure_c.h"

#ifndef _Matrix_double
Define_Matrix(double);
#define _Matrix_double 1
#endif


//*
// KttI:       [nVol]*[nVol]
// Kut:     [nne*nsd]*[nVol]
// KutKttI: [nne*nsd]*[nVol]
// Kpt:       [npres]*[nVol]
// KptKttI:   [npres]*[nVol]
	
//*/


void condense_Fupt_to_Fup(double *fe, int nne, int nsd, int npres, int nVol,
                   Matrix(double) fu, Matrix(double) ft, Matrix(double) fp, 
                   Matrix(double) Kut, Matrix(double) Ktt, Matrix(double) Kpt)
{
  int ndofn = nsd + 1; 
  Matrix(double) KttI, KutKttI, KptKttI, _fu, _fp;

  Matrix_construct_redim(double,KttI,   nVol,   nVol);  
	Matrix_construct_redim(double,KutKttI,nne*nsd,nVol);
	Matrix_construct_redim(double,KptKttI,npres,  nVol);
	Matrix_construct_redim(double,_fu,    nne*nsd,1);
	Matrix_construct_redim(double,_fp,    npres,  1);    
  
	Matrix_inv(Ktt, KttI);
         
  // Matrix_AxB(C, a, b, A, AT, B, BT) <-- C = aAxB + bC
  Matrix_AxB(KutKttI, 1.0, 0.0,     Kut, 0, KttI, 0);
  Matrix_AxB(KptKttI, 1.0, 0.0,     Kpt, 0, KttI, 0);  
  Matrix_AxB(    _fu, 1.0, 0.0, KutKttI, 0,   ft, 0);
  Matrix_AxB(    _fp, 1.0, 0.0, KptKttI, 0,   ft, 0);  
                        
	for(long a = 0; a<nne; a++)
	{
	  for(long b=0; b<ndofn; b++)
	  {
       if(b<nsd)
         fe[a*ndofn + b] = Mat_v(fu,a*nsd+b+1,1) - Mat_v(_fu, a*nsd+b+1, 1);
       else
         fe[a*ndofn + b] = Mat_v(fp,a+1,1) - Mat_v(_fp, a, 1);
	  }
	}  
  
  Matrix_cleanup(_fu);
  Matrix_cleanup(_fp);
  
  Matrix_cleanup(KttI);
  Matrix_cleanup(KutKttI);
  Matrix_cleanup(KptKttI);      
  
}

void condense_Fupt_to_Fu(double *fe, int nne, int nsd, int npres, int nVol,
                   Matrix(double) fu, Matrix(double) ft, Matrix(double) fp, 
                   Matrix(double) Kup, Matrix(double) Ktp, Matrix(double) Ktt, Matrix(double) Kpt)                               
{
  
	int m = nVol;
	
	Matrix(double) KptI, KtpI, _fu;
  Matrix_construct_redim(double,KptI,npres,nVol); 
	Matrix_construct_redim(double,KtpI,npres,nVol);
	Matrix_construct_init(double,_fu, nne*nsd,1,0.0);

  Matrix_inv(Ktp,KtpI);
  Matrix_inv(Kpt,KptI);

  Matrix(double) KupKtpI, KupKtpIKtt, KptIFp;
  Matrix_construct_redim(double,KupKtpI   ,nne*nsd,nVol); 
  Matrix_construct_redim(double,KupKtpIKtt,nne*nsd,nVol); 
  Matrix_construct_redim(double,KptIFp    ,nVol,   1   );     
  
  Matrix_AxB(KupKtpI,     1.0, 0.0, Kup,        0, KtpI,     0);
  Matrix_AxB(KupKtpIKtt,  1.0, 0.0, KupKtpI,    0, Ktt,      0);
  Matrix_AxB(KptIFp,      1.0, 0.0, KptI,       0, fp,       0);
  Matrix_AxB(_fu,         1.0, 0.0, KupKtpI,    0, ft,       0);
  Matrix_AxB(_fu,        -1.0, 1.0, KupKtpIKtt, 0, KptIFp,   0);  
              
  for(int a=0; a<nne*nsd; a++)
    fe[a] =  Mat_v(fu,a+1,1) -  Mat_v(_fu,a+1,1); 

  Matrix_cleanup(KptI);
  Matrix_cleanup(KtpI);
  Matrix_cleanup(KupKtpI);
  Matrix_cleanup(KupKtpIKtt);
  Matrix_cleanup(KptIFp); 
  Matrix_cleanup(_fu);  
}

void condense_K2_to_K1(double *K11, int nne, int nsd, int npres,
                   Matrix(double) Kuu, Matrix(double) Kup, Matrix(double) Kpu, Matrix(double) Kpp)                               
{
  Matrix(double) KppI, KupKppI, _Kuu;
  Matrix_construct_redim(double, KppI   ,npres,  npres  );
  Matrix_construct_redim(double, KupKppI,nne*nsd,npres  );
  Matrix_construct_redim(double, _Kuu   ,nne*nsd,nne*nsd);

  Matrix_inv(Kpp, KppI);

  Matrix_AxB(KupKppI, 1.0, 0.0, Kup,     0, KppI, 0);
  Matrix_AxB(_Kuu,    1.0, 0.0, KupKppI, 0, Kpu,  0);  
  
  for(long a=0;a<nne*nsd*nne*nsd; a++)
    K11[a] = Kuu.m_pdata[a] - _Kuu.m_pdata[a];
    
  Matrix_cleanup(KppI);    
  Matrix_cleanup(KupKppI);
  Matrix_cleanup(_Kuu);          
  
}

void condense_K3_to_K2(Matrix(double) K11, Matrix(double) K12, Matrix(double) K21, Matrix(double) K22, 
                   int nne, int nsd, int npres, int nVol,
                   Matrix(double) Kuu, Matrix(double) Kut, Matrix(double) Kup,
                   Matrix(double) Ktu, Matrix(double) Ktt, Matrix(double) Ktp,
                   Matrix(double) Kpu, Matrix(double) Kpt, Matrix(double) Kpp)                               
{
	Matrix(double) KttI, KutKttI, KptKttI;
	
	Matrix_construct_redim(double, KttI   ,nVol,   nVol);
	Matrix_construct_redim(double, KutKttI,nne*nsd,nVol);
	Matrix_construct_redim(double, KptKttI,npres,  nVol);
	
  Matrix_inv(Ktt, KttI);
   
  Matrix(double) _Kuu, _Kup, _Kpu, _Kpp;
	Matrix_construct_redim(double, _Kuu,nne*nsd,nne*nsd);
	Matrix_construct_redim(double, _Kup,nne*nsd,npres  );
	Matrix_construct_redim(double, _Kpu,npres,  nne*nsd);  
  Matrix_construct_redim(double, _Kpp,npres,  npres  );
  
  Matrix_AxB(KutKttI, 1.0, 0.0, Kut,     0, KttI, 0);
  Matrix_AxB(KptKttI, 1.0, 0.0, Kpt,     0, KttI, 0);
  Matrix_AxB(_Kuu,    1.0, 0.0, KutKttI, 0, Ktu,  0);
  Matrix_AxB(_Kup,    1.0, 0.0, KutKttI, 0, Ktp,  0);
  Matrix_AxB(_Kpu,    1.0, 0.0, KptKttI, 0, Ktu,  0);
  Matrix_AxB(_Kpp,    1.0, 0.0, KptKttI, 0, Ktp,  0);  
  
  // C = aA + bB
  // Matrix_AplusB(C, a, A, b, B)  
  Matrix_AplusB(K11, 1.0, Kuu, -1.0, _Kuu);
  Matrix_AplusB(K12, 1.0, Kup, -1.0, _Kup);
  Matrix_AplusB(K21, 1.0, Kpu, -1.0, _Kpu);
  Matrix_AplusB(K22, 1.0, Kpp, -1.0, _Kpp);  
        
  Matrix_cleanup(KttI);  
  Matrix_cleanup(KutKttI);
  Matrix_cleanup(KptKttI);    
  
  Matrix_cleanup(_Kuu);
  Matrix_cleanup(_Kup);
  Matrix_cleanup(_Kpu);
  Matrix_cleanup(_Kpp);  
} 

void condense_Kupt_to_Ku(double *Ks, int nne, int nsd, int npres, int nVol,
                   Matrix(double) Kuu, Matrix(double) Kut, Matrix(double) Kup,
                   Matrix(double) Ktu, Matrix(double) Ktt, Matrix(double) Ktp,
                   Matrix(double) Kpu, Matrix(double) Kpt, Matrix(double) Kpp)
{

	Matrix(double) KptI, KtpI, Kuu_add;
	Matrix_construct_redim(double, KptI   ,npres,  nVol   );
	Matrix_construct_redim(double, KtpI   ,nVol,   npres  );
	Matrix_construct_redim(double, Kuu_add,nne*nsd,nne*nsd);	
	
	Matrix_inv(Kpt, KptI);
	Matrix_inv(Ktp, KtpI);	
	

  Matrix(double) KupKtpI, KupKtpIKtt, KptIKpu;
  Matrix_construct_redim(double,KupKtpI,   nne*nsd,nVol   );
  Matrix_construct_redim(double,KupKtpIKtt,nne*nsd,nVol   );
  Matrix_construct_redim(double,KptIKpu,   nVol,   nne*nsd);
           
  // Matrix_AxB(C, a, b, A, AT, B, BT) <-- C = aAxB + bC
  Matrix_AxB(KupKtpI,    1.0, 0.0, Kup,        0, KtpI,    0);
  Matrix_AxB(KupKtpIKtt, 1.0, 0.0, KupKtpI,    0, Ktt,     0);  
  Matrix_AxB(KptIKpu,    1.0, 0.0, KptI,       0, Kpu,     0);  
  Matrix_AxB(Kuu_add,    1.0, 0.0, KupKtpIKtt, 0, KptIKpu, 0);    
               
  for(int a=0; a<nne*nsd*nne*nsd; a++)
    Ks[a] = Kuu.m_pdata[a] + Kuu_add.m_pdata[a];
    
  Matrix_cleanup(KptI);
  Matrix_cleanup(KtpI);
  Matrix_cleanup(Kuu_add);
  Matrix_cleanup(KupKtpI);
  Matrix_cleanup(KupKtpIKtt);
  Matrix_cleanup(KptIKpu);     
}

void condense_Kupt_to_Kup(double *Ks, int nne, int nsd, int npres, int nVol,
                   Matrix(double) Kuu, Matrix(double) Kut, Matrix(double) Kup,
                   Matrix(double) Ktu, Matrix(double) Ktt, Matrix(double) Ktp,
                   Matrix(double) Kpu, Matrix(double) Kpt, Matrix(double) Kpp)                               
{
  int ndofn = nsd + 1;
  Matrix(double) K11, K12, K21, K22;
  
  Matrix_construct_redim(double,K11,nne*nsd,nne*nsd);
  Matrix_construct_redim(double,K12,nne*nsd,npres  );
  Matrix_construct_redim(double,K21,npres,  nne*nsd);
  Matrix_construct_redim(double,K22,npres,  npres  );
          
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
	          Ks[idx] = K11.m_pdata[idx_uu];
	        }
	        else if(b<nsd && d==nsd)
	        {
	          int idx_up = idx_K_gen(a,b,c,0,nne,nsd,npres,1);
	          Ks[idx] = K12.m_pdata[idx_up];
	        }
	        else if(b==nsd && d<nsd)
	        {
	          int idx_pu = idx_K_gen(a,0,c,d,npres,1,nne,nsd);
	          Ks[idx] = K21.m_pdata[idx_pu];
	        }
	        else if(b==nsd && d==nsd)
	        {
	          int idx_pp = idx_K_gen(a,0,c,0,npres,1,npres,1);
	          Ks[idx] = K22.m_pdata[idx_pp];
	        }
	      }
	    }
	  }
	} 
  Matrix_cleanup(K11);
  Matrix_cleanup(K12);
  Matrix_cleanup(K21);
  Matrix_cleanup(K22);
} 

void condense_K_out(double *Ks, int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp)                               
{
  
  Matrix(double) Kuu_temp, Kut_temp, Kup_temp;
  Matrix(double) Ktu_temp, Ktt_temp, Ktp_temp;  
  Matrix(double) Kpu_temp, Kpt_temp, Kpp_temp;  
  
	Matrix_construct(double, Kuu_temp);
	Matrix_construct(double, Ktu_temp);
	Matrix_construct(double, Kpu_temp);

	Matrix_construct(double, Kut_temp);
	Matrix_construct(double, Ktt_temp);
	Matrix_construct(double, Kpt_temp);
	
  Matrix_construct(double, Kup_temp);
	Matrix_construct(double, Ktp_temp);
	Matrix_construct(double, Kpp_temp);	
                            
	Matrix_init_w_array(Kuu_temp, nne*nsd, nne*nsd, Kuu);
	Matrix_init_w_array(Ktu_temp, nVol,    nne*nsd, Ktu);
	Matrix_init_w_array(Kpu_temp, npres,   nne*nsd, Kpu);
	Matrix_init_w_array(Kut_temp, nne*nsd, nVol,    Kut);
	Matrix_init_w_array(Ktt_temp, nVol,    nVol,    Ktt);
	Matrix_init_w_array(Kpt_temp, npres,   nVol,    Kpt);	
	Matrix_init_w_array(Kup_temp, nne*nsd, npres,   Kup);
	Matrix_init_w_array(Ktp_temp, nVol,    npres,   Ktp);
	Matrix_init_w_array(Kpp_temp, npres,   npres,   Kpp);	
   
  switch(npres)
  {
    case 1:
      condense_Kupt_to_Ku(Ks,nne,nsd,npres,nVol,
                     Kuu_temp,Kut_temp,Kup_temp,
                     Ktu_temp,Ktt_temp,Ktp_temp,
                     Kpu_temp,Kpt_temp,Kpp_temp);
      break;
    case 4:                  
      condense_Kupt_to_Kup(Ks,nne,nsd,npres,nVol,
                     Kuu_temp,Kut_temp,Kup_temp,
                     Ktu_temp,Ktt_temp,Ktp_temp,
                     Kpu_temp,Kpt_temp,Kpp_temp);
      break;                     
  }
  
  Matrix_cleanup(Kuu_temp);
  Matrix_cleanup(Ktu_temp);
  Matrix_cleanup(Kpu_temp);
  Matrix_cleanup(Kut_temp);
  Matrix_cleanup(Ktt_temp);
  Matrix_cleanup(Kpt_temp);
  Matrix_cleanup(Kup_temp);
  Matrix_cleanup(Ktp_temp);
  Matrix_cleanup(Kpp_temp);
  
  if(npres!=1 && npres!=4)
  {    
    printf("codensing with pressure # %d is not supported\n", npres);
    exit(0);
  }                                             
}

void condense_F_out(double *fe, int nne, int nsd, int npres, int nVol,
                   double *fu, double *ft, double *fp, 
                   double *Kut, double *Kup, double *Ktp, double *Ktt,double *Kpt)
{
  Matrix(double) fu_temp,ft_temp,fp_temp;
	Matrix_construct(double, fu_temp);	  
	Matrix_construct(double, ft_temp);
	Matrix_construct(double, fp_temp);
	
	Matrix_init_w_array(fu_temp, nne*nsd, 1, fu);	
	Matrix_init_w_array(ft_temp, nVol,    1, ft);	
	Matrix_init_w_array(fp_temp, npres,   1, fp);	
	
  Matrix(double) Kut_temp,Kup_temp,Ktp_temp,Ktt_temp,Kpt_temp;
	Matrix_construct(double, Kut_temp);	
	Matrix_construct(double, Kup_temp);		
	Matrix_construct(double, Ktp_temp);	
	Matrix_construct(double, Ktt_temp);		
	Matrix_construct(double, Kpt_temp);

	Matrix_init_w_array(Kut_temp, nne*nsd, nVol,  Kut);	
	Matrix_init_w_array(Kup_temp, nne*nsd, npres, Kup);	
	Matrix_init_w_array(Ktp_temp, nVol,    npres, Ktp);
	Matrix_init_w_array(Ktt_temp, nVol,    nVol,  Ktt); 	
	Matrix_init_w_array(Kpt_temp, npres,   nVol,  Kpt);		
	
  if(npres==1)
    condense_Fupt_to_Fu(fe,nne,nsd,npres,nVol,
                   fu_temp,ft_temp,fp_temp,Kup_temp,Ktp_temp,Ktt_temp,Kpt_temp);     	
  if(npres==4)	
    condense_Fupt_to_Fup(fe,nne,nsd,npres,nVol, 
                   fu_temp,ft_temp,fp_temp,Kut_temp,Ktt_temp,Kpt_temp);

  Matrix_cleanup(fu_temp);
  Matrix_cleanup(ft_temp);
  Matrix_cleanup(fp_temp);

  Matrix_cleanup(Kut_temp);
  Matrix_cleanup(Kup_temp);
  Matrix_cleanup(Ktp_temp);  
  Matrix_cleanup(Ktt_temp);
  Matrix_cleanup(Kpt_temp);      
}