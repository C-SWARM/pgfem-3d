#include "condense.h"

#include "data_structure.h"
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
	
//*/


void condense_Fupt_to_Fup(double *fe, int nne, int nsd, int npres, int nVol,
                   Matrix<double> &fu, Matrix<double> &ft, Matrix<double> &fp, 
                   Matrix<double> &Kut, Matrix<double> &Ktt, Matrix<double> &Kpt)
{
  int ndofn = nsd + 1; 
  Matrix<double> KttI, KutKttI, KptKttI, _fu, _fp;

     KttI.initialization(nVol,   nVol);  
	KutKttI.initialization(nne*nsd,nVol);
	KptKttI.initialization(npres,  nVol);
	    _fu.initialization(nne*nsd,1);
	    _fp.initialization(npres,  1);    
  
  KttI.inv(Ktt);
         
  // Matrix_AxB(C, a, b, A, AT, B, BT) <-- C = aAxB + bC
  KutKttI.prod(    Kut, KttI);
  KptKttI.prod(    Kpt, KttI);  
      _fu.prod(KutKttI,   ft);
      _fp.prod(KptKttI,   ft);  
                        
	for(long a = 0; a<nne; a++)
	{
	  for(long b=0; b<ndofn; b++)
	  {
       if(b<nsd)
         fe[a*ndofn + b] = fu.m_pdata[a*nsd+b] - _fu.m_pdata[a*nsd+b];
       else
         fe[a*ndofn + b] = fp.m_pdata[a] - _fp.m_pdata[a-1];
	  }
	}    
}

void condense_Fupt_to_Fu(double *fe, int nne, int nsd, int npres, int nVol,
                   Matrix<double> &fu, Matrix<double> &ft, Matrix<double> &fp, 
                   Matrix<double> &Kup, Matrix<double> &Ktp, Matrix<double> &Ktt, Matrix<double> &Kpt)                               
{  
	Matrix<double> KptI, KtpI, fu_add;

  KtpI.inv(Ktp);
  KptI.inv(Kpt);

  Matrix<double> KttKptIFp, ft_KttKptIFp;
  
  fu_add.prod(Kup, KtpI);
  KttKptIFp.prod(Ktt,KptI,fp);
  ft_KttKptIFp.sub(ft,KttKptIFp);

  fu_add.prod(ft_KttKptIFp);
  
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);	
  
  if(myrank==-1)
  { 
	  printf("-------------------------------------\n");	     
	  Kup.print("Kup");
    Ktp.print("Ktp");
	   ft.print("ft");
	  Ktt.print("Ktt"); 
	  Kpt.print("Kpt"); 
	   fp.print("fp"); 
	  fu_add.print("fu");	  	  
	  printf("-------------------------------------\n");	 
	   	  	  
	}  
              
  for(int a=0; a<nne*nsd; a++)
    fe[a] =  fu.m_pdata[a] - fu_add.m_pdata[a];
    
}

void condense_K2_to_K1(double *K11, int nne, int nsd, int npres,
                   Matrix<double> &Kuu, Matrix<double> &Kup, Matrix<double> &Kpu, Matrix<double> &Kpp)                               
{
  Matrix<double> KppI, KupKppI, _Kuu;
     KppI.initialization(npres,  npres  );
  KupKppI.initialization(nne*nsd,npres  );
     _Kuu.initialization(nne*nsd,nne*nsd);

  KppI.inv(Kpp);

  KupKppI.prod(Kup, KppI);
  _Kuu.prod(KupKppI,Kpu);  
  
  for(long a=0;a<nne*nsd*nne*nsd; a++)
    K11[a] = Kuu.m_pdata[a] - _Kuu.m_pdata[a];  
}

void condense_K3_to_K2(Matrix<double> &K11, Matrix<double> &K12, Matrix<double> &K21, Matrix<double> &K22, 
                   int nne, int nsd, int npres, int nVol,
                   Matrix<double> &Kuu, Matrix<double> &Kut, Matrix<double> &Kup,
                   Matrix<double> &Ktu, Matrix<double> &Ktt, Matrix<double> &Ktp,
                   Matrix<double> &Kpu, Matrix<double> &Kpt, Matrix<double> &Kpp)                               
{
	Matrix<double> KttI, KutKttI, KptKttI;
	
	   KttI.initialization(nVol,   nVol);
	KutKttI.initialization(nne*nsd,nVol);
	KptKttI.initialization(npres,  nVol);
	
  KttI.inv(Ktt);
   
  Matrix<double> _Kuu, _Kup, _Kpu, _Kpp;
	_Kuu.initialization(nne*nsd,nne*nsd);
	_Kup.initialization(nne*nsd,npres  );
	_Kpu.initialization(npres,  nne*nsd);  
  _Kpp.initialization(npres,  npres  );
  
  KutKttI.prod(Kut,     KttI);
  KptKttI.prod(Kpt,     KttI);
     _Kuu.prod(KutKttI, Ktu );
     _Kup.prod(KutKttI, Ktp );
     _Kpu.prod(KptKttI, Ktu );
     _Kpp.prod(KptKttI, Ktp );  
  
  K11.sub(Kuu,_Kuu);
  K12.sub(Kup,_Kup);
  K21.sub(Kpu,_Kpu);
  K22.sub(Kpp,_Kpp);        
} 

void condense_Kupt_to_Ku(double *Ks, int nne, int nsd, int npres, int nVol,
                   Matrix<double> &Kuu, Matrix<double> &Kut, Matrix<double> &Kup,
                   Matrix<double> &Ktu, Matrix<double> &Ktt, Matrix<double> &Ktp,
                   Matrix<double> &Kpu, Matrix<double> &Kpt, Matrix<double> &Kpp)
{

	Matrix<double> KptI, KtpI, Kuu_add;
	   KptI.initialization(npres,  nVol   );
	   KtpI.initialization(nVol,   npres  );
	Kuu_add.initialization(nne*nsd,nne*nsd);	
	
	KptI.inv(Kpt);
	KtpI.inv(Ktp);

  Matrix<double> KupKtpI, KupKtpIKtt, KptIKpu;
     KupKtpI.initialization(nne*nsd,nVol   );
  KupKtpIKtt.initialization(nne*nsd,nVol   );
     KptIKpu.initialization(nVol,   nne*nsd);
           
  // Matrix_AxB(C, a, b, A, AT, B, BT) <-- C = aAxB + bC
     KupKtpI.prod(Kup,KtpI);
  KupKtpIKtt.prod(KupKtpI,    Ktt    );  
     KptIKpu.prod(KptI,       Kpu    );  
     Kuu_add.prod(KupKtpIKtt, KptIKpu);    


  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  int myrank;
  MPI_Comm_rank(mpi_comm,&myrank);	
  
  if(myrank==-1)
  { 
	  printf("-------------------------------------\n");
	  Kuu.print("Kuu");	  
	  Kup.print("Kup");
	  Ktp.print("Ktp");
	  Ktt.print("Ktt");	  	  
	  Kpt.print("Kpt");	  	  
	  Kpu.print("Kpu");	  	  	  
	  Kuu_add.print("Kuu_add");	  	  	  	  
	  printf("-------------------------------------\n");
	}  

  for(int a=0; a<nne*nsd*nne*nsd; a++)
    Ks[a] = Kuu.m_pdata[a] + Kuu_add.m_pdata[a];    
}

void condense_Kupt_to_Kup(double *Ks, int nne, int nsd, int npres, int nVol,
                   Matrix<double> &Kuu, Matrix<double> &Kut, Matrix<double> &Kup,
                   Matrix<double> &Ktu, Matrix<double> &Ktt, Matrix<double> &Ktp,
                   Matrix<double> &Kpu, Matrix<double> &Kpt, Matrix<double> &Kpp)                               
{
  int ndofn = nsd + 1;
  Matrix<double> K11, K12, K21, K22;
  
  K11.initialization(nne*nsd,nne*nsd);
  K12.initialization(nne*nsd,npres  );
  K21.initialization(npres,  nne*nsd);
  K22.initialization(npres,  npres  );
          
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
} 

void condense_K_out(double *Ks, int nne, int nsd, int npres, int nVol,
                   double *Kuu, double *Kut, double *Kup,
                   double *Ktu, double *Ktt, double *Ktp,
                   double *Kpu, double *Kpt, double *Kpp)                               
{
  
  Matrix<double> Kuu_temp, Kut_temp, Kup_temp;
  Matrix<double> Ktu_temp, Ktt_temp, Ktp_temp;  
  Matrix<double> Kpu_temp, Kpt_temp, Kpp_temp;  
  
	Kuu_temp.use_reference(nne*nsd, nne*nsd, Kuu);
	Ktu_temp.use_reference(nVol,    nne*nsd, Ktu);
	Kpu_temp.use_reference(npres,   nne*nsd, Kpu);
	Kut_temp.use_reference(nne*nsd, nVol,    Kut);
	Ktt_temp.use_reference(nVol,    nVol,    Ktt);
	Kpt_temp.use_reference(npres,   nVol,    Kpt);	
	Kup_temp.use_reference(nne*nsd, npres,   Kup);
	Ktp_temp.use_reference(nVol,    npres,   Ktp);
	Kpp_temp.use_reference(npres,   npres,   Kpp);	
   
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
  Matrix<double> fu_temp,ft_temp,fp_temp;
	fu_temp.use_reference(nne*nsd, 1, fu);	
	ft_temp.use_reference(nVol,    1, ft);	
	fp_temp.use_reference(npres,   1, fp);	
	
  Matrix<double> Kut_temp,Kup_temp,Ktp_temp,Ktt_temp,Kpt_temp;
  Kut_temp.use_reference(nne*nsd, nVol,  Kut);	
  Kup_temp.use_reference(nne*nsd, npres, Kup);	
  Ktp_temp.use_reference(nVol,    npres, Ktp);
  Ktt_temp.use_reference(nVol,    nVol,  Ktt); 	
  Kpt_temp.use_reference(npres,   nVol,  Kpt);		
	
  if(npres==1)
    condense_Fupt_to_Fu(fe,nne,nsd,npres,nVol,
                   fu_temp,ft_temp,fp_temp,Kup_temp,Ktp_temp,Ktt_temp,Kpt_temp);     	
  if(npres==4)	
    condense_Fupt_to_Fup(fe,nne,nsd,npres,nVol, 
                   fu_temp,ft_temp,fp_temp,Kut_temp,Ktt_temp,Kpt_temp);
}