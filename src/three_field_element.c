#include "three_field_element.h"
#include "mkl_cblas.h"
#include "index_macros.h"
#include "utils.h"
#include "allocation.h"

void add_3F_Kuu_ip(double *K,
        int nne, double *ST, double *F, double jj, double wt, double Pn, double Tn,
        double dt, double alpha)
{
  int nsd = 3;
  double AA[9], BB[9], CC[9];
  
  double F_I[9];
  double Jn = det3x3(F);
  inverse(F,3,F_I);
  
  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const double* const ptrST_ab = &ST[idx_4_gen(a,b,0,0,
              nne,nsd,nsd,nsd)];
                    
      /* AA = F_I Grad(del u) */
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              3,3,3,1.0,F_I,3,ptrST_ab,3,0.0,AA,3);
      double trAA = 0.0;
      for(int m=0;m<nsd;m++)
        trAA += AA[m*nsd + m];
      
      for(int w=0; w<nne; w++)
      {
        for(int g=0; g<nsd; g++)
        {
          const double * const ptrST_wg = &ST[idx_4_gen(w,g,0,0,
                  nne,nsd,nsd,nsd)];
                                    
          /* BB = F_I Grad(d u) */
          cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                  3,3,3,1.0,F_I,3,ptrST_wg,3,0.0,BB,3);
          
          double trBB = 0.0;
          for(int m=0;m<nsd;m++)
            trBB += BB[m*nsd + m];
          
          /* CC = BB*AA */
          cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                  3,3,3,1.0,BB,3,AA,3,0.0,CC,3);
          
          double trCC = 0.0;
          for(int m=0;m<nsd;m++)
            trCC += CC[m*nsd + m];
          
          const int K_idx = idx_K(a,b,w,g,nne,nsd);
          K[K_idx] += -alpha*(1.0-alpha)*dt*jj*wt*Pn*Jn*(trAA*trBB - trCC);
          
        }
      }
    }
  }
}

void add_3F_Kut_ip(double *K,
        int nne, int nVol, double *ST, double *F, double jj, double wt, double *Nt,
        double dt, double alpha)
{
	int nsd = 3;
  return;
}
void add_3F_Kup_ip(double *K,
        int nne, int npres, double *ST, double *F, double jj, double wt, double *Np,
        double dt, double alpha)
{
	int nsd = 3;

  double *AA  = aloc1(9);
  double *C   = aloc1(9);
  double *C_I = aloc1(9);
  
  double Jn = det3x3(F);
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,F,3,F,3,0.0,C,3);
  
  inverse(C,3,C_I);
  
  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const double* const ptrST_ab = &ST[idx_4_gen(a,b,0,0,nne,nsd,nsd,nsd)];
      
      /* AA = F'Grad(del u) */
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              3,3,3,1.0,F,3,ptrST_ab,3,0.0,AA,3);
            
      /* C_IFdu = C_I:AA*/
      double C_IFdu = cblas_ddot(9,C_I,1,AA,1);
      
      for(int w=0; w<npres; w++)
      {
        int idx_up = idx_K_gen(a,b,w,0,nne,nsd,npres,1);
        K[idx_up] += -alpha*(1.0-alpha)*dt*jj*wt*Jn*C_IFdu*Np[w];
      }
    }
  }
  
  free(AA);
  free(C);
  free(C_I);
}

void add_3F_Ktu_ip(double *K,
        int nne, int nVol,
        double *ST, double *F, double jj, double wt, double *Nt,
        double dt, double alpha)
{
	int nsd = 3;
  return;
}

void add_3F_Ktt_ip(double *K, int nVol, double jj, double wt, double *Nt, double Upp,
        double dt, double alpha)
{
  for(int a=0; a<nVol; a++)
  {
    for(int b=0; b<nVol; b++)
      K[idx_K(a,0,b,0,nVol,1)] += -alpha*(1.0-alpha)*dt*Upp*(jj*wt*Nt[a]*Nt[b]);
  }
}


void add_3F_Ktp_ip(double *K,
        int nVol, int npres, double jj, double wt, double *Nt, double *Np,
        double dt, double alpha)
{

  for(int a=0; a<nVol; a++)
  {
    for(int b=0; b<npres; b++)
      K[idx_K_gen(a,0,b,0,nVol,1,npres,1)] += alpha*(1.0-alpha)*dt*jj*wt*Nt[a]*Np[b];
  }
}

void add_3F_Kpu_ip(double *K,
        int nne, int npres, double *ST, double *F, double jj, double wt, double *Np,
        double dt, double alpha)
{
	int nsd = 3;
  double *AA  = aloc1(9);
  double *C   = aloc1(9);
  double *C_I = aloc1(9);
  double *F_I = aloc1(9);
  
  double Jn = det3x3(F);
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,F,3,F,3,0.0,C,3);
  
  inverse(F,3,F_I);
  inverse(C,3,C_I);
  
  for(int a=0; a<npres; a++)
  {
    int b = 0;
    
    for(int w=0;w<nne; w++)
    {
      
      for(int g=0;g<nsd;g++)
      {
        const double * const ptrST_wg = &ST[idx_4_gen(w,g,0,0,
                nne,nsd,nsd,nsd)];
        
        /* AA = F'Grad(du) */
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                3,3,3,1.0,F,3,ptrST_wg,3,0.0,AA,3);
        
        double C_IFdu = 0.0;
        double F_I_Tdu = 0.0;
        /* C_IFdu = C_I:AA*/
        /* F_I_Tdu = F_I':du*/
        for(int i = 0; i<3; i++)
        {
          for(int j=0; j<3; j++)
          {
            C_IFdu += C_I[i*3+j]*AA[i*3+j];
            F_I_Tdu += F_I[j*3+i]*ptrST_wg[i*3+j];
          }
        }
        
        int idx_pu = idx_K_gen(a,b,w,g,npres,1,nne,nsd);
        K[idx_pu] += -alpha*(1.0-alpha)*dt*jj*wt*Jn*Np[a]*(F_I_Tdu*0.0 + C_IFdu);
      }
    }
  }
  
  free(AA);
  free(C);
  free(C_I);
  free(F_I);
}

void add_3F_Kpt_ip(double *K,
        int nVol, int npres, double jj, double wt, double *Nt, double *Np,
        double dt, double alpha)
{
		
  for(int a=0; a<npres; a++)
  {
    for(int b=0; b<nVol; b++)
      K[idx_K_gen(a,0,b,0,npres,1,nVol,1)] += alpha*(1.0-alpha)*dt*jj*wt*Np[a]*Nt[b];
  }
}

void resid_w_inertia_Ru_ip(double *fu,
        int nne, double *ST, double *F, double jj, double wt, double Pn)
{
	int nsd = 3;

  double *AA  = aloc1(9);
  double *C   = aloc1(9);
  double *C_I = aloc1(9);
  
  double Jn = det3x3(F);
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,F,3,F,3,0.0,C,3);
  
  inverse(C,3,C_I);
  
  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const double* const ptrST_ab = &ST[idx_4_gen(a,b,0,0,nne,nsd,nsd,nsd)];
      
      /* AA = F'Grad(del u) */
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              3,3,3,1.0,F,3,ptrST_ab,3,0.0,AA,3);
      
      /* C_IFdu = C_I:AA*/
      double C_IFdu = cblas_ddot(9,C_I,1,AA,1);      
            
      fu[a*nsd + b] += jj*wt*Jn*Pn*C_IFdu;
    }
  }
  
  free(AA);
  free(C);
  free(C_I);
}

void resid_w_inertia_Rt_ip(double *ft, int nVol, double jj, double wt, double *Nt, double Pn, double Up)
{

  for(int a=0; a<nVol; a++)
    ft[a] += jj*wt*Nt[a]*(Pn-Up);
}

void resid_w_inertia_Rp_ip(double *fp, int npres, double *F, double jj, double wt, double *Np, double Tn)
{
  double Jn = det3x3(F);
	for(int a=0; a<npres; a++)
		fp[a] += jj*wt*Np[a]*(Jn - Tn);
}    

