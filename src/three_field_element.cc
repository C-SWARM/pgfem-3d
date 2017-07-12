#include "three_field_element.h"
#include "mkl_cblas.h"
#include "index_macros.h"
#include "utils.h"
#include "allocation.h"
#include "femlib.h"

#include "Hu_Washizu_element.h"
#include "displacement_based_element.h"

#include "condense.h"
#include "cast_macros.h"
#include "def_grad.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "enumerations.h"
#include "dynamics.h"

#define ndn 3
#define USE_HW_FUNCS 0
#define INTG_ORDER 1

void add_3F_Kuu_ip_disp(double *K, FEMLIB *fe, 
        double *ST, Matrix<double> &F, Matrix<double> &S, double *L,
        double dt_alpha_1_minus_alpha)
{
  int nne = fe->nne;
  int nsd = fe->nsd;
  Matrix<double> ST_ab, ST_wg, AA, BB, CC;
  Matrix<double> sAA, sBB, sCC;
  AA.initialization(3,3);
  BB.initialization(3,3);
  CC.initialization(3,3);
  sAA.initialization(3,3);
  sBB.initialization(3,3);
  sCC.initialization(3,3);
  
  
  double LsBB[9];
    
  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const double* const ptrST_ab = &ST[idx_4_gen(a,b,0,0,
              nne,nsd,nsd,nsd)];
      ST_ab.use_reference(3,3,ptrST_ab);

      // AA = F^T Grad(del u)
      AA.prod(F,1,ST_ab,0); 
      symmetric_part(sAA.m_pdata,AA.m_pdata,3);
                          
      for(int w=0; w<nne; w++)
      {
        for(int g=0; g<nsd; g++)
        {
          const double * const ptrST_wg = &ST[idx_4_gen(w,g,0,0,
                  nne,nsd,nsd,nsd)];
          ST_wg.use_reference(3,3,ptrST_wg);
                                    
          // BB = F^T Grad(du) 
          BB.prod(F,1,ST_wg,0);
          symmetric_part(sBB.m_pdata,BB.m_pdata,3);
          
          // CC = Grad(del u)^T Grad(du) 
          CC.prod(ST_ab,1,ST_wg,0);
          symmetric_part(sCC.m_pdata,CC.m_pdata,3);
          
          // Compute L:sBB 

          for (int i=0; i<nsd; i++){
            for (int j=0; j<nsd; j++){
              LsBB[idx_2(i,j)] = cblas_ddot(9,&L[idx_4(i,j,0,0)],1,
          			      sBB.m_pdata,1);
            }
          }            
          const int K_idx = idx_K(a,b,w,g,nne,nsd);
          K[K_idx] += -dt_alpha_1_minus_alpha*fe->detJxW*(cblas_ddot(9,sAA.m_pdata,1,LsBB,1) +
                                                          cblas_ddot(9,sCC.m_pdata,1,S.m_pdata,1));          
        }
      }
    }
  }    
}

void add_3F_Kuu_ip(double *K,
        int nne, double *ST, double *F, double jj, double wt, double Pn,
        double dt_alpha_1_minus_alpha)
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
          K[K_idx] += -dt_alpha_1_minus_alpha*jj*wt*Pn*Jn*(trAA*trBB - trCC);
          
        }
      }
    }
  }
}
void add_3F_Kup_ip(double *K,
        int nne, int npres, double *ST, double *F, double jj, double wt, double *Np,
        double dt_alpha_1_minus_alpha)
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
        K[idx_up] += -dt_alpha_1_minus_alpha*jj*wt*Jn*C_IFdu*Np[w];
      }
    }
  }
  
  free(AA);
  free(C);
  free(C_I);
}

void add_3F_Ktt_ip(double *K, int nVol, double jj, double wt, double *Nt, double Upp,
        double dt_alpha_1_minus_alpha)
{
  for(int a=0; a<nVol; a++)
  {
    for(int b=0; b<nVol; b++)
      K[idx_K(a,0,b,0,nVol,1)] += -dt_alpha_1_minus_alpha*Upp*(jj*wt*Nt[a]*Nt[b]);
  }
}
void add_3F_Ktp_ip(double *K,
        int nVol, int npres, double jj, double wt, double *Nt, double *Np,
        double dt_alpha_1_minus_alpha)
{

  for(int a=0; a<nVol; a++)
  {
    for(int b=0; b<npres; b++)
      K[idx_K_gen(a,0,b,0,nVol,1,npres,1)] += dt_alpha_1_minus_alpha*jj*wt*Nt[a]*Np[b];
  }
}
void add_3F_Kpu_ip(double *K,
        int nne, int npres, double *ST, double *F, double jj, double wt, double *Np,
        double dt_alpha_1_minus_alpha)
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
        K[idx_pu] += -dt_alpha_1_minus_alpha*jj*wt*Jn*Np[a]*(F_I_Tdu*0.0 + C_IFdu);
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
        double dt_alpha_1_minus_alpha)
{
		
  for(int a=0; a<npres; a++)
  {
    for(int b=0; b<nVol; b++)
      K[idx_K_gen(a,0,b,0,npres,1,nVol,1)] += dt_alpha_1_minus_alpha*jj*wt*Np[a]*Nt[b];
  }
}
void resid_w_inertia_Ru_ip(double *fu,
        int nne, double *ST, double *F, double *S, double jj, double wt, double Pn)
{
	int nsd = 3;

  double *AA  = aloc1(9);
  double *sAA = aloc1(9);
  double *C   = aloc1(9);
  double *C_I = aloc1(9);
  double *S_temp = aloc1(9);
  
  double Jn = det3x3(F);
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
          3,3,3,1.0,F,3,F,3,0.0,C,3);
  
  inverse(C,3,C_I);
//  for(int a = 0; a<9; a++)
//    S_temp[a] = Jn*Pn*C_I[a] + S[a];
//  cblas_daxpy(9,Jn*Pn,C_I,1,S,1);
  
  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const double* const ptrST_ab = &ST[idx_4_gen(a,b,0,0,nne,nsd,nsd,nsd)];
      
      // AA = F'Grad(del u);
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              3,3,3,1.0,F,3,ptrST_ab,3,0.0,AA,3);
      symmetric_part(sAA,AA,3);
      //sAA:S
      double sAAS = cblas_ddot(9,sAA,1,S,1);            
      //C_IFdu = C_I:AA
      double C_IFdu = cblas_ddot(9,C_I,1,AA,1);      
            
      fu[a*nsd + b] += jj*wt*(sAAS + Jn*Pn*C_IFdu);
//fu[a*nsd + b] += cblas_ddot(9,S_temp,1,AA,1)*jj*wt;
      
    }
  }
  
  free(AA);
  free(sAA);
  free(C);
  free(C_I);
  free(S_temp);
}
void resid_w_inertia_Rt_ip(double *ft, int nVol, double jj, double wt, double *Nt, double Pn, double Up)
{

  for(int a=0; a<nVol; a++)
    ft[a] += -jj*wt*Nt[a]*(Pn-Up);
}
void resid_w_inertia_Rp_ip(double *fp, int npres, double *F, double jj, double wt, double *Np, double Tn)
{
  double Jn = det3x3(F);
 // printf("%e, %e, %e\n", Jn, Tn, Jn - Tn); 
	for(int a=0; a<npres; a++)
		fp[a] += jj*wt*Np[a]*(Jn - Tn);
}    
void TF_Ru_ip(Matrix<double> &f, FEMLIB *fe, 
        Matrix<double> &F, Matrix<double> &S, double Pn)
{
  int nne = fe->nne;
  if(USE_HW_FUNCS)
  {
    double Fn[9];
    Fn[0] = Fn[4] = Fn[8] = 1.0;
    Fn[1] = Fn[2] = Fn[3] = Fn[5] = Fn[6] = Fn[7] = 0.0; 
    double F_I[9];
    inverse(F.m_pdata,3,F_I);       
    HW_Ru_at_ip(f.m_pdata,nne,nne,fe->ST,Fn,F.m_pdata,F_I,
    		1.0,0.0,1.0,1.0,S.m_pdata,Pn,fe->detJ,fe->temp_v.w_ip,0);
  }  		
  else
    resid_w_inertia_Ru_ip(f.m_pdata,nne,fe->ST,F.m_pdata,S.m_pdata,fe->detJ,fe->temp_v.w_ip, Pn);
}
void TF_Rt_ip(Matrix<double> &f, FEMLIB *fe, 
        Matrix<double> &F, int nVol, Matrix<double> &Nt, double Pn, double Up)
{
  if(USE_HW_FUNCS)
    HW_Rt_at_ip(f.m_pdata,nVol,Nt.m_pdata,1.0,1.0,Pn,1.0,Up,fe->detJ,fe->temp_v.w_ip);
  else
    resid_w_inertia_Rt_ip(f.m_pdata,nVol,fe->detJ,fe->temp_v.w_ip, Nt.m_pdata, Pn, Up);
}
void TF_Rp_ip(Matrix<double> &f, FEMLIB *fe, 
        Matrix<double> &F, int npres, Matrix<double> &Np, double Tn)
{
  if(USE_HW_FUNCS)
    HW_Rp_at_ip(f.m_pdata,npres,Np.m_pdata,1.0,1.0,1.0,1.0,fe->detJ,fe->temp_v.w_ip);
  else
    resid_w_inertia_Rp_ip(f.m_pdata, npres, F.m_pdata, fe->detJ,fe->temp_v.w_ip, Np.m_pdata, Tn);
}
void TF_Kuu_ip(Matrix<double> &K, FEMLIB *fe, 
        Matrix<double> &F, Matrix<double> &S, double *L, 
        double Pn, double dt_alpha_1_minus_alpha)
{
  int nne = fe->nne;
  int nsd = fe->nsd;
  
  if(USE_HW_FUNCS){
    double Fn[9];
    Fn[0] = Fn[4] = Fn[8] = 1.0;
    Fn[1] = Fn[2] = Fn[3] = Fn[5] = Fn[6] = Fn[7] = 0.0; 
    double F_I[9];
    inverse(F.m_pdata,3,F_I);             
    HW_Kuu_at_ip(K.m_pdata,nne,nne,fe->ST,Fn,F.m_pdata,F_I,1.0,1.0,1.0,Pn,S.m_pdata,L,fe->detJ,fe->temp_v.w_ip,0);
  }
  else{
    add_3F_Kuu_ip(K.m_pdata,nne ,fe->ST,F.m_pdata,fe->detJ,fe->temp_v.w_ip,Pn,dt_alpha_1_minus_alpha);  
    add_3F_Kuu_ip_disp(K.m_pdata, fe, fe->ST, F, S, L, dt_alpha_1_minus_alpha);    
  }    
}
void TF_Kup_ip(Matrix<double> &K, FEMLIB *fe, 
        Matrix<double> &F, int npres, Matrix<double> &Np, double dt_alpha_1_minus_alpha)
{
  int nne = fe->nne;
  if(USE_HW_FUNCS)
    HW_Kup_at_ip(K.m_pdata,nne,nne,npres,Np.m_pdata,fe->ST,F.m_pdata,1.0,1.0,1.0,fe->detJ,fe->temp_v.w_ip,0);
  else
    add_3F_Kup_ip(K.m_pdata,nne,npres,fe->ST,F.m_pdata,fe->detJ,fe->temp_v.w_ip,Np.m_pdata,dt_alpha_1_minus_alpha);    
}
void TF_Kut_ip(Matrix<double> &K, FEMLIB *fe, 
        Matrix<double> &F, int nVol, Matrix<double> &Nt, double dt_alpha_1_minus_alpha)
{
  return;
}
void TF_Kpt_ip(Matrix<double> &K, FEMLIB *fe, 
        Matrix<double> &F, int npres, int nVol, Matrix<double> &Np, Matrix<double> &Nt,
        double dt_alpha_1_minus_alpha)
{
  int nne = fe->nne;
  if(USE_HW_FUNCS)
    HW_Kpt_at_ip(K.m_pdata,npres,Np.m_pdata,nVol,Nt.m_pdata,1.0,1.0,fe->detJ,fe->temp_v.w_ip);
  else
    add_3F_Kpt_ip(K.m_pdata,nVol,npres,fe->detJ,fe->temp_v.w_ip,Nt.m_pdata,Np.m_pdata,dt_alpha_1_minus_alpha);    
}
void TF_Ktt_ip(Matrix<double> &K, FEMLIB *fe, 
        Matrix<double> &F, int nVol, Matrix<double> &Nt, double Upp,
        double dt_alpha_1_minus_alpha)
{
  int nne = fe->nne;
  if(USE_HW_FUNCS)
    HW_Ktt_at_ip(K.m_pdata,nVol,Nt.m_pdata,1.0,1.0,1.0,Upp,fe->detJ,fe->temp_v.w_ip);             
  else
    add_3F_Ktt_ip(K.m_pdata,nVol,fe->detJ,fe->temp_v.w_ip,Nt.m_pdata,Upp,dt_alpha_1_minus_alpha); 
}

void stiffmat_3f_el(double *Ks,
        const int ii,
        const int ndofn,
        const int nne,
        int npres,
        int nVol,
        int nsd,
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const long *nod,
        const NODE *node,
        double dt,
        SIG *sig,
        EPS *eps,
        const SUPP sup,
        double alpha, double *r_e)
{
  double alpha_1;
  double alpha_2;
  double dt_alpha_1_minus_alpha;

  if(alpha<0)
  {  
    alpha = 1.0;
    alpha_1 = 0.0;
    alpha_2 = alpha;
    dt_alpha_1_minus_alpha = -1.0; 
  }
  else
  {
    alpha_1 = 1.0 - alpha;
    alpha_2 = alpha;
    dt_alpha_1_minus_alpha = dt*alpha_1*alpha_2;     
  } 
   
  double *u, *P;
  u = aloc1(nne*nsd);        
  P = aloc1(npres);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
    	u[a*nsd+b] = r_e[a*ndofn+b];
    	
    if(npres==nne)
      P[a] = r_e[a*ndofn+nsd];
  }
      
  if(npres==1)
    P[0] = alpha_1*eps[ii].d_T[1] + alpha_2*eps[ii].d_T[0]; 
     
  int ndofe = nne*ndofn;
  memset(Ks,0,ndofe*ndofe*sizeof(double));
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));  
  double rho = hommat[mat].density;
  
  dUdJFuncPtr          UP   = getDUdJFunc(1, &hommat[mat]);        
  d2UdJ2FuncPtr       UPP   = getD2UdJ2Func(1, &hommat[mat]);  
  devStressFuncPtr Stress   = getDevStressFunc(1,&hommat[mat]);
  matStiffFuncPtr Stiffness = getMatStiffFunc(1,&hommat[mat]);

  int err = 0;
  
  double *L = aloc1(81);
  Matrix<double> F,C,S;
  Matrix<double> Kuu,Kut,Kup;
  Matrix<double> Ktu,Ktt,Ktp;
  Matrix<double> Kpu,Kpt,Kpp;

  F.initialization(3,3,0.0);
  C.initialization(3,3,0.0);
  S.initialization(3,3,0.0);
  
  Kuu.initialization(nne*nsd,nne*nsd,0.0);
  Kut.initialization(nne*nsd,nVol   ,0.0);
  Kup.initialization(nne*nsd,npres  ,0.0);
  Ktu.initialization(nVol   ,nne*nsd,0.0);
  Ktt.initialization(nVol   ,nVol   ,0.0);
  Ktp.initialization(nVol   ,npres  ,0.0);    
  Kpu.initialization(npres  ,nne*nsd,0.0);
  Kpt.initialization(npres  ,nVol   ,0.0);
  Kpp.initialization(npres  ,npres  ,0.0);  


  FEMLIB fe(ii, elem, node, INTG_ORDER,1);
  
  Matrix<double> Np, Nt;
  Np.initialization(npres,1,0.0);
  Nt.initialization(nVol, 1,0.0);    
                         
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    fe.elem_basis_V(ip);  
    fe.elem_shape_function(ip,npres,Np.m_pdata);
    fe.elem_shape_function(ip,nVol, Nt.m_pdata);
    fe.update_shape_tensor();
    fe.update_deformation_gradient(ndofn,u,F.m_pdata);    
              
    double Tn = 0.0; 
    double Pn = 0.0;

    for(int a=0; a<nVol; a++)
      Tn += Nt(a+1)*(alpha_1*eps[ii].T[a*3+1] + alpha_2*eps[ii].T[a*3+0]);
       
    for(int a=0; a<npres; a++)
      Pn += Np(a+1)*P[a];  
      
    double Jn = det3x3(F.m_pdata);   
    
    C.prod(F,1,F,0);  
    double Upp = 0.0;          
    UPP(Tn,&hommat[mat],&Upp);
    Upp = Upp*kappa;
    double Up = 0.0;
    UP(Tn,&hommat[mat],&Up);
    Up = Up*kappa;
    
    Stress(C.m_pdata,&hommat[mat],S.m_pdata);       
    Stiffness(C.m_pdata,&hommat[mat],L);            

    TF_Kuu_ip(Kuu,&fe,F,S,L,Pn,dt_alpha_1_minus_alpha);
    TF_Kup_ip(Kup,&fe,F,npres,Np,dt_alpha_1_minus_alpha);     
    TF_Kpt_ip(Kpt,&fe,F,npres,nVol,Np,Nt,dt_alpha_1_minus_alpha);
    TF_Ktt_ip(Ktt,&fe,F,nVol,Nt,Upp,dt_alpha_1_minus_alpha);
  } 

  Kpu.trans(Kup);
  Ktp.trans(Kpt);
                
  condense_K_out(Ks,nne,nsd,npres,nVol,
                 Kuu.m_pdata,Kut.m_pdata,Kup.m_pdata,
                 Ktu.m_pdata,Ktt.m_pdata,Ktp.m_pdata,
                 Kpu.m_pdata,Kpt.m_pdata,Kpp.m_pdata);
    
  free(P); free(u);     
}

void residuals_3f_el(double *f,
        const int ii,
        const int ndofn,
        const int nne,
        const int npres,
        const int nVol,
        const int nsd,
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const long *nod,
        const NODE *node,
        double dt,
        SIG *sig,
        EPS *eps,
        const SUPP sup,
        double *r_e)
{
  double *u, *P;
  u = aloc1(nne*nsd);
  P = aloc1(npres);		

  for(int a=0;a<nne;a++)
  {
  	for(int b=0; b<nsd;b++)
  		u[a*nsd+b] = r_e[a*ndofn+b];	

  	if(npres==nne)
  	  P[a] = r_e[a*ndofn+nsd];
  }	
    					
  if(npres==1)
    P[0] =  eps[ii].d_T[0]; 
                        			  
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.0*(1.0-2.0*hommat[mat].nu));
      
  dUdJFuncPtr UP = getDUdJFunc(1, &hommat[mat]);
  d2UdJ2FuncPtr UPP = getD2UdJ2Func(1, &hommat[mat]);
  devStressFuncPtr Stress = getDevStressFunc(1,&hommat[mat]);

  Matrix<double> F,S,C;
  F.initialization(3,3,0.0);
  C.initialization(3,3,0.0);  
  S.initialization(3,3,0.0);  
  
  Matrix<double> fu,Kut;
  Matrix<double> fp,Kpt;
  Matrix<double> ft,Ktt;
  Matrix<double> Kup,Ktp;
  
  int ndofe = nne*ndofn;

   fu.initialization( nne*nsd,1,    0.0);
   fp.initialization( npres,  1,    0.0);
   ft.initialization( nVol,   1,    0.0);
  Kut.initialization(nne*nsd,nVol, 0.0);
  Kpt.initialization(npres,  nVol, 0.0);
  Ktt.initialization(nVol,   nVol, 0.0);    
  Kup.initialization(nne*nsd,npres,0.0);
  Ktp.initialization(nVol,   npres,0.0);
  
  FEMLIB fe(ii, elem, node, INTG_ORDER,1);
                          
  Matrix<double> Np, Nt;
  Np.initialization(npres,1,0.0);
  Nt.initialization(nVol, 1,0.0);    
                         
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    fe.elem_basis_V(ip);  
    fe.elem_shape_function(ip,npres,Np.m_pdata);
    fe.elem_shape_function(ip,nVol, Nt.m_pdata);
    fe.update_shape_tensor();
    fe.update_deformation_gradient(ndofn,u,F.m_pdata);
                  
    double Tn = 0.0; 
    double Pn = 0.0;
    
    for(int a=0; a<nVol; a++)
      Tn += Nt(a+1)*eps[ii].T[a*3+0];
    
    for(int a=0; a<npres; a++)
      Pn += Np(a+1)*P[a];
      
                
    double Upp = 0.0;
    double Up = 0.0;
    UPP(Tn,&hommat[mat],&Upp);
    Upp = Upp*kappa;

    UP(Tn, &hommat[mat], &Up);
    Up = Up*kappa;    
    
    C.prod(F,1,F,0);
    Stress(C.m_pdata,&hommat[mat],S.m_pdata);
    		
    TF_Kpt_ip(Kpt,&fe,F,npres,nVol,Np,Nt,-1.0);
    TF_Ktt_ip(Ktt,&fe,F,nVol,Nt,Upp,-1.0);        

    if(npres==nne)
      TF_Kut_ip(Kut,&fe,F,nVol,Nt,-1.0);
      
		if(npres==1)
      TF_Kup_ip(Kup,&fe,F,npres,Np,-1.0);
      
    TF_Ru_ip(fu,&fe,F,S,Pn);
    TF_Rp_ip(fp,&fe,F,npres,Np,Tn);
    TF_Rt_ip(ft,&fe,F,nVol,Nt,Pn,Up);
  }

  Ktp.trans(Kpt);

  condense_F_out(f,nne,nsd,npres,nVol,fu.m_pdata,ft.m_pdata,fp.m_pdata,
                    Kut.m_pdata,Kup.m_pdata,Ktp.m_pdata,Ktt.m_pdata,Kpt.m_pdata);
  
  free(P); free(u);   
}

void residuals_3f_w_inertia_el(double *f,
        const int ii,
        const int ndofn,
        const int nne,
        const int npres,
        const int nVol,
        const int nsd,
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const NODE *node,
        const double *dts,
        SIG *sig,
        EPS *eps,
        double alpha, double *r_n_a, double *r_n_1_a)
{
  double alpha_1;
  double alpha_2;
  double dt_alpha_1_minus_alpha;
  double dt_alpha_1;
  double dt_alpha_2;

  if(alpha<0)
  {  
    alpha = 1.0;
    alpha_1 = 0.0;
    alpha_2 = alpha;
    dt_alpha_1_minus_alpha = -1.0; 
    dt_alpha_1 = 1.0;
    dt_alpha_2 = 0.0;
  }
  else
  {
    alpha_1 = 1.0 - alpha;
    alpha_2 = alpha;
    dt_alpha_1_minus_alpha = dts[DT_NP1]*alpha_1*alpha_2; 
    dt_alpha_1 = -dts[DT_NP1]*alpha_1;  
    dt_alpha_2 = -dts[DT_N]*alpha_2;  
  } 
  
  double *u1 = aloc1(nne*nsd);
  double *u2 = aloc1(nne*nsd);  
  double *P1 = aloc1(npres);		
  double *P2 = aloc1(npres);		  

  for(int a=0;a<nne;a++)
  {
  	for(int b=0; b<nsd;b++)
  	{
  		u1[a*nsd+b] = r_n_1_a[a*ndofn+b];	
  		u2[a*nsd+b] = r_n_a[a*ndofn+b];	  		
  	}

  	if(npres==nne)
  	{
  	  P1[a] = r_n_1_a[a*ndofn+nsd];
  	  P2[a] = r_n_a[a*ndofn+nsd];  	  
  	}
  }	
    					
  if(npres==1)
  {  
    P1[0] = alpha_1*eps[ii].d_T[2] + alpha_2*eps[ii].d_T[1]; 
    P2[0] = alpha_1*eps[ii].d_T[1] + alpha_2*eps[ii].d_T[0];
  }                      			  
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.0*(1.0-2.0*hommat[mat].nu));
      
  dUdJFuncPtr          UP = getDUdJFunc(1, &hommat[mat]);
  d2UdJ2FuncPtr       UPP = getD2UdJ2Func(1, &hommat[mat]);
  devStressFuncPtr Stress = getDevStressFunc(1,&hommat[mat]);

  Matrix<double> F1,F2,C1,C2,S1,S2;
  F1.initialization(3,3,0.0);
  F2.initialization(3,3,0.0);  
  C1.initialization(3,3,0.0);
  C2.initialization(3,3,0.0);    
  S1.initialization(3,3,0.0);
  S2.initialization(3,3,0.0); 
  
  Matrix<double> fu,fu1,fu2,Kut;
  Matrix<double> fp,fp1,fp2,Kpt;
  Matrix<double> ft,ft1,ft2,Ktt;
  Matrix<double> Kup,Ktp;
  
  int ndofe = nne*ndofn;

   fu.initialization(nne*nsd,1,    0.0);
  fu1.initialization(nne*nsd,1,    0.0);
  fu2.initialization(nne*nsd,1,    0.0);  
   fp.initialization(npres,  1,    0.0);
  fp1.initialization(npres,  1,    0.0);
  fp2.initialization(npres,  1,    0.0);  
   ft.initialization(nVol,   1,    0.0);
  ft1.initialization(nVol,   1,    0.0);
  ft2.initialization(nVol,   1,    0.0);  
  Kut.initialization(nne*nsd,nVol, 0.0);
  Kpt.initialization(npres,  nVol, 0.0);
  Ktt.initialization(nVol,   nVol, 0.0);    
  Kup.initialization(nne*nsd,npres,0.0);
  Ktp.initialization(nVol,   npres,0.0);
  
  FEMLIB fe(ii, elem, node, INTG_ORDER,1);
                        
  Matrix<double> Np, Nt;
  Np.initialization(npres,1,0.0);
  Nt.initialization(nVol, 1,0.0);    
                         
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    fe.elem_basis_V(ip);
    fe.elem_shape_function(ip,npres,Np.m_pdata);
    fe.elem_shape_function(ip,nVol, Nt.m_pdata);
    fe.update_shape_tensor();
    fe.update_deformation_gradient(ndofn,u1,F1.m_pdata);
    fe.update_deformation_gradient(ndofn,u2,F2.m_pdata);    
      
    double Tn1 = 0.0; // 1: n-1+alpha
    double Tn2 = 0.0; // 2: n+alpha
    double Pn1 = 0.0;
    double Pn2 = 0.0;

    for(int a=0; a<nVol; a++)
    {
      Tn1 += Nt(a+1)*(alpha_1*eps[ii].T[a*3+2] + alpha_2*eps[ii].T[a*3+1]);
      Tn2 += Nt(a+1)*(alpha_1*eps[ii].T[a*3+1] + alpha_2*eps[ii].T[a*3+0]);      
    }
    for(int a=0; a<npres; a++)
    {
      Pn1 += Np(a+1)*P1[a];
      Pn2 += Np(a+1)*P2[a];
    }
      
                
    double Upp2 = 0.0;
    double Up1 = 0.0;
    double Up2 = 0.0;    
    UPP(Tn2,&hommat[mat],&Upp2);
    Upp2 = Upp2*kappa;

    UP(Tn1, &hommat[mat], &Up1);
    Up1 = Up1*kappa;    
    
    UP(Tn2, &hommat[mat], &Up2);
    Up2 = Up2*kappa;    
    
    C1.prod(F1,1,F1,0);
    Stress(C1.m_pdata,&hommat[mat],S1.m_pdata);
    
    C2.prod(F2,1,F2,0);
    Stress(C2.m_pdata,&hommat[mat],S2.m_pdata);    
    		
    TF_Kpt_ip(Kpt,&fe,F2,npres,nVol,Np,Nt,dt_alpha_1_minus_alpha);
    TF_Ktt_ip(Ktt,&fe,F2,nVol,Nt,Upp2,dt_alpha_1_minus_alpha);        

    if(npres==nne)
      TF_Kut_ip(Kut,&fe,F2,nVol,Nt,dt_alpha_1_minus_alpha);
      
		if(npres==1)
      TF_Kup_ip(Kup,&fe,F2,npres,Np,dt_alpha_1_minus_alpha);     

    
    TF_Ru_ip(fu1,&fe,F1,S1,Pn1);
    TF_Rp_ip(fp1,&fe,F1,npres,Np,Tn1);
    TF_Rt_ip(ft1,&fe,F1,nVol,Nt,Pn1,Up1);
    
    TF_Ru_ip(fu2,&fe,F2,S2,Pn2);
    TF_Rp_ip(fp2,&fe,F2,npres,Np,Tn2);
    TF_Rt_ip(ft2,&fe,F2,nVol,Nt,Pn2,Up2);    
  }
  Ktp.trans(Kpt);

  for(int a=0; a<nne*nsd; a++)
  	fu.m_pdata[a] = dt_alpha_1*fu2.m_pdata[a] + dt_alpha_2*fu1.m_pdata[a];  	

  for(int a=0; a<nVol; a++)
  	ft.m_pdata[a] = dt_alpha_1*ft2.m_pdata[a] + dt_alpha_2*ft1.m_pdata[a];

  for(int a=0; a<npres; a++)
  	fp.m_pdata[a] = dt_alpha_1*fp2.m_pdata[a] + dt_alpha_2*fp1.m_pdata[a];


  condense_F_out(f,nne,nsd,npres,nVol,fu.m_pdata,ft.m_pdata,fp.m_pdata,
                    Kut.m_pdata,Kup.m_pdata,Ktp.m_pdata,Ktt.m_pdata,Kpt.m_pdata);

  free(P1); free(P2); free(u1); free(u2);
}

void update_3f_state_variables_ip(int ii, int ip, 
        const ELEMENT *elem,
        const HOMMAT *hommat,
        SIG *sig,
        EPS *eps,
        double *F,
        double pressure,
        double theta,
        double volume,
        double detJxW)
{
  if(ip==1)
  { 
    memset(sig[ii].el.o,0,6*sizeof(double));
    memset(eps[ii].el.o,0,6*sizeof(double));
  }    
  
  const int mat = elem[ii].mat[2];
    
  double F_total[9], C[9], C_I[9], S[9], AA[9];
  int err = 0;
  double J = getJacobian(F,ii,&err);
  double alpha = pow(theta/J,1.0/3.0);
  
  for(int a = 0; a<9; a++)
    F_total[a] = alpha*F[a];

  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
  		3,3,3,1.0,F_total,3,F_total,3,0.0,C,3);
  
  inverse(C,3,C_I);
  
  // Get Deviatoric 2 P-K stress
  devStressFuncPtr Stress;
  Stress = getDevStressFunc(1,&hommat[mat]);
  Stress(C,&hommat[mat],S);
    
  /* Compute total stress S = dev(S) + p*Tn*Tr*C_I */
  cblas_daxpy(9,pressure*theta,C_I,1,S,1);  
    
	/* store S at ip */
	sig[ii].il[ip-1].o[0] = S[idx_2(0,0)]; /* XX */
	sig[ii].il[ip-1].o[1] = S[idx_2(1,1)]; /* YY */
	sig[ii].il[ip-1].o[2] = S[idx_2(2,2)]; /* ZZ */
	sig[ii].il[ip-1].o[3] = S[idx_2(1,2)]; /* YZ */
	sig[ii].il[ip-1].o[4] = S[idx_2(0,2)]; /* XZ */
	sig[ii].il[ip-1].o[5] = S[idx_2(0,1)]; /* XY */	   
  
  /* Compute Cauchy stress (theta)^-1 F S F' */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
  		3,3,3,1.0,F_total,3,S,3,0.0,AA,3);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
  		3,3,3,1.0/theta,AA,3,F_total,3,0.0,S,3);  
    
  sig[ii].el.o[0] += detJxW/volume*S[idx_2(0,0)];
  sig[ii].el.o[1] += detJxW/volume*S[idx_2(1,1)];
  sig[ii].el.o[2] += detJxW/volume*S[idx_2(2,2)];

  sig[ii].el.o[3] += detJxW/volume*S[idx_2(1,2)];
  sig[ii].el.o[4] += detJxW/volume*S[idx_2(0,2)];
  sig[ii].el.o[5] += detJxW/volume*S[idx_2(0,1)];	
  
  double E[9], F_total_I[9];
  
  /* Compute G-L strain E = 0.5(C - 1) */
  for(int a=0; a<9; a++)
    E[a] = 0.5*(C[a] - 1.0*(a==0 || a==4 || a==8));

	/* Elastic Green Lagrange strain */
	eps[ii].il[ip-1].o[0] = E[idx_2(0,0)];
	eps[ii].il[ip-1].o[1] = E[idx_2(1,1)];
	eps[ii].il[ip-1].o[2] = E[idx_2(2,2)];
	eps[ii].il[ip-1].o[3] = E[idx_2(1,2)]*2.0;
	eps[ii].il[ip-1].o[4] = E[idx_2(0,2)]*2.0;
	eps[ii].il[ip-1].o[5] = E[idx_2(0,1)]*2.0;

  
  /* Compute logarithmic strain */
  inverse(F_total,3,F_total_I);
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
  		3,3,3,1.0,F_total_I,3,E,3,0.0,AA,3);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
  		3,3,3,1.0,AA,3,F_total_I,3,0.0,E,3);  
  
  eps[ii].el.o[0] += detJxW/volume*E[idx_2(0,0)];
  eps[ii].el.o[1] += detJxW/volume*E[idx_2(1,1)];
  eps[ii].el.o[2] += detJxW/volume*E[idx_2(2,2)];  
  eps[ii].el.o[3] += detJxW/volume*E[idx_2(1,2)]*2.0;
  eps[ii].el.o[4] += detJxW/volume*E[idx_2(0,2)]*2.0;
  eps[ii].el.o[5] += detJxW/volume*E[idx_2(0,1)]*2.0;    
}  


void update_3f_state_variables_el(const int ii,
        const int ndofn,
        const int nne,
        const int npres,
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const NODE *node,
        double *u,
        const double *P,
        double dt,
        SIG *sig,
        EPS *eps)
{

  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));   
  const double volume = Tetra_V(x,y,z);       			    
  
  const int nsd = 3;
  const int nVol = N_VOL_TREE_FIELD;
  
  Matrix<double> F(3,3,0.0);
  
  FEMLIB fe(ii, elem, node, 1,1);
    
  Matrix<double> Np(npres,1,0.0), Nt(nVol, 1,0.0);       
    
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    fe.elem_basis_V(ip);  
    fe.elem_shape_function(ip,npres,Np.m_pdata);
    fe.elem_shape_function(ip,nVol, Nt.m_pdata);
    fe.update_shape_tensor();
    fe.update_deformation_gradient(ndofn,u,F.m_pdata);
                      
    double Tn = 0.0; 
    double Pn = 0.0;
    
    for(int a=0; a<nVol; a++)
      Tn += Nt(a+1)*eps[ii].T[a*3+0];
    
    for(int a=0; a<npres; a++)
      Pn += Np(a+1)*P[a];
                      
    update_3f_state_variables_ip(ii,ip,elem,hommat,sig,eps,F.m_pdata,Pn,Tn,volume,fe.detJxW);

  }
}


void evaluate_PT_el(const int ii,
        const int ndofn,
        const int nne,
        const int npres,
        const int nVol,
        const int nsd,        
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const NODE *node,
        const double *r_e,
        const double *du,
        double dt,
        SIG *sig,
        EPS *eps)
{
  double *P, *u;
  u = aloc1(nne*nsd);
	P = aloc1(npres);
		
  for(int a=0;a<nne;a++)
  {
  	for(int b=0; b<nsd;b++)
  		u[a*nsd+b] = r_e[a*ndofn+b];
 
  	if(npres==nne)  		
  	  P[a] = r_e[a*ndofn+nsd];
  }
  if(npres==1)
    P[0] = eps[ii].d_T[0];
  
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));      			    
  
  dUdJFuncPtr UP = getDUdJFunc(1, &hommat[mat]);
  d2UdJ2FuncPtr UPP = getD2UdJ2Func(1, &hommat[mat]);
  devStressFuncPtr Stress = getDevStressFunc(1,&hommat[mat]);
  
  int ndofe = nne*ndofn;
  
  Matrix<double> F(3,3,0.0);
	
  Matrix<double> fp(npres, 1, 0.0), ft(nVol,  1, 0.0);
		
  Matrix<double> Ktu,Ktt,Ktp;
  Matrix<double> Kup,Kpu,Kpt;
  
  Ktu.initialization(nVol   ,nne*nsd,0.0);
  Ktt.initialization(nVol   ,nVol   ,0.0);
  Kup.initialization(nne*nsd,npres  ,0.0);      
  Kpu.initialization(npres  ,nne*nsd,0.0);
  Kpt.initialization(npres  ,nVol   ,0.0);
  Ktp.initialization(nVol   ,npres  ,0.0);
  
  FEMLIB fe(ii, elem, node, INTG_ORDER,1);
  
  Matrix<double> Np(npres,1,0.0), Nt(nVol, 1,0.0);    
                         
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    fe.elem_basis_V(ip);  
    fe.elem_shape_function(ip,npres,Np.m_pdata);
    fe.elem_shape_function(ip,nVol, Nt.m_pdata);
    fe.update_shape_tensor();
    fe.update_deformation_gradient(ndofn,u,F.m_pdata);
              
    double Tn = 0.0; 
    double Pn = 0.0;
    
    for(int a=0; a<nVol; a++)
      Tn += Nt(a+1)*eps[ii].T[a*3+0];
    
    for(int a=0; a<npres; a++)
      Pn += Np(a+1)*P[a];
      
                
    double Upp = 0.0;
    double Up = 0.0;
    UPP(Tn,&hommat[mat],&Upp);
    Upp = Upp*kappa;

    UP(Tn, &hommat[mat], &Up);
    Up = Up*kappa;
                  
    TF_Kup_ip(Kup,&fe,F,npres,Np,-1.0);     
    
    TF_Kpt_ip(Kpt,&fe,F,npres,nVol,Np,Nt,-1.0);
    TF_Ktt_ip(Ktt,&fe,F,nVol,Nt,Upp,-1.0);        

    TF_Rt_ip(ft,&fe,F,nVol,Nt,Pn,Up);
    TF_Rp_ip(fp,&fe,F,npres,Np,Tn);				
  }

  free(u);
  free(P);
  
  Kpu.trans(Kup);
  Ktp.trans(Kpt);
  
  Matrix<double> theta(nVol,1,0.0),KptI(nVol,nVol,0.0);  
  KptI.inv(Kpt);
  
  Matrix<double> uu;
  uu.use_reference(nne*ndofn, 1, du);
    
  Matrix<double> Kpu_uu;
  Kpu_uu.prod(Kpu,uu);
  fp.add(Kpu_uu);
  
  theta.prod(KptI,fp);
  theta.prod(-1.0);
  
	for(int a = 0; a<nVol; a++)
	  eps[ii].T[a*3+0] += theta.m_pdata[a];

  Matrix<double> press(npres,1,0.0),KtpI(nVol,npres,0.0);  
  KtpI.inv(Ktp);
  
  Matrix<double> Ktt_theta;
  Ktt_theta.prod(Ktt,theta);
  ft.add(Ktt_theta);
  
  press.prod(KtpI,ft);
  press.prod(-1.0);  
  
  for(int a = 0; a<npres; a++)
    eps[ii].d_T[a*3+0] += press.m_pdata[a]; 
}

void evaluate_PT_w_inertia_el(const int ii,
        const int ndofn,
        const int nne,
        const int npres,
        const int nVol,
        const int nsd,
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const NODE *node,
        const double *du,
        double dt,
        SIG *sig,
        EPS *eps,
        double alpha, double *r_n_a, double *r_n_1_a)
{ 
  double alpha_1;
  double alpha_2;
  double dt_alpha_1_minus_alpha;
  double dt_alpha_1;
  double dt_alpha_2;

  if(alpha<0)
  {  
    alpha = 1.0;
    alpha_1 = 0.0;
    alpha_2 = alpha;
    dt_alpha_1_minus_alpha = -1.0; 
    dt_alpha_1 = 1.0;
    dt_alpha_2 = 0.0;
  }
  else
  {
    alpha_1 = 1.0 - alpha;
    alpha_2 = alpha;
    dt_alpha_1_minus_alpha = dt*alpha_1*alpha_2; 
    dt_alpha_1 = -dt*alpha_1;  
    dt_alpha_2 = -dt*alpha_2;  
  } 
  
  double *u2 = aloc1(nne*nsd); // 2: n+alpha
  double *u1 = aloc1(nne*nsd); // 1: n-1+alpha
  double *P2 = aloc1(npres);   // 1: n+alpha
  double *P1 = aloc1(npres);   // 1: n-1+alpha  
  
  for(int a=0;a<nne;a++)
  {
  	for(int b=0; b<nsd;b++)
  	{
  		u1[a*nsd+b] = r_n_1_a[a*ndofn+b]; 		  			
  		u2[a*nsd+b] = r_n_a[a*ndofn+b];
  	}
  	if(npres==nne)
  	{
  	  P1[a] = r_n_1_a[a*ndofn+nsd];
  	  P2[a] = r_n_a[a*ndofn+nsd];
  	}
  }	
  if(npres==1)
  {  
    P1[0] = alpha_1*eps[ii].d_T[2] + alpha_2*eps[ii].d_T[1]; 
    P2[0] = alpha_1*eps[ii].d_T[1] + alpha_2*eps[ii].d_T[0];
  }  
  
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.0*(1.0-2.0*hommat[mat].nu));
  
  dUdJFuncPtr          UP = getDUdJFunc(1, &hommat[mat]);
  d2UdJ2FuncPtr       UPP = getD2UdJ2Func(1, &hommat[mat]);
  devStressFuncPtr Stress = getDevStressFunc(1,&hommat[mat]);
             
  int ndofe = nne*ndofn;


  Matrix<double> F1(3,3,0.0),F2(3,3,0.0);      
	
  Matrix<double> fp, fp1, fp2;
  Matrix<double> ft, ft1, ft2;
  
   fp.initialization(npres,  1,    0.0);
  fp1.initialization(npres,  1,    0.0);
  fp2.initialization(npres,  1,    0.0);
   ft.initialization(nVol,   1,    0.0);  
  ft1.initialization(nVol,   1,    0.0);  
  ft2.initialization(nVol,   1,    0.0);
		
  Matrix<double> Ktu,Ktt,Ktp;
  Matrix<double> Kup,Kpu,Kpt;
  
  Ktu.initialization(nVol   ,nne*nsd,0.0);
  Ktt.initialization(nVol   ,nVol   ,0.0);
  Kup.initialization(nne*nsd,npres  ,0.0);      
  Kpu.initialization(npres  ,nne*nsd,0.0);
  Kpt.initialization(npres  ,nVol   ,0.0);
  Ktp.initialization(nVol   ,npres  ,0.0);   

  FEMLIB fe(ii, elem, node, INTG_ORDER,1);
    
  Matrix<double> Np(npres,1,0.0),Nt(nVol, 1,0.0);    
                             
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    fe.elem_basis_V(ip);  
    fe.elem_shape_function(ip,npres,Np.m_pdata);
    fe.elem_shape_function(ip,nVol, Nt.m_pdata);   
    fe.update_shape_tensor();
    fe.update_deformation_gradient(ndofn,u1,F1.m_pdata);
    fe.update_deformation_gradient(ndofn,u2,F2.m_pdata);    
    
    double Tn1 = 0.0; // 1: n-1+alpha
    double Tn2 = 0.0; // 2: n+alpha
    double Pn1 = 0.0;
    double Pn2 = 0.0;
    
    for(int a=0; a<nVol; a++)
    {
      Tn1 += Nt(a+1)*(alpha_1*eps[ii].T[a*3+2] + alpha_2*eps[ii].T[a*3+1]);
      Tn2 += Nt(a+1)*(alpha_1*eps[ii].T[a*3+1] + alpha_2*eps[ii].T[a*3+0]);      
    }
    for(int a=0; a<npres; a++)
    {
      Pn1 += Np(a+1)*P1[a];
      Pn2 += Np(a+1)*P2[a];
    }
        
    double Upp2 = 0.0;
    double Up1 = 0.0;
    double Up2 = 0.0;    
    UPP(Tn2,&hommat[mat],&Upp2);
    Upp2 = Upp2*kappa;

    UP(Tn1, &hommat[mat], &Up1);
    Up1 = Up1*kappa;    
    
    UP(Tn2, &hommat[mat], &Up2);
    Up2 = Up2*kappa;        
    
    TF_Kup_ip(Kup,&fe,F2,npres,Np,dt_alpha_1_minus_alpha);
    
    TF_Kpt_ip(Kpt,&fe,F2,npres,nVol,Np,Nt,dt_alpha_1_minus_alpha);
    TF_Ktt_ip(Ktt,&fe,F2,nVol,Nt,Upp2,dt_alpha_1_minus_alpha); 
    
    TF_Rt_ip(ft1,&fe,F1,nVol,Nt,Pn1,Up1);
    TF_Rt_ip(ft2,&fe,F2,nVol,Nt,Pn2,Up2);
        
    TF_Rp_ip(fp1,&fe,F1,npres,Np,Tn1);	
    TF_Rp_ip(fp2,&fe,F2,npres,Np,Tn2);	
  }
  
  Kpu.trans(Kup);
  Ktp.trans(Kpt); 

  for(int a=0; a<nVol; a++)
  	ft.m_pdata[a] = dt_alpha_1*ft2.m_pdata[a] + dt_alpha_2*ft1.m_pdata[a];

  for(int a=0; a<npres; a++)
  	fp.m_pdata[a] = dt_alpha_1*fp2.m_pdata[a] + dt_alpha_2*fp1.m_pdata[a];
  	    
  free(P1); free(P2); free(u1); free(u2);
  
  Matrix<double> theta(nVol,1,0.0),KptI(nVol,npres,0.0);  
  KptI.inv(Kpt);
  
  
  Matrix<double> uu;
  uu.use_reference(nne*ndofn, 1, du);
    
  Matrix<double> Kpu_uu;
  Kpu_uu.prod(Kpu,uu);
  fp.add(Kpu_uu);
  
  theta.prod(KptI,fp);
  theta.prod(-1.0);
  
	for(int a = 0; a<nVol; a++)
	  eps[ii].T[a*3+0] += theta.m_pdata[a];

  Matrix<double> press(npres,1,0.0),KtpI(nVol,npres,0.0);  
  KtpI.inv(Ktp);
  
  Matrix<double> Ktt_theta;
  Ktt_theta.prod(Ktt,theta);
  ft.add(Ktt_theta);
  
  press.prod(KtpI,ft);
  press.prod(-1.0);  
  
  for(int a = 0; a<npres; a++)
    eps[ii].d_T[a*3+0] += press.m_pdata[a];   
}

void evaluate_theta_el(const int ii,
        const int ndofn,
        const int nne,
        const int npres,
        const int nVol,
        const int nsd,
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const NODE *node,
        const double *r_e,
        const double *du,
        double dt,
        SIG *sig,
        EPS *eps)
{
  double *P, *u;
  u = aloc1(nne*nsd);
	P = aloc1(npres);
		
  for(int a=0;a<nne;a++)
  {
  	for(int b=0; b<nsd;b++)
  		u[a*nsd+b] = r_e[a*ndofn+b];
 
  	if(npres==nne)  		
  			P[a] = r_e[a*ndofn+nsd];
  }
  if(npres==1)
    P[0] = eps[ii].d_T[0];  
  
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  
  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UP = getDUdJFunc(1, &hommat[mat]);  
  UPP = getD2UdJ2Func(1, &hommat[mat]);
  
  const double Vol = Tetra_V (x,y,z);
  int err = 0;  
  int ndofe = nne*ndofn;

  double *ft  = aloc1(nVol);   
  double *Ktu = aloc1(nne*nsd*nVol);
  double *Ktp = aloc1(npres*nVol);
  double *Ktt = aloc1(nVol*nVol);

  memset(ft, 0,        nVol*sizeof(double));      
	memset(Ktu,0,nne*nsd*nVol*sizeof(double));
	memset(Ktp,0,  nVol*npres*sizeof(double));    
	memset(Ktt,0,     nVol*nVol*sizeof(double));	
  
  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int itg_order = nne;
  int_point(nne,&npt_z);
  
  double *int_pt_ksi, *int_pt_eta, *int_pt_zet, *weights;
  int_pt_ksi = aloc1(npt_z);
  int_pt_eta = aloc1(npt_z);
  int_pt_zet = aloc1(npt_z);
  weights = aloc1(npt_z);
  
  /* allocate space for the shape functions, derivatives etc */
  double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  double Upp, Up;
 
	Na = aloc1(nne); 
  Np = aloc1(npres);
  Nt = aloc1(nVol);
  N_x = aloc1(nne);
  N_y = aloc1(nne);
  N_z = aloc1(nne);
  ST_tensor = aloc4(3,3,nsd,nne);
  ST = aloc1(3*3*nsd*nne);
  
  /*=== INTEGRATION LOOP ===*/
  integrate(nne,&npt_x,&npt_y,&npt_z,
          int_pt_ksi,int_pt_eta,int_pt_zet,
          weights);
  
  double **F_mat,*F;
  
  F = aloc1(9);
  F_mat = aloc2(3,3);

  for(long i=0; i<npt_x; i++)
  {
    for(long j=0; j<npt_y; j++)
    {
      for(long k=0; k<npt_z; k++)
      {     	
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nne,Na);
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],npres,Np);
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nVol,Nt);
        
        double detJ = deriv(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nne,x,y,z,N_x,N_y,N_z);
        double wt = weights[k];
        
        double Tn = 0.0;
        double Pn = 0.0;
       
        for(int a=0; a<nVol; a++)
          Tn += Nt[a]*eps[ii].T[a*3+0];          
        
        for(int a=0; a<npres; a++)
          Pn += Np[a]*P[a];          
        
    
        UPP(Tn,&hommat[mat],&Upp);
        Upp = Upp*kappa;

        UP(Tn, &hommat[mat], &Up);
        Up = Up*kappa;
                
        shape_tensor(nne,ndn,N_x,N_y,N_z,ST_tensor);
        shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
        
        def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u,F_mat);
        mat2array(F,CONST_2(double) F_mat,3,3);
        
//        add_3F_Ktu_ip(Ktu,nne,nVol,  ST,F,detJ,wt,Nt,-1.0);        
				add_3F_Ktp_ip(Ktp,nVol,npres,detJ,wt,Nt,Np,-1.0);
				add_3F_Ktt_ip(Ktt,nVol,detJ,wt,Nt,Upp,-1.0);	

				resid_w_inertia_Rt_ip(ft, nVol, detJ, wt, Nt, Pn, Up);
								
      }
    }
  }

  free(u);
  free(P);
     	  	  	 
  dealoc4(ST_tensor,3,3,nsd);
  free(ST);
  
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(Np);
  free(Nt);
  free(N_x);
  free(N_y);
  free(N_z);
  
  free(F);
  dealoc2(F_mat,3);  
  
  	
	double *KttI, *KttIKtu, *KttIKtp;
  KttI = aloc1(nVol*nVol);
  KttIKtu = aloc1(nne*nsd*nVol);
  KttIKtp = aloc1(npres*nVol);
				
	inverse(Ktt,nVol,KttI);
     
  double *theta1 = aloc1(nVol);
  double *theta2 = aloc1(nVol);
  double *theta3 = aloc1(nVol);    
    
     
	// KttI:    [nVol]*[nVol]
	// Ktu:     [nVol]*[nne*nsd]
	// KttIKtu: [nVol]*[nne*nsd]
	// Ktp:     [nVol]*[npres]
	// KttIKpt: [nVol]*[npres]
        
  // A[m,k] x B[k,n] = C[m,n]
  // cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,A,k,B,n,0.0,C,n);      
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,nne*nsd,nVol,1.0,KttI,nVol,Ktu,nne*nsd,0.0,KttIKtu,nne*nsd);
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,npres,nVol,1.0,KttI,nVol,Ktp,npres,0.0,KttIKtp,npres); 
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nne*nsd,1.0,KttIKtu,nne*nsd,u,1,0.0,theta1,1);
                            
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,npres,1.0,KttIKtp,npres,P,1,0.0,theta2,1);
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nVol,1.0,KttI,nVol,ft,1,0.0,theta3,1);              
                                                         
	for(long a = 0; a<nVol; a++)
	  eps[ii].T[a*3+0] = eps[ii].T[a*3+0] - (theta1[a] + theta2[a] + theta3[a]);
  
  free(ft);
    
  free(Ktu); 
  free(Ktp); 
  free(Ktt); free(KttI);
  free(KttIKtu);
  free(KttIKtp);
  free(theta1); free(theta2); free(theta3);
}

void evaluate_theta_w_inertia_el(const int ii,
        const int ndofn,
        const int nne,
        const int npres,
        const int nVol,
        const int nsd,
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const NODE *node,
        const double *du,
        double dt,
        SIG *sig,
        EPS *eps,
        double alpha, double *r_n_a, double *r_n_1_a)
{
  double *u2 = aloc1(nne*nsd); // 2: n+alpha
  double *u1 = aloc1(nne*nsd); // 1: n-1+alpha
  double *P2 = aloc1(npres);   // 1: n+alpha
  double *P1 = aloc1(npres);   // 1: n-1+alpha  
  
  for(int a=0;a<nne;a++)
  {
  	for(int b=0; b<nsd;b++)
  	{
  		u1[a*nsd+b] = r_n_1_a[a*ndofn+b]; 		  			
  		u2[a*nsd+b] = r_n_a[a*ndofn+b];
  	}
  	if(npres==nne)
  	{
  	  P1[a] = r_n_1_a[a*ndofn+nsd];
  	  P2[a] = r_n_a[a*ndofn+nsd];
  	}
  }	
  if(npres==1)
  {
    P1[0] = (1.0-alpha)*eps[ii].d_T[2] + alpha*eps[ii].d_T[1]; 
    P2[0] = (1.0-alpha)*eps[ii].d_T[1] + alpha*eps[ii].d_T[0];             					  
  } 
    
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  
  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UP = getDUdJFunc(1, &hommat[mat]);  
  UPP = getD2UdJ2Func(1, &hommat[mat]);
  
  const double Vol = Tetra_V (x,y,z);
  int err = 0;  
  int ndofe = nne*ndofn;

  double *ft  = aloc1(nVol);   
  double *ft1 = aloc1(nVol);
  double *ft2 = aloc1(nVol);  
  double *Ktu = aloc1(nne*nsd*nVol);
  double *Ktp = aloc1(npres*nVol);
  double *Ktt = aloc1(nVol*nVol);

  memset(ft, 0,        nVol*sizeof(double));      
  memset(ft1,0,        nVol*sizeof(double));    
  memset(ft2,0,        nVol*sizeof(double));      
	memset(Ktu,0,nne*nsd*nVol*sizeof(double));
	memset(Ktp,0,  nVol*npres*sizeof(double));    
	memset(Ktt,0,     nVol*nVol*sizeof(double));	
  
  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(nne,&npt_z);
  
  double *int_pt_ksi, *int_pt_eta, *int_pt_zet, *weights;
  int_pt_ksi = aloc1(npt_z);
  int_pt_eta = aloc1(npt_z);
  int_pt_zet = aloc1(npt_z);
  weights = aloc1(npt_z);
  
  /* allocate space for the shape functions, derivatives etc */
  double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  double Upp2, Up1, Up2;
 
	Na = aloc1(nne); 
  Np = aloc1(npres);
  Nt = aloc1(nVol);
  N_x = aloc1(nne);
  N_y = aloc1(nne);
  N_z = aloc1(nne);
  ST_tensor = aloc4(3,3,nsd,nne);
  ST = aloc1(3*3*nsd*nne);
  
  /*=== INTEGRATION LOOP ===*/
  integrate(nne,&npt_x,&npt_y,&npt_z,
          int_pt_ksi,int_pt_eta,int_pt_zet,
          weights);
  
  double **F_mat1, **F_mat2, *F1, *F2;
  
  F1 = aloc1(9);
  F2 = aloc1(9);
  F_mat1 = aloc2(3,3);
  F_mat2 = aloc2(3,3);

  for(long i=0; i<npt_x; i++)
  {
    for(long j=0; j<npt_y; j++)
    {
      for(long k=0; k<npt_z; k++)
      {     	
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nne,Na);
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],npres,Np);
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nVol,Nt);
        
        double detJ = deriv(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nne,x,y,z,N_x,N_y,N_z);
        double wt = weights[k];
        
        double Tn1 = 0.0; // 1: n-1+alpha
        double Tn2 = 0.0; // 2: n+alpha
        double Pn1 = 0.0;
        double Pn2 = 0.0;
       
        for(int a=0; a<nVol; a++)
        {
          Tn1 += Nt[a]*((alpha-1.0)*eps[ii].T[a*3+2] + alpha*eps[ii].T[a*3+1]);
          Tn2 += Nt[a]*((alpha-1.0)*eps[ii].T[a*3+1] + alpha*eps[ii].T[a*3+0]);          
        }
        
        for(int a=0; a<npres; a++)
        {
          Pn1 += Np[a]*P1[a];
          Pn2 += Np[a]*P2[a];          
        }
    
        UPP(Tn2,&hommat[mat],&Upp2);
        Upp2 = Upp2*kappa;

        UP(Tn1, &hommat[mat], &Up1);
        UP(Tn2, &hommat[mat], &Up2);
        Up1 = Up1*kappa;
        Up2 = Up2*kappa;        
                
        shape_tensor(nne,ndn,N_x,N_y,N_z,ST_tensor);
        shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
        
        def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u1,F_mat1);
        mat2array(F1,CONST_2(double) F_mat1,3,3);
        
        def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u2,F_mat2);
        mat2array(F2,CONST_2(double) F_mat2,3,3);        

//        add_3F_Ktu_ip(Ktu,nne,nVol,  ST,F2,detJ,wt,Nt,dt*(1.0-alpha)*alpha);        
				add_3F_Ktp_ip(Ktp,nVol,npres,detJ,wt,Nt,Np,dt*(1.0-alpha)*alpha);
				add_3F_Ktt_ip(Ktt,nVol,detJ,wt,Nt,Upp2,dt*(1.0-alpha)*alpha);	

				resid_w_inertia_Rt_ip(ft1, nVol, detJ, wt, Nt, Pn1, Up1);
				resid_w_inertia_Rt_ip(ft2, nVol, detJ, wt, Nt, Pn2, Up2);
								
      }
    }
  }
  
  for(int a=0; a<nVol; a++)
  	ft[a] = -(1.0 - alpha)*dt*ft2[a] - alpha*dt*ft1[a];
  	  	  	 
  free(P1); free(P2); free(u1); free(u2);
  dealoc4(ST_tensor,3,3,nsd);
  free(ST);
  
  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(Np);
  free(Nt);
  free(N_x);
  free(N_y);
  free(N_z);
  
  free(F1);
  dealoc2(F_mat1,3);  
  free(F2);
  dealoc2(F_mat2,3);   
  
  	
	double *KttI, *KttIKtu, *KttIKtp;
  KttI = aloc1(nVol*nVol);
  KttIKtu = aloc1(nne*nsd*nVol);
  KttIKtp = aloc1(npres*nVol);
				
	inverse(Ktt,nVol,KttI);
     
  double *theta1 = aloc1(nVol);
  double *theta2 = aloc1(nVol);
  double *theta3 = aloc1(nVol);    
    
     
	// KttI:    [nVol]*[nVol]
	// Ktu:     [nVol]*[nne*nsd]
	// KttIKtu: [nVol]*[nne*nsd]
	// Ktp:     [nVol]*[npres]
	// KttIKpt: [nVol]*[npres]
        
  // A[m,k] x B[k,n] = C[m,n]
  // cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,A,k,B,n,0.0,C,n);      
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,nne*nsd,nVol,1.0,KttI,nVol,Ktu,nne*nsd,0.0,KttIKtu,nne*nsd);
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,npres,nVol,1.0,KttI,nVol,Ktp,npres,0.0,KttIKtp,npres); 
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nne*nsd,1.0,KttIKtu,nne*nsd,u2,1,0.0,theta1,1);
                            
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,npres,1.0,KttIKtp,npres,P2,1,0.0,theta2,1);
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nVol,1.0,KttI,nVol,ft,1,0.0,theta3,1);              
                                                         
	for(long a = 0; a<nVol; a++)
	  eps[ii].T[a*3+0] = eps[ii].T[a*3+0] - (theta1[a] + theta2[a] + theta3[a]);
  
  free(ft);
  free(ft1); free(ft2);
    
  free(Ktu); 
  free(Ktp); 
  free(Ktt); free(KttI);
  free(KttIKtu);
  free(KttIKtp);
  free(theta1); free(theta2); free(theta3);
}

void update_3f(long ne, long ndofn, long npres, double *d_r, double *r, double *rr,
               NODE *node, ELEMENT *elem, HOMMAT *hommat, SUPP sup, EPS *eps, SIG *sig, double dt, double t,
		           MPI_Comm mpi_comm, const PGFem3D_opt *opts, double alpha, double *r_n, double *r_n_1,
		           const int mp_id)
{
  const int mat = elem[0].mat[2];
  double rho = hommat[mat].density;
  long include_inertia = 1;

  const int nsd = 3;
  const int nVol = N_VOL_TREE_FIELD;
  
  if(fabs(rho)<MIN_DENSITY)
    include_inertia = 0;

//////////////////////////////////////////////////////////////////////    
//  include_inertia = 0;    
//////////////////////////////////////////////////////////////////////

            
  double err = 0.0;
      
  for (int i=0;i<ne;i++)
  {

    int nne = elem[i].toe;
    long *nod = aloc1l (nne);
    elemnodes (i,nne,nod,elem);
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);

    double *r_e = aloc1 (ndofe);
    double *dr_e = aloc1 (ndofe);

    double *x,*y,*z;
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);

    long *cn = aloc1l (ndofe);

    nodecoord_total(nne,nod,node,x,y,z);

    /* code numbers on element */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);
    
    /* deformation on element */
    def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
    def_elem(cn,ndofe,rr,elem,node,dr_e,sup,2);
    
    double*du;
    du = aloc1(nne*nsd);
	  
    for(int a=0;a<nne;a++)
    {
    	for(int b=0; b<nsd;b++)
    		du[a*nsd+b] = dr_e[a*ndofn+b];
    }	  
	  
	  if(include_inertia && (0<alpha && alpha<1))
	  {
	    double *r0, *r0_;
  	  r0      = aloc1(ndofe);
	    r0_     = aloc1(ndofe);
	
  	  for (long I=0;I<nne;I++)
	    {
	      for(long J=0; J<nsd; J++)
	      {
	        r0[I*ndofn + J] =   r_n[nod[I]*ndofn + J];
	        r0_[I*ndofn + J] = r_n_1[nod[I]*ndofn + J];            
	      }
	    }	  
	  
   		double *r_n_a, *r_n_1_a;
   		r_n_a = aloc1(ndofe);
   		r_n_1_a = aloc1(ndofe);
 			
  	  mid_point_rule(r_n_1_a,r0_,r0,  alpha, ndofe); 
  	  mid_point_rule(r_n_a,  r0, r_e, alpha, ndofe); 		
   			   	  		
  		if(npres==1)
  		  evaluate_PT_w_inertia_el(i,ndofn,nne,npres,nVol,nsd,x,y,z,elem,hommat,node,du,dt,sig,eps,alpha,r_n_a,r_n_1_a);
  		else
  		  evaluate_theta_w_inertia_el(i,ndofn,nne,npres,nVol,nsd,x,y,z,elem,hommat,node,du,dt,sig,eps,alpha,r_n_a,r_n_1_a); 

	    free(r0);
	    free(r0_);
		  free(r_n_a);
		  free(r_n_1_a);
    }
    else
    {
  		if(npres==1)
  		  evaluate_PT_el(i,ndofn,nne,npres,nVol,nsd,x,y,z,elem,hommat,node,r_e,du,dt,sig,eps);
  		else
  		  evaluate_theta_el(i,ndofn,nne,npres,nVol,nsd,x,y,z,elem,hommat,node,r_e,du,dt,sig,eps); 
  	}

 	  free(du);  		  	
   	free(nod);
   	free(cn);
    free(x);
    free(y);
    free(z);  	
  	free(r_e);
  	free(dr_e);
	}
}

void update_3f_state_variables(long ne, long ndofn, long npres, double *d_r, double *r,
               NODE *node, ELEMENT *elem, HOMMAT *hommat, SUPP sup, EPS *eps, SIG *sig, double dt, double t,
		           MPI_Comm mpi_comm,const int mp_id)
{
  const int mat = elem[0].mat[2];
  double rho = hommat[mat].density;
  long include_inertia = 1;
  
  if(fabs(rho)<MIN_DENSITY)
    include_inertia = 0;
            
  double err = 0.0;
  
  for (int i=0;i<ne;i++)
  {

    int nne = elem[i].toe;
    long *nod = aloc1l (nne);
    elemnodes (i,nne,nod,elem);
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);

    double *r_e = aloc1 (ndofe);

    double *x,*y,*z;
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);

    long *cn = aloc1l (ndofe);

    nodecoord_total(nne,nod,node,x,y,z);

    // code numbers on element
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);
    
    // deformation on element
    def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
    

  	int nsd = 3;
	  
		double *u, *P;
		u = aloc1(nne*nsd);
		P = aloc1(npres);

		for(int a=0;a<nne;a++)
		{
			for(int b=0; b<nsd;b++)
			{
				u[a*nsd+b] = r_e[a*ndofn+b];
			}
			if(npres==nne)  		
  			P[a] = r_e[a*ndofn+nsd];
		}
    if(npres==1)
      P[0] = eps[i].d_T[0];
		
    update_3f_state_variables_el(i,ndofn,nne,npres,x,y,z,elem,hommat,node,u,P,dt,sig,eps);

	  free(u);
		free(P);
  		  	
   	free(nod);
   	free(cn);
    free(x);
    free(y);
    free(z);  	
  	free(r_e);
	}
}