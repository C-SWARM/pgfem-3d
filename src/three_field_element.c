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


#define ndn 3
#define N_VOL_TF 1

void add_3F_Kuu_ip_disp(double *K, FEMLIB *fe, 
        double *ST, Matrix(double) F, Matrix(double) S, double *L,
        double dt_alpha_1_minus_alpha)
{
  int nne = fe->nne;
  int nsd = fe->nsd;
  Matrix(double) ST_ab, ST_wg, AA, BB, CC;
  Matrix(double) sAA, sBB, sCC;  
  Matrix_construct(double,ST_ab);
  Matrix_construct(double,ST_wg);
  Matrix_construct_redim(double,AA,3,3);
  Matrix_construct_redim(double,BB,3,3);
  Matrix_construct_redim(double,CC,3,3);
  
  Matrix_construct_redim(double,sAA,3,3);
  Matrix_construct_redim(double,sBB,3,3);
  Matrix_construct_redim(double,sCC,3,3);
  
  
  double LsBB[9];
    
  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const double* const ptrST_ab = &ST[idx_4_gen(a,b,0,0,
              nne,nsd,nsd,nsd)];
      Matrix_init_w_array(ST_ab,3,3,ptrST_ab);

      // AA = F^T Grad(del u)      
      Matrix_AxB(AA,1.0,0.0,F,1,ST_ab,0); 
      symmetric_part(sAA.m_pdata,AA.m_pdata,3);
                          
      for(int w=0; w<nne; w++)
      {
        for(int g=0; g<nsd; g++)
        {
          const double * const ptrST_wg = &ST[idx_4_gen(w,g,0,0,
                  nne,nsd,nsd,nsd)];
          Matrix_init_w_array(ST_wg,3,3,ptrST_wg);        
                                    
          // BB = F^T Grad(du) 
          Matrix_AxB(BB,1.0,0.0,F,1,ST_wg,0);
          symmetric_part(sBB.m_pdata,BB.m_pdata,3);
          
          // CC = Grad(del u)^T Grad(du) 
          Matrix_AxB(CC,1.0,0.0,ST_ab,1,ST_wg,0);
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
  Matrix_cleanup(ST_ab);
  Matrix_cleanup(ST_wg);
  Matrix_cleanup(AA);
  Matrix_cleanup(BB);
  Matrix_cleanup(CC);
  Matrix_cleanup(sAA);
  Matrix_cleanup(sBB);
  Matrix_cleanup(sCC);    
}

void add_3F_Kuu_ip_(double *K, FEMLIB *fe, 
        double *ST, Matrix(double) F, double Pn, double Tn,
        double dt_alpha_1_minus_alpha)
{
  int nne = fe->nne;
  int nsd = fe->nsd;
  Matrix(double) F_I, ST_ab, ST_wg, AA, BB, CC;
  Matrix_construct(double,F_I);  Matrix_redim(F_I,3,3);
  Matrix_construct(double,ST_ab);
  Matrix_construct(double,ST_wg);
  Matrix_construct(double,AA);  Matrix_redim(AA,3,3);
  Matrix_construct(double,BB);  Matrix_redim(BB,3,3);
  Matrix_construct(double,CC);  Matrix_redim(CC,3,3);      
  
  Matrix_inv(F, F_I);
  double Jn;
  Matrix_det(F, Jn);
    
  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const double* const ptrST_ab = &ST[idx_4_gen(a,b,0,0,
              nne,nsd,nsd,nsd)];
      Matrix_init_w_array(ST_ab,3,3,ptrST_ab);

      /* AA = F_I Grad(del u) */      
      Matrix_AxB(AA,1.0,0.0,F_I,0,ST_ab,0);        
                    
      double trAA = 0;
      Matrix_trace(AA,trAA);
      
      for(int w=0; w<nne; w++)
      {
        for(int g=0; g<nsd; g++)
        {
          const double * const ptrST_wg = &ST[idx_4_gen(w,g,0,0,
                  nne,nsd,nsd,nsd)];
          Matrix_init_w_array(ST_wg,3,3,ptrST_wg);        
                                    
          /* BB = F_I Grad(d u) */
          Matrix_AxB(BB,1.0,0.0,F_I,0,ST_wg,0);
          
          double trBB = 0;
          Matrix_trace(BB,trBB);
                          
          /* CC = BB*AA */
          Matrix_AxB(CC,1.0,0.0,AA,0,BB,0);
          
          double trCC = 0;
          Matrix_trace(CC,trCC);
          
          const int K_idx = idx_K(a,b,w,g,nne,nsd);
          K[K_idx] += -dt_alpha_1_minus_alpha*fe->detJxW*Pn*Jn*(trAA*trBB - trCC);          
        }
      }
    }
  }
  Matrix_cleanup(F_I);
  Matrix_cleanup(ST_ab);
  Matrix_cleanup(ST_wg);
  Matrix_cleanup(AA);
  Matrix_cleanup(BB);
  Matrix_cleanup(CC);
}
void add_3F_Kuu_ip(double *K,
        int nne, double *ST, double *F, double jj, double wt, double Pn, double Tn,
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

void add_3F_Kut_ip(double *K,
        int nne, int nVol, double *ST, double *F, double jj, double wt, double *Nt,
        double dt_alpha_1_minus_alpha)
{
	int nsd = 3;
  return;
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

void add_3F_Ktu_ip(double *K,
        int nne, int nVol,
        double *ST, double *F, double jj, double wt, double *Nt,
        double dt_alpha_1_minus_alpha)
{
	int nsd = 3;
  return;
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
      symmetric_part(sAA,AA,3);
      // sAA:S
      double sAAS = cblas_ddot(9,sAA,1,S,1);            
      /* C_IFdu = C_I:AA*/
      double C_IFdu = cblas_ddot(9,C_I,1,AA,1);      
            
      fu[a*nsd + b] += jj*wt*(sAAS + Jn*Pn*C_IFdu);
    }
  }
  
  free(AA);
  free(sAA);
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
        double *r_e)
{ 
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
    P[0] =  eps[ii].d_T[0];	        
              
  int ndofe = nne*ndofn;
  memset(Ks,0,ndofe*ndofe*sizeof(double));
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.0*(1.0-2.0*hommat[mat].nu));
          
  d2UdJ2FuncPtr UPP = getD2UdJ2Func(1, &hommat[mat]);  
  devStressFuncPtr Stress = getDevStressFunc(1,&hommat[mat]);
  matStiffFuncPtr Stiffness = getMatStiffFunc(1,&hommat[mat]);

  int err = 0;
  
  double *L;
  L = aloc1(81);
  Matrix(double) F,C,S;
  Matrix(double) Kuu_add;
  Matrix(double) Kuu,Kut,Kup;
  Matrix(double) Ktu,Ktt,Ktp;
  Matrix(double) Kpu,Kpt,Kpp;

  Matrix_construct_init(double,F,3,3,0.0);
  Matrix_construct_init(double,C,3,3,0.0);
  Matrix_construct_init(double,S,3,3,0.0);
  
  Matrix_construct_init(double,Kuu_add,nne*nsd,nne*nsd,0.0);
  Matrix_construct_init(double,Kuu,nne*nsd,nne*nsd,0.0);
  Matrix_construct_init(double,Kut,nne*nsd,nVol   ,0.0);
  Matrix_construct_init(double,Kup,nne*nsd,npres  ,0.0);
  Matrix_construct_init(double,Ktu,nVol   ,nne*nsd,0.0);
  Matrix_construct_init(double,Ktt,nVol   ,nVol   ,0.0);
  Matrix_construct_init(double,Ktp,nVol   ,npres  ,0.0);    
  Matrix_construct_init(double,Kpu,npres  ,nne*nsd,0.0);
  Matrix_construct_init(double,Kpt,npres  ,nVol   ,0.0);
  Matrix_construct_init(double,Kpp,npres  ,npres  ,0.0); 

// \/this is for test ////////////////////////////////////////////////////////////			  
double *Kuu_ = aloc1(nne*nsd*nne*nsd);
memset(Kuu_,0,ndofe*ndofe*sizeof(double));        

double *Kup_ = aloc1(nne*nsd*npres);
memset(Kup_,0,nne*nsd*npres*sizeof(double));   

double *Ktt_ = aloc1(nVol*nVol);
memset(Ktt_,0,nVol*nVol*sizeof(double));  

double *Ktp_ = aloc1(nVol*npres);
memset(Ktp_,0,nVol*npres*sizeof(double)); 

        double Fn[9];
        Fn[0] = Fn[4] = Fn[8] = 1.0;
        Fn[1] = Fn[2] = Fn[3] = Fn[5] = Fn[6] = Fn[7] = 0.0; 
        double F_I[9]; 
   
// /\this is for test ////////////////////////////////////////////////////////////

  MPI_Comm mpi_comm = MPI_COMM_WORLD;    
  int myrank;
  MPI_Comm_rank(mpi_comm,&myrank);

  FEMLIB fe;
  Matrix(double) xe;  
  Matrix_construct_init(double,xe,nne,3,0.0);

  for(int a=0; a<nne; a++)
  {
    Mat_v(xe, a+1, 1) = x[a];  
    Mat_v(xe, a+1, 2) = y[a];  
    Mat_v(xe, a+1, 3) = z[a];  
  }        

  int itg_order = nne;
  if(nne==4)
    itg_order = nne + 1; 
        
  double ****ST_tensor, *ST;
  double Upp;
  ST_tensor = aloc4(3,3,nsd,nne);
  ST = aloc1(3*3*nsd*nne);
    
  double **F_mat;
  
  F_mat = aloc2(3,3);
  
  Matrix(double) Np, Nt;
  Matrix_construct_init(double,Np,npres,1,0.0);
  Matrix_construct_init(double,Nt,nVol, 1,0.0);    
                         
  FEMLIB_initialization(&fe, itg_order, 1, nne);
  FEMLIB_set_element(&fe, xe, ii);
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);  
    FEMLIB_elem_shape_function(&fe, ip,npres,Np);
    FEMLIB_elem_shape_function(&fe, ip,nVol, Nt);
                    
    double Tn = 0.0; 
    double Pn = 0.0;
    
    for(int a=0; a<nVol; a++)
      Tn += Mat_v(Nt,a+1,1)*eps[ii].T[a*3+0];
    
    for(int a=0; a<npres; a++)
      Pn += Mat_v(Np,a+1,1)*P[a];       
            
    shape_tensor(nne,nsd,fe.temp_v.N_x.m_pdata,fe.temp_v.N_y.m_pdata,fe.temp_v.N_z.m_pdata,ST_tensor);
    shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
                    
    def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u,F_mat);
    mat2array(F.m_pdata,CONST_2(double) F_mat,3,3);

    Matrix_AxB(C,1.0,0.0,F,1,F,0);            
    UPP(Tn,&hommat[mat],&Upp);
    Upp = Upp*kappa;

    Stress(C.m_pdata,&hommat[mat],S.m_pdata);       
    Stiffness(C.m_pdata,&hommat[mat],L);     
                       
    add_3F_Kuu_ip(Kuu.m_pdata,nne ,ST,F.m_pdata,fe.detJ,fe.temp_v.w_ip,Pn,Tn,-1.0);    
    add_3F_Kup_ip(Kup.m_pdata,nne ,npres, ST,F.m_pdata,fe.detJ,fe.temp_v.w_ip,Np.m_pdata,-1.0);
    add_3F_Kpu_ip(Kpu.m_pdata,nne ,npres, ST,F.m_pdata,fe.detJ,fe.temp_v.w_ip,Np.m_pdata,-1.0);        
    add_3F_Kut_ip(Kut.m_pdata,nne ,nVol,  ST,F.m_pdata,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,-1.0);
    add_3F_Ktu_ip(Ktu.m_pdata,nne ,nVol,  ST,F.m_pdata,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,-1.0);        
    add_3F_Ktp_ip(Ktp.m_pdata,nVol,npres,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Np.m_pdata,-1.0);
    add_3F_Kpt_ip(Kpt.m_pdata,nVol,npres,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Np.m_pdata,-1.0);				
    add_3F_Ktt_ip(Ktt.m_pdata,nVol,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Upp,-1.0); 

    add_3F_Kuu_ip_disp(Kuu_add.m_pdata, &fe, ST, F, S, L, -1.0);        				    
    
    inverse(F.m_pdata,3,F_I);
    double Jn = det3x3(F.m_pdata);
//    HW_Kuu_at_ip(Kuu.m_pdata,nne,nne,ST,Fn,F.m_pdata,F_I,1.0,Tn,Jn,Pn,S.m_pdata,L,fe.detJ,fe.temp_v.w_ip,0);    
//    HW_Kup_at_ip(Kup.m_pdata,nne,nne,npres,Np.m_pdata,ST,F.m_pdata,1.0,Tn,Jn,fe.detJ,fe.temp_v.w_ip,0);   
//    HW_Kup_at_ip(Kpu.m_pdata,nne,nne,npres,Np.m_pdata,ST,F.m_pdata,1.0,Tn,Jn,fe.detJ,fe.temp_v.w_ip,0);       
    HW_Ktt_at_ip(Ktt_,nVol,Nt.m_pdata,Tn,Jn,kappa,Upp/kappa,fe.detJ,fe.temp_v.w_ip);             
    HW_Kpt_at_ip(Ktp_,npres,Np.m_pdata,nVol,Nt.m_pdata,Tn,Jn,fe.detJ,fe.temp_v.w_ip);
  } 
  
  for(int a=0; a<nne; a++)
  {
    for(int b=1; b<=nsd; b++)
    {
      for(int c=0; c<nne; c++)
      {
        for(int d=1; d<=nsd; d++)
          Mat_v(Kuu, a*nsd+b, c*nsd+d) += Mat_v(Kuu_add, a*nsd+b, c*nsd+d);         
      }
    }
  }

// this is for test ////////////////////////////////////////////////////////////

//  for(int a=0; a<ndofe*ndofe; a++)
//    printf("%e, %e, %e\n", Kuu_[a], Kuu.m_pdata[a], Kuu_[a]-Kuu.m_pdata[a]);
//
//  for(int a=0; a<ndofe*npres; a++)
//    printf("%e, %e, %e\n", Kup_[a], Kup.m_pdata[a], Kup_[a]-Kup.m_pdata[a]);
//
//  for(int a=0; a<nVol*nVol; a++)
//    printf("%e, %e, %e\n", Ktt_[a], Ktt.m_pdata[a], Ktt_[a]-Ktt.m_pdata[a]);

//  for(int a=0; a<nVol*npres; a++)
//    printf("%e, %e, %e\n", Ktp_[a], Ktp.m_pdata[a], Ktp_[a]-Ktp.m_pdata[a]);    

dealoc1(Kuu_);
dealoc1(Kup_);
dealoc1(Ktt_);
dealoc1(Ktp_);    
// this is for test ////////////////////////////////////////////////////////////  

          
  dealoc4(ST_tensor,3,3,nsd);
  free(ST);
    
  dealoc2(F_mat,3);
  condense_K_out(Ks,nne,nsd,npres,nVol,
                 Kuu.m_pdata,Kut.m_pdata,Kup.m_pdata,
                 Ktu.m_pdata,Ktt.m_pdata,Ktp.m_pdata,
                 Kpu.m_pdata,Kpt.m_pdata,Kpp.m_pdata);
  
  FEMLIB_destruct(&fe);
  Matrix_cleanup(xe); 

  Matrix_cleanup(F);
  Matrix_cleanup(C);
  Matrix_cleanup(S);  
  Matrix_cleanup(Kuu_add);
  Matrix_cleanup(Kuu);
  Matrix_cleanup(Kut);
  Matrix_cleanup(Kup);
  Matrix_cleanup(Ktu);
  Matrix_cleanup(Ktt);
  Matrix_cleanup(Ktp);
  Matrix_cleanup(Kpu);
  Matrix_cleanup(Kpt);
  Matrix_cleanup(Kpp);
  
  Matrix_cleanup(Np);
  Matrix_cleanup(Nt);   
  
  free(P); free(u);    
}

void stiffmat_3f_w_inertia_el(double *Ks,
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
        double alpha, double *r_n, double *r_e)
{   
  int ndofe = nne*ndofn;
  memset(Ks,0,ndofe*ndofe*sizeof(double));
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  
  double rho = hommat[mat].density;

  double alpha_1 = 1.0 - alpha;
  double alpha_2 = alpha;
  double dt_alpha_1_minus_alpha = dt*alpha_1*alpha_2;
          
  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UPP = getD2UdJ2Func(1, &hommat[mat]);

  int err = 0;
  
  Matrix(double) Kuu_I;
  Matrix(double) Kuu_add;
  Matrix(double) Kuu,Kut,Kup;
  Matrix(double) Ktu,Ktt,Ktp;
  Matrix(double) Kpu,Kpt,Kpp;

  Matrix_construct_init(double,Kuu_I,nne*nsd,nne*nsd,0.0);
  Matrix_construct_init(double,Kuu_add,nne*ndofn,nne*ndofn,0.0);
  Matrix_construct_init(double,Kuu,nne*nsd,nne*nsd,0.0);
  Matrix_construct_init(double,Kut,nne*nsd,nVol   ,0.0);
  Matrix_construct_init(double,Kup,nne*nsd,npres  ,0.0);
  Matrix_construct_init(double,Ktu,nVol   ,nne*nsd,0.0);
  Matrix_construct_init(double,Ktt,nVol   ,nVol   ,0.0);
  Matrix_construct_init(double,Ktp,nVol   ,npres  ,0.0);    
  Matrix_construct_init(double,Kpu,npres  ,nne*nsd,0.0);
  Matrix_construct_init(double,Kpt,npres  ,nVol   ,0.0);
  Matrix_construct_init(double,Kpp,npres  ,npres  ,0.0);  

//needs to be verified-------------------------------------->
  DISP_stiffmat_el(Kuu_add.m_pdata,ii,ndofn,nne,x,y,z,elem,
    hommat,nod,node,eps,sig,sup,r_e);                       
//needs to be verified<-------------------------------------
  FEMLIB fe;
  Matrix(double) xe;  
  Matrix_construct_init(double,xe,nne,3,0.0);

//get nodal values---------------------------------------->  
  double *u, *P, *r0, *r_mid;
  u = aloc1(nne*nsd);        
  P = aloc1(npres);
  r0 = aloc1(ndofe);
  r_mid = aloc1(ndofe);
  
  for (int I=0;I<nne;I++)
  {
    for(int J=0; J<ndofn; J++)
      r0[I*ndofn + J] = r_n[nod[I]*ndofn + J];
  }

  mid_point_rule(r_mid, r0, r_e, alpha, ndofe);

  for (int I=0;I<nne;I++)
  {
    for(int J=0; J<nsd; J++)
      u[I*nsd + J] = r_mid[I*ndofn + J];
    if(npres==nne)  
      P[I] = r_mid[I*ndofn + nsd];
  }
  
  if(npres==1)          
      P[0] =  alpha_1*eps[ii].d_T[1] + alpha_2*eps[ii].d_T[0]; 
      
//get nodal values<----------------------------------------

//start integration----------------------------------------> 
  for(int a=0; a<nne; a++)
  {
    Mat_v(xe, a+1, 1) = x[a];  
    Mat_v(xe, a+1, 2) = y[a];  
    Mat_v(xe, a+1, 3) = z[a];  
  }        

  int itg_order = nne;
  if(nne==4)
    itg_order = nne + 1; 
        
  double ****ST_tensor, *ST;
  double Upp;
  ST_tensor = aloc4(3,3,nsd,nne);
  ST = aloc1(3*3*nsd*nne);
    
  double F[9];  
  double **F_mat;
  
  F_mat = aloc2(3,3);
  
  Matrix(double) Np, Nt;
  Matrix_construct_init(double,Np,npres,1,0.0);
  Matrix_construct_init(double,Nt,nVol, 1,0.0);    
                         
  FEMLIB_initialization(&fe, itg_order, 1, nne);
  FEMLIB_set_element(&fe, xe, ii);
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);  
    FEMLIB_elem_shape_function(&fe, ip,npres,Np);
    FEMLIB_elem_shape_function(&fe, ip,nVol, Nt);
              
    if(npres==1)
      Mat_v(Nt,1,1) = 1.0;
    
    if(nVol==1)
      Mat_v(Np,1,1) = 1.0;  
      
    double Tn = 0.0; 
    double Pn = 0.0;

    for(int a=0; a<nVol; a++)
      Tn += Mat_v(Nt,a+1,1)*(alpha_1*eps[ii].T[a*3+1] + alpha_2*eps[ii].T[a*3+0]);
       
    for(int a=0; a<npres; a++)
      Pn += Mat_v(Np,a+1,1)*P[a];      
                        
    UPP(Tn,&hommat[mat],&Upp);
    Upp = Upp*kappa;
    
    shape_tensor(nne,nsd,fe.temp_v.N_x.m_pdata,fe.temp_v.N_y.m_pdata,fe.temp_v.N_z.m_pdata,ST_tensor);
    shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
                    
    def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u,F_mat);
    mat2array(F,CONST_2(double) F_mat,3,3);        
                       
    add_3F_Kuu_ip(Kuu.m_pdata,nne ,ST,F,fe.detJ,fe.temp_v.w_ip,Pn,Tn,dt_alpha_1_minus_alpha);    
    add_3F_Kup_ip(Kup.m_pdata,nne ,npres, ST,F,fe.detJ,fe.temp_v.w_ip,Np.m_pdata,dt_alpha_1_minus_alpha);
    add_3F_Kpu_ip(Kpu.m_pdata,nne ,npres, ST,F,fe.detJ,fe.temp_v.w_ip,Np.m_pdata,dt_alpha_1_minus_alpha);        
    add_3F_Kut_ip(Kut.m_pdata,nne ,nVol,  ST,F,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,dt_alpha_1_minus_alpha);
    add_3F_Ktu_ip(Ktu.m_pdata,nne ,nVol,  ST,F,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,dt_alpha_1_minus_alpha);        
    add_3F_Ktp_ip(Ktp.m_pdata,nVol,npres,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Np.m_pdata,dt_alpha_1_minus_alpha);
    add_3F_Kpt_ip(Kpt.m_pdata,nVol,npres,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Np.m_pdata,dt_alpha_1_minus_alpha);				
    add_3F_Ktt_ip(Ktt.m_pdata,nVol,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Upp,dt_alpha_1_minus_alpha);
    
    for(long a = 0; a<nne; a++)
    {
      for(long c=0; c<nne; c++)
      {
        for(long b=1; b<=nsd; b++)
          Mat_v(Kuu_I,a*nsd+b,c*nsd+b) += rho/dt*Mat_v(fe.N,a+1,1)*Mat_v(fe.N,c+1,1)*fe.detJxW;
      }
	  }             				
  } 
  
  for(int a=0; a<nne; a++)
  {
    for(int b=1; b<=nsd; b++)
    {
      for(int c=0; c<nne; c++)
      {
        for(int d=1; d<=nsd; d++)
          Mat_v(Kuu, a*nsd+b, c*nsd+d) -= Mat_v(Kuu_I, a*nsd+b, c*nsd+d)+dt_alpha_1_minus_alpha*Mat_v(Kuu_add, a*ndofn+b, c*ndofn+d);         
      }
    }
  }
          
  dealoc4(ST_tensor,3,3,nsd);
  free(ST);
    
  dealoc2(F_mat,3);
  condense_K_out(Ks,nne,nsd,npres,nVol,
                 Kuu.m_pdata,Kut.m_pdata,Kup.m_pdata,
                 Ktu.m_pdata,Ktt.m_pdata,Ktp.m_pdata,
                 Kpu.m_pdata,Kpt.m_pdata,Kpp.m_pdata);
  
  FEMLIB_destruct(&fe); 

  Matrix_cleanup(Kuu_I);
  Matrix_cleanup(Kuu_add);
  Matrix_cleanup(Kuu);
  Matrix_cleanup(Kut);
  Matrix_cleanup(Kup);
  Matrix_cleanup(Ktu);
  Matrix_cleanup(Ktt);
  Matrix_cleanup(Ktp);
  Matrix_cleanup(Kpu);
  Matrix_cleanup(Kpt);
  Matrix_cleanup(Kpp);   
  
  free(r_mid);
  free(r0);
  free(u);
  free(P);    
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

  Matrix(double) F,S,C;
  Matrix_construct_init(double,F,3,3,0.0);
  Matrix_construct_init(double,C,3,3,0.0);  
  Matrix_construct_init(double,S,3,3,0.0);  
  
  //Matrix(double) fu_add;
  Matrix(double) fu,Kut;
  Matrix(double) fp,Kpt;
  Matrix(double) ft,Ktt;
  Matrix(double) Kup,Ktp;
  
  int ndofe = nne*ndofn;

//  Matrix_construct_init(double,fu_add, nne*ndofn,1,0.0);
  Matrix_construct_init(double,fu, nne*nsd,1,    0.0);
  Matrix_construct_init(double,fp, npres,  1,    0.0);
  Matrix_construct_init(double,ft, nVol,   1,    0.0);
  Matrix_construct_init(double,Kut,nne*nsd,nVol, 0.0);
  Matrix_construct_init(double,Kpt,npres,  nVol, 0.0);
  Matrix_construct_init(double,Ktt,nVol,   nVol, 0.0);    
  Matrix_construct_init(double,Kup,nne*nsd,npres,0.0);
  Matrix_construct_init(double,Ktp,nVol,   npres,0.0);

//residuals_disp_el(fu_add.m_pdata,ii,ndofn,nne,nsd,x,y,z,elem,
//	         hommat,node,r_e);
//  DISP_resid_el(fu_add.m_pdata,ii,ndofn,nne,x,y,z,elem,
//	         hommat,nod,node,eps,sig,sup,r_e);  


// \/this is for test ////////////////////////////////////////////////////////////			  
//
  double *Ru = aloc1(nne*nsd);
  double *Rp = aloc1(nne*nsd);
  double *Rt = aloc1(nne*nsd);    
  memset(Ru,0,nne*nsd*sizeof(double));  
  memset(Rp,0,npres*sizeof(double));        
  memset(Rt,0,nVol*sizeof(double));    
  
  double F_I[9];

  double Fn[9];
  Fn[0] = Fn[4] = Fn[8] = 1.0;
  Fn[1] = Fn[2] = Fn[3] = Fn[5] = Fn[6] = Fn[7] = 0.0;  
  
  
// 
// /\this is for test ////////////////////////////////////////////////////////////			


  FEMLIB fe;
  Matrix(double) xe;  
  Matrix_construct_init(double,xe,nne,3,0.0);

  for(int a=0; a<nne; a++)
  {
    Mat_v(xe, a+1, 1) = x[a];  
    Mat_v(xe, a+1, 2) = y[a];  
    Mat_v(xe, a+1, 3) = z[a];  
  }        

  int itg_order = nne;
  if(nne==4)
    itg_order = nne + 1; 
    
    
  double ****ST_tensor, *ST;
   
  ST_tensor = aloc4(3,3,nsd,nne);
  ST = aloc1(3*3*nsd*nne);

  double **F_mat;

  F_mat = aloc2(3,3);
                        
  Matrix(double) Np, Nt;
  Matrix_construct_init(double,Np,npres,1,0.0);
  Matrix_construct_init(double,Nt,nVol, 1,0.0);    
                         
  FEMLIB_initialization(&fe, itg_order, 1, nne);
  FEMLIB_set_element(&fe, xe, ii);
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip);  
    FEMLIB_elem_shape_function(&fe, ip,npres,Np);
    FEMLIB_elem_shape_function(&fe, ip,nVol, Nt);
              
    double Tn = 0.0; 
    double Pn = 0.0;
    
    for(int a=0; a<nVol; a++)
      Tn += Mat_v(Nt,a+1,1)*eps[ii].T[a*3+0];
    
    for(int a=0; a<npres; a++)
      Pn += Mat_v(Np,a+1,1)*P[a];
      
                
    double Upp = 0.0;
    double Up = 0.0;
    UPP(Tn,&hommat[mat],&Upp);
    Upp = Upp*kappa;

    UP(Tn, &hommat[mat], &Up);
    Up = Up*kappa;
    
    shape_tensor(nne,nsd,fe.temp_v.N_x.m_pdata,fe.temp_v.N_y.m_pdata,fe.temp_v.N_z.m_pdata,ST_tensor);
    shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
        
    def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u,F_mat);
    mat2array(F.m_pdata,CONST_2(double) F_mat,3,3);        

    Matrix_AxB(C,1.0,0.0,F,1,F,0);
    Stress(C.m_pdata,&hommat[mat],S.m_pdata);
    
		add_3F_Kpt_ip(Kpt.m_pdata,nVol,npres,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Np.m_pdata,-1.0);	
		add_3F_Ktt_ip(Ktt.m_pdata,nVol,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Upp,-1.0);	

    if(npres==nne)
    {
      add_3F_Kut_ip(Kut.m_pdata,nne,nVol,ST,F.m_pdata,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,-1.0);
    }	
    
		if(npres==1)
		{
      add_3F_Kup_ip(Kup.m_pdata,nne,npres,ST,F.m_pdata,fe.detJ,fe.temp_v.w_ip,Np.m_pdata,-1.0);
	  	add_3F_Ktp_ip(Ktp.m_pdata,nVol,npres,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Np.m_pdata,-1.0);			
	  }	
		resid_w_inertia_Rp_ip(fp.m_pdata, npres, F.m_pdata, fe.detJ,fe.temp_v.w_ip, Np.m_pdata, Tn);
//		resid_w_inertia_Rt_ip(ft.m_pdata, nVol, fe.detJ,fe.temp_v.w_ip, Nt.m_pdata, Pn, Up);
		resid_w_inertia_Ru_ip(fu.m_pdata, nne, ST, F.m_pdata, S.m_pdata, fe.detJ,fe.temp_v.w_ip, Pn);


// \/this is for test ////////////////////////////////////////////////////////////			   
   double Jn = det3x3(F.m_pdata);
   inverse(F.m_pdata,3,F_I);  

    HW_Ru_at_ip(Ru,nne,nne,ST,Fn,F.m_pdata,F_I,
    		1.0,0.0,Jn,Tn,S.m_pdata,Pn,fe.detJ,fe.temp_v.w_ip,0);
    HW_Rp_at_ip(Rp,npres,Np.m_pdata,1.0,Jn,1.0,Tn,fe.detJ,fe.temp_v.w_ip);
    HW_Rt_at_ip(ft.m_pdata,nVol,Nt.m_pdata,Tn,Jn,Pn,kappa,Up/kappa,fe.detJ,fe.temp_v.w_ip);	
// /\this is for test ////////////////////////////////////////////////////////////												
 
  }


// this is for test ////////////////////////////////////////////////////////////
//for(int a=0; a<ndofe; a++)
//  printf("%e %e %e\n", Ru[a], fu.m_pdata[a], Ru[a]-fu.m_pdata[a]);
  
//for(int a=0; a<nVol; a++)
//  printf("%e %e %e\n", Rt[a], ft.m_pdata[a], Rt[a]-ft.m_pdata[a]);
//  
//for(int a=0; a<npres; a++)
//  printf("%e %e %e\n", Rp[a], fp.m_pdata[a], Rp[a]-fp.m_pdata[a]);  
    
  dealoc1(Ru);
  dealoc1(Rp);
  dealoc1(Rt);  
// this is for test ////////////////////////////////////////////////////////////  


  dealoc4(ST_tensor,3,3,nsd);
  free(ST);
  
  dealoc2(F_mat,3);  
  
  condense_F_out(f,nne,nsd,npres,nVol,fu.m_pdata,ft.m_pdata,fp.m_pdata,
                    Kut.m_pdata,Kup.m_pdata,Ktp.m_pdata,Ktt.m_pdata,Kpt.m_pdata);

  Matrix_cleanup(F);
  Matrix_cleanup(C);
  Matrix_cleanup(S); 
                        
  FEMLIB_destruct(&fe);                   
  Matrix_cleanup(xe);
  Matrix_cleanup(Np);  
  Matrix_cleanup(Nt);
  
  Matrix_cleanup(fu);
  Matrix_cleanup(fp); 
  Matrix_cleanup(ft); 
  Matrix_cleanup(Kut);
  Matrix_cleanup(Kpt);
  Matrix_cleanup(Ktt);
  Matrix_cleanup(Kup);
  Matrix_cleanup(Ktp); 
  
  free(P); free(u);   
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

void evaluate_PT_el(const int ii,
        const int ndofn,
        const int nne,
        const int npres,
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const NODE *node,
        const double *r_e,
        const double *u,
        const double *P,
        double dt,
        SIG *sig,
        EPS *eps)
{
  Matrix(double) uu;
	Matrix_construct(double, uu);	  
	
	Matrix_init_w_array(uu, nne*ndofn, 1, u);	  
  
  
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));  
  
  
  double volume;
  if(nne == 4){ /* linear tet */
    volume = Tetra_V(x,y,z);
  } else if(nne == 10){ /* quadradic tet */
    volume = Tetra_qv_V(nne,ndofn,x,y,z);
  } else if(nne == 8){ /* trilinear hex */
    volume = Hexa_V(x,y,z);
  } else volume = 0.0;  
    			    
  
  dUdJFuncPtr UP = getDUdJFunc(1, &hommat[mat]);
  d2UdJ2FuncPtr UPP = getD2UdJ2Func(1, &hommat[mat]);
  devStressFuncPtr Stress = getDevStressFunc(1,&hommat[mat]);
  
  const int nsd = 3;
  const int nVol = N_VOL_TF;
  int err = 0;  
  int ndofe = nne*ndofn;
	
  Matrix(double) fp;
  Matrix(double) ft;
  
//  Matrix_construct_init(double,fu_add, nne*ndofn,1,0.0);
  Matrix_construct_init(double,fp, npres,  1,    0.0);
  Matrix_construct_init(double,ft, nVol,   1,    0.0);
		
  Matrix(double) Ktu,Ktt,Ktp;
  Matrix(double) Kpu,Kpt;
  
  Matrix_construct_init(double,Ktu,nVol   ,nne*nsd,0.0);
  Matrix_construct_init(double,Ktt,nVol   ,nVol   ,0.0);
  Matrix_construct_init(double,Ktp,nVol   ,npres  ,0.0);    
  Matrix_construct_init(double,Kpu,npres  ,nne*nsd,0.0);
  Matrix_construct_init(double,Kpt,npres  ,nVol   ,0.0);
  
  double F_I[9];

  double Fn[9];
  Fn[0] = Fn[4] = Fn[8] = 1.0;
  Fn[1] = Fn[2] = Fn[3] = Fn[5] = Fn[6] = Fn[7] = 0.0;    
	  
  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int itg_order = nne;
  if(nne==4)
    itg_order = nne + 1;  
  int_point(itg_order,&npt_z);
  
  double *int_pt_ksi, *int_pt_eta, *int_pt_zet, *weights;
  int_pt_ksi = aloc1(npt_z);
  int_pt_eta = aloc1(npt_z);
  int_pt_zet = aloc1(npt_z);
  weights = aloc1(npt_z);
  
  /* allocate space for the shape functions, derivatives etc */
  double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
 
	Na = aloc1(nne); 
  Np = aloc1(npres);
  Nt = aloc1(nVol);
  N_x = aloc1(nne);
  N_y = aloc1(nne);
  N_z = aloc1(nne);
  ST_tensor = aloc4(3,3,nsd,nne);
  ST = aloc1(3*3*nsd*nne);
  
  /*=== INTEGRATION LOOP ===*/
  integrate(itg_order,&npt_x,&npt_y,&npt_z,
          int_pt_ksi,int_pt_eta,int_pt_zet,
          weights);
  
  double **F_mat, *F;
  
  F = aloc1(9);
  F_mat = aloc2(3,3);

  int ip = 1;
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
        
        double Tn = 0.0; // 1: n-1+alpha
        double Pn = 0.0;
       
        if(npres==1)
          Nt[0] = 1.0;
        
        if(nVol==1)
          Np[0] = 1.0;
                 
        for(int a=0; a<nVol; a++)
          Tn += Nt[a]*eps[ii].T[a*3+0];          
        
        for(int a=0; a<npres; a++)
          Pn += Np[a]*P[a];
          
        double Upp = 0.0;
        double Up = 0.0;          
    
        UPP(Tn,&hommat[mat],&Upp);
        Upp = Upp*kappa;

        UP(Tn, &hommat[mat], &Up);
        Up = Up*kappa;
                
        shape_tensor(nne,ndn,N_x,N_y,N_z,ST_tensor);
        shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
        
        def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u,F_mat);
        mat2array(F,CONST_2(double) F_mat,3,3);
        
        add_3F_Ktu_ip(Ktu.m_pdata,nne,nVol,  ST,F,detJ,wt,Nt,-1.0);        
				add_3F_Ktp_ip(Ktp.m_pdata,nVol,npres,detJ,wt,Nt,Np,-1.0);
				add_3F_Ktt_ip(Ktt.m_pdata,nVol,detJ,wt,Nt,Upp,-1.0);	
        add_3F_Kpu_ip(Kpu.m_pdata,nne,npres, ST,F,detJ,wt,Np,-1.0);        
				add_3F_Kpt_ip(Kpt.m_pdata,nVol,npres,detJ,wt,Nt,Np,-1.0);		
double Jn = det3x3(F);						
//    HW_Kup_at_ip(Kpu.m_pdata,nne,nne,npres,Np,ST,F,1.0,Tn,Jn,detJ,wt,0); 
//				resid_w_inertia_Rt_ip(ft.m_pdata, nVol, detJ, wt, Nt, Pn, Up);							
				resid_w_inertia_Rp_ip(fp.m_pdata, npres, F, detJ, wt, Np, Tn);
    HW_Rt_at_ip(ft.m_pdata,nVol,Nt,Tn,Jn,Pn,kappa,Up/kappa,detJ,wt);	
				
//        update_3f_state_variables_ip(ii,ip,elem,hommat,sig,eps,F,Pn,Tn,volume,detJ*wt);
        ip++;				
								
      }
    }
  }
    	  	  	 
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


  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  int myrank;
  MPI_Comm_rank(mpi_comm,&myrank);	
  
  Matrix(double) theta, KptI;
  Matrix_construct_init(double,theta, nVol,1,0.0);
  Matrix_construct_init(double,KptI, nVol,nVol,0.0);  
  Matrix_inv(Kpt,KptI);
  
  Matrix_AxB(fp, 1.0, 1.0, Kpu, 0, uu, 0);        
  Matrix_AxB(theta, -1.0, 0.0, KptI, 0, fp, 0); 
  
	for(int a = 0; a<nVol; a++)
	  eps[ii].T[a*3+0] += theta.m_pdata[a];

  Matrix(double) press, KtpI;
  Matrix_construct_init(double,press, npres,1,0.0);
  Matrix_construct_init(double,KtpI, nVol,npres,0.0);  
  Matrix_inv(Ktp,KtpI);
  
  Matrix_AxB(ft, 1.0, 1.0, Ktt, 0, theta, 0);
  Matrix_AxB(press, -1.0, 0.0, KtpI, 0, ft, 0);  
  
  for(int a = 0; a<npres; a++)
    eps[ii].d_T[a*3+0] += press.m_pdata[a]; 
  
  Matrix_cleanup(uu); 
  Matrix_cleanup(fp);  
  Matrix_cleanup(ft);  
  Matrix_cleanup(Ktu);
  Matrix_cleanup(Ktp);
  Matrix_cleanup(Ktt);
  Matrix_cleanup(Kpu);
  Matrix_cleanup(Kpt);
  Matrix_cleanup(theta); 
  Matrix_cleanup(KptI); 
  Matrix_cleanup(press);
  Matrix_cleanup(KtpI);
}

void evaluate_PT_w_inertia_el(const int ii,
        const int ndofn,
        const int nne,
        const int npres,
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const NODE *node,
        const double *u,
        const double *u2, // 2: n+alpha
        const double *u1, // 1: n-1+alpha
        const double *P2, // 1: n+alpha
        const double *P1, // 1: n-1+alpha
        double dt,
        SIG *sig,
        EPS *eps,
        double alpha)
{
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  
  double rho = hommat[mat].density;
  long include_inertia = 1;

  double alpha_1 = 1.0 - alpha;
  double alpha_2 = alpha;
  double dt_alpha_1_minus_alpha = dt*alpha_1*alpha_2;
      
  if(fabs(rho)<1.0e-15)
    include_inertia = 0;

  if(include_inertia==0)
  { 
    dt_alpha_1_minus_alpha = -1.0;
    alpha_1 = 0.0;
    alpha_2 = 1.0;
  }  
  
  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UP = getDUdJFunc(1, &hommat[mat]);  
  UPP = getD2UdJ2Func(1, &hommat[mat]);
  
  const int nsd = 3;
  const int nVol = N_VOL_TF;
//  const int npres = nVol;
  const double Vol = Tetra_V (x,y,z);
  int err = 0;  
  int ndofe = nne*ndofn;

  double *fp  = aloc1(npres);   
  double *fp1 = aloc1(npres);
  double *fp2 = aloc1(npres);  
  
  double *ft  = aloc1(nVol);   
  double *ft1 = aloc1(nVol);
  double *ft2 = aloc1(nVol); 
    
  double *Ktu = aloc1(nne*nsd*nVol);
  double *Ktp = aloc1(npres*nVol);
  double *Ktt = aloc1(nVol*nVol);
  
  double *Kpu = aloc1(npres*nne*nsd);
  double *Kpt = aloc1(npres*nVol);  

  memset(fp, 0,        npres*sizeof(double));      
  memset(fp1,0,        npres*sizeof(double));    
  memset(fp2,0,        npres*sizeof(double));      

  memset(ft, 0,        nVol*sizeof(double));      
  memset(ft1,0,        nVol*sizeof(double));    
  memset(ft2,0,        nVol*sizeof(double));      
  
	memset(Ktu,0,nne*nsd*nVol*sizeof(double));
	memset(Ktp,0,  nVol*npres*sizeof(double));    
	memset(Ktt,0,     nVol*nVol*sizeof(double));	
	
	memset(Kpu,0,npres*nne*nsd*sizeof(double));
	memset(Kpt,0,  npres*nVol*sizeof(double));    

  
  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int itg_order = nne;
  if(nne==4)
    itg_order = nne + 1;  
  int_point(itg_order,&npt_z);
  
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
  integrate(itg_order,&npt_x,&npt_y,&npt_z,
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
       
        if(npres==1)
          Nt[0] = 1.0;
        
        if(nVol==1)
          Np[0] = 1.0;
                 
        for(int a=0; a<nVol; a++)
        {
          Tn1 += Nt[a]*((1.0-alpha)*eps[ii].T[a*3+2] + alpha*eps[ii].T[a*3+1]);
          Tn2 += Nt[a]*((1.0-alpha)*eps[ii].T[a*3+1] + alpha*eps[ii].T[a*3+0]);          
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

        add_3F_Ktu_ip(Ktu,nne,nVol,  ST,F2,detJ,wt,Nt,dt*(1.0-alpha)*alpha);        
				add_3F_Ktp_ip(Ktp,nVol,npres,detJ,wt,Nt,Np,dt*(1.0-alpha)*alpha);
				add_3F_Ktt_ip(Ktt,nVol,detJ,wt,Nt,Upp2,dt*(1.0-alpha)*alpha);	
        add_3F_Kpu_ip(Kpu,nne,npres, ST,F2,detJ,wt,Np,dt*(1.0-alpha)*alpha);        
				add_3F_Kpt_ip(Kpt,nVol,npres,detJ,wt,Nt,Np,dt*(1.0-alpha)*alpha);				

				resid_w_inertia_Rt_ip(ft1, nVol, detJ, wt, Nt, Pn1, Up1);
				resid_w_inertia_Rt_ip(ft2, nVol, detJ, wt, Nt, Pn2, Up2);
							
				resid_w_inertia_Rp_ip(fp1, npres, F1, detJ, wt, Np, Tn1);
				resid_w_inertia_Rp_ip(fp2, npres, F2, detJ, wt, Np, Tn2);			
								
      }
    }
  }
  
  for(int a=0; a<nVol; a++)
  	ft[a] = (1.0 - alpha)*dt*ft2[a] + alpha*dt*ft1[a];

  for(int a=0; a<npres; a++)
  	fp[a] = -(1.0 - alpha)*dt*fp2[a] - alpha*dt*fp1[a];
  	  	  	 
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
    
  free(fp1);
  free(fp2);
  free(ft1);
  free(ft2);
  
  double *theta   = aloc1(nVol);
  double *KptI    = aloc1(nVol*nVol);
  
  inverse(Kpt,nVol,KptI);
  
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nne*nsd,1.0,Kpu,nne*nsd,u,1,1.0,fp,1);  
              
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nVol,-1.0,KptI,nVol,fp,1,0.0,theta,1);    


	for(int a = 0; a<nVol; a++)
	  eps[ii].T[a*3+0] += theta[a];

  double *press   = aloc1(nVol);  
  double *KtpI = aloc1(nVol*nVol);
  
  inverse(Ktp,nVol,KtpI);        
         
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nVol,1.0,Ktt,nVol,theta,1,1.0,ft,1);               
              
//	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//              nVol,1,nne*nsd,1.0,Ktu,nne*nsd,u,1,1.0,ft,1);               

//	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
//              nVol,1,nVol,-1.0,KtpI,nVol,ft,1,0.0,press,1);                                  	                                                          

	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nVol,-1.0,KtpI,nVol,ft,1,0.0,press,1);                                  	                                                          

  for(int a = 0; a<npres; a++)
    eps[ii].d_T[a*3+0] += press[a];      

  free(fp);  
  free(ft);  
  free(Ktu);
  free(Ktp);
  free(Ktt);
  free(Kpu);
  free(Kpt);

  free(theta); 
  free(KptI); 
  free(press);
  free(KtpI);
}

void evaluate_theta_el(const int ii,
        const int ndofn,
        const int nne,
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const NODE *node,
        const double *u,
        const double *P,
        double dt,
        SIG *sig,
        EPS *eps)
{
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  
  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UP = getDUdJFunc(1, &hommat[mat]);  
  UPP = getD2UdJ2Func(1, &hommat[mat]);
  
  const int nsd = 3;
  const int nVol = N_VOL_TF;
  const int npres = nne;
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
  if(nne==4)
    itg_order = nne + 1;  
  int_point(itg_order,&npt_z);
  
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
  integrate(itg_order,&npt_x,&npt_y,&npt_z,
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
        
        add_3F_Ktu_ip(Ktu,nne,nVol,  ST,F,detJ,wt,Nt,-1.0);        
				add_3F_Ktp_ip(Ktp,nVol,npres,detJ,wt,Nt,Np,-1.0);
				add_3F_Ktt_ip(Ktt,nVol,detJ,wt,Nt,Upp,-1.0);	

				resid_w_inertia_Rt_ip(ft, nVol, detJ, wt, Nt, Pn, Up);
								
      }
    }
  }
   	  	  	 
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
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const NODE *node,
        const double *u2, // 2: n+alpha
        const double *u1, // 1: n-1+alpha
        const double *P2, // 1: n+alpha
        const double *P1, // 1: n-1+alpha
        double dt,
        SIG *sig,
        EPS *eps,
        double alpha)
{
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  
  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UP = getDUdJFunc(1, &hommat[mat]);  
  UPP = getD2UdJ2Func(1, &hommat[mat]);
  
  const int nsd = 3;
  const int nVol = N_VOL_TF;
  const int npres = nne;
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
  int itg_order = nne;
  if(nne==4)
    itg_order = nne + 1;  
  int_point(itg_order,&npt_z);
  
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
  integrate(itg_order,&npt_x,&npt_y,&npt_z,
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

        add_3F_Ktu_ip(Ktu,nne,nVol,  ST,F2,detJ,wt,Nt,dt*(1.0-alpha)*alpha);        
				add_3F_Ktp_ip(Ktp,nVol,npres,detJ,wt,Nt,Np,dt*(1.0-alpha)*alpha);
				add_3F_Ktt_ip(Ktt,nVol,detJ,wt,Nt,Upp2,dt*(1.0-alpha)*alpha);	

				resid_w_inertia_Rt_ip(ft1, nVol, detJ, wt, Nt, Pn1, Up1);
				resid_w_inertia_Rt_ip(ft2, nVol, detJ, wt, Nt, Pn2, Up2);
								
      }
    }
  }
  
  for(int a=0; a<nVol; a++)
  	ft[a] = -(1.0 - alpha)*dt*ft2[a] - alpha*dt*ft1[a];
  	  	  	 
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

void update_3f(long ne,
		  long ndofn,
		  long npres,
		  double *d_r,
		  double *r,
		  NODE *node,
		  ELEMENT *elem,
		  HOMMAT *hommat,
		  SUPP sup,
		  EPS *eps,
		  SIG *sig,
		  double dt,
		  double t,
		  MPI_Comm mpi_comm,
		  const PGFem3D_opt *opts,
		  double alpha, double *r_n, double *r_n_1)
{
  const int mat = elem[0].mat[2];
  double rho = hommat[mat].density;
  long include_inertia = 1;
  
  if(fabs(rho)<1.0e-15)
    include_inertia = 0;
            
  int err = 0;
      
  for (int i=0;i<ne;i++)
  {

    int nne = elem[i].toe;
    long *nod = aloc1l (nne);
    elemnodes (i,nne,nod,elem);
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node);

    double *r_e = aloc1 (ndofe);

    double *x,*y,*z;
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);

    long *cn = aloc1l (ndofe);

    nodecoord_total(nne,nod,node,x,y,z);

    /* code numbers on element */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);
    
    /* deformation on element */
    def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
    
  	int nsd = 3;
	  double *r_en;
	
    r_en = aloc1(ndofe);

//    def_elem(cn,ndofe,d_r,elem,node,r_e,sup,0);
    def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
	  
	  if(include_inertia)
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

   			   	
  		double *u1, *u2, *P1, *P2;
  		u1 = aloc1(nne*nsd);
  		u2 = aloc1(nne*nsd);
  		P1 = aloc1(npres);		  			
  		P2 = aloc1(npres);
		
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
        P1[0] = (1.0-alpha)*eps[i].d_T[2] + alpha*eps[i].d_T[1]; 
        P2[0] = (1.0-alpha)*eps[i].d_T[1] + alpha*eps[i].d_T[0];
      }					
  		
  		if(npres==1)
  		  evaluate_PT_w_inertia_el(i,ndofn,nne,npres,x,y,z,elem,hommat,node,r_e,u2,u1,P2,P1,dt,sig,eps,alpha);
  		else
  		  evaluate_theta_w_inertia_el(i,ndofn,nne,x,y,z,elem,hommat,node,u2,u1,P2,P1,dt,sig,eps,alpha); 

	    free(r_en);
	    free(r0);
	    free(r0_);
		  free(r_n_a);
		  free(r_n_1_a);
		  free(P1); free(P2); free(u1); free(u2);
    }
    else
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
        P[0] = eps[i].d_T[0];
  		
  		if(npres==1)
  		  evaluate_PT_el(i,ndofn,nne,npres,x,y,z,elem,hommat,node,r_e,u,P,dt,sig,eps);
  		else
  		  evaluate_theta_el(i,ndofn,nne,x,y,z,elem,hommat,node,u,P,dt,sig,eps);
 
  	  free(u);
  		free(P);
  		
  	}
   	free(nod);
   	free(cn);
    free(x);
    free(y);
    free(z);  	
  	free(r_e);
  	free(r_en);  	
	}	        							  											   		
}
