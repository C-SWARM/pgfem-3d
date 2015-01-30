/* HEADER */
/**
 * AUTHORS:
 * Matthew Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 * Karel Matous, University of Notre Dame, kmatous [at] nd.edu
 */
#include "stiffmat_fd.h"
#include "enumerations.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "elem3d.h"
#include "allocation.h"
#include "PLoc_Sparse.h"
#include "stabilized.h"
#include "stiffmatel_fd.h"
#include "utils.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "displacement_based_element.h"
#include "matice.h"

#include "three_field_element.h"
#include "condense.h"
#include "new_potentials.h"
#include "tensors.h"
#include "cast_macros.h"
#include "index_macros.h"
#include "def_grad.h"

#ifndef PFEM_DEBUG
#define PFEM_DEBUG 0
#endif

#define ndn 3
#define N_VOL_TF 1

static const int periodic = 0;

void stiffmat_w_inertia_el(double *Ks,
         const int ii,
         const int ndofn,
         const int nne,
         const double *x,
         const double *y,
         const double *z,		     
         const ELEMENT *elem,
         const HOMMAT *hommat,
		     const NODE *node, double dt)
{
  int err = 0;
  const int mat = elem[ii].mat[2];
  double rho = hommat[mat].density;
  
  int ndofe = nne*ndofn;

  /* make sure the stiffenss matrix contains all zeros */
  memset(Ks,0,ndofe*ndofe*sizeof(double));

  FEMLIB fe;
  Matrix(double) xe, du;
  
  Matrix_construct_redim(double,xe,nne,3);
  
  for(int a=0; a<nne; a++)
  {
    Mat_v(xe, a+1, 1) = x[a];  
    Mat_v(xe, a+1, 2) = y[a];  
    Mat_v(xe, a+1, 3) = z[a];  
  }        

  int itg_order = nne;
  if(nne==4)
    itg_order = nne + 1; 
    
               
  FEMLIB_initialization(&fe, itg_order, 1, nne);
  FEMLIB_set_element(&fe, xe, ii);      
  for(int ip = 1; ip<=fe.nint; ip++)
  {
    FEMLIB_elem_basis_V(&fe, ip); 
    
    for(long a = 0; a<nne; a++)
    {
      for(long b=0; b<ndofn; b++)
      {
        for(long c=0; c<nne; c++)
        {
          for(long d = 0; d<=ndofn; d++)
          {
            if(b==d)
            {
               const int K_idx = idx_K(a,b,c,d,nne,ndofn);                  
               Ks[K_idx] += rho/dt*Mat_v(fe.N,a+1,1)*Mat_v(fe.N,c+1,1)*fe.detJxW;
	          }
	        } 
	      }		     
      }
	  }
  }
  Matrix_cleanup(xe);  
  FEMLIB_destruct(&fe);
}

void stiffmat_3f_w_inertia_el(double *Ks,
        const int ii,
        double *KsI,
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
        const NODE *node,
        const double *u,
        const double *P,        
        double dt,
        SIG *sig,
        EPS *eps,
        double alpha)
{ 
  int ndofe = nne*ndofn;
  memset(Ks,0,ndofe*ndofe*sizeof(double));
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
  UPP = getD2UdJ2Func(1, &hommat[mat]);

  int err = 0;  
  

  double *Kuu = aloc1(nne*nsd*nne*nsd);
  double *Kut = aloc1(nne*nsd*nVol);
  double *Kup = aloc1(nne*nsd*npres);
  
  double *Ktu = aloc1(nVol*nne*nsd);
  double *Ktt = aloc1(nVol*nVol);
  double *Ktp = aloc1(nVol*npres);
  
  double *Kpu = aloc1(npres*nsd*nne);
  double *Kpt = aloc1(npres*nVol);
  double *Kpp = aloc1(npres*npres);
  
	memset(Kuu,0,nne*nsd*nne*nsd*sizeof(double));  	
	memset(Kut,0,   nne*nsd*nVol*sizeof(double));
	memset(Kup,0,  nne*nsd*npres*sizeof(double));	
	memset(Ktu,0,   nne*nsd*nVol*sizeof(double));			
	memset(Ktt,0,      nVol*nVol*sizeof(double));		
	memset(Ktp,0,     nVol*npres*sizeof(double));		
	memset(Kpu,0,  nne*nsd*npres*sizeof(double));			
	memset(Kpt,0,     nVol*npres*sizeof(double));  
	memset(Kpp,0,    npres*npres*sizeof(double));

  // make sure the stiffenss matrix contains all zeros

  
  // INTEGRATION
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
  
  // allocate space for the shape functions, derivatives etc 
  double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  double Upp;
 
	Na = aloc1(nne); 
  Np = aloc1(npres);
  Nt = aloc1(nVol);
  N_x = aloc1(nne);
  N_y = aloc1(nne);
  N_z = aloc1(nne);
  ST_tensor = aloc4(3,3,nsd,nne);
  ST = aloc1(3*3*nsd*nne);
  
  //=== INTEGRATION LOOP ===
  integrate(itg_order,&npt_x,&npt_y,&npt_z,
          int_pt_ksi,int_pt_eta,int_pt_zet,
          weights);
  
  double **F_mat, *F;
  
  F = aloc1(9);
  F_mat = aloc2(3,3);
  
// \/this is for test ////////////////////////////////////////////////////////////			  
//  double *Kuu_ = aloc1(nne*nsd*nne*nsd);
//  double *Ku_ = aloc1(nne*nsd*nne*nsd);
//  memset(Kuu_,0,ndofe*ndofe*sizeof(double));        
//  memset(Ku_,0,ndofe*ndofe*sizeof(double));  
//  
//  double *Kup_ = aloc1(nne*nsd*npres);
//  memset(Kup_,0,nne*nsd*npres*sizeof(double));   
//  
//  double *Ktt_ = aloc1(nVol*nVol);
//  memset(Ktt_,0,nVol*nVol*sizeof(double));  
//  
//  double *Ktp_ = aloc1(nVol*npres);
//  memset(Ktp_,0,nVol*npres*sizeof(double));    
// /\this is for test ////////////////////////////////////////////////////////////			
  
  int ip = 0;
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
       
        if(npres==1)
          Nt[0] = 1.0;
        
        if(nVol==1)
          Np[0] = 1.0;

        for(int a=0; a<nVol; a++)
          Tn += Nt[a]*(alpha_1*eps[ii].T[a*3+1] + alpha_2*eps[ii].T[a*3+0]);
           
        for(int a=0; a<npres; a++)
          Pn += Np[a]*P[a];
          
        UPP(Tn,&hommat[mat],&Upp);
        Upp = Upp*kappa;
        shape_tensor(nne,ndn,N_x,N_y,N_z,ST_tensor);
        shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
        
        def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u,F_mat);
        mat2array(F,CONST_2(double) F_mat,3,3);        
       
        add_3F_Kuu_ip(Kuu,nne,ST,F,detJ,wt,Pn,Tn,dt_alpha_1_minus_alpha);
    
        if(include_inertia==0 || (alpha>1.0e-15 && fabs(alpha-1.0)>1.0e-15))
        {       
          add_3F_Kup_ip(Kup,nne,npres, ST,F,detJ,wt,Np,dt_alpha_1_minus_alpha);
          add_3F_Kpu_ip(Kpu,nne,npres, ST,F,detJ,wt,Np,dt_alpha_1_minus_alpha);        
          add_3F_Kut_ip(Kut,nne,nVol,  ST,F,detJ,wt,Nt,dt_alpha_1_minus_alpha);
          add_3F_Ktu_ip(Ktu,nne,nVol,  ST,F,detJ,wt,Nt,dt_alpha_1_minus_alpha);        
				  add_3F_Ktp_ip(Ktp,nVol,npres,detJ,wt,Nt,Np,dt_alpha_1_minus_alpha);
				  add_3F_Kpt_ip(Kpt,nVol,npres,detJ,wt,Nt,Np,dt_alpha_1_minus_alpha);				
				  add_3F_Ktt_ip(Ktt,nVol,detJ,wt,Nt,Upp,dt_alpha_1_minus_alpha); 
				}
				
// \/this is for test ////////////////////////////////////////////////////////////				
//        double F_I[9];
//        double Jn = det3x3(F);
//        inverse(F,3,F_I);
//  
//
//        double Fn[9], S[9], C[9], L[81];
//        Fn[0] = Fn[4] = Fn[8] = 1.0;
//        Fn[1] = Fn[2] = Fn[3] = Fn[5] = Fn[6] = Fn[7] = 0.0;  
//
//        devStressFuncPtr Stress;
//        matStiffFuncPtr Stiffness;
//        
//        const int mat = elem[ii].mat[2];
//        Stress = getDevStressFunc(1,&hommat[mat]);
//        Stiffness = getMatStiffFunc(1,&hommat[mat]);
//
//        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//        3,3,3,1.0,F,3,F,3,0.0,C,3);
//		

//        Stress(C,&hommat[mat],S);
//        Stiffness(C,&hommat[mat],L);
//
//        HW_Kuu_at_ip(Kuu_,nne,nne,ST,Fn,F,F_I,1.0,Tn,Jn,Pn,S,L,detJ,wt,0);
//        HW_Kup_at_ip(Kup_,nne,nne,npres,Np,ST,F,1.0,Tn,Jn,detJ,wt,0);   
//        HW_Ktt_at_ip(Ktt_,nVol,Nt,Tn,Jn,kappa,Upp/kappa,detJ,wt);             
//        HW_Kpt_at_ip(Ktp_,npres,Np,nVol,Nt,Tn,Jn,detJ,wt);
//        const damage *ptrDam = &(eps[ii].dam[ip]);
//        material_geometric_sitff_ip(Ku_,nne,ST,F,S,L,ptrDam,detJ,wt);
//        ip++;
//
//        if(i==0)
//        {
//          for(int a=0; a<=nVol*nVol; a++)
//            printf("%e, %e, %e\n", Ktp[a] + alpha*(1.0-alpha)*dt*Ktp_[a], Ktp[a], -alpha*(1.0-alpha)*dt*Ktp_[a]);
//        }
//
// /\this is for test ////////////////////////////////////////////////////////////				
				
       }
    }
  } 
// this is for test ////////////////////////////////////////////////////////////
//  dealoc1(Kuu_);
//  dealoc1(Ku_);
//  dealoc1(Kup_);
//  dealoc1(Ktt_);
//  dealoc1(Ktp_);  
// this is for test ////////////////////////////////////////////////////////////  
          
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

  
  for(int a=0; a<nne*nne*nsd*nsd; a++)
  	Kuu[a] += KsI[a];

  condense_K_out(Ks,nne,nsd,npres,nVol,
                     Kuu,Kut,Kup,Ktu,Ktt,Ktp,Kpu,Kpt,Kpp);


  free(Kuu);
  free(Kut);
  free(Kup);
  
  free(Ktu);
  free(Ktt);
  free(Ktp);
  
  free(Kpu);
  free(Kpt);  
  free(Kpp);        
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
    	u[a*nsd+b] = r_e[a*ndofn+b]; 		  			
    
    if(npres==nne)
      P[a] = r_e[a*ndofn+nsd];
  }	        
              
  int ndofe = nne*ndofn;
  memset(Ks,0,ndofe*ndofe*sizeof(double));
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
          
  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UPP = getD2UdJ2Func(1, &hommat[mat]);

  int err = 0;
  
  Matrix(double) Kuu_add;
  Matrix(double) Kuu,Kut,Kup;
  Matrix(double) Ktu,Ktt,Ktp;
  Matrix(double) Kpu,Kpt,Kpp;

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

  DISP_stiffmat_el(Kuu_add.m_pdata,ii,ndofn,nne,x,y,z,elem,
    hommat,nod,node,eps,sig,sup,r_e);                       

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
      Tn += Mat_v(Nt,a+1,1)*eps[ii].T[a*3+0];
    
    for(int a=0; a<npres; a++)
      Pn += Mat_v(Np,a+1,1)*P[a];       
            
    UPP(Tn,&hommat[mat],&Upp);
    Upp = Upp*kappa;
    
    shape_tensor(nne,ndn,fe.temp_v.N_x.m_pdata,fe.temp_v.N_y.m_pdata,fe.temp_v.N_z.m_pdata,ST_tensor);
    shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
                    
    def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u,F_mat);
    mat2array(F,CONST_2(double) F_mat,3,3);        
                       
    add_3F_Kuu_ip(Kuu.m_pdata,nne ,ST,F,fe.detJ,fe.temp_v.w_ip,Pn,Tn,-1.0);    
    add_3F_Kup_ip(Kup.m_pdata,nne ,npres, ST,F,fe.detJ,fe.temp_v.w_ip,Np.m_pdata,-1.0);
    add_3F_Kpu_ip(Kpu.m_pdata,nne ,npres, ST,F,fe.detJ,fe.temp_v.w_ip,Np.m_pdata,-1.0);        
    add_3F_Kut_ip(Kut.m_pdata,nne ,nVol,  ST,F,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,-1.0);
    add_3F_Ktu_ip(Ktu.m_pdata,nne ,nVol,  ST,F,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,-1.0);        
    add_3F_Ktp_ip(Ktp.m_pdata,nVol,npres,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Np.m_pdata,-1.0);
    add_3F_Kpt_ip(Kpt.m_pdata,nVol,npres,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Np.m_pdata,-1.0);				
    add_3F_Ktt_ip(Ktt.m_pdata,nVol,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Upp,-1.0);         				
  } 
  
  for(int a=0; a<nne; a++)
  {
    for(int b=1; b<=nsd; b++)
    {
      for(int c=0; c<nne; c++)
      {
        for(int d=1; d<=nsd; d++)
          Mat_v(Kuu, a*nsd+b, c*nsd+d) += Mat_v(Kuu_add, a*ndofn+b, c*ndofn+d);         
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
  
  free(P); free(u);    
}



void add_3F_el(double *Ks,
        const int ii,
        double *KsI,
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
        const NODE *node,
        const double *u,
        const double *P,
        double dt,
        SIG *sig,
        EPS *eps,
        double alpha)
{ 
  int ndofe = nne*ndofn;
  memset(Ks,0,ndofe*ndofe*sizeof(double));
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
  UPP = getD2UdJ2Func(1, &hommat[mat]);

  int err = 0;  
  

  double *Kuu = aloc1(nne*nsd*nne*nsd);
  double *Kut = aloc1(nne*nsd*nVol);
  double *Kup = aloc1(nne*nsd*npres);
  
  double *Ktu = aloc1(nVol*nne*nsd);
  double *Ktt = aloc1(nVol*nVol);
  double *Ktp = aloc1(nVol*npres);
  
  double *Kpu = aloc1(npres*nsd*nne);
  double *Kpt = aloc1(npres*nVol);
  double *Kpp = aloc1(npres*npres);
  
	memset(Kuu,0,nne*nsd*nne*nsd*sizeof(double));  	
	memset(Kut,0,   nne*nsd*nVol*sizeof(double));
	memset(Kup,0,  nne*nsd*npres*sizeof(double));	
	memset(Ktu,0,   nne*nsd*nVol*sizeof(double));			
	memset(Ktt,0,      nVol*nVol*sizeof(double));		
	memset(Ktp,0,     nVol*npres*sizeof(double));		
	memset(Kpu,0,  nne*nsd*npres*sizeof(double));			
	memset(Kpt,0,     nVol*npres*sizeof(double));  
	memset(Kpp,0,    npres*npres*sizeof(double));

  // make sure the stiffenss matrix contains all zeros

  
  // INTEGRATION
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
  
  // allocate space for the shape functions, derivatives etc 
  double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST, J;
  double Upp;
 
	Na = aloc1(nne); 
  Np = aloc1(npres);
  Nt = aloc1(nVol);
  N_x = aloc1(nne);
  N_y = aloc1(nne);
  N_z = aloc1(nne);
  ST_tensor = aloc4(3,3,nsd,nne);
  ST = aloc1(3*3*nsd*nne);
  
  //=== INTEGRATION LOOP ===
  integrate(itg_order,&npt_x,&npt_y,&npt_z,
          int_pt_ksi,int_pt_eta,int_pt_zet,
          weights);
  
  double **F_mat, *F;
  
  F = aloc1(9);
  F_mat = aloc2(3,3);
  
// \/this is for test ////////////////////////////////////////////////////////////			  
//  double *Kuu_ = aloc1(nne*nsd*nne*nsd);
//  double *Ku_ = aloc1(nne*nsd*nne*nsd);
//  memset(Kuu_,0,ndofe*ndofe*sizeof(double));        
//  memset(Ku_,0,ndofe*ndofe*sizeof(double));  
//  
//  double *Kup_ = aloc1(nne*nsd*npres);
//  memset(Kup_,0,nne*nsd*npres*sizeof(double));   
//  
//  double *Ktt_ = aloc1(nVol*nVol);
//  memset(Ktt_,0,nVol*nVol*sizeof(double));  
//  
//  double *Ktp_ = aloc1(nVol*npres);
//  memset(Ktp_,0,nVol*npres*sizeof(double));    
// /\this is for test ////////////////////////////////////////////////////////////			
  
  int ip = 0;
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
       
        if(npres==1)
          Nt[0] = 1.0;
        
        if(nVol==1)
          Np[0] = 1.0;

        for(int a=0; a<nVol; a++)
          Tn += Nt[a]*(alpha_1*eps[ii].T[a*3+1] + alpha_2*eps[ii].T[a*3+0]);
           
        for(int a=0; a<npres; a++)
          Pn += Np[a]*P[a];
          
        UPP(Tn,&hommat[mat],&Upp);
        Upp = Upp*kappa;
        shape_tensor(nne,ndn,N_x,N_y,N_z,ST_tensor);
        shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
        
        def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u,F_mat);
        mat2array(F,CONST_2(double) F_mat,3,3);        
       
        add_3F_Kuu_ip(Kuu,nne,ST,F,detJ,wt,Pn,Tn,dt_alpha_1_minus_alpha);
    
        if(include_inertia==0 || (alpha>1.0e-15 && fabs(alpha-1.0)>1.0e-15))
        {       
          add_3F_Kup_ip(Kup,nne,npres, ST,F,detJ,wt,Np,dt_alpha_1_minus_alpha);
          add_3F_Kpu_ip(Kpu,nne,npres, ST,F,detJ,wt,Np,dt_alpha_1_minus_alpha);        
          add_3F_Kut_ip(Kut,nne,nVol,  ST,F,detJ,wt,Nt,dt_alpha_1_minus_alpha);
          add_3F_Ktu_ip(Ktu,nne,nVol,  ST,F,detJ,wt,Nt,dt_alpha_1_minus_alpha);        
				  add_3F_Ktp_ip(Ktp,nVol,npres,detJ,wt,Nt,Np,dt_alpha_1_minus_alpha);
				  add_3F_Kpt_ip(Kpt,nVol,npres,detJ,wt,Nt,Np,dt_alpha_1_minus_alpha);				
				  add_3F_Ktt_ip(Ktt,nVol,detJ,wt,Nt,Upp,dt_alpha_1_minus_alpha); 
				}
				
// \/this is for test ////////////////////////////////////////////////////////////				
//        double F_I[9];
//        double Jn = det3x3(F);
//        inverse(F,3,F_I);
//  
//
//        double Fn[9], S[9], C[9], L[81];
//        Fn[0] = Fn[4] = Fn[8] = 1.0;
//        Fn[1] = Fn[2] = Fn[3] = Fn[5] = Fn[6] = Fn[7] = 0.0;  
//
//        devStressFuncPtr Stress;
//        matStiffFuncPtr Stiffness;
//        
//        const int mat = elem[ii].mat[2];
//        Stress = getDevStressFunc(1,&hommat[mat]);
//        Stiffness = getMatStiffFunc(1,&hommat[mat]);
//
//        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//        3,3,3,1.0,F,3,F,3,0.0,C,3);
//		

//        Stress(C,&hommat[mat],S);
//        Stiffness(C,&hommat[mat],L);
//
//        HW_Kuu_at_ip(Kuu_,nne,nne,ST,Fn,F,F_I,1.0,Tn,Jn,Pn,S,L,detJ,wt,0);
//        HW_Kup_at_ip(Kup_,nne,nne,npres,Np,ST,F,1.0,Tn,Jn,detJ,wt,0);   
//        HW_Ktt_at_ip(Ktt_,nVol,Nt,Tn,Jn,kappa,Upp/kappa,detJ,wt);             
//        HW_Kpt_at_ip(Ktp_,npres,Np,nVol,Nt,Tn,Jn,detJ,wt);
//        const damage *ptrDam = &(eps[ii].dam[ip]);
//        material_geometric_sitff_ip(Ku_,nne,ST,F,S,L,ptrDam,detJ,wt);
//        ip++;
//
//        if(i==0)
//        {
//          for(int a=0; a<=nVol*nVol; a++)
//            printf("%e, %e, %e\n", Ktp[a] + alpha*(1.0-alpha)*dt*Ktp_[a], Ktp[a], -alpha*(1.0-alpha)*dt*Ktp_[a]);
//        }
//
// /\this is for test ////////////////////////////////////////////////////////////				
				
       }
    }
  } 
// this is for test ////////////////////////////////////////////////////////////
//  dealoc1(Kuu_);
//  dealoc1(Ku_);
//  dealoc1(Kup_);
//  dealoc1(Ktt_);
//  dealoc1(Ktp_);  
// this is for test ////////////////////////////////////////////////////////////  
          
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

  
  for(int a=0; a<nne*nne*nsd*nsd; a++)
  	Kuu[a] += KsI[a];

  condense_K_out(Ks,nne,nsd,npres,nVol,
                     Kuu,Kut,Kup,Ktu,Ktt,Ktp,Kpu,Kpt,Kpp);


  free(Kuu);
  free(Kut);
  free(Kup);
  
  free(Ktu);
  free(Ktt);
  free(Ktp);
  
  free(Kpu);
  free(Kpt);  
  free(Kpp);        
}

/* This function may not be used outside this file */
static int el_stiffmat(int i, /* Element ID */
			double **Lk,
			int *Ap,
			int *Ai,
			long ndofn,
			ELEMENT *elem,
			NODE *node,
			HOMMAT *hommat,
			MATGEOM matgeom,
			SIG *sig,
			EPS *eps,
			double *d_r,
			double *r,
			long npres,
			SUPP sup,
			long iter,
			double nor_min,
			double dt,
			CRPL *crpl,
			double stab,
			long FNR,
			double lm,
			double *f_u,
			int myrank,
			int nproc,
			long GDof,
			COMMUN comm,
			int *Ddof,
			int interior,
			const int analysis,
			PGFEM_HYPRE_solve_info *PGFEM_hypre,
			double alpha, double *r_n, double *r_n_1)
{
/* make decision to include ineria*/
  
  const int mat = elem[i].mat[2];
  double rho = hommat[mat].density;
  long include_inertia = 1;
  int nsd = 3;
  
  if(fabs(rho)<1.0e-15)
  {
    include_inertia = 0;
  }
/* decision end*/ 
  
  
  int err = 0;
  long j,l,nne,ndofe,*cnL,*cnG,*nod,II;
  double *lk,*x,*y,*z,*r_e,*sup_def,*fe;
  long kk;
  
  /* Number of element nodes */
  nne = elem[i].toe;
    
  /* Nodes on element */
  nod = aloc1l (nne);
  elemnodes (i,nne,nod,elem);
    
  /* Element Dof */
  ndofe = get_ndof_on_elem_nodes(nne,nod,node);

  /* allocation */
  cnL = aloc1l (ndofe);
  cnG = aloc1l (ndofe);
  lk = aloc1 ((ndofe*ndofe)); 

  const int nne_t = nne + elem[i].n_bub;

  if (analysis == MINI
      || analysis == MINI_3F){ /* P1+B/P1 element */
    x = aloc1 (nne_t);
    y = aloc1 (nne_t);
    z = aloc1 (nne_t);
  } else {
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);
  }
  r_e = aloc1 (ndofe);
  fe = aloc1 (ndofe);
  if(sup->npd>0){
    sup_def = aloc1(sup->npd);
  } else {
    sup_def = NULL;
  }
    
  /* Coordinates of nodes */
  /*
   * The coordinates define the configuration for computing gradients. 
   */
  if(sup->multi_scale){
    /* multi-scale analysis is TOTAL LAGRANGIAN for all element
       formulations */
    nodecoord_total (nne,nod,node,x,y,z);
  } else {
    switch(analysis){
    case DISP:
      nodecoord_total (nne,nod,node,x,y,z);
      break;
    default:
      nodecoord_updated (nne,nod,node,x,y,z);
      break;
    }
  }

  /* if P1+B/P1, get element centroid coords */
  if (analysis == MINI || analysis == MINI_3F){
    element_center(nne,x,y,z);
  }
    
  /* code numbers on element */
  get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cnL);
  get_dof_ids_on_elem_nodes(1,nne,ndofn,nod,node,cnG);
    
  /*=== deformation on element ===*/

  /* set the increment of applied def=0 on first iter */
  if (iter == 0) {
    for (j=0;j<sup->npd;j++){
      sup_def[j] = sup->defl_d[j];
      sup->defl_d[j] = 0.0;
    }
  }

  /* get the deformation on the element */
  if(sup->multi_scale){
    /* multi-scale analysis is TOTAL LAGRANGIAN for all element
       formulations */
    def_elem_total(cnL,ndofe,r,d_r,elem,node,sup,r_e);
  } else {
    switch(analysis){
    case DISP: /* TOTAL LAGRANGIAN */
      def_elem_total(cnL,ndofe,r,d_r,elem,node,sup,r_e);
      break;
    default:
      def_elem (cnL,ndofe,d_r,elem,node,r_e,sup,0);
      break;
    }
  }

  /* recover thei increment of applied def on first iter */
  if (iter == 0){
    for (j=0;j<sup->npd;j++)
      sup->defl_d[j] = sup_def[j];
   }
    
  nulld (lk,ndofe*ndofe);  
  if(include_inertia)
  {
    switch(analysis){
      case DISP:
      {
        /* Get TOTAL deformation on element; r_e already contains
         * INCREMENT of deformation, add the deformation from previous. */
///////////////////////////////////////////////////////////////////////////////////        
        double *r_en, *r_mid, *r0;
        double *lk_k, *lk_i;
        lk_k = aloc1(ndofe*ndofe);
        lk_i = aloc1(ndofe*ndofe);
        
        r_en  = aloc1(ndofe);
        r_mid = aloc1(ndofe);
        r0    = aloc1(ndofe);
              
        for (long I=0;I<nne;I++)
        {
          for(long J=0; J<ndofn; J++)
            r0[I*ndofn + J] = r_n[nod[I]*ndofn + J];
        }
        
        def_elem (cnL,ndofe,r,elem,node,r_en,sup,1);
        vvplus(r_e,r_en,ndofe);
        
        mid_point_rule(r_mid, r0, r_e, alpha, ndofe);
        
        err = DISP_stiffmat_el(lk_k,i,ndofn,nne,x,y,z,elem,
                hommat,nod,node,eps,sig,sup,r_mid);
        
        stiffmat_w_inertia_el(lk_i,i,ndofn,nne,x,y,z,elem,hommat,node,dt);
        
        memset(lk,0,ndofe*ndofe*sizeof(double));
        for(long a = 0; a<ndofe*ndofe; a++)
            lk[a] = -lk_i[a]-alpha*(1-alpha)*dt*lk_k[a];

        free(r_en);
        free(r_mid);
        free(r0);
        free(lk_k);
        free(lk_i);
        break;
      }
      case TF: 
      {
        /* Get TOTAL deformation on element; r_e already contains
         * INCREMENT of deformation, add the deformation from previous. */
///////////////////////////////////////////////////////////////////////////////////
        double *r_en, *r_mid, *r0;
        double *lk_k, *lk_i;
        double *U, *P, *T;
        int nVol = N_VOL_TF;
          
        lk_k = aloc1(ndofe*ndofe);
        lk_i = aloc1(nne*nne*nsd*nsd);
        
        r_en  = aloc1(ndofe);
        r_mid = aloc1(ndofe);
        r0    = aloc1(ndofe);
        U     = aloc1(nne*nsd);
        P     = aloc1(npres);
        T     = aloc1(nVol);        

        for (long I=0;I<nne;I++)
        {
          for(long J=0; J<ndofn; J++)
            r0[I*ndofn + J] = r_n[nod[I]*ndofn + J];
        }
        def_elem(cnL,ndofe,r,elem,node,r_en,sup,1);
        vvplus(r_e,r_en,ndofe);

        mid_point_rule(r_mid, r0, r_e, alpha, ndofe);

        for (long I=0;I<nne;I++)
        {
          for(long J=0; J<nsd; J++)
            U[I*nsd + J] = r_mid[I*ndofn + J];
          if(npres==nne)  
            P[I] = r_mid[I*ndofn + nsd];
        }
        if(npres==1)          
            P[0] =  (1.0-alpha)*eps[i].d_T[1] + alpha*eps[i].d_T[0]; 

        err = DISP_stiffmat_el(lk_k,i,ndofn,nne,x,y,z,elem,
                hommat,nod,node,eps,sig,sup,r_mid);

        stiffmat_w_inertia_el(lk_i,i,nsd,nne,x,y,z,elem,hommat,node,dt);

        for(int a=0; a<nne*nsd*nne*nsd; a++)
        	lk_i[a] = -lk_i[a] - alpha*(1-alpha)*dt*lk_k[a];
        	
        add_3F_el(lk,i,lk_i, ndofn,nne,npres,nVol,nsd,
                  x,y,z,elem,hommat,node,U,P,dt,sig,eps,alpha);   
                                           
///////////////////////////////////////////////////////////////////////////////////          
        free(r_en);
        free(r_mid);
        free(r0);
        free(U);
        free(P);
        free(T);        
        free(lk_k);
        free(lk_i);

        break;
      }        
    } /* switch (analysis) */
  }
  else
  {        
  switch(analysis){
  case STABILIZED:
    err = stiffmatel_st (i,ndofn,nne,x,y,z,elem,hommat,nod,node,sig,eps,
			 sup,r_e,npres,nor_min,lk,dt,stab,FNR,lm,fe);
    break;
  case MINI:
    err = MINI_stiffmat_el(lk,i,ndofn,nne,x,y,z,elem,
			   hommat,nod,node,eps,sig,r_e);
    break;
  case MINI_3F:
    err = MINI_3f_stiffmat_el(lk,i,ndofn,nne,x,y,z,elem,
			      hommat,nod,node,eps,sig,r_e);
    break;
  case DISP:
    err = DISP_stiffmat_el(lk,i,ndofn,nne,x,y,z,elem,
			   hommat,nod,node,eps,sig,sup,r_e);
    break;
  case TF:
    {
      
        int nVol = N_VOL_TF;
        stiffmat_3f_el(lk,i,ndofn,nne,npres,nVol,nsd,
                  x,y,z,elem,hommat,nod,node,dt,sig,eps,sup,r_e);
        break;
      }        
  default:
    err = stiffmatel_fd (i,ndofn,nne,nod,x,y,z,elem,matgeom,
			 hommat,node,sig,eps,r_e,npres,
			 nor_min,lk,dt,crpl,FNR,lm,fe,analysis);
    break;
  } /* switch (analysis) */
  } /* if(include_inertia) */
    
  if (PFEM_DEBUG){
    char filename[50];
    switch(analysis){
    case STABILIZED:
      sprintf(filename,"stab_stiff_%d.log",myrank);
      break;
    case MINI:
      sprintf(filename,"MINI_stiff_%d.log",myrank);
      break;
    case MINI_3F:
      sprintf(filename,"MINI_3f_stiff_%d.log",myrank);
      break;
    default:
      sprintf(filename,"stiff_%d.log",myrank);
      break;
    }

    FILE *output;
    output = fopen(filename,"a");
    print_array_d(output,lk,ndofe*ndofe,ndofe,ndofe);
    fclose(output);
  }

  /* Localization of TANGENTIAL LOAD VECTOR */
  if (periodic == 1 && (FNR == 2 || FNR == 3)){
    for (l=0;l<nne;l++){
      for (kk=0;kk<node[nod[l]].ndofn;kk++){
	II = node[nod[l]].id[kk]-1;
	if (II < 0)  continue;
	f_u[II] += fe[l*node[nod[l]].ndofn+kk];
      }/*end l */
    }/*end kk */
  }/* end periodic */

  /* Assembly */
  PLoc_Sparse (Lk,lk,Ai,Ap,cnL,cnG,ndofe,Ddof,GDof,
	       myrank,nproc,comm,interior,PGFEM_hypre,analysis);

  /*  dealocation  */
  free (cnL);
  free (cnG);
  free (nod);
  free (lk);
  free (x);
  free (y);
  free (z);
  free (fe);
  free (r_e);
  free (sup_def);

  return err;

} /* ELEMENT STIFFNESS */

/* This function may not be used outside of this file */
static void coel_stiffmat(int i, /* coel ID */
			  double **Lk,
			  int *Ap,
			  int *Ai,
			  long ndofc,
			  ELEMENT *elem,
			  NODE *node,
			  EPS *eps,
			  double *d_r,
			  double *r,
			  long npres,
			  SUPP sup,
			  long iter,
			  double nor_min,
			  double dt,
			  CRPL *crpl,
			  double stab,
			  COEL *coel,
			  long FNR,
			  double lm,
			  double *f_u,
			  int myrank,
			  int nproc,
			  long *DomDof,
			  long GDof,
			  COMMUN comm,
			  int *Ddof,
			  int interior,
			  const int analysis,
			  PGFEM_HYPRE_solve_info *PGFEM_hypre)
{
  long j,l,nne,ndofe,*cnL,*cnG,*nod,P,R,II;
  double *lk,*x,*y,*z,*r_e,*sup_def,*fe, *X, *Y;
  long kk;

  nne = coel[i].toe/2;
  ndofe = coel[i].toe*ndofc;
      
  /* Alocation */
  nod = aloc1l (coel[i].toe);
  lk = aloc1 (ndofe*ndofe);
  x = aloc1 (coel[i].toe); 
  y = aloc1 (coel[i].toe);
  z = aloc1 (coel[i].toe);
  r_e = aloc1 (ndofe);
  fe = aloc1 (ndofe);
  cnL = aloc1l (ndofe);
  cnG = aloc1l (ndofe);

  if(sup->npd > 0){
    sup_def = aloc1(sup->npd);
  } else {
    sup_def = NULL;
  }

  /* Element Node */
  for (j=0;j<coel[i].toe;j++)
    nod[j] = coel[i].nod[j];

  /* Coordinates */
  /* 
   * NOTE: get updated coordinates for all formulations, computes jump
   * based on both coordinates and deformation.
   */
  nodecoord_updated (coel[i].toe,nod,node,x,y,z);
      
  /* code numbers on element */
  get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nod,node,cnL);
  get_dof_ids_on_elem_nodes(1,coel[i].toe,ndofc,nod,node,cnG);

  /* deformation on element */
  if (iter == 0){
    for (j=0;j<sup->npd;j++){
      sup_def[j] = sup->defl_d[j];
      sup->defl_d[j] = 0.0;
    }
  }

  def_elem (cnL,ndofe,d_r,elem,node,r_e,sup,0);
  if (iter == 0)
    for (j=0;j<sup->npd;j++)
      sup->defl_d[j] = sup_def[j];
      
  if (periodic == 1){/* Periodic */

    X = aloc1 (3);
    Y = aloc1 (3);
	
    for (j=0;j<coel[i].toe;j++){
      X[0] = x[j];
      X[1] = y[j];
      X[2] = z[j];
      for (P=0;P<3;P++){
	Y[P] = 0.0;
	for (R=0;R<3;R++){
	  Y[P] += eps[0].F[P][R]*X[R];
	}
      }
      for (P=0;P<3;P++){
	II = node[nod[j]].id[P];
	    
	if (P == 0)
	  x[j] = Y[P];
	if (P == 1)
	  y[j] = Y[P];
	if (P == 2)
	  z[j] = Y[P];
	    
	if (II > 0){
	  if (P == 0)
	    x[j] += r[II-1];
	  if (P == 1)
	    y[j] += r[II-1];
	  if (P == 2)
	    z[j] += r[II-1];
	}
      }/* P < 3 */
    }/* j < coel[i].toe */

    dealoc1 (X);
    dealoc1 (Y);
	
  }/* end periodic */

  /**************************/
  /*     FOR PERIODIC       */
  /*                        */
  /* xn+1 = Fn+1*X + un + u */
  /*    x = Fn+1*X + un     */
  /**************************/

  nulld (lk,ndofe*ndofe); 
  stiff_mat_coh (i,ndofc,nne,nod,x,y,z,coel,r_e,lk,
		 nor_min,eps,FNR,lm,fe,myrank);
      
  /* Assembly */
  PLoc_Sparse (Lk,lk,Ai,Ap,cnL,cnG,ndofe,Ddof,
	       GDof,myrank,nproc,comm,interior,PGFEM_hypre,analysis);

  /* Localization of TANGENTIAL LOAD VECTOR */
  if (periodic == 1 && (FNR == 2 || FNR == 3)){
    for (l=0;l<coel[i].toe;l++){
      for (kk=0;kk<ndofc;kk++){
	II = node[nod[l]].id[kk]-1;
	if (II < 0)  continue;
	f_u[II] += fe[l*ndofc+kk];
      }/*end l */
    }/*end kk */
  }/* end periodic */

  /*  dealocation  */
  dealoc1l (cnL);
  dealoc1l (cnG);
  dealoc1l (nod);
  dealoc1 (lk); 
  dealoc1 (x);
  dealoc1(y); 
  dealoc1 (z);
  dealoc1 (r_e);
  dealoc1 (fe);
  free(sup_def);

} /* COHESIVE ELEMENT STIFFNESS */

static int bnd_el_stiffmat(int belem_id,
			   double **Lk,
			   int *Ap,
			   int *Ai,
			   long ndofn,
			   ELEMENT *elem,
			   BOUNDING_ELEMENT *b_elems,
			   NODE *node,
			   HOMMAT *hommat,
			   MATGEOM matgeom,
			   SIG *sig,
			   EPS *eps,
			   double *d_r,
			   double *r,
			   long npres,
			   SUPP sup,
			   long iter,
			   double nor_min,
			   double dt,
			   CRPL *crpl,
			   double stab,
			   long FNR,
			   double lm,
			   double *f_u,
			   int myrank,
			   int nproc,
			   long GDof,
			   COMMUN comm,
			   int *Ddof,
			   int interior,
			   const int analysis,
			   PGFEM_HYPRE_solve_info *PGFEM_hypre)
{
  int err = 0;
  const BOUNDING_ELEMENT *ptr_be = &b_elems[belem_id];
  const ELEMENT *ptr_ve = &elem[ptr_be->vol_elem_id];
  const long *ptr_vnodes = ptr_ve->nod;
  const int nn_ve = ptr_ve->toe;

  /* get coordinated for BOUNDING ELEMENT */
  double *x = aloc1(nn_ve);
  double *y = aloc1(nn_ve);
  double *z = aloc1(nn_ve);
  switch(analysis){
  case DISP:
    nodecoord_total(nn_ve,ptr_vnodes,node,x,y,z);
    break;
  default:
    nodecoord_updated(nn_ve,ptr_vnodes,node,x,y,z);
    break;
  }

  /* get the local and global dof id's */
  int ndof_ve = get_ndof_on_bnd_elem(node,ptr_be,elem);

  long *cn_ve = aloc1l(ndof_ve);
  long *Gcn_ve = aloc1l(ndof_ve);

  get_dof_ids_on_bnd_elem(0,ndofn,node,ptr_be,elem,cn_ve);
  get_dof_ids_on_bnd_elem(1,ndofn,node,ptr_be,elem,Gcn_ve);

  /* compute the deformation on the element */
  double *v_disp = aloc1(ndof_ve);

  if(iter == 0){
    /* on iter == 0, null increment of deflection */
    double *sup_def = aloc1(sup->npd);
    for (int j=0;j<sup->npd;j++){
      sup_def[j] = sup->defl_d[j];
      sup->defl_d[j] = 0.0;
    }

    def_elem(cn_ve,ndof_ve,d_r,NULL,NULL,v_disp,sup,0);

    for (int j=0;j<sup->npd;j++){
      sup->defl_d[j] = sup_def[j];
    }
    free(sup_def);
  } else {
    def_elem(cn_ve,ndof_ve,d_r,NULL,NULL,v_disp,sup,0);
  }

  if(analysis == DISP){ /* TOTAL LAGRANGIAN formulation */
    double *ve_n = aloc1(ndof_ve);
    def_elem(cn_ve,ndof_ve,r,NULL,NULL,ve_n,sup,1);
    vvplus(v_disp,ve_n,ndof_ve);
    free(ve_n);
  }

  /* compute the local stiffness matrix */
  double *lk = aloc1(ndof_ve*ndof_ve);
  if(analysis == DISP){
    err += DISP_stiffmat_bnd_el(lk,belem_id,ndofn,ndof_ve,
				x,y,z,b_elems,elem,hommat,node,eps,
				sig,sup,v_disp);
  } else {
    /* Not implemented, do nothing */
  }

  /* only assemble to global stiffness if no error */
  if(err == 0){
    /* PLoc_Sparse_rec(Lk,lk,Ai,Ap,Gcn_be,Gcn_ve,ndof_be,ndof_ve,Ddof, */
    /* 		   GDof,myrank,nproc,comm,interior); */
    PLoc_Sparse(Lk,lk,Ai,Ap,cn_ve,Gcn_ve,ndof_ve,Ddof,
		GDof,myrank,nproc,comm,interior,PGFEM_hypre,analysis);
  }


  free(x);
  free(y);
  free(z);

  free(cn_ve);
  free(Gcn_ve);

  free(v_disp);
  free(lk);

  return err;
} /* Bounding element stiffnes matrix */

/* This is the re-written function which computes elem stiffness on
   boundaries first, then interior elem stiffnesses before
   assembly. */
int stiffmat_fd(int *Ap,
		 int *Ai,
		 long ne,
		 int n_be,
		 long ndofn,
		 ELEMENT *elem,
		 BOUNDING_ELEMENT *b_elems,
		 long nbndel,
		 long *bndel,
		 NODE *node,
		 HOMMAT *hommat,
		 MATGEOM matgeom,
		 SIG *sig,
		 EPS *eps,
		 double *d_r,
		 double *r,
		 long npres,
		 SUPP sup,
		 long iter,
		 double nor_min,
		 double dt,
		 CRPL *crpl,
		 double stab,
		 long nce,
		 COEL *coel,
		 long FNR,
		 double lm,
		 double *f_u,
		 int myrank,
		 int nproc,
		 long *DomDof,
		 long GDof,
		 COMMUN comm,
		 MPI_Comm mpi_comm,
		 PGFEM_HYPRE_solve_info *PGFEM_hypre,
		 const PGFem3D_opt *opts,double alpha, double *r_n, double *r_n_1)
{
  int err = 0;
  long i,ndofc;
  int *Ddof;
  double **Lk,**recieve;
  MPI_Status *sta_s,*sta_r;
  MPI_Request *req_s,*req_r;

  /* interior element counters */
  int idx = 0;
  int skip = 0;

  if(opts->solverpackage == 0){
    PGFEM_printf("BlockSolve no longer supported\n");
    abort();
  }

  err += init_and_post_stiffmat_comm(&Lk,&recieve,&req_r,&sta_r,
				     mpi_comm,comm);

  /* Lk = (double**) PGFEM_calloc (nproc,sizeof(double*)); */
  /* for (i=0;i<nproc;i++) { */
  /*   if (myrank == i || comm->S[i] == 0) */
  /*     k = 1; */
  /*   else  */
  /*     k = comm->AS[i]; */
  /*   Lk[i] = (double*) PGFEM_calloc (k,sizeof(double)); */
  /* } */
  /* if (Lk == NULL){ */
  /*   PGFEM_printf ("\n Memory is full. %s:%s:%d\n",__func__,__FILE__,__LINE__); */
  /*   fflush(stdout);  */
  /*   PGFEM_Comm_code_abort (mpi_comm,i); */
  /* } */
  
  /* /\* Allocate recieve *\/ */
  /* recieve = (double**) PGFEM_calloc (nproc,sizeof(double*)); */
  /* for (i=0;i<nproc;i++) { */
  /*   if (comm->AR[i] == 0) */
  /*     KK = 1; */
  /*   else */
  /*     KK = comm->AR[i]; */
  /*   recieve[i] = (double*) PGFEM_calloc (KK,sizeof(double)); */
  /* } */
  
  /* /\* Allocate request fields *\/ */
  /* if (comm->Nr == 0) */
  /*   KK = 1; */
  /* else */
  /*   KK = comm->Nr; */
  /* sta_r = (MPI_Status*) PGFEM_calloc (KK,sizeof(MPI_Status)); */
  /* req_r = (MPI_Request*) PGFEM_calloc (KK,sizeof(MPI_Request)); */
  
  /* /\* Receive data *\/ */
  /* for (i=0;i<comm->Nr;i++){ */
  /*   KK = comm->Nrr[i]; */
  /*   MPI_Irecv (recieve[KK],comm->AR[KK],MPI_DOUBLE,KK, */
  /* 	       MPI_ANY_TAG,mpi_comm,&req_r[i]); */
  /* }/\* end i < nproc *\/ */
  
  /* Allocate */
  Ddof = aloc1i (nproc);
  
  /* Set Ddof */
  Ddof[0] = DomDof[0];
  for (i=1;i<nproc;i++)
    Ddof[i] = Ddof[i-1] + DomDof[i];
  
  /***** COMM BOUNDARY ELEMENTS *****/
  for(i=0; i<nbndel; i++){
    err += el_stiffmat(bndel[i],Lk,Ap,Ai,ndofn,elem,node,hommat,
		       matgeom,sig,eps,d_r,r,npres,sup,iter,nor_min,
		       dt,crpl,stab,FNR,lm,f_u,myrank,nproc,GDof,comm,
		       Ddof,0,opts->analysis_type,PGFEM_hypre,alpha,r_n,r_n_1);

    /* If there is an error, complete communication and exit */
    if(err != 0) goto send;
  }

  /***** COHESIVE ELEMENTS *****/

  /* Need to split into boundary and interior parts as with the
     regular elements */
  if (opts->cohesive == 1){/* COHESIVE STIFFNESS */
    ndofc = 3; 
    
    /* WHY IS THIS HERE */
    if (iter == 0) FNR = 0;
    
    if (nor_min < 1.e-10)
      nor_min = 1.e-10;
    
    for (i=0;i<nce;i++){
      coel_stiffmat(i,Lk,Ap,Ai,ndofc,elem,node,eps,
		    d_r,r,npres,sup,iter,nor_min,dt,crpl,
		    stab,coel,FNR,lm,f_u,myrank,nproc,DomDof,
		    GDof,comm,Ddof,0,opts->analysis_type,PGFEM_hypre);
    }
  }


  /**** BOUNDING ELEMENTS ******/
  /* In the future, these elements will be listed as with the
     volumetric elements to properly overlay computation and
     communication. For now, this is the best place for them as they
     are a proportinally smaller group than the interior volume
     elements for typical problems and are guarenteed to be on the
     communication boundary for periodic domains. */

  /* temporary for compile testing */
  /* int n_be = 0; */
  /* BOUNDING_ELEMENT *b_elems = NULL; */
  for(i=0; i<n_be; i++){
    err += bnd_el_stiffmat(i,Lk,Ap,Ai,ndofn,elem,b_elems,node,hommat,
			   matgeom,sig,eps,d_r,r,npres,sup,iter,nor_min,
			   dt,crpl,stab,FNR,lm,f_u,myrank,nproc,GDof,
			   comm,Ddof,0,opts->analysis_type,PGFEM_hypre);

    /* If there is an error, complete communication and exit */
    if(err != 0) goto send;
  }

  if (PFEM_DEBUG){
    char ofile[50];
    switch(opts->analysis_type){
    case STABILIZED:
      sprintf(ofile,"stab_el_stiff_send_%d.log",myrank);
      break;
    case MINI:
      sprintf(ofile,"MINI_el_stiff_send_%d.log",myrank);
      break;
    case MINI_3F:
      sprintf(ofile,"MINI_3f_el_stiff_send_%d.log",myrank);
      break;
    default:
      sprintf(ofile,"el_stiff_send_%d.log",myrank);
      break;
    }

    FILE *out;
    out = fopen(ofile,"a");
    for(int send_proc=0; send_proc<nproc; send_proc++){
      if(send_proc==myrank || comm->S[send_proc] ==0) continue;
      for(int n_data=0; n_data<comm->AS[send_proc]; n_data++){
	PGFEM_fprintf(out,"%12.12e ",Lk[send_proc][n_data]);
      }
      PGFEM_fprintf(out,"\n");
    }
    PGFEM_fprintf(out,"\n");
    fclose(out);
  }

  /**********************************/
  /***** SEND BOUNDARY AND COEL *****/
  /**********************************/
 send:
  err += send_stiffmat_comm(&sta_s,&req_s,Lk,mpi_comm,comm);

  /* If error, complete communication and exit */
  if(err != 0) goto wait;

   /**********************************/
  /*****    CONTINUE WORKING    *****/
  /**********************************/

  /***** INTERIOR ELEMENTS *****/

  if(nbndel > 0){/* this is 99% of the time */
    for(i=0; i<ne; i++){
      if(idx < nbndel-1){
	if(i == 0 && idx == 0 && bndel[idx] == 0){
	  idx++;
	  skip++;
	  continue;
	} else if(i == bndel[idx]){
	  idx++;
	  skip++;
	  continue;
	} else if (idx == 0 && i < bndel[idx]){
	  err = el_stiffmat(i,Lk,Ap,Ai,ndofn,elem,node,hommat,matgeom,
			    sig,eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,
			    stab,FNR,lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			    opts->analysis_type,PGFEM_hypre,alpha,r_n,r_n_1);
	} else if (idx > 0 && bndel[idx-1] < i && i < bndel[idx]){
	  err = el_stiffmat(i,Lk,Ap,Ai,ndofn,elem,node,hommat,matgeom,sig,
			    eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,stab,
			    FNR,lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			    opts->analysis_type,PGFEM_hypre,alpha,r_n,r_n_1);
	} else {
	  PGFEM_printf("[%d]ERROR: problem in determining if element %ld"
		 " is on interior.\n", myrank, i);
	}
      } else {
	if(i != bndel[nbndel-1]){
	  err = el_stiffmat(i,Lk,Ap,Ai,ndofn,elem,node,hommat,matgeom,sig,
			    eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,stab,
			    FNR,lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			    opts->analysis_type,PGFEM_hypre,alpha,r_n,r_n_1);
	}
      }

      /* If there is an error, complete communication and exit */
      if(err != 0) goto wait;
    }

    /* Check to make sure I got all of them */
    if(skip != nbndel - 1){
      PGFEM_printf("[%d]WARNING: number skipped elem != nbndel, check code\n",myrank);
    }
  } else { /* communication by coheisve elements only, nbndel = 0 */
    for(i=0; i<ne; i++){
      err = el_stiffmat(i,Lk,Ap,Ai,ndofn,elem,node,hommat,matgeom,sig,
			eps,d_r,r,npres,sup,iter,nor_min,dt,crpl,stab,FNR,
			lm,f_u,myrank,nproc,GDof,comm,Ddof,1,
			opts->analysis_type,PGFEM_hypre,alpha,r_n,r_n_1);
      /* If there is an error, complete communication and exit */
      if(err != 0) goto wait;
    }
  }

  /**********************************/
  /*****    RECEIVE AND ADD     *****/
  /**********************************/
  
 wait:
  err += assemble_nonlocal_stiffmat(comm,sta_r,req_r,PGFEM_hypre,recieve);
  err += finalize_stiffmat_comm(sta_s,sta_r,req_s,req_r,comm);
    
 /* exit_function: */
  /* Deallocate recieve */
  for (i=0;i<nproc;i++)
    free (recieve[i]);
  free (recieve);
  
  /*  dealocation  */
  for (i=0;i<nproc;i++)
    free(Lk[i]);
  
  free (Lk);
  free (Ddof);
  free (sta_s);
  free (sta_r);
  free (req_s);
  free (req_r);

  return err;
}


/** Assemble non-local parts as they arrive */
int assemble_nonlocal_stiffmat(const COMMUN pgfem_comm,
			       MPI_Status *sta_r,
			       MPI_Request *req_r,
			       PGFEM_HYPRE_solve_info *PGFEM_hypre,
			       double **recv)
{
  int err = 0;
  int comm_idx = 0;
  int n_received = 0;
  while (n_received < pgfem_comm->Nr){
    /* get the communication index */
    err += MPI_Waitany(pgfem_comm->Nr,req_r,&comm_idx,sta_r);

    /* convert communication index to proc_id */
    const int proc = pgfem_comm->Nrr[comm_idx];

    /* get number of rows */
    const int nrows = pgfem_comm->R[proc];

    /* allocate rows and cols to receive */
    int *row_idx = PGFEM_calloc(nrows,sizeof(int));
    int *ncols = PGFEM_calloc(nrows,sizeof(int));
    int *col_idx = PGFEM_calloc(pgfem_comm->AR[proc],sizeof(int));

    /* get row and column ids */
    int idx = 0;
    for(int j=0; j<pgfem_comm->R[proc]; j++){
      row_idx[j] = pgfem_comm->RGID[proc][j];
      ncols[j] = pgfem_comm->RAp[proc][j];
      for(int k=0; k<ncols[j]; k++){
	col_idx[idx] = pgfem_comm->RGRId[proc][idx];
	++idx;
      }
    }

    /* assemble to local part of global stiffness */
    err += HYPRE_IJMatrixAddToValues(PGFEM_hypre->hypre_k,
				     nrows,ncols,row_idx,col_idx,
				     recv[proc]);

    /* free memory */
    free(row_idx);
    free(ncols);
    free(col_idx);

    /* increment counter */
    n_received++;
  }
  return err;
}
