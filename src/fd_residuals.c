#include "fd_residuals.h"
#include "enumerations.h"
#include "utils.h"
#include "allocation.h"
#include "resid_on_elem.h"
#include "stabilized.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "displacement_based_element.h"
#include "matice.h"
#include "elem3d.h"
#include <math.h>
#include <string.h>
#include "three_field_element.h"
#include "condense.h"
#include "new_potentials.h"
#include "tensors.h"
#include "cast_macros.h"
#include "index_macros.h"
#include "def_grad.h"
#include "mkl_cblas.h"
#include "femlib.h"

#define ndn 3
#define N_VOL_TF 1

void residuals_disp_el(double *f,
        const int ii,
        const int ndofn,
        const int nne,
        const int nsd,
        const double *x,
        const double *y,
        const double *z,
        const ELEMENT *elem,
        const HOMMAT *hommat,
        const NODE *node,
        double *r_e)
{
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
      
  Matrix(double) F,S,P,C,C_I;

  int ndofe = nne*ndofn;

  Matrix_construct_init(double,F,3,3,0.0);
  Matrix_construct_init(double,C,3,3,0.0);  
  Matrix_construct_init(double,C_I,3,3,0.0);    
  Matrix_construct_init(double,S,3,3,0.0);
  Matrix_construct_init(double,P,3,3,0.0);

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
                             
  devStressFuncPtr Stress = getDevStressFunc(0,&hommat[mat]);
  dUdJFuncPtr DUDJ = getDUdJFunc(0,hommat);
  
  double AA[9];
                             
  FEMLIB_initialization(&fe, itg_order, 1, nne);
  FEMLIB_set_element(&fe, xe, ii);
  for(int ip = 1; ip<=fe.nint; ip++)
  {   
    FEMLIB_elem_basis_V(&fe, ip);  
         
    double dUdJ = 0.0;

    shape_tensor(nne,ndn,fe.temp_v.N_x.m_pdata,fe.temp_v.N_y.m_pdata,fe.temp_v.N_z.m_pdata,ST_tensor);
    shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
        
    def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,r_e,F_mat);
    mat2array(F.m_pdata,CONST_2(double) F_mat,3,3);
    
    Matrix_AxB(C,1.0,0.0,F,1,F,0);
    double J = 0.0;
    Matrix_det(F, J);
    Matrix_inv(C, C_I);
    		            
    Stress(C.m_pdata,&hommat[mat],S.m_pdata);
    DUDJ(J,&hommat[mat],&dUdJ);
    
    cblas_daxpy(9,kappa*J*dUdJ,C_I.m_pdata,1,S.m_pdata,1);    
            
    Matrix_AxB(P,1.0,0.0,F,0,S,0);
        
    for(int a=0; a<fe.nne; a++)
    {
      for(int b=0; b<fe.nsd; b++)
      {
        for(int c=1; c<=3; c++)
          f[a*ndofn + b] += Mat_v(P,b+1,c)*Mat_v(fe.dN,a+1,c)*fe.detJxW;
      }
    }     
  }
      
  FEMLIB_destruct(&fe);                   
  Matrix_cleanup(F);
  Matrix_cleanup(C);
  Matrix_cleanup(C_I);  
  Matrix_cleanup(S);
  Matrix_cleanup(P);
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
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
      
  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UP = getDUdJFunc(1, &hommat[mat]);  
  UPP = getD2UdJ2Func(1, &hommat[mat]);

  Matrix(double) fu_add;
  Matrix(double) fu,Kut;
  Matrix(double) fp,Kpt;
  Matrix(double) ft,Ktt;
  Matrix(double) Kup,Ktp;
  
  int ndofe = nne*ndofn;

  Matrix_construct_init(double,fu_add, nne*ndofn,1,0.0);
  Matrix_construct_init(double,fu, nne*nsd,1,    0.0);
  Matrix_construct_init(double,fp, npres,  1,    0.0);
  Matrix_construct_init(double,ft, nVol,   1,    0.0);
  Matrix_construct_init(double,Kut,nne*nsd,nVol, 0.0);
  Matrix_construct_init(double,Kpt,npres,  nVol, 0.0);
  Matrix_construct_init(double,Ktt,nVol,   nVol, 0.0);    
  Matrix_construct_init(double,Kup,nne*nsd,npres,0.0);
  Matrix_construct_init(double,Ktp,nVol,   npres,0.0);
          
  DISP_resid_el(fu_add.m_pdata,ii,ndofn,nne,x,y,z,elem,
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
  double Upp, Up;
 
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

    UP(Tn, &hommat[mat], &Up);
    Up = Up*kappa;
            
    shape_tensor(nne,ndn,fe.temp_v.N_x.m_pdata,fe.temp_v.N_y.m_pdata,fe.temp_v.N_z.m_pdata,ST_tensor);
    shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
        
    def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u,F_mat);
    mat2array(F,CONST_2(double) F_mat,3,3);        

		add_3F_Kpt_ip(Kpt.m_pdata,nVol,npres,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Np.m_pdata,-1.0);	
		add_3F_Ktt_ip(Ktt.m_pdata,nVol,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Upp,-1.0);	

    if(npres==nne)
    {
      add_3F_Kut_ip(Kut.m_pdata,nne,nVol,ST,F,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,-1.0);
    }	
    
		if(npres==1)
		{
      add_3F_Kup_ip(Kup.m_pdata,nne,npres,ST,F,fe.detJ,fe.temp_v.w_ip,Np.m_pdata,-1.0);
	  	add_3F_Ktp_ip(Ktp.m_pdata,nVol,npres,fe.detJ,fe.temp_v.w_ip,Nt.m_pdata,Np.m_pdata,-1.0);			
	  }	
		resid_w_inertia_Rp_ip(fp.m_pdata, npres, F, fe.detJ,fe.temp_v.w_ip, Np.m_pdata, Tn);
		resid_w_inertia_Rt_ip(ft.m_pdata, nVol, fe.detJ,fe.temp_v.w_ip, Nt.m_pdata, Pn, Up);
		resid_w_inertia_Ru_ip(fu.m_pdata, nne, ST, F, fe.detJ,fe.temp_v.w_ip, Pn);
  }
  
  for(int a = 0; a<nne; a++)
  {
    for(int b = 0; b<nsd; b++)
      Mat_v(fu, nsd*a+b+1, 1) += Mat_v(fu_add, ndofn*a+b+1, 1); 
  }    

  dealoc4(ST_tensor,3,3,nsd);
  free(ST);
  
  dealoc2(F_mat,3);  
  
  condense_F_out(f,nne,nsd,npres,nVol,fu.m_pdata,ft.m_pdata,fp.m_pdata,
                    Kut.m_pdata,Kup.m_pdata,Ktp.m_pdata,Ktt.m_pdata,Kpt.m_pdata);
                        
  FEMLIB_destruct(&fe);                   
  Matrix_cleanup(xe);
  Matrix_cleanup(Np);  
  Matrix_cleanup(Nt);
  
  Matrix_cleanup(fu_add);
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


void residuals_3f_w_inertia_el(double *f,
        const int ii,
        double *fuI,
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

  Matrix(double) fu,fu1,fu2,Kut;
  Matrix(double) fp,fp1,fp2,Kpt;
  Matrix(double) ft,ft1,ft2,Ktt;
  Matrix(double) Kup,Ktp;
  
  int ndofe = nne*ndofn;

  Matrix_construct_init(double,fu, nne*nsd,1,    0.0);
  Matrix_construct_init(double,fu1,nne*nsd,1,    0.0);
  Matrix_construct_init(double,fu2,nne*nsd,1,    0.0);  
  Matrix_construct_init(double,fp, npres,  1,    0.0);
  Matrix_construct_init(double,fp1,npres,  1,    0.0);
  Matrix_construct_init(double,fp2,npres,  1,    0.0);  
  Matrix_construct_init(double,ft, nVol,   1,    0.0);
  Matrix_construct_init(double,ft1,nVol,   1,    0.0);
  Matrix_construct_init(double,ft2,nVol,   1,    0.0);  
  Matrix_construct_init(double,Kut,nne*nsd,nVol, 0.0);
  Matrix_construct_init(double,Kpt,npres,  nVol, 0.0);
  Matrix_construct_init(double,Ktt,nVol,   nVol, 0.0);    
  Matrix_construct_init(double,Kup,nne*nsd,npres,0.0);
  Matrix_construct_init(double,Ktp,nVol,   npres,0.0);
      	
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
  double Upp2, Up1, Up2;
 
  ST_tensor = aloc4(3,3,nsd,nne);
  ST = aloc1(3*3*nsd*nne);
  
  double F1[9], F2[9]; 
  double **F_mat1, **F_mat2;

  F_mat1 = aloc2(3,3);
  F_mat2 = aloc2(3,3);
                        
  Matrix(double) Np, Nt;
  Matrix_construct_init(double,Np,npres,1,0.0);
  Matrix_construct_init(double,Nt,nVol, 1,0.0);    
                         
  FEMLIB_initialization(&fe, itg_order, 1, nne);
  FEMLIB_set_element(&fe, xe, ii);
  for(int a = 1; a<=fe.nint; a++)
  {
    FEMLIB_elem_basis_V(&fe, a);  
    FEMLIB_elem_shape_function(&fe, a,npres,Np);
    FEMLIB_elem_shape_function(&fe, a,nVol, Nt);    
         
    if(npres==1)
      Mat_v(Nt,1,1) = 1.0;
    
    if(nVol==1)
      Mat_v(Np,1,1) = 1.0;
    
    
    double Tn1 = 0.0; // 1: n-1+alpha
    double Tn2 = 0.0; // 2: n+alpha
    double Pn1 = 0.0;
    double Pn2 = 0.0;
    
    for(int a=0; a<nVol; a++)
    {
      Tn1 += Mat_v(Nt,a,1)*((1.0-alpha)*eps[ii].T[a*3+2] + alpha*eps[ii].T[a*3+1]);
      Tn2 += Mat_v(Nt,a,1)*((1.0-alpha)*eps[ii].T[a*3+1] + alpha*eps[ii].T[a*3+0]);          
    }
    
    for(int a=0; a<npres; a++)
    {
      Pn1 += Mat_v(Np,a,1)*P1[a];
      Pn2 += Mat_v(Np,a,1)*P2[a];          
    }
    
    UPP(Tn2,&hommat[mat],&Upp2);
    Upp2 = Upp2*kappa;

    UP(Tn1, &hommat[mat], &Up1);
    UP(Tn2, &hommat[mat], &Up2);
    Up1 = Up1*kappa;
    Up2 = Up2*kappa;        
            
    shape_tensor(nne,ndn,fe.temp_v.N_x.m_pdata,fe.temp_v.N_y.m_pdata,fe.temp_v.N_z.m_pdata,ST_tensor);
    shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);
    
    def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u1,F_mat1);
    mat2array(F1,CONST_2(double) F_mat1,3,3);
    
    def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u2,F_mat2);
    mat2array(F2,CONST_2(double) F_mat2,3,3);        

/*
				add_3F_Kpt_ip(Kpt,nVol,npres,detJ,wt,Nt,Np,dt,alpha);	
				add_3F_Ktt_ip(Ktt,nVol,detJ,wt,Nt,Upp2,dt,alpha);	

        if(npres==nne)
        {
          add_3F_Kut_ip(Kut,nne,nVol,ST,F2,detJ,wt,Nt,dt,alpha);
        }	
        
				if(npres==1)
				{
          add_3F_Kup_ip(Kup,nne,npres,ST,F2,detJ,wt,Np,dt,alpha);
	  			add_3F_Ktp_ip(Ktp,nVol,npres,detJ,wt,Nt,Np,dt,alpha);			
	  		}	

				resid_w_inertia_Rp_ip(fp1, npres, F1, detJ, wt, Np, Tn1);
				resid_w_inertia_Rt_ip(ft1, nVol, detJ, wt, Nt, Pn1, Up1);

				resid_w_inertia_Rp_ip(fp2, npres, F2, detJ, wt, Np, Tn2);
				resid_w_inertia_Rt_ip(ft2, nVol, detJ, wt, Nt, Pn2, Up2);
				
				resid_w_inertia_Ru_ip(fu1, nne, ST, F1, detJ, wt, Pn1);				
				resid_w_inertia_Ru_ip(fu2, nne, ST, F2, detJ, wt, Pn2);			
	*/			
// \/this is for test ////////////////////////////////////////////////////////////			  
//        double F_I[9];
//        double Jn = det3x3(F2);
//        inverse(F2,3,F_I);  
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
//        3,3,3,1.0,F2,3,F2,3,0.0,C,3);
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
//
//        const damage *ptrDam = &(eps[ii].dam[ip]);
//        disp_based_resid_at_ip(Ru_,nne,ST,F2,S,ptrDam,detJ,wt);
//        HW_Ru_at_ip(Ru,nne,nne,ST,Fn,F2,F_I,
//        		1.0,0.0,Jn,Tn2,S,Pn2,detJ,wt,0);
//		    HW_Rp_at_ip(Rp,npres,Np,1.0,Jn,1.0,Tn2,detJ,wt);
//        HW_Rt_at_ip(Rt,nVol,Nt,Tn2,Jn,Pn2,kappa,Up2/kappa,detJ,wt);		    
//		    
//        if(i==0)
//        {
//          for(int a=0; a<=nVol; a++)
//            printf("%e, %e, %e\n", ft2[a], Rt[a],ft2[a]-Rt[a]);
//        }		
//		
//		    ip++;
//				
// /\this is for test ////////////////////////////////////////////////////////////												
								
  }
   	
   	
// this is for test ////////////////////////////////////////////////////////////
//  dealoc1(Ru);
//  dealoc1(Ru_);
//  dealoc1(Rp);
//  dealoc1(Rt);  
//
// this is for test ////////////////////////////////////////////////////////////  
 
/*    	
  for(int a=0; a<nne*nsd; a++)
  	fu[a] = fuI[a] + (-(1.0 - alpha)*dt*fu2[a] - alpha*dt*fu1[a]);  	
  	
  for(int a=0; a<nVol; a++)
  	ft[a] = ((1.0 - alpha)*dt*ft2[a] + alpha*dt*ft1[a]);

  for(int a=0; a<npres; a++)
  	fp[a] = (-(1.0 - alpha)*dt*fp2[a] - alpha*dt*fp1[a]);
  */  	  	  	 
  dealoc4(ST_tensor,3,3,nsd);
  free(ST);
  
  dealoc2(F_mat1,3);  
  dealoc2(F_mat2,3); 
  
//  condense_F_out(fe,nne,nsd,npres,nVol,fu,ft,fp,Kut,Kup,Ktp,Ktt,Kpt);

  FEMLIB_destruct(&fe);                   
  Matrix_cleanup(xe);
  Matrix_cleanup(Np);  
  Matrix_cleanup(Nt);
  
  Matrix_cleanup(fu);
  Matrix_cleanup(fu1);
  Matrix_cleanup(fu2);
  Matrix_cleanup(fp); 
  Matrix_cleanup(fp1);
  Matrix_cleanup(fp2);
  Matrix_cleanup(ft); 
  Matrix_cleanup(ft1);
  Matrix_cleanup(ft2);
  Matrix_cleanup(Kut);
  Matrix_cleanup(Kpt);
  Matrix_cleanup(Ktt);
  Matrix_cleanup(Kup);
  Matrix_cleanup(Ktp);
}

void resid_w_inertia_3f_el(double *fe,
        const int ii,
        double *fuI,
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
  
  double alpha_1 = 1.0 - alpha;
  double alpha_2 = alpha;
  double dt_alpha_1_minus_alpha = dt*alpha_1*alpha_2;  
    
  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UP = getDUdJFunc(1, &hommat[mat]);  
  UPP = getD2UdJ2Func(1, &hommat[mat]);
  
  int err = 0;  
  int ndofe = nne*ndofn;

  double *fu  = aloc1(nne*nsd);
  double *fu1 = aloc1(nne*nsd);
  double *fu2 = aloc1(nne*nsd);  
  double *fp  = aloc1(npres);
  double *fp1 = aloc1(npres);
  double *fp2 = aloc1(npres);
  double *ft  = aloc1(nVol);   
  double *ft1 = aloc1(nVol);
  double *ft2 = aloc1(nVol);  
  double *Kut = aloc1(nne*nsd*nVol);
  double *Kpt = aloc1(npres*nVol);
  double *Ktt = aloc1(nVol*nVol);
  
  double *Kup = aloc1(nne*nsd*npres);
  double *Ktp = aloc1(nVol*npres); 

  memset(fu, 0,      nne*nsd*sizeof(double));
  memset(fp, 0,        npres*sizeof(double));
  memset(ft, 0,         nVol*sizeof(double));      
  memset(fu1,0,      nne*nsd*sizeof(double));
  memset(fp1,0,        npres*sizeof(double));
  memset(ft1,0,         nVol*sizeof(double));    
  memset(fu2,0,      nne*nsd*sizeof(double));
  memset(fp2,0,        npres*sizeof(double));
  memset(ft2,0,         nVol*sizeof(double));      
	memset(Kut,0, nne*nsd*nVol*sizeof(double));
	memset(Kpt,0,   npres*nVol*sizeof(double));    
	memset(Ktt,0,    nVol*nVol*sizeof(double));	
	
	memset(Kup,0,nne*nsd*npres*sizeof(double));
	memset(Ktp,0,   nVol*npres*sizeof(double)); 	
  
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
  integrate(itg_order ,&npt_x,&npt_y,&npt_z,
          int_pt_ksi,int_pt_eta,int_pt_zet,
          weights);
  
  double **F_mat1, **F_mat2, *F1, *F2;
  
  F1 = aloc1(9);
  F2 = aloc1(9);
  F_mat1 = aloc2(3,3);
  F_mat2 = aloc2(3,3);
  
  
// \/this is for test ////////////////////////////////////////////////////////////			  
//
//  double *Ru = aloc1(nne*nsd);
//  double *Ru_ = aloc1(nne*nsd);
//  double *Rp = aloc1(nne*nsd);
//  double *Rt = aloc1(nne*nsd);    
//  memset(Ru_,0,nne*nsd*sizeof(double));        
//  memset(Ru,0,nne*nsd*sizeof(double));  
//  memset(Rp,0,npres*sizeof(double));        
//  memset(Rt,0,nVol*sizeof(double));    
//  int ip = 0;  
// 
// /\this is for test ////////////////////////////////////////////////////////////			
   
  for(long i=0; i<npt_x; i++)
  {
    for(long j=0; j<npt_y; j++)
    {
      for(long k=0; k<npt_z; k++)
      {     	
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nne,Na);
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],npres,Np);
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nVol,Nt);
        
        if(npres==1)
          Nt[0] = 1.0;
        
        if(nVol==1)
          Np[0] = 1.0;
        
        double detJ = deriv(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nne,x,y,z,N_x,N_y,N_z);
        double wt = weights[k];
        
        double Tn1 = 0.0; // 1: n-1+alpha
        double Tn2 = 0.0; // 2: n+alpha
        double Pn1 = 0.0;
        double Pn2 = 0.0;
       
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


        if(alpha<1.0e-15 || fabs(alpha-1.0)<1.0e-15)
        {        
				  resid_w_inertia_Rp_ip(fp1, npres, F1, detJ, wt, Np, Tn1);
				  resid_w_inertia_Rt_ip(ft1, nVol, detJ, wt, Nt, Pn1, Up1);
				  resid_w_inertia_Ru_ip(fu1, nne, ST, F1, detJ, wt, Pn1);
				}
				else
				{  				  
  				add_3F_Kpt_ip(Kpt,nVol,npres,detJ,wt,Nt,Np,dt_alpha_1_minus_alpha);	
	  			add_3F_Ktt_ip(Ktt,nVol,detJ,wt,Nt,Upp2,dt_alpha_1_minus_alpha);	

          if(npres==nne)
          {
            add_3F_Kut_ip(Kut,nne,nVol,ST,F2,detJ,wt,Nt,dt_alpha_1_minus_alpha);
          }	
        
				  if(npres==1)
				  {
            add_3F_Kup_ip(Kup,nne,npres,ST,F2,detJ,wt,Np,dt_alpha_1_minus_alpha);
	  			  add_3F_Ktp_ip(Ktp,nVol,npres,detJ,wt,Nt,Np,dt_alpha_1_minus_alpha);			
	  		  }				  
				  
  				resid_w_inertia_Rp_ip(fp2, npres, F2, detJ, wt, Np, Tn2);
	  			resid_w_inertia_Rt_ip(ft2, nVol, detJ, wt, Nt, Pn2, Up2);				
		  		resid_w_inertia_Ru_ip(fu2, nne, ST, F2, detJ, wt, Pn2);
				  resid_w_inertia_Rp_ip(fp1, npres, F1, detJ, wt, Np, Tn1);
				  resid_w_inertia_Rt_ip(ft1, nVol, detJ, wt, Nt, Pn1, Up1);
				  resid_w_inertia_Ru_ip(fu1, nne, ST, F1, detJ, wt, Pn1);		  		
		  	}			
				
// \/this is for test ////////////////////////////////////////////////////////////			  
//        double F_I[9];
//        double Jn = det3x3(F2);
//        inverse(F2,3,F_I);  
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
//        3,3,3,1.0,F2,3,F2,3,0.0,C,3);
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
//
//        const damage *ptrDam = &(eps[ii].dam[ip]);
//        disp_based_resid_at_ip(Ru_,nne,ST,F2,S,ptrDam,detJ,wt);
//        HW_Ru_at_ip(Ru,nne,nne,ST,Fn,F2,F_I,
//        		1.0,0.0,Jn,Tn2,S,Pn2,detJ,wt,0);
//		    HW_Rp_at_ip(Rp,npres,Np,1.0,Jn,1.0,Tn2,detJ,wt);
//        HW_Rt_at_ip(Rt,nVol,Nt,Tn2,Jn,Pn2,kappa,Up2/kappa,detJ,wt);		    
//		    
//        if(i==0)
//        {
//          for(int a=0; a<=nVol; a++)
//            printf("%e, %e, %e\n", ft2[a], Rt[a],ft2[a]-Rt[a]);
//        }		
//		
//		    ip++;
//				
// /\this is for test ////////////////////////////////////////////////////////////												
								
      }
    }
  }
   	
   	
// this is for test ////////////////////////////////////////////////////////////
//  dealoc1(Ru);
//  dealoc1(Ru_);
//  dealoc1(Rp);
//  dealoc1(Rt);  
//
// this is for test ////////////////////////////////////////////////////////////  
 
    	
  for(int a=0; a<nne*nsd; a++)
  	fu[a] = fuI[a] + (-(1.0 - alpha)*dt*fu2[a] - alpha*dt*fu1[a]);  	
  	
  for(int a=0; a<nVol; a++)
  	ft[a] = ((1.0 - alpha)*dt*ft2[a] + alpha*dt*ft1[a]);

  for(int a=0; a<npres; a++)
  	fp[a] = (-(1.0 - alpha)*dt*fp2[a] - alpha*dt*fp1[a]);
    	  	  	 
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
  
  free(fu1); free(fu2);
  free(fp1); free(fp2);
  free(ft1); free(ft2);

  condense_F_out(fe,nne,nsd,npres,nVol,fu,ft,fp,Kut,Kup,Ktp,Ktt,Kpt);
                   
  free(fu);
  free(fp);
  free(ft);
  free(Kup);
  free(Kut); 
  free(Ktp);  
  free(Kpt); 
  free(Ktt);
}

void DISP_resid_w_inertia_el(double *f,
         const int ii,
         const int ndofn,
         const int nne,
         const double *x,
         const double *y,
         const double *z,		     
         const ELEMENT *elem,
         const HOMMAT *hommat,
		     const NODE *node, double dt, double t,
		     double *r_2, double* r_1, double *r_0, double alpha)
{
  const int mat = elem[ii].mat[2];
  double rho = hommat[mat].density;
  int ndofe = nne*ndofn;
    
  /* make sure the f vector contains all zeros */
  memset(f,0,ndofe*sizeof(double));
  
  FEMLIB fe;
  Matrix(double) xe, du;
  
  Matrix_construct_redim(double,xe,nne,3);
  Matrix_construct_redim(double,du,3,1);  
  
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
  for(int a = 1; a<=fe.nint; a++)
  {
    FEMLIB_elem_basis_V(&fe, a); 

    Matrix_init(du, 0.0);  
    for(long a = 0; a<nne; a++)
    {                              
      for(long b = 0; b<ndofn; b++)
      {
        long id = a*ndofn + b;
        Mat_v(du,b+1,1) += Mat_v(fe.N,a+1,1)*(r_2[id]-2.0*r_1[id]+r_0[id]);
      }
    }
    
	  for(long a = 0; a<nne; a++)
	  {
      for(long b=0; b<ndofn; b++)
	    {
	      long id = a*ndofn + b;
	      f[id] += rho/dt*Mat_v(fe.N,a+1,1)*Mat_v(du,b+1,1)*fe.detJxW;	      	      
	    }
	  }
	          
  }
        
  Matrix_cleanup(xe);  
  Matrix_cleanup(du); 
  FEMLIB_destruct(&fe);
}

int residuals_w_inertia_el(double *fe, int i, 
			int nne, long ndofn, long ndofe, double *r_e,                               
		  NODE *node, ELEMENT *elem, HOMMAT *hommat, SUPP sup, EPS *eps, SIG *sig,
		  long* nod, long *cn, double *x, double *y, double *z,                                
		  double dt, double t, double alpha, double *r_n, double *r_n_1, double *r_n_1_a, double *r_n_a)
{
	int nsd = 3;
	int err = 0;
	double *r0, *r0_;
	
	double *f_i     = aloc1(ndofe);
	double *f_n_a   = aloc1(ndofe);
	double *f_n_1_a = aloc1(ndofe);
	
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
		
	mid_point_rule(r_n_1_a,r0_,r0,  alpha, ndofe); 
	mid_point_rule(r_n_a,  r0, r_e, alpha, ndofe); 
	  
	err =  DISP_resid_el(f_n_1_a,i,ndofn,nne,x,y,z,elem,
	         hommat,nod,node,eps,sig,sup,r_n_1_a);                                             
	       
	err =  DISP_resid_el(f_n_a,i,ndofn,nne,x,y,z,elem,
	        hommat,nod,node,eps,sig,sup,r_n_a);

	DISP_resid_w_inertia_el(f_i,i,ndofn,nne,x,y,z,elem,hommat,node,dt,t,r_e, r0, r0_, alpha); 
	
	for(long a = 0; a<ndofe; a++)
	    fe[a] = -f_i[a] - (1.0-alpha)*dt*f_n_a[a] - alpha*dt*f_n_1_a[a];
	
	free(r0);
	free(r0_);
	free(f_i);
	free(f_n_a);
	free(f_n_1_a);		
	return err;
} 

int fd_residuals (double *f_u,
		  long ne,
		  int n_be,
		  long ndofn,
		  long npres,
		  double *d_r,
		  double *r,
		  NODE *node,
		  ELEMENT *elem,
		  BOUNDING_ELEMENT *b_elems,
		  MATGEOM matgeom,
		  HOMMAT *hommat,
		  SUPP sup,
		  EPS *eps,
		  SIG *sig,
		  double nor_min,
		  CRPL *crpl,
		  double dt,
		  double t,
		  double stab,
		  long nce,
		  COEL *coel,
		  MPI_Comm mpi_comm,
		  const PGFem3D_opt *opts,
		  double alpha, double *r_n, double *r_n_1)
{
/* make decision to include ineria*/
  const int mat = elem[0].mat[2];
  double rho = hommat[mat].density;
  long include_inertia = 1;
  
  if(fabs(rho)<1.0e-15)
    include_inertia = 0;
   
/* decision end*/   
  int err = 0;
  /* long i,j,nne,ndofe,ndofc,k,kk,II,*nod,P,R,*cn; */
  /* double *r_e,*x,*y,*z,*fe,*X,*Y; */

  /* int nne_t; */

  int myrank,nproc;
  MPI_Comm_size(mpi_comm,&nproc);
  MPI_Comm_rank(mpi_comm,&myrank);
  
  for (int i=0;i<ne;i++){
    /* Number of element nodes */
    int nne = elem[i].toe;
    int nne_t = nne + elem[i].n_bub;
    /* Nodes on element */
    long *nod = aloc1l (nne);
    elemnodes (i,nne,nod,elem);
    /* Element Dof */
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node);

    /* allocation */
    double *r_e = aloc1 (ndofe);

    double *x,*y,*z;
    if(opts->analysis_type == MINI
       || opts->analysis_type == MINI_3F){/* P1+B/P1 */
      x = aloc1 (nne_t);
      y = aloc1 (nne_t);
      z = aloc1 (nne_t);
    } else {
      x = aloc1 (nne);
      y = aloc1 (nne);
      z = aloc1 (nne);
    }

    double *fe = aloc1 (ndofe);
    long *cn = aloc1l (ndofe);
    /* code numbers on element */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn);
    
    /* coordinates */
    if(sup->multi_scale)
    {
	    nodecoord_total(nne,nod,node,x,y,z);
	    def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
    } 
    else
    {
      switch(opts->analysis_type)
      {
      case DISP: 
      case TF:
	      nodecoord_total(nne,nod,node,x,y,z);
	      def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
	      break;
      default:
	      nodecoord_updated(nne,nod,node,x,y,z);
	      /* deformation on element */
	      def_elem (cn,ndofe,d_r,elem,node,r_e,sup,0);    
	      break;
      }
    }

    if(opts->analysis_type == MINI
       || opts->analysis_type == MINI_3F){/* P1+B/P1 */
      element_center(nne,x,y,z);
    }

		int nVol = N_VOL_TF;
    int nsd = 3;

    if(include_inertia)
    {
   		double *r_n_a, *r_n_1_a;
   		r_n_a = aloc1(ndofe);
   		r_n_1_a = aloc1(ndofe);
   		      
    	switch(opts->analysis_type)
    	{   		
      	case DISP:
      	{
					/* Get TOTAL deformation on element; r_e already contains
	   			INCREMENT of deformation, add the deformation from previous. */
///////////////////////////////////////////////////////////////////////////////////
					residuals_w_inertia_el(fe,i,nne,ndofn,ndofe,r_e,node,elem,hommat,sup,eps,sig,
		  			nod,cn,x,y,z,dt,t,alpha,r_n,r_n_1, r_n_1_a, r_n_a);		
      	}
      	break;
      	case TF:
      	{
					double *u1, *u2, *P1, *P2;
					double *f_i;
					f_i = aloc1(ndofe);
					u1 = aloc1(nne*nsd);
					u2 = aloc1(nne*nsd);
		      P1 = aloc1(npres);		  			
		      P2 = aloc1(npres);		      
		      
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
        	
        	mid_point_rule(r_n_1_a,r0_,r0,  alpha, ndofe); 
	        mid_point_rule(r_n_a,  r0, r_e, alpha, ndofe);         	
				      
		      DISP_resid_w_inertia_el(f_i,i,ndofn,nne,x,y,z,elem,hommat,node,dt,t,r_e, r0, r0_, alpha);
		          		  
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
            P1[0] =  (1.0-alpha)*eps[i].d_T[2] + alpha*eps[i].d_T[1]; 
            P2[0] =  (1.0-alpha)*eps[i].d_T[1] + alpha*eps[i].d_T[0];             					  
				  }   
					resid_w_inertia_3f_el(fe,i,f_i,ndofn,nne,npres,nVol,nsd,
					                      x,y,z,elem,hommat,node,u2,u1,P2,P1,dt,sig,eps,alpha);        							  								
					free(r_n_a);
					free(r_n_1_a);
					free(P1); free(P2); free(u1); free(u2);
					free(f_i);
    		} 
    		break;
    		default:
      		err = resid_on_elem(i,ndofn,nne,nod,elem,node,matgeom,
			   											hommat,x,y,z,eps,sig,r_e,npres,
			   											nor_min,fe,crpl,dt,opts->analysis_type);
      	break;    		  
      }        
    }
    else
    {    
    /* Residuals on element */
    switch(opts->analysis_type){
    case STABILIZED:
      err = resid_st_elem (i,ndofn,nne,elem,nod,node,hommat,
			   x,y,z,eps,sig,sup,r_e,nor_min,fe,dt,stab);
      break;
    case MINI:
      err = MINI_resid_el(fe,i,ndofn,nne,x,y,z,elem,
			  nod,node,hommat,eps,sig,r_e);
      break;
    case MINI_3F:
      err = MINI_3f_resid_el(fe,i,ndofn,nne,x,y,z,elem,
			     nod,node,hommat,eps,sig,r_e);
      break;
    case DISP:
      err =  DISP_resid_el(fe,i,ndofn,nne,x,y,z,elem,
			   hommat,nod,node,eps,sig,sup,r_e);
      break;
    case TF:
{

      err =  DISP_resid_el(fe,i,ndofn,nne,x,y,z,elem,
			   hommat,nod,node,eps,sig,sup,r_e);
	         			         
//       residuals_3f_el(fe,i,ndofn,nne,npres,nVol,nsd,
//			                      x,y,z,elem,hommat,nod,node,dt,sig,eps,sup,r_e);
			break;                      
}      
    default:
      err = resid_on_elem (i,ndofn,nne,nod,elem,node,matgeom,
			   hommat,x,y,z,eps,sig,r_e,npres,
			   nor_min,fe,crpl,dt,opts->analysis_type);
      break;
    }
  }

    /* Assembly */
    {
      int j = 0;
      int II = 0;
      for (int k=0;k<nne;k++){
	for (int kk=0;kk<node[nod[k]].ndofn;kk++){
	  II = node[nod[k]].id[kk]-1;
	  if (II < 0) continue;
	  f_u[II] += fe[j+kk];
	}/*end kk*/
	j += node[nod[k]].ndofn;
      }/*end k*/
    }

    dealoc1l (nod);
    dealoc1 (r_e);
    dealoc1 (x);
    dealoc1 (y);
    dealoc1 (z);
    dealoc1 (fe);
    dealoc1l (cn);

    /*** RETURN on error ***/
    if(err != 0) return err;

  }/* end i < ne*/

  /**** COHESIVE ELEMENT RESIDUALS ****/
  if (opts->cohesive == 1){
    int ndofc = 3;

    for (int i=0;i<nce;i++){
      
      int nne = coel[i].toe/2;
      int ndofe = coel[i].toe*ndofc;
      
      long *nod = aloc1l (coel[i].toe);
      double *r_e = aloc1 (ndofe);
      double *x = aloc1 (coel[i].toe);
      double *y = aloc1 (coel[i].toe);
      double *z = aloc1 (coel[i].toe);
      double *fe = aloc1 (ndofe);
      long *cn = aloc1l (ndofe);
      
      for (int j=0;j<coel[i].toe;j++)
	nod[j] = coel[i].nod[j];

      nodecoord_updated (coel[i].toe,nod,node,x,y,z);
      
      /* code numbers on element */
      get_dof_ids_on_elem_nodes(0,coel[i].toe,ndofc,nod,node,cn);
      
      /* deformation on element */
      def_elem (cn,ndofe,d_r,elem,node,r_e,sup,0);

      /* Residuals on element */
      resid_co_elem (i,ndofc,nne,nod,x,y,z,coel,r_e,fe,nor_min,myrank);
      
      /* Localization */
      {
	int II = 0;
	for (int k=0;k<coel[i].toe;k++){
	  for (int kk=0;kk<ndofc;kk++){
	    II = node[nod[k]].id[kk]-1;
	    if (II < 0)
	      continue;
	    f_u[II] += fe[k*ndofc+kk];
	  }/*end kk*/
	}/*end k*/
      }
      
      dealoc1l (nod);
      dealoc1 (r_e);
      dealoc1 (x);
      dealoc1 (y);
      dealoc1 (z);
      dealoc1 (fe);
      dealoc1l (cn);
    }/* end i < nce */
  }/* end coh == 1 */

  /*===============================================
    |             BOUNDARY ELEMENTS               |
    ===============================================*/
    
  for (int i=0; i<n_be; i++){

    /* get the coordinates and dof id's on the element nodes */
    const long *ptr_vnodes = elem[b_elems[i].vol_elem_id].nod; /* --"-- */
    const ELEMENT *ptr_elem = &elem[b_elems[i].vol_elem_id]; /* --"-- */
    const BOUNDING_ELEMENT *ptr_be = &b_elems[i];
    const int nne = ptr_elem->toe;

    double *x = aloc1(nne);
    double *y = aloc1(nne);
    double *z = aloc1(nne);

    switch(opts->analysis_type){
    case DISP:
      nodecoord_total(nne,ptr_vnodes,node,x,y,z);
      break;
    default:
      nodecoord_updated(nne,ptr_vnodes,node,x,y,z);
      break;
    }

    int ndofe = get_ndof_on_bnd_elem(node,ptr_be,elem);
    double *r_e = aloc1(ndofe);
    double *fe = aloc1(ndofe);
    long *cn = aloc1l(ndofe);
    long *Gcn = aloc1l(ndofe);

    get_dof_ids_on_bnd_elem(0,ndofn,node,ptr_be,elem,cn);
    get_dof_ids_on_bnd_elem(1,ndofn,node,ptr_be,elem,Gcn);

    /* compute the deformation on the element */
    def_elem(cn,ndofe,d_r,NULL,NULL,r_e,sup,0);

    /* TOTAL LAGRANGIAN formulation */
    if(opts->analysis_type == DISP){
      double *r_en = aloc1(ndofe);
      def_elem(cn,ndofe,r,NULL,NULL,r_en,sup,1);
      vvplus(r_e,r_en,ndofe);
      free(r_en);
    }

    /* for debugging */
    double *RR = aloc1(ndofe);
    def_elem(cn,ndofe,f_u,NULL,NULL,RR,sup,2);

    if(opts->analysis_type == DISP){
      err += DISP_resid_bnd_el(fe,i,ndofn,ndofe,x,y,z,b_elems,elem,
			       hommat,node,eps,sig,sup,r_e);
    } else {
      /* not implemented/needed so do nothing */
    }

    /* Assembly to local part of the residual vector */
    {
      for(int j=0; j<ndofe; j++){
	int II = cn[j] - 1;
	if(II >= 0){
	  f_u[II] += fe[j];
	}
      }
    }

    def_elem(cn,ndofe,f_u,NULL,NULL,RR,sup,2);

    free(x);
    free(y);
    free(z);

    free(r_e);
    free(fe);
    free(cn);
    free(Gcn);
    free(RR);

    /*** RETURN on error ***/
    if(err != 0) return err;
  } /* for each bounding element */
  return err;
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
    			    
  
  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UP = getDUdJFunc(1, &hommat[mat]);  
  UPP = getD2UdJ2Func(1, &hommat[mat]);
  
  const int nsd = 3;
  const int nVol = N_VOL_TF;
  int err = 0;  
  int ndofe = nne*ndofn;

  double *fp  = aloc1(npres);   
  
  double *ft  = aloc1(nVol);   
    
  double *Ktu = aloc1(nne*nsd*nVol);
  double *Ktp = aloc1(npres*nVol);
  double *Ktt = aloc1(nVol*nVol);
  
  double *Kpu = aloc1(npres*nne*nsd);
  double *Kpt = aloc1(npres*nVol);  

  memset(fp, 0,        npres*sizeof(double));      

  memset(ft, 0,        nVol*sizeof(double));      
  
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
        add_3F_Kpu_ip(Kpu,nne,npres, ST,F,detJ,wt,Np,-1.0);        
				add_3F_Kpt_ip(Kpt,nVol,npres,detJ,wt,Nt,Np,-1.0);				

				resid_w_inertia_Rt_ip(ft, nVol, detJ, wt, Nt, Pn, Up);							
				resid_w_inertia_Rp_ip(fp, npres, F, detJ, wt, Np, Tn);
				
        update_3f_state_variables_ip(ii,ip,elem,hommat,sig,eps,F,Pn,Tn,volume,detJ*wt);
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
    
  double *theta   = aloc1(nVol);
  double *KptI    = aloc1(nVol*nVol);
  
  inverse(Kpt,nVol,KptI);
  
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nne*nsd,1.0,Kpu,nne*nsd,r_e,1,1.0,fp,1);  
              
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
//              nVol,1,nne*nsd,1.0,Ktu,nne*nsd,r_e,1,1.0,ft,1);               

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
    
  if(opts->analysis_type != TF)
    return;
        
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
    def_elem (cn,ndofe,d_r,elem,node,r_e,sup,0);
    
  	int nsd = 3;
	  double *r_en;
	
    r_en = aloc1(ndofe);

	  def_elem(cn,ndofe,r,elem,node,r_en,sup,1);
	  vvplus(r_e,r_en,ndofe);
	  
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
