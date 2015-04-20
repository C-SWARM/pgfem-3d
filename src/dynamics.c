
#include "femlib.h"
#include "dynamics.h"
#include "utils.h"
#include "enumerations.h"
#include "allocation.h"
#include "displacement_based_element.h"
#include "three_field_element.h"

void MMS_body_force(double *b, HOMMAT const * hommat, double t, double X, double Y, double Z)
{
  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;
}

void stiffmat_disp_w_inertia_el(double *Ks,
         const int ii,
         const int ndofn,
         const int nne, const int npres, const int nVol, const int nsd,
         const double *x, const double *y, const double *z,		     
         const ELEMENT *elem, const HOMMAT *hommat, const long *nod, const NODE *node, double dt,
         SIG *sig, EPS *eps, const SUPP sup, const int analysis,		     
		     double alpha, double *r_n, double *r_e)
{
  int err = 0;
  const int mat = elem[ii].mat[2];
  double rho = hommat[mat].density;  
   
  int ndofe = nne*ndofn;

  Matrix(double) Kuu_K,Kuu_I, u, u_n;
  Matrix_construct_init(double,Kuu_I,ndofe,ndofe,0.0);
  Matrix_construct_init(double,Kuu_K,ndofe,ndofe,0.0);  
  Matrix_construct_init(double,u,ndofe,1,0.0);
  Matrix_construct_init(double,u_n,ndofe,1,0.0);      

  /* make sure the stiffenss matrix contains all zeros */  
  memset(Ks,0,ndofe*ndofe*sizeof(double));
        
  for (long I=0;I<nne;I++)
  {
    for(long J=0; J<ndofn; J++)
      Mat_v(u_n,I*ndofn+J+1,1) = r_n[nod[I]*ndofn + J];  
  }
   
  mid_point_rule(u.m_pdata, u_n.m_pdata, r_e, alpha, ndofe); 
  
  FEMLIB fe;
  Matrix(double) xe;
  
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

  if(analysis == DISP || analysis == TF)    
  {       
    FEMLIB_initialization(&fe, itg_order, 1, nne);
    FEMLIB_set_element(&fe, xe, ii);      
    for(int ip = 1; ip<=fe.nint; ip++)
    {
      FEMLIB_elem_basis_V(&fe, ip); 
    
      for(long a = 0; a<nne; a++)
      {
        for(long c=0; c<nne; c++)
        {
          for(long b=1; b<=ndofn; b++)
            Mat_v(Kuu_I,a*ndofn+b,c*ndofn+b) += rho/dt*Vec_v(fe.N,a+1)*Vec_v(fe.N,c+1)*fe.detJxW;
        }
	    } 
    }
    
    Matrix_cleanup(xe);  
    FEMLIB_destruct(&fe);
  }
  
  switch(analysis)
  {
    case DISP:
      err = DISP_stiffmat_el(Kuu_K.m_pdata,ii,ndofn,nne,x,y,z,elem,
                             hommat,nod,node,eps,sig,sup,u.m_pdata);

      for(long a = 0; a<ndofe*ndofe; a++)
        Ks[a] = -Kuu_I.m_pdata[a]-alpha*(1.0-alpha)*dt*Kuu_K.m_pdata[a];                                
      
      break;

    case TF:
      if(0<alpha && alpha<1.0)
      {
        stiffmat_3f_el(Kuu_K.m_pdata,ii,ndofn,nne,npres,nVol,nsd,x,y,z,
                                 elem,hommat,nod,node,dt,
                                 sig,eps,sup,alpha,u.m_pdata);                                                              
      }                          
      for(long a = 0; a<ndofe*ndofe; a++)
        Ks[a] = -Kuu_I.m_pdata[a] + Kuu_K.m_pdata[a];                                

      break;
      
    default:
      printf("Only displacement based element and three field element are supported\n");
      break;                         
  }              

  Matrix_cleanup(Kuu_I); 
  Matrix_cleanup(Kuu_K);  
  Matrix_cleanup(u); 
  Matrix_cleanup(u_n);  
}


void DISP_resid_body_force_el(double *f,
         const int ii,
         const int ndofn,
         const int nne,
         const double *x,
         const double *y,
         const double *z,		     
         const ELEMENT *elem,
         const HOMMAT *hommat,
		     const NODE *node, double dt, double t)
{
  const int mat = elem[ii].mat[2];
  double rho = hommat[mat].density;
  int ndofe = nne*ndofn;
    
  /* make sure the f vector contains all zeros */
  memset(f,0,ndofe*sizeof(double));
  
  FEMLIB fe;
  Matrix(double) xe;
  
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
    
  double *bf = aloc1(ndofn);
                 
  FEMLIB_initialization(&fe, itg_order, 1, nne);
  FEMLIB_set_element(&fe, xe, ii);      
  for(int a = 1; a<=fe.nint; a++)
  {
    FEMLIB_elem_basis_V(&fe, a); 
    double X[3];
    X[0] = X[1] = X[2] = 0.0;
    
    memset(bf,    0,ndofn*sizeof(double));
              
    for(long a = 0; a<nne; a++)
    {
      X[0] += Vec_v(fe.N,a+1)*x[a];
      X[1] += Vec_v(fe.N,a+1)*y[a];          
      X[2] += Vec_v(fe.N,a+1)*z[a];                                    
    }
    
    MMS_body_force(bf, &hommat[mat], t,  X[0], X[1], X[2]);
        
	  for(long a = 0; a<nne; a++)
	  {
      for(long b=0; b<3; b++)
	    {
	      long id = a*ndofn + b;
        f[id] += bf[b]*Vec_v(fe.N,a+1)*fe.detJxW;	      	      	      
	    }
	  }
	          
  }
        
  dealoc1(bf);
        
  Matrix_cleanup(xe);  
  FEMLIB_destruct(&fe);
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
    
  double *bf0, *bf1, *bf2, *bf_n1a, *bf;
  bf0 = aloc1(ndofn);
  bf1 = aloc1(ndofn);
  bf2 = aloc1(ndofn);
  bf_n1a = aloc1(ndofn);        
  bf     = aloc1(ndofn);
  
               
  FEMLIB_initialization(&fe, itg_order, 1, nne);
  FEMLIB_set_element(&fe, xe, ii);      
  for(int a = 1; a<=fe.nint; a++)
  {
    FEMLIB_elem_basis_V(&fe, a); 

    Matrix_init(du, 0.0);
    double X[3];
    X[0] = X[1] = X[2] = 0.0;
    
    memset(bf0,   0,ndofn*sizeof(double));
    memset(bf1,   0,ndofn*sizeof(double));
    memset(bf2,   0,ndofn*sizeof(double));
    memset(bf_n1a,0,ndofn*sizeof(double));
    memset(bf,    0,ndofn*sizeof(double));
              
    for(long a = 0; a<nne; a++)
    {
      X[0] += Vec_v(fe.N,a+1)*x[a];
      X[1] += Vec_v(fe.N,a+1)*y[a];          
      X[2] += Vec_v(fe.N,a+1)*z[a];                                    
      for(long b = 0; b<3; b++)
      {
        long id = a*ndofn + b;
        Vec_v(du,b+1) += Vec_v(fe.N,a+1)*(r_2[id]-2.0*r_1[id]+r_0[id]);
      }
    }
    
    double t1 = t-dt;
    double t0 = t - dt - dt;
    
    if(t0>0) 
      MMS_body_force(bf0, &hommat[mat], t0, X[0], X[1], X[2]);

    if(t1>0)  
      MMS_body_force(bf1, &hommat[mat], t1, X[0], X[1], X[2]); 
    
    MMS_body_force(bf2, &hommat[mat], t,  X[0], X[1], X[2]);
        
    mid_point_rule(bf_n1a, bf0, bf1, alpha, ndofn);
    mid_point_rule(bf, bf1, bf2, alpha, ndofn);	    
    
	  for(long a = 0; a<nne; a++)
	  {
      for(long b=0; b<3; b++)
	    {
	      long id = a*ndofn + b;
	      f[id] += rho/dt*Vec_v(fe.N,a+1)*Vec_v(du,b+1)*fe.detJxW;
        f[id] -= (1.0-alpha)*dt*bf[b]*Vec_v(fe.N,a+1)*fe.detJxW;
        f[id] -= alpha*dt*bf_n1a[b]*Vec_v(fe.N,a+1)*fe.detJxW;	      	      	      
	    }
	  }	          
  }
        
  dealoc1(bf0);
  dealoc1(bf1);
  dealoc1(bf2);
  dealoc1(bf_n1a);        
  dealoc1(bf);         
        
  Matrix_cleanup(xe);  
  Matrix_cleanup(du); 
  FEMLIB_destruct(&fe);
}

int residuals_w_inertia_el(double *fe, int i, 
			int nne, long ndofn, long npres, long nVol,long ndofe, double *r_e,                               
		  NODE *node, ELEMENT *elem, HOMMAT *hommat, SUPP sup, EPS *eps, SIG *sig,
		  long* nod, long *cn, double *x, double *y, double *z,                                
		  double dt, double t, const PGFem3D_opt *opts, double alpha, double *r_n, double *r_n_1)
{
	int nsd = 3;
	int err = 0;
	
  double *r_n_a   = aloc1(ndofe);
  double *r_n_1_a = aloc1(ndofe);   			
	double *r0      = aloc1(ndofe);
	double *r0_     = aloc1(ndofe); 	
	double *f_i     = aloc1(ndofe);
	memset(fe, 0, sizeof(double)*ndofe);
	memset(f_i, 0, sizeof(double)*ndofe);
		
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
	
  if(opts->analysis_type == DISP || opts->analysis_type == TF) 
    DISP_resid_w_inertia_el(f_i,i,ndofn,nne,x,y,z,elem,hommat,node,dt,t,r_e, r0, r0_, alpha);   

  switch(opts->analysis_type)
  {   		
    case DISP:
    {   
	    double *f_n_a   = aloc1(ndofe);
	    double *f_n_1_a = aloc1(ndofe);      
       
	    err =  DISP_resid_el(f_n_1_a,i,ndofn,nne,x,y,z,elem,
	         hommat,nod,node,eps,sig,sup,r_n_1_a);                                             
	       
	    err =  DISP_resid_el(f_n_a,i,ndofn,nne,x,y,z,elem,
	        hommat,nod,node,eps,sig,sup,r_n_a);
            	
	    for(long a = 0; a<ndofe; a++)
	      fe[a] = -f_i[a] - (1.0-alpha)*dt*f_n_a[a] - alpha*dt*f_n_1_a[a];
	      
	    free(f_n_a);
	    free(f_n_1_a);		      
	      
	    break;
	  }  
	  case TF:
	  {  
      if(0<alpha && alpha<1.0)
      {	    
	      residuals_3f_w_inertia_el(fe,i,ndofn,nne,npres,nVol,nsd,x,y,z,elem,hommat,node,
	                                dt,sig,eps,alpha,r_n_a,r_n_1_a);
	    }
	                              
	    for(long a = 0; a<ndofe; a++)
	      fe[a] = -f_i[a];
	      
	    break;
	  }  
    default:
      printf("Only displacement based element and three field element are supported\n");
  	  break;
  }

	free(r_n_a);
	free(r_n_1_a);	
	free(r0);
	free(r0_);
	free(f_i);	
	return err;
} 
