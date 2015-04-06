
#include "femlib.h"
#include "dynamics.h"
#include "utils.h"
#include "enumerations.h"
#include "displacement_based_element.h"
#include "three_field_element.h"

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
