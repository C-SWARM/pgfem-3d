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
      
      DISP_resid_w_inertia_el(f_i,i,ndofn,nne,x,y,z,elem,hommat,node,dt,t,r_e, r0, r0_, alpha); 
	
	    for(long a = 0; a<ndofe; a++)
	      fe[a] = -f_i[a] - (1.0-alpha)*dt*f_n_a[a] - alpha*dt*f_n_1_a[a];
	      
	    free(f_n_a);
	    free(f_n_1_a);		      
	      
	    break;
	  }  
	  case TF:	    
	    DISP_resid_w_inertia_el(f_i,i,ndofn,nne,x,y,z,elem,hommat,node,dt,t,r_e, r0, r0_, alpha);
	    residuals_3f_w_inertia_el(fe,i,f_i,ndofn,nne,npres,nVol,nsd,x,y,z,elem,hommat,node,
	                              dt,sig,eps,alpha,r_n_a,r_n_1_a);
	    break;
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
      residuals_w_inertia_el(fe,i,nne,ndofn,npres,nVol,ndofe,r_e,node,elem,hommat,sup,eps,sig,
		                         nod,cn,x,y,z,dt,t,opts,alpha,r_n,r_n_1);		
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
       residuals_3f_el(fe,i,ndofn,nne,npres,nVol,nsd,
			       x,y,z,elem,hommat,nod,node,dt,sig,eps,sup,r_e);			       			       
			break;                      
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
