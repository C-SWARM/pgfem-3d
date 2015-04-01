#include "femlib.h"
#include "elem3d.h"

#include "data_structure_c.h"
#include "utils.h"
#include "cast_macros.h"
#include "allocation.h"
#include "tensors.h"
#include "def_grad.h"

#ifndef _Matrix_double
Define_Matrix(double);
#define _Matrix_double 1
#endif

#ifndef _Matrix_int
Define_Matrix(int);
#define _Matrix_int 1
#endif


#define LINE           1
#define TRIANGLE       3
#define QUADRILATERAL  4
#define TETRAHEDRON    5
#define HEXAHEDRAL     8
#define QTETRAHEDRON   10

long FEMLIB_determine_integration_type(int e_type, int i_order)
{
  if(i_order == 0)
  {
  	switch(e_type)
  	{
      case LINE:
      	return 1;
      case TRIANGLE:
      	return 1;
      case QUADRILATERAL:
      	return 1;
      case TETRAHEDRON:
      	return 4;
      case HEXAHEDRAL:
      	return 4;
      default:
        return 4;
    }
  }
  
  if(i_order == 1)
  {
  	switch(e_type)
  	{  	
      case LINE:
      	return 2;
      case TRIANGLE:
      	return 2;
      case QUADRILATERAL:
      	return 2;
      case TETRAHEDRON:
      	return 5;
      case HEXAHEDRAL:
      	return 8;
      case QTETRAHEDRON:
      	return 10;      	
      default:
        return 5;
    }        
  }
  return 4;
}

void FEMLIB_set_variable_size(TEMP_VARIABLES *v, int nne)
{
  Matrix_construct_redim(double,v->x  ,nne,1);
  Matrix_construct_redim(double,v->y  ,nne,1);
  Matrix_construct_redim(double,v->z  ,nne,1);
  Matrix_construct_redim(double,v->N_x,nne,1);
  Matrix_construct_redim(double,v->N_y,nne,1);
  Matrix_construct_redim(double,v->N_z,nne,1);         
}

void FEMLIB_destruct_variable(TEMP_VARIABLES *v)
{
  Matrix_cleanup(v->x);  
  Matrix_cleanup(v->y);  
  Matrix_cleanup(v->z);  
  Matrix_cleanup(v->N_x);
  Matrix_cleanup(v->N_y);
  Matrix_cleanup(v->N_z);         
}

void FEMLIB_initialization(FEMLIB *fe, int e_type, int i_order, int nne)
{
  int nsd = 3;
  long nint;
  
	fe->nsd = nsd;
	fe->nne = nne;
  fe->elem_type = e_type;
  fe->intg_order = i_order;

  long intg_type = FEMLIB_determine_integration_type(e_type, i_order);
  int_point(intg_type, &nint);
    
  fe->nint = nint;

  Matrix_construct_redim(double,fe->ksi       ,nint,1);
  Matrix_construct_redim(double,fe->eta       ,nint,1);    
  Matrix_construct_redim(double,fe->zet       ,nint,1);    
  Matrix_construct_redim(double,fe->weights   ,nint,1);
  Matrix_construct_redim(double,fe->N         ,nne ,1);
  Matrix_construct_redim(double,fe->dN        ,nne ,nsd);
  Matrix_construct_redim(double,fe->x_ip      ,nsd ,1);  
  Matrix_construct_redim(double,fe->node_coord,nne ,nsd);
  Matrix_construct_redim(long,fe->node_id,nne,1);
      
  FEMLIB_set_variable_size(&(fe->temp_v), nne);      
    
  long npt_x, npt_y, npt_z;
  integrate(intg_type, &npt_x, &npt_y, &npt_z,
          fe->ksi.m_pdata, fe->eta.m_pdata, fe->zet.m_pdata,
          fe->weights.m_pdata);

  Matrix_construct(int, fe->itg_ids); 
  Matrix_redim(fe->itg_ids, npt_x*npt_y*npt_z, nsd);

  long cnt = 0;
  for(int a=1; a<=npt_x; a++)
  {
  	for(int b=1; b<=npt_y; b++)
    {
    	for(int c=1; c<=npt_z; c++)
    	{
    		cnt++;
    		Mat_v(fe->itg_ids, cnt, 1) = a;
    		Mat_v(fe->itg_ids, cnt, 2) = b;
    		Mat_v(fe->itg_ids, cnt, 3) = c;    		    		
      }
    }
  } 
  
  fe->ST_tensor = aloc4(3,3,nsd,nne);
  fe->ST = aloc1(3*3*nsd*nne);  
}

void FEMLIB_initialization_by_elem(FEMLIB *fe, int e, ELEMENT *elem, NODE *node)
{
  int nne = elem[e].toe;
  
  int itg_order = nne;
  if(nne==4)
    itg_order = nne + 1;   
  
  FEMLIB_initialization(fe, itg_order, 1, nne);
  
  long *nod = aloc1l (nne);
  elemnodes(e,nne,nod,elem);
  
  double *x,*y,*z;
  x = aloc1(nne);
  y = aloc1(nne);
  z = aloc1(nne);
  nodecoord_total(nne,nod,node,x,y,z);
  
  Matrix(double) xe;  
  Matrix_construct_init(double,xe,nne,3,0.0); 
  
  for(int a=0; a<nne; a++)
  {
    Vec_v(fe->node_id, a+1) = nod[a];
    Mat_v(xe, a+1, 1) = x[a];  
    Mat_v(xe, a+1, 2) = y[a];  
    Mat_v(xe, a+1, 3) = z[a];  
  }
  FEMLIB_set_element(fe, xe, e);  
  
  dealoc1(nod);
  
  dealoc1(x);
  dealoc1(y);
  dealoc1(z);
  Matrix_cleanup(xe);       
}

void FEMLIB_set_element(FEMLIB *fe, Matrix(double) x, int eid)
{
  Matrix_AeqB(fe->node_coord,1.0,x);	
  fe->curt_elem_id = eid;
  for(int a = 1; a<=fe->nne; a++)
  {
  	Mat_v(fe->temp_v.x, a, 1) = Mat_v(x, a, 1);
  	Mat_v(fe->temp_v.y, a, 1) = Mat_v(x, a, 2);
  	Mat_v(fe->temp_v.z, a, 1) = Mat_v(x, a, 3);  	  	
  }  
}

void FEMLIB_elem_shape_function(FEMLIB *fe, long ip, int nne, Matrix(double) N)
{
  double ksi_, eta_, zet_, wt;	
               
  int itg_id_ip_1 = Mat_v(fe->itg_ids, ip, 1);
  int itg_id_ip_2 = Mat_v(fe->itg_ids, ip, 2);
  int itg_id_ip_3 = Mat_v(fe->itg_ids, ip, 3);    

  if(fe->elem_type == HEXAHEDRAL)
  {/* hexahedron */
    ksi_ = Mat_v(fe->ksi, itg_id_ip_1, 1);
    eta_ = Mat_v(fe->eta, itg_id_ip_2, 1);
    zet_ = Mat_v(fe->zet, itg_id_ip_3, 1);
           
    wt = Mat_v(fe->weights, itg_id_ip_1, 1)*Mat_v(fe->weights, itg_id_ip_2, 1)*Mat_v(fe->weights, itg_id_ip_3, 1);
  } 
  else 
  { /* tetrahedron type */
    ksi_ = Mat_v(fe->ksi, itg_id_ip_3, 1);
    eta_ = Mat_v(fe->eta, itg_id_ip_3, 1);
    zet_ = Mat_v(fe->zet, itg_id_ip_3, 1);
    wt = Mat_v(fe->weights, itg_id_ip_3, 1);
  } 
  
  shape_func(ksi_, eta_, zet_, nne, N.m_pdata);
}  

void FEMLIB_elem_basis_V(FEMLIB *fe, long ip)
{
  double ksi_, eta_, zet_, wt;	
               
  int itg_id_ip_1 = Mat_v(fe->itg_ids, ip, 1);
  int itg_id_ip_2 = Mat_v(fe->itg_ids, ip, 2);
  int itg_id_ip_3 = Mat_v(fe->itg_ids, ip, 3);    

  if(fe->elem_type == HEXAHEDRAL)
  {/* hexahedron */
    ksi_ = Mat_v(fe->ksi, itg_id_ip_1, 1);
    eta_ = Mat_v(fe->eta, itg_id_ip_2, 1);
    zet_ = Mat_v(fe->zet, itg_id_ip_3, 1);
           
    wt = Mat_v(fe->weights, itg_id_ip_1, 1)*Mat_v(fe->weights, itg_id_ip_2, 1)*Mat_v(fe->weights, itg_id_ip_3, 1);
  } 
  else 
  { /* tetrahedron type */
    ksi_ = Mat_v(fe->ksi, itg_id_ip_3, 1);
    eta_ = Mat_v(fe->eta, itg_id_ip_3, 1);
    zet_ = Mat_v(fe->zet, itg_id_ip_3, 1);
    wt = Mat_v(fe->weights, itg_id_ip_3, 1);
  } 
  
  fe->temp_v.ksi_ip = ksi_;
  fe->temp_v.eta_ip = eta_;
  fe->temp_v.zet_ip = zet_;
  fe->temp_v.w_ip   = wt;
  
  shape_func(ksi_, eta_, zet_, fe->nne, fe->N.m_pdata);
    
  fe->detJ = deriv(ksi_, eta_, zet_, fe->nne, 
               fe->temp_v.x.m_pdata,   fe->temp_v.y.m_pdata,   fe->temp_v.z.m_pdata, 
               fe->temp_v.N_x.m_pdata, fe->temp_v.N_y.m_pdata, fe->temp_v.N_z.m_pdata);    
  
  Matrix_init(fe->x_ip, 0.0); 
  for(int a = 1; a<=fe->nne; a++)
  {
  	Mat_v(fe->dN, a, 1) = Mat_v(fe->temp_v.N_x, a, 1);
  	Mat_v(fe->dN, a, 2) = Mat_v(fe->temp_v.N_y, a, 1);
  	Mat_v(fe->dN, a, 3) = Mat_v(fe->temp_v.N_z, a, 1);
  	  	  	
  	Mat_v(fe->x_ip,1,1) += Mat_v(fe->N,a,1)*Mat_v(fe->node_coord,a,1);
  	Mat_v(fe->x_ip,2,1) += Mat_v(fe->N,a,1)*Mat_v(fe->node_coord,a,2);
  	Mat_v(fe->x_ip,3,1) += Mat_v(fe->N,a,1)*Mat_v(fe->node_coord,a,3);  	  	  	
  }
  fe->detJxW = fe->detJ*wt;
} 

void FEMLIB_update_shape_tensor(FEMLIB *fe)
{  
  shape_tensor(fe->nne,fe->nsd,fe->temp_v.N_x.m_pdata,
                               fe->temp_v.N_y.m_pdata,
                               fe->temp_v.N_z.m_pdata,fe->ST_tensor);
  shapeTensor2array(fe->ST,CONST_4(double) fe->ST_tensor,fe->nne);
}

void FEMLIB_update_deformation_gradient(FEMLIB *fe, const int ndofn, double *u, Matrix(double) F)
{
  double **F_mat;  
  F_mat = aloc2(3,3);                 
  def_grad_get(fe->nne,ndofn,CONST_4(double) fe->ST_tensor,u,F_mat);
  mat2array(F.m_pdata,CONST_2(double) F_mat,3,3);
  dealoc2(F_mat,3);  
}

double FEMLIB_elem_volume(FEMLIB *fe)
{  
  double Volume;
  if(fe->nne == 4){ /* linear tet */
    Volume = Tetra_V(fe->temp_v.x.m_pdata,fe->temp_v.y.m_pdata,fe->temp_v.z.m_pdata);
  } else if(fe->nne == 10){ /* quadradic tet */
    Volume = Tetra_qv_V(10,3,fe->temp_v.x.m_pdata,fe->temp_v.y.m_pdata,fe->temp_v.z.m_pdata);
  } else if(fe->nne == 8){ /* trilinear hex */
    Volume = Hexa_V(fe->temp_v.x.m_pdata,fe->temp_v.y.m_pdata,fe->temp_v.z.m_pdata);
  } else Volume = 0.0; 
  return Volume;
}

void FEMLIB_destruct(FEMLIB *fe)
{
  FEMLIB_destruct_variable(&(fe->temp_v));

  Matrix_cleanup(fe->ksi);     
  Matrix_cleanup(fe->eta);     
  Matrix_cleanup(fe->zet);     
  Matrix_cleanup(fe->weights);                  
  Matrix_cleanup(fe->N);       
  Matrix_cleanup(fe->dN);      
  Matrix_cleanup(fe->x_ip);    
  Matrix_cleanup(fe->node_coord);
  Matrix_cleanup(fe->itg_ids);
  Matrix_cleanup(fe->node_id);
  
  dealoc4(fe->ST_tensor,3,3,fe->nsd);
  free(fe->ST);
}
