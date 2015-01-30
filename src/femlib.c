#include "femlib.h"
#include "elem3d.h"

#include "data_structure_c.h"
#include "utils.h"

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
      
  FEMLIB_set_variable_size(&(fe->temp_v), nne);      
    
  long npt_x, npt_y, npt_z;
  integrate(intg_type, &npt_x, &npt_y, &npt_z,
          fe->ksi.m_pdata, fe->eta.m_pdata, fe->zet.m_pdata,
          fe->weights.m_pdata);

  Matrix_construct(int, fe->itg_ids); Matrix_redim(fe->itg_ids, npt_x*npt_y*npt_z, nsd);

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
}


void FEMLIB_set_element(FEMLIB *fe, Matrix(double) x, int eid)
{
  Matrix_AeqB(fe->node_coord,1,x);	
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
}
