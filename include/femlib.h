#ifndef _PGFem3D_FEMLIB_H_
#define _PGFem3D_FEMLIB_H_

#include "data_structure_c.h"
#include "element.h"
#include "node.h"

#ifndef _Matrix_double
Define_Matrix(double);
#define _Matrix_double 1
#endif

#ifndef _Matrix_int
Define_Matrix(int);
#define _Matrix_int 1
#endif

#ifndef _Matrix_long
Define_Matrix(long);
#define _Matrix_long 1
#endif

typedef struct {
  Matrix(double) N_x, N_y, N_z, x, y, z;
  double ksi_ip, eta_ip, zet_ip, w_ip;
} TEMP_VARIABLES;

 
typedef struct {
  Matrix(double) N;
  Matrix(double) dN;
  
	Matrix(double) ksi, eta, zet, weights;
	Matrix(int) itg_ids;
	Matrix(long) node_id;

	long intg_order;
		  
  Matrix(double) node_coord;
  
  double detJ, detJxW, normal;
  Matrix(double) x_ip;
  long elem_type, nint, nne, nsd, curt_elem_id, curt_itg_id;  
  TEMP_VARIABLES temp_v;
  double ****ST_tensor;
  double *ST;
         
 } FEMLIB;

long FEMLIB_determine_integration_type(int e_type, int i_order);

//  FEMLIB(long elem_type, long intg_order = 1)
//  {
//  	integration_points(elem_type, intg_order);
//  };

void FEMLIB_set_variable_size(TEMP_VARIABLES *v, int nne);
void FEMLIB_initialization(FEMLIB *fe, int e_type, int i_order, int nne);
void FEMLIB_initialization_by_elem(FEMLIB *fe, int e, const ELEMENT *elem, const NODE *node, int i_order);
void FEMLIB_set_element(FEMLIB *fe, Matrix(double) x, int eid);
void FEMLIB_elem_shape_function(FEMLIB *fe, long ip, int nne, Matrix(double) N);
void FEMLIB_elem_basis_V(FEMLIB *fe, long ip);
void FEMLIB_update_shape_tensor(FEMLIB *fe);
void FEMLIB_update_deformation_gradient(FEMLIB *fe, const int ndofn, double *u, Matrix(double) F);
double FEMLIB_elem_volume(FEMLIB *fe);
void FEMLIB_destruct(FEMLIB *fe);

/*  
void FEMLIB_integration_points(FEMLIB *fe, long e_type, long i_order)
{ 
	nsd = 3;
  elem_type = e_type;
  intg_order = i_order;

  long intg_type = FEMLIB_determine_integration_type(e_type, i_order);

  int_point(intg_type, &nint);
  
  fe->ksi.redim(nint);
  eta.redim(nint);
  zet.redim(nint);
  weights.redim(nint);  
  
  long npt_x, npt_y, npt_z;
  integrate(e_type, &npt_x, &npt_y, &npt_z,
          ksi.m_pdata, eta.m_pdata, zet.m_pdata,
          weights.m_pdata);

  long cnt = 0;
  itg_ids.redim(npt_x*npt_y*npt_z, nsd);
  for(int a=1; a<=npt_x; a++)
  {
  	for(int b=1; b<=npt_y; b++)
    {
    	for(int c=1; c<=npt_z; c++)
    	{
    		cnt++;
        itg_ids(cnt, 1) = a; itg_ids(cnt, 2) = b;  itg_ids(cnt, 3) = c;    		
       }
     }
   }           
}

void FEMLIB::initialize_integration(long nne_, long nsd_)
{
	nne = nne_;
	nsd = nsd_;
  N.redim(nne);	
  dN.redim(nne, nsd);
  temp_v.set_variable_size(nne);
  x_ip.redim(nsd);
}

void FEMLIB::set_element(Matrix<double> &x, long eid)
{
  node_coord = x;	
  curt_elem_id = eid;
  for(int a = 1; a<=nne; a++)
  {
  	temp_v.x(a) = x(a, 1);
  	temp_v.y(a) = x(a, 2);
  	temp_v.z(a) = x(a, 3);  	  	
  }  
}

void FEMLIB::elem_basis_V(long ip)
{
  
  double ksi_, eta_, zet_, wt;	
  int err;
               
        
  if(elem_type == HEXAHEDRAL)
  {// hexahedron 
    ksi_ = ksi(itg_ids(ip, 1));
    eta_ = eta(itg_ids(ip, 2));
    zet_ = zet(itg_ids(ip, 3));
    wt = weights(itg_ids(ip, 1))*weights(itg_ids(ip, 2))*weights(itg_ids(ip, 3));
  } 
  else 
  { // tetrahedron type
    ksi_ = ksi(itg_ids(ip, 3));
    eta_ = eta(itg_ids(ip, 3));
    zet_ = zet(itg_ids(ip, 3));
    wt = weights(itg_ids(ip, 3));
  }
    	
  shape_func(ksi_, eta_, zet_, nne, N.m_pdata);
  detJ = deriv(ksi_, eta_, zet_, nne, temp_v.x.m_pdata, temp_v.y.m_pdata, temp_v.z.m_pdata, 
               temp_v.N_x.m_pdata, temp_v.N_y.m_pdata, temp_v.N_z.m_pdata);

  x_ip.redim(nsd);
  for(int a = 1; a<=nne; a++)
  {
  	dN(a, 1) = temp_v.N_x(a);
  	dN(a, 2) = temp_v.N_y(a);
  	dN(a, 3) = temp_v.N_z(a);
  	x_ip(1) += N(a)*node_coord(a, 1);
  	x_ip(2) += N(a)*node_coord(a, 2);
  	x_ip(3) += N(a)*node_coord(a, 3);  	  	
  }

  detJxW = detJ*wt;
}
	*/
#endif
