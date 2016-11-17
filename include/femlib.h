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


struct FEMLIB;
#ifndef TYPE_FEMLIB
#define TYPE_FEMLIB
typedef struct FEMLIB FEMLIB;
#endif
 
struct FEMLIB {
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
         
 };

long FEMLIB_determine_integration_type(int e_type, int i_order);

/// set FEM libray by element
///
/// \param[out] fe self of FEMLIB
/// \param[in] e element id
/// \param[in] elem list of ELEMENT
/// \param[in] node list of NODE
/// \param[in] i_order integration order, 0: linear, 1: quadratic, and 2: higher
/// \param[in] is_total if 1: total Lagrangian, 0: updated Lagrangian
/// \return void
void FEMLIB_initialization_by_elem(FEMLIB *fe, int e, const ELEMENT *elem, const NODE *node, int i_order, int is_total);

/// set FEM libray by element with bubble
///
/// Set the FEM library same as FEMLIB_initialization_by_elem, but if bubble is added, it computes
/// center node in the element (MINI and MINI_3F)
///
/// \param[out] fe self of FEMLIB
/// \param[in] e element id
/// \param[in] elem list of ELEMENT
/// \param[in] node list of NODE
/// \param[in] i_order integration order, 0: linear, 1: quadratic, and 2: higher
/// \param[in] is_total if 1: total Lagrangian, 0: updated Lagrangian
/// \return void
void FEMLIB_initialization_by_elem_w_bubble(FEMLIB *fe, 
                                            int e, 
                                            const ELEMENT *elem, 
                                            const NODE *node, 
                                            int i_order, 
                                            int is_total);
                                            
void FEMLIB_elem_shape_function(FEMLIB *fe, long ip, int nne, Matrix(double) *N);
void FEMLIB_elem_basis_V(FEMLIB *fe, long ip);
void FEMLIB_update_shape_tensor(FEMLIB *fe);
void FEMLIB_update_deformation_gradient(FEMLIB *fe, const int ndofn, double *u, Matrix(double) *F);
double FEMLIB_elem_volume(FEMLIB *fe);
void FEMLIB_destruct(FEMLIB *fe);
#endif
