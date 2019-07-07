#ifndef FEMLIB_FEMLIB_H
#define FEMLIB_FEMLIB_H

#include "data_structure.h"
#include "element.h"
#include "node.h"

class TEMP_VARIABLES
{
 public:
  gcm::Matrix<double> x, y, z;
  TEMP_VARIABLES(){};
  TEMP_VARIABLES(int nne_t, int nne)
  {
    set_variable_size(nne_t, nne);
  };
  ~TEMP_VARIABLES(){};
  void set_variable_size(int nne_t, int nne);
};

class FEMLIB
{
 public:

  gcm::Matrix<double> N;          //! shape function
  gcm::Matrix<double> dN;         //! dNdX

  gcm::Matrix<double> ksi, eta, zet, weights; //! isoparametric coordinates and weights
                                              //! for the quadrature rules
  gcm::Matrix<int> itg_ids;       //! list of integration IDs of ksi, eta, and zet
  gcm::Matrix<long> node_id;      //! list of nodes in element under FEM integration
  gcm::Matrix<double> u0;         //! nodal variable values

  int intg_order;                 //! order of FEM integration

  gcm::Matrix<double> node_coord; //! nodal coordinates, [node id, xyz id = 0,1,2]

  double detJ;                    //! determinant of Jaccobian, det(J)
  double detJxW;                  //! det(J)*Weight(ip)
  gcm::Matrix<double> x_ip;       //! coordiate at integration points
    
  int elem_type, nint, nne, nsd, curt_elem_id, curt_itg_id;
  int bnd_elem_no;
  
  TEMP_VARIABLES temp_v;
  double *ST;

  FEMLIB(){
    ST = NULL;
  };
  FEMLIB(const int nne,
         const int nsd,
         const int i_order){
    initialization(nne, nsd, i_order);
  };
  FEMLIB(int e,
         const Element *elem,
         const Node *node,
         int i_order,
         int is_total,
         bool add_bubble = false)
  {
    initialization(e,elem,node,i_order,is_total, NULL, add_bubble);
  };

  ~FEMLIB();

  int determine_integration_type(const int e_type,
                                  const int e_order,
                                  const int i_order);

  int number_of_boundary_elements(const int e_type);
  int number_of_boundary_elements(const int nne,
                                  const int nsd);

  void initialization(const int nne,
                      const int nsd,
                      const int i_order);

  /// set FEM libray by element
  void initialization(int e,
                      const Element *elem,
                      const Node *node,
                      int i_order,
                      int is_total,
                      double *pF0I,
                      bool add_bubble = false);

  void elem_shape_function(long ip,
                           int nne,
                           double *N);
  double compute_integration_weight(const int ip);
  void elem_basis_V(long ip);
  void update_shape_tensor(void);
  void update_deformation_gradient(const int ndofn,
                                           double *u,
                                           double *F);
  void update_deformation_gradient(const int ndofn, 
                                   const double *u, 
                                   double *F,
                                   double *pF0I);

  double elem_volume(void);
};

class FemLibBoundary : public FEMLIB{
public:
  FEMLIB *feVol;                  //! Volumetric FEMLIB
  const int *Volume2Boundary;     //! Volumetric to boundary map of node IDs
  const int *kez_map;             //! Volumetric to boundary map of
                                  //! isoparametric coordinates
  int face_id;                    //! face id of element
  double normal[3];               //! normal vector
  int nne_bnd;                    //! number of nodes of boundary element

  FemLibBoundary(){};
  FemLibBoundary(FEMLIB *fe,
                 const int face_id,
                 const int i_order){
    initialization(fe, face_id, i_order);
  };

  void initialization(FEMLIB *fe,
                      const int face_id,
                      const int i_order);

  void set_volume_to_boundary_map(void);
  double compute_integration_weight(const int ip);
  void elem_basis_S(const int ip);
};

#endif // #define FEMLIB_FEMLIB_H
