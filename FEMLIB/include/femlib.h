#ifndef FEMLIB_FEMLIB_H
#define FEMLIB_FEMLIB_H

#include "data_structure.h"
#include "element.h"
#include "node.h"

class TEMP_VARIABLES
{
 public:
  gcm::Matrix<double> N_x, N_y, N_z, x, y, z;
  double ksi_ip, eta_ip, zet_ip, w_ip;
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

  gcm::Matrix<double> N;
  gcm::Matrix<double> dN;

  gcm::Matrix<double> ksi, eta, zet, weights;
  gcm::Matrix<int> itg_ids;
  gcm::Matrix<long> node_id;

  long intg_order;

  gcm::Matrix<double> node_coord;

  double detJ, detJxW, normal;
  gcm::Matrix<double> x_ip;
  long elem_type, nint, nne, nsd, curt_elem_id, curt_itg_id;
  TEMP_VARIABLES temp_v;
  double ****ST_tensor;
  double *ST;

  FEMLIB(){};
  FEMLIB(int e_type,
         int i_order,
         int nne)
  {
    initialization(e_type, i_order, nne);
  };
  FEMLIB(int e,
         const Element *elem,
         const Node *node,
         int i_order,
         int is_total,
         bool add_bubble = false)
  {
    initialization(e,elem,node,i_order,is_total, add_bubble);
  };

  ~FEMLIB();

  long determine_integration_type(int e_type,
                                  int i_order);

  void initialization(int e_type,
                      int i_order,
                      int nne);

  /// set FEM libray by element
  void initialization(int e,
                      const Element *elem,
                      const Node *node,
                      int i_order,
                      int is_total,
                      bool add_bubble = false);

  void elem_shape_function(long ip,
                           int nne,
                           double *N);
  void elem_basis_V(long ip);
  void update_shape_tensor(void);
  void update_deformation_gradient(const int ndofn,
                                   double *u,
                                   double *F);
  double elem_volume(void);
};

#endif // #define FEMLIB_FEMLIB_H
