#ifndef FEMLIB_FEMLIB_H
#define FEMLIB_FEMLIB_H

#include "data_structure.h"
#include "element.h"
#include "node.h"

static const constexpr int LinearElement    = 0;
static const constexpr int QuadraticElement = 1;

static const constexpr int FemDim0D = 0;
static const constexpr int FemDim1D = 1;
static const constexpr int FemDim2D = 2;
static const constexpr int FemDim3D = 3;

static const constexpr int Tet2Tri[4][3] = {{0,2,1},{0,1,3},{1,2,3},{0,3,2}};
static const constexpr int QTetQTri[4][6] = {{0,2,1,6,5,4},{0,1,3,4,8,7},{1,2,3,5,9,8},{0,3,2,7,9,6}};
static const constexpr int Hex2Quad[6][4] = {{3,2,1,0},{7,4,5,6},{0,1,5,4},{2,6,5,1},{3,7,6,2},{3,0,4,7}};

static const constexpr int element_type(const int nne,
                                        const int dim){
  return dim*10 + nne;
}

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
  gcm::Matrix<double> u0;

  long intg_order;

  gcm::Matrix<double> node_coord;

  double detJ, detJxW, normal;
  gcm::Matrix<double> x_ip;
  long elem_type, nint, nne, nsd, curt_elem_id, curt_itg_id;
  TEMP_VARIABLES temp_v;
  double ****ST_tensor;
  double *ST;

  FEMLIB(){};
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

  long determine_integration_type(const int e_type,
                                  const int e_order,
                                  const int i_order);

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
  void elem_basis_V(long ip);
  void update_shape_tensor(void);
  void update_deformation_gradient(const int ndofn,
                                   double *u,
                                   double *F);
  void update_deformation_gradient(const int ndofn, 
                                   double *u, 
                                   double *F,
                                   double *pF0I);                                   
  double elem_volume(void);
};

class FemLib2D : public FEMLIB{
public:
  const FEMLIB *fe3D;

  FemLib2D(){};
  FemLib2D(const FEMLIB *fe,
           const int face_id,
           const int i_order){
    initialization(fe, face_id, i_order);
  };

  void initialization(const FEMLIB *fe,
                      const int face_id,
                      const int i_order);

};

#endif // #define FEMLIB_FEMLIB_H
