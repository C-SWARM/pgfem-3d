/// Define finite element integration rules
///
/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
///
/// [1] University of Notre Dame, Notre Dame, IN

// In this file, Vol and BND stands for Volume and boundary element, respectively.
//               nsd is number of spatial dimesions.

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "femlib.h"
#include "allocation.h"
#include "cast_macros.h"
#include "data_structure.h"
#include "def_grad.h"
#include "elem3d.h"
#include "tensors.h"
#include "utils.h"

// element order
static const constexpr int LinearElement    = 0;
static const constexpr int QuadraticElement = 1;

// finite element demesion
static const constexpr int FemDim0D = 0;
static const constexpr int FemDim1D = 1;
static const constexpr int FemDim2D = 2;
static const constexpr int FemDim3D = 3;

/// compute element type by number of nodes per element and spatial dimension
/// spatial dimension increase by 10 in setting element type 
///
/// \param[in] 
static const constexpr int element_type(const int nne,
                                        const int dim){
  return dim*10 + nne;
}


// define element types
static const constexpr int          LINE = element_type( 2, FemDim1D);
static const constexpr int         QLINE = element_type( 3, FemDim1D);
static const constexpr int      TRIANGLE = element_type( 3, FemDim2D);
static const constexpr int     QTRIANGLE = element_type( 6, FemDim2D);
static const constexpr int QUADRILATERAL = element_type( 4, FemDim2D);
static const constexpr int   TETRAHEDRON = element_type( 4, FemDim3D);
static const constexpr int  QTETRAHEDRON = element_type(10, FemDim3D);
static const constexpr int    HEXAHEDRAL = element_type( 8, FemDim3D);
static const constexpr int        WEDGE  = element_type( 6, FemDim3D);

// define mapping from volume to boundary elements by face id
// Volume2Boundary[face_id][nne_BND]
static const constexpr int Tet2Tri[4][3] = {{0,2,1},{0,1,3},{1,2,3},{0,3,2}};
static const constexpr int QTet2QTri[4][6] = {{0,2,1,6,5,4},{0,1,3,4,8,7},{1,2,3,5,9,8},{0,3,2,7,9,6}};
static const constexpr int Hex2Quad[6][4] = {{3,2,1,0},{7,4,5,6},{0,1,5,4},{2,6,5,1},{3,7,6,2},{3,0,4,7}};

// define mapping from nsd_Vol-1 integration rules to nsd_Vol for boundary integration in Volume element
static const constexpr int kez_map_Hex[6][4] = {{1, 0, 2, -1},  // ksi2D(0) -> eta3D(1), eta2D(1) -> ksi3D(0), zet3D(2) = -1
                                                {0, 1, 2,  1},  // ksi2D(0) -> ksi3D(0), eta2D(1) -> eta3D(1), zet3D(2) =  1
                                                {1, 2, 0,  1},  // ksi2D(0) -> eta3D(1), eta2D(1) -> zet3D(2), ksi3D(0) =  1
                                                {2, 0, 1,  1},  // ksi2D(0) -> zet3D(2), eta2D(1) -> ksi3D(0), eta2D(1) =  1
                                                {2, 1, 0, -1},  // ksi2D(0) -> zet3D(2), eta2D(1) -> eta3D(1), ksi3D(0) = -1
                                                {0, 2, 1, -1}}; // ksi2D(0) -> ksi3D(0), eta2D(1) -> zet3D(2), eta2D(1) = -1

static const constexpr int kez_map_Tet[4][4] = {{0,1, 2, 0},  // ksi2D(0) -> ksi3D(0), eta2D(1) -> eta3D(1), zet3D(2) = 0
                                                {0,2, 1, 0},  // ksi2D(0) -> ksi3D(0), eta2D(1) -> zet3D(2), eta2D(1) = 0
                                                {0,1, 2, 0},  // ksi2D(0) -> ksi3D(0), eta2D(1) -> eta3D(1), zet3D(2) = 1 - ksi_3D - eta_3D;
                                                {1,2, 0, 0}}; // ksi2D(0) -> eta3D(1), eta2D(1) -> zet3D(2), ksi3D(0) = 0
                                                  
/// 
void print_func_file_line_and_abort(const char *fname,
                                    const char *filename,
                                    const int line){
 PGFEM_printerr("(%s:%s:%d)\n", fname, filename, line);
 pgfem3d::PGFEM_Abort();
}

void FEM_shape_functions(const int element_type, 
                         const double ksi,
                         const double eta,
                         const double zet,
                         double *N){
  switch(element_type){
    case LINE:
      N[0] = 0.5*(1.0 - ksi);
      N[1] = 0.5*(1.0 + ksi);
      break;
    case QLINE:
      N[0] = -0.5*(1.0 - ksi)*ksi;
      N[1] = (1.0 - ksi)*(1.0 + ksi);
      N[2] = 0.5*(1.0 + ksi)*ksi;
      break;
    case TRIANGLE:
      N[0]= 0.25*(1.0 - ksi)*(1.0 - eta);
      N[1]= 0.25*(1.0 + ksi)*(1.0 - eta);
      N[2]= 0.50*(1.0 + eta);
      break;
    case QUADRILATERAL:
      N[0] = 0.25*(1.0 - ksi)*(1.0 - eta);
      N[1] = 0.25*(1.0 + ksi)*(1.0 - eta);
      N[2] = 0.25*(1.0 + ksi)*(1.0 + eta);
      N[3] = 0.25*(1.0 - ksi)*(1.0 + eta);
      break;
    case QTRIANGLE:
      break;
    case TETRAHEDRON:
      N[0] = (1.0 - ksi - eta - zet);
      N[1] = ksi;
      N[2] = eta;
      N[3] = zet;
      break;
    case QTETRAHEDRON:
    {
      double ksi_eta_zet = 1.0 - ksi - eta - zet;
      N[0] = (2.0*ksi_eta_zet - 1.0)*ksi_eta_zet;
      N[1] = (2.0*ksi-1.0)*ksi;
      N[2] = (2.0*eta-1.0)*eta;
      N[3] = (2.0*zet-1.0)*zet;
      N[4] = 4.0*ksi_eta_zet*ksi;
      N[5] = 4.0*ksi*eta;
      N[6] = 4.0*ksi_eta_zet*eta;
      N[7] = 4.0*ksi_eta_zet*zet;
      N[8] = 4.0*ksi*zet;
      N[9] = 4.0*eta*zet;
      break;
    }
    case HEXAHEDRAL:
      N[0] = 1.0/8.0*(1.0 - ksi)*(1.0 - eta)*(1.0 - zet);
      N[1] = 1.0/8.0*(1.0 + ksi)*(1.0 - eta)*(1.0 - zet);
      N[2] = 1.0/8.0*(1.0 + ksi)*(1.0 + eta)*(1.0 - zet);
      N[3] = 1.0/8.0*(1.0 - ksi)*(1.0 + eta)*(1.0 - zet);  
      N[4] = 1.0/8.0*(1.0 - ksi)*(1.0 - eta)*(1.0 + zet);
      N[5] = 1.0/8.0*(1.0 + ksi)*(1.0 - eta)*(1.0 + zet);
      N[6] = 1.0/8.0*(1.0 + ksi)*(1.0 + eta)*(1.0 + zet);
      N[7] = 1.0/8.0*(1.0 - ksi)*(1.0 + eta)*(1.0 + zet);
      break;
    case WEDGE:
      N[0] = 0.5*(1.0 - zet)*ksi;
      N[1] = 0.5*(1.0 - zet)*eta;
      N[2] = 0.5*(1.0 - zet)*(1.0-ksi-eta);
      N[3] = 0.5*(1.0 + zet)*ksi;
      N[4] = 0.5*(1.0 + zet)*eta;
      N[5] = 0.5*(1.0 + zet)*(1.0-ksi-eta);
    default:
      PGFEM_printerr("ERROR: Sahpe function for element %d is not implemented.\n", element_type);
      print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
  }
}

static const constexpr double one_over_sqrt_3 = 0.57735026918962584; // 1.0/sqrt(3.0);

typedef struct {
  const double weights[1]  = {2.0};
  const double gk[1] = {0.0};
  const double ge[1] = {0.0};
  const double gz[1] = {0.0};
  int gpno = 1;
} GaussIntegrationPoints1;

typedef struct {
  const double weights[2]  = {1.0,1.0};
  const double gk[2] = {-one_over_sqrt_3, one_over_sqrt_3};
  const double ge[2] = {-one_over_sqrt_3, one_over_sqrt_3};
  const double gz[2] = {-one_over_sqrt_3, one_over_sqrt_3};
  int gpno = 2;
} GaussIntegrationPoints2;

typedef struct {
  const double weights[3]  = {0.555555555555555, 0.888888888888888, 0.555555555555555};
  const double gk[3] = {-0.774596669241483, 0.0, 0.774596669241483};
  const double ge[3] = {-0.774596669241483, 0.0, 0.774596669241483};
  const double gz[3] = {-0.774596669241483, 0.0, 0.774596669241483};
  int gpno = 3;
} GaussIntegrationPoints3;

typedef struct {
    const double weights[1]  = {1.0/6.0};
    const double gk[1] = {1.0/4.0};
    const double ge[1] = {1.0/4.0};
    const double gz[1] = {1.0/4.0};
    
    int gpno = 1;
} TetIntegrationPoints1;

typedef struct {
    const double weights[4]  = {(1.0/4.0)/6.0, (1.0/4.0)/6.0, (1.0/4.0)/6.0, (1.0/4.0)/6.0};
    const double gk[4] = {0.58541020, 0.13819660, 0.13819660, 0.13819660};
    const double ge[4] = {0.13819660, 0.58541020, 0.13819660, 0.13819660};
    const double gz[4] = {0.13819660, 0.13819660, 0.58541020, 0.13819660};
    
    int gpno = 4;
} TetIntegrationPoints4;

typedef struct {
    const double weights[5]  = {(-4.0/5.0)/6.0, (9.0/20.0)/6.0, (9.0/20.0)/6.0, (9.0/20.0)/6.0};
    const double gk[5] = {1.0/4.0, 1.0/2.0, 1.0/6.0, 1.0/6.0, 1.0/6.0};
    const double ge[5] = {1.0/4.0, 1.0/6.0, 1.0/2.0, 1.0/6.0, 1.0/6.0};
    const double gz[5] = {1.0/4.0, 1.0/6.0, 1.0/6.0, 1.0/2.0, 1.0/6.0};
    
    int gpno = 5;
    
} TetIntegrationPoints5;

typedef struct {
    const double weights[11]  = {-0.01315555556,  0.007622222222, 0.007622222222, 0.007622222222,
                                  0.007622222222, 0.02488888889,  0.02488888889,  0.02488888889,
                                  0.02488888889,  0.02488888889,  0.02488888889};
    const double gk[11] = {1.0/4.0,       0.07142857143, 0.7857142857,  0.7857142857,
                           0.7857142857,  0.3994035762,  0.3994035762,  0.3994035762,
                           0.3994035762,  0.3994035762,  0.3994035762};
    const double ge[11] = {1.0/4.0,       0.07142857143, 0.07142857143, 0.07142857143,
                           0.07142857143, 0.3994035762,  0.1005964238,  0.1005964238,
                           0.3994035762,  0.3994035762,  0.1005964238};
    const double gz[11] = {1.0/4.0,       0.07142857143, 0.07142857143, 0.07142857143,
                           0.07142857143, 0.1005964238,  0.1005964238,  0.1005964238,
                           0.1005964238, 0.1005964238, 0.1005964238};
    
    int gpno = 11;
    
} TetIntegrationPoints11;

typedef struct {
    const double weights[1]  = {0.5};
    const double gk[1] = {1.0/3.0};
    const double ge[1] = {1.0/3.0};
    const double gz[1] = {};
    
    int gpno = 1;
} TriIntegrationPoints1;

typedef struct {
    const double weights[3]  = {1.0/6.0, 1.0/6.0, 1.0/6.0};
    const double gk[3] = {2.0/3.0, 1.0/6.0, 1.0/6.0};
    const double ge[3] = {1.0/6.0, 1.0/6.0, 2.0/3.0};
    const double gz[3] = {};
        
    int gpno = 3;
} TriIntegrationPoints3;

typedef struct {
    const double weights[4]  = {-9.0/32.0, 25.0/96.0, 25.0/96.0, 25.0/96.0};
    const double gk[4] = {1.0/3.0, 0.6, 0.2, 0.2};
    const double ge[4] = {1.0/3.0, 0.2, 0.2, 0.6};
    const double gz[4] = {};
    
    int gpno = 4;
} TriIntegrationPoints4;

class QuadratureRule{
  public:
    QuadratureRule(){};

    template <class GP>
    void get_gauss_integration_point_ids_3D(const GP &gp,
                                            Matrix<int> &itg_ids){
      const int nsd = 3;
      itg_ids.initialization(gp.gpno*gp.gpno*gp.gpno, nsd);
          
      long cnt = 0;
      for(int ia=0; ia<gp.gpno; ++ia){
        for(int ib=0; ib<gp.gpno; ++ib){
          for(int ic=0; ic<gp.gpno; ++ic){
            itg_ids(cnt, 0) = ia;
            itg_ids(cnt, 1) = ib;
            itg_ids(cnt, 2) = ic;
            ++cnt;
          }
        }
      }
    }
    
    template <class GP>
    void get_gauss_integration_point_ids_2D(const GP &gp,
                                            Matrix<int> &itg_ids){
      const int nsd = 2;
      itg_ids.initialization(gp.gpno*gp.gpno, nsd);
          
      long cnt = 0;
      for(int ia=0; ia<gp.gpno; ++ia){
        for(int ib=0; ib<gp.gpno; ++ib){
          itg_ids(cnt, 0) = ia;
          itg_ids(cnt, 1) = ib;
          ++cnt;
        }
      }
    }
    
    template <class GP>
    void get_integration_point_ids_listed(const GP &gp,
                                          Matrix<int> &itg_ids,
                                          const int nsd){
      itg_ids.initialization(gp.gpno, nsd);
          
      for(int ia=0; ia<gp.gpno; ++ia)
        for(int ib=0; ib<nsd; ++ib)
          itg_ids(ia, ib) = ia;
    }
    
    template <class GP>
    void get_integration_point_ids_listed_Boundary(const GP &gp,
                                                   Matrix<int> &itg_ids,
                                                   const int volume_nsd,
                                                   const int face_id){
      itg_ids.initialization(gp.gpno, volume_nsd);
          
      for(int ia=0; ia<gp.gpno; ++ia)
        for(int ib=0; ib<volume_nsd; ++ib)
          itg_ids(ia, ib) = ia;
    }

    template <class GP, class T>
    void get_gauss_integration_point_ids_3D_Boundary(const GP &gp,
                                                     Matrix<int> &itg_ids, 
                                                     const int face_id,
                                                     const T kez_map){
      const int nsd = 3;
      itg_ids.initialization(gp.gpno*gp.gpno, nsd);
      
      long cnt = 0;
      for(int ia=0; ia<gp.gpno; ++ia){
        for(int ib=0; ib<gp.gpno; ++ib){
            itg_ids(cnt, kez_map[face_id][0]) = ia;
            itg_ids(cnt, kez_map[face_id][1]) = ib;
            itg_ids(cnt, kez_map[face_id][2]) = 0;
            ++cnt;
        }
      }
    }
    
    template <class GP, class T>
    void get_gauss_integration_point_ids_2D_Boundary(const GP &gp,
                                                     Matrix<int> &itg_ids, 
                                                     const int face_id,
                                                     const T kez_map){
      const int nsd = 2;
      itg_ids.initialization(gp.gpno, nsd);
      
      for(int ia=0; ia<gp.gpno; ++ia){
        itg_ids(ia, kez_map[face_id][0]) = ia;
        itg_ids(ia, kez_map[face_id][1]) = 0;
      }
    }
};

template <class GP> class QuadratureRuleVolume : public QuadratureRule{
  public:
    GP gp;
    QuadratureRuleVolume(){};
    void get_integration_points(Matrix<double> &ksi, Matrix<double> &eta, Matrix<double> &zet, Matrix<double> &weights){
      for(int ia=0; ia<gp.gpno; ++ia){
        ksi(ia) = gp.gk[ia];
        eta(ia) = gp.ge[ia];
        zet(ia) = gp.gz[ia];
        weights(ia) = gp.weights[ia];
      }
    }
    void get_integration_point_ids(Matrix<int> &itg_ids, const int nsd, bool is_Gaussian = true){
      if(is_Gaussian){
        switch(nsd){
          case 1:
            get_integration_point_ids_listed(gp, itg_ids, nsd);
            return;
          case 2:
            get_gauss_integration_point_ids_2D(gp, itg_ids);
            return;
          case 3:
            get_gauss_integration_point_ids_3D(gp, itg_ids);
            return;
          default:
            PGFEM_printerr("ERROR: number of spatial dimension %d is not supported.\n", nsd);
            print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
        }              
      } 
      else
        get_integration_point_ids_listed(gp, itg_ids, nsd);
    }    
};

template <class GP> class QuadratureRuleBoundary : public QuadratureRule{
  public:
    GP gp;
    QuadratureRuleBoundary(){};
    void set_hex2quad(Matrix<double> &ksi,
                      Matrix<double> &eta,
                      Matrix<double> &zet,
                      Matrix<double> &weights,
                      const int face_id){
      double *kez3D[3] = {ksi.m_pdata, eta.m_pdata, zet.m_pdata};
      
      for(int ia=0; ia<gp.gpno; ++ia){
        kez3D[kez_map_Hex[face_id][0]][ia] = gp.gk[ia];
        kez3D[kez_map_Hex[face_id][1]][ia] = gp.ge[ia];
        kez3D[kez_map_Hex[face_id][2]][ia] = 1.0*kez_map_Hex[face_id][3];
      }
    }
    
    void set_tet2tri(Matrix<double> &ksi,
                     Matrix<double> &eta,
                     Matrix<double> &zet,
                     Matrix<double> &weights,
                     const int face_id){
      double *kez3D[3] = {ksi.m_pdata, eta.m_pdata, zet.m_pdata};
      
      for(int ia=0; ia<gp.gpno; ++ia){
        kez3D[kez_map_Tet[face_id][0]][ia] = gp.gk[ia];
        kez3D[kez_map_Tet[face_id][1]][ia] = gp.ge[ia];
        if(face_id==2)
          zet(ia) = 1.0 - gp.gk[ia] - gp.ge[ia];
      }
    }
         
    void get_integration_points(Matrix<double> &ksi,
                                Matrix<double> &eta,
                                Matrix<double> &zet,
                                Matrix<double> &weights,
                                const int face_id,
                                const int e_type_volume){
      if(e_type_volume == HEXAHEDRAL)
        set_hex2quad(ksi, eta, zet, weights, face_id);
        
      if(e_type_volume == TETRAHEDRON || e_type_volume == QTETRAHEDRON)
        set_tet2tri(ksi, eta, zet, weights, face_id);        

    }
    
    void get_integration_point_ids(Matrix<int> &itg_ids, 
                                   const int volume_nsd, 
                                   const int face_id, 
                                   bool is_Gaussian = true){     
      if(is_Gaussian){
        switch(volume_nsd){
          case 3:
            get_gauss_integration_point_ids_3D_Boundary(gp, itg_ids, face_id, kez_map_Hex);
            return;
          default:
            PGFEM_printerr("ERROR: number of spatial dimension %d is not implemented.\n", volume_nsd);
            print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
        }              
      } else{
        switch(volume_nsd){
          case 3:
            get_integration_point_ids_listed_Boundary(gp, itg_ids, volume_nsd, face_id);
            return;
          default:
            PGFEM_printerr("ERROR: number of spatial dimension %d is not implemented.\n", volume_nsd);
            print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
        }
      }
    }
};


void
TEMP_VARIABLES::set_variable_size(int nne_t, int nne)
{
  this->x.initialization(nne_t, 1);
  this->y.initialization(nne_t, 1);
  this->z.initialization(nne_t, 1);
  this->N_x.initialization(nne, 1);
  this->N_y.initialization(nne, 1);
  this->N_z.initialization(nne, 1);
}

long number_of_integration_points_line(const int order){
  return order + 1;
}

long number_of_integration_points_tri(const int order){
  switch(order){
    case 0:
    case 1:
      return 1;
    case 2:
      return 3;
    case 3:
      return 4;
    default:
      PGFEM_printerr("ERROR: number of integration points: order %d is not implemented for TRIANGLE element.\n", order);
      print_func_file_line_and_abort(__func__,__FILE__,__LINE__); 
  }
  return 0;
}

long number_of_integration_points_quad(const int order){
  switch(order){
    case 0:
      return 1;
    case 1:
      return 4;
    case 2:
      return 9;
    default:
      PGFEM_printerr("ERROR: number of integration points: order %d is not implemented for QUADRILATERAL element.\n", order);
      print_func_file_line_and_abort(__func__,__FILE__,__LINE__);       
  }
  return 0;
}

long number_of_integration_points_tet(const int order){
  switch(order){
    case 0:
      return 1;
    case 1:
      return 4;
    case 2:
      return 5;
    case 3:
      return 11;
    case 4:
      return 24;
    default:
      PGFEM_printerr("ERROR: number of integration points: "
                     " order %d is not implemented for TETRAHEDRON element.\n", order);
      print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
  }
  return 0;
}

long number_of_integration_points_hex(const int order){
  switch(order){
    case 0:
      return 1;
    case 1:
      return 8; // 4th order accuracy
    case 2:
      return 14; // 6th order accuracy
    case 3:
      return 27; // 
    default:
      PGFEM_printerr("ERROR: number of integration points: "
                     " order %d is not implemented for HEXAHEDRAL element.\n", order);
      print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
  }
  return 0;
}

long
FEMLIB::determine_integration_type(const int e_type,
                                   const int e_order,
                                   const int i_order)
{
  switch(e_type){
    case LINE:
      return number_of_integration_points_line(i_order);
    case TRIANGLE:
      return number_of_integration_points_tri(i_order);
    case QUADRILATERAL:
      return number_of_integration_points_quad(i_order);
    case TETRAHEDRON:
    case QTETRAHEDRON:
      return number_of_integration_points_tet(i_order);
    case HEXAHEDRAL:
      return number_of_integration_points_hex(i_order);
    default:
      PGFEM_printerr("ERROR: element type %d is not implemented.\n", e_type);
      print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
  }
  return 0;
}

void
FEMLIB::initialization(const int nne,
                       const int nsd,
                       const int i_order)
{
  this->nsd = nsd;
  this->nne = nne;
  this->elem_type = element_type(nne, nsd);
  
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   
  if(this->curt_elem_id==0 && myrank == 0)
  {  
    FemLibBoundary fes(this, 1, i_order);
    std::cout << "Element type " << this->elem_type << " has boundary element " << fes.elem_type << endl;    
  }

  int e_order = LinearElement;

  if(this->elem_type == QTETRAHEDRON)
    e_order = QuadraticElement;  

  this->intg_order = i_order;
  if(i_order == 0){
    // these elements need at least 1st order accuracy
    if(this->elem_type == QUADRILATERAL || this->elem_type == HEXAHEDRAL    ||
       this->elem_type == QTRIANGLE     || this->elem_type == QTETRAHEDRON)
      this->intg_order = 1;
  }  
      
  this->nint = this->determine_integration_type(this->elem_type, e_order, this->intg_order);

  this->ksi.initialization(this->nint,1);
  this->eta.initialization(this->nint,1);
  this->zet.initialization(this->nint,1);
  this->weights.initialization(this->nint,1);

  this->N.initialization(nne ,1);
  this->dN.initialization(nne ,nsd);
  this->x_ip.initialization(nsd ,1);
  this->node_coord.initialization(nne ,nsd);
  this->node_id.initialization(nne,1);

  //currently supports for TETRAHEDRON
  switch(this->elem_type)
  {
    case TETRAHEDRON:  // intended to flow to next
    case QTETRAHEDRON:
    {
      bool no_Gaussian = false;
      switch(this->intg_order){
        case 0:
        {
          QuadratureRuleVolume<TetIntegrationPoints1> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd, no_Gaussian);
          break;
        }
        case 1:
        {
          QuadratureRuleVolume<TetIntegrationPoints4> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd, no_Gaussian);
          break;
        }
        case 2:
        {
          QuadratureRuleVolume<TetIntegrationPoints5> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd, no_Gaussian);
          break;
        }
        case 3:
        {
          QuadratureRuleVolume<TetIntegrationPoints5> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd, no_Gaussian);
          break;
        }
        default:
          PGFEM_printerr("ERROR: Integration order %d is not supported.\n", this->intg_order);
          print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
      }
      
      break;
    }
    case QUADRILATERAL: // intended to flow
    case HEXAHEDRAL:
    {
      const int gpno = this->intg_order + 1; // number of Gaussian points
      switch(gpno){
        case 1:
        {
          QuadratureRuleVolume<GaussIntegrationPoints1> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd);
          break;
        }
        case 2:
        {
          QuadratureRuleVolume<GaussIntegrationPoints2> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd);
          break;
        }
        case 3:
        {
          QuadratureRuleVolume<GaussIntegrationPoints3> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd);
          break;
        }
        default:
          PGFEM_printerr("ERROR: Number of Gaussian integration points %d is not implemented.\n", gpno);
          print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
      }
      break;                 
    }
  }
  
  this->ST_tensor = (double ****) aloc4(3,3,nsd,nne);
  this->ST = (double *) aloc1(3*3*nsd*nne);
}

/// set FEM libray by element
///
/// Set the FEM libray: shape function, deriviative of the shape functions, weight,
/// Gaussian points, node ids, element id, and node coordinates
/// number of node on an element and integration order determine size of member varialbes (Memory will be allocated)
/// FEMLIB_destruct should be called to cleanup this structure.
///
/// \param[in] e element id
/// \param[in] elem list of Element
/// \param[in] node list of Node
/// \param[in] i_order integration order, 0: linear, 1: quadratic, and 2: higher
/// \param[in] is_total if 1: total Lagrangian, 0: updated Lagrangian
/// \param[in] add_bubble if ture, set FEM libray by element with bubble
/// \return void
void
FEMLIB::initialization(int e,
                       const Element *elem,
                       const Node *node,
                       int i_order,
                       int is_total,
                       double *pF0I,
                       bool add_bubble)
{
  int nne = elem[e].toe;
  int nne_t = nne;
  int nsd = 3;

  if(add_bubble)
    nne_t = nne + elem[e].n_bub;

  initialization(nne, nsd, i_order);

  long *nod = (this->node_id).m_pdata;  // no memory is allocated
  elemnodes(e,nne,nod,elem);

  this->temp_v.set_variable_size(nne_t, nne);
  
  double *x = (this->temp_v).x.m_pdata; // no memory is allocated
  double *y = (this->temp_v).y.m_pdata;
  double *z = (this->temp_v).z.m_pdata;

  if(is_total)
  {
    if(pF0I != NULL)
    {
      this->u0.initialization(nne_t, 3, 0.0);
      Matrix<double> X(nne_t, 1), Y(nne_t, 1), Z(nne_t, 1);
      nodecoord_total(nne,nod,node,X.m_pdata,Y.m_pdata,Z.m_pdata);
      for(int ia=0; ia<nne; ++ia)
      {
        this->u0(ia,0) = X(ia) - (pF0I[0]*X(ia) + pF0I[1]*Y(ia) + pF0I[2]*Z(ia));
        this->u0(ia,1) = Y(ia) - (pF0I[3]*X(ia) + pF0I[4]*Y(ia) + pF0I[5]*Z(ia));
        this->u0(ia,2) = Z(ia) - (pF0I[6]*X(ia) + pF0I[7]*Y(ia) + pF0I[8]*Z(ia));
        (this->temp_v).x(ia) = X(ia) - this->u0(ia, 0);
        (this->temp_v).y(ia) = Y(ia) - this->u0(ia, 1);
        (this->temp_v).z(ia) = Z(ia) - this->u0(ia, 2);
      }
    }
    else
      nodecoord_total(nne,nod,node,x,y,z);
  }
  else
    nodecoord_updated(nne,nod,node,x,y,z);

  if(add_bubble)
    element_center(nne,x,y,z);

  for(int a=0; a<nne; a++)
  {
    this->node_coord(a, 0) = x[a];
    this->node_coord(a, 1) = y[a];
    this->node_coord(a, 2) = z[a];
  }

  this->curt_elem_id = e;
}

void
FEMLIB::elem_shape_function(long ip, int nne, double *N)
{
  if(nne == 1){ // constant
    N[0] = 1.0;
    return;
  }
  
  double ksi_ = {};
  double eta_ = {};
  double zet_ = {};
  switch(this->nsd){
    case 3: // intended to flow
      zet_ = this->zet(this->itg_ids(ip, 2));
    case 2: // intended to flow
      eta_ = this->eta(this->itg_ids(ip, 1));
    case 1: 
      ksi_ = this->ksi(this->itg_ids(ip, 0));
  }
  FEM_shape_functions(this->elem_type, ksi_, eta_, zet_, N);
}

double
FEMLIB::compute_integration_weight(const int ip){
  double wt = 0.0;
  switch(this->elem_type){
    case TRIANGLE:    // intended to flow for triangle type of element
    case QTRIANGLE:
    case TETRAHEDRON:
    case QTETRAHEDRON:
      wt = this->weights(this->itg_ids(ip, 0));
      break;
    case LINE:
    case QLINE:
      wt = this->weights(this->itg_ids(ip, 0));
      break;      
    case QUADRILATERAL:
      wt = this->weights(this->itg_ids(ip, 0))
          *this->weights(this->itg_ids(ip, 1));
      break;
    case HEXAHEDRAL:
      wt = this->weights(this->itg_ids(ip, 0))
          *this->weights(this->itg_ids(ip, 1))
          *this->weights(this->itg_ids(ip, 2));
      break;
    default:
      PGFEM_printerr("ERROR: FEM integration for element %d is not implemented.\n", this->elem_type);
      print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
  }
  return wt;
}

void
FEMLIB::elem_basis_V(long ip)
{
  this->curt_itg_id = ip;

  double ksi_ = {}, eta_ = {}, zet_ = {};
  switch(this->nsd){
    case 3: // intended to flow
      zet_ = this->zet(this->itg_ids(ip, 2));
    case 2: // intended to flow
      eta_ = this->eta(this->itg_ids(ip, 1));
    case 1: 
      ksi_ = this->ksi(this->itg_ids(ip, 0));
  }
  FEM_shape_functions(this->elem_type, ksi_, eta_, zet_, this->N.m_pdata);
  double wt = this->compute_integration_weight(ip);

  this->temp_v.ksi_ip = ksi_;
  this->temp_v.eta_ip = eta_;
  this->temp_v.zet_ip = zet_;
  this->temp_v.w_ip   = wt;

  this->detJ = deriv(ksi_, eta_, zet_, this->nne,
                     this->temp_v.x.m_pdata,   this->temp_v.y.m_pdata,   this->temp_v.z.m_pdata,
                     this->temp_v.N_x.m_pdata, this->temp_v.N_y.m_pdata, this->temp_v.N_z.m_pdata);

  this->x_ip.set_values(0.0);
  for(int a = 0; a<this->nne; a++)
  {
    this->dN(a, 0) = this->temp_v.N_x(a);
    this->dN(a, 1) = this->temp_v.N_y(a);
    this->dN(a, 2) = this->temp_v.N_z(a);

    this->x_ip(0) += this->N(a)*this->node_coord(a,0);
    this->x_ip(1) += this->N(a)*this->node_coord(a,1);
    this->x_ip(2) += this->N(a)*this->node_coord(a,2);
  }
  this->detJxW = this->detJ*wt;
}

void
FEMLIB::update_shape_tensor(void)
{
  shape_tensor(this->nne,this->nsd,this->temp_v.N_x.m_pdata,
               this->temp_v.N_y.m_pdata,
               this->temp_v.N_z.m_pdata,this->ST_tensor);
  shapeTensor2array(this->ST,CONST_4(double) this->ST_tensor,this->nne);
}

void
FEMLIB::update_deformation_gradient(const int ndofn, double *u, double *F)
{
  double **F_mat;
  F_mat = aloc2(3,3);
  def_grad_get(this->nne,ndofn,CONST_4(double) this->ST_tensor,u,F_mat);
  mat2array(F,CONST_2(double) F_mat,3,3);
  dealoc2(F_mat,3);
}

void
FEMLIB::update_deformation_gradient(const int ndofn, double *u, double *F, double *pF0I)
{
  if(pF0I==NULL)
    this->update_deformation_gradient(ndofn, u, F);
  else
  {
    Matrix<double> r_e(ndofn*this->nne, 1, 0.0);          
    for(int ia=0; ia<this->nne; ++ia)
    {
      r_e.m_pdata[ia*ndofn+0] = u[ia*ndofn+0] + this->u0(ia,0);
      r_e.m_pdata[ia*ndofn+1] = u[ia*ndofn+1] + this->u0(ia,1);
      r_e.m_pdata[ia*ndofn+2] = u[ia*ndofn+2] + this->u0(ia,2);
    }
    this->update_deformation_gradient(ndofn, r_e.m_pdata, F);
  }
}

double
FEMLIB::elem_volume(void)
{
  double Volume;
  if(this->nne == 4){ // linear tet
    Volume = Tetra_V(this->temp_v.x.m_pdata,this->temp_v.y.m_pdata,this->temp_v.z.m_pdata);
  } else if(this->nne == 10){ // quadradic tet
    Volume = Tetra_qv_V(10,3,this->temp_v.x.m_pdata,this->temp_v.y.m_pdata,this->temp_v.z.m_pdata);
  } else if(this->nne == 8){ // trilinear hex
    Volume = Hexa_V(this->temp_v.x.m_pdata,this->temp_v.y.m_pdata,this->temp_v.z.m_pdata);
  } else Volume = 0.0;
  return Volume;
}

FEMLIB::~FEMLIB()
{
  if(ST_tensor != NULL)
    dealoc4(ST_tensor,3,3,nsd);
  if(ST != NULL)
    PGFEM_free(ST);
}

int nne_boundary(const int elem_type){
  switch(elem_type){
  case TRIANGLE:      // tri  -> line
  case QUADRILATERAL: // quad -> line  
    return 2;
  case QTRIANGLE:     // qtri -> qline
    return 3;
  case TETRAHEDRON:   // Tet --> Tri
    return 3;
  case QTETRAHEDRON:  // qTet ->> qTri
    return 6;
  case HEXAHEDRAL:    // Hex --> quad
    return 4;
  default:
    PGFEM_printerr("WARNING: unrecognized element type for 2D femlib: %s:%s:%d\n",
                     __func__,__FILE__,__LINE__);
      break;
  }
  return 0;
}

void
FemLibBoundary::set_volume_to_boundary_map(void){
  switch(this->feVol->elem_type){
    case TETRAHEDRON:
      this->Volume2Boundary = Tet2Tri[this->face_id];
      this->kez_map         = kez_map_Tet[this->face_id];
      break;
    case QTETRAHEDRON:
      this->Volume2Boundary = QTet2QTri[this->face_id];
      this->kez_map         = kez_map_Tet[this->face_id];
      break;
    case HEXAHEDRAL:
      this->Volume2Boundary = Hex2Quad[this->face_id];
      this->kez_map         = kez_map_Hex[this->face_id];
      break;
    default:
      PGFEM_printerr("ERROR: Boundary integration for "
                     "%d element is not implemented.\n", this->feVol->elem_type);
      print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
  }
}

void
FemLibBoundary::initialization(const FEMLIB *fe,
                               const int face_id,
                               const int i_order)
{
  this->feVol = fe;
  this->face_id = face_id;
  this->nsd = fe->nsd;
  this->nne = nne_boundary(fe->elem_type);
  this->curt_elem_id = fe->curt_elem_id;
  this->elem_type = element_type(this->nne, 2);
  
  int e_order = LinearElement;

  if(this->elem_type == QTRIANGLE)
    e_order = QuadraticElement; 
   
  this->intg_order = i_order;
  if(i_order == 0){
    // these elements need at least 1st order accuracy
    if(this->elem_type == QUADRILATERAL || this->elem_type == QTRIANGLE)
      this->intg_order = 1;
  }
  
  // set maps for node id and integration points from volum to boundary element
  // volume femlib(fe) needs to be set prior to this
  this->set_volume_to_boundary_map();
  
  this->nint = this->determine_integration_type(this->elem_type, e_order, this->intg_order);
    
  this->ksi.initialization(this->nint,1);
  this->eta.initialization(this->nint,1);
  this->zet.initialization(this->nint,1);
  this->weights.initialization(this->nint,1);
  
  this->N.initialization(this->feVol->nne, 1);
  this->dN.initialization(this->feVol->nne, this->feVol->nsd);
  this->x_ip.initialization(this->feVol->nsd ,1);
  this->node_coord.initialization(this->nne ,this->nsd);
  this->node_id.initialization(this->nne,1);
  
  switch(fe->elem_type){
    case TRIANGLE:
    case QUADRILATERAL:
      PGFEM_printerr("ERROR: Boundary integration for elelement %d is not implemented.\n", fe->elem_type);
      print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
    case TETRAHEDRON:  // intended to flow
    case QTETRAHEDRON:
    {
      bool no_Gaussian = false;
      if(this->intg_order > 3){
        PGFEM_printerr("ERROR: Integration order %d is not supported for element %d\n", this->intg_order, this->elem_type);
        print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
      }    
      switch(this->intg_order){
        case 0:
        {
          QuadratureRuleBoundary<TriIntegrationPoints1> qr;
          qr.get_integration_points(ksi, eta, zet, weights, face_id, fe->elem_type);
          qr.get_integration_point_ids(this->itg_ids, this->nsd, face_id, no_Gaussian);
          break;
        }
        case 1:
        {
          QuadratureRuleBoundary<TriIntegrationPoints3> qr;
          qr.get_integration_points(ksi, eta, zet, weights, face_id, fe->elem_type);
          qr.get_integration_point_ids(this->itg_ids, this->nsd, face_id, no_Gaussian);
          break;
        }
        case 2:
        {
          QuadratureRuleBoundary<TriIntegrationPoints4> qr;
          qr.get_integration_points(ksi, eta, zet, weights, face_id, fe->elem_type);
          qr.get_integration_point_ids(this->itg_ids, this->nsd, face_id, no_Gaussian);
          break;
        }                
      }
      break;
    }
    case HEXAHEDRAL:
    {
      const int gpno = this->intg_order + 1;
      if(gpno > 3){
        PGFEM_printerr("ERROR: Number of Gaussian integration points %d is not supported.\n", gpno);
        print_func_file_line_and_abort(__func__,__FILE__,__LINE__);        
      }
      
      switch(gpno){
        case 1:
        {
          QuadratureRuleBoundary<GaussIntegrationPoints1> qr;
          qr.get_integration_points(ksi, eta, zet, weights, face_id, fe->elem_type);
          qr.get_integration_point_ids(this->itg_ids, this->nsd, face_id);
          break;
        }
        case 2:
        {
          QuadratureRuleBoundary<GaussIntegrationPoints2> qr;
          qr.get_integration_points(ksi, eta, zet, weights, face_id, fe->elem_type);
          qr.get_integration_point_ids(this->itg_ids, this->nsd, face_id);
          break;
        }
        case 3:
        {
          QuadratureRuleBoundary<GaussIntegrationPoints3> qr;
          qr.get_integration_points(ksi, eta, zet, weights, face_id, fe->elem_type);
          qr.get_integration_point_ids(this->itg_ids, this->nsd, face_id);
          break;
        }
      }
      break;
    }
    default:
      PGFEM_printerr("ERROR: Boundary integration for "
                     "%d element is supported.\n", fe->elem_type);
      print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
      break;
  }
  
  if(this->nsd == 3){
    this->ST_tensor = (double ****) aloc4(3,3,nsd,nne);
    this->ST = (double *) aloc1(3*3*nsd*nne);
  }

  this->temp_v.set_variable_size(this->nne, this->nne);
  
  double *x = (this->temp_v).x.m_pdata; // no memory is allocated
  double *y = (this->temp_v).y.m_pdata;
  double *z = (this->temp_v).z.m_pdata;  

  for(int ia=0; ia<this->nne; ++ia){
    int id_V2B = this->Volume2Boundary[ia];
    this->node_id(ia) = this->feVol->node_id.m_pdata[id_V2B];
    x[ia] = fe->temp_v.x.m_pdata[id_V2B];
    y[ia] = fe->temp_v.y.m_pdata[id_V2B];
    z[ia] = fe->temp_v.z.m_pdata[id_V2B];        

    this->node_coord(ia, 0) = x[ia];
    this->node_coord(ia, 1) = y[ia];
    this->node_coord(ia, 2) = z[ia];
  }
}

double FemLibBoundary::compute_integration_weight(const int ip){
  double wt = 0.0;
  const int itg_ids_0 = this->itg_ids(ip, kez_map[0]);  

  switch(this->feVol->elem_type){
    case TRIANGLE:    // intended to flow for triangle type of element
    case QTRIANGLE:
    case TETRAHEDRON:
    case QTETRAHEDRON:
    case QUADRILATERAL:
    case LINE:
    case QLINE:      
      wt = this->weights(itg_ids_0);
      break;      
    case HEXAHEDRAL:
    {
      const int itg_ids_1 = this->itg_ids(ip, kez_map[1]);
      wt = this->weights(itg_ids_0)
          *this->weights(itg_ids_1);
      break;
    }
    default:
      PGFEM_printerr("ERROR: FEM integration for element %d is not implemented.\n", this->elem_type);
      print_func_file_line_and_abort(__func__,__FILE__,__LINE__);
  }
  return wt;  
}

void
FemLibBoundary::elem_basis_S(const int ip)
{
  this->curt_itg_id = ip;

  double ksi_ = {}, eta_ = {}, zet_ = {};
  switch(this->nsd){
    case 3: // intended to flow
      zet_ = this->zet(this->itg_ids(ip, 2));
    case 2: // intended to flow
      eta_ = this->eta(this->itg_ids(ip, 1));
    case 1: 
      ksi_ = this->ksi(this->itg_ids(ip, 0));
  }
  FEM_shape_functions(this->feVol->elem_type, ksi_, eta_, zet_, this->N.m_pdata);
  double wt = this->compute_integration_weight(ip);

  this->temp_v.ksi_ip = ksi_;
  this->temp_v.eta_ip = eta_;
  this->temp_v.zet_ip = zet_;
  this->temp_v.w_ip   = wt;
}