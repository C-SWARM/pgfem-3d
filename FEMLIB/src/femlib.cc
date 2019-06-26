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
#include "integrate_surface.h"
#include "quadrature_rules.h"

static const constexpr int          LINE    = element_type( 2, FemDim1D);
static const constexpr int      TRIANGLE    = element_type( 3, FemDim2D);
static const constexpr int     QTRIANGLE    = element_type( 6, FemDim2D);
static const constexpr int QUADRILATERAL    = element_type( 4, FemDim2D);
static const constexpr int   TETRAHEDRON    = element_type( 4, FemDim3D);
static const constexpr int  QTETRAHEDRON    = element_type(10, FemDim3D);
static const constexpr int    HEXAHEDRAL    = element_type( 8, FemDim3D);

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



template <class GP> class QuadratureRule{
  public:
    GP gp;
    QuadratureRule(){};
    void get_integration_points(Matrix<double> &ksi, Matrix<double> &eta, Matrix<double> &zet, Matrix<double> &weights){
      for(int ia=0; ia<gp.gpno; ++ia){
        ksi(ia) = gp.gk[ia];
        eta(ia) = gp.ge[ia];
        zet(ia) = gp.gz[ia];
        weights(ia) = gp.weights[ia];
      }
    }

    void get_integration_point_ids_3D(Matrix<int> &itg_ids, 
                                      const int npt_x,
                                      const int npt_y,
                                      const int npt_z){
      const int nsd = 3;
      itg_ids.initialization(npt_x*npt_y*npt_z, nsd);
          
      long cnt = 0;
      for(int ia=0; ia<npt_x; ia++){
        for(int ib=0; ib<npt_y; ib++){
          for(int ic=0; ic<npt_z; ic++){
            itg_ids(cnt, 0) = ia;
            itg_ids(cnt, 1) = ib;
            itg_ids(cnt, 2) = ic;
            ++cnt;
          }
        }
      }
    }
    
    void get_integration_point_ids_2D(Matrix<int> &itg_ids, 
                                      const int npt_x,
                                      const int npt_y){
      const int nsd = 2;
      itg_ids.initialization(npt_x*npt_y, nsd);
      
      long cnt = 0;
      for(int ia=0; ia<npt_x; ia++){
        for(int ib=0; ib<npt_y; ib++){
            itg_ids(cnt, 0) = ia;
            itg_ids(cnt, 1) = ib;
            ++cnt;
        }
      }
    }
    void get_integration_point_ids_1D(Matrix<int> &itg_ids, 
                                      const int npt_x){

      const int nsd = 1;
      itg_ids.initialization(npt_x, nsd);
      
      for(int ia=0; ia<npt_x; ia++){
        itg_ids(ia, 0) = ia;
      }
    }    

    
    void get_integration_point_ids(Matrix<int> &itg_ids, const int nsd, bool is_Gaussian = true){
      int npt[3] = {};
      if(is_Gaussian){
        for(int ia=0; ia<nsd; ++ia)
          npt[ia] = this->gp.gpno;          
      } else {
        npt[0] = this->gp.gpno;
        for(int ia=1; ia<nsd; ++ia)
          npt[ia] = 1;
      }
      switch(nsd){
        case 1:
          get_integration_point_ids_1D(itg_ids, npt[0]);
          return;
        case 2:
          get_integration_point_ids_2D(itg_ids, npt[0], npt[1]);
          return;
        case 3:
          get_integration_point_ids_3D(itg_ids, npt[0], npt[1], npt[2]);
          return;
        default:
          PGFEM_printerr("ERROR: number of spatial dimension %d is not supported.\n", nsd);
          pgfem3d::PGFEM_Abort();
      }      
    }
};

template <class GP> class QuadratureRule3Dto2D{
  public:
    GP gp;
    QuadratureRule3Dto2D(){};
    void get_integration_points(Matrix<double> &ksi, Matrix<double> &eta, Matrix<double> &zet, Matrix<double> &weights){
      for(int ia=0; ia<gp.gpno; ++ia){
        ksi(ia) = gp.gk[ia];
        eta(ia) = gp.ge[ia];
        zet(ia) = gp.gz[ia];
        weights(ia) = gp.weights[ia];
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
      PGFEM_printerr("ERROR: number of integration points: "
                     " order %d is not yet implemented for TRIANGLE element.\n", order);
      pgfem3d::PGFEM_Abort();
      break;
  }
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
      PGFEM_printerr("ERROR: number of integration points: "
                     " order %d is not yet implemented for QUADRILATERAL element.\n", order);
      pgfem3d::PGFEM_Abort();
      break;
  }
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
                     " order %d is not yet implemented for TETRAHEDRON element.\n", order);
      pgfem3d::PGFEM_Abort();
      break;
  }
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
                     " order %d is not yet implemented for HEXAHEDRAL element.\n", order);
      pgfem3d::PGFEM_Abort();
      break;
  }
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
      PGFEM_printerr("ERROR: element type %d is not yet implemented for finite element integration.\n", e_type);
      pgfem3d::PGFEM_Abort();
  }
}

void
FEMLIB::initialization(const int nne,
                       const int nsd,
                       const int i_order)
{
  this->nsd = nsd;
  this->nne = nne;
  this->elem_type = element_type(nne, nsd);

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
      if(this->intg_order > 3){
        PGFEM_printerr("ERROR: Integration order %d is not supported.\n", this->intg_order);
        pgfem3d::PGFEM_Abort();        
      }

      switch(this->intg_order){
        case 0:
        {
          QuadratureRule<TetIntegrationPoints1> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd);
          break;
        }
        case 1:
        {
          QuadratureRule<TetIntegrationPoints4> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd);
          break;
        }
        case 2:
        {
          QuadratureRule<TetIntegrationPoints5> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd);
          break;
        }
        case 3:
        {
          QuadratureRule<TetIntegrationPoints5> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd);
          break;
        }
      }
      
      break;
    }
    case QUADRILATERAL: // intended to flow
    case HEXAHEDRAL:
    {
      const int gpno = this->intg_order + 1;
      if(gpno > 3){
        PGFEM_printerr("ERROR: Number of Gaussian integration points %d is not supported.\n", gpno);
        pgfem3d::PGFEM_Abort();        
      }
      
      switch(gpno){
        case 1:
        {
          QuadratureRule<GaussIntegrationPoints1> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd);
          break;
        }
        case 2:
        {
          QuadratureRule<GaussIntegrationPoints2> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd);
          break;
        }
        case 3:
        {
          QuadratureRule<GaussIntegrationPoints3> qr;
          qr.get_integration_points(ksi, eta, zet, weights);
          qr.get_integration_point_ids(this->itg_ids, this->nsd);
          break;
        }
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
      for(int ia=0; ia<nne; ia++)
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
  double ksi_, eta_, zet_;

  int itg_id_ip_1 = this->itg_ids(ip, 0);
  int itg_id_ip_2 = this->itg_ids(ip, 1);
  int itg_id_ip_3 = this->itg_ids(ip, 2);

  if(this->elem_type == HEXAHEDRAL)
  {// hexahedron
    ksi_ = this->ksi(itg_id_ip_1);
    eta_ = this->eta(itg_id_ip_2);
    zet_ = this->zet(itg_id_ip_3);
  }
  else
  { // tetrahedron type
    ksi_ = this->ksi(itg_id_ip_3);
    eta_ = this->eta(itg_id_ip_3);
    zet_ = this->zet(itg_id_ip_3);
  }

  shape_func(ksi_, eta_, zet_, nne, N);
}

void
FEMLIB::elem_basis_V(long ip)
{
  this->curt_itg_id = ip;
  double ksi_, eta_, zet_, wt;

  int itg_id_ip_1 = this->itg_ids(ip, 0);
  int itg_id_ip_2 = this->itg_ids(ip, 1);
  int itg_id_ip_3 = this->itg_ids(ip, 2);

  if(this->elem_type == HEXAHEDRAL)
  {// hexahedron
    ksi_ = this->ksi(itg_id_ip_1);
    eta_ = this->eta(itg_id_ip_2);
    zet_ = this->zet(itg_id_ip_3);

    wt = this->weights(itg_id_ip_1)*this->weights(itg_id_ip_2)*this->weights(itg_id_ip_3);
  }
  else
  { // tetrahedron type
    ksi_ = this->ksi(itg_id_ip_3);
    eta_ = this->eta(itg_id_ip_3);
    zet_ = this->zet(itg_id_ip_3);
    wt = this->weights(itg_id_ip_3);
  }

  this->temp_v.ksi_ip = ksi_;
  this->temp_v.eta_ip = eta_;
  this->temp_v.zet_ip = zet_;
  this->temp_v.w_ip   = wt;

  shape_func(ksi_, eta_, zet_, this->nne, this->N.m_pdata);

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
    for(int ia=0; ia<this->nne; ia++)
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

int nne_2D(const int nne_3D){
  switch(nne_3D){
  case TETRAHEDRON:  // Tet --> Tri
    return 3;
  case QTETRAHEDRON: // qTet ->> qTri
    return 6;
  case HEXAHEDRAL: // Hex --> quad
    return 4;
  default:
    PGFEM_printerr("WARNING: unrecognized element type for 2D femlib: %s:%s:%d\n",
                     __func__,__FILE__,__LINE__);
      break;
  }
  return 0;
}


void
FemLib2D::initialization(const FEMLIB *fe,
                         const int face_id,
                         const int i_order)
{
  this->nsd = 3;
  this->nne = nne_2D(fe->nne);
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

  this->nint = this->determine_integration_type(this->elem_type, e_order, this->intg_order);
    
  this->ksi.initialization(this->nint,1);
  this->eta.initialization(this->nint,1);
  this->zet.initialization(this->nint,1);
  this->weights.initialization(this->nint,1);
  
  this->N.initialization(this->nne, 1);
  this->dN.initialization(this->nne, this->nsd);
  this->x_ip.initialization(this->nsd ,1);
  this->node_coord.initialization(this->nne ,this->nsd);
  this->node_id.initialization(this->nne,1);
  
//  int tmp;
  switch(this->elem_type){
    case TRIANGLE:
//      get_tria_quadrature_rule(this->intg_order, &tmp,
//                               this->ksi.m_pdata, 
//                               this->eta.m_pdata,
//                               this->weights.m_pdata);
      break;
    case QTETRAHEDRON:
//      get_tria_quadrature_rule(this->intg_order, &tmp,
//                               this->ksi.m_pdata, 
//                               this->eta.m_pdata,
//                               this->weights.m_pdata);
    case HEXAHEDRAL:
    default:
      PGFEM_printerr("ERROR: 3D surface integration: "
                     "%d element is supported.\n", fe->elem_type);
      pgfem3d::PGFEM_Abort();
      break;
  }
  /*
  this->fe3D = fe;

  int int_order = 1;

  if(this->fe3D->nne == 10)
    int_order = 2;

  this->nne = nne_2D(this->fe3D->nne);

  Matrix<double> ksi_3D;
  integrate_surface(this->fe3D->nne, face_id, int_order,
                    &(this->nint),&ksi_3D,&eta_3D,&zet_3D,
                               &ksi_2D,&eta_2D,&wt_2D,
                               &nne_2D,&nod_2D);
*/
}
