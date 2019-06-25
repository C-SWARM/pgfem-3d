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

static const constexpr int          LINE    = 2;
static const constexpr int      TRIANGLE    = 3;
static const constexpr int QUADRILATERAL    = 4;
static const constexpr int   TETRAHEDRON    = 4;
static const constexpr int    HEXAHEDRAL    = 8;
static const constexpr int LinearElement    = 0;
static const constexpr int QuadraticElement = 1;

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
      return 1;
    case 1:
      return 3;
    case 2:
      return 4;
    default:
      PGFEM_printerr("ERROR: number of integration points: "
                     " order %d is not yet implemented for TRIANGLE element.\n", order);
      PGFEM_Abort();
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
      PGFEM_Abort();
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
      PGFEM_Abort();
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
      PGFEM_Abort();
      break;
  }
}

long
FEMLIB::determine_integration_type(const int e_type,
                                   const int e_order,
                                   const int order)
{
  int i_order = order;

  if(order<1 && e_order>0)
    i_order = 1;

  if(e_type == HEXAHEDRAL || e_type == QUADRILATERAL)
    i_order = 1;

  switch(e_type){
    case LINE:
      return number_of_integration_points_line(i_order);
    case TRIANGLE:
      return number_of_integration_points_tri(i_order);
    case QUADRILATERAL:
      return number_of_integration_points_quad(i_order);
    case TETRAHEDRON:
      return number_of_integration_points_tet(i_order);
    case HEXAHEDRAL:
      return number_of_integration_points_hex(i_order);
    default:
      PGFEM_printerr("ERROR: element type %d is not yet implemented for finite element integration.\n", e_type);
      PGFEM_Abort();
  }
}

void
FEMLIB::initialization(int e_type, int i_order, int nne)
{
  int nsd = 3;
  long nint;

  this->nsd = nsd;
  this->nne = nne;
  this->elem_type = e_type;
  this->intg_order = i_order;

  int e_order = 0;

  if(nne==QTETRAHEDRON)
    e_order = 1;

  this->nint = this->determine_integration_type(e_type, e_order, i_order);

  this->ksi.initialization(this->nint,1);
  this->eta.initialization(this->nint,1);
  this->zet.initialization(this->nint,1);
  this->weights.initialization(this->nint,1);

  this->N.initialization(nne ,1);
  this->dN.initialization(nne ,nsd);
  this->x_ip.initialization(nsd ,1);
  this->node_coord.initialization(nne ,nsd);
  this->node_id.initialization(nne,1);

  long npt_x = 0, npt_y = 0, npt_z = 0;

  //currently supports for TETRAHEDRON
  switch(e_type)
  {
   case TETRAHEDRON:
   case QTETRAHEDRON:
     {
       switch(i_order)
       {

        case 0:
         int_tetra_1(this->ksi.m_pdata, this->eta.m_pdata, this->zet.m_pdata,
                     this->weights.m_pdata);
         npt_x = 1;
         npt_y = 1;
         npt_z = 1;
         break;
        case 1:
         int_tetra_4(this->ksi.m_pdata, this->eta.m_pdata, this->zet.m_pdata,
                     this->weights.m_pdata);
         npt_x = 1;
         npt_y = 1;
         npt_z = 4;
         break;
        case 2:
         int_tetra_5(this->ksi.m_pdata, this->eta.m_pdata, this->zet.m_pdata,
                     this->weights.m_pdata);
         npt_x = 1;
         npt_y = 1;
         npt_z = 5;
         break;
        case 3:
         int_tetra_11(this->ksi.m_pdata, this->eta.m_pdata, this->zet.m_pdata,
                      this->weights.m_pdata);
         npt_x = 1;
         npt_y = 1;
         npt_z = 11;

         break;
        default:
         integrate(e_type, &npt_x, &npt_y, &npt_z,
                   this->ksi.m_pdata, this->eta.m_pdata, this->zet.m_pdata,
                   this->weights.m_pdata);
       }

       this->itg_ids.initialization(npt_x*npt_y*npt_z, nsd);
       break;
     }
   case HEXAHEDRAL:
     {
       double gk[5]; // maximum integration order = 5
       double w[5];  // maximum integration order = 5
       int n = 0;
       switch(i_order)
       {
        case 1:
         intpoints_2(gk,w);
         n = 2;
         break;
        case 2:
         intpoints_3(gk,w);
         n=3;
         break;
        case 3:
        default:
         integrate(e_type, &npt_x, &npt_y, &npt_z, gk, NULL, NULL, w);
       }
       npt_x = npt_y = npt_z = n;
       for(int ia=0; ia<n; ia++)
       {
         this->ksi.m_pdata[ia] = gk[ia];
         this->eta.m_pdata[ia] = gk[ia];
         this->zet.m_pdata[ia] = gk[ia];
         this->weights.m_pdata[ia] = w[ia];
       }

       this->itg_ids.initialization(npt_x*npt_y*npt_z, nsd);
       break;
     }

  }


  long cnt = 0;
  for(int a=0; a<npt_x; a++)
  {
    for(int b=0; b<npt_y; b++)
    {
      for(int c=0; c<npt_z; c++)
      {
        this->itg_ids(cnt, 0) = a;
        this->itg_ids(cnt, 1) = b;
        this->itg_ids(cnt, 2) = c;
        ++cnt;
      }
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

  if(add_bubble)
    nne_t = nne + elem[e].n_bub;

  initialization(nne, i_order, nne);

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
  case 4:  // Tet --> Tri
    return 3;
  case 10: // qTet ->> qTri
    return 6;
  case 8: // Hex --> quad
    return 4;
  default:
    PGFEM_printerr("WARNING: unrecognized element type for 2D femlib: %s:%s:%d\n",
                     __func__,__FILE__,__LINE__);
      break;
  }
}

void
FemLib2D::initialization(const FEMLIB *fe,
                         const int face_id)
{
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

}
