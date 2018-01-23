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

static const constexpr int          LINE = 1;
static const constexpr int      TRIANGLE = 3;
static const constexpr int QUADRILATERAL = 4;
static const constexpr int   TETRAHEDRON = 4;
static const constexpr int    HEXAHEDRAL = 8;
static const constexpr int  QTETRAHEDRON = 10;

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

long
FEMLIB::determine_integration_type(int e_type, int i_order)
{
  switch(i_order)
  {
   case 0:
    switch(e_type)
    {
     case LINE:
      return 1;
     case TRIANGLE:
      return 1;
     case TETRAHEDRON: // QUADRILATERAL
      return 1;
     default:
      return 1;
    }
    break;
   case 1:
    switch(e_type)
    {
     case LINE:
      return 2;
     case TRIANGLE:
      return 3;
     case TETRAHEDRON: //QUADRILATERAL
      return 4;
     case HEXAHEDRAL:
      return 8;
     case QTETRAHEDRON:
      return 4;
     default:
      return 4;
    }
    break;
   case 2:
    switch(e_type)
    {
     case TETRAHEDRON:
      return 5;
     case QTETRAHEDRON:
      return 5;
     default:
      return 5;
    }
    break;
   case 3:
    switch(e_type)
    {
     case TETRAHEDRON:
      return 11;
     case QTETRAHEDRON:
      return 11;
     default:
      return 11;
    }
    break;
  }
  return 4;
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

  if(i_order<1)
  {
    if(nne==QTETRAHEDRON)
      i_order = 1;
    if(nne==HEXAHEDRAL)
      i_order = 1;
  }

  nint = this->determine_integration_type(e_type, i_order);
  this->nint = nint;

  this->ksi.initialization(nint,1);
  this->eta.initialization(nint,1);
  this->zet.initialization(nint,1);
  this->weights.initialization(nint,1);
  this->N.initialization(nne ,1);
  this->dN.initialization(nne ,nsd);
  this->x_ip.initialization(nsd ,1);
  this->node_coord.initialization(nne ,nsd);
  this->node_id.initialization(nne,1);

  long npt_x, npt_y, npt_z;

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
  for(int a=1; a<=npt_x; a++)
  {
    for(int b=1; b<=npt_y; b++)
    {
      for(int c=1; c<=npt_z; c++)
      {
        cnt++;
        this->itg_ids(cnt, 1) = a;
        this->itg_ids(cnt, 2) = b;
        this->itg_ids(cnt, 3) = c;
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
    nodecoord_total(nne,nod,node,x,y,z);
  else
    nodecoord_updated(nne,nod,node,x,y,z);

  if(add_bubble)
    element_center(nne,x,y,z);

  for(int a=0; a<nne; a++)
  {
    this->node_coord(a+1, 1) = x[a];
    this->node_coord(a+1, 2) = y[a];
    this->node_coord(a+1, 3) = z[a];
  }

  this->curt_elem_id = e;
}

void
FEMLIB::elem_shape_function(long ip, int nne, double *N)
{
  double ksi_, eta_, zet_;

  int itg_id_ip_1 = this->itg_ids(ip, 1);
  int itg_id_ip_2 = this->itg_ids(ip, 2);
  int itg_id_ip_3 = this->itg_ids(ip, 3);

  if(this->elem_type == HEXAHEDRAL)
  {// hexahedron
    ksi_ = this->ksi(itg_id_ip_1, 1);
    eta_ = this->eta(itg_id_ip_2, 1);
    zet_ = this->zet(itg_id_ip_3, 1);
  }
  else
  { // tetrahedron type
    ksi_ = this->ksi(itg_id_ip_3, 1);
    eta_ = this->eta(itg_id_ip_3, 1);
    zet_ = this->zet(itg_id_ip_3, 1);
  }

  shape_func(ksi_, eta_, zet_, nne, N);
}

void
FEMLIB::elem_basis_V(long ip)
{
  this->curt_itg_id = ip;
  double ksi_, eta_, zet_, wt;

  int itg_id_ip_1 = this->itg_ids(ip, 1);
  int itg_id_ip_2 = this->itg_ids(ip, 2);
  int itg_id_ip_3 = this->itg_ids(ip, 3);

  if(this->elem_type == HEXAHEDRAL)
  {// hexahedron
    ksi_ = this->ksi(itg_id_ip_1, 1);
    eta_ = this->eta(itg_id_ip_2, 1);
    zet_ = this->zet(itg_id_ip_3, 1);

    wt = this->weights(itg_id_ip_1, 1)*this->weights(itg_id_ip_2, 1)*this->weights(itg_id_ip_3, 1);
  }
  else
  { // tetrahedron type
    ksi_ = this->ksi(itg_id_ip_3, 1);
    eta_ = this->eta(itg_id_ip_3, 1);
    zet_ = this->zet(itg_id_ip_3, 1);
    wt = this->weights(itg_id_ip_3, 1);
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
  for(int a = 1; a<=this->nne; a++)
  {
    this->dN(a, 1) = this->temp_v.N_x(a, 1);
    this->dN(a, 2) = this->temp_v.N_y(a, 1);
    this->dN(a, 3) = this->temp_v.N_z(a, 1);

    this->x_ip(1,1) += this->N(a,1)*this->node_coord(a,1);
    this->x_ip(2,1) += this->N(a,1)*this->node_coord(a,2);
    this->x_ip(3,1) += this->N(a,1)*this->node_coord(a,3);
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
    free(ST);
}
