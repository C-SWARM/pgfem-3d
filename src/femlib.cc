#include "femlib.h"
#include "elem3d.h"

#include "data_structure_c.h"
#include "utils.h"
#include "cast_macros.h"
#include "allocation.h"
#include "tensors.h"
#include "def_grad.h"

#define LINE           1
#define TRIANGLE       3
#define QUADRILATERAL  4
#define TETRAHEDRON    4
#define HEXAHEDRAL     8
#define QTETRAHEDRON   10

long FEMLIB_determine_integration_type(int e_type, int i_order)
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

void FEMLIB_set_variable_size(TEMP_VARIABLES *v, int nne_t, int nne)
{
  Matrix_construct_redim(double,v->x  ,nne_t,1);
  Matrix_construct_redim(double,v->y  ,nne_t,1);
  Matrix_construct_redim(double,v->z  ,nne_t,1);
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

  if(i_order<1)
  {
    if(nne==QTETRAHEDRON)
      i_order = 1;
    if(nne==HEXAHEDRAL)
      i_order = 1;
  }

  nint = FEMLIB_determine_integration_type(e_type, i_order);
  fe->nint = nint;

  Matrix_construct_redim(double,fe->ksi       ,nint,1); memset((fe->ksi).m_pdata,0,sizeof(double)*nint);
  Matrix_construct_redim(double,fe->eta       ,nint,1); memset((fe->eta).m_pdata,0,sizeof(double)*nint);
  Matrix_construct_redim(double,fe->zet       ,nint,1); memset((fe->zet).m_pdata,0,sizeof(double)*nint);
  Matrix_construct_redim(double,fe->weights   ,nint,1);
  Matrix_construct_redim(double,fe->N         ,nne ,1);
  Matrix_construct_redim(double,fe->dN        ,nne ,nsd);
  Matrix_construct_redim(double,fe->x_ip      ,nsd ,1);
  Matrix_construct_redim(double,fe->node_coord,nne ,nsd);
  Matrix_construct_redim(long,fe->node_id,nne,1);

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
          int_tetra_1(fe->ksi.m_pdata, fe->eta.m_pdata, fe->zet.m_pdata,
                      fe->weights.m_pdata);
          npt_x = 1;
          npt_y = 1;
          npt_z = 1;
          break;
        case 1:
          int_tetra_4(fe->ksi.m_pdata, fe->eta.m_pdata, fe->zet.m_pdata,
                      fe->weights.m_pdata);
          npt_x = 1;
          npt_y = 1;
          npt_z = 4;
          break;
        case 2:
          int_tetra_5(fe->ksi.m_pdata, fe->eta.m_pdata, fe->zet.m_pdata,
                      fe->weights.m_pdata);
          npt_x = 1;
          npt_y = 1;
          npt_z = 5;
          break;
        case 3:
          int_tetra_11(fe->ksi.m_pdata, fe->eta.m_pdata, fe->zet.m_pdata,
                      fe->weights.m_pdata);
          npt_x = 1;
          npt_y = 1;
          npt_z = 11;

          break;
        default:
          integrate(e_type, &npt_x, &npt_y, &npt_z,
                  fe->ksi.m_pdata, fe->eta.m_pdata, fe->zet.m_pdata,
                  fe->weights.m_pdata);
      }
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
        fe->ksi.m_pdata[ia] = gk[ia];
        fe->eta.m_pdata[ia] = gk[ia];
        fe->zet.m_pdata[ia] = gk[ia];
        fe->weights.m_pdata[ia] = w[ia];
      }
      break;
    }

  }

  Matrix_construct(int, fe->itg_ids);
  Matrix_redim(fe->itg_ids, npt_x*npt_y*npt_z, nsd);

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

  fe->ST_tensor = aloc4(3,3,nsd,nne);
  fe->ST = aloc1(3*3*nsd*nne);
}

/// set FEM libray by element
///
/// Set the FEM libray: shape function, deriviative of the shape functions, weight,
/// Gaussian points, node ids, element id, and node coordinates
/// number of node on an element and integration order determine size of member varialbes (Memory will be allocated)
/// FEMLIB_destruct should be called to cleanup this structure.
///
/// \param[out] fe self of FEMLIB
/// \param[in] e element id
/// \param[in] elem list of ELEMENT
/// \param[in] node list of NODE
/// \param[in] i_order integration order, 0: linear, 1: quadratic, and 2: higher
/// \param[in] is_total if 1: total Lagrangian, 0: updated Lagrangian
/// \return void
void FEMLIB_initialization_by_elem(FEMLIB *fe, int e, const ELEMENT *elem, const NODE *node, int i_order, int is_total)
{
  int nne = elem[e].toe;

  FEMLIB_initialization(fe, nne, i_order, nne);

  long *nod = (fe->node_id).m_pdata;  // no memory is allocated
  elemnodes(e,nne,nod,elem);

  FEMLIB_set_variable_size(&(fe->temp_v), nne,nne);

  double *x = (fe->temp_v).x.m_pdata; // no memory is allocated
  double *y = (fe->temp_v).y.m_pdata;
  double *z = (fe->temp_v).z.m_pdata;

  if(is_total)
    nodecoord_total(nne,nod,node,x,y,z);
  else
    nodecoord_updated(nne,nod,node,x,y,z);

  for(int a=0; a<nne; a++)
  {
    Mat_v(fe->node_coord, a+1, 1) = x[a];
    Mat_v(fe->node_coord, a+1, 2) = y[a];
    Mat_v(fe->node_coord, a+1, 3) = z[a];
  }

  fe->curt_elem_id = e;
}

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
                                            int is_total)
{
  int nne = elem[e].toe;
  int nne_t = nne + elem[e].n_bub;

  FEMLIB_initialization(fe, nne, i_order, nne);

  long *nod = (fe->node_id).m_pdata;  // no memory is allocated
  elemnodes(e,nne,nod,elem);

  FEMLIB_set_variable_size(&(fe->temp_v), nne_t, nne);

  double *x = (fe->temp_v).x.m_pdata; // no memory is allocated
  double *y = (fe->temp_v).y.m_pdata;
  double *z = (fe->temp_v).z.m_pdata;

  if(is_total)
    nodecoord_total(nne,nod,node,x,y,z);
  else
    nodecoord_updated(nne,nod,node,x,y,z);

  element_center(nne,x,y,z);

  for(int a=0; a<nne; a++)
  {
    Mat_v(fe->node_coord, a+1, 1) = x[a];
    Mat_v(fe->node_coord, a+1, 2) = y[a];
    Mat_v(fe->node_coord, a+1, 3) = z[a];
  }

  fe->curt_elem_id = e;
}

void FEMLIB_elem_shape_function(FEMLIB *fe, long ip, int nne, Matrix<double> *N)
{
  double ksi_, eta_, zet_;

  int itg_id_ip_1 = Mat_v(fe->itg_ids, ip, 1);
  int itg_id_ip_2 = Mat_v(fe->itg_ids, ip, 2);
  int itg_id_ip_3 = Mat_v(fe->itg_ids, ip, 3);

  // @todo Removed as dead code. @cp should review (see also comments below). LD
  // double wt;

  if(fe->elem_type == HEXAHEDRAL)
  {/* hexahedron */
    ksi_ = Mat_v(fe->ksi, itg_id_ip_1, 1);
    eta_ = Mat_v(fe->eta, itg_id_ip_2, 1);
    zet_ = Mat_v(fe->zet, itg_id_ip_3, 1);

    // wt = Mat_v(fe->weights, itg_id_ip_1, 1)*Mat_v(fe->weights, itg_id_ip_2, 1)*Mat_v(fe->weights, itg_id_ip_3, 1);
  }
  else
  { /* tetrahedron type */
    ksi_ = Mat_v(fe->ksi, itg_id_ip_3, 1);
    eta_ = Mat_v(fe->eta, itg_id_ip_3, 1);
    zet_ = Mat_v(fe->zet, itg_id_ip_3, 1);
    // wt = Mat_v(fe->weights, itg_id_ip_3, 1);
  }

  shape_func(ksi_, eta_, zet_, nne, N->m_pdata);
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

void FEMLIB_update_shape_tensor(FEMLIB *fe)
{
  shape_tensor(fe->nne,fe->nsd,fe->temp_v.N_x.m_pdata,
                               fe->temp_v.N_y.m_pdata,
                               fe->temp_v.N_z.m_pdata,fe->ST_tensor);
  shapeTensor2array(fe->ST,CONST_4(double) fe->ST_tensor,fe->nne);
}

void FEMLIB_update_deformation_gradient(FEMLIB *fe, const int ndofn, double *u, Matrix<double> *F)
{
  double **F_mat;
  F_mat = aloc2(3,3);
  def_grad_get(fe->nne,ndofn,CONST_4(double) fe->ST_tensor,u,F_mat);
  mat2array(F->m_pdata,CONST_2(double) F_mat,3,3);
  dealoc2(F_mat,3);
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
  Matrix_cleanup(fe->node_id);

  dealoc4(fe->ST_tensor,3,3,fe->nsd);
  free(fe->ST);
}
