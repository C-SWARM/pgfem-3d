/// Define energy equation: function for computing stiffness matrix and residual vector
///
/// Authors:
///  Sangmin Lee, [1], <slee43@nd.edu>
///  Aaron Howell, [1], <ahowell3@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "PGFem3D_data_structure.h"
#include "PLoc_Sparse.h"
#include "data_structure.h"
#include "constitutive_model.h"
#include "energy_equation.h"
#include "femlib.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "hyperelasticity.h"     // <= constitutive model elasticity
#include "material_properties.h" // <= constitutive model material properties
#include "utils.h"
#include <ttl/ttl.h>
#include <cmath>

namespace {
using namespace gcm;
using namespace ttl;

using pgfem3d::Solver;

const constexpr int DIM_3 = 3;
const constexpr int DIM_3x3 = 9;
const constexpr int DIM_3x3x3 = 27;
const constexpr int DIM_3x3x3x3 = 81;

const constexpr double TOL_FHS = 1.0e-6;
const constexpr double PLASTIC_HEAT_FACTOR = 0.8;

const constexpr Index<'A'> A;
const constexpr Index<'B'> B;
const constexpr Index<'C'> C;
const constexpr Index<'D'> D;
const constexpr Index<'I'> I;
const constexpr Index<'J'> J;
const constexpr Index<'K'> K;
const constexpr Index<'L'> L;
const constexpr Index<'M'> M;
const constexpr Index<'N'> N;
const constexpr Index<'O'> O;
const constexpr Index<'P'> P;
const constexpr Index<'Q'> Q;
const constexpr Index<'X'> X;
const constexpr Index<'Y'> Y;

template <class T1, class T2>
inline int inv(T1 &A, T2 &AI)
{
  int err = inv3x3(A.data, AI.data);
  return err;
}

const Tensor<2,3,double> TD = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
}

/*
/// compute derivative of PK1 w.r.t F
///
/// dePdeF(I,J,K,L) = delta(I,K)*S(L,J) + F(I,M)*C(M,J,P,Q)*dEdF(P,Q,K,L)
/// dEdF(P,Q,K,L) = delta(P,L)*F(K,Q)+F(K,P)*delta(Q,L);
///
/// \param[out] dePdeF coumputed 4th order tensor
/// \param[in] S PK2 stress
/// \param[in] dWdE elasticity tensor
/// \param[in] F deformation gradient tensor
/// \return non-zero on internal error
int compute_dePdeF(double *_dePdeF,
                   double *_S,
                   double *_dWdE,
                   double *_F)
{
  int err = 0;

  TensorA<4> dePdeF(_dePdeF);
  TensorA<2> S(_S);
  TensorA<4> dWdE(_dWdE);
  TensorA<2> F(_F);

  Tensor<4> dEdF = {};
  Tensor<2,3,double> delta;
  delta(I,J) = identity(I,J);
  dEdF(P,Q,K,L) = delta(P,L)*F(K,Q)+F(K,P)*delta(Q,L);
  dePdeF(I,J,K,L) = delta(I,K)*S(L,J) + F(I,M)*dWdE(M,J,P,Q)*dEdF(P,Q,K,L);

  return err;
}*/
/*
/// compute derivative of PK1 w.r.t F
///
/// d2PdF2(I,J,K,L,A,B) = delta(I,K)*dWdE(L,J,M,X) + dEdF(M,X,A,B)
///                    + delta(I,A)*dWdE(B,J,P,Q) + dEdF(P,Q,K,L)
///                    + F(I,M)*dCdE(M,J,P,Q,X,Y)*dEdF(X,Y,A,B)*dEdF(P,Q,K,L)
///                    + F(I,M)*dWdE(M,J,P,Q)*d2PdF2(P,Q,K,L,A,B)
/// dEdF(P,Q,K,L) = delta(P,L)*F(K,Q)+F(K,P)*delta(Q,L);
///
/// \param[out] d2PdF2 coumputed 6th order tensor
/// \param[in] S PK2 stress
/// \param[in] dWdE elasticity tensor
/// \param[in] dCdE 6th order dCdE tensor
/// \param[in] F deformation gradient tensor
/// \return non-zero on internal error
int compute_d2ePdeF2(double *_d2PdF2,
                     double *_S,
                     double *_dWdE,
                     double *_dCdE,
                     double *_F)
{
  int err = 0;

  TensorA<6, 3, double*> d2PdF2(_d2PdF2);
  TensorA<2, 3, double*> S(_S);
  TensorA<4, 3, double*> dWdE(_dWdE);
  TensorA<6, 3, double*> dCdE(_dCdE);
  TensorA<2, 3, double*> F(_F);

  Tensor<2, 3, double> delta = identity(I,J);

  Tensor<4, 3, double> dEdF = delta(P,L)*F(K,Q)+F(K,P)*delta(Q,L);

  Tensor<6, 3, double> temp_a = delta(I,K)*dWdE(L,J,M,X)*dEdF(M,X,A,B)
                 + delta(I,A)*dWdE(B,J,P,Q)*dEdF(P,Q,K,L)
                 + F(I,M)*dWdE(M,J,P,Q)*(delta(P,L)*delta(K,A)*delta(Q,B) + delta(K,A)*delta(P,B)*delta(Q,L));

  Tensor<6,3,double> temp_b = dCdE(M,J,P,Q,X,Y)*dEdF(X,Y,A,B);
  Tensor<6,3,double> temp_c = temp_b(M,J,P,Q,A,B)*dEdF(P,Q,K,L);
  Tensor<6,3,double> temp_d = F(I,M)*temp_c(M,J,A,B,K,L);

  d2PdF2(I,J,K,L,A,B) = temp_a(I,J,K,L,A,B) + temp_d(I,J,K,L,A,B);

  return err;
} */

/// compute deformation gradient due to heat expansion ttl version
///
/// \param[out] hF deformation gradient due to heat expansion
/// \param[in] dT temperature difference
/// \param[in] mat MaterialProperty object
/// \param[in] mat_id material id
/// \param[in] diff_order, if 0 deformation gradient
///                        if 1 1st order of differentiation of hF w.r.t temperature
///                        if 2 2nd order of differentiation of hF w.r.t temperature
/// \return non-zero with interal error
int compute_hF_ttl(Tensor<2, DIM_3, double> &hF,
                   double dT,
                   const MaterialProperty *mat,
                   const int mat_id,
                   const int diff_order)
{
  int err = 0.0;
  // compute thermal part of deformation gradient
  double ax = mat->mater[mat_id].ax;
  double ay = mat->mater[mat_id].ay;
  double az = mat->mater[mat_id].az;

  hF = 0.0*TD(I,J);
  switch(diff_order)
  {
    case 0:
      hF[0][0] = (1.0 + ax*dT);
      hF[1][1] = (1.0 + ay*dT);
      hF[2][2] = (1.0 + az*dT);
      break;
    case 1:
      hF[0][0] = ax*dT;
      hF[1][1] = ay*dT;
      hF[2][2] = az*dT;
      break;
    case 2:
      break;
  }
  return err;
}
/*
/// compute differentiation eF w.r.t hF
///
/// \param[out] dF computed 4th order tensor
/// \param[in] F total deformation tensor
/// \param[in] pFI inverse of the plastic part deformation gradient
/// \param[in] hFI inverse of the thermal part deformation gradient
/// \return non-zero on interal error
int compute_deF_over_dhF(Matrix<double> *dF,
                         Matrix<double> *F,
                         Matrix<double> *pFI,
                         Matrix<double> *hFI)
{
  int err = 0;
  for(int I=1; I<=DIM_3; I++)
  {
    for(int J=1; J<=DIM_3; J++)
    {
      for(int K=1; K<=DIM_3; K++)
      {
        for(int L=1; L<=DIM_3; L++)
        {
          Tns4_v(*dF,I,J,K,L) = 0.0;
          for(int M=1; M<=DIM_3; M++)
          {
            for(int N=1; N<=DIM_3; N++)
              Tns4_v(*dF,I,J,K,L) -= Mat_v(*F,I,M)*Mat_v(*hFI,M,K)*Mat_v(*hFI,L,N)*Mat_v(*pFI,N,J);
          }
        }
      }
    }
  }
  return err;
}*/
/*
/// compute differentiation eF w.r.t pF
///
/// \param[out] dF computed 4th order tensor
/// \param[in] F total deformation tensor
/// \param[in] pFI inverse of the plastic part deformation gradient
/// \param[in] hFI inverse of the thermal part deformation gradient
/// \return non-zero on interal error
int compute_deF_over_dpF(Matrix<double> *dF,
                         Matrix<double> *F,
                         Matrix<double> *pFI,
                         Matrix<double> *hFI)
{
  int err = 0;
  for(int I=1; I<=DIM_3; I++)
  {
    for(int J=1; J<=DIM_3; J++)
    {
      for(int K=1; K<=DIM_3; K++)
      {
        for(int L=1; L<=DIM_3; L++)
        {
          Tns4_v(*dF,I,J,K,L) = 0.0;
          for(int M=1; M<=DIM_3; M++)
          {
            for(int N=1; N<=DIM_3; N++)
              Tns4_v(*dF,I,J,K,L) -= Mat_v(*F,I,M)*Mat_v(*hFI,M,N)*Mat_v(*pFI,N,K)*Mat_v(*pFI,L,J);
          }
        }
      }
    }
  }
  return err;
}*/

/*
/// compute derivative of PK1 w.r.t F
///
/// d2PdF2(I,J,K,L,A,B) = delta(I,K)*dWdE(L,J,M,X) + dEdF(M,X,A,B)
///                    + delta(I,A)*dWdE(B,J,P,Q) + dEdF(P,Q,K,L)
///                    + F(I,M)*dCdE(M,J,P,Q,X,Y)*dEdF(X,Y,A,B)*dEdF(P,Q,K,L)
///                    + F(I,M)*dWdE(M,J,P,Q)*d2PdF2(P,Q,K,L,A,B)
/// dEdF(P,Q,K,L) = delta(P,L)*F(K,Q)+F(K,P)*delta(Q,L);
///
/// \param[out] d2PdF2 coumputed 6th order tensor
/// \param[in] S PK2 stress
/// \param[in] dWdE elasticity tensor
/// \param[in] dCdE 6th order dCdE tensor
/// \param[in] F deformation gradient tensor
/// \return non-zero on internal error
int compute_dPdhF(Matrix<double> *dPdhF_in,
                  Matrix<double> *deFdhF_in,
                  Matrix<double> *dePdeF_in,
                  Matrix<double> *eF_in,
                  Matrix<double> *eP_in,
                  Matrix<double> *pFI_in,
                  Matrix<double> *hFI_in)
{
  int err = 0;

  Tensor<4, 3, double*> dPdhF(dPdhF_in->m_pdata);
  Tensor<4, 3, double*> deFdhF(deFdhF_in->m_pdata);
  Tensor<4, 3, double*> dePdeF(dePdeF_in->m_pdata);
  Tensor<2, 3, double*> eF(eF_in->m_pdata);
  Tensor<2, 3, double*> eP(eP_in->m_pdata);
  Tensor<2, 3, double*> pFI(pFI_in->m_pdata);
  Tensor<2, 3, double*> hFI(hFI_in->m_pdata);

  dPdhF(I,J,K,L) = (deFdhF(I,J,K,M)*eP(M,N) + eF(I,M)*dePdeF(M,J,A,B)*deFdhF(A,B,K,N))*pFI(O,N)*hFI(L,O)
                 + eF(I,M)*eP(M,N)*pFI(O,N)*hFI(O,K)*hFI(L,J);

  return err;
}*/
/// compute heat generation due to mechanical (reference configureation)
///
/// \param[out] Qe thermal source due to mechanical work (elastic part)
/// \param[out] Qp thermal source due to mechanical work (plastic part)
/// \param[out] DQ tangent of heat generation w.r.t temperature (computed only if compute_tangent = 1)
/// \param[in] mat MaterialProperty object
/// \param[in] fv_m mechanical FieldVariables object
/// \param[in] T temperature
/// \param[in] dT temperature changes
/// \param[in] dt time step size
/// \param[in] eid element id
/// \param[in] ip integration point id
/// \param[in] mat_id material id
/// \return non-zero on internal error
int compute_mechanical_heat_gen(double *Qe,
                                double *Qp,
                                double *DQ,
                                const MaterialProperty *mat,
                                const FieldVariables *fv_m,
                                const double T,
                                const double deltaT_np1,
                                const double deltaT_n,
                                const double dt,
                                const int eid,
                                const int ip,
                                const int mat_id,
                                const int compute_tangent)
{
  int err = 0;
  int compute_stiffness = 1;

  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  double hJ = 0.0;
  double pJ = 0.0;

  // 1. compute deformation gradient of thermal expansions
  int diff_order = 0; // if 0 deformation gradient
                      // if 1 1st order of differentiation of hF w.r.t temperature
                      // if 2 2nd order of differentiation of hF w.r.t temperature
  Tensor<2,DIM_3,double> hFn,hF,hFp,hFpp,hFI,hFnI;

  err += compute_hF_ttl(hFn,  deltaT_n,   mat, mat_id, diff_order);
  err += compute_hF_ttl(hF,   deltaT_np1, mat, mat_id, diff_order);
  err += compute_hF_ttl(hFp,  deltaT_np1, mat, mat_id, 1);
  err += compute_hF_ttl(hFpp, deltaT_np1, mat, mat_id, 2);

  hFI   = inverse(hF);
  hFnI  = inverse(hFn);

  // 2. obtain deformation gradient from mechanical part
  Constitutive_model *m = &(fv_m->eps[eid].model[ip-1]);
  //Model_parameters *func = m->param;
  ELASTICITY *elast = (m->param)->cm_elast;

  Tensor<2,DIM_3,double> F,Fn,pF,pFn;

  err += m->param->get_F(  m,F.data, 1); // this brings  F(t(n+1))
  err += m->param->get_F( m,Fn.data, 0); // this brings  F(t(n))
  err += m->param->get_pF( m,pF.data,1); //             pF(t(n+1))
  err += m->param->get_pF(m,pFn.data,0); //             pF(t(n))

  // 3. compute eFs, dot_xFs, det(xFs)
  Tensor<2,DIM_3,double> eF,eFn,pFI,pFnI,deF,dpF,dhF;
  pFI  = inverse(pF);
  pFnI = inverse(pFn);

  eF  = F(I,K)*hFI(K,L)*pFI(L,J);
  eFn = F(I,K)*hFI(K,L)*pFI(L,J);

  hJ = det(hF);
  pJ = det(pF);

  deF = (eF(I,J) - eFn(I,J))/dt;
  dpF = (pF(I,J) - pFn(I,J))/dt;
  dhF = (hF(I,J) - hFn(I,J))/dt;
  double Tdot = deltaT_np1/dt;

  // 4. compute stress and elasticity
  elast->update_elasticity(elast,eF.data,compute_stiffness);

  const Tensor<2,DIM_3,double*> S{elast->S};
  const Tensor<4,DIM_3,double*> dWdE{elast->L};
  Tensor<2,DIM_3,double> dWdeF = eF(I,K)*S(K,J);

  // 5. compute heat gen of elastic part:
  *Qe = 0.0;
  *Qp = 0.0;
  *DQ = 0.0;

  Tensor<2,DIM_3,double> eFpF   = eF(I,K)*pF(K,J);
  Tensor<2,DIM_3,double> hFIpFI = hFI(I,K)*pFI(K,J);

  Tensor<4,DIM_3,double> d2WdeF2;

  Tensor<4, 3, double> dEdF = {};
  dEdF(P,Q,K,L) = TD(P,L)*F(K,Q)+F(K,P)*TD(Q,L);
  d2WdeF2(I,J,K,L) = TD(I,K)*S(L,J) + eF(I,M)*dWdE(M,J,P,Q)*dEdF(P,Q,K,L);

  Tensor<2,DIM_3,double> Qe_1 = eFpF(A,K)*hFp(K,B)*d2WdeF2(A,B,I,M)*hFIpFI(M,J); //Qe_1(I,J)
  Tensor<2,DIM_3,double> Qe_2 = dWdeF(I,K)*hFIpFI(M,K)*hFp(L,M)*pF(J,L);         //Qe_2(I,J)

  double qe1 = (Qe_1(I,J) + Qe_2(I,J))*deF(I,J);

  double factor = dWdeF(I,K)*hFIpFI(J,K)*eFpF(I,L)*hFp(L,J);
  Tensor<2,DIM_3,double> d2WdTdeF = {};
  Tensor<2,DIM_3,double> Qe_3 = Tdot*eFpF(K,I)*d2WdTdeF(K,L)*hFIpFI(J,L);
  Tensor<2,DIM_3,double> Qe_4 = -hFI(K,I)*hFp(L,K)*eFpF(M,L)*dWdeF(M,N)*hFIpFI(J,N);

  double qe2 = (factor*hFI(J,I) + Qe_3(I,J) + Qe_4(I,J))*dhF(I,J);

  double qe3 = Tdot*dWdeF(I,K)*hFIpFI(J,K)*eFpF(I,L)*hFpp(L,J);;
  *Qe = -T*hJ*pJ*(qe1 + qe2 + qe3);

  // 6. compute heat gen of plastic part:

  Tensor<2,DIM_3,double> Qp_1 = -hFIpFI(K,I)*hFp(L,K)*eFpF(M,L)*dWdeF(M,N)*pFI(J,N);
  Tensor<2,DIM_3,double> Qp_2 = eF(K,I)*dWdeF(K,L)*hFIpFI(M,L)*hFp(J,M);
  Tensor<2,DIM_3,double> Qp_3 = -eF(K,I)*dWdeF(K,L)*pFI(J,L);

  *Qp = -PLASTIC_HEAT_FACTOR*hJ*pJ*(T*factor*hFI(J,I) + T*Qp_1(I,J) + T*Qp_2(I,J) + Qp_3(I,J))*dpF(I,J);
  return err;
}

/// compute temperature on node in a element
///
/// \param[in] cn local DOF ids for field variables
/// \param[in] ndofe number of degree of freedom for element
/// \param[in] T nodal temperature for all node
/// \param[in] dT temperature increment
/// \param[in] elem Element object
/// \param[in] node Node object
/// \param[in] sup struct for BC's
/// \param[out] T_e computed nodal temperature for current element
/// \param[in] T0 reference temperature
/// \return non-zero on internal error
int get_temperature_elem(const long *cn,
                         const long ndofe,
                         const double *T,
                         const double *dT,
                         const Element *elem,
                         const Node *node,
                         const SUPP sup,
                         double *T_e,
                         double T0)
{
  int err = 0;
  for(int i=0; i< ndofe; i++){
    const int id = cn[i];
    const int aid = abs(id) - 1;

    if (id == 0){
      T_e[i] = T0;
    } else if (id > 0){
      T_e[i] = T[aid] + dT[aid];
    } else {
      T_e[i] = T0 + sup->defl[aid] + sup->defl_d[aid];
    }
  }
  return err;
}

/// assemble residuals for heat conduction problem
///
/// element-wize residuals are merged into global residual vector
///
/// \param[in] fe container of finite element resources
/// \param[in] fi local residual(element-wize)
/// \param[in] grid a mesh object
/// \param[in,out] fv field variable object
/// \param[in] mp_id mutiphysics id
/// \return non-zero on internal error
int energy_equation_residuals_assemble(FEMLIB *fe,
                                       double *fi,
                                       Grid *grid,
                                       FieldVariables *fv,
                                       const int mp_id)
{
  int err = 0;
  long *nod = fe->node_id.m_pdata;

  for(int ia = 0; ia<fe->nne; ia++)
  {
    for(int ib = 0; ib<fv->ndofn; ib++)
    {
      int II = grid->node[nod[ia]].id_map[mp_id].id[ib] - 1;
      if (II < 0) continue;
      fv->f_u[II] += fi[ia*(fv->ndofn) + ib];
    }
  }
  return err;
}

/// compute residuals for heat conduction problem at the element level
///
/// Actual computation for constructing residual vector takes place here.
/// Memory of the residual vector should be allocated before this function calls.
///
/// \param[in] fe container of finite element resources
/// \param[out] fi_in local residual vector (element level)
/// \param[in] du temperature increment
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in] FV array of field variable object
/// \param[in] load object for loading
/// \param[in] mp_id mutiphysics id
/// \param[in] dt_in time step size
/// \param[in] total_Lagrangian flag denotes total_Lagrangian= 1: Total Lagrangian formulation
///                                          total_Lagrangian= 0: Updated Lagrangian formulation
/// \return non-zero on internal error
int energy_equation_compute_residuals_elem(FEMLIB *fe,
                                           double *fi,
                                           double *du,
                                           Grid *grid,
                                           MaterialProperty *mat,
                                           FieldVariables *fv,
                                           LoadingSteps *load,
                                           int mp_id,
                                           double dt_in,
                                           int total_Lagrangian)

{
  int err = 0;
  double Jn = 1.0; // Total Lagrangian = 1.0; // it is updated if mechanical is coupled.

  int eid = fe->curt_elem_id;
  const int mat_id = (grid->element[eid]).mat[0];
  double rho_0 = mat->density[mat_id];

  int is_it_couple_w_mechanical  = -1;
  int is_it_couple_w_chemical    = -1;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_MECHANICAL)
      is_it_couple_w_mechanical = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  if(is_it_couple_w_chemical>=0)
  {
    // compute for chemical
  }

  const MaterialThermal&     thermal = mat->thermal[mat_id];
  const Tensor<2,3,const double*> k0 = thermal.k;
  Tensor<2,3,double>               k = k0(I,J);

  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  double cp = thermal.cp;

  double dt = 1.0; // for the quasi steady state
  if(rho_0>0)
    dt = dt_in;

  SUPP sup = load->sups[mp_id];

  long *nod = fe->node_id.m_pdata; // use only address, no need deallication
  int ndofn = fv->ndofn;
  long ndofe = (fe->nne)*ndofn;

  Matrix<long> cnL(ndofe, 1);

  get_dof_ids_on_elem_nodes(0,fe->nne,ndofn,nod,grid->node,cnL.m_pdata,mp_id);

  Matrix<double> q(grid->nsd,1), Tnp1(fe->nne,1), Tn(fe->nne,1);

  //compute nodal value
  get_temperature_elem(cnL.m_pdata,ndofe,fv->u_np1,du,grid->element,grid->node,sup,Tnp1.m_pdata,fv->u0);

  for(int ia=0; ia<fe->nne; ia++)
    Tn(ia+1) = fv->u_n[nod[ia]];

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    double Q = 0.0;
    // Udate basis functions at the integration points.
    fe->elem_basis_V(ip);
    fe->update_shape_tensor();

    double Temp = 0.0;
    double dT   = 0.0;
    double deltaT_np1 = 0.0;
    double deltaT_n   = 0.0;

    // compute varialbes at the integration point
    for(int ia=1; ia<=fe->nne; ia++)
    {
      Temp += fe->N(ia)*Tnp1(ia);
      dT   += fe->N(ia)*(Tnp1(ia)-Tn(ia));
      deltaT_np1 += fe->N(ia)*(Tnp1(ia) - fv->u0);
      deltaT_n   += fe->N(ia)*(Tn(ia) - fv->u0);
    }

    if(is_it_couple_w_mechanical>=0)
    {
      FieldVariables *fv_m = fv->fvs[is_it_couple_w_mechanical];
      int compute_tangent = 0;
      double DQ = 0.0;
      double Qe = 0.0;
      double Qp = 0.0;

      Tensor<2,3,double> F;
      Constitutive_model *m = &(fv_m->eps[eid].model[ip-1]);
      const Model_parameters *func = m->param;
      err += func->get_F(m,F.data,1);        // this brings  F(t(n+1))

      if(total_Lagrangian)
      {
        double detF = det(F);

        Tensor<2,3,double> FI;
        inv(F,FI);
        k = detF * FI(I,K) * k0(K,L) * FI(J,L); // k = J*FI*k0*FIT
      }
      else
        Jn = det(F); //updated Lagrangian

      if (thermal.FHS_MW > TOL_FHS)
      {
        err += compute_mechanical_heat_gen(&Qe,&Qp,&DQ,mat,fv_m,Temp,deltaT_np1,deltaT_np1,dt,eid,ip,mat_id,compute_tangent);
        Q += thermal.FHS_MW * (Qe + Qp);

      }
    }

    q.set_values(0.0);

    // compute heat flux
    for(int ia=1; ia<=fe->nne; ia++)
    {
      // k = [nsd, nsd], dN = [nne, nsd], Tnp1  = [nne, 1],
      // q = [nsd, 1] = k*dN'*T
      for(int ib=1; ib<=grid->nsd; ib++)
      {
        for(int ic=1; ic<=grid->nsd; ic++) {
          q(ib) += k(ib-1,ic-1)      // NB: TTL does not support 1-based indices
                   *fe->dN(ia,ic)*Tnp1(ia);
        }
      }
    }

    // R = rho*cp*dT + dt*grad.q - dt*Q = 0;
    // rho = rho_0/Jn
    // Q = Q_0/Jn;
    for(int ia=1; ia<=fe->nne; ia++)
    {
      fi[ia-1] += fe->N(ia)*(rho_0*cp*dT - dt*Q)*(fe->detJxW)/Jn;
      for(int ib=1; ib<=grid->nsd; ib++)
        fi[ia-1] += dt*fe->dN(ia,ib)*q(ib)*(fe->detJxW)/Jn;
    }
  }

  return err;
}

/// compute residuals for heat conduction problem
///
/// Every element will be visited to compute residuals and local residuals are
/// assembled to global residual vector. The actual computation for constructing
/// residuals takes place in energy_equation_compute_residuals_elem function.
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv field variable object
/// \param[in] load object for loading
/// \param[in] mp_id mutiphysics id
/// \param[in] use_updated if use_updated=1, compute residuals updated temperature
///                           use_updated=0, compute residuals using temporal temperature
/// \param[in] dt time step size
/// \return non-zero on internal error
int energy_equation_compute_residuals(Grid *grid,
                                      MaterialProperty *mat,
                                      FieldVariables *fv,
                                      LoadingSteps *load,
                                      const int mp_id,
                                      int use_updated,
                                      double dt)
{
  int err = 0;
  int total_Lagrangian = 1;
  int intg_order = 0;

  int is_it_couple_w_mechanical  = -1;
  int is_it_couple_w_chemical    = -1;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_MECHANICAL)
      is_it_couple_w_mechanical = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  if(is_it_couple_w_chemical>=0)
  {
    // compute for chemical
  }

  // save original pointer to access mechanical part
  double *u_n = NULL;
  double *u_nm1 = NULL;
  State_variables *statv_list = NULL;

  if(is_it_couple_w_mechanical>=0)
  {
    FieldVariables *fv_m = fv->fvs[is_it_couple_w_mechanical];
    u_n   = fv_m->u_n;
    u_nm1 = fv_m->u_nm1;
    statv_list = fv_m->statv_list;

    fv_m->u_n   = fv_m->temporal->u_n;
    fv_m->u_nm1 = fv_m->temporal->u_nm1;
    fv_m->statv_list = fv_m->temporal->var;
  }

  double *du;
  if(use_updated)
    du = fv->f;
  else
    du = fv->d_u;

  for(int eid=0; eid<grid->ne; eid++)
  {
    // Construct finite element library.
    // It provide element wise integration info
    // such as basis function, weights, ...
    FEMLIB fe(eid,grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element,
    // fe needs to be updated by integration points
    Matrix<double> fi(fe.nne, 1, 0.0);
    err += energy_equation_compute_residuals_elem(&fe,fi.m_pdata,du,grid,mat,fv,load,mp_id,dt,total_Lagrangian);
    err += energy_equation_residuals_assemble(&fe,fi.m_pdata,grid,fv,mp_id);
  }

  if(is_it_couple_w_mechanical>=0)
  {
    FieldVariables *fv_m = fv->fvs[is_it_couple_w_mechanical];
    fv_m->u_n   = u_n;
    fv_m->u_nm1 = u_nm1;
    fv_m->statv_list = statv_list;
  }
  return err;
}

/// compute stiffness for heat conduction problem at the element level
///
/// Actual computation for constructing stiffness matrix takes place here.
/// Local memory for the element level stiffness matrix is created and
/// merged into memory for assembling and communications.
///
/// \param[in] fe container of finite element resources
/// \param[out] Lk local stiffness matrix for assmebling and communication
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in] fv field variable object
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com object for communications
/// \param[in] Ddof number of degree of freedoms accumulated through processes
/// \param[in] myrank current process rank
/// \param[in] interior indentifier to distinguish element in interior or
///             communication boundary
/// \param[in] opts structure of PGFem3D options
/// \param[in] mp_id mutiphysics id
/// \param[in] do_assemble if yes, assmeble computed element stiffness matrix into
///            sparce solver matrix or compute element stiffness only
/// \param[in] compute_load4pBCs, if yes, compute external flux due to Dirichlet BCs
///                               if no, no compute external flux due to Dirichlet BCs
/// \param[in] dt_in time step size
/// \param[in] total_Lagrangian flag denotes total_Lagrangian= 1: Total Lagrangian formulation
///                                          total_Lagrangian= 0: Updated Lagrangian formulation
/// \return non-zero on internal error
int energy_equation_compute_stiffness_elem(FEMLIB *fe,
                                           double **Lk,
                                           Grid *grid,
                                           MaterialProperty *mat,
                                           FieldVariables *fv,
                                           Solver *sol,
                                           LoadingSteps *load,
                                           CommunicationStructure *com,
                                           int *Ddof,
                                           int myrank,
                                           int interior,
                                           const PGFem3D_opt *opts,
                                           int mp_id,
                                           int do_assemble,
                                           int compute_load4pBCs,
                                           double dt_in,
                                           int total_Lagrangian)

{
  int err = 0;
  int eid = fe->curt_elem_id;
  double Jn = 1.0; // Total Lagrangian = 1.0; // it is updated if mechanical is coupled.

  const int mat_id = (grid->element[eid]).mat[0];
  double rho_0 = mat->density[mat_id];

  int is_it_couple_w_mechanical  = -1;
  int is_it_couple_w_chemical    = -1;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_MECHANICAL)
      is_it_couple_w_mechanical = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  if(is_it_couple_w_chemical>=0)
  {
    // compute for chemical
  }

  const MaterialThermal& thermal = mat->thermal[mat_id];
  const Tensor<2,3,const double*> k0 = thermal.k;
  Tensor<2,3,double> k = k0(I,J);

  double cp = thermal.cp;

  double dt = 1.0; // for the quasi steady state
  if(rho_0>0)
    dt = dt_in;

  SUPP sup = load->sups[mp_id];

  long *nod = fe->node_id.m_pdata;
  int ndofn = fv->ndofn;
  long ndofe = (fe->nne)*ndofn;
  Matrix<long> cnL(ndofe, 1);
  Matrix<long> cnG(ndofe, 1);

  get_dof_ids_on_elem_nodes(0,fe->nne,ndofn,nod,grid->node,cnL.m_pdata,mp_id);
  get_dof_ids_on_elem_nodes(1,fe->nne,ndofn,nod,grid->node,cnG.m_pdata,mp_id);

  Matrix<double> lk(ndofe,ndofe,0.0);
  Matrix<double> Tnp1(fe->nne,1,0.0), Tn(fe->nne,1,0.0);

  //compute nodal value
  get_temperature_elem(cnL.m_pdata,ndofe,fv->u_np1,fv->f,grid->element,grid->node,sup,Tnp1.m_pdata,fv->u0);

  for(int ia=0; ia<fe->nne; ia++)
    Tn(ia+1) = fv->u_n[nod[ia]];

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    // Udate basis functions at the integration points.
    fe->elem_basis_V(ip);
    fe->update_shape_tensor();

    double Temp = 0.0;
    double dT   = 0.0;
    double deltaT_np1 = 0.0;
    double deltaT_n   = 0.0;

    // compute varialbes at the integration point
    for(int ia=1; ia<=fe->nne; ia++)
    {
      Temp += fe->N(ia)*Tnp1(ia);
      dT   += fe->N(ia)*(Tnp1(ia)-Tn(ia));
      deltaT_np1 += fe->N(ia)*(Tnp1(ia) - fv->u0);
      deltaT_n   += fe->N(ia)*(Tn(ia) - fv->u0);
    }

    double DQ = 0.0;

    if(is_it_couple_w_mechanical>=0)
    {
      FieldVariables *fv_m = fv->fvs[is_it_couple_w_mechanical];
      int compute_tangent = 1;
      double Qe = 0.0;
      double Qp = 0.0;

      Tensor<2,3,double> F = {};
      Constitutive_model *m = &(fv_m->eps[eid].model[ip-1]);
      const Model_parameters *func = m->param;
      err += func->get_F(m,F.data,1);        // this brings  F(t(n+1))

      if(total_Lagrangian)
      {
        double detF = det(F);

        Tensor<2,3,double> FI;
        inv(F,FI);
        k = detF*FI(I,K)*k0(K,L)*FI(J,L); // k = J*FI*k0*FIT
      }
      else
        Jn = det(F); //updated Lagrangian

      if (thermal.FHS_MW > TOL_FHS)
      {
        err += compute_mechanical_heat_gen(&Qe,&Qp,&DQ,mat,fv_m,Temp,deltaT_np1,deltaT_n,dt,eid,ip,mat_id,compute_tangent);
        DQ = thermal.FHS_MW * DQ;
      }
    }

    // R = rho_0*cp*(Tnp1-Tn) + dt*grad.q - dt*Q = 0;
    // DR = rho_0*cp*dT + dt*D[grad.q]dT - dt*D[Q]dT = 0;
    for(int ia=1; ia<=fe->nne; ia++)
    {
      for(int ib=1; ib<=fe->nne; ib++)
      {
        lk(ia,ib) += (rho_0*cp - dt*DQ)*fe->N(ia)*fe->N(ib)*(fe->detJxW)/Jn;
        for(int im = 1; im<=grid->nsd; im++)
        {
          for(int in = 1; in<=grid->nsd; in++)
            lk(ia,ib) += dt*fe->dN(ia,in)*
                         k(in-1, im-1) // NB: TTL does not support 1-based indices
                         *fe->dN(ib,im)*(fe->detJxW)/Jn;
        }
      }
    }
  }

  // Assemble
  if(do_assemble)
  {
    PLoc_Sparse(Lk,lk.m_pdata,
                com->Ai,
                com->Ap,
                cnL.m_pdata,cnG.m_pdata,ndofe,Ddof,
                com->GDof,
                myrank,
                com->nproc,
                com->comm,
                interior,
                sol->system,
                opts->analysis_type);
  }

  if(compute_load4pBCs)
  {
    int ndofn = 1;
    Matrix<double> u((fe->nne)*ndofn,1,0.0), f_loc((fe->nne)*ndofn,1,0.0);

    // get the bc increment
    for(int ia=0; ia<fe->nne; ia++)
    {
      for(int ib=0; ib<ndofn; ib++)
      {
        int id = ia*ndofn + ib;
        if(cnL.m_pdata[id] <= -1)
          u.m_pdata[id] = load->sups[mp_id]->defl_d[abs(cnL.m_pdata[id])-1];
        else
          u.m_pdata[id] = 0.0;
      }
    }
    f_loc.prod(lk,u);

    // element -> localization
    for(int ia=0; ia<fe->nne; ia++)
    {
      for(int ib=0; ib<ndofn; ib++)
      {
        int id_e = ia*ndofn + ib;
        int id_l = grid->node[nod[ia]].id_map[mp_id].id[ib]-1;
        if (id_l < 0)  continue;
          fv->f_defl[id_l] += f_loc.m_pdata[id_e];
      }
    }
  }

  return err;
}

/// compute stiffness for heat conduction problem
///
/// In building global stiffness matrix, elements on communcation boundary are
/// first computed and interior elements are visited later. As soon as stiffness is
/// computed on the communcation boundary, communication is established and the interior stiffness
/// is commputed. In this way, the communcation and the computation are overlaid.
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv field variable object
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] com object for communications
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int energy_equation_compute_stiffness(Grid *grid,
                                      MaterialProperty *mat,
                                      FieldVariables *fv,
                                      Solver *sol,
                                      LoadingSteps *load,
                                      CommunicationStructure *com,
                                      MPI_Comm mpi_comm,
                                      int myrank,
                                      const PGFem3D_opt *opts,
                                      const int mp_id,
                                      double dt)
{
  int err = 0;
  int total_Lagrangian  = 1;
  int intg_order        = 0;
  int do_assemble       = 1; // udate stiffness matrix
  int compute_load4pBCs = 0; // if 1, compute load due to Dirichlet BCs.
                             // if 0, update only stiffness

  int is_it_couple_w_mechanical  = -1;
  int is_it_couple_w_chemical    = -1;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_MECHANICAL)
      is_it_couple_w_mechanical = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  if(is_it_couple_w_chemical>=0)
  {
    // compute for chemical
  }

  // save original pointer to access mechanical part
  double *u_n = NULL;
  double *u_nm1 = NULL;
  State_variables *statv_list = NULL;

  if(is_it_couple_w_mechanical>=0)
  {
    FieldVariables *fv_m = fv->fvs[is_it_couple_w_mechanical];
    u_n   = fv_m->u_n;
    u_nm1 = fv_m->u_nm1;
    statv_list = fv_m->statv_list;

    fv_m->u_n   = fv_m->temporal->u_n;
    fv_m->u_nm1 = fv_m->temporal->u_nm1;
    fv_m->statv_list = fv_m->temporal->var;
  }

  double **Lk,**recieve;
  MPI_Status *sta_s,*sta_r;
  MPI_Request *req_s,*req_r;

  err += init_and_post_stiffmat_comm(&Lk,&recieve,&req_r,&sta_r,
                                     mpi_comm,com->comm);

  Matrix<int> Ddof(com->nproc,1);

  Ddof.m_pdata[0] = com->DomDof[0];
  for (int ia=1; ia<com->nproc; ia++)
    Ddof.m_pdata[ia] = Ddof.m_pdata[ia-1] + com->DomDof[ia];

  for(int eid=0; eid<com->nbndel; eid++)
  {
    // construct finite element library
    // it provide element wise integration info
    // such as basis function, weights, ...
    FEMLIB fe(com->bndel[eid],grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    int interior = 0;
    err += energy_equation_compute_stiffness_elem(&fe,Lk,grid,mat,fv,sol,load,com,Ddof.m_pdata,
                                                  myrank,interior,opts,mp_id,
                                                  do_assemble,compute_load4pBCs,dt,total_Lagrangian);

    if(err != 0)
      break;
  }
  err += send_stiffmat_comm(&sta_s,&req_s,Lk,mpi_comm,com->comm);

  int skip = 0;
  int idx  = 0;

  for(int eid=0; eid<grid->ne; eid++)
  {
    int is_it_in = is_element_interior(eid,&idx,&skip,com->nbndel,com->bndel,myrank);

    if(is_it_in==-1)
    {
      err = 1;
      break;
    }

    if(is_it_in==0)
      continue;

    FEMLIB fe(eid,grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    int interior = 1;

    err += energy_equation_compute_stiffness_elem(&fe,Lk,grid,mat,fv,sol,load,com,Ddof.m_pdata,
                                                  myrank,interior,opts,mp_id,
                                                  do_assemble,compute_load4pBCs,dt,total_Lagrangian);

    if(err != 0)
      break;
  }

  err += assemble_nonlocal_stiffmat(com->comm,sta_r,req_r,sol->system,recieve);
  err += finalize_stiffmat_comm(sta_s,sta_r,req_s,req_r,com->comm);

  if(is_it_couple_w_mechanical>=0)
  {
    FieldVariables *fv_m = fv->fvs[is_it_couple_w_mechanical];
    fv_m->u_n   = u_n;
    fv_m->u_nm1 = u_nm1;
    fv_m->statv_list = statv_list;
  }

  // stiffnes build is completed
  // deallocate memory
  for(int ia=0; ia<com->nproc; ia++)
  {
    free (recieve[ia]);
    free(Lk[ia]);
  }
  free (recieve);
  free (Lk);
  free (sta_s);
  free (sta_r);
  free (req_s);
  free (req_r);

  return err;
}

/// compute flux due to Dirichlet BCs
///
/// Compute flux vector for prescribed BCs(Dirichlet)
/// This compute load, f, as below:
/// [Kii Kio]<ui>   <bi>
/// [Koi Koo]<uo> = <bo>
/// [Kii][ui] = <bi> - [Kio]<uo>
/// where f = [Kio]<uo>, uo is Drichlet BCs
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv field variable object
/// \param[in] sol object for solution scheme
/// \param[in] load object for loading
/// \param[in] mpi_comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \param[in] opts structure PGFem3D option
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \return non-zero on internal error
int energy_equation_compute_load4pBCs(Grid *grid,
                                      MaterialProperty *mat,
                                      FieldVariables *fv,
                                      Solver *sol,
                                      LoadingSteps *load,
                                      int myrank,
                                      const PGFem3D_opt *opts,
                                      const int mp_id,
                                      double dt)
{
  int err = 0;
  int total_Lagrangian = 1;
  int intg_order       = 0;
  int interior         = 1;
  int do_assemble      = 0;
  int compute_load4pBCs= 1;

  int is_it_couple_w_mechanical  = -1;
  int is_it_couple_w_chemical    = -1;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_MECHANICAL)
      is_it_couple_w_mechanical = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  if(is_it_couple_w_chemical>=0)
  {
    // compute for chemical
  }

  // save original pointer to access mechanical part
  double *u_n = NULL;
  double *u_nm1 = NULL;
  State_variables *statv_list = NULL;

  if(is_it_couple_w_mechanical>=0)
  {
    FieldVariables *fv_m = fv->fvs[is_it_couple_w_mechanical];
    u_n   = fv_m->u_n;
    u_nm1 = fv_m->u_nm1;
    statv_list = fv_m->statv_list;

    fv_m->u_n   = fv_m->temporal->u_n;
    fv_m->u_nm1 = fv_m->temporal->u_nm1;
    fv_m->statv_list = fv_m->temporal->var;
  }

  for(int ia=0; ia<load->sups[mp_id]->nde; ia++)
  {
    int eid = load->sups[mp_id]->lepd[ia];
    FEMLIB fe(eid,grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    err += energy_equation_compute_stiffness_elem(&fe,NULL,grid,mat,fv,sol,load,NULL,NULL,
                                                  myrank,interior,opts,mp_id,
                                                  do_assemble,compute_load4pBCs,dt,total_Lagrangian);
    if(err != 0)
      break;
  }

  if(is_it_couple_w_mechanical>=0)
  {
    FieldVariables *fv_m = fv->fvs[is_it_couple_w_mechanical];
    fv_m->u_n   = u_n;
    fv_m->u_nm1 = u_nm1;
    fv_m->statv_list = statv_list;
  }
  return err;
}

/// update for for print
///
/// compute and store values (e.g flux) for print
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv a field variable object
/// \param[in] dt time step size
/// \return non-zero on internal error
int update_thermal_flux4print(Grid *grid,
                              MaterialProperty *mat,
                              FieldVariables *fv,
                              double dt)
{
  int err = 0;
  int total_Lagrangian = 1;
  int intg_order = 0;

  EPS *eps = fv->eps;

  int is_it_couple_w_mechanical  = -1;
  int is_it_couple_w_chemical    = -1;

  for(int ia=0; ia<fv->n_coupled; ia++)
  {
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_MECHANICAL)
      is_it_couple_w_mechanical = ia;
    if(fv->coupled_physics_ids[ia] == MULTIPHYSICS_CHEMICAL)
      is_it_couple_w_chemical = ia;
  }

  if(is_it_couple_w_chemical>=0)
  {
    // compute for chemical
  }

  // save original pointer to access mechanical part
  double *u_n = NULL;
  double *u_nm1 = NULL;
  State_variables *statv_list = NULL;

  if(is_it_couple_w_mechanical>=0)
  {
    FieldVariables *fv_m = fv->fvs[is_it_couple_w_mechanical];
    u_n   = fv_m->u_n;
    u_nm1 = fv_m->u_nm1;
    statv_list = fv_m->statv_list;

    fv_m->u_n   = fv_m->temporal->u_n;
    fv_m->u_nm1 = fv_m->temporal->u_nm1;
    fv_m->statv_list = fv_m->temporal->var;
  }

  if(is_it_couple_w_chemical>=0)
  {
    // compute for chemical
  }

  int myrank = 0;
  MPI_Comm_rank (MPI_COMM_WORLD,&myrank);

  for(int eid=0; eid<grid->ne; eid++)
  {
    FEMLIB fe(eid,grid->element,grid->node,intg_order,total_Lagrangian);
    memset(eps[eid].el.o,0,6*sizeof(double));

    // get material constants (parameters)
    const int mat_id = (grid->element[eid]).mat[0];

    const MaterialThermal& thermal = mat->thermal[mat_id];
    const Tensor<2,3,const double*> k0 = thermal.k;
    Tensor<2,3,double> k = k0(I,J);

    // compute noda values
    Matrix<double> q(grid->nsd,1), Tnp1(fe.nne,1), Tn(fe.nne,1);

    long *nod = fe.node_id.m_pdata;
    for(int ia=0; ia<fe.nne; ia++)
    {
      Tnp1.m_pdata[ia] = fv->u_n[nod[ia]];
      Tn.m_pdata[ia] = fv->u_nm1[nod[ia]];
    }

    // do integration point loop
    double volume = 0.0;
    for(int ip = 1; ip<=fe.nint; ip++)
    {
      fe.elem_basis_V(ip);
      fe.update_shape_tensor();

      double Temp = 0.0;
      double dT   = 0.0;
      double deltaT_np1 = 0.0;
      double deltaT_n   = 0.0;

      // compute values at the integration point
      for(int ia=1; ia<=fe.nne; ia++)
      {
        Temp += fe.N(ia,1)*Tnp1(ia,1);
        dT   += fe.N(ia,1)*(Tnp1(ia,1)-Tn(ia,1));
        deltaT_np1 += fe.N(ia,1)*(Tnp1(ia,1) - fv->u0);
        deltaT_n   += fe.N(ia,1)*(Tn(  ia,1) - fv->u0);
      }

      // compute heat sources
      double Qe = 0.0;
      double Qp = 0.0;
      double DQ = 0.0;
      double Jn = 1.0;
      if(is_it_couple_w_mechanical>=0)
      {
        FieldVariables *fv_m = fv->fvs[is_it_couple_w_mechanical];
        int compute_tangent = 0;

        Tensor<2,3,double> F;
        Constitutive_model *m = &(fv_m->eps[eid].model[ip-1]);
        const Model_parameters *func = m->param;
        err += func->get_F(m,F.data,1);        // this brings  F(t(n+1))

        if(total_Lagrangian)
        {
          double detF = det(F);

          Tensor<2,3,double> FI;
          inv(F,FI);
          k = detF*FI(I,K)*k0(K,L)*FI(J,L); // k = J*FI*k0*FIT
        }
        else
          Jn = det(F); //updated Lagrangian

        if(thermal.FHS_MW > TOL_FHS)
        {
          err += compute_mechanical_heat_gen(&Qe,&Qp,&DQ,mat,fv_m,Temp,deltaT_np1,deltaT_n,dt,eid,ip,mat_id,compute_tangent);
          Qe = thermal.FHS_MW * Qe;
          Qp = thermal.FHS_MW * Qp;
        }
      }

      q.set_values(0.0);

      // compute heat flux
      for(int ia=1; ia<=fe.nne; ia++)
      {
        // k = [nsd, nsd], dN = [nne, nsd], Tnp1  = [nne, 1],
        // q = [nsd, 1] = k*dN'*T
        for(int ib=1; ib<=grid->nsd; ib++)
        {
          for(int ic=1; ic<=grid->nsd; ic++) {
            q(ib) -= k(ib-1, ic-1)   // NB: TTL does not support 1-based indices
                     *fe.dN(ia,ic)*Tnp1(ia);
          }
        }
      }

      // save the values
      fv->eps[eid].el.o[0] += fe.detJxW*q(1);
      fv->eps[eid].el.o[1] += fe.detJxW*q(2);
      fv->eps[eid].el.o[2] += fe.detJxW*q(3);
      fv->eps[eid].el.o[3] += fe.detJxW*Qe/Jn;
      fv->eps[eid].el.o[4] += fe.detJxW*Qp/Jn;
      volume += fe.detJxW;
    }
    for(int ia=0; ia<6; ia++)
      eps[eid].el.o[ia] = eps[eid].el.o[ia]/volume;

  }

  if(is_it_couple_w_mechanical>=0)
  {
    FieldVariables *fv_m = fv->fvs[is_it_couple_w_mechanical];
    fv_m->u_n   = u_n;
    fv_m->u_nm1 = u_nm1;
    fv_m->statv_list = statv_list;
  }

  return err;
}

