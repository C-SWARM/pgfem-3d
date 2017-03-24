/// Define energy equation: function for computing stiffness matrix and residual vector
/// 
/// Authors:
///  Sangmin Lee, [1], <slee43@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN

#include "energy_equation.h"
#include "femlib.h"
#include "PGFem3D_data_structure.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "PLoc_Sparse.h"
#include <math.h>
#include "constitutive_model.h"
#include "utils.h"
#include "material_properties.h" // <= constitutive model material properties
#include "hyperelasticity.h"     // <= constitutive model elasticity

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

#define TOL_FHS 1.0e-6
#define PLASTIC_HEAT_FACTOR 0.8

#define Tns6_v(p, I,J,K,L,M,N) (p).m_pdata[DIM_3x3x3x3*3*(I-1)+DIM_3x3x3x3*(J-1)+DIM_3x3x3*(K-1)+DIM_3x3*(L-1)+DIM_3*(M-1)+(N-1)]

/// compute effective rate strain (Von Mises strain)
///
/// I = delta_ij; ed = e-tr(e)/3*I; eff = sqrt(2/3*ed:ed)
///
/// \param[out] eff computed effective strain
/// \param[in] F_in deformation gradient at t(n+1)
/// \param[in] Fn_in deformation gradient at t(n)
/// \param[in] dt time step size  
/// \return non-zero on internal error 
int compute_effective_dot_strain(double *eff, 
                                 double *F_in, 
                                 double *Fn_in, 
                                 double dt)
{
  int err = 0;
  
  Matrix(double) b, b_I, bn_I, dot_e, F;
  Matrix_construct_init(double, b,    DIM_3,DIM_3,0.0);
  Matrix_construct_init(double, bn_I, DIM_3,DIM_3,0.0);
  Matrix_construct_init(double, b_I,  DIM_3,DIM_3,0.0);
    
  Matrix_construct_init(double, dot_e,  DIM_3,DIM_3,0.0);
    
  F.m_row = F.m_col = DIM_3;
  F.m_pdata = F_in;
    
  Matrix_AxB(b,1.0,0.0,F,0,F,1);
  Matrix_inv(b,b_I);


  F.m_pdata = Fn_in;
  Matrix_AxB(b,1.0,0.0,F,0,F,1);
  Matrix_inv(b,bn_I);
    
  for(int ia=1; ia<=DIM_3; ia++)
  {
    for(int ib=1; ib<=DIM_3; ib++)
      Mat_v(dot_e,ia,ib) = -0.5*(Mat_v(b_I,ia,ib) - Mat_v(bn_I,ia,ib))/dt;
  }
  
  double tr = (Mat_v(dot_e,1,1) + Mat_v(dot_e,2,2) + Mat_v(dot_e,3,3))/3.0;

  Mat_v(dot_e,1,1) -= tr;
  Mat_v(dot_e,2,2) -= tr;
  Mat_v(dot_e,3,3) -= tr;
  
  double e_dd_e = 0.0;
  for(int ia=0; ia<DIM_3x3; ia++)
    e_dd_e += dot_e.m_pdata[ia]*dot_e.m_pdata[ia];
  
  *eff = sqrt(2.0/3.0*e_dd_e);
  
  Matrix_cleanup(b);
  Matrix_cleanup(b_I);
  Matrix_cleanup(bn_I);
  Matrix_cleanup(dot_e);
  
  return err;
}

/// compute derivative of PK1 w.r.t F
///
/// dPdF(I,J,K,L) = delta(I,K)*S(L,J) + F(I,M)*C(M,J,P,Q)*dEdF(P,Q,K,L)
/// dEdF(P,Q,K,L) = delt(P,L)*F(K,Q)+F(K,P)*delta(Q,L);
///
/// \param[out] dPdF coumputed 4th order tensor
/// \param[in] S PK2 stress
/// \param[in] dWdE elasticity tensor
/// \param[in] F deformation gradient tensor
/// \return non-zero on internal error
int compute_dPdF(Matrix(double) *dPdF,
                 Matrix(double) *S,
                 Matrix(double) *dWdE,
                 Matrix(double) *F)
{
  int err = 0;
  Matrix(double) delta;
  Matrix_construct_redim(double, delta,DIM_3,DIM_3);
  Matrix_eye(delta,DIM_3);
  
  for(int I=1; I<=DIM_3; I++)
  {
    for(int J=1; J<=DIM_3; J++)
    {
      for(int K=1; K<=DIM_3; K++)
      {
        for(int L=1; L<=DIM_3; L++)
        {
          Tns4_v(*dPdF,I,J,K,L) = 0.0;
          for(int M=1; M<=DIM_3; M++)
          {
            Tns4_v(*dPdF,I,J,K,L) += Mat_v(delta,I,K)*Mat_v(*S,M,J);
            for(int P=1; P<=DIM_3; P++)
            {
              for(int Q=1; Q<=DIM_3; Q++)
              {
                Tns4_v(*dPdF,I,J,K,L) += Mat_v(*F,I,M)*Tns4_v(*dWdE,M,J,P,Q)*
                                         (Mat_v(delta,P,L)*Mat_v(*F,K,Q)+Mat_v(*F,K,P)*Mat_v(*F,Q,L));
              }
            }
          }                  
        }
      }
    }
  }
  Matrix_cleanup(delta);
  return err;
}                 


/// compute derivative of PK1 w.r.t F
///
/// d2PdF2(I,J,K,L,AB) = delta(I,K)*dWdE(L,J,M,X) + dEdF(M,X,A,B)
///                    + delta(I,A)*dWdE(B,J,P,Q) + dEdF(P,Q,K,L)
///                    + F(I,M)*dCdE(M,J,P,Q,X,Y)*dEdF(X,Y,A,B)*dEdF(P,Q,K,L)
///                    + F(I,M)*dWdE(M,J,P,Q)*d2EdF2(P,Q,K,L,A,B)
/// dEdF(P,Q,K,L) = delt(P,L)*F(K,Q)+F(K,P)*delta(Q,L);
///
/// \param[out] d2PdF2 coumputed 6th order tensor
/// \param[in] S PK2 stress
/// \param[in] dWdE elasticity tensor
/// \param[in] dCdE 6th order dCdE tensor
/// \param[in] F deformation gradient tensor
/// \return non-zero on internal error
int compute_d2PdF2(Matrix(double) *d2PdF2,
                   Matrix(double) *S,
                   Matrix(double) *dWdE,
                   Matrix(double) *dCdE,
                   Matrix(double) *F)
{
  int err = 0;
  Matrix(double) delta, dEdF;
  Matrix_construct_redim(double, delta,DIM_3,DIM_3);
  Matrix_construct_redim(double, dEdF,DIM_3x3x3x3,1);  
  Matrix_eye(delta,DIM_3);
  
  for(int K=1; K<=DIM_3; K++)
  {
    for(int L=1; L<=DIM_3; L++)
    {

      for(int P=1; P<=DIM_3; P++)
      {
        for(int Q=1; Q<=DIM_3; Q++)
         Tns4_v(dEdF,P,Q,K,L) = Mat_v(delta,P,L)*Mat_v(*F,K,Q) + Mat_v(*F,K,P)*Mat_v(delta,Q,L);
      }
    }
  }
                
  
  for(int I=1; I<=DIM_3; I++)
  {
    for(int J=1; J<=DIM_3; J++)
    {
      for(int K=1; K<=DIM_3; K++)
      {
        for(int L=1; L<=DIM_3; L++)
        {
          for(int A=1; A<=DIM_3; A++)
          {
            for(int B=1; B<=DIM_3; B++)
            {
              Tns6_v(*d2PdF2,I,J,K,L,A,B) = 0.0;
              
              for(int P=1; P<=DIM_3; P++)
              {
                for(int Q=1; Q<=DIM_3; Q++)
                {
                  double d2EdF2_PQKLAB = Mat_v(delta,P,L)*Mat_v(delta,K,A)*Mat_v(delta,Q,B)
                                       + Mat_v(delta,K,A)*Mat_v(delta,P,B)*Mat_v(delta,Q,L);
                                       
                  Tns6_v(*d2PdF2,I,J,K,L,A,B) += Mat_v(delta,I,A)*Tns4_v(*dWdE,B,J,P,Q)*Tns4_v(dEdF,P,Q,K,L);
                  for(int M=1; M<=DIM_3; M++)
                  {
                    Tns6_v(*d2PdF2,I,J,K,L,A,B) += Mat_v(*F,I,M)*Tns4_v(*dWdE,M,J,P,Q)*d2EdF2_PQKLAB;
                    for(int X=1; X<=DIM_3; X++)
                    {
                      Tns6_v(*d2PdF2,I,J,K,L,A,B) += Mat_v(delta,I,K)*Tns4_v(*dWdE,L,J,M,X)*Tns4_v(dEdF,M,X,A,B);
                      for(int Y=1; Y<=DIM_3; Y++)
                        Tns6_v(*d2PdF2,I,J,K,L,A,B) += Mat_v(*F,I,M)*Tns6_v(*dCdE,M,J,P,Q,X,Y)*Tns4_v(dEdF,X,Y,A,B)*Tns4_v(dEdF,P,Q,K,L);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }    

  Matrix_cleanup(delta);
  Matrix_cleanup(dEdF);
  return err;
}

/// compute 2nd derivative of deF/dpF w.r.t hF
///
/// d2eFdhFdpF(I,J,K,L,O,P) = F(I,M)*hFI(M,O)*hFI(P,N)*pFI(N,K)*pFI(L,J)
///
/// \param[out] dF computed 6th order tensor 
/// \param[in] F total deformation tensor
/// \param[in] pFI inverse of the plastic part deformation gradient
/// \param[in] hFI inverse of the thermal part deformation gradient
/// \return non-zero on interal error
int compute_d2eFdhFdpF(Matrix(double) *dF,
                       Matrix(double) *F,
                       Matrix(double) *pFI, 
                       Matrix(double) *hFI)
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
          for(int O=1; O<=DIM_3; O++)
          {
            for(int P=1; P<=DIM_3; P++)
            {
              Tns6_v(*dF,I,J,K,L,O,P) = 0.0;
              for(int M=1; M<=DIM_3; M++)
              {
                for(int N=1; N<=DIM_3; N++)
                  Tns6_v(*dF,I,J,K,L,O,P) += Mat_v(*F,I,M)*Mat_v(*hFI,M,O)*Mat_v(*hFI,P,N)*Mat_v(*pFI,N,K)*Mat_v(*pFI,L,J);
              }
            }
          }
        }
      }
    }
  }
  return err;
}                       
 
/// compute ABC = A:B:C, A,B,C are all 4th order tensors
///
/// ABC(I,J,K,L) = A(I,J,M,N)*B(M,N,O,P)*C(O,P,K,L);
///
/// \param[out] ABC computed 4th order tensor
/// \param[in] A 1st 4th order tensor
/// \param[in] B 2nd 4th order tensor
/// \param[in] C 2nd 4th order tensor
/// \return non-zero on internal error
int compute_Ten4_A_dd_B_dd_C(Matrix(double) *ABC,
                             Matrix(double) *A,
                             Matrix(double) *B,
                             Matrix(double) *C)
{
  int err = 0;
  double temp[DIM_3x3x3x3];
  Matrix(double) T;
  T.m_pdata = temp;
  T.m_row = DIM_3x3x3x3;
  T.m_col = 1;
  
  for(int I=1; I<=DIM_3; I++)
  {
    for(int J=1; J<=DIM_3; J++)
    {
      for(int K=1; K<=DIM_3; K++)
      {
        for(int L=1; L<=DIM_3; L++)
        {
          Tns4_v(T,I,J,K,L) = 0.0;
          for(int M=1; M<=DIM_3; M++)
          {
            for(int N=1; N<=DIM_3; N++)
              Tns4_v(T,I,J,K,L) += Tns4_v(*A,I,J,M,N)*Tns4_v(*B,N,N,K,L);
          }                  
        }
      }
    }
  }
  
  for(int I=1; I<=DIM_3; I++)
  {
    for(int J=1; J<=DIM_3; J++)
    {
      for(int K=1; K<=DIM_3; K++)
      {
        for(int L=1; L<=DIM_3; L++)
        {
          Tns4_v(*ABC,I,J,K,L) = 0.0;
          for(int M=1; M<=DIM_3; M++)
          {
            for(int N=1; N<=DIM_3; N++)
              Tns4_v(*ABC,I,J,K,L) += Tns4_v(T,I,J,M,N)*Tns4_v(*C,N,N,K,L);
          }                  
        }
      }
    }
  }  
  
  return err;
}                             


/// compute deformation gradient due to heat expansion
///
/// \param[out] hF deformation gradient due to heat expansion
/// \param[in] dT temperature difference
/// \param[in] mat MATERIAL_PROPERTY object
/// \param[in] mat_id material id
/// \param[in] diff_order, if 0 deformation gradient
///                        if 1 1st order of differentiation of hF w.r.t temperature
///                        if 2 2nd order of differentiation of hF w.r.t temperature
/// \return non-zero with interal error
int compute_hF(Matrix(double) *hF, 
               double dT,
               const MATERIAL_PROPERTY *mat,
               const int mat_id,
               const int diff_order)
{
  int err = 0.0;
  // compute thermal part of deformation gradient    
  double ax = mat->mater[mat_id].ax;
  double ay = mat->mater[mat_id].ay;
  double az = mat->mater[mat_id].az;  
  
  Matrix_init(*hF, 0.0);
  switch(diff_order)
  {
    case 0:
      Mat_v(*hF, 1,1) = 1.0 + ax*dT;  
      Mat_v(*hF, 2,2) = 1.0 + ay*dT;
      Mat_v(*hF, 3,3) = 1.0 + az*dT;
      break;
    case 1:
      Mat_v(*hF, 1,1) = ax;
      Mat_v(*hF, 2,2) = ay;
      Mat_v(*hF, 3,3) = az;
      break;
    case 2:
      Mat_v(*hF, 1,1) = 0.0;  
      Mat_v(*hF, 2,2) = 0.0;
      Mat_v(*hF, 3,3) = 0.0;
      break;
  }    
  return err;
}                

/// compute differentiation eF w.r.t hF
///
/// \param[out] dF computed 4th order tensor 
/// \param[in] F total deformation tensor
/// \param[in] pFI inverse of the plastic part deformation gradient
/// \param[in] hFI inverse of the thermal part deformation gradient
/// \return non-zero on interal error
int compute_deF_over_dhF(Matrix(double) *dF,
                         Matrix(double) *F,
                         Matrix(double) *pFI, 
                         Matrix(double) *hFI)
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
}                         

/// compute differentiation eF w.r.t pF
///
/// \param[out] dF computed 4th order tensor 
/// \param[in] F total deformation tensor
/// \param[in] pFI inverse of the plastic part deformation gradient
/// \param[in] hFI inverse of the thermal part deformation gradient
/// \return non-zero on interal error
int compute_deF_over_dpF(Matrix(double) *dF,
                         Matrix(double) *F,
                         Matrix(double) *pFI, 
                         Matrix(double) *hFI)
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
}

/// compute heat generation due to mechanical (reference configureation)
///
/// \param[out] Qe thermal source due to mechanical work (elastic part)
/// \param[out] Qp thermal source due to mechanical work (plastic part)
/// \param[out] DQ tangent of heat generation w.r.t temperature (computed only if compute_tangent = 1)
/// \param[in] mat MATERIAL_PROPERTY object
/// \param[in] fv_m mechanical FIELD_VARIABLES object
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
                                const MATERIAL_PROPERTY *mat,
                                const FIELD_VARIABLES *fv_m,
                                const double T,
                                const double dT,
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
  
  double Tdot = dT/dt;
  double hJ = 0.0;
  double pJ = 0.0;
    
  // construct 4th order tensors
  enum {dWdE,dPdF,d2fdhF2,deFdhF,d2eFdhF2,deFdpF,hB,F4_temp,F4end};
  Matrix(double) *F4 = malloc(F4end*sizeof(Matrix(double)));
  for (int a = 0; a < F4end; a++)
    Matrix_construct_init(double, F4[a],DIM_3x3x3x3,1,0.0);
  
  // construct 2nd order tensors    
  enum {F,eF,pF,pFI,pFn,pFdot,hF,hFI,hFp,hFpp,eS,eP,dfdhF,F2_temp,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++)
    Matrix_construct_init(double, F2[a],DIM_3,DIM_3 ,0.0);

  // 1. compute deformation gradient of thermal expansions
  int diff_order = 0; // if 0 deformation gradient
                      // if 1 1st order of differentiation of hF w.r.t temperature
                      // if 2 2nd order of differentiation of hF w.r.t temperature
  err += compute_hF(F2+hF,   dT, mat, mat_id, diff_order);
  err += compute_hF(F2+hFp,  dT, mat, mat_id, 1);
  err += compute_hF(F2+hFpp, dT, mat, mat_id, 2);
  
  Matrix_inv(F2[hF],F2[hFI]);
  Matrix_det(F2[hF],hJ);
  
  // 2. obtain deformation gradient from mechanical part 
  Constitutive_model *m = &(fv_m->eps[eid].model[ip-1]);
  const Model_parameters *func = m->param;
  ELASTICITY *elast = (m->param)->cm_elast;
  double *tempS = elast->S; // temporal pointer to update F4[dWdE], and F2[S] using elast
  double *tempL = elast->L;   
  elast->S = F2[eS].m_pdata;
  elast->L = F4[dWdE].m_pdata;
  
  int stepno = 1; // 0 = time step = n-1
                  // 1 = time step = n
                  // 2 = time step = n+1
  err += m->param->get_eF_of_hF(m,F2+eF,F2+hFI,stepno);
  
  elast->update_elasticity(elast,F2[eF].m_pdata,compute_stiffness);
  Matrix_AxB(F2[eP],1.0,0.0,F2[eF],0,F2[eS],0);
  
  // When compute deformation gradients using m->param->get_xF functions
  // all xF(t(n)) are updated from xF(t(n+1)) such that if xFn is needed
  // temporal field variable should be used for coupled problem.
  err += func->get_F(  m, F2+F);   // this brings  F(t(n+1))       
  err += func->get_pF( m, F2+pF);  //             pF(t(n+1))
  err += func->get_pFn(m, F2+pFn); //             pF(t(n))  
  
  Matrix_det(F2[pF],pJ);

  for(int ia=0; ia<9; ia++)
    F2[pFdot].m_pdata[ia] = (F2[pF].m_pdata[ia] - F2[pFn].m_pdata[ia])/dt;   

  // 4. start computations
  // compute dPdF
  err += compute_dPdF(F4+dPdF,F2+eS,F4+dWdE,F2+eF);
  // compute deFdhF
  Matrix_inv(F2[pF], F2[pFI]);
  err += compute_deF_over_dhF(F4+deFdhF,F2+F,F2+pFI,F2+hFI);
  // compute dfdhF, f=Psi(potential function) 
  Matrix_Tns2_dd_Tns4(F2[dfdhF],F2[eP],F4[deFdhF]);
  // compute d2fdhF2 
  err += compute_Ten4_A_dd_B_dd_C(F4+d2fdhF2, F4+dPdF,F4+deFdhF,F4+deFdhF);

  // compute deFdpF
  err += compute_deF_over_dpF(F4+deFdpF,F2+F,F2+pFI,F2+hFI);
  // compute heat generations
  // xi = dhFdT:d2fdhF2:dhFdT + dfdhF:hFpp
  
  // Heat gen by plastic field variabls is not used (commented out). Instead, 
  // lumped into overall heat gen by plasticity and take 80%(PLASTIC_HEAT_FACTOR)from it
  // 
  // double eff   = 0.0;
  // double g   = 0.0;
  // err += compute_effective_dot_strain(&eff, F2[pF].m_pdata, F2[pFn].m_pdata, dt);
  // err += func->get_hardening(m,&g);
  // Q_p = hJ*g*eff; 
  
  double xi = 0.0;
  double Q_p = 0.0;
  for(int I = 1; I<=DIM_3; I++)
  {
    for(int J = 1; J<=DIM_3; J++)
    {
      xi += Mat_v(F2[dfdhF],I,J)*Mat_v(F2[hFpp],I,J);
      for(int K = 1; K<=DIM_3; K++)
      {
        for(int L = 1; L<=DIM_3; L++)
        {
          xi += Mat_v(F2[hFp],I,J)*Tns4_v(F4[d2fdhF2],I,J,K,L)*Mat_v(F2[hFp],K,L);
          Q_p -= Mat_v(F2[eP],I,J)*Tns4_v(F4[deFdpF],I,J,K,L)*Mat_v(F2[pFdot],K,L);
        }
      }
    }
  } 
  
  *Qe = hJ*pJ*xi*T*Tdot;  
  *Qp = PLASTIC_HEAT_FACTOR*hJ*pJ*Q_p;
    
  // compute tangent of Qe
  double DQe =0.0;
  double DQp = 0.0;
  if(compute_tangent)
  {
    Matrix(double) d3WdC3, d2PdF2, df3dhF3, d2eFdhFdpF;
    Matrix_construct_redim(double,d3WdC3,    DIM_3x3x3x3*DIM_3x3,1);
    Matrix_construct_redim(double,d2PdF2,    DIM_3x3x3x3*DIM_3x3,1);
    Matrix_construct_redim(double,df3dhF3,   DIM_3x3x3x3*DIM_3x3,1);
    Matrix_construct_redim(double,d2eFdhFdpF,DIM_3x3x3x3*DIM_3x3,1);
        
    err += elast->compute_d3W_dC3(elast,F2[eF].m_pdata,d3WdC3.m_pdata);
    err += compute_d2PdF2(&d2PdF2,F2+eS,F4+dWdE,&d3WdC3,F2+eF);
    err += compute_d2eFdhFdpF(&d2eFdhFdpF,F2+F,F2+pFI,F2+hFI);
    
    for(int I=1; I<=DIM_3; I++)
    {
      for(int J=1; J<=DIM_3; J++)
      {
        for(int K=1; K<=DIM_3; K++)
        {
          for(int L=1; L<=DIM_3; L++)
          {
            for(int M=1; M<=DIM_3; M++)
            {
              for(int N=1; N<=DIM_3; N++)
              {
                Tns6_v(df3dhF3,I,J,K,L,M,N) = 0.0;
                DQp -= Mat_v(F2[eP],I,J)*Tns6_v(d2eFdhFdpF,I,J,K,L,M,N)*Mat_v(F2[hFp],M,N)*Mat_v(F2[pFdot],K,L); 
                for(int O=1; O<=DIM_3; O++)
                {
                  for(int P=1; P<=DIM_3; P++)
                  {
                    DQp -= Tns4_v(d2PdF2,I,J,K,L)*Tns4_v(F4[deFdhF],K,L,M,N)*Mat_v(F2[hFp],M,N)*Tns4_v(F4[deFdpF],I,J,O,P)*Mat_v(F2[pFdot],O,P);                    
                    for(int A=1; A<=DIM_3; A++)
                    {
                      for(int B=1; B<=DIM_3; B++)
                      {
                        for(int C=1; C<=DIM_3; C++)
                        {
                          for(int D=1; D<=DIM_3; D++)
                            Tns6_v(df3dhF3,I,J,K,L,M,N) += Tns6_v(d2PdF2,I,J,K,L,O,P)*Tns4_v(F4[deFdhF],O,P,A,B)*Tns4_v(F4[deFdhF],A,B,C,D)*Tns4_v(F4[deFdhF],C,D,M,N);                          
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    double Dxi = 0.0;
    
    for(int I=1; I<=DIM_3; I++)
    {
      for(int J=1; J<=DIM_3; J++)
      {
        Dxi += Mat_v(F2[dfdhF],I,J)*Mat_v(F2[hFpp],I,J);
        for(int K=1; K<=DIM_3; K++)
        {
          for(int L=1; L<=DIM_3; L++)
          {
            Dxi += Mat_v(F2[hFpp],I,J)*Tns4_v(F4[d2fdhF2],I,J,K,L)*Mat_v(F2[hFp], K,L); 
                +  Mat_v(F2[hFp], I,J)*Tns4_v(F4[d2fdhF2],I,J,K,L)*Mat_v(F2[hFpp],K,L);
                +  Tns4_v(F4[d2fdhF2],I,J,K,L)*Mat_v(F2[hFp],K,L)*Mat_v(F2[hFpp], I,J);
                +  Mat_v(F2[dfdhF],I,J)*Mat_v(F2[hFpp],I,J);

            for(int M=1; M<=DIM_3; M++)
            {
              for(int N=1; N<=DIM_3; N++)
                Dxi += Mat_v(F2[hFp],I,J)*Tns6_v(df3dhF3,I,J,K,L,M,N)*Mat_v(F2[hFp],M,N)*Mat_v(F2[hFp],K,L);
            } 
          }
        }
      }
    }
 
    DQe = xi*(dT/dt + T/dt) + Dxi*dT/dt*T;
    Matrix_cleanup(d3WdC3);
    Matrix_cleanup(d2PdF2);
    Matrix_cleanup(df3dhF3);
    Matrix_cleanup(d2eFdhFdpF);
  }
      
  *DQ = hJ*pJ*(DQe + PLASTIC_HEAT_FACTOR*DQp);

  elast->S = tempS;
  elast->L = tempL;

  for(int a = 0; a < F4end; a++)
    Matrix_cleanup(F4[a]);   
  free(F4);
  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);   
  free(F2);

  return err;
}     
                                      
/// compute temperature on node in a element
///
/// \param[in] cn local DOF ids for field variables
/// \param[in] ndofe number of degree of freedom for element
/// \param[in] T nodal temperature for all node
/// \param[in] dT temperature increment
/// \param[in] elem ELEMENT object
/// \param[in] node NODE object
/// \param[in] sup struct for BC's
/// \param[out] T_e computed nodal temperature for current element
/// \param[in] T0 reference temperature
/// \return non-zero on internal error
int get_temperature_elem(const long *cn,
                         const long ndofe,
                         const double *T,
                         const double *dT,
                         const ELEMENT *elem,
                         const NODE *node,
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
                                       GRID *grid,
                                       FIELD_VARIABLES *fv,
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
                                           double *fi_in,
                                           double *du,
                                           GRID *grid,
                                           MATERIAL_PROPERTY *mat,
                                           FIELD_VARIABLES *fv,
                                           LOADING_STEPS *load,
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
  
  MATERIAL_THERMAL *thermal = (mat->thermal) + mat_id;  
  Matrix(double) k0,k;
  k0.m_pdata = thermal->k;
  k0.m_row = k0.m_col = DIM_3;
  
  Matrix_construct_init(double, k,DIM_3,DIM_3,0.0);
  Matrix_AeqB(k,1.0,k0);
  
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  
  double cp = thermal->cp;
  double Q = 0.0;
  
  double dt = 1.0; // for the quasi steady state
  if(rho_0>0)
    dt = dt_in;  
  
  SUPP sup = load->sups[mp_id];

  long *nod = fe->node_id.m_pdata; // use only address, no need deallication
  int ndofn = fv->ndofn;  
  long ndofe = (fe->nne)*ndofn;
  long *cnL = aloc1l(ndofe);
  
  Matrix(double) fi;
  fi.m_pdata = fi_in;
  fi.m_row = ndofe;  
  
  get_dof_ids_on_elem_nodes(0,fe->nne,ndofn,nod,grid->node,cnL,mp_id); 
  
  Matrix(double) q, Tnp1, Tn;  
  Matrix_construct_redim(double, q,grid->nsd,1);
  Matrix_construct_init(double, Tnp1,fe->nne,1,0.0); 
  Matrix_construct_init(double, Tn,  fe->nne,1,0.0);
  
  //compute nodal value
  get_temperature_elem(cnL,ndofe,fv->u_np1,du,grid->element,grid->node,sup,Tnp1.m_pdata,fv->u0);
  
  for(int ia=0; ia<fe->nne; ia++)
    Vec_v(Tn,   ia+1) = fv->u_n[nod[ia]];    

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    // Udate basis functions at the integration points.
    FEMLIB_elem_basis_V(fe, ip);
    FEMLIB_update_shape_tensor(fe);
        
    double Temp = 0.0;
    double dT   = 0.0;
    Matrix_init( q,0.0);
            
    // compute varialbes at the integration point
    for(int ia=1; ia<=fe->nne; ia++)
    { 
      Temp += Vec_v(fe->N,ia)*Vec_v(Tnp1, ia);
      dT   += Vec_v(fe->N,ia)*(Vec_v(Tnp1, ia)-Vec_v(Tn, ia));
    }
    
    if(is_it_couple_w_mechanical>=0)
    {  
      FIELD_VARIABLES *fv_m = fv->fvs[is_it_couple_w_mechanical];
      int compute_tangent = 0;
      double DQ = 0.0;
      double Qe = 0.0;
      double Qp = 0.0;

      Matrix(double) F;      
      Matrix_construct_init(double, F, DIM_3,DIM_3,0.0);
                  
      Constitutive_model *m = &(fv_m->eps[eid].model[ip-1]);
      const Model_parameters *func = m->param;
      err += func->get_F(m, &F);        // this brings  F(t(n+1))
      
      if(total_Lagrangian)
      {
        double detF = 1.0;
        Matrix_det(F, detF);
        
        Matrix(double) FI, FIT;
        Matrix_construct_init(double, FI,  DIM_3,DIM_3,0.0);
        Matrix_construct_init(double, FIT, DIM_3,DIM_3,0.0);
          
        Matrix_inv(F,FI);
        Matrix_AeqBT(FIT,1.0,FI);
        Matrix_Tns2_AxBxC(k,detF,0.0,FI,k0,FIT); // k = J*FI*k0*FIT
        
        Matrix_cleanup(FI);
        Matrix_cleanup(FIT);        
      }
      else
        Matrix_det(F, Jn); // updated Lagrangian
          
      Matrix_cleanup(F);
      
      if(thermal->FHS_MW>TOL_FHS)
      {   
        err += compute_mechanical_heat_gen(&Qe,&Qp,&DQ,mat,fv_m,Temp,dT,dt,eid,ip,mat_id,compute_tangent);
        Q += thermal->FHS_MW*(Qe + Qp);
        
      }
    }
    Matrix_init( q,0.0);
            
    // compute heat flux
    for(int ia=1; ia<=fe->nne; ia++)
    { 
      // k = [nsd, nsd], dN = [nne, nsd], Tnp1  = [nne, 1],
      // q = [nsd, 1] = k*dN'*T     
      for(int ib=1; ib<=grid->nsd; ib++)
      {
        for(int ic=1; ic<=grid->nsd; ic++)        
          Vec_v(q,ib) += Mat_v(k, ib, ic)*Mat_v(fe->dN,ia,ic)*Vec_v(Tnp1, ia);
      }        
    }    

    // R = rho*cp*dT + dt*grad.q - dt*Q = 0;
    // rho = rho_0/Jn
    // Q = Q_0/Jn;
    for(int ia=1; ia<=fe->nne; ia++)
    {      
      Vec_v(fi,ia) += Vec_v(fe->N,ia)*(rho_0*cp*dT - dt*Q)*(fe->detJxW)/Jn;
      for(int ib=1; ib<=grid->nsd; ib++)
        Vec_v(fi,ia) += dt*Mat_v(fe->dN,ia,ib)*Vec_v(q,ib)*(fe->detJxW);
    }  
  }
  
  free(cnL);
  Matrix_cleanup(k);
  Matrix_cleanup(q);
  Matrix_cleanup(Tnp1);
  Matrix_cleanup(Tn);
  
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
int energy_equation_compute_residuals(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      LOADING_STEPS *load,
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
  
  // save original pointer to access mechanical part
  double *u_n;
  double *u_nm1;
  State_variables *statv_list = NULL;

  if(is_it_couple_w_mechanical>=0)
  {
    FIELD_VARIABLES *fv_m = fv->fvs[is_it_couple_w_mechanical];    
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
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);
    
    // do volume integration at an element, 
    // fe needs to be updated by integration points 
    Matrix(double) fi;
    Matrix_construct_init(double, fi, fe.nne, 1, 0.0);
    err += energy_equation_compute_residuals_elem(&fe,fi.m_pdata,du,grid,mat,fv,load,mp_id,dt,total_Lagrangian);
    err += energy_equation_residuals_assemble(&fe,fi.m_pdata,grid,fv,mp_id);
    
    Matrix_cleanup(fi);
    FEMLIB_destruct(&fe);
  }
  
  if(is_it_couple_w_mechanical>=0)
  {
    FIELD_VARIABLES *fv_m = fv->fvs[is_it_couple_w_mechanical];    
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
                                           GRID *grid,
                                           MATERIAL_PROPERTY *mat,
                                           FIELD_VARIABLES *fv,
                                           SOLVER_OPTIONS *sol,
                                           LOADING_STEPS *load,
                                           COMMUNICATION_STRUCTURE *com,
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
  
  MATERIAL_THERMAL *thermal = (mat->thermal) + mat_id;  
  Matrix(double) k0,k;
  k0.m_pdata = thermal->k;
  k0.m_row = k0.m_col = DIM_3;
  
  Matrix_construct_init(double, k,DIM_3,DIM_3,0.0);
  Matrix_AeqB(k,1.0,k0);
  
  double cp = thermal->cp;
  
  double dt = 1.0; // for the quasi steady state
  if(rho_0>0)
    dt = dt_in;
  
  SUPP sup = load->sups[mp_id];
  
  long *nod = fe->node_id.m_pdata;
  int ndofn = fv->ndofn;  
  long ndofe = (fe->nne)*ndofn;
  long *cnL = aloc1l(ndofe);
  long *cnG = aloc1l(ndofe);
  
  get_dof_ids_on_elem_nodes(0,fe->nne,ndofn,nod,grid->node,cnL,mp_id);
  get_dof_ids_on_elem_nodes(1,fe->nne,ndofn,nod,grid->node,cnG,mp_id);
  
  Matrix(double) lk;
  Matrix_construct_init(double,lk,ndofe,ndofe,0.0);
  
  Matrix(double) Tnp1, Tn;  
  Matrix_construct_init(double, Tnp1,fe->nne,1,0.0); 
  Matrix_construct_init(double, Tn,  fe->nne,1,0.0);
    
  //compute nodal value
  get_temperature_elem(cnL,ndofe,fv->u_np1,fv->f,grid->element,grid->node,sup,Tnp1.m_pdata,fv->u0);

  for(int ia=0; ia<fe->nne; ia++)
    Vec_v(Tn,   ia+1) = fv->u_n[nod[ia]];    

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    // Udate basis functions at the integration points.    
    FEMLIB_elem_basis_V(fe, ip);
    FEMLIB_update_shape_tensor(fe);
    
    double Temp = 0.0;
    double dT   = 0.0;
    
    // compute varialbes at the integration point
    for(int ia=1; ia<=fe->nne; ia++)
    { 
      Temp += Vec_v(fe->N,ia)*Vec_v(Tnp1, ia);
      dT   += Vec_v(fe->N,ia)*(Vec_v(Tnp1, ia)-Vec_v(Tn, ia));
    }
    
    double DQ = 0.0;
    
    if(is_it_couple_w_mechanical>=0)
    {  
      FIELD_VARIABLES *fv_m = fv->fvs[is_it_couple_w_mechanical];
      int compute_tangent = 1;
      double Qe = 0.0;
      double Qp = 0.0;

      Matrix(double) F;
      Matrix_construct_init(double, F, DIM_3,DIM_3,0.0);
      Constitutive_model *m = &(fv_m->eps[eid].model[ip-1]);
      const Model_parameters *func = m->param;
      err += func->get_F(m, &F);        // this brings  F(t(n+1))

      if(total_Lagrangian)
      {
        double detF = 1.0;
        Matrix_det(F, detF);
        
        Matrix(double) FI, FIT;
        Matrix_construct_init(double, FI,  DIM_3,DIM_3,0.0);
        Matrix_construct_init(double, FIT, DIM_3,DIM_3,0.0);
          
        Matrix_inv(F,FI);
        Matrix_AeqBT(FIT,1.0,FI);
        Matrix_Tns2_AxBxC(k,detF,0.0,FI,k0,FIT); // k = J*FI*k0*FIT
        
        Matrix_cleanup(FI);
        Matrix_cleanup(FIT);        
      }
      else
        Matrix_det(F, Jn); // updated Lagrangian
        
      Matrix_cleanup(F);
            
      if(thermal->FHS_MW>TOL_FHS)
      {   
        err += compute_mechanical_heat_gen(&Qe,&Qp,&DQ,mat,fv_m,Temp,dT,dt,eid,ip,mat_id,compute_tangent);
        DQ = thermal->FHS_MW*DQ;
      }
    }
    
    // R = rho_0*cp*(Tnp1-Tn) + dt*grad.q - dt*Q = 0;
    // DR = rho_0*cp*dT + dt*D[grad.q]dT - dt*D[Q]dT = 0;    
    for(int ia=1; ia<=fe->nne; ia++)
    {
      for(int ib=1; ib<=fe->nne; ib++)
      {
        Mat_v(lk,ia,ib) += (rho_0*cp - dt*DQ)*Vec_v(fe->N,ia)*Vec_v(fe->N,ib)*(fe->detJxW)/Jn;
        for(int im = 1; im<=grid->nsd; im++)
        {
          for(int in = 1; in<=grid->nsd; in++)
            Mat_v(lk,ia,ib) += dt*Mat_v(fe->dN,ia,in)*Mat_v(k, in, im)*Mat_v(fe->dN,ib,im)*(fe->detJxW);
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
                cnL,cnG,ndofe,Ddof,
                com->GDof,
                myrank,
                com->nproc,
                com->comm,
                interior,
                sol->PGFEM_hypre,
                opts->analysis_type);                               
  }
  
  if(compute_load4pBCs)
  {
    int ndofn = 1;
    int k = 0;
    int jj = 0;
    Matrix(double) u, f_loc;    
    Matrix_construct_init(double,u    ,(fe->nne)*ndofn,1,0.0);
    Matrix_construct_init(double,f_loc,(fe->nne)*ndofn,1,0.0);
    
    // get the bc increment
    for(int ia=0; ia<fe->nne; ia++)
    {
      for(int ib=0; ib<ndofn; ib++)
      {
        int id = ia*ndofn + ib;
        if(cnL[id] <= -1)
          u.m_pdata[id] = load->sups[mp_id]->defl_d[abs(cnL[id])-1];
        else
          u.m_pdata[id] = 0.0;
      }
    }
    Matrix_AxB(f_loc,1.0,0.0,lk,0,u,0);
    
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
    
    Matrix_cleanup(u);
    Matrix_cleanup(f_loc);        
  }
  
  Matrix_cleanup(lk);
    
  free(cnL);
  free(cnG);
  Matrix_cleanup(k);
  Matrix_cleanup(Tnp1);
  Matrix_cleanup(Tn);  
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
int energy_equation_compute_stiffness(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      SOLVER_OPTIONS *sol,
                                      LOADING_STEPS *load,
                                      COMMUNICATION_STRUCTURE *com,
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
  
  // save original pointer to access mechanical part
  double *u_n;
  double *u_nm1;
  State_variables *statv_list = NULL;

  if(is_it_couple_w_mechanical>=0)
  {
    FIELD_VARIABLES *fv_m = fv->fvs[is_it_couple_w_mechanical];    
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

  Matrix(int) Ddof;
  Matrix_construct_redim(int, Ddof,com->nproc,1);
  
  Ddof.m_pdata[0] = com->DomDof[0];
  for (int ia=1; ia<com->nproc; ia++)
    Ddof.m_pdata[ia] = Ddof.m_pdata[ia-1] + com->DomDof[ia];
  
  for(int eid=0; eid<com->nbndel; eid++)
  {
    // construct finite element library
    // it provide element wise integration info 
    // such as basis function, weights, ... 
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe,com->bndel[eid],grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    int interior = 0;
    err += energy_equation_compute_stiffness_elem(&fe,Lk,grid,mat,fv,sol,load,com,Ddof.m_pdata,
                                                  myrank,interior,opts,mp_id,
                                                  do_assemble,compute_load4pBCs,dt,total_Lagrangian);
    
    FEMLIB_destruct(&fe);
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
    
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);
    
    // do volume integration at an element
    int interior = 1;
    
    err += energy_equation_compute_stiffness_elem(&fe,Lk,grid,mat,fv,sol,load,com,Ddof.m_pdata,
                                                  myrank,interior,opts,mp_id,
                                                  do_assemble,compute_load4pBCs,dt,total_Lagrangian);
      
    FEMLIB_destruct(&fe);
    if(err != 0)
      break;
  }        

  err += assemble_nonlocal_stiffmat(com->comm,sta_r,req_r,sol->PGFEM_hypre,recieve);
  err += finalize_stiffmat_comm(sta_s,sta_r,req_s,req_r,com->comm);

  if(is_it_couple_w_mechanical>=0)
  {
    FIELD_VARIABLES *fv_m = fv->fvs[is_it_couple_w_mechanical];    
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
  
  Matrix_cleanup(Ddof);
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
int energy_equation_compute_load4pBCs(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      SOLVER_OPTIONS *sol,
                                      LOADING_STEPS *load,
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
  
  // save original pointer to access mechanical part
  double *u_n;
  double *u_nm1;
  State_variables *statv_list = NULL;

  if(is_it_couple_w_mechanical>=0)
  {
    FIELD_VARIABLES *fv_m = fv->fvs[is_it_couple_w_mechanical];    
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
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);

    // do volume integration at an element
    err += energy_equation_compute_stiffness_elem(&fe,NULL,grid,mat,fv,sol,load,NULL,NULL,
                                                  myrank,interior,opts,mp_id,
                                                  do_assemble,compute_load4pBCs,dt,total_Lagrangian);    
    FEMLIB_destruct(&fe);
    if(err != 0)
      break;      
  }
  
  if(is_it_couple_w_mechanical>=0)
  {
    FIELD_VARIABLES *fv_m = fv->fvs[is_it_couple_w_mechanical];    
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
int update_thermal_flux4print(GRID *grid,
                              MATERIAL_PROPERTY *mat,
                              FIELD_VARIABLES *fv,
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
    
  // save original pointer to access mechanical part
  double *u_n;
  double *u_nm1;
  State_variables *statv_list = NULL;

  if(is_it_couple_w_mechanical>=0)
  {
    FIELD_VARIABLES *fv_m = fv->fvs[is_it_couple_w_mechanical];    
    u_n   = fv_m->u_n;
    u_nm1 = fv_m->u_nm1;     
    statv_list = fv_m->statv_list;
  
    fv_m->u_n   = fv_m->temporal->u_n;
    fv_m->u_nm1 = fv_m->temporal->u_nm1;
    fv_m->statv_list = fv_m->temporal->var;
  }  
  
  int myrank = 0;
  MPI_Comm_rank (MPI_COMM_WORLD,&myrank);
  
  for(int eid=0; eid<grid->ne; eid++)
  {    
    FEMLIB fe;
    FEMLIB_initialization_by_elem(&fe,eid,grid->element,grid->node,intg_order,total_Lagrangian);

    memset(eps[eid].el.o,0,6*sizeof(double));

    // get material constants (parameters)
    const int mat_id = (grid->element[eid]).mat[0];  
    double rho_0 = mat->density[mat_id];
    
    MATERIAL_THERMAL *thermal = (mat->thermal) + mat_id;  
    Matrix(double) k0,k;
    k0.m_pdata = thermal->k;
    k0.m_row = k0.m_col = DIM_3;
  
    Matrix_construct_init(double, k,DIM_3,DIM_3,0.0);
    Matrix_AeqB(k,1.0,k0);

    double cp = thermal->cp;
      
    // compute noda values
    Matrix(double) q, Tnp1, Tn;  
    Matrix_construct_redim(double, q,grid->nsd,1);    
    Matrix_construct_init(double, Tnp1,fe.nne,1,0.0); 
    Matrix_construct_init(double, Tn,  fe.nne,1,0.0);
    
    long *nod = fe.node_id.m_pdata;
    for(int ia=0; ia<fe.nne; ia++)
    {
      Vec_v(Tnp1, ia+1) = fv->u_n[nod[ia]];
      Vec_v(Tn,   ia+1) = fv->u_nm1[nod[ia]];
    }
    
    // do integration point loop
    double volume = 0.0;
    for(int ip = 1; ip<=fe.nint; ip++)
    {
      FEMLIB_elem_basis_V(&fe, ip);
      FEMLIB_update_shape_tensor(&fe);
    
      double Temp = 0.0;
      double dT   = 0.0;
      
      // compute values at the integration point 
      for(int ia=1; ia<=fe.nne; ia++)
      {         
        Temp += Vec_v(fe.N,ia)*Vec_v(Tnp1, ia);
        dT   += Vec_v(fe.N,ia)*(Vec_v(Tnp1, ia)-Vec_v(Tn, ia));
      }
      
      // compute heat sources 
      double Qe = 0.0;
      double Qp = 0.0;
      double DQ = 0.0;
      double Jn = 1.0;      
      if(is_it_couple_w_mechanical>=0)
      {
        FIELD_VARIABLES *fv_m = fv->fvs[is_it_couple_w_mechanical];
        int compute_tangent = 0;

        Matrix(double) F;
        Matrix_construct_init(double, F, DIM_3,DIM_3,0.0);
        Constitutive_model *m = &(fv_m->eps[eid].model[ip-1]);
        const Model_parameters *func = m->param;
        err += func->get_F(m, &F);        // this brings  F(t(n+1))
  
        if(total_Lagrangian)
        {
          double detF = 1.0;
          Matrix_det(F, detF);
          
          Matrix(double) FI, FIT;
          Matrix_construct_init(double, FI,  DIM_3,DIM_3,0.0);
          Matrix_construct_init(double, FIT, DIM_3,DIM_3,0.0);
            
          Matrix_inv(F,FI);
          Matrix_AeqBT(FIT,1.0,FI);
          Matrix_Tns2_AxBxC(k,detF,0.0,FI,k0,FIT); // k = J*FI*k0*FIT
          
          Matrix_cleanup(FI);
          Matrix_cleanup(FIT);        
        }
        else
          Matrix_det(F, Jn); // updated Lagrangian

        Matrix_cleanup(F);
        
        if(thermal->FHS_MW>TOL_FHS)
        { 
          err += compute_mechanical_heat_gen(&Qe,&Qp,&DQ,mat,fv_m,Temp,dT,dt,eid,ip,mat_id,compute_tangent);
          Qe = thermal->FHS_MW*Qe;
          Qp = thermal->FHS_MW*Qp;
        }
      }
      
      Matrix_init( q,0.0);
              
      // compute heat flux
      for(int ia=1; ia<=fe.nne; ia++)
      { 
        // k = [nsd, nsd], dN = [nne, nsd], Tnp1  = [nne, 1],
        // q = [nsd, 1] = k*dN'*T     
        for(int ib=1; ib<=grid->nsd; ib++)
        {
          for(int ic=1; ic<=grid->nsd; ic++)        
            Vec_v(q,ib) -= Mat_v(k, ib, ic)*Mat_v(fe.dN,ia,ic)*Vec_v(Tnp1, ia);
        }        
      }       
      
      // save the values
      fv->eps[eid].el.o[0] += fe.detJxW*Vec_v(q,1);
      fv->eps[eid].el.o[1] += fe.detJxW*Vec_v(q,2);
      fv->eps[eid].el.o[2] += fe.detJxW*Vec_v(q,3);
      fv->eps[eid].el.o[3] += fe.detJxW*Qe/Jn;
      fv->eps[eid].el.o[4] += fe.detJxW*Qp/Jn;
      volume += fe.detJxW;
    }
    for(int ia=0; ia<6; ia++)
      eps[eid].el.o[ia] = eps[eid].el.o[ia]/volume;
    
    Matrix_cleanup(k);
    Matrix_cleanup(q);
    Matrix_cleanup(Tnp1);
    Matrix_cleanup(Tn);
  
    FEMLIB_destruct(&fe);    
  }
  
  if(is_it_couple_w_mechanical>=0)
  {
    FIELD_VARIABLES *fv_m = fv->fvs[is_it_couple_w_mechanical];    
    fv_m->u_n   = u_n;
    fv_m->u_nm1 = u_nm1;
    fv_m->statv_list = statv_list;
  }
   
  return err;
}
