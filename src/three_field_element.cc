#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "three_field_element.h"
#include "Hu_Washizu_element.h"
#include "allocation.h"
#include "condense.h"
#include "cast_macros.h"
#include "def_grad.h"
#include "displacement_based_element.h"
#include "dynamics.h"
#include "enumerations.h"
#include "femlib.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "index_macros.h"
#include "utils.h"
#include <mkl_cblas.h>

static constexpr int ndn = 3;
static constexpr int USE_HW_FUNCS = 0;
static constexpr int INTG_ORDER  = 1;
static constexpr int DIM_3       = 3;
static constexpr int DIM_3x3     = 9;
static constexpr int DIM_3x3x3x3 = 81;

#include <ttl/ttl.h>

//ttl declarations
namespace {
  template<int R, int D = 3, class S = double>
  using Tensor = ttl::Tensor<R, D, S>;

  template<int R, int D = 3, class S = double *>
  using TensorA = ttl::Tensor<R, D, S>;
    
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'o'> o;
  static constexpr ttl::Index<'p'> p; 
    
  template<class T1, class T2> int inv(T1 &A, T2 &AI)
  {
    int err = inv3x3(A.data, AI.data);
    return err;
  }
  template<class T1, class T2> int symm(T1 &A, T2 &Asym)
  {
    Asym(i,j) = 0.5*(A(i,j)+A(j,i));
    return 0;
  }
  template<class T1> double tr(T1 &A)
  {
    double v = A(i,i);
    return v;
  }        
}

/// Variables at integration point for Total Lagrangian 
/// and displacement based three field mixed
class Var_Carrier
{
  public:
    dUdJFuncPtr      UP;
    d2UdJ2FuncPtr    UPP;
    devStressFuncPtr Stress;
    matStiffFuncPtr  Stiffness;
    double kappa;
    double factor;
    HOMMAT *e_hommat;
    double mF[DIM_3x3], mFI[DIM_3x3], mS[DIM_3x3], mL[DIM_3x3x3x3];
    double J;
    double theta;
    double P;
    double Upp;
    double Up;

    Var_Carrier()
    {
      UP        = NULL;
      UPP       = NULL;
      Stress    = NULL;
      Stiffness = NULL;
    };

    ~Var_Carrier()
    {
      UP        = NULL;
      UPP       = NULL;
      Stress    = NULL;
      Stiffness = NULL;
    };
    void set_elasticity_functions(MaterialProperty *mat,
                                  Grid *grid, 
                                  int eid)
    {
      const int matid = grid->element[eid].mat[2];
      e_hommat = mat->hommat + matid;
      kappa = e_hommat->E/(3.0*(1.0-2.0*(e_hommat->nu)));

      UP = getDUdJFunc(1, e_hommat);
      UPP = getD2UdJ2Func(1, e_hommat);
      Stress = getDevStressFunc(1, e_hommat);
      Stiffness = getMatStiffFunc(1, e_hommat);

    }
    void set_variables(double theta_in, double P_in, bool compute_elasticity_tenosr = true)
    {
      TensorA<2> F(mF), FI(mFI);
      inv(F,FI);
      J = ttl::det(F);
      theta  = theta_in;
      P      = P_in;
      factor = pow(theta_in/J, 2.0/3.0);
      
      Tensor<2> C = F(k,i)*F(k,j);
      
      Upp = compute_d2UdJ2(theta_in);
      Up  = compute_dUdJ(theta_in);
      compute_PK2(C.data, mS);
      
      if(compute_elasticity_tenosr)
        compute_L(C.data, mL);            
    }
    
    double compute_dUdJ(double J)
    {
      double Up = 0.0;
      UP(J, e_hommat, &Up);
      return  kappa*Up;
    }

    double compute_d2UdJ2(double J)
    {
      double Up = 0.0;
      UPP(J, e_hommat, &Up);
      return kappa*Up;
    }
    void compute_PK2(double *C, double *S)
    {
      Stress(C,e_hommat,S);
    }
    void compute_L(double *C, double *L)
    {
      Stiffness(C,e_hommat,L);
    }

    template<class T1, class T2> void compute_A_of_du(T1 &A, T2 &du)
    {
      TensorA<2> F(mF), FI(mFI);
      A = {};
      
      Tensor<2> duTF = du(k,i)*F(k,j);
      Tensor<2> sduTF;
      
      symm(duTF, A);
    }
};

int compute_d_theta_dP_test(Matrix<double> &d_theta,
                            Matrix<double> &dP, 
                            Matrix<double> &du,
                            Matrix<double> &Rt, 
                            Matrix<double> &Rp, 
                            Matrix<double> &Kpu, 
                            Matrix<double> &Ktu, 
                            Matrix<double> &Ktp, 
                            Matrix<double> &Ktt, 
                            Matrix<double> &Kpt);

/// compute stiffness: Kuu at ip
///
/// \param[out] Kuu                    computed Kuu part
/// \param[in]  fe                     fem library object
/// \param[in]  var                    3f related variable object
/// \param[in]  dt_alpha_1_minus_alpha coefficient for inertia or qusi-steady state
void TF_Kuu_ip(Matrix<double> &Kuu,
               FEMLIB *fe,
               Var_Carrier &var,
               double dt_alpha_1_minus_alpha)
{
  int nne = fe->nne;
  int nsd = fe->nsd;
  TensorA<2> S(var.mS), F(var.mF), FI(var.mFI);
  TensorA<4> L(var.mL);
  
  for(int ia=0; ia<nne; ia++)
  {
    for(int ib=0; ib<nsd; ib++)
    {
      int id_ab = ia*nsd + ib + 1;
      TensorA<2> Grad_du(fe->ST + idx_4_gen(ia,ib,0,0,nne,nsd,nsd,nsd));
      Tensor<2> A_ab;
      var.compute_A_of_du(A_ab, Grad_du);
      double FIT_Grad_du = FI(j,i)*Grad_du(i,j);
      
      for(int iw=0; iw<nne; iw++)
      {
        for(int ig=0; ig<nsd; ig++)
        { 
          int id_wg = iw*nsd + ig + 1;
          TensorA<2> Grad_tu(fe->ST + idx_4_gen(iw,ig,0,0,nne,nsd,nsd,nsd));

          // eq. 1
          Tensor<2> A_wg;
          var.compute_A_of_du(A_wg, Grad_tu);         
          double ALA = A_ab(i,j)*L(i,j,k,l)*A_wg(k,l);
          
          // eq. 2
          Tensor<2> duTtu = Grad_du(k,i)*Grad_tu(k,j);
          Tensor<2> sduTtu;
          symm(duTtu, sduTtu);
          double sduTtuS = S(i,j)*sduTtu(i,j);

          // eq. 3          
          Tensor<2> FIT_Grad_tu_FIT = FI(k,i)*Grad_tu(l,k)*FI(j,l);
          double FIT_Grad_tu_FIT_du = FIT_Grad_tu_FIT(i,j)*Grad_du(i,j);
                                          
          // eq. 4
          double FIT_Grad_tu = FI(j,i)*Grad_tu(i,j);

          Kuu(id_ab,id_wg) += -dt_alpha_1_minus_alpha*fe->detJxW*(ALA + sduTtuS
                                     + var.P*var.J*(FIT_Grad_du*FIT_Grad_tu - FIT_Grad_tu_FIT_du));
        }
      }
    }
  }
  
}

/// compute stiffness: Kup at ip
///
/// \param[out] Kup                    computed Kup part
/// \param[in]  fe                     fem library object
/// \param[in]  var                    3f related variable object
/// \param[in]  Pno                    number of pressure variables
/// \param[in]  Np                     shape function for pressure
/// \param[in]  dt_alpha_1_minus_alpha coefficient for inertia or qusi-steady state
void TF_Kup_ip(Matrix<double> &Kup,
               FEMLIB *fe,
               Var_Carrier &var,
               const int Pno,
               Matrix<double> &Np,
               double dt_alpha_1_minus_alpha)
{
  TensorA<2>FI(var.mFI);
  
  for(int ia=0; ia<fe->nne; ia++)
  {
    for(int ib=0; ib<fe->nsd; ib++)
    {
      const int id_ab = idx_4_gen(ia,ib,0,0,fe->nne,fe->nsd,fe->nsd,fe->nsd);
      TensorA<2> Grad_du((fe->ST)+id_ab);

      for(int iw=0; iw<Pno; iw++)
      {
        int idx_up = idx_K_gen(ia,ib,iw,0,fe->nne,fe->nsd,Pno,1);
        Kup.m_pdata[idx_up] += -dt_alpha_1_minus_alpha*fe->detJxW*var.J*FI(j,i)*Grad_du(i,j)*Np(iw+1);
      }
    }
  }
}

/// compute stiffness: Kpt at ip
///
/// \param[out] Kpt                    computed Kpt part
/// \param[in]  fe                     fem library object
/// \param[in]  var                    3f related variable object
/// \param[in]  Pno                    number of pressure variable
/// \param[in]  Np                     shape function for pressure
/// \param[in]  Vno                    number of volume variables
/// \param[in]  Nt                     shape function for volume
/// \param[in]  dt_alpha_1_minus_alpha coefficient for inertia or qusi-steady state
void TF_Kpt_ip(Matrix<double> &Kpt,
               FEMLIB *fe,
               Var_Carrier &var,
               const int Pno,
               Matrix<double> &Np,               
               const int Vno,
               Matrix<double> &Nt,
               double dt_alpha_1_minus_alpha)
{
  for(int ia=1; ia<=Pno; ia++)
  {
    for(int ib=1; ib<=Vno; ib++)
      Kpt(ia,ib) += dt_alpha_1_minus_alpha*fe->detJxW*Np(ia)*Nt(ib);
  }
}

/// compute stiffness: Kpt at ip
///
/// \param[out] Kpt                    computed Kpt part
/// \param[in]  fe                     fem library object
/// \param[in]  var                    3f related variable object
/// \param[in]  Vno                    number of volume variables
/// \param[in]  Nt                     shape function for volume
/// \param[in]  dt_alpha_1_minus_alpha coefficient for inertia or qusi-steady state
void TF_Ktt_ip(Matrix<double> &Ktt,
               FEMLIB *fe,
               Var_Carrier &var,               
               int Vno,
               Matrix<double> &Nt,
               double dt_alpha_1_minus_alpha)
{
  for(int ia=1; ia<=Vno; ia++)
  {
    for(int iw=1; iw<=Vno; iw++)
    { 
      Ktt(ia,iw) += -dt_alpha_1_minus_alpha*fe->detJxW*Nt(ia)*Nt(iw)*var.Upp;
    }
  }  
}

/// compute residual: Ru at ip
///
/// \param[out] Ru  computed Rt part
/// \param[in]  fe  fem library object
/// \param[in]  var 3f related variable object
void TF_Ru_ip(Matrix<double> &Ru,
              FEMLIB *fe,
              Var_Carrier &var)
{
  TensorA<2> F(var.mF), FI(var.mFI), S(var.mS);
  Tensor<2> C, CI;
  
  C  = F(k,i)*F(k,j);
  CI = FI(i,k)*FI(j,k);
 
  for(int ia=0; ia<fe->nne; ia++)
  {
    for(int ib=0; ib<fe->nsd; ib++)
    {
      const int id_ab = idx_4_gen(ia,ib,0,0,fe->nne,fe->nsd,fe->nsd,fe->nsd);
      TensorA<2> Grad_du(fe->ST + id_ab);
      Tensor<2> AA, sAA;
      AA = F(k,i)*Grad_du(k,j);
      symm(AA,sAA);
      
      double sAAS = sAA(i,j)*S(i,j);
      double C_IFdu = CI(i,j)*AA(i,j);

      Ru.m_pdata[ia*fe->nsd + ib] += fe->detJxW*(sAAS + var.J*var.P*C_IFdu);
    }
  }    
}

/// compute residual: Rt at ip
///
/// \param[out] Rt  computed Rt part
/// \param[in]  fe  fem library object
/// \param[in]  var 3f related variable object
/// \param[in]  Vno number of volume variables
/// \param[in]  Nt  shape function for volume
void TF_Rt_ip(Matrix<double> &Rt,
              FEMLIB *fe,
              Var_Carrier &var,
              int Vno,
              Matrix<double> &Nt)
{
  for(int ia=1; ia<=Vno; ia++)
    Rt(ia) += fe->detJxW*Nt(ia)*(var.Up - var.P);
}

/// compute residual: Rp at ip
///
/// \param[out] Rp  computed Rp part
/// \param[in]  fe  fem library object
/// \param[in]  var 3f related variable object
/// \param[in]  Pno number of pressure variable
/// \param[in]  Np  shape function for pressure
void TF_Rp_ip(Matrix<double> &Rp,
              FEMLIB *fe,
              Var_Carrier &var,
              int Pno,
              Matrix<double> &Np)
{
  for(int ia=1; ia<=Pno; ia++)
    Rp(ia) += fe->detJxW*Np(ia)*(var.J - var.theta);
}              

/// compute element stiffness matrix in quasi steady state
///
/// Total Lagrangian based three-field mixed method for Hyperelasticity.
///
/// \param[in]  fe    finite element helper object
/// \param[out] lk    computed element stiffness matrix
/// \param[in]  r_e   nodal variabls(displacements) on the current element
/// \param[in]  grid  a mesh object
/// \param[in]  mat   a material object
/// \param[in]  fv    object for field variables
/// \param[in]  alpha mid point alpha
/// \param[in]  dt    time step size
/// \return non-zero on internal error
int stiffmat_3f_el(FEMLIB *fe,
                   double *lk,
                   double *r_e,
                   Grid *grid,
                   MaterialProperty *mat,
                   FieldVariables *fv,
                   double alpha,
                   double dt)
{
  int err = 0;

  double alpha_1;
  double alpha_2;
  double dt_alpha_1_minus_alpha;
  
  if(alpha<0)
  {
    alpha = 1.0;
    alpha_1 = 0.0;
    alpha_2 = alpha;
    dt_alpha_1_minus_alpha = -1.0;
  }
  else
  {
    alpha_1 = 1.0 - alpha;
    alpha_2 = alpha;
    dt_alpha_1_minus_alpha = dt*alpha_1*alpha_2;
  }  

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  int ndofe = nne*ndofn;
  int Pno   = fv->npres;
  int Vno   = fv->nVol;
  
  Matrix<double> Kuu(nne*nsd,nne*nsd,0.0);
  Matrix<double> Kut(nne*nsd,Vno    ,0.0);
  Matrix<double> Kup(nne*nsd,Pno    ,0.0);
  Matrix<double> Ktu(Vno    ,nne*nsd,0.0);
  Matrix<double> Ktt(Vno    ,Vno    ,0.0);
  Matrix<double> Ktp(Vno    ,Pno    ,0.0);
  Matrix<double> Kpu(Pno    ,nne*nsd,0.0);
  Matrix<double> Kpt(Pno    ,Vno    ,0.0);
  Matrix<double> Kpp(Pno    ,Pno    ,0.0);  

  Var_Carrier var;
  var.set_elasticity_functions(mat,grid,eid);
  
  Matrix<double> u(nne*nsd, 1), P(Pno, 1);
  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u.m_pdata[a*nsd+b] = r_e[a*ndofn+b];

    if(Pno==nne)
      P.m_pdata[a] = r_e[a*ndofn+nsd];
  }
  if(Pno==1)
  {  
    P(1) = alpha_1*fv->tf.P_n(eid+1,1) + 
           alpha_2*(fv->tf.P_np1(eid+1,1) + fv->tf.dP(eid+1,1));
  }  
  memset(lk,0,ndofe*ndofe*sizeof(double)); 
   
  Matrix<double> Np(Pno,1,0.0);
  Matrix<double> Nt(Vno,1,0.0);  
       
  for(int ip = 1; ip<=fe->nint; ip++)
  {
    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    fe->update_deformation_gradient(ndofn,u.m_pdata,var.mF);

    fe->elem_shape_function(ip,Pno, Np.m_pdata);
    fe->elem_shape_function(ip,Vno, Nt.m_pdata);
    
    double Tn = 0.0;
    double Pn = 0.0;

    for(int ia=1; ia<=Vno; ia++)
    {
      Tn += Nt(ia)*(alpha_1*fv->tf.V_n(eid+1, ia) + 
            alpha_2*(fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia)));
    }
    for(int ia=1; ia<=Pno; ia++)
      Pn += Np(ia)*P(ia);
      
    var.set_variables(Tn, Pn);      

    TF_Kup_ip(Kup, fe, var, Pno, Np, dt_alpha_1_minus_alpha);
    TF_Kpt_ip(Kpt, fe, var, Pno, Np, Vno, Nt, dt_alpha_1_minus_alpha);

    TF_Kuu_ip(Kuu, fe, var, dt_alpha_1_minus_alpha);
    TF_Ktt_ip(Ktt, fe, var, Vno,Nt, dt_alpha_1_minus_alpha);
  }
  
  Kpu.trans(Kup);
  Ktp.trans(Kpt);
  
  err += condense_K_3F_to_1F(lk, nne, nsd, Pno, Vno,
                             Kuu.m_pdata, Kut.m_pdata, Kup.m_pdata,
                             Ktu.m_pdata, Ktt.m_pdata, Ktp.m_pdata,
                             Kpu.m_pdata, Kpt.m_pdata, NULL);

  // check diagonal for zeros/nans
  for (int a = 0; a < nne; a++) {
    for (int b = 0; b < nsd; b++) {
      if ( !isnormal(lk[idx_K(a,b,a,b,nne,nsd)]) ) err++;
    }
  }
  
  return err;
}

/// compute element residual vector in quasi steady state.
/// Total Lagrangian based three-field mixed method for Hyperelasticity.
///
/// \param[in]  fe    finite element helper object
/// \param[out] f     computed element residual vector
/// \param[in]  r_e   nodal variabls(displacements) on the current element
/// \param[in]  grid  a mesh object
/// \param[in]  mat   a material object
/// \param[in]  fv    object for field variables
/// \return non-zero on internal error
int residuals_3f_el(FEMLIB *fe,
                    double *f,
                    double *r_e,
                    Grid *grid,
                    MaterialProperty *mat,
                    FieldVariables *fv)
{
  int err = 0;

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  int ndofe = nne*ndofn;
  int Pno   = fv->npres;
  int Vno   = fv->nVol;
  
  Matrix<double> Kut(nne*nsd,Vno, 0.0);
  Matrix<double> Kup(nne*nsd,Pno, 0.0);
  Matrix<double> Ktt(Vno    ,Vno, 0.0);
  Matrix<double> Ktp(Vno    ,Pno, 0.0);
  Matrix<double> Kpt(Pno    ,Vno, 0.0);  
  Matrix<double> Ru(nne*nsd,1,    0.0);
  Matrix<double> Rp(Pno,  1,    0.0);
  Matrix<double> Rt(Vno,   1,    0.0);  

  Var_Carrier var;
  var.set_elasticity_functions(mat,grid,eid);
  
  Matrix<double> u(nne*nsd, 1), P(Pno, 1);
  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u.m_pdata[a*nsd+b] = r_e[a*ndofn+b];

    if(Pno==nne)
      P.m_pdata[a] = r_e[a*ndofn+nsd];
  }
  if(Pno==1)
    P.m_pdata[0] = fv->tf.P_np1(eid+1,1) + fv->tf.dP(eid+1,1);
    
  memset(f,0,ndofe*sizeof(double)); 
   
  Matrix<double> Np(Pno,1,0.0);
  Matrix<double> Nt(Vno,1,0.0);  
       
  for(int ip = 1; ip<=fe->nint; ip++)
  {
    fe->elem_basis_V(ip);
    fe->update_shape_tensor();
    fe->update_deformation_gradient(ndofn,u.m_pdata,var.mF);

    fe->elem_shape_function(ip,Pno, Np.m_pdata);
    fe->elem_shape_function(ip,Vno, Nt.m_pdata);
    
    double Tn = 0.0;
    double Pn = 0.0;

    for(int ia=1; ia<=Vno; ia++)
      Tn += Nt(ia)*(fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia));

    for(int ia=1; ia<=Pno; ia++)
      Pn += Np(ia)*P(ia);
      
    var.set_variables(Tn, Pn);      
    
    TF_Kup_ip(Kup, fe, var, Pno, Np, -1.0);
    TF_Kpt_ip(Kpt, fe, var, Pno, Np, Vno, Nt, -1.0);
    TF_Ktt_ip(Ktt, fe, var, Vno, Nt, -1.0);
        
    TF_Ru_ip(Ru,fe,var);                  
    TF_Rp_ip(Rp,fe,var,Pno,Np);
    TF_Rt_ip(Rt,fe,var,Vno,Nt);    
  }
  
  Ktp.trans(Kpt);  
  
  err += condense_F_3F_to_1F(f, nne, nsd, Pno, Vno,
                             Ru.m_pdata, Rt.m_pdata, Rp.m_pdata,
                             Kut.m_pdata, Kup.m_pdata, Ktp.m_pdata, Ktt.m_pdata, Kpt.m_pdata);
  return err;
}

/// compute element residual vector in transient .
/// Total Lagrangian based three-field mixed method for Hyperelasticity.
///
/// \param[in]  fe    finite element helper object
/// \param[out] f     computed element stiffness matrix
/// \param[in]  r_e   nodal variabls(displacements) on the current element
/// \param[in]  grid  a mesh object
/// \param[in]  mat   a material object
/// \param[in]  fv    object for field variables
/// \param[in]  alpha mid point alpha
/// \param[in]  dts   time step size at t(n), t(n+1); dts[DT_N]   = t(n)   - t(n-1)
///                                                   dts[DT_NP1] = t(n+1) - t(n)
/// \param[in]  re_nma nodal variabls(displacements) on the current element: (1-alpha)*u(n-1) + alpha*u(n)
/// \param[in]  re_npa nodal variabls(displacements) on the current element: (1-alpha)*u(n)   + alpha*u(n+1)
/// \return non-zero on internal error
int residuals_3f_w_inertia_el(FEMLIB *fe,
                              double *f,
                              double *r_e,
                              Grid *grid,
                              MaterialProperty *mat,
                              FieldVariables *fv,
                              Solver *sol,
                              const double *dts,
                              double *re_nma,
                              double *re_npa)
{

  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  int Pno   = fv->npres;
  int Vno   = fv->nVol;  
  
  double alpha_1;
  double alpha_2;
  double dt_alpha_1_minus_alpha;
  double dt_alpha_1;
  double dt_alpha_2;

  if(sol->alpha<0)
  {
    alpha_1 = 0.0;
    alpha_2 = sol->alpha;
    dt_alpha_1_minus_alpha = -1.0;
    dt_alpha_1 = 1.0;
    dt_alpha_2 = 0.0;
  }
  else
  {
    alpha_1 = 1.0 - sol->alpha;
    alpha_2 = sol->alpha;
    dt_alpha_1_minus_alpha = dts[DT_NP1]*alpha_1*alpha_2;
    dt_alpha_1 = -dts[DT_NP1]*alpha_1;
    dt_alpha_2 = -dts[DT_N]*alpha_2;
  }

  Matrix<double> u_npa(nne*nsd, 1), u_nma(nne*nsd, 1), Pe_npa(Pno, 1), Pe_nma(Pno, 1);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
    {
      u_npa.m_pdata[a*nsd+b] = re_npa[a*ndofn+b];
      u_nma.m_pdata[a*nsd+b] = re_nma[a*ndofn+b];      
    }

    if(Pno==nne)
    {
      Pe_npa.m_pdata[a] = re_npa[a*ndofn+nsd];      
      Pe_nma.m_pdata[a] = re_nma[a*ndofn+nsd];
    }
  }

  if(Pno==1)
  {
    Pe_npa(1) = alpha_1*fv->tf.P_n(eid+1,1) + 
            alpha_2*( fv->tf.P_np1(eid+1,1) + fv->tf.dP(eid+1,1));
    Pe_nma(1) = alpha_1*fv->tf.P_nm1(eid+1,1) + alpha_2*fv->tf.P_n(eid+1,1);
            
  }
  
  Var_Carrier var_npa, var_nma;
  var_npa.set_elasticity_functions(mat,grid,eid);
  var_nma.set_elasticity_functions(mat,grid,eid);
  
  Matrix<double> Ru(  nne*nsd,1, 0.0);
  Matrix<double> Ru_npa(nne*nsd,1, 0.0);
  Matrix<double> Ru_nma(nne*nsd,1, 0.0);
  Matrix<double> Rp(  Pno,    1, 0.0);
  Matrix<double> Rp_npa(Pno,    1, 0.0);
  Matrix<double> Rp_nma(Pno,    1, 0.0);
  Matrix<double> Rt(  Vno,    1, 0.0);
  Matrix<double> Rt_npa(Vno,    1, 0.0);
  Matrix<double> Rt_nma(Vno,    1, 0.0);
  Matrix<double> Kut(nne*nsd,Vno, 0.0);
  Matrix<double> Kpt(Pno,    Vno, 0.0);
  Matrix<double> Ktt(Vno,    Vno, 0.0);
  Matrix<double> Kup(nne*nsd,Pno, 0.0);
  Matrix<double> Ktp(Vno,    Pno, 0.0);

  Matrix<double> Np(Pno,1,0.0), Nt(Vno, 1,0.0);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    fe->elem_basis_V(ip);
    fe->elem_shape_function(ip,Pno, Np.m_pdata);
    fe->elem_shape_function(ip,Vno, Nt.m_pdata);
    fe->update_shape_tensor();
    fe->update_deformation_gradient(ndofn,u_npa.m_pdata,var_npa.mF);
    fe->update_deformation_gradient(ndofn,u_nma.m_pdata,var_nma.mF);

    double Tnma = 0.0; // 1: n-1+alpha
    double Tnpa = 0.0; // 2: n+alpha
    double Pnma = 0.0;
    double Pnpa = 0.0;

    for(int ia=1; ia<=Vno; ia++)
    {
      Tnpa += Nt(ia)*(alpha_1*fv->tf.V_n(eid+1, ia) + 
             alpha_2*(fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia)));
             
      Tnma += Nt(ia)*(alpha_1*fv->tf.V_nm1(eid+1, ia) + alpha_2*fv->tf.V_n(eid+1, ia));             
    }
    for(int ia=1; ia<=Pno; ia++)
    {
      Pnpa += Np(ia)*Pe_npa(ia);
      Pnma += Np(ia)*Pe_nma(ia);
    }
        
    var_npa.set_variables(Tnpa, Pnpa);
    var_nma.set_variables(Tnma, Pnma);    
    
    TF_Kup_ip(Kup, fe, var_npa, Pno, Np, dt_alpha_1_minus_alpha);
    TF_Kpt_ip(Kpt, fe, var_npa, Pno, Np, Vno, Nt, dt_alpha_1_minus_alpha);
    TF_Ktt_ip(Ktt, fe, var_npa, Vno, Nt, dt_alpha_1_minus_alpha);

    TF_Ru_ip(Ru_npa,fe,var_npa);
    TF_Rp_ip(Rp_npa,fe,var_npa,Pno,Np);
    TF_Rt_ip(Rt_npa,fe,var_npa,Vno,Nt);
            
    TF_Ru_ip(Ru_nma,fe,var_nma);
    TF_Rp_ip(Rp_nma,fe,var_nma,Pno,Np);
    TF_Rt_ip(Rt_nma,fe,var_nma,Vno,Nt);
  }
  Ktp.trans(Kpt);

  for(int a=0; a<nne*nsd; a++)
    Ru.m_pdata[a] = dt_alpha_1*Ru_npa.m_pdata[a] + dt_alpha_2*Ru_nma.m_pdata[a];

  for(int a=0; a<Vno; a++)
    Rt.m_pdata[a] = dt_alpha_1*Rt_npa.m_pdata[a] + dt_alpha_2*Rt_nma.m_pdata[a];

  for(int a=0; a<Pno; a++)
    Rp.m_pdata[a] = dt_alpha_1*Rp_npa.m_pdata[a] + dt_alpha_2*Rp_nma.m_pdata[a];

  Ktp.trans(Kpt);  
  condense_F_3F_to_1F(f, nne, nsd, Pno, Vno,
                      Ru.m_pdata, Rt.m_pdata, Rp.m_pdata,
                      Kut.m_pdata, Kup.m_pdata, Ktp.m_pdata, Ktt.m_pdata, Kpt.m_pdata);
  return 0;
} 

/// compute and update increments of prssure and volume for qusi-steady state
/// in an element
///
/// \param[in]  fe       finite element helper object
/// \param[in]  grid     a mesh object
/// \param[in]  material a material object
/// \param[in]  fv       object for field variables, fv->tf will be updated
/// \param[in]  r_e      nodal variabls(displacements + pressure (if ndofn = 4)) 
///                      on the current element
/// \param[in]  d_u      increments of displacement (updated by NR before calling this function)
void evaluate_PT_el(FEMLIB *fe,
                    Grid *grid,
                    MaterialProperty *material,
                    FieldVariables *fv,
                    Matrix<double> &r_e,
                    Matrix<double> &du)
{
  int eid = fe->curt_elem_id;
  Var_Carrier var;
  var.set_elasticity_functions(material,grid,eid);
  
  int ndofn = fv->ndofn;
  int Pno   = fv->npres;
  int Vno   = fv->nVol;
  int nsd   = fe->nsd;
  int nne   = fe->nne;
  
  Matrix<double> P(Pno, 1, 0.0);
  Matrix<double> u(nne*nsd, 1, 0.0);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u.m_pdata[a*nsd+b] = r_e.m_pdata[a*ndofn+b];

    if(Pno==nne)
      P.m_pdata[a] = r_e.m_pdata[a*ndofn+nsd];
  }
  if(Pno==1)
    P.m_pdata[0] = fv->tf.P_np1(eid+1,1) + fv->tf.dP(eid+1,1);

  Matrix<double> Rp(Pno, 1, 0.0), Rt(Vno,  1, 0.0);

  Matrix<double> Ktu(Vno, nne*nsd,0.0); // will be transposed: Kut -> Ktu
  Matrix<double> Ktt(Vno   ,Vno  ,0.0);
  Matrix<double> Kpu(nne*nsd,Pno ,0.0);  // will be transposed: Kpu -> Kup 
  Matrix<double> Kpt(Pno  ,Vno   ,0.0);
  Matrix<double> Ktp(Vno   ,Pno  ,0.0);

  Matrix<double> Np(Pno,1,0.0), Nt(Vno, 1,0.0);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    fe->elem_basis_V(ip);
    fe->elem_shape_function(ip,Pno,Np.m_pdata);
    fe->elem_shape_function(ip,Vno, Nt.m_pdata);
    fe->update_shape_tensor();
    fe->update_deformation_gradient(ndofn,u.m_pdata,var.mF);    

    double Tn = 0.0;
    double Pn = 0.0;
    
    for(int ia=1; ia<=Vno; ia++)
      Tn += Nt(ia)*(fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia));

    for(int ia=1; ia<=Pno; ia++)
      Pn += Np(ia)*P(ia);

    var.set_variables(Tn, Pn);
    Matrix<double> S,F;
    F.use_reference(DIM_3,DIM_3,var.mF);
    S.use_reference(DIM_3,DIM_3,var.mS);         

    TF_Kup_ip(Kpu,fe,var,Pno,Np,-1.0);
    TF_Kpt_ip(Kpt,fe,var,Pno,Np,Vno,Nt,-1.0);
    TF_Ktt_ip(Ktt,fe,var,Vno,Nt,-1.0);
    
    TF_Rp_ip(Rp,fe,var,Pno,Np);
    TF_Rt_ip(Rt,fe,var,Vno,Nt);        
  }

  Kpu.trans();
  Ktp.trans(Kpt);

  Matrix<double> d_theta(Vno, 1, 0.0), dP(Pno, 1, 0.0);
  
  compute_d_theta_dP_test(d_theta,dP, 
                            du,Rt,Rp,Kpu,Ktu,Ktp,Ktt,Kpt);

  for(int ia=1; ia<=Pno; ia++)
    fv->tf.ddP(eid+1,ia) = dP(ia);

  for(int ia=1; ia<=Vno; ia++)
    fv->tf.ddV(eid+1,ia) = d_theta(ia);
}

/// compute and update increments of prssure and volume for transient
/// in an element
///
/// \param[in]  fe     finite element helper object
/// \param[in]  grid   a mesh object
/// \param[in]  mat    a material object
/// \param[in]  fv     object for field variables, fv->tf will be updated
/// \param[in]  r_e    nodal variabls(displacements + pressure (if ndofn = 4)) 
///                    on the current element
/// \param[in]  d_u    increments of displacement (updated by NR before calling this function)
/// \param[in]  dts   time step size at t(n), t(n+1); dts[DT_N]   = t(n)   - t(n-1)
///                                                   dts[DT_NP1] = t(n+1) - t(n)
/// \param[in]  re_nma nodal variabls(displacements) on the current element: (1-alpha)*u(n-1) + alpha*u(n)
/// \param[in]  re_npa nodal variabls(displacements) on the current element: (1-alpha)*u(n)   + alpha*u(n+1)
/// \param[in]  alpha  mid point alpha
void evaluate_PT_w_inertia_el(FEMLIB *fe,
                              Grid *grid,
                              MaterialProperty *mat,
                              FieldVariables *fv,
                              Matrix<double> &r_e,                              
                              Matrix<double> &du,
                              const double *dts,
                              Matrix<double> &re_nma,
                              Matrix<double> &re_npa,
                              double alpha)
{
  int eid = fe->curt_elem_id;
  int nsd = fe->nsd;
  int nne = fe->nne;
  int ndofn = fv->ndofn;
  int Pno   = fv->npres;
  int Vno   = fv->nVol;  
    
  double alpha_1;
  double alpha_2;
  double dt_alpha_1_minus_alpha;
  double dt_alpha_1;
  double dt_alpha_2;

  if(alpha<0)
  {
    alpha = 1.0;
    alpha_1 = 0.0;
    alpha_2 = alpha;
    dt_alpha_1_minus_alpha = -1.0;
    dt_alpha_1 = 1.0;
    dt_alpha_2 = 0.0;
  }
  else
  {
    alpha_1 = 1.0 - alpha;
    alpha_2 = alpha;    
    dt_alpha_1_minus_alpha = dts[DT_NP1]*alpha_1*alpha_2;
    dt_alpha_1 = -dts[DT_NP1]*alpha_1;
    dt_alpha_2 = -dts[DT_N]*alpha_2;    
  }
  
  Matrix<double> u_npa(nne*nsd, 1), u_nma(nne*nsd, 1), Pe_npa(Pno, 1), Pe_nma(Pno, 1);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
    {
      u_npa.m_pdata[a*nsd+b] = re_npa.m_pdata[a*ndofn+b];
      u_nma.m_pdata[a*nsd+b] = re_nma.m_pdata[a*ndofn+b];
    }

    if(Pno==nne)
    {
      Pe_npa.m_pdata[a] = re_npa.m_pdata[a*ndofn+nsd];
      Pe_nma.m_pdata[a] = re_nma.m_pdata[a*ndofn+nsd];
    }
  }

  if(Pno==1)
  {
    Pe_npa(1) = alpha_1*fv->tf.P_n(eid+1,1) + 
            alpha_2*( fv->tf.P_np1(eid+1,1) + fv->tf.dP(eid+1,1));
    Pe_nma(1) = alpha_1*fv->tf.P_nm1(eid+1,1) + alpha_2*fv->tf.P_n(eid+1,1);
            
  }
  
  Var_Carrier var_npa, var_nma;
  var_npa.set_elasticity_functions(mat,grid,eid);
  var_nma.set_elasticity_functions(mat,grid,eid);
  
  Matrix<double> Rp(    Pno, 1, 0.0);
  Matrix<double> Rp_npa(Pno, 1, 0.0);
  Matrix<double> Rp_nma(Pno, 1, 0.0);
  Matrix<double> Rt(    Vno, 1, 0.0);
  Matrix<double> Rt_npa(Vno, 1, 0.0);
  Matrix<double> Rt_nma(Vno, 1, 0.0);
  
  Matrix<double> Ktu(Vno    ,nne*nsd,0.0);
  Matrix<double> Ktt(Vno    ,Vno    ,0.0);
  Matrix<double> Kup(nne*nsd,Pno    ,0.0);
  Matrix<double> Kpu(Pno    ,nne*nsd,0.0);
  Matrix<double> Kpt(Pno    ,Vno    ,0.0);
  Matrix<double> Ktp(Vno    ,Pno    ,0.0); 

  Matrix<double> Np(Pno,1,0.0),Nt(Vno, 1,0.0);

  for(int ip = 1; ip<=fe->nint; ip++)
  {
    fe->elem_basis_V(ip);
    fe->elem_shape_function(ip,Pno, Np.m_pdata);
    fe->elem_shape_function(ip,Vno, Nt.m_pdata);
    fe->update_shape_tensor();
    fe->update_deformation_gradient(ndofn,u_npa.m_pdata,var_npa.mF);
    fe->update_deformation_gradient(ndofn,u_nma.m_pdata,var_nma.mF);

    double Tnpa = 0.0; // 1: n-1+alpha
    double Tnma = 0.0; // 2: n+alpha
    double Pnpa = 0.0;
    double Pnma = 0.0;

    for(int ia=1; ia<=Vno; ia++)
    {
      Tnpa += Nt(ia)*(alpha_1*fv->tf.V_n(eid+1, ia) + 
             alpha_2*(fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia)));
      Tnma += Nt(ia)*(alpha_1*fv->tf.V_nm1(eid+1, ia) + alpha_2*fv->tf.V_n(eid+1, ia));             
    }
    for(int ia=1; ia<=Pno; ia++)
    {
      Pnpa += Np(ia)*Pe_npa(ia);
      Pnma += Np(ia)*Pe_nma(ia);
    }
    
    var_npa.set_variables(Tnpa, Pnpa);
    var_nma.set_variables(Tnma, Pnma);

    TF_Kup_ip(Kup, fe, var_npa, Pno, Np, dt_alpha_1_minus_alpha);
    TF_Kpt_ip(Kpt, fe, var_npa, Pno, Np, Vno, Nt, dt_alpha_1_minus_alpha);
    TF_Ktt_ip(Ktt, fe, var_npa, Vno, Nt, dt_alpha_1_minus_alpha);

    TF_Rp_ip(Rp_npa,fe,var_npa,Pno,Np);
    TF_Rt_ip(Rt_npa,fe,var_npa,Vno,Nt);

    TF_Rp_ip(Rp_nma,fe,var_nma,Pno,Np);
    TF_Rt_ip(Rt_nma,fe,var_nma,Vno,Nt);
  }

  Kpu.trans(Kup);
  Ktp.trans(Kpt);

  for(int a=0; a<Vno; a++)
    Rt.m_pdata[a] = dt_alpha_1*Rt_nma.m_pdata[a] + dt_alpha_2*Rt_npa.m_pdata[a];

  for(int a=0; a<Pno; a++)
    Rp.m_pdata[a] = dt_alpha_1*Rp_nma.m_pdata[a] + dt_alpha_2*Rp_npa.m_pdata[a];

  Matrix<double> d_theta(Vno, 1, 0.0), dP(Pno, 1, 0.0);
  
  compute_d_theta_dP_test(d_theta,dP, 
                          du,Rt,Rp,Kpu,Ktu,Ktp,Ktt,Kpt);


  for(int ia=1; ia<=Pno; ia++)
    fv->tf.ddP(eid+1,ia) = dP(ia);

  for(int ia=1; ia<=Vno; ia++)
    fv->tf.ddV(eid+1,ia) = d_theta(ia);
}

/// compute increments of prssure and volume in an element
///
/// \param[out] d_theta computed increment of volume
/// \param[out] dP      computed increment of pressure
/// \param[in]  du      increments of displacement (updated by NR before calling this function)
/// \param[in]  Rt      residual for volume     
/// \param[in]  Rp      residual for pressure
/// \param[in]  Kpu     stiffness for pu
/// \param[in]  Ktu     stiffness for tu
/// \param[in]  Ktp     stiffness for tp
/// \param[in]  Ktt     stiffness for tt
/// \param[in]  Kpt     stiffness for pt
/// \return non-zero on internal error 
int compute_d_theta_dP_test(Matrix<double> &d_theta,
                            Matrix<double> &dP, 
                            Matrix<double> &du,
                            Matrix<double> &Rt, 
                            Matrix<double> &Rp, 
                            Matrix<double> &Kpu, 
                            Matrix<double> &Ktu, 
                            Matrix<double> &Ktp, 
                            Matrix<double> &Ktt, 
                            Matrix<double> &Kpt)
{
  int err = 0;
    	  
	Matrix<double> KptI;
  KptI.inv(Kpt);


//	Matrix<double> KptIFp, KptIKpu_du;
//  KptIFp.prod(KptI,Rp);
//  KptIKpu_du.prod(KptI, Kpu, du);
//                
//  for(int ia=1; ia<=d_theta.m_row; ia++)
//    d_theta(ia) = -KptIFp(ia) - KptIKpu_du(ia);   

  Matrix<double> Kpu_du;
  Kpu_du.prod(Kpu,du);
  Kpu_du.add(Rp);

  d_theta.prod(KptI,Kpu_du);
  d_theta.prod(-1.0);   

        
	Matrix<double> KtpI, KtpIFt, KtpIKtu_du;
	Matrix<double> KtpIKtt_d_theta;
  KtpI.inv(Ktp);
  
  KtpIFt.prod(KtpI, Rt);
  KtpIKtu_du.prod(KtpI, Ktu,du);

  KtpIKtt_d_theta.prod(KtpI, Ktt, d_theta);  
  
  for(int ia=1; ia<=dP.m_row; ia++)
    dP(ia) = -KtpIFt(ia) - KtpIKtu_du(ia) - KtpIKtt_d_theta(ia);
    
  return err;
}

/// compute ouput variables e.g. effective stress and strain
///
/// Visit each element and compute output variables according to the element model type.
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv array of field variable object
/// \param[in] load object for loading
/// \param[in] mp_id mutiphysics id
/// \param[in] dts   time step size at t(n), t(n+1); dts[DT_N]   = t(n)   - t(n-1)
///                                                  dts[DT_NP1] = t(n+1) - t(n)
/// \param[in] alpha mid point rule alpha
/// \return non-zero on internal error
int update_3f_NR(Grid *grid,
                 MaterialProperty *mat,
                 FieldVariables *fv,
                 LoadingSteps *load,
                 const PGFem3D_opt *opts,
                 int mp_id,
                 const double *dts,
                 double alpha)
{

  const int mat_id = grid->element[0].mat[2];
  double rho = mat->hommat[mat_id].density;
  long include_inertia = 1;

  if(fabs(rho)<MIN_DENSITY)
    include_inertia = 0;
    

  int total_Lagrangian = 1;
  
  Node *node = grid->node;
  Element *elem = grid->element;

  int ndofn = fv->ndofn;
  SUPP sup = load->sups[mp_id];
      

  for (int eid=0;eid<grid->ne;eid++)
  {
    FEMLIB fe(eid,elem,node,1,total_Lagrangian);
    int nne   = fe.nne;
    int nsd   = fe.nsd;
    int ndofe = nne*ndofn;

    Matrix<long> cn(ndofe,1);
    long *nod = fe.node_id.m_pdata;

    Matrix<double> dr_e(ndofe,1), re_np1(ndofe, 1), du(nne*nsd,1), u(nne*nsd, 1);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn.m_pdata,mp_id);
    
    // get the deformation on the element
    def_elem_total(cn.m_pdata,ndofe,fv->u_np1,fv->d_u,elem,node,sup,re_np1.m_pdata);    
    def_elem(cn.m_pdata,ndofe,fv->dd_u,elem,node,dr_e.m_pdata,sup,2);

    for(int a=0;a<nne;a++)
    {
      for(int b=0; b<nsd;b++)
        du.m_pdata[a*nsd+b] = dr_e.m_pdata[a*ndofn+b];
    }
    
    if(include_inertia && (0<alpha && alpha<1))
    {
      Matrix<double> u_n(ndofe,1), u_nm1(ndofe,1);
      Matrix<double> re_npa(ndofe,1), re_nma(ndofe,1);
      
      for (long I=0;I<nne;I++)
      {
        for(long J=0; J<nsd; J++)
        {
            u_n.m_pdata[I*ndofn + J] =   fv->u_n[nod[I]*ndofn + J];
          u_nm1.m_pdata[I*ndofn + J] = fv->u_nm1[nod[I]*ndofn + J];
        }
      }

      mid_point_rule(re_nma.m_pdata, u_nm1.m_pdata, u_n.m_pdata, alpha, ndofe);
      mid_point_rule(re_npa.m_pdata, u_n.m_pdata,re_np1.m_pdata, alpha, ndofe);      
      evaluate_PT_w_inertia_el(&fe, grid, mat, fv, re_np1, du, dts, re_nma, re_npa, alpha);      
    }
    else
      evaluate_PT_el(&fe,grid,mat,fv,re_np1,du);
  }
  return 0;  
}


/// Update variables during Newton Raphson iterations
///
/// \param[in] grid  a mesh object
/// \param[in] mat   a material object
/// \param[in] fv    object for field variables
/// \param[in] load  object for loading
/// \param[in] mp_id mutiphysics id
/// \param[in] dt    time step size
/// \param[in] t     time
/// \param[in] mpi_comm MPI_COMM_WORLD
void update_3f_output_variables(Grid *grid,
                                MaterialProperty *mat,
                                FieldVariables *fv,
                                LoadingSteps *load,
                                const int mp_id,
                                const double dt,
                                const double t,
                                MPI_Comm mpi_comm)
{
  int ndofn = fv->ndofn;
  int Pno   = fv->npres;
  int Vno   = fv->nVol;
  
  SIG *sig = fv->sig;
  EPS *eps = fv->eps;
  Node *node = grid->node;
  SUPP sup = load->sups[mp_id];
  Element *elem = grid->element;
  
  int idx[6];
  idx[0] = idx_2(0,0); //XX
  idx[1] = idx_2(1,1); //YY
  idx[2] = idx_2(2,2); //ZZ
  idx[3] = idx_2(1,2); //YZ
  idx[4] = idx_2(0,2); //XZ
  idx[5] = idx_2(0,1); //XY
  
  double idnty[9] = {1.0,0.0,0.0,
                     0.0,1.0,0.0,
                     0.0,0.0,1.0};
  TensorA<2> I(idnty);
  
  for (int eid=0;eid<grid->ne;eid++)
  {
    Var_Carrier var;
    var.set_elasticity_functions(mat,grid,eid);
    
    memset(sig[eid].el.o,0,6*sizeof(double));
    memset(eps[eid].el.o,0,6*sizeof(double));
      
    FEMLIB fe(eid,elem,node,0,1);
    int nne   = fe.nne;
    int nsd   = fe.nsd;
    int ndofe = nne*ndofn;

    Matrix<long> cn(ndofe,1);
    long *nod = fe.node_id.m_pdata;

    Matrix<double> r_e(ndofe, 1);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn.m_pdata,mp_id);
    def_elem_total(cn.m_pdata,ndofe,fv->u_np1,fv->d_u,elem,node,sup,r_e.m_pdata);

    Matrix<double> P(Pno, 1, 0.0);
    Matrix<double> u(nne*nsd, 1, 0.0);

    for(int a=0;a<nne;a++)
    {
      for(int b=0; b<nsd;b++)
        u.m_pdata[a*nsd+b] = r_e.m_pdata[a*ndofn+b];

      if(Pno==nne)
        P.m_pdata[a] = r_e.m_pdata[a*ndofn+nsd];
    }
    if(Pno==1)
      P.m_pdata[0] = fv->tf.P_np1(eid+1,1) + fv->tf.dP(eid+1,1);
      
    Matrix<double> Np(Pno,1,0.0), Nt(Vno, 1,0.0);
    
    double V = 0.0;
    for(int ip = 1; ip<=fe.nint; ip++)
    {
      fe.elem_basis_V(ip);
      fe.elem_shape_function(ip,Pno,Np.m_pdata);
      fe.elem_shape_function(ip,Vno, Nt.m_pdata);
      fe.update_shape_tensor();
      fe.update_deformation_gradient(ndofn,u.m_pdata,var.mF);
      
      V += fe.detJxW;
      
      double Tn = 0.0;
      double Pn = 0.0;

      for(int ia=1; ia<=Vno; ia++)
        Tn += Nt(ia)*(fv->tf.V_np1(eid+1, ia) + fv->tf.dV(eid+1, ia));

      for(int ia=1; ia<=Pno; ia++)
        Pn += Np(ia)*P(ia);
        
      var.set_variables(Tn, Pn);
      
      TensorA<2> F(var.mF), FI(var.mFI);  
      Tensor<2> C = var.factor*F(k,i)*F(k,j);
      Tensor<2> b = var.factor*F(i,k)*F(j,k);
      Tensor<2> CI, bI, S_bar, S, sigma;
      
      inv(C, CI);
      inv(b, bI);

      // Elastic Green Lagrange strain
      Tensor<2> E = 0.5*(C(i,j) - I(i,j));
      
      // Compute the logarithmic strain e = 1/2(I - inv(FF'))
      Tensor<2> e = 0.5*(I(i,j) - bI(i,j));
     
      // compute stresses
      double Up = var.compute_dUdJ(var.theta);
      var.compute_PK2(C.data, S_bar.data);
  
      // compute PKII    
      S = S_bar(i,j) + Up*var.theta*CI(i,j);

      // Compute Cauchy stress (theta)^-1 F S F'
      sigma = 1.0/var.theta*F(i,k)*var.factor*S(k,l)*F(j,l);
      
      for(int ia=0; ia<6; ia++)
      {
        int id = idx[ia];
        sig[eid].il[ip-1].o[ia] = S.data[id];
        eps[eid].il[ip-1].o[ia] = E.data[id];

        sig[eid].el.o[ia] += fe.detJxW*sigma.data[id];
        
        // engineering strain
        if(ia<3)
          eps[eid].el.o[ia] += fe.detJxW*e.data[id];
        else      
          eps[eid].el.o[ia] += 2.0*fe.detJxW*e.data[id]; 
      }
    }

    for(int ia=0; ia<6; ia++)
    {
      sig[eid].el.o[ia] /=V;
      eps[eid].el.o[ia] /=V;
    }
  }                                                        
}

/// compute and set initial conditions for three field mixed method
///
/// \param[in] grid an object containing all mesh info
/// \param[in] mat a material object
/// \param[in,out] fv array of field variable object
/// \return non-zero on internal error
void compute_3f_initial_conditions(Grid *grid,
                                   MaterialProperty *mat,
                                   FieldVariables *fv)
{

  int total_Lagrangian = 1;  
  Node *node = grid->node;
  Element *elem = grid->element;

  int ndofn = fv->ndofn;      

  for (int eid=0;eid<grid->ne;eid++)
  {
    FEMLIB fe(eid,elem,node,1,total_Lagrangian);
    int nne   = fe.nne;
    int nsd   = fe.nsd;
    
    long *nod = fe.node_id.m_pdata;

    Var_Carrier var;
    var.set_elasticity_functions(mat,grid,eid);    

    Matrix<double> u_nm1(nne*nsd, 1), u_n(nne*nsd, 1);
    
    for (long I=0;I<nne;I++){
      for(long J=0; J<nsd; J++){
        u_n.m_pdata[I*ndofn + J] =   fv->u_n[nod[I]*ndofn + J];
        u_nm1.m_pdata[  I*ndofn + J] = fv->u_nm1[nod[I]*ndofn + J];
      }
    }

    double J_n    = 0.0;
    double J_nm1  = 0.0;
    double Up_n   = 0.0;
    double Up_nm1 = 0.0;
    double V = 0.0;
    for(int ip=1; ip<=fe.nint; ip++)
    {
      Tensor<2> Fn;
      Tensor<2> Fnm1;
      fe.elem_basis_V(ip);
      fe.update_shape_tensor();
      fe.update_deformation_gradient(nsd,u_n.m_pdata,    Fn.data);
      fe.update_deformation_gradient(nsd,u_nm1.m_pdata,Fnm1.data);
    
      double Jip_n   = ttl::det(Fn);
      double Jip_nm1 = ttl::det(Fnm1);
        
      V      += fe.detJxW;
      J_n    += fe.detJxW*Jip_n;
      J_nm1  += fe.detJxW*Jip_nm1;
      Up_n   += fe.detJxW*var.compute_dUdJ(Jip_n);
      Up_nm1 += fe.detJxW*var.compute_dUdJ(Jip_nm1);      
    }
    if(fv->npres == 1)
    {  
      fv->tf.P_np1(eid+1, 1) = fv->tf.P_n(eid+1, 1) = Up_n/V;
      fv->tf.P_nm1(eid+1, 1) = Up_nm1/V;
    }
    
    fv->tf.V_np1(eid+1, 1) = fv->tf.V_n(eid+1, 1) = J_n/V;
    fv->tf.V_nm1(eid+1, 1) = J_nm1/V;
  }
}