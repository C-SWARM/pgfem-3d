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

    };
    void set_variables(double theta_in, double P_in, bool compute_elasticity_tenosr = true)
    {
      TensorA<2> F(mF), FI(mFI);
      inv(F,FI);
      J = ttl::det(F);
      theta  = theta_in;
      P      = P_in;
      factor = pow(theta_in/J, 2.0/3.0);
      
      Tensor<2> C;
      C = factor*F(k,i)*F(k,j);
      
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
    };

    double compute_d2UdJ2(double J)
    {
      double Up = 0.0;
      UPP(J, e_hommat, &Up);
      return kappa*Up;
    };
    void compute_PK2(double *C, double *S)
    {
      Stress(C,e_hommat,S);
    };
    void compute_L(double *C, double *L)
    {
      Stiffness(C,e_hommat,L);
    };

    template<class T1, class T2> void compute_A_of_du(T1 &A, T2 &du)
    {
      TensorA<2> F(mF), FI(mFI);
      A = {};
      
      Tensor<2> duTF = du(k,i)*F(k,j);
      Tensor<2> sduTF;
      
      symm(duTF, sduTF);
      double FITdu = FI(j,i)*du(i,j); 
      A(i,j) = sduTF(i,j) - 1.0/3.0*FITdu*F(k,i)*F(k,j);
    };
    
    template<class T1> double compute_B_of_du(T1 &du)
    {
      TensorA<2> FI(mFI);
      double B = -1.0/3.0*factor*du(j,i)*FI(i,j);
      return B;
    };
};

void print_diff(Matrix<double> &A, Matrix<double> &B)
{
  for(int ia=0; ia<A.m_row*A.m_col; ia++)
    printf("%d: %e %e %e\n", ia, A.m_pdata[ia], B.m_pdata[ia], A.m_pdata[ia] - B.m_pdata[ia]); 
};

int compute_d_theta_dP_test(Matrix<double> &d_theta,
                            Matrix<double> &dP, 
                            Matrix<double> &du,
                            Matrix<double> &ft, 
                            Matrix<double> &fp, 
                            Matrix<double> &Kpu, 
                            Matrix<double> &Ktu, 
                            Matrix<double> &Ktp, 
                            Matrix<double> &Ktt, 
                            Matrix<double> &Kpt);
                            
              
void TF_Kuu_ip(Matrix<double> &Kuu,
               FEMLIB *fe,
               Var_Carrier &var,
               double dt_alpha_1_minus_alpha)
{
  int nne = fe->nne;
  int nsd = fe->nsd;
  TensorA<2> S(var.mS), F(var.mF), FI(var.mFI);
  TensorA<4> L(var.mL);
  double factor = var.factor;
  
  for(int ia=0; ia<nne; ia++)
  {
    for(int ib=0; ib<nsd; ib++)
    {
      int id_ab = ia*nsd + ib + 1;
      TensorA<2> Grad_du(fe->ST + idx_4_gen(ia,ib,0,0,nne,nsd,nsd,nsd));
      Tensor<2> A_ab;
      var.compute_A_of_du(A_ab, Grad_du);
      double SA = S(i,j)*A_ab(i,j);
      double FIT_Grad_du = FI(j,i)*Grad_du(i,j);
      
      for(int iw=0; iw<nne; iw++)
      {
        for(int ig=0; ig<nsd; ig++)
        { 
          int id_wg = iw*nsd + ig + 1;
          // eq. 1
          TensorA<2> Grad_tu(fe->ST + idx_4_gen(iw,ig,0,0,nne,nsd,nsd,nsd));
          double B = var.compute_B_of_du(Grad_tu);

          // eq. 2
          Tensor<2> A_wg;
          var.compute_A_of_du(A_wg, Grad_tu);         
          double ALA = factor*factor*A_ab(i,j)*L(i,j,k,l)*A_wg(k,l);
          
          // eq. 3 - 1
          Tensor<2> duTtu = Grad_du(k,i)*Grad_tu(k,j);
          Tensor<2> sduTtu;
          symm(duTtu, sduTtu);

          // eq. 3 - 2          
          Tensor<2> FIT_Grad_tu_FIT = FI(k,i)*Grad_tu(l,k)*FI(j,l);
          double FIT_Grad_tu_FIT_du = FIT_Grad_tu_FIT(i,j)*Grad_du(i,j);
          
          Tensor<2> FF = 1.0/3.0*FIT_Grad_tu_FIT_du*F(k,i)*F(k,j);
          
          // eq. 3 - 3
          Tensor<2> tuTF = Grad_tu(k,i)*F(k,j);
          Tensor<2> stuTF;
          symm(tuTF,stuTF);

          // eq. 3
          double S_all = factor*S(i,j)*(sduTtu(i,j) + FF(i,j)
                                          -2.0/3.0*FIT_Grad_du*stuTF(i,j));
                                          
          // eq. 4
          double FIT_Grad_tu = FI(j,i)*Grad_tu(i,j);

          Kuu(id_ab,id_wg) += -dt_alpha_1_minus_alpha*fe->detJxW*(2.0*B*SA + ALA + S_all
                                     + var.P*var.J*(FIT_Grad_du*FIT_Grad_tu - FIT_Grad_tu_FIT_du));
        }
      }
    }
  }
  
}

void TF_Kut_ip(Matrix<double> &Kut,
               FEMLIB *fe,
               Var_Carrier &var,               
               double dt_alpha_1_minus_alpha,
               int Vno,
               Matrix<double> &Nt)
{
  int nne = fe->nne;
  int nsd = fe->nsd;
  TensorA<2> S(var.mS), F(var.mF);
  TensorA<4> L(var.mL);
  double factor = var.factor;
  Tensor<2> FFL = F(o,i)*F(o,j)*L(i,j,k,l);
  
  for(int ia=0; ia<nne; ia++)
  {
    for(int ib=0; ib<nsd; ib++)
    {
      int id_ab = ia*nsd + ib + 1;
      TensorA<2> Grad_du(fe->ST + idx_4_gen(ia,ib,0,0,nne,nsd,nsd,nsd));
      Tensor<2> A_ab;
      var.compute_A_of_du(A_ab, Grad_du);
      double SA = S(i,j)*A_ab(i,j);
      double FFLA = FFL(i,j)*A_ab(i,j);
      
      for(int iw=1; iw<=Vno; iw++)
      { 
        Kut(id_ab,iw) += -dt_alpha_1_minus_alpha*fe->detJxW*Nt(iw)*factor/var.theta/3.0*(2.0*SA + factor*FFLA);
      }
    }
  }
  
}

void TF_Ktt_ip(Matrix<double> &Ktt,
               FEMLIB *fe,
               Var_Carrier &var,               
               double dt_alpha_1_minus_alpha,
               int Vno,
               Matrix<double> &Nt)
{
  TensorA<2> S(var.mS), F(var.mF);
  TensorA<4> L(var.mL);
  double factor = var.factor;
  Tensor<2> C = F(k,i)*F(k,j);
  double SFF = S(i,j)*C(i,j);
  double FFLFF = C(i,j)*L(i,j,k,l)*C(k,l);
  
  for(int ia=1; ia<=Vno; ia++)
  {
    for(int iw=1; iw<=Vno; iw++)
    { 
      Ktt(ia,iw) += -dt_alpha_1_minus_alpha*fe->detJxW*Nt(ia)*Nt(iw)*
                      (1.0/9.0*factor/var.theta/var.theta*(-SFF + factor*FFLFF) + var.Upp);
    }
  }  
}

void TF_Rt_ip(double *Ktt,
              FEMLIB *fe,
              Var_Carrier &var,
              int Vno,
              Matrix<double> &Nt)
{
}              

void add_3F_Kup_ip(double *K,
                   int nne, int npres, double *ST, double *F, double jj, double wt, double *Np,
                   double dt_alpha_1_minus_alpha)
{
  int nsd = 3;

  double *AA  = aloc1(9);
  double *C   = aloc1(9);
  double *C_I = aloc1(9);

  double Jn = det3x3(F);
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              3,3,3,1.0,F,3,F,3,0.0,C,3);

  inverse(C,3,C_I);

  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      const double* const ptrST_ab = &ST[idx_4_gen(a,b,0,0,nne,nsd,nsd,nsd)];

      /* AA = F'Grad(del u) */
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                  3,3,3,1.0,F,3,ptrST_ab,3,0.0,AA,3);

      /* C_IFdu = C_I:AA*/
      double C_IFdu = cblas_ddot(9,C_I,1,AA,1);

      for(int w=0; w<npres; w++)
      {
        int idx_up = idx_K_gen(a,b,w,0,nne,nsd,npres,1);
        K[idx_up] += -dt_alpha_1_minus_alpha*jj*wt*Jn*C_IFdu*Np[w];
      }
    }
  }

  free(AA);
  free(C);
  free(C_I);
}

void add_3F_Ktt_ip(double *K, int nVol, double jj, double wt, double *Nt, double Upp,
                   double dt_alpha_1_minus_alpha)
{
  for(int a=0; a<nVol; a++)
  {
    for(int b=0; b<nVol; b++)
      K[idx_K(a,0,b,0,nVol,1)] += -dt_alpha_1_minus_alpha*Upp*(jj*wt*Nt[a]*Nt[b]);
  }
}
void add_3F_Ktp_ip(double *K,
                   int nVol, int npres, double jj, double wt, double *Nt, double *Np,
                   double dt_alpha_1_minus_alpha)
{

  for(int a=0; a<nVol; a++)
  {
    for(int b=0; b<npres; b++)
      K[idx_K_gen(a,0,b,0,nVol,1,npres,1)] += dt_alpha_1_minus_alpha*jj*wt*Nt[a]*Np[b];
  }
}
void add_3F_Kpu_ip(double *K,
                   int nne, int npres, double *ST, double *F, double jj, double wt, double *Np,
                   double dt_alpha_1_minus_alpha)
{
  int nsd = 3;
  double *AA  = aloc1(9);
  double *C   = aloc1(9);
  double *C_I = aloc1(9);
  double *F_I = aloc1(9);

  double Jn = det3x3(F);
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              3,3,3,1.0,F,3,F,3,0.0,C,3);

  inverse(F,3,F_I);
  inverse(C,3,C_I);

  for(int a=0; a<npres; a++)
  {
    int b = 0;
    for(int w=0;w<nne; w++)
    {

      for(int g=0;g<nsd;g++)
      {
        const double * const ptrST_wg = &ST[idx_4_gen(w,g,0,0,
                                                      nne,nsd,nsd,nsd)];

        /* AA = F'Grad(du) */
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                    3,3,3,1.0,F,3,ptrST_wg,3,0.0,AA,3);

        double C_IFdu = 0.0;
        double F_I_Tdu = 0.0;
        /* C_IFdu = C_I:AA*/
        /* F_I_Tdu = F_I':du*/
        for(int i = 0; i<3; i++)
        {
          for(int j=0; j<3; j++)
          {
            C_IFdu += C_I[i*3+j]*AA[i*3+j];
            F_I_Tdu += F_I[j*3+i]*ptrST_wg[i*3+j];
          }
        }

        int idx_pu = idx_K_gen(a,b,w,g,npres,1,nne,nsd);
        K[idx_pu] += -dt_alpha_1_minus_alpha*jj*wt*Jn*Np[a]*(F_I_Tdu + C_IFdu);
      }
    }
  }

  free(AA);
  free(C);
  free(C_I);
  free(F_I);
}
void add_3F_Kpt_ip(double *K,
                   int nVol, int npres, double jj, double wt, double *Nt, double *Np,
                   double dt_alpha_1_minus_alpha)
{

  for(int a=0; a<npres; a++)
  {
    for(int b=0; b<nVol; b++)
      K[idx_K_gen(a,0,b,0,npres,1,nVol,1)] += dt_alpha_1_minus_alpha*jj*wt*Np[a]*Nt[b];
  }
}
void resid_w_inertia_Ru_ip(double *fu,
                           int nne, double *ST, double *F_in, double *S_in, double jj, double wt, double Pn, double theta)
{
  int nsd = 3;
  
  Matrix<double> AA(DIM_3,DIM_3), sAA(DIM_3,DIM_3), C(DIM_3,DIM_3), S_temp(DIM_3,DIM_3);
  Matrix<double> F, S;
  F.use_reference(DIM_3,DIM_3,F_in);
  S.use_reference(DIM_3,DIM_3,S_in);  

  C.prod(F,1,F,0);
  double Jn = det3x3(F_in);
  
  Matrix<double> C_I(DIM_3,DIM_3), F_I(DIM_3,DIM_3), FIduC(DIM_3,DIM_3);  
  C_I.inv(C);
  F_I.inv(F);

  double factor = pow(theta/Jn, 2.0/3.0);
  
  for(int a=0; a<nne; a++)
  {
    for(int b=0; b<nsd; b++)
    {
      Matrix<double> Grad_du;
      Grad_du.use_reference(DIM_3,DIM_3,&ST[idx_4_gen(a,b,0,0,nne,nsd,nsd,nsd)]);
      AA.prod(F,1,Grad_du,0);

      symmetric_part(sAA.m_pdata,AA.m_pdata,DIM_3);
      
      double FIdu = 0.0;      
      for(int p=1; p<=3; p++)
      {
        for(int q=1; q<=3; q++)
        {
          FIdu += F_I(q,p)*Grad_du(p,q);
          FIduC(p,q) = C(p,q);
        }
      }

      FIduC.prod(FIdu);
      double FIduCS = FIduC.ddot(S);
      
      //sAA:S
      double sAAS = sAA.ddot(S);
      //C_IFdu = C_I:AA
      double C_IFdu = C_I.ddot(AA);

      fu[a*nsd + b] += jj*wt*(factor*(sAAS - 1.0/3.0*FIduCS) + Jn*Pn*C_IFdu);
    }
  }
}
void resid_w_inertia_Rt_ip(double *ft, int nVol, double *F_in, double *S_in, double jj, double wt, double *Nt, double theta, double Pn, double Up)
{
  Matrix<double> F, S, C(3,3);  
  F.use_reference(DIM_3,DIM_3,F_in);
  S.use_reference(DIM_3,DIM_3,S_in);  

  double Jn = det3x3(F_in);
  C.prod(F,1,F,0);
  
  double factor = pow(theta/Jn, 2.0/3.0);
  
  double SC = S.ddot(C);
  for(int a=0; a<nVol; a++)
    ft[a] += jj*wt*Nt[a]*(factor*1.0/3.0/theta*SC + Up - Pn);
}
void resid_w_inertia_Rp_ip(double *fp, int npres, double *F, double jj, double wt, double *Np, double Tn)
{
  double Jn = det3x3(F);

  for(int a=0; a<npres; a++)
    fp[a] += jj*wt*Np[a]*(Jn - Tn);
}
void TF_Ru_ip(Matrix<double> &f, FEMLIB *fe,
              Matrix<double> &F, Matrix<double> &S, double Pn, double theta)
{
  int nne = fe->nne;
  if(USE_HW_FUNCS)
  {
    double Fn[9];
    Fn[0] = Fn[4] = Fn[8] = 1.0;
    Fn[1] = Fn[2] = Fn[3] = Fn[5] = Fn[6] = Fn[7] = 0.0;
    double F_I[9];
    inverse(F.m_pdata,3,F_I);
    HW_Ru_at_ip(f.m_pdata,nne,nne,fe->ST,Fn,F.m_pdata,F_I,
                1.0,0.0,1.0,1.0,S.m_pdata,Pn,fe->detJ,fe->temp_v.w_ip,0);
  }
  else
    resid_w_inertia_Ru_ip(f.m_pdata,nne,fe->ST,F.m_pdata,S.m_pdata,fe->detJ,fe->temp_v.w_ip, Pn, theta);
}
void TF_Rt_ip(Matrix<double> &f, FEMLIB *fe,
              Matrix<double> &F, 
              Matrix<double> &S, int nVol, Matrix<double> &Nt, double theta, double Pn, double Up)
{
  if(USE_HW_FUNCS)
    HW_Rt_at_ip(f.m_pdata,nVol,Nt.m_pdata,1.0,1.0,Pn,1.0,Up,fe->detJ,fe->temp_v.w_ip);
  else
    resid_w_inertia_Rt_ip(f.m_pdata,nVol,F.m_pdata, S.m_pdata, fe->detJ,fe->temp_v.w_ip, Nt.m_pdata, theta, Pn, Up);
}
void TF_Rp_ip(Matrix<double> &f, FEMLIB *fe,
              Matrix<double> &F, int npres, Matrix<double> &Np, double Tn)
{
  if(USE_HW_FUNCS)
    HW_Rp_at_ip(f.m_pdata,npres,Np.m_pdata,1.0,1.0,1.0,1.0,fe->detJ,fe->temp_v.w_ip);
  else
    resid_w_inertia_Rp_ip(f.m_pdata, npres, F.m_pdata, fe->detJ,fe->temp_v.w_ip, Np.m_pdata, Tn);
}


void TF_Kup_ip(Matrix<double> &K, FEMLIB *fe,
               Matrix<double> &F, int npres, Matrix<double> &Np, double dt_alpha_1_minus_alpha)
{
  int nne = fe->nne;
  if(USE_HW_FUNCS)
    HW_Kup_at_ip(K.m_pdata,nne,nne,npres,Np.m_pdata,fe->ST,F.m_pdata,1.0,1.0,1.0,fe->detJ,fe->temp_v.w_ip,0);
  else
    add_3F_Kup_ip(K.m_pdata,nne,npres,fe->ST,F.m_pdata,fe->detJ,fe->temp_v.w_ip,Np.m_pdata,dt_alpha_1_minus_alpha);
}
void TF_Kut_ip(Matrix<double> &K, FEMLIB *fe,
               Matrix<double> &F, int nVol, Matrix<double> &Nt, double dt_alpha_1_minus_alpha)
{
  return;
}
void TF_Kpt_ip(Matrix<double> &K, FEMLIB *fe,
               Matrix<double> &F, int npres, int nVol, Matrix<double> &Np, Matrix<double> &Nt,
               double dt_alpha_1_minus_alpha)
{
  if(USE_HW_FUNCS)
    HW_Kpt_at_ip(K.m_pdata,npres,Np.m_pdata,nVol,Nt.m_pdata,1.0,1.0,fe->detJ,fe->temp_v.w_ip);
  else
    add_3F_Kpt_ip(K.m_pdata,nVol,npres,fe->detJ,fe->temp_v.w_ip,Nt.m_pdata,Np.m_pdata,dt_alpha_1_minus_alpha);
}
void TF_Ktt_ip(Matrix<double> &K, FEMLIB *fe,
               Matrix<double> &F, int nVol, Matrix<double> &Nt, double Upp,
               double dt_alpha_1_minus_alpha)
{
  if(USE_HW_FUNCS)
    HW_Ktt_at_ip(K.m_pdata,nVol,Nt.m_pdata,1.0,1.0,1.0,Upp,fe->detJ,fe->temp_v.w_ip);
  else
    add_3F_Ktt_ip(K.m_pdata,nVol,fe->detJ,fe->temp_v.w_ip,Nt.m_pdata,Upp,dt_alpha_1_minus_alpha);
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
    P.m_pdata[0] = alpha_1*fv->eps[eid].d_T[1] + alpha_2*fv->eps[eid].d_T[0];
    
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

    for(int a=0; a<Vno; a++)
      Tn += Nt(a+1)*(alpha_1*fv->eps[eid].T[a*3+1] + alpha_2*fv->eps[eid].T[a*3+0]);

    for(int a=1; a<=Pno; a++)
      Pn += Np(a)*P(a);
      
    var.set_variables(Tn, Pn);      
    
    Matrix<double> S,F;
    F.use_reference(DIM_3,DIM_3,var.mF);
    S.use_reference(DIM_3,DIM_3,var.mS);
            
    TF_Kup_ip(Kup,fe,F,Pno,Np,dt_alpha_1_minus_alpha);
    TF_Kpt_ip(Kpt,fe,F,Pno,Vno,Np,Nt,dt_alpha_1_minus_alpha);

    TF_Kuu_ip(Kuu, fe, var, dt_alpha_1_minus_alpha);
    TF_Kut_ip(Kut, fe, var, dt_alpha_1_minus_alpha,Vno,Nt);   
    TF_Ktt_ip(Ktt, fe, var, dt_alpha_1_minus_alpha,Vno,Nt);
  }
  
  Kpu.trans(Kup);
  Ktp.trans(Kpt);
  
  Ktu.trans(Kut);

  err += condense_K_3F_to_1F(lk, nne, nsd, Pno, Vno,
                             Kuu.m_pdata, Kut.m_pdata, Kup.m_pdata,
                             Ktu.m_pdata, Ktt.m_pdata, Ktp.m_pdata,
                             Kpu.m_pdata, Kpt.m_pdata, NULL);

if(eid==-1)
{
  printf("this is running\n");

Kuu.print("Kuu");
Kut.print("Kut");
Kup.print("Kup");
Ktu.print("Ktu");
Ktt.print("Ktt");
Ktp.print("Ktp");
Kpu.print("Kpu");
Kpt.print("Kpt");

//  Matrix<double> k;
//  k.use_reference(nne*nsd,nne*nsd, lk);
//  k.print("k");
  exit(0);
}

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
/// \param[out] lk    computed element stiffness matrix
/// \param[in]  r_e   nodal variabls(displacements) on the current element
/// \param[in]  grid  a mesh object
/// \param[in]  mat   a material object
/// \param[in]  fv    object for field variables
/// \param[in]  alpha mid point alpha
/// \param[in]  dt    time step size
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
  
  Matrix<double> fu(nne*nsd,1,    0.0);
  Matrix<double> fp(Pno,  1,    0.0);
  Matrix<double> ft(Vno,   1,    0.0);  

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
    P.m_pdata[0] = fv->eps[eid].d_T[0];
    
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

    for(int a=0; a<Vno; a++)
      Tn += Nt(a+1)*fv->eps[eid].T[a*3+0];

    for(int a=1; a<=Pno; a++)
      Pn += Np(a)*P(a);
      
    var.set_variables(Tn, Pn);      
    
    Matrix<double> S,F;
    F.use_reference(DIM_3,DIM_3,var.mF);
    S.use_reference(DIM_3,DIM_3,var.mS);
            
    TF_Kup_ip(Kup,fe,F,Pno,Np,-1.0);
    TF_Kpt_ip(Kpt,fe,F,Pno,Vno,Np,Nt,-1.0);

    TF_Kut_ip(Kut, fe, var, -1.0,Vno,Nt);   
    TF_Ktt_ip(Ktt, fe, var, -1.0,Vno,Nt);
    
    TF_Ru_ip(fu,fe,F,S,Pn,Tn);
    TF_Rp_ip(fp,fe,F,Pno,Np,Tn);
    TF_Rt_ip(ft,fe,F,S,Vno,Nt,Tn,Pn,var.Up);    
  }
  
  Ktp.trans(Kpt);
  
  err += condense_F_3F_to_1F(f, nne, nsd, Pno, Vno,
                             fu.m_pdata, ft.m_pdata, fp.m_pdata,
                             Kut.m_pdata, Kup.m_pdata, Ktp.m_pdata, Ktt.m_pdata, Kpt.m_pdata);
  
  if(eid==-1)
  {

    u.print("u");
    fu.print("Ru");
    ft.print("Rt");
    fp.print("Rp");
    Matrix<double> ff;
    ff.use_reference(nne*nsd, 1, f);
    ff.print("f");
    exit(0);
  }
    
  return err;
}
void residuals_3f_w_inertia_el(double *f,
                               const int ii,
                               const int ndofn,
                               const int nne,
                               const int npres,
                               const int nVol,
                               const int nsd,
                               const double *x,
                               const double *y,
                               const double *z,
                               const Element *elem,
                               const HOMMAT *hommat,
                               const Node *node,
                               const double *dts,
                               SIG *sig,
                               EPS *eps,
                               double alpha, double *r_n_a, double *r_n_1_a)
{
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

  double *u1 = aloc1(nne*nsd);
  double *u2 = aloc1(nne*nsd);
  double *P1 = aloc1(npres);
  double *P2 = aloc1(npres);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
    {
      u1[a*nsd+b] = r_n_1_a[a*ndofn+b];
      u2[a*nsd+b] = r_n_a[a*ndofn+b];
    }

    if(npres==nne)
    {
      P1[a] = r_n_1_a[a*ndofn+nsd];
      P2[a] = r_n_a[a*ndofn+nsd];
    }
  }

  if(npres==1)
  {
    P1[0] = alpha_1*eps[ii].d_T[2] + alpha_2*eps[ii].d_T[1];
    P2[0] = alpha_1*eps[ii].d_T[1] + alpha_2*eps[ii].d_T[0];
  }
  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.0*(1.0-2.0*hommat[mat].nu));

  dUdJFuncPtr          UP = getDUdJFunc(1, &hommat[mat]);
  d2UdJ2FuncPtr       UPP = getD2UdJ2Func(1, &hommat[mat]);
  devStressFuncPtr Stress = getDevStressFunc(1,&hommat[mat]);

  Matrix<double> F1,F2,C1,C2,S1,S2;
  F1.initialization(3,3,0.0);
  F2.initialization(3,3,0.0);
  C1.initialization(3,3,0.0);
  C2.initialization(3,3,0.0);
  S1.initialization(3,3,0.0);
  S2.initialization(3,3,0.0);

  Matrix<double> fu,fu1,fu2,Kut;
  Matrix<double> fp,fp1,fp2,Kpt;
  Matrix<double> ft,ft1,ft2,Ktt;
  Matrix<double> Kup,Ktp;

  fu.initialization(nne*nsd,1,    0.0);
  fu1.initialization(nne*nsd,1,    0.0);
  fu2.initialization(nne*nsd,1,    0.0);
  fp.initialization(npres,  1,    0.0);
  fp1.initialization(npres,  1,    0.0);
  fp2.initialization(npres,  1,    0.0);
  ft.initialization(nVol,   1,    0.0);
  ft1.initialization(nVol,   1,    0.0);
  ft2.initialization(nVol,   1,    0.0);
  Kut.initialization(nne*nsd,nVol, 0.0);
  Kpt.initialization(npres,  nVol, 0.0);
  Ktt.initialization(nVol,   nVol, 0.0);
  Kup.initialization(nne*nsd,npres,0.0);
  Ktp.initialization(nVol,   npres,0.0);

  FEMLIB fe(ii, elem, node, INTG_ORDER,1);

  Matrix<double> Np, Nt;
  Np.initialization(npres,1,0.0);
  Nt.initialization(nVol, 1,0.0);

  for(int ip = 1; ip<=fe.nint; ip++)
  {
    fe.elem_basis_V(ip);
    fe.elem_shape_function(ip,npres,Np.m_pdata);
    fe.elem_shape_function(ip,nVol, Nt.m_pdata);
    fe.update_shape_tensor();
    fe.update_deformation_gradient(ndofn,u1,F1.m_pdata);
    fe.update_deformation_gradient(ndofn,u2,F2.m_pdata);

    double Tn1 = 0.0; // 1: n-1+alpha
    double Tn2 = 0.0; // 2: n+alpha
    double Pn1 = 0.0;
    double Pn2 = 0.0;

    for(int a=0; a<nVol; a++)
    {
      Tn1 += Nt(a+1)*(alpha_1*eps[ii].T[a*3+2] + alpha_2*eps[ii].T[a*3+1]);
      Tn2 += Nt(a+1)*(alpha_1*eps[ii].T[a*3+1] + alpha_2*eps[ii].T[a*3+0]);
    }
    for(int a=0; a<npres; a++)
    {
      Pn1 += Np(a+1)*P1[a];
      Pn2 += Np(a+1)*P2[a];
    }


    double Upp2 = 0.0;
    double Up1 = 0.0;
    double Up2 = 0.0;
    UPP(Tn2,&hommat[mat],&Upp2);
    Upp2 = Upp2*kappa;

    UP(Tn1, &hommat[mat], &Up1);
    Up1 = Up1*kappa;

    UP(Tn2, &hommat[mat], &Up2);
    Up2 = Up2*kappa;

    C1.prod(F1,1,F1,0);
    Stress(C1.m_pdata,&hommat[mat],S1.m_pdata);

    C2.prod(F2,1,F2,0);
    Stress(C2.m_pdata,&hommat[mat],S2.m_pdata);

    TF_Kpt_ip(Kpt,&fe,F2,npres,nVol,Np,Nt,dt_alpha_1_minus_alpha);
    TF_Ktt_ip(Ktt,&fe,F2,nVol,Nt,Upp2,dt_alpha_1_minus_alpha);

    if(npres==nne)
      TF_Kut_ip(Kut,&fe,F2,nVol,Nt,dt_alpha_1_minus_alpha);

    if(npres==1)
      TF_Kup_ip(Kup,&fe,F2,npres,Np,dt_alpha_1_minus_alpha);


    TF_Ru_ip(fu1,&fe,F1,S1,Pn1,Tn1);
    TF_Rp_ip(fp1,&fe,F1,npres,Np,Tn1);
    TF_Rt_ip(ft1,&fe,F1,S1,nVol,Nt,Tn1,Pn1,Up1);

    TF_Ru_ip(fu2,&fe,F2,S2,Pn2,Tn2);
    TF_Rp_ip(fp2,&fe,F2,npres,Np,Tn2);
    TF_Rt_ip(ft2,&fe,F2,S2,nVol,Nt,Tn2,Pn2,Up2);
  }
  Ktp.trans(Kpt);

  for(int a=0; a<nne*nsd; a++)
    fu.m_pdata[a] = dt_alpha_1*fu2.m_pdata[a] + dt_alpha_2*fu1.m_pdata[a];

  for(int a=0; a<nVol; a++)
    ft.m_pdata[a] = dt_alpha_1*ft2.m_pdata[a] + dt_alpha_2*ft1.m_pdata[a];

  for(int a=0; a<npres; a++)
    fp.m_pdata[a] = dt_alpha_1*fp2.m_pdata[a] + dt_alpha_2*fp1.m_pdata[a];


  condense_F_out(f,nne,nsd,npres,nVol,fu.m_pdata,ft.m_pdata,fp.m_pdata,
                 Kut.m_pdata,Kup.m_pdata,Ktp.m_pdata,Ktt.m_pdata,Kpt.m_pdata);

  free(P1); free(P2); free(u1); free(u2);
}

void update_3f_state_variables_ip(int ii, int ip,
                                  const Element *elem,
                                  const HOMMAT *hommat,
                                  SIG *sig,
                                  EPS *eps,
                                  double *F,
                                  double pressure,
                                  double theta,
                                  double volume,
                                  double detJxW)
{
  if(ip==1)
  {
    memset(sig[ii].el.o,0,6*sizeof(double));
    memset(eps[ii].el.o,0,6*sizeof(double));
  }

  const int mat = elem[ii].mat[2];

  double F_total[9], C[9], C_I[9], S[9], AA[9];
  int err = 0;
  double J = getJacobian(F,ii,&err);
  double alpha = pow(theta/J,1.0/3.0);

  for(int a = 0; a<9; a++)
    F_total[a] = alpha*F[a];

  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              3,3,3,1.0,F_total,3,F_total,3,0.0,C,3);

  inverse(C,3,C_I);

  // Get Deviatoric 2 P-K stress
  devStressFuncPtr Stress;
  Stress = getDevStressFunc(1,&hommat[mat]);
  Stress(C,&hommat[mat],S);

  /* Compute total stress S = dev(S) + p*Tn*Tr*C_I */
  cblas_daxpy(9,pressure*theta,C_I,1,S,1);

  /* store S at ip */
  sig[ii].il[ip-1].o[0] = S[idx_2(0,0)]; /* XX */
  sig[ii].il[ip-1].o[1] = S[idx_2(1,1)]; /* YY */
  sig[ii].il[ip-1].o[2] = S[idx_2(2,2)]; /* ZZ */
  sig[ii].il[ip-1].o[3] = S[idx_2(1,2)]; /* YZ */
  sig[ii].il[ip-1].o[4] = S[idx_2(0,2)]; /* XZ */
  sig[ii].il[ip-1].o[5] = S[idx_2(0,1)]; /* XY */

  /* Compute Cauchy stress (theta)^-1 F S F' */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              3,3,3,1.0,F_total,3,S,3,0.0,AA,3);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
              3,3,3,1.0/theta,AA,3,F_total,3,0.0,S,3);

  sig[ii].el.o[0] += detJxW/volume*S[idx_2(0,0)];
  sig[ii].el.o[1] += detJxW/volume*S[idx_2(1,1)];
  sig[ii].el.o[2] += detJxW/volume*S[idx_2(2,2)];

  sig[ii].el.o[3] += detJxW/volume*S[idx_2(1,2)];
  sig[ii].el.o[4] += detJxW/volume*S[idx_2(0,2)];
  sig[ii].el.o[5] += detJxW/volume*S[idx_2(0,1)];

  double E[9], F_total_I[9];

  /* Compute G-L strain E = 0.5(C - 1) */
  for(int a=0; a<9; a++)
    E[a] = 0.5*(C[a] - 1.0*(a==0 || a==4 || a==8));

  /* Elastic Green Lagrange strain */
  eps[ii].il[ip-1].o[0] = E[idx_2(0,0)];
  eps[ii].il[ip-1].o[1] = E[idx_2(1,1)];
  eps[ii].il[ip-1].o[2] = E[idx_2(2,2)];
  eps[ii].il[ip-1].o[3] = E[idx_2(1,2)]*2.0;
  eps[ii].il[ip-1].o[4] = E[idx_2(0,2)]*2.0;
  eps[ii].il[ip-1].o[5] = E[idx_2(0,1)]*2.0;


  /* Compute logarithmic strain */
  inverse(F_total,3,F_total_I);
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              3,3,3,1.0,F_total_I,3,E,3,0.0,AA,3);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              3,3,3,1.0,AA,3,F_total_I,3,0.0,E,3);

  eps[ii].el.o[0] += detJxW/volume*E[idx_2(0,0)];
  eps[ii].el.o[1] += detJxW/volume*E[idx_2(1,1)];
  eps[ii].el.o[2] += detJxW/volume*E[idx_2(2,2)];
  eps[ii].el.o[3] += detJxW/volume*E[idx_2(1,2)]*2.0;
  eps[ii].el.o[4] += detJxW/volume*E[idx_2(0,2)]*2.0;
  eps[ii].el.o[5] += detJxW/volume*E[idx_2(0,1)]*2.0;
}


void update_3f_state_variables_el(const int ii,
                                  const int ndofn,
                                  const int nne,
                                  const int npres,
                                  const double *x,
                                  const double *y,
                                  const double *z,
                                  const Element *elem,
                                  const HOMMAT *hommat,
                                  const Node *node,
                                  double *u,
                                  const double *P,
                                  double dt,
                                  SIG *sig,
                                  EPS *eps)
{

  // const int mat = elem[ii].mat[2];
  // const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));
  const double volume = Tetra_V(x,y,z);

  const int nVol = N_VOL_TREE_FIELD;

  Matrix<double> F(3,3,0.0);

  FEMLIB fe(ii, elem, node, 1,1);

  Matrix<double> Np(npres,1,0.0), Nt(nVol, 1,0.0);

  for(int ip = 1; ip<=fe.nint; ip++)
  {
    fe.elem_basis_V(ip);
    fe.elem_shape_function(ip,npres,Np.m_pdata);
    fe.elem_shape_function(ip,nVol, Nt.m_pdata);
    fe.update_shape_tensor();
    fe.update_deformation_gradient(ndofn,u,F.m_pdata);

    double Tn = 0.0;
    double Pn = 0.0;

    for(int a=0; a<nVol; a++)
      Tn += Nt(a+1)*eps[ii].T[a*3+0];

    for(int a=0; a<npres; a++)
      Pn += Np(a+1)*P[a];

    update_3f_state_variables_ip(ii,ip,elem,hommat,sig,eps,F.m_pdata,Pn,Tn,volume,fe.detJxW);

  }
}


void evaluate_PT_el(const int ii,
                    const int ndofn,
                    const int nne,
                    const int npres,
                    const int nVol,
                    const int nsd,
                    const double *x,
                    const double *y,
                    const double *z,
                    const Element *elem,
                    const HOMMAT *hommat,
                    const Node *node,
                    const double *r_e,
                    double *du,
                    double dt,
                    SIG *sig,
                    EPS *eps)
{
  double *P, *u;
  u = aloc1(nne*nsd);
  P = aloc1(npres);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];

    if(npres==nne)
      P[a] = r_e[a*ndofn+nsd];
  }
  if(npres==1)
    P[0] = eps[ii].d_T[0];

  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));

  dUdJFuncPtr UP = getDUdJFunc(1, &hommat[mat]);
  d2UdJ2FuncPtr UPP = getD2UdJ2Func(1, &hommat[mat]);
  devStressFuncPtr Stress   = getDevStressFunc(1,&hommat[mat]);

  Matrix<double> F(3,3,0.0);

  Matrix<double> fp(npres, 1, 0.0), ft(nVol,  1, 0.0);

  Matrix<double> Ktu,Ktt,Ktp;
  Matrix<double> Kup,Kpu,Kpt;

  Ktu.initialization(nVol   ,nne*nsd,0.0);
  Ktt.initialization(nVol   ,nVol   ,0.0);
  Kup.initialization(nne*nsd,npres  ,0.0);
  Kpu.initialization(npres  ,nne*nsd,0.0);
  Kpt.initialization(npres  ,nVol   ,0.0);
  Ktp.initialization(nVol   ,npres  ,0.0);

  FEMLIB fe(ii, elem, node, INTG_ORDER,1);

  Matrix<double> Np(npres,1,0.0), Nt(nVol, 1,0.0);

  for(int ip = 1; ip<=fe.nint; ip++)
  {
    fe.elem_basis_V(ip);
    fe.elem_shape_function(ip,npres,Np.m_pdata);
    fe.elem_shape_function(ip,nVol, Nt.m_pdata);
    fe.update_shape_tensor();
    fe.update_deformation_gradient(ndofn,u,F.m_pdata);

    double Tn = 0.0;
    double Pn = 0.0;

    for(int a=0; a<nVol; a++)
      Tn += Nt(a+1)*eps[ii].T[a*3+0];

    for(int a=0; a<npres; a++)
      Pn += Np(a+1)*P[a];


    double Upp = 0.0;
    double Up = 0.0;
    UPP(Tn,&hommat[mat],&Upp);
    Upp = Upp*kappa;

    UP(Tn, &hommat[mat], &Up);
    Up = Up*kappa;

    TF_Kup_ip(Kup,&fe,F,npres,Np,-1.0);

    TF_Kpt_ip(Kpt,&fe,F,npres,nVol,Np,Nt,-1.0);
    TF_Ktt_ip(Ktt,&fe,F,nVol,Nt,Upp,-1.0);
    
    Matrix<double> S(DIM_3,DIM_3), C(DIM_3,DIM_3);
    C.prod(F,1,F,0);
    Stress(C.m_pdata,&hommat[mat],S.m_pdata);
    TF_Rt_ip(ft,&fe,F,S,nVol,Nt,Tn,Pn,Up);
    TF_Rp_ip(fp,&fe,F,npres,Np,Tn);
  }

  free(u);
  free(P);

  Kpu.trans(Kup);
  Ktp.trans(Kpt);

  Matrix<double> theta(nVol,1,0.0),KptI(nVol,nVol,0.0);
  KptI.inv(Kpt);

  Matrix<double> uu;
  uu.use_reference(nne*ndofn, 1, du);


  Matrix<double> Kpu_uu;
  Kpu_uu.prod(Kpu,uu);
  fp.add(Kpu_uu);
    

  theta.prod(KptI,fp);
  theta.prod(-1.0); 

  for(int a = 0; a<nVol; a++)
    eps[ii].T[a*3+0] += theta.m_pdata[a];
    
    

  Matrix<double> press(npres,1,0.0),KtpI(nVol,npres,0.0);
  KtpI.inv(Ktp);

  Matrix<double> Ktt_theta;
  Ktt_theta.prod(Ktt,theta);
  ft.add(Ktt_theta);

  press.prod(KtpI,ft);
  press.prod(-1.0);

  for(int a = 0; a<npres; a++)
    eps[ii].d_T[a*3+0] += press.m_pdata[a];
    
  if(ii==-1)
        press.print("B");     
}

void evaluate_PT_el_test(FEMLIB *fe,
                         Grid *grid,
                         MaterialProperty *material,
                         FieldVariables *fv,
                         Matrix<double> &r_e,
                         Matrix<double> &du)
{
  int eid = fe->curt_elem_id;
  Var_Carrier var;
  var.set_elasticity_functions(material,grid,eid);
  EPS *eps = fv->eps;
  
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
    P.m_pdata[0] = eps[eid].d_T[0];

  Matrix<double> fp(Pno, 1, 0.0), ft(Vno,  1, 0.0);

  Matrix<double> Ktu(nne*nsd,Vno, 0.0); // will be transposed: Kut -> Ktu
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

    for(int a=0; a<Vno; a++)
      Tn += Nt(a+1)*eps[eid].T[a*3+0];

    for(int ia=1; ia<=Pno; ia++)
      Pn += Np(ia)*P(ia);

    var.set_variables(Tn, Pn);
    Matrix<double> S,F;
    F.use_reference(DIM_3,DIM_3,var.mF);
    S.use_reference(DIM_3,DIM_3,var.mS);    

    TF_Kup_ip(Kpu,fe,F,Pno,Np,-1.0);
    TF_Kpt_ip(Kpt,fe,F,Pno,Vno,Np,Nt,-1.0);
    TF_Ktt_ip(Ktt,fe,F,Vno,Nt,var.Upp,-1.0);
    TF_Kut_ip(Ktu, fe, var, -1.0,Vno,Nt);    
    TF_Rt_ip(ft,fe,F,S,Vno,Nt,Tn,Pn,var.Up);
    TF_Rp_ip(fp,fe,F,Pno,Np,Tn);
  }

  Kpu.trans();
  Ktp.trans(Kpt);
  Ktu.trans();

  Matrix<double> d_theta(Vno, 1, 0.0), dP(Pno, 1, 0.0);
  
  compute_d_theta_dP_test(d_theta,dP, 
                            du,ft,fp,Kpu,Ktu,Ktp,Ktt,Kpt);

  for(int ia = 0; ia<Pno; ia++)
    eps[eid].d_T[ia*3+0] += dP.m_pdata[ia];

  for(int ia=0; ia<Vno; ia++)
    eps[eid].T[ia*3+0] += d_theta.m_pdata[ia];
}

void evaluate_PT_w_inertia_el(const int ii,
                              const int ndofn,
                              const int nne,
                              const int npres,
                              const int nVol,
                              const int nsd,
                              const double *x,
                              const double *y,
                              const double *z,
                              const Element *elem,
                              const HOMMAT *hommat,
                              const Node *node,
                              double *du,
                              double dt,
                              SIG *sig,
                              EPS *eps,
                              double alpha, double *r_n_a, double *r_n_1_a)
{
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
    dt_alpha_1_minus_alpha = dt*alpha_1*alpha_2;
    dt_alpha_1 = -dt*alpha_1;
    dt_alpha_2 = -dt*alpha_2;
  }

  double *u2 = aloc1(nne*nsd); // 2: n+alpha
  double *u1 = aloc1(nne*nsd); // 1: n-1+alpha
  double *P2 = aloc1(npres);   // 1: n+alpha
  double *P1 = aloc1(npres);   // 1: n-1+alpha

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
    {
      u1[a*nsd+b] = r_n_1_a[a*ndofn+b];
      u2[a*nsd+b] = r_n_a[a*ndofn+b];
    }
    if(npres==nne)
    {
      P1[a] = r_n_1_a[a*ndofn+nsd];
      P2[a] = r_n_a[a*ndofn+nsd];
    }
  }
  if(npres==1)
  {
    P1[0] = alpha_1*eps[ii].d_T[2] + alpha_2*eps[ii].d_T[1];
    P2[0] = alpha_1*eps[ii].d_T[1] + alpha_2*eps[ii].d_T[0];
  }

  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.0*(1.0-2.0*hommat[mat].nu));

  dUdJFuncPtr          UP = getDUdJFunc(1, &hommat[mat]);
  d2UdJ2FuncPtr       UPP = getD2UdJ2Func(1, &hommat[mat]);
  devStressFuncPtr Stress   = getDevStressFunc(1,&hommat[mat]);

  Matrix<double> F1(3,3,0.0),F2(3,3,0.0);

  Matrix<double> fp, fp1, fp2;
  Matrix<double> ft, ft1, ft2;

  fp.initialization(npres,  1,    0.0);
  fp1.initialization(npres,  1,    0.0);
  fp2.initialization(npres,  1,    0.0);
  ft.initialization(nVol,   1,    0.0);
  ft1.initialization(nVol,   1,    0.0);
  ft2.initialization(nVol,   1,    0.0);

  Matrix<double> Ktu,Ktt,Ktp;
  Matrix<double> Kup,Kpu,Kpt;

  Ktu.initialization(nVol   ,nne*nsd,0.0);
  Ktt.initialization(nVol   ,nVol   ,0.0);
  Kup.initialization(nne*nsd,npres  ,0.0);
  Kpu.initialization(npres  ,nne*nsd,0.0);
  Kpt.initialization(npres  ,nVol   ,0.0);
  Ktp.initialization(nVol   ,npres  ,0.0);

  FEMLIB fe(ii, elem, node, INTG_ORDER,1);

  Matrix<double> Np(npres,1,0.0),Nt(nVol, 1,0.0);

  for(int ip = 1; ip<=fe.nint; ip++)
  {
    fe.elem_basis_V(ip);
    fe.elem_shape_function(ip,npres,Np.m_pdata);
    fe.elem_shape_function(ip,nVol, Nt.m_pdata);
    fe.update_shape_tensor();
    fe.update_deformation_gradient(ndofn,u1,F1.m_pdata);
    fe.update_deformation_gradient(ndofn,u2,F2.m_pdata);

    double Tn1 = 0.0; // 1: n-1+alpha
    double Tn2 = 0.0; // 2: n+alpha
    double Pn1 = 0.0;
    double Pn2 = 0.0;

    for(int a=0; a<nVol; a++)
    {
      Tn1 += Nt(a+1)*(alpha_1*eps[ii].T[a*3+2] + alpha_2*eps[ii].T[a*3+1]);
      Tn2 += Nt(a+1)*(alpha_1*eps[ii].T[a*3+1] + alpha_2*eps[ii].T[a*3+0]);
    }
    for(int a=0; a<npres; a++)
    {
      Pn1 += Np(a+1)*P1[a];
      Pn2 += Np(a+1)*P2[a];
    }

    double Upp2 = 0.0;
    double Up1 = 0.0;
    double Up2 = 0.0;
    UPP(Tn2,&hommat[mat],&Upp2);
    Upp2 = Upp2*kappa;

    UP(Tn1, &hommat[mat], &Up1);
    Up1 = Up1*kappa;

    UP(Tn2, &hommat[mat], &Up2);
    Up2 = Up2*kappa;

    TF_Kup_ip(Kup,&fe,F2,npres,Np,dt_alpha_1_minus_alpha);

    TF_Kpt_ip(Kpt,&fe,F2,npres,nVol,Np,Nt,dt_alpha_1_minus_alpha);
    TF_Ktt_ip(Ktt,&fe,F2,nVol,Nt,Upp2,dt_alpha_1_minus_alpha);

    Matrix<double> S1(DIM_3,DIM_3), C1(DIM_3,DIM_3);
    Matrix<double> S2(DIM_3,DIM_3), C2(DIM_3,DIM_3);    
    C1.prod(F1,1,F1,0);
    C2.prod(F2,1,F2,0);
    Stress(C1.m_pdata,&hommat[mat],S2.m_pdata);
    Stress(C2.m_pdata,&hommat[mat],S2.m_pdata);    
    
    TF_Rt_ip(ft1,&fe,F1,S1,nVol,Nt,Tn1,Pn1,Up1);
    TF_Rt_ip(ft2,&fe,F2,S2,nVol,Nt,Tn2,Pn2,Up2);

    TF_Rp_ip(fp1,&fe,F1,npres,Np,Tn1);
    TF_Rp_ip(fp2,&fe,F2,npres,Np,Tn2);
  }

  Kpu.trans(Kup);
  Ktp.trans(Kpt);

  for(int a=0; a<nVol; a++)
    ft.m_pdata[a] = dt_alpha_1*ft2.m_pdata[a] + dt_alpha_2*ft1.m_pdata[a];

  for(int a=0; a<npres; a++)
    fp.m_pdata[a] = dt_alpha_1*fp2.m_pdata[a] + dt_alpha_2*fp1.m_pdata[a];

  free(P1); free(P2); free(u1); free(u2);

  Matrix<double> theta(nVol,1,0.0),KptI(nVol,npres,0.0);
  KptI.inv(Kpt);


  Matrix<double> uu;
  uu.use_reference(nne*ndofn, 1, du);

  Matrix<double> Kpu_uu;
  Kpu_uu.prod(Kpu,uu);
  fp.add(Kpu_uu);

  theta.prod(KptI,fp);
  theta.prod(-1.0);

  for(int a = 0; a<nVol; a++)
    eps[ii].T[a*3+0] += theta.m_pdata[a];

  Matrix<double> press(npres,1,0.0),KtpI(nVol,npres,0.0);
  KtpI.inv(Ktp);

  Matrix<double> Ktt_theta;
  Ktt_theta.prod(Ktt,theta);
  ft.add(Ktt_theta);

  press.prod(KtpI,ft);
  press.prod(-1.0);

  for(int a = 0; a<npres; a++)
    eps[ii].d_T[a*3+0] += press.m_pdata[a];
}

void evaluate_theta_el(const int ii,
                       const int ndofn,
                       const int nne,
                       const int npres,
                       const int nVol,
                       const int nsd,
                       const double *x,
                       const double *y,
                       const double *z,
                       const Element *elem,
                       const HOMMAT *hommat,
                       const Node *node,
                       const double *r_e,
                       const double *du,
                       double dt,
                       SIG *sig,
                       EPS *eps)
{
  double *P, *u;
  u = aloc1(nne*nsd);
  P = aloc1(npres);

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
      u[a*nsd+b] = r_e[a*ndofn+b];

    if(npres==nne)
      P[a] = r_e[a*ndofn+nsd];
  }
  if(npres==1)
    P[0] = eps[ii].d_T[0];

  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));

  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UP = getDUdJFunc(1, &hommat[mat]);
  UPP = getD2UdJ2Func(1, &hommat[mat]);

  double *ft  = aloc1(nVol);
  double *Ktu = aloc1(nne*nsd*nVol);
  double *Ktp = aloc1(npres*nVol);
  double *Ktt = aloc1(nVol*nVol);

  memset(ft, 0,        nVol*sizeof(double));
  memset(Ktu,0,nne*nsd*nVol*sizeof(double));
  memset(Ktp,0,  nVol*npres*sizeof(double));
  memset(Ktt,0,     nVol*nVol*sizeof(double));

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(nne,&npt_z);

  double *int_pt_ksi, *int_pt_eta, *int_pt_zet, *weights;
  int_pt_ksi = aloc1(npt_z);
  int_pt_eta = aloc1(npt_z);
  int_pt_zet = aloc1(npt_z);
  weights = aloc1(npt_z);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST;
  double Upp, Up;

  Na = aloc1(nne);
  Np = aloc1(npres);
  Nt = aloc1(nVol);
  N_x = aloc1(nne);
  N_y = aloc1(nne);
  N_z = aloc1(nne);
  ST_tensor = aloc4(3,3,nsd,nne);
  ST = aloc1(3*3*nsd*nne);

  /*=== INTEGRATION LOOP ===*/
  integrate(nne,&npt_x,&npt_y,&npt_z,
            int_pt_ksi,int_pt_eta,int_pt_zet,
            weights);

  double **F_mat,*F;

  F = aloc1(9);
  F_mat = aloc2(3,3);

  for(long i=0; i<npt_x; i++)
  {
    for(long j=0; j<npt_y; j++)
    {
      for(long k=0; k<npt_z; k++)
      {
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nne,Na);
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],npres,Np);
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nVol,Nt);

        double detJ = deriv(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nne,x,y,z,N_x,N_y,N_z);
        double wt = weights[k];

        double Tn = 0.0;
        double Pn = 0.0;

        for(int a=0; a<nVol; a++)
          Tn += Nt[a]*eps[ii].T[a*3+0];

        for(int a=0; a<npres; a++)
          Pn += Np[a]*P[a];


        UPP(Tn,&hommat[mat],&Upp);
        Upp = Upp*kappa;

        UP(Tn, &hommat[mat], &Up);
        Up = Up*kappa;

        shape_tensor(nne,ndn,N_x,N_y,N_z,ST_tensor);
        shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);

        def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u,F_mat);
        mat2array(F,CONST_2(double) F_mat,3,3);

        //        add_3F_Ktu_ip(Ktu,nne,nVol,  ST,F,detJ,wt,Nt,-1.0);
        add_3F_Ktp_ip(Ktp,nVol,npres,detJ,wt,Nt,Np,-1.0);
        add_3F_Ktt_ip(Ktt,nVol,detJ,wt,Nt,Upp,-1.0);

        Matrix<double> S(3,3,0.0);
        resid_w_inertia_Rt_ip(ft, nVol, F, S.m_pdata, detJ, wt, Nt, Tn, Pn, Up);

      }
    }
  }

  free(u);
  free(P);

  dealoc4(ST_tensor,3,3,nsd);
  free(ST);

  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(Np);
  free(Nt);
  free(N_x);
  free(N_y);
  free(N_z);

  free(F);
  dealoc2(F_mat,3);


  double *KttI, *KttIKtu, *KttIKtp;
  KttI = aloc1(nVol*nVol);
  KttIKtu = aloc1(nne*nsd*nVol);
  KttIKtp = aloc1(npres*nVol);

  inverse(Ktt,nVol,KttI);

  double *theta1 = aloc1(nVol);
  double *theta2 = aloc1(nVol);
  double *theta3 = aloc1(nVol);


  // KttI:    [nVol]*[nVol]
  // Ktu:     [nVol]*[nne*nsd]
  // KttIKtu: [nVol]*[nne*nsd]
  // Ktp:     [nVol]*[npres]
  // KttIKpt: [nVol]*[npres]

  // A[m,k] x B[k,n] = C[m,n]
  // cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,A,k,B,n,0.0,C,n);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,nne*nsd,nVol,1.0,KttI,nVol,Ktu,nne*nsd,0.0,KttIKtu,nne*nsd);

  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,npres,nVol,1.0,KttI,nVol,Ktp,npres,0.0,KttIKtp,npres);

  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nne*nsd,1.0,KttIKtu,nne*nsd,u,1,0.0,theta1,1);

  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,npres,1.0,KttIKtp,npres,P,1,0.0,theta2,1);

  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nVol,1.0,KttI,nVol,ft,1,0.0,theta3,1);

  for(long a = 0; a<nVol; a++)
    eps[ii].T[a*3+0] = eps[ii].T[a*3+0] - (theta1[a] + theta2[a] + theta3[a]);

  free(ft);

  free(Ktu);
  free(Ktp);
  free(Ktt); free(KttI);
  free(KttIKtu);
  free(KttIKtp);
  free(theta1); free(theta2); free(theta3);
}

void evaluate_theta_w_inertia_el(const int ii,
                                 const int ndofn,
                                 const int nne,
                                 const int npres,
                                 const int nVol,
                                 const int nsd,
                                 const double *x,
                                 const double *y,
                                 const double *z,
                                 const Element *elem,
                                 const HOMMAT *hommat,
                                 const Node *node,
                                 const double *du,
                                 double dt,
                                 SIG *sig,
                                 EPS *eps,
                                 double alpha, double *r_n_a, double *r_n_1_a)
{
  double *u2 = aloc1(nne*nsd); // 2: n+alpha
  double *u1 = aloc1(nne*nsd); // 1: n-1+alpha
  double *P2 = aloc1(npres);   // 1: n+alpha
  double *P1 = aloc1(npres);   // 1: n-1+alpha

  for(int a=0;a<nne;a++)
  {
    for(int b=0; b<nsd;b++)
    {
      u1[a*nsd+b] = r_n_1_a[a*ndofn+b];
      u2[a*nsd+b] = r_n_a[a*ndofn+b];
    }
    if(npres==nne)
    {
      P1[a] = r_n_1_a[a*ndofn+nsd];
      P2[a] = r_n_a[a*ndofn+nsd];
    }
  }
  if(npres==1)
  {
    P1[0] = (1.0-alpha)*eps[ii].d_T[2] + alpha*eps[ii].d_T[1];
    P2[0] = (1.0-alpha)*eps[ii].d_T[1] + alpha*eps[ii].d_T[0];
  }

  const int mat = elem[ii].mat[2];
  const double kappa = hommat[mat].E/(3.*(1.-2.*hommat[mat].nu));

  dUdJFuncPtr UP;
  d2UdJ2FuncPtr UPP;
  UP = getDUdJFunc(1, &hommat[mat]);
  UPP = getD2UdJ2Func(1, &hommat[mat]);

  double *ft  = aloc1(nVol);
  double *ft1 = aloc1(nVol);
  double *ft2 = aloc1(nVol);
  double *Ktu = aloc1(nne*nsd*nVol);
  double *Ktp = aloc1(npres*nVol);
  double *Ktt = aloc1(nVol*nVol);

  memset(ft, 0,        nVol*sizeof(double));
  memset(ft1,0,        nVol*sizeof(double));
  memset(ft2,0,        nVol*sizeof(double));
  memset(Ktu,0,nne*nsd*nVol*sizeof(double));
  memset(Ktp,0,  nVol*npres*sizeof(double));
  memset(Ktt,0,     nVol*nVol*sizeof(double));

  /* INTEGRATION */
  long npt_x, npt_y, npt_z;
  int_point(nne,&npt_z);

  double *int_pt_ksi, *int_pt_eta, *int_pt_zet, *weights;
  int_pt_ksi = aloc1(npt_z);
  int_pt_eta = aloc1(npt_z);
  int_pt_zet = aloc1(npt_z);
  weights = aloc1(npt_z);

  /* allocate space for the shape functions, derivatives etc */
  double *Na, *Np, *Nt, *N_x, *N_y, *N_z, ****ST_tensor, *ST;
  double Upp2, Up1, Up2;

  Na = aloc1(nne);
  Np = aloc1(npres);
  Nt = aloc1(nVol);
  N_x = aloc1(nne);
  N_y = aloc1(nne);
  N_z = aloc1(nne);
  ST_tensor = aloc4(3,3,nsd,nne);
  ST = aloc1(3*3*nsd*nne);

  /*=== INTEGRATION LOOP ===*/
  integrate(nne,&npt_x,&npt_y,&npt_z,
            int_pt_ksi,int_pt_eta,int_pt_zet,
            weights);

  double **F_mat1, **F_mat2, *F1, *F2;

  F1 = aloc1(9);
  F2 = aloc1(9);
  F_mat1 = aloc2(3,3);
  F_mat2 = aloc2(3,3);

  for(long i=0; i<npt_x; i++)
  {
    for(long j=0; j<npt_y; j++)
    {
      for(long k=0; k<npt_z; k++)
      {
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nne,Na);
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],npres,Np);
        shape_func(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nVol,Nt);

        double detJ = deriv(int_pt_ksi[k],int_pt_eta[k],int_pt_zet[k],nne,x,y,z,N_x,N_y,N_z);
        double wt = weights[k];

        double Tn1 = 0.0; // 1: n-1+alpha
        double Tn2 = 0.0; // 2: n+alpha
        double Pn1 = 0.0;
        double Pn2 = 0.0;

        for(int a=0; a<nVol; a++)
        {
          Tn1 += Nt[a]*((alpha-1.0)*eps[ii].T[a*3+2] + alpha*eps[ii].T[a*3+1]);
          Tn2 += Nt[a]*((alpha-1.0)*eps[ii].T[a*3+1] + alpha*eps[ii].T[a*3+0]);
        }

        for(int a=0; a<npres; a++)
        {
          Pn1 += Np[a]*P1[a];
          Pn2 += Np[a]*P2[a];
        }

        UPP(Tn2,&hommat[mat],&Upp2);
        Upp2 = Upp2*kappa;

        UP(Tn1, &hommat[mat], &Up1);
        UP(Tn2, &hommat[mat], &Up2);
        Up1 = Up1*kappa;
        Up2 = Up2*kappa;

        shape_tensor(nne,ndn,N_x,N_y,N_z,ST_tensor);
        shapeTensor2array(ST,CONST_4(double) ST_tensor,nne);

        def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u1,F_mat1);
        mat2array(F1,CONST_2(double) F_mat1,3,3);

        def_grad_get(nne,ndofn,CONST_4(double) ST_tensor,u2,F_mat2);
        mat2array(F2,CONST_2(double) F_mat2,3,3);

        //        add_3F_Ktu_ip(Ktu,nne,nVol,  ST,F2,detJ,wt,Nt,dt*(1.0-alpha)*alpha);
        add_3F_Ktp_ip(Ktp,nVol,npres,detJ,wt,Nt,Np,dt*(1.0-alpha)*alpha);
        add_3F_Ktt_ip(Ktt,nVol,detJ,wt,Nt,Upp2,dt*(1.0-alpha)*alpha);

        Matrix<double> S1(3,3,0.0);
        Matrix<double> S2(3,3,0.0);
        resid_w_inertia_Rt_ip(ft1, nVol, F1, S1.m_pdata, detJ, wt, Nt, Tn1, Pn1, Up1);
        resid_w_inertia_Rt_ip(ft2, nVol, F2, S2.m_pdata, detJ, wt, Nt, Tn2, Pn2, Up2);

      }
    }
  }

  for(int a=0; a<nVol; a++)
    ft[a] = -(1.0 - alpha)*dt*ft2[a] - alpha*dt*ft1[a];

  free(P1); free(P2); free(u1); free(u2);
  dealoc4(ST_tensor,3,3,nsd);
  free(ST);

  free(int_pt_ksi);
  free(int_pt_eta);
  free(int_pt_zet);
  free(weights);
  free(Na);
  free(Np);
  free(Nt);
  free(N_x);
  free(N_y);
  free(N_z);

  free(F1);
  dealoc2(F_mat1,3);
  free(F2);
  dealoc2(F_mat2,3);


  double *KttI, *KttIKtu, *KttIKtp;
  KttI = aloc1(nVol*nVol);
  KttIKtu = aloc1(nne*nsd*nVol);
  KttIKtp = aloc1(npres*nVol);

  inverse(Ktt,nVol,KttI);

  double *theta1 = aloc1(nVol);
  double *theta2 = aloc1(nVol);
  double *theta3 = aloc1(nVol);


  // KttI:    [nVol]*[nVol]
  // Ktu:     [nVol]*[nne*nsd]
  // KttIKtu: [nVol]*[nne*nsd]
  // Ktp:     [nVol]*[npres]
  // KttIKpt: [nVol]*[npres]

  // A[m,k] x B[k,n] = C[m,n]
  // cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,A,k,B,n,0.0,C,n);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,nne*nsd,nVol,1.0,KttI,nVol,Ktu,nne*nsd,0.0,KttIKtu,nne*nsd);

  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,npres,nVol,1.0,KttI,nVol,Ktp,npres,0.0,KttIKtp,npres);

  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nne*nsd,1.0,KttIKtu,nne*nsd,u2,1,0.0,theta1,1);

  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,npres,1.0,KttIKtp,npres,P2,1,0.0,theta2,1);

  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nVol,1,nVol,1.0,KttI,nVol,ft,1,0.0,theta3,1);

  for(long a = 0; a<nVol; a++)
    eps[ii].T[a*3+0] = eps[ii].T[a*3+0] - (theta1[a] + theta2[a] + theta3[a]);

  free(ft);
  free(ft1); free(ft2);

  free(Ktu);
  free(Ktp);
  free(Ktt); free(KttI);
  free(KttIKtu);
  free(KttIKtp);
  free(theta1); free(theta2); free(theta3);
}

void update_3f(long ne, long ndofn, long npres, double *d_r, double *r, double *rr,
               Node *node, Element *elem, HOMMAT *hommat, SUPP sup, EPS *eps, SIG *sig, double dt, double t,
               MPI_Comm mpi_comm, const PGFem3D_opt *opts, double alpha, double *r_n, double *r_n_1,
               const int mp_id)
{
  const int mat = elem[0].mat[2];
  double rho = hommat[mat].density;
  long include_inertia = 1;

  const int nsd = 3;
  const int nVol = N_VOL_TREE_FIELD;

  if(fabs(rho)<MIN_DENSITY)
    include_inertia = 0;

  //////////////////////////////////////////////////////////////////////
  //  include_inertia = 0;
  //////////////////////////////////////////////////////////////////////

  for (int i=0;i<ne;i++)
  {

    int nne = elem[i].toe;
    long *nod = aloc1l (nne);
    elemnodes (i,nne,nod,elem);
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);

    double *r_e = aloc1 (ndofe);
    double *dr_e = aloc1 (ndofe);

    double *x,*y,*z;
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);

    long *cn = aloc1l (ndofe);

    nodecoord_total(nne,nod,node,x,y,z);

    /* code numbers on element */
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    /* deformation on element */
    def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);
    def_elem(cn,ndofe,rr,elem,node,dr_e,sup,2);

    double*du;
    du = aloc1(nne*nsd);

    for(int a=0;a<nne;a++)
    {
      for(int b=0; b<nsd;b++)
        du[a*nsd+b] = dr_e[a*ndofn+b];
    }

    if(include_inertia && (0<alpha && alpha<1))
    {
      double *r0, *r0_;
      r0      = aloc1(ndofe);
      r0_     = aloc1(ndofe);

      for (long I=0;I<nne;I++)
      {
        for(long J=0; J<nsd; J++)
        {
          r0[I*ndofn + J] =   r_n[nod[I]*ndofn + J];
          r0_[I*ndofn + J] = r_n_1[nod[I]*ndofn + J];
        }
      }

      double *r_n_a, *r_n_1_a;
      r_n_a = aloc1(ndofe);
      r_n_1_a = aloc1(ndofe);

      mid_point_rule(r_n_1_a,r0_,r0,  alpha, ndofe);
      mid_point_rule(r_n_a,  r0, r_e, alpha, ndofe);

      if(npres==1)
        evaluate_PT_w_inertia_el(i,ndofn,nne,npres,nVol,nsd,x,y,z,elem,hommat,node,du,dt,sig,eps,alpha,r_n_a,r_n_1_a);
      else
        evaluate_theta_w_inertia_el(i,ndofn,nne,npres,nVol,nsd,x,y,z,elem,hommat,node,du,dt,sig,eps,alpha,r_n_a,r_n_1_a);

      free(r0);
      free(r0_);
      free(r_n_a);
      free(r_n_1_a);
    }
    else
    {
      if(npres==1)
        evaluate_PT_el(i,ndofn,nne,npres,nVol,nsd,x,y,z,elem,hommat,node,r_e,du,dt,sig,eps);
      else
        evaluate_theta_el(i,ndofn,nne,npres,nVol,nsd,x,y,z,elem,hommat,node,r_e,du,dt,sig,eps);
    }

    free(du);
    free(nod);
    free(cn);
    free(x);
    free(y);
    free(z);
    free(r_e);
    free(dr_e);
  }
}

int compute_d_theta_dP_test(Matrix<double> &d_theta,
                            Matrix<double> &dP, 
                            Matrix<double> &du,
                            Matrix<double> &ft, 
                            Matrix<double> &fp, 
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
//  KptIFp.prod(KptI,fp);
//  KptIKpu_du.prod(KptI, Kpu, du);
//                
//  for(int ia=1; ia<=d_theta.m_row; ia++)
//    d_theta(ia) = -KptIFp(ia) - KptIKpu_du(ia);   

  Matrix<double> Kpu_du;
  Kpu_du.prod(Kpu,du);
  Kpu_du.add(fp);

  d_theta.prod(KptI,Kpu_du);
  d_theta.prod(-1.0);   

        
	Matrix<double> KtpI, KtpIFt, KtpIKtu_du;
	Matrix<double> KtpIKtt_d_theta;
  KtpI.inv(Ktp);
  
  KtpIFt.prod(KtpI, ft);
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
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] dt time step size
/// \param[in] alpha mid point rule alpha
/// \return non-zero on internal error
int update_3f(Grid *grid,
              MaterialProperty *mat,
              FieldVariables *fv,
              LoadingSteps *load,
              const PGFem3D_opt *opts,
              Multiphysics *mp,
              int mp_id,
              const double dt,
              double alpha)
{

  int total_Lagrangian = 1;
  
  Node *node = grid->node;
  Element *elem = grid->element;

  int ndofn = fv->ndofn;
  SUPP sup = load->sups[mp_id];
      

  //////////////////////////////////////////////////////////////////////
  //  include_inertia = 0;
  //////////////////////////////////////////////////////////////////////

  for (int eid=0;eid<grid->ne;eid++)
  {
    FEMLIB fe(eid,elem,node,1,total_Lagrangian);
    int nne   = fe.nne;
    int nsd   = fe.nsd;
    int ndofe = nne*ndofn;

    Matrix<long> cn(ndofe,1);
    long *nod = fe.node_id.m_pdata;

    Matrix<double> dr_e(ndofe,1), r_e(ndofe, 1), du(nne*nsd,1), u(nne*nsd, 1);
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn.m_pdata,mp_id);
    def_elem_total(cn.m_pdata,ndofe,fv->u_np1,fv->d_u,elem,node,sup,r_e.m_pdata);
    
    def_elem(cn.m_pdata,ndofe,fv->dd_u,elem,node,dr_e.m_pdata,sup,2);

    for(int a=0;a<nne;a++)
    {
      for(int b=0; b<nsd;b++)
        du.m_pdata[a*nsd+b] = dr_e.m_pdata[a*ndofn+b];
    }

    evaluate_PT_el_test(&fe,grid,mat,fv,r_e,du);
  }
return 0;  
}

void update_3f_state_variables(long ne, long ndofn, long npres, double *d_r, double *r,
                               Node *node, Element *elem, HOMMAT *hommat, SUPP sup, EPS *eps, SIG *sig, double dt, double t,
                               MPI_Comm mpi_comm,const int mp_id)
{
  // @todo Commented out as dead code. @cp shoud review. LD
  // long include_inertia = 1;
  //
  // if(fabs(rho)<MIN_DENSITY)
  //   include_inertia = 0;

  for (int i=0;i<ne;i++)
  {

    int nne = elem[i].toe;
    long *nod = aloc1l (nne);
    elemnodes (i,nne,nod,elem);
    int ndofe = get_ndof_on_elem_nodes(nne,nod,node,ndofn);

    double *r_e = aloc1 (ndofe);

    double *x,*y,*z;
    x = aloc1 (nne);
    y = aloc1 (nne);
    z = aloc1 (nne);

    long *cn = aloc1l (ndofe);

    nodecoord_total(nne,nod,node,x,y,z);

    // code numbers on element
    get_dof_ids_on_elem_nodes(0,nne,ndofn,nod,node,cn,mp_id);

    // deformation on element
    def_elem_total(cn,ndofe,r,d_r,elem,node,sup,r_e);


    int nsd = 3;

    double *u, *P;
    u = aloc1(nne*nsd);
    P = aloc1(npres);

    for(int a=0;a<nne;a++)
    {
      for(int b=0; b<nsd;b++)
      {
        u[a*nsd+b] = r_e[a*ndofn+b];
      }
      if(npres==nne)
        P[a] = r_e[a*ndofn+nsd];
    }
    if(npres==1)
      P[0] = eps[i].d_T[0];

    update_3f_state_variables_el(i,ndofn,nne,npres,x,y,z,elem,hommat,node,u,P,dt,sig,eps);

    free(u);
    free(P);

    free(nod);
    free(cn);
    free(x);
    free(y);
    free(z);
    free(r_e);
  }
}
