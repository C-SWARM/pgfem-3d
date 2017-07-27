
/// Authors:
///  Sangmin Lee, [1], <slee43@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN


#include "constitutive_model_3f.h"
#include "cm_placeholder_functions.h"

#include "hommat.h"
#include "index_macros.h"
#include "utils.h"

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
}

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

#define MAX(a, b) ((a) >= (b)? (a) : (b))
  

void Var_Carrier::compute_factor(void)
{
  factor = pow(theta_r/tJr, 2.0/3.0);
};

void Var_Carrier::set_tenosrs(double *tFr,
                              double *eFn,
                              double *M,
                              double *pFnp1,
                              double *eSd,
                              double *Ld)
{
  tFr_in   = tFr;
  eFn_in   = eFn;
  M_in     = M;
  pFnp1_in = pFnp1;
  eSd_in   = eSd;
  Ld_in    = Ld;
  compute_Z();
};
    
void Var_Carrier::set_scalars(double theta_r_in,
                              double theta_n_in,
                              double tJn_in,
                              double Jn_in,
                              double P_in,
                              double dUd_theta_in,
                              double d2Ud_theta2_in)
{
  theta_r     = theta_r_in;
  theta_n     = theta_n_in;
  tJn         = tJn_in;
  Jn          = Jn_in;
  P           = P_in;
  dUd_theta   = dUd_theta_in;
  d2Ud_theta2 = d2Ud_theta2_in;
  
  TensorA<2> tFr(tFr_in), M(M_in), eFn(eFn_in);
  eJn = ttl::det(eFn);
  JM  = ttl::det(M);
  tJr = ttl::det(tFr);
  compute_factor();
};

void Var_Carrier::compute_Z(void)
{
  TensorA<2> tFr(tFr_in), M(M_in), eFn(eFn_in);    
  Tensor<2> FrFeM = tFr(i,k)*eFn(k,l)*M(l,j);
  Z = FrFeM(k,i)*FrFeM(k,j);
};

template<class T1> void Var_Carrier::compute_Phi(T1 &Phi, double *Grad_beta_in)
{
  TensorA<2> tFr(tFr_in), eFn(eFn_in); 
  Tensor<2> A, Asym,tFrI;
  TensorA<2> Grad_beta(Grad_beta_in);
  A = Grad_beta(k,i)*tFr(k,j);
  symm(A, Asym);
  
  inv(tFr,tFrI);
  double Gb_tFrI = Grad_beta(j,i)*tFrI(i,j);
   
  Phi = eFn(k,i)*(Asym(k,l) - 1.0/3.0*Gb_tFrI*tFr(o,k)*tFr(o,l))*eFn(l,j);
};

template<class T1, class T2> void Var_Carrier::compute_Psi(T1 &Psi, T2 &Phi)    
{
  TensorA<2> M(M_in); 
  Psi = factor*M(k,i)*Phi(k,l)*M(l,j);
};   
 
template<class T1, class T2> void Var_Carrier::compute_Gamma(T1 &Gamma,
                                                T2 &dM)
{
  TensorA<2> tFr(tFr_in), M(M_in), eFn(eFn_in); 
  Tensor<2> tFreFn = tFr(i,k)*eFn(k,j);
  Tensor<2> G = dM(k,i)*tFreFn(l,k)*tFreFn(l,o)*M(o,j);
  symm(G,Gamma);
};

template<class T> void Var_Carrier::compute_Lambda(double *Lambda,
                                      T &dM)
{
  TensorA<2> M(M_in), pFnp1(pFnp1_in); 
  Tensor<2> pFnp1I,MI;
  inv(pFnp1, pFnp1I);
  inv(M,MI);
  *Lambda = pFnp1I(k,i)*MI(k,j)*dM(i,o)*pFnp1(o,j); 
}; 
   
template<class T> void Var_Carrier::compute_DPsi(T &DPsi,
                                    double *Grad_du_in,
                                    double *Grad_tu_in,
                                    double *dMdu_in,
                                    double *Phi_du_in)
{
  TensorA<2> tFr(tFr_in), M(M_in), eFn(eFn_in);
  
  TensorA<2> Grad_du(Grad_du_in), Grad_tu(Grad_tu_in), dMdu(dMdu_in), Phi_du(Phi_du_in);
  Tensor<2> eFnM = eFn(i,k)*M(k,j);
  Tensor<2> dMPhiM = dMdu(i,k)*Phi_du(k,l)*M(l,j);
  Tensor<2> sdMPhiM;
  symm(dMPhiM, sdMPhiM);

  Tensor<2> GradGrad = Grad_du(k,i)*Grad_tu(k,j);
  Tensor<2> sGradGrad, tFrI, Grad_duFr, sGrad_duFr, FrFr, Grad_tuFr, sGrad_tuFr;

  FrFr = tFr(k,i)*tFr(k,j);
  symm(GradGrad, sGradGrad);
  inv(tFr, tFrI);
  double Grad_tutFrI = Grad_tu(j,i)*tFrI(i,j);
  double Grad_dutFrI = Grad_du(j,i)*tFrI(i,j);
  Grad_duFr = Grad_du(k,i)*tFr(k,j);
  symm(Grad_duFr,sGrad_duFr);

  Grad_tuFr = Grad_tu(k,i)*tFr(k,j);
  symm(Grad_tuFr,sGrad_tuFr);

  double GduFrIGtuFrI = Grad_tu(l,i)*tFrI(i,j)*Grad_tu(j,p)*tFrI(p,l);
        
  Tensor<2> zeta = sGradGrad(i,j) - 2.0/3.0*Grad_tutFrI*sGrad_duFr(i,j)
                 + 2.0/9.0*Grad_tutFrI*Grad_dutFrI*FrFr(i,j)
                 - 1.0/3.0*GduFrIGtuFrI*FrFr(i,j) - 2.0/3.0*Grad_dutFrI*sGrad_tuFr(i,j);

  DPsi = factor*eFnM(k,i)*zeta(k,l)*eFnM(l,j) + 2.0*factor*sdMPhiM(i,j);
};

int compute_Ru(FEMLIB *fe,
               double *Ru,
               Var_Carrier &vc)
{  
  int err = 0;
  
  TensorA<2> tFr(vc.tFr_in), eSd(vc.eSd_in);
  
  Tensor<2> tFrI;
  err += inv(tFr, tFrI);

  for(int ia=0; ia<fe->nne; ia++)
  {
    for(int ib=0; ib<fe->nsd; ib++)
    {
      
      const int id_ab = idx_4_gen(ia,ib,0,0,fe->nne,fe->nsd,fe->nsd,fe->nsd);
      TensorA<2> Grad_du((fe->ST)+id_ab);
      Tensor<2> Phi, Psi;
      vc.compute_Phi(Phi, Grad_du.data);
      vc.compute_Psi(Psi, Phi);
      
      int Ru_id = ia*fe->nsd + ib;              
      Ru[Ru_id] += (Psi(i,j)*eSd(i,j) + 
                   vc.P*vc.tJr*vc.tJn*Grad_du(j,i)*tFrI(i,j))/vc.Jn*fe->detJxW;
    }
  }

  return err;
}

int compute_Rp(FEMLIB *fe,
               double *Rp,
               int Pno,
               double *Np,
               Var_Carrier &vc)
{
  int err = 0;
  
  for(int ia=0; ia<Pno; ia++)
    Rp[ia] += Np[ia]*(vc.tJr*vc.tJn - vc.theta_r*vc.theta_n)/vc.Jn*fe->detJxW;
  
  return err;
}

int compute_Rt(FEMLIB *fe,
               double *Rt,
               int Vno,
               double *Nt,
               Var_Carrier &vc)
{
  int err = 0;
  
  TensorA<2> eSd(vc.eSd_in);
  double factor = 2.0/3.0*pow(1.0/vc.tJr/vc.tJr/vc.theta_r, 2.0/3.0);
  double ZeSd = vc.Z(i,j)*eSd(i,j);
  
  for(int ia=0; ia<Vno; ia++)
  {
    Rt[ia] += Nt[ia]*(factor*ZeSd + vc.dUd_theta*vc.eJn*vc.JM 
                                  - vc.P*vc.theta_n)/vc.Jn*fe->detJxW;
  }  
  return err;
}

int compute_Kuu(FEMLIB *fe,
                double *Kuu,
                double *dMdu_all,
                Var_Carrier &vc)
{
  int err = 0;
  
  TensorA<2> tFr(vc.tFr_in), eSd(vc.eSd_in), eFn(vc.eFn_in); 
  TensorA<4> Ld(vc.Ld_in);
  for(int ia=0; ia<fe->nne; ia++)
  {
    for(int ib=0; ib<fe->nsd; ib++)
    {
      const int id_ab = idx_4_gen(ia,ib,0,0,fe->nne,fe->nsd,fe->nsd,fe->nsd);
      TensorA<2> Grad_du((fe->ST)+id_ab);
      Tensor<2> Phi_du, Psi_du;
      vc.compute_Phi(Phi_du, Grad_du.data);
      vc.compute_Psi(Psi_du, Phi_du);

      for(int iw=0; iw<fe->nne; iw++)
      {
        for(int ig=0; ig<fe->nsd; ig++)
        { 
          const int id_wg = idx_4_gen(iw,ig,0,0,fe->nne,fe->nsd,fe->nsd,fe->nsd);
          TensorA<2> Grad_tu((fe->ST)+id_wg);
          Tensor<2> Phi_tu, Psi_tu,Gamma_tu;
          vc.compute_Phi(Phi_tu, Grad_tu.data);
          vc.compute_Psi(Psi_tu, Phi_tu);

          TensorA<2> dMdu(dMdu_all + id_wg);
          vc.compute_Gamma(Gamma_tu,dMdu);
          double PsiCPsi = Psi_du(i,j)*Ld(i,j,k,l)*(Psi_tu(k,l) + vc.factor*Gamma_tu(k,l));
          
          Tensor<2> tFrI;
          err += inv(tFr,tFrI);
          double Grad_du_tFrI = Grad_du(j,i)*tFrI(i,j);
          double pJrJn = vc.P*vc.tJr*vc.tJn*(Grad_du_tFrI*Grad_tu(j,i)*tFrI(i,j) - Grad_du(j,i)*tFrI(i,k)*Grad_tu(k,l)*tFrI(l,j));

          Tensor<2> DPsi;
          vc.compute_DPsi(DPsi,Grad_du.data,Grad_tu.data,dMdu.data,Phi_du.data);                                                   
          double DPsi_eSd = DPsi(i,j)*eSd(i,j);

          double Lambda;
          vc.compute_Lambda(&Lambda,dMdu);
          double Lambda_Psi = Lambda*(Psi_du(i,j)*eSd(i,j) + vc.P*vc.tJr*vc.tJn*Grad_du(i,j)*tFrI(i,j));

          const int lk_idx = idx_K(ia,ib,iw,ig,fe->nne,fe->nsd);                      
          Kuu[lk_idx] += 1.0/vc.Jn*fe->detJxW*(PsiCPsi + pJrJn + DPsi_eSd + Lambda_Psi);
        }
      }
    }
  }  
  return err;
}

int compute_Kup(FEMLIB *fe,
                double *Kup,
                Var_Carrier &vc,
                int Pno,
                double *Np)
{
  int err = 0;

  TensorA<2> tFr(vc.tFr_in);
  Tensor<2> tFrI;
  err += inv(tFr, tFrI);
  
  for(int ia=0; ia<fe->nne; ia++)
  {
    for(int ib=0; ib<fe->nsd; ib++)
    {
      const int id_ab = idx_4_gen(ia,ib,0,0,fe->nne,fe->nsd,fe->nsd,fe->nsd);
      TensorA<2> Grad_du((fe->ST)+id_ab);
            
      for(int iw=0; iw<Pno; iw++)
      {
        int idx_up = idx_K_gen(ia,ib,iw,0,fe->nne,fe->nsd,Pno,1);
        Kup[idx_up] += 1.0/vc.Jn*fe->detJxW*Grad_du(i,j)*tFrI(i,j)*Np[iw];
      }
    }
  }
  return err;
}

int compute_Kut(FEMLIB *fe,
                double *Kut,
                double *dMdt_all,
                Var_Carrier &vc,
                int Vno,
                double *Nt)
{
  int err = 0;

  TensorA<2> tFr(vc.tFr_in), M(vc.M_in), eSd(vc.eSd_in);
  TensorA<4> Ld(vc.Ld_in);
  Tensor<2> tFrI;  
  err += inv(tFr, tFrI);
         
  double Jtheta = 2.0/3.0*pow(vc.tJr*vc.tJr*vc.theta_r, -1.0/3.0);
  
  for(int ia=0; ia<fe->nne; ia++)
  {
    for(int ib=0; ib<fe->nsd; ib++)
    {
      const int id_ab = idx_4_gen(ia,ib,0,0,fe->nne,fe->nsd,fe->nsd,fe->nsd);
      TensorA<2> Grad_du((fe->ST)+id_ab);
      Tensor<2> Phi_du, Psi_du;
      vc.compute_Phi(Phi_du, Grad_du.data);
      vc.compute_Psi(Psi_du, Phi_du);       
                            
      for(int iw=0; iw<Vno; iw++)
      {
        const int id_wg = idx_4_gen(iw,0,0,0,Vno,1,fe->nsd,fe->nsd);
        TensorA<2> dMdt(dMdt_all + id_wg);
        Tensor<2> MPhidMdt = M(k,i)*Phi_du(k,l)*dMdt(l,j);
        
        Tensor<2> DPsi = Jtheta*Psi_du(i,j)/vc.factor
                       + 2.0*vc.factor*0.5*(MPhidMdt(i,j) + MPhidMdt(j,i));
        
        double DPsi_eSd = DPsi(i,j)*eSd(i,j);
        
        Tensor<2> Gamma_tu;
        vc.compute_Gamma(Gamma_tu,dMdt);
        Tensor<2> DeSd = 0.5*Ld(i,j,k,l)*(Jtheta*vc.Z(k,l) + 2.0*vc.factor*Gamma_tu(k,l));
        
        double PsiDeSd = Psi_du(i,j)*DeSd(i,j);
        
        double Lambda;
        vc.compute_Lambda(&Lambda, dMdt);
        
        double Lambda_Psi = Lambda*(Psi_du(i,j)*eSd(i,j) 
                          + vc.P*vc.tJr*vc.tJn*Grad_du(j,i)*tFrI(i,j));
        int idx_ut = idx_K_gen(ia,ib,iw,0,fe->nne,fe->nsd,Vno,1);
        Kut[idx_ut] += 1.0/vc.Jn*fe->detJxW*(DPsi_eSd + PsiDeSd - Lambda_Psi)*Nt[iw];
      }
    }
  }
  return err;
}

int compute_Ktu(FEMLIB *fe,
                double *Ktu,
                double *dMdu_all,
                Var_Carrier &vc,
                int Vno,
                double *Nt)
{
  int err = 0;

  TensorA<2> tFr(vc.tFr_in), M(vc.M_in), eSd(vc.eSd_in);
  TensorA<4> Ld(vc.Ld_in);
  
  Tensor<2> tFrI;
  err += inv(tFr, tFrI);
  
  Tensor<2> MI;
  err += inv(M, MI);
  
  double Jtheta = 2.0/3.0*pow(vc.tJr*vc.tJr*vc.theta_r, -1.0/3.0);
  double eJnJM = vc.eJn*vc.JM;
  
  for(int ia=0; ia<Vno; ia++)
  {
    for(int iw=0; iw<fe->nne; iw++)
    {
      for(int ig=0; ig<fe->nsd; ig++)
      { 
        const int id_wg = idx_4_gen(iw,ig,0,0,fe->nne,fe->nsd,fe->nsd,fe->nsd);

        TensorA<2> Grad_tu((fe->ST)+id_wg);
        Tensor<2> Phi_tu, Psi_tu,Gamma_tu;
        vc.compute_Phi(Phi_tu, Grad_tu.data);
        vc.compute_Psi(Psi_tu, Phi_tu);

        TensorA<2> dMdu(dMdu_all + id_wg);
        vc.compute_Gamma(Gamma_tu,dMdu);
        
        double DZ_eSd = (2.0*Jtheta/vc.theta_r*Psi_tu(i,j) + 2.0*Jtheta*Gamma_tu(i,j))*eSd(i,j);
        double ZCPsi = Jtheta*vc.Z(i,j)*Ld(i,j,k,l)*(Psi_tu(k,l) + Jtheta*Gamma_tu(k,l));
        double JMdMdu = eJnJM*(vc.d2Ud_theta2*vc.theta_r*eJnJM + vc.dUd_theta)*MI(j,i)*dMdu(i,j);

        double Lambda;
        vc.compute_Lambda(&Lambda, dMdu);        

        double Lmabda_ZeSd = Lambda*Jtheta*vc.Z(i,j)*eSd(i,j);
        Lmabda_ZeSd += Lambda*(vc.dUd_theta*eJnJM - vc.P*vc.theta_n);        
        
        int idx_tu = idx_K_gen(ia,0,iw,0,Vno,1,fe->nne,fe->nsd);
        Ktu[idx_tu] += Nt[ia]/vc.Jn*fe->detJxW*(DZ_eSd + ZCPsi + JMdMdu - Lmabda_ZeSd);
      }
    }
  }
  return err;
}


int compute_Ktp(FEMLIB *fe,
                double *Ktp,
                Var_Carrier &vc,
                int Vno,
                double *Nt,
                int Pno,
                double *Np)
{
  int err = 0;  
    
  for(int ia=0; ia<Vno; ia++)
  {
    for(int iw=0; iw<Pno; iw++)
    {
      int idx_tp = idx_K_gen(ia,0,iw,0,Vno,1,Pno,1);
      Ktp[idx_tp] -= Nt[ia]*Np[iw]/vc.Jn*fe->detJxW*vc.theta_n;
    }
  }
  return err;
}

int compute_Ktt(FEMLIB *fe,
                double *Ktt,
                double *dMdt_all,
                Var_Carrier &vc,
                int Vno,
                double *Nt)
{
  int err = 0;
  
  TensorA<2> tFr(vc.tFr_in), M(vc.M_in), eSd(vc.eSd_in);
  TensorA<4> Ld(vc.Ld_in);
  Tensor<2> tFrI;
  err += inv(tFr, tFrI);
  
  Tensor<2> MI;
  err += inv(M, MI);
         
  double Jtheta = 2.0/3.0*pow(vc.tJr*vc.tJr*vc.theta_r, -1.0/3.0);
  double eJnJM = vc.eJn*vc.JM;
  
  for(int ia=0; ia<Vno; ia++)
  {
    for(int iw=0; iw<Vno; iw++)
    {
      const int id_wg = idx_4_gen(iw,0,0,0,Vno,1,fe->nsd,fe->nsd);
      TensorA<2> dMdt(dMdt_all + id_wg);
      Tensor<2> Gamma_tt;
      vc.compute_Gamma(Gamma_tt,dMdt);
      double ZGamm_eSd = Jtheta*(-1.0/3.0*vc.Z(i,j) + 2.0*Gamma_tt(i,j))*eSd(i,j);
      
      Tensor<2> DeSd = 0.5*Ld(i,j,k,l)*(Jtheta*vc.Z(k,l) + 2.0*vc.factor*Gamma_tt(k,l));
      double ZDeSd = Jtheta*vc.Z(i,j)*DeSd(i,j);
      double JUdMdt = eJnJM*vc.dUd_theta*MI(j,i)*dMdt(i,j);
      JUdMdt += eJnJM*vc.d2Ud_theta2*vc.theta_r*eJnJM;
      
      double Lambda;
      vc.compute_Lambda(&Lambda, dMdt);
      
      double Lambda_ZeSd = Lambda*Jtheta*vc.Z(i,j)*eSd(i,j);
      Lambda_ZeSd += Lambda*(vc.dUd_theta*eJnJM - vc.P*vc.theta_n);
      
      int idx_tt = idx_K_gen(ia,0,iw,0,Vno,1,Vno,1);
      Ktt[idx_tt] += 1.0/vc.Jn*fe->detJxW*Nt[ia]*Nt[iw]*(ZGamm_eSd + ZDeSd + JUdMdt - Lambda_ZeSd);
    }
  }
  return err;
}


int condense_K_3F_to_1F(double *K, int nne, int nsd, int Pno, int Vno,
                        double *Kuu_in, double *Kut_in, double *Kup_in,
                        double *Ktu_in, double *Ktt_in, double *Ktp_in,
                        double *Kpu_in, double *Kpt_in, double *Kpp_in)                               
{
  int err = 0;
  Matrix<double> Kuu, Kut, Kup;
  Matrix<double> Ktu, Ktt, Ktp;  
  Matrix<double> Kpu, Kpt, Kpp;  
  
	Kuu.use_reference(nne*nsd, nne*nsd, Kuu_in);
	Ktu.use_reference(Vno,     nne*nsd, Ktu_in);
	Kpu.use_reference(Pno,     nne*nsd, Kpu_in);
	Kut.use_reference(nne*nsd, Vno,     Kut_in);
	Ktt.use_reference(Vno,     Vno,     Ktt_in);
	Kpt.use_reference(Pno,     Vno,     Kpt_in);	
	Kup.use_reference(nne*nsd, Pno,     Kup_in);
	Ktp.use_reference(Vno,     Pno,     Ktp_in);
	Kpp.use_reference(Pno,     Pno,     Kpp_in);	
   

	Matrix<double> KptI(Vno, Pno), KtpI(Pno, Vno), Kuu_add(nne*nsd,nne*nsd);
	
	KptI.inv(Kpt);
	KtpI.inv(Ktp);

  Matrix<double> KupKtpI(nne*nsd,Vno), KupKtpIKtt(nne*nsd,Vno), KptIKpu(Vno,nne*nsd);
           
     KupKtpI.prod(Kup,KtpI);
  KupKtpIKtt.prod(KupKtpI,    Ktt    );  
     KptIKpu.prod(KptI,       Kpu    );  
     Kuu_add.prod(KupKtpIKtt, KptIKpu);

  Matrix<double> KupKtpIKtu(nne*nsd,nne*nsd);
  KupKtpIKtu.prod(KupKtpI, Ktu); 
  
  Matrix<double> KutKptIKpu(nne*nsd,nne*nsd);
  KutKptIKpu.prod(Kut,KptIKpu);
  
  for(int a=0; a<nne*nsd*nne*nsd; a++)
    K[a] = Kuu.m_pdata[a] - KupKtpIKtu.m_pdata[a] + Kuu_add.m_pdata[a] - KutKptIKpu.m_pdata[a];
    
  return err;  
}

int condense_F_3F_to_1F(double *fe, int nne, int nsd, int Pno, int Vno,
                        double *fu_in, double *ft_in, double *fp_in, 
                        double *Kut_in, double *Kup_in, double *Ktp_in, double *Ktt_in,double *Kpt_in)
{
  int err = 0;
  
  Matrix<double> fu,ft,fp;
	fu.use_reference(nne*nsd, 1, fu_in);	
	ft.use_reference(Vno,     1, ft_in);	
	fp.use_reference(Pno,     1, fp_in);	
	
  Matrix<double> Kut,Kup,Ktp,Ktt,Kpt;
  Kut.use_reference(nne*nsd, Vno,  Kut_in);	
  Kup.use_reference(nne*nsd, Pno,  Kup_in);	
  Ktp.use_reference(Vno,     Pno,  Ktp_in);
  Ktt.use_reference(Vno,     Vno,  Ktt_in); 	
  Kpt.use_reference(Pno,     Vno,  Kpt_in);

	Matrix<double> KptI(Vno,Pno), KtpI(Pno,Vno), fu_add(nne*nsd,1,0.0);

  KtpI.inv(Ktp);
  KptI.inv(Kpt);

  Matrix<double> KupKtpI(nne*nsd,Vno), KupKtpIKtt(nne*nsd,Vno), KptIFp(Vno,1), KupKtpIFt(nne*nsd,1), KutKptIFp;
  
     KupKtpI.prod(Kup,        KtpI  );
  KupKtpIKtt.prod(KupKtpI,    Ktt   );
      KptIFp.prod(KptI,       fp    );
      fu_add.prod(KupKtpIKtt, KptIFp);
   KupKtpIFt.prod(KupKtpI,    ft);
   KutKptIFp.prod(Kut,        KptIFp);
              
  for(int a=0; a<nne*nsd; a++)
    fe[a] -=  (-fu.m_pdata[a] + fu_add.m_pdata[a] - KupKtpIFt.m_pdata[a] - KutKptIFp.m_pdata[a]);

  return err;
}


int compute_d_theta_dP(double *d_theta, double *dP, double *du_in, 
                       int nne, int nsd, int Pno, int Vno,
                       double *fu_in, double *ft_in, double *fp_in, 
                       double *Kpu_in, double *Ktu_in, double *Ktp_in, double *Ktt_in,double *Kpt_in)
{
  int err = 0;
  
  Matrix<double> du;
  du.use_reference(nne*nsd, 1, du_in);	
  	
  Matrix<double> fu,ft,fp;
	fu.use_reference(nne*nsd, 1, fu_in);	
	ft.use_reference(Vno,     1, ft_in);	
	fp.use_reference(Pno,     1, fp_in);
  
  Matrix<double> Kpu,Ktu,Ktp,Ktt,Kpt;
  Kpu.use_reference(Pno, nne*nsd, Kpu_in);
  Ktu.use_reference(Vno, nne*nsd, Ktu_in);  
  Ktp.use_reference(Vno,     Pno, Ktp_in);
  Ktt.use_reference(Vno,     Vno, Ktt_in);   
  Kpt.use_reference(Pno,     Vno, Kpt_in);

	Matrix<double> KptI(Vno,Pno), KptIFp(Vno,1), KptIKpu_du(Vno,1);
  KptI.inv(Kpt);
  
  KptIFp.prod(KptI,fp);
  KptIKpu_du.prod(KptI, Kpu, du);
              
  for(int ia=0; ia<Vno; ia++)
    d_theta[ia] = -KptIFp.m_pdata[ia] - KptIKpu_du.m_pdata[ia];
        
	Matrix<double> KtpI(Pno,Vno), KtpIFt(Pno, 1), KtpIKtu_du(Pno,1);
	Matrix<double> KtpIKtt_d_theta(Pno, 1);
  KtpI.inv(Ktp);
  
  KtpIFt.prod(KtpI, ft);
  KtpIKtu_du.prod(KtpI, Ktu,du);

  Matrix<double> d_t;
  d_t.use_reference(Vno, 1, d_theta);
  KtpIKtt_d_theta.prod(KtpI, Ktt, d_t);  
  
  for(int ia=0; ia<Pno; ia++)
    dP[ia] = -KtpIFt.m_pdata[ia] - KtpIKtu_du.m_pdata[ia] - KtpIKtt_d_theta.m_pdata[ia];
    
  return err;
}