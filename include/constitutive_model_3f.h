/// Define the interface to generalized constitutive models with
/// their associated data structure(s). The user is responsible for
/// defining and linking the functions associated with a particular
/// model. NEED MORE USAGE DETAILS.
/// 
/// Authors:
///  Sangmin Lee, [1], <slee43@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN

#ifndef H__CONSTITUTIVE_MODEL_3F__H
#define H__CONSTITUTIVE_MODEL_3F__H

#include "femlib.h"
#include <ttl/ttl.h>

class Var_Carrier
{
  public:
    double *tFr_in, *eFn_in, *M_in, *pFnp1_in, *eSd_in, *Ld_in;
    ttl::Tensor<2,3,double> Z;    

    double tJr, tJn, Jn, theta_r, theta_n, eJn, JM;
    double P;        // pressure    
    double factor;   // factor = pow(theta_r/tJr, 2.0/3.0);
    double dUd_theta;
    double d2Ud_theta2;

    Var_Carrier()
    {
      tFr_in   = NULL;
      eFn_in   = NULL;
      M_in     = NULL;
      pFnp1_in = NULL;
      eSd_in   = NULL;
      Ld_in    = NULL;
    };

    ~Var_Carrier()
    {
      tFr_in   = NULL;
      eFn_in   = NULL;
      M_in     = NULL;
      pFnp1_in = NULL;
      eSd_in   = NULL;
      Ld_in    = NULL;
    };

    void compute_factor(void);
    void set_tenosrs(double *tFr,
                     double *eFn,
                     double *M,
                     double *pFnp1,
                     double *eSd,
                     double *Ld);
    
    void set_scalars(double theta_r_in,
                     double theta_n_in,
                     double tJn_in,
                     double Jn_in,
                     double P_in,
                     double dUd_theta_in,
                     double d2Ud_theta2_in);
        
    void compute_Z(void);
    template<class T1> void compute_Phi(T1 &Phi, double *Grad_beta_in);
    template<class T1, class T2> void compute_Psi(T1 &Psi, T2 &Phi);
    template<class T1, class T2> void compute_Gamma(T1 &Gamma,
                                                    T2 &dM);
    template<class T> void compute_Lambda(double *Lambda,
                                          T &dM);
    template<class T> void compute_DPsi(T &DPsi,
                                        double *Grad_du_in,
                                        double *Grad_tu_in,
                                        double *dMdu_in,
                                        double *Phi_du_in);
};

int compute_Ru(FEMLIB *fe,
               double *Ru,
               Var_Carrier &vc);

int compute_Rp(FEMLIB *fe,
               double *Rp,
               int Pno,
               double *Np,               
               Var_Carrier &vc);

int compute_Rt(FEMLIB *fe,
               double *Rt,
               int Vno,
               double *Nt,               
               Var_Carrier &vc);
               
int compute_Kuu(FEMLIB *fe,
                double *Kuu,
                double *dMdu_all,
                Var_Carrier &vc);

int compute_Kup(FEMLIB *fe,
                double *Kup,
                Var_Carrier &vc,
                int Pno,
                double *Np);

int compute_Kut(FEMLIB *fe,
                double *Kut,
                double *dMdt_all,
                Var_Carrier &vc,
                int Vno,
                double *Nt);

int compute_Ktu(FEMLIB *fe,
                double *Ktu,
                double *dMdu_all,
                Var_Carrier &vc,
                int Vno,
                double *Nt);


int compute_Ktp(FEMLIB *fe,
                double *Ktp,
                Var_Carrier &vc,
                int Vno,
                double *Nt,
                int Pno,
                double *Np);

int compute_Ktt(FEMLIB *fe,
                double *Ktt,
                double *dMdt_all,
                Var_Carrier &vc,
                int Vno,
                double *Nt);

int condense_K_3F_to_1F(double *Ks, int nne, int nsd, int Pno, int Vno,
                        double *Kuu_in, double *Kut_in, double *Kup_in,
                        double *Ktu_in, double *Ktt_in, double *Ktp_in,
                        double *Kpu_in, double *Kpt_in, double *Kpp_in);

int condense_F_3F_to_1F(double *fe, int nne, int nsd, int Pno, int Vno,
                        double *fu_in, double *ft_in, double *fp_in, 
                        double *Kut_in, double *Kup_in, double *Ktp_in, double *Ktt_in,double *Kpt_in);

int compute_d_theta_dP(double *d_theta, double *dP, double *du_in, 
                       int nne, int nsd, int Pno, int Vno,
                       double *fu_in, double *ft_in, double *fp_in, 
                       double *Kpu_in, double *Ktu_in, double *Ktp_in, double *Ktt_in,double *Kpt_in);
#endif
