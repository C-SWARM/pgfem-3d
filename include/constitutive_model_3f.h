/// Define the functions for computing stiffness matrix
/// and residulas for constitutive models with three field mixed method.
/// 
/// Authors:
///  Sangmin Lee, [1], <slee43@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN

#ifndef H__CONSTITUTIVE_MODEL_3F__H
#define H__CONSTITUTIVE_MODEL_3F__H

#include "data_structure.h"
#include "femlib.h"
#include <ttl/ttl.h>

/// Variables at integration point for Updated Lagrangian 
/// and three field mixed method formulation
class CM_ThreeField
{
  public:
    FEMLIB *fe; // reference of femlib object
    int Vno;    // number of volume variables
    int Pno;    // number of pressrue variables
    double *Nt; // finite element shape function for volume
    double *Np; // finite element shape function for pressure
    
    double *tFr_in, *eFn_in, *M_in, *pFnp1_in, *eSd_in;
    double *Ld_in;

    double tJr, tJn, Jn, theta_r, theta_n, eJn, JM;
    double P;        // pressure    
    double factor;   // factor = pow(theta_r/tJr, 2.0/3.0);
    double dUd_theta;
    double d2Ud_theta2;
    double dt_alpha_1_minus_alpha;

    CM_ThreeField()
    {
      Vno = 1;
      Pno = 1;
      dt_alpha_1_minus_alpha = -1.0;
    };

    ~CM_ThreeField()
    {
      fe = NULL;
      Nt = NULL;
      Np = NULL;
    };

    void compute_factor(void);
    void set_tenosrs(double *tFr_in,
                     double *eFn_in,
                     double *M_in,
                     double *pFnp1_in,
                     double *eSd_in,
                     double *Ld_in);
    
    void set_scalars(double theta_r_in,
                     double theta_n_in,
                     double tJn_in,
                     double Jn_in,
                     double P_in,
                     double dUd_theta_in,
                     double d2Ud_theta2_in,
                     double dt_alpha_1_minus_alpha_in = -1.0);
                     
    void set_femlib(FEMLIB *fe_in,
                    int     Vno_in,
                    int     Pno_in,
                    double *Nt_in,
                    double *Np_in);                
        
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
                                        
    /// compute residual: Ru at ip
    int compute_Ru(double *Ru);
    
    /// compute residual: Rp at ip
    int compute_Rp(double *Rp);
    
    /// compute residual: Rt at ip
    int compute_Rt(double *Rt);
    
    /// compute stiffness: Kpt at ip               
    int compute_Kuu(double *Kuu,
                    double *dMdu_all);
    
    /// compute stiffness: Kpt at ip
    int compute_Kup(double *Kup);
                    
    /// compute stiffness: Kpt at ip
    int compute_Kut(double *Kut,
                    double *dMdt_all);
    
    /// compute stiffness: Kpt at ip
    int compute_Ktu(double *Ktu,
                    double *dMdu_all);
    
    /// compute stiffness: Kpt at ip
    int compute_Ktp(double *Ktp);
    
    /// compute stiffness: Kpt at ip
    int compute_Ktt(double *Ktt,
                    double *dMdt_all);                                        
};

class ThreeFieldStiffness
{
  public:
    gcm::Matrix<double> Kuu, Kut, Kup, Ktu, Ktt, Ktp, Kpu, Kpt, Kpp;
    bool is_for_residual;
    
    ThreeFieldStiffness(void){ is_for_residual = false; };
    ThreeFieldStiffness(FEMLIB *fe, 
                        int Vno, 
                        int Pno,
                        bool for_residual = false)
    { initialization(fe, Vno, Pno, for_residual);};

    void initialization(FEMLIB *fe, 
                        int Vno, 
                        int Pno, 
                        bool for_residual = false)
    {
      is_for_residual = for_residual;
      int nsd = fe->nsd;
      int nne = fe->nne;
      if(is_for_residual)
      {
        Kut.initialization(nne*nsd,Vno    ,0.0);
        Ktu.initialization(Vno,nne*nsd    ,0.0);
        Kup.initialization(nne*nsd,Pno    ,0.0);
        Ktt.initialization(Vno    ,Vno    ,0.0);
        Ktp.initialization(Vno    ,Pno    ,0.0);
        Kpt.initialization(Pno    ,Vno    ,0.0);
      }  
      else
      {    
        Kuu.initialization(nne*nsd,nne*nsd,0.0);
        Kut.initialization(nne*nsd,Vno    ,0.0);
        Kup.initialization(nne*nsd,Pno    ,0.0);
        Ktu.initialization(Vno    ,nne*nsd,0.0);
        Ktt.initialization(Vno    ,Vno    ,0.0);
        Ktp.initialization(Vno    ,Pno    ,0.0);
        Kpu.initialization(Pno    ,nne*nsd,0.0);
        Kpt.initialization(Pno    ,Vno    ,0.0);
        Kpp.initialization(Pno    ,Pno    ,0.0);
      }    
    };
    int compute_stiffness(CM_ThreeField &cm_tf, 
                          gcm::Matrix<double> &dMdu, 
                          gcm::Matrix<double> &dMdt)
    {
      int err = 0;
      if(is_for_residual)
      {
        err += cm_tf.compute_Kut(Kut.m_pdata, dMdt.m_pdata);
        err += cm_tf.compute_Ktu(Ktu.m_pdata, dMdu.m_pdata);
        err += cm_tf.compute_Kup(Kup.m_pdata);
        err += cm_tf.compute_Ktp(Ktp.m_pdata);
        err += cm_tf.compute_Ktt(Ktt.m_pdata, dMdt.m_pdata);
      }
      else
      {
        err += cm_tf.compute_Kuu(Kuu.m_pdata, dMdu.m_pdata);
        err += cm_tf.compute_Kut(Kut.m_pdata, dMdt.m_pdata);
        err += cm_tf.compute_Kup(Kup.m_pdata);
        err += cm_tf.compute_Ktu(Ktu.m_pdata, dMdu.m_pdata);
        err += cm_tf.compute_Ktp(Ktp.m_pdata);
        err += cm_tf.compute_Ktt(Ktt.m_pdata, dMdt.m_pdata);
      }        
      return err;      
    };
};

class ThreeFieldResidual
{
  public:
    gcm::Matrix<double> Ru, Rp, Rt;
    bool is_for_PT_update;

    ThreeFieldResidual(){ is_for_PT_update = false; };
    ThreeFieldResidual(FEMLIB *fe, 
                       int Vno, 
                       int Pno,
                       bool for_PT_update = false)
    { initialization(fe, Vno, Pno, for_PT_update); };
                           
    void initialization(FEMLIB *fe, 
                        int Vno, 
                        int Pno,
                        bool for_PT_update = false)
    {
      is_for_PT_update = for_PT_update;
      int nsd = fe->nsd;
      int nne = fe->nne;
      if(is_for_PT_update)
      {
        Rt.initialization(Vno,1, 0.0);
        Rp.initialization(Pno,1, 0.0);
      }
      else
      {
        Ru.initialization(nne*nsd,1,0.0);
        Rt.initialization(Vno,    1, 0.0);
        Rp.initialization(Pno,    1, 0.0);
      }
    };
    int compute_residual(CM_ThreeField &cm_tf)
    {
      int err = 0;
      if(is_for_PT_update)
      {      
        err += cm_tf.compute_Rp(Rp.m_pdata);
        err += cm_tf.compute_Rt(Rt.m_pdata);
      }
      else
      {
        err += cm_tf.compute_Ru(Ru.m_pdata);
        err += cm_tf.compute_Rp(Rp.m_pdata);
        err += cm_tf.compute_Rt(Rt.m_pdata);        
      }
      return err;
    };                             
};

/// compute and update increments of prssure and volume for transient
int compute_d_theta_dP(double *d_theta_in, double *dP_in, double *du_in, 
                       int nne, int nsd, int Pno, int Vno,
                       double *fu_in, double *ft_in, double *fp_in, 
                       double *Kpu_in, double *Ktu_in, double *Ktp_in, double *Ktt_in,double *Kpt_in);
#endif
