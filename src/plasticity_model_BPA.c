/**
 * This file defines the implementation for the BPA plasticity model.
 * REFERENCES:
 *
 * Mary C. BOYCE, David M. PARKS, Ali S. ARGON (1988). LARGE INELASTIC
 * DEFORMATION OF GLASSY POLYMERS.PART I: RATE DEPENDENT CONSTITUTIVE
 * MODEL. Mechanics of Materials, 7:15-33.
 *
 * Holopainen, S. (2013). Modeling of the mechanical behavior of
 * amorphous glassy polymers under variable loadings and comparison
 * with state-of-the-art model predictions. Mechanics of Materials,
 * 66:35–58.
 *
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */

#include "plasticity_model_BPA.h"
#include "_plasticity_model_BPA.h"
#include "constitutive_model.h"
#include "state_variables.h"
#include "new_potentials.h"
#include "data_structure_c.h"
#include "utils.h"
#include "index_macros.h"

/* Define constant dimensions. Note cannot use `static const` with
   initialization list */
#define dim  3
#define tensor 9
#define tensor4 81
#define tangent 19

Define_Matrix(double);

/* constants/enums */
static const int _n_Fs = 6;
static const int _n_vars = 4;
enum {_M,_W,_Fe,_M_n,_W_n,_Fe_n};
enum {_s,_lam,_s_n,_lam_n};
static const double eye[tensor] = {1.0,0,0, 0,1.0,0, 0,0,1.0};

/* material parameters (constant for testing purposes, need to be
   migrated to Model_parameters in a programatic/general way)*/
/* parameter values from S. Holopanien, Mech. of Mat. (2013) */
static const double param_A = 289;
static const double param_T = 240;
static const double param_N = 2.15;
static const double param_Cr = 12.8;
static const double param_alpha = 00.8;
static const double param_gdot0 = 2.0e15;
static const double param_h = 500;
static const double param_s0 = 97;
static const double param_s_ss = 76.6;

/*
 * Purely static functions
 */
static double bpa_compute_bulk_mod(const HOMMAT *mat)
{
  return ( (2* mat->G * (1 + mat->nu)) / (3 * (1 - 2 * mat->nu)) );
}

static void bpa_compute_Fe(const double * restrict F,
                           const double * restrict M,
                           double * restrict Fe)
{
  memset(Fe, 0, tensor * sizeof(*Fe));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
       for (int k = 0; k < dim; k++) {
        Fe[idx_2(i,j)] += F[idx_2(i,k)] * M[idx_2(k,j)];
      }
    }
  }
}

static double bpa_compute_plam(const double *Cp)
{
  return sqrt(1./3. * (Cp[0] + Cp[4] + Cp[8]));
}

static void bpa_compute_Cp(double * restrict Cp,
                           const double * restrict Fp)
{
  memset(Cp, 0, tensor * sizeof(*Cp));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        Cp[idx_2(i,j)] += Fp[idx_2(i,k)] * Fp[idx_2(j,k)];
      }
    }
  }
}

static void bpa_compute_Cpdev(double * restrict Cpdev,
                              const double * restrict Cp)
{
  const double Jdev = pow(det3x3(Cp),-1./3.);
  for (int i = 0; i < tensor; i++) {
    Cpdev[i] = Jdev * Cp[i];
  }
}

static void bpa_compute_Fe_Ce(const double * restrict F,
                              const double * restrict M,
                              double * restrict Fe,
                              double * restrict Ce)
{
  bpa_compute_Fe(F,M,Fe);
  memset(Ce, 0, tensor * sizeof(*Ce));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        Ce[idx_2(i,j)] += Fe[idx_2(k,i)] * Fe[idx_2(k,j)];
      }
    }
  }
}

static void bpa_compute_Sdev(const double *Ce,
                             const HOMMAT *p_hmat,
                             double *Sdev)
{
  devStressFuncPtr Stress = getDevStressFunc(-1,p_hmat);
  Stress(Ce,p_hmat,Sdev);
}

static void bpa_compute_Ldev(const double *Ce,
                             const HOMMAT *p_hmat,
                             double *Ldev)
{
  matStiffFuncPtr Tangent = getMatStiffFunc(-1,p_hmat);
  Tangent(Ce,p_hmat,Ldev);
}

static void bpa_compute_dudj(const double Je,
                             const HOMMAT *p_hmat,
                             double *dudj)
{
  dUdJFuncPtr Pressure = getDUdJFunc(-1,p_hmat);
  Pressure(Je,p_hmat,dudj);
}

static void bpa_compute_d2udj2(const double Je,
                               const HOMMAT *p_hmat,
                               double *d2udj2)
{
  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,p_hmat);
  D_Pressure(Je,p_hmat,d2udj2);
}

static int bpa_compute_step1_terms(double *gdot,
                                   double *s_s,
                                   double *tau,
                                   double *eq_sig_dev,
                                   double *normal,
                                   const double s,
                                   const double kappa,
                                   const double *F,
                                   const double *M,
                                   const HOMMAT *p_hmat)
{
  int err = 0;
  double Fp[tensor] = {};
  double Fe[tensor] = {};
  double Ce[tensor] = {};
  err += inv3x3(M,Fp);
  bpa_compute_Fe_Ce(F,M,Fe,Ce);

  /* compute the elastic stress */
  double Sdev[tensor] = {};
  bpa_compute_Sdev(Ce,p_hmat,Sdev);

  /* compute the pressure */
  double pressure = 0.0;
  double Je = det3x3(Fe);
  bpa_compute_dudj(Je,p_hmat,&pressure);
  pressure *= kappa / 2.0;

  /* compute the pressure-dependent athermal shear stress */
  *s_s = s + param_alpha * pressure;

  /* compute the plastic backstress */
  double Bdev[tensor] = {};
  err += BPA_compute_Bdev(Bdev,Fp);
  err += BPA_compute_loading_dir(normal,eq_sig_dev,tau,Sdev,Bdev,Fe);
  err += BPA_compute_gdot(gdot,*tau,*s_s);
  return err;
}

/**
 * Compute an initial guess for the solution vector (M,Wp,lam). The
 * initial guess for the plastic deformation (M) is that the increment
 * of deformation (Fn -> Fn+1) is all plastic. This is reasonable as
 * the elastic deformations are assumed to be small.
 *
 * Note that we use the most up-to-date estimates (i.e., from previous
 * iterations in the PDE solve) to compute the initial guess.
 */
int bpa_int_alg_initial_guess(const double *F,
                              const Constitutive_model *m,
                              double * restrict M,
                              double *Wp,
                              double *lam,
                              double *s)
{
  int err = 0;

  /* M ~= inv(F) Fe */
  memset(M, 0, tensor * sizeof(*M));
  double invF[tensor] = {};
  err += inv3x3(F, invF);
  const double * restrict Fen = m->vars.Fs[_Fe].m_pdata;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        M[idx_2(i,j)] += invF[idx_2(i,k)] * Fen[idx_2(k,j)];
      }
    }
  }

  /* Wp ~= Wp */
  memcpy(Wp, m->vars.Fs[_W].m_pdata, tensor * sizeof(*Wp));

  /* s ~= s, lam ~= lam */
  *s = m->vars.state_vars->m_pdata[_s];
  *lam = m->vars.state_vars->m_pdata[_lam];

  return err;
}

int bpa_update_state_variables(const double *F,
                               const double *M,
                               const double *Wp,
                               const double lam,
                               const double s,
                               Constitutive_model *m)
{
  int err = 0;
  bpa_compute_Fe(F,M,m->vars.Fs[_Fe].m_pdata);
  memcpy(m->vars.Fs[_M].m_pdata, M, tensor * sizeof(*M));
  memcpy(m->vars.Fs[_W].m_pdata, Wp, tensor * sizeof(*Wp));
  m->vars.state_vars->m_pdata[_lam] = lam;
  m->vars.state_vars->m_pdata[_s] = s;
  return err;
}

/**
 * Update the current value of the solution with the increment
 *
 */
int bpa_update_solution(const double * restrict increment,
                        double * restrict M,
                        double * restrict Wp,
                        double * restrict lam)
{
  int err = 0;
  for (int i = 0; i < tensor; i++) {
    M[i] += increment[i];
    Wp[i] += increment[i + tensor];
  }
  *lam += increment[2 * tensor];
  return err;
}

/*
 * Private interface for the model
 */
int BPA_int_alg(Constitutive_model *m,
                const void *ctx)
{
  static double TOL1 = 1.0e-5;
  static double TOL2 = 1.0e-5;
  static int maxit1 = 10;
  static int maxit2 = 5;

  int err = 0;
  const BPA_ctx *CTX = ctx;
  const HOMMAT *p_hmat = m->param->p_hmat;
  double Fpn[tensor] = {};
  err += inv3x3(m->vars.Fs[_M_n].m_pdata, Fpn);
  const double kappa = bpa_compute_bulk_mod(p_hmat);
  const double s_n = m->vars.state_vars->m_pdata[_s_n];

  /* initial guess for the solution */
  double M[tensor] = {};
  double Wp[tensor] = {};
  double lam = 0;
  double s = 0;
  err += bpa_int_alg_initial_guess(CTX->F, m, M, Wp, &lam, &s);

  /* internal variables */
  double gdot = 0;
  double s_s = 0;
  double tau = 0;
  double eq_sig_dev[tensor] = {};
  double normal[tensor] = {};

  double TAN[tangent * tangent] = {};
  double RES[tangent] = {};
  double norm1 = 1.0;
  double norm2 = 1.0;
  int iter2 = 0;
  err += bpa_compute_step1_terms(&gdot, &s_s, &tau, eq_sig_dev, normal, /*< out */
                                 s, kappa, CTX->F, M, p_hmat);          /*< in */
  err += BPA_int_alg_res(RES, CTX->dt, gdot, lam, Fpn,
                         normal, CTX->F, M, Wp);

  while ((norm1 > TOL1 || norm2 > TOL2) && iter2 < maxit2) {

    /* compute plstic deformation */
    int iter1 = 0;
    while (norm1 > TOL1 && iter1 < maxit1) {
      err += BPA_int_alg_tan(TAN, CTX->dt, gdot, lam, tau, s_s, Fpn,
                             normal, eq_sig_dev, CTX->F, M, Wp, p_hmat);

      /* solve for increment: note solution in RES on exit */
      err += solve_Ax_b(tangent,TAN,RES);

      /* update deformation */
      err += bpa_update_solution(RES, M, Wp, &lam);

      /* update vars (gdot, s_s, tau, eq_sig_dev, normal) */
      err += bpa_compute_step1_terms(&gdot, &s_s, &tau, eq_sig_dev, normal, /*< out */
                                     s, kappa, CTX->F, M, p_hmat);          /*< in */

      /* compute residual */
      err += BPA_int_alg_res(RES, CTX->dt, gdot, lam, Fpn,
                             normal, CTX->F, M, Wp);
      norm1 = cblas_dnrm2(tangent,RES,1);
      iter1++;
    }

    /* update s */
    const double s_k = s;
    s = s_s * (s_n + CTX->dt * param_h * gdot) / (s_s + CTX->dt * param_h * gdot);
    norm2 = fabs(s - s_k) / param_s_ss;

    /* compute the residual with updated value of s */
    err += bpa_compute_step1_terms(&gdot, &s_s, &tau, eq_sig_dev, normal, /*< out */
                                   s, kappa, CTX->F, M, p_hmat);          /*< in */
    err += BPA_int_alg_res(RES, CTX->dt, gdot, lam, Fpn,
                           normal, CTX->F, M, Wp);
    norm1 = cblas_dnrm2(tangent,RES,1);

    /* output of the iterative procedure */
    printf("[%d] R1 = %6e (%d) || R2 = %6e\n", iter2, norm1, iter1, norm2);
  }

  /* Update state variables with converged values */
  err += bpa_update_state_variables(CTX->F,M, Wp, lam, s, m);
  return err;
}

int BPA_int_alg_res_terms(double * restrict Rm,
                          double * restrict Rw,
                          double * restrict Rlam,
                          const double dt,
                          const double gdot,
                          const double lam,
                          const double * restrict Fpn,
                          const double * restrict n,
                          const double * restrict F,
                          const double * restrict M,
                          const double * restrict Wp)
{
  int err = 0;

  /* compute intermediate terms (mat mults) */
  double FpnM_dt[tensor] = {[0] = 1.0, [4] = 1.0, [8] = 1.0};
  double Fe[tensor] = {};
  double Isym[tensor] = {};
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        FpnM_dt[idx_2(i,j)] -= Fpn[idx_2(i,k)] * M[idx_2(k,j)];
        Fe[idx_2(i,j)] += F[idx_2(i,k)] * M[idx_2(k,j)];
      }
      FpnM_dt[idx_2(i,j)] /= dt;
      Isym[idx_2(i,j)] = Fe[idx_2(i,j)] - Fe[idx_2(j,i)];
    }
  }

  /* compute the residual terms */
  memset(Rm, 0, tensor * sizeof(*Rm));
  const double coef = gdot / sqrt(2);
  for (int k = 0; k < dim; k++) {
    for (int L = 0; L < dim; L++) {
      Rw[idx_2(k,L)] = (dt * (coef * n[idx_2(k,L)] + Wp[idx_2(k,L)] - FpnM_dt[idx_2(k,L)])
                        + 2.0 * (Wp[idx_2(k,L)] + Wp[idx_2(L,k)]) );
      for (int i = 0; i < dim; i++) {
        Rm[idx_2(k,L)] += (Fpn[idx_2(i,k)] * (coef * n[idx_2(i,L)] + Wp[idx_2(i,L)] - FpnM_dt[idx_2(i,L)])
                           + 2.0 * lam * (F[idx_2(i,k)] * Fe[idx_2(i,L)] - Fe[idx_2(L,i)] * F[idx_2(i,k)]) );
      }
    }
  }
  *Rlam = 0.5 * cblas_ddot(tensor, Isym, 1, Isym, 1);
  return err;
}

int BPA_compute_DFe_DM(double * restrict DFe_DM,
                       const double * restrict F)
{
  int err = 0;
  memset(DFe_DM, 0, tensor4 * sizeof(*DFe_DM));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        /* NOTE: no loop over L idx since we have I(L,j) which reduces
           to a loop over only j */
        const int ijkl = idx_4(i,j,k,j);
        DFe_DM[ijkl] = F[idx_2(i,k)]; /* eye[idx_2(L,j)] */
      }
    }
  }
  return err;
}

int BPA_compute_Dgdot_Dtau(double *Dgdot_Dtau,
                           const double s_s,
                           const double tau)
{
  int err = 0;
  double pow_term = pow(tau/s_s,5.0/6.0);
  *Dgdot_Dtau = ((5 * param_A * s_s * param_gdot0) / (6 * param_T * tau) * pow_term
                 * exp(-param_A * s_s * (1.0 - pow_term) / param_T));
  return err;
}

int BPA_compute_Dtau_Dsig(double * restrict Dtau_Dsig,
                          const double * restrict sig)
{
  int err = 0;
  double mag = 0;
  for (int i = 0; i < tensor; i++){
    mag += sig[i] * sig[i];
  }
  mag = sqrt(2.0 * mag);

  /* corner case, sig = 0 (should never happen, produce error code) */
  if (mag == 0) {
    err++;
    mag = 1.0;
  }

  for (int i = 0; i < tensor; i++){
    Dtau_Dsig[i] = sig[i] / mag;
  }
  return err;
}

int BPA_compute_Dn_Dsig(double * restrict Dn_Dsig,
                        const double * restrict n,
                        const double tau)
{
  int err = 0;
  double coef = 1.0 / (2*sqrt(2));
  if (tau > 0) coef /= tau;
  else err++;

  memset(Dn_Dsig, 0, tensor4 * sizeof(*Dn_Dsig));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          Dn_Dsig[idx_4(i,j,k,l)] = coef * (eye[idx_2(i,k)] * eye[idx_2(l,j)]
                                            + eye[idx_2(i,l)] * eye[idx_2(k,j)]
                                            - n[idx_2(i,j)] * n[idx_2(k,l)]);
        }
      }
    }
  }
  return err;
}

int BPA_compute_Dplam_DM(double * restrict Dplam_DM,
                         const double * restrict Fp,
                         const double * restrict Cp)
{
  int err = 0;
  const double coef = -1.0 / sqrt(3 * (Cp[0] + Cp[4] + Cp[8]));
  memset(Dplam_DM, 0, tensor * sizeof(*Dplam_DM));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        Dplam_DM[idx_2(i,j)] += coef * Fp[idx_2(k,i)] * Cp[idx_2(k,j)];
      }
    }
  }
  return err;
}

int BPA_compute_DCpdev_DM(double * restrict DCpdev_DM,
                          const double * restrict Fp,
                          const double * restrict Cpdev)
{
  int err = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {

          const int ijkl = idx_4(i,j,k,l);
          DCpdev_DM[ijkl] = 2./3. * Cpdev[idx_2(i,j)] * Fp[idx_2(l,k)];

          for (int p = 0; p < dim; p++) {
            for(int q = 0; q < dim; q++) {
              DCpdev_DM[ijkl] -= (0.5 * (eye[idx_2(i,k)] *  eye[idx_2(l,j)]
                                         + eye[idx_2(i,l)] *  eye[idx_2(k,j)])
                                  * (Fp[idx_2(p,k)] * Cpdev[idx_2(l,q)] 
                                     + Fp[idx_2(p,l)] * Cpdev[idx_2(k,q)])
                                  );
            }
          }
        }
      }
    }
  }
  return err;
}

int BPA_compute_DBdev_DM(double * restrict DB_DM,
                         const double * restrict Fp)
{
  int err = 0;
  double Cp[tensor] = {};
  double Cpdev[tensor] = {};
  bpa_compute_Cp(Cp,Fp);
  bpa_compute_Cpdev(Cpdev,Cp);
  const double plam = bpa_compute_plam(Cp);

  /* compute chain rule terms */
  double der_inv_lang = 0;
  double Dplam_DM[tensor] = {};
  double DCpdev_DM[tensor4] = {};
  err += BPA_der_inverse_langevin(plam / sqrt(param_N), &der_inv_lang);
  err += BPA_compute_Dplam_DM(Dplam_DM,Fp,Cp);
  err += BPA_compute_DCpdev_DM(DCpdev_DM,Fp,Cpdev);

  /* coefficients */
  double inv_lang = 0;
  err += BPA_inverse_langevin(plam / sqrt(param_N), &inv_lang);
  const double coef1 = (-param_Cr * sqrt(param_N) * inv_lang / (3 * plam * plam)
                        + param_Cr * der_inv_lang / (3 * plam));
  const double coef2 = param_Cr * sqrt(param_N) * inv_lang / (3 * plam);

  /* compute result */
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          DB_DM[idx_4(i,j,k,l)] = (coef1 * Cpdev[idx_2(i,j)] * Dplam_DM[idx_2(k,l)] 
                                   + coef2 * DCpdev_DM[idx_4(i,j,k,l)]);
        }
      }
    }
  }

  return err;
}

int BPA_compute_DSdev_DM(double * restrict DSdev_DM,
                         const double * restrict F,
                         const double * restrict Fe,
                         const double * restrict Ce,
                         const HOMMAT *p_hmat)
{
  int err = 0;

  double FeF[tensor] = {};
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        FeF[idx_2(i,j)] += Fe[idx_2(k,i)] * F[idx_2(k,j)];
      }
    }
  }

  /* compute Ldev */
  double Ldev[tensor4] = {};
  bpa_compute_Ldev(Ce,p_hmat,Ldev);

  /* compute result */
  memset(DSdev_DM, 0, tensor4 * sizeof(*DSdev_DM));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          for (int p = 0; p < dim; p++) {
            DSdev_DM[idx_4(i,j,k,l)] += 0.5 * (Ldev[idx_4(i,j,l,p)] * FeF[idx_2(p,k)]
                                               + Ldev[idx_4(i,j,p,l)] * FeF[idx_2(p,k)]);
          }
        }
      }
    }
  }
  return err;
}

int BPA_compute_Dsig_DM(double * restrict Dsig_DM,
                        const double * restrict F,
                        const double * restrict M,
                        const HOMMAT *p_hmat)
{
  int err = 0;

  /* compute deformation tensors */
  double Fe[tensor] = {};
  double Ce[tensor] = {};
  double Fp[tensor] = {};
  bpa_compute_Fe_Ce(F,M,Fe,Ce);
  err += inv3x3(M,Fp);
  const double Je_inv = 1.0 / det3x3(Fe);

  /* compute stress tensors */
  double Sdev[tensor] = {};
  double Bdev[tensor] = {};
  bpa_compute_Sdev(Ce,p_hmat,Sdev);
  err += BPA_compute_Bdev(Bdev,Fp);

  /* compute chain rule tensors */
  double DSdev_DM[tensor4] = {};
  double DBdev_DM[tensor4] = {};
  err += BPA_compute_DSdev_DM(DSdev_DM,F,Fe,Ce,p_hmat);
  err += BPA_compute_DBdev_DM(DBdev_DM,Fp);

  /* compute result */
  memset(Dsig_DM, 0, tensor4 * sizeof(*Dsig_DM));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          const int ijkl = idx_4(i,j,k,l);
          for (int p = 0; p < dim; p++) {
            Dsig_DM[ijkl] += (Je_inv * (Sdev[idx_2(l,p)] - Bdev[idx_2(l,p)])
                              *(F[idx_2(i,k)] * Fe[idx_2(j,p)] + Fe[idx_2(i,p)] * F[idx_2(j,k)]));
            for (int q = 0; q < dim; q++) {
              const int pqkl = idx_4(p,q,k,l);
              Dsig_DM[ijkl] += (Je_inv
                                * (DSdev_DM[pqkl] - DBdev_DM[pqkl]
                                   - Fp[idx_2(l,k)] * (Sdev[idx_2(l,p)] - Bdev[idx_2(l,q)]))
                                * Fe[idx_2(i,p)] * F[idx_2(j,q)]);
            }
          }
        }
      }
    }
  }
  return err;
}

int BPA_compute_Dn_DM(double * restrict Dn_DM,
                      const double * restrict Dsig_DM,
                      const double * restrict n,
                      const double tau)
{
  int err = 0;

  double Dn_Dsig[tensor4] = {};
  err += BPA_compute_Dn_Dsig(Dn_Dsig,n,tau);

  memset(Dn_DM, 0, tensor4 * sizeof(*Dn_DM));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          const int ijkl = idx_4(i,j,k,l);
          for (int p = 0; p < dim; p++) {
            for (int q = 0; q < dim; q++) {
              Dn_DM[ijkl] += Dn_Dsig[idx_4(i,j,p,q)] * Dsig_DM[idx_4(p,q,k,l)];
            }
          }
        }
      }
    }
  }

  return err;
}

int BPA_compute_Dgdot_DM(double *restrict Dgdot_DM,
                         const double * restrict Dsig_DM,
                         const double * restrict sig,
                         const double tau,
                         const double s_s)
{
  int err = 0;
  double Dgdot_Dtau = 0;
  double Dtau_Dsig[tensor] = {};
  err += BPA_compute_Dgdot_Dtau(&Dgdot_Dtau,s_s,tau);
  err += BPA_compute_Dtau_Dsig(Dtau_Dsig,sig);

  memset(Dgdot_DM, 0, tensor * sizeof(*Dgdot_DM));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          Dgdot_DM[idx_2(k,l)] += Dgdot_Dtau * Dtau_Dsig[idx_2(i,j)] * Dsig_DM[idx_4(i,j,k,l)];
        }
      }
    }
  }

  return err;
}

int BPA_int_alg_tan_terms(double * restrict DM_M,
                          double * restrict DM_W,
                          double * restrict DM_lam,
                          double * restrict DW_M,
                          double * restrict DW_W,
                          const double dt,
                          const double gdot,
                          const double lam,
                          const double tau,
                          const double s_s,
                          const double * restrict Fpn,
                          const double * restrict n,
                          const double * restrict sig,
                          const double * restrict F,
                          const double * restrict M,
                          const double * restrict Wp,
                          const HOMMAT *p_hmat)
{
  int err = 0;
  double Dsig_DM[tensor4] = {};
  double Dgdot_DM[tensor] = {};
  double Dn_DM[tensor4] = {};
  double Fe[tensor] = {};
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
     Fe[idx_2(i,j)] = 0.0;
      for (int k = 0; k < dim; k++) {
        Fe[idx_2(i,j)] += F[idx_2(i,k)] * M[idx_2(k,j)];
      }
    }
  }

  err += BPA_compute_Dsig_DM(Dsig_DM,F,M,p_hmat);
  err += BPA_compute_Dgdot_DM(Dgdot_DM,Dsig_DM,sig,tau,s_s);
  err += BPA_compute_Dn_DM(Dn_DM,Dsig_DM,n,tau);

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      DM_lam[idx_2(i,j)] = 0.0;
      for (int k = 0; k < dim; k++) {
        DM_lam[idx_2(i,j)] += 2.0 * (F[idx_2(k,i)] * Fe[idx_2(k,j)] - Fe[idx_2(j,k)] * F[idx_2(k,i)]);
        for (int L = 0; L < dim; L++) {
          const int ijkl = idx_4(i,j,k,L);
          DM_M[ijkl] = 0.0;
          for (int p = 0; p < dim; p++) {
            DM_M[ijkl] += (Fpn[idx_2(p,i)] * (Dgdot_DM[idx_2(k,L)] * n[idx_2(p,j)]
                                              + gdot * Dn_DM[idx_4(p,j,k,L)]
                                              + Fpn[idx_2(p,k)] * eye[idx_2(L,j)])
                           + 2.0 * lam * (F[idx_2(p,i)] * F[idx_2(p,k)] * eye[idx_2(L,j)]
                                          - F[idx_2(j,k)] * F[idx_2(L,i)]) );
          }
          DM_W[ijkl] = Fpn[idx_2(k,i)] * eye[idx_2(L,j)];

          DW_M[ijkl] = dt * (Dgdot_DM[idx_2(k,L)] * n[idx_2(i,j)]
                             + gdot * Dn_DM[ijkl] + Fpn[idx_2(i,k)] * eye[idx_2(L,i)] / dt);
          DW_W[ijkl] = (dt + 2.0) * eye[idx_2(i,k)] * eye[idx_2(L,j)] + 2.0 * eye[idx_2(j,k)] * eye[idx_2(L,i)];

        }
      }
    }
  }

  return err;
}

int BPA_int_alg_tan(double * restrict TAN,
                    const double dt,
                    const double gdot,
                    const double lam,
                    const double tau,
                    const double s_s,
                    const double *Fpn,
                    const double *n,
                    const double *sig,
                    const double *F,
                    const double *M,
                    const double *Wp,
                    const HOMMAT *p_hmat)
{
  int err = 0;

  /* compute terms */
  double DM_M[tensor4] = {};
  double DM_W[tensor4] = {};
  double DM_lam[tensor] = {};
  double DW_M[tensor4] = {};
  double DW_W[tensor4] = {};
  err += BPA_int_alg_tan_terms(DM_M, DM_W, DM_lam, DW_M, DW_W,  /*< out */
                               dt, gdot, lam, tau, s_s, Fpn, n, /*< in (and next) */
                               sig, F, M, Wp, p_hmat);

  /* Assemble tangent blocks, TAN [19x19]
     TAN = T-----------------T
           | DM_M | DM_W | DM|
           |      |      | _l|
           |------|------|---|
           | DW_M | DW_W | 0 |
           |      |      |   |
           |------|------|---|
           | DM_l |   0  | 0 |
           L-----------------J
   */
  const int row2_start = tensor * tangent;
  const int row3_start = 2 * row2_start;
  for (int i = 0; i < tensor; i++) {
    const int blk1_start = i * tangent;
    const int blk2_start = blk1_start + tensor;
    const int blk3_start = blk2_start + tensor;

    /* 1st row */
    memcpy(TAN + blk1_start, DM_M + i * tensor, tensor * sizeof(*TAN));
    memcpy(TAN + blk2_start, DM_W + i * tensor, tensor * sizeof(*TAN));
    TAN[blk3_start] = DM_lam[i];

    /* 2nd row */
    memcpy(TAN + row2_start + blk1_start, DW_M + i * tensor, tensor * sizeof(*TAN));
    memcpy(TAN + row2_start + blk2_start, DW_W + i * tensor, tensor * sizeof(*TAN));

    /* 3rd row */
    TAN[row3_start + i] = DM_lam[i];
  }

  return err;
}

int BPA_int_alg_res(double *RES,
                    const double dt,
                    const double gdot,
                    const double lam,
                    const double * restrict Fpn,
                    const double * restrict n,
                    const double * restrict F,
                    const double * restrict M,
                    const double * restrict Wp)
{
  int err = 0;
  /*
   * OPTIMIZATION NOTE: After debugging/error-checking with separate
   * arrays, can pass in appropriately indexed pointers to locations
   * in RES to remove additional assembly/memory overhead.
   */
  double Rm[tensor] = {};
  double Rw[tensor] = {};
  double Rlam;
  err += BPA_int_alg_res_terms(Rm, Rw, &Rlam, dt, gdot, lam, Fpn, n, F, M, Wp);
  memcpy(RES, Rm, tensor * sizeof(*RES));
  memcpy(RES + tensor, Rw, tensor * sizeof(*RES));
  RES[tangent - 1] = Rlam;

  /* OPTIMIZED */
  /* err += BPA_int_alg_res_terms(RES, RES + tensor, RES + tangent - 1, */
  /*                              dt,gdot,lam,Fpn,F,M,Wp); */
  return err;
}

int BPA_dev_stress(const Constitutive_model *m,
                   const void *ctx,
                   Matrix_double *dev_stress)
{
  int err = 0;
  const BPA_ctx *CTX = ctx;
  double Fe[tensor] = {}, Ce[tensor] = {};
  bpa_compute_Fe_Ce(CTX->F,m->vars.Fs[_M].m_pdata,Fe,Ce);
  bpa_compute_Sdev(Ce,m->param->p_hmat,dev_stress->m_pdata);
  return err;
}

int BPA_dudj(const Constitutive_model *m,
             const void *ctx,
             double *dudj)
{
  int err = 0;
  const BPA_ctx *CTX = ctx;
  double Fe[tensor] = {};
  bpa_compute_Fe(CTX->F,m->vars.Fs[_M].m_pdata,Fe);
  const double Je = det3x3(Fe);
  bpa_compute_dudj(Je,m->param->p_hmat,dudj);
  return err;
}

int BPA_dev_tangent(const Constitutive_model *m,
                    const void *ctx,
                    Matrix_double *dev_tangent)
{
  int err = 0;
  const BPA_ctx *CTX = ctx;
  double Fe[tensor] = {}, Ce[tensor] = {};
  bpa_compute_Fe_Ce(CTX->F,m->vars.Fs[_M].m_pdata,Fe,Ce);
  bpa_compute_Ldev(Ce,m->param->p_hmat,dev_tangent->m_pdata);
  return err;
}

int BPA_d2udj2(const Constitutive_model *m,
               const void *ctx,
               double *d2udj2)
{
  int err = 0;
  const BPA_ctx *CTX = ctx;
  double Fe[tensor] = {};
  bpa_compute_Fe(CTX->F,m->vars.Fs[_M].m_pdata,Fe);
  const double Je = det3x3(Fe);
  bpa_compute_d2udj2(Je,m->param->p_hmat,d2udj2);
  return err;
}

int BPA_update_vars(Constitutive_model *m)
{
  int err = 0;
  Matrix_copy(m->vars.Fs[_M_n],m->vars.Fs[_M]);
  Matrix_copy(m->vars.Fs[_W_n],m->vars.Fs[_W]);
  Matrix_copy(m->vars.Fs[_Fe_n],m->vars.Fs[_Fe]);

  /* alias */
  Vector_double *vars = m->vars.state_vars;
  Vec_v(*vars,_s_n + 1) = Vec_v(*vars,_s + 1);
  Vec_v(*vars,_lam_n + 1) = Vec_v(*vars,_lam + 1);
  return err;
}

int BPA_reset_vars(Constitutive_model *m)
{
  int err = 0;
  Matrix_copy(m->vars.Fs[_M],m->vars.Fs[_M_n]);
  Matrix_copy(m->vars.Fs[_W],m->vars.Fs[_W_n]);
  Matrix_copy(m->vars.Fs[_Fe],m->vars.Fs[_Fe_n]);

  /* alias */
  Vector_double *vars = m->vars.state_vars;
  Vec_v(*vars,_s + 1) = Vec_v(*vars,_s_n + 1);
  Vec_v(*vars,_lam + 1) = Vec_v(*vars,_lam_n + 1);
  return err;
}

int BPA_model_info(Model_var_info **info)
{
  int err = 0;

  /* make sure I don't leak memory */
  if (*info != NULL) err += model_var_info_destroy(info);

  /* allocate pointers */
  (*info) = malloc(sizeof(**info));
  (*info)->n_Fs = _n_Fs;
  (*info)->n_vars = _n_vars;
  (*info)->F_names = malloc(_n_Fs * sizeof( ((*info)->F_names) ));
  (*info)->var_names = malloc( _n_vars * sizeof( ((*info)->var_names) ));

  /* allocate/copy strings */
  (*info)->F_names[_M] = strdup("M");
  (*info)->F_names[_M_n] = strdup("M_n");
  (*info)->F_names[_W] = strdup("Wp");
  (*info)->F_names[_W_n] = strdup("Wp_n");
  (*info)->F_names[_Fe] = strdup("Fe");
  (*info)->F_names[_Fe_n] = strdup("Fe_n");
  (*info)->var_names[_s_n] = strdup("s_n");
  (*info)->var_names[_s] = strdup("s");
  (*info)->var_names[_lam_n] = strdup("lam_n");
  (*info)->var_names[_lam] = strdup("lam");

  return err;
}

int BPA_inverse_langevin(const double y,
                         double *inv_lang)
{
  /*
   * REFERENCE:
   *   Radoslaw Jedynak, Rheol Acta (2015) 54:29–39, DOI 10.1007/s00397-014-0802-2
   *
   * See rounded NPA[3/2] formulation, eq. 21
   */
  int err = 0;
  if (y < 0 || y > 1) err++;
  *inv_lang = y * (3.0 - 2.6 * y + 0.7 * y * y) / (1.0 - 0.9 * y - 0.1 * y * y);
  return err;
}

int BPA_der_inverse_langevin(const double y,
                             double  *der_inv_lang)
{
 /*
   * REFERENCE:
   *   Radoslaw Jedynak, Rheol Acta (2015) 54:29–39, DOI 10.1007/s00397-014-0802-2
   *
   * See rounded NPA[3/2] formulation, eq. 30
   */
  int err = 0;
  if (y < 0 || y > 1) err++;
  *der_inv_lang = -7.0 + 1.0 / ((1.0 - y) * (1.0 - y)) + 900.0 / ((y + 10.0) * (y + 10.0));
  return err;
}

int BPA_compute_gdot(double *gdot,
                     const double tau,
                     const double s_s)
{
  int err = 0;
  /* NOTE: Makes use of global variables for material properties!! */
  *gdot = param_gdot0 * exp( - param_A * s_s / param_T
                             * (1.0 - pow(tau / s_s, 5.0 / 6.0)) );
  return err;
}

int BPA_compute_sdot(double *sdot,
                     const double s,
                     const double gdot)
{
  /* NOTE: makes use of (file) global parameters! */
  int err = 0;
  *sdot = param_h * (1.0 - s / param_s_ss) * gdot;
  return err;
}

int BPA_compute_Bdev(double *Bdev,
                     const double *Fp)
{
  /* NOTE: makes use of (file) global parameters! */
  int err = 0;

  /* compute the deviatoric plastic deformation */
  double Cp[tensor] = {};
  double Cpdev[tensor] = {};
  bpa_compute_Cp(Cp,Fp);
  bpa_compute_Cpdev(Cpdev,Cp);

  /* compute the backstess coefficient */
  const double lam_p = bpa_compute_plam(Cp);
  double coeff = 0;
  err += BPA_inverse_langevin(lam_p / sqrt(param_N), &coeff);
  coeff *= param_Cr * sqrt(param_N) / (3 * lam_p);

  /* compute the deviatoric back stress tensor */
  for (int i = 0; i < tensor; i++) {
    Bdev[i] = coeff * Cpdev[i];
  }

  return err;
}

int BPA_compute_loading_dir(double * restrict normal,
                            double * restrict eq_sig_dev,
                            double * restrict tau,
                            const double * restrict Sdev,
                            const double * restrict Bdev,
                            const double * restrict Fe)
{
  int err = 0;

  /* clear the result buffer(s) */
  memset(eq_sig_dev, 0, tensor * sizeof(*eq_sig_dev));
  memset(normal, 0, tensor * sizeof(*normal));

  /* compute temporary buffers */
  const double inv_Je = 1.0 / det3x3(Fe);
  double tmp[tensor] = {}, tmp2[tensor] = {};
  for (int i = 0; i < tensor; i++) {
    tmp[i] = inv_Je * (Sdev[i] - Bdev[i]);
    tmp2[i] = 0;
  }

  /* compute the intermediate result */
  for (int i = 0; i < dim; i++) {
    for(int j = 0; j < dim; j++) {
      for(int k = 0; k < dim; k++) {
        tmp2[idx_2(i,j)] += Fe[idx_2(i,k)] * tmp[idx_2(k,j)];
      }
    }
  }

  /* comptue the final result */
  for (int i = 0; i < dim; i++) {
    for(int j = 0; j < dim; j++) {
      for(int k = 0; k < dim; k++) {
        eq_sig_dev[idx_2(i,j)] += tmp2[idx_2(i,k)] * Fe[idx_2(j,k)];
      }
    }
  }

  /* compute equivalent plastic stress and the normal to the loading
     direction */
  *tau = sqrt(0.5 * cblas_ddot(tensor, eq_sig_dev, 1, eq_sig_dev, 1) );
  cblas_daxpy(tensor,1.0 / (sqrt(2.0) * (*tau)), eq_sig_dev, 1, normal, 1);

  return err;
}

/*
 * Public interface for the BPA model
 */
int plasticity_model_BPA_initialize(Model_parameters *p)
{
  int err = 0;
  p->integration_algorithm = BPA_int_alg;
  p->compute_dev_stress = BPA_dev_stress;
  p->compute_dudj = BPA_dudj;
  p->compute_dev_tangent = BPA_dev_tangent;
  p->compute_d2udj2 = BPA_d2udj2;
  p->update_state_vars = BPA_update_vars;
  p->reset_state_vars = BPA_reset_vars;
  p->get_var_info = BPA_model_info;
  return err;
}

int plasticity_model_BPA_set_initial_values(Constitutive_model *m)
{
  int err = 0;
  /* Algorithmic Wp is 0 at start */
  memset(m->vars.Fs[_W].m_pdata, 0, tensor * sizeof(double));
  memset(m->vars.Fs[_W_n].m_pdata, 0, tensor * sizeof(double));

  /* s is s0 at start */
  m->vars.state_vars->m_pdata[_s] = param_s0;
  m->vars.state_vars->m_pdata[_s_n] = param_s0;

  /* lam is 0 at start */
  m->vars.state_vars->m_pdata[_lam] = 0;
  m->vars.state_vars->m_pdata[_lam_n] = 0;

  return err;
}

int plasticity_model_BPA_ctx_build(void **ctx,
                                   const double *F,
                                   const double dt)
{
  int err = 0;
  BPA_ctx *t_ctx = malloc(sizeof(*t_ctx));
  *ctx = t_ctx;
  t_ctx->dt = dt;
  memcpy(t_ctx->F, F, tensor * sizeof(*F));
  return err;
}

int plasticity_model_BPA_ctx_destroy(void **ctx)
{
  int err = 0;
  free(*ctx);
  *ctx = NULL;
  return err;
}
