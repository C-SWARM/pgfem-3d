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
#include "allocation.h"
#include "constitutive_model.h"
#include "cm_placeholder_functions.h"
#include "mkl_cblas.h"
#include "mkl_lapack.h"
#include "index_macros.h"
#include "new_potentials.h"
#include "state_variables.h"
#include "utils.h"
#include <assert.h>

/* Define constant dimensions. Note cannot use `static const` with
   initialization list */
namespace {
constexpr int     dim = 3;
constexpr int  tensor = 9;
constexpr int tensor4 = 81;
constexpr int tan_row = 10;
constexpr int tan_col = 10;

/* Set to value > 0 for extra diagnostics/printing */
constexpr int BPA_PRINT_LEVEL = 0;

/* constants/enums */
constexpr int    _n_Fs = 6;
constexpr int  _n_vars = 4;
constexpr int _n_flags = 0;
enum {_Fe,_Fp,_F,_Fe_n,_Fp_n,_F_n};
enum {_s,_lam,_s_n,_lam_n};
constexpr double eye[tensor] = {1.0,0,0, 0,1.0,0, 0,0,1.0};

/* enumerations for indexing into the list of model parameters */
enum {mcA, mcAlpha, mcCr, mcGdot0, mcH, mcN, mcS0, mcSss, mcT, N_PARAM};
}

static double bpa_compute_bulk_mod(const HOMMAT *mat) {
  return ( (2* mat->G * (1 + mat->nu)) / (3 * (1 - 2 * mat->nu)) );
}

static void bpa_compute_Ce(double *Ce, const double *Fe) {
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              dim,dim,dim,1.0,Fe,dim,Fe,dim,
              0.0,Ce,dim);
}

static void bpa_compute_Cp(double * __restrict Cp, const double * __restrict Fp) {
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
              dim,dim,dim,1.0,Fp,dim,Fp,dim,
              0.0,Cp,dim);
}

static double bpa_compute_plam(const double *Cp) {
  return sqrt((Cp[0] + Cp[4] + Cp[8]) / 3.0);
}

static void bpa_compute_Cpdev(double * __restrict Cpdev,
                              const double * __restrict Cp)
{
  const double lam = (Cp[0] + Cp[4] + Cp[8]) / 3.0;
  memcpy(Cpdev, Cp, tensor * sizeof(*Cp));
  Cpdev[0] -= lam;
  Cpdev[4] -= lam;
  Cpdev[8] -= lam;
}

static void bpa_compute_Fp(double * __restrict Fp, const double * __restrict F,
                           const double * __restrict Fe)
{
  double invFe[tensor] = {};
  inv3x3(Fe,invFe);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              dim,dim,dim,1.0,invFe,dim,F,dim,
              0.0,Fp,dim);
}

static void bpa_compute_Sdev(const double *Ce, const HOMMAT *p_hmat,
                             double *Sdev)
{
  devStressFuncPtr Stress = getDevStressFunc(-1,p_hmat);
  Stress(Ce,p_hmat,Sdev);
}

static void bpa_compute_Ldev(const double *Ce, const HOMMAT *p_hmat,
                             double *Ldev)
{
  matStiffFuncPtr Tangent = getMatStiffFunc(-1,p_hmat);
  Tangent(Ce,p_hmat,Ldev);
}

static void bpa_compute_dudj(const double Je, const HOMMAT *p_hmat,
                             double *dudj)
{
  dUdJFuncPtr Pressure = getDUdJFunc(-1,p_hmat);
  Pressure(Je,p_hmat,dudj);
}

static void bpa_compute_d2udj2(const double Je, const HOMMAT *p_hmat,
                               double *d2udj2)
{
  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,p_hmat);
  D_Pressure(Je,p_hmat,d2udj2);
}

int bpa_compute_loading_dir(double * __restrict normal,
                            double * __restrict eq_sig_dev,
                            double * __restrict tau,
                            const double * __restrict Sdev,
                            const double * __restrict Bdev,
                            const double * __restrict Fe)
{
  int err = 0;

  /* compute temporary buffers */
  const double Je = det3x3(Fe);
  double SmB[tensor] = {};
  double tmp[tensor] = {};
  for (int i = 0; i < tensor; i++) {
    SmB[i] = Sdev[i] - Bdev[i];
  }

  /* compute the intermediate result */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dim, dim, dim, 1.0 / Je, Fe, dim, SmB, dim,
              0.0, tmp, dim);

  /* comptue eq. stress */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
              dim, dim, dim, 1.0, tmp, dim, Fe, dim,
              0.0, eq_sig_dev, dim);

  /* compute eff. stress */
  *tau = sqrt(0.5 * cblas_ddot(tensor, eq_sig_dev, 1, eq_sig_dev, 1) );

  /* compute the loading direction (normal) */
  memset(normal, 0, tensor * sizeof(*normal));
  if (*tau != 0) {
    cblas_daxpy(tensor, 1.0 / (sqrt(2.0) * (*tau)), eq_sig_dev, 1, normal, 1);
  }

  return err;
}

static int bpa_compute_gdot(double *gdot,
                            const double param_gdot0,
                            const double param_A,
                            const double param_T,
                            const double tau,
                            const double s_s)
{
  int err = 0;
  *gdot = param_gdot0 * exp( - param_A * s_s / param_T
                             * (1.0 - pow(tau / s_s, 5.0 / 6.0)) );
  return err;
}

static int bpa_compute_res_vec(double * __restrict RES,
                               const double dt,
                               const double gdot,
                               const double lam,
                               const double Jp,
                               const double * __restrict n,
                               const double * __restrict Mn,
                               const double * __restrict Wp,
                               const double * __restrict F,
                               const double * __restrict Fe)
{
  int err = 0;

  /* compute intermediate terms */
  double invFe[tensor] = {};
  double FMn[tensor] = {};
  double Lp[tensor] = {};
  double dtFMnLp[tensor] = {};
  err += inv3x3(Fe,invFe);
  for (int i = 0; i < tensor; i++) {
    Lp[i] = gdot * n[i] + Wp[i];
  }
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dim, dim, dim, 1.0, F, dim, Mn, dim,
              0.0, FMn, dim);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dim, dim, dim, dt, FMn, dim, Lp, dim,
              0.0, dtFMnLp, dim);

  /* compute residual vector */
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      const int ij = idx_2(i,j);
      const int ji = idx_2(j,i);
      RES[ij] = -(2 * Fe[ij] - Fe[ji] - FMn[ij] + dtFMnLp[ij]
                 - lam * Jp * invFe[ji]);
    }
  }
  RES[tensor] = -(Jp - 1.0);

  return err;
}

int bpa_inverse_langevin(const double y,
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

int bpa_der_inverse_langevin(const double y,
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

static int bpa_compute_Bdev(double *Bdev,
                            const double param_N,
                            const double param_Cr,
                            const double *Fp)
{
  int err = 0;

  /* compute the deviatoric plastic deformation */
  double Cp[tensor] = {};
  double Cpdev[tensor] = {};
  bpa_compute_Cp(Cp,Fp);
  bpa_compute_Cpdev(Cpdev,Cp);

  /* compute the backstess coefficient */
  const double lam_p = bpa_compute_plam(Cp);
  double coeff = 0;
  err += bpa_inverse_langevin(lam_p / sqrt(param_N), &coeff);
  coeff *= param_Cr * sqrt(param_N) / (3.0 * lam_p);

  /* compute the deviatoric back stress tensor */
  for (int i = 0; i < tensor; i++) {
    Bdev[i] = coeff * Cpdev[i];
  }

  return err;
}

static int bpa_compute_DBdev_DFp(double * __restrict DB_DFp,
                                 const double * __restrict Fp,
                                 const double param_N,
                                 const double param_Cr)
{
  int err = 0;
  double Cp[tensor] = {};
  double Cpdev[tensor] = {};
  double invFp[tensor] = {};
  bpa_compute_Cp(Cp,Fp);
  bpa_compute_Cpdev(Cpdev,Cp);
  err += inv3x3(Fp,invFp);
  const double plam = bpa_compute_plam(Cp);
  const double inv_lang_arg = plam / sqrt(param_N);
  // const double Jp23 = pow(det3x3(Fp),-2.0 / 3.0);
  double inv_lang = 0;
  double inv_lang_p = 0;
  err += bpa_inverse_langevin(inv_lang_arg, &inv_lang);
  err += bpa_der_inverse_langevin(inv_lang_arg, &inv_lang_p);
  const double coef_1 = param_Cr * sqrt(param_N) * inv_lang / (3.0 * plam);
  const double coef_2 = ((param_Cr * inv_lang_p / (3.0 * plam)
                         - param_Cr * sqrt(param_N) / (3.0 * plam * plam))
                         / sqrt(3.0 * (Cp[0] + Cp[4] + Cp[8])) );

  /* compute result */
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          DB_DFp[idx_4(i,j,k,l)] = (coef_1 * (eye[idx_2(i,k)] * Fp[idx_2(j,l)]
                                              + Fp[idx_2(i,l)] * eye[idx_2(j,k)]
                                              -2./3. * eye[idx_2(i,j)] * Fp[idx_2(k,l)])
                                    + coef_2 * Cpdev[idx_2(i,j)] * Fp[idx_2(k,l)]);
        }
      }
    }
  }

  return err;
}

static int bpa_compute_DSdev_DFe(double * __restrict DSdev_DFe,
                                 const double * __restrict Ce,
                                 const double * __restrict Fe,
                                 const HOMMAT *p_hmat)
{
  int err = 0;
  double Ldev[tensor4] = {};
  bpa_compute_Ldev(Ce,p_hmat,Ldev);

  memset(DSdev_DFe, 0, tensor4 * sizeof(*DSdev_DFe));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          for (int r = 0; r < dim; r++) {
            DSdev_DFe[idx_4(i,j,k,l)] += 0.5 * (Ldev[idx_4(i,j,l,r)] * Fe[idx_2(k,r)]
                                                + Ldev[idx_4(i,j,r,l)] * Fe[idx_2(k,r)]);
          }
        }
      }
    }
  }

  return err;
}

static int bpa_compute_DFp_DFe(double * __restrict DFp_DFe,
                               const double * __restrict invFe,
                               const double * __restrict Fp)
{
  int err = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          DFp_DFe[idx_4(i,j,k,l)] = -invFe[idx_2(i,k)] * Fp[idx_2(l,j)];
        }
      }
    }
  }
  return err;
}

static int bpa_compute_Dsig_DFe(double * __restrict Dsig_DFe,
                                const double * __restrict F,
                                const double * __restrict Fe,
                                const double * __restrict Fp,
                                const HOMMAT *p_hmat,
                                const double param_N,
                                const double param_Cr)
{
  int err = 0;

  /* compute deformation */
  const double inv_Je = 1.0 / det3x3(Fe);
  double invFe[tensor] = {};
  err += inv3x3(Fe,invFe);
  double Ce[tensor] = {};
  bpa_compute_Ce(Ce,Fe);

  /* compute stress tensors */
  double Sdev[tensor] = {};
  double Bdev[tensor] = {};
  bpa_compute_Sdev(Ce,p_hmat,Sdev);
  err += bpa_compute_Bdev(Bdev, param_N, param_Cr, Fp);

  /* compute derivatives */
  double DB_DFp[tensor4] = {};
  double DSdev_DFe[tensor4] = {};
  double DFp_DFe[tensor4] = {};
  err += bpa_compute_DBdev_DFp(DB_DFp,Fp,param_N, param_Cr);
  err += bpa_compute_DSdev_DFe(DSdev_DFe,Ce,Fe,p_hmat);
  err += bpa_compute_DFp_DFe(DFp_DFe,invFe,Fp);

  /* compute other terms */
  double SmB[tensor] = {};
  for(int i = 0; i < tensor; i++){
    SmB[i] = Sdev[i] - Bdev[i];
  }

  memset(Dsig_DFe, 0, tensor4 * sizeof(*Dsig_DFe));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          for (int p = 0; p < dim; p++) {
            for (int q = 0; q < dim; q++) {
              Dsig_DFe[idx_4(i,j,k,l)] += ((-inv_Je * invFe[idx_2(l,k)] * Fe[idx_2(i,p)] * SmB[idx_2(p,q)] * Fe[idx_2(j,q)])
                                           + (inv_Je * eye[idx_2(i,k)] * eye[idx_2(l,p)] * SmB[idx_2(p,q)] * Fe[idx_2(j,q)])
                                           + (inv_Je * Fe[idx_2(i,p)] * SmB[idx_2(p,q)] * eye[idx_2(j,k)] * eye[idx_2(l,q)])
                                           + (inv_Je * Fe[idx_2(i,p)] * DSdev_DFe[idx_4(p,q,k,l)] * Fe[idx_2(j,q)]));
              for (int r = 0; r < dim; r++) {
                for (int s = 0; s < dim; s++) {
                  Dsig_DFe[idx_4(i,j,k,l)] += - (inv_Je * Fe[idx_2(i,p)] * DB_DFp[idx_4(p,q,r,s)] * DFp_DFe[idx_4(r,s,k,l)] * Fe[idx_2(j,q)]);
                }
              }
            }
          }
        }
      }
    }
  }

  return err;
}

static int bpa_compute_Dgdot_Dtau(double *Dgdot_Dtau,
                                  const double param_gdot0,
                                  const double param_A,
                                  const double param_T,
                                  const double s_s,
                                  const double tau)
{
  int err = 0;
  const double pow_term = pow(tau / s_s, 1.0 / 6.0);
  double gdot = 0;
  err += bpa_compute_gdot(&gdot, param_gdot0, param_A, param_T, tau, s_s);
  *Dgdot_Dtau = 5.0 * param_A * gdot / (6.0 * param_T * pow_term);
  return err;
}

static int bpa_compute_Dgdot_Ds_s(double *Dgdot_Ds,
                                  const double param_gdot0,
                                  const double param_A,
                                  const double param_T,
                                  const double tau,
                                  const double s_s)
{
  int err = 0;
  const double pow_term = pow(tau / s_s, 5. / 6.);
  double gdot = 0.0;
  err += bpa_compute_gdot(&gdot,param_gdot0,param_A,param_T,tau,s_s);
  *Dgdot_Ds = gdot * param_A * (pow_term - 6.) / (6. * param_T);
  return err;
}

static int bpa_compute_Dgdot_DFe(double * __restrict Dgdot_DFe,
                                 const double param_gdot0,
                                 const double param_A,
                                 const double param_T,
                                 const double param_alpha,
                                 const double s_s,
                                 const double tau,
                                 const double * __restrict Dsig_DFe,
                                 const double * __restrict normal,
                                 const double * __restrict Dp_DFe)
{
  int err = 0;
  double Dgdot_Dtau = 0;
  if(tau != 0){
    err += bpa_compute_Dgdot_Dtau(&Dgdot_Dtau, param_gdot0, param_A, param_T, s_s, tau);
  }
  const double rt2 = sqrt(2.0);
  double Dgdot_Dp = 0;
  err += bpa_compute_Dgdot_Ds_s(&Dgdot_Dp, param_gdot0, param_A, param_T, tau, s_s);
  Dgdot_Dp *= -param_alpha;

  memset(Dgdot_DFe, 0, tensor * sizeof(*Dgdot_DFe));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      Dgdot_DFe[idx_2(i,j)] += Dgdot_Dp * Dp_DFe[idx_2(i,j)];
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          Dgdot_DFe[idx_2(i,j)] += Dgdot_Dtau / rt2 * normal[idx_2(k,l)] * Dsig_DFe[idx_4(k,l,i,j)];
        }
      }
    }
  }

  return err;
}

int bpa_compute_Dn_Dsig(double * __restrict Dn_Dsig,
                        const double * __restrict n,
                        const double tau)
{
  int err = 0;
  const double coef = 1.0 / (sqrt(2) * tau);
  if( !isfinite(coef) ) {
    memset(Dn_Dsig, 0, tensor4 * sizeof(*Dn_Dsig));
    return err;
  }

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          Dn_Dsig[idx_4(i,j,k,l)] = coef * (0.5 * (eye[idx_2(i,k)] * eye[idx_2(l,j)]
                                                   + eye[idx_2(i,l)] * eye[idx_2(k,j)])
                                            - n[idx_2(i,j)] * n[idx_2(k,l)]);
        }
      }
    }
  }
  return err;
}

static int bpa_compute_Dn_DFe(double * __restrict Dn_DFe,
                              const double tau,
                              const double * __restrict Dsig_DFe,
                              const double * __restrict n)
{
  int err = 0;
  double Dn_Dsig[tensor4] = {};
  err += bpa_compute_Dn_Dsig(Dn_Dsig,n,tau);
  memset(Dn_DFe, 0, tensor4 * sizeof(*Dn_DFe));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          for (int r = 0; r < dim; r++) {
            for (int s = 0; s < dim; s++) {
              Dn_DFe[idx_4(i,j,k,l)] += Dn_Dsig[idx_4(i,j,r,s)] * Dsig_DFe[idx_4(r,s,k,l)];
            }
          }
        }
      }
    }
  }

  return err;
}

static int bpa_compute_Dp_DFe(double * __restrict Dp_DFe,
                              const double * __restrict Fe,
                              const HOMMAT *p_hmat)
{
  int err = 0;
  const double kappa = bpa_compute_bulk_mod(p_hmat);
  const double Je = det3x3(Fe);
  double invFe[tensor] = {};
  err += inv3x3(Fe,invFe);

  double d2udj2 = 0;
  bpa_compute_d2udj2(Je,p_hmat,&d2udj2);

  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      Dp_DFe[idx_2(i,j)] = 0.5 * kappa * d2udj2 * Je * invFe[idx_2(j,i)];
    }
  }
  return err;
}

static int bpa_compute_tan_Fe_Fe(double * __restrict tan,
                                 const double param_gdot0,
                                 const double param_A,
                                 const double param_T,
                                 const double param_N,
                                 const double param_Cr,
                                 const double param_alpha,
                                 const double dt,
                                 const double gdot,
                                 const double tau,
                                 const double s_s,
                                 const double lam,
                                 const double * __restrict eff_sig,
                                 const double * __restrict F,
                                 const double * __restrict n,
                                 const double * __restrict Mn,
                                 const double * __restrict Fe,
                                 const double * __restrict Fp,
                                 const HOMMAT *p_hmat)
{
  int err = 0;
  double invFe[tensor] = {};
  double FMn[tensor] = {};
  double Deff_sig_DFe[tensor4] = {};
  double Dp_DFe[tensor] = {};
  double Dgdot_DFe[tensor] = {};
  double Dn_DFe[tensor4] = {};

  const double Jp = det3x3(Fp);
  err += inv3x3(Fe,invFe);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              dim, dim, dim, 1.0, F, dim, Mn, dim,
              0.0, FMn, dim);
  err += bpa_compute_Dsig_DFe(Deff_sig_DFe, F, Fe, Fp, p_hmat, param_N, param_Cr);
  err += bpa_compute_Dp_DFe(Dp_DFe,Fe,p_hmat);
  err += bpa_compute_Dgdot_DFe(Dgdot_DFe, param_gdot0, param_A, param_T,
                               param_alpha, s_s, tau, Deff_sig_DFe, n, Dp_DFe);
  err += bpa_compute_Dn_DFe(Dn_DFe,tau,Deff_sig_DFe,n);

  memset(tan, 0, tensor4 * sizeof(*tan));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          const int ijkl = idx_4(i,j,k,l);
          tan[ijkl] += (2 * eye[idx_2(i,k)] * eye[idx_2(l,j)]
                        - eye[idx_2(j,k)] * eye[idx_2(l,i)]
                        + lam * Jp * (invFe[idx_2(j,i)] * invFe[idx_2(l,k)]
                                      + invFe[idx_2(j,k)] * invFe[idx_2(l,i)]));
          for (int q = 0; q < dim; q++) {
            tan[ijkl] += (dt * FMn[idx_2(i,q)]
                          * (Dgdot_DFe[idx_2(k,l)] * n[idx_2(q,j)]
                             + gdot * Dn_DFe[idx_4(q,j,k,l)]));
          }
        }
      }
    }
  }

  return err;
}

static int bpa_compute_tan_Fe_lam(double *tan,
                                  const double *Fe,
                                  const double *Fp)
{
  int err = 0;
  const double Jp = det3x3(Fp);
  double invFe[tensor] = {};
  err += inv3x3(Fe,invFe);
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      tan[idx_2(i,j)] = -Jp * invFe[idx_2(j,i)];
    }
  }

  return err;
}

static int bpa_compute_tan(double * __restrict tan,
                           const double param_gdot0,
                           const double param_A,
                           const double param_T,
                           const double param_N,
                           const double param_Cr,
                           const double param_alpha,
                           const double dt,
                           const double gdot,
                           const double tau,
                           const double s_s,
                           const double lam,
                           const double * __restrict eff_sig,
                           const double * __restrict n,
                           const double * __restrict Mn,
                           const double * __restrict F,
                           const double * __restrict Fe,
                           const double * __restrict Fp,
                           const HOMMAT *p_hmat)
{
  int err = 0;

  /* compute terms */
  double tan_Fe_Fe[tensor4] = {};
  double tan_Fe_lam[tensor] = {};
  err += bpa_compute_tan_Fe_Fe(tan_Fe_Fe,
                               param_gdot0, param_A, param_T, param_N,
                               param_Cr, param_alpha, dt, gdot, tau, s_s, lam,
                               eff_sig, F, n, Mn, Fe, Fp,
                               p_hmat);
  err += bpa_compute_tan_Fe_lam(tan_Fe_lam, Fe, Fp);

  /* assemble */
  const int last_row = tensor * tan_col;
  for (int i = 0; i < tensor; i++) {
    const int row_start = i * tan_col;
    const int last_col = row_start + tensor;
    memcpy(tan + i * tan_col, tan_Fe_Fe + i * tensor, tensor * sizeof(*tan));
    tan[last_col] = tan[last_row + i] = tan_Fe_lam[i];
  }
  tan[tan_row * tan_col - 1] = 0;

  return err;
}

static int bpa_compute_step1_terms(double *gdot,
                                   double *s_s,
                                   double *tau,
                                   double *Jp,
                                   double *eq_sig_dev,
                                   double *normal,
                                   double *Fp,
                                   const double param_gdot0,
                                   const double param_A,
                                   const double param_T,
                                   const double param_N,
                                   const double param_Cr,
                                   const double param_alpha,
                                   const double s,
                                   const double kappa,
                                   const double *F,
                                   const double *Fe,
                                   const HOMMAT *p_hmat)
{
  int err = 0;
  double Ce[tensor] = {};
  bpa_compute_Fp(Fp,F,Fe);
  *Jp = det3x3(Fp);
  bpa_compute_Ce(Ce,Fe);

  /* compute the pressure */
  double pressure = 0.0;
  double Je = det3x3(Fe);
  bpa_compute_dudj(Je,p_hmat,&pressure);
  pressure *= kappa / 2.0;

  /* compute the pressure-dependent athermal shear stress */
  *s_s = s - param_alpha * pressure;

  /* compute the plastic backstress */
  double Sdev[tensor] = {};
  double Bdev[tensor] = {};
  bpa_compute_Sdev(Ce,p_hmat,Sdev);
  err += bpa_compute_Bdev(Bdev, param_N, param_Cr, Fp);
  err += bpa_compute_loading_dir(normal, eq_sig_dev, tau, Sdev, Bdev, Fe);
  err += bpa_compute_gdot(gdot, param_gdot0, param_A, param_T, *tau, *s_s);
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
static int bpa_int_alg_initial_guess(const double *F,
                                     const Constitutive_model *m,
                                     double *Fe,
                                     double *s,
                                     double *lam)
{
  int err = 0;

  /* copy tensors */
  memcpy(Fe, m->vars_list[0][m->model_id].Fs[_Fe].m_pdata, tensor * sizeof(*Fe));

  /* copy s */
  *s = m->vars_list[0][m->model_id].state_vars->m_pdata[_s];
  *lam = m->vars_list[0][m->model_id].state_vars->m_pdata[_lam];

  /* store the imposed deformation gradient */
  memcpy(m->vars_list[0][m->model_id].Fs[_F].m_pdata, F, tensor * sizeof(*F));

  return err;
}


static int bpa_update_state_variables(const double *F,
                                      const double *Fe,
                                      const double *Fp,
                                      const double s,
                                      const double lam,
                                      Constitutive_model *m)
{
  int err = 0;
  memcpy(m->vars_list[0][m->model_id].Fs[_Fe].m_pdata, Fe, tensor * sizeof(*Fe));
  memcpy(m->vars_list[0][m->model_id].Fs[_Fp].m_pdata, Fp, tensor * sizeof(*Fp));
  memcpy(m->vars_list[0][m->model_id].Fs[_F].m_pdata, F, tensor * sizeof(*F));
  m->vars_list[0][m->model_id].state_vars->m_pdata[_s] = s;
  m->vars_list[0][m->model_id].state_vars->m_pdata[_lam] = lam;
  return err;
}

/**
 * Update the current value of the solution with the increment
 *
 */
static int bpa_update_solution(const double * __restrict increment,
                               double * __restrict Fe,
                               double * __restrict lam)
{
  int err = 0;
  for (int i = 0; i < tensor; i++) {
    Fe[i] += increment[i];
  }
  *lam += increment[tan_row - 1];
  return err;
}

static int bpa_compute_vel_grad(const double * __restrict F,
                                const double * __restrict Fn,
                                const double dt,
                                double * __restrict L,
                                double * __restrict d,
                                double * __restrict ome)
{
  int err = 0;

  /* compute Fdot by BW Euler */
  double Fdot[tensor] = {};
  for (int i = 0; i < tensor; i++) {
    Fdot[i] = (F[i] - Fn[i]) / dt;
  }

  double invF[tensor] = {};
  err += inv3x3(F,invF);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              dim,dim,dim,1.0,Fdot,dim,invF,dim,
              0.0,L,dim);

  /* compute spin and stretch */
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      d[idx_2(i,j)] = 0.5 * (L[idx_2(i,j)] + L[idx_2(j,i)]);
      ome[idx_2(i,j)] = 0.5 * (L[idx_2(i,j)] - L[idx_2(j,i)]);
    }
  }

  return err;
}

static int bpa_compute_Wp(double *Wp,
                          const double gdot,
                          const double * __restrict n,
                          const double * __restrict ome,
                          const double * __restrict d,
                          const double * __restrict Fe)
{
  int err = 0;
  memset(Wp, 0, tensor * sizeof(*Wp));
  double A[tensor4] = {};
  for (int i = 0; i < dim; i++) {
    for (int k = 0; k < dim; k++) {
      for (int p = 0; p < dim; p++) {
        Wp[idx_2(i,k)] += (Fe[idx_2(i,p)] * (ome[idx_2(p,k)] - d[idx_2(p,k)] - gdot * n[idx_2(p,k)])
                           + Fe[idx_2(p,k)] * (ome[idx_2(i,p)] + d[idx_2(i,p)] + gdot * n[idx_2(i,p)]));
        for (int q = 0; q < dim; q++) {
          A[idx_4(i,k,p,q)] = Fe[idx_2(i,p)] * eye[idx_2(k,q)] + eye[idx_2(i,p)] * Fe[idx_2(k,q)];
        }
      }
    }
  }

  err += solve_Ax_b(tensor,tensor,A,Wp);
  if(err) PGFEM_printerr("WARNING: received error (%d) from 'solve_Ax_b'\n",err);
  return err;
}

static int bpa_get_state_at_n(const Constitutive_model *m,
                              double *Mn,
                              double *s_n)
{
  int err = 0;
  err += inv3x3(m->vars_list[0][m->model_id].Fs[_Fp_n].m_pdata,Mn);
  *s_n = m->vars_list[0][m->model_id].state_vars->m_pdata[_s_n];
  return err;
}

static int bpa_compute_Dtau_DFe(double * __restrict Dtau_DFe,
                                const double * __restrict n,
                                const double * __restrict Dsig_DFe)
{
  int err = 0;
  const double rt2 = sqrt(2.0);
  memset(Dtau_DFe, 0, tensor * sizeof(*Dtau_DFe));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          Dtau_DFe[idx_2(i,j)] += n[idx_2(k,l)] / rt2 * Dsig_DFe[idx_4(k,l,i,j)];
        }
      }
    }
  }
  return err;
}

static int bpa_compute_Ds_Dgdot(double * __restrict Ds_Dgdot,
                                const double param_h,
                                const double param_s_ss,
                                const double s_n,
                                const double dt,
                                const double gdot)
{
  int err = 0;
  *Ds_Dgdot =  gdot * (dt * param_h * param_s_ss * (param_s_ss - s_n)
                       / pow(dt * param_h + param_s_ss, 2));
  return err;
}

static int bpa_compute_Ds_DFe(double * __restrict Ds_DFe,
                              const double param_h,
                              const double param_s_ss,
                              const double s_n,
                              const double dt,
                              const double gdot,
                              const double Dgdot_Dss,
                              const double * __restrict Dgdot_DFe)
{
  int err = 0;
  double Ds_Dgdot = 0.0;
  err += bpa_compute_Ds_Dgdot(&Ds_Dgdot, param_h, param_s_ss, s_n, dt, gdot);

  const double coeff = Ds_Dgdot / (1.0 - Ds_Dgdot * Dgdot_Dss);

  memset(Ds_DFe, 0, tensor * sizeof(*Ds_DFe));
  cblas_daxpy(tensor, coeff, Dgdot_DFe, 1, Ds_DFe, 1);
  return err;
}

static int bpa_compute_DM_DFe(double * __restrict DM_DFe,
                              const double param_gdot0,
                              const double param_A,
                              const double param_T,
                              const double param_N,
                              const double param_Cr,
                              const double param_alpha,
                              const double param_h,
                              const double param_s_ss,
                              const double dt,
                              const double s_n,
                              const double s,
                              const double * __restrict Fe,
                              const double * __restrict F,
                              const HOMMAT *p_hmat)
{
  int err = 0;

  /* compute needed terms for derivatives */
  const double kappa = bpa_compute_bulk_mod(p_hmat);
  double gp = 0.0;
  double s_s = 0.0;
  double tau = 0.0;
  double Jp = 0.0;
  double eff_sig[tensor] = {};
  double n[tensor] = {};
  double Fp[tensor] = {};
  err += bpa_compute_step1_terms(&gp, &s_s, &tau, &Jp, eff_sig, n, Fp,
                                 param_gdot0, param_A, param_T, param_N,
                                 param_Cr, param_alpha,
                                 s, kappa, F, Fe, p_hmat);

  /* compute derivative terms */
  double Dgp_DFe[tensor] = {}; /* does not include terms from s */
  double Dsig_DFe[tensor4] = {};
  double Dgp_Dss = 0.0;
  double Dtau_DFe[tensor] = {};
  double Ds_DFe[tensor] = {};
  double Dp_DFe[tensor] = {};
  double Dn_DFe[tensor4] = {};
  err += bpa_compute_Dsig_DFe(Dsig_DFe, F, Fe, Fp, p_hmat, param_N, param_Cr);
  err += bpa_compute_Dp_DFe(Dp_DFe, Fe, p_hmat);
  err += bpa_compute_Dgdot_DFe(Dgp_DFe, param_gdot0, param_A, param_T,
                               param_alpha,s_s, tau, Dsig_DFe, n, Dp_DFe);
  err += bpa_compute_Dgdot_Ds_s(&Dgp_Dss, param_gdot0, param_A, param_T, tau, s_s);
  err += bpa_compute_Dtau_DFe(Dtau_DFe, n, Dsig_DFe);
  err += bpa_compute_Ds_DFe(Ds_DFe, param_h, param_s_ss, s_n, dt, gp, Dgp_Dss, Dgp_DFe);
  err += bpa_compute_Dn_DFe(Dn_DFe, tau, Dsig_DFe, n);

  /* compute full Dgp_DFe */
  cblas_daxpy(tensor, Dgp_Dss, Ds_DFe, 1, Dgp_DFe, 1);

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          const int ijkl = idx_4(i,j,k,l);
          DM_DFe[ijkl] = dt * eye[idx_2(i,j)] * Dgp_DFe[idx_2(k,l)];
          for (int p = 0; p < dim; p++) {
            DM_DFe[ijkl] -= dt * gp * eye[idx_2(i,p)] * Dn_DFe[idx_4(p,j,k,l)];
          }
        }
      }
    }
  }

  return err;
}

/*
 * Private interface for the model
 */
int BPA_int_alg(Constitutive_model *m,
                CM_Ctx &cm_ctx)
{
  int err = 0;
  const HOMMAT *p_hmat = m->param->p_hmat;
  const double kappa = bpa_compute_bulk_mod(p_hmat);
  const double *Fn = m->vars_list[0][m->model_id].Fs[_F_n].m_pdata;
  const double *params = m->param->model_param;

  double s_n = 0;
  double Mn[tensor] = {};
  err += bpa_get_state_at_n(m,Mn,&s_n);

  /* initial guess for the solution */
  double Fe[tensor] = {};
  double lam = 0;
  double s = 0;
  err += bpa_int_alg_initial_guess(cm_ctx.F,m,Fe,&s,&lam);

  /* internal variables */
  double gdot = 0;
  double s_s = 0;
  double tau = 0;
  double Jp = 0;
  double eq_sig_dev[tensor] = {};
  double normal[tensor] = {};
  double Fp[tensor] = {};
  err += bpa_compute_step1_terms(&gdot, &s_s, &tau, &Jp, eq_sig_dev, normal, Fp,         /*< out */
                                 params[mcGdot0], params[mcA], params[mcT], params[mcN], /*< in */
                                 params[mcCr], params[mcAlpha], s, kappa, cm_ctx.F,
                                 Fe, p_hmat);

  /* compute plastic spin */
  double L[tensor] = {};
  double d[tensor] = {};
  double ome[tensor] = {};
  double Wp[tensor] = {};
  err += bpa_compute_vel_grad(cm_ctx.F,Fn,cm_ctx.dt,L,d,ome);
  err += bpa_compute_Wp(Wp,gdot,normal,ome,d,Fe);

  double TAN[tan_row * tan_col] = {};
  double RES[tan_row] = {};
  double norm = 1.0;
  double normWp = 1.0;
  int iter = 0;
  int iterWp = 0;
  int total_it = 0;
  static const double TOL = 1.0e-5;
  static const int maxit = 10;

  /* COMPUTE THE RESIDUAL */
  err += bpa_compute_res_vec(RES,cm_ctx.dt,gdot,lam,Jp,normal,Mn,Wp,cm_ctx.F,Fe);

  while ((normWp > TOL || norm > TOL) && iterWp < maxit) {
    iter = 0;
    while (norm > TOL  && iter < maxit) {
      /* COMPUTE THE TANGENT */
      err += bpa_compute_tan(TAN, params[mcGdot0], params[mcA], params[mcT], params[mcN],
                             params[mcCr], params[mcAlpha], cm_ctx.dt, gdot, tau, s_s, lam,
                             eq_sig_dev, normal, Mn, cm_ctx.F, Fe,
                             Fp, p_hmat);

      if(BPA_PRINT_LEVEL > 1){
        print_array_d(stdout, TAN, tan_row * tan_col, tan_row, tan_col);
        print_array_d(stdout, RES, tan_row, 1, tan_row);
      }

      /* solve for increment: note solution in RES on exit */
      int err_s = 0;
      err_s += solve_Ax_b(tan_row, tan_row, TAN, RES);
      if(err_s) PGFEM_printerr("WARNING: received error (%d) from 'solve_Ax_b'\n",err_s);
      err += err_s;
      if(err_s) goto exit_err;

      /* update deformation */
      err += bpa_update_solution(RES, Fe, &lam);

      /* update vars (gdot, s_s, tau, eq_sig_dev, normal) */
      err += bpa_compute_step1_terms(&gdot, &s_s, &tau, &Jp, eq_sig_dev, normal, Fp,         /*< out */
                                     params[mcGdot0], params[mcA], params[mcT], params[mcN], /*< in */
                                     params[mcCr], params[mcAlpha], s, kappa, cm_ctx.F,
                                     Fe, p_hmat);

      /* COMPUTE RESIDUAL */
      err += bpa_compute_res_vec(RES,cm_ctx.dt,gdot,lam,Jp,normal,Mn,Wp,cm_ctx.F,Fe);
      norm = cblas_dnrm2(tan_row,RES,1);
      if (BPA_PRINT_LEVEL > 0) {
        PGFEM_printf("\tR1 = %6e (%d)\n", norm, iter);
        /* print_array_d(stdout, RES, tan_row, 1, tan_row); */
      }
      iter++;
    }
    total_it += iter;

    double s_k = s;
    s = (s_n + params[mcH] * gdot * (cm_ctx.dt)) / (1 + params[mcH] * gdot * (cm_ctx.dt)/ params[mcSss]);
    normWp = pow((s-s_k) / params[mcSss],2);

    /* compute updated Wp and residual norm */
    err += bpa_compute_Wp(Wp,gdot,normal,ome,d,Fe);
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        normWp += pow(Wp[idx_2(i,j)] + Wp[idx_2(j,i)],2);
      }
    }
    normWp = sqrt(normWp);

    /* COMPUTE NEW RESIDUAL AND NORM */
    err += bpa_compute_step1_terms(&gdot, &s_s, &tau, &Jp, eq_sig_dev, normal, Fp,         /*< out */
                                   params[mcGdot0], params[mcA], params[mcT], params[mcN], /*< in */
                                   params[mcCr], params[mcAlpha], s, kappa, cm_ctx.F,
                                   Fe, p_hmat);

    err += bpa_compute_res_vec(RES,cm_ctx.dt,gdot,lam,Jp,normal,Mn,Wp,cm_ctx.F,Fe);
    norm = cblas_dnrm2(tan_row,RES,1);
    /* output of the iterative procedure */
    if (BPA_PRINT_LEVEL > 0) {
      PGFEM_printf("[%d] R = %6e (%d) || RWp = %6e\n", iterWp, norm, iter - 1, normWp);
    }
    iterWp++;
  }

  /* induce subdivision of PDE if too many iterations or failed to
     converge */
  if (total_it >= 2*maxit
      || norm > TOL
      || normWp > TOL){ err ++;}

  /* Update state variables with converged values */
  err += bpa_update_state_variables(cm_ctx.F,Fe, Fp, s, lam, m);
 exit_err:
  return err;
}

int BPA_PARAM::compute_dev_stress(const Constitutive_model *m,
                                  CM_Ctx &cm_ctx,
                                  double *stress)
const
{
  int err = 0;
  double Ce[tensor] = {};
  bpa_compute_Ce(Ce,m->vars_list[0][m->model_id].Fs[_Fe].m_pdata);
  bpa_compute_Sdev(Ce,m->param->p_hmat,stress);
  return err;
}

int BPA_PARAM::compute_dudj(const Constitutive_model *m,
                            CM_Ctx &cm_ctx,
                            double *dudj)
const
{
  int err = 0;
  const double Je = det3x3(m->vars_list[0][m->model_id].Fs[_Fe].m_pdata);
  bpa_compute_dudj(Je,m->param->p_hmat,dudj);
  return err;
}

int BPA_PARAM::compute_dev_tangent(const Constitutive_model *m,
                                   CM_Ctx &cm_ctx,
                                   double *L)
const
{
  int err = 0;
  double Ce[tensor] = {};
  bpa_compute_Ce(Ce,m->vars_list[0][m->model_id].Fs[_Fe].m_pdata);
  bpa_compute_Ldev(Ce,m->param->p_hmat, L);
  return err;
}

int BPA_PARAM::compute_d2udj2(const Constitutive_model *m,
                              CM_Ctx &cm_ctx,
                              double *d2udj2)
const
{
  int err = 0;
  const double Je = det3x3(m->vars_list[0][m->model_id].Fs[_Fe].m_pdata);
  bpa_compute_d2udj2(Je,m->param->p_hmat,d2udj2);
  return err;
}

int BPA_PARAM::update_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  memcpy(Fs[_Fe].m_pdata, Fs[_Fe_n].m_pdata, sizeof(double)*tensor);
  memcpy(Fs[_Fp].m_pdata, Fs[_Fp_n].m_pdata, sizeof(double)*tensor);
  memcpy(Fs[_F ].m_pdata, Fs[_F_n ].m_pdata, sizeof(double)*tensor);

  /* alias */
  double *vars = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  vars[_s_n]   = vars[_s];
  vars[_lam_n] = vars[_lam];
  return err;
}

int BPA_PARAM::reset_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  memcpy(Fs[_Fe_n].m_pdata, Fs[_Fe].m_pdata, sizeof(double)*tensor);
  memcpy(Fs[_Fp_n].m_pdata, Fs[_Fp].m_pdata, sizeof(double)*tensor);
  memcpy(Fs[_F_n ].m_pdata, Fs[_F ].m_pdata, sizeof(double)*tensor);

  /* alias */
  double *vars = m->vars_list[0][m->model_id].state_vars[0].m_pdata;
  vars[_s  ] = vars[_s_n  ];
  vars[_lam] = vars[_lam_n];
  return err;
}

int BPA_PARAM::get_var_info(Model_var_info &info)
const
{
  int err = 0;

  info.n_Fs = _n_Fs;
  info.n_vars = _n_vars;
  info.n_flags = _n_flags;
  info.F_names = PGFEM_malloc<char*>(_n_Fs);
  info.var_names = PGFEM_malloc<char*>(_n_vars);
  info.flag_names = PGFEM_malloc<char*>(_n_flags);

  /* allocate/copy strings */
  info.F_names[_Fe]      = strdup("Fe");
  info.F_names[_Fe_n]    = strdup("Fe_n");
  info.F_names[_Fp]      = strdup("Fp");
  info.F_names[_Fp_n]    = strdup("Fp_n");
  info.F_names[_F]       = strdup("F");
  info.F_names[_F_n]     = strdup("F_n");
  info.var_names[_s]     = strdup("s");
  info.var_names[_s_n]   = strdup("s_n");
  info.var_names[_lam]   = strdup("lam");
  info.var_names[_lam_n] = strdup("lam_n");
  return err;
}

int BPA_PARAM::get_pF(const Constitutive_model *m,
                      double *F,
                      const int stepno)
const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  switch(stepno)
  {
    case 0: // n-1
      memcpy(F,Fs[_Fp_n].m_pdata, tensor*sizeof(double));
      break;
    case 1: // n
      memcpy(F,Fs[_Fp_n].m_pdata, tensor*sizeof(double));
      break;
    case 2: // n+1
      memcpy(F,Fs[_Fp].m_pdata, tensor*sizeof(double));
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;
}



int BPA_PARAM::get_F(const Constitutive_model *m,
                     double *F,
                     const int stepno)
const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  switch(stepno)
  {
    case 0: // n-1
      memcpy(F,Fs[_F_n].m_pdata,tensor*sizeof(double));
      break;
    case 1: // n
      memcpy(F,Fs[_F_n].m_pdata,  tensor*sizeof(double));
      break;
    case 2: // n+1
      memcpy(F,Fs[_F].m_pdata,tensor*sizeof(double));
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;
}

int BPA_PARAM::get_eF(const Constitutive_model *m,
                      double *eF_in,
                      const int stepno)
const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  switch(stepno)
  {
    case 0: // n-1
      memcpy(eF_in,Fs[_Fe_n].m_pdata,tensor*sizeof(double));
      break;
    case 1: // n
      memcpy(eF_in,Fs[_Fe_n].m_pdata,  tensor*sizeof(double));
      break;
    case 2: // n+1
      memcpy(eF_in,Fs[_Fe].m_pdata,tensor*sizeof(double));
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;
}


int BPA_PARAM::get_hardening(const Constitutive_model *m,
                             double *var,
                             const int stepno)
const
{
  int err = 0;
  double *s_var = m->vars_list[0][m->model_id].state_vars->m_pdata;
  switch(stepno)
  {
    case 0: // n-1
      *var = s_var[_s_n];
      break;
    case 1: // n
      *var = s_var[_s_n];
      break;
    case 2: // n+1
      *var = s_var[_s];
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;
}

int BPA_PARAM::compute_dMdu(const Constitutive_model *m,
                            CM_Ctx &cm_ctx,
                            double *Grad_op,
                            const int nne,
                            const int ndofn,
                            double *dM_du)
const
{
  int err = 0;
  const double *F = cm_ctx.F;
  const double dt = cm_ctx.dt;
  const double *Fe = m->vars_list[0][m->model_id].Fs[_Fe].m_pdata;
  const double *Fp = m->vars_list[0][m->model_id].Fs[_Fp].m_pdata;
  const double *Fp_n = m->vars_list[0][m->model_id].Fs[_Fp_n].m_pdata;
  const double s = m->vars_list[0][m->model_id].state_vars->m_pdata[_s];
  const double s_n = m->vars_list[0][m->model_id].state_vars->m_pdata[_s_n];
  const double *p = m->param->model_param;


  /* compute 4th order tensor operators */
  double M[tensor] = {};
  double Mn[tensor] = {};
  double DM_DFe[tensor4] = {};
  // double DFe_DF[tensor4] = {};
  // double DFe_DM[tensor4] = {};
  err += inv3x3(Fp,M);
  err += inv3x3(Fp_n,Mn);
  err += bpa_compute_DM_DFe(DM_DFe, p[mcGdot0], p[mcA], p[mcT], p[mcN],
                            p[mcCr], p[mcAlpha], p[mcH], p[mcSss],
                            dt, s_n, s, Fe, F, m->param->p_hmat);

  /* compute the tangent operator */
  double U[tensor4] = {};
  double B[tensor4] = {};
  for (int i = 0; i < dim; i++) {
    for(int j = 0; j < dim; j++) {
      for(int m = 0; m < dim; m++) {
        for(int n = 0; n < dim; n++) {
          const int ijmn = idx_4(i,j,m,n);
          U[ijmn] = eye[idx_2(i,m)] * eye[idx_2(n,j)];
          for(int l = 0; l < dim; l++) {
            U[ijmn] -= DM_DFe[idx_4(i,j,l,n)] * F[idx_2(l,m)];
            B[ijmn] += DM_DFe[idx_4(i,j,m,l)] * M[idx_2(n,l)];
          }
        }
      }
    }
  }

  /* NOTE that we assemble the multiple RHS (and the resulting
     solution) as rows instead of columns. Therefore, we do not need
     to transpose the RHS before passing to LAPACK, nor the solution
     after returning from LAPACK */
  for (int a = 0; a < nne; a++) {
    for (int b = 0; b < ndofn; b++) {
      const int idx_ab = idx_4_gen(a,b,0,0,nne,ndofn,dim,dim);
      for (int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          const int ijmn = idx_4(i,j,0,0);
          const int ij = idx_2(i,j);
          *(dM_du + idx_ab + ij) = cblas_ddot(tensor, B + ijmn, 1, Grad_op + idx_ab, 1);
        }
      }
    }
  }

  /* LAPACK solution procedure:
   * 1) Transpose U
   * 2) Allocate workspace
   * 3) Solve
   * 4) Clean up
   */
  {
    int l_err = 0;
    double Ut[tensor4] = {};
    transpose(Ut,U,tensor,tensor);
    int *IPIV = PGFEM_malloc<int>(tensor);
    int NRHS = nne * ndofn;
    int DIM = tensor;
#ifdef ARCH_BGQ
    dgesv(DIM, NRHS, Ut, DIM, IPIV, dM_du, DIM, &l_err);
#else
    dgesv(&DIM, &NRHS, Ut, &DIM, IPIV, dM_du, &DIM, &l_err);
#endif
    assert(l_err >= 0);
    err += l_err;
    free(IPIV);
  }

  return err;
}

int BPA_PARAM::read_param(FILE *in)
const
{
  int err = 0;

  /* get pointer to parameter data */
  double *param = this->model_param;
  assert(param != NULL); // check the pointer

  /* scan to non-blank/comment line */
  err += scan_for_valid_line(in);

  /* READ PROPERTIES IN ALPHABETICAL ORDER */
  int match = fscanf(in, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                     param + mcA, param + mcAlpha, param + mcCr,
                     param + mcGdot0, param + mcH, param + mcN,
                     param + mcS0, param + mcSss, param + mcT);
  if (match != N_PARAM) err++;
  assert(match == N_PARAM && "Did not read expected number of parameters");

  /* scan past any other comment/blank lines in the block */
  err += scan_for_valid_line(in);

  /* not expecting EOF, check and return error if encountered */
  if (feof(in)) err ++;
  assert(!feof(in) && "EOF reached prematurely");

  return err;
}

int BPA_PARAM::set_init_vals(Constitutive_model *m)
const
{
  int err = 0;

  /* s is s0 at start */
  m->vars_list[0][m->model_id].state_vars->m_pdata[_s] = m->param->model_param[mcS0];
  m->vars_list[0][m->model_id].state_vars->m_pdata[_s_n] = m->param->model_param[mcS0];

  m->vars_list[0][m->model_id].state_vars->m_pdata[_lam] = 0;
  m->vars_list[0][m->model_id].state_vars->m_pdata[_lam_n] = 0;

  return err;
}

int BPA_PARAM::write_restart(FILE *out, const Constitutive_model *m)
const
{
  /* write all state variables at n */
  int err = 0;
  const double *Fen = m->vars_list[0][m->model_id].Fs[_Fe_n].m_pdata;
  const double *Fpn = m->vars_list[0][m->model_id].Fs[_Fp_n].m_pdata;
  const double *Fn = m->vars_list[0][m->model_id].Fs[_F_n].m_pdata;
  const double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  if(fprintf(out,"%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
             Fen[0], Fen[1], Fen[2], Fen[3], Fen[4],
             Fen[5], Fen[6], Fen[7], Fen[8]) < 0) err ++;
  if(fprintf(out,"%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
             Fpn[0], Fpn[1], Fpn[2], Fpn[3], Fpn[4],
             Fpn[5], Fpn[6], Fpn[7], Fpn[8]) < 0) err ++;
  if(fprintf(out,"%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
             Fn[0], Fn[1], Fn[2], Fn[3], Fn[4],
             Fn[5], Fn[6], Fn[7], Fn[8]) < 0) err ++;
  if(fprintf(out,"%.17e %.17e\n",vars[_s_n], vars[_lam_n]) < 0) err++;
  return err;
}

int BPA_PARAM::read_restart(FILE *in, Constitutive_model *m)
const
{
  /* read all state variables at n and set all vars at n+1 = n */
  int err = 0;
  double *Fen = m->vars_list[0][m->model_id].Fs[_Fe_n].m_pdata;
  double *Fpn = m->vars_list[0][m->model_id].Fs[_Fp_n].m_pdata;
  double *Fn = m->vars_list[0][m->model_id].Fs[_F_n].m_pdata;
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;

  if(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
            Fen, Fen + 1, Fen + 2, Fen + 3, Fen + 4,
            Fen + 5, Fen + 6, Fen + 7, Fen + 8) != tensor) err++;
  if(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
            Fpn, Fpn + 1, Fpn + 2, Fpn + 3, Fpn + 4,
            Fpn + 5, Fpn + 6, Fpn + 7, Fpn + 8) != tensor) err++;
  if(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
            Fn, Fn + 1, Fn + 2, Fn + 3, Fn + 4,
            Fn + 5, Fn + 6, Fn + 7, Fn + 8) != tensor) err++;
  if(fscanf(in,"%lf %lf", vars + _s_n, vars + _lam_n) != 2) err++;
  err += this->reset_state_vars(m);
  return err;
}

int BPA_PARAM::update_elasticity(const Constitutive_model *m,
                                 CM_Ctx &cm_ctx,
                                 double *L,
                                 double *S,
                                 const bool compute_stiffness)
const
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  err += constitutive_model_default_update_elasticity(m, Fs[_Fe].m_pdata, L, S, compute_stiffness);
  return err;
}

/*
 * Public interface for the BPA model
 */
int BPA_PARAM::model_dependent_initialization(void)
{
  int err = 0;
  this->type              = BPA_PLASTICITY;
  this->n_param           = N_PARAM;
  this->model_param       = new double[N_PARAM]();
  return err;
}

int BPA_PARAM::integration_algorithm(Constitutive_model *m,
                                     CM_Ctx &cm_ctx)
const
{
  return BPA_int_alg(m, cm_ctx);
}
