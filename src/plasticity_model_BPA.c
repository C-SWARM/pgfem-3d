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
#define tangent 10

Define_Matrix(double);

/* Set to value > 0 for extra diagnostics/printing */
static const int BPA_PRINT_LEVEL = 1;

/* constants/enums */
static const int _n_Fs = 8;
static const int _n_vars = 2;
enum {_M,_W,_Fe,_F,_M_n,_W_n,_Fe_n,_F_n};
enum {_s,_s_n};
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

static void bpa_compute_Ce(double *Ce,
                           const double *Fe)
{
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
              dim,dim,dim,1.0,Fe,dim,Fe,dim,
              0.0,Ce,dim);
}

static double bpa_compute_bulk_mod(const HOMMAT *mat)
{
  return ( (2* mat->G * (1 + mat->nu)) / (3 * (1 - 2 * mat->nu)) );
}

static double bpa_compute_plam(const double *Cp)
{
  return sqrt(1./3. * (Cp[0] + Cp[4] + Cp[8]));
}

static void bpa_compute_Cp(double * restrict Cp,
                           const double * restrict Fp)
{
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
              dim,dim,dim,1.0,Fp,dim,Fp,dim,
              0.0,Cp,dim);
}

static void bpa_compute_Cpdev(double * restrict Cpdev,
                              const double * restrict Cp)
{
  const double Jdev = pow(det3x3(Cp),-1./3.);
  for (int i = 0; i < tensor; i++) {
    Cpdev[i] = Jdev * Cp[i];
  }
}

static void bpa_compute_Fp(double * restrict Fp,
                           const double * restrict F,
                           const double * restrict Fe)
{
  double invFe[tensor] = {};
  inv3x3(Fe,invFe);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              dim,dim,dim,1.0,invFe,dim,F,dim,
              0.0,Fp,dim);
}

/**
 * Compute chain rule term DM_DFe
 */
static void bpa_compute_DM_DFe(double * restrict DM_DFe,
                               const double * restrict F,
                               const double * restrict Fe)
{
  double invF[tensor] = {};
  inv3x3(F,invF);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              dim,dim,dim,1.0,invF,dim,Fe,dim,
              0.0,DM_DFe,dim);
}

static void bpa_compute_M(double * restrict M,
                          const double * restrict F,
                          const double * restrict Fe)
{
  double invF[tensor] = {};
  inv3x3(F,invF);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim,
              1.0, invF, dim, Fe, dim, 0.0, M, dim);
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
                                   double *sdot,
                                   double *tau,
                                   double *eq_sig_dev,
                                   double *normal,
                                   const double s,
                                   const double kappa,
                                   const double *F,
                                   const double *Fe,
                                   const HOMMAT *p_hmat)
{
  int err = 0;
  double Fp[tensor] = {};
  double Ce[tensor] = {};
  bpa_compute_Fp(Fp,F,Fe);
  bpa_compute_Ce(Ce,Fe);

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
  err += BPA_compute_sdot(sdot,s,*gdot);
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
                                     double *s)
{
  int err = 0;

  /* copy tensors */
  memcpy(Fe, m->vars.Fs[_Fe].m_pdata, tensor * sizeof(*Fe));

  /* copy s */
  *s = m->vars.state_vars->m_pdata[_s];

  /* store the imposed deformation gradient */
  memcpy(m->vars.Fs[_F].m_pdata, F, tensor * sizeof(*F));

  return err;
}


static int bpa_update_state_variables(const double *F,
                                      const double *Fe,
                                      const double *Wp,
                                      const double s,
                                      Constitutive_model *m)
{
  int err = 0;
  double M[tensor] = {};
  bpa_compute_M(M,F,Fe);
  memcpy(m->vars.Fs[_M].m_pdata, M, tensor * sizeof(*M));
  memcpy(m->vars.Fs[_Fe].m_pdata, Fe, tensor * sizeof(*Fe));
  memcpy(m->vars.Fs[_W].m_pdata, Wp, tensor * sizeof(*Wp));
  m->vars.state_vars->m_pdata[_s] = s;
  return err;
}

/**
 * Update the current value of the solution with the increment
 *
 */
static int bpa_update_solution(const double * restrict increment,
                               double * restrict Fe,
                               double * restrict s)
{
  int err = 0;
  for (int i = 0; i < tensor; i++) {
    Fe[i] += increment[i];
  }
  *s += increment[tangent - 1];
  return err;
}

static int bpa_compute_vel_grad(const double * restrict F,
                                const double * restrict Fn,
                                const double dt,
                                double * restrict L,
                                double * restrict d,
                                double * restrict ome)
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
                          const double * restrict n,
                          const double * restrict ome,
                          const double * restrict d,
                          const double * restrict Fe)
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
  assert(err == 0);
  return err;
}

/*
 * Private interface for the model
 */
int BPA_int_alg(Constitutive_model *m,
                const void *ctx)
{
  int err = 0;
  const BPA_ctx *CTX = ctx;
  const HOMMAT *p_hmat = m->param->p_hmat;
  const double kappa = bpa_compute_bulk_mod(p_hmat);
  const double s_n = m->vars.state_vars->m_pdata[_s_n];
  const double *Mn = m->vars.Fs[_M_n].m_pdata; /* alias */
  const double *Fn = m->vars.Fs[_F_n].m_pdata; /* alias */

  /* initial guess for the solution */
  double Fe[tensor] = {};
  double s = 0;
  err += bpa_int_alg_initial_guess(CTX->F,m,Fe,&s);

  /* internal variables */
  double gdot = 0;
  double s_s = 0;
  double tau = 0;
  double sdot = 0;
  double eq_sig_dev[tensor] = {};
  double normal[tensor] = {};
  err += bpa_compute_step1_terms(&gdot, &s_s, &sdot, &tau, eq_sig_dev, normal, /*< out */
                                 s, kappa, CTX->F, Fe, p_hmat);          /*< in */

  /* compute plastic spin */
  double L[tensor] = {};
  double d[tensor] = {};
  double ome[tensor] = {};
  double Wp[tensor] = {};
  err += bpa_compute_vel_grad(CTX->F,Fn,CTX->dt,L,d,ome);
  err += bpa_compute_Wp(Wp,gdot,normal,ome,d,Fe);

  double TAN[tangent * tangent] = {};
  double RES[tangent] = {};
  double norm = 1.0;
  double normWp = 1.0;
  int iter = 0;
  int iterWp = 0;
  static const double TOL = 1.0e-5;
  static const int maxit = 15;

  /* compute the residual */
  err += BPA_int_alg_res(RES, CTX->dt, gdot, sdot, s_n, s,
                         normal, CTX->F, Fe, Wp, Mn);

  while (normWp > TOL && iterWp < maxit) {
    while (norm > TOL  && iter < maxit) {

      /* compute the tangent */
      err += BPA_int_alg_tan(TAN, CTX->dt, gdot, tau, s, s_s,
                             eq_sig_dev, Mn, normal, CTX->F, Fe, p_hmat);

      if(BPA_PRINT_LEVEL > 1){
        print_array_d(stdout,TAN,tangent * tangent, tangent,tangent);
        print_array_d(stdout,RES,tangent,1,tangent);
      }

      /* solve for increment: note solution in RES on exit */
      int err_s = 0;
      err_s += solve_Ax_b(tangent,tangent,TAN,RES);
      if(err_s) PGFEM_printerr("WARNING: received error (%d) from 'solve_Ax_b'\n",err_s);
      assert(err_s == 0);
      err += err_s;

      /* update deformation */
      err += bpa_update_solution(RES, Fe, &s);

      /* update vars (gdot, s_s, tau, eq_sig_dev, normal) */
      err += bpa_compute_step1_terms(&gdot, &s_s, &sdot, &tau, eq_sig_dev, normal, /*< out */
                                     s, kappa, CTX->F, Fe, p_hmat);          /*< in */

      /* compute residual */
      err += BPA_int_alg_res(RES, CTX->dt, gdot, sdot, s_n, s,
                             normal, CTX->F, Fe, Wp, Mn);

      norm = cblas_dnrm2(tangent,RES,1);
      if (BPA_PRINT_LEVEL > 0) {
        printf("\tR1 = %6e (%d)\n", norm, iter);
      }
      iter++;
    }
    assert(norm < TOL);

    /* compute updated Wp and residual norm */
    err += bpa_compute_Wp(Wp,gdot,normal,ome,d,Fe);
    normWp = 0.0;
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        normWp += pow(Wp[idx_2(i,j)] + Wp[idx_2(j,i)],2);
      }
    }
    normWp = sqrt(normWp);

    /* compute new residual and norm */
    err += BPA_int_alg_res(RES, CTX->dt, gdot, sdot, s_n, s,
                           normal, CTX->F, Fe, Wp, Mn);
    norm = cblas_dnrm2(tangent,RES,1);
    /* output of the iterative procedure */
    if (BPA_PRINT_LEVEL > 0) {
      printf("[%d] R = %6e (%d) || RWp = %6e\n", iterWp, norm, iter - 1, normWp);
    }
    iterWp++;
  }
  assert(normWp < TOL);

  /* Update state variables with converged values */
  err += bpa_update_state_variables(CTX->F,Fe, Wp, s, m);
  return err;
}

int BPA_int_alg_res_terms(double * restrict RFe,
                          double * restrict Rs,
                          const double dt,
                          const double gdot,
                          const double sdot,
                          const double s_n,
                          const double s,
                          const double * restrict n,
                          const double * restrict F,
                          const double * restrict Fe,
                          const double * restrict Wp,
                          const double * restrict Mn)
{
  int err = 0;
  memset(RFe, 0, tensor * sizeof(*RFe));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      RFe[idx_2(i,j)] += 2 * Fe[idx_2(i,j)] - Fe[idx_2(j,i)];
      for (int k = 0; k < dim; k++) {
        RFe[idx_2(i,j)] -= F[idx_2(i,k)] * Mn[idx_2(k,j)];
        for (int p = 0; p < dim; p++) {
          RFe[idx_2(i,j)] += dt * F[idx_2(i,k)] * Mn[idx_2(k,p)] * (gdot * n[idx_2(p,j)] + Wp[idx_2(p,j)]);
        }
      }
    }
  }
  *Rs = (s - s_n - sdot * dt) / param_s_ss;
  return err;
}

static int bpa_compute_DSdev_DFe(double * restrict DSdev_DFe,
                                 const double * restrict Ce,
                                 const double * restrict Fe,
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

  return 0;
}

static int bpa_compute_Dgdot_DFe(double * restrict Dgdot_DFe,
                                 const double s_s,
                                 const double tau,
                                 const double * restrict Dsig_DFe,
                                 const double * restrict sig)
{
  int err = 0;
  double Dgdot_Dtau = 0;
  err += BPA_compute_Dgdot_Dtau(&Dgdot_Dtau,s_s,tau);
  double Dtau_Dsig[tensor] = {};
  err += BPA_compute_Dtau_Dsig(Dtau_Dsig,sig);

  memset(Dgdot_DFe, 0, tensor * sizeof(*Dgdot_DFe));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          Dgdot_DFe[idx_2(i,j)] += Dgdot_Dtau * Dtau_Dsig[idx_2(k,l)] * Dsig_DFe[idx_4(k,l,i,j)];
        }
      }
    }
  }

  return err;
}

static int bpa_compute_Dn_DFe(double * restrict Dn_DFe,
                              const double tau,
                              const double * restrict Dsig_DFe,
                              const double * restrict n)
{
  int err = 0;
  double Dn_Dsig[tensor4] = {};
  err += BPA_compute_Dn_Dsig(Dn_Dsig,n,tau);
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


static int bpa_compute_Dsig_DFe(double * restrict Dsig_DFe,
                                const double * restrict F,
                                const double * restrict Fe,
                                const HOMMAT *p_hmat)
{
  int err = 0;

  /* compute deformation */
  const double inv_Je = 1.0 / det3x3(Fe);
  double invFe[tensor] = {};
  err += inv3x3(Fe,invFe);
  double Fp[tensor] = {};
  double Ce[tensor] = {};
  bpa_compute_Ce(Ce,Fe);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              dim,dim,dim,1.0,invFe,dim,Fe,dim,
              0.0,Fp,dim);

  /* compute stress tensors */
  double Sdev[tensor] = {};
  double Bdev[tensor] = {};
  bpa_compute_Sdev(Ce,p_hmat,Sdev);
  err += BPA_compute_Bdev(Bdev,Fp);

  /* compute derivatives */
  double DB_DM[tensor4] = {};
  double DSdev_DFe[tensor4] = {};
  double DM_DFe[tensor4] = {};
  err += BPA_compute_DBdev_DM(DB_DM,Fp);
  err += bpa_compute_DSdev_DFe(DSdev_DFe,Ce,Fe,p_hmat);
  bpa_compute_DM_DFe(DM_DFe,F,Fe);

  /* compute other terms */
  double SmB[tensor] = {};
  for(int i = 0; i < tensor; i++){
    SmB[i] = Sdev[i] - Bdev[i];
  }

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
                  Dsig_DFe[idx_4(i,j,k,l)] += - (inv_Je * Fe[idx_2(i,p)] * DB_DM[idx_4(p,q,r,s)] * DM_DFe[idx_4(r,s,k,l)] * Fe[idx_2(j,q)]);
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

static int bpa_compute_Dp_DFe(double * restrict Dp_DFe,
                              const double * restrict Fe,
                              const HOMMAT *p_hmat)
{
  int err = 0;
  const double kappa = bpa_compute_bulk_mod(p_hmat);
  const double Je = det3x3(Fe);
  double invFe[tensor] = {};
  err += inv3x3(Fe,invFe);

  double d2udj2 = 0;
  bpa_compute_d2udj2(Je,p_hmat,&d2udj2);

  for(int i = 0; i < tensor; i++){
    Dp_DFe[i] = 0.5 * kappa * d2udj2 * invFe[i];
  }
  return err;
}

int BPA_compute_Dgdot_Dtau(double *Dgdot_Dtau,
                           const double s_s,
                           const double tau)
{
  int err = 0;
  const double pow_term = pow(tau / s_s, 5.0 / 6.0);
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
  if (mag <= 0) {
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
  const double coef = 1.0 / (sqrt(2) * tau);
  assert(isfinite(coef));

  memset(Dn_Dsig, 0, tensor4 * sizeof(*Dn_Dsig));
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
              DCpdev_DM[ijkl] -= (0.5 * (eye[idx_2(i,p)] *  eye[idx_2(q,j)]
                                         + eye[idx_2(i,q)] *  eye[idx_2(p,j)])
                                  * (Fp[idx_2(p,k)] * Cpdev[idx_2(l,q)] 
                                     + Cpdev[idx_2(p,l)] * Fp[idx_2(q,k)])
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

int BPA_int_alg_tan(double * restrict TAN,
                    const double dt,
                    const double gdot,
                    const double tau,
                    const double s,
                    const double s_s,
                    const double *sig,
                    const double *Mn,
                    const double *n,
                    const double *F,
                    const double *Fe,
                    const HOMMAT *p_hmat)
{
  int err = 0;

  /* compute terms */
  double DRFe_Fe[tensor4] = {};
  double DRFe_s[tensor] = {};
  double DRs_Fe[tensor] = {};
  double DRs_s = 0;
  err += BPA_int_alg_tan_terms(DRFe_Fe, DRFe_s, DRs_Fe, &DRs_s,
                               dt, gdot, tau, s, s_s, sig, Mn, n, F, Fe, p_hmat);

  /* !!! REWRITE !!! */
  /* Assemble tangent blocks, TAN [19x19]
     TAN = T------------------T
           | DRFe_Fe | DRFe_s |
           |         |        |
           |---------|--------|
           | DRs_Fe  |  DRs_s |
           L------------------J
   */
  const int row2_start = tensor * tangent;
  for (int i = 0; i < tensor; i++) {
    const int blk1_start = i * tangent;
    const int blk2_start = blk1_start + tensor;

    /* 1st row */
    memcpy(TAN + blk1_start, DRFe_Fe + i * tensor, tensor * sizeof(*TAN));
    TAN[blk2_start] = DRFe_s[i];

    /* 2nd row */
    TAN[row2_start + i] = DRs_Fe[i];
  }
  TAN[tangent * tangent - 1] = DRs_s;

  return err;
}

int BPA_int_alg_res(double *RES,
                    const double dt,
                    const double gdot,
                    const double sdot,
                    const double s_n,
                    const double s,
                    const double * restrict n,
                    const double * restrict F,
                    const double * restrict Fe,
                    const double * restrict Wp,
                    const double * restrict Mn)
{
  int err = 0;
  /*
   * OPTIMIZATION NOTE: After debugging/error-checking with separate
   * arrays, can pass in appropriately indexed pointers to locations
   * in RES to remove additional assembly/memory overhead.
   */
  double RFe[tensor] = {};
  double Rs = 0;
  err += BPA_int_alg_res_terms(RFe,&Rs,dt,gdot,sdot,s_n,s,
                               n,F,Fe,Wp,Mn);
  memcpy(RES, RFe, tensor * sizeof(*RES));
  RES[tangent - 1] = Rs;

  /* OPTIMIZED */
  /* err += BPA_int_alg_res_terms(RES, RES + tangent - 1, */
  /*                              dt,gdot,sdot,s_n,s,n,F,Fe,Wp,Mn); */

  for(int i = 0; i < tangent; i++){
    RES[i] *= -1.0;
  }
  return err;
}

int BPA_dev_stress(const Constitutive_model *m,
                   const void *ctx,
                   Matrix_double *dev_stress)
{
  int err = 0;
  double Ce[tensor] = {};
  bpa_compute_Ce(Ce,m->vars.Fs[_Fe].m_pdata);
  bpa_compute_Sdev(Ce,m->param->p_hmat,dev_stress->m_pdata);
  return err;
}

int BPA_dudj(const Constitutive_model *m,
             const void *ctx,
             double *dudj)
{
  int err = 0;
  const double Je = det3x3(m->vars.Fs[_Fe].m_pdata);
  bpa_compute_dudj(Je,m->param->p_hmat,dudj);
  return err;
}

int BPA_dev_tangent(const Constitutive_model *m,
                    const void *ctx,
                    Matrix_double *dev_tangent)
{
  int err = 0;
  double Ce[tensor] = {};
  bpa_compute_Ce(Ce,m->vars.Fs[_Fe].m_pdata);
  bpa_compute_Ldev(Ce,m->param->p_hmat,dev_tangent->m_pdata);
  return err;
}

int BPA_d2udj2(const Constitutive_model *m,
               const void *ctx,
               double *d2udj2)
{
  int err = 0;
  const double Je = det3x3(m->vars.Fs[_Fe].m_pdata);
  bpa_compute_d2udj2(Je,m->param->p_hmat,d2udj2);
  return err;
}

int BPA_update_vars(Constitutive_model *m)
{
  int err = 0;
  Matrix_copy(m->vars.Fs[_M_n],m->vars.Fs[_M]);
  Matrix_copy(m->vars.Fs[_W_n],m->vars.Fs[_W]);
  Matrix_copy(m->vars.Fs[_Fe_n],m->vars.Fs[_Fe]);
  Matrix_copy(m->vars.Fs[_F_n],m->vars.Fs[_F]);

  /* alias */
  Vector_double *vars = m->vars.state_vars;
  Vec_v(*vars,_s_n + 1) = Vec_v(*vars,_s + 1);
  return err;
}

int BPA_reset_vars(Constitutive_model *m)
{
  int err = 0;
  Matrix_copy(m->vars.Fs[_M],m->vars.Fs[_M_n]);
  Matrix_copy(m->vars.Fs[_W],m->vars.Fs[_W_n]);
  Matrix_copy(m->vars.Fs[_Fe],m->vars.Fs[_Fe_n]);
  Matrix_copy(m->vars.Fs[_F],m->vars.Fs[_F_n]);

  /* alias */
  Vector_double *vars = m->vars.state_vars;
  Vec_v(*vars,_s + 1) = Vec_v(*vars,_s_n + 1);
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
  (*info)->F_names[_F] = strdup("F");
  (*info)->F_names[_F_n] = strdup("F_n");
  (*info)->var_names[_s_n] = strdup("s_n");
  (*info)->var_names[_s] = strdup("s");
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

static int bpa_compute_Dgdot_Ds(double *Dgdot_Ds,
                                const double tau,
                                const double s_s)
{
  int err = 0;
  const double pow_term = pow(tau / s_s, 5. / 6.);
  double gdot = 0.0;
  err += BPA_compute_gdot(&gdot,tau,s_s);
  *Dgdot_Ds = gdot * param_A * (pow_term - 6.) / (6. * param_T);
  return err;
}

static int bpa_compute_Dsdot_Dgdot(double *Dsdot_Dgdot,
                                   const double s)
{
  int err = 0;
  *Dsdot_Dgdot = param_h * (1 - s / param_s_ss);
  return err;
}

static int bpa_compute_Dsdot_Ds(double *Dsdot_Ds,
                                const double s,
                                const double gdot,
                                const double Dgdot_Ds)
{
  int err = 0;
  *Dsdot_Ds = param_h * (1 - s / param_s_ss) * Dgdot_Ds - param_h * gdot / param_s_ss;
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

int BPA_int_alg_tan_terms(double * restrict DRFe_Fe,
                          double * restrict DRFe_s,
                          double * restrict DRs_Fe,
                          double * restrict DRs_s,
                          const double dt,
                          const double gdot,
                          const double tau,
                          const double s,
                          const double s_s,
                          const double * restrict sig,
                          const double * restrict Mn,
                          const double * restrict n,
                          const double * restrict F,
                          const double * restrict Fe,
                          const HOMMAT *p_hmat)
{
  int err = 0;

  /* declare dervatives */
  double sig_Fe[tensor4] = {};
  double n_Fe[tensor4] = {};
  double gdot_Fe[tensor] = {};
  double gdot_s = 0;
  double sdot_gdot = 0;
  double sdot_s = 0;

  /* compute derivatives */
  err += bpa_compute_Dsig_DFe(sig_Fe,F,Fe,p_hmat);
  err += bpa_compute_Dn_DFe(n_Fe,tau,sig_Fe,n);
  err += bpa_compute_Dgdot_DFe(gdot_Fe,s_s,tau,sig_Fe,sig);
  err += bpa_compute_Dgdot_Ds(&gdot_s, tau, s_s);
  err += bpa_compute_Dsdot_Dgdot(&sdot_s,s);
  err += bpa_compute_Dsdot_Ds(&sdot_s,s,gdot,gdot_s);

  /* compute intermediate terms */
  double FMn[tensor] = {};
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              dim,dim,dim,1.0,F,dim,Mn,dim,
              0.0,FMn,dim);

  /* comptue tangent terms */
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      DRs_Fe[idx_2(i,j)] -= dt * sdot_gdot * gdot_Fe[idx_2(i,j)] / param_s_ss;
      for (int k = 0; k < dim; k++) {
        DRFe_s[idx_2(i,j)] += dt * gdot_s * FMn[idx_2(i,k)] * n[idx_2(k,j)];
        for (int l = 0; l < dim; l++) {
          DRFe_Fe[idx_4(i,j,k,l)] += 2 * eye[idx_2(i,k)] * eye[idx_2(l,j)] - eye[idx_2(i,l)] * eye[idx_2(k,j)];;
          for (int p = 0; p < dim; p++) {
            DRFe_Fe[idx_4(i,j,k,l)] += dt * FMn[idx_2(i,p)] * (gdot_Fe[idx_2(k,l)] * n[idx_2(p,j)]
                                                               + gdot * n_Fe[idx_4(p,j,k,l)]);
          }
        }
      }
    }
  }
  *DRs_s = (1 - dt * sdot_s) / param_s_ss;
  return err;
}
