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
static const int _n_Fs = 4;
static const int _n_vars = 4;
enum {_M,_W,_M_n,_W_n};
enum {_s,_lam,_s_n,_lam_n};
static const double eye[tensor] = {1.0,0,0, 0,1.0,0, 0,0,1.0};

/* material parameters */
static double param_A;
static double param_T;
static double param_N;
static double param_Cr;
static double param_alpha;
static double param_gdot0;
static double param_h;
static double param_s_ss;

/*
 * Purely static functions
 */
static double compute_bulk_mod(const HOMMAT *mat)
{
  return ( (2* mat->G * (1 + mat->nu)) / (3 * (1 - 2 * mat->nu)) );
}

static void compute_Fe(const double * restrict F,
                       const double * restrict M,
                       double * restrict Fe)
{
  memset(Fe, 0, tensor * dim * sizeof(*Fe));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
       for (int k = 0; k < dim; k++) {
        Fe[idx_2(i,j)] = F[idx_2(i,k)] * M[idx_2(k,j)];
      }
    }
  }
}

static double compute_plam(const double *Cp)
{
  return sqrt(1./3. * (Cp[0] + Cp[4] + Cp[8]));
}

static void compute_Cp(double * restrict Cp,
                       const double * restrict Fp)
{
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        Cp[idx_2(i,j)] = Fp[idx_2(i,k)] * Fp[idx_2(j,k)];
      }
    }
  }
}

static void compute_Cpdev(double * restrict Cpdev,
                          const double * restrict Cp)
{
  const double Jdev = pow(det3x3(Cp),-1./3.);
  for (int i = 0; i < tensor; i++) {
    Cpdev[i] = Jdev * Cp[i];
  }
}

static void compute_Fe_Ce(const double * restrict F,
                          const double * restrict M,
                          double * restrict Fe,
                          double * restrict Ce)
{
  compute_Fe(F,M,Fe);
  memset(Ce, 0, tensor * sizeof(*Ce));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        Ce[idx_2(i,j)] = Fe[idx_2(k,i)] * Fe[idx_2(k,j)];
      }
    }
  }
}

static void compute_Sdev(const double *Ce,
                         const HOMMAT *p_hmat,
                         double *Sdev)
{
  devStressFuncPtr Stress = getDevStressFunc(-1,p_hmat);
  Stress(Ce,p_hmat,Sdev);
}

static void compute_Ldev(const double *Ce,
                         const HOMMAT *p_hmat,
                         double *Ldev)
{
  matStiffFuncPtr Tangent = getMatStiffFunc(-1,p_hmat);
  Tangent(Ce,p_hmat,Ldev);
}

static void compute_dudj(const double Je,
                         const HOMMAT *p_hmat,
                         double *dudj)
{
  dUdJFuncPtr Pressure = getDUdJFunc(-1,p_hmat);
  Pressure(Je,p_hmat,dudj);
}

static void compute_d2udj2(const double Je,
                           const HOMMAT *p_hmat,
                           double *d2udj2)
{
  d2UdJ2FuncPtr D_Pressure = getD2UdJ2Func(-1,p_hmat);
  D_Pressure(Je,p_hmat,d2udj2);
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

  /* compute the deformation */
  double Fe[tensor], Ce[tensor];
  compute_Fe_Ce(CTX->F,m->vars.Fs[_M].m_pdata,Fe,Ce);
  double Fp[tensor];
  inv3x3(m->vars.Fs[_M].m_pdata,Fp);
  double Je = det3x3(Fe);

  /* compute the elastic stress */
  double Sdev[tensor];
  compute_Sdev(Ce,p_hmat,Sdev);

  /* compute the pressure */
  const double kappa = compute_bulk_mod(m->param->p_hmat);
  double pressure = 0.0;
  compute_dudj(Je,p_hmat,&pressure);
  pressure *= kappa / 2.0;

  /* compute the pressure-dependent athermal shear stress */
  double s_s = m->vars.state_vars->m_pdata[_s] + param_alpha * pressure;

  /* compute the plastic backstress */
  double Bdev[tensor];
  err += BPA_compute_Bdev(Bdev,Fp);

  /* compute the loading direction */
  double normal[tensor], eq_sig_dev[tensor];
  double tau = 0;
  err += BPA_compute_loading_dir(normal,eq_sig_dev,&tau,Sdev,Bdev,Fe);

  /* compute gdot */
  double gdot = 0.0;
  err += BPA_compute_gdot(&gdot,tau,s_s);


  return err;
}

int BPA_dev_stress(const Constitutive_model *m,
                   const void *ctx,
                   Matrix_double *dev_stress)
{
  int err = 0;
  const BPA_ctx *CTX = ctx;
  double Fe[tensor], Ce[tensor];
  compute_Fe_Ce(CTX->F,m->vars.Fs[_M].m_pdata,Fe,Ce);
  compute_Sdev(Ce,m->param->p_hmat,dev_stress->m_pdata);
  return err;
}

int BPA_dudj(const Constitutive_model *m,
             const void *ctx,
             double *dudj)
{
  int err = 0;
  const BPA_ctx *CTX = ctx;
  double Fe[tensor];
  compute_Fe(CTX->F,m->vars.Fs[_M].m_pdata,Fe);
  double Je = det3x3(Fe);
  compute_dudj(Je,m->param->p_hmat,dudj);
  return err;
}

int BPA_dev_tangent(const Constitutive_model *m,
                    const void *ctx,
                    Matrix_double *dev_tangent)
{
  int err = 0;
  const BPA_ctx *CTX = ctx;
  double Fe[tensor], Ce[tensor];
  compute_Fe_Ce(CTX->F,m->vars.Fs[_M].m_pdata,Fe,Ce);
  compute_Ldev(Ce,m->param->p_hmat,dev_tangent->m_pdata);
  return err;
}

int BPA_d2udj2(const Constitutive_model *m,
               const void *ctx,
               double *d2udj2)
{
  int err = 0;
  const BPA_ctx *CTX = ctx;
  double Fe[tensor];
  compute_Fe(CTX->F,m->vars.Fs[_M].m_pdata,Fe);
  double Je = det3x3(Fe);
  compute_d2udj2(Je,m->param->p_hmat,d2udj2);
  return err;
}

int BPA_update_vars(Constitutive_model *m)
{
  int err = 0;
  Matrix_copy(m->vars.Fs[_M_n],m->vars.Fs[_M]);
  Matrix_copy(m->vars.Fs[_W_n],m->vars.Fs[_W]);

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
  double Cp[tensor];
  double Cpdev[tensor];
  compute_Cp(Cp,Fp);
  compute_Cpdev(Cpdev,Cp);

  /* compute the backstess coefficient */
  const double lam_p = compute_plam(Cp);
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
  double tmp[tensor], tmp2[tensor];
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
  cblas_daxpy(tensor,1.0 / (*tau), eq_sig_dev, 1, normal, 1);

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

int plasticity_model_BPA_ctx_build(void **ctx
                                   //...
                                   )
{
  int err = 0;
  BPA_ctx *t_ctx = malloc(sizeof(*t_ctx));
  *ctx = t_ctx;
  return err;
}

int plasticity_model_BPA_ctx_destroy(void **ctx)
{
  int err = 0;
  free(*ctx);
  *ctx = NULL;
  return err;
}
