/* HEADER */
/**
 * This file defines the implementation for finite strain J2
 * plasticity with damage employing the stress split. This
 * formulations was chosen as it allows for simple computation of the
 * exact algorithmic tangent for the coupled plasticity and
 * damage. This is because the evaluation of the damage parameter is
 * not a function of any of the plastic variables.
 *
 * REFERENCES:
 *
 * Simo, J. C., and J. W. Ju. "On continuum damage-elastoplasticity at
 * finite strains." Computational Mechanics 5.5 (1989): 375-400.
 *
 * Mosby, Matthew and K. Matous. "On mechanics and material length scales of failure
 * in heterogeneous interfaces using a finite strain high performance
 * solver." Modelling and Simulation in Materials Science and
 * Engineering 23.8 (2015): 085014.
 *
 * Authors:
 *  Matt Mosby, University of Notre Dame, <mmosby1@nd.edu>
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "allocation.h"
#include "cm_j2_plasticity.h"
#include "cm_iso_viscous_damage.h"
#include "constitutive_model.h"
#include "index_macros.h"
#include "utils.h"
#include "new_potentials.h"
#include <math.h>
#include <string.h> 
#include <assert.h>

/* Define constant dimensions. Note cannot use `static const` with
   initialization list */
namespace {
constexpr int                dim = 3;
constexpr int             tensor = 9;
constexpr int            tensor4 = 81;
constexpr double   DAMAGE_THRESH = 0.9999;
constexpr double j2d_int_alg_tol = 1.0e-10;
//static const double eye[tensor] = {[0] = 1.0, [4] = 1.0, [8] = 1.0};
constexpr  double eye[tensor] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

/* macros for easy access to Constitutive_model structure */
#define cm_Fs(m) ((m)->vars_list[0][m->model_id].Fs)
#define cm_Fs_data(m, idx) ((m)->vars_list[0][m->model_id].Fs[(idx)].m_pdata)
#define cm_vars(m) ((m)->vars_list[0][m->model_id].state_vars->m_pdata)
#define cm_flags(m) ((m)->vars_list[0][m->model_id].flags)
#define cm_func(m) ((m)->param)
#define cm_hmat(m) ((m)->param->p_hmat)
#define cm_param(m) ((m)->param->model_param)

/* private context structure */
struct j2d_ctx {
  double F[tensor];
  double dt;
};

enum {FN, FNP1, SPN, SP, NUM_Fs};
enum {epn, ep, gam_n, gam, wn, w, Xn, X, Hn, H, NUM_vars};
enum {damaged_n, damaged, NUM_flags};
enum {G, nu, beta, hp, k0, mu, ome_max, p1, p2, Yin, NUM_param};
}

int CM_J2P_PARAM::get_var_info(Model_var_info &info)
const 
{
  int err = 0;

  info.n_Fs = NUM_Fs;
  info.n_vars = NUM_vars;
  info.n_flags = NUM_flags;
  info.F_names = PGFEM_malloc<char*>(NUM_Fs);
  info.var_names = PGFEM_malloc<char*>(NUM_vars);
  info.flag_names = PGFEM_malloc<char*>(NUM_flags);

  /* allocate/copy strings */
  info.F_names[FNP1] = strdup("F");
  info.F_names[FN]   = strdup("Fn");
  info.F_names[SPN]   = strdup("spn");
  info.F_names[SP]   = strdup("sp");

  info.var_names[epn] = strdup("epn");
  info.var_names[ep] = strdup("ep");
  info.var_names[gam_n] = strdup("gam_n");
  info.var_names[gam] = strdup("gam");

  info.var_names[wn] = strdup("wn");
  info.var_names[w] = strdup("w");
  info.var_names[Xn] = strdup("Xn");
  info.var_names[X] = strdup("X");
  info.var_names[Hn] = strdup("Hn");
  info.var_names[H] = strdup("H");

  info.flag_names[damaged_n] = strdup("damaged_n");
  info.flag_names[damaged] = strdup("damaged");

  return err;
}

int CM_J2P_PARAM::set_init_vals(Constitutive_model *m)
const
{
  int err = 0;
  /* tensors default init to identity and vars default to 0. Only need
     to ensure that the stress tensors (spn, sp) are all zeros */
  memset(cm_Fs_data(m, SPN), 0, tensor * sizeof(double));
  memset(cm_Fs_data(m, SP), 0, tensor * sizeof(double));
  return err;
}

int CM_J2P_PARAM::reset_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  double *vars = cm_vars(m);
  int *flags = cm_flags(m);
  Matrix<double> *Fs = cm_Fs(m);  
  Fs[FNP1]= Fs[FN];
  Fs[SP]  = Fs[SPN];
  vars[ep] = vars[epn];
  vars[gam] = vars[gam_n];
  vars[w] = vars[wn];
  vars[X] = vars[Xn];
  vars[H] = vars[Hn];
  flags[damaged] = flags[damaged_n];
  return err;
}
int CM_J2P_PARAM::update_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  double *vars = cm_vars(m);
  int *flags = cm_flags(m);
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  Fs[FN] = Fs[FNP1];
  Fs[SPN]= Fs[SP];
  vars[epn] = vars[ep];
  vars[gam_n] = vars[gam];
  vars[wn] = vars[w];
  vars[Xn] = vars[X];
  vars[Hn] = vars[H];
  flags[damaged_n] = flags[damaged];
  return err;
}

int CM_J2P_PARAM::read_param(FILE *in)
const
{
  int err = 0;
  /* get pointer to parameter data */
  double *param = this->model_param;
  assert(param != NULL); // check the pointer

  /* There are three (3) sets of parameters:

     1. elastic properties (G, nu) **NOTE**: These are currently read
        from the HOMMAT object.
     2. plastic properties (hp, beta, k0)
     3. damage properties (mu, p1, p2, Yin)

     Each set of parameters is read in alphabetical order */

  /* ELASTIC PROPERTIES */
  /* err += scan_for_valid_line(in); */
  /* int match = fscanf(in, "%lf %lf", */
  /*                    param + G, param + nu); */

  /*
    Get the elastic properties from the HOMMAT object.

    This is to maintain consistency with the other implementations and
    to ensure that the correct bulk modulous is computed elsewhere.
    Note that we **DO NOT** use the deviatoric constitutive law
    prescribed in the HOMMAT object, but use a spatial formulation of
    the Neo-Hookean model to maintain consistency with the J2+damage
    formulation in Simo and Ju.
  */
  int match = 2;
  param[G] = this->p_hmat->G;
  param[nu] = this->p_hmat->nu;

  /* PLASTIC PROPERTIES */
  err += scan_for_valid_line(in);
  match += fscanf(in, "%lf %lf %lf",
                  param + beta, param + hp, param + k0);

  err += scan_for_valid_line(in);

  /* DAMAGE PROPERTIES */
  match += fscanf(in, "%lf %lf %lf %lf %lf",
                  param + mu, param + ome_max, param + p1,
                  param + p2, param + Yin);

  if (match != NUM_param) err++;
  assert(match == NUM_param && "Did not read expected number of parameters");

  /* ome_max in [0, 1) */
  if (param[ome_max] >= 1) param[ome_max] = DAMAGE_THRESH;

  /* scan past any other comment/blank lines in the block */
  err += scan_for_valid_line(in);

  /* not expecting EOF, check and return error if encountered */
  if (feof(in)) err ++;
  assert(!feof(in) && "EOF reached prematurely");

  return err;
}

int CM_J2P_PARAM::destroy_ctx(void **ctx)
const
{
  int err = 0;
  free(*ctx);
  *ctx = NULL;
  return err;
}

/* bbar = J^-(2/3) F F' */
static void j2d_bbar(const double * __restrict F,
                     double * __restrict bbar)
{
  memset(bbar, 0, tensor * sizeof(*bbar));
  const double J23 = pow(det3x3(F), -2./3.);
  for (int i = 0; i < dim; i++){
    for (int j = 0; j < dim; j++) {
      for(int k = 0; k < dim; k++) {
        bbar[idx_2(i,j)] += J23 * F[idx_2(i,k)] * F[idx_2(j,k)];
      }
    }
  }
}

/* dev(a) = a - 1/3 tr(a) i */
static void j2d_dev(const double * __restrict a,
                    double * __restrict dev_a)
{
  const double tra = 1./3. * (a[0] + a[4] + a[8]);
  memcpy(dev_a, a, tensor * sizeof(*a));
  dev_a[0] -= tra;
  dev_a[4] -= tra;
  dev_a[8] -= tra;
}

/* s0 = G * dev(bbar) */
static void j2d_compute_s0(const double G,
                           const double * __restrict bbar,
                           double * __restrict s0)
{
  double devbbar[tensor] = {0};
  j2d_dev(bbar, devbbar);
  for (int i = 0; i < tensor; i++) {
    s0[i] = G * devbbar[i];
  }
}

/* compute Fubar */
static int j2d_compute_Fubar(const double * __restrict F,
                             const double * __restrict Fn,
                             double * __restrict Fubar)
{
  int err = 0;
  double FnI[tensor] = {0};
  err += inv3x3(Fn, FnI);
  for (int i = 0; i < dim; i++){
    for (int j = 0; j < dim; j++) {
      for(int k = 0; k < dim; k++) {
        Fubar[idx_2(i,j)] += F[idx_2(i,k)] * FnI[idx_2(k,j)];
      }
    }
  }
  const double Ju13 = pow(det3x3(Fubar), -1./3.);
  for (int i = 0; i < tensor; i++) Fubar[i] *= Ju13;

  return err;
}

/* push-forward operation a = F A F'. */
static void j2d_push_forward(const double * __restrict F,
                             const double * __restrict A,
                             double * __restrict a)
{
  memset(a, 0, tensor * sizeof(*a));
  double wkspc[tensor] = {0};
  for (int i = 0; i < dim; i++){
    for (int j = 0; j < dim; j++) {
      for(int k = 0; k < dim; k++) {
        wkspc[idx_2(i,j)] += F[idx_2(i,k)] * A[idx_2(k,j)];
      }
    }
  }
  for (int i = 0; i < dim; i++){
    for (int j = 0; j < dim; j++) {
      for(int k = 0; k < dim; k++) {
        a[idx_2(i,j)] += wkspc[idx_2(i,k)] * F[idx_2(j,k)];
      }
    }
  }
}

static int j2d_pull_back(const double *F,
                         const double *a,
                         double *A)
{
  double FI[tensor] = {0};
  int err = inv3x3(F, FI);
  j2d_push_forward(FI, a, A);
  return err;
}

/* compute the trial plastic stress */
static int j2d_compute_sp_tr(const double *F,
                             const double *Fn,
                             const double *spn,
                             double *sp_tr)
{
  int err = 0;
  double Fubar[tensor] = {0};
  err += j2d_compute_Fubar(F, Fn, Fubar);
  j2d_push_forward(Fubar, spn, sp_tr);
  const double tr = (sp_tr[0] + sp_tr[4] + sp_tr[8]) / 3.0;
  sp_tr[0] -= tr;
  sp_tr[4] -= tr;
  sp_tr[8] -= tr;

  return err;
}

static double j2d_compute_normal(const double * __restrict s_tr,
                                 const double * __restrict sp_tr,
                                 const double * __restrict param,
                                 double  * __restrict n)
{
  const double coef = param[hp] / (3. * param[G]) * (1. - param[beta]);
  double nrm = 0;

  /* compute ksi_tr and ||ksi_tr|| simultaneously */
  for(int i = 0; i < tensor; i++){
    n[i] = s_tr[i] - coef * sp_tr[i];
    nrm += n[i] * n[i];
  }
  nrm = sqrt(nrm);

  /* compute normal: n = ksi_tr / ||ksi_tr|| */
  if(nrm > 0) {
    for(int i = 0; i < tensor; i++) n[i] /= nrm;
  }
  return nrm;
}

/* Y = G / 2 *(tr(bbar) - 3) + k U */
static int j2d_compute_Y0(const HOMMAT *p_hmat,
                          const double *bbar,
                          const double J,
                          const double G,
                          double *Y0)
{
  const double kappa = hommat_get_kappa(p_hmat);
  double U = 0.0;
  new_pot_compute_U(J, p_hmat, &U);
  *Y0 = 0.5 * G * (bbar[0] + bbar[4] + bbar[8] - 3.0) + kappa * U;
  return 0;
}

/* see box 4 of Simo and Ju (1989) */
int CM_J2P_PARAM::integration_algorithm(Constitutive_model *m,
                                        const void *CTX)
const
{
  int err = 0;
  auto ctx = (j2d_ctx *) CTX;

  /* store the current (n + 1) deformation gradient in the CM object. */
  memcpy(cm_Fs_data(m, FNP1), ctx->F, tensor * sizeof(*(ctx->F)));

  /* get pointers to tensor variables */
  const double *param = cm_param(m);
  const double *F = cm_Fs_data(m, FNP1);
  const double *Fn = cm_Fs_data(m, FN);
  const double *spn = cm_Fs_data(m, SPN);
  double *sp = cm_Fs_data(m, SP);

  /* get pointer to scalar state variables */
  double *vars = cm_vars(m);

  /* compute bbar at n + 1 */
  double bbar[tensor] = {0};
  j2d_bbar(F, bbar);

  /* compute mu_bar */
  const double J = det3x3(F);
  const double J23 = pow(J, -2./3.);
  const double mu_bar = param[G] * J23 * (bbar[0] + bbar[4] + bbar[8]) / 3;

  /* compute sp_tr */
  double sp_tr[tensor] = {0};
  j2d_compute_sp_tr(F, Fn, spn, sp_tr);

  /* compute s_tr */
  double s0[tensor] = {0};
  double s_tr[tensor] = {0};
  j2d_compute_s0(param[G], bbar, s0);
  for (int i = 0; i < tensor; i++) s_tr[i] = s0[i] - sp_tr[i];

  /* compute ksi_tr and the normal of plastic loading */
  double n[tensor] = {0};
  const double ksi_nrm = j2d_compute_normal(s_tr, sp_tr, param, n);

  /* yield function */
  const double phi = ksi_nrm - sqrt(2./3.) * (param[k0] + param[beta] * param[hp] * vars[epn]);

  if (phi <= j2d_int_alg_tol) {
    vars[gam] = 0.0;
    vars[ep] = vars[epn];
    memcpy(sp, sp_tr, tensor * sizeof(*sp));
  } else {
    vars[gam] = phi / (2. * mu_bar
                         * (1. + param[hp] / (3. * param[G]) * (1. - param[beta])
                            + param[beta] * param[hp] / (3. * mu_bar) )
                         );
    const double tmp = 2 * mu_bar * vars[gam];
    vars[ep] = vars[epn] + sqrt(2./3.) * vars[gam];
    for (int i = 0; i < tensor; i++) sp[i] = sp_tr[i] + tmp * n[i];
  }

  /* damage integration algorithm */
  double Y0 = 0.0;
  err += j2d_compute_Y0(cm_hmat(m), bbar, J, param[G], &Y0);
  err += ivd_public_int_alg(&vars[w], &vars[X], &vars[H], &cm_flags(m)[damaged],
                            vars[wn], vars[Xn], ctx->dt, Y0, param[mu],
                            param[ome_max], param[p1], param[p2], param[Yin]);
  return err;
}

/* compute the deviatoric initial/unloading tangent in the reference
   configuration */
static int j2d_unloading_Aep_dev(const Constitutive_model *m,
                                 const void *CTX,
                                 double * __restrict Aep_dev)
{
  int err = 0;
  memset(Aep_dev, 0, tensor4 * sizeof(*Aep_dev));

  //const j2d_ctx *ctx = CTX;
  const double *F = cm_Fs_data(m,FNP1);
  const double *Fn = cm_Fs_data(m, FN);
  const double *spn = cm_Fs_data(m, SPN);

  /* mu is shear modulus (G) in this context!!! */
  const double mu = cm_param(m)[G];

  /* compute C and related terms */
  double C[tensor] = {0};
  double CI[tensor] = {0};
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
              dim, dim, dim, 1.0, F, dim, F, dim, 0.0, C, dim);
  err += inv3x3(C,CI);
  const double Cpp = C[0] + C[4] + C[8];
  const double J23 = pow(det3x3(C), -1./3.);

  /* compute pull-back of spn */
  double Spn[tensor] = {0};
  err += j2d_pull_back(Fn, spn, Spn);

  double CSp = 0.0;
  for (int i = 0; i < tensor; i++) CSp += C[i] * Spn[i];

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      const int ij = idx_2(i,j);
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          const int ijkl = idx_4(i,j,k,l);
          const int kl = idx_2(k,l);
          Aep_dev[ijkl] = ((2./3. * J23 * (mu * Cpp - CSp)
                            * (CI[idx_2(i,k)] * CI[idx_2(j,l)] + CI[ij] * CI[kl] / 3.))
                           - 2./3. * J23 *(CI[ij] * (mu * eye[kl] - Spn[kl])
                                           + (mu * eye[ij] - Spn[ij]) * CI[kl]));
        }
      }
    }
  }

  return err;
}

static int j2d_pull_back4(const double * __restrict FI,
                          const double * __restrict aep,
                          double * __restrict Aep)
{
  int err = 0;
  memset(Aep, 0, tensor4 * sizeof(*Aep));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        for (int l = 0; l < dim; l++) {
          const int ijkl = idx_4(i,j,k,l);
          for (int m = 0; m < dim; m++) {
            for (int n = 0; n < dim; n++) {
              for (int o = 0; o < dim; o++) {
                for (int p = 0; p < dim; p++) {
                  const int mnop = idx_4(m,n,o,p);
                  Aep[ijkl] += (FI[idx_2(i,m)] * FI[idx_2(j,n)]
                                * FI[idx_2(k,o)] * FI[idx_2(l,p)]
                                * aep[mnop]);
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

/* compute the deviatoric plastic loading tangent in the reference
   configuration */
static int j2d_loading_Aep_dev(const Constitutive_model *m,
                               const void *CTX,
                               double * __restrict Aep_dev)
{
  int err = 0;
  const auto ctx = (j2d_ctx *) CTX;
  const double *param = cm_param(m);
  const double *Fn = cm_Fs_data(m, FN);
  const double *spn = cm_Fs_data(m, SPN);

  /* get pointer to scalar state variables */
  const double *vars = cm_vars(m);

  /* compute bbar at n + 1 */
  const double J23 = pow(det3x3(ctx->F), -2./3.);
  double bbar[tensor] = {0};
  double devbbar[tensor] = {0};
  double FI[tensor] = {0};
  j2d_bbar(ctx->F, bbar);
  j2d_dev(bbar, devbbar);
  err += inv3x3(ctx->F, FI);

  /* compute Itr = G tr(bbar) - tr(Fubar spn Fubar') */
  double sp_tr[tensor] = {0};
  double Itr = 0;
  {
    double Fubar[tensor] = {0};
    err += j2d_compute_Fubar(ctx->F, Fn, Fubar);
    j2d_push_forward(Fubar, spn, sp_tr);
    const double trace_sp = sp_tr[0] + sp_tr[4] + sp_tr[8];
    sp_tr[0] -= trace_sp / 3.0;
    sp_tr[4] -= trace_sp / 3.0;
    sp_tr[8] -= trace_sp / 3.0;
    Itr = param[G] * (bbar[0] + bbar[4] + bbar[8]) - trace_sp;
  }

  /* compute s_tr, normal and ||s_tr|| */
  double s_tr[tensor] = {0};
  double normal[tensor] = {0};
  double normal2[tensor] = {0};
  j2d_compute_s0(param[G], bbar, s_tr);
  double norm_s_tr = 0.0;
  for (int i = 0; i < tensor; i++) {
    s_tr[i] -= sp_tr[i];
    norm_s_tr += s_tr[i] * s_tr[i];
  }
  norm_s_tr = sqrt(norm_s_tr);

  j2d_compute_normal(s_tr, sp_tr, param, normal);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim,
              1.0, normal, dim, normal, dim, 0.0, normal2, dim);

  /* compute factors */
  const double mu_bar = param[G] * J23 * (bbar[0] + bbar[4] + bbar[8]) / 3;
  const double f0 = 1.0 - 2.0 * mu_bar * vars[gam] / norm_s_tr;
  const double del0 = 1.0 + param[hp] / (3.0 * mu_bar);
  const double f1 = (1.0 / del0 - 1.0 + f0);
  const double del1 = f1 * 2./3. * Itr;
  const double del2 = (1.0 / del0 - 1.0) * 4./3. * mu_bar * vars[gam] - 2./3. * norm_s_tr * f1;
  const double del3 = 2.0 * norm_s_tr * f1;
  const double del4 = (1.0 / del0  - 1.0) * 4./3. * param[G] * vars[gam] * J23;

  /* compute aep according to box 5 in Simo and Ju 1989*/
  double aep[tensor4] = {0};
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      const int ij = idx_2(i,j);
      for (int k = 0; k < dim; k++) {
        const int ik = idx_2(i,k);
        for (int l = 0; l < dim; l++) {
          const int lj = idx_2(l,j);
          const int kl = idx_2(k,l);
          const int ijkl = idx_4(i,j,k,l);
          aep[ijkl] = (f0 * (2./3. * Itr * (eye[ik] * eye[lj] - eye[ij] * eye[kl] / 3.)
                             - 2./3. * (s_tr[ij] * eye[kl] + eye[ij] * s_tr[kl]))
                       - del1 * normal[ij] * normal[kl]
                       - del2 * (normal[ij] * eye[kl] + normal[kl] * eye[ij])* 0.5
                       - del3 * (normal[ij] * normal2[kl] + normal[kl] * normal2[ij])* 0.5
                       - del4 * (normal[ij] * devbbar[kl] + normal[kl] * devbbar[ij])* 0.5);
        }
      }
    }
  }

  err += j2d_pull_back4(FI, aep, Aep_dev);
  return err;
}

static int j2d_compute_Lbar(const Constitutive_model *m,
                            const void *CTX,
                            double * __restrict Lbar)
{
  int err = 0;
  const auto ctx = (j2d_ctx *) CTX;
  const double kappa = hommat_get_kappa(cm_hmat(m));
  const double J = det3x3(ctx->F);
  double C[tensor] = {0};
  double C_I[tensor] = {0};
  double dudj = 0.0;
  double d2udj2 = 0.0;

  /* compute C, CI */
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
              dim, dim, dim, 1.0, ctx->F, dim, ctx->F, dim, 0.0, C, dim);
  err += inv3x3(C, C_I);

  /* compute deviatoric tangent */
  if (cm_vars(m)[gam] > 0) err += j2d_loading_Aep_dev(m, CTX, Lbar);
  else err += j2d_unloading_Aep_dev(m, CTX, Lbar);

  new_pot_compute_dudj(J, cm_hmat(m), &dudj);
  new_pot_compute_d2udj2(J, cm_hmat(m), &d2udj2);

  for (int i=0; i < dim; i++) {
    for (int j=0; j < dim; j++) {
      for (int k=0; k < dim; k++) {
    for (int l=0; l < dim; l++) {
      const int idx4 = idx_4(i,j,k,l);
      /* Deviatoric + Volumetric stiffness */
      Lbar[idx4] += ((kappa * J * (dudj + J * d2udj2)
                          * C_I[idx_2(i,j)] * C_I[idx_2(k,l)])
             - (2. * kappa * J * dudj
                            * C_I[idx_2(i,k)] * C_I[idx_2(l,j)])
                         );
    }
      }
    }
  }

  return err;
}

static int j2d_compute_S0_Sbar(const Constitutive_model *m,
                               const void *CTX,
                               double *S0,
                               double *Sbar)
{
  int err = 0;
  auto ctx = (j2d_ctx *) CTX;

  /* compute the current configuration deviatoric stresses */
  double s0[tensor] = {0};
  double sbar[tensor] = {0};
  double bbar[tensor] = {0};
  j2d_bbar(ctx->F, bbar);
  j2d_compute_s0(cm_param(m)[G], bbar, s0);
  for (int i = 0; i < tensor; i++) sbar[i] = s0[i] - cm_Fs_data(m, SP)[i];

  /* perform pull-back */
  err += j2d_pull_back(ctx->F, s0, S0);
  err += j2d_pull_back(ctx->F, sbar, Sbar);

  /* compute volumetric stress */
  double C[tensor] = {0};
  double CI[tensor] = {0};
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim, dim, dim,
              1.0, ctx->F, dim, ctx->F, dim, 0.0, C, dim);
  err += inv3x3(C, CI);
  const double J = det3x3(ctx->F);
  const double kappa = hommat_get_kappa(cm_hmat(m));
  double dudj = 0.0;
  new_pot_compute_dudj(J, cm_hmat(m), &dudj);
  for (int i = 0; i < tensor; i++) {
    S0[i] += kappa * J * dudj * CI[i];
    Sbar[i] += kappa * J * dudj * CI[i];
  }

  return err;
}

static int j2d_modify_AST(const Constitutive_model *m,
                          const void *CTX,
                          double *L)
{
  int err = 0;
  auto ctx = (j2d_ctx *) CTX;
  const double dmu = ctx->dt * cm_param(m)[mu];
  const double evo = dmu * cm_vars(m)[H] / (1.0 + dmu);
  double S0[tensor] = {0};
  double Sbar[tensor] = {0};
  err += j2d_compute_S0_Sbar(m, CTX, S0, Sbar);

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
    for (int l = 0; l < dim; l++) {
      const int idx4 = idx_4(i,j,k,l);
      L[idx4] -= evo * 0.5 * (Sbar[idx_2(i,j)] * S0[idx_2(k,l)]
                                  + S0[idx_2(i,j)] * Sbar[idx_2(k,l)]);
        }
      }
    }
  }

  return err;
}

static int j2d_compute_AST(const Constitutive_model *m,
                           const void *ctx,
                           double *L)
{
  int err = 0;
  err += j2d_compute_Lbar(m, ctx, L);

  /* scale by damage parameter */
  const double dam = 1.0 - cm_vars(m)[w];
  for (int i = 0; i < tensor4; i++) {
    L[i] *= dam;
  }

  /* if evolving damage, modify tangent */
  if (cm_flags(m)[damaged]) {
    err += j2d_modify_AST(m, ctx, L);
  }

  return err;
}

/* compute the deviatoric part of the effective (undamaged) Kirckhoff stress */
static int j2d_compute_sbar(const double *F,
                            const double *sp,
                            const double G,
                            double * __restrict sbar)
{
  int err = 0;
  double bbar[tensor] = {0};
  j2d_bbar(F, bbar);
  j2d_compute_s0(G, bbar, sbar);
  for (int i = 0; i < tensor; i++) sbar[i] -= sp[i];
  return err;
}

/* compute the deviatoric PK2 stress */
int CM_J2P_PARAM::compute_dev_stress(const Constitutive_model *m,
                                     const void *CTX,
                                     double *stress)
const 
{
  int err = 0;
  auto ctx = (j2d_ctx *) CTX;
  double sbar[tensor] = {0};
  err += j2d_compute_sbar(ctx->F,
                          cm_Fs_data(m, SP),
                          cm_param(m)[G],
                          sbar);
  err += j2d_pull_back(ctx->F, sbar, stress);
  const double dam = 1. - cm_vars(m)[w];
  for(int i = 0; i < tensor; i++) stress[i] *= dam;
  return err;
}

/* compute the volumetric stress. NOTE we are still using the actual
   HOMMAT object here for convenience, but the only used information
   is for the flag */
int CM_J2P_PARAM::compute_dudj(const Constitutive_model *m,
                               const void *ctx,
                               double *dudj)
const
{
  int err = 0;
  auto CTX = (j2d_ctx *) ctx;
  double J = det3x3(CTX->F);
  new_pot_compute_dudj(J, cm_hmat(m), dudj);

  /* scale by damage variable */
  *dudj *= (1.0 - cm_vars(m)[w]);
  return err;
}

int CM_J2P_PARAM::get_hardening(const Constitutive_model *m,
                                double *var,
                                const int stepno)
const 
{
  int err = 0;
  double *s_var = cm_vars(m);
  switch(stepno)
  {
    case 0: // n-1
      *var = s_var[wn];
      break;
    case 1: // n
      *var = s_var[wn];
      break;
    case 2: // n+1
      *var = s_var[w];
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;  
}

int CM_J2P_PARAM::get_plast_strain_var(const Constitutive_model *m,
                                       double *ep)
const
{
  int err = 0;
  *ep = cm_vars(m)[epn];
  return err;
}

static int j2d_identity_tensor(const Constitutive_model *m,
                               double *F)
{
  int err = 0;
  double I[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
  memcpy(F,I,tensor*sizeof(double));
  return err;
}

int CM_J2P_PARAM::get_pF(const Constitutive_model *m,
                         double *F,
                         const int stepno)
const
{
  return j2d_identity_tensor(m, F);
}

int CM_J2P_PARAM::get_F(const Constitutive_model *m,
                        double *F,
                        const int stepno)
const 
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  switch(stepno)
  {
    case 0: // n-1
      memcpy(F,Fs[FN].m_pdata, tensor*sizeof(double));
      break;
    case 1: // n
      memcpy(F,Fs[FN].m_pdata, tensor*sizeof(double));
      break;
    case 2: // n+1
      memcpy(F,Fs[FNP1].m_pdata, tensor*sizeof(double));
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;
}

int CM_J2P_PARAM::get_eF(const Constitutive_model *m,
                         double *eF_in,
                         const int stepno)
const
{
  return this->get_F(m,eF_in,stepno);
}

static int j2d_read_tensor(FILE *in,
                           double *FF)
{
  return (fscanf(in, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &FF[0], &FF[1], &FF[2], &FF[3],
                 &FF[4], &FF[5], &FF[6], &FF[7], &FF[8]) != tensor);
}

int CM_J2P_PARAM::read_restart(FILE *in, Constitutive_model *m)
const
{
  int err = 0;
  err += j2d_read_tensor(in, cm_Fs_data(m, FN));
  err += j2d_read_tensor(in, cm_Fs_data(m, SPN));
  double *vars = cm_vars(m);
  int *flags = cm_flags(m);
  if( fscanf(in, "%lf %lf %lf %lf %lf %d",
             &vars[epn], &vars[gam_n],
             &vars[wn], &vars[Xn], &vars[Hn],
             &flags[damaged_n]) != 6 ) err ++;
  err += this->reset_state_vars(m);
  return err;
}

static int j2d_write_tensor(FILE *out,
                            const double *FF)
{
  return (fprintf(out, "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
                  FF[0], FF[1], FF[2], FF[3], FF[4], FF[5], FF[6], FF[7], FF[8]) < 0);
}

int CM_J2P_PARAM::write_restart(FILE *out, const Constitutive_model *m)
const
{
  int err = 0;
  err += j2d_write_tensor(out, cm_Fs_data(m, FN));
  err += j2d_write_tensor(out, cm_Fs_data(m, SPN));
  const double *vars = cm_vars(m);
  const int *flags = cm_flags(m);
  err += (fprintf(out, "%.17e %.17e %.17e %.17e %.17e %d\n",
                  vars[epn], vars[gam_n],
                  vars[wn], vars[Xn], vars[Hn],
                  flags[damaged_n]) < 0);
  return err;
}

int CM_J2P_PARAM::get_subdiv_param(const Constitutive_model *m,
                                   double *subdiv_param,
                                   double dt)
const
{
  int err = 0;

  /* compute the plastic subdivision parameter */
  double plast_param = 0.0;
  /* no subdivision scheme yet... */

  /* compute the damage subdivision parameter */
  double damage_param = 0.0;
  err += ivd_public_subdiv_param(cm_vars(m)[wn], cm_vars(m)[w], &damage_param);

  /* return the maximum subdivision parameter */
  *subdiv_param = (plast_param >= damage_param)? plast_param : damage_param;

  return err;
}

int CM_J2P_PARAM::update_elasticity(const Constitutive_model *m,
                                    const void *ctx,
                                    double *L,
                                    double *S,
                                    const int compute_stiffness)
const
{
  int err = 0;

  double S0[tensor] = {0};
  err += j2d_compute_S0_Sbar(m,ctx,S0,S);

  double dam = (1.0 - cm_vars(m)[w]);
  for(int a=0; a<tensor; a++)
    S[a] *= dam;

  if(compute_stiffness)
    err += j2d_compute_AST(m, ctx, L); //compute stiffness

  return err;
}

int CM_J2P_PARAM::model_dependent_initialization(void)
{
  int err = 0;

  /* reset counters/flags */
  this->type = J2_PLASTICITY_DAMAGE;

  /* allocate room for parameters */
  this->n_param = NUM_param;
  this->model_param = PGFEM_calloc(double, NUM_param);

  return err;
}

int j2d_plasticity_model_ctx_build(void **ctx,
                                   const double *F,
                                   const double dt)
{
  j2d_ctx *CTX = PGFEM_malloc<j2d_ctx>();
  //CTX = (j2d_ctx *) malloc(sizeof(*CTX));
  memcpy(CTX->F, F, tensor * sizeof(*F));
  CTX->dt = dt;
  *ctx = CTX;
  return 0;
}
