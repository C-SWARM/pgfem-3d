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

#include "cm_j2_plasticity.h"
#include "constitutive_model.h"
#include <math.h>
#include <string.h>
#include "utils.h"
#include "index_macros.h"
#include "data_structure_c.h"
#include "new_potentials.h"

/* Define constant dimensions. Note cannot use `static const` with
   initialization list */
#define dim  3
#define tensor 9
#define tensor4 81

static const double j2d_int_alg_tol = 1.0e-10;

Define_Matrix(double);

/* macros for easy access to Constitutive_model structure */
#define cm_Fs(m) ((m)->vars.Fs)
#define cm_Fs_data(m, idx) ((m)->vars.Fs[(idx)].m_pdata)
#define cm_vars(m) ((m)->vars.state_vars->m_pdata)
#define cm_flags(m) ((m)->vars.flags)
#define cm_func(m) ((m)->param)
#define cm_hmat(m) ((m)->param->p_hmat)
#define cm_param(m) ((m)->param->model_param)

/* private context structure */
typedef struct {
  double F[tensor];
  double dt;
} j2d_ctx;

enum {FN, FNP1, SPN, SP, NUM_Fs};
enum {epn, ep, gam_n, gam, wn, w, Xn, X, Hn, H, NUM_vars};
enum {damaged_n, damaged, NUM_flags};
enum {G, nu, hp, beta, k0, mu, p1, p2, Yin, NUM_param};

static int j2d_get_info(Model_var_info **info)
{
  int err = 0;
  /* make sure I don't leak memory */
  if( *info != NULL) err += model_var_info_destroy(info);

  /* allocate pointers */
  (*info) = malloc(sizeof(**info));
  (*info)->n_Fs = NUM_Fs;
  (*info)->n_vars = NUM_vars;
  (*info)->n_flags = NUM_flags;
  (*info)->F_names = malloc(NUM_Fs * sizeof( ((*info)->F_names) ));
  (*info)->var_names = malloc( NUM_vars * sizeof( ((*info)->var_names) ));
  (*info)->flag_names = malloc( NUM_flags * sizeof( ((*info)->flag_names) ));

  /* allocate/copy strings */
  (*info)->F_names[FNP1] = strdup("F");
  (*info)->F_names[FN]   = strdup("Fn");
  (*info)->F_names[SPN]   = strdup("spn");
  (*info)->F_names[SP]   = strdup("sp");

  (*info)->var_names[epn] = strdup("epn");
  (*info)->var_names[ep] = strdup("ep");
  (*info)->var_names[gam_n] = strdup("gam_n");
  (*info)->var_names[gam] = strdup("gam");

  (*info)->var_names[wn] = strdup("wn");
  (*info)->var_names[w] = strdup("w");
  (*info)->var_names[Xn] = strdup("Xn");
  (*info)->var_names[X] = strdup("X");
  (*info)->var_names[Hn] = strdup("Hn");
  (*info)->var_names[H] = strdup("H");

  (*info)->flag_names[damaged_n] = strdup("damaged_n");
  (*info)->flag_names[damaged] = strdup("damaged");

  return err;
}

static int j2d_set_init_vals(Constitutive_model *m)
{
  int err = 0;
  /* tensors default init to identity and vars default to 0. Only need
     to ensure that the stress tensors (spn, sp) are all zeros */
  memset(cm_Fs_data(m, SPN), 0, tensor * sizeof(double));
  memset(cm_Fs_data(m, SP), 0, tensor * sizeof(double));
  return err;
}

static int j2d_reset(Constitutive_model *m)
{
  int err = 0;
  double *vars = cm_vars(m);
  int *flags = cm_flags(m);
  Matrix_copy(cm_Fs(m)[FNP1], cm_Fs(m)[FN]);
  Matrix_copy(cm_Fs(m)[SP], cm_Fs(m)[SPN]);
  vars[ep] = vars[epn];
  vars[gam] = vars[gam_n];
  vars[w] = vars[wn];
  vars[X] = vars[Xn];
  vars[H] = vars[Hn];
  flags[damaged] = flags[damaged_n];
  return err;
}

static int j2d_update(Constitutive_model *m)
{
  int err = 0;
  double *vars = cm_vars(m);
  int *flags = cm_flags(m);
  Matrix_copy(cm_Fs(m)[FN], cm_Fs(m)[FNP1]);
  Matrix_copy(cm_Fs(m)[SPN], cm_Fs(m)[SP]);
  vars[epn] = vars[ep];
  vars[gam_n] = vars[gam];
  vars[wn] = vars[w];
  vars[Xn] = vars[X];
  vars[Hn] = vars[H];
  flags[damaged_n] = flags[damaged];
  return err;
}

static int j2d_read_param(Model_parameters *p,
                          FILE *in)
{
  int err = 0;
  /* get pointer to parameter data */
  double *param = p->model_param;
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
  param[G] = p->p_hmat->G;
  param[nu] = p->p_hmat->nu;

  /* PLASTIC PROPERTIES */
  err += scan_for_valid_line(in);
  match += fscanf(in, "%lf %lf %lf",
                  param + beta, param + hp, param + k0);


  /* DAMAGE PROPERTIES */
  match += fscanf(in, "%lf %lf %lf %lf",
                  param + mu, param + p1,
                  param + p2, param + Yin);

  if (match != NUM_param) err++;
  assert(match == NUM_param && "Did not read expected number of parameters");

  /* scan past any other comment/blank lines in the block */
  err += scan_for_valid_line(in);

  /* not expecting EOF, check and return error if encountered */
  if (feof(in)) err ++;
  assert(!feof(in) && "EOF reached prematurely");

  return err;
}

static int j2d_plasticity_model_ctx_destroy(void **ctx)
{
  int err = 0;
  free(*ctx);
  *ctx = NULL;
  return err;
}

/* bbar = J^-(2/3) F F' */
static void j2d_bbar(const double * restrict F,
                     double * restrict bbar)
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
static void j2d_dev(const double * restrict a,
                    double * restrict dev_a)
{
  const double tra = 1./3. * a[0] + a[4] + a[8];
  memcpy(dev_a, a, tensor * sizeof(*a));
  dev_a[0] -= tra;
  dev_a[4] -= tra;
  dev_a[8] -= tra;
}

/* s0 = G * dev(bbar) */
static void j2d_compute_s0(const double G,
                           const double * restrict bbar,
                           double * restrict s0)
{
  double devbbar[tensor] = {0};
  j2d_dev(bbar, devbbar);
  for (int i = 0; i < tensor; i++) {
    s0[i] = G * devbbar[i];
  }
}

/* compute Fubar */
static int j2d_compute_Fubar(const double * restrict F,
                             const double * restrict Fn,
                             double * restrict Fubar)
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
static void j2d_push_forward(const double * restrict F,
                             const double * restrict A,
                             double * restrict a)
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

static double j2d_compute_normal(const double * restrict s_tr,
                                 const double * restrict sp_tr,
                                 const double * restrict param,
                                 double  * restrict n)
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
  for(int i = 0; i < tensor; i++) n[i] /= nrm;
  return nrm;
}

/* see box 4 of Simo and Ju (1989) */
static int j2d_int_alg(Constitutive_model *m,
                       const void *CTX)
{
  int err = 0;
  const j2d_ctx *ctx = CTX;

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
  double J23 = pow(det3x3(F), -2./3.);
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

  /* testing after the integration algorithm */
  {
    for (int i = 0; i < tensor; i++) s_tr[i] = s0[i] - sp[i];
    const double ksi_nrm2 = j2d_compute_normal(s_tr, sp, param, n);
    const double phi2 = ksi_nrm2 - sqrt(2./3.) * (param[k0] + param[beta] * param[hp] * vars[epn]);
    if (phi2 > j2d_int_alg_tol) {
      printf("UH OH... phi > 0 after integration algorithm...\n");
      abort();
    }
  }

  return err;
}

/* compute the deviatoric part of the effective (undamaged) Kirckhoff stress */
static int j2d_compute_sbar(const double *F,
                            const double *sp,
                            const double G,
                            double * restrict sbar)
{
  int err = 0;
  double bbar[tensor] = {0};
  j2d_bbar(F, bbar);
  j2d_compute_s0(G, bbar, sbar);
  for (int i = 0; i < tensor; i++) sbar[i] -= sp[i];
  return err;
}

/* compute the deviatoric PK2 stress */
static int j2d_Sdev(const Constitutive_model *m,
                    const void *CTX,
                    Matrix_double *Sdev)
{
  int err = 0;
  const j2d_ctx *ctx = CTX;
  double sbar[tensor] = {0};
  err += j2d_compute_sbar(ctx->F,
                          cm_Fs_data(m, SP),
                          cm_param(m)[G],
                          sbar);
  err += j2d_pull_back(ctx->F, sbar, Sdev->m_pdata);
  const double dam = 1. - cm_vars(m)[w];
  for(int i = 0; i < tensor; i++) Sdev->m_pdata[i] *= dam;
  return err;
}

/* compute the volumetric stress. NOTE we are still using the actual
   HOMMAT object here for convenience, but the only used information
   is for the flag */
static int j2d_dudj(const Constitutive_model *m,
                    const void *ctx,
                    double *dudj)
{
  int err = 0;
  const j2d_ctx *CTX = ctx;
  double J = det3x3(CTX->F);
  new_pot_compute_dudj(J, cm_hmat(m), dudj);

  /* scale by damage variable */
  *dudj *= (1.0 - cm_vars(m)[w]);
  return err;
}


int j2d_plasticity_model_initialize(Model_parameters *p)
{
  int err = 0;

  /* set function pointers */
  p->integration_algorithm = j2d_int_alg;
  p->compute_dev_stress = j2d_Sdev;
  p->compute_dudj = j2d_dudj;
  p->compute_dev_tangent = NULL;
  p->compute_d2udj2 = NULL;
  p->compute_AST = NULL;
  p->update_state_vars = j2d_update;
  p->reset_state_vars = j2d_reset;
  p->get_var_info = j2d_get_info;
  p->get_Fn = NULL;
  p->get_Fnm1 = NULL;
  p->get_pF = NULL;
  p->get_pFn = NULL;
  p->get_pFnm1 = NULL;
  p->get_eF = NULL;
  p->get_eFn = NULL;
  p->get_eFnm1 = NULL;

  p->get_hardening = NULL;
  p->get_hardening_nm1 = NULL;

  p->write_restart = NULL;
  p->read_restart = NULL;

  p->destroy_ctx = j2d_plasticity_model_ctx_destroy;
  p->compute_dMdu = NULL;

  p->set_init_vals = j2d_set_init_vals;
  p->read_param = j2d_read_param;

  p->get_size = NULL;
  p->pack = NULL;
  p->unpack = NULL;

  /* reset counters/flags */
  p->type = J2_PLASTICITY_DAMAGE;
  p->N_SYS = 0;

  /* allocate room for parameters */
  p->n_param = NUM_param;
  p->model_param = calloc(NUM_param, sizeof(*(p->model_param)));

  return err;
}

int j2d_plasticity_model_ctx_build(void **ctx,
                                   const double *F,
                                   const double dt)
{
  j2d_ctx *CTX = malloc(sizeof(*CTX));
  memcpy(CTX->F, F, tensor * sizeof(*F));
  CTX->dt = dt;
  *ctx = CTX;
  return 0;
}
