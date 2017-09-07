/* HEADER */
/**
 * This file defines the implementation for the isotropic viscous
 * damage model.
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
#include "cm_iso_viscous_damage.h"
#include "constitutive_model.h"
#include "index_macros.h"
#include "new_potentials.h"
#include "utils.h"
#include <math.h>
#include <string.h>
#include <assert.h>

/* Define constant dimensions. Note cannot use `static const` with
   initialization list */
namespace {
constexpr int              dim = 3;
constexpr int           tensor = 9;
constexpr int          tensor4 = 81;
constexpr double DAMAGE_THRESH = 0.9999;
constexpr double   DELTA_W_MAX = 0.05;

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* private context structure */
struct ivd_ctx {
  double F[tensor];
  double dt;
};

enum {Fn, F, NUM_Fs};
enum {wn, w, Xn, X, Hn, H, NUM_vars};
enum {damaged_n, damaged, NUM_flags};
enum {mu, ome_max, p1, p2, Yin, NUM_param};
}

/**
 * Matrix multiplication b = a'a, dim(a) = [3 3]
 */
static void ata(const double *a, double *b) {
  memset(b, 0, tensor * sizeof(*b));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for(int k = 0; k < dim; k++) {
        b[idx_2(i,j)] += a[idx_2(k,i)] * a[idx_2(k,j)];
      }
    }
  }
}

static double ivd_weibull_function(const double Ybar, const double *params) {
  if(Ybar <= params[Yin]) return 0.0;
  return (params[ome_max] - params[ome_max]
      * exp(- pow((Ybar - params[Yin]) / (params[p1] * params[Yin]),
            params[p2])
           )
      );
}

static double ivd_weibull_evolution(const double Ybar, const double *params) {
  if(Ybar <= params[Yin]) return 0.0;
  return (params[ome_max] * params[p2] / (params[p1] * params[Yin])
      * exp(- pow((Ybar - params[Yin]) / (params[p1] * params[Yin]), params[p2]) )
      * pow((Ybar - params[Yin]) / (params[p1] * params[Yin]), params[p2] - 1.0)
      );
}

static int ivd_private_damage_int_alg(double *vars,
                                      int *flags,
                                      const double *params,
                                      const double Ybar,
                                      const double dt)
{
  int err = 0;
  const double dmu = dt * params[mu];
  double G = ivd_weibull_function(Ybar, params);
  double g = G - vars[Xn];

  if (g > 0.0){ /* Damage propagation */
    /* flag integration point */
    flags[damaged] = 1;

    /* update damage parameter */
    vars[w] = vars[wn] + dmu / (1 + dmu) * g;

    /* update softening parameter */
    vars[X] = MAX(vars[Xn],(vars[Xn] + dmu * G) / (1 + dmu));

    /* update evolution parameter */
    vars[H] = ivd_weibull_evolution(Ybar, params);

  } else { /* no damage propagation */
    flags[damaged] = 0;
    vars[w] = vars[wn];
    vars[X] = vars[Xn];
    vars[H] = 0.0;
  }

  return err;
}

int ivd_public_int_alg(double *var_w,
                       double *var_X,
                       double *var_H,
                       int *flag_damaged,
                       const double var_wn,
                       const double var_Xn,
                       const double dt,
                       const double Ybar,
                       const double param_mu,
                       const double param_ome_max,
                       const double param_p1,
                       const double param_p2,
                       const double param_Yin)
{
  int err = 0;
  double *params = PGFEM_calloc(double, NUM_param);
  double *vars = PGFEM_calloc(double, NUM_vars);
  int *flags = PGFEM_calloc(int, NUM_flags);

  /* pack state at n */
  params[mu] = param_mu;
  params[ome_max] = param_ome_max;
  params[p1] = param_p1;
  params[p2] = param_p2;
  params[Yin] = param_Yin;
  vars[wn] = var_wn;
  vars[Xn] = var_Xn;

  /* run the integration algorithm */
  err +=  ivd_private_damage_int_alg(vars, flags, params, Ybar, dt);

  /* unpack the state at n+1 */
  *var_w = vars[w];
  *var_X = vars[X];
  *var_H = vars[H];
  *flag_damaged = flags[damaged];

  /* cleanup and exit */
  free(params);
  free(vars);
  free(flags);
  return err;
}

static int ivd_compute_Sbar(const Constitutive_model *m,
                            const void *ctx,
                            double *Sbar)
{
  int err = 0;
  const double kappa = hommat_get_kappa(m->param->p_hmat);
  auto CTX = (ivd_ctx *) ctx;
  const double J = det3x3(CTX->F);
  double C[tensor] = {0};
  double CI[tensor] = {0};
  ata(CTX->F,C);
  err += inv3x3(C, CI);
  new_pot_compute_Sdev(C, m->param->p_hmat,Sbar);

  double dudj = 0;
  new_pot_compute_dudj(J, m->param->p_hmat, &dudj);

  for (int i = 0; i < tensor; i++) {
    Sbar[i] += kappa * J * dudj * CI[i];
  }

  return err;
}

static int ivd_compute_Lbar(const Constitutive_model *m,
                            const void *ctx,
                            double *Lbar)
{
  int err = 0;
  const double kappa = hommat_get_kappa(m->param->p_hmat);
  auto CTX = (ivd_ctx *) ctx;
  const double J = det3x3(CTX->F);
  double C[tensor] = {0};
  double C_I[tensor] = {0};
  double dudj = 0.0;
  double d2udj2 = 0.0;

  ata(CTX->F,C);
  err += inv3x3(C, C_I);
  new_pot_compute_Ldev(C, m->param->p_hmat, Lbar);
  new_pot_compute_dudj(J, m->param->p_hmat, &dudj);
  new_pot_compute_d2udj2(J, m->param->p_hmat, &d2udj2);

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

static int ivd_modify_AST(const Constitutive_model *m,
                          const void *ctx,
                          double *L)
{
  int err = 0;
  auto CTX = (ivd_ctx *) ctx;
  const double dmu = CTX->dt * m->param->model_param[mu];
  const double evo = dmu * m->vars_list[0][m->model_id].state_vars->m_pdata[H] / (1.0 + dmu);
  double Sbar[tensor] = {0};
  err += ivd_compute_Sbar(m, ctx, Sbar);

  for (int i=0; i < dim; i++) {
    for (int j=0; j < dim; j++) {
      for (int k=0; k < dim; k++) {
    for (int l=0; l < dim; l++) {
      const int idx4 = idx_4(i,j,k,l);
      /* Deviatoric + Volumetric stiffness */
      L[idx4] -= evo * Sbar[idx_2(i,j)] * Sbar[idx_2(k,l)];
        }
      }
    }
  }

  return err;
}
/* function stubs for CM interface */
static int ivd_compute_AST(const Constitutive_model *m,
                           const void *ctx,
                           double *L)
{
  int err = 0;
  /* compute virgin material tangent */
  err += ivd_compute_Lbar(m, ctx, L);

  /* scale by damage parameter */
  const double dam = 1.0 - m->vars_list[0][m->model_id].state_vars->m_pdata[w];
  for (int i = 0; i < tensor4; i++) {
    L[i] *= dam;
  }

  /* if evolving damage, modify tangent */
  if (m->vars_list[0][m->model_id].flags[damaged]) {
    err += ivd_modify_AST(m, ctx, L);
  }

  return err;
}

int CM_IVD_PARAM::integration_algorithm(Constitutive_model *m,
                                        const void *ctx)
const
{
  int err = 0;
  auto CTX = (ivd_ctx *) ctx;
  const double *params = m->param->model_param;
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  int *flags = m->vars_list[0][m->model_id].flags;

  double Wdev = 0.0;
  double U = 0.0;
  double C[tensor] = {0};
  ata(CTX->F,C);
  const double J = det3x3(CTX->F);
  new_pot_compute_Wdev(C, m->param->p_hmat, &Wdev);
  new_pot_compute_U(J, m->param->p_hmat, &U);
  const double Ybar =  Wdev + U * hommat_get_kappa(m->param->p_hmat);

  /* integration algorithm */
  err += ivd_private_damage_int_alg(vars, flags, params, Ybar, CTX->dt);

  /* store the deformation gradient */
  memcpy(m->vars_list[0][m->model_id].Fs[F].m_pdata, CTX->F, tensor * sizeof(double));

  return err;
}

int CM_IVD_PARAM::compute_dev_stress(const Constitutive_model *m,
                                     const void *ctx,
                                     double *stress)
const
{
  int err = 0;
  auto CTX = (ivd_ctx *) ctx;
  double C[tensor] = {0};
  ata(CTX->F,C);
  new_pot_compute_Sdev(C, m->param->p_hmat, stress);

  /* scale by damage variable */
  const double dam = 1.0 - m->vars_list[0][m->model_id].state_vars->m_pdata[w];
  for(int i = 0; i < tensor; i++) stress[i] *= dam;
  return err;
}

int CM_IVD_PARAM::compute_dudj(const Constitutive_model *m,
                               const void *ctx,
                               double *dudj)
const
{
  int err = 0;
  auto CTX = (ivd_ctx *) ctx;
  double J = det3x3(CTX->F);
  new_pot_compute_dudj(J, m->param->p_hmat, dudj);

  /* scale by damage variable */
  *dudj *= (1.0 - m->vars_list[0][m->model_id].state_vars->m_pdata[w]);
  return err;
}


int CM_IVD_PARAM::compute_dev_tangent(const Constitutive_model *m,
                                      const void *ctx,
                                      double *L)
const
{
  int err = 0;
  auto CTX = (ivd_ctx *) ctx;
  double C[tensor] = {0};
  ata(CTX->F,C);
  new_pot_compute_Ldev(C, m->param->p_hmat, L);

  /* scale by damage variable */
  const double dam = 1.0 - m->vars_list[0][m->model_id].state_vars->m_pdata[w];
  for(int i = 0; i < tensor4; i++) L[i] *= dam;

  return err;
}

int CM_IVD_PARAM::compute_d2udj2(const Constitutive_model *m,
                                 const void *ctx,
                                 double *d2udj2)
const
{
  int err = 0;
  auto CTX = (ivd_ctx *) ctx;
  double J = det3x3(CTX->F);
  new_pot_compute_d2udj2(J, m->param->p_hmat, d2udj2);

  /* scale by damage variable */
  *d2udj2 *= (1.0 - m->vars_list[0][m->model_id].state_vars->m_pdata[w]);
  return err;
}

int CM_IVD_PARAM::update_state_vars(Constitutive_model *m)
const
{
  int err = 0;
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  int *flags = m->vars_list[0][m->model_id].flags;

  /* deformation gradients */
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  memcpy(Fs[F].m_pdata, Fs[Fn].m_pdata, sizeof(double)*tensor);

  /* state variables */
  vars[wn] = vars[w];
  vars[Xn] = vars[X];
  vars[Hn] = vars[H];

  /* flags */
  flags[damaged_n] = flags[damaged];

  return err;
}

int CM_IVD_PARAM::reset_state_vars(Constitutive_model *m)
const 
{
  int err = 0;
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  int *flags = m->vars_list[0][m->model_id].flags;

  /* deformation gradients */
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  memcpy(Fs[Fn].m_pdata, Fs[F].m_pdata, sizeof(double)*tensor);  

  /* state variables */
  vars[w] = vars[wn];
  vars[X] = vars[Xn];
  vars[H] = vars[Hn];

  /* flags */
  flags[damaged] = flags[damaged_n];

  return err;
}

int CM_IVD_PARAM::get_var_info(Model_var_info &info)
const
{
  int err = 0;

  info.n_Fs = NUM_Fs;
  info.n_vars = NUM_vars;
  info.n_flags    = NUM_flags;
  info.F_names    = PGFEM_malloc<char*>(NUM_Fs);
  info.var_names  = PGFEM_malloc<char*>(NUM_vars);
  info.flag_names = PGFEM_malloc<char*>(NUM_flags);

  /* allocate/copy strings */
  info.F_names[F] = strdup("F");
  info.F_names[Fn]   = strdup("Fn");
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

int ivd_get_F(const Constitutive_model *m,
              double *F_out,
              const int stepno)
{
  int err = 0;
  Matrix<double> *Fs = m->vars_list[0][m->model_id].Fs;
  switch(stepno)
  {
    case 0: // n-1
      memcpy(F_out, Fs[Fn].m_pdata, tensor*sizeof(double));
      break;
    case 1: // n
      memcpy(F_out, Fs[Fn].m_pdata, tensor*sizeof(double));
      break;
    case 2: // n+1
      memcpy(F_out, Fs[F].m_pdata, tensor*sizeof(double));
      break;
    default:
      PGFEM_printerr("ERROR: Unrecognized step number (%zd)\n",stepno);
      err++;
  }
  assert(err == 0);
  return err;    
}


int CM_IVD_PARAM::get_F(const Constitutive_model *m,
                        double *F_out,
                        const int stepno)
const
{
  return ivd_get_F(m, F_out, stepno);
}


int CM_IVD_PARAM::get_eF(const Constitutive_model *m,
                         double *eF_in,
                         const int stepno)
const
{
  return ivd_get_F(m, eF_in, stepno);
}

/// Return the identity tensor for the deformation gradient. Used for ,
/// e.g., Fp
///
int CM_IVD_PARAM::get_pF(const Constitutive_model *m,
                         double *F_out,
                         const int stepno)
const
{
  double I[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
  memcpy(F_out,I,tensor*sizeof(double));
  return 0;
}

int CM_IVD_PARAM::get_hardening(const Constitutive_model *m,
                                double *var,
                                const int stepno)
const
{
  int err = 0;
  double *s_var = m->vars_list[0][m->model_id].state_vars->m_pdata;
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

int CM_IVD_PARAM::get_plast_strain_var(const Constitutive_model *m,
                                       double *chi)
const                                   
{
  *chi = m->vars_list[0][m->model_id].state_vars->m_pdata[Xn];
  return 0;
}

int CM_IVD_PARAM::write_restart(FILE *out, const Constitutive_model *m)
const
{
  int err = 0;
  const double *FF = m->vars_list[0][m->model_id].Fs[Fn].m_pdata;
  const double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  const int *flags = m->vars_list[0][m->model_id].flags;
  if(fprintf(out, "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
             FF[0], FF[1], FF[2], FF[3], FF[4], FF[5], FF[6], FF[7], FF[8]) < 0) err++;
  if(fprintf(out, "%.17e %.17e %.17e %d\n", vars[wn], vars[Xn], vars[Hn], flags[damaged_n]) < 0) err++;
  return err;
}

int CM_IVD_PARAM::read_restart(FILE *in, Constitutive_model *m)
const
{
  int err = 0;
  double *FF = m->vars_list[0][m->model_id].Fs[Fn].m_pdata;
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  int *flags = m->vars_list[0][m->model_id].flags;
  if(fscanf(in, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &FF[0], &FF[1], &FF[2], &FF[3],
            &FF[4], &FF[5], &FF[6], &FF[7], &FF[8]) != tensor) err++;
  if(fscanf(in, "%lf %lf %lf %d", &vars[wn], &vars[Xn], &vars[Hn], &flags[damaged_n]) != 4) err++;
  err += this->reset_state_vars(m);
  return err;
}

int CM_IVD_PARAM::set_init_vals(Constitutive_model *m)
const
{
  int err = 0;
  double I[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
  memcpy(m->vars_list[0][m->model_id].Fs[Fn].m_pdata,I, sizeof(double)*tensor);
  double *vars = m->vars_list[0][m->model_id].state_vars->m_pdata;
  vars[wn] = vars[Xn] = vars[Hn] = 0.0;
  m->vars_list[0][m->model_id].flags[damaged_n] = 0;
  err += this->reset_state_vars(m);
  return err;
}

int CM_IVD_PARAM::read_param(FILE *in)
const
{
  int err = 0;
  /* get pointer to parameter data */
  double *param     = this->model_param;
  assert(param != NULL); // check the pointer

  /* scan to non-blank/comment line */
  err += scan_for_valid_line(in);

  /* READ PROPERTIES IN ALPHABETICAL ORDER */
  int match = fscanf(in, "%lf %lf %lf %lf %lf",
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

int CM_IVD_PARAM::get_subdiv_param(const Constitutive_model *m,
                                   double *subdiv_param,
                                   double dt)
const
{
  return ivd_public_subdiv_param(m->vars_list[0][m->model_id].state_vars->m_pdata[wn],
                                 m->vars_list[0][m->model_id].state_vars->m_pdata[w],
                                 subdiv_param);
}

int CM_IVD_PARAM::update_elasticity(const Constitutive_model *m,
                                    const void *ctx_in,
                                    double *L,
                                    double *S,
                                    const int compute_stiffness)
const
{
  int err = 0;

  err += ivd_compute_Sbar(m,ctx_in,S);
  double dam = (1.0 - m->vars_list[0][m->model_id].state_vars->m_pdata[w]);
  for(int a=0; a<tensor; a++)
    S[a] *= dam;

  if(compute_stiffness)
    err += ivd_compute_AST(m, ctx_in, L); //compute stiffness

  return err;
}

/* API functions */

int CM_IVD_PARAM::model_dependent_initialization(void)
{
  int err = 0;
  this->type = ISO_VISCOUS_DAMAGE;
  this->n_param = NUM_param;
  this->model_param = PGFEM_calloc(double, NUM_param);

  return err;
}

int iso_viscous_damage_model_ctx_build(void **ctx,
                                       const double *F,
                                       const double dt)
{
  ivd_ctx *CTX = PGFEM_malloc<ivd_ctx>();
  memcpy(CTX->F, F, tensor * sizeof(*F));
  CTX->dt = dt;
  *ctx = CTX;
  return 0;
}

int CM_IVD_PARAM::destroy_ctx(void **ctx)
const
{
  int err = 0;
  free(*ctx);
  *ctx = NULL;
  return err;
}

int ivd_public_subdiv_param(const double var_wn,
                            const double var_w,
                            double *subdiv_param)
{
  *subdiv_param = (var_w - var_wn) / DELTA_W_MAX;
  return 0;
}
