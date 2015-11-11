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

#include "cm_iso_viscous_damage.h"
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
#define DAMAGE_THRESH 0.9999
#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

Define_Matrix(double);

/* private context structure */
typedef struct {
  double F[tensor];
  double dt;
} ivd_ctx;

enum {Fn, F, NUM_Fs};
enum {wn, w, Xn, X, Hn, H, NUM_vars};
enum {damaged_n, damaged, NUM_flags};
enum {mu, p1, p2, Yin, NUM_param};

/**
 * Matrix multiplication b = a'a, dim(a) = [3 3]
 */
static void ata(const double *a,
                double *b)
{
  memset(b, 0, tensor * sizeof(*b));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for(int k = 0; k < dim; k++) {
        b[idx_2(i,j)] += a[idx_2(k,i)] * a[idx_2(k,j)];
      }
    }
  }
}

static double ivd_weibull_function(const double Ybar,
                                   const double *params)
{
  if(Ybar <= params[Yin]) return 0.0;
  return (DAMAGE_THRESH - DAMAGE_THRESH
	  * exp(- pow((Ybar - params[Yin]) / (params[p1] * params[Yin]),
		    params[p2])
	       )
	  );
}

static double ivd_weibull_evolution(const double Ybar,
                                    const double *params)
{
  if(Ybar <= params[Yin]) return 0.0;
  return (DAMAGE_THRESH * params[p2] / (params[p1] * params[Yin])
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

static int ivd_compute_Sbar(const Constitutive_model *m,
                            const void *ctx,
                            double *Sbar)
{
  int err = 0;
  const double kappa = hommat_get_kappa(m->param->p_hmat);
  const ivd_ctx *CTX = ctx;
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
  const ivd_ctx *CTX = ctx;
  const double J = det3x3(CTX->F);
  double C[tensor] = {0};
  double C_I[tensor] = {0};
  double dudj = 0.0;
  double d2udj2 = 0.0;

  ata(CTX->F,C);
  err += inv3x3(C, C_I);
  new_pot_compute_Ldev(C, m->param->p_hmat, Lbar);
  new_pot_compute_dudj(J, m->param->p_hmat, &dudj);
  new_pot_compute_d2udj2(J, m->param->p_hmat, &dudj);

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
  const ivd_ctx *CTX = ctx;
  const double dmu = CTX->dt * m->param->model_param[mu];
  const double evo = dmu * m->vars.state_vars->m_pdata[H] / (1.0 + dmu);
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
                           Matrix_double *L)
{
  int err = 0;
  err += ivd_compute_Lbar(m, ctx, L->m_pdata);

  if (m->vars.flags[damaged]) {
    err += ivd_modify_AST(m, ctx, L->m_pdata);
  }

  return err;
}

static int ivd_int_alg(Constitutive_model *m,
                       const void *ctx)
{
  int err = 0;
  const ivd_ctx *CTX = ctx;
  const double *params = m->param->model_param;
  double *vars = m->vars.state_vars->m_pdata;
  int *flags = m->vars.flags;

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
  memcpy(m->vars.Fs[F].m_pdata, CTX->F, tensor * sizeof(double));

  return err;
}

static int ivd_dev_stress(const Constitutive_model *m,
                          const void *ctx,
                          Matrix_double *devS)
{
  int err = 0;
  const ivd_ctx *CTX = ctx;
  double C[tensor] = {0};
  ata(CTX->F,C);
  new_pot_compute_Sdev(C, m->param->p_hmat, devS->m_pdata);

  /* scale by damage variable */
  const double dam = 1.0 - m->vars.state_vars->m_pdata[w];
  for(int i = 0; i < tensor; i++) devS->m_pdata[i] *= dam;
  return err;
}

static int ivd_dudj(const Constitutive_model *m,
                    const void *ctx,
                    double *dudj)
{
  int err = 0;
  const ivd_ctx *CTX = ctx;
  double J = det3x3(CTX->F);
  new_pot_compute_dudj(J, m->param->p_hmat, dudj);

  /* scale by damage variable */
  *dudj *= (1.0 - m->vars.state_vars->m_pdata[w]);
  return err;
}

static int ivd_dev_tangent(const Constitutive_model *m,
                           const void *ctx,
                           Matrix_double *devL)
{
  int err = 0;
  const ivd_ctx *CTX = ctx;
  double C[tensor] = {0};
  ata(CTX->F,C);
  new_pot_compute_Ldev(C, m->param->p_hmat, devL->m_pdata);

  /* scale by damage variable */
  const double dam = 1.0 - m->vars.state_vars->m_pdata[w];
  for(int i = 0; i < tensor4; i++) devL->m_pdata[i] *= dam;

  return err;
}

static int ivd_d2udj2(const Constitutive_model *m,
                      const void *ctx,
                      double *d2udj2)
{
  int err = 0;
  const ivd_ctx *CTX = ctx;
  double J = det3x3(CTX->F);
  new_pot_compute_d2udj2(J, m->param->p_hmat, d2udj2);

  /* scale by damage variable */
  *d2udj2 *= (1.0 - m->vars.state_vars->m_pdata[w]);
  return err;
}

static int ivd_update(Constitutive_model *m)
{
  int err = 0;
  double *vars = m->vars.state_vars->m_pdata;
  int *flags = m->vars.flags;

  /* deformation gradients */
  Matrix_copy(m->vars.Fs[Fn], m->vars.Fs[F]);

  /* state variables */
  vars[wn] = vars[w];
  vars[Xn] = vars[X];
  vars[Hn] = vars[H];

  /* flags */
  flags[damaged_n] = flags[damaged];

  return err;
}

static int ivd_reset(Constitutive_model *m)
{
  int err = 0;
  double *vars = m->vars.state_vars->m_pdata;
  int *flags = m->vars.flags;

  /* deformation gradients */
  Matrix_copy(m->vars.Fs[F], m->vars.Fs[Fn]);

  /* state variables */
  vars[w] = vars[wn];
  vars[X] = vars[Xn];
  vars[H] = vars[Hn];

  /* flags */
  flags[damaged] = flags[damaged_n];

  return err;
}

static int ivd_get_info(Model_var_info **info)
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
  (*info)->F_names[F] = strdup("F");
  (*info)->F_names[Fn]   = strdup("Fn");

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

static int ivd_get_Fn(const Constitutive_model *m,
                      Matrix_double *FF)
{
  int err = 0;
  Matrix_copy(*FF, m->vars.Fs[Fn]);
  return err;
}

static int ivd_get_F(const Constitutive_model *m,
                     Matrix_double *FF)
{
  int err = 0;
  Matrix_copy(*FF, m->vars.Fs[F]);
  return err;
}

/**
 * Return the identity tensor for the deformation gradient. Used for ,
 * e.g., Fp
 */
static int ivd_get_F_I(const Constitutive_model *m,
                       Matrix_double *F)
{
  int err = 0;
  Matrix_eye(*F, dim);
  return err;
}

static int ivd_get_damage(const Constitutive_model *m,
                          double *w)
{
  *w = m->vars.state_vars->m_pdata[wn];
  return 0;
}

static int ivd_write_restart(FILE *out,
                             const Constitutive_model *m)
{
  int err = 0;
  const double *FF = m->vars.Fs[Fn].m_pdata;
  const double *vars = m->vars.state_vars->m_pdata;
  const int *flags = m->vars.flags;
  if(fprintf(out, "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
             FF[0], FF[1], FF[2], FF[3], FF[4], FF[5], FF[6], FF[7], FF[8]) < 0) err++;
  if(fprintf(out, "%.17e %.17e %.17e %d\n", vars[wn], vars[Xn], vars[Hn], flags[damaged_n]) < 0) err++;

  /* do I need to write out the damaged flag(s)??? */

  return err;
}

static int ivd_read_restart(FILE *in,
                            Constitutive_model *m)
{
  int err = 0;
  double *FF = m->vars.Fs[Fn].m_pdata;
  double *vars = m->vars.state_vars->m_pdata;
  int *flags = m->vars.flags;
  if(fscanf(in, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &FF[0], &FF[1], &FF[2], &FF[3],
            &FF[4], &FF[5], &FF[6], &FF[7], &FF[8]) != tensor) err++;
  if(fscanf(in, "%lf %lf %lf %d", &vars[wn], &vars[Xn], &vars[Hn], &flags[damaged_n]) != 4) err++;
  err += ivd_reset(m);
  return err;
}

static int ivd_compute_dMdu(const Constitutive_model *m,
                            const void *ctx,
                            const double *Grad_op,
                            const int nne,
                            const int ndofn,
                            double *dM_du)
{
  memset(dM_du, 0, tensor * nne *dim * sizeof(*dM_du));
  return 0;
}

static int ivd_set_init_vals(Constitutive_model *m)
{
  int err = 0;
  Matrix_eye(m->vars.Fs[Fn], dim);
  double *vars = m->vars.state_vars->m_pdata;
  vars[wn] = vars[Xn] = vars[Hn] = 0.0;
  m->vars.flags[damaged_n] = 0;
  err += ivd_reset(m);
  return err;
}

static int ivd_read_param(Model_parameters *p,
                          FILE *in)
{
  int err = 0;
  /* get pointer to parameter data */
  double *param = p->model_param;
  assert(param != NULL); // check the pointer

  /* scan to non-blank/comment line */
  err += scan_for_valid_line(in);

  /* READ PROPERTIES IN ALPHABETICAL ORDER */
  int match = fscanf(in, "%lf %lf %lf %lf",
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

static size_t ivd_get_size(const Constitutive_model *m)
{
  return state_variables_get_packed_size(&(m->vars));
}

static int ivd_pack(const Constitutive_model *m,
                    char *buffer,
                    size_t *pos)
{
  return state_variables_pack(&(m->vars), buffer, pos);
}

static int ivd_unpack(Constitutive_model *m,
                      const char *buffer,
                      size_t *pos)
{
  return state_variables_unpack(&(m->vars), buffer, pos);
}

/* API functions */
int iso_viscous_damage_model_initialize(Model_parameters *p)
{
  int err = 0;

  /* set function pointers */
  p->integration_algorithm = ivd_int_alg;
  p->compute_dev_stress = ivd_dev_stress;
  p->compute_dudj = ivd_dudj;
  p->compute_dev_tangent = ivd_dev_tangent;
  p->compute_d2udj2 = ivd_d2udj2;
  p->compute_AST = ivd_compute_AST;
  /* p->compute_AST = NULL; */
  p->update_state_vars = ivd_update;
  p->reset_state_vars = ivd_reset;
  p->get_var_info = ivd_get_info;
  p->get_Fn = ivd_get_Fn;
  p->get_Fnm1 = NULL;
  p->get_pF = ivd_get_F_I;
  p->get_pFn = ivd_get_F_I;
  p->get_pFnm1 = NULL;
  p->get_eF = ivd_get_F;
  p->get_eFn = ivd_get_Fn;
  p->get_eFnm1 = NULL;

  p->get_hardening = ivd_get_damage;
  p->get_hardening_nm1 = NULL;

  p->write_restart = ivd_write_restart;
  p->read_restart = ivd_read_restart;

  p->destroy_ctx = iso_viscous_damage_model_ctx_destroy;
  p->compute_dMdu = ivd_compute_dMdu;

  p->set_init_vals = ivd_set_init_vals;
  p->read_param = ivd_read_param;

  p->get_size = ivd_get_size;
  p->pack = ivd_pack;
  p->unpack = ivd_unpack;

  /* reset counters/flags */
  p->type = ISO_VISCOUS_DAMAGE;
  p->N_SYS = 0;

  /* allocate room for parameters */
  p->n_param = NUM_param;
  p->model_param = calloc(NUM_param, sizeof(*(p->model_param)));

  return err;
}

int iso_viscous_damage_model_ctx_build(void **ctx,
                                       const double *F,
                                       const double dt)
{
  ivd_ctx *CTX = malloc(sizeof(*CTX));
  memcpy(CTX->F, F, tensor * sizeof(*F));
  CTX->dt = dt;
  *ctx = CTX;
  return 0;
}

int iso_viscous_damage_model_ctx_destroy(void **ctx)
{
  int err = 0;
  free(*ctx);
  *ctx = NULL;
  return err;
}
