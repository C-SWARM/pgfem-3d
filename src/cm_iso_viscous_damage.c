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
#include "data_structure_c.h"

/* Define constant dimensions. Note cannot use `static const` with
   initialization list */
#define dim  3
#define tensor 9
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
enum {Yin, p1, p2, mu, NUM_param};

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

static int ivd_private_damage_int_alg(double *dw,
                                      double *vars,
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

  *dw = g;

  return err;
}

/* function stubs for CM interface */
static int ivd_int_alg(Constitutive_model *m,
                       const void *ctx)
{
  int err = 0;

  return err;
}

static int ivd_dev_stress(const Constitutive_model *m,
                          const void *ctx,
                          Matrix_double *devS)
{
  int err = 0;

  return err;
}

static int ivd_dudj(const Constitutive_model *m,
                    const void *ctx,
                    double *dudj)
{
  int err = 0;

  return err;
}

static int ivd_dev_tangent(const Constitutive_model *m,
                           const void *ctx,
                           Matrix_double *devL)
{
  int err = 0;

  return err;
}

static int ivd_d2udj2(const Constitutive_model *m,
                      const void *ctx,
                      double *d2udj2)
{
  int err = 0;

  return err;
}

static int ivd_update(Constitutive_model *m)
{
  int err = 0;
  double *vars = m->vars.state_vars->m_pdata;
  Matrix_copy(m->vars.Fs[Fn], m->vars.Fs[F]);
  vars[wn] = vars[w];
  vars[Xn] = vars[X];
  vars[Hn] = vars[H];

  /* flags? */
  return err;
}

static int ivd_reset(Constitutive_model *m)
{
  int err = 0;
  double *vars = m->vars.state_vars->m_pdata;
  Matrix_copy(m->vars.Fs[F], m->vars.Fs[Fn]);
  vars[w] = vars[wn];
  vars[X] = vars[Xn];
  vars[H] = vars[Hn];

  /* flags? */
  return err;
}

static int ivd_get_info(Model_var_info **info)
{
  int err = 0;

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

static int ivd_get_damage_nm1(const Constitutive_model *m,
                              double *w)
{
  *w = 0;
  return 0;
}

static int ivd_write_restart(FILE *out,
                             const Constitutive_model *m)
{
  int err = 0;
  const double *FF = m->vars.Fs[Fn].m_pdata;
  const double *vars = m->vars.state_vars->m_pdata;
  if(fprintf(out, "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
             FF[0], FF[1], FF[2], FF[3], FF[4], FF[5], FF[6], FF[7], FF[8]) < 0) err++;
  if(fprintf(out, "%.17e %.17e %.17e\n", vars[wn], vars[Xn], vars[Hn]) < 0) err++;

  /* do I need to write out the damaged flag(s)??? */

  return err;
}

static int ivd_read_restart(FILE *in,
                            Constitutive_model *m)
{
  int err = 0;
  double *FF = m->vars.Fs[Fn].m_pdata;
  double *vars = m->vars.state_vars->m_pdata;
  if(fscanf(in, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &FF[0], &FF[1], &FF[2], &FF[3],
            &FF[4], &FF[5], &FF[6], &FF[7], &FF[8]) != tensor) err++;
  if(fscanf(in, "%lf %lf %lf", &vars[wn], &vars[Xn], &vars[Hn]) != 2) err++;
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
  int err = 0;

  return err;
}

static int ivd_set_init_vals(Constitutive_model *m)
{
  int err = 0;
  Matrix_eye(m->vars.Fs[Fn], dim);
  double *vars = m->vars.state_vars->m_pdata;
  vars[wn] = vars[Xn] = vars[Hn] = 0.0;
  /* flags? */
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
  return 0;
}

static int ivd_pack(const Constitutive_model *m,
                    char *buffer,
                    size_t *pos)
{
  int err = 0;

  return err;
}

static int ivd_unpack(Constitutive_model *m,
                      const char *buffer,
                      size_t *pos)
{
  int err = 0;

  return err;
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
