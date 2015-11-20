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

// #include "cm_j2_plasticity.h"
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

Define_Matrix(double);

/* macros for easy access to Constitutive_model structure */
#define cm_Fs(m) ((m)->vars.Fs)
#define cm_vars(m) ((m)->vars.state_variables->m_pdata)
#define cm_func(m) ((m)->param)
#define cm_hmat(m) ((m)->param->p_hmat)

/* macro for accessing data of Fs */
#define Fs_data(Fs,idx) (Fs)[(idx)].m_pdata

/* private context structure */
typedef struct {
  double F[tensor];
  double dt;
} j2d_ctx;

enum {Fn, F, sp_n, sp, NUM_Fs};
enum {ep_n, ep, wn, w, Xn, X, Hn, H, NUM_vars};
enum {damaged_n, damaged, NUM_flags};
enum {G, nu, hp, beta, k0, mu, p1, p2, Yin, NUM_param};

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

/* DEV(A) = A - 1/3 (A:C) C^-1 */
static void j2d_DEV(const double * restrict A,
                    const double * restrict C,
                    double * restrict DEV_A)
{
  double CI[tensor] = {0};
  inv3x3(C,CI);
  double AC3 = 0;
  for (int i = 0; i < tensor; i++) {
    AC3 += A[i] * C[i] / 3.;
  }
  for (int i = 0; i < tensor; i++) {
    DEV_A[i] = A[i] - AC3 * CI[i];
  }
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

/* see box 4 of Simo and Ju (1989) */
static void j2d_radial_return(const double *F,
                              const double *Fn,
                              const double *param,
                              Constitutive_model *m)
{

}
