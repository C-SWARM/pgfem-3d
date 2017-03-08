#ifndef _H_DYNAMICS_H_
#define _H_DYNAMICS_H_

#include "element.h"
#include "node.h"
#include "hommat.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "PGFem3D_options.h"
#include "PGFem3D_data_structure.h"

#define MIN_DENSITY 1.0e-16
#define DT_NP1 0
#define DT_N   1

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */
		     
void DISP_resid_body_force_el(double *f,
         const int ii,
         const int ndofn,
         const int nne,
         const double *x,
         const double *y,
         const double *z,		     
         const ELEMENT *elem,
         const HOMMAT *hommat,
		     const NODE *node, double dt, double t);		     
		     
struct FEMLIB;
#ifndef TYPE_FEMLIB
#define TYPE_FEMLIB
typedef struct FEMLIB FEMLIB;
#endif

int residual_with_inertia(FEMLIB *fe,
                          double *be,
                          double *r_e,
                          GRID *grid,
                          MATERIAL_PROPERTY *mat,
                          FIELD_VARIABLES *fv,
                          SOLVER_OPTIONS *sol,
                          LOADING_STEPS *load,
                          CRPL *crpl,
                          const PGFem3D_opt *opts,
                          MULTIPHYSICS *mp,
                          int mp_id,
                          double *dts,
                          double t);

int stiffness_with_inertia(FEMLIB *fe,
                           double *Ks,
                           double *r_e,
                           GRID *grid,
                           MATERIAL_PROPERTY *mat,
                           FIELD_VARIABLES *fv,
                           SOLVER_OPTIONS *sol,
                           LOADING_STEPS *load,
                           CRPL *crpl,
                           const PGFem3D_opt *opts,
                           MULTIPHYSICS *mp,
                           int mp_id,
                           double dt);
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef _H_DYNAMICS_H_ */
