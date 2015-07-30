#ifndef _H_DYNAMICS_H_
#define _H_DYNAMICS_H_

#include "element.h"
#include "node.h"
#include "hommat.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "PGFem3D_options.h"

#define MIN_DENSITY 1.0e-16

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * Computes element stiffness matrices for transient terms
   */
void stiffmat_disp_w_inertia_el(double *Ks,
         const int ii,
         const int ndofn,
         const int nne, const int npres, const int nVol, const int nsd,
         const double *x, const double *y, const double *z,		     
         const ELEMENT *elem, const HOMMAT *hommat, const long *nod, const NODE *node, double dt,
         SIG *sig, EPS *eps, const SUPP sup, const int analysis,		     
		     double alpha, double *r_n, double *r_e);
		     
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
		     
void DISP_resid_w_inertia_el(double *f,
         const int ii,
         const int ndofn,
         const int nne,
         const double *x,
         const double *y,
         const double *z,		     
         const ELEMENT *elem,
         const HOMMAT *hommat,
		     const NODE *node, double dt, double t,
		     double *r_2, double* r_1, double *r_0, double alpha);	     
		     
int residuals_w_inertia_el(double *fe, int i, 
			int nne, long ndofn, long npres, long nVol,long ndofe, double *r_e,                               
		  NODE *node, ELEMENT *elem, HOMMAT *hommat, SUPP sup, EPS *eps, SIG *sig,
		  long* nod, long *cn, double *x, double *y, double *z,                                
		  double dt, double t, const PGFem3D_opt *opts, double alpha, double *r_n, double *r_n_1);
		  
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef _H_DYNAMICS_H_ */
