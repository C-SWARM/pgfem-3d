#ifndef _H_DYNAMICS_H_
#define _H_DYNAMICS_H_

#include "element.h"
#include "node.h"
#include "hommat.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"


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

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef _H_DYNAMICS_H_ */
