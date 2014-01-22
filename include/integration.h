#ifndef INTEGRATION_H
#define INTEGRATION_H

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifndef CRPL_H
#include "crpl.h"
#endif

#ifndef MATGEOM_H
#include "matgeom.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef SUPP_H
#include "supp.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Explicit hardening integration. */
  double int_har_explicit (long ii,
			   long ip,
			   long mat,
			   SIG *sig,
			   CRPL *crpl,
			   double dt);

  double int_alg_res (double lam,
		      double dt,
		      long mat,
		      long nss,
		      double *TA,
 		      double *GA,
		      CRPL *crpl,
		      double **UU,
		      double **Fp,
		      double **BB);

  long int_UU (long ii,
	       long ip,
	       double NOR_MIN,
	       double dt,
	       long mat,
	       double *NOR,
	       long nss,
	       CRPL *crpl,
	       double **Fr,
	       double **FnB,
	       double **Fp,
	       double **BB,
	       double *DD,
	       double **S,
	       double L[3][3][3][3],
	       double *TA,
	       double *GA,
 	       double Tr,
	       double Jr,
	       double Har,
	       double *LAM,
	       double **UU,
	       long GAMA,
	       long iter);

  long int_HA (long ii,
	       long ip,
	       long nss,
	       long mat,
	       double dt,
	       CRPL *crpl,
	       double *GA,
 	       double Har,
	       double *HAR,
	       EPS *eps,
	       const PGFem3D_opt *opts);

  long integration_alg (long ne,
			long ndofn,
			long ndofd,
			long npres,
			CRPL *crpl,
			ELEMENT *elem,
			NODE *node,
			double *d_r,
			double *rr,
			SUPP sup,
			MATGEOM matgeom,
			HOMMAT *hommat,
			EPS *eps,
			SIG *sig,
			long tim,
			long iter,
			double dt,
			double nor_min,
			long STEP,
			long GAMA,
			const PGFem3D_opt *opts);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef INTEGRATION_H */
