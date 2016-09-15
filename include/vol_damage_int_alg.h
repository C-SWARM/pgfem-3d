/* HEADER */
#ifndef VOL_DAMAGE_INT_ALG_H
#define VOL_DAMAGE_INT_ALG_H

#include "PGFEM_mpi.h"

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef SUPP_H
#include "supp.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Compute new damage parameters */
  int vol_damage_int_alg(const int ne,
			 const int ndofn,
			 const double *d_r,
			 const double *r,
			 const ELEMENT *elem,
			 const NODE *node,
			 const HOMMAT *hommat,
			 const SUPP sup,
			 const double dt,
			 const int iter,
			 const MPI_Comm mpi_comm,
			 EPS *eps,
			 SIG *sig,
			 double *max_omega,
			 double *dissipation,
			 const int analysis,
			 const int mp_id);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef VOL_DAMAGE_INT_ALG_H */
