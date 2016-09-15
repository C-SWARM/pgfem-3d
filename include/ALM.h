#ifndef ALM_H
#define ALM_H

#include "PGFEM_mpi.h"

#ifndef MATGEOM_H
#include "matgeom.h"
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

#ifndef CRPL_H
#include "crpl.h"
#endif

#ifndef PGFEM_COMM_H
#include "pgfem_comm.h"
#endif

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef BOUNDING_ELEMENT_H
#include "bounding_element.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  double D_lam_ALM (long ndofd,
		    double *BS_rr,
		    double *BS_d_r,
		    double *BS_D_R,
		    double *BS_R,
		    double *BS_DK,
		    double dlm,
		    double dAL,
		    long *DomDof,
		    MPI_Comm mpi_comm);

  /** Returns 0. */
  double d_ALM2 (long ndofd,
		 double *rr,
		 double *R,
		 double *DK,
		 double d_lm);

  /** Returns 0. */
  double d_lam_ALM2 (long ndofd,
		     double *rr,
		     double *R,
		     double *DK,
		     double dAL,
		     double DET,
		     double DET0,
		     double dlm0,
		     double nor_min,
		     double *dR);

  double D_lam_ALM2 (double *BS_rr,
		     double *BS_D_R,
		     double *BS_R,
		     double *BS_DK,
		     double dlm,
		     double lm,
		     double dAL,
		     long ne,
		     int n_be,
		     long ndofd,
		     long npres,
		     double *BS_d_r,
		     double *r,
		     NODE *node,
		     ELEMENT *elem,
		     BOUNDING_ELEMENT *b_elems,
		     MATGEOM matgeom,
		     HOMMAT *hommat,
		     SUPP sup,
		     EPS *eps,
		     SIG *sig,
		     double nor_min,
		     CRPL *crpl,
		     double dt,
		     double stab,
		     long nce,
		     COEL *coel,
		     long *DomDof,
		     int GDof,
		     COMMUN comm 
		     /*,
		       GNOD *gnod,
		       GEEL *geel*/  ,
		     MPI_Comm mpi_comm,
		     const PGFem3D_opt *opts,
		     const int mp_id);

  double d_ALM4 (long ndofd,
		 double *BS_rr,
		 double *BS_DK,
		 double dlm,
		 long *DomDof,
		 MPI_Comm mpi_comm);

  double d_lam_ALM4 (long ndofd,
		     double *BS_rr,
		     double *BS_DK,
		     double *BS_dR,
		     double dAL,
		     long *DomDof,
		     MPI_Comm mpi_comm);

  double D_lam_ALM4 (long ndofd,
		     double *BS_rr,
		     double *BS_d_r,
		     double *BS_D_R,
		     double *BS_DK,
		     double dlm,
		     double dAL,
		     long *DomDof,
		     MPI_Comm mpi_comm);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef ALM_H */
