#ifndef FD_RESIDUALS_H
#define FD_RESIDUALS_H

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
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

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifndef CRPL_H
#include "crpl.h"
#endif

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef BOUNDING_ELEMENT_H
#include "bounding_element.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  int fd_residuals (double *f_u,
		    long ne,
		    int n_be,
		    long ndofn,
		    long npres,
		    double *d_r,
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
		    MPI_Comm mpi_comm,
		    const PGFem3D_opt *opts);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef FD_RESIDUALS_H */
