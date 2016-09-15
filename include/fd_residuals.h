/**
 * AUTHORS:
 * Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 * Karel Matous, University of Notre Dame, kmatous [at] nd.
 */
#pragma once
#ifndef FD_RESIDUALS_H
#define FD_RESIDUALS_H

#include "PGFEM_mpi.h"
#include "element.h"
#include "node.h"
#include "matgeom.h"
#include "hommat.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "cohesive_element.h"
#include "bounding_element.h"
#include "PGFem3D_options.h"

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
		    double *dts,
		    double t,
		    double stab,
		    long nce,
		    COEL *coel,
		    MPI_Comm mpi_comm,
		    const PGFem3D_opt *opts,
		    double alpha, double *r_n, double *r_n_1,
		    const int mp_id);

  /* compute the reaction force for each magnitude of prescribed
     deflection. CAVEATS: Does not include contributions from cohesive
     or boundary elements. */
  int fd_res_compute_reactions(const long ndofn,
                               const long npres,
                               const double *d_r,
                               const double *r,
                               ELEMENT *elem,
                               NODE *node,
                               MATGEOM matgeom,
                               HOMMAT *hommat,
                               SUPP sup,
                               EPS *eps,
                               SIG *sig,
                               const double nor_min,
                               CRPL *crpl,
                               const double *dts,
                               const double t,
                               const double stab,
                               MPI_Comm mpi_comm,
                               const PGFem3D_opt *opts,
                               const double alpha,
                               double *r_n,
                               double *r_n_1,
                               const int mp_id);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef FD_RESIDUALS_H */
