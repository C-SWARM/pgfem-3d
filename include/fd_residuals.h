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
