/* HEADER */

#pragma once
#ifndef FD_INCREMENT_H
#define FD_INCREMENT_H

#include "PGFEM_mpi.h"
#include "element.h"
#include "cohesive_element.h"
#include "matgeom.h"
#include "hommat.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "PGFem3D_options.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * Increment the 'fd' formulation elements after a completed step.
   */
  void fd_increment (long ne,
		     long nn,
		     long ndofn,
		     long npres,
		     MATGEOM matgeom,
		     HOMMAT *hommat,
		     ELEMENT *elem,
		     NODE *node,
		     SUPP sup,
		     EPS *eps,
		     SIG *sig,
		     double *d_r,
		     double *r,
		     double nor_min,
		     CRPL *crpl,
		     double dt,
		     long nce,
		     COEL *coel,
		     double *pores,
		     MPI_Comm mpi_comm,
		     const double VVolume,
		     const PGFem3D_opt *opts);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef FD_INCREMENT_H */
