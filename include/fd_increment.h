#ifndef FD_INCREMENT_H
#define FD_INCREMENT_H

#include "PGFEM_mpi.h"

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef MATGEOM_H
#include "matgeom.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
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

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
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
