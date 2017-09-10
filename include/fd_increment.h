#ifndef PGFEM3D_FD_INCREMENT_H
#define PGFEM3D_FD_INCREMENT_H

#include "PGFEM_mpi.h"
#include "PGFem3D_options.h"
#include "cohesive_element.h"
#include "crpl.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "matgeom.h"
#include "sig.h"

/**
 * Increment the 'fd' formulation elements after a completed step.
 */
void fd_increment (long ne,
                   long nn,
                   long ndofn,
                   long npres,
                   MATGEOM matgeom,
                   HOMMAT *hommat,
                   Element *elem,
                   Node *node,
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
                   const PGFem3D_opt *opts,
                   const int mp_id);

#endif /* #define PGFEM3D_FD_INCREMENT_H */
