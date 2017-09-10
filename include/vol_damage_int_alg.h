/* HEADER */
#ifndef PGFEM3D_VOL_DAMAGE_INT_ALG_H
#define PGFEM3D_VOL_DAMAGE_INT_ALG_H

#include "PGFEM_mpi.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "node.h"
#include "sig.h"
#include "supp.h"

/** Compute new damage parameters */
int vol_damage_int_alg(const int ne,
                       const int ndofn,
                       const double *d_r,
                       const double *r,
                       const Element *elem,
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

#endif /* #define PGFEM3D_VOL_DAMAGE_INT_ALG_H */
