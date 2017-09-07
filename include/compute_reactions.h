/* HEADER */
/**
 * AUTHORS:
 * Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 */
#pragma once
#ifndef COMPUTE_REACTIONS_H
#define COMPUTE_REACTIONS_H

#include "data_structure.h"
#include "PGFEM_mpi.h"
#include "element.h"
#include "node.h"
#include "matgeom.h"
#include "hommat.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"

  /** After convergence, compute the reaction forces at nodes with
      prescribed deflections. This is achieved by computing the
      residuals on elements which contain node with prescribed
      deflections, filtering the result by deflection ID and summing
      together. NOTE: Only computes the reaction in the direction of
      the prescribed deflection!!! Since dV_u U dV_t = 0, there is no
      need to subtract external forces */

int compute_reactions(long ne,
              long ndofn,
              long npres,
              double *r,
              NODE *node,
              ELEMENT *elem,
              MATGEOM matgeom,
              HOMMAT *hommat,
              SUPP sup,
              EPS *eps,
              SIG *sig,
              double nor_min,
              CRPL *crpl,
              double dt,
              double stab,
              MPI_Comm mpi_comm,
              const int analysis,
              const int mp_id);

#endif /* #ifndef COMPUTE_REACTIONS_H */
