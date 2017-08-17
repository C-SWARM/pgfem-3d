/* HEADER */
#pragma once
#ifndef MINI_ELEMENT_H
#define MINI_ELEMENT_H

#include "PGFEM_mpi.h"
#include "element.h"
#include "node.h"
#include "hommat.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "pgfem_comm.h"

  /** Reset the internal variables on the elements */
  void MINI_reset(ELEMENT *elem,
          const int nelem,
          const int npres,
          SIG *sig);

  /** Compute the element stiffness matrix for the linear bubble
      element */
  int MINI_stiffmat_el(double *Ks,            /**< Element stiffmat */
               const int ii,          /**< element id */
               const int ndofn,
               const int nne,
               const double *x,
               const double *y,
               const double *z,
               const ELEMENT *elem,
               const HOMMAT *hommat,
               const long *nod,
               const NODE *node,
               const EPS *eps,
               const SIG *sig,
               const double *r_e);    /**< dof values on elem */

  /** Compute the element residual vector for the linear bubble
      element */
  int MINI_resid_el(double *Res,         /**< Element residual */
            const int ii,        /**< element id */
            const int ndofn,
            const int nne,
            const double *x,
            const double *y,
            const double *z,
            const ELEMENT *elem,
            const long *nod,
            const NODE *node,
            const HOMMAT *hommat,
            const EPS *eps,
            const SIG *sig,
            const double *r_e);    /**< dof values on elem */

  /** Update the bubble dofs on an element */
  int MINI_update_bubble_el(ELEMENT *elem,
                const int ii, /* id of element working on */
                const int nne,
                const NODE *node,
                const int ndofn,
                const double *x,
                const double *y,
                const double *z,
                const long *nod, /* list of node ids on elem */
                const EPS *eps,
                const SIG *sig,
                const HOMMAT *hommat,
                const double *sol_e, /* accum. incr on elem */
                const double *dsol_e); /* increment on element */

  /** Update the bubble dofs on ALL elements */
  int MINI_update_bubble(ELEMENT *elem,
             const int nelem,
             const NODE *node,
             const int ndofn,
             SUPP sup,
             const EPS *eps,
             const SIG *sig,
             const HOMMAT *hommat,
             const double *sol, /* accum. solution  on incr */
             const double *dsol, /* sol from current iter */
             const int iter,
             const int mp_id);

  /** Increment an element after a converged solution is obtained */
  void MINI_increment_el(ELEMENT *elem,
             const int ii, /* id of element working on */
             const int nne,
             const NODE *node,
             const long *nod,
             const int ndofn,
             const double *x,
             const double *y,
             const double *z,
             EPS *eps,
             SIG *sig,
             const HOMMAT *hommat,
             const double *sol_e);

  /** Increment ALL elements after a converged solution is obtained */
  void MINI_increment(ELEMENT *elem,
              const int nelem,
              NODE *node,
              const int nnodes,
              const int ndofn,
              SUPP sup,
              EPS *eps,
              SIG *sig,
              const HOMMAT *hommat,
              const double *sol,
              const MPI_Comm mpi_comm,
              const int mp_id);

  void MINI_check_resid(const int ndofn,
            const int ne,
            const ELEMENT *elem,
            const NODE *node,
            const HOMMAT *hommat,
            const EPS *eps,
            const SIG *sig,
            const double *d_r,
            SUPP sup,
            const double *RR,
            const long *DomDof,
            const int ndofd,
            const int GDof,
            const COMMUN comm,
            const MPI_Comm mpi_comm,
            const int mp_id);

#endif /* #ifndef  MINI_ELEMENT_H */
