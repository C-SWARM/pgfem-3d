/**
 * @file This file declares the functions for implimenting the three
 * field Hu-Washizu extension of the MINI element.
 *
 * AUTHORS:
 * Matthew Mosby
 */
#ifndef PGFEM3D_MINI_3F_ELEMENT_H
#define PGFEM3D_MINI_3F_ELEMENT_H

#include "PGFem3D_data_structure.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "node.h"
#include "sig.h"
#include "supp.h"

/** Reset the internal variables for the MINI_3f element */
void MINI_3f_reset(Element *elem,
                   const int nelem,
                   const int npres,
                   const int nvol,
                   SIG *sig,
                   EPS *eps);

/** Compute the stiffness matrix on the element */
int MINI_3f_stiffmat_el(double *Ks,            /**< Element stiffmat */
                        const int ii,          /**< element id */
                        const int ndofn,
                        const int nne,
                        const double *x,
                        const double *y,
                        const double *z,
                        const Element *elem,
                        const HOMMAT *hommat,
                        const long *nod,
                        const Node *node,
                        const EPS *eps,
                        const SIG *sig,
                        const double *r_e);    /**< dof values on elem */

int MINI_3f_resid_el(double *Res,         /**< Element residual */
                     const int ii,        /**< element id */
                     const int ndofn,
                     const int nne,
                     const double *x,
                     const double *y,
                     const double *z,
                     const Element *elem,
                     const long *nod,
                     const Node *node,
                     const HOMMAT *hommat,
                     const EPS *eps,
                     const SIG *sig,
                     const double *r_e);    /**< dof values on elem */

int MINI_3f_update_bubble_el(Element *elem,
                             const int ii, /* id of element working on */
                             const int nne,
                             const Node *node,
                             const int ndofn,
                             const double *x,
                             const double *y,
                             const double *z,
                             const long *nod, /* list of node ids on elem */
                             const EPS *eps,
                             const SIG *sig,
                             const HOMMAT *hommat,
                             const double *sol_e, /* accum. incr on elem */
                             const double *dsol_e); /*increment on element*/

int MINI_3f_update_bubble(Element *elem,
                          const int nelem,
                          const Node *node,
                          const int ndofn,
                          const SUPP sup,
                          const EPS *eps,
                          const SIG *sig,
                          const HOMMAT *hommat,
                          const double *sol, /* accum. solution  on incr */
                          const double *dsol, /* sol from current iter */
                          const int iter,
                          const int mp_id);

void MINI_3f_increment_el(Element *elem,
                          const int ii, /* id of element working on */
                          const int nne,
                          const Node *node,
                          const long *nod,
                          const int ndofn,
                          const double *x,
                          const double *y,
                          const double *z,
                          EPS *eps,
                          SIG *sig,
                          const HOMMAT *hommat,
                          const double *sol_e);

void MINI_3f_increment(Element *elem,
                       const int nelem,
                       Node *node,
                       const int nnodes,
                       const int ndofn,
                       const SUPP sup,
                       EPS *eps,
                       SIG *sig,
                       const HOMMAT *hommat,
                       const double *sol,
		       const pgfem3d::CommunicationStructure *com,
                       const int mp_id);

void MINI_3f_check_resid(const int ndofn,
                         const int ne,
                         const Element *elem,
                         const Node *node,
                         const HOMMAT *hommat,
                         const EPS *eps,
                         const SIG *sig,
                         const double *d_r,
                         const SUPP sup,
                         const double *RR,
                         const long *DomDof,
                         const int ndofd,
			 const pgfem3d::CommunicationStructure *com,
                         const int mp_id);

#endif // #define PGFEM3D_MINI_3F_ELEMENT_H
