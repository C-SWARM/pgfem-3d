/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 *  Karel Matous, University of Notre Dame, <kmatous [at] nd.edu>
 */
#pragma once
#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

#include "PGFem3D_data_structure.h"
#include "PGFEM_mpi.h"
#include "element.h"
#include "matgeom.h"
#include "hommat.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "cohesive_element.h"
#include "bounding_element.h"
#include "PGFem3D_options.h"
#include "pgfem_comm.h"

#include "hypre_global.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

double Newton_Raphson_test(const int print_level,
                           GRID *grid,
                           MATERIAL_PROPERTY *mat,
                           FIELD_VARIABLES *variables,
                           SOLVER_OPTIONS *sol,
                           LOADING_STEPS *load,
                           COMMUNICATION_STRUCTURE *com, 
                           CRPL *crpl, /**< Crystal plasticity stuff */
                           double GNOR, /**< should be local variable. */
                           double nor1, /**< should be local variable. */
                           long nt, /**< _DEPRECATED_ */
                           MPI_Comm mpi_comm,
                           const double VVolume, /**< original volume of the domain */
                           const PGFem3D_opt *opts /**< structure of options */);
			 
  /**
   * \brief Newton-Raphson solution algorithm.
   *
   * This is the Newton-Raphson nonlinear solution algorithm. It
   * includes a globally convergent line-search algorithm and
   * subdivision procedure. If at any point the linear solution does
   * not converge to the specified tolerance (+ an additional
   * tolerance for "wiggle room"), the load is subdivided and we
   * attempt to solve the prescribed deformation in multiple
   * increments. The prescribed load is likewise subdivided if the
   * solution is not obtained within the specified number of nonlinear
   * iterations or if the line-search algorithm fails to converge. The
   * relative residual is checked for convergence. The error norm is
   * also checked for convergence to ensure that the solution
   * progresses in the case that the relative residual is not less
   * than the tolerance but the total residual is very small (err^2).
   *
   * \return time in linear solver (seconds).
   *
   * Side effects:
   *  Multiple calls to linear solver.
   *  Multiple calls to collective communication routines.
   *  Abort on unrecoverable errors (e.g., too many substeps)
   */   
  double Newton_Raphson (const int print_level,
			 int *n_step, /**< returns the number of nonlinear steps taken to solve the given increment */
			 long ne, /**< number of ELEMENT */
			 int n_be, /**< number of BOUNDING_ELEMENT */
			 long nn, /**< number of NODE */
			 long ndofn, /**< number of DOFs/node */
			 long ndofd, /**< number of DOFs on the domain */
			 long npres, /**< number of pressure DOFs/element */
			 long tim, /**< index of the current time in times */
			 double *times, /**< array of times in load history. Modified for subdivision */
			 double nor_min, /**< nonlinearr convergence tolerance */
			 double dt, /**< time increment (times[tim] - times[tim-1]) */
			 ELEMENT *elem, /**< list of volumetric elements */
			 BOUNDING_ELEMENT *b_elems, /**< _DEPRECATED_ list of bounding elements */
			 NODE *node, /**< list of nodes */
			 SUPP sup, /**< list of Dirichlet boundary conditions */
			 double *sup_defl, /**< sum of Dirichlet BC increments to step n */
			 HOMMAT *hommat, /**< list of material properites */
			 MATGEOM matgeom, /**< information related to material geometry (for crystal plasticity) */
			 SIG *sig_e, /**< stress(es) in each ELEMENT */
			 EPS *eps, /**< strain(s) in each ELEMENT */ 
			 int *Ap, /**< n_cols in each owned row of global stiffness matrix */
			 int *Ai, /**< column ids for owned rows of global stiffness matrix */
			 double *r, /**< local total solution vector to time n */
			 double *f, /**< workspace for local residual */
			 double *d_r, /**< workspace for local increment of the solution n->n+1 */
			 double *rr, /**< workspace for local _iterative_ increment of the solution */
			 double *R, /**< [in] vector of Neumann loads */
			 double *f_defl, /**< workspace for the load vector due to derichlet conditions */
			 double *RR, /**< [out] total Neumann load */
			 double *f_u, /**< workspace for load due to body force */
			 double *RRn,/**< [in] Neumann load to time n */
			 CRPL *crpl, /**< Crystal plasticity stuff */
			 double stab, /**< stabilization parameter (for -st branch) */
			 long nce, /**< number of COEL (cohesive elements) */
			 COEL *coel, /**< list of cohesive elements */
			 long FNR, /**< "Full Newton-Raphson" == 0, only compute tangent on 1st iteration */
			 double *pores, /**< [out] opening volume of failed cohesive interfaces */
			 PGFEM_HYPRE_solve_info *PGFEM_hypre, /**< custom HYPRE solver object */
			 double *BS_x, /**< workspace for the locally owned part of the global solution 'rr'. */
			 double *BS_f, /**< Global part of 'f'. */
			 double *BS_RR, /**< Global part of 'RR'. */
			 double gama, /**< related to linesearch, but is modified internally... */
			 double GNOR, /**< should be local variable. */
			 double nor1, /**< should be local variable. */
			 double err, /**< linear solve tolerance */
			 double *BS_f_u, /**< Global part of 'f_u'. */
			 long *DomDof, /**< number of global DOFs on each domain */
			 COMMUN comm, /**< sparse communication structure */
			 int GDof, /**< maximum id of locally owned global DOF */
			 long nt, /**< _DEPRECATED_ */
			 long iter_max, /**< maximum number of nonlinear iterations */
			 double *NORM, /**< [out] residual of first iteration (tim = 0, iter = 0). */
			 long nbndel, /**< number of ELEMENT on the communication boundary */
			 long *bndel, /**< ELEMENT ids on the communication boundary */
			 MPI_Comm mpi_comm,
			 const double VVolume, /**< original volume of the domain */
			 const PGFem3D_opt *opts, /**< structure of options */
			 void *microscale, /**< Container of microscale information. */
			 double alpha_alpha, /**< mid_point_rule alpha */
			 double *r_n, /**< local total solution vector to times[tim-1] */
			 double *r_n_1 /**< local total solution vector to times[tim-2] */			 
			 );

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef NEWTON_RAPHSON_H */
