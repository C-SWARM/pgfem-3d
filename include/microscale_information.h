/* HEADER */
/**
 * \file This header defines the structure(s) for the microscale solutions.
 *
 * AUTHORS:
 *    Matthew Mosby, University of Notre Dame, <mmosby1 (at) nd.edu>
 */
#pragma once
#ifndef MICROSCALE_INFORMATION_H
#define MICROSCALE_INFORMATION_H

#include "PGFEM_mpi.h"
#include "element.h"
#include "node.h"
#include "pgfem_comm.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "matgeom.h"
#include "hommat.h"
#include "cohesive_element.h"
#include "PGFem3D_options.h"
#include "hypre_global.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** This is the structure of microscale information that is
      identical for all microstructures. */
  typedef struct COMMON_MICROSCALE{
    /* options /solver information */
    PGFEM_HYPRE_solve_info *SOLVER;
    double lin_err;
    double lim_zero;
    int maxit_nl;
    int *Ap;
    int *Ai;

    /* communication information */
    COMMUN pgfem_comm;
    MPI_Comm mpi_comm;
    long nbndel; /* no. elements on comm boundary */
    long *bndel; /* bnd elem ids */
    long ndofd; /* dof on dom */
    long *DomDof; /* global dof on dom */
    int GDof; /* first global dof id on dom */

    /* mesh info */
    long nn; /**< no. nodes */
    long ne; /**< no. elements */
    long nce; /**< no. cohesive elements */
    long ndofn; /**< no. dof on each node */
    long npres; /**< no. pressure nodes on elem */
    double VVolume; /**< original volume */
    NODE *node;
    ELEMENT *elem; /* NOTE: state/solution information is copied from
		      solution structure */
    COEL *coel; /* NOTE: state/solution information is copied from
		   solution structure */
    long n_orient;
    MATGEOM matgeom; /**< !pointer */
    long nhommat;
    HOMMAT *hommat;
    ENSIGHT ensight; /**< !pointer */
    SUPP supports; /**< !pointer */
    int n_co_props;
    cohesive_props *co_props;

    /* mixed tangent matrices. */
    void *K_01;
    void *K_10;

    /* buffer for solution stuff */
    void *solution_buffer;

  } COMMON_MICROSCALE;

  /** This structure contains the information for the
      history-dependent solution. */
  typedef struct MICROSCALE_SOLUTION{
    /* stress/strain/state information */
    SIG *sig_e;
    SIG *sig_n;
    EPS *eps;
    CRPL *crpl;
    long npres;

    /* solution information */
    double *r; /**< current solution r at macro time n+1 */

    /* State variables at time n */
    size_t packed_state_var_len;
    char *packed_state_var_n;

    /** The following are pointers to a shared buffer elsewhere! The
	buffers are used purely as a workspace. */
    /* local vectors */
    double *f;
    double *d_r;
    double *rr;
    double *D_R;
    double *R;
    double *RR;
    double *f_u;
    double *f_defl;
    double *RRn;
    double *U;
    double *DK;
    double *dR;

    /* global vectors */
    double *BS_f;
    double *BS_x;
    double *BS_RR;
    double *BS_f_u;
    double *BS_d_r;
    double *BS_D_R;
    double *BS_rr;
    double *BS_R;
    double *BS_U;
    double *BS_DK;
    double *BS_dR;

    /* convergence info */
    double dt;
    double *times;
    int tim;
    int p_tim;
    double NORM;
  } MICROSCALE_SOLUTION;

  typedef struct{
    size_t size;
    int *map;
  } sol_idx_map;
  void sol_idx_map_build(sol_idx_map *map,
			 const size_t size);
  void sol_idx_map_destroy(sol_idx_map *map);
  void sol_idx_map_sort_id(sol_idx_map *map);
  void sol_idx_map_sort_idx(sol_idx_map *map);
  int sol_idx_map_id_get_idx(const sol_idx_map *map,
			     const int id);
  int sol_idx_map_idx_get_id(const sol_idx_map *map,
			     const int idx);
  void sol_idx_map_idx_set_id(sol_idx_map *map,
			      const int idx,
			      const int id);
  /** structure to contain all microscale information */
  typedef struct MICROSCALE{
    PGFem3D_opt *opts;
    COMMON_MICROSCALE *common;
    sol_idx_map idx_map;
    MICROSCALE_SOLUTION *sol;
  } MICROSCALE;

  /** Instantiate a MICROSCALE. Some space is allocated, but full
      allocation is deffered to the build_MICROSCALE
      function. Unallocated information is set to NULL/default
      values. */
  void initialize_MICROSCALE(MICROSCALE **microscale);

  /** build the full MICROSCALE given by the list of command-line
      style options */
  void build_MICROSCALE(MICROSCALE *microscale,
			MPI_Comm mpi_comm,
			const int argc,
			char **argv);

  /** build n solutions to compute on the scale */
  void build_MICROSCALE_solutions(const int n_solutions,
				  MICROSCALE *microscale);

  void destroy_MICROSCALE(MICROSCALE *microscale);

  /** resets a single MICROSCALE_SOLUTION to time (n) */
  int reset_MICROSCALE_SOLUTION(MICROSCALE_SOLUTION *sol,
				const MICROSCALE *micro);

  /** updates a single MICROSCALE_SOLUTION to time (n+1) */
  int update_MICROSCALE_SOLUTION(MICROSCALE_SOLUTION *sol,
				 const MICROSCALE *micro);

  /**=== Aliases for MACROSCALE ===*/
  typedef COMMON_MICROSCALE COMMON_MACROSCALE;
  typedef MICROSCALE_SOLUTION MACROSCALE_SOLUTION;
  typedef MICROSCALE MACROSCALE;
#define initialize_MACROSCALE(macro) initialize_MICROSCALE(macro)
#define build_MACROSCALE(macro,comm,argc,argv)	\
  build_MICROSCALE(macro,comm,argc,argv)
#define build_MACROSCALE_solution(macro) \
  build_MICROSCALE_solutions(1,macro)
#define destroy_MACROSCALE(macro) destroy_MICROSCALE(macro)

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */

