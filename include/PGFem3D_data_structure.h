//////////////////////////////////////////////////////////////////////
/// Define the PGFem3D data structs
/// 
/// Authors:
///   Sangmin Lee, University of Notre Dame <slee43 [at] nd.edu>
//////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _H_PGFEM3D_DATA_STRUCTURE_H_
#define _H_PGFEM3D_DATA_STRUCTURE_H_

#include "node.h"
#include "element.h"
#include "cohesive_element.h"
#include "bounding_element.h"
#include "sig.h"
#include "eps.h"
#include "matgeom.h"
#include "hommat.h"
#include "hypre_global.h"
#include "pgfem_comm.h"

/// Time stepping struct
/// Has time stepping information
typedef struct {
  int nt;        /// total number of times
  double *times; /// list of of time
} PGFem3D_TIME_STEPPING;

/// Mesh  
/// Has all mesh data
typedef struct {
  long ne;                   /// number of elements
  long nn;                   /// number of nodes
  int n_be;                  /// number of bounding elements
  long nce;                  /// number of COEL (cohesive elements) 
  NODE *node;                /// list of node
  ELEMENT *element;          /// list of element
  BOUNDING_ELEMENT *b_elems; /// list of bounding element
  COEL *coel;                /// list of cohesive elements 
} GRID;

/// struct for field variables
typedef struct {
  long ndofn;     /// number of degree of freedom on a node
  long ndofd;     /// number of degree of freedom in the domain
  long npres;     /// number of pressure per element
  double *u_np1;  /// displacement at n+1
  double *u_n;    /// displacement at n
  double *u_nm1;  /// displacement at n-1
  double *d_u;    /// workspace for local increment of the solution n->n+1 
  double *dd_u;   /// workspace for local _iterative_ increment of the solution 
  double *f;      /// workspace for local residual 
  double *R;      /// [in] vector of Neumann loads 
  double *f_defl; /// workspace for the load vector due to derichlet conditions 
  double *RR;     /// [out] total Neumann load 
  double *f_u;    /// workspace for load due to body force 
  double *RRn;    /// [in] Neumann load to time n 
  double *pores;  /// [out] opening volume of failed cohesive interfaces 
  double *BS_x;   /// workspace for the locally owned part of the global solution 'rr'. 
  double *BS_f;   /// Global part of 'f'. 
  double *BS_f_u; /// Global part of 'f_u'.
  double *BS_RR;  /// Global part of 'RR'. 
  double *NORM;   /// [out] residual of first iteration (tim = 0, iter = 0).
  SIG *sig;       /// pointer for the stress
  EPS *eps;       /// pointer for strain
} FIELD_VARIABLES;

/// struct for material properties
typedef struct {
  HOMMAT *hommat;  /// list of material properites 
  MATGEOM matgeom; /// information related to material geometry (for crystal plasticity) 
} MATERIAL_PROPERTY;

/// struct for solution scheme
typedef struct {
  long tim;         /// current time step number
  double *times;    /// list of time
  double dt_n;      /// dt at n
  double dt_np1;    /// dt at n+1
  int *n_step;      /// 
  double nor_min;   /// 
  long iter_max;    /// maximum number of iterations
  double alpha;     /// midpoint rule alpha
  void *microscale;
  double stab;      /// stabilization parameter (for -st branch) 
  PGFEM_HYPRE_solve_info *PGFEM_hypre; /// custom HYPRE solver object 
  long FNR;         /// "Full Newton-Raphson" == 0, only compute tangent on 1st iteration 
  double gama;      /// related to linesearch, but is modified internally... 
  double err;       /// linear solve tolerance 
} SOLVER_OPTIONS;

/// struct for the boundary conditions
typedef struct {
  SUPP sup;
  double *sup_defl;
} LOADING_STEPS;

/// struct for the communication
typedef struct {
  int *Ap;      /// n_cols in each owned row of global stiffness matrix 
	int *Ai;      /// column ids for owned rows of global stiffness matrix 
  long *DomDof; /// number of global DOFs on each domain 
  long nbndel;  /// number of ELEMENT on the communication boundary 
  long *bndel;  /// ELEMENT ids on the communication boundary 
  COMMUN comm;  /// sparse communication structure 
  int GDof;     /// maximum id of locally owned global DOF 
} COMMUNICATION_STRUCTURE;


#endif