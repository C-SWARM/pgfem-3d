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
#include "material.h"
#include "matgeom.h"
#include "hommat.h"
#include "mesh_load.h"
#include "hypre_global.h"
#include "pgfem_comm.h"
#include "comm_hints.h"

/// Time stepping struct
/// Has time stepping information
typedef struct {
  long nt;          /// total number of times
  long tim;         /// current time step number
  double *times;    /// list of time
  double dt_n;      /// dt at n
  double dt_np1;    /// dt at n+1
  long *print;      /// step numbers to be printed  
} PGFem3D_TIME_STEPPING;

/// Mesh  
/// Has all mesh data
typedef struct {
  long Gnn;                  /// global number of nodes
  long Gne;                  /// global number of elements
  long Gnbndel;              /// global number of boundary elements
  long Gn_be;                /// global number of bounding elements
  long ne;                   /// number of elements
  long nn;                   /// number of nodes
  int nsd;                   /// number of spatial dimension
  int n_be;                  /// number of bounding elements
  long nce;                  /// number of COEL (cohesive elements)
  long Gnce;                 /// global number of cohesive elements
  NODE *node;                /// list of node
  ELEMENT *element;          /// list of element
  BOUNDING_ELEMENT *b_elems; /// list of bounding element
  COEL *coel;                /// list of cohesive elements 
} GRID;

/// struct for field variables
typedef struct {
  long Gndof;   /// total number of degree freedom
  long ndofn;     /// number of degree of freedom on a node
  long ndofd;     /// number of degree of freedom in the domain
  long npres;     /// number of pressure per element
  long n_concentrations; /// number of concentrations
  double *u_np1;  /// displacement at n+1
  double *u_n;    /// displacement at n
  double *u_nm1;  /// displacement at n-1
  double *d_u;    /// workspace for local increment of the solution n->n+1 
  double *dd_u;   /// workspace for local _iterative_ increment of the solution 
  double *f;      /// workspace for local residual 
  double *R;      /// [in] vector of Neumann loads (Incramental forces) 
  double *f_defl; /// workspace for the load vector due to derichlet conditions 
  double *RR;     /// [out] total Neumann load (Total forces for sudivided increment)
  double *f_u;    /// workspace for load due to body force 
  double *RRn;    /// [in] Neumann load to time n (Total force after equiblirium)
  double pores;   /// [out] opening volume of failed cohesive interfaces 
  double *BS_x;   /// workspace for the locally owned part of the global solution 'rr'
  double *BS_f;   /// Global part of 'f' 
  double *BS_f_u; /// Global part of 'f_u'
  double *BS_RR;  /// Global part of 'RR' 
  double NORM;    /// [out] residual of first iteration (tim = 0, iter = 0).
  SIG *sig;       /// pointer for the stress
  EPS *eps;       /// pointer for strain
  SIG *sig_n;     /// smoothed stress
} FIELD_VARIABLES;

/// struct for field variables
typedef struct {
  long Gndof;   /// total number of degree freedom
  long ndofn;     /// number of degree of freedom on a node
  long ndofd;     /// number of degree of freedom in the domain
  double *T_np1;  /// displacement at n+1
  double *T_n;    /// displacement at n
  double *T_nm1;  /// displacement at n-1
  double *dT;     /// workspace for local increment of the temperature solution n->n+1 
  double *dd_T;   /// workspace for local _iterative_ increment of the solution 
  double NORM;    /// [out] residual of first iteration (tim = 0, iter = 0).
} FIELD_VARIABLES_THERMAL;

/// struct for material properties
typedef struct {
  double *density; /// list of material density
  MATERIAL *mater; /// list of material properites
  HOMMAT *hommat;  /// list of homogeneous material properites 
  MATGEOM matgeom; /// information related to material geometry (for crystal plasticity)
  long nhommat;    /// number of homogeneous materials
  long nmat;       /// number of materials
  long n_orient;   /// number of orientations
  int n_co_props;  /// number of cohesive material properites;
  cohesive_props *co_props; /// list of cohesive material properites
} MATERIAL_PROPERTY;

/// struct for solution scheme
typedef struct {
  int *n_step;      /// the number of nonlinear steps taken to solve the given increment
  double nor_min;   /// nonlinearr convergence tolerance for Newton Raphson
  long iter_max;    /// maximum number of iterations for Newton Raphson
  double alpha;     /// midpoint rule alpha
  void *microscale; /// Container of microscale information
  PGFEM_HYPRE_solve_info *PGFEM_hypre; /// HYPRE solver object 
  long FNR;         /// "Full Newton-Raphson" == 0, only compute tangent on 1st iteration 
  double gama;      /// related to linesearch, but is modified internally... 
  double err;       /// linear solve tolerance 
  long iter_max_sol;/// maximum number of iterations for linear solver 
  double computer_zero; /// computer zero   
} SOLVER_OPTIONS;

/// struct for the boundary conditions
typedef struct {
  SUPP sup;          /// list of Dirichlet boundary conditions
  double *sup_defl;  /// sum of Dirichlet BC increments to step n
  long nln;       ;  /// number of nodes with loads
  long nle_s;        /// number of surface element with loads
  long nle_v;        /// number of volume element with loads
  ZATNODE *znod;     /// list of nodes with loads
  ZATELEM *zele_s;   /// list of surface element with loads
  ZATELEM *zele_v;   /// list of volume element with loads
  long *tim_load;    /// list of time steps to be saved
  FILE *solver_file; /// file pointer for reading loads increments
} LOADING_STEPS;

/// struct for the communication
typedef struct {
  int nproc;         /// number of mpi processes
  int *Ap;           /// n_cols in each owned row of global stiffness matrix 
	int *Ai;           /// column ids for owned rows of global stiffness matrix 
  long *DomDof;      /// number of global DOFs on each domain 
  long nbndel;       /// number of ELEMENT on the communication boundary 
  long *bndel;       /// ELEMENT ids on the communication boundary 
  COMMUN comm;       /// sparse communication structure 
  int GDof;          /// maximum id of locally owned global DOF
  long NBN;          /// Number of nodes on domain interfaces
  Comm_hints *hints; /// Comm_hints structure
} COMMUNICATION_STRUCTURE;

/// for setting physics ids
enum {MULTIPHYSICS_MECHANICAL, 
      MULTIPHYSICS_THERMAL,
      MULTIPHYSICS_CHEMICAL, 
      MULTIPHYSICS_NO} MULTIPHYSICS_ANALYSIS;

/// struct for setting multiphysics
typedef struct {
  int physicsno;      /// number of physics 
  char **physicsname; /// physics names
  int *physics_ids;   /// physics ids
} MULTIPHYSICS;

/// initialize time stepping variable
///
/// \param[in, out] ts an object for time stepping
/// \return non-zero on internal error
int time_stepping_initialization(PGFem3D_TIME_STEPPING *ts);

/// destruct time stepping variable
///
/// \param[in, out] ts an object for time stepping
/// \return non-zero on internal error
int destruct_time_stepping(PGFem3D_TIME_STEPPING *ts);

/// initialize mesh object
///
/// \param[in, out] grid an object containing all mesh data
/// \return non-zero on internal error
int grid_initialization(GRID *grid);

/// destruct of mesh
///
/// \param[in, out] grid an object containing all mesh data
/// \return non-zero on internal error
int destruct_grid(GRID *grid, 
                  const PGFem3D_opt *opts,
                  MULTIPHYSICS *mp);


/// initialize field variables
/// 
/// \param[in, out] fv an object containing all field variables
/// \return non-zero on internal error
int field_varialbe_initialization(FIELD_VARIABLES *fv);

/// construct field variables
/// 
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \param[in] com an object for communication
/// \param[in] opts structure PGFem3D option
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int construct_field_varialbe(FIELD_VARIABLES *fv, 
                             GRID *grid,
                             COMMUNICATION_STRUCTURE *com,
                             const PGFem3D_opt *opts,
                             int myrank);

/// destruct field variables
/// 
/// \param[in, out] fv an object containing all field variables
/// \param[in] grid an object containing all mesh data
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int destruct_field_varialbe(FIELD_VARIABLES *fv, 
                            GRID *grid,
                            const PGFem3D_opt *opts); 

/// initialize field variables thermal part
/// 
/// \param[in, out] fv an object containing all field variables for thermal
/// \return non-zero on internal error
int thermal_field_varialbe_initialization(FIELD_VARIABLES_THERMAL *fv);

/// initialize material properties
///
/// \param[in, out] mat an object containing all material parameters
/// \return non-zero on internal error
int material_initialization(MATERIAL_PROPERTY *mat);

/// destruct material properties
///
/// \param[in, out] mat an object containing all material parameters
/// \return non-zero on internal error
int destruct_material(MATERIAL_PROPERTY *mat,
                      const PGFem3D_opt *opts);                            


/// initialize iterative solver object
///
/// \param[in, out] sol an object containing data for linear solver
/// \return non-zero on internal error
int solution_scheme_initialization(SOLVER_OPTIONS *sol);

/// initialize loading steps object
///
/// \param[in, out] load an object containing boundary increments
/// \return non-zero on internal error
int loading_steps_initialization(LOADING_STEPS *load);

/// destruct loading steps object
///
/// \param[in, out] load an object containing boundary increments
/// \return non-zero on internal error
int destruct_loading_steps(LOADING_STEPS *load);

/// initialize communication structures
///
/// \param[in, out] com an object for communication
/// \return non-zero on internal error
int communication_structure_initialization(COMMUNICATION_STRUCTURE *com);


/// destruct communication structures
///
/// \param[in, out] com an object for communication
/// \return non-zero on internal error
int destruct_communication_structure(COMMUNICATION_STRUCTURE *com);
#endif