/* HEADER */
#pragma once
#ifndef READ_INPUT_FILE_H
#define READ_INPUT_FILE_H

#include "PGFEM_mpi.h"
#include "PGFem3D_options.h"
#include "node.h"
#include "element.h"
#include "material.h"
#include "matgeom.h"
#include "supp.h"
#include "mesh_load.h"
#include "PGFem3D_data_structure.h"
#include "Arc_length.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Function for reading the entire input file. All required space
      is allocated within this function */
  int read_input_file(const PGFem3D_opt *opts,
		      MPI_Comm comm,
		      long *nn,
		      long *Gnn,
		      long *ndofn,
		      long *ne,
		      long *lin_maxit,
		      double *lin_err,
		      double *lim_zero,
		      long *nmat,
		      long *n_concentrations,
		      long *n_orient,
		      NODE **node,
		      ELEMENT **elem,
		      MATERIAL **material,
		      MATGEOM *matgeom,
		      SUPP *sup,
		      long *nln,
		      ZATNODE **znod,
		      long *nel_s,
		      ZATELEM **zelem_s,
		      long *nel_v,
		      ZATELEM **zelem_v);
		      
  /// read mesh file
  ///
  /// \param[out] grid a mesh object
  /// \param[out] mat a material object
  /// \param[out] variables object for field variables
  /// \param[out] sol object for solution scheme
  /// \param[out] load object for loading
  /// \param[in] comm MPI_COMM_WORLD
  /// \param[in] opts structure PGFem3D option
  /// \return non-zero on internal error
  int read_mesh_file(GRID *grid, 
                     MATERIAL_PROPERTY *mat,
                     FIELD_VARIABLES *variables,
                     SOLVER_OPTIONS *sol,
                     LOADING_STEPS *load,
                     MPI_Comm mpi_comm,
                     const PGFem3D_opt *opts);

  int read_solver_file(PGFem3D_TIME_STEPPING *ts,
                       MATERIAL_PROPERTY *mat,
                       FIELD_VARIABLES *variables,
                       SOLVER_OPTIONS *sol,
                       LOADING_STEPS *load,
                       ARC_LENGTH_VARIABLES *arc,
                       CRPL *crpl,
                       const PGFem3D_opt *opts,
                       int myrank);

int read_initial_values(GRID *grid,
                        MATERIAL_PROPERTY *mat,
                        FIELD_VARIABLES *fv,
                        SOLVER_OPTIONS *sol,
                        LOADING_STEPS *load,
                        PGFem3D_TIME_STEPPING *ts,
                        PGFem3D_opt *opts,
                        int *restart, 
                        double *tnm1, 
                        int myrank);
                       
int read_and_apply_load_increments(GRID *grid,
                                   FIELD_VARIABLES *variables,
                                   LOADING_STEPS *load, 
                                   long tim, 
                                   MPI_Comm mpi_comm,
                                   int myrank);

int read_cohesive_elements(GRID *grid,
                           MATERIAL_PROPERTY *mat,
                           const PGFem3D_opt *opts,
                           ENSIGHT ensight,
                           MPI_Comm mpi_comm,
                           int myrank);                           
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
