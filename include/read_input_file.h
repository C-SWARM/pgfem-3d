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
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
