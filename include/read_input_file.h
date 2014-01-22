/* HEADER */
#ifndef READ_INPUT_FILE_H
#define READ_INPUT_FILE_H

#include "PGFEM_mpi.h"

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef MATERIAL_H
#include "material.h"
#endif

#ifndef MATGEOM_H
#include "matgeom.h"
#endif

#ifndef SUPP_H
#include "supp.h"
#endif

#ifndef MESH_LOAD_H
#include "mesh_load.h"
#endif

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
