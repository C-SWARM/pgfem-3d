/* HEADER */
#ifndef _PGFEM_TO_VTK_HPP_
#define _PGFEM_TO_VTK_HPP_

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef SUPP_H
#include "supp.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

#ifndef INTERFACE_MACRO_H
#include "interface_macro.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

/** Create a vtkUnstructuredGrid object from a PGFem3D mesh and
    associated information. ***CURRENTLY ONLY SUPPORT LINEAR TETRAS
    FOR DAMAGE*** */
void* PGFem3D_to_vtkUnstructuredGrid(const int nnode,
				     const int nelems,
				     const NODE *nodes,
				     const ELEMENT *elems,
				     const SUPP supports,
				     const SIG *stress,
				     const EPS *strain,
				     const double *dofs,
				     const int analysis_type);

void PGFem3D_destroy_vtkUnstructuredGrid(void *grid);

/** Write a vtkUnstructuredGrid object to a specified file in binary
    or ASCII format. It does not matter how the VTK object is
    constructed, it will be printed properly. */
void PGFem3D_write_vtkUnstructuredGrid(const char* filename,
				       const void *grid,
				       const int ascii);
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
