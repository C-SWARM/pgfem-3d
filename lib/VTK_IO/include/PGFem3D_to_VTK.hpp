/* HEADER */
#pragma once
#ifndef _PGFEM_TO_VTK_HPP_
#define _PGFEM_TO_VTK_HPP_

// pre include mpi. This is necessary so that there is not a conflict
// with c/c++ when wrapped with extern c
#include "mpi.h"

#include "element.h"
#include "node.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "enumerations.h"
#include "interface_macro.h"

/** Create a vtkUnstructuredGrid object from a PGFem3D mesh and
    associated information. ***CURRENTLY ONLY SUPPORT LINEAR TETRAS
    FOR DAMAGE*** */
void* PGFem3D_to_vtkUnstructuredGrid(const int nnode,
                     const int nelems,
                     const NODE *nodes,
                     const Element *elems,
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

int read_VTK_file(char fn[], double *r);
int read_VTK_file4TF(char fn[], double *r, double *P, double *V);

#endif
