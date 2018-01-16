/* HEADER */
#ifndef PGFEM3D_PGFEM3D_TO_VTK_HPP
#define PGFEM3D_PGFEM3D_TO_VTK_HPP

#include "element.h"
#include "enumerations.h"
#include "eps.h"
#include "interface_macro.h"
#include "node.h"
#include "sig.h"
#include "supp.h"

/** Create a vtkUnstructuredGrid object from a PGFem3D mesh and
    associated information. ***CURRENTLY ONLY SUPPORT LINEAR TETRAS
    FOR DAMAGE*** */
void* PGFem3D_to_vtkUnstructuredGrid(const int nnode,
                                     const int nelems,
                                     const Node *nodes,
                                     const Element *elems,
                                     const SUPP supports,
                                     const SIG *stress,
                                     const EPS *strain,
                                     const double *dofs,
                                     const double *Pnp1,
                                     const double *Vnp1,                                                                          
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

#endif // #define PGFEM3D_PGFEM3D_TO_VTK_HPP
