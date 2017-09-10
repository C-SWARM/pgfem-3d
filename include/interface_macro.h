/**
 * \file functions to compute the macroscopic deformation gradient and
 *  macroscopic displacements
 *
 * AUTHORS:
 *    Matthew Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */
#ifndef PGFEM3D_INTERFACE_MACRO_H
#define PGFEM3D_INTERFACE_MACRO_H

#include "node.h"
#include "supp.h"

/** read the macroscopic normal and cell thickness (lc) from a file
    called "normal.in" placed in the input directory.
    Format: [lc Nx Ny Nz] */
int read_interface_macro_normal_lc(char *in_dir,
                                   SUPP sup);

/** computes the macroscopic jump across the interface. The jump is
    computed according to the first six (6) prescribed displacements
    assuming the order [u+ v+ w+ u- v- w-]. If there are not at least
    six prescribed displacements, the jump is filled with zeros and an
    error is returned. */
int compute_interface_macro_jump_u(double *jump_u,
                                   const SUPP sup,
                                   const int analysis);

/** Compute the macroscopic interface deformation gradient
    contribution grad(^0u) */
int compute_interface_macro_grad_u(double *F_0,
                                   const double lc,
                                   const double *jump_u,
                                   const double *normal);

/** Compute the macroscopic displacements at a particular node */
int compute_interface_macro_disp_at_node(double *u_0,
                                         const Node *ptrNode,
                                         const double *F_0,
                                         const int analysis);

/** Create the maacroscopic grad(u) for either interfaces or
    bulk. Bulk vs. interface is determined by the number of
    prescribed displacements: >= 9 --> bulk, >= 8 --> interface */
int compute_macro_grad_u(double *F0,
                         const SUPP sup,
                         const int analysis);

#endif /* #define PGFEM3D_INTERFACE_MACRO_H */
