/* HEADER */

/**
 * AUTHORS:
 *    Matthew Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 *    Karel Matous, University of Notre Dame, <kmatous [at] nd.edu>
 *    Jaroslav Kruis, Czech Technical University
 */

#pragma once
#ifndef IN_H
#define IN_H
#include "PGFEM_io.h"
#include "element.h"
#include "node.h"
#include "PGFem3D_options.h"
#include "material.h"
#include "matgeom.h"
#include "mesh_load.h"
#include "supp.h"

/* Function reads parameters of supports */
SUPP read_supports (FILE *in,
		    long nn,
		    long ndofn,
		    NODE *node);

/* Function reads parameters of materials */
void read_material (FILE *in,
		    long nmat,
		    MATERIAL *mater,
		    int legacy);

/* Function gives stiffnesses matrix of materials */
void read_matgeom (FILE *in,
		   long nc,
		   long np,
		   MATGEOM matgeom);

/* Function reads nodal load */
void read_nodal_load (FILE *in,
		      long nln,
		      long ndofn,
		      ZATNODE *znod);

/* Function reads element surface load */
void read_elem_surface_load (FILE *in,
			     long nle_s,
			     long ndofn,
			     ELEMENT *elem,
			     ZATELEM *zele_s);

/** Function reads a specified file from the command line options to
    override the prescribed displacements used on the first step of
    the simulation. A non-zero value is returned if there is an error
    reading the file. */
int override_prescribed_displacements(SUPP sup,
				      const PGFem3D_opt *opt);
#endif
