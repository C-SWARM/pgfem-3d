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

#include <stdlib.h>

/* Function reads parameters of supports */
SUPP read_supports (FILE *in,
		    long nn,
		    long ndofn,
		    NODE *node);

/**
 * Read material property listing for material mat_id [0,nmat).
 * 
 * \param[in] in, File to read from
 * \param[in] mat_id, Index to 'mater'
 * \param[in,out] mater, Allocated material list
 * \param[in] legacy, format flag for how to read the listing
 *
 * \return non-zero on internal error, i.e., I/O errors.
 */
int read_material (FILE *in,
                   const size_t mat_id,
                   MATERIAL *mater,
                   const int legacy);

/**
 * Read and override material properties as specified in the options.
 *
 * If a file containing material property overrides is specified on
 * the command line, parse the file and override the specified
 * material properties.
 *
 * \param[in] nmat, number of material properties allocated in 'mater'
 * \param[in] opt, program options container
 * \param[in,out] mater, list of material properties
 *
 * \return non-zero on internal error, i.e., I/O or bounds.
 */
int override_material_properties(const size_t nmat,
                                 const PGFem3D_opt *opt,
                                 MATERIAL *mater);


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
