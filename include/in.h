/*******************************************
 * Program FEM3d ver. 2.0                  *
 * FEM - 3D analysis                       *
 * CTU, Department of Structural Mechanics *
 * Karel Matous & Jaroslav Kruis           *
 *******************************************/

/*****************/
/* November 2000 */
/*****************/

#ifndef IN_H
#define IN_H

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifndef MATERIAL_H
#include "material.h"
#endif

#ifndef MATGEOM_H
#include "matgeom.h"
#endif

#ifndef MESH_LOAD_H
#include "mesh_load.h"
#endif

#ifndef SUPP_H
#include "supp.h"
#endif

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
