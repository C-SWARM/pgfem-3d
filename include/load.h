/*******************************************
 * Program FEM3d ver. 2.0                  *
 * FEM - 3D analysis                       *
 * CTU, Department of Structural Mechanics *
 * Karel Matous & Jaroslav Kruis           *
 *******************************************/

/*****************/
/* November 2000 */
/*****************/

#ifndef LOAD_H
#define LOAD_H

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef SUPP_H
#include "supp.h"
#endif

#ifndef MATGEOM_H
#include "matgeom.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef MESH_LOAD_H
#include "mesh_load.h"
#endif

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef BOUNDING_ELEMENT_H
#include "bounding_element.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifndef CRPL_H
#include "crpl.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

long* compute_times_load (FILE *in1,
			  long nt,
			  long nlod_tim);

void load_vec_node (double *f,
		    long nln,
		    long ndofn,
		    ZATNODE *znode,
		    NODE *node);

int load_vec_node_defl (double *f,
			long ne,
			long ndofn,
			ELEMENT *elem,
			BOUNDING_ELEMENT *b_elems,
			NODE *node,
			HOMMAT *hommat,
			MATGEOM matgeom,
			SUPP sup,
			long npres,
			double nor_min,
			SIG *sig,
			EPS *eps,
			double dt,
			CRPL *crpl,
			double stab,
			double *r,
			const PGFem3D_opt *opts);

void load_vec_elem_sur (double *f,
			long nle_s,
			long ndofn,
			ELEMENT *elem,
			ZATELEM *zele_s);

#endif
