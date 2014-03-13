/* HEADER */

#pragma once
#ifndef LOAD_H
#define LOAD_H

#include "element.h"
#include "node.h"
#include "supp.h"
#include "matgeom.h"
#include "hommat.h"
#include "mesh_load.h"
#include "cohesive_element.h"
#include "bounding_element.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "PGFem3D_options.h"

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
