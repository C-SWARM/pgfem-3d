/* HEADER */
/**
 * AUTHORS:
 * Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 * Karel Matous, University of Notre Dame, kmatous [at] nd.edu
 */

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

/**
 * \brief Get the list of times to increment the load from the file.
 */
long* compute_times_load (FILE *in1,
			  const long nt,
			  const long nlod_tim);

/**
 * \brief Get the load from the nodes with prescribed force.
 */
void load_vec_node (double *f,
		    const long nln,
		    const long ndofn,
		    const ZATNODE *znode,
		    const NODE *node,
		    const int mp_id);

/**
 * \brief Compute the load vector from the prescribed boundary
 * conditions.
 */
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
			double *r_n,
			const PGFem3D_opt *opts,double alpha,
			const int mp_id);

/**
 * \brief Compute the load vector from the elements with surface load
 * [NOT IMPLEMENTED].
 */
void load_vec_elem_sur (double *f,
			const long nle_s,
			const long ndofn,
			const ELEMENT *elem,
			const ZATELEM *zele_s);
			
/// Compute load vector for prescribed BCs(Dirichlet)
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv field variable object
/// \param[in] sol solution scheme object
/// \param[in] load object for loading
/// \param[in] dt time step
/// \param[in] crpl object for lagcy crystal plasticity
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[in] mp_id mutiphysics id
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int compute_load_vector_for_prescribed_BC(GRID *grid,
                                          MATERIAL_PROPERTY *mat,
                                          FIELD_VARIABLES *fv,
                                          SOLVER_OPTIONS *sol,
                                          LOADING_STEPS *load,
                                          double dt,                                          
                                          CRPL *crpl,
                                          const PGFem3D_opt *opts,
                                          MULTIPHYSICS *mp,
                                          int mp_id,
                                          int myrank);
#endif
