/**
 * AUTHORS:
 * Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 * Karel Matous, University of Notre Dame, kmatous [at] nd.edu
 */
#ifndef PGFEM3D_LOAD_H
#define PGFEM3D_LOAD_H

#include "PGFem3D_options.h"
#include "bounding_element.h"
#include "cohesive_element.h"
#include "crpl.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "macro_micro_functions.h"
#include "matgeom.h"
#include "mesh_load.h"
#include "node.h"
#include "sig.h"
#include "supp.h"

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
                    const Node *node,
                    const int mp_id);

/**
 * \brief Compute the load vector from the elements with surface load
 * [NOT IMPLEMENTED].
 */
void load_vec_elem_sur (double *f,
                        const long nle_s,
                        const long ndofn,
                        const Element *elem,
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

/// Multiscale simulation interface to compute load vector due to BCs(Dirichlet)
///
/// s->f_defl will be updated
///
/// \param[in] c structure of macroscale information
/// \param[in,out] s contains the information for the history-dependent solution
/// \param[in] opts structure PGFem3D option
/// \param[in] nor_min nonlinear convergence tolerance
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int compute_load_vector_for_prescribed_BC_multiscale(COMMON_MACROSCALE *c,
                                                     MACROSCALE_SOLUTION *s,
                                                     const PGFem3D_opt *opts,
                                                     double nor_min,
                                                     int myrank);

#endif // #define PGFEM3D_LOAD_H
