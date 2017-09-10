/* HEADER */
/**
 * AUTHORS:
 *    Matthew Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */
#ifndef PGFEM3D_INITIALIZE_DAMAGE_H
#define PGFEM3D_INITIALIZE_DAMAGE_H

#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"

/** Initialize the damage part of the strain object */
void initialize_damage(const int ne,
                       const Element *elem,
                       const HOMMAT *hommat,
                       EPS *eps,
                       const int analysis);

#endif /* #define PGFEM3D_INITIALIZE_DAMAGE_H */
