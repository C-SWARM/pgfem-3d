/* HEADER */
/**
 * AUTHORS:
 *    Matthew Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
 */

#pragma once
#ifndef INITIALIZE_DAMAGE_H
#define INITIALIZE_DAMAGE_H

#include "element.h"
#include "hommat.h"
#include "eps.h"

  /** Initialize the damage part of the strain object */
  void initialize_damage(const int ne,
             const ELEMENT *elem,
             const HOMMAT *hommat,
             EPS *eps,
             const int analysis);

#endif /* #ifndef INITIALIZE_DAMAGE_H */
