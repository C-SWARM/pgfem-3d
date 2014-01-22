/* HEADER */
#ifndef INITIALIZE_DAMAGE_H
#define INITIALIZE_DAMAGE_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

  /** Initialize the damage part of the strain object */
  void initialize_damage(const int ne,
			 const ELEMENT *elem,
			 const HOMMAT *hommat,
			 EPS *eps,
			 const int analysis);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef INITIALIZE_DAMAGE_H */
