#ifndef ROWLENGTH_H
#define ROWLENGTH_H

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Computes lengths of the rows. */
  void rowlength (long *adr,
		  long ne,
		  long ndofn,
		  long ndof,
		  NODE *node,
		  ELEMENT *elem,
		  long gr4);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef ROWLENGTH_H */
