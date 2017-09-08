/* HEADER */
#pragma once
#ifndef ROWLENGTH_H
#define ROWLENGTH_H

#include "data_structure.h"
#include "element.h"
#include "node.h"

  /** Computes lengths of the rows. */
  void rowlength (long *adr,
          long ne,
          long ndofn,
          long ndof,
          NODE *node,
          ELEMENT *elem,
          long gr4,
          const int mp_id);

#endif /* #ifndef ROWLENGTH_H */
