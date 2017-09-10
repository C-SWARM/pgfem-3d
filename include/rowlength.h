/* HEADER */
#pragma once
#ifndef PGFEM3D_ROWLENGTH_H
#define PGFEM3D_ROWLENGTH_H

#include "data_structure.h"
#include "element.h"
#include "node.h"

/** Computes lengths of the rows. */
void rowlength (long *adr,
                long ne,
                long ndofn,
                long ndof,
                NODE *node,
                Element *elem,
                long gr4,
                const int mp_id);

#endif /* #define PGFEM3D_ROWLENGTH_H */
