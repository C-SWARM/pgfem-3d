#ifndef PGFEM3D_RES_FINI_DEF_H
#define PGFEM3D_RES_FINI_DEF_H

#include "crpl.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "sig.h"

/** */
void res_fini_def (long ne,
                   long npres,
                   Element *elem,
                   EPS *eps,
                   SIG *sig,
                   CRPL *crpl,
                   const int analysis);

#endif /* #define PGFEM3D_RES_FINI_DEF_H */
