/* HEADER */
#pragma once
#ifndef PGFEM3D_SET_FINI_DEF
#define PGFEM3D_SET_FINI_DEF

#include "crpl.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "sig.h"

/** */
void set_fini_def (long ne,
                   long npres,
                   Element *elem,
                   EPS *eps,
                   SIG *sig,
                   const int analysis);

void set_fini_def_pl (long ne,
                      long npres,
                      Element *elem,
                      EPS *eps,
                      SIG *sig,
                      CRPL *crpl,
                      const int analysis,
                      const int plc);

#endif /* #define PGFEM3D_SET_FINI_DEF */
