/* HEADER */
#pragma once
#ifndef SET_FINI_DEF
#define SET_FINI_DEF

#include "element.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"

  /** */
  void set_fini_def (long ne,
             long npres,
             ELEMENT *elem,
             EPS *eps,
             SIG *sig,
             const int analysis);

  void set_fini_def_pl (long ne,
            long npres,
            ELEMENT *elem,
            EPS *eps,
            SIG *sig,
            CRPL *crpl,
            const int analysis,
            const int plc);

#endif /* #ifndef SET_FINI_DEF */
