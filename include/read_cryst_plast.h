#ifndef READ_CRYST_PLAST_H
#define READ_CRYST_PLAST_H

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef CRPL_H
#include "crpl.h"
#endif

  /** Read in crystal plasticity information. */
  void read_cryst_plast (FILE *in1,
             long nmat,
             CRPL *crpl,
             const int plc);

#endif /* #ifndef READ_CRYST_PLAST_H */
