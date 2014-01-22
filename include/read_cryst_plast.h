#ifndef READ_CRYST_PLAST_H
#define READ_CRYST_PLAST_H

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef CRPL_H
#include "crpl.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Read in crystal plasticity information. */
  void read_cryst_plast (FILE *in1,
			 long nmat,
			 CRPL *crpl,
			 const int plc);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef READ_CRYST_PLAST_H */
