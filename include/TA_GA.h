#ifndef TA_GA_H
#define TA_GA_H

#include "data_structure.h"

#ifndef CRPL_H
#include "crpl.h"
#endif

  /** */
  void TA_GA (long nss,
          long mat,
          CRPL *crpl,
          double Har,
          double **eFn,
          double **S,
          double *TA,
          double *GA);

#endif /* #ifndef TA_GA_H */
