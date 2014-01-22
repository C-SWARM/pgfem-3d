#ifndef PLC_H
#define PLC_H

#ifndef CRPL_H
#include "crpl.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  double PLC_diff_g_ga (long ii,
			long ip,
			long nss,
			long mat,
			double dt,
			CRPL *crpl,
			double *GA,
			EPS *eps);

  /** */
  long PLC_Inaccessible (long ii,
			 long ip,
			 double nor_min,
			 long nss,
			 long mat,
			 double dt,
			 CRPL *crpl,
			 double *GA,
			 double *GA1,
			 EPS *eps,
			 double **UU,
			 double HAR);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef PLC_H */
