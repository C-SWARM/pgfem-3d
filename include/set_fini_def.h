#ifndef SET_FINI_DEF
#define SET_FINI_DEF

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifndef CRPL_H
#include "crpl.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

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

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef SET_FINI_DEF */
