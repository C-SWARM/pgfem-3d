#ifndef RES_FINI_DEF_H
#define RES_FINI_DEF_H

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
  void res_fini_def (long ne,
		     long npres,
		     ELEMENT *elem,
		     EPS *eps,
		     SIG *sig,
		     CRPL *crpl,
		     const int analysis);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef RES_FINI_DEF_H */
