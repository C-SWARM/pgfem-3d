/* HEADER */
#ifndef SIG_H
#define SIG_H

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure of stresses SIG */
  typedef struct { /* Inelastic stress in all integration points */
    double *o,
      *f,
      *m,
      *d,
      *a;
    /* crystal plasticity */
    double *Tau,Har,Har1;
  }IL0_sig;

  /** Structure of stresses SIG */
  typedef struct { /* Inelastic stress in all integration points */
    double *o,*f,*m,*d,*a;
  }IL1_sig;

  struct SIG{
    /** Elastic stress */
    struct {
      /** Overall stress. */
      double *o;
      /** Stress in the fibre */
      double *f;
      /** Stress in the matrix */
      double *m;
      double *d,eq,eq_m;
    }el;
    IL0_sig *il;
    IL1_sig *d_il;
    /** crystal plasticity */
    double *p,*d_p,*pn_1;
  };
  typedef struct SIG SIG;

  SIG* build_sig_el(const long ne);

  void destroy_sig_el(SIG* sig,
		      const long ne);

  SIG* build_sig_il (const long ne,
		     const int analysis,
		     ELEMENT *elem);

  /*** MUST be called before destroy_elem */
  void destroy_sig_il(SIG* sig,
		      const ELEMENT *elem,
		      const long ne,
		      const int analysis);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
