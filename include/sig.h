/* HEADER */
#ifndef SIG_H
#define SIG_H

#include "element.h"

/** Structure of stresses SIG */
struct IL0_sig { /* Inelastic stress in all integration points */
  double *o,
    *f,
    *m,
    *d,
    *a;
  /* crystal plasticity */
  double *Tau,Har,Har1;
};

/** Structure of stresses SIG */
struct IL1_sig { /* Inelastic stress in all integration points */
  double *o,*f,*m,*d,*a;
};

struct SIG {
  /** Elastic stress */
  struct {
    /** Overall stress. */
    double *o;
    /** Stress in the fibre */
    double *f;
    /** Stress in the matrix */
    double *m;
    double *d,eq,eq_m;
  } el;

  IL0_sig *il;
  IL1_sig *d_il;
  /** crystal plasticity */
  double *p,*d_p,*pn_1;
};

SIG* build_sig_el(const long ne);

void destroy_sig_el(SIG* sig,
                    const long ne);

SIG* build_sig_il (const long ne,
                   const int analysis,
                   Element *elem);

/*** MUST be called before destroy_elem */
void destroy_sig_il(SIG* sig,
                    const Element *elem,
                    const long ne,
                    const int analysis);

#endif /* #ifndef  */
