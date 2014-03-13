/* HEADER */
#pragma once
#ifndef COHESIVE_ELEMENT_H
#define COHESIVE_ELEMENT_H

/**
 * @file This file describes the cohesive element and generic
 * functions regarding its use.
 */

#include "PGFEM_io.h"
#include "node.h"
#include "cohesive_potentials.h"
#include "ensight.h"
#include "eps.h"
#include "supp.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure of COHESIVE ELEMENTS */
  struct COEL {
    long toe,*nod,pr,typ;
    double Sc,Xc,b,Jjn,*e1,*e2,*n,*x,*y,
      *z,*Xmax,*tmax,
      *Xi,*ti,Xxi,txi,tn,ts,Xn,Xs,vo,k;
    const cohesive_props *props; /* this is a pointer to a particular
				    set of material properties. */

    /* internal cohesive state variables */
    int nvar;
    double **vars; /* variables at each ip */
  };
  typedef struct COEL COEL;

  void destroy_coel(COEL* coel, long nce);

  COEL* read_cohe_elem (FILE* in1,
			long ncom,
			long ndofn,
			long nn,
			NODE *node,
			long *NCE,
			double **comat,
			ENSIGHT ensight,
			long gr2,
			int myrank,
			const cohesive_props *co_props);

  void stiff_mat_coh (long ii,
		      long ndofn,
		      long nne,
		      long *nod,
		      double *x,
		      double *y,
		      double *z,
		      COEL *coel,
		      double *r_e,
		      double *Kch,
		      double nor_min,
		      EPS *eps,
		      long FNR,
		      double lm,
		      double *fe,
		      int myrank);

  void resid_co_elem (long ii,
		      long ndofn,
		      long nne,
		      long *nod,
		      double *x,
		      double *y,
		      double *z,
		      COEL *coel,
		      double *r_e,
		      double *fe,
		      double nor_min,
		      int myrank);

  int increment_cohesive_elements(const int nce,
				  COEL *coel,
				  double *pores,
				  const NODE *node,
				  const SUPP sup,
				  const double *d_r);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
