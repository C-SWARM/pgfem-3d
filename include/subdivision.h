/* HEADER */
#pragma once
#ifndef SUBDIVISION_H
#define SUBDIVISION_H

#include "PGFEM_mpi.h"
#include "element.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  void subdivision (long INFO,
		    double *dt,
		    long *STEP,
		    long *DIV,
		    long tim,
		    double *times,
		    long *ST,
		    long ne,
		    long ndofn,
		    long ndofd,
		    long npres,
		    ELEMENT *elem,
		    CRPL *crpl,
		    EPS *eps,
		    SIG *sig,
		    SUPP sup,
		    double *sup_defl,
		    double *rr,
		    double *d_r,
		    double *f_defl,
		    double *f,
		    double *RRn,
		    double *R,
		    long *GAMA,
		    double *DT,
		    long *OME,
		    double stab,
		    int iter,
		    int max_iter,
		    double alpha,
		    MPI_Comm mpi_comm,
		    const int analysis);

  double subdiv_arc (long INFO,
		     double *dt,
		     double dt0,
		     long *STEP,
		     long *DIV,
		     long tim,
		     double *times,
		     long *ST,
		     long ne,
		     long ndofd,
		     long npres,
		     ELEMENT *elem,
		     CRPL *crpl,
		     EPS *eps,
		     SIG *sig,
		     SUPP sup,
		     double *sup_defl,
		     double *rr,
		     double *d_r,
		     double *D_R,
		     double *f_defl,
		     double *f,
		     long *GAMA,
		     double *DT,
		     long *OME,
		     double stab,
		     double dAL0,
		     double dAL,
		     double dALMAX,
		     double nor_min,
		     double dlm0,
		     long *ITT,
		     long iter,
		     long iter_max,
		     long TYPE,
		     MPI_Comm mpi_comm,
		     const int analysis);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef SUBDIVISION_H */
