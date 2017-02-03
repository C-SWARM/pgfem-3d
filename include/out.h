/* HEADER */
#pragma once
#ifndef OUT_H
#define OUT_H

#include "PGFEM_mpi.h"
#include "PGFEM_io.h"
#include "element.h"
#include "node.h"
#include "supp.h"
#include "cohesive_element.h"
#include "eps.h"
#include "sig.h"
#include "ensight.h"
#include "PGFem3D_options.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

void logo (FILE *out);

void coordinates (FILE *out,
		  NODE *node,
		  long nn);

void deform (FILE *out,
	     NODE *node,
	     ELEMENT *elem,
	     long nn,
	     long ne,
	     long ndofn,
	     SUPP sup,
	     double *r);

void stress_out (FILE *out,
		 long ne,
		 long nn,
		 ELEMENT *elem,
		 SIG *sig_e,
		 SIG *sig_n,
		 long gr4);

void strain_out (FILE *out,
		 long ne,
		 ELEMENT *elem,
		 EPS *eps,
		 const PGFem3D_opt *opts);

void deform_grad_out (FILE *out,
		      long ne,
		      ELEMENT *elem,
		      EPS *eps);

void macro_fields_out (FILE *out,
		       EPS *eps,
		       const PGFem3D_opt *opts);

void cohesive_out (FILE *out,
		   long nce,
		   COEL *coel);

void damage_out(FILE *out,
		const long ne,
		const ELEMENT *elem,
		const EPS *eps);


void elixir (char jmeno[50],
	     long nn,
	     long ne,
	     long ndofn,
	     NODE *node,
	     ELEMENT *elem,
	     SUPP sup,
	     double *r,
	     SIG *sig_e,
	     SIG *sig_n,
	     EPS *eps,
	     long gr4,
	     long nce,
	     COEL *coel,
	     const PGFem3D_opt *opts);

void EnSight (char jmeno[500],
	      long tim,
	      long nt,
	      long nn,
	      long ne,
	      long ndofn,
	      NODE *node,
	      ELEMENT *elem,
	      SUPP sup,
	      double *r,
	      SIG *sig_e,
	      SIG *sig_n,
	      EPS *eps,
	      long gr4,
	      long nce,
	      COEL *coel,
	      /*long nge,
		GEEL *geel,
		long ngn,
		GNOD *gnod,
	      */long FNR,
	      double lm,
	      ENSIGHT ensight,
	      MPI_Comm mpi_comm,
	      const PGFem3D_opt *opts);

void ASCII_output(const PGFem3D_opt *opts,
		  MPI_Comm comm,
		  long tim,
		  double *times,
		  long  Gnn,
		  long nn,
		  long ne,
		  long nce,
		  long ndofd,
		  long *DomDof,
		  int *Ap,
		  long FNR,
		  double lm,
		  double pores,
		  double VVolume,
		  NODE *node,
		  ELEMENT *elem,
		  SUPP sup,
		  double *r,
		  EPS *eps,
		  SIG *sig_e,
		  SIG *sig_n,
		  COEL *coel);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
