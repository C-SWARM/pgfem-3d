/*******************************************
 * Program FEM3d ver. 2.0                  *
 * FEM  - 3D analysis                      *
 * CTU, Department of Structural Mechanics *
 * Karel Matous & Jaroslav Kruis           *
 *******************************************/

/*****************/
/* November 2000 */
/*****************/

#ifndef OUT_H
#define OUT_H

#include "PGFEM_mpi.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef SUPP_H
#include "supp.h"
#endif

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef ENSIGHT_H
#include "ensight.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

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
#endif
