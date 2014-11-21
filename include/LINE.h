#ifndef LINE_H
#define LINE_H

#include "PGFEM_mpi.h"

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef MATGEOM_H
#include "matgeom.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef SUPP_H
#include "supp.h"
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

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef BOUNDING_ELEMENT_H
#include "bounding_element.h"
#endif

#ifndef PGFEM_COMM_H
#include "pgfem_comm.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  long LINE_S1 (double *nor,
		double *gama,
		double nor1,
		double NOR,
		long iter,
		double *f_u,
		long ne,
		int n_be,
		long ndofd,
		long ndofn,
		long npres,
		double *d_r,
		double *r,
		NODE *node,
		ELEMENT *elem,
		BOUNDING_ELEMENT *b_elems,
		MATGEOM matgeom,
		HOMMAT *hommat,
		SUPP sup,
		EPS *eps,
		SIG *sig_e,
		double nor_min,
		CRPL *crpl,
		double dt,
		double stab,
		long nce,
		COEL *coel,
		double *f,
		double *rr,
		double *RR,
		long tim,
		/*GNOD *gnod,
		  GEEL *geel,
		*/double *BS_f,
		double *BS_RR,
		double *BS_f_u,
		long *DomDof,
		long STEP,
		COMMUN comm,
		int GDof,
		MPI_Comm mpi_comm,
		double *max_damage,
		double *dissipation,
		const PGFem3D_opt *opts);

  long LINE_S3 (double *nor,
		double *nor2,
		double *gama,
		double nor1,
		double NOR,
		double LS1,
		long iter,
		double *f_u,
		long ne,
		int n_be,
		long ndofd,
		long ndofn,
		long npres,
		double *d_r,
		double *r,
		NODE *node,
		ELEMENT *elem,
		BOUNDING_ELEMENT *b_elems,
		MATGEOM matgeom,
		HOMMAT *hommat,
		SUPP sup,
		EPS *eps,
		SIG *sig_e,
		double nor_min,
		CRPL *crpl,
		double dt,
		double t,
		double stab,
		long nce,
		COEL *coel,
		double *f,
		double *rr,
		double *RR,
		long tim,
		double *BS_f,
		double *BS_RR,
		double *BS_f_u,
		long *DomDof,
		COMMUN comm,
		int GDof,
		long STEP,
		MPI_Comm mpi_comm,
		double *max_damage,
		double *dissipation,
		const PGFem3D_opt *opts,double alpha, double *r_n, double *r_n_1);

  long ALINE_S3 (long ARC,
		 double *DLM,
		 double *nor,
		 double *nor2,
		 double *gama,
		 double nor1,
		 double LS1,
		 long iter,
		 long ne,
		 int n_be,
		 long ndofd,
		 long ndofn,
		 long npres,
		 long tim,
		 double nor_min,
		 double dt,
		 double stab,
		 long nce,
		 double dlm,
		 double lm,
		 double dAL,
		 double *d_r,
		 double *r,
		 double *D_R,
		 NODE *node,
		 ELEMENT *elem,
		 BOUNDING_ELEMENT *b_elems,
		 MATGEOM matgeom,
		 HOMMAT *hommat,
		 SUPP sup,
		 EPS *eps,
		 SIG *sig_e,
		 CRPL *crpl,
		 COEL *coel,
		 double *f_u,
		 double *f,
		 double *R/*,
			    GNOD *gnod,
			    GEEL *geel*/,
 		 double *BS_f,
		 double *BS_R,
		 double *BS_D_R,
		 double *BS_d_r,
		 double *BS_DK,
		 double *BS_U,
		 double *BS_rr,
		 long *DomDof,
 		 int GDof,
		 COMMUN comm,
		 long STEP,
		 MPI_Comm mpi_comm ,
		 double *max_damage,
		 double *dissipation,
		 const PGFem3D_opt *opts);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef LINE_H */
