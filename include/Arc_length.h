#ifndef ARC_LENGTH_H
#define ARC_LENGTH_H

#include "BSprivate.h"

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

#ifndef HYPRE_GLOBAL_H
#include "hypre_global.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef SUPP_H
#include "supp.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef MATGEOM_H
#include "matgeom.h"
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

#ifndef PGFEM_COMM_H
#include "pgfem_comm.h"
#endif

#ifndef BOUNDING_ELEMENT_H
#include "bounding_element.h"
#endif

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  double Arc_length (long ne,
		     int n_be,
		     long nn,
		     long ndofn,
		     long ndofd,
		     long npres,
		     long nt,
		     long tim,
		     double *times,
		     double nor_min,
		     long iter_max,
		     double dt,
		     double dt0,
		     ELEMENT *elem,
		     BOUNDING_ELEMENT *b_elems,
		     long nbndel,
		     long *bndel,
		     NODE *node,
		     SUPP sup,
		     double *sup_defl,
		     HOMMAT *hommat,
		     MATGEOM matgeom,
		     SIG *sig_e,
		     EPS *eps,
		     int *Ap,
		     int *Ai,
		     BSspmat *k,
		     PGFEM_HYPRE_solve_info *PGFEM_hypre,
		     double *RRn,
		     double *f_defl,
		     CRPL *crpl,
		     double stab,
		     long nce,
		     COEL *coel,
		     double *r,
		     double *f,
		     double *d_r,
		     double *D_R,
		     double *rr,
		     double *R,
		     double *RR,
		     double *f_u,
		     double *U,
		     double *DK,
		     double *dR,
		     double *BS_f,
		     double *BS_d_r,
		     double *BS_D_R,
		     double *BS_rr,
		     double *BS_R,
		     double *BS_RR,
		     double *BS_f_u,
		     double *BS_U,
		     double *BS_DK,
		     double *BS_dR,
		     long FNR,
		     double lm,
		     double dAL0,
		     double *DET0,
		     double *DLM0,
		     double *dlmdlm,
		     long gr2,
		     long gr4,
		     SIG *sig_n,
		     char *out_dat,
		     long *print,
		     long *AT,
		     long ARC,
		     double dALMAX,
		     long *ITT,
		     double *DAL,
		     double *pores,
		     /*long nge,
		       GEEL *geel,
		       long ngn,
		       GNOD *gnod*/
		     long *DomDof,
		     int GDof,
		     COMMUN comm,
		     BSprocinfo *BSinfo,
		     BSpar_mat **pk,
		     BSpar_mat **f_pk,
		     double err,
		     double *NORM,
		     MPI_Comm mpi_comm,
		     const double VVolume,
		     const PGFem3D_opt *opts);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef ARC_LENGTH_H */
