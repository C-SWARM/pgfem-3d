#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

#include "BSprivate.h"

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
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

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifndef COMPUTE_FORCE_ON_MODEL_ENT
#include "compute_force_on_model_ent.h"
#endif

#ifndef PGFEM_COMM_H
#include "pgfem_comm.h"
#endif

#ifndef HYPRE_GLOBAL_H
#include "hypre_global.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  double Newton_Raphson (const int print_level,
			 long ne,
			 int n_be,
			 long nn,
			 long ndofn,
			 long ndofd,
			 long npres,
			 long tim,
			 double *times,
			 double nor_min,
			 double dt,
			 ELEMENT *elem,
			 BOUNDING_ELEMENT *b_elems,
			 NODE *node,
			 SUPP sup,
			 double *sup_defl,
			 HOMMAT *hommat,
			 MATGEOM matgeom,
			 SIG *sig_e,
			 EPS *eps,
			 int *Ap,
			 int *Ai,
			 double *r,
			 double *f,
			 double *d_r,
			 double *rr,
			 double *R,
			 double *f_defl,
			 double *RR,
			 double *f_u,
			 double *RRn,
			 CRPL *crpl,
			 double stab,
			 long nce,
			 COEL *coel,
			 long FNR,
			 double *pores,
			 PGFEM_HYPRE_solve_info *PGFEM_hypre,
			 BSprocinfo *BSinfo,
			 BSspmat *k,
			 BSpar_mat **pk,
			 BSpar_mat **f_pk,
			 double *BS_x,
			 double *BS_f,
			 double *BS_RR,
			 double gama,
			 double GNOR,
			 double nor1,
			 double err,
			 double *BS_f_u,
			 long *DomDof,
			 COMMUN comm,
			 int GDof,
			 long nt,
			 long iter_max,
			 double *NORM,
			 long nbndel,
			 long *bndel,
			 MPI_Comm mpi_comm,
			 const double VVolume,
			 const PGFem3D_opt *opts,
			 MODEL_ENTITY *me,
			 double *forces,
			 void *ms_job_list,
			 void *microscale);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef NEWTON_RAPHSON_H */
