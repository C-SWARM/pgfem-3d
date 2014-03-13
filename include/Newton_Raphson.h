#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

#include "PGFEM_mpi.h"
#include "element.h"
#include "matgeom.h"
#include "hommat.h"
#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "crpl.h"
#include "cohesive_element.h"
#include "bounding_element.h"
#include "PGFem3D_options.h"
#include "compute_force_on_model_ent.h"
#include "pgfem_comm.h"

#include "hypre_global.h"
#include "blocksolve_interface.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**
   * Newton-Raphson solutionn algorithm.
   */
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
