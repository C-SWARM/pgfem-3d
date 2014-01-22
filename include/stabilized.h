/****************************
 *  Program FEM3d ver. 3.0  *
 *  FEM - 3D analysis       *
 *  Karel Matous            *
 ****************************/

#ifndef STABILIZED_H
#define STABILIZED_H

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef COHESIVE_ELEMENT_H
#include "cohesive_element.h"
#endif

#ifndef MATGEOM_H
#include "matgeom.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

#ifndef SUPP_H
#include "supp.h"
#endif

/** Compute the material potential at an integration point */
int stab_get_material_potential(double *Wbar,
				const double kappa,
				const double Un_1, /* U(n-1) */
				const double Jn_1, /* J(n-1) */
				const double Jnp1, /* J)n+1) */
				const double Pn_1, /* p(n-1) */
				const double Pn,
				const double P,    /* p(n+1) */
				const double *C,
				const HOMMAT *mat);

void res_stab_def (long ne,
		   long npres,
		   ELEMENT *elem,
		   EPS *eps,
		   SIG *sig,
		   double stab);

int resid_st_elem (long ii,
		   long ndofn,
		   long nne,
		   ELEMENT *elem,
		   long *nod,
		   NODE *node,
		   HOMMAT *hommat,
		   double *x,
		   double *y,
		   double *z,
		   EPS *eps,
		   SIG *sig,
		   double *r_e,
		   double nor_min,
		   double *fe,
		   double dt,
		   double stab);

int stiffmatel_st (long ii,
		   long ndofn,
		   long nne,
		   double *x,
		   double *y,
		   double *z,
		   ELEMENT *elem,
		   HOMMAT *hommat,
		   long *nod,
		   NODE *node,
		   SIG *sig,
		   EPS *eps,
		   double *r_e,
		   long npres,
		   double nor_min,
		   double *Ks,
		   double dt,
		   double stab,
		   long FNR,
		   double lm,
		   double *fe);

int st_increment (long ne,
		  long nn,
		  long ndofn,
		  long ndofd,
		  MATGEOM matgeom,
		  HOMMAT *hommat,
		  ELEMENT *elem,
		  NODE *node,
		  SUPP sup,
		  EPS *eps,
		  SIG *sig,
		  double *d_r,
		  double *r,
		  double nor_min,
		  double stab,
		  double dt,
		  long nce,
		  COEL *coel,
		  double *pores,
		  MPI_Comm mpi_comm,
		  const int coh);

#endif
