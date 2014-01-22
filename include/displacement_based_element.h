/* HEADER */

/** This file declares the functions for the displacement-based TOTAL
    LAGRANGIAN element. Note that there are several helper functions
    which are defined in the definitions file which are not listed
    here.  This is on purpose as these functions should ONLY be used
    within the diplacement based element context. */

#ifndef DISP_BASED_ELEM_H
#define DISP_BASED_ELEM_H

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef BOUNDING_ELEMENT_H
#include "bounding_element.h"
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

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** get the material potential */
  int DISP_get_material_potential(const double kappa,
				  const HOMMAT *hommat, /* pointer to single hommat */
				  const double *C,
				  const double J,
				  double *Wbar);

  /** Compute the element stiffness for a volumetric element. */
  int DISP_stiffmat_el(double *Ks,
		       const int elem_id,
		       const int ndofn,
		       const int nne,
		       const double *x,
		       const double *y,
		       const double *z,
		       const ELEMENT *elem,
		       const HOMMAT *hommat,
		       const long *nod,
		       const NODE *node,
		       const EPS *eps,
		       const SIG *sig,
		       const SUPP sup,
		       const double *disp);

  /** Compute the element residual for a volumetric element */
  int DISP_resid_el(double *R,
		    const int elem_id,
		    const int ndofn,
		    const int nne,
		    const double *x,
		    const double *y,
		    const double *z,
		    const ELEMENT *elem,
		    const HOMMAT *hommat,
		    const long *nod,
		    const NODE *node,
		    const EPS *eps,
		    const SIG *sig,
		    const SUPP sup,
		    const double *disp);

  /** Compute the residual on a boundary element (Lagrange multiplier
      formulation) */
  int DISP_resid_bnd_el(double *R,
			const int b_el_id,
			const int ndofn,
			const int ndof_ve,
			const double *x_ve,
			const double *y_ve,
			const double *z_ve,
			const BOUNDING_ELEMENT *b_elems,
			const ELEMENT *elem,
			const HOMMAT *hommat,
			const NODE *node,
			const EPS *eps,
			const SIG *sig,
			const SUPP sup,
			const double *vol_disp);

  /** Compute the tangent on a boundary element (Lagrange multiplier
      formulation). */
  int DISP_stiffmat_bnd_el(double *Ks,
			   const int b_el_id,
			   const int ndofn,
			   const int ndof_ve,
			   const double *x_ve,
			   const double *y_ve,
			   const double *z_ve,
			   const BOUNDING_ELEMENT *b_elems,
			   const ELEMENT *elem,
			   const HOMMAT *hommat,
			   const NODE *node,
			   const EPS *eps,
			   const SIG *sig,
			   const SUPP sup,
			   const double *vol_disp);

  /** Compute the stress and strain for visualization. */
  void DISP_increment_el(const ELEMENT *elem,
			 const int elem_id,
			 const int nne,
			 const NODE *node,
			 const long *nod,
			 const int ndofn,
			 const double *x,
			 const double *y,
			 const double *z,
			 EPS *eps,
			 SIG *sig,
			 const SUPP sup,
			 const HOMMAT *hommat,
			 const double *disp);

  /** Compute quantities for visualization volume in current
      configuration */
  void DISP_increment(const ELEMENT *elem,
		      const int nelem,
		      NODE *node,
		      const int nnodes,
		      const int ndofn,
		      SUPP sup,
		      EPS *eps,
		      SIG *sig,
		      const HOMMAT *hommat,
		      const double *sol_incr,
		      const double *sol,
		      MPI_Comm mpi_comm);

  /** Compute element contributions to macro and mixed tangents in
      multiscale modeling of coheisive interfaces */
  int DISP_cohe_micro_terms_el(double *K_00_e,
			       double *K_01_e,
			       double *K_10_e,
			       double *traction_res_e,
			       double *traction_e,
			       /* macro information */
			       const int macro_nnode,
			       const int macro_ndofn,
			       const double macro_int_wt,
			       const double *macro_shape_func,
			       const double *macro_normal,
			       const double layer_thickness,
			       /* micro information */
			       const double micro_volume0,
			       const int elem_id,
			       const int ndofn,
			       const int nne,
			       const double *x,
			       const double *y,
			       const double *z,
			       const ELEMENT *elem,
			       const HOMMAT *hommat,
			       const long *nod,
			       const NODE *node,
			       const EPS *eps,
			       const SIG *sig,
			       const SUPP sup,
			       const double *disp);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef DISP_BASED_ELEM_H */
