/* HEADER */
#ifndef APPLIED_TRACTION_H
#define APPLIED_TRACTION_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef NODE_H
#include "node.h"
#endif

#ifndef SIG_H
#include "sig.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

  typedef struct SURFACE_TRACTION_ELEM{
    int elem_id; /**< index to get elem from array */
    int n_faces; /**< number of faces on the element assoc. w/ feat*/
    int *faces; /**< face ids */
    int *feat_num; /** feat id from original list */
  } SUR_TRAC_ELEM;
 
  /** Reads the list of features and loads to apply from a
      file. Returns allocated arrays and error flag (SUCCESS = 0). */
  int read_applied_surface_tractions_fname(char *fname,
					   int *n_feats,
					   int **feat_type,
					   int **feat_id,
					   double **loads);

  int read_applied_surface_tractions(FILE *in,
				     int *n_feats,
				     int **feat_type,
				     int **feat_id,
				     double **loads);

  /** Get the list of elements and their faces which have applied
      tractions. Done once for the entire simulation. feat_type,
      feat_id and load_id are all of length n_feats */
  int generate_applied_surface_traction_list(const int ne,
					     const ELEMENT *elem,
					     const int n_feats,
					     const int *feat_type,
					     const int *feat_id,
					     int *n_sur_trac_elem,
					     SUR_TRAC_ELEM **sur_trac_elem);

  int destroy_applied_surface_traction_list(const int n_sur_trac_elem,
					    SUR_TRAC_ELEM *sur_trac_elem);

  int compute_applied_traction_res(const int ndofn,
				   const NODE *nodes,
				   const ELEMENT *elem,
				   const int n_ste,
				   const SUR_TRAC_ELEM *ste,
				   const int n_feats,
				   const double *loads,
				   double *res,
				   const int mp_id);

  /** integrate the force on the marked boundaries of the LOCAL DOMAIN
      in the LAGRANGIAN FRAME. forces vector is [n_feats x ndim]
      non-collective */
  int compute_resultant_force(const int n_feats,
			      const int n_ste,
			      const SUR_TRAC_ELEM *ste,
			      const NODE *nodes,
			      const ELEMENT *elem,
			      const SIG *sig,
			      const EPS *eps,
			      double *forces);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
