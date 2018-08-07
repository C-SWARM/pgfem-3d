#ifndef PGFEM3D_APPLIED_TRACTION_H
#define PGFEM3D_APPLIED_TRACTION_H

#include "PGFEM_io.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "node.h"
#include "sig.h"

struct SURFACE_TRACTION_ELEM {
  int elem_id; /**< index to get elem from array */
  int n_faces; /**< number of faces on the element assoc. w/ feat*/
  int *faces; /**< face ids */
  int *feat_num; /** feat id from original list */
};

/** Reads the list of features and loads to apply from a
    file. Returns allocated arrays and error flag (SUCCESS = 0). */
int read_applied_surface_tractions_fname(char *fname,
                                         int *n_feats,
                                         int **feat_type,
                                         int **feat_id,
                                         double **loads,
					 int myrank);

int read_applied_surface_tractions(FILE *in,
                                   int *n_feats,
                                   int **feat_type,
                                   int **feat_id,
                                   double **loads);

/** Get the list of elements and their faces which have applied
    tractions. Done once for the entire simulation. feat_type,
    feat_id and load_id are all of length n_feats */
int generate_applied_surface_traction_list(const int ne,
                                           const Element *elem,
                                           const int n_feats,
                                           const int *feat_type,
                                           const int *feat_id,
                                           int *n_sur_trac_elem,
                                           SURFACE_TRACTION_ELEM **sur_trac_elem);

int destroy_applied_surface_traction_list(const int n_sur_trac_elem,
                                          SURFACE_TRACTION_ELEM *sur_trac_elem);

int compute_applied_traction_res(const int ndofn,
                                 const Node *nodes,
                                 const Element *elem,
                                 const int n_ste,
                                 const SURFACE_TRACTION_ELEM *ste,
                                 const int n_feats,
                                 const double *loads,
                                 const EPS *eps,
                                 double *res,
                                 const int mp_id);

/** integrate the force on the marked boundaries of the LOCAL DOMAIN
    in the LAGRANGIAN FRAME. forces vector is [n_feats x ndim]
    non-collective */
int compute_resultant_force(const int n_feats,
                            const int n_ste,
                            const SURFACE_TRACTION_ELEM *ste,
                            const Node *nodes,
                            const Element *elem,
                            const SIG *sig,
                            const EPS *eps,
                            double *forces);

#endif /* #define PGFEM3D_APPLIED_TRACTION_H  */
