/* HEADER */
#pragma once
#ifndef ELEMENT_H
#define ELEMENT_H

#include "PGFEM_io.h"
#include "supp.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure of element properties */
  struct ELEMENT{
  
    long pr, /**< Property of element */

      toe, /**< Type of element. This is synonomous to the number of
	      nodes on the element. */

      *nod, /**< Which nodes are on the element. The nodes are
	       identified by their local numbers. */

      *mat, /**< Pointer to material. mat[0] = matrix | mat[1] =
	       fibre | mat[2] = homogeneous medium */

      *hom; /**< Pointer to volume fraction (cf), fibre
	       orientation (psi). hom[0] = cf | hom[1] = psi */

    int *bnd_type; /**< t3d feature type on face */
    int *bnd_id; /**< t3d feature id on face */

    /** DOFs which cannot be condensed out of the global system. Note
	that the values associated with these DOFs will be treated in
	the same way as nodeal DOFs and are not stored with the element
	as with bubble dofs. */
    int n_dofs;
    long *G_dof_ids;
    long *L_dof_ids;

    /** Bounding elements */
    int n_be;
    long *be_ids;

    int n_bub;            /**< number of bubble nodes in element */
    int n_bub_dofs;       /**< number of dofs per bubble node */
    double *bub_dofs;     /**< Array of bubble dof values */
    double *d_bub_dofs;   /**< Bubble dof increment */

    double **L;
    long *LO;
  };
  typedef struct ELEMENT ELEMENT;

  ELEMENT* build_elem(FILE *in,
		      const long ne,
		      const int analysis);

  void destroy_elem(ELEMENT *elem,
		    const long ne);

  /* Function reads parameters of elements */
  void read_elem (FILE *in,
		  long ne,
		  ELEMENT *elem,
		  SUPP sup,
		  const int legacy);

  /** */
  void write_element_fname(const char *filename,
			   const int nelem,
			   const ELEMENT *elems);

  void write_element(FILE *ofile,
		     const int nelem,
		     const ELEMENT *elems);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef ELEMENT_H */
