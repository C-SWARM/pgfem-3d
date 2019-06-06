/* HEADER */
#pragma once
#ifndef PGFEM3D_ELEMENT_H
#define PGFEM3D_ELEMENT_H

#include "PGFEM_io.h"
#include "bounding_element.h"
#include "node.h"
#include "supp.h"

#include <cassert>
#include <iostream>
#include <vector>

/** Structure of element properties */
struct Element {

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

  void init();
  void fini();

  /// Get the total number of degrees of freedom for this element.
  ///
  /// @param      ndofn The number of degrees of freedom per node.
  /// @param        bes The bounding element array.
  ///
  /// @returns          The total number of degrees of freedom for this
  ///                   element.
  size_t getNDof(size_t ndofn, const BoundingElement* bes) const
  {
    size_t n = toe * ndofn + n_dofs;
    for (size_t i = 0, e = n_be; i < e; ++i) {
      n += bes[be_ids[i]].n_dofs;
    }
    return n;
  }

  /// Apply the function operator to each ID in nodes as defined by the
  /// element's nod map.
  ///
  /// @tparam        Op The operator type to apply (called with a node id).
  ///
  /// @param        pid The multiphysics id that we're dealing with.
  /// @param      ndofn The number of degrees of freedom per node (we don't seem
  ///                   to store this anywhere).
  /// @param      nodes The array of nodes (we have the nod[] map to access
  ///                   them).
  /// @param         op The operator to apply.
  template <class Op>
  void forEachLinkedDof(int pid, size_t ndofn, const Node* nodes, Op&& op) const
  {
    for (size_t i = 0, e = toe; i < e; ++i) {
      for (size_t j = 0, e = ndofn; j < e; ++j) {
        auto lid = nodes[nod[i]].id_map[pid].id[j];
        auto gid = nodes[nod[i]].id_map[pid].Gid[j];
        op(lid, gid);
      }
    }
  }

  /// Apply the function operator to each DOF ID in the node.
  ///
  /// @tparam        Op The operator type to apply (called with a node id).
  ///
  /// @param         op The operator to apply.
  template <class Op>
  void forEachDof(Op&& op) const
  {
    for (size_t i = 0, e = n_dofs; i < e; ++i) {
      op(L_dof_ids[i], G_dof_ids[i]);
    }
  }

  /// Apply the function operator to each DOF ID in the bounding elements in
  /// this node (as mapped by be_ids).
  ///
  /// @tparam        Op The operator type to apply (called with a node id).
  ///
  /// @param        bes The global array of bounding elements.
  /// @param         op The operator to apply.
  template <class Op>
  void forEachBoundingDof(const BoundingElement* bes, Op&& op) const
  {
    for (size_t i = 0, e = n_be; i < e; ++i) {
      auto& be = bes[be_ids[i]];
      for (size_t j = 0, e = be.n_dofs; j < e; ++j) {
        op(be.L_dof_ids[j], be.G_dof_ids[j]);
      }
    }
  }

  /// Apply the function operator to each DOF ID in the mapped nodes, this
  /// element, and the mapped bounding elements.
  ///
  /// @tparam        Op The operator type to apply (called with a node id).
  ///
  /// @param        pid The multiphysics id that we're dealing with.
  /// @param      ndofn The number of degrees of freedom per node (we don't seem
  ///                   to store this anywhere).
  /// @param      nodes The array of nodes (we have the nod[] map to access
  ///                   them).
  /// @param        bes The global array of bounding elements.
  /// @param         op The operator to apply.
  template <class Op>
  void forEachDof(int pid, size_t ndofn, const Node* nodes,
                  const BoundingElement* bes, Op&& op) const
  {
    forEachLinkedDof(pid, ndofn, nodes, std::forward<Op>(op));
    forEachDof(std::forward<Op>(op));
    forEachBoundingDof(bes, std::forward<Op>(op));
  }
};

Element* build_elem(FILE *in,
                    const long ne,
                    const int analysis);

void destroy_elem(Element *elem,
                  const long ne);

/* Function reads parameters of elements */
void read_elem (FILE *in,
                long ne,
                Element *elem,
                SUPP sup,
                const int legacy);

/** */
void write_element_fname(const char *filename,
                         const int nelem,
                         const Element *elems);

void write_element(FILE *ofile,
                   const int nelem,
                   const Element *elems);

#endif /* #ifndef PGFEM3D_ELEMENT_H */
