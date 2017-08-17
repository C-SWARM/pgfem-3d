/* HEADER */
#ifndef BOUNDING_ELEMENT_H
#define BOUNDING_ELEMENT_H


#include "PGFEM_mpi.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

  /** Define the bounding element type */
  typedef struct BOUNDING_ELEMENT{
    /* data */
    int id;
    int vol_elem_id;
    int nnodes;
    long *nodes;
    long *loc_nodes; /**< element local node id */
    /* double *normal; */

    int periodic;
    int master;
    int master_dom;
    int master_be_id;
    int slave_dom;
    int slave_be_id;
    double other_val;

    int n_dofs;
    long *G_dof_ids;
    long *L_dof_ids;
  } BOUNDING_ELEMENT;

  /** read the bounding element information from a file and allocate
      the data structure. */
  int read_bounding_elements(FILE *ifile,
                 const int ndof_be,
                 int *nbe,
                 BOUNDING_ELEMENT **b_elems,
                 MPI_Comm mpi_comm);

  /** read the bounding element information from a file with provided
      filename */
  int read_bounding_elements_fname(char *filename,
                   const int ndof_be,
                   int *nbe,
                   BOUNDING_ELEMENT **b_elems,
                   MPI_Comm mpi_comm);

  /** Write the bounding element information to a file. */
  int write_bounding_elements(FILE *ofile,
                  int nbe,
                  BOUNDING_ELEMENT *b_elems);

  /** Write the bounding element information to a file with a provided
      filename */
  int write_bounding_elements_fname(char *filename,
                    int nbe,
                    BOUNDING_ELEMENT *b_elems);

  /** allocate the bounding element data structure given the number of
      elements and a list of number of nodes per element */
  int construct_bounding_elements(int nbe,
                  int *nn_per_elem,
                  BOUNDING_ELEMENT **b_elems);

  /** Safely destroy the bounding element structure */
  int destroy_bounding_elements(int nbe,
                BOUNDING_ELEMENT *b_elems);

#endif /* #ifndef  */

