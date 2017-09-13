/**
 * @file This file describes the cohesive element and generic functions
 * regarding its use.
 */
#ifndef PGFEM3D_COHESIVE_ELEMENT_H
#define PGFEM3D_COHESIVE_ELEMENT_H

#include "PGFEM_io.h"
#include "cohesive_potentials.h"
#include "data_structure.h"
#include "ensight.h"
#include "eps.h"
#include "node.h"
#include "supp.h"

/** Structure of COHESIVE ELEMENTS */
struct COEL {
  long toe;  /**< Number of nodes on the element */
  long *nod; /**< Node ids on the element */
  long pr;   /**< property flag */
  long typ;  /**< type flag */
  long mat; /**< material id */

  /** depricated material properties (should not be used) */
  double Sc,Xc,b,k;

  /** Transformation for updated Lagrangian formulation. Note that
      cohesive elements are *always* total Lagrangian, thus Jjn = 1
      always. */
  double Jjn;

  /** buffers for basis vectors (updated every time) */
  double *e1,*e2,*n;

  /** buffers for mean-map coordinates (updated every time) */
  double *x,*y,*z;

  /** depricated internal state variables */
  double *Xmax,*tmax;

  /** average element  values for visualization */
  double *Xi,*ti,Xxi,txi,tn,ts,Xn,Xs;

  /** unknown, unused? */
  double vo;

  /** pointer to a particular set of material properties. */
  const cohesive_props *props;

  /** internal cohesive state variables */
  int nvar;
  double **vars; /* variables at each ip */
};

void destroy_coel(COEL* coel, long nce);

/**
 * Reset the coel properties accorting to p_coel->mat.
 */
void reset_coel_props(const cohesive_props *co_props,
                      COEL *p_coel);

COEL* read_cohe_elem (FILE* in1,
                      long ncom,
                      long ndofn,
                      long nn,
                      Node *node,
                      long *NCE,
                      double **comat,
                      Ensight *ensight,
                      long gr2,
                      int myrank,
                      const cohesive_props *co_props);

void stiff_mat_coh (long ii,
                    long ndofn,
                    long nne,
                    long *nod,
                    double *x,
                    double *y,
                    double *z,
                    COEL *coel,
                    double *r_e,
                    double *Kch,
                    double nor_min,
                    EPS *eps,
                    long FNR,
                    double lm,
                    double *fe,
                    int myrank);

void resid_co_elem (long ii,
                    long ndofn,
                    long nne,
                    long *nod,
                    double *x,
                    double *y,
                    double *z,
                    COEL *coel,
                    double *r_e,
                    double *fe,
                    double nor_min,
                    int myrank);

int increment_cohesive_elements(const int nce,
                                COEL *coel,
                                double *pores,
                                const Node *node,
                                const SUPP sup,
                                const double *d_r,
                                const int mp_id);

/**
 * Get the storage size (in bytes) of the internal state variables
 * for all cohesive elements in the domain.
 */
size_t coel_list_get_state_length_bytes(const int nce,
                                        const COEL *coel);

/**
 * Pack the cohesive element state variables into a buffer.
 *
 * BUFFER is the buffer to copy data to. BUF_POS is the current
 * insertion position in BUFFER. \return modified BUFFER and updated
 * insertion location BUF_POS.
 */
void coel_list_pack_state(const int nce,
                          const COEL *coel,
                          char *buffer,
                          size_t *buf_pos);

/**
 * Unpack the cohesive element state variables from a buffer.
 *
 * BUFFER is the buffer to copy data from. BUF_POS is the current
 * copy position in BUFFER. \return modified cohesive elements COEL
 * and updated insertion/copy location BUF_POS.
 */
void coel_list_unpack_state(const int nce,
                            COEL *coel,
                            const cohesive_props *co_props,
                            const char *buffer,
                            size_t *buf_pos);

#endif // #define PGFEM3D_COHESIVE_ELEMENT_H
