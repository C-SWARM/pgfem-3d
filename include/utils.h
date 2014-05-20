/* HEADER */

/**
 * @file This file declares utility functions for PGFem3D.
 */

#pragma once
#ifndef UTILS_H
#define UTILS_H

#include "PGFEM_mpi.h"
#include "element.h"
#include "node.h"
#include "material.h"
#include "matgeom.h"
#include "hommat.h"
#include "bounding_element.h"
#include "pgfem_comm.h"
#include "sig.h"
#include "eps.h"

/** Dynamically allocate and populate a formated string */
int alloc_sprintf(char **str,
		  const char *format,
		  ...);

/** Pack data into contiguous array of char. pos is incremented to
    point to the next insertion point in buffer*/
void pack_data(const void *src,
	       char *buffer,
	       size_t *pos,
	       const size_t n_el,
	       const size_t size);

/** Unpack data from contiguous array of char. pos is incremented to
    point to the next extraction point in buffer */
void unpack_data(const char *buffer,
		 void *dest,
		 size_t *pos,
		 const size_t n_el,
		 const size_t size);

/** Copies the data from a general double-pointer matrix into a
    continuous arry row-major matrix of the dimension nrows x ncols */
void mat2array(double *array,
	       const double **mat,
	       const unsigned int nrows,
	       const unsigned int ncols);

/** Does the reverse mat2array. */
void array2mat(const double *array,
	       double **tensor,
	       const unsigned int I,
	       const unsigned int J);

/** Copies the data from a general triple-pointer tensor into a
    contiguous array tensor (I x J x K) with row-major sub-matrices 
    (J x K). */
void tensor3_2array(double *array,
		    const double ***tensor,
		    const unsigned int I,
		    const unsigned int J,
		    const unsigned int K);

/** Does the reverse of tensor3_2array. */
void array2tensor3(const double *array,
		   double ***tensor,
		   const unsigned int I,
		   const unsigned int J,
		   const unsigned int K);

/** Analogous to tensor3_2array but for a general quad-pointer
    tensor */
void tensor4_2array(double *array,
		    const double ****tensor,
		    const unsigned int I,
		    const unsigned int J,
		    const unsigned int K,
		    const unsigned int L);

/** Does the reverse of tensor4_2array. */
void array2tensor4(const double *array,
		   double ****tensor,
		   const unsigned int I,
		   const unsigned int J,
		   const unsigned int K,
		   const unsigned int L);

/** Copies the specially ordered Shape Tensor object into a contiguous
    array which may be accessed using the idx_4_gen function with the
    index order (node,dof,i,j). NOTE 1: (ij) is always 3 x 3. NOTE 2:
    the sub-matrices (ij) are contiguous in memory for convenient use
    with BLAS/LAPACK. */
void shapeTensor2array(double *array,
		       const double ****ST,
		       const unsigned int nne);

/** Copies the array form of the Shape Tensor object back to a
    quad-pointer format */
void array2shapeTensor(const double *array,
		       double ****ST,
		       const unsigned int nne);

/** Compute the determinant of a 3 x 3 row-major matrix. */
double det2x2(const double *mat);
double det3x3(const double *mat);
double det4x4(const double *mat);

double getJacobian(const double *mat,
		   const int elem_id,
		   int *err);

/** Compute the inverse of a 3 x 3 row-major matrix. */
int inv2x2(const double *mat,
	   double *mat_inv);
int inv3x3(const double *mat,
	   double *mat_inv);
int inv4x4(const double *mat,
	   double *mat_inv);

/** Compute the inverse of a N x N row-major matrix by factorization
    (LAPACK) and return 0 on success. */
int inverse(double const* A,
	    const int M,
	    double *A_I);

/** Copy the transpose of 'mat' into 'mat_t'. */
void transpose(double *mat_t,
	       const double *mat,
	       const int mat_row,
	       const int mat_col);

/** Compute sym = 1/2 (mat + mat'). NOTE: reuires square matrix. */
void symmetric_part(double *sym,
		    const double *mat,
		    const int dim);

/** Print the coordinates of an element to a file. */
void print_coords(FILE *out,
		  const int nne,
		  const double *x,
		  const double *y,
		  const double *z);

/** Print a double type array to a file in block format. */
void print_array_d(FILE *out,
		   const double *array,
		   const int length,
		   const int nrow,
		   const int ncol);

/** Print an integer type array to a file in block format. */
void print_array_i(FILE *out,
		   const int *array,
		   const int length,
		   const int nrow,
		   const int ncol);

/** Print an long type array to a file in block format. */
void print_array_l(FILE *out,
		   const long *array,
		   const int length,
		   const int nrow,
		   const int ncol);

/** Print the elements of a MATERIAL element to a file */
void print_material(FILE *out,
		    const MATERIAL *mat);

/** Update the bubble dofs on the elements */
void update_elem_bub_dofs(const long ne,
			  ELEMENT *const elem);

/* Compute the displacement gradient from a specified node 
   (0 <= node < nne) or from all nodes (node < 0 || node > nne) */
void compute_disp_grad(const int nne,
		       const double *ST,
		       const double *disp,
		       double *grad,
		       const int node);

/** Detemine the number of COMMUNICATION boundary elements and create
    a list of their indices */
long* list_boundary_el(const long ne,
		       const ELEMENT *elem,
		       const long nn,
		       const NODE *node,
		       const long myrank,
		       long *nbndel);

/** Get the output times from the input file. */
long* times_print (FILE *in1,
		   const long nt,
		   const long n_p);

/**
 * Get the global partition (Gf) of the data (f)for the process.
 *
 * \param[out] Gf Contains the global DOF values owned by the domin in Gid-order
 *
 * Side effects: point-to-point communication according to comm.
 */
void LToG (const double *f,
	   double *Gf,
	   const int myrank,
	   const int nproc,
	   const long ndofd,
	   const long *DomDof,
	   const long GDof,
	   const COMMUN comm,
	   const MPI_Comm mpi_comm);

/**
 * Get the local part (r) of the global data (Gr).
 *
 * \param[out] r Contains the DOF values on the domain, including
 * information from other domains.
 *
 * Side effects: point-to-point communication according to comm.
 */
void GToL (const double *Gr,
	   double *r,
	   const int myrank,
	   const int nproc,
	   const long ndofd,
	   const long *DomDof,
	   const long GDof,
	   const COMMUN comm,
	   const MPI_Comm mpi_comm);

MPI_Comm* CreateGraph (int nproc,
		       int myrank,
		       long nn,
		       NODE *node,
		       MPI_Comm mpi_comm);

/** Pause for t seconds */
void pause_time(int t);

/** Increases the length of a vector by a specified amount.  Used
    primarily in RNPsparse_ApAi. NOTE: frees the pointer 'orig'.
    Usage: original vector = increase_length(original vector, original
    length, new length) */
long* change_length(long *orig,
		    const long old_len,
		    const long new_len);

/** Checks whether the array is empty and exits */
void null_quit(void *array,
	       int error);

/**************/
/*** SINGLE ***/
/**************/

long num_fib (long nmat,
	      long ne,
	      ELEMENT *elem);

long num_matr (long nmat,
	       long ne,
	       ELEMENT *elem);

long list (long ***a,
	   long ne,
	   long nmat,
	   long nc,
	   ELEMENT *elem);

/** Compute the volume of a linear tetrahedron */
double Tetra_V (const double *x,
		const double *y,
		const double *z);

/** Compute the volume of a quadradic tetrahedron */
double Tetra_qv_V (const long nne,
		   const long ndofn,
		   const double *x,
		   const double *y,
		   const double *z);

/** Compute the volume of a hexahedron */
double Hexa_V (const double *x,
	       const double *y,
	       const double *z);

/** Returns the deformation on the element (using r and sup) in r_e. */
void def_elem (const long *cn,
	       const long ndofe,
	       const double *r,
	       const ELEMENT *elem,
	       const NODE *node,
	       double *r_e,
	       const SUPP sup,
	       const long TYPE);

/** Returns the local node numbers in a given element in nod[]. */
void elemnodes (const long ii,
		const long nne,
		long *nod,
		const ELEMENT *elem);

/** Returns the coordinates of the nodes on the element in element
    connectivity order. NOTE: returns undeformed configuration
    if(periodic == 1 || analysis == DISP). */
/* void nodecoord (const long nne, */
/* 		const long *nod, */
/* 		const NODE *node, */
/* 		double *x, */
/* 		double *y, */
/* 		double *z); */

/** returns node coords for total Lagrangian formulation
    (i.e. undeformed) */
void nodecoord_total (const long nne,
		      const long *nod,
		      const NODE *node,
		      double *x,
		      double *y,
		      double *z);

/** returns node coords for updated Lagrangian formulation
    (i.e. deformed config) */
void nodecoord_updated (const long nne,
			const long *nod,
			const NODE *node,
			double *x,
			double *y,
			double *z);

/** Determines the elements which have prescribed nodes and adds the
    indices to a list in the SUPP object.*/
void list_el_prescribed_def(SUPP sup,
			    const NODE *node,
			    const ELEMENT *elem,
			    const BOUNDING_ELEMENT *b_elems,
			    const long ne,
			    const int n_be,
			    const long nn);

void fun_eps (double *r_e,
	      long ndofe,
	      double **B_T,
	      double *eps);

void eps_element(long nne,
		 long ndofn,
		 double V,
		 double *r_e,
		 double *x,
		 double *y,
		 double *z,
		 double *EPSi);

void stress (long ne,
	     long ndofn,
	     NODE *node,
	     ELEMENT *elem,
	     MATGEOM matgeom,
	     HOMMAT *hommat,
	     double *r,
	     SIG *sig,
	     EPS *eps,
	     SUPP sup,
	     const int analysis);

void Mises_sig (long ne,
		SIG *sig,
		long TYPE);

void Mises_eps (long ne,
		EPS *eps,
		long TYPE);

/********************************************/
/****** SMOOTHING OF STRESSES TO NODES ******/
/********************************************/

void str_elem_matrix (long kk,
		      long nne,
		      long ndofn,
		      double *x,
		      double *y,
		      double *z,
		      double *K);

void str_proj_matrix (long *adr,
		      long ne,
		      long ndofn,
		      ELEMENT *elem,
		      NODE *node,
		      HOMMAT *hommat,
		      double *k,
		      const int analysis);

void stress_projector (long nne,
		       double *N,
		       double *P);

void str_solve (double *r,
		double *k,
		double *s,
		double *f,
		long *adr,
		long smo,
		long ne,
		long nn,
		long ndofn,
		NODE *node,
		ELEMENT *elem,
		HOMMAT *hommat,
		SIG *sig_e,
		SIG *sig_n,
		SUPP sup,
		const int analysis);

void str_prj_load (long ii,
		   long kk,
		   long nne,
		   long ndofn,
		   double *r_e,
		   double **D,
		   double *x,
		   double *y,
		   double *z,
		   SIG *sig_e,
		   double *f_e,
		   const int analysis);


/************************************************************************/

double eq_M_sig (long i,
		 long ip,
		 SIG *sig,
		 long TYPE);

double eq_M_eps (long i,
		 long ip,
		 EPS *eps,
		 double **deps_i,
		 long TYPE);

void eps_e_in (long nne,
	       long ndofn,
	       double *r_e,
	       double *x,
	       double *y,
	       double *z,
	       double **EPSi);

void unequal_forces (long ii,
		     double *x,
		     double *y,
		     double *z,
		     long nne,
		     long ndofn,
		     ELEMENT *elem,
		     double **dsig,
		     double *fe);

void aver_stress (long ii,
		  long nne,
		  long ndofn,
		  double *x,
		  double *y,
		  double *z,
		  SIG *sig,
		  EPS *eps);

/*******************************************************************************************/
/*****************************  ARC-LENGTH PROCEDURES  *************************************/
/*******************************************************************************************/

long diag_K (double *k,
	     long *adr,
	     long ndofd);

double det_K (double *k,
	      long *adr,
	      long ndofd);

/*
  double d_lam_ALM (long ndofd,
  double *rr,
  double *R,
  double dAL,
  double DET,
  double DET0,
  double dlm0,
  long PD,
  long PD0);

  double D_lam_ALM (long ndofd,
  double *rr,
  double *d_r,
  double *D_R,
  double *R,
  double dlm,
  double dAL);
*/

double new_arc_length (long iter,
		       long iter_des,
		       double dAL,
		       double dAL0);

/***************************************************************/

void check_equi (double *fu,
		 long ne,
		 long ndofd,
		 long ndofn,
		 ELEMENT *elem,
		 NODE *node,
		 MATGEOM matgeom,
		 SIG *sig,
		 const int analysis);

double* Energy_functional (long ne,
			   long ndofn,
			   long ndofd,
			   ELEMENT *elem,
			   NODE *node,
			   SIG *sig,
			   EPS *eps,
			   MATGEOM matgeom,
			   double *f,
			   double *r,
			   const int analysis);

/*********************  NONSYMMETRIC SPARSE SOLVER  ****************/

long* sparse_ApAi (long ne,
		   long ndofd,
		   long ndofn,
		   ELEMENT *elem,
		   NODE *node,
		   long *Ap);

/****************************************************************/

void tensor_9x9 (double **K,
		 double A[3][3][3][3],
		 long pom);


double equivalent_Mises (long i,
			 SIG *sig);

double equivalent_M_eps (long i,
			 EPS *eps);

double equivalent_M_eps_pl (long i,
			    EPS *eps);

void Mises (long ne,
	    SIG *sig,
	    EPS *eps,
	    const int analysis);

double T_VOLUME (const long ne,
		 const long ndofn,
		 const ELEMENT *elem,
		 const NODE *node);

double area (long nne,
	     double *x,
	     double *y);

void Logarithmic_strain (double **F,
			 double **EL);

#endif
