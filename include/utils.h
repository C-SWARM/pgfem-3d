/**
 * @file This file declares utility functions for PGFem3D.
 */
#ifndef PGFEM3D_UTILS_H
#define PGFEM3D_UTILS_H

#include "PGFem3D_data_structure.h"
#include "bounding_element.h"
#include "data_structure.h"
#include "element.h"
#include "hommat.h"
#include "material.h"
#include "matgeom.h"
#include "node.h"
#include "sig.h"
#include <cstdio>

using pgfem3d::solvers::SparseSystem;
struct EPS;                                     // defined in eps.h

/// compute Eulerian Almansi strain from a given F
void compute_Eulerian_Almansi_strain(double *e,
                                     double *F);

/// compute Equivalent (Von Mises) Eulerian Almansi strain from a given F
double compute_Equivalent_strain(double *F_in);

/// Lots of system headers require that the return value from scanf be
/// checked. This small utility simplifies this process by counting the number
/// of arguments that are expected automatically.
#define CHECK_SCANF(stream, ...) do {                                   \
    check_scanf(__FILE__, __LINE__, __func__, (stream), __VA_ARGS__);   \
  } while (0)

template <class... Args>
static inline void check_scanf(const char* file, int line, const char* func,
                               FILE *stream, Args... args)
{
  constexpr int N = (sizeof...(args) - 1);   // the number of expected arguments
  const int n = fscanf(stream, args...);     // the number of actual arguments
  if (n != N) {
    std::string e;
    e += file;
    e += ":";
    e += line;
    e += " (";
    e += func;
    e += ") ";
    e += "Expected " + std::to_string(N) + "args, saw " + std::to_string(n);
    throw std::runtime_error(std::move(e));
  }
}

/// Simple utility to convince the compiler to tell us the length of an array.
template <class T, int N>
constexpr int size(T (&)[N]) {
  return N;
}

/// Simple clock macro to get the current time as a double
static inline double CLOCK() {
  timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec/1e9;
}

/**
 * Scan the file for a valid line (non-blank and does not start with a
 * '#').  This function may be called multiple times on the same file.
 *
 * \param[in,out] in, File to scan
 *
 * \return non-zero on error. Upon successful completion, 'in' is
 * returned with the file position set to the beginning of the valid
 * line.
 */
int scan_for_valid_line(FILE *in);

void pack_2mat(const void **src,
               const int nrow,
               const int ncol,
               const size_t elem_size,
               char *buffer,
               size_t *pos);

void unpack_2mat(void **dest,
                 const int nrow,
                 const int ncol,
                 const size_t elem_size,
                 const char *buffer,
                 size_t *pos);

void pack_3mat(const void ***src,
               const int n_1,
               const int n_2,
               const int n_3,
               const size_t elem_size,
               char *buffer,
               size_t *pos);

void unpack_3mat(void ***dest,
                 const int n_1,
                 const int n_2,
                 const int n_3,
                 const size_t elem_size,
                 const char *buffer,
                 size_t *pos);

void pack_4mat(const void ****src,
               const int n_1,
               const int n_2,
               const int n_3,
               const int n_4,
               const size_t elem_size,
               char *buffer,
               size_t *pos);

void unpack_4mat(void ****dest,
                 const int n_1,
                 const int n_2,
                 const int n_3,
                 const int n_4,
                 const size_t elem_size,
                 const char *buffer,
                 size_t *pos);

/** src and dest must be unique */
void copy_2mat(void **dest,
               const void **src,
               const int nrow,
               const int ncol,
               const size_t elem_size);

/** src and dest must be unique */
void copy_3mat(void ***dest,
               const void ***src,
               const int n_1,
               const int n_2,
               const int n_3,
               const size_t elem_size);

/** src and dest must be unique */
void copy_4mat(void ****dest,
               const void ****src,
               const int n_1,
               const int n_2,
               const int n_3,
               const int n_4,
               const size_t elem_size);

/**
 * \brief Determine the number of duplicate values in an array.
 *
 * A copy of arr is sorted according to the compare function. The
 * sorted array is then checked for duplicates using the compare
 * function again.
 *
 * \return number of duplicate values.
 */
int number_of_duplicates(const void *arr,
                         const size_t n_elem,
                         const size_t size,
                         int (*compare)(const void *a, const void *b));

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

/**
 * Compute the solution to Ax=b for an NxN system using LAPACK and
 * return 0 on success.
 *
 * \param[in] n_eq, the number of equations to solve
 * \param[in] mat_dim, the size of the (square) matrix
 * \param[in] A, the coefficient matrix
 * \param[in,out] b_x, on entry: RHS (b), on exit: solution (x)
 * \return non-zero on recoverable error, abort on logic error
 *
 * CAVEATES: The matrix A is copied/transposed for migration to
 * FORTRAN LAPACK routine.
 */
int solve_Ax_b(const int n_eq,
               const int mat_dim,
               const double *A,
               double *b_x);

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
                   const Ai_t *array,
                   const int length,
                   const int nrow,
                   const int ncol);

/** Print an long type array to a file in block format. */
void print_array_l(FILE *out,
                   const long *array,
                   const int length,
                   const int nrow,
                   const int ncol);

/** Print the elements of a Material element to a file */
void print_material(FILE *out,
                    const Material *mat);

/** Update the bubble dofs on the elements */
void update_elem_bub_dofs(const long ne,
                          Element *const elem);

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
                       const Element *elem,
                       const long nn,
                       const Node *node,
                       const int myrank,
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
       const long ndofd,
       const pgfem3d::CommunicationStructure *com);

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
       const long ndofd,
       const pgfem3d::CommunicationStructure *com);

multiscale::net::MSNET_Comm*
CreateGraph (int nproc,
	     int myrank,
	     long nn,
	     Node *node,
	     multiscale::net::MSNET_Comm comm);

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
              Element *elem);

long num_matr (long nmat,
               long ne,
               Element *elem);

long list (long ***a,
           long ne,
           long nmat,
           long nc,
           Element *elem);

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
               const Element *elem,
               const Node *node,
               double *r_e,
               const SUPP sup,
               const long TYPE);

/** Compute the TOTAL deformation on an element using r, d_r, and
    sup->{defl,defl_d} and store in r_e */
void def_elem_total (const long *cn,
                     const long ndofe,
                     const double *r,
                     const double *d_r,
                     const Element *elem,
                     const Node *node,
                     const SUPP sup,
                     double *r_e);

/// compute value of nodal variables
///
/// \param[in] cn id of nodal values
/// \param[in] ndofe number of degree of freedom on an element
/// \param[in] r nodal variables at n+1
/// \param[in] d_r nodal variable increments at n+1
/// \param[in] elem Element object
/// \param[in] node Node object
/// \param[out] r_e computed nodal variables for an element
/// \param[in] reference nodal value
/// \return non-zero on interal error
int def_elem_with_reference(const long *cn,
                            const long ndofe,
                            const double *r,
                            const double *d_r,
                            const Element *elem,
                            const Node *node,
                            const SUPP sup,
                            double *r_e,
                            double r0);

/** Returns the local node numbers in a given element in nod[]. */
void elemnodes (const long ii,
                const long nne,
                long *nod,
                const Element *elem);

/** returns node coords for total Lagrangian formulation
    (i.e. undeformed) */
void nodecoord_total (const long nne,
                      const long *nod,
                      const Node *node,
                      double *x,
                      double *y,
                      double *z);

/** returns node coords for updated Lagrangian formulation
    (i.e. deformed config) */
void nodecoord_updated (const long nne,
                        const long *nod,
                        const Node *node,
                        double *x,
                        double *y,
                        double *z);

/** Determines the elements which have prescribed nodes and adds the
    indices to a list in the SUPP object.*/
void list_el_prescribed_def(SUPP sup,
                            const Node *node,
                            const Element *elem,
                            const BoundingElement *b_elems,
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
             Node *node,
             Element *elem,
             MATGEOM matgeom,
             HOMMAT *hommat,
             double *r,
             SIG *sig,
             EPS *eps,
             SUPP sup,
             const int analysis,
             const int mp_id);

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
                      Element *elem,
                      Node *node,
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
                Node *node,
                Element *elem,
                HOMMAT *hommat,
                SIG *sig_e,
                SIG *sig_n,
                SUPP sup,
                const int analysis,
                const int mp_id);

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
                     Element *elem,
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

/***************************************************************/
/*************  ARC-LENGTH PROCEDURES  *************************/
/***************************************************************/

long diag_K (double *k,
             long *adr,
             long ndofd);

double det_K (double *k,
              long *adr,
              long ndofd);

double new_arc_length (long iter,
                       long iter_des,
                       double dAL,
                       double dAL0);

/***************************************************************/

void check_equi (double *fu,
                 long ne,
                 long ndofd,
                 long ndofn,
                 Element *elem,
                 Node *node,
                 MATGEOM matgeom,
                 SIG *sig,
                 const int analysis,
                 const int mp_id);

double* Energy_functional (long ne,
                           long ndofn,
                           long ndofd,
                           Element *elem,
                           Node *node,
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
                   Element *elem,
                   Node *node,
                   long *Ap,
                   const int mp_id);

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
                 const Element *elem,
                 const Node *node);

double area (long nne,
             double *x,
             double *y);

void Logarithmic_strain (double **F,
                         double **EL);
void mid_point_rule(double *v, double *w, double *x, double alpha, long n_row);

/// determine whether the element is on communication boundary or in interior
///
/// \parma[in] eid element id
/// \param[in,out] idx id of bndel (communication boundary element)
/// \param[in,out] skip count element on communication boundary
/// \param[in] nbndel number of elements on communication boundary
/// \param[in] bndel list of elements on communcation boundary
/// \param[in] myrank current process rank
/// \return return 1 if the element is interior or 0 if the element on the communication boundary
int is_element_interior(int eid, int *idx, int *skip, long nbndel, long *bndel, int myrank);

double compute_volumes_from_coordinates(double *x,
                                        double *y,
                                        double *z,
                                        long nne);

/// find roots of cubic equations(a*x^3+ b*x^2 + c*x + d = 0) numerically.
void compute_root_of_cubic_euqation(double *x,
                                    double a,
                                    double b,
                                    double c,
                                    double d,
                                    double x0 = 0.0,
                                    bool print_cnvg = false);

/// compute sting math expression as a function of time
double string_function_of_time(const std::string &expr,
                               double t);
                               
/// compute sting math expression as a function of position (x, y, z)
double string_function_of_xyz(const std::string &expr,
                              double x,
                              double y,
                              double z);
                              
/// compute sting math expression as a function of time and position (x, y, z)
double string_function_of_txyz(const std::string &expr,
                               double t,
                               double x,
                               double y,
                               double z);

#endif // #define PGFEM3D_UTILS_H
