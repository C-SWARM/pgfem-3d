/* HEADER */
/**
 * AUTHORS:
 * Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 */
#pragma once
#ifndef PGFEM_COMM_H
#define PGFEM_COMM_H

#include "PGFEM_mpi.h"

  /** Structure for MPI communication */
  struct COMMUN_1{
    long *S, /**< Contains how many rows to send to each processor
        (nproc) */
      *R, /**< Contains how many rows to recieve from each processor
         (nproc) */
      *AS, /**< Amount to send (nproc) */
      *AR, /**< Amount to recieve (nproc)*/
      **SLID, /**< Local ID of communicated rows (nproc)(S[rank]) */
      **RGID, /**< Global ID of rows to recieve (nproc)(R[rank]) */
      *LG, /**< Local-to-global (ndofd) */
      *GL, /**< Global-to-local (num global on dom) */
      **SAp, /**< nnz in row to send (nproc)(S[rank]) */
      **SGRId, /**< Global row IDs to be sent (nproc)() */
      **RAp, /**< nnz of each row to recieve (nproc)(R[rank]) */
      **RGRId, /**< Global row IDs to receive (nproc)() */
      Ns, /**< Number of procs to send to*/
      Nr, /**< Number of procs to recieve from*/
      *Nss, /**< Which procs to send to (Ns) */
      *Nrr; /**< Which procs to receive from (Nr) */

    void *fast_LG_map; /**< fast map of the LG indices */
    void *fast_GL_map; /**< fast map of the GL indices */
  };
  typedef struct COMMUN_1 COMMUN_1;
  typedef struct COMMUN_1 *COMMUN;

  /** Set initial values for a COMMUN structure */
  void initialize_commun(COMMUN comm);

  /** Destroy the communication graph object */
  void destroy_commun(COMMUN comm,long nproc);

  /** builds the off-process buffer for assembling the gloabal stiffness
      matrix. Posts non-blocking recieve */
  int init_and_post_stiffmat_comm(double ***Lk,
                  double ***receive,
                  MPI_Request **req_r,
                  MPI_Status **sta_r,
                  const MPI_Comm mpi_comm,
                  const COMMUN pgfem_comm);

  /** Send the off-process data for stiffness assembly */
  int send_stiffmat_comm(MPI_Status **sta_s,
             MPI_Request **req_s,
             double **Lk,
             const MPI_Comm mpi_comm,
             const COMMUN pgfem_comm);

  /** Finalize the stiffmat communication. Simple wrapper for the
      waitall commands. */
  int finalize_stiffmat_comm(MPI_Status *sta_s,
                 MPI_Status *sta_r,
                 MPI_Request *req_s,
                 MPI_Request *req_r,
                 const COMMUN pgfem_comm);

    typedef struct PGFEM_HYPRE_solve_info PGFEM_HYPRE_solve_info;

  /** Assemble non-local parts as they arrive */
  int assemble_nonlocal_stiffmat(const COMMUN pgfem_comm,
                   MPI_Status *sta_r,
                   MPI_Request *req_r,
                   PGFEM_HYPRE_solve_info *PGFEM_hypre,
                   double **recv);

  /**
   * Build the fast LG and GL index maps.
   *
   * \param[in/out] pgfem_comm Communicator object. Accesses and
   * modifies internal values.
   * \praram[in] ndofd Number of degrees of freedom on the domain.
   * \param[in] ngdof_owned Number of OWNED global dofs on the domain
   * (DomDof[myrank])
   * \param[in] start_gdof_id Starting global dof id owned by the domain (GDof).
   *
   * \return non-zero on error
   */
  int pgfem_comm_build_fast_maps(COMMUN pgfem_comm,
                 const long ndofd,
                 const long ngdof_owned,
                 const long start_gdof_id);

  /**
   * Set the locally owned GLOBAL degrees of freedom VALUES from the local dof array.
   *
   * \param[in] pgfem_comm Commuication structure containing
   * information on the mapping.
   * \param[in] local_dofs Array of local dof values (size ndofd)
   * \param[out] global_dofs Array of global dof values (size
   * DomDof[myrank]). Contains locally owned values on return, others
   * must be communicated.
   *
   * \return non-zero on error
   *
   */
  int pgfem_comm_get_owned_global_dof_values(const COMMUN pgfem_comm,
                         const double *local_dofs,
                         double *global_dofs);

  /**
   * Set the LOCAL dof VALUES from the global dof array.
   *
   * \param[in] pgfem_comm Commuication structure containing
   * information on the mapping.
   * \param[in] global_dofs Array of global dof values (size DomDof[myrank])
   * \param[out] local_dofs Array of loacal dof values (size  ndofd).
   * Contains locally owned values on return, others must be communicated.
   *
   * \return non-zero on error
   *
   */
  int pgfem_comm_get_local_dof_values_from_global(const COMMUN pgfem_comm,
                          const double *global_dofs,
                          double *local_dofs);

#endif /* #ifndef  */
