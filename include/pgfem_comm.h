/* HEADER */
#ifndef PGFEM_COMM_H
#define PGFEM_COMM_H

#ifndef PGFEM_MPI_H
#include "PGFEM_mpi.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

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
  };
  typedef struct COMMUN_1 COMMUN_1;
  typedef struct COMMUN_1 *COMMUN;

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

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
