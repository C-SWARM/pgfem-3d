#ifndef PGFEM3D_COMMUNICATION_H
#define PGFEM3D_COMMUNICATION_H

/// @brief This file defines the system-wide Communication abstractions
#include "pgfem3d/Boot.hpp"
#include "pgfem3d/Network.hpp"
#include "pgfem3d/Solver.hpp"
#include "datatype.hpp"
#include <cassert>
#include <cstdio>
#include <string>
#include <vector>

namespace pgfem3d {

/**
 * Structure for containing an index map.
 *
 * \param n_map The number of mapped indices.
 * \param maps The index mapping stored in two consecutive
 * elements. Example: the i-th map is INDEX maps[i*2] maps to INDEX
 * maps[i*2+1].
 */
struct fast_map {
  int n_maps = 0;             /**< number of maps */
  size_t *maps = nullptr;     /**< mapping stored in two consecutive elements */
};

/// The communication hints class.
///
/// In PGFem3d, communication hints are lists of send and receive values that
/// are provided as per-rank input files. The CommHints class reads those hints
/// and provides simplified interfaces to the resulting vectors.
class CommHints {
 public:
  /// Create the name of the comm hints file for a rank.
  static std::string BuildFilename(const char *ipath,
                                   const char* basefname,
                                   int rank);

  /// Construct a comm hints class with the passed path and base filename.
  ///
  /// This will throw a system_error exception if it fails to find a file.
  CommHints(const char *ipath, const char* basefname, int rank);

  size_t nSend() const {
    return send_.size();
  }

  const std::vector<int>& sends() const {
    return send_;
  }

  const int* get_send_list() const {
    return &send_[0];
  }

  size_t nRecv() const {
    return recv_.size();
  }

  const std::vector<int>& recvs() const {
    return recv_;
  }

  const int* get_recv_list() const {
    return &recv_[0];
  }

 private:
  std::vector<int> send_ = {};
  std::vector<int> recv_ = {};
};

// The Sparse Communication class
class SparseComm {
 public:
  virtual ~SparseComm();

  static SparseComm* Create(net::Network *net, net::PGFem3D_Comm);

  virtual void initialize() = 0;
  virtual void post_stiffmat(double ***Lk, double ***receive) = 0;
  virtual void send_stiffmat() = 0;
  virtual void finalize_stiffmat() = 0;
  virtual void assemble_nonlocal_stiffmat(pgfem3d::solvers::SparseSystem* system) = 0;

  virtual void LToG(const double *f, double *Gf, const long ndofd,
                    const long *DomDof, const long GDof) = 0;
  virtual void GToL(const double *Gr, double *r, const long ndofd,
                    const long GDof) = 0;

  void get_owned_global_dof_values(const double *local_dofs,
                                   double *global_dofs);
  void get_local_dof_values_from_global(const double *global_dofs,
                                        double *local_dofs);

  int build_fast_maps(const long ndofd,
                      const long ngdof_owned,
                      const long start_gdof_id);


  long *S,   /**< Contains how many rows to send to each processor
                (nproc) */
    *R,      /**< Contains how many rows to recieve from each processor
                (nproc) */
    *AS,     /**< Amount to send (nproc) */
    *AR,     /**< Amount to receive (nproc)*/
    **SLID,  /**< Local ID of communicated rows (nproc)(S[rank]) */
    **RGID,  /**< Global ID of rows to recieve (nproc)(R[rank]) */
    *LG,     /**< Local-to-global (ndofd) */
    *GL,     /**< Global-to-local (num global on dom) */
    **SAp,   /**< nnz in row to send (nproc)(S[rank]) */
    **SGRId, /**< Global row IDs to be sent (nproc)() */
    **RAp,   /**< nnz of each row to recieve (nproc)(R[rank]) */
    **RGRId, /**< Global row IDs to receive (nproc)() */
    Ns,      /**< Number of procs to send to*/
    Nr,      /**< Number of procs to recieve from*/
    *Nss,    /**< Which procs to send to (Ns) */
    *Nrr;    /**< Which procs to receive from (Nr) */

 protected:
  static void get_mapped_values(const fast_map *map, const double *A, double *B)
  {
    for(int i=0; i<map->n_maps; i++){
      B[map->maps[2*i+1]] = A[map->maps[2*i]];
    }
  }

  int build_fast_LG_map(const long ndofd,
                        const long ngdof_owned,
                        const long start_gdof_id);
  int build_fast_GL_map(const long ndofd,
                        const long ngdof_owned);

  double **send;            //!< The local buffer to send
  double **recv;            //!< The receive buffer

  fast_map *fast_LG_map;    //!< fast map of the LG indices */
  fast_map *fast_GL_map;    //!< fast map of the GL indices */

  int nproc;                //!< number of processes
  int myrank;               //!< this process rank

  net::Network *net;        //!< Network handle for data exchange
  net::PGFem3D_Comm comm;   //!< The PGFem3D network communicator handle
};


// The primary Communication structure
struct CommunicationStructure {
  virtual ~CommunicationStructure();

  int nproc;                //!< number of processes
  int rank;                 //!< this process rank
  int *Ap;                  //!< n_cols in each owned row of global stiffness matrix
  Ai_t *Ai;                 //!< column ids for owned rows of global stiffness matrix
  long *DomDof;             //!< number of global DOFs on each domain
  long nbndel;              //!< number of Element on the communication boundary
  long *bndel;              //!< Element ids on the communication boundary
  long GDof;                //!< maximum id of locally owned global DOF
  long NBN;                 //!< Number of nodes on domain interfaces

  net::PGFem3D_Comm comm;   //!< The PGFem3D network communicator handle

  net::Network *net;        //!< Network handle for data exchange
  net::Boot *boot;          //!< Boot handle for PMI info
  SparseComm *spc;          //!< Sparse communication handle
  CommHints *hints;         //!< Communication hints handle
};

/** Wrappers for network abort. Allows for simple change
    throughout code for different abort mechanism */
[[noreturn]] void PGFEM_Comm_code_abort(const CommunicationStructure *com, int code);
[[noreturn]] void PGFEM_Comm_abort(const CommunicationStructure *com);
[[noreturn]] void PGFEM_Abort();

/** Get the rank on COMM_WORLD for error messages */
int PGFEM_Error_rank(int* myrank);

} // namespace pgfem3d

#endif
