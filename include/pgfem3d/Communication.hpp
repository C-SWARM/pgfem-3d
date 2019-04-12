#ifndef PGFEM3D_COMMUNICATION_H
#define PGFEM3D_COMMUNICATION_H

/// @brief This file defines the system-wide Communication abstractions
#include "pgfem3d/Boot.hpp"
#include "pgfem3d/Network.hpp"
#include "pgfem3d/Solver.hpp"
#include <cstdio>

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
  
// Tracks communication hints
class CommHints {
public:
  CommHints() : filestr(0), send(0), recv(0), nsend(0), nrecv(0) {}
  CommHints(const char *ipath, const char* basefname, const int rank);
  ~CommHints();
  
  void read(FILE* in);
  void write(FILE *out);
  void read_filename(const char *fn);
  void write_filename(const char *fn);

  char* get_filename() const { return filestr; };
  int get_nsend() const {return nsend;};
  int get_nrecv() const {return nrecv;};
  int* get_send_list() const {return send;};
  int* get_recv_list() const {return recv;};
  
private:
  char *filestr;
  int *send;
  int *recv;
  int nsend;
  int nrecv;
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
		    const long *DomDof, const int GDof) = 0;
  virtual void GToL(const double *Gr, double *r, const long ndofd,
		    const int GDof) = 0;
  
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
  int *Ai;                  //!< column ids for owned rows of global stiffness matrix
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
