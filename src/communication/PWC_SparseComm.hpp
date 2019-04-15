#include "pgfem3d/Communication.hpp"

using pgfem3d::SparseComm;

namespace pgfem3d {
  
class PWC_SparseComm : public SparseComm {
public:
  PWC_SparseComm(multiscale::net::PWCNetwork *n, multiscale::net::MSNET_Comm c);
  ~PWC_SparseComm();
  
  void initialize();
  void post_stiffmat(double ***Lk, double ***receive);
  void send_stiffmat();
  void finalize_stiffmat();
  void assemble_nonlocal_stiffmat(pgfem3d::solvers::SparseSystem* system);

  void LToG(const double *f, double *Gf, const long ndofd,
	    const long *DomDof, const long GDof);
  void GToL(const double *Gr, double *r, const long ndofd,
	    const long GDof);
  
private:
  // stiffmat
  size_t total_ssz;
  size_t total_rsz;

  double *backing_s;
  double *backing_r;

  multiscale::net::Buffer *local_s;
  multiscale::net::Buffer *local_r;
  multiscale::net::Buffer *remote;
  
  // LToG and GToL
  size_t total_lg_ssz;
  size_t total_lg_rsz;

  double **lg_S;
  double **lg_R;
  
  double *backing_lg_S;
  double *backing_lg_R;

  multiscale::net::Buffer *local_lg_S;
  multiscale::net::Buffer *local_lg_R;
  multiscale::net::Buffer *remote_lg_S;
  multiscale::net::Buffer *remote_lg_R;
  
  void exchange(multiscale::net::Buffer *in,
		multiscale::net::Buffer *out,
		long nsend, long nrecv, long *idxs);
  
  multiscale::net::PWCNetwork *net;
  static constexpr multiscale::net::CID LOCAL_ID = 0xcafebabe;
};
  
} // namespace pgfem3d
