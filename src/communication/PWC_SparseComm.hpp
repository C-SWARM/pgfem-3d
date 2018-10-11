#include "pgfem3d/Communication.hpp"

using pgfem3d::SparseComm;

namespace pgfem3d {
  
class PWC_SparseComm : public SparseComm {
public:
  PWC_SparseComm(pgfem3d::net::PWCNetwork *n, pgfem3d::net::PGFem3D_Comm c);
  ~PWC_SparseComm();
  
  void initialize();
  void post_stiffmat(double ***Lk, double ***receive);
  void print_stiffmat();
  void send_stiffmat();
  void finalize_stiffmat();
  void assemble_nonlocal_stiffmat(pgfem3d::solvers::SparseSystem* system);

  void LToG(const double *f, double *Gf, const long ndofd,
	    const long *DomDof, const int GDof);
  void GToL(const double *Gr, double *r, const long ndofd,
	    const int GDof);
  
private:
  // stiffmat
  size_t total_ssz;
  size_t total_rsz;

  double *backing_s;
  double *backing_r;

  net::Buffer *local_s;
  net::Buffer *local_r;
  net::Buffer *remote;
  
  // LToG and GToL
  size_t total_lg_ssz;
  size_t total_lg_rsz;

  double **lg_S;
  double **lg_R;
  
  double *backing_lg_S;
  double *backing_lg_R;

  net::Buffer *local_lg_S;
  net::Buffer *local_lg_R;
  net::Buffer *remote_lg_S;
  net::Buffer *remote_lg_R;
  
  void exchange(net::Buffer *in, net::Buffer *out, long nsend,
		long nrecv, long *idxs);
  
  net::PWCNetwork *net;
  static constexpr net::CID LOCAL_ID = 0xcafebabe;
};
  
} // namespace pgfem3d
