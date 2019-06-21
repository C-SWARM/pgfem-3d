#include "pgfem3d/Communication.hpp"

using pgfem3d::SparseComm;

namespace pgfem3d {
  
class ISIR_SparseComm : public SparseComm {
public:
  ISIR_SparseComm(pgfem3d::net::ISIRNetwork *n, pgfem3d::net::PGFem3D_Comm c);
  ~ISIR_SparseComm();
  
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
  net::ISIRNetwork *net;
  net::Request *req_s, *req_r;
  net::Status *sta_s, *sta_r;
};
  
} // namespace pgfem3d
