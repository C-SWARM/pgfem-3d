#ifndef _THREE_FIELD_ELEMENT_H_
#define _THREE_FIELD_ELEMENT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "femlib.h"
#include "elem3d.h"
#include "new_potentials.h"
//#include "supp.h"
#include "sig.h"
#include "eps.h"
#include "tensors.h"
#define N_VOL_TREE_FIELD 1
void resid_w_inertia_Ru_ip(double *fu,
        int nne, double *ST, double *F, double *S, double jj, double wt, double Pn);

void resid_w_inertia_Rt_ip(double *ft, int nVol, double jj, double wt, double *Nt, double Pn, double Up);

void resid_w_inertia_Rp_ip(double *fp, int npres, double *F, double jj, double wt, double *Np, double Tn);

void stiffmat_3f_el(double *Ks, const int ii, const int ndofn, const int nne, int npres, int nVol, int nsd,
              const double *x, const double *y, const double *z, const ELEMENT *elem, const HOMMAT *hommat,
              const long *nod, const NODE *node, double dt, SIG *sig, EPS *eps, const SUPP sup, double alpha, double *r_e);

void residuals_3f_el(double *f, const int ii, const int ndofn, const int nne, const int npres, const int nVol, const int nsd,
        const double *x, const double *y, const double *z, 
        const ELEMENT *elem, const HOMMAT *hommat, const long *nod, const NODE *node,
        double dt, SIG *sig, EPS *eps,  const SUPP sup, double *r_e); 

void residuals_3f_w_inertia_el(double *f,const int ii,
        const int ndofn,const int nne,const int npres,const int nVol,const int nsd,
        const double *x,const double *y,const double *z,
        const ELEMENT *elem,const HOMMAT *hommat,const NODE *node,
        const double *dts,SIG *sig,EPS *eps,double alpha, double *r_n_a, double *r_n_1_a);
void update_3f(long ne, long ndofn, long npres, double *d_r, double *r, double *rr,
               NODE *node, ELEMENT *elem, HOMMAT *hommat, SUPP sup, EPS *eps, SIG *sig, double dt, double t,
		           MPI_Comm mpi_comm, const PGFem3D_opt *opts, double alpha, double *r_n, double *r_n_1);
		           
void update_3f_state_variables(long ne, long ndofn, long npres, double *d_r, double *r,
               NODE *node, ELEMENT *elem, HOMMAT *hommat, SUPP sup, EPS *eps, SIG *sig, double dt, double t,
		           MPI_Comm mpi_comm);
		           		           
void compute_stress(double *GS, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm,int analysis);		           
#ifdef __cplusplus
}
#endif

#endif