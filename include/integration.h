#ifndef PGFEM3D_INTEGRATION_H
#define PGFEM3D_INTEGRATION_H

#include "PGFem3D_options.h"
#include "crpl.h"
#include "data_structure.h"
#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "matgeom.h"
#include "node.h"
#include "sig.h"
#include "supp.h"

/** Explicit hardening integration. */
double int_har_explicit (long ii,
                         long ip,
                         long mat,
                         SIG *sig,
                         CRPL *crpl,
                         double dt);

double int_alg_res (double lam,
                    double dt,
                    long mat,
                    long nss,
                    double *TA,
                    double *GA,
                    CRPL *crpl,
                    double **UU,
                    double **Fp,
                    double **BB);

long int_UU (long ii,
             long ip,
             double NOR_MIN,
             double dt,
             long mat,
             double *NOR,
             long nss,
             CRPL *crpl,
             double **Fr,
             double **FnB,
             double **Fp,
             double **BB,
             double *DD,
             double **S,
             double L[3][3][3][3],
             double *TA,
             double *GA,
             double Tr,
             double Jr,
             double Har,
             double *LAM,
             double **UU,
             long GAMA,
             long iter);

long int_HA (long ii,
             long ip,
             long nss,
             long mat,
             double dt,
             CRPL *crpl,
             double *GA,
             double Har,
             double *HAR,
             EPS *eps,
             const PGFem3D_opt *opts);

long integration_alg (long ne,
                      long ndofn,
                      long ndofd,
                      long npres,
                      CRPL *crpl,
                      Element *elem,
                      Node *node,
                      double *d_r,
                      double *rr,
                      SUPP sup,
                      MATGEOM matgeom,
                      HOMMAT *hommat,
                      EPS *eps,
                      SIG *sig,
                      long tim,
                      long iter,
                      double dt,
                      double nor_min,
                      long STEP,
                      long GAMA,
                      const PGFem3D_opt *opts,
                      const int mp_id);

#endif /* #define PGFEM3D_INTEGRATION_H */
