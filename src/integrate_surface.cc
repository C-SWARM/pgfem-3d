/* HEADER */
#include "integrate_surface.h"
#include <string.h>
#include <math.h>
#include "mkl_cblas.h"

#include "PGFEM_mpi.h"
#include "PGFEM_io.h"
#include "quadrature_rules.h"
#include "cohesive_element_utils.h"
#include "allocation.h"
#include "utils.h"

int integrate_surface(const int nne,
              const int face_id,
              const int int_order,
              int *n_ip,
              double **ksi_3D,
              double **eta_3D,
              double **zet_3D,
              double **ksi_2D,
              double **eta_2D,
              double **wt_2D,
              int *nne_2D,
              int **nod_2D)
{
  int err = 0;
  /* initialize variables */
  (*n_ip) = 0;
  (*ksi_3D) = NULL;
  (*eta_3D) = NULL;
  (*zet_3D) = NULL;
  (*ksi_2D) = NULL;
  (*eta_2D) = NULL;
  (*wt_2D) = NULL;
  (*nod_2D) = NULL;

  /* NOTE: face numbering is according to T3d user manual */
  switch(nne){
  case 4: /* Tet --> tria */
    /* get 2D integration points */
    err = get_tria_quadrature_rule(int_order,n_ip,ksi_2D,eta_2D,
                   wt_2D);
    *nne_2D = 3;
    *ksi_3D = PGFEM_calloc(double, *n_ip);
    *eta_3D = PGFEM_calloc(double, *n_ip);
    *zet_3D = PGFEM_calloc(double, *n_ip);
    *nod_2D = PGFEM_calloc(int, *nne_2D);
    /* assign 3D integration points */
    switch(face_id){
    case 0:
      memcpy(*ksi_3D,*ksi_2D,(*n_ip)*sizeof(double));
      memcpy(*eta_3D,*eta_2D,(*n_ip)*sizeof(double));
      /* zet_3D = 0 */
      (*nod_2D)[0] = 0;
      (*nod_2D)[1] = 2;
      (*nod_2D)[2] = 1;
      break;
    case 1:
      memcpy(*ksi_3D,*ksi_2D,(*n_ip)*sizeof(double));
      memcpy(*zet_3D,*eta_2D,(*n_ip)*sizeof(double));
      /* eta_3D = 0 */
      (*nod_2D)[0] = 0;
      (*nod_2D)[1] = 1;
      (*nod_2D)[2] = 3;
      break;
    case 2:
      memcpy(*ksi_3D,*ksi_2D,(*n_ip)*sizeof(double));
      memcpy(*eta_3D,*eta_2D,(*n_ip)*sizeof(double));
      for(int i=0; i<*n_ip; i++){
        *zet_3D[i] = 1 - *ksi_3D[i] - *eta_3D[i];
      }
      (*nod_2D)[0] = 1;
      (*nod_2D)[1] = 2;
      (*nod_2D)[2] = 3;
      break;
    case 3:
      memcpy(*eta_3D,*ksi_2D,(*n_ip)*sizeof(double));
      memcpy(*zet_3D,*eta_2D,(*n_ip)*sizeof(double));
      /* ksi_3D = 0 */
      (*nod_2D)[0] = 0;
      (*nod_2D)[1] = 3;
      (*nod_2D)[2] = 2;
      break;
    }/* end switch face */
    break;
  case 10:

    err = get_tria_quadrature_rule(int_order,n_ip,ksi_2D,eta_2D,
                   wt_2D);
    *nne_2D = 6;
    *ksi_3D = PGFEM_calloc(double, *n_ip);
    *eta_3D = PGFEM_calloc(double, *n_ip);
    *zet_3D = PGFEM_calloc(double, *n_ip);
    *nod_2D = PGFEM_calloc(int, *nne_2D);
    /* assign 3D integration points */
    switch(face_id){
    case 0:
      memcpy(*ksi_3D,*ksi_2D,(*n_ip)*sizeof(double));
      memcpy(*eta_3D,*eta_2D,(*n_ip)*sizeof(double));
      (*nod_2D)[0] = 0;
      (*nod_2D)[1] = 2;
      (*nod_2D)[2] = 1;
      (*nod_2D)[3] = 6;
      (*nod_2D)[4] = 5;
      (*nod_2D)[5] = 4;
      break;
    case 1:
      memcpy(*ksi_3D,*ksi_2D,(*n_ip)*sizeof(double));
      memcpy(*zet_3D,*eta_2D,(*n_ip)*sizeof(double));
      /* eta_3D = 0 */
      (*nod_2D)[0] = 0;
      (*nod_2D)[1] = 1;
      (*nod_2D)[2] = 3;
      (*nod_2D)[3] = 4;
      (*nod_2D)[4] = 8;
      (*nod_2D)[5] = 7;
      break;
    case 2:
      memcpy(*ksi_3D,*ksi_2D,(*n_ip)*sizeof(double));
      memcpy(*eta_3D,*eta_2D,(*n_ip)*sizeof(double));
      for(int i=0; i<*n_ip; i++){
    *zet_3D[i] = 1 - *ksi_3D[i] - *eta_3D[i];
      }
      (*nod_2D)[0] = 1;
      (*nod_2D)[1] = 2;
      (*nod_2D)[2] = 3;
      (*nod_2D)[3] = 5;
      (*nod_2D)[4] = 9;
      (*nod_2D)[5] = 8;
      break;
    case 3:
      memcpy(*eta_3D,*ksi_2D,(*n_ip)*sizeof(double));
      memcpy(*zet_3D,*eta_2D,(*n_ip)*sizeof(double));
      /* ksi_3D = 0 */
      (*nod_2D)[0] = 0;
      (*nod_2D)[1] = 3;
      (*nod_2D)[2] = 2;
      (*nod_2D)[3] = 7;
      (*nod_2D)[4] = 9;
      (*nod_2D)[5] = 6;
      break;
    }/* end switch face */
    break;
  case 8: // Hex --> quad
    err = get_quad_quadrature_rule(int_order,n_ip,ksi_2D,eta_2D,wt_2D);
    *nne_2D = 4;
    *ksi_3D = PGFEM_calloc(double, *n_ip);
    *eta_3D = PGFEM_calloc(double, *n_ip);
    *zet_3D = PGFEM_calloc(double, *n_ip);
    *nod_2D = PGFEM_calloc(int, *nne_2D);    
    switch(face_id)
    {
      case 0:
        memcpy(*eta_3D,*ksi_2D,(*n_ip)*sizeof(double));
        memcpy(*ksi_3D,*eta_2D,(*n_ip)*sizeof(double));
        (*nod_2D)[0] = 3; // 4
        (*nod_2D)[1] = 2; // 3
        (*nod_2D)[2] = 1; // 2
        (*nod_2D)[3] = 0; // 1
        break;
      case 1:
        memcpy(*ksi_3D,*ksi_2D,(*n_ip)*sizeof(double));
        memcpy(*eta_3D,*eta_2D,(*n_ip)*sizeof(double));
        (*nod_2D)[0] = 7; // 8
        (*nod_2D)[1] = 4; // 5
        (*nod_2D)[2] = 5; // 6
        (*nod_2D)[3] = 6; // 7
        break;
      case 2:
        memcpy(*eta_3D,*ksi_2D,(*n_ip)*sizeof(double));
        memcpy(*zet_3D,*eta_2D,(*n_ip)*sizeof(double));
        (*nod_2D)[0] = 0; // 1
        (*nod_2D)[1] = 1; // 2
        (*nod_2D)[2] = 5; // 6
        (*nod_2D)[3] = 4; // 5
        break;
      case 3:
        memcpy(*zet_3D,*ksi_2D,(*n_ip)*sizeof(double));
        memcpy(*ksi_3D,*eta_2D,(*n_ip)*sizeof(double));
        (*nod_2D)[0] = 2; // 3
        (*nod_2D)[1] = 6; // 7
        (*nod_2D)[2] = 5; // 6
        (*nod_2D)[3] = 1; // 2
        break;
       case 4:
        memcpy(*zet_3D,*ksi_2D,(*n_ip)*sizeof(double));
        memcpy(*eta_3D,*eta_2D,(*n_ip)*sizeof(double));
        (*nod_2D)[0] = 3; // 4
        (*nod_2D)[1] = 7; // 8
        (*nod_2D)[2] = 6; // 7
        (*nod_2D)[3] = 2; // 3
        break;
       case 5:
        memcpy(*ksi_3D,*ksi_2D,(*n_ip)*sizeof(double));
        memcpy(*zet_3D,*eta_2D,(*n_ip)*sizeof(double));
        (*nod_2D)[0] = 3; // 4
        (*nod_2D)[1] = 0; // 1
        (*nod_2D)[2] = 4; // 5
        (*nod_2D)[3] = 7; // 8
        break;
      }
      break;
    default:
    {
      int myrank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
      PGFEM_printerr("[%d] WARNING: unrecognized element type!"
          " Surface will not be integrated! %s:%s:%d\n",
          myrank,__func__,__FILE__,__LINE__);
      break;
    }
  }/* end switch type (nne) */

  return err;
}

double compute_surface_jacobian(const int nne_2D,
                const int *nod_2D,
                const double *x,
                const double *y,
                const double *z,
                const double ksi_2D,
                const double eta_2D)
{
  double J = 0.0;
  double *Nk = PGFEM_calloc(double, nne_2D);
  double *Ne = PGFEM_calloc(double, nne_2D);

  /* compute derivative of 2D shape functions with respect to natural
     coordinates */
  dN_2D(nne_2D,ksi_2D,eta_2D,Nk,Ne);

  double *jac = PGFEM_calloc(double, 6);
  double *jtj = PGFEM_calloc(double, 4);

  /* compute the Jacobian matrix */
  for(int i=0; i<nne_2D; i++){
    jac[0] += Nk[i]*x[nod_2D[i]];
    jac[1] += Ne[i]*x[nod_2D[i]];
    jac[2] += Nk[i]*y[nod_2D[i]];
    jac[3] += Ne[i]*y[nod_2D[i]];
    jac[4] += Nk[i]*z[nod_2D[i]];
    jac[5] += Ne[i]*z[nod_2D[i]];
  }

  /* Compute the Jacobian of the transformation */
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,2,2,3,
          1.0,jac,2,jac,2,0.0,jtj,2);
  J = sqrt(det2x2(jtj));

  free(Nk);
  free(Ne);
  free(jac);
  free(jtj);

  return J;
}
