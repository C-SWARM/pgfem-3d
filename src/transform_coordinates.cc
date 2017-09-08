#include "transform_coordinates.h"
#include <string.h>
#include <stdlib.h>
#include "mkl_cblas.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

/* TESTED MM 1/21/2013 */
void transform_coordinates(const int n_pts,
               const double *x,
               const double *y,
               const double *z,
               const double *e1,
               const double *e2,
               const double *n,
               const int TRANS, /* boolean to use transpose */
               double *tx,
               double *ty,
               double *tz)
{
  /*
    TRANS = 0 ==> Global to local coordinates
    TRANS = 1 ==> Local to global coordinates
   */
  double *TM = PGFEM_calloc(double, 9);
  transformation_matrix(e1,e2,n,TM);

  double *result_mat = PGFEM_calloc(double, 3*n_pts);
  double *coords_mat = PGFEM_calloc(double, 3*n_pts);
  /* coords_mat = { x , y , z }' */
  memcpy(coords_mat,             x,n_pts*sizeof(double));
  memcpy(coords_mat + n_pts,     y,n_pts*sizeof(double));
  memcpy(coords_mat + (2*n_pts), z,n_pts*sizeof(double));

  if(TRANS){
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
        3,n_pts,3,1.0,TM,3,coords_mat,n_pts,
        0.0,result_mat,n_pts);
  } else {
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
        3,n_pts,3,1.0,TM,3,coords_mat,n_pts,
        0.0,result_mat,n_pts);
  }

  /* result_mat = { tx , ty , tz }' */
  memcpy(tx, result_mat,             n_pts*sizeof(double));
  memcpy(ty, result_mat + n_pts,     n_pts*sizeof(double));
  memcpy(tz, result_mat + (2*n_pts), n_pts*sizeof(double));

  free(TM);
  free(coords_mat);
  free(result_mat);
}/* transform_coordinates() */

/* TESTED MM 1/21/2013 */
void transformation_matrix(const double *e1,
               const double *e2,
               const double *n,
               double *tm)
{
  /* T = { e1 , e2 , n }' */
  memcpy(tm,e1,3*sizeof(double));
  memcpy(tm+3,e2,3*sizeof(double));
  memcpy(tm+6,n,3*sizeof(double));
}/* transformation_matrix() */


/*=============================================================================
  DRIVER PROGRAM
  compile command:

  icc -Wall -O0 -std=c99 -I../include \
  -I/opt/crc/intel/12.0/mkl/include -DTRANSFORM_COORDINATES_DRIVER -o \
  tc_driver transform_coordinates.c -Wl,--start-group \
  /opt/crc/intel/12.0/mkl/lib/intel64/libmkl_intel_lp64.a \
  /opt/crc/intel/12.0/mkl/lib/intel64/libmkl_sequential.a \
  /opt/crc/intel/12.0/mkl/lib/intel64/libmkl_core.a \
  /opt/crc/intel/12.0/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
  -Wl,--end-group -lpthread -limf -lm
  =============================================================================*/
#ifdef TRANSFORM_COORDINATES_DRIVER
#include <math.h>
#include <stdio.h>
static void cross_product(const double *a,
              const double *b,
              double *c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

static void print_array(const int len,
            const double *a)
{
  for(int i=0; i<len; i++){
    printf("%f ",a[i]);
  }
  printf("\n");
}

int main()
{
  constexpr int n_pts = 3; /* here we MUST use three to guarantee plane */
  double *x = PGFEM_calloc(double, n_pts);
  double *y = PGFEM_calloc(double, n_pts);
  double *z = PGFEM_calloc(double, n_pts);
  double *tx = PGFEM_calloc(double, n_pts);
  double *ty = PGFEM_calloc(double, n_pts);
  double *tz = PGFEM_calloc(double, n_pts);

  double *coords = PGFEM_calloc(double, n_pts*3);
  double *e1 = PGFEM_calloc(double, n_pts);
  double *e2 = PGFEM_calloc(double, n_pts);
  double *n = PGFEM_calloc(double, n_pts);

  /* get random coordinates for the three points */
  for(int i=0; i<n_pts; i++){
    x[i] = ((double) rand())/((double) (rand()+1));
    y[i] = ((double) rand())/((double) (rand()+1));
    z[i] = ((double) rand())/((double) (rand()+1));
    coords[i*3 + 0] = x[i];
    coords[i*3 + 1] = y[i];
    coords[i*3 + 2] = z[i];
  }

  /* compute basis */
  /* e1 = c1 - c0 || e2 = c2 - c0  */
  memcpy(e1,coords + 3,3*sizeof(double));
  memcpy(e2,coords + 6,3*sizeof(double));

  cblas_daxpy(3,-1.0,coords,1,e1,1);
  cblas_daxpy(3,-1.0,coords,1,e2,1);

  cblas_dscal(3,1./cblas_dnrm2(3,e1,1),e1,1);
  cblas_dscal(3,1./cblas_dnrm2(3,e2,1),e2,1);

  cross_product(e1,e2,n);
  cblas_dscal(3,1./cblas_dnrm2(3,n,1),n,1);

  cross_product(e1,n,e2);
  cblas_dscal(3,1./cblas_dnrm2(3,e2,1),e2,1);

  /* print coords */
  printf("x || y || z\n");
  print_array(n_pts,x);
  print_array(n_pts,y);
  print_array(n_pts,z);

  /* print basis */
  printf("e1 || e2 || n\n");
  print_array(3,e1);
  print_array(3,e2);
  print_array(3,n);

  transform_coordinates(n_pts,x,y,z,e1,e2,n,0,tx,ty,tz);

  printf("tx || ty || tz\n");
  print_array(n_pts,tx);
  print_array(n_pts,ty);
  print_array(n_pts,tz);

  memset(x,0,n_pts*sizeof(double));
  memset(y,0,n_pts*sizeof(double));
  memset(z,0,n_pts*sizeof(double));

  transform_coordinates(n_pts,tx,ty,tz,e1,e2,n,1,x,y,z);

  printf("x || y || z\n");
  print_array(n_pts,x);
  print_array(n_pts,y);
  print_array(n_pts,z);

  free(x);
  free(y);
  free(z);
  free(tx);
  free(ty);
  free(tz);
  free(coords);
  free(e1);
  free(e2);
  free(n);
}
#endif
