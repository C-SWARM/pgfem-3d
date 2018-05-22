/* HEADER */
#include "pgfem3d/Communication.hpp"
#include "quadrature_rules.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

using namespace pgfem3d;

/* REFERENCES:
   [1] T.J.R. Hughes, "The Finite Element Method", Dover, 2000
*/

/*=========== TRIANGLE QUADRATURE RULES ==================*/
static int tria_1(int *n_ip,
          double **ksi,
          double **eta,
          double **wt)
{
  /* integration at centroid */

  int err = 0;
  *n_ip = 1;
  *ksi = PGFEM_calloc(double, *n_ip);
  *eta = PGFEM_calloc(double, *n_ip);
  *wt = PGFEM_calloc(double, *n_ip);

  (*wt)[0] = 0.5;
  (*ksi)[0] = 1./3.;
  (*eta)[0] = 1./3.;
  return err;
} /*  tria_1() */

static int tria_2(int *n_ip,
          double **ksi,
          double **eta,
          double **wt)
{
  /* 3-pt formula from [1] Table 3.1.1 pg 173 */
  int err = 0;
  *n_ip = 3;
  *ksi = PGFEM_calloc(double, *n_ip);
  *eta = PGFEM_calloc(double, *n_ip);
  *wt = PGFEM_calloc(double, *n_ip);

  (*wt)[0] = (*wt)[1] = (*wt)[2] = 1./6.;

  (*ksi)[0] = 2./3.; (*eta)[0] = 1./6.;
  (*ksi)[1] = 1./6.; (*eta)[1] = 1./6.;
  (*ksi)[2] = 1./6.; (*eta)[2] = 2./3.;

  return err;
} /*  tria_2() */

static int tria_3(int *n_ip,
          double **ksi,
          double **eta,
          double **wt)
{
  /* 4-pt formula from [1] Table 3.1.1 pg 173 */
  int err = 0;
  *n_ip = 4;
  *ksi = PGFEM_calloc(double, *n_ip);
  *eta = PGFEM_calloc(double, *n_ip);
  *wt = PGFEM_calloc(double, *n_ip);

  (*wt)[0] = -9./32.;
  (*ksi)[0] = (*eta)[0] = 1./3.;

  (*wt)[1] = (*wt)[2] = (*wt)[3] = 25./96.;
  (*ksi)[1] = 0.6; (*eta)[1] = 0.2;
  (*ksi)[2] = 0.2; (*eta)[2] = 0.2;
  (*ksi)[3] = 0.2; (*eta)[3] = 0.6;

  return err;
} /*  tria_3() */

/*=========== QUADRILATERAL QUADRATURE RULES =============*/

/*=========== TETRAHEDRON QUADRATURE RULES ===============*/

/*=========== HEXAHEDRON QUADRATURE RULES ================*/

/*=========== WEDGE QUADRATURE RULES =====================*/

/*=========== PYRAMID QUADRATURE RULES ===================*/


/*=========== API FUNCTIONS ==============================*/
int get_tria_quadrature_rule(const int int_order,
                 int *n_ip,
                 double **ksi,
                 double **eta,
                 double **wt)
{
  int err = 0;
  switch(int_order){
  case 0:
  case 1:
    err = tria_1(n_ip,ksi,eta,wt);
    break;
  case 2:
    err = tria_2(n_ip,ksi,eta,wt);
    break;
  case 3:
    err = tria_3(n_ip,ksi,eta,wt);
    break;
  default:
    {
      int err_rank = 0;
      PGFEM_Error_rank(&err_rank);
      PGFEM_printerr("[%d] WARNING: triangle integration rule for "
             "%d(st/nd/rd/th) order not implemented! %s:%s:%d\n",
             err_rank,int_order,__func__,__FILE__,__LINE__);
      err++;
      break;
    }
  }
  return err;
}

/// 4-pt formula from [1] pg 145
///
/// \param[out] n_in  number of integration point
/// \param[out] ksi   Gaussian quadrature
/// \param[out] eta   Gaussian quadrature
/// \param[out] wt    weight
/// \return non-zero on internal error
static int quad_1(int *n_ip,
                  double **ksi,
                  double **eta,
                  double **wt)
{
  
  int err = 0;
  *n_ip = 4;
  *ksi = PGFEM_calloc(double, *n_ip);
  *eta = PGFEM_calloc(double, *n_ip);
  *wt =  PGFEM_calloc(double, *n_ip);

  double one_over_sqrt_3 = 0.57735026918962584; // 1.0/sqrt(3.0);
  
  (*wt)[0] = (*wt)[1] = (*wt)[2] = (*wt)[3] = 1.0;
  (*ksi)[1] = -one_over_sqrt_3; (*eta)[1] = -one_over_sqrt_3;
  (*ksi)[1] =  one_over_sqrt_3; (*eta)[1] = -one_over_sqrt_3;
  (*ksi)[2] =  one_over_sqrt_3; (*eta)[2] = one_over_sqrt_3;
  (*ksi)[3] = -one_over_sqrt_3; (*eta)[3] = one_over_sqrt_3;

  return err;
}	              
		  
int get_quad_quadrature_rule(const int int_order,
                             int *n_ip,
                             double **ksi,
                             double **eta,
                             double **wt)
{
  int err = 0;
  switch(int_order){
  case 0:
  case 1:
    err = quad_1(n_ip,ksi,eta,wt);
    break;
  case 2:
    err = quad_1(n_ip,ksi,eta,wt);
    break;
  case 3:
    err = quad_1(n_ip,ksi,eta,wt);
    break;
  default:
    {
      int err_rank = 0;
      PGFEM_Error_rank(&err_rank);
      PGFEM_printerr("[%d] WARNING: quadrilateral integration rule for "
		     "%d(st/nd/rd/th) order not implemented! %s:%s:%d\n",
		     err_rank,int_order,__func__,__FILE__,__LINE__);
      err++;
      break;
    }
  }
  return err;
}

int get_tet_quadrature_rule(const int int_order,
                int *n_ip,
                double **ksi,
                double **eta,
                double **zet,
                double **wt)
{
  int err = 0;

  /* function not yet implemented */
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);
  PGFEM_printerr("[%d] ERROR: Function not yet implemented! %s:%s:%d\n",
         err_rank,__func__,__FILE__,__LINE__);
  PGFEM_Abort();
  err++;
  return err;
}

int get_hex_quadrature_rule(const int int_order,
                int *n_ip,
                double **ksi,
                double **eta,
                double **zet,
                double **wt)
{
  int err = 0;

  /* function not yet implemented */
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);
  PGFEM_printerr("[%d] ERROR: Function not yet implemented! %s:%s:%d\n",
         err_rank,__func__,__FILE__,__LINE__);
  PGFEM_Abort();
  err++;
  return err;
}

int get_wedge_quadrature_rule(const int int_order,
                  int *n_ip,
                  double **ksi,
                  double **eta,
                  double **zet,
                  double **wt)
{
  int err = 0;

  /* function not yet implemented */
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);
  PGFEM_printerr("[%d] ERROR: Function not yet implemented! %s:%s:%d\n",
         err_rank,__func__,__FILE__,__LINE__);
  PGFEM_Abort();
  err++;
  return err;
}

int get_pyram_quadrature_rule(const int int_order,
                  int *n_ip,
                  double **ksi,
                  double **eta,
                  double **zet,
                  double **wt)
{
  int err = 0;

  /* function not yet implemented */
  int err_rank = 0;
  PGFEM_Error_rank(&err_rank);
  PGFEM_printerr("[%d] ERROR: Function not yet implemented! %s:%s:%d\n",
         err_rank,__func__,__FILE__,__LINE__);
  PGFEM_Abort();
  err++;
  return err;
}
