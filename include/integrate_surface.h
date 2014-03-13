/* HEADER */
#pragma once
#ifndef INTEGRATE_SURFACE_H
#define INTEGRATE_SURFACE_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

 /** get the number of integration points and integration point
    locations for the surface integral. NOTE: ksi, eta and zet are
    allocated within the function. Memory must be free'd by calling
    function */
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
		      int **nod_2D);

  /** compute the Jacobian of the transformation for the surface
      integration */
  double compute_surface_jacobian(const int nne_2D,
				  const int *nod_2D,
				  const double *x,
				  const double *y,
				  const double *z,
				  const double ksi_2D,
				  const double eta_2D);
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
