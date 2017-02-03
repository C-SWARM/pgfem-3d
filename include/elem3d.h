/* HEADER */

/**
 * @file @todo Describe remaining functions and migrate integration
 * rules to @see quadrature_rules.c
 */

#pragma once
#ifndef ELEM3D_H
#define ELEM3D_H

#include "element.h"
#include "node.h"
#include "hommat.h"
#include "PGFem3D_options.h"


#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

/**
 * Get the number of integration points for the element.
 */
void int_point (const long nne,
		long *II);

/** 
 * Get the number of integration points in each direction, parent
 * coordinates and weights for integration.
 */
void integrate (const long nne,
		long *II,
		long *JJ,
		long *KK,
		double *gk,
		double *ge,
		double *gz,
		double *w);

void intpoints_2 (double *gs,
		  double *w);

void intpoints_3 (double *gs,
		  double *w);

void int_tetra_1 (double *gk,
		  double *ge,
		  double *gz,
		  double *w);

void int_tetra_4 (double *gk,
		  double *ge,
		  double *gz,
		  double *w);

void int_tetra_5 (double *gk,
		  double *ge,
		  double *gz,
		  double *w);

void int_tetra_11 (double *gk,
		   double *ge,
		   double *gz,
		   double *w);

void int_tetra_24(double *gk,
		  double *ge,
		  double *gz,
		  double *w);

void shape_tria (const double ksi,
		 const double eta,
		 const double zet,
		 double *N);

void shape_func(const double ksi,
		const double eta,
		const double zet,
		const long nne,
		double *N);

void shape_2D (const long nne,
	       const double ksi,
	       const double eta,
	       double *N);

void dN_kez(const double ksi,
	    const double eta,
	    const double zet,
	    const long nne,
	    double *N_ksi,
	    double *N_eta,
	    double *N_zet);

void dxyz_kez (const double ksi,
	       const double eta,
	       const double zet,
	       const long nne,
	       const double *x,
	       const double *y,
	       const double *z,
	       const double *N_ksi,
	       const double *N_eta,
	       const double *N_zet,
	       double *dx,
	       double *dy,
	       double *dz);

double Jacobi(const double ksi,
	      const double eta,
	      const double zet,
	      const double *x,
	      const double *y,
	      const double *z,
	      const double *dx,
	      const double *dy,
	      const double *dz);

double deriv(const double ksi,
	     const double eta,
	     const double zet,
	     const long nne,
	     const double *x,
	     const double *y,
	     const double *z,
	     double *N_x,
	     double *N_y,
	     double *N_z);

/** Get the element nodal coordinates in the parent frame. Useful for
    integrating over the boundary of a volume element */
void get_element_node_parent_coords(const int nne,
				    double *ksi,
				    double *eta,
				    double *zet);

/** Compute the center of an element and put the coordinates as the
    nne+1 element of x,y,z */
void element_center(const int nne,
		    double *x,
		    double *y,
		    double *z);

/** Compute the (x,y,z) coordinates from the center of the element in
    the natural coordinates.  Put results in x_i[nne] */
int element_center_kez(const int nne,
		       double *x,
		       double *y,
		       double *z);

/** Get the gradient of the bubble function */
void get_bubble_grad(const int nne_t, /* The bubble is the last node */
		     const double ksi,
		     const double eta,
		     const double zet,
		     const double *x,
		     const double *y,
		     const double *z,
		     double *N_x,
		     double *N_y,
		     double *N_z);

double Bmat(const double ksi,
	    const double eta,
	    const double zet,
	    const long nne,
	    const double *x,
	    const double *y,
	    const double *z,
	    double **B);

void B_BAR (double **B,
	    const long nne,
	    const double *x,
	    const double *y,
	    const double *z);

void stiffmat (long *adr,
	       long ne,
	       long ndofn,
	       ELEMENT *elem,
	       NODE *node,
	       HOMMAT *hommat,
	       long *ci,
	       long typsolveru,
	       double *k,
	       const PGFem3D_opt *opts,
	       const int mp_id);

void stiffmatel (long ii,
		 double *x,
		 double *y,
		 double *z,
		 long nne,
		 long ndofn,
		 ELEMENT *elem,
		 HOMMAT *hommat,
		 NODE *node,
		 double *K,
		 const PGFem3D_opt *opts);


#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
