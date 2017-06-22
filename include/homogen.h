/*******************************************
 * HOMOGENIZATION of Composite media    3D *
 * Karel Matous                            *
 * December 2000                       (c) *
 *******************************************/

#ifndef HOMOGEN_H
#define HOMOGEN_H

#include "data_structure.h"

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef MATERIAL_H
#include "material.h"
#endif

#ifndef MATGEOM_H
#include "matgeom.h"
#endif

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef TFA_H
#include "tfa.h"
#endif

#ifndef PGFEM_OPTIONS_H
#include "PGFem3D_options.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

void Mat_3D_orthotropic (const long nmat,
			 MATERIAL *mater,
			 const int analysis);

void Mat_trans_isotropic (long nmat,
			  MATERIAL *mater);

void hom_matrices (long ***a,
		   const long ne,
		   const long nmat,
		   const long nc,
		   ELEMENT *elem,
		   MATERIAL *mater,
		   MATGEOM matgeom,
		   HOMMAT *hommat,
		   const long SHAPE,
		   const int analysis);

/** Assign the homogenized material to the list of elements */
void assign_homogenized_material(const long ne,
				 ELEMENT *elem,
				 long ***a,
				 const int analysis);

void funkce_Wf (long jj,
		MATERIAL *mater,
		double **Lf,
		double **Mf,
		double **Lm,
		double **Mm,
		double **Wf);

void funkce_Bf (long kk,
		MATGEOM matgeom,
		double **Wf,
		double **Bf);

void funkce_Bm (long kk,
		MATGEOM matgeom,
		double **Wf,
		double **Bm);

void mat_vrs (long nn,
	      long kk,
	      HOMMAT *hommat,
	      MATERIAL *mater,
	      MATGEOM matgeom,
	      double **Bf,
	      double **Bm,
	      double **Mf,
	      double **Mm);

void mori_tanaka (long ii,
		  long jj,
		  long kk,
		  long nn,
		  MATERIAL *mater,
		  MATGEOM matgeom,
		  HOMMAT *hommat);

void Stiffness_Matrix_3D (long ii,
			  long ip,
			  ELEMENT *elem,
			  HOMMAT *hommat,
			  double **D,
			  long TYPE);

TFA* build_tfa (long i);

void Overall_Mat (long ii,
		  long jj,
		  long kk,
		  long nn,
		  MATERIAL *mater,
		  MATGEOM matgeom,
		  HOMMAT *hommat,
		  long TYPE);

void TFA_tensors (long ***a,
		  long ne,
		  long nmat,
		  long nc,
		  ELEMENT *elem,
		  MATERIAL *mater,
		  MATGEOM matgeom,
		  HOMMAT *hommat,
		  TFA *tfa,
		  long SHAPE,
		  long TYPE);

void A_D_tensors (long ii,
		  long jj,
		  long kk,
		  long nn,
		  MATERIAL *mater,
		  MATGEOM matgeom,
		  HOMMAT *hommat,
		  TFA *tfa,
		  long SHAPE,
		  long TYPE);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
