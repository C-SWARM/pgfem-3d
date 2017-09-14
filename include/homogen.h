/*******************************************
 * HOMOGENIZATION of Composite media    3D *
 * Karel Matous                            *
 * December 2000                       (c) *
 *******************************************/
#ifndef PGFEM3D_HOMOGEN_H
#define PGFEM3D_HOMOGEN_H

#include "PGFem3D_options.h"
#include "data_structure.h"
#include "element.h"
#include "hommat.h"
#include "material.h"
#include "matgeom.h"
#include "tfa.h"

void Mat_3D_orthotropic (const long nmat,
                         Material *mater,
                         const int analysis);

void Mat_trans_isotropic (long nmat,
                          Material *mater);

void hom_matrices (long ***a,
                   const long ne,
                   const long nmat,
                   const long nc,
                   Element *elem,
                   Material *mater,
                   MATGEOM matgeom,
                   HOMMAT *hommat,
                   const long SHAPE,
                   const int analysis);

/** Assign the homogenized material to the list of elements */
void assign_homogenized_material(const long ne,
                                 Element *elem,
                                 long ***a,
                                 const int analysis);

void funkce_Wf (long jj,
                Material *mater,
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
              Material *mater,
              MATGEOM matgeom,
              double **Bf,
              double **Bm,
              double **Mf,
              double **Mm);

void mori_tanaka (long ii,
                  long jj,
                  long kk,
                  long nn,
                  Material *mater,
                  MATGEOM matgeom,
                  HOMMAT *hommat);

void Stiffness_Matrix_3D (long ii,
                          long ip,
                          Element *elem,
                          HOMMAT *hommat,
                          double **D,
                          long TYPE);

TFA* build_tfa (long i);

void Overall_Mat (long ii,
                  long jj,
                  long kk,
                  long nn,
                  Material *mater,
                  MATGEOM matgeom,
                  HOMMAT *hommat,
                  long TYPE);

void TFA_tensors (long ***a,
                  long ne,
                  long nmat,
                  long nc,
                  Element *elem,
                  Material *mater,
                  MATGEOM matgeom,
                  HOMMAT *hommat,
                  TFA *tfa,
                  long SHAPE,
                  long TYPE);

void A_D_tensors (long ii,
                  long jj,
                  long kk,
                  long nn,
                  Material *mater,
                  MATGEOM matgeom,
                  HOMMAT *hommat,
                  TFA *tfa,
                  long SHAPE,
                  long TYPE);

#endif // PGFEM3D_HOMOGEN_H
