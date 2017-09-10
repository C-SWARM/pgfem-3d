/*****************************
 * Program FEM3d ver. 2.0    *
 * FEM - 3D analysis         *
 * Karel Matous              *
 *****************************/

/*=== NO MATCHING SOURCE...? === */

#ifndef PGFEM3D_INELASTIC_H
#define PGFEM3D_INELASTIC_H

#include "element.h"
#include "eps.h"
#include "hommat.h"
#include "material.h"
#include "matgeom.h"
#include "sig.h"
#include "tfa.h"

void inelastic_material (long ne,
                         Element *elem,
                         HOMMAT *hommat,
                         MATGEOM matgeom,
                         MATERIAL *mater,
                         SIG *sig,
                         EPS *eps,
                         TFA *tfa);

void LmaL (long ii,
           long ip,
           long II,
           Element *elem,
           HOMMAT *hommat,
           MATGEOM matgeom,
           MATERIAL *mater,
           SIG *sig,
           EPS *eps,
           TFA *tfa,
           double **LmL);

void LfiL (long ii,
           long ip,
           long II,
           Element *elem,
           HOMMAT *hommat,
           MATGEOM matgeom,
           MATERIAL *mater,
           SIG *sig,
           EPS *eps,
           TFA *tfa,
           double **LmL,
           double **LfL);

double H_modulus (long i,
                  MATGEOM matgeom,
                  MATERIAL *mater,
                  double Eeq);

void plast_normal (long ii,
                   long ip,
                   SIG *sig,
                   double **dsig_m,
                   double **dsig_a,
                   double *m);

double plast_functi (long ii,
                     long ip,
                     Element *elem,
                     MATERIAL *mater,
                     SIG *sig,
                     double **dsig_m,
                     double **dsig_a,
                     long TYPE);

double plast_step (long ii,
                   long ip,
                   Element *elem,
                   MATERIAL *mater,
                   SIG *sig,
                   double **dsig_m,
                   long TYPE);

void plast_dam_material (long ii,
                         long ip,
                         long II,
                         Element *elem,
                         HOMMAT *hommat,
                         MATGEOM matgeom,
                         MATERIAL *mater,
                         SIG *sig,
                         EPS *eps,
                         TFA *tfa);

void sig_e_in (long ii,
               long nne,
               Element *elem,
               HOMMAT *hommat,
               MATGEOM matgeom,
               MATERIAL *mater,
               TFA *tfa,
               SIG *sig,
               double **dsig_m,
               double **dsig_f,
               double **dsig_a,
               EPS *eps,
               double **deps,
               double **deps_i,
               double nor_min);

void eps_phas (long ii,
               long nne,
               Element *elem,
               HOMMAT *hommat,
               MATGEOM matgeom,
               MATERIAL *mater,
               TFA *tfa,
               double **dsig_m,
               double **dsig_f,
               double **deps_m,
               double **deps_f,
               double **deps_i);

void inelastic_increment (double *f_u,
                          long ne,
                          long ndofn,
                          long ndofd,
                          double *d_r,
                          Element *elem,
                          NODE *node,
                          MATERIAL *mater,
                          HOMMAT *hommat,
                          MATGEOM matgeom,
                          SIG *sig,
                          EPS *eps,

                          SUPP sup,
                          double nor_min,
                          TFA *tfa);

void check_load (long ne,
                 Element *elem,
                 SIG *sig);

void increment (long ne,
                long ndofn,
                NODE *node,
                Element *elem,
                SIG *sig,
                EPS *eps);

void set_back (long ne,
               Element *elem,
               SIG *sig,
               EPS *eps);

#endif // PGFEM3D_INELASTIC_H
