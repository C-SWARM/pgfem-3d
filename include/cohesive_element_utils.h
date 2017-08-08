/* HEADER */
/**
 * AUTHORS:
 * Matt Mosby, University of Notre Dame, mmosby1 [at] nd.edu
 */
#pragma once
#ifndef COHESIVE_ELEMENT_UTILS_H
#define COHESIVE_ELEMENT_UTILS_H

  long int_pointC (long nne);

  void integrateC (long nne,
           long *II,
           long *JJ,
           double *gk,
           double *ge,
           double *w);

  void shape_2DC (long nne,
          double ksi,
          double eta,
          double *N);

  void dN_2D (long nne,
          double ksi,
          double eta,
          double *N_ksi,
          double *N_eta);

  void dN_dxb (long nne,
           double ksi,
           double eta,
           double *e1,
           double *e2,
           double *n,
           double ***Nxb);

  double dN3_xy (double ksi,
         double eta,
         long nne,
         double *x,
         double *y,
         double *z,
         double *N_x,
         double *N_y);

  void mean_map (long nne,
         double *x,
         double *y,
         double *z,
         double *r_u,
         double *xb,
         double *yb,
         double *zb);

  void base_vec (const long nne,
         const double ksi,
         const double eta,
         const double *x,
         const double *y,
         const double *z,
         double *e1,
         double *e2,
         double *e2h,
         double *n,
         const int myrank);

  void tran_coord (long nne,
           double *x,
           double *y,
           double *z,
           double *e1,
           double *e2,
           double *n,
           double *xl,
           double *yl,
           double *zl,
           long TO);

  int get_jump(const long nne,
           const double *x,
           const double *y,
           const double *z,
           const double *r_u,
           const double *N,
           double *jump);

  int get_jump_nt(double *jump_n,
          double *jump_t,
          const double *jump,
          const double *normal);

  int get_eff_jump(double *eff_jump,
           const double jump_n,
           const double jump_t,
           const double beta);

  double get_Xxi (const double bet,
          const double *Xi,
          const double *n);

#endif
