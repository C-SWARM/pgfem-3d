#ifndef TENSORS_H
#define TENSORS_H

#include "data_structure.h"
#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef EPS_H
#include "eps.h"
#endif

/** */
  void shape_tensor (long nne,
             long ndofn,
             double *N_x,
             double *N_y,
             double *N_z,
             double ****AA);

  /** Function convert D from matrix to tensor form : D(6,
      6) -> Dijkl. */
  void matrix_tensor_3D (long mat,
             HOMMAT *hommat,
             double A[3][3][3][3]);

  /** Fn -> Fn-BAR || LOOK OUT. */
  void tensors_FF_f (long nne,
             long ndn,
             double **Fn,
             double **Fr,
             double **Fr_I,
             double Jr,
             double ****ST,
             double ****FF,
             double **f);

  /** Fn -> Fn-BAR || LOOK OUT. */
  void tensors_FFla (long ii,
             long ip,
             long nne,
             long ndn,
             EPS *eps,
             double **Fn,
             double **Fr,
             double **Fr_I,
             double Jr,
             double **f,
             double **UU,
             double **Flam,
             double **FFlam,
             double **pGplam,
             const int analysis);

  /** Fn -> Fn-BAR || LOOK OUT. */
  void tensors_aa_bb_dd_mm (long nne,
                long ndn,
                long npres,
                double ai,
                double aj,
                double ak,
                double J,
                double *Psi,
                double L[3][3][3][3],
                double **Fn,
                double **Fr,
                double **Fr_I,
                double Jn,
                double Jr,
                double Tn,
                double Tr,
                double ****ST,
                double ****FF,
                double **f,
                double ***AA,
                double ***aa,
                double ***BB,
                double ***bb,
                double **DD,
                double **dd,
                double **MM,
                double **mm);

  /** */
  void tensors_GG (long ii,
           long ip,
           long nne,
           long ndn,
           double ****ST,
           double ****FF,
           double **f,
           double **UU,
           EPS *eps,
           double ****GT,
           double gT[3][3],
           double ****gg);

  /** */
  void tensors_aa_bb_dd_mm_plast (long ii,
                  long ip,
                  long nne,
                  long ndn,
                  long npres,
                  double ai,
                  double aj,
                  double ak,
                  double J,
                  double *Psi,
                  double L[3][3][3][3],
                  double **Fn,
                  double **Fr,
                  double **Fr_I,
                  double Jn,
                  double Jr,
                  double Tn,
                  double Tr,
                  double ****ST,
                  double ****FF,
                  double **f,
                  double ***AA,
                  double ***aa,
                  double ***BB,
                  double ***bb,
                  double **DD,
                  double **dd,
                  double **MM,
                  double **mm,
                  double **S,
                  double ****pFFp,
                  double **pfp,
                  double **UU,
                  EPS *eps,
                  double ***BB1,
                  double ***bb1,
                  double ***BB2,
                  double ***bb2,
                  double **MMp,
                  double **mmp,
                  double ****pGp);

#endif /* #ifndef TENSORS_H */
