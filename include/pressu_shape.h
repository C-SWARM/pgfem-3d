#ifndef PRESSU_SHAPE_H
#define PRESSU_SHAPE_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Pressure shape functions. */
  void pressu_shape (long npres,double ksi,double eta,double zet,double *Psi);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef PRESSU_SHAPE_H */
