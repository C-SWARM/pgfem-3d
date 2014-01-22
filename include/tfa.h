/* HEADER */
#ifndef TFA_H
#define TFA_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure of material properties */
  struct TFA{
    double *Dmm,
      *Dmd,
      *Dmb,
      *Dbm,
      *Dbd,
      *Dbb,
      *Ddm,
      *Ddd,
      *Ddb,
      *Am,
      *A2;
  };
  typedef struct TFA TFA;

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
