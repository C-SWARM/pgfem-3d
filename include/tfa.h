/* HEADER */
#ifndef TFA_H
#define TFA_H

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


#endif /* #ifndef  */
