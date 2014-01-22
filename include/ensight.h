/* HEADER */
#ifndef ENSIGHT_H
#define ENSIGHT_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure of ENSIGHT */
  struct ENSIGHT_1{
    long ncn,*Sm,*Sp,*No,NVp,NCp,*Vp,*Cp;
  };
  typedef struct ENSIGHT_1 ENSIGHT_1;
  typedef struct ENSIGHT_1 *ENSIGHT;

  void destroy_ensight(ENSIGHT ensight);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
