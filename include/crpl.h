/* HEADER */
#ifndef CRPL_H
#define CRPL_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure of CRYSTAL PLASTICITY material properties */
  struct CRPL {
    long nss;
    double a,m,TH,To,Tso,Go,mm,**P,Tn;
    /* PLC */
    double b,fo,nu,l1,l2,c1,c2,c3,to;
  };
  typedef struct CRPL CRPL;

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
