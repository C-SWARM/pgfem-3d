/* HEADER */
#ifndef CRPL_H
#define CRPL_H

/** Structure of CRYSTAL PLASTICITY material properties */
struct CRPL {
  long nss;
  double a,m,TH,To,Tso,Go,mm,**P,Tn;
  /* PLC */
  double b,fo,nu,l1,l2,c1,c2,c3,to;
};

#endif /* #ifndef  */
