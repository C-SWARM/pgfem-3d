/* HEADER */
#ifndef MATGEOM_H
#define MATGEOM_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /**  Structures of geometric properties */
  struct MATGEOM_1{
    /* Volume fraction (cf), Fibre orientation [vector base ee =
       (e1,e2,e3)] */
    double *cf,*cm,*cd,a1,a2,**ee,*H2;
    long SH,*H1;
  };
  typedef struct MATGEOM_1 MATGEOM_1;
  typedef struct MATGEOM_1 *MATGEOM;

  MATGEOM build_matgeom(long nc,long np);
  void destroy_matgeom(MATGEOM mg, long np);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
