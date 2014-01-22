/* HEADER */
#ifndef MATERIAL_H
#define MATERIAL_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure of material properties */
  struct MATERIAL{
    double Ex,Ey,Ez,Gyz,Gxz,Gxy,nyz,nxz,nxy,ax,ay,az,sig;
    /*Elastic stiffness*/
    double L[9];/* Local coordinate system */
    double M[9];/* Local coordinate system */

    /* potential function flags */
    int devPotFlag;
    int volPotFlag;
  };
  typedef struct MATERIAL MATERIAL;

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
