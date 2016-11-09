/* HEADER */
#pragma once
#ifndef MATERIAL_H
#define MATERIAL_H

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

#ifndef TYPE_MATERIAL
#define TYPE_MATERIAL
typedef struct MATERIAL MATERIAL;
#endif

typedef struct {
  double k[9]; /// heat conductivity
  double cp;   /// heat capacity
} MATERIAL_THERMAL;

#endif /* #ifndef  */
