#include "allocation.h"
#include "utils.h"
#include "MMS.h"
#include<math.h>
#include "hommat.h"

#include "cm_MMS.h"
#include "hyperelasticity.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef ELEM3D_H
#include "elem3d.h"
#endif


//void MMS4cm_pressure_volume(double *P, double *V, ELASTICITY *elast, double t, double X, double Y, double Z);
//void MMS4cm_body_force(double *b, HOMMAT const * hommat, ELASTICITY *elast, double t, double X, double Y, double Z);

void MMS_HE_displacement(double *u, double t, double X, double Y, double Z);
void MMS_HE_initial_velocity(double *v, double X, double Y, double Z);
void MMS_HE_pressure_volume(double *P, double *V, HOMMAT const * hommat, double t, double X, double Y, double Z);

void MMS_HE_body_force(double *b, HOMMAT const * hommat, double t, double X, double Y, double Z);


void MMS_displacement(double *u, double t, double X, double Y, double Z, const bool is4cm){
  if(is4cm)
    MMS4cm_displacement(u, t, X, Y, Z);
  else
    MMS_HE_displacement(u, t, X, Y, Z);
}

void MMS_initial_velocity(double *v, double X, double Y, double Z, const bool is4cm){
  if(is4cm)
    MMS4cm_initial_velocity(v, X, Y, Z);
  else
    MMS_HE_initial_velocity(v, X, Y, Z);
  
}
void MMS_pressure_volume(double *P, double *V, const HOMMAT *hommat, ELASTICITY *elast, double t, double X, double Y, double Z, const bool is4cm){
  if(is4cm)
    MMS4cm_pressure_volume(P, V, elast, t, X, Y, Z);
  else
    MMS_HE_pressure_volume(P, V, hommat, t, X, Y, Z);
}

void MMS_body_force(double *b, const HOMMAT *hommat, ELASTICITY *elast, double t, double X, double Y, double Z, const bool is4cm){
  if(is4cm)
    MMS4cm_body_force(b, hommat, elast, t, X, Y, Z);
  else
    MMS_HE_body_force(b, hommat, t, X, Y, Z);
}
