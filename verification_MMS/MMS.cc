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

static constexpr double MMS4cm_c[6] = {10.0, 0.5, 1000.0, 1.0, 0.7, -50.0};
static constexpr double MMS_HE_c[6] = {10.0, 0.5, 1000.0, 1.0, 0.0, -50.0};

void MMS_displacement(double *u, double t, double X, double Y, double Z, const bool is4cm){
  if(t<0.0)
    u[0] = u[1] = u[2] = 0.0;
    
  if(is4cm)
    MMS4cm_displacement(u, t, X, Y, Z, MMS4cm_c);
  else
    MMS4cm_displacement(u, t, X, Y, Z, MMS_HE_c);
}

void MMS_initial_velocity(double *v, double X, double Y, double Z, const bool is4cm){
  if(is4cm)
    MMS4cm_initial_velocity(v, X, Y, Z, MMS4cm_c);
  else
    MMS4cm_initial_velocity(v, X, Y, Z, MMS_HE_c);
  
}
void MMS_pressure_volume(double *P, double *V, ELASTICITY *elast, double t, double X, double Y, double Z, const bool is4cm){
  if(is4cm)
    MMS4cm_pressure_volume(P, V, elast, t, X, Y, Z, MMS4cm_c);
  else
    MMS4cm_pressure_volume(P, V, elast, t, X, Y, Z, MMS_HE_c);
}

void MMS_body_force(double *b, const HOMMAT *hommat, ELASTICITY *elast, double t, double X, double Y, double Z, const bool is4cm){
  if(t<0.0)
    b[0] = b[1] = b[2] = 0.0;
  if(is4cm)
    MMS4cm_body_force(b, hommat, elast, t, X, Y, Z, MMS4cm_c);
  else
    MMS4cm_body_force(b, hommat, elast, t, X, Y, Z, MMS_HE_c);
}
