#ifndef __MMS_MMS_h__
#define __MMS_MMS_h__

#include<math.h>
#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef ELEM3D_H
#include "elem3d.h"
#endif

void MMS_displacement(double *u, double t, double X, double Y, double Z);
void MMS_initial_velocity(double *v, double X, double Y, double Z);
void MMS_pressure_volume(double *P, double *V, HOMMAT const * hommat, double t, double X, double Y, double Z);

void MMS_body_force(double *b, HOMMAT const * hommat, double t, double X, double Y, double Z);

#endif
