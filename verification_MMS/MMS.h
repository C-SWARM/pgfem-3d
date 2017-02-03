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

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */


void MMS_displacement(double *u, double t, double X, double Y, double Z);
void MMS_velocity(double *v, double t, double X, double Y, double Z);
void MMS_initial_velocity(double *v, double X, double Y, double Z);
double MMS_DC(double t, double X, double Y, double Z, double u_t);
void MMS_pressure_volume(double *P, double *V, HOMMAT const * hommat, double t, double X, double Y, double Z);
double MMS_S1_1(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_S1_2(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_S1_3(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_S2_1(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_S2_2(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_S2_3(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_S3_1(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_S3_2(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_S3_3(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dS1_1(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dS1_2(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dS1_3(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dS2_1(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dS2_2(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dS2_3(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dS3_1(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dS3_2(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dS3_3(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dSdDC1_1(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dSdDC1_2(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dSdDC1_3(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dSdDC2_1(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dSdDC2_2(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dSdDC2_3(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dSdDC3_1(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dSdDC3_2(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);
double MMS_dSdDC3_3(double m10, double m01, double kappa, double mu, double X, double Y, double Z, double u_t, double DC);

void MMS_body_force(double *b, HOMMAT const * hommat, double t, double X, double Y, double Z);


#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
