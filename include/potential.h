#ifndef POTENTIAL_H
#define POTENTIAL_H

#ifndef HOMMAT_H
#include "hommat.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** These functions contain the calculations of important tensors
      relating to the deviatoric and volumetric energy potential
      functions.

      NOTE: the bulk modulus has been factored out to
      normalize the residual vector and thus allow for the l_2 norm
      convergence test. */

  typedef void (*deviatoricStressFunctionPtr)(const double * const *C,
					      const HOMMAT *material,
					      double **S);

  typedef void (*materialStiffnessFunctionPtr)(const double * const *C,
					       const HOMMAT *material,
					       const double * const *Fn,
					       const double * const *Fr,
					       double ****L);

  typedef void (*volumetricPressureFunctionPtr)(const double Jn,
						const double Jr,
						const HOMMAT *material,
						double *pres);

  typedef void (*linearizedPressureFunctionPtr)(const double Jn,
						const double Jr,
						const double *const *Fn,
						const double *const *Fr,
						const HOMMAT *material,
						double **linPres);

  typedef void (*d2UdJ2FunctionPtr)(const double Jn,
				    const double Jr,
				    const HOMMAT *material,
				    double *d2UdJ2);

  typedef void (*totalUpFunctPtr)(const double Tn,
				  const double Tr,
				  const HOMMAT *material,
				  double *Up);

  typedef void (*totalUppFunctPtr)(const double Tn,
				   const double Tr,
				   const HOMMAT *material,
				   double *Upp);

  /*=========================================================*/
  /* Retrieve functions                                      */
  /*=========================================================*/
  deviatoricStressFunctionPtr getDeviatoricStressFunction(const int flag,
							  const HOMMAT *material);

  materialStiffnessFunctionPtr getMaterialStiffnessFunction(const int flag,
							    const HOMMAT *material);

  volumetricPressureFunctionPtr getVolumetricPressureFunction(const int flag,
							      const HOMMAT *material);

  linearizedPressureFunctionPtr getLinearizedPressureFunction(const int flag,
							      const HOMMAT *material);

  d2UdJ2FunctionPtr getd2UdJ2Function(const int flag,
				      const HOMMAT *material);

  totalUpFunctPtr getTotalUpFunction(const int flag,
				     const HOMMAT *material);

  totalUppFunctPtr getTotalUppFunction(const int flag,
				       const HOMMAT *material);


  /*=========================================================*/
  /* Stress functions                                        */
  /*=========================================================*/

  /*** Second Piola-Kirchhoff Stress Functions ***/
  void deviatoricStress_MooneyRivlin(const double * const *C,
				     const HOMMAT *material,
				     double **S);

  void deviatoricStress_Linear(const double * const *C,
			       const HOMMAT *material,
			       double **S);

  /*=========================================================*/
  /* Tangent functions                                       */
  /*=========================================================*/

  /*** Lagrangian Stiffness Tensor Functions ***/
  void materialStiffness_MooneyRivlin(const double * const *C,
				      const HOMMAT *material,
				      const double * const *Fn,
				      const double * const *Fr,
				      double ****L);

  void materialStiffness_Linear(const double * const *C,
				const HOMMAT *material,
				const double * const *Fn,
				const double * const *Fr,
				double ****L);

  /*=========================================================*/
  /* Pressure functions                                      */
  /*=========================================================*/

  /*** Volumetric Pressure Function ***/
  void volumetricPressure_Common(const double Jn,
				 const double Jr,
				 const HOMMAT *material,
				 double *pres);

  void totalUp_Common(const double Tn,
		      const double Tr,
		      const HOMMAT *material,
		      double *Up);

  void volumetricPressure_DollSchweizerhof_7(const double Jn,
					     const double Jr,
					     const HOMMAT *material,
					     double *pres);

  void totalUp_DS7(const double Tn,
		   const double Tr,
		   const HOMMAT *material,
		   double *Up);

  void volumetricPressure_DollSchweizerhof_8(const double Jn,
					     const double Jr,
					     const HOMMAT *material,
					     double *pres);

  void totalUp_DS8(const double Tn,
		   const double Tr,
		   const HOMMAT *material,
		   double *Up);


  /*=========================================================*/
  /* Linearized Pressure functions                           */
  /*=========================================================*/

  /*** Linearized Volumetric Pressure Functions ***/
  void linearizedPressure_Common(const double Jn,
				 const double Jr,
				 const double *const *Fn,
				 const double *const *Fr,
				 const HOMMAT *material,
				 double **linPres);

  void linearizedPressure_DollSchweizerhof_7(const double Jn,
					     const double Jr,
					     const double *const *Fn,
					     const double *const *Fr,
					     const HOMMAT *material,
					     double **linPres);

  void linearizedPressure_DollSchweizerhof_8(const double Jn,
					     const double Jr,
					     const double *const *Fn,
					     const double *const *Fr,
					     const HOMMAT *material,
					     double **linPres);

  /*=========================================================*/
  /* Linearized Pressure functions (no deformation gradient) */
  /*=========================================================*/

  void d2UdJ2_Common(const double Jn,
		     const double Jr,
		     const HOMMAT *material,
		     double *d2UdJ2);

  void totalUpp_Common(const double Tn,
		       const double Tr,
		       const HOMMAT *material,
		       double *Upp);

  void d2UdJ2_DollSchweizerhof_7(const double Jn,
				 const double Jr,
				 const HOMMAT *material,
				 double *d2UdJ2);

  void totalUpp_DS7(const double Tn,
		    const double Tr,
		    const HOMMAT *material,
		    double *Upp);

  void d2UdJ2_DollSchweizerhof_8(const double Jn,
				 const double Jr,
				 const HOMMAT *material,
				 double *d2UdJ2);

  void totalUpp_DS8(const double Tn,
		    const double Tr,
		    const HOMMAT *material,
		    double *Upp);
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */
#endif /* #ifndef POTENTIAL_H */
