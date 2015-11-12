/* HEADER */
#ifndef NEW_POTENTIALS_H
#define NEW_POTENTIALS_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#include "hommat.h"

  /* This file contains functions to return stress, pressure and their
     tangents.  The tangents are the derivatives with respect to their
     native variable, i.e. it is up to the user to take care of chain
     rule. e.g. dev_S = 2*d (dev_W)/dC, dev_L = 2*d (dev_S)/dC */

  /* Quick drivers */
  void new_pot_compute_Wdev(const double *Ce,
                            const HOMMAT *p_hmat,
                            double *Wdev);

  void new_pot_compute_U(const double J,
                         const HOMMAT *p_hmat,
                         double *U);

  void new_pot_compute_Sdev(const double *Ce,
                            const HOMMAT *p_hmat,
                            double *Sdev);

  void new_pot_compute_dudj(const double Je,
                            const HOMMAT *p_hmat,
                            double *dudj);

  void new_pot_compute_Ldev(const double *Ce,
                            const HOMMAT *p_hmat,
                            double *Ldev);

  void new_pot_compute_d2udj2(const double Je,
                              const HOMMAT *p_hmat,
                              double *d2udj2);

  /*===== Function types ====*/
  typedef void (*devPotentialFuncPtr)(double const *C,
				      HOMMAT const *mat,
				      double *W);

  typedef void (*devStressFuncPtr)(double const *C,
				   HOMMAT const *mat,
				   double *S);

  typedef void (*matStiffFuncPtr)(double const *C,
				  HOMMAT const *mat,
				  double *L);

  typedef void (*matSensFuncPtr)(double const *C,
				 HOMMAT const *mat,
				 double *K);

  typedef void (*UFuncPtr)(double const J,
			   HOMMAT const *mat,
			   double *U);

  typedef void (*dUdJFuncPtr)(double const J,
			      HOMMAT const *mat,
			      double *dUdJ);

  typedef void (*d2UdJ2FuncPtr)(double const J,
				HOMMAT const *mat,
				double *d2UdJ2);

  typedef void (*d3UdJ3FuncPtr)(double const J,
				HOMMAT const *mat,
				double *d3UdJ3);

  /*==== Get functions pointers ====*/
  devPotentialFuncPtr getDevPotentialFunc(const int flag,
					  HOMMAT const * mat);

  devStressFuncPtr getDevStressFunc(const int flag,
				    HOMMAT const *mat);

  matStiffFuncPtr getMatStiffFunc(const int flag,
				  HOMMAT const *mat);

  matSensFuncPtr getMatSensFunc(const int flag,
				HOMMAT const *mat);

  UFuncPtr getUFunc(const int flag,
		    HOMMAT const *mat);

  dUdJFuncPtr getDUdJFunc(const int flag,
			  HOMMAT const *mat);

  d2UdJ2FuncPtr getD2UdJ2Func(const int flag,
			      HOMMAT const *mat);

  d3UdJ3FuncPtr getD3UdJ3Func(const int flag,
			      HOMMAT const *mat);

  /*==== Deviatoric potential functions ====*/
  void devPotential_Mooney_Rivlin(double const *C,
				  HOMMAT const *mat,
				  double *W);

  void devPotential_Linear(double const *C,
			   HOMMAT const *mat,
			   double *W);

  /*==== Deviatoric stress functions ====*/
  void devStress_Mooney_Rivlin(double const *C,
			       HOMMAT const *mat,
			       double *S);

  void devStress_Linear(double const *C,
			HOMMAT const *mat,
			double *S);

  /*==== Material stiffness functions ====*/
  void matStiff_Mooney_Rivlin(double const *C,
			      HOMMAT const *mat,
			      double *L);

  void matStiff_Linear(double const *C,
		       HOMMAT const *mat,
		       double *L);

  /*==== Material sensitivity functions ====*/
  void matSens_Mooney_Rivlin(double const *C,
			     HOMMAT const *mat,
			     double *L);

  void matSens_Linear(double const *C,
		      HOMMAT const *mat,
		      double *L);

  /*==== Volumetric Potetntial Functions ====*/
  void U_Common(double const J,
		HOMMAT const *mat,
		double *U);

  void U_Doll_Schweizerhof_7(double const J,
			     HOMMAT const *mat,
			     double *U);

  void U_Doll_Schweizerhof_8(double const J,
			     HOMMAT const *mat,
			     double *U);

  /*==== dUdJ functions ====*/
  void dUdJ_Common(double const J,
		   HOMMAT const *mat,
		   double *dUdJ);

  void dUdJ_Doll_Schweizerhof_7(double const J,
				HOMMAT const *mat,
				double *dUdJ);

  void dUdJ_Doll_Schweizerhof_8(double const J,
				HOMMAT const *mat,
				double *dUdJ);

  /*==== d2UdJ2 functions ===*/
  void d2UdJ2_Common_new(double const J,
			 HOMMAT const *mat,
			 double *dUdJ);

  void d2UdJ2_Doll_Schweizerhof_7(double const J,
				  HOMMAT const *mat,
				  double *d2UdJ2);

  void d2UdJ2_Doll_Schweizerhof_8(double const J,
				  HOMMAT const *mat,
				  double *d2UdJ2);

  /*==== d3UdJ3 functions ===*/
  void d3UdJ3_Common_new(double const J,
			 HOMMAT const *mat,
			 double *d3UdJ3);

  void d3UdJ3_Doll_Schweizerhof_7(double const J,
				  HOMMAT const *mat,
				  double *d3UdJ3);

  void d3UdJ3_Doll_Schweizerhof_8(double const J,
				  HOMMAT const *mat,
				  double *d3UdJ3);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
