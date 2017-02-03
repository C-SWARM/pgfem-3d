/* HEADER */
#pragma once
#ifndef VOLUMETRIC_DAMAGE_H
#define VOLUMETRIC_DAMAGE_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Define volumetric damage variables */
  typedef struct DAMAGE_PARAMS{
    /* Weibull distribution parameters */
    double Yin, p1, p2;
    /* Viscous regularization parameter */
    double mu;
  } damage_params;

  typedef double (*vd_fun)(const double Ybar, const damage_params *params);

  typedef struct DAMAGE{
    /* *** damage parameters*** */
    damage_params params;

    /* *** damage variables *** */
    double wn,w; /* isotropic damage variable */
    double Hn,H; /* isotropic damage tangent */
    double Hpn,Hp; /* isotropic damage tangent rate (H') */
    double Xn,X; /* softening parameter */
    double dmu; /* dt*params.mu so I dont have to pass around dt */
    int damaged_n,damaged;  /* damage propagation flag */
    int broken;
    /* *** functions *** */
    int eq_flag;
    vd_fun function;
    vd_fun evolution;
    vd_fun evolution_rate;
  } damage;

  /** initialize the damage structure */
  void init_damage(damage *dam,
		   vd_fun function,
		   vd_fun evolution,
		   vd_fun evolution_rate,
		   const double *params,
		   const int len);

  /** initialize the damage structure using flags to set the function
      and evolution equations. */
  void init_damagef(damage *dam,
		    const int eq_flag,
		    const double *params,
		    const int len);

  /** Re-assign the damage functions according to eq_flag */
  void reset_damage_functions(damage *dam,
			      const int eq_flag);

  /** Copy the data from src to dest. The functions are reset
      according to src->eq_flag */
  void copy_damage(damage *dest,
		   const damage *src);

  /** Get the size of a damage object */
  size_t sizeof_damage(const damage *dam);

  /** Pack a damage object into a buffer */
  void pack_damage(const damage *src,
		   char *buffer,
		   size_t *pos);

  /** Unpack a damage object from a buffer */
  void unpack_damage(damage *dest,
		     const char *buffer,
		     size_t *pos);

  /** Reset the damage variables to n (e.g. after restart). */
  void reset_damage(damage *dam);

  /** Update the damage variables to n+1 (e.g. after completed step). */
  void update_damage(damage *dam);
 
  /** Damage integration algorithm */
  double damage_int_alg(damage *dam,
			const double Ybar,
			const double dt);

  /* *** Damage Equations *** */
  double weibull_function(const double Ybar,
			  const damage_params *params);

  /* *** Damage Evolution Equations *** */
  double weibull_evolution(const double Ybar,
			   const damage_params *params);

  /* *** Damage Evolution Rate Equations *** */
  double weibull_evolution_rate(const double Ybar,
				const damage_params *params);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef VOLUMETRIC_DAMAGE_H */
