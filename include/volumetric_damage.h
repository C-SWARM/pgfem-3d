/* HEADER */
#pragma once
#ifndef VOLUMETRIC_DAMAGE_H
#define VOLUMETRIC_DAMAGE_H

#include <stdlib.h>

  /** Define volumetric damage variables */
  struct DAMAGE_PARAMS{
    /* Weibull distribution parameters */
    double Yin = 0.0;
    double p1 = 0.0;
    double p2 = 0.0;
    /* Viscous regularization parameter */
    double mu = 0.0;
  };
  using DamageParams = DAMAGE_PARAMS;

  typedef double (*vd_fun)(const double Ybar, const DamageParams *params);

  struct DAMAGE{
    /* *** damage parameters*** */
    DamageParams params = {};

    /* *** damage variables *** */
    double wn = 0.0, w = 0.0; /* isotropic damage variable */
    double Hn = 0.0, H = 0.0; /* isotropic damage tangent */
    double Hpn = 0.0, Hp = 0.0; /* isotropic damage tangent rate (H') */
    double Xn = 0.0, X = 0.0; /* softening parameter */
    double dmu = 0.0; /* dt*params.mu so I dont have to pass around dt */
    int damaged_n,damaged = 0;  /* damage propagation flag */
    int broken = 0;
    /* *** functions *** */
    int eq_flag = 0;
    vd_fun function = nullptr;
    vd_fun evolution = nullptr;
    vd_fun evolution_rate = nullptr;
  };
  using Damage = DAMAGE;

  /** initialize the damage structure */
  void init_damage(Damage *dam,
           vd_fun function,
           vd_fun evolution,
           vd_fun evolution_rate,
           const double *params,
           const int len);

  /** initialize the damage structure using flags to set the function
      and evolution equations. */
  void init_damagef(Damage *dam,
            const int eq_flag,
            const double *params,
            const int len);

  /** Re-assign the damage functions according to eq_flag */
  void reset_damage_functions(Damage *dam,
                  const int eq_flag);

  /** Copy the data from src to dest. The functions are reset
      according to src->eq_flag */
  void copy_damage(Damage *dest,
           const Damage *src);

  /** Get the size of a damage object */
  size_t sizeof_damage(const Damage *dam);

  /** Pack a damage object into a buffer */
  void pack_damage(const Damage *src,
           char *buffer,
           size_t *pos);

  /** Unpack a damage object from a buffer */
  void unpack_damage(Damage *dest,
             const char *buffer,
             size_t *pos);

  /** Reset the damage variables to n (e.g. after restart). */
  void reset_damage(Damage *dam);

  /** Update the damage variables to n+1 (e.g. after completed step). */
  void update_damage(Damage *dam);

  /** Damage integration algorithm */
  double damage_int_alg(Damage *dam,
            const double Ybar,
            const double dt);

  /* *** Damage Equations *** */
  double weibull_function(const double Ybar,
              const DamageParams *params);

  /* *** Damage Evolution Equations *** */
  double weibull_evolution(const double Ybar,
               const DamageParams *params);

  /* *** Damage Evolution Rate Equations *** */
  double weibull_evolution_rate(const double Ybar,
                const DamageParams *params);

#endif /* #ifndef VOLUMETRIC_DAMAGE_H */
