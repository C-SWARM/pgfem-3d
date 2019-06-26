/* HEADER */
#pragma once
#ifndef PGFEM3D_EPS_H
#define PGFEM3D_EPS_H

#include "data_structure.h"
#include "element.h"
#include "state_variables.h"
#include "volumetric_damage.h"
#include <stdio.h>

class Constitutive_model;

/** Structure of strains EPS */
struct IL0_eps { /* Inelastic strain in all integration points */
  double *o,*f,*m,*d,*i;

  /* volumetric potential for damage */
  double Un_1, Un, Jn_1;

  /* virgin potential */
  double Y;

  /** crystal plasticity */
  double *F,
    lam,
    *Fp,
    *UU,
    ****dUU_Fr,
    **dUU_Tr,
    ****dUU_Fr_n,
    **dUU_Tr_n,
    *Fe,
    *Fe1,
    eff,
    *GA,
    *GA1,
    GAMA;
  long *PLC_B;
  /* For plasticity F denotes Fe => F = FeFp */
};

/** Structure of strains EPS */
struct IL1_eps { /* Inelastic strain in all integration points */
  double *o,*f,*m,*d,*i;
};

/** Structure of strains EPS */
struct IL2_eps { /* Stabilized */
  double *Fpp = nullptr;
  Damage dam = {}; /* damage at extra int points */

  /* volumetric potential for damage */
  double Un_1 = 0.0, Un = 0.0, Jn_1 = 0.0;
};

/** ???.
    o - Overall strain
    f - Strain in the fibre
    m - Strain in the matrix
*/
struct EPS{
  /** Elastic strain */
  struct {
    double *o = nullptr,
      *f = nullptr,
      *m = nullptr,
      *i = nullptr,
      *d = nullptr,
      eq = 0.0, eq_m = 0.0, eq_i = 0.0;
  } el = {};

  /** Plastic strain */
  struct {
    double *o = nullptr, eq[2] = {0.0, 0.0};
  } pl = {};

  IL0_eps *il = nullptr;
  IL1_eps *d_il = nullptr;

  /** Stabilized */
  IL2_eps *st = nullptr;

  /** Volumetric Damage */
  Damage *dam = nullptr;

  /** Generalized constitutive modeling interface */
  Constitutive_model *model = nullptr;

  /** Crystal plasticity */
  double *T = nullptr, *d_T = nullptr, GD = 0.0;

  /** HOMOGENIZATION */
  double   **F = nullptr;
  double  **Fn = nullptr;
  double   **P = nullptr;
  double   **S = nullptr;
  double  **Fe = nullptr;
  double  **Fp = nullptr;
  double  **FB = nullptr;
  double load1 = 0.0;
  double  load = 0.0;
  double  **Dp = nullptr;
  double   eff = 0.0;
  long    type = 0;
};

EPS* build_eps_il (const long ne,
                   const Element *elem,
                   const int analysis,
                   State_variables **statv_list);

/**
 * Copy EPS for all elements. Additionally copies the rate of
 * plastic strain (which is stored only on the 0-th element). The
 * lists src and dest should be constructed by identical calls to
 * build_eps_il.
 */
void copy_eps_list(EPS *dest,
                   const EPS *src,
                   const long ne,
                   const Element *elem,
                   const int analysis);

/**
 * return size of EPS for all elements in bytes
 */
size_t sizeof_eps_list(const EPS *src,
                       const long ne,
                       const Element *elem,
                       const int analysis);

/**
 * Pack EPS for all elements into a buffer.
 */
void pack_eps_list(const EPS *src,
                   const long ne,
                   const Element *elem,
                   const int analysis,
                   char *buffer,
                   size_t *pos);

/**
 * Unpack EPS for all elements from a buffer.
 */
void unpack_eps_list(EPS *dest,
                     const long ne,
                     const Element *elem,
                     const int analysis,
                     const char *buffer,
                     size_t *pos);

/*** MUST be called before destroy_elem */
void destroy_eps_il(EPS* eps,
                    const Element *elem,
                    const long ne,
                    const int analysis);

#endif /* #define PGFEM3D_EPS_H  */
