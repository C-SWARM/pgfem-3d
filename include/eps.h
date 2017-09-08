/* HEADER */
#pragma once
#ifndef EPS_H
#define EPS_H

#include "data_structure.h"
#include "element.h"
#include "volumetric_damage.h"
#include <stdio.h>
#include "state_variables.h"

#ifndef TYPE_CONSTITUTIVE_MODEL
#define TYPE_CONSTITUTIVE_MODEL
typedef struct Constitutive_model Constitutive_model;
#endif

  /** Structure of strains EPS */
  typedef struct { /* Inelastic strain in all integration points */
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
  }IL0_eps;

  /** Structure of strains EPS */
  typedef struct { /* Inelastic strain in all integration points */
    double *o,*f,*m,*d,*i;
  }IL1_eps;

  /** Structure of strains EPS */
  typedef struct { /* Stabilized */
    double *Fpp;
    damage dam; /* damage at extra int points */

    /* volumetric potential for damage */
    double Un_1, Un, Jn_1;
  }IL2_eps;

  /** ???.
      o - Overall strain
      f - Strain in the fibre
      m - Strain in the matrix
  */
  struct EPS{
    /** Elastic strain */
    struct {
      double *o,*f,*m,*i,
    *d,eq,eq_m,eq_i;
    }el;

    /** Plastic strain */
    struct {
      double *o,eq[2];
    }pl;

    IL0_eps *il;
    IL1_eps *d_il;

    /** Stabilized */
    IL2_eps *st;

    /** Volumetric Damage */
    damage *dam;

    /** Generalized constitutive modeling interface */
    Constitutive_model *model;

    /** Crystal plasticity */
    double *T,*d_T,GD;

    /** HOMOGENIZATION */
    double **F,
      **Fn,
      **P,
      **S,
      **Fe,
      **Fp,
      **FB,
      load1,load,
      **Dp,eff;
    long type;
  };

  EPS* build_eps_il (const long ne,
             const ELEMENT *elem,
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
             const int ne,
             const ELEMENT *elem,
             const int analysis);

  /**
   * return size of EPS for all elements in bytes
   */
  size_t sizeof_eps_list(const EPS *src,
             const int ne,
             const ELEMENT *elem,
             const int analysis);

  /**
   * Pack EPS for all elements into a buffer.
   */
  void pack_eps_list(const EPS *src,
             const int ne,
             const ELEMENT *elem,
             const int analysis,
             char *buffer,
             size_t *pos);

  /**
   * Unpack EPS for all elements from a buffer.
   */
  void unpack_eps_list(EPS *dest,
               const int ne,
               const ELEMENT *elem,
               const int analysis,
               const char *buffer,
               size_t *pos);

  /*** MUST be called before destroy_elem */
  void destroy_eps_il(EPS* eps,
              const ELEMENT *elem,
              const long ne,
              const int analysis);

#endif /* #ifndef  */
