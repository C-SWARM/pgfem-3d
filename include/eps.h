/* HEADER */
#ifndef EPS_H
#define EPS_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

#ifndef ELEMENT_H
#include "element.h"
#endif

#ifndef VOLUMETRIC_DAMAGE_H
#include "volumetric_damage.h"
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
  typedef struct EPS EPS;

  EPS* build_eps_el (const long ne);

  EPS* build_eps_il (const long ne,
		     const ELEMENT *elem,
		     const int analysis);

  /*** MUST be called before destroy_elem */
  void destroy_eps_il(EPS* eps,
		      const ELEMENT *elem,
		      const long ne,
		      const int analysis);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
