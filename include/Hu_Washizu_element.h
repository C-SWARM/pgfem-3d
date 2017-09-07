/* HEADER */

/**
 * @file This file declares the element tangent and residual martices
 * and vectors for updated lagrangian formulation of three-field
 * Hu-Washizu hyperelasticity.
 *
 * AUTHORS:
 * Matthew Mosby
*/

#pragma once
#ifndef HU_WASHIZU_ELEMENT_H
#define HU_WASHIZU_ELEMENT_H

  /*======================================*/
  /*            TANGENTS                  */
  /*======================================*/

  /** Tangent stiffness matrix for displacement dofs.  TYPE allows for
      bubble enriched displacement elements such as the MINI element.
      TYPE = 0 => Kuu
      TYPE = 1 => Kbb
      TYPE = 2 => Kub
  */
  void HW_Kuu_at_ip(double *Kuu,
            const int nne,
            const int nne_t,
            const double *ST,
            const double *Fn,
            const double *Fr,
            const double *Fr_I,
            const double Jn,
            const double Tn,
            const double Jr,
            const double pres,
            const double *S,
            const double *L,
            const double jj,
            const double wt,
            const int TYPE);

  /** Tangent stiffness matrix for displacement-pressure terms.TYPE
      allows for bubble enriched displacement elements. Kup = (Kpu)'.
      TYPE = 0 => Kup
      TYPE = 1 => Kbp
   */
  void HW_Kup_at_ip(double *Kup,
            const int nne,
            const int nne_t,
            const int nPres,
            const double *Np,
            const double *ST,
            const double *Fr_I,
            const double Jn,
            const double Tn,
            const double Jr,
            const double jj,
            const double wt,
            const int TYPE);

  /** Tangent stiffness for volume DOFs. */
  void HW_Ktt_at_ip(double *Ktt,
            const int nVol,
            const double *Nt,
            const double Tn,
            const double Jn,
            const double kappa,
            const double Upp,
            const double jj,
            const double wt);

  /** Tangent stiffness for pressure-volume DOFs. Kpt = (Ktp)' */
  void HW_Kpt_at_ip(double *Kpt,
            const int nPres,
            const double *Np,
            const int nVol,
            const double *Nt,
            const double Tn,
            const double Jn,
            const double jj,
            const double wt);

  /*======================================*/
  /*            RESIDUALS                 */
  /*======================================*/

  /** Element residual for the displacement DOFs.  TYPE allows for
      bubble enriched displacement elements such as the MINI element.
      TYPE = 0 => Ru
      TYPE = 1 => Rb
   */
  void HW_Ru_at_ip(double *Ru,
           const int nne,
           const int nne_t,
           const double *ST,
           const double *Fn,
           const double *Fr,
           const double *Fr_I,
           const double Jn,
           const double Tn,
           const double Jr,
           const double Tr,
           const double *S,
           const double pres,
           const double jj,
           const double wt,
           const int TYPE);

  void HW_Ru1_at_ip(double *Ru,
            const int nne,
            const int nne_t,
            const double *ST,
            const double *Fn,
            const double *Fr,
            const double *Fr_I,
            const double Jn,
            const double Tn,
            const double Jr,
            const double Tr,
            const double *S,
            const double pres,
            const double jj,
            const double wt,
            const int TYPE);

  void HW_Ru2_at_ip(double *Ru,
            const int nne,
            const int nne_t,
            const double *ST,
            const double *Fn,
            const double *Fr,
            const double *Fr_I,
            const double Jn,
            const double Tn,
            const double Jr,
            const double *S,
            const double pres,
            const double jj,
            const double wt,
            const int TYPE);

  /** Element residual for the pressure DOFs. */
  void HW_Rp_at_ip(double *Rp,
           const int nPres,
           const double *Np,
           const double Jn,
           const double Jr,
           const double Tn,
           const double Tr,
           const double jj,
           const double wt);

  /** Element residual for the Volume DOFs */
  void HW_Rt_at_ip(double *Rt,
           const int nVol,
           const double *Nt,
           const double Tn,
           const double Jn,
           const double pres,
           const double kappa,
           const double Up,
           const double jj,
           const double wt);

  void debug_NH_Ru_at_ip(double *Ru,
             const int nne,
             const int nne_t,
             const double *ST,
             const double *Fn,
             const double *Fr,
             const double *Fr_I,
             const double Jn,
             const double Tn,
             const double Jr,
             const double Tr,
             const double *S,
             const double pres,
             const double G,
             const double jj,
             const double wt,
             const int TYPE);

#endif
