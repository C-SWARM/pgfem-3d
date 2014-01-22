/* HEADER */

/** This file declares the element tangent and residual martices and
    vectors for updated lagrangian formulation of hyperelasticity. */

#ifndef TWO_FIELD_ELEMENT_H
#define TWO_FIELD_ELEMENT_H

#ifdef __cplusplus
extern "C" {
#endif

  /*======================================*/
  /*            TANGENTS                  */
  /*======================================*/

  /** Tangent stiffness matrix for displacement dofs.  TYPE allows for
      bubble enriched displacement elements such as the MINI element.
      TYPE = 0 => Kuu
      TYPE = 1 => Kbb
      TYPE = 2 => Kub
  */
  void UL_Kuu_at_ip(double *Kuu,
		    const int nne,
		    const int nne_t,
		    const double *ST,
		    const double *Fn,
		    const double *Fr,
		    const double *Fr_I,
		    const double Jn,
		    const double Jr,
		    const double pres,
		    const double *S,
		    const double *dSdF,
		    const double jj,
		    const double wt,
		    const int TYPE);

  /** Tangent stiffness matrix for displacement-pressure terms.  In
      general, Kup != (Kpu)'.  TYPE allows for bubble enriched
      displacement elements such as the MINI element.
      TYPE = 0 => Kup
      TYPE = 1 => Kbp
   */
  void UL_Kup_at_ip(double *Kup,
		    const int nne,
		    const int nne_t,
		    const double *Na,
		    const double *ST,
		    const double *Fr_I,
		    const double Jr, 
		    const double jj,
		    const double wt,
		    const int TYPE);

  /** Tangent stiffness matrix for pressure-displacement terms.  In
      general, Kpu != (Kup)'.  TYPE allows for bubble enriched
      displacement elements such as the MINI element.
      TYPE = 0 => Kpu
      TYPE = 1 => Kpb
  */
  void UL_Kpu_at_ip(double *Kpu,
		    const int nne,
		    const int nne_t,
		    const double *Na,
		    const double *ST,
		    const double *Fr_I,
		    const double Jn,
		    const double Jr,
		    const double Upp, 
		    const double jj,
		    const double wt,
		    const int TYPE);

  /** compute Kpu with damage */
  void damage_UL_Kpu_at_ip(double *Kpu,
			   const int nne,
			   const int nne_t,
			   const double *Na,
			   const double *ST,
			   const double *Fr_I,
			   const double Jn,
			   const double Jr,
			   const double kappa,
			   const double Upp,
			   const double UP, /* n+1 */
			   const double P,
			   const double *SS,
			   const double omega, 
			   const double jj,
			   const double wt,
			   const int TYPE);

  /** Tangent stiffness matrix for the pressure terms. */
  void UL_Kpp_at_ip(double *Kpp,
		    const int nne,
		    const int nne_t,
		    const double *Na,
		    const double kappa,
		    const double jj,
		    const double wt);

  void damage_UL_Kpp_at_ip(double *Kpp,
			   const int nne,
			   const int nne_t,
			   const double *Na,
			   const double kappa,
			   const double omega,
			   const double HH,
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
  void UL_Ru_at_ip(double *Ru,
		   const int nne,
		   const int nne_t,
		   const double *ST,
		   const double *Fn,
		   const double *Fr,
		   const double *Fr_I,
		   const double Jn,
		   const double Jr,
		   const double *S,
		   const double pres,
		   const double jj,
		   const double wt,
		   const int TYPE);

  /** Element residual for the pressure DOFs. */
  void UL_Rp_at_ip(double *Rp,
		   const int nne,
		   const double *Na,
		   const double kappa,
		   const double Up,
		   const double pres,
		   const double jj,
		   const double wt);
#ifdef __cplusplus
}
#endif

#endif
