#ifndef DEF_GRAD_H
#define DEF_GRAD_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  void def_grad_get (const long nne,
		     const long ndofn,
		     const double ****AA,
		     const double *r_e,
		     double **F);
  
  void def_grad_inv (const double *const *F,
		     double **F_I);
  
  double def_grad_det (const double *const *F);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef DEF_GRAD_H */
