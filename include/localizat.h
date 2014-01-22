#ifndef LOCALIZAT_H
#define LOCALIZAT_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Localizes the matrices of one layer on one element to the
      stiffness matrix of subdomain.
      %%%%%%%%%%%%%%%% TESTED 19.05.2000 %%%%%%%%%%%%%%%%%
  */
  void localizat (double *k, double *lk, long *adr, long *cn, long ndofe);

  /** Localizes the matrices of one layer on one element to the
      stiffness matrix of subdomain.
      
      %%%%%%%%%%%%%%%% TESTED 19.05.2000 %%%%%%%%%%%%%%%%%
  */
  void localizat_nonsym (double *kL, double *kU, double *lk,
			 long *adr, long *cn, long ndofe);

  /** Sparse nonsymmetric column storage format.
       
            (i)
          xxxxxx
      (j) xxxxxx
          xxxxxx
  */
  void localizat_sparse (double *K, double *lk, long *Ai,
			 long *Ap, long *cn, long ndofe);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef LOCALIZAT_H */
