/* HEADER */
#ifndef HOMMAT_H
#define HOMMAT_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** Structure of material matrices */
  struct HOMMAT{
    /* Elastic stiffness matrix */
    double *M,*L;
    /* Mooney - Rivlin material parameters */
    double m10,m01,E,G,nu;

    /* extra params from file */
    double e1,e2,e3,e4;

    /* potential function flags */
    int devPotFlag;
    int volPotFlag;
  };
  typedef struct HOMMAT HOMMAT;

  HOMMAT* build_hommat(long i);

  void destroy_hommat(HOMMAT* hm, long nm);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef  */
