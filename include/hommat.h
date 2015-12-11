/* HEADER */
#pragma once
#ifndef HOMMAT_H
#define HOMMAT_H

/** Structure of material matrices */
struct HOMMAT{
  /* Elastic stiffness matrix */
  double *M,*L;
  /* Mooney - Rivlin material parameters */
  double m10,m01,E,G,nu,density;

  /* extra params from file */
  double e1,e2,e3,e4;

  /* potential function flags */
  int devPotFlag;
  int volPotFlag;
  int mat_id; // material ID
};
#ifndef TYPE_HOMMAT
#define TYPE_HOMMAT
typedef struct HOMMAT HOMMAT;
#endif

HOMMAT* build_hommat(long i);

void destroy_hommat(HOMMAT* hm, long nm);

/**
 * \return bulk modulus (kappa) computed from linear elastic
 * properties.
 */
double hommat_get_kappa(const HOMMAT *mat);

#endif /* #ifndef  */
