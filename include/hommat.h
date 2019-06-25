/* HEADER */
#pragma once
#ifndef HOMMAT_H
#define HOMMAT_H

class Model_parameters;
/* #ifndef TYPE_Model_parameters */
/* #define TYPE_Model_parameters */
/* typedef class Model_parameters Model_parameters; */
/* #endif */

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

  Model_parameters *param; // model parameters constitutive model dependent
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
