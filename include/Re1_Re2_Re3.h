#ifndef RE1_RE2_RE3_H
#define RE1_RE2_RE3_H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

  /** */
  void Re1_Re2_Re3 (long nne,
		    long ndn,
		    long npres,
		    double ai,
		    double aj,
		    double ak,
		    double J,
		    double *Psi,
		    double Jr,
		    double Jn,
		    double Tr,
		    double Tn,
		    double **Fr_I,
		    double ****ST,
		    double ****FF,
		    double **S,
		    double **f,
		    double pp,
		    double **Re1,
		    double **re1,
		    double *Re2,
		    double *re2,
		    double *Re3,
double *re3);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef RE1_RE2_RE3_H */
