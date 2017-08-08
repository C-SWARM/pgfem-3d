#ifndef STRESS_STRAIN_H
#define STRESS_STRAIN_H

  /** Function returns Green Langange strain tensor E[3][3].
      Fn -> Fn-BAR || LOOK OUT */
  void get_GL_strain (double **Fn,
              double **Fr,
              double Jr,
              double Tr,
              double **E);

  /** Function returns Second Piola Kirchhoff stress tensor S[3][3]. */
  void get_SPK_stress (double L[3][3][3][3],
               double **E,
               double **S);

#endif /* #ifndef STRESS_STRAIN_H */
