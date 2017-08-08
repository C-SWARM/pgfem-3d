/* HEADER */
#ifndef TRANSFORM_COORDINATES_H
#define TRANSFORM_COORDINATES_H

  /** transform the coordinates (x,y,z) using a transformation matrix
      (or its transpose) built from the orthonormal basis (e1,e2,n) */
  void transform_coordinates(const int n_pts,
                 const double *x,
                 const double *y,
                 const double *z,
                 const double *e1,
                 const double *e2,
                 const double *n,
                 const int TRANS, /* boolean to use transpose */
                 double *tx,
                 double *ty,
                 double *tz);

  /** construct the transformation matrix from the orthonormal basis
      (e1,e2,n) as T = { e1 , e2 , n }' */
  void transformation_matrix(const double *e1,
                 const double *e2,
                 const double *n,
                 double *tm);

#endif /* #ifndef TRANSFORM_COORDINATES_H */
