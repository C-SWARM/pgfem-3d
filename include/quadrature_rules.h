/* HEADER */

/* This is a set of functions for retrieving quadrature rules of a
   given order for a given element type. Given the desired quadrature
   order, return the number of integration points and the
   isoparametric coordinates of the integration points. Note that ksi,
   eta, zet and wt are allocated internally and must be deallocated by
   the controlling function. */

#ifndef PGFEM_QUADRRATURE_RULES_H
#define PGFEM_QUADRRATURE_RULES_H

  int get_tria_quadrature_rule(const int int_order,
                   int *n_ip,
                   double **ksi,
                   double **eta,
                   double **wt);

  int get_quad_quadrature_rule(const int int_order,
                   int *n_ip,
                   double **ksi,
                   double **eta,
                   double **wt);

  int get_tet_quadrature_rule(const int int_order,
                  int *n_ip,
                  double **ksi,
                  double **eta,
                  double **zet,
                  double **wt);

  int get_hex_quadrature_rule(const int int_order,
                  int *n_ip,
                  double **ksi,
                  double **eta,
                  double **zet,
                  double **wt);

  int get_wedge_quadrature_rule(const int int_order,
                int *n_ip,
                double **ksi,
                double **eta,
                double **zet,
                double **wt);

  int get_pyram_quadrature_rule(const int int_order,
                int *n_ip,
                double **ksi,
                double **eta,
                double **zet,
                double **wt);

#endif /* #ifndef  */
