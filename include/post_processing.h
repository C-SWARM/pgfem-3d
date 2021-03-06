#ifndef PGFEM3D_POST_PROCESSING_H
#define PGFEM3D_POST_PROCESSING_H

#include "pgfem3d/Communication.hpp"
#include "PGFem3D_options.h"
#include "eps.h"
#include "femlib.h"
#include "new_potentials.h"
#include "PGFem3D_data_structure.h"

int read_from_VTK(const PGFem3D_opt *opts, int myrank, int step, double *u);

void post_processing_compute_stress_disp_ip(FEMLIB *fe, int e, double *S,
                                            HOMMAT *hommat, Element *elem,
                                            double *F, double Pn);

void post_processing_compute_stress(double *GS, Element *elem, HOMMAT *hommat,
                                    long ne, int npres, Node *node, EPS *eps,
                                    double* r, double *Vnp1, int ndofn,
				    const pgfem3d::CommunicationStructure *com,
                                    const PGFem3D_opt *opts);

void post_processing_deformation_gradient(double *GF, Element *elem,
                                          HOMMAT *hommat, long ne, int npres,
                                          Node *node, EPS *eps, double* r,
                                          int ndofn, const pgfem3d::CommunicationStructure *com,
                                          const PGFem3D_opt *opts);

void post_processing_deformation_gradient_elastic_part(double *GF,
                                                       Element *elem,
                                                       HOMMAT *hommat,
                                                       long ne, int npres,
                                                       Node *node, EPS *eps,
                                                       double* r, int ndofn,
                                                       const pgfem3d::CommunicationStructure *com,
                                                       const PGFem3D_opt *opts);

void post_processing_plastic_hardness(double *G_gn, Element *elem,
                                      HOMMAT *hommat, long ne, int npres,
                                      Node *node, EPS *eps, double* r,
                                      int ndofn, const pgfem3d::CommunicationStructure *com,
                                      const PGFem3D_opt *opts);

void post_processing_potential_energy(double *GE, Element *elem, HOMMAT *hommat,
                                      long ne, int npres, Node *node, EPS *eps,
                                      double* r, int ndofn,
				      const pgfem3d::CommunicationStructure *com,
                                      const PGFem3D_opt *opts);

void post_processing_deformed_volume(double *GV, Element *elem, long ne,
                                     Node *node, EPS *eps, double* r, int ndofn,
                                     const pgfem3d::CommunicationStructure *com,
                                     const PGFem3D_opt *opts);

/// compute and print maximum element pressure
void post_processing_max_pressure(const Grid &grid,
                                  const MaterialProperty &mat,
                                  const pgfem3d::CommunicationStructure *com,
                                  const FieldVariables &fv);

/// compute surface area of face on which 
/// Neumann boundary (NB) conditions are applied
void post_processing_compute_NBE_area(Grid &grid,
                                      const pgfem3d::CommunicationStructure *com);

#endif // #define PGFEM3D_POST_PROCESSING_H
