#ifndef _H_POST_PROCESSING_H_
#define _H_POST_PROCESSING_H_

#include "new_potentials.h"
#include "femlib.h"
#include "eps.h"
#include "PGFem3D_options.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

int read_from_VTK(const PGFem3D_opt *opts, int myrank, int step, double *u);

void post_processing_compute_stress_disp_ip(FEMLIB *fe, int e, Matrix(double) S, HOMMAT *hommat, ELEMENT *elem, 
                          Matrix(double) F, double Pn);
void post_processing_compute_stress(double *GS, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts);
void post_processing_deformation_gradient(double *GF, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts);   
void post_processing_deformation_gradient_elastic_part(double *GF, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts);                                     
void post_processing_plastic_hardness(double *G_gn, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts);
void post_processing_potential_energy(double *GE, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts);
void post_processing_deformed_volume(double *GV, ELEMENT *elem, long ne, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, const PGFem3D_opt *opts);                    

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
