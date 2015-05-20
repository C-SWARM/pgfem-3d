#ifndef _H_POST_PROCESSING_H_
#define _H_POST_PROCESSING_H_

#include "new_potentials.h"
#include "femlib.h"
#include "eps.h"

void post_processing_compute_stress_disp_ip(FEMLIB *fe, int e, Matrix(double) S, HOMMAT *hommat, ELEMENT *elem, 
                          Matrix(double) F, double Pn);
void post_processing_compute_stress(double *GS, ELEMENT *elem, HOMMAT *hommat, long ne, int npres, NODE *node, EPS *eps,
                    double* r, int ndofn, MPI_Comm mpi_comm, int analysis);
#endif