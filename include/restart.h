#ifndef __H_PGEM3D_RESTART_H__
#define __H_PGEM3D_RESTART_H__

#include "PGFem3D_options.h"

#ifndef NO_VTK_LIB
#include "PGFem3D_to_VTK.hpp"
#endif

int read_initial_from_VTK(const PGFem3D_opt *opts, int myrank, int *restart, double *u0, double *u1);

int read_restart(double *u0, double *u1, const PGFem3D_opt *opts, 
                 ELEMENT *elem, NODE *node, SIG * sig_e, EPS *eps, SUPP sup,
                 int myrank, int elemno, int nodeno, int nsd, int *stepno);
                 
int write_restart(double *u0, double *u1, const PGFem3D_opt *opts, 
                  ELEMENT *elem, NODE *node, SIG * sig_e, EPS *eps, SUPP sup,                  
                  int myrank, int elemno, int nodeno, int ndofn, int ndofd, int stepno);
#endif