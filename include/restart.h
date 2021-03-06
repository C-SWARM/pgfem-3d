#ifndef PGEM3D_RESTART_H
#define PGEM3D_RESTART_H

#include "PGFem3D_options.h"
#include "PGFem3D_data_structure.h"

#ifndef NO_VTK_LIB
#include "PGFem3D_to_VTK.hpp"
#endif

/// read restart files for multiphysics problem
///
/// \param[in] grid a mesh object
/// \param[in, out] fv array of field variable object
/// \param[in, out] time_steps object for time stepping
/// \param[in] load object for loading
/// \param[in] opts PGFem3D commend line options
/// \param[in] mp multiphysics object
/// \param[out] tnm1 times at t(n-1), t(n)
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_restart(Grid *grid,
                 FieldVariables *fv,
                 TimeStepping *time_steps,
                 LoadingSteps *load,
                 const PGFem3D_opt *opts,
                 const Multiphysics& mp,
                 double *tnm1,
                 int myrank);

/// write restart files for mechanical part
///
/// \param[in] grid   a mesh object
/// \param[in, out]   fv array of field variable object
/// \param[in] opts   PGFem3D commend line options
/// \param[in] mp     multiphysics object
/// \param[in] tns    times  at n for multiple physics
/// \param[in] tnm1   times at t(n-1), t(n), t(n+1)
/// \param[in] myrank current process rank
/// \param[in] stepno current time step number
/// \return non-zero on internal error
int write_restart(Grid *grid,
                  FieldVariables *fv,
                  LoadingSteps *load,
                  const PGFem3D_opt *opts,
                  const Multiphysics& mp,
                  const double *tns,
                  const double *tnm1,
                  int myrank,
                  int stepno);

int read_initial_from_VTK(PGFem3D_opt *opts, int myrank, int *restart, double *u0, double *u1);

int read_restart(double *u0, double *u1, const PGFem3D_opt *opts,
                 Element *elem, Node *node, SIG * sig_e, EPS *eps, SUPP sup,
                 int myrank, int elemno, int nodeno, int nsd, int *stepno, double *tnm1, double *NORM);

int write_restart(double *u0, double *u1, const PGFem3D_opt *opts,
                  Element *elem, Node *node, SIG * sig_e, EPS *eps, SUPP sup,
                  int myrank, int elemno, int nodeno, int ndofn, int ndofd, int stepno, double *times, double NORM);

#endif // #define PGFEM3D_RESTART_H
