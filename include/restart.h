#ifndef __H_PGEM3D_RESTART_H__
#define __H_PGEM3D_RESTART_H__

#include "PGFem3D_options.h"
#include "PGFem3D_data_structure.h"

#ifndef NO_VTK_LIB
#include "PGFem3D_to_VTK.hpp"
#endif

/// read restart files for multiphysics problem
///
/// \param[in] grid a mesh object
/// \param[in, out] fv array of field variable object
/// \param[in] load object for loading
/// \param[in] opts PGFem3D commend line options
/// \param[in] mp multiphysics object
/// \param[out] tnm1 times at t(n-1), t(n) 
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_restart(GRID *grid,
                  FIELD_VARIABLES *fv,
                  LOADING_STEPS *load,
                  const PGFem3D_opt *opts,
                  MULTIPHYSICS *mp,
                  double *tnm1,
                  int myrank);

/// write restart files for mechanical part
///
/// \param[in] grid a mesh object
/// \param[in, out] fv array of field variable object
/// \param[in] load object for loading
/// \param[in] time_steps object for time stepping
/// \param[in] opts PGFem3D commend line options
/// \param[in] mp multiphysics object
/// \param[in] myrank current process rank
/// \param[in] mp_id multiphysics id
/// \param[in] stepno current time step number
/// \param[in] rs_path directory path for restart files
/// \return non-zero on internal error                                   
int write_restart(GRID *grid,
                  FIELD_VARIABLES *fv,
                  LOADING_STEPS *load,
                  PGFem3D_TIME_STEPPING *time_steps,
                  const PGFem3D_opt *opts,
                  MULTIPHYSICS *mp,
                  int stepno,
                  int myrank);
#endif