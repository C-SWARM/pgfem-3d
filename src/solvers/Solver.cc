// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/// This file provides a home for the Solver abstract class' vtable, as well as
/// the constructor dispatch functionality.
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "pgfem3d/Solver.hpp"
#include "PGFEM_io.h"

pgfem3d::Solver::~Solver() {
}

int solve_system_check_error(FILE *out, const SOLVER_INFO info)
{
  switch(info.err){
   case SOLVE_SUCCESS:
    return 0;
   case SOLVE_ERROR:
    PGFEM_fprintf(out,"solve system returned error code SOLVE_ERROR.\n");
    return 1;
   case UNREC_SOLVER:
    PGFEM_fprintf(out,"solve system returned error code UNREC_SOLVER.\n");
    return 1;
   case UNREC_SOLVER_TYPE:
    PGFEM_fprintf(out,"solve system returned error code UNREC_SOLVER_TYPE.\n");
    return 1;
   case UNREC_PRECOND:
    PGFEM_fprintf(out,"solve system returned error code UNREC_PRECOND.\n");
    return 1;
   case BAD_PRECOND:
    PGFEM_fprintf(out,"solve system returned error code BAD_PRECOND.\n");
    return 1;
   default:
    PGFEM_fprintf(out,"solve system returned unrecognized error code (%d).\n",
                  info.err);
    return 1;
  }
}

