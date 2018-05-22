/* HEADER */
/**
 * Routines for checkpointing/restart of FE2 simulations. Under the
 * FE2 framework, the original input files are used to perform the
 * necessary allocations and the restart files are read only to
 * populate the current state variables (including solution).
 *
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame <mmosby1@nd.edu>
 */
#pragma once
#ifndef PGF_FE2_RESTART_H
#define PGF_FE2_RESTART_H

#include "pgfem3d/MultiscaleCommon.hpp"
#include <stdlib.h> /* for size_t */

/**
 * Print a restart file for the macroscale domain.
 *
 * Note, calls update_M*CROSCALE_SOLUTION as this is not typically
 * necessary under the current implementation of the FE2 solver at the
 * macroscale. \return non-zero on error.
 */
int pgf_FE2_restart_print_macro(pgfem3d::Macroscale *macro);

/**
 * Print a restart file for the microscale domain.
 *
 * \return non-zero on error.
 */
int pgf_FE2_restart_print_micro(const pgfem3d::Microscale *micro,
				const size_t cell_id);

/**
 * Read a restart file for the macroscale domain and reset the
 * solution/state.
 *
 * \return non-zero on error.
 */
int pgf_FE2_restart_read_macro(pgfem3d::Macroscale *macro,
			       const size_t step,
			       const int mp_id);

/**
 * Read a restart file for a microscale domain and reset the
 * solution/state.
 *
 * \return non-zero on error.
 */
int pgf_FE2_restart_read_micro(pgfem3d::Microscale *micro,
			       const size_t step,
			       const size_t cell_id);

#endif /* #ifndef  */
