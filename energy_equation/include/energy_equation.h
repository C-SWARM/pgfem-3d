#pragma once
#ifndef H__PGFEM3D_ENERGY_EQUATION__H
#define H__PGFEM3D_ENERGY_EQUATION__H

#include "PGFem3D_data_structure.h"

int energy_equation_compute_residuals(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      const int mp_id);

int energy_equation_compute_stiffness(GRID *grid,
                                      MATERIAL_PROPERTY *mat,
                                      FIELD_VARIABLES *fv,
                                      SOLVER_OPTIONS *sol,
                                      COMMUNICATION_STRUCTURE *com,
                                      MPI_Comm mpi_comm,
                                      int myrank,
                                      const PGFem3D_opt *opts,
                                      const int mp_id);

#endif
