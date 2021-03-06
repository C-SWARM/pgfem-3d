/* HEADER */
/**
 * Declare functions for reading/manipulating the solver
 * file. Introduce a structure for storing/updating the data. All
 * functions return non-zero on (detected) error.
 *
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame <mmosby1@nd.edu>
 */
#pragma once
#ifndef SOLVER_FILE_H
#define SOLVER_FILE_H

#include <stdio.h>
#include <stdlib.h>

struct SOLVER_FILE {
  FILE *file; /**< file to read data from */
  double nonlin_tol;
  size_t max_nonlin_iter;
  size_t n_pressure_nodes;
  size_t nonlin_method;
  size_t n_nonlin_method_opts;
  double *nonlin_method_opts;
  size_t n_step;
  double *times; /**< length of times, # computed steps + 1 */
  long   *print_steps; /**< T/F print on step */
  long   *load_steps; /**< T/F increment load on step */
};

enum {UNDEFINED_METHOD = -1, /**< gives overflow on assignment to nonlin_method */
      NEWTON_METHOD = 1,
      ARC_LENGTH_METHOD = 2,
      AUX_ARC_LENGTH_METHOD = 3,
      MULTIPHYSICS_NEWTON_METHOD = 4,
      TAYLOR_MODEL = 5};

void solver_file_init_values(SOLVER_FILE *sf);

/**
 * Read a load increment from the solver file.
 */
int solver_file_read_load(SOLVER_FILE *sf,
              const size_t step,
              const size_t len_load,
              double *incr_load);

/**
 * Scan the load section of the file to the specified step. The
 * specified loads are accumulated in incr_load. Note: incr_load is
 * *not* reset upon entry
 */
int solver_file_scan_to_step(SOLVER_FILE *sf,
                 const size_t step,
                 const size_t len_load,
                 double *incr_load);

#endif /* #ifndef  */
