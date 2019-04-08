/* HEADER */
/**
 * Declare functions for reading/manipulating the solver
 * file. Introduce a structure for storing/updating the data. All
 * functions return non-zero on (detected) error.
 *
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame <mmosby1@nd.edu>
 */

#include "solver_file.h"
#include "allocation.h"
#include "PGFEM_io.h"
#include "utils.h"
#include <string.h>

void solver_file_init_values(SOLVER_FILE *sf)
{
  sf->file = NULL;
  sf->nonlin_tol = 0;
  sf->max_nonlin_iter = 0;
  sf->n_pressure_nodes = 0;
  sf->nonlin_method = UNDEFINED_METHOD;
  sf->n_nonlin_method_opts = 0;
  sf->nonlin_method_opts = NULL;
  sf->n_step = 0;
  sf->times = NULL;
  sf->print_steps = NULL;
  sf->load_steps = NULL;
}

int solver_file_read_load(SOLVER_FILE *sf,
                          const size_t step,
                          const size_t len_load,
                          double *incr_load)
{
  int err = 0;
  if(step == 0 && sf->load_steps[step]){
    err++;
  }
  if(sf->load_steps[step]){
    for(size_t i = 0; i < len_load; i++){
      CHECK_SCANF(sf->file,"%lf",incr_load + i);
    }
  } else {
    memset(incr_load,0,len_load*sizeof(*incr_load));
  }
  return err;
}

int solver_file_scan_to_step(SOLVER_FILE *sf,
                             const size_t step,
                             const size_t len_load,
                             double *__restrict incr_load)
{
  int err = 0;

  double *__restrict tmp = PGFEM_malloc<double>(len_load);
  /* We start from i = 1 since the initial increment comes from a
     different file. */
  for(size_t i = 1; i <= step; i++){
    solver_file_read_load(sf,i,len_load,tmp);
    for(size_t j = 0; j < len_load; j++){
      incr_load[j] += tmp[j];
    }
  }
  free(tmp);
  return err;
}
