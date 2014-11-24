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

#include <string.h>
#include "PGFEM_io.h"

static void solver_file_init_values(SOLVER_FILE *sf)
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

int solver_file_open(const char *filename,
		     SOLVER_FILE **sf)
{
  int err = 0;
  *sf = malloc(sizeof(**sf));
  solver_file_init_values(*sf);
  (*sf)->file = PGFEM_fopen(filename,"r");
  return err;
}

int solver_file_close(SOLVER_FILE *sf)
{
  int err = 0;
  PGFEM_fclose(sf->file);
  free(sf->nonlin_method_opts);
  free(sf->times);
  free(sf->print_steps);
  free(sf->load_steps);
  free(sf);
  return err;
}

int solver_file_read_header(SOLVER_FILE *sf)
{
  int err = 0;

  /* read first line */
  fscanf(sf->file,"%lf %ld %ld %ld",
	 &(sf->nonlin_tol),
	 &(sf->max_nonlin_iter),
	 &(sf->n_pressure_nodes),
	 &(sf->nonlin_method));

  switch(sf->nonlin_method){
  case NEWTON_METHOD:
    sf->n_nonlin_method_opts = 0;
    sf->nonlin_method_opts = NULL;
    break;

  case ARC_LENGTH_METHOD:
  case AUX_ARC_LENGTH_METHOD:
    /* read additional options */
    sf->n_nonlin_method_opts = 2;
    sf->nonlin_method_opts = malloc((sf->n_nonlin_method_opts)
				    *sizeof(*(sf->nonlin_method_opts)));
    fscanf(sf->file,"%lf %lf",
	   sf->nonlin_method_opts,
	   sf->nonlin_method_opts + 1);
    break;

  default:
    PGFEM_printerr("Undefined nonlinear solution method!\n");
    sf->nonlin_method = UNDEFINED_METHOD;
    sf->n_nonlin_method_opts = 0;
    sf->nonlin_method_opts = NULL;
    err++;
    break;
  }

  /* read number of steps and allocate */
  fscanf(sf->file,"%ld",&(sf->n_step)); /* # computed steps */
  sf->print_steps = calloc(sf->n_step,sizeof(*(sf->print_steps)));
  sf->load_steps = calloc(sf->n_step,sizeof(*(sf->load_steps)));

  /* allocate/read times to read */
  sf->times = malloc((sf->n_step + 1) * sizeof(*(sf->times)));
  for(size_t i = 0, e = sf->n_step + 1; i < e; i++){
    fscanf(sf->file,"%lf",(sf->times) + i);
  }

  /* read print steps */
  size_t n_step = 0; /* NOT sf->n_step!! */
  fscanf(sf->file,"%ld",&n_step);
  for(size_t i = 0; i < n_step; i++){
    size_t idx = 0;
    fscanf(sf->file,"%ld",&idx);
    sf->print_steps[idx] = 1;
  }

  /* read load steps */
  fscanf(sf->file,"%ld",&n_step);
  for(size_t i = 0; i < n_step; i++){
    size_t idx = 0;
    fscanf(sf->file,"%ld",&idx);
    sf->load_steps[idx] = 1;
  }

  return err;
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
      fscanf(sf->file,"%lf",incr_load + i);
    }
  } else {
    memset(incr_load,0,len_load*sizeof(*incr_load));
  }
  return err;
}

int solver_file_scan_to_step(SOLVER_FILE *sf,
			     const size_t step,
			     const size_t len_load,
			     double *restrict accum_load,
			     double *restrict incr_load)
{
  int err = 0;

  /* start with empty increment */
  memset(incr_load,0,len_load*sizeof(*incr_load));

  /* note that we increment to i = step so that we get increment for
     current step. Also, we start from i = 1 since the initial
     increment comes from a different file. */
  for(size_t i = 1; i <= step; i++){
    for(size_t j = 0; j < len_load; j++){
      accum_load[j] += incr_load[j];
    }
    solver_file_read_load(sf,i,len_load,incr_load);
  }
  return err;
}
