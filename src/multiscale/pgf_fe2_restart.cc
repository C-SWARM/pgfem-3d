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

#include "pgf_fe2_restart.h"
#include "gen_path.h"
#include "utils.h"
#include "PGFEM_io.h"

#include "enumerations.h"

#include "fd_increment.h"
#include "stabilized.h"
#include "MINI_element.h"
#include "MINI_3f_element.h"
#include "displacement_based_element.h"
#include "cohesive_element.h"

#include <stdio.h>
#include <assert.h>

using namespace pgfem3d;
using namespace multiscale::net;

/**
 * Generate the restart path string and ensure that the path exists,
 * creating directories as needed.
 *
 * User must ensure that restart_path_str is available for allocation
 * upon entry. User is responsible for releasing memory pointed to by
 * restart_path_str on exit. \return non-zero on error.
 */
static int generate_restart_path(const char *opath,
				 const size_t cell_id,
				 const size_t step,
				 char **restart_path_str)
{
  int err = 0;
  alloc_sprintf(restart_path_str,"%s/restart/STEP_%.5ld/cell_%.5ld",opath,step,cell_id);
  err += make_path(*restart_path_str,DIR_MODE);
  return (err == -1); /* make_path returns -1 on error */
}

/**
 * Generate the restart filename given the cell id and step.
 *
 * User must ensure that restart_fname is available for allocation
 * upon entry. User is responsible for releasing memory pointed to
 * restart_fname upon exit. rank refers to the MPI rank on the
 * macro/micro communicators, NOT rank on MPI_COMM_WORLD. \return
 * non-zero on error.
 */
static int generate_restart_fname(const char *opath,
				  const size_t rank,
				  const size_t cell_id,
				  const size_t step,
				  char **restart_fname)
{
  int err = 0;
  char *path = NULL;
  err += generate_restart_path(opath,cell_id,step,&path);
  alloc_sprintf(restart_fname,"%s/restart.%.4ld.%ld.%ld",path,rank,cell_id,step);
  free(path);
  return err;
}

int pgf_FE2_restart_print_macro(Macroscale *macro)
{
  int err = 0;
  int rank = macro->rank;

  /* update dt from any macroscale subdivision. */
  {
    int tim = macro->sol->tim;
    const double *times = macro->sol->times;
    macro->sol->dt = times[tim + 1] - times[tim];
  }

  /* The current implementation of the FE2 solver does ont use the
     state vector for updating/resetting the solution at the
     macroscale. To populate the state vector with current values we
     need to call update_M*CROSCALE_SOLUTION */
  err += update_MULTISCALE_SOLUTION(macro->sol,macro);

  /* generate the restart filename and write the file */
  char *restart_fname = NULL;
  err += generate_restart_fname(macro->opts->opath,rank,0,
				macro->sol->tim,&restart_fname);
  FILE *out = PGFEM_fopen(restart_fname,"w");
  err += dump_MULTISCALE_SOLUTION_state(macro->sol,out);
  PGFEM_fclose(out);
  free(restart_fname);

  return err;
}

int pgf_FE2_restart_print_micro(const Microscale *micro,
				const size_t cell_id)
{
  int err = 0;

  /* get the index to the solution for cell_id */
  int idx = -1; /* poison */
  idx = sol_idx_map_id_get_idx(&(micro->idx_map),cell_id);
  assert(idx >= 0);

  /* get the MPI rank on the microscale communicator */
  int rank = micro->rank;

  /* alias to the correct solution */
  const MULTISCALE_SOLUTION *s = micro->sol + idx;

  /* generate the restart filename and write the file */
  char *restart_fname = NULL;
  err += generate_restart_fname(micro->opts->opath,rank,cell_id,s->p_tim,
				&restart_fname);


  FILE *out = PGFEM_fopen(restart_fname,"w");
  err += dump_MULTISCALE_SOLUTION_state(s,out);
  PGFEM_fclose(out);
  free(restart_fname);

  return err;
}

int pgf_FE2_restart_read_macro(Macroscale *macro,
			       const size_t step,
			       const int mp_id)
{
  int err = 0;

  /* get rank on macroscale communicator */
  int rank = macro->rank;

  /* generate the restart filename and read the file */
  char *restart_fname = NULL;
  err += generate_restart_fname(macro->opts->opath,rank,0,
				step,&restart_fname);
  FILE *in = PGFEM_fopen(restart_fname,"r");
  err += read_MULTISCALE_SOLUTION_state(macro->sol,in);
  PGFEM_fclose(in);
  free(restart_fname);

  /* need to reset the M*CROSCALE_SOLUTION because current
     implementation of the FE2 solver does not typically use the
     packed vector at the macroscale. */
  reset_MULTISCALE_SOLUTION(macro->sol,macro);

  /* aliases */
  MultiscaleCommon *c = macro;
  MULTISCALE_SOLUTION *s = macro->sol;

  /* set time increment to be consistent w/ previous solve for
     substeping */
  s->times[step] = s->times[step+1] - s->dt;

  /* push r -> d_r and r <- 0 */
  memcpy(s->d_r,s->r,c->ndofd*sizeof(double));
  memset(s->r,0,c->ndofd*sizeof(double));

  /* increment the solution */
  double pores = 0;
  double nor_min = 1e-5;
  if(macro->opts->cohesive){
    increment_cohesive_elements(c->nce,c->coel,&pores,
				c->node,c->supports,s->d_r,mp_id);
  }

  /* Finite deformations increment */
  switch(macro->opts->analysis_type){
  case FS_CRPL:
  case FINITE_STRAIN:
    fd_increment (c->ne,c->nn,c->ndofn,c->npres,c->matgeom,c->hommat,
		  c->elem,c->node,c->supports,s->eps,s->sig_e,s->d_r,s->r,
		  nor_min,s->crpl,s->dt,c->nce,c->coel,&pores,c,
		  c->VVolume,macro->opts,mp_id);
    break;
  case STABILIZED:
    st_increment (c->ne,c->nn,c->ndofn,c->ndofd,c->matgeom,c->hommat,
		  c->elem,c->node,c->supports,s->eps,s->sig_e,s->d_r,s->r,
		  nor_min,macro->opts->stab,s->dt,c->nce,c->coel,&pores,
		  c,macro->opts->cohesive,mp_id);
    break;
  case MINI:
    MINI_increment(c->elem,c->ne,c->node,c->nn,c->ndofn,
		   c->supports,s->eps,s->sig_e,c->hommat,s->d_r,c,mp_id);
    break;
  case MINI_3F:
    MINI_3f_increment(c->elem,c->ne,c->node,c->nn,c->ndofn,
		      c->supports,s->eps,s->sig_e,c->hommat,s->d_r,c,mp_id);
    break;
  case DISP:
    DISP_increment(c->elem,c->ne,c->node,c->nn,c->ndofn,c->supports,s->eps,
		   s->sig_e,c->hommat,s->d_r,s->r,c,mp_id);
    break;
  default: break;
  }

  /* set displacement vector and clear increment */
  memcpy(s->r,s->d_r,c->ndofd*sizeof(double));
  memset(s->d_r,0,c->ndofd*sizeof(double));

  /* increment the supports */
  memcpy(c->supports->defl,c->supports->defl_d,c->supports->npd * sizeof(double));
  memset(c->supports->defl_d,0,c->supports->npd * sizeof(double));

  return err;
}

int pgf_FE2_restart_read_micro(Microscale *micro,
			       const size_t step,
			       const size_t cell_id)
{
  int err = 0;

  /* get the index to the solution for cell_id */
  int idx = -1; /* poison value */
  idx = sol_idx_map_id_get_idx(&(micro->idx_map),cell_id);
  assert(idx >= 0);
  MULTISCALE_SOLUTION *s = micro->sol + idx; /* alias */
  
  /* get rank on microscale communicator */
  int rank = micro->rank;

  /* generate the restart filename and read the file */
  char *restart_fname = NULL;
  err += generate_restart_fname(micro->opts->opath,rank,cell_id,
				step,&restart_fname);

#ifndef NDEBUG
  if(rank == 0){
    PGFEM_printerr("Reading restart file: %s\n",restart_fname);
  }
#endif

  FILE *in = PGFEM_fopen(restart_fname,"r");
  err += read_MULTISCALE_SOLUTION_state(s,in);
  PGFEM_fclose(in);
  free(restart_fname);

  /* Do not need to reset the solution because this is done at the
     beginning of every FE2 job at the microscale. */

  return err;
}

