/* HEADER */
/**
 * AUTHORS:
 *  Matthew Mosby, University of Notre Dame <mmosby1[at]nd.edu>
 */

#include "pgf_fe2_compute_max_n_jobs.h"
#include "PGFEM_io.h"
#include <assert.h>
#include <math.h>

int pgf_FE2_compute_max_n_jobs(const MACROSCALE *macro,
			       const PGFEM_mpi_comm *mpi_comm,
			       int *max_n_jobs)
{
  int err = 0;
  static const double factor = 1.1;
  if ( macro == NULL ){ /* microscale server(s) */

    /* Matching call to Bcast from Macroscale */
    err += MPI_Bcast(max_n_jobs,1,MPI_INT,0,mpi_comm->world);
    assert(*max_n_jobs > 0);

  } else { /* Macroscale */

    /* get rank and sizes */
    int rank = -1;
    int n_macro = -1;
    int n_server = -1;
    err += MPI_Comm_rank(mpi_comm->world,&rank);
    err += MPI_Comm_size(mpi_comm->mm_inter,&n_server);
    err += MPI_Comm_size(mpi_comm->macro,&n_macro);
    n_server -= n_macro;
    assert(n_server > 0);

    /* Get number of cohesive elements on macro domain */
    int total_macro_nce = macro->common->nce;

    if(rank == 0){
      /* We use MPI_Reduce rather than MPI_Allreduce since we will
	 call MPI_Bcast from rank 0 anyway. */
      err += MPI_Reduce(MPI_IN_PLACE,&total_macro_nce,1,
			MPI_INT,MPI_SUM,0,mpi_comm->macro);

      /* compute minimum number of jobs */
      int min_n_jobs = total_macro_nce / n_server;
      assert( min_n_jobs > 0 );
      if ( (total_macro_nce % n_server) > 0 ) min_n_jobs ++;

      /* compute/override maximum number of jobs */
      *max_n_jobs = (int) ceil(factor * min_n_jobs);
      if (!macro->opts->no_migrate
	  && macro->opts->max_n_jobs > 0 ){
	if( macro->opts->max_n_jobs >= min_n_jobs ) {
	  *max_n_jobs = macro->opts->max_n_jobs;
	} else {
	  PGFEM_printerr("WARNING: specified max. server jobs too small! Using default.\n");
	}
      }

      /* if we are not migrating cells, use min */
      if(macro->opts->no_migrate) *max_n_jobs = min_n_jobs;

      PGFEM_printf("Max. No. Jobs/Server: %d\n",*max_n_jobs);
    } else {
      /* recv buffer for reduce. Contents invalid after Reduce */
      int junk = -1;
      err += MPI_Reduce(&total_macro_nce,&junk,1,
			MPI_INT,MPI_SUM,0,mpi_comm->macro);
    }

    /* broadcacst on macroscale processes in world communicator */
    err += MPI_Bcast(max_n_jobs,1,MPI_INT,0,mpi_comm->world);
  }
  return err;
}
