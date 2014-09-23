/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "pgf_fe2_macro_client.h"
#include "PGFEM_mpi.h"
#include "microscale_information.h"
#include <stdlib.h>

/* fully define the macro client structure */
struct pgf_FE2_macro_client{
  size_t n_jobs;
  /* macroscale job(s) */
  /* communication buffers */
};

void pgf_FE2_macro_client_init(pgf_FE2_macro_client **client)
{
  *client = malloc(sizeof(**client));
  /* other initialization stuff */
}

void pgf_FE2_macro_client_destroy(pgf_FE2_macro_client *client)
{
  /* destroy internal objects */
  free(client);
  client = NULL;
}

void pgf_FE2_macro_client_create_job_list(pgf_FE2_macro_client *client
					  /* TBD */)
{
  /* create list of jobs */
}

void pgf_FE2_macro_client_assign_initial_servers(pgf_FE2_macro_client *client
						 /* TBD */)
{

}

void pgf_FE2_macro_client_rebalance_servers(pgf_FE2_macro_client *client
					    /* TBD */)
{

}

void pgf_FE2_macro_client_send_jobs(pgf_FE2_macro_client *client
				    /* TBD */)
{

}

void pgf_FE2_macro_client_recv_jobs(pgf_FE2_macro_client *client
				    /* TBD */)
{

}

void pgf_FE2_macro_client_send_exit(pgf_FE2_macro_client *client
				     /* TBD */)
{

}

