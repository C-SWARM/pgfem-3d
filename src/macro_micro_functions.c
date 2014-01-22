#include "macro_micro_functions.h"

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef COMPUTE_MS_COHE_JOB_H
#include "compute_ms_cohe_job.h"
#endif

#ifndef PLOC_SPARSE_H
#include "PLoc_Sparse.h"
#endif

#ifndef QUADRATURE_RULES_H
#include "quadrature_rules.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef COHESIVE_ELEMENT_UTILS_H
#include "cohesive_element_utils.h"
#endif

#ifndef STIFFMAT_FD_H
#include "stiffmat_fd.h"
#endif

static const int ndim = 3;

int compute_n_job_and_job_sizes(const COMMON_MACROSCALE *c,
				int *n_jobs,
				int **job_buff_sizes)
{
  int err = 0;

  /* initialize */
  *n_jobs = 0;
  *job_buff_sizes = NULL;

  /* currently only compute jobs for cohesive elements */
  *n_jobs = c->nce;
  *job_buff_sizes = NULL;
  if(*n_jobs > 0)*job_buff_sizes = PGFEM_calloc(*n_jobs,sizeof(int));
  int *buff_sizes = *job_buff_sizes; /*alias */
  MS_COHE_JOB_INFO *tmp = PGFEM_calloc(1,sizeof(MS_COHE_JOB_INFO));
  for(int i=0; i<c->nce; i++){
    /* intitialize JOB_INFO */
    err += build_MS_COHE_JOB_INFO(tmp, c->coel[i].toe);

    /* compute size */
    buff_sizes[i] = compute_MS_COHE_JOB_INFO_size(tmp);

    /* destroy JOB_INFO */
    destroy_MS_COHE_JOB_INFO(tmp);
  }

  /* cleanup */
  free(tmp);
  return err;
}

int start_microscale_server(const PGFEM_mpi_comm *mpi_comm,
			    const PGFEM_ms_job_intercomm *ic,
			    MICROSCALE *microscale)
{
  int err = 0;
  /* error check */
  if(!mpi_comm->valid_micro) return ++err;
  if(mpi_comm->valid_mm_inter && ic == NULL) return ++err;

  /* begin serving */
  if(mpi_comm->valid_mm_inter){
    /* set up server contexts */
    PGFEM_server_ctx *send = PGFEM_calloc(1,sizeof(PGFEM_server_ctx));
    PGFEM_server_ctx *recv = PGFEM_calloc(1,sizeof(PGFEM_server_ctx));
    err += initialize_PGFEM_server_ctx(send);
    err += initialize_PGFEM_server_ctx(recv);
    err += build_PGFEM_server_ctx_from_PGFEM_comm_info(ic->send_info,send);
    err += build_PGFEM_server_ctx_from_PGFEM_comm_info(ic->recv_info,recv);

    /* post receives */
    for(int i=0; i<recv->n_comms; i++){
      err += MPI_Irecv(recv->buffer[i],recv->sizes[i],MPI_CHAR,
		       recv->procs[i],MPI_ANY_TAG,ic->comm,
		       &(recv->req[i]));
    }

    /* server loop */
    int exit_server = 0;
    while(!exit_server){
      int idx = 0;
      /* wait for ANY job to show up */
      err += MPI_Waitany(recv->n_comms,recv->req,
			 &idx,recv->stat);

      /* got job, compute job */
      err += micro_job_master(mpi_comm,idx,recv->sizes[idx],
			      recv->buffer[idx],send->buffer[idx],
			      microscale,&exit_server);

      /* wait for any previous send of this process to complete. May
	 need to do some status checking here to avoid blocking in the
	 future */
      err += MPI_Wait(&(send->req[idx]),&(send->stat[idx]));
      
      /* post send of job */
      err += MPI_Isend(send->buffer[idx],send->sizes[idx],MPI_CHAR,
		       send->procs[idx],idx /*tag*/,ic->comm,
		       &(send->req[idx]));

      /* post recv for next job */
      err += MPI_Irecv(recv->buffer[idx],recv->sizes[idx],MPI_CHAR,
		       recv->procs[idx],MPI_ANY_TAG,ic->comm,
		       &(recv->req[idx]));
    }

    /* cancel any pending communication */
    for(int i=0; i<recv->n_comms; i++){
      err += MPI_Cancel(&(recv->req[i]));
    }
    for(int i=0; i<send->n_comms; i++){
      err += MPI_Cancel(&(send->req[i]));
    }

    /* cleanup */
    err += destroy_PGFEM_server_ctx(send);
    err += destroy_PGFEM_server_ctx(recv);
    free(send);
    free(recv);

  } else {
   err += micro_job_slave(mpi_comm,microscale);
  }

  return err;
}
			    
int micro_job_master(const PGFEM_mpi_comm *mpi_comm,
		     const int idx,
		     const int buff_size,
		     char *in_buffer,
		     char *out_buffer,
		     MICROSCALE *micro,
		     int *exit_server)
{
  int err = 0;
  /* exit if not master on microscale */
  if(!mpi_comm->valid_micro || mpi_comm->rank_micro != 0) return ++err;

  /* broadcast job information */
  int job_init_info[2] = {0,0};
  job_init_info[0] = idx;
  job_init_info[1] = buff_size;
  err += MPI_Bcast(job_init_info,2,MPI_INT,0,mpi_comm->micro);
  err += MPI_Bcast(in_buffer,buff_size,MPI_CHAR,0,mpi_comm->micro);

  /* copy the in buffer to the out buffer */
  memcpy(out_buffer,in_buffer,buff_size*sizeof(char));

  /* compute the job */
  err += microscale_compute_job(idx,buff_size,
				out_buffer,micro,
				exit_server);
  return err;
}

int micro_job_slave(const PGFEM_mpi_comm *mpi_comm,
		    MICROSCALE *micro)
{
  int err = 0;
  int exit_server = 0;

  /* exit if not slave on microscale */
  if(!mpi_comm->valid_micro && mpi_comm->rank_micro <= 0) return ++err;

  while(!exit_server){
    /* get job information from master */
    int job_init_info[2] = {0,0};
    err += MPI_Bcast(job_init_info,2,MPI_INT,0,mpi_comm->micro);
    char *buffer = PGFEM_calloc(job_init_info[1],sizeof(char));
    err += MPI_Bcast(buffer,job_init_info[1],MPI_CHAR,0,mpi_comm->micro);

    /* compute the job */
    err += microscale_compute_job(job_init_info[0],
				  job_init_info[1],
				  buffer,micro,&exit_server);

    /* cleanup */
    free(buffer);
  }
  return err;
}

int microscale_compute_job(const int idx,
			   const int buff_len,
			   char *buffer,
			   MICROSCALE *micro,
			   int *exit_server)
{
  int err = 0;

  /* construct the job info from the buffer */
  MS_COHE_JOB_INFO *job = PGFEM_calloc(1,sizeof(MS_COHE_JOB_INFO));
  err += build_MS_COHE_JOB_INFO_buffer(buff_len*sizeof(char),buffer,job);

  /* check the exit_server status and exit early */
  if(job->job_type == JOB_EXIT){
    *exit_server = 1;
    memset(buffer,0,buff_len*sizeof(char));
    goto exit_function;
  }

  /* compute the job */
  err += compute_ms_cohe_job(idx,job,micro);

  /* pack the job back into buffer */
  err += pack_MS_COHE_JOB_INFO(job,buff_len,buffer);

 exit_function:
  destroy_MS_COHE_JOB_INFO(job);
  free(job);
  return err;
}

/******** Main MACROSCALE interface functions *********/
int start_macroscale_compute_jobs(const PGFEM_ms_job_intercomm *ic,
				  const MACROSCALE *macro,
				  const int job_type,
				  const double *loc_sol,
				  MS_COHE_JOB_INFO *job_list,
				  PGFEM_server_ctx *send,
				  PGFEM_server_ctx *recv)
{
  int err = 0;
  /* check that things are allocated and return controll */
  if(ic == NULL || macro == NULL) return ++err;

  /* Wait to complete any pending communcication. The error flag is
     incremented as this should only occur by a logical/programming
     error. Note that buffers may be overwritten. If the communication
     is completed but the flags have not been reset to 0, the
     MPI_Wait* commands should return immediately */
  if(send->in_process || recv->in_process){
    PGFEM_printerr("WARNING: communication in progress!(%s:%s:%d)\n",
		   __func__,__FILE__,__LINE__);
    err++;
    err += MPI_Waitall(recv->n_comms,recv->req,recv->stat);
    err += MPI_Waitall(send->n_comms,send->req,send->stat);
    recv->in_process = 0;
    send->in_process = 0;
  }


  /* post recieves (from the running server) */
  for(int i=0; i<recv->n_comms; i++){
    err += MPI_Irecv(recv->buffer[i],recv->sizes[i],MPI_CHAR,
		     recv->procs[i],MPI_ANY_TAG,ic->comm,
		     &(recv->req[i]));
  }

  for(int i=0; i<send->n_comms; i++){
    /* update the job information according to job_type */
    err += macroscale_update_job_info(macro,job_type,loc_sol,job_list+i);

    /* pack the job info into the buffer to send */
    err += pack_MS_COHE_JOB_INFO(job_list + i,send->sizes[i],
				 send->buffer[i]);

    /* post the send (to the running server) */
    err += MPI_Isend(send->buffer[i],send->sizes[i],MPI_CHAR,
		     send->procs[i],i /*tag*/,ic->comm,
		     &(send->req[i]));
  }

  /* set in_process flags to true (1) */
  recv->in_process = 1;
  send->in_process = 1;

  /* cancel communication if JOB_EXIT */
  if(job_type == JOB_EXIT){
    for(int i=0; i<recv->n_comms; i++){
      err += MPI_Cancel(&(recv->req[i]));
    }
    /* have canceled the communication, not in_process */
    recv->in_process = 0;
  }
  return err;
}

int finish_macroscale_compute_jobs(MS_COHE_JOB_INFO *job_list,
				   MACROSCALE *macro,
				   PGFEM_server_ctx *send,
				   PGFEM_server_ctx *recv)
{
  int err = 0;
  COMMON_MACROSCALE *c = macro->common;
  MACROSCALE_SOLUTION *s = macro->sol;
  int rank_macro = 0;
  int nproc_macro = 0;

  /* exit early if !*->in_process */
  if(!recv->in_process && !send->in_process) return err;

  err += MPI_Comm_rank(c->mpi_comm,&rank_macro);
  err += MPI_Comm_size(c->mpi_comm,&nproc_macro);

  /* if expecting to receive buffers */
  if(recv->in_process){
    /* set up the stiffness matrix communication */
    double **Lk = NULL;
    double **receive = NULL;
    MPI_Request *req_r = NULL;
    MPI_Status *sta_r = NULL;
    err += init_and_post_stiffmat_comm(&Lk,&receive,&req_r,&sta_r,
				       c->mpi_comm,c->pgfem_comm);

    /* assemble jobs as they are received */
    for(int i=0; i<recv->n_comms; i++){
      int idx = 0;
      err += MPI_Waitany(recv->n_comms,recv->req,&idx,recv->stat);
      MS_COHE_JOB_INFO *job = job_list + idx;
      err += unpack_MS_COHE_JOB_INFO(job,recv->sizes[idx],
				     recv->buffer[idx]);

      /* finish jobs based on job_type */
      switch(job->job_type){
      case JOB_COMPUTE_EQUILIBRIUM:
	/* assemble residual (local) */
	for(int j=0; j<job->ndofe; j++){
	  int dof_id = job->loc_dof_ids[j] - 1;
	  if(dof_id < 0) continue; /* boundary condition */
	  s->f_u[dof_id] += job->traction_res[j];
	}
	/*** Deliberate drop through ***/
      case JOB_NO_COMPUTE_EQUILIBRIUM:
	/* assemble tangent to local and off-proc buffers */
	PLoc_Sparse(NULL,Lk,job->K_00_contrib,NULL,NULL,NULL,job->g_dof_ids,
		    job->ndofe,NULL,c->GDof,rank_macro,nproc_macro,
		    c->pgfem_comm,0,c->SOLVER,macro->opts->analysis_type);
	break;
      case JOB_UPDATE:
	/* update cohesive elements */
	err += macroscale_update_coel(job,macro);
	break;
      default: /* do nothing */ break;
      }
    }

    /* send/finalize communication of the stiffness matrix */
    MPI_Status *sta_s = NULL;
    MPI_Request *req_s = NULL;
    err += send_stiffmat_comm(&sta_s,&req_s,Lk,c->mpi_comm,c->pgfem_comm);

    err += assemble_nonlocal_stiffmat(c->pgfem_comm,sta_r,req_r,
				      c->SOLVER,receive);

    err += finalize_stiffmat_comm(sta_s,sta_r,req_s,req_r,c->pgfem_comm);

    /* re-initialize preconditioner ? */
    /* err += PGFEM_HYPRE_create_preconditioner(c->SOLVER,c->mpi_comm); */

    /* clean up memory */
    for(int i=0; i<nproc_macro; i++){
      if(Lk != NULL) free(Lk[i]);
      if(receive != NULL) free(receive[i]);
    }
    free(Lk);
    free(receive);
    free(sta_r);
    free(req_r);
    free(sta_s);
    free(req_s);

    /* set in_process to false (0) */
    recv->in_process = 0;
  }/* end if(recv->in_process) */

  /* wait for any pending send communication to finish */
  if(send->in_process){
    MPI_Waitall(send->n_comms,send->req,send->stat);
    send->in_process = 0;
  }

  return err;
}

int macroscale_update_job_info(const MACROSCALE *macro,
			       const int job_type,
			       const double *loc_sol,
			       MS_COHE_JOB_INFO *job)
{
  int err = 0;
  const COMMON_MACROSCALE *c = macro->common;
  const MACROSCALE_SOLUTION *s = macro->sol;
  const COEL *cel = c->coel + job->elem_id;
  const int nne = cel->toe;
  const int nne_2D = nne/2;
  double *x = PGFEM_calloc(nne,sizeof(double));
  double *y = PGFEM_calloc(nne,sizeof(double));
  double *z = PGFEM_calloc(nne,sizeof(double));
  double *disp = PGFEM_calloc(nne*ndim,sizeof(double));

  /* get nodal coordinates and displacements on the element */
  nodecoord_updated(nne,cel->nod,c->node,x,y,z);
  def_elem(job->loc_dof_ids,nne*ndim,loc_sol,
	   NULL,c->node,disp,c->supports,0);

  /* set up integration SINGLE INTEGRATION POINT */
  double *gk = NULL;
  double *ge = NULL;
  double *w = NULL;
  int n_ip = 0;
  switch(nne_2D){
  case 3: err += get_tria_quadrature_rule(0,&n_ip,&gk,&ge,&w); break;
  case 4: err += get_quad_quadrature_rule(0,&n_ip,&gk,&ge,&w); break;
  }

  /* error checking */
  if(n_ip != 1){
    PGFEM_printerr("ERROR: currently only support single"
		   " integration point for macro-micro!\n");
    PGFEM_Abort();
  }
  
  /* set the jump */
  get_jump(nne_2D,x,y,z,disp,job->shape + nne_2D,job->jump);

  /* scale the jump */
  /* for(int i=0; i<ndim; i++){ */
  /*   /\* mm to micron *\/ */
  /*   job->jump[i] *= 1000; */
  /* } */

  /* set the job type */
  job->job_type = job_type;

  /* clear the tangent and res */
  memset(job->traction_res,0,(job->ndofe)*sizeof(double));
  memset(job->K_00_contrib,0,(job->ndofe)*(job->ndofe)*sizeof(double));

  /* set job->tim and job->times */
  job->tim = s->tim;
  if(s->tim == 0){
    memcpy(job->times,s->times,3*sizeof(double));
  } else {
    memcpy(job->times,s->times + (s->tim - 1),3*sizeof(double));
  }

  free(x);
  free(y);
  free(z);
  free(disp);
  free(gk);
  free(ge);
  free(w);

  return err;
}

int macroscale_update_coel(const MS_COHE_JOB_INFO *job,
			   MACROSCALE *macro)
{
  int err = 0;
  COEL *coel = macro->common->coel + job->elem_id;

  /* set macroscale traction on the element from the microscale job */
  memcpy(coel->ti,job->traction,ndim*sizeof(double));

  /* That's all folks! */
  return err;
}
