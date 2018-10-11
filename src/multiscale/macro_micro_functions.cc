#include "macro_micro_functions.h"
#include "allocation.h"
#include "compute_ms_cohe_job.h"
#include "PLoc_Sparse.h"
#include "quadrature_rules.h"
#include "utils.h"
#include "cohesive_element_utils.h"
#include "stiffmat_fd.h"

using namespace pgfem3d;
using namespace pgfem3d::net;

static const int ndim = 3;

int compute_n_job_and_job_sizes(const MultiscaleCommon *c,
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
  if(*n_jobs > 0)*job_buff_sizes = PGFEM_calloc(int, *n_jobs);
  int *buff_sizes = *job_buff_sizes; /*alias */
  MS_COHE_JOB_INFO *tmp = PGFEM_calloc(MS_COHE_JOB_INFO, 1);
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

int start_microscale_server(const MultiscaleComm *mscom,//deprecated
                const PGFEM_ms_job_intercomm *ic,
                Microscale *microscale,
                const int mp_id)
{
  ISIRNetwork *net = static_cast<ISIRNetwork*>(microscale->net);
  int err = 0;
  /* error check */
  if(!mscom->valid_micro_1) return ++err;
  if(mscom->valid_mm_inter && ic == NULL) return ++err;

  /* begin serving */
  if(mscom->valid_mm_inter){
    /* set up server contexts */
    MultiscaleServerContext *send = new MultiscaleServerContext(net);
    MultiscaleServerContext *recv = new MultiscaleServerContext(net);
    send->initialize(ic->send_info);
    recv->initialize(ic->recv_info);

    /* post receives */
    for(int i=0; i<recv->n_comms; i++){
      net->irecv(recv->buffer[i],recv->sizes[i],NET_DT_CHAR,
		 recv->procs[i],NET_ANY_TAG,ic->comm,
		 &(recv->req[i]));
    }

    /* server loop */
    int exit_server = 0;
    while(!exit_server){
      int idx = 0;
      /* wait for ANY job to show up */
      net->waitany(recv->n_comms,recv->req,
		   &idx,recv->stat);

      /* got job, compute job */
      err += micro_job_master(mscom,idx,recv->sizes[idx],
			      recv->buffer[idx],send->buffer[idx],
			      microscale,&exit_server,mp_id);

      /* wait for any previous send of this process to complete. May
     need to do some status checking here to avoid blocking in the
     future */
      net->wait(&(send->req[idx]),&(send->stat[idx]));

      /* post send of job */
      net->isend(send->buffer[idx],send->sizes[idx],NET_DT_CHAR,
		 send->procs[idx],idx /*tag*/,ic->comm,
		 &(send->req[idx]));

      /* post recv for next job */
      net->irecv(recv->buffer[idx],recv->sizes[idx],NET_DT_CHAR,
		 recv->procs[idx],NET_ANY_TAG,ic->comm,
		 &(recv->req[idx]));
    }

    /* cancel any pending communication */
    for(int i=0; i<recv->n_comms; i++){
      net->cancel(&(recv->req[i]));
    }
    for(int i=0; i<send->n_comms; i++){
      if(send->req[i].getData() != NULL) {
	net->cancel(&(send->req[i]));
      }
    }
    net->waitall(recv->n_comms,recv->req,NET_STATUS_IGNORE);
    net->waitall(send->n_comms,send->req,NET_STATUS_IGNORE);

    /* cleanup */
    delete send;
    delete recv;
  } else {
    err += micro_job_slave(mscom,microscale,mp_id);
  }

  return err;
}

int micro_job_master(const MultiscaleComm *mscom, //deprecated
		     const int idx,
		     const int buff_size,
		     char *in_buffer,
		     char *out_buffer,
		     Microscale *micro,
		     int *exit_server,
		     const int mp_id)
{
  int err = 0;
  /* exit if not master on microscale */
<<<<<<< HEAD
<<<<<<< HEAD
  if(!mscom->valid_micro || mscom->rank_micro != 0) return ++err;

=======
 int micro_model = 1;//this is unused code anyways
 if(!mscom->valid_micro_1 || mscom->rank_micro != 0) return ++err;
>>>>>>> fixed all mpi/memory errors. beginning tests
=======
 int micro_model = 1;//this is unused code anyways
 if(!mscom->valid_micro_1 || mscom->rank_micro != 0) return ++err;
>>>>>>> 56768dcd05fd9525ebaf1db9f08b1889daef65bf
  /* broadcast job information */
  int job_init_info[2] = {0,0};
  job_init_info[0] = idx;
  job_init_info[1] = buff_size;
  micro->net->bcast(job_init_info,2,NET_DT_INT,0,mscom->micro);
  micro->net->bcast(in_buffer,buff_size,NET_DT_CHAR,0,mscom->micro);

  /* copy the in buffer to the out buffer */
  memcpy(out_buffer,in_buffer,buff_size*sizeof(char));

  /* compute the job */
  err += microscale_compute_job(idx,buff_size,
				out_buffer,micro,
				exit_server,mp_id,micro_model);
  return err;
}

int micro_job_slave(const MultiscaleComm *mscom, //deprecated
		    Microscale *micro,const int mp_id)
{
  int err = 0;
  int exit_server = 0;
  int micro_model = 1;//this is unused code anyways
  /* exit if not slave on microscale */
  if(!mscom->valid_micro_1 && mscom->rank_micro <= 0) return ++err;

  while(!exit_server){
    /* get job information from master */
    int job_init_info[2] = {0,0};
    micro->net->bcast(job_init_info,2,NET_DT_INT,0,mscom->micro);
    char *buffer = PGFEM_calloc(char, job_init_info[1]);
    micro->net->bcast(buffer,job_init_info[1],NET_DT_CHAR,0,mscom->micro);
    
    /* compute the job */
    err += microscale_compute_job(job_init_info[0],
				  job_init_info[1],
				  buffer,micro,&exit_server,mp_id,micro_model);

    /* cleanup */
    free(buffer);
  }
  return err;
}

int microscale_compute_job(const int idx,
               const int buff_len,
               char *buffer,
               Microscale *micro,
               int *exit_server,
               const int mp_id,
                int micro_model)
{
  int err = 0;

  /* construct the job info from the buffer */
  MS_COHE_JOB_INFO *job = PGFEM_calloc(MS_COHE_JOB_INFO, 1);
  err += build_MS_COHE_JOB_INFO_buffer(buff_len*sizeof(char),buffer,job);

  /* check the exit_server status and exit early */
  if(job->job_type == JOB_EXIT){
    *exit_server = 1;
    memset(buffer,0,buff_len*sizeof(char));
    goto exit_function;
  }

  /* compute the job */
  err += compute_ms_cohe_job(idx,job,micro,mp_id,micro_model);

  /* pack the job back into buffer */
  err += pack_MS_COHE_JOB_INFO(job,buff_len,buffer);

 exit_function:
  destroy_MS_COHE_JOB_INFO(job);
  free(job);
  return err;
}

/******** Main Macroscale interface functions *********/
int start_macroscale_compute_jobs(const PGFEM_ms_job_intercomm *ic,//deprecated
				  const Macroscale *macro,
				  const int job_type,
				  const double *loc_sol,
				  MS_COHE_JOB_INFO *job_list,
				  MultiscaleServerContext *send,
				  MultiscaleServerContext *recv)
{
  int err = 0;
  ISIRNetwork *net = static_cast<ISIRNetwork*>(macro->net);
  /* check that things are allocated and return controll */
  if(ic == NULL || macro == NULL) return ++err;

  /* Wait to complete any pending communcication. The error flag is
     incremented as this should only occur by a logical/programming
     error. Note that buffers may be overwritten. If the communication
     is completed but the flags have not been reset to 0, the
     NET_Wait* commands should return immediately */
  if(send->in_process || recv->in_process){
    PGFEM_printerr("WARNING: communication in progress!(%s:%s:%d)\n",
           __func__,__FILE__,__LINE__);
    err++;
    net->waitall(recv->n_comms,recv->req,recv->stat);
    net->waitall(send->n_comms,send->req,send->stat);
    recv->in_process = 0;
    send->in_process = 0;
  }


  /* post recieves (from the running server) */
  for(int i=0; i<recv->n_comms; i++){
    net->irecv(recv->buffer[i],recv->sizes[i],NET_DT_CHAR,
	       recv->procs[i],NET_ANY_TAG,ic->comm,
	       &(recv->req[i]));
  }

  for(int i=0; i<send->n_comms; i++){
    /* update the job information according to job_type */
    err += macroscale_update_job_info(macro,job_type,loc_sol,job_list+i);

    /* pack the job info into the buffer to send */
    err += pack_MS_COHE_JOB_INFO(job_list + i,send->sizes[i],
                 send->buffer[i]);

    /* post the send (to the running server) */
    net->isend(send->buffer[i],send->sizes[i],NET_DT_CHAR,
	       send->procs[i],i /*tag*/,ic->comm,
	       &(send->req[i]));
  }

  /* set in_process flags to true (1) */
  recv->in_process = 1;
  send->in_process = 1;

  /* cancel communication if JOB_EXIT */
  if(job_type == JOB_EXIT){
    for(int i=0; i<recv->n_comms; i++){
      net->cancel(&(recv->req[i]));
    }
    /* have canceled the communication, not in_process */
    recv->in_process = 0;
  }
  return err;
}

int finish_macroscale_compute_jobs(MS_COHE_JOB_INFO *job_list,//deprecated
				   Macroscale *macro,
				   MultiscaleServerContext *send,
				   MultiscaleServerContext *recv)
{
  int err = 0;
  ISIRNetwork *net = static_cast<ISIRNetwork*>(macro->net);
  MultiscaleCommon *c = macro;
  MULTISCALE_SOLUTION *s = macro->sol;
  int rank_macro = macro->rank;
  int nproc_macro = macro->nproc;

  /* exit early if !*->in_process */
  if(!recv->in_process && !send->in_process) return err;
  
  /* if expecting to receive buffers */
  if(recv->in_process){
    /* set up the stiffness matrix communication */
    double **Lk = NULL;
    double **receive = NULL;
    c->spc->post_stiffmat(&Lk,&receive);
    
    /* assemble jobs as they are received */
    for(int i=0; i<recv->n_comms; i++){
      int idx = 0;
      net->waitany(recv->n_comms,recv->req,&idx,recv->stat);
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
	PLoc_Sparse(Lk,job->K_00_contrib,NULL,NULL,NULL,job->g_dof_ids,
		    job->ndofe,NULL,c->GDof,rank_macro,nproc_macro,
		    c->spc,0,c->SOLVER,macro->opts->analysis_type);
	break;
      case JOB_UPDATE:
	/* update cohesive elements */
	err += macroscale_update_coel(job,macro);
	break;
      default: /* do nothing */ break;
      }
    }
    
    /* send/finalize communication of the stiffness matrix */
    c->spc->send_stiffmat();
    c->spc->assemble_nonlocal_stiffmat(c->SOLVER);
    c->spc->finalize_stiffmat();
    
    /* re-initialize preconditioner ? */
    /* err += PGFEM_HYPRE_create_preconditioner(c->SOLVER,c->mscom); */

    /* set in_process to false (0) */
    recv->in_process = 0;
  }/* end if(recv->in_process) */

  /* wait for any pending send communication to finish */
  if(send->in_process){
    net->waitall(send->n_comms,send->req,send->stat);
    send->in_process = 0;
  }
  
  return err;
}

int macroscale_update_job_info(const Macroscale *macro,
                   const int job_type,
                   const double *loc_sol,
                   MS_COHE_JOB_INFO *job)
{
  int err = 0;
  const MultiscaleCommon *c = macro;
  const MULTISCALE_SOLUTION *s = macro->sol;
  const COEL *cel = c->coel + job->global_job_id;
  const int nne = cel->toe;
  const int nne_2D = nne/2;
  double *x = PGFEM_calloc(double, nne);
  double *y = PGFEM_calloc(double, nne);
  double *z = PGFEM_calloc(double, nne);
  double *disp = PGFEM_calloc(double, nne*ndim);

  /* get nodal coordinates and displacements on the element */
  nodecoord_updated(nne,cel->nod,c->node,x,y,z);
  def_elem(job->loc_dof_ids,nne*ndim,loc_sol,
       NULL,c->node,disp,c->supports,0);

  /* set up integration SINGLE INTEGRATION POINT */
  double *gk = NULL;
  double *ge = NULL;
  double *w = NULL;
  int n_ip = 0;
  int ip = 0;
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

  if(macro->opts->restart >= 0){
    memcpy(job->jump_n,job->jump,ndim*sizeof(double));
    job->max_traction = cel->vars[ip][0];
    job->max_jump = cel->vars[ip][1];
    memcpy(job->traction_n,cel->vars[ip] + 2,ndim*sizeof(double));
  }

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
               Macroscale *macro)
{
  int err = 0;
  COEL *coel = macro->coel + job->global_job_id;  //this particular macroscale cohesive element

  /* set traction for output */
  memcpy(coel->ti,job->traction_n,ndim*sizeof(double));

  /* set the state variable(s) */
  coel->vars[job->int_pt][0] = job->max_traction;
  coel->vars[job->int_pt][1] = job->max_jump;
  memcpy(coel->vars[job->int_pt] + 2,job->traction_n,ndim*sizeof(double));

  /* That's all folks! */
  return err;
}
