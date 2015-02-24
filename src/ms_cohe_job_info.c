/* HEADER */
#include "ms_cohe_job_info.h"
#include <string.h>
#include "allocation.h"
#include "utils.h"
#include "pgf_fe2_job.h"

#ifndef DEBUG_MS_JOB_INFO
#define DEBUG_MS_JOB_INFO 0
#endif

static const int ndim = 3;

inline size_t compute_MS_COHE_JOB_INFO_size(const MS_COHE_JOB_INFO *info)
{
  size_t result = (7*sizeof(int)
		   + sizeof(double)
		   + 5*ndim*sizeof(double)
		   + sizeof(double)
		   + info->nnode*sizeof(double)
		   + info->ndofe*sizeof(double)
		   + (info->ndofe*info->ndofe)*sizeof(double)
		   + 2*info->ndofe*sizeof(long)
		   + sizeof(int)
		   + 3*sizeof(double)
		   + sizeof(int));
  return result;
}	 

int build_MS_COHE_JOB_INFO(MS_COHE_JOB_INFO *info,
			   const int nnode)
{
  int err = 0;
  /* set variables */
  info->nnode = nnode;
  info->ndofe = nnode*ndim;
  info->elem_id = -1; /* poisoned value */
  info->proc_id = -1; /* poisoned value */
  info->int_pt = -1; /* poisoned value */
  info->job_type = -1; /* poisoned value */
  info->print_flag = 0; /* print output */

  info->int_wt = 0.0;

  /* allocate len = ndim */
  info->jump = PGFEM_calloc(ndim,sizeof(double));
  info->jump_n = PGFEM_calloc(ndim,sizeof(double));
  info->normal = PGFEM_calloc(ndim,sizeof(double));
  info->traction = PGFEM_calloc(ndim,sizeof(double));
  info->traction_n = PGFEM_calloc(ndim,sizeof(double));

  /* allocate len = nnode */
  info->shape = PGFEM_calloc(info->nnode,sizeof(double));

  /* allocate len = ndofe */
  info->traction_res = PGFEM_calloc(info->ndofe,sizeof(double));
  info->loc_dof_ids = PGFEM_calloc(info->ndofe,sizeof(double));
  info->g_dof_ids = PGFEM_calloc(info->ndofe,sizeof(double));

  /* allocate len = ndofe*ndofe */
  info->K_00_contrib = PGFEM_calloc(info->ndofe*info->ndofe,sizeof(double));

  /* allocate other */
  info->times = PGFEM_calloc(3,sizeof(double));
  info->n_step = 1;

  return err;
}

int build_MS_COHE_JOB_INFO_buffer(const size_t buff_len,
				  const char *buffer,
				  MS_COHE_JOB_INFO *info)
{
  int err = 0;

  /* unpack the number of nodes from the buffer */
  size_t pos = 0;
  unpack_data(buffer,&info->nnode,&pos,1,sizeof(int));

  /* allocate the object */
  err += build_MS_COHE_JOB_INFO(info,info->nnode);

  /* unpack information from the buffer */
  err += unpack_MS_COHE_JOB_INFO(info,buff_len,buffer);

  return err;
}

int set_MS_COHE_JOB_INFO(MS_COHE_JOB_INFO *info,
			 const double *normal,
			 const double *jump,
			 const double *shape,
			 const long *loc_dof_ids,
			 const long *g_dof_ids)
{
  int err = 0;
  memcpy(info->jump,jump,ndim*sizeof(double));
  memcpy(info->normal,normal,ndim*sizeof(double));
  memcpy(info->shape,shape,info->nnode*sizeof(double));
  memcpy(info->loc_dof_ids,loc_dof_ids,info->nnode*ndim*sizeof(long));
  memcpy(info->g_dof_ids,g_dof_ids,info->nnode*ndim*sizeof(long));
  return err;
}

void destroy_MS_COHE_JOB_INFO(MS_COHE_JOB_INFO *info)
{
  if(info != NULL){
    free(info->jump);
    free(info->jump_n);
    free(info->normal);
    free(info->traction);
    free(info->traction_n);
    free(info->shape);
    free(info->traction_res);
    free(info->K_00_contrib);
    free(info->loc_dof_ids);
    free(info->g_dof_ids);
    free(info->times);
  }
}
			     
int pack_MS_COHE_JOB_INFO(const MS_COHE_JOB_INFO *info,
			  const size_t buffer_len,
			  char *buffer)
{
  int err = 0;
  size_t pos = 0;

  /* pack scalars */
  pack_data(&info->nnode,buffer,&pos,1,sizeof(int));
  pack_data(&info->ndofe,buffer,&pos,1,sizeof(int));
  pack_data(&info->elem_id,buffer,&pos,1,sizeof(int));
  pack_data(&info->proc_id,buffer,&pos,1,sizeof(int));
  pack_data(&info->int_pt,buffer,&pos,1,sizeof(int));
  pack_data(&info->job_type,buffer,&pos,1,sizeof(int));
  pack_data(&info->print_flag,buffer,&pos,1,sizeof(int));
  pack_data(&info->int_wt,buffer,&pos,1,sizeof(double));
  pack_data(&info->tim,buffer,&pos,1,sizeof(int));
  pack_data(&info->n_step,buffer,&pos,1,sizeof(int));
  pack_data(&info->max_traction,buffer,&pos,1,sizeof(double));

  /* pack arrays */
  pack_data(info->jump,buffer,&pos,ndim,sizeof(double));
  pack_data(info->jump_n,buffer,&pos,ndim,sizeof(double));
  pack_data(info->normal,buffer,&pos,ndim,sizeof(double));
  pack_data(info->traction,buffer,&pos,ndim,sizeof(double));
  pack_data(info->traction_n,buffer,&pos,ndim,sizeof(double));
  pack_data(info->shape,buffer,&pos,info->nnode,sizeof(double));
  pack_data(info->traction_res,buffer,&pos,info->ndofe,sizeof(double));
  pack_data(info->K_00_contrib,buffer,&pos,
	    info->ndofe*info->ndofe,sizeof(double));
  pack_data(info->loc_dof_ids,buffer,&pos,info->ndofe,sizeof(long));
  pack_data(info->g_dof_ids,buffer,&pos,info->ndofe,sizeof(long));
  pack_data(info->times,buffer,&pos,3,sizeof(double));

  /* error checking */
  if(pos > buffer_len) err++;
  return err;
}

int unpack_MS_COHE_JOB_INFO(MS_COHE_JOB_INFO *info,
			    const size_t buffer_len,
			    const char *buffer)
{
  int err = 0;

  /* Check that we have enough space in info */
  {
    size_t loc_buffer_len = compute_MS_COHE_JOB_INFO_size(info);
    if(loc_buffer_len < buffer_len) return (err++);
  }

  size_t pos = 0;
  /* unpack nnode and error check */
  {
    int nnode = info->nnode;
    unpack_data(buffer,&info->nnode,&pos,1,sizeof(int));
    if(nnode != info->nnode) return (err++);
  }

  /* unpack remaining scalars */
  unpack_data(buffer,&info->ndofe,&pos,1,sizeof(int));
  unpack_data(buffer,&info->elem_id,&pos,1,sizeof(int));
  unpack_data(buffer,&info->proc_id,&pos,1,sizeof(int));
  unpack_data(buffer,&info->int_pt,&pos,1,sizeof(int));
  unpack_data(buffer,&info->job_type,&pos,1,sizeof(int));
  unpack_data(buffer,&info->print_flag,&pos,1,sizeof(int));
  unpack_data(buffer,&info->int_wt,&pos,1,sizeof(double));
  unpack_data(buffer,&info->tim,&pos,1,sizeof(int));
  unpack_data(buffer,&info->n_step,&pos,1,sizeof(int));
  unpack_data(buffer,&info->max_traction,&pos,1,sizeof(double));

  /* unpack arrays */
  unpack_data(buffer,info->jump,&pos,ndim,sizeof(double));
  unpack_data(buffer,info->jump_n,&pos,ndim,sizeof(double));
  unpack_data(buffer,info->normal,&pos,ndim,sizeof(double));
  unpack_data(buffer,info->traction,&pos,ndim,sizeof(double));
  unpack_data(buffer,info->traction_n,&pos,ndim,sizeof(double));
  unpack_data(buffer,info->shape,&pos,info->nnode,sizeof(double));
  unpack_data(buffer,info->traction_res,&pos,info->ndofe,sizeof(double));
  unpack_data(buffer,info->K_00_contrib,&pos,
	      info->ndofe*info->ndofe,sizeof(double));
  unpack_data(buffer,info->loc_dof_ids,&pos,info->ndofe,sizeof(long));
  unpack_data(buffer,info->g_dof_ids,&pos,info->ndofe,sizeof(long));
  unpack_data(buffer,info->times,&pos,3,sizeof(double));

  /* error checking */
  if(pos > buffer_len) err++;
  return err;
}

static void job_type_str(const int job_type,
			 char **str)
{
  switch(job_type){
  case JOB_NO_COMPUTE_EQUILIBRIUM:
    alloc_sprintf(str,"JOB_NO_COMPUTE_EQUILIBRIUM");
    break;
  case JOB_COMPUTE_EQUILIBRIUM:
    alloc_sprintf(str,"JOB_COMPUTE_EQUILIBRIUM");
    break;
  case JOB_UPDATE:
    alloc_sprintf(str,"JOB_UPDATE");
    break;
  case JOB_PRINT:
     alloc_sprintf(str,"JOB_PRINT");
    break;
  case JOB_EXIT:
     alloc_sprintf(str,"JOB_EXIT");
    break;
  default:
     alloc_sprintf(str,"JOB_UNDEFINED");
    break;
  }
}

int print_MS_COHE_JOB_INFO(FILE *out,
			   const MS_COHE_JOB_INFO *info)
{
  int err = 0;
  const int cell_id = pgf_FE2_job_compute_encoded_id(info->proc_id,
						     info->elem_id,
						     info->int_pt);
  char *job_str = NULL;
  PGFEM_fprintf(out,"===== START JOB INFO =====\n");
  PGFEM_fprintf(out,"NNODE:   %d\n",info->nnode);
  PGFEM_fprintf(out,"NDOFE:   %d\n",info->ndofe);
  PGFEM_fprintf(out,"ELEM_ID: %d\n",info->elem_id);
  PGFEM_fprintf(out,"PROC_ID: %d\n",info->proc_id);
  PGFEM_fprintf(out,"INT_PT:  %d\n",info->int_pt);
  PGFEM_fprintf(out,"CELL_ID: %d\n",cell_id);

  job_type_str(info->job_type,&job_str);
  PGFEM_fprintf(out,"JOB TYPE: %s\n",job_str);
  free(job_str); job_str = NULL;

  PGFEM_fprintf(out,"PRINT:   %d\n",info->print_flag);
  PGFEM_fprintf(out,"TIM:     %d\n",info->tim);
  PGFEM_fprintf(out,"TIMES:   %3.5e %3.5e %3.5e\n",
		info->times[0],info->times[1],info->times[2]);
  PGFEM_fprintf(out,"N_STEP:  %d\n",info->n_step);
  PGFEM_fprintf(out,"JUMP(n): %3.5e %3.5e %3.5e\n",
		info->jump_n[0],info->jump_n[1],info->jump_n[2]);
  PGFEM_fprintf(out,"JUMP:    %3.5e %3.5e %3.5e\n",
		info->jump[0],info->jump[1],info->jump[2]);
  PGFEM_fprintf(out,"NORMAL:  %3.5e %3.5e %3.5e\n",
		info->normal[0],info->normal[1],info->normal[2]);
  PGFEM_fprintf(out,"INT_WT:  %3.5e\n",info->int_wt);

  PGFEM_fprintf(out,"SHAPE:   ");
  for(int i=0; i<info->nnode; i++){
    PGFEM_fprintf(out,"%3.5e ",info->shape[i]);
  }
  PGFEM_fprintf(out,"\n");

  PGFEM_fprintf(out,"TRAC(n): %3.5e %3.5e %3.5e\n",
		info->traction_n[0],info->traction_n[1],info->traction_n[2]);
  PGFEM_fprintf(out,"TRAC:    %3.5e %3.5e %3.5e\n",
		info->traction[0],info->traction[1],info->traction[2]);
  PGFEM_fprintf(out,"LID: ");
  print_array_l(out,info->loc_dof_ids,info->ndofe,1,info->ndofe);
  PGFEM_fprintf(out,"GID: ");
  print_array_l(out,info->g_dof_ids,info->ndofe,1,info->ndofe);

  /* print tangent and residual if debugging */
  if(DEBUG_MS_JOB_INFO){
    PGFEM_fprintf(out,"RES: ");
    print_array_d(out,info->traction_res,info->ndofe,1,info->ndofe);
    PGFEM_fprintf(out,"TANGENT:\n");
    print_array_d(out,info->K_00_contrib,info->ndofe*info->ndofe,
		  info->ndofe,info->ndofe);
  }

  PGFEM_fprintf(out,"====== END JOB INFO ======\n");

  return err;
}

