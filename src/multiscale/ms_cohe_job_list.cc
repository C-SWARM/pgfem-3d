#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "PLoc_Sparse.h"
#include "allocation.h"
#include "cohesive_element_utils.h"
#include "compute_ms_cohe_job.h"
#include "get_dof_ids_on_elem.h"
#include "incl.h"
#include "ms_cohe_job_list.h"
#include "quadrature_rules.h"
#include "utils.h"

using namespace pgfem3d;
using namespace pgfem3d::net;
using pgfem3d::solvers::SparseSystem;

static constexpr int ndim = 3;

/*=== LOC utils ===*/
/** print diagnostic message and return code for further action */
static int
CHECK_WARNING(int err_code, int rank, const char *func, const char *file, int line)
{
  if(err_code){
    PGFEM_printerr("[%d]WARINING: received error code %d! %s:%s:%d\n",
                   rank,err_code,func,file,line);
  }
  return err_code;
}

#define check_warning(err_code,rank) CHECK_WARNING(err_code,rank,__func__,__FILE__,__LINE__)

/*==== STATIC FUNCTION PROTOTYPES ====*/

/** compute the number of multiscale coheisve jobs on the domain. */
static int compute_local_ms_cohe_n_jobs(const long nce,
                                        const COEL *coel,
                                        const Node *node,
                                        long *n_jobs);

/** Loop through the cohesive elements and pre-compute all of the
    job information that will remain constant throughout the
    analysis. Creates a job for each integration point. */
static int create_local_ms_cohe_job_list(const long nce,
                                         const COEL *coel,
                                         const Node *node,
                                         const int group_id,
                                         const long n_jobs,
                                         MS_COHE_JOB_INFO *job_list,
                                         int *local_buffer_size,
					 int myrank_macro,
                                         const int mp_id);

/** update the local job list displacement jumps */
static int update_loc_ms_cohe_job_list(const int nce,
                                       const COEL *coel,
                                       const Node *node,
                                       const SUPP sup,
                                       const double *sol,
                                       MS_COHE_JOB_INFO *loc_job_list,
                                       int *local_buffer_size);

/** allocates n_job_dom for the number of processors in ms_comm and
    polulates with the number of jobs on each proc. Sums the number of
    jobs into Gn_job and computes start idx. Contains collective
    communication on ms_comm */
static int compute_loc_job_list_metadata(const long nce,
                                         const COEL *coel,
                                         const Node *node,
                                         long **n_job_dom,
                                         long *Gn_job,
                                         long *start_id,
					 Network *net,
                                         PGFem3D_Comm ms_comm);

/** Distribute the job list to all procs on the communicator. Contains
    collective communication on ms_comm */
static int  distribute_group_ms_cohe_job_list(MS_COHE_JOB_INFO *job_list,
                                              const int start_id,
                                              /* const */ int buff_size,
                                              const long *n_job_dom,
					      Network *net,
                                              PGFem3D_Comm ms_comm);

/*==== API FUNCTION DEFINITIONS ====*/

int create_group_ms_cohe_job_list(const long nce,
                                  const COEL *coel,
                                  const Node *node,
                                  const PGFem3D_Comm macro_comm,
                                  const PGFem3D_Comm ms_comm,
                                  const int group_id,
                                  long *Gn_jobs,
                                  long **n_job_dom,
                                  MS_COHE_JOB_INFO **job_list,
				  Network *net,
                                  const int mp_id)
{
  int err = 0;
  int myrank = 0;
  int macro_rank = 0;
  int nproc = 0;
  net->comm_rank(ms_comm,&myrank);
  net->comm_size(ms_comm,&nproc);
  net->comm_rank(macro_comm,&macro_rank);
  
  int *buff_sizes = NULL;
  int *buff_starts = NULL;
  char *buffer = NULL;
  long job_id_start = 0;

  *job_list = NULL;
  err += compute_loc_job_list_metadata(nce,coel,node,n_job_dom,
                                       Gn_jobs,&job_id_start,net,ms_comm);

  /* check error status */
  if(check_warning(err,myrank)) goto exit_function;

  /* exit early if no jobs on communication group */
  if(*Gn_jobs <= 0) goto exit_function;

  /* allocate list */
  *job_list = PGFEM_calloc(MS_COHE_JOB_INFO, *Gn_jobs);

  buff_sizes = PGFEM_calloc(int, nproc);
  buff_sizes[myrank] = 0;
  err += create_local_ms_cohe_job_list(nce,coel,node,group_id,
                                       (*n_job_dom)[myrank],
                                       *job_list + job_id_start,
                                       &buff_sizes[myrank],macro_rank,
				       mp_id);

  /* check error status */
  if(check_warning(err,myrank)) goto exit_function;

  /* compute size of global buffer */
  buff_starts = PGFEM_calloc(int, nproc);
  net->allgather(NET_IN_PLACE,1,NET_DT_INT,buff_sizes,1,NET_DT_INT,ms_comm);
  {
    size_t g_buff_size = 0;
    for(int i=0; i<nproc; i++){
      buff_starts[i] = g_buff_size;
      g_buff_size += buff_sizes[i];
    }
    buffer = PGFEM_calloc(char, g_buff_size);
  }

  /* pack the local job info */
  {
    size_t pos = buff_starts[myrank];
    for(long i=0; i<(*n_job_dom)[myrank]; i++){
      MS_COHE_JOB_INFO *info = *job_list + job_id_start + i;
      size_t len = compute_MS_COHE_JOB_INFO_size(info);
      err += pack_MS_COHE_JOB_INFO(info,len,buffer+pos);
      pos += len;
    }
  }

  /* gather on all processes in group */
  net->allgatherv(NET_IN_PLACE,buff_sizes[myrank],NET_DT_CHAR,
		  buffer,buff_sizes,buff_starts,NET_DT_CHAR,ms_comm);
  
  /* check error status */
  if(check_warning(err,myrank)) goto exit_function;

  /* unpack data from other processes */
  for(int i=0; i<nproc; i++){
    /* compute start */
    if(i == 0) job_id_start = 0;
    else job_id_start += (*n_job_dom)[i-1];

    if(i == myrank)continue;

    size_t pos = buff_starts[i];
    for(long j=0; j<(*n_job_dom)[i]; j++){
      MS_COHE_JOB_INFO *info = *job_list + job_id_start + j;
      /* peek the number of nodes from the buffer at pos */
      {
        int nnode = 0;
        size_t tmp = 0;
        unpack_data(buffer+pos,&nnode,&tmp,1,sizeof(int));
        err += build_MS_COHE_JOB_INFO(info,nnode);
      }
      size_t len_job = compute_MS_COHE_JOB_INFO_size(info);
      err += unpack_MS_COHE_JOB_INFO(info,len_job,buffer+pos);
      pos += len_job;
      buff_sizes[i] -= len_job;
    }
    if(buff_sizes[i] != 0) err++;

    /* check error status */
    if(check_warning(err,myrank)) goto exit_function;
  }

 exit_function:
  free(buffer);
  free(buff_sizes);
  free(buff_starts);
  return err;
} /* create_group_ms_cohe_job_list() */

int update_group_ms_cohe_job_list(const long nce,
                                  const COEL *coel,
                                  const Node *node,
                                  const SUPP sup,
                                  const double *sol,
				  Network *net,
                                  PGFem3D_Comm ms_comm,
                                  MS_COHE_JOB_INFO *job_list)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  net->comm_rank(ms_comm,&myrank);
  net->comm_size(ms_comm,&nproc);


  long *n_job_dom = NULL;
  long Gn_jobs =  0;
  long start_id = 0;
  int buff_size = 0;
  err += compute_loc_job_list_metadata(nce,coel,node,&n_job_dom,
                                       &Gn_jobs,&start_id,net,ms_comm);

  /* check error status/zero jobs */
  if(check_warning(err,myrank) || Gn_jobs <= 0) goto exit_function;

  /* update local jobs and compute buffer size */
  err += update_loc_ms_cohe_job_list(nce,coel,node,sup,sol,
                                     job_list + start_id,
                                     &buff_size);

  /* check error status */
  if(check_warning(err,myrank)) goto exit_function;

  /* distribute the job list */
  err += distribute_group_ms_cohe_job_list(job_list,
                                           start_id,buff_size,
                                           n_job_dom,net,ms_comm);
  /* check error status */
  if(check_warning(err,myrank));

 exit_function:
  free(n_job_dom);

  return err;
}/* update_group_ms_cohe_job_list() */

int compute_ms_cohe_tan_res(const int compute_micro_eq, //deprecated
                            const CommunicationStructure *com,
                            MS_COHE_JOB_INFO *job_list,
                            SparseSystem *macro_solver,
                            Microscale *microscale,
                            const int mp_id)
{
  int err = 0;

  MultiscaleCommon         *c = microscale;
  SparseSystem  *micro_solver = c->SOLVER;
  const int            n_sols = microscale->idx_map.size;
  const int          analysis = microscale->opts->analysis_type;

  /* get MPI ranks */
  int  macro_rank = com->rank;
  int macro_nproc = com->nproc;

  /* redirect the I/O to micro logging */
  PGFEM_redirect_io_micro();

  /* destroy the preconditioner object(s) and re-initialize the
     microscale preconditioner */
  micro_solver->resetPreconditioner();

  /* setup the stiffness matrix communication */
  double        **Lk = nullptr;
  double   **receive = nullptr;
  com->spc->post_stiffmat(&Lk, &receive);

  /* for each solution, compute job */
  for (int i = 0; i < n_sols; ++i) {
    MS_COHE_JOB_INFO *job = job_list + i;
    err += compute_ms_cohe_job(i, job, microscale, mp_id);

    /* assemble tangent if owning process */
    if (macro_rank == job->proc_id) {
      PLoc_Sparse(Lk, job->K_00_contrib, nullptr, nullptr, nullptr,
                  job->g_dof_ids, job->ndofe, nullptr, 0,
                  macro_rank, macro_nproc, com->spc,
                  0/* interior (t/f) */, macro_solver, analysis);
    }
  }

  /*=== OPTIMIZATION NOTES: overlay computation and communication as
    in the volumetric elements in stiffmat_fd */

  /* send/finalize the communication */
  com->spc->send_stiffmat();
  com->spc->finalize_stiffmat();

  /* assemble to macro tangent on this process */
  {
    SparseSystem::sp_idx *row_idx = nullptr;
    SparseSystem::sp_idx   *ncols = nullptr;
    SparseSystem::sp_idx *col_idx = nullptr;
    for (int i = 0; i < com->spc->Nr; ++i) {             //ksaha
      const int proc = com->spc->Nrr[i];
      const int nrows = com->spc->R[proc];

      /* allocate rows and cols to receive */
      row_idx = PGFEM_calloc(SparseSystem::sp_idx, nrows);
      ncols   = PGFEM_calloc(SparseSystem::sp_idx, nrows);
      col_idx = PGFEM_calloc(SparseSystem::sp_idx, com->spc->AR[proc]);

      /* get row and column ids */
      SparseSystem::sp_idx idx = 0;
      for (int j = 0, e = com->spc->R[proc]; j < e; ++j) {
        row_idx[j] = com->spc->RGID[proc][j];
        ncols[j] = com->spc->RAp[proc][j];
        for (SparseSystem::sp_idx k = 0, e = ncols[j]; k < e; ++k) {
          col_idx[idx] = com->spc->RGRId[proc][idx];
          ++idx;
        }
      }

      /* assemble to local part of global stiffness */
      // err += HYPRE_IJMatrixAddToValues(macro_solver->hypre_k, nrows, ncols,
      //                                  row_idx, col_idx, receive[proc]);
      macro_solver->add(nrows, ncols, row_idx, col_idx, receive[proc]);
      free(row_idx);
      free(ncols);
      free(col_idx);
    }
  } /* finish assembly */

  /* reset I/O to macro logging */
  PGFEM_redirect_io_macro();

  /* re-initialize macroscale preconditioner */
  macro_solver->resetPreconditioner();

  /* exit function */
  return err;
}/* compute_ms_cohe_tan_res() */

int assemble_ms_cohe_res(const Microscale *micro,
                         const MS_COHE_JOB_INFO *jobs,
                         const PGFem3D_Comm macro_mpi_comm,
                         double *macro_loc_res,
                         const int mp_id)
{
  int err = 0;
  int macro_rank;
  micro->net->comm_rank(macro_mpi_comm, &macro_rank);
  
  /* redirect I/O to microscale */
  PGFEM_redirect_io_micro();

  /* assemble the residual from each job on the local part of the
     macroscale residual. Assembly takes place on the macro owning
     process only!! */
  for(int i=0, e = micro->idx_map.size; i < e; i++){
    err += assemble_ms_cohe_job_res(i,jobs+i,
                                    micro->rank,
				    macro_rank,
                                    macro_loc_res,
                                    mp_id);
  }

  /* redirect I/O to macroscale */
  PGFEM_redirect_io_macro();
  return err;
}/* assemble_ms_cohe_res() */


void destroy_ms_cohe_job_list(const long Gn_job,
                              MS_COHE_JOB_INFO *job_list)
{
  if(job_list != NULL){
    for(long i=0; i<Gn_job; i++){
      destroy_MS_COHE_JOB_INFO(job_list + i);
    }
  }
  free(job_list);
}/* destroy_ms_cohe_job_list() */

/*==== STATIC FUNCTION DEFINITIONS ====*/

static int compute_local_ms_cohe_n_jobs(const long nce,
                                        const COEL *coel,
                                        const Node *node,
                                        long *n_jobs)
{
  int err = 0;
  /* compute number of jobs I need to allocate */
  *n_jobs = nce;
  /* for(int i=0; i<nce; i++){ */
  /*   *n_jobs += int_pointC(coel[i].toe/2); */
  /* } */

  if(*n_jobs != nce){
    PGFEM_printerr("ERROR: Only support one integration point"
                   " per macro element!\n");
    PGFEM_Abort();
  }
  return err;
}

static int create_local_ms_cohe_job_list(const long nce,
                                         const COEL *coel,
                                         const Node *node,
                                         const int group_id,
                                         const long n_jobs,
                                         MS_COHE_JOB_INFO *job_list,
                                         int *local_buffer_size,
					 int myrank_macro,
                                         const int mp_id)
{
  int err = 0;
  /* exit early if there are no jobs on this domain */
  if(n_jobs <= 0)return err;

  *local_buffer_size = 0;
  int job_id = 0;
  double *normal = PGFEM_calloc(double, ndim);
  double *jump = PGFEM_calloc(double, ndim);

  for(int i=0; i<nce; i++){
    /* information that is constant per element */
    const COEL *cel = &coel[i];
    const int nne = cel->toe;
    const int nne_2D = nne/2;
    double *shape_2D = PGFEM_calloc(double, nne_2D);
    double *N_x = PGFEM_calloc(double, nne_2D);
    double *N_y = PGFEM_calloc(double, nne_2D);
    double *shape = PGFEM_calloc(double, nne);
    long *loc_dof_ids = PGFEM_calloc(long, nne*ndim);
    long *g_dof_ids = PGFEM_calloc(long, nne*ndim);
    get_dof_ids_on_elem_nodes(0,nne,ndim,cel->nod,node,loc_dof_ids,mp_id);
    get_dof_ids_on_elem_nodes(1,nne,ndim,cel->nod,node,g_dof_ids  ,mp_id);

    double *x = PGFEM_calloc(double, nne);
    double *y = PGFEM_calloc(double, nne);
    double *z = PGFEM_calloc(double, nne);
    double *disp = PGFEM_calloc(double, nne*ndim);
    nodecoord_updated(nne,cel->nod,node,x,y,z);

    double *xl = PGFEM_calloc(double, nne_2D);
    double *yl = PGFEM_calloc(double, nne_2D);
    double *zl = PGFEM_calloc(double, nne_2D);
    tran_coord(nne_2D,cel->x,cel->y,cel->z,cel->e1,
               cel->e2,cel->n,xl,yl,zl,1);

    double *xb = PGFEM_calloc(double, nne_2D);
    double *yb = PGFEM_calloc(double, nne_2D);
    double *zb = PGFEM_calloc(double, nne_2D);
    mean_map(nne_2D,x,y,z,disp,xb,yb,zb);

    /* set up integration SINGLE INTEGRATION POINT */
    int n_ip = 0;
    double *gk = NULL;
    double *ge = NULL;
    double *w = NULL;
    switch(nne_2D){
     case 3: err += get_tria_quadrature_rule(0,&n_ip,&gk,&ge,&w); break;
     case 4: err += get_quad_quadrature_rule(0,&n_ip,&gk,&ge,&w); break;
    }

    double *e1 = PGFEM_calloc(double, ndim);
    double *e2 = PGFEM_calloc(double, ndim);
    double *e2h = PGFEM_calloc(double, ndim);

    /* for each integration point */
    for(int ip = 0; ip<n_ip; ip++){
      const double ksi = gk[ip];
      const double eta = ge[ip];
      double wt = w[ip];

      /* build the job for the integration point */
      MS_COHE_JOB_INFO *job = &job_list[job_id];
      err += build_MS_COHE_JOB_INFO(job,nne);

      /* reset information */
      nulld(jump,ndim);
      nulld(e1,ndim);
      nulld(e2,ndim);
      nulld(e2h,ndim);
      nulld(normal,ndim);
      nulld(shape,nne);
      nulld(shape_2D,nne_2D);
      nulld(N_x,nne_2D);
      nulld(N_y,nne_2D);

      /* compute normal */
      base_vec(nne_2D,ksi,eta,xb,yb,zb,e1,e2,e2h,normal,0);

      /* compute integral transformation and multiply with wt */
      wt *= dN3_xy(ksi,eta,nne_2D,xl,yl,zl,N_x,N_y);

      /* compute 2D shape functions */
      shape_2DC(nne_2D,ksi,eta,shape_2D);

      /* push shape functions into array for 3D interface element */
      for(int k=0; k<nne_2D; k++){
        shape[k] = -shape_2D[k]; /* - */
        shape[k+nne_2D] = shape_2D[k]; /* + */
      }

      /* compute the jump accross the interface */
      get_jump(nne_2D,x,y,z,disp,shape_2D,jump);

      /* set the job information */
      err += set_MS_COHE_JOB_INFO(job,normal,jump,shape,
                                  loc_dof_ids,g_dof_ids);
      job->int_wt = wt;
      job->elem_id = i;
      job->proc_id = myrank_macro;
      job->int_pt = ip;
      job->job_type = JOB_NO_COMPUTE_EQUILIBRIUM;

      /* set print flag equal to property. Non-zero value will cause
         microscale output. */
      job->print_flag = cel->pr;

      (*local_buffer_size) += compute_MS_COHE_JOB_INFO_size(job);

      job_id++;
    }/* each integration point/job */

    free(shape_2D);
    free(N_x);
    free(N_y);
    free(shape);
    free(g_dof_ids);
    free(loc_dof_ids);
    free(x);
    free(y);
    free(z);
    free(disp);
    free(xl);
    free(yl);
    free(zl);
    free(xb);
    free(yb);
    free(zb);
    free(e1);
    free(e2);
    free(e2h);
    free(gk);
    free(ge);
    free(w);
  }/* each cohesive element */

  free(normal);
  free(jump);
  return err;
}

static int update_loc_ms_cohe_job_list(const int nce,
                                       const COEL *coel,
                                       const Node *node,
                                       const SUPP sup,
                                       const double *sol,
                                       MS_COHE_JOB_INFO *loc_job_list,
                                       int *local_buffer_size)
{
  int err = 0;
  /* NOTE: jobs are sequential by element. Thus, the coodinates and
     displacement on subsequent jobs associated with the same element
     are identical. */

  int job_id = 0;
  for(int i=0; i<nce; i++){
    const COEL *cel = coel + i;
    const int nne = cel->toe;
    const int nne_2D = nne/2;
    double *x = PGFEM_calloc(double, nne);
    double *y = PGFEM_calloc(double, nne);
    double *z = PGFEM_calloc(double, nne);
    double *disp = PGFEM_calloc(double, nne*ndim);

    /* get nodal coordinates and displacements on the element */
    nodecoord_updated(nne,cel->nod,node,x,y,z);
    def_elem(loc_job_list[job_id].loc_dof_ids,nne*ndim,
             sol,NULL,node,disp,sup,0);

    /* set up integration SINGLE INTEGRATION POINT */
    double *gk = NULL;
    double *ge = NULL;
    double *w = NULL;
    int n_ip = 0;
    switch(nne_2D){
     case 3: err += get_tria_quadrature_rule(0,&n_ip,&gk,&ge,&w); break;
     case 4: err += get_quad_quadrature_rule(0,&n_ip,&gk,&ge,&w); break;
    }

    /* use shape functions stored in the job list to compute the new
       jump. ONLY THE JUMP IS UPDATED */
    for(int ip=0; ip<n_ip; ip++){
      MS_COHE_JOB_INFO *job = &loc_job_list[job_id];
      get_jump(nne_2D,x,y,z,disp,job->shape + nne_2D,job->jump);
      (*local_buffer_size) += compute_MS_COHE_JOB_INFO_size(job);
      job_id++;
    }

    free(x);
    free(y);
    free(z);
    free(disp);
    free(gk);
    free(ge);
    free(w);
  }
  return err;
}

static int compute_loc_job_list_metadata(const long nce,
                                         const COEL *coel,
                                         const Node *node,
                                         long **n_job_dom,
                                         long *Gn_job,
                                         long *start_id,
					 Network *net,
                                         PGFem3D_Comm ms_comm)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  net->comm_rank(ms_comm,&myrank);
  net->comm_size(ms_comm,&nproc);

  *n_job_dom = PGFEM_calloc(long, nproc);
  err += compute_local_ms_cohe_n_jobs(nce,coel,node,*n_job_dom + myrank);
  net->allgather(NET_IN_PLACE,1,NET_DT_LONG,
		 *n_job_dom,1,NET_DT_LONG,ms_comm);

  *Gn_job = 0;
  for(int i=0; i<nproc; i++){
    *Gn_job += (*n_job_dom)[i];
  }

  *start_id = 0;
  for(int i=0; i<myrank; i++){
    start_id += (*n_job_dom)[i];
  }

  return err;
}


static int  distribute_group_ms_cohe_job_list(MS_COHE_JOB_INFO *job_list,
                                              const int start_id,
                                              /* const */ int buff_size,
                                              const long *n_job_dom,
					      Network *net,
                                              PGFem3D_Comm ms_comm)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  net->comm_rank(ms_comm,&myrank);
  net->comm_size(ms_comm,&nproc);

  /* exit early if no communication is required */
  if(nproc <= 1) return err;

  /* allocate space for the buffer sizes and starts. */
  int *buff_sizes = PGFEM_calloc(int, nproc);
  int *buff_starts = PGFEM_calloc(int, nproc);

  net->allgather(&buff_size,1,NET_DT_INT,
		 buff_sizes,1,NET_DT_INT,ms_comm);
  
  /* compute total buffer size and starts */
  size_t g_buff_size = 0;
  for(int i=0; i<nproc; i++){
    buff_starts[i] = g_buff_size;
    g_buff_size += buff_sizes[i];
  }
  char * buffer = PGFEM_calloc(char, g_buff_size);

  /* pack the local job info */
  {
    size_t pos = buff_starts[myrank];
    for(long i=0; i<n_job_dom[myrank]; i++){
      const MS_COHE_JOB_INFO *info = job_list + start_id + i;
      const size_t len = compute_MS_COHE_JOB_INFO_size(info);
      err += pack_MS_COHE_JOB_INFO(info,len,buffer + pos);
      pos += len;
    }
  }

  /* gather buffer on all procs */
  net->allgatherv(NET_IN_PLACE,buff_sizes[myrank],NET_DT_CHAR,
		  buffer,buff_sizes,buff_starts,NET_DT_CHAR,ms_comm);
  
  /* unpack data from other processes */
  int loc_start = 0;
  for(int i=0; i<nproc; i++){
    if(i == 0) loc_start = 0;
    else loc_start += n_job_dom[i-1];
    if(i == myrank) continue;

    size_t pos = buff_starts[i];
    for(long j=0; j<n_job_dom[i]; j++){
      MS_COHE_JOB_INFO *info = job_list + loc_start + j;
      size_t len_job = compute_MS_COHE_JOB_INFO_size(info);
      err += unpack_MS_COHE_JOB_INFO(info,len_job,buffer+pos);
      pos += len_job;
      buff_sizes[i] -= len_job;
    }
    if(buff_sizes[i] != 0) err++;
  }

  free(buff_sizes);
  free(buff_starts);
  free(buffer);
  return err;
}
