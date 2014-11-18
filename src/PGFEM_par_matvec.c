/* HEADER */
#include "PGFEM_par_matvec.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mkl_cblas.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

#ifndef INDEX_MACROS_H
#include "index_macros.h"
#endif

/*=== ENTRY ===*/
typedef struct entry{
  int row;
  int col;
  double val;
} entry;
static inline void build_entry(entry *e);
static inline void set_entry(const int row,
			     const int col,
			     const double val,
			     entry *e);
static inline void destroy_entry(entry *e);
static inline void destroy_entries(const int n_entry,
				   entry *e);
static int compare_entry_row(const void *a, const void *b);
static int get_sorted_list_of_entries(const int n_entries,
				      const int *row_idx,
				      const int *col_idx,
				      const double *vals,
				      entry **e);

/*=== COMM_INFO ===*/
typedef struct comm_info{
  int nproc; /* nproc to communicate with */
  int *proc; /* proc list */
  int *n_info; /* number of rows to s/r with each proc */
} comm_info;
static int build_comm_info(const int nproc,
			   const int myrank,
			   const int *proc_comm_info,
			   comm_info *info);
static void destroy_comm_info(comm_info *info);
static void print_comm_info_graph(const int myrank,
				  comm_info *send,
				  comm_info *recv);

/*=== OFF_PROC_ENTRIES ===*/
typedef struct off_proc_entries{
  int *row_loc_ids; /* loc row id on owning domain */
  double *data; /* com_info.n_info[proc]*n_col */
}off_proc_entries;
static int build_off_proc_entries(const int n_rows,
				  const int n_col,
				  off_proc_entries *info);
static void destroy_off_proc_entries(off_proc_entries *info);

/*=== COMM_CTX ===*/
typedef struct comm_ctx{
  MPI_Request *r_req;
  MPI_Request *s_req;
  MPI_Status *r_stat;
  MPI_Status *s_stat;
} comm_ctx;
static int build_comm_ctx(const int n_send,
			  const int n_recv,
			  comm_ctx *COMM);
static void destroy_comm_ctx(comm_ctx *COMM);

/*=== LOC UTILS ===*/
#ifndef NDEBUG
#define my_assert(arg) \
  do {							    \
    if(!(arg)){						    \
      PGFEM_printerr("Assertion ("#arg") failed! %s:%s:%d\n",	\
		     __func__,__FILE__,__LINE__);		\
      PGFEM_Abort(); } } while(0)
#else 
#define my_assert(arg) do{/* nothing */}
#endif

#ifndef LOGGING
#define LOGGING 0
#endif

#if LOGGING == 1
void LOG_MSG(const char *msg)
{
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  PGFEM_printerr("%s\t%s\n",msg,asctime(timeinfo));
  /* free(timeinfo); */
}
#else
void LOG_MSG(const char *msg){}
#endif

static int PGFEM_par_matvec_compare_int(const void *a, const void *b);
static int set_local_value(const int row_start,
			   const int n_row,
			   const int n_col,
			   const entry *e,
			   double *data);
static int set_off_proc_value(const int row_start,
			      const int n_row,
			      const int n_col,
			      const entry *e,
			      const int *row_ids,
			      double *data);
static int add_to_local_value(const int row_start,
			      const int n_row,
			      const int n_col,
			      const entry *e,
			      double *data);
static int add_to_off_proc_value(const int row_start,
				 const int n_row,
				 const int n_col,
				 const entry *e,
				 const int *row_ids,
				 double *data);
static int add_assemble_matrix(PGFEM_par_matrix *mat);
static int set_assemble_matrix(PGFEM_par_matrix *mat);

/*==== API FUNCTIONS ====*/
int initialize_PGFEM_par_matrix(const int n_rows,
				const int n_cols,
				const int n_own_rows,
				const int n_entries,
				const int *row_idx,
				const int *col_idx,
				MPI_Comm mpi_comm,
				PGFEM_par_matrix **mat)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  MPI_Comm_rank(mpi_comm,&myrank);
  MPI_Comm_size(mpi_comm,&nproc);

  if(myrank == 0) LOG_MSG("Entering initialization.");

  /* allocate the object and get alias */
  *mat = PGFEM_calloc(1,sizeof(PGFEM_par_matrix));
  PGFEM_par_matrix *m = *mat;

  /* set internal variables */
  /* duplicate communicator to get unique context. Allows simultaneous
     assembly of multiple matrices without comm errors */
  MPI_Comm_dup(mpi_comm,&m->mpi_comm);
  m->assembled = 0;
  m->add_values = 0;
  m->n_rows = n_rows;
  m->n_cols = n_cols;
  m->n_own_rows = PGFEM_calloc(nproc,sizeof(int));
  m->idx_starts = PGFEM_calloc(nproc+1,sizeof(int));
  m->data = PGFEM_calloc(n_cols*n_own_rows,sizeof(double));
  m->send_info = PGFEM_calloc(1,sizeof(comm_info));
  m->recv_info = PGFEM_calloc(1,sizeof(comm_info));
  m->s_rows = NULL;
  m->r_rows = NULL;

  /* get aliases to send info */
  comm_info *send_info = (comm_info*) m->send_info;
  comm_info *recv_info = (comm_info*) m->recv_info;

  /* post receive of information */
  int *r_off_proc_rows = PGFEM_calloc(nproc,sizeof(int));
  MPI_Request *r_req = NULL;
  MPI_Status *r_stat = NULL;

  if(myrank == 0) LOG_MSG("\tPosting non-blocking receives.");
  if(nproc > 1){
    r_req = PGFEM_calloc(nproc-1,sizeof(MPI_Request));
    r_stat = PGFEM_calloc(nproc-1,sizeof(MPI_Status));
    {
      int idx = 0;
      int req_idx = 0;
      for(int i=0; i<nproc; i++){
	if(i == myrank) continue;
	MPI_Irecv(&r_off_proc_rows[i],1,MPI_INT,i,MPI_ANY_TAG,
		  m->mpi_comm,&r_req[req_idx]);
	req_idx++;
      }
    }
  }

  /* gather n_own_rows on all processes and compute the idx_starts */
  m->n_own_rows[myrank] = n_own_rows;
  MPI_Allgather(MPI_IN_PLACE,1,MPI_INT,m->n_own_rows,
		1,MPI_INT,m->mpi_comm);
  m->idx_starts[0] = 0;
  for(int i=1; i<=nproc; i++){
    m->idx_starts[i] = m->idx_starts[i-1] + m->n_own_rows[i-1];
  }

  /* push row/col info into entry, sort */
  entry *e = PGFEM_calloc(n_entries,sizeof(entry));
  for(int i=0; i<n_entries; i++){
    build_entry(e+i);
    set_entry(row_idx[i],col_idx[i],0.0,e+i);
  }

  qsort(e,n_entries,sizeof(entry),compare_entry_row);

  /* unique rows */
  int n_dup = 0;
  {
    int idx = 0;
    for(int i=1; i<n_entries; i++){
      if(e[i].row <= e[idx].row){
	e[i].row = -1;
	n_dup++;
      } else {
	idx = i;
      }
    }
  }

  /* sort entries, duplicates are first */
  qsort(e,n_entries,sizeof(entry),compare_entry_row);

  /* check that last duplicate is valid entry */
  if(e[n_dup].row < 0) n_dup++;

  /* count number of rows for each process */
  int *s_off_proc_rows = PGFEM_calloc(nproc,sizeof(int));
  {
    int idx = n_dup;
    int total_rows = 0;
    for(int i=0; i<nproc; i++){
      while(idx < n_entries){
	/* break loop if row owned by other dom */
	if(e[idx].row >= m->idx_starts[i+1]) break;
	s_off_proc_rows[i] ++;
	idx ++;
      }
    }
  }

  /* broadcast how much I am sending to everyone */
  if(myrank == 0) LOG_MSG("\tPost non-blocking sends.");
  MPI_Request *s_req = NULL;
  MPI_Status *s_stat = NULL;
  if(nproc > 1){
    s_req = PGFEM_calloc(nproc-1,sizeof(MPI_Request));
    s_stat = PGFEM_calloc(nproc-1,sizeof(MPI_Status));
    int req_idx = 0;
    for(int i=0; i<nproc; i++){
      if(i == myrank) continue;
      MPI_Isend(s_off_proc_rows + i,1,MPI_INT,i,i /*tag*/,
		m->mpi_comm,&s_req[req_idx]);
      req_idx++;
    }
  }

  /* build send communication information structure */
  if(myrank == 0) LOG_MSG("\tBuild send info.");
  err += build_comm_info(nproc,myrank,s_off_proc_rows,send_info);

  /* build s_rows and get alias */
  if(myrank == 0) LOG_MSG("\tBuild off-proc send data.");
  if(send_info->nproc > 0){
    m->s_rows = PGFEM_calloc(send_info->nproc,sizeof(off_proc_entries));
  }
  off_proc_entries *s_rows = (off_proc_entries*) m->s_rows;
  for(int i=0; i<send_info->nproc; i++){
    err += build_off_proc_entries(send_info->n_info[i],
				  m->n_cols,s_rows+i);
  }

  /* set the row_loc_ids */
  {
    int ent_idx = n_dup;
    for(int i=0; i<send_info->nproc; i++){
      /* scan entry rows until the first off-process start */
      int off_proc_start = m->idx_starts[send_info->proc[i]];
      while(e[ent_idx].row < off_proc_start) ent_idx++;
      for(int j=0; j<send_info->n_info[i]; j++){
	s_rows[i].row_loc_ids[j] = e[ent_idx].row - off_proc_start;
	ent_idx++;
      }
    }
  }

  /* wait for communication to finish */
  if(nproc > 1){
    MPI_Waitall(nproc-1,r_req,r_stat);
    MPI_Waitall(nproc-1,s_req,s_stat);
  }

  /* free and reset requests and statuses */
  free(r_req); r_req = NULL;
  free(s_req); s_req = NULL;
  free(r_stat); r_stat = NULL;
  free(s_stat); s_stat = NULL;

  /* build recv_info and r_rows */
  if(myrank == 0) LOG_MSG("\tBuild off-proc recv data.");
  err += build_comm_info(nproc,myrank,r_off_proc_rows,recv_info);

  if(LOGGING){
    print_comm_info_graph(myrank,send_info,recv_info);
  }

  if(recv_info->nproc > 0){
    m->r_rows = PGFEM_calloc(recv_info->nproc,sizeof(off_proc_entries));
  }

  off_proc_entries *r_rows = (off_proc_entries*) m->r_rows;
  for(int i=0; i<recv_info->nproc; i++){
    err += build_off_proc_entries(recv_info->n_info[i],
				  m->n_cols,r_rows+i);
  }

  /* communicate local row ids for off-proc entries */
  if(myrank == 0) LOG_MSG("\tCommunicate off-process row ids (nonblocking).");
  if(send_info->nproc > 0){
    s_req = PGFEM_calloc(send_info->nproc,sizeof(MPI_Request));
    s_stat = PGFEM_calloc(send_info->nproc,sizeof(MPI_Status));
  }else{
    PGFEM_printerr("[%d]WARNING: no send_info in %s!\n",myrank,__func__);
  }

  if(recv_info->nproc > 0){
    r_req = PGFEM_calloc(recv_info->nproc,sizeof(MPI_Request));
    r_stat = PGFEM_calloc(recv_info->nproc,sizeof(MPI_Status));
  }else{
    PGFEM_printerr("[%d]WARNING: no recv_info in %s!\n",myrank,__func__);
  }

  /* post receive */
  for(int i=0; i<recv_info->nproc; i++){
    MPI_Irecv(r_rows[i].row_loc_ids,recv_info->n_info[i],MPI_INT,
	      recv_info->proc[i],MPI_ANY_TAG,m->mpi_comm,&r_req[i]);
  }

  /* post send */
  for(int i=0; i<send_info->nproc; i++){
    MPI_Isend(s_rows[i].row_loc_ids,send_info->n_info[i],MPI_INT,
	      send_info->proc[i],i /*tag*/,m->mpi_comm,&s_req[i]);
  }

  /* complete communication */
  if(myrank == 0) LOG_MSG("\tWait all II.");
  MPI_Waitall(recv_info->nproc,r_req,r_stat);
  if(LOGGING){
    PGFEM_printerr("\t\t[%d] Finished RECV II\n",myrank);
    MPI_Barrier(m->mpi_comm);
    if(myrank == 0) LOG_MSG("\t\tFinished Receives.");
  }
  MPI_Waitall(send_info->nproc,s_req,s_stat);

  if(myrank == 0) LOG_MSG("\tBarrier II.");
  if(LOGGING){
    PGFEM_printerr("\t\t[%d] Finished Comm II\n",myrank);
    MPI_Barrier(m->mpi_comm);
  }

  /* deallocate local arrays/structures */
  free(r_req);
  free(s_req);
  free(r_stat);
  free(s_stat);
  destroy_entries(n_entries,e);
  free(s_off_proc_rows);
  free(r_off_proc_rows);

  /* exit function */
  if(myrank == 0) LOG_MSG("Exiting initialization");
  return err;
}

void destroy_PGFEM_par_matrix(PGFEM_par_matrix *mat)
{
  if(mat != NULL){
    free(mat->n_own_rows);
    free(mat->idx_starts);
    free(mat->data);

    /* destroy send information */
    comm_info *s = (comm_info *) mat->send_info;
    off_proc_entries *sr = (off_proc_entries *) mat->s_rows;
    for(int i=0; i<s->nproc; i++){
      destroy_off_proc_entries(sr+i);
    }
    free(sr);
    destroy_comm_info(s);

    /* destroy receive information */
    comm_info *r = (comm_info *) mat->recv_info;
    off_proc_entries *rr = (off_proc_entries *) mat->r_rows;
    for(int i=0; i<r->nproc; i++){
      destroy_off_proc_entries(rr+i);
    }
    free(rr);
    destroy_comm_info(r);

    /* destroy duplicated communicator */
    MPI_Comm_free(&(mat->mpi_comm));
  }
  free(mat);
}

int PGFEM_par_matrix_set_values(const int n_entries,
				const int *row_idx,
				const int *col_idx,
				const double *values,
				PGFEM_par_matrix *mat)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  MPI_Comm_rank(mat->mpi_comm,&myrank);
  MPI_Comm_size(mat->mpi_comm,&nproc);
  mat->assembled = 0;
  mat->add_values = 0;
  /* exit early if no entries */
  if(n_entries <= 0) return 0;

  entry *e = NULL;
  err += get_sorted_list_of_entries(n_entries,row_idx,col_idx,values,&e);

  /* get aliases to send information */
  comm_info *send = (comm_info *) mat->send_info;
  off_proc_entries *s_rows = (off_proc_entries *) mat->s_rows;

  /* get proc ids and sort */
  int *proc_order = PGFEM_calloc(send->nproc+1,sizeof(int));
  for(int i=0; i<send->nproc; i++){
    proc_order[i] = send->proc[i];
  }
  proc_order[send->nproc] = myrank;

  /* sort proc ids. Note that the communication order in 'send' is
     already ordered by proc id. Therefore indexing through
     'proc_order' will maintain the order in 'send' */
  qsort(proc_order,send->nproc+1,sizeof(int),
	PGFEM_par_matvec_compare_int);

  /* check to make sure that bounds of row values are valid */
  if(e[0].row < mat->idx_starts[proc_order[0]] ||
     e[n_entries-1].row >= mat->idx_starts[proc_order[send->nproc]+1]){
    PGFEM_printerr("[%d]ERROR: invalid row! %s:%s:%d\n",
	    myrank,__func__,__FILE__,__LINE__);
    PGFEM_Comm_code_abort(mat->mpi_comm,0);
  }

  /* loop through entries and set values */
  int ent_idx = 0;
  int send_idx = 0; /* see previous note */
  for(int i=0; i<=send->nproc; i++){
    int proc = proc_order[i];
    int row_start = mat->idx_starts[proc];
    int row_end = mat->idx_starts[proc+1];
    int n_rows = mat->n_own_rows[proc];
    if(proc == myrank){
      while(ent_idx < n_entries){
	/* break loop if owned by other domain */
	if(e[ent_idx].row >= row_end) break;
	err += set_local_value(row_start,n_rows,mat->n_cols,
			       e+ent_idx,mat->data);
	ent_idx++;
      }
    } else {
      while(ent_idx < n_entries){
	/* break loop if owned by other domain */
	if(e[ent_idx].row >= row_end) break;
	err += set_off_proc_value(row_start,send->n_info[send_idx],
				  mat->n_cols,e+ent_idx,
				  s_rows[send_idx].row_loc_ids,
				  s_rows[send_idx].data);
	ent_idx++;
      }
      send_idx++;
    }
  }

  if(ent_idx != n_entries){
    PGFEM_printerr("[%d]ERROR: did not reach end of entry list!\n"
	    "%s:%s:%d\n",myrank,__func__,__FILE__,__LINE__);
    PGFEM_Comm_code_abort(mat->mpi_comm,0);
  }

  free(proc_order);
  destroy_entries(n_entries,e);
  return err;
}

int PGFEM_par_matrix_add_to_values(const int n_entries,
				   const int *row_idx,
				   const int *col_idx,
				   const double *values,
				   PGFEM_par_matrix *mat)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  MPI_Comm_rank(mat->mpi_comm,&myrank);
  MPI_Comm_size(mat->mpi_comm,&nproc);
  mat->assembled = 0;
  mat->add_values = 1;
  /* exit early if no entries */
  if(n_entries <= 0) return 0;

  entry *e = NULL;
  err += get_sorted_list_of_entries(n_entries,row_idx,
				    col_idx,values,&e);

  /* get aliases to send information */
  comm_info *send = (comm_info *) mat->send_info;
  off_proc_entries *s_rows = (off_proc_entries *) mat->s_rows;

  /* get proc ids and sort */
  int *proc_order = PGFEM_calloc(send->nproc+1,sizeof(int));
  for(int i=0; i<send->nproc; i++){
    proc_order[i] = send->proc[i];
  }
  proc_order[send->nproc] = myrank;

  /* sort proc ids. Note that the communication order in 'send' is
     already ordered by proc id. Therefore indexing through
     'proc_order' will maintain the order in 'send' */
  qsort(proc_order,send->nproc+1,sizeof(int),
	PGFEM_par_matvec_compare_int);

  /* check to make sure that bounds of row values are valid */
  if(e[0].row < mat->idx_starts[proc_order[0]] ||
     e[n_entries-1].row >= mat->idx_starts[proc_order[send->nproc]+1]){
    PGFEM_printerr("[%d]ERROR: invalid row! %s:%s:%d\n",
	    myrank,__func__,__FILE__,__LINE__);
    PGFEM_Comm_code_abort(mat->mpi_comm,0);
  }

  /* loop through entries and set values */
  int ent_idx = 0;
  int send_idx = 0; /* see previous note */
  for(int i=0; i<=send->nproc; i++){
    int proc = proc_order[i];
    int row_start = mat->idx_starts[proc];
    int row_end = mat->idx_starts[proc+1];
    int n_rows = mat->n_own_rows[proc];
    if(proc == myrank){
       while(ent_idx < n_entries){
	/* break loop if owned by other domain */
	if(e[ent_idx].row >= row_end) break;
	err += add_to_local_value(row_start,n_rows,mat->n_cols,
				  e+ent_idx,mat->data);
	ent_idx++;
      }
    } else {
       while(ent_idx < n_entries){
	/* break loop if owned by other domain */
	if(e[ent_idx].row >= row_end) break;
	err += add_to_off_proc_value(row_start,send->n_info[send_idx],
				     mat->n_cols,e+ent_idx,
				     s_rows[send_idx].row_loc_ids,
				     s_rows[send_idx].data);
	ent_idx++;
      }
      send_idx++;
    }
  }

  if(ent_idx != n_entries){
    PGFEM_printerr("[%d]ERROR: did not reach end of entry list!\n"
	    "%s:%s:%d\n",myrank,__func__,__FILE__,__LINE__);
    PGFEM_Comm_code_abort(mat->mpi_comm,0);
  }

  free(proc_order);
  destroy_entries(n_entries,e);
  return err;
}

int PGFEM_par_matrix_zero_values(PGFEM_par_matrix *mat)
{
  int err = 0;
  int myrank = 0;
  comm_info *send = (comm_info *) mat->send_info;
  comm_info *recv = (comm_info *) mat->recv_info;
  off_proc_entries *s_rows = (off_proc_entries *) mat->s_rows;
  off_proc_entries *r_rows = (off_proc_entries *) mat->r_rows;

  MPI_Comm_rank(mat->mpi_comm,&myrank);
  /* zero local data */
  memset(mat->data,0,mat->n_cols*mat->n_own_rows[myrank]*sizeof(double));

  /* zero send data */
  for(int i=0; i<send->nproc; i++){
    memset(s_rows[i].data,0,mat->n_cols*send->n_info[i]*sizeof(double));
  }

  /* zero receive data */
  for(int i=0; i<recv->nproc; i++){
    memset(r_rows[i].data,0,mat->n_cols*recv->n_info[i]*sizeof(double));
  }

  /* set assembled flag */
  mat->assembled = 1;
  return err;
}

int PGFEM_par_matrix_start_assembly(PGFEM_par_matrix *mat,
				    PGFEM_par_matrix_comm *comm)
{
  int err = 0;

  /* allocate a communication context */
  *comm = PGFEM_calloc(1,sizeof(comm_ctx));
  comm_info *send = (comm_info *) mat->send_info;
  comm_info *recv = (comm_info *) mat->recv_info;
  comm_ctx *COMM = (comm_ctx*) *comm;
  build_comm_ctx(send->nproc,recv->nproc,COMM);

  /* post the recieves */
  off_proc_entries *r_rows = (off_proc_entries*) mat->r_rows;
  for(int i=0; i<recv->nproc; i++){
    int len = recv->n_info[i]*mat->n_cols;
    MPI_Irecv(r_rows[i].data,len,MPI_DOUBLE,recv->proc[i],
	      MPI_ANY_TAG,mat->mpi_comm,&COMM->r_req[i]);
  }

  /* post the sends */
  off_proc_entries *s_rows = (off_proc_entries*) mat->s_rows;
  for(int i=0; i<send->nproc; i++){
    int len = send->n_info[i]*mat->n_cols;
    MPI_Isend(s_rows[i].data,len,MPI_DOUBLE,send->proc[i],
	      i /*tag*/,mat->mpi_comm,&COMM->s_req[i]);
  }

  return err;
}

int PGFEM_par_matrix_end_assembly(PGFEM_par_matrix *mat,
				  PGFEM_par_matrix_comm comm)
{
  int err = 0;
  comm_ctx *COMM = (comm_ctx*) comm;
  comm_info *send = (comm_info *) mat->send_info;
  comm_info *recv = (comm_info *) mat->recv_info;

  /* finish communication */
  MPI_Waitall(recv->nproc,COMM->r_req,COMM->r_stat);
  MPI_Waitall(send->nproc,COMM->s_req,COMM->s_stat);

  if(mat->add_values){
    err += add_assemble_matrix(mat);
  } else {
    err += set_assemble_matrix(mat);
  }

  destroy_comm_ctx(COMM);

  /* set the assembled flag */
  mat->assembled = 1;
  return err;
}

int PGFEM_par_matrix_get_column(const PGFEM_par_matrix *mat,
				const int col_idx,
				double *col)
{
  int err = 0;
  /* return error code if matrix is not assembled or specified column
     is outside dimension of matrix */
  if(!mat->assembled || col_idx >= mat->n_cols) return 1;

  /* compute number of rows/cols */
  int myrank = 0;
  MPI_Comm_rank(mat->mpi_comm,&myrank);
  const int n_loc_row = mat->n_own_rows[myrank];
  const int ncol = mat->n_cols;

  /* copy column */
  for(int i=0; i<n_loc_row; i++){
    col[i] = mat->data[idx_2_gen(i,col_idx,n_loc_row,ncol)];
  }

  return err;
}

int PGFEM_par_vec_dot(const int len,
		      const double *vec_a,
		      const double *vec_b,
		      MPI_Comm comm,
		      double *result)
{
  int err = 0;
  *result = cblas_ddot(len,vec_a,1,vec_b,1);
  err += MPI_Allreduce(MPI_IN_PLACE,result,1,MPI_DOUBLE,MPI_SUM,comm);
  return err;
}

/*===== STATIC FUNCTIONS ====*/
static int PGFEM_par_matvec_compare_int(const void *a, const void *b)
{
  return (*((const int*) a) - *((const int*) b));
}

static int set_local_value(const int row_start,
			   const int n_row,
			   const int n_col,
			   const entry *e,
			   double *data)
{
  int err = 0;
  int loc_row_id = e->row - row_start;
  int idx = idx_2_gen(loc_row_id,e->col,n_row,n_col);
  my_assert(loc_row_id >= 0 && loc_row_id < n_row);
  data[idx] = e->val;
  return err;
}

static int set_off_proc_value(const int row_start,
			      const int n_row,
			      const int n_col,
			      const entry *e,
			      const int *row_ids,
			      double *data)
{
  int err = 0;
  int loc_row_id = e->row - row_start;
  int *row = bsearch(&loc_row_id,row_ids,n_row,sizeof(int),
		     PGFEM_par_matvec_compare_int);
  my_assert(row != NULL);
  int row_id = (row-row_ids);///sizeof(int);
  int idx = idx_2_gen(row_id,e->col,
		      n_row,n_col);
  data[idx] = e->val;
  return err;
}
static int add_to_local_value(const int row_start,
			      const int n_row,
			      const int n_col,
			      const entry *e,
			      double *data)
{
  int err = 0;
  int loc_row_id = e->row - row_start;
  int idx = idx_2_gen(loc_row_id,e->col,n_row,n_col);
  my_assert(loc_row_id >= 0 && loc_row_id < n_row);
  data[idx] += e->val;
  return err;
}

static int add_to_off_proc_value(const int row_start,
				 const int n_row,
				 const int n_col,
				 const entry *e,
				 const int *row_ids,
				 double *data)
{
  int err = 0;
  int loc_row_id = e->row - row_start;
  int *row = bsearch(&loc_row_id,row_ids,n_row,sizeof(int),
		     PGFEM_par_matvec_compare_int);
  my_assert(row != NULL);
  int row_id = (row-row_ids);///sizeof(int);
  int idx = idx_2_gen(row_id,e->col,
		      n_row,n_col);
  data[idx] += e->val;
  return err;
}

static int add_assemble_matrix(PGFEM_par_matrix *mat)
{
  int err = 0;
  int myrank = 0;
  comm_info *recv = (comm_info *) mat->recv_info;
  off_proc_entries *r_rows = (off_proc_entries*) mat->r_rows;

  MPI_Comm_rank(mat->mpi_comm,&myrank);
  int n_col = mat->n_cols;
  int n_row = mat->n_own_rows[myrank];

  for(int i=0; i<recv->nproc; i++){ /* from each proc */
    int n_rec_row = recv->n_info[i];
    for(int j=0; j<n_rec_row; j++){ /* each row on proc */
      int loc_idx = idx_2_gen(r_rows[i].row_loc_ids[j],0,n_row,n_col);
      int rec_idx = idx_2_gen(j,0,n_rec_row,n_col);
      cblas_daxpy(n_col,1.0,&r_rows[i].data[rec_idx],1,
		  &mat->data[loc_idx],1);
    }
  }
  return err;
}

static int set_assemble_matrix(PGFEM_par_matrix *mat)
{
  int err = 0;
  int myrank = 0;
  comm_info *recv = (comm_info *) mat->recv_info;
  off_proc_entries *r_rows = (off_proc_entries*) mat->r_rows;

  MPI_Comm_rank(mat->mpi_comm,&myrank);
  int n_col = mat->n_cols;
  int n_row = mat->n_own_rows[myrank];
  int row_len = n_col*sizeof(double);

  for(int i=0; i<recv->nproc; i++){ /* from each proc */
    int n_rec_row = recv->n_info[i];
    for(int j=0; j<n_rec_row; j++){ /* each row on proc */
      int loc_idx = idx_2_gen(r_rows[i].row_loc_ids[j],0,n_row,n_col);
      int rec_idx = idx_2_gen(j,0,n_rec_row,n_col);
      memcpy(&mat->data[loc_idx],&r_rows[i].data[rec_idx],row_len);
    }
  }
  return err;
}

static inline void build_entry(entry *e){}
static inline void set_entry(const int row,
			     const int col,
			     const double val,
			     entry *e)
{
  e->row = row;
  e->col = col;
  e->val = val;
}

static inline void destroy_entry(entry *e){};
static inline void destroy_entries(const int n_entry,
				   entry *e)
{
  /* for(int i=0; i<n_entry) destroy_entry(e+i); */
  free(e);
}
static int compare_entry_row(const void *a, const void *b)
{
  return (((const entry*) a)->row - ((const entry*) b)->row);
}
static int get_sorted_list_of_entries(const int n_entries,
				      const int *row_idx,
				      const int *col_idx,
				      const double *vals,
				      entry **e)
{
  int err = 0;
  /* push data into entry list and sort by rows */
  *e = PGFEM_calloc(n_entries,sizeof(entry));
  for(int i=0; i<n_entries; i++){
    build_entry(*e+i);
    set_entry(row_idx[i],col_idx[i],vals[i],*e+i);
  }
  qsort(*e,n_entries,sizeof(entry),compare_entry_row);
  return err;
}


static int build_comm_info(const int nproc,
			   const int myrank,
			   const int *proc_comm_info,
			   comm_info *info)
{
  int err = 0;
  info->nproc = 0;
  info->proc = NULL;
  info->n_info = NULL;

  /* exit early if there is no communication to occur */
  if(nproc <= 1) return err;

  for(int i=0; i<nproc; i++){
    if(i == myrank) continue;
    if(proc_comm_info[i] > 0) info->nproc ++;
  }

  info->proc = PGFEM_calloc(info->nproc,sizeof(int));
  info->n_info = PGFEM_calloc(info->nproc,sizeof(int));

  int idx = 0;
  for(int i=0; i<nproc; i++){
    if(i == myrank) continue;
    if(proc_comm_info[i] > 0){
      info->proc[idx] = i;
      info->n_info[idx] = proc_comm_info[i];
      idx++;
    }
  }

  return err;
}

static void destroy_comm_info(comm_info *info)
{
  if(info != NULL){
    free(info->proc);
    free(info->n_info);
  }
  free(info);
}

static int build_off_proc_entries(const int n_rows,
				  const int n_col,
				  off_proc_entries *info)
{
  int err = 0;
  info->row_loc_ids = PGFEM_calloc(n_rows,sizeof(int));
  info->data = PGFEM_calloc(n_rows*n_col,sizeof(double));
  return err;
}

static void destroy_off_proc_entries(off_proc_entries *info)
{
  if(info != NULL){
    free(info->row_loc_ids);
    free(info->data);
  }
}

static int build_comm_ctx(const int n_send,
			  const int n_recv,
			  comm_ctx *COMM)
{
  int err = 0;
  COMM->r_req = NULL;
  COMM->r_stat = NULL;
  if(n_recv > 0){
    COMM->r_req = PGFEM_calloc(n_recv,sizeof(MPI_Request));
    COMM->r_stat = PGFEM_calloc(n_recv,sizeof(MPI_Status));
  }

  COMM->s_req = NULL;
  COMM->s_stat = NULL;
  if(n_send > 0){
    COMM->s_req = PGFEM_calloc(n_send,sizeof(MPI_Request));
    COMM->s_stat = PGFEM_calloc(n_send,sizeof(MPI_Status));
  }
  return err;
}

static void destroy_comm_ctx(comm_ctx *COMM)
{
  if(COMM != NULL){
    free(COMM->r_req);
    free(COMM->s_req);
    free(COMM->r_stat);
    free(COMM->s_stat);
  }
  free(COMM);
  COMM = NULL;
}
static void print_comm_info_graph(const int myrank,
				  comm_info *send,
				  comm_info *recv)
{
  char fname[500];
  comm_info *info = NULL;
  FILE *out = NULL;
  for(int i=0; i<2; i++){
    switch(i){
    case 0: /* send */
      sprintf(fname,"send_info_%0.3d.log",myrank);
      out = PGFEM_fopen(fname,"w");
      info = send;
      break;
    case 1: /* recv */
      sprintf(fname,"recv_info_%0.3d.log",myrank);
      out = PGFEM_fopen(fname,"w");
      info = recv;
      break;
    }
    PGFEM_fprintf(out,"%% FORMAT: %% nproc, rank||procs||count\n");
    PGFEM_fprintf(out,"%% %d\n",info->nproc);
    for(int j=0; j<info->nproc; j++){
      PGFEM_fprintf(out,"%d %d %d\n",myrank,info->proc[j],
		    info->n_info[j]);
    }
    PGFEM_fprintf(out,"\n");
    PGFEM_fclose(out);
  }
}
