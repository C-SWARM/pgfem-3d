/* HEADER */
/**
 * AUTHORS:
 *  Matt Mosby, University of Notre Dame
 */

#include "pgf_fe2_micro_server.h"

#include "utils.h"
#include "PGFEM_io.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

using namespace pgfem3d;
using namespace pgfem3d::net;

void pgf_FE2_micro_server_stats_print(const pgf_FE2_micro_server *server)
{

  const pgf_FE2_micro_server_stats *stat = server->stats;
  PGFEM_printf("SERVER STATS: njob total avg std min max\n");
  PGFEM_printf("%ld %0.5e %0.5e %0.5e %0.5e %0.5e\n\n",
           server->n_jobs,
           stat->total,stat->avg,
           stat->std,stat->min,stat->max);
}

void pgf_FE2_micro_server_init(pgf_FE2_micro_server **server)
{
  static constexpr auto bytes = sizeof(pgf_FE2_micro_server);
  *server = static_cast<pgf_FE2_micro_server*>(malloc(bytes));

  (*server)->n_jobs = 0;
  (*server)->jobs = NULL;
  (*server)->stats = NULL;
}

void pgf_FE2_micro_server_build(pgf_FE2_micro_server *server,
                const pgf_FE2_server_rebalance *rebal)
{
  const size_t keep = pgf_FE2_server_rebalance_n_keep(rebal);
  const size_t recv = pgf_FE2_server_rebalance_n_recv(rebal);
  server->n_jobs = keep + recv;
  
  server->jobs = static_cast<pgf_FE2_job*>(malloc(server->n_jobs *
                                                  sizeof(pgf_FE2_job)));
  server->stats = static_cast<pgf_FE2_micro_server_stats*>(malloc(sizeof(pgf_FE2_micro_server_stats)));

  /* initialize the jobs that do not need rebalancing */
  const int *job_ids_keep = pgf_FE2_server_rebalance_keep_buf(rebal);
  for(size_t i=0; i<keep; i++){
    pgf_FE2_job_init((server->jobs) + i,job_ids_keep[i],FE2_STATE_NEED_INFO);
  }

  /* initialize the jobs that do not need rebalancing */
  const int *job_ids_recv = pgf_FE2_server_rebalance_recv_buf(rebal);
  for(size_t i=0; i<recv; i++){
    pgf_FE2_job_init((server->jobs) + keep + i,job_ids_recv[i],
             FE2_STATE_NEED_INFO_REBALANCE);
  }
}

void pgf_FE2_micro_server_destroy(pgf_FE2_micro_server *server)
{
  for(size_t i=0,e=server->n_jobs; i<e; i++){
    pgf_FE2_job_destroy((server->jobs) + i);
  }
  free(server->jobs);
  free(server->stats);
  free(server);
  server = NULL;
}

/**
 * Check to see if all jobs are finished. Returns true or false.
 */
static int pgf_FE2_micro_server_done(pgf_FE2_micro_server *server,
				     const Microscale *micro)
{
  pgf_FE2_job *__restrict jobs = server->jobs; /* alias */
  for(size_t i=0, n=server->n_jobs; i<n; i++){
    if(pgf_FE2_job_complete(jobs + i,micro) != FE2_STATE_DONE) return 0;
  }
  return 1;
}

/**
 * Attempt to get information from the microscale.
 */
static void pgf_FE2_micro_server_get_info(pgf_FE2_micro_server *server,
					  const MultiscaleComm *mscom,
					  const Microscale *micro)
{
  pgf_FE2_job *__restrict jobs = server->jobs; /* alias */
  for(size_t i=0, n=server->n_jobs; i<n; i++){
    pgf_FE2_job_get_info(jobs + i,mscom,micro);
  }
}

/**
 * Attempt to compute ready jobs.
 */
static void pgf_FE2_micro_server_compute_ready(pgf_FE2_micro_server *server,
                           Microscale *micro,
                           const MultiscaleComm *mscom,
                           const int mp_id)
{
  pgf_FE2_job *__restrict jobs = server->jobs; /* alias */
  for(size_t i=0, n=server->n_jobs; i<n; i++){
    if (mscom->valid_micro_1) {
      pgf_FE2_job_compute(jobs + i,micro,mscom,mp_id,1);
    } else { //ROM
      pgf_FE2_job_compute(jobs + i,micro,mscom,mp_id,2);
    }
  }
}


/**
 * Busy loop looking for message to start the server cycle.
 */
static void pgf_FE2_micro_server_probe_start(const MultiscaleComm *mscom,
					     const Microscale *micro,
					     Status *stat)
{
  ISIRNetwork *net = static_cast<ISIRNetwork*>(micro->net);
  int msg_waiting = 0;
  while (1){
    if (mscom->valid_mm_inter){
      /* exit */
      net->iprobe(NET_ANY_SOURCE,FE2_MICRO_SERVER_EXIT,
		  mscom->mm_inter,&msg_waiting,stat);
      if(msg_waiting) break;
    
      /* rebalance */
      net->iprobe(NET_ANY_SOURCE,FE2_MICRO_SERVER_REBALANCE,
  		mscom->mm_inter,&msg_waiting,stat);
      if(msg_waiting) break;
    }

    if (mscom->valid_mm_inter_ROM){
      /* exit */
      net->iprobe(NET_ANY_SOURCE,FE2_MICRO_SERVER_EXIT,
      mscom->mm_inter_ROM,&msg_waiting,stat);
      if(msg_waiting) break;

      /* rebalance */
      net->iprobe(NET_ANY_SOURCE,FE2_MICRO_SERVER_REBALANCE,
      mscom->mm_inter_ROM,&msg_waiting,stat);
      if(msg_waiting) break;
    }

    /* other ... */
  }

  /* this assert will fail if we accidentally overlap with
     NET_ANY_TAG. FYI: often NET_ANY_TAG = -1.*/
  assert(stat->NET_TAG == FE2_MICRO_SERVER_EXIT
	 || stat->NET_TAG == FE2_MICRO_SERVER_REBALANCE);
}


/**
 * On the MASTER server process look for info from the macroscale and
 * propogate to the workers.
 */
static void pgf_FE2_micro_server_start_cycle(const MultiscaleComm *mscom,
					     const Microscale *micro,
					     pgf_FE2_micro_server *server,
					     pgf_FE2_server_rebalance **rebal,
					     int *exit_server)
{
  ISIRNetwork *net = static_cast<ISIRNetwork*>(micro->net);
  static const size_t n_info = 2;
  long info[n_info]; /* = {0 /\* tag *\/, */
            /*    0 /\* buffer_len *\/}; */
  *exit_server = 0;
  Status stat;

  /* probe for incomming message on comm->mm_inter */
  pgf_FE2_micro_server_probe_start(mscom,micro,&stat);

  /* analyze source of the message and its content */
#ifndef NDEBUG
  int n_micro_proc = 0;
  int n_macro_proc = 0;
  if (mscom->valid_mm_inter) {
    net->comm_size(mscom->worker_inter,&n_micro_proc);
    net->comm_size(mscom->mm_inter,&n_macro_proc);
  } else {
    net->comm_size(mscom->worker_inter_ROM,&n_micro_proc);
    net->comm_size(mscom->mm_inter_ROM,&n_macro_proc);
  }
  n_macro_proc -= n_micro_proc;
  assert(stat.NET_SOURCE < n_macro_proc);
#endif
  int buf_len = 0;
  net->get_status_count(&stat,NET_DT_CHAR,&buf_len);
  info[0] = stat.NET_TAG;
  info[1] = buf_len;

  /* allocate buffer for receive */
  char *buf = static_cast<char*>(malloc(buf_len));

  /* post non-blocking receive matching the probed message */
  Request req;
  if (mscom->valid_mm_inter) {
    net->irecv(buf,buf_len,NET_DT_CHAR,
		    stat.NET_SOURCE,stat.NET_TAG,
		    mscom->mm_inter,&req);
  } else {
    net->irecv(buf,buf_len,NET_DT_CHAR,
        stat.NET_SOURCE,stat.NET_TAG,
        mscom->mm_inter_ROM,&req);
  }


  /* broadcast information to workers !!Signature must match that in
     worker busy loop!!*/
  if (mscom->valid_mm_inter) {
    net->bcast(info,n_info,NET_DT_LONG,0,mscom->micro);
  } else {
    net->bcast(info,n_info,NET_DT_LONG,0,mscom->micro_ROM);
  }
  /* complete communication w/ macroscale */
  net->wait(&req,NET_STATUS_IGNORE);

  switch(stat.NET_TAG){
  case FE2_MICRO_SERVER_EXIT:
    /* set flag and move on to exit function */
    *exit_server = 1;

    /* free the buffer */
    free(buf);
    break;

 case FE2_MICRO_SERVER_REBALANCE:
    /* broadcast information to microscale (could change to
       non-blocking?) */
  if (mscom->valid_mm_inter) {
    net->bcast(buf,buf_len,NET_DT_CHAR,0,mscom->micro);
  } else {
    net->bcast(buf,buf_len,NET_DT_CHAR,0,mscom->micro_ROM);
  }  
   /* build rebalance */
   pgf_FE2_server_rebalance_build_from_buffer(rebal,(void**) &buf);
   
   /* build server list */
   pgf_FE2_micro_server_build(server,*rebal);
   
   /* do not free the buffer. Maintain access through rebal. Buffer
      is free'd when rebal is destroyed */
   break;
  default: assert(0); /* should never get here */ break;
  }
}

static void pgf_FE2_micro_server_compute_stats(pgf_FE2_micro_server *server)
{
  /* aliases */
  pgf_FE2_micro_server_stats *s = server->stats;
  pgf_FE2_job *j = server->jobs;

  /* initialization */
  const size_t njob = server->n_jobs;
  s->total = s->avg = s->std = 0.0;

  /* total, min, max, and average */
  double total = 0.0;
  double min = INT_MAX;
  double max = 0.0;
  for(size_t i=0; i<njob; i++){
    double cur = (double) j[i].time;
    total += cur;
    if(min > cur) min = cur;
    if(max < cur) max = cur;
  }
  const double avg = total / njob;

  /* standard deviation */
  double std = 0.0;
  for(size_t i=0; i<njob; i++){
    double cur = (double) j[i].time;
    std += (cur-avg)*(cur-avg);
  }
  s->total = total;
  s->avg = avg;
  s->std = sqrt(1./njob*std);
  s->min = min;
  s->max = max;
}

/**
 * Make a fairly shallow copy of the job data into a communication
 * buffer.
 */
static void pgf_FE2_micro_server_pack_summary(pgf_FE2_micro_server *server,
					      size_t *buf_len,
					      char **buf)
{
  /* compute stats */
  pgf_FE2_micro_server_compute_stats(server);

  /* compute buffer length and allocate */
  const size_t n_jobs = server->n_jobs;
  static const size_t size_job = sizeof(*(server->jobs));
  static const size_t size_stats = sizeof(*(server->stats));
  *buf_len = sizeof(server->n_jobs) + n_jobs*size_job + size_stats;
  (*buf) = static_cast<char*>(malloc(*buf_len));

  /* pack the server job summary */
  size_t pos = 0;
  pack_data(&n_jobs,*buf,&pos,1,sizeof(n_jobs));
  pack_data(server->stats,*buf,&pos,1,size_stats);
  pack_data(server->jobs,*buf,&pos,n_jobs,size_job);
}

static void pgf_FE2_micro_server_finish_cycle(const MultiscaleComm *mscom,
					      const Microscale *micro,
					      pgf_FE2_micro_server *server)
{
  ISIRNetwork *net = static_cast<ISIRNetwork*>(micro->net);
  int n_micro_proc = 0;
  int n_macro_proc = 0;
  if (mscom->valid_micro_1) {
    net->comm_size(mscom->worker_inter,&n_micro_proc);
    net->comm_size(mscom->mm_inter,&n_macro_proc);
  } else {
    net->comm_size(mscom->worker_inter_ROM,&n_micro_proc);
    net->comm_size(mscom->mm_inter_ROM,&n_macro_proc);
  }

  n_macro_proc -= n_micro_proc;

  Request *req;
  net->allocRequestArray(n_macro_proc, &req);
  
  char *buf = NULL;
  size_t buf_len = 0;
  pgf_FE2_micro_server_pack_summary(server,&buf_len,&buf);
  if (mscom->valid_micro_1) {
    for (int i=0; i<n_macro_proc; i++) {
      net->isend(buf, buf_len, NET_DT_CHAR, i, FE2_MICRO_SERVER_REBALANCE,
		      mscom->mm_inter, req+i);
    }
  } else {
    for (int i=0; i<n_macro_proc; i++) {
      net->isend(buf, buf_len, NET_DT_CHAR, i, FE2_MICRO_SERVER_REBALANCE,
          mscom->mm_inter_ROM, req+i);
    }
  }
  net->waitall(n_macro_proc,req,NET_STATUS_IGNORE);

//  net->waitall(1,req,NET_STATUS_IGNORE);
  free(buf);
  delete [] req;
}

void
pgf_FE2_micro_server_unpack_summary(pgf_FE2_micro_server **Server,
                                    const char *buf)
{
  pgf_FE2_micro_server_init(Server);
  pgf_FE2_micro_server *server = *Server; /* alias */
  static constexpr size_t size_jobs = sizeof(*(server->jobs));
  static constexpr size_t size_stats = sizeof(*(server->stats));

  server->stats = static_cast<pgf_FE2_micro_server_stats*>(malloc(size_stats));
  size_t pos = 0;
  unpack_data(buf,&(server->n_jobs),&pos,1,sizeof(server->n_jobs));
  server->jobs = static_cast<pgf_FE2_job*>(malloc(server->n_jobs*size_jobs));
  unpack_data(buf,server->stats,&pos,1,size_stats);
  unpack_data(buf,server->jobs,&pos,server->n_jobs,size_jobs);
  for(int i=0,e=server->n_jobs; i<e; i++){
    server->jobs[i].comm_buf = NULL;
  }
}

static int
pgf_FE2_micro_server_master(const MultiscaleComm *mscom,
                            Microscale *micro,
                            const int mp_id)
{
  int err = 0;
  int exit_server = 0;

  while(1){
    /* begin server loop */
    pgf_FE2_micro_server *server = NULL;
    pgf_FE2_server_rebalance *rebal = NULL;

    pgf_FE2_micro_server_init(&server);
    pgf_FE2_micro_server_start_cycle(mscom,micro,server,&rebal,&exit_server);
    if(exit_server){
      pgf_FE2_micro_server_destroy(server);
      pgf_FE2_server_rebalance_destroy(rebal);
      break;
    }

    /* swap microscale information according to rebal */
    /* Future implementation/improvements will overlay this with
       communication of job information from macroscale and computation
       of microstructures that are not to be rebalanced. */
    pgf_FE2_server_rebalance_post_exchange(rebal,mscom,micro);
    /* blocks until all exchanges are done. Could improve
       async. comm here. */
    pgf_FE2_server_rebalance_finalize_exchange(rebal,mscom,micro);

    /* for now, go ahead and mark all jobs as needing information */
    for(size_t i=0,e=server->n_jobs; i<e; i++){
      pgf_FE2_job_set_state(server->jobs + i,FE2_STATE_NEED_INFO);
    }

    while(!pgf_FE2_micro_server_done(server,micro)) {
      /* Get info from macroscale */
      pgf_FE2_micro_server_get_info(server,mscom,micro);

      /* in improved implementation, advance rebalancing here. I.e.,
     check for rebalancing jobs where the communication has
     completed and update their state appropriately. */

      /* compute the jobs that are ready. Posts sends */
      pgf_FE2_micro_server_compute_ready(server,micro,mscom,mp_id);

    }

    /* send server/job statistics to macroscale for rebalancing */
    pgf_FE2_micro_server_finish_cycle(mscom,micro,server);
    pgf_FE2_micro_server_stats_print(server);

    pgf_FE2_micro_server_destroy(server);
    pgf_FE2_server_rebalance_destroy(rebal);
  }

  return err;
}

/**
 * Worker busy loop. Some logic as to which jobs to initiate.
 */
static int pgf_FE2_micro_server_worker(const MultiscaleComm *mscom,
				       Microscale *micro,
				       const int mp_id,
               int micro_model)
{
  int err = 0;
  int exit_server = 0;
  static const size_t n_info = 2;
  long info[n_info]; /* = {0 /\*id*\/, */
             /*   0 /\*msg_size in bytes *\/}; */
  while( !exit_server ){
    /* get information from master on server */
    micro->net->bcast(info,n_info,NET_DT_LONG,0,micro->comm);

    switch(info[0]){
    case FE2_MICRO_SERVER_EXIT:
      exit_server = 1;
      break;

    case FE2_MICRO_SERVER_REBALANCE:
      /* rebalance on workers */
      {
	/* allocate buffer for receive */
	int buf_len = info[1];
	char *buf = static_cast<char*>(malloc(buf_len));
	pgf_FE2_server_rebalance *rebal = NULL;
	
	/* get information from MASTER */
  if (mscom->valid_micro_1) {
	  micro->net->bcast(buf,buf_len,NET_DT_CHAR,0,mscom->micro);
	} else {
    micro->net->bcast(buf,buf_len,NET_DT_CHAR,0,mscom->micro_ROM);// should it also be bcast_ROM?
  }
	/* build rebalancing data structure */
	pgf_FE2_server_rebalance_build_from_buffer(&rebal,(void**) &buf);
	
	/* perform rebalancing */
	pgf_FE2_server_rebalance_post_exchange(rebal,mscom,micro);
	/* blocks until all exchanges are done. Could improve
	   async. comm here. */
	pgf_FE2_server_rebalance_finalize_exchange(rebal,mscom,micro);
	
	/* implicitly free's buf */
	pgf_FE2_server_rebalance_destroy(rebal);
      }
      break;
      
    default: /* valid job, compute work */
        if (mscom->valid_micro_1) {
          pgf_FE2_job_compute_worker(info[0],info[1],micro,mp_id,1);
        } else { //ROM
          pgf_FE2_job_compute_worker(info[0],info[1],micro,mp_id,2);
        }
      break;
    }
  }

  return err;
}

int pgf_FE2_micro_server_START(const MultiscaleComm *mscom,
			       Microscale *micro,
			       const int mp_id)
{
  int err = 0;
//  assert(mscom->valid_micro);
  if (mscom->valid_mm_inter || mscom->valid_mm_inter_ROM){
   err += pgf_FE2_micro_server_master(mscom,micro,mp_id);  
  } else {
    if (mscom->valid_micro_1) {
      err += pgf_FE2_micro_server_worker(mscom,micro,mp_id,1);
    } else { //ROM
      err += pgf_FE2_micro_server_worker(mscom,micro,mp_id,2);
    }
  }
  return err;
}
