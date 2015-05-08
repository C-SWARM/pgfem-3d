/* HEARER */
#include "compute_ms_cohe_job.h"
#include <math.h>
#include <string.h>
#include <assert.h>
#include "mkl_cblas.h"

#include "vtk_output.h"
#include "allocation.h"
#include "enumerations.h"
#include "utils.h"
#include "index_macros.h"
#include "Newton_Raphson.h"
#include "matice.h"
#include "PGFEM_par_matvec.h"
#include "stiffmat_fd.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "displacement_based_element.h"
#include "interface_macro.h"
#include "solve_system.h"
#include "pgf_fe2_restart.h"

#ifndef JOB_LOGGING
#define JOB_LOGGING 1
#endif

static const int ndim = 3;

/*==== STATIC FUNCTION PROTOTYPES ====*/

/** Wrapper for Newton Raphson. */
static int ms_cohe_job_nr(COMMON_MICROSCALE *c,
			  MICROSCALE_SOLUTION *s,
			  const PGFem3D_opt *opts,
			  int *n_step);

/** Set the job supports appropriately from the jump at (n) and
    (n+1). Also set the normal to the interface. */
static int set_job_supports(const MS_COHE_JOB_INFO *p_job,
			    SUPP sup);


/** compute and initialize the dense matrix structure for the mixed
    tangents */
static int initialize_ms_cohe_job_mixed_tangents
(const int n_cols,
 const COMMON_MICROSCALE *common);

/** Call stiffmat_fd using information for the microscale job. The
    microscale tangent reset and assembled. Contains (blocking)
    collective communication */ 
static int ms_cohe_job_compute_micro_tangent(COMMON_MICROSCALE *c,
					     MICROSCALE_SOLUTION *s,
					     PGFem3D_opt *o);

/** Compute all of the microscale terms for the macroscale tangent and
    residual. The mixed tangents are assembled and the cohesive
    tangent at the macroscale is computed and distributed to all
    processors. Note that only DISP analysis type is currently
    supported. */
static int compute_ms_cohe_job_micro_terms(const MS_COHE_JOB_INFO *job,
					   const COMMON_MICROSCALE *c,
					   const MICROSCALE_SOLUTION *s,
					   const PGFem3D_opt *o,
					   PGFEM_par_matrix *K_01,
					   PGFEM_par_matrix *K_10,
					   double *K_00_contrib);

/** Compute the contibutions to the microscale terms from a single
    element */
static int compute_elem_micro_terms(const int elem_id,
				    const MS_COHE_JOB_INFO *job,
				    const COMMON_MICROSCALE *c,
				    const MICROSCALE_SOLUTION *s,
				    const PGFem3D_opt *o,
				    PGFEM_par_matrix *K_01,
				    PGFEM_par_matrix *K_10,
				    double *K_00_contrib);

/** Compute and assemble the contribution to the macroscale
    tangent. Involves multiple linear solves of the microscale */
static int compute_ms_cohe_job_tangent(const int macro_ndof,
				       const long *macro_dof_ids,
				       const PGFEM_par_matrix *K_01,
				       const PGFEM_par_matrix *K_10,
				       const PGFem3D_opt *o,
				       const int print_level,
				       COMMON_MICROSCALE *c,
				       double *K_00_contrib);

static int print_ms_cohe_job(const MS_COHE_JOB_INFO *job,
			     const COMMON_MICROSCALE *c,
			     const MICROSCALE_SOLUTION *s,
			     const PGFem3D_opt *o);

/**
 * Update the job information from n <-- n+1.
 *
 * Modifies the job, particularly updates traction_n, jump_n,
 * max_traction, and max_jump.
 *
 * \return non-zero if cell failure is detected.
 */
static int update_job_information(MS_COHE_JOB_INFO *job)
{
  /* hard-code scaling for failure detection */
  static const double small_val = 0.05;

  int cell_failure_detected = 0;
  const double tn = cblas_dnrm2(ndim,job->traction_n,1);
  const double tnp1 = cblas_dnrm2(ndim,job->traction,1);
  const double jnp1 = cblas_dnrm2(ndim,job->jump,1);

  /* LOADING CONDITION: eff. jump is greater than any previous */
  if (jnp1 > job->max_jump) {
    /* update state variable */
    job->max_jump = jnp1;

    /* update max traction */
    if (tnp1 > job->max_traction) job->max_traction = tnp1;

    /* 
     * FAILURE CONDITION: traction is small and decreasing under
     * LOADING CONDITIO
     */
    if( (tnp1 < (small_val * job->max_traction)) && (tnp1 < tn) ){
      cell_failure_detected++;
    }
  }

  /* UPDATE: n <-- n+1 */
  memcpy(job->jump_n,job->jump,ndim*sizeof(double));
  memcpy(job->traction_n,job->traction,ndim*sizeof(double));

  /* RESET: traction, residual and tangent */
  memset(job->traction,0,ndim*sizeof(double));
  memset(job->traction_res,0,job->ndofe*sizeof(double));
  memset(job->K_00_contrib,0,job->ndofe*job->ndofe*sizeof(double));

  return cell_failure_detected;
}

/*==== API FUNCTION DEFINITIONS ====*/

/** Given a microscale job, compute the solution on the time
    increment, store the traction and assemble the tangent */
int compute_ms_cohe_job(const int job_id,
			MS_COHE_JOB_INFO *p_job,
			MICROSCALE *microscale)
{
  int err = 0;
  const int print_level = 1;
  COMMON_MICROSCALE *common = microscale->common;
  MICROSCALE_SOLUTION *sol = microscale->sol + job_id;

  int myrank = 0;
  MPI_Comm_rank(common->mpi_comm,&myrank);
  if(myrank == 0){
    PGFEM_printf("=== MICROSCALE cell %d of %d ===\n",
		 job_id+1,microscale->idx_map.size);
  }

  /* switch set up job */
  switch(p_job->job_type){
  case JOB_UPDATE: case JOB_PRINT:/* do nothing */ break;

  default: /* reset state to macro time (n) */
    /* reset the microscale solution to macro time (n) */
    err += reset_MICROSCALE_SOLUTION(sol,microscale);

    /* reset the supports to contain the previous jump and current
       increment */
    err += set_job_supports(p_job,common->supports);
    break;
  }

  /* copy the solve time from the job */
  sol->p_tim = p_job->tim;
  sol->tim = (p_job->tim > 0);
  memcpy(sol->times,p_job->times,3*sizeof(double));
  sol->dt = sol->times[sol->tim + 1] - sol->times[sol->tim];

  /* swtich compute job */
  switch(p_job->job_type){
  case JOB_COMPUTE_EQUILIBRIUM:
    if(sol->failed) break;
    /* print time step information */
    if(myrank == 0){
      PGFEM_printf("=== EQUILIBRIUM SOLVE ===\n");
      PGFEM_printf("\nFinite deformations time step %ld) "
		   " Time %f | dt = %10.10f\n",
		   p_job->tim,sol->times[sol->tim+1],sol->dt);
    }

    /* compute the microscale equilibrium. */
    err += ms_cohe_job_nr(common,sol,microscale->opts,&(p_job->n_step));

    /*=== INTENTIONAL DROP THROUGH ===*/
  case JOB_NO_COMPUTE_EQUILIBRIUM:
    if(sol->failed) break;
    /* Do not compute the microscale equilibrium, just assemble the
       microscale tangent */
    {

      if(JOB_LOGGING && myrank == 0) PGFEM_printf("=== MICROSCALE TANGENT ===\n");
      /*  compute the equilibriated microscale tangent */
      err += ms_cohe_job_compute_micro_tangent(common,sol,microscale->opts);
      if(myrank == 0) PGFEM_printf("\n");

      static int init_mixed = 0;
      if(!init_mixed){
	/* initialize matrices K_01 and K_10 */
	/* allocate for quad macro elements */
	const int n_cols = 8*ndim;
	/* const int n_cols = p_job->ndofe; */
      if(JOB_LOGGING && myrank == 0) PGFEM_printf("=== INIT MIXED TANGENT ===\n");
	err += initialize_ms_cohe_job_mixed_tangents(n_cols,common);
	init_mixed++;
      }

      /* get pointers to matrix objects and number of active columns */
      const int n_cols = p_job->ndofe;
      PGFEM_par_matrix *K_01 = (PGFEM_par_matrix *) common->K_01;
      PGFEM_par_matrix *K_10 = (PGFEM_par_matrix *) common->K_10;
      double *K_00_contrib = p_job->K_00_contrib;

      if(JOB_LOGGING && myrank == 0) PGFEM_printf("=== COMPUTE MICRO TERMS ===\n");
      /* compute microscale terms */
      err += compute_ms_cohe_job_micro_terms(p_job,common,sol,
					     microscale->opts,
					     K_01,K_10,K_00_contrib);

      if(JOB_LOGGING && myrank == 0) PGFEM_printf("=== MACROSCALE TANGENT/SCHUR ===\n");
      /* compute macroscale tangent contribution */
      err += compute_ms_cohe_job_tangent(n_cols,p_job->g_dof_ids,K_01,K_10,
					 microscale->opts,print_level,
					 microscale->common,K_00_contrib);
    }
    break;

  case JOB_UPDATE:
    /* update job information and set cell failure condition. */
    {
      int failure_detected = update_job_information(p_job);
      if (!sol->failed && failure_detected) sol->failed = 1;
    }

    /* update the solution and state variables n <- n+1 */
    err += update_MICROSCALE_SOLUTION(sol,microscale);

    /* Set the supports for the new n-state and null the increment */
    err += set_job_supports(p_job,common->supports);

    break;

  case JOB_PRINT:
    /* output the job based on the print flag to the file specified by
       the options and solution step id */
    err += print_ms_cohe_job(p_job,common,sol,microscale->opts);

    /* print the restart file regardless of output parameters */
    {
      int cell_id = sol_idx_map_idx_get_id(&(microscale->idx_map),job_id);
      assert(cell_id >= 0);
      err += pgf_FE2_restart_print_micro(microscale,cell_id);
    }					 
    break;

  case JOB_EXIT: /* do nothing */ break;
  default:
    PGFEM_printerr("Have not implemented these job types yet!\n");
    PGFEM_Abort();
    break;
  }

  if(myrank == 0){
    print_MS_COHE_JOB_INFO(PGFEM_stdout,p_job);
    PGFEM_printf("=====================================\n\n");
  }

  return err;
}/* compute_ms_cohe_job() */

int assemble_ms_cohe_job_res(const int job_id,
			     const MS_COHE_JOB_INFO *p_job,
			     const MPI_Comm micro_comm,
			     const MPI_Comm macro_comm,
			     double *loc_res)
{
  int err = 0;
  int micro_rank = 0;
  int macro_rank = 0;
  MPI_Comm_rank(micro_comm,&micro_rank);
  MPI_Comm_rank(macro_comm,&macro_rank);

  /* exit early without doing anything if not the owning domain */
  if(macro_rank != p_job->proc_id) return err;

  for(int i=0; i<p_job->nnode*ndim; i++){
    const int idx = p_job->loc_dof_ids[i] - 1;
    if(idx >= 0) loc_res[idx] += p_job->traction_res[i];
  }

  return err;
}

/*=====================================*/
/*==== STATIC FUNCTION DEFINITIONS ====*/
/*=====================================*/

static int ms_cohe_job_nr(COMMON_MICROSCALE *c,
			  MICROSCALE_SOLUTION *s,
			  const PGFem3D_opt *opts,
			  int *n_step)
{
  int err = 0;
  double time = 0.0;
  double nl_err = c->lin_err;
  int full_NR = 1; /* 0 is modified NR */
  double pores = 0.0;
  const int print_level = 0;
  *n_step = 0;

  /* copy of load increment */
  double *sup_defl = PGFEM_calloc(c->supports->npd,sizeof(double));
  memcpy(sup_defl,c->supports->defl_d,c->supports->npd*sizeof(double));

  time += Newton_Raphson(print_level,n_step,c->ne,0,c->nn, c->ndofn,
			 c->ndofd,c->npres,s->tim,s->times,
			 nl_err,s->dt,c->elem,NULL,
			 c->node,c->supports,sup_defl,c->hommat,
			 c->matgeom,s->sig_e,s->eps,c->Ap,
			 c->Ai,s->r,s->f,s->d_r,
			 s->rr,s->R,s->f_defl,s->RR,
			 s->f_u,s->RRn,s->crpl,opts->stab,
			 c->nce,c->coel,full_NR,&pores,
			 c->SOLVER,s->BS_x,s->BS_f,s->BS_RR,
			 0.0,0.0,0.0,c->lin_err,
			 s->BS_f_u,c->DomDof,c->pgfem_comm,c->GDof,
			 1,c->maxit_nl,&s->NORM,c->nbndel,
			 c->bndel,c->mpi_comm,c->VVolume,opts,NULL, 0, NULL, NULL);
	
  free(sup_defl);
  return err;
}/* ms_cohe_job_nr() */

static int set_job_supports(const MS_COHE_JOB_INFO *p_job,
			    SUPP sup)
{
  /* NOTE: This is written for TOTAL LAGRANGIAN */
  int err = 0;
  if(!sup->multi_scale || sup->npd != 6) err++;
  memcpy(sup->N0,p_job->normal,ndim*sizeof(double));
  for(int i=0; i<ndim; i++){
    sup->defl[i] = 0.5*p_job->jump_n[i];
    sup->defl[i+ndim] = -0.5*p_job->jump_n[i];
    sup->defl_d[i] = 0.5*(p_job->jump[i] - p_job->jump_n[i]);
    sup->defl_d[i+ndim] = -0.5*(p_job->jump[i] - p_job->jump_n[i]);
  }

  /* compute the macroscopic deformation gradient. This is
     particularly important for if microscale equilibrium is not
     computed. */
  err += compute_macro_grad_u(sup->F0,sup,DISP);
  return err;
}/* set_job_supports() */

static int initialize_ms_cohe_job_mixed_tangents
(const int n_cols,
 const COMMON_MICROSCALE *common)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  if(!(common->K_01 == NULL && common->K_10 == NULL)) return err; 
  MPI_Comm_rank(common->mpi_comm,&myrank);
  MPI_Comm_size(common->mpi_comm,&nproc);

  /* compute global number of rows/columns and how many are owned on
     this process */
  int n_rows = 0;
  for(int i=0; i<nproc; i++){
    n_rows += common->DomDof[i];
  }
  const int n_own_rows = common->DomDof[myrank];

  /* compute the number of entries and list the global row
     numbers. NOTE: only nodal dofs are counted currently */
  const int n_entries = common->ndofn*common->nn;
  int *rows = PGFEM_calloc(n_entries,sizeof(int));
  int *cols = PGFEM_calloc(n_entries,sizeof(int));

  /* get the row numbers (G dof ids). If the dof is a BC, set to
     -1. This makes the duplicate counting work properly in the
     initialization */
  {
    int idx = 0;
    for(int i=0; i<common->nn; i++){
      for(int j=0; j<common->ndofn; j++){
	rows[idx] = common->node[i].Gid[j] - 1;
	if(rows[idx] < 0) rows[idx] = -1;
	idx++;
      }
    }
  }
  /* sort list of entries so multiple calls to
     initialize_PGFEM_par_matrix will be faster */
  qsort(rows,n_entries,sizeof(int),compare_int);

  err += initialize_PGFEM_par_matrix(n_rows,n_cols,n_own_rows,
				     n_entries,rows,cols,
				     common->mpi_comm,
				     (PGFEM_par_matrix**) &common->K_01);

  err += initialize_PGFEM_par_matrix(n_rows,n_cols,n_own_rows,
				     n_entries,rows,cols,
				     common->mpi_comm,
				     (PGFEM_par_matrix**) &common->K_10);
  free(rows);
  free(cols);
  return err;
}

static int ms_cohe_job_compute_micro_tangent(COMMON_MICROSCALE *c,
					     MICROSCALE_SOLUTION *s,
					     PGFem3D_opt *o)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  double nor_min = c->lin_err;
  MPI_Comm_rank(c->mpi_comm,&myrank);
  MPI_Comm_size(c->mpi_comm,&nproc);

  /* reset the microscale tangent to zeros */
  ZeroHypreK(c->SOLVER,c->Ai,c->DomDof[myrank]);

  /* assemble to the microscale tangent matrix */
  err += stiffmat_fd(c->Ap,c->Ai,c->ne,0,c->ndofn,
		     c->elem,NULL,c->nbndel,c->bndel,
		     c->node,c->hommat,c->matgeom,s->sig_e,
		     s->eps,s->d_r,s->r,c->npres,c->supports,
		     /*iter*/0,nor_min,s->dt,s->crpl,o->stab,
		     c->nce,c->coel,0,0.0,s->f_u,myrank,nproc,
		     c->DomDof,c->GDof,c->pgfem_comm,c->mpi_comm,
		     c->SOLVER,o,0,NULL,NULL);

  /* finalize the microscale tangent matrix assembly */
  err += HYPRE_IJMatrixAssemble(c->SOLVER->hypre_k);

  return err;
}

static int compute_ms_cohe_job_micro_terms(const MS_COHE_JOB_INFO *job,
					   const COMMON_MICROSCALE *c,
					   const MICROSCALE_SOLUTION *s,
					   const PGFem3D_opt *o,
					   PGFEM_par_matrix *K_01,
					   PGFEM_par_matrix *K_10,
					   double *K_00_contrib)
{
  int err = 0;

  /* error check analysis type */
  switch(o->analysis_type){
  case DISP: break;
  default:
    PGFEM_printerr("ERROR: analysis type %d is not supported!\n"
	    "%s:%s:%d",o->analysis_type,__func__,__FILE__,__LINE__);
    PGFEM_Comm_code_abort(c->mpi_comm,0);
  }

  /* zero the matrices/vectors */
  const int macro_ndof = job->ndofe;
  err += PGFEM_par_matrix_zero_values(K_01);
  err += PGFEM_par_matrix_zero_values(K_10);
  nulld(K_00_contrib,macro_ndof*macro_ndof);
  nulld(job->traction_res,macro_ndof);
  nulld(job->traction,ndim);


  /* volume elements */
  for(int i=0; i<c->ne; i++){
    err += compute_elem_micro_terms(i,job,c,s,o,K_01,
				    K_10,K_00_contrib);
  }

  /* cohesive elements */
  for(int i=0; i<c->nce; i++){
    /* compute and assemble cohesive element contributions */

    /*=== NOT CURRENTLY IMPLEMENTED ===*/
    /* need to formulate fully-coupled multiscale with contributions
       from cohesive elements at microscale */
  }

  /*=== ASSEMBLE MATRICES from all domains ===*/

  /* start assembly of the distributed matrices */
  PGFEM_par_matrix_comm K_01_comm = NULL;
  PGFEM_par_matrix_comm K_10_comm = NULL;
  err += PGFEM_par_matrix_start_assembly(K_01,&K_01_comm);
  err += PGFEM_par_matrix_start_assembly(K_10,&K_10_comm);

  /* assemble K_00_contrib from all domains */
  err += MPI_Allreduce(MPI_IN_PLACE,K_00_contrib,
		       job->ndofe*job->ndofe,
		       MPI_DOUBLE,MPI_SUM,c->mpi_comm);

  /* assemble traction residual from all domains */
  err += MPI_Allreduce(MPI_IN_PLACE,job->traction_res,job->ndofe,
		       MPI_DOUBLE,MPI_SUM,c->mpi_comm);

  /* assemble traction from all domains */
  err += MPI_Allreduce(MPI_IN_PLACE,job->traction,ndim,
		       MPI_DOUBLE,MPI_SUM,c->mpi_comm);

  /* finish assembly of the distributed matrices */
  err += PGFEM_par_matrix_end_assembly(K_01,K_01_comm);
  err += PGFEM_par_matrix_end_assembly(K_10,K_10_comm);

  return err;
}

static int compute_elem_micro_terms(const int elem_id,
				    const MS_COHE_JOB_INFO *job,
				    const COMMON_MICROSCALE *c,
				    const MICROSCALE_SOLUTION *s,
				    const PGFem3D_opt *o,
				    PGFEM_par_matrix *K_01,
				    PGFEM_par_matrix *K_10,
				    double *K_00_contrib)
{
  int err = 0;
  int myrank = 0;
  MPI_Comm_rank(c->mpi_comm,&myrank);

  /* important macro quantities */
  const int macro_nnode = job->nnode;
  const int macro_ndofn = ndim;
  const int macro_ndof = job->ndofe;
  const double *macro_shape_func = job->shape;
  const double *macro_normal = job->normal;
  const double macro_int_wt = job->int_wt;
  const double layer_thickness = c->supports->lc;

  /* microscale element information */
  const ELEMENT *elem = c->elem + elem_id;
  const int ndofn = c->ndofn;
  const int nne = elem->toe;
  const long *node_ids = elem->nod;
  const int ndofe = get_ndof_on_elem_nodes(nne,node_ids,c->node);

  /* allocate microscale information */
  long *local_dof_ids = PGFEM_calloc(ndofe,sizeof(long));
  long *global_dof_ids = PGFEM_calloc(ndofe,sizeof(long));
  double *K_01_e = PGFEM_calloc(macro_ndof*ndofe,sizeof(double));
  double *K_10_e = PGFEM_calloc(ndofe*macro_ndof,sizeof(double));
  double *K_00_e = PGFEM_calloc(macro_ndof*macro_ndof,sizeof(double));
  double *traction_res_e = PGFEM_calloc(macro_ndof,sizeof(double));
  double *traction_e = PGFEM_calloc(ndim,sizeof(double));
  double *disp = PGFEM_calloc(ndofe,sizeof(double));
  double *x = PGFEM_calloc(nne,sizeof(double));
  double *y = PGFEM_calloc(nne,sizeof(double));
  double *z = PGFEM_calloc(nne,sizeof(double));


  /* microscale dof ids */
  get_dof_ids_on_elem_nodes(0,nne,ndofn,node_ids,c->node,local_dof_ids);
  get_dof_ids_on_elem_nodes(1,nne,ndofn,node_ids,c->node,global_dof_ids);

  /* microscale node coordinates and displacements */
  switch(o->analysis_type){
  case DISP:
    nodecoord_total(nne,node_ids,c->node,x,y,z);
    /*==== microscale displacement ====*/
    /* The increment of displacement is zero after equilibrium,
       therefore displacement is zero except for total Lagrangian
       formulation. */
    def_elem(local_dof_ids,ndofe,s->r,NULL,c->node,disp,c->supports,1);
    break;
  default: nodecoord_updated(nne,node_ids,c->node,x,y,z); break;
  }

  /* compute the microscale contributions based on analysis type */
  switch(o->analysis_type){
  case DISP:
    err += DISP_cohe_micro_terms_el(K_00_e,K_01_e,K_10_e,traction_res_e,
				    traction_e,
				    macro_nnode,macro_ndofn,macro_int_wt,
				    macro_shape_func,macro_normal,
				    layer_thickness,c->VVolume,
				    elem_id,ndofn,nne,x,y,z,c->elem,
				    c->hommat,node_ids,c->node,
				    s->eps,s->sig_e,c->supports,disp);
    if(err){
      PGFEM_printerr("[%d]ERROR: DISP_cohe_micro_terms_el returned error status!\n",myrank);
      PGFEM_Abort();
    }
    break;
  default:
    /*=== NOT IMPLEMENTED ===*/
    PGFEM_printerr("ERROR: multiscale cohesive modeling not "
	    "yet implemented for this analysis type (%d)!\n"
	    "%s:%s:%d\n",o->analysis_type,__func__,__FILE__,__LINE__);
    PGFEM_Abort();
    break;
  }

  /*=== ASSEMBLY ===*/
  /* K_00_contrib. Add Contribution from LOCAL element */
  cblas_daxpy(macro_ndof*macro_ndof,1.0,K_00_e,1,K_00_contrib,1);

  /* traction.Add conttribution from LOCAL element */
  cblas_daxpy(macro_ndof,1.0,traction_res_e,1,job->traction_res,1);
  cblas_daxpy(ndim,1.0,traction_e,1,job->traction,1);

  /* K_01 & K_10 */
  /* { */
  /*   for(int a=0; a<macro_nnode; a++){ */
  /*     for(int b=0; b<macro_ndofn; b++){ */
  /* 	for(int w=0; w<nne; w++){ */
  /* 	  for(int g=0; g<ndim; g++){ */
  /* 	    /\* column idx is simply the index of the macro dof *\/ */
  /* 	    const int col_idx = idx_2_gen(a,b,macro_nnode,macro_ndofn); */

  /* 	    /\* row idx is the global dof id *\/ */
  /* 	    const int row_idx = global_dof_ids[idx_2_gen(w,g,nne,ndim)] - 1; */

  /* 	    /\* skip if boundary condition *\/  */
  /* 	    if(row_idx < 0) continue; */

  /* 	    /\* get indices *\/ */
  /* 	    const int idx_01 = idx_K_gen(a,b,w,g,macro_nnode, */
  /* 					 macro_ndofn,nne,ndim); */
  /* 	    const int idx_10 = idx_K_gen(w,g,a,b,nne,ndim, */
  /* 					 macro_nnode,macro_ndofn); */

  /* 	    /\* add values to distributed matrix *\/ */
  /* 	    err += PGFEM_par_matrix_add_to_values(1,&row_idx,&col_idx, */
  /* 						  K_01_e + idx_01,K_01); */
  /* 	    err += PGFEM_par_matrix_add_to_values(1,&row_idx,&col_idx, */
  /* 						  K_10_e + idx_10,K_10); */
  /* 	  } */
  /* 	} */
  /*     } */
  /*   } */
  /* } */

  {
    /* allocate enough space for full matrix */
    double *val_01 = PGFEM_calloc(macro_nnode*macro_ndofn*nne*ndim,sizeof(double));
    double *val_10 = PGFEM_calloc(macro_nnode*macro_ndofn*nne*ndim,sizeof(double));
    int *row = PGFEM_calloc(macro_nnode*macro_ndofn*nne*ndim,sizeof(int));
    int *col = PGFEM_calloc(macro_nnode*macro_ndofn*nne*ndim,sizeof(int));
    int idx = 0;

    /* get list */
    for(int a=0; a<macro_nnode; a++){
      for(int b=0; b<macro_ndofn; b++){
  	for(int w=0; w<nne; w++){
  	  for(int g=0; g<ndim; g++){
  	    /* row idx is the global dof id */
  	    row[idx] = global_dof_ids[idx_2_gen(w,g,nne,ndim)] - 1;

  	    /* skip if boundary condition */
  	    if(row[idx] < 0) continue;
	    else {
	      /* column idx is simply the index of the macro dof */
	      col[idx] = idx_2_gen(a,b,macro_nnode,macro_ndofn);

	      /* get values */
	      val_01[idx] = *(K_01_e + idx_K_gen(a,b,w,g,macro_nnode,
						 macro_ndofn,nne,ndim));
	      val_10[idx] = *(K_10_e + idx_K_gen(w,g,a,b,nne,ndim,
						 macro_nnode,macro_ndofn));
	      /* increment counter */
	      idx++;
	    }
	  }
	}
      }
    }

    /* add list to matrices */
    err += PGFEM_par_matrix_add_to_values(idx,row,col,val_01,K_01);
    err += PGFEM_par_matrix_add_to_values(idx,row,col,val_10,K_10);

    /* free memory */
    free(val_01);
    free(val_10);
    free(row);
    free(col);
  }

  /* deallocate function-scope information */
  free(local_dof_ids);
  free(global_dof_ids);
  free(K_01_e);
  free(K_10_e);
  free(K_00_e);
  free(traction_res_e);
  free(traction_e);
  free(disp);
  free(x);
  free(y);
  free(z);

  return err;
}

/* Compute and assemble the contribution to the macroscale
    tangent. Involves multiple linear solves of the microscale */
static int compute_ms_cohe_job_tangent(const int macro_ndof,
				       const long *macro_dof_ids,
				       const PGFEM_par_matrix *K_01,
				       const PGFEM_par_matrix *K_10,
				       const PGFem3D_opt *o,
				       const int print_level,
				       COMMON_MICROSCALE *c,
				       double *K_00_contrib)
{
  int err = 0;
  int myrank = 0;
  int nproc = 0;
  MPI_Comm_rank(c->mpi_comm,&myrank);
  MPI_Comm_size(c->mpi_comm,&nproc);

  /* allocate stuff for solver loops */
  const int n_loc_rows = c->DomDof[myrank];
  double *loc_rhs = PGFEM_calloc(n_loc_rows,sizeof(double));
  double *loc_sol = PGFEM_calloc(n_loc_rows,sizeof(double));
  SOLVER_INFO *info = PGFEM_calloc(1,sizeof(SOLVER_INFO));

  /* KK = K_01 K_11^{-1} K_10 */
  double *KK = PGFEM_calloc(macro_ndof*macro_ndof,sizeof(double));
  int total_it = 0;
  double max_res = 0.0;

  /* compute the (++) block */
  const int block_ndof = macro_ndof/2;
  double *KK_block = PGFEM_calloc(block_ndof*block_ndof,sizeof(double));
  for(int i=0; i<block_ndof; i++){
    /* null the local vectors */
    nulld(loc_rhs,n_loc_rows);
    nulld(loc_sol,n_loc_rows);

    /* extract the i-th column from K_10 */
    err += PGFEM_par_matrix_get_column(K_10,i,loc_rhs);

    /* compute loc_sol = K_11^{-1} K_10{i} */
    if(i == 0){
      /* only setup solver for first time through */
      solve_system(o,loc_rhs,loc_sol,1,1,c->DomDof,info,
		   c->SOLVER,c->mpi_comm);
    } else {
      solve_system_no_setup(o,loc_rhs,loc_sol,1,1,c->DomDof,info,
			    c->SOLVER,c->mpi_comm);
    }

    /* check solver error status and print solve information */
    err += solve_system_check_error(PGFEM_stderr,*info);
    total_it += info->n_iter;
    max_res = (max_res > info->res_norm) ? max_res : (info->res_norm);

    for(int j=0; j<block_ndof; j++){
      /* extract row (column) from K_01 */
      err += PGFEM_par_matrix_get_column(K_01,j,loc_rhs);

      /* compute KK(j,i) = K_01(j,:) {K_11^-1 K_10}(:,i). NOTE: The
	 dot product is gathered to all procs, thus KK is identical on
	 all processes. NO ASSEMBLY OF KK IS REQUIRED! */
      const int idx = idx_2_gen(j,i,block_ndof,block_ndof);
      err += PGFEM_par_vec_dot(n_loc_rows,loc_rhs,loc_sol,
			       c->mpi_comm,KK_block + idx);
    }
  }

  /* copy (++) block to others */
  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      const int sign = ((i==j)? 1:-1);
      for(int k=0; k<block_ndof; k++){
  	for(int m=0; m<block_ndof; m++){
  	  const int idx = idx_K(i,k,j,m,2,block_ndof);
  	  *(KK + idx) = sign*KK_block[idx_2_gen(k,m,block_ndof,block_ndof)];
  	}
      }
    }
  }

  /* compute final K_00_contrib */
  /* NOTE: need to scale for K_11 */
  double scale = c->supports->lc/c->VVolume;
  cblas_daxpy(macro_ndof*macro_ndof,-1.0/scale,KK,1,K_00_contrib,1);
  
  /* if(myrank == 0){ */
  /*   PGFEM_printf("KK_block:\n"); */
  /*   print_array_d(PGFEM_stdout,KK_block,block_ndof*block_ndof, */
  /* 		  block_ndof,block_ndof); */

  /*   PGFEM_printf("K_00:\n"); */
  /*   print_array_d(PGFEM_stdout,K_00_contrib,macro_ndof*macro_ndof, */
  /* 		  macro_ndof,macro_ndof); */

  /*   PGFEM_printf("KK:\n"); */
  /*   print_array_d(PGFEM_stdout,KK,macro_ndof*macro_ndof, */
  /* 		  macro_ndof,macro_ndof); */
  /* } */

  /* print summary */
  if(myrank == 0){
    PGFEM_printf("=== TANGENT CONTRIBUTION SUMMARY ===\n");
    PGFEM_printf("Macro Ndof: %d || Total iter: %d || Max. Res.: %11.5e\n",
		 macro_ndof,total_it,max_res);
  }

  free(loc_rhs);
  free(loc_sol);
  free(info);
  free(KK);
  free(KK_block);
  return err;
}

static int print_ms_cohe_job(const MS_COHE_JOB_INFO *job,
			     const COMMON_MICROSCALE *c,
			     const MICROSCALE_SOLUTION *s,
			     const PGFem3D_opt *o)
{
  int err = 0;
  /* exit early if not marked for output */
  if(!job->print_flag) return err;
  if(o->vis_format == VIS_NONE) return err;

  int myrank = 0;
  int nproc = 0;
  err += MPI_Comm_rank(c->mpi_comm,&myrank);
  err += MPI_Comm_size(c->mpi_comm,&nproc);

  /* compute Mises stresses and strains */
  Mises(c->ne,s->sig_e,s->eps,o->analysis_type);

  /* compute the output filename for the job */
  int len = snprintf(NULL,0,"%s_p%.5d_cel%.5d",
		     o->ofname,job->proc_id,job->elem_id) + 1;
  char *ofname = PGFEM_calloc(len,sizeof(char));
  sprintf(ofname,"%s_p%.5d_cel%.5d",
	  o->ofname,job->proc_id,job->elem_id);

  switch(o->vis_format){
  case VIS_NONE: /* do nothing */ break;
  case VIS_VTK:
    if(myrank == 0){
      VTK_print_master(o->opath,ofname,s->p_tim,nproc,o);
    }
    VTK_print_vtu(o->opath,ofname,s->p_tim,myrank,c->ne,
		  c->nn,c->node,c->elem,c->supports,
		  s->r,s->sig_e,s->eps,o);

    if(o->cohesive){
      if(myrank == 0){
	VTK_print_cohesive_master(o->opath,ofname,s->p_tim,nproc,o);
      }
      VTK_print_cohesive_vtu(o->opath,ofname,s->p_tim,myrank,c->nce,c->node,
			     c->coel,c->supports,s->r,c->ensight,o);
    }
    break;
  default:
    if(myrank == 0){
      PGFEM_printerr("WARNING: output format not supported at microscale!\n"
		     "No output generated!\n");
    }
    break;
  }

  free(ofname);
  return err;
}
