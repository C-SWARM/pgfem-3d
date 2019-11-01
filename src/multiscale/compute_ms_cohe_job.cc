#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "Newton_Raphson.h"
#include "PGFEM_par_matvec.h"
#include "allocation.h"
#include "compute_ms_cohe_job.h"
#include "displacement_based_element.h"
#include "enumerations.h"
#include "get_ndof_on_elem.h"
#include "get_dof_ids_on_elem.h"
#include "incl.h"
#include "index_macros.h"
#include "interface_macro.h"
#include "matice.h"
#include "pgf_fe2_restart.h"
#include "stiffmat_fd.h"
#include "utils.h"
#include "vtk_output.h"
#include <mkl_cblas.h>
#include <cassert>
#include <cmath>
#include <cstring>

using namespace pgfem3d;
using namespace multiscale::net;

#ifndef JOB_LOGGING
#define JOB_LOGGING 1
#endif

static constexpr int NDIM = 3;

/*==== STATIC FUNCTION PROTOTYPES ====*/

/** Wrapper for Newton Raphson. */
static int ms_cohe_job_nr(MultiscaleCommon *c,
                          MULTISCALE_SOLUTION *s,
                          const PGFem3D_opt *opts,
                          int *n_step,
                          const int mp_id,int tag_TP);

/** Set the job supports appropriately from the jump at (n) and
    (n+1). Also set the normal to the interface. */
static int set_job_supports(const MS_COHE_JOB_INFO *p_job,
                            SUPP sup);


/** compute and initialize the dense matrix structure for the mixed
    tangents */
static int initialize_ms_cohe_job_mixed_tangents
(const int n_cols,
 MultiscaleCommon *common,
 const int mp_id);

/** Call stiffmat_fd using information for the microscale job. The
    microscale tangent reset and assembled. Contains (blocking)
    collective communication */
static int ms_cohe_job_compute_micro_tangent(MultiscaleCommon *c,
                                             MULTISCALE_SOLUTION *s,
                                             PGFem3D_opt *o,
                                             const int mp_id);

/** Compute all of the microscale terms for the macroscale tangent and
    residual. The mixed tangents are assembled and the cohesive
    tangent at the macroscale is computed and distributed to all
    processors. Note that only DISP analysis type is currently
    supported. */
static int compute_ms_cohe_job_micro_terms(const MS_COHE_JOB_INFO *job,
                                           const MultiscaleCommon *c,
                                           const MULTISCALE_SOLUTION *s,
                                           const PGFem3D_opt *o,
                                           PGFEM_par_matrix *K_01,
                                           PGFEM_par_matrix *K_10,
                                           double *K_00_contrib,
                                           const int mp_id);

/** Compute the contibutions to the microscale terms from a single
    element */
static int compute_elem_micro_terms(const int elem_id,
                                    const MS_COHE_JOB_INFO *job,
                                    const MultiscaleCommon *c,
                                    const MULTISCALE_SOLUTION *s,
                                    const PGFem3D_opt *o,
                                    PGFEM_par_matrix *K_01,
                                    PGFEM_par_matrix *K_10,
                                    double *K_00_contrib,
                                    const int mp_id);

/** Compute and assemble the contribution to the macroscale
    tangent. Involves multiple linear solves of the microscale */
static int compute_ms_cohe_job_tangent(const int macro_ndof,
                                       const long *macro_dof_ids,
                                       const PGFEM_par_matrix *K_01,
                                       const PGFEM_par_matrix *K_10,
                                       const PGFem3D_opt *o,
                                       const int print_level,
                                       MultiscaleCommon *c,
                                       double *K_00_contrib);

static int print_ms_cohe_job(const MS_COHE_JOB_INFO *job,
                             const MultiscaleCommon *c,
                             const MULTISCALE_SOLUTION *s,
                             const PGFem3D_opt *o,
                             const int mp_id);

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
  const double tn = cblas_dnrm2(NDIM,job->traction_n,1);
  const double tnp1 = cblas_dnrm2(NDIM,job->traction,1);
  const double jnp1 = cblas_dnrm2(NDIM,job->jump,1);

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
  memcpy(job->jump_n,job->jump,NDIM*sizeof(double));
  memcpy(job->traction_n,job->traction,NDIM*sizeof(double));

  /* RESET: traction, residual and tangent */
  memset(job->traction,0,NDIM*sizeof(double));
  memset(job->traction_res,0,job->ndofe*sizeof(double));
  memset(job->K_00_contrib,0,job->ndofe*job->ndofe*sizeof(double));

  return cell_failure_detected;
}

//Frobenius of the deflection matrix
int find_norm(MultiscaleCommon *c, MULTISCALE_SOLUTION *s,int myrank) {
  int i,j;
  double sum = 0;
  int high_norm = 0;
  for (i = 0; i < 3 ; i++) {
    for(j = 0; j < 3; j++) {
      if(i==j) {
        sum += (1 + c->supports->defl[i + 3*j])*(1 + c->supports->defl[i + 3*j]);
      } else {
        sum += c->supports->defl[i + 3*j]*c->supports->defl[i + 3*j];
      }
    }
  }

  sum = sqrt(sum);

  if (sum > 2.01 || sum < 1.65 ) {
    high_norm = 1;
  }

return high_norm;
}


/*==== API FUNCTION DEFINITIONS ====*/

/** Given a microscale job, compute the solution on the time
    increment, store the traction and assemble the tangent */
int compute_ms_cohe_job(const int job_id,
                        MS_COHE_JOB_INFO *p_job,
                        Microscale *microscale,
                        const int mp_id)
{
  int err = 0;
  const int print_level = 1;
  MultiscaleCommon *common = microscale;
  MULTISCALE_SOLUTION *sol = microscale->sol + job_id;

  int myrank = 0;
  microscale->net->comm_rank(common->comm,&myrank);
  if(myrank == 0){
    PGFEM_printf("=== Microscale cell %d of %d ===\n",
                 job_id+1,microscale->idx_map.size-1);
  }

  /* switch set up job */
  switch(p_job->job_type){
   case JOB_UPDATE: case JOB_PRINT:/* do nothing */ break;

   default: /* reset state to macro time (n) */
    /* reset the microscale solution to macro time (n) */
    err += reset_MULTISCALE_SOLUTION(sol,microscale);

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

    /* consider file_.msm for PDE or Taylor */
    if(microscale->opts->custom_micro) 
    {
      int tag_TP = microscale->opts->methods[p_job->elem_id + common->nce*0];
      PGFEM_printf("\n ---> simulation_method=%ld", tag_TP);
      
      if (tag_TP)
      {
        PGFEM_printf("\n CUSTOM: ms_cohe_job_nr for PDE called!!! \n");
        err += ms_cohe_job_nr(common,sol,microscale->opts,&(p_job->n_step), mp_id, tag_TP);
      } 
      else {
        PGFEM_printf("\n CUSTOM: ms_cohe_job_nr for Taylor called!!! \n");
        err += ms_cohe_job_nr(common,sol,microscale->opts,&(p_job->n_step), mp_id, tag_TP);
      }
    }  
    //else {
      //PGFEM_printf("\n No CUSTOM: NR default for PDE called!!! \n");
      //err += ms_cohe_job_nr(common,sol,microscale->opts,&(p_job->n_step), mp_id, 1); 
    //}

    /*=== INTENTIONAL DROP THROUGH ===*/
   case JOB_NO_COMPUTE_EQUILIBRIUM:
    if(sol->failed) break;
    /* Do not compute the microscale equilibrium, just assemble the
       microscale tangent */
    {

      if(JOB_LOGGING && myrank == 0) PGFEM_printf("=== Microscale TANGENT ===\n");
      /*  compute the equilibriated microscale tangent */
      err += ms_cohe_job_compute_micro_tangent(common,sol,microscale->opts,mp_id);
      if(myrank == 0) PGFEM_printf("\n");

      static int init_mixed = 0;
      if(!init_mixed){
        /* initialize matrices K_01 and K_10 */
        /* allocate for quad macro elements */
        const int n_cols = 8*NDIM;
        /* const int n_cols = p_job->ndofe; */
        if(JOB_LOGGING && myrank == 0) PGFEM_printf("=== INIT MIXED TANGENT ===\n");
        err += initialize_ms_cohe_job_mixed_tangents(n_cols,common,mp_id);
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
                                             K_01,K_10,K_00_contrib,mp_id);

      if(JOB_LOGGING && myrank == 0) PGFEM_printf("=== MACROSCALE TANGENT/SCHUR ===\n");
      /* compute macroscale tangent contribution */
      err += compute_ms_cohe_job_tangent(n_cols,p_job->g_dof_ids,K_01,K_10,
                                         microscale->opts,print_level,
                                         common,K_00_contrib);
    }
    break;

   case JOB_UPDATE:
    /* update job information and set cell failure condition. */
     {
       int failure_detected = update_job_information(p_job);
       if (!sol->failed && failure_detected) sol->failed = 1;
     }

     /* update the solution and state variables n <- n+1 */
     err += update_MULTISCALE_SOLUTION(sol,microscale);

     /* Set the supports for the new n-state and null the increment */
     err += set_job_supports(p_job,common->supports);

     break;

   case JOB_PRINT:
    /* output the job based on the print flag to the file specified by
       the options and solution step id */
    err += print_ms_cohe_job(p_job,common,sol,microscale->opts,mp_id);

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
                             int micro_rank,
			     int macro_rank,
                             double *loc_res,
                             const int mp_id)
{
  int err = 0;
  /* exit early without doing anything if not the owning domain */
  if(macro_rank != p_job->proc_id) return err;

  for(int i=0; i<p_job->nnode*NDIM; i++){
    const int idx = p_job->loc_dof_ids[i] - 1;
    if(idx >= 0) loc_res[idx] += p_job->traction_res[i];
  }

  return err;
}

/*=====================================*/
/*==== STATIC FUNCTION DEFINITIONS ====*/
/*=====================================*/

static int ms_cohe_job_nr(MultiscaleCommon *c,
                          MULTISCALE_SOLUTION *s,
                          const PGFem3D_opt *opts,
                          int *n_step,const int mp_id,int tag_TP)
{
  int err = 0;
  double pores = 0.0;
  const int print_level = 0;
  *n_step = 0;

  /* copy of load increment */
  double *sup_defl = PGFEM_calloc(double, c->supports->npd);
  memcpy(sup_defl,c->supports->defl_d,c->supports->npd*sizeof(double));

  Newton_Raphson_multiscale(print_level,c,s,NULL,NULL,opts,sup_defl,&pores,n_step,tag_TP);

  free(sup_defl);
  return err;
}/* ms_cohe_job_nr() */

static int set_job_supports(const MS_COHE_JOB_INFO *p_job,
                            SUPP sup)
{
  /* NOTE: This is written for TOTAL LAGRANGIAN */
  int err = 0;
  if(!sup->multi_scale || sup->npd != 6) err++;
  memcpy(sup->N0,p_job->normal,NDIM*sizeof(double));
  for(int i=0; i<NDIM; i++){
    sup->defl[i] = 0.5*p_job->jump_n[i];
    sup->defl[i+NDIM] = -0.5*p_job->jump_n[i];
    sup->defl_d[i] = 0.5*(p_job->jump[i] - p_job->jump_n[i]);
    sup->defl_d[i+NDIM] = -0.5*(p_job->jump[i] - p_job->jump_n[i]);
  }

  /* compute the macroscopic deformation gradient. This is
     particularly important for if microscale equilibrium is not
     computed. */
  err += compute_macro_grad_u(sup->F0,sup,DISP);
  return err;
}/* set_job_supports() */

static int initialize_ms_cohe_job_mixed_tangents
(const int n_cols,
 MultiscaleCommon *common,
 const int mp_id)
{
  int err = 0;
  int myrank = common->rank;
  int nproc = common->nproc;
  if(!(common->K_01 == NULL && common->K_10 == NULL)) return err;

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
  int *rows = PGFEM_calloc(int, n_entries);
  int *cols = PGFEM_calloc(int, n_entries);

  /* get the row numbers (G dof ids). If the dof is a BC, set to
     -1. This makes the duplicate counting work properly in the
     initialization */
  {
    int idx = 0;
    for(int i=0; i<common->nn; i++){
      for(int j=0; j<common->ndofn; j++){
        rows[idx] = common->node[i].id_map[mp_id].Gid[j] - 1;
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
                                     common,
                                     (PGFEM_par_matrix**) &common->K_01);

  err += initialize_PGFEM_par_matrix(n_rows,n_cols,n_own_rows,
                                     n_entries,rows,cols,
                                     common,
                                     (PGFEM_par_matrix**) &common->K_10);
  free(rows);
  free(cols);
  return err;
}

static int ms_cohe_job_compute_micro_tangent(MultiscaleCommon *c,
                                             MULTISCALE_SOLUTION *s,
                                             PGFem3D_opt *o,
                                             const int mp_id)
{
  int err = 0;
  int myrank = c->rank;
  int nproc = c->nproc;
  double nor_min = c->lin_err;

  /* reset the microscale tangent to zeros */
  c->SOLVER->zero();

  err += stiffmat_fd_multiscale(c,s,o,0,nor_min,0,myrank,nproc);

  /* finalize the microscale tangent matrix assembly */
  c->SOLVER->assemble();

  return err;
}

static int compute_ms_cohe_job_micro_terms(const MS_COHE_JOB_INFO *job,
                                           const MultiscaleCommon *c,
                                           const MULTISCALE_SOLUTION *s,
                                           const PGFem3D_opt *o,
                                           PGFEM_par_matrix *K_01,
                                           PGFEM_par_matrix *K_10,
                                           double *K_00_contrib,
                                           const int mp_id)
{
  int err = 0;

  /* error check analysis type */
  switch(o->analysis_type){
   case DISP: break;
   case CM: if (o->cm == DISP) break;
    /* deliberate drop through */
   default:
    PGFEM_printerr("ERROR: analysis type %d is not supported!\n"
                   "%s:%s:%d",o->analysis_type,__func__,__FILE__,__LINE__);
    PGFEM_Comm_code_abort(c, 0);
  }

  /* zero the matrices/vectors */
  const int macro_ndof = job->ndofe;
  err += PGFEM_par_matrix_zero_values(K_01);
  err += PGFEM_par_matrix_zero_values(K_10);
  nulld(K_00_contrib,macro_ndof*macro_ndof);
  nulld(job->traction_res,macro_ndof);
  nulld(job->traction,NDIM);


  /* volume elements */
  for(int i=0; i<c->ne; i++){
    err += compute_elem_micro_terms(i,job,c,s,o,K_01,
                                    K_10,K_00_contrib,mp_id);
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
  c->net->allreduce(NET_IN_PLACE,K_00_contrib,
		    job->ndofe*job->ndofe,
		    NET_DT_DOUBLE,NET_OP_SUM,c->comm);

  /* assemble traction residual from all domains */
  c->net->allreduce(NET_IN_PLACE,job->traction_res,job->ndofe,
		    NET_DT_DOUBLE,NET_OP_SUM,c->comm);

  /* assemble traction from all domains */
  c->net->allreduce(NET_IN_PLACE,job->traction,NDIM,
		    NET_DT_DOUBLE,NET_OP_SUM,c->comm);
  
  /* finish assembly of the distributed matrices */
  err += PGFEM_par_matrix_end_assembly(K_01,K_01_comm);
  err += PGFEM_par_matrix_end_assembly(K_10,K_10_comm);

  return err;
}

static int compute_elem_micro_terms(const int elem_id,
                                    const MS_COHE_JOB_INFO *job,
                                    const MultiscaleCommon *c,
                                    const MULTISCALE_SOLUTION *s,
                                    const PGFem3D_opt *o,
                                    PGFEM_par_matrix *K_01,
                                    PGFEM_par_matrix *K_10,
                                    double *K_00_contrib,
                                    const int mp_id)
{
  int err = 0;
  int myrank = c->rank;

  /* important macro quantities */
  const int macro_nnode = job->nnode;
  const int macro_ndofn = NDIM;
  const int macro_ndof = job->ndofe;
  const double *macro_shape_func = job->shape;
  const double *macro_normal = job->normal;
  const double macro_int_wt = job->int_wt;
  const double layer_thickness = c->supports->lc;

  /* microscale element information */
  const Element *elem = c->elem + elem_id;
  const int ndofn = c->ndofn;
  const int nne = elem->toe;
  const long *node_ids = elem->nod;
  const int ndofe = get_ndof_on_elem_nodes(nne,node_ids,c->node,ndofn);

  /* allocate microscale information */
  long *local_dof_ids = PGFEM_calloc(long, ndofe);
  long *global_dof_ids = PGFEM_calloc(long, ndofe);
  double *K_01_e = PGFEM_calloc(double, macro_ndof*ndofe);
  double *K_10_e = PGFEM_calloc(double, ndofe*macro_ndof);
  double *K_00_e = PGFEM_calloc(double, macro_ndof*macro_ndof);
  double *traction_res_e = PGFEM_calloc(double, macro_ndof);
  double *traction_e = PGFEM_calloc(double, NDIM);
  double *disp = PGFEM_calloc(double, ndofe);
  double *x = PGFEM_calloc(double, nne);
  double *y = PGFEM_calloc(double, nne);
  double *z = PGFEM_calloc(double, nne);


  /* microscale dof ids */
  get_dof_ids_on_elem_nodes(0,nne,ndofn,node_ids,c->node,local_dof_ids ,mp_id);
  get_dof_ids_on_elem_nodes(1,nne,ndofn,node_ids,c->node,global_dof_ids,mp_id);

  /* microscale node coordinates and displacements */
  switch(o->analysis_type){
   case CM:
    if (o->cm != DISP) {
      PGFEM_printerr("ERROR: multiscale cohesive modeling not "
                     "yet implemented for this analysis type (%d)!\n"
                     "%s:%s:%d\n",o->analysis_type,__func__,__FILE__,__LINE__);
      PGFEM_Abort();
      break;
    }
    /* deliberate drop through */
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
   case CM:
    if (o->cm != DISP) {
      PGFEM_printerr("ERROR: multiscale cohesive modeling not "
                     "yet implemented for this analysis type (%d)!\n"
                     "%s:%s:%d\n",o->analysis_type,__func__,__FILE__,__LINE__);
      PGFEM_Abort();
      break;
    }
    /* deliberate drop through */
   case DISP:
    err += DISP_cohe_micro_terms_el(K_00_e,K_01_e,K_10_e,traction_res_e,
                                    traction_e,
                                    macro_nnode,macro_ndofn,macro_int_wt,
                                    macro_shape_func,macro_normal,
                                    layer_thickness,c->VVolume,
                                    elem_id,ndofn,nne,x,y,z,c->elem,
                                    c->hommat,node_ids,c->node,
                                    s->eps,s->sig_e,c->supports,disp,s->dt);
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
  cblas_daxpy(NDIM,1.0,traction_e,1,job->traction,1);

  {
    /* allocate enough space for full matrix */
    double *val_01 = PGFEM_calloc(double, macro_nnode*macro_ndofn*nne*NDIM);
    double *val_10 = PGFEM_calloc(double, macro_nnode*macro_ndofn*nne*NDIM);
    int *row = PGFEM_calloc(int, macro_nnode*macro_ndofn*nne*NDIM);
    int *col = PGFEM_calloc(int, macro_nnode*macro_ndofn*nne*NDIM);
    int idx = 0;

    /* get list */
    for(int a=0; a<macro_nnode; a++){
      for(int b=0; b<macro_ndofn; b++){
        for(int w=0; w<nne; w++){
          for(int g=0; g<NDIM; g++){
            /* row idx is the global dof id */
            row[idx] = global_dof_ids[idx_2_gen(w,g,nne,NDIM)] - 1;

            /* skip if boundary condition */
            if(row[idx] < 0) continue;
            else {
              /* column idx is simply the index of the macro dof */
              col[idx] = idx_2_gen(a,b,macro_nnode,macro_ndofn);

              /* get values */
              val_01[idx] = *(K_01_e + idx_K_gen(a,b,w,g,macro_nnode,
                                                 macro_ndofn,nne,NDIM));
              val_10[idx] = *(K_10_e + idx_K_gen(w,g,a,b,nne,NDIM,
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
                                       MultiscaleCommon *c,
                                       double *K_00_contrib)
{
  int err = 0;
  int myrank = c->rank;

  /* allocate stuff for solver loops */
  const int n_loc_rows = c->DomDof[myrank];
  double *loc_rhs = PGFEM_calloc(double, n_loc_rows);
  double *loc_sol = PGFEM_calloc(double, n_loc_rows);
  SOLVER_INFO *info = PGFEM_calloc(SOLVER_INFO, 1);

  /* KK = K_01 K_11^{-1} K_10 */
  double *KK = PGFEM_calloc(double, macro_ndof*macro_ndof);
  int total_it = 0;
  double max_res = 0.0;

  /* compute the (++) block */
  const int block_ndof = macro_ndof/2;
  double *KK_block = PGFEM_calloc(double, block_ndof*block_ndof);
  for(int i=0; i<block_ndof; i++){
    /* null the local vectors */
    nulld(loc_rhs,n_loc_rows);
    nulld(loc_sol,n_loc_rows);

    /* extract the i-th column from K_10 */
    err += PGFEM_par_matrix_get_column(K_10,i,loc_rhs);

    /* compute loc_sol = K_11^{-1} K_10{i} */
    if(i == 0){
      /* only setup solver for first time through */
      c->SOLVER->solveSystem(o, loc_rhs, loc_sol, 1, 1, c->DomDof, info);
    } else {
      c->SOLVER->solveSystemNoSetup(o, loc_rhs, loc_sol, 1, 1, c->DomDof, info);
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
                               c,KK_block + idx);
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
  /*          block_ndof,block_ndof); */

  /*   PGFEM_printf("K_00:\n"); */
  /*   print_array_d(PGFEM_stdout,K_00_contrib,macro_ndof*macro_ndof, */
  /*          macro_ndof,macro_ndof); */

  /*   PGFEM_printf("KK:\n"); */
  /*   print_array_d(PGFEM_stdout,KK,macro_ndof*macro_ndof, */
  /*          macro_ndof,macro_ndof); */
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
                             const MultiscaleCommon *c,
                             const MULTISCALE_SOLUTION *s,
                             const PGFem3D_opt *o,const int mp_id)
{
  int err = 0;
  /* exit early if not marked for output */
  if(!job->print_flag) return err;
  if(o->vis_format == VIS_NONE) return err;

  int myrank = c->rank;
  int nproc = c->nproc;

  /* compute Mises stresses and strains */
  Mises(c->ne,s->sig_e,s->eps,o->analysis_type);

  /* compute the output filename for the job */
  int len = snprintf(NULL,0,"%s_p%.5d_cel%.5d",
                     o->ofname,job->proc_id,job->elem_id) + 1;
  char *ofname = PGFEM_calloc(char, len);
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
                  s->r,NULL,NULL,s->sig_e,s->eps,o,mp_id);

    if(o->cohesive){
      if(myrank == 0){
        VTK_print_cohesive_master(o->opath,ofname,s->p_tim,nproc,o);
      }
      VTK_print_cohesive_vtu(o->opath,ofname,s->p_tim,myrank,c->nce,c->node,
                             c->coel,c->supports,s->r,c->ensight,o,mp_id);
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
