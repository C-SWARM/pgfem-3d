/* HEADER */
#include "microscale_information.h"

#include "allocation.h"
#include "enumerations.h"
#include "gen_path.h"
#include "mesh_load.h"
#include "material.h"
#include "read_input_file.h"
#include "interface_macro.h"
#include "utils.h"
#include "homogen.h"
#include "set_fini_def.h"
#include "generate_dof_ids.h"
#include "Psparse_ApAi.h"
#include "initialize_damage.h"
#include "in.h"
#include "PGFEM_par_matvec.h"

static const int ndim = 3;

/*==== STATIC FUNCTION PROTOTYPES ====*/
static void initialize_COMMON_MICROSCALE(COMMON_MICROSCALE *common);
static void build_COMMON_MICROSCALE(const PGFem3D_opt *opts,
				    MPI_Comm mpi_comm,
				    COMMON_MICROSCALE *common);
static void destroy_COMMON_MICROSCALE(COMMON_MICROSCALE *common);

static void initialize_MICROSCALE_SOLUTION(MICROSCALE_SOLUTION *sol);
static void build_MICROSCALE_SOLUTION(MICROSCALE_SOLUTION *sol,
				      const COMMON_MICROSCALE *common,
				      const int analysis);
static void destroy_MICROSCALE_SOLUTION(MICROSCALE_SOLUTION *sol,
					const COMMON_MICROSCALE *common,
					const int analysis);

static void build_MICROSCALE_SOLUTION_BUFFERS(void *buffer,
					      const int local_len,
					      const int global_len);

static void destroy_MICROSCALE_SOLUTION_BUFFERS(void *buffer);

/**
 * \brief Private data type for storing common solution buffers.
 *
 * This portion is allocated and held by COMMON_MICROSCALE and
 * MICROSCALE_SOLUTION simply holds pointers to the various buffers.
 */
typedef struct MICROSCALE_SOLUTION_BUFFERS{
    /* local vectors */
    double *f;
    double *d_r;
    double *rr;
    double *D_R;
    double *R;
    double *RR;
    double *f_u;
    double *f_defl;
    double *RRn;
    double *U;
    double *DK;
    double *dR;

    /* global vectors */
    double *BS_f;
    double *BS_x;
    double *BS_RR;
    double *BS_f_u;
    double *BS_d_r;
    double *BS_D_R;
    double *BS_rr;
    double *BS_R;
    double *BS_U;
    double *BS_DK;
    double *BS_dR;
} MICROSCALE_SOLUTION_BUFFERS;


/*==== API FUNCTIONS ====*/
void initialize_MICROSCALE(MICROSCALE **microscale)
{
  *microscale = PGFEM_calloc(1,sizeof(MICROSCALE));
  (*microscale)->opts = PGFEM_calloc(1,sizeof(PGFem3D_opt));
  (*microscale)->common = PGFEM_calloc(1,sizeof(COMMON_MICROSCALE));
  (*microscale)->n_solutions = 0;
  (*microscale)->sol = NULL;

  set_default_options((*microscale)->opts);
  initialize_COMMON_MICROSCALE((*microscale)->common);
}/* initialize_MICROSCALE */

void build_MICROSCALE(MICROSCALE *microscale,
		      MPI_Comm mpi_comm,
		      const int argc,
		      char **argv)
{
  int myrank = 0;
  int nproc = 0;
  MPI_Comm_rank(mpi_comm,&myrank);
  MPI_Comm_size(mpi_comm,&nproc);

  /* parse the command-line-style options */
  parse_command_line(argc,argv,myrank,microscale->opts);

  /* error check the command line */
  if (microscale->opts->solverpackage != HYPRE){
    if(myrank == 0)
      PGFEM_printerr("ERROR: only the HYPRE solver are supported!"
	"%s:%s:%d\n",__func__,__FILE__,__LINE__);
  PGFEM_Comm_code_abort(mpi_comm,0);
}

  if (microscale->opts->analysis_type != DISP){
    if(myrank == 0)
      PGFEM_printerr("ERROR: only DISP analysis is supported!"
	      "%s:%s:%d\n",__func__,__FILE__,__LINE__);
    PGFEM_Comm_code_abort(mpi_comm,0);
  }

  /* attempt to write to output directory */
  if(make_path(microscale->opts->opath,DIR_MODE) != 0){
    if(myrank == 0){
      PGFEM_printerr("ERROR: could not write to %s!\n",
	      microscale->opts->opath);
    }
    PGFEM_Comm_code_abort(mpi_comm,0);
  }

  /*=== BUILD COMMON ===*/
  build_COMMON_MICROSCALE(microscale->opts,mpi_comm,microscale->common);
  microscale->common->supports->multi_scale = microscale->opts->multi_scale;

}/* build_MICROSCALE */

void build_MICROSCALE_solutions(const int n_solutions,
				MICROSCALE *microscale)
{
  /*=== BUILD SOLUTIONS ===*/
  microscale->n_solutions = n_solutions;
  microscale->sol = PGFEM_calloc(n_solutions,sizeof(MICROSCALE_SOLUTION));
  for(int i=0; i<n_solutions; i++){
    initialize_MICROSCALE_SOLUTION(microscale->sol + i);
    build_MICROSCALE_SOLUTION(microscale->sol +i,microscale->common,
			      microscale->opts->analysis_type);
  }
}

void destroy_MICROSCALE(MICROSCALE *microscale)
{
  if(microscale != NULL){
    for(int i=0; i<microscale->n_solutions; i++){
      destroy_MICROSCALE_SOLUTION(microscale->sol+i,
				  microscale->common,
				  microscale->opts->analysis_type);
    }
    free(microscale->sol);
    destroy_COMMON_MICROSCALE(microscale->common);
    free(microscale->common);
    free(microscale->opts);
  }
  free(microscale);
} /* destroy_MICROSCALE */

int reset_MICROSCALE(MICROSCALE *m)
{
  int err = 0;
  int myrank = 0;
  COMMON_MICROSCALE *c = m->common;
  err += MPI_Comm_rank(c->mpi_comm,&myrank);

  const long loc_ndof = c->ndofd;
  const long g_ndof = c->DomDof[myrank];

  for(int i=0; i<m->n_solutions; i++){
    MICROSCALE_SOLUTION *sol = m->sol+i;
    err += reset_MICROSCALE_SOLUTION(sol,loc_ndof,g_ndof);
  }

  return err;
}/* reset_MICROSCALE_SOLUTION */

int update_MICROSCALE(MICROSCALE *m)
{
  int err = 0;
  int myrank = 0;
  COMMON_MICROSCALE *c = m->common;
  err += MPI_Comm_rank(c->mpi_comm,&myrank);

  const long loc_ndof = c->ndofd;
  const long g_ndof = c->DomDof[myrank];

  for(int i=0; i<m->n_solutions; i++){
    MICROSCALE_SOLUTION *sol = m->sol+i;
    err += update_MICROSCALE_SOLUTION(sol,loc_ndof,g_ndof);
  }

  return err;
} /* update_MICROSCALE_SOLUTION */

int reset_MICROSCALE_SOLUTION(MICROSCALE_SOLUTION *sol,
			      const int loc_ndof,
			      const int g_ndof)
{
  int err = 0;
  /* reset displacement (solution) vector */
  memcpy(sol->r,sol->rn,loc_ndof*sizeof(double));

  /* null all of the other local vectors */
  nulld(sol->f,loc_ndof);
  nulld(sol->d_r,loc_ndof);
  nulld(sol->rr,loc_ndof);
  nulld(sol->D_R,loc_ndof);
  nulld(sol->R,loc_ndof);
  nulld(sol->RR,loc_ndof);
  nulld(sol->f_u,loc_ndof);
  nulld(sol->f_defl,loc_ndof);
  nulld(sol->RRn,loc_ndof);
  nulld(sol->U,loc_ndof);
  nulld(sol->DK,loc_ndof);
  nulld(sol->dR,loc_ndof);

  /* null all of the "global" vectors */
  nulld(sol->BS_f,g_ndof);
  nulld(sol->BS_x,g_ndof);
  nulld(sol->BS_RR,g_ndof);
  nulld(sol->BS_f_u,g_ndof);
  nulld(sol->BS_d_r,g_ndof);
  nulld(sol->BS_D_R,g_ndof);
  nulld(sol->BS_rr,g_ndof);
  nulld(sol->BS_R,g_ndof);
  nulld(sol->BS_U,g_ndof);
  nulld(sol->BS_DK,g_ndof);
  nulld(sol->BS_dR,g_ndof);

  /* reset state variables */
  /*** NOT YET IMPLEMENTED ***/

  return err;
}/* reset_MICROSCALE_SOLUTION */

int update_MICROSCALE_SOLUTION(MICROSCALE_SOLUTION *sol,
			       const int loc_ndof,
			       const int g_ndof)
{
  int err = 0;
  /* copy r -> rn  */
  memcpy(sol->rn,sol->r,loc_ndof*sizeof(double));

  /* leave other solution vectors alone */

  /* update state variables */
  /*** NOT IMPLEMENTED YET ***/

  return err;
}/* update_MICROSCALE_SOLUTION */

/*==== STATIC FUNCTION DEFINITIONS ===*/

static void initialize_COMMON_MICROSCALE(COMMON_MICROSCALE *common)
{
  /* options /solver information */
  common->SOLVER = NULL;
  common->Ap = NULL;
  common->Ai = NULL;

  /* communication information */
  common->pgfem_comm = NULL;
  common->mpi_comm = MPI_COMM_WORLD;
  common->nbndel = 0;
  common->bndel = NULL;
  common->ndofd = 0;
  common->DomDof = NULL;
  common->GDof = 0;

  /* mesh info */
  common->nn = 0;
  common->ne = 0;
  common->nce = 0;
  common->ndofn = 0;
  common->npres = 0;
  common->VVolume = 0.0;
  common->node = NULL;
  common->elem = NULL;
  common->coel = NULL;
  common->n_orient = 0;
  common->matgeom = NULL;
  common->nhommat = 0;
  common->hommat = NULL;
  common->ensight = NULL;
  common->supports = NULL;
  common->n_co_props = 0;
  common->co_props = NULL;

  /* mixed tangents */
  common->K_01 = NULL;
  common->K_10 = NULL;

  /* solution buffers */
  common->solution_buffer = NULL;
}

static void build_COMMON_MICROSCALE(const PGFem3D_opt *opts,
				    MPI_Comm mpi_comm,
				    COMMON_MICROSCALE *common)
{
  int myrank = 0;
  int nproc = 0;
  MPI_Comm_rank(mpi_comm,&myrank);
  MPI_Comm_size(mpi_comm,&nproc);
  /* initialize the solver information */
  initialize_PGFEM_HYPRE_solve_info(&common->SOLVER);
  common->SOLVER->solver_type = opts->solver;
  common->SOLVER->precond_type = opts->precond;
  common->mpi_comm = mpi_comm;
  common->maxit_nl = 5;

  switch(opts->vis_format){
  case VIS_ENSIGHT: case VIS_VTK:
    common->ensight = PGFEM_calloc (1,sizeof(ENSIGHT_1));
    break;
  default: break;
  }

  /*=== READ MICROSCALE INPUT FILES ===*/
  char *in_fname = PGFEM_calloc(500,sizeof(char));
  sprintf(in_fname,"%s/%s",opts->ipath,opts->ifname);

  long Gnn = 0;
  /* read *.in input files */
  {
    int err = 0;

    /* supports that should not exist for microscale */
    long nln = 0;
    long nle_s = 0;
    long nle_v = 0;
    ZATNODE *znod = NULL;
    ZATELEM *zele_s = NULL;
    ZATELEM *zele_v = NULL;

    /* other variables that will be ignored/un-saved */
    long ni = 0;
    long nmat = 0;
    long n_con = 0;
    MATERIAL *mater = NULL;

    err = read_input_file(opts,mpi_comm,&common->nn, &Gnn,&common->ndofn,
			  &common->ne, &ni,&common->lin_err,
			  &common->lim_zero,&nmat,&n_con,
			  &common->n_orient,&common->node,
			  &common->elem,&mater,&common->matgeom,
			  &common->supports,&nln,&znod,&nle_s,&zele_s,
			  &nle_v,&zele_v);

    /* error reading file(s) */
    if(err){
      PGFEM_printerr("[%d]ERROR: incorrectly formatted input file!\n",
	      myrank);
      PGFEM_Comm_code_abort(mpi_comm,0);
    }

    /* error in input */
    if(nln > 0 || nle_s > 0 || nle_v > 0){
      PGFEM_printerr("[%d]ERROR: applied loads are not consistent w/"
	      " multiscale analysis!\n",myrank);
      PGFEM_Comm_code_abort(mpi_comm,0);
    }

    destroy_zatnode(znod,nln);
    destroy_zatelem(zele_s,nle_s);
    destroy_zatelem(zele_v,nle_v);

    /* override prescribed displacements */
    if(opts->override_pre_disp){
      if(override_prescribed_displacements(common->supports,opts) != 0){
	PGFEM_printerr("[%d]ERROR: an error was encountered when"
		       " reading the displacement override file.\n"
		       "Be sure that there are enough prescribed"
		       " displacements in the file.\n",myrank);
	PGFEM_Abort();
      }
    }

    /* read microscale normal/thickness */
    if(opts->multi_scale){
      common->supports->multi_scale = opts->multi_scale;
      err = read_interface_macro_normal_lc(opts->ipath,common->supports);
      if(err != 0){
	PGFEM_printerr("[%d] ERROR: could not read normal from file!\n"
		       "Check that the file \"%s/normal.in\""
		       " exists and try again.\n",
		       myrank,opts->ipath);
	PGFEM_Comm_code_abort(mpi_comm,0);
      }
    }

    /* homogenized material properties */
    Mat_3D_orthotropic (nmat,mater,opts->analysis_type);
    long ***a = aloc3l (nmat,nmat,n_con);
    common->nhommat = list(a,common->ne,nmat,n_con,common->elem);
    common->hommat= build_hommat( common->nhommat);
    hom_matrices(a,common->ne,nmat,n_con,common->elem,mater,common->matgeom,
		 common->hommat,common->matgeom->SH,opts->analysis_type);

    dealoc3l(a,nmat,nmat);
    free(mater);
  }/* end reading *.in */

  /* read cohesive stuff */
  if (opts->cohesive){
    char *filename = PGFEM_calloc(500,sizeof(char));
    sprintf(filename,"%s%d.in.co_props",in_fname,myrank);
    FILE *in1 = PGFEM_fopen(filename,"r");
    read_cohesive_properties(in1,&common->n_co_props,
			     &common->co_props,mpi_comm);
    PGFEM_fclose(in1);
      
    /* read coheisve elements */
    sprintf (filename,"%s%d.in.co",in_fname,myrank);
    in1 = PGFEM_fopen(filename,"r");

    /* temporary leftovers from old file format */ 
    long ncom = 0;     
    fscanf (in1,"%ld\n",&ncom); 
    double **comat = aloc2 (ncom,4);
      
    /* read the cohesive element info */
    common->coel = read_cohe_elem (in1,ncom,ndim,common->nn,common->node,
				   &common->nce,comat,common->ensight,
				   opts->vis_format,myrank,
				   common->co_props);
    dealoc2 (comat,ncom);
    PGFEM_fclose (in1);
    free(filename);
  }

  /*=== END FILE READING ===*/
  list_el_prescribed_def(common->supports,common->node,common->elem,
			 NULL,common->ne,0,common->nn);

  common->bndel = list_boundary_el(common->ne,common->elem,common->nn,
				   common->node,myrank,&common->nbndel);

  common->DomDof = PGFEM_calloc(nproc,sizeof(long));

  common->ndofd = generate_local_dof_ids(common->ne,common->nce,common->nn,
					 common->ndofn,common->node,
					 common->elem,common->coel,NULL,
					 mpi_comm);

  common->DomDof[myrank] =
    generate_global_dof_ids(common->ne,common->nce,common->nn,
			    common->ndofn,common->node,common->elem,
			    common->coel,NULL,mpi_comm);

  MPI_Allgather(MPI_IN_PLACE,1,MPI_LONG,common->DomDof,
		1,MPI_LONG,mpi_comm);

  {
    long Gndof = 0;
    for(int i=0; i<nproc; i++)  Gndof += common->DomDof[i];
    set_HYPRE_row_col_bounds(common->SOLVER,Gndof,common->DomDof,myrank);
  }

  renumber_global_dof_ids(common->ne,common->nce,0,common->nn,
			  common->ndofn,common->DomDof,common->node,
			  common->elem,common->coel,NULL,mpi_comm);

  long NBN = distribute_global_dof_ids(common->ne,common->nce,
				       0,common->nn,
				       common->ndofn,ndim,
				       common->node,common->elem,
				       common->coel,NULL,mpi_comm);

  /* global stiffness pattern and communication structure */
  common->Ap = PGFEM_calloc(common->DomDof[myrank]+1,sizeof(int));
  common->pgfem_comm = PGFEM_calloc (1,sizeof(COMMUN_1));
  initialize_commun(common->pgfem_comm);
  common->Ai = Psparse_ApAi(nproc,myrank,common->ne,0,common->nn,
			    common->ndofn,common->ndofd,common->elem,
			    NULL,common->node,common->Ap,common->nce,
			    common->coel,common->DomDof,&common->GDof,
			    common->pgfem_comm,mpi_comm,opts->cohesive);
  pgfem_comm_build_fast_maps(common->pgfem_comm,common->ndofd,
			     common->DomDof[myrank],common->GDof);

  hypre_initialize(common->Ap,common->Ai,common->DomDof[myrank],
		   opts->maxit,common->lin_err,common->SOLVER,opts,
		   mpi_comm);

  if(!common->supports->multi_scale){
  common->VVolume = T_VOLUME (common->ne,ndim,common->elem,common->node);
  MPI_Allreduce(MPI_IN_PLACE,&common->VVolume,1,MPI_DOUBLE,
		MPI_SUM,mpi_comm);
  } else {
    common->VVolume = common->supports->v0;
  }

  /* allocate solution_buffers */
  common->solution_buffer = PGFEM_calloc(1,sizeof(MICROSCALE_SOLUTION_BUFFERS));
  build_MICROSCALE_SOLUTION_BUFFERS(common->solution_buffer,
				    common->ndofd,common->DomDof[myrank]);


  free(in_fname);

  /* compute/print summary information */
  long mesh_info[7];
  mesh_info[0] = Gnn;
  mesh_info[1] = NBN;
  mesh_info[2] = common->ne;
  mesh_info[3] = common->nce;
  mesh_info[4] = common->nbndel;
  mesh_info[5] = common->DomDof[myrank];
  mesh_info[6] = common->Ap[common->DomDof[myrank]];
  MPI_Allreduce(MPI_IN_PLACE,mesh_info+2,5,MPI_LONG,MPI_SUM,mpi_comm);

  if(myrank == 0){
    print_interpreted_options(opts);
    if(opts->multi_scale){
      if(common->supports->npd >= 9){
	PGFEM_printf("\n*** BULK Multiscale Modelling ***\n");
      } else {
	PGFEM_printf("\n*** INTERFACE Multiscale Modelling ***\n");
      }
    }

    PGFEM_printf ("\n");
    PGFEM_printf ("Number of total nodes                    : %ld\n",mesh_info[0]);
    PGFEM_printf ("Number of nodes on domain interfaces     : %ld\n",mesh_info[1]);
    PGFEM_printf ("Total number of elements                 : %ld\n",mesh_info[2]);
    if (opts->cohesive == 1){
      PGFEM_printf ("Number of cohesive elements              : %ld\n",mesh_info[3]);
    }
    PGFEM_printf ("Number of elems on the COMM interfaces   : %ld\n",mesh_info[4]);
    PGFEM_printf ("Total number of degrees of freedom       : %ld\n",mesh_info[5]);
    PGFEM_printf ("Total number of nonzeros in the matrix   : %ld\n",mesh_info[6]);
    PGFEM_printf ("\n");
    PGFEM_printf ("Volume: %f\n\n",common->VVolume);
  }
}

static void destroy_COMMON_MICROSCALE(COMMON_MICROSCALE *common)
{
  int nproc = 0;
  MPI_Comm_size(common->mpi_comm,&nproc);
  destroy_PGFEM_HYPRE_solve_info(common->SOLVER);
  free(common->Ap);
  free(common->Ai);
  destroy_commun(common->pgfem_comm,nproc);
  free(common->bndel);
  free(common->DomDof);
  destroy_node(common->nn,common->node);
  destroy_elem(common->elem,common->ne);
  destroy_coel(common->coel,common->nce);
  destroy_matgeom(common->matgeom,common->n_orient);
  destroy_hommat(common->hommat,common->nhommat);
  destroy_ensight(common->ensight);
  destroy_supp(common->supports);
  destroy_cohesive_props(common->n_co_props,common->co_props);
  destroy_PGFEM_par_matrix((PGFEM_par_matrix *) common->K_01);
  destroy_PGFEM_par_matrix((PGFEM_par_matrix *) common->K_10);
  destroy_MICROSCALE_SOLUTION_BUFFERS(common->solution_buffer);
  free(common->solution_buffer);
}

static void initialize_MICROSCALE_SOLUTION(MICROSCALE_SOLUTION *sol)
{
  /* stress/strain/state information */
  sol->sig_e = NULL;
  sol->sig_n = NULL;
  sol->eps = NULL;
  sol->crpl = NULL;
  sol->npres = 0;
  sol->elem_state_info = NULL;
  sol->coel_state_info = NULL;

  /* solution information */
  /* local vectors */
  sol->r = NULL;
  sol->f = NULL;
  sol->d_r = NULL;
  sol->rr = NULL;
  sol->D_R = NULL;
  sol->R = NULL;
  sol->f_defl = NULL;
  sol->RR = NULL;
  sol->f_u = NULL;
  sol->RRn = NULL;
  sol->U = NULL;
  sol->DK = NULL;
  sol->dR = NULL;

  /* global vectors */
  sol->BS_f = NULL;
  sol->BS_f_u = NULL;
  sol->BS_x = NULL;
  sol->BS_RR = NULL;
  sol->BS_d_r = NULL;
  sol->BS_D_R = NULL;
  sol->BS_rr = NULL;
  sol->BS_R = NULL;
  sol->BS_U = NULL;
  sol->BS_DK = NULL;
  sol->BS_dR = NULL;

  /* convergence info */
  sol->dt = 0.0;
  sol->times = NULL;
  sol->tim = 0;
  sol->p_tim = 0;
  sol->NORM = 0.0;
}

static void build_MICROSCALE_SOLUTION(MICROSCALE_SOLUTION *sol,
				      const COMMON_MICROSCALE *common,
				      const int analysis)
{
  int myrank = 0;
  int nproc = 0;
  MPI_Comm_rank(common->mpi_comm,&myrank);
  MPI_Comm_size(common->mpi_comm,&nproc);

  const long local_len = common->ndofd;
  const long global_len = common->DomDof[myrank];
  const size_t len_double = sizeof(double);

  sol->sig_e = build_sig_il(common->ne,analysis,common->elem);
  sol->eps = build_eps_il(common->ne,common->elem,analysis);
  initialize_damage(common->ne,common->elem,common->hommat,
		    sol->eps,analysis);

  /* need to figure out elem/coel_state_info indexing */

  switch(analysis){
  case STABILIZED: case MINI: case MINI_3F: sol->npres = 4; break;
  case DISP: sol->npres = 0; break;
  default: sol->npres = 1; break;
  }
  build_pressure_nodes(common->ne,sol->npres,common->elem,
		       sol->sig_e,sol->eps,analysis);
  set_fini_def (common->ne,sol->npres,common->elem,
		sol->eps,sol->sig_e,analysis);

  /*=== crystal plasticity is not currently supported ===*/
  /* if (analysis == FS_CRPL) { */
  /*   sol->crpl = PGFEM_calloc (common->nhommat,sizeof(CRPL)); */
  /*   read_cryst_plast (in1,common->nhommat,sol->crpl,plc); */
  /*   build_crystal_plast (common->ne,common->elem,sol->sig_e,sol->eps, */
  /* 			 sol->crpl,analysis,plc); */
  /*   set_fini_def_pl(common->ne,sol->npres,common->elem, */
  /* 		    sol->eps,sol->sig_e,sol->crpl,analysis,plc); */
  /* } */

  /* local solution vectors */
  sol->rn = PGFEM_calloc(local_len,len_double);
  sol->r = PGFEM_calloc(local_len,len_double);

  /* sol->f = PGFEM_calloc(local_len,len_double); */
  /* sol->d_r = PGFEM_calloc(local_len,len_double); */
  /* sol->rr = PGFEM_calloc(local_len,len_double); */
  /* sol->D_R = PGFEM_calloc(local_len,len_double); */
  /* sol->R = PGFEM_calloc(local_len,len_double); */
  /* sol->f_defl = PGFEM_calloc(local_len,len_double); */
  /* sol->RR = PGFEM_calloc(local_len,len_double); */
  /* sol->f_u = PGFEM_calloc(local_len,len_double); */
  /* sol->RRn = PGFEM_calloc(local_len,len_double); */
  /* sol->U = PGFEM_calloc(local_len,len_double); */
  /* sol->DK = PGFEM_calloc(local_len,len_double); */
  /* sol->dR = PGFEM_calloc(local_len,len_double); */

  /* sol->BS_f = PGFEM_calloc(global_len,len_double); */
  /* sol->BS_f_u = PGFEM_calloc(global_len,len_double); */
  /* sol->BS_x = PGFEM_calloc(global_len,len_double); */
  /* sol->BS_RR = PGFEM_calloc(global_len,len_double); */
  /* sol->BS_d_r = PGFEM_calloc(global_len,len_double); */
  /* sol->BS_D_R = PGFEM_calloc(global_len,len_double); */
  /* sol->BS_rr = PGFEM_calloc(global_len,len_double); */
  /* sol->BS_R = PGFEM_calloc(global_len,len_double); */
  /* sol->BS_U = PGFEM_calloc(global_len,len_double); */
  /* sol->BS_DK = PGFEM_calloc(global_len,len_double); */
  /* sol->BS_dR = PGFEM_calloc(global_len,len_double); */

  /* Get pointers to the shared solution workspace */
  {
    MICROSCALE_SOLUTION_BUFFERS *buff =
      (MICROSCALE_SOLUTION_BUFFERS *) common->solution_buffer;

    sol->f	= buff->f     ;
    sol->d_r	= buff->d_r   ;
    sol->rr	= buff->rr    ;
    sol->D_R	= buff->D_R   ;
    sol->R	= buff->R     ;
    sol->f_defl = buff->f_defl;
    sol->RR	= buff->RR    ;
    sol->f_u	= buff->f_u   ;
    sol->RRn	= buff->RRn   ;
    sol->U	= buff->U     ;
    sol->DK	= buff->DK    ;
    sol->dR	= buff->dR    ;

    sol->BS_f	= buff->BS_f  ;
    sol->BS_f_u = buff->BS_f_u;
    sol->BS_x	= buff->BS_x  ;
    sol->BS_RR	= buff->BS_RR ;
    sol->BS_d_r = buff->BS_d_r;
    sol->BS_D_R = buff->BS_D_R;
    sol->BS_rr	= buff->BS_rr ;
    sol->BS_R	= buff->BS_R  ;
    sol->BS_U	= buff->BS_U  ;
    sol->BS_DK	= buff->BS_DK ;
    sol->BS_dR	= buff->BS_dR ;
  }

  sol->dt = 0.0;
  sol->times = PGFEM_calloc(3,len_double);
  sol->tim = 0;
  sol->NORM = 0.0;
}

static void destroy_MICROSCALE_SOLUTION(MICROSCALE_SOLUTION *sol,
					const COMMON_MICROSCALE *common,
					const int analysis)
{
  free(sol->r);
  free(sol->rn);
  /* free(sol->f); */
  /* free(sol->d_r); */
  /* free(sol->rr); */
  /* free(sol->D_R); */
  /* free(sol->R); */
  /* free(sol->f_defl); */
  /* free(sol->RR); */
  /* free(sol->f_u); */
  /* free(sol->RRn); */
  /* free(sol->U); */
  /* free(sol->DK); */
  /* free(sol->dR); */

  /* free(sol->BS_f); */
  /* free(sol->BS_f_u); */
  /* free(sol->BS_x); */
  /* free(sol->BS_RR); */
  /* free(sol->BS_d_r); */
  /* free(sol->BS_D_R); */
  /* free(sol->BS_rr); */
  /* free(sol->BS_R); */
  /* free(sol->BS_U); */
  /* free(sol->BS_DK); */
  /* free(sol->BS_dR); */
  free(sol->times);

  destroy_sig_il(sol->sig_e,common->elem,common->ne,analysis);
  /* destroy sig_n */
  destroy_eps_il(sol->eps,common->elem,common->ne,analysis);
  //destroy_crpl
  free(sol->elem_state_info);
  free(sol->coel_state_info);
}

static void build_MICROSCALE_SOLUTION_BUFFERS(void *buffer,
					      const int local_len,
					      const int global_len)
{
  MICROSCALE_SOLUTION_BUFFERS *buff = (MICROSCALE_SOLUTION_BUFFERS*) buffer;
  const unsigned long len_double = sizeof(double);

 /* local solution vectors */
  buff->f = PGFEM_calloc(local_len,len_double);
  buff->d_r = PGFEM_calloc(local_len,len_double);
  buff->rr = PGFEM_calloc(local_len,len_double);
  buff->D_R = PGFEM_calloc(local_len,len_double);
  buff->R = PGFEM_calloc(local_len,len_double);
  buff->f_defl = PGFEM_calloc(local_len,len_double);
  buff->RR = PGFEM_calloc(local_len,len_double);
  buff->f_u = PGFEM_calloc(local_len,len_double);
  buff->RRn = PGFEM_calloc(local_len,len_double);
  buff->U = PGFEM_calloc(local_len,len_double);
  buff->DK = PGFEM_calloc(local_len,len_double);
  buff->dR = PGFEM_calloc(local_len,len_double);

  /* global solution vectors */
  buff->BS_f = PGFEM_calloc(global_len,len_double);
  buff->BS_f_u = PGFEM_calloc(global_len,len_double);
  buff->BS_x = PGFEM_calloc(global_len,len_double);
  buff->BS_RR = PGFEM_calloc(global_len,len_double);
  buff->BS_d_r = PGFEM_calloc(global_len,len_double);
  buff->BS_D_R = PGFEM_calloc(global_len,len_double);
  buff->BS_rr = PGFEM_calloc(global_len,len_double);
  buff->BS_R = PGFEM_calloc(global_len,len_double);
  buff->BS_U = PGFEM_calloc(global_len,len_double);
  buff->BS_DK = PGFEM_calloc(global_len,len_double);
  buff->BS_dR = PGFEM_calloc(global_len,len_double);
}

static void destroy_MICROSCALE_SOLUTION_BUFFERS(void *buffer)
{
  MICROSCALE_SOLUTION_BUFFERS *buff = (MICROSCALE_SOLUTION_BUFFERS*) buffer;
  free(buff->f);
  free(buff->d_r);
  free(buff->rr);
  free(buff->D_R);
  free(buff->R);
  free(buff->f_defl);
  free(buff->RR);
  free(buff->f_u);
  free(buff->RRn);
  free(buff->U);
  free(buff->DK);
  free(buff->dR);

  free(buff->BS_f);
  free(buff->BS_f_u);
  free(buff->BS_x);
  free(buff->BS_RR);
  free(buff->BS_d_r);
  free(buff->BS_D_R);
  free(buff->BS_rr);
  free(buff->BS_R);
  free(buff->BS_U);
  free(buff->BS_DK);
  free(buff->BS_dR);
}
