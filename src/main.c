/*
 * In this implmentation, initial displacement, velocity, and material density are spcicified in file_base_*.initial files.
 * NON zero density triggers adding inertia. If density is equal to zero, no inertia and transient terms are evaluated.
 */
 /*
 Restart is added which is dependent on read_VTK_file(char fn[], double *r) which can be found lib/VTK_IO/src (written in C++)
 For the restart you have to set #define SAVE_RESTART_FILE 1.
 When you restart the simulation, change -1 to your time step number in the file_base_0.inital
 Nov. 17/2014, Sangmin Lee
 */

/* HEADER */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "PFEM3d.h"
#ifndef ENUMERATIONS_H
#include "enumerations.h"
#endif

/* Standard headers/libs */
#include <time.h>
#include <stdlib.h>
#include <sys/time.h> 
#include <sys/resource.h>
#include <assert.h>

/* Extra libs */
#include "renumbering.h"

/*=== PFEM3d headers ===*/
#include "PGFEM_io.h"
#include "allocation.h"
#include "Arc_length.h"
#include "build_distribution.h"
#include "homogen.h"
#include "hypre_global.h"
#include "in.h"
#include "load.h"
#include "matice.h"
#include "matrix_printing.h"
#include "Newton_Raphson.h"
#include "out.h"
#include "Printing.h"
#include "print_dist.h"
#include "profiler.h"
#include "Psparse_ApAi.h"
#include "read_cryst_plast.h"
#include "renumber_ID.h"
#include "RNPsparse_ApAi.h"
#include "set_fini_def.h"
#include "skyline.h"
#include "utils.h"
#include "SetGlobalNodeNumbers.h"
#include "interface_macro.h"
#include "computeMacroF.h"
#include "computeMacroS.h"
#include "vtk_output.h"
#include "PGFem3D_options.h"
#include "gen_path.h"
#include "initialize_damage.h"
#include "bounding_element.h"
#include "bounding_element_utils.h"
#include "element.h"
#include "node.h"
#include "generate_dof_ids.h"
#include "applied_traction.h"
#include "read_input_file.h"

#include "three_field_element.h"

static const int periodic = 0;
static const int ndim = 3;

/*****************************************************/
/*           BEGIN OF THE COMPUTER CODE              */
/*****************************************************/
#define SAVE_RESTART_FILE 1

#ifndef NO_VTK_LIB
#include "PGFem3D_to_VTK.hpp"


int read_initial_from_VTK(const PGFem3D_opt *opts, int myrank, int *restart, double *u0, double *u1)
{
  int err = 0;
  char filename[1024];
  
  sprintf(filename,"%s/restart/VTK/STEP_%.5d/%s_%d_%d.vtu",opts->opath,*restart,opts->ofname,myrank, *restart);   
  err += read_VTK_file(filename, u0);      
  sprintf(filename,"%s/VTK/STEP_%.5d/%s_%d_%d.vtu",opts->opath,*restart,opts->ofname,myrank, *restart);   
  err += read_VTK_file(filename, u1);
    
  return err;
}      

#else
int read_initial_from_VTK(const PGFem3D_opt *opts, int myrank, int *restart, double *u0, double *u1s)
{
  if(myrank==0)
  {
    PGFEM_printerr("Restart with VTK is not supported!\n");
    PGFEM_printerr("Enforce to turn off restart!\n");
  }
  
  *restart = -1;  
  return 0;
}
#endif

int write_restart_disp(double *u0, double *u1, const PGFem3D_opt *opts, int myrank, int nodeno, int nsd, int stepno)
{
  
  char restart_path[1024];
  sprintf(restart_path, "%s/restart/STEP_%.5d", opts->opath,stepno);
                
  if(make_path(restart_path,DIR_MODE) != 0)
  {
    PGFEM_printf("Directory (%s) not created!\n",restart_path);
    abort();                   
  }  
 
  char filename[1024];
  sprintf(filename,"%s/restart/STEP_%.5d/%s_%d_%d.res",opts->opath,stepno,opts->ofname,myrank, stepno);   
  FILE *fp = fopen(filename,"w");

  if(fp == NULL)
  {    
    printf("Fail to open file [%s]. finishing\n", filename);
    exit(0);  
  }
  
  for(int a=0; a<nodeno; a++)
  {
    for(int b=0; b<nsd; b++)
      fprintf(fp, "%e %e ", u0[a*nsd + b], u1[a*nsd + b]);
    
    fprintf(fp, "\n");    
  }
  
  fclose(fp);
  return 0;
}


int read_restart_disp(double *u0, double *u1, const PGFem3D_opt *opts, int myrank, int nodeno, int nsd, int stepno)
{
 
  char filename[1024];
  sprintf(filename,"%s/restart/STEP_%.5d/%s_%d_%d.res",opts->opath,stepno,opts->ofname,myrank, stepno);   
  FILE *fp = fopen(filename,"r");

  if(fp == NULL)
  {    
    printf("Fail to open file [%s]. finishing\n", filename);      
    exit(0);  
  }
  
  for(int a=0; a<nodeno; a++)
  {
    for(int b=0; b<nsd; b++)
      fscanf(fp, "%lf %lf", u0+a*nsd + b, u1+a*nsd + b);    
  }
  
  fclose(fp);
  return 0;
}

double read_initial_values(double *u0, double *u1, double *rho, const PGFem3D_opt *opts, int myrank, int nodeno, int nmat, double dt, int *restart)
{
  char filename[1024];
  char line[1024];
  double alpha = 0.5;
  for(int a=0; a<nmat; a++)
    rho[a] = 0.0;
    
  sprintf(filename,"%s/%s%d.initial",opts->ipath,opts->ifname,0);
  FILE *fp_0 = fopen(filename,"r");

  *restart = -1;
  if(fp_0 != NULL)
  {  
    while(fgets(line, 1024, fp_0)!=NULL) 
	  {
	    if(line[0]=='#')
	      continue;
        
	    sscanf(line, "%d", restart);
	    break;
	  }
    fclose(fp_0);
  }
  
  if(*restart>0)
  { 
    int nsd = 3;
    if(opts->analysis_type==DISP)
      read_restart_disp(u0, u1, opts, myrank, nodeno, nsd, *restart);
    else
      read_initial_from_VTK(opts, myrank, restart, u0, u1);
  } 
  
  sprintf(filename,"%s/%s%d.initial",opts->ipath,opts->ifname,myrank);
  FILE *fp = fopen(filename,"r");

       
  if(fp == NULL)
  {    
    if(myrank==0)
      printf("Fail to open file [%s]. Quasi steady state\n", filename);
    return alpha;
  }

  if(myrank==0)
  {
    while(fgets(line, 1024, fp)!=NULL) 
	  {
	    if(line[0]=='#')
	      continue;
      
	    double temp;  
	    sscanf(line, "%lf", &temp);
	    break;
	  }
  }
  
  while(fgets(line, 1024, fp)!=NULL) 
  {
    if(line[0]=='#')
	    continue;
        
    sscanf(line, "%lf", &alpha);
    break;
  }
     

  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
	    continue;
    for(int a=0; a<nmat; a++)
  	{    
	    sscanf(line, "%lf", rho+a);
	    if(a<nmat-1)
	      fgets(line, 1024, fp);      
	  }
    break;
  } 
  if(*restart>0)
  {
    fclose(fp);
      return alpha;  
  }
                   
  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
	    continue;
        
    long nid;
    double u[3], v[3];        
    sscanf(line, "%ld %lf %lf %lf %lf %lf %lf", &nid, u+0, u+1, u+2, v+0, v+1, v+2);

    u1[nid*3+0] = u[0];
    u1[nid*3+1] = u[1];
    u1[nid*3+2] = u[2];    
    u0[nid*3+0] = u[0]-dt*v[0];
    u0[nid*3+1] = u[1]-dt*v[1];
    u0[nid*3+2] = u[2]-dt*v[2];
  }              
    
  fclose(fp);
  return alpha;       
}

int single_scale_main(int argc,char *argv[])
{
  /* Create MPI communicator. Currently aliased to MPI_COMM_WORLD but
     may change */
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  struct rusage usage;
  long nn = 0;
  long ndofn = 0;
  long ne = 0;
  long nmat = 0;
  long nln = 0;
  long nc = 0;
  long np = 0;
  long i = 0;
  long ndofd = 0;
  long nle_s = 0;
  long nle_v = 0;
  long ni = 0;
  long ***a = NULL;
  long nhommat = 0;
  long iter_max = 0;
  long FNR = 0;
  long nt = 0;
  long tim = 0;
  long AT = 0;
  double *r = NULL;
  double *f = NULL;
  double err = 0.0;
  double limit = 0.0;
  double nor_min = 0.0;
  double *d_r = NULL;
  double *rr = NULL;
  double *D_R = NULL;
  double *R = NULL;
  double *f_defl = NULL;
  double *f_u = NULL;
  double nor1 = 0.0;
  double DET = 0.0;
  double dlm = 0.0;
  double lm = 0.0;
  double dlm0 = 0.0;
  double DLM = 0.0;
  double dAL = 0.0;
  double *DK = NULL;
  double NORM = 0.0;
  double VVolume = 0.0;

  FILE *in1 = NULL;
  FILE *out = NULL;
  NODE *node = NULL;

  ELEMENT *elem = NULL;
  MATERIAL *mater = NULL;
  MATGEOM matgeom = NULL;
  HOMMAT *hommat = NULL;
  SIG *sig_e = NULL;
  SIG *sig_n = NULL;
  EPS *eps = NULL;
  ZATNODE *znod = NULL;
  ZATELEM *zele_s = NULL;
  ZATELEM *zele_v = NULL;
  SUPP sup = NULL;
  char filename[500],in_dat[500],out_dat[500];

  /* boundary elements */
  int n_be = 0;
  BOUNDING_ELEMENT *b_elems = NULL;

  long gem = 0;
  long temp_int = 0;
  long sky = 0;
  
  /* CRYSTAL PLASTICITY */
  CRPL *crpl = NULL;
  double *RR = NULL;
  double *RRn = NULL;
  double *sup_defl = NULL; 
  double *times = NULL;
  double dt = 0.0;
  double gama = 0.0;
  double GNOR = 0.0;
  long npres = 0;
  long n_p = 0;
  long *print = NULL;
  long *tim_load = NULL;
  long nlod_tim = 0;

  /* ARC LENGTH METHOD */
  double dAL0 = 0.0;
  double dALMAX = 0.0;

  /*
    ARC == 0 :: Crisfield
    ARC == 1 :: Simo ONLY ARC LENG METHOD FOR PARALELL
  */
  
  long ARC = 1;
  
  /* MI */
  double GVolume = 0.0;
  int namelen = 0;
  int *Ap = NULL;
  int *Ai = NULL;
  int GDof = 0;
  int APP = 0;
  long *DomDof = NULL;
  long Gndof=0;
  long *DomNe = NULL;
  long Gne=0;
  long *DomNn = NULL;
  long Gnn = 0;
  long NBN = 0;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  COMMUN comm;
  
  /* COHESIVE ELEMENTS */
  long nce = 0;
  long ncom = 0;
  long ITT =0;
  long Gnce = 0;
  double **comat = NULL;
  double dt0 = 0.0;
  double *U = NULL;
  double *dR = NULL;
  double pores = 0.0;
  COEL *coel = NULL;
  int n_co_props = 0;
  cohesive_props *co_props = NULL;

  double  *BS_x = NULL;
  double *BS_f = NULL;
  double *BS_RR = NULL;
  double *BS_d_r = NULL;
  double *BS_D_R = NULL;
  double *BS_rr = NULL;
  double *BS_R = NULL;
  double *BS_f_u = NULL;
  double *BS_U = NULL;
  double *BS_DK = NULL;
  double *BS_dR = NULL;
  
  /* HYPRE specific */
  PGFEM_HYPRE_solve_info *PGFEM_hypre = NULL;
  
  /* Ensight */
  ENSIGHT ensight;

  int *dist = NULL;

  /* original volume */
  double oVolume = 0.0;

  /* macro F S */
  double *GF = NULL;
  double *GS = NULL;
  double *GP = NULL;

  double hypre_time = 0.0;

  /* ***** Set up debug log ***** */
  FILE *debug_log = NULL;
  /* debug_log = fopen("debug.log","w"); */
  debug_log = stdout;
  /* debug_log = stderr; */

  int myrank = 0;
  int nproc = 0;

  /*=== END INITIALIZATION === */

  int flag_MPI_Init;
  MPI_Initialized(&flag_MPI_Init);
  if(!flag_MPI_Init)
  {
    MPI_Init (&argc,&argv);  
  }

  MPI_Comm_rank (mpi_comm,&myrank);
  MPI_Comm_size (mpi_comm,&nproc);
  MPI_Get_processor_name (processor_name,&namelen);
  PGFEM_initialize_io(NULL,NULL);

  if(myrank == 0){
    PGFEM_printf("=== SINGLE SCALE ANALYSIS ===\n\n");
  }

  double total_time = 0.0;
  /* MPI_Barrier(mpi_comm); */
  total_time -= MPI_Wtime();

  MPI_Errhandler_set(mpi_comm,MPI_ERRORS_ARE_FATAL);

  if(myrank == 0) {
    PGFEM_printf("\n\nInitializing PFEM3d\n\n");
  }

 
  /*=== Parse the command line for options ===*/
  PGFem3D_opt options;
  if (argc <= 2){
    if(myrank == 0){
      print_usage(stdout);
    }
    exit(0);
  }
  set_default_options(&options);
  re_parse_command_line(myrank,2,argc,argv,&options);
  if(myrank == 0){
    print_options(stdout,&options);
  }

  /* set up solver variables */
  if(options.solverpackage == HYPRE){
    initialize_PGFEM_HYPRE_solve_info(&PGFEM_hypre);
    PGFEM_hypre->solver_type = options.solver;
    PGFEM_hypre->precond_type = options.precond;
  } else {
    if(myrank == 0){
      PGFEM_printerr("ERROR: Only HYPRE solvers are supported.\n");
    }
    PGFEM_Comm_code_abort(mpi_comm,0);
  }
  
  /* for attaching debugger */
  while (options.debug);

  /* visualization */
  switch(options.vis_format){
  case VIS_ENSIGHT: case VIS_VTK:
    ensight = PGFEM_calloc (1,sizeof(ENSIGHT_1));
    break;
  default: ensight = NULL; break;
  }

  /* abort early if unrecognized analysis type */
  if(options.analysis_type < 0
     || options.analysis_type >= ANALYSIS_MAX){
    if(myrank == 0){
      PGFEM_printerr("ERROR: Unregognized analysis type given (%d)!"
	      " Please provide an analysis type (see the help menu).\n",
	      options.analysis_type);
    }
    PGFEM_Abort();
  }

  if(options.restart){
    if(myrank == 0){
      PGFEM_printerr("WARNING: Restart is not currently supported!\n"
	      "Starting the simulation from the beginning with "
	      "given options.\n");
    }
  }

  sprintf(in_dat,"%s/%s",options.ipath,options.ifname);
  sprintf(out_dat,"%s/%s",options.opath,options.ofname);

  if(make_path(options.opath,DIR_MODE) != 0){
    if(myrank == 0){
      PGFEM_printf("Could not create path (%s)!\n"
	     "Please check input and try again.\n\n",options.opath);
      print_usage(stdout);
    }
    PGFEM_Comm_code_abort(mpi_comm,-1);
  }

  /*=== READ *.in INPUT FILES ===*/
  {
    int in_err = 0;
    in_err = read_input_file(&options,mpi_comm,&nn,&Gnn,&ndofn,
			     &ne,&ni,&err,&limit,&nmat,&nc,&np,&node,
			     &elem,&mater,&matgeom,&sup,&nln,&znod,
			     &nle_s,&zele_s,&nle_v,&zele_v);
    if(in_err){
      PGFEM_printerr("[%d]ERROR: incorrectly formatted input file!\n",
	      myrank);
      PGFEM_Abort();
    }

    /* TEMP reset ndofn for testing. In many/most cases ndofn should
       actually be ndim */
    /* ndofn = 3; */
  }

  /*=== OVERRIDE PRESCRIBED DISPLACEMENTS ===*/
  if(options.override_pre_disp){
    if(override_prescribed_displacements(sup,&options) != 0){
      PGFEM_printerr("[%d]ERROR: an error was encountered when"
	      " reading the displacement override file.\n"
	      "Be sure that there are enough prescribed"
	      " displacements in the file.\n",myrank);
      PGFEM_Abort();
    }
  }

  /*=== MULTISCALE INFORMATION ===*/
  if(options.multi_scale){
    sup->multi_scale = options.multi_scale;
    int ms_err = read_interface_macro_normal_lc(options.ipath,sup);
    if(ms_err != 0){
      PGFEM_printerr("[%d] ERROR: could not read normal from file!\n"
	      "Check that the file \"%s/normal.in\""
	      " exists and try again.\n",
	      myrank,options.ipath);
      PGFEM_Abort();
    }
  }

  /*=== BOUNDING ELEMENTS ===*/
  /* NOTE: These might be ripped out... */
  /* ADDED/tested 12/18/2012 MM */
  {
    char bnd_file[500];
    sprintf(bnd_file,"%s%d.in.bnd",in_dat,myrank);
    read_bounding_elements_fname(bnd_file,3,&n_be,&b_elems,mpi_comm);
    bounding_element_set_local_ids(n_be,b_elems,elem);
    bounding_element_reverse_mapping(n_be,b_elems,elem);
  }

  /*==== ADDITIONAL SETUP ===*/

  /* list of elements with prescribed deflection */
  list_el_prescribed_def (sup,node,elem,b_elems,ne,n_be,nn);

  /* list of elements on the COMMUNICATION boundary */
  long *bndel = NULL;
  long nbndel = 0;
  bndel = list_boundary_el(ne,elem,nn,node,myrank,&nbndel);
  
  /*  material matrices of the phases  */
  Mat_3D_orthotropic (nmat,mater,options.analysis_type);
  
  a = aloc3l (nmat,nmat,nc);
  nhommat = list (a,ne,nmat,nc,elem); 
  
  /*  alocation of the material matrices  */
  hommat = build_hommat (nhommat);
  
  /* creates material matrices of the homogeneous medium : LOCAL
     COORDINATE SYSTEM */
  hom_matrices (a,ne,nmat,nc,elem,mater,matgeom,
		hommat,matgeom->SH,options.analysis_type);

  dealoc3l(a,nmat,nmat);
  free(mater);
  
  /* Create Graph for communication */
  /* GrComm = CreateGraph (nproc,myrank,nn,node); */
  
  DomDof = aloc1l (nproc);
  DomNe = aloc1l (nproc);
  DomNn = aloc1l (nproc);
  DomNe[myrank] = ne;
  DomNn[myrank] = nn;
  
  /* Read cohesive elements */
  if (options.cohesive == 1){

    /* read cohesive properties */
    sprintf(filename,"%s%d.in.co_props",in_dat,myrank);
    in1 = PGFEM_fopen(filename,"r");
    read_cohesive_properties(in1,&n_co_props,&co_props,mpi_comm);
    PGFEM_fclose(in1);
      
    /* read coheisve elements */
    sprintf (filename,"%s%d.in.co",in_dat,myrank);
    in1 = PGFEM_fopen(filename,"r");

    /* temporary leftovers from old file format */      
    fscanf (in1,"%ld\n",&ncom); 

    /* to silence warning message. need to pull this legacy bit of
       code out completely. Cohesive porperties provided in separate
       file. This leads to *very* small memory leak */
    if(ncom <= 0) comat = aloc2 (1,4);
    else     comat = aloc2 (ncom,4);

      
    /* read the cohesive element info */
    coel = read_cohe_elem (in1,ncom,ndim,nn,node,&nce,
			   comat,ensight,options.vis_format,
			   myrank,co_props);
    dealoc2 (comat,ncom);
    PGFEM_fclose (in1);

    /* Global number of cohesive elements */
    MPI_Allreduce(&nce,&Gnce,1,MPI_LONG,MPI_SUM,mpi_comm);
  }

  /* use new functions to get code numbers */
  ndofd = generate_local_dof_ids(ne,nce,nn,ndofn,node,
				 elem,coel,b_elems,mpi_comm);
  DomDof[myrank] = generate_global_dof_ids(ne,nce,nn,ndofn,node,
					   elem,coel,b_elems,mpi_comm);

  /* Gather degrees of freedom from all domains */
  MPI_Allgather (MPI_IN_PLACE,1,MPI_LONG,DomDof,1,MPI_LONG,mpi_comm);

  /* Make integer copy of DomDof.  May eventually switch everything to
     integer since 64 bit */
  dist = aloc1i(nproc+1);
  /* build_dist(DomDof,dist,nproc); */
  build_distribution(DomDof,dist,mpi_comm);

  /* Gather number of element from all domains */
  MPI_Gather (&ne,1,MPI_LONG,DomNe,1,MPI_LONG,0,mpi_comm);

  /* Total number of boundary elements */
  long Gnbndel;
  MPI_Reduce(&nbndel,&Gnbndel,1,MPI_LONG,MPI_SUM,0,mpi_comm);

  /* Total number of bounding elements */
  int Gn_be;
  MPI_Reduce(&n_be,&Gn_be,1,MPI_INT,MPI_SUM,0,mpi_comm);

  /* Gather number of nodes from all domains */
  MPI_Gather (&nn,1,MPI_LONG,DomNn,1,MPI_LONG,0,mpi_comm);

  if (myrank == 0 && PFEM_DEBUG) { 
    PGFEM_printf(" Done.\nRedistributing information...");
  }

  for (i=0;i<nproc;i++) {
    Gndof += DomDof[i];
    Gne += DomNe[i];
    /*Gnn += DomNn[i];*/
  }

  /* Compute global matrix row partitioning */
  set_HYPRE_row_col_bounds(PGFEM_hypre,Gndof,DomDof,myrank);

  renumber_global_dof_ids(ne,nce,n_be,nn,ndofn,DomDof,node,
			  elem,coel,b_elems,mpi_comm);
  NBN = distribute_global_dof_ids(ne,nce,n_be,nn,ndofn,ndim,node,
				  elem,coel,b_elems,mpi_comm);

  if (myrank == 0){
    PrintTitleV1();
    for (i = 0; i < argc; i++){
      PGFEM_printf("%s ",argv[i]);
    }
    PGFEM_printf("\n\n");

    switch(options.analysis_type){
    case ELASTIC:
      PGFEM_printf ("ELASTIC ANALYSIS\n");
      break;
    case TP_ELASTO_PLASTIC:
      PGFEM_printf ("TWO PHASE COMPOSITE SYSTEM : ELASTO-PLASTIC ANALYSIS\n");
      break;
    case FS_CRPL:
      PGFEM_printf ("FINITE STRAIN CRYSTAL ELASTO-PLASTICITY\n");
      break;
    case FINITE_STRAIN:
      if (options.cohesive == 0) {
	PGFEM_printf ("FINITE STRAIN ELASTICITY\n");
      } else {
	PGFEM_printf ("FINITE STRAIN ELASTICITY WITH COHESIVE FRACTURE\n");
      }
      break;
    case STABILIZED:
      if (options.cohesive == 0 && gem == 0) {
	PGFEM_printf ("FINITE STRAIN STABILIZED FORMULATION : stb = %12.5e\n",
		options.stab);
      } else if( options.cohesive == 1) {
	PGFEM_printf ("FINITE STRAIN STABILIZED FORMULATION"
		" WITH COHESIVE FRACTURE : stb = %12.5e\n",
		options.stab);
      } else if ( gem == 1) {
	PGFEM_printf ("GENERALIZED FINITE ELEMENT METHOD\n");
	PGFEM_printf ("FINITE STRAIN STABILIZED FORMULATION : stb = %12.5e\n",
		options.stab);
      }
      break;
    case MINI:
      PGFEM_printf("FINITE STRAIN HYPERELASTICITY W/ MINI ELEMENT\n");
      break;
    case MINI_3F:
      PGFEM_printf("FINITE STRAIN HYPERELASTICITY W/ MINI 3 FIELD ELEMENT\n");
      break;
    case DISP:
      PGFEM_printf("FINITE STRAIN DAMAGE HYPERELASTICITY:\n"
	     "TOTAL LAGRANGIAN DISPLACEMENT-BASED ELEMENT\n");
      break;
    case TF:
        PGFEM_printf("FINITE STRAIN TREE FIELDS HYPERELASTICITY:\n"
                "TOTAL LAGRANGIAN TREE FIELDS-BASED ELEMENT\n"); 
      break;                            

    default:
      PGFEM_printerr("ERROR: unrecognized analysis type!\n");
      PGFEM_Abort();
      break;
    }

    if(options.multi_scale){
      if(sup->npd >= 9){
	PGFEM_printf("*** BULK Multiscale Modelling ***\n");
      } else {
	PGFEM_printf("*** INTERFACE Multiscale Modelling ***\n");
      }
    }

    PGFEM_printf ("\n");
    PGFEM_printf ("SolverPackage: ");
    assert(options.solverpackage == HYPRE);
    switch(options.solver){
    case HYPRE_GMRES: PGFEM_printf ("HYPRE - GMRES\n"); break;
    case HYPRE_BCG_STAB: PGFEM_printf ("HYPRE - BiCGSTAB\n"); break;
    case HYPRE_AMG: PGFEM_printf ("HYPRE - BoomerAMG\n"); break;
    case HYPRE_FLEX: PGFEM_printf ("HYPRE - FlexGMRES\n"); break;
    case HYPRE_HYBRID: PGFEM_printf ("HYPRE - Hybrid (GMRES)\n"); break;
    default:
      PGFEM_printerr("Unrecognized solver package!\n");
      PGFEM_Abort();
      break;
    }

    PGFEM_printf("Preconditioner: ");
    switch(options.precond){
    case PARA_SAILS: PGFEM_printf ("HYPRE - PARASAILS\n"); break;
    case PILUT: PGFEM_printf ("HYPRE - PILUT\n"); break;
    case EUCLID: PGFEM_printf ("HYPRE - EUCLID\n"); break;
    case BOOMER: PGFEM_printf ("HYPRE - BoomerAMG\n"); break;
    case NONE: PGFEM_printf ("PGFEM3D - NONE\n"); break;
    case DIAG_SCALE: PGFEM_printf ("PGFEM3D - DIAGONAL SCALE\n"); break;
    case JACOBI: PGFEM_printf ("PGFEM3D - JACOBI\n"); break;
    }
    PGFEM_printf ("\n");
    PGFEM_printf ("Number of total nodes                    : %ld\n",Gnn);
    PGFEM_printf ("Number of nodes on domain interfaces     : %ld\n",NBN);
    PGFEM_printf ("Total number of elements                 : %ld\n",Gne);
    PGFEM_printf ("Number of elems on the COMM interfaces   : %ld\n",Gnbndel);
    PGFEM_printf ("Total number of bounding (surf) elems    : %d\n",Gn_be);
    PGFEM_printf ("Total number of degrees of freedom       : %ld\n",Gndof); 
  }

  {
    /* ALlocate Ap, Ai */
    Ap = aloc1i (DomDof[myrank]+1);
    comm = (COMMUN) PGFEM_calloc (1,sizeof(COMMUN_1));
    initialize_commun(comm);

    Ai = Psparse_ApAi (nproc,myrank,ne,n_be,nn,ndofn,ndofd,
		       elem,b_elems,node,Ap,nce,coel,DomDof,
		       &GDof,comm,mpi_comm,options.cohesive);

    pgfem_comm_build_fast_maps(comm,ndofd,DomDof[myrank],GDof);

    /* Total number of nonzeros and skyline */
    MPI_Reduce (&Ap[DomDof[myrank]],&APP,1,MPI_INT,MPI_SUM,0,mpi_comm);
    temp_int = skyline((int) DomDof[myrank],Ap,Ai,dist[myrank]);
    MPI_Reduce (&temp_int,&sky,1,MPI_INT,MPI_SUM,0,mpi_comm);

    if (myrank == 0){
      PGFEM_printf ("Total number of nonzeros in the matrix   : %d\n",APP);
      if (options.cohesive == 1){
	PGFEM_printf ("Number of cohesive elements              : %ld\n",Gnce);
      }
      PGFEM_printf ("Symmetric skyline (including diagonal)   : %ld\n",sky);
    }

    /*=== NO RENUMBERING === */

    /* allocate vectors */
    r = aloc1 (ndofd);
    f = aloc1 (ndofd);
    d_r = aloc1 (ndofd);
    rr = aloc1 (ndofd);
    D_R = aloc1(ndofd);
    R = aloc1 (ndofd);
    f_defl = aloc1 (ndofd);
    RR = aloc1 (ndofd); 
    f_u = aloc1 (ndofd);
    RRn = aloc1 (ndofd);

    if(sup->npd > 0){
      sup_defl = aloc1 (sup->npd);
    } else {
      sup_defl = NULL;
    }

    U = aloc1 (ndofd);
    DK = aloc1 (ndofd);
    dR = aloc1 (ndofd);
    BS_f = aloc1 (DomDof[myrank]);
    BS_x = aloc1 (DomDof[myrank]);
    BS_RR = aloc1 (DomDof[myrank]);
    BS_d_r = aloc1 (DomDof[myrank]);
    BS_D_R = aloc1 (DomDof[myrank]);
    BS_rr = aloc1 (DomDof[myrank]);
    BS_R = aloc1 (DomDof[myrank]); 
    BS_f_u = aloc1 (DomDof[myrank]);
    BS_U = aloc1 (DomDof[myrank]);
    BS_DK = aloc1 (DomDof[myrank]);
    BS_dR = aloc1 (DomDof[myrank]);

    /*=== TESTING ===*/
    double *nodal_forces = PGFEM_calloc(ndofd,sizeof(double));
    int n_feats = 0;
    int n_sur_trac_elem = 0;
    SUR_TRAC_ELEM *ste = NULL;
    {
      int *feat_type = NULL;
      int *feat_id = NULL;
      double *loads = NULL;

      char *trac_fname = NULL;
      alloc_sprintf(&trac_fname,"%s/traction.in",options.ipath);

      read_applied_surface_tractions_fname(trac_fname,&n_feats,
    					   &feat_type,&feat_id,&loads);

      generate_applied_surface_traction_list(ne,elem,
    					     n_feats,feat_type,
    					     feat_id,&n_sur_trac_elem,
    					     &ste);

      compute_applied_traction_res(ndofn,node,elem,
    				   n_sur_trac_elem,ste,
    				   n_feats,loads,
    				   nodal_forces);

      double tmp_sum = 0.0;
      for(int i=0; i<ndofd; i++){
    	tmp_sum += nodal_forces[i];
      }
 
      MPI_Allreduce(MPI_IN_PLACE,&tmp_sum,1,MPI_DOUBLE,
    		    MPI_SUM,mpi_comm);

      if(myrank == 0){
    	PGFEM_printf("Total load from surface tractions: %.8e\n\n",tmp_sum);
      }

      free(feat_type);
      free(feat_id);
      free(loads);
      free(trac_fname);
    }

    /* push nodal_forces to s->R */
    vvplus  (R,nodal_forces,ndofd);

    /*=== READ SOLVER FILE ===*/
    /* override the default solver file with one specified
       at commandline */
    if(options.override_solver_file){
      if(myrank == 0){
	PGFEM_printf("Overriding the default solver file with:\n%s\n",
	       options.solver_file);
      }
      in1 = PGFEM_fopen(options.solver_file,"r");
    } else {
      /* use the default file/filename */
      sprintf (filename,"%s%d.in.st",in_dat,myrank);
      in1 = PGFEM_fopen(filename,"r");
    }

    fscanf (in1,"%lf %ld %ld %ld",&nor_min,&iter_max,&npres,&FNR);
    if (FNR == 2 || FNR == 3){
      fscanf (in1,"%lf %lf",&dAL0,&dALMAX);
    }
    
    /* Nonlinear solver */
    if (myrank == 0) {
      if (FNR == 0 || FNR == 1) {
	PGFEM_printf ("\nNONLINEAR SOLVER : NEWTON-RAPHSON METHOD\n");
      }
      if ((FNR == 2 || FNR == 3) && ARC == 0){
	PGFEM_printf ("\nNONLINEAR SOLVER : ARC-LENGTH METHOD - Crisfield\n");
      }
      if ((FNR == 2 || FNR == 3) && ARC == 1) {
	PGFEM_printf ("\nNONLINEAR SOLVER : ARC-LENGTH METHOD - Simo\n");
      }
    }

    /* HYPRE INITIALIZATION ROUTINES */
    if(options.solverpackage == HYPRE){ /* HYPRE */
      /* Initialize HYPRE */
      hypre_initialize(Ap,Ai,DomDof[myrank],ni,err,PGFEM_hypre,
		       &options,mpi_comm);
    }

    /*=== CRYSTAL PLASTICITY ===*/
    if (options.analysis_type == FS_CRPL) {
      crpl = (CRPL*) PGFEM_calloc (nmat,sizeof(CRPL));
      read_cryst_plast (in1,nmat,crpl,options.plc);
    }

    /* read number of computational times */
    fscanf (in1,"%ld",&nt);
    
    /* Compute times */
    times = aloc1 (nt+1);
    for (i=0;i<nt+1;i++){
      fscanf (in1,"%lf",&times[i]);
    }
    
    /* read times for output */
    fscanf (in1,"%ld",&n_p);
    
    /* Times for printing */
    print = times_print (in1,nt,n_p);
    
    fscanf (in1,"%ld",&nlod_tim);
    
    /* read times dependent load */
    tim_load = compute_times_load (in1,nt,nlod_tim);
    
    /* alocation of the sigma vector */
    sig_e = build_sig_il (ne,options.analysis_type,elem);

    /* alocation of the sigma vector */
    if (options.smoothing == 0) {
      sig_n = build_sig_el (nn);
    }
    
    /* alocation of the eps vector */
    eps = build_eps_il (ne,elem,options.analysis_type);
    initialize_damage(ne,elem,hommat,eps,options.analysis_type);
  
    /* alocation of pressure variables */
    int nVol = 1;
    switch(options.analysis_type){
      case TF:
        nVol = 1;
        if(elem[0].toe==10 && ndofn==3)
        {  
          npres = 1;
          nVol = 1;
        }          
        break;      
    case STABILIZED: case MINI: case MINI_3F:
      if(npres != 4){
	npres = 4;
	if(myrank == 0){
	  PGFEM_printf("WARNING: Incorrect pressure nodes input, should be 4.\n"
		 "Re-setting to 4 and continuing...\n");
	}
      }
      break;
    case DISP:
      if(npres != 0){
	npres = 0;
	if (myrank == 0) {
	  PGFEM_printf("WARNING: Incorrect pressure nodes input, should be 0.\n"
		 "Re-setting to 0 and continuing...\n");
	}
      }
      break;
    default:
      if(npres != 1){
	npres = 1;
	if (myrank == 0) {
	  PGFEM_printf("WARNING: Incorrect pressure nodes input, should be 1.\n"
		 "Re-setting to 1 and continuing...\n");
	}
      }
      break;
    }/* switch */
    build_pressure_nodes (ne,npres,elem,sig_e,eps,options.analysis_type);
    build_crystal_plast (ne,elem,sig_e,eps,crpl,
			 options.analysis_type,options.plc);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /* \/ initialized element varialbes */    
    if(options.analysis_type==TF)
    {		            
    	for (int e=0;e<ne;e++)
    	{
    	  if(npres==1)
    	  {
        	eps[e].d_T   = (double *) PGFEM_calloc(3,sizeof(double));
      	  for(int a=0; a<3; a++)
      		  eps[e].d_T[a] = 0.0;
        }
    
       	eps[e].T   = (double *) PGFEM_calloc(nVol*3,sizeof(double));	
      	for(int a=0; a<nVol*3; a++)
      		eps[e].T[a] = 1.0;
      }		  
    }    
    /* /\ initialized element varialbes */  			 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

    /* set finite deformations variables */
    set_fini_def (ne,npres,elem,eps,sig_e,options.analysis_type);
    if (options.analysis_type == FS_CRPL){
      set_fini_def_pl (ne,npres,elem,eps,sig_e,crpl,
		       options.analysis_type,options.plc);
    }

    /*  NODE (PRESCRIBED DEFLECTION)- SUPPORT COORDINATES generation
	of the load vector  */
    dt = dt0 = times[1] - times[0];
    if (dt == 0.0){
      if (myrank == 0){
	PGFEM_printf("Incorrect dt\n");
      }
      PGFEM_Comm_code_abort(mpi_comm,0);
    }
    
    load_vec_node_defl (f_defl,ne,ndofn,elem,b_elems,node,hommat,
			matgeom,sup,npres,nor_min,sig_e,eps,dt,
			crpl,options.stab,r,&options);
    
    /*  NODE - generation of the load vector  */
    load_vec_node (R,nln,ndim,znod,node);
    /*  ELEMENT - generation of the load vector  */
    load_vec_elem_sur (R,nle_s,ndim,elem,zele_s);

    /* R   -> Incramental forces 
       RR  -> Total forces for sudivided increment 
       RRn -> Total force after equiblirium */

    vvplus  (f,R,ndofd);
    vvplus  (RR,f,ndofd);
    vvminus (f,f_defl,ndofd);

    /* Transform LOCAL load vector to GLOBAL */
    LToG (R,BS_R,myrank,nproc,ndofd,DomDof,GDof,comm,mpi_comm);

    /* Prescribed deflection */
    for (i=0;i<sup->npd;i++){
      sup_defl[i] = sup->defl_d[i];
    }

    /*=== NO PERIODIC ===*/

    lm = dlm = DLM = DET = 0.0;
    /*PD = 1;*/ 
    dAL = dlm0 = dAL0;
    tim = AT = ITT = 0;

    /* compute un-deformed volume */
    oVolume = 0;
    GVolume = T_VOLUME (ne,ndim,elem,node);
    MPI_Allreduce (&GVolume,&oVolume,1,MPI_DOUBLE,MPI_SUM,mpi_comm);   
    if (myrank == 0){
      PGFEM_printf ("oVolume = %12.12f\n",oVolume);
    }
    VVolume = oVolume;

    /*=== BEGIN SOLVE ===*/
    double *sup_check = NULL;
    if(sup->npd > 0){
      sup_check = aloc1(sup->npd);
    }

    /*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
    
    /* this is for inertia */
    /* material density*/
    double *rho;
    double alpha = 0.5;    /* mid point rule alpha */
    
    double *r_n   = NULL; /* displacement at time is t_n*/
    double *r_n_1 = NULL; /* displacement at time is t_n-1*/
    double *r_n_dof = NULL;
    
    r_n   = aloc1(nn*ndofn);
    r_n_1 = aloc1(nn*ndofn);
    r_n_dof = aloc1(ndofd);
        
    rho = malloc(sizeof(double)*nmat);    
    int restart_tim = 0;
    
    alpha = read_initial_values(r_n_1, r_n, rho, &options, myrank, nn, nmat, times[1] - times[0], &restart_tim);
    for(long a = 0; a<nn; a++)
    {
      for(long b = 0; b<ndofn; b++)
      {
        long id = node[a].id[b];
        if(id>0)
          r[id-1] = r_n[a*ndofn + b];
      }   
    }

    for(int a = 0; a<nmat; a++)
      hommat[a].density = rho[a];      

/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
  

    while (nt > tim){
      dt = times[tim+1] - times[tim];
      if (dt <= 0.0){
	if (myrank == 0) {
	  PGFEM_printf("Incorrect dt\n");
	}
	PGFEM_Comm_code_abort(mpi_comm,0);
      }

      if (myrank == 0){
	PGFEM_printf("\nFinite deformations time step %ld) "
	       " Time %e | dt = %e\n",
	       tim,times[tim+1],dt);
      }

      /*=== NEWTON RAPHSON ===*/
      if (FNR == 0 || FNR == 1){
	if (tim_load[tim] == 1 && tim == 0) {
	  if (myrank == 0){
	    PGFEM_printf ("Incorrect load input for Time = 0\n");
	  }
	  PGFEM_Comm_code_abort (mpi_comm,0);
	}
	if (tim_load[tim] == 1 && tim != 0) {
	  /*  read nodal prescribed deflection */
	  for (i=0;i<sup->npd;i++){
	    fscanf (in1,"%lf",&sup->defl_d[i]);
	    sup_defl[i] = sup_check[i] = sup->defl_d[i];
	  }
	  /* read nodal load in the subdomain */
	  read_nodal_load (in1,nln,ndim,znod);
	  /* read elem surface load */
	  read_elem_surface_load (in1,nle_s,ndim,elem,zele_s);
	  /*  NODE - generation of the load vector  */
	  load_vec_node (R,nln,ndim,znod,node);
	  /*  ELEMENT - generation of the load vector  */
	  load_vec_elem_sur (R,nle_s,ndim,elem,zele_s);

	  /*
	    R   -> Incramental forces 
	    RR  -> Total forces for sudivided increment 
	    RRn -> Total force after equiblirium
	  */

	  /* push nodal_forces to s->R */
	  vvplus  (R,nodal_forces,ndofd);

	} /* end load increment */
	int n_step = 0;
	
/////////////////////////////////////////////////////////////////////////////////////////////////
        if(tim<restart_tim+1)
        {
          tim++;
          continue;
        }
/////////////////////////////////////////////////////////////////////////////////////////////////
		
	hypre_time += Newton_Raphson ( 1,&n_step,ne,n_be,nn,ndofn,ndofd,npres,tim,
				       times,nor_min,dt,elem,b_elems,node,
				       sup,sup_defl,hommat,matgeom,sig_e,
				       eps,Ap,Ai,r,f,d_r,rr,R,f_defl,RR,
				       f_u,RRn,crpl,options.stab,
				       nce,coel,FNR,&pores,PGFEM_hypre,
				       BS_x,BS_f,BS_RR,gama,GNOR,nor1,err,
				       BS_f_u,DomDof,comm,GDof,nt,iter_max,
				       &NORM,nbndel,bndel,mpi_comm,VVolume,
				       &options,NULL,alpha, r_n, r_n_1);

	/* Null global vectors */
	for (i=0;i<ndofd;i++){
	  RRn[i] += R[i];
	  RR[i] = RRn[i];
	  R[i] = 0.0;
	}

	/* null the prescribed displacement increment */
	nulld(sup_defl,sup->npd);
      }/* end NR */

      /*=== ARC LENGTH ===*/
      if (FNR == 2 || FNR == 3){
	dlm = Arc_length ( ne,n_be,nn,ndofn,ndofd,npres,nt,tim,times,
			   nor_min,iter_max,dt,dt0,elem,b_elems,nbndel,
			   bndel,node,sup,sup_defl,hommat,matgeom,sig_e,
			   eps,Ap,Ai,PGFEM_hypre,RRn,f_defl,crpl,
			   options.stab,
			   nce,coel,r,f,d_r,D_R,rr,R,RR,f_u,U,DK,dR,
			   BS_f,BS_d_r,BS_D_R,BS_rr,BS_R,BS_RR,BS_f_u,
			   BS_U,BS_DK,BS_dR,FNR,lm,dAL0,&DET,&dlm0,&DLM,
			   options.vis_format,options.smoothing,
			   sig_n,out_dat,print,&AT,ARC,
			   (times[tim+1]-times[tim])/dt0*dALMAX,&ITT,
			   &dAL,&pores,DomDof,GDof,comm,err,&NORM,
			   mpi_comm,VVolume,&options);

	/* Load multiplier */
	lm += dlm;
	dlm = 0.0;
	
	/* Total force vector */
	for (i=0;i<ndofd;i++){
	  RR[i] = lm*R[i];
	}
      }/* end AL */

      /*=== OUTPUT ===*/
      /* Calculating equvivalent Mises stresses and strains vectors */
      Mises (ne,sig_e,eps,options.analysis_type);

      /* print tractions on marked features */
      {
	double *sur_forces = NULL;
	if(n_feats > 0){
	  sur_forces = PGFEM_calloc(n_feats*ndim,sizeof(double));
	  compute_resultant_force(n_feats,n_sur_trac_elem,
				  ste,node,elem,
				  sig_e,eps,sur_forces);
	  MPI_Allreduce(MPI_IN_PLACE,sur_forces,n_feats*ndim,
			MPI_DOUBLE,MPI_SUM,mpi_comm);
	  if(myrank == 0){
	    PGFEM_printf("Forces on marked features:\n");
	    print_array_d(PGFEM_stdout,sur_forces,n_feats*ndim,
			  n_feats,ndim);
	  }
	}
	free(sur_forces);
      }

      /* Calculate macro deformation gradient */
      GF = computeMacroF(elem,ne,node,nn,eps,oVolume,mpi_comm);
      GS = computeMacroS(elem,ne,node,nn,sig_e,oVolume,mpi_comm);
      GP = computeMacroP(elem,ne,node,nn,sig_e,eps,oVolume,mpi_comm);           

      /* print GF & GS to file */
      if(myrank == 0){                
	sprintf(filename,"%s_macro.out.%ld",out_dat,tim);
	out = fopen(filename,"w");
	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GF[0],GF[1],GF[2]);
	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GF[3],GF[4],GF[5]);
	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GF[6],GF[7],GF[8]);
	PGFEM_fprintf(out,"\n");
	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\t",GS[0],GS[1],GS[2]);
	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GS[3],GS[4],GS[5]);
	PGFEM_fprintf(out,"\n");
	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GP[0],GP[1],GP[2]);
	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GP[3],GP[4],GP[5]);
	PGFEM_fprintf(out,"%8.8e\t%8.8e\t%8.8e\n",GP[6],GP[7],GP[8]);
	fclose(out);
      }

      free(GF);
      free(GS);
      free(GP);

      if (print[tim] == 1 && options.vis_format != VIS_NONE ) {
	if(options.ascii){
	  ASCII_output(&options,mpi_comm,tim,times,Gnn,nn,ne,nce,ndofd,
		       DomDof,Ap,FNR,lm,pores,VVolume,node,elem,sup,
		       r,eps,sig_e,sig_n,coel);
	} /* End ASCII output */

	switch(options.vis_format){
	case VIS_ELIXIR:/* Print to elix file */
	  sprintf (filename,"%s_%d.elx%ld",out_dat,myrank,tim);
	  elixir (filename,nn,ne,ndim,node,elem,sup,r,sig_e,
		  sig_n,eps,options.smoothing,nce,coel,&options);
	  break;
	case VIS_ENSIGHT:/* Print to EnSight files */
	  sprintf (filename,"%s",out_dat);
	  EnSight (filename,tim,nt,nn,ne,ndim,node,elem,sup,
		   r,sig_e,sig_n,eps,options.smoothing,nce,coel,
		   /*nge,geel,ngn,gnod,*/FNR,lm,ensight,mpi_comm,
		   &options);
	  break;
	case VIS_VTK:/* Print to VTK files */
	  if(myrank == 0){
	    VTK_print_master(options.opath,options.ofname,
			     tim,nproc,&options);
	  }

	  VTK_print_vtu(options.opath,options.ofname,tim,
			myrank,ne,nn,node,elem,sup,r,sig_e,eps,
			&options);

///////////////////////////////////////////////////////////////////////////////////      
///////////////////////////////////////////////////////////////////////////////////                          
// this is only for the restart
              if(SAVE_RESTART_FILE)
              {
                for(long a = 0; a<nn; a++)
                {              
                  for(long b = 0; b<ndofn; b++)
                  {
                    long id = node[a].id[b];
                    if(id>0)
                      r_n_dof[id-1] = r_n_1[a*ndofn + b];
                  }
                }
                              
                char restart_path[1024];
                sprintf(restart_path, "%s/restart", options.opath);
                
                if(options.analysis_type==DISP)
                {    
                  if(make_path(restart_path,DIR_MODE) != 0)
                  {
                    PGFEM_printf("Directory (%s) not created!\n",restart_path);
                    abort();                   
                  }
                  write_restart_disp(r_n_1, r_n, &options, myrank, nn, ndofn, tim);
                }                
                else
                {
                  VTK_print_vtu(restart_path,options.ofname,tim,
                      myrank,ne,nn,node,elem,sup,r_n_dof,sig_e,eps,
                      &options);
                }                
              } 
///////////////////////////////////////////////////////////////////////////////////      
///////////////////////////////////////////////////////////////////////////////////  

	  if (options.cohesive == 1){
	    if(myrank == 0){
	      VTK_print_cohesive_master(options.opath,
					options.ofname,tim,nproc,
					&options);
	    }

	    VTK_print_cohesive_vtu(options.opath,options.ofname,
				   tim,myrank,nce,node,coel,sup,r,ensight,
				   &options);
	  }
	  break;
	default: /* no output */ break;
	}/* switch(format) */

      }/* end output */
 
      if (myrank == 0){
	PGFEM_printf("\n");
	PGFEM_printf("*********************************************\n");
	PGFEM_printf("*********************************************\n");
      }
      
      tim++;
    }/* end while */
    
    /*=== FREE MEMORY ===*/
    free(sup_check);
    fclose (in1);
    free(Ap);
    dealoc1i (Ai);
    dealoc1 (f);
    dealoc1 (d_r);
    dealoc1 (rr);
    dealoc1 (R);
    dealoc1 (f_defl);
    dealoc1 (RR);
    dealoc1 (f_u);
    dealoc1 (RRn);
    dealoc1 (sup_defl);
    dealoc1 (times);
    dealoc1l (tim_load);
    dealoc1l (print);
    dealoc1 (U);
    dealoc1 (dR);
    dealoc1 (DK);
    dealoc1 (D_R);
    free(bndel);
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
    dealoc1(r_n);
    dealoc1(r_n_1);
    dealoc1(r_n_dof);
////////////////////////////////////////////////////////////////////////////////////////////////////////
    /* HYPRE */
    if (options.solverpackage == HYPRE){
      destroy_PGFEM_HYPRE_solve_info(PGFEM_hypre);
    }
    
    dealoc1 (BS_x);
    dealoc1 (BS_f);
    dealoc1 (BS_RR);
    dealoc1 (BS_d_r); 
    dealoc1 (BS_D_R);
    dealoc1 (BS_rr);
    dealoc1 (BS_R);
    dealoc1 (BS_f_u);
    dealoc1 (BS_U);
    dealoc1 (BS_DK);
    dealoc1 (BS_dR);
    destroy_applied_surface_traction_list(n_sur_trac_elem,ste);
  }

  /* Deallocate */
  dealoc1 (r);
  dealoc1l (DomDof);
  dealoc1l (DomNe);
  dealoc1l (DomNn);
  free(dist);

  destroy_zatnode(znod,nln);
  destroy_zatelem(zele_s,nle_s);
  destroy_zatelem(zele_v,nle_v);
  destroy_matgeom(matgeom,np);
  destroy_hommat(hommat,nhommat);

  destroy_eps_il(eps,elem,ne,options.analysis_type);
  destroy_sig_il(sig_e,elem,ne,options.analysis_type);

  if(options.cohesive == 1){
    destroy_coel(coel,nce);
    destroy_cohesive_props(n_co_props,co_props);
  }

  destroy_supp(sup);
  destroy_commun(comm,nproc);
  destroy_ensight(ensight);

  destroy_elem(elem,ne);
  destroy_bounding_elements(n_be,b_elems);
  destroy_node(nn,node);

  /*=== PRINT TIME OF ANALYSIS ===*/
  total_time += MPI_Wtime();
  if (myrank == 0){
    getrusage (RUSAGE_SELF, &usage);
    PGFEM_printf ("\n");  
    PGFEM_printf ("Time of analysis on processor [%d] - "
	    " System %ld.%ld, User %ld.%ld.\n",
	    myrank,usage.ru_stime.tv_sec,usage.ru_stime.tv_usec,
	    usage.ru_utime.tv_sec,usage.ru_utime.tv_usec);

    PGFEM_printf ("Total HYPRE solve time on processor [%d] - %f\n",
	    myrank,hypre_time);
    PGFEM_printf("Total time (no MPI_Init()) - %f\n",total_time);
  }

  /*=== FINALIZE AND EXIT ===*/
  PGFEM_finalize_io();
  
  int flag_MPI_finalized;
  MPI_Finalized(&flag_MPI_finalized);
  if(!flag_MPI_finalized)
  {
    if(myrank==0)
      printf("MPI finalizing\n");  

    MPI_Finalize();  
  }
    
  return(0);
}
