#include "allocation.h"
#include "homogen.h"
#include "utils.h"

#include "read_input_file.h"
#include "L2norm.h"
#include "PGFem3D_to_VTK.hpp"
#include "post_processing.h"
/*****************************************************/
/*           BEGIN OF THE COMPUTER CODE              */
/*****************************************************/

int main(int argc,char *argv[])
{
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  int myrank = 0;
  int nproc = 0;

  /*=== END INITIALIZATION === */
  MPI_Init (&argc,&argv);
  MPI_Comm_rank (mpi_comm,&myrank);
  MPI_Comm_size (mpi_comm,&nproc);
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int namelen = 0;  
  MPI_Get_processor_name (processor_name,&namelen);
  PGFEM_initialize_io(NULL,NULL);  
    
  PGFem3D_opt options;
  if (argc <= 2){
    if(myrank == 0){
      print_usage(stdout);
    }
    exit(0);
  }
  set_default_options(&options);
  re_parse_command_line(myrank,2,argc,argv,&options); 

  long nn = 0;
  long Gnn = 0;
  long ndofn = 0;
  long ne = 0;
  long ni = 0;
  double err = 0.0;
  double limit = 0.0;  
  long nmat = 0;  
  long nc = 0;
  long np = 0;
  NODE *node = NULL;  
  ELEMENT *elem = NULL; 
  MATERIAL *mater = NULL;   
  MATGEOM matgeom = NULL;
  SUPP sup = NULL;
  long nln = 0;  
  ZATNODE *znod = NULL;   
  long nle_s = 0;
  ZATELEM *zele_s = NULL;
  long nle_v = 0;
  ZATELEM *zele_v = NULL;    
  
  
  int in_err = 0;
  int physicsno = 1;
  int ndim = 3;
  in_err = read_input_file(&options,mpi_comm,&nn,&Gnn,&ndofn,
         &ne,&ni,&err,&limit,&nmat,&nc,&np,&node,
         &elem,&mater,&matgeom,&sup,&nln,&znod,
         &nle_s,&zele_s,&nle_v,&zele_v,physicsno,&ndim);
  if(in_err){
    PGFEM_printerr("[%d]ERROR: incorrectly formatted input file!\n",
	    myrank);
    PGFEM_Abort();
  }   

  HOMMAT *hommat = NULL;
  Mat_3D_orthotropic (nmat,mater,options.analysis_type);
  
  long ***a = NULL;
  a = aloc3l (nmat,nmat,nc);
  long nhommat = list(a,ne,nmat,nc,elem); 
  
  /*  alocation of the material matrices  */
  hommat = build_hommat(nhommat);

  hom_matrices(a,ne,nmat,nc,elem,mater,matgeom,
		hommat,matgeom->SH,options.analysis_type);

  dealoc3l(a,nmat,nmat);
  free(mater);  
  
  EPS *eps = NULL;    
  eps = build_eps_il(ne,elem,options.analysis_type);

  // read time steps
  char filename[1024];
  double nor_min;  
  long iter_max, npres, FNR, nt;
  sprintf (filename,"%s/%s%d.in.st",options.ipath,options.ifname,myrank);
  FILE *in_st = fopen(filename,"r");
  fscanf (in_st,"%lf %ld %ld %ld",&nor_min,&iter_max,&npres,&FNR);
  fscanf (in_st,"%ld",&nt);
  
  /* Compute times */
  double *times = (double *) malloc((nt+1)*sizeof(double));
  for(int a=0;a<nt+1;a++){
      fscanf (in_st,"%lf",&times[a]);
  }
  
  long n_p;
  /* read times for output */
  fscanf (in_st,"%ld",&n_p);
  
  /* Times for printing */
  long *print = times_print(in_st,nt,n_p);  
  
  fclose(in_st);  
            
  double GL2_err[3];
  double t = 0.0;
  
  int tim = nt+1;
  while(print[tim]!=1)
  {
    tim--;
  }  

  /////////////////////////////////////////////////////////////////////////////////////
  // read inputs
  double *u = aloc1(nn*ndofn);
  double *Ph = aloc1(ne);
  double *Vh = aloc1(ne);

  sprintf(filename,"%s/VTK/STEP_%.5d/%s_%d_%d.vtu",options.opath,tim,options.ofname,myrank,tim);

  read_VTK_file4TF(filename, u,Ph,Vh);

  /////////////////////////////////////////////////////////////////////////////////////
  // compute errors    
  compute_L2_error(GL2_err, elem, ne, node, u, Ph, Vh, times[tim+1], mpi_comm, &options, hommat);
  if(myrank==0)
  {
    FILE *f = fopen("error.txt", "w");
    printf("time step = %d, time = %e\n", tim, times[tim+1]);
    printf("error = %e, %e, %e\n", sqrt(GL2_err[0]), sqrt(GL2_err[1]), sqrt(GL2_err[2]));
    fprintf(f , "%e, %e, %e\n", sqrt(GL2_err[0]), sqrt(GL2_err[1]), sqrt(GL2_err[2]));
  }

  free(u);
  free(times);
  free(print);    

  destroy_zatnode(znod,nln);
  destroy_zatelem(zele_s,nle_s);
  destroy_zatelem(zele_v,nle_v);
  destroy_matgeom(matgeom,np);
  destroy_hommat(hommat,nhommat);
  destroy_eps_il(eps,elem,ne,options.analysis_type);
  destroy_supp(sup);
  destroy_elem(elem,ne);
  destroy_node_multi_physics(nn,node,physicsno);

  /*=== FINALIZE AND EXIT ===*/
  PGFEM_finalize_io();
  MPI_Finalize(); 
  return(0);
}
