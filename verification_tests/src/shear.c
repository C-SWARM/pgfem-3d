#include "allocation.h"
#include "homogen.h"
#include "utils.h"

#include "read_input_file.h"
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
  if(myrank == 0){
    print_options(stdout,&options);
  }  

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
  in_err = read_input_file(&options,mpi_comm,&nn,&Gnn,&ndofn,
		     &ne,&ni,&err,&limit,&nmat,&nc,&np,&node,
		     &elem,&mater,&matgeom,&sup,&nln,&znod,
		     &nle_s,&zele_s,&nle_v,&zele_v);
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
  
/////////////////////////////////////////////////////////////////////////////////////
// read inputs
  double *u = aloc1(nn*ndofn); 
  read_from_VTK(&options, myrank, 0, u);
  
  int npres = 0;
  double *GS = aloc1(9);    
  post_processing_compute_stress(GS,elem,hommat,ne,npres,node,eps,u,ndofn,mpi_comm, options.analysis_type);            

  printf("computed stress\n");
    
  FILE *fp = fopen("stress.out", "w");  
  fprintf(fp, "%e\n", GS[1]);
  fclose(fp);  
  free(u);    
  free(GS);

  destroy_zatnode(znod,nln);
  destroy_zatelem(zele_s,nle_s);
  destroy_zatelem(zele_v,nle_v);
  destroy_matgeom(matgeom,np);
  destroy_hommat(hommat,nhommat);
  destroy_eps_il(eps,elem,ne,options.analysis_type);
  destroy_supp(sup);
  destroy_elem(elem,ne);
  destroy_node(nn,node);

  /*=== FINALIZE AND EXIT ===*/
  PGFEM_finalize_io();
  MPI_Finalize(); 
  return(0);
}
