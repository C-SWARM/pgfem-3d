#include "allocation.h"
#include "homogen.h"
#include "utils.h"
#include "constitutive_model.h"

#include "read_input_file.h"
#include "post_processing.h"
#include "enumerations.h"
#include "PGFem3D_to_VTK.hpp"

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
  Model_parameters *param_list = NULL;   
  
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
    
  if (options.analysis_type == CM) {
    /* parameter list and initialize const. model at int points.
     * NOTE: should catch/handle returned error flag...
     */
    char *cm_filename = NULL;
    alloc_sprintf(&cm_filename,"%s/model_params.in",options.ipath);
    FILE *cm_in = PGFEM_fopen(cm_filename, "r");
    read_model_parameters_list(&param_list, nhommat, hommat, cm_in);
    free(cm_filename);
    fclose(cm_in);
    init_all_constitutive_model(eps,ne,elem,nhommat,param_list);
  } 
  

  FILE *in_st;
  char filename[1024];

  if(options.override_solver_file)
  {
    if(myrank == 0)
	    printf("Overriding the default solver file with:\n%s\n", options.solver_file);

    in_st = fopen(options.solver_file,"r");
  }
  else 
  {
    sprintf (filename,"%s/%s%d.in.st",options.ipath,options.ifname,myrank);
    in_st = fopen(filename,"r");  
  }
  
  double nor_min;
  long iter_max, npres, FNR, nt;
  
  fscanf (in_st,"%lf %ld %ld %ld",&nor_min,&iter_max,&npres,&FNR);
  fscanf (in_st,"%ld",&nt);
  
  fclose(in_st);  
/////////////////////////////////////////////////////////////////////////////////////
// read inputs  
  sprintf(filename,"%s/VTK/STEP_%.5ld/%s_%d_%ld.vtu",options.opath,nt-1,options.ofname,myrank,nt-1);   
  
  
  double *u = aloc1(nn*ndofn);
  double *P, *V;

  int nVol = 1;

  if(options.analysis_type==TF)
  {
    if(elem[0].toe==10 && ndofn==3)
    {  
      npres = 1;
      nVol = 1;
    }
    
    P = (double *) malloc(sizeof(double)*ne*npres);
    V = (double *) malloc(sizeof(double)*ne*nVol);
    read_VTK_file4TF(filename, u, P, V);        
    for (int e=0;e<ne;e++)
    {
      if(npres==1)
      {
      	eps[e].d_T   = (double *) PGFEM_calloc(3,sizeof(double));
      	eps[e].d_T[0] = P[e];
      }
    
     	eps[e].T   = (double *) PGFEM_calloc(nVol*3,sizeof(double));	
  		eps[e].T[0] = V[e];
    }
                              
    free(P);
    free(V);          
  }
  else
    read_VTK_file(filename, u);
  
  enum{MAX_H, MAX_DISP, MAX_END};
  double value_max[MAX_END], Gvalue_max[MAX_END];
  
  for(int a=0; a<MAX_END; a++)
  {
    Gvalue_max[a] = -1.0e-15;
     value_max[a] = -1.0e-15;
  }  
  for(int a=0; a<nn; a++)
  {
    double x0 = node[a].x1_fd;
    double y0 = node[a].x2_fd;
    double z0 = node[a].x3_fd;    
    
    double x = u[a*3+0] + x0;
    double y = u[a*3+1] + y0;
    double z = u[a*3+2] + z0;
    
    double disp = sqrt(x*x + y*y);
    if(fabs(x0)<1.0e-12 && z>value_max[MAX_H])
      value_max[MAX_H] = z;
      
    if(disp>value_max[MAX_DISP])
      value_max[MAX_DISP] = disp;
  }
    
  Matrix(double) PK2,sigma,eFeff,sigma_dev,eFeffPK2;

  Matrix_construct_init(double, PK2, 3,3,0.0);
  Matrix_construct_init(double, sigma, 3,3,0.0);
  Matrix_construct_init(double, eFeff, 3,3,0.0);
  Matrix_construct_init(double, sigma_dev, 3,3,0.0);
  Matrix_construct_init(double, eFeffPK2, 3,3,0.0);
            
  post_processing_compute_stress(PK2.m_pdata,elem,hommat,ne,npres,node,eps,u,ndofn,mpi_comm, &options);
  post_processing_deformation_gradient_elastic_part(eFeff.m_pdata,elem,hommat,ne,npres,node,eps,u,ndofn,mpi_comm, &options);              

  double PE = 0.0;
  post_processing_potential_energy(&PE,elem,hommat,ne,npres,node,eps,u,ndofn,mpi_comm, &options);

  MPI_Allreduce(value_max,Gvalue_max,MAX_END,MPI_DOUBLE,MPI_MAX,mpi_comm);

  if(myrank==0)
  {  
    double det_Fe, tr_sigma;
    Matrix_det(eFeff, det_Fe);
    Matrix_AxB(eFeffPK2,1.0,0.0,eFeff,0,PK2,0);
    Matrix_AxB(sigma,1.0/det_Fe,0.0,eFeffPK2,0,eFeff,1);            
    Matrix_trace(sigma,tr_sigma);
    Matrix_eye(sigma_dev, 3);
    Matrix_AplusB(sigma_dev, 1.0, sigma, -tr_sigma/3.0, sigma_dev);  

    double norm_sigma;
    Matrix_ddot(sigma_dev,sigma_dev,norm_sigma);  
    double sigma_eff=sqrt(3.0/2.0*norm_sigma);
        
            
    FILE *fp = fopen("maximum_disp.out", "w");  
    fprintf(fp, "%e %e %e %e\n", Gvalue_max[MAX_DISP], sigma_eff, Gvalue_max[MAX_H], PE);
    fclose(fp);
  } 


  Matrix_cleanup(PK2);
  Matrix_cleanup(sigma);
  Matrix_cleanup(eFeff);
  Matrix_cleanup(sigma_dev);
  Matrix_cleanup(eFeffPK2);
 
  free(u);    
  destroy_zatnode(znod,nln);
  destroy_zatelem(zele_s,nle_s);
  destroy_zatelem(zele_v,nle_v);
  destroy_matgeom(matgeom,np);
  destroy_hommat(hommat,nhommat);
  destroy_model_parameters_list(nhommat,param_list);  
  destroy_eps_il(eps,elem,ne,options.analysis_type);
  destroy_supp(sup);
  destroy_elem(elem,ne);
  destroy_node(nn,node);

  /*=== FINALIZE AND EXIT ===*/
  PGFEM_finalize_io();
  MPI_Finalize(); 
  return(0);
}
