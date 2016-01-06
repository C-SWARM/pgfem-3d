#include "allocation.h"
#include "homogen.h"
#include "utils.h"
#include "constitutive_model.h"

#include "read_input_file.h"
#include "post_processing.h"
#include "enumerations.h"
#include "PGFem3D_to_VTK.hpp"
#include "restart.h"
#include "math.h"

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

  
  double *u0 = aloc1(nn*ndofn);
  double *u1 = aloc1(nn*ndofn);
  double *sigma_r = aloc1(1000);  
  
  double temp[6];
  if(myrank==0)
  {
    constitutive_model_test(hommat, NULL, -1); // F(t)
    FILE *fp_r = fopen("single_crystal_results.txt", "r");
    for(int a=0; a<1000; a++)
    {
      fscanf(fp_r, "%lf %lf %lf %lf %lf %lf", temp+0, temp+1, temp+2, temp+3, temp+4, temp+5);
      sigma_r[a] = temp[1];
    }
    fclose(fp_r);
  }
  int nsd = 3;
  int npres=0;
  double dt = 0.1;
  
  double G_gn = 0.0;
  Matrix(double) PK2,sigma,Feff,Eeff,eFeff,E,PK2dev,sigma_dev,eFeffPK2;
  double Err_of_stress = 0.0;  

  for(int istep=0; istep<1000; istep++)
  {
    read_restart(u0,u1,&options,elem,node,NULL,eps,sup,
                 myrank,ne,nn,nsd,&istep); 

    Matrix_construct_init(double, PK2, 3,3,0.0);
    Matrix_construct_init(double, sigma, 3,3,0.0);
    Matrix_construct_init(double, Feff, 3,3,0.0);
    Matrix_construct_init(double, Eeff, 3,3,0.0);              
    Matrix_construct_init(double, eFeff, 3,3,0.0);
    Matrix_construct_init(double, E, 3,3,0.0);
    Matrix_construct_init(double, PK2dev, 3,3,0.0);
    Matrix_construct_init(double, sigma_dev, 3,3,0.0);
    Matrix_construct_init(double, eFeffPK2, 3,3,0.0);
            
    post_processing_compute_stress(PK2.m_pdata,elem,hommat,ne,npres,node,eps,u1,ndofn,mpi_comm, &options);
    post_processing_deformation_gradient(Feff.m_pdata,elem,hommat,ne,npres,node,eps,u1,ndofn,mpi_comm, &options);
    post_processing_deformation_gradient_elastic_part(eFeff.m_pdata,elem,hommat,ne,npres,node,eps,u1,ndofn,mpi_comm, &options);              
    post_processing_plastic_hardness(&G_gn,elem,hommat,ne,npres,node,eps,u1,ndofn,mpi_comm, &options); 
                           
    if(myrank==0)
    { 
      Matrix_eye(Eeff, 3);
      Matrix_AxB(Eeff,0.5,-0.5,Feff,1,Feff,0);
      double det_Fe;
      Matrix_det(eFeff, det_Fe);
      Matrix_AxB(eFeffPK2,1.0,0.0,eFeff,0,PK2,0);
      Matrix_AxB(sigma,1.0/det_Fe,0.0,eFeffPK2,0,eFeff,1);        
                      
      double trPK2, tr_sigma;
      
      Matrix_trace(PK2,trPK2);
      Matrix_trace(sigma,tr_sigma);
      Matrix_eye(PK2dev, 3);
      Matrix_eye(sigma_dev, 3);  
      
      Matrix_AplusB(PK2dev,    1.0, PK2,      -trPK2/3.0, PK2dev);
      Matrix_AplusB(sigma_dev, 1.0, sigma, -tr_sigma/3.0, sigma_dev);    

      double norm_sigma, norm_PK2;
      Matrix_ddot(PK2dev,PK2dev,norm_PK2);    
      Matrix_ddot(sigma_dev,sigma_dev,norm_sigma);

      double sigma_eff=sqrt(3.0/2.0*norm_sigma);
      double PK2_eff = sqrt(3.0/2.0*norm_PK2);    
                    
      FILE *fp_ss;
      if(istep==0) 
        fp_ss = fopen("strain_stress_crystal_plasticity.txt", "w");
      else
        fp_ss = fopen("strain_stress_crystal_plasticity.txt", "a");
        
      fprintf(fp_ss,"%e %e %e %e %e %e\n",dt*(istep+1),sigma_eff,PK2_eff, G_gn, Mat_v(Eeff,1,1), Mat_v(PK2,1,1));                                    

      fclose(fp_ss);
      
      // compute error
      double Err_pt = sqrt((sigma_r[istep] - sigma_eff)*(sigma_r[istep] - sigma_eff))/sigma_r[istep];
      if(Err_pt > Err_of_stress)
        Err_of_stress = Err_pt;
    }                  
  }

  Matrix_cleanup(PK2);
  Matrix_cleanup(sigma);
  Matrix_cleanup(Feff);
  Matrix_cleanup(Eeff);              
  Matrix_cleanup(eFeff);
  Matrix_cleanup(E);
  Matrix_cleanup(PK2dev);
  Matrix_cleanup(sigma_dev);
  Matrix_cleanup(eFeffPK2);                  
  
  free(u0);    
  free(u1);      
  free(sigma_r);
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
  if(myrank==0)
  { 
    FILE *fp_err = fopen("order_of_error.txt", "w");
    fprintf(fp_err,"%d\n", (int) log10(Err_of_stress/10.0));
    fclose(fp_err);
  }
  return(0);
}
