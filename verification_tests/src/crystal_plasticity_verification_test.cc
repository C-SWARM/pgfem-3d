#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "crystal_plasticity_integration.h"
#include "PGFem3D_data_structure.h"
#include "PGFem3D_to_VTK.hpp"
#include "allocation.h"
#include "constitutive_model.h"
#include "enumerations.h"
#include "homogen.h"
#include "math_help.h"
#include "post_processing.h"
#include "read_input_file.h"
#include "restart.h"
#include "utils.h"
#include "elem3d.h"
#include <ttl/ttl.h>
#include <ttl/Library/matrix.h>
#include <cmath>

#define DIM_3       3

using namespace ttl;

namespace {
  const constexpr Index<'i'> i;
  const constexpr Index<'j'> j;
  const constexpr Index<'k'> k;
  const constexpr Index<'l'> l;

  template<class T1, class T2> int inv(T1 &A, T2 &AI)
  {
    int err = inverse(A.data, DIM_3, AI.data);
    return err;
  }
}

/*****************************************************/
/*           BEGIN OF THE COMPUTER CODE              */
/*****************************************************/


void test_crystal_plasticity_single_crystal(void)
{
  //double lame1 = 75600.0;
  //double lame2     = 26100.0;
  double E = 70.0e+3;
  double nu = 0.25;
  
  double gamma_dot_0 = 1.0;
  double gamma_dot_s = 50.0e+9;
  double m           = 0.1;  
  double g0          = 210.0;
  double G0          = 200.0;
  double gs_0        = 330.0;
  double w           = 0.005;
  
  int max_itr_stag      = 100;
  int max_itr_hardening = 5;
  int max_itr_M         = 100;
  double tol_hardening  = 1.0e-6;
  double tol_M          = 1.0e-6;
  double computer_zero  = 1.0e-15;


  // create material properties: Elasticity
  MATERIAL_ELASTICITY mat_e;  
  set_properties_using_E_and_nu(&mat_e,E,nu); 
  // or you can use : set_properties_using_Lame_constants(&mat_e,lame1,lame2);
  //print_material_property_elasticity(&mat_e); // <= this is optional
  
  // create slip system : 0 for FCC
  // it is needed to setup plasticity
  SLIP_SYSTEM slip;
  construct_slip_system(&slip,0); 
  
  // create material properties: Plasticity
  MATERIAL_CRYSTAL_PLASTICITY mat_p;
  set_properties_crystal_plasticity(&mat_p,&slip,gamma_dot_0,gamma_dot_s, 
                                     m,g0,G0,gs_0,w);
  //print_material_property_crystal_plasticity(&mat_p);  // <= this is optional 

  // create material plasticity: it needs material properties for elasticity and plasticity
  MATERIAL_CONSTITUTIVE_MODEL mat;
  set_properties_constitutive_model(&mat,&mat_e,&mat_p);
  
  // create solver info: criteria for numerical iterations
  CRYSTAL_PLASTICITY_SOLVER_INFO solver_info;
  set_crystal_plasticity_solver_info(&solver_info,max_itr_stag,
                                                  max_itr_hardening,
                                                  max_itr_M,
                                                  tol_hardening,
                                                  tol_M,
                                                  computer_zero);  
  //print_crystal_plasticity_solver_info(&solver_info); // <= this is optional
  
  // create elasticity object for integration
  // this creates memory for stress and elasticity tensor s.t. requres destructor
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);  

  // set variables for integration
  ttl::Tensor<2,DIM_3,double> M,MI,pFn,pFnp1,pFnp1_I,eFnp1,Fn,Fnp1 = {},sigma,PK2dev,sigma_dev;
  
  pFn   = ttl::identity(i,j);
  pFnp1 = ttl::identity(i,j);
  eFnp1 = ttl::identity(i,j);
  Fn    = ttl::identity(i,j);
  Fnp1  = ttl::identity(i,j);
      
  double g_n,g_np1;
  g_n = g_np1 = mat_p.g0;
  
  double dt = 0.1;
  
  // start integration  
  ttl::Tensor<2,DIM_3,double*> PK2(elast.S);
  
  FILE *fp = fopen("single_crystal_results.txt", "w");
  
  for(int a = 1; a<=1000; a++)
  {
    double lambda = 0.0;
    double t = a*dt;
    
    Fnp1[0][0] = 1.0 - t*1.0e-3;
    Fnp1[1][1] = Fnp1[2][2] = 1.0 + t*0.5*1.0e-3;
    
    staggered_Newton_Rapson(pFnp1.data,M.data, &g_np1, &lambda, 
                            pFn.data, Fn.data,Fnp1.data, 
                            g_n, dt, &mat, &elast, &solver_info);
                            
    pFn = pFnp1(i,j);
    Fn  = Fnp1(i,j);
    inv(pFnp1,pFnp1_I);
    eFnp1 = Fnp1(i,k)*pFnp1_I(k,j);
       
    g_n = g_np1;
    
    
    // print result at time t
    double sigma_eff;
    double PK2_eff;
    double det_pF = ttl::det(pFnp1);
        
    elast.update_elasticity(&elast,eFnp1.data,0);
    elast.compute_PK2_eff(&elast,&PK2_eff);
    elast.compute_Cauchy_eff(&elast,&sigma_eff,eFnp1.data);   

    fprintf(fp, "%e %e %e %e %e %e\n",t,sigma_eff,PK2_eff, g_np1, det_pF, PK2[0][0]);
  }    
  
  fclose(fp);  
  destruct_elasticity(&elast);
  destruct_slip_system(&slip);
}


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
  Node *node = NULL;
  Element *elem = NULL;
  Material *mater = NULL;
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
  int fv_ndofn = ndim;

  in_err = read_input_file(&options,mpi_comm,&nn,&Gnn,&ndofn,
                           &ne,&ni,&err,&limit,&nmat,&nc,&np,&node,
                           &elem,&mater,&matgeom,&sup,&nln,&znod,
                           &nle_s,&zele_s,&nle_v,&zele_v,&fv_ndofn,physicsno,&ndim,NULL);
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

  int n_state_varialbles = 0;
  for(int eid=0; eid<ne; eid++)
  {
    int nne = elem[eid].toe;
    long nint = 0;
    int_point(nne,&nint);
    n_state_varialbles += nint;
  }

  State_variables *statv_list = (State_variables *) malloc(sizeof(State_variables)*n_state_varialbles);

  eps = build_eps_il(ne,elem,options.analysis_type,&statv_list);

  if (options.analysis_type == CM) {
    /* parameter list and initialize const. model at int points.
     * NOTE: should catch/handle returned error flag...
     */
    char *cm_filename = NULL;
    alloc_sprintf(&cm_filename,"%s/model_params.in",options.ipath);
    FILE *cm_in = PGFEM_fopen(cm_filename, "r");
    read_model_parameters_list(nhommat, hommat, cm_in);
    free(cm_filename);
    fclose(cm_in);
    init_all_constitutive_model(eps,ne,elem,nhommat,hommat);
  }


  double *u0 = aloc1(nn*ndofn);
  double *u1 = aloc1(nn*ndofn);
  double *sigma_r = aloc1(1000);

  double temp[6];
  if(myrank==0)
  {
    test_crystal_plasticity_single_crystal(); // F(t)
    FILE *fp_r = fopen("single_crystal_results.txt", "r");
    for(int a=0; a<1000; a++)
    {
      fscanf(fp_r, "%lf %lf %lf %lf %lf %lf", temp+0, temp+1, temp+2, temp+3, temp+4, temp+5);
      sigma_r[a] = temp[2];
    }
    fclose(fp_r);
  }
  int nsd = 3;
  int npres=0;
  double dt = 0.1;

  double G_gn = 0.0;
  Tensor<2, 3, double> I,PK2,sigma,Feff,Eeff,eFeff,E,PK2dev,sigma_dev;
  I = identity(i,j);

  double Err_of_stress = 0.0;

  double tnm1[2] = {-1.0,-1.0};
  double NORM = 0.0;

  // initialize and define multiphysics
  Multiphysics mp;
  int id = MULTIPHYSICS_MECHANICAL;
  int write_no = 0;
  int *coupled_ids = (int *) malloc(sizeof(int));
  char *physicsname = (char *) malloc(sizeof(char)*1024);
  {
    coupled_ids[0] = 0;
    sprintf(physicsname, "Mechanical");

    mp.physicsno      = 1;
    mp.physicsname    = &physicsname;
    mp.physics_ids    = &id;
    mp.ndim           = &ndim;
    mp.write_no       = &write_no;
    mp.write_ids      = NULL;
    mp.coupled_ids    = &coupled_ids;
    mp.total_write_no = 0;
  }

  Grid grid;
  grid_initialization(&grid);
  {
    grid.ne          = ne;
    grid.nn          = nn;
    grid.nsd         = nsd;
    grid.element     = elem;
    grid.node        = node;
  }

  // initialize and define field variables
  FieldVariables fv;
  {
    field_varialbe_initialization(&fv);
    fv.ndofn  = ndofn;
    fv.npres  = npres;
    fv.eps    = eps;
    fv.u_nm1  = u0;
    fv.u_n    = u1;
    fv.NORM   = NORM;
    fv.statv_list = statv_list;
  }

  double tns[2];
  TimeStepping ts;
  {
    time_stepping_initialization(&ts);
    ts.tns    = tns;
  }

  // initialize and define loading steps object
  LoadingSteps load;
  {
    loading_steps_initialization(&load);
    load.sups = &sup;
  }

  for(int istep=0; istep<1000; istep++)
  {
    options.restart = istep;
    read_restart(&grid,&fv,&ts,&load,&options,&mp,tnm1,myrank);

    post_processing_compute_stress(PK2.data,elem,hommat,ne,npres,node,eps,u1,ndofn,mpi_comm, &options);
    post_processing_deformation_gradient(Feff.data,elem,hommat,ne,npres,node,eps,u1,ndofn,mpi_comm, &options);
    post_processing_deformation_gradient_elastic_part(eFeff.data,elem,hommat,ne,npres,node,eps,u1,ndofn,mpi_comm, &options);
    post_processing_plastic_hardness(&G_gn,elem,hommat,ne,npres,node,eps,u1,ndofn,mpi_comm, &options);

    if(myrank==0)
    {
      Eeff = 0.5*(Feff(k,i)*Feff(k,j) - I(i,j));
      double det_Fe = det(eFeff);

      sigma = 1.0/det_Fe*eFeff(i,k)*PK2(k,l)*eFeff(j,l);

      double trPK2    = PK2(i,i);
      double tr_sigma = sigma(i,i);

      PK2dev    = PK2(i,j) - trPK2/3.0*I(i,j);
      sigma_dev = sigma(i,j) - tr_sigma/3.0*I(i,j);

      double norm_sigma = sigma_dev(i,j)*sigma_dev(i,j);
      double norm_PK2   = PK2dev(i,j)*PK2dev(i,j);

      double sigma_eff=sqrt(3.0/2.0*norm_sigma);
      double PK2_eff = sqrt(3.0/2.0*norm_PK2);

      FILE *fp_ss;
      if(istep==0)
        fp_ss = fopen("strain_stress_crystal_plasticity.txt", "w");
      else
        fp_ss = fopen("strain_stress_crystal_plasticity.txt", "a");

      fprintf(fp_ss,"%e %e %e %e %e %e\n",dt*(istep+1),sigma_eff,PK2_eff, G_gn, Eeff.data[0], PK2.data[0]);

      fclose(fp_ss);

      // compute error
      double Err_pt = sqrt((sigma_r[istep] - PK2_eff)*(sigma_r[istep] - PK2_eff))/sigma_r[istep];
      if(Err_pt > Err_of_stress)
        Err_of_stress = Err_pt;
    }
  }

  free(u0);
  free(u1);
  free(sigma_r);
  destroy_zatnode(znod,nln);
  destroy_zatelem(zele_s,nle_s);
  destroy_zatelem(zele_v,nle_v);
  destroy_matgeom(matgeom,np);
  destroy_hommat(hommat,nhommat);
  destroy_eps_il(eps,elem,ne,options.analysis_type);
  destroy_supp(sup);
  destroy_elem(elem,ne);
  destroy_node_multi_physics(nn,node,physicsno);

  free(coupled_ids);
  free(physicsname);

  if(myrank==0)
  {
    FILE *fp_err = fopen("order_of_error.txt", "w");
    fprintf(fp_err,"%d\n", (int) log10(Err_of_stress/10.0));
    fclose(fp_err);
  }

  /*=== FINALIZE AND EXIT ===*/
  PGFEM_finalize_io();
  MPI_Finalize();

  return(0);
}
