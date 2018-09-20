#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "L2norm.h"
#include "PGFem3D_to_VTK.hpp"
#include "allocation.h"
#include "homogen.h"
#include "post_processing.h"
#include "read_input_file.h"
#include "utils.h"
#include "pgfem3d/Communication.hpp"
#include "constitutive_model.h"

using namespace pgfem3d;
using namespace pgfem3d::net;

/*****************************************************/
/*           BEGIN OF THE COMPUTER CODE              */
/*****************************************************/

int main(int argc,char *argv[])
{  
  Boot *boot = new Boot();
  int myrank = boot->get_rank();
  int nproc = boot->get_nproc();

  PGFem3D_opt options;
  if (argc <= 2){
    if(myrank == 0){
      print_usage(stdout);
    }
    exit(0);
  }
  set_default_options(&options);
  re_parse_command_line(myrank,2,argc,argv,&options);
  
  // Create the desired network
  Network *net = Network::Create(options);
  
  CommunicationStructure *com = new CommunicationStructure();
  com->rank = myrank;
  com->nproc = nproc;
  com->net = net;
  com->comm = NET_COMM_WORLD;

  char processor_name[NET_MAX_PROCESSOR_NAME];
  int namelen = 0;
  boot->get_processor_name(processor_name, &namelen);
  
  PGFEM_initialize_io(NULL,NULL);

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
  
  char physicsname[1024] = "Mechanical";
  char *p = physicsname;
  in_err = read_input_file(&options,com,&nn,&Gnn,&ndofn,
                           &ne,&ni,&err,&limit,&nmat,&nc,&np,&node,
                           &elem,&mater,&matgeom,&sup,&nln,&znod,
                           &nle_s,&zele_s,&nle_v,&zele_v,&fv_ndofn,
			                     physicsno,&ndim,&p);
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
  eps = build_eps_il(ne,elem,options.analysis_type,NULL);

  // read time steps
  char filename[2048];
  char in_dat[1024];
  FILE *fp = NULL;  
  double nor_min;
  long iter_max, npres, FNR, nt;
  
  if(options.override_solver_file)
  {
    if(myrank == 0)
      PGFEM_printf("Overriding the default solver file with:\n%s\n", options.solver_file);

    fp = fopen(options.solver_file,"r");
  }
  else
  {
    // use the default file/filename
    sprintf(in_dat,"%s/%s",options.ipath,options.ifname);
    sprintf(filename,"%s%d.in.st",in_dat,myrank);
    fp = fopen(filename,"r");
    if(fp==NULL)
    {
      sprintf(filename,"%s%d.in.st",in_dat,0);
      fp = fopen(filename,"r");
    }
  }  
  
  CHECK_SCANF(fp,"%lf %ld %ld %ld", &nor_min,&iter_max,&npres,&FNR);
  CHECK_SCANF(fp,"%ld", &nt);

  /* Compute times */
  double *times = (double *) malloc((nt+1)*sizeof(double));
  for(int a=0;a<nt+1;a++){
    CHECK_SCANF(fp,"%lf",times+a);
  }

  long n_p;
  /* read times for output */
  CHECK_SCANF(fp,"%ld",&n_p);

  /* Times for printing */
  long *print = times_print(fp,nt,n_p);

  fclose(fp);

  const int error_no = 3;
  double sum_err[error_no] = {};
  
  double *u = aloc1(nn*ndofn);
  double *Ph = aloc1(ne);
  double *Vh = aloc1(ne);

  bool is4cm = false;
  
  if(options.analysis_type == CM || options.analysis_type == CM3F){
    char cm_filename[2048];
    sprintf(cm_filename, "%s/model_params.in", options.ipath);
    FILE *cm_in = PGFEM_fopen(cm_filename, "r"); 
              
    int num_entries = -1;
    in_err += scan_for_valid_line(cm_in);
    CHECK_SCANF(cm_in, "%d", &num_entries);

    int mat_id = 0;
    int model_type = -1;
    err += scan_for_valid_line(cm_in);

    CHECK_SCANF(cm_in, "%d %d", &mat_id, &model_type);
    fclose(cm_in);
    
    if(model_type == MANUFACTURED_SOLUTIONS)
      is4cm = true;
  }       
  
  FILE *f = NULL;
  if(myrank==0)
    f = fopen("error.txt", "w");

  for(int tim=0; tim<nt; tim++){
    if(print[tim]){
      double err[error_no] = {};
      sprintf(filename,"%s/VTK/STEP_%.6d/%s_%d_%d.vtu",options.opath,tim,options.ofname,myrank,tim);
      read_VTK_file4TF(filename, u,Ph,Vh);

      compute_L2_error(err, elem, ne, node, u, Ph, Vh, times[tim+1], com, &options, hommat, is4cm);

      if(myrank==0)
      {
        PGFEM_printf("error = %e, %e, %e\n", sqrt(err[0]), sqrt(err[1]), sqrt(err[2]));
        PGFEM_fprintf(f, "%e %e %e\n", sqrt(err[0]), sqrt(err[1]), sqrt(err[2]));
      }
      for(int ib=0; ib<error_no; ib++)
        sum_err[ib] += err[ib];
    }
  }

  if(myrank==0)
  {
    PGFEM_printf("total error = %e, %e, %e\n", sqrt(sum_err[0]), sqrt(sum_err[1]), sqrt(sum_err[2]));
    PGFEM_fprintf(f, "%e %e %e\n", sqrt(sum_err[0]), sqrt(sum_err[1]), sqrt(sum_err[2]));
    fclose(f);
  }

  free(u);
  free(Ph);
  free(Vh);
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
  net->finalize();

  delete com;
  delete net;
  delete boot;
  
  return(0);
}
