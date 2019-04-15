#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "PGFem3D_to_VTK.hpp"
#include "allocation.h"
#include "constitutive_model.h"
#include "enumerations.h"
#include "homogen.h"
#include "post_processing.h"
#include "read_input_file.h"
#include "utils.h"
#include "pgfem3d/Communication.hpp"

using namespace pgfem3d;
using namespace multiscale;
using namespace multiscale::net;

/*****************************************************/
/*           BEGIN OF THE COMPUTER CODE              */
/*****************************************************/

int main(int argc,char *argv[])
{
  Boot *boot = Boot::Create(BOOT_DEFAULT);
  int myrank = boot->get_rank();
  int nproc = boot->get_ranks();

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
  Network *net = Network::Create(NET_DEFAULT);
  
  CommunicationStructure *com = new CommunicationStructure();
  com->rank = myrank;
  com->nproc = nproc;
  com->boot = boot;
  com->net = net;
  com->comm = NET_COMM_WORLD;

  char processor_name[MAX_PROCESSOR_NAME];
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

  in_err = read_input_file(&options,com,&nn,&Gnn,&ndofn,
                           &ne,&ni,&err,&limit,&nmat,&nc,&np,&node,
                           &elem,&mater,&matgeom,&sup,&nln,&znod,
                           &nle_s,&zele_s,&nle_v,&zele_v,&fv_ndofn,
			   physicsno,&ndim,NULL);
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
    init_all_constitutive_model(eps,ne,elem,nhommat,hommat,myrank);
  }
  /////////////////////////////////////////////////////////////////////////////////////
  // read inputs
  char filename[1024];
  sprintf(filename,"%s/VTK/STEP_%.6d/%s_%d_%d.vtu",options.opath,0,options.ofname,myrank,0);


  double *u = aloc1(nn*ndofn);
  double *P = NULL;
  double *V = NULL;

  int nVol = 1;
  int npres = 1;

  if(options.analysis_type==TF)
  {
    P = (double *) malloc(sizeof(double)*ne*npres);
    V = (double *) malloc(sizeof(double)*ne*nVol);
    read_VTK_file4TF(filename, u, P, V);
  }
  else
    read_VTK_file(filename, u);

  double *GS = aloc1(9);
  post_processing_compute_stress(GS,elem,hommat,ne,npres,node,eps,u,V,ndofn,com, &options);

  if(options.analysis_type==TF)
  {
    free(P);
    free(V);
  }

  FILE *fp = fopen("stress_tension_pressure.out", "w");
  fprintf(fp, "%e\n", GS[0]);
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
  destroy_node_multi_physics(nn,node,physicsno);

  /*=== FINALIZE AND EXIT ===*/
  PGFEM_finalize_io();
  net->finalize();

  delete com;
  delete net;
  delete boot;
  
  return(0);
}
