#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "allocation.h"
#include "constitutive_model.h"
#include "homogen.h"
#include "post_processing.h"
#include "read_input_file.h"
#include "utils.h"
#include "pgfem3d/Communication.hpp"
#include <unistd.h>

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


  // Load_Type = 0: Uniaxial_tension
  // Load_Type = 1: Uniaxial_compression
  // Load_Type = 2: Simple_shear
  // Load_Type = 3: Plain_strain_compression
  // Load_Type = 4: Cyclic_loading
  // Load_Type = 5: Stress_relaxation
  // Load_Type = 6: Uniaxial_tension

  if(myrank==0)
  {
    for(int Load_Type = 0; Load_Type<6; Load_Type++)
      constitutive_model_test(hommat, NULL, Load_Type);

    constitutive_model_test(hommat, NULL, -1); // F(t)

    Matrix<double> L;
    L.initialization(3,3,0.0);
    L(0.0) = -1.001;
    L(1,1) = L(2,2) = 0.4998;

    constitutive_model_test(hommat, L.m_pdata, -1); // user defined L

  }

  destroy_zatnode(znod,nln);
  destroy_zatelem(zele_s,nle_s);
  destroy_zatelem(zele_v,nle_v);
  destroy_matgeom(matgeom,np);
  destroy_hommat(hommat,nhommat);
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
