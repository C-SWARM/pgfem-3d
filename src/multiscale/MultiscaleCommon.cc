/* HEADER */
#include "pgfem3d/MultiscaleCommon.hpp"
#include "allocation.h"
#include "constitutive_model.h"
#include "elem3d.h"
#include "enumerations.h"
#include "gen_path.h"
#include "generate_dof_ids.h"
#include "homogen.h"
#include "in.h"
#include "incl.h"
#include "initialize_damage.h"
#include "interface_macro.h"
#include "mesh_load.h"
#include "material.h"
#include "PGFEM_par_matvec.h"
#include "Psparse_ApAi.h"
#include "read_input_file.h"
#include "set_fini_def.h"
#include "utils.h"
#include <stdlib.h>
#include <search.h>
#include <assert.h>
#include <stdio.h>

using namespace pgfem3d;
using namespace pgfem3d::net;

namespace pgfem3d {

static const int ndim = 3;

MultiscaleCommon::MultiscaleCommon()
{
  opts = PGFEM_calloc(PGFem3D_opt, 1);
  sol = NULL;

  set_default_options(opts);

  /* options /solver information */
  SOLVER = NULL;
  Ap = NULL;
  Ai = NULL;

  nbndel = 0;
  bndel = NULL;
  ndofd = 0;
  DomDof = NULL;
  GDof = 0;

  /* mesh info */
  nn = 0;
  ne = 0;
  nce = 0;
  ndofn = 0;
  npres = 0;
  VVolume = 0.0;
  node = NULL;
  elem = NULL;
  coel = NULL;
  n_orient = 0;
  matgeom = NULL;
  nhommat = 0;
  hommat = NULL;
  ensight = NULL;
  supports = NULL;
  n_co_props = 0;
  co_props = NULL;

  /* mixed tangents */
  K_01 = NULL;
  K_10 = NULL;

  /* solution buffers */
  solution_buffer = NULL;
}

MultiscaleCommon::~MultiscaleCommon()
{
  for(int i = 0, e = idx_map.size;
      i < e; i++){
    destroy_MULTISCALE_SOLUTION(sol+i,
                this,
                opts->analysis_type);
  }
  free(sol);
  free(opts);
  destroy_common();
  sol_idx_map_destroy(&(idx_map));
}

void MultiscaleCommon::build_solutions(const int n_solutions)
{
  /*=== BUILD SOLUTIONS ===*/
  sol_idx_map_build(&idx_map, n_solutions);
  sol = PGFEM_calloc(MULTISCALE_SOLUTION, n_solutions);
  for(int i=0; i<n_solutions; i++){
    initialize_MULTISCALE_SOLUTION(sol + i);
    build_MULTISCALE_SOLUTION(sol +i, this, opts->analysis_type);
  }
}

void MultiscaleCommon::initialize(const int argc,
                  char **argv,
                  const CommunicationStructure *com,
                  const int mp_id,
          Multiphysics& mp)
{
  int myrank = com->rank;
  /* parse the command-line-style options */
  parse_command_line(argc, argv, myrank, opts);

  /* error check the command line */
  if (opts->solverpackage != HYPRE && opts->solverpackage != TRILINOS){
    if(myrank == 0)
      PGFEM_printerr("ERROR: only the HYPRE and TRILINOS solvers are supported!"
                     "%s:%s:%d\n",__func__,__FILE__,__LINE__);
    PGFEM_Comm_code_abort(com, 0);
  }

  switch (opts->analysis_type) {
  case DISP: break;
  case CM:
    if (opts->cm == DISP) break;
    /* deliberate drop through */
  default:
    if(myrank == 0)
      PGFEM_printerr("ERROR: only DISP analysis is supported!"
                     "%s:%s:%d\n",__func__,__FILE__,__LINE__);
    PGFEM_Comm_code_abort(com, 0);
    break;
  }

  /* attempt to write to output directory */
  if (make_path(opts->opath,DIR_MODE) != 0){
    if (myrank == 0){
      PGFEM_printerr("ERROR: could not write to %s!\n",
                     opts->opath);
    }
    PGFEM_Comm_code_abort(com, 0);
  }

  /*=== BUILD COMMON ===*/
  build_common(com, mp_id,mp);
  supports->multi_scale = opts->multi_scale;
}

void MultiscaleCommon::build_common(const CommunicationStructure *com,
                    const int mp_id,
            Multiphysics& mp)
{
  // save initial communication properties from parent com
  rank = com->rank;    // my rank
  nproc = com->nproc;  // number of processes
  boot = com->boot;    // boot handle
  net = com->net;      // network handle
  comm = com->comm;    // communicator
  hints = com->hints;  // communication hints

  /* initialize the solver information */
  SOLVER = nullptr;
  maxit_nl = 5;

  switch(opts->vis_format){
   case VIS_ENSIGHT: case VIS_VTK:
    ensight = new Ensight{};
    break;
   default: break;
  }

  /*=== READ MULTISCALE INPUT FILES ===*/
  char *in_fname = PGFEM_calloc(char, 500);
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
    Material *mater = NULL;

    int fv_ndofn = ndim;

    switch(opts->analysis_type)
    {
     case STABILIZED: //intented to over flow.
     case MINI:
     case MINI_3F:
      fv_ndofn = 4;
      break;
    }

    int physicsno = 1; // currently support only mechanical part
    err = read_input_file(opts,this,&nn, &Gnn,&ndofn,
                          &ne, &ni,&lin_err,
                          &lim_zero,&nmat,&n_con,
                          &n_orient,&node,
                          &elem,&mater,&matgeom,
                          &supports,&nln,&znod,&nle_s,&zele_s,
                          &nle_v,&zele_v, &fv_ndofn,physicsno,&ndim, mp.physicsname);

    /* error reading file(s) */
    if(err){
      PGFEM_printerr("[%d]ERROR: incorrectly formatted input file!\n",
                     rank);
      PGFEM_Comm_code_abort(com, 0);
    }

    /* error in input */
    if(nln > 0 || nle_s > 0 || nle_v > 0){
      PGFEM_printerr("[%d]ERROR: applied loads are not consistent w/"
                     " multiscale analysis!\n",rank);
      PGFEM_Comm_code_abort(com, 0);
    }

    destroy_zatnode(znod,nln);
    destroy_zatelem(zele_s,nle_s);
    destroy_zatelem(zele_v,nle_v);

    /* override prescribed displacements */
    if(opts->override_pre_disp){
      auto fn = opts->pre_disp_file;
      if (override_prescribed_displacements(*supports, fn) != 0) {
        PGFEM_printerr("[%d]ERROR: an error was encountered when"
                       " reading the displacement override file.\n"
                       "Be sure that there are enough prescribed"
                       " displacements in the file.\n",rank);
        PGFEM_Abort();
      }
    }

    /* read microscale normal/thickness */
    if(opts->multi_scale){
      supports->multi_scale = opts->multi_scale;
      err = read_interface_macro_normal_lc(opts->ipath,supports);
      if(err != 0){
        PGFEM_printerr("[%d] ERROR: could not read normal from file!\n"
                       "Check that the file \"%s/normal.in\""
                       " exists and try again.\n",
                       rank,opts->ipath);
        PGFEM_Comm_code_abort(com, 0);
      }
    }

    /* homogenized material properties */
    Mat_3D_orthotropic (nmat,mater,opts->analysis_type);
    long ***a = aloc3l (nmat,nmat,n_con);
    nhommat = list(a,ne,nmat,n_con,elem);
    hommat= build_hommat(nhommat);
    hom_matrices(a,ne,nmat,n_con,elem,mater,matgeom,
                 hommat,matgeom->SH,opts->analysis_type);

    dealoc3l(a,nmat,nmat);
    free(mater);

    if (opts->analysis_type == CM) {
      char *cm_filename = NULL;
      alloc_sprintf(&cm_filename,"%s/model_params.in",opts->ipath);
      FILE *cm_in = PGFEM_fopen(cm_filename, "r");
      read_model_parameters_list(nhommat,
                                 hommat, cm_in);
      free(cm_filename);
      fclose(cm_in);
    }
  }/* end reading *.in */

  /* read cohesive stuff */
  if (opts->cohesive){
    char *filename = PGFEM_calloc(char, 500);
    sprintf(filename,"%s%d.in.co_props",in_fname,rank);
    FILE *in1 = PGFEM_fopen(filename,"r");
    read_cohesive_properties(in1,&n_co_props,
                             &co_props,com);
    PGFEM_fclose(in1);

    /* read coheisve elements */
    sprintf (filename,"%s%d.in.co",in_fname,rank);
    in1 = PGFEM_fopen(filename,"r");

    /* temporary leftovers from old file format */
    long ncom = 0;
    CHECK_SCANF (in1,"%ld\n",&ncom);
    double **comat = aloc2 (ncom,4);

    /* read the cohesive element info */
    coel = read_cohe_elem (in1,ncom,ndim,nn,node,
               &nce,comat,ensight,
               opts->vis_format,rank,
               co_props);
    dealoc2 (comat,ncom);
    PGFEM_fclose (in1);
    free(filename);
  }

  /*=== END FILE READING ===*/
  list_el_prescribed_def(supports,node,elem,
                         NULL,ne,0,nn);

  bndel = list_boundary_el(ne,elem,nn,
               node,rank,&nbndel);

  DomDof = PGFEM_calloc(long, nproc);

  ndofd = generate_local_dof_ids(ne,nce,nn,
                 ndofn,node,
                 elem,coel,NULL,
                 this,mp_id);

  DomDof[rank] =
  generate_global_dof_ids(ne,nce,nn,
                          ndofn,node,elem,
                          coel,NULL,this,mp_id);

  net->allgather(NET_IN_PLACE,1,NET_DT_LONG,DomDof,
         1,NET_DT_LONG,comm);

  renumber_global_dof_ids(ne,nce,0,nn,
                          ndofn,DomDof,node,
                          elem,coel,NULL,this,mp_id);

  long NBN = distribute_global_dof_ids(ne,nce,
                                       0,nn,
                                       ndofn,ndim,
                                       node,elem,
                                       coel,NULL,this,mp_id);

  /* global stiffness pattern and communication structure */
  Ap = PGFEM_calloc(int, DomDof[rank]+1);

  // now create a new sparse structure
  spc = pgfem3d::SparseComm::Create(net, comm);

  Ai =  Psparse_ApAi(ne,0,nn, ndofn,ndofd,elem,
            NULL,node,nce, coel, this,
            opts->cohesive,mp_id);

  // initialize SparseComm buffers after pattern determined
  spc->initialize();
  spc->build_fast_maps(ndofd, DomDof[rank], GDof);

  SOLVER = pgfem3d::solvers::SparseSystem::Create(*opts,
                          comm,
                          Ap,
                          Ai,
                          DomDof,
                          opts->maxit,
                          lin_err);

  if (!supports->multi_scale) {
    VVolume = T_VOLUME (ne,ndim,elem,node);
    net->allreduce(NET_IN_PLACE, &VVolume, 1, NET_DT_DOUBLE,
           NET_OP_SUM, comm);
  } else {
    VVolume = supports->v0;
  }

  /* allocate solution_buffers */
  build_MULTISCALE_SOLUTION_BUFFERS(&solution_buffer,
                                    ndofd, DomDof[rank]);

  free(in_fname);

  /* compute/print summary information */
  long mesh_info[7];
  mesh_info[0] = Gnn;
  mesh_info[1] = NBN;
  mesh_info[2] = ne;
  mesh_info[3] = nce;
  mesh_info[4] = nbndel;
  mesh_info[5] = DomDof[rank];
  mesh_info[6] = Ap[DomDof[rank]];
  net->allreduce(NET_IN_PLACE, mesh_info+2, 5, NET_DT_LONG, NET_OP_SUM, comm);

  if (rank == 0) {
    print_interpreted_options(opts);
    if (opts->multi_scale) {
      if (supports->npd >= 9){
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
    PGFEM_printf ("Volume: %f\n\n",VVolume);
  }
}

void MultiscaleCommon::destroy_common()
{
  int mp_id = 0; // id of mutiphysics. Supported only Mechanical part (=0)
  delete SOLVER;
  delete spc;
  free(Ap);
  free(Ai);
  free(bndel);
  free(DomDof);
  destroy_node_multi_physics(nn,node,mp_id);
  destroy_elem(elem,ne);
  destroy_coel(coel,nce);
  destroy_matgeom(matgeom,n_orient);
  destroy_hommat(hommat,nhommat);
  delete ensight;
  destroy_supp(supports);
  destroy_cohesive_props(n_co_props,co_props);
  destroy_PGFEM_par_matrix((PGFEM_par_matrix *) K_01);
  destroy_PGFEM_par_matrix((PGFEM_par_matrix *) K_10);
  destroy_MULTISCALE_SOLUTION_BUFFERS(solution_buffer);
  free(solution_buffer);
}

} // end namespace pgfem3d
