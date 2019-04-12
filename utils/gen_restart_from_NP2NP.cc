#ifdef HAVE_CONFIG_H
# include "config.h"
#endif


#include "allocation.h"
#include "homogen.h"
#include "post_processing.h"
#include "read_input_file.h"
#include "utils.h"
#include "pgfem3d/Communication.hpp"
#include "constitutive_model.h"
#include "restart.h"
#include "elem3d.h"

using namespace pgfem3d;
using namespace pgfem3d::net;

// define from t3d2fsifel  
#define ELM_ID   0 
#define ELM_Type 1
#define ELM_Mat  2

#define NODE_ID  0
  
template<class T>  
T max(Matrix<T> A,
      const int column_id){
  
  T max_num = -1e15;
  for(long ia=0; ia<A.m_row; ++ia)
    if(max_num<A(ia, column_id))
      max_num = A(ia, column_id);
  
  return max_num;
}

void read_a_map(Matrix<long> &node,
                Matrix<long> &elem,
                const int myrank,
                const int processno,
                const PGFem3D_opt &options){
  long nodeno;
  long elemno;
  char filename[2048];
  sprintf(filename, "%s/map.%d/%s%d.map", options.ipath, processno, options.ifname, myrank);
  FILE *in = fopen(filename, "r");
  
  if(in==NULL){
    printf("Cannot open file [%s]\n", filename);
    return;
  }
    
  CHECK_SCANF(in, "%ld %ld", &nodeno, &elemno);
  
  node.initialization(nodeno, 1);
  elem.initialization(elemno, 3);
  
  long nid = {};
  for(long ia=0; ia<nodeno; ++ia){
    CHECK_SCANF(in, "%ld", &nid);
    node(ia) = nid;
  }
  
  long eid = {};
  int  mat_id = {}, etype = {};
  
  for(long ia=0; ia<elemno; ++ia){
    CHECK_SCANF(in, "%ld %d %d", &eid, &mat_id, &etype);
    elem(ia, 0) = eid;
    elem(ia, 1) = mat_id;
    elem(ia, 2) = etype;
  }  
  
  fclose(in);
}

class Locals2Global{
  public:
    Matrix<Matrix<long>> N, E;
    long nodeno;
    long elemno;
    long matno;
    int  np_from;
    
    Locals2Global(){
      nodeno  = 0;
      elemno  = 0;
      matno   = 0;
      np_from = 0;
    }

    // read decomposed global node and element IDs    
    void read_maps(const int np_from_in,
                   const PGFem3D_opt &options){
      np_from = np_from_in;
      N.initialization(np_from, 1);
      E.initialization(np_from, 1);

      for(int ia=0; ia<np_from; ++ia){
        read_a_map(N(ia), E(ia), ia, np_from, options);
        long max_num = max(N(ia), NODE_ID);
        if(nodeno < max_num)
          nodeno = max_num;
        
        max_num = max(E(ia), ELM_ID);
        
        if(elemno < max_num)
          elemno = max_num;
        
        max_num = max(E(ia), ELM_Mat);
        
        if(matno < max_num) // count material IDs
          matno = max_num;        
      }
  
      ++matno;
    }
};

void read_a_restart_file(Matrix<Matrix<double>> &gU_nm1,
                         Matrix<Matrix<double>> &gU_n,
                         Matrix<Matrix<Constitutive_model>> &gM,
                         Locals2Global &L2G,
                         const Multiphysics &mp,
                         int myrank,
                         const PGFem3D_opt &options){

  double tnm1[2] = {-1.0,-1.0};
  
  long nodeno = L2G.N(myrank).m_row;
  long elemno = L2G.E(myrank).m_row;
  
  Grid grid;
  grid.nn  = nodeno;
  grid.ne  = elemno;
  grid.nsd = 3;
  
  Matrix<Element> elem(elemno, 1);
  for(long ia=0; ia<elemno; ++ia)
    elem(ia).toe = L2G.E.m_pdata[myrank](ia, 2);
  
  grid.element = elem.m_pdata;
  
  //set fv
  Matrix<FieldVariables> fv(mp.physicsno, 1);
  Matrix<Matrix<double>> u_nm1, u_n;
  EPS *eps = NULL;
  
  for(int mp_id=0; mp_id<mp.physicsno; ++mp_id){
    int ndofn = mp.ndim[mp_id];
    u_nm1(mp_id).initialization(nodeno, ndofn, 0.0);
      u_n(mp_id).initialization(nodeno, ndofn, 0.0);
       fv(mp_id).u_nm1 = u_nm1(mp_id).m_pdata;
       fv(mp_id).u_n   =   u_n(mp_id).m_pdata;
    if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL){
      if(options.analysis_type == CM || options.analysis_type == CM3F){
        eps = (EPS *) malloc(sizeof(EPS)*elemno);
        for(long eid = 0; eid<elemno; ++eid)
          eps[eid].model = gM(L2G.E(myrank)(eid, ELM_ID)).m_pdata;

        fv(mp_id).eps;
      }
    }
  }
  
  TimeStepping ts;
  LoadingSteps *load = NULL;
  read_restart(&grid, fv.m_pdata, &ts, load, &options, mp, tnm1, myrank);
  
  if(NULL != eps)
    free(eps);

}

class GlobalRestartValues{
  public:
    Matrix<Matrix<double>> gU_nm1, gU_n;
    Matrix<Matrix<Constitutive_model>> gM;
    State_variables *statv_list;
    Matrix<Model_parameters> params;
    bool isMechanicalCMActive;

    GlobalRestartValues(){
      isMechanicalCMActive = false;
    }
    
    ~GlobalRestartValues(){
      if(NULL != statv_list)
        delete[] statv_list;
    }
    
    // initialization. Resize gU_x, and gM
    void initialization(Locals2Global &L2G,
                        const Multiphysics &mp,
                        const PGFem3D_opt &options){

      gU_nm1.initialization(mp.physicsno, 1);
        gU_n.initialization(mp.physicsno, 1);
                  
      for(int mp_id=0; mp_id<mp.physicsno; ++mp_id){
        int ndofn = mp.ndim[mp_id];
        gU_nm1(mp_id).initialization(L2G.nodeno, ndofn);
        gU_n(  mp_id).initialization(L2G.nodeno, ndofn);
        
        if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL){
          if(options.analysis_type == CM || options.analysis_type == CM3F){
            isMechanicalCMActive = true;
            gM.initialization(L2G.elemno, 1);
          }          
        }
      }
    }
    
    void read_restart_files(Matrix<int> &model_types,
                            const PGFem3D_opt &options){
      
    }
    
    void construct_model_parameters(Matrix<int> &model_types){
      params.initialization(model_types.m_row, 1);
      HOMMAT h = {};
      h.devPotFlag = 1;
      h.volPotFlag = 2;

      for(int ia=0; ia<model_types.m_row; ++ia){
        params(ia).initialization(&h, model_types(ia));
      }
    } 
    void construct_constitutive_models(Locals2Global &L2G,
                                       const Multiphysics &mp){
      if(isMechanicalCMActive){
        long n_state_variables = 0;
        for(int ia=0; ia<L2G.np_from; ++ia){
          for(long ib=0; ib<L2G.E(ia).m_row; ++ib){
            long nint = {};
            int_point(L2G.E(ia)(ib, ELM_Type), &nint);
            n_state_variables += nint;
          }
        }
        statv_list = new State_variables[n_state_variables];
        long svid = 0;
        for(int ia=0; ia<L2G.np_from; ++ia){
          for(long ib=0; ib<L2G.E(ia).m_row; ++ib){
            long eid = L2G.E(ia)(ib, ELM_ID) - 1; // element ID starts from 1
            int  mid = L2G.E(ia)(ib, ELM_Mat);
            long nint = {};
            int_point(L2G.E(ia)(ib, ELM_Type), &nint);
            gM(eid).initialization(nint, 1);
            for(int ic=0; ic<nint; ++ic){
              gM(eid)(ic).model_id = svid;
              gM(eid)(ic).vars_list = &statv_list;              
              gM(eid)(ic).initialization(&params(mid));              
              ++svid;
            }
          }
        }
      }
    }
};

int read_model_params(Matrix<int> &model_types,
                      const PGFem3D_opt &options){

  int err = 0;
  int num_entries = {};
  
  char cm_filename[2048];
  sprintf(cm_filename,"%s/model_params.in",options.ipath);
  FILE *in = PGFEM_fopen(cm_filename, "r");
            
  err += scan_for_valid_line(in);
  CHECK_SCANF(in, "%d", &num_entries);
  
  char line[2048];
  
  for(int ia=0; ia<num_entries; ia++){
    int mat_id = {};
    int type = {};
    err += scan_for_valid_line(in);
    if(feof(in)) break;
    
    CHECK_SCANF(in, "%d %d", &mat_id, &type);  

    if(mat_id>0)
      model_types(mat_id) = type;
    
    while(fgetc(in) != '}'){
      err += scan_for_valid_line(in);      
      fgets(line, 2048, in);
    }

    
    if(feof(in)) break;
  }
    
  fclose(in);
  
  return err;
}



void get_options(int argc, 
                 char *argv[],
                 int &np_from,
                 int &np_to,
                 PGFem3D_opt *options,
                 const int myrank){
  if (argc <= 4) {
    if (myrank == 0) {
      printf("gen_restart_from_NP2NP [np(from)] [np(to)] [PGFem3D options]\n"); 
      print_usage(stdout);
    }
    exit(0);
  }
  
  sscanf(argv[1], "%d", &np_from);
  sscanf(argv[2], "%d", &np_to);
  
  set_default_options(options);
  re_parse_command_line(myrank, 4, argc, argv, options);
}


                          
int main(int argc,char *argv[])
{  
  Boot *boot = new Boot();
  int myrank = boot->get_rank();
  int nproc = boot->get_nproc();

  // read options
  PGFem3D_opt options;
  int np_from, np_to;
  get_options(argc, argv, np_from, np_to, &options, myrank);  

  // read decomposed global node and element IDs
  Locals2Global L2G;
  L2G.read_maps(np_from, options);
  
  // read model parameters
  Matrix<int> model_types(L2G.matno, 1, 0);
  read_model_params(model_types, options);  
  
  for(int ia=0; ia<L2G.matno; ia++)
    printf("%d -> %d: models(%d}: %d %d\n", np_from, np_to, myrank, ia, model_types(ia));  
     
  Multiphysics mp;
  read_multiphysics_settings(mp,&options,myrank);
  
  GlobalRestartValues grv;
  grv.initialization(L2G, mp, options);
  grv.construct_model_parameters(model_types);
  grv.construct_constitutive_models(L2G, mp);
  
  // set constitutive models
//  Matrix<Matrix<State_variables>> sv;
//  Matrix<Matrix<Constitutive_model>> gM;
  
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
  

  
  Matrix<Matrix<double>> Un(mp.physicsno, 1), Unm1(mp.physicsno, 1);

//  if(options.analysis == CM || options.analysis == CM3F || options.analysis == TF)
    
  
  net->barrier(com->comm);
  printf("%ld, %ld %ld\n", L2G.nodeno, L2G.elemno, L2G.matno);
        
  net->barrier(com->comm);
  
  destruct_multiphysics(mp);
                      
  exit(0);


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
  if(0)
    printf("%d\n", (int) is4cm);
  
  FILE *f = NULL;
  if(myrank==0)
    f = fopen("error.txt", "w");

  for(int tim=0; tim<nt; tim++){
    if(print[tim]){
      double err[error_no] = {};
      sprintf(filename,"%s/VTK/STEP_%.6d/%s_%d_%d.vtu",options.opath,tim,options.ofname,myrank,tim);
      //read_VTK_file4TF(filename, u,Ph,Vh);

      //compute_L2_error(err, elem, ne, node, u, Ph, Vh, times[tim+1], com, &options, hommat, is4cm);

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
