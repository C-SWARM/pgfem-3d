/// A util generating restart files from a prior run of different decomposition.
/// Maps from t3d2psifel are used to map restart file from a mesh decomposition to
/// another mesh decomposition.
/// 
/// How to use: gen_restart_from_NP2NP [np(from)] [np(to)] [PGFem3D options]
///
/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN


#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "local2global_id_mapping.h"
#include "allocation.h"
#include "homogen.h"
#include "post_processing.h"
#include "read_input_file.h"
#include "utils.h"
#include "pgfem3d/Communication.hpp"
#include "constitutive_model.h"
#include "restart.h"
#include "elem3d.h"
#include "three_field_element.h"
#include "plasticity_model.h"

using namespace pgfem3d;
using namespace pgfem3d::net;
  
typedef struct {
  int e_id;
  int n_ip;
  int mat_id;
  Matrix<int> ip_ids;
} IP_ID_LIST;


// global restart values
class GlobalCOValues{
  public:
    Matrix<double> angles; // global orietation angles
    Matrix<IP_ID_LIST> ip_ids;

    GlobalCOValues(){
    }
    
    ~GlobalCOValues(){
    }
    
    /// initialization of members
    void initialization(Locals2Global &L2G_from){

      angles.initialization(L2G_from.elemno, 1);
      ip_ids.initialization(L2G_from.elemno, 1);

      long n_ips = 0;      
      for(int ia=0; ia<L2G_from.np; ++ia){
        for(int ib=0; ib<L2G_from.E(ia).m_row; ++ib){
          long nint = {};
          long g_eid = L2G_from.E(ia)(ib, ELM_ID);
          int_point(L2G_from.E(ia)(ib, ELM_Type), &nint);
          ip_ids(g_eid).ip_ids.initialization(nint, 1);
          ip_ids(g_eid).n_ip = nint;
          ip_ids(g_eid).e_id = g_eid;
          ip_ids(g_eid).mat_id = L2G_from.E(ia)(ib, ELM_Mat);
                    
          for(int ic=0; ic<nint; ++ic)
            ip_ids(g_eid).ip_ids(ic) = n_ips + ic;
          
          n_ips += nint;
        }
      }
    }
    
    void read_a_orientation_file(Locals2Global &L2G_from,
                                 const char *CO,
                                 const char *co,
                                 const int myrank){      
      char fn[FILE_NAME_SIZE];
      sprintf(fn, "%s/%s_%d.in", CO, co, myrank);
      
      int ne = L2G_from.E(myrank).m_row;
      Matrix<int> e_ids(ne, 1);
      
      Matrix<double> angles(ne, 3);
      //IP_ID_LIST *elm_ip_map = NULL;
      //char fn_in[1024];
      //plasticity_model_read_orientations(e_ids, angles, elm_ip_map, fn_in, myrank, ne);
      
    }
    
    void read_orientation_files(Locals2Global &L2G_from,
                                const char *CO,
                                const char *co){
                                  
      for(int ia=0; ia<L2G_from.np; ++ia)
        read_a_orientation_file(L2G_from, CO, co, ia);

    }     
};
                        
int main(int argc, char *argv[])
{ /* 
  Boot *boot = new Boot();
  int myrank = boot->get_rank();
  int nproc = boot->get_nproc();

  // read options
  PGFem3D_opt options;
  
  MapInfo mapF; // map info from
  MapInfo mapT; // map info to
  
  char rs_path[FILE_NAME_SIZE];
  
  try{
    get_options(argc, argv, mapF, mapT, rs_path, &options, myrank);
  }
  catch(const char *msg){
    if(myrank==0)
      std::cerr << msg;
    return 0;
  }
  
  // read Multiphysics info
  Multiphysics mp;
  read_multiphysics_settings(mp,&options,myrank);
  
  // read decomposed global node and element IDs
  Locals2Global L2G_from, L2G_to;  
  L2G_from.read_maps(mapF);
  L2G_to.read_maps(mapT);

  // object to save global restart values
  GlobalCOValues grv;
  grv.read_a_orientation();
  grv.initialization(L2G_from, mp, options);
    
  // read and build model parameters
  try{
    read_model_params(grv, L2G_from.matno, options);
  }catch(const char *msg){
    if(myrank==0)
      std::cerr << msg;
    return 0;
  }
      
  for(int ia=0; ia<L2G_from.matno; ia++){
    std::cout << mapF.np << " -> " << mapT.np << ": models(rank = " << myrank << "): "
              << grv.params(ia)->mat_id << " " << grv.params(ia)->type << endl;
  }

      
  grv.construct_constitutive_models(L2G_from, mp);        
  grv.read_restart_files(L2G_from, mp, options);

  const char *tmp = options.opath;
  options.opath = rs_path;

  grv.write_restart_files(L2G_to, mp, options);
  
  options.opath = tmp;


    
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
    
  
  net->barrier(com->comm);
  printf("%ld, %ld %ld\n", L2G_from.nodeno, L2G_from.elemno, L2G_from.matno);
        
  net->barrier(com->comm);
  
  destruct_multiphysics(mp);

  delete com;
  //delete net;
  delete boot;
  */
  return(0);
}
