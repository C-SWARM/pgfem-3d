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

#include <sstream>

using namespace pgfem3d;
using namespace pgfem3d::net;

// defined from t3d2fsifel
#define ELM_ID   0  // element ID
#define ELM_Type 1  // element Type
#define ELM_Mat  2  // element material

#define NODE_ID  0  // node ID

constexpr const int FILE_NAME_SIZE = 2048;

/// maximum function
template<class T>
T max(Matrix<T> A,
      const int column_id){

  T max_num = -1e15;
  for(long ia=0; ia<A.m_row; ++ia)
    if(max_num<A(ia, column_id))
      max_num = A(ia, column_id);

  return max_num;
}

/// mapping info of number of partitions and mapping filenames
class MapInfo{
  public:
    int np;
    char filenames[FILE_NAME_SIZE];

    MapInfo(){
      np=0;
    }

    void print(void){
      std::cout << "number of partitions: " << np << endl;
      std::cout << "filebase: " << filenames << endl;
    }
};

/// read a local to global map. This function will re-size
/// the node and elem to store list of nodes and elements
///
/// \param[out] node    list of local node IDs to global IDs
/// \param[out] elem    list of local element IDs to global IDs
/// \param[in]  minfo   mapping file info to read
/// \param[in]  myrank  current process rank (partition ID of the mesh)
void read_a_map(Matrix<long> &node,
                Matrix<long> &elem,
                const MapInfo &minfo,
                const int myrank){
  long nodeno;
  long elemno;
  char filename[FILE_NAME_SIZE+1024];
  sprintf(filename, "%s_%d.map", minfo.filenames, myrank);
  FILE *in = fopen(filename, "r");

  if(in==NULL){
    printf("Cannot open file [%s]\n", filename);
    return;
  }

  CHECK_SCANF(in, "%ld %ld", &nodeno, &elemno);

  // re-size list of IDs
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

/// Local to global mapping class
class Locals2Global{
  public:

    Matrix<Matrix<long>> N; // list of local to global node IDs.
                            // size: {np}x[nodeno x 1]
    Matrix<Matrix<long>> E; // list of local to global element IDs.
                            // size: {np}x[elemno x 3(ELEM_ID, ELEM_type, ELEM_mat)]
    long nodeno; // total number of nodes
    long elemno; // total number of elements
    long matno;  // total number of material
    int  np;     // number of processes (number partitions)

    Locals2Global(){
      nodeno  = 0;
      elemno  = 0;
      matno   = 0;
      np = 0;
    }

    /// get global element ID of local element ID(index) for a process(myrank)
    long eid(const int myrank,
             const int index){
      return E(myrank)(index, ELM_ID) - 1; // ID from partitioner starts from 0
    }

    /// get global node ID of local node ID(index) for a process(myrank)
    long nid(const int myrank,
             const int index){
      return N(myrank)(index, NODE_ID) - 1; // ID from partitioner starts from 0
    }


    // read decomposed global node and element IDs
    void read_maps(const MapInfo minfo){
      np = minfo.np;
      // re-size members to number of partitions
      N.initialization(np, 1);
      E.initialization(np, 1);

      // read maps for each decomposed mesh
      for(int ia=0; ia<minfo.np; ++ia){
        // read a map for (ia)th map
        read_a_map(N(ia), E(ia), minfo, ia);

        // get total number of nodes
        long max_num = max(N(ia), NODE_ID);
        if(nodeno < max_num)
          nodeno = max_num;

        // get total number of elements
        max_num = max(E(ia), ELM_ID);
        if(elemno < max_num)
          elemno = max_num;

        // get total number of materials
        max_num = max(E(ia), ELM_Mat);

        if(matno < max_num) // count material IDs
          matno = max_num;
      }

      ++matno; // material ID starts from 0 such that number of materials = mat_id + 1
    }
};

// global restart values
class GlobalRestartValues{
  public:
    Matrix<Matrix<double>> gU_nm1, gU_n; // global nodal values for each physics
                                         // at t(n-1) and t(n)
    Matrix<double> gV_nm1, gV_n;         // global 3f volumes at t(n-1) and t(n)
    Matrix<double> gP_nm1, gP_n;         // global 3f pressures at t(n-1) and t(n)

    Matrix<Matrix<Constitutive_model>> gM; // global constitutive models CM or CM3F
    State_variables *statv_list;           // global state variables for CM or CM3F
    Matrix<Model_parameters *> params;     // list of model parameters for CM or CM3F
    bool isMechanicalCMActive;             // identify use of CM or CM3F
    Matrix<double> tns;  // t(n) for each physics
    double tnm1[3];      // t(n-1), t(n), t(n+1) at the sychronization
    Matrix<double> norm; // norm of residuals saved for each physics

    GlobalRestartValues(){
      isMechanicalCMActive = false;
      statv_list = NULL;
    }

    ~GlobalRestartValues(){
      if(statv_list != NULL)
        delete[] statv_list;
      if(isMechanicalCMActive && params.m_row*params.m_col > 0){
        for(int ia=0; ia<params.m_row*params.m_col; ++ia)
          params(ia)->finalization();
      }
    }

    /// initialization of members. re-size to number of physics
    void initialization(Locals2Global &L2G_from,
                        const Multiphysics &mp,
                        const PGFem3D_opt &options){

      // re-size to number of physics
      gU_nm1.initialization(mp.physicsno, 1); // nodal value at t(n-1)
        gU_n.initialization(mp.physicsno, 1); // nodal value at t(n)

         tns.initialization(mp.physicsno, 1); // t(n)
        norm.initialization(mp.physicsno, 1); // norm of residuals

      // re-size variables for each physics
      for(int mp_id=0; mp_id<mp.physicsno; ++mp_id){
        int ndofn = mp.ndim[mp_id]; // degree of freedom of a physics
        // set nodal values. size : nodeno x ndofn
        gU_nm1(mp_id).initialization(L2G_from.nodeno, ndofn);
        gU_n(  mp_id).initialization(L2G_from.nodeno, ndofn);

        // set only for momentum equation
        if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL){
          // build constitutive models if CM or CM3F is used
          if(options.analysis_type == CM || options.analysis_type == CM3F){
            isMechanicalCMActive = true;
            gM.initialization(L2G_from.elemno, 1);
          }
          // build 3f members if TF or CM3F is used
          if(options.analysis_type == TF || options.analysis_type == CM3F){
            long npres = {}, nVol = {};
            get_3f_pressure_volume_number(&npres, &nVol, options, 0);
            gV_nm1.initialization(L2G_from.elemno, nVol);
              gV_n.initialization(L2G_from.elemno, nVol);
            gP_nm1.initialization(L2G_from.elemno, npres);
              gP_n.initialization(L2G_from.elemno, npres);
          }
        }
      }
    }

    /// read restart files (n number of files for individual physics) for a partitioned mesh (myrank)
    ///
    /// \param[in] L2G_from local to global map
    /// \param[in] mp       mutiphysics information object
    /// \param[in] myrank   partion ID (current process rank)
    /// \param[in] options  PGFem3D option object
    void read_a_restart_file(Locals2Global &L2G_from,
                             const Multiphysics &mp,
                             int myrank,
                             const PGFem3D_opt &options){

      // set dummy PGFem3D object in order to used reading restart file function in PGFem3D
      //--------------------->
      // mesh object

      // materials
      Matrix<HOMMAT> hmat(L2G_from.matno, 1);
      Matrix<Matrix<long>> mat(L2G_from.matno, 1);
      for(int ia=0; ia<L2G_from.matno; ++ia){
        mat(ia).initialization(3,1);
        mat(ia)(0) = ia;

        hmat(ia).devPotFlag = 1;
        hmat(ia).volPotFlag = 2;
        hmat(ia).param      = params.m_pdata[ia];
        hmat(ia).mat_id     = ia;
      }

      // grid
      long nodeno = L2G_from.N(myrank).m_row;
      long elemno = L2G_from.E(myrank).m_row;

      Grid grid;
      grid.nn  = nodeno;
      grid.ne  = elemno;
      grid.nsd = 3;

      Matrix<Element> elem(elemno, 1);
      for(long ia=0; ia<elemno; ++ia){
        elem(ia).toe = L2G_from.E(myrank)(ia, ELM_Type);
        elem(ia).mat = mat(L2G_from.E(myrank)(ia, ELM_Mat)).m_pdata;
      }

      grid.element = elem.m_pdata;

      //set fv
      Matrix<FieldVariables> fv(mp.physicsno, 1);
      Matrix<Matrix<double>> u_nm1(mp.physicsno, 1), u_n(mp.physicsno, 1); // temporal nodal values
      EPS *eps = NULL;
      if(isMechanicalCMActive){
        eps = (EPS *) malloc(sizeof(EPS)*elemno);
        for(long eid = 0; eid<elemno; ++eid)
          eps[eid].model = gM(L2G_from.eid(myrank, eid)).m_pdata;

      }
      for(int mp_id=0; mp_id<mp.physicsno; ++mp_id){
        int ndofn = mp.ndim[mp_id];
        u_nm1(mp_id).initialization(nodeno, ndofn, 0.0);
        u_n(  mp_id).initialization(nodeno, ndofn, 0.0);
        fv(mp_id).u_nm1 = u_nm1(mp_id).m_pdata;
        fv(mp_id).u_n   = u_n(  mp_id).m_pdata;

        if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL){
          fv(mp_id).eps = eps;

          // nothing will happen if there is no use of the crystal plasticity model
          plasticity_model_set_orientations(eps, elemno, elem.m_pdata, L2G_from.matno, hmat.m_pdata, myrank);

          // construct 3f variables
          if(options.analysis_type == TF || options.analysis_type == CM3F){
            get_3f_pressure_volume_number(&(fv(mp_id).npres), &(fv(mp_id).nVol), options, myrank);
            fv(mp_id).tf.construct(elemno, fv(mp_id).nVol, fv(mp_id).npres);
          }
        }
      }

      // set time stepping object
      Matrix<double> tns_tmp(mp.physicsno, 1);
      TimeStepping ts;
      ts.tns = tns_tmp.m_pdata;

      double tnm1_tmp[3] = {};

      LoadingSteps *load = NULL; // in reading restart, load is not used
      // set dummy objects
      //<---------------------

      // read restart files using PGFem3D function

      read_restart(&grid, fv.m_pdata, &ts, load, &options, mp, tnm1, myrank);

      // save read temporal values into global
      // NOTE: global constitutive models are alread set because their references are used.
      for(int mp_id=0; mp_id<mp.physicsno; ++mp_id){
        int ndofn = mp.ndim[mp_id];
        for(long ia=0; ia<nodeno; ++ia){
          for(int ib=0; ib<ndofn; ++ib){
            gU_nm1(mp_id)(L2G_from.nid(myrank, ia), ib) = u_nm1(mp_id).m_pdata[ia*ndofn + ib];
            gU_n(  mp_id)(L2G_from.nid(myrank, ia), ib) = u_n(  mp_id).m_pdata[ia*ndofn + ib];
          }
        }
        if(options.analysis_type == TF || options.analysis_type == CM3F){
          get_3f_pressure_volume_number(&(fv(mp_id).npres), &(fv(mp_id).nVol), options, myrank);
          for(long ia=0; ia<elemno; ++ia){
            for(int ib=0; ib<fv(mp_id).nVol; ++ib){
              gV_nm1(L2G_from.eid(myrank, ia), ib) = fv(mp_id).tf.V_nm1(ia, ib);
              gV_n(  L2G_from.eid(myrank, ia), ib) = fv(mp_id).tf.V_n(  ia, ib);
            }
            for(int ib=0; ib<fv(mp_id).npres; ++ib){
              gP_nm1(L2G_from.eid(myrank, ia), ib) = fv(mp_id).tf.P_nm1(ia, ib);
              gP_n(  L2G_from.eid(myrank, ia), ib) = fv(mp_id).tf.P_n(  ia, ib);
            }
          }
        }
      }

      if(myrank==0){
        for(int mp_id=0; mp_id<mp.physicsno; ++mp_id){
          tns(mp_id)  = tns_tmp(mp_id);
          norm(mp_id) = fv(mp_id).NORM;
        }
        tnm1[0] = tnm1_tmp[0];
        tnm1[1] = tnm1_tmp[1];
        tnm1[2] = tnm1_tmp[2];
      }

      if(eps != NULL)
        free(eps);

    }

    /// read all restart files
    ///
    /// \param[in] L2G_from local to global map
    /// \param[in] mp       mutiphysics information object
    /// \param[in] options  PGFem3D option object
    void read_restart_files(Locals2Global &L2G_from,
                            const Multiphysics &mp,
                            const PGFem3D_opt &options){

      for(int ia=0; ia<L2G_from.np; ++ia)
        read_a_restart_file(L2G_from, mp, ia, options);
    }

    /// write restart files (n number of files for individual physics) for a partitioned mesh (myrank)
    ///
    /// \param[in] L2G_from local to global map
    /// \param[in] mp       mutiphysics information object
    /// \param[in] myrank   partion ID (current process rank)
    /// \param[in] options  PGFem3D option object
    void write_a_restart_file(Locals2Global &L2G_from,
                             const Multiphysics &mp,
                             int myrank,
                             const PGFem3D_opt &options){
      // set dummy PGFem3D object in order to used reading restart file function in PGFem3D
      //--------------------->
      // mesh object
      long nodeno = L2G_from.N(myrank).m_row;
      long elemno = L2G_from.E(myrank).m_row;

      Grid grid;
      grid.nn  = nodeno;
      grid.ne  = elemno;
      grid.nsd = 3;

      Matrix<Element> elem(elemno, 1);
      for(long ia=0; ia<elemno; ++ia)
        elem(ia).toe = L2G_from.E(myrank)(ia, ELM_Type);

      grid.element = elem.m_pdata;

      //set fv
      Matrix<FieldVariables> fv(mp.physicsno, 1);
      Matrix<Matrix<double>> u_nm1(mp.physicsno, 1), u_n(mp.physicsno, 1); // temporal nodal values
      EPS *eps = NULL;

      for(int mp_id=0; mp_id<mp.physicsno; ++mp_id)
        fv(mp_id).NORM = norm(mp_id);

      if(isMechanicalCMActive){
        eps = (EPS *) malloc(sizeof(EPS)*elemno);
        for(long eid = 0; eid<elemno; ++eid)
          eps[eid].model = gM(L2G_from.eid(myrank, eid)).m_pdata;
      }

      for(int mp_id=0; mp_id<mp.physicsno; ++mp_id){
        int ndofn = mp.ndim[mp_id];
        u_nm1(mp_id).initialization(nodeno, ndofn, 0.0);
        u_n(  mp_id).initialization(nodeno, ndofn, 0.0);

        for(long ia=0; ia<nodeno; ++ia){
          for(int ib=0; ib<ndofn; ++ib){
            u_nm1(mp_id).m_pdata[ia*ndofn + ib] = gU_nm1(mp_id)(L2G_from.nid(myrank, ia), ib);
            u_n(  mp_id).m_pdata[ia*ndofn + ib] = gU_n(  mp_id)(L2G_from.nid(myrank, ia), ib);
          }
        }

        fv(mp_id).u_nm1 = u_nm1(mp_id).m_pdata;
        fv(mp_id).u_n   = u_n(  mp_id).m_pdata;

        if(mp.physics_ids[mp_id] == MULTIPHYSICS_MECHANICAL){
          fv(mp_id).eps = eps;
          // construct 3f variables
          if(options.analysis_type == TF || options.analysis_type == CM3F){
            get_3f_pressure_volume_number(&(fv(mp_id).npres), &(fv(mp_id).nVol), options, myrank);
            fv(mp_id).tf.construct(elemno, fv(mp_id).nVol, fv(mp_id).npres);

            for(long ia=0; ia<elemno; ++ia){
              for(int ib=0; ib<fv(mp_id).nVol; ++ib){
                fv(mp_id).tf.V_nm1(ia, ib) = gV_nm1(L2G_from.eid(myrank, ia), ib);
                fv(mp_id).tf.V_n(  ia, ib) = gV_n(  L2G_from.eid(myrank, ia), ib);
              }
              for(int ib=0; ib<fv(mp_id).npres; ++ib){
                fv(mp_id).tf.P_nm1(ia, ib) = gP_nm1(L2G_from.eid(myrank, ia), ib);
                fv(mp_id).tf.P_n(  ia, ib) = gP_n(  L2G_from.eid(myrank, ia), ib);
              }
            }
          }
        }
      }

      LoadingSteps *load = NULL; // in reading restart, load is not used
      // set dummy objects
      //<---------------------

      // writerestart files using PGFem3D function
      write_restart(&grid,fv.m_pdata,load,&options,mp,tns.m_pdata,tnm1,myrank,options.restart);

      if(NULL != eps)
        free(eps);
    }

    /// write all restart files
    ///
    /// \param[in] L2G_from local to global map
    /// \param[in] mp       mutiphysics information object
    /// \param[in] options  PGFem3D option object
    void write_restart_files(Locals2Global &L2G_to,
                            const Multiphysics &mp,
                            const PGFem3D_opt &options){

      for(int ia=0; ia<L2G_to.np; ++ia)
        write_a_restart_file(L2G_to, mp, ia, options);
    }

    /// build list of constitutive models and their state variables
    /// for every element if momentum equation and CM or CM3F are used
    ///
    /// \param[in] L2G_from local to global map
    /// \param[in] mp       mutiphysics information object
    void construct_constitutive_models(Locals2Global &L2G_from,
                                       const Multiphysics &mp){

      // do only if momentum equation and CM or CM3F are used
      if(isMechanicalCMActive){
        // count number of state variables
        long n_state_variables = 0;
        for(int ia=0; ia<L2G_from.np; ++ia){
          for(long ib=0; ib<L2G_from.E(ia).m_row; ++ib){
            long nint = {};
            int_point(L2G_from.E(ia)(ib, ELM_Type), &nint);
            n_state_variables += nint;
          }
        }
        // allocate state variables as many as counted
        statv_list = new State_variables[n_state_variables];
        long svid = 0;

        // create constitutive models for an element
        // and assign state variables
        for(int ia=0; ia<L2G_from.np; ++ia){
          for(long ib=0; ib<L2G_from.E(ia).m_row; ++ib){
            long eid = L2G_from.eid(ia,ib);
            int  mid = L2G_from.E(ia)(ib, ELM_Mat);
            long nint = {}; // number of integration point
            int_point(L2G_from.E(ia)(ib, ELM_Type), &nint);
            // set size of constitutive model as many as nint
            gM(eid).initialization(nint, 1);
            for(int ic=0; ic<nint; ++ic){
              gM(eid)(ic).model_id = svid;
              gM(eid)(ic).vars_list = &statv_list;
              gM(eid)(ic).initialization(params(mid));
              ++svid;
            }
          }
        }
      }
    }
};

/// read model parameters from the model_param.in file
/// Reading restart file needs only model ID and constitutive model type.
/// This function skips all other constitutive model parameters.
///
/// \param[out] grv     object to save global restart values
/// \param[in]  matno   number of materials
/// \param[in]  options PGFem3D option object
void read_model_params(GlobalRestartValues &grv,
                       const int matno,
                       const PGFem3D_opt &options){
  if(grv.isMechanicalCMActive){

    // create Model_parameters
    grv.params.initialization(matno, 1);

    // temporal homat
    Matrix<HOMMAT> hmat(matno, 1);
    for(int ia=0; ia<matno; ++ia){
      hmat(ia).devPotFlag = 1;
      hmat(ia).volPotFlag = 2;
      hmat(ia).mat_id     = ia;
    }

    char cm_filename[FILE_NAME_SIZE];
    sprintf(cm_filename,"%s/model_params.in", options.ipath);
    FILE *in = PGFEM_fopen(cm_filename, "r");

    if(in == NULL){
      std::stringstream ss;
      ss << "Error: Cannot open [" << cm_filename << "]\n";
      throw ss.str();
    }

    if(read_model_parameters_list(matno, hmat.m_pdata, in)>0){
      fclose(in);
      std::stringstream ss;
      ss << "Error: Cannot read [" << cm_filename << "] properly\n";
      throw ss.str();
    }

    for(int ia=0; ia<matno; ++ia)
      grv.params.m_pdata[ia] = hmat(ia).param;

    fclose(in);
  }
}


/// get runtime options
/// The option [np_from] and [np_to] are expected first, and
/// PGFem3D run time options used in the simulation which generated restart
/// file is expected to be following.
///
/// \param[in]  argc    number of arguments
/// \param[in]  argv    argument values
/// \param[out] from    mapping local to global to be read
/// \param[out] to      mapping local to global to be written
/// \param[out] out     filebase of restart files to be written
/// \param[out] options PGFem3D option object
/// \param[in]  myrank  current process rank
void get_options(int argc,
                 char *argv[],
                 MapInfo &from,
                 MapInfo &to,
                 char *out,
                 PGFem3D_opt *options,
                 const int myrank){
  if (argc <= 3)
    throw "How to use: gen_restart_from_NP2NP [input_file_path] [PGFem3D options]\n";


  FILE *fp = fopen(argv[1], "r");

  if(fp==NULL){
    char err_msg[FILE_NAME_SIZE];
    sprintf(err_msg, "Error: Cannot open file [%s].\n", argv[1]);
    throw err_msg;
  }

  int err = 0;

  // map from
  err += scan_for_valid_line(fp);
  CHECK_SCANF(fp, "%d", &from.np);
  err += scan_for_valid_line(fp);
  CHECK_SCANF(fp, "%s", &from.filenames);

  // map to
  err += scan_for_valid_line(fp);
  CHECK_SCANF(fp, "%d", &to.np);
  err += scan_for_valid_line(fp);
  CHECK_SCANF(fp, "%s", &to.filenames);

  // new restart file base
  err += scan_for_valid_line(fp);
  CHECK_SCANF(fp, "%s", out);

  if(myrank==0){
    std::cout << "Map info from:" << endl;
    from.print();
    std::cout << "Map info to:" << endl;
    to.print();
    std::cout << "restart out: " << out << endl;
  }

  if(err>0){
    char err_msg[FILE_NAME_SIZE];
    sprintf(err_msg, "Error: Cannot read file [%s].\n", argv[1]);
    throw err_msg;
  }

  set_default_options(options);
  re_parse_command_line(myrank, 3, argc, argv, options);
}

int main(int argc, char *argv[])
{
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
  GlobalRestartValues grv;
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

  return(0);
}
