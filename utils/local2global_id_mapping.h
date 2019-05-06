#ifndef H__H__LOCAL2GLOBAL_ID_MAPPING__H__H
#define H__H__LOCAL2GLOBAL_ID_MAPPING__H__H
/// A util generating restart files from a prior run of different decomposition.
/// Maps from t3d2psifel are used to map restart file from a mesh decomposition to
/// another mesh decomposition.
/// 
/// How to use: gen_restart_from_NP2NP [np(from)] [np(to)] [PGFem3D options]
///
/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN


#include "data_structure.h"
#include "utils.h"

using namespace gcm;

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

#endif
