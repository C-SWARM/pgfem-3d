/// A util generating crystal orientaion (CO)files from inputs for a different decomposition.
/// Maps from t3d2psifel are used to map CO files from a mesh decomposition to
/// another mesh decomposition.
/// 
/// How to use: gen_co_from_NP2NP [map info file] [CO_path/co_filebase]
///
/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN


#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "local2global_id_mapping.h"
#include "elem3d.h"
#include "gen_path.h"

using namespace pgfem3d;
using namespace pgfem3d::net;

enum CO_TYPE {
  CO_ELME_BY_ELEM, // set crystal orientation element by element
  CO_MAT_BY_MAT,   // set crystal orientation material by material
  CO_NO_CO,        // no crystal orientation
  CO_TYPE_end
};

// object for header of the crystal orientation file
typedef struct {
  char text[1024];
} HeaderText;

// global restart values
class GlobalCOValues{
  public:
    Matrix<Matrix<double>> angles; // global orietation angles
    Matrix<int> isSet;             // identify CO is set for element
    Matrix<HeaderText> COhead;     // header 
    Matrix<double> mat_angles;     // global CO set material by material
    bool mat_by_mat;               // if true, CO is set material by material
                                   // default is false

    GlobalCOValues(){
      mat_by_mat = false;
    }
    
    ~GlobalCOValues(){
    }
    
    /// initialization of members    
    void initialization(Locals2Global &L2G_from){
      
      // re-size members to total number of global elements
      angles.initialization(L2G_from.elemno, 1);
      isSet.initialization(L2G_from.elemno, 1, CO_NO_CO);
      mat_angles.initialization(L2G_from.matno, 1);
      
      // prepare each angle of an element as many as number of integration points 
      // looping all decomposition
      for(int ia=0; ia<L2G_from.np; ++ia){
        // looping all element in a decomposed mesh
        for(int ib=0; ib<L2G_from.E(ia).m_row; ++ib){
          long nint = {};
          long g_eid = L2G_from.eid(ia, ib);
          
          // get number of integration points
          int_point(L2G_from.E(ia)(ib, ELM_Type), &nint);
          // set angle size
          angles(g_eid).initialization(nint, 3); // [nint] x [3] matrix
        }
      }
    }
    
    /// read a CO file for a partitioned mesh (myrank)
    /// 
    /// \param[in] L2G_from local to global map
    /// \param[in] co_in    CO_path/co_filebase
    /// \param[in] myrank   partion ID (current process rank)
    void read_a_orientation_file(Locals2Global &L2G_from,
                                 const char *co_in,
                                 const int myrank){
      // set CO file name for the decomposition_myrank
      char fn[FILE_NAME_SIZE];
      sprintf(fn, "%s_%d.in", co_in, myrank);
      
      // try to open the CO file
      FILE *fp = fopen(fn, "r");
      if(fp==NULL)
      {
        std::cout << "fail to read [" << fn << "]" << endl;
        abort();
      }
      
      char line[1024];
      
      if(myrank>0){
        // read header from file_base_0.in file
        // when write orientation file, this header will be used for every rank.
        
        // count number of headers
        int headerno = 0;
        while(fgets(line, 1024, fp)!=NULL){
          if(line[0]=='#')
            ++headerno;
          else
            break;
        }
        rewind(fp); // go back to beginning
        
        // set CO headers by reading file_base_0.in file
        COhead.initialization(headerno, 1);
        
        for(int ia=0; ia<headerno; ia++)
          fgets(COhead(ia).text, 1024, fp); 

        rewind(fp); // go back to beginning in order to let file_base_0.in use
                    // same routines as other files [file_base_#.in]
      }      

      // read COs
      while(fgets(line, 1024, fp)!=NULL){
        // skip comments
        if(line[0]=='#')
          continue;
    
        
        int e, ip;
        double x1, x2, x3;
        sscanf(line, "%d %d %lf %lf %lf", &e, &ip, &x1, &x2, &x3);
        
        if(e<0){ 
          // orientation is set material by material
          mat_by_mat = true;
          mat_angles(ip, 0) = x1;
          mat_angles(ip, 1) = x2;
          mat_angles(ip, 2) = x3;
        }else{
          // orientation is set element by element
          long g_eid = L2G_from.eid(myrank, e);
          angles(g_eid)(ip, 0) = x1;
          angles(g_eid)(ip, 1) = x2;
          angles(g_eid)(ip, 2) = x3;
          isSet(g_eid) = CO_ELME_BY_ELEM;
        }
      }
      fclose(fp);      
    }
    
    /// read all CO files
    /// 
    /// \param[in] L2G_from local to global map
    /// \param[in] co_in    CO_path/co_filebase
    void read_orientation_files(Locals2Global &L2G_from,
                                const char *co_in){
                                  
      for(int ia=0; ia<L2G_from.np; ++ia)
        read_a_orientation_file(L2G_from, co_in, ia);

    }

    /// write a CO file for a partitioned mesh (myrank)
    /// 
    /// \param[in] L2G_from local to global map
    /// \param[in] co_out    CO_path_out/co_filebase
    /// \param[in] myrank   partion ID (current process rank)    
    void write_a_orientation_file(Locals2Global &L2G_to,
                                  const char *co_out,
                                  const int myrank){                                    
      char fn[FILE_NAME_SIZE];
      sprintf(fn, "%s_%d.in", co_out, myrank);
      
      FILE *fp = fopen(fn, "w");
      if(fp==NULL)
      {
        std::cout << "fail to create [" << fn << "]" << endl;
        abort();
      }
      
      // write heads first
      for(int ia=0; ia<COhead.m_row; ++ia)          
        fprintf(fp, "%s", COhead(ia).text);
      
      // write COs if set material by material
      if(mat_by_mat){  
        for(int mat_id=0; mat_id<mat_angles.m_row; ++mat_id){
          fprintf(fp, "-1 %d %e %e %e\n", mat_id, mat_angles(mat_id, 0),
                                                  mat_angles(mat_id, 1),
                                                  mat_angles(mat_id, 2));
        }
      }
      
      // looping all element and write angles if CO is set element by element
      for(int eid=0; eid<L2G_to.E(myrank).m_row; ++eid){
        long g_eid = L2G_to.eid(myrank, eid);
          
        if(isSet(g_eid)==CO_ELME_BY_ELEM){
          for(int ip=0; ip<angles(g_eid).m_row; ++ip)
            fprintf(fp, "%d %d %e %e %e\n", eid, ip, angles(g_eid)(ip, 0), 
                                                     angles(g_eid)(ip, 1), 
                                                     angles(g_eid)(ip, 2));
        }
      }
      
      fclose(fp);      
    }

    /// write all CO files. 
    /// co_out_dir/co_out_fb_#.in file will be ceated
    /// 
    /// \param[in] L2G_from   local to global map
    /// \param[in] co_out_dir CO out directory path, if not exist, 
    ///                       this folder will be created
    /// \param[in] co_out_fb  CO out filebase 
    void write_orientation_files(Locals2Global &L2G_to,
                                 const char *co_out_dir,
                                 const char *co_out_fb){
    // check default CO directory exists
    if(make_path(co_out_dir,DIR_MODE) != 0){
      std::cout << "Cannot create directory [" << co_out_dir << "]" << endl;
      abort();
    } 

    char co_out[FILE_NAME_SIZE];
    sprintf(co_out, "%s/%s", co_out_dir, co_out_fb);
    
    for(int ia=0; ia<L2G_to.np; ++ia)
      write_a_orientation_file(L2G_to, co_out, ia);
    }

};

/// get runtime options
/// The option [map info file] is expected first, and [CO_path/co_filebase]
/// is expected to be following.
///
/// \param[in]  argc       number of arguments
/// \param[in]  argv       argument values
/// \param[out] from       mapping local to global to be read
/// \param[out] to         mapping local to global to be written
/// \param[out] co_in      filebase of CO files to be read
/// \param[out] co_out_dir CO directory path to be written
/// \param[out] co_out_fb  filebase of CO files to be written
void get_options(int argc, 
                 char *argv[],
                 MapInfo &from,
                 MapInfo &to,
                 char *co_in,
                 char *co_out_dir,
                 char *co_out_fb){                  
  if (argc < 3){
    std::cout << "How to use: gen_co_from_NP2NP [input_file_path] [CO dir/filebase]" << endl;
    std::cout << "----------------------"      << endl;
    std::cout << "input_file_path format"      << endl;
    std::cout << "----------------------"      << endl;
    std::cout << "# : this is a comment"       << endl;
    std::cout << "# number of partitions from" << endl;
    std::cout << "5040"                        << endl;
    std::cout << "# filebases of map from"     << endl;
    std::cout << "./map.5040/case_1"           << endl;
    std::cout << "# number of partitions to"   << endl;
    std::cout << "192"                         << endl;
    std::cout << "# filebases of map to"       << endl;                     
    std::cout << "./map.192/case_1"            << endl;                     
    std::cout << "# orientation out directory" << endl;
    std::cout << "./CO.192"                    << endl;
    std::cout << "# orientation out filebase"  << endl;
    std::cout << "co"                          << endl;
    abort();    
  }

  
  FILE *fp = fopen(argv[1], "r");
  
  if(fp==NULL){
    std::cout << "Error: Cannot open file [" << argv[1] << "]." << endl;
    abort();
  }
  
  sprintf(co_in, argv[2]);
  
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
  
  // new CO directory path
  err += scan_for_valid_line(fp);
  CHECK_SCANF(fp, "%s", co_out_dir);
  
  // new CO file base 
  err += scan_for_valid_line(fp);
  CHECK_SCANF(fp, "%s", co_out_fb);  
  
  // print read options
  std::cout << "Map info from:" << endl;
  from.print();
  std::cout << "Map info to:" << endl;
  to.print();
  std::cout << "Read CO from: " << co_in << endl;
  std::cout << "Write CO to: " << co_out_dir << "/" << co_out_fb<< endl;    
    
  if(err>0){
    std::cout << "Error: Cannot read file [" << argv[1] << "]" << endl;
    abort();
  }    
}
                        
int main(int argc, char *argv[])
{

  MapInfo mapF; // map info from
  MapInfo mapT; // map info to
  char co_in[FILE_NAME_SIZE];
  char co_out_dir[FILE_NAME_SIZE];
  char co_out_fb[FILE_NAME_SIZE];
    
  get_options(argc, argv, mapF, mapT, co_in, co_out_dir, co_out_fb);
  
  // read decomposed global node and element IDs
  Locals2Global L2G_from, L2G_to;  
  L2G_from.read_maps(mapF);
  L2G_to.read_maps(mapT);

  // object to save global crystal orientation values
  GlobalCOValues gco;
  gco.initialization(L2G_from);
    
  // read orientation files
  gco.read_orientation_files(L2G_from, co_in);
  
  // write orientation files
  gco.write_orientation_files(L2G_to, co_out_dir, co_out_fb);
  
  printf("All files are successfully mapped.\n");
  return(0);
}
