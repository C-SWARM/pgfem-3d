#include "allocation.h"
#include "homogen.h"
#include "utils.h"

#include "read_input_file.h"
#include "post_processing.h"
#include "PGFem3D_to_VTK.hpp"
#include <stdlib.h>
/*****************************************************/
/*           BEGIN OF THE COMPUTER CODE              */
/*****************************************************/
int single_scale_main(int argc, char **argv);

void copy_filename(char fn_from[], char f_to[]) 
{
  int c = 0;
 
  while (fn_from[c] != '\0') 
  {
    f_to[c] = fn_from[c];
    c++;
  }
  f_to[c-1] = '\0';
}   

int read_from_VTK(const PGFem3D_opt *opts, int myrank, int step, double *u)
{
  int err = 0;
  char filename[1024];
  sprintf(filename,"%s/VTK/STEP_%.5d/%s_%d_%d.vtu",opts->opath,step,opts->ofname,myrank, step);   
  err += read_VTK_file(filename, u);      
  return err;
}      

int change_material_properties(int argc, char **argv)
{
  /*=== END INITIALIZATION === */
  
  MPI_Comm mpi_comm = MPI_COMM_WORLD;    
  int myrank = 0;
  int nprocs = 0;
  MPI_Comm_rank(mpi_comm, &myrank); 
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs); 
  
  if(myrank==0)
  {  
    FILE *fp_mat = fopen("material_properties.in", "r");
    double E;
    fscanf(fp_mat, "%lf", &E);
    fclose(fp_mat);
  
    printf("%e is read\n", E);
  
    PGFem3D_opt options;

    if (argc <= 2)
    {
      if(myrank == 0){
        print_usage(stdout);
      }
      exit(0);
    }
    
    double nu = 0.25;
    double mu = E/2.0/(1.0+nu);
    double C01 = 0.5*mu;
    double C10 = 0.0;
    
    set_default_options(&options);
    re_parse_command_line(myrank,2,argc,argv,&options);

    char filename[1024];  
    char filename_symbol[1024];
    char filename_header[1024];    
    char system_cmd[2048];
    char system_cmd_meshing[2048];    
   
    copy_filename(options.ifname, filename);
      
    sprintf(filename_symbol,"%s/../%s.out.header.symbol",options.ipath,filename);
    sprintf(filename_header,"%s/../%s.out.header",       options.ipath,filename);    

    sprintf(system_cmd, "sed -e \"s|Exx|%e|\" -e \"s|C01|%e|\" -e \"s|C10|%e|\" -e \"s|mu|%e|\"  -e \"s|nu|%f|\" %s > %s", 
                                    E,               C01,             C10,             mu,              nu, filename_symbol, filename_header);
                                    
    system(system_cmd);
    sprintf(system_cmd_meshing, "./makeset.pl -np %d", nprocs);
    system(system_cmd_meshing);

  }  
  
  return 0;  
}

int main(int argc, char **argv)
{

  char **argv_4func;
  argv_4func = (char**)malloc((argc + 1) * sizeof(char*));
         
  for(long a = 0; a < argc; a++)
  {
    long n = strlen(argv[a]) + 1;
                argv_4func[a] = (char*)malloc(n);
                strcpy(argv_4func[a], argv[a]);
  }  
  
  int flag_MPI_Init;
  int err = MPI_Initialized(&flag_MPI_Init);
  if(!flag_MPI_Init)
    MPI_Init (&argc,&argv);  

  err += change_material_properties(argc,argv);
  err += single_scale_main(argc,argv_4func);
  return err;
}