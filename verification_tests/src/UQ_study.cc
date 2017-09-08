#include "allocation.h"
#include "homogen.h"
#include "utils.h"

#include "read_input_file.h"
#include "post_processing.h"

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

int change_material_properties(int argc, char **argv, char *filename_out)
{
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  int myrank = 0;
  int nprocs = 0;
  MPI_Comm_rank(mpi_comm, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if(myrank==0)
  {
    FILE *fp_mat = fopen("material_properties.in", "r");
    int matno = 0;
    char line[1024];
    fgets(line, 1024, fp_mat);

    sscanf(line, "%d", &matno);

    printf("material properties to be read: %d\n", matno);

    double *mat = (double *) malloc(matno*sizeof(double));
    for(int a=0; a<matno; a++)
    {
      fgets(line, 1024, fp_mat);
      sscanf(line, "%lf", mat+a);
      printf("%d: %e\n", a, mat[a]);
    }
    fclose(fp_mat);

    PGFem3D_opt options;

    if (argc <= 2)
    {
      if(myrank == 0){
        print_usage(stdout);
      }
      exit(0);
    }

    set_default_options(&options);
    re_parse_command_line(myrank,2,argc,argv,&options);
    char filename[1024];
    char filename_symbol[1024];
    char filename_header[1024];
    char system_cmd[32768];
    char system_cmd_meshing[2048];

    copy_filename(options.ifname, filename);
    sprintf(filename_symbol,"%s/../%s.out.header.symbol",options.ipath,filename);
    sprintf(filename_header,"%s/../%s.out.header",       options.ipath,filename);
    sprintf(filename_out,"%s/%s_macro.out.%d",options.opath,options.ofname,0);


    // read center positions of grains and store then is X
    FILE *fp_grains = fopen("center_postion_of_grains.in", "r");
    int grainno = 0;
    int layerno = 0;
    fgets(line, 1024, fp_grains); //skip comment
    fgets(line, 1024, fp_grains); //skip comment
    fgets(line, 1024, fp_grains); //skip comment
    fgets(line, 1024, fp_grains); //skip comment
    fgets(line, 1024, fp_grains); //skip comment
    fscanf(fp_grains, "%d%d", &grainno,&layerno);
    printf("Total number of grains: %d\n", grainno);
    printf("Total number of layers: %d\n", layerno);

    double *X = (double *) malloc(sizeof(double)*3*grainno);
    int *ngl = (int *) malloc(sizeof(int)*layerno);
    for(int a=0; a<layerno; a++)
    {
      fscanf(fp_grains, "%d", ngl+a);
      printf("%d ", ngl[a]);
    }
    printf("\n");
    for(int a=0; a<grainno; a++)
    {
      fscanf(fp_grains, "%lf%lf%lf", X+a*3+0, X+a*3+1, X+a*3+2);
      //printf("%e %e %e\n", X[a*3+0], X[a*3+1], X[a*3+2]);
    }

    fclose(fp_grains);

    // below is just for this current developing version. Later, these need to be updated in order for changing
    // material properties properly using correlation fuctions

    sprintf(system_cmd, "sed -e \"s|MAT_VALUES|");

    double E0[2] = {100.0e-3,211.0e-3};
    double sigma[2] = {16.0e-3,13.0e-3};

    int cnt = 0;
    int lid = 0;
    int mat_id = 0; // if mat_id==0: Al
                    // if mat_id==1: Ni

    for(int a=0; a<grainno; a++)
    {
      if(a==ngl[lid]+cnt)
      {
        cnt += ngl[lid];
        lid++;
        printf("switch materials: %d -> ", mat_id);
        mat_id = (mat_id == 0);
        printf("%d\n", mat_id);
      }

      //use x,y,z coordinate for the functions
      double x = X[a*3+0];
      double y = X[a*3+1];
      double z = X[a*3+2];

      double xi = mat[mat_id];
      double E = E0[mat_id] + sigma[mat_id]*sin(2*x)*cos(2*y)*sin(z)*xi;

      double nu = 0.25;
      double mu = E/2.0/(1.0+nu);
      double C01 = 0.5*mu;
      // double C10 = 0.0;

      sprintf(system_cmd, "%s %e %e 0.0 %e 0.0 0.0 0.25 0.0 0.0 1.0 1.0 1.0 1.7e+03 1 2\\n", system_cmd, E, C01, mu);
    }

    sprintf(system_cmd, "%s|\" %s > %s", system_cmd, filename_symbol, filename_header);
    //printf("[%s] is performed\n",  system_cmd);

    system(system_cmd);

    sprintf(system_cmd_meshing, "cd %s/..; ls; ./gen_meshes.pl -f %s -np %d", options.ipath, filename, nprocs);
    system(system_cmd_meshing);
    free(mat);
    free(X);
    free(ngl);
  }
  int G_temp = 0;
  int n_temp = myrank;
  MPI_Allreduce(&G_temp,&n_temp,1,MPI_INT,MPI_SUM,mpi_comm);

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

  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  int myrank = 0;
  // int nprocs = 0;
  MPI_Comm_rank(mpi_comm, &myrank);

  char filename_out[1024];
  err += change_material_properties(argc,argv,filename_out);
  err += single_scale_main(argc,argv_4func);

  if(myrank==0)
  {
    FILE *fp_in = fopen(filename_out, "r");
    char line[1024];
    double stress[6];
    fgets(line, 1024, fp_in);
    fgets(line, 1024, fp_in);
    fgets(line, 1024, fp_in);
    fgets(line, 1024, fp_in);
    fgets(line, 1024, fp_in);

    sscanf(line, "%lf %lf %lf %lf %lf %lf",
                  stress+0, stress+1, stress+2,
                  stress+3, stress+4, stress+5);
    fclose(fp_in);

    FILE *fp_out = fopen("stress.out", "w");
    fprintf(fp_out, "%e\n", stress[5]);
    fclose(fp_out);
  }

  return err;
}
