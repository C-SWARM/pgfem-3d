/// Functions are defined for reading input files
/// 
/// Authors:
///   Matt Mosby, University of Notre Dame, <mmosby1 [at] nd.edu>
///   Karel Matous, University of Notre Dame, <kmatous [at] nd.edu>
///   Sangmin Lee, University of Notre Dame <slee43 [at] nd.edu>


/* HEADER */
#include "read_input_file.h"

#include "enumerations.h"
#include "allocation.h"
#include "in.h"
#include <string.h>
#include "PGFem3D_data_structure.h"
#include "utils.h"
#include "Arc_length.h"
#include "read_cryst_plast.h"
#include "load.h"
#include "constitutive_model.h"
#include "restart.h"
#include "gen_path.h"

/// read mechanical part of material properties
///
/// \param[in] fp file pointer for reading mechanical part of material properties
/// \param[in,out] mat material property object
/// \param[in] opts PGFem3D options
/// \return non-zero on interal error
int read_material_for_Mechanical(FILE *fp,
                                 MATERIAL_PROPERTY *mat,
                                 const PGFem3D_opt *opts)
{
  int err = 0;  
  for(int ia=0; ia<mat->nmat; ia++)
  {
    scan_for_valid_line(fp);
    if(read_material(fp,ia,mat->mater,opts->legacy))
      PGFEM_Abort();
  }  
  return err;
}

/// read mechanical part of material properties
///
/// \param[in] fp file pointer for reading mechanical part of material properties
/// \param[in,out] mat material property object
/// \param[in] opts PGFem3D options
/// \return non-zero on interal error
int read_material_for_Thermal(FILE *fp,
                              MATERIAL_PROPERTY *mat,
                              const PGFem3D_opt *opts)
{
  int err = 0;  
  int param_in = 10;
  
  MATERIAL_THERMAL *thermal = (MATERIAL_THERMAL *) malloc(sizeof(MATERIAL_THERMAL)*(mat->nmat));
  
  for(int ia=0; ia<mat->nmat; ia++)
  {
    double cp;
    double k[9];

    int match = 0;
    scan_for_valid_line(fp);
    match += fscanf(fp, "%lf", &cp);
     
    scan_for_valid_line(fp);
    match += fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", k+0, k+1, k+2
                                                             , k+3, k+4, k+5
                                                             , k+6, k+7, k+8);
    if(match != param_in)                                                         
      PGFEM_Abort();
      
    thermal[ia].cp = cp;
    for(int ib=0; ib<9; ib++)
      thermal[ia].k[ib] = k[ib];  
  }
  
  mat->thermal = thermal; 
  return err;
}
                                                                
/// read material properties for multiphysics problem
///
/// \param[in,out] mat material property object
/// \param[in] opts PGFem3D options
/// \param[in] mp multiphysics object
/// \return non-zero on interal error
int read_multiphysics_material_properties(MATERIAL_PROPERTY *mat,
                                          const PGFem3D_opt *opts,
                                          const MULTIPHYSICS *mp)
{
  int err = 0;

  int myrank = 0; 
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  char dirname[1024], fn[1024];
  sprintf(dirname,"%s/Material",opts->ipath);
      
  for(int ia=0; ia<mp->physicsno; ia++)
  {
    sprintf(fn,"%s/%s.mat",dirname,mp->physicsname[ia]);
    
    FILE *fp = NULL;
    fp = fopen(fn, "r");
    
    if(fp==NULL)
    { 
      if(myrank==0) 
        printf("No [%s] exists.\n", fn);

      continue;      
    }
    
    switch(mp->physics_ids[ia])
    {
      case MULTIPHYSICS_MECHANICAL:
        err += read_material_for_Mechanical(fp,mat,opts);
        break;
      case MULTIPHYSICS_THERMAL:
        err += read_material_for_Thermal(fp,mat,opts);        
        break;
      case MULTIPHYSICS_CHEMICAL:
        break;
      default:
        break;  
    }
    fclose(fp);        
  }
  
  // read and set general material properties e.g. density
  
  mat->density = (double *) malloc(sizeof(double)*(mat->nmat));
  double *d = mat->density;
  for(int ia=0; ia<mat->nmat; ia++)
    d[ia] = 0.0;  
    
  sprintf(fn,"%s/material.mat",dirname);
  FILE *fp = NULL;
  fp = fopen(fn, "r");
  
  if(fp == NULL)
  {
    if(myrank==0)
      printf("No [%s] exists. \nDensity is set to zero.\n", fn);
    
    return err;  
  }
  
  int match = 0;    
  for(int ia=0; ia<mat->nmat; ia++)
  {
    scan_for_valid_line(fp);
    match += fscanf(fp, "%lf", d+ia);
  }

  fclose(fp);   
  
  // check number of densities that is read.
  if(match != mat->nmat)
  {
    if(myrank==0)
      printf("Material density is not read as many as number of materials.\n");      

    PGFEM_Abort();
  }
  else
  {
    for(int ia = 0; ia<mat->nhommat; ia++)
    {
      (mat->hommat[ia]).density = d[(mat->hommat[ia]).mat_id];
      if(myrank==0)
        printf("Density(%d), %e\n", ia, (mat->hommat[ia]).density);
    }      
  }
  return err;
}                                          

/// count number of ranges that are seperated by comma
/// 
/// \param[in] str a string containing ranges
/// \return a integer number (counted number of ranges)
int count_number_of_ranges(char str[])
{
  int charno = strlen(str);
  int rangeno = 1; 

  for(int a=0; a<=charno; a++)
  {
    if(str[a]==',')
      rangeno++;
  }
  
  return rangeno;
}

/// get list of ranges
/// \param[in] str_in a string containing ranges
///      str_in format: 
///      ---------------------------------------------------------------
///      srt_in = "0:0.001:0.1"
///         or
///      str_in = "0:0.001:0.1,0.105:0.005:0.2"
///      ---------------------------------------------------------------
///      range format: start_value:step_size:end_value
///      ex) 1:1:10 => 1,2,3,4,5,6,7,8,9,10
///      ex) 0:0.2:1 => 0,0.2,0.4,0.6,0.8,1
/// \param[out] ranges double array of ranges read from the input string
/// \return non-zero on internal error
int interpret_ranges(double *ranges, char str_in[])
{
  int err = 0;
  int rangeno = count_number_of_ranges(str_in); 
  int charno = strlen(str_in);
    
  char str[1024], *pch;
  memcpy(str, str_in, charno+1);
  
  pch = strtok (str,","); 
  for(int a=0; a<rangeno; a++)
  {
    sscanf(pch, "%lf:%lf:%lf", ranges+a*3+0, ranges+a*3+1, ranges+a*3+2);
    pch = strtok(NULL, ",");
  }      
  
  return err;
}
                   
int read_input_file(const PGFem3D_opt *opts,
		    MPI_Comm comm,
		    long *nn,
		    long *Gnn,
		    long *ndofn,
		    long *ne,
		    long *lin_maxit,
		    double *lin_err,
		    double *lim_zero,
		    long *nmat,
		    long *n_concentrations,
		    long *n_orient,
		    NODE **node,
		    ELEMENT **elem,
		    MATERIAL **material,
		    MATGEOM *matgeom,
		    SUPP *sup,
		    long *nln,
		    ZATNODE **znod,
		    long *nel_s,
		    ZATELEM **zelem_s,
		    long *nel_v,
		    ZATELEM **zelem_v,
		    const int physicsno,
		    const int *ndim)
{
  int err = 0;
  int myrank = 0;
  MPI_Comm_rank(comm,&myrank);

  /* compute filename and open file */
  char *filename = PGFEM_calloc(500,sizeof(char));
  sprintf(filename,"%s/%s%d.in",opts->ipath,opts->ifname,myrank);
  FILE *in = PGFEM_fopen(filename,"r");

  /* read header lines */
  fscanf (in,"%ld %ld %ld",nn,ndofn,ne);
  fscanf (in,"%ld %lf %lf",lin_maxit,lin_err,lim_zero);
  fscanf (in,"%ld %ld %ld",nmat,n_concentrations,n_orient);

  /* Set ndofn according to analysis type */
  switch(opts->analysis_type){
  case STABILIZED: case MINI: case MINI_3F: *ndofn = 4; break;
  default: *ndofn = 3; break;
  }

  (*node) = build_node_multi_physics(*nn,ndim,physicsno);
  (*elem) = build_elem(in,*ne,opts->analysis_type);
  (*material) = PGFEM_calloc(*nmat,sizeof(MATERIAL));
  (*matgeom) = build_matgeom(*n_concentrations,*n_orient);

  *Gnn = read_nodes(in,*nn,*node,opts->legacy,comm);
  /* NOTE: Supports assume only ndim supported dofs per node! */
  
  for(int ia=0; ia<physicsno; ia++)
    sup[ia] = read_supports(in,*nn,ndim[ia],*node, ia);    
 
  read_elem(in,*ne,*elem,*sup,opts->legacy);

  for(int i=0, e=*nmat; i<e; i++){
    if ( read_material(in,i,*material,opts->legacy) ){
      PGFEM_Abort();
    }
  }
  
  if (feof(in)) {
    PGFEM_printerr("ERROR: prematurely reached EOF in %s(%s)\n",
                   __func__,__FILE__);
    PGFEM_Abort();
  }

  if ( override_material_properties(*nmat,opts,*material) ) {
    PGFEM_Abort();
  }

  read_matgeom(in,*n_concentrations,*n_orient,*matgeom);

  /* NOTE: Node/Element loading assumes forces only in ndim
     directions */
  /* node */
  fscanf(in,"%ld",nln);
  *znod = build_zatnode (*ndofn,*nln);
  read_nodal_load (in,*nln,*ndofn,*znod);
  /* surface */
  fscanf (in,"%ld",nel_s);
  *zelem_s = build_zatelem (*ndofn,*nel_s);
  read_elem_surface_load (in,*nel_s,*ndofn,*elem,*zelem_s);
  /* volume */
  fscanf (in,"%ld",nel_v);
  *zelem_v = build_zatelem (*ndofn,*nel_v);

  /* check the ferror bit */
  if(ferror(in)) err++;

  /* free local memory and close file */
  free(filename);
  fclose(in);
  return err;
}

/// This function read mesh info, boundary conditions, and material properties 
/// from main input files (*.in). While reading inputs, node, element, material, 
/// and support objects are constructed. This function still utilizes 
/// lagacy function (read_input_file). Later, this function and read_input_file need to
/// be combined.
///
/// \param[out] grid a mesh object
/// \param[out] mat a material object
/// \param[out] FV array of field variable object
/// \param[out] SOL array of solution scheme object
/// \param[out] load object for loading
/// \param[in] mp multiphysics object
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int read_mesh_file(GRID *grid, 
                   MATERIAL_PROPERTY *mat,
                   FIELD_VARIABLES *FV,
                   SOLVER_OPTIONS *SOL,
                   LOADING_STEPS *load,
                   MULTIPHYSICS *mp,
                   MPI_Comm mpi_comm,
                   const PGFem3D_opt *opts)
{
  long ndofn;
  int err = read_input_file(opts,
                            mpi_comm,
                            &(grid->nn),
                            &(grid->Gnn),
                            &ndofn,
                            &(grid->ne),
                            &(SOL[0].iter_max_sol),
                            &(SOL[0].err),
                            &(SOL[0].computer_zero),
                            &(mat->nmat),
                            &(FV[0].n_concentrations),
                            &(mat->n_orient),
                            &(grid->node),
                            &(grid->element),
                            &(mat->mater),
                            &(mat->matgeom),
                            load->sups,
                            &(load->nln),
                            &(load->znod),
                            &(load->nle_s),
                            &(load->zele_s),
                            &(load->nle_v),
                            &(load->zele_v),
                            mp->physicsno,
                            mp->ndim);
  // read multiphysics material properties
  err += read_multiphysics_material_properties(mat,opts,mp);

  // update numerical solution scheme parameters
  for(int iA=1; iA<mp->physicsno; iA++)
  {
    SOL[iA].iter_max_sol  = SOL[0].iter_max_sol;
    SOL[iA].err           = SOL[0].err; 
    SOL[iA].computer_zero = SOL[0].computer_zero;
    FV[iA].n_concentrations = FV[0].n_concentrations;
    if(mp->physics_ids[iA]==MULTIPHYSICS_MECHANICAL)
      FV[iA].ndofn = ndofn;   
  }                            
  // need to update number of elements that have prescribed BCs (supported)                             
  for (long ia=0;ia<grid->ne;ia++)
  {
    for(int iA = 1; iA<mp->physicsno; iA++) // iA = 0 is alreaded accounted in read_elem in read_input_file
    {
      const long *nod = grid->element[ia].nod;
      const int nne = grid->element[ia].toe;
      int is_it_supp = 0;
      for(long ja=0; ja<load->sups[iA]->ndn; ja++)
      {
        for (int ka=0; ka<nne; ka++)
        {
          if(load->sups[iA]->lnpd[ja] == nod[ka])
          {
            is_it_supp = 1;
            break;
          }
        }
        if(is_it_supp)
          break;        
      }
      if(is_it_supp)    
      (load->sups[iA]->nde)++;
    }
  }                        
  return err;
}                   

int read_time_steps(FILE *fp, PGFem3D_TIME_STEPPING *ts)
{
  int err = 0;
  
  // read number of computational times
  fscanf(fp,"%ld",&(ts->nt));
  
  // read times
  ts->times = aloc1(ts->nt+1);
  for (long ia=0; ia<ts->nt+1; ia++)
    fscanf(fp,"%lf",(ts->times)+ia);
  
  long n_p = 0;  
  // read times for output
  fscanf(fp,"%ld",&n_p);
  
  //Times for printing
  ts->print = times_print(fp,ts->nt,n_p);      
      
  return err;
} 

/// Read solver file for time stepping.
/// If command line includes override solver file option, all solver files will be overrided.
/// At the end of this function, file pointer is stored in LOADING_STEPS 
/// in order to read load increments as time elapses.
/// Detailed slover file format can be found at the following link:
/// https://wiki-cswarm.crc.nd.edu/foswiki/pub/Main/CodeDevelopment/PGFem3DQuickStarts/generate_input_file.pdf
///
/// \param[out] time_steps object for time stepping
/// \param[out] mat a material object
/// \param[out] FV array of field variable object
/// \param[out] SOL array of solution scheme object
/// \param[out] load object for loading
/// \param[out] arc an object for Arc length scheme
/// \param[out] crpl object for lagcy crystal plasticity
/// \param[in] mp multiphysics object
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_solver_file(PGFem3D_TIME_STEPPING *ts,
                     MATERIAL_PROPERTY *mat,
                     FIELD_VARIABLES *FV,
                     SOLVER_OPTIONS *SOL,
                     LOADING_STEPS *load,
                     ARC_LENGTH_VARIABLES *arc,                    
                     CRPL *crpl,
                     MULTIPHYSICS *mp,
                     const PGFem3D_opt *opts,
                     int myrank)
{
  int err = 0;
  // READ SOLVER FILE
  // override the default solver file with one specified
  // at commandline
  char filename[1024];
  char in_dat[1024];
  FILE *fp = NULL;
  if(opts->override_solver_file)
  {
    if(myrank == 0)
      printf("Overriding the default solver file with:\n%s\n", opts->solver_file);

    fp = fopen(opts->solver_file,"r");
  } 
  else
  {
    // use the default file/filename 
    sprintf(in_dat,"%s/%s",opts->ipath,opts->ifname);    
    sprintf(filename,"%s%d.in.st",in_dat,myrank);
    fp = fopen(filename,"r");
    if(fp==NULL)
    {
      sprintf(filename,"%s%d.in.st",in_dat,0);
      fp = fopen(filename,"r");      
    }      
  }
  
  long npres = 0;
  
  fscanf (fp,"%lf %ld %ld %ld",&(SOL[0].nor_min),&(SOL[0].iter_max),&npres,&(SOL[0].FNR));
  if (SOL[0].FNR == 2 || SOL[0].FNR == 3)
    fscanf (fp,"%lf %lf",&(arc->dAL0),&(arc->dALMAX));

  // CRYSTAL PLASTICITY
  if(opts->analysis_type == FS_CRPL) {
    crpl = (CRPL*) PGFEM_calloc (mat->nmat,sizeof(CRPL));
    read_cryst_plast(fp,mat->nmat,crpl,opts->plc);
  }
  
  err += read_time_steps(fp, ts);
  
  // loading history exists in load directory
  char load_path[1024];
  char load_fn[1024];
  sprintf(load_path,"%s/load",opts->ipath);

  int is_load_exist = 0;
  for(int ia=0; ia<mp->physicsno; ia++)
  {
    sprintf(load_fn,"%s/%s.load",load_path,mp->physicsname[ia]);
    load->solver_file[ia] = NULL;
    load->solver_file[ia] = fopen(load_fn, "r"); // Load increments are needed to be read
                                             // while time is elapsing.
                                             // This file point needs to be freed end of the simulation
                                             // by calling destruction of the LOADING_STEPS 
    
    if(load->solver_file[ia]==NULL)
      continue;
    
    is_load_exist = 1;
    long nlod_tim = 0;
    fscanf(load->solver_file[ia],"%ld",&nlod_tim);
    
    // read times dependent load
    load->tim_load[ia] = NULL;
    load->tim_load[ia] = compute_times_load(load->solver_file[ia],ts->nt,nlod_tim);
  }

  // if no load directory exists, 
  // loading history will be read and saved only in the 1st physics
  // from the lagacy solver file.
  if(is_load_exist==0)
  { 
    long nlod_tim = 0;
    fscanf (fp,"%ld",&nlod_tim);
    /* read times dependent load */
    load->tim_load[0] = compute_times_load(fp,ts->nt,nlod_tim);
    load->solver_file[0] = fp; // load increments are still need to be read 
                               // while time is elapsing
                               // this file point needs to be freed end of the simulation
                               // by calling destruction of the LOADING_STEPS
  }

  for(int ia=0; ia<mp->physicsno; ia++)
  {
    // if loading is not defined, set default (all zeros)
    if(load->solver_file[ia] == NULL)
    {  
      if(ts->nt == 0)
        load->tim_load[ia] = aloc1l(1);
      else
        load->tim_load[ia] = aloc1l(ts->nt);
    }      
  }
  
  // update numerical solution scheme parameters and field variables
  for(int iA=1; iA<mp->physicsno; iA++)
  {
    // apply only mechanical part
    if(mp->physics_ids[iA] == MULTIPHYSICS_MECHANICAL)
      FV[iA].npres = npres;
          
    if(iA==0)
      continue;
      
    SOL[iA].nor_min  = SOL[0].nor_min;
    SOL[iA].iter_max = SOL[0].iter_max; 
    SOL[iA].FNR      = SOL[0].FNR;
  }
                            
  return 0;
}

/// Read initial conditions from lagcy format.
///
/// If no restart, this function reads initial conditions from *.initial files. 
/// The file format can be found at the following link:
/// https://gitlab-cswarm.crc.nd.edu/pgfem_3d/pgfem_3d/wikis/how-to-set-initial-values
/// If restart is set from the commend line, *.inital files are used only for reading 
/// mid point rule and material densities. Initial conditions are set by reading restart files.   
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] fv object for field variables
/// \param[out] sol object for solution scheme
/// \param[out] load object for loading
/// \param[out] ts object for time stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[out] tnm1 if restart, read time step info from the previous run
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_initial_values_lagcy(GRID *grid,
                              MATERIAL_PROPERTY *mat,
                              FIELD_VARIABLES *fv,
                              SOLVER_OPTIONS *sol,
                              LOADING_STEPS *load,
                              PGFem3D_TIME_STEPPING *ts,
                              PGFem3D_opt *opts,
                              MULTIPHYSICS *mp,
                              double *tnm1, 
                              int myrank)
{
  int err = 0;
  int mp_id = 0;
  
  char filename[1024];
  char line[1024];
  double dt = ts->times[1] - ts->times[0];;
    
  sprintf(filename,"%s/%s%d.initial",opts->ipath,opts->ifname,0);
  
  // restart option from command line is -1
  // check restart form initial file.
  if(opts->restart < 0)
  {
    FILE *fp_0 = fopen(filename,"r");
    
    if(fp_0 != NULL)
    {
      while(fgets(line, 1024, fp_0)!=NULL)
      {
        if(line[0]=='#')
          continue;
        
        sscanf(line, "%d", &(opts->restart));
        break;
      }
      fclose(fp_0);
    }
  }
  
  // check restart and read values
  if(opts->restart >= 0)
    err += read_restart(grid,fv,load,opts,mp,tnm1,myrank);
  
  sprintf(filename,"%s/%s%d.initial",opts->ipath,opts->ifname,myrank);
  FILE *fp = fopen(filename,"r");
  
  if(fp == NULL)
  {
    if(myrank==0)
      printf("Fail to open file [%s]. Quasi steady state\n", filename);
    return 0;
  }
  else
  {
    if(opts->analysis_type == CM && opts->cm == UPDATED_LAGRANGIAN)
    {
      opts->cm = TOTAL_LAGRANGIAN;
      if(myrank==0)
      {
        printf("Updated Lagrangian is currently unavailable with inertia.\n");
        printf("Forced to Total Lagrangian (-cm = %d)\n", TOTAL_LAGRANGIAN);
      }
    }
  }
  
  if(myrank==0)
  {
    while(fgets(line, 1024, fp)!=NULL)
    {
      if(line[0]=='#')
        continue;
      
      double temp;
      sscanf(line, "%lf", &temp);
      break;
    }
  }
  
  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
      continue;
    
    sscanf(line, "%lf", &(sol[mp_id].alpha));
    break;
  }
  
  // read material density
  double *rho = (double *) malloc(sizeof(double)*mat->nmat);  
  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
      continue;
    for(int a=0; a<mat->nmat; a++)
    {
      sscanf(line, "%lf", rho+a);
      if(a<mat->nmat-1)
        fgets(line, 1024, fp);
    }
    break;
  }

  for(int ia = 0; ia<mat->nhommat; ia++)
  {
    (mat->hommat[ia]).density = rho[(mat->hommat[ia]).mat_id];
    if(myrank==0)
      printf("Density(%d), %e\n", ia, rho[(mat->hommat[ia]).mat_id]);
  }  
  
  free(rho);
  
  if(opts->restart<0)
  {  
    while(fgets(line, 1024, fp)!=NULL)
    {
      if(line[0]=='#')
        continue;
    
      long nid;
      double u[3], v[3];
      sscanf(line, "%ld %lf %lf %lf %lf %lf %lf", &nid, u+0, u+1, u+2, v+0, v+1, v+2);
    
      fv[mp_id].u_n[nid*3+0] = u[0];
      fv[mp_id].u_n[nid*3+1] = u[1];
      fv[mp_id].u_n[nid*3+2] = u[2];
      fv[mp_id].u_nm1[nid*3+0] = u[0]-dt*v[0];
      fv[mp_id].u_nm1[nid*3+1] = u[1]-dt*v[1];
      fv[mp_id].u_nm1[nid*3+2] = u[2]-dt*v[2];
    }
  }
  
  for(long idx_a = 0; idx_a<grid->nn; idx_a++)
  {
    for(long idx_b = 0; idx_b<fv[mp_id].ndofn; idx_b++)
    {
      long id = grid->node[idx_a].id_map[mp_id].id[idx_b];
      if(id>0)
        fv[mp_id].u_np1[id-1] = fv[mp_id].u_n[idx_a*fv[mp_id].ndofn + idx_b];
    }
  }  
  
  fclose(fp);
  return err;
}


int read_initial_for_Mechanical(FILE *fp,
                                GRID *grid,
                                MATERIAL_PROPERTY *mat, 
                                FIELD_VARIABLES *fv,
                                SOLVER_OPTIONS *sol,
                                PGFem3D_TIME_STEPPING *ts,
                                PGFem3D_opt *opts,
                                int myrank,
                                int mp_id)
{
  int err = 0;
  char line[1024];
  double dt = ts->times[1] - ts->times[0];;
    
  if(opts->analysis_type == CM && opts->cm == UPDATED_LAGRANGIAN)
  {
    opts->cm = TOTAL_LAGRANGIAN;
    if(myrank==0)
    {
      printf("Updated Lagrangian is currently unavailable with inertia.\n");
      printf("Forced to Total Lagrangian (-cm = %d)\n", TOTAL_LAGRANGIAN);
    }
  }
  
  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
      continue;
    
    sscanf(line, "%lf", &(sol->alpha));
    break;
  }
  
  // read material density density
  double *rho = (double *) malloc(sizeof(double)*mat->nmat);  
  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
      continue;
    for(int a=0; a<mat->nmat; a++)
    {
      sscanf(line, "%lf", rho+a);
      if(a<mat->nmat-1)
        fgets(line, 1024, fp);
    }
    break;
  }

  for(int ia = 0; ia<mat->nhommat; ia++)
  {
    (mat->hommat[ia]).density = rho[(mat->hommat[ia]).mat_id];
    if(myrank==0)
      printf("Density(%d), %e\n", ia, rho[(mat->hommat[ia]).mat_id]);
  }  
  
  free(rho);

  if(opts->restart < 0)
  {      
    while(fgets(line, 1024, fp)!=NULL)
    {
      if(line[0]=='#')
        continue;
    
      long nid;
      double u[3], v[3];
      sscanf(line, "%ld %lf %lf %lf %lf %lf %lf", &nid, u+0, u+1, u+2, v+0, v+1, v+2);
    
      fv->u_n[nid*3+0] = u[0];
      fv->u_n[nid*3+1] = u[1];
      fv->u_n[nid*3+2] = u[2];
      fv->u_nm1[nid*3+0] = u[0]-dt*v[0];
      fv->u_nm1[nid*3+1] = u[1]-dt*v[1];
      fv->u_nm1[nid*3+2] = u[2]-dt*v[2];
    }
  }
  
  for(long idx_a = 0; idx_a<grid->nn; idx_a++)
  {
    for(long idx_b = 0; idx_b<fv->ndofn; idx_b++)
    {
      long id = grid->node[idx_a].id_map[mp_id].id[idx_b];
      if(id>0)
        fv->u_np1[id-1] = fv->u_n[idx_a*fv->ndofn + idx_b];
    }
  }  
  
  return err;
}


int read_initial_for_Thermal(FILE *fp,
                             GRID *grid,
                             MATERIAL_PROPERTY *mat, 
                             FIELD_VARIABLES *fv,
                             SOLVER_OPTIONS *sol,
                             PGFem3D_TIME_STEPPING *ts,
                             PGFem3D_opt *opts,
                             int myrank,
                             int mp_id)
{
  int err = 0;
  char line[1024];

  double T0 = 300.0;
  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
      continue;
    
    sscanf(line, "%lf", &T0);
    fv->u0 = T0; // set reference temperature
    if(myrank==0)
      printf("Default initial temperature: %e\n", T0);

    break;
  }  

  if(opts->restart < 0)
  {      
    // set default
    for(int ia=0; ia<grid->nn; ia++)
    {
      fv->u_nm1[ia] = T0;
      fv->u_n[ia] = T0;
    }
      
    while(fgets(line, 1024, fp)!=NULL)
    {
      if(line[0]=='#')
        continue;
    
      long nid;
      double u;
      sscanf(line, "%ld %lf", &nid, &u);
    
      fv->u_n[nid] = u;
      fv->u_nm1[nid] = u;
    }
  }
  
  for(int ia = 0; ia<grid->nn; ia++)
  {
    long id = grid->node[ia].id_map[mp_id].id[0];
    if(id>0)
      fv->u_np1[id-1] = fv->u_n[ia];    
  }  
  
  return err;
}

/// Read initial conditions from lagcy format.
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[out] SOL array of solution scheme object
/// \param[out] load object for loading
/// \param[out] ts object for time stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[out] tnm1 if restart, read time step info from the previous run
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_initial_values_IC(GRID *grid,
                           MATERIAL_PROPERTY *mat,
                           FIELD_VARIABLES *FV,
                           SOLVER_OPTIONS *SOL,
                           LOADING_STEPS *load,
                           PGFem3D_TIME_STEPPING *ts,
                           PGFem3D_opt *opts,
                           MULTIPHYSICS *mp,
                           double *tnm1, 
                           int myrank)
{
  int err = 0;
  int mp_id = 0;
  
  // check restart and read restart values
 if(opts->restart >= 0)
    err += read_restart(grid,FV,load,opts,mp,tnm1,myrank);  
  
  char IC[1024];
  sprintf(IC,"%s/IC",opts->ipath);
  
  char fn_0[1024], fn[1024];  
  
  for(int ia=0; ia<mp->physicsno; ia++)
  {
    sprintf(fn_0,"%s/%s_0.initial",IC,mp->physicsname[ia]);
    sprintf(fn  ,"%s/%s_%d.initial",IC,mp->physicsname[ia], myrank);
    
    FILE *fp = NULL;
    fp = fopen(fn, "r");
    
    if(fp==NULL)
    {  
      fp = fopen(fn_0, "r");
      if(fp==NULL)
      { 
        if(myrank==0) 
          printf("No [%s] exists. Use default ICs.\n", fn_0);

        continue;
      }
    }
    switch(mp->physics_ids[ia])
    {
      case MULTIPHYSICS_MECHANICAL:
        err += read_initial_for_Mechanical(fp,grid,mat,FV+ia,SOL+ia,ts,opts,myrank,ia);
        break;
      case MULTIPHYSICS_THERMAL:
        err += read_initial_for_Thermal(fp,grid,mat,FV+ia,SOL+ia,ts,opts,myrank,ia);        
        break;
      case MULTIPHYSICS_CHEMICAL:
        break;
      default:
        break;  
    }
    
    fclose(fp);
  }
  return err;
}



/// Read initial conditions.
///
/// \param[in] grid a mesh object
/// \param[in] mat a material object
/// \param[in,out] FV array of field variable object
/// \param[out] SOL array of solution scheme object
/// \param[out] load object for loading
/// \param[out] ts object for time stepping
/// \param[in] opts structure PGFem3D option
/// \param[in] mp mutiphysics object
/// \param[out] tnm1 if restart, read time step info from the previous run
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_initial_values(GRID *grid,
                        MATERIAL_PROPERTY *mat,
                        FIELD_VARIABLES *FV,
                        SOLVER_OPTIONS *SOL,
                        LOADING_STEPS *load,
                        PGFem3D_TIME_STEPPING *ts,
                        PGFem3D_opt *opts,
                        MULTIPHYSICS *mp,
                        double *tnm1, 
                        int myrank)
{
  int err = 0;
  char IC[1024];
  sprintf(IC,"%s/IC",opts->ipath);
  
  if(is_directory_exist(IC))
  { 
    if(myrank==0) 
      printf("IC directory exists, read initial conditions from IC\n");
    err += read_initial_values_IC(grid,mat,FV,SOL,load,ts,opts,mp,tnm1,myrank);
  }     
  else
  {
    if(myrank==0)  
      printf("No IC directory exists, read inital conditions from *.initial\n");

    err += read_initial_values_lagcy(grid,mat,FV+0,SOL+0,load,ts,opts,mp,tnm1,myrank);
  }
  return err;
}

/// Read loads increments.
/// As time is elapsing, loads increments are read from solver file which 
/// file pointer is saved in LOADING_STEPS. Prior to run this function, 
/// read_initial_values function shold be called which open and the solver file pointer.
/// The file pointer will be freed when LOADING_STEPS object is destoryed.
/// The number of loads increments should be exact as read before in read_initial_values.
///
/// \param[in] grid a mesh object
/// \param[in] fv object for field variables
/// \param[out] load object for loading
/// \param[in] mp multiphysics object
/// \param[in] tim time step ID
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \return non-zero on internal error 
int read_and_apply_load_increments(GRID *grid,
                                   FIELD_VARIABLES *fv,
                                   LOADING_STEPS *load,
                                   MULTIPHYSICS *mp, 
                                   long tim, 
                                   MPI_Comm mpi_comm,
                                   int myrank)
{
  int err = 0;

  //  read nodal prescribed boundary values
  for(int mp_id=0; mp_id<mp->physicsno; mp_id++)
  {
    if(load->solver_file[mp_id]==NULL)
      continue;
      
    if(load->tim_load[mp_id][tim] == 1 && tim == 0)
    {
      if (myrank == 0)
        PGFEM_printf ("Incorrect load input for Time = 0\n");

      PGFEM_Comm_code_abort(mpi_comm,0);
    }
    
    if(load->tim_load[mp_id][tim] == 1 && tim != 0)
    {  
      for(long ia=0;ia<load->sups[mp_id]->npd;ia++)
      {
        fscanf(load->solver_file[mp_id],"%lf",(load->sups[mp_id])->defl_d + ia);
        (load->sup_defl[mp_id])[ia] = load->sups[mp_id]->defl_d[ia];
      }
      if(mp->physics_ids[mp_id]==MULTIPHYSICS_MECHANICAL)
      {  
        // read nodal load in the subdomain
        read_nodal_load(load->solver_file[mp_id],load->nln,grid->nsd,load->znod);
        // read elem surface load */
        read_elem_surface_load(load->solver_file[mp_id],load->nle_s,grid->nsd,grid->element,load->zele_s);
        //  NODE - generation of the load vector 
        load_vec_node(fv[mp_id].R,load->nln,grid->nsd,load->znod,grid->node,MULTIPHYSICS_MECHANICAL);
        //  ELEMENT - generation of the load vector
        load_vec_elem_sur(fv[mp_id].R,load->nle_s,grid->nsd,grid->element,load->zele_s);
      }
    }
  }
          
  return err;
}

/// Read read cohesive elements.
/// This function was part of the main function and extracted to 
/// Modularized reading cohesive elements and reorganize the main function structures.
/// The function call read not only cohesive elements but also cohesive properties (*.in.co_props).
/// The functin dependency includes lagacy code (read_cohe_elem) to read *.in.co files, too.
///
/// \param[out] grid a mesh object
/// \param[out] mat a material object
/// \param[in] opts structure PGFem3D option
/// \param[in] ensight ENSIGHT object
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] myrank current process rank
/// \return non-zero on internal error  
int read_cohesive_elements(GRID *grid,
                           MATERIAL_PROPERTY *mat,
                           const PGFem3D_opt *opts,
                           ENSIGHT ensight,
                           MPI_Comm mpi_comm,
                           int myrank)
{
  int err = 0;
  char in_dat[1024], filename[1024];
  
  sprintf(in_dat,"%s/%s",opts->ipath,opts->ifname);

  FILE *fp;
  // read cohesive properties
  sprintf(filename,"%s%d.in.co_props",in_dat,myrank);
  fp = fopen(filename,"r");
  read_cohesive_properties(fp,&(mat->n_co_props),&(mat->co_props),mpi_comm);
  fclose(fp);
  
  /* read coheisve elements */
  sprintf(filename,"%s%d.in.co",in_dat,myrank);
  fp = fopen(filename,"r");
  
  long ncom = 0;
  double **comat = NULL;
  
  /* temporary leftovers from old file format */
  fscanf(fp,"%ld\n",&ncom);
  
  /* to silence warning message. need to pull this legacy bit of
   * code out completely. Cohesive porperties provided in separate
   * file. This leads to *very* small memory leak */
  if(ncom <= 0) 
    comat = aloc2(1,4);
  else
    comat = aloc2(ncom,4);
  
  
  /* read the cohesive element info */
  grid->coel = read_cohe_elem (fp,ncom,grid->nsd,grid->nn,grid->node,&(grid->nce),
          comat,ensight,opts->vis_format,
          myrank,mat->co_props);
  if(ncom <= 0)
    dealoc2(comat,ncom);
  else
    dealoc2(comat,1);  

  fclose (fp);
  
  /* Global number of cohesive elements */
  MPI_Allreduce(&(grid->nce),&(grid->Gnce),1,MPI_LONG,MPI_SUM,mpi_comm);
  return err;
}   