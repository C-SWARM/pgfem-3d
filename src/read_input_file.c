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

static const int ndim = 3;

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
		    ZATELEM **zelem_v)
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
  default: *ndofn = ndim; break;
  }

  (*node) = build_node(*nn,*ndofn);
  (*elem) = build_elem(in,*ne,opts->analysis_type);
  (*material) = PGFEM_calloc(*nmat,sizeof(MATERIAL));
  (*matgeom) = build_matgeom(*n_concentrations,*n_orient);

  *Gnn = read_nodes(in,*nn,*node,opts->legacy,comm);
  /* NOTE: Supports assume only ndim supported dofs per node! */
  *sup = read_supports(in,*nn,ndim,*node);

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
  *znod = build_zatnode (ndim,*nln);
  read_nodal_load (in,*nln,ndim,*znod);
  /* surface */
  fscanf (in,"%ld",nel_s);
  *zelem_s = build_zatelem (ndim,*nel_s);
  read_elem_surface_load (in,*nel_s,ndim,*elem,*zelem_s);
  /* volume */
  fscanf (in,"%ld",nel_v);
  *zelem_v = build_zatelem (ndim,*nel_v);

  /* check the ferror bit */
  if(ferror(in)) err++;

  /* free local memory and close file */
  free(filename);
  fclose(in);
  return err;
}

/// read mesh file
///
/// \param[out] grid a mesh object
/// \param[out] mat a material object
/// \param[out] variables a mesh object
/// \param[out] sol a mesh object
/// \param[out] load a mesh object
/// \param[in] comm MPI_COMM_WORLD
/// \param[in] opts structure PGFem3D option
/// \return non-zero on internal error
int read_mesh_file(GRID *grid, 
                   MATERIAL_PROPERTY *mat,
                   FIELD_VARIABLES *variables,
                   SOLVER_OPTIONS *sol,
                   LOADING_STEPS *load,
                   MPI_Comm mpi_comm,
                   const PGFem3D_opt *opts)
{
  int err = read_input_file(opts,
                            mpi_comm,
                            &(grid->nn),
                            &(grid->Gnn),
                            &(variables->ndofn),
                            &(grid->ne),
                            &(sol->iter_max_sol),
                            &(sol->err),
                            &(sol->computer_zero),
                            &(mat->nmat),
                            &(variables->n_concentrations),
                            &(mat->n_orient),
                            &(grid->node),
                            &(grid->element),
                            &(mat->mater),
                            &(mat->matgeom),
                            &(load->sup),
                            &(load->nln),
                            &(load->znod),
                            &(load->nle_s),
                            &(load->zele_s),
                            &(load->nle_v),
                            &(load->zele_v));
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

/// read solver file for time stepping
///
///////////////////////////////////////////////////////////////
// read input file name
int read_solver_file(PGFem3D_TIME_STEPPING *ts,
                     MATERIAL_PROPERTY *mat,
                     FIELD_VARIABLES *variables,
                     SOLVER_OPTIONS *sol,
                     LOADING_STEPS *load,
                     ARC_LENGTH_VARIABLES *arc,
                     CRPL *crpl,
                     const PGFem3D_opt *opts,
                     int myrank)
{
  int err = 0;
  // READ SOLVER FILE
  // override the default solver file with one specified
  // at commandline
  char filename[1024];
  char in_dat[1024];
  FILE *fp;
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
  }
  
  fscanf (fp,"%lf %ld %ld %ld",&(sol->nor_min),&(sol->iter_max),&(variables->npres),&(sol->FNR));
  if (sol->FNR == 2 || sol->FNR == 3)
    fscanf (fp,"%lf %lf",&(arc->dAL0),&(arc->dALMAX));

  // CRYSTAL PLASTICITY
  if(opts->analysis_type == FS_CRPL) {
    crpl = (CRPL*) PGFEM_calloc (mat->nmat,sizeof(CRPL));
    read_cryst_plast(fp,mat->nmat,crpl,opts->plc);
  }
  
  err += read_time_steps(fp, ts);
  
  long nlod_tim = 0;
  fscanf (fp,"%ld",&nlod_tim);
  /* read times dependent load */
  load->tim_load = compute_times_load(fp,ts->nt,nlod_tim);
  load->solver_file = fp; // load increments are still need to be read 
                          // while time is elapsing
                          // this file point needs to be freed end of the simulation
                          // by calling destruction of the LOADING_STEPS
  return 0;
}

int read_initial_values(GRID *grid,
                        MATERIAL_PROPERTY *mat,
                        FIELD_VARIABLES *fv,
                        SOLVER_OPTIONS *sol,
                        LOADING_STEPS *load,
                        PGFem3D_TIME_STEPPING *ts,
                        PGFem3D_opt *opts,
                        int *restart, 
                        double *tnm1, 
                        int myrank)
{
  int err = 0;
  
  char filename[1024];
  char line[1024];
  double dt = ts->times[1] - ts->times[0];;
    
  sprintf(filename,"%s/%s%d.initial",opts->ipath,opts->ifname,0);
  
  // restart option from command line is -1
  // check restart form initial file.
  if(*restart < 0)
  {
    FILE *fp_0 = fopen(filename,"r");
    
    if(fp_0 != NULL)
    {
      while(fgets(line, 1024, fp_0)!=NULL)
      {
        if(line[0]=='#')
          continue;
        
        sscanf(line, "%d", restart);
        break;
      }
      fclose(fp_0);
    }
  }
  
  // check restart and read values
  if(*restart >= 0)
  {
    int nsd = 3;
    read_restart(fv->u_nm1,fv->u_n,opts,grid->element,grid->node,fv->sig,fv->eps,load->sup,
            myrank,grid->ne,grid->nn,nsd,restart,tnm1, &(fv->NORM));
  }
  
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
  
  if(*restart>0)
  {
    fclose(fp);
    return 0;
  }
  
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
  
  fclose(fp);
  return err;
}

int read_and_apply_load_increments(GRID *grid,
                                   FIELD_VARIABLES *variables,
                                   LOADING_STEPS *load, 
                                   long tim, 
                                   MPI_Comm mpi_comm,
                                   int myrank)
{
  int err = 0;
  
  if (load->tim_load[tim] == 1 && tim == 0)
  {
    if (myrank == 0)
      PGFEM_printf ("Incorrect load input for Time = 0\n");

    PGFEM_Comm_code_abort(mpi_comm,0);
  }
  if (load->tim_load[tim] == 1 && tim != 0)
  {
    //  read nodal prescribed deflection
    for(long ia=0;ia<load->sup->npd;ia++)
    {
      fscanf (load->solver_file,"%lf",((load->sup)->defl_d) + ia);
      (load->sup_defl)[ia] = (load->sup)->defl_d[ia];
    }
    // read nodal load in the subdomain
    read_nodal_load(load->solver_file,load->nln,grid->nsd,load->znod);
    // read elem surface load */
    read_elem_surface_load(load->solver_file,load->nle_s,grid->nsd,grid->element,load->zele_s);
    //  NODE - generation of the load vector 
    load_vec_node(variables->R,load->nln,ndim,load->znod,grid->node);
    //  ELEMENT - generation of the load vector
    load_vec_elem_sur(variables->R,load->nle_s,grid->nsd,grid->element,load->zele_s);
    
  } /* end load increment */
          
  return err;
}


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