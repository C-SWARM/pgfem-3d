#include "restart.h"
#include "PGFem3D_options.h"
#include "gen_path.h"
#include "element.h"
#include "node.h"
#include "supp.h"
#include "vtk_output.h"
#include "constitutive_model.h"
#include "data_structure_c.h"
#include "elem3d.h"
#include "PGFem3D_data_structure.h"


#ifndef NO_VTK_LIB
#include "PGFem3D_to_VTK.hpp"

/// read inital values from VTK files, if VTK library is used,
/// it return -1 restart number, such that run will start from t(n=0)
///
/// \param[in] opts structure PGFem3D option 
/// \param[in] myrank current process rank
/// \param[out] u0 displacement at t(n-1)
/// \param[out] u1 displacement at t(n)
/// \param[in] rs_path directory path for restart files
/// \return non-zero on internal error 
int read_initial_from_VTK(PGFem3D_opt *opts, int myrank, double *u0, double *u1, char *rs_path)
{
  int err = 0;
  char filename[1024];
  int stepno = opts->restart;
  
  sprintf(filename,"%s/VTK/STEP_%.6d/%s_%d_%d.vtu",rs_path,stepno,opts->ofname,myrank,stepno);   
  err += read_VTK_file(filename, u0);      
  sprintf(filename,"%s/VTK/STEP_%.6d/%s_%d_%d.vtu",rs_path,stepno,opts->ofname,myrank,stepno);   
  err += read_VTK_file(filename, u1);
    
  return err;
}      

#else
// in case VTK library is not used.
int read_initial_from_VTK(PGFem3D_opt *opts, int myrank, double *u0, double *u1, char *rs_path)
{
  if(myrank==0)
  {
    PGFEM_printerr("Restart with VTK is not supported!\n");
    PGFEM_printerr("Enforce to turn off restart!\n");
  }
  
  opts->restart = -1;  
  return 0;
}
#endif

/// read time stepping information
///
/// When subdivision is made, t(n-1) and t(n) are different than read values from
/// solver file, so that restart requires subdivision information in order to start the run
/// using saved restart values at t(n). Otherwise, artificial acceleration will be generated due to
/// large time step size.
///
/// \param[in,out] fv array of field variable object, fv.NORM will be updated
/// \param[in,out] time_steps object for time stepping, time_steps.tns are updated
/// \param[in] opts PGFem3D commend line options
/// \param[in] mp multiphysics object
/// \param[out] tnm1 times at t(n-1), t(n) 
/// \param[in] myrank current process rank
/// \return non-zero on internal error 
int read_time_step_info(FIELD_VARIABLES *fv,
                        PGFem3D_TIME_STEPPING *time_steps,
                        const PGFem3D_opt *opts,
                        MULTIPHYSICS *mp,
                        double *tnm1,
                        int myrank)
{ 
  int err = 0;
   
  char fn[1024];
  sprintf(fn, "%s/restart/time_step_info_%.6d.res",opts->opath,opts->restart);

  double t[3];
  t[0] = t[1] = t[2] = -1.0;
  
  FILE *fp = fopen(fn, "r");    
  
  if(fp != NULL)
  { 
    fscanf(fp, "%lf %lf %lf", t+0, t+1, t+2);
    tnm1[0] = t[0];
    tnm1[1] = t[1];

    if(myrank==0)
      printf("read time stpe info t(n-1)=%e, t(n)=%e, t(n+1) = %e\n", t[0], t[1], t[2]); 
    
    for(int ia=0; ia<mp->physicsno; ia++)
    {
      fscanf(fp, "%lf %lf", &(fv[ia].NORM), time_steps->tns+ia);    
      if(myrank==0)
        printf("\t\t%s: NORM = %e, t(n) = %e\n",mp->physicsname[ia], fv[ia].NORM, time_steps->tns[ia]); 
    }
    if(myrank==0)
      printf("\n");
      
    fclose(fp);
  }
  else
  {
    if(myrank==0)
      printf("WARNING: cannot read time steps info [%s] \n", fn);
  }
  return err;
}

/// write time stepping information
///
/// A file will be generated about subdivision; t(n-1) and t(n) are different than read values from
/// solver file, so that restart requires subdivision information in order to start the run
/// using saved restart values at t(n). Otherwise, artificial acceleration will be generated due to
/// large time step size.
//
/// \param[in] fv array of field variable object
/// \param[in] time_steps object for time stepping
/// \param[in] opts PGFem3D commend line options
/// \param[in] mp multiphysics object
/// \param[in] myrank current process rank
/// \param[in] stepno current time step number
/// \return non-zero on internal error 
int write_time_step_info(FIELD_VARIABLES *fv,
                         PGFem3D_TIME_STEPPING *time_steps,
                         const PGFem3D_opt *opts,
                         MULTIPHYSICS *mp, 
                         int myrank,
                         int stepno)
{
  int err = 0;
  if(myrank==0)
  {
    double *times = time_steps->times;

    // write time stepping info
    char fn[1024];
    sprintf(fn, "%s/restart/time_step_info_%.6d.res",opts->opath,stepno);
      
    FILE *fp = fopen(fn, "w");
    if(fp==NULL)
    {
      printf("Cannot create a file [%s]\n", fn);
      printf("Anyway continue ...\n");
    }
    else
    {  
      if(stepno>0)
        fprintf(fp, "%.17e %.17e %.17e ", times[stepno-1], times[stepno], times[stepno+1]);  
      else
        fprintf(fp, "0.0 %.17e %.17e ", times[stepno], times[stepno+1]);
      
      for(int ia=0; ia<mp->physicsno; ia++)
        fprintf(fp, "\n%.17e %.17e",fv[ia].NORM, time_steps->tns[ia]);

      fprintf(fp, "\n");
                          
      fclose(fp);
    }
  }
  return err;  
}

/// write restart file for constitutive model interface
///
/// First, desplacements are written, but afterwards depending on the constitutive model, 
/// restart file is written differently.
/// Each file format is defined in each constitutive model. 
///
/// \param[in] grid a mesh object
/// \param[in] fv array of field variable object
/// \param[in] time_steps object for time stepping
/// \param[in] opts PGFem3D commend line options
/// \param[in] mp multiphysics object
/// \param[in] myrank current process rank
/// \param[in] mp_id multiphysics id
/// \param[in] stepno current time step number
/// \param[in] rs_path directory path for restart files
/// \return non-zero on internal error 
int write_restart_constitutive_model(GRID *grid,
                                     FIELD_VARIABLES *fv,
                                     const PGFem3D_opt *opts,
                                     MULTIPHYSICS *mp, 
                                     int myrank,
                                     int mp_id,
                                     int stepno,
                                     char rs_path[1024])                                     
{
  int err = 0;

  char restart_path[1024];
  sprintf(restart_path, "%s/STEP_%.6d", rs_path,stepno);
                
  if(make_path(restart_path,DIR_MODE) != 0)
  {
    PGFEM_printf("Directory (%s) not created!\n",restart_path);
    abort();                   
  }  
 
  char filename[1024];
  sprintf(filename,"%s/STEP_%.6d/%s_%d_%d.res",rs_path,stepno,opts->ofname,myrank, stepno);
  FILE *fp = fopen(filename,"w");

  if(fp == NULL)
  {    
    printf("Fail to open file [%s]. finishing\n", filename);
    exit(0);  
  }
  
  for(int a=0; a<grid->nn; a++)
  {
    for(int b=0; b<grid->nsd; b++) {
      fprintf(fp, "%.17e %.17e ", fv[mp_id].u_nm1[a*(grid->nsd) + b], fv[mp_id].u_n[a*(grid->nsd) + b]);
    }
    fprintf(fp, "\n");    
  }
  
  if(opts->analysis_type==CM)
  {  
    for (int e = 0; e < grid->ne; e++)
    {
      const ELEMENT *p_el = grid->element + e;
      long n_ip = 0;
      int_point(p_el->toe,&n_ip);
      fprintf(fp, "%ld\n", n_ip);
      for (int ip = 0; ip < n_ip; ip++)
      {
        Constitutive_model *m = &(fv[mp_id].eps[e].model[ip]);
        err += m->param->write_restart(fp, m);
      }
    }
  }

  fclose(fp);
  return err;
}

/// read restart file for constitutive model interface
///
/// Depending on constitutive mode, reading restart file has different format except displacement.
/// Each file format is defined in each constitutive model. 
///
/// \param[in] grid a mesh object
/// \param[in, out] fv array of field variable object
/// \param[in] opts PGFem3D commend line options
/// \param[in] mp multiphysics object
/// \param[in] myrank current process rank
/// \param[in] mp_id multiphysics id
/// \param[in] rs_path directory path for restart files
/// \return non-zero on internal error 
int read_restart_constitutive_model(GRID *grid,
                                    FIELD_VARIABLES *fv,
                                    const PGFem3D_opt *opts,
                                    MULTIPHYSICS *mp, 
                                    int myrank,
                                    int mp_id,
                                    char rs_path[1024]) 
{
  int err = 0;
  int stepno = opts->restart;
    
  char filename[1024];
  sprintf(filename,"%s/STEP_%.6d/%s_%d_%d.res",rs_path,stepno,opts->ofname,myrank, stepno);   
  FILE *fp = fopen(filename,"r");

  if(fp == NULL)
  {    
    printf("Fail to open file [%s]. finishing\n", filename);      
    exit(0);  
  }
  
  for(int a=0; a<grid->nn; a++)
  {
    for(int b=0; b<grid->nsd; b++)
      fscanf(fp, "%lf %lf", (fv[mp_id].u_nm1)+a*(grid->nsd) + b, (fv[mp_id].u_n)+a*(grid->nsd) + b);    
  }
  
  if(opts->analysis_type==CM)
  {      
    for (int e = 0; e < grid->ne; e++) 
    {
      const ELEMENT *p_el = grid->element + e;
      long n_ip = 0;
      long n_ip_read = 0;
      int_point(p_el->toe,&n_ip);
      fscanf(fp, "%ld", &n_ip_read); 
      if(n_ip!=n_ip_read)
      {  
        printf("Error: restart file has wrong integration number (PN: %d, elem: %d)\n", myrank, e);
        return -1;
      }
      for (int ip = 0; ip < n_ip; ip++)     
      {
        Constitutive_model *m = &(fv[mp_id].eps[e].model[ip]);
        err += m->param->read_restart(fp, m);
      }
    }
  }  
        
  fclose(fp);  
  return err;  
}

/// read restart files for mechanical part
///
/// Read restart values at t(n-1) and t(n) in the order as: displacement, 
/// element values according to the analysis type.
///
/// \param[in] grid a mesh object
/// \param[in, out] fv array of field variable object
/// \param[in] load object for loading
/// \param[in] opts PGFem3D commend line options
/// \param[in] mp multiphysics object
/// \param[out] tnm1 times at t(n-1), t(n) 
/// \param[in] myrank current process rank
/// \param[in] mp_id multiphysics id
/// \param[in] rs_path directory path for restart files
/// \return non-zero on internal error
int read_restart_mechanical(GRID *grid,
                            FIELD_VARIABLES *fv,
                            LOADING_STEPS *load,
                            PGFem3D_opt *opts,
                            MULTIPHYSICS *mp,                            
                            double *tnm1,
                            int myrank,
                            int mp_id,
                            char rs_path[1024])
{
  int err = 0;
    
  switch(opts->analysis_type)
  {
    case DISP: // intended to flow
    case CM:
      err += read_restart_constitutive_model(grid,fv,opts,mp,myrank,mp_id,rs_path); 
      break;
    default:
      read_initial_from_VTK(opts, myrank, fv[mp_id].u_nm1, fv[mp_id].u_n, rs_path);
      break;
  }
      
  return err;
}

/// read restart files for thermal part
///
/// Read temperatures at t(n-1) and t(n).
///
/// \param[in] grid a mesh object
/// \param[in, out] fv array of field variable object
/// \param[in] opts PGFem3D commend line options
/// \param[in] myrank current process rank
/// \param[in] mp_id multiphysics id
/// \param[in] rs_path directory path for restart files
/// \return non-zero on internal error
int read_restart_thermal(GRID *grid,
                         FIELD_VARIABLES *fv,
                         const PGFem3D_opt *opts,
                         int myrank,
                         int mp_id,
                         char *rs_path)
{
  int err = 0;
  int stepno = opts->restart;
  
  char filename[1024];
  sprintf(filename,"%s/STEP_%.6d/%s_%d_%d.res",rs_path,stepno,opts->ofname,myrank, stepno);   
  FILE *fp = fopen(filename,"r");

  if(fp == NULL)
  {    
    printf("Fail to open file [%s]. finishing\n", filename);      
    exit(0);  
  }
  
  for(int ia=0; ia<grid->nn; ia++)
    fscanf(fp, "%lf %lf\n", (fv[mp_id].u_nm1)+ia, (fv[mp_id].u_n)+ia);

  fclose(fp);
  return err;
}

/// read restart files for multiphysics problem
///
/// Read field variables that are needed to restart at t(n-1) and t(n).
/// By going through all physics, restart files for multiple physics are read 
/// as many as the number of physics.
///
/// \param[in] grid a mesh object
/// \param[in, out] fv array of field variable object
/// \param[in, out] time_steps object for time stepping
/// \param[in] load object for loading
/// \param[in] opts PGFem3D commend line options
/// \param[in] mp multiphysics object
/// \param[out] tnm1 times at t(n-1), t(n) 
/// \param[in] myrank current process rank
/// \return non-zero on internal error
int read_restart(GRID *grid,
                 FIELD_VARIABLES *fv,
                 PGFem3D_TIME_STEPPING *time_steps,
                 LOADING_STEPS *load,
                 const PGFem3D_opt *opts,
                 MULTIPHYSICS *mp,
                 double *tnm1,
                 int myrank)
{
  int err = 0;
  
  char rs_path[1024]; 

  for(int ia=0; ia<mp->physicsno; ia++)
  {
    sprintf(rs_path, "%s/restart/%s", opts->opath,mp->physicsname[ia]);
    if(make_path(rs_path,DIR_MODE) != 0)
    {
      PGFEM_printf("Directory (%s) not created!\n",rs_path);
      abort();                   
    }    
    switch(mp->physics_ids[ia])
    {
      case MULTIPHYSICS_MECHANICAL:
        err += read_restart_mechanical(grid,fv,load,opts,mp,tnm1,myrank,ia,rs_path);
        break;
      case MULTIPHYSICS_THERMAL:
        err += read_restart_thermal(grid,fv,opts,myrank,ia,rs_path);
        break;
      case MULTIPHYSICS_CHEMICAL:
        // not yet implemented
        break;
      default:
        err += read_restart_mechanical(grid,fv,load,opts,mp,tnm1,myrank,ia,rs_path);
    }
  }
  // read time stepping info  
  err += read_time_step_info(fv,time_steps,opts,mp,tnm1,myrank);
    
  return err;
}

/// write restart files for mechanical part
///
/// Write restart values at t(n-1) and t(n) in the order as: displacement, 
/// element values according to the analysis type.
///
/// \param[in] grid a mesh object
/// \param[in, out] fv array of field variable object
/// \param[in] load object for loading
/// \param[in] time_steps object for time stepping
/// \param[in] opts PGFem3D commend line options
/// \param[in] mp multiphysics object
/// \param[in] myrank current process rank
/// \param[in] mp_id multiphysics id
/// \param[in] stepno current time step number
/// \param[in] rs_path directory path for restart files
/// \return non-zero on internal error
int write_restart_mechanical(GRID *grid,
                             FIELD_VARIABLES *fv,
                             LOADING_STEPS *load,
                             PGFem3D_TIME_STEPPING *time_steps,
                             const PGFem3D_opt *opts,
                             MULTIPHYSICS *mp,
                             int myrank,
                             int mp_id,
                             int stepno,
                             char *rs_path)
{
  int err = 0;
  
  switch(opts->analysis_type)
  {
    case DISP: // intended to flow
    case CM:
      err += write_restart_constitutive_model(grid,fv,opts,mp,myrank,mp_id,stepno,rs_path); 
      break;
    default:
    {  
      double *r_n_dof = (double *) malloc(sizeof(double)*(fv[mp_id].ndofd));
      for(long a = 0; a<grid->nn; a++)
      {              
        for(long b = 0; b<fv[mp_id].ndofn; b++)
        {
          long id = grid->node[a].id_map[mp_id].id[b];
          if(id>0)
          r_n_dof[id-1] = fv[mp_id].u_nm1[a*(fv[mp_id].ndofn) + b];
        }
      }            
      VTK_print_vtu(rs_path,opts->ofname,stepno,
                    myrank,grid->ne,grid->nn,grid->node,grid->element,load->sups[mp_id],
                    r_n_dof,fv[mp_id].sig,fv[mp_id].eps,
                    opts, mp_id);
      free(r_n_dof);              
      break;
    }
  }  
  return err;
}

/// write restart files for thermal part
///
/// Write temperatures at t(n-1) and t(n)
///
/// \param[in] grid a mesh object
/// \param[in] fv array of field variable object
/// \param[in] opts PGFem3D commend line options
/// \param[in] myrank current process rank
/// \param[in] mp_id multiphysics id
/// \param[in] stepno current time step number
/// \param[in] rs_path directory path for restart files
/// \return non-zero on internal error
int write_restart_thermal(GRID *grid,
                          FIELD_VARIABLES *fv,
                          const PGFem3D_opt *opts,
                          int myrank,
                          int mp_id,
                          int stepno,
                          char *rs_path)
{
  int err = 0;

  char restart_path[1024];
  sprintf(restart_path, "%s/STEP_%.6d", rs_path,stepno);
                
  if(make_path(restart_path,DIR_MODE) != 0)
  {
    PGFEM_printf("Directory (%s) not created!\n",restart_path);
    abort();                   
  }  
 
  char filename[1024];
  sprintf(filename,"%s/%s_%d_%d.res",restart_path,opts->ofname,myrank, stepno);
  FILE *fp = fopen(filename,"w");

  if(fp == NULL)
  {    
    printf("Fail to open file [%s]. finishing\n", filename);
    exit(0);  
  }
  
  for(int ia=0; ia<grid->nn; ia++)
    fprintf(fp, "%.17e %.17e\n", fv[mp_id].u_nm1[ia], fv[mp_id].u_n[ia]);
  
  fclose(fp);
  return err;
}

/// write restart files for mechanical part
///
/// Write field variables that are needed to restart at t(n-1) and t(n).
/// By going through all physics, restart files for multiple physics are written 
/// as many as the number of physics.
///
/// \param[in] grid a mesh object
/// \param[in, out] fv array of field variable object
/// \param[in] load object for loading
/// \param[in] time_steps object for time stepping
/// \param[in] opts PGFem3D commend line options
/// \param[in] mp multiphysics object
/// \param[in] myrank current process rank
/// \param[in] mp_id multiphysics id
/// \param[in] stepno current time step number
/// \param[in] rs_path directory path for restart files
/// \return non-zero on internal error
int write_restart(GRID *grid,
                  FIELD_VARIABLES *fv,
                  LOADING_STEPS *load,
                  PGFem3D_TIME_STEPPING *time_steps,
                  const PGFem3D_opt *opts,
                  MULTIPHYSICS *mp,
                  int myrank,
                  int stepno)

{
  int err = 0;

  char rs_path[1024];
    
  for(int ia=0; ia<mp->physicsno; ia++)
  {
    sprintf(rs_path, "%s/restart/%s", opts->opath,mp->physicsname[ia]);
    if(make_path(rs_path,DIR_MODE) != 0)
    {
      PGFEM_printf("Directory (%s) not created!\n",rs_path);
      abort();                   
    }
            
    switch(mp->physics_ids[ia])
    {
      case MULTIPHYSICS_MECHANICAL:
        err += write_restart_mechanical(grid,fv,load,time_steps,opts,mp,myrank,ia,stepno,rs_path);
        break;
      case MULTIPHYSICS_THERMAL:
        err += write_restart_thermal(grid,fv,opts,myrank,ia,stepno,rs_path);
        break;
      case MULTIPHYSICS_CHEMICAL:
        // not yet implemented
        break;
      default:
        err += write_restart_mechanical(grid,fv,load,time_steps,opts,mp,myrank,ia,stepno,rs_path);
    }
  }
  
  err += write_time_step_info(fv,time_steps,opts,mp,myrank,stepno);  
  return err;
}
