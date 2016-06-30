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

#ifndef _Matrix_double
Define_Matrix(double);
#define _Matrix_double 1
#endif

#ifndef NO_VTK_LIB
#include "PGFem3D_to_VTK.hpp"

int read_initial_from_VTK(const PGFem3D_opt *opts, int myrank, int *restart, double *u0, double *u1)
{
  int err = 0;
  char filename[1024];
  
  sprintf(filename,"%s/restart/VTK/STEP_%.5d/%s_%d_%d.vtu",opts->opath,*restart,opts->ofname,myrank, *restart);   
  err += read_VTK_file(filename, u0);      
  sprintf(filename,"%s/VTK/STEP_%.5d/%s_%d_%d.vtu",opts->opath,*restart,opts->ofname,myrank, *restart);   
  err += read_VTK_file(filename, u1);
    
  return err;
}      

#else
int read_initial_from_VTK(const PGFem3D_opt *opts, int myrank, int *restart, double *u0, double *u1s)
{
  if(myrank==0)
  {
    PGFEM_printerr("Restart with VTK is not supported!\n");
    PGFEM_printerr("Enforce to turn off restart!\n");
  }
  
  *restart = -1;  
  return 0;
}
#endif

int read_restart_disp(double *u0, double *u1, const PGFem3D_opt *opts, 
                      int myrank, int nodeno, int nsd, int stepno)
{
 
  char filename[1024];
  sprintf(filename,"%s/restart/STEP_%.5d/%s_%d_%d.res",opts->opath,stepno,opts->ofname,myrank, stepno);   
  FILE *fp = fopen(filename,"r");

  if(fp == NULL)
  {    
    printf("Fail to open file [%s]. finishing\n", filename);      
    exit(0);  
  }
  
  for(int a=0; a<nodeno; a++)
  {
    for(int b=0; b<nsd; b++)
      fscanf(fp, "%lf %lf", u0+a*nsd + b, u1+a*nsd + b);    
  }
  
  fclose(fp);
  return 0;
}

int write_restart_disp(double *u0, double *u1, const PGFem3D_opt *opts, 
                       int myrank, int nodeno, int nsd, int stepno)
{
  
  char restart_path[1024];
  sprintf(restart_path, "%s/restart/STEP_%.5d", opts->opath,stepno);
                
  if(make_path(restart_path,DIR_MODE) != 0)
  {
    PGFEM_printf("Directory (%s) not created!\n",restart_path);
    abort();                   
  }  
 
  char filename[1024];
  sprintf(filename,"%s/restart/STEP_%.5d/%s_%d_%d.res",opts->opath,stepno,opts->ofname,myrank, stepno);   
  FILE *fp = fopen(filename,"w");

  if(fp == NULL)
  {    
    printf("Fail to open file [%s]. finishing\n", filename);
    exit(0);  
  }
  
  for(int a=0; a<nodeno; a++)
  {
    for(int b=0; b<nsd; b++)
      fprintf(fp, "%e %e ", u0[a*nsd + b], u1[a*nsd + b]);
    
    fprintf(fp, "\n");    
  }
  
  fclose(fp);
  return 0;
}

int write_restart_plasticity(double *u0, double *u1, const PGFem3D_opt *opts,
                             EPS *eps, const ELEMENT *elem,
                             int myrank, int elemno, int nodeno, int nsd, int stepno)
{
  int err = 0;

  char restart_path[1024];
  sprintf(restart_path, "%s/restart/STEP_%.5d", opts->opath,stepno);
                
  if(make_path(restart_path,DIR_MODE) != 0)
  {
    PGFEM_printf("Directory (%s) not created!\n",restart_path);
    abort();                   
  }  
 
  char filename[1024];
  sprintf(filename,"%s/restart/STEP_%.5d/%s_%d_%d.res",opts->opath,stepno,opts->ofname,myrank, stepno);
  FILE *fp = fopen(filename,"w");

  if(fp == NULL)
  {    
    printf("Fail to open file [%s]. finishing\n", filename);
    exit(0);  
  }
  
  for(int a=0; a<nodeno; a++)
  {
    for(int b=0; b<nsd; b++) {
      fprintf(fp, "%.17e %.17e ", u0[a*nsd + b], u1[a*nsd + b]);
    }
    fprintf(fp, "\n");    
  }
  
  for (int e = 0; e < elemno; e++)
  {
    const ELEMENT *p_el = elem + e;
    long n_ip = 0;
    int_point(p_el->toe,&n_ip);
    fprintf(fp, "%ld\n", n_ip);
    for (int ip = 0; ip < n_ip; ip++)
    {
      Constitutive_model *m = &(eps[e].model[ip]);
      err += m->param->write_restart(fp, m);
    }
  }

  fclose(fp);
  return err;
}

int read_restart_plasticity(double *u0, double *u1, const PGFem3D_opt *opts,
                             EPS *eps, const ELEMENT *elem,
                             int myrank, int elemno, int nodeno, int nsd, int stepno)
{
  int err = 0;
    
  char filename[1024];
  sprintf(filename,"%s/restart/STEP_%.5d/%s_%d_%d.res",opts->opath,stepno,opts->ofname,myrank, stepno);   
  FILE *fp = fopen(filename,"r");

  if(fp == NULL)
  {    
    printf("Fail to open file [%s]. finishing\n", filename);      
    exit(0);  
  }
  
  for(int a=0; a<nodeno; a++)
  {
    for(int b=0; b<nsd; b++)
      fscanf(fp, "%lf %lf", u0+a*nsd + b, u1+a*nsd + b);    
  }
      
  for (int e = 0; e < elemno; e++) 
  {
    const ELEMENT *p_el = elem + e;
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
      Constitutive_model *m = &(eps[e].model[ip]);
      err += m->param->read_restart(fp, m);
    }
  }
      
  fclose(fp);  
  return err;  
}


int read_restart(double *u0, double *u1, const PGFem3D_opt *opts, 
                 ELEMENT *elem, NODE *node, SIG * sig_e, EPS *eps, SUPP sup,
                 int myrank, int elemno, int nodeno, int nsd, int *stepno, double *tnm1)
{
  int err = 0;
  
  // write time stepping info
  char fn[1024];
  sprintf(fn, "%s/restart/STEP_%.5d/time_step_info.res",opts->opath,*stepno);

  double t[3];
  t[0] = t[1] = t[2] = -1.0;
  
  FILE *fp = fopen(fn, "r");    
  
  if(fp != NULL)
  { 
    fscanf(fp, "%lf %lf %lf\n", t+0, t+1, t+2);  
    *tnm1 = t[1];
    
    if(myrank==0)
      printf("read time stpe info %e %e %e\n", t[0], t[1], t[2]); 
    
    fclose(fp);
  }
  else
  {
    if(myrank==0)
      printf("WARNING: cannot read time steps info [%s] \n", fn);
  }
    
  switch(opts->analysis_type)
  {
    case DISP:
      err += read_restart_disp(u0,u1,opts,myrank,nodeno,nsd,*stepno);
      break;
    case CM:
      err += read_restart_plasticity(u0,u1,opts,eps,elem,
                                     myrank,elemno,nodeno,nsd,*stepno);
      break;
    default:
      read_initial_from_VTK(opts, myrank, stepno, u0, u1);
      break;
  }  
  return err;
}

int write_restart(double *u0, double *u1, const PGFem3D_opt *opts, 
                  ELEMENT *elem, NODE *node, SIG * sig_e, EPS *eps, SUPP sup,                  
                  int myrank, int elemno, int nodeno, int ndofn, int ndofd, int stepno, double *times)
{
  int err = 0;
  int nsd = 3;
  char restart_path[1024];
  sprintf(restart_path, "%s/restart", opts->opath);
  if(make_path(restart_path,DIR_MODE) != 0)
  {
    PGFEM_printf("Directory (%s) not created!\n",restart_path);
    abort();                   
  }
  
  switch(opts->analysis_type)
  {
    case DISP:
      write_restart_disp(u0,u1,opts,myrank,nodeno,ndofn,stepno);
      break;
  case CM:
    err += write_restart_plasticity(u0,u1,opts,eps,elem,
                                    myrank,elemno,nodeno,nsd,stepno);
    break;
    default:
    {  
      double *r_n_dof = (double *) malloc(sizeof(double)*ndofd);
      for(long a = 0; a<nodeno; a++)
      {              
        for(long b = 0; b<ndofn; b++)
        {
          long id = node[a].id[b];
          if(id>0)
          r_n_dof[id-1] = u0[a*ndofn + b];
        }
      }            
      VTK_print_vtu(restart_path,opts->ofname,stepno,
                    myrank,elemno,nodeno,node,elem,sup,r_n_dof,sig_e,eps,
                    opts);
      free(r_n_dof);              
      break;
    }
  }
  
  if(myrank==0)
  {
    // write time stepping info
    char fn[1024];
    sprintf(fn, "%s/STEP_%.5d/time_step_info.res", restart_path,stepno);
    
    FILE *fp = fopen(fn, "w");
    if(fp==NULL)
    {
      printf("Cannot create a file [%s]\n", fn);
      printf("Anyway continue ...\n");
    }
    else
    {  
      if(stepno>0)
        fprintf(fp, "%e %e %e\n", times[stepno-1], times[stepno], times[stepno+1]);  
      else
        fprintf(fp, "0.0 %e %e\n", times[stepno], times[stepno+1]);
                  
      fclose(fp);
    }
  }
      
  return err;
}
